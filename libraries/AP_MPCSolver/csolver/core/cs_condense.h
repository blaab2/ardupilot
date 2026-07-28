/* Full condensing of the r5 stage QP onto the dense control+time space
 * (CUSTOM-SOLVER-PLAN.md sections 4.1-4.4).
 *
 * QP variable (SCALED, see cs_scaling.h):
 *     zeta = [du_0/su .. du_{N-1}/su, ds_T, s_lo[npin], s_up[npin]]
 * nz0 = N*nu + 1 real columns + ns = 2*npin elastic terminal slack columns.
 * State increments are eliminated exactly:
 *     dx_{k+1} = A_k dx_k + B_k du_k + bT_k dT + c_k,   dx_0 = x0 - X_0
 *  => dx_k = E_k zeta + e_k    (E_k: nx x nz0 col-major, stored all k).
 *
 * Objective (the acados twin's, exactly): terminal 0.5*T^2 (GN-exact min-T),
 * stage (rho/N)*0.5*||u_k||^2 (acados' default tf/N cost scaling), plus the
 * EXACT full condensing of acados' Levenberg-Marquardt damping
 * 0.5*lm*(||dx~_k||^2 + ||du_k||^2) per stage with the augmented state
 * x~ = [x; T]: lm*(E'E + T_ref^2 on s_T) into H and lm*E'e into f. LM only
 * shapes the iteration path — at a fixed point dz = 0 is stationary
 * regardless — but the path is what the measured "LM 1e-1 fixed point"
 * finding is about, so it is replicated, not approximated.
 *
 * Rows (two-sided, DAQP layout: simple bounds on all nz variables first,
 * then mA general rows, A row-major):
 *   stage 0:        [jerk node, jerk end-of-interval]           (x0 pinned;
 *                   the x-only h rows are constants there — omitted to avoid
 *                   zero rows, the acados stage-0 pins make them vacuous)
 *   stage 1..N-1:   [5 h rows (gated per _gate_corridor), jerk end, box rows
 *                   for the both-bounds-finite states (acados idxbx filter)]
 *   terminal:       npin elastic equality pins with +-1 slack columns.
 * The conservative end-of-interval jerk row enforces
 *     adot^2 + (a + adot*T/N)^2 thetadot^2 <= j_max^2
 * — a(t) is linear inside the interval (u piecewise constant), so jerk^2 is
 * convex in t and node row + end row bound the WHOLE interval (the N=10
 * between-node audit fix, MEM-SWEEP). Linearized analytically here (it is a
 * solver-layer row, not a model h row).
 *
 * All storage arena-carved at init; build/expand allocate nothing. */
#ifndef CS_CONDENSE_H
#define CS_CONDENSE_H

#include "cs_arena.h"
#include "cs_scaling.h"
#include "cs_types.h"

typedef struct {
    /* dims (layout-dependent: recomputed by cs_condense_set_pattern; the
     * default pattern — entry fully pinned, all nx terminal pins — is the
     * M1 layout exactly) */
    int N, nx, nu, nh, npin;
    int nz0, ns, nz, mA, m;
    int nsc;                   /* elastic soft-state slack columns (M1.2):
                                * CS_NSC + CS_PD_NBX when soft_corr, else 0.
                                * A shared nonneg slack per gated corridor
                                * BLOCK (outbound / crossing / return — the
                                * gate relaxes each block uniformly, so one
                                * slack absorbing the block's worst residual
                                * frees every row) plus one per state BOX
                                * softening its UPPER bound (v<=V_INF,
                                * a<=A_CAP ride their ceilings near the exit /
                                * apex; the shrinking horizon overshoots them
                                * by cm). Laid out after the 2*npin terminal
                                * slacks: [out, cross, ret, box_v, box_a].   */
#define CS_NSC 3               /* shared corridor-block slacks (M1.2)       */
    int n_out, n_ret;
    int iT;                    /* s_T column index (= N*nu)               */

    /* boundary pin pattern (M2 free-boundary support, basin_map):
     * free0 lists the FREE entry states (each becomes one extra scaled
     * dx0 decision column between s_T and the slacks); pinN flags which
     * terminal states carry an elastic pin. Free terminal states get
     * gated x-only h rows + box rows (xn_lo/xn_hi, default = the full
     * state box; override via cs_solver_set_xn_box — e.g. the P1 exit
     * v >= 0.999 V_INF inequality). */
    int free0[6], nf0;
    int pinN[6];
    cs_real xn_lo[6], xn_hi[6];
    int nz_max, mA_max, m_max; /* worst-case layout (arena sizing)        */

    /* options (defaults in cs_condense_init; setters via cs_solver) */
    cs_real lm;        /* Levenberg-Marquardt damping, default 1e-1 */
    cs_real rho;       /* control regularization, default 1e-3      */
    cs_real slack_z;   /* terminal slack L1 weight (exact penalty)  */
    cs_real slack_Z;   /* terminal slack L2 weight                  */
    cs_real slack_zc;  /* corridor slack L1 weight (M1.2, acados    */
    cs_real slack_Zc;  /* rate1 corridor: zl 1e2 / Zl 1e3)          */
    cs_real relax;     /* ungated-row bound (acados _RELAX = 1e9)   */
    int jerk_end;      /* conservative end-of-interval jerk row on  */
    int soft_corr;     /* elastic corridor rows in rh mode (M1.2)   */

    /* progress-aware corridor gating (M1.1, rh mode). gate_tau0 < 0
     * (default) keeps the original full-horizon INDEX gating (n_out/n_ret).
     * gate_tau0 in [0, 1]: the iterate is a shrinking-horizon restretch
     * whose stage k had ORIGINAL progress tau_k = tau0 + (1 - tau0)*k/N —
     * gate by tau_k against the same fractions, and relax each gated bound
     * by the caller-measured violation + margin (softening-to-the-measured-
     * state; the acados twin needed exactly this experiment-side, see
     * sitl_closed_loop.py RHCorridor / SITL-CLOSED-LOOP-NOTES.md). */
    cs_real gate_tau0;
    cs_real gate_rlx_out, gate_rlx_ret, gate_rlx_xi;

    /* runtime row-spacing shift (G0.2, AUTO/MpcTurn P0): row_off = d - 14.1.
     * The corridor d-generalization is ROW-RELATIVE (D-GENERALIZATION-
     * NOTES.md): the approach interpolant is d-independent; the RETURN-leg
     * rows and the eta box upper shift RIGIDLY by row_off. Applied at build
     * time to (a) the reI gate bounds (return row becomes
     * eta - reI(xi) >= row_off when gated, batch and rh alike), (b) the
     * interior eta box-row upper (CS_PD_BX eta entry + row_off) and (c) the
     * free-entry eta column upper (M2 pattern). Stage-N free-terminal boxes
     * (xn_lo/xn_hi, M2 experiments only) are NOT auto-shifted — override
     * via cs_solver_set_xn_box if a free-terminal-eta problem needs it.
     * The eta SCALE stays the canonical 14.1 (conditioning-only; < 8 %
     * variation over d in [13, 15]). Default 0 = the canonical problem
     * bit-unchanged. */
    cs_real row_off;

    /* R5-proximity reference tracking (turn-2 conditioning fix). Adds a
     * Gauss-Newton state cost w_ref * sum_k ||x_k - x_ref_k||^2 to the QP,
     * x_ref = the reference R5 seed restretched to the current progress:
     *   H += w_ref * E~' diag(lm_w) E~   (curvature in the flat bulge<->time
     *                                     mode -> raises the small diagonal,
     *                                     lowers the kappa proxy)
     *   f += w_ref * E~' diag(lm_w) r~,  r~_k = (X_k - x_ref_k)/sx  (the pull
     *                                     back toward R5 that damps the
     *                                     stochastic over-bulge)
     * Same syrk/gemv structure as the LM term (which anchors the scaled
     * defect); this anchors the deviation-from-reference instead. ref_mode:
     * 0 off (default, canonical problem bit-unchanged), 1 L2 (H+f), 2 L1-
     * surrogate (f += w_ref*E~'diag(lm_w) sign(r~) only, no curvature -- the
     * conditioning-null comparison for the L1-vs-L2 sweep), 3 NORMAL-only
     * path tube (ABI 10): per node, the signed distance from (xi_k, eta_k)
     * to the nearest point of the x_ref (xi, eta) POLYLINE, penalized L2
     * rank-1 along the local path normal only. Unlike modes 1/2 the residual
     * cannot be re-anchored away by a rigid along-path translation (the
     * flare-ratchet null direction): tangential progress and timing stay
     * free, normal migration of the flare is priced in metres. x_ref/w_ref/
     * ref_mode set via cs_solver_set_ref. */
    cs_real w_ref;
    int ref_mode;

    /* R5-anchored CROP floor (one-sided, soft): xi_k >= min(0,
     * env(eta_k)) - crop_margin over the whole transfer, env = the
     * reference path's minimum xi at the stage's current row-transfer
     * position eta_k (crop_env_xi — eta-aligned, NOT node-aligned: a
     * cutting plan is ahead of the reference in eta, so a per-node bound
     * never binds). Holds the plan out of the inter-row crop (the area
     * between the straights and the corner line) no deeper than R5 cuts,
     * WITHOUT pushing a bulge (capped at 0 where R5 sits in the headland).
     * Loose near the rows where R5 itself flares. Softened by the shared
     * c_cross slack wherever it lifts lo[4] above -gth. Replaces the
     * relaxable vx crossing gate when on. Uses x_ref (the same restretched
     * R5 seed as the w_ref proximity), set by cs_rh_step whenever
     * crop_bound || w_ref>0. */
    int crop_bound;
    cs_real crop_margin;
    cs_real crop_v_xi, crop_v_eta;  /* vehicle state at build (x0[2], x0[3]):
                                     * stages within CS_CROP_VEH_WIN of the
                                     * vehicle's eta admit the measured xi
                                     * (hard-pin accommodation, LOCAL only —
                                     * a global relaxation would re-open the
                                     * whole crossing, the old vx bug) */
    cs_real crop_eta_max;           /* ref eta max (build) for the row
                                     * dead-zone: floor off within
                                     * CS_CROP_ETA_DEAD of either row */
    cs_real crop_xi_up;             /* soft xi upper cap (build): ref bulge
                                     * + CS_CROP_UP_MARGIN — damps the flat
                                     * bulge<->time manifold wander that the
                                     * validator otherwise rejects (a cap,
                                     * not a pull) */

    cs_scales sc;
    /* LM state-metric weights over the SCALED E rows (M2): lm_w[i] = sx_i^2
     * reproduces the physical-metric acados damping EXACTLY; lm_w[i] = 1 is
     * the scaled-metric variant. Default set in cs_condense_init (measured
     * choice — see csolver/README.md M2). */
    cs_real lm_w[6];

    /* dense QP data (DAQP conventions) */
    cs_real *H;        /* nz x nz col-major, symmetric              */
    cs_real *f;        /* nz                                        */
    cs_real *A;        /* mA x nz ROW-major                         */
    cs_real *bl, *bu;  /* m = nz + mA (simple bounds then rows)     */
    int *sense;        /* m (zeroed each build; solver adds bits)   */

    /* condensing state */
    cs_real *E;        /* (N+1) blocks of nx*nz0, col-major         */
    cs_real *e;        /* (N+1) * nx                                */
    cs_real *x_ref;    /* (N+1)*nx R5 reference states (phys units) */
    cs_real *rref;     /* nx scratch: scaled iterate-vs-ref residual */
    cs_real *vref;     /* nz_max scratch: mode-3 normal-row Jacobian */

    /* scratch */
    cs_real *Ak, *Bk, *bTk, *xnext, *hk, *Jxk, *Juk, *jrow;

    /* diagnostics of the last build */
    cs_real defect_max;
} cs_qp;

/* Carve all buffers from the arena (sized for the WORST-CASE pattern
 * unless CS_PINNED_ONLY is defined — the target build's lever: the flight
 * problem is fully pinned and pays no free-boundary reserve). Installs the
 * default fully-pinned pattern. Returns CS_OK / CS_ERR_ARENA. */
int cs_condense_init(cs_qp *qp, cs_arena *arena, int N);

/* Install a boundary pattern: free0/pinN are per-state flags (nonzero =
 * entry state FREE / terminal state PINNED). Recomputes the column/row
 * layout inside the init-time allocation; resets xn_lo/xn_hi to the full
 * state box. Returns CS_ERR_ARG if the pattern exceeds the allocated
 * worst case (only possible under CS_PINNED_ONLY). */
int cs_condense_set_pattern(cs_qp *qp, const int *free0_flags,
                            const int *pinN_flags);

/* Toggle the elastic corridor slacks (M1.2). Recomputes ns/nz/m in place
 * (nsc = on ? N-1 : 0); mA is unchanged (slacks are columns, not rows).
 * Returns CS_ERR_ARG if the new nz exceeds the allocated worst case. */
int cs_condense_set_soft(cs_qp *qp, int on);
int cs_condense_set_crop(cs_qp *qp, int on, cs_real margin);

/* Linearize at (X, U, T) with pins (x0, xN) and fill H/f/A/bl/bu/sense,
 * E/e, defect_max. Calls the cs_model layer (must be initialized). */
int cs_condense_build(cs_qp *qp, const cs_real *X, const cs_real *U,
                      cs_real T, const cs_real *x0, const cs_real *xN);

/* Expand the QP solution zeta (length nz) back: X/U/T updated IN PLACE with
 * the full (RTI) step; returns the step norm max(|dX|, |dU|, |dT|). */
cs_real cs_condense_expand(const cs_qp *qp, const cs_real *zeta,
                           cs_real *X, cs_real *U, cs_real *T);

#endif /* CS_CONDENSE_H */
