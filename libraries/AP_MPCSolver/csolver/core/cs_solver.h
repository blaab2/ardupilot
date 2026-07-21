/* Public solver API — one cs_solver_iterate() call = ONE SQP-RTI iteration
 * (linearize via the generated model, full condensing, warm-started dense
 * DAQP solve, expansion, LM-damped full step), mirroring the acados
 * backend's per-call semantics (CUSTOM-SOLVER-PLAN.md sections 2.2, 4.2).
 *
 * Memory contract: cs_solver_init() receives ONE caller-provided block and
 * carves everything (model call buffers, condensing storage, dense QP data,
 * DAQP workspace) with the bump arena; every pointer is fixed after init and
 * no call path allocates (CS_NO_LIBC_ALLOC discipline). The required size is
 * static arithmetic: cs_solver_min_arena(N) — cross-checked at runtime by
 * cs_solver_high_water(). */
#ifndef CS_SOLVER_H
#define CS_SOLVER_H

#include "cs_types.h"

typedef struct cs_solver cs_solver;

/* Static workspace size for horizon N (must be a compiled model variant). */
size_t cs_solver_min_arena(int N);

/* Carve a solver out of mem[0..size). Also (re)initializes the cs_model call
 * buffers from the same block. Returns NULL on failure (*rc = CS_ERR_*). */
cs_solver *cs_solver_init(void *mem, size_t size, int N, int *rc);

/* Entry / terminal pins (nx each; r5 pins all six physical states).
 * M2 free-boundary support: a NaN component is FREE — entry NaNs become
 * decision columns, terminal NaNs drop their pin and gain the gated
 * boundary h/box rows (basin_map's partial-pin problem). A pattern change
 * re-dims the QP layout and clears the QP warm start. */
int cs_solver_set_pins(cs_solver *s, const cs_real *x0, const cs_real *xN);

/* Override the terminal box row of FREE terminal state idx (defaults: the
 * full state box; reset whenever the pin pattern changes). E.g. the P1
 * exit-speed inequality v_N in [0.999*V_INF, V_INF]. */
int cs_solver_set_xn_box(cs_solver *s, int idx, cs_real lo, cs_real hi);

/* Load the iterate: X (N+1)*nx (stage-major: X[k*nx+i]), U N*nu, T.
 * clip_rollout != 0 applies the plan-4.3 dynamics-consistent seeding: clip
 * U into the control boxes, clamp T into its box, then roll X out from X_0
 * through the transcription's own discrete map (~machine-zero defects). */
int cs_solver_set_seed(cs_solver *s, const cs_real *X, const cs_real *U,
                       const cs_real *T, int clip_rollout);

/* One RTI iteration. Returns CS_OK when the iteration ran (inspect
 * cs_solver_get_info for the QP status — a failed QP freezes the iterate,
 * the acados-backend convention), or CS_ERR_* on misuse. */
int cs_solver_iterate(cs_solver *s);

int cs_solver_get_iterate(const cs_solver *s, cs_real *X, cs_real *U,
                          cs_real *T);

/* step_norm: max(|dX|inf, |dU|inf, |dT|) of the last step (-1 if the last
 *            QP failed and the iterate was frozen);
 * defect:    max dynamics defect |F(X_k,U_k,T) - X_{k+1}| at linearization;
 * slack_max: largest elastic terminal slack (0 at a healthy r5 iterate);
 * kappa:     max/min diagonal ratio of the scaled dense H (cheap proxy,
 *            plan-4.4 instrumentation);
 * qp_status: DAQP exit flag (1 optimal, 2 soft-optimal, <=0 failure) or
 *            update_ldp error; qp_iters: DAQP active-set iterations. */
int cs_solver_get_info(const cs_solver *s, cs_real *step_norm,
                       cs_real *defect, cs_real *slack_max, cs_real *kappa,
                       int *qp_status, int *qp_iters);

int cs_solver_get_dims(const cs_solver *s, int *N, int *nx, int *nu,
                       int *nz, int *m);

/* Options (defaults: lm 1e-1, rho 1e-3, slack z 1e3 / Z 1e4, jerk_end 1). */
int cs_solver_set_opts(cs_solver *s, cs_real lm, cs_real rho,
                       cs_real slack_z, cs_real slack_Z, int jerk_end);

/* M1.4b: weights of the SHARED corridor-block + box-upper elastic slacks
 * (L1 zc / L2 Zc exact penalty; defaults 1e2/1e3 = the acados rate1 scale).
 * Step-E SITL showed re-plans buying 6-15 m of slack against the min-time
 * objective at those defaults - raise to make corridor violation expensive
 * relative to seconds of T. ABI 5. */
int cs_solver_set_soft_weights(cs_solver *s, cs_real zc, cs_real Zc);

/* Progress-aware corridor gating for the shrinking-horizon rh mode (M1.1).
 * tau0 < 0 (the init default) restores the batch full-horizon index gating.
 * tau0 in [0, 1]: stage k is gated by its ORIGINAL progress
 * tau_k = tau0 + (1 - tau0)*k/N, and each gated bound is relaxed by the
 * caller-measured corridor violation + margin (rlx_* >= 0, meters):
 * outbound  eta - apI(xi) <= rlx_out,  return  eta - reI(xi) >= -rlx_ret,
 * crossing  xi >= xi_cross_min - rlx_xi. Takes effect at the next
 * cs_solver_iterate() (bounds are rebuilt every linearization).
 *
 * CALLER CONTRACT (M1.3): each relaxation must be LEG-GATED by the vehicle's
 * own progress before it is passed in — rlx_out = violation only while
 * tau0 <= FRAC_OUT, rlx_ret only while tau0 >= FRAC_RET, rlx_xi only while
 * FRAC_OUT < tau0 < FRAC_RET (all + margin otherwise just the margin). The
 * relaxations exist ONLY to keep the hard x0 pin feasible while the vehicle
 * is on that leg; passing the raw violation un-gated lets a vehicle on the
 * outbound row relax the return-leg bound by the full row spacing, and the
 * re-plans legally cut meters inside the return envelope (see
 * RHCorridor.regate in experiments/sitl_closed_loop.py — the same rule). */
int cs_solver_set_gate(cs_solver *s, cs_real tau0, cs_real rlx_out,
                       cs_real rlx_ret, cs_real rlx_xi);

/* Runtime row-spacing shift for AUTO-detected turns (G0.2, ABI 6):
 * d_minus_14p1 = d - 14.1 (meters; the family band d in [13, 15] gives
 * [-1.1, +0.9]). The corridor d-generalization is ROW-RELATIVE
 * (D-GENERALIZATION-NOTES.md): NO re-codegen per d — the baked approach
 * interpolant is d-independent; this call rigidly shifts (a) the RETURN-leg
 * gate bounds (the reI row becomes eta - reI(xi) >= d_minus_14p1, batch and
 * rh gating alike) and (b) the eta box UPPER (the AP lateral peak,
 * CS_PD_BX eta entry + shift; also the free-entry eta column upper).
 * Terminal pins are caller data and must already be in the shifted frame.
 * Takes effect at the next cs_solver_iterate (bounds are rebuilt every
 * linearization); no QP re-dim, warm start kept. Default 0 = the canonical
 * problem bit-unchanged. NOT shifted: the eta scale (conditioning-only)
 * and the M2 free-terminal xn boxes (override via cs_solver_set_xn_box). */
int cs_solver_set_row_offset(cs_solver *s, cs_real d_minus_14p1);

/* Elastic (soft) corridor rows for the shrinking-horizon rh mode (M1.2).
 * on != 0 attaches one nonnegative L1+L2-penalized slack column to the ONE
 * gated corridor row of each interior stage (outbound eta<=apI, return
 * eta>=reI, or crossing xi>=0 — whichever the progress gate selects), so the
 * condensed QP is ALWAYS primal-feasible and DAQP never returns -1; the slack
 * magnitude (cs_solver_get_info slack_max) reports the residual corridor
 * violation. Off (default) = hard corridor rows, the M1/batch behavior
 * bit-unchanged. Re-dims the QP (like a pin-pattern change); pair with the
 * progress gate (cs_solver_set_gate tau0 >= 0). */
int cs_solver_set_soft_corr(cs_solver *s, int on);

/* R5-proximity reference tracking (ABI 7). Adds a Gauss-Newton state cost
 * w * sum_k ||x_k - x_ref_k||^2 to every subsequent QP build (see the
 * cs_condense.h w_ref/ref_mode comment): mode 1 = L2 (Hessian curvature ->
 * lowers the kappa proxy of the flat bulge<->time mode + a pull toward the
 * reference, the recommended damping of the stochastic turn-2 over-bulge);
 * mode 2 = L1 surrogate (gradient only, no curvature -- the conditioning-null
 * comparison); w = 0 or mode = 0 disables it (canonical problem bit-unchanged).
 * x_ref (may be NULL to update only w/mode) is the reference R5 seed states,
 * (N+1)*nx PHYSICAL units, restretched by the caller to the current progress
 * (the rh wrapper / Python rh_step restretches the init seed each cycle).
 * Copied in; not retained. */
int cs_solver_set_ref(cs_solver *s, const cs_real *x_ref, cs_real w, int mode);
/* R5-anchored one-sided soft CROP floor: xi_k >= min(0, env(eta_k)) - margin
 * over the whole transfer, env = the reference path's minimum xi at the
 * stage's row-transfer position eta_k (eta-aligned lookup on qp->x_ref).
 * Holds the plan out of the inter-row crop no deeper than the restretched R5
 * reference cuts, with no bulge push. Reads qp->x_ref, so the caller must
 * keep it fresh (cs_solver_set_ref / cs_rh_step). on=0 restores the legacy
 * vx crossing gate; margin >= 0. */
int cs_solver_set_crop_bound(cs_solver *s, int on, cs_real margin);

/* Largest M1.2 soft-state slack (corridor-block + box-upper) of the last
 * iterate — the shrinking-horizon residual corridor/box violation, in
 * physical units (0 = the plan is inside the gated corridor and boxes).
 * Distinct from cs_solver_get_info's slack (= the terminal elastic pins). */
int cs_solver_get_soft_slack(const cs_solver *s, cs_real *soft);

/* Diagnostic dump of the CURRENT condensed QP (valid after cs_solver_iterate).
 * Copies the mA general rows A (row-major, mA*nz), the full two-sided bounds
 * bl/bu (length m = nz simple bounds then mA general rows) and sense (length
 * m). Any out pointer may be NULL. Sizes come from cs_solver_get_dims
 * (mA = m - nz). Used by the offline rh probe to attribute infeasibility to a
 * constraint group. */
int cs_solver_get_qp(const cs_solver *s, cs_real *A, cs_real *bl, cs_real *bu,
                     int *sense);

/* Arena bytes actually carved (== high-water mark; nothing allocates after
 * init). The measured M1 memory-gate number. */
size_t cs_solver_high_water(const cs_solver *s);

#endif /* CS_SOLVER_H */
