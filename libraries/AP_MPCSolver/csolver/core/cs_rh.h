/* cs_rh — the embedded receding-horizon wrapper (D0, ABI 6): portable-C port
 * of the HOST-side rh cycle logic that so far lived in Python, namely
 *   solvers/csolver_solver.py rh_step        (fractional polyline-projection
 *                                             restretch, tau0 bookkeeping,
 *                                             leg-membership-gated corridor
 *                                             relaxations)
 *   experiments/sitl_closed_loop.py          (plan validation, accept/reject
 *                                             with iterate restore + tau0
 *                                             rollback, SpiralDetector DEEP/
 *                                             TSTALL recovery, T_rem aging)
 * so the firmware MpcTurn feature drives ONE C entry point per MPC cycle.
 * Behavioral parity with the Python twin is gated by tests/rh_parity.py
 * (identical accept/reject/tau0/T sequences over a recorded walk).
 *
 * Memory contract: cs_rh is a plain caller-allocated struct (static or
 * arena) with FIXED internal arrays (N <= CS_RH_NMAX = 30, the flight
 * horizon); no call path allocates (CS_NO_LIBC_ALLOC discipline). The
 * wrapped cs_solver keeps its own arena and is NOT owned.
 *
 * Per-cycle protocol (cs_rh_step):
 *  1. SpiralDetector update from the measured state, BEFORE the solve
 *     (recovery must solve from the reference-family seed, not the degraded
 *     iterate):
 *       DEEP    reference-speed deficit > 1.5 m/s for 3 consecutive cycles,
 *               progress tau in (0.05, 0.97)  -> reseed from the init seed
 *               restretched to tau (event CS_RH_EVENT_RESEED)
 *       TSTALL  accepted-plan T_rem drop < 0.03 s for 8 consecutive cycles,
 *               T_rem <= 1.5, tau >= 0.6      -> the vehicle is at/past the
 *               window exit: report CS_RH_EVENT_EXIT_FORCED and DO NOT
 *               solve — the caller performs the handback (both rules force
 *               the exit when tau > 0.95 or xi is past the seed exit).
 *       cooldown 2 s, max 2 recoveries per turn.
 *  2. shrinking-horizon restretch: project x_meas onto the plan's XY
 *     polyline -> FRACTIONAL tau*; advance tau0 += (1-tau0)*tau*, restretch
 *     the iterate onto [0,1], T *= (1-tau*) (resolution-independent — the
 *     integer nearest-node rule froze tau0 at the embedded N).
 *  3. leg-membership-gated corridor relaxations (M1.3): measured violations
 *     from cs_h rows 2/3 at x_meas WITH the row offset applied
 *     (v_out = max(0, h[2]) only while tau0 <= FRAC_OUT, v_ret =
 *     max(0, row_off - h[3]) only while tau0 >= FRAC_RET, crossing between;
 *     each + the 0.05 m margin), via cs_solver_set_gate.
 *  4. pin x_meas, run the RTI iterations (base k, escalated to >= 5 on a
 *     > 1 m / > 1 m/s innovation vs the nearest plan node, >= 10 on a
 *     recovery cycle — the M1.4/B2 lever).
 *  5. validate the plan: finite, 0.05 < T < 15, max xi < 5,
 *     max eta < 17 + row_off, min v > 0, DAQP status > 0, soft corridor
 *     slack <= 1.0. Accept -> it becomes the streamed plan (T_rem restarts).
 *     Reject -> restore the ACCEPTED plan as the iterate (terminal pin
 *     follows the accepted plan's last node, the Python rh_reset semantics)
 *     and roll tau0 back (the restretch advance of a rejected cycle would
 *     otherwise re-add progress every cycle of a rejection burst).
 *  6. T_rem = accepted T - time since accept; the caller exits MPC below
 *     its T_EXIT and on CS_RH_EVENT_EXIT_FORCED. */
#ifndef CS_RH_H
#define CS_RH_H

#include "cs_types.h"
#include "cs_solver.h"

#define CS_RH_NMAX 30            /* flight horizon (N = 30 fixed is fine)  */
#define CS_RH_D_CANON ((cs_real)14.1)   /* canonical row spacing [m]       */
#define CS_RH_MARGIN ((cs_real)0.05)    /* gate softening margin [m]       */
#define CS_RH_MAX_EVENTS 8

/* cs_rh_step *event out */
enum {
    CS_RH_EVENT_NONE = 0,
    CS_RH_EVENT_RESEED = 1,      /* recovery reseed happened this cycle    */
    CS_RH_EVENT_EXIT_FORCED = 2, /* caller must hand back NOW (no solve)   */
};

/* spiral rules (cs_rh_event.rule) */
enum { CS_RH_RULE_DEEP = 1, CS_RH_RULE_TSTALL = 2 };

/* last_reject bitmask (diagnosability — the audit-B lesson: log WHY) */
enum {
    CS_RH_REJ_NONFINITE = 1,
    CS_RH_REJ_T = 2,
    CS_RH_REJ_XI = 4,
    CS_RH_REJ_ETA = 8,
    CS_RH_REJ_VMIN = 16,
    CS_RH_REJ_QP = 32,
    CS_RH_REJ_SLACK = 64,
};

typedef struct {
    cs_real t;                   /* wrapper-relative trigger time [s]      */
    int rule;                    /* CS_RH_RULE_*                           */
    int action;                  /* the event code it produced             */
    cs_real tau, deficit, v, dist;
} cs_rh_event;

/* All state in the open struct so the caller can place it statically;
 * treat the fields as read-only diagnostics — mutate only through the API. */
typedef struct {
    cs_solver *s;                /* wrapped solver (not owned)             */
    int N, nx, nu;
    int k_iters;                 /* base RTI iterations per cycle          */
    cs_real d_off;               /* row offset d - 14.1 (cs_rh_set_frame_d)*/

    /* the init seed (family member, embedded N+1 nodes) */
    cs_real S[(CS_RH_NMAX + 1) * 6];
    cs_real T_seed;
    cs_real engage_xi;           /* = S[0][2]: engage the turn at xi >= it */
    cs_real xN[6];               /* current terminal pin (see step 5)      */

    /* accepted plan (the compact snapshot the reject path restores) */
    cs_real Xa[(CS_RH_NMAX + 1) * 6];
    cs_real Ua[CS_RH_NMAX * 2];
    cs_real Ta;
    cs_real age;                 /* seconds since the accept               */
    int age_valid;

    cs_real tau0;                /* cumulative shrinking-horizon progress  */

    /* SpiralDetector state */
    int c_deep, c_stall;
    cs_real prev_trem;
    int prev_trem_valid;
    int recoveries;
    cs_real cool_until;
    cs_real sp_deficit, sp_tau, sp_dist;   /* last-update diagnostics      */

    cs_real t_now;               /* wrapper-relative clock (sum of dt)     */
    int n_cycles, n_fallbacks;
    int last_reject;             /* CS_RH_REJ_* mask of the last cycle     */

    cs_rh_event events[CS_RH_MAX_EVENTS];
    int n_events;

    /* R5-proximity reference tracking (turn-2 over-bulge fix,
     * cs_rh_set_ref_tracking). Each step the init seed rh->S is restretched
     * to the current progress tau0 and applied via cs_solver_set_ref as the
     * L2 tracking reference (ref_w>0). val_slack_max is the runtime slack
     * reject cap (VAL_SLACK_MAX default); the fix pairs w=3 with a loosened
     * 1.5 cap (the over-bulge and the corridor-slack are two facets of the
     * same lag). Defaults ref_w=0 / cap=VAL_SLACK_MAX = the canonical rh
     * bit-unchanged (rh_parity holds at defaults). */
    cs_real ref_w;
    int ref_mode;
    cs_real val_slack_max;
    cs_real Xref[(CS_RH_NMAX + 1) * 6];

    /* R5-anchored one-sided soft CROP floor (cs_rh_set_crop_bound): each step
     * the restretched Xref is installed as the solver reference (like the
     * proximity, but the pull weight may be 0) and cs_solver_set_crop_bound
     * applies xi_k >= min(0, Xref_xi_k) - crop_margin. Holds the plan out of
     * the inter-row crop no deeper than R5 cuts, without a bulge push. Off by
     * default (crop_bound=0) => canonical rh, rh_parity bit-unchanged. */
    int crop_bound;
    cs_real crop_margin;

    /* scratch (restretch / iterate readback) */
    cs_real Xs[(CS_RH_NMAX + 1) * 6];
    cs_real Us[CS_RH_NMAX * 2];
    cs_real Xi_[(CS_RH_NMAX + 1) * 6];
    cs_real Ui_[CS_RH_NMAX * 2];
} cs_rh;

/* Bind the wrapper to an initialized solver and load the seed member:
 * S = embedded seed states ((N+1)*6 stage-major, e.g. a cs_seed_family.h
 * row), T_seed its turn time, N the horizon (must equal the solver's N,
 * <= CS_RH_NMAX), k_iters the base RTI iterations per cycle (>= 1).
 * Arms the mission: pins (S_0, S_N), elastic corridor ON, progress gate
 * tau0 = 0 at the margin, seed loaded as the iterate, accepted plan = the
 * seed with finite-difference controls, row offset RESET to canonical
 * (call cs_rh_set_frame_d afterwards for a detected d != 14.1).
 * Seed preconvergence is the caller's: run cs_solver_iterate() during row
 * transit (the harness's <= 30-iteration RUN_UP budget) before engaging. */
int cs_rh_init(cs_rh *rh, cs_solver *s, const cs_real *S, cs_real T_seed,
               int N, int k_iters);

/* Set the detected row spacing d (ABSOLUTE meters; family band [13, 15]).
 * Stores d - 14.1 and forwards it to cs_solver_set_row_offset (G0.2), so
 * the return-leg gate bounds, the eta box upper, the measured-violation
 * relaxations and the eta_max validation all live in the shifted frame.
 * Call after cs_rh_init and before the first cs_rh_step. Handedness is NOT
 * handled here: the caller mirrors a right-handed turn's states into the
 * canonical frame before every call and mirrors plans back (G0.4). */
int cs_rh_set_frame_d(cs_rh *rh, cs_real d);

/* Arm the R5-proximity reference tracking (turn-2 over-bulge fix). w = the L2
 * reference weight (0 disables, the canonical rh); mode 1 = L2 / 2 = L1
 * surrogate; slack_max = the runtime slack reject cap (<= 0 keeps the current
 * value). Installs the current turn's seed reference immediately so the raw
 * ARMED/READY cs_solver_iterate() preconvergence uses the same objective as
 * ENGAGED cs_rh_step(); cs_rh_step restretches it thereafter. Call after
 * cs_rh_init (which clears the previous turn's solver-side reference). */
int cs_rh_set_ref_tracking(cs_rh *rh, cs_real w, int mode, cs_real slack_max);
int cs_rh_set_crop_bound(cs_rh *rh, int on, cs_real margin);

/* READY-phase reference alignment: install the reference restretched to the
 * VEHICLE's absolute progress along the seed path, allowing pre-path states
 * (negative tau by exact backward extrapolation of the constant-cruise
 * approach segment). Call before every live-pinned tracking solve (the
 * firmware's thread_track): a tau0=0-anchored reference there is node-
 * misaligned and teaches the iterate a premature braking profile. No-op
 * when reference tracking is off. Engaged cycles keep their own per-cycle
 * restretch in cs_rh_step (tau0 >= 0 by the engage gate). */
int cs_rh_arm_ready_ref(cs_rh *rh, const cs_real x_meas[6]);

/* The engage gate from the seed: enter the MPC turn at xi >= this value. */
cs_real cs_rh_engage_xi(const cs_rh *rh);

/* One MPC cycle (the full protocol above). x_meas = measured embedded state
 * [v, psi, xi, eta, a, theta] in the CANONICAL frame (caller applies the
 * harness clips and the plan-carried a/theta/psi-branch); dt = seconds since
 * the previous cs_rh_step (0 on the first). Outputs (any may be NULL):
 * X_plan/U_plan/T_plan = the ACCEPTED plan (the one to stream/track),
 * T_rem = its remaining time, accepted = 1 if THIS cycle's solve was
 * accepted (0 = fallback kept the previous plan; see cs_rh_last_reject),
 * event = CS_RH_EVENT_*. On CS_RH_EVENT_EXIT_FORCED no solve ran and the
 * caller must hand control back to its stock guidance. */
int cs_rh_step(cs_rh *rh, const cs_real *x_meas, cs_real dt,
               cs_real *X_plan, cs_real *U_plan, cs_real *T_plan,
               cs_real *T_rem, int *accepted, int *event);

/* Copy up to max_events recovery events; *n = number copied. */
int cs_rh_events(const cs_rh *rh, cs_rh_event *out, int max_events, int *n);

cs_real cs_rh_progress(const cs_rh *rh);     /* cumulative tau0            */
int cs_rh_last_reject(const cs_rh *rh);      /* CS_RH_REJ_* mask (0 = ok)  */
size_t cs_rh_sizeof(void);                   /* sizeof(cs_rh) for bindings */

#endif /* CS_RH_H */
