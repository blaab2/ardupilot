/* cs_rh.c — embedded receding-horizon wrapper (D0). Behavioral port of
 * solvers/csolver_solver.py rh_step + the sitl_closed_loop.py csolver TURN
 * branch (SpiralDetector, plan validation, accept/reject, T_rem aging).
 * Parity with the Python twin is gated by tests/rh_parity.py. */
#include "cs_rh.h"

#include "cs_linalg.h"
#include "cs_model.h"

#include "../model/generated/turn_r5_probdata.h"

#include <math.h>

/* SpiralDetector tuning (sitl_closed_loop.py class SpiralDetector, CAND-2) */
#define SP_DEEP_DFC ((cs_real)1.5)
#define SP_DEEP_N 3
#define SP_STALL_DT ((cs_real)0.03)
#define SP_STALL_N 8
#define SP_STALL_TREM ((cs_real)1.5)
#define SP_STALL_TAU ((cs_real)0.6)
#define SP_TAU_LO ((cs_real)0.05)
#define SP_TAU_HI ((cs_real)0.97)
#define SP_COOLDOWN ((cs_real)2.0)
#define SP_MAX_RECOV 2

/* plan validation (the harness's reject terms) */
#define VAL_T_LO ((cs_real)0.05)
#define VAL_T_HI ((cs_real)15.0)
#define VAL_XI_MAX ((cs_real)5.0)
#define VAL_ETA_MAX ((cs_real)17.0)     /* + row offset at runtime */
#define VAL_SLACK_MAX ((cs_real)1.0)

/* forced-exit + iteration escalation (M1.4/B2, CAND-2) */
#define FE_TAU ((cs_real)0.95)
#define FE_XI_PAD ((cs_real)0.3)
#define INNOV_GATE ((cs_real)1.0)
#define K_INNOV 5
#define K_RECOV 10
/* READY-phase leg-gate relaxation (cs_rh_arm_ready_ref): large enough to
 * make the corridor leg gates inert for the pre-path warm start, far below
 * the gth slack trigger (relax/2 = 5e8) so no elastic machinery engages. */
#define CS_RH_READY_RLX ((cs_real)50.0)

/* np.interp over the uniform node grid linspace(0, 1, n+1); col strided
 * (stage-major state columns). Clamps outside [0, 1] like np.interp —
 * pre-path (tq < 0) reference nodes are NOT segment-extrapolated: the seed
 * decelerates from node 0 (braking starts AT the engage gate), so linear
 * extrapolation would inflate the backward reference's speed; the READY
 * extension (cs_rh_arm_ready_ref) instead holds the node-0 state and
 * extends only xi linearly (cruise-hold along the row). */
static cs_real interp_node(const cs_real *col, int stride, int n, cs_real tq)
{
    cs_real pos = tq * (cs_real)n;
    cs_real w;
    int j;
    if (pos <= (cs_real)0)
        return col[0];
    if (pos >= (cs_real)n)
        return col[(size_t)n * stride];
    j = (int)pos;
    if (j > n - 1)
        j = n - 1;
    w = pos - (cs_real)j;
    return col[(size_t)j * stride]
        + (col[(size_t)(j + 1) * stride] - col[(size_t)j * stride]) * w;
}

/* np.interp over the CONTROL node grid (positions j/N, j = 0..N-1): the
 * rh_step restretch's tu = tau[:-1] — queries beyond (N-1)/N clamp. */
static cs_real interp_unode(const cs_real *col, int stride, int N, cs_real tq)
{
    cs_real pos = tq * (cs_real)N;
    cs_real w;
    int j;
    if (pos <= (cs_real)0)
        return col[0];
    if (pos >= (cs_real)(N - 1))
        return col[(size_t)(N - 1) * stride];
    j = (int)pos;
    if (j > N - 2)
        j = N - 2;
    w = pos - (cs_real)j;
    return col[(size_t)j * stride]
        + (col[(size_t)(j + 1) * stride] - col[(size_t)j * stride]) * w;
}

/* fractional projection of (xi, eta) onto a plan's XY polyline (the
 * resolution-independent progress rule of rh_step / SpiralDetector.project):
 * per-segment clipped projection, first argmin, tau = (j + s_j) / N. */
static void project_xy(const cs_real *X, int nx, int N, cs_real xi,
                       cs_real eta, cs_real *tau_out, cs_real *dist_out)
{
    cs_real best = (cs_real)0, sbest = (cs_real)0;
    int j, jn = 0, first = 1;
    for (j = 0; j < N; ++j) {
        const cs_real px = X[(size_t)j * nx + 2];
        const cs_real py = X[(size_t)j * nx + 3];
        const cs_real sx = X[(size_t)(j + 1) * nx + 2] - px;
        const cs_real sy = X[(size_t)(j + 1) * nx + 3] - py;
        cs_real L2 = sx * sx + sy * sy;
        cs_real s, dx, dy, d2;
        if (L2 < (cs_real)1e-12)
            L2 = (cs_real)1e-12;
        s = ((xi - px) * sx + (eta - py) * sy) / L2;
        if (s < (cs_real)0)
            s = (cs_real)0;
        if (s > (cs_real)1)
            s = (cs_real)1;
        dx = xi - (px + s * sx);
        dy = eta - (py + s * sy);
        d2 = dx * dx + dy * dy;
        if (first || d2 < best) {
            best = d2;
            jn = j;
            sbest = s;
            first = 0;
        }
    }
    *tau_out = ((cs_real)jn + sbest) / (cs_real)N;
    *dist_out = (cs_real)sqrt((double)best);
}

/* restretch a plan from fractional progress tau_r onto [0, 1]: X columns
 * over the node grid, U columns over the control grid (rh_step's interp). */
static void restretch(const cs_real *X, const cs_real *U, int N, int nx,
                      int nu, cs_real tau_r, cs_real *Xo, cs_real *Uo)
{
    int k, j;
    for (k = 0; k <= N; ++k) {
        const cs_real tq = tau_r + ((cs_real)1 - tau_r)
            * (cs_real)k / (cs_real)N;
        for (j = 0; j < nx; ++j)
            Xo[(size_t)k * nx + j] = interp_node(X + j, nx, N, tq);
        if (k < N)
            for (j = 0; j < nu; ++j)
                Uo[(size_t)k * nu + j] = interp_unode(U + j, nu, N, tq);
    }
}

static void default_controls(const cs_real *X, int N, int nx, int nu,
                             cs_real T, cs_real *U)
{
    const cs_real dt = T / (cs_real)N;
    int k, j;
    for (k = 0; k < N; ++k)
        for (j = 0; j < nu; ++j)
            U[(size_t)k * nu + j] =
                (X[(size_t)(k + 1) * nx + (nx - nu + j)]
                 - X[(size_t)k * nx + (nx - nu + j)]) / dt;
}

int cs_rh_init(cs_rh *rh, cs_solver *s, const cs_real *S, cs_real T_seed,
               int N, int k_iters)
{
    int sN = 0, snx = 0, snu = 0, i, rc;
    if (!rh || !s || !S || N < 3 || N > CS_RH_NMAX || k_iters < 1
        || !(T_seed > (cs_real)0))
        return CS_ERR_ARG;
    rc = cs_solver_get_dims(s, &sN, &snx, &snu, (int *)0, (int *)0);
    if (rc != CS_OK)
        return rc;
    if (sN != N || snx != 6 || snu != 2)
        return CS_ERR_ARG;
    rh->s = s;
    rh->N = N;
    rh->nx = snx;
    rh->nu = snu;
    rh->k_iters = k_iters;
    rh->d_off = (cs_real)0;
    cs_copy(S, rh->S, (N + 1) * snx);
    rh->T_seed = T_seed;
    rh->engage_xi = S[2];
    cs_copy(S + (size_t)N * snx, rh->xN, snx);
    /* accepted plan = the seed with finite-difference controls (the
     * harness's pre-first-accept fallback) */
    cs_copy(S, rh->Xa, (N + 1) * snx);
    default_controls(S, N, snx, snu, T_seed, rh->Ua);
    rh->Ta = T_seed;
    rh->age = (cs_real)0;
    rh->age_valid = 0;
    rh->ref_w = (cs_real)0;              /* R5-proximity off = canonical rh */
    rh->ref_mode = 1;
    rh->val_slack_max = VAL_SLACK_MAX;
    rh->crop_bound = 0;                  /* R5-anchored crop floor off       */
    rh->crop_margin = (cs_real)0;
    rh->tau0 = (cs_real)0;
    rh->c_deep = rh->c_stall = 0;
    rh->prev_trem = (cs_real)0;
    rh->prev_trem_valid = 0;
    rh->recoveries = 0;
    rh->cool_until = (cs_real)-1e9;
    rh->sp_deficit = rh->sp_tau = rh->sp_dist = (cs_real)0;
    rh->t_now = (cs_real)0;
    rh->n_cycles = rh->n_fallbacks = 0;
    rh->last_reject = 0;
    rh->n_events = 0;
    for (i = 0; i < CS_RH_MAX_EVENTS; ++i) {
        rh->events[i].t = rh->events[i].tau = (cs_real)0;
        rh->events[i].deficit = rh->events[i].v = (cs_real)0;
        rh->events[i].dist = (cs_real)0;
        rh->events[i].rule = rh->events[i].action = 0;
    }
    /* arm the mission on the solver (fly_closed_loop's csolver setup +
     * rh_gate_reset): pins, elastic corridor, progress gate at the margin,
     * seed as iterate, canonical row offset. */
    rc = cs_solver_set_row_offset(s, (cs_real)0);
    if (rc == CS_OK)
        rc = cs_solver_set_pins(s, S, S + (size_t)N * snx);
    if (rc == CS_OK)
        rc = cs_solver_set_soft_corr(s, 1);
    if (rc == CS_OK)
        rc = cs_solver_set_gate(s, (cs_real)0, CS_RH_MARGIN, CS_RH_MARGIN,
                                CS_RH_MARGIN);
    if (rc == CS_OK)
        rc = cs_solver_set_seed(s, S, rh->Ua, &T_seed, 0);
    /* Clear any R5-proximity reference left on the SOLVER by a previous turn's
     * cs_rh_step. set_seed resets the iterate + QP warm start but NOT
     * qp->w_ref/x_ref, and the caller's preconvergence iterates the solver
     * DIRECTLY (bypassing cs_rh_step), so a stale reference from the prior turn
     * would pull the fresh preconvergence off the new seed and stall it (the
     * back-to-back second-turn collapse). rh_step re-installs the correct,
     * progress-restretched reference once armed (ref_w > 0). */
    if (rc == CS_OK)
        rc = cs_solver_set_ref(s, (const cs_real *)0, (cs_real)0, 0);
    /* likewise the crop floor: it reads the same qp->x_ref, so a stale
     * solver-side arm would bound the fresh preconvergence against the
     * PREVIOUS turn's reference. rh_step re-arms it per cycle. */
    if (rc == CS_OK)
        rc = cs_solver_set_crop_bound(s, 0, (cs_real)0);
    return rc;
}

int cs_rh_set_frame_d(cs_rh *rh, cs_real d)
{
    int rc;
    if (!rh || !rh->s || !(d == d))
        return CS_ERR_ARG;
    rc = cs_solver_set_row_offset(rh->s, d - CS_RH_D_CANON);
    if (rc != CS_OK)
        return rc;
    rh->d_off = d - CS_RH_D_CANON;
    return CS_OK;
}

int cs_rh_set_ref_tracking(cs_rh *rh, cs_real w, int mode, cs_real slack_max)
{
    if (!rh || w < (cs_real)0 || mode < 0 || mode > 3)
        return CS_ERR_ARG;
    rh->ref_w = w;
    rh->ref_mode = mode;
    if (slack_max > (cs_real)0)
        rh->val_slack_max = slack_max;
    /* Install this turn's reference NOW. The firmware preconverges through
     * raw cs_solver_iterate() calls before the first cs_rh_step(), so merely
     * storing ref_w leaves that phase reference-free (or, before the init
     * clear above, exposed to the previous turn's restretched x_ref). Keep
     * the objective continuous from ARMED through ENGAGED by arming the full
     * seed at tau0=0, or the appropriately restretched seed if this setter is
     * used later in a turn. cs_rh_step keeps it fresh thereafter. */
    if (w > (cs_real)0 && mode != 0) {
        restretch(rh->S, rh->Ua, rh->N, rh->nx, rh->nu, rh->tau0,
                  rh->Xref, rh->Us);
        return cs_solver_set_ref(rh->s, rh->Xref, w, mode);
    }
    return cs_solver_set_ref(rh->s, (const cs_real *)0, (cs_real)0, 0);
}

int cs_rh_set_crop_bound(cs_rh *rh, int on, cs_real margin)
{
    if (!rh || margin < (cs_real)0)
        return CS_ERR_ARG;
    rh->crop_bound = on ? 1 : 0;
    rh->crop_margin = margin;
    return CS_OK;
}

int cs_rh_arm_ready_ref(cs_rh *rh, const cs_real x_meas[6])
{
    cs_real tau_v, dist;
    if (!rh || !rh->s || !x_meas)
        return CS_ERR_ARG;
    if (rh->ref_w <= (cs_real)0 || rh->ref_mode == 0)
        return CS_OK;
    /* the vehicle's absolute progress along the seed path, ALLOWING
     * pre-path states: the READY tracking phase pins x0 up to 12 m before
     * the seed's first node, where a tau0=0-anchored reference is node-
     * misaligned — the per-index pull then drags every plan node toward a
     * further-along, SLOWER reference state and the tracked iterate learns
     * a premature braking profile (flown: both turns at ~3 m/s six metres
     * before the corner; doc/ready_ref_misalignment.pdf). The seed's entry
     * segment is the constant-cruise approach straight, so the backward
     * extension by first-segment extrapolation (negative tau through the
     * extrapolating interp_node/interp_unode) is exact: reference node 0
     * sits AT the vehicle, node-aligned with the tracked plan's span. */
    if (x_meas[2] >= rh->S[2]) {
        project_xy(rh->S, rh->nx, rh->N, x_meas[2], x_meas[3],
                   &tau_v, &dist);
    } else {
        const cs_real dxi = rh->S[(size_t)rh->nx + 2] - rh->S[2];
        if (!(dxi > (cs_real)1e-6) || !(x_meas[2] == x_meas[2]))
            return CS_ERR_ARG;
        tau_v = (x_meas[2] - rh->S[2]) / ((cs_real)rh->N * dxi);
    }
    if (tau_v < (cs_real)-0.5)
        tau_v = (cs_real)-0.5;
    if (tau_v > (cs_real)0.99)
        tau_v = (cs_real)0.99;
    restretch(rh->S, rh->Ua, rh->N, rh->nx, rh->nu, tau_v,
              rh->Xref, rh->Us);
    /* pre-path nodes (tq < 0, clamped to S[0] by the interp): hold the
     * node-0 state (cruise on the row — the seed BRAKES from node 0, so
     * extrapolating its first segment would inflate v backward) and extend
     * only xi linearly with the first-segment slope; at tq = tau_v this
     * lands exactly on the vehicle. */
    if (tau_v < (cs_real)0) {
        const cs_real dxi0 = (cs_real)rh->N
            * (rh->S[(size_t)rh->nx + 2] - rh->S[2]);
        int k;
        for (k = 0; k <= rh->N; ++k) {
            const cs_real tq = tau_v
                + ((cs_real)1 - tau_v) * (cs_real)k / (cs_real)rh->N;
            if (tq >= (cs_real)0)
                break;
            rh->Xref[(size_t)k * rh->nx + 2] = rh->S[2] + tq * dxi0;
        }
    }
    /* Align the WHOLE tracked problem, not just the reference: the
     * progress gates are index-anchored too (stage k's original progress
     * tk = tau0 + (1-tau0)k/N assumes node 0 = the seed start), so with
     * the pin up to 12 m earlier the corridor legs are enforced ~2-3
     * nodes too soon and the tracked iterate gets squeezed into an early-
     * braking shape REGARDLESS of the reference (host A/B: approach speed
     * collapses to ~5-8 m/s with the pull on or off; the engage chaining
     * then inherits that shape). The gate API's batch sentinel occupies
     * tau0 < 0, so a negative-progress gate cannot be expressed; instead
     * relax the leg gates to inert for the READY warm start — it is never
     * flown, and the aligned w-pull holds it on the seed. The first
     * engaged cs_rh_step restores real gating (it calls set_gate every
     * cycle). */
    {
        const int rc = cs_solver_set_gate(rh->s, (cs_real)0,
                                          (cs_real)CS_RH_READY_RLX,
                                          (cs_real)CS_RH_READY_RLX,
                                          (cs_real)CS_RH_READY_RLX);
        if (rc != CS_OK)
            return rc;
    }
    return cs_solver_set_ref(rh->s, rh->Xref, rh->ref_w, rh->ref_mode);
}

cs_real cs_rh_engage_xi(const cs_rh *rh)
{
    return rh ? rh->engage_xi : (cs_real)0;
}

cs_real cs_rh_progress(const cs_rh *rh)
{
    return rh ? rh->tau0 : (cs_real)0;
}

int cs_rh_last_reject(const cs_rh *rh)
{
    return rh ? rh->last_reject : 0;
}

size_t cs_rh_sizeof(void)
{
    return sizeof(cs_rh);
}

int cs_rh_events(const cs_rh *rh, cs_rh_event *out, int max_events, int *n)
{
    int i, m;
    if (!rh || !out || !n || max_events < 0)
        return CS_ERR_ARG;
    m = rh->n_events < max_events ? rh->n_events : max_events;
    for (i = 0; i < m; ++i)
        out[i] = rh->events[i];
    *n = m;
    return CS_OK;
}

/* SpiralDetector.update: feed one TURN cycle's measured state + the current
 * accepted plan's T_rem; returns the fired rule (0 = none), tau via out. */
static int spiral_update(cs_rh *rh, cs_real v, cs_real xi, cs_real eta,
                         cs_real trem, cs_real *tau_out)
{
    cs_real tau, dist, vref;
    int in_win, hit, rule = 0;
    project_xy(rh->S, rh->nx, rh->N, xi, eta, &tau, &dist);
    vref = interp_node(rh->S, rh->nx, rh->N, tau);
    rh->sp_deficit = vref - v;
    rh->sp_tau = tau;
    rh->sp_dist = dist;
    in_win = (tau > SP_TAU_LO && tau < SP_TAU_HI);
    rh->c_deep = (in_win && rh->sp_deficit > SP_DEEP_DFC)
        ? rh->c_deep + 1 : 0;
    hit = (rh->prev_trem_valid && trem == trem
           && rh->prev_trem - trem < SP_STALL_DT
           && trem <= SP_STALL_TREM && tau >= SP_STALL_TAU);
    rh->c_stall = hit ? rh->c_stall + 1 : 0;
    if (trem == trem) {
        rh->prev_trem = trem;
        rh->prev_trem_valid = 1;
    }
    if (rh->t_now >= rh->cool_until && rh->recoveries < SP_MAX_RECOV) {
        if (rh->c_deep >= SP_DEEP_N)
            rule = CS_RH_RULE_DEEP;
        else if (rh->c_stall >= SP_STALL_N)
            rule = CS_RH_RULE_TSTALL;
    }
    if (!rule)
        return 0;
    ++rh->recoveries;
    rh->c_deep = rh->c_stall = 0;
    rh->cool_until = rh->t_now + SP_COOLDOWN;
    if (rh->n_events < CS_RH_MAX_EVENTS) {
        cs_rh_event *e = &rh->events[rh->n_events++];
        e->t = rh->t_now;
        e->rule = rule;
        e->action = 0;
        e->tau = tau;
        e->deficit = rh->sp_deficit;
        e->v = v;
        e->dist = dist;
    }
    *tau_out = tau;
    return rule;
}

int cs_rh_step(cs_rh *rh, const cs_real *x_meas, cs_real dt,
               cs_real *X_plan, cs_real *U_plan, cs_real *T_plan,
               cs_real *T_rem, int *accepted, int *event)
{
    const int N = rh ? rh->N : 0, nx = rh ? rh->nx : 0,
        nu = rh ? rh->nu : 0;
    cs_real trem_now, tau_star, dproj, Tit, sl;
    cs_real h[CS_PD_NH];
    const cs_real u0[CS_PD_NU] = {(cs_real)0, (cs_real)0};
    cs_real tau0_prev, v_out, v_ret, vx;
    int rc, i, k, rule, recov = 0, k_eff, k0, mask, qp_status, ev;
    cs_real tau_r = (cs_real)0;
    if (!rh || !rh->s || !x_meas)
        return CS_ERR_ARG;
    ++rh->n_cycles;
    rh->t_now += dt;
    if (rh->age_valid)
        rh->age += dt;
    trem_now = rh->age_valid ? rh->Ta - rh->age : (cs_real)NAN;
    ev = CS_RH_EVENT_NONE;
    rh->last_reject = 0;

    /* ---- 1. spiral detection BEFORE the solve (CAND-2) ---- */
    rule = spiral_update(rh, x_meas[0], x_meas[2], x_meas[3], trem_now,
                         &tau_r);
    if (rule) {
        if (tau_r > FE_TAU
            || x_meas[2] < rh->S[(size_t)N * nx + 2] + FE_XI_PAD) {
            /* LATE trigger: at/past the window exit — a tau~1 reseed would
             * place targets BEHIND the vehicle. No solve; the caller hands
             * back to stock guidance. */
            rh->events[rh->n_events - 1].action = CS_RH_EVENT_EXIT_FORCED;
            if (X_plan)
                cs_copy(rh->Xa, X_plan, (N + 1) * nx);
            if (U_plan)
                cs_copy(rh->Ua, U_plan, N * nu);
            if (T_plan)
                *T_plan = rh->Ta;
            if (!rh->age_valid) {
                rh->age = (cs_real)0;
                rh->age_valid = 1;
            }
            if (T_rem)
                *T_rem = rh->Ta - rh->age;
            if (accepted)
                *accepted = 0;
            if (event)
                *event = CS_RH_EVENT_EXIT_FORCED;
            return CS_OK;
        }
        /* DEEP/TSTALL reseed: solve from the reference-family seed
         * restretched to the vehicle's progress (the seed-policy lesson:
         * the solver only converges from R5-family seeds). The reseed
         * becomes the accepted plan NOW so a rejected recovery solve falls
         * back to IT, not to the degraded pre-recovery plan. */
        rh->events[rh->n_events - 1].action = CS_RH_EVENT_RESEED;
        restretch(rh->S, rh->Ua, N, nx, nu, tau_r, rh->Xs, rh->Us);
        rh->Ta = rh->T_seed * ((cs_real)1 - tau_r);
        default_controls(rh->Xs, N, nx, nu, rh->Ta, rh->Us);
        cs_copy(rh->Xs, rh->Xa, (N + 1) * nx);
        cs_copy(rh->Us, rh->Ua, N * nu);
        cs_copy(rh->Xa + (size_t)N * nx, rh->xN, nx);
        rc = cs_solver_set_seed(rh->s, rh->Xa, rh->Ua, &rh->Ta, 0);
        if (rc != CS_OK)
            return rc;
        rh->tau0 = tau_r;
        rh->age = (cs_real)0;
        rh->age_valid = 1;
        recov = 1;
        ev = CS_RH_EVENT_RESEED;
    }
    tau0_prev = rh->tau0;                  /* rollback point (audit B3) */

    /* ---- 2. innovation-based iteration escalation (M1.4/B2) ---- */
    rc = cs_solver_get_iterate(rh->s, rh->Xi_, rh->Ui_, &Tit);
    if (rc != CS_OK)
        return rc;
    k0 = 0;
    {
        cs_real best = (cs_real)0;
        for (k = 0; k <= N; ++k) {
            const cs_real dx = rh->Xi_[(size_t)k * nx + 2] - x_meas[2];
            const cs_real dy = rh->Xi_[(size_t)k * nx + 3] - x_meas[3];
            const cs_real d2 = dx * dx + dy * dy;
            if (k == 0 || d2 < best) {
                best = d2;
                k0 = k;
            }
        }
    }
    k_eff = rh->k_iters;
    {
        const cs_real dx = x_meas[2] - rh->Xi_[(size_t)k0 * nx + 2];
        const cs_real dy = x_meas[3] - rh->Xi_[(size_t)k0 * nx + 3];
        const cs_real dv = x_meas[0] - rh->Xi_[(size_t)k0 * nx + 0];
        if ((cs_real)sqrt((double)(dx * dx + dy * dy)) > INNOV_GATE
            || (dv < (cs_real)0 ? -dv : dv) > INNOV_GATE)
            k_eff = rh->k_iters > K_INNOV ? rh->k_iters : K_INNOV;
    }
    if (recov)
        k_eff = k_eff > K_RECOV ? k_eff : K_RECOV;

    /* ---- 3. shrinking-horizon restretch (rh_step, fractional) ----
     * ABSOLUTE progress along the immutable SEED path, not the iterate's
     * own polyline. Iterate-anchored progress was self-referential: a plan
     * that drifted one node early along the path measured the vehicle one
     * node further, re-licensed itself through the restretch reseed, and
     * the loop held a self-consistent solution exactly one node spacing
     * (~1.5 m) field-side — on EVERY non-first turn (the 3-turn fork:
     * first-turn approaches are pristine, post-handback approaches carry
     * the small eta settling transient that kicks the loop into the
     * shifted basin; neither the w=10 reference pull nor its removal
     * changed it, because a node-slide along the same path is near-tangent
     * to the pull's metric). The seed path is restretch-invariant and is
     * what the spiral detector and the Xref pull already anchor to, so
     * tau_abs is true row progress and the shifted equilibrium stops
     * being self-consistent. The iterate restretch keeps its incremental
     * form: tau_star = (tau_abs - tau0)/(1 - tau0); jitter behind
     * (tau_star <= 0) skips the advance — monotonicity as before. */
    {
        cs_real tau_abs = (cs_real)0;
        project_xy(rh->S, nx, N, x_meas[2], x_meas[3], &tau_abs, &dproj);
        tau_star = (rh->tau0 < (cs_real)0.999)
            ? (tau_abs - rh->tau0) / ((cs_real)1 - rh->tau0)
            : (cs_real)0;
    }
    if (tau_star > (cs_real)0
        && tau_star < (cs_real)(N - 2) / (cs_real)N) {
        cs_real Tn = Tit * ((cs_real)1 - tau_star);
        rh->tau0 += ((cs_real)1 - rh->tau0) * tau_star;
        restretch(rh->Xi_, rh->Ui_, N, nx, nu, tau_star, rh->Xs, rh->Us);
        rc = cs_solver_set_seed(rh->s, rh->Xs, rh->Us, &Tn, 0);
        if (rc != CS_OK)
            return rc;
    }

    /* ---- 4. leg-membership-gated relaxations (M1.3), via cs_h with the
     * row offset applied to the return row ---- */
    rc = cs_h(x_meas, u0, h);
    if (rc != CS_OK)
        return rc;
    v_out = h[2] > (cs_real)0 ? h[2] : (cs_real)0;
    v_ret = rh->d_off - h[3] > (cs_real)0 ? rh->d_off - h[3] : (cs_real)0;
    if (rh->tau0 > (cs_real)CS_PD_FRAC_OUT)
        v_out = (cs_real)0;              /* vehicle past the outbound leg */
    if (rh->tau0 < (cs_real)CS_PD_FRAC_RET)
        v_ret = (cs_real)0;              /* not yet on the return leg     */
    vx = (cs_real)0;
    if (rh->tau0 > (cs_real)CS_PD_FRAC_OUT
        && rh->tau0 < (cs_real)CS_PD_FRAC_RET) {
        vx = (cs_real)CS_PD_XI_CROSS_MIN - x_meas[2];
        if (vx < (cs_real)0)
            vx = (cs_real)0;
    }
    rc = cs_solver_set_gate(rh->s, rh->tau0, v_out + CS_RH_MARGIN,
                            v_ret + CS_RH_MARGIN, vx + CS_RH_MARGIN);
    if (rc != CS_OK)
        return rc;

    /* ---- 4b. R5-proximity reference (turn-2 over-bulge fix) ---- */
    /* restretch the init seed to the current progress tau0 (progress-aligned
     * with the plan nodes, like the Python rh_step) and apply as the L2
     * tracking reference; rh->Us is free scratch after the warm-start
     * set_seed above. Off (ref_w=0) leaves the solver's ref weight at 0. */
    if (rh->ref_w > (cs_real)0 || rh->crop_bound) {
        /* the crop floor reads x_ref too, so restretch/install whenever
         * EITHER the proximity pull (ref_w>0) OR the crop floor is armed;
         * ref_w may be 0 (pure crop floor, no bulge push). */
        restretch(rh->S, rh->Ua, N, nx, nu, rh->tau0, rh->Xref, rh->Us);
        rc = cs_solver_set_ref(rh->s, rh->Xref, rh->ref_w, rh->ref_mode);
        if (rc != CS_OK)
            return rc;
    }
    rc = cs_solver_set_crop_bound(rh->s, rh->crop_bound, rh->crop_margin);
    if (rc != CS_OK)
        return rc;

    /* ---- 5. pin the measured state, run the RTI iterations ---- */
    rc = cs_solver_set_pins(rh->s, x_meas, rh->xN);
    if (rc != CS_OK)
        return rc;
    for (i = 0; i < k_eff; ++i) {
        rc = cs_solver_iterate(rh->s);
        if (rc != CS_OK)
            return rc;
    }

    /* ---- 6. plan validation (the harness's reject terms) ---- */
    rc = cs_solver_get_iterate(rh->s, rh->Xi_, rh->Ui_, &Tit);
    if (rc != CS_OK)
        return rc;
    rc = cs_solver_get_info(rh->s, (cs_real *)0, (cs_real *)0, (cs_real *)0,
                            (cs_real *)0, &qp_status, (int *)0);
    if (rc != CS_OK)
        return rc;
    sl = (cs_real)0;
    rc = cs_solver_get_soft_slack(rh->s, &sl);
    if (rc != CS_OK)
        return rc;
    mask = 0;
    {
        int finite = 1;
        cs_real xi_max = rh->Xi_[2], eta_max = rh->Xi_[3],
            v_min = rh->Xi_[0];
        for (k = 0; k <= N && finite; ++k)
            for (i = 0; i < nx; ++i)
                if (!isfinite((double)rh->Xi_[(size_t)k * nx + i])) {
                    finite = 0;
                    break;
                }
        for (k = 0; k < N && finite; ++k)
            for (i = 0; i < nu; ++i)
                if (!isfinite((double)rh->Ui_[(size_t)k * nu + i])) {
                    finite = 0;
                    break;
                }
        if (!finite)
            mask |= CS_RH_REJ_NONFINITE;
        else {
            for (k = 1; k <= N; ++k) {
                const cs_real xi_k = rh->Xi_[(size_t)k * nx + 2];
                const cs_real eta_k = rh->Xi_[(size_t)k * nx + 3];
                const cs_real v_k = rh->Xi_[(size_t)k * nx + 0];
                if (xi_k > xi_max)
                    xi_max = xi_k;
                if (eta_k > eta_max)
                    eta_max = eta_k;
                if (v_k < v_min)
                    v_min = v_k;
            }
            if (!(Tit > VAL_T_LO && Tit < VAL_T_HI))
                mask |= CS_RH_REJ_T;
            if (!(xi_max < VAL_XI_MAX))
                mask |= CS_RH_REJ_XI;
            if (!(eta_max < VAL_ETA_MAX + rh->d_off))
                mask |= CS_RH_REJ_ETA;
            if (!(v_min > (cs_real)0))
                mask |= CS_RH_REJ_VMIN;
        }
        if (!(qp_status > 0))
            mask |= CS_RH_REJ_QP;
        if (sl > rh->val_slack_max)
            mask |= CS_RH_REJ_SLACK;
    }
    rh->last_reject = mask;

    /* ---- 7. accept / reject with iterate restore + tau0 rollback ---- */
    if (mask == 0) {
        cs_copy(rh->Xi_, rh->Xa, (N + 1) * nx);
        cs_copy(rh->Ui_, rh->Ua, N * nu);
        rh->Ta = Tit;
        rh->age = (cs_real)0;
        rh->age_valid = 1;
    } else {
        /* anytime fallback: restore the accepted plan as the iterate
         * (rh_reset semantics: the terminal pin follows the accepted
         * plan's last node) and roll the progress advance back. */
        ++rh->n_fallbacks;
        rc = cs_solver_set_seed(rh->s, rh->Xa, rh->Ua, &rh->Ta, 0);
        if (rc != CS_OK)
            return rc;
        cs_copy(rh->Xa + (size_t)N * nx, rh->xN, nx);
        rc = cs_solver_set_pins(rh->s, x_meas, rh->xN);
        if (rc != CS_OK)
            return rc;
        rh->tau0 = tau0_prev;
    }
    if (!rh->age_valid) {         /* first cycle rejected: age the seed */
        rh->age = (cs_real)0;
        rh->age_valid = 1;
    }

    if (X_plan)
        cs_copy(rh->Xa, X_plan, (N + 1) * nx);
    if (U_plan)
        cs_copy(rh->Ua, U_plan, N * nu);
    if (T_plan)
        *T_plan = rh->Ta;
    if (T_rem)
        *T_rem = rh->Ta - rh->age;
    if (accepted)
        *accepted = (mask == 0) ? 1 : 0;
    if (event)
        *event = ev;
    return CS_OK;
}
