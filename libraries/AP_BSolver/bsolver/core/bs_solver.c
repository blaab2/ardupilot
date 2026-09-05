/* bsolver core — see bs_solver.h.
 *
 * Correspondence with the certificates, which is the whole reason this file
 * exists (MixedHorizon in wp5_selfconcordance.py / wp5_interleaved_certificate.py):
 *
 *   Phi[0] = I,  Phi[t+1] = T(k+t+1) A Phi[t]
 *   Abar[6t:6t+6]              = Phi[t+1]
 *   Gam[6t:6t+6, 3t:3t+3]      = T(k+t+1) B
 *   Gam[6t:6t+6, :3t]          = T(k+t+1) A Gam[6(t-1):6t, :3t]
 *   Wt = kron(I_N, Q) with the last block replaced by P[fam[N]]
 *   J  = xi'Q xi + U'Rbar U + Xf'Wt Xf
 *        + sum_t eps * [ b(z) - b(m) + b'(m)(m - z) ],   z = m - g
 *   g  = 2(Rbar U + Gam' Wt Xf) + sum_t Rmat[t]' ( eps (b'(m) - b'(z)) )
 *   H  = 2(Rbar + Gam' Wt Gam) + sum_t Rmat[t]' diag( eps / max(z,DR)^2 ) Rmat[t]
 *
 * with b(z) = -log z for z >= DR and the quadratic extension below it, and
 * Rmat[t] = Hr[t] Gam[6(t-1):6t] (zero for t = 0) with Er[t] added into the
 * t-th input block.  The Hessian barrier weight uses the CLIPPED slack
 * max(z, DR) — the extension's curvature — exactly as the reference does.
 */
#include "bs_solver.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

#define NX BS_NX
#define NU BS_NU
#define NN BS_N
#define NV BS_NV
#define NR BS_NROW
#define NXT (BS_N * BS_NX)

#ifndef BS_TABLES_RUNTIME
#ifdef BS_N_PHASE
/* The 59-phase periodic parity clock lives only in the record's data header
 * -- it is an "a"/"t" object and the corner-online header retires family "a".
 * Guarded, NOT deleted: the record build is still the certified one and B1
 * is stated on this clock. */
const bs_schedule bs_sched_periodic = {
    bs_periodic_family, bs_periodic_rot, NULL, BS_N_PHASE, 1,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL
};
#endif
const bs_schedule bs_sched_mission = {
    bs_mission_family, bs_mission_rot,
#ifdef BS_NOFF
    bs_mission_off,
#else
    NULL,                      /* the record's mission has no affine part */
#endif
    BS_N_MISSION, 0,
    NULL, NULL, NULL,          /* no RAM-extra tables on a flash schedule */
    NULL, NULL, NULL, NULL, NULL
};
#endif /* !BS_TABLES_RUNTIME */

/* The affine part of a schedule that has none.  off_of() hands this back for
 * schedule->offset == NULL, so every plant recursion is written once. */
static const bs_real bs_off_zero[NX] = { 0.0 };

static int tick_of(const bs_schedule *schedule, int t)
{
    if (schedule->periodic) {
        int wrapped = t % schedule->length;
        if (wrapped < 0) wrapped += schedule->length;
        return wrapped;
    }
    if (t < 0) t = 0;
    if (t >= schedule->length) t = schedule->length - 1;
    if (schedule->ring_mask) t &= schedule->ring_mask;
    return t;
}

static const bs_real *rot_of(const bs_schedule *schedule, int t)
{
    const int idx = schedule->rotation[tick_of(schedule, t)];
    if (idx < 0) return schedule->rot_extra;      /* runtime seam (ingress) */
    if (schedule->rot_tab) return &schedule->rot_tab[(size_t)idx * NX * NX];
#ifndef BS_TABLES_RUNTIME
    return &bs_rot[(size_t)idx * NX * NX];
#else
    return schedule->rot_extra;   /* unreachable by construction */
#endif
}

/* The affine seam offset applied at DESTINATION tick t, i.e. the c(t) of
 *     xi_t = T(t) (A xi_{t-1} + B u_{t-1}) + c(t).
 * Mirrors rot_of(): same clamped/wrapped tick, same table-plus-index
 * encoding, slot 0 the zero vector as slot 0 of bs_rot is the identity. */
static const bs_real *off_of(const bs_schedule *schedule, int t)
{
    if (schedule->offset == NULL) return bs_off_zero;
    {
        const int idx = schedule->offset[tick_of(schedule, t)];
        if (idx < 0) return schedule->off_extra;  /* runtime seam (ingress) */
        if (schedule->off_tab) return &schedule->off_tab[(size_t)idx * NX];
#ifdef BS_NOFF
        return &bs_off[(size_t)idx * NX];
#else
        return bs_off_zero;                       /* no flash table to index */
#endif
    }
}

static int family_of(const bs_schedule *schedule, int t)
{
    return schedule->family[tick_of(schedule, t)];
}

static const bs_real *rows_of(const bs_schedule *schedule, int t)
{
    const int family = schedule->family[tick_of(schedule, t)];
    if (family < 0) return schedule->rows_extra;  /* runtime family (ingress) */
    if (schedule->rows_tab)
        return &schedule->rows_tab[(size_t)family * NR * 10];
#ifndef BS_TABLES_RUNTIME
    return &bs_rows[(size_t)family * NR * 10];
#else
    return schedule->rows_extra;  /* unreachable by construction */
#endif
}

/* Terminal-law family lookup: bs_P / bs_K have no RAM-extra entry, and by
 * the ingress driver's construction (prefix capped at BS_N ticks) a negative
 * family never reaches a horizon-terminal tick.  Clamp defensively to the
 * trim family rather than index flash with a negative — wrong-but-bounded
 * beats out-of-bounds if the invariant is ever broken. */
static int terminal_family_of(const bs_schedule *schedule, int t)
{
    const int family = family_of(schedule, t);
    return (family < 0) ? 0 : family;
}

/* c = a * b, (m x k) * (k x n), all row-major. */
static void mat_mul(const bs_real *a, const bs_real *b, bs_real *c,
                    int m, int k, int n)
{
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) c[(size_t)i * n + j] = 0.0;
        for (int p = 0; p < k; ++p) {
            const bs_real aip = a[(size_t)i * k + p];
            if (aip == 0.0) continue;
            const bs_real *brow = &b[(size_t)p * n];
            bs_real *crow = &c[(size_t)i * n];
            for (int j = 0; j < n; ++j) crow[j] += aip * brow[j];
        }
    }
}

size_t bs_problem_size(void)
{
    return sizeof(bs_problem);
}

size_t bs_workspace_size(void)
{
    /* Gam is NOT stored.  Its block rows obey the one-step recursion
     *   Gam_t = [ T(k+t+1) A Gam_(t-1)(:, <3t) | T(k+t+1) B | 0 ]
     * and every consumer walks t in increasing order, so a ping-pong pair of
     * NX x NV block rows suffices.  Storing it costs NXT x NV = 230,400 B on
     * an f64 target, which is 55% of the working set and the difference
     * between fitting an STM32H743 and missing it by 29,664 B.  The predicted
     * states are taken from the plant recursion instead of from Gam, so no
     * consumer needs a block row out of order. */
    return (size_t)NV * NV           /* Hquad                         */
         + (size_t)NV * NX           /* Cquad                         */
         + (size_t)NV                /* Dquad, the affine drift       */
         + (size_t)NR * NV           /* Rmat scratch (eval)           */
         + (size_t)NXT               /* predicted states (eval)       */
         + (size_t)NR                /* deferred barrier weights      */
         + 2 * (size_t)NX * NV       /* Gam block row, ping-pong      */
         + (size_t)NX * NV           /* W*Gam block row (init)        */
         + (size_t)NX * NX           /* W*Abar block (init)           */
         + (size_t)NV * NV           /* Newton Hessian copy           */
         + 3 * (size_t)NV;           /* grad, direction, trial        */
}

/* Build the phase-dependent objects, accumulating the quadratic half only on
 * the rows a >= NU*npin.
 *
 * WHY THIS IS EXACT, not an approximation.  Hquad and Cquad are read in
 * exactly two places (eval_impl): the gradient loop, which visits rows
 * i >= NU*npin, and the Hessian copy, which visits rows a >= NU*npin.  The
 * OBJECTIVE VALUE does not touch them at all -- it is summed from the plant
 * recursion's predicted states and the stage costs directly.  So on the
 * pinned path the rows below NU*npin are dead storage, and not writing them
 * changes no computed quantity.  Measured, not argued: tail gradient and
 * tail Hessian against the fully built problem, 0.000e+00.
 *
 * THE CONTRACT IT CHANGES, stated so nobody trips over it: a problem built
 * with npin = p may only be evaluated at npin >= p.  eval_impl enforces it
 * rather than trusting it -- reading an unbuilt row would silently return
 * zeros, which is the kind of defect that passes every smoke test.
 *
 * bs_problem_init() is this at npin = 0 and is bit-identical to the version
 * that had no npin at all.
 *
 * THE AFFINE DRIFT, Dquad.  With affine seams the predicted states carry a
 * constant Cbar_t = T A Cbar_(t-1) + c(t) on top of Abar xi + Gam U, so the
 * reference gradient
 *     g = 2 ( Rbar U + Gam' Wt (Abar xi + Gam U + Cbar) )
 * gains a third piece, Dquad = sum_t 2 Gam_t' W_t Cbar_t, which multiplies
 * neither xi nor U and therefore cannot hide inside Hquad or Cquad.  It is
 * accumulated here, in the same stage loop and under the SAME npin
 * restriction (rows a >= NU*npin only), because it is read in exactly one
 * place: eval_impl's gradient loop, which visits those rows.  The HESSIAN is
 * untouched -- a constant added to the states shifts the gradient, not the
 * curvature.  On a schedule without offsets Cbar is identically zero and
 * Dquad is a dead vector of zeros. */
bs_status bs_problem_init_pinned(bs_problem *problem,
                                 const bs_schedule *schedule, int phase,
                                 int npin, bs_real *memory, size_t memory_len)
{
    if (memory_len < bs_workspace_size()) return BS_ERR_ARENA;
    if (schedule == NULL || schedule->length <= 0) return BS_ERR_PHASE;
    if (npin < 0 || npin > NN) return BS_ERR_PHASE;
    const int qpin = NU * npin;     /* first quadratic row that is BUILT */

    bs_real *cursor = memory;
    problem->schedule = schedule;
    problem->phase = phase;
    problem->npin_built = npin;
    problem->Hquad = cursor; cursor += (size_t)NV * NV;
    problem->Cquad = cursor; cursor += (size_t)NV * NX;
    problem->Dquad = cursor; cursor += (size_t)NV;
    problem->gam = cursor; cursor += (size_t)NX * NV;
    problem->gam_next = cursor; cursor += (size_t)NX * NV;
    problem->nt_hess = cursor; cursor += (size_t)NV * NV;
    problem->nt_grad = cursor; cursor += (size_t)NV;
    problem->nt_dir = cursor; cursor += (size_t)NV;
    problem->nt_trial = cursor; cursor += (size_t)NV;
    problem->scratch = cursor;

    const bs_real *Pterm = schedule->P_tab
        ? &schedule->P_tab[(size_t)terminal_family_of(schedule, phase + NN)
                           * NX * NX]
#ifndef BS_TABLES_RUNTIME
        : &bs_P[(size_t)terminal_family_of(schedule, phase + NN) * NX * NX];
#else
        : schedule->P_tab;        /* unreachable by construction */
#endif
    bs_real *WG = problem->scratch;            /* NX x NV */
    bs_real *WA = WG + (size_t)NX * NV;        /* NX x NX */
    bs_real phi[NX * NX], next[NX * NX], TA[NX * NX], TB[NX * NU];
    /* Cbar_t, the free drift of the affine plant: the state recursion run at
     * xi = 0, U = 0.  Rolled forward alongside Phi and never stored. */
    bs_real cbar[NX], cnext[NX], Wc[NX];

    for (int i = 0; i < NX * NX; ++i) phi[i] = (i % (NX + 1) == 0) ? 1.0 : 0.0;
    for (int i = 0; i < NX; ++i) cbar[i] = 0.0;
    memset(&problem->Hquad[(size_t)qpin * NV], 0,
           sizeof(bs_real) * (size_t)(NV - qpin) * NV);
    memset(&problem->Cquad[(size_t)qpin * NX], 0,
           sizeof(bs_real) * (size_t)(NV - qpin) * NX);
    memset(&problem->Dquad[qpin], 0,
           sizeof(bs_real) * (size_t)(NV - qpin));
    memset(problem->gam, 0, sizeof(bs_real) * (size_t)NX * NV);

    for (int t = 0; t < NN; ++t) {
        const bs_real *T = rot_of(schedule, phase + t + 1);
        mat_mul(T, bs_A, TA, NX, NX, NX);
        mat_mul(T, bs_B, TB, NX, NX, NU);

        /* Free response block: Abar_t = T A Phi[t] = Phi[t+1].  Kept in the
         * local `next` and consumed below in this same iteration, so the
         * NXT x NX array it used to occupy is not needed. */
        mat_mul(TA, phi, next, NX, NX, NX);

        /* Drift block: Cbar_t = T(k+t+1) A Cbar_(t-1) + c(k+t+1), the same
         * recursion the predicted states run, seeded at zero.  Verbatim
         * MixedHorizonAff.__init__: `cprev = Tt @ A6 @ cprev` then
         * `+ off[t]`, so the T*A product is formed FIRST here too. */
        {
            const bs_real *c = off_of(schedule, phase + t + 1);
            for (int i = 0; i < NX; ++i) {
                bs_real acc = 0.0;
                for (int j = 0; j < NX; ++j)
                    acc += TA[(size_t)i * NX + j] * cbar[j];
                cnext[i] = acc + c[i];
            }
            memcpy(cbar, cnext, sizeof(cbar));
        }

        /* Gam block row t, from block row t-1 (ping-pong, never stored) */
        bs_real *here = (t == 0) ? problem->gam : problem->gam_next;
        if (t >= 1) {
            const bs_real *prev = problem->gam;
            for (int i = 0; i < NX; ++i) {
                bs_real *out = &here[(size_t)i * NV];
                for (int j = 0; j < NU * t; ++j) {
                    bs_real acc = 0.0;
                    for (int p = 0; p < NX; ++p)
                        acc += TA[(size_t)i * NX + p] * prev[(size_t)p * NV + j];
                    out[j] = acc;
                }
                for (int j = NU * t; j < NV; ++j) out[j] = 0.0;
            }
        }
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NU; ++j)
                here[(size_t)i * NV + (NU * t + j)] = TB[(size_t)i * NU + j];
        if (t >= 1) {
            bs_real *swap = problem->gam;
            problem->gam = problem->gam_next;
            problem->gam_next = swap;
        }

        /* Accumulate the quadratic half against this block row only:
         *   Hquad += 2 Gam_t' W_t Gam_t,   Cquad += 2 Gam_t' W_t Abar_t. */
        const bs_real *W = (t == NN - 1) ? Pterm : bs_Q;
        const bs_real *Gt = problem->gam;
        const bs_real *At = next;      /* = Abar_t, this iteration only */
        const int width = NU * (t + 1);      /* causality: zero beyond this */

        for (int i = 0; i < NX; ++i) {
            bs_real *wg = &WG[(size_t)i * NV];
            for (int j = 0; j < width; ++j) wg[j] = 0.0;
            for (int p = 0; p < NX; ++p) {
                const bs_real w = W[(size_t)i * NX + p];
                if (w == 0.0) continue;
                const bs_real *g = &Gt[(size_t)p * NV];
                for (int j = 0; j < width; ++j) wg[j] += w * g[j];
            }
            bs_real *wa = &WA[(size_t)i * NX];
            for (int j = 0; j < NX; ++j) {
                bs_real acc = 0.0;
                for (int p = 0; p < NX; ++p)
                    acc += W[(size_t)i * NX + p] * At[(size_t)p * NX + j];
                wa[j] = acc;
            }
            bs_real wc = 0.0;
            for (int p = 0; p < NX; ++p)
                wc += W[(size_t)i * NX + p] * cbar[p];
            Wc[i] = wc;
        }
        for (int p = 0; p < NX; ++p) {
            const bs_real *grow = &Gt[(size_t)p * NV];
            const bs_real *wg = &WG[(size_t)p * NV];
            const bs_real *wa = &WA[(size_t)p * NX];
            const bs_real wcp = Wc[p];
            for (int a = qpin; a < width; ++a) {
                const bs_real ga = grow[a];
                if (ga == 0.0) continue;
                const bs_real scaled = 2.0 * ga;
                bs_real *hrow = &problem->Hquad[(size_t)a * NV];
                for (int b = 0; b < width; ++b) hrow[b] += scaled * wg[b];
                bs_real *crow = &problem->Cquad[(size_t)a * NX];
                for (int b = 0; b < NX; ++b) crow[b] += scaled * wa[b];
                problem->Dquad[a] += scaled * wcp;
            }
        }
#ifdef BS_SPEED_TILT_W
        /* gradient of the speed tilt: d/dU of -w delta_t is
         * -w * (row 0 of Gam_t); a constant, so it lives in Dquad
         * beside the affine drift. */
        for (int a = qpin; a < width; ++a)
            problem->Dquad[a] -= (bs_real)BS_SPEED_TILT_W * Gt[a];
#endif
        memcpy(phi, next, sizeof(phi));    /* advance after Abar_t is used */
    }
    /* + 2 Rbar (block diagonal over the input stages), built rows only */
    for (int t = npin; t < NN; ++t)
        for (int i = 0; i < NU; ++i)
            for (int j = 0; j < NU; ++j)
                problem->Hquad[(size_t)(t * NU + i) * NV + (t * NU + j)]
                    += 2.0 * bs_R[(size_t)i * NU + j];

    return BS_OK;
}

bs_status bs_problem_init(bs_problem *problem, const bs_schedule *schedule,
                          int phase, bs_real *memory, size_t memory_len)
{
    return bs_problem_init_pinned(problem, schedule, phase, 0,
                                  memory, memory_len);
}

static bs_real barrier_value(bs_real z)
{
    if (z >= BS_DR) return -log(z);
    const bs_real u = (z - 2.0 * BS_DR) / BS_DR;
    return 0.5 * (u * u - 1.0) - log(BS_DR);
}

static bs_real barrier_slope(bs_real z)
{
    if (z >= BS_DR) return -1.0 / z;
    return (z - 2.0 * BS_DR) / (BS_DR * BS_DR);
}

/* Objective, gradient and Hessian with the first `npin` input stages treated
 * as PINNED DATA.  The value returned is the FULL-horizon objective; the
 * gradient and Hessian returned are the TAIL BLOCK of the full-horizon
 * objects at (U, xi),
 *
 *     grad[j - 3*npin] = g[j],           j  >= 3*npin
 *     hess[(a-p)*nf + (b-p)] = H[a][b],  a, b >= 3*npin,  p = 3*npin
 *
 * packed to nf = NV - 3*npin.  npin = 0 recovers bs_eval() on the same code
 * path, bit for bit.
 *
 * THIS IS NOT A TRUNCATION.  By causality (bsolver.tex Prop. 1) the row map
 * R_t is zero beyond column NU(t+1), so a stage t < npin cannot reach a
 * column >= 3*npin, and Gam_(npin-1) restricted to those columns is zero.
 * The tail block therefore receives contributions from stages t >= npin only,
 * and skipping the earlier stages reproduces it TERM FOR TERM AND IN THE SAME
 * ORDER — bit-identical to assembling the full objects and slicing, which
 * tests/test_parity_pinned.py checks directly (max |delta| exactly 0).
 *
 * That restriction IS the horizon-(N - npin) tail problem at phase k + npin
 * with boundary state x_npin, which is how wp5_h1_interleaved_q.py formulates
 * the interleaved certificate; the two agree by the boundary-state identity.
 *
 * The value is deliberately NOT restricted: the reference's guard (ii) and
 * Armijo test are stated on the full objective (newton_q_pinned in
 * wp5_anytime_sim.py calls self.J on the full plan), and the pinned prefix
 * contributes a constant to both sides of the Armijo inequality — equal in
 * exact arithmetic, not equal in the last bits.  Matching the reference costs
 * one full value sweep and is worth it.
 */
static bs_status eval_impl(const bs_problem *problem, const bs_real *U,
                           const bs_real *xi, int npin, bs_real *value,
                           bs_real *grad, bs_real *hess)
{
    const int phase = problem->phase;
    const bs_schedule *schedule = problem->schedule;
    bs_real *Rmat = problem->scratch;                       /* NR x NV */
    bs_real *states = problem->scratch + (size_t)NR * NV;   /* NXT      */
    bs_real *weight = states + (size_t)NXT;                 /* NR       */
    const int want_derivatives = (grad != NULL) || (hess != NULL);
    const int pin = NU * npin;               /* first FREE column          */
    const int nf = NV - pin;                 /* free-tail dimension        */
    const int t_start = npin;                /* first contributing stage   */
    const int t_adv = (npin > 0) ? npin : 1; /* first Gam advance          */

    if (npin < 0 || npin > NN) return BS_ERR_PHASE;
    /* A problem built with the quadratic half restricted to rows >= 3*p can
     * only answer evaluations at npin >= p; below that the rows this reads
     * were never written.  Refuse rather than return zeros. */
    if (want_derivatives && npin < problem->npin_built) return BS_ERR_PHASE;

    /* Predicted states by the PLANT recursion, not through Gam:
     *   x_0 = xi,  x_(t+1) = T(k+t+1) (A x_t + B u_t) + c(k+t+1).
     * This is what frees Gam from being stored — see bs_workspace_size.  The
     * affine part c is added AFTER the rotation and indexed by the
     * DESTINATION tick, exactly as SchedAff and both closed-loop runners in
     * model/sim_corner_online.py apply it; it is the zero vector on a
     * schedule without offsets, so the record path is unchanged.  Carrying
     * the drift here is also what makes the barrier rows right: they read
     * states[] and so pick up Hr Cbar without a separate constant. */
    {
        bs_real x[NX], ax[NX];
        memcpy(x, xi, sizeof(x));
        for (int t = 0; t < NN; ++t) {
            const bs_real *T = rot_of(schedule, phase + t + 1);
            const bs_real *c = off_of(schedule, phase + t + 1);
            for (int i = 0; i < NX; ++i) {
                bs_real acc = 0.0;
                for (int j = 0; j < NX; ++j)
                    acc += bs_A[(size_t)i * NX + j] * x[j];
                for (int j = 0; j < NU; ++j)
                    acc += bs_B[(size_t)i * NU + j] * U[t * NU + j];
                ax[i] = acc;
            }
            for (int i = 0; i < NX; ++i) {
                bs_real acc = 0.0;
                for (int j = 0; j < NX; ++j)
                    acc += T[(size_t)i * NX + j] * ax[j];
                acc += c[i];
                states[t * NX + i] = acc;
                x[i] = acc;
            }
        }
    }

    bs_real quad = 0.0;
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j)
            quad += xi[i] * bs_Q[(size_t)i * NX + j] * xi[j];
    for (int t = 0; t < NN; ++t)
        for (int i = 0; i < NU; ++i)
            for (int j = 0; j < NU; ++j)
                quad += U[t * NU + i] * bs_R[(size_t)i * NU + j] * U[t * NU + j];
    const bs_real *Pterm = schedule->P_tab
        ? &schedule->P_tab[(size_t)terminal_family_of(schedule, phase + NN)
                           * NX * NX]
#ifndef BS_TABLES_RUNTIME
        : &bs_P[(size_t)terminal_family_of(schedule, phase + NN) * NX * NX];
#else
        : schedule->P_tab;        /* unreachable by construction */
#endif
    for (int t = 0; t < NN; ++t) {
        const bs_real *W = (t == NN - 1) ? Pterm : bs_Q;
        const bs_real *x = &states[t * NX];
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NX; ++j)
                quad += x[i] * W[(size_t)i * NX + j] * x[j];
#ifdef BS_SPEED_TILT_W
        /* SPEED TILT (experimental, compiled out by default): each stage
         * cost gains -w * delta_t, a linear reward on the tangential
         * pace, so the optimum moves from "match the published pace" to
         * "ride the forward face".  Linear in the state: the Hessian --
         * and with it the contraction machinery -- is untouched, exactly
         * like the affine drift.  The matching gradient constant is
         * accumulated into Dquad at build time. */
        quad += -(bs_real)BS_SPEED_TILT_W * x[0];
#endif
    }

    if (grad) {
        /* The quadratic half of the gradient is Hquad U + Cquad xi + Dquad.
         * Dquad is the AFFINE DRIFT term and is easy to miss: the reference
         * writes g = 2(Rbar U + Gam' Wt (Abar xi + Gam U + Cbar)), and the
         * Cbar summand depends on NEITHER xi nor U, so it never appears in a
         * matrix that multiplies one of them.  Dropping it leaves a constant
         * error in every gradient entry -- the objective and Hessian stay
         * right, and the solve quietly converges to the wrong point. */
        for (int i = pin; i < NV; ++i) {
            bs_real acc = 0.0;
            for (int j = 0; j < NV; ++j)
                acc += problem->Hquad[(size_t)i * NV + j] * U[j];
            for (int j = 0; j < NX; ++j)
                acc += problem->Cquad[(size_t)i * NX + j] * xi[j];
            grad[i - pin] = acc + problem->Dquad[i];
        }
    }
    if (hess) {
        for (int a = pin; a < NV; ++a)
            memcpy(&hess[(size_t)(a - pin) * nf],
                   &problem->Hquad[(size_t)a * NV + pin],
                   sizeof(bs_real) * (size_t)nf);
    }

    /* Gam block row, rolled forward alongside the stage loop.  Only the
     * derivative path needs it: a value-only evaluation (every Armijo trial)
     * reads the rows straight off the states and skips Gam and Rmat
     * entirely.  With npin > 0 the roll starts at t = npin from Gam_(npin-1)
     * restricted to columns >= 3*npin, which is IDENTICALLY ZERO, so no seed
     * is needed and no prefix stage is ever touched. */
    bs_real *gam = problem->gam;
    bs_real *gam_next = problem->gam_next;
    if (want_derivatives && t_start == 0) {
        memset(gam, 0, sizeof(bs_real) * (size_t)NX * NV);
        const bs_real *T0 = rot_of(schedule, phase + 1);
        bs_real TB0[NX * NU];
        mat_mul(T0, bs_B, TB0, NX, NX, NU);
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NU; ++j)
                gam[(size_t)i * NV + j] = TB0[(size_t)i * NU + j];
    }

    bs_real barrier = 0.0;
    for (int t = 0; t < NN; ++t) {
        const bs_real *rows = rows_of(schedule, phase + t);
        const bs_real *state = (t == 0) ? xi : &states[(t - 1) * NX];
        const int width = NU * (t + 1);
        const int deriv = want_derivatives && (t >= t_start);

        if (deriv) {
            /* Rmat[t] = Hr Gam_(t-1) with Er in the t-th input block, over
             * the FREE columns only; at t = t_start the Gam term vanishes
             * (t = 0: no Gam; t = npin: Gam_(npin-1) has no free columns). */
            for (int i = 0; i < NR; ++i) {
                bs_real *out = &Rmat[(size_t)i * NV];
                for (int j = pin; j < width; ++j) out[j] = 0.0;
                if (t > t_start) {
                    for (int p = 0; p < NX; ++p) {
                        const bs_real h = rows[(size_t)i * 10 + p];
                        if (h == 0.0) continue;
                        const bs_real *grow = &gam[(size_t)p * NV];
                        for (int j = pin; j < NU * t; ++j)
                            out[j] += h * grow[j];
                    }
                }
                for (int j = 0; j < NU; ++j)
                    out[t * NU + j] += rows[(size_t)i * 10 + 6 + j];
            }
        }

        for (int i = 0; i < NR; ++i) {
            const bs_real *row = &rows[(size_t)i * 10];
            bs_real g_i = 0.0;
            for (int p = 0; p < NX; ++p) g_i += row[p] * state[p];
            for (int p = 0; p < NU; ++p) g_i += row[6 + p] * U[t * NU + p];
            const bs_real m = row[9];
            const bs_real z = m - g_i;

            barrier += BS_EPS * (barrier_value(z) - barrier_value(m)
                                 + barrier_slope(m) * (m - z));
            if (grad && deriv) {
                const bs_real w = BS_EPS * (barrier_slope(m)
                                            - barrier_slope(z));
                const bs_real *rm = &Rmat[(size_t)i * NV];
                for (int j = pin; j < width; ++j) grad[j - pin] += w * rm[j];
            }
            if (hess && deriv) {
                /* Deferred: the weight is stored and the whole stage is
                 * accumulated as one tiled block update below.  See the
                 * comment there for why this is not a micro-optimization. */
                const bs_real clipped = (z >= BS_DR) ? z : BS_DR;
                weight[i] = BS_EPS / (clipped * clipped);
            }
        }

        if (hess && deriv) {
            /* BLOCKED SYMMETRIC RANK-k UPDATE
             *   H[p:w, p:w] += R_t[:, p:w]' diag(weight) R_t[:, p:w]
             * with p = 3*npin and w = NU (t+1).
             *
             * Mathematically identical to the 50 rank-one updates it
             * replaces -- the same sum, reassociated.  Causality and symmetry
             * are still respected: the loops stop at `width` and at the
             * diagonal, and now also start at the commitment boundary. */
            for (int a0 = pin; a0 < width; a0 += BS_HESS_TILE) {
                const int a1 = (a0 + BS_HESS_TILE < width)
                             ? a0 + BS_HESS_TILE : width;
                for (int b0 = pin; b0 <= a0; b0 += BS_HESS_TILE) {
                    const int b1 = (b0 + BS_HESS_TILE < width)
                                 ? b0 + BS_HESS_TILE : width;
                    for (int i = 0; i < NR; ++i) {
                        const bs_real w = weight[i];
                        if (w == 0.0) continue;
                        const bs_real *rm = &Rmat[(size_t)i * NV];
                        for (int a = a0; a < a1; ++a) {
                            const bs_real scaled = w * rm[a];
                            if (scaled == 0.0) continue;
                            bs_real *hrow = &hess[(size_t)(a - pin) * nf];
                            const int bend = (b1 < a + 1) ? b1 : a + 1;
                            for (int b = b0; b < bend; ++b)
                                hrow[b - pin] += scaled * rm[b];
                        }
                    }
                }
            }
        }

        /* Advance the rolling block row.  INDEXING, because it is easy to
         * get wrong and the parity gate is what caught it: at stage t the row
         * map needs Gam_(t-1), so `gam` must hold Gam_(t-1) on entry to stage
         * t.  For npin = 0 it is seeded with Gam_0 before the loop (used
         * first at t = 1); for npin > 0 the advance at t = npin builds
         * Gam_npin from an empty free part.  Either way the advance at the
         * end of stage t uses rotation T(k+t+1), multiplies over the free
         * columns Gam_(t-1) occupies, and places the new B block at column
         * NU*t. */
        if (want_derivatives && t >= t_adv && t + 1 < NN) {
            const bs_real *T = rot_of(schedule, phase + t + 1);
            bs_real TA[NX * NX], TB[NX * NU];
            mat_mul(T, bs_A, TA, NX, NX, NX);
            mat_mul(T, bs_B, TB, NX, NX, NU);
            for (int i = 0; i < NX; ++i) {
                bs_real *out = &gam_next[(size_t)i * NV];
                for (int j = pin; j < NU * t; ++j) {
                    bs_real acc = 0.0;
                    for (int p = 0; p < NX; ++p)
                        acc += TA[(size_t)i * NX + p] * gam[(size_t)p * NV + j];
                    out[j] = acc;
                }
                for (int j = 0; j < NU; ++j)
                    out[NU * t + j] = TB[(size_t)i * NU + j];
            }
            bs_real *swap = gam;
            gam = gam_next;
            gam_next = swap;
        }
    }

    if (hess) {
        /* mirror the lower triangle accumulated above */
        for (int a = 0; a < nf; ++a)
            for (int b = 0; b < a; ++b)
                hess[(size_t)b * nf + a] = hess[(size_t)a * nf + b];
    }
    if (value) *value = quad + barrier;
    return BS_OK;
}

bs_status bs_eval(const bs_problem *problem, const bs_real *U,
                  const bs_real *xi, bs_real *value,
                  bs_real *grad, bs_real *hess)
{
    return eval_impl(problem, U, xi, 0, value, grad, hess);
}

bs_status bs_eval_pinned(const bs_problem *problem, const bs_real *U,
                         const bs_real *xi, int npin, bs_real *value,
                         bs_real *grad, bs_real *hess)
{
    return eval_impl(problem, U, xi, npin, value, grad, hess);
}

void bs_terminal_state(const bs_problem *problem, const bs_real *U,
                       const bs_real *xi, bs_real *xi_terminal)
{
    /* xi_N is simply the last predicted state, so the plant recursion gives
     * it without touching Gam. */
    const bs_schedule *schedule = problem->schedule;
    bs_real x[NX], ax[NX];
    memcpy(x, xi, sizeof(x));
    for (int t = 0; t < NN; ++t) {
        const bs_real *T = rot_of(schedule, problem->phase + t + 1);
        const bs_real *c = off_of(schedule, problem->phase + t + 1);
        for (int i = 0; i < NX; ++i) {
            bs_real acc = 0.0;
            for (int j = 0; j < NX; ++j)
                acc += bs_A[(size_t)i * NX + j] * x[j];
            for (int j = 0; j < NU; ++j)
                acc += bs_B[(size_t)i * NU + j] * U[t * NU + j];
            ax[i] = acc;
        }
        for (int i = 0; i < NX; ++i) {
            bs_real acc = 0.0;
            for (int j = 0; j < NX; ++j)
                acc += T[(size_t)i * NX + j] * ax[j];
            x[i] = acc + c[i];
        }
    }
    memcpy(xi_terminal, x, sizeof(x));
}

void bs_shift_append(const bs_problem *problem, const bs_real *U,
                     const bs_real *xi, bs_real *U_shifted)
{
    bs_real terminal[NX];
    bs_terminal_state(problem, U, xi, terminal);
    for (int i = 0; i < NV - NU; ++i) U_shifted[i] = U[i + NU];
    const bs_real *K = problem->schedule->K_tab
        ? &problem->schedule->K_tab[
              (size_t)terminal_family_of(problem->schedule,
                                         problem->phase + NN) * NU * NX]
#ifndef BS_TABLES_RUNTIME
        : &bs_K[(size_t)terminal_family_of(problem->schedule,
                                           problem->phase + NN) * NU * NX];
#else
        : problem->schedule->K_tab;  /* unreachable by construction */
#endif
    for (int i = 0; i < NU; ++i) {
        bs_real acc = 0.0;
        for (int j = 0; j < NX; ++j) acc += K[(size_t)i * NX + j] * terminal[j];
        U_shifted[NV - NU + i] = acc;
    }
}

bs_status bs_chol_factor(bs_real *A, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            bs_real acc = A[(size_t)i * n + j];
            for (int p = 0; p < j; ++p)
                acc -= A[(size_t)i * n + p] * A[(size_t)j * n + p];
            if (i == j) {
                if (!(acc > 0.0)) return BS_ERR_FACTOR;
                A[(size_t)i * n + j] = sqrt(acc);
            } else {
                A[(size_t)i * n + j] = acc / A[(size_t)j * n + j];
            }
        }
        for (int j = i + 1; j < n; ++j) A[(size_t)i * n + j] = 0.0;
    }
    return BS_OK;
}

void bs_chol_solve(const bs_real *L, int n, bs_real *b)
{
    for (int i = 0; i < n; ++i) {
        bs_real acc = b[i];
        for (int p = 0; p < i; ++p) acc -= L[(size_t)i * n + p] * b[p];
        b[i] = acc / L[(size_t)i * n + i];
    }
    for (int i = n - 1; i >= 0; --i) {
        bs_real acc = b[i];
        for (int p = i + 1; p < n; ++p) acc -= L[(size_t)p * n + i] * b[p];
        b[i] = acc / L[(size_t)i * n + i];
    }
}

static bs_status newton_impl(const bs_problem *problem, bs_real *U,
                             const bs_real *xi, int q, int npin,
                             bs_real *lambda0, bs_newton_stats *stats)
{
    /* Mirrors MixedHorizon.newton_q() in the certificates, GUARDS INCLUDED.
     * The guards are not defensive garnish — they are what the reference
     * does, and two of them decide behaviour over most of a maintained
     * mission:
     *   (a) no descent direction        -> take NO step;
     *   (b) predicted decrease below    -> take NO step.  Without this a
     *       numerical noise                converged iterate fails Armijo 40
     *                                      times (the objective change is
     *                                      pure rounding) and then moves by a
     *                                      meaningless sliver;
     *   (c) Armijo exhausted            -> take NO step, and count it.
     * The B2 mission gate is what exposed their absence; B1's random states
     * never reach the converged regime where they bind. */
    bs_real *const hess = problem->nt_hess;
    bs_real *const grad = problem->nt_grad;
    bs_real *const direction = problem->nt_dir;
    bs_real *const trial = problem->nt_trial;
    const bs_real c1 = 0.1, rho_bt = 0.5;
    const int pin = NU * npin;      /* first FREE entry of the plan */
    const int nf = NV - pin;        /* free-tail dimension          */

    if (stats) { stats->backtracks = 0; stats->armijo_fail = 0;
                 stats->factor_fail = 0; stats->no_step = 0; }
    if (npin < 0 || npin > NN) return BS_ERR_PHASE;
    if (lambda0) *lambda0 = 0.0;
    if (nf == 0) return BS_OK;      /* nothing free: the plan is all data */

    for (int iter = 0; iter < q; ++iter) {
        bs_real value = 0.0;
        bs_status st = eval_impl(problem, U, xi, npin, &value, grad, hess);
        if (st != BS_OK) return st;

        for (int i = 0; i < nf; ++i) direction[i] = grad[i];
        if (bs_chol_factor(hess, nf) != BS_OK) {
            if (stats) ++stats->factor_fail;
            break;
        }
        bs_chol_solve(hess, nf, direction);      /* direction = H_ff^-1 g_f */

        /* gd = g . d with d supported on the free tail and d_f = -H_ff^-1 g_f,
         * i.e. gd = -(g_f' H_ff^-1 g_f) <= 0.  The pinned entries of d are
         * zero, so the full inner product reduces to the free one exactly —
         * which is what newton_q_pinned() computes as gr @ d. */
        bs_real quadratic = 0.0;
        for (int i = 0; i < nf; ++i) quadratic += grad[i] * direction[i];
        if (iter == 0 && lambda0)
            *lambda0 = sqrt(quadratic > 0.0 ? quadratic : 0.0) / sqrt(BS_EPS);

        const bs_real gd = -quadratic;
        if (!(gd < -1e-300)) {                              /* (a) */
            if (stats) ++stats->no_step;
            break;
        }
        const bs_real magnitude = fabs(value) > 1.0 ? fabs(value) : 1.0;
        if (-gd / 2.0 <= 1e-13 * magnitude) {               /* (b) */
            if (stats) ++stats->no_step;
            break;
        }

        bs_real alpha = 1.0;
        int accepted = 0;
        for (int attempt = 0; attempt < 60; ++attempt) {
            for (int i = 0; i < pin; ++i) trial[i] = U[i];
            for (int i = 0; i < nf; ++i)
                trial[pin + i] = U[pin + i] - alpha * direction[i];
            bs_real candidate = 0.0;
            if (eval_impl(problem, trial, xi, npin, &candidate, NULL, NULL)
                    == BS_OK
                && candidate <= value + c1 * alpha * gd) {
                accepted = 1;
                break;
            }
            alpha *= rho_bt;
            if (stats) ++stats->backtracks;
        }
        if (!accepted) {                                    /* (c) */
            if (stats) ++stats->armijo_fail;
            break;
        }
        for (int i = 0; i < nf; ++i) U[pin + i] -= alpha * direction[i];
    }
    return BS_OK;
}

bs_status bs_newton(const bs_problem *problem, bs_real *U, const bs_real *xi,
                    int q, bs_real *lambda0, bs_newton_stats *stats)
{
    return newton_impl(problem, U, xi, q, 0, lambda0, stats);
}

/* The INTERLEAVED-cadence operator: solve every tick, but the first `npin`
 * input stages are committed DATA — replayed, not revised — and only the free
 * tail moves.  Mirrors MixedHorizon.newton_q_pinned() in wp5_anytime_sim.py
 * line for line, guards included; it is the object the interleaved
 * certificate wp5_h1_interleaved_q.md is stated over (uniform q >= 1,
 * varrho = 4.220230e-2). */
bs_status bs_newton_pinned(const bs_problem *problem, bs_real *U,
                           const bs_real *xi, int q, int npin,
                           bs_real *lambda0, bs_newton_stats *stats)
{
    return newton_impl(problem, U, xi, q, npin, lambda0, stats);
}
