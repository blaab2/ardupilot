#include "cs_solver.h"
#include "cs_arena.h"
#include "cs_condense.h"
#include "cs_linalg.h"
#include "cs_model.h"

#include "../model/generated/turn_r5_probdata.h"

/* vendored DAQP kernel subset (qp/daqp, see PROVENANCE) */
#include "../qp/daqp/daqp.h"
#include "../qp/daqp/utils.h"

/* the QP kernel and the solver must agree on the scalar type */
typedef char cs_real_matches_c_float
    [(sizeof(cs_real) == sizeof(c_float)) ? 1 : -1];

static const cs_real k_lbu_s[CS_PD_NU] = CS_PD_LBU;
static const cs_real k_ubu_s[CS_PD_NU] = CS_PD_UBU;

/* DAQP workspace arrays, held BY VALUE in the solver struct so the sizing
 * pass (fake arena base) never dereferences an arena pointer — the real init
 * wires them into the arena-carved DAQPWorkspace afterwards. */
typedef struct {
    c_float *Rinv, *M, *scaling, *v, *dupper, *dlower;
    c_float *lam, *lam_star, *D, *xldl, *zldl, *L, *u, *xold;
    int *sense, *WS;
} cs_daqp_arrays;

struct cs_solver {
    int N, nx, nu;
    cs_qp qp;

    /* iterate */
    cs_real *X;                /* (N+1)*nx, stage-major */
    cs_real *U;                /* N*nu                  */
    cs_real T;
    cs_real *x0, *xN;
    int have_pins, have_seed;

    /* DAQP */
    DAQPWorkspace *work;
    DAQPProblem *dqp;
    DAQPSettings *settings;
    cs_daqp_arrays da;

    /* warm-start memory (active set of the previous QP) */
    int *ws_ids;
    unsigned char *ws_lower;
    int ws_n;

    /* last-iteration info */
    cs_real step_norm, defect, slack_max, slack_soft, kappa;
    int qp_status, qp_iters;
    size_t used_bytes;
};

/* All allocations in ONE place, shared by the sizing pass (fake base, no
 * stores into arena memory — only pointer bookkeeping in *s, which lives on
 * the caller's stack during sizing) and the real init. */
static int cs_carve(cs_solver *s, cs_arena *arena, int N,
                    void **model_mem, size_t *model_bytes)
{
    int rc, nz, m, ldl;
    *model_bytes = cs_model_min_arena();
    *model_mem = cs_arena_alloc_default(arena, *model_bytes);
    if (!*model_mem)
        return CS_ERR_ARENA;
    rc = cs_condense_init(&s->qp, arena, N);
    if (rc != CS_OK)
        return rc;
    /* DAQP workspace sized for the WORST-CASE pattern (cs_condense_init's
     * nz_max/m_max) so cs_solver_set_pins can re-dim without reallocation */
    nz = s->qp.nz_max;
    m = s->qp.m_max;
    ldl = nz + 1;
    s->N = N;
    s->nx = s->qp.nx;
    s->nu = s->qp.nu;
    s->X = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)(N + 1) * s->qp.nx * sizeof(cs_real));
    s->U = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)N * s->qp.nu * sizeof(cs_real));
    s->x0 = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)s->qp.nx * sizeof(cs_real));
    s->xN = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)s->qp.nx * sizeof(cs_real));
    s->work = (DAQPWorkspace *)cs_arena_alloc_default(
        arena, sizeof(DAQPWorkspace));
    s->dqp = (DAQPProblem *)cs_arena_alloc_default(arena, sizeof(DAQPProblem));
    s->settings = (DAQPSettings *)cs_arena_alloc_default(
        arena, sizeof(DAQPSettings));
    s->ws_ids = (int *)cs_arena_alloc_default(arena, (size_t)ldl * sizeof(int));
    s->ws_lower = (unsigned char *)cs_arena_alloc_default(arena, (size_t)ldl);
    if (!s->X || !s->U || !s->x0 || !s->xN || !s->work || !s->dqp ||
        !s->settings || !s->ws_ids || !s->ws_lower)
        return CS_ERR_ARENA;
    /* DAQP workspace arrays (allocate_daqp_workspace + setup_daqp_ldp sizes,
     * ns = 0: elastic columns are explicit QP columns, not DAQP softs) */
    {
        cs_daqp_arrays *d = &s->da;
        d->Rinv = (c_float *)cs_arena_alloc_default(
            arena, (size_t)nz * (nz + 1) / 2 * sizeof(c_float));
        d->M = (c_float *)cs_arena_alloc_default(
            arena, (size_t)nz * s->qp.mA_max * sizeof(c_float));
        d->scaling = (c_float *)cs_arena_alloc_default(
            arena, (size_t)m * sizeof(c_float));
        d->v = (c_float *)cs_arena_alloc_default(
            arena, (size_t)nz * sizeof(c_float));
        d->dupper = (c_float *)cs_arena_alloc_default(
            arena, (size_t)m * sizeof(c_float));
        d->dlower = (c_float *)cs_arena_alloc_default(
            arena, (size_t)m * sizeof(c_float));
        d->sense = (int *)cs_arena_alloc_default(arena, (size_t)m * sizeof(int));
        d->lam = (c_float *)cs_arena_alloc_default(
            arena, (size_t)ldl * sizeof(c_float));
        d->lam_star = (c_float *)cs_arena_alloc_default(
            arena, (size_t)ldl * sizeof(c_float));
        d->WS = (int *)cs_arena_alloc_default(arena, (size_t)ldl * sizeof(int));
        d->D = (c_float *)cs_arena_alloc_default(
            arena, (size_t)ldl * sizeof(c_float));
        d->xldl = (c_float *)cs_arena_alloc_default(
            arena, (size_t)ldl * sizeof(c_float));
        d->zldl = (c_float *)cs_arena_alloc_default(
            arena, (size_t)ldl * sizeof(c_float));
        d->L = (c_float *)cs_arena_alloc_default(
            arena, (size_t)ldl * (ldl + 1) / 2 * sizeof(c_float));
        d->u = (c_float *)cs_arena_alloc_default(
            arena, (size_t)nz * sizeof(c_float));
        d->xold = (c_float *)cs_arena_alloc_default(
            arena, (size_t)nz * sizeof(c_float));
        if (!d->Rinv || !d->M || !d->scaling || !d->v || !d->dupper ||
            !d->dlower || !d->sense || !d->lam || !d->lam_star || !d->WS ||
            !d->D || !d->xldl || !d->zldl || !d->L || !d->u || !d->xold)
            return CS_ERR_ARENA;
    }
    return CS_OK;
}

size_t cs_solver_min_arena(int N)
{
    /* sizing pass: fake aligned base, solver struct on OUR stack — the arena
     * pointers are computed but never dereferenced. +64 alignment slack for
     * the struct carve of the real init. */
    cs_solver dummy;
    cs_arena a;
    void *mm;
    size_t mb;
    int ok = 0;
    {   /* is N a compiled variant? */
        int i;
        for (i = 0; i < cs_model_num_variants(); ++i)
            if (cs_model_variant_N(i) == N)
                ok = 1;
    }
    if (!ok)
        return 0;
    cs_arena_init(&a, (void *)(uintptr_t)64, (size_t)1 << 30);
    cs_arena_alloc_default(&a, sizeof(cs_solver));
    if (cs_carve(&dummy, &a, N, &mm, &mb) != CS_OK)
        return 0;
    return cs_arena_high_water(&a) + 64;
}

cs_solver *cs_solver_init(void *mem, size_t size, int N, int *rc)
{
    cs_arena arena;
    cs_solver *s;
    void *model_mem;
    size_t model_bytes;
    int st, i;
    if (rc)
        *rc = CS_ERR_ARG;
    if (!mem || cs_solver_min_arena(N) == 0)
        return (cs_solver *)0;
    cs_arena_init(&arena, mem, size);
    s = (cs_solver *)cs_arena_alloc_default(&arena, sizeof(cs_solver));
    if (!s) {
        if (rc)
            *rc = CS_ERR_ARENA;
        return (cs_solver *)0;
    }
    st = cs_carve(s, &arena, N, &model_mem, &model_bytes);
    if (st == CS_OK)
        st = cs_model_init(model_mem, model_bytes);
    if (st != CS_OK) {
        if (rc)
            *rc = st;
        return (cs_solver *)0;
    }
    /* values */
    s->T = (cs_real)CS_PD_T_REF;
    s->have_pins = s->have_seed = 0;
    s->ws_n = 0;
    s->step_norm = (cs_real)-1;
    s->defect = s->slack_max = s->slack_soft = s->kappa = (cs_real)0;
    s->qp_status = 0;
    s->qp_iters = 0;
    cs_zero(s->X, (size_t)(N + 1) * s->nx);
    cs_zero(s->U, (size_t)N * s->nu);
    /* DAQP wiring (the api.c roles, arena edition) */
    s->dqp->n = s->qp.nz;
    s->dqp->m = s->qp.m;
    s->dqp->ms = s->qp.nz;          /* simple bounds on ALL variables */
    s->dqp->H = s->qp.H;
    s->dqp->f = s->qp.f;
    s->dqp->A = s->qp.A;
    s->dqp->bupper = s->qp.bu;
    s->dqp->blower = s->qp.bl;
    s->dqp->sense = s->qp.sense;
    s->dqp->bin_ids = (int *)0;
    s->dqp->nb = 0;
    /* daqp_default_settings (api.c, replicated — api.c is not vendored) */
    s->settings->primal_tol = (c_float)DEFAULT_PRIM_TOL;
    s->settings->dual_tol = (c_float)DEFAULT_DUAL_TOL;
    s->settings->zero_tol = (c_float)DEFAULT_ZERO_TOL;
    s->settings->pivot_tol = (c_float)DEFAULT_PIVOT_TOL;
    s->settings->progress_tol = (c_float)DEFAULT_PROG_TOL;
#ifdef CS_REAL_SINGLE
    /* float-honest QP tolerances (plan 4.4: the ~1e-5/1e-6 design scale).
     * DAQP's double-scale defaults (dual 1e-12, zero 1e-11, progress 1e-14)
     * sit BELOW float epsilon — the singularity/feasibility guards would
     * compare against pure rounding noise and can never trigger. Measured
     * A/B on repro_warmstart N=20/30 (see csolver/README.md M2): identical
     * QP path either way on this problem — these values change nothing
     * measured and make the guards meaningful in float. */
    s->settings->primal_tol = (c_float)1e-5;
    s->settings->dual_tol = (c_float)1e-7;
    s->settings->zero_tol = (c_float)1e-6;
    s->settings->pivot_tol = (c_float)1e-4;
    s->settings->progress_tol = (c_float)1e-7;
#endif
    s->settings->cycle_tol = DEFAULT_CYCLE_TOL;
    s->settings->iter_limit = DEFAULT_ITER_LIMIT;
    s->settings->fval_bound = DAQP_INF;
    s->settings->eps_prox = 0;
    s->settings->eta_prox = (c_float)DEFAULT_ETA;
    s->settings->rho_soft = (c_float)DEFAULT_RHO_SOFT;
    s->settings->rel_subopt = 0;
    s->settings->abs_subopt = 0;
    s->work->qp = s->dqp;
    s->work->n = s->qp.nz;
    s->work->m = s->qp.m;
    s->work->ms = s->qp.nz;
    s->work->Rinv = s->da.Rinv;
    s->work->M = s->da.M;
    s->work->scaling = s->da.scaling;
    s->work->v = s->da.v;
    s->work->dupper = s->da.dupper;
    s->work->dlower = s->da.dlower;
    s->work->sense = s->da.sense;
    s->work->lam = s->da.lam;
    s->work->lam_star = s->da.lam_star;
    s->work->WS = s->da.WS;
    s->work->D = s->da.D;
    s->work->xldl = s->da.xldl;
    s->work->zldl = s->da.zldl;
    s->work->L = s->da.L;
    s->work->u = s->da.u;
    s->work->x = s->da.u;           /* DAQP aliases x onto u */
    s->work->xold = s->da.xold;
    s->work->settings = s->settings;
    s->work->bnb = (DAQPBnB *)0;
    s->work->soft_slack = 0;
    for (i = 0; i < s->qp.m; ++i)
        s->work->sense[i] = 0;
    reset_daqp_workspace(s->work);
    s->used_bytes = cs_arena_high_water(&arena);
    if (rc)
        *rc = CS_OK;
    return s;
}

int cs_solver_set_pins(cs_solver *s, const cs_real *x0, const cs_real *xN)
{
    int free0[6], pinN[6], f0cur[6] = {0, 0, 0, 0, 0, 0};
    int i, changed = 0, rc;
    if (!s || !x0 || !xN)
        return CS_ERR_ARG;
    /* NaN = free (M2 free-boundary support): a NaN entry component becomes
     * a decision column, a NaN terminal component drops its pin (gaining
     * the gated boundary h/box rows instead). All-finite pins = the M1
     * fully-pinned problem, layout-identical. */
    for (i = 0; i < s->nx; ++i) {
        free0[i] = (x0[i] != x0[i]) ? 1 : 0;
        pinN[i] = (xN[i] != xN[i]) ? 0 : 1;
    }
    for (i = 0; i < s->qp.nf0; ++i)
        f0cur[s->qp.free0[i]] = 1;
    for (i = 0; i < s->nx; ++i)
        if (f0cur[i] != free0[i] || s->qp.pinN[i] != pinN[i])
            changed = 1;
    if (changed) {
        rc = cs_condense_set_pattern(&s->qp, free0, pinN);
        if (rc != CS_OK)
            return rc;
        /* re-dim the (max-sized) DAQP workspace onto the new layout */
        s->dqp->n = s->qp.nz;
        s->dqp->m = s->qp.m;
        s->dqp->ms = s->qp.nz;
        s->work->n = s->qp.nz;
        s->work->m = s->qp.m;
        s->work->ms = s->qp.nz;
        for (i = 0; i < s->qp.m; ++i)
            s->work->sense[i] = 0;
        reset_daqp_workspace(s->work);
        s->ws_n = 0;
    }
    cs_copy(x0, s->x0, s->nx);
    cs_copy(xN, s->xN, s->nx);
    s->have_pins = 1;
    return CS_OK;
}

int cs_solver_set_xn_box(cs_solver *s, int idx, cs_real lo, cs_real hi)
{
    if (!s || idx < 0 || idx >= s->nx || !(lo <= hi))
        return CS_ERR_ARG;
    s->qp.xn_lo[idx] = lo;
    s->qp.xn_hi[idx] = hi;
    return CS_OK;
}

int cs_solver_set_seed(cs_solver *s, const cs_real *X, const cs_real *U,
                       const cs_real *T, int clip_rollout)
{
    int k, j, rc;
    if (!s || !X || !U || !T)
        return CS_ERR_ARG;
    cs_copy(X, s->X, (s->N + 1) * s->nx);
    cs_copy(U, s->U, s->N * s->nu);
    s->T = *T;
    if (clip_rollout) {
        if (s->T < (cs_real)CS_PD_T_MIN)
            s->T = (cs_real)CS_PD_T_MIN;
        if (s->T > (cs_real)CS_PD_T_MAX)
            s->T = (cs_real)CS_PD_T_MAX;
        for (k = 0; k < s->N; ++k)
            for (j = 0; j < s->nu; ++j) {
                cs_real *u = s->U + (size_t)k * s->nu + j;
                if (*u < k_lbu_s[j])
                    *u = k_lbu_s[j];
                if (*u > k_ubu_s[j])
                    *u = k_ubu_s[j];
            }
        for (k = 0; k < s->N; ++k) {
            rc = cs_step(s->N, s->X + (size_t)k * s->nx,
                         s->U + (size_t)k * s->nu, &s->T, s->qp.wind,
                         s->X + (size_t)(k + 1) * s->nx);
            if (rc != CS_OK)
                return rc;
        }
    }
    /* fresh problem => fresh QP warm start */
    s->ws_n = 0;
    s->have_seed = 1;
    return CS_OK;
}

/* diag-ratio kappa proxy of the scaled dense H (plan-4.4 instrumentation) */
static cs_real kappa_proxy(const cs_qp *qp)
{
    int i;
    cs_real lo = (cs_real)0, hi = (cs_real)0;
    for (i = 0; i < qp->nz; ++i) {
        cs_real d = qp->H[(size_t)i * qp->nz + i];
        if (i == 0 || d < lo)
            lo = d;
        if (i == 0 || d > hi)
            hi = d;
    }
    return (lo > (cs_real)0) ? hi / lo : (cs_real)-1;
}

int cs_solver_iterate(cs_solver *s)
{
    DAQPWorkspace *work;
    int rc, i, flag;
    if (!s)
        return CS_ERR_ARG;
    if (!s->have_pins || !s->have_seed)
        return CS_ERR_INIT;
    work = s->work;

    rc = cs_condense_build(&s->qp, s->X, s->U, s->T, s->x0, s->xN);
    if (rc != CS_OK)
        return rc;
    s->defect = s->qp.defect_max;
    s->kappa = kappa_proxy(&s->qp);

    /* sense bits: terminal pins are IMMUTABLE equalities, always in the
     * working set (the dense_qp_daqp convention for equality rows — with
     * the +-1 slack columns they can never be linearly dependent); the
     * previous active set is the RTI warm start. */
    for (i = 0; i < s->qp.npin; ++i)
        s->qp.sense[s->qp.m - s->qp.npin + i] = ACTIVE + LOWER + IMMUTABLE;
    for (i = 0; i < s->ws_n; ++i)
        s->qp.sense[s->ws_ids[i]] |= ACTIVE + (s->ws_lower[i] ? LOWER : 0);

    /* two-step update: sense FIRST so normalize_M's zero-row IMMUTABLE
     * marking (vendored utils.c) is not clobbered by the sense copy */
    update_ldp(UPDATE_sense, work);
    rc = update_ldp(UPDATE_Rinv + UPDATE_M + UPDATE_v + UPDATE_d, work);
    if (rc < 0) {
        s->qp_status = rc;
        s->qp_iters = 0;
        s->step_norm = (cs_real)-1;
        s->ws_n = 0;
        return CS_OK;                 /* iterate frozen, status reported */
    }
    rc = activate_constraints(work);
    if (rc < 0) {
        /* overdetermined warm start: strip the mutable actives, keep pins */
        reset_daqp_workspace(work);
        for (i = 0; i < s->qp.m; ++i)
            if (!(work->sense[i] & IMMUTABLE))
                work->sense[i] &= ~ACTIVE;
        rc = activate_constraints(work);
        if (rc < 0) {
            s->qp_status = EXIT_OVERDETERMINED_INITIAL;
            s->qp_iters = 0;
            s->step_norm = (cs_real)-1;
            s->ws_n = 0;
            return CS_OK;
        }
    }

    flag = daqp_ldp(work);
    s->qp_status = flag;
    s->qp_iters = work->iterations;
    if (flag == EXIT_OPTIMAL || flag == EXIT_SOFT_OPTIMAL) {
        ldp2qp_solution(work);
        s->step_norm = cs_condense_expand(&s->qp, work->x, s->X, s->U, &s->T);
        /* slack_max = terminal elastic pins only (M1/M2 semantics: 0 at the
         * r5 fixed point); slack_soft = the M1.2 corridor-block + box-upper
         * slacks (the shrinking-horizon residual violation reporter). */
        s->slack_max = cs_norm_inf(work->x + s->qp.nz0, 2 * s->qp.npin);
        s->slack_soft = (s->qp.nsc > 0)
            ? cs_norm_inf(work->x + s->qp.nz0 + 2 * s->qp.npin, s->qp.nsc)
            : (cs_real)0;
        /* save the active set for the next QP (skip the immutable pins —
         * they are re-activated structurally every iteration) */
        s->ws_n = 0;
        for (i = 0; i < work->n_active; ++i) {
            int id = work->WS[i];
            if (work->sense[id] & IMMUTABLE)
                continue;
            s->ws_ids[s->ws_n] = id;
            s->ws_lower[s->ws_n] = (work->sense[id] & LOWER) ? 1 : 0;
            ++s->ws_n;
        }
    } else {
        s->step_norm = (cs_real)-1;   /* QP failure froze the iterate */
        s->slack_max = s->slack_soft = (cs_real)0;
        s->ws_n = 0;
    }
    return CS_OK;
}

int cs_solver_get_iterate(const cs_solver *s, cs_real *X, cs_real *U,
                          cs_real *T)
{
    if (!s || !X || !U || !T)
        return CS_ERR_ARG;
    cs_copy(s->X, X, (s->N + 1) * s->nx);
    cs_copy(s->U, U, s->N * s->nu);
    *T = s->T;
    return CS_OK;
}

int cs_solver_get_info(const cs_solver *s, cs_real *step_norm,
                       cs_real *defect, cs_real *slack_max, cs_real *kappa,
                       int *qp_status, int *qp_iters)
{
    if (!s)
        return CS_ERR_ARG;
    if (step_norm) *step_norm = s->step_norm;
    if (defect)    *defect = s->defect;
    if (slack_max) *slack_max = s->slack_max;
    if (kappa)     *kappa = s->kappa;
    if (qp_status) *qp_status = s->qp_status;
    if (qp_iters)  *qp_iters = s->qp_iters;
    return CS_OK;
}

int cs_solver_get_dims(const cs_solver *s, int *N, int *nx, int *nu,
                       int *nz, int *m)
{
    if (!s)
        return CS_ERR_ARG;
    if (N)  *N = s->N;
    if (nx) *nx = s->nx;
    if (nu) *nu = s->nu;
    if (nz) *nz = s->qp.nz;
    if (m)  *m = s->qp.m;
    return CS_OK;
}

int cs_solver_set_opts(cs_solver *s, cs_real lm, cs_real rho,
                       cs_real slack_z, cs_real slack_Z, int jerk_end)
{
    if (!s || lm < (cs_real)0 || rho < (cs_real)0 ||
        slack_z < (cs_real)0 || slack_Z <= (cs_real)0)
        return CS_ERR_ARG;
    s->qp.lm = lm;
    s->qp.rho = rho;
    s->qp.slack_z = slack_z;
    s->qp.slack_Z = slack_Z;
    s->qp.jerk_end = jerk_end ? 1 : 0;
    return CS_OK;
}

int cs_solver_set_soft_weights(cs_solver *s, cs_real zc, cs_real Zc)
{
    if (!s || zc < (cs_real)0 || Zc <= (cs_real)0)
        return CS_ERR_ARG;
    s->qp.slack_zc = zc;
    s->qp.slack_Zc = Zc;
    return CS_OK;
}

int cs_solver_set_ref(cs_solver *s, const cs_real *x_ref, cs_real w, int mode)
{
    if (!s || w < (cs_real)0 || mode < 0 || mode > 3)
        return CS_ERR_ARG;
    if (x_ref)
        cs_copy(x_ref, s->qp.x_ref, (size_t)(s->N + 1) * s->nx);
    s->qp.w_ref = w;
    s->qp.ref_mode = mode;
    return CS_OK;
}

int cs_solver_set_crop_bound(cs_solver *s, int on, cs_real margin)
{
    if (!s)
        return CS_ERR_ARG;
    return cs_condense_set_crop(&s->qp, on, margin);
}

int cs_solver_set_gate(cs_solver *s, cs_real tau0, cs_real rlx_out,
                       cs_real rlx_ret, cs_real rlx_xi)
{
    if (!s)
        return CS_ERR_ARG;
    if (tau0 >= (cs_real)0 &&
        (tau0 > (cs_real)1 || rlx_out < (cs_real)0 ||
         rlx_ret < (cs_real)0 || rlx_xi < (cs_real)0))
        return CS_ERR_ARG;
    s->qp.gate_tau0 = tau0;
    s->qp.gate_rlx_out = rlx_out;
    s->qp.gate_rlx_ret = rlx_ret;
    s->qp.gate_rlx_xi = rlx_xi;
    return CS_OK;
}

int cs_solver_set_row_offset(cs_solver *s, cs_real d_minus_14p1)
{
    if (!s || !(d_minus_14p1 == d_minus_14p1))      /* NaN guard */
        return CS_ERR_ARG;
    s->qp.row_off = d_minus_14p1;
    return CS_OK;
}

int cs_solver_set_wind(cs_solver *s, cs_real w_xi, cs_real w_eta)
{
    if (!s || !(w_xi == w_xi) || !(w_eta == w_eta)) /* NaN guard */
        return CS_ERR_ARG;
    s->qp.wind[0] = w_xi;
    s->qp.wind[1] = w_eta;
    return CS_OK;
}

int cs_solver_set_soft_corr(cs_solver *s, int on)
{
    int i, rc;
    if (!s)
        return CS_ERR_ARG;
    rc = cs_condense_set_soft(&s->qp, on);
    if (rc != CS_OK)
        return rc;
    /* re-dim the (max-sized) DAQP workspace onto the new column count */
    s->dqp->n = s->qp.nz;
    s->dqp->m = s->qp.m;
    s->dqp->ms = s->qp.nz;
    s->work->n = s->qp.nz;
    s->work->m = s->qp.m;
    s->work->ms = s->qp.nz;
    for (i = 0; i < s->qp.m; ++i)
        s->work->sense[i] = 0;
    reset_daqp_workspace(s->work);
    s->ws_n = 0;
    return CS_OK;
}

int cs_solver_get_soft_slack(const cs_solver *s, cs_real *soft)
{
    if (!s || !soft)
        return CS_ERR_ARG;
    *soft = s->slack_soft;
    return CS_OK;
}

int cs_solver_get_qp(const cs_solver *s, cs_real *A, cs_real *bl, cs_real *bu,
                     int *sense)
{
    int i;
    if (!s)
        return CS_ERR_ARG;
    if (A)
        cs_copy(s->qp.A, A, s->qp.mA * s->qp.nz);
    if (bl)
        cs_copy(s->qp.bl, bl, s->qp.m);
    if (bu)
        cs_copy(s->qp.bu, bu, s->qp.m);
    if (sense)
        for (i = 0; i < s->qp.m; ++i)
            sense[i] = s->qp.sense[i];
    return CS_OK;
}

size_t cs_solver_high_water(const cs_solver *s)
{
    return s ? s->used_bytes : 0;
}
