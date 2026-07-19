#include "cs_condense.h"
#include "cs_linalg.h"
#include "cs_model.h"

#include "../model/generated/turn_r5_probdata.h"

#define CS_BIG ((cs_real)1e20)   /* "no bound" for simple bounds / one-sided rows */

static const int k_bx_idx[CS_PD_NBX] = CS_PD_BX_IDX;
static const cs_real k_bx_lb[CS_PD_NBX] = CS_PD_BX_LB;
static const cs_real k_bx_ub[CS_PD_NBX] = CS_PD_BX_UB;
static const cs_real k_lbu[CS_PD_NU] = CS_PD_LBU;
static const cs_real k_ubu[CS_PD_NU] = CS_PD_UBU;
/* full per-state boxes (one-sided kept; +-1e20 = unbounded) */
static const cs_real k_x_lb[6] = CS_PD_X_LB;
static const cs_real k_x_ub[6] = CS_PD_X_UB;
/* structural state dependencies of the h rows [jerk, chi, ap, re, xi]
 * (row 0 also depends on u) — which boundary h rows are non-constant
 * under a partial pin pattern */
static const int k_hdep[5] = { 1 << 4, (1 << 1) | (1 << 5),
                               (1 << 2) | (1 << 3), (1 << 2) | (1 << 3),
                               1 << 2 };

/* rows a boundary pattern adds at stage 0 (h rows made non-constant by the
 * free entry states) and stage N (h rows + box rows for free terminal
 * states; box rows only where some bound is finite) */
static void pattern_rows(const int *free0_flags, const int *pinN_flags,
                         int *h0_extra, int *hN, int *boxN)
{
    int i, fm0 = 0, fmN = 0;
    for (i = 0; i < 6; ++i) {
        if (free0_flags[i])
            fm0 |= 1 << i;
        if (!pinN_flags[i])
            fmN |= 1 << i;
    }
    *h0_extra = 0;
    *hN = 0;
    *boxN = 0;
    for (i = 1; i < 5; ++i) {          /* x-only h rows (jerk row is u-row) */
        if (k_hdep[i] & fm0)
            ++*h0_extra;
        if (k_hdep[i] & fmN)
            ++*hN;
    }
    for (i = 0; i < 6; ++i)
        if (!pinN_flags[i] &&
            (k_x_lb[i] > -CS_BIG || k_x_ub[i] < CS_BIG))
            ++*boxN;
}

int cs_condense_set_pattern(cs_qp *qp, const int *free0_flags,
                            const int *pinN_flags)
{
    int i, h0_extra, hN, boxN;
    qp->nf0 = 0;
    qp->npin = 0;
    for (i = 0; i < qp->nx; ++i) {
        if (free0_flags[i])
            qp->free0[qp->nf0++] = i;
        qp->pinN[i] = pinN_flags[i] ? 1 : 0;
        if (qp->pinN[i])
            ++qp->npin;
        qp->xn_lo[i] = k_x_lb[i];
        qp->xn_hi[i] = k_x_ub[i];
    }
    qp->iT = qp->N * qp->nu;
    qp->nz0 = qp->iT + 1 + qp->nf0;
    qp->nsc = qp->soft_corr ? (CS_NSC + CS_PD_NBX) : 0;
    qp->ns = 2 * qp->npin + qp->nsc;
    qp->nz = qp->nz0 + qp->ns;
    pattern_rows(free0_flags, pinN_flags, &h0_extra, &hN, &boxN);
    qp->mA = (2 + h0_extra)
        + (qp->N - 1) * (qp->nh + 1 + CS_PD_NBX)
        + hN + boxN + qp->npin;
    qp->m = qp->nz + qp->mA;
    if (qp->nz > qp->nz_max || qp->mA > qp->mA_max || qp->m > qp->m_max)
        return CS_ERR_ARG;
    return CS_OK;
}

int cs_condense_set_soft(cs_qp *qp, int on)
{
    qp->soft_corr = on ? 1 : 0;
    qp->nsc = qp->soft_corr ? (CS_NSC + CS_PD_NBX) : 0;
    qp->ns = 2 * qp->npin + qp->nsc;
    qp->nz = qp->nz0 + qp->ns;
    qp->m = qp->nz + qp->mA;
    if (qp->nz > qp->nz_max || qp->m > qp->m_max)
        return CS_ERR_ARG;
    return CS_OK;
}

int cs_condense_init(cs_qp *qp, cs_arena *arena, int N)
{
    int nx = CS_PD_NX, nu = CS_PD_NU, nh = CS_PD_NH;
    int all_pinned[6] = {1, 1, 1, 1, 1, 1};
    int none_free[6] = {0, 0, 0, 0, 0, 0};
    qp->N = N; qp->nx = nx; qp->nu = nu; qp->nh = nh;
    /* worst-case layout for the allocation: all entry states free + all
     * terminal states pinned (columns), all boundary h/box rows (rows), PLUS
     * the M1.2 elastic soft-state reserve (CS_NSC corridor-block slacks +
     * CS_PD_NBX box-upper slacks) — the flight rh path turns it on, so
     * CS_PINNED_ONLY reserves it too. CS_PINNED_ONLY: no free-boundary. */
#ifdef CS_PINNED_ONLY
    qp->nz_max = N * nu + 1 + 2 * nx + CS_NSC + CS_PD_NBX;
    qp->mA_max = 2 + (N - 1) * (nh + 1 + CS_PD_NBX) + nx;
#else
    {
        int all_free[6] = {1, 1, 1, 1, 1, 1};
        int h0_extra, hN, boxN;
        qp->nz_max = N * nu + 1 + nx + 2 * nx + CS_NSC + CS_PD_NBX;
        pattern_rows(all_free, none_free, &h0_extra, &hN, &boxN);
        qp->mA_max = (2 + h0_extra)
            + (N - 1) * (nh + 1 + CS_PD_NBX) + hN + boxN + nx;
    }
#endif
    qp->m_max = qp->nz_max + qp->mA_max;
    qp->n_out = (int)(CS_PD_FRAC_OUT * N);
    qp->n_ret = (int)(CS_PD_FRAC_RET * N);

    qp->lm = (cs_real)1e-1;       /* measured fixed point (acados LM sweep) */
    qp->rho = (cs_real)1e-3;
    qp->slack_z = (cs_real)1e3;   /* exact-penalty L1, physical units */
    qp->slack_Z = (cs_real)1e4;
    qp->slack_zc = (cs_real)1e2;  /* corridor slack (acados rate1: zl 1e2) */
    qp->slack_Zc = (cs_real)1e3;  /* corridor slack L2 (acados: Zl 1e3)   */
    qp->soft_corr = 0;            /* hard corridor by default (M1 parity) */
    qp->relax = (cs_real)1e9;     /* the acados _RELAX */
    qp->jerk_end = 1;
    qp->gate_tau0 = (cs_real)-1;  /* full-horizon index gating (batch) */
    qp->gate_rlx_out = qp->gate_rlx_ret = qp->gate_rlx_xi = (cs_real)0;
    qp->row_off = (cs_real)0;     /* canonical d = 14.1 (G0.2) */
    cs_scaling_default(&qp->sc);
    {   /* LM state metric: physical (lm_w = sx^2) — the M1/acados-parity
         * damping; E/e themselves are stored scaled (M2). */
        int i;
        for (i = 0; i < nx; ++i)
            qp->lm_w[i] = qp->sc.x[i] * qp->sc.x[i];
    }

    qp->H = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)qp->nz_max * qp->nz_max * sizeof(cs_real));
    qp->f = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)qp->nz_max * sizeof(cs_real));
    qp->A = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)qp->mA_max * qp->nz_max * sizeof(cs_real));
    qp->bl = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)qp->m_max * sizeof(cs_real));
    qp->bu = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)qp->m_max * sizeof(cs_real));
    qp->sense = (int *)cs_arena_alloc_default(
        arena, (size_t)qp->m_max * sizeof(int));
    qp->E = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)(N + 1) * nx * (qp->nz_max - 2 * nx) * sizeof(cs_real));
    qp->e = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)(N + 1) * nx * sizeof(cs_real));
    qp->Ak = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)nx * nx * sizeof(cs_real));
    qp->Bk = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)nx * nu * sizeof(cs_real));
    qp->bTk = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)nx * sizeof(cs_real));
    qp->xnext = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)nx * sizeof(cs_real));
    qp->hk = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)nh * sizeof(cs_real));
    qp->Jxk = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)nh * nx * sizeof(cs_real));
    qp->Juk = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)nh * nu * sizeof(cs_real));
    qp->jrow = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)nx * sizeof(cs_real));
    if (!qp->H || !qp->f || !qp->A || !qp->bl || !qp->bu || !qp->sense ||
        !qp->E || !qp->e || !qp->Ak || !qp->Bk || !qp->bTk || !qp->xnext ||
        !qp->hk || !qp->Jxk || !qp->Juk || !qp->jrow)
        return CS_ERR_ARENA;
    /* default pattern = the M1 layout: entry fully pinned, nx terminal pins */
    return cs_condense_set_pattern(qp, none_free, all_pinned);
}

/* one h-type constraint row: A[r,:] = Jxrow' * E_k (nz0 cols) + Ju entries at
 * the u_k columns; bounds shifted by h + Jxrow'e_k.
 * E/e are stored in SCALED state units (M2 staged state scaling): the
 * physical Jx row is folded with the state scales, Jx~ = Jx * Sx, so
 * Jx~' * E~ and Jx~' * e~ reproduce the exact physical row/shift values —
 * the constraint set is IDENTICAL to the unscaled M1 build. */
static void fill_h_row(cs_qp *qp, int r, int k,
                       const cs_real *Jxrow, const cs_real *Jurow,
                       cs_real hval, cs_real lo, cs_real hi)
{
    const cs_real *Ek = qp->E + (size_t)k * qp->nx * qp->nz0;
    const cs_real *ek = qp->e + (size_t)k * qp->nx;
    cs_real *row = qp->A + (size_t)r * qp->nz;
    cs_real jxs[6];
    cs_acc shift = (cs_acc)hval;
    int j;
    for (j = 0; j < qp->nx; ++j)
        jxs[j] = Jxrow[j] * qp->sc.x[j];
    cs_row_gemv(qp->nx, qp->nz0, jxs, Ek, row);     /* row[0..nz0) */
    if (k < qp->N && Jurow)
        for (j = 0; j < qp->nu; ++j)
            row[k * qp->nu + j] += Jurow[j] * qp->sc.u[j];
    for (j = 0; j < qp->nx; ++j)
        shift += (cs_acc)jxs[j] * (cs_acc)ek[j];
    qp->bl[qp->nz + r] = (cs_real)(lo - shift);
    qp->bu[qp->nz + r] = (cs_real)(hi - shift);
}

/* conservative end-of-interval jerk row at stage k (analytic linearization) */
static void fill_jerk_end_row(cs_qp *qp, int r, int k,
                              const cs_real *X, const cs_real *U, cs_real T)
{
    const cs_real a = X[(size_t)k * qp->nx + 4];
    const cs_real ad = U[(size_t)k * qp->nu + 0];
    const cs_real td = U[(size_t)k * qp->nu + 1];
    const cs_real hN = T / (cs_real)qp->N;
    const cs_real ae = a + ad * hN;
    const cs_real g = ad * ad + ae * ae * td * td;
    const cs_real dg_da = (cs_real)2 * ae * td * td;
    cs_real Jxrow[6] = {0, 0, 0, 0, 0, 0};
    cs_real Jurow[2];
    cs_real *row;
    Jxrow[4] = dg_da;
    Jurow[0] = (cs_real)2 * ad + dg_da * hN;
    Jurow[1] = (cs_real)2 * ae * ae * td;
    fill_h_row(qp, r, k, Jxrow, Jurow, g,
               qp->jerk_end ? -qp->relax : -CS_BIG,
               qp->jerk_end ? (cs_real)(CS_PD_J_MAX * CS_PD_J_MAX) : CS_BIG);
    /* direct T dependence (beyond the E-propagated a(T) part) */
    row = qp->A + (size_t)r * qp->nz;
    row[qp->iT] += dg_da * (ad / (cs_real)qp->N) * qp->sc.T_ref;
}

/* gated bounds of the nh h rows at stage k (batch index gating or the
 * M1.1 rh progress gating) — shared by the interior and boundary stages.
 * G0.2: the RETURN row (h[3] = eta - reI(xi)) is row-relative to row 2 —
 * its gate bound shifts rigidly by row_off = d - 14.1 (the shifted problem
 * is eta >= reI(xi) + row_off, i.e. h[3] >= row_off). Approach/crossing
 * are row-1-relative and d-independent. row_off = 0 is bit-unchanged. */
static void gate_bounds(const cs_qp *qp, int k, cs_real *lo, cs_real *hi)
{
    lo[0] = -qp->relax;            hi[0] = (cs_real)(CS_PD_J_MAX * CS_PD_J_MAX);
    lo[1] = (cs_real)CS_PD_CHI_LO; hi[1] = (cs_real)CS_PD_CHI_HI;
    lo[2] = -qp->relax;            hi[2] = qp->relax;
    lo[3] = -qp->relax;            hi[3] = qp->relax;
    lo[4] = -qp->relax;            hi[4] = qp->relax;
    if (qp->gate_tau0 < (cs_real)0) {
        /* batch: the original full-horizon index gating */
        if (k <= qp->n_out)
            hi[2] = (cs_real)0;
        if (k >= qp->n_ret)
            lo[3] = qp->row_off;
        if (k > qp->n_out && k < qp->n_ret)
            lo[4] = (cs_real)CS_PD_XI_CROSS_MIN;
    } else {
        /* rh: gate by the stage's ORIGINAL progress (M1.1) */
        const cs_real tk = qp->gate_tau0
            + ((cs_real)1 - qp->gate_tau0)
              * (cs_real)k / (cs_real)qp->N;
        if (tk <= (cs_real)CS_PD_FRAC_OUT)
            hi[2] = qp->gate_rlx_out;
        if (tk >= (cs_real)CS_PD_FRAC_RET)
            lo[3] = qp->row_off - qp->gate_rlx_ret;
        if (tk > (cs_real)CS_PD_FRAC_OUT &&
            tk < (cs_real)CS_PD_FRAC_RET)
            lo[4] = (cs_real)CS_PD_XI_CROSS_MIN - qp->gate_rlx_xi;
    }
}

int cs_condense_build(cs_qp *qp, const cs_real *X, const cs_real *U,
                      cs_real T, const cs_real *x0, const cs_real *xN)
{
    const int N = qp->N, nx = qp->nx, nu = qp->nu, nh = qp->nh;
    const int nz0 = qp->nz0, nz = qp->nz, npin = qp->npin;
    const int iT = qp->iT, nf0 = qp->nf0;
    const cs_real Tref = qp->sc.T_ref;
    int free0_flag[6] = {0, 0, 0, 0, 0, 0};
    int k, i, j, r, rc;
    for (i = 0; i < nf0; ++i)
        free0_flag[qp->free0[i]] = 1;

    /* ---- 1. sweep: sensitivities, defects, E/e propagation ----
     * M2 staged state scaling (plan 4.4): E/e carry the state increments in
     * SCALED units, dx~ = Sx^-1 dx, so the whole condensing recursion runs
     * on O(1) quantities:
     *   A~_k = Sx^-1 A_k Sx,  B~_k = Sx^-1 B_k Su,  bT~_k = Sx^-1 bT_k Tref,
     *   e~_{k+1} = A~_k e~_k + Sx^-1 c_k.
     * defect_max stays PHYSICAL (measured before scaling). */
    cs_zero(qp->E, (size_t)(N + 1) * nx * nz0);
    qp->defect_max = (cs_real)0;
    /* dx_0: pinned comps enter e_0; FREE comps are decision columns
     * (unit entries in E_0 at cols iT+1..iT+nf0, scaled by sx) */
    for (i = 0; i < nx; ++i)
        qp->e[i] = free0_flag[i] ? (cs_real)0
                                 : (x0[i] - X[i]) / qp->sc.x[i];
    for (i = 0; i < nf0; ++i)
        qp->E[(size_t)(iT + 1 + i) * nx + qp->free0[i]] = (cs_real)1;
    for (k = 0; k < N; ++k) {
        const cs_real *Xk = X + (size_t)k * nx;
        const cs_real *Uk = U + (size_t)k * nu;
        const cs_real *Ek = qp->E + (size_t)k * nx * nz0;
        const cs_real *ek = qp->e + (size_t)k * nx;
        cs_real *En = qp->E + (size_t)(k + 1) * nx * nz0;
        cs_real *en = qp->e + (size_t)(k + 1) * nx;
        rc = cs_step_jac(N, Xk, Uk, &T, qp->xnext, qp->Ak, qp->Bk, qp->bTk);
        if (rc != CS_OK)
            return rc;
        /* e~_{k+1} = A~_k e~_k + Sx^-1 c_k,  c_k = F(X_k,U_k,T) - X_{k+1} */
        for (i = 0; i < nx; ++i) {
            cs_real c = qp->xnext[i] - X[(size_t)(k + 1) * nx + i];
            cs_real d = c < (cs_real)0 ? -c : c;
            if (d > qp->defect_max)
                qp->defect_max = d;
            en[i] = c / qp->sc.x[i];
        }
        /* scale the stage sensitivities in place (col-major nx x nx / nu) */
        for (j = 0; j < nx; ++j)
            for (i = 0; i < nx; ++i)
                qp->Ak[(size_t)j * nx + i] *= qp->sc.x[j] / qp->sc.x[i];
        for (j = 0; j < nu; ++j)
            for (i = 0; i < nx; ++i)
                qp->Bk[(size_t)j * nx + i] *= qp->sc.u[j] / qp->sc.x[i];
        for (i = 0; i < nx; ++i)
            qp->bTk[i] *= Tref / qp->sc.x[i];
        cs_gemm_acc(nx, 1, nx, qp->Ak, ek, en);
        /* E~_{k+1}: u-prefix cols [0, k*nu), the s_T col + the nf0 free-dx0
         * cols (contiguous block [iT, nz0)), then the new u_k cols = B~_k.
         * Cols in between stay zero (E zeroed above). */
        if (k > 0)
            cs_gemm(nx, k * nu, nx, qp->Ak, Ek, En);
        cs_gemm(nx, 1 + nf0, nx, qp->Ak, Ek + (size_t)iT * nx,
                En + (size_t)iT * nx);
        for (i = 0; i < nx; ++i)
            En[(size_t)iT * nx + i] += qp->bTk[i];
        for (j = 0; j < nu; ++j)
            for (i = 0; i < nx; ++i)
                En[(size_t)(k * nu + j) * nx + i] = qp->Bk[(size_t)j * nx + i];
    }

    /* ---- 2. objective: H / f ----
     * The LM state term is lm * E~' diag(lm_w) E~ over the SCALED E: with
     * lm_w = sx^2 this is the EXACT physical-metric acados damping (M1
     * parity); with lm_w = 1 it is the scaled-metric variant (see the
     * cs_condense_init comment for the measured choice). */
    cs_zero(qp->H, (size_t)nz * nz);
    cs_zero(qp->f, nz);
    for (k = 0; k < N; ++k)
        for (j = 0; j < nu; ++j) {
            const int c = k * nu + j;
            const cs_real su = qp->sc.u[j];
            qp->H[(size_t)c * nz + c] +=
                (qp->rho / (cs_real)N + qp->lm) * su * su;
            qp->f[c] += (qp->rho / (cs_real)N) * U[(size_t)k * nu + j] * su;
        }
    for (k = 0; k <= N; ++k) {
        const cs_real *Ek = qp->E + (size_t)k * nx * nz0;
        const cs_real *ek = qp->e + (size_t)k * nx;
        cs_syrk_w_acc(nz0, nx, qp->lm, qp->lm_w, Ek, qp->H, nz);
        cs_gemv_tw_acc(nx, nz0, qp->lm, qp->lm_w, Ek, ek, qp->f);
        qp->H[(size_t)iT * nz + iT] += qp->lm * Tref * Tref;
    }
    qp->H[(size_t)iT * nz + iT] += Tref * Tref;                 /* 0.5 T^2 */
    qp->f[iT] += T * Tref;
    for (i = 0; i < qp->ns; ++i) {
        const int c = nz0 + i;
        const int cor = (i >= 2 * npin);   /* corridor slack vs terminal pin */
        qp->H[(size_t)c * nz + c] += cor ? qp->slack_Zc : qp->slack_Z;
        qp->f[c] += cor ? qp->slack_zc : qp->slack_z;
    }

    /* ---- 3. simple bounds on all nz variables ---- */
    for (k = 0; k < N; ++k)
        for (j = 0; j < nu; ++j) {
            const int c = k * nu + j;
            const cs_real u = U[(size_t)k * nu + j], su = qp->sc.u[j];
            qp->bl[c] = (k_lbu[j] - u) / su;
            qp->bu[c] = (k_ubu[j] - u) / su;
        }
    qp->bl[iT] = ((cs_real)CS_PD_T_MIN - T) / Tref;
    qp->bu[iT] = ((cs_real)CS_PD_T_MAX - T) / Tref;
    for (i = 0; i < nf0; ++i) {
        /* free-dx0 columns: state box relative to the current iterate,
         * in scaled units (one-sided boxes keep their infinite side) */
        const int f = qp->free0[i];
        const cs_real sx = qp->sc.x[f];
        qp->bl[iT + 1 + i] = (k_x_lb[f] <= -CS_BIG)
            ? -CS_BIG : (k_x_lb[f] - X[f]) / sx;
        qp->bu[iT + 1 + i] = (k_x_ub[f] >= CS_BIG)
            ? CS_BIG : (k_x_ub[f] + (f == 3 ? qp->row_off : (cs_real)0)
                        - X[f]) / sx;                     /* G0.2 eta shift */
    }
    for (i = 0; i < qp->ns; ++i) {
        qp->bl[nz0 + i] = (cs_real)0;
        qp->bu[nz0 + i] = CS_BIG;
    }

    /* ---- 4. general rows ---- */
    cs_zero(qp->A, (size_t)qp->mA * nz);
    r = 0;
    /* stage 0: jerk node row (u_0) + jerk end row; with FREE entry states,
     * additionally the x-only h rows they make non-constant (all-pinned
     * pattern: constants — omitted to avoid zero rows, the M1 layout) */
    rc = cs_h_jac(X, U, qp->hk, qp->Jxk, qp->Juk);
    if (rc != CS_OK)
        return rc;
    {
        cs_real lo[5], hi[5];
        int fm0 = 0;
        for (j = 0; j < nx; ++j)
            if (free0_flag[j])
                fm0 |= 1 << j;
        gate_bounds(qp, 0, lo, hi);
        for (i = 0; i < nh; ++i) {
            cs_real Jxrow[6], Jurow[2];
            if (i > 0 && !(k_hdep[i] & fm0))
                continue;
            for (j = 0; j < nx; ++j)
                Jxrow[j] = qp->Jxk[i + (size_t)nh * j];
            for (j = 0; j < nu; ++j)
                Jurow[j] = qp->Juk[i + (size_t)nh * j];
            fill_h_row(qp, r++, 0, Jxrow, Jurow, qp->hk[i], lo[i], hi[i]);
        }
    }
    fill_jerk_end_row(qp, r++, 0, X, U, T);
    /* stages 1..N-1: gated h rows + jerk end + box rows */
    for (k = 1; k < N; ++k) {
        const cs_real *Xk = X + (size_t)k * nx;
        cs_real lo[5], hi[5];
        const cs_real gth = qp->relax * (cs_real)0.5;  /* gated < gth < relax */
        /* shared per-block corridor slack columns: [out, cross, ret] */
        const int c_out = qp->soft_corr ? (nz0 + 2 * npin + 0) : -1;
        const int c_cross = qp->soft_corr ? (nz0 + 2 * npin + 1) : -1;
        const int c_ret = qp->soft_corr ? (nz0 + 2 * npin + 2) : -1;
        rc = cs_h_jac(Xk, U + (size_t)k * nu, qp->hk, qp->Jxk, qp->Juk);
        if (rc != CS_OK)
            return rc;
        gate_bounds(qp, k, lo, hi);
        for (i = 0; i < nh; ++i) {
            cs_real Jxrow[6], Jurow[2];
            for (j = 0; j < nx; ++j)
                Jxrow[j] = qp->Jxk[i + (size_t)nh * j];
            for (j = 0; j < nu; ++j)
                Jurow[j] = qp->Juk[i + (size_t)nh * j];
            fill_h_row(qp, r, k, Jxrow, Jurow, qp->hk[i], lo[i], hi[i]);
            /* M1.2: soften the ONE gated corridor row of this stage with the
             * shared nonnegative slack of its block. Outbound is an UPPER
             * bound (coeff -1: A z - s <= hi[2]); return/crossing are LOWER
             * bounds (coeff +1: A z + s >= lo[3|4]). Exactly one fires per
             * stage (the tk partition); a single unbounded-above slack per
             * block absorbs that block's worst residual, so every corridor
             * row is feasible => the condensed QP is always primal-feasible
             * on the corridor and DAQP never returns -1 mid-turn. */
            if (qp->soft_corr) {
                cs_real *row = qp->A + (size_t)r * nz;
                if (i == 2 && hi[2] < gth)
                    row[c_out] = (cs_real)-1;
                else if (i == 3 && lo[3] > -gth)
                    row[c_ret] = (cs_real)1;
                else if (i == 4 && lo[4] > -gth)
                    row[c_cross] = (cs_real)1;
            }
            ++r;
        }
        fill_jerk_end_row(qp, r++, k, X, U, T);
        for (i = 0; i < CS_PD_NBX; ++i) {
            const int s = k_bx_idx[i];
            /* G0.2: the eta box upper (the AP lateral peak, measured at
             * d = 14.1) is row-2-relative — shift by row_off = d - 14.1 */
            const cs_real ub = k_bx_ub[i]
                + (s == 3 ? qp->row_off : (cs_real)0);
            cs_real Jxrow[6] = {0, 0, 0, 0, 0, 0};
            Jxrow[s] = (cs_real)1;
            fill_h_row(qp, r, k, Jxrow, (const cs_real *)0,
                       Xk[s], k_bx_lb[i], ub);
            /* M1.2: soften the box UPPER bound with a shared nonneg slack
             * (coeff -1: x_s - s <= ub). v<=V_INF and a<=A_CAP ride their
             * ceilings near the exit/apex, and the shrinking-horizon
             * linearization overshoots them by cm — the LOWER bounds
             * (v>=V_MIN huge margin, a>=A_TRIM the load-bearing floor) stay
             * HARD. Column: after the CS_NSC corridor-block slacks. */
            if (qp->soft_corr)
                (qp->A + (size_t)r * nz)[nz0 + 2 * npin + CS_NSC + i] =
                    (cs_real)-1;
            ++r;
        }
    }
    /* stage N (free-boundary patterns only): the x-only h rows made
     * non-constant by the free terminal states (e.g. the chi row the R5
     * exit rides at -0.2), then box rows for the free terminal states
     * (xn_lo/xn_hi — the full state box unless overridden, e.g. the P1
     * exit-speed inequality v_N >= 0.999 V_INF) */
    {
        int fmN = 0;
        for (j = 0; j < nx; ++j)
            if (!qp->pinN[j])
                fmN |= 1 << j;
        if (fmN) {
            cs_real lo[5], hi[5];
            rc = cs_h_jac(X + (size_t)N * nx, U + (size_t)(N - 1) * nu,
                          qp->hk, qp->Jxk, qp->Juk);
            if (rc != CS_OK)
                return rc;
            gate_bounds(qp, N, lo, hi);
            for (i = 1; i < nh; ++i) {          /* x-only rows */
                cs_real Jxrow[6];
                if (!(k_hdep[i] & fmN))
                    continue;
                for (j = 0; j < nx; ++j)
                    Jxrow[j] = qp->Jxk[i + (size_t)nh * j];
                fill_h_row(qp, r++, N, Jxrow, (const cs_real *)0,
                           qp->hk[i], lo[i], hi[i]);
            }
            for (i = 0; i < nx; ++i) {
                cs_real Jxrow[6] = {0, 0, 0, 0, 0, 0};
                if (qp->pinN[i] ||
                    (k_x_lb[i] <= -CS_BIG && k_x_ub[i] >= CS_BIG))
                    continue;
                Jxrow[i] = (cs_real)1;
                fill_h_row(qp, r++, N, Jxrow, (const cs_real *)0,
                           X[(size_t)N * nx + i],
                           qp->xn_lo[i] <= -CS_BIG ? -CS_BIG : qp->xn_lo[i],
                           qp->xn_hi[i] >= CS_BIG ? CS_BIG : qp->xn_hi[i]);
            }
        }
    }
    /* terminal: npin elastic equality pins, slack columns +-1 */
    {
        int p = 0;
        for (i = 0; i < nx; ++i) {
            cs_real Jxrow[6] = {0, 0, 0, 0, 0, 0};
            cs_real *row;
            cs_real target;
            if (!qp->pinN[i])
                continue;
            target = xN[i] - X[(size_t)N * nx + i];
            Jxrow[i] = (cs_real)1;
            fill_h_row(qp, r, N, Jxrow, (const cs_real *)0, (cs_real)0,
                       target, target);
            row = qp->A + (size_t)r * nz;
            row[nz0 + p] = (cs_real)1;             /* s_lo */
            row[nz0 + npin + p] = (cs_real)-1;     /* s_up */
            ++p;
            ++r;
        }
    }
    if (r != qp->mA)
        return CS_ERR_INIT;

    /* sense: fresh each build; the solver ORs its warm-start bits in */
    for (i = 0; i < qp->m; ++i)
        qp->sense[i] = 0;
    return CS_OK;
}

cs_real cs_condense_expand(const cs_qp *qp, const cs_real *zeta,
                           cs_real *X, cs_real *U, cs_real *T)
{
    const int N = qp->N, nx = qp->nx, nu = qp->nu, nz0 = qp->nz0;
    cs_real step = (cs_real)0, dT;
    int k, i, j;
    for (k = 0; k <= N; ++k) {
        const cs_real *Ek = qp->E + (size_t)k * nx * nz0;
        const cs_real *ek = qp->e + (size_t)k * nx;
        cs_acc dx[6];                       /* double accumulation (plan 4.4) */
        for (i = 0; i < nx; ++i)
            dx[i] = (cs_acc)ek[i];
        for (j = 0; j < nz0; ++j) {
            const cs_acc z = (cs_acc)zeta[j];
            if (z == (cs_acc)0)
                continue;
            for (i = 0; i < nx; ++i)
                dx[i] += (cs_acc)Ek[(size_t)j * nx + i] * z;
        }
        for (i = 0; i < nx; ++i) {
            cs_real dxi = (cs_real)(dx[i] * (cs_acc)qp->sc.x[i]);  /* unscale */
            cs_real d = dxi < (cs_real)0 ? -dxi : dxi;
            if (d > step)
                step = d;
            X[(size_t)k * nx + i] += dxi;
        }
    }
    for (k = 0; k < N; ++k)
        for (j = 0; j < nu; ++j) {
            cs_real du = zeta[k * nu + j] * qp->sc.u[j];
            cs_real d = du < (cs_real)0 ? -du : du;
            if (d > step)
                step = d;
            U[(size_t)k * nu + j] += du;
        }
    dT = zeta[qp->iT] * qp->sc.T_ref;
    *T += dT;
    if (dT < (cs_real)0)
        dT = -dT;
    if (dT > step)
        step = dT;
    return step;
}
