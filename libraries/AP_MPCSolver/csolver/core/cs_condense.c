#include "cs_condense.h"
#include "cs_linalg.h"
#include "cs_model.h"

#include <math.h>

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

int cs_condense_set_crop(cs_qp *qp, int on, cs_real margin)
{
    if (!qp || margin < (cs_real)0)
        return CS_ERR_ARG;
    qp->crop_bound = on ? 1 : 0;
    qp->crop_margin = margin;
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
    qp->x_ref = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)(N + 1) * nx * sizeof(cs_real));
    qp->rref = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)nx * sizeof(cs_real));
    qp->vref = (cs_real *)cs_arena_alloc_default(
        arena, (size_t)qp->nz_max * sizeof(cs_real));
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
        !qp->E || !qp->e || !qp->x_ref || !qp->rref || !qp->vref ||
        !qp->Ak || !qp->Bk || !qp->bTk || !qp->xnext || !qp->hk ||
        !qp->Jxk || !qp->Juk || !qp->jrow)
        return CS_ERR_ARENA;
    qp->w_ref = (cs_real)0;      /* reference tracking off (canonical); */
    qp->ref_mode = 0;            /* x_ref filled by cs_solver_set_ref    */
    qp->crop_bound = 0;          /* R5-anchored crop floor off (canonical) */
    qp->crop_margin = (cs_real)0;
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
/* R5 crop envelope: the deepest (minimum) xi the reference path x_ref
 * reaches at row-transfer position eta — the crop boundary "as R5 cuts it".
 * Aligned by ETA, not node index: a plan that cuts the corner is AHEAD of
 * the reference in eta, so at equal node index the reference is still on
 * the approach (xi << 0) and a per-node bound never binds. The (xi, eta)
 * path is restretch-invariant, so the envelope is stable as tau0 advances.
 * Multivalued regions (eta flat at 0 / d on the rows, both flares) take the
 * minimum xi over all crossings = the most permissive, matching "no deeper
 * than R5 anywhere at this eta". Outside the reference's eta range, falls
 * back to the eta-nearest node's xi. O(N) scan, N=30. */
#define CS_GATE_FAR_MARGIN ((cs_real)0.05) /* = CS_RH_MARGIN: bare softening
                                            * for stages beyond the vehicle */
#define CS_GATE_VLOC ((cs_real)0.06)       /* ~2 nodes of tk: the pin's
                                            * neighborhood for measured-
                                            * violation relaxations */
#define CS_CROP_ETA_DEAD ((cs_real)0.5) /* floor off within this of a row  */
#define CS_CROP_VEH_WIN ((cs_real)0.3)  /* vehicle pin-accommodation window */
#define CS_CROP_UP_MARGIN ((cs_real)1.0) /* soft bulge cap: ref bulge + this */

/* mode-3 normal tube dead radius: cm-noise, the 0.2-0.5 m mid-turn
 * tracking offsets, and the apex projection flicker stay (mostly) out of
 * the objective; only migration beyond it is priced, via the hinge
 * 0.5*w*(|d_n|-r0)^2. Flight-tuned across three sentinel flights: full
 * residual at w=30 storms the slack validator from ordinary tracking
 * offsets (tube3t); capping the tube to the outbound stages instead frees
 * the apex into the over-bulge branch of the flat bulge<->time manifold
 * (tube3u, bulge 3.8 vs seed 1.98); splicing the index-matched pull onto
 * the tail re-anchors TIMING by node index against the tube's tangential
 * freedom and back-propagates as braking + a re-seeded field-side shift
 * (tube3v). One geometric hinge over the whole horizon is the composition
 * that works: flat-manifold migration (metres) is decided by the hinge,
 * tracking noise is not. */
#define CS_REF_TUBE_R0 ((cs_real)0.25)

/* EXACT segment interpolation, no eta window: the window (a 0.5 m min-
 * neighborhood, added for the flat row segments) softened the steep entry
 * flare by ~0.9 m — R5's xi(eta) rises 2.56 m over 1.6 m of eta there, and
 * the flown flare rode the softened bound ~1 m deep / ~2 m eta-late
 * (measured, PSCN target == actual). The flats' problem is instead handled
 * by the CS_CROP_ETA_DEAD dead-zone in gate_bounds: within 0.5 m of either
 * row, "off the row" is measurement noise, not crop, and the floor is off
 * — which also covers the multivalued terminal-settle band. */
static cs_real crop_env_xi(const cs_qp *qp, cs_real eta)
{
    const cs_real *xr = qp->x_ref;
    const int nx = qp->nx;
    cs_real fl = (cs_real)0, best = CS_BIG, xin = (cs_real)0;
    int j, hit = 0;
    for (j = 0; j <= qp->N; ++j) {
        const cs_real ej = xr[(size_t)j * nx + 3];
        const cs_real xj = xr[(size_t)j * nx + 2];
        cs_real de = ej - eta;
        if (de < (cs_real)0)
            de = -de;
        if (de < best) {
            best = de;
            xin = xj;
        }
        if (j < qp->N) {
            const cs_real e1 = xr[(size_t)(j + 1) * nx + 3];
            const cs_real elo = ej < e1 ? ej : e1;
            const cs_real ehi = ej < e1 ? e1 : ej;
            if (eta >= elo && eta <= ehi &&
                ehi - elo > (cs_real)1e-9) {
                const cs_real x1 = xr[(size_t)(j + 1) * nx + 2];
                const cs_real xi =
                    xj + (eta - ej) / (e1 - ej) * (x1 - xj);
                if (!hit || xi < fl) {
                    fl = xi;
                    hit = 1;
                }
            }
        }
    }
    return hit ? fl : xin;
}

static void gate_bounds(const cs_qp *qp, int k, cs_real eta_k,
                        cs_real *lo, cs_real *hi)
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
        /* LOCAL measured-violation relaxations (the corridor-ratchet fix):
         * the v_out/v_ret relaxations exist ONLY to keep the hard x0 pin
         * feasible when the VEHICLE violates its leg's envelope — applied
         * to the whole leg they are a self-licensing loop: a +2 cm entry
         * cross-track (every post-handback approach settles slightly
         * row-high; a first turn's pristine leg rides row-low, v_out = 0)
         * relaxes the outbound gate everywhere, the replan flares early,
         * the replay flies it, the measured violation grows, and the loop
         * ratchets ~7 cm/cycle into a rigid one-node (-1.5 m) field-side
         * shift on EVERY non-first turn (3-turn fork + MPCC telemetry;
         * seed-anchored tau0 alone did not break it). Stages beyond the
         * vehicle's progress neighborhood keep the bare softening margin. */
        const cs_real far = (cs_real)CS_GATE_FAR_MARGIN;
        const cs_real near_v = tk - qp->gate_tau0 <= (cs_real)CS_GATE_VLOC
            ? (cs_real)1 : (cs_real)0;
        if (tk <= (cs_real)CS_PD_FRAC_OUT) {
            hi[2] = near_v > (cs_real)0 ? qp->gate_rlx_out
                : (qp->gate_rlx_out < far ? qp->gate_rlx_out : far);
        }
        if (tk >= (cs_real)CS_PD_FRAC_RET) {
            const cs_real r = near_v > (cs_real)0 ? qp->gate_rlx_ret
                : (qp->gate_rlx_ret < far ? qp->gate_rlx_ret : far);
            lo[3] = qp->row_off - r;
        }
        if (qp->crop_bound) {
          /* soft bulge cap (whole horizon, eta-independent): see the
           * crop_xi_up build comment — damps the flat-manifold wander */
          if (qp->crop_xi_up < hi[4])
              hi[4] = qp->crop_xi_up;
          if (eta_k > CS_CROP_ETA_DEAD &&
              eta_k < qp->crop_eta_max - CS_CROP_ETA_DEAD) {
            /* R5-anchored crop floor over the transfer (replaces the
             * relaxable vx crossing gate): xi_k >= min(0, env(eta_k)) -
             * crop_margin, env = the reference path's xi at the stage's
             * CURRENT row-transfer position eta_k (crop_env_xi, exact
             * interpolation). The row dead-zone (0.5 m of either row)
             * excludes the flat row segments and the terminal settle,
             * where eta is noise, not crop. Capped at 0 so it never pushes
             * a bulge where R5 sits in the headland; follows R5's own
             * flare depth near the rows; tight through the mid-transfer
             * (floor ~ -crop_margin). The i==4 softening (shared c_cross
             * slack) fires wherever this raises lo[4] above -gth, so the
             * crop bound is elastic. */
            cs_real fl = crop_env_xi(qp, eta_k);
            cs_real dev = eta_k - qp->crop_v_eta;
            if (fl > (cs_real)0)
                fl = (cs_real)0;
            fl -= qp->crop_margin;
            /* local pin accommodation: stages near the VEHICLE's eta admit
             * the measured xi (the hard x0 anchor must stay feasible when
             * the vehicle is already in the crop). LOCAL by design — the
             * old vx gate relaxed the whole crossing by the measured
             * violation and let plans legally sit in the crop; here only
             * the vehicle's eta neighborhood softens, the rest of the
             * transfer holds the R5 line. NaN vehicle (free entry) fails
             * the compare and disables it. */
            if (dev < (cs_real)0)
                dev = -dev;
            if (dev <= CS_CROP_VEH_WIN &&
                qp->crop_v_xi - qp->crop_margin < fl)
                fl = qp->crop_v_xi - qp->crop_margin;
            if (fl > lo[4])
                lo[4] = fl;
          }
          /* dead-zone stages: no bound (crop mode fully replaces the
           * legacy vx crossing gate — no fall-through) */
        } else if (tk > (cs_real)CS_PD_FRAC_OUT &&
                   tk < (cs_real)CS_PD_FRAC_RET) {
            lo[4] = (cs_real)CS_PD_XI_CROSS_MIN - qp->gate_rlx_xi;
        }
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
    /* vehicle (xi, eta) for the crop floor's local pin accommodation
     * (gate_bounds); NaN (free entry) disables it via failed compares */
    qp->crop_v_xi = x0[2];
    qp->crop_v_eta = x0[3];
    /* reference eta max (crop dead-zone row test) and reference bulge
     * (the soft xi UPPER cap): without the two-sided proximity pull the
     * bulge<->time manifold is flat/undamped, and the exact crop floor
     * tips early replans into 7-8 m wide arcs that the flight validator
     * then rejects (REJ_XI bursts -> stale fallbacks). The cap bounds the
     * wander INSIDE the QP — a cap, not a pull: it pushes nothing, plans
     * at/below the reference bulge are untouched. */
    if (qp->crop_bound) {
        cs_real em = qp->x_ref[3], xm = qp->x_ref[2];
        for (k = 1; k <= N; ++k) {
            const cs_real e = qp->x_ref[(size_t)k * nx + 3];
            const cs_real x = qp->x_ref[(size_t)k * nx + 2];
            if (e > em)
                em = e;
            if (x > xm)
                xm = x;
        }
        qp->crop_eta_max = em;
        qp->crop_xi_up = (xm > (cs_real)0 ? xm : (cs_real)0)
            + CS_CROP_UP_MARGIN;
    }

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
    /* R5-proximity: w_ref * sum_k ||x_k - x_ref_k||^2 (Gauss-Newton state
     * cost). Same syrk/gemv machinery as the LM term, but the residual is the
     * iterate-vs-reference deviation r~_k = (X_k - x_ref_k)/sx (scaled), so H
     * gains curvature in the flat bulge<->time mode (mode 1) and f gains the
     * pull toward the R5 reference (mode 1) or its sign (mode 2, the
     * curvature-null L1 surrogate). w_ref = 0 / mode 0 -> QP bit-unchanged. */
    if (qp->w_ref > (cs_real)0 &&
        (qp->ref_mode == 1 || qp->ref_mode == 2)) {
        for (k = 0; k <= N; ++k) {
            const cs_real *Ek = qp->E + (size_t)k * nx * nz0;
            const cs_real *xrk = qp->x_ref + (size_t)k * nx;
            for (i = 0; i < nx; ++i) {
                const cs_real r =
                    (X[(size_t)k * nx + i] - xrk[i]) / qp->sc.x[i];
                qp->rref[i] = (qp->ref_mode == 2)
                    ? (r > (cs_real)0 ? (cs_real)1
                       : (r < (cs_real)0 ? (cs_real)-1 : (cs_real)0))
                    : r;
            }
            if (qp->ref_mode == 1)      /* L2 only adds Hessian curvature  */
                cs_syrk_w_acc(nz0, nx, qp->w_ref, qp->lm_w, Ek, qp->H, nz);
            cs_gemv_tw_acc(nx, nz0, qp->w_ref, qp->lm_w, Ek, qp->rref, qp->f);
        }
    }
    /* mode 3: NORMAL-only path tube. Per node, project (xi_k, eta_k) onto
     * the x_ref (xi, eta) polyline and penalize 0.5*w_ref*(|d_n|-r0)^2
     * beyond the CS_REF_TUBE_R0 dead radius, d_n = the SIGNED distance
     * along the local path normal (rank-1 Gauss-Newton in position space;
     * v, psi, a, th, T untouched). The index-matched modes above re-anchor
     * with a rigid along-path shift (tau0 re-projection slides the
     * reference with the vehicle -> the flare-ratchet's xi translation is
     * a null direction of their pull); the tube residual is geometric, so
     * tangential progress and timing stay free while normal flare
     * migration is priced in metres. Vertex-clamped projections keep the
     * SEGMENT normal, so overshooting the path's END tangentially is also
     * free (n'(p-c) discards the tangential part). sqrt-free: the
     * unnormalized segment normal (-dy, dx) is folded into alpha = w/L2.
     * Applied over the WHOLE horizon — see the CS_REF_TUBE_R0 comment for
     * why the dead-radius hinge (not a progress cap, not an index-pull
     * splice) is the composition that survived the sentinel flights. */
    if (qp->w_ref > (cs_real)0 && qp->ref_mode == 3) {
        for (k = 0; k <= N; ++k) {
            {
            const cs_real *Ek = qp->E + (size_t)k * nx * nz0;
            const cs_real px = X[(size_t)k * nx + 2];
            const cs_real py = X[(size_t)k * nx + 3];
            cs_real best_d2 = (cs_real)0, bnx = (cs_real)0,
                    bny = (cs_real)0, br = (cs_real)0, bL2 = (cs_real)0;
            int found = 0, j;
            for (j = 0; j < N; ++j) {
                const cs_real ax = qp->x_ref[(size_t)j * nx + 2];
                const cs_real ay = qp->x_ref[(size_t)j * nx + 3];
                const cs_real dxs = qp->x_ref[(size_t)(j + 1) * nx + 2] - ax;
                const cs_real dys = qp->x_ref[(size_t)(j + 1) * nx + 3] - ay;
                const cs_real L2 = dxs * dxs + dys * dys;
                cs_real t, cx, cy, ex, ey, d2;
                if (!(L2 > (cs_real)1e-9))       /* degenerate segment */
                    continue;
                t = ((px - ax) * dxs + (py - ay) * dys) / L2;
                if (t < (cs_real)0)
                    t = (cs_real)0;
                else if (t > (cs_real)1)
                    t = (cs_real)1;
                cx = ax + t * dxs;
                cy = ay + t * dys;
                ex = px - cx;
                ey = py - cy;
                d2 = ex * ex + ey * ey;
                if (!found || d2 < best_d2) {
                    found = 1;
                    best_d2 = d2;
                    bnx = -dys;              /* unnormalized normal */
                    bny = dxs;
                    br = bnx * ex + bny * ey;   /* = L * d_n (signed) */
                    bL2 = L2;
                }
            }
            if (found &&
                best_d2 > (cs_real)(CS_REF_TUBE_R0 * CS_REF_TUBE_R0)) {
                /* hinge beyond the dead radius: br holds L*d_n, so the
                 * shrunk residual L*(|d_n|-r0)*sign(d_n) = br*(1-r0/|d_n|);
                 * the Gauss-Newton curvature of 0.5w(|d|-r0)^2 is w — the
                 * full rank-1 H stays. */
                const cs_real dn = (cs_real)sqrt((double)best_d2);
                const cs_real shr =
                    (cs_real)1 - (cs_real)CS_REF_TUBE_R0 / dn;
                const cs_real alpha = qp->w_ref / bL2;
                const cs_real s2 = qp->sc.x[2], s3 = qp->sc.x[3];
                cs_real *v = qp->vref;
                int j1, j2;
                for (j = 0; j < nz0; ++j)   /* physical d(n~'p)/dz row */
                    v[j] = bnx * s2 * Ek[(size_t)j * nx + 2]
                         + bny * s3 * Ek[(size_t)j * nx + 3];
                for (j2 = 0; j2 < nz0; ++j2) {
                    const cs_real a2 = alpha * v[j2];
                    if (a2 != (cs_real)0) {
                        for (j1 = 0; j1 <= j2; ++j1) {
                            const cs_real hv = a2 * v[j1];
                            qp->H[(size_t)j2 * nz + j1] += hv;
                            if (j1 != j2)
                                qp->H[(size_t)j1 * nz + j2] += hv;
                        }
                    }
                    qp->f[j2] += a2 * (br * shr);
                }
            }
            }
        }
    }
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
        gate_bounds(qp, 0, X[3], lo, hi);
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
        gate_bounds(qp, k, Xk[3], lo, hi);
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
                /* crop bulge cap: UPPER side of the crossing row; shares
                 * the outbound block slack (coeff -1: A z - s <= hi[4]) —
                 * c_out is otherwise unused on non-outbound stages, and a
                 * shared slack legally absorbs the block's worst residual
                 * across rows (M1.2). May coexist with the lower-side
                 * c_cross on the same row. */
                if (i == 4 && hi[4] < gth)
                    row[c_out] = (cs_real)-1;
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
            gate_bounds(qp, N, (X + (size_t)N * nx)[3], lo, hi);
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
