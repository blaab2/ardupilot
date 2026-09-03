/* bs_mission_build — see bs_mission_build.h. */
#include "bs_mission_build.h"

#include <math.h>
#include <string.h>

#include "bs_dare.h"

/* The mission-independent constants.  On the firmware target this is the
 * same header bs_solver.h consumed (BS_DATA_HEADER); on host builds the
 * flash model header was included instead, so pull the runtime constants
 * explicitly when they are absent. */
#ifndef BS_TABLES_RUNTIME
#include "../model/generated/bs_model_runtime_const.h"
#endif

#define MB_MAX_SEAM 64
#define MB_PI 3.14159265358979323846
/* cornering budget: how much of the corridor the rounding arc may use,
 * and the share of the lateral accel face it may ride (the rest is the
 * closed loop's tracking margin) */
/* d = 2.0 is the corridor face itself: the honest ceiling, because
 * beyond it the law promises an arc the barrier will not allow and the
 * closed loop starts violating (measured: harsh-angle missions break
 * first, at d = 2.6-3.2).  kappa = 1.1 lets the corner ride slightly
 * past the LATERAL face because the anisotropic polygon's braking
 * direction (5.609) is what a decelerating corner actually loads.
 * Measured at these values: zero face violations and zero re-timings on
 * the SN77 waypoints and on a mixed-angle synthetic at 3-11.8 m/s. */
#ifndef MB_CORNER_D
#define MB_CORNER_D 2.0
#endif
#ifndef MB_CORNER_KAPPA
#define MB_CORNER_KAPPA 1.1
#endif
#define MB_A_LAT 5.203

/* HOST CALIBRATION HOOK.  Zero on the target (never written); the host
 * calibration driver sets it to bisect the admissible crossing pace per
 * junction angle in the CLOSED LOOP, which is what produces the
 * calibrated table this law is then checked against. */
double bs_mb_vj_cap = 0.0;

/* ------------------------------------------------------------------ util */
static double wrap180(double a)
{
    while (a > 180.0) a -= 360.0;
    while (a < -180.0) a += 360.0;
    return a;
}

static double v_adm_interp(double chi_deg)
{
    const double c = fabs(chi_deg);
    if (c <= bs_rc_adm_chi[0]) return bs_rc_adm_v[0];
    for (int i = 1; i < BS_NADM; ++i) {
        if (c <= bs_rc_adm_chi[i]) {
            const double t = (c - bs_rc_adm_chi[i - 1])
                           / (bs_rc_adm_chi[i] - bs_rc_adm_chi[i - 1]);
            return bs_rc_adm_v[i - 1]
                 + t * (bs_rc_adm_v[i] - bs_rc_adm_v[i - 1]);
        }
    }
    return bs_rc_adm_v[BS_NADM - 1];
}

/* The cold-injection admissible speed (slow_junction_resim V7) is the
 * right budget for an INGRESS junction — a state handed to the chart
 * with no preparation — and the WRONG one mid-mission, where the closed
 * loop pre-banks into the seam.  Exported for the ingress builder;
 * mission seams use the corridor-radius law below. */
double bs_mb_v_adm(double chi_deg) { return v_adm_interp(chi_deg); }

/* TRIED AND REJECTED (2026-09-03): holding the corner pace across the
 * rounding span T = r tan(chi/2) either side of each vertex, so the
 * tracker would finish decelerating before it rounds.  Measured on the
 * SN77 waypoints it cost 38 ticks (1072 -> 1110) and introduced a face
 * violation -- about 1.5 s per corner spent at corner pace on straight
 * ground, for a tracking benefit the model cannot see.  The profile
 * therefore brakes to the corner pace AT the vertex. */

/* junction crossing speed for a seam of angle chi (deg): the kink-
 * rounding law fitted in flight, capped by the state-level admissible
 * speed and floored. */
static double v_junction(double chi_deg, double v_trim)
{
    /* THE CORRIDOR IS THE RADIUS.  An arc tangent to both legs deviates
     * from the vertex by d = r (sec(chi/2) - 1), so the corridor budget
     * d_eff the barrier already allows IS a cornering radius
     *     r = d_eff / (sec(chi/2) - 1),
     * and the fastest crossing the chart's own faces admit is
     *     v = sqrt(kappa a_lat r).
     * This is not a new constraint on the optimizer: it is the frame pace
     * at which the optimizer's OWN constraint set is feasible, which it
     * cannot choose itself (the lag band caps how far the published track
     * may fall behind the frame).  At the deployed faces it reproduces
     * the hand-tuned SN77 window pace: 6.6 m/s at 72 deg. */
    const double half = fabs(chi_deg) * MB_PI / 360.0;
    double v = v_trim;
    /* The geometric ceiling.  THE CORRIDOR IS DISTANCE TO THE NEAREST
     * LEG, not to the vertex: an inscribed arc of radius r departs the
     * legs by r(1 - cos(chi/2)), while its distance from the VERTEX is
     * the much larger r(sec(chi/2) - 1).  Using the vertex form (as an
     * earlier revision did) understates the admissible radius by 70 % at
     * 108 deg and priced that corner at 3.85 m/s where the flown record
     * achieves 5.51. */
    if (half > 1e-4) {
        const double one_m_cos = 1.0 - cos(half);
        if (one_m_cos > 1e-9) {
            const double r = MB_CORNER_D / one_m_cos;
            const double vc = sqrt(MB_CORNER_KAPPA * MB_A_LAT * r);
            if (vc < v) v = vc;
        }
    }
    /* the CLOSED-LOOP calibration governs: the law prices the rounding
     * arc but not the along/cross coupling, and is optimistic in the
     * mid range (measured: 5.26 vs 3.78 m/s admissible at 90 deg). */
    {
        const double c = fabs(chi_deg);
        double vs;
        if (c <= bs_rc_seam_chi[0]) {
            vs = bs_rc_seam_v[0];
        } else if (c >= bs_rc_seam_chi[BS_NSEAM - 1]) {
            vs = bs_rc_seam_v[BS_NSEAM - 1];
        } else {
            int i = 1;
            while (i < BS_NSEAM - 1 && c > bs_rc_seam_chi[i]) i++;
            const double t = (c - bs_rc_seam_chi[i - 1])
                           / (bs_rc_seam_chi[i] - bs_rc_seam_chi[i - 1]);
            vs = bs_rc_seam_v[i - 1]
               + t * (bs_rc_seam_v[i] - bs_rc_seam_v[i - 1]);
        }
        if (vs < v) v = vs;
    }
    if (bs_mb_vj_cap > 0.0) v = bs_mb_vj_cap;   /* calibration override */
    if (v < BS_VJ_FLOOR) v = BS_VJ_FLOOR;
    if (v > v_trim) v = v_trim;
    return v;
}

/* -------------------------------------------------------------- families */
/* family t: rows_trim_aniso(m = 0.1) with the runtime VLEG; family h:
 * rows_hover floors.  Emission order is the constructors' order, which
 * the flash generator mirrors. */
static void build_rows_t(double *rows, double v_leg)
{
    int n = 0;
    for (int k = 0; k < 20; ++k) {
        const double th = 18.0 * (double)k * MB_PI / 180.0;
        const double c = cos(th), s = sin(th);
        double *r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[0] = c; r[4] = s;
        r[9] = v_leg - (v_leg - BS_M_NU_FWD) * c;
        n++;
        r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[1] = c; r[5] = s;
        r[9] = bs_rc_aniso[k];
        n++;
    }
    for (int sg = 0; sg < 2; ++sg) {
        const double sgn = (sg == 0) ? 1.0 : -1.0;
        double *r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[2] = sgn; r[9] = BS_M_LAG; n++;
        r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[3] = sgn; r[9] = BS_M_CORR; n++;
        r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[6] = sgn; r[9] = BS_JBOX; n++;
        r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[8] = sgn; r[9] = BS_JBOX; n++;
    }
    double *r = &rows[(size_t)n * 10];
    memset(r, 0, sizeof(double) * 10);
    r[7] = 1.0; r[9] = BS_M_NU_FWD; n++;
    r = &rows[(size_t)n * 10];
    memset(r, 0, sizeof(double) * 10);
    r[7] = -1.0; r[9] = v_leg - BS_M_NU_FWD; n++;
    /* n == 50 == BS_NROW by construction */
}

static void build_rows_h(double *rows)
{
    int n = 0;
    for (int k = 0; k < 20; ++k) {
        const double th = 18.0 * (double)k * MB_PI / 180.0;
        const double c = cos(th), s = sin(th);
        double *r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[0] = c; r[4] = s; r[9] = BS_HOV_SPEED; n++;
        r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[1] = c; r[5] = s; r[9] = BS_HOV_ACCEL; n++;
    }
    for (int sg = 0; sg < 2; ++sg) {
        const double sgn = (sg == 0) ? 1.0 : -1.0;
        double *r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[2] = sgn; r[9] = BS_HOV_LAG; n++;
        r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[3] = sgn; r[9] = BS_HOV_CORR; n++;
        r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[6] = sgn; r[9] = BS_HOV_JERK; n++;
        r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[8] = sgn; r[9] = BS_HOV_JERK; n++;
        r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[7] = sgn; r[9] = BS_HOV_NU; n++;
    }
    /* n == 50 == BS_NROW */
}

/* the terminal pair of the family's rows: dare_pair of record. */
static int family_dare(const double *rows, double *P6, double *K36)
{
    double Mx[36], Mu[9];
    memset(Mx, 0, sizeof(Mx));
    memset(Mu, 0, sizeof(Mu));
    for (int i = 0; i < BS_NROW; ++i) {
        const double *r = &rows[(size_t)i * 10];
        const double w = 2.0 * BS_EPS / (r[9] * r[9]);
        for (int a = 0; a < 6; ++a) {
            for (int b = 0; b < 6; ++b) {
                Mx[a * 6 + b] += w * r[a] * r[b];
            }
        }
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                Mu[a * 3 + b] += w * r[6 + a] * r[6 + b];
            }
        }
    }
    double QL[9] = { BS_W_QD, 0, 0, 0, BS_W_QA, 0, 0, 0, BS_W_RL };
    double QT[9] = { BS_W_RC, 0, 0, 0, BS_W_REG, 0, 0, 0, BS_W_REG };
    double RLm[4] = { BS_W_RJ, 0, 0, BS_W_RNU };
    double RTm[1] = { BS_W_RJ };
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            QL[a * 3 + b] += Mx[a * 6 + b];
            QT[a * 3 + b] += Mx[(a + 3) * 6 + (b + 3)];
        }
    }
    RLm[0] += Mu[0]; RLm[1] += Mu[1]; RLm[2] += Mu[3]; RLm[3] += Mu[4];
    RTm[0] += Mu[8];

    double PL[9], KL[6], PT[9], KT[3];
    if (bs_dare_solve(bs_rc_AL, bs_rc_BL, QL, RLm, 3, 2, PL, KL,
                      5000, 1e-13) != 0) return 1;
    if (bs_dare_solve(bs_rc_AT, bs_rc_BT, QT, RTm, 3, 1, PT, KT,
                      5000, 1e-13) != 0) return 1;
    memset(P6, 0, sizeof(double) * 36);
    memset(K36, 0, sizeof(double) * 18);
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            P6[a * 6 + b] = PL[a * 3 + b];
            P6[(a + 3) * 6 + (b + 3)] = PT[a * 3 + b];
        }
    }
    for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 3; ++b) {
            K36[a * 6 + b] = KL[a * 3 + b];
        }
    }
    for (int b = 0; b < 3; ++b) {
        K36[2 * 6 + 3 + b] = KT[b];
    }
    return 0;
}

/* T_switch(chi) in the core's layout: pairs (0,4), (1,5), (2,3). */
static void seam_rot(double *T, double chi_deg)
{
    const double c = cos(chi_deg * MB_PI / 180.0);
    const double s = sin(chi_deg * MB_PI / 180.0);
    static const int pa[3] = { 0, 1, 2 }, pb[3] = { 4, 5, 3 };
    memset(T, 0, sizeof(double) * 36);
    for (int k = 0; k < 3; ++k) {
        T[pa[k] * 6 + pa[k]] = c;
        T[pa[k] * 6 + pb[k]] = s;
        T[pb[k] * 6 + pa[k]] = -s;
        T[pb[k] * 6 + pb[k]] = c;
    }
}

static void seam_off(double *c6, double chi_deg, double v_trim)
{
    const double c = cos(chi_deg * MB_PI / 180.0);
    const double s = sin(chi_deg * MB_PI / 180.0);
    memset(c6, 0, sizeof(double) * 6);
    c6[0] = v_trim * (c - 1.0);
    c6[4] = -v_trim * s;
}

/* ---------------------------------------------------------------- sizing */
size_t bs_mission_size(int n_wp, double len_m, double v_cap_ms)
{
    double v_trim = v_cap_ms * BS_COS_PI20 - BS_VTRIM_STANDOFF;
    if (v_trim < 1.0) v_trim = 1.0;
    /* ramps cost ~44 extra ticks per breakpoint; cruise covers CELL_T
     * per tick.  Generous by design: this sizes an allocation, the
     * build bounds-checks against it. */
    int n_est = (int)(len_m / (BS_TS * v_trim * 0.35)) + 44 * (n_wp + 2)
              + BS_MB_HOLD + BS_MB_PAD + 64;
    if (n_est > BS_MB_MAX_TICKS + BS_MB_HOLD + BS_MB_PAD) {
        n_est = BS_MB_MAX_TICKS + BS_MB_HOLD + BS_MB_PAD;
    }
    return (size_t)n_est * (5 * sizeof(double) + 1 + 3 * sizeof(int))
         + (size_t)(2 * BS_NROW * 10 + 2 * 36 + 2 * 18) * sizeof(double)
         + (size_t)MB_MAX_SEAM * (36 + 6) * sizeof(double)
         + 1024;
}

/* ----------------------------------------------------------------- build */
bs_mb_status bs_mission_build(const double *vx, const double *vy, int n_wp,
                              const bs_mission_params *pp,
                              void *mem, size_t mem_len,
                              bs_mission_tables *out,
                              bs_mission_report *rep)
{
    memset(rep, 0, sizeof(*rep));
#define FAIL(st, g, d) do { rep->status = (st); rep->gate = (g); \
                            rep->detail = (d); return (st); } while (0)

    if (n_wp < 2 || n_wp > BS_MB_MAX_WP) FAIL(BS_MB_ERR_NWP, 0, n_wp);
    /* Tolerance on the ceiling: the configured speed arrives as a
     * single-precision parameter, and 11.8f promotes to 11.80000019 --
     * one ulp over the double literal.  Refusing the airframe's own
     * cruise on that is a boundary bug, measured once. */
    if (!(pp->v_cap_ms >= BS_VCAP_MIN - 1e-3 &&
          pp->v_cap_ms <= BS_VCAP_MAX + 1e-3)) {
        FAIL(BS_MB_ERR_SPEED, 6, pp->v_cap_ms);
    }

    /* speed system */
    const double v_leg = pp->v_cap_ms * BS_COS_PI20;
    const double v_trim = v_leg - BS_VTRIM_STANDOFF;
    const double cell_t = BS_TS * v_trim;
    /* The re-timing threshold scales with the cell, FLOORED at 1.8 m:
     * the natural e_l excursions (junction transients, engage sprint)
     * sit near 1.4 m independent of speed, and an unfloored threshold
     * dives below them at low cruise (measured: 18 spurious retards on
     * a 4 m/s synthetic mission, 0 with the floor).  At the record's
     * cell the floor is inactive (2.106 > 1.8), so the flown config is
     * unchanged. */
    double hyst = BS_HYST_RATIO * cell_t;
    if (hyst < 1.8) hyst = 1.8;
    const double cell_min = hyst / 1.5;
    if (v_trim < BS_NU_MIN + 0.4) FAIL(BS_MB_ERR_SPEED, 6, v_trim);

    /* ---- vertex cleanup: merge short legs, drop straight-through ---- */
    double px[BS_MB_MAX_WP], py[BS_MB_MAX_WP];
    int src[BS_MB_MAX_WP];          /* cleaned index -> input index */
    int n_v = 0;
    for (int i = 0; i < n_wp; ++i) {
        if (n_v > 0) {
            const double dx = vx[i] - px[n_v - 1];
            const double dy = vy[i] - py[n_v - 1];
            if (sqrt(dx * dx + dy * dy) < 1.0) {
                continue;                     /* merge into predecessor */
            }
        }
        px[n_v] = vx[i];
        py[n_v] = vy[i];
        src[n_v] = i;
        n_v++;
    }
    /* drop straight-through vertices (tiny angle AND tiny chord miss) */
    for (int pass = 0; pass < BS_MB_MAX_WP; ++pass) {
        int dropped = 0;
        for (int i = 1; i + 1 < n_v; ++i) {
            const double h0 = atan2(py[i] - py[i - 1], px[i] - px[i - 1]);
            const double h1 = atan2(py[i + 1] - py[i], px[i + 1] - px[i]);
            const double chi = wrap180((h1 - h0) * 180.0 / MB_PI);
            if (fabs(chi) >= 0.5) continue;
            const double l1 = hypot(px[i] - px[i - 1], py[i] - py[i - 1]);
            const double l2 = hypot(px[i + 1] - px[i], py[i + 1] - py[i]);
            const double miss = (l1 * l2 / (l1 + l2))
                              * fabs(sin(chi * MB_PI / 180.0));
            if (miss > 0.5) continue;
            for (int j = i; j + 1 < n_v; ++j) {
                px[j] = px[j + 1];
                py[j] = py[j + 1];
                src[j] = src[j + 1];
            }
            n_v--;
            dropped = 1;
            break;
        }
        if (!dropped) break;
    }
    if (n_v < 2) FAIL(BS_MB_ERR_GEOM, 0, n_v);
    if (n_v - 2 >= MB_MAX_SEAM) FAIL(BS_MB_ERR_GEOM, 0, n_v);

    /* legs and junction angles */
    double leg_len[BS_MB_MAX_WP], leg_hd[BS_MB_MAX_WP];
    double chi[BS_MB_MAX_WP];       /* chi[i] at cleaned vertex i (1..) */
    double total = 0.0;
    for (int i = 0; i + 1 < n_v; ++i) {
        leg_len[i] = hypot(px[i + 1] - px[i], py[i + 1] - py[i]);
        leg_hd[i] = atan2(py[i + 1] - py[i], px[i + 1] - px[i]);
        total += leg_len[i];
    }
    for (int i = 1; i + 1 < n_v; ++i) {
        chi[i] = wrap180((leg_hd[i] - leg_hd[i - 1]) * 180.0 / MB_PI);
    }

    /* ---- breakpoint speeds ---- */
    double vb[BS_MB_MAX_WP + 1];    /* pace at each cleaned vertex */
    vb[0] = 0.0;                    /* leg 1 published from rest */
    for (int i = 1; i + 1 < n_v; ++i) {
        vb[i] = v_junction(chi[i], v_trim);
        /* THE ARC MUST FIT BETWEEN ITS NEIGHBOURS.  The rounding arc is
         * tangent at T = r tan(chi/2) either side of the vertex, so with
         * adjacent legs of length L the radius is bounded by
         * r <= (L/2) / tan(chi/2) -- the half-leg is shared with the
         * neighbouring corner.  The seam-pace table is calibrated on
         * ISOLATED corners (140 m legs); without this term a mission
         * whose corners sit 15 m apart inherits an isolated corner's
         * pace and breaks (measured: 0.149 of face violation on the
         * mixed-angle synthetic). */
        const double half = fabs(chi[i]) * MB_PI / 360.0;
        const double tn = tan(half);
        if (tn > 1e-6) {
            double Lmin = leg_len[i - 1];
            if (leg_len[i] < Lmin) Lmin = leg_len[i];
            const double r_fit = 0.5 * Lmin / tn;
            const double v_fit = sqrt(MB_CORNER_KAPPA * MB_A_LAT * r_fit);
            if (v_fit < vb[i]) vb[i] = v_fit;
            if (vb[i] < BS_VJ_FLOOR) vb[i] = BS_VJ_FLOOR;
        }
    }
    /* hover splice: r_h fixed point on the final leg */
    const int last = n_v - 2;       /* index of the final leg */
    double r_h = 1.0;
    for (int it = 0; it < 6; ++it) {
        const double d_h = -bs_rc_hov_x[0] * r_h;    /* 20*r_h */
        double avail = leg_len[last] - d_h;
        if (avail < 0.0) avail = 0.0;
        const double v_in = (last == 0) ? vb[0] : vb[last];
        double v_reach = sqrt(v_in * v_in + 2.0 * BS_ARC_A_EFF * avail);
        if (v_reach > v_trim) v_reach = v_trim;
        double r_new = v_reach / BS_HOV_PACE_REF;
        if (r_new > 1.0) r_new = 1.0;
        r_h = r_new;
    }
    const double v_entry = r_h * BS_HOV_PACE_REF;
    if (v_entry < 1.0) FAIL(BS_MB_ERR_HOVER, 8, v_entry);
    vb[n_v - 1] = v_entry;          /* target at the splice breakpoint */
    const double d_hov = -bs_rc_hov_x[0] * r_h;

    /* forward + backward reachability over breakpoints (splice-aware:
     * the final leg's cruise portion ends d_hov early) */
    for (int i = 1; i < n_v; ++i) {
        const double L = leg_len[i - 1] - ((i == n_v - 1) ? d_hov : 0.0);
        const double vmax = sqrt(vb[i - 1] * vb[i - 1]
                                 + 2.0 * BS_ARC_A_EFF * (L > 0 ? L : 0));
        if (vb[i] > vmax) vb[i] = vmax;
    }
    for (int i = n_v - 2; i >= 0; --i) {
        const double L = leg_len[i] - ((i == n_v - 2) ? d_hov : 0.0);
        const double vmax = sqrt(vb[i + 1] * vb[i + 1]
                                 + 2.0 * BS_ARC_A_EFF * (L > 0 ? L : 0));
        if (vb[i] > vmax) vb[i] = vmax;
    }

    /* ---- allocate the block ---- */
    const int cap_ticks = BS_MB_MAX_TICKS + BS_MB_HOLD + BS_MB_PAD;
    /* carve the block in DESCENDING alignment order (doubles, ints,
     * chars) so every pointer is naturally aligned; the integer
     * round-trip keeps -Wcast-align=strict quiet about the (provably
     * aligned) conversions. */
    unsigned long long cur = (unsigned long long)(size_t)mem;
    size_t need = 0;
#define TAKE(ptr, ty, count) do { \
        need += sizeof(ty) * (size_t)(count); \
        if (need > mem_len) FAIL(BS_MB_ERR_MEM, 9, (double)need); \
        (ptr) = (ty *)(size_t)cur; \
        cur += sizeof(ty) * (size_t)(count); \
    } while (0)

    memset(out, 0, sizeof(*out));
    double *path, *psi, *angt, *rows, *P, *K, *rot, *off;
    signed char *fam;
    int *sfam, *srot, *soff;
    /* worst-case per-tick arrays are sized after counting; count first by
     * integrating the profile without writes. */

    /* ---- pass 1: integrate the published profile, count ticks ---- */
    int g_tick[BS_MB_MAX_WP + 1];   /* tick of each cleaned vertex */
    int n_tick = 0;                 /* ticks BEFORE the hover splice */
    {
        double v = vb[0], a = 0.0;
        g_tick[0] = 0;
        for (int i = 0; i + 1 < n_v; ++i) {
            const double L = leg_len[i] - ((i == n_v - 2) ? d_hov : 0.0);
            const double vt = vb[i + 1];
            double arc = 0.0;
            int nt = 0;
            while (arc < L && nt < cap_ticks) {
                double an = a + BS_ARC_J_EFF * BS_TS;
                if (an > BS_ARC_A_EFF) an = BS_ARC_A_EFF;
                double vn = v + 0.5 * (a + an) * BS_TS;
                if (vn > v_trim) { vn = v_trim; an = 0.0; }
                const double rem_b = L - arc;
                const double rem = (rem_b > 0.0) ? rem_b : 0.0;
                const double vbrk = sqrt(vt * vt
                                         + 2.0 * BS_ARC_A_EFF * rem);
                if (vn > vbrk) { vn = vbrk; an = 0.0; }
                arc += 0.5 * (v + vn) * BS_TS;
                v = vn;
                a = an;
                nt++;
                if (arc <= 0.0 && nt > 8) {
                    FAIL(BS_MB_ERR_GEOM, 2, (double)i);  /* stalled */
                }
            }
            if (arc < L) FAIL(BS_MB_ERR_TICKS, 9, (double)n_tick + nt);
            n_tick += nt;
            g_tick[i + 1] = n_tick;    /* vertex i+1 (or splice) tick */
            if (n_tick > BS_MB_MAX_TICKS) {
                FAIL(BS_MB_ERR_TICKS, 9, (double)n_tick);
            }
        }
    }
    const int node = n_tick + BS_HOV_TICKS;
    const int n_end = node + BS_MB_HOLD;
    const int n_clk = n_end + BS_MB_PAD;
    const int n_path = node + 1;
    const int hov_in = node - BS_HOV_TICKS;
    if (node > BS_MB_MAX_TICKS) FAIL(BS_MB_ERR_TICKS, 9, (double)node);

    TAKE(path, double, 2 * n_path);
    TAKE(psi, double, n_path);
    TAKE(angt, double, n_clk);
    TAKE(rows, double, 2 * BS_NROW * 10);
    TAKE(P, double, 2 * 36);
    TAKE(K, double, 2 * 18);
    TAKE(rot, double, MB_MAX_SEAM * 36);
    TAKE(off, double, MB_MAX_SEAM * 6);
    TAKE(sfam, int, n_clk);
    TAKE(srot, int, n_clk);
    TAKE(soff, int, n_clk);
    TAKE(fam, signed char, n_clk);
#undef TAKE

    /* ---- pass 2: fill the per-tick tables ---- */
    memset(angt, 0, sizeof(double) * (size_t)n_clk);
    for (int t = 0; t < n_clk; ++t) {
        fam[t] = (signed char)((t >= hov_in) ? 1 : 0);
        sfam[t] = (t >= hov_in) ? BS_MB_FAM_H : BS_MB_FAM_T;
        srot[t] = 0;
        soff[t] = 0;
    }

    /* seam stores: slot 0 identity/zero, then one slot per distinct
     * (chi, crossing pace).  THE KICK IS AT THE PUBLISHED PACE OF THE
     * CROSSING, off_vec(v_j, v_j, chi) -- the certified tables always
     * carried the reference's pace at the seam, and a trim-pace kick
     * at a slow published crossing injects a phantom cross-rate the
     * model then has to fight (measured on the SN77 gate: 4.0 of face
     * violation at the first 72-deg seam; with the crossing-pace kick
     * the v_adm cap makes the injected |v_j sin chi| sit inside the
     * certified cold-injection basin by construction). */
    int n_seam = 1;
    memset(rot, 0, sizeof(double) * 36);
    for (int d = 0; d < 6; ++d) rot[d * 6 + d] = 1.0;
    memset(off, 0, sizeof(double) * 6);
    for (int i = 1; i + 1 < n_v; ++i) {
        double c6[6];
        seam_off(c6, chi[i], vb[i]);
        int slot = -1;
        for (int k = 1; k < n_seam; ++k) {
            if (fabs(off[k * 6 + 0] - c6[0]) < 1e-12 &&
                fabs(off[k * 6 + 4] - c6[4]) < 1e-12) {
                slot = k;
                break;
            }
        }
        if (slot < 0) {
            if (n_seam >= MB_MAX_SEAM) FAIL(BS_MB_ERR_GEOM, 5, n_seam);
            slot = n_seam++;
            seam_rot(&rot[(size_t)slot * 36], chi[i]);
            memcpy(&off[(size_t)slot * 6], c6, sizeof(c6));
        }
        const int gt = g_tick[i];
        angt[gt] = chi[i];
        srot[gt] = slot;
        soff[gt] = slot;
    }

    /* published arc + chart heading, leg by leg (re-integrated with the
     * same recursion as pass 1, positions scaled so each breakpoint
     * lands exactly on its tick) */
    {
        double v = vb[0], a = 0.0;
        int t0 = 0;
        path[0] = px[0];
        path[1] = py[0];
        for (int i = 0; i + 1 < n_v; ++i) {
            const double L = leg_len[i] - ((i == n_v - 2) ? d_hov : 0.0);
            const double vt = vb[i + 1];
            const int nt = g_tick[i + 1] - g_tick[i];
            /* integrate; the unscaled arc is STAGED in the path array
             * (no stack buffer: the solver thread's stack is 8 kB) */
            double arc = 0.0;
            for (int k = 0; k < nt; ++k) {
                double an = a + BS_ARC_J_EFF * BS_TS;
                if (an > BS_ARC_A_EFF) an = BS_ARC_A_EFF;
                double vn = v + 0.5 * (a + an) * BS_TS;
                if (vn > v_trim) { vn = v_trim; an = 0.0; }
                const double rem_b = L - arc;
                const double rem = (rem_b > 0.0) ? rem_b : 0.0;
                const double vbrk = sqrt(vt * vt
                                         + 2.0 * BS_ARC_A_EFF * rem);
                if (vn > vbrk) { vn = vbrk; an = 0.0; }
                arc += 0.5 * (v + vn) * BS_TS;
                v = vn;
                a = an;
                path[2 * (t0 + k + 1)] = arc;
            }
            const double scale = (arc > 0.0) ? (L / arc) : 0.0;
            const double ch = cos(leg_hd[i]), sh = sin(leg_hd[i]);
            for (int k = 0; k < nt; ++k) {
                const double s_al = path[2 * (t0 + k + 1)] * scale;
                path[2 * (t0 + k + 1)] = px[i] + s_al * ch;
                path[2 * (t0 + k + 1) + 1] = py[i] + s_al * sh;
            }
            /* chart heading: this leg owns ticks (t0, t0+nt]; the seam's
             * destination tick carries the NEW leg heading, so psi at
             * tick t0 (vertex i) belongs to leg i already. */
            for (int k = 0; k < nt; ++k) {
                psi[t0 + k] = leg_hd[i];
            }
            t0 += nt;
        }
        /* hover splice: the certified anchor profile, scaled by r_h,
         * along the final heading, ending at rest on the last vertex */
        const double ch = cos(leg_hd[n_v - 2]), sh = sin(leg_hd[n_v - 2]);
        for (int k = 0; k <= BS_HOV_TICKS; ++k) {
            const double s_al = bs_rc_hov_x[k] * r_h;   /* -20r .. 0 */
            path[2 * (hov_in + k)] = px[n_v - 1] + s_al * ch;
            path[2 * (hov_in + k) + 1] = py[n_v - 1] + s_al * sh;
            if (hov_in + k < n_path) psi[hov_in + k] = leg_hd[n_v - 2];
        }
        psi[node] = leg_hd[n_v - 2];
    }

    /* ---- families ---- */
    build_rows_t(&rows[0], v_leg);
    build_rows_h(&rows[(size_t)BS_NROW * 10]);
    if (family_dare(&rows[0], &P[0], &K[0]) != 0) {
        FAIL(BS_MB_ERR_DARE, 10, 0.0);
    }
    if (family_dare(&rows[(size_t)BS_NROW * 10], &P[36], &K[18]) != 0) {
        FAIL(BS_MB_ERR_DARE, 10, 1.0);
    }
    /* family h cross-check against the emitted reference: a free sanity
     * gate on the on-target DARE solver. */
    {
        double dmax = 0.0;
        for (int i = 0; i < 36; ++i) {
            const double d = fabs(P[36 + i] - bs_rc_Ph[i]);
            if (d > dmax) dmax = d;
        }
        for (int i = 0; i < 18; ++i) {
            const double d = fabs(K[18 + i] - bs_rc_Kh[i]);
            if (d > dmax) dmax = d;
        }
        if (dmax > 1e-7) FAIL(BS_MB_ERR_DARE, 10, dmax);
    }

    /* ---- acceptance gates ---- */
    /* G1 no NaN; G2 monotone arc with bounded step; G4 on-polyline */
    for (int t = 0; t < n_path; ++t) {
        if (path[2 * t] != path[2 * t] || path[2 * t + 1] != path[2 * t + 1]
            || psi[t] != psi[t]) {
            FAIL(BS_MB_ERR_GATE, 1, (double)t);
        }
    }
    for (int t = 1; t < n_path; ++t) {
        const double st = hypot(path[2 * t] - path[2 * (t - 1)],
                                path[2 * t + 1] - path[2 * (t - 1) + 1]);
        if (st > cell_t * (1.0 + 1e-9) + 1e-12) {
            FAIL(BS_MB_ERR_GATE, 2, st);
        }
    }
    /* G4: every cleaned vertex lands exactly on its tick */
    for (int i = 0; i < n_v; ++i) {
        const int gt = (i == n_v - 1) ? node : g_tick[i];
        const double dv = hypot(path[2 * gt] - px[i],
                                path[2 * gt + 1] - py[i]);
        if (dv > 1e-6) FAIL(BS_MB_ERR_GATE, 4, dv);
    }
    /* G6: published pace below the face */
    for (int t = 1; t < n_path; ++t) {
        const double st = hypot(path[2 * t] - path[2 * (t - 1)],
                                path[2 * t + 1] - path[2 * (t - 1) + 1]);
        if (st / BS_TS > v_leg - 0.05) FAIL(BS_MB_ERR_GATE, 6, st / BS_TS);
    }

    /* ---- publish ---- */
    out->n_path = n_path;
    out->n_clk = n_clk;
    out->node = node;
    out->n_end = n_end;
    out->hov_in = hov_in;
    out->n_seam = n_seam;
    out->v_leg = v_leg;
    out->v_trim = v_trim;
    out->cell_t = cell_t;
    out->cell_min = cell_min;
    out->hyst = hyst;
    out->v0 = vb[0];
    out->length_m = total;
    out->path = path;
    out->psi = psi;
    out->ang = angt;
    out->fam = fam;
    out->sfam = sfam;
    out->srot = srot;
    out->soff = soff;
    out->rows = rows;
    out->P = P;
    out->K = K;
    out->rot = rot;
    out->off = off;
    out->n_wp_in = n_wp;
    for (int i = 0; i < BS_MB_MAX_WP; ++i) out->wp_tick[i] = -1;
    {
        /* input waypoint -> vertex tick; merged/dropped inputs inherit
         * the surviving predecessor's tick, and the final waypoint maps
         * to `node` (the hover terminal AT that waypoint). */
        int ci = 0;
        for (int i = 0; i < n_wp; ++i) {
            while (ci + 1 < n_v && src[ci + 1] <= i) ci++;
            out->wp_tick[i] = (ci == n_v - 1) ? node : g_tick[ci];
        }
        out->wp_tick[n_wp - 1] = node;
    }

    out->sched.family = sfam;
    out->sched.rotation = srot;
    out->sched.offset = soff;
    out->sched.length = n_clk;
    out->sched.periodic = 0;
    out->sched.rows_extra = NULL;
    out->sched.rot_extra = NULL;
    out->sched.off_extra = NULL;
    out->sched.rows_tab = rows;
    out->sched.P_tab = P;
    out->sched.K_tab = K;
    out->sched.rot_tab = rot;
    out->sched.off_tab = off;

    rep->status = BS_MB_OK;
    rep->n_ticks = n_clk;
    rep->n_seam = n_seam;
    rep->n_vert = n_v;
    rep->build_bytes = (double)need;
    return BS_MB_OK;
#undef FAIL
}
