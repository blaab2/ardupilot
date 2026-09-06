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
/* Lag-row half-width.  2 m is the deployed value; the host study opens
 * it so the barrier stops fighting a vehicle that slows for a corner
 * (the band, not the speed gain, is what pins the corner speed --
 * measured: worst face 3.36 -> 1.04 as the band goes 2 -> 20 m). */
/* 2 m is the flown value and it is the right one: a 5 m band was tried
 * (with the threshold tracking it, so the ledger did not claw the
 * freedom back) and the host saw it as neutral -- same 1072 ticks, same
 * zero violations, corridor peaks marginally lower -- but IN FLIGHT it
 * cost 16 s (288.9 s against 272.9) and the corridor got slightly worse
 * (2.09 against 2.06).  The extra freedom is spent letting the plan sit
 * further behind the schedule, and the tracker follows it there. */
/* Published cruise standoff below the speed face: V_TRIM = VLEG - this.
 * Smaller = a faster plan.  The record value is 0.2888; shrinking it
 * lengthens every cell (CELL_T = Ts*V_TRIM) so the mission takes fewer
 * ticks, at the cost of transient recovery authority -- the tracker
 * chases a reference that sits closer to the airframe's ceiling and is
 * expected to trail it slightly (owner decision 2026-09-04: valid). */
#ifndef MB_TRIM_STANDOFF
#define MB_TRIM_STANDOFF 0.15   /* was 0.2888; host: -9 ticks, clean */
#endif

#ifndef MB_LAG_BAND
#define MB_LAG_BAND 2.0
#endif
/* Absolute re-timing threshold override (0 = the scaled default).  It
 * belongs with the lag band: widening the band alone hands the loop
 * freedom the hysteresis takes straight back, because the ledger starts
 * retarding at 1.8-2.1 m and every retard costs a tick. */
#ifndef MB_HYST_ABS
#define MB_HYST_ABS 0.0
#endif
/* 1 = publish a uniform trim arc and let the loop find the corner
 * speed itself (the band + ledger doing the pacing). */
/* Above this junction angle the measured table governs.  DEFAULT 0 =
 * everywhere, and that is not timidity: confining it to the reversal
 * regime recovers 104 ticks on SN77 (1072 -> 968) with no violation
 * there, but every relaxation tried breaks a DIFFERENT cell of the
 * (mission, speed) sweep -- threshold 120 breaks the mixed-angle
 * synthetic at 11.8 m/s by 0.31, and pulling kappa back to 0.95 fixes
 * that only to break the same mission at 7.0 m/s by 0.31.  The
 * violations move rather than vanish because what they measure is
 * corner-PAIR interaction (in the synthetic, a 75 deg vertex 35 m
 * before a 150 deg reversal), which no per-angle rule can see.  Closing
 * that ~25 s honestly needs pair-aware pacing -- a profile over the
 * whole path -- not another constant. */
#ifndef MB_CAL_CHI
#define MB_CAL_CHI 0.0
#endif
#ifndef MB_NO_CORNER_PACE
#define MB_NO_CORNER_PACE 0
#endif

/* HOST CALIBRATION HOOK.  Zero on the target (never written); the host
 * calibration driver sets it to bisect the admissible crossing pace per
 * junction angle in the CLOSED LOOP, which is what produces the
 * calibrated table this law is then checked against. */
double bs_mb_vj_cap = 0.0;
double bs_mb_vj_scale = 1.0;

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

/* DRAG-AWARE FORWARD AUTHORITY.  The airframe's net forward
 * acceleration falls with speed -- drag consumes (k/m) v^2 of the tilt
 * budget -- so a ramp integrated at constant A_EFF demands ~3x the real
 * authority near cruise, and the tracker trails it (measured: ~3 s per
 * straight, 21 s per SN77 mission, lean p95 23.8 deg riding the limit).
 * The published ACCELERATION side therefore uses
 *     a_fwd(v) = 0.8 (ACCF - (k/m) v^2),   floored at 0.15,
 * while the braking side keeps constant A_EFF: drag HELPS deceleration,
 * so that side was always feasible. */
/* DEFAULT OFF, and not from timidity -- measured in SITL: the honest
 * (drag-aware) ramp made the plan 7 s slower on SN77 and the FLOWN time
 * worse by the same 7 s, with the tracker's ~21 s overhead UNCHANGED
 * (286.8 -> 294.2 s, corridor 1.92 -> 1.91).  The tracker's deficit is
 * plan-shape-independent, so honest-but-slower plans cost flown time
 * one-for-one.  The infeasibility diagnosis was falsified; whatever the
 * ~21 s is, it is not ramp authority. */
#ifndef MB_DRAG_RAMPS
#define MB_DRAG_RAMPS 0
#endif

static double arc_a_fwd(double v)
{
#if MB_DRAG_RAMPS
    double a = 0.8 * (BS_ACCF - BS_DRAG_KM * v * v);
    if (a < 0.15) a = 0.15;
    return a;
#else
    (void)v;
    return BS_ARC_A_EFF;
#endif
}

/* forward reach with drag: dv^2/ds = 2 a_fwd(v) is linear in w = v^2,
 * so w(s) = w_eq + (w0 - w_eq) exp(-c s), w_eq = ACCF/(k/m). */
static double arc_reach_fwd(double v_in, double dist)
{
#if MB_DRAG_RAMPS
    const double w_eq = BS_ACCF / BS_DRAG_KM;
    const double c = 1.6 * BS_DRAG_KM;
    const double w = w_eq + (v_in * v_in - w_eq) * exp(-c * dist);
    return (w > 0.0) ? sqrt(w) : 0.0;
#else
    return sqrt(v_in * v_in + 2.0 * BS_ARC_A_EFF * dist);
#endif
}

#include "bs_arc_step.h"

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
    /* THE CLOSED-LOOP CALIBRATION GOVERNS THE REVERSAL REGIME ONLY.
     * It was measured on ISOLATED corners against a 1.45 m model
     * corridor target, which is stricter than "clean" (the face is at
     * 2.0), so below the reversal regime it is simply conservative:
     * applying it everywhere costs SN77 113 ticks (1072 vs 959) with no
     * violation either way.  Above ~120 deg the arc model itself breaks
     * down -- the cross-kick v sin(chi) dies away and the corner
     * becomes an along-track brake -- and there the measured table is
     * the only thing that protects a harsh mission (the mixed-angle
     * synthetic breaks by 4.01 without it). */
    if (fabs(chi_deg) >= MB_CAL_CHI) {
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
    v *= bs_mb_vj_scale;                        /* sweep knob (tests/sweep_junction_pace.py) */
#if MB_NO_CORNER_PACE
    v = v_trim;                    /* uniform arc: the loop paces itself */
#endif
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
        {
            /* forward margin capped at the standoff: the face sits at
             * trim + margin and must stay inside the drag disc VLEG */
            double mfwd = BS_M_NU_FWD;
            if (mfwd > MB_TRIM_STANDOFF) mfwd = MB_TRIM_STANDOFF;
            r[9] = v_leg - (v_leg - mfwd) * c;
        }
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
        r[2] = sgn; r[9] = MB_LAG_BAND; n++;
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
#if BS_NROW >= 60
    /* DRAG-AWARE FORWARD-ACCEL CHORDS (Tier 1) -- BLOCKED, kept for the
     * pace-carrying chart.  Measured 2026-09-04 with BS_NROW=60: the
     * closed loop BREAKS at cruise (SN77 viol 4.38, corridor 2.31;
     * mixed-angle synthetic 2.57), because the chart's speed is
     * FICTIONAL exactly where drag matters: mid-corner the published
     * arc is at 4-6 m/s but the chart still carries ~trim (the
     * chart/publish divergence, measured 38 seams: authored 4 vs chart
     * 10.6), so these rows clamp corner-recovery authority the airframe
     * physically has.  Correct math, wrong coordinate: the chords
     * become deployable exactly when the chart carries the physical
     * pace (the nu-driven-publish follow-up in the papers).  Until
     * then BS_NROW stays 50 and this block compiles out.
     *
     * The geometry (valid once the precondition holds): forward
     * authority along a direction with tangential component cos(th) is
     *     A_k(v) = aniso_k - (k/m) cos(th) v^2,
     * a CONCAVE curve in v, so chords between speed knots lie BELOW it:
     * each chord is a conservative affine row coupling accel and delta,
     * and the problem stays in the certified affine class.  Applied to
     * the five most-forward polygon directions (|th| <= 36 deg), two
     * chords each over [0, 7] and [7, VLEG+0.25].  The constant aniso
     * row remains as the v = 0 bound; braking directions are already
     * conservative (drag helps them) and stay untouched. */
    {
        const double vt = v_leg - MB_TRIM_STANDOFF;   /* chart trim */
        const double kn[3] = { 0.0, 7.0, v_leg + 0.25 };
        const int dir[5] = { 0, 1, 19, 2, 18 };
        for (int d = 0; d < 5; ++d) {
            const double th = 18.0 * (double)dir[d] * MB_PI / 180.0;
            const double ct = cos(th), st = sin(th);
            for (int c2 = 0; c2 < 2; ++c2) {
                const double va = kn[c2], vb2 = kn[c2 + 1];
                const double Aa = bs_rc_aniso[dir[d]]
                                - BS_DRAG_KM * ct * va * va;
                const double Ab = bs_rc_aniso[dir[d]]
                                - BS_DRAG_KM * ct * vb2 * vb2;
                const double sl = (Ab - Aa) / (vb2 - va);
                r = &rows[(size_t)n * 10];
                memset(r, 0, sizeof(double) * 10);
                r[1] = ct; r[5] = st;
                r[0] = -sl;                       /* -sl > 0: tighter fast */
                r[9] = Aa + sl * (vt - va);
                n++;
            }
        }
    }
#endif
    /* n == BS_NROW by construction */
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
#if BS_NROW >= 60
    /* pad to BS_NROW with inert rows (h = e = 0, m = 1): each adds a
     * constant to the objective, zero gradient, zero Mx contribution */
    while (n < BS_NROW) {
        double *r = &rows[(size_t)n * 10];
        memset(r, 0, sizeof(double) * 10);
        r[9] = 1.0; n++;
    }
#endif
    /* n == BS_NROW */
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
#ifdef BS_STAGE2_QD
    double QL[9] = { BS_STAGE2_QD, 0, 0, 0, BS_W_QA, 0, 0, 0, BS_W_RL };
#else
    double QL[9] = { BS_W_QD, 0, 0, 0, BS_W_QA, 0, 0, 0, BS_W_RL };
#endif
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
    double v_trim = v_cap_ms * BS_COS_PI20 - MB_TRIM_STANDOFF;
    if (v_trim < 1.0) v_trim = 1.0;
    /* ramps cost ~44 extra ticks per breakpoint; cruise covers CELL_T
     * per tick.  Generous by design: this sizes an allocation, the
     * build bounds-checks against it.
     *
     * Saturate BEFORE the double->int cast.  The worst-case caller
     * (AP_BSolver::init, len_m = 1e9) used to overflow that cast: the
     * result is UB, went to INT_MIN on both targets, the size_t product
     * wrapped, and on a 64-bit host the driver's reservation search
     * cycled forever (SITL hang, 2026-09-05).  The int path below is
     * byte-identical to the old one for every in-range input. */
    const int cap = BS_MB_MAX_TICKS + BS_MB_HOLD + BS_MB_PAD;
    const double q = len_m / (BS_TS * v_trim * 0.35);
    int n_est;
    if (!(q < (double)cap)) {          /* huge, inf or NaN: the cap */
        n_est = cap;
    } else {
        n_est = (int)q + 44 * (n_wp + 2) + BS_MB_HOLD + BS_MB_PAD + 64;
        if (n_est > cap) {
            n_est = cap;
        }
    }
    return (size_t)n_est * (5 * sizeof(double) + 1 + 3 * sizeof(int))
         + (size_t)(2 * BS_NROW * 10 + 2 * 36 + 2 * 18) * sizeof(double)
         + (size_t)MB_MAX_SEAM * (36 + 6) * sizeof(double)
         + 1024;
}

/* ------------------------------------------------------------------ plan */
/* Build the vertex-level plan: everything whole-mission (cleanup, speed
 * system, junction paces with the leg-fit cap, hover splice fixed
 * point, forward/backward reachability, pass-1 tick counting with the
 * integrator carry recorded at every vertex, deduped seam stores,
 * families with their DARE pairs, wp map, clock geometry).  The batch
 * builder and the streaming renderer both consume THIS object, so the
 * two fills cannot drift: single-sourcing is the structural guarantee
 * behind the stream==batch byte-equality gate (tests/stream_gate.c).
 *
 * max_ticks: BS_MB_MAX_TICKS reproduces the classic batch refusal
 * semantics exactly -- the per-leg integration guard runs to
 * max_ticks + BS_MB_HOLD + BS_MB_PAD and the total-tick and node
 * refusals fire at max_ticks -- while a streaming caller passes its
 * own, larger cap (mission length is then bounded by AP_Mission
 * storage and battery, not by RAM). */
bs_mb_status bs_mission_plan_build(const double *vx, const double *vy,
                                   int n_wp, const bs_mission_params *pp,
                                   int max_ticks, bs_mission_plan *pl,
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
    const double v_trim = v_leg - MB_TRIM_STANDOFF;
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
    /* THE RE-TIMING THRESHOLD TRACKS THE LAG BAND.  They are one
     * mechanism: the band says how far behind the loop may fall, the
     * threshold says when the reference waits for it.  Widening the
     * band alone hands back freedom the ledger immediately reclaims --
     * measured at band 5 with the 2.1 m threshold, 16 retards and 16
     * ticks lost for nothing.  At the deployed 2 m band 0.8*band is
     * below the 1.8 m floor, so the flown configuration is unchanged. */
    if (0.8 * MB_LAG_BAND > hyst) hyst = 0.8 * MB_LAG_BAND;
    /* the explicit override outranks the tracking rule (it exists for
     * experiments that place the threshold INSIDE the band) */
    if (MB_HYST_ABS > 0.0) hyst = MB_HYST_ABS;
    const double cell_min = hyst / 1.5;
    if (v_trim < BS_NU_MIN + 0.4) FAIL(BS_MB_ERR_SPEED, 6, v_trim);

    /* ---- vertex cleanup: merge short legs, drop straight-through ---- */
    double *px = pl->px, *py = pl->py;
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
    double *leg_len = pl->leg_len, *leg_hd = pl->leg_hd, *chi = pl->chi;
    for (int i = 0; i + 1 < n_v; ++i) {
        leg_len[i] = hypot(px[i + 1] - px[i], py[i + 1] - py[i]);
        leg_hd[i] = atan2(py[i + 1] - py[i], px[i + 1] - px[i]);
    }
    for (int i = 1; i + 1 < n_v; ++i) {
        chi[i] = wrap180((leg_hd[i] - leg_hd[i - 1]) * 180.0 / MB_PI);
    }

    /* ---- breakpoint speeds ---- */
    double *vb = pl->vb;            /* pace at each cleaned vertex */
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
        double v_reach = arc_reach_fwd(v_in, avail);
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
        const double vmax = arc_reach_fwd(vb[i - 1], (L > 0 ? L : 0));
        if (vb[i] > vmax) vb[i] = vmax;
    }
    for (int i = n_v - 2; i >= 0; --i) {
        const double L = leg_len[i] - ((i == n_v - 2) ? d_hov : 0.0);
        const double vmax = arc_reach_fwd(vb[i + 1], (L > 0 ? L : 0));
        if (vb[i] > vmax) vb[i] = vmax;
    }

    /* ---- pass 1: integrate the published profile, count ticks, and
     * record the integrator carry at every vertex (the renderer re-runs
     * the same recursion per leg from this carry) ---- */
    const int cap_ticks = max_ticks + BS_MB_HOLD + BS_MB_PAD;
    int *g_tick = pl->g_tick;       /* tick of each cleaned vertex */
    int n_tick = 0;                 /* ticks BEFORE the hover splice */
    {
        double v = vb[0], a = 0.0;
        g_tick[0] = 0;
        pl->carry_v[0] = v;
        pl->carry_a[0] = a;
        for (int i = 0; i + 1 < n_v; ++i) {
            const double L = leg_len[i] - ((i == n_v - 2) ? d_hov : 0.0);
            const double vt = vb[i + 1];
            double arc = 0.0;
            int nt = 0;
            while (arc < L && nt < cap_ticks) {
                bs_arc_step(&v, &a, &arc, L, vt, v_trim);
                nt++;
                if (arc <= 0.0 && nt > 8) {
                    FAIL(BS_MB_ERR_GEOM, 2, (double)i);  /* stalled */
                }
            }
            if (arc < L) FAIL(BS_MB_ERR_TICKS, 9, (double)n_tick + nt);
            n_tick += nt;
            g_tick[i + 1] = n_tick;    /* vertex i+1 (or splice) tick */
            pl->carry_v[i + 1] = v;
            pl->carry_a[i + 1] = a;
            if (n_tick > max_ticks) {
                FAIL(BS_MB_ERR_TICKS, 9, (double)n_tick);
            }
        }
    }
    pl->node = n_tick + BS_HOV_TICKS;
    pl->n_end = pl->node + BS_MB_HOLD;
    pl->n_clk = pl->n_end + BS_MB_PAD;
    pl->n_path = pl->node + 1;
    pl->hov_in = pl->node - BS_HOV_TICKS;
    if (pl->node > max_ticks) FAIL(BS_MB_ERR_TICKS, 9, (double)pl->node);

    /* seam stores: slot 0 identity/zero, then one slot per distinct
     * (chi, crossing pace), with each interior vertex's slot in
     * vslot[].  THE KICK IS AT THE PUBLISHED PACE OF THE CROSSING,
     * off_vec(v_j, v_j, chi) -- the certified tables always carried
     * the reference's pace at the seam, and a trim-pace kick at a slow
     * published crossing injects a phantom cross-rate the model then
     * has to fight (measured on the SN77 gate: 4.0 of face violation
     * at the first 72-deg seam; with the crossing-pace kick the v_adm
     * cap makes the injected |v_j sin chi| sit inside the certified
     * cold-injection basin by construction). */
    pl->n_seam = 1;
    memset(pl->rot, 0, sizeof(double) * 36);
    for (int d = 0; d < 6; ++d) pl->rot[d * 6 + d] = 1.0;
    memset(pl->off, 0, sizeof(double) * 6);
    memset(pl->vslot, 0, sizeof(pl->vslot));
    for (int i = 1; i + 1 < n_v; ++i) {
        double c6[6];
        seam_off(c6, chi[i], vb[i]);
        int slot = -1;
        for (int k = 1; k < pl->n_seam; ++k) {
            if (fabs(pl->off[k * 6 + 0] - c6[0]) < 1e-12 &&
                fabs(pl->off[k * 6 + 4] - c6[4]) < 1e-12) {
                slot = k;
                break;
            }
        }
        if (slot < 0) {
            if (pl->n_seam >= BS_MB_SEAM_CAP) {
                FAIL(BS_MB_ERR_GEOM, 5, pl->n_seam);
            }
            slot = pl->n_seam++;
            seam_rot(&pl->rot[(size_t)slot * 36], chi[i]);
            memcpy(&pl->off[(size_t)slot * 6], c6, sizeof(c6));
        }
        pl->vslot[i] = slot;
    }

    /* ---- families ---- */
    build_rows_t(&pl->rows[0], v_leg);
    build_rows_h(&pl->rows[(size_t)BS_NROW * 10]);
    if (family_dare(&pl->rows[0], &pl->P[0], &pl->K[0]) != 0) {
        FAIL(BS_MB_ERR_DARE, 10, 0.0);
    }
    if (family_dare(&pl->rows[(size_t)BS_NROW * 10], &pl->P[36],
                    &pl->K[18]) != 0) {
        FAIL(BS_MB_ERR_DARE, 10, 1.0);
    }
    /* family h cross-check against the emitted reference: a free sanity
     * gate on the on-target DARE solver. */
    {
        double dmax = 0.0;
        for (int i = 0; i < 36; ++i) {
            const double d = fabs(pl->P[36 + i] - bs_rc_Ph[i]);
            if (d > dmax) dmax = d;
        }
        for (int i = 0; i < 18; ++i) {
            const double d = fabs(pl->K[18 + i] - bs_rc_Kh[i]);
            if (d > dmax) dmax = d;
        }
#ifndef BS_STAGE2_QD
        if (dmax > 1e-7) FAIL(BS_MB_ERR_DARE, 10, dmax);
#else
        (void)dmax;   /* Stage-2 Q_delta override: the pair cannot reproduce the record's table by construction */
#endif
    }

    /* input waypoint -> vertex tick; merged/dropped inputs inherit the
     * surviving predecessor's tick, and the final waypoint maps to
     * `node` (the hover terminal AT that waypoint). */
    pl->n_wp_in = n_wp;
    {
        int ci = 0;
        for (int i = 0; i < n_wp; ++i) {
            while (ci + 1 < n_v && src[ci + 1] <= i) ci++;
            pl->wp_tick[i] = (ci == n_v - 1) ? pl->node : g_tick[ci];
        }
        pl->wp_tick[n_wp - 1] = pl->node;
    }

    pl->n_v = n_v;
    pl->v_leg = v_leg;
    pl->v_trim = v_trim;
    pl->cell_t = cell_t;
    pl->cell_min = cell_min;
    pl->hyst = hyst;
    pl->r_h = r_h;
    pl->d_hov = d_hov;
    rep->status = BS_MB_OK;
    rep->n_ticks = pl->n_clk;
    rep->n_seam = pl->n_seam;
    rep->n_vert = n_v;
    return BS_MB_OK;
#undef FAIL
}

/* ---------------------------------------------------------------- render */
/* Render one unit of the per-tick tables from the plan's stored carry.
 * Unit i covers leg i: path ticks g[i]+1..g[i+1] (plus tick 0 when
 * i == 0), psi and schedule metadata for ticks g[i]..g[i+1]-1, seam
 * metadata at the destination tick g[i] for interior vertices.  The
 * LAST unit (i == n_v-2) additionally renders the hover splice, the
 * hold and the pad through n_clk-1.  Array indices wrap with `mask`
 * (0x7fffffff = flat arrays, W-1 = a W-slot ring).
 *
 * The integration re-runs the pass-1 recursion (bs_arc_step.h) from
 * the stored carry, so the doubles are BIT-IDENTICAL to the counting
 * pass; positions are scaled so each breakpoint lands exactly on its
 * tick.
 *
 * ALIASING CONTRACT: the unscaled arc is staged in scratch[0..nt-1]
 * and consumed BACKWARD (k = nt-1 .. 0).  Backward consumption makes
 * scratch = &path[2*(g[i]+1)] legal: iteration k writes path cells at
 * scratch-relative offsets 2k and 2k+1, both strictly beyond every
 * still-unread scratch cell (offset k' < k <= 2k), and the staging
 * region [2*(g[i]+1), 2*(g[i]+1)+nt-1] lies wholly inside the leg's
 * own output span.  That is how the batch builder stages inside its
 * path array instead of a leg-sized stack buffer -- the solver thread
 * has 8 kB of stack.  A disjoint scratch (the streaming driver's) is
 * trivially safe. */
int bs_mission_render_unit(const bs_mission_plan *pl, int leg_i,
                           double *path, double *psi, double *ang,
                           signed char *fam, int *sfam, int *srot,
                           int *soff, int mask,
                           double *scratch, int scratch_len)
{
    const int t0 = pl->g_tick[leg_i];
    const int nt = pl->g_tick[leg_i + 1] - t0;
    if (nt > scratch_len) return -1;
    const double L = pl->leg_len[leg_i]
                   - ((leg_i == pl->n_v - 2) ? pl->d_hov : 0.0);
    const double vt = pl->vb[leg_i + 1];
    double v = pl->carry_v[leg_i], a = pl->carry_a[leg_i], arc = 0.0;
    for (int k = 0; k < nt; ++k) {
        bs_arc_step(&v, &a, &arc, L, vt, pl->v_trim);
        scratch[k] = arc;
    }
    const double scale = (arc > 0.0) ? (L / arc) : 0.0;
    const double ch = cos(pl->leg_hd[leg_i]), sh = sin(pl->leg_hd[leg_i]);
    if (leg_i == 0) {
        path[2 * (0 & mask)] = pl->px[0];
        path[2 * (0 & mask) + 1] = pl->py[0];
    }
    for (int k = nt - 1; k >= 0; --k) {      /* backward: see contract */
        const double s_al = scratch[k] * scale;
        const int s = (t0 + k + 1) & mask;
        path[2 * s] = pl->px[leg_i] + s_al * ch;
        path[2 * s + 1] = pl->py[leg_i] + s_al * sh;
    }
    /* chart heading + schedule metadata: this leg owns ticks
     * [g[i], g[i+1]-1]; the seam's destination tick carries the NEW leg
     * heading, so psi at tick g[i] (vertex i) belongs to leg i already. */
    for (int k = 0; k < nt; ++k) {
        const int t = t0 + k;
        const int s = t & mask;
        psi[s] = pl->leg_hd[leg_i];
        fam[s] = (signed char)((t >= pl->hov_in) ? 1 : 0);
        sfam[s] = (t >= pl->hov_in) ? BS_MB_FAM_H : BS_MB_FAM_T;
        if (k == 0 && leg_i >= 1) {          /* seam at vertex i */
            ang[s] = pl->chi[leg_i];
            srot[s] = pl->vslot[leg_i];
            soff[s] = pl->vslot[leg_i];
        } else {
            ang[s] = 0.0;
            srot[s] = 0;
            soff[s] = 0;
        }
    }
    if (leg_i + 2 >= pl->n_v) {
        /* LAST unit: hover splice -- the certified anchor profile,
         * scaled by r_h, along the final heading, ending at rest on the
         * last vertex -- then hold and pad metadata through n_clk-1. */
        const double hch = cos(pl->leg_hd[pl->n_v - 2]);
        const double hsh = sin(pl->leg_hd[pl->n_v - 2]);
        for (int k = 0; k <= BS_HOV_TICKS; ++k) {
            const int t = pl->hov_in + k;
            const int s = t & mask;
            const double s_al = bs_rc_hov_x[k] * pl->r_h;   /* -20r .. 0 */
            path[2 * s] = pl->px[pl->n_v - 1] + s_al * hch;
            path[2 * s + 1] = pl->py[pl->n_v - 1] + s_al * hsh;
            if (t < pl->n_path) psi[s] = pl->leg_hd[pl->n_v - 2];
        }
        for (int t = pl->hov_in; t < pl->n_clk; ++t) {
            const int s = t & mask;
            fam[s] = (signed char)1;         /* t >= hov_in throughout */
            sfam[s] = BS_MB_FAM_H;
            srot[s] = 0;
            soff[s] = 0;
            ang[s] = 0.0;
        }
    }
    return 0;
}

/* ----------------------------------------------------------------- build */
bs_mb_status bs_mission_build(const double *vx, const double *vy, int n_wp,
                              const bs_mission_params *pp,
                              void *mem, size_t mem_len,
                              bs_mission_tables *out,
                              bs_mission_report *rep)
{
    /* ONE plan serves the whole build.  It is ~30 kB -- far beyond the
     * solver thread's 8 kB stack -- and it cannot be carved from `mem`
     * without changing bs_mission_size()'s refusal boundary, so it
     * lives in a file-scope static.  Builds are serialized by
     * construction: bs_mission_build runs only in solver-thread context
     * while the solver is idle (see AP_BSMissionBuilder). */
    static bs_mission_plan plan;
    bs_mission_plan *const pl = &plan;

    const bs_mb_status pst =
        bs_mission_plan_build(vx, vy, n_wp, pp, BS_MB_MAX_TICKS, pl, rep);
    if (pst != BS_MB_OK) return pst;

#define FAIL(st, g, d) do { rep->status = (st); rep->gate = (g); \
                            rep->detail = (d); return (st); } while (0)

    const int n_v = pl->n_v;
    const int node = pl->node;
    const int n_clk = pl->n_clk;
    const int n_path = pl->n_path;

    /* ---- allocate the block ---- */
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

    /* ---- copy the plan's whole-mission tables into the block ---- */
    memcpy(rows, pl->rows, sizeof(pl->rows));
    memcpy(P, pl->P, sizeof(pl->P));
    memcpy(K, pl->K, sizeof(pl->K));
    memcpy(rot, pl->rot, sizeof(double) * 36 * (size_t)pl->n_seam);
    memcpy(off, pl->off, sizeof(double) * 6 * (size_t)pl->n_seam);

    /* ---- render the per-tick tables, unit by unit, flat mask.  The
     * scratch is the path array itself (see the renderer's aliasing
     * contract): no leg-sized buffer ever touches the 8 kB stack. ---- */
    for (int i = 0; i + 1 < n_v; ++i) {
        const int nt = pl->g_tick[i + 1] - pl->g_tick[i];
        if (bs_mission_render_unit(pl, i, path, psi, angt, fam, sfam,
                                   srot, soff, 0x7fffffff,
                                   &path[2 * (pl->g_tick[i] + 1)], nt)
            != 0) {
            FAIL(BS_MB_ERR_GEOM, 2, (double)i);      /* unreachable */
        }
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
        if (st > pl->cell_t * (1.0 + 1e-9) + 1e-12) {
            FAIL(BS_MB_ERR_GATE, 2, st);
        }
    }
    /* G4: every cleaned vertex lands exactly on its tick */
    for (int i = 0; i < n_v; ++i) {
        const int gt = (i == n_v - 1) ? node : pl->g_tick[i];
        const double dv = hypot(path[2 * gt] - pl->px[i],
                                path[2 * gt + 1] - pl->py[i]);
        if (dv > 1e-6) FAIL(BS_MB_ERR_GATE, 4, dv);
    }
    /* G6: published pace below the face */
    for (int t = 1; t < n_path; ++t) {
        const double st = hypot(path[2 * t] - path[2 * (t - 1)],
                                path[2 * t + 1] - path[2 * (t - 1) + 1]);
        if (st / BS_TS > pl->v_leg - 0.05) {
            FAIL(BS_MB_ERR_GATE, 6, st / BS_TS);
        }
    }

    /* ---- publish ---- */
    double total = 0.0;
    for (int i = 0; i + 1 < n_v; ++i) total += pl->leg_len[i];
    out->n_path = n_path;
    out->n_clk = n_clk;
    out->node = node;
    out->n_end = pl->n_end;
    out->hov_in = pl->hov_in;
    out->n_seam = pl->n_seam;
    out->v_leg = pl->v_leg;
    out->v_trim = pl->v_trim;
    out->cell_t = pl->cell_t;
    out->cell_min = pl->cell_min;
    out->hyst = pl->hyst;
    out->v0 = pl->vb[0];
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
    for (int i = 0; i < n_wp; ++i) out->wp_tick[i] = pl->wp_tick[i];

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
    rep->n_seam = pl->n_seam;
    rep->n_vert = n_v;
    rep->build_bytes = (double)need;
    return BS_MB_OK;
#undef FAIL
}
