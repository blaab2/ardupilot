/* bs_mission_build — the RUNTIME mission builder for the bsolver anytime
 * MPC.  Turns an uploaded waypoint polyline into the solver's working
 * tables, on the aircraft, at mission-upload time.  Pure C, no libc
 * allocation (caller provides one block), host-testable.
 *
 * ARCHITECTURE ("a junction everywhere", owner decision 2026-09-02):
 * there is no turn-window construction.  Every mission vertex is a bare
 * SEAM — chart rotation T_switch(chi) plus the affine kick
 * off_vec(v_j, v_j, chi) at the CROSSING pace, at the destination tick —
 * with plain row-family "t" cells on both sides.  All pace variation
 * lives in the PUBLISHED reference arc: a jerk-shaped profile that
 * decelerates to v_j(chi) into each vertex and re-accelerates after (the
 * flight-verified BSLV_INGRESS ramp, generalized), starts leg 1 from
 * rest, and
 * ends with the certified hover-anchor deceleration scaled to the
 * final-leg entry pace.  The model chart is pace-agnostic — uniform
 * cells, families t then h, exactly the structure of the flown
 * SN77 tables (whose hover transition carries no offset either).
 *
 * The speed system derives from the configured cruise speed at build
 * time: VLEG = v_cap*cos(pi/20) (inscribed 20-gon), V_TRIM =
 * VLEG - 0.2888, CELL_T = Ts*V_TRIM, HYST = 0.74116*CELL_T.  Families
 * t (speed margins from VLEG) and h (constant floors) are rebuilt per
 * mission; their terminal pairs come from bs_dare_solve, with family
 * h cross-checked against the emitted reference as a solver sanity
 * gate.
 *
 * THE CORNER PACE IS DERIVED AND THEN MEASURED.  The corridor the
 * barrier already allows IS a cornering radius: an arc tangent to both
 * legs departs THE LEGS (which is what the corridor measures) by
 * d = r (1 - cos(chi/2)), so the geometric ceiling is
 *     v_j(chi) = sqrt(kappa a_lat d / (1 - cos(chi/2))).
 * That ceiling prices the rounding arc but not the along/cross coupling
 * the closed loop actually pays, so it is capped by a table measured IN
 * THE CLOSED LOOP (tests/calibrate_seam_pace.py): per junction angle,
 * the largest crossing pace with zero face violations, zero re-timings
 * and a model corridor peak inside the tracking budget.
 * The frame pace has to be authored at all because the optimizer cannot
 * slow its own reference by more than the lag band; its VALUE is the
 * speed at which the optimizer's own constraint set is feasible, and the
 * rounding within the band is then the optimizer's to find.
 */
#ifndef BS_MISSION_BUILD_H
#define BS_MISSION_BUILD_H

#include <stddef.h>

#include "bs_solver.h"

#define BS_MB_MAX_WP 64            /* input waypoints (run) cap          */
#define BS_MB_MAX_TICKS 1800       /* mission clock cap                  */
#define BS_MB_HOLD 12              /* terminal hold ticks (of record)    */
#define BS_MB_PAD (BS_N + 8)       /* clock padding beyond the hold      */

typedef struct {
    double v_cap_ms;               /* configured cruise (WPNAV_SPEED)    */
} bs_mission_params;

typedef enum {
    BS_MB_OK = 0,
    BS_MB_ERR_NWP,                 /* < 2 or > cap waypoints             */
    BS_MB_ERR_SPEED,               /* v_cap outside [3.0, 11.8]          */
    BS_MB_ERR_GEOM,                /* degenerate after cleanup           */
    BS_MB_ERR_TICKS,               /* > BS_MB_MAX_TICKS                  */
    BS_MB_ERR_MEM,                 /* block too small                    */
    BS_MB_ERR_DARE,                /* terminal pair failed               */
    BS_MB_ERR_HOVER,               /* final-leg entry pace < 1 m/s       */
    BS_MB_ERR_GATE                 /* an acceptance gate failed          */
} bs_mb_status;

typedef struct {
    bs_mb_status status;
    int gate;                      /* failing gate id (0 = none)         */
    double detail;                 /* gate-specific value                */
    int n_ticks, n_seam, n_vert;
    double build_bytes;
} bs_mission_report;

typedef struct {
    /* clock geometry */
    int n_path;                    /* path samples = node + 1            */
    int n_clk;                     /* schedule length                    */
    int node;                      /* first hold tick (rest at last WP)  */
    int n_end;                     /* node + BS_MB_HOLD (mission end)    */
    int hov_in;                    /* first hover-family tick (node-14)  */
    int n_seam;                    /* seam slots incl. slot 0 = identity */
    /* speed system */
    double v_leg, v_trim, cell_t, cell_min, hyst;
    double v0;                     /* leg-1 published start pace (rest)  */
    double length_m;
    /* per-tick tables */
    double *path;                  /* [2*n_path] published NE            */
    double *psi;                   /* [n_path] chart heading, rad        */
    double *ang;                   /* [n_clk] seam angle at DEST tick, deg */
    signed char *fam;              /* [n_clk] 0 = t, 1 = h               */
    int *sfam, *srot, *soff;       /* [n_clk] core schedule indices      */
#if BS_STAGE2_PACE_TAB
    double *pace;                  /* [n_clk] published pace, m/s (2c)   */
#endif
    /* family tables (t = 0, h = 1) */
    double *rows;                  /* [2*BS_NROW*10]                     */
    double *P;                     /* [2*36]                             */
    double *K;                     /* [2*18]                             */
    /* seam stores; slot 0 = identity / zero                             */
    double *rot;                   /* [n_seam*36]                        */
    double *off;                   /* [n_seam*6]                         */
    /* tau of each INPUT waypoint's vertex (merged ones inherit the
     * survivor's tick); [BS_MB_MAX_WP] */
    int wp_tick[BS_MB_MAX_WP];
    int n_wp_in;
    /* the schedule the core consumes (points at the arrays above)       */
    bs_schedule sched;
} bs_mission_tables;

/* Upper bound on the block size for a mission of n_wp waypoints whose
 * polyline length is len_m at cruise v_cap.  Call before allocating. */
size_t bs_mission_size(int n_wp, double len_m, double v_cap_ms);

/* host calibration hook: >0 caps every junction pace (target: unused) */
extern double bs_mb_vj_cap;
extern double bs_mb_vj_scale;   /* host experiments: scale the junction pace (1.0 = shipped) */

/* ---- streaming (two-level) build ----------------------------------
 * PLAN AT VERTEX LEVEL GLOBALLY, RENDER TO TICKS LOCALLY.  The plan is
 * everything whole-mission (cleanup, junction paces, reachability,
 * hover splice, per-vertex tick indices and integrator carry, deduped
 * seam stores, families) and is small regardless of mission length.
 * The per-tick tables render leg by leg from the stored carry -- into
 * flat arrays (the batch builder, mask = 0x7fffffff) or into a W-slot
 * ring (the streaming driver, mask = W-1).  One renderer serves both:
 * the stream==batch byte-equality gate (tests/stream_gate.c) holds by
 * construction and is verified anyway. */
#define BS_MB_SEAM_CAP 48              /* deduped (chi, pace) seam slots */

typedef struct {
    int n_v;
    double px[BS_MB_MAX_WP], py[BS_MB_MAX_WP];
    double leg_len[BS_MB_MAX_WP], leg_hd[BS_MB_MAX_WP];
    double chi[BS_MB_MAX_WP], vb[BS_MB_MAX_WP + 1];
    double carry_v[BS_MB_MAX_WP + 1], carry_a[BS_MB_MAX_WP + 1];
    int g_tick[BS_MB_MAX_WP + 1];
    int vslot[BS_MB_MAX_WP + 1];
    int node, n_end, n_clk, n_path, hov_in, n_seam;
    double v_leg, v_trim, cell_t, cell_min, hyst, r_h, d_hov;
    double rows[2 * BS_NROW * 10], P[2 * 36], K[2 * 18];
    double rot[BS_MB_SEAM_CAP * 36], off[BS_MB_SEAM_CAP * 6];
    int wp_tick[BS_MB_MAX_WP], n_wp_in;
#if BS_STAGE2_PACE_TAB
    int pace_t0;                   /* first tick of leg 0 at cruise pace */
#endif
#if BS_STAGE2_PAIR || BS_STAGE2_REPACE
    /* Stage 2 (c): per-leg pace ceiling for the integrator (v_trim on
     * every leg unless a corner pair holds its pace across the short
     * leg between the two corners) */
    double leg_cap[BS_MB_MAX_WP];
#endif
#if BS_STAGE2_REPACE
    /* the authored (base) paces and ceilings the online re-pace scales,
     * the current scale, and the input-waypoint -> cleaned-vertex map
     * the tick map is refreshed from after a re-pace */
    double vb_base[BS_MB_MAX_WP + 1], cap_base[BS_MB_MAX_WP];
    double pace_scale;
    int wp_ci[BS_MB_MAX_WP];
#endif
} bs_mission_plan;

/* Stage 2 (c) re-pace tuning + the verdict helpers (always visible so the
 * host gate's accumulators compile on every build; the API below is
 * gated) */
#ifndef BS_STAGE2_REPACE_MAX
#define BS_STAGE2_REPACE_MAX 1.30   /* pace scale ceiling                */
#endif
#ifndef BS_STAGE2_REPACE_GAIN
#define BS_STAGE2_REPACE_GAIN 0.15  /* scale step per unit verdict         */
#endif
#ifndef BS_STAGE2_REPACE_DB
#define BS_STAGE2_REPACE_DB 0.10    /* verdict dead band: |u| below it holds */
#endif
#ifndef BS_STAGE2_REPACE_RIDE
#define BS_STAGE2_REPACE_RIDE 0.8   /* band-riding: e_l above this x band   */
#endif
#ifndef BS_STAGE2_REPACE_MIN_N
#define BS_STAGE2_REPACE_MIN_N 4    /* slow ticks needed for a verdict   */
#endif
#ifndef BS_STAGE2_REPACE_NU_REF
#define BS_STAGE2_REPACE_NU_REF 0.5  /* mean |nu| that cancels a full-band lead */
#endif
#ifndef BS_STAGE2_REPACE_VIOL_REF
#define BS_STAGE2_REPACE_VIOL_REF 1.0 /* face violation that cancels a full-band lead */
#endif
#ifndef BS_STAGE2_REPACE_VIOL_DZ
#define BS_STAGE2_REPACE_VIOL_DZ 1.0  /* dead zone: the relaxed barrier's own
                                       * small excursions (0.1-0.55 under the
                                       * Stage-2 (a) cost) are not a pace signal */
#endif
/* THE RE-PACE VERDICT on the slow ticks since the previous leg render:
 *     u = ride - mean|nu|/nu_ref - max(0, viol_max - dz)/viol_ref
 * `ride` is the fraction of those ticks on which the plan sits at the
 * lag band (e_l > 0.8 band): a frame the model can follow is led by
 * 0.8-0.9 m on average but the band is touched only briefly (SN77 pair
 * pace: ride 0.05), while a frame that is too slow has the plan parked
 * at the band (the record's finding: every frontier rides the 3 m
 * band).  The MEAN lead is not the signal -- it is +0.3 band at every
 * pace under the Stage-2 (a) cost and drove a 1.0 <-> 1.07 flip-flop
 * against the health terms.  The other two terms are the model's HEALTH: a
 * frame the model cannot follow does not show up as lag -- the clock
 * input nu (cheap, R = 0.01) absorbs it and e_l stays positive while
 * the published track exceeds the disc (measured at scale 1.3 on SN77:
 * mean|nu| 0.08 -> 0.9, published v_max 13.2 > VLEG, jerk faces
 * violated).  So mean |nu| and the worst face violation of the applied
 * input (seam ticks excluded: the frame rotation re-expresses the
 * acceleration polygon there) vote AGAINST the pace.  The update is
 * multiplicative and bounded to [1, BS_STAGE2_REPACE_MAX]. */
static inline double bs_mission_repace_verdict(double ride, double nu_mean,
                                               double viol_max)
{
    double vx = viol_max - BS_STAGE2_REPACE_VIOL_DZ;
    if (vx < 0.0) vx = 0.0;
    return ride
         - nu_mean / BS_STAGE2_REPACE_NU_REF
         - vx / BS_STAGE2_REPACE_VIOL_REF;
}
static inline double bs_mission_repace_scale(double s, double u)
{
    if (u > 1.0) u = 1.0;
    if (u < -1.0) u = -1.0;
    if (u > -BS_STAGE2_REPACE_DB && u < BS_STAGE2_REPACE_DB) return s;
    s *= 1.0 + BS_STAGE2_REPACE_GAIN * u;
    if (s < 1.0) s = 1.0;
    if (s > BS_STAGE2_REPACE_MAX) s = BS_STAGE2_REPACE_MAX;
    return s;
}

#if BS_STAGE2_REPACE
/* Re-pace every not-yet-rendered leg (>= leg_i) and vertex (>= leg_i+1)
 * to `scale` x the authored pace, then re-integrate the clock from the
 * stored carry at vertex leg_i.  Call BEFORE rendering leg leg_i; every
 * tick < g_tick[leg_i] and the seam at vertex leg_i are untouched.  On
 * refusal (seam-slot table full, clock cap) the plan is restored and a
 * negative code returned. */
int bs_mission_plan_repace(bs_mission_plan *pl, int leg_i, double scale,
                           int max_ticks);
/* the batch builder's plan (host gates emulate the streaming driver
 * on the batch tables) */
bs_mission_plan *bs_mission_batch_plan(void);
#endif

#if BS_STAGE2_WIN
/* the margin of row i of the pace family at a tick, as the core reads it
 * (bs_solver.c row_margin, chart frame): shared by the driver's and the
 * host gate's violation ledgers */
static inline double bs_mb_win_margin(const bs_schedule *sch, int fam_tau,
                                      int tau, double pace_tau, int i,
                                      const double *r)
{
    if (tau >= sch->pace_t0 && fam_tau == sch->pace_fam) {
        if (i < 40 && (i & 1) == 0) return sch->pace_vleg - pace_tau * r[0];
        if ((i == 40 || i == 44) && pace_tau < sch->pace_vslow)
            return BS_STAGE2_WIN_BAND;
    }
    return r[9];
}
#endif

#if BS_STAGE2_PACE_TAB
/* Stage 2 (c): the published pace per tick of unit leg_i, derived from
 * the rendered path (|p(t+1) - p(t)| / Ts; 0 at rest), same mask/ring as
 * the path.  Call after bs_mission_render_unit(leg_i). */
void bs_mission_render_pace(const bs_mission_plan *pl, int leg_i,
                            const double *path, double *pace, int mask);
#endif

/* Build the vertex-level plan.  max_ticks: pass BS_MB_MAX_TICKS for the
 * batch semantics (refusal at the classic cap); larger for streaming. */
bs_mb_status bs_mission_plan_build(const double *vx, const double *vy,
                                   int n_wp, const bs_mission_params *pp,
                                   int max_ticks, bs_mission_plan *pl,
                                   bs_mission_report *rep);

/* Render one unit of the per-tick tables: unit i = leg i for
 * i < n_v-2; unit n_v-2 = the last leg PLUS the hover splice, hold and
 * pad (finalizes through n_clk-1).  Array indices are wrapped with
 * `mask` (0x7fffffff = flat).  scratch stages one leg's unscaled arc;
 * returns 0, or -1 if the leg exceeds scratch_len ticks. */
int bs_mission_render_unit(const bs_mission_plan *pl, int leg_i,
                           double *path, double *psi, double *ang,
                           signed char *fam, int *sfam, int *srot,
                           int *soff, int mask,
                           double *scratch, int scratch_len);

/* cold-injection admissible junction speed (the ingress budget) */
double bs_mb_v_adm(double chi_deg);

/* Build.  `mem` is one caller-owned block of `mem_len` bytes; `out`
 * points into it after success.  On failure nothing is published and
 * `rep` carries the reason. */
bs_mb_status bs_mission_build(const double *vx, const double *vy, int n_wp,
                              const bs_mission_params *pp,
                              void *mem, size_t mem_len,
                              bs_mission_tables *out,
                              bs_mission_report *rep);

#endif /* BS_MISSION_BUILD_H */
