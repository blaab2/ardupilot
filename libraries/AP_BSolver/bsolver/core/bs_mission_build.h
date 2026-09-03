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
