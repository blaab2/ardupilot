/* bs_arc — the ARC-CELL CLOCK runtime (fillet chart) on top of the bsolver
 * core.  The C port of model/stage2_arc_clock.py + stage2_arc_clock_mission.py
 * (--fillet), the controller papers/p4-anytime-racing-mpc/tex_journal is
 * written around, so that it can be flown (SITL, Sim-on-Hardware) and timed
 * on the flight processor.  Reference archive: data/derived/p4-current/
 * (tooling/current_evidence.py), which this layer must reproduce.
 *
 * WHAT THE CLOCK IS.  The mission polyline is cut into CELLS: straight cells
 * on the legs (corner pace 8 m/s within 20 m before / 16 m after a vertex,
 * trim pace elsewhere) and CHORD cells along a fillet at every vertex (the
 * largest arc inside the capsule union, r_f = w / (1 - cos(chi/2)), chords
 * <= 18 deg and <= one corner-pace cell).  The state at a stage assigned to
 * cell k is the phase quotient relative to that cell,
 *     xi = (delta = v - pace_k, a, e_l = s - s_start(k), e_c, edot_c, a_n),
 * and the cell advances only when the plan has EARNED it: after one plant
 * step the raw e_l is the position relative to the NEXT cell's start; a
 * crossing (by position past the interior bisector, or by heading past
 * the bisector inside the vertex box) advances the cell and re-expresses
 * the state through the seam map T_switch(chi), off_vec(pace_k, pace_k+1,
 * chi); otherwise the frame point stays and the transition carries +l_k on
 * e_l.  Nothing advances because time has elapsed; there is no ledger and
 * no clock input (v1: nu is pinned to |nu| <= 0.01; v2: nu does not exist).
 *
 * MODEL VERSIONS.  The layer is dimension-generic and follows the generated
 * header selected at compile time (BS_DATA_HEADER):
 *   v1  bs_model_runtime_const.h     u = (j_l, nu, j_c), BS_NU 3, BS_NV 90,
 *                                    BS_NROW 50, forward-Euler e_l
 *                                    (the P4 archive data/derived/p4-current)
 *   v2  bs_model_runtime_const_v2.h  u = (j_l, j_c), BS_NU 2, BS_NV 60,
 *                                    BS_NROW 48, exact zero-order-hold e_l
 *                                    (e_l+ = e_l + Ts delta + Ts^2/2 a_l +
 *                                    Ts^3/6 j_l), BS_MODEL_V2 defined
 * Everything else -- chart, crossing rule, seam maps, margins, trust
 * region, reward, engage, rounds, shift-append -- is identical; the v1
 * build is bit-identical to the pre-v2 layer.
 *
 * ENGAGE SEED (v2b).  The engage's three fixed-point rounds start from a
 * seed U.  With U = 0 the first rollout assigns every stage to cell 0 and
 * the terminal e_l weight makes the plan back up 2-3 m from rest before
 * accelerating (the archives' pre-arrival corridor maxima are that
 * back-up; doc/ARC-FLIGHTS-2026-09-11.md section 8).  The driver's
 * `seed_ramp` (1 by default = v2b) seeds the jerk-limited acceleration
 * ramp to the first cell's pace instead (bs_arc_ramp_seed); 0 keeps
 * U = 0, which reproduces the v1 (p4-current) and v2 (p4-v2) archives.
 *
 * WHAT A SOLVE IS.  Before each solve the shifted warm start is rolled out
 * with the same rule to produce the stage -> cell assignment, the composed
 * per-stage transitions (T_t, c_t) and the rollout's within-cell progress
 * e_l(t) (the anchor of the 12 m schedule trust region).  Frozen, those
 * data make the tick's problem exactly the scheduled-LTV relaxed-barrier
 * problem the core solves: the arc layer RENDERS them into a bs_schedule
 * (family, rotation/offset slots, per-stage margins, the tip family's
 * terminal pair) and calls bs_problem_init_pinned + bs_newton_pinned.  Every
 * family shares one (h, e) row template -- families differ only in their
 * margins (speed disc about the cell's own pace, iso/aniso acceleration
 * polygon, corridor rows, the fillet's asymmetric corridor) -- so the
 * schedule carries a per-stage MARGIN table instead of per-stage row blocks.
 *
 * The core is compiled with -DBS_ARC=1 for this layer: that switches on the
 * per-stage margin table in bs_schedule, the (a) re-priced stage cost
 * (Q_delta = Q_el = 0; the DARE keeps RL, exactly as the Python does) and
 * the 1.5 progress reward.  The default build is byte-identical.
 */
#ifndef BS_ARC_H
#define BS_ARC_H

#include <stddef.h>

#include "bs_solver.h"
#include "bs_dare.h"

#ifndef BS_ARC_MAX_WP
#define BS_ARC_MAX_WP 64
#endif

/* controller constants of record (tooling/current_evidence.py) */
#define BS_ARC_CORNER_PACE 8.0       /* chart pace of the corner cells, m/s */
#define BS_ARC_W 2.0                 /* corridor half-width, m              */
#define BS_ARC_BAND 1000.0           /* lag band: none (trust region instead) */
#define BS_ARC_TRUST 12.0            /* schedule trust region, m            */
#define BS_ARC_ROUNDS 3              /* assignment/solve rounds per tick    */
#define BS_ARC_TAIL 60.0             /* trailing straight after the last WP */
#define BS_ARC_MAX_DPSI 18.0         /* chord angle cap, deg                */
#define BS_ARC_VERTEX_BOX 2.0        /* capsule disc at a vertex, m         */
#define BS_ARC_ACCC 5.66             /* iso acceleration tier, m/s^2        */
#define BS_ARC_M_NU 0.01             /* v1 only: clock input pin, both sides */
#define BS_ARC_SEED_MAXIT 300        /* engage: unpinned Newton cap         */
#define BS_ARC_SEED_TOL 1e-7         /* engage: gradient tolerance          */
#define BS_ARC_SEED_ROUNDS 3         /* engage: assignment fixed-point rounds */
#define BS_ARC_WIN_IN 20.0           /* corner-pace straight before a vertex */
#define BS_ARC_WIN_OUT 16.0          /* corner-pace straight after a vertex  */

enum { BS_ARC_TA = 0, BS_ARC_CA = 1, BS_ARC_FILLET = 2 };

typedef struct {
    double px, py;        /* start point (N, E)                              */
    double thx, thy;      /* unit heading; n_hat = (-thy, thx)               */
    double l;             /* cell length; chart pace = l / BS_TS             */
    double chi;           /* seam angle at the cell START, deg (0 = none)    */
    double d_in;          /* fillet cell: chord-midpoint depth into the corner */
    signed char kind;     /* BS_ARC_TA / BS_ARC_CA / BS_ARC_FILLET           */
    signed char sense;    /* fillet: +1 / -1 (turn sense), else 0           */
    short leg;            /* polyline leg index                              */
} bs_arc_cell;

typedef struct {
    bs_arc_cell *cells;
    int n_cells;
    int k_end;            /* arrival cell = first trailing cell              */
    int n_v;              /* vertices kept                                   */
    double vx[BS_ARC_MAX_WP], vy[BS_ARC_MAX_WP];
    double chi_v[BS_ARC_MAX_WP];       /* turn angle at each vertex, deg     */
    int leg_first[BS_ARC_MAX_WP];      /* first cell of leg i's straight     */
    int leg_last[BS_ARC_MAX_WP];       /* last cell (incl. fillet) of leg i  */
    double v_trim, v_leg, corner_pace, w;
} bs_arc_chart;

/* Upper bound on the cell count / block size for n_wp waypoints of total
 * polyline length len_m (call before allocating). */
int bs_arc_chart_cells_bound(int n_wp, double len_m);
size_t bs_arc_chart_size(int n_wp, double len_m);

/* Build the fillet chart.  `mem` holds the cells (>= bs_arc_chart_size).
 * Returns 0, or a negative code (-1 waypoints, -2 memory, -3 geometry). */
int bs_arc_chart_build(const double *vx, const double *vy, int n_wp,
                       void *mem, size_t mem_len, bs_arc_chart *out);

/* world position / velocity / acceleration of a quotient state in cell k */
void bs_arc_world(const bs_arc_chart *ch, int k, const double *xi,
                  double *pn, double *pe, double *vn, double *ve,
                  double *an, double *ae);

/* ---- the per-solve schedule ---------------------------------------- */
#define BS_ARC_NSLOT (BS_N + 1)
#define BS_ARC_DARE_CACHE 6

typedef struct {
    int family[BS_ARC_NSLOT], rotation[BS_ARC_NSLOT], offset[BS_ARC_NSLOT];
    double rows[BS_NROW * 10];             /* the (h, e) template + margins */
    double mstk[BS_N * BS_NROW];           /* per-stage margins (trust in)  */
    double rot[BS_ARC_NSLOT * BS_NX * BS_NX];
    double off[BS_ARC_NSLOT * BS_NX];
    double P[BS_NX * BS_NX], K[BS_NU * BS_NX];  /* the tip family's pair    */
    int n_rot, n_off;
    bs_schedule sched;
    /* the assignment the schedule was rendered from */
    int ks[BS_N + 1];
    double els[BS_N + 1];
    double xiN[BS_NX];                     /* rollout terminal state        */
    /* terminal-pair cache: key = cell family key (see bs_arc_fam_key) */
    struct { int key; double P[BS_NX * BS_NX], K[BS_NU * BS_NX]; }
        dare[BS_ARC_DARE_CACHE];
    int dare_n, dare_next;
    int dare_fail;
} bs_arc_sched;

/* the family key of a cell: 0 = ta, 1 = ca, 2 + k for a fillet cell k */
int bs_arc_fam_key(const bs_arc_chart *ch, int k);

/* margins of cell k's family (no trust), BS_NROW entries in template order */
void bs_arc_margins(const bs_arc_chart *ch, int k, double *m);

/* One plant step + the clock rule: y = A xi + B u, then (k_new, T, c) with
 * xi_new = T y + c.  T (6x6 row-major) and c may be NULL when only the
 * state is wanted.  Returns k_new. */
int bs_arc_transition(const bs_arc_chart *ch, int k, const double *y,
                      double *T, double *c);
void bs_arc_plant(const double *xi, const double *u, double *y);

/* Roll U out from (k0, xi0): the assignment ks[0..N], els[0..N], the
 * composed per-stage (T, c) rendered into S, the terminal state S->xiN. */
void bs_arc_assign(const bs_arc_chart *ch, bs_arc_sched *S, int k0,
                   const double *xi0, const double *U);

/* Render the schedule (margins with trust, terminal pair) after
 * bs_arc_assign.  Returns 0, or -1 if the tip family's DARE failed. */
int bs_arc_render(const bs_arc_chart *ch, bs_arc_sched *S, double trust);

/* ---- the driver ------------------------------------------------------ */
typedef struct {
    const bs_arc_chart *chart;
    bs_arc_sched *S;
    double *arena;
    size_t arena_len;
    bs_problem prob;
    int q, rounds;
    double trust;
    int seed_ramp;            /* engage seed: 1 = jerk-limited ramp (v2b,
                               * the default), 0 = U = 0 (v1 / v2 archives) */
    /* loop state */
    int k;                    /* the vehicle's cell                        */
    int kw;                   /* applied ticks so far                      */
    double xi[BS_NX];
    double U[BS_NV];
    double committed[BS_NU * BS_M_COMMIT];   /* inputs for ticks kw..kw+M-1 */
    int engaged, finished;
    /* per-tick diagnostics */
    double lambda0;
    int bt, armijo_fail, no_step, factor_fail;
    int reassign;             /* rounds whose assignment changed           */
    int consistent;           /* final plan's re-assignment == its schedule */
    int holds, skips;
    double viol;              /* applied row excess of the last tick       */
    int seed_iters;
    /* optional platform hooks: yield between seed iterations (the engage
     * runs for seconds), microsecond clock for the per-tick split */
    void (*yield)(void);
    unsigned long (*micros)(void);
    /* per-tick split of the last tick, us (0 without a clock hook):
     * assignment + render (rollouts, margins, DARE), problem init, Newton */
    unsigned long us_assign, us_init, us_newton;
} bs_arc_driver;

/* Bind the driver to a chart, a schedule block and a solver arena
 * (>= bs_workspace_size() doubles). */
void bs_arc_driver_init(bs_arc_driver *d, const bs_arc_chart *ch,
                        bs_arc_sched *S, double *arena, size_t arena_len,
                        int q, int rounds, double trust);

/* The v2b engage seed (doc/ARC-FLIGHTS-2026-09-11.md section 8): the
 * jerk-limited acceleration ramp from rest to the first cell's chart pace
 * v_target = cells[0].l / BS_TS -- j_l = +BS_JBOX until a_l reaches the
 * forward face AF of cell 0's family (bs_arc_margins at raster index 0:
 * bs_rc_aniso[0] = 3.119 for a trim cell, BS_ARC_ACCC = 5.66 for a
 * corner / fillet cell), hold, j_l = -BS_JBOX so that a_l returns to 0 as
 * the speed arrives at v_target; every other input 0.  Writes all BS_NV
 * entries of U. */
void bs_arc_ramp_seed(const bs_arc_chart *ch, double *U);

/* The engage: at rest on the first vertex, three unpinned fixed-point
 * rounds (assignment + full Newton, <= 300 steps, tol 1e-7) from the
 * seed U (d->seed_ramp: the ramp; else U = 0), commits the first M
 * inputs.  Returns 0 on success. */
int bs_arc_engage(bs_arc_driver *d);

/* One maintained tick, in the order of ac.run: publish (the caller lifts
 * before calling), plant step with the clock rule, shift-append, `rounds`
 * x (assign, render, pinned solve), commit stage M-1.  u0 receives the
 * applied input.  Returns 0 (ticked), 1 (arrived: k >= k_end, nothing
 * applied), or a negative code. */
int bs_arc_tick(bs_arc_driver *d, double *u0);

/* Positions of the committed prefix's nodes (M + 1) from the current
 * state, in the plane (N, E); also the cell of each node. */
void bs_arc_lift(const bs_arc_driver *d, double *pn, double *pe, int *kn);

#endif /* BS_ARC_H */
