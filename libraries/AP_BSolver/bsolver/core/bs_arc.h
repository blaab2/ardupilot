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
 * ROW FAMILIES (v2c: cell-centred rows).  The chart's `cell_rows` (1 by
 * default: bs_arc_chart_build; bs_arc_chart_build_rows takes it as an
 * argument) makes EVERY cell its own row family: the twenty velocity
 * rows are written about the cell's own chart pace l_k / Ts,
 * m_i = VLEG - pace_k cos(theta_i), so that the velocity polygon is the
 * world disc |v| <= VLEG on every cell; the acceleration faces stay by
 * cell kind (the anisotropic vector on trim straights, BS_ARC_ACCC on
 * corner-pace straights and chords), the corridor / band / jerk rows as
 * before, and the terminal DARE pair is per cell (key 2 + k for every
 * cell; the cache is a plain ring).  So that no cell's pace exceeds its
 * nominal pace (V_TRIM on a trim straight, the corner pace on a
 * corner-pace straight; a chord is <= one corner-pace cell by its
 * rule), the straight cell counts are rounded UP under cell rows,
 * n = ceil(L / (v_nom Ts)), so l = L / n <= v_nom Ts and the forward
 * velocity margin VLEG - pace_k is >= 0.289 (trim) / >= 3.65 (corner)
 * on every cell -- the family rows' round() (Python round-half-even)
 * put 4 trim cells of the Turns mission at 12.34 m/s and 14 of the
 * Mixed mission at 11.75, above VLEG.  The gate checks every margin of
 * every cell is positive under cell rows.
 * `cell_rows = 0` restores the FAMILY-centred rows AND the round()
 * cell counts of the p4-current / p4-v2 / p4-v2b archives: 'ta' about
 * V_TRIM, 'ca' about the corner pace, one family per fillet chord
 * about its own pace, with the DARE cache keeping the two straight
 * families resident (the gate option --family-rows; tests/arc_regress.sh
 * runs the archives with it).
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
 *
 * THE STOP AT THE LAST WAYPOINT (v2d, doc/V2D-HOVER-CELL-PLAN.md D1-D9;
 * brake onset rule G1 of doc/V2D-G1-DESIGN-2026-09-25.md, which replaced
 * the stop_profile enumeration of rule G0, doc/V2D-P1-W1-2026-09-22.md,
 * on 2026-09-25; the Python of record is
 * model/arc_v2.py: has_hover .. install_stop and _run).  A chart built
 * with `hover` = 1 (bs_arc_chart_build: BS_ARC_HOVER, the configuration
 * of record since v2d; bs_arc_chart_build_rows takes the flag) has NO
 * trailing straight: the chart ends in ONE cell of kind BS_ARC_HV at the
 * last waypoint (l = 0, pace 0, chi 0, the last leg's tangent), k_end =
 * its index.  Its row family (bs_arc_margins): the velocity disc of
 * radius BS_HV_SPEED_ARC about pace 0, the isotropic BS_ARC_ACCC faces,
 * +-e_l at BS_HV_BAND (the hover zone), +-e_c at BS_HV_EC_BAND, the jerk
 * box -- 48 rows in the template order; its terminal pair is the EMITTED
 * bs_rc_arc_Ph / bs_rc_arc_Kh of the generated header (never solved on
 * target, counted in bs_arc_sched.dare_hover).  The progress reward is
 * per stage: BS_HV_W (0) on a stage whose next state is assigned to the
 * hover cell, BS_SPEED_TILT_W elsewhere (bs_arc_sched.tilt_w ->
 * bs_schedule.tilt_w; a horizon without a hover stage keeps today's
 * arithmetic bit for bit).
 *   GATE   closed from engage until the first tick at which rule G1
 *          yields a brake seed (or the vehicle is already in the hover
 *          cell; doc/V2D-G1-DESIGN-2026-09-25.md, variant cdB).  While
 *          closed the hover cell is NOT assignable: bs_arc_assign rolls
 *          out on the chart without it and the last straight cell holds
 *          beyond its end with T = I, c = 0 once the rollout state has
 *          e_l >= 0 (the 'virtual tail': exactly the v2c tail behaviour).
 *   SEED   every tick while the gate is closed (after the shift-append,
 *          before the rounds): the natural rollout of the shifted warm
 *          start and ONE scan of its stages t = N-1 down to M (closed
 *          form, O(N), no enumeration): v_t = pace + delta, d_t = the
 *          distance to WP_N; stages with v_t <= 0.3 or d_t <= 0 skipped;
 *          a_req = v_t^2 / (2 d_t), rest = t + max(2, ceil(v_t / (a_req
 *          Ts))).  A stage qualifies when a_req <= a_nom = v_c^2 / (2
 *          d_stop) (the stop of record, bs_arc_stop_model) and rest <= N
 *          - BS_HV_GATE_MARGIN; t_b = the latest one.  The seed from
 *          (v_b, a_b, d_b) at t_b is the constant deceleration a_s entered
 *          and left at the jerk box (a trapezoid in acceleration, a_s the
 *          smaller root of one quadratic, re-solved on the Ts grid so the
 *          seed ends at rest); it overwrites the warm start's j_l from t_b
 *          (stages < N) and the latch is set to min(N, its rest stage).
 *          No admissible root at t_b: the gate stays closed this tick
 *          (bs_arc_driver.no_root), no fallback seed.
 *   LATCH  a stage index (-1 none): every stage >= latch whose natural
 *          assignment is short of the hover cell is FORCED into it
 *          (force_to_end: chain through the remaining cells, e_l = the
 *          signed shortfall to WP_N in the hover frame, so the +-e_l band
 *          rows demand the plan to be there by that stage); the latch
 *          moves earlier on natural entry (the first stage > 0 whose
 *          rollout cell is >= k_end), never later.
 *   ARRIVAL = zone entry: the first post-step state with e_l >= -BS_HV_BAND
 *          in the hover frame (the state itself in the hover cell, e_l -
 *          l_{k_end-1} in the last straight cell).  SETTLED = the first
 *          tick whose last solved plan has X0..XM in the hover box (speed
 *          <= BS_HV_SPEED_ARC, |e_l| and |e_c| <= BS_HV_BAND in the hover
 *          frame, |a| <= BS_HV_A_SETTLE): bs_arc_tick returns 2 on that
 *          tick (the handback); the driver never sets `finished` by cell
 *          index on a hover chart.  The caller (tests/arc_gate.c --hover,
 *          the firmware) runs M more ticks after settle, and treats
 *          BS_HV_SETTLE_CAP ticks after arrival unsettled as a failed run.
 * Charts built with hover = 0 reach none of these branches: the v1 / v2 /
 * v2b / v2c regression charts are bit-identical to before.
 *
 * TUBE FILLETS (v2e = v2d + tube fillets, 2026-09-27): a chart built with
 * tube_delta > 0 (bs_arc_chart_build: BS_TUBE_DELTA = 0.4 m when
 * BS_ARC_TUBE, the configuration flown since v2e; bs_arc_chart_build_tube
 * takes the delta) sizes every fillet with that inside clearance and
 * gives its chords the inner-corner / outer-boundary corridor rows, all
 * computed once at the chart build (see bs_arc_chart_build below).  The
 * tick path is unchanged: bs_arc_margins reads the two stored margins.
 * tube_delta = 0 (bs_arc_chart_build_rows, and every chart before v2e)
 * is bit-identical to the v1 .. v2d charts.
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
#define BS_ARC_HOVER 1               /* v2d: bs_arc_chart_build ends the chart
                                      * in the hover cell (no tail)          */
#ifndef BS_ARC_TUBE
#define BS_ARC_TUBE 1                /* v2e: bs_arc_chart_build sizes the
                                      * fillets with the inside clearance
                                      * BS_TUBE_DELTA and gives the chords
                                      * the tube rows (0: the v2d chart; a
                                      * firmware build takes --define
                                      * BS_ARC_TUBE=0 at configure)        */
#endif
/* the v2d block (BS_HV_*, bs_rc_arc_Ph / Kh) is emitted into the v2
 * header only; the v1 header (byte-identical to its record) has none, so
 * a v1 build refuses a hover chart (bs_arc_chart_build_rows -> -4) and
 * bs_arc_chart_build stays the tail chart there */
#ifdef BS_HV_SPEED_ARC
#define BS_ARC_HV_AVAILABLE 1
#else
#define BS_ARC_HV_AVAILABLE 0
#endif
/* the v2e constant BS_TUBE_DELTA (the owner's delta = 0.4 m) is likewise
 * emitted into the v2 header only; without it bs_arc_chart_build builds
 * no tube (bs_arc_chart_build_tube still takes any delta explicitly) */
#ifdef BS_TUBE_DELTA
#define BS_ARC_TUBE_AVAILABLE 1
#else
#define BS_ARC_TUBE_AVAILABLE 0
#endif

enum { BS_ARC_TA = 0, BS_ARC_CA = 1, BS_ARC_FILLET = 2, BS_ARC_HV = 3 };

typedef struct {
    double px, py;        /* start point (N, E)                              */
    double thx, thy;      /* unit heading; n_hat = (-thy, thx)               */
    double l;             /* cell length; chart pace = l / BS_TS             */
    double chi;           /* seam angle at the cell START, deg (0 = none)    */
    double d_in;          /* fillet cell: chord-midpoint depth into the corner */
    double tube_in;       /* v2e fillet chord (chart tube on): b_k, the
                           * inward corridor margin s e_c <= b_k (the
                           * distance from the chord line to the corridor's
                           * inner corner C_in); 0 otherwise                */
    double tube_out;      /* v2e fillet chord: o_k, the outward margin
                           * -s e_c <= o_k (closed-form exit to the
                           * corridor's outer boundary, the smaller over
                           * the chord's two endpoints); 0 otherwise        */
    signed char kind;     /* BS_ARC_TA / BS_ARC_CA / BS_ARC_FILLET / BS_ARC_HV */
    signed char sense;    /* fillet: +1 / -1 (turn sense), else 0           */
    short leg;            /* polyline leg index                              */
} bs_arc_cell;

typedef struct {
    bs_arc_cell *cells;
    int n_cells;
    int k_end;            /* arrival cell = first trailing cell (tail chart)
                           * / the hover cell = the last cell (hover chart)  */
    int n_v;              /* vertices kept                                   */
    double vx[BS_ARC_MAX_WP], vy[BS_ARC_MAX_WP];
    double chi_v[BS_ARC_MAX_WP];       /* turn angle at each vertex, deg     */
    int leg_first[BS_ARC_MAX_WP];      /* first cell of leg i's straight     */
    int leg_last[BS_ARC_MAX_WP];       /* last cell (incl. fillet) of leg i  */
    double v_trim, v_leg, corner_pace, w;
    int cell_rows;        /* 1 (default): every cell its own row family, the
                           * velocity rows about the cell's own pace, straight
                           * cell counts rounded up (v2c); 0: the family-
                           * centred rows and round() counts of the
                           * p4-current / p4-v2 / p4-v2b archives           */
    int hover;            /* 1: no tail, the chart ends in the hover cell
                           * (kind BS_ARC_HV) at the last waypoint (v2d);
                           * 0: the 60 m trailing straight (v1 .. v2c)      */
    int tube;             /* 1: tube fillets (v2e): the chords carry
                           * tube_in / tube_out and bs_arc_margins uses
                           * them; 0: the fillets and rows of v1 .. v2d     */
    double tube_delta;    /* the inside clearance the fillets were sized
                           * with, m (0 when tube = 0)                      */
} bs_arc_chart;

/* Upper bound on the cell count / block size for n_wp waypoints of total
 * polyline length len_m (call before allocating). */
int bs_arc_chart_cells_bound(int n_wp, double len_m);
size_t bs_arc_chart_size(int n_wp, double len_m);

/* Build the fillet chart.  `mem` holds the cells (>= bs_arc_chart_size).
 * bs_arc_chart_build builds the configuration flown: cell_rows = 1,
 * hover = BS_ARC_HOVER (v2d: the chart ends in the hover cell) and the
 * tube fillets at BS_TUBE_DELTA when BS_ARC_TUBE (v2e; the v2 header
 * carries the constant, the v1 header does not: no tube there);
 * bs_arc_chart_build_rows takes the two flags (cell_rows 0 = the archives'
 * family rows and cell counts; hover 0 = the 60 m tail of the v1 .. v2c
 * archives) with the tube OFF (the v1 .. v2d charts, bit for bit);
 * bs_arc_chart_build_tube takes the inside clearance as well (0 = off).
 * The cell count depends on the mode, so it is a build-time choice, not a
 * field to flip afterwards.  Returns 0, or a negative code (-1 waypoints,
 * -2 memory, -3 geometry, -4 hover without cell rows, -5 tube: delta not
 * in [0, BS_ARC_W), a chord endpoint outside the corridor, or a
 * non-positive tube margin -- stage2_arc_clock_mission.build_cells_fillet
 * raises ValueError on the same charts).
 *
 * TUBE FILLETS (v2e = v2d + tube fillets, owner decision 2026-09-27,
 * study tests/fix1b_tube_2026-09-27; Python of record
 * stage2_arc_clock_mission.build_cells_fillet(tube_delta) and
 * tube_outer_reach, arc_v2.fillet_rows).  tube_delta > 0: every fillet
 * radius is r_f = (w - tube_delta) / (1 - cos h), h = |chi| / 2 (the
 * tangent-length cap and the chord-count rule unchanged), and every chord
 * cell k carries, computed here once per chart (closed form, no search):
 *   tube_in  b_k = n_k . (C_in - p_k), n_k = sense * (-thy, thx) the
 *            inward chord normal, p_k the chord start, C_in = V + (nA +
 *            nB) * (w / (1 + tA . tB)) the corner of the two legs' inner
 *            offset lines at w (not w - delta);
 *   tube_out o_k = the smaller over the chord's two endpoints of the exit
 *            distance along -n_k from the union of the two leg capsules
 *            of radius w (the outer offset lines joined by the radius-w
 *            cap about the vertex): per capsule the hull of the ray's
 *            slab interval on the rectangle and its two disc intervals
 *            (t = -beta +- sqrt(beta^2 - gamma)).
 * bs_arc_margins then writes s e_c <= b_k, -s e_c <= o_k on the chord in
 * place of max(w - d_in, 0.05) / w + d_in.  Straight cells, the hover
 * cell, all other rows, the reward, the solver and rule G1 are
 * unchanged.  tube_delta = 0: r_f = (w - 0) / (1 - cos h) = w / (1 - cos
 * h) bit for bit and no chord carries a tube margin. */
int bs_arc_chart_build(const double *vx, const double *vy, int n_wp,
                       void *mem, size_t mem_len, bs_arc_chart *out);
int bs_arc_chart_build_rows(const double *vx, const double *vy, int n_wp,
                            void *mem, size_t mem_len, bs_arc_chart *out,
                            int cell_rows, int hover);
int bs_arc_chart_build_tube(const double *vx, const double *vy, int n_wp,
                            void *mem, size_t mem_len, bs_arc_chart *out,
                            int cell_rows, int hover, double tube_delta);

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
    /* terminal-pair cache: key = cell family key (see bs_arc_fam_key);
     * a plain ring under cell_rows, else slots 0/1 stay resident */
    struct { int key; double P[BS_NX * BS_NX], K[BS_NU * BS_NX]; }
        dare[BS_ARC_DARE_CACHE];
    int dare_n, dare_next;
    int dare_fail;
    int dare_solves;      /* DARE pairs solved (cache misses), cumulative */
    int dare_hover;       /* terminal pairs served from the emitted hover
                           * pair (bs_rc_arc_Ph / Kh), cumulative; never a
                           * cache miss                                    */
    /* D8: the on-target DARE solves of the current tick, microseconds by
     * the driver's clock hook (0 without one); the driver copies `micros`
     * in at the start of every tick and reads us_dare out at its end */
    unsigned long (*micros)(void);
    unsigned long us_dare;
    /* v2d: the per-stage progress reward weight of the rendered horizon,
     * w_t for the term on delta of X_{t+1}: BS_HV_W when ks[t + 1] is the
     * hover cell, BS_SPEED_TILT_W otherwise (bs_schedule.tilt_w) */
    double tilt_w[BS_N];
    /* the world-frame jerk-change penalty weight (arc_v2.W_DJ, runtime
     * port 2026-11-13): copied into bs_schedule.w_dj by every
     * bs_arc_render.  bs_arc_driver_init zeroes it (v2e bit for bit); the
     * caller sets d->S->w_dj AFTER bs_arc_driver_init and before
     * bs_arc_engage (the firmware: BSLV_WDJ).  BS_NU == 2 only. */
    double w_dj;
    /* v2d: the stop state the assignment consults (arc_v2.HoverStop):
     * gate_open / latch are the rule's state (the latch a running minimum
     * over EVERY assign call, -1 = none); forced, t_nat, n_vtail describe
     * the LAST bs_arc_assign call (stages forced into the hover cell, the
     * first stage > 0 whose natural cell is >= k_end or -1, virtual-tail
     * stages).  All zero / -1 on a tail chart. */
    struct { int gate_open, latch, forced, t_nat, n_vtail; } hv;
} bs_arc_sched;

/* the family key of a cell: under cell_rows 2 + k for every cell; else
 * 0 = ta, 1 = ca, 2 + k for a fillet cell k */
int bs_arc_fam_key(const bs_arc_chart *ch, int k);

/* margins of cell k's family (no trust), BS_NROW entries in template
 * order (the velocity rows about the cell's own pace under cell_rows,
 * about the family pace otherwise; every margin is positive on every
 * chart the builder produces -- the gate asserts it) */
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

/* Render the schedule (margins with trust, terminal pair, the per-stage
 * reward weights) after bs_arc_assign.  Returns 0, or -1 if the tip
 * family's DARE failed. */
int bs_arc_render(const bs_arc_chart *ch, bs_arc_sched *S, double trust);

/* ---- the driver ------------------------------------------------------ */
/* the brake seed of record (arc_v2.brake_seed's info, rule G1): written
 * once, at the tick the gate opens with a seed; tick = -1 when none (the
 * gate opened by the vehicle entering the hover cell, or is still
 * closed).  t_b = the onset stage, n_stop = n_r + n_h + n_r2 the seed's
 * stages (ramp in at -jerk, hold at -a_s, ramp out), rest_stage = t_b +
 * n_stop (may exceed N: stages >= N are not written); d_stop = the seed's
 * own rest distance, rest_short = dist_wp - d_stop; a_req = v_b^2 / (2
 * dist_wp), a_nom = the stop of record's mean deceleration, a_s_cont =
 * the continuous-time root, a_s the grid re-solve, a_min = -a_s, j_max =
 * the seed's largest |j_l|. */
typedef struct {
    int tick, t_b, n_stop, rest_stage, n_r, n_h, n_r2;
    double d_stop, v_b, a_b, dist_wp, rest_short;
    double a_req, a_nom, a_s_cont, a_s, a_min, j_max;
} bs_arc_seed;

/* brake_seed / tick_gate outcomes (the no-root reasons of
 * arc_v2.trapezoid_closed_form; LONG is a C-side guard on a non-finite or
 * absurd seed duration, unreachable on a qualifying stage) */
enum {
    BS_ARC_SEED_NONE = 0,     /* no stage qualifies: gate stays closed   */
    BS_ARC_SEED_OK = 1,       /* seeded: gate opens                      */
    BS_ARC_SEED_DPRIME = 2,   /* d' <= 0                                 */
    BS_ARC_SEED_DISC = 3,     /* the quadratic's discriminant < 0        */
    BS_ARC_SEED_AB = 4,       /* a_s > AB                                */
    BS_ARC_SEED_RAMPUP = 5,   /* a_s < -a_b (a ramp up would be needed)  */
    BS_ARC_SEED_LONG = 6      /* duration beyond 4 N stages (C guard)    */
};

/* the stop model of a hover chart (arc_v2.stop_params, computed ONCE at
 * engage as HoverStop.__init__ does): the last straight cell's pace v_c,
 * its family's rear face AB (the bound on rule G1's plateau a_s) and
 * forward face AF, the symmetric stop of record from v_c
 * (brake_profile_exact: n_stop stages, distance d_stop, jerk J, peak
 * deceleration a_max, residual speed v_residual = 0) and rule G1's
 * nominal deceleration a_nom = v_c^2 / (2 d_stop).  The Python
 * refuses a chart without that stop (AssertionError); bs_arc_engage
 * returns -4 on it.  valid = 1 once computed (hover chart only).
 * dist_cp / dist_stride: checkpoints of rule G1's distance table (arc_v2.
 * dist_suffix, S[c] = S[c + 1] + l_c accumulated backward from WP_N):
 * dist_cp[i] = S[i * dist_stride] for i * dist_stride < k_end, dist_stride
 * = max(1, ceil(k_end / (BS_HV_DIST_CP - 1))), built once at engage.  The
 * onset scan walks down from the checkpoint at or above a stage's cell,
 * the same additions in the same order as the Python's table, so a
 * closed tick costs at most dist_stride - 1 + (the cells the horizon
 * spans) additions (O(N) for charts up to (BS_HV_DIST_CP - 1) N cells;
 * verts: 967 cells, stride 16) with 8 B per checkpoint, not per cell. */
#define BS_HV_DIST_CP 64
typedef struct {
    int valid, n_stop;
    double AB, AF, v_c, d_stop, J, a_max, v_residual, a_nom;
    int dist_stride;
    double dist_cp[BS_HV_DIST_CP];
} bs_arc_stop_model;

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
    /* ---- v2d: the stop at the last waypoint (hover chart only; every
     * field reset at engage, -1 / 0 on a tail chart) ---- */
    int gate_open;            /* the gate of rule G1                       */
    int latch;                /* the latch as the LAST SOLVED plan of the
                               * tick saw it (the value of the last round's
                               * assignment, -1 none); S->hv.latch is the
                               * running state                             */
    int open_tick;            /* tick index at which the gate opened, -1  */
    double open_speed;        /* HoverStop.open_speed: the vehicle's speed
                               * (pace + delta) at the gate tick, m/s      */
    double open_dist;         /* HoverStop.open_dist = the seed record's
                               * dist_vehicle: the vehicle's distance to
                               * WP_N along the chart at the gate tick, m */
    int arrival_tick;         /* index of the first state sample in the
                               * hover zone, -1: the engage state (0) when
                               * it is already in the zone (arc_v2._run
                               * tests in_zone on sample 0), else the
                               * post-step state kw + 1                    */
    int settle_tick;          /* the tick whose plan first settled, -1     */
    int settled;              /* 1 from the settle tick on                 */
    int forced;               /* stages forced into the hover cell by the
                               * last solved plan's assignment             */
    int n_hv;                 /* stages (of 0..N) of the last solved plan's
                               * schedule in the hover cell                */
    bs_arc_seed seed;         /* the brake seed record                     */
    int no_root;              /* rule G1: closed-gate ticks whose scan
                               * qualified a stage but whose quadratic had
                               * no admissible root (HoverStop.no_root)    */
    int no_root_tick, no_root_tb, no_root_reason;  /* the first such tick,
                               * its t_b and BS_ARC_SEED_* reason, -1 / 0  */
    bs_arc_stop_model stop;   /* the stop model of record (at engage)     */
    unsigned us_dare;         /* D8: on-target DARE solves this tick, us
                               * (copied from S->us_dare on EVERY chart,
                               * tail charts included)                     */
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
 * inputs.  On a hover chart the stop model of record (d->stop) is
 * computed first.  Returns 0 on success, -1 when the engage's render
 * failed, -4 on a hover chart without the symmetric stop of record
 * (arc_v2.stop_params' AssertionError), -5 on an invalid jerk-change
 * weight d->S->w_dj (negative / NaN, or non-zero on the three-input v1
 * model). */
int bs_arc_engage(bs_arc_driver *d);

/* One maintained tick, in the order of ac.run / arc_v2._run: publish (the
 * caller lifts before calling), plant step with the clock rule (arrival =
 * zone entry recorded on the post-step state), shift-append, [hover
 * chart: the gate test once, HoverStop.tick_gate], `rounds` x (assign,
 * render, pinned solve), [hover chart: the settle predicate on the last
 * solved plan's X0..XM], commit stage M-1.  u0 receives the applied
 * input.  Returns
 *    0  ticked (an input applied);
 *    1  tail chart only: arrived, k >= k_end, nothing applied (finished);
 *    2  hover chart only: ticked, and this is the tick at which `settled`
 *       first became 1 (later ticks return 0 with settled staying 1);
 *   <0  failure (-1 not engaged, -2 DARE / render, -3 solve). */
int bs_arc_tick(bs_arc_driver *d, double *u0);

/* Positions of the committed prefix's nodes (M + 1) from the current
 * state, in the plane (N, E); also the cell of each node. */
void bs_arc_lift(const bs_arc_driver *d, double *pn, double *pe, int *kn);

#endif /* BS_ARC_H */
