#include "AP_BSolver.h"

#if AP_BSOLVER_ENABLED

#include <AP_AHRS/AP_AHRS.h>
#include <AP_Mission/AP_Mission.h>
#include <GCS_MAVLink/GCS.h>
#include <AP_Logger/AP_Logger.h>
#include <stdlib.h>
#include <string.h>

// RUNTIME TABLES (2026-09-02): the flight build bakes NO mission.  The
// model constants come from bs_model_runtime_const.h (selected as
// BS_DATA_HEADER in wscript) and the mission tables are built on the
// aircraft at upload time by AP_BSMissionBuilder / bs_mission_build.c.
extern "C" {
#include "bsolver/core/bs_solver.h"
#include "bsolver/core/bs_mission_build.h"
}
#include "AP_BSMissionBuilder.h"

extern const AP_HAL::HAL &hal;

// Arena, sized from the core's own requirement rather than guessed.  The
// first version reserved 64*1024 doubles (524,288 B) — about 2.6x what
// bs_workspace_size() asks for, and more static RAM than the board has heap.
// ONE arena for the whole library, taken from the HEAP at init rather than
// reserved in .bss.  Static arrays for this overflowed the board's RAM region
// at link (0x300751e0 against a 0x30040000 boundary) — csolver puts its own
// ~288 kB arena on the heap for the same reason.  The bench and the live
// solve never run together, so they share it.
static double *bs_arena;
static size_t bs_arena_words;

// THE PUBLISHED MISSION TABLES.  Written by the builder on the solver
// thread while the solver is inactive, read by everything below.  NULL
// until the first successful build; the READY gate requires it.
static const bs_mission_tables *mt;

// THE MISSION BLOCK, reserved ONCE at init beside the arena.  Measured
// on the CubeOrange+: a per-upload calloc of the worst-case table block
// (129 kB) FAILS at runtime -- 330 kB are free but fragmented across
// the H7's RAM regions, and no single region still holds it once the
// system is up.  At init, right after the arena, the big region is
// intact.  The builder and the bench both build into this reservation.
static void *bs_mb_block;
static size_t bs_mb_block_len;
static size_t bs_mb_largest;      // largest allocatable block, probed at init (bounded by worst + margin)
static AP_BSMissionBuilder bs_mb;
static AP_Mission *bs_mb_mission;      // captured at poll time

// The ONE bs_problem the driver owns.  It is file-scope so ensure_problem()
// can keep it across ticks: the reference loop assembles at tau at the top of
// a tick and at (k+1) - o_off at the bottom, and those are the SAME phase one
// tick apart, so caching turns two assemblies per tick into one.  On target
// bs_problem_init is 16.8 % of a pinned cycle, so this is not cosmetic.
static bs_problem bs_driver_problem;

// ---------------------------------------------------------------- ingress
//
// BSLV_INGRESS: the approach to the path start becomes the FIRST LEG of the
// anytime problem (the negative-tau design of
// scripts/mpc/bsolver/doc/INGRESS-CAMPAIGN.md, with the probe-verified
// simplifications: plain row-family cells at trim pace, corridor pair
// opened to 5 m, no ramp cells, no geometric refusals).  Everything here is
// RUNTIME-BUILT at the ingress engage: one straight leg of ing_n cells from
// the hover point to the path start, the junction seam (rotation + affine
// kick, exactly off_vec(V, V, chi) of the host model) at destination tick
// 0, and the first N body schedule entries copied so the horizon can span
// the seam.  For tau >= 0 every lookup reads the untouched flash tables:
// the body problem sequence is the certified one BY CONSTRUCTION, which is
// ground rule 3 of the campaign.
//
// The prefix is capped at BS_N cells (~85 m): beyond that the horizon tail
// could land on the leg, whose runtime family has no bs_P / bs_K entry.
static const int32_t ING_MAX = BS_N;
static const int32_t ING_LEN = ING_MAX + BS_N + 1;
static int32_t ing_n;                       // 0 = classic rest-engage
static double ing_anchor_n, ing_anchor_e;   // p_sched(-ing_n), NED xy
static double ing_dir_n, ing_dir_e;         // leg tangent (unit)
static double ing_ang_deg;                  // junction seam angle chi
static int ing_fam[ING_LEN], ing_rot[ING_LEN], ing_off[ING_LEN];
static double ing_rows[BS_NROW * 10];       // family t, corridor pair at 5 m
static double ing_arc[ING_MAX + 1];         // published arc per leg tick
static double ing_rotmat[BS_NX * BS_NX];    // T_switch(chi)
static double ing_offvec[BS_NX];            // off_vec(V_TRIM, V_TRIM, chi)
static bs_schedule bs_sched_ingress;
static const double ING_COR_WIDE = 5.0;     // leg corridor half-width, m
static uint32_t ing_climb_ok_ms;            // climb-complete timer start
static uint32_t ing_txt_ms;                 // refusal statustext throttle

// ---------------------------------------------------------------- WCET bench
// Compile-gated, default OFF, mirroring MPC_BENCH_WCET in AP_MPCSolver.cpp.
// NEVER leave enabled for a flight build.
//
// IT MUST BE SEEDED AT AN ENGAGE ITERATE.  On the maintained domain the
// iterate is already at its fixed point, the "predicted decrease below noise"
// guard fires, no step is taken, and the solve costs the same regardless of q
// — 873 of the mission's 931 ticks look like that.  Measured on the host:
// engage q=1/2/3 = 1.173/2.222/3.187 ms (scaling as it should), maintained
// q=1/2/3 = 1.161/1.161/1.159 ms (all guard).  Timing the maintained state
// would understate the q=3 worst case by 2.75x.

// The committed prefix length, M.  PLAN_NODES is M + 1 and the core's
// BS_M_COMMIT is the same number arriving from the generated model table;
// they are checked against each other rather than kept in step by hand.
static const int BS_M = BS_M_COMMIT;
static_assert(AP_BSolver::PLAN_NODES == BS_M_COMMIT + 1,
              "PLAN_NODES must be the committed prefix plus its terminal node");

// The driver's fixed-size members are declared in AP_BSolver.h, which does
// NOT see the generated model header (see bs_defs.h for why).  This TU does,
// so it is where the mirror is checked.  A model/reference pair that does not
// match the mirror -- e.g. BS_DATA_HEADER switched to the corner tables while
// BS_DRV_NV still says 120 -- fails HERE rather than overrunning _U at run
// time.
static_assert(BS_DRV_NV == BS_NV,
              "BS_DRV_NV (AP_BSolver.h buffers) must equal BS_NV of the "
              "selected BS_DATA_HEADER");
static_assert(BS_DRV_M == BS_M_COMMIT,
              "BS_DRV_M must equal BS_M_COMMIT of the selected "
              "BS_DATA_HEADER");
// The old model/reference cross-header asserts are gone with the flash
// reference header; the equivalent invariants are now runtime gates in
// bs_mission_build (G1-G10), checked per uploaded mission.

// The airframe drag coefficient the replay's lean map feeds forward,
// k/m in 1/m — the same SN77 number AP_MPCSolver streams (MPC_K_OVER_M).
// Handing 0 instead leaves the ~2.3 m/s^2 of cruise drag entirely to the
// velocity PID's integrator, which is a lag the plan does not model.
static const float BS_DRAG_K = 0.018017f;

const AP_Param::GroupInfo AP_BSolver::var_info[] = {
    // @Param: ENABLE
    // @DisplayName: Enable the onboard certified anytime MPC
    // @Values: 0:Disabled,1:Enabled
    // @User: Advanced
    AP_GROUPINFO_FLAGS("ENABLE", 1, AP_BSolver, _enable, 0, AP_PARAM_FLAG_ENABLE),

    // @Param: Q
    // @DisplayName: Newton iterations per 4 Hz solve
    // @Range: 1 3
    // @User: Advanced
    // q = 2 IS THE CADENCE OF RECORD for the corner-online configuration,
    // and the default is 2 for that reason and not as margin.  On this
    // design q = 1 is MEASURED UNUSABLE: the free tail does not keep up with
    // the re-timing kick through the corner windows, the corridor is missed
    // by 2-4x and the in-window |e_c| runs to 2.69-3.74 m against the 2.0 m
    // band.  At q = 2 the host C mission closes it (|e_c| 1.9618, |e_l|
    // 2.7721) and the WCET fits: 132.2 ms on the Cube against the 250 ms
    // tick (worst of two boots on the deployed anisotropic configuration;
    // 39-47 % margin after the x1.16 whole-tick calibration).
    //
    // The record's N = 40 ANCHORED configuration read the other way round
    // (uniform q >= 1 certified on the interleaved cadence,
    // wp5_h1_interleaved_q.md, varrho = 4.220230e-2).  That statement is
    // about a different clock and a different horizon; do not carry it over.
    AP_GROUPINFO("Q", 2, AP_BSolver, _q, 2),

    // @Param: QE
    // @DisplayName: Cap on the engage seed's full-horizon Newton steps
    // @Range: 1 200
    // @User: Advanced
    // The interleaved loop of record seeds from an EXACT solve at tau = 0
    // (newton(..., maxit=300), gradient tolerance 1e-7); there is no
    // commitment to pin before it.  The firmware iterates the same way and
    // stops at the same tolerance, but caps the work so a non-converging
    // seed cannot hang the mission clock.  Measured on the host: the exact
    // seed converges in 32 steps on the CORNER-ONLINE tables (27 on the
    // record's), and the default 60 clears both with room.  Measured on the
    // corner host C mission at q = 2, that seed reproduces 983 ticks,
    // arrival tick 971 (242.75 s), offset +55, 55 re-timings.
    AP_GROUPINFO("QE", 3, AP_BSolver, _qe, 60),

    // @Param: IR
    // @DisplayName: Restrict the per-tick setup to the pinned rows
    // @Values: 0:Build the whole quadratic half,1:Build only rows >= 3*M
    // @User: Advanced
    // bs_problem_init does NOT shrink when the SOLVE is pinned -- it still
    // assembles the full BS_NV x BS_NV quadratic half -- on target that setup
    // is a third of the pinned cycle.  Rows below 3*M are dead on the pinned
    // path (the value never reads Hquad/Cquad and the derivatives visit only
    // rows >= 3*M), so not writing them is EXACT: measured 0.000e+00 on the
    // tail gradient and the tail Hessian over 21 draws / 7 phases, with the
    // core refusing an npin = 0 evaluation of a restricted problem rather
    // than silently returning zeros.  Kept switchable because it is the one
    // thing here that changes bs_problem_init's contract.
    AP_GROUPINFO("IR", 4, AP_BSolver, _ir, 1),

    // (param index 5 was BSLV_INGRESS, retired 2026-09-02: the ingress
    // leg is now the unconditional AUTO approach — do not reuse the slot)

    AP_GROUPEND
};

AP_BSolver::AP_BSolver()
{
    AP_Param::setup_object_defaults(this, var_info);
}

// ---------------------------------------------------------------- geometry
//
// Everything below reads GENERATED mission data only.  There is no projection
// and no measured position anywhere in this file's control path: the
// reference is a function of the mission clock and the model state, which is
// what premise A4 requires.

// p_sched(tau): the reference position AT schedule tick tau.  A table read,
// not an arc-length interpolation.
//
// The predecessor of this function searched bs_ref_sched for an arc length
// s_veh and interpolated LINEARLY between the two path samples bracketing it.
// That is not the object the quotient is defined against, and on the
// corner-online tables it is not even close: the published in-window
// reference is the centreline polyline sampled at 25 stages, so consecutive
// samples turn by 72 deg and then 108 deg, and an interpolation that walks
// along those chords reproduces the kinks verbatim.  Measured on the flown
// mission, that lift demanded 70.86 deg of lean (p95 58.66) against a 35 deg
// authority and put the vehicle 3.11 m outside the 2 m corridor.
//
// The reconstruction below is instead the one the certificates and the
// successor simulator state the quotient in -- sim_corner_online.py's
// window_track / published_track -- and it needs only the sample at tau.
static void ref_sample(int32_t tau, double *x, double *y)
{
    if (mt == nullptr) { *x = 0.0; *y = 0.0; return; }
    if (tau < 0) {
        if (ing_n > 0) {
            // The ingress leg publishes a JERK-SHAPED REST-TO-TRIM RAMP, not
            // uniform trim cells.  p_sched is a table read that never waits,
            // so a trim-paced leg publishes a reference running at >= 11.4
            // m/s from the engage instant while the vehicle hovers -- the
            // classic engage has 79.6 m of straight row to chase that down
            // (measured: the SITL vehicle overspeeds to ~18 m/s and
            // converges before the first window), but a short leg puts the
            // JUNCTION CORNER inside the chase window (measured: 5-9 m of
            // corner cut).  Ramping the published arc keeps node 0 at rest
            // and the chase closed; the model chart is untouched (its own
            // slow junction crossing is a certified-probe result), and the
            // velocity/acceleration feed-forward stays consistent because
            // the lift differentiates these positions.
            const int32_t k = (tau < -ing_n) ? 0 : (int32_t)(ing_n + tau);
            const double a = ing_arc[k];
            *x = ing_anchor_n + a * ing_dir_n;
            *y = ing_anchor_e + a * ing_dir_e;
            return;
        }
        tau = 0;
    }
    int32_t i = (tau > mt->n_path - 1 ? mt->n_path - 1 : tau);
    if (mt->sched.ring_mask) i &= mt->sched.ring_mask;
    *x = mt->path[2 * i];
    *y = mt->path[2 * i + 1];
}

// CHART frame at schedule tick tau: the heading the quotient's deviations are
// expressed in.  It is NOT the finite-difference tangent of the polyline —
// the chart turns only at the frame seams, which is exactly why the seam
// angles exist as their own table.
static void ref_frame(int32_t tau, double *tx, double *ty)
{
    if (mt == nullptr) { *tx = 1.0; *ty = 0.0; return; }
    if (tau < 0) {
        if (ing_n > 0) {              // the leg's chart heading is constant
            *tx = ing_dir_n;
            *ty = ing_dir_e;
            return;
        }
        tau = 0;
    }
    int32_t i = (tau > mt->n_path - 1 ? mt->n_path - 1 : tau);
    if (mt->sched.ring_mask) i &= mt->sched.ring_mask;
    const double a = mt->psi[i];
    *tx = cos(a);
    *ty = sin(a);
}

static double ang_at(int32_t tau)
{
    if (mt == nullptr) return 0.0;
    if (ing_n > 0) {
        // the junction seam sits at DESTINATION tick 0: the transition from
        // the last leg tick onto the body chart.  bs_ref_ang[0] is 0 on this
        // mission (it starts on a row), so serving chi here rewrites nothing
        // certified; it also makes ledger_step seam-suppress the junction
        // tick, exactly as the host probe's schedule does.
        if (tau == 0) return ing_ang_deg;
        if (tau < 0) return 0.0;      // the leg is straight
    }
    int32_t i = (tau < 0) ? 0
              : (tau > mt->n_clk - 1 ? mt->n_clk - 1 : tau);
    if (mt->sched.ring_mask) i &= mt->sched.ring_mask;
    return mt->ang[i];
}

#if BS_STAGE2_PACE_TAB
// Stage 2 (c): the published pace at a mission tick (ring-aware, like
// ang_at); negative on the ingress leg, whose pace is the ramp's own.
static double pace_at(int32_t tau)
{
    if (mt == nullptr) return -1.0;
    if (ing_n > 0 && tau < 0) return -1.0;
    int32_t i = (tau < 0) ? 0
              : (tau > mt->n_clk - 1 ? mt->n_clk - 1 : tau);
    if (mt->sched.ring_mask) i &= mt->sched.ring_mask;
    return mt->pace[i];
}
#endif

static int fam_at(int32_t tau)
{
    // The ingress family is deliberately NOT the trim family: ledger_step's
    // family gate then cannot fire on the leg (the campaign's D4 gating,
    // confirmed as the intended mechanics by the wide-leg probe).
    if (mt == nullptr) return BS_MB_FAM_T;
    if (tau < 0 && ing_n > 0) return -1;
    int32_t i = (tau < 0) ? 0
              : (tau > mt->n_clk - 1 ? mt->n_clk - 1 : tau);
    if (mt->sched.ring_mask) i &= mt->sched.ring_mask;
    return (int)mt->sfam[i];
}

// THE AFFINE PART OF THE SEAM MAP, c(tau), applied at the DESTINATION tick.
//
// The corner-online reference pace is piecewise constant across a window, so
// the transition on the phase quotient is affine rather than linear:
//
//     xi_{t+1} = T(t+1) (A xi_t + B u_t) + c(t+1)
//
// with the offset added AFTER the rotation.  This is off_of() in
// bsolver/core/bs_solver.c, written out for the driver's own plant step: same
// tables (bs_off / bs_mission_off, visible here through bs_solver.h), same
// destination-tick indexing, same clamp as the core's tick_of() for a
// non-periodic schedule.
//
// IT IS NOT OPTIONAL AND IT IS NOT SMALL.  On the flown clock 76 of the 966
// ticks carry a non-zero c, 38 of them with an IDENTITY rotation — so a
// rotation-only plant step misses them entirely, and the largest entry is
// 11.27 m/s of along-track pace.  The record's mission has no affine part at
// all (bs_mission_off is not emitted), which is why the predecessor of this
// function was correct there; BS_NOFF is what distinguishes the two.
static const double *off_at(int32_t tau)
{
    static const double zero[BS_NX] = { 0.0 };
    if (mt == nullptr) return zero;
    if (ing_n > 0) {
        if (tau == 0) return ing_offvec;   // the junction's affine kick
        if (tau < 0) return zero;          // constant pace on the leg
    }
    int32_t i = (tau < 0) ? 0
              : (tau > mt->n_clk - 1 ? mt->n_clk - 1 : tau);
    if (mt->sched.ring_mask) i &= mt->sched.ring_mask;
    return &mt->off[(size_t)mt->soff[i] * BS_NX];
}

// One exact-ZOH plant step followed by the chart re-expression at the solve
// boundary:  xi <- T(tau+1) (A xi + B u).  This is A6/B6 of the certificates
// written out (verified entry for entry against the exported matrices), and
// T_switch's rotation of the three tangential/normal pairs (0,4), (1,5),
// (2,3) by the seam angle.
static void plant_step(double *xi, const double *u, double seam_deg,
                       const double *off)
{
    const double Ts = BS_TS;
    const double d0 = xi[0], a0 = xi[1], el0 = xi[2];
    const double ec0 = xi[3], ed0 = xi[4], an0 = xi[5];
    const double j = u[0], nu = u[1], jn = u[2];

    double y0 = d0 + Ts * a0 + 0.5 * Ts * Ts * j;
    double y1 = a0 + Ts * j;
    double y2 = el0 + Ts * d0 - Ts * nu;
    double y3 = ec0 + Ts * ed0 + 0.5 * Ts * Ts * an0 + (Ts * Ts * Ts / 6.0) * jn;
    double y4 = ed0 + Ts * an0 + 0.5 * Ts * Ts * jn;
    double y5 = an0 + Ts * jn;

    // Exact zero test, written as |x| > 0 because the tree builds with
    // -Werror=float-equal.  The seam table holds exact 0 / +-72 / +-108, so
    // this IS the reference's `if a == 0.0: return I`, not a tolerance.
    if (fabs(seam_deg) > 0.0) {
        const double c = cos(radians(seam_deg)), s = sin(radians(seam_deg));
        double t;
        t = y0; y0 = c * t + s * y4; y4 = -s * t + c * y4;
        t = y1; y1 = c * t + s * y5; y5 = -s * t + c * y5;
        t = y2; y2 = c * t + s * y3; y3 = -s * t + c * y3;
    }
    // ...and the affine part, AFTER the rotation, indexed by the destination
    // tick.  All-zero for a schedule without one, so the record build takes
    // the same path and gets the same answer.
    xi[0] = y0 + off[0]; xi[1] = y1 + off[1]; xi[2] = y2 + off[2];
    xi[3] = y3 + off[3]; xi[4] = y4 + off[4]; xi[5] = y5 + off[5];
}

// THE ASSIGNMENT LEDGER OF RECORD, verbatim from wp5_anytime_sim.py's
// run_mission() (`assignment: rigid clock, one integer offset per solve,
// frozen trigger 1.5 cell_min = 2.106 m, deferred at frame seams`):
//
//     p_phase += Ts * xi[0]                        (done by the caller)
//     if |p_phase| > HYST and ang(tau) == 0 and ang(tau+1) == 0:
//         cell = CELL_T if fam(tau) == 't' else CELL_MIN
//         p_phase < 0 ? (o += 1, p_phase += cell) : (o -= 1, p_phase -= cell)
//
// It is driven by the MODEL state, which is the only thing that makes it
// consistent: the accumulator credits a fixed trim cell per re-timing, and a
// trim cell is what the MODEL clock advances by.  Fed a MEASURED arc length
// instead it is not consistent — the mission's true mean cell is 2.4073 m,
// which over 916 intervals is a 140-tick advance deficit no accumulator of
// Ts*(v - v_trim) can produce, and the loop deadlocks.  That is a finding
// about the model-to-geometry mapping, not a licence to change the rule.
//
// FINDING F-LEDGER, and the one semantic change the corner-online tables
// force here.  The rule above is the RECORD's, and its gate is the seam test
// ALONE — the family only picks the cell.  On the record's anchored mission
// that is harmless: measured, it and run_mission_interleaved()'s one-sided
// trim-only rule produce the IDENTICAL mission (937 ticks, arrival 925,
// offset +9, 9 shifts, at q = 1, 2 and 3), because every shift that fires is
// a retard on the trim family during the engage acceleration.
//
// On the corner-online mission it is NOT harmless.  The corner windows carry
// interior ticks whose seam angle is zero (the seams sit at in-window stages
// 7 and 16 only), so the seam-only gate lets the ledger re-time INSIDE a
// window, against a reference whose pace there is not the trim cell it
// credits.  Measured with the C core on the corner tables at q = 2: 473
// re-timings, offset -377, the mission "arriving" at tick 546 with the model
// state diverged (in-window |e_c| 5327 m against the 2.0 m band).  It does
// not fly.
//
// The fix is to gate on the CURRENT family, which is what the certified
// interleaved rule does (run_mission_interleaved in
// model/sim_corner_online.py: `p_phase < -HYST and ang(tau) == 0 and
// ang(tau+1) == 0 and fam(tau) == "t"`).  Written symbolically through
// BS_REF_FAM_TRIM, so the same source is correct for either table pair —
// the constant is 2 in the record header and 0 in the corner header.
//
// The two-sided form is KEPT.  Measured against the certified one-sided rule
// on the corner tables at q = 2, the two are identical in every invariant
// (983 ticks, arrival tick 971, offset +55, 55 re-timings, identical tau
// axis): with the family gate in place the advance branch never fires, so
// the extra arm is unreachable rather than merely unused.  CELL_MIN likewise
// becomes unreachable and is kept only so the rule still reads as the
// record's.
// TRIGGER QUANTITY (2026-08-22).  The record triggers on *phase, its raw
// integral of the tangential pace deviation.  The plant makes that
// e_l + sum Ts*nu + credits, so on the corner-online families it is
// dominated by nu -- the optimizer's own free clock input, which absorbs
// the segment pace schedule -- and is nearly uncorrelated with any lag
// (measured on the flown trace: corr(phase, e_l) = 0.11, corr(phase,
// vehicle along-track) = 0.30).  Consequence, measured: all 55 retards
// fired while the vehicle was AHEAD (median +3.5 m), each stalling the
// published position for one tick at 10.7 m/s while the feed-forward
// still commanded 11.2 m/s -- an internally inconsistent PVA triple that
// injected +2.4 m of lead per event and owned 95.7 % of the squared
// tracking error.  `trig` is therefore the model's own longitudinal
// deviation e_l = xi[2], the axis the corridor faces are already written
// on.  Host measurement at q=2 (face 1.60): 986 ticks / 58 retards /
// 243.50 s / realized violation 0.0646  ->  928 / 0 / 229.00 / 0.0000.
#if BS_STAGE2_WIN
// a slow tick of the pace family past the engage ramp: the window band
// and the absolute disc apply there, and the re-timing ledger does not
// (F-LEDGER's family gate, in the runtime's per-tick form)
static bool slow_at(int32_t tau)
{
    if (mt == nullptr) return false;
    if (tau < mt->sched.pace_t0) return false;
    if (fam_at(tau) != BS_MB_FAM_T) return false;
    return pace_at(tau) < mt->sched.pace_vslow;
}
#endif

static bool ledger_step(double *phase, int32_t *offset, int32_t tau,
                        double trig)
{
    if (mt == nullptr) return false;
    if (!(fabs(trig) > mt->hyst)) return false;
#if BS_STAGE2_WIN
    if (slow_at(tau)) return false;                          // Stage 2 (c)
#endif
    // exact zero test (see plant_step): the seam table holds exact values
    if (fabs(ang_at(tau)) > 0.0 || fabs(ang_at(tau + 1)) > 0.0) return false;
    if (fam_at(tau) != BS_MB_FAM_T) return false;          // F-LEDGER
    const double cell = (fam_at(tau) == BS_MB_FAM_T) ? mt->cell_t
                                                     : mt->cell_min;
    if (trig < 0.0) { *offset += 1; *phase += cell; }
    else            { *offset -= 1; *phase -= cell; }
    return true;
}

// Read #1 of the ingress engage: shape the leg from the measured position.
// Returns false when the geometry cannot be served (farther than the
// ING_MAX-cell cap, or an unexpected face layout); the READY gate then keeps
// waiting.  There is NO angle refusal and NO minimum length: the
// slow-junction / existing-family probes certify (at the state level) that
// every junction geometry is admissible when the crossing is slow, and the
// horizon schedules the crossing speed itself.
static bool ing_build(const Vector3p &p)
{
    if (mt == nullptr) return false;
    const double dn = mt->path[0] - (double)p.x;
    const double de = mt->path[1] - (double)p.y;
    const double L = sqrt(dn * dn + de * de);
    if (!(L > 0.1)) return false;
    // Tick count from the published ramp: rest to trim under the
    // feedback-preserving forward-face share (0.8 * 3.119) with the
    // tangential jerk face, capped at trim -- the D3/D4 profile, applied at
    // the PUBLISH layer.  n = first tick whose cumulative arc covers L.
    static const double ING_A = 0.8 * 3.119;     // m/s^2
    static const double ING_J = 3.5355339;       // m/s^3, jerk face
    ing_dir_n = dn / L;
    ing_dir_e = de / L;
    // junction angle: the heading change from the leg chart onto row 1's
    // chart, in the seam table's own convention (applied by plant_step /
    // T_switch as R(-chi) on the state pairs).
    {
        const double psi_leg = atan2(ing_dir_e, ing_dir_n);
        double chi = degrees(mt->psi[0] - psi_leg);
        while (chi > 180.0) chi -= 360.0;
        while (chi < -180.0) chi += 360.0;
        ing_ang_deg = chi;
    }
    // Junction-speed cap: the position controller rounds a published kink
    // by ~K * v^2 * tan(chi/2) (fitted on measured SITL misses: 1.65 m at
    // 11.1 m/s / 45 deg and 2.53 m at 8.7 m/s / 90 deg), so the ramp's
    // tail decelerates to the speed whose rounding stays inside the
    // corridor.  chi <= ~20 deg is uncapped; a reversal (chi -> 180) goes
    // to the floor.
    double v_j = mt->v_trim;
    {
        const double tn = tan(radians(fabs(ing_ang_deg)) * 0.5);
        if (tn > 1e-3) {
            v_j = sqrt(1.2 / (0.035 * tn));
        }
        // the mission's own leg-1 published arc starts from rest, so the
        // ingress junction must arrive slow regardless of chi
        if (v_j > 2.0) v_j = 2.0;
        if (v_j < 1.5) v_j = 1.5;
    }
    double arc = 0.0, vv = 0.0, aa = 0.0;
    int32_t n = 0;
    ing_arc[0] = 0.0;
    while (arc < L && n < ING_MAX) {
        // one Ts step of the jerk-shaped ramp (accelerate, cap, coast),
        // bounded by backward reachability of the junction speed under the
        // same face share: v <= sqrt(v_j^2 + 2 A (L - arc)).
        double aa_n = aa + ING_J * BS_TS;
        if (aa_n > ING_A) aa_n = ING_A;
        double vv_n = vv + 0.5 * (aa + aa_n) * BS_TS;
        if (vv_n > mt->v_trim) { vv_n = mt->v_trim; aa_n = 0.0; }
        const double rem = (L - arc > 0.0) ? (L - arc) : 0.0;
        const double v_brk = sqrt(v_j * v_j + 2.0 * ING_A * rem);
        if (vv_n > v_brk) { vv_n = v_brk; aa_n = 0.0; }
        arc += 0.5 * (vv + vv_n) * BS_TS;
        vv = vv_n;
        aa = aa_n;
        n++;
        ing_arc[n] = arc;
    }
    if (arc < L) return false;                   // beyond the ramp's reach
    // normalize so the last sample lands exactly on the path start: the
    // residual (< one cell) is scaled across the profile rather than
    // stepped at the junction.
    const double scale = L / arc;
    for (int32_t k = 0; k <= n; ++k) {
        ing_arc[k] *= scale;
    }
    // anchor = the ramp's rest point: exactly L before the path start, so
    // p_sched(-n) is the hover point and p_sched(0) the junction vertex.
    ing_anchor_n = mt->path[0] - L * ing_dir_n;
    ing_anchor_e = mt->path[1] - L * ing_dir_e;

    // family "tc": the row family with the corridor pair opened to 5 m.
    // The e_c rows are identified structurally (h = +-e3, no input part)
    // rather than by index, and exactly two must exist.
    memcpy(ing_rows, &mt->rows[(size_t)BS_MB_FAM_T * BS_NROW * 10],
           sizeof(ing_rows));
    int patched = 0;
    for (int i = 0; i < BS_NROW; ++i) {
        double *row = &ing_rows[(size_t)i * 10];
        bool ec_only = fabs(row[3]) > 0.0;
        for (int j = 0; j < 9 && ec_only; ++j) {
            if (j != 3 && fabs(row[j]) > 0.0) {
                ec_only = false;
            }
        }
        if (ec_only) {
            row[9] = ING_COR_WIDE;
            patched++;
        }
    }
    if (patched != 2) return false;
    // ang_at(0)/off_at(0) override the flash values at the junction tick;
    // that is sound only while the mission's own tick 0 is seam-free.
    // Refuse rather than silently rewrite if a future table breaks this.
    if (fabs(mt->ang[0]) > 0.0) return false;

    // junction seam: T_switch(chi) in bs_rot's own layout (pairs (0,4),
    // (1,5), (2,3): y_a' = c y_a + s y_b, y_b' = -s y_a + c y_b) and the
    // affine kick off_vec(V_TRIM, V_TRIM, chi).
    const double cr = cos(radians(ing_ang_deg));
    const double sr = sin(radians(ing_ang_deg));
    memset(ing_rotmat, 0, sizeof(ing_rotmat));
    const int pa[3] = { 0, 1, 2 }, pb[3] = { 4, 5, 3 };
    for (int k = 0; k < 3; ++k) {
        ing_rotmat[pa[k] * BS_NX + pa[k]] = cr;
        ing_rotmat[pa[k] * BS_NX + pb[k]] = sr;
        ing_rotmat[pb[k] * BS_NX + pa[k]] = -sr;
        ing_rotmat[pb[k] * BS_NX + pb[k]] = cr;
    }
    // the ingress crossing is capped at ~2 m/s (leg-1 starts from rest),
    // so the seam kick carries the crossing pace, not trim
    memset(ing_offvec, 0, sizeof(ing_offvec));
    ing_offvec[0] = v_j * (cr - 1.0);
    ing_offvec[4] = -v_j * sr;

    // the RAM schedule: n leg cells, the junction at destination tick n,
    // then the first BS_N body entries so the horizon spans the seam.  For
    // tau >= 0 the driver switches to bs_sched_mission (ensure_problem), so
    // entries beyond index n + BS_N are never read.
    for (int32_t i = 0; i < n; ++i) {
        ing_fam[i] = -1;                    // rows_extra: the wide leg family
        ing_rot[i] = 0;                     // identity
        ing_off[i] = 0;                     // zero vector
    }
    ing_fam[n] = mt->sfam[0];
    ing_rot[n] = -1;                        // rot_extra: T_switch(chi)
    ing_off[n] = -1;                        // off_extra: the junction kick
    for (int32_t i = 1; i <= BS_N; ++i) {
        ing_fam[n + i] = mt->sfam[i];
        ing_rot[n + i] = mt->srot[i];
        ing_off[n + i] = mt->soff[i];
    }
    bs_sched_ingress.family = ing_fam;
    bs_sched_ingress.rotation = ing_rot;
    bs_sched_ingress.offset = ing_off;
    bs_sched_ingress.length = (int)(n + BS_N + 1);
    bs_sched_ingress.periodic = 0;
    bs_sched_ingress.rows_extra = ing_rows;
    bs_sched_ingress.rot_extra = ing_rotmat;
    bs_sched_ingress.off_extra = ing_offvec;
    bs_sched_ingress.rows_tab = mt->sched.rows_tab;
    bs_sched_ingress.P_tab = mt->sched.P_tab;
    bs_sched_ingress.K_tab = mt->sched.K_tab;
    bs_sched_ingress.rot_tab = mt->sched.rot_tab;
    bs_sched_ingress.off_tab = mt->sched.off_tab;
    ing_n = n;
    return true;
}

// Measured first-order lag of the horizontal-acceleration response,
// seconds.  0 disables the pre-compensation.  See seam_lag_attribution.md.
#ifndef BS_LAG_TAU
#define BS_LAG_TAU 0.305
#endif

// ------------------------------------------------------------------ driver
void AP_BSolver::init()
{
    if (_inited || _enable == 0) {
        return;
    }
#ifdef BS_STAGE2_VJ_SCALE
    bs_mb_vj_scale = BS_STAGE2_VJ_SCALE;   // Stage-2 experiment builds only
#endif
    bs_arena_words = bs_workspace_size();
    bs_arena = (double *)calloc(bs_arena_words, sizeof(double));
    if (bs_mb_block == nullptr) {
        // worst case over speed and geometry: the tick cap dominates
        // (bs_mission_size saturates at BS_MB_MAX_TICKS: 129442 B, gated
        // by tests/size_gate.c)
        const size_t worst = bs_mission_size(BS_MB_MAX_WP, 1.0e9, 3.0);
        const size_t margin = 16u * 1024u;
        // Measure what IS available -- the largest contiguous block, probed
        // up to worst + margin -- and reserve min(largest - margin, worst):
        // the margin stays free for the rest of boot, every mission whose
        // tables fit is served, and the builder refuses the rest with a
        // number.  Measured on the CubeOrange+: largest 127.0 kB after the
        // arena -> ~111 kB reserved (unchanged by this form to within the
        // 1 kB probe granularity); a 64-bit host reserves the worst case.
        //   `volatile`: the probe is a MEASUREMENT.  GCC deletes a
        //   calloc/free pair whose result is only null-tested and folds the
        //   test as "succeeded", which made this search a call-free loop
        //   that could only end by converging -- and with the old
        //   overflowed bound (> 2^63) `(lo + hi) / 2` wrapped and it never
        //   did: SITL's main thread spun here at 100 % (2026-09-05).
        //   `lo + (hi - lo) / 2` cannot overflow.
        size_t lo = 0, hi = worst + margin;
        {   // the top of the range first: if worst + margin allocates, the
            // search is moot and the reservation is exactly `worst` (the
            // bisection alone stops up to 1 kB short of hi and could never
            // reach it -- found in review)
            void *volatile t = calloc(1, hi);
            if (t) { free(t); lo = hi; }
        }
        while (hi - lo > 1024) {
            const size_t mid = lo + (hi - lo) / 2;
            void *volatile t = calloc(1, mid);
            if (t) { free(t); lo = mid; } else { hi = mid; }
        }
        bs_mb_largest = lo;
        bs_mb_block_len = (lo > margin) ? MIN(lo - margin, worst) : 0;
        bs_mb_block = bs_mb_block_len
                    ? calloc(1, bs_mb_block_len) : nullptr;
        if (bs_mb_block == nullptr) {
            bs_mb_block_len = 0;
        }
        bs_mb.set_reserved(bs_mb_block, bs_mb_block_len);
    }
    if (bs_arena == nullptr) {
        // Turn the failure into a measurement: MEMINFO reports total free
        // across regions, but a single arena needs one CONTIGUOUS block, and
        // that is the number that actually constrains this core.  Probe it.
        static size_t largest;
        static uint32_t last_ms;
        if (largest == 0) {
            size_t lo = 0, hi = 512u * 1024u;
            while (hi - lo > 1024) {
                const size_t mid = lo + (hi - lo) / 2;
                void *volatile t = calloc(1, mid);   // real probe, see above
                if (t) { free(t); lo = mid; } else { hi = mid; }
            }
            largest = lo ? lo : 1;
        }
        // RE-BROADCAST: a one-shot here is lost to USB re-enumeration, which
        // is the whole reason the bench repeats itself.
        const uint32_t now_ms = AP_HAL::millis();
        if (now_ms - last_ms >= 5000) {
            last_ms = now_ms;
            GCS_SEND_TEXT(MAV_SEVERITY_ERROR,
                          "BSLV: need %u B, largest block %u B",
                          (unsigned)(bs_arena_words * sizeof(double)),
                          (unsigned)largest);
        }
        return;
    }
    _inited = true;
    _active = false;
    _seeded = false;
    _finished = false;
    _tick = 0;
    _offset = 0;
    _phase_ledger = 0.0;
    _problem_phase = INT32_MIN;
    _problem_npin = -1;
    _prefix_moved = 0.0f;
    _overruns = 0;
    _render_n = 0;
    _render_us_max = _render_us_sum = _render_noop_us_max = 0;
    _bcast_ms = 0;
    _plan_valid = false;
    ing_n = 0;
    ing_climb_ok_ms = 0;
    memset(_xi, 0, sizeof(_xi));
    memset(_U, 0, sizeof(_U));
    memset(_committed, 0, sizeof(_committed));
    // the rest-engage IC is set at seed time, from the mission tables
    GCS_SEND_TEXT(MAV_SEVERITY_INFO, "BSLV: init (%u B arena)",
                  (unsigned)(bs_workspace_size() * sizeof(double)));
#if BS_BENCH_WCET
    if (!hal.scheduler->thread_create(
            FUNCTOR_BIND_MEMBER(&AP_BSolver::bench_thread, void),
            "BSB", 8192, AP_HAL::Scheduler::PRIORITY_IO, 0)) {
        GCS_SEND_TEXT(MAV_SEVERITY_ERROR, "BSB: thread create failed");
    }
#else
    // The solve gets its own thread; see the declaration for why a scheduler
    // slot cannot host it.
    if (!_thread_created) {
        if (!hal.scheduler->thread_create(
                FUNCTOR_BIND_MEMBER(&AP_BSolver::solver_thread, void),
                "BSLV", 8192, AP_HAL::Scheduler::PRIORITY_IO, 0)) {
            GCS_SEND_TEXT(MAV_SEVERITY_ERROR, "BSLV: thread create failed");
            return;
        }
        _thread_created = true;
    }
#endif
}

void AP_BSolver::engage()
{
    if (!_inited) init();
    if (!_inited) return;
    _active = true;
    _next_solve_ms = AP_HAL::millis();
    GCS_SEND_TEXT(MAV_SEVERITY_INFO, "BSLV: engaged (q=%d)", (int)_q);
}

void AP_BSolver::disengage()
{
    _active = false;
    WITH_SEMAPHORE(_plan_sem);
    _plan_valid = false;
}

void AP_BSolver::mission_poll(AP_Mission &mission, float v_cap_ms)
{
    if (_enable == 0) {
        return;
    }
    bs_mb_mission = &mission;
    bs_mb.poll(mission, v_cap_ms);
}

void AP_BSolver::update(bool ingress_ready)
{
    if (_enable == 0) {
        return;
    }
    if (!_inited) {
        init();
    }
    // PERIODIC RE-BROADCAST.  Boot statustexts are unrecoverable over USB:
    // the port re-enumerates in ~0.4 s and opens at ~0.7 s, by which time
    // "BSLV: init (N B arena)" is long gone.  The WCET bench solved this by
    // repeating itself every 5 s; the mission build needs the same, because
    // on Sim-on-Hardware the ONLY evidence that the arena was allocated and
    // the solver thread exists is that one boot line.  Repeating it also
    // carries the live tick/offset/overrun state, so a ground station that
    // attaches mid-mission still learns where the loop is.
    const uint32_t now_bcast = AP_HAL::millis();
    if (now_bcast - _bcast_ms >= 5000) {
        _bcast_ms = now_bcast;
        if (!_inited) {
            GCS_SEND_TEXT(MAV_SEVERITY_WARNING, "BSLV: NOT INITED (en %d)",
                          (int)_enable);
        } else if (_finished) {
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSLV: done %d ticks off %+d ovr %u",
                          (int)_tick, (int)_offset, (unsigned)_overruns);
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSLV: rnd n %u max %lu us mean %lu noop %lu",
                          (unsigned)_render_n, (unsigned long)_render_us_max,
                          (unsigned long)(_render_n ? _render_us_sum / _render_n : 0),
                          (unsigned long)_render_noop_us_max);
#if BS_STAGE2_REPACE
            GCS_SEND_TEXT(MAV_SEVERITY_INFO, "BSLV: repace n %u ref %u s %.3f end %ld",
                          (unsigned)bs_mb.repace_events(), (unsigned)bs_mb.repace_refused(),
                          (double)bs_mb.repace_scale(), (long)(mt ? mt->n_end : -1));
#endif
        } else if (_active) {
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSLV: t %d off %+d %.0f ms ovr %u",
                          (int)_tick, (int)_offset, (double)_solve_ms,
                          (unsigned)_overruns);
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSLV: rnd n %u max %lu us mean %lu noop %lu",
                          (unsigned)_render_n, (unsigned long)_render_us_max,
                          (unsigned long)(_render_n ? _render_us_sum / _render_n : 0),
                          (unsigned long)_render_noop_us_max);
#if BS_STAGE2_REPACE
            GCS_SEND_TEXT(MAV_SEVERITY_INFO, "BSLV: repace n %u ref %u s %.3f end %ld",
                          (unsigned)bs_mb.repace_events(), (unsigned)bs_mb.repace_refused(),
                          (double)bs_mb.repace_scale(), (long)(mt ? mt->n_end : -1));
#endif
        } else {
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSLV: armed-wait q%d ir%d arena %u B",
                          (int)_q, (int)_ir,
                          (unsigned)(bs_workspace_size() * sizeof(double)));
        }
    }
    // ENGAGEMENT MOVED TO THE OWNING MODE (AUTO's BSLV_MISSION submode
    // calls engage_try() while its run is armed).  This slot keeps only
    // init and the re-broadcast; the solve is on the solver thread.
    (void)ingress_ready;
}

// THE READY GATE, relocated from update(): the only place this class
// reads the AHRS.  Called by the owning mode (AUTO, staged toward the
// run's first waypoint) every loop while it wants engagement.  Two ways
// in: CLASSIC rest-engage within READY_RADIUS_M of the path start, or
// the INGRESS leg from any near-hover within the prefix cap — the
// approach then becomes the first leg of the anytime problem.  Returns
// true when the gate passed and any ingress leg is built; the CALLER
// then runs the engage contract in order: mpc_replay.start() (captures
// the pre-engage desired state) followed by engage().
bool AP_BSolver::engage_try()
{
    if (_enable == 0 || !_inited || _active || _finished ||
        mt == nullptr) {
        return false;
    }
    Vector3p p;
    Vector3f v;
    if (!(hal.util->get_soft_armed() &&
          AP::ahrs().get_relative_position_NED_origin(p) &&
          AP::ahrs().get_velocity_NED(v) && (-p.z) > 2.0)) {
        ing_climb_ok_ms = 0;
        return false;
    }
    const double dn = (double)p.x - mt->path[0];
    const double de = (double)p.y - mt->path[1];
    const float speed = v.xy().length();
    if (dn * dn + de * de <= (double)READY_RADIUS_M * READY_RADIUS_M &&
        speed <= READY_SPEED_MS) {
        ing_n = 0;                     // classic rest-engage at the start
        return true;                   // caller: replay.start(); engage()
    }
    if (speed <= READY_SPEED_MS) {
        const uint32_t now_ms = AP_HAL::millis();
        if (fabsf(v.z) < 0.5f) {
            if (ing_climb_ok_ms == 0) {
                ing_climb_ok_ms = now_ms;
            }
        } else {
            ing_climb_ok_ms = 0;
        }
        if (ing_climb_ok_ms != 0 && now_ms - ing_climb_ok_ms >= 1000) {
            if (ing_build(p)) {
                _problem_phase = INT32_MIN;    // RAM tables changed
                _problem_npin = -1;
                GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                              "BSLV: ingress %ld cells chi %+d deg",
                              (long)ing_n, (int)lround(ing_ang_deg));
                return true;           // caller: replay.start(); engage()
            }
            if (now_ms - ing_txt_ms >= 5000) {
                ing_txt_ms = now_ms;
                GCS_SEND_TEXT(MAV_SEVERITY_WARNING,
                              "BSLV: ingress too far for %d cells",
                              (int)ING_MAX);
            }
        }
    } else {
        ing_climb_ok_ms = 0;
    }
    return false;
}

// Reset the mission clock and model state for a fresh engage.  SAFE ONLY
// WHILE INACTIVE: the solver thread touches this state exclusively when
// _active (the idle branch only rebuilds mission tables), so a
// main-thread reset while disengaged is race-free.  This is what makes
// the sticky _finished re-armable per AUTO entry.
void AP_BSolver::reset_mission()
{
    if (_active) {
        return;
    }
    _seeded = false;
    _finished = false;
    _tick = 0;
    _offset = 0;
    _phase_ledger = 0.0;
    _problem_phase = INT32_MIN;
    _problem_npin = -1;
    _prefix_moved = 0.0f;
    _overruns = 0;
    ing_n = 0;
    ing_climb_ok_ms = 0;
    bs_mb.ring_reset();
    {
        WITH_SEMAPHORE(_plan_sem);
        _plan_valid = false;
    }
    memset(_xi, 0, sizeof(_xi));
    memset(_U, 0, sizeof(_U));
    memset(_committed, 0, sizeof(_committed));
    _tau_pub = INT32_MIN;
}

#if !BS_BENCH_WCET
void AP_BSolver::solver_thread()
{
    while (true) {
        hal.scheduler->delay(10);
        if (!_active) {
            // rebuild the mission tables while the solver is idle
            if (bs_mb_mission != nullptr) {
                const bool ran = bs_mb.build_pending(*bs_mb_mission);
                if (ran) {
                    mt = bs_mb.tables();
                }
            }
            continue;
        }
        if (!_seeded) {
            _anchor_ms = AP_HAL::millis();
            seed_once();
            _next_solve_ms = AP_HAL::millis();
            continue;
        }
        const uint32_t now = AP_HAL::millis();
        if (now < _next_solve_ms) {
            continue;
        }
        // Pace from the deadline, not from completion: if a solve overruns
        // the 4 Hz period the next one starts immediately rather than
        // compounding the lag.
        //
        // AND FROM THE DEADLINE MEANS THE OLD DEADLINE, NOT now: advancing
        // from `now` re-anchors every tick to the thread's wake-up time, so
        // each tick inherits the wake-up latency and the mission clock runs
        // slow by that latency FOREVER (measured in SITL: 254.8 ms mean
        // cadence, +5.2 s over the SN77 mission).  The deadline advances on
        // its own 250 ms grid; only a true overrun re-anchors.  With the N = 30 corner tables the pinned
        // solve FITS -- 135.7 ms at q = 2 on the Cube against 250 ms -- so
        // unlike the record's N = 40 build (469.5 ms at q = 1, B3) this path
        // should stay cold; _overruns says whether it did.
        _next_solve_ms += (uint32_t)(BS_TS * 1000.0);
        if (now >= _next_solve_ms) {
            _next_solve_ms = now;              // overrun: restart, don't chase
        }
        _anchor_ms = now;
        // WALL TIME OF THE TICK.  micros() around the whole of solve_once():
        // the problem init, the shift-append warm start, the pinned Newton
        // solve and the lift.  That is the bench's quantity plus O(BS_NV)
        // bookkeeping, so 135.7 ms (idle board, q = 2, N = 30) is directly
        // comparable.  Under Sim-on-Hardware the physics model and the full
        // vehicle stack share this processor, so what this measures is the
        // solve AND every preemption it suffers -- which is the point.
        const uint32_t tick_us0 = AP_HAL::micros();
        const bool ok = solve_once();
        const uint32_t tick_us = AP_HAL::micros() - tick_us0;
        if (ok) {
            _solve_ms = (float)tick_us * 1e-3f;
            if (tick_us > (uint32_t)(BS_TS * 1e6)) {
                // Not a failure on its own: the M-input committed prefix is
                // 2.5 s of standing plan and replay keeps publishing from it.
                // Counted so the log says how often, and when.
                if (_overruns < 65535) _overruns++;
            }
            write_log();
        }
    }
}
#endif

bool AP_BSolver::ensure_problem(int32_t tau_raw)
{
    // Schedule selection is the negative-tau index switch of the ingress
    // design: tau < 0 solves on the RUNTIME ingress schedule at phase
    // tau + ing_n; tau >= 0 solves on the flash mission schedule at phase
    // tau, exactly as before -- the certified body problem sequence is
    // untouched by construction.
    if (mt == nullptr) return false;
    const bs_schedule *sched = &mt->sched;
    int32_t phase, key;
    if (tau_raw < 0 && ing_n > 0) {
        phase = tau_raw + ing_n;
        if (phase < 0) phase = 0;
        sched = &bs_sched_ingress;
        key = phase - ing_n;         // negative: cannot collide with body keys
    } else {
        phase = (tau_raw < 0) ? 0 : tau_raw;
        if (phase > mt->n_end) phase = mt->n_end;
        key = phase;
    }
    // The engage seed is a FULL-horizon solve, so its problem must carry the
    // whole quadratic half; every maintained solve is pinned at npin = M and
    // reads only the rows at or above it.
    const int32_t npin = (_ir != 0 && _seeded) ? BS_M : 0;
    if (_problem_phase == key && _problem_npin == npin) {
        return true;
    }
    if (bs_arena == nullptr) {
        return false;
    }
    if (bs_problem_init_pinned(&bs_driver_problem, sched,
                               (int)phase, (int)npin,
                               bs_arena, bs_arena_words) != BS_OK) {
        // INT32_MIN, not -1: negative keys are legitimate during an ingress
        _problem_phase = INT32_MIN;
        return false;
    }
    _problem_phase = key;
    _problem_npin = npin;
    return true;
}

// ------------------------------------------------------------- engage seed
//
// The interleaved loop of record starts from an EXACT solve at tau = 0:
//   U, _ = MX(0).newton(np.zeros(3N), xi, maxit=300)     (tolerance 1e-7)
//   committed[t] = U[3t:3t+3] for t in 0..M-1
// There is no commitment to pin yet, so this ONE solve is full-horizon
// (npin = 0) and it commits the whole prefix at once; every solve after it is
// pinned at npin = M.  That is the entire transition rule.
//
// F-ENG-PRECONV IS OPEN and this is engineering under it: the certificates
// say nothing about the pre-convergence stretch, and capping the iteration
// (BSLV_QE) is a firmware safety choice, not a certified schedule.
void AP_BSolver::seed_once()
{
    // With an ingress leg the seed problem sits at tau = -ing_n on the
    // runtime schedule; classic keeps tau = 0 on the flash schedule.  The
    // NOMINAL IC is the same rest-engage state either way -- read #2 below
    // injects the measured deviations afterwards (the D6 order: shape,
    // seed against nominal, then set the actual state).
    if (!ensure_problem(-(int32_t)ing_n)) {
        return;
    }
    if (mt == nullptr) {
        return;
    }
    memset(_U, 0, sizeof(_U));
    memset(_xi, 0, sizeof(_xi));
    _xi[0] = -mt->v_trim;
    _tick = 0;
    _offset = 0;
    _phase_ledger = 0.0;

    const uint32_t t0 = AP_HAL::micros();
    int it = 0;
    const int cap = (int)_qe > 0 ? (int)_qe : 1;
#if BS_STAGE2_HB
    bs_hb_buckets(&bs_driver_problem, _U, _xi, mt->v_trim, BS_MB_FAM_T);
#endif
    static double grad[BS_NV];
    for (; it < cap; ++it) {
        double value = 0.0;
        if (bs_eval(&bs_driver_problem, _U, _xi, &value, grad, nullptr)
            != BS_OK) {
            break;
        }
        double gmax = 0.0;
        for (int i = 0; i < BS_NV; ++i) {
            const double a = fabs(grad[i]);
            if (a > gmax) gmax = a;
        }
        if (gmax <= 1e-7) break;
        double lam = 0.0;
        bs_newton_stats st = {};
        if (bs_newton(&bs_driver_problem, _U, _xi, 1, &lam, &st) != BS_OK) {
            break;
        }
        _lambda0 = (float)lam;
        if (st.armijo_fail || st.factor_fail || st.no_step) break;
        hal.scheduler->delay(1);        // yield: this loop runs for seconds
    }
    _solve_ms = (AP_HAL::micros() - t0) * 1e-3f;

    // READ #2 of the ingress engage.  The seed ran for seconds; re-read the
    // estimator and set the ACTUAL deviations on the leg chart.  The leg
    // geometry itself stays anchored to the flash path start -- drift lands
    // in the deviation state, not in a table translation, so the published
    // junction stays exact (equivalent to D6's translation for the cm-dm
    // drifts a hover produces, without a step at the junction).
    // Implausible reads REFUSE into the classic READY wait.
    if (ing_n > 0) {
        Vector3p p2;
        Vector3f v2;
        bool ok;
        {
            // this runs on the SOLVER thread: unlike the READY gate's reads
            // (scheduler context), the AHRS must be locked here
            WITH_SEMAPHORE(AP::ahrs().get_semaphore());
            ok = AP::ahrs().get_relative_position_NED_origin(p2) &&
                 AP::ahrs().get_velocity_NED(v2);
        }
        if (ok) {
            const double rn = (double)p2.x - ing_anchor_n;
            const double re = (double)p2.y - ing_anchor_e;
            const double el = rn * ing_dir_n + re * ing_dir_e;
            const double ec = -rn * ing_dir_e + re * ing_dir_n;
            const double vt = (double)v2.x * ing_dir_n +
                              (double)v2.y * ing_dir_e;
            const double vc = -(double)v2.x * ing_dir_e +
                              (double)v2.y * ing_dir_n;
            if (fabs(el) <= 2.0 && fabs(ec) <= 2.0 &&
                fabs(vt) <= 2.0 && fabs(vc) <= 2.0) {
                _xi[0] = -mt->v_trim + vt;
                _xi[2] = el;
                _xi[3] = ec;
                _xi[4] = vc;
            } else {
                ok = false;
            }
        }
        if (!ok) {
            GCS_SEND_TEXT(MAV_SEVERITY_WARNING,
                          "BSLV: ingress read refused, classic READY");
            // ORDER: everything else first, _active LAST — the gate only
            // runs while !_active, so this closes the window in which a
            // preempting gate could rebuild the leg mid-reset.
            _seeded = false;
            ing_n = 0;
            ing_climb_ok_ms = 0;
            _problem_phase = INT32_MIN;
            _problem_npin = -1;
            _active = false;
            return;
        }
    }

    memcpy(_committed, _U, sizeof(double) * 3 * BS_M);
    _seeded = true;
    GCS_SEND_TEXT(MAV_SEVERITY_INFO, "BSLV: seed %d steps, %.0f ms", it,
                  (double)_solve_ms);
}

// ------------------------------------------------------------ maintained tick
//
// One tick of run_mission_interleaved(), with the environment (plant step and
// rigid clock) and the controller (pinned solve) in the order the reference
// runs them.  Nothing here reads the vehicle.
bool AP_BSolver::solve_once()
{
    // The ingress prefix shifts the clock start: at engage, tau = -ing_n.
    // No re-timing can fire on the leg (fam_at returns the non-trim
    // sentinel there), so tau is non-decreasing and crosses 0 exactly once,
    // at the junction.
    if (mt == nullptr) {
        return false;
    }
    const int32_t tau_raw = _tick - _offset - ing_n;
    if (tau_raw >= mt->n_end) {
        if (!_finished) {
            _finished = true;
            _active = false;
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSLV: mission complete, %d ticks, off %+d",
                          (int)_tick, (int)_offset);
        }
        return false;
    }
    const int32_t tau = (ing_n > 0) ? tau_raw
                       : ((tau_raw < 0) ? 0 : tau_raw);

    {   // render slice, timed: the one per-tick cost the bench ladder
        // does not cover (leg renders happen ~n_v times per mission)
        const int32_t r_before = bs_mb.rendered_to();
        const uint32_t r0 = AP_HAL::micros();
        const bool r_ok = bs_mb.render_to(tau + BS_N + 2);
        const uint32_t r_us = AP_HAL::micros() - r0;
        if (bs_mb.rendered_to() != r_before) {
            _render_n++;
            _render_us_sum += r_us;
            if (r_us > _render_us_max) _render_us_max = r_us;
        } else if (r_us > _render_noop_us_max) {
            _render_noop_us_max = r_us;
        }
        if (!r_ok) {
            GCS_SEND_TEXT(MAV_SEVERITY_WARNING, "BSLV: ring starved at %ld",
                          (long)tau);
            _active = false;
            return false;
        }
    }

    if (!ensure_problem(tau)) {
        return false;
    }

    // (a) REPLAY the committed prefix.  The pinned solve provably cannot have
    //     moved it (causality: no stage t < M reaches a column >= 3M), so
    //     this is an identity — measure it rather than trust it.
    double moved = 0.0;
    for (int i = 0; i < 3 * BS_M; ++i) {
        const double dev = fabs(_U[i] - _committed[i]);
        if (dev > moved) moved = dev;
        _U[i] = _committed[i];
    }
    if ((float)moved > _prefix_moved) _prefix_moved = (float)moved;

    _tau_pub = tau_raw;
#if BS_STAGE2_REPACE
    // Stage 2 (c): feed the re-pace verdict on slow ticks -- the model's
    // lead, its clock input and the worst face violation of the applied
    // input (chart-frame margins as the host gate reads them; seam ticks
    // excluded, the frame rotation re-expresses the polygon there).
    if (ing_n == 0 && slow_at(tau)) {
        double viol = 0.0;
        if (!(fabs(ang_at(tau)) > 0.0)) {
            const double *rows = &mt->rows[(size_t)BS_MB_FAM_T * BS_NROW * 10];
            const double pc = pace_at(tau);
            for (int i = 0; i < BS_NROW; ++i) {
                const double *r = &rows[(size_t)i * 10];
                double g = 0.0;
                for (int j = 0; j < 6; ++j) g += r[j] * _xi[j];
                for (int j = 0; j < 3; ++j) g += r[6 + j] * _U[j];
                const double v = g - bs_mb_win_margin(&mt->sched, BS_MB_FAM_T,
                                                      tau, pc, i, r);
                if (v > viol) viol = v;
            }
        }
        bs_mb.repace_feed(_xi[2], _U[1], viol);
    }
#endif

    // (b) PUBLISH the standing commitment.  Node i of this plan is node i+1
    //     of the previous one, exactly: same inputs, same recursion, one tick
    //     on.  Nothing is re-anchored at a measurement.
    {
        WITH_SEMAPHORE(_plan_sem);
        lift_prefix(_plan);
        _plan_valid = true;
    }

    // (c) warm start: shift-1 + terminal-law append, at the OLD phase and the
    //     PRE-plant state (the reference's `KK[fam(tau+N)] @ xiN`).
    static double U_shift[BS_NV];
    bs_shift_append(&bs_driver_problem, _U, _xi, U_shift);

    // (d) plant step and the rigid clock, on the MODEL state
    _phase_ledger += BS_TS * _xi[0];
    const double trig_el = _xi[2];        // e_l at the ledger instant
    plant_step(_xi, _U, ang_at(tau + 1), off_at(tau + 1));
    ledger_step(&_phase_ledger, &_offset, tau, trig_el);

    // (e) the PINNED solve at the next phase.  npin = M: the committed prefix
    //     is data, only the free tail is optimized.  This is the cadence the
    //     interleaved certificate (wp5_h1_interleaved_q.md, uniform q >= 1 at
    //     varrho = 4.220230e-2) is stated over, and the object premise A4
    //     describes.  Re-committing every input every tick — which is what
    //     calling bs_newton here would do — is a different controller.
    const int32_t tau_next = _tick + 1 - _offset - ing_n;
    if (!ensure_problem(tau_next)) {
        return false;
    }
#if BS_STAGE2_HB
    bs_hb_buckets(&bs_driver_problem, U_shift, _xi, mt->v_trim, BS_MB_FAM_T);
#endif
    double lambda0 = 0.0;
    bs_newton_stats stats = {};
    if (bs_newton_pinned(&bs_driver_problem, U_shift, _xi, (int)_q,
                         BS_M, &lambda0, &stats) != BS_OK) {
        return false;
    }
    _lambda0 = (float)lambda0;
    _bt = (uint8_t)MIN(stats.backtracks, 255);
    _armijo_fail = (uint8_t)MIN(stats.armijo_fail, 255);
    _no_step = (uint8_t)MIN(stats.no_step, 255);
    memcpy(_U, U_shift, sizeof(_U));

    // (f) exactly ONE input becomes data per tick: stage M-1 of the new plan,
    //     which is the input that will be flown M ticks from now.  The window
    //     slides by one.
    memmove(_committed, _committed + 3, sizeof(double) * 3 * (BS_M - 1));
    memcpy(_committed + 3 * (BS_M - 1), _U + 3 * (BS_M - 1),
           sizeof(double) * 3);

    _tick++;
    return true;
}

// -------------------------------------------------------------------- lift
//
// Quotient -> plane: THE CHART-FRAME RECONSTRUCTION, which is the map the
// quotient is defined by and the one the successor simulator publishes
// (model/sim_corner_online.py :: window_track, used by published_track and
// by the affine-identity gate):
//
//     p(tau) = p_sched(tau) + e_l * t_hat(tau) + e_c * n_hat(tau)
//
// with e_l = xi[2] positive AHEAD of the reference sample, e_c = xi[3]
// positive LEFT of the chart heading, and (t_hat, n_hat) the CHART frame --
// bs_ref_psi, which turns only at the frame seams -- not the local tangent of
// the polyline.  Host gate L1 replays the flown q = 2 mission through exactly
// this arithmetic and reproduces the simulator's published_track over all 19
// windows to 5.9e-14 m.
//
// THREE POINTS, EACH SETTLED BY MEASUREMENT RATHER THAN TASTE.
//
// (1) p_phase DOES NOT APPEAR.  Its predecessor used s_veh = SCHED[tau] +
//     p_phase as the along-track position and omitted e_l entirely, on the
//     argument that p_phase gains a cell at a re-timing and so keeps the
//     track continuous.  That is a double book-keeping of the along-track
//     axis: p_phase is the LEDGER accumulator (a raw integral of Ts*delta
//     plus the cells it has spent) and reached 9.40 m on the flown mission,
//     while the model's actual longitudinal deviation e_l stays inside
//     2.77 m.  On a straight row the two lifts differ by exactly
//     (p_phase - e_l) * t_hat -- verified to 4.6e-13 m over all 482 row
//     ticks -- i.e. by the double-booked term and nothing else.
//
//     There is also nothing to add back for the re-timing.  tau = k - o
//     ALREADY carries every re-timing the ledger has performed, so a retard
//     is expressed by tau repeating for one tick, and p_phase's only job is
//     to decide WHEN that happens.  It is clock state, not geometry.
//
//     The consequence is real and is not hidden: when tau repeats, the
//     published sample repeats, so the track advances by |delta e_l| (median
//     0.23 m) instead of a cell (median 2.76 m) on each of the 55 re-timing
//     ticks.  That is the reference waiting for a vehicle that is by
//     construction >= 2.1 m behind, and it is the correct position command;
//     what it costs is a one-tick notch in the DIFFERENTIATED feed-forward
//     below.  Measured on the published track: every tick whose lean demand
//     exceeds the 35 deg authority is adjacent to one of those 55 holds
//     (peak 67.77 deg); with them excluded the peak is 31.04 deg and the
//     in-window p95 is 30.80 deg.
//
// (2) The velocity and acceleration handed over are the DERIVATIVES OF THIS
//     POSITION TRACK (central differences over the node set), not the model's
//     own (pace + delta, a).  The position track is the object of record, so
//     the feed-forward has to be its own derivative or the position
//     controller is fed a contradiction.  With the chart reconstruction the
//     differentiator has a far easier job than it did with the arc lift: the
//     in-window peak lean demand falls from 70.86 deg to 38.51 deg and the
//     ticks over the authority from 194 to 100 -- all 100 of them now the
//     clock holds of point (1) rather than the corner geometry.
//
// (3) On the straight parts this reduces to the predecessor's form, as it
//     must: see the identity in (1).  It is only the along-track scalar that
//     changes there; through the corner windows the two are different objects
//     entirely, because the arc lift interpolated across the polyline kinks
//     and this one never leaves the published samples.
void AP_BSolver::lift_prefix(Plan &plan) const
{
    plan.n = PLAN_NODES;
    plan.dt_ms = (uint16_t)(BS_TS * 1000.0);
    plan.anchor_ms = _anchor_ms;
    plan.drag_k = BS_DRAG_K;

    // Roll the model state, the ledger and the clock forward through the
    // committed inputs, exactly as the maintained tick will.  `phase` is
    // carried because ledger_step needs it to decide when to re-time; it does
    // NOT enter the geometry -- see (1) above.
    double xi[6];
    memcpy(xi, _xi, sizeof(xi));
    double phase = _phase_ledger;
    int32_t off = _offset;
    int32_t tick = _tick;

    double pn[PLAN_NODES], pe[PLAN_NODES];
    int32_t tau_of[PLAN_NODES];
    for (uint8_t i = 0; i < PLAN_NODES; ++i) {
        int32_t tau = tick - off - ing_n;
        if (tau < -ing_n) tau = -ing_n;
        if (ing_n == 0 && tau < 0) tau = 0;
        if (mt != nullptr && tau > mt->n_end) tau = mt->n_end;
        tau_of[i] = tau;

        double px, py, tx, ty;
        ref_sample(tau, &px, &py);
        ref_frame(tau, &tx, &ty);
        pn[i] = px + xi[2] * tx + xi[3] * (-ty);
        pe[i] = py + xi[2] * ty + xi[3] * (tx);

        if (i + 1 >= PLAN_NODES) break;

        phase += BS_TS * xi[0];
        const double trig_el_i = xi[2];
        plant_step(xi, &_U[3 * i], ang_at(tau + 1), off_at(tau + 1));
        ledger_step(&phase, &off, tau, trig_el_i);
        tick += 1;
    }

    // ---------------------------------------------- feed-forward derivation
    //
    // THE POSITION TARGETS ABOVE ARE FINAL and are not touched here: they are
    // the certified published geometry and they sit inside the corridor at
    // every tick.  What follows derives only the velocity and acceleration
    // feed-forward from them, and it is HOLD-AWARE.
    //
    // WHY.  A re-timing holds the mission clock: tau repeats for one wall
    // tick, so two consecutive nodes carry the SAME reference sample and the
    // published track advances by |delta e_l| (median 0.23 m) instead of a
    // cell (median 2.76 m).  Differentiating that in WALL TICK k reads the
    // repeat as a near-stop and hands the position controller a brake at
    // exactly the moment the vehicle is, by the ledger's own trigger, at
    // least HYST = 2.106 m behind and should be closing the gap.  Measured on
    // the flown mission with the plain k-stencil: every published tick whose
    // lean demand exceeded the 35 deg authority was adjacent to one of the 55
    // holds (peak 67.77 deg against 31.04 deg everywhere else), the flown
    // track left the 2 m band by 0.61 m through the corner windows, and 8.2 %
    // of samples rode the authority.
    //
    // THE STENCIL.  Differentiate in MISSION CLOCK tau rather than in wall
    // tick k.  A hold contributes zero tau-advance, so the repeated node is
    // removed from the stencil rather than being allowed to contribute a
    // zero-velocity sample; the derivative at that tau is then taken from its
    // tau-adjacent neighbours.  Concretely: compress the node list to its
    // DISTINCT mission ticks (keeping the first node of each run, which is
    // the same representative sim_corner_online.published_track picks), apply
    // the existing central-difference rule -- unchanged, one-sided at the
    // ends, spacing Ts -- to that compressed sequence, and scatter the result
    // back so both nodes of a held pair publish the tau-rate.
    //
    // Behaviour at the four cases that exist on this mission:
    //   isolated hold        - one node dropped; both nodes of the pair get
    //                          the central difference over tau-1 .. tau+1.
    //   two holds one tick apart - both dropped; the stencil stays central
    //                          because it is built on the compressed axis.
    //   the engage cluster   - 9 holds inside the first 11 ticks, so L can
    //                          fall to about 6 of 11 nodes; the same rule
    //                          applies and only the one-sided ends move.
    //   the mission end      - tau clamps at BS_REF_NEND, so the trailing
    //                          nodes collapse to one entry.  That is correct
    //                          rather than special: the reference really is
    //                          stationary at the hover point, and the last
    //                          distinct pair still supplies the deceleration.
    //   L == 1               - no distinct pair at all; publish zero rather
    //                          than divide by a zero span.
    //
    // ON A HOLD-FREE PLAN THIS IS THE OLD CODE.  Every tau is distinct, so
    // L == PLAN_NODES, keep[m] == m and rep[i] == i, and the arithmetic below
    // is performed on the same values in the same order.  Gated numerically
    // on the host: max feed-forward difference on hold-free ticks is exactly
    // 0.
    uint8_t keep[PLAN_NODES], rep[PLAN_NODES], nkeep = 0;
    for (uint8_t i = 0; i < PLAN_NODES; ++i) {
        if (i == 0 || tau_of[i] != tau_of[i - 1]) {
            keep[nkeep++] = i;
        }
        rep[i] = (uint8_t)(nkeep - 1);
    }

    const double inv = 1.0 / BS_TS;
    double kn[PLAN_NODES], ke[PLAN_NODES];
    double kvn[PLAN_NODES], kve[PLAN_NODES];
    double kan[PLAN_NODES], kae[PLAN_NODES];
    for (uint8_t m = 0; m < nkeep; ++m) {
        kn[m] = pn[keep[m]];
        ke[m] = pe[keep[m]];
    }
    if (nkeep < 2) {
        kvn[0] = kve[0] = kan[0] = kae[0] = 0.0;
    } else {
        for (uint8_t m = 0; m < nkeep; ++m) {
            if (m == 0) {
                // second-order one-sided: the instantaneous velocity at the
                // node, not the mean over the first interval (which differs
                // by a Ts/2 * acc bias exactly where the replay reads).
                if (nkeep >= 3) {
                    kvn[m] = (-3.0 * kn[0] + 4.0 * kn[1] - kn[2]) * 0.5 * inv;
                    kve[m] = (-3.0 * ke[0] + 4.0 * ke[1] - ke[2]) * 0.5 * inv;
                } else {
                    kvn[m] = (kn[1] - kn[0]) * inv;
                    kve[m] = (ke[1] - ke[0]) * inv;
                }
            } else if (m == nkeep - 1) {
                kvn[m] = (kn[m] - kn[m - 1]) * inv;
                kve[m] = (ke[m] - ke[m - 1]) * inv;
            } else {
                kvn[m] = (kn[m + 1] - kn[m - 1]) * 0.5 * inv;
                kve[m] = (ke[m + 1] - ke[m - 1]) * 0.5 * inv;
            }
        }
        // ACCELERATION: the SECOND DIFFERENCE OF POSITION, not the first
        // difference of the published velocity.
        //
        // The previous form applied the velocity rule twice.  That is exact
        // in the interior, where both passes are central, but it is wrong at
        // the head of the plan, and the head is the only part the replay
        // interpolates between.  With kvn[0] one-sided and kvn[1] central,
        //
        //   kan[0] = (kvn[1] - kvn[0])/Ts
        //          = [ (p2 - p0)/2Ts - (p1 - p0)/Ts ] / Ts
        //          = (p2 - 2 p1 + p0) / (2 Ts^2),
        //
        // i.e. exactly HALF the second difference the same three points
        // carry, and kan[1] is three quarters of it.  Measured consequence
        // on the flown mission: the position controller was handed 0.36-0.50
        // of the published track's own curvature through every corner
        // window, and the missing 2.5-5 m/s^2 had to be produced by the
        // position-P / velocity-PI loop, which can only produce it by
        // standing off the plan.  That accounted for 42 % of the corridor
        // excursion; the remaining 58 % is the airframe's attitude lag
        // (measured tau = 0.305 s), which this stencil does not address.
        // See scripts/mpc/bsolver/doc/seam_lag_attribution.md.
        //
        // The interior stays a central second difference on the compressed
        // axis; the two ends take the second-order one-sided form so that
        // node 0 carries the curvature it actually has.  Positions are
        // untouched -- they are the certified published geometry.
        for (uint8_t m = 0; m < nkeep; ++m) {
            const double h2 = inv * inv;
            if (nkeep < 3) {
                kan[m] = kae[m] = 0.0;
            } else if (m == 0) {
                kan[m] = (kn[0] - 2.0 * kn[1] + kn[2]) * h2;
                kae[m] = (ke[0] - 2.0 * ke[1] + ke[2]) * h2;
            } else if (m == nkeep - 1) {
                kan[m] = (kn[m] - 2.0 * kn[m - 1] + kn[m - 2]) * h2;
                kae[m] = (ke[m] - 2.0 * ke[m - 1] + ke[m - 2]) * h2;
            } else {
                kan[m] = (kn[m + 1] - 2.0 * kn[m] + kn[m - 1]) * h2;
                kae[m] = (ke[m + 1] - 2.0 * ke[m] + ke[m - 1]) * h2;
            }
        }
    }
    // LAG PRE-COMPENSATION.  The airframe's horizontal acceleration follows
    // its demand through a first-order lag; the measured constant on this
    // platform is tau = 0.305 s (fitted against the position controller's own
    // acceleration target over all nineteen corners; agrees with the
    // configured attitude-loop constant 1/ATC_ANG_RLL_P = 0.222 s and the
    // interface-inclusive 0.32 s).  Publishing the plan's acceleration
    // therefore asks for a tilt the vehicle reaches one lag later, and the
    // position loop can only make up the difference by standing off the plan.
    //
    // Inverting the lag to first order, the demand whose ACHIEVED response is
    // the plan's acceleration is  a + tau * da/dt.  da/dt is taken from the
    // same compressed axis as the rest of the stencil.  This buys the lag
    // without adding the two actuator states an augmented model would need
    // (see papers/.../notes/attitude-lag-model.tex for that alternative).
    // Positions and velocities are untouched.
    double kjn[PLAN_NODES], kje[PLAN_NODES];
    for (uint8_t m = 0; m < nkeep; ++m) {
        if (nkeep < 2) {
            kjn[m] = kje[m] = 0.0;
        } else if (m == 0) {
            kjn[m] = (kan[1] - kan[0]) * inv;
            kje[m] = (kae[1] - kae[0]) * inv;
        } else if (m == nkeep - 1) {
            kjn[m] = (kan[m] - kan[m - 1]) * inv;
            kje[m] = (kae[m] - kae[m - 1]) * inv;
        } else {
            kjn[m] = (kan[m + 1] - kan[m - 1]) * 0.5 * inv;
            kje[m] = (kae[m + 1] - kae[m - 1]) * 0.5 * inv;
        }
    }
    for (uint8_t i = 0; i < PLAN_NODES; ++i) {
        const uint8_t m = rep[i];
        plan.pos_n[i] = (float)pn[i];
        plan.pos_e[i] = (float)pe[i];
        plan.vel_n[i] = (float)kvn[m];
        plan.vel_e[i] = (float)kve[m];
        plan.acc_n[i] = (float)(kan[m] + BS_LAG_TAU * kjn[m]);
        plan.acc_e[i] = (float)(kae[m] + BS_LAG_TAU * kje[m]);
    }
}

#if BS_BENCH_WCET
// Runs on its own thread: the solve costs hundreds of milliseconds on this
// target and must not be serviced from a scheduler slot.  (The arming refusal
// once blamed on that slot was in fact the scheduler task table being out of
// priority order -- see the note in AP_BSolver.h.)
//
// WHAT IS TIMED IS THE MAINTAINED CADENCE: bs_problem_init + shift-append +
// bs_newton_pinned at npin = M, which is what solve_once() runs 936 times out
// of 937.  The ENGAGE SEED (the one full-horizon solve) is timed separately
// as q = 0 in the report, because it happens once and off the deadline.
#if BS_BENCH_WCET
/* the demo mission (SN77 rows), embedded so the bench is SELF-HOSTING:
 * on a bench boot no ground station ever uploads a mission, so the
 * builder must be fed here or mt stays null and the ladder cannot run.
 * Building on-target is itself the new bench item: the mission compile
 * (including both DARE pairs) has only ever been timed on the host. */
static const double BVX[40] = {
    0.000, 0.000, 14.100, 14.100, 28.200, 28.200, 42.300, 42.300, 56.400, 56.400, 70.500, 70.500, 84.600, 84.600, 98.700, 98.700, 112.800, 112.800, 126.900, 126.900, 141.000, 141.000, 155.100, 155.100, 169.200, 169.200, 183.300, 183.300, 197.400, 197.400, 211.500, 211.500, 225.600, 225.600, 239.700, 239.700, 253.800, 253.800, 267.900, 267.900 };
static const double BVY[40] = {
    0.000, -99.561, -104.143, -3.684, 0.897, -99.561, -104.143, -3.684, 0.897, -99.561, -104.143, -3.684, 0.897, -99.561, -104.143, -3.684, 0.897, -99.561, -104.143, -3.684, 0.897, -99.561, -104.143, -3.684, 0.897, -99.561, -104.143, -3.684, 0.897, -99.561, -104.143, -3.684, 0.897, -99.561, -104.143, -3.684, 0.897, -99.561, -104.143, -3.684 };
#endif

void AP_BSolver::bench_thread()
{
    double *bench_mem = bs_arena;
    static double U[BS_NV];
    const size_t arena_words = bs_arena_words;
    static bs_problem problem;

    const double POISON = -1.2345678901234e-300;
    // [variant][q]: variant 0 = full setup, variant 1 = setup restricted to
    // the rows the pinned path reads.  BOTH are measured in one run, because
    // the tiling episode showed the host:target scaling can invert.
    struct { uint32_t mn, mx; uint64_t sum; uint16_t n; } r[2][4] = {};
    uint32_t high_water = 0;
    int status_ok = 1;

    /* ---- on-target mission build: the upload-time cost, measured ----
     * v_cap 11.8 is the deployment case; 6.0 is the DARE stress case
     * (rho(A_cl) rises as cruise falls; the fixed point runs longest
     * there).  The 11.8 tables are kept and published to mt for the
     * ladder, so the ladder runs on exactly what a real upload yields. */
    static bs_mission_tables bench_mt;
    static bs_mission_report bench_rep;
    uint32_t t_build[2] = {0, 0};
    int32_t n_build[2] = {0, 0};
    int8_t st_build[2] = {-1, -1};
    uint32_t fail_need = 0;
    int8_t fail_where = 0;          /* 1 calloc, 2 build, 3 rebuild, 4 ladder */
    int16_t fail_gate = 0;
    float fail_det = 0.0f;
    {
        double len = 0.0;
        for (int i = 0; i + 1 < 40; ++i) {
            const double dx = BVX[i + 1] - BVX[i];
            const double dy = BVY[i + 1] - BVY[i];
            len += sqrt(dx * dx + dy * dy);
        }
        const double caps[2] = {11.8, 6.0};
        for (int c = 0; c < 2 && status_ok; ++c) {
            /* hand the builder the whole reservation: its TAKE macro
             * bounds-checks the ACTUAL carve, and the loose sizing bound
             * (129 kB) overstates the real need (82-99 kB) by more than
             * the H7 has to spare */
            void *blk = bs_mb_block;
            const size_t need = bs_mb_block_len;
            if (blk == nullptr) {
                status_ok = 0; fail_where = 1;
                fail_need = (uint32_t)bs_mission_size(40, len, caps[c]);
                break;
            }
            memset(blk, 0, need);
            bs_mission_params pp;
            pp.v_cap_ms = caps[c];
            const uint32_t b0 = AP_HAL::micros();
            const bs_mb_status st = bs_mission_build(BVX, BVY, 40, &pp,
                                                     blk, need,
                                                     &bench_mt, &bench_rep);
            t_build[c] = AP_HAL::micros() - b0;
            n_build[c] = bench_rep.n_ticks;
            st_build[c] = (int8_t)st;
            if (st != BS_MB_OK) {
                status_ok = 0; fail_where = 2;
                fail_gate = (int16_t)bench_rep.gate;
                fail_det = (float)bench_rep.detail;
                break;
            }
        }
        /* rebuild at 11.8 LAST so mt points at the deployment tables
         * (the c-loop freed only the 6.0 block; bench_mt now holds the
         * 6.0 build, so build 11.8 again into a fresh block) */
        if (status_ok) {
            const size_t need = bs_mb_block_len;
            bs_mission_params pp;
            pp.v_cap_ms = 11.8;
            memset(bs_mb_block, 0, need);
            if (bs_mission_build(BVX, BVY, 40, &pp, bs_mb_block, need,
                                 &bench_mt, &bench_rep) != BS_MB_OK) {
                status_ok = 0;
                fail_where = 3;
                fail_need = (uint32_t)need;
                fail_gate = (int16_t)bench_rep.gate;
                fail_det = (float)bench_rep.detail;
            } else {
                mt = &bench_mt;
            }
        }
    }

    const uint8_t REPS = 25;
    for (uint8_t rep = 0; rep < REPS && status_ok; ++rep) {
        for (uint8_t v = 0; v < 2 && status_ok; ++v) {
            for (uint8_t q = 1; q <= 3; ++q) {
                // ENGAGE seed: the mission's own rest-engage IC at tau = 0
                // with a zero plan — the state the loop actually starts from,
                // and the one where the Newton iteration does real work.
                if (mt == nullptr) { status_ok = 0; break; }
                double xi[BS_NX] = { -mt->v_trim, 0, 0, 0, 0, 0 };
                for (uint16_t i = 0; i < BS_NV; ++i) U[i] = 0.0;
                if (rep == 0 && v == 0 && q == 1) {
                    for (uint32_t i = 0; i < arena_words; ++i)
                        bench_mem[i] = POISON;
                }

                const uint32_t t0 = AP_HAL::micros();
                if (bs_problem_init_pinned(&problem, &mt->sched, 0,
                                           v ? BS_M : 0, bench_mem,
                                           arena_words) != BS_OK) {
                    status_ok = 0; fail_where = 4; break;
                }
                double U_shift[BS_NV];
                bs_shift_append(&problem, U, xi, U_shift);
                double lam = 0.0;
                bs_newton_stats st;
                if (bs_newton_pinned(&problem, U_shift, xi, (int)q, BS_M,
                                     &lam, &st) != BS_OK) {
                    status_ok = 0; break;
                }
                const uint32_t dt = AP_HAL::micros() - t0;

                if (rep == 0 && v == 0 && q == 1) {
                    for (uint32_t i = arena_words; i > 0; --i) {
                        if (bench_mem[i - 1] != POISON) { high_water = i; break; }
                    }
                }
                if (r[v][q].n == 0 || dt < r[v][q].mn) r[v][q].mn = dt;
                if (dt > r[v][q].mx) r[v][q].mx = dt;
                r[v][q].sum += dt; r[v][q].n++;
                hal.scheduler->delay(1);
            }
        }
    }

    // Re-broadcast forever at 5 s: boot statustexts are unrecoverable over USB
    // (re-enumeration 0.4 s, port open 0.7 s — still too late to catch them).
    uint8_t line = 0;
    while (true) {
        if (!status_ok) {
            if (line & 1) {
                GCS_SEND_TEXT(MAV_SEVERITY_INFO, "BSB: FAILED");
            } else {
                GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                              "BSB fd w%d g%d d%.3g nd%lu av%lu lg%lu r%lu",
                              (int)fail_where, (int)fail_gate,
                              (double)fail_det,
                              (unsigned long)fail_need,
                              (unsigned long)hal.util->available_memory(),
                              (unsigned long)bs_mb_largest,
                              (unsigned long)bs_mb_block_len);
            }
        } else if (line == 0) {
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSB: PINNED npin=%u, engage seed, n=%u/q, f64",
                          (unsigned)BS_M, (unsigned)r[0][1].n);
        } else if (line <= 6) {
            const uint8_t v = (uint8_t)((line - 1) / 3);
            const uint8_t q = (uint8_t)(((line - 1) % 3) + 1);
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSB IR%u q=%u min %lu avg %lu max %lu us",
                          (unsigned)v, (unsigned)q,
                          (unsigned long)r[v][q].mn,
                          (unsigned long)(r[v][q].sum /
                                          (r[v][q].n ? r[v][q].n : 1)),
                          (unsigned long)r[v][q].mx);
        } else if (line == 7) {
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSB arena hw %lu B of %lu, bss %u B",
                          (unsigned long)(high_water * sizeof(double)),
                          (unsigned long)(bs_workspace_size() * sizeof(double)),
                          (unsigned)BS_STATIC_BSS);
        } else {
            const uint8_t c = (uint8_t)(line - 8);
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSB build v%u %lu us %ld ticks st%d",
                          (unsigned)(c == 0 ? 118 : 60),
                          (unsigned long)t_build[c],
                          (long)n_build[c], (int)st_build[c]);
        }
        line = (line + 1) % 10;
        hal.scheduler->delay(5000 / 8);
    }
}
#endif // BS_BENCH_WCET

#if HAL_LOGGING_ENABLED
// One record per accepted solve.  This is the evidence channel: without it
// the only sign the certified solver ran at all is a single engage statustext,
// which cannot distinguish "solving at rate" from "engaged once and silent".
// The fields are exactly the quantities the WP5 certificates bound -- the
// Newton decrement lambda_0, the iteration count q, the Armijo backtracks and
// the guarded exits -- plus the clock the mission gate checks (tick, offset,
// the ledger accumulator) and the MODEL deviations the plan drives to zero.
//
// El/Ec/Del are MODEL states, not tracking errors: nothing in this class
// measures the vehicle.  The tracking error is the vehicle against the
// published reference and is reconstructed from the flight log.
void AP_BSolver::write_log() const
{
    AP::logger().WriteStreaming(
        "BSLV",
        "TimeUS,Tick,Off,Q,Lam,Ms,El,Ec,Del,Ph,BT,AF,NS,Ovr,Rn,Re",
        "s-----mmnm----mm", "F-----0000----00", "QiiBffffffBBBHff",
        AP_HAL::micros64(),
        (int32_t)_tick, (int32_t)_offset, (uint8_t)_q,
        (float)_lambda0, (float)_solve_ms,
        (float)_xi[2], (float)_xi[3], (float)_xi[0], (float)_phase_ledger,
        _bt, _armijo_fail, _no_step, _overruns,
        _plan.pos_n[0], _plan.pos_e[0]);
}
#endif

uint16_t AP_BSolver::run_first_idx() const
{
    return bs_mb.run_first_idx();
}

uint16_t AP_BSolver::run_last_idx() const
{
    return bs_mb.run_last_idx();
}

int32_t AP_BSolver::tau_vertex(uint16_t wp_idx) const
{
    return bs_mb.tau_vertex(wp_idx);
}

bool AP_BSolver::mission_ready() const
{
    return mt != nullptr;
}

bool AP_BSolver::take_plan(Plan &out)
{
    WITH_SEMAPHORE(_plan_sem);
    if (!_plan_valid) {
        return false;
    }
    out = _plan;
    _plan_valid = false;
    return true;
}

#endif // AP_BSOLVER_ENABLED
