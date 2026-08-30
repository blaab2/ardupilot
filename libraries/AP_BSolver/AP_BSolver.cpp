#include "AP_BSolver.h"

#if AP_BSOLVER_ENABLED

#include <AP_AHRS/AP_AHRS.h>
#include <GCS_MAVLink/GCS.h>
#include <AP_Logger/AP_Logger.h>
#include <stdlib.h>
#include <string.h>

// The mission REFERENCE tables are a build-time choice, mirroring the one
// bsolver/core/bs_solver.h makes for the MODEL tables with BS_DATA_HEADER.
// The record pair (N = 40) and the corner-online pair (N = 30) define the
// same bs_ref_* symbols and cannot coexist in one translation unit, so this
// is a substitution, not an addition; wscript sets both flags together and
// AP_BSolver.cpp static_asserts that the two halves agree.  Default = the
// record, so the file stays buildable with no flags at all.
#ifndef BS_REF_HEADER
#define BS_REF_HEADER "bsolver/model/generated/bs_reference_data.h"
#endif

extern "C" {
#include "bsolver/core/bs_solver.h"
#include BS_REF_HEADER
}

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

// The ONE bs_problem the driver owns.  It is file-scope so ensure_problem()
// can keep it across ticks: the reference loop assembles at tau at the top of
// a tick and at (k+1) - o_off at the bottom, and those are the SAME phase one
// tick apart, so caching turns two assemblies per tick into one.  On target
// bs_problem_init is 16.8 % of a pinned cycle, so this is not cosmetic.
static bs_problem bs_driver_problem;

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
// ... and the MODEL half and the REFERENCE half must be the same pair.  They
// are selected by two independent flags, so a half-switched build is a real
// hazard, and it is silent: both halves define every symbol this file uses.
// Three checks, each catching a different mix.
//  (a) the clamped lookups here index bs_ref_ang / bs_ref_fam over a clock
//      the core's schedule must cover.  Record 930/976, corner 966/966.
static_assert(BS_REF_NCLK <= BS_N_MISSION,
              "the reference clock must not outrun the model's mission "
              "schedule: BS_DATA_HEADER and BS_REF_HEADER are mismatched");
//  (b) the family enum is renumbered between the pairs (record a=0,h=1,t=2;
//      corner t=0,c1=1,c2=2,c3=3,h=4).  BS_FAM_T exists only in the corner
//      model header, so this fires on corner-model + record-reference.
#ifdef BS_FAM_T
static_assert(BS_REF_FAM_TRIM == BS_FAM_T,
              "the reference's trim-family code disagrees with the model's: "
              "BS_DATA_HEADER and BS_REF_HEADER are mismatched");
#endif
//  (c) the corner model's clock arrays are the FULL padded mission and its
//      reference header matches them exactly; the record's do not.
#ifdef BS_NOFF
static_assert(BS_REF_NCLK == BS_N_MISSION,
              "the corner model header needs the corner reference header: "
              "BS_DATA_HEADER and BS_REF_HEADER are mismatched");
#endif

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
    // 2.7721) and the WCET fits: 135.7 ms on the Cube against the 250 ms
    // tick, 46 % margin.
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
    const int32_t i = (tau < 0) ? 0
                    : (tau > BS_REF_NPATH - 1 ? BS_REF_NPATH - 1 : tau);
    *x = bs_ref_path[2 * i];
    *y = bs_ref_path[2 * i + 1];
}

// CHART frame at schedule tick tau: the heading the quotient's deviations are
// expressed in.  It is NOT the finite-difference tangent of the polyline —
// the chart turns only at the frame seams, which is exactly why the seam
// angles exist as their own table.
static void ref_frame(int32_t tau, double *tx, double *ty)
{
    const int32_t i = (tau < 0) ? 0
                    : (tau > BS_REF_NPSI - 1 ? BS_REF_NPSI - 1 : tau);
    const double a = bs_ref_psi[i];
    *tx = cos(a);
    *ty = sin(a);
}

static double ang_at(int32_t tau)
{
    const int32_t i = (tau < 0) ? 0
                    : (tau > BS_REF_NCLK - 1 ? BS_REF_NCLK - 1 : tau);
    return bs_ref_ang[i];
}

static int fam_at(int32_t tau)
{
    const int32_t i = (tau < 0) ? 0
                    : (tau > BS_REF_NCLK - 1 ? BS_REF_NCLK - 1 : tau);
    return (int)bs_ref_fam[i];
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
#ifdef BS_NOFF
    const int32_t i = (tau < 0) ? 0
                    : (tau > BS_N_MISSION - 1 ? BS_N_MISSION - 1 : tau);
    return &bs_off[(size_t)bs_mission_off[i] * BS_NX];
#else
    static const double zero[BS_NX] = { 0.0 };
    (void)tau;
    return zero;
#endif
}

// One exact-ZOH plant step followed by the chart re-expression at the solve
// boundary:  xi <- T(tau+1) (A xi + B u).  This is A6/B6 of the certificates
// written out (verified entry for entry against the exported matrices), and
// T_switch's rotation of the three tangential/normal pairs (0,4), (1,5),
// (2,3) by the seam angle.
static void plant_step(double *xi, const double *u, double seam_deg,
                       const double *off)
{
    const double Ts = BS_REF_TS;
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
static bool ledger_step(double *phase, int32_t *offset, int32_t tau,
                        double trig)
{
    if (!(fabs(trig) > BS_REF_HYST)) return false;
    // exact zero test (see plant_step): the seam table holds exact values
    if (fabs(ang_at(tau)) > 0.0 || fabs(ang_at(tau + 1)) > 0.0) return false;
    if (fam_at(tau) != BS_REF_FAM_TRIM) return false;      // F-LEDGER
    const double cell = (fam_at(tau) == BS_REF_FAM_TRIM) ? BS_REF_CELLT
                                                         : BS_REF_CELLMIN;
    if (trig < 0.0) { *offset += 1; *phase += cell; }
    else            { *offset -= 1; *phase -= cell; }
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
    bs_arena_words = bs_workspace_size();
    bs_arena = (double *)calloc(bs_arena_words, sizeof(double));
    if (bs_arena == nullptr) {
        // Turn the failure into a measurement: MEMINFO reports total free
        // across regions, but a single arena needs one CONTIGUOUS block, and
        // that is the number that actually constrains this core.  Probe it.
        static size_t largest;
        static uint32_t last_ms;
        if (largest == 0) {
            size_t lo = 0, hi = 512u * 1024u;
            while (hi - lo > 1024) {
                const size_t mid = (lo + hi) / 2;
                void *t = calloc(1, mid);
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
    _problem_phase = -1;
    _problem_npin = -1;
    _prefix_moved = 0.0f;
    _overruns = 0;
    _bcast_ms = 0;
    _plan_valid = false;
    memset(_xi, 0, sizeof(_xi));
    memset(_U, 0, sizeof(_U));
    memset(_committed, 0, sizeof(_committed));
    _xi[0] = -BS_REF_VTRIM;          // rest-engage initial condition
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
        } else if (_active) {
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSLV: t %d off %+d %.0f ms ovr %u",
                          (int)_tick, (int)_offset, (double)_solve_ms,
                          (unsigned)_overruns);
        } else {
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSLV: armed-wait q%d ir%d arena %u B",
                          (int)_q, (int)_ir,
                          (unsigned)(bs_workspace_size() * sizeof(double)));
        }
    }
    // THE READY GATE, and the only place this class reads the AHRS.  The
    // mission clock starts at tau = 0 from the reference's rest-engage IC,
    // whose lifted node 0 is the path start; engaging anywhere else hands the
    // position controller a step it did not ask for.  So: armed, airborne,
    // the ingress actually open, and the vehicle station-keeping within
    // READY_RADIUS_M of the path start.
    //
    // The ingress term is not decoration.  Measured before it existed, the
    // solver self-engaged at 2 m during the GUIDED takeoff — while replay
    // still refuses plans — and burned 35 ticks at zero along-track speed,
    // slipping the assignment offset by one EVERY tick before tracking began.
    //
    // Deliberately NOT triggered over MAVLink: the point of this architecture
    // is that nothing in the control path crosses a link.
    if (!_active) {
        Vector3p p;
        Vector3f v;
        if (ingress_ready && hal.util->get_soft_armed() &&
            AP::ahrs().get_relative_position_NED_origin(p) &&
            AP::ahrs().get_velocity_NED(v) && (-p.z) > 2.0) {
            const double dn = (double)p.x - bs_ref_path[0];
            const double de = (double)p.y - bs_ref_path[1];
            const float speed = v.xy().length();
            if (dn * dn + de * de <= (double)READY_RADIUS_M * READY_RADIUS_M &&
                speed <= READY_SPEED_MS) {
                engage();
            }
        }
        return;
    }
    // Nothing else happens here.  The solve is on the solver thread; this
    // function runs in a scheduler slot and must stay cheap.
}

#if !BS_BENCH_WCET
void AP_BSolver::solver_thread()
{
    while (true) {
        hal.scheduler->delay(10);
        if (!_active) {
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
        // compounding the lag.  With the N = 30 corner tables the pinned
        // solve FITS -- 135.7 ms at q = 2 on the Cube against 250 ms -- so
        // unlike the record's N = 40 build (469.5 ms at q = 1, B3) this path
        // should stay cold; _overruns says whether it did.
        _next_solve_ms = now + (uint32_t)(BS_REF_TS * 1000.0);
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
            if (tick_us > (uint32_t)(BS_REF_TS * 1e6)) {
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

bool AP_BSolver::ensure_problem(int32_t phase)
{
    if (phase < 0) phase = 0;
    if (phase > BS_REF_NEND) phase = BS_REF_NEND;
    // The engage seed is a FULL-horizon solve, so its problem must carry the
    // whole quadratic half; every maintained solve is pinned at npin = M and
    // reads only the rows at or above it.
    const int32_t npin = (_ir != 0 && _seeded) ? BS_M : 0;
    if (_problem_phase == phase && _problem_npin == npin) {
        return true;
    }
    if (bs_arena == nullptr) {
        return false;
    }
    if (bs_problem_init_pinned(&bs_driver_problem, &bs_sched_mission,
                               (int)phase, (int)npin,
                               bs_arena, bs_arena_words) != BS_OK) {
        _problem_phase = -1;
        return false;
    }
    _problem_phase = phase;
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
    if (!ensure_problem(0)) {
        return;
    }
    memset(_U, 0, sizeof(_U));
    memset(_xi, 0, sizeof(_xi));
    _xi[0] = -BS_REF_VTRIM;
    _tick = 0;
    _offset = 0;
    _phase_ledger = 0.0;

    const uint32_t t0 = AP_HAL::micros();
    int it = 0;
    const int cap = (int)_qe > 0 ? (int)_qe : 1;
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
    const int32_t tau_raw = _tick - _offset;
    if (tau_raw >= BS_REF_NEND) {
        if (!_finished) {
            _finished = true;
            _active = false;
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSLV: mission complete, %d ticks, off %+d",
                          (int)_tick, (int)_offset);
        }
        return false;
    }
    const int32_t tau = (tau_raw < 0) ? 0 : tau_raw;

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
    _phase_ledger += BS_REF_TS * _xi[0];
    const double trig_el = _xi[2];        // e_l at the ledger instant
    plant_step(_xi, _U, ang_at(tau + 1), off_at(tau + 1));
    ledger_step(&_phase_ledger, &_offset, tau, trig_el);

    // (e) the PINNED solve at the next phase.  npin = M: the committed prefix
    //     is data, only the free tail is optimized.  This is the cadence the
    //     interleaved certificate (wp5_h1_interleaved_q.md, uniform q >= 1 at
    //     varrho = 4.220230e-2) is stated over, and the object premise A4
    //     describes.  Re-committing every input every tick — which is what
    //     calling bs_newton here would do — is a different controller.
    const int32_t tau_next = _tick + 1 - _offset;
    if (!ensure_problem(tau_next)) {
        return false;
    }
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
    plan.dt_ms = (uint16_t)(BS_REF_TS * 1000.0);
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
        int32_t tau = tick - off;
        if (tau < 0) tau = 0;
        if (tau > BS_REF_NEND) tau = BS_REF_NEND;
        tau_of[i] = tau;

        double px, py, tx, ty;
        ref_sample(tau, &px, &py);
        ref_frame(tau, &tx, &ty);
        pn[i] = px + xi[2] * tx + xi[3] * (-ty);
        pe[i] = py + xi[2] * ty + xi[3] * (tx);

        if (i + 1 >= PLAN_NODES) break;

        phase += BS_REF_TS * xi[0];
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

    const double inv = 1.0 / BS_REF_TS;
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

    const uint8_t REPS = 25;
    for (uint8_t rep = 0; rep < REPS && status_ok; ++rep) {
        for (uint8_t v = 0; v < 2 && status_ok; ++v) {
            for (uint8_t q = 1; q <= 3; ++q) {
                // ENGAGE seed: the mission's own rest-engage IC at tau = 0
                // with a zero plan — the state the loop actually starts from,
                // and the one where the Newton iteration does real work.
                double xi[BS_NX] = { -BS_REF_VTRIM, 0, 0, 0, 0, 0 };
                for (uint16_t i = 0; i < BS_NV; ++i) U[i] = 0.0;
                if (rep == 0 && v == 0 && q == 1) {
                    for (uint32_t i = 0; i < arena_words; ++i)
                        bench_mem[i] = POISON;
                }

                const uint32_t t0 = AP_HAL::micros();
                if (bs_problem_init_pinned(&problem, &bs_sched_mission, 0,
                                           v ? BS_M : 0, bench_mem,
                                           arena_words) != BS_OK) {
                    status_ok = 0; break;
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
            GCS_SEND_TEXT(MAV_SEVERITY_INFO, "BSB: FAILED");
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
        } else {
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSB arena hw %lu B of %lu, bss %u B",
                          (unsigned long)(high_water * sizeof(double)),
                          (unsigned long)(bs_workspace_size() * sizeof(double)),
                          (unsigned)BS_STATIC_BSS);
        }
        line = (line + 1) % 8;
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
