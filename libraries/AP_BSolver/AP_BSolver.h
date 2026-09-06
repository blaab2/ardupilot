/*
 * AP_BSolver — the certified anytime MPC, solved ONBOARD.
 *
 * Wraps the vendored bsolver C core (bsolver/, f64 — see
 * scripts/mpc/bsolver/sync_firmware.sh and doc/bsolver.tex) behind the same
 * plan-producing contract AP_MPCSolver uses, so MPCTrajReplay can commit its
 * output through the existing onboard ingress.
 *
 * WHY ONBOARD RATHER THAN STREAMED. The controller commits M inputs per solve
 * at 4 Hz, which is 2.5 s of committed input, and the inputs are jerks plus a
 * speed reference — so integrating the committed prefix through the exact-ZOH
 * plant yields a locally known, jerk-continuous position/velocity/acceleration
 * trajectory. Streaming setpoints instead puts the link and the position loop
 * inside the control path in BOTH directions; measured, that diverges (66 m
 * cross-track, the assignment clock re-timing 79 times in 200 ticks against
 * the mission gate's +3). There is also a correctness reason: the certified
 * premise A4 is the EXACT committed-prefix predictor, and a controller that
 * re-plans between the solver and the vehicle weakens the premise the proofs
 * rest on.
 *
 * Frames. The core's state is the 6-dim phase quotient
 *   xi = (delta, a, e_l, e_c, ec_dot, a_n)
 * with no absolute position: delta is speed deviation from trim, a the
 * along-track acceleration, e_l/e_c the along/cross-track errors, ec_dot the
 * cross-track rate, a_n the normal acceleration. The mission geometry that
 * turns that back into a plane — the composed path, its arc-length schedule,
 * the per-tick chart heading, and the ledger's seam angle / row family — is
 * generated data (the reference header BS_REF_HEADER selects — the flight
 * configuration is bs_reference_corner_data.h), not reconstructed here.  The
 * published plan is the CHART-FRAME reconstruction
 *   p(tau) = p_sched(tau) + e_l t_hat(tau) + e_c n_hat(tau),
 * the map the quotient is defined by; see lift_prefix() for why the ledger
 * accumulator p_phase is clock state and never appears in it.
 *
 * THE PLANNER IS OPEN-LOOP ON ITS OWN MODEL STATE, and that is the
 * architecture, not an omission.  Premise A4 is the EXACT committed-prefix
 * predictor: the certificates are stated over a loop whose state evolves by
 * the model recursion xi <- T(tau+1) (A xi + B u_0) under the committed
 * inputs and is NEVER re-measured (run_mission / run_mission_interleaved in
 * wp5_anytime_sim.py do exactly this).  So:
 *
 *   - the solver's state comes from the model, never from the AHRS;
 *   - the assignment clock runs the certificate's own rigid-clock ledger on
 *     that model state (accumulator, seam gate, per-family cell);
 *   - the committed prefix is REPLAYED, and the published plan is the
 *     standing commitment extended by one node per tick — never a fresh plan
 *     re-anchored at the measured position;
 *   - the AHRS is read ONLY by the supervision layer (the READY gate that
 *     decides when to start the mission clock) and downstream by the position
 *     controller that tracks the published reference.
 *
 * The predecessor of this file did the opposite: it projected the measured
 * NED position onto the path every tick and fed the measured gap into the
 * solver and into the clock.  Measured, that deadlocked — the ledger credited
 * a fixed trim cell per re-timing against a mission whose mean cell is
 * 2.4073 m, e_l ratcheted ~ +19.9 m per turn window against a 2.0 m barrier
 * row, the plan braked, and the vehicle crept until the harness cut.  The
 * deficit was structural: 140 ticks of advance the accumulator could never
 * produce.  Closing the loop at the SOLVER was the wrong repair; the loop
 * closes at the POSITION CONTROLLER, which is where a tracking error belongs.
 *
 * CONFIGURATION OF RECORD (corner-online, 2026-08-22).  The tables this
 * builds against are the N = 30 corner-online pair, selected in wscript
 * through BS_DATA_HEADER / BS_REF_HEADER.  Three things about them change
 * this file's assumptions and are handled explicitly below:
 *   - the seam map on the phase quotient is AFFINE,
 *     xi <- T(tau+1) (A xi + B u) + c(tau+1), so plant_step() carries the
 *     offset term.  76 of the 966 mission ticks have a non-zero c, and 38 of
 *     those have an IDENTITY rotation, so dropping it is not a small error;
 *   - the family enum is t=0, c1=1, c2=2, c3=3, h=4 (the record's was
 *     a=0, h=1, t=2), so every family test goes through BS_REF_FAM_TRIM;
 *   - the periodic 59-phase parity clock does not exist in this
 *     configuration (it was an "a"/"t" object).  Nothing here used it.
 *
 * SCOPE OF THIS RUNG, stated so it is not mistaken for more:
 *  - the solve runs on a dedicated thread at a 4 Hz DEADLINE.  At N = 30 the
 *    pinned solve FITS: 135.7 ms at q = 2 on the Cube against a 250 ms
 *    period.  (The record's N = 40 full solve was 469.5 ms at q = 1, 1.88x
 *    over budget — that measurement is what N = 30 was chosen against.)
 *  - the ENGAGE seed is a full-horizon solve iterated to the reference's own
 *    gradient tolerance (or BSLV_QE steps, whichever comes first).  That is
 *    engineering under the OPEN item F-ENG-PRECONV, not a certified staging.
 */
#pragma once

#include <AP_HAL/AP_HAL_Boards.h>

#include "bs_defs.h"

// The two cores cannot co-reside on the real target (614 kB against 507 kB
// free heap), so they are mutually exclusive at build time.
#if AP_BSOLVER_ENABLED && AP_MPCSOLVER_ENABLED
#error "AP_BSolver and AP_MPCSolver cannot both be enabled: they do not fit"
#endif

#if AP_BSOLVER_ENABLED

#include <AP_Common/AP_Common.h>
#include <AP_HAL/AP_HAL.h>
#include <AP_Math/AP_Math.h>
#include <AP_Param/AP_Param.h>

class AP_BSolver {
public:
    AP_BSolver();
    CLASS_NO_COPY(AP_BSolver);

    static const struct AP_Param::GroupInfo var_info[];

    // one committed prefix, lifted into the plane
    static const uint8_t PLAN_NODES = 11;      // M + 1, the committed prefix

    struct Plan {
        uint16_t n;
        uint16_t dt_ms;
        uint32_t anchor_ms;
        float drag_k;
        float pos_n[PLAN_NODES], pos_e[PLAN_NODES];
        float vel_n[PLAN_NODES], vel_e[PLAN_NODES];
        float acc_n[PLAN_NODES], acc_e[PLAN_NODES];
    };

    void init();
    bool enabled() const { return _enable != 0; }
    bool active() const { return _active; }

    // Mission intake: cheap change detection against the uploaded
    // mission and the configured cruise speed; the rebuild itself runs
    // on the solver thread while the solver is idle.  Scheduler context.
    void mission_poll(class AP_Mission &mission, float v_cap_ms);

    // Scheduler service.  `ingress_ready` must be the vehicle's answer to
    // "would a plan handed over right now actually be accepted by replay?".
    // It gates ENGAGE, and it has to: the mission clock starts at engage and
    // advances every solve, so engaging while the ingress is still refusing
    // plans (during the GUIDED takeoff, say) spends real ticks against a
    // vehicle that cannot move yet.  Measured, engaging at 2 m during takeoff
    // burned 35 ticks at zero along-track speed and drove the rigid-clock
    // assignment offset up by one EVERY tick before tracking even began.
    void update(bool ingress_ready);

    // Distance / speed the READY gate allows at the path start.  The mission
    // clock starts from the reference's own rest-engage IC at tau = 0, whose
    // lifted node 0 IS the path start, so engaging anywhere else hands the
    // position controller a step.  3 m is the record's own D_PROJ_MAX.
    static constexpr float READY_RADIUS_M = 3.0f;
    static constexpr float READY_SPEED_MS = 1.0f;

    // consume a freshly staged plan (MPCTrajReplay ingress)
    bool take_plan(Plan &out);

    // start/stop the mission clock (the calling mode owns replay entry)
    void engage();
    void disengage();

    // The READY gate (classic rest at the path start OR the ingress leg
    // from a near-hover).  Called by the owning mode each loop while it
    // wants engagement; true exactly when engage() ran.
    bool engage_try();
    // reset clock/model/plan state for a fresh engage (only while
    // inactive) — makes the sticky finished() re-armable per AUTO entry
    void reset_mission();
    bool finished() const { return _finished; }
    // mission-clock position of the last published plan (INT32_MIN
    // before the first publish); aligned 32-bit read, any thread
    int32_t mission_tau() const { return _tau_pub; }
    // the compiled run's bookkeeping (for the AUTO submode)
    uint16_t run_first_idx() const;
    uint16_t run_last_idx() const;
    int32_t tau_vertex(uint16_t wp_idx) const;
    bool mission_ready() const;

    // telemetry
    int32_t tick() const { return _tick; }
    int32_t offset() const { return _offset; }
    float last_lambda0() const { return _lambda0; }
    float last_solve_ms() const { return _solve_ms; }
    uint16_t overruns() const { return _overruns; }

private:
#if BS_BENCH_WCET
    void bench_thread();
#endif
    // The solve runs on its OWN thread, not a scheduler slot: it costs
    // hundreds of milliseconds on this target (135.7 ms at q = 2 for the
    // N = 30 pinned solve; 469.5 ms at q = 1 for the record's N = 40 full
    // solve) against a 250 ms period, so no scheduler budget can hold it.
    // AP_MPCSolver uses the same pattern.
    //
    // HISTORICAL CORRECTION, kept because the wrong version was believed for
    // two sessions: the arming refusal that prompted this thread was NOT
    // caused by the solve overrunning its slot.  It was internal error
    // 0x100000 (flow_of_control) raised at BOOT by AP_Scheduler.cpp:161,
    // which requires the vehicle task table to be sorted by ascending
    // priority -- the task had been inserted at priority 82 immediately
    // BEFORE one_hz_loop at 81.  The latched internal error is what refused
    // the arm.  The thread is still correct and still required; it simply was
    // not the fix.  Do not re-derive "task overrun" from an arming refusal
    // without first reading the internal-error line number.
    void solver_thread();
    // one full-horizon solve from U = 0 at tau = 0: the mission-start
    // initialization the interleaved loop of record seeds from.
    void seed_once();
    bool solve_once();
    // (re)build the phase-dependent objects, skipping the rebuild when the
    // phase has not moved.  The reference loop inits at tau at the top of a
    // tick and at (k+1) - o_off at the bottom; those two are the SAME phase
    // one tick apart, so caching turns two assemblies per tick into one.
    bool ensure_problem(int32_t phase);
#if HAL_LOGGING_ENABLED
    void write_log() const;
#else
    void write_log() const {}
#endif
    void lift_prefix(Plan &plan) const;

    AP_Int8  _enable;
    AP_Int8  _q;                 // Newton iterations per maintained solve
    AP_Int16 _qe;                // cap on the engage seed's Newton steps
    AP_Int8  _ir;                // restrict the setup to the pinned rows

    bool _inited;
    bool _thread_created;
    bool _active;
    bool _seeded;                // the engage full solve has run
    bool _finished;              // tau reached the end of the mission clock
    int32_t _tick;               // absolute clock index k
    int32_t _offset;             // rigid-clock assignment offset
    double _phase_ledger;        // along-track ledger driving the offset
    int32_t _problem_phase;      // phase the built problem belongs to (-1 none)
    int32_t _problem_npin;       // quadratic-row restriction it was built at
    uint32_t _next_solve_ms;
    uint32_t _anchor_ms;
    volatile int32_t _tau_pub = INT32_MIN;   // last published plan's tau

    double _xi[6];               // MODEL quotient state (never measured)
    double _U[BS_DRV_NV];        // the plan (decision vector)
    // The committed prefix, keyed by position in the window [k, k+M-1].  The
    // pinned solve provably cannot move it (B1p/B2i measure 0.0e+00), so
    // replaying it into _U is an identity — this array exists so the identity
    // is CHECKED each tick (_prefix_moved) rather than assumed.
    double _committed[3 * BS_DRV_M];
    float _prefix_moved;
    float _lambda0;
    // WALL TIME OF THE WHOLE TICK, not of the Newton call alone: micros()
    // around init + shift-append + pinned solve + lift, i.e. exactly the
    // quantity the on-target bench reports (bs_problem_init_pinned +
    // bs_shift_append + bs_newton_pinned, 135.7 ms worst case at q = 2 on an
    // IDLE board with the N = 30 corner tables; the record's N = 40 number
    // was 245.8 ms at q = 1) plus the O(BS_NV) bookkeeping.  On
    // Sim-on-Hardware the same core competes with the physics model and the
    // whole vehicle stack for one 480 MHz processor, so this number minus
    // the bench number IS the preemption cost, and it is the only way to
    // answer whether 4 Hz holds.
    float _solve_ms;
    // render-slice instrumentation (solver thread): time spent in
    // AP_BSMissionBuilder::render_to per tick.  _render_n / _max / _sum
    // count only calls that rendered a leg; _render_noop_us_max is the
    // check-only cost of the other calls.  Broadcast every 5 s.
    uint32_t _render_us_max, _render_us_sum, _render_noop_us_max;
    uint16_t _render_n;
    // Ticks whose wall time exceeded the 250 ms period.  Cumulative, logged
    // every tick, so the record carries both the event and its ordinal.  An
    // overrun is NOT a failure by itself: the committed prefix is M = 10
    // inputs = 2.5 s of standing plan and MPCTrajReplay keeps sampling it
    // while the late solve lands.  It is a failure only if the replay runs
    // off the end of that plan.
    uint16_t _overruns;
    // last periodic re-broadcast (boot statustexts are unrecoverable over USB)
    uint32_t _bcast_ms;
    // Line-search health of the last solve.  Held as plain counters rather
    // than a bs_newton_stats so this header need not pull in the C core.
    uint8_t _bt, _armijo_fail, _no_step;

    bool _plan_valid;
    Plan _plan;
    // guards the plan handoff: written on the solver thread, read from the
    // main loop through take_plan()
    HAL_Semaphore _plan_sem;
};

#endif // AP_BSOLVER_ENABLED
