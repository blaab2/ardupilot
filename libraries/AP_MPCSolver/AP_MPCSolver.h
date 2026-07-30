/*
 * AP_MPCSolver — onboard time-optimal turn MPC for the AUTO/MpcTurn feature.
 *
 * Wraps the vendored csolver C core (csolver/, f32 flight build — see
 * mpc_cs_defs.h and scripts/mpc/csolver/sync_firmware.sh) behind a small
 * vehicle-facing API:
 *
 *   boot:    init()                      — arena malloc + "MPC" thread (only
 *                                          when MPC_ENABLE=1; before arming)
 *   mission: arm_turn(TurnRecord)        — IDLE  -> ARMED  (seed loads +
 *                                          preconverges on the MPC thread)
 *            ready()                     — ARMED -> READY reached?
 *            engage()                    — READY -> ENGAGED (10 Hz rh cycles)
 *            t_rem() / exit_forced()     — exit decision inputs for the mode
 *            disarm()                    — any   -> IDLE
 *   plumbing (main loop): update()       — publish the AHRS world-NED state
 *                                          snapshot (seqlock; solver thread
 *                                          consumes it)
 *            take_plan(Plan&)            — consume a newly accepted plan in
 *                                          world NED (feed the TrajBuffer;
 *                                          MPCTrajReplay::commit_pending())
 *
 * Threading: ALL csolver calls happen on the dedicated "MPC" thread
 * (PRIORITY_IO, 8k stack). The main loop only touches the seqlock snapshot,
 * the staged-plan buffer (semaphore) and the small control flags.
 *
 * Frames: the solver works in the canonical turn frame ([v, psi, xi, eta,
 * a, theta]; xi along-row toward the corner line at xi=0, eta from row 1
 * toward row 2). TurnRecord carries the SE(2) world<->canonical transform +
 * a mirror flag; mirrored turns solve canonically with eta/psi/theta
 * sign-flipped on the way in AND out. Streamed accelerations are the
 * KINEMATIC a_kin = a_thr - (k/m) v^2 v_hat (the Step-E interface of
 * record); positions are absolute NED (EKF-origin meters).
 */
#pragma once

#include <AP_HAL/AP_HAL_Boards.h>

#ifndef AP_MPCSOLVER_ENABLED
#define AP_MPCSOLVER_ENABLED 1
#endif

#if AP_MPCSOLVER_ENABLED

#include <AP_Common/AP_Common.h>
#include <AP_HAL/AP_HAL.h>
#include <AP_Math/AP_Math.h>
#include <AP_Param/AP_Param.h>

class AP_MPCSolver {
public:
    AP_MPCSolver();

    CLASS_NO_COPY(AP_MPCSolver);

    // parameter table (Agent C wires this as an MPC_ subgroup of ParametersG2:
    //   AP_SUBGROUPINFO(mpc_solver, "MPC_", <idx>, ParametersG2, AP_MPCSolver)
    // giving MPC_ENABLE / MPC_TURN_ANG / MPC_TURN_DMIN / MPC_TURN_DMAX)
    static const struct AP_Param::GroupInfo var_info[];

    // solver horizon: N = 30 intervals -> 31-node plans
    static const uint8_t PLAN_N = 30;
    static const uint8_t PLAN_NODES = PLAN_N + 1;

    // world->canonical SE(2)+mirror frame of one detected turn (built by
    // ModeAuto's turn detector, consumed by arm_turn)
    struct TurnRecord {
        float d;                // detected row spacing [m] (clamped into the seed-family band [13, 15])
        Vector2f origin_ne;     // world NED position of the canonical origin (row-1 end corner, xi=0/eta=0) [m]
        float cos_h, sin_h;     // unit vector of the +xi axis (row direction toward the corner line) in world NE
        bool mirror;            // true: the turn curves to the -eta side; solved canonically with eta/psi/theta sign-flipped in AND out
        uint16_t entry_wp_idx;  // mission bookkeeping (opaque to the solver)
        uint16_t exit_wp_idx;   //   "
        float engage_xi;        // engage gate override [m]; NaN or 0 = use the seed family's gate (recommended)
    };

    // one accepted plan in world NED, ready for the replay TrajBuffer
    struct Plan {
        uint16_t n;             // number of nodes (== PLAN_NODES)
        uint16_t dt_ms;         // node spacing [ms]
        uint32_t anchor_ms;     // vehicle-clock timestamp of the solve's state snapshot (TrajBuffer commit_ms anchor)
        float drag_k;           // k/m [1/m] of the model (engages the firmware drag-FF lean map)
        float pos_n[PLAN_NODES], pos_e[PLAN_NODES];   // absolute NED [m]
        float vel_n[PLAN_NODES], vel_e[PLAN_NODES];   // [m/s]
        float acc_n[PLAN_NODES], acc_e[PLAN_NODES];   // KINEMATIC a_kin [m/s^2]
    };

    // boot-time init: one arena allocation + the "MPC" solver thread.
    // No-op unless MPC_ENABLE=1. Call before arming (arm_turn() also calls
    // it lazily, but only accepts that path while still disarmed).
    void init();

    bool enabled() const { return _enable != 0; }

    // ---- ModeAuto (Agent C) API — call order: arm_turn -> poll ready ->
    //      engage -> poll t_rem/exit_forced -> disarm ----
    // returns false if disabled, not initializable, or the record is invalid
    bool arm_turn(const TurnRecord &rec);
    bool ready() const;         // seed preconverged; engage() will be honored
    void engage();              // start the 10 Hz receding-horizon cycles
    bool engaged() const;
    void disarm();              // abort/finish: back to IDLE (plans stop)
    // accepted-plan remaining time [s] (NaN unless ENGAGED). Exit the turn
    // (hand back to waypoint nav) when this drops below the T_EXIT the
    // caller uses (harness value: 0.8 s).
    float t_rem() const;
    // TSTALL latch (plans re-adding T ~ the restretch shrink near the exit):
    // once true the caller must force the handback NOW (SpiralDetector rule)
    bool exit_forced() const;
    // effective engage gate [m]: record override or the seed family's xi
    float engage_xi() const;
    // vehicle xi [m] in the armed canonical frame from the latest snapshot
    // (NaN while IDLE or before the first update()) — engage when
    // current_xi() >= engage_xi()
    float current_xi() const;

    // true when the vehicle is settled on the RETURN row: canonical eta
    // within d_tol of the row (d) AND the cross-row rate |eta_dot| below
    // rate_tol. The turn-exit gate uses this so handover happens on the
    // straight with ~zero cross-track motion - otherwise the resumed WPNav
    // projects the settling-oscillation velocity into the inter-row area.
    bool settled_on_row(float d_tol, float rate_tol) const;

    // ---- main-loop plumbing ----
    // publish the AHRS world-NED pos/vel snapshot for the solver thread.
    // Call at >= 50 Hz while non-IDLE (MPCTrajReplay::run() calls it; before
    // the replay is active the controlling mode must call it directly).
    void update();
    // consume a newly staged accepted plan (false if none since last call)
    bool take_plan(Plan &out);
    // replay-side feed for plan-time pinning (MPC_PLAN_TIME_PIN): the
    // replay's CURRENT commanded NE position/velocity target and its
    // vehicle-clock timestamp. Called from MPCTrajReplay::run() every
    // Tier-1 loop; consumed by the engaged solve to pin x0 at the plan's
    // predicted state instead of the measured one (freshness + divergence
    // gated in thread_cycle).
    void set_plan_pin_target(const Vector2f &pos_ne_m,
                             const Vector2f &vel_ne_ms, uint32_t time_ms);

    // turn-detector parameters (read by Agent C's ModeAuto scanner)
    float turn_ang_deg() const { return _turn_ang_deg; }
    float turn_dmin_m() const { return _turn_dmin_m; }
    float turn_dmax_m() const { return _turn_dmax_m; }

private:
    enum class State : uint8_t {
        IDLE = 0,       // no turn armed
        ARMED,          // record loaded; seed selected; preconverging
        READY,          // preconverged; waiting for engage()
        ENGAGED,        // 10 Hz rh cycles; plans staged for the replay
    };

    // world-NED state snapshot published by update() (single writer: main
    // loop; single reader: MPC thread) — classic seqlock, odd seq = writing
    struct Snapshot {
        uint32_t time_ms;       // vehicle clock at sampling (anchor domain)
        float pos_n, pos_e, pos_d;
        float vel_n, vel_e, vel_d;
    };

    void thread_main();
    void thread_rearm();        // ARMED entry work (seed load) on the thread
    void thread_preconverge();  // ARMED -> READY
    void thread_track();        // READY: live-state tracking preconvergence
    void thread_cycle();        // one ENGAGED 10 Hz rh cycle
    bool build_x_meas(const Snapshot &snap, float x_meas[6], float &xi_out);
    bool read_snapshot(Snapshot &out) const;
    bool stage_plan(const float *X, float T, const Snapshot &snap,
                    uint32_t arm_generation);

    // canonical<->world helpers of the armed frame (_rec)
    void world_to_canonical(const Snapshot &snap, float &v, float &psi_w,
                            float &xi, float &eta) const;

    // parameters
    AP_Int8  _enable;
    AP_Float _turn_ang_deg;
    AP_Float _turn_dmin_m;
    AP_Float _turn_dmax_m;

    // init
    bool _initialised;
    bool _thread_created;
    void *_arena;               // csolver arena (malloc'd once at init)
    size_t _arena_bytes;
    struct cs_solver *_solver;  // carved from _arena at init
    void *_rh;                  // cs_rh block (malloc'd at init; the cs_rh
                                // typedef names an anonymous struct, so the
                                // pointer stays untyped here)

    // control state (guarded by _sem where written cross-thread)
    HAL_Semaphore _sem;
    volatile State _state;
    bool _rearm_pending;        // arm_turn asked the thread to (re)load
    uint32_t _arm_generation;   // increments on every arm; rejects results
                                // completed by a previous turn's cycle
    TurnRecord _rec;            // armed frame (stable while non-IDLE)
    float _engage_xi;           // effective gate (filled on the thread)
    volatile bool _exit_forced; // TSTALL latch
    // accepted-plan timeline (thread writes under _sem; t_rem reads)
    float _accepted_T;
    uint32_t _accept_ms;        // snapshot time of the accepted cycle
    bool _have_accepted;

    // state snapshot seqlock (main-loop writer, thread reader)
    volatile uint32_t _snap_seq;
    Snapshot _snap;
    bool _snap_valid;           // one-way latch: first good publish

    // staged plan hand-off (thread writer, main-loop reader; _plan_sem)
    HAL_Semaphore _plan_sem;
    Plan _staged;
    uint32_t _staged_seq;       // bumped per staged plan
    uint32_t _taken_seq;

    // solver-thread locals (no cross-thread access)
    float _X[(PLAN_N + 1) * 6]; // iterate / accepted-plan scratch
    float _U[PLAN_N * 2];
    uint32_t _next_cycle_ms;
    uint32_t _last_cycle_snap_ms;   // snapshot time of the previous rh cycle
    // plan-time pin feed (set_plan_pin_target; guarded by _sem)
    Vector2f _pin_tgt_pos_ne_m;
    Vector2f _pin_tgt_vel_ne_ms;
    uint32_t _pin_tgt_ms;
    bool _first_cycle;          // engage bootstrap: stage the first plan
                                // regardless of accept (harness `sent_traj`)
    uint8_t _pre_iters;         // preconverge bookkeeping
    float _track_xN[6];         // seed terminal pin for READY tracking
    uint32_t _last_track_ms;    // READY tracking cadence
    uint32_t _track_snap_ms;    // snapshot time of the last READY tracking
                                // pin (diagnostic)
    float _seed_xi0;            // seed node-0 xi — READY tracking pins the
                                // along-row coordinate HERE (the handover
                                // point), never at the live xi
    float _pre_step;
};

#endif  // AP_MPCSOLVER_ENABLED
