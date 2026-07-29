/*
 * MPCTrajReplay — copter-level replay of streamed/onboard MPC trajectories
 * (the Step-E TrajStream interface of record, moved out of ModeGuided so any
 * mode can run it: GUIDED SubMode::TrajStream delegates here; AUTO's MpcTurn
 * stage calls the same start()/run()).
 *
 * Double-buffer chunk assembly (TrajBuffer), clock-anchored replay with
 * linear interpolation, drag-FF lean map engage/refresh, Tier-1 (hold final
 * node until a fresh plan) / Tier-2 (duration-aware station-keep) staleness
 * fallback, fence check at commit, GUIP logging at 10 Hz.
 *
 * Plans arrive by two ingress paths that share one commit pipeline:
 *  - MAVLink MPC_TRAJECTORY chunks (GCS_MAVLink_Copter.cpp) -> handle_chunk()
 *    (auto-starts GUIDED SubMode::TrajStream on the committing chunk);
 *  - the onboard AP_MPCSolver -> commit_pending() polls take_plan() and
 *    commits with the solve snapshot's timestamp as the replay anchor (no
 *    submode auto-start: the calling mode owns its own replay entry).
 *
 * CONTRACT for a mode driving this helper:
 *  1. on entering its replay submode: call start() (sizes the NE/D
 *     controller envelopes, captures the entry altitude, engages the drag-FF
 *     lean map from the committed drag_k, points auto_yaw at HOLD);
 *  2. every loop of that submode: do the mode's own disarm/landed ground
 *     handling FIRST, then call run() (which also services
 *     AP_MPCSolver::update() + commit_pending());
 *  3. on leaving the submode: call stop() (clears the active flag used for
 *     fence-check altitude selection of later commits).
 * GUIDED keeps its bit-for-bit Step-E behavior through exactly this wiring
 * (mode_guided.cpp trajstream_control_start/run).
 */
#pragma once

#include <AP_Common/AP_Common.h>
#include <AP_HAL/AP_HAL_Boards.h>
#include <AP_Math/AP_Math.h>

class MPCTrajReplay {
public:
    MPCTrajReplay() {}

    CLASS_NO_COPY(MPCTrajReplay);

    // ingest one chunk of an onboard-MPC trajectory (NE-only, local NED).
    // Chunks sharing a traj_id assemble one trajectory; the chunk containing
    // the final node commits it and (if not already there) enters GUIDED
    // SubMode::TrajStream. Returns true when the trajectory commits.
    bool handle_chunk(uint16_t traj_id, uint16_t num_points, uint16_t start_index, uint8_t count, uint16_t dt_ms,
                      const float* pos_n, const float* pos_e,
                      const float* vel_n, const float* vel_e,
                      const float* acc_n, const float* acc_e,
                      uint32_t time_boot_ms = 0, float drag_k = 0.0f);

    // poll the onboard solver for a newly accepted plan and commit it into
    // the replay buffer (anchor = the plan's solve-snapshot timestamp).
    // Called from run(); safe to call from anywhere on the main loop.
    void commit_pending();

    // replay controller entry: call once when the mode enters its replay
    // submode (auto-called for GUIDED by the committing MAVLink chunk)
    void start();

    // replay controller: call at the loop rate from the owning mode's run()
    // AFTER the mode's own disarmed/landed ground handling
    void run();

    // mark the replay controller inactive (owning mode left its submode).
    // GUIDED does not need this (its activity test is the guided submode);
    // AUTO's MpcTurn stage must pair start()/stop().
    void stop();

    // targets of the last run() sample, for mode telemetry (wp reporting,
    // POSITION_TARGET_LOCAL_NED)
    const Vector3p &target_pos_NED_m() const { return _target_pos_ned_m; }
    const Vector3f &target_vel_NED_ms() const { return _target_vel_ned_ms; }
    const Vector3f &target_accel_NED_mss() const { return _target_accel_ned_mss; }

private:
    // Onboard-MPC trajectory replay buffer. Chunks fill pending[]; the
    // committing chunk copies pending->active and re-bases commit_ms, so
    // run() can sample by true elapsed time. A single main-loop thread
    // touches this (MAVLink ingestion, commit_pending and the mode run()
    // are all on the main loop), so no locking is needed - the payload
    // (active[]) is written before commit_ms/active_n are published.
    struct TrajBuffer {
        struct Sample {
            float t_s;      // node time from trajectory start (s)
            float pos_n;    // north position (m, local NED)
            float pos_e;    // east position (m, local NED)
            float vel_n;    // north velocity (m/s)
            float vel_e;    // east velocity (m/s)
            float acc_n;    // north acceleration (m/s/s)
            float acc_e;    // east acceleration (m/s/s)
        };
        static const uint16_t TRAJ_MAX = 128;   // 61-node plan + headroom

        Sample active[TRAJ_MAX];    // trajectory currently being replayed
        Sample pending[TRAJ_MAX];   // trajectory being assembled from chunks
        uint16_t active_n {0};      // number of valid nodes in active[]
        uint16_t pending_n {0};     // number of valid nodes in pending[]
        uint16_t pending_next {0};  // next expected node index (chunks must arrive in order, gap-free)
        bool pending_valid {false}; // false once a gap/out-of-order chunk poisons the assembly
        uint32_t commit_ms {0};     // replay time base: sender anchor time_boot_ms if given, else arrival millis()
        uint16_t active_id {0};     // traj_id of active[]
        uint16_t pending_id {0};    // traj_id of pending[]
        uint16_t dt_ms {0};         // nominal node spacing (ms)
        float drag_k {0};           // k/m (1/m) of the committed trajectory's airframe model (0 = no drag FF)
        float pending_drag_k {0};   // drag_k of the trajectory being assembled

        // begin assembling a new pending trajectory (dt must be non-zero, checked by the caller)
        void begin(uint16_t id, uint16_t n, uint16_t dt, float drag_k_perm);
        // store one node into pending[index]; returns false if index out of range
        bool put(uint16_t index, float pos_n, float pos_e, float vel_n, float vel_e, float acc_n, float acc_e);
        // publish pending[] as active[] and re-base the replay clock to base_ms
        void commit(uint32_t base_ms);
        // replay duration of the active trajectory in ms (0 if none)
        uint32_t duration_ms() const;
        // sample the active trajectory at elapsed_s (linear interpolation).
        // Before start clamps to the first node; past the end holds the final
        // position with zero velocity/acceleration (Tier-1 replay-until-new).
        // Returns false if there is no active trajectory.
        bool sample(float elapsed_s, float &posN, float &posE, float &velN, float &velE, float &accN, float &accE) const;
    } _traj;

    // shared tail of both ingress paths: fence-check pending[], re-base the
    // replay clock to base_ms (with plausibility fallback to now), publish,
    // log the anchor node. auto_start_guided: enter GUIDED SubMode::TrajStream
    // if not already there (MAVLink path only). Returns false if fence-rejected.
    bool commit_assembled(uint32_t time_boot_ms, bool auto_start_guided);

    // true while the replay controller is the active controller (fence-check
    // altitude source selection: entry alt vs current estimate)
    bool replay_active() const;

    // altitude (NED down) captured at start(); held while replaying (NE-only first cut)
    postype_t _entry_pos_d_m {0};

    // Tier-2 station-keep: position captured once when the failsafe hold engages
    // (re-pinning to the live estimate every loop would drift with the estimate)
    bool _hold_active {false};
    Vector2p _hold_pos_ne_m;

    // Engage-handover bridge: the pre-start() desired state (captured BEFORE
    // the controller init re-seeds the targets from the estimate) is
    // continued at constant velocity until the first plan commit, instead of
    // the zero-velocity Tier-2 hold that kicked the velocity PID with the
    // full cruise speed for the pre-commit tick(s). Expires into the normal
    // Tier-2 hold if no plan arrives in time.
    bool _bridge_active {false};
    uint32_t _bridge_start_ms {0};
    Vector2p _bridge_pos_ne_m;
    Vector2f _bridge_vel_ne_ms;
    Vector2f _bridge_accel_ne_mss;

    // Engage-handover accel-FF jerk blend: the WPNav->plan source switch
    // steps the kinematic accel FF discontinuously (the plan's jerk bound
    // governs only within a plan); the commanded FF slews toward the sampled
    // target at PSC's shaping jerk until caught up, then feeds raw again
    bool _blend_active {false};
    Vector2f _blend_accel_ne_mss;

    // rate limiter for GUIP dataflash logging of the replayed target
    uint32_t _last_log_ms {0};

    // Accepted-plan diagnostic snapshots: start()/each AUTO turn resets the
    // schedule; commit_pending logs all plan nodes at 0, 3, 6, ... seconds.
    uint32_t _plan_snapshot_start_ms {0};
    uint32_t _plan_snapshot_next_ms {0};
    uint8_t _plan_snapshot_id {0};

    // set by start()/stop(): a non-GUIDED mode (AUTO MpcTurn) is driving the replay
    bool _external_active {false};

    // last sampled targets (telemetry mirrors of the old guided_* statics)
    Vector3p _target_pos_ned_m;
    Vector3f _target_vel_ned_ms;
    Vector3f _target_accel_ned_mss;
};
