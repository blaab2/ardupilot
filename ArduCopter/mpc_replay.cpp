#include "Copter.h"

#if MODE_GUIDED_ENABLED

// DIAGNOSTIC: replay time dilation. 1.0f = normal. 2.0f flies the SAME
// committed plan geometry at HALF speed: the timeline is sampled at
// elapsed/SCALE, velocities scale by 1/SCALE and accelerations by
// 1/SCALE^2 (chain rule); the drag-FF needs no change (it reads measured
// velocity). Pair of this define lives in AP_MPCSolver.cpp t_rem() — the
// exit countdown must tick at 1/SCALE or AUTO hands back mid-turn. Keep
// both in sync. Affects EVERY replayed trajectory (incl. the MAVLink
// path) — experiments only.
#define MPC_REPLAY_TIME_SCALE 1.0f

// TRUE tilt-rate (thrust-vector angular velocity) feedforward. 1 = derive the
// plan's jerk by finite-differencing the sampled KINEMATIC plan acceleration
// over a +0.1 s (real-time) window, map it to the angular velocity of the
// commanded specific-force direction f = (ax, ay, -g), and feed that rate to
// the attitude controller's roll/pitch command model so the shaped attitude
// target tracks the moving thrust vector without input-TC lag (~0.3-0.5 m
// cross-track measured without it). The tilt ANGLE target is unchanged — it
// still comes only from PSC's thrust vector (PID corrections and drag-FF
// included); ONLY the rate comes from the plan, so nothing is double-counted.
// The rate is clamped to 1.5 rad/s and forced to zero in the Tier-2 hold
// branch and for trajectories with fewer than 2 nodes. Works identically for
// the MAVLink/offboard TrajStream path (same buffer, same run()).
// 0 = bit-identical to the original input_thrust_vector_heading() hand-off.
#define MPC_TILT_RATE_FF 0

/*
 * MPCTrajReplay — implementation. This is the Step-E GUIDED TrajStream code
 * moved verbatim out of mode_guided.cpp (semantics preserved bit-for-bit;
 * see mpc_replay.h for the ownership/entry-point contract). The only
 * additions are the onboard-solver ingress (commit_pending) and the
 * AP_MPCSolver::update() service call at the top of run().
 */

// ingest one chunk of an onboard-MPC trajectory (NE-only, local NED)
// returns true when the final chunk commits the trajectory
bool MPCTrajReplay::handle_chunk(uint16_t traj_id, uint16_t num_points, uint16_t start_index, uint8_t count, uint16_t dt_ms,
                                 const float* pos_n, const float* pos_e,
                                 const float* vel_n, const float* vel_e,
                                 const float* acc_n, const float* acc_e,
                                 uint32_t time_boot_ms, float drag_k)
{
    // reject implausible sizes; dt_ms == 0 would collapse every node onto t=0 and
    // make sample() hold the final node (a raw full-plan position step)
    if (num_points == 0 || num_points > TrajBuffer::TRAJ_MAX || count == 0 || count > 6 || dt_ms == 0) {
        return false;
    }
    if ((uint32_t)start_index + count > num_points) {
        return false;
    }

    // never yank the vehicle out of a guided takeoff still in progress
    if (copter.mode_guided.is_taking_off()) {
        return false;
    }

    // the first chunk of a solve starts a fresh pending trajectory
    if (start_index == 0) {
        _traj.begin(traj_id, num_points, dt_ms, drag_k);
    }

    // chunks must belong to the pending trajectory and arrive in order, gap-free:
    // committing around a dropped chunk would publish stale (or, on the first
    // solve, zero-initialized) nodes in the hole. A stray-id chunk is dropped
    // without poisoning the assembly; an out-of-order same-id chunk poisons it
    // (recovered by the next start_index==0 chunk).
    if (!_traj.pending_valid || traj_id != _traj.pending_id) {
        return false;
    }
    if (start_index != _traj.pending_next) {
        _traj.pending_valid = false;
        return false;
    }

    // store this chunk's nodes into pending[]
    for (uint8_t i = 0; i < count; i++) {
        _traj.put(start_index + i, pos_n[i], pos_e[i], vel_n[i], vel_e[i], acc_n[i], acc_e[i]);
    }
    _traj.pending_next = start_index + count;

    // commit when every node up to the final one has been stored
    if (_traj.pending_next == _traj.pending_n) {
        return commit_assembled(time_boot_ms, true);
    }

    return false;
}

// poll the onboard solver for a newly accepted plan and commit it into the
// replay buffer. The plan's anchor_ms is the vehicle-clock timestamp of the
// state snapshot its solve started from — committing on it preserves the
// latency-anchored replay semantics of the MAVLink path.
void MPCTrajReplay::commit_pending()
{
#if AP_MPCSOLVER_ENABLED
    if (!copter.mpc_solver.enabled()) {
        return;
    }
    AP_MPCSolver::Plan plan;
    if (!copter.mpc_solver.take_plan(plan)) {
        return;
    }
    if (plan.n == 0 || plan.n > TrajBuffer::TRAJ_MAX || plan.dt_ms == 0) {
        return;
    }
    // same takeoff guard as the MAVLink ingress — but only while GUIDED is
    // actually the flying mode: ModeGuided's submode variable idles at
    // TakeOff (with the static takeoff_complete flag still false) in a
    // flight that never entered GUIDED, which would otherwise drop every
    // onboard plan of an AUTO MpcTurn (found in AUTO SITL: the replay
    // Tier-2 station-kept at the engage point because no plan ever
    // committed)
    if (copter.flightmode == &copter.mode_guided && copter.mode_guided.is_taking_off()) {
        return;
    }
    _traj.begin(uint16_t(_traj.pending_id + 1), plan.n, plan.dt_ms, plan.drag_k);
    for (uint16_t i = 0; i < plan.n; i++) {
        _traj.put(i, plan.pos_n[i], plan.pos_e[i], plan.vel_n[i], plan.vel_e[i],
                  plan.acc_n[i], plan.acc_e[i]);
    }
    _traj.pending_next = plan.n;
#if HAL_LOGGING_ENABLED
    // Capture the COMPLETE accepted trajectory every three seconds, rather
    // than only GUIP's current replay setpoint. These node sets are the exact
    // trajectory fan needed to see how successive replans move relative to
    // R5 and the flown track.
    if (_external_active && plan.anchor_ms >= _plan_snapshot_start_ms &&
        plan.anchor_ms - _plan_snapshot_start_ms >= _plan_snapshot_next_ms) {
        for (uint16_t i = 0; i < plan.n; i++) {
            copter.Log_Write_MPC_Plan(_plan_snapshot_id, uint8_t(i),
                                      plan.dt_ms, plan.pos_n[i], plan.pos_e[i]);
        }
        _plan_snapshot_id++;
        _plan_snapshot_next_ms += 3000U;
    }
#endif
    // no guided auto-start: the mode driving the onboard solver (AUTO
    // MpcTurn) owns its own replay submode entry (start()/run()/stop())
    commit_assembled(plan.anchor_ms, false);
#endif
}

// shared commit tail of both ingress paths (see mpc_replay.h)
bool MPCTrajReplay::commit_assembled(uint32_t time_boot_ms, bool auto_start_guided)
{
#if AP_FENCE_ENABLED
    // parity with the stock guided setters: reject a plan that leaves the fence
    const float node_alt_d_m = replay_active() ?
        (float)_entry_pos_d_m : copter.pos_control->get_pos_estimate_NED_m().z;
    for (uint16_t i = 0; i < _traj.pending_n; i++) {
        const Vector3p node_ned_m {(postype_t)_traj.pending[i].pos_n, (postype_t)_traj.pending[i].pos_e, (postype_t)node_alt_d_m};
        const Location node_loc = Location::from_ekf_offset_NED_m(node_ned_m, Location::AltFrame::ABOVE_ORIGIN);
        if (!copter.fence.check_location_within_fence(node_loc)) {
            LOGGER_WRITE_ERROR(LogErrorSubsystem::NAVIGATION, LogErrorCode::DEST_OUTSIDE_FENCE);
            _traj.pending_valid = false;
            return false;
        }
    }
#endif
    // replay time base: the sender's anchor timestamp (node 0 = the state at
    // the solve snapshot, stamped in the vehicle clock domain) compensates
    // solve+transport latency; fall back to arrival time if absent/implausible
    const uint32_t now_ms = millis();
    uint32_t base_ms = time_boot_ms;
    if (base_ms == 0 || base_ms > now_ms || now_ms - base_ms > 1000) {
        base_ms = now_ms;
    }
    _traj.commit(base_ms);
#if HAL_LOGGING_ENABLED
    // log the committed trajectory's anchor node (parity with the stock setters)
    copter.Log_Write_Guided_Position_Target(ModeGuided::SubMode::TrajStream,
        Vector3p{(postype_t)_traj.active[0].pos_n, (postype_t)_traj.active[0].pos_e, _entry_pos_d_m}, false,
        Vector3f{_traj.active[0].vel_n, _traj.active[0].vel_e, 0.0f},
        Vector3f{_traj.active[0].acc_n, _traj.active[0].acc_e, 0.0f});
#endif
    // auto-start the guided submode on first commit (mirrors how the PVA path
    // switches submode); subsequent commits just re-base the replay cursor
    if (auto_start_guided &&
        copter.mode_guided.submode() != ModeGuided::SubMode::TrajStream) {
        copter.mode_guided.trajstream_control_start();
    }
    return true;
}

// true while the replay controller is the active position controller: either
// GUIDED is in SubMode::TrajStream (the variable the pre-refactor code
// tested) or a non-guided mode has claimed the replay via start()/stop()
bool MPCTrajReplay::replay_active() const
{
    return _external_active ||
           copter.mode_guided.submode() == ModeGuided::SubMode::TrajStream;
}

// replay controller entry — mirrors pva_control_start() but sizes the NE/D
// limits to the plan envelope, holds heading (crab) and feeds the buffered
// trajectory RAW (no S-curve)
void MPCTrajReplay::start()
{
    // size the horizontal speed/accel envelope to cover the MPC plan
    // (>=12 m/s, accel >= g*tan(35deg) ~ 6.9 m/s/s) so the raw feed is not clamped
    const float ne_speed_ms = MAX(copter.wp_nav->get_default_speed_NE_ms(), 12.0f);
    const float ne_accel_mss = MAX(copter.wp_nav->get_wp_acceleration_mss(), GRAVITY_MSS * tanf(radians(35.0f)));
    copter.pos_control->NE_set_max_speed_accel_m(ne_speed_ms, ne_accel_mss);
    copter.pos_control->NE_set_correction_speed_accel_m(ne_speed_ms, ne_accel_mss);

    // vertical speeds and acceleration (altitude is held in this first cut)
    copter.pos_control->D_set_max_speed_accel_m(copter.wp_nav->get_default_speed_down_ms(), copter.wp_nav->get_default_speed_up_ms(), copter.wp_nav->get_accel_D_mss());
    copter.pos_control->D_set_correction_speed_accel_m(copter.wp_nav->get_default_speed_down_ms(), copter.wp_nav->get_default_speed_up_ms(), copter.wp_nav->get_accel_D_mss());

    // initialise the position controllers
    copter.pos_control->D_init_controller();
    copter.pos_control->NE_init_controller();

    // capture the entry altitude to hold while replaying (NE-only trajectory)
    _entry_pos_d_m = copter.pos_control->get_pos_estimate_NED_m().z;

    // drag-aware lean map: engage the k/m*|v|*v thrust-side feed-forward. The
    // velocity-PID I term was just seeded from the current attitude, which at
    // cruise already contains the drag share - preserve_output shifts the I term
    // down by the new FF so the commanded acceleration is continuous and drag is
    // not double-counted. From here the PID corrects only unmodeled effects.
    {
        const Vector2f vel_ne_ms = copter.pos_control->get_vel_estimate_NED_ms().xy();
        copter.pos_control->set_drag_accel_ff_NE_mss(vel_ne_ms * (_traj.drag_k * vel_ne_ms.length()), true);
    }

    // reset the Tier-2 station-keep latch
    _hold_active = false;

    // accepted-plan trajectory-fan logging schedule (0, 3, 6, ... s)
    _plan_snapshot_start_ms = millis();
    _plan_snapshot_next_ms = 0;
    _plan_snapshot_id = 0;

    // crab: hold current heading, do not turn the nose along the track
    Mode::auto_yaw.set_mode(Mode::AutoYaw::Mode::HOLD);

    // a non-guided owner (AUTO MpcTurn) claims the replay until stop();
    // for GUIDED the activity test is the guided submode itself
    _external_active = !copter.flightmode->in_guided_mode();
}

void MPCTrajReplay::stop()
{
    _external_active = false;
}

// replay controller — call from the owning mode's run() at the main loop
// rate (~400 Hz), after the mode's own disarmed/landed ground handling
void MPCTrajReplay::run()
{
    // service the onboard solver: publish the state snapshot it solves from
    // and commit any newly accepted plan into the replay buffer (both no-ops
    // for the MAVLink/offboard path)
#if AP_MPCSOLVER_ENABLED
    copter.mpc_solver.update();
#endif
    commit_pending();

    // set motors to full range
    copter.motors->set_desired_spool_state(AP_Motors::DesiredSpoolState::THROTTLE_UNLIMITED);

    Vector3p pos_ned_m;
    Vector3f vel_ned_ms;
    Vector3f accel_ned_mss;

    // hold the entry altitude on the D axis (first cut: NE-only trajectory)
    pos_ned_m.z = _entry_pos_d_m;

    const uint32_t now_ms = millis();
    const uint32_t age_ms = now_ms - _traj.commit_ms;

#if MPC_TILT_RATE_FF
    // thrust-vector angular velocity feedforward (NED, rad/s); stays zero in
    // the Tier-2 hold branch and for degenerate (<2 node) trajectories
    Vector3f tilt_rate_ff_ned_rads;
#endif

    float posN, posE, velN, velE, accN, accE;
    if (_traj.active_n == 0 ||
        age_ms > (uint32_t)(_traj.duration_ms() * MPC_REPLAY_TIME_SCALE) +
                 copter.mode_guided.get_timeout_ms()) {
        // Tier-2 failsafe: no trajectory yet, or the plan finished replaying and
        // no fresh solution arrived within the timeout on top of the plan's own
        // duration (a valid long plan must never be cut mid-replay). Station-keep
        // at the position captured when the hold engaged - re-pinning to the live
        // estimate every loop would drift with the estimate instead of holding.
        // Stay in-submode so a fresh trajectory resumes replay instantly.
        if (!_hold_active) {
            _hold_active = true;
            _hold_pos_ne_m = copter.pos_control->get_pos_estimate_NED_m().xy();
        }
        pos_ned_m.xy() = _hold_pos_ne_m;
        vel_ned_ms.zero();
        accel_ned_mss.zero();
    } else {
        // Tier-1: sample the buffered trajectory at true elapsed time (holds the
        // final node past the end until a new trajectory arrives).
        _hold_active = false;
        const float elapsed_s = age_ms * 0.001f / MPC_REPLAY_TIME_SCALE;
        _traj.sample(elapsed_s, posN, posE, velN, velE, accN, accE);
        pos_ned_m.x = posN;
        pos_ned_m.y = posE;
        vel_ned_ms.x = velN / MPC_REPLAY_TIME_SCALE;
        vel_ned_ms.y = velE / MPC_REPLAY_TIME_SCALE;
        accel_ned_mss.x = accN / (MPC_REPLAY_TIME_SCALE * MPC_REPLAY_TIME_SCALE);
        accel_ned_mss.y = accE / (MPC_REPLAY_TIME_SCALE * MPC_REPLAY_TIME_SCALE);

#if MPC_TILT_RATE_FF
        if (_traj.active_n >= 2) {
            // Jerk FF: finite-difference the SCALED commanded acceleration over
            // a REAL +h window. The commanded accel at real age t is
            // a_plan(t/S)/S^2, so (a_cmd(t+h) - a_cmd(t))/h = j_plan/S^3 — the
            // correctly time-dilated jerk with no extra scale factor. Only the
            // KINEMATIC plan accel enters here (not the drag-FF, not the PID).
            // CENTRAL difference: the first flight of the forward-FD variant
            // showed the ~h/2 phase LEAD as a ~0.3 m inside-shift of the whole
            // turn (fits -0.4..-0.7 vs -0.05..-0.5 without FF, at 5 m/s flare
            // speed) — the FF must carry zero phase to feed rate without
            // moving the flare. sample() clamps before node 0, so the early
            // half-window degrades to forward-FD only in the first 50 ms.
            const float h_s = 0.1f;     // window in REAL time
            float pN2, pE2, vN2, vE2, aN2, aE2, aN1, aE1;
            _traj.sample(elapsed_s + 0.5f * h_s / MPC_REPLAY_TIME_SCALE,
                         pN2, pE2, vN2, vE2, aN2, aE2);
            _traj.sample(elapsed_s - 0.5f * h_s / MPC_REPLAY_TIME_SCALE,
                         pN2, pE2, vN2, vE2, aN1, aE1);
            const float inv_scale2 = 1.0f / (MPC_REPLAY_TIME_SCALE * MPC_REPLAY_TIME_SCALE);
            const Vector3f jerk_ned_msss{
                (aN2 - aN1) * inv_scale2 / h_s,
                (aE2 - aE1) * inv_scale2 / h_s,
                0.0f};

            // Map jerk to the angular velocity of the commanded specific-force
            // direction: f = (ax, ay, -g) (NED, matching what PSC's
            // get_thrust_vector() builds from its accel target), b = f/|f|,
            // db/dt = (I - b b^T) j / |f|, omega = b x db/dt. omega has no
            // component along b — no spin about the thrust axis (yaw is
            // commanded separately and held).
            const Vector3f f_ned{accel_ned_mss.x, accel_ned_mss.y, -GRAVITY_MSS};
            const float f_len = f_ned.length();     // >= g, never zero
            const Vector3f b_ned = f_ned / f_len;
            const Vector3f db_dt = (jerk_ned_msss - b_ned * (b_ned * jerk_ned_msss)) / f_len;
            tilt_rate_ff_ned_rads = b_ned % db_dt;

            // safety clamp on the feedforward magnitude
            const float w_rads = tilt_rate_ff_ned_rads.length();
            if (w_rads > 1.5f) {
                tilt_rate_ff_ned_rads *= 1.5f / w_rads;
            }
        }
#endif
    }

    // drag-aware lean map: refresh the thrust-side drag feed-forward from the
    // measured velocity (k/m from the committed trajectory; 0 disables). The
    // streamed accelerations stay kinematic (= velocity-reference derivative);
    // the modeled drag share enters after the velocity PID so the PID's own
    // correction path is untouched and free for unmodeled disturbances.
    {
        const Vector2f vel_ne_ms = copter.pos_control->get_vel_estimate_NED_ms().xy();
        copter.pos_control->set_drag_accel_ff_NE_mss(vel_ne_ms * (_traj.drag_k * vel_ne_ms.length()), false);
    }

    // feed the controller RAW - no S-curve reshaping (this is the clean seam,
    // identical to how AC_WPNav feeds AC_PosControl)
    copter.pos_control->set_pos_vel_accel_NED_m(pos_ned_m, vel_ned_ms, accel_ned_mss);

#if AP_MPCSOLVER_ENABLED
    // plan-time pinning feed: the commanded target just fed to PSC, but only
    // while actually replaying a trajectory (Tier-2 hold must NOT pin the
    // solver to a station-keep target)
    if (!_hold_active && _traj.active_n > 0) {
        copter.mpc_solver.set_plan_pin_target(pos_ned_m.xy().tofloat(),
                                              vel_ned_ms.xy(), now_ms);
    }
#endif

    // run position controllers (add only P (pos->vel) + PID (vel->accel) correction)
    copter.pos_control->NE_update_controller();
    copter.pos_control->D_update_controller();

    // record the commanded targets for telemetry / wp reporting
    _target_pos_ned_m = pos_ned_m;
    _target_vel_ned_ms = vel_ned_ms;
    _target_accel_ned_mss = accel_ned_mss;

#if HAL_LOGGING_ENABLED
    // log the replayed target at 10 Hz (parity with the stock guided setters,
    // which log per set; here the "set" is the 400 Hz replay sampling)
    if (now_ms - _last_log_ms >= 100) {
        _last_log_ms = now_ms;
        copter.Log_Write_Guided_Position_Target(ModeGuided::SubMode::TrajStream, _target_pos_ned_m, false, _target_vel_ned_ms, _target_accel_ned_mss);
    }
#endif

    // call attitude controller with auto yaw (HOLD = crab)
#if MPC_TILT_RATE_FF
    {
        const Vector3f thrust_vector = copter.pos_control->get_thrust_vector();
        const AC_AttitudeControl::HeadingCommand heading = Mode::auto_yaw.get_heading();
        if (heading.heading_mode == AC_AttitudeControl::HeadingMode::Rate_Only) {
            // AutoYaw HOLD (set in start()) reports Rate_Only with zero rate, so
            // input_thrust_vector_heading() would dispatch to
            // input_thrust_vector_rate_heading_rads(thrust_vector, 0) — call the
            // feedforward overload of that same entry point directly. Rotate the
            // NED thrust-vector angular velocity into the body frame of the
            // controller's attitude TARGET (the frame thrust_vector_rotation_angles
            // expresses the roll/pitch error in, i.e. the frame the command-model
            // rates live in; _attitude_target carries the held heading). The
            // feedforward's body-z component is zeroed: heading is held, and the
            // yaw command model keeps its own rate input (0 in HOLD).
            Vector3f ang_vel_ff_body_rads =
                copter.attitude_control->get_attitude_target_quat().inverse() * tilt_rate_ff_ned_rads;
            ang_vel_ff_body_rads.z = 0.0f;
            copter.attitude_control->input_thrust_vector_rate_heading_rads(
                thrust_vector, heading.yaw_rate_rads, ang_vel_ff_body_rads, true);
        } else {
            // non-crab heading modes: original hand-off, no tilt-rate feedforward
            copter.attitude_control->input_thrust_vector_heading(thrust_vector, heading);
        }
    }
#else
    copter.attitude_control->input_thrust_vector_heading(copter.pos_control->get_thrust_vector(), Mode::auto_yaw.get_heading());
#endif
}

// TrajBuffer::begin - start assembling a new pending trajectory
void MPCTrajReplay::TrajBuffer::begin(uint16_t id, uint16_t n, uint16_t dt, float drag_k_perm)
{
    pending_id = id;
    pending_n = n;          // caller guarantees n <= TRAJ_MAX
    pending_next = 0;
    pending_valid = true;
    dt_ms = dt;
    pending_drag_k = drag_k_perm;
}

// TrajBuffer::put - store one node into pending[index]
bool MPCTrajReplay::TrajBuffer::put(uint16_t index, float pos_n, float pos_e, float vel_n, float vel_e, float acc_n, float acc_e)
{
    if (index >= TRAJ_MAX) {
        return false;
    }
    Sample &s = pending[index];
    s.t_s = index * dt_ms * 0.001f;
    s.pos_n = pos_n;
    s.pos_e = pos_e;
    s.vel_n = vel_n;
    s.vel_e = vel_e;
    s.acc_n = acc_n;
    s.acc_e = acc_e;
    return true;
}

// TrajBuffer::commit - publish pending[] as active[] and re-base the time cursor.
// Writes the payload (active[]) before publishing active_n/commit_ms so a
// concurrent sample() (single main-loop thread here) never sees a torn buffer.
void MPCTrajReplay::TrajBuffer::commit(uint32_t base_ms)
{
    for (uint16_t i = 0; i < pending_n; i++) {
        active[i] = pending[i];
    }
    active_id = pending_id;
    active_n = pending_n;
    drag_k = pending_drag_k;
    commit_ms = base_ms;
}

// TrajBuffer::duration_ms - replay duration of the active trajectory
uint32_t MPCTrajReplay::TrajBuffer::duration_ms() const
{
    if (active_n == 0) {
        return 0;
    }
    return (uint32_t)(active[active_n - 1].t_s * 1000.0f);
}

// TrajBuffer::sample - linearly interpolate the active trajectory at elapsed_s
bool MPCTrajReplay::TrajBuffer::sample(float elapsed_s, float &posN, float &posE, float &velN, float &velE, float &accN, float &accE) const
{
    if (active_n == 0) {
        return false;
    }

    // before the trajectory start: hold the first node
    if (elapsed_s <= active[0].t_s) {
        const Sample &s = active[0];
        posN = s.pos_n; posE = s.pos_e;
        velN = s.vel_n; velE = s.vel_e;
        accN = s.acc_n; accE = s.acc_e;
        return true;
    }

    // past the end (Tier-1): hold the final node position, zero vel/accel
    const Sample &last = active[active_n - 1];
    if (elapsed_s >= last.t_s) {
        posN = last.pos_n; posE = last.pos_e;
        velN = 0.0f; velE = 0.0f;
        accN = 0.0f; accE = 0.0f;
        return true;
    }

    // binary search for the interval [lo].t_s <= elapsed_s < [hi].t_s
    uint16_t lo = 0;
    uint16_t hi = active_n - 1;
    while (hi - lo > 1) {
        const uint16_t mid = (lo + hi) / 2;
        if (active[mid].t_s <= elapsed_s) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    const Sample &a = active[lo];
    const Sample &b = active[hi];
    const float span = b.t_s - a.t_s;
    const float u = (span > 1.0e-6f) ? (elapsed_s - a.t_s) / span : 0.0f;
    posN = a.pos_n + u * (b.pos_n - a.pos_n);
    posE = a.pos_e + u * (b.pos_e - a.pos_e);
    velN = a.vel_n + u * (b.vel_n - a.vel_n);
    velE = a.vel_e + u * (b.vel_e - a.vel_e);
    accN = a.acc_n + u * (b.acc_n - a.acc_n);
    accE = a.acc_e + u * (b.acc_e - a.acc_e);
    return true;
}

#endif  // MODE_GUIDED_ENABLED
