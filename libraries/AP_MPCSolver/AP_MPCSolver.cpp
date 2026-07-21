#include "AP_MPCSolver.h"

#if AP_MPCSOLVER_ENABLED

#include <AP_AHRS/AP_AHRS.h>
#include <AP_HAL/AP_HAL.h>
#include <GCS_MAVLink/GCS.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include "mpc_cs_defs.h"
#include "csolver/core/cs_solver.h"
#include "csolver/core/cs_model.h"
#include "csolver/core/cs_rh.h"
// the committed seed family (10 members, d in [13, 15]; G0.6). Included in
// this one TU only (static data).
#include "csolver/model/generated/cs_seed_family.h"
}

extern const AP_HAL::HAL &hal;

// ---- embedded-formulation constants (scripts/mpc/formulation.py; the plan's
// a_kin subtraction and the firmware drag-FF lean map MUST share K_OVER_M) ----
// k/m [1/m] — SN77 wind-corrected drag constant (log-analysis fit; the same
// 0.018017 the Step-E harness streams as MPC_TRAJECTORY.drag_k)
static const float MPC_K_OVER_M = 0.018017f;
static const float MPC_V_MIN = 0.5f;         // formulation V_MIN [m/s]
static const float MPC_V_INF = 11.8f;        // formulation V_INF [m/s]
static const float MPC_A_TRIM = 2.5060086f;  // g*tan(phi_trim) [m/s^2]
static const float MPC_A_CAP = 5.6638061f;   // g*tan(30 deg) [m/s^2]
static const float MPC_ENV_TOL = 0.01f;      // formulation ENV_TOL

// receding-horizon cycle period (the Step-E harness re-solves per
// LOCAL_POSITION_NED sample ~ 10 Hz)
static const uint32_t MPC_CYCLE_MS = 100;
// seed preconvergence budget (harness: capped at 30 — longer budgets let the
// flat u-reg objective drift the plan off R5; measured in SITL)
static const uint8_t MPC_PRECONV_MAX_ITERS = 30;
static const float MPC_PRECONV_TOL = 1.0e-4f;
static const uint8_t MPC_PRECONV_ITERS_PER_TICK = 5;

// base RTI iterations per rh cycle (the committed Step-E config is k=1;
// cs_rh escalates internally on a large innovation / recovery cycle)
static const int MPC_K_ITERS = 1;

const AP_Param::GroupInfo AP_MPCSolver::var_info[] = {
    // @Param: ENABLE
    // @DisplayName: Onboard turn-MPC enable
    // @Description: Enable the onboard time-optimal turn MPC (AUTO MpcTurn). Requires a reboot to allocate the solver.
    // @Values: 0:Disabled,1:Enabled
    // @RebootRequired: True
    // @User: Advanced
    AP_GROUPINFO_FLAGS("ENABLE", 1, AP_MPCSolver, _enable, 0, AP_PARAM_FLAG_ENABLE),

    // @Param: TURN_ANG
    // @DisplayName: Turn-detection heading tolerance
    // @Description: Maximum deviation of each corner from 90 degrees for a headland turn to be flown by the MPC (the pair must sum to a 180 degree row reversal).
    // @Units: deg
    // @Range: 0 45
    // @User: Advanced
    AP_GROUPINFO("TURN_ANG", 2, AP_MPCSolver, _turn_ang_deg, 5),

    // @Param: TURN_DMIN
    // @DisplayName: Minimum MPC turn row spacing
    // @Description: Detected turns with row spacing below this are not flown by the MPC (seed family band).
    // @Units: m
    // @Range: 10 20
    // @User: Advanced
    AP_GROUPINFO("TURN_DMIN", 3, AP_MPCSolver, _turn_dmin_m, 13),

    // @Param: TURN_DMAX
    // @DisplayName: Maximum MPC turn row spacing
    // @Description: Detected turns with row spacing above this are not flown by the MPC (seed family band).
    // @Units: m
    // @Range: 10 20
    // @User: Advanced
    AP_GROUPINFO("TURN_DMAX", 4, AP_MPCSolver, _turn_dmax_m, 15),

    AP_GROUPEND
};

// unwrap `angle` onto the branch nearest `ref` (embedded psi/theta are
// continuous; atan2 of the EKF velocity is not) — harness wrap_near
static float wrap_near(float angle, float ref)
{
    return ref + wrap_PI(angle - ref);
}

AP_MPCSolver::AP_MPCSolver() :
    _initialised(false),
    _thread_created(false),
    _arena(nullptr),
    _arena_bytes(0),
    _solver(nullptr),
    _rh(nullptr),
    _state(State::IDLE),
    _rearm_pending(false),
    _arm_generation(0),
    _engage_xi(0),
    _exit_forced(false),
    _accepted_T(0),
    _accept_ms(0),
    _have_accepted(false),
    _snap_seq(0),
    _snap_valid(false),
    _staged_seq(0),
    _taken_seq(0),
    _next_cycle_ms(0),
    _last_cycle_snap_ms(0),
    _first_cycle(false),
    _pre_iters(0),
    _pre_step(0)
{
    AP_Param::setup_object_defaults(this, var_info);
}

void AP_MPCSolver::init()
{
    if (_initialised || !enabled()) {
        return;
    }

    // one arena allocation at boot (before arming) — the csolver core is
    // allocation-free after this (CS_NO_LIBC_ALLOC); the cs_rh block is a
    // plain fixed-size struct (re-initialized per turn by cs_rh_init)
    _arena_bytes = cs_solver_min_arena(PLAN_N);
    if (_arena_bytes == 0) {
        GCS_SEND_TEXT(MAV_SEVERITY_ERROR, "MPC: no N=%u model variant", unsigned(PLAN_N));
        return;
    }
    _arena = calloc(1, _arena_bytes);
    _rh = calloc(1, cs_rh_sizeof());
    if (_arena == nullptr || _rh == nullptr) {
        GCS_SEND_TEXT(MAV_SEVERITY_ERROR, "MPC: arena alloc failed (%u B)",
                      unsigned(_arena_bytes + cs_rh_sizeof()));
        free(_arena); _arena = nullptr;
        free(_rh); _rh = nullptr;
        return;
    }

    int rc = 0;
    _solver = cs_solver_init(_arena, _arena_bytes, PLAN_N, &rc);
    if (_solver == nullptr) {
        GCS_SEND_TEXT(MAV_SEVERITY_ERROR, "MPC: solver init failed (rc %d)", rc);
        return;
    }

    if (!_thread_created) {
        if (!hal.scheduler->thread_create(FUNCTOR_BIND_MEMBER(&AP_MPCSolver::thread_main, void),
                                          "MPC", 8192, AP_HAL::Scheduler::PRIORITY_IO, 0)) {
            GCS_SEND_TEXT(MAV_SEVERITY_ERROR, "MPC: thread create failed");
            return;
        }
        _thread_created = true;
    }

    _initialised = true;
    GCS_SEND_TEXT(MAV_SEVERITY_INFO, "MPC: init (%u B arena)", unsigned(_arena_bytes));
}

// ---------------------------------------------------------------- mode API
bool AP_MPCSolver::arm_turn(const TurnRecord &rec)
{
    if (!enabled()) {
        return false;
    }
    if (!_initialised) {
        // lazy path: only while still disarmed (malloc + thread creation)
        if (hal.util->get_soft_armed()) {
            return false;
        }
        init();
        if (!_initialised) {
            return false;
        }
    }
    // record sanity: finite frame, unit-ish heading vector
    if (isnan(rec.d) || isnan(rec.origin_ne.x) || isnan(rec.origin_ne.y) ||
        isnan(rec.cos_h) || isnan(rec.sin_h)) {
        return false;
    }
    const float n2 = rec.cos_h * rec.cos_h + rec.sin_h * rec.sin_h;
    if (n2 < 0.99f || n2 > 1.01f) {
        return false;
    }
    {
        WITH_SEMAPHORE(_sem);
        ++_arm_generation;
        _rec = rec;
        _exit_forced = false;
        _have_accepted = false;
        _rearm_pending = true;      // the MPC thread performs the seed load
        _state = State::ARMED;
    }
    return true;
}

bool AP_MPCSolver::ready() const
{
    return _state == State::READY;
}

void AP_MPCSolver::engage()
{
    WITH_SEMAPHORE(_sem);
    if (_state == State::READY) {
        _state = State::ENGAGED;
    }
}

bool AP_MPCSolver::engaged() const
{
    return _state == State::ENGAGED;
}

void AP_MPCSolver::disarm()
{
    {
        WITH_SEMAPHORE(_sem);
        _state = State::IDLE;
        _rearm_pending = false;
        _have_accepted = false;
        _exit_forced = false;
    }
    {
        // drop any staged-but-untaken plan
        WITH_SEMAPHORE(_plan_sem);
        _taken_seq = _staged_seq;
    }
}

float AP_MPCSolver::t_rem() const
{
    WITH_SEMAPHORE(const_cast<HAL_Semaphore&>(_sem));
    if (_state != State::ENGAGED || !_have_accepted) {
        return nanf("");
    }
    return _accepted_T - (AP_HAL::millis() - _accept_ms) * 1.0e-3f;
}

bool AP_MPCSolver::exit_forced() const
{
    return _exit_forced;
}

float AP_MPCSolver::engage_xi() const
{
    return _engage_xi;
}

float AP_MPCSolver::current_xi() const
{
    if (_state == State::IDLE) {
        return nanf("");
    }
    Snapshot snap;
    if (!read_snapshot(snap)) {
        return nanf("");
    }
    float v, psi, xi, eta;
    world_to_canonical(snap, v, psi, xi, eta);
    return xi;
}

bool AP_MPCSolver::settled_on_row(float d_tol, float rate_tol) const
{
    if (_state != State::ENGAGED) {
        return false;
    }
    Snapshot snap;
    if (!read_snapshot(snap)) {
        return false;
    }
    // canonical cross-row position eta and rate eta_dot (world_to_canonical
    // math, extended to return the cross-row velocity component)
    const float dx = snap.pos_n - _rec.origin_ne.x;
    const float dy = snap.pos_e - _rec.origin_ne.y;
    float eta = -_rec.sin_h * dx + _rec.cos_h * dy;
    float eta_dot = -_rec.sin_h * snap.vel_n + _rec.cos_h * snap.vel_e;
    if (_rec.mirror) {
        eta = -eta;
        eta_dot = -eta_dot;
    }
    // fire at/above the row (headland side = the turn-around zone, safe),
    // NOT below it (the crop): the plan overshoots to the headland then
    // settles THROUGH the row toward the terminal, so the only clean
    // handover is the overshoot apex (eta >= d, cross-rate ~0). d_tol is the
    // small BELOW-row slack; above-row is bounded loosely so a normal
    // overshoot (~0.7 m) qualifies but a still-climbing crossing does not.
    return eta > _rec.d - d_tol && eta < _rec.d + 1.0f &&
           fabsf(eta_dot) < rate_tol;
}

// ------------------------------------------------------ main-loop plumbing
void AP_MPCSolver::update()
{
    if (!_initialised || _state == State::IDLE) {
        return;
    }
    Vector3f pos_ned, vel_ned;
    auto &ahrs = AP::ahrs();
    if (!ahrs.get_relative_position_NED_origin_float(pos_ned) ||
        !ahrs.get_velocity_NED(vel_ned)) {
        return;
    }
    // seqlock publish (single writer: the vehicle main loop)
    const uint32_t seq = _snap_seq + 1;             // odd: write in progress
    __atomic_store_n(&_snap_seq, seq, __ATOMIC_RELEASE);
    _snap.time_ms = AP_HAL::millis();
    _snap.pos_n = pos_ned.x;
    _snap.pos_e = pos_ned.y;
    _snap.pos_d = pos_ned.z;
    _snap.vel_n = vel_ned.x;
    _snap.vel_e = vel_ned.y;
    _snap.vel_d = vel_ned.z;
    __atomic_store_n(&_snap_seq, seq + 1, __ATOMIC_RELEASE);
    _snap_valid = true;
}

bool AP_MPCSolver::read_snapshot(Snapshot &out) const
{
    if (!_snap_valid) {
        return false;
    }
    for (uint8_t tries = 0; tries < 8; tries++) {
        const uint32_t s1 = __atomic_load_n(&_snap_seq, __ATOMIC_ACQUIRE);
        if (s1 & 1U) {
            continue;                               // write in progress
        }
        memcpy(&out, &_snap, sizeof(out));
        __atomic_thread_fence(__ATOMIC_ACQUIRE);
        const uint32_t s2 = __atomic_load_n(&_snap_seq, __ATOMIC_ACQUIRE);
        if (s1 == s2) {
            return true;
        }
    }
    return false;
}

bool AP_MPCSolver::take_plan(Plan &out)
{
    WITH_SEMAPHORE(_plan_sem);
    if (_staged_seq == _taken_seq) {
        return false;
    }
    out = _staged;
    _taken_seq = _staged_seq;
    return true;
}

// ------------------------------------------------------- frame transforms
void AP_MPCSolver::world_to_canonical(const Snapshot &snap, float &v,
                                      float &psi, float &xi, float &eta) const
{
    const float dx = snap.pos_n - _rec.origin_ne.x;
    const float dy = snap.pos_e - _rec.origin_ne.y;
    xi = _rec.cos_h * dx + _rec.sin_h * dy;
    float eta_u = -_rec.sin_h * dx + _rec.cos_h * dy;
    float vx_c = _rec.cos_h * snap.vel_n + _rec.sin_h * snap.vel_e;
    float vy_c = -_rec.sin_h * snap.vel_n + _rec.cos_h * snap.vel_e;
    if (_rec.mirror) {
        // mirrored problems solve canonically with the cross-row axis (and
        // hence all headings) sign-flipped in AND out
        eta_u = -eta_u;
        vy_c = -vy_c;
    }
    eta = eta_u;
    v = norm(vx_c, vy_c);
    psi = atan2f(vy_c, vx_c);
}

// ------------------------------------------------------- solver thread side
void AP_MPCSolver::thread_main()
{
    State last_state = State::IDLE;
    while (true) {
        hal.scheduler->delay(5);
        if (!_initialised) {
            continue;
        }
        if (_rearm_pending) {
            thread_rearm();
        }
        const State st = _state;
        const uint32_t now = AP_HAL::millis();
        switch (st) {
        case State::IDLE:
            break;
        case State::ARMED:
            thread_preconverge();
            break;
        case State::READY:
            thread_track();
            break;
        case State::ENGAGED:
            if (last_state != State::ENGAGED) {
                // engage: the first cycle runs immediately and stages its
                // returned accepted plan regardless of the accept flag (the
                // harness's `not sent_traj` first-commit bootstrap)
                _first_cycle = true;
                _last_cycle_snap_ms = 0;
                _next_cycle_ms = now;
                GCS_SEND_TEXT(MAV_SEVERITY_INFO, "MPC: engaged");
            }
            if (int32_t(now - _next_cycle_ms) >= 0) {
                _next_cycle_ms += MPC_CYCLE_MS;
                if (int32_t(now - _next_cycle_ms) > 0) {
                    _next_cycle_ms = now + MPC_CYCLE_MS;    // fell behind: re-pace
                }
                thread_cycle();
            }
            break;
        }
        last_state = st;
    }
}

void AP_MPCSolver::thread_rearm()
{
    TurnRecord rec;
    {
        WITH_SEMAPHORE(_sem);
        rec = _rec;
        _rearm_pending = false;
    }
    cs_rh *rh = (cs_rh *)_rh;
    // clamp the detected spacing into the seed-family band and select the
    // NEAREST family member as the seed (the runtime row-offset shift below
    // moves the corridor to the exact d; D-GENERALIZATION-NOTES)
    const float d = constrain_float(rec.d, cs_seed_d[0],
                                    cs_seed_d[CS_SEED_FAMILY_COUNT - 1]);
    uint8_t sel = 0;
    for (uint8_t i = 1; i < CS_SEED_FAMILY_COUNT; i++) {
        if (fabsf(cs_seed_d[i] - d) < fabsf(cs_seed_d[sel] - d)) {
            sel = i;
        }
    }
    int rc = cs_rh_init(rh, _solver, cs_seed_X[sel], cs_seed_T[sel],
                        PLAN_N, MPC_K_ITERS);
    memcpy(_track_xN, &cs_seed_X[sel][PLAN_N * 6], sizeof(_track_xN));
    _last_track_ms = 0;
    if (rc == CS_OK) {
        rc = cs_rh_set_frame_d(rh, d);
    }
    if (rc == CS_OK) {
        // R5-proximity reference tracking (w=10 L2) — hold the onboard
        // re-plan on the restretched R5 seed through the flat bulge<->time
        // min-time manifold. Crop floor DISARMED while the second-turn
        // rigid xi-shift is under diagnosis: the 2x2 translation fit showed
        // the "collapse" is R5's exact shape displaced -1.5 m in xi
        // (second-turn-specific, handedness-independent), so the crop
        // symptoms were the shift, not solver greed. Floor machinery stays
        // available via cs_rh_set_crop_bound once the shift is fixed.
        rc = cs_rh_set_ref_tracking(rh, 10.0f, 1, 1.5f);
    }
    if (rc != CS_OK) {
        GCS_SEND_TEXT(MAV_SEVERITY_ERROR, "MPC: frame load failed (rc %d)", rc);
        WITH_SEMAPHORE(_sem);
        _state = State::IDLE;
        return;
    }
    // effective engage gate: explicit record override, else the seed family's
    _engage_xi = (!isnan(rec.engage_xi) && !is_zero(rec.engage_xi))
                     ? rec.engage_xi : float(cs_rh_engage_xi(rh));
    _pre_iters = 0;
    _pre_step = FLT_MAX;
    GCS_SEND_TEXT(MAV_SEVERITY_INFO, "MPC: armed d=%.2f (seed %.2f)%s engage_xi=%.1f",
                  (double)d, (double)cs_seed_d[sel],
                  rec.mirror ? " (mirror)" : "", (double)_engage_xi);
}

void AP_MPCSolver::thread_preconverge()
{
    // live seed policy (harness RUN_UP): converge the turn solution before
    // engage; budget CAPPED at 30 iterations — longer budgets let the flat
    // u-reg objective drift the plan off R5 (measured in SITL)
    uint8_t burst = 0;
    while (_pre_iters < MPC_PRECONV_MAX_ITERS && _pre_step > MPC_PRECONV_TOL &&
           burst < MPC_PRECONV_ITERS_PER_TICK) {
        if (cs_solver_iterate(_solver) != CS_OK) {
            break;
        }
        cs_real step = 0, defect = 0, slack = 0, kappa = 0;
        int qp_status = 0, qp_iters = 0;
        if (cs_solver_get_info(_solver, &step, &defect, &slack, &kappa,
                               &qp_status, &qp_iters) == CS_OK && step >= 0) {
            _pre_step = step;
        }
        _pre_iters++;
        burst++;
    }
    if (_pre_iters >= MPC_PRECONV_MAX_ITERS || _pre_step <= MPC_PRECONV_TOL) {
        {
            WITH_SEMAPHORE(_sem);
            if (_state == State::ARMED) {
                _state = State::READY;
            }
        }
        GCS_SEND_TEXT(MAV_SEVERITY_INFO, "MPC: ready (%u iters, step %.1e)",
                      unsigned(_pre_iters), (double)_pre_step);
    }
}

bool AP_MPCSolver::stage_plan(const float *X, float T, const Snapshot &snap,
                              uint32_t arm_generation)
{
    // A solve from the previous turn can finish after ModeAuto has disarmed
    // and immediately armed a back-to-back turn. Never publish that stale
    // plan into the new turn's replay queue.
    WITH_SEMAPHORE(_sem);
    if (_state != State::ENGAGED || _arm_generation != arm_generation) {
        return false;
    }
    WITH_SEMAPHORE(_plan_sem);
    _staged.n = PLAN_NODES;
    _staged.dt_ms = MAX(1, uint16_t(lrintf(T / float(PLAN_N) * 1000.0f)));
    _staged.anchor_ms = snap.time_ms;   // solve-snapshot anchor (latency comp)
    _staged.drag_k = MPC_K_OVER_M;
    for (uint16_t k = 0; k < PLAN_NODES; k++) {
        const float v = X[k * 6 + 0];
        float psi = X[k * 6 + 1];
        const float xi = X[k * 6 + 2];
        float eta = X[k * 6 + 3];
        const float a = X[k * 6 + 4];
        float th = X[k * 6 + 5];
        if (_rec.mirror) {
            eta = -eta;
            psi = -psi;
            th = -th;
        }
        // canonical -> world SE(2)
        _staged.pos_n[k] = _rec.origin_ne.x + _rec.cos_h * xi - _rec.sin_h * eta;
        _staged.pos_e[k] = _rec.origin_ne.y + _rec.sin_h * xi + _rec.cos_h * eta;
        const float vx_c = v * cosf(psi);
        const float vy_c = v * sinf(psi);
        _staged.vel_n[k] = _rec.cos_h * vx_c - _rec.sin_h * vy_c;
        _staged.vel_e[k] = _rec.sin_h * vx_c + _rec.cos_h * vy_c;
        // KINEMATIC accel a_kin = a_thr - (k/m) v^2 v_hat (the Step-E
        // interface of record: the streamed accel FF is the velocity-
        // reference derivative; the firmware lean map adds the drag share
        // after the velocity PID from drag_k — streaming a_thr would
        // double-count drag and under-brake, measured dmax 1.35 -> 3.33 m)
        const float ax_c = a * cosf(th) - MPC_K_OVER_M * v * v * cosf(psi);
        const float ay_c = a * sinf(th) - MPC_K_OVER_M * v * v * sinf(psi);
        _staged.acc_n[k] = _rec.cos_h * ax_c - _rec.sin_h * ay_c;
        _staged.acc_e[k] = _rec.sin_h * ax_c + _rec.cos_h * ay_c;
    }
    _staged_seq++;
    return true;
}

// x_meas the harness way (plan-carried a/theta + psi branch at the nearest
// node of the current iterate; a MEASURED accel IC was tested and REGRESSED
// the closed loop — sitl_closed_loop.py TURN branch). Shared by the ENGAGED
// cycle and the READY tracking preconvergence.
bool AP_MPCSolver::build_x_meas(const Snapshot &snap, float x_meas[6],
                                float &xi_out)
{
    float v_raw, psi_raw, xi, eta;
    world_to_canonical(snap, v_raw, psi_raw, xi, eta);
    float T_cur = 0;
    if (cs_solver_get_iterate(_solver, _X, _U, &T_cur) != CS_OK) {
        return false;
    }
    uint16_t k0 = 0;
    float best_d2 = FLT_MAX;
    for (uint16_t k = 0; k <= PLAN_N; k++) {
        const float ddx = _X[k * 6 + 2] - xi;
        const float ddy = _X[k * 6 + 3] - eta;
        const float d2 = ddx * ddx + ddy * ddy;
        if (d2 < best_d2) {
            best_d2 = d2;
            k0 = k;
        }
    }
    const float psi_m = wrap_near(psi_raw, _X[k0 * 6 + 1]);
    const float a_pin = constrain_float(_X[k0 * 6 + 4],
                                        MPC_A_TRIM - MPC_ENV_TOL, MPC_A_CAP);
    const float chi = constrain_float(_X[k0 * 6 + 5] - psi_m,
                                      -0.2f + 1.0e-3f, M_PI - 1.0e-3f);
    x_meas[0] = constrain_float(v_raw, MPC_V_MIN + 1.0e-3f,
                                MPC_V_INF - 5.0e-3f);
    x_meas[1] = psi_m;
    x_meas[2] = xi;
    x_meas[3] = MAX(eta, -5.0f);
    x_meas[4] = a_pin;
    x_meas[5] = psi_m + chi;
    xi_out = xi;
    return true;
}

// READY-state tracking preconvergence (the turn-2 fix): the ARMED
// preconvergence converged the SEED's nominal entry, but the vehicle can
// arrive laterally offset (measured: turn 2 engages at eta ~ -0.15 after
// turn 1's handback leaves row 2 a shallow diagonal) and the engaged
// solver then pays a ~1 s slack-reject burst plus tau0/restretch re-adds
// re-converging. Anchor the iterate at the LIVE state (solve-only, no
// plan publishing, no cs_rh bookkeeping) in the final approach so engage
// starts from a converged-to-reality warm start whatever the offset.
void AP_MPCSolver::thread_track()
{
    const uint32_t now = AP_HAL::millis();
    if (now - _last_track_ms < 150) {
        return;
    }
    Snapshot snap;
    if (!read_snapshot(snap)) {
        return;
    }
    float x_meas[6], xi = 0;
    if (!build_x_meas(snap, x_meas, xi)) {
        return;
    }
    // only inside the final-approach zone: far from the gate the min-time
    // problem from the live state has a much longer horizon and would
    // reshape the warm start away from the turn solution
    if (xi < _engage_xi - 12.0f || xi >= _engage_xi) {
        return;
    }
    _last_track_ms = now;
    if (cs_solver_set_pins(_solver, x_meas, _track_xN) != CS_OK) {
        return;
    }
    for (uint8_t i = 0; i < 2; i++) {
        if (cs_solver_iterate(_solver) != CS_OK) {
            break;
        }
    }
}

void AP_MPCSolver::thread_cycle()
{
    uint32_t arm_generation;
    {
        WITH_SEMAPHORE(_sem);
        if (_state != State::ENGAGED) {
            return;
        }
        arm_generation = _arm_generation;
    }
    Snapshot snap;
    if (!read_snapshot(snap)) {
        return;
    }
    cs_rh *rh = (cs_rh *)_rh;
    float x_meas[6], xi_unused = 0;
    if (!build_x_meas(snap, x_meas, xi_unused)) {
        return;
    }

    // ---- one full rh cycle: spiral detection, restretch, leg-gated
    // corridor relaxations, pin, iterate, validate, accept/fallback and
    // T_rem aging all live in cs_rh (see cs_rh.h protocol) ----
    const float dt_s = (_last_cycle_snap_ms == 0)
        ? 0.0f : (snap.time_ms - _last_cycle_snap_ms) * 1.0e-3f;
    _last_cycle_snap_ms = snap.time_ms;
    float T_plan = 0, T_rem_out = 0;
    int accepted = 0, event = CS_RH_EVENT_NONE;
    const int rc = cs_rh_step(rh, x_meas, dt_s, _X, _U, &T_plan,
                              &T_rem_out, &accepted, &event);
    if (rc != CS_OK) {
        GCS_SEND_TEXT(MAV_SEVERITY_ERROR, "MPC: rh_step failed (rc %d)", rc);
        return;
    }

    bool report_forced_exit = false;
    {
        WITH_SEMAPHORE(_sem);
        /* ModeAuto can finish one turn and arm the next while this cycle is
         * solving. Before the generation guard, the old cycle then wrote its
         * nearly-expired T_rem into the new arm: at the second engage,
         * t_rem() was immediately -5.336 s and AUTO skipped the whole turn.
         * State alone is insufficient because the next turn can already be
         * ENGAGED; bind every result to the arm generation that started it. */
        if (_state != State::ENGAGED || _arm_generation != arm_generation) {
            return;
        }
        if (event == CS_RH_EVENT_EXIT_FORCED) {
            // TSTALL past the window exit: no solve ran; the mode must hand
            // control back to its stock guidance NOW (SpiralDetector rule)
            report_forced_exit = !_exit_forced;
            _exit_forced = true;
        } else {
            // t_rem() extrapolation base for the mode's exit decision
            _accepted_T = T_rem_out;
            _accept_ms = snap.time_ms;
            _have_accepted = true;
        }
    }
    if (event == CS_RH_EVENT_EXIT_FORCED) {
        if (report_forced_exit) {
            GCS_SEND_TEXT(MAV_SEVERITY_WARNING, "MPC: forced exit (T_rem %.2f)",
                          (double)T_rem_out);
        }
        return;
    }

    // TEMP xi-shift diagnosis: per-cycle canonical telemetry for the first
    // 25 engaged cycles — measured (xi, eta), the plan's canonical bulge
    // (max xi) and flare depth at eta=1.5 (interp), tau0, accept flag.
    {
        static uint8_t dbg_n;
        if (_first_cycle) {
            dbg_n = 0;
        }
        if (dbg_n < 25) {
            dbg_n++;
            float bmax = -99, xfl = 99;
            for (uint8_t k = 0; k <= PLAN_N; k++) {
                const float xk = _X[k * 6 + 2], ek = _X[k * 6 + 3];
                if (xk > bmax) {
                    bmax = xk;
                }
                if (k > 0) {
                    const float e0 = _X[(k - 1) * 6 + 3];
                    if ((e0 - 1.5f) * (ek - 1.5f) <= 0.0f &&
                        fabsf(ek - e0) > 1.0e-6f) {
                        const float t = (1.5f - e0) / (ek - e0);
                        const float xi_i = _X[(k - 1) * 6 + 2]
                            + t * (xk - _X[(k - 1) * 6 + 2]);
                        if (xi_i < xfl) {
                            xfl = xi_i;
                        }
                    }
                }
            }
            // reference (Xref) flare depth at eta=1.5: if this creeps with
            // tau0 the w=10 pull is the creep engine; if it stays R5-true
            // while the plan flare creeps, the corridor side is
            float rfl = 99;
            for (uint8_t k = 1; k <= PLAN_N; k++) {
                const float e0 = rh->Xref[(k - 1) * 6 + 3];
                const float ek = rh->Xref[k * 6 + 3];
                if ((e0 - 1.5f) * (ek - 1.5f) <= 0.0f &&
                    fabsf(ek - e0) > 1.0e-6f) {
                    const float t = (1.5f - e0) / (ek - e0);
                    const float xi_i = rh->Xref[(k - 1) * 6 + 2]
                        + t * (rh->Xref[k * 6 + 2]
                               - rh->Xref[(k - 1) * 6 + 2]);
                    if (xi_i < rfl) {
                        rfl = xi_i;
                    }
                }
            }
            GCS_SEND_TEXT(MAV_SEVERITY_DEBUG,
                          "MPCC %u xi%.2f eta%.2f b%.2f f%.2f rf%.2f t%.3f a%d",
                          unsigned(dbg_n), (double)x_meas[2],
                          (double)x_meas[3], (double)bmax, (double)xfl,
                          (double)rfl, (double)rh->tau0, accepted);
        }
    }

    // stage ONLY a genuinely new accepted plan (or the engage bootstrap /
    // a recovery reseed — the harness's streaming rule: re-sending an
    // unchanged plan would re-base the replay clock -> hold at the start)
    if (accepted || _first_cycle || event == CS_RH_EVENT_RESEED) {
        if (stage_plan(_X, T_plan, snap, arm_generation)) {
            _first_cycle = false;
        }
    } else {
        GCS_SEND_TEXT(MAV_SEVERITY_DEBUG, "MPC: reject (mask 0x%x)",
                      unsigned(cs_rh_last_reject(rh)));
    }
}

#endif  // AP_MPCSOLVER_ENABLED
