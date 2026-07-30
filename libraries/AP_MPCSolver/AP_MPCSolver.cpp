#include "AP_MPCSolver.h"

#if AP_MPCSOLVER_ENABLED

#include <AP_AHRS/AP_AHRS.h>
#include <AP_HAL/AP_HAL.h>
#include <GCS_MAVLink/GCS.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
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

// receding-horizon cycle period. The Step-E harness re-solved per
// LOCAL_POSITION_NED sample (~10 Hz), but the replay tracks the accepted
// plan between solves, so the cycle only needs to refresh faster than the
// plan goes stale. 250 ms = what the H743 can actually sustain (measured
// WCET: rh cycle avg 240.6 / max 255.8 ms at N=30, f32, -Os).
static const uint32_t MPC_CYCLE_MS = 250;
// seed preconvergence budget (harness: capped at 30 — longer budgets let the
// flat u-reg objective drift the plan off R5; measured in SITL)
static const uint8_t MPC_PRECONV_MAX_ITERS = 30;
static const float MPC_PRECONV_TOL = 1.0e-4f;
static const uint8_t MPC_PRECONV_ITERS_PER_TICK = 5;

// PLAN-TIME PINNING: pin each engaged re-solve at the REPLAY'S CURRENT
// COMMANDED state (position/velocity) instead of the measured one, when the
// command is fresh (<= 300 ms) and the vehicle is within the divergence
// gate (2.0 m). Rationale (frozen-plan experiments, 07-24): the 10 Hz loop
// re-anchors every plan at the lag-carrying measured state, converting
// ordinary tracking lag into commanded slow-down (loop tempo gain < 1) —
// the mechanism behind the non-first-turn shift/crop-cut. Pinning at the
// commanded state gives the loop feedback in SPACE (divergence beyond the
// gate re-anchors at the measured state — wind-class disturbances still
// correct) without re-anchoring in TIME, and a plan solved from the
// commanded state re-bases the replay timeline WITHOUT a jump (clock
// continuity for free). a/chi were already plan-derived in build_x_meas;
// only pos/vel change source. NOTE: while plan-pinned, the spiral/TSTALL
// detector watches the commanded state, weakening forced-exit protection —
// evaluate before real flights.
#define MPC_PLAN_TIME_PIN 1
#define MPC_PLAN_PIN_MAX_AGE_MS 300
#define MPC_PLAN_PIN_MAX_DIV_M 2.0f

// DIAGNOSTIC: freeze the plan one engaged cycle after engage — the first
// cycle stages the bootstrap/accepted plan (its READY-converged iterate,
// live-pinned at the engage state), then ALL further cycling stops and the
// replay flies that single plan open-loop for the whole turn. t_rem() stays
// valid (wall-clock extrapolation from the first accept). Separates the
// per-cycle replan loop from execution/entry-state as the corruptor of
// non-first turns. Spiral/forced-exit protection is OFF while frozen —
// SITL experiments only.
#define MPC_FREEZE_PLAN 0

// base RTI iterations per rh cycle (the committed Step-E config is k=1;
// cs_rh escalates internally on a large innovation / recovery cycle)
static const int MPC_K_ITERS = 1;

// WCET bench build (M5, v2): at boot, on the MPC thread, run ARMED-style
// preconvergence then replay the RECORDED flight cycle sequence
// (mpc_bench_seq.h — a real turn's per-cycle x_meas/dt from the MPCC/MPCP
// telemetry, so cycle times cover a realistic work profile) plus a few
// large-innovation STRESS cycles driving the escalation/reject/recovery
// paths a nominal flight never exercises. EVERY sample is retained and
// dumped as chunked MPCB2 statustexts, re-broadcast every 10 s (one-shot
// boot statustexts are emitted before a USB host can connect and are
// lost). Capture host-side with scripts/mpc/experiments/cube_bench_v2.py
// (multi-boot). The dirtied solver state is irrelevant: arm_turn re-seeds
// via thread_rearm. NEVER leave enabled for flight builds.
#ifndef MPC_BENCH_WCET
#define MPC_BENCH_WCET 0
#endif

// LATENCY-INJECTION diagnostic (validation-ladder step B): pad every
// solver call on this thread to the MEASURED H743 cost (v2 bench:
// ~290 ms/iteration + ~60 ms rh overhead) so host SITL exhibits the
// target's timing behavior - escalated-cycle gaps, engage-commit
// latency, protection-window timing - with host-fast numerics. SITL
// experiments only; NEVER for flight or SimOnHW builds.
#ifndef MPC_SIM_SOLVE_LATENCY
#define MPC_SIM_SOLVE_LATENCY 0
#endif
#if MPC_SIM_SOLVE_LATENCY
static const uint32_t MPC_SIM_ITER_MS = 290;
static const uint32_t MPC_SIM_OVERHEAD_MS = 60;
#endif
#if MPC_BENCH_WCET
#include "mpc_bench_seq.h"
// stress cycles appended after the recorded sequence: growing canonical
// innovation from the last recorded state (the flight sequence is
// all-accepted; these force the expensive fallback/recovery machinery)
#define MPC_BENCH_STRESS_N 5
#endif

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
    _pin_tgt_ms(0),
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
    // keep in sync with MPC_REPLAY_TIME_SCALE in mpc_replay.cpp: the replay
    // dilates the plan timeline by SCALE, so plan-time remaining ticks at
    // 1/SCALE of wall clock (diagnostic; 1.0f = normal).
    const float scale = 1.0f;
    WITH_SEMAPHORE(const_cast<HAL_Semaphore&>(_sem));
    if (_state != State::ENGAGED || !_have_accepted) {
        return nanf("");
    }
    return _accepted_T - (AP_HAL::millis() - _accept_ms) * 1.0e-3f / scale;
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

void AP_MPCSolver::set_plan_pin_target(const Vector2f &pos_ne_m,
                                       const Vector2f &vel_ne_ms,
                                       uint32_t time_ms)
{
    WITH_SEMAPHORE(_sem);
    _pin_tgt_pos_ne_m = pos_ne_m;
    _pin_tgt_vel_ne_ms = vel_ne_ms;
    _pin_tgt_ms = time_ms;
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
#if MPC_BENCH_WCET
        {
            // v2 bench (see the gate comment at the top of this file):
            // per-sample retention, recorded flight sequence, stress tail.
            // INCREMENTAL: one setup step / one iterate / one cycle per
            // thread pass, so the paced dump below keeps streaming even if
            // a solver call hangs - everything recorded so far still gets
            // out, and the header's counts identify the hanging call.
            static uint8_t bench_phase;      // 0 setup, 1 pre, 2 cyc,
                                             // 3 done, 4 failed
            static uint8_t n_pre, n_cyc;
            static uint32_t t_pre[MPC_PRECONV_MAX_ITERS];
            static uint32_t t_cyc[BENCH_SEQ_N + MPC_BENCH_STRESS_N];
            // per-cycle code: bit0 accepted, bits1.. event (CS_RH_EVENT_*);
            // 0xE = cs_rh_step returned an error on this cycle
            static uint8_t c_cyc[BENCH_SEQ_N + MPC_BENCH_STRESS_N];
            static uint32_t bench_sent_ms;
            cs_rh *rh = (cs_rh *)_rh;
            if (bench_phase == 0) {
                uint8_t sel = 0;
                for (uint8_t i = 1; i < CS_SEED_FAMILY_COUNT; i++) {
                    if (fabsf(cs_seed_d[i] - BENCH_SEQ_D) <
                        fabsf(cs_seed_d[sel] - BENCH_SEQ_D)) {
                        sel = i;
                    }
                }
                int rc = cs_rh_init(rh, _solver, cs_seed_X[sel],
                                    cs_seed_T[sel], PLAN_N, MPC_K_ITERS);
                if (rc == CS_OK) {
                    rc = cs_rh_set_frame_d(rh, BENCH_SEQ_D);
                }
                if (rc == CS_OK) {
                    rc = cs_rh_set_ref_tracking(rh, 3.0f, 1, 1.5f);
                }
                bench_phase = (rc == CS_OK) ? 1 : 4;
            } else if (bench_phase == 1) {
                // ARMED preconvergence, one iterate per pass (sample 0 is
                // the cold cache - kept, it is a real boot cost)
                const uint32_t t0 = AP_HAL::micros();
                if (cs_solver_iterate(_solver) != CS_OK) {
                    bench_phase = 4;
                } else {
                    t_pre[n_pre++] = AP_HAL::micros() - t0;
                    if (n_pre >= MPC_PRECONV_MAX_ITERS) {
                        bench_phase = 2;
                    }
                }
            } else if (bench_phase == 2) {
                // one recorded-sequence or stress cycle per pass
                const uint8_t k = n_cyc;
                float x_meas[6];
                float dt_s;
                if (k < BENCH_SEQ_N) {
                    memcpy(x_meas, bench_seq_xmeas[k], sizeof(x_meas));
                    dt_s = bench_seq_dt[k];
                } else {
                    // stress tail: growing canonical innovation from the
                    // last recorded state (clamps mirror build_x_meas)
                    const uint8_t j = k - BENCH_SEQ_N + 1;
                    memcpy(x_meas, bench_seq_xmeas[BENCH_SEQ_N - 1],
                           sizeof(x_meas));
                    x_meas[0] = constrain_float(x_meas[0] - 0.5f * j,
                                                MPC_V_MIN + 1.0e-3f,
                                                MPC_V_INF - 5.0e-3f);
                    x_meas[2] -= 0.8f * j;                  // xi back
                    x_meas[3] = MAX(x_meas[3] +
                                    ((j & 1) ? 0.6f : -0.6f) * j,
                                    -5.0f);                 // eta zigzag
                    dt_s = 0.25f;
                }
                float T_plan = 0, T_rem_b = 0;
                int accepted = 0, event = CS_RH_EVENT_NONE;
                const uint32_t t0 = AP_HAL::micros();
                const int rc = cs_rh_step(rh, x_meas,
                                          (k == 0) ? 0.0f : dt_s,
                                          _X, _U, &T_plan, &T_rem_b,
                                          &accepted, &event);
                t_cyc[n_cyc] = AP_HAL::micros() - t0;
                c_cyc[n_cyc] = (rc != CS_OK)
                    ? 0xE
                    : uint8_t((accepted ? 1 : 0) |
                              (unsigned(event) << 1));
                n_cyc++;
                if (n_cyc >= BENCH_SEQ_N + MPC_BENCH_STRESS_N) {
                    bench_phase = 3;
                }
            }
            // paced round-robin dump: ONE line per 300 ms tick (a burst
            // overflows the GCS statustext queue), cycling header -> pre
            // chunks -> cyc chunks -> code string; runs from the first
            // pass on, so partial data streams even mid-bench
            const uint32_t bench_now = AP_HAL::millis();
            if (bench_now - bench_sent_ms >= 300) {
                bench_sent_ms = bench_now;
                static uint8_t dump_idx;
                char line[50];
                if (bench_phase == 4) {
                    GCS_SEND_TEXT(MAV_SEVERITY_INFO, "MPCB2 failed");
                } else {
                    const uint8_t n_pl = (n_pre + 3) / 4;
                    const uint8_t n_cl = (n_cyc + 3) / 4;
                    const uint8_t n_lines = 1 + n_pl + n_cl + 1;
                    const uint8_t li = dump_idx++ % n_lines;
                    if (li == 0) {
                        snprintf(line, sizeof(line),
                                 "MPCB2 h f%u p%u c%u d%.1f a%u r%u",
                                 unsigned(bench_phase), unsigned(n_pre),
                                 unsigned(n_cyc), (double)BENCH_SEQ_D,
                                 unsigned(_arena_bytes),
                                 unsigned(cs_rh_sizeof()));
                    } else if (li <= n_pl + n_cl) {
                        const bool is_pre = li <= n_pl;
                        const uint8_t base = (is_pre ? (li - 1)
                                                     : (li - 1 - n_pl)) * 4;
                        const uint8_t n = is_pre ? n_pre : n_cyc;
                        const uint32_t *ts = is_pre ? t_pre : t_cyc;
                        int len = snprintf(line, sizeof(line),
                                           "MPCB2 %c%02u",
                                           is_pre ? 'p' : 'c',
                                           unsigned(base));
                        for (uint8_t i = base; i < MIN(base + 4, n); i++) {
                            len += snprintf(line + len, sizeof(line) - len,
                                            " %u", unsigned(ts[i]));
                        }
                    } else {
                        int len = snprintf(line, sizeof(line), "MPCB2 e ");
                        for (uint8_t i = 0; i < n_cyc; i++) {
                            len += snprintf(line + len, sizeof(line) - len,
                                            "%X", unsigned(c_cyc[i]));
                        }
                        (void)len;
                    }
                    GCS_SEND_TEXT(MAV_SEVERITY_INFO, "%s", line);
                }
            }
        }
#endif
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
                _first_cycle = true;
                _last_cycle_snap_ms = 0;
                _next_cycle_ms = now;
                GCS_SEND_TEXT(MAV_SEVERITY_INFO, "MPC: engaged");
                // ENGAGE-STAGE (re-anchor diagnosis, CUBE-BENCH-NOTES):
                // stage the READY-converged iterate as a flyable plan NOW.
                // On target the first engaged solve can take >1 s; waiting
                // for it (the old first-commit bootstrap) let the bridge
                // expire into the Tier-2 zero-velocity hold and re-created
                // the engage transient at speed. Staging also feeds
                // plan-time pinning from cycle one, so the first solve
                // sees a small innovation instead of an escalation trigger.
                //
                // ANCHOR BY PROJECTION, not by timestamps: the anchor's
                // only job is to make the replay sample the plan AT the
                // vehicle. Timestamp guesses encode assumptions about
                // where the READY iterate's node 0 sits, and both guesses
                // tried were wrong on target (pin-time anchoring sampled
                // 1.1 s AHEAD - a measured 9.5 m command spike - because
                // the iterate stays on the seed frame at the engage gate).
                // Project the live snapshot onto the plan's canonical
                // polyline instead and back the anchor off by tau* T.
                float T_cur = 0;
                Snapshot snap_now;
                if (cs_solver_get_iterate(_solver, _X, _U, &T_cur) == CS_OK &&
                    T_cur > 0.5f && read_snapshot(snap_now)) {
                    bool sane = true;
                    for (uint16_t k = 0; k <= PLAN_N && sane; k++) {
                        const float v = _X[k * 6 + 0];
                        const float eta = _X[k * 6 + 3];
                        sane = isfinite(v) && v >= 0.0f && v < 20.0f &&
                               isfinite(eta) && fabsf(eta) < 40.0f;
                    }
                    // project the vehicle's canonical position onto the
                    // plan polyline (fractional node index -> tau*)
                    float vv, ppsi, pxi, peta;
                    world_to_canonical(snap_now, vv, ppsi, pxi, peta);
                    float best_d2 = FLT_MAX, tau_star = 0.0f;
                    for (uint16_t k = 0; k < PLAN_N; k++) {
                        const float ax = _X[k * 6 + 2], ay = _X[k * 6 + 3];
                        const float bx = _X[(k + 1) * 6 + 2];
                        const float by = _X[(k + 1) * 6 + 3];
                        const float dx = bx - ax, dy = by - ay;
                        const float L2 = MAX(dx * dx + dy * dy, 1.0e-9f);
                        const float sseg = constrain_float(
                            ((pxi - ax) * dx + (peta - ay) * dy) / L2,
                            0.0f, 1.0f);
                        const float ex = pxi - (ax + sseg * dx);
                        const float ey = peta - (ay + sseg * dy);
                        const float d2 = ex * ex + ey * ey;
                        if (d2 < best_d2) {
                            best_d2 = d2;
                            tau_star = (float(k) + sseg) / float(PLAN_N);
                        }
                    }
                    tau_star = constrain_float(tau_star, 0.0f, 0.9f);
                    // if the vehicle is not actually near the READY plan,
                    // staging it would command a jump - decline and let
                    // the bridge + first-cycle bootstrap handle engage
                    if (best_d2 > 9.0f) {
                        sane = false;
                    }
                    // corridor guard: never stage a plan whose row-to-row
                    // crossing sits over the crop (xi < 0 mid-crossing).
                    // The gate-anchored tracking redesign makes this
                    // unreachable by construction; it stays as the safety
                    // net for the invariant "only flight-legal shapes are
                    // flown" (v4 turn 3 flew a gate-relaxed diagonal at
                    // xi -7: CUBE-BENCH-NOTES correction 2026-07-30)
                    for (uint16_t k = 0; k <= PLAN_N && sane; k++) {
                        const float eta = _X[k * 6 + 3];
                        if (eta > 2.0f && eta < _rec.d - 2.0f &&
                            _X[k * 6 + 2] < -0.5f) {
                            sane = false;
                        }
                    }
                    uint32_t gen;
                    {
                        WITH_SEMAPHORE(_sem);
                        gen = _arm_generation;
                    }
                    Snapshot anchor;
                    memset(&anchor, 0, sizeof(anchor));
                    anchor.time_ms = snap_now.time_ms -
                        uint32_t(tau_star * T_cur * 1000.0f);
                    if (sane && stage_plan(_X, T_cur, anchor, gen)) {
                        _first_cycle = false;   // this IS the bootstrap
                        {
                            WITH_SEMAPHORE(_sem);
                            _accepted_T = T_cur;
                            _accept_ms = anchor.time_ms;
                            _have_accepted = true;
                        }
                        GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                                      "MPC: engage-staged (T %.1f tau %.2f "
                                      "d %.1f)", (double)T_cur,
                                      (double)tau_star,
                                      (double)sqrtf(best_d2));
                    }
                }
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
    _seed_xi0 = cs_seed_X[sel][2];
    _last_track_ms = 0;
    _track_snap_ms = 0;
    if (rc == CS_OK) {
        rc = cs_rh_set_frame_d(rh, d);
    }
    if (rc == CS_OK) {
        // Index-matched L2 pull, mode 1, w=3 (was 10 through the SimOnHW
        // campaign). The three-way lambda map (cert_lambda_map.py,
        // CUBE-BENCH-NOTES) measured w=10 NON-contracting at the mid-family
        // seed members (period-2 limit cycle, |lambda|>=0.97 at
        // d=13.75..14.5, 36>30 iters even at 14.1) while w=3 contracts at
        // ALL 10 members in 18-21 iters with |lambda|<=0.87 — and w~3 was
        // already the over-bulge sweep's sufficient weight. Tube
        // alternative (ABI 10 mode 3, hinge r0=0.25) is measured inert for
        // preconvergence (dead-band contains the seed => identical to
        // w=0). Crop floor still DISARMED (cs_rh_set_crop_bound).
        rc = cs_rh_set_ref_tracking(rh, 3.0f, 1, 1.5f);
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
#if MPC_SIM_SOLVE_LATENCY
        hal.scheduler->delay(MPC_SIM_ITER_MS);
#endif
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
    // GATE-ANCHORED tracking (redesign, 2026-07-30): the engage trigger IS
    // the xi crossing, so at the only moment this iterate will ever be
    // used the vehicle's along-row position is the handover point by
    // definition. Converge the DEVIATIONS that are genuinely unknown at
    // engage (eta settle residual, speed, heading branch, plan-carried
    // a/theta) but pin xi at the SEED START. This keeps tau0 = 0 so the
    // corridor leg gates stay REAL (no cs_rh_arm_ready_ref, no negative-
    // tau restretch, no CS_RH_READY_RLX gate relaxation), the plan is
    // never backward-stretched (T stays ~seed T), and the tracked iterate
    // remains a real-gated seed-frame solution BY CONSTRUCTION - the
    // engage-staged plan is then flight-legal without further checks.
    // Predecessor design pinned at the live xi (-30..-18): the stretched,
    // gate-relaxed re-solve migrated toward the min-time diagonal crop
    // cut and v4 flew it (CUBE-BENCH-NOTES correction 2026-07-30).
    x_meas[2] = _seed_xi0;
    if (cs_solver_set_pins(_solver, x_meas, _track_xN) != CS_OK) {
        return;
    }
    _track_snap_ms = snap.time_ms;
    for (uint8_t i = 0; i < 2; i++) {
        if (cs_solver_iterate(_solver) != CS_OK) {
            break;
        }
#if MPC_SIM_SOLVE_LATENCY
        hal.scheduler->delay(MPC_SIM_ITER_MS);
#endif
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
#if MPC_FREEZE_PLAN
        if (!_first_cycle && _have_accepted) {
            return;     // plan frozen: first engaged cycle staged + accepted
        }
#endif
        arm_generation = _arm_generation;
    }
    Snapshot snap;
    if (!read_snapshot(snap)) {
        return;
    }
#if MPC_PLAN_TIME_PIN
    {
        // pin at the replay's commanded state when fresh and close (see the
        // MPC_PLAN_TIME_PIN comment); otherwise fall through to the measured
        // snapshot (event-triggered re-anchor)
        Vector2f tp, tv;
        uint32_t tms;
        {
            WITH_SEMAPHORE(_sem);
            tp = _pin_tgt_pos_ne_m;
            tv = _pin_tgt_vel_ne_ms;
            tms = _pin_tgt_ms;
        }
        if (tms != 0 && snap.time_ms - tms <= MPC_PLAN_PIN_MAX_AGE_MS) {
            const float div = Vector2f(snap.pos_n - tp.x,
                                       snap.pos_e - tp.y).length();
            if (div <= MPC_PLAN_PIN_MAX_DIV_M) {
                snap.pos_n = tp.x;
                snap.pos_e = tp.y;
                snap.vel_n = tv.x;
                snap.vel_e = tv.y;
            } else {
                GCS_SEND_TEXT(MAV_SEVERITY_DEBUG,
                              "MPC: pin re-anchor (d %.1f)", (double)div);
            }
        }
    }
#endif
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
#if MPC_SIM_SOLVE_LATENCY
    hal.scheduler->delay(MPC_SIM_OVERHEAD_MS +
                         MPC_SIM_ITER_MS * uint32_t(cs_rh_last_keff(rh)));
#endif
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
        static uint32_t dbg_gen;
        // keyed to the arm generation, NOT _first_cycle: engage-staging
        // consumes the bootstrap flag before the first cycle runs, which
        // silently disabled this telemetry for every non-first turn
        if (dbg_gen != arm_generation) {
            dbg_gen = arm_generation;
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
            // e3 = the plan's node-3 eta (~0.3 s ahead = what the replay
            // will command next): discriminates target-led vs vehicle-led
            // in the approach eta drift (STATUSTEXT budget: drop rf/tau0)
            GCS_SEND_TEXT(MAV_SEVERITY_DEBUG,
                          "MPCC %u xi%.2f eta%.3f e3:%.3f b%.2f f%.2f a%d",
                          unsigned(dbg_n), (double)x_meas[2],
                          (double)x_meas[3], (double)_X[3 * 6 + 3],
                          (double)bmax, (double)xfl, accepted);
            // MPCP: the FULL pinned state, so the flight's exact x_meas
            // stream can be replayed through the host lib (residual-delta
            // hunt: flight plans lead the flare ~7 m earlier than the host
            // twin from nominal pins). Second line: STATUSTEXT ~50-char cap.
            GCS_SEND_TEXT(MAV_SEVERITY_DEBUG,
                          "MPCP %u v%.2f p%.3f a%.3f th%.3f dt%.3f",
                          unsigned(dbg_n), (double)x_meas[0],
                          (double)x_meas[1], (double)x_meas[4],
                          (double)x_meas[5], (double)dt_s);
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
