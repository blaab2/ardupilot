#include "AP_MPCTurnPlanner.h"

#if AP_MPCSOLVER_ENABLED

#include <AP_Mission/AP_Mission.h>

#include <math.h>

bool AP_MPCTurnPlanner::update(AP_Mission &mission, const AP_MPCSolver &solver)
{
    if (!solver.enabled()) {
        _num_turns = 0;
        _valid = false;
        return false;
    }
    const uint32_t change_ms = mission.last_change_time_ms();
    if (_valid && change_ms == _scan_change_ms) {
        return false;
    }
    scan(mission, solver);
    _scan_change_ms = change_ms;
    _valid = true;
    return true;
}

const AP_MPCSolver::TurnRecord *AP_MPCTurnPlanner::find_entry(uint16_t wp_index) const
{
    if (!_valid) {
        return nullptr;
    }
    for (uint8_t i = 0; i < _num_turns; i++) {
        if (_turns[i].entry_wp_idx == wp_index) {
            return &_turns[i];
        }
    }
    return nullptr;
}

void AP_MPCTurnPlanner::scan(AP_Mission &mission, const AP_MPCSolver &solver)
{
    _num_turns = 0;

    const float tol_deg = constrain_float(solver.turn_ang_deg(), 0.0f, 45.0f);
    const float dmin_m = solver.turn_dmin_m();
    const float dmax_m = solver.turn_dmax_m();

    // rolling window of the last four consecutive plain nav waypoints
    struct NavWP {
        uint16_t idx;       // mission index
        uint16_t p1;        // hold time [s] (corner guard)
        Vector3f ned_m;     // NED meters relative to the EKF origin
    } w[4];
    uint8_t have = 0;

    AP_Mission::Mission_Command cmd;
    uint16_t search = 1;                    // cmd 0 is the home item
    while (_num_turns < MAX_TURNS && mission.get_next_nav_cmd(search, cmd)) {
        search = cmd.index + 1;

        Vector3f ned_m;
        if (cmd.id != MAV_CMD_NAV_WAYPOINT ||
            !cmd.content.location.get_vector_from_origin_NED_m(ned_m)) {
            // a non-plain-waypoint nav item (takeoff, land, spline, loiter,
            // ...) or an unresolvable location (no EKF origin, terrain alt
            // without data) breaks the consecutive-waypoint chain
            have = 0;
            continue;
        }
        if (have == 4) {
            w[0] = w[1]; w[1] = w[2]; w[2] = w[3];
            have = 3;
        }
        w[have].idx = cmd.index;
        w[have].p1 = cmd.p1;
        w[have].ned_m = ned_m;
        have++;
        if (have < 4) {
            continue;
        }

        // ---- candidate window A,B,C,D ----
        const NavWP &A = w[0], &B = w[1], &C = w[2], &D = w[3];

        // strictly increasing mission indices (protects the DO-guard loop
        // below against DO_JUMP-reordered iteration)
        if (!(A.idx < B.idx && B.idx < C.idx && C.idx < D.idx)) {
            continue;
        }

        // altitudes ignored beyond a flat-altitude guard
        const float alt_hi = MAX(MAX(A.ned_m.z, B.ned_m.z), MAX(C.ned_m.z, D.ned_m.z));
        const float alt_lo = MIN(MIN(A.ned_m.z, B.ned_m.z), MIN(C.ned_m.z, D.ned_m.z));
        if (alt_hi - alt_lo > 0.5f) {
            continue;
        }

        const Vector2f ab = (B.ned_m - A.ned_m).xy();
        const Vector2f bc = (C.ned_m - B.ned_m).xy();
        const Vector2f cd = (D.ned_m - C.ned_m).xy();

        // degenerate rows (the connector length is checked by the d band)
        if (ab.length() < 1.0f || cd.length() < 1.0f) {
            continue;
        }

        // row spacing = the corner connector |BC| (perpendicular by the
        // 90-deg corner checks, so |BC| IS the projected row spacing)
        const float d_m = bc.length();
        if (d_m < dmin_m || d_m > dmax_m) {
            continue;
        }

        // corner angles: |heading change| at B and C within 90 +/- tol
        const float head_ab = atan2f(ab.y, ab.x);
        const float head_bc = atan2f(bc.y, bc.x);
        const float head_cd = atan2f(cd.y, cd.x);
        const float turn_b_deg = degrees(fabsf(wrap_PI(head_bc - head_ab)));
        const float turn_c_deg = degrees(fabsf(wrap_PI(head_cd - head_bc)));
        if (fabsf(turn_b_deg - 90.0f) > tol_deg || fabsf(turn_c_deg - 90.0f) > tol_deg) {
            continue;
        }

        // one continuous reversal: the same turning sign at both corners
        const float cross_b = ab % bc;
        const float cross_c = bc % cd;
        if (cross_b * cross_c <= 0.0f) {
            continue;
        }

        // entry and return rows anti-parallel within 2*tol
        const float anti_deg = degrees(fabsf(wrap_PI(head_cd - head_ab - M_PI)));
        if (anti_deg > 2.0f * tol_deg) {
            continue;
        }

        // conservative corner guards: no hold time on either corner and no
        // DO_/CONDITION_ command between B and D (resuming at D must not
        // skip pump/camera triggers)
        if (B.p1 > 0 || C.p1 > 0) {
            continue;
        }
        bool clean = true;
        for (uint16_t idx = B.idx + 1; idx < D.idx && clean; idx++) {
            AP_Mission::Mission_Command between;
            if (!mission.read_cmd_from_storage(idx, between) ||
                !AP_Mission::is_nav_cmd(between)) {
                clean = false;
            }
        }
        if (!clean) {
            continue;
        }

        // ---- detected: build the canonical-frame TurnRecord ----
        AP_MPCSolver::TurnRecord &rec = _turns[_num_turns];
        const Vector2f h = ab.normalized();     // +xi axis (entry row direction)
        rec.d = d_m;
        rec.origin_ne = B.ned_m.xy();           // first corner = canonical origin
        rec.cos_h = h.x;
        rec.sin_h = h.y;
        // eta of the second corner along (-sin_h, cos_h): positive = the
        // return row is on the canonical (+eta) side, negative = mirrored
        const float eta_c = -h.y * bc.x + h.x * bc.y;
        rec.mirror = (eta_c < 0.0f);
        rec.entry_wp_idx = B.idx;
        rec.exit_wp_idx = D.idx;
        rec.engage_xi = 0.0f;                   // 0 = the seed family's gate
        _num_turns++;
    }
}

#endif  // AP_MPCSOLVER_ENABLED
