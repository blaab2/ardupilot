#include "AP_BSMissionBuilder.h"

#if AP_BSOLVER_ENABLED

#include <AP_Mission/AP_Mission.h>
#include <GCS_MAVLink/GCS.h>

#include <stdlib.h>
#include <string.h>

void AP_BSMissionBuilder::poll(const AP_Mission &mission, float v_cap_ms)
{
    const uint32_t stamp = mission.last_change_time_ms();
    if (_valid || _failed) {
        if (stamp == _stamp_ms && fabsf(v_cap_ms - _stamp_speed) < 0.01f) {
            return;
        }
        _valid = false;
        _failed = false;
    }
    _req = true;
    _stamp_ms = stamp;
    _stamp_speed = v_cap_ms;
}

// walk the mission into the FIRST maximal contiguous plain-waypoint run
bool AP_BSMissionBuilder::scan(AP_Mission &mission)
{
    _n_wp = 0;
    _run_first = _run_last = 0;
    AP_Mission::Mission_Command cmd;
    uint16_t search = 1;                       // item 0 is home
    bool in_run = false;
    float alt_lo = 0.0f, alt_hi = 0.0f;

    while (mission.get_next_nav_cmd(search, cmd)) {
        search = cmd.index + 1;
        Vector3f ned_m;
        const bool plain =
            cmd.id == MAV_CMD_NAV_WAYPOINT &&
            cmd.p1 == 0 &&
            cmd.content.location.get_vector_from_origin_NED_m(ned_m);
        if (!plain) {
            if (in_run) break;                 // the run ends here
            continue;                          // still before the run
        }
        // altitude band check (planar solver; replay pins altitude)
        if (!in_run) {
            alt_lo = alt_hi = ned_m.z;
        } else {
            if (ned_m.z < alt_lo) alt_lo = ned_m.z;
            if (ned_m.z > alt_hi) alt_hi = ned_m.z;
            if (alt_hi - alt_lo > 2.0f) break; // band exceeded: run ends
        }
        // DO_/CONDITION_ items between this and the previous run item
        // end the run before this waypoint
        if (in_run) {
            bool breaker = false;
            for (uint16_t j = _wp_idx[_n_wp - 1] + 1; j < cmd.index; ++j) {
                AP_Mission::Mission_Command between;
                if (mission.read_cmd_from_storage(j, between) &&
                    !AP_Mission::is_nav_cmd(between)) {
                    breaker = true;
                    break;
                }
            }
            if (breaker) break;
        }
        if (_n_wp >= BS_MB_MAX_WP) break;
        _vx[_n_wp] = (double)ned_m.x;
        _vy[_n_wp] = (double)ned_m.y;
        _wp_idx[_n_wp] = cmd.index;
        _n_wp++;
        if (!in_run) {
            _run_first = cmd.index;
            in_run = true;
        }
        _run_last = cmd.index;
    }
    if (_n_wp < 2) {
        return false;
    }
    _run_alt_d = 0.5f * (alt_lo + alt_hi);
    return true;
}

bool AP_BSMissionBuilder::build_pending(AP_Mission &mission)
{
    if (!_req) {
        return false;
    }
    _req = false;
    _valid = false;
    _failed = true;                            // until proven otherwise

    if (!scan(mission)) {
        // no eligible run: not an error, stock AUTO simply flies it
        if (AP_HAL::millis() - _txt_ms >= 10000) {
            _txt_ms = AP_HAL::millis();
            GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                          "BSLV mission: no waypoint run");
        }
        return true;
    }

    // total length for the sizing bound
    double len = 0.0;
    for (int i = 0; i + 1 < _n_wp; ++i) {
        const double dx = _vx[i + 1] - _vx[i], dy = _vy[i + 1] - _vy[i];
        len += sqrt(dx * dx + dy * dy);
    }
    const size_t need = bs_mission_size(_n_wp, len, (double)_v_cap);
    if (_block == nullptr || _block_len < need) {
        free(_block);
        _block = calloc(1, need);
        _block_len = _block ? need : 0;
        if (_block == nullptr) {
            GCS_SEND_TEXT(MAV_SEVERITY_WARNING,
                          "BSLV mission: no memory (%u B)", (unsigned)need);
            return true;
        }
    }

    bs_mission_params pp;
    pp.v_cap_ms = (double)_stamp_speed;
    _v_cap = _stamp_speed;
    const uint32_t t0 = AP_HAL::micros();
    const bs_mb_status st = bs_mission_build(_vx, _vy, _n_wp, &pp,
                                             _block, _block_len,
                                             &_mt, &_rep);
    const uint32_t dt_us = AP_HAL::micros() - t0;
    if (st != BS_MB_OK) {
        GCS_SEND_TEXT(MAV_SEVERITY_WARNING,
                      "BSLV mission: refused st%d G%d %.3g",
                      (int)st, _rep.gate, _rep.detail);
        return true;
    }
    _failed = false;
    _valid = true;
    GCS_SEND_TEXT(MAV_SEVERITY_INFO,
                  "BSLV mission: %d ticks %d seams wp %u-%u %.1fms",
                  _rep.n_ticks, _rep.n_seam - 1,
                  (unsigned)_run_first, (unsigned)_run_last,
                  (double)dt_us * 1e-3);
    return true;
}

int32_t AP_BSMissionBuilder::tau_vertex(uint16_t wp_idx) const
{
    if (!_valid) {
        return -1;
    }
    for (int i = 0; i < _n_wp; ++i) {
        if (_wp_idx[i] == wp_idx) {
            return (int32_t)_mt.wp_tick[i];
        }
    }
    return -1;
}

#endif // AP_BSOLVER_ENABLED
