/* AP_BSMissionBuilder — the AP-side half of the runtime mission builder.
 *
 * Walks the uploaded AP_Mission into a waypoint polyline, feeds it to the
 * pure-C bs_mission_build core (bsolver/core/bs_mission_build.h), owns the
 * one heap block the tables live in, and caches on
 * mission.last_change_time_ms() + the configured cruise speed (the
 * AP_MPCTurnPlanner::update idiom).  Nothing here runs in the per-tick
 * control path: rebuild_poll() is a cheap timestamp compare from a
 * scheduler slot, and the actual build (sub-millisecond measured on host,
 * budgeted 200 ms) runs on the SOLVER THREAD while the solver is idle.
 *
 * RUN EXTRACTION (v1): the FIRST maximal contiguous run of plain
 * MAV_CMD_NAV_WAYPOINT items (no hold time, resolvable location, altitude
 * within a 2 m band, no DO/CONDITION movers in between).  Commands before
 * and after the run fly stock AUTO; anything that breaks the run simply
 * ENDS it — enabling the bsolver never makes a flyable mission unflyable.
 */
#pragma once

#include "bs_defs.h"

#if AP_BSOLVER_ENABLED

#include <AP_Common/AP_Common.h>
#include <AP_HAL/AP_HAL.h>

extern "C" {
#include "bsolver/core/bs_mission_build.h"
}

class AP_Mission;

class AP_BSMissionBuilder {
public:
    // cheap change detection; sets the rebuild request when the mission
    // or the configured speed moved.  Scheduler context.
    void poll(const class AP_Mission &mission, float v_cap_ms);

    // perform a pending rebuild.  SOLVER-THREAD context, call only while
    // the solver is inactive.  Returns true when a build ran (pass or
    // refuse).
    bool build_pending(class AP_Mission &mission);

    // the published tables; nullptr until a build succeeded for the
    // CURRENT mission/speed stamps
    const bs_mission_tables *tables() const {
        return _valid ? &_mt : nullptr;
    }

    // STREAMING: render the per-tick ring far enough that every read at
    // mission tick `upto` (horizon tip inclusive) hits a valid slot.
    // Solver-thread context, microseconds per call.  Returns false only
    // on the impossible starvation case (ring smaller than the horizon).
    bool render_to(int32_t upto);
    int32_t rendered_to() const { return _cursor; }
    void ring_reset();
#if BS_STAGE2_REPACE
    // Stage 2 (c): one slow tick's verdict inputs (model lead, clock input,
    // worst face violation of the applied input; solver-thread context)
    void repace_feed(double el, double nu, double viol) {
        _rp_n++;
        _rp_nu_sum += fabs(nu);
        if (el > BS_STAGE2_REPACE_RIDE * BS_STAGE2_WIN_BAND) _rp_ride++;
        if (viol > _rp_viol) _rp_viol = viol;
    }
    float repace_scale() const { return (float)_rp_scale; }
    uint16_t repace_events() const { return _rp_events; }
    uint16_t repace_refused() const { return _rp_refused; }
#endif
    bool build_failed() const { return _failed; }

    void invalidate() { _valid = false; _failed = false; _req = true; }

    // init-time memory reservation (see AP_BSolver::init): builds go into
    // this block instead of a per-upload calloc, which fragments out on
    // the flight processor
    void set_reserved(void *blk, size_t len) {
        _reserved = blk;
        _reserved_len = len;
    }

    // mission bookkeeping for the AUTO layer
    uint16_t run_first_idx() const { return _run_first; }
    uint16_t run_last_idx() const { return _run_last; }
    // schedule tick at which the given mission item's vertex is passed
    // (-1 when the item is not part of the compiled run)
    int32_t tau_vertex(uint16_t wp_idx) const;
    float run_alt_D_m() const { return _run_alt_d; }

private:
    bool scan(class AP_Mission &mission);

    bs_mission_tables _mt;
    bs_mission_report _rep;
    // streaming state: the vertex plan and the ring live inside the
    // boot-time reservation (see AP_BSolver::init); W is a power of two
    bs_mission_plan *_plan;
    double *_r_path, *_r_psi, *_r_ang;
    int *_r_sfam, *_r_srot, *_r_soff;
    signed char *_r_fam;
#if BS_STAGE2_PACE_TAB
    double *_r_pace;                  // published pace per tick (ring)
#endif
#if BS_STAGE2_REPACE
    double _rp_scale;                 // current pace scale, [1, MAX]
    double _rp_nu_sum, _rp_viol;      // accumulators since the last render
    int _rp_n, _rp_ride;
    uint16_t _rp_events, _rp_refused;
    void refresh_clock();             // plan -> tables after a re-pace
#endif
    int _ring_w, _ring_mask;
    int32_t _cursor;
    int _next_leg;
    void *_block;
    size_t _block_len;
    void *_reserved;
    size_t _reserved_len;
    bool _valid;
    bool _failed;
    bool _req;
    uint32_t _stamp_ms;            // mission.last_change_time_ms of build
    float _stamp_speed;
    uint32_t _txt_ms;              // refusal statustext throttle
    uint8_t _fail_st;              // latched refusal (0 = none)
    uint8_t _fail_gate;
    float _fail_detail;

    // scan output
    double _vx[BS_MB_MAX_WP], _vy[BS_MB_MAX_WP];
    uint16_t _wp_idx[BS_MB_MAX_WP];   // mission index per polyline vertex
    int _n_wp;
    uint16_t _run_first, _run_last;
    float _run_alt_d;
    float _v_cap;
};

#endif // AP_BSOLVER_ENABLED
