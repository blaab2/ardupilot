/*
 * AP_MPCTurnPlanner — AUTO-mission scanner for 180-deg boustrophedon
 * headland turns flyable by the onboard turn MPC (AP_MPCSolver).
 *
 * Scans the mission's nav waypoints for four consecutive plain waypoints
 * A,B,C,D forming a row reversal:
 *   - corner angle at B (between AB and BC) and at C (between BC and CD)
 *     each within 90 +/- MPC_TURN_ANG degrees;
 *   - the same turning sign at B and C (one continuous 180-deg reversal);
 *   - row spacing |BC| within [MPC_TURN_DMIN, MPC_TURN_DMAX] (the solver's
 *     seed-family band);
 *   - AB and CD anti-parallel within 2*MPC_TURN_ANG degrees;
 *   - all four waypoints at the same altitude (|dAlt| <= 0.5 m);
 *   - conservative guards: no DO_/CONDITION_ command between B and D
 *     (skipping the corner waypoints must not skip pump/camera triggers)
 *     and no hold time (p1 > 0) on either corner.
 *
 * Each detection is stored as a ready-to-arm AP_MPCSolver::TurnRecord in
 * the canonical R5 frame: entry along +xi = unit(B - A), corner line at
 * xi = 0 through origin_ne = B, eta positive toward the return row along
 * (-sin_h, cos_h); mirror = true when the return row lies on the -eta side.
 * entry_wp_idx = B's mission index (arm when that leg starts),
 * exit_wp_idx = D's mission index (resume the mission there after the turn).
 *
 * Scans are cached on AP_Mission::last_change_time_ms(): update() is a
 * timestamp compare unless the mission changed since the last scan.
 */
#pragma once

#include "AP_MPCSolver.h"

#if AP_MPCSOLVER_ENABLED

class AP_Mission;

class AP_MPCTurnPlanner {
public:
    AP_MPCTurnPlanner() {}

    CLASS_NO_COPY(AP_MPCTurnPlanner);

    // fixed table of detected turns
    static const uint8_t MAX_TURNS = 8;

    // rescan the mission if it changed since the last scan (otherwise a
    // cheap timestamp compare). Returns true when a scan actually ran
    // (i.e. the detection table is fresh — the caller logs DETECT events).
    bool update(AP_Mission &mission, const AP_MPCSolver &solver);

    // force the next update() to rescan
    void invalidate() { _valid = false; }

    // detected turn whose entry corner (B) has this mission index, or
    // nullptr if the leg ending at wp_index is not a detected turn entry
    const AP_MPCSolver::TurnRecord *find_entry(uint16_t wp_index) const;

    // detection table access (DETECT logging)
    uint8_t num_turns() const { return _num_turns; }
    const AP_MPCSolver::TurnRecord &turn(uint8_t i) const { return _turns[i]; }

private:
    void scan(AP_Mission &mission, const AP_MPCSolver &solver);

    AP_MPCSolver::TurnRecord _turns[MAX_TURNS];
    uint8_t _num_turns;
    uint32_t _scan_change_ms;   // mission.last_change_time_ms() of the scan
    bool _valid;
};

#endif  // AP_MPCSOLVER_ENABLED
