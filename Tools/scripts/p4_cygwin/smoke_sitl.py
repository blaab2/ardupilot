#!/usr/bin/env python3
"""Smoke test of the P4 Windows SITL build (CI, windows runner).

Expects ArduCopter.elf.exe already running as Mission Planner would start
it with the P4 extra command line:
    --model +:SN77.json --defaults @ROMFS/models/p4_sitl.parm
    --home 50.2917,10.3986,0,0 -w
and MAVLink on tcp:127.0.0.1:5760.

Checks:
  1. BSLV_ENABLE/BSLV_Q/BSLV_K/BSLV_WDJ come up as 1/1/3/1.0 from the
     ROMFS parameter file (no parameter is set over MAVLink),
  2. the uploaded turns mission is built on board ('BSLV mission: arc'),
  3. after GUIDED takeoff + AUTO: 'BSLV: engaged (q=1 K=3 wdj=1.00)',
     'BSLV: arc wdj 1.000', 'BSLV: mission complete, N ticks, settled ...'
     and 'BSLV: run complete'.
Writes smoke_result.json and exits non-zero on any failed check.

Usage: smoke_sitl.py <mission.waypoints> [out_dir] [conn]
"""
import json
import math
import os
import re
import sys
import time

from pymavlink import mavutil

WP_FILE = sys.argv[1]
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "."
CONN = sys.argv[3] if len(sys.argv) > 3 else "tcp:127.0.0.1:5760"
TAKEOFF_ALT = 20.0
MODE_GUIDED = 4
MODE_AUTO = 3

result = {"checks": {}, "statustexts": [], "params": {}}
T0 = time.time()


def log(msg):
    print("[smoke %7.1f] %s" % (time.time() - T0, msg), flush=True)


def finish(code):
    ok = all(result["checks"].values()) and code == 0
    result["pass"] = bool(ok)
    with open(os.path.join(OUT_DIR, "smoke_result.json"), "w") as fh:
        json.dump(result, fh, indent=1)
    log("RESULT %s %s" % ("PASS" if ok else "FAIL",
                          json.dumps(result["checks"])))
    sys.exit(0 if ok else 1)


def connect():
    deadline = time.time() + 120
    while time.time() < deadline:
        try:
            m = mavutil.mavlink_connection(CONN, source_system=255)
            hb = m.wait_heartbeat(timeout=10)
            if hb is not None:
                return m
        except Exception as ex:  # connection refused while SITL boots
            log("connect: %s" % ex)
        time.sleep(3)
    return None


def note_text(msg):
    t = msg.text
    result["statustexts"].append([round(time.time() - T0, 2), t])
    log("statustext: " + t)


def pump(m, seconds, types=None):
    """read for `seconds`, logging statustexts; returns list of other msgs"""
    out = []
    end = time.time() + seconds
    while time.time() < end:
        msg = m.recv_match(blocking=True, timeout=0.5)
        if msg is None:
            continue
        if msg.get_type() == "STATUSTEXT":
            note_text(msg)
        elif types and msg.get_type() in types:
            out.append(msg)
    return out


def seen(pattern):
    rx = re.compile(pattern)
    for _, t in result["statustexts"]:
        mt = rx.search(t)
        if mt:
            return mt
    return None


def get_param(m, name):
    for _ in range(5):
        m.mav.param_request_read_send(m.target_system, m.target_component,
                                      name.encode(), -1)
        end = time.time() + 3
        while time.time() < end:
            msg = m.recv_match(blocking=True, timeout=0.5)
            if msg is None:
                continue
            if msg.get_type() == "STATUSTEXT":
                note_text(msg)
            elif msg.get_type() == "PARAM_VALUE" and msg.param_id == name:
                return msg.param_value
    return None


def read_waypoints(path):
    items = []
    with open(path) as fh:
        head = fh.readline()
        assert head.startswith("QGC WPL"), head
        for line in fh:
            f = line.split()
            if len(f) < 12:
                continue
            items.append(dict(seq=int(f[0]), current=int(f[1]),
                              frame=int(f[2]), cmd=int(f[3]),
                              p=[float(x) for x in f[4:8]],
                              lat=float(f[8]), lon=float(f[9]),
                              alt=float(f[10]), auto=int(f[11])))
    return items


def upload(m, items):
    m.mav.mission_count_send(m.target_system, m.target_component, len(items),
                             mavutil.mavlink.MAV_MISSION_TYPE_MISSION)
    sent = set()
    end = time.time() + 30
    ack = None
    while time.time() < end:
        msg = m.recv_match(type=["MISSION_REQUEST", "MISSION_REQUEST_INT",
                                 "MISSION_ACK", "STATUSTEXT"],
                           blocking=True, timeout=2)
        if msg is None:
            continue
        if msg.get_type() == "STATUSTEXT":
            note_text(msg)
            continue
        if msg.get_type() == "MISSION_ACK":
            ack = msg
            break
        it = items[msg.seq]
        m.mav.mission_item_int_send(
            m.target_system, m.target_component, it["seq"], it["frame"],
            it["cmd"], it["current"], it["auto"], *it["p"],
            int(round(it["lat"] * 1e7)), int(round(it["lon"] * 1e7)),
            it["alt"], mavutil.mavlink.MAV_MISSION_TYPE_MISSION)
        sent.add(msg.seq)
    ok = ack is not None and ack.type == 0 and len(sent) == len(items)
    log("mission upload %d/%d items, ack %s" % (len(sent), len(items), ack))
    return ok


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    m = connect()
    result["checks"]["heartbeat"] = m is not None
    if m is None:
        finish(1)
    log("heartbeat from system %d" % m.target_system)

    for name in ("BSLV_ENABLE", "BSLV_Q", "BSLV_K", "BSLV_WDJ", "BSLV_SEED",
                 "WP_SPD", "PSC_ANGLE_MAX", "SIM_BATT_VOLTAGE"):
        result["params"][name] = get_param(m, name)
    log("params %s" % result["params"])
    p = result["params"]
    result["checks"]["params_from_romfs"] = (
        p["BSLV_ENABLE"] == 1 and p["BSLV_Q"] == 1 and p["BSLV_K"] == 3
        and p["BSLV_WDJ"] is not None and abs(p["BSLV_WDJ"] - 1.0) < 1e-6)

    for mid, rate in ((32, 10), (30, 10)):
        m.mav.command_long_send(m.target_system, m.target_component,
                                mavutil.mavlink.MAV_CMD_SET_MESSAGE_INTERVAL,
                                0, mid, int(1e6 / rate), 0, 0, 0, 0, 0)

    log("waiting for the EKF to use GPS")
    end = time.time() + 120
    while time.time() < end and not seen(r"EKF3 IMU0 is using GPS"):
        pump(m, 1.0)
    pump(m, 5.0)

    items = read_waypoints(WP_FILE)
    result["checks"]["mission_upload"] = upload(m, items)
    pump(m, 3.0)
    result["checks"]["mission_built"] = seen(r"BSLV mission: arc") is not None

    m.mav.set_mode_send(m.target_system,
                        mavutil.mavlink.MAV_MODE_FLAG_CUSTOM_MODE_ENABLED,
                        MODE_GUIDED)
    pump(m, 1.5)
    armed = False
    for attempt in range(30):
        m.arducopter_arm()
        for hb in pump(m, 3.0, types=("HEARTBEAT",)):
            if hb.base_mode & mavutil.mavlink.MAV_MODE_FLAG_SAFETY_ARMED:
                armed = True
        if armed:
            break
        log("arm attempt %d not yet accepted" % (attempt + 1))
    result["checks"]["armed"] = armed
    if not armed:
        finish(1)
    m.mav.command_long_send(m.target_system, m.target_component,
                            mavutil.mavlink.MAV_CMD_NAV_TAKEOFF,
                            0, 0, 0, 0, 0, 0, 0, TAKEOFF_ALT)
    at_alt = False
    end = time.time() + 120
    while time.time() < end and not at_alt:
        for lp in pump(m, 1.0, types=("LOCAL_POSITION_NED",)):
            if -lp.z > TAKEOFF_ALT - 1.0:
                at_alt = True
    result["checks"]["takeoff"] = at_alt
    log("at altitude: %s" % at_alt)

    end = time.time() + 40
    while time.time() < end:
        m.mav.set_position_target_local_ned_send(
            0, m.target_system, m.target_component,
            mavutil.mavlink.MAV_FRAME_LOCAL_NED, 0b0000111111111000,
            0.0, 0.0, -TAKEOFF_ALT, 0, 0, 0, 0, 0, 0, 0, 0)
        lps = pump(m, 1.0, types=("LOCAL_POSITION_NED",))
        if lps and math.hypot(lps[-1].x, lps[-1].y) < 1.0 and \
                math.hypot(lps[-1].vx, lps[-1].vy) < 0.4:
            break

    m.mav.set_mode_send(m.target_system,
                        mavutil.mavlink.MAV_MODE_FLAG_CUSTOM_MODE_ENABLED,
                        MODE_AUTO)
    log("switched to AUTO")
    end = time.time() + 300
    while time.time() < end and not seen(r"BSLV: run complete"):
        pump(m, 1.0)
    pump(m, 3.0)

    eng = seen(r"BSLV: engaged \(q=(\d+) K=(\d+) wdj=([0-9.]+)\)")
    result["engaged"] = eng.group(0) if eng else None
    result["checks"]["engaged_q1_k3_wdj1"] = bool(
        eng and eng.group(1) == "1" and eng.group(2) == "3"
        and eng.group(3) == "1.00")
    result["checks"]["arc_wdj_1"] = seen(r"BSLV: arc wdj 1\.000") is not None
    done = seen(r"BSLV: mission complete, (\d+) ticks, settled t(\d+)")
    result["mission_complete"] = done.group(0) if done else None
    result["ticks"] = int(done.group(1)) if done else None
    result["checks"]["mission_complete_settled"] = done is not None
    result["checks"]["run_complete"] = seen(r"BSLV: run complete") is not None
    ovr = seen(r"BSLV: done (\d+) ticks off ([+-]?\d+) ovr (\d+)")
    result["done"] = ovr.group(0) if ovr else None
    finish(0)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as ex:
        result["exception"] = repr(ex)
        log("exception %r" % ex)
        result["checks"]["no_exception"] = False
        finish(1)
