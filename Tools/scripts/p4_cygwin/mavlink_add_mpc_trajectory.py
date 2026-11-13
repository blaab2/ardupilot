#!/usr/bin/env python3
"""Insert the MPC_TRAJECTORY (id 11061) message into an upstream
modules/mavlink checkout.

This fork's modules/mavlink gitlink (0dafe8cab2) is a local commit that
adds exactly one message to message_definitions/v1.0/ardupilotmega.xml on
top of upstream ArduPilot/mavlink 13f2f7351a.  That commit is not on
GitHub, so CI checks out the upstream parent and re-applies the message
from mavlink_mpc_trajectory.xml (the verbatim block, extracted from the
local commit).  The insertion goes immediately before the closing
</messages>, which is where the local commit put it, so the resulting
file is byte-identical to the local one (modulo line endings).

Usage: mavlink_add_mpc_trajectory.py <path/to/ardupilotmega.xml>
Idempotent: a file that already defines MPC_TRAJECTORY is left alone.
"""
import os
import sys


def main():
    xml_path = sys.argv[1]
    frag_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "mavlink_mpc_trajectory.xml")
    with open(xml_path, "rb") as fh:
        data = fh.read()
    with open(frag_path, "rb") as fh:
        frag = fh.read().replace(b"\r\n", b"\n")
    if b'name="MPC_TRAJECTORY"' in data:
        print("MPC_TRAJECTORY already present in %s" % xml_path)
        return 0
    crlf = b"\r\n" in data
    data = data.replace(b"\r\n", b"\n")
    tail = b"  </messages>\n</mavlink>"
    idx = data.rfind(tail)
    if idx < 0:
        print("closing </messages></mavlink> not found in %s" % xml_path)
        return 1
    data = data[:idx] + frag + data[idx:]
    if crlf:
        data = data.replace(b"\n", b"\r\n")
    with open(xml_path, "wb") as fh:
        fh.write(data)
    print("inserted MPC_TRAJECTORY (id 11061) into %s" % xml_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
