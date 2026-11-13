#!/usr/bin/env bash

# P4 variant of Tools/scripts/cygwin_build.sh: builds ONLY ArduCopter SITL
# for Mission Planner (Windows, Cygwin x86_64) with the arc-cell clock
# planner (AP_BSolver, BSLV_ARC=1 = revision v2e + BSLV_WDJ) and packages
#   ArduCopter.elf.exe  (rename to ArduCopter.exe for Mission Planner)
#   cyg*.dll            (every Cygwin DLL cygcheck lists for the exe)
#   SN77.json           (SITL frame model; also embedded in ROMFS)
#   p4_sitl.parm        (P4 SITL parameter set; also embedded in ROMFS)
#   git.txt, build.txt, romfs.txt
#
# The frame model and the parameter file are embedded in the binary's ROMFS
# (Tools/autotest/models/*.json|*.parm are embedded by boards.py at
# configure time), so the exe needs no files in its working directory:
#   --model +:SN77.json            -> SIM_Frame falls back to @ROMFS/models/SN77.json
#   --defaults @ROMFS/models/p4_sitl.parm
#
# Run from the root of the ardupilot tree.

set -x
set -e

TOOLCHAIN=x86_64-pc-cygwin
GPP_COMPILER="${TOOLCHAIN}-g++"
P4DIR=Tools/scripts/p4_cygwin

$GPP_COMPILER --version
SYS_ROOT=$($GPP_COMPILER -print-sysroot)
echo "SYS_ROOT=$SYS_ROOT"

rm -rf artifacts
mkdir artifacts

if [ -z "$GITHUB_ACTIONS" ] || [ "$GITHUB_ACTIONS" != "true" ]; then
WAF_OPTIONS="-j8"
fi

# ROMFS payload: the SN77 frame model (the SimOnHardWare hwdef already
# carries the reference copy) and the P4 SITL parameter set.  Copied, not
# committed, so the Linux SITL tree keeps its upstream ROMFS.
cp -v libraries/AP_HAL_ChibiOS/hwdef/CubeOrangePlus-SimOnHardWare/SN77.json Tools/autotest/models/SN77.json
cp -v $P4DIR/p4_sitl.parm Tools/autotest/models/p4_sitl.parm

# modules/mavlink is NOT at its gitlink (that commit is local-only; CI
# checks out its upstream parent + MPC_TRAJECTORY), so waf must not try to
# "update" it back: --no-submodule-update.
export BSLV_ARC=1
(
    python ./waf --color yes --toolchain $TOOLCHAIN --board sitl configure \
        --no-submodule-update \
        --define AP_BSOLVER_ENABLED=1 --define AP_MPCSOLVER_ENABLED=0 2>&1
    python ./waf copter $WAF_OPTIONS 2>&1
) | tee artifacts/build.txt

grep -q "AP_BSolver arc-cell clock runtime" artifacts/build.txt

cp -v build/sitl/bin/arducopter artifacts/ArduCopter.elf.exe

# the arc build and the W_DJ statustexts must be in the binary
grep -a -q "BSLV: arc wdj" artifacts/ArduCopter.elf.exe
grep -a -q "BSLV: engaged (q=%d K=%d wdj=%.2f)" artifacts/ArduCopter.elf.exe

for exe in artifacts/*.exe; do
    echo $exe
    cygcheck $exe | grep -oP 'cyg[^\s\\/]+\.dll' | sort -u | while read -r line; do
      cp -v /usr/bin/$line artifacts/
    done
done

cp -v Tools/autotest/models/SN77.json artifacts/SN77.json
cp -v $P4DIR/p4_sitl.parm artifacts/p4_sitl.parm

# what went into ROMFS (configure-time list)
grep -a -o "models/[A-Za-z0-9_.-]*" build/sitl/ap_romfs_embedded.h 2>/dev/null | sort -u > artifacts/romfs.txt || true
grep -q "models/SN77.json" artifacts/romfs.txt
grep -q "models/p4_sitl.parm" artifacts/romfs.txt

git log -1 > artifacts/git.txt
echo "BSLV_ARC=1 --define AP_BSOLVER_ENABLED=1 --define AP_MPCSOLVER_ENABLED=0 --no-submodule-update" >> artifacts/git.txt
(cd modules/mavlink && echo "modules/mavlink: $(git rev-parse HEAD) + MPC_TRAJECTORY (Tools/scripts/p4_cygwin/mavlink_mpc_trajectory.xml)") >> artifacts/git.txt || true
ls -l artifacts/
