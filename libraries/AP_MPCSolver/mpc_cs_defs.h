/*
 * mpc_cs_defs.h — the f32 FLIGHT build configuration of the vendored csolver
 * core (csolver/, synced from scripts/mpc/csolver by sync_firmware.sh).
 *
 * MUST be included FIRST by every TU that compiles or calls the vendored
 * sources (the mpc_cs_*.c build wrappers and AP_MPCSolver.cpp), so cs_real,
 * casadi_real and the DAQP c_float all agree on float. This mirrors the
 * host Makefile's DEFS_F32 (scripts/mpc/csolver/host/Makefile):
 *   -DCS_REAL_SINGLE -DDAQP_SINGLE_PRECISION -Dcasadi_real=float
 *   -DCS_NO_LIBC_ALLOC (+ CS_PINNED_ONLY, the flight arena lever)
 */
#pragma once

#define CS_REAL_SINGLE 1
#define DAQP_SINGLE_PRECISION 1
/* flight problems are fully pinned — drops the free-boundary reserve from
 * the arena (csolver/README.md "CS_PINNED_ONLY") */
#define CS_PINNED_ONLY 1

/* Layer-2 chi-box symmetrization (csolver ABI 11 compile gate). Default 0 =
 * the flown [-0.2, pi] gauge box, bit-identical to the pre-ABI-11 behavior.
 * 1 widens the h2 box to [-pi, pi] (brake-and-crab quadrant) and compiles the
 * cs_rh entry-point branch discipline (branch_commit: nearest-branch wrap of
 * measured psi/theta onto the accepted plan + hysteresis at +-pi) PLUS the
 * firmware-side twins in AP_MPCSolver.cpp (build_x_meas chi clamp, READY-pin
 * branch wrap). Overridable from the command line for build checks
 * (waf configure --define CS_CHI_SYM=1); MUST be consistent across every TU
 * including this header. DO NOT fly =1 before the sym re-certification stage
 * (SOSC/multipliers/kappa/lambda-map on the sym formulation) has run —
 * Layer-2 roadmap, CUBE-BENCH-NOTES. */
#ifndef CS_CHI_SYM
#define CS_CHI_SYM 0
#endif

/* the generated CasADi code guards its scalar typedef with #ifndef */
#ifndef casadi_real
#define casadi_real float
#endif

#ifndef __cplusplus
/* Core no-allocation contract (poisons malloc/free via cs_types.h) — for the
 * vendored C TUs only. The C++ side (AP_MPCSolver.cpp) allocates the arena
 * once at boot, which is exactly the caller role the contract expects. */
#define CS_NO_LIBC_ALLOC 1

/* The vendored sources are pinned upstream code (DAQP) and generated code
 * (CasADi) — their style warnings are not ours to fix (same waivers as the
 * host Makefile's -Wno-* set, widened for ArduPilot's -Werror=* selections). */
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#endif
