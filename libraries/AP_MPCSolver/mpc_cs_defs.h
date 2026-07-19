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
