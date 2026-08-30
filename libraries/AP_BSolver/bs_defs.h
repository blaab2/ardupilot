/* bs_defs.h — C-safe build configuration shared by the vendored core wrapper
 * and AP_BSolver.  Kept separate from AP_BSolver.h because the wrapper TU is
 * plain C.
 *
 * The core is f64 unconditionally: doc/bsolver.tex section 3 measures
 * cond(H) up to 5.2e8, at which single precision leaves the Newton direction
 * with no correct digits.  There is therefore no precision lever here, unlike
 * the csolver side's mpc_cs_defs.h. */
#pragma once

#ifndef AP_BSOLVER_ENABLED
#define AP_BSOLVER_ENABLED 0
#endif

/* Compile-gated on-target WCET bench (B3).  Default OFF; never for flight.
 * Lives here rather than in the .cpp so AP_BSolver.h sees it when it guards
 * the member declaration. */
#ifndef BS_BENCH_WCET
#define BS_BENCH_WCET 0
#endif

/* Static (.bss) footprint of the core, measured with `size` on the
 * cross-compiled object.  It is now ZERO: the Newton workspace that used to
 * sit here (hess[120*120] plus three [120] vectors, 118,080 B) moved into the
 * caller-provided arena, so bs_workspace_size() reports the core's WHOLE
 * footprint and the one-arena contract holds.  Kept as a named constant so
 * the bench keeps asserting it rather than assuming it.
 *
 * The previous value 118080u was left frozen here after that move and the
 * bench printed it as if live — a stale-constant reporting defect.  Verified
 * against `arm-none-eabi-size` on the cross-compiled core: bss = 0. */
#define BS_STATIC_BSS 0u

/* ------------------------------------------------- driver buffer dimensions
 * AP_BSolver.h declares the decision vector and the committed prefix as
 * fixed-size members, and it must NOT include the 94 kB generated model
 * header to learn their length -- every TU that includes AP_BSolver.h would
 * pay for it.  So the two dimensions are MIRRORED here and CHECKED in
 * AP_BSolver.cpp, which does see the real header:
 *
 *     static_assert(BS_DRV_NV == BS_NV, ...)
 *     static_assert(BS_DRV_M  == BS_M_COMMIT, ...)
 *
 * A stale mirror is therefore a compile error, not a silent buffer overrun.
 * The default is the RECORD configuration (N = 40 -> NV = 120); the
 * corner-online flight build overrides BS_DRV_NV to 90 on the command line,
 * next to BS_DATA_HEADER, in wscript.  M is 10 in both configurations, so it
 * has no override -- only the check. */
#ifndef BS_DRV_NV
#define BS_DRV_NV 120
#endif
#ifndef BS_DRV_M
#define BS_DRV_M 10
#endif
