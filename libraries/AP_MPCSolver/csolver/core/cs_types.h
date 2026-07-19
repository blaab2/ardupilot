/* csolver core types — portable C99, no deps, no malloc (CUSTOM-SOLVER-PLAN.md
 * sections 1, 2.1). cs_real is float on target (build with -DCS_REAL_SINGLE,
 * matching -Dcasadi_real=float for the generated model) and double for the
 * host reference twin. */
#ifndef CS_TYPES_H
#define CS_TYPES_H

#include <stddef.h>
#include <stdint.h>

#ifdef CS_REAL_SINGLE
typedef float cs_real;
#else
typedef double cs_real;
#endif

/* Accumulation type for residual/defect/step bookkeeping and the dense dot
 * products (CUSTOM-SOLVER-PLAN.md section 4.4 "mixed precision where it is
 * cheap"): the M7's FPU is double-capable (fpv5-d16), so dot products and
 * the feasibility-oracle quantities accumulate in double while ALL storage
 * and the QP solve stay cs_real. In the f64 build cs_acc == cs_real and
 * nothing changes. */
typedef double cs_acc;

/* problem dimensions (fixed structure; N is a per-build/model-file choice) */
typedef struct {
    int nx;   /* states   (r5: 6)  */
    int nu;   /* controls (r5: 2)  */
    int nh;   /* h rows   (r5: 5)  */
    int N;    /* shooting intervals */
} cs_dims;

/* status codes */
enum {
    CS_OK        = 0,
    CS_ERR_ARENA = 1,   /* arena exhausted / too small          */
    CS_ERR_INIT  = 2,   /* module used before successful init   */
    CS_ERR_BADN  = 3,   /* no generated model for this N        */
    CS_ERR_MODEL = 4,   /* generated function returned nonzero  */
    CS_ERR_ARG   = 5,   /* NULL pointer / bad argument          */
};

/* CS_NO_LIBC_ALLOC: the no-allocation contract, enforced at compile time.
 * Core code must never allocate — every buffer comes from the caller-provided
 * arena (cs_arena.h), fixed after init. The poison below makes any direct
 * use of the libc allocator in core translation units a compile error; the
 * host Makefile's `check-noalloc` target additionally asserts via `nm` that
 * no core object references the allocator symbols (belt and braces).
 * <stdlib.h>/<string.h> are included FIRST so their declarations pre-date the
 * poison; include cs_types.h before any further libc headers in core code. */
#if defined(CS_NO_LIBC_ALLOC) && defined(__GNUC__)
#include <stdlib.h>
#include <string.h>
#pragma GCC poison malloc calloc realloc free
#endif

#endif /* CS_TYPES_H */
