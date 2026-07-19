/* Small dense kernels for the condensing layer (CUSTOM-SOLVER-PLAN.md
 * section 1). Portable C99, no BLAS, no allocation; every matrix is dense
 * COLUMN-MAJOR (the CasADi/cs_model convention). Sizes are tiny (nx=6,
 * nz<=~75) — plain triple loops beat any tiling here and keep the flash
 * cost near zero (plan section 6). */
#ifndef CS_LINALG_H
#define CS_LINALG_H

#include "cs_types.h"

void cs_zero(cs_real *x, int n);
void cs_copy(const cs_real *src, cs_real *dst, int n);

/* C = A*B      A: m x k, B: k x n, C: m x n (all col-major) */
void cs_gemm(int m, int n, int k,
             const cs_real *A, const cs_real *B, cs_real *C);

/* C += A*B */
void cs_gemm_acc(int m, int n, int k,
                 const cs_real *A, const cs_real *B, cs_real *C);

/* y += alpha * A' * x    A: m x n col-major, x: m, y: n */
void cs_gemv_t_acc(int m, int n, cs_real alpha,
                   const cs_real *A, const cs_real *x, cs_real *y);

/* C += alpha * A' * A    A: k x n col-major, C: n x n block with leading
 * dimension ldc >= n, col-major (full, both triangles written — DAQP reads
 * the upper triangle of dense H) */
void cs_syrk_acc(int n, int k, cs_real alpha, const cs_real *A, cs_real *C,
                 int ldc);

/* row-weighted variants: C += alpha * A' diag(w) A and
 * y += alpha * A' diag(w) x, w: k (resp. m) — the LM metric weights over the
 * SCALED state rows of E (cs_condense M2 staged state scaling). */
void cs_syrk_w_acc(int n, int k, cs_real alpha, const cs_real *w,
                   const cs_real *A, cs_real *C, int ldc);
void cs_gemv_tw_acc(int m, int n, cs_real alpha, const cs_real *w,
                    const cs_real *A, const cs_real *x, cs_real *y);

/* row' = a' * B          a: m, B: m x n col-major, row: n  (row = B'a) */
void cs_row_gemv(int m, int n, const cs_real *a, const cs_real *B,
                 cs_real *row);

cs_real cs_norm_inf(const cs_real *x, int n);

#endif /* CS_LINALG_H */
