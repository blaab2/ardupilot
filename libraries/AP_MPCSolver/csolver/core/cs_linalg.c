#include "cs_linalg.h"

void cs_zero(cs_real *x, int n)
{
    int i;
    for (i = 0; i < n; ++i)
        x[i] = (cs_real)0;
}

void cs_copy(const cs_real *src, cs_real *dst, int n)
{
    int i;
    for (i = 0; i < n; ++i)
        dst[i] = src[i];
}

void cs_gemm(int m, int n, int k,
             const cs_real *A, const cs_real *B, cs_real *C)
{
    cs_zero(C, m * n);
    cs_gemm_acc(m, n, k, A, B, C);
}

/* All dot products below accumulate in cs_acc (double) and round to cs_real
 * only when stored — the plan-4.4 double-accumulation discipline. In the f64
 * build this is a no-op; in f32 it removes the O(sqrt(k))*eps_f32 summation
 * noise from H, the condensed rows, and the E/e propagation. */

void cs_gemm_acc(int m, int n, int k,
                 const cs_real *A, const cs_real *B, cs_real *C)
{
    int i, j, l;
    for (j = 0; j < n; ++j) {
        const cs_real *bj = B + (size_t)j * k;
        cs_real *cj = C + (size_t)j * m;
        for (i = 0; i < m; ++i) {
            cs_acc sum = (cs_acc)cj[i];
            for (l = 0; l < k; ++l)
                sum += (cs_acc)A[(size_t)l * m + i] * (cs_acc)bj[l];
            cj[i] = (cs_real)sum;
        }
    }
}

void cs_gemv_t_acc(int m, int n, cs_real alpha,
                   const cs_real *A, const cs_real *x, cs_real *y)
{
    int i, j;
    for (j = 0; j < n; ++j) {
        const cs_real *aj = A + (size_t)j * m;
        cs_acc sum = (cs_acc)0;
        for (i = 0; i < m; ++i)
            sum += (cs_acc)aj[i] * (cs_acc)x[i];
        y[j] = (cs_real)((cs_acc)y[j] + (cs_acc)alpha * sum);
    }
}

void cs_syrk_w_acc(int n, int k, cs_real alpha, const cs_real *w,
                   const cs_real *A, cs_real *C, int ldc)
{
    int i, j, l;
    for (j = 0; j < n; ++j) {
        const cs_real *aj = A + (size_t)j * k;
        for (i = 0; i <= j; ++i) {
            const cs_real *ai = A + (size_t)i * k;
            cs_acc sum = (cs_acc)0;
            cs_real s;
            for (l = 0; l < k; ++l)
                sum += (cs_acc)w[l] * (cs_acc)ai[l] * (cs_acc)aj[l];
            s = (cs_real)((cs_acc)alpha * sum);
            C[(size_t)j * ldc + i] += s;
            if (i != j)
                C[(size_t)i * ldc + j] += s;
        }
    }
}

void cs_gemv_tw_acc(int m, int n, cs_real alpha, const cs_real *w,
                    const cs_real *A, const cs_real *x, cs_real *y)
{
    int i, j;
    for (j = 0; j < n; ++j) {
        const cs_real *aj = A + (size_t)j * m;
        cs_acc sum = (cs_acc)0;
        for (i = 0; i < m; ++i)
            sum += (cs_acc)w[i] * (cs_acc)aj[i] * (cs_acc)x[i];
        y[j] = (cs_real)((cs_acc)y[j] + (cs_acc)alpha * sum);
    }
}

void cs_syrk_acc(int n, int k, cs_real alpha, const cs_real *A, cs_real *C,
                 int ldc)
{
    int i, j, l;
    for (j = 0; j < n; ++j) {
        const cs_real *aj = A + (size_t)j * k;
        for (i = 0; i <= j; ++i) {
            const cs_real *ai = A + (size_t)i * k;
            cs_acc sum = (cs_acc)0;
            cs_real s;
            for (l = 0; l < k; ++l)
                sum += (cs_acc)ai[l] * (cs_acc)aj[l];
            s = (cs_real)((cs_acc)alpha * sum);
            C[(size_t)j * ldc + i] += s;
            if (i != j)
                C[(size_t)i * ldc + j] += s;
        }
    }
}

void cs_row_gemv(int m, int n, const cs_real *a, const cs_real *B,
                 cs_real *row)
{
    int i, j;
    for (j = 0; j < n; ++j) {
        const cs_real *bj = B + (size_t)j * m;
        cs_acc sum = (cs_acc)0;
        for (i = 0; i < m; ++i)
            sum += (cs_acc)a[i] * (cs_acc)bj[i];
        row[j] = (cs_real)sum;
    }
}

cs_real cs_norm_inf(const cs_real *x, int n)
{
    int i;
    cs_real m = (cs_real)0;
    for (i = 0; i < n; ++i) {
        cs_real v = x[i] < (cs_real)0 ? -x[i] : x[i];
        if (v > m)
            m = v;
    }
    return m;
}
