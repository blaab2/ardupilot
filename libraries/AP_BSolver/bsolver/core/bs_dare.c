/* bs_dare — see bs_dare.h. */
#include "bs_dare.h"

#include <math.h>
#include <string.h>

#define NMAX 3
#define MMAX 2

/* c = a*b, row-major, (ra x ca) * (ca x cb) */
static void mul(const double *a, const double *b, double *c,
                int ra, int ca, int cb)
{
    for (int i = 0; i < ra; ++i) {
        for (int j = 0; j < cb; ++j) {
            double acc = 0.0;
            for (int k = 0; k < ca; ++k) {
                acc += a[i * ca + k] * b[k * cb + j];
            }
            c[i * cb + j] = acc;
        }
    }
}

/* c = a'*b, a is (ra x ca) so a' is (ca x ra); b is (ra x cb) */
static void mul_tn(const double *a, const double *b, double *c,
                   int ra, int ca, int cb)
{
    for (int i = 0; i < ca; ++i) {
        for (int j = 0; j < cb; ++j) {
            double acc = 0.0;
            for (int k = 0; k < ra; ++k) {
                acc += a[k * ca + i] * b[k * cb + j];
            }
            c[i * cb + j] = acc;
        }
    }
}

/* invert a symmetric positive-definite m x m (m <= 2) in place-of `inv`;
 * returns 0 on success. */
static int spd_inv(const double *S, double *inv, int m)
{
    if (m == 1) {
        if (!(S[0] > 0.0)) return 1;
        inv[0] = 1.0 / S[0];
        return 0;
    }
    /* m == 2 */
    const double a = S[0], b = S[1], c = S[2], d = S[3];
    const double det = a * d - b * c;
    if (!(a > 0.0) || !(det > 0.0)) return 1;
    inv[0] = d / det;
    inv[1] = -b / det;
    inv[2] = -c / det;
    inv[3] = a / det;
    return 0;
}

int bs_dare_solve(const double *A, const double *B, const double *Q,
                  const double *R, int n, int m,
                  double *P, double *K, int max_iter, double tol)
{
    if (n < 1 || n > NMAX || m < 1 || m > MMAX) return 1;

    double PA[NMAX * NMAX], PB[NMAX * MMAX];
    double BtPA[MMAX * NMAX], BtPB[MMAX * MMAX], Sinv[MMAX * MMAX];
    double G[MMAX * NMAX];                    /* (R+B'PB)^-1 B'PA */
    double AtPA[NMAX * NMAX], AtPB_G[NMAX * NMAX], Pn[NMAX * NMAX];
    double AtP[NMAX * NMAX];

    memcpy(P, Q, sizeof(double) * (size_t)(n * n));

    int converged = 0;
    for (int it = 0; it < max_iter; ++it) {
        mul(P, A, PA, n, n, n);
        mul(P, B, PB, n, n, m);
        mul_tn(B, PA, BtPA, n, m, n);         /* B'PA: m x n */
        mul_tn(B, PB, BtPB, n, m, m);         /* B'PB: m x m */
        for (int i = 0; i < m * m; ++i) BtPB[i] += R[i];
        if (spd_inv(BtPB, Sinv, m) != 0) return 2;
        mul(Sinv, BtPA, G, m, m, n);          /* G: m x n */

        mul_tn(A, PA, AtPA, n, n, n);         /* A'PA */
        /* A'PB * G = (B'PA)' G  (n x n) */
        mul_tn(BtPA, G, AtPB_G, m, n, n);

        double step = 0.0, scale = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                const double v = Q[i * n + j] + AtPA[i * n + j]
                               - AtPB_G[i * n + j];
                Pn[i * n + j] = v;
            }
        }
        /* symmetrize, measure the step */
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                const double v = 0.5 * (Pn[i * n + j] + Pn[j * n + i]);
                const double d = fabs(v - P[i * n + j]);
                if (d > step) step = d;
                const double av = fabs(v);
                if (av > scale) scale = av;
                P[i * n + j] = v;
            }
        }
        if (step <= tol * (scale > 1.0 ? scale : 1.0)) {
            converged = 1;
            break;
        }
    }
    if (!converged) return 3;

    /* K = -(R + B'PB)^-1 B'PA at the converged P */
    mul(P, A, PA, n, n, n);
    mul(P, B, PB, n, n, m);
    mul_tn(B, PA, BtPA, n, m, n);
    mul_tn(B, PB, BtPB, n, m, m);
    for (int i = 0; i < m * m; ++i) BtPB[i] += R[i];
    if (spd_inv(BtPB, Sinv, m) != 0) return 2;
    mul(Sinv, BtPA, G, m, m, n);
    for (int i = 0; i < m * n; ++i) K[i] = -G[i];

    /* residual acceptance: || Q + A'PA - A'PB G - P ||_max small */
    mul_tn(A, PA, AtPA, n, n, n);
    mul_tn(BtPA, G, AtPB_G, m, n, n);
    (void)AtP;
    double res = 0.0, scale = 0.0;
    for (int i = 0; i < n * n; ++i) {
        const double v = Q[i] + AtPA[i] - AtPB_G[i] - P[i];
        if (fabs(v) > res) res = fabs(v);
        if (fabs(P[i]) > scale) scale = fabs(P[i]);
    }
    if (res > 1e-8 * (scale > 1.0 ? scale : 1.0)) return 4;
    /* P must be positive definite: diagonal positivity is a cheap
     * necessary check; the builder's gate does the full Cholesky. */
    for (int i = 0; i < n; ++i) {
        if (!(P[i * n + i] > 0.0)) return 5;
    }
    return 0;
}
