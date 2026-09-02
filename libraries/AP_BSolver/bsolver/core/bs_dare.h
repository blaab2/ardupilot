/* bs_dare — small dense discrete algebraic Riccati solver for the runtime
 * mission builder.  Fixed-point iteration
 *
 *     P <- Q + A'PA - A'PB (R + B'PB)^-1 B'PA
 *
 * with symmetrization each step.  Sized for the terminal pair of record
 * (n <= 3 states, m <= 2 inputs; wp5_anytime_sim.py::dare_pair): the
 * closed loops contract at rho ~ 0.7-0.9, so a few hundred iterations
 * converge far below the acceptance tolerance.  Runs at mission-build
 * time on the solver thread, never in the per-tick path.
 */
#ifndef BS_DARE_H
#define BS_DARE_H

/* Solve the DARE for P (n x n, row-major) and the LQR gain
 * K = -(R + B'PB)^-1 B'PA (m x n).  A: n x n, B: n x m, Q: n x n,
 * R: m x m.  Returns 0 on success; nonzero when the iteration fails to
 * converge to `tol` (relative step) within `max_iter`, when the inner
 * m x m system is not positive definite, or when the residual check
 * fails.  P is seeded from Q. */
int bs_dare_solve(const double *A, const double *B, const double *Q,
                  const double *R, int n, int m,
                  double *P, double *K, int max_iter, double tol);

#endif /* BS_DARE_H */
