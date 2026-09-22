/* bs_dare — small dense discrete algebraic Riccati solver of the bsolver
 * runtime.  Fixed-point iteration
 *
 *     P <- Q + A'PA - A'PB (R + B'PB)^-1 B'PA
 *
 * seeded from Q and symmetrized each step (bs_dare.c), sized for the
 * terminal pair of record (n <= 3 states, m <= 2 inputs: the
 * longitudinal and the lateral DARE of dare_pair, wp5_anytime_sim.py /
 * arc_v2.py).
 *
 * WHERE IT RUNS (the contract since model revision v2c, 2026-09-13):
 *   - the arc-cell clock runtime solves LAZILY IN THE TICK PATH.
 *     bs_arc.c terminal_pair() needs the pair of the schedule's tip cell
 *     at the shift-append and at every render; the pair is keyed per
 *     cell (bs_arc_fam_key: 2 + k for every cell under cell_rows) and
 *     kept in a ring of BS_ARC_DARE_CACHE slots in bs_arc_sched, so a
 *     solve happens on a cache miss -- while the tip walks the chart,
 *     typically once per new cell -- with max_iter 20000, tol 1e-14
 *     (bs_arc.c arc_dare).  A nonzero return here is a FAILED TICK with
 *     no fallback: bs_arc_render returns -1 (the engage path,
 *     bs_arc_engage -> -1) and bs_arc_tick returns -2; the firmware
 *     (AP_BSolver.cpp) then reports the code and deactivates
 *     (_active = false).  bs_arc_sched counts the misses (dare_solves)
 *     and the failures (dare_fail).
 *   - the hover family of the v2d stop (cell kind BS_ARC_HV) is NEVER
 *     solved on target: its pair is mission-independent and is emitted
 *     by model/gen_runtime_const.py --v2 as bs_rc_arc_Ph / bs_rc_arc_Kh
 *     (arc_v2.dare_pair(arc_v2.hover_rows())); terminal_pair returns
 *     the emitted pair for that cell.
 *   - the v1 mission builder (bs_mission_build.c, stubbed under
 *     BS_NU == 2) solved once per family at mission-build time on the
 *     solver thread with max_iter 5000, tol 1e-13; that is where the
 *     old "never in the per-tick path" wording came from.
 *
 * ITERATION FACTS (bs_dare.c): P is seeded from Q; each step is
 * symmetrized, 0.5 (Pn + Pn'); convergence is the max-abs step
 * <= tol * max(1, max|P|) (a relative step) within max_iter; at the
 * converged P the gain K = -(R + B'PB)^-1 B'PA is formed and the
 * residual || Q + A'PA - A'PB G - P ||_max must be <= 1e-8 * max(1,
 * max|P|); finally every diagonal entry of P must be positive (a cheap
 * necessary check for P > 0; the full Cholesky is the host gate's).
 * The moving cells' closed loops contract at rho ~ 0.7-0.9 and converge
 * in a few hundred iterations; the hover pair (rho 0.991) converges
 * with a 13x iteration margin (doc/V2D-P0-2026-09-22.md section 3.2).
 */
#ifndef BS_DARE_H
#define BS_DARE_H

/* Solve the DARE for P (n x n, row-major) and the LQR gain
 * K = -(R + B'PB)^-1 B'PA (m x n).  A: n x n, B: n x m, Q: n x n,
 * R: m x m, all row-major.  P is seeded from Q.  Returns 0 on success,
 * else one of the codes bs_dare.c implements:
 *   1  dimensions out of range (n not in 1..3 or m not in 1..2);
 *   2  the inner m x m system R + B'PB is not positive definite (at an
 *      iteration step or when forming the final gain);
 *   3  no convergence to `tol` (relative max-abs step) within
 *      `max_iter`;
 *   4  the residual acceptance failed, || Q + A'PA - A'PB G - P ||_max
 *      > 1e-8 * max(1, max|P|) at the converged P;
 *   5  a diagonal entry of P is not positive.
 * On a nonzero return the outputs must not be used: P holds the last
 * iterate once the iteration has started (codes 2-5), K is written only
 * after the gain is formed (codes 4 and 5). */
int bs_dare_solve(const double *A, const double *B, const double *Q,
                  const double *R, int n, int m,
                  double *P, double *K, int max_iter, double tol);

#endif /* BS_DARE_H */
