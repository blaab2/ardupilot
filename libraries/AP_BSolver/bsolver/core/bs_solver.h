/* bsolver — barrier-Newton core for the certified anytime MPC.
 *
 * The object of record is the relaxed-barrier economic MPC the WP5 theorems
 * are stated over: horizon N with an M-of-N committed prefix, solved by
 * damped Newton with Armijo backtracking at 4 Hz.  There is no QP anywhere
 * in this file, and that is the point — see ../README.md section 0.
 *
 * Precision is double by design (README section 3): the certified objects are
 * f64, the record documents an f64 resolution floor and a condensed-Hessian
 * conditioning ceiling, and at 4 Hz the budget affords it.
 */
#ifndef BS_SOLVER_H
#define BS_SOLVER_H

#include <stddef.h>

/* The model tables are a BUILD-TIME choice, not a source edit.  The default
 * is the record's N = 40 header; the corner-online flight configuration is
 * the N = 30 drop-in replacement, selected by defining BS_DATA_HEADER on the
 * compiler command line (host/Makefile target `corner`).  The two headers
 * define the same symbols and cannot coexist in one translation unit, which
 * is why this is a substitution and not an addition. */
#ifndef BS_DATA_HEADER
#define BS_DATA_HEADER "../model/generated/bs_model_data.h"
#endif
#include BS_DATA_HEADER

typedef double bs_real;

/* Output-tile edge for the blocked Hessian accumulation.
 *
 * MEASURED ON TARGET 2026-08-21, and the result refutes the hypothesis the
 * tiling was built for.  The B3 diagnosis held that the accumulation was
 * memory-bound — a 115 kB Hessian swept 50 times per stage against a 16 kB
 * D-cache — so tiling should have paid.  It does not.  Cube Orange+, engage
 * seed, 25 reps, q = 1 worst case:
 *
 *     rank-one (pre-tiling)      469.5 ms
 *     TILE =   8                 600.2 ms
 *     TILE =  16                 538.9 ms
 *     TILE =  32                 509.3 ms
 *     TILE = 128 (single tile)   475.1 ms   (min 459.9 vs 459.7 — identical)
 *
 * Monotone: every halving of the tile costs more, and the degenerate single
 * tile recovers the original exactly.  So the restructuring itself is free
 * and the tiling is pure overhead — panel re-reads and loop control buying
 * nothing, because the accumulation was never cache-bound.  The 400x
 * host:target ratio is therefore arithmetic, not stall: scalar double on the
 * FPv5-D16 against SIMD double on the host.
 *
 * The default is set to a no-op tile so the target keeps the faster path.
 * The parameter is retained because it is exact at any value (B1 and B2 are
 * bit-identical across the sweep) and another target with a different cache
 * may want it.
 */
#ifndef BS_HESS_TILE
#define BS_HESS_TILE 128
#endif

typedef enum {
    BS_OK = 0,
    BS_ERR_PHASE = 1,      /* phase index out of range                     */
    BS_ERR_DOMAIN = 2,     /* a row slack left the barrier domain          */
    BS_ERR_FACTOR = 3,     /* Hessian not positive definite to working tol */
    BS_ERR_ARENA = 4       /* workspace too small                          */
} bs_status;

/* A clock: family, chart rotation and affine seam offset per tick.  The
 * 59-phase instance wraps (it is periodic); the SN77 mission clamps at its
 * ends, exactly as Sched / SchedAff in wp5_anytime_sim.py and
 * model/sim_corner_online.py do.
 *
 * `offset` is the per-tick index into bs_off, the affine part of the seam
 * map.  The corner-online reference pace is piecewise constant across a
 * window, so the pace jumps at four ticks per window and the transition on
 * the phase quotient becomes AFFINE,
 *
 *     xi_{t+1} = T(t+1) (A xi_t + B u_t) + c(t+1),
 *
 * the offset applied AFTER the rotation and indexed by the DESTINATION tick,
 * clamped not wrapped -- the same indexing `rotation` already uses.  A
 * schedule with no affine part sets `offset` to NULL, which off_of() reads as
 * the zero vector; that is the record's plain linear step. */
typedef struct {
    const int *family;
    const int *rotation;
    const int *offset;     /* NULL == no affine part, i.e. all-zero */
    int length;
    int periodic;
} bs_schedule;

#ifdef BS_N_PHASE
/* The 59-phase periodic parity clock is an "a"/"t" object and exists only in
 * the record's data header; the corner-online header retires family "a" and
 * does not emit it. */
extern const bs_schedule bs_sched_periodic;
#endif
extern const bs_schedule bs_sched_mission;

/* Everything that depends only on the phase, built once per solve. */
typedef struct {
    const bs_schedule *schedule;
    int phase;
    /* Quadratic-half restriction: Hquad/Cquad rows below NU*npin_built were
     * NOT accumulated (see bs_problem_init_pinned).  An evaluation at
     * npin < npin_built is refused. */
    int npin_built;
    /* Abar is NOT stored: its block row is consumed in the same loop
     * iteration that produces it (the Cquad accumulation), and every other
     * consumer now takes predicted states from the plant recursion. */
    bs_real *gam;          /* BS_NX x BS_NV, the CURRENT block row */
    bs_real *gam_next;     /* BS_NX x BS_NV, ping-pong partner     */
    bs_real *Hquad;        /* BS_NV x BS_NV, the quadratic half */
    bs_real *Cquad;        /* BS_NV x BS_NX                     */
    /* Affine drift of the quadratic half: Dquad = sum_t 2 Gam_t' W_t Cbar_t,
     * the xi- and U-INDEPENDENT part of the gradient that the affine seam
     * offsets create.  Zero for a schedule without offsets, so the record
     * build carries it as a dead NV-vector of zeros.  See
     * bs_problem_init_pinned. */
    bs_real *Dquad;        /* BS_NV                             */
    bs_real *scratch;      /* >= BS_NROW*BS_NV + BS_NV*BS_NV    */
    /* Newton workspace.  These were file-scope statics, which put a third of
     * the core's RAM outside the caller's block and broke the one-arena
     * contract this header claims; on the target it also overflowed the
     * linker's RAM region.  They live in the arena now, so
     * bs_workspace_size() reports the core's whole footprint. */
    bs_real *nt_hess;      /* BS_NV x BS_NV */
    bs_real *nt_grad;      /* BS_NV */
    bs_real *nt_dir;       /* BS_NV */
    bs_real *nt_trial;     /* BS_NV */
} bs_problem;

/* Workspace sizing, in bs_real units. */
size_t bs_workspace_size(void);

/* sizeof(bs_problem), for callers that allocate it opaquely (the Python
 * harnesses do).  Exported because guessing it is a silent buffer overflow:
 * the struct grew from 64 to 88 bytes when the Newton workspace moved in, and
 * the harnesses kept allocating 64 and kept passing. */
size_t bs_problem_size(void);

/* Build the phase-dependent objects into `memory` (>= bs_workspace_size()). */
bs_status bs_problem_init(bs_problem *problem, const bs_schedule *schedule,
                          int phase, bs_real *memory, size_t memory_len);

/* Same, but accumulate the quadratic half (Hquad, Cquad) ONLY on the rows the
 * pinned path can read, a >= NU*npin.  Exact for evaluations at npin' >= npin
 * -- the objective value never touches Hquad/Cquad, and the gradient and
 * Hessian visit only those rows -- and refused for npin' < npin.  This is the
 * setup half of the interleaved cadence: bs_problem_init does not shrink when
 * the SOLVE does, and on target it is a third of the pinned cycle.
 * npin = 0 is bs_problem_init(), bit for bit. */
bs_status bs_problem_init_pinned(bs_problem *problem,
                                 const bs_schedule *schedule, int phase,
                                 int npin, bs_real *memory, size_t memory_len);

/* Objective, gradient and Hessian at (U, xi).  grad may be NULL; hess may be
 * NULL.  Mirrors MixedHorizon.gh() in the certificates term for term. */
bs_status bs_eval(const bs_problem *problem, const bs_real *U,
                  const bs_real *xi, bs_real *value,
                  bs_real *grad, bs_real *hess);

/* Same, with the first `npin` input stages PINNED as data.  `value` is still
 * the FULL objective (the reference's guards and Armijo test are stated on
 * it); `grad` and `hess` are the TAIL BLOCK of the full-horizon objects,
 * packed to nf = BS_NV - 3*npin (grad: nf entries; hess: nf x nf, row stride
 * nf).  This is exact, not a truncation: by causality a stage t < npin cannot
 * reach a column >= 3*npin, so the tail block is assembled term for term and
 * in the same order as the full one — bit-identical, gated in
 * tests/test_parity_pinned.py.  npin = 0 is bs_eval(). */
bs_status bs_eval_pinned(const bs_problem *problem, const bs_real *U,
                         const bs_real *xi, int npin, bs_real *value,
                         bs_real *grad, bs_real *hess);

typedef struct {
    int backtracks;     /* Armijo halvings taken                     */
    int armijo_fail;    /* line searches that never satisfied Armijo */
    int factor_fail;    /* Cholesky failures                         */
    int no_step;        /* guarded exits: no descent, or below noise */
} bs_newton_stats;

/* q damped-Newton steps in place.  Reports the first Newton decrement
 * lambda_0 = sqrt(g' H^-1 g) / sqrt(eps) — the quantity the certificates
 * bound at 59/59 arrivals — and the number of Armijo backtracks taken. */
bs_status bs_newton(const bs_problem *problem, bs_real *U, const bs_real *xi,
                    int q, bs_real *lambda0, bs_newton_stats *stats);

/* q damped-Newton steps on the FREE TAIL only: the first `npin` input stages
 * are committed data, replayed and not revised, and the Newton system is the
 * tail block H[3*npin:, 3*npin:] d_f = -g[3*npin:] of the full-horizon
 * objects.  This is the CONTROLLER OF RECORD's cadence — solve every tick
 * with the M-input committed prefix pinned (run_mission_interleaved /
 * MixedHorizon.newton_q_pinned in wp5_anytime_sim.py) — and the object the
 * interleaved certificate wp5_h1_interleaved_q.md is stated over.
 *
 * Two reasons it exists rather than bs_newton with a post-hoc overwrite of
 * the prefix.  (1) Correctness: premise A4 is the exact committed-prefix
 * predictor, so a solve that revises committed inputs is a different
 * controller.  (2) Cost: the tail is the horizon-(N - npin) problem, and the
 * barrier-Hessian work scales as sum_{t>=npin} [3(t+1) - 3*npin]^2 over
 * sum_t [3(t+1)]^2 = 0.427 at npin = 10, with the Cholesky at (90/120)^3.
 *
 * lambda0 reports sqrt(g_f' H_ff^-1 g_f)/sqrt(eps) at the first iteration —
 * the tail decrement, which is what the interleaved certificate bounds.
 * npin = 0 is bs_newton(). */
bs_status bs_newton_pinned(const bs_problem *problem, bs_real *U,
                           const bs_real *xi, int q, int npin,
                           bs_real *lambda0, bs_newton_stats *stats);

/* Terminal state of the plan: the last predicted state, from the plant
 * recursion. */
void bs_terminal_state(const bs_problem *problem, const bs_real *U,
                       const bs_real *xi, bs_real *xi_terminal);

/* Warm start for the next solve: shift by one and append the terminal-law
 * control of the family at the OLD horizon tip,
 *     U_shift = [ U[3:] , K[fam(phase + N)] xi_N ].
 * The append family is an identity, not a convention — it is gated upstream
 * (G6/G9 in the certificates) and must not be "tidied" to the new phase. */
void bs_shift_append(const bs_problem *problem, const bs_real *U,
                     const bs_real *xi, bs_real *U_shifted);

/* Dense Cholesky on a symmetric positive-definite n x n matrix, in place
 * (lower triangle), then solve.  Exposed for the parity harness. */
bs_status bs_chol_factor(bs_real *A, int n);
void bs_chol_solve(const bs_real *L, int n, bs_real *b);

#endif /* BS_SOLVER_H */
