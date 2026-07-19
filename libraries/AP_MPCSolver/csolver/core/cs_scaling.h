/* Physical-range diagonal scaling (CUSTOM-SOLVER-PLAN.md section 4.4).
 *
 * The dense QP is formed in the SCALED variable
 *     zeta = S^-1 * [dU_0 .. dU_{N-1}, ds_T, slacks]
 * with S = diag(scale_u per control column, 1 for s_T, 1 for slacks):
 *   - controls are scaled by their box ranges (CS_PD_SCALE_U), applied
 *     INSIDE condensing when the B / Ju columns are written, so the
 *     generated model functions keep physical units (no regeneration churn);
 *   - the free end time enters as s_T = T / T_ref (T_ref = 8, the documented
 *     cold-start scale) — the s_T column of the sensitivities is bT * T_ref,
 *     which is the plan's "s_T rather than raw T" fix for the worst-offender
 *     T column;
 *   - terminal slack columns stay in physical state units (they are ~0 at
 *     any healthy iterate; weights are quoted in physical units).
 *
 * State scales (CS_PD_SCALE_X) are PROVIDED here for the M2 equilibration /
 * kappa instrumentation; the M1 f64 build keeps states physical because
 *   (a) the LM damping replicates the acados twin's physical-space
 *       lambda*||dx||^2 exactly (parity first), and
 *   (b) DAQP row-normalizes every constraint row internally (normalize_M /
 *       normalize_Rinv), so external row equilibration is redundant for the
 *       QP — column scaling is what shapes kappa(H). Instrumented, not
 *       speculative: cs_solver reports a kappa proxy each iteration. */
#ifndef CS_SCALING_H
#define CS_SCALING_H

#include "cs_types.h"

typedef struct {
    cs_real x[6];    /* physical ranges of the states (M2 equilibration) */
    cs_real u[2];    /* control-column scales, applied in condensing     */
    cs_real T_ref;   /* s_T = T / T_ref                                  */
} cs_scales;

/* Fill from the generated problem data (turn_r5_probdata.h). */
void cs_scaling_default(cs_scales *s);

#endif /* CS_SCALING_H */
