/* bs_arc_step — THE published-profile integration step, shared verbatim
 * by the batch builder's two passes and the streaming renderer.  One
 * body, one expression order: the stream==batch byte-equality gate
 * depends on both paths computing bit-identical doubles, so the step
 * lives here and nowhere else.
 *
 * Semantics (unchanged from the flown builder): jerk-limited ramp
 * toward the (drag-aware when enabled) forward acceleration cap,
 * clipped at trim, braked so the leg can still reach the target pace
 * vt at its end with constant BS_ARC_A_EFF authority.
 */
#ifndef BS_ARC_STEP_H
#define BS_ARC_STEP_H

#include <math.h>

/* the caller's translation unit defines arc_a_fwd() and the BS_ARC_*
 * constants before including this header */
static inline void bs_arc_step(double *v, double *a, double *arc,
                               double L, double vt, double v_trim)
{
    double an = *a + BS_ARC_J_EFF * BS_TS;
    {
        const double cap = arc_a_fwd(*v);
        if (an > cap) an = cap;
    }
    double vn = *v + 0.5 * (*a + an) * BS_TS;
    if (vn > v_trim) { vn = v_trim; an = 0.0; }
    const double rem_b = L - *arc;
    const double rem = (rem_b > 0.0) ? rem_b : 0.0;
    const double vbrk = sqrt(vt * vt + 2.0 * BS_ARC_A_EFF * rem);
    if (vn > vbrk) { vn = vbrk; an = 0.0; }
    *arc += 0.5 * (*v + vn) * BS_TS;
    *v = vn;
    *a = an;
}

#endif /* BS_ARC_STEP_H */
