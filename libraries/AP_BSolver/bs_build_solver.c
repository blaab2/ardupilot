/* build wrapper: compiles the vendored bsolver core TU. One wrapper per
 * source file — the vendored tree itself is not globbed by waf. The core is
 * f64 by design (doc/bsolver.tex section 3: cond(H) reaches 5.2e8, so
 * single precision would leave the Newton direction with no correct digits),
 * so unlike the csolver wrapper there is no precision macro to set.
 * DO NOT add code here. */
#include "bs_defs.h"

#if AP_BSOLVER_ENABLED
#include "bsolver/core/bs_solver.c"
#endif
