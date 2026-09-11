/* build wrapper for the vendored arc-cell clock layer (see
 * bs_build_solver.c).  Compiles to nothing unless BS_ARC=1. */
#include "bs_defs.h"

#if AP_BSOLVER_ENABLED && defined(BS_ARC) && BS_ARC
#include "bsolver/core/bs_arc.c"
#endif
