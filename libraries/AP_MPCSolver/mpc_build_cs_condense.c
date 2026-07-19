/* build wrapper: compiles the vendored csolver core TU with the f32 flight
 * config (mpc_cs_defs.h). One wrapper per source file — the vendored tree
 * itself is not globbed by waf, and per-TU wrapping keeps static symbols
 * (e.g. CasADi helpers) from colliding. DO NOT add code here. */
#include "mpc_cs_defs.h"
#include "csolver/core/cs_condense.c"
