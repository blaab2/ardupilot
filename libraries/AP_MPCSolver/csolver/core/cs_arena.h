/* Bump allocator over ONE caller-provided memory block (CUSTOM-SOLVER-PLAN.md
 * section 2.1: static allocation, no cycle-time malloc). Every pointer is
 * carved at init and fixed afterwards; per-cycle work touches only the arena
 * and the stack. Exhaustion returns NULL (caller maps it to CS_ERR_ARENA) —
 * there is no fallback allocation path by design. */
#ifndef CS_ARENA_H
#define CS_ARENA_H

#include "cs_types.h"

typedef struct {
    unsigned char *base;
    size_t         size;
    size_t         off;         /* current bump offset            */
    size_t         high_water;  /* max off ever reached (report)  */
} cs_arena;

/* Wrap a caller-provided block. mem may have any alignment; allocations are
 * aligned internally. Returns CS_OK or CS_ERR_ARG. */
int cs_arena_init(cs_arena *a, void *mem, size_t size);

/* Carve `bytes` aligned to `align` (power of two). NULL on exhaustion. */
void *cs_arena_alloc(cs_arena *a, size_t bytes, size_t align);

/* Convenience: alignof(max_align_t)-aligned carve. */
void *cs_arena_alloc_default(cs_arena *a, size_t bytes);

size_t cs_arena_high_water(const cs_arena *a);

#endif /* CS_ARENA_H */
