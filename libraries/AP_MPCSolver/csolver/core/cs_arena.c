#include "cs_arena.h"

#define CS_ARENA_DEFAULT_ALIGN 16u  /* covers max_align_t on M7 and x86-64 */

int cs_arena_init(cs_arena *a, void *mem, size_t size)
{
    if (!a || (!mem && size))
        return CS_ERR_ARG;
    a->base = (unsigned char *)mem;
    a->size = size;
    a->off = 0;
    a->high_water = 0;
    return CS_OK;
}

void *cs_arena_alloc(cs_arena *a, size_t bytes, size_t align)
{
    size_t start, end;
    uintptr_t p;
    if (!a || !a->base || align == 0 || (align & (align - 1)) != 0)
        return (void *)0;
    p = (uintptr_t)(a->base + a->off);
    start = a->off + (size_t)((align - (p & (align - 1))) & (align - 1));
    end = start + bytes;
    if (end < start || end > a->size)   /* overflow or exhausted */
        return (void *)0;
    a->off = end;
    if (end > a->high_water)
        a->high_water = end;
    return a->base + start;
}

void *cs_arena_alloc_default(cs_arena *a, size_t bytes)
{
    return cs_arena_alloc(a, bytes, CS_ARENA_DEFAULT_ALIGN);
}

size_t cs_arena_high_water(const cs_arena *a)
{
    return a ? a->high_water : 0;
}
