#include "cs_model.h"
#include "cs_arena.h"

#include "../model/generated/turn_r5_common.h"
#include "../model/generated/turn_r5_N15.h"
#include "../model/generated/turn_r5_N20.h"
#include "../model/generated/turn_r5_N30.h"
#include "../model/generated/turn_r5_N60.h"

/* the generated model and the core must be built with the same scalar */
typedef char cs_real_matches_casadi_real
    [(sizeof(cs_real) == sizeof(casadi_real)) ? 1 : -1];
typedef char cs_int_matches_casadi_int
    [(sizeof(int) == sizeof(casadi_int)) ? 1 : -1];

typedef int (*cs_fn)(const casadi_real **, casadi_real **, casadi_int *,
                     casadi_real *, int);
typedef int (*cs_work_fn)(casadi_int *, casadi_int *, casadi_int *,
                          casadi_int *);

/* every generated function this library calls — the _work() maxima are taken
 * over THIS table, so adding a function here automatically resizes buffers */
typedef struct {
    cs_fn fn;
    cs_work_fn work;
} cs_entry;

static const cs_entry k_funcs[] = {
    { turn_r5_expl_ode_fun,      turn_r5_expl_ode_fun_work },
    { turn_r5_expl_vde_forw,     turn_r5_expl_vde_forw_work },
    { turn_r5_constr_h,          turn_r5_constr_h_work },
    { turn_r5_constr_h_jac,      turn_r5_constr_h_jac_work },
    { turn_r5_N15_disc_step,     turn_r5_N15_disc_step_work },
    { turn_r5_N15_disc_step_jac, turn_r5_N15_disc_step_jac_work },
    { turn_r5_N20_disc_step,     turn_r5_N20_disc_step_work },
    { turn_r5_N20_disc_step_jac, turn_r5_N20_disc_step_jac_work },
    { turn_r5_N30_disc_step,     turn_r5_N30_disc_step_work },
    { turn_r5_N30_disc_step_jac, turn_r5_N30_disc_step_jac_work },
    { turn_r5_N60_disc_step,     turn_r5_N60_disc_step_work },
    { turn_r5_N60_disc_step_jac, turn_r5_N60_disc_step_jac_work },
};
#define CS_N_FUNCS ((int)(sizeof(k_funcs) / sizeof(k_funcs[0])))

static const int k_variants[] = { 15, 20, 30, 60 };
static const cs_fn k_step[]     = { turn_r5_N15_disc_step,
                                    turn_r5_N20_disc_step,
                                    turn_r5_N30_disc_step,
                                    turn_r5_N60_disc_step };
static const cs_fn k_step_jac[] = { turn_r5_N15_disc_step_jac,
                                    turn_r5_N20_disc_step_jac,
                                    turn_r5_N30_disc_step_jac,
                                    turn_r5_N60_disc_step_jac };
#define CS_N_VARIANTS ((int)(sizeof(k_variants) / sizeof(k_variants[0])))

/* shared call buffers — carved ONCE from the caller's arena, fixed after */
static struct {
    const casadi_real **arg;
    casadi_real       **res;
    casadi_int         *iw;
    casadi_real        *w;
    casadi_int sz_arg, sz_res, sz_iw, sz_w;
    int ready;
} g;

static void work_maxima(casadi_int *a, casadi_int *r, casadi_int *i,
                        casadi_int *w)
{
    casadi_int sa, sr, si, sw;
    int k;
    *a = *r = *i = *w = 0;
    for (k = 0; k < CS_N_FUNCS; ++k) {
        k_funcs[k].work(&sa, &sr, &si, &sw);
        if (sa > *a) *a = sa;
        if (sr > *r) *r = sr;
        if (si > *i) *i = si;
        if (sw > *w) *w = sw;
    }
}

int cs_abi_version(void) { return CS_ABI_VERSION; }
int cs_real_size(void) { return (int)sizeof(cs_real); }

void cs_model_dims(int *nx, int *nu, int *nh)
{
    if (nx) *nx = 6;
    if (nu) *nu = 2;
    if (nh) *nh = 5;
}

int cs_model_num_variants(void) { return CS_N_VARIANTS; }
int cs_model_variant_N(int i)
{
    return (i >= 0 && i < CS_N_VARIANTS) ? k_variants[i] : -1;
}

size_t cs_model_min_arena(void)
{
    casadi_int sa, sr, si, sw;
    work_maxima(&sa, &sr, &si, &sw);
    return (size_t)sa * sizeof(const casadi_real *)
         + (size_t)sr * sizeof(casadi_real *)
         + (size_t)si * sizeof(casadi_int)
         + (size_t)sw * sizeof(casadi_real)
         + 4 * 16;                     /* per-carve alignment slack */
}

int cs_model_init(void *mem, size_t size)
{
    cs_arena arena;
    casadi_int sa, sr, si, sw;
    int rc;
    g.ready = 0;
    rc = cs_arena_init(&arena, mem, size);
    if (rc != CS_OK)
        return rc;
    work_maxima(&sa, &sr, &si, &sw);
    /* THE work-array rule: buffers sized per _work() maxima, not n_in/n_out */
    g.arg = (const casadi_real **)
        cs_arena_alloc_default(&arena, (size_t)sa * sizeof(*g.arg));
    g.res = (casadi_real **)
        cs_arena_alloc_default(&arena, (size_t)sr * sizeof(*g.res));
    g.iw = (casadi_int *)
        cs_arena_alloc_default(&arena, (size_t)si * sizeof(*g.iw));
    g.w = (casadi_real *)
        cs_arena_alloc_default(&arena, (size_t)sw * sizeof(*g.w));
    if ((sa && !g.arg) || (sr && !g.res) || (si && !g.iw) || (sw && !g.w))
        return CS_ERR_ARENA;
    g.sz_arg = sa; g.sz_res = sr; g.sz_iw = si; g.sz_w = sw;
    g.ready = 1;
    return CS_OK;
}

int cs_model_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w)
{
    if (!g.ready)
        return CS_ERR_INIT;
    if (sz_arg) *sz_arg = (int)g.sz_arg;
    if (sz_res) *sz_res = (int)g.sz_res;
    if (sz_iw)  *sz_iw = (int)g.sz_iw;
    if (sz_w)   *sz_w = (int)g.sz_w;
    return CS_OK;
}

/* call helper: load declared ins/outs, leave scratch slots to the wrapper */
static int cs_call(cs_fn fn, const cs_real *const *in, int n_in,
                   cs_real *const *out, int n_out)
{
    int k;
    if (!g.ready)
        return CS_ERR_INIT;
    for (k = 0; k < n_in; ++k)
        g.arg[k] = in[k];
    for (k = 0; k < n_out; ++k)
        g.res[k] = out[k];
    return fn(g.arg, g.res, g.iw, g.w, 0) == 0 ? CS_OK : CS_ERR_MODEL;
}

int cs_ode(const cs_real *x, const cs_real *u, const cs_real *T,
           cs_real *xdot)
{
    const cs_real *in[3];
    cs_real *out[1];
    if (!x || !u || !T || !xdot)
        return CS_ERR_ARG;
    in[0] = x; in[1] = u; in[2] = T;
    out[0] = xdot;
    return cs_call(turn_r5_expl_ode_fun, in, 3, out, 1);
}

int cs_vde_forw(const cs_real *x, const cs_real *Sx, const cs_real *Sw,
                const cs_real *u, const cs_real *T,
                cs_real *xdot, cs_real *Sxdot, cs_real *Swdot)
{
    const cs_real *in[5];
    cs_real *out[3];
    if (!x || !Sx || !Sw || !u || !T || !xdot || !Sxdot || !Swdot)
        return CS_ERR_ARG;
    in[0] = x; in[1] = Sx; in[2] = Sw; in[3] = u; in[4] = T;
    out[0] = xdot; out[1] = Sxdot; out[2] = Swdot;
    return cs_call(turn_r5_expl_vde_forw, in, 5, out, 3);
}

int cs_h(const cs_real *x, const cs_real *u, cs_real *h)
{
    const cs_real *in[2];
    cs_real *out[1];
    if (!x || !u || !h)
        return CS_ERR_ARG;
    in[0] = x; in[1] = u;
    out[0] = h;
    return cs_call(turn_r5_constr_h, in, 2, out, 1);
}

int cs_h_jac(const cs_real *x, const cs_real *u,
             cs_real *h, cs_real *Jx, cs_real *Ju)
{
    const cs_real *in[2];
    cs_real *out[3];
    if (!x || !u || !h || !Jx || !Ju)
        return CS_ERR_ARG;
    in[0] = x; in[1] = u;
    out[0] = h; out[1] = Jx; out[2] = Ju;
    return cs_call(turn_r5_constr_h_jac, in, 2, out, 3);
}

static int variant_index(int N)
{
    int k;
    for (k = 0; k < CS_N_VARIANTS; ++k)
        if (k_variants[k] == N)
            return k;
    return -1;
}

int cs_step(int N, const cs_real *x, const cs_real *u, const cs_real *T,
            cs_real *x_next)
{
    const cs_real *in[3];
    cs_real *out[1];
    int v = variant_index(N);
    if (v < 0)
        return CS_ERR_BADN;
    if (!x || !u || !T || !x_next)
        return CS_ERR_ARG;
    in[0] = x; in[1] = u; in[2] = T;
    out[0] = x_next;
    return cs_call(k_step[v], in, 3, out, 1);
}

int cs_step_jac(int N, const cs_real *x, const cs_real *u, const cs_real *T,
                cs_real *x_next, cs_real *A, cs_real *B, cs_real *bT)
{
    const cs_real *in[3];
    cs_real *out[4];
    int v = variant_index(N);
    if (v < 0)
        return CS_ERR_BADN;
    if (!x || !u || !T || !x_next || !A || !B || !bT)
        return CS_ERR_ARG;
    in[0] = x; in[1] = u; in[2] = T;
    out[0] = x_next; out[1] = A; out[2] = B; out[3] = bT;
    return cs_call(k_step_jac[v], in, 3, out, 4);
}
