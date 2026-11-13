/* bs_arc — see bs_arc.h.  Every formula below names the Python line it
 * mirrors (model/stage2_arc_clock.py = ac, stage2_arc_clock_mission.py = am,
 * stage2_host_anytime.py = ha, sim_corner_online.py = so). */
#include "bs_arc.h"

#include <math.h>
#include <string.h>
#ifdef BS_ARC_HOST_DEBUG
#include <stdio.h>
#endif

#define NX BS_NX
#define NU BS_NU
#define NN BS_N
#define NV BS_NV
#define NR BS_NROW
#define MM BS_M_COMMIT

/* MODEL DIMENSIONS, decided by the generated header at compile time.
 *   v1 (bs_model_runtime_const.h):    u = (j_l, nu, j_c), NU 3, 50 rows
 *                                     (the two |nu| <= 0.01 clock rows last)
 *   v2 (bs_model_runtime_const_v2.h): u = (j_l, j_c),     NU 2, 48 rows
 *                                     (exact-ZOH e_l, no clock input)
 * The row template is one layout: the (speed, accel) pairs, then +e_l,
 * +e_c, +j_l, +j_c, -e_l, -e_c, -j_l, -j_c (indices 40..47), and the
 * nu pair (48, 49) only when the input exists.  j_l is input 0 and j_c
 * input NU-1 in both.  The row stride the core reads is 10 in both
 * (h[0..5], e[6..6+NU-1], margin at [9]); under NU == 2 slot 8 is
 * zero. */
#if NU == 3
#if NR != 50
#error "bs_arc: NU == 3 (v1, clock input nu) needs BS_NROW == 50"
#endif
#elif NU == 2
#if NR != 48
#error "bs_arc: NU == 2 (v2, no clock input) needs BS_NROW == 48"
#endif
#else
#error "bs_arc: BS_NU must be 3 (v1) or 2 (v2)"
#endif
#if NV != NU * NN
#error "bs_arc: BS_NV must equal BS_NU * BS_N"
#endif
#define ARC_JL 0            /* input index of j_l */
#define ARC_JC (NU - 1)     /* input index of j_c */

#if !BS_ARC_HV_AVAILABLE
/* v1 header: the hover branches are unreachable (a hover chart is refused
 * at build, -4); placeholders so that the translation unit compiles */
#define BS_HV_SPEED_ARC 0.0
#define BS_HV_BAND 0.0
#define BS_HV_EC_BAND 0.0
#define BS_HV_W 0.0
#define BS_HV_GATE_MARGIN 0
#define BS_HV_A_SETTLE 0.0
#endif

/* rule G1's two scan literals (arc_v2.V_ONSET_MIN / N_REST_MIN): stages
 * with v_t <= BS_HV_V_ONSET_MIN are skipped by the onset scan; the rest
 * count is n = max(BS_HV_N_REST_MIN, ceil(v_t / (a_req Ts))).  Not
 * emitted by gen_runtime_const.py (model constants of the rule, not of
 * the plant); a header that does define them wins. */
#ifndef BS_HV_V_ONSET_MIN
#define BS_HV_V_ONSET_MIN 0.3
#endif
#ifndef BS_HV_N_REST_MIN
#define BS_HV_N_REST_MIN 2
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* numpy's degrees / deg2rad: x * (180/pi), x * (pi/180) */
static double deg_of(double rad) { return rad * (180.0 / M_PI); }
static double rad_of(double deg) { return deg * (M_PI / 180.0); }

/* the model's two speed constants (so.VLEG, so.V_TRIM of the certificate) */
static double v_leg(void) { return 11.8 * cos(M_PI / 20.0); }
static double v_trim(void) { return v_leg() - 0.2888; }

/* ((deg + 540) % 360) - 180, Python's floored modulo */
static double wrap180(double deg)
{
    double a = fmod(deg + 540.0, 360.0);
    if (a < 0.0) a += 360.0;
    return a - 180.0;
}

/* ------------------------------------------------------------- the chart */

typedef struct {
    bs_arc_cell *cells;      /* NULL: count only */
    int n, cap;
    int overflow;
} cell_sink;

static void add_cell(cell_sink *s, double px, double py, double thx,
                     double thy, double l, int kind, double chi, int leg,
                     double d_in, int sense, double tube_in, double tube_out)
{
    if (s->cells != NULL) {
        if (s->n >= s->cap) { s->overflow = 1; return; }
        bs_arc_cell *c = &s->cells[s->n];
        c->px = px; c->py = py; c->thx = thx; c->thy = thy;
        c->l = l; c->chi = chi; c->d_in = d_in;
        c->tube_in = tube_in; c->tube_out = tube_out;
        c->kind = (signed char)kind; c->sense = (signed char)sense;
        c->leg = (short)leg;
    }
    s->n++;
}

/* ---- v2e: the corridor's outer boundary in closed form
 * (stage2_arc_clock_mission._ray_slab / _ray_disc / _ray_capsule /
 * tube_outer_reach, the same operations in the same order).  An interval
 * of the ray parameter t is (lo, hi); "none" is a return of 0. */

/* _ray_slab: {t : lo <= s0 + t su <= hi}; su == 0 and s0 inside: the
 * whole line (-inf, inf).  fpclassify(su) == FP_ZERO is exactly the
 * Python's su == 0.0 (both signed zeros; false on NaN) without the
 * firmware's -Wfloat-equal. */
static int ray_slab(double s0, double su, double lo, double hi,
                    double *t_lo, double *t_hi)
{
    if (fpclassify(su) == FP_ZERO) {
        if (lo <= s0 && s0 <= hi) { *t_lo = -HUGE_VAL; *t_hi = HUGE_VAL; return 1; }
        return 0;
    }
    const double t1 = (lo - s0) / su;
    const double t2 = (hi - s0) / su;
    if (t1 <= t2) { *t_lo = t1; *t_hi = t2; }
    else          { *t_lo = t2; *t_hi = t1; }
    return 1;
}

/* _ray_disc: {t : |p + t u - c| <= w} for a unit u, the two roots of
 * t^2 + 2 beta t + gamma = 0 */
static int ray_disc(double px, double py, double ux, double uy, double cx,
                    double cy, double w, double *t_lo, double *t_hi)
{
    const double dx = px - cx;
    const double dy = py - cy;
    const double beta = dx * ux + dy * uy;
    const double gamma = (dx * dx + dy * dy) - w * w;
    const double disc = beta * beta - gamma;
    if (disc < 0.0) return 0;
    const double sq = sqrt(disc);
    *t_lo = -beta - sq;
    *t_hi = -beta + sq;
    return 1;
}

/* _ray_capsule: {t : dist(p + t u, [a, b]) <= w}, the hull of the
 * rectangle's slab interval (0 <= s <= L, |c| <= w, n = (-ty, tx)) and
 * the two disc intervals about a and b, in that order */
static int ray_capsule(double px, double py, double ux, double uy,
                       double ax, double ay, double bx, double by,
                       double tx, double ty, double L, double w,
                       double *t_lo, double *t_hi)
{
    const double rx = px - ax;
    const double ry = py - ay;
    const double s0 = rx * tx + ry * ty;
    const double su = ux * tx + uy * ty;
    const double c0 = rx * (-ty) + ry * tx;
    const double cu = ux * (-ty) + uy * tx;
    int have = 0;
    double lo = 0.0, hi = 0.0;
    double s_lo, s_hi, c_lo, c_hi;
    if (ray_slab(s0, su, 0.0, L, &s_lo, &s_hi)
        && ray_slab(c0, cu, -w, w, &c_lo, &c_hi)) {
        const double a = (s_lo >= c_lo) ? s_lo : c_lo;
        const double b = (s_hi <= c_hi) ? s_hi : c_hi;
        if (a <= b) { lo = a; hi = b; have = 1; }
    }
    for (int e = 0; e < 2; ++e) {
        double a, b;
        if (ray_disc(px, py, ux, uy, e ? bx : ax, e ? by : ay, w, &a, &b)) {
            if (!have) { lo = a; hi = b; have = 1; }
            else {
                if (a < lo) lo = a;
                if (b > hi) hi = b;
            }
        }
    }
    if (!have) return 0;
    *t_lo = lo; *t_hi = hi;
    return 1;
}

/* tube_outer_reach: the first t >= 0 at which p + t u leaves the union of
 * the capsules of radius w about the legs [V[k-1], V[k]] and [V[k],
 * V[k+1]].  Returns 0 (and *e) or -5 when p lies in neither capsule. */
static int tube_outer_reach(double px, double py, double ux, double uy,
                            const double *vx, const double *vy,
                            const double *tx, const double *ty,
                            const double *Ls, int k, double w, double *e_out)
{
    double a0 = 0.0, a1 = 0.0, b0 = 0.0, b1 = 0.0;
    const int hA = ray_capsule(px, py, ux, uy, vx[k - 1], vy[k - 1], vx[k], vy[k],
                               tx[k - 1], ty[k - 1], Ls[k - 1], w, &a0, &a1);
    const int hB = ray_capsule(px, py, ux, uy, vx[k], vy[k], vx[k + 1], vy[k + 1],
                               tx[k], ty[k], Ls[k], w, &b0, &b1);
    const int inA = hA && a0 <= 0.0 && 0.0 <= a1;
    const int inB = hB && b0 <= 0.0 && 0.0 <= b1;
    double e;
    if (inA && inB) {
        e = (a1 >= b1) ? a1 : b1;
    } else if (inA) {
        e = a1;
        if (hB && b0 <= e && b1 > e) e = b1;
    } else if (inB) {
        e = b1;
        if (hA && a0 <= e && a1 > e) e = a1;
    } else {
        return -5;
    }
    *e_out = e;
    return 0;
}

/* The cell count of a straight of length Lp at nominal pace `pace`
 * (am.build_cells_fillet :: n_cells).  Cell rows (v2c) round UP,
 * n = ceil(Lp / (pace Ts)), so that every cell's own pace l / Ts is <=
 * the nominal pace its velocity rows are written about (the forward
 * margin VLEG - pace_k stays >= VLEG - v_nom on every cell); family
 * rows keep the archives' int(round(.)), Python's round-half-even,
 * which is nearbyint under the default rounding mode.  At least one
 * cell either way.  The ceil is taken on the same double the Python
 * sees (np.ceil of Lp / (pace * Ts)), no slack: an exact integer ratio
 * (a 16 m corner-pace part at 8 m/s = 8.0) must give that integer in
 * both. */
static int straight_count(double Lp, double pace, int cell_rows)
{
    const double r = Lp / (pace * BS_TS);
    int n = cell_rows ? (int)ceil(r) : (int)nearbyint(r);
    if (n < 1) n = 1;
    return n;
}

/* am.build_cells_fillet :: straight(): corner-pace cells within 16 m
 * after / 20 m before a vertex, trim cells between. */
static void straight(cell_sink *s, double p0x, double p0y, double thx,
                     double thy, double L, int leg, double chi_at, int n_v,
                     double corner_pace, double vtrim, int cell_rows)
{
    const double a = (leg > 0) ? BS_ARC_WIN_OUT : 0.0;
    const double b = (leg < n_v - 2) ? BS_ARC_WIN_IN : 0.0;
    double Lp[3], pace[3];
    int kind[3], np = 0;
    Lp[np] = (a < L) ? a : L; pace[np] = corner_pace; kind[np] = BS_ARC_CA; np++;
    if (a + b < L) {
        Lp[np] = L - a - b; pace[np] = vtrim; kind[np] = BS_ARC_TA; np++;
        Lp[np] = b; pace[np] = corner_pace; kind[np] = BS_ARC_CA; np++;
    } else if (L > a) {
        Lp[np] = L - a; pace[np] = corner_pace; kind[np] = BS_ARC_CA; np++;
    }
    double sacc = 0.0;
    int first = 1;
    for (int p = 0; p < np; ++p) {
        if (Lp[p] <= 1e-9) continue;
        const int n = straight_count(Lp[p], pace[p], cell_rows);
        const double l = Lp[p] / n;
        for (int j = 0; j < n; ++j) {
            add_cell(s, p0x + thx * sacc, p0y + thy * sacc, thx, thy, l,
                     kind[p], (first && j == 0) ? chi_at : 0.0, leg, 0.0, 0,
                     0.0, 0.0);
            sacc += l;
        }
        first = 0;
    }
}

/* per-vertex scratch of the chart pass: static rather than on the stack
 * (7 x 64 doubles = 3.6 kB; the firmware's solver thread has 8 kB) --
 * the builder runs on one thread at a time */
static struct {
    double tx[BS_ARC_MAX_WP], ty[BS_ARC_MAX_WP], Ls[BS_ARC_MAX_WP];
    double hd[BS_ARC_MAX_WP], chi[BS_ARC_MAX_WP];
    double T[BS_ARC_MAX_WP], rf[BS_ARC_MAX_WP];
} cp;

static int chart_pass(const double *vx, const double *vy, int n_v,
                      cell_sink *s, bs_arc_chart *out, int cell_rows,
                      int hover, double tube_delta)
{
    const double corner_pace = BS_ARC_CORNER_PACE, w = BS_ARC_W;
    const double vtrim = v_trim();
    double *const tx = cp.tx, *const ty = cp.ty, *const Ls = cp.Ls;
    double *const hd = cp.hd, *const chi = cp.chi;
    double *const T = cp.T, *const rf = cp.rf;
    for (int i = 0; i + 1 < n_v; ++i) {
        const double dx = vx[i + 1] - vx[i], dy = vy[i + 1] - vy[i];
        const double L = sqrt(dx * dx + dy * dy);
        if (!(L > 1e-9)) return -3;
        tx[i] = dx / L; ty[i] = dy / L; Ls[i] = L;
        hd[i] = atan2(ty[i], tx[i]);
    }
    chi[0] = 0.0;
    for (int k = 1; k + 1 < n_v; ++k)
        chi[k] = wrap180(deg_of(hd[k] - hd[k - 1]));
    chi[n_v - 1] = 0.0;
    for (int k = 0; k < n_v; ++k) { T[k] = 0.0; rf[k] = 0.0; }
    for (int k = 1; k + 1 < n_v; ++k) {
        const double h = 0.5 * rad_of(fabs(chi[k]));
        if (h < 1e-6) continue;
        /* v2e: (w - delta) / (1 - cos h); tube_delta 0: w / (1 - cos h)
         * bit for bit (w - 0.0 == w) */
        rf[k] = (w - tube_delta) / (1.0 - cos(h));
        T[k] = rf[k] * tan(h);
    }
    for (int k = 1; k + 1 < n_v; ++k) {
        const double ra = Ls[k - 1] * ((k - 1 > 0) ? 0.5 : 1.0);
        const double rb = Ls[k] * ((k + 1 < n_v - 1) ? 0.5 : 1.0);
        const double room = ((ra < rb) ? ra : rb) - 0.5;
        if (T[k] > room) {
            T[k] = (room > 0.2) ? room : 0.2;
            const double h = 0.5 * rad_of(fabs(chi[k]));
            rf[k] = T[k] / tan(h);
        }
    }
    if (out != NULL) {
        for (int k = 0; k < n_v; ++k) out->chi_v[k] = chi[k];
    }
    double chi_carry = 0.0;
    for (int i = 0; i + 1 < n_v; ++i) {
        const double t0 = tx[i], t1 = ty[i];
        const double sx = vx[i] + t0 * T[i], sy = vy[i] + t1 * T[i];
        const double ex = vx[i + 1] - t0 * T[i + 1], ey = vy[i + 1] - t1 * T[i + 1];
        const double L = sqrt((ex - sx) * (ex - sx) + (ey - sy) * (ey - sy));
        if (out != NULL) out->leg_first[i] = s->n;
        straight(s, sx, sy, t0, t1, L, i, chi_carry, n_v, corner_pace, vtrim,
                 cell_rows);
        const int k = i + 1;
        if (k < n_v - 1 && T[k] > 0.0) {
            const double sgn = (chi[k] > 0.0) ? 1.0 : ((chi[k] < 0.0) ? -1.0 : 0.0);
            const double h = rad_of(fabs(chi[k]));
            const double nAx = -t1 * sgn, nAy = t0 * sgn;      /* toward the centre */
            const double Cx = ex + nAx * rf[k], Cy = ey + nAy * rf[k];
            const double arc = rf[k] * h;
            int n_ch = 2;
            {
                const int n1 = (int)ceil(deg_of(h) / BS_ARC_MAX_DPSI);
                const int n2 = (int)ceil(arc / (corner_pace * BS_TS));
                if (n1 > n_ch) n_ch = n1;
                if (n2 > n_ch) n_ch = n2;
            }
            const double dpsi = sgn * h / n_ch;
            const double ang0 = atan2(ey - Cy, ex - Cx);
            double prevx = t0, prevy = t1;
            const double tBx = tx[i + 1], tBy = ty[i + 1];
            /* v2e: the inner corner C_in of the corridor (the legs' inner
             * offset lines at w, not w - delta), algebraic form:
             * C_in = V + (nA + nB) * (w / (1 + tA . tB)) */
            double cinx = 0.0, ciny = 0.0;
            if (tube_delta > 0.0) {
                const double nBx = -tBy * sgn, nBy = tBx * sgn;
                const double f_in = w / (1.0 + (t0 * tBx + t1 * tBy));
                cinx = vx[k] + (nAx + nBx) * f_in;
                ciny = vy[k] + (nAy + nBy) * f_in;
            }
            for (int j = 0; j < n_ch; ++j) {
                const double a0 = ang0 + j * dpsi, a1 = ang0 + (j + 1) * dpsi;
                const double p0x = Cx + rf[k] * cos(a0), p0y = Cy + rf[k] * sin(a0);
                const double p1x = Cx + rf[k] * cos(a1), p1y = Cy + rf[k] * sin(a1);
                const double dx = p1x - p0x, dy = p1y - p0y;
                const double l = sqrt(dx * dx + dy * dy);
                const double thx = dx / l, thy = dy / l;
                const double dchi = wrap180(deg_of(atan2(thy, thx) - atan2(prevy, prevx)));
                const double mx = 0.5 * (p0x + p1x), my = 0.5 * (p0y + p1y);
                const double dA = fabs((mx - vx[k]) * (-t1) + (my - vy[k]) * t0);
                const double dB = fabs((mx - vx[k]) * (-tBy) + (my - vy[k]) * tBx);
                double b_k = 0.0, o_k = 0.0;
                if (tube_delta > 0.0 && s->cells != NULL) {
                    /* v2e: the chord's inward margin b_k (to C_in) and
                     * outward margin o_k (to the outer boundary), on the
                     * filling pass only (the counting pass needs none) */
                    const double nix = -thy * sgn, niy = thx * sgn;   /* n_k = sense * nh */
                    b_k = nix * (cinx - p0x) + niy * (ciny - p0y);
                    double o0, o1;
                    if (tube_outer_reach(p0x, p0y, -nix, -niy, vx, vy, tx, ty, Ls, k, w, &o0) != 0
                        || tube_outer_reach(p1x, p1y, -nix, -niy, vx, vy, tx, ty, Ls, k, w, &o1) != 0)
                        return -5;
                    o_k = (o0 <= o1) ? o0 : o1;
                    if (!(b_k > 0.0 && o_k > 0.0)) return -5;
                }
                add_cell(s, p0x, p0y, thx, thy, l, BS_ARC_FILLET, dchi, i,
                         (dA < dB) ? dA : dB, (int)sgn, b_k, o_k);
                prevx = thx; prevy = thy;
            }
            chi_carry = wrap180(deg_of(atan2(tBy, tBx) - atan2(prevy, prevx)));
        } else {
            chi_carry = 0.0;
        }
        if (out != NULL) out->leg_last[i] = s->n - 1;
    }
    const int k_end = s->n;
    {   /* the trailing straight (trim cells): tail / n with n = round(tail
         * / (V_TRIM Ts)) = 21 under family rows, ceil(.) = 22 under cell
         * rows.  am.build_cells_fillet: th = cells[-1]['th'] (the last
         * cell's tangent = the last leg's), p0 = V[-1]. */
        const double thx = (s->cells != NULL && s->n > 0) ? s->cells[s->n - 1].thx : tx[n_v - 2];
        const double thy = (s->cells != NULL && s->n > 0) ? s->cells[s->n - 1].thy : ty[n_v - 2];
        if (hover) {
            /* v2d (tail = 0): no tail -- the chart ends in ONE hover cell
             * at WP_N (kind BS_ARC_HV, l = 0, pace 0, chi 0, leg n_v - 2);
             * k_end is its index (am.build_cells_fillet tail <= 0). */
            add_cell(s, vx[n_v - 1], vy[n_v - 1], thx, thy, 0.0, BS_ARC_HV,
                     0.0, n_v - 2, 0.0, 0, 0.0, 0.0);
        } else {
            const int n = straight_count(BS_ARC_TAIL, vtrim, cell_rows);
            const double l = BS_ARC_TAIL / n;
            for (int j = 0; j < n; ++j)
                add_cell(s, vx[n_v - 1] + thx * (j * l), vy[n_v - 1] + thy * (j * l),
                         thx, thy, l, BS_ARC_TA, 0.0, n_v - 2, 0.0, 0, 0.0, 0.0);
        }
    }
    if (out != NULL) {
        out->k_end = k_end;
        out->n_cells = s->n;
        out->n_v = n_v;
        out->v_trim = vtrim;
        out->v_leg = v_leg();
        out->corner_pace = corner_pace;
        out->w = w;
        out->cell_rows = cell_rows;
        out->hover = hover;
        out->tube = (tube_delta > 0.0) ? 1 : 0;
        out->tube_delta = tube_delta;
    }
    return 0;
}

int bs_arc_chart_cells_bound(int n_wp, double len_m)
{
    /* loose: 0.75 len + 13 per vertex + tail; the builder counts exactly */
    return (int)(0.75 * len_m) + 13 * n_wp + 32;
}

size_t bs_arc_chart_size(int n_wp, double len_m)
{
    return (size_t)bs_arc_chart_cells_bound(n_wp, len_m) * sizeof(bs_arc_cell);
}

int bs_arc_chart_build_tube(const double *vx, const double *vy, int n_wp,
                            void *mem, size_t mem_len, bs_arc_chart *out,
                            int cell_rows, int hover, double tube_delta)
{
    if (n_wp < 2 || n_wp > BS_ARC_MAX_WP) return -1;
    memset(out, 0, sizeof(*out));
    for (int i = 0; i < n_wp; ++i) { out->vx[i] = vx[i]; out->vy[i] = vy[i]; }
    cell_rows = cell_rows ? 1 : 0;
    hover = hover ? 1 : 0;
    /* arc_v2.build_chart: the hover cell needs cell rows; and the v2d
     * constants / emitted pair, which only the v2 header carries */
    if (hover && (!cell_rows || !BS_ARC_HV_AVAILABLE)) return -4;
    /* arc_v2.configure: 0 <= tube_delta < EC_BAND (= w); NaN refused */
    if (!(tube_delta >= 0.0 && tube_delta < BS_ARC_W)) return -5;
    cell_sink cnt = { NULL, 0, 0, 0 };
    int rc = chart_pass(vx, vy, n_wp, &cnt, NULL, cell_rows, hover, tube_delta);
    if (rc != 0) return rc;
    if ((size_t)cnt.n * sizeof(bs_arc_cell) > mem_len) return -2;
    cell_sink fill = { (bs_arc_cell *)mem, 0, cnt.n, 0 };
    rc = chart_pass(vx, vy, n_wp, &fill, out, cell_rows, hover, tube_delta);
    if (rc != 0) return rc;
    if (fill.overflow || fill.n != cnt.n) return -2;
    out->cells = (bs_arc_cell *)mem;
    return 0;
}

int bs_arc_chart_build_rows(const double *vx, const double *vy, int n_wp,
                            void *mem, size_t mem_len, bs_arc_chart *out,
                            int cell_rows, int hover)
{
    return bs_arc_chart_build_tube(vx, vy, n_wp, mem, mem_len, out,
                                   cell_rows, hover, 0.0);
}

/* the flown configuration: v2e under the v2 header (cell rows, the hover
 * cell, tube fillets at BS_TUBE_DELTA); BS_ARC_TUBE = 0 builds v2d */
#if BS_ARC_TUBE_AVAILABLE && BS_ARC_TUBE
#define BS_ARC_TUBE_DEFAULT BS_TUBE_DELTA
#else
#define BS_ARC_TUBE_DEFAULT 0.0
#endif
int bs_arc_chart_build(const double *vx, const double *vy, int n_wp,
                       void *mem, size_t mem_len, bs_arc_chart *out)
{
    return bs_arc_chart_build_tube(vx, vy, n_wp, mem, mem_len, out, 1,
                                   BS_ARC_HV_AVAILABLE ? BS_ARC_HOVER : 0,
                                   BS_ARC_TUBE_DEFAULT);
}

void bs_arc_world(const bs_arc_chart *ch, int k, const double *xi,
                  double *pn, double *pe, double *vn, double *ve,
                  double *an, double *ae)
{
    const bs_arc_cell *c = &ch->cells[k];
    const double nhx = -c->thy, nhy = c->thx;
    const double pace = c->l / BS_TS;
    if (pn) *pn = c->px + xi[2] * c->thx + xi[3] * nhx;
    if (pe) *pe = c->py + xi[2] * c->thy + xi[3] * nhy;
    if (vn) *vn = (pace + xi[0]) * c->thx + xi[4] * nhx;
    if (ve) *ve = (pace + xi[0]) * c->thy + xi[4] * nhy;
    if (an) *an = xi[1] * c->thx + xi[5] * nhx;
    if (ae) *ae = xi[1] * c->thy + xi[5] * nhy;
}

/* ------------------------------------------------------------ families */

int bs_arc_fam_key(const bs_arc_chart *ch, int k)
{
    const bs_arc_cell *c = &ch->cells[k];
    if (ch->cell_rows || c->kind == BS_ARC_FILLET) return 2 + k;
    return (c->kind == BS_ARC_TA) ? 0 : 1;
}

/* the pace the velocity rows are written about */
static double row_pace(const bs_arc_chart *ch, const bs_arc_cell *c)
{
    /* cell rows (v2c): the cell's own chart pace on every cell.  Family
     * rows: V_TRIM for ta, the corner pace for ca (ac.register_families),
     * the cell's own pace only on a fillet chord (am.fillet_rows). */
    if (ch->cell_rows || c->kind == BS_ARC_FILLET) return c->l / BS_TS;
    return (c->kind == BS_ARC_TA) ? ch->v_trim : ch->corner_pace;
}

/* ha.rows_corner_m layout: (speed_a, accel_a) pairs for a = 0..19, then
 * +e_l, +e_c, +j_l, +j_c, -e_l, -e_c, -j_l, -j_c, and (v1 only) +nu, -nu.
 * Every family of the fillet chart has all NR rows: rows_corner_m drops
 * a velocity row with margin < 1e-9 (pace cos(theta) >= VLEG), which no
 * cell reaches -- the family paces are below VLEG, and under cell rows
 * the rounded-up cell counts keep every cell pace at or below its
 * nominal pace (the gate asserts every margin positive). */
static void template_rows(double *rows)
{
    memset(rows, 0, sizeof(double) * NR * 10);
    for (int a = 0; a < 20; ++a) {
        const double th = rad_of(18.0 * a);
        const double c = cos(th), s = sin(th);
        double *r = &rows[(size_t)(2 * a) * 10];
        r[0] = c; r[4] = s; r[9] = 1.0;
        r = &rows[(size_t)(2 * a + 1) * 10];
        r[1] = c; r[5] = s; r[9] = 1.0;
    }
    const int sgn[2] = { 1, -1 };
    for (int q = 0; q < 2; ++q) {
        double *r = &rows[(size_t)(40 + 4 * q) * 10];
        r[2] = sgn[q]; r[9] = 1.0;                 /* +-e_l */
        r = &rows[(size_t)(41 + 4 * q) * 10];
        r[3] = sgn[q]; r[9] = 1.0;                 /* +-e_c */
        r = &rows[(size_t)(42 + 4 * q) * 10];
        r[6 + ARC_JL] = sgn[q]; r[9] = 1.0;        /* +-j_l */
        r = &rows[(size_t)(43 + 4 * q) * 10];
        r[6 + ARC_JC] = sgn[q]; r[9] = 1.0;        /* +-j_c */
    }
#if NU == 3
    rows[(size_t)48 * 10 + 7] = 1.0;  rows[(size_t)48 * 10 + 9] = 1.0;
    rows[(size_t)49 * 10 + 7] = -1.0; rows[(size_t)49 * 10 + 9] = 1.0;
#endif
}

void bs_arc_margins(const bs_arc_chart *ch, int k, double *m)
{
    const bs_arc_cell *c = &ch->cells[k];
    const double pace = row_pace(ch, c);
    /* arc_v2.rows_corner_m(r_disc): VLEG on every moving cell, the hover
     * disc BS_HV_SPEED_ARC on the hover cell (pace 0: all 20 margins) */
    const double VLEG = (c->kind == BS_ARC_HV) ? BS_HV_SPEED_ARC : ch->v_leg;
    for (int a = 0; a < 20; ++a) {
        const double th = rad_of(18.0 * a);
        m[2 * a] = VLEG - pace * cos(th);
        /* the acceleration faces stay by cell KIND in both row modes
         * (the hover cell: the isotropic tier, D4) */
        m[2 * a + 1] = (c->kind == BS_ARC_TA) ? bs_rc_aniso[a] : BS_ARC_ACCC;
    }
    if (c->kind == BS_ARC_HV) {
        /* arc_v2.hover_rows: +-e_l at BAND_H, +-e_c at EC_BAND_H (D3) */
        m[40] = BS_HV_BAND; m[44] = BS_HV_BAND;
        m[41] = BS_HV_EC_BAND; m[45] = BS_HV_EC_BAND;
    } else {
    m[40] = BS_ARC_BAND; m[44] = BS_ARC_BAND;
    if (c->kind == BS_ARC_FILLET) {
        /* am.fillet_rows: s e_c <= max(w - d_in, 0.05), -s e_c <= w + d_in;
         * a tube chart (v2e, arc_v2.fillet_rows): s e_c <= b_k, -s e_c <=
         * o_k, built once with the chart */
        double in, outm;
        if (ch->tube) {
            in = c->tube_in;
            outm = c->tube_out;
        } else {
            in = ch->w - c->d_in;
            if (in < 0.05) in = 0.05;
            outm = ch->w + c->d_in;
        }
        if (c->sense > 0) { m[41] = in; m[45] = outm; }
        else              { m[41] = outm; m[45] = in; }
    } else {
        m[41] = ch->w; m[45] = ch->w;
    }
    }
    m[42] = BS_JBOX; m[43] = BS_JBOX; m[46] = BS_JBOX; m[47] = BS_JBOX;
#if NU == 3
    m[48] = BS_ARC_M_NU; m[49] = BS_ARC_M_NU;
#endif
}

/* so.dare_pair on the template rows with the given margins, QD = 0 (the
 * Python zeroes NS["QD"]; RL stays).  The pair splits into the
 * longitudinal DARE on (delta, a_l, e_l) with the ML = NU-1 longitudinal
 * inputs -- (j_l, nu) with RL_ = diag(RJ, RNU) + Mu[0:2,0:2] under v1,
 * j_l alone with RL_ = RJ + Mu[0,0] under v2 -- and the lateral DARE on
 * (e_c, v_c, a_c) with j_c, RT_ = RJ + Mu[NU-1,NU-1].  K is NU x 6:
 * rows 0..NU-2 carry KL on columns 0..2, row NU-1 carries KT on 3..5. */
#define ARC_ML (NU - 1)
static int arc_dare(const double *rows, const double *m, double *P6,
                    double *K36)
{
    double Mx[36], Mu[NU * NU];
    memset(Mx, 0, sizeof(Mx));
    memset(Mu, 0, sizeof(Mu));
    for (int i = 0; i < NR; ++i) {
        const double *r = &rows[(size_t)i * 10];
        const double wgt = 2.0 * BS_EPS / (m[i] * m[i]);
        for (int a = 0; a < 6; ++a)
            for (int b = 0; b < 6; ++b)
                Mx[a * 6 + b] += wgt * r[a] * r[b];
        for (int a = 0; a < NU; ++a)
            for (int b = 0; b < NU; ++b)
                Mu[a * NU + b] += wgt * r[6 + a] * r[6 + b];
    }
    double QL[9] = { 0.0, 0, 0, 0, BS_W_QA, 0, 0, 0, BS_W_RL };
    double QT[9] = { BS_W_RC, 0, 0, 0, BS_W_REG, 0, 0, 0, BS_W_REG };
#if NU == 3
    double RLm[4] = { BS_W_RJ, 0, 0, BS_W_RNU };
#else
    double RLm[1] = { BS_W_RJ };
#endif
    double RTm[1] = { BS_W_RJ };
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            QL[a * 3 + b] += Mx[a * 6 + b];
            QT[a * 3 + b] += Mx[(a + 3) * 6 + (b + 3)];
        }
    }
    for (int a = 0; a < ARC_ML; ++a)
        for (int b = 0; b < ARC_ML; ++b) RLm[a * ARC_ML + b] += Mu[a * NU + b];
    RTm[0] += Mu[(NU - 1) * NU + (NU - 1)];
    double PL[9], KL[3 * ARC_ML], PT[9], KT[3];
    if (bs_dare_solve(bs_rc_AL, bs_rc_BL, QL, RLm, 3, ARC_ML, PL, KL, 20000, 1e-14) != 0)
        return -1;
    if (bs_dare_solve(bs_rc_AT, bs_rc_BT, QT, RTm, 3, 1, PT, KT, 20000, 1e-14) != 0)
        return -1;
    memset(P6, 0, sizeof(double) * 36);
    memset(K36, 0, sizeof(double) * NU * NX);
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            P6[a * 6 + b] = PL[a * 3 + b];
            P6[(a + 3) * 6 + (b + 3)] = PT[a * 3 + b];
        }
    for (int a = 0; a < ARC_ML; ++a)
        for (int b = 0; b < 3; ++b) K36[a * 6 + b] = KL[a * 3 + b];
    for (int b = 0; b < 3; ++b) K36[(NU - 1) * 6 + 3 + b] = KT[b];
    return 0;
}

static int terminal_pair(const bs_arc_chart *ch, bs_arc_sched *S, int k,
                         double *P6, double *K36)
{
    if (ch->cells[k].kind == BS_ARC_HV) {
        /* D7: the hover family's pair is mission-independent and emitted
         * by gen_runtime_const.py --v2 (arc_v2.dare_pair(hover_rows()));
         * never solved on target, never a cache miss */
#if BS_ARC_HV_AVAILABLE
        memcpy(P6, bs_rc_arc_Ph, sizeof(double) * 36);
        memcpy(K36, bs_rc_arc_Kh, sizeof(double) * NU * NX);
        S->dare_hover++;
        return 0;
#else
        return -1;                     /* the v1 header carries no v2d pair */
#endif
    }
    const int key = bs_arc_fam_key(ch, k);
    for (int i = 0; i < S->dare_n; ++i) {
        if (S->dare[i].key == key) {
            memcpy(P6, S->dare[i].P, sizeof(double) * 36);
            memcpy(K36, S->dare[i].K, sizeof(double) * NU * NX);
            return 0;
        }
    }
    double m[NR];
    bs_arc_margins(ch, k, m);
    S->dare_solves++;
    {   /* D8: the on-target solve, timed by the driver's clock hook */
        const unsigned long t0 = S->micros ? S->micros() : 0UL;
        const int rc = arc_dare(S->rows, m, P6, K36);
        if (S->micros) S->us_dare += S->micros() - t0;
        if (rc != 0) { S->dare_fail++; return -1; }
    }
    const int slot = (S->dare_n < BS_ARC_DARE_CACHE) ? S->dare_n++ : S->dare_next;
    S->dare_next = (slot + 1) % BS_ARC_DARE_CACHE;
    /* family rows: the two straight families (the first two pairs solved)
     * stay resident, the fillet pairs cycle through slots 2..5.  Cell
     * rows: every cell is a family and the tip walks the chart, so the
     * cache is a plain ring (the tip cell of the shift, the K rounds and
     * usually the next tick's rounds is one of the last few solved). */
    if (!ch->cell_rows && S->dare_next < 2) S->dare_next = 2;
    S->dare[slot].key = key;
    memcpy(S->dare[slot].P, P6, sizeof(double) * 36);
    memcpy(S->dare[slot].K, K36, sizeof(double) * NU * NX);
    return 0;
}

/* ---------------------------------------------------------- transitions */

void bs_arc_plant(const double *xi, const double *u, double *y)
{
    for (int i = 0; i < NX; ++i) {
        double acc = 0.0;
        for (int j = 0; j < NX; ++j) acc += bs_A[(size_t)i * NX + j] * xi[j];
        for (int j = 0; j < NU; ++j) acc += bs_B[(size_t)i * NU + j] * u[j];
        y[i] = acc;
    }
}

/* so.T_switch: pairs (0,4), (1,5), (2,3) */
static void t_switch(double chi_deg, double *T)
{
    const double c = cos(rad_of(chi_deg)), s = sin(rad_of(chi_deg));
    static const int pa[3] = { 0, 1, 2 }, pb[3] = { 4, 5, 3 };
    for (int i = 0; i < NX * NX; ++i) T[i] = (i % (NX + 1) == 0) ? 1.0 : 0.0;
    for (int p = 0; p < 3; ++p) {
        const int ip = pa[p], iq = pb[p];
        T[ip * NX + ip] = c; T[ip * NX + iq] = s;
        T[iq * NX + ip] = -s; T[iq * NX + iq] = c;
    }
}

/* ac.crossed */
static int crossed(const bs_arc_chart *ch, int kn, const double *st,
                   double pace_prev)
{
    const bs_arc_cell *c = &ch->cells[kn];
    const double el = st[2], ec = st[3];
    if (!(fabs(c->chi) > 0.0)) return el >= 0.0;
    const double cc = cos(rad_of(c->chi)), sn = sin(rad_of(c->chi));
    const double h = 0.5 * rad_of(c->chi);
    const double va = pace_prev + st[0], vc = st[4];
    const double zone = BS_ARC_VERTEX_BOX;      /* fillet cells carry no zone */
    const int by_pos = (el >= -BS_ARC_VERTEX_BOX) && (el * cos(h) + ec * sin(h) >= 0.0);
    const int by_hdg = (fabs(el) <= zone + 0.5) && (va * va + vc * vc >= 0.3 * 0.3)
                       && ((va * cc + vc * sn) >= va);
    return by_pos || by_hdg;
}

static void mat6_mul(const double *a, const double *b, double *c)
{
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j) {
            double acc = 0.0;
            for (int p = 0; p < NX; ++p) acc += a[i * NX + p] * b[p * NX + j];
            c[i * NX + j] = acc;
        }
}

static void mat6_vec(const double *a, const double *x, double *y)
{
    for (int i = 0; i < NX; ++i) {
        double acc = 0.0;
        for (int p = 0; p < NX; ++p) acc += a[i * NX + p] * x[p];
        y[i] = acc;
    }
}

/* ac.transition: given the raw stepped state y (e_l relative to cell k+1's
 * start), hold or advance through every crossed boundary.  `n_cells` is
 * the chart length the rule sees: ch->n_cells (the vehicle, the natural
 * rollout), or k_end on a hover chart with the gate closed (arc_v2.
 * HoverStop.assign's `trunc = cells[:k_end]`: the hover cell is not
 * assignable, the last straight cell holds at the chart's end). */
static int transition_n(const bs_arc_chart *ch, int n_cells, int k,
                        const double *y, double *T, double *c)
{
    double Tl[NX * NX], cl[NX];
    for (int i = 0; i < NX * NX; ++i) Tl[i] = (i % (NX + 1) == 0) ? 1.0 : 0.0;
    for (int i = 0; i < NX; ++i) cl[i] = 0.0;
    const double pace_k = ch->cells[k].l / BS_TS;
    if (k + 1 >= n_cells || !crossed(ch, k + 1, y, pace_k)) {
        cl[2] = ch->cells[k].l;                         /* HOLD */
        if (T) memcpy(T, Tl, sizeof(Tl));
        if (c) memcpy(c, cl, sizeof(cl));
        return k;
    }
    int kn = k + 1;
    double state[NX];
    memcpy(state, y, sizeof(state));
    for (;;) {
        const bs_arc_cell *cn = &ch->cells[kn];
        const double chi = cn->chi;
        const double pace_old = ch->cells[kn - 1].l / BS_TS;
        const double pace_new = cn->l / BS_TS;
        double ck[NX] = { 0, 0, 0, 0, 0, 0 };
        {   /* so.off_vec */
            const double ca = cos(rad_of(chi)), sa = sin(rad_of(chi));
            ck[0] = pace_old * ca - pace_new;
            ck[4] = -pace_old * sa;
        }
        if (fabs(chi) > 0.0) {
            double Tk[NX * NX], tmp[NX * NX], tv[NX];
            t_switch(chi, Tk);
            mat6_mul(Tk, Tl, tmp); memcpy(Tl, tmp, sizeof(Tl));
            mat6_vec(Tk, cl, tv);
            for (int i = 0; i < NX; ++i) cl[i] = tv[i] + ck[i];
            mat6_vec(Tk, state, tv);
            for (int i = 0; i < NX; ++i) state[i] = tv[i] + ck[i];
        } else {
            for (int i = 0; i < NX; ++i) { cl[i] += ck[i]; state[i] += ck[i]; }
        }
        double st2[NX];
        memcpy(st2, state, sizeof(st2));
        st2[2] -= cn->l;
        if (kn + 1 < n_cells && crossed(ch, kn + 1, st2, pace_new)) {
            cl[2] -= cn->l;                    /* re-base onto cell kn+1 */
            state[2] -= cn->l;
            kn += 1;
            continue;
        }
        break;
    }
    if (T) memcpy(T, Tl, sizeof(Tl));
    if (c) memcpy(c, cl, sizeof(cl));
    return kn;
}

int bs_arc_transition(const bs_arc_chart *ch, int k, const double *y,
                      double *T, double *c)
{
    return transition_n(ch, ch->n_cells, k, y, T, c);
}

/* ------------------------------------------ v2d: the stop at WP_N */

/* arc_v2.force_to_end: the forced transition of a latched stage -- from
 * cell k through every remaining cell into the hover cell (rotations
 * and offsets applied, each intermediate cell's length subtracted), so
 * the returned e_l is the signed shortfall to WP_N in the hover frame. */
static int force_to_end(const bs_arc_chart *ch, int k, const double *y,
                        double *T, double *c)
{
    double Tl[NX * NX], cl[NX], state[NX];
    for (int i = 0; i < NX * NX; ++i) Tl[i] = (i % (NX + 1) == 0) ? 1.0 : 0.0;
    for (int i = 0; i < NX; ++i) cl[i] = 0.0;
    memcpy(state, y, sizeof(state));
    const int k_end = ch->n_cells - 1;
    int kn = k;
    while (kn < k_end) {
        kn += 1;
        const bs_arc_cell *cn = &ch->cells[kn];
        const double chi = cn->chi;
        const double pace_old = ch->cells[kn - 1].l / BS_TS;
        const double pace_new = cn->l / BS_TS;
        double ck[NX] = { 0, 0, 0, 0, 0, 0 };
        {   /* so.off_vec */
            const double ca = cos(rad_of(chi)), sa = sin(rad_of(chi));
            ck[0] = pace_old * ca - pace_new;
            ck[4] = -pace_old * sa;
        }
        if (fabs(chi) > 0.0) {
            double Tk[NX * NX], tmp[NX * NX], tv[NX];
            t_switch(chi, Tk);
            mat6_mul(Tk, Tl, tmp); memcpy(Tl, tmp, sizeof(Tl));
            mat6_vec(Tk, cl, tv);
            for (int i = 0; i < NX; ++i) cl[i] = tv[i] + ck[i];
            mat6_vec(Tk, state, tv);
            for (int i = 0; i < NX; ++i) state[i] = tv[i] + ck[i];
        } else {
            for (int i = 0; i < NX; ++i) { cl[i] += ck[i]; state[i] += ck[i]; }
        }
        if (kn < k_end) {
            cl[2] -= cn->l;
            state[2] -= cn->l;
        }
    }
    if (T) memcpy(T, Tl, sizeof(Tl));
    if (c) memcpy(c, cl, sizeof(cl));
    return kn;
}

/* arc_v2.hover_frame: (e_l, e_c) of a state in the hover cell's frame
 * (frame point WP_N, pace 0, the last leg's tangent): the state itself
 * in the hover cell, shifted by the last straight cell's length one cell
 * before it, none further back.  Returns 0 when the frame is undefined. */
static int hover_frame(const bs_arc_chart *ch, int k, const double *xi,
                       double *el, double *ec)
{
    const int k_end = ch->n_cells - 1;
    if (k >= k_end) { *el = xi[2]; *ec = xi[3]; return 1; }
    if (k == k_end - 1) { *el = xi[2] - ch->cells[k].l; *ec = xi[3]; return 1; }
    return 0;
}

/* arc_v2.in_zone: arrival -- e_l >= -BAND_H in the hover frame */
static int in_zone(const bs_arc_chart *ch, int k, const double *xi)
{
    double el, ec;
    return hover_frame(ch, k, xi, &el, &ec) && el >= -BS_HV_BAND;
}

/* arc_v2.in_box: the settle predicate -- speed <= R_H, |e_l| <= BAND_H
 * and |e_c| <= BAND_H in the hover frame, |a| <= A_SETTLE */
static int in_box(const bs_arc_chart *ch, int k, const double *xi)
{
    double el, ec;
    if (!hover_frame(ch, k, xi, &el, &ec)) return 0;
    const double speed = hypot(ch->cells[k].l / BS_TS + xi[0], xi[4]);
    const double acc = hypot(xi[1], xi[5]);
    return speed <= BS_HV_SPEED_ARC && fabs(el) <= BS_HV_BAND
           && fabs(ec) <= BS_HV_BAND && acc <= BS_HV_A_SETTLE;
}

/* arc_v2.dist_to_go (the gate distance of record, open_dist; the G1
 * onset scan uses the backward suffix sums instead): the distance along the
 * chart from a state in cell k0 to WP_N (negative past it); the cell
 * lengths summed in chart order from k0 as the Python's sum() does. */
static double dist_to_go(const bs_arc_chart *ch, int k0, const double *xi0)
{
    const int k_end = ch->n_cells - 1;
    if (k0 >= k_end) return -xi0[2];
    double sum = 0.0;
    for (int c = k0; c < k_end; ++c) sum += ch->cells[c].l;
    return sum - xi0[2];
}

/* arc_v2.stop_profile::integrate, the ZOH integration in the Python's
 * operation order; returns v_final, a_final, the distance and the
 * minima, and writes js[0..n1+n2+n3-1].  Jp from a_final = 0. */
static void stop_integrate(double v0, double a0, int n1, int n2, int n3,
                           double J, double *js, double *x_out,
                           double *v_out, double *a_out, double *Jp_out,
                           double *amin_out, double *vmin_out)
{
    const double Ts = BS_TS;
    const double Ts3 = Ts * Ts * Ts;          /* Ts ** 3, exact (a power of two) */
    const double Jp = (J * n1 * Ts - a0) / (n3 * Ts);
    const int n = n1 + n2 + n3;
    for (int i = 0; i < n1; ++i) js[i] = -J;
    for (int i = n1; i < n1 + n2; ++i) js[i] = 0.0;
    for (int i = n1 + n2; i < n; ++i) js[i] = Jp;
    double v = v0, a = a0, x = 0.0, amin = a, vmin = v;
    for (int i = 0; i < n; ++i) {
        const double j = js[i];
        x += v * Ts + 0.5 * a * Ts * Ts + j * Ts3 / 6;
        v += a * Ts + 0.5 * j * Ts * Ts; a += j * Ts;
        if (a < amin) amin = a;
        if (v < vmin) vmin = v;
    }
    *x_out = x; *v_out = v; *a_out = a; *Jp_out = Jp;
    *amin_out = amin; *vmin_out = vmin;
}

/* arc_v2.stop_profile (G1: used only by brake_profile_exact, once per
 * chart, for the stop of record): js = [-J] n1 + [0] n2 + [+Jp] n3 from (v0, a0) to
 * EXACT rest on the Ts grid -- Jp fixed by a_final = 0, J by v_final = 0
 * (v_final affine in J: J = -vA / (vB - vA) from the integrations at
 * J = 0 and J = 1).  Returns 0 when infeasible (J <= 0, J > JB, |Jp| >
 * JB, a below -AB, v reversing, rest inexact), else 1 with js, d_stop
 * (x), n_stop, J, Jp, a_min. */
static int stop_profile(double v0, double a0, double AB, double JB,
                        int n1, int n2, int n3, double *js, double *d,
                        int *n_stop, double *J_out, double *Jp_out,
                        double *amin_out)
{
    double x, vA, vB, a, Jp, amin, vmin;
    stop_integrate(v0, a0, n1, n2, n3, 0.0, js, &x, &vA, &a, &Jp, &amin, &vmin);
    stop_integrate(v0, a0, n1, n2, n3, 1.0, js, &x, &vB, &a, &Jp, &amin, &vmin);
    if (fabs(vB - vA) < 1e-12) return 0;
    const double J = -vA / (vB - vA);
    double vf, af;
    stop_integrate(v0, a0, n1, n2, n3, J, js, &x, &vf, &af, &Jp, &amin, &vmin);
    if (J <= 0.0 || J > JB + 1e-9 || fabs(Jp) > JB + 1e-9 || amin < -AB - 1e-9
        || vmin < -1e-6 || fabs(vf) > 1e-9 || fabs(af) > 1e-9)
        return 0;
    *d = x; *n_stop = n1 + n2 + n3; *J_out = J; *Jp_out = Jp; *amin_out = amin;
    return 1;
}

/* arc_v2.plan_states: the natural-clock rollout of U from (k, xi) --
 * transition() on the full chart, the hover cell assignable -- cells
 * ks[0..N] and states xs[0..N]. */
static void plan_states(const bs_arc_chart *ch, int k, const double *xi,
                        const double *U, int *ks, double *xs)
{
    double x[NX], y[NX], T[NX * NX], c[NX], xn[NX];
    memcpy(x, xi, sizeof(x));
    ks[0] = k;
    memcpy(&xs[0], xi, sizeof(x));
    int kk = k;
    for (int t = 0; t < NN; ++t) {
        bs_arc_plant(x, &U[NU * t], y);
        kk = bs_arc_transition(ch, kk, y, T, c);
        mat6_vec(T, y, xn);
        for (int i = 0; i < NX; ++i) x[i] = xn[i] + c[i];
        ks[t + 1] = kk;
        memcpy(&xs[(size_t)(t + 1) * NX], x, sizeof(x));
    }
}

/* arc_v2.seed_integrate over the grid trapezoid of arc_v2.trapezoid_grid
 * (js_of(s) = [j1] * n_r + [0.0] * n_h + [j3] * n_r2): the exact ZOH
 * integration from (v0, a0) in the Python's operation order, the jerk of
 * stage i generated in place (no array: the profile may run past N).
 * Returns (x, v) at its end. */
static void seed_integrate(double v0, double a0, int n_r, int n_h, int n_r2,
                           double j1, double j3, double *x_out, double *v_out)
{
    const double Ts = BS_TS;
    const double Ts3 = Ts * Ts * Ts;          /* Ts ** 3, exact (a power of two) */
    const int n = n_r + n_h + n_r2;
    double v = v0, a = a0, x = 0.0;
    for (int i = 0; i < n; ++i) {
        const double j = (i < n_r) ? j1 : ((i < n_r + n_h) ? 0.0 : j3);
        x += v * Ts + 0.5 * a * Ts * Ts + j * Ts3 / 6;
        v += a * Ts + 0.5 * j * Ts * Ts; a += j * Ts;
    }
    *x_out = x; *v_out = v;
}

/* arc_v2.trapezoid_closed_form: the plateau deceleration a_s of the
 * jerk-limited trapezoid in acceleration from (v, a) resting at distance
 * d (ramp a -> -a_s at -JB, hold, ramp -a_s -> 0 at +JB), the smaller
 * root of (v' / (2 JB)) a_s^2 - d' a_s + v'^2 / 2 = 0 after the pre-phase
 * to a = 0 on the same jerk line.  Returns BS_ARC_SEED_OK with a_s and
 * the continuous duration, else the no-root reason (d' <= 0, disc < 0,
 * a_s > AB, a_s < -a), in the Python's operation order.  The powers are
 * explicit products in the Python's order (tau0 * tau0, tau0 * tau0 *
 * tau0, vp * vp * vp): no libm pow, so the result does not depend on the
 * target's libm (glibc on the host, newlib on the Cube). */
static int trapezoid_closed_form(double v, double a, double d, double JB, double AB,
                                 double *s_out, double *dur_out)
{
    const double tau0 = fabs(a) / JB;
    const double vp = v + a * a / (2.0 * JB);
    const double x0 = (a >= 0)
        ? (v * tau0 + a * (tau0 * tau0) / 2 - JB * (tau0 * tau0 * tau0) / 6)
        : -(vp * tau0 - JB * (tau0 * tau0 * tau0) / 6);
    const double dp = d - x0;
    if (dp <= 0.0) return BS_ARC_SEED_DPRIME;
    const double disc = dp * dp - (vp * vp * vp) / JB;
    if (disc < 0.0) return BS_ARC_SEED_DISC;
    const double s = JB * (dp - sqrt(disc)) / vp;
    if (s > AB) return BS_ARC_SEED_AB;
    if (s < -a) return BS_ARC_SEED_RAMPUP;
    const double T_h = vp / s - s / JB;
    const double tau1 = (a + s) / JB, tau2 = s / JB;
    *s_out = s;
    *dur_out = tau1 + T_h + tau2;
    return BS_ARC_SEED_OK;
}

/* the onset scan's distance table (arc_v2.dist_suffix): S[kk] = sum_{c =
 * kk}^{k_end - 1} l_c accumulated BACKWARD from WP_N (S[k_end] = 0, S[c]
 * = S[c + 1] + l_c), from the engage checkpoints (bs_arc_stop_model.
 * dist_cp, every dist_stride-th S): a lookup starts at the checkpoint at
 * or above kk (k_end with S = 0 past the last one) and walks down, the
 * Python's additions in the Python's order, so the value is bit-identical
 * to its table.  The scan visits the stages t = N-1 down to M, whose
 * cells never increase (transition_n holds k or advances past it, never
 * lowers it), so one walk serves the whole scan: at most dist_stride - 1
 * + (ks[N-1] - ks[M]) additions per tick.  The re-start branch (kk above
 * the walk's position) is defensive only; it cannot be reached from the
 * scan. */
typedef struct { int c; double s; } dist_walk;

static void dist_checkpoints(const bs_arc_chart *ch, bs_arc_stop_model *sm)
{
    const int k_end = ch->n_cells - 1;
    int stride = (k_end + (BS_HV_DIST_CP - 1) - 1) / (BS_HV_DIST_CP - 1);
    if (stride < 1) stride = 1;
    sm->dist_stride = stride;
    double s = 0.0;
    for (int c = k_end - 1; c >= 0; --c) {
        s = s + ch->cells[c].l;
        if (c % stride == 0) sm->dist_cp[c / stride] = s;
    }
}

static void dist_walk_start(const bs_arc_chart *ch, const bs_arc_stop_model *sm,
                            dist_walk *w, int kk)
{
    const int k_end = ch->n_cells - 1;
    const int stride = sm->dist_stride;
    const int c_up = ((kk + stride - 1) / stride) * stride;   /* checkpoint at or above kk */
    if (c_up >= k_end) { w->c = k_end; w->s = 0.0; }
    else { w->c = c_up; w->s = sm->dist_cp[c_up / stride]; }
}

static double suffix_at(const bs_arc_chart *ch, const bs_arc_stop_model *sm,
                        dist_walk *w, int kk)
{
    if (w->c < 0 || kk > w->c) dist_walk_start(ch, sm, w, kk);
    while (w->c > kk) { w->c--; w->s = w->s + ch->cells[w->c].l; }
    return w->s;
}

/* arc_v2.brake_seed (rule G1, doc/V2D-G1-DESIGN-2026-09-25.md, variant
 * cdB), once per tick while the gate is closed: one natural rollout of
 * the warm start U from (k, xi), one onset scan over t = N-1 down to M
 * (v_t = pace + delta, d_t = S[kk] - e_l, or -e_l in the hover cell;
 * skip v_t <= BS_HV_V_ONSET_MIN or d_t <= 0; a_req = v_t^2 / (2 d_t),
 * n = max(BS_HV_N_REST_MIN, ceil(v_t / (a_req Ts))), rest = t + n; the
 * first stage of the downward scan with a_req <= a_nom and rest <= N -
 * margin is t_b); at t_b the jerk-box trapezoid at the constant
 * deceleration a_s (trapezoid_closed_form), on the Ts grid with a_s
 * re-solved so the seed ends at rest (v_final affine in a_s: two ZOH
 * integrations), written into j_l from t_b (stages at or beyond N not
 * written).  Returns BS_ARC_SEED_OK (U seeded, the record in *sd),
 * BS_ARC_SEED_NONE when no stage qualifies (U untouched), or a no-root
 * reason BS_ARC_SEED_DPRIME .. BS_ARC_SEED_LONG (U untouched; *sd carries
 * t_b, v_b, a_b, dist_wp, a_req, a_nom): the gate stays closed this tick
 * (owner policy 2026-09-25, no fallback seed).
 * The rollout scratch is static (31 x 6 doubles + 31 ints), off the
 * 8 kB solver-thread stack; the seed itself needs no array. */
/* The margin guard, re-derived for G1.  Memory is safe for ANY margin:
 * the seed write is clamped at stage N (t_b + i < N, as the Python's
 * `if t_b + i < N`) and the latch is min(N, rest stage).  What the range
 * protects is the rule: a qualifying stage needs t >= M and t + n <= N -
 * margin with n >= BS_HV_N_REST_MIN, so margin > N - M - BS_HV_N_REST_MIN
 * leaves the scan nothing to find (the gate would open only on entry into
 * the hover cell), and a negative margin admits stops that rest beyond
 * the horizon (the latch then sits at N while the plan is still moving). */
#if BS_HV_N_REST_MIN < 1
#error "BS_HV_N_REST_MIN < 1: the rest stage must lie after the onset stage"
#endif
#if BS_HV_GATE_MARGIN < 0 || BS_HV_GATE_MARGIN > BS_N - BS_M_COMMIT - BS_HV_N_REST_MIN
#error "BS_HV_GATE_MARGIN out of range: need 0 <= margin <= N - M - BS_HV_N_REST_MIN (a qualifying stage t >= M with t + n <= N - margin must exist; a negative margin admits stops resting past the horizon)"
#endif
/* the grid trapezoid's stage count is bounded by the scan (a qualifying
 * stage has 2 d_t / (v_t Ts) <= N, and the trapezoid's duration is that
 * of a stop over the same distance, about 2 d' / v'); this cap only
 * refuses a non-finite or absurd duration before it is cast to int -- a
 * C-side guard with no Python counterpart, unreachable on the scan's
 * stages (reason BS_ARC_SEED_LONG, the gate stays closed) */
#define BS_HV_SEED_STAGES_CAP (4 * NN)
static struct {
    int ks[NN + 1];
    double xs[(NN + 1) * NX];
} bsd;

static int brake_seed(const bs_arc_chart *ch, const bs_arc_stop_model *sm,
                      int k, const double *xi, double *U, double AB,
                      double a_nom, int margin, bs_arc_seed *sd)
{
    const int M = MM;
    const int k_end = ch->n_cells - 1;
    const double Ts = BS_TS;
    const double JB = BS_JBOX;
    int *const ks = bsd.ks;
    double *const xs = bsd.xs;
    plan_states(ch, k, xi, U, ks, xs);
    /* arc_v2.onset_scan */
    dist_walk w = { -1, 0.0 };               /* started at the first lookup */
    int t_b = -1;
    double a_req = 0.0, v_b = 0.0, a_b = 0.0, d_b = 0.0;
    for (int t = NN - 1; t >= M; --t) {
        const int kk = ks[t];
        const double *xt = &xs[(size_t)t * NX];
        const double v_t = ch->cells[kk].l / BS_TS + xt[0];
        const double d_t = (kk >= k_end) ? (-xt[2]) : (suffix_at(ch, sm, &w, kk) - xt[2]);
        if (v_t <= BS_HV_V_ONSET_MIN || d_t <= 0.0) continue;
        const double ar = v_t * v_t / (2.0 * d_t);
        /* n and rest as doubles: exact integers, no int overflow on a far
         * stage (v_t / (a_req Ts) = 2 d_t / (v_t Ts)) */
        double n = ceil(v_t / (ar * Ts) - 1e-12);
        if (n < BS_HV_N_REST_MIN) n = BS_HV_N_REST_MIN;
        const double rest = t + n;
        if (ar <= a_nom && rest <= NN - margin) {
            t_b = t; a_req = ar; v_b = v_t; a_b = xt[1]; d_b = d_t;
            break;
        }
    }
    if (t_b < 0) return BS_ARC_SEED_NONE;
    sd->t_b = t_b; sd->v_b = v_b; sd->a_b = a_b; sd->dist_wp = d_b;
    sd->a_req = a_req; sd->a_nom = a_nom;
    double s_c = 0.0, dur = 0.0;
    const int why = trapezoid_closed_form(v_b, a_b, d_b, JB, AB, &s_c, &dur);
    if (why != BS_ARC_SEED_OK) return why;
    if (!(dur / Ts <= BS_HV_SEED_STAGES_CAP)) return BS_ARC_SEED_LONG;
    /* arc_v2.trapezoid_grid */
    int n_r = (int)ceil((a_b + s_c) / (JB * Ts) - 1e-9);
    if (n_r < 1) n_r = 1;
    int n_r2 = (int)ceil(s_c / (JB * Ts) - 1e-9);
    if (n_r2 < 1) n_r2 = 1;
    int n_tot = (int)rint(dur / Ts);          /* Python's round(): half to even */
    if (n_tot < n_r + n_r2) n_tot = n_r + n_r2;
    const int n_h = n_tot - n_r - n_r2;
    double x, vA, vB;
    seed_integrate(v_b, a_b, n_r, n_h, n_r2, -(a_b + 0.0) / (n_r * Ts), 0.0 / (n_r2 * Ts), &x, &vA);
    seed_integrate(v_b, a_b, n_r, n_h, n_r2, -(a_b + 1.0) / (n_r * Ts), 1.0 / (n_r2 * Ts), &x, &vB);
    const double s = -vA / (vB - vA);
    const double j1 = -(a_b + s) / (n_r * Ts), j3 = s / (n_r2 * Ts);
    /* the third integration serves the record only (d_stop, rest_short:
     * the Python's x_r, compared by the parity), once per mission on the
     * gate tick; the seed write below does not use it */
    double x_r, v_r;
    seed_integrate(v_b, a_b, n_r, n_h, n_r2, j1, j3, &x_r, &v_r);
    for (int i = 0; i < n_tot && t_b + i < NN; ++i)
        U[NU * (t_b + i) + ARC_JL] = (i < n_r) ? j1 : ((i < n_r + n_h) ? 0.0 : j3);
    sd->n_stop = n_tot; sd->rest_stage = t_b + n_tot;
    sd->d_stop = x_r; sd->rest_short = d_b - x_r;
    sd->a_s_cont = s_c; sd->a_s = s; sd->a_min = -s;
    sd->n_r = n_r; sd->n_h = n_h; sd->n_r2 = n_r2;
    sd->j_max = (fabs(j1) >= fabs(j3)) ? fabs(j1) : fabs(j3);
    return BS_ARC_SEED_OK;
}

/* arc_v2.brake_profile_exact: the symmetric stop of record from (v0, a0 =
 * 0): n1 = round(AB / JB / Ts) stages at -J (Python's round = half to
 * even = rint under the default rounding mode), the smallest n2 with J <=
 * JB at 0, n1 at +J, through stop_profile.  Returns 1 with d_stop, n_stop,
 * J, a_max (= -a_min) and v_residual (0 by stop_profile's exact rest), 0
 * when no such stop exists or its profile exceeds the scratch (2 n1 + n2
 * > BS_HV_STOP_MAX stages, a pace beyond 70 m/s at the record's faces;
 * the Python has no such cap -- both refuse the chart). */
#define BS_HV_STOP_MAX 64
static int brake_profile_exact(double v0, double a0, double AB, double JB,
                               double *d_stop, int *n_stop, double *J,
                               double *a_max, double *v_res)
{
    const double Ts = BS_TS;
    int n1 = (int)rint(AB / JB / Ts);
    if (n1 < 1) n1 = 1;
    int n2 = (int)ceil(v0 / (JB * n1 * Ts * Ts) - n1 - 1e-9);
    if (n2 < 0) n2 = 0;
    if (2 * n1 + n2 > BS_HV_STOP_MAX) return 0;
    double js[BS_HV_STOP_MAX], Jp, amin;
    if (!stop_profile(v0, a0, AB, JB, n1, n2, n1, js, d_stop, n_stop, J, &Jp, &amin))
        return 0;
    *a_max = -amin; *v_res = 0.0;
    return 1;
}

/* arc_v2.stop_params (HoverStop.__init__): the stop model of a hover
 * chart -- the last straight cell's pace v_c, its family's rear face AB
 * (arc_v2.rear_face: the acceleration row at raster index 10, template
 * row 21; 5.609 on a trim cell, 5.66 on a corner-pace one) and forward
 * face AF (forward_face: raster index 0, template row 1), and the
 * symmetric stop of record from v_c, and rule G1's nominal deceleration
 * a_nom = v_c^2 / (2 d_stop).  Returns 0 when that stop does not exist
 * (the Python's AssertionError). */
static int stop_params(const bs_arc_chart *ch, bs_arc_stop_model *sm)
{
    const int k_end = ch->n_cells - 1;
    double m[NR];
    memset(sm, 0, sizeof(*sm));
    bs_arc_margins(ch, k_end - 1, m);
    sm->AB = m[2 * 10 + 1];
    sm->AF = m[1];
    sm->v_c = ch->cells[k_end - 1].l / BS_TS;
    if (!brake_profile_exact(sm->v_c, 0.0, sm->AB, BS_JBOX, &sm->d_stop, &sm->n_stop,
                             &sm->J, &sm->a_max, &sm->v_residual))
        return 0;
    /* rule G1's nominal deceleration (HoverStop.__init__: a_nom =
     * v_c * v_c / (2.0 * d_stop), the same product, no libm pow) */
    sm->a_nom = sm->v_c * sm->v_c / (2.0 * sm->d_stop);
    sm->valid = 1;
    return 1;
}

/* ------------------------------------------------------- the schedule */

static int rot_slot(bs_arc_sched *S, const double *T, int is_identity)
{
    if (is_identity) return 0;
    const int slot = S->n_rot++;
    memcpy(&S->rot[(size_t)slot * NX * NX], T, sizeof(double) * NX * NX);
    return slot;
}

static int off_slot(bs_arc_sched *S, const double *c)
{
    int nz = 0;
    for (int i = 0; i < NX; ++i) if (fabs(c[i]) > 0.0) nz = 1;
    if (!nz) return 0;
    const int slot = S->n_off++;
    memcpy(&S->off[(size_t)slot * NX], c, sizeof(double) * NX);
    return slot;
}

void bs_arc_assign(const bs_arc_chart *ch, bs_arc_sched *S, int k0,
                   const double *xi0, const double *U)
{
    double xi[NX], y[NX], T[NX * NX], c[NX];
    memcpy(xi, xi0, sizeof(xi));
    int k = k0;
    S->ks[0] = k0;
    S->els[0] = xi0[2];
    /* slot 0 = identity / zero, the rest allocated per stage */
    for (int i = 0; i < NX * NX; ++i) S->rot[i] = (i % (NX + 1) == 0) ? 1.0 : 0.0;
    for (int i = 0; i < NX; ++i) S->off[i] = 0.0;
    S->n_rot = 1; S->n_off = 1;
    S->family[0] = 0; S->rotation[0] = 0; S->offset[0] = 0;
    /* v2d (arc_v2.HoverStop.assign): gate closed -> the virtual tail on
     * the chart without the hover cell; gate open -> the natural
     * transition with every stage >= latch short of the hover cell
     * forced into it.  A tail chart takes neither branch. */
    const int k_end = ch->n_cells - 1;
    const int vtail = ch->hover && !S->hv.gate_open && k0 < k_end;
    const int gated = ch->hover && !vtail;
    const int t_from = S->hv.latch;            /* -1: none */
    int forced = 0, n_vt = 0;
    for (int t = 0; t < NN; ++t) {
        bs_arc_plant(xi, &U[NU * t], y);
        int kn;
        if (vtail) {
            if (k == k_end - 1 && y[2] >= 0.0) {
                /* the last straight cell repeats beyond WP_N: T = I, c = 0 */
                for (int i = 0; i < NX * NX; ++i) T[i] = (i % (NX + 1) == 0) ? 1.0 : 0.0;
                for (int i = 0; i < NX; ++i) c[i] = 0.0;
                kn = k_end - 1; n_vt++;
            } else {
                kn = transition_n(ch, k_end, k, y, T, c);
            }
        } else {
            kn = bs_arc_transition(ch, k, y, T, c);
            if (gated && kn < k_end && t_from >= 0 && (t + 1) >= t_from) {
                kn = force_to_end(ch, k, y, T, c);
                forced++;
            }
        }
        double xn[NX];
        mat6_vec(T, y, xn);
        for (int i = 0; i < NX; ++i) xi[i] = xn[i] + c[i];
        /* identity exactly when no seam was crossed (the composition
         * starts from I and is only touched at a rotating boundary) */
        int ident = 1;
        for (int i = 0; i < NX * NX && ident; ++i) {
            const double e = (i % (NX + 1) == 0) ? 1.0 : 0.0;
            if (fabs(T[i] - e) > 0.0) ident = 0;
        }
        S->rotation[t + 1] = rot_slot(S, T, ident);
        S->offset[t + 1] = off_slot(S, c);
        S->family[t + 1] = 0;
        k = kn;
        S->ks[t + 1] = k;
        S->els[t + 1] = xi[2];
    }
    memcpy(S->xiN, xi, sizeof(xi));
    if (ch->hover) {
        /* the latch moves earlier on natural entry (the first stage > 0
         * whose cell is >= k_end), never later */
        int t_first = -1;
        if (!vtail) {
            for (int t = 0; t <= NN; ++t)
                if (S->ks[t] >= k_end) { t_first = t; break; }
            if (t_first > 0 && (S->hv.latch < 0 || t_first < S->hv.latch))
                S->hv.latch = t_first;
        }
        S->hv.forced = forced; S->hv.t_nat = t_first; S->hv.n_vtail = n_vt;
    }
}

int bs_arc_render(const bs_arc_chart *ch, bs_arc_sched *S, double trust)
{
    for (int t = 0; t < NN; ++t) {
        double *m = &S->mstk[(size_t)t * NR];
        bs_arc_margins(ch, S->ks[t], m);
        /* ac.ArcHorizon: the schedule trust region on the e_l rows */
        const double up = S->els[t] + trust, dn = -S->els[t] + trust;
        if (up < m[40]) m[40] = up;
        if (dn < m[44]) m[44] = dn;
    }
    if (terminal_pair(ch, S, S->ks[NN], S->P, S->K) != 0) return -1;
    /* v2d (arc_v2.TiltedHorizon): w_t by the family of X_{t+1} -- BS_HV_W
     * on the hover cell, BS_SPEED_TILT_W elsewhere.  On a tail chart every
     * entry is BS_SPEED_TILT_W and the core's product is the same. */
    for (int t = 0; t < NN; ++t)
        S->tilt_w[t] = (ch->hover && S->ks[t + 1] >= ch->n_cells - 1)
                       ? (double)BS_HV_W : (double)BS_SPEED_TILT_W;
    bs_schedule *sc = &S->sched;
    memset(sc, 0, sizeof(*sc));
    sc->family = S->family;
    sc->rotation = S->rotation;
    sc->offset = S->offset;
    sc->length = BS_ARC_NSLOT;
    sc->periodic = 0;
    sc->rows_tab = S->rows;
    sc->P_tab = S->P;
    sc->K_tab = S->K;
    sc->rot_tab = S->rot;
    sc->off_tab = S->off;
    sc->ring_mask = 0;
#if BS_ARC
    sc->mstk = S->mstk;
    sc->tilt_w = S->tilt_w;
    sc->w_dj = S->w_dj;
#endif
    return 0;
}

/* ------------------------------------------------------------ the driver */

void bs_arc_driver_init(bs_arc_driver *d, const bs_arc_chart *ch,
                        bs_arc_sched *S, double *arena, size_t arena_len,
                        int q, int rounds, double trust)
{
    memset(d, 0, sizeof(*d));
    d->chart = ch;
    d->S = S;
    d->arena = arena;
    d->arena_len = arena_len;
    d->q = q;
    d->rounds = rounds;
    d->trust = trust;
    d->seed_ramp = 1;                         /* v2b: the ramp seed */
    memset(S, 0, sizeof(*S));
    template_rows(S->rows);
    S->hv.latch = -1; S->hv.t_nat = -1;        /* no latch before engage */
    d->latch = -1; d->open_tick = -1; d->arrival_tick = -1; d->settle_tick = -1;
    d->seed.tick = d->seed.t_b = d->seed.n_stop = d->seed.rest_stage = -1;
}

/* doc section 8, the ramp verbatim: JB = BS_JBOX, AF = the forward face
 * of cell 0's family (template row 1 = the acceleration row at raster
 * index 0), v_target = cells[0].l / BS_TS; only j_l (input 0) is set. */
void bs_arc_ramp_seed(const bs_arc_chart *ch, double *U)
{
    const double JB = BS_JBOX, Ts = BS_TS;
    const double v_target = ch->cells[0].l / Ts;
    double m[NR];
    bs_arc_margins(ch, 0, m);
    const double AF = m[1];
    memset(U, 0, sizeof(double) * NV);
    double v = 0.0, a = 0.0;
    for (int t = 0; t < NN; ++t) {
        double j;
        if (v + a * a / (2.0 * JB) < v_target - 1e-9) {
            j = (a < AF - 1e-9) ? JB : 0.0;
            if (a + j * Ts > AF) j = (AF - a) / Ts;
        } else {
            j = (a > 1e-9) ? -JB : 0.0;
            if (a + j * Ts < 0.0) j = -a / Ts;
        }
        U[NU * t + ARC_JL] = j;
        v += a * Ts + 0.5 * j * Ts * Ts;
        a += j * Ts;
    }
}

static unsigned long clock_us(const bs_arc_driver *d)
{
    return d->micros ? d->micros() : 0UL;
}

static int solve_pinned(bs_arc_driver *d, int npin, int q,
                        bs_newton_stats *st, double *lam)
{
    const unsigned long t0 = clock_us(d);
    if (bs_problem_init_pinned(&d->prob, &d->S->sched, 0, npin, d->arena,
                               d->arena_len) != BS_OK) return -1;
    const unsigned long t1 = clock_us(d);
    if (bs_newton_pinned(&d->prob, d->U, d->xi, q, npin, lam, st) != BS_OK)
        return -2;
    const unsigned long t2 = clock_us(d);
    d->us_init += t1 - t0;
    d->us_newton += t2 - t1;
    return 0;
}

/* reset the stop state (driver and schedule) -- at engage */
static void stop_reset(bs_arc_driver *d)
{
    bs_arc_sched *S = d->S;
    S->hv.gate_open = 0; S->hv.latch = -1;
    S->hv.forced = 0; S->hv.t_nat = -1; S->hv.n_vtail = 0;
    S->us_dare = 0; S->micros = d->micros;
    d->gate_open = 0; d->latch = -1; d->open_tick = -1;
    d->open_speed = 0.0; d->open_dist = 0.0;
    d->arrival_tick = -1; d->settle_tick = -1; d->settled = 0;
    d->forced = 0; d->n_hv = 0; d->us_dare = 0;
    memset(&d->seed, 0, sizeof(d->seed));
    d->seed.tick = -1; d->seed.t_b = -1; d->seed.n_stop = -1; d->seed.rest_stage = -1;
    d->no_root = 0; d->no_root_tick = -1; d->no_root_tb = -1; d->no_root_reason = 0;
    memset(&d->stop, 0, sizeof(d->stop));
}

int bs_arc_engage(bs_arc_driver *d)
{
    const bs_arc_chart *ch = d->chart;
    d->k = 0;
    d->kw = 0;
    stop_reset(d);
    /* v2d: the stop model of record, once (HoverStop.__init__ ->
     * stop_params); a hover chart without the symmetric stop from the
     * last straight cell's pace is refused as the Python refuses it */
    if (ch->hover && !stop_params(ch, &d->stop)) return -4;
    if (!(d->S->w_dj >= 0.0)) return -5;       /* negative / NaN weight */
#if BS_NU != 2
    /* the jerk-change penalty is defined for the two-input model only
     * (arc_v2.jerk_change_gram raises at NU != 2; bs_solver ignores it) */
    if (d->S->w_dj > 0.0) return -5;
#endif
    /* rule G1's distance checkpoints (HoverStop.__init__ -> dist_suffix) */
    if (ch->hover) dist_checkpoints(ch, &d->stop);
    memset(d->xi, 0, sizeof(d->xi));
    d->xi[0] = -(ch->cells[0].l / BS_TS);          /* at rest on vertex 0 */
    /* v2d: arc_v2._run tests in_zone on EVERY state sample, the engage
     * state (sample 0) included -- a chart whose first cell is the last
     * straight cell (two waypoints closer than BS_HV_BAND) arrives at
     * tick 0; the settle cap then counts from 0 */
    if (ch->hover && in_zone(ch, d->k, d->xi)) d->arrival_tick = 0;
    if (d->seed_ramp) bs_arc_ramp_seed(ch, d->U);  /* v2b */
    else memset(d->U, 0, sizeof(d->U));            /* the archives' U = 0 */
    d->seed_iters = 0;
    for (int r = 0; r < BS_ARC_SEED_ROUNDS; ++r) {
        bs_arc_assign(ch, d->S, d->k, d->xi, d->U);
        if (bs_arc_render(ch, d->S, d->trust) != 0) return -1;
        if (bs_problem_init_pinned(&d->prob, &d->S->sched, 0, 0, d->arena,
                                   d->arena_len) != BS_OK) return -2;
        /* so.MixedHorizonAff.newton(U, xi, tol=1e-7, maxit=300) */
        double *grad = d->prob.nt_grad;
#ifdef BS_ARC_HOST_DEBUG
        {
            double v0 = 0.0;
            bs_eval(&d->prob, d->U, d->xi, &v0, NULL, NULL);
            printf("seed round %d: J0 %.10g ks[N] %d n_rot %d n_off %d\n", r, v0, d->S->ks[NN], d->S->n_rot, d->S->n_off);
            printf("  P diag %.6g %.6g %.6g %.6g %.6g %.6g  P01 %.6g\n", d->S->P[0], d->S->P[7], d->S->P[14], d->S->P[21], d->S->P[28], d->S->P[35], d->S->P[1]);
            printf("  mstk[0][0..4] %.6g %.6g %.6g %.6g %.6g  [40..47] %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g",
                   d->S->mstk[0], d->S->mstk[1], d->S->mstk[2], d->S->mstk[3], d->S->mstk[4],
                   d->S->mstk[40], d->S->mstk[41], d->S->mstk[42], d->S->mstk[43], d->S->mstk[44], d->S->mstk[45], d->S->mstk[46], d->S->mstk[47]);
            for (int i = 48; i < NR; ++i) printf(" %.6g", d->S->mstk[i]);
            printf("\n");
            printf("  off[1] %.6g %.6g %.6g %.6g %.6g %.6g  xi0 %.6g\n", d->S->off[6], d->S->off[7], d->S->off[8], d->S->off[9], d->S->off[10], d->S->off[11], d->xi[0]);
        }
        int iters_r = 0;
#endif
        for (int it = 0; it < BS_ARC_SEED_MAXIT; ++it) {
            double value = 0.0;
            if (bs_eval(&d->prob, d->U, d->xi, &value, grad, NULL) != BS_OK)
                return -3;
            double gmax = 0.0;
            for (int i = 0; i < NV; ++i) {
                const double a = fabs(grad[i]);
                if (a > gmax) gmax = a;
            }
            if (gmax <= BS_ARC_SEED_TOL) break;
            bs_newton_stats st = {0, 0, 0, 0};
            double lam = 0.0;
            if (bs_newton(&d->prob, d->U, d->xi, 1, &lam, &st) != BS_OK)
                return -3;
            d->seed_iters++;
#ifdef BS_ARC_HOST_DEBUG
            iters_r++;
#endif
            if (st.armijo_fail || st.factor_fail || st.no_step) break;
            if (d->yield) d->yield();
        }
#ifdef BS_ARC_HOST_DEBUG
        {
            double v1 = 0.0;
            bs_eval(&d->prob, d->U, d->xi, &v1, NULL, NULL);
            printf("  -> %d iters, J1 %.10g, U[0] %.8g U[NU-1] %.8g, U[NV-NU] %.8g\n", iters_r, v1, d->U[0], d->U[NU - 1], d->U[NV - NU]);
        }
#endif
    }
    memcpy(d->committed, d->U, sizeof(double) * NU * MM);
    d->engaged = 1;
    d->finished = 0;
    return 0;
}

/* arc_v2.HoverStop.tick_gate: once per tick after the shift-append and
 * before the rounds.  The gate opens at the first tick rule G1 yields a
 * brake seed (or the vehicle is in the hover cell); the seed is written
 * into U and the latch set to min(N, its rest stage).  A tick whose scan
 * qualifies a stage but whose quadratic has no admissible root leaves
 * the gate closed (HoverStop.no_root: counted, the first one recorded).
 * `tick` = the Python's stop.tick (kw + 1). */
static void tick_gate(bs_arc_driver *d, int tick)
{
    const bs_arc_chart *ch = d->chart;
    bs_arc_sched *S = d->S;
    const int k_end = ch->n_cells - 1;
    if (S->hv.gate_open) return;
    /* HoverStop.open_dist / open_speed: the vehicle's distance to WP_N
     * along the chart and its speed at the tick the gate opens (the seed
     * record's dist_vehicle).  Computed on the opening tick only, as the
     * Python does (dist_to_go is a forward sum over k..k_end-1, O(n_cells);
     * brake_seed changes neither d->k nor d->xi, so the value is the same
     * before or after it). */
    if (d->k >= k_end) {
        S->hv.gate_open = 1; d->open_tick = tick;
        d->open_dist = dist_to_go(ch, d->k, d->xi);
        d->open_speed = ch->cells[d->k].l / BS_TS + d->xi[0];
        return;
    }
    /* AB = the last straight family's rear acceleration face and a_nom =
     * v_c^2 / (2 d_stop), both from the stop model computed at engage
     * (HoverStop.tick_gate passes self.params['AB'] and self.a_nom) */
    bs_arc_seed sd;
    memset(&sd, 0, sizeof(sd));
    const int r = brake_seed(ch, &d->stop, d->k, d->xi, d->U, d->stop.AB,
                             d->stop.a_nom, BS_HV_GATE_MARGIN, &sd);
    if (r == BS_ARC_SEED_NONE) return;
    if (r != BS_ARC_SEED_OK) {                 /* no admissible root: stays closed */
        if (d->no_root == 0) {
            d->no_root_tick = tick; d->no_root_tb = sd.t_b; d->no_root_reason = r;
        }
        d->no_root++;
        return;
    }
    S->hv.gate_open = 1; d->open_tick = tick;
    d->open_dist = dist_to_go(ch, d->k, d->xi);
    d->open_speed = ch->cells[d->k].l / BS_TS + d->xi[0];
    sd.tick = tick;
    d->seed = sd;
    S->hv.latch = (sd.rest_stage < NN) ? sd.rest_stage : NN;
}

/* the settle predicate (arc_v2.HoverStop.prefix_settled) on the plan
 * just solved: X0..XM over the schedule's (ks, T, c) all in the box */
static int prefix_settled(const bs_arc_chart *ch, const bs_arc_sched *S,
                          const double *xi0, const double *U)
{
    double x[NX], y[NX], xn[NX];
    memcpy(x, xi0, sizeof(x));
    if (!in_box(ch, S->ks[0], x)) return 0;
    for (int t = 0; t < MM; ++t) {
        bs_arc_plant(x, &U[NU * t], y);
        const double *T = &S->rot[(size_t)S->rotation[t + 1] * NX * NX];
        const double *c = &S->off[(size_t)S->offset[t + 1] * NX];
        mat6_vec(T, y, xn);
        for (int i = 0; i < NX; ++i) x[i] = xn[i] + c[i];
        if (!in_box(ch, S->ks[t + 1], x)) return 0;
    }
    return 1;
}

int bs_arc_tick(bs_arc_driver *d, double *u0)
{
    const bs_arc_chart *ch = d->chart;
    bs_arc_sched *S = d->S;
    if (!d->engaged) return -1;
    if (!ch->hover && d->k >= ch->k_end) { d->finished = 1; return 1; }
    S->us_dare = 0; S->micros = d->micros;
    /* replay the committed prefix (an identity on the pinned path) */
    memcpy(d->U, d->committed, sizeof(double) * NU * MM);
    memcpy(u0, d->U, sizeof(double) * NU);
    /* applied row excess (diagnostic, the family's own margins) */
    {
        double m[NR];
        bs_arc_margins(ch, d->k, m);
        double worst = 0.0;
        for (int i = 0; i < NR; ++i) {
            const double *r = &S->rows[(size_t)i * 10];
            double g = 0.0;
            for (int j = 0; j < NX; ++j) g += r[j] * d->xi[j];
            for (int j = 0; j < NU; ++j) g += r[6 + j] * u0[j];
            if (g - m[i] > worst) worst = g - m[i];
        }
        d->viol = worst;
    }
    const int k_prev = d->k;
    double xi_prev[NX];
    memcpy(xi_prev, d->xi, sizeof(xi_prev));
    /* plant step with the vehicle's own clock rule */
    {
        double y[NX], T[NX * NX], c[NX], xn[NX];
        bs_arc_plant(d->xi, u0, y);
        const int kn = bs_arc_transition(ch, d->k, y, T, c);
        mat6_vec(T, y, xn);
        for (int i = 0; i < NX; ++i) d->xi[i] = xn[i] + c[i];
        if (kn == d->k) d->holds++;
        if (kn > d->k + 1) d->skips++;
        d->k = kn;
        /* v2d: arrival = zone entry, on the post-step state (its index is
         * kw + 1 = the Python's arrival tick).  arc_v2._run tests in_zone
         * at the top of the NEXT loop iteration, so the sample its last
         * iteration produces (index max_ticks) is never tested there,
         * while it is tested here; a run whose first zone entry is that
         * sample fails in both (no settle), differing only in the
         * recorded arrival -- no archive run comes near it. */
        if (ch->hover && d->arrival_tick < 0 && in_zone(ch, d->k, d->xi))
            d->arrival_tick = d->kw + 1;
    }
    /* warm start: shift + terminal law of the previous assignment's tip
     * (ac.run: assign from the PRE-step state along the final plan) */
    const unsigned long t_shift = clock_us(d);
    {
        bs_arc_assign(ch, S, k_prev, xi_prev, d->U);
        double P6[36], K36[NU * NX];
        if (terminal_pair(ch, S, S->ks[NN], P6, K36) != 0) return -2;
        double Us[NV];
        for (int i = 0; i < NV - NU; ++i) Us[i] = d->U[i + NU];
        for (int i = 0; i < NU; ++i) {
            double acc = 0.0;
            for (int j = 0; j < NX; ++j) acc += K36[i * NX + j] * S->xiN[j];
            Us[NV - NU + i] = acc;
        }
        memcpy(d->U, Us, sizeof(Us));
    }
    /* v2d: the gate test, once, after the shift-append and before the rounds */
    if (ch->hover) tick_gate(d, d->kw + 1);
    /* rounds: assignment from the iterate, pinned solve */
    d->reassign = 0;
    d->bt = d->armijo_fail = d->no_step = d->factor_fail = 0;
    d->us_init = d->us_newton = 0;
    d->us_assign = clock_us(d) - t_shift;
    int ks_prev[NN + 1];
    int have_prev = 0;
    for (int r = 0; r < d->rounds; ++r) {
        const unsigned long ta = clock_us(d);
        bs_arc_assign(ch, S, d->k, d->xi, d->U);
        if (have_prev && memcmp(ks_prev, S->ks, sizeof(ks_prev)) != 0) d->reassign++;
        memcpy(ks_prev, S->ks, sizeof(ks_prev));
        have_prev = 1;
        if (bs_arc_render(ch, S, d->trust) != 0) return -2;
        d->us_assign += clock_us(d) - ta;
        if (ch->hover && r == d->rounds - 1) {
            /* the stop state as the last solved plan sees it (the archive's
             * solve_latch_stage / solve_forced_stages of the last round) */
            d->latch = S->hv.latch; d->forced = S->hv.forced;
            int n_hv = 0;
            for (int t = 0; t <= NN; ++t) if (S->ks[t] >= ch->n_cells - 1) n_hv++;
            d->n_hv = n_hv;
        }
        bs_newton_stats st = {0, 0, 0, 0};
        double lam = 0.0;
        if (solve_pinned(d, MM, d->q, &st, &lam) != 0) return -3;
        if (r == 0) d->lambda0 = lam;
        d->bt += st.backtracks;
        d->armijo_fail += st.armijo_fail;
        d->no_step += st.no_step;
        d->factor_fail += st.factor_fail;
    }
    /* v2d: the settle predicate on the last solved plan (before the
     * consistency re-assignment below rewrites the schedule slots) */
    int settled_now = 0;
    if (ch->hover && !d->settled && d->rounds > 0
        && prefix_settled(ch, S, d->xi, d->U)) {
        d->settled = 1; d->settle_tick = d->kw + 1; settled_now = 1;
    }
    /* the fixed point of the assignment (diagnostic) */
    {
        int ks_solved[NN + 1];
        memcpy(ks_solved, S->ks, sizeof(ks_solved));
        /* re-assign over the schedule itself: bs_arc_assign rewrites the
         * slots, which are not read again this tick.  (An earlier draft
         * declared a scratch bs_arc_sched here and never used it -- a
         * 29 kB automatic that only -O elides, on an 8 kB thread stack.) */
        bs_arc_assign(ch, S, d->k, d->xi, d->U);
        d->consistent = (memcmp(ks_solved, S->ks, sizeof(ks_solved)) == 0);
    }
    /* commit stage M-1: the input flown M ticks from now */
    memmove(d->committed, d->committed + NU, sizeof(double) * NU * (MM - 1));
    memcpy(d->committed + NU * (MM - 1), d->U + NU * (MM - 1), sizeof(double) * NU);
    d->kw++;
    d->us_dare = (unsigned)S->us_dare;       /* D8, on every chart */
    if (ch->hover) {
        d->gate_open = S->hv.gate_open;
        return settled_now ? 2 : 0;          /* never finished by cell index */
    }
    if (d->k >= ch->k_end) d->finished = 1;
    return 0;
}

void bs_arc_lift(const bs_arc_driver *d, double *pn, double *pe, int *kn)
{
    const bs_arc_chart *ch = d->chart;
    double xi[NX], y[NX], T[NX * NX], c[NX], xn[NX];
    memcpy(xi, d->xi, sizeof(xi));
    int k = d->k;
    for (int i = 0; i <= MM; ++i) {
        bs_arc_world(ch, k, xi, &pn[i], &pe[i], NULL, NULL, NULL, NULL);
        if (kn) kn[i] = k;
        if (i == MM) break;
        bs_arc_plant(xi, &d->committed[NU * i], y);
        k = bs_arc_transition(ch, k, y, T, c);
        mat6_vec(T, y, xn);
        for (int j = 0; j < NX; ++j) xi[j] = xn[j] + c[j];
    }
}
