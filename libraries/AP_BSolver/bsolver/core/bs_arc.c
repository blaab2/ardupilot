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
                     double d_in, int sense)
{
    if (s->cells != NULL) {
        if (s->n >= s->cap) { s->overflow = 1; return; }
        bs_arc_cell *c = &s->cells[s->n];
        c->px = px; c->py = py; c->thx = thx; c->thy = thy;
        c->l = l; c->chi = chi; c->d_in = d_in;
        c->kind = (signed char)kind; c->sense = (signed char)sense;
        c->leg = (short)leg;
    }
    s->n++;
}

/* am.build_cells_fillet :: straight(): corner-pace cells within 16 m
 * after / 20 m before a vertex, trim cells between. */
static void straight(cell_sink *s, double p0x, double p0y, double thx,
                     double thy, double L, int leg, double chi_at, int n_v,
                     double corner_pace, double vtrim)
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
        /* int(round(Lp / (pace Ts))): Python rounds half to even */
        int n = (int)nearbyint(Lp[p] / (pace[p] * BS_TS));
        if (n < 1) n = 1;
        const double l = Lp[p] / n;
        for (int j = 0; j < n; ++j) {
            add_cell(s, p0x + thx * sacc, p0y + thy * sacc, thx, thy, l,
                     kind[p], (first && j == 0) ? chi_at : 0.0, leg, 0.0, 0);
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
                      cell_sink *s, bs_arc_chart *out)
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
        rf[k] = w / (1.0 - cos(h));
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
        straight(s, sx, sy, t0, t1, L, i, chi_carry, n_v, corner_pace, vtrim);
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
                add_cell(s, p0x, p0y, thx, thy, l, BS_ARC_FILLET, dchi, i,
                         (dA < dB) ? dA : dB, (int)sgn);
                prevx = thx; prevy = thy;
            }
            chi_carry = wrap180(deg_of(atan2(tBy, tBx) - atan2(prevy, prevx)));
        } else {
            chi_carry = 0.0;
        }
        if (out != NULL) out->leg_last[i] = s->n - 1;
    }
    const int k_end = s->n;
    {   /* the trailing straight: tail / round(tail / (V_TRIM Ts)) cells */
        const double thx = (s->cells != NULL && s->n > 0) ? s->cells[s->n - 1].thx : tx[n_v - 2];
        const double thy = (s->cells != NULL && s->n > 0) ? s->cells[s->n - 1].thy : ty[n_v - 2];
        int n = (int)nearbyint(BS_ARC_TAIL / (vtrim * BS_TS));
        if (n < 1) n = 1;
        const double l = BS_ARC_TAIL / n;
        for (int j = 0; j < n; ++j)
            add_cell(s, vx[n_v - 1] + thx * (j * l), vy[n_v - 1] + thy * (j * l),
                     thx, thy, l, BS_ARC_TA, 0.0, n_v - 2, 0.0, 0);
    }
    if (out != NULL) {
        out->k_end = k_end;
        out->n_cells = s->n;
        out->n_v = n_v;
        out->v_trim = vtrim;
        out->v_leg = v_leg();
        out->corner_pace = corner_pace;
        out->w = w;
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

int bs_arc_chart_build(const double *vx, const double *vy, int n_wp,
                       void *mem, size_t mem_len, bs_arc_chart *out)
{
    if (n_wp < 2 || n_wp > BS_ARC_MAX_WP) return -1;
    memset(out, 0, sizeof(*out));
    for (int i = 0; i < n_wp; ++i) { out->vx[i] = vx[i]; out->vy[i] = vy[i]; }
    cell_sink cnt = { NULL, 0, 0, 0 };
    int rc = chart_pass(vx, vy, n_wp, &cnt, NULL);
    if (rc != 0) return rc;
    if ((size_t)cnt.n * sizeof(bs_arc_cell) > mem_len) return -2;
    cell_sink fill = { (bs_arc_cell *)mem, 0, cnt.n, 0 };
    rc = chart_pass(vx, vy, n_wp, &fill, out);
    if (rc != 0) return rc;
    if (fill.overflow || fill.n != cnt.n) return -2;
    out->cells = (bs_arc_cell *)mem;
    return 0;
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
    if (c->kind == BS_ARC_FILLET) return 2 + k;
    return (c->kind == BS_ARC_TA) ? 0 : 1;
}

/* ha.rows_corner_m layout: (speed_a, accel_a) pairs for a = 0..19, then
 * +e_l, +e_c, +j_l, +j_c, -e_l, -e_c, -j_l, -j_c, and (v1 only) +nu, -nu.
 * Every family of the fillet chart has all NR rows (the disc row with
 * margin < 1e-9 would need pace >= VLEG, which no cell has). */
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
    /* the disc is about the FAMILY's nominal pace: V_TRIM for ta, the
     * corner pace for ca (ac.register_families), the cell's own pace only
     * on a fillet chord (am.fillet_rows) */
    const double pace = (c->kind == BS_ARC_FILLET) ? c->l / BS_TS
                      : ((c->kind == BS_ARC_TA) ? ch->v_trim : ch->corner_pace);
    const double VLEG = ch->v_leg;
    for (int a = 0; a < 20; ++a) {
        const double th = rad_of(18.0 * a);
        m[2 * a] = VLEG - pace * cos(th);
        m[2 * a + 1] = (c->kind == BS_ARC_TA) ? bs_rc_aniso[a] : BS_ARC_ACCC;
    }
    m[40] = BS_ARC_BAND; m[44] = BS_ARC_BAND;
    if (c->kind == BS_ARC_FILLET) {
        /* am.fillet_rows: s e_c <= max(w - d_in, 0.05), -s e_c <= w + d_in */
        double in = ch->w - c->d_in;
        if (in < 0.05) in = 0.05;
        const double outm = ch->w + c->d_in;
        if (c->sense > 0) { m[41] = in; m[45] = outm; }
        else              { m[41] = outm; m[45] = in; }
    } else {
        m[41] = ch->w; m[45] = ch->w;
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
    if (arc_dare(S->rows, m, P6, K36) != 0) { S->dare_fail++; return -1; }
    const int slot = (S->dare_n < BS_ARC_DARE_CACHE) ? S->dare_n++ : S->dare_next;
    S->dare_next = (slot + 1) % BS_ARC_DARE_CACHE;
    if (S->dare_next < 2) S->dare_next = 2;    /* keep ta / ca resident */
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
 * start), hold or advance through every crossed boundary. */
int bs_arc_transition(const bs_arc_chart *ch, int k, const double *y,
                      double *T, double *c)
{
    double Tl[NX * NX], cl[NX];
    for (int i = 0; i < NX * NX; ++i) Tl[i] = (i % (NX + 1) == 0) ? 1.0 : 0.0;
    for (int i = 0; i < NX; ++i) cl[i] = 0.0;
    const double pace_k = ch->cells[k].l / BS_TS;
    if (k + 1 >= ch->n_cells || !crossed(ch, k + 1, y, pace_k)) {
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
        if (kn + 1 < ch->n_cells && crossed(ch, kn + 1, st2, pace_new)) {
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
    for (int t = 0; t < NN; ++t) {
        bs_arc_plant(xi, &U[NU * t], y);
        const int kn = bs_arc_transition(ch, k, y, T, c);
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
    memset(S, 0, sizeof(*S));
    template_rows(S->rows);
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

int bs_arc_engage(bs_arc_driver *d)
{
    const bs_arc_chart *ch = d->chart;
    d->k = 0;
    d->kw = 0;
    memset(d->xi, 0, sizeof(d->xi));
    d->xi[0] = -(ch->cells[0].l / BS_TS);          /* at rest on vertex 0 */
    memset(d->U, 0, sizeof(d->U));
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

int bs_arc_tick(bs_arc_driver *d, double *u0)
{
    const bs_arc_chart *ch = d->chart;
    bs_arc_sched *S = d->S;
    if (!d->engaged) return -1;
    if (d->k >= ch->k_end) { d->finished = 1; return 1; }
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
        bs_newton_stats st = {0, 0, 0, 0};
        double lam = 0.0;
        if (solve_pinned(d, MM, d->q, &st, &lam) != 0) return -3;
        if (r == 0) d->lambda0 = lam;
        d->bt += st.backtracks;
        d->armijo_fail += st.armijo_fail;
        d->no_step += st.no_step;
        d->factor_fail += st.factor_fail;
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
