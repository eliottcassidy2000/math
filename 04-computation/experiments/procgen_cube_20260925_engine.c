/*
 * procgen_cube_20260925_engine.c -- fast engine for the strategy cube (session
 * collatz-procgen-20260922, cube lane, 2026-09-25).  Cross-checked against the exact
 * pure-Python reference procgen_cube_20260925_core.py by procgen_cube_20260925_run.py.
 *
 * A level-k sign strategy sigma: odd residues mod 2^k -> {+1,-1};
 *   T(n) = n/2 (n even), (3n + sigma(n mod 2^k))/2 (n odd).
 * Mask encoding (k <= 7): bit i <-> residue 2i+1, bit 1 = minus.
 *
 * Per strategy it computes
 *   - the parity-transition graph G on Z/2^k (s -> two lifts of T(s) mod 2^(k-1)),
 *     node weight w(s) = log(3/2) (odd) / -log 2 (even);
 *   - SCCs (Tarjan), bottom SCCs, Karp max / min cycle means (global and per bottom SCC;
 *     only for k <= KARPMAX, otherwise signs by Bellman-Ford);
 *   - stationary law on each bottom SCC (Gauss for K <= 64, else power iteration),
 *     drift = pi_odd log 3 - log 2, and the Haar basin measure of each bottom SCC;
 *   - |Bad_L| for L = 1..LMAX by the Terras DP over (node, a) (exact uint64 counts), and the
 *     first L with |Bad_L| = 0 (searched to LZMAX when the graph has no expanding cycle);
 *   - exceptional dimension D = inf_{beta>=0} log2 rho(A(beta)), A(beta)_{s,t} = e^{beta w(s)};
 *   - cycle search over positive starts 1..N (memoised fates; value cap CAP; Brent detection).
 *
 * Modes:
 *   engine all   k N          -- every mask at level k (k <= 5 sensible)
 *   engine list  k N          -- masks (hex) read from stdin, one per line
 *   engine rand  k N count seed  -- random strategies at level k (sigma bits from splitmix64)
 *   engine rsign N count seed    -- k = infinity model: independent random sign per odd integer;
 *                                    reports cycles among starts 1..N only
 *   engine ball  k r which       -- exhaustive Hamming ball around Collatz (k <= 7): which=1 class (i), 2 class (ii)
 *   engine ballr k r which       -- the same restricted to env BALL_ALLOWED flips (exploratory; upper bounds only)
 *   engine expcyc k M            -- expanding simple cycles of G_sigma (for the lazy MaxSAT of the boundary script)
 *   engine cycsub k M sign       -- cycles of a given sign inside a node subset (for the (ii) boundary)
 * Output: one line per strategy, key=value fields (parsed by the Python driver).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define LMAX 40
#define LZMAX 400
#define KARPMAX 10
#define CAP ((int64_t)1 << 60)
#define STEPMAX 200000
#define MAXCYC 64

static int k, K, H;
static int8_t *sig;          /* sig[r] for r mod K (odd r) */
static int *s0, *s1;          /* successors */
static double LG2, LG3, WODD, WEVEN;

static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void build_graph(void) {
    for (int s = 0; s < K; s++) {
        int t;
        if (s % 2 == 0) t = s / 2;
        else t = (3 * s + sig[s]) / 2;
        int t0 = H ? t % H : 0;
        s0[s] = t0;
        s1[s] = t0 + H;
    }
}
static inline double wt(int s) { return (s & 1) ? WODD : WEVEN; }

/* ---------------------------------------------------------------- Tarjan (iterative) */
static int *comp_of, ncomp;
static void tarjan(void) {
    int *idx = malloc(K * sizeof(int)), *low = malloc(K * sizeof(int));
    int *onst = calloc(K, sizeof(int)), *st = malloc(K * sizeof(int)), sp = 0;
    int *wv = malloc(K * sizeof(int)), *we = malloc(K * sizeof(int)), wp;
    int counter = 0;
    for (int i = 0; i < K; i++) idx[i] = -1;
    ncomp = 0;
    for (int r = 0; r < K; r++) {
        if (idx[r] >= 0) continue;
        wp = 0; wv[wp] = r; we[wp] = 0; wp++;
        idx[r] = low[r] = counter++; st[sp++] = r; onst[r] = 1;
        while (wp) {
            int v = wv[wp - 1];
            if (we[wp - 1] < 2) {
                int w = we[wp - 1] == 0 ? s0[v] : s1[v];
                we[wp - 1]++;
                if (idx[w] < 0) {
                    idx[w] = low[w] = counter++; st[sp++] = w; onst[w] = 1;
                    wv[wp] = w; we[wp] = 0; wp++;
                } else if (onst[w] && idx[w] < low[v]) low[v] = idx[w];
            } else {
                wp--;
                if (wp) { int u = wv[wp - 1]; if (low[v] < low[u]) low[u] = low[v]; }
                if (low[v] == idx[v]) {
                    int w;
                    do { w = st[--sp]; onst[w] = 0; comp_of[w] = ncomp; } while (w != v);
                    ncomp++;
                }
            }
        }
    }
    free(idx); free(low); free(onst); free(st); free(wv); free(we);
}

/* ---------------------------------------------------------------- Karp on a node subset */
/* returns max cycle mean over cycles inside mask (sign = +1) or min (sign = -1); -INF if acyclic */
static double karp(const char *in, int sign) {
    int n = 0;
    for (int v = 0; v < K; v++) if (in[v]) n++;
    if (n == 0) return -INFINITY * sign;
    double *D = malloc((size_t)(n + 1) * K * sizeof(double));
    const double NEG = -1e300;
    for (int v = 0; v < K; v++) D[v] = in[v] ? 0.0 : NEG;
    for (int j = 1; j <= n; j++) {
        double *prev = D + (size_t)(j - 1) * K, *cur = D + (size_t)j * K;
        for (int v = 0; v < K; v++) cur[v] = NEG;
        for (int u = 0; u < K; u++) {
            if (!in[u] || prev[u] <= NEG / 2) continue;
            double val = prev[u] + sign * wt(u);
            int a = s0[u], b = s1[u];
            if (in[a] && val > cur[a]) cur[a] = val;
            if (in[b] && val > cur[b]) cur[b] = val;
        }
    }
    double best = NEG;
    double *Dn = D + (size_t)n * K;
    for (int v = 0; v < K; v++) {
        if (!in[v] || Dn[v] <= NEG / 2) continue;
        double worst = 1e300;
        for (int j = 0; j < n; j++) {
            double dj = D[(size_t)j * K + v];
            if (dj <= NEG / 2) continue;
            double q = (Dn[v] - dj) / (n - j);
            if (q < worst) worst = q;
        }
        if (worst < 1e299 && worst > best) best = worst;
    }
    free(D);
    if (best <= NEG / 2) return -INFINITY * sign;
    return sign * best;
}

/* Bellman-Ford: does the subset contain a cycle with sign*weight > 0 ? */
static int bf_poscycle(const char *in, int sign) {
    double *d = malloc(K * sizeof(double));
    int n = 0;
    for (int v = 0; v < K; v++) { d[v] = 0.0; if (in[v]) n++; }
    int changed = 1;
    for (int it = 0; it <= n + 1 && changed; it++) {
        changed = 0;
        for (int u = 0; u < K; u++) {
            if (!in[u]) continue;
            double val = d[u] + sign * wt(u);
            int a = s0[u], b = s1[u];
            if (in[a] && val > d[a] + 1e-12) { d[a] = val; changed = 1; }
            if (in[b] && val > d[b] + 1e-12) { d[b] = val; changed = 1; }
        }
    }
    free(d);
    return changed;
}

/* ---------------------------------------------------------------- stationary law */
static void stationary(const char *in, double *pi) {
    int n = 0;
    int *ix = malloc(K * sizeof(int)), *nodes = malloc(K * sizeof(int));
    for (int v = 0; v < K; v++) { ix[v] = -1; if (in[v]) { ix[v] = n; nodes[n++] = v; } }
    if (n <= 64) {
        double *A = calloc((size_t)n * (n + 1), sizeof(double));
        for (int i = 0; i < n; i++) {
            int v = nodes[i];
            A[(size_t)ix[s0[v]] * (n + 1) + i] += 0.5;
            A[(size_t)ix[s1[v]] * (n + 1) + i] += 0.5;
            A[(size_t)i * (n + 1) + i] -= 1.0;
        }
        for (int i = 0; i < n; i++) A[(size_t)(n - 1) * (n + 1) + i] = 1.0;
        A[(size_t)(n - 1) * (n + 1) + n] = 1.0;
        for (int c = 0; c < n; c++) {
            int p = c;
            for (int r = c + 1; r < n; r++) if (fabs(A[(size_t)r * (n + 1) + c]) > fabs(A[(size_t)p * (n + 1) + c])) p = r;
            if (p != c) for (int j = 0; j <= n; j++) { double t = A[(size_t)c * (n + 1) + j]; A[(size_t)c * (n + 1) + j] = A[(size_t)p * (n + 1) + j]; A[(size_t)p * (n + 1) + j] = t; }
            double pv = A[(size_t)c * (n + 1) + c];
            for (int j = 0; j <= n; j++) A[(size_t)c * (n + 1) + j] /= pv;
            for (int r = 0; r < n; r++) if (r != c) {
                double f = A[(size_t)r * (n + 1) + c];
                if (f != 0) for (int j = 0; j <= n; j++) A[(size_t)r * (n + 1) + j] -= f * A[(size_t)c * (n + 1) + j];
            }
        }
        for (int v = 0; v < K; v++) pi[v] = 0;
        for (int i = 0; i < n; i++) pi[nodes[i]] = A[(size_t)i * (n + 1) + n];
        free(A);
    } else {
        double *p = malloc(K * sizeof(double)), *q = malloc(K * sizeof(double));
        for (int v = 0; v < K; v++) p[v] = in[v] ? 1.0 / n : 0.0;
        for (int it = 0; it < 6000; it++) {
            for (int v = 0; v < K; v++) q[v] = 0.5 * p[v];
            for (int v = 0; v < K; v++) if (in[v]) { q[s0[v]] += 0.25 * p[v]; q[s1[v]] += 0.25 * p[v]; }
            double *t = p; p = q; q = t;
        }
        for (int v = 0; v < K; v++) pi[v] = p[v];
        free(p); free(q);
    }
    free(ix); free(nodes);
}

/* Haar basin measure of each SCC id (only bottom ones get mass in the limit) */
static void basins(double *mass_by_comp) {
    double *p = malloc(K * sizeof(double)), *q = malloc(K * sizeof(double));
    for (int v = 0; v < K; v++) p[v] = 1.0 / K;
    for (int it = 0; it < 4000; it++) {
        for (int v = 0; v < K; v++) q[v] = 0.5 * p[v];
        for (int v = 0; v < K; v++) { q[s0[v]] += 0.25 * p[v]; q[s1[v]] += 0.25 * p[v]; }
        double *t = p; p = q; q = t;
    }
    for (int c = 0; c < ncomp; c++) mass_by_comp[c] = 0;
    for (int v = 0; v < K; v++) mass_by_comp[comp_of[v]] += p[v];
    free(p); free(q);
}

/* ---------------------------------------------------------------- exceptional counts */
static int amin_tab[LZMAX + 2];
static void init_amin(void) {
    for (int j = 0; j <= LZMAX + 1; j++) {
        /* least a with 3^a > 2^j ; exact via long double logs with a safety check */
        int a = (int)floorl((long double)j * logl(2.0L) / logl(3.0L)) + 1;
        amin_tab[j] = a;
    }
}
/* counts[L] for L = 0..Lcap ; returns first L with zero count (<= Lcap) or -1 */
static int bad_dp(int Lcap, uint64_t *counts_out) {
    int A = Lcap + 2;
    uint64_t *cur = calloc((size_t)K * A, sizeof(uint64_t)), *nxt = calloc((size_t)K * A, sizeof(uint64_t));
    for (int s = 0; s < K; s++) cur[(size_t)s * A + 0] = 1;
    if (counts_out) counts_out[0] = K;
    int zero = -1;
    for (int j = 0; j < Lcap; j++) {
        memset(nxt, 0, (size_t)K * A * sizeof(uint64_t));
        uint64_t tot = 0;
        for (int s = 0; s < K; s++) {
            int odd = s & 1;
            for (int a = 0; a <= j; a++) {
                uint64_t c = cur[(size_t)s * A + a];
                if (!c) continue;
                int a2 = a + odd;
                if (a2 >= amin_tab[j + 1]) {
                    nxt[(size_t)s0[s] * A + a2] += c;
                    nxt[(size_t)s1[s] * A + a2] += c;
                    tot += 2 * c;
                }
            }
        }
        uint64_t *t = cur; cur = nxt; nxt = t;
        if (counts_out) counts_out[j + 1] = tot;
        if (tot == 0) { zero = j + 1; break; }
    }
    free(cur); free(nxt);
    return zero;
}

/* existence version (no counts, no overflow) for the long search: first L <= Lcap with Bad_L empty */
static int bad_zero_bool(int Lcap) {
    int A = Lcap + 2;
    char *cur = calloc((size_t)K * A, 1), *nxt = calloc((size_t)K * A, 1);
    for (int s = 0; s < K; s++) cur[(size_t)s * A + 0] = 1;
    int zero = -1;
    for (int j = 0; j < Lcap; j++) {
        memset(nxt, 0, (size_t)K * A);
        int any = 0;
        for (int s = 0; s < K; s++) {
            int odd = s & 1;
            for (int a = 0; a <= j; a++) {
                if (!cur[(size_t)s * A + a]) continue;
                int a2 = a + odd;
                if (a2 >= amin_tab[j + 1]) { nxt[(size_t)s0[s] * A + a2] = 1; nxt[(size_t)s1[s] * A + a2] = 1; any = 1; }
            }
        }
        char *t = cur; cur = nxt; nxt = t;
        if (!any) { zero = j + 1; break; }
    }
    free(cur); free(nxt);
    return zero;
}

/* ---------------------------------------------------------------- dimension */
static double logrho(double beta) {
    /* rho(A(beta)) via power iteration on A + I (aperiodic), A_{s,t} = e^{beta w(s)} for s->t */
    double eo = exp(beta * WODD), ee = exp(beta * WEVEN);
    double *x = malloc(K * sizeof(double)), *y = malloc(K * sizeof(double));
    for (int v = 0; v < K; v++) x[v] = 1.0;
    double lam = 0;
    for (int it = 0; it < 300; it++) {
        for (int v = 0; v < K; v++) y[v] = x[v];
        for (int u = 0; u < K; u++) { double f = ((u & 1) ? eo : ee) * x[u]; y[s0[u]] += f; y[s1[u]] += f; }
        double m = 0; for (int v = 0; v < K; v++) if (y[v] > m) m = y[v];
        for (int v = 0; v < K; v++) x[v] = y[v] / m;
        lam = m;
    }
    free(x); free(y);
    return log(lam - 1.0);   /* rho(A) = rho(A+I) - 1 */
}
static double exc_dimension(void) {
    /* inf over beta >= 0 of log2 rho(A(beta)); golden section on [0, 40] (convex) */
    double a = 0, b = 40, g = (sqrt(5.0) - 1) / 2;
    double c = b - g * (b - a), d = a + g * (b - a);
    double fc = logrho(c), fd = logrho(d);
    for (int it = 0; it < 40; it++) {
        if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = logrho(c); }
        else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = logrho(d); }
    }
    double f0 = logrho(0.0);
    double best = fc < fd ? fc : fd;
    if (f0 < best) best = f0;
    return best / LG2;
}

/* ---------------------------------------------------------------- cycle search */
typedef struct { int64_t mn; int len, odd; int64_t basin; } cyc_t;
static int32_t *fate;   /* 0 unknown, >0 index+1 of cycle, -1 escape, -2 unresolved */
static int64_t *pathbuf;
static cyc_t cycs[MAXCYC];
static int ncyc, ncyc_total;

static inline int64_t Tstep(int64_t x) {
    if (!(x & 1)) return x >> 1;
    return (3 * x + sig[x & (K - 1)]) / 2;
}
static inline int64_t Tstep_rs(int64_t x, uint64_t seed) {
    if (!(x & 1)) return x >> 1;
    uint64_t h = seed ^ ((uint64_t)x * 0x9e3779b97f4a7c15ULL);
    uint64_t z = splitmix64(&h);
    int s = (z & 1) ? -1 : 1;
    return (3 * x + s) / 2;
}

static int find_or_add_cycle(int64_t x, int rs, uint64_t seed) {
    int64_t mn = x, y = x; int len = 0, odd = 0;
    do { if (y < mn) mn = y; odd += (int)(y & 1); y = rs ? Tstep_rs(y, seed) : Tstep(y); len++; } while (y != x);
    for (int i = 0; i < ncyc; i++) if (cycs[i].mn == mn) return i;
    ncyc_total++;
    if (ncyc < MAXCYC) { cycs[ncyc].mn = mn; cycs[ncyc].len = len; cycs[ncyc].odd = odd; cycs[ncyc].basin = 0; return ncyc++; }
    return -3;   /* overflow of the table (reported) */
}

static void cycle_search(int64_t N, int rs, uint64_t seed, int64_t *esc, int64_t *unres) {
    memset(fate, 0, (size_t)(N + 1) * sizeof(int32_t));
    ncyc = ncyc_total = 0; *esc = *unres = 0;
    for (int64_t n = 1; n <= N; n++) {
        if (fate[n]) { if (fate[n] > 0) cycs[fate[n] - 1].basin++; else if (fate[n] == -1) (*esc)++; else (*unres)++; continue; }
        int64_t x = n; int np = 0; int32_t res = 0;
        int64_t tort = x; long power = 1, lam = 0;
        while (1) {
            if (x <= N && fate[x]) { res = fate[x]; break; }
            if (x > CAP || x < 1) { res = -1; break; }
            if (np >= STEPMAX) { res = -2; break; }
            pathbuf[np++] = x;
            x = rs ? Tstep_rs(x, seed) : Tstep(x);
            lam++;
            if (x == tort) {   /* Brent: cycle through tort */
                int ci = find_or_add_cycle(x, rs, seed);
                res = ci >= 0 ? ci + 1 : -2;
                break;
            }
            if (power == lam) { tort = x; power *= 2; lam = 0; }
        }
        for (int i = 0; i < np; i++) if (pathbuf[i] <= N) fate[pathbuf[i]] = res;
        fate[n] = res;
        if (res > 0) cycs[res - 1].basin++; else if (res == -1) (*esc)++; else (*unres)++;
    }
}

/* ---------------------------------------------------------------- per-strategy analysis */
static void analyse(const char *label, int64_t N, int do_dim) {
    build_graph();
    tarjan();
    char *in = malloc(K);
    /* global expanding-cycle test */
    double mumax_glob, mumin_glob;
    memset(in, 1, K);
    if (k <= KARPMAX) { mumax_glob = karp(in, +1); mumin_glob = karp(in, -1); }
    else { mumax_glob = bf_poscycle(in, +1) ? 1.0 : -1.0; mumin_glob = bf_poscycle(in, -1) ? -1.0 : 1.0; }
    int cls_i = mumax_glob < 0;
    /* bottom SCCs */
    int *isbot = malloc(ncomp * sizeof(int)), *csize = calloc(ncomp, sizeof(int));
    for (int c = 0; c < ncomp; c++) isbot[c] = 1;
    for (int v = 0; v < K; v++) {
        csize[comp_of[v]]++;
        if (comp_of[s0[v]] != comp_of[v] || comp_of[s1[v]] != comp_of[v]) isbot[comp_of[v]] = 0;
    }
    double *mass = malloc(ncomp * sizeof(double));
    basins(mass);
    double *pi = malloc(K * sizeof(double));
    int cls_ii = 0, nbot = 0, npos = 0;
    char botbuf[8192]; int bl = 0; botbuf[0] = 0;
    for (int c = 0; c < ncomp; c++) {
        if (!isbot[c]) continue;
        nbot++;
        for (int v = 0; v < K; v++) in[v] = comp_of[v] == c;
        double mn, mx;
        if (k <= KARPMAX) { mn = karp(in, -1); mx = karp(in, +1); }
        else { mn = bf_poscycle(in, -1) ? -1.0 : 1.0; mx = bf_poscycle(in, +1) ? 1.0 : -1.0; }
        stationary(in, pi);
        double podd = 0; for (int v = 0; v < K; v++) if (in[v] && (v & 1)) podd += pi[v];
        double drift = podd * LG3 - LG2;
        if (mn > 0) cls_ii = 1;
        if (drift > 0) npos++;
        int minnode = -1; for (int v = 0; v < K; v++) if (in[v]) { minnode = v; break; }
        if (bl < 7800)
            bl += snprintf(botbuf + bl, sizeof(botbuf) - bl, "%s%d,%d,%.9f,%.9f,%.12f,%.9f,%.9f", nbot > 1 ? ";" : "",
                           minnode, csize[c], mn, mx, podd, drift, mass[c]);
    }
    /* exceptional counts */
    uint64_t counts[LZMAX + 2];
    int lz = bad_dp(LMAX, counts);
    if (lz < 0 && cls_i) { lz = bad_zero_bool(LZMAX); if (lz < 0) lz = -9; }
    double D = -1e9;   /* -1e9 = not computed */
    if (cls_i) D = -INFINITY;
    else if (npos > 0 && do_dim < 2) D = 1.0;          /* theory: a positive-drift bottom SCC gives D = 1 */
    else if (do_dim) D = exc_dimension();
    /* cycles */
    int64_t esc = 0, unres = 0;
    if (N > 0) cycle_search(N, 0, 0, &esc, &unres);
    printf("S %s k=%d prim=%d mumax=%.9f mumin=%.9f ci=%d cii=%d nscc=%d nbot=%d npos=%d bots=%s lz=%d",
           label, k, 1, mumax_glob, mumin_glob, cls_i, cls_ii, ncomp, nbot, npos, botbuf, lz);
    printf(" bad=");
    for (int L = 1; L <= LMAX; L++) {
        uint64_t v = (lz > 0 && L >= lz) ? 0 : counts[L];
        printf("%s%llu", L > 1 ? "," : "", (unsigned long long)v);
    }
    printf(" D=%.6f", D);
    if (N > 0) {
        printf(" N=%lld esc=%lld unres=%lld ncyc=%d cyc=", (long long)N, (long long)esc, (long long)unres, ncyc_total);
        for (int i = 0; i < ncyc; i++) printf("%s%lld:%d:%d:%lld", i ? ";" : "", (long long)cycs[i].mn, cycs[i].len, cycs[i].odd, (long long)cycs[i].basin);
    }
    printf("\n");
    free(in); free(isbot); free(csize); free(mass); free(pi);
}

static void setup(int kk, int64_t N) {
    k = kk; K = 1 << k; H = K >> 1;
    sig = calloc(K, 1); s0 = malloc(K * sizeof(int)); s1 = malloc(K * sizeof(int));
    comp_of = malloc(K * sizeof(int));
    if (N > 0) { fate = malloc((size_t)(N + 1) * sizeof(int32_t)); pathbuf = malloc((size_t)(STEPMAX + 2) * sizeof(int64_t)); }
}

int main(int argc, char **argv) {
    LG2 = log(2.0); LG3 = log(3.0); WODD = LG3 - LG2; WEVEN = -LG2;
    init_amin();
    if (argc < 2) { fprintf(stderr, "usage\n"); return 1; }
    const char *mode = argv[1];
    if (!strcmp(mode, "all") || !strcmp(mode, "list") || !strcmp(mode, "rand")) {
        int kk = atoi(argv[2]); int64_t N = atoll(argv[3]);
        int do_dim = 1;
        if (getenv("CUBE_DIM")) do_dim = atoi(getenv("CUBE_DIM"));   /* 0 none, 1 default, 2 always numeric */
        setup(kk, N);
        if (!strcmp(mode, "all")) {
            uint64_t nm = 1ULL << (K / 2);
            for (uint64_t m = 0; m < nm; m++) {
                for (int i = 0; i < K / 2; i++) sig[2 * i + 1] = ((m >> i) & 1) ? -1 : 1;
                char lab[64]; snprintf(lab, sizeof lab, "mask=%llx", (unsigned long long)m);
                analyse(lab, N, do_dim);
                fflush(stdout);
            }
        } else if (!strcmp(mode, "list")) {
            char line[70000];
            while (fgets(line, sizeof line, stdin)) {
                int len = (int)strlen(line); while (len && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = 0;
                if (!len) continue;
                /* hex string, least significant digit last; bit i <-> residue 2i+1 */
                for (int i = 0; i < K / 2; i++) {
                    int dig = i / 4, bit = i % 4;
                    int pos = len - 1 - dig;
                    int v = 0;
                    if (pos >= 0) { char ch = line[pos]; v = (ch >= '0' && ch <= '9') ? ch - '0' : (ch | 32) - 'a' + 10; }
                    sig[2 * i + 1] = ((v >> bit) & 1) ? -1 : 1;
                }
                char lab[80]; snprintf(lab, sizeof lab, "mask=%.60s", line);
                analyse(lab, N, do_dim);
                fflush(stdout);
            }
        } else {
            int count = atoi(argv[4]); uint64_t seed = strtoull(argv[5], 0, 10);
            if (!getenv("CUBE_DIM")) do_dim = (k <= 8);
            for (int c = 0; c < count; c++) {
                uint64_t st = seed * 1000003ULL + c;
                for (int i = 0; i < K / 2; i++) {
                    uint64_t z = splitmix64(&st);
                    sig[2 * i + 1] = (z & 1) ? -1 : 1;
                }
                /* label: hex of the first up-to-60 hex digits is not needed; store sample id */
                char lab[64]; snprintf(lab, sizeof lab, "sample=%d", c);
                /* also print the sign word compactly for k <= 8 */
                if (k <= 8) {
                    char w[200]; int p = 0;
                    for (int i = 0; i < K / 2; i++) w[p++] = sig[2 * i + 1] > 0 ? '+' : '-';
                    w[p] = 0;
                    printf("W sample=%d %s\n", c, w);
                }
                analyse(lab, N, do_dim);
                fflush(stdout);
            }
        }
    } else if (!strcmp(mode, "ball")) {
        /* engine ball k r which : all masks (k <= 7) at Hamming distance <= r from Collatz (mask 0);
           which = 1: print those in class (i) [prefilter sigma(1)=+, sigma(-1)=-];
           which = 2: print those in class (ii).  Prints "B mask=<hex> d=<dist>" per hit, then a summary. */
        int kk = atoi(argv[2]); int r = atoi(argv[3]); int which = atoi(argv[4]);
        setup(kk, 0);
        char *in = malloc(K);
        int Hn = K / 2;
        long long tested = 0, hits = 0;
        int best = -1;
        for (int d = 0; d <= r && d <= Hn; d++) {
            /* iterate over all d-subsets of {0..Hn-1} via Gosper's hack on uint64 */
            if (d == 0) { /* mask 0 */ }
            uint64_t m = d ? ((d == 64) ? ~0ULL : ((1ULL << d) - 1)) : 0;
            uint64_t lim = (Hn == 64) ? 0 : (1ULL << Hn);
            while (1) {
                int ok = 1;
                if (which == 1) {
                    if (m & 1ULL) ok = 0;                          /* sigma(1) must be + */
                    if (!((m >> (Hn - 1)) & 1ULL)) ok = 0;         /* sigma(-1) must be - */
                }
                if (ok) {
                    for (int i = 0; i < Hn; i++) sig[2 * i + 1] = ((m >> i) & 1ULL) ? -1 : 1;
                    build_graph();
                    tested++;
                    int hit = 0;
                    if (which == 1) {
                        memset(in, 1, K);
                        hit = (k <= KARPMAX ? karp(in, +1) < 0 : !bf_poscycle(in, +1));
                    } else {
                        tarjan();
                        for (int c = 0; c < ncomp && !hit; c++) {
                            int bottom = 1;
                            for (int v = 0; v < K; v++) if (comp_of[v] == c && (comp_of[s0[v]] != c || comp_of[s1[v]] != c)) bottom = 0;
                            if (!bottom) continue;
                            for (int v = 0; v < K; v++) in[v] = comp_of[v] == c;
                            if (k <= KARPMAX ? karp(in, -1) > 0 : !bf_poscycle(in, -1)) hit = 1;
                        }
                    }
                    if (hit) {
                        hits++;
                        if (best < 0) best = d;
                        printf("B mask=%llx d=%d\n", (unsigned long long)m, d);
                    }
                }
                if (d == 0) break;
                /* Gosper */
                uint64_t c = m & (~m + 1), rr = m + c;
                if (rr == 0) break;
                m = (((rr ^ m) >> 2) / c) | rr;
                if (lim && m >= lim) break;
            }
            if (best >= 0) { printf("BALL k=%d which=%d first_hit_distance=%d hits_at_that_distance=%lld tested=%lld\n", k, which, best, hits, tested); fflush(stdout); break; }
            printf("BALLD k=%d which=%d d=%d none tested=%lld\n", k, which, d, tested); fflush(stdout);
        }
        if (best < 0) printf("BALL k=%d which=%d none_within=%d tested=%lld\n", k, which, r, tested);
    } else if (!strcmp(mode, "ballr")) {
        /* engine ballr k r which : like ball, but flips (relative to Collatz) are restricted to the residues
           listed in env BALL_ALLOWED (comma-separated odd residues mod 2^k, at most 64 of them) and the residues
           in env BALL_FORCED are always flipped.  Prints "B flips=r1,r2,... d=<dist>" per hit (upper bounds only). */
        int kk = atoi(argv[2]); int r = atoi(argv[3]); int which = atoi(argv[4]);
        setup(kk, 0);
        char *in = malloc(K);
        int allowed[64], na = 0, forced[64], nf = 0;
        char *ea = getenv("BALL_ALLOWED"), *ef = getenv("BALL_FORCED");
        if (ea) { char *t = strdup(ea); for (char *q = strtok(t, ","); q && na < 64; q = strtok(0, ",")) allowed[na++] = atoi(q); }
        if (ef) { char *t = strdup(ef); for (char *q = strtok(t, ","); q && nf < 64; q = strtok(0, ",")) forced[nf++] = atoi(q); }
        long long tested = 0, hits = 0; int best = -1;
        for (int d = 0; d + nf <= r && d <= na; d++) {
            uint64_t m = d ? ((d == 64) ? ~0ULL : ((1ULL << d) - 1)) : 0;
            uint64_t lim = (na == 64) ? 0 : (1ULL << na);
            while (1) {
                for (int i = 0; i < K / 2; i++) sig[2 * i + 1] = 1;
                for (int i = 0; i < nf; i++) sig[forced[i]] = -1;
                for (int i = 0; i < na; i++) if ((m >> i) & 1ULL) sig[allowed[i]] = -sig[allowed[i]];
                build_graph(); tested++;
                int hit = 0;
                if (which == 1) { memset(in, 1, K); hit = (k <= KARPMAX ? karp(in, +1) < 0 : !bf_poscycle(in, +1)); }
                else {
                    tarjan();
                    for (int c = 0; c < ncomp && !hit; c++) {
                        int bottom = 1;
                        for (int v = 0; v < K; v++) if (comp_of[v] == c && (comp_of[s0[v]] != c || comp_of[s1[v]] != c)) bottom = 0;
                        if (!bottom) continue;
                        for (int v = 0; v < K; v++) in[v] = comp_of[v] == c;
                        if (k <= KARPMAX ? karp(in, -1) > 0 : !bf_poscycle(in, -1)) hit = 1;
                    }
                }
                if (hit) {
                    hits++; if (best < 0) best = d + nf;
                    printf("B flips=");
                    int first = 1;
                    for (int v = 1; v < K; v += 2) if (sig[v] < 0) { printf("%s%d", first ? "" : ",", v); first = 0; }
                    printf(" d=%d\n", d + nf);
                }
                if (d == 0) break;
                uint64_t c = m & (~m + 1), rr = m + c;
                if (rr == 0) break;
                m = (((rr ^ m) >> 2) / c) | rr;
                if (lim && m >= lim) break;
            }
            if (best >= 0) { printf("BALL k=%d which=%d restricted first_hit_distance=%d hits=%lld tested=%lld\n", k, which, best, hits, tested); fflush(stdout); break; }
            printf("BALLD k=%d which=%d restricted d=%d none tested=%lld\n", k, which, d + nf, tested); fflush(stdout);
        }
        if (best < 0) printf("BALL k=%d which=%d restricted none_within=%d tested=%lld\n", k, which, r, tested);
    } else if (!strcmp(mode, "expcyc")) {
        /* engine expcyc k M : for each hex mask on stdin, print up to M distinct expanding simple cycles of G_sigma
           (Bellman-Ford longest path with predecessors; after each cycle one of its odd-node edges is banned).
           Output per mask: "E n=<count>" then lines "C v1,v2,...". */
        int kk = atoi(argv[2]); int M = atoi(argv[3]);
        setup(kk, 0);
        double *d = malloc(K * sizeof(double)); int *pred = malloc(K * sizeof(int));
        char *ban0 = malloc(K), *ban1 = malloc(K);
        int *cyc = malloc((K + 1) * sizeof(int));
        char line[70000];
        while (fgets(line, sizeof line, stdin)) {
            int len = (int)strlen(line); while (len && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = 0;
            if (!len) continue;
            for (int i = 0; i < K / 2; i++) {
                int dig = i / 4, bit = i % 4, pos = len - 1 - dig, v = 0;
                if (pos >= 0) { char ch = line[pos]; v = (ch >= '0' && ch <= '9') ? ch - '0' : (ch | 32) - 'a' + 10; }
                sig[2 * i + 1] = ((v >> bit) & 1) ? -1 : 1;
            }
            build_graph();
            memset(ban0, 0, K); memset(ban1, 0, K);
            int nfound = 0;
            char *outbuf = malloc((size_t)M * (K * 6 + 8) + 64); int ol = 0;
            for (int rep = 0; rep < M; rep++) {
                for (int v = 0; v < K; v++) { d[v] = 0; pred[v] = -1; }
                int last = -1;
                for (int it = 0; it <= K; it++) {
                    last = -1;
                    for (int u = 0; u < K; u++) {
                        double val = d[u] + wt(u);
                        if (!ban0[u] && val > d[s0[u]] + 1e-12) { d[s0[u]] = val; pred[s0[u]] = u; last = s0[u]; }
                        if (!ban1[u] && val > d[s1[u]] + 1e-12) { d[s1[u]] = val; pred[s1[u]] = u; last = s1[u]; }
                    }
                    if (last < 0) break;
                }
                if (last < 0) break;
                int v = last;
                for (int i = 0; i < K; i++) v = pred[v];
                int n = 0; cyc[n++] = v; int u = pred[v];
                while (u != v && n <= K) { cyc[n++] = u; u = pred[u]; }
                /* reverse to forward order */
                for (int i = 0; i < n / 2; i++) { int t = cyc[i]; cyc[i] = cyc[n - 1 - i]; cyc[n - 1 - i] = t; }
                int a = 0; for (int i = 0; i < n; i++) a += cyc[i] & 1;
                if (!((long double)a * logl(3.0L) > (long double)n * logl(2.0L))) break;   /* not expanding: stop */
                nfound++;
                ol += sprintf(outbuf + ol, "C ");
                for (int i = 0; i < n; i++) ol += sprintf(outbuf + ol, "%s%d", i ? "," : "", cyc[i]);
                ol += sprintf(outbuf + ol, "\n");
                for (int i = 0; i < n; i++) if (cyc[i] & 1) { int nx = cyc[(i + 1) % n]; if (s0[cyc[i]] == nx) ban0[cyc[i]] = 1; else ban1[cyc[i]] = 1; break; }
            }
            printf("E n=%d\n%s", nfound, ol ? outbuf : "");
            free(outbuf);
            fflush(stdout);
        }
    } else if (!strcmp(mode, "cycsub")) {
        /* engine cycsub k M sign : queries of two lines on stdin: a hex mask, then a 0/1 string of length 2^k
           (the node subset Y).  Prints up to M distinct simple cycles inside Y whose weight has the given sign
           (+1: expanding, -1: contracting): "E n=<count>" then "C v1,v2,...". */
        int kk = atoi(argv[2]); int M = atoi(argv[3]); int sg = atoi(argv[4]);
        setup(kk, 0);
        double *d = malloc(K * sizeof(double)); int *pred = malloc(K * sizeof(int));
        char *ban0 = malloc(K), *ban1 = malloc(K), *in = malloc(K);
        int *cyc = malloc((K + 1) * sizeof(int));
        char *line = malloc(200000), *sub = malloc(200000);
        while (fgets(line, 200000, stdin) && fgets(sub, 200000, stdin)) {
            int len = (int)strlen(line); while (len && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = 0;
            for (int i = 0; i < K / 2; i++) {
                int dig = i / 4, bit = i % 4, pos = len - 1 - dig, v = 0;
                if (pos >= 0) { char ch = line[pos]; v = (ch >= '0' && ch <= '9') ? ch - '0' : (ch | 32) - 'a' + 10; }
                sig[2 * i + 1] = ((v >> bit) & 1) ? -1 : 1;
            }
            for (int v = 0; v < K; v++) in[v] = sub[v] == '1';
            build_graph();
            memset(ban0, 0, K); memset(ban1, 0, K);
            int nfound = 0;
            char *outbuf = malloc((size_t)M * (K * 6 + 8) + 64); int ol = 0;
            for (int rep = 0; rep < M; rep++) {
                for (int v = 0; v < K; v++) { d[v] = 0; pred[v] = -1; }
                int last = -1;
                for (int it = 0; it <= K; it++) {
                    last = -1;
                    for (int u = 0; u < K; u++) {
                        if (!in[u]) continue;
                        double val = d[u] + sg * wt(u);
                        if (!ban0[u] && in[s0[u]] && val > d[s0[u]] + 1e-12) { d[s0[u]] = val; pred[s0[u]] = u; last = s0[u]; }
                        if (!ban1[u] && in[s1[u]] && val > d[s1[u]] + 1e-12) { d[s1[u]] = val; pred[s1[u]] = u; last = s1[u]; }
                    }
                    if (last < 0) break;
                }
                if (last < 0) break;
                int v = last;
                for (int i = 0; i < K; i++) v = pred[v];
                int n = 0; cyc[n++] = v; int u = pred[v];
                while (u != v && n <= K) { cyc[n++] = u; u = pred[u]; }
                for (int i = 0; i < n / 2; i++) { int t = cyc[i]; cyc[i] = cyc[n - 1 - i]; cyc[n - 1 - i] = t; }
                int a = 0; for (int i = 0; i < n; i++) a += cyc[i] & 1;
                long double wsum = (long double)a * logl(3.0L) - (long double)n * logl(2.0L);
                if (!(sg * wsum > 0)) break;
                nfound++;
                ol += sprintf(outbuf + ol, "C ");
                for (int i = 0; i < n; i++) ol += sprintf(outbuf + ol, "%s%d", i ? "," : "", cyc[i]);
                ol += sprintf(outbuf + ol, "\n");
                /* ban an edge of the cycle: prefer one leaving an odd node, else the first */
                int bi = 0; for (int i = 0; i < n; i++) if (cyc[i] & 1) { bi = i; break; }
                int nx = cyc[(bi + 1) % n]; if (s0[cyc[bi]] == nx) ban0[cyc[bi]] = 1; else ban1[cyc[bi]] = 1;
            }
            printf("E n=%d\n%s", nfound, ol ? outbuf : "");
            free(outbuf);
            fflush(stdout);
        }
    } else if (!strcmp(mode, "rsign")) {
        int64_t N = atoll(argv[2]); int count = atoi(argv[3]); uint64_t seed = strtoull(argv[4], 0, 10);
        k = 1; K = 2; setup(1, N);
        for (int c = 0; c < count; c++) {
            uint64_t sd = seed * 7919ULL + (uint64_t)c * 104729ULL + 12345ULL;
            int64_t esc, unres;
            cycle_search(N, 1, sd, &esc, &unres);
            printf("R sample=%d N=%lld esc=%lld unres=%lld ncyc=%d cyc=", c, (long long)N, (long long)esc, (long long)unres, ncyc_total);
            for (int i = 0; i < ncyc; i++) printf("%s%lld:%d:%d:%lld", i ? ";" : "", (long long)cycs[i].mn, cycs[i].len, cycs[i].odd, (long long)cycs[i].basin);
            printf("\n");
            fflush(stdout);
        }
    }
    return 0;
}
