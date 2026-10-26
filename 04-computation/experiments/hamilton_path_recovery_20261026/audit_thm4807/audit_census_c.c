/* Auditor's C census for THM-4807 at n = 9 (all classes) and n = 10 (strong classes).
 *
 * Reads gentourng output (upper-triangle strings, '1' at (i,j), i<j, means i->j) on stdin.
 * For every tournament:
 *   - subset DP for first/last vertices of Hamilton paths of every T[S];
 *   - traced arcs;  strong components (own reachability routine);
 *   - cut positions computed from the realisable prefix sets (DP form of (e) step 2);
 *   - three decoders from the trace support / path data: endpoint (THM-4807 (e), DP form),
 *     two-step count (G1), deleted-vertex one-way reachability (Session70 G4);
 *   - for strong T: (a) both halves, (b) |C_k|>=3, HP u->v exists, no HP v->u, the recipe
 *     path uP' and the cycle uP'u in T_D for every set D of untraced arcs containing u->v,
 *     (c) on every arc, Busch Cor 4.8 consequence d+(u) >= d+(v), two-step min <= 1.
 * No code from the lane is used.  Compile: cc -O2 -o audit_census_c audit_census_c.c
 * Usage: gentourng -q 9 | ./audit_census_c 9      gentourng -q -c 10 | ./audit_census_c 10
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 12
typedef uint32_t u32;
typedef uint64_t u64;

static int n;
static u32 out_[MAXN], inn_[MAXN], tr[MAXN], trin[MAXN];
static u32 *F, *B, *G;
static u64 *CNT;

static long long n_tour = 0, n_strong = 0, n_untraced = 0, n_strong_with_untraced = 0;
static long long c_arcs = 0, c_exact = 0, fails = 0;
static long long dec_E = 0, dec_G1 = 0, dec_G4 = 0, dec_agree = 0, cut_ok = 0;
static long long intra_pairs = 0, b_flipsets = 0;
static long long twostep[16][16];   /* [true][false] */
static long long hpuv[64];          /* number of Hamilton paths u->v */
static long long dgap[16];          /* d+(u)-d+(v) */
static long long kdist[16], ckdist[16], c1dist[16], elig[16];
static long long untr_hist[64];   /* strong classes by number of untraced arcs */

#define FAIL(msg) do { fails++; if (fails <= 20) fprintf(stderr, "FAILURE %s (tournament #%lld)\n", msg, n_tour); } while (0)

static int popc(u32 x) { return __builtin_popcount(x); }

static void dp_tables(void) {
    u32 N = 1u << n;
    F[0] = B[0] = 0;
    for (int v = 0; v < n; v++) { F[1u << v] = 1u << v; B[1u << v] = 1u << v; }
    for (u32 S = 1; S < N; S++) {
        if ((S & (S - 1)) == 0) continue;
        u32 f = 0, b = 0, x = S;
        while (x) {
            int v = __builtin_ctz(x); u32 low = x & (~x + 1); x ^= low;
            u32 S2 = S ^ low;
            if (F[S2] & inn_[v]) f |= low;
            if (B[S2] & out_[v]) b |= low;
        }
        F[S] = f; B[S] = b;
    }
}

static void traced_arcs(void) {
    u32 full = (1u << n) - 1;
    memset(tr, 0, sizeof tr);
    for (u32 S = 1; S < full; S++) {
        u32 f = F[S]; if (!f) continue;
        u32 b = B[full ^ S]; if (!b) continue;
        while (f) { int u = __builtin_ctz(f); f &= f - 1; tr[u] |= out_[u] & b; }
    }
    memset(trin, 0, sizeof trin);
    for (int a = 0; a < n; a++) { u32 x = tr[a]; while (x) { int b = __builtin_ctz(x); x &= x - 1; trin[b] |= 1u << a; } }
}

static u32 reach(const u32 *adj, int v, u32 allowed) {
    u32 seen = 1u << v, fr = seen;
    while (fr) {
        u32 nx = 0, x = fr;
        while (x) { int w = __builtin_ctz(x); x &= x - 1; nx |= adj[w]; }
        nx &= allowed & ~seen; seen |= nx; fr = nx;
    }
    return seen;
}

/* ordered strong components of adj restricted to allowed; returns k; sets *chain */
static int sccs(const u32 *adj, u32 allowed, u32 *comps, int *chain) {
    u32 R[MAXN]; u32 x = allowed;
    while (x) { int v = __builtin_ctz(x); x &= x - 1; R[v] = reach(adj, v, allowed); }
    int k = 0; u32 assigned = 0; int rs[MAXN]; int rep[MAXN];
    x = allowed;
    while (x) {
        int v = __builtin_ctz(x); x &= x - 1;
        if ((assigned >> v) & 1) continue;
        u32 c = 0, y = R[v];
        while (y) { int w = __builtin_ctz(y); y &= y - 1; if ((R[w] >> v) & 1) c |= 1u << w; }
        comps[k] = c; rs[k] = popc(R[v]); rep[k] = v; assigned |= c; k++;
    }
    for (int i = 1; i < k; i++) {      /* insertion sort by reach size, descending */
        u32 c = comps[i]; int r = rs[i], p = rep[i], j = i - 1;
        while (j >= 0 && rs[j] < r) { comps[j + 1] = comps[j]; rs[j + 1] = rs[j]; rep[j + 1] = rep[j]; j--; }
        comps[j + 1] = c; rs[j + 1] = r; rep[j + 1] = p;
    }
    *chain = 1;
    for (int t = 0; t + 1 < k; t++)
        if (R[rep[t]] != (comps[t] | R[rep[t + 1]])) *chain = 0;
    return k;
}

/* ends of Hamilton paths of T[A] that start at s */
static u32 ends_from(int s, u32 A) {
    u32 m = A & ~(1u << s), sub = 0;
    G[1u << s] = 1u << s;
    for (;;) {
        sub = (sub - m) & m;           /* next subset of m in increasing order */
        if (sub == 0) break;
        u32 S = sub | (1u << s), g = 0, x = sub;
        while (x) { int v = __builtin_ctz(x); x &= x - 1; if (G[S ^ (1u << v)] & inn_[v]) g |= 1u << v; }
        G[S] = g;
    }
    return G[A];
}

static u64 count_paths(int s, int t, u32 A) {   /* # Hamilton paths of T[A] from s to t */
    u32 m = A & ~(1u << s), sub = 0;
    u32 N = 1u << n;
    for (u32 S = 0; S < N; S++) for (int v = 0; v < n; v++) CNT[(u64)S * n + v] = 0;
    CNT[(u64)(1u << s) * n + s] = 1;
    for (;;) {
        sub = (sub - m) & m;
        if (sub == 0) break;
        u32 S = sub | (1u << s), x = sub;
        while (x) {
            int v = __builtin_ctz(x); x &= x - 1;
            u32 S2 = S ^ (1u << v), y = S2 & inn_[v];
            u64 c = 0;
            while (y) { int w = __builtin_ctz(y); y &= y - 1; c += CNT[(u64)S2 * n + w]; }
            CNT[(u64)S * n + v] = c;
        }
    }
    return CNT[(u64)A * n + t];
}

static int path_start(u32 S, int start, int *buf) {
    int len = 0; buf[len++] = start; u32 rem = S & ~(1u << start); int cur = start;
    while (rem) { u32 cand = out_[cur] & B[rem]; if (!cand) return -1; int w = __builtin_ctz(cand); buf[len++] = w; rem &= ~(1u << w); cur = w; }
    return len;
}
static int path_end(u32 S, int end, int *buf) {
    int tmp[MAXN], len = 0; tmp[len++] = end; u32 rem = S & ~(1u << end); int cur = end;
    while (rem) { u32 cand = inn_[cur] & F[rem]; if (!cand) return -1; int w = __builtin_ctz(cand); tmp[len++] = w; rem &= ~(1u << w); cur = w; }
    for (int i = 0; i < len; i++) buf[i] = tmp[len - 1 - i];
    return len;
}

static int legal(const u32 *adj, const int *w, int len) {
    for (int i = 0; i + 1 < len; i++) if (!((adj[w[i]] >> w[i + 1]) & 1)) return 0;
    return 1;
}

static void analyse(void) {
    u32 full = (1u << n) - 1;
    for (int i = 0; i < n; i++) inn_[i] = 0;
    for (int i = 0; i < n; i++) { u32 x = out_[i]; while (x) { int j = __builtin_ctz(x); x &= x - 1; inn_[j] |= 1u << i; } }
    dp_tables();
    traced_arcs();
    u32 comps[MAXN]; int chain;
    int k = sccs(out_, full, comps, &chain);
    if (!chain) FAIL("tournament condensation not a chain");
    /* cut positions from realisable prefix sets */
    {
        int cnt[MAXN + 1]; memset(cnt, 0, sizeof cnt);
        cnt[0] = 1; cnt[n] = 1;
        for (u32 S = 1; S < full; S++) {
            u32 f = F[S]; if (!f) continue;
            u32 b = B[full ^ S]; if (!b) continue;
            int ok = 0;
            while (f) { int x = __builtin_ctz(f); f &= f - 1; if (out_[x] & b) { ok = 1; break; } }
            if (ok) cnt[popc(S)]++;
        }
        int bounds[MAXN + 1], nb = 0, acc = 0; bounds[nb++] = 0;
        for (int t = 0; t < k; t++) { acc += popc(comps[t]); bounds[nb++] = acc; }
        int cuts[MAXN + 1], nc = 0;
        for (int m = 0; m <= n; m++) if (cnt[m] == 1) cuts[nc++] = m;
        int same = (nc == nb);
        for (int i = 0; same && i < nc; i++) if (cuts[i] != bounds[i]) same = 0;
        if (same) cut_ok++; else FAIL("cut positions != component boundaries");
    }
    /* decoders */
    u32 dcomps[MAXN]; int dchain;
    int dk = sccs(tr, full, dcomps, &dchain);
    if (!dchain || dk != k) FAIL("trace-support components differ");
    for (int t = 0; t < k && t < dk; t++) if (dcomps[t] != comps[t]) FAIL("trace-support component mismatch");
    int idx[MAXN];
    for (int t = 0; t < dk; t++) { u32 x = dcomps[t]; while (x) { int v = __builtin_ctz(x); x &= x - 1; idx[v] = t; } }
    u32 dE[MAXN], d1[MAXN], d4[MAXN];
    int okE = 1, ok1 = 1, ok4 = 1;
    for (int x = 0; x < n; x++) {
        u32 fw = 0;
        for (int y = 0; y < n; y++) if (y != x && idx[x] < idx[y]) fw |= 1u << y;
        dE[x] = d1[x] = d4[x] = tr[x] | fw;
    }
    for (int t = 0; t < dk; t++) {
        u32 C = dcomps[t];
        for (int x = 0; x < n; x++) {
            if (!((C >> x) & 1)) continue;
            for (int y = x + 1; y < n; y++) {
                if (!((C >> y) & 1)) continue;
                if (((tr[x] >> y) & 1) || ((tr[y] >> x) & 1)) continue;
                intra_pairs++;
                int truexy = (out_[x] >> y) & 1;
                /* G1 */
                int cxy = popc(tr[x] & trin[y] & C), cyx = popc(tr[y] & trin[x] & C);
                int tc = truexy ? cxy : cyx, fc = truexy ? cyx : cxy;
                twostep[tc < 15 ? tc : 15][fc < 15 ? fc : 15]++;
                if (!(tc >= 2 && fc <= 1)) FAIL("two-step rule");
                if ((cxy >= 2) == (cyx >= 2)) ok1 = 0;
                else if (cxy >= 2) d1[x] |= 1u << y; else d1[y] |= 1u << x;
                /* G4 */
                int dirs = 0, ne = 0;  /* bit0: x->y seen, bit1: y->x seen */
                for (int v = 0; v < n; v++) {
                    if (!((C >> v) & 1) || v == x || v == y) continue;
                    u32 A = C & ~(1u << v);
                    int xy = (reach(tr, x, A) >> y) & 1, yx = (reach(tr, y, A) >> x) & 1;
                    if (xy != yx) { ne++; dirs |= xy ? 1 : 2; }
                }
                elig[ne < 15 ? ne : 15]++;
                if (dirs == 1) d4[x] |= 1u << y; else if (dirs == 2) d4[y] |= 1u << x; else ok4 = 0;
                /* endpoint decoder (DP form): Hamilton path of T[C] from x to y / from y to x */
                int exy = (ends_from(x, C) >> y) & 1, eyx = (ends_from(y, C) >> x) & 1;
                if (exy == eyx) okE = 0;
                else if (exy) dE[x] |= 1u << y; else dE[y] |= 1u << x;
            }
        }
    }
    int eqE = okE, eq1 = ok1, eq4 = ok4, agree = 1;
    for (int x = 0; x < n; x++) {
        if (dE[x] != out_[x]) eqE = 0;
        if (d1[x] != out_[x]) eq1 = 0;
        if (d4[x] != out_[x]) eq4 = 0;
        if (dE[x] != d1[x] || d1[x] != d4[x]) agree = 0;
    }
    dec_E += eqE; dec_G1 += eq1; dec_G4 += eq4; dec_agree += agree;
    if (!eqE) FAIL("endpoint decoder");
    if (!eq1) FAIL("G1 decoder");
    if (!eq4) FAIL("G4 decoder");

    if (k != 1) return;
    n_strong++;
    /* untraced arcs of strong T */
    int ua[64][2], nu = 0;
    for (int u = 0; u < n; u++) { u32 x = out_[u] & ~tr[u]; while (x) { int v = __builtin_ctz(x); x &= x - 1; ua[nu][0] = u; ua[nu][1] = v; nu++; } }
    n_untraced += nu; if (nu) n_strong_with_untraced++;
    untr_hist[nu < 63 ? nu : 63]++;
    for (int a = 0; a < nu; a++) {
        int u = ua[a][0], v = ua[a][1];
        u32 cu[MAXN], cv[MAXN]; int ch1, ch2;
        int ku = sccs(out_, full & ~(1u << u), cu, &ch1);
        int kv = sccs(out_, full & ~(1u << v), cv, &ch2);
        if (ku < 2 || !((cu[ku - 1] >> v) & 1)) FAIL("(a) tail half");
        if (kv < 2 || !((cv[0] >> u) & 1)) FAIL("(a) dual half");
        kdist[ku < 15 ? ku : 15]++;
        ckdist[popc(cu[ku - 1])]++; c1dist[popc(cv[0])]++;
        if (popc(cu[ku - 1]) < 3 || popc(cv[0]) < 3) FAIL("(b) |C_k| >= 3");
        if (!((ends_from(u, full) >> v) & 1)) FAIL("(b) no Hamilton path u->v");
        if ((ends_from(v, full) >> u) & 1) FAIL("(b) Hamilton path v->u exists");
        u64 c = count_paths(u, v, full); hpuv[c < 63 ? c : 63]++;
        /* recipe path P' */
        u32 cset = out_[u] & cu[0];
        if (!cset) { FAIL("(b) no c"); continue; }
        int cc = __builtin_ctz(cset);
        int w[MAXN + 1], len = 0, buf[MAXN];
        w[len++] = u;
        for (int t = 0; t < ku; t++) {
            int l = (t == 0) ? path_start(cu[t], cc, buf) : (t == ku - 1) ? path_end(cu[t], v, buf)
                    : path_start(cu[t], __builtin_ctz(B[cu[t]]), buf);
            if (l < 0) { FAIL("(b) recipe path"); l = 0; }
            for (int i = 0; i < l; i++) w[len++] = buf[i];
        }
        if (len != n || !legal(out_, w, len) || w[n - 1] != v) FAIL("(b) uP' not a Hamilton path u->v");
        /* every D: set of untraced arcs containing u->v */
        for (u32 sub = 0; sub < (1u << nu); sub++) {
            if (!((sub >> a) & 1)) continue;
            u32 TD[MAXN]; memcpy(TD, out_, sizeof(u32) * n);
            for (int b2 = 0; b2 < nu; b2++) if ((sub >> b2) & 1) {
                int p = ua[b2][0], q = ua[b2][1];
                TD[p] &= ~(1u << q); TD[q] |= 1u << p;
            }
            int cyc[MAXN + 1];
            for (int i = 0; i < n; i++) cyc[i] = w[i];
            cyc[n] = u;
            if (!legal(TD, cyc, n + 1)) FAIL("(b) uP'u not a Hamiltonian cycle of T_D");
            if (legal(out_, cyc + 1, n)) FAIL("(b) P'u legal in T");
            b_flipsets++;
        }
        int gap = popc(out_[u]) - popc(out_[v]);
        if (gap < 0) FAIL("d+(u) < d+(v)");
        dgap[gap >= 0 ? (gap < 15 ? gap : 15) : 0]++;
    }
    /* (c) on every arc */
    for (int u = 0; u < n; u++) {
        u32 cu[MAXN]; int ch;
        int ku = sccs(out_, full & ~(1u << u), cu, &ch);
        u32 x = out_[u];
        while (x) {
            int v = __builtin_ctz(x); x &= x - 1;
            c_arcs++;
            int pred = 0;
            if (ku >= 2 && ((cu[ku - 1] >> v) & 1)) {
                u32 Ck = cu[ku - 1];
                int cond2 = (cu[ku - 2] & ~out_[u]) == 0;
                int cond3 = 1;
                u32 rest = Ck & ~(1u << v), A = rest;
                while (A) {
                    if ((F[A] & inn_[u]) && ((B[Ck & ~A] >> v) & 1)) { cond3 = 0; break; }
                    A = (A - 1) & rest;
                }
                pred = cond2 && cond3;
            }
            int actual = !((tr[u] >> v) & 1);
            if (pred == actual) c_exact++; else FAIL("(c) characterisation");
        }
    }
}

int main(int argc, char **argv) {
    n = (argc > 1) ? atoi(argv[1]) : 9;
    if (n < 2 || n > MAXN) { fprintf(stderr, "bad n\n"); return 1; }
    F = malloc(sizeof(u32) << n); B = malloc(sizeof(u32) << n); G = malloc(sizeof(u32) << n);
    CNT = malloc(sizeof(u64) * ((size_t)1 << n) * n);
    char line[256];
    int m = n * (n - 1) / 2;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '>') continue;
        int L = (int)strlen(line);
        while (L > 0 && (line[L - 1] == '\n' || line[L - 1] == '\r' || line[L - 1] == ' ')) line[--L] = 0;
        if (L == 0) continue;
        if (L != m) { fprintf(stderr, "bad line length %d\n", L); return 1; }
        memset(out_, 0, sizeof out_);
        int t = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (line[t] == '1') out_[i] |= 1u << j; else out_[j] |= 1u << i;
            t++;
        }
        n_tour++;
        analyse();
    }
    printf("n = %d: tournaments read %lld, strong %lld, untraced arcs in strong %lld (in %lld strong)\n",
           n, n_tour, n_strong, n_untraced, n_strong_with_untraced);
    printf("    cut positions = component boundaries: %lld/%lld\n", cut_ok, n_tour);
    printf("    decoders exact: endpoint %lld, two-step G1 %lld, deleted-reachability G4 %lld, all agree %lld\n",
           dec_E, dec_G1, dec_G4, dec_agree);
    printf("    (c) exact on %lld/%lld arcs of strong tournaments; (b) flip sets %lld; untraced intra-block pairs %lld\n",
           c_exact, c_arcs, b_flipsets, intra_pairs);
    printf("    two-step (true,false):");
    for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) if (twostep[a][b]) printf(" (%d,%d):%lld", a, b, twostep[a][b]);
    printf("\n    G4 eligible separators per pair:");
    for (int a = 0; a < 16; a++) if (elig[a]) printf(" %d:%lld", a, elig[a]);
    printf("\n    k=#comp(T-u):");
    for (int a = 0; a < 16; a++) if (kdist[a]) printf(" %d:%lld", a, kdist[a]);
    printf("\n    |C_k|:");
    for (int a = 0; a < 16; a++) if (ckdist[a]) printf(" %d:%lld", a, ckdist[a]);
    printf("\n    |C'_1|:");
    for (int a = 0; a < 16; a++) if (c1dist[a]) printf(" %d:%lld", a, c1dist[a]);
    printf("\n    #Hamilton paths u->v:");
    for (int a = 0; a < 64; a++) if (hpuv[a]) printf(" %d:%lld", a, hpuv[a]);
    printf("\n    d+(u)-d+(v):");
    for (int a = 0; a < 16; a++) if (dgap[a]) printf(" %d:%lld", a, dgap[a]);
    printf("\n    strong classes by number of untraced arcs:");
    for (int a = 0; a < 64; a++) if (untr_hist[a]) printf(" %d:%lld", a, untr_hist[a]);
    printf("  (Busch Thm 4.26 bound floor((n^2-4n+3)/8) = %d)", (n * n - 4 * n + 3) / 8);
    printf("\n    FAILURES: %lld\n", fails);
    return fails ? 2 : 0;
}
