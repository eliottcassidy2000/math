/* audit_labelled.c -- independent audit (S10 auditor) of THM-4868 by exhaustive enumeration of ALL labelled
 * tournaments on n <= 7 vertices (2^21 at n = 7).  No generation lemma, no structural assumption.
 *
 * Per labelled tournament: strong connectivity, non-critical set, Hamiltonian cycles counted exactly by DFS
 * (and re-counted by a subset DP whenever hc = 1).  Every labelled UH tournament is checked against THM-4868
 * (A)-(E) directly from the definitions; classes are formed with a canonical form that sorts vertices by score
 * and minimises the adjacency code by branch and bound (validated against full n! brute force at n <= 6).
 *
 * Side checks over ALL tournaments: Moon (strong, n >= 4 => >= 2 non-critical vertices); Las Vergnas as quoted
 * (strong with exactly two non-critical vertices => isomorphic to A_n); the "equivalently" form of (B)
 * (strong, min out-degree >= 2 => hc >= 2); HYP-9459's t_k(n), u_k(n) for k <= 5.
 *
 * Build: gcc -O2 -fopenmp -o audit_labelled audit_labelled.c ;  run: ./audit_labelled 7
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXN 8
static int N, NP;
static int PI_[32], PJ_[32];

static void build(uint32_t code, uint8_t *adj) {
    memset(adj, 0, MAXN);
    for (int b = 0; b < NP; b++) {
        if (code >> b & 1u) adj[PI_[b]] |= (uint8_t)(1u << PJ_[b]);
        else adj[PJ_[b]] |= (uint8_t)(1u << PI_[b]);
    }
}

/* strong connectivity of a general digraph (out-masks) on vertex set S; in-masks computed explicitly */
static int strong_dg(const uint8_t *out, int n, uint32_t S) {
    if (__builtin_popcount(S) <= 1) return 1;
    uint32_t in[MAXN] = {0};
    for (int x = 0; x < n; x++) if (S >> x & 1u) { uint32_t o = out[x] & S; while (o) { int y = __builtin_ctz(o); o &= o - 1; in[y] |= 1u << x; } }
    int s = __builtin_ctz(S);
    uint32_t seen = 1u << s, fr = seen;
    while (fr) { uint32_t nf = 0; for (uint32_t f = fr; f; f &= f - 1) nf |= out[__builtin_ctz(f)]; nf &= S & ~seen; seen |= nf; fr = nf; }
    if (seen != S) return 0;
    seen = 1u << s; fr = seen;
    while (fr) { uint32_t nf = 0; for (uint32_t f = fr; f; f &= f - 1) nf |= in[__builtin_ctz(f)]; nf &= S & ~seen; seen |= nf; fr = nf; }
    return seen == S;
}
static uint32_t reach_fwd(const uint8_t *out, uint32_t S, int s) {
    uint32_t seen = 1u << s, fr = seen;
    while (fr) { uint32_t nf = 0; for (uint32_t f = fr; f; f &= f - 1) nf |= out[__builtin_ctz(f)]; nf &= S & ~seen; seen |= nf; fr = nf; }
    return seen;
}

/* Hamiltonian cycles of T[S] by DFS from min(S); exact count; first cycle recorded (as a vertex sequence) */
typedef struct { const uint8_t *adj; uint32_t S; int s, target; long cnt; int path[MAXN], first[MAXN]; } dfs_t;
static void dfs(dfs_t *d, int v, uint32_t vis, int depth) {
    d->path[depth] = v;
    if (depth == d->target - 1) {
        if (d->adj[v] >> d->s & 1u) { if (d->cnt == 0) memcpy(d->first, d->path, sizeof d->path); d->cnt++; }
        return;
    }
    uint32_t nx = d->adj[v] & d->S & ~vis;
    while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dfs(d, w, vis | (1u << w), depth + 1); }
}
static long hc_dfs(const uint8_t *adj, uint32_t S, int *cyc) {
    dfs_t d; d.adj = adj; d.S = S; d.s = __builtin_ctz(S); d.target = __builtin_popcount(S); d.cnt = 0;
    if (d.target < 3) return 0;
    dfs(&d, d.s, 1u << d.s, 0);
    if (cyc && d.cnt) memcpy(cyc, d.first, sizeof(int) * d.target);
    return d.cnt;
}
/* second, independent count: subset DP */
static uint64_t hc_dp(const uint8_t *adj, int n, uint32_t S) {
    if (__builtin_popcount(S) < 3) return 0;
    int s = __builtin_ctz(S);
    uint64_t dp[1 << MAXN][MAXN];
    memset(dp, 0, sizeof(uint64_t) * ((size_t)1 << n) * MAXN);
    dp[1u << s][s] = 1;
    for (uint32_t m = 0; m < (1u << n); m++) {
        if ((m & ~S) || !(m >> s & 1u)) continue;
        for (int v = 0; v < n; v++) { uint64_t c = dp[m][v]; if (!c) continue;
            uint32_t nx = adj[v] & S & ~m;
            while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[m | (1u << w)][w] += c; } }
    }
    uint64_t t = 0;
    for (int v = 0; v < n; v++) if (v != s && (S >> v & 1u) && (adj[v] >> s & 1u)) t += dp[S][v];
    return t;
}

/* canonical form: positions sorted by out-degree, code = bits b(j,i) (j -> i) for i = 1..n-1, j = 0..i-1,
 * most significant first; minimised by branch and bound; cnt = number of minimising labellings = |Aut| */
typedef struct { const uint8_t *adj; int n, NPn; int deg[MAXN], slot[MAXN], perm[MAXN]; uint32_t used; uint64_t best; int have; long cnt; } canon_t;
static void canon_rec(canon_t *c, int i, uint64_t pc) {
    if (i == c->n) {
        if (!c->have || pc < c->best) { c->best = pc; c->have = 1; c->cnt = 1; }
        else if (pc == c->best) c->cnt++;
        return;
    }
    int L = (i + 1) * i / 2;
    for (int v = 0; v < c->n; v++) {
        if (c->used >> v & 1u) continue;
        if (c->deg[v] != c->slot[i]) continue;
        uint64_t bits = 0;
        for (int j = 0; j < i; j++) bits = (bits << 1) | ((c->adj[c->perm[j]] >> v) & 1u);
        uint64_t npc = (pc << i) | bits;
        if (c->have && npc > (c->best >> (c->NPn - L))) continue;
        c->perm[i] = v; c->used |= 1u << v;
        canon_rec(c, i + 1, npc);
        c->used &= ~(1u << v);
    }
}
static uint64_t canon_ss(const uint8_t *adj, int n, long *aut) {
    canon_t c; memset(&c, 0, sizeof c); c.adj = adj; c.n = n; c.NPn = n * (n - 1) / 2;
    for (int v = 0; v < n; v++) c.deg[v] = __builtin_popcount(adj[v]);
    int tmp[MAXN]; memcpy(tmp, c.deg, sizeof tmp);
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (tmp[b] < tmp[a]) { int t = tmp[a]; tmp[a] = tmp[b]; tmp[b] = t; }
    memcpy(c.slot, tmp, sizeof tmp);
    canon_rec(&c, 0, 0);
    if (aut) *aut = c.cnt;
    return c.best;
}
/* full n! brute force (validation only) */
static uint64_t canon_full(const uint8_t *adj, int n, long *aut) {
    int p[MAXN]; for (int i = 0; i < n; i++) p[i] = i;
    uint64_t best = UINT64_MAX; long cnt = 0;
    while (1) {
        uint64_t code = 0;
        for (int i = 1; i < n; i++) for (int j = 0; j < i; j++) code = (code << 1) | ((adj[p[j]] >> p[i]) & 1u);
        if (code < best) { best = code; cnt = 1; } else if (code == best) cnt++;
        int k = n - 2; while (k >= 0 && p[k] > p[k + 1]) k--;
        if (k < 0) break;
        int l = n - 1; while (p[l] < p[k]) l--;
        int t = p[k]; p[k] = p[l]; p[l] = t;
        for (int a = k + 1, b = n - 1; a < b; a++, b--) { t = p[a]; p[a] = p[b]; p[b] = t; }
    }
    if (aut) *aut = cnt;
    return best;
}

/* failure bits */
enum { F_A1, F_A2, F_B1, F_B2, F_C1, F_C2, F_C3, F_C4, F_C5, F_D1, F_D2, F_D3, F_L1, F_HC2, F_AUT, NFAIL };
static const char *fname[NFAIL] = { "A1_roundness", "A2_propagation", "B1_no_sole_exit_or_entry", "B2_gt2_or_nonadjacent",
    "C1_not_one_run_or_r_range", "C2_run_ends", "C3_succ_of_sole_exit_noncritical", "C4_unique_nc_sole_entry_exit",
    "C5_T_minus_y_not_UH_or_parent", "D1_cut_predicate", "D2_cut_plus_r", "D3_components_at_critical_w",
    "L1_lemma1_T_minus_v", "HC_dfs_vs_dp", "E_aut_nontrivial" };

static unsigned check_uh(const uint8_t *adj, int n, const int *cyc, int *o_r, int *o_k, int *o_se, int *o_sen) {
    unsigned f = 0; uint32_t F = (1u << n) - 1;
    int pos[MAXN], suc[MAXN], pre[MAXN], od[MAXN], id[MAXN];
    for (int i = 0; i < n; i++) pos[cyc[i]] = i;
    for (int v = 0; v < n; v++) { suc[v] = cyc[(pos[v] + 1) % n]; pre[v] = cyc[(pos[v] + n - 1) % n];
        od[v] = __builtin_popcount(adj[v]); id[v] = n - 1 - od[v]; }
    for (int v = 0; v < n; v++) if (!(adj[v] >> suc[v] & 1u)) f |= 1u << F_HC2;   /* sanity: cycle arcs exist */
    uint32_t NC = 0;
    for (int v = 0; v < n; v++) if (strong_dg(adj, n, F & ~(1u << v))) NC |= 1u << v;
    int r = __builtin_popcount(NC);
    /* (A1), (A2) */
    for (int v = 0; v < n; v++) if (NC >> v & 1u) {
        uint32_t want = 0; for (int s = 1; s <= od[v]; s++) want |= 1u << cyc[(pos[v] + s) % n];
        if (adj[v] != want) f |= 1u << F_A1;
        if (od[v] >= 2 && !(NC >> suc[v] & 1u)) f |= 1u << F_A2;
        if (id[v] >= 2 && !(NC >> pre[v] & 1u)) f |= 1u << F_A2;
    }
    /* (B) */
    int se = 0, sen = 0; int sev[MAXN], senv[MAXN];
    for (int v = 0; v < n; v++) { if (od[v] == 1) sev[se++] = v; if (id[v] == 1) senv[sen++] = v; }
    if (se < 1 || sen < 1) f |= 1u << F_B1;
    if (se > 2 || sen > 2) f |= 1u << F_B2;
    if (se == 2 && !(suc[sev[0]] == sev[1] || suc[sev[1]] == sev[0])) f |= 1u << F_B2;
    if (sen == 2 && !(suc[senv[0]] == senv[1] || suc[senv[1]] == senv[0])) f |= 1u << F_B2;
    /* (C) */
    int starts = 0, y = -1;
    for (int v = 0; v < n; v++) if ((NC >> v & 1u) && !(NC >> pre[v] & 1u)) { starts++; y = v; }
    if (starts != 1 || r < 2 || r > n - 1) f |= 1u << F_C1;
    int x = -1;
    if (starts == 1) {
        x = y; while (NC >> suc[x] & 1u) x = suc[x];
        if (id[y] != 1 || od[x] != 1) f |= 1u << F_C2;
    }
    for (int v = 0; v < n; v++) {
        if (od[v] == 1 && (NC >> suc[v] & 1u)) f |= 1u << F_C3;
        if (id[v] == 1 && (NC >> pre[v] & 1u)) f |= 1u << F_C3;
    }
    { int a = 0, b = 0, ya = -1, xb = -1;
      for (int v = 0; v < n; v++) if (NC >> v & 1u) { if (id[v] == 1) { a++; ya = v; } if (od[v] == 1) { b++; xb = v; } }
      if (a != 1 || b != 1 || ya != y || xb != x) f |= 1u << F_C4; }
    if (y >= 0) {
        if (hc_dfs(adj, F & ~(1u << y), NULL) != 1) f |= 1u << F_C5;
        if ((F & ~(uint32_t)adj[y] & ~(1u << y)) != (1u << pre[y])) f |= 1u << F_C5;
    } else f |= 1u << F_C5;
    /* Lemma 1 at every non-critical v: T - v has exactly one Hamiltonian cycle and it is C - v + (v- -> v+) */
    for (int v = 0; v < n; v++) if (NC >> v & 1u) {
        int c2[MAXN]; long h = hc_dfs(adj, F & ~(1u << v), c2);
        if (h != 1) { f |= 1u << F_L1; continue; }
        for (int i = 0; i < n - 1; i++) { int a = c2[i], b = c2[(i + 1) % (n - 1)];
            int expect = (suc[a] == v) ? suc[v] : suc[a];
            if (b != expect) f |= 1u << F_L1; }
    }
    /* (D) cut arcs by direct deletion */
    int k = 0;
    for (int v = 0; v < n; v++) {
        uint8_t a2[MAXN]; memcpy(a2, adj, MAXN); a2[v] &= (uint8_t)~(1u << suc[v]);
        int cut = !strong_dg(a2, n, F);
        k += cut;
        int pred = !((NC >> v & 1u) && (NC >> suc[v] & 1u));
        if (cut != pred) f |= 1u << F_D1;
    }
    if (k + r != n + 1) f |= 1u << F_D2;
    /* (D) intermediate: components of T - w at critical w */
    for (int w = 0; w < n; w++) if (!(NC >> w & 1u)) {
        uint32_t S = F & ~(1u << w), Q1 = 0, Qm = 0;
        uint8_t rev[MAXN]; memset(rev, 0, MAXN);
        for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (adj[a] >> b & 1u) rev[b] |= (uint8_t)(1u << a);
        for (int u = 0; u < n; u++) if (S >> u & 1u) {
            if (reach_fwd(adj, S, u) == S) Q1 |= 1u << u;
            if (reach_fwd(rev, S, u) == S) Qm |= 1u << u;
        }
        if (!(Q1 >> suc[w] & 1u) || !(Qm >> pre[w] & 1u)) f |= 1u << F_D3;
        if ((adj[w] & Q1) != (1u << suc[w])) f |= 1u << F_D3;
        if ((rev[w] & Qm) != (1u << pre[w])) f |= 1u << F_D3;
        if (Q1 == S || Qm == S) f |= 1u << F_D3;
    }
    *o_r = r; *o_k = k; *o_se = se; *o_sen = sen;
    return f;
}

static long fact(int n) { long f = 1; for (int i = 2; i <= n; i++) f *= i; return f; }
static long long fib(int k) { if (k < 0) return (k % 2 == 0) ? -fib(-k) : fib(-k); long long a = 0, b = 1; for (int i = 0; i < k; i++) { long long c = a + b; a = b; b = c; } return a; }

static void run(int n, FILE *cls) {
    double t0 = omp_get_wtime();
    N = n; NP = n * (n - 1) / 2;
    { int b = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { PI_[b] = i; PJ_[b] = j; b++; } }
    uint32_t TOT = 1u << NP;
    uint32_t *uh_cnt = calloc(TOT, sizeof(uint32_t));
    int32_t *inv = malloc(TOT * sizeof(int32_t)); for (uint32_t c = 0; c < TOT; c++) inv[c] = -1;
    uint8_t *lv = calloc(TOT, 1);
    uint8_t *tk = calloc(TOT, 1), *uk = calloc(TOT, 1);   /* bit k-1 set if class has hc = k (k<=5); uk: and min out >= 2 */
    long n_strong = 0, n_uh = 0, moon_fail = 0, eq_fail = 0, h9459_fail = 0, inv_incons = 0;
    long failc[NFAIL]; memset(failc, 0, sizeof failc);
    long rt5_hc = -1;

    #pragma omp parallel
    {
        long l_strong = 0, l_uh = 0, l_moon = 0, l_eq = 0, l_h = 0, l_fail[NFAIL]; memset(l_fail, 0, sizeof l_fail);
        #pragma omp for schedule(dynamic, 4096)
        for (long code = 0; code < (long)TOT; code++) {
            uint8_t adj[MAXN]; build((uint32_t)code, adj);
            uint32_t F = (1u << n) - 1;
            if (!strong_dg(adj, n, F)) continue;
            l_strong++;
            uint32_t NC = 0;
            for (int v = 0; v < n; v++) if (strong_dg(adj, n, F & ~(1u << v))) NC |= 1u << v;
            int r = __builtin_popcount(NC);
            if (n >= 4 && r < 2) l_moon++;
            if (n >= 4 && r == 2) { uint64_t cc = canon_ss(adj, n, NULL);
                #pragma omp atomic write
                lv[cc] = 1;
            }
            int cyc[MAXN]; long hc = hc_dfs(adj, F, cyc);
            int minout = n; for (int v = 0; v < n; v++) { int d = __builtin_popcount(adj[v]); if (d < minout) minout = d; }
            if (minout >= 2 && hc < 2) l_eq++;
            if (minout >= 2 && n >= 6 && hc < 3) l_h++;
            if (hc >= 1 && hc <= 5) {
                uint64_t cc = canon_ss(adj, n, NULL);
                #pragma omp atomic
                tk[cc] |= (uint8_t)(1u << (hc - 1));
                if (minout >= 2) {
                    #pragma omp atomic
                    uk[cc] |= (uint8_t)(1u << (hc - 1));
                }
            }
            if (hc != 1) continue;
            l_uh++;
            int rr, kk, se, sen; unsigned f = check_uh(adj, n, cyc, &rr, &kk, &se, &sen);
            if (hc_dp(adj, n, F) != 1) f |= 1u << F_HC2;
            long aut; uint64_t cc = canon_ss(adj, n, &aut);
            if (aut != 1) f |= 1u << F_AUT;
            for (int b = 0; b < NFAIL; b++) if (f >> b & 1u) l_fail[b]++;
            int32_t pack = rr | (kk << 5) | (se << 10) | (sen << 13);
            #pragma omp critical(uhbook)
            {
                uh_cnt[cc]++;
                if (inv[cc] == -1) inv[cc] = pack; else if (inv[cc] != pack) inv_incons++;
            }
        }
        #pragma omp critical(red)
        {
            n_strong += l_strong; n_uh += l_uh; moon_fail += l_moon; eq_fail += l_eq; h9459_fail += l_h;
            for (int b = 0; b < NFAIL; b++) failc[b] += l_fail[b];
        }
    }
    /* RT5 */
    if (n == 5) { uint8_t a[MAXN] = {0}; for (int i = 0; i < 5; i++) a[i] = (uint8_t)((1u << ((i + 1) % 5)) | (1u << ((i + 2) % 5))); rt5_hc = hc_dfs(a, 31, NULL); }
    /* A_n: i -> i+1, j -> i for j >= i+2 */
    uint8_t An[MAXN] = {0};
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) { if (j == i + 1 || i >= j + 2) An[i] |= (uint8_t)(1u << j); }
    uint64_t can_An = canon_ss(An, n, NULL); long hc_An = hc_dfs(An, (1u << n) - 1, NULL);

    long classes = 0, badsize = 0, Nr[MAXN + 2] = {0}, Tk[MAXN + 2] = {0}, oi[4][4] = {{0}}, rooted = 0;
    for (uint32_t c = 0; c < TOT; c++) if (uh_cnt[c]) {
        classes++;
        if ((long)uh_cnt[c] != fact(n)) badsize++;
        int r = inv[c] & 31, k = (inv[c] >> 5) & 31, se = (inv[c] >> 10) & 7, sen = (inv[c] >> 13) & 7;
        Nr[r]++; Tk[k]++; if (se <= 3 && sen <= 3) oi[se][sen]++; rooted += se;
        if (cls) fprintf(cls, "%d %u\n", n, c);
    }
    long nlv = 0; int lv_is_An = 1; for (uint32_t c = 0; c < TOT; c++) if (lv[c]) { nlv++; if (c != can_An) lv_is_An = 0; }
    long tcount[6] = {0}, ucount[6] = {0};
    for (uint32_t c = 0; c < TOT; c++) for (int k = 1; k <= 5; k++) { if (tk[c] >> (k - 1) & 1u) tcount[k]++; if (uk[c] >> (k - 1) & 1u) ucount[k]++; }

    printf("n=%d  labelled=%u strong=%ld UH_labelled=%ld UH_classes=%ld (Garey F(2n-6)=%lld) class_size!=n!:%ld invariant_inconsistent=%ld\n",
           n, TOT, n_strong, n_uh, classes, n >= 3 ? fib(2 * n - 6) : 0, badsize, inv_incons);
    printf("   failures:"); for (int b = 0; b < NFAIL; b++) printf(" %s=%ld", fname[b], failc[b]); printf("\n");
    printf("   N(n,r) r=2..n-1:"); for (int r = 2; r <= n - 1; r++) printf(" %ld", Nr[r]);
    printf("   | r outside [2,n-1]: %ld\n", Nr[0] + Nr[1] + Nr[n]);
    printf("   cut arcs k=2..n-1:"); for (int k = 2; k <= n - 1; k++) printf(" %ld", Tk[k]); printf("\n");
    printf("   (#sole exits,#sole entries): (1,1)=%ld (1,2)=%ld (2,1)=%ld (2,2)=%ld other=%ld ; formula F(2n-10)+1=%lld F(2n-9)-1=%lld F(2n-8)+1=%lld\n",
           oi[1][1], oi[1][2], oi[2][1], oi[2][2], classes - oi[1][1] - oi[1][2] - oi[2][1] - oi[2][2],
           fib(2 * n - 10) + 1, fib(2 * n - 9) - 1, fib(2 * n - 8) + 1);
    printf("   one sole exit=%ld (F(2n-8)=%lld) two sole exits=%ld (F(2n-7)=%lld) rooted (T,sole exit)=%ld (F(2n-5)=%lld)\n",
           oi[1][1] + oi[1][2], fib(2 * n - 8), oi[2][1] + oi[2][2], fib(2 * n - 7), rooted, fib(2 * n - 5));
    if (n >= 4) printf("   Moon fails=%ld ; strong classes with exactly 2 non-critical: %ld, all = A_n: %s ; A_n hc=%ld\n", moon_fail, nlv, lv_is_An ? "yes" : "NO", hc_An);
    printf("   strong & min out>=2 & hc<2: %ld ; (n>=6) strong & min out>=2 & hc<3: %ld", eq_fail, h9459_fail);
    if (n == 5) printf(" ; RT5 hc=%ld", rt5_hc);
    printf("\n   classes t_k (exactly k Ham cycles) k=1..5:"); for (int k = 1; k <= 5; k++) printf(" %ld", tcount[k]);
    printf(" ; u_k (and min out>=2):"); for (int k = 1; k <= 5; k++) printf(" %ld", ucount[k]);
    printf("\n   time %.2f s\n", omp_get_wtime() - t0);
    fflush(stdout);
    free(uh_cnt); free(inv); free(lv); free(tk); free(uk);
}

/* validate canon_ss against full n! brute force on ALL labelled tournaments of order n (partition + |Aut|) */
static void validate_canon(int n) {
    N = n; NP = n * (n - 1) / 2;
    { int b = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { PI_[b] = i; PJ_[b] = j; b++; } }
    uint32_t TOT = 1u << NP;
    int64_t *f2s = malloc(sizeof(int64_t) << NP), *s2f = malloc(sizeof(int64_t) << NP);
    for (uint32_t c = 0; c < TOT; c++) f2s[c] = s2f[c] = -1;
    long bad = 0, autbad = 0, ncls = 0;
    for (uint32_t code = 0; code < TOT; code++) {
        uint8_t adj[MAXN]; build(code, adj);
        long a1, a2; uint64_t cs = canon_ss(adj, n, &a1), cf = canon_full(adj, n, &a2);
        if (a1 != a2) autbad++;
        if (f2s[cf] == -1) { f2s[cf] = (int64_t)cs; ncls++; } else if (f2s[cf] != (int64_t)cs) bad++;
        if (s2f[cs] == -1) s2f[cs] = (int64_t)cf; else if (s2f[cs] != (int64_t)cf) bad++;
    }
    printf("canon validation n=%d: classes=%ld partition_mismatch=%ld aut_mismatch=%ld\n", n, ncls, bad, autbad);
    free(f2s); free(s2f);
}

int main(int argc, char **argv) {
    int nmax = argc > 1 ? atoi(argv[1]) : 7;
    double t0 = omp_get_wtime();
    for (int n = 3; n <= 6 && n <= nmax; n++) validate_canon(n);
    FILE *cls = fopen("audit_labelled_uh_classes.txt", "w");
    for (int n = 3; n <= nmax; n++) run(n, cls);
    if (cls) fclose(cls);
    printf("total time %.2f s\n", omp_get_wtime() - t0);
    return 0;
}
