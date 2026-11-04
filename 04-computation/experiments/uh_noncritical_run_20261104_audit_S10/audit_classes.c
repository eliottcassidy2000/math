/* audit_classes.c -- independent audit (S10 auditor) of THM-4868 at the level of isomorphism classes.
 *
 * PART I (no structural assumption).  All tournament classes of orders 1..8 by plain one-vertex extension
 *   (every T of order n is an extension of a class representative of T - v), canonical form = score-sorted
 *   branch-and-bound minimum (same definition as audit_labelled.c).  The UH classes of order 8 are read off the
 *   6880 classes; the UH classes of order 9 come from extending all 6880 classes of order 8 by all 2^8 out-sets
 *   and counting Hamiltonian cycles exactly (subset DP).  Same sweep: Las Vergnas (strong, exactly two
 *   non-critical vertices => A_n) at orders 8, 9; Moon at 8; HYP-9459's t_k, u_k (k <= 5) at orders 8, 9.
 *
 * PART II (Shelor-style insertion, NOT the author's interval extension).  From every UH class T' of order n-1 and
 *   every parent p: add y with p -> y and y -> every other vertex; keep iff hc(T) = 1, where hc is an exact subset
 *   DP on T itself (n <= DIRECT); for n > DIRECT, hc(T) = e_p(T') (Hamiltonian paths of T' ending at p; elementary)
 *   is used only as a filter and every accepted T is re-counted by the exact DP on T.  Canonical form: maximum over
 *   the n rotations of the unique Hamiltonian cycle of a 120-bit upper-triangle code (valid because hc = 1 is
 *   verified exactly, so every isomorphism maps the cycle to the cycle).  Completeness is NOT assumed: it is
 *   cross-checked against PART I (orders <= 9), audit_labelled.c (orders <= 7) and Garey's F(2n-6).
 *   Per class: (A1) roundness, (A2) propagation, (B) degree-1 counts and adjacency, (C) one run, ends, Lemma 2,
 *   unique non-critical sole entry/exit, Shelor parent T - y found among the order-(n-1) classes, (D) cut arcs by
 *   direct deletion and the component structure at critical w, Lemma 1 (n <= L1MAX), (E) rotations, (F) the
 *   two-sole-exit bijection checked directly (hits every rooted pair (T', sole exit) of order n-1 exactly once),
 *   Las Vergnas class = A_n, and the tables N(n,r), T(n,k), (#sole exits, #sole entries); P7 recurrence.
 *
 * Build: gcc -O2 -march=native -fopenmp -o audit_classes audit_classes.c ; run: ./audit_classes NMAX [DIRECT=13] [L1MAX=12]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXN 17
typedef unsigned __int128 u128;
static int popc(uint32_t x) { return __builtin_popcount(x); }
static int ctz(uint32_t x) { return __builtin_ctz(x); }
static long long fib(int k) { if (k < 0) return (k % 2 == 0) ? -fib(-k) : fib(-k); long long a = 0, b = 1; for (int i = 0; i < k; i++) { long long c = a + b; a = b; b = c; } return a; }

/* ---------------- digraph helpers ---------------- */
static int strong_dg(const uint32_t *out, int n, uint32_t S) {
    if (popc(S) <= 1) return 1;
    uint32_t in[MAXN] = {0};
    for (int x = 0; x < n; x++) if (S >> x & 1u) { uint32_t o = out[x] & S; while (o) { int y = ctz(o); o &= o - 1; in[y] |= 1u << x; } }
    int s = ctz(S);
    uint32_t seen = 1u << s, fr = seen;
    while (fr) { uint32_t nf = 0; for (uint32_t f = fr; f; f &= f - 1) nf |= out[ctz(f)]; nf &= S & ~seen; seen |= nf; fr = nf; }
    if (seen != S) return 0;
    seen = 1u << s; fr = seen;
    while (fr) { uint32_t nf = 0; for (uint32_t f = fr; f; f &= f - 1) nf |= in[ctz(f)]; nf &= S & ~seen; seen |= nf; fr = nf; }
    return seen == S;
}
static uint32_t reach(const uint32_t *out, uint32_t S, int s) {
    uint32_t seen = 1u << s, fr = seen;
    while (fr) { uint32_t nf = 0; for (uint32_t f = fr; f; f &= f - 1) nf |= out[ctz(f)]; nf &= S & ~seen; seen |= nf; fr = nf; }
    return seen;
}

/* exact Hamiltonian-cycle count of the whole tournament (paths from vertex 0); dp has 2^(n-1)*n entries.
 * If the count is 1 and cyc != NULL, the unique cycle is written as cyc[0] = 0, cyc[1], ..., cyc[n-1]. */
static uint64_t hc_exact(const uint32_t *adj, int n, uint64_t *dp, int *cyc) {
    if (n < 3) return 0;
    uint32_t M = 1u << n;
    #define IDX(m, v) ((size_t)((m) >> 1) * (size_t)n + (size_t)(v))
    memset(dp, 0, sizeof(uint64_t) * (size_t)(M >> 1) * n);
    dp[IDX(1u, 0)] = 1;
    for (uint32_t m = 1; m < M; m += 2) {
        for (int v = 0; v < n; v++) { uint64_t c = dp[IDX(m, v)]; if (!c) continue;
            uint32_t nx = adj[v] & ~m & (M - 1);
            while (nx) { int w = ctz(nx); nx &= nx - 1; dp[IDX(m | (1u << w), w)] += c; } }
    }
    uint64_t tot = 0;
    for (int v = 1; v < n; v++) if (adj[v] & 1u) tot += dp[IDX(M - 1, v)];
    if (cyc && tot == 1) {
        uint32_t mask = M - 1; int v = -1;
        for (int u = 1; u < n; u++) if ((adj[u] & 1u) && dp[IDX(mask, u)]) v = u;
        for (int pos = n - 1; pos >= 1; pos--) {
            cyc[pos] = v; uint32_t pm = mask & ~(1u << v); int pu = -1, cntu = 0;
            for (int u = 0; u < n; u++) if ((pm >> u & 1u) && (adj[u] >> v & 1u) && dp[IDX(pm, u)]) { pu = u; cntu++; }
            if (cntu != 1) { fprintf(stderr, "backtrack ambiguity\n"); exit(2); }
            mask = pm; v = pu;
        }
        cyc[0] = v; if (v != 0) { fprintf(stderr, "backtrack did not end at 0\n"); exit(2); }
    }
    return tot;
    #undef IDX
}
/* Hamiltonian cycles of T[S]: compress and count */
static uint64_t hc_sub(const uint32_t *adj, int n, uint32_t S, uint64_t *dp) {
    int map[MAXN], k = 0; for (int v = 0; v < n; v++) map[v] = (S >> v & 1u) ? k++ : -1;
    uint32_t a2[MAXN] = {0};
    for (int v = 0; v < n; v++) if (map[v] >= 0) { uint32_t o = adj[v] & S; while (o) { int w = ctz(o); o &= o - 1; a2[map[v]] |= 1u << map[w]; } }
    return hc_exact(a2, k, dp, NULL);
}
/* Hamiltonian-path end counts of T (any start), saturated at 2; dp has 2^m * m bytes */
static void endpoints(const uint32_t *adj, int m, uint8_t *dp, uint8_t *e) {
    size_t M = (size_t)1 << m;
    memset(dp, 0, M * m);
    for (int v = 0; v < m; v++) dp[((size_t)1 << v) * m + v] = 1;
    for (size_t mask = 1; mask < M; mask++)
        for (int v = 0; v < m; v++) { uint8_t c = dp[mask * m + v]; if (!c) continue;
            uint32_t nx = adj[v] & ~(uint32_t)mask & (uint32_t)(M - 1);
            while (nx) { int w = ctz(nx); nx &= nx - 1; size_t id = (mask | ((size_t)1 << w)) * m + w; int s = dp[id] + c; dp[id] = (uint8_t)(s > 2 ? 2 : s); } }
    for (int v = 0; v < m; v++) e[v] = dp[(M - 1) * m + v];
}

/* ---------------- score-sorted canonical form (identical definition to audit_labelled.c) ---------------- */
typedef struct { const uint32_t *adj; int n, NPn; int deg[MAXN], slot[MAXN], perm[MAXN]; uint32_t used; uint64_t best; int have; long cnt; } canon_t;
static void canon_rec(canon_t *c, int i, uint64_t pc) {
    if (i == c->n) { if (!c->have || pc < c->best) { c->best = pc; c->have = 1; c->cnt = 1; } else if (pc == c->best) c->cnt++; return; }
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
static uint64_t canon_ss(const uint32_t *adj, int n, long *aut) {
    canon_t c; memset(&c, 0, sizeof c); c.adj = adj; c.n = n; c.NPn = n * (n - 1) / 2;
    for (int v = 0; v < n; v++) c.deg[v] = popc(adj[v]);
    int tmp[MAXN]; memcpy(tmp, c.deg, sizeof tmp);
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (tmp[b] < tmp[a]) { int t = tmp[a]; tmp[a] = tmp[b]; tmp[b] = t; }
    memcpy(c.slot, tmp, sizeof tmp);
    canon_rec(&c, 0, 0);
    if (aut) *aut = c.cnt;
    return c.best;
}
static void decode_ss(uint64_t code, int n, uint32_t *adj) {
    int NPn = n * (n - 1) / 2; memset(adj, 0, sizeof(uint32_t) * MAXN);
    for (int i = 1; i < n; i++) for (int j = 0; j < i; j++) {
        int p = i * (i - 1) / 2 + j;
        if (code >> (NPn - 1 - p) & 1u) adj[j] |= 1u << i; else adj[i] |= 1u << j;
    }
}
static int cmp_u64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return x < y ? -1 : x > y; }

/* ---------------- rotation canonical form (UH tournaments in cycle coordinates: i -> i+1 mod n) ---------------- */
static u128 code_rot(const uint32_t *adj, int n, int s) {
    u128 c = 0;
    for (int a = 0; a < n; a++) { uint32_t row = adj[(a + s) % n];
        for (int b = a + 1; b < n; b++) c = (c << 1) | (u128)((row >> ((b + s) % n)) & 1u); }
    return c;
}
static u128 canon_rot(const uint32_t *adj, int n, int *nfix, int *arg) {
    u128 c0 = code_rot(adj, n, 0), best = c0; int fix = 1, a = 0;
    for (int s = 1; s < n; s++) { u128 c = code_rot(adj, n, s); if (c == c0) fix++; if (c > best) { best = c; a = s; } }
    if (nfix) *nfix = fix; if (arg) *arg = a; return best;
}
static void rotate_adj(const uint32_t *adj, int n, int s, uint32_t *out) {
    for (int a = 0; a < n; a++) { uint32_t row = adj[(a + s) % n], r = 0;
        for (int b = 0; b < n; b++) if (row >> ((b + s) % n) & 1u) r |= 1u << b; out[a] = r; }
}
/* relabel T so that the Hamiltonian cycle cyc[0] -> cyc[1] -> ... becomes 0 -> 1 -> ... */
static void to_cycle_coords(const uint32_t *adj, int n, const int *cyc, uint32_t *out) {
    int pos[MAXN]; for (int i = 0; i < n; i++) pos[cyc[i]] = i;
    for (int i = 0; i < n; i++) { uint32_t row = adj[cyc[i]], r = 0; while (row) { int w = ctz(row); row &= row - 1; r |= 1u << pos[w]; } out[i] = r; }
}
/* delete vertex v of a UH tournament in cycle coordinates; result in cycle coordinates of C - v + (v- -> v+),
 * starting at v+.  Returns 0 if the closing arc v- -> v+ is missing. */
static int delete_in_cycle(const uint32_t *adj, int n, int v, uint32_t *out, int *newlabel_of) {
    int m = n - 1, old[MAXN];
    for (int i = 0; i < m; i++) old[i] = (v + 1 + i) % n;
    for (int u = 0; u < n; u++) newlabel_of[u] = -1;
    for (int i = 0; i < m; i++) newlabel_of[old[i]] = i;
    for (int i = 0; i < m; i++) { uint32_t row = adj[old[i]], r = 0; while (row) { int w = ctz(row); row &= row - 1; if (w != v) r |= 1u << newlabel_of[w]; } out[i] = r; }
    return (adj[(v + n - 1) % n] >> ((v + 1) % n)) & 1u;
}

/* ---------------- levels ---------------- */
typedef struct { u128 code; uint32_t adj[MAXN]; } cls_t;
typedef struct { int n; long cnt; cls_t *v; } level_t;
static int cmp_cls(const void *a, const void *b) { u128 x = ((const cls_t *)a)->code, y = ((const cls_t *)b)->code; return x < y ? -1 : x > y; }
static long find_cls(const level_t *L, u128 code) {
    long lo = 0, hi = L->cnt - 1;
    while (lo <= hi) { long mid = (lo + hi) / 2; if (L->v[mid].code == code) return mid; if (L->v[mid].code < code) lo = mid + 1; else hi = mid - 1; }
    return -1;
}

enum { F_A1, F_A2, F_B1, F_B2, F_C1, F_C2, F_C3, F_C4, F_C5, F_D1, F_D2, F_D3, F_L1, F_E, F_HC, NFAIL };
static const char *fname[NFAIL] = { "A1", "A2", "B1", "B2", "C1", "C2", "C3", "C4", "C5_shelor_parent", "D1", "D2", "D3", "L1", "E", "hc" };

typedef struct { long fail[NFAIL]; long Nr[MAXN + 2], Tk[MAXN + 2], oi[4][4], rooted; } stats_t;

/* checks on one UH class in cycle coordinates; returns r; fills y, x, se, sen */
static int check_class(const uint32_t *adj, int n, const level_t *prev, int L1MAX, int HCMAX_TY, uint64_t *dp, stats_t *st,
                       int *o_y, int *o_x) {
    unsigned f = 0; uint32_t F = (n >= 32) ? 0xffffffffu : ((1u << n) - 1);
    int od[MAXN], id[MAXN];
    for (int v = 0; v < n; v++) { od[v] = popc(adj[v]); id[v] = n - 1 - od[v]; if (!(adj[v] >> ((v + 1) % n) & 1u)) f |= 1u << F_HC; }
    uint32_t NC = 0;
    for (int v = 0; v < n; v++) if (strong_dg(adj, n, F & ~(1u << v))) NC |= 1u << v;
    int r = popc(NC);
#define SUC(v) (((v) + 1) % n)
#define PRE(v) (((v) + n - 1) % n)
    for (int v = 0; v < n; v++) if (NC >> v & 1u) {
        uint32_t want = 0; for (int s = 1; s <= od[v]; s++) want |= 1u << ((v + s) % n);
        if (adj[v] != want) f |= 1u << F_A1;
        if (od[v] >= 2 && !(NC >> SUC(v) & 1u)) f |= 1u << F_A2;
        if (id[v] >= 2 && !(NC >> PRE(v) & 1u)) f |= 1u << F_A2;
    }
    int se = 0, sen = 0, sev[MAXN], senv[MAXN];
    for (int v = 0; v < n; v++) { if (od[v] == 1) sev[se++] = v; if (id[v] == 1) senv[sen++] = v; }
    if (se < 1 || sen < 1) f |= 1u << F_B1;
    if (se > 2 || sen > 2) f |= 1u << F_B2;
    if (se == 2 && !(SUC(sev[0]) == sev[1] || SUC(sev[1]) == sev[0])) f |= 1u << F_B2;
    if (sen == 2 && !(SUC(senv[0]) == senv[1] || SUC(senv[1]) == senv[0])) f |= 1u << F_B2;
    int starts = 0, y = -1, x = -1;
    for (int v = 0; v < n; v++) if ((NC >> v & 1u) && !(NC >> PRE(v) & 1u)) { starts++; y = v; }
    if (starts != 1 || r < 2 || r > n - 1) f |= 1u << F_C1;
    if (starts == 1) { x = y; while (NC >> SUC(x) & 1u) x = SUC(x); if (id[y] != 1 || od[x] != 1) f |= 1u << F_C2; }
    for (int v = 0; v < n; v++) { if (od[v] == 1 && (NC >> SUC(v) & 1u)) f |= 1u << F_C3; if (id[v] == 1 && (NC >> PRE(v) & 1u)) f |= 1u << F_C3; }
    { int a = 0, b = 0, ya = -1, xb = -1;
      for (int v = 0; v < n; v++) if (NC >> v & 1u) { if (id[v] == 1) { a++; ya = v; } if (od[v] == 1) { b++; xb = v; } }
      if (a != 1 || b != 1 || ya != y || xb != x) f |= 1u << F_C4; }
    /* Shelor parent: y's only in-neighbour is y-, and T - y (along C - y + (y- -> y+)) is a known UH class */
    if (y >= 0) {
        if ((F & ~adj[y] & ~(1u << y)) != (1u << PRE(y))) f |= 1u << F_C5;
        uint32_t sub[MAXN]; int nl[MAXN];
        if (!delete_in_cycle(adj, n, y, sub, nl)) f |= 1u << F_C5;
        else if (prev && prev->n == n - 1) { u128 c = canon_rot(sub, n - 1, NULL, NULL); if (find_cls(prev, c) < 0) f |= 1u << F_C5; }
        if (n <= HCMAX_TY && hc_sub(adj, n, F & ~(1u << y), dp) != 1) f |= 1u << F_C5;
    } else f |= 1u << F_C5;
    /* Lemma 1 at every non-critical v (n <= L1MAX): hc(T - v) = 1 */
    if (n <= L1MAX) for (int v = 0; v < n; v++) if (NC >> v & 1u) if (hc_sub(adj, n, F & ~(1u << v), dp) != 1) f |= 1u << F_L1;
    /* (D) cut arcs */
    int k = 0;
    for (int v = 0; v < n; v++) {
        uint32_t a2[MAXN]; memcpy(a2, adj, sizeof(uint32_t) * n); a2[v] &= ~(1u << SUC(v));
        int cut = !strong_dg(a2, n, F); k += cut;
        int pred = !((NC >> v & 1u) && (NC >> SUC(v) & 1u));
        if (cut != pred) f |= 1u << F_D1;
    }
    if (k + r != n + 1) f |= 1u << F_D2;
    {   uint32_t rev[MAXN] = {0};
        for (int a = 0; a < n; a++) { uint32_t o = adj[a]; while (o) { int b = ctz(o); o &= o - 1; rev[b] |= 1u << a; } }
        for (int w = 0; w < n; w++) if (!(NC >> w & 1u)) {
            uint32_t S = F & ~(1u << w), Q1 = 0, Qm = 0;
            for (int u = 0; u < n; u++) if (S >> u & 1u) { if (reach(adj, S, u) == S) Q1 |= 1u << u; if (reach(rev, S, u) == S) Qm |= 1u << u; }
            if (!(Q1 >> SUC(w) & 1u) || !(Qm >> PRE(w) & 1u)) f |= 1u << F_D3;
            if ((adj[w] & Q1) != (1u << SUC(w)) || (rev[w] & Qm) != (1u << PRE(w))) f |= 1u << F_D3;
            if (Q1 == S || Qm == S) f |= 1u << F_D3;
        } }
    { int fix; canon_rot(adj, n, &fix, NULL); if (fix != 1) f |= 1u << F_E; }
    for (int b = 0; b < NFAIL; b++) if (f >> b & 1u) st->fail[b]++;
    st->Nr[r]++; st->Tk[k]++; if (se <= 3 && sen <= 3) st->oi[se][sen]++; st->rooted += se;
    *o_y = y; *o_x = x;
    return r;
#undef SUC
#undef PRE
}

/* ===================================================================================================== */
static uint64_t *uh_ss[MAXN + 1]; static long uh_ss_cnt[MAXN + 1];   /* PART I UH classes (score-sorted codes) */

static void part1(void) {
    double t0 = omp_get_wtime();
    uint64_t *cur = malloc(sizeof(uint64_t)); long ncur = 1; cur[0] = 0;   /* order 1 */
    long expect[9] = {0, 1, 1, 2, 4, 12, 56, 456, 6880};
    for (int n = 2; n <= 8; n++) {
        int m = n - 1; long cap = ncur << m; uint64_t *nx = malloc(sizeof(uint64_t) * cap);
        #pragma omp parallel for schedule(dynamic, 8)
        for (long k = 0; k < ncur; k++) {
            uint32_t base[MAXN]; decode_ss(cur[k], m, base);
            for (uint32_t O = 0; O < (1u << m); O++) {
                uint32_t adj[MAXN]; memcpy(adj, base, sizeof adj);
                for (int u = 0; u < m; u++) if (!(O >> u & 1u)) adj[u] |= 1u << m;
                adj[m] = O;
                nx[(k << m) | O] = canon_ss(adj, n, NULL);
            }
        }
        qsort(nx, cap, sizeof(uint64_t), cmp_u64);
        long u = 0; for (long i = 0; i < cap; i++) if (i == 0 || nx[i] != nx[i - 1]) nx[u++] = nx[i];
        free(cur); cur = nx; ncur = u;
        printf("PART I: all tournament classes of order %d: %ld (A000568: %ld) %s\n", n, ncur, expect[n], ncur == expect[n] ? "OK" : "MISMATCH");
    }
    /* order 8: UH classes, Las Vergnas, Moon, t_k, u_k */
    {
        int n = 8; uint32_t F = 255; long uh = 0, lvc = 0, moon = 0, tk[6] = {0}, uk[6] = {0}, eqf = 0, strongc = 0;
        uint64_t *uhl = malloc(sizeof(uint64_t) * ncur); uint64_t lvcode = 0; int lv_ok = 1;
        uint32_t An[MAXN] = {0}; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && (j == i + 1 || i >= j + 2)) An[i] |= 1u << j;
        uint64_t canAn = canon_ss(An, n, NULL);
        uint64_t *dp = malloc(sizeof(uint64_t) * (1u << (n - 1)) * n);
        for (long k = 0; k < ncur; k++) {
            uint32_t adj[MAXN]; decode_ss(cur[k], n, adj);
            if (!strong_dg(adj, n, F)) continue;
            strongc++;
            int r = 0; for (int v = 0; v < n; v++) r += strong_dg(adj, n, F & ~(1u << v));
            if (r < 2) moon++;
            if (r == 2) { lvc++; lvcode = cur[k]; if (cur[k] != canAn) lv_ok = 0; }
            uint64_t h = hc_exact(adj, n, dp, NULL);
            int mo = n; for (int v = 0; v < n; v++) if (popc(adj[v]) < mo) mo = popc(adj[v]);
            if (mo >= 2 && h < 3) eqf++;
            if (h >= 1 && h <= 5) { tk[h]++; if (mo >= 2) uk[h]++; }
            if (h == 1) uhl[uh++] = cur[k];
        }
        (void)lvcode;
        uh_ss[8] = uhl; uh_ss_cnt[8] = uh;
        printf("PART I order 8: strong classes %ld; UH classes %ld (Garey 55); Moon fails %ld; exactly-2-noncritical classes %ld, = A_8: %s; strong & min out>=2 & hc<3: %ld\n",
               strongc, uh, moon, lvc, lv_ok ? "yes" : "NO", eqf);
        printf("   t_k(8) k=1..5: %ld %ld %ld %ld %ld (HYP-9459: 55 60 124 91 138) ; u_k(8): %ld %ld %ld %ld %ld (HYP-9459: 0 0 13 10 9)\n",
               tk[1], tk[2], tk[3], tk[4], tk[5], uk[1], uk[2], uk[3], uk[4], uk[5]);
        free(dp);
    }
    /* order 9: extend all 6880 classes of order 8 by every out-set */
    {
        int n = 9, m = 8; uint32_t F = 511;
        uint64_t *uhl = NULL; long nuh = 0, capuh = 0, nlv = 0, caplv = 0; uint64_t *lvl = NULL;
        uint64_t *tkl = NULL; long ntk = 0, captk = 0;
        long eqf = 0, ext_total = 0;
        #pragma omp parallel
        {
            uint64_t *dp = malloc(sizeof(uint64_t) * (1u << (n - 1)) * n);
            long l_eq = 0, l_tot = 0;
            #pragma omp for schedule(dynamic, 4)
            for (long k = 0; k < ncur; k++) {
                uint32_t base[MAXN]; decode_ss(cur[k], m, base);
                for (uint32_t O = 0; O < (1u << m); O++) {
                    uint32_t adj[MAXN]; memcpy(adj, base, sizeof adj);
                    for (int u = 0; u < m; u++) if (!(O >> u & 1u)) adj[u] |= 1u << m;
                    adj[m] = O; l_tot++;
                    if (!strong_dg(adj, n, F)) continue;
                    int r = 0; for (int v = 0; v < n; v++) r += strong_dg(adj, n, F & ~(1u << v));
                    if (r == 2) { uint64_t c = canon_ss(adj, n, NULL);
                        #pragma omp critical(lvadd)
                        { if (nlv == caplv) { caplv = caplv ? 2 * caplv : 1024; lvl = realloc(lvl, sizeof(uint64_t) * caplv); } lvl[nlv++] = c; } }
                    uint64_t h = hc_exact(adj, n, dp, NULL);
                    int mo = n; for (int v = 0; v < n; v++) if (popc(adj[v]) < mo) mo = popc(adj[v]);
                    if (mo >= 2 && h < 3) l_eq++;
                    if (h >= 1 && h <= 5) { uint64_t c = canon_ss(adj, n, NULL);
                        uint64_t tagged = (c << 4) | (uint64_t)(h << 1) | (uint64_t)(mo >= 2);
                        #pragma omp critical(tkadd)
                        { if (ntk == captk) { captk = captk ? 2 * captk : 4096; tkl = realloc(tkl, sizeof(uint64_t) * captk); } tkl[ntk++] = tagged; } }
                    if (h == 1) { uint64_t c = canon_ss(adj, n, NULL);
                        #pragma omp critical(uhadd)
                        { if (nuh == capuh) { capuh = capuh ? 2 * capuh : 1024; uhl = realloc(uhl, sizeof(uint64_t) * capuh); } uhl[nuh++] = c; } }
                }
            }
            #pragma omp critical(red9)
            { eqf += l_eq; ext_total += l_tot; }
            free(dp);
        }
        qsort(uhl, nuh, sizeof(uint64_t), cmp_u64); long u = 0; for (long i = 0; i < nuh; i++) if (i == 0 || uhl[i] != uhl[i - 1]) uhl[u++] = uhl[i];
        long uh_inst = nuh; nuh = u; uh_ss[9] = uhl; uh_ss_cnt[9] = nuh;
        qsort(lvl, nlv, sizeof(uint64_t), cmp_u64); long lu = 0; for (long i = 0; i < nlv; i++) if (i == 0 || lvl[i] != lvl[i - 1]) lvl[lu++] = lvl[i];
        uint32_t An[MAXN] = {0}; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && (j == i + 1 || i >= j + 2)) An[i] |= 1u << j;
        uint64_t canAn = canon_ss(An, n, NULL);
        qsort(tkl, ntk, sizeof(uint64_t), cmp_u64); long tk[6] = {0}, uk[6] = {0};
        for (long i = 0; i < ntk; i++) if (i == 0 || tkl[i] != tkl[i - 1]) { int h = (int)((tkl[i] >> 1) & 7); tk[h]++; if (tkl[i] & 1u) uk[h]++; }
        printf("PART I order 9: %ld extensions; UH instances %ld -> UH classes %ld (Garey 144); exactly-2-noncritical classes %ld, = A_9: %s; strong & min out>=2 & hc<3: %ld\n",
               ext_total, uh_inst, nuh, lu, (lu == 1 && lvl[0] == canAn) ? "yes" : "NO", eqf);
        printf("   t_k(9) k=1..5: %ld %ld %ld %ld %ld (HYP-9459: 144 191 390 315 466) ; u_k(9): %ld %ld %ld %ld %ld (HYP-9459: 0 0 34 26 23)\n",
               tk[1], tk[2], tk[3], tk[4], tk[5], uk[1], uk[2], uk[3], uk[4], uk[5]);
        free(lvl); free(tkl);
    }
    free(cur);
    printf("PART I time %.1f s\n", omp_get_wtime() - t0);
    fflush(stdout);
}

/* ===================================================================================================== */
/* published tables for comparison */
static const long HYP9458[14][12] = {   /* T(n,k), k = 2..n-1, rows n = 5..13 */
    [5] = {1, 1, 1}, [6] = {2, 3, 2, 1}, [7] = {4, 7, 6, 3, 1}, [8] = {8, 16, 16, 10, 4, 1}, [9] = {16, 36, 41, 30, 15, 5, 1},
    [10] = {32, 80, 102, 85, 50, 21, 6, 1}, [11] = {64, 176, 248, 231, 155, 77, 28, 7, 1},
    [12] = {128, 384, 592, 608, 456, 259, 112, 36, 8, 1}, [13] = {256, 832, 1392, 1560, 1289, 819, 406, 156, 45, 9, 1} };
static const long HYP9458_16[14] = {2048, 8192, 16896, 23552, 24600, 20256, 13504, 7380, 3312, 1210, 352, 78, 12, 1};
static const long P7pred[19][17] = {     /* HYP-9460 P7 predicted N(n,r), r = 2..n-1 */
    [14] = {1, 10, 55, 210, 606, 1372, 2471, 3530, 3920, 3232, 1792, 512},
    [15] = {1, 11, 66, 275, 870, 2178, 4396, 7183, 9420, 9680, 7424, 3840, 1024},
    [16] = {1, 12, 78, 352, 1210, 3312, 7380, 13504, 20256, 24600, 23552, 16896, 8192, 2048},
    [17] = {1, 13, 91, 442, 1639, 4862, 11826, 23868, 40081, 55692, 63072, 56576, 38144, 17408, 4096},
    [18] = {1, 14, 105, 546, 2171, 6930, 18238, 40140, 74313, 115598, 149856, 159168, 134400, 85504, 36864, 8192} };
static const long NOTE_N[14][12] = {    /* the note's section-4 table N(n,r), r = 2..n-1 */
    [5] = {1, 1, 1}, [6] = {1, 2, 3, 2}, [7] = {1, 3, 6, 7, 4}, [8] = {1, 4, 10, 16, 16, 8}, [9] = {1, 5, 15, 30, 41, 36, 16},
    [10] = {1, 6, 21, 50, 85, 102, 80, 32}, [11] = {1, 7, 28, 77, 155, 231, 248, 176, 64},
    [12] = {1, 8, 36, 112, 259, 456, 608, 592, 384, 128}, [13] = {1, 9, 45, 156, 406, 819, 1289, 1560, 1392, 832, 256} };

static long NN[MAXN + 2][MAXN + 2];   /* my N(n,r) */

/* negative controls: the per-class checks must fire on inputs that violate the theorem's hypotheses */
static void selftest(void) {
    uint64_t *dp = malloc(sizeof(uint64_t) * (1u << 15) * 16);
    struct { const char *name; int n; uint32_t adj[MAXN]; } cases[3];
    /* RT5 = Circ(Z5,{1,2}), two Hamiltonian cycles, in cycle coordinates */
    cases[0].name = "RT5 (hc=2)"; cases[0].n = 5; for (int i = 0; i < 5; i++) cases[0].adj[i] = (1u << ((i + 1) % 5)) | (1u << ((i + 2) % 5));
    /* A_7 (UH) with the chord 6 -> 0 ... reversed into 0 -> 6? that is the closing arc; flip chord 5 -> 2 to 2 -> 5 */
    { int n = 7; cases[1].name = "A_7 with chord 5->2 reversed"; cases[1].n = n; memset(cases[1].adj, 0, sizeof(cases[1].adj));
      for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && (j == i + 1 || i >= j + 2)) cases[1].adj[i] |= 1u << j;
      cases[1].adj[5] &= ~(1u << 2); cases[1].adj[2] |= 1u << 5; }
    /* A_7 relabelled so that i -> i+1 is NOT its Hamiltonian cycle (swap labels 2 and 4) */
    { int n = 7; uint32_t a[MAXN] = {0}; cases[2].name = "A_7 with labels 2,4 swapped"; cases[2].n = n; memset(cases[2].adj, 0, sizeof(cases[2].adj));
      for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && (j == i + 1 || i >= j + 2)) a[i] |= 1u << j;
      int pi[7] = {0, 1, 4, 3, 2, 5, 6};
      for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (a[i] >> j & 1u) cases[2].adj[pi[i]] |= 1u << pi[j]; }
    for (int c = 0; c < 3; c++) {
        stats_t st; memset(&st, 0, sizeof st); int y, x;
        uint64_t h = hc_exact(cases[c].adj, cases[c].n, dp, NULL);
        check_class(cases[c].adj, cases[c].n, NULL, 12, 15, dp, &st, &y, &x);
        printf("selftest %-32s hc=%llu  fired:", cases[c].name, (unsigned long long)h);
        int any = 0; for (int b = 0; b < NFAIL; b++) if (st.fail[b]) { printf(" %s", fname[b]); any = 1; }
        printf("%s\n", any ? "" : " NONE (check is vacuous!)");
    }
    free(dp);
}

int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "selftest")) { selftest(); return 0; }
    int NMAX = argc > 1 ? atoi(argv[1]) : 13;
    int DIRECT = argc > 2 ? atoi(argv[2]) : 13;
    int L1MAX = argc > 3 ? atoi(argv[3]) : 12;
    int HCMAX_TY = argc > 4 ? atoi(argv[4]) : 15;
    if (NMAX > 16) { fprintf(stderr, "NMAX <= 16 (128-bit codes)\n"); return 1; }
    double T0 = omp_get_wtime();
    printf("# audit_classes NMAX=%d DIRECT=%d L1MAX=%d HCMAX_TY=%d threads=%d\n", NMAX, DIRECT, L1MAX, HCMAX_TY, omp_get_max_threads());
    part1();

    /* read audit_labelled's UH classes (score-sorted codes), orders <= 7 */
    uint64_t *lab[MAXN + 1] = {0}; long labc[MAXN + 1] = {0};
    { FILE *fp = fopen("audit_labelled_uh_classes.txt", "r");
      if (fp) { int nn; unsigned long long c; while (fscanf(fp, "%d %llu", &nn, &c) == 2) { lab[nn] = realloc(lab[nn], sizeof(uint64_t) * (labc[nn] + 1)); lab[nn][labc[nn]++] = c; } fclose(fp);
                for (int i = 0; i <= MAXN; i++) if (labc[i]) qsort(lab[i], labc[i], sizeof(uint64_t), cmp_u64); }
      else printf("(audit_labelled_uh_classes.txt not found: skipping labelled cross-check)\n"); }

    level_t prev; prev.n = 3; prev.cnt = 1; prev.v = calloc(1, sizeof(cls_t));
    { uint32_t c3[MAXN] = {2, 4, 1}; memcpy(prev.v[0].adj, c3, sizeof c3); prev.v[0].code = canon_rot(c3, 3, NULL, NULL); }
    NN[3][0] = 1;
    for (int n = 4; n <= NMAX; n++) {
        double t0 = omp_get_wtime();
        int m = n - 1;
        cls_t *cand = NULL; long ncand = 0, capc = 0;
        long lemma_mismatch = 0, hc_reject_after_filter = 0, considered = 0;
        /* THM-4868 (G), (H) checks on every parent U of order m = n - 1 (m >= 4) */
        long G_bad = 0, adm2 = 0, adm3 = 0, admother = 0, H_bad = 0, Hfact_bad = 0, serule_bad = 0, children = 0;
        #pragma omp parallel
        {
            uint64_t *dp = malloc(sizeof(uint64_t) * ((size_t)1 << (n - 1)) * n);
            uint8_t *dpe = malloc(((size_t)1 << m) * m);
            long l_lm = 0, l_rej = 0, l_cons = 0, l_G = 0, l_a2 = 0, l_a3 = 0, l_ao = 0, l_H = 0, l_Hf = 0, l_se = 0, l_ch = 0;
            #pragma omp for schedule(dynamic, 4)
            for (long k = 0; k < prev.cnt; k++) {
                const uint32_t *pa = prev.v[k].adj; uint8_t e[MAXN];
                endpoints(pa, m, dpe, e);
                /* parent structure: NC(U), run start yU, run end xU, sole entries */
                uint32_t FU = (1u << m) - 1, NCU = 0; int yU = -1, xU = -1, rU = 0;
                for (int v = 0; v < m; v++) if (strong_dg(pa, m, FU & ~(1u << v))) NCU |= 1u << v;
                rU = popc(NCU);
                for (int v = 0; v < m; v++) if ((NCU >> v & 1u) && !(NCU >> ((v + m - 1) % m) & 1u)) yU = v;
                int okU = (m >= 4 && yU >= 0 && rU < m);
                if (m >= 4 && !okU) l_G++;              /* no proper run: would already contradict (C) */
                if (okU) { xU = yU; while (NCU >> ((xU + 1) % m) & 1u) xU = (xU + 1) % m; } else yU = 0;
                int yUm = (yU + m - 1) % m, yUmm = (yU + m - 2) % m;
                int yUm_sole_entry = (m - 1 - popc(pa[yUm])) == 1;
                if (okU) {
                    uint32_t adm = 0, pred = (1u << yU) | (1u << yUm) | (yUm_sole_entry ? (1u << yUmm) : 0u);
                    for (int p = 0; p < m; p++) if (e[p] == 1) adm |= 1u << p;
                    if (adm != pred) l_G++;
                    if (popc(adm) == 2) l_a2++; else if (popc(adm) == 3) l_a3++; else l_ao++;
                }
                for (int p = 0; p < m; p++) {
                    l_cons++;
                    uint32_t adj[MAXN] = {0}; memcpy(adj, pa, sizeof(uint32_t) * m);
                    adj[p] |= 1u << m; adj[m] = ((1u << m) - 1) & ~(1u << p);
                    uint64_t h; int cyc[MAXN];
                    if (n <= DIRECT) { h = hc_exact(adj, n, dp, cyc); if ((h == 1) != (e[p] == 1)) l_lm++; }
                    else { if (e[p] != 1) continue; h = hc_exact(adj, n, dp, cyc); if (h != 1) { l_rej++; continue; } }
                    if (h != 1) continue;
                    if (okU) {   /* (H): run length of the child and the three facts of its proof; sole-entry rule */
                        l_ch++;
                        uint32_t FT = (1u << n) - 1, NCT = 0;
                        for (int v = 0; v < n; v++) if (strong_dg(adj, n, FT & ~(1u << v))) NCT |= 1u << v;
                        int delta = (p == yU) ? 0 : (p == yUm) ? 1 : (p == yUmm && yUm_sole_entry) ? 2 : -100;
                        if (popc(NCT) != rU + delta) l_H++;
                        if (((NCU & ~(1u << p)) & ~NCT) != 0) l_Hf++;          /* NC(U) - p stays non-critical */
                        if (NCT >> p & 1u) l_Hf++;                              /* p becomes critical */
                        if (!(NCT >> m & 1u)) l_Hf++;                           /* the new vertex is non-critical */
                        if (xU == p || popc(adj[xU]) != 1) l_Hf++;              /* x_U still a sole exit */
                        if (NCT >> ((xU + 1) % m) & 1u) l_Hf++;                 /* x_U+ critical */
                        if (delta == 2 && !(NCT >> yUm & 1u)) l_Hf++;           /* y_U- becomes non-critical */
                        int seT = 0; for (int v = 0; v < n; v++) if (n - 1 - popc(adj[v]) == 1) seT++;
                        int pSoleU = (m - 1 - popc(pa[p])) == 1;
                        if (seT != 1 + pSoleU) l_se++;
                    }
                    uint32_t cc[MAXN] = {0}, can[MAXN] = {0}; to_cycle_coords(adj, n, cyc, cc);
                    int arg; u128 code = canon_rot(cc, n, NULL, &arg); rotate_adj(cc, n, arg, can);
                    #pragma omp critical(candadd)
                    { if (ncand == capc) { capc = capc ? 2 * capc : 4096; cand = realloc(cand, sizeof(cls_t) * capc); }
                      cand[ncand].code = code; memset(cand[ncand].adj, 0, sizeof(cand[ncand].adj)); memcpy(cand[ncand].adj, can, sizeof(uint32_t) * n); ncand++; }
                }
            }
            #pragma omp critical(redg)
            { lemma_mismatch += l_lm; hc_reject_after_filter += l_rej; considered += l_cons;
              G_bad += l_G; adm2 += l_a2; adm3 += l_a3; admother += l_ao; H_bad += l_H; Hfact_bad += l_Hf; serule_bad += l_se; children += l_ch; }
            free(dp); free(dpe);
        }
        qsort(cand, ncand, sizeof(cls_t), cmp_cls);
        long u = 0, dups = 0; for (long i = 0; i < ncand; i++) { if (i == 0 || cand[i].code != cand[i - 1].code) cand[u++] = cand[i]; else dups++; }
        level_t cur; cur.n = n; cur.cnt = u; cur.v = cand;
        double tgen = omp_get_wtime() - t0;

        /* per-class checks */
        stats_t tot; memset(&tot, 0, sizeof tot);
        long *exit_hit = NULL; long two_exit_classes = 0, bij_bad = 0;
        /* rooted pairs of the previous level: index k*(n-1)+pos */
        exit_hit = calloc((size_t)prev.cnt * m, sizeof(long));
        long r2cls = -1, r2count = 0;
        #pragma omp parallel
        {
            uint64_t *dp = malloc(sizeof(uint64_t) * ((size_t)1 << (n - 1)) * n);
            stats_t st; memset(&st, 0, sizeof st);
            long l_two = 0, l_bad = 0;
            #pragma omp for schedule(dynamic, 16)
            for (long k = 0; k < cur.cnt; k++) {
                const uint32_t *adj = cur.v[k].adj; int y, x;
                int r = check_class(adj, n, &prev, L1MAX, HCMAX_TY, dp, &st, &y, &x);
                if (r == 2) {
                    #pragma omp critical(r2)
                    { r2count++; r2cls = k; }
                }
                /* (F) bijection: two sole exits z -> z+ ; map to (T - z, z+) */
                int z = -1, nexit = 0;
                for (int v = 0; v < n; v++) if (popc(adj[v]) == 1) { nexit++; if (popc(adj[(v + 1) % n]) == 1) z = v; }
                if (nexit == 2) {
                    l_two++;
                    if (z < 0) { l_bad++; continue; }
                    uint32_t sub[MAXN], can2[MAXN]; int nl[MAXN];
                    if (!delete_in_cycle(adj, n, z, sub, nl)) { l_bad++; continue; }
                    int arg; u128 c = canon_rot(sub, m, NULL, &arg); rotate_adj(sub, m, arg, can2);
                    long j = find_cls(&prev, c);
                    if (j < 0) { l_bad++; continue; }
                    int ypos = (nl[(z + 1) % n] - arg + m) % m;     /* label of z+ in the canonical labelling of T - z */
                    if (popc(prev.v[j].adj[ypos]) != 1) { l_bad++; continue; }
                    #pragma omp atomic
                    exit_hit[j * m + ypos]++;
                }
            }
            #pragma omp critical(redc)
            { for (int b = 0; b < NFAIL; b++) tot.fail[b] += st.fail[b];
              for (int i = 0; i < MAXN + 2; i++) { tot.Nr[i] += st.Nr[i]; tot.Tk[i] += st.Tk[i]; }
              for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) tot.oi[i][j] += st.oi[i][j];
              tot.rooted += st.rooted; two_exit_classes += l_two; bij_bad += l_bad; }
            free(dp);
        }
        /* every rooted pair (T', sole exit) of order n-1 hit exactly once? */
        long rooted_prev = 0, hit_once = 0, hit_other = 0;
        for (long j = 0; j < prev.cnt; j++) for (int pos = 0; pos < m; pos++) {
            int isexit = popc(prev.v[j].adj[pos]) == 1;
            if (isexit) { rooted_prev++; if (exit_hit[j * m + pos] == 1) hit_once++; else hit_other++; }
            else if (exit_hit[j * m + pos]) hit_other++;
        }
        free(exit_hit);
        /* Las Vergnas class = A_n ? */
        int lv_ok = 0; uint64_t hcAn = 0;
        { uint32_t An[MAXN] = {0}; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && (j == i + 1 || i >= j + 2)) An[i] |= 1u << j;
          uint64_t *dp = malloc(sizeof(uint64_t) * ((size_t)1 << (n - 1)) * n); int cyc[MAXN];
          hcAn = hc_exact(An, n, dp, cyc); free(dp);
          if (hcAn == 1 && r2count == 1) { uint32_t cc[MAXN]; to_cycle_coords(An, n, cyc, cc); if (canon_rot(cc, n, NULL, NULL) == cur.v[r2cls].code) lv_ok = 1; } }
        /* cross-check with PART I and the labelled enumeration (score-sorted canonical codes) */
        const char *x1 = "n/a", *x2 = "n/a";
        if (n <= 9) {
            uint64_t *ss = malloc(sizeof(uint64_t) * cur.cnt);
            for (long k = 0; k < cur.cnt; k++) ss[k] = canon_ss(cur.v[k].adj, n, NULL);
            qsort(ss, cur.cnt, sizeof(uint64_t), cmp_u64);
            if (uh_ss[n]) x1 = (uh_ss_cnt[n] == cur.cnt && !memcmp(ss, uh_ss[n], sizeof(uint64_t) * cur.cnt)) ? "IDENTICAL" : "DIFFERENT";
            if (lab[n]) x2 = (labc[n] == cur.cnt && !memcmp(ss, lab[n], sizeof(uint64_t) * cur.cnt)) ? "IDENTICAL" : "DIFFERENT";
            free(ss);
        }
        for (int r = 0; r <= n; r++) NN[n][r] = tot.Nr[r];

        printf("n=%2d Shelor candidates %ld, accepted %ld, distinct classes %ld (duplicates %ld) Garey F(2n-6)=%lld %s | hc vs e_p lemma mismatches %ld, filter-accepted but hc!=1: %ld\n",
               n, considered, ncand, cur.cnt, dups, fib(2 * n - 6), cur.cnt == fib(2 * n - 6) ? "OK" : "MISMATCH", lemma_mismatch, hc_reject_after_filter);
        printf("   same class set as PART I (no-assumption): %s ; as labelled enumeration: %s\n", x1, x2);
        if (m >= 4)
            printf("   (G) on all %ld parents of order %d: adm != {y,y-,(y-- iff y- sole entry)}: %ld ; |adm|=2: %ld (F(2m-8)=%lld) |adm|=3: %ld (F(2m-7)=%lld) other: %ld\n"
                   "   (H) on all %ld children: run length != r(U)+0/1/2: %ld ; proof facts violated: %ld ; sole-entry rule violated: %ld\n",
                   prev.cnt, m, G_bad, adm2, fib(2 * m - 8), adm3, fib(2 * m - 7), admother, children, H_bad, Hfact_bad, serule_bad);
        printf("   failures:"); for (int b = 0; b < NFAIL; b++) printf(" %s=%ld", fname[b], tot.fail[b]); printf("\n");
        printf("   N(n,r) r=2..n-1:"); for (int r = 2; r <= n - 1; r++) printf(" %ld", tot.Nr[r]); printf("  | outside: %ld\n", tot.Nr[0] + tot.Nr[1] + tot.Nr[n]);
        printf("   cut arcs k=2..n-1:"); for (int k = 2; k <= n - 1; k++) printf(" %ld", tot.Tk[k]);
        if (n >= 5 && n <= 13) { int ok = 1; for (int k = 2; k <= n - 1; k++) if (tot.Tk[k] != HYP9458[n][k - 2]) ok = 0; printf("  [HYP-9458 table row: %s]", ok ? "equal" : "DIFFERENT"); }
        if (n == 16) { int ok = 1; for (int k = 2; k <= n - 1; k++) if (tot.Tk[k] != HYP9458_16[k - 2]) ok = 0; printf("  [S74 order-16 row: %s]", ok ? "equal" : "DIFFERENT"); }
        printf("\n");
        { int ok = 1; for (int k = 2; k <= n - 1; k++) if (tot.Tk[k] != tot.Nr[n + 1 - k]) ok = 0; printf("   T(n,k) = N(n,n+1-k): %s", ok ? "yes" : "NO"); }
        if (n >= 5 && n <= 13) { int ok = 1; for (int r = 2; r <= n - 1; r++) if (tot.Nr[r] != NOTE_N[n][r - 2]) ok = 0; printf(" ; note's section-4 row: %s", ok ? "equal" : "DIFFERENT"); }
        if (n >= 14) { int ok = 1; for (int r = 2; r <= n - 1; r++) if (tot.Nr[r] != P7pred[n][r - 2]) ok = 0; printf(" ; HYP-9460 P7 predicted row: %s", ok ? "equal" : "DIFFERENT"); }
        if (n >= 6) { int ok = 1; for (int r = 2; r <= n - 1; r++) { long rhs = NN[n - 1][r] + 2 * NN[n - 1][r - 1] - NN[n - 2][r - 1]; if (NN[n][r] != rhs) ok = 0; }
                      printf(" ; P7 recurrence: %s", ok ? "holds" : "FAILS"); }
        printf("\n");
        printf("   (#se,#sen): (1,1)=%ld (1,2)=%ld (2,1)=%ld (2,2)=%ld other=%ld ; formula %lld %lld %lld %lld\n",
               tot.oi[1][1], tot.oi[1][2], tot.oi[2][1], tot.oi[2][2], cur.cnt - tot.oi[1][1] - tot.oi[1][2] - tot.oi[2][1] - tot.oi[2][2],
               fib(2 * n - 10) + 1, fib(2 * n - 9) - 1, fib(2 * n - 9) - 1, fib(2 * n - 8) + 1);
        printf("   one/two sole exits: %ld (F(2n-8)=%lld) / %ld (F(2n-7)=%lld) ; rooted (T, sole exit) %ld (F(2n-5)=%lld)\n",
               tot.oi[1][1] + tot.oi[1][2], fib(2 * n - 8), tot.oi[2][1] + tot.oi[2][2], fib(2 * n - 7), tot.rooted, fib(2 * n - 5));
        printf("   (F) bijection: two-sole-exit classes %ld, map failures %ld, rooted pairs of order %d: %ld, hit exactly once %ld, other %ld\n",
               two_exit_classes, bij_bad, m, rooted_prev, hit_once, hit_other);
        printf("   Las Vergnas/top diagonal: N(n,2)=%ld, that class = A_n: %s (hc(A_n)=%llu)\n", r2count, lv_ok ? "yes" : "NO", (unsigned long long)hcAn);
        printf("   time: generation %.1f s, checks %.1f s\n", tgen, omp_get_wtime() - t0 - tgen);
        fflush(stdout);
        free(prev.v); prev = cur;
    }
    printf("total time %.1f s\n", omp_get_wtime() - T0);
    return 0;
}
