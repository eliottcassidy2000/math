/* Independent audit code for THM-4870 (A), (B) - kind-pasteur-2026-11-04-S10 audit.
 *
 * Written from scratch for the audit; shares no code with the author's programs.
 *
 * Modes
 *   full  NMAX            : all tournament classes of order 1..NMAX (NMAX <= 9) by vertex extension +
 *                           canonical form (individualisation-refinement, IR).  Per class: H (Hamiltonian
 *                           paths, subset DP), hc (Hamiltonian cycles), strong?, |Aut| by IR leaf count AND by an
 *                           independent backtracking enumerator (cross-checked), cyclicity of Aut, and for
 *                           uniquely Hamiltonian classes r (#non-critical) and k (#cut arcs of the cycle).
 *                           Writes one line per STRONG class to strong_classes_n<=NMAX.txt.
 *   strongB NMAX B        : strong classes of order 3..NMAX with H <= B only, generated as strong one-vertex
 *                           extensions of strong classes with H <= B (valid by Moon 1966 + H(T) >= H(T-v)).
 *   uh NMAX               : uniquely Hamiltonian classes of order 3..NMAX, generated as all one-vertex
 *                           extensions of UH classes (valid by Moon 1966 + hc(T) >= hc(T-v) when T-v strong).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 12
typedef uint64_t u64;
typedef uint32_t u32;

/* ---------------- basic tournament functions ---------------- */
static __thread u64 *DP;

static u64 hpaths(int n, const u32 *out) {
    int full = (1 << n) - 1;
    memset(DP, 0, sizeof(u64) * ((size_t)1 << n) * n);
    for (int v = 0; v < n; v++) DP[((size_t)1 << v) * n + v] = 1;
    for (int m = 1; m <= full; m++) {
        for (int v = 0; v < n; v++) {
            u64 c = DP[(size_t)m * n + v];
            if (!c) continue;
            u32 av = out[v] & ~(u32)m;
            while (av) {
                int w = __builtin_ctz(av); av &= av - 1;
                DP[(size_t)(m | (1 << w)) * n + w] += c;
            }
        }
    }
    u64 s = 0;
    for (int v = 0; v < n; v++) s += DP[(size_t)full * n + v];
    return s;
}

static u64 hcycles(int n, const u32 *out) { /* directed Hamiltonian cycles, n >= 3 */
    if (n < 3) return 0;
    int full = (1 << n) - 1;
    memset(DP, 0, sizeof(u64) * ((size_t)1 << n) * n);
    DP[(size_t)1 * n + 0] = 1;
    for (int m = 1; m <= full; m += 2) { /* masks containing vertex 0 */
        for (int v = 0; v < n; v++) {
            u64 c = DP[(size_t)m * n + v];
            if (!c) continue;
            u32 av = out[v] & ~(u32)m;
            while (av) {
                int w = __builtin_ctz(av); av &= av - 1;
                DP[(size_t)(m | (1 << w)) * n + w] += c;
            }
        }
    }
    u64 s = 0;
    for (int v = 1; v < n; v++) if (out[v] & 1u) s += DP[(size_t)full * n + v];
    return s;
}

/* strong connectivity of a digraph given by out-sets on vertex set 'vs' (bitmask) */
static int strong_on(int n, const u32 *out, u32 vs) {
    if (!vs) return 1;
    int s = __builtin_ctz(vs);
    u32 r = 1u << s, prev = 0;
    while (r != prev) { prev = r; u32 t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; r |= out[v] & vs; } }
    if (r != vs) return 0;
    /* backward */
    u32 inn[MAXN];
    for (int v = 0; v < n; v++) inn[v] = 0;
    for (int u = 0; u < n; u++) if (vs >> u & 1) { u32 t = out[u] & vs; while (t) { int v = __builtin_ctz(t); t &= t - 1; inn[v] |= 1u << u; } }
    r = 1u << s; prev = 0;
    while (r != prev) { prev = r; u32 t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; r |= inn[v] & vs; } }
    return r == vs;
}

/* ---------------- canonical form by individualisation-refinement ---------------- */
static __thread int cmp_sig_n;
typedef struct { int key[MAXN + 2]; int v; } Sig;
static int sigcmp(const void *a, const void *b) {
    const Sig *x = a, *y = b;
    for (int i = 0; i < cmp_sig_n; i++) if (x->key[i] != y->key[i]) return x->key[i] < y->key[i] ? -1 : 1;
    return 0;
}
/* refine colouring col[] (values 0..k-1) to the coarsest equitable refinement; returns #colours */
static int refine(int n, const u32 *out, int *col) {
    int k = 0;
    for (int v = 0; v < n; v++) if (col[v] + 1 > k) k = col[v] + 1;
    for (;;) {
        Sig s[MAXN];
        cmp_sig_n = k + 1;
        for (int v = 0; v < n; v++) {
            s[v].v = v;
            for (int i = 0; i <= k; i++) s[v].key[i] = 0;
            s[v].key[0] = col[v];
            u32 t = out[v];
            while (t) { int w = __builtin_ctz(t); t &= t - 1; s[v].key[1 + col[w]]++; }
        }
        qsort(s, n, sizeof(Sig), sigcmp);
        int nc = 0, newcol[MAXN];
        for (int i = 0; i < n; i++) {
            if (i > 0 && sigcmp(&s[i - 1], &s[i]) != 0) nc++;
            newcol[s[i].v] = nc;
        }
        nc++;
        for (int v = 0; v < n; v++) col[v] = newcol[v];
        if (nc == k) return k;
        k = nc;
    }
}
static inline int pidx(int n, int i, int j) { return i * n - i * (i + 1) / 2 + (j - i - 1); }
static u64 code_of(int n, const u32 *out, const int *lab /* vertex -> label */) {
    int inv[MAXN];
    for (int v = 0; v < n; v++) inv[lab[v]] = v;
    u64 c = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++)
        if (out[inv[i]] >> inv[j] & 1u) c |= (u64)1 << pidx(n, i, j);
    return c;
}
static __thread u64 best_code; static __thread int best_cnt; static __thread int best_lab[MAXN];
static void ir_search(int n, const u32 *out, const int *col0) {
    int col[MAXN];
    memcpy(col, col0, sizeof(int) * n);
    int k = refine(n, out, col);
    if (k == n) {
        u64 c = code_of(n, out, col);
        if (best_cnt == 0 || c < best_code) { best_code = c; best_cnt = 1; memcpy(best_lab, col, sizeof(int) * n); }
        else if (c == best_code) best_cnt++;
        return;
    }
    /* target cell: smallest colour with >= 2 vertices */
    int sz[MAXN] = {0};
    for (int v = 0; v < n; v++) sz[col[v]]++;
    int tc = -1;
    for (int c = 0; c < k; c++) if (sz[c] >= 2) { tc = c; break; }
    for (int v = 0; v < n; v++) if (col[v] == tc) {
        int c2[MAXN];
        /* individualise v: v gets colour tc, rest of the cell tc+1, colours above shift by 1 */
        for (int u = 0; u < n; u++) {
            if (col[u] < tc) c2[u] = col[u];
            else if (col[u] == tc) c2[u] = (u == v) ? tc : tc + 1;
            else c2[u] = col[u] + 1;
        }
        ir_search(n, out, c2);
    }
}
/* canonical code and |Aut| (= number of leaves with the minimal code, full tree, no pruning) */
static u64 canon(int n, const u32 *out, int *aut_ir) {
    int col[MAXN];
    for (int v = 0; v < n; v++) col[v] = 0;
    best_cnt = 0;
    ir_search(n, out, col);
    if (aut_ir) *aut_ir = best_cnt;
    return best_code;
}
static void decode(int n, u64 c, u32 *out) {
    for (int v = 0; v < n; v++) out[v] = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        if (c >> pidx(n, i, j) & 1) out[i] |= 1u << j; else out[j] |= 1u << i;
    }
}

/* ---------------- independent automorphism enumerator (backtracking) ---------------- */
#define MAXAUT 4096
static int AUTS[MAXAUT][MAXN]; static int naut;
static int bt_map[MAXN], bt_used;
static void bt(int n, const u32 *out, int i) {
    if (i == n) { if (naut < MAXAUT) memcpy(AUTS[naut], bt_map, sizeof(int) * n); naut++; return; }
    int di = __builtin_popcount(out[i]);
    for (int c = 0; c < n; c++) {
        if (bt_used >> c & 1) continue;
        if (__builtin_popcount(out[c]) != di) continue;
        int ok = 1;
        for (int j = 0; j < i && ok; j++) {
            int a = (out[i] >> j) & 1, b = (out[c] >> bt_map[j]) & 1;
            if (a != b) ok = 0;
        }
        if (!ok) continue;
        bt_map[i] = c; bt_used |= 1 << c;
        bt(n, out, i + 1);
        bt_used &= ~(1 << c);
    }
}
static int aut_bt(int n, const u32 *out) { naut = 0; bt_used = 0; bt(n, out, 0); return naut; }
static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static int perm_order(int n, const int *p) {
    int seen = 0, o = 1;
    for (int i = 0; i < n; i++) if (!(seen >> i & 1)) {
        int L = 0, j = i;
        while (!(seen >> j & 1)) { seen |= 1 << j; j = p[j]; L++; }
        o = o / gcd(o, L) * L;
    }
    return o;
}
static int aut_is_cyclic(int n) { /* uses AUTS from the last aut_bt call */
    int m = naut < MAXAUT ? naut : MAXAUT, best = 1;
    for (int a = 0; a < m; a++) { int o = perm_order(n, AUTS[a]); if (o > best) best = o; }
    return best == naut;
}

/* ---------------- hash set of u64 codes ---------------- */
typedef struct { u64 *t; size_t cap, cnt; } HSet;
static void hs_init(HSet *h, size_t cap) { h->cap = cap; h->cnt = 0; h->t = calloc(cap, sizeof(u64)); }
static int hs_add(HSet *h, u64 key) { /* returns 1 if newly inserted */
    u64 k = key + 1; /* 0 = empty */
    size_t i = (size_t)((k * 0x9E3779B97F4A7C15ULL) >> 20) & (h->cap - 1);
    while (h->t[i]) { if (h->t[i] == k) return 0; i = (i + 1) & (h->cap - 1); }
    h->t[i] = k; h->cnt++;
    if (h->cnt * 2 > h->cap) { fprintf(stderr, "hash too full\n"); exit(1); }
    return 1;
}

/* ---------------- uniquely Hamiltonian analysis ---------------- */
static int cyc[MAXN], found_cyc;
static void find_cycle_dfs(int n, const u32 *out, int pos, u32 used) {
    if (found_cyc) return;
    if (pos == n) { if (out[cyc[n - 1]] & 1u) found_cyc = 1; return; }
    u32 av = out[cyc[pos - 1]] & ~used;
    while (av && !found_cyc) { int w = __builtin_ctz(av); av &= av - 1; cyc[pos] = w; find_cycle_dfs(n, out, pos + 1, used | (1u << w)); }
}
/* returns r, k; sets *onerun = 1 if NC is one cyclic run */
static void uh_rk(int n, const u32 *out, int *r, int *k, int *onerun) {
    found_cyc = 0; cyc[0] = 0; find_cycle_dfs(n, out, 1, 1u);
    if (!found_cyc) { fprintf(stderr, "no cycle?\n"); exit(1); }
    int full = (1 << n) - 1, nc[MAXN];
    *r = 0;
    for (int i = 0; i < n; i++) { nc[i] = strong_on(n, out, full & ~(1u << cyc[i])); *r += nc[i]; }
    *k = 0;
    for (int i = 0; i < n; i++) {
        u32 o2[MAXN]; memcpy(o2, out, sizeof(u32) * n);
        o2[cyc[i]] &= ~(1u << cyc[(i + 1) % n]);
        if (!strong_on(n, o2, full)) (*k)++;
    }
    /* one run: number of i with nc[i] && !nc[i-1] equals 1 (and not all NC) */
    int starts = 0;
    for (int i = 0; i < n; i++) if (nc[i] && !nc[(i + n - 1) % n]) starts++;
    *onerun = (starts == 1);
}

/* ---------------- modes ---------------- */
static void mode_full(int NMAX) {
    u64 *cur = malloc(sizeof(u64)), *nxt; size_t ncur = 1; cur[0] = 0; /* order 1 */
    FILE *fs = fopen("strong_classes.txt", "w");
    fprintf(fs, "# n H hc aut cyclic code(hex)   [strong classes only, orders 1..%d]\n", NMAX);
    fprintf(fs, "1 1 0 1 1 0\n");
    printf("order 1: classes 1\n");
    for (int n = 2; n <= NMAX; n++) {
        HSet hs; hs_init(&hs, (size_t)1 << 21);
        size_t cap = 1024; nxt = malloc(sizeof(u64) * cap); size_t nn = 0;
        for (size_t c = 0; c < ncur; c++) {
            u32 o[MAXN]; decode(n - 1, cur[c], o);
            for (u32 p = 0; p < (1u << (n - 1)); p++) {
                u32 t[MAXN];
                for (int v = 0; v < n - 1; v++) t[v] = o[v] | (((p >> v) & 1) ? 0 : (1u << (n - 1)));
                t[n - 1] = p;
                u64 cd = canon(n, t, NULL);
                if (hs_add(&hs, cd)) { if (nn == cap) { cap *= 2; nxt = realloc(nxt, sizeof(u64) * cap); } nxt[nn++] = cd; }
            }
        }
        free(hs.t); free(cur); cur = nxt; ncur = nn;
        /* statistics */
        size_t nstrong = 0, autmismatch = 0, uh = 0, uh_bad = 0, uh_asym = 0, autdivH_fail = 0;
        long autcount[64] = {0}; long autcount_strong[64] = {0};
        u64 minHstrong = (u64)-1;
        long Nr[MAXN + 1] = {0}; long sumr_num = 0;
        for (size_t c = 0; c < ncur; c++) {
            u32 t[MAXN]; decode(n, cur[c], t);
            int air; u64 cc = canon(n, t, &air);
            if (cc != cur[c]) { fprintf(stderr, "canon not idempotent\n"); exit(1); }
            int ab = aut_bt(n, t);
            if (ab != air) autmismatch++;
            int cyclic = aut_is_cyclic(n);
            u64 H = hpaths(n, t);
            if (H % (u64)ab) autdivH_fail++;
            int st = (n == 1) ? 1 : strong_on(n, t, (1u << n) - 1);
            if (ab < 64) autcount[ab]++;
            if (st) {
                nstrong++; if (ab < 64) autcount_strong[ab]++;
                if (H < minHstrong) minHstrong = H;
                u64 hc = hcycles(n, t);
                fprintf(fs, "%d %llu %llu %d %d %llx\n", n, (unsigned long long)H, (unsigned long long)hc, ab, cyclic,
                        (unsigned long long)cur[c]);
                if (hc == 1 && n >= 3) {
                    uh++;
                    if (ab == 1) uh_asym++;
                    if (n >= 4) {
                        int r, k, one; uh_rk(n, t, &r, &k, &one);
                        if (r + k != n + 1 || !one || r < 2 || r > n - 1) uh_bad++;
                        Nr[r]++; sumr_num += r;
                    }
                }
            }
        }
        printf("order %d: classes %zu  strong %zu  min H(strong) %llu  |Aut| IR-vs-backtrack mismatches %zu  "
               "|Aut| does not divide H: %zu\n", n, ncur, nstrong, (unsigned long long)minHstrong, autmismatch, autdivH_fail);
        printf("   |Aut| distribution (all classes):");
        for (int a = 1; a < 64; a++) if (autcount[a]) printf(" %d:%ld", a, autcount[a]);
        printf("\n   |Aut| distribution (strong classes):");
        for (int a = 1; a < 64; a++) if (autcount_strong[a]) printf(" %d:%ld", a, autcount_strong[a]);
        printf("\n   UH classes (hc=1): %zu, asymmetric %zu", uh, uh_asym);
        if (n >= 4) {
            printf(", violations of [r+k=n+1, one run, 2<=r<=n-1]: %zu\n   N(n,r) r=2..n-1:", uh_bad);
            for (int r = 2; r <= n - 1; r++) printf(" %ld", Nr[r]);
            printf("   E[r] = %ld/%zu = %.12f", sumr_num, uh, uh ? (double)sumr_num / uh : 0.0);
        }
        printf("\n");
        fflush(stdout);
    }
    fclose(fs);
}

/* strong classes with H <= B, by strong extension of strong classes */
static void mode_strongB(int NMAX, u64 B) {
    u64 *cur = malloc(sizeof(u64)); size_t ncur = 1;
    { u32 c3[3] = {2, 4, 1}; cur[0] = canon(3, c3, NULL); }
    FILE *fs = fopen("strongB_classes.txt", "w");
    fprintf(fs, "# n H aut cyclic code(hex)  strong classes with H <= %llu, orders 3..%d (plus order 1)\n", (unsigned long long)B, NMAX);
    fprintf(fs, "1 1 1 1 0\n");
    { u32 c3[3]; decode(3, cur[0], c3); int a = aut_bt(3, c3); fprintf(fs, "3 %llu %d %d %llx\n", (unsigned long long)hpaths(3, c3), a, aut_is_cyclic(3), (unsigned long long)cur[0]); }
    printf("order 3: strong classes with H<=%llu: 1\n", (unsigned long long)B);
    for (int n = 4; n <= NMAX; n++) {
        HSet hs; hs_init(&hs, (size_t)1 << 22);
        size_t cap = 1024, nn = 0; u64 *nxt = malloc(sizeof(u64) * cap);
        for (size_t c = 0; c < ncur; c++) {
            u32 o[MAXN]; decode(n - 1, cur[c], o);
            for (u32 p = 1; p < (1u << (n - 1)) - 1; p++) { /* p = 0 or all: new vertex is a sink/source, not strong */
                u32 t[MAXN];
                for (int v = 0; v < n - 1; v++) t[v] = o[v] | (((p >> v) & 1) ? 0 : (1u << (n - 1)));
                t[n - 1] = p;
                if (!strong_on(n, t, (1u << n) - 1)) continue;
                u64 H = hpaths(n, t);
                if (H > B) continue;
                u64 cd = canon(n, t, NULL);
                if (hs_add(&hs, cd)) { if (nn == cap) { cap *= 2; nxt = realloc(nxt, sizeof(u64) * cap); } nxt[nn++] = cd; }
            }
        }
        free(hs.t); free(cur); cur = nxt; ncur = nn;
        u64 minH = (u64)-1; long autc[64] = {0}; size_t mism = 0;
        for (size_t c = 0; c < ncur; c++) {
            u32 t[MAXN]; decode(n, cur[c], t);
            int air; canon(n, t, &air);
            int a = aut_bt(n, t); if (a != air) mism++;
            int cy = aut_is_cyclic(n);
            u64 H = hpaths(n, t); if (H < minH) minH = H;
            if (a < 64) autc[a]++;
            fprintf(fs, "%d %llu %d %d %llx\n", n, (unsigned long long)H, a, cy, (unsigned long long)cur[c]);
        }
        printf("order %d: strong classes with H<=%llu: %zu  (min H %llu)  |Aut| mismatches %zu  |Aut| dist:",
               n, (unsigned long long)B, ncur, ncur ? (unsigned long long)minH : 0ULL, mism);
        for (int a = 1; a < 64; a++) if (autc[a]) printf(" %d:%ld", a, autc[a]);
        printf("\n"); fflush(stdout);
        if (ncur == 0) { printf("no strong classes with H<=%llu at order %d; by Moon + monotonicity none at any larger order\n", (unsigned long long)B, n); break; }
    }
    fclose(fs);
}

/* UH classes by extension of UH classes */
static void mode_uh(int NMAX) {
    u64 *cur = malloc(sizeof(u64)); size_t ncur = 1;
    { u32 c3[3] = {2, 4, 1}; cur[0] = canon(3, c3, NULL); }
    FILE *fu = fopen("uh_classes.txt", "w");
    printf("order 3: UH classes 1 (C3, |Aut| = 3)\n");
    for (int n = 4; n <= NMAX; n++) {
        HSet hs; hs_init(&hs, (size_t)1 << 20);
        size_t cap = 1024, nn = 0; u64 *nxt = malloc(sizeof(u64) * cap);
        for (size_t c = 0; c < ncur; c++) {
            u32 o[MAXN]; decode(n - 1, cur[c], o);
            for (u32 p = 1; p < (1u << (n - 1)) - 1; p++) {
                u32 t[MAXN];
                for (int v = 0; v < n - 1; v++) t[v] = o[v] | (((p >> v) & 1) ? 0 : (1u << (n - 1)));
                t[n - 1] = p;
                if (hcycles(n, t) != 1) continue;
                u64 cd = canon(n, t, NULL);
                if (hs_add(&hs, cd)) { if (nn == cap) { cap *= 2; nxt = realloc(nxt, sizeof(u64) * cap); } nxt[nn++] = cd; }
            }
        }
        free(hs.t); free(cur); cur = nxt; ncur = nn;
        long Nr[MAXN + 1] = {0}; long sumr = 0; size_t bad = 0, asym = 0;
        for (size_t c = 0; c < ncur; c++) {
            u32 t[MAXN]; decode(n, cur[c], t);
            int a = aut_bt(n, t); if (a == 1) asym++;
            int r, k, one; uh_rk(n, t, &r, &k, &one);
            if (r + k != n + 1 || !one || r < 2 || r > n - 1) bad++;
            Nr[r]++; sumr += r;
            fprintf(fu, "%d %llu %d %d %d %llx\n", n, (unsigned long long)hpaths(n, t), a, r, k, (unsigned long long)cur[c]);
        }
        printf("order %d: UH classes %zu  asymmetric %zu  violations of [r+k=n+1, one run, 2<=r<=n-1] %zu\n   N(n,r) r=2..n-1:",
               n, ncur, asym, bad);
        for (int r = 2; r <= n - 1; r++) printf(" %ld", Nr[r]);
        printf("   E[r] = %ld/%zu = %.15f\n", sumr, ncur, (double)sumr / ncur);
        fflush(stdout);
    }
    fclose(fu);
}

/* ---------------- parallel strongB (this file only) ---------------- */
static int cmp_u64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
static void mode_strongB_par(int NMAX, u64 B) {
    u64 *cur = malloc(sizeof(u64)); size_t ncur = 1;
    { u32 c3[3] = {2, 4, 1}; cur[0] = canon(3, c3, NULL); }
    FILE *fs = fopen("strongB1125_classes.txt", "w");
    fprintf(fs, "# n H aut cyclic code(hex)  strong classes with H <= %llu, orders 3..%d (plus order 1)\n", (unsigned long long)B, NMAX);
    fprintf(fs, "1 1 1 1 0\n");
    { u32 c3[3]; decode(3, cur[0], c3); int a = aut_bt(3, c3); fprintf(fs, "3 %llu %d %d %llx\n", (unsigned long long)hpaths(3, c3), a, aut_is_cyclic(3), (unsigned long long)cur[0]); }
    printf("order 3: strong classes with H<=%llu: 1\n", (unsigned long long)B);
    for (int n = 4; n <= NMAX; n++) {
        size_t tot = 0, cap_all = 1 << 20; u64 *all = malloc(sizeof(u64) * cap_all);
        #pragma omp parallel
        {
            if (!DP) DP = malloc(sizeof(u64) * ((size_t)1 << MAXN) * MAXN);
            size_t lcap = 1 << 14, ln = 0; u64 *loc = malloc(sizeof(u64) * lcap);
            #pragma omp for schedule(dynamic, 16)
            for (long c = 0; c < (long)ncur; c++) {
                u32 o[MAXN]; decode(n - 1, cur[c], o);
                for (u32 p = 1; p < (1u << (n - 1)) - 1; p++) {
                    u32 t[MAXN];
                    for (int v = 0; v < n - 1; v++) t[v] = o[v] | (((p >> v) & 1) ? 0 : (1u << (n - 1)));
                    t[n - 1] = p;
                    if (!strong_on(n, t, (1u << n) - 1)) continue;
                    if (hpaths(n, t) > B) continue;
                    u64 cd = canon(n, t, NULL);
                    if (ln == lcap) { lcap *= 2; loc = realloc(loc, sizeof(u64) * lcap); }
                    loc[ln++] = cd;
                }
                if (ln > (1u << 16)) { /* compact locally */
                    qsort(loc, ln, sizeof(u64), cmp_u64);
                    size_t k = 0; for (size_t i = 0; i < ln; i++) if (i == 0 || loc[i] != loc[i - 1]) loc[k++] = loc[i];
                    ln = k;
                }
            }
            #pragma omp critical
            {
                if (tot + ln > cap_all) { while (tot + ln > cap_all) cap_all *= 2; all = realloc(all, sizeof(u64) * cap_all); }
                memcpy(all + tot, loc, sizeof(u64) * ln); tot += ln;
            }
            free(loc);
        }
        qsort(all, tot, sizeof(u64), cmp_u64);
        size_t k = 0; for (size_t i = 0; i < tot; i++) if (i == 0 || all[i] != all[i - 1]) all[k++] = all[i];
        free(cur); cur = all; ncur = k;
        u64 minH = (u64)-1; long autc[128] = {0}; size_t mism = 0;
        for (size_t c = 0; c < ncur; c++) {
            u32 t[MAXN]; decode(n, cur[c], t);
            int air; canon(n, t, &air);
            int a = aut_bt(n, t); if (a != air) mism++;
            int cy = aut_is_cyclic(n);
            u64 H = hpaths(n, t); if (H < minH) minH = H;
            if (a < 128) autc[a]++;
            fprintf(fs, "%d %llu %d %d %llx\n", n, (unsigned long long)H, a, cy, (unsigned long long)cur[c]);
        }
        printf("order %d: strong classes with H<=%llu: %zu  (min H %llu)  |Aut| IR/backtrack mismatches %zu  |Aut| dist:",
               n, (unsigned long long)B, ncur, ncur ? (unsigned long long)minH : 0ULL, mism);
        for (int a = 1; a < 128; a++) if (autc[a]) printf(" %d:%ld", a, autc[a]);
        printf("\n"); fflush(stdout);
    }
    fclose(fs);
}

int main(int argc, char **argv) {
    DP = malloc(sizeof(u64) * ((size_t)1 << MAXN) * MAXN);
    if (argc >= 4 && !strcmp(argv[1], "strongBpar")) mode_strongB_par(atoi(argv[2]), strtoull(argv[3], 0, 10));
    else { fprintf(stderr, "usage: strongBpar NMAX B\n"); return 1; }
    return 0;
}
