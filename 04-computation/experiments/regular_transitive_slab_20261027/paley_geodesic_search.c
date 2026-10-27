/* paley_geodesic_search.c -- mac-mini-2026-10-27-slab (orders up to 26; mode 2 = the Paley tournament P_n, n prime = 3 mod 4)
 * Is a given regular tournament R of odd order n exactly c = (n^2-1)/8 native moves from the transitive class?
 * Equivalently (THM-4819): is there a chain of c L-descending moves from R to a transitive tournament?  Each
 * L-descending move from a slab state stays in the slab.  A move lowers L iff its endpoint pair (a, b) has s(a) above
 * and s(b) below its sorted target in some optimal matching (THM-4815 B1, B2), so a state's candidate moves are the
 * Hamiltonian paths joining such pairs (endpoint relation by subset DP; paths sampled by backtracking).
 * Depth-first search with lookahead: candidate images are ranked by (no source/sink first, then the number of
 * descending pairs their own endpoint relation realises); failures are memoised by canonical form.  A found word is
 * replayed and checked; a failure within the node budget is reported as UNDECIDED, never as a counterexample.
 * Instances: mode 0 = random regular tournaments (rotational circulant, then many random 3-cycle reversals, which
 * preserve scores); mode 1 = every circulant on Z_n (one per multiplier class is enough, all are tried).
 * Build: cc -O3 -march=native -I/opt/homebrew/opt/nauty/include regular_geodesic_search.c /opt/homebrew/opt/nauty/lib/libnauty.a -o regular_geodesic_search
 * Run:   ./regular_geodesic_search n mode count seed [node_budget] [cands]
 */
#define MAXN 32
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n, m, C;
static uint64_t rs;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
typedef struct { uint32_t o[32]; } Tour;
static uint32_t *dp;
static void fill_dp(const uint32_t *o) {
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(uint32_t) * (size_t)(full + 1) * n);
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1u << v;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { uint32_t st = dp[(size_t)S * n + v]; if (!st) continue;
        uint32_t nx = o[v] & ~S; while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[(size_t)(S | 1 << w) * n + w] |= st; } }
}
static int sample_path(const uint32_t *o, int a, int b, int *Q) {
    int full = (1 << n) - 1; if (!(dp[(size_t)full * n + b] >> a & 1)) return 0;
    int S = full, v = b, k = n - 1; Q[k] = b;
    while (k > 0) { int S2 = S & ~(1 << v); int cand[32], nc = 0;
        for (int u = 0; u < n; u++) if ((S2 >> u & 1) && (dp[(size_t)S2 * n + u] >> a & 1) && (o[u] >> v & 1)) cand[nc++] = u;
        if (!nc) return 0; int u = cand[rnd() % nc]; Q[--k] = u; S = S2; v = u; }
    return 1;
}
static void apply_move(const uint32_t *o, const int *Q, uint32_t *r) { uint32_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) r[u] = (~o[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; r[u] |= 1u << v; r[v] &= ~(1u << u); } }
static int is_path(const uint32_t *o, const int *Q) { uint32_t seen = 0; for (int i = 0; i < n; i++) { if (seen >> Q[i] & 1) return 0; seen |= 1u << Q[i]; if (i && !(o[Q[i-1]] >> Q[i] & 1)) return 0; } return seen == (1u << n) - 1; }
static int Lof(const int *s0) { int s[32]; memcpy(s, s0, n * sizeof(int));
    for (int i = 1; i < n; i++) { int v = s[i], j = i - 1; while (j >= 0 && s[j] > v) { s[j + 1] = s[j]; j--; } s[j + 1] = v; }
    int l = 0; for (int i = 0; i < n; i++) l += abs(s[i] - i); return l / 2; }
static void scores(const uint32_t *o, int *s) { for (int v = 0; v < n; v++) s[v] = __builtin_popcount(o[v]); }
/* descending pairs realised by the current dp: list (a,b) */
static int desc_pairs(const uint32_t *o, int *pa, int *pb) {
    int s[32], t[32]; scores(o, s); int L0 = Lof(s), np = 0, full = (1 << n) - 1;
    for (int b = 0; b < n; b++) { uint32_t st = dp[(size_t)full * n + b];
        while (st) { int a = __builtin_ctz(st); st &= st - 1;
            for (int v = 0; v < n; v++) t[v] = n - 1 - s[v]; t[a]++; t[b]--;
            if (Lof(t) == L0 - 1) { pa[np] = a; pb[np] = b; np++; } } }
    return np;
}
typedef struct { uint64_t w[4]; } key128;   /* 256 bits: n <= 23 */
static key128 canon(const uint32_t *q) {
    graph g[MAXN], cg[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_DIGRAPH(options); statsblk stats; options.getcanon = TRUE;
    EMPTYGRAPH(g, 1, n);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (q[i] >> j & 1) ADDONEARC(g, i, j, 1);
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, cg);
    key128 k; memset(&k, 0, sizeof k); int pos = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (ISELEMENT(GRAPHROW(cg, i, 1), j)) k.w[pos >> 6] |= 1ULL << (pos & 63); pos++; }
    return k;
}
#define HC (1 << 20)
static key128 *hk; static unsigned char *hu;
static uint64_t khash(key128 k) { return ((k.w[0] * 0x9E3779B97F4A7C15ULL) ^ (k.w[1] * 0xC2B2AE3D27D4EB4FULL) ^ (k.w[2] * 0x165667B19E3779F9ULL) ^ (k.w[3] * 0x27D4EB2F165667C5ULL)) >> 40; }
static int fail_known(key128 k) { uint64_t h = khash(k) & (HC - 1);
    while (hu[h]) { if (!memcmp(&hk[h], &k, sizeof k)) return 1; h = (h + 1) & (HC - 1); } return 0; }
static void fail_add(key128 k) { uint64_t h = khash(k) & (HC - 1);
    while (hu[h]) { if (!memcmp(&hk[h], &k, sizeof k)) return; h = (h + 1) & (HC - 1); } hu[h] = 1; hk[h] = k; }
static long nodes, budget; static int cands;
static int word[128][32]; static int wlen;
static int has_srcsink(const uint32_t *o) { for (int v = 0; v < n; v++) { int s = __builtin_popcount(o[v]); if (s == 0 || s == n - 1) return 1; } return 0; }
static int search(const uint32_t *o, int depth) {
    int s[32]; scores(o, s);
    if (Lof(s) == 0) return 1;
    if (++nodes > budget) return -1;
    key128 k = canon(o); if (fail_known(k)) return 0;
    fill_dp(o);
    static int pa[1024], pb[1024]; int np = desc_pairs(o, pa, pb);
    if (!np) { fail_add(k); return 0; }
    /* candidates: images of sampled paths for random descending pairs, ranked by lookahead */
    Tour cand[64]; int qs[64][32]; int score[64]; int nc = 0;
    for (int t = 0; t < cands && nc < 64; t++) {
        int c = rnd() % np; int Q[32];
        if (!sample_path(o, pa[c], pb[c], Q)) continue;
        apply_move(o, Q, cand[nc].o); memcpy(qs[nc], Q, sizeof Q); nc++;
    }
    for (int i = 0; i < nc; i++) {   /* lookahead score */
        int ss = 0; if (!has_srcsink(cand[i].o)) ss += 1000000;
        int s2[32]; scores(cand[i].o, s2);
        if (Lof(s2) == 0) ss += 100000000;
        else { fill_dp(cand[i].o); static int pa2[1024], pb2[1024]; ss += desc_pairs(cand[i].o, pa2, pb2); }
        score[i] = ss;
    }
    for (int r = 0; r < nc; r++) {
        int best = -1; for (int i = 0; i < nc; i++) if (score[i] >= 0 && (best < 0 || score[i] > score[best])) best = i;
        if (best < 0) break; score[best] = -1;
        memcpy(word[depth], qs[best], sizeof(int) * 32);
        int res = search(cand[best].o, depth + 1);
        if (res == 1) { if (depth + 1 > wlen) wlen = depth + 1; return 1; }
        if (res < 0) return -1;
    }
    fail_add(k); return 0;
}
static void setarc(uint32_t *o, int u, int v) { o[u] |= 1u << v; o[v] &= ~(1u << u); }
int main(int argc, char **argv) {
    n = atoi(argv[1]); int mode = atoi(argv[2]); int count = atoi(argv[3]); rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atol(argv[4]);
    budget = argc > 5 ? atol(argv[5]) : 20000; cands = argc > 6 ? atoi(argv[6]) : 12;
    m = (n - 1) / 2; C = (n * n - 1) / 8;
    dp = malloc(sizeof(uint32_t) * ((size_t)1 << n) * n);
    hk = calloc(HC, sizeof(key128)); hu = calloc(HC, 1);
    int found = 0, undecided = 0, tested = 0;
    int nconn = 0; uint32_t conns[4096];
    if (mode == 1) {   /* all connection sets S with S ∩ -S = empty, |S| = m (circulants) */
        for (uint32_t mask = 0; mask < (1u << m); mask++) { uint32_t S = 0; for (int d = 1; d <= m; d++) S |= ((mask >> (d - 1) & 1) ? (1u << d) : (1u << (n - d))); conns[nconn++] = S; }
        count = nconn;
    }
    for (int it = 0; it < count; it++) {
        uint32_t o[32] = {0};
        if (mode == 0) {
            for (int i = 0; i < n; i++) for (int d = 1; d <= m; d++) setarc(o, i, (i + d) % n);
            for (int r = 0; r < 60 * n * n; r++) { int a = rnd() % n, b = rnd() % n, c2 = rnd() % n; if (a == b || b == c2 || a == c2) continue;
                if ((o[a] >> b & 1) && (o[b] >> c2 & 1) && (o[c2] >> a & 1)) { setarc(o, b, a); setarc(o, c2, b); setarc(o, a, c2); } }
        } else if (mode == 1) {
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && (conns[it] >> ((j - i + n) % n) & 1)) o[i] |= 1u << j;
        } else if (mode == 3) {   /* explicit rows (hex out-masks) from the environment variable ROWS */
            const char *env = getenv("ROWS"); char buf[1024]; strncpy(buf, env, 1023); char *tok = strtok(buf, " ");
            for (int i = 0; i < n && tok; i++, tok = strtok(NULL, " ")) o[i] = (uint32_t)strtoul(tok, NULL, 16);
        } else {   /* Paley: i -> j iff j - i is a nonzero square mod n */
            int qr[64] = {0}; for (int x = 1; x < n; x++) qr[(x * x) % n] = 1;
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && qr[((j - i) % n + n) % n]) o[i] |= 1u << j;
        }
        tested++; nodes = 0; wlen = 0;
        memset(hu, 0, HC);   /* fresh failure memo per instance (it would also be valid to keep it) */
        int res = search(o, 0);
        if (res == 1) {
            uint32_t q[32]; memcpy(q, o, sizeof q); int ok = 1;
            for (int i = 0; i < C; i++) { if (!is_path(q, word[i])) { ok = 0; break; } uint32_t r2[32]; apply_move(q, word[i], r2); memcpy(q, r2, sizeof q); }
            int s[32]; scores(q, s); if (!ok || Lof(s) != 0) { fprintf(stderr, "REPLAY FAILED\n"); return 5; }
            found++;
        } else {
            undecided++;
            printf("UNDECIDED (budget %ld, result %d):", budget, res); for (int i = 0; i < n; i++) printf(" %06x", o[i]); printf("\n");
        }
        fflush(stdout);
    }
    printf("n=%d mode=%d: regular tournaments tested %d; geodesic words of length %d found and replayed %d; undecided %d; search nodes (last) %ld\n", n, mode, tested, C, found, undecided, nodes);
    if (mode == 2 && found) { printf("word:"); for (int i = 0; i < C; i++) { printf(" "); for (int j = 0; j < n; j++) printf("%d%s", word[i][j], j + 1 < n ? "," : ""); printf(";"); } printf("\n"); }
    return 0;
}
