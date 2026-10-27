/* slab_sample_search.c -- mac-mini-2026-10-27-slab.  Test of HYP-9426 (1): every slab class (L + R = (n^2-1)/8) with
 * R <= m-1 is exactly L moves from the transitive class.  Samples slab classes NOT built from TT: a random regular
 * tournament (as in regular_geodesic_search.c mode 0) followed by k random L-descending moves (each lands in the slab
 * with R = k, THM-4819 (A)); then searches a chain of exactly L(S) L-descending moves from S to a transitive tournament
 * (the search of regular_geodesic_search.c), replays and checks it.  UNDECIDED means the node budget ran out.
 * Build: as regular_geodesic_search.c.  Run: ./slab_sample_search n count seed kmax [node_budget] [cands]
 * With ROWS="hex rows" in the environment, the single given slab tournament is searched instead (count 1).
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
typedef struct { uint32_t o[16]; } Tour;
static uint16_t *dp;
static void fill_dp(const uint32_t *o) {
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(uint16_t) * (size_t)(full + 1) * n);
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1 << v;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { uint16_t st = dp[(size_t)S * n + v]; if (!st) continue;
        uint32_t nx = o[v] & ~S; while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[(size_t)(S | 1 << w) * n + w] |= st; } }
}
static int sample_path(const uint32_t *o, int a, int b, int *Q) {
    int full = (1 << n) - 1; if (!(dp[(size_t)full * n + b] >> a & 1)) return 0;
    int S = full, v = b, k = n - 1; Q[k] = b;
    while (k > 0) { int S2 = S & ~(1 << v); int cand[16], nc = 0;
        for (int u = 0; u < n; u++) if ((S2 >> u & 1) && (dp[(size_t)S2 * n + u] >> a & 1) && (o[u] >> v & 1)) cand[nc++] = u;
        if (!nc) return 0; int u = cand[rnd() % nc]; Q[--k] = u; S = S2; v = u; }
    return 1;
}
static void apply_move(const uint32_t *o, const int *Q, uint32_t *r) { uint32_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) r[u] = (~o[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; r[u] |= 1u << v; r[v] &= ~(1u << u); } }
static int is_path(const uint32_t *o, const int *Q) { uint32_t seen = 0; for (int i = 0; i < n; i++) { if (seen >> Q[i] & 1) return 0; seen |= 1u << Q[i]; if (i && !(o[Q[i-1]] >> Q[i] & 1)) return 0; } return seen == (1u << n) - 1; }
static int Lof(const int *s0) { int s[16]; memcpy(s, s0, n * sizeof(int));
    for (int i = 1; i < n; i++) { int v = s[i], j = i - 1; while (j >= 0 && s[j] > v) { s[j + 1] = s[j]; j--; } s[j + 1] = v; }
    int l = 0; for (int i = 0; i < n; i++) l += abs(s[i] - i); return l / 2; }
static void scores(const uint32_t *o, int *s) { for (int v = 0; v < n; v++) s[v] = __builtin_popcount(o[v]); }
/* descending pairs realised by the current dp: list (a,b) */
static int desc_pairs(const uint32_t *o, int *pa, int *pb) {
    int s[16], t[16]; scores(o, s); int L0 = Lof(s), np = 0, full = (1 << n) - 1;
    for (int b = 0; b < n; b++) { uint16_t st = dp[(size_t)full * n + b];
        while (st) { int a = __builtin_ctz(st); st &= st - 1;
            for (int v = 0; v < n; v++) t[v] = n - 1 - s[v]; t[a]++; t[b]--;
            if (Lof(t) == L0 - 1) { pa[np] = a; pb[np] = b; np++; } } }
    return np;
}
typedef struct { uint64_t lo, hi; } key128;
static key128 canon(const uint32_t *q) {
    graph g[MAXN], cg[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_DIGRAPH(options); statsblk stats; options.getcanon = TRUE;
    EMPTYGRAPH(g, 1, n);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (q[i] >> j & 1) ADDONEARC(g, i, j, 1);
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, cg);
    key128 k = {0, 0}; int pos = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (ISELEMENT(GRAPHROW(cg, i, 1), j)) { if (pos < 64) k.lo |= 1ULL << pos; else k.hi |= 1ULL << (pos - 64); } pos++; }
    return k;
}
#define HC (1 << 22)
static key128 *hk; static unsigned char *hu;
static int fail_known(key128 k) { uint64_t h = (k.lo * 0x9E3779B97F4A7C15ULL ^ k.hi * 0xC2B2AE3D27D4EB4FULL) >> 40; h &= HC - 1;
    while (hu[h]) { if (hk[h].lo == k.lo && hk[h].hi == k.hi) return 1; h = (h + 1) & (HC - 1); } return 0; }
static void fail_add(key128 k) { uint64_t h = (k.lo * 0x9E3779B97F4A7C15ULL ^ k.hi * 0xC2B2AE3D27D4EB4FULL) >> 40; h &= HC - 1;
    while (hu[h]) { if (hk[h].lo == k.lo && hk[h].hi == k.hi) return; h = (h + 1) & (HC - 1); } hu[h] = 1; hk[h] = k; }
static long nodes, budget; static int cands;
static int word[64][16]; static int wlen;
static int has_srcsink(const uint32_t *o) { for (int v = 0; v < n; v++) { int s = __builtin_popcount(o[v]); if (s == 0 || s == n - 1) return 1; } return 0; }
static int search(const uint32_t *o, int depth) {
    int s[16]; scores(o, s);
    if (Lof(s) == 0) return 1;
    if (++nodes > budget) return -1;
    key128 k = canon(o); if (fail_known(k)) return 0;
    fill_dp(o);
    int pa[512], pb[512]; int np = desc_pairs(o, pa, pb);
    if (!np) { fail_add(k); return 0; }
    /* candidates: images of sampled paths for random descending pairs, ranked by lookahead */
    Tour cand[64]; int qs[64][16]; int score[64]; int nc = 0;
    for (int t = 0; t < cands && nc < 64; t++) {
        int c = rnd() % np; int Q[16];
        if (!sample_path(o, pa[c], pb[c], Q)) continue;
        apply_move(o, Q, cand[nc].o); memcpy(qs[nc], Q, sizeof Q); nc++;
    }
    for (int i = 0; i < nc; i++) {   /* lookahead score */
        int ss = 0; if (!has_srcsink(cand[i].o)) ss += 1000000;
        int s2[16]; scores(cand[i].o, s2);
        if (Lof(s2) == 0) ss += 100000000;
        else { fill_dp(cand[i].o); int pa2[512], pb2[512]; ss += desc_pairs(cand[i].o, pa2, pb2); }
        score[i] = ss;
    }
    for (int r = 0; r < nc; r++) {
        int best = -1; for (int i = 0; i < nc; i++) if (score[i] >= 0 && (best < 0 || score[i] > score[best])) best = i;
        if (best < 0) break; score[best] = -1;
        memcpy(word[depth], qs[best], sizeof(int) * 16);
        int res = search(cand[best].o, depth + 1);
        if (res == 1) { if (depth + 1 > wlen) wlen = depth + 1; return 1; }
        if (res < 0) return -1;
    }
    fail_add(k); return 0;
}
static void setarc(uint32_t *o, int u, int v) { o[u] |= 1u << v; o[v] &= ~(1u << u); }
int main(int argc, char **argv) {
    n = atoi(argv[1]); int count = atoi(argv[2]); rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atol(argv[3]); int kmax = atoi(argv[4]);
    budget = argc > 5 ? atol(argv[5]) : 20000; cands = argc > 6 ? atoi(argv[6]) : 12;
    m = (n - 1) / 2;
    dp = malloc(sizeof(uint16_t) * ((size_t)1 << n) * n);
    hk = calloc(HC, sizeof(key128)); hu = calloc(HC, 1);
    long found[16] = {0}, undec[16] = {0}, tested[16] = {0}, srcsink[16] = {0};
    for (int it = 0; it < count; it++) {
        uint32_t o[16] = {0};
        const char *env = getenv("ROWS");   /* optional: one explicit slab tournament (hex out-masks), no random walk */
        if (env) { char buf[1024]; strncpy(buf, env, 1023); buf[1023] = 0; char *tok = strtok(buf, " "); for (int i = 0; i < n && tok; i++, tok = strtok(NULL, " ")) o[i] = (uint32_t)strtoul(tok, NULL, 16);
            int sx[16]; scores(o, sx); int Rx = 0; for (int v = 0; v < n; v++) Rx += abs(sx[v] - m); Rx /= 2; kmax = Rx > kmax ? Rx : kmax; goto have_instance; }
        for (int i = 0; i < n; i++) for (int d = 1; d <= m; d++) setarc(o, i, (i + d) % n);
        for (int r = 0; r < 60 * n * n; r++) { int a = rnd() % n, b = rnd() % n, c2 = rnd() % n; if (a == b || b == c2 || a == c2) continue;
            if ((o[a] >> b & 1) && (o[b] >> c2 & 1) && (o[c2] >> a & 1)) { setarc(o, b, a); setarc(o, c2, b); setarc(o, a, c2); } }
        int k = 1 + (int)(rnd() % kmax), ok = 1;
        for (int step = 0; step < k && ok; step++) {   /* one random L-descending move */
            fill_dp(o); int pa[512], pb[512]; int np = desc_pairs(o, pa, pb); if (!np) { ok = 0; break; }
            int Q[16], got = 0; for (int t = 0; t < 50 && !got; t++) { int c = rnd() % np; got = sample_path(o, pa[c], pb[c], Q); }
            if (!got) { ok = 0; break; } uint32_t r2[16]; apply_move(o, Q, r2); memcpy(o, r2, sizeof r2);
        }
        if (!ok) continue;
        have_instance:;
        int s0[16]; scores(o, s0); int L0 = Lof(s0); int R0 = 0; for (int v = 0; v < n; v++) R0 += abs(s0[v] - m); R0 /= 2;
        if (L0 + R0 != (n * n - 1) / 8 || (!env && R0 != k)) { fprintf(stderr, "slab bookkeeping failed\n"); return 6; }
        tested[R0]++; if (has_srcsink(o)) { srcsink[R0]++; continue; }   /* R <= m-1 forbids a source or sink, so this never fires */
        C = L0; nodes = 0; wlen = 0; memset(hu, 0, HC);
        int res = search(o, 0);
        if (res == 1) {
            uint32_t q[16]; memcpy(q, o, sizeof q); int good = 1;
            for (int i = 0; i < C; i++) { if (!is_path(q, word[i])) { good = 0; break; } uint32_t r2[16]; apply_move(q, word[i], r2); memcpy(q, r2, sizeof q); }
            int s[16]; scores(q, s); if (!good || Lof(s) != 0) { fprintf(stderr, "REPLAY FAILED\n"); return 5; }
            found[R0]++;
        } else { undec[R0]++; printf("UNDECIDED R=%d (budget %ld, result %d):", R0, budget, res); for (int i = 0; i < n; i++) printf(" %04x", o[i]); printf("\n"); }
        fflush(stdout);
    }
    for (int R = 1; R <= kmax; R++) printf("n=%d slab samples with R=%d (L=%d): %ld; exactly L moves from TT (word found and replayed): %ld; undecided: %ld; with source/sink: %ld\n",
        n, R, (n * n - 1) / 8 - R, tested[R], found[R], undec[R], srcsink[R]);
    return 0;
}
