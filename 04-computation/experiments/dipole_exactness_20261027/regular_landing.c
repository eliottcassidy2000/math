/* regular_landing.c -- mac-mini-2026-10-27-dipole
 * Which regular classes are exactly (n^2-1)/8 native moves from the transitive class?  By THM-4815 B8 every R-descending
 * chain from the transitive tournament TT_n has exactly R(TT_n) = (n^2-1)/8 moves and ends at a regular tournament, and
 * moves are involutions, so every regular class reached this way is at class distance exactly (n^2-1)/8 = L from the
 * transitive class (lower bound THM-4810 Theorem L).  Random R-descending chains (each move along a random Hamiltonian
 * path joining a random pair (a, b) with s(a) > (n-1)/2 > s(b)); the landing classes are canonicalised with nauty and
 * counted.  Build: cc -O3 -I/opt/homebrew/opt/nauty/include regular_landing.c /opt/homebrew/opt/nauty/lib/libnauty.a -o regular_landing
 * Run: ./regular_landing n trials seed [listfile] */
#define MAXN 32
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n; static uint16_t *dp; static uint64_t rs;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
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
static int is_ham_path(const uint32_t *o, const int *Q) { uint32_t seen = 0; for (int i = 0; i < n; i++) { if (seen >> Q[i] & 1) return 0; seen |= 1u << Q[i]; if (i && !(o[Q[i-1]] >> Q[i] & 1)) return 0; } return seen == (1u << n) - 1; }
static void apply_move(uint32_t *o, const int *Q) { uint32_t full = (1u << n) - 1, r[16];
    for (int u = 0; u < n; u++) r[u] = (~o[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; r[u] |= 1u << v; r[v] &= ~(1u << u); }
    memcpy(o, r, n * sizeof(uint32_t)); }
#define HC (1 << 20)
static key128 hk[HC]; static unsigned char hu[HC]; static long nseen = 0;
static int ins(key128 k) { uint64_t h = (k.lo * 0x9E3779B97F4A7C15ULL ^ k.hi * 0xC2B2AE3D27D4EB4FULL) >> 40; h &= HC - 1;
    while (hu[h]) { if (hk[h].lo == k.lo && hk[h].hi == k.hi) return 0; h = (h + 1) & (HC - 1); } hu[h] = 1; hk[h] = k; nseen++; return 1; }
int main(int argc, char **argv) {
    n = atoi(argv[1]); long trials = atol(argv[2]); rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atol(argv[3]);
    FILE *fl = argc > 4 ? fopen(argv[4], "w") : NULL;
    dp = malloc(sizeof(uint16_t) * ((size_t)1 << n) * n);
    int m = (n - 1) / 2, target = (n * n - 1) / 8;
    long milestones[] = {10, 100, 1000, 10000, 100000, 1000000};
    for (long t = 1; t <= trials; t++) {
        uint32_t o[16]; for (int i = 0; i < n; i++) { o[i] = 0; for (int j = i + 1; j < n; j++) o[i] |= 1u << j; }
        int steps = 0;
        for (;;) {
            int s[16], R = 0; for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(o[v]); R += abs(s[v] - m); }
            if (R == 0) break;
            fill_dp(o); int full = (1 << n) - 1; int pa[256], pb[256], np = 0;
            for (int b = 0; b < n; b++) if (s[b] < m) { uint16_t st = dp[(size_t)full * n + b];
                while (st) { int a = __builtin_ctz(st); st &= st - 1; if (s[a] > m && np < 256) { pa[np] = a; pb[np] = b; np++; } } }
            if (!np) { fprintf(stderr, "B7/B8 violated\n"); return 2; }
            int c = rnd() % np, Q[16]; if (!sample_path(o, pa[c], pb[c], Q)) { fprintf(stderr, "sample failed\n"); return 3; }
            if (!is_ham_path(o, Q)) { fprintf(stderr, "not a Hamiltonian path\n"); return 5; }
            apply_move(o, Q); steps++;
        }
        if (steps != target) { fprintf(stderr, "steps %d\n", steps); return 4; }
        { int s2[16]; for (int v = 0; v < n; v++) { s2[v] = __builtin_popcount(o[v]); if (s2[v] != m) { fprintf(stderr, "not regular\n"); return 6; } } }
        key128 k = canon(o);
        if (ins(k) && fl) { fprintf(fl, "%016llx%016llx\n", (unsigned long long)k.hi, (unsigned long long)k.lo); }
        for (int q = 0; q < 6; q++) if (t == milestones[q]) printf("  after %ld chains: %ld distinct regular classes\n", t, nseen);
    }
    printf("n=%d: %ld random R-descending chains from the transitive tournament (each %d moves) land on %ld distinct regular classes\n", n, trials, target, nseen);
    if (fl) fclose(fl);
    return 0;
}
