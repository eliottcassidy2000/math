/* regular_L_geodesics.c -- mac-mini-2026-10-27-dipole
 * For each regular tournament read from gentourng (e.g. gentourng -d5D5 11), search for a word of exactly
 * L = (n^2-1)/8 native moves to a transitive tournament, i.e. a chain along which L = (1/2) sum |s_(i) - i| drops by
 * one per move (THM-4815: a move lowers L iff its path joins a descending endpoint pair).  Randomised greedy descent
 * with restarts: at each state compute the endpoint relation Pi(T) by a subset DP, list the descending pairs it
 * realises, pick one at random, sample a random Hamiltonian path between them, apply the move.  A found word proves
 * d(R, transitive class) = L(R) = (n^2-1)/8 for that class (the lower bound is THM-4810 Theorem L); every word found is
 * replayed and checked.  Build: cc -O3 regular_L_geodesics.c -o regular_L_geodesics
 * Run: gentourng -q -d5D5 11 | ./regular_L_geodesics 11 restarts seed */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n; static uint32_t out_[16];
static uint16_t *dp;    /* dp[S*n+v] = starts a of paths covering S ending at v (bitmask) */
static uint64_t rs;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static void fill_dp(const uint32_t *o) {
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(uint16_t) * (size_t)(full + 1) * n);
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1 << v;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { uint16_t st = dp[(size_t)S * n + v]; if (!st) continue;
        uint32_t nx = o[v] & ~S; while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[(size_t)(S | 1 << w) * n + w] |= st; } }
}
static int Lof(const int *s) { int t[16]; memcpy(t, s, n * sizeof(int));
    for (int i = 1; i < n; i++) { int v = t[i], j = i - 1; while (j >= 0 && t[j] > v) { t[j + 1] = t[j]; j--; } t[j + 1] = v; }
    int l = 0; for (int i = 0; i < n; i++) l += abs(t[i] - i); return l / 2; }
/* sample a uniformly-random-ish Hamiltonian path from a to b (backtrack through dp, which records starts) */
static int sample_path(const uint32_t *o, int a, int b, int *Q) {
    int full = (1 << n) - 1;
    if (!(dp[(size_t)full * n + b] >> a & 1)) return 0;
    int S = full, v = b, k = n - 1; Q[k] = b;
    while (k > 0) {
        int S2 = S & ~(1 << v); int cand[16], nc = 0;
        for (int u = 0; u < n; u++) if ((S2 >> u & 1) && (dp[(size_t)S2 * n + u] >> a & 1) && (o[u] >> v & 1)) cand[nc++] = u;
        if (!nc) return 0;
        int u = cand[rnd() % nc]; Q[--k] = u; S = S2; v = u;
    }
    return Q[0] == a;
}
static void apply_move(uint32_t *o, const int *Q) {
    uint32_t full = (1u << n) - 1, r[16];
    for (int u = 0; u < n; u++) r[u] = (~o[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; r[u] |= 1u << v; r[v] &= ~(1u << u); }
    memcpy(o, r, n * sizeof(uint32_t));
}
static int is_path(const uint32_t *o, const int *Q) { uint32_t seen = 0; for (int i = 0; i < n; i++) { if (seen >> Q[i] & 1) return 0; seen |= 1u << Q[i]; if (i && !(o[Q[i-1]] >> Q[i] & 1)) return 0; } return seen == (1u << n) - 1; }
int main(int argc, char **argv) {
    n = atoi(argv[1]); int restarts = atoi(argv[2]); rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atol(argv[3]);
    dp = malloc(sizeof(uint16_t) * ((size_t)1 << n) * n);
    char line[256]; long classes = 0, found = 0; long tries_hist[64] = {0};
    int target = (n * n - 1) / 8;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        classes++;
        uint32_t T0[16] = {0}; int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T0[i] |= 1u << j; else T0[j] |= 1u << i; k++; }
        int ok = 0, tries = 0;
        for (tries = 1; tries <= restarts && !ok; tries++) {
            uint32_t o[16]; memcpy(o, T0, sizeof o);
            int word[64][16]; int len = 0;
            for (;;) {
                int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(o[v]);
                int L0 = Lof(s); if (L0 == 0) { ok = 1; break; }
                fill_dp(o);
                int full = (1 << n) - 1; int pa[256], pb[256], np = 0, t[16];
                for (int b = 0; b < n; b++) { uint16_t st = dp[(size_t)full * n + b];
                    while (st) { int a = __builtin_ctz(st); st &= st - 1;
                        for (int v = 0; v < n; v++) t[v] = n - 1 - s[v]; t[a]++; t[b]--;
                        if (Lof(t) == L0 - 1 && np < 256) { pa[np] = a; pb[np] = b; np++; } } }
                if (!np) break;   /* local minimum: restart */
                int c = rnd() % np; int Q[16];
                if (!sample_path(o, pa[c], pb[c], Q)) { fprintf(stderr, "sampling failed\n"); return 2; }
                memcpy(word[len++], Q, sizeof Q); apply_move(o, Q);
            }
            if (ok) {   /* replay and verify */
                if (len != target) { fprintf(stderr, "length %d != %d\n", len, target); return 3; }
                uint32_t o2[16]; memcpy(o2, T0, sizeof o2);
                for (int i = 0; i < len; i++) { if (!is_path(o2, word[i])) { fprintf(stderr, "replay failed\n"); return 4; } apply_move(o2, word[i]); }
                int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(o2[v]); if (Lof(s) != 0) { fprintf(stderr, "not transitive\n"); return 5; }
            }
        }
        if (ok) { found++; tries_hist[tries - 1 < 63 ? tries - 1 : 63]++; }
        else printf("NO WORD FOUND within %d restarts: %s", restarts, line);
    }
    printf("n=%d regular classes %ld; L-geodesic words of length %d found (and replayed) for %ld; restarts needed:", n, classes, target, found);
    for (int i = 1; i < 64; i++) if (tries_hist[i]) printf(" %d:%ld", i, tries_hist[i]); printf("\n");
    return 0;
}
