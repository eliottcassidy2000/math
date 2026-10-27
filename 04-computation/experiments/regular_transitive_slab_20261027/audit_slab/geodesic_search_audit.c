/* geodesic_search_audit.c -- the auditor's own search for THM-4819 (C): a word of exactly c = (n^2-1)/8 native moves
 * from a regular tournament to a transitive one (L drops by one per move; the lower bound d >= L = c is THM-4810).
 * Written independently of regular_geodesic_search.c / paley_geodesic_search.c (no nauty; own endpoint DP, own
 * canonical form for the failure memo, own ranking: images without a source or sink first (THM-4819 (A4)), then by the
 * number of L-descending endpoint pairs they realise).  Prints every word found as  WORD <rows hex> : <path>;<path>;...
 * so that replay_words.py can re-check it.  Input lines: n hex out-masks (one tournament per line).
 * Usage: ./geodesic_search_audit n budget cands seed < instances.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, c; static uint32_t FULL;
static uint64_t rs = 0x2545F4914F6CDD1Dull;
static uint32_t rnd(uint32_t k) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (uint32_t)(rs % k); }
typedef struct { uint16_t o[16]; } Tour;

static int Lof(const int *s0) { int s[16]; memcpy(s, s0, sizeof(int) * n); for (int i = 1; i < n; i++) for (int j = i; j > 0 && s[j] < s[j-1]; j--) { int t = s[j]; s[j] = s[j-1]; s[j-1] = t; } int l = 0; for (int i = 0; i < n; i++) l += abs(s[i] - i); return l / 2; }
static void scores(const Tour *T, int *s) { for (int v = 0; v < n; v++) s[v] = __builtin_popcount(T->o[v]); }
static Tour native(const Tour *T, const int *Q) {
    Tour R; for (int u = 0; u < n; u++) R.o[u] = (uint16_t)(~T->o[u] & FULL & ~(1u << u));
    for (int i = 0; i + 1 < n; i++) { int a = Q[i], b = Q[i + 1]; R.o[a] |= (uint16_t)(1u << b); R.o[b] &= (uint16_t)~(1u << a); }
    return R;
}
/* endpoint DP: dp[S*n + v] = set of a with a Hamilton path of T[S] from a to v */
static uint16_t *dp;
static void fill(const Tour *T) {
    for (uint32_t S = 1; S <= FULL; S++) {
        for (int v = 0; v < n; v++) {
            if (!(S >> v & 1)) { dp[S * n + v] = 0; continue; }
            uint32_t S2 = S & ~(1u << v);
            if (!S2) { dp[S * n + v] = (uint16_t)(1u << v); continue; }
            uint16_t acc = 0;
            for (uint32_t x = S2; x; x &= x - 1) { int u = __builtin_ctz(x); if (T->o[u] >> v & 1) acc |= dp[S2 * n + u]; }
            dp[S * n + v] = acc;
        }
    }
}
static int sample(const Tour *T, int a, int b, int *Q) {
    if (!(dp[FULL * n + b] >> a & 1)) return 0;
    uint32_t S = FULL; int cur = b; Q[n - 1] = b;
    for (int k = n - 2; k >= 0; k--) {
        uint32_t S2 = S & ~(1u << cur); int cand[16], nc = 0;
        for (uint32_t x = S2; x; x &= x - 1) { int u = __builtin_ctz(x); if ((T->o[u] >> cur & 1) && (dp[S2 * n + u] >> a & 1)) cand[nc++] = u; }
        if (!nc) return 0;
        cur = cand[rnd(nc)]; Q[k] = cur; S = S2;
    }
    return Q[0] == a;
}
static int desc_pairs(const Tour *T, int *pa, int *pb) {   /* uses the current dp of T */
    int s[16]; scores(T, s); int L0 = Lof(s), np = 0;
    for (int b = 0; b < n; b++) for (uint32_t x = dp[FULL * n + b]; x; x &= x - 1) {
        int a = __builtin_ctz(x); if (a == b) continue;
        int t[16]; for (int v = 0; v < n; v++) t[v] = n - 1 - s[v]; t[a]++; t[b]--;
        if (Lof(t) == L0 - 1) { pa[np] = a; pb[np] = b; np++; }
    }
    return np;
}
/* canonical form (individualisation-refinement), 128-bit code */
typedef struct { uint64_t w[2]; } Key;
static Key bestk; static int bestset;
static int keyless(Key a, Key b) { return a.w[1] != b.w[1] ? a.w[1] < b.w[1] : a.w[0] < b.w[0]; }
static void refine(const Tour *T, int *col) {
    int k_old = -1;
    for (;;) {
        int k = 0; for (int v = 0; v < n; v++) if (col[v] + 1 > k) k = col[v] + 1;
        if (k == k_old) return; k_old = k;
        int sig[16][17];
        for (int v = 0; v < n; v++) { sig[v][0] = col[v]; for (int j = 0; j < k; j++) sig[v][1 + j] = 0; for (uint32_t x = T->o[v]; x; x &= x - 1) sig[v][1 + col[__builtin_ctz(x)]]++; }
        int order[16]; for (int v = 0; v < n; v++) order[v] = v;
        for (int i = 1; i < n; i++) for (int j = i; j > 0; j--) { int a = order[j-1], b = order[j], cmp = 0; for (int t = 0; t <= k && !cmp; t++) cmp = (sig[a][t] > sig[b][t]) - (sig[a][t] < sig[b][t]); if (cmp > 0) { order[j-1] = b; order[j] = a; } else break; }
        int nc[16], r = 0;
        for (int i = 0; i < n; i++) { if (i > 0) { int a = order[i-1], b = order[i], same = 1; for (int t = 0; t <= k; t++) if (sig[a][t] != sig[b][t]) { same = 0; break; } if (!same) r++; } nc[order[i]] = r; }
        memcpy(col, nc, sizeof(int) * n);
    }
}
static void csearch(const Tour *T, int *col) {
    refine(T, col);
    int k = 0; for (int v = 0; v < n; v++) if (col[v] + 1 > k) k = col[v] + 1;
    if (k == n) {
        int at[16]; for (int v = 0; v < n; v++) at[col[v]] = v;
        Key K = {{0, 0}}; int b = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (T->o[at[i]] >> at[j] & 1) K.w[b >> 6] |= 1ull << (b & 63); b++; }
        if (!bestset || keyless(K, bestk)) { bestk = K; bestset = 1; }
        return;
    }
    int cnt[16] = {0}; for (int v = 0; v < n; v++) cnt[col[v]]++;
    int target = 0; while (cnt[target] < 2) target++;
    for (int w = 0; w < n; w++) {
        if (col[w] != target) continue;
        int c2[16]; for (int v = 0; v < n; v++) c2[v] = 2 * col[v] + 1; c2[w] = 2 * col[w];
        int used[40] = {0}; for (int v = 0; v < n; v++) used[c2[v]] = 1; int rk[40], r = 0; for (int x = 0; x < 40; x++) { rk[x] = r; if (used[x]) r++; }
        for (int v = 0; v < n; v++) c2[v] = rk[c2[v]];
        csearch(T, c2);
    }
}
static Key canon(const Tour *T) {
    int col[16]; for (int v = 0; v < n; v++) col[v] = __builtin_popcount(T->o[v]);
    int used[17] = {0}; for (int v = 0; v < n; v++) used[col[v]] = 1; int rk[17], r = 0; for (int x = 0; x < 17; x++) { rk[x] = r; if (used[x]) r++; }
    for (int v = 0; v < n; v++) col[v] = rk[col[v]];
    bestset = 0; csearch(T, col); return bestk;
}
#define HC (1 << 20)
static Key *hk; static unsigned char *hu;
static uint32_t kh(Key k) { return (uint32_t)((k.w[0] * 0x9E3779B97F4A7C15ull ^ k.w[1] * 0xC2B2AE3D27D4EB4Full) >> 44); }
static int known(Key k) { uint32_t h = kh(k) & (HC - 1); while (hu[h]) { if (hk[h].w[0] == k.w[0] && hk[h].w[1] == k.w[1]) return 1; h = (h + 1) & (HC - 1); } return 0; }
static void addk(Key k) { uint32_t h = kh(k) & (HC - 1); while (hu[h]) { if (hk[h].w[0] == k.w[0] && hk[h].w[1] == k.w[1]) return; h = (h + 1) & (HC - 1); } hu[h] = 1; hk[h] = k; }
static long nodes, budget; static int cands;
static int word[64][16];
static int srcsink(const Tour *T) { for (int v = 0; v < n; v++) { int s = __builtin_popcount(T->o[v]); if (s == 0 || s == n - 1) return 1; } return 0; }
static int search(const Tour *T, int depth) {
    int s[16]; scores(T, s); if (Lof(s) == 0) return 1;
    if (++nodes > budget) return -1;
    Key K = canon(T); if (known(K)) return 0;
    fill(T);
    static int pa[64][256], pb[64][256];
    int np = desc_pairs(T, pa[depth], pb[depth]);
    if (!np) { addk(K); return 0; }
    Tour im[32]; int qs[32][16], rank[32], ni = 0;
    for (int t = 0; t < cands && ni < 32; t++) { int k = rnd(np); int Q[16]; if (!sample(T, pa[depth][k], pb[depth][k], Q)) continue; im[ni] = native(T, Q); memcpy(qs[ni], Q, sizeof Q); ni++; }
    for (int i = 0; i < ni; i++) {
        int s2[16]; scores(&im[i], s2);
        if (Lof(s2) == 0) { rank[i] = 1 << 30; continue; }
        rank[i] = srcsink(&im[i]) ? 0 : (1 << 20);
        static int qa[256], qb[256]; fill(&im[i]); rank[i] += desc_pairs(&im[i], qa, qb);
    }
    for (int r = 0; r < ni; r++) {
        int best = -1; for (int i = 0; i < ni; i++) if (rank[i] >= 0 && (best < 0 || rank[i] > rank[best])) best = i;
        if (best < 0) break; rank[best] = -1;
        memcpy(word[depth], qs[best], sizeof(int) * 16);
        int res = search(&im[best], depth + 1);
        if (res != 0) return res;
    }
    addk(K); return 0;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); budget = atol(argv[2]); cands = atoi(argv[3]); rs ^= (uint64_t)atoll(argv[4]) * 0x9E3779B97F4A7C15ull;
    c = (n * n - 1) / 8; FULL = (1u << n) - 1;
    dp = malloc(sizeof(uint16_t) * ((size_t)FULL + 1) * n); hk = calloc(HC, sizeof(Key)); hu = calloc(HC, 1);
    char line[4096]; int inst = 0, found = 0, und = 0;
    while (fgets(line, sizeof line, stdin)) {
        Tour T; memset(&T, 0, sizeof T); char *tok = strtok(line, " \n"); int k = 0;
        while (tok && k < n) { T.o[k++] = (uint16_t)strtoul(tok, NULL, 16); tok = strtok(NULL, " \n"); }
        if (k < n) continue;
        int s[16]; scores(&T, s); for (int v = 0; v < n; v++) if (s[v] != (n - 1) / 2) { fprintf(stderr, "not regular\n"); return 2; }
        inst++; nodes = 0; memset(hu, 0, HC);
        int res = search(&T, 0);
        if (res == 1) {
            found++; printf("WORD"); for (int v = 0; v < n; v++) printf(" %x", T.o[v]); printf(" :");
            for (int j = 0; j < c; j++) { printf(" "); for (int i = 0; i < n; i++) printf("%d%s", word[j][i], i + 1 < n ? "," : ""); printf(";"); }
            printf("  [nodes %ld]\n", nodes);
        } else { und++; printf("UNDECIDED"); for (int v = 0; v < n; v++) printf(" %x", T.o[v]); printf("  [result %d]\n", res); }
        fflush(stdout);
    }
    printf("n=%d: instances %d; words found %d; undecided %d\n", n, inst, found, und);
    return 0;
}
