/* hratio5.c (klein-2026-11-15-age) -- history-refined feature caps for literal span-3 ribbon words.
   usage: hratio5 mats3.txt frag.txt L K [b0]
   Depth-first search over all words of length <= L (restricted to first letter b0 if given) on the exact 119-state
   integer transfer (sparse).  For every word W with H(W) > 0 and every history length k = 0..K, the key is the last k
   letters of W (words shorter than k are skipped for that k); for each key the maximum over words of
     f_j(x_W) / H(W)
   is recorded for the features j = 0..N-1 (the coordinates x_s), N..N+2 (sums of x_s over states with 2, 3, 4
   fragments), N+3..N+10 (Z_l = sum over two-fragment states of x_s (h A_l)_s, l = 0..7), N+11..N+18 (the full next-step counts
   H(W l) = (h A_l) . x, l = 0..7, so that cap/1 is the per-letter growth ratio), as exact fractions.
   Output lines: "k key feature num den" for every recorded (k, key, feature) with a positive maximum. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef __int128 i128;
static int N, L, K, b0 = -1, NF;
static ll *Ad[8]; static ll *v0, *h; static int *nfrag;
static int nnz[8]; static int *ri[8], *ci[8]; static ll *val[8];
static ll hA[8][256];
static ll *F;                 /* NF x N feature vectors */
static ll *bnum, *bden;       /* [k][key][feature] */
static int word[64];
static size_t idx(int k, int key, int j) { size_t base = 0; for (int t = 0; t < k; t++) base += (size_t)1 << (3 * t); return ((base + key) * NF) + j; }
static ll *fv;                /* feature values of the current node */
static void record(int depth, const ll *x) {
    ll Hx = 0; for (int j = 0; j < N; j++) Hx += h[j] * x[j];
    if (Hx <= 0) return;
    for (int j = 0; j < N; j++) fv[j] = x[j];                       /* coordinate features */
    for (int j = N; j < NF; j++) { ll v = 0; const ll *f = F + (size_t)j * N; for (int s = 0; s < N; s++) v += f[s] * x[s]; fv[j] = v; }
    for (int k = 0; k <= K && k <= depth; k++) {
        int key = 0; for (int t = depth - k; t < depth; t++) key = key * 8 + word[t];
        size_t base = idx(k, key, 0);
        for (int j = 0; j < NF; j++) {
            ll v = fv[j]; if (v <= 0) continue;
            size_t id = base + j;
            if (bden[id] == 0 || (i128)v * bden[id] > (i128)bnum[id] * Hx) { bnum[id] = v; bden[id] = Hx; }
        }
    }
}
static void dfs(int depth, const ll *x) {
    record(depth, x);
    if (depth == L) return;
    ll *y = malloc(sizeof(ll) * N);
    int lo = 0, hi = 7; if (depth == 0 && b0 >= 0) lo = hi = b0;
    for (int b = lo; b <= hi; b++) {
        memset(y, 0, sizeof(ll) * N);
        for (int t = 0; t < nnz[b]; t++) y[ri[b][t]] += val[b][t] * x[ci[b][t]];
        word[depth] = b; dfs(depth + 1, y);
    }
    free(y);
}
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); if (!f) return 1; if (fscanf(f, "%d", &N) != 1) return 1; L = atoi(argv[3]); K = atoi(argv[4]); if (argc > 5) b0 = atoi(argv[5]); if (b0 > 7) b0 = -1;
    for (int b = 0; b < 8; b++) { Ad[b] = malloc(sizeof(ll) * N * N); for (int i = 0; i < N * N; i++) if (fscanf(f, "%lld", &Ad[b][i]) != 1) return 2; }
    v0 = malloc(sizeof(ll) * N); h = malloc(sizeof(ll) * N);
    for (int i = 0; i < N; i++) if (fscanf(f, "%lld", &v0[i]) != 1) return 3;
    for (int i = 0; i < N; i++) if (fscanf(f, "%lld", &h[i]) != 1) return 4;
    fclose(f);
    FILE *g = fopen(argv[2], "r"); if (!g) return 5; nfrag = malloc(sizeof(int) * N); for (int i = 0; i < N; i++) if (fscanf(g, "%d", &nfrag[i]) != 1) return 6; fclose(g);
    for (int b = 0; b < 8; b++) {
        int m = 0; for (int i = 0; i < N * N; i++) if (Ad[b][i]) m++;
        nnz[b] = m; ri[b] = malloc(sizeof(int) * m); ci[b] = malloc(sizeof(int) * m); val[b] = malloc(sizeof(ll) * m);
        int t = 0; for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (Ad[b][(size_t)i * N + j]) { ri[b][t] = i; ci[b][t] = j; val[b][t] = Ad[b][(size_t)i * N + j]; t++; }
        for (int j = 0; j < N; j++) { ll s = 0; for (int i = 0; i < N; i++) s += h[i] * Ad[b][(size_t)i * N + j]; hA[b][j] = s; }
    }
    NF = N + 3 + 8 + 8; F = calloc((size_t)NF * N, sizeof(ll));
    for (int s = 0; s < N; s++) F[(size_t)s * N + s] = 1;
    for (int s = 0; s < N; s++) for (int c = 2; c <= 4; c++) if (nfrag[s] == c) F[(size_t)(N + c - 2) * N + s] = 1;
    for (int b = 0; b < 8; b++) for (int s = 0; s < N; s++) if (nfrag[s] == 2) F[(size_t)(N + 3 + b) * N + s] = hA[b][s];
    for (int b = 0; b < 8; b++) for (int s = 0; s < N; s++) F[(size_t)(N + 11 + b) * N + s] = hA[b][s];
    size_t total = 0; for (int k = 0; k <= K; k++) total += (size_t)1 << (3 * k); total *= NF;
    bnum = calloc(total, sizeof(ll)); bden = calloc(total, sizeof(ll)); fv = malloc(sizeof(ll) * NF);
    ll *x = malloc(sizeof(ll) * N); memcpy(x, v0, sizeof(ll) * N);
    dfs(0, x);
    printf("# hratio5 N=%d L=%d K=%d b0=%d features=%d (0..%d coordinates, %d..%d fragment-class sums 2,3,4, %d..%d Z_l for l=0..7, %d..%d H(W l) for l=0..7)\n", N, L, K, b0, NF, N - 1, N, N + 2, N + 3, N + 10, N + 11, N + 18);
    for (int k = 0; k <= K; k++) for (int key = 0; key < (1 << (3 * k)); key++) for (int j = 0; j < NF; j++) {
        size_t id = idx(k, key, j); if (bden[id]) printf("%d %d %d %lld %lld\n", k, key, j, bnum[id], bden[id]);
    }
    return 0;
}
