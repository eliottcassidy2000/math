/* hratio.c (klein-2026-11-15) -- the one-vertex growth ratio of literal span-3 ribbon words.
   usage: hratio mats3.txt L
   Reads the 119-state profile transfer (eight letter matrices A_0..A_7, the initial vector v0 = A_0 e0 (the first
   vertex) and the final functional h) and enumerates every word W of length 1..L over the eight letters by depth-first
   search on state vectors x_W = A_W v0 (64-bit; H < 2^63 for L <= 25).  For each length it reports
     max over W of length n and letters l of H(W l) / H(W),
   the maximising word and letter, and the number of words.  Also the maximum over all lengths <= L. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
static int N; static ll *A[8]; static ll *v0, *h;
static int L; static double best[64]; static int bestword[64][64], bestl[64]; static ll bestnum[64], bestden[64];
static ll cnt[64]; static int *nfrag; static ll *hA[8]; static double envA[64]; static ll bestZ[64], bestZH[64]; static double bestZr[64];
#define EMAX 70000
static ll env[EMAX];
static void apply(int b, const ll *x, ll *y) {  /* y = A_b x */
    const ll *M = A[b];
    for (int i = 0; i < N; i++) { ll s = 0; const ll *row = M + (size_t)i * N; for (int j = 0; j < N; j++) s += row[j] * x[j]; y[i] = s; }
}
static ll dot(const ll *a, const ll *b) { ll s = 0; for (int j = 0; j < N; j++) s += a[j] * b[j]; return s; }
static int word[64];
static void dfs(int depth, const ll *x) {
    /* x = state of the word of length depth (>= 1) */
    ll Hx = dot(h, x); cnt[depth]++;
    ll *y = malloc(sizeof(ll) * N);
    for (int b = 0; b < 8; b++) {
        apply(b, x, y); ll Hy = dot(h, y);
        if (Hx > 0) {
            ll two = 0; for (int j = 0; j < N; j++) if (nfrag[j] == 2) two += hA[b][j] * x[j];
            double zr = (double)two / (double)Hx; if (zr > bestZr[depth]) { bestZr[depth] = zr; bestZ[depth] = two; bestZH[depth] = Hx; }
            if (Hx < EMAX && Hy > env[Hx]) env[Hx] = Hy;
            double r = (double)Hy / (double)Hx;
            if (r > best[depth] + 1e-12 || (r > best[depth] - 1e-12 && (bestden[depth] == 0))) {
                best[depth] = r; bestnum[depth] = Hy; bestden[depth] = Hx; bestl[depth] = b; memcpy(bestword[depth], word, sizeof(int) * depth);
            }
        }
        if (depth < L) { word[depth] = b; dfs(depth + 1, y); }
    }
    free(y);
}
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); L = atoi(argv[2]); if (fscanf(f, "%d", &N) != 1) return 1;
    for (int b = 0; b < 8; b++) { A[b] = malloc(sizeof(ll) * N * N); for (int i = 0; i < N * N; i++) if (fscanf(f, "%lld", &A[b][i]) != 1) return 2; }
    v0 = malloc(sizeof(ll) * N); h = malloc(sizeof(ll) * N);
    for (int i = 0; i < N; i++) if (fscanf(f, "%lld", &v0[i]) != 1) return 3;
    for (int i = 0; i < N; i++) if (fscanf(f, "%lld", &h[i]) != 1) return 4;
    fclose(f);
    nfrag = malloc(sizeof(int) * N); FILE *g = fopen("frag.txt", "r"); for (int i = 0; i < N; i++) if (fscanf(g, "%d", &nfrag[i]) != 1) return 5; fclose(g);
    for (int b = 0; b < 8; b++) { hA[b] = malloc(sizeof(ll) * N); for (int j = 0; j < N; j++) { ll s = 0; for (int i = 0; i < N; i++) s += h[i] * A[b][(size_t)i * N + j]; hA[b][j] = s; } }
    /* words of length n: the first vertex is v0 (length 0 state); a word of length n applies n letters */
    ll *x = malloc(sizeof(ll) * N); memcpy(x, v0, sizeof(ll) * N);
    /* depth 0 = empty word (one vertex): ratio for the first letter */
    dfs(0, x);
    double overall = 0;
    for (int n = 0; n <= L; n++) {
        printf("length %d: words %lld, max H(Wl)/H(W) = %.6f = %lld/%lld, letter %d, word", n, cnt[n], best[n], bestnum[n], bestden[n], bestl[n]);
        for (int i = 0; i < n; i++) printf(" %d", bestword[n][i]); printf("\n");
        if (best[n] > overall) overall = best[n];
    }
    printf("overall max ratio over lengths 0..%d: %.6f\n", L, overall);
    for (int n = 0; n <= L; n++) printf("length %d: max two-fragment part / H(W) = %.4f = %lld/%lld\n", n, bestZr[n], bestZ[n], bestZH[n]);
    printf("envelope env[H] = max H(Wl) over words with H(W) = H (all lengths <= L), with env[H]/H and env[H] - 5H:\n");
    for (ll H = 1; H < EMAX; H++) if (env[H]) { static int shown = 0; if (shown < 60 || env[H] > 5 * H + 30) { printf("  H=%lld env=%lld ratio=%.3f env-5H=%lld\n", H, env[H], (double)env[H] / H, env[H] - 5 * H); shown++; } }
    return 0;
}
