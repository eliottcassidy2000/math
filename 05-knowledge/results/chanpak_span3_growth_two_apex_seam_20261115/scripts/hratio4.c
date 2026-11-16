/* hratio4.c (klein-2026-11-15) -- exact one-vertex growth maxima of literal span-3 ribbon words, with the two-fragment
   part Z and exact tie counts.  usage: hratio4 mats3.txt frag.txt L [b0] [print]   (b0 = -1 for all first letters; print lists every maximising pair)
   Reads the 119-state profile transfer (eight letter matrices A_0..A_7, v0, h) as hratio.c does, plus frag.txt (the
   number of fragments of each state).  Enumerates every word W of length 0..L (all 8^n words; with b0 given, only the
   words whose first letter is b0, so that eight processes partition the work) by depth-first search on the exact integer
   state vectors x_W = A_W v0, using the sparse form of the letter matrices.  For each length n it reports, over all words
   W of length n and letters l,
     max H(W l)/H(W)  and  max Z_l(W)/H(W),  Z_l(W) = sum over two-fragment states s of x_s (h A_l)_s
   as exact fractions (128-bit cross-multiplication), the number of pairs (W, l) attaining each maximum, and the first
   maximiser found.  H(W l) = e + s + Z_l (THM-5243 (1)); Z_l is the bridged two-path-cover count. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef __int128 i128;
static int N, L, b0 = -1;
static ll *Ad[8]; static ll *v0, *h; static int *nfrag;
static int nnz[8]; static int *ri[8], *ci[8]; static ll *val[8];   /* sparse triples per letter */
static ll hA[8][256], twoA[8][256];
typedef struct { ll num, den; ll ties; int letter; int word[64]; } rec;
static rec bestH[64], bestZ[64]; static ll cnt[64];
static int word[64]; static int printing = 0;
static void show(const char *what, int depth, int b, ll num, ll den) { printf("  %s length %d: letter %d, word", what, depth, b); for (int i = 0; i < depth; i++) printf(" %d", word[i]); printf("  (%lld/%lld)\n", num, den); }
static void consider(rec *r, ll num, ll den, int depth, int b) {
    if (den <= 0) return;
    if (r->den == 0) { r->num = num; r->den = den; r->ties = 1; r->letter = b; memcpy(r->word, word, sizeof(int) * depth); return; }
    i128 lhs = (i128)num * r->den, rhs = (i128)r->num * den;
    if (lhs > rhs) { r->num = num; r->den = den; r->ties = 1; r->letter = b; memcpy(r->word, word, sizeof(int) * depth); }
    else if (lhs == rhs) r->ties++;
}
static void dfs(int depth, const ll *x) {
    ll Hx = 0; for (int j = 0; j < N; j++) Hx += h[j] * x[j];
    if (!printing) cnt[depth]++;
    ll *y = malloc(sizeof(ll) * N);
    int lo = 0, hi = 7; if (depth == 0 && b0 >= 0) { lo = hi = b0; }
    for (int b = lo; b <= hi; b++) {
        memset(y, 0, sizeof(ll) * N);
        const int *R = ri[b], *C = ci[b]; const ll *V = val[b]; int m = nnz[b];
        for (int k = 0; k < m; k++) y[R[k]] += V[k] * x[C[k]];
        ll Hy = 0, Z = 0; for (int j = 0; j < N; j++) { Hy += h[j] * y[j]; Z += twoA[b][j] * x[j]; }
        if (Hx > 0 && !printing) { consider(&bestH[depth], Hy, Hx, depth, b); consider(&bestZ[depth], Z, Hx, depth, b); }
        if (Hx > 0 && printing) { if ((i128)Hy * bestH[depth].den == (i128)bestH[depth].num * Hx) show("maxH", depth, b, Hy, Hx);
                                  if ((i128)Z * bestZ[depth].den == (i128)bestZ[depth].num * Hx) show("maxZ", depth, b, Z, Hx); }
        if (depth < L) { word[depth] = b; dfs(depth + 1, y); }
    }
    free(y);
}
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); if (!f) return 1; if (fscanf(f, "%d", &N) != 1) return 1; L = atoi(argv[3]); if (argc > 4) b0 = atoi(argv[4]); if (b0 > 7) b0 = -1;
    for (int b = 0; b < 8; b++) { Ad[b] = malloc(sizeof(ll) * N * N); for (int i = 0; i < N * N; i++) if (fscanf(f, "%lld", &Ad[b][i]) != 1) return 2; }
    v0 = malloc(sizeof(ll) * N); h = malloc(sizeof(ll) * N);
    for (int i = 0; i < N; i++) if (fscanf(f, "%lld", &v0[i]) != 1) return 3;
    for (int i = 0; i < N; i++) if (fscanf(f, "%lld", &h[i]) != 1) return 4;
    fclose(f);
    FILE *g = fopen(argv[2], "r"); if (!g) return 5; nfrag = malloc(sizeof(int) * N); for (int i = 0; i < N; i++) if (fscanf(g, "%d", &nfrag[i]) != 1) return 6; fclose(g);
    for (int b = 0; b < 8; b++) {
        int m = 0; for (int i = 0; i < N * N; i++) if (Ad[b][i]) m++;
        nnz[b] = m; ri[b] = malloc(sizeof(int) * m); ci[b] = malloc(sizeof(int) * m); val[b] = malloc(sizeof(ll) * m);
        int k = 0; for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (Ad[b][(size_t)i * N + j]) { ri[b][k] = i; ci[b][k] = j; val[b][k] = Ad[b][(size_t)i * N + j]; k++; }
        for (int j = 0; j < N; j++) { ll s = 0; for (int i = 0; i < N; i++) s += h[i] * Ad[b][(size_t)i * N + j]; hA[b][j] = s; twoA[b][j] = (nfrag[j] == 2) ? s : 0; }
    }
    ll *x = malloc(sizeof(ll) * N); memcpy(x, v0, sizeof(ll) * N);
    dfs(0, x);
    if (argc > 5 && strcmp(argv[5], "print") == 0) { printing = 1; printf("all maximising pairs (W, l):\n"); dfs(0, x); printing = 0; }
    for (int n = 0; n <= L; n++) {
        rec *r = &bestH[n], *z = &bestZ[n];
        printf("length %d: words %lld; max H(Wl)/H(W) = %lld/%lld = %.6f, attained by %lld pairs (W,l), first: letter %d, word", n, cnt[n], r->num, r->den, (double)r->num / (double)r->den, r->ties, r->letter);
        for (int i = 0; i < n; i++) printf(" %d", r->word[i]);
        printf("; max Z_l/H = %lld/%lld = %.6f, attained by %lld pairs, first: letter %d, word", z->num, z->den, (double)z->num / (double)z->den, z->ties, z->letter);
        for (int i = 0; i < n; i++) printf(" %d", z->word[i]); printf("\n");
    }
    return 0;
}
