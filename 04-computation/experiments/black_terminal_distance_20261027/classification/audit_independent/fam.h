/* fam.h -- bottom layers and E_t(S, x0) built from the definition (shared by famE.c and perturb.c) */
#ifndef FAM_H
#define FAM_H
#include "at.h"
typedef struct { int a; int q; int sz[16]; int bit[16]; } blayer;   /* bit[i] = 1 if x0 -> R_i */

static int nBL = 0; static blayer BL[4096];
static void gen_comp(int rem, int *parts, int k, int a) {
    if (rem == 0) {
        if (k == 0) return;                                  /* need q - 1 >= 1 blocks before R_q */
        int q = k + 1; int nb = q >= 3 ? q - 2 : 0;
        for (int bits = 0; bits < (1 << nb); bits++) {
            blayer B; memset(&B, 0, sizeof B); B.a = a; B.q = q;
            for (int i = 0; i < k; i++) B.sz[i] = parts[i];
            B.sz[k] = 3;
            B.bit[0] = 1;                                     /* x0 => R_1 */
            for (int i = 1; i <= q - 2; i++) B.bit[i] = (bits >> (i - 1)) & 1;
            B.bit[q - 1] = 0;                                 /* R_q => x0 */
            BL[nBL++] = B;
        }
        return;
    }
    if (rem >= 1) { parts[k] = 1; gen_comp(rem - 1, parts, k + 1, a); }
    if (rem >= 3) { parts[k] = 3; gen_comp(rem - 3, parts, k + 1, a); }
}

/* S on a vertices: x0 = 0, blocks consecutive; returns block start indices */
static void build_S(const blayer *B, tour *S, int *bstart) {
    int a = B->a; t_clear(S, a);
    int idx = 1;
    for (int i = 0; i < B->q; i++) { bstart[i] = idx; idx += B->sz[i]; }
    bstart[B->q] = idx;
    if (idx != a) { fprintf(stderr, "size mismatch\n"); exit(1); }
    for (int i = 0; i < B->q; i++) {
        int b0 = bstart[i];
        if (B->sz[i] == 3) { t_arc(S, b0, b0 + 1); t_arc(S, b0 + 1, b0 + 2); t_arc(S, b0 + 2, b0); }
        for (int j = i + 1; j < B->q; j++) for (int u = bstart[i]; u < bstart[i + 1]; u++) for (int v = bstart[j]; v < bstart[j + 1]; v++) t_arc(S, u, v);
        for (int u = bstart[i]; u < bstart[i + 1]; u++) { if (B->bit[i]) t_arc(S, 0, u); else t_arc(S, u, 0); }
    }
    t_fin(S);
}

/* E_t(S, x0=0): bottom 0..a-1, middles c_k = a+k-1 (k=1..t-1), top sigma(u) = a+t-1+u */
static void build_E(const tour *S, int t, tour *T, int *sigma) {
    int a = S->n, n = 2 * a + t - 1; t_clear(T, n);
    int top0 = a + t - 1;
    /* default: higher layer -> lower layer; order of layers: bottom (0), middles by k, top */
    int layer[MAXV];
    for (int u = 0; u < a; u++) layer[u] = 0;
    for (int k = 1; k <= t - 1; k++) layer[a + k - 1] = k;
    for (int u = 0; u < a; u++) layer[top0 + u] = t;
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) {
        if (layer[u] != layer[v]) { if (layer[u] > layer[v]) t_arc(T, u, v); else t_arc(T, v, u); }
    }
    for (int u = 0; u < a; u++) for (int v = 0; v < a; v++) if (u != v && t_has(S, u, v)) { t_arc(T, u, v); t_arc(T, top0 + v, top0 + u); }
    /* jumps */
    if (t == 1) t_arc(T, 0, top0 + 0);
    else { t_arc(T, 0, a); for (int k = 1; k + 1 <= t - 1; k++) t_arc(T, a + k - 1, a + k); t_arc(T, a + t - 2, top0 + 0); }
    t_fin(T);
    for (int u = 0; u < a; u++) { sigma[u] = top0 + u; sigma[top0 + u] = u; }
    for (int k = 1; k <= t - 1; k++) sigma[a + k - 1] = a + (t - k) - 1;
}

#endif
