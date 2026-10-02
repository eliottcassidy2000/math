/*
 * procgen_tcpc_20261001_karp4.c  (TCPC lane, collatz-procgen-20260922, 2026-10-01)
 *
 * Same inclusion-exclusion parity computation as procgen_tcpc_20261001_karp3.c (see its header),
 * but the per-orbit cross-correlations
 *        C_t(S) = XOR_k parity( f_k & Src_t & rot(b_(N-2-k), delta_t) )
 * are obtained for ALL offsets delta at once by carry-less multiplication over GF(2)
 * (ARM64 PMULL, vmull_p64):  sum_v F[v] B[v + delta] is the coefficient of x^(delta + N - 1)
 * in clmul(B, rev_N(F)) folded modulo x^N - 1.  The two source sheets (even/odd vertices) are
 * correlated separately.  Requires an ARM64 CPU with the crypto extension (build with
 * -march=armv8-a+crypto).
 *
 * stdin/stdout and the optional "part total" arguments are as in karp3.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <arm_neon.h>

#define MAXN 62
#define MAXK 64
static int N, K, r;
static int A[MAXN][MAXN];
static uint64_t Out[MAXN], In[MAXN];
static int ru[MAXK], rv[MAXK], dl[MAXK], sh[MAXK];
static int NCH;
static uint64_t TabO[8][256], TabI[8][256];

static void build_tab(uint64_t *rows, uint64_t tab[8][256]) {
    for (int c = 0; c < NCH; c++)
        for (int pat = 0; pat < 256; pat++) {
            uint64_t x = 0;
            for (int b = 0; b < 8; b++) {
                int v = 8 * c + b;
                if (v < N && (pat >> b & 1)) x ^= rows[v];
            }
            tab[c][pat] = x;
        }
}

static inline uint64_t mv(uint64_t tab[8][256], uint64_t f) {
    uint64_t x = 0;
    for (int c = 0; c < NCH; c++) x ^= tab[c][(f >> (8 * c)) & 0xFF];
    return x;
}

static inline void clmul_acc(uint64_t a, uint64_t b, uint64_t *lo, uint64_t *hi) {
    poly128_t p = vmull_p64((poly64_t)a, (poly64_t)b);
    *lo ^= (uint64_t)p;
    *hi ^= (uint64_t)(p >> 64);
}

int main(int argc, char **argv) {
    uint64_t part = 0, total = 1;
    if (argc >= 3) { part = strtoull(argv[1], NULL, 10); total = strtoull(argv[2], NULL, 10); }
    if (scanf("%d", &N) != 1 || N < 4 || N > MAXN || (N & 1)) { fprintf(stderr, "bad N\n"); return 1; }
    r = N / 2;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
        if (scanf("%d", &A[i][j]) != 1) { fprintf(stderr, "bad matrix\n"); return 1; }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
        if (A[i][j] != A[(i + 2) % N][(j + 2) % N]) { fprintf(stderr, "not tau-invariant\n"); return 1; }
    if (scanf("%d", &K) != 1 || K < 0 || K > MAXK) { fprintf(stderr, "bad K\n"); return 1; }
    uint64_t FULL = (1ULL << N) - 1;
    uint64_t EVEN = 0;
    for (int v = 0; v < N; v += 2) EVEN |= 1ULL << v;
    uint64_t ODD = FULL ^ EVEN;
    for (int t = 0; t < K; t++) {
        if (scanf("%d %d", &ru[t], &rv[t]) != 2 || !A[ru[t]][rv[t]]) { fprintf(stderr, "bad arc\n"); return 1; }
        dl[t] = ((rv[t] - ru[t]) % N + N) % N;
        sh[t] = ru[t] & 1;
    }
    for (int v = 0; v < N; v++) {
        Out[v] = 0; In[v] = 0;
        for (int w = 0; w < N; w++) {
            if (v != w && A[v][w]) Out[v] |= 1ULL << w;
            if (v != w && A[w][v]) In[v] |= 1ULL << w;
        }
    }
    NCH = (N + 7) / 8;
    build_tab(Out, TabO);
    build_tab(In, TabI);
    int Hpar = 0;
    /* global XOR accumulators of the folded correlations, per sheet */
    uint64_t GE = 0, GO = 0;
    uint64_t reps = 0, idx = 0;
    uint64_t f[MAXN], b[MAXN];
    const int RS = 64 - N;
    for (uint64_t S = 1; S <= FULL; S++) {
        uint64_t R = S;
        int canon = 1;
        for (int j = 1; j < r; j++) {
            R = ((R << 2) | (R >> (N - 2))) & FULL;
            if (R < S) { canon = 0; break; }
            if (R == S) break;
        }
        if (!canon) continue;
        reps++;
        if (total > 1 && (idx++ % total) != part) continue;
        f[0] = S; b[0] = S;
        for (int k = 1; k < N; k++) {
            f[k] = mv(TabO, f[k - 1]) & S;
            b[k] = mv(TabI, b[k - 1]) & S;
        }
        Hpar ^= __builtin_popcountll(f[N - 1]) & 1;
        uint64_t eLo = 0, eHi = 0, oLo = 0, oHi = 0;
        for (int k = 0; k <= N - 2; k++) {
            uint64_t B = b[N - 2 - k];
            uint64_t FE = __builtin_bitreverse64(f[k] & EVEN) >> RS;   /* rev_N */
            uint64_t FO = __builtin_bitreverse64(f[k] & ODD) >> RS;
            if (FE) clmul_acc(B, FE, &eLo, &eHi);
            if (FO) clmul_acc(B, FO, &oLo, &oHi);
        }
        /* fold modulo x^N - 1: bit j of fold = bit j xor bit j + N of the 128-bit product */
        uint64_t foldE = (eLo ^ ((eLo >> N) | (eHi << (64 - N)))) & FULL;
        uint64_t foldO = (oLo ^ ((oLo >> N) | (oHi << (64 - N)))) & FULL;
        GE ^= foldE;
        GO ^= foldO;
        if (S == FULL) break;
    }
    printf("H %d\n", Hpar);
    for (int t = 0; t < K; t++) {
        int pos = (dl[t] + N - 1) % N;
        uint64_t G = sh[t] ? GO : GE;
        printf("C %d %d %d\n", ru[t], rv[t], (int)((G >> pos) & 1));
    }
    printf("REPS %llu\n", (unsigned long long)reps);
    return 0;
}
