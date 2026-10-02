/*
 * procgen_tcpc_20261001_karp2.c  (TCPC lane, collatz-procgen-20260922, 2026-10-01)
 *
 * Arc-HP parities of a TWO-SHEET CLOCK tournament by inclusion-exclusion over walks
 * (Karp / Bax), over the ring GF(2)[e_1..e_K]/(e_i e_j), with O(1) memory.
 * Algorithmically independent of the subset DPs in procgen_tcpc_20261001_hp.c / _par26.c.
 *
 *   Sum over S subset of V of  (#walks with N vertices inside S, arc weights (1 + sum_t e_t [arc in O_t]))
 *   = Sum over Hamiltonian paths P of prod_(arcs) (1 + ...)       (mod 2; inclusion-exclusion signs vanish)
 *   = H + sum_t e_t * (sum_(e in O_t) c(e))                       (mod 2).
 * O_t is the orbit of an arc under tau: (2j, 2j+1) -> (2j+2, 2j+3) (indices mod N), |O_t| = r odd,
 * so sum_(e in O_t) c(e) = r c(e_t) = c(e_t) (mod 2).  Since every weight is tau-invariant and every
 * tau-orbit of subsets has odd size (it divides r, r odd), the sum over S may be restricted to one
 * representative per tau-orbit (mod 2).
 *
 * stdin: N (even, <= 40), the N x N 0/1 matrix, K, then K lines "u v" (orbit representatives).
 * stdout: "H <H mod 2>" and "C u v <c(u->v) mod 2>" for the K representatives, then "REPS <count>".
 * The input digraph MUST be invariant under tau (checked).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXN 40
#define MAXK 40
static int N, K, r;
static int A[MAXN][MAXN];
static uint64_t Out[MAXN];
static uint64_t OutO[MAXK][MAXN];
static int ru[MAXK], rv[MAXK];
static int NCH;                       /* number of 8-bit chunks */
static uint64_t TabA[5][256];         /* XOR of Out rows over an 8-bit chunk pattern */
static uint64_t TabO[MAXK][5][256];

static inline uint64_t rotS(uint64_t S, uint64_t FULL) { return ((S << 2) | (S >> (N - 2))) & FULL; }

static void build_tab(uint64_t *rows, uint64_t tab[5][256]) {
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

static inline uint64_t matvec(uint64_t tab[5][256], uint64_t f) {
    uint64_t x = 0;
    for (int c = 0; c < NCH; c++) x ^= tab[c][(f >> (8 * c)) & 0xFF];
    return x;
}

int main(void) {
    if (scanf("%d", &N) != 1 || N < 4 || N > MAXN || (N & 1)) { fprintf(stderr, "bad N\n"); return 1; }
    r = N / 2;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
        if (scanf("%d", &A[i][j]) != 1) { fprintf(stderr, "bad matrix\n"); return 1; }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
        if (A[i][j] != A[(i + 2) % N][(j + 2) % N]) { fprintf(stderr, "not tau-invariant\n"); return 1; }
    if (scanf("%d", &K) != 1 || K < 0 || K > MAXK) { fprintf(stderr, "bad K\n"); return 1; }
    for (int t = 0; t < K; t++) {
        if (scanf("%d %d", &ru[t], &rv[t]) != 2 || !A[ru[t]][rv[t]]) { fprintf(stderr, "bad arc\n"); return 1; }
    }
    uint64_t FULL = (N == 64) ? ~0ULL : ((1ULL << N) - 1);
    for (int v = 0; v < N; v++) {
        Out[v] = 0;
        for (int w = 0; w < N; w++) if (v != w && A[v][w]) Out[v] |= 1ULL << w;
    }
    for (int t = 0; t < K; t++) {
        for (int v = 0; v < N; v++) OutO[t][v] = 0;
        for (int j = 0; j < r; j++) {
            int u = (ru[t] + 2 * j) % N, w = (rv[t] + 2 * j) % N;
            OutO[t][u] |= 1ULL << w;
        }
    }
    NCH = (N + 7) / 8;
    build_tab(Out, TabA);
    for (int t = 0; t < K; t++) build_tab(OutO[t], TabO[t]);
    int Hpar = 0;
    int Cpar[MAXK];
    memset(Cpar, 0, sizeof(Cpar));
    uint64_t reps = 0;
    uint64_t g[MAXK], g2[MAXK];
    for (uint64_t S = 1; S <= FULL; S++) {
        /* canonical representative of the tau-orbit: S minimal among its rotations */
        uint64_t R = S;
        int canon = 1;
        for (int j = 1; j < r; j++) {
            R = rotS(R, FULL);
            if (R < S) { canon = 0; break; }
            if (R == S) break;          /* period found; all rotations seen */
        }
        if (!canon) continue;
        reps++;
        uint64_t f = S;
        for (int t = 0; t < K; t++) g[t] = 0;
        for (int step = 1; step < N; step++) {
            uint64_t f2 = matvec(TabA, f) & S;
            for (int t = 0; t < K; t++) g2[t] = (matvec(TabA, g[t]) ^ matvec(TabO[t], f)) & S;
            f = f2;
            for (int t = 0; t < K; t++) g[t] = g2[t];
            if (!f) {
                int any = 0;
                for (int t = 0; t < K; t++) if (g[t]) { any = 1; break; }
                if (!any) break;
            }
        }
        Hpar ^= __builtin_popcountll(f) & 1;
        for (int t = 0; t < K; t++) Cpar[t] ^= __builtin_popcountll(g[t]) & 1;
        if (S == FULL) break;
    }
    printf("H %d\n", Hpar);
    for (int t = 0; t < K; t++) printf("C %d %d %d\n", ru[t], rv[t], Cpar[t]);
    printf("REPS %llu\n", (unsigned long long)reps);
    return 0;
}
