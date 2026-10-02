/*
 * procgen_tcpc_20261001_karp3.c  (TCPC lane, collatz-procgen-20260922, 2026-10-01)
 *
 * Faster variant of procgen_tcpc_20261001_karp2.c (same inclusion-exclusion principle, O(1) memory):
 * for every subset S (one per tau-orbit) it computes forward walk-parity vectors f_k (walks with k+1
 * vertices inside S ending at v) and backward vectors b_k (walks with k+1 vertices inside S starting
 * at w), and gets, for every tau-orbit O_t of arcs (u_t + 2j -> v_t + 2j), the parity of the number of
 * N-vertex walks inside S that use a marked arc of O_t:
 *        C_t(S) = XOR_k parity( f_k & Src_t & rot(b_(N-2-k), delta_t) ),   delta_t = v_t - u_t,
 * where Src_t is the sheet of u_t.  Summed over S (mod 2) this is sum_(e in O_t) c(e) = r c(e_t)
 * = c(e_t) (mod 2), r = N/2 odd (inclusion-exclusion kills non-Hamiltonian walks).
 *
 * stdin: N (even, <= 62), the N x N 0/1 matrix (tau-invariant: tau = +2 mod N), K, K lines "u v".
 * stdout: "H <H mod 2>", "C u v <c(u->v) mod 2>" per representative, "REPS <count>".
 * Optional argv[1], argv[2] = "part total": only process canonical subsets S with (S % total == part)
 * (for splitting; parities of the parts XOR together).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXN 62
#define MAXK 64
static int N, K, r;
static int A[MAXN][MAXN];
static uint64_t Out[MAXN], In[MAXN];
static int ru[MAXK], rv[MAXK], dl[MAXK];
static uint64_t Src[MAXK];
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
    for (int t = 0; t < K; t++) {
        if (scanf("%d %d", &ru[t], &rv[t]) != 2 || !A[ru[t]][rv[t]]) { fprintf(stderr, "bad arc\n"); return 1; }
        dl[t] = ((rv[t] - ru[t]) % N + N) % N;
        Src[t] = (ru[t] % 2 == 0) ? EVEN : (FULL ^ EVEN);
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
    int Cpar[MAXK];
    memset(Cpar, 0, sizeof(Cpar));
    uint64_t reps = 0;
    uint64_t f[MAXN], b[MAXN];
    for (uint64_t S = 1; S <= FULL; S++) {
        uint64_t R = S;
        int canon = 1;
        for (int j = 1; j < r; j++) {
            R = ((R << 2) | (R >> (N - 2))) & FULL;
            if (R < S) { canon = 0; break; }
            if (R == S) break;
        }
        if (!canon) continue;
        if (total > 1 && (reps++ % total) != part) continue;
        if (total == 1) reps++;
        f[0] = S; b[0] = S;
        for (int k = 1; k < N; k++) {
            f[k] = mv(TabO, f[k - 1]) & S;      /* walks extended forward: union of out-rows */
            b[k] = mv(TabI, b[k - 1]) & S;      /* walks extended backward: union of in-rows */
        }
        Hpar ^= __builtin_popcountll(f[N - 1]) & 1;
        for (int t = 0; t < K; t++) {
            int d = dl[t];
            uint64_t acc = 0;
            for (int k = 0; k <= N - 2; k++) {
                uint64_t bb = b[N - 2 - k];
                /* rot(bb, d)[v] = bb[v + d] */
                uint64_t rb = d ? (((bb >> d) | (bb << (N - d))) & FULL) : bb;
                acc ^= f[k] & Src[t] & rb;
            }
            Cpar[t] ^= __builtin_popcountll(acc) & 1;
        }
        if (S == FULL) break;
    }
    printf("H %d\n", Hpar);
    for (int t = 0; t < K; t++) printf("C %d %d %d\n", ru[t], rv[t], Cpar[t]);
    printf("REPS %llu\n", (unsigned long long)reps);
    return 0;
}
