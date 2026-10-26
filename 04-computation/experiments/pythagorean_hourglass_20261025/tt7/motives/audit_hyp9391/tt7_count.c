/* Independent audit counter for TT7(P_q), q = 3 (mod 4) prime (audit of HYP-9391).
 *
 * Method (different from the lane's paley_tt7.py, which fixes the source and counts TT6 in the QR tournament):
 * the affine maps x -> a x + b with a a nonzero square act simply transitively on the arcs of P_q, and every
 * transitive 7-set has a unique (source, second vertex) arc, so
 *     TT7(P_q) = (q(q-1)/2) * T5(D),   D = N+(0) cap N+(1) = { x : chi(x) = chi(x-1) = 1 },
 * where T5(D) is the number of transitive 5-subsets of the subtournament on D.  T5 is counted as
 *     T5 = sum over 3-chains v1 -> v2 -> v3 in D of C(d, 2),  d = |D cap N+(v1) cap N+(v2) cap N+(v3)|
 * (any 2-subset of the common out-set is a 2-chain).  Output: "q T5".  A brute-force mode ("brute q") counts
 * transitive 7-subsets of P_q directly for small q.
 * Usage: tt7_count q1 q2 ...   |   tt7_count brute q
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int *chi;
static void make_chi(int q) {
    chi = calloc(q, sizeof(int));
    for (int x = 1; x < q; x++) chi[x] = -1;
    for (long x = 1; x < q; x++) chi[(x * x) % q] = 1;
}

static uint64_t T5_of_D(int q) {
    int *D = malloc(q * sizeof(int)); int n = 0;
    for (int x = 2; x < q; x++) if (chi[x] == 1 && chi[x - 1] == 1) D[n++] = x;
    int W = (n + 63) / 64;
    uint64_t *A = calloc((size_t)n * W, sizeof(uint64_t));      /* A[i] = out-set of D[i] inside D */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j && chi[((D[j] - D[i]) % q + q) % q] == 1) A[(size_t)i * W + j / 64] |= 1ULL << (j % 64);
    uint64_t *S1 = malloc(W * 8), *S2 = malloc(W * 8);
    uint64_t total = 0;
    for (int v1 = 0; v1 < n; v1++) {
        memcpy(S1, A + (size_t)v1 * W, W * 8);
        for (int w1 = 0; w1 < W; w1++) {
            uint64_t b1 = S1[w1];
            while (b1) {
                int v2 = w1 * 64 + __builtin_ctzll(b1); b1 &= b1 - 1;
                const uint64_t *A2 = A + (size_t)v2 * W;
                for (int k = 0; k < W; k++) S2[k] = S1[k] & A2[k];
                for (int w2 = 0; w2 < W; w2++) {
                    uint64_t b2 = S2[w2];
                    while (b2) {
                        int v3 = w2 * 64 + __builtin_ctzll(b2); b2 &= b2 - 1;
                        const uint64_t *A3 = A + (size_t)v3 * W;
                        uint64_t d = 0;
                        for (int k = 0; k < W; k++) d += __builtin_popcountll(S2[k] & A3[k]);
                        total += d * (d - 1) / 2;
                    }
                }
            }
        }
    }
    free(D); free(A); free(S1); free(S2);
    return total;
}

/* brute force: transitive 7-subsets of P_q by ordered chains (each transitive set has exactly one) */
static uint64_t brute_count;
static int q_g;
static void ext(int depth, int *cand, int nc) {
    if (depth == 7) { brute_count++; return; }
    int *nxt = malloc(nc * sizeof(int));
    for (int i = 0; i < nc; i++) {
        int v = cand[i], m = 0;
        for (int j = 0; j < nc; j++) if (chi[((cand[j] - v) % q_g + q_g) % q_g] == 1) nxt[m++] = cand[j];
        ext(depth + 1, nxt, m);
    }
    free(nxt);
}

int main(int argc, char **argv) {
    if (argc >= 3 && strcmp(argv[1], "brute") == 0) {
        int q = atoi(argv[2]); q_g = q; make_chi(q);
        int *all = malloc(q * sizeof(int)); for (int i = 0; i < q; i++) all[i] = i;
        brute_count = 0; ext(0, all, q);
        printf("%d %llu\n", q, (unsigned long long)brute_count);
        return 0;
    }
    for (int a = 1; a < argc; a++) {
        int q = atoi(argv[a]); make_chi(q);
        uint64_t t5 = T5_of_D(q);
        printf("%d %llu\n", q, (unsigned long long)t5); fflush(stdout);
        free(chi);
    }
    return 0;
}
