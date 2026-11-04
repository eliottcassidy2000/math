/* uh_graph_downset_20261104.c -- exhaustive census of uniquely Hamiltonian graphs on the fixed cycle 0..n-1.
 * Removing a chord never creates a Hamiltonian cycle, so the chord sets S for which C_n + S has exactly one
 * Hamiltonian cycle form a down-set. We enumerate it by depth-first search, adding chords in increasing index and
 * pruning as soon as a second cycle appears. Every UH graph of order n is isomorphic to C_n + S for such an S, since it
 * can be labelled along its cycle. For each UH graph we test whether two consecutive cycle vertices are both
 * non-critical (G - v Hamiltonian), and record the minimum number of edges among witnesses.
 * Usage: uh_graph_downset_20261104 NMIN NMAX   (OpenMP over the first chord) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

static uint64_t hc_dir(const uint32_t *adj, uint32_t mask, uint64_t cap, uint64_t *dp) {
    int vs[32], k = 0; for (int v = 0; v < 32; v++) if (mask >> v & 1u) vs[k++] = v;
    if (k < 3) return 0;
    uint32_t A[32]; for (int i = 0; i < k; i++) { uint32_t r = 0; for (int j = 0; j < k; j++) if (adj[vs[i]] >> vs[j] & 1u) r |= 1u << j; A[i] = r; }
    size_t M = (size_t)1 << k; memset(dp, 0, M * k * sizeof(uint64_t)); dp[k] = 1;   /* mask 1, end 0 */
    for (size_t m = 1; m < M; m += 2) for (int e = 0; e < k; e++) { uint64_t c = dp[m * k + e]; if (!c) continue;
        uint32_t nx = A[e] & ~(uint32_t)m & (uint32_t)(M - 1);
        while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; uint64_t *d = &dp[(m | ((size_t)1 << w)) * k + w]; *d += c; if (*d > cap) *d = cap; } }
    uint64_t t = 0; for (int e = 1; e < k; e++) if (A[e] & 1u) { t += dp[(M - 1) * k + e]; if (t > cap) t = cap; }
    return t;
}

typedef struct { long uh, witness, minwitness_edges; uint32_t wit[32]; int have; long ncdist[33]; } acc_t;
static int N, NP, PI[512], PJ[512];

static void visit(uint32_t *adj, int next, int edges, acc_t *a, uint64_t *dp) {
    uint32_t full = (1u << N) - 1;
    /* current graph is UH (checked by caller) */
    a->uh++;
    uint32_t NC = 0;
    for (int v = 0; v < N; v++) if (hc_dir(adj, full & ~(1u << v), 1, dp) > 0) NC |= 1u << v;
    a->ncdist[__builtin_popcount(NC)]++;
    int bad = 0; for (int i = 0; i < N; i++) if ((NC >> i & 1u) && (NC >> ((i + 1) % N) & 1u)) bad = 1;
    if (bad) { a->witness++; if (!a->have || edges < a->minwitness_edges) { a->have = 1; a->minwitness_edges = edges; memcpy(a->wit, adj, sizeof(uint32_t) * N); } }
    for (int q = next; q < NP; q++) {
        adj[PI[q]] |= 1u << PJ[q]; adj[PJ[q]] |= 1u << PI[q];
        if (hc_dir(adj, full, 3, dp) == 2) visit(adj, q + 1, edges + 1, a, dp);
        adj[PI[q]] &= ~(1u << PJ[q]); adj[PJ[q]] &= ~(1u << PI[q]);
    }
}

int main(int argc, char **argv) {
    int lo = atoi(argv[1]), hi = atoi(argv[2]);
    for (N = lo; N <= hi; N++) {
        NP = 0; for (int i = 0; i < N; i++) for (int j = i + 2; j < N; j++) { if (i == 0 && j == N - 1) continue; PI[NP] = i; PJ[NP] = j; NP++; }
        acc_t tot; memset(&tot, 0, sizeof tot);
        double t0 = omp_get_wtime();
        /* root (no chord) handled serially; subtrees by first chord in parallel */
        {
            uint64_t *dp = malloc(((size_t)1 << N) * N * sizeof(uint64_t));
            uint32_t adj[32] = {0}; for (int i = 0; i < N; i++) { int j = (i + 1) % N; adj[i] |= 1u << j; adj[j] |= 1u << i; }
            uint32_t full = (1u << N) - 1; uint32_t NC = 0;
            for (int v = 0; v < N; v++) if (hc_dir(adj, full & ~(1u << v), 1, dp) > 0) NC |= 1u << v;
            tot.uh++; tot.ncdist[__builtin_popcount(NC)]++;
            free(dp);
        }
        #pragma omp parallel
        {
            acc_t a; memset(&a, 0, sizeof a);
            uint64_t *dp = malloc(((size_t)1 << N) * N * sizeof(uint64_t));
            #pragma omp for schedule(dynamic, 1)
            for (int q = 0; q < NP; q++) {
                uint32_t adj[32] = {0}; for (int i = 0; i < N; i++) { int j = (i + 1) % N; adj[i] |= 1u << j; adj[j] |= 1u << i; }
                adj[PI[q]] |= 1u << PJ[q]; adj[PJ[q]] |= 1u << PI[q];
                if (hc_dir(adj, (1u << N) - 1, 3, dp) == 2) visit(adj, q + 1, N + 1, &a, dp);
            }
            #pragma omp critical
            {
                tot.uh += a.uh; tot.witness += a.witness;
                for (int c = 0; c <= 32; c++) tot.ncdist[c] += a.ncdist[c];
                if (a.have && (!tot.have || a.minwitness_edges < tot.minwitness_edges)) { tot.have = 1; tot.minwitness_edges = a.minwitness_edges; memcpy(tot.wit, a.wit, sizeof a.wit); }
            }
            free(dp);
        }
        printf("n=%2d labelled UH graphs on the fixed cycle=%ld  with two consecutive non-critical vertices=%ld  (%.1fs)  |NC| dist:", N, tot.uh, tot.witness, omp_get_wtime() - t0);
        for (int c = 0; c <= N; c++) if (tot.ncdist[c]) printf(" %d:%ld", c, tot.ncdist[c]);
        if (tot.have) { printf("\n   a witness with the fewest edges (%ld edges), rows:", tot.minwitness_edges); for (int i = 0; i < N; i++) printf(" %x", tot.wit[i]); }
        printf("\n"); fflush(stdout);
    }
    return 0;
}
