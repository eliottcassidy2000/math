/* uh_graph_noncritical_20261104.c -- non-critical vertices of uniquely Hamiltonian GRAPHS (undirected).
 *
 * A graph G is uniquely Hamiltonian (UH) if it has exactly one Hamiltonian cycle C; v is non-critical if G - v is
 * Hamiltonian.  Question (HYP-9462 (B)): can two consecutive vertices of C both be non-critical?  For tournaments
 * the non-critical vertices form ONE run of C (THM-4868 (C)); for graphs the conjecture is the opposite: never two
 * consecutive ones.
 *
 * Sampler: G = C_n (vertices 0..n-1 in cycle order) plus each non-cycle pair independently with probability p.
 * Keep G if it has exactly one Hamiltonian cycle (exact subset DP, undirected cycles = directed / 2).
 * Rejection loop is capped (MAXTRIALS) and reports its yield, so an infeasible parameter set cannot hang.
 * For every kept G: compute the non-critical set NC and test
 *   consec : some i with c_i and c_(i+1) both in NC   (the conjecture says never)
 *   adj    : some two NC vertices adjacent in G       (allowed; reported)
 * Usage: uh_graph_noncritical_20261104 NMIN NMAX TARGET_PER_N SEED P1 [P2 ...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static uint64_t rng_s;
static uint64_t rnd(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }
static double unif(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

/* number of directed Hamiltonian cycles of the graph induced on `mask` (adjacency rows adj), capped at cap */
static uint64_t *dpbuf; static size_t dpcap;
static uint64_t hcycles_dir(const uint32_t *adj, uint32_t mask, uint64_t cap) {
    int vs[32], k = 0; for (int v = 0; v < 32; v++) if (mask >> v & 1u) vs[k++] = v;
    if (k < 3) return 0;
    uint32_t A[32]; for (int i = 0; i < k; i++) { uint32_t r = 0; for (int j = 0; j < k; j++) if (adj[vs[i]] >> vs[j] & 1u) r |= 1u << j; A[i] = r; }
    size_t M = (size_t)1 << k; size_t need = M * k;
    if (need > dpcap) { free(dpbuf); dpbuf = malloc(need * sizeof(uint64_t)); dpcap = need; }
    memset(dpbuf, 0, need * sizeof(uint64_t));
    dpbuf[1 * k + 0] = 1;
    for (size_t m = 1; m < M; m += 2) for (int e = 0; e < k; e++) { uint64_t c = dpbuf[m * k + e]; if (!c) continue;
        uint32_t nx = A[e] & ~(uint32_t)m & (uint32_t)(M - 1);
        while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; uint64_t *d = &dpbuf[(m | ((size_t)1 << w)) * k + w]; *d += c; if (*d > cap) *d = cap; } }
    uint64_t tot = 0; for (int e = 1; e < k; e++) if (A[e] & 1u) { tot += dpbuf[(M - 1) * k + e]; if (tot > cap) tot = cap; }
    return tot;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: NMIN NMAX TARGET_PER_N SEED P1 [P2 ...]\n"); return 1; }
    int nmin = atoi(argv[1]), nmax = atoi(argv[2]); long target = atol(argv[3]); rng_s = strtoull(argv[4], NULL, 10) | 1ULL;
    int np = argc - 5; double ps[16]; for (int i = 0; i < np && i < 16; i++) ps[i] = atof(argv[5 + i]);
    const long MAXTRIALS = 20000000;
    long tot_uh = 0, tot_consec = 0, tot_adj = 0; int maxnc = 0;
    printf("# uh_graph_noncritical_20261104 n=%d..%d target/n=%ld seed=%s p:", nmin, nmax, target, argv[4]);
    for (int i = 0; i < np; i++) printf(" %.3f", ps[i]); printf("\n");
    for (int n = nmin; n <= nmax; n++) {
        long kept = 0, trials = 0, consec = 0, adjnc = 0; long ncdist[33]; memset(ncdist, 0, sizeof ncdist);
        uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
        /* exhaustive mode (TARGET_PER_N == 0): every chord set on the fixed cycle 0..n-1; every UH graph of order n
           is isomorphic to one of these (label it along its unique cycle) */
        int npairs = 0, pi[512], pj[512];
        for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { if (i == 0 && j == n - 1) continue; pi[npairs] = i; pj[npairs] = j; npairs++; }
        int exhaustive = (target == 0);
        if (exhaustive && npairs > 24) { printf("n=%d exhaustive mode refused (2^%d chord sets)\n", n, npairs); continue; }
        uint64_t nsets = exhaustive ? ((uint64_t)1 << npairs) : 0;
        while (exhaustive ? ((uint64_t)trials < nsets) : (kept < target && trials < MAXTRIALS)) {
            uint64_t code = (uint64_t)trials;
            trials++;
            double p = ps[trials % np];
            uint32_t adj[32]; memset(adj, 0, sizeof adj);
            for (int i = 0; i < n; i++) { int j = (i + 1) % n; adj[i] |= 1u << j; adj[j] |= 1u << i; }
            if (exhaustive) { for (int q = 0; q < npairs; q++) if (code >> q & 1ULL) { adj[pi[q]] |= 1u << pj[q]; adj[pj[q]] |= 1u << pi[q]; } }
            else for (int q = 0; q < npairs; q++) if (unif() < p) { adj[pi[q]] |= 1u << pj[q]; adj[pj[q]] |= 1u << pi[q]; }
            if (hcycles_dir(adj, full, 3) != 2) continue;          /* exactly one undirected Hamiltonian cycle */
            kept++;
            uint32_t NC = 0;
            for (int v = 0; v < n; v++) if (hcycles_dir(adj, full & ~(1u << v), 1) > 0) NC |= 1u << v;
            int c = __builtin_popcount(NC); ncdist[c]++; if (c > maxnc) maxnc = c;
            int bad = 0, ad = 0;
            for (int i = 0; i < n; i++) if ((NC >> i & 1u) && (NC >> ((i + 1) % n) & 1u)) bad = 1;
            for (int u = 0; u < n; u++) if (NC >> u & 1u) if (adj[u] & NC) ad = 1;
            if (bad) { consec++; if (consec <= 3) { printf("COUNTEREXAMPLE n=%d rows:", n); for (int i = 0; i < n; i++) printf(" %x", adj[i]); printf(" NC=%x\n", NC); } }
            if (ad) adjnc++;
        }
        printf("n=%2d kept=%ld trials=%ld yield=%.4f consecutive_NC_pairs=%ld adjacent_NC_pairs=%ld  |NC| dist:", n, kept, trials, trials ? (double)kept / trials : 0.0, consec, adjnc);
        for (int c = 0; c <= n; c++) if (ncdist[c]) printf(" %d:%ld", c, ncdist[c]);
        if (kept < target) printf("  [CAP REACHED: infeasible or low-yield parameters]");
        printf("\n"); fflush(stdout);
        tot_uh += kept; tot_consec += consec; tot_adj += adjnc;
    }
    printf("TOTAL UH graphs=%ld consecutive_NC_pairs=%ld adjacent_NC_pairs=%ld max|NC|=%d\n", tot_uh, tot_consec, tot_adj, maxnc);
    return 0;
}
