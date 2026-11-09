/* d2small.c - audit D2 (THM-5056): exhaustive enumeration of ALL labelled graphs on n <= NMAX vertices (default 7),
 * written from scratch (no census input, no canonical forms).
 *
 * For every labelled connected graph: d(G) by subset enumeration, and the complete automorphism group by
 * backtracking, with the sign (-1)^{#edges a<b with g(a)>g(b)} of every element.  The number of isomorphism classes
 * with given (n, d, parity) is (sum of |Aut(G)| over the labelled graphs in that bucket) / n!  (orbit-stabiliser:
 * a class with automorphism group A has n!/|A| labelled copies).
 *
 * Output:
 *   "C n d parity classes"  for every (n, d, parity) with connected graphs;
 *   "MIN n mind"            the least d of a connected graph on n vertices;
 *   "REP n d parity |Aut| edges"  one canonical representative (lexicographically least adjacency word over all
 *                           relabellings) of every connected class with d in the watch list (default 5, 37), n <= 6.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define NMAX 7
static int n; static uint32_t adj[NMAX];
static int img[NMAX], used[NMAX], dg[NMAX];
static long long naut; static int anyodd;
static int edge(int a, int b) { return (adj[a] >> b) & 1; }
static void bt(int v) {
    if (v == n) {
        naut++;
        int rev = 0;
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (edge(a, b) && img[a] > img[b]) rev++;
        if (rev & 1) anyodd = 1;
        return;
    }
    for (int c = 0; c < n; c++) {
        if (used[c] || dg[c] != dg[v]) continue;
        int ok = 1;
        for (int j = 0; j < v; j++) if (edge(v, j) != edge(c, img[j])) { ok = 0; break; }
        if (!ok) continue;
        img[v] = c; used[c] = 1; bt(v + 1); used[c] = 0;
    }
}
static int dcount(void) {
    int c = 0; uint32_t full = (1u << n) - 1;
    for (uint32_t S = 1; S < (1u << n); S++) {
        uint32_t dom = 0;
        for (int v = 0; v < n; v++) if ((S >> v) & 1) dom |= adj[v] | (1u << v);
        if (dom == full) c++;
    }
    return c;
}
static int connected(void) {
    uint32_t seen = 1, frontier = 1;
    while (frontier) { uint32_t nf = 0; for (int v = 0; v < n; v++) if ((frontier >> v) & 1) nf |= adj[v]; frontier = nf & ~seen; seen |= nf; }
    return seen == (1u << n) - 1;
}
/* canonical word: least, over all permutations, of the upper-triangle bit string read in a fixed order */
static uint64_t best; static int perm[NMAX], pused[NMAX];
static void canon_rec(int i) {
    if (i == n) {
        uint64_t w = 0; int bit = 0;
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) { if (edge(perm[a], perm[b])) w |= 1ull << bit; bit++; }
        if (w < best) best = w;
        return;
    }
    for (int c = 0; c < n; c++) if (!pused[c]) { pused[c] = 1; perm[i] = c; canon_rec(i + 1); pused[c] = 0; }
}
#define DMAX 128
static long long sumaut[NMAX + 1][DMAX + 1][2];
int main(int argc, char **argv) {
    int nmax = argc > 1 ? atoi(argv[1]) : 7;
    int watch[8] = {5, 37}; int nwatch = 2;
    if (nmax > NMAX) nmax = NMAX;
    for (n = 1; n <= nmax; n++) {
        int E = n * (n - 1) / 2; int ea[32], eb[32], t = 0;
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) { ea[t] = a; eb[t] = b; t++; }
        int mind = 1 << 30;
        static uint64_t reps[4096]; static int nrep; nrep = 0;
        for (uint64_t mask = 0; mask < (1ull << E); mask++) {
            memset(adj, 0, sizeof adj);
            for (int e = 0; e < E; e++) if ((mask >> e) & 1) { adj[ea[e]] |= 1u << eb[e]; adj[eb[e]] |= 1u << ea[e]; }
            if (!connected()) continue;
            int d = dcount();
            if (d < mind) mind = d;
            for (int v = 0; v < n; v++) { dg[v] = __builtin_popcount(adj[v]); used[v] = 0; }
            naut = 0; anyodd = 0; bt(0);
            sumaut[n][d][anyodd] += naut;
            int w = 0; for (int i = 0; i < nwatch; i++) if (d == watch[i]) w = 1;
            if (w && n <= 6) {
                best = ~0ull; memset(pused, 0, sizeof pused); canon_rec(0);
                int f = 0; for (int i = 0; i < nrep; i++) if (reps[i] == best) { f = 1; break; }
                if (!f) {
                    reps[nrep++] = best;
                    /* print the canonical representative */
                    uint32_t save[NMAX]; memcpy(save, adj, sizeof adj);
                    memset(adj, 0, sizeof adj); int bit = 0;
                    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) { if ((best >> bit) & 1) { adj[a] |= 1u << b; adj[b] |= 1u << a; } bit++; }
                    for (int v = 0; v < n; v++) { dg[v] = __builtin_popcount(adj[v]); used[v] = 0; }
                    naut = 0; anyodd = 0; bt(0);
                    printf("REP %d %d %s %lld |", n, d, anyodd ? "O" : "E", naut);
                    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (edge(a, b)) printf(" %d-%d", a, b);
                    printf("\n");
                    memcpy(adj, save, sizeof adj);
                }
            }
        }
        long long fact = 1; for (int i = 2; i <= n; i++) fact *= i;
        for (int d = 0; d <= DMAX; d++) for (int p = 0; p < 2; p++) if (sumaut[n][d][p]) {
            if (sumaut[n][d][p] % fact) { printf("ERROR non-integral class count n=%d d=%d\n", n, d); return 1; }
            printf("C %d %d %s %lld\n", n, d, p ? "O" : "E", sumaut[n][d][p] / fact);
        }
        printf("MIN %d %d\n", n, mind);
        fflush(stdout);
    }
    return 0;
}
