/* a7_anneal.c (AUDIT, independent) - search for a tournament on the 15 points of PG(3,2) with FEWER than 21 even
 * zero-sum 4-sets.  The cost is computed directly: for each zero-sum 4-set the number of cyclic triangles inside
 * it is counted from the adjacency matrix (even = 0 or 2 triangles).  No linear algebra, no planes.
 * Simulated annealing with restarts + greedy polishing; several schedules.  usage: a7_anneal restarts steps seed */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned long long u64;
static int adj[16][16];
static int quads[105][4], nq;
static int qof[16][16][6], nqof[16][16];
static u64 rs;
static inline u64 rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline int cyc(int a, int b, int c) { return (adj[a][b] && adj[b][c] && adj[c][a]) || (adj[b][a] && adj[c][b] && adj[a][c]); }
static inline int even(int q) { int *p = quads[q]; int t = cyc(p[0], p[1], p[2]) + cyc(p[0], p[1], p[3]) + cyc(p[0], p[2], p[3]) + cyc(p[1], p[2], p[3]); return t != 1; }
int main(int argc, char **argv) {
    int restarts = atoi(argv[1]); long steps = atol(argv[2]); rs = 88172645463325252ULL ^ (u64)atoll(argv[3]) * 0x9E3779B97F4A7C15ULL;
    for (int a = 1; a < 16; a++) for (int b = a + 1; b < 16; b++) for (int c = b + 1; c < 16; c++) { int d = a ^ b ^ c; if (d > c) { quads[nq][0] = a; quads[nq][1] = b; quads[nq][2] = c; quads[nq][3] = d; nq++; } }
    if (nq != 105) { printf("nq = %d\n", nq); return 1; }
    for (int q = 0; q < nq; q++) for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) if (i != j) { int u = quads[q][i], v = quads[q][j]; qof[u][v][nqof[u][v]++] = q; }
    for (int u = 1; u < 16; u++) for (int v = 1; v < 16; v++) if (u != v && nqof[u][v] != 6) { printf("bad incidence\n"); return 1; }
    long hist[106]; memset(hist, 0, sizeof hist); int best_all = 999; int bestadj[16][16];
    for (int rep = 0; rep < restarts; rep++) {
        for (int u = 1; u < 16; u++) for (int v = u + 1; v < 16; v++) { int b = (int)(rnd() >> 33) & 1; adj[u][v] = b; adj[v][u] = !b; }
        int cur = 0; for (int q = 0; q < nq; q++) cur += even(q);
        int best = cur; double T = 1.2 + 0.6 * (rep % 4), cool = pow(0.03 / T, 1.0 / (double)steps);
        for (long st = 0; st < steps; st++) {
            int u = 1 + (int)(rnd() % 15), v = 1 + (int)(rnd() % 15); if (u == v) continue;
            int before = 0, after = 0; for (int i = 0; i < 6; i++) before += even(qof[u][v][i]);
            adj[u][v] ^= 1; adj[v][u] ^= 1;
            for (int i = 0; i < 6; i++) after += even(qof[u][v][i]);
            int d = after - before;
            if (d <= 0 || (double)(rnd() >> 11) * (1.0 / 9007199254740992.0) < exp(-d / T)) { cur += d; if (cur < best) { best = cur; if (best < best_all) { best_all = best; memcpy(bestadj, adj, sizeof adj); } } }
            else { adj[u][v] ^= 1; adj[v][u] ^= 1; }
            T *= cool;
        }
        hist[best]++;
    }
    printf("restarts %d, steps %ld: least number of even 4-sets per restart -> count:", restarts, steps);
    for (int e = 0; e < 106; e++) if (hist[e]) printf(" %d:%ld", e, hist[e]); printf("\noverall least: %d\n", best_all);
    /* recount the best from scratch and print it */
    memcpy(adj, bestadj, sizeof adj); int c = 0; for (int q = 0; q < nq; q++) c += even(q);
    printf("best tournament recounted: %d even 4-sets; out-masks (points 1..15):", c);
    for (int u = 1; u < 16; u++) { int m = 0; for (int v = 1; v < 16; v++) if (adj[u][v]) m |= 1 << (v - 1); printf(" %d", m); } printf("\n");
    return 0;
}
