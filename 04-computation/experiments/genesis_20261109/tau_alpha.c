/* tau_alpha.c (= tau_share.c + ALPHA dump of first orders for N <= 120) - pre-registered Q3 (genesis lane): the spanning-tree analogue of the first-creator share.
 * For m = 2..MMAX, over ALL labelled simple graphs on m vertices: the set S_tau(m) of spanning-tree counts of
 * connected graphs, and the labelled share of connected graphs whose count is NEW (not in S_tau(m-1)).
 * tau(G) = det of the reduced Laplacian (Kirchhoff), computed exactly by fraction-free Bareiss elimination.
 * S_tau(m-1) is contained in S_tau(m) (attach a pendant vertex), like S(n-1) in S(n) for tournaments.
 * usage: tau_share MMAX   (MMAX <= 8)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
#define MAXV 300000   /* tau <= m^(m-2) = 262144 for m = 8 */
static unsigned char seenPrev[MAXV], seenCur[MAXV];

static ll det_bareiss(ll a[8][8], int n) {
    ll prev = 1; int sign = 1;
    for (int k = 0; k < n - 1; k++) {
        if (a[k][k] == 0) {
            int r = -1; for (int i = k + 1; i < n; i++) if (a[i][k]) { r = i; break; }
            if (r < 0) return 0;
            for (int j = 0; j < n; j++) { ll t = a[k][j]; a[k][j] = a[r][j]; a[r][j] = t; }
            sign = -sign;
        }
        for (int i = k + 1; i < n; i++) for (int j = k + 1; j < n; j++)
            a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) / prev;
        prev = a[k][k];
    }
    return sign * a[n - 1][n - 1];
}

int main(int argc, char **argv) {
    int MMAX = argc > 1 ? atoi(argv[1]) : 7;
    memset(seenPrev, 0, sizeof seenPrev);
    seenPrev[1] = 1;   /* S_tau(1) = {1} (the single vertex) */
    for (int m = 2; m <= MMAX; m++) {
        int E = m * (m - 1) / 2;
        int eu[28], ev[28], e = 0;
        for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { eu[e] = i; ev[e] = j; e++; }
        memset(seenCur, 0, sizeof seenCur);
        ll connected = 0, isnew = 0;
        for (uint32_t mask = 0; mask < (1u << E); mask++) {
            /* connectivity by union-find-free BFS on adjacency bitmasks */
            uint32_t adj[8] = {0};
            for (int k = 0; k < E; k++) if (mask >> k & 1) { adj[eu[k]] |= 1u << ev[k]; adj[ev[k]] |= 1u << eu[k]; }
            uint32_t r = 1, prev = 0;
            while (r != prev) { prev = r; for (int v = 0; v < m; v++) if (r >> v & 1) r |= adj[v]; }
            if (r != (1u << m) - 1) continue;
            connected++;
            ll L[8][8];
            for (int i = 1; i < m; i++) for (int j = 1; j < m; j++) {
                if (i == j) L[i - 1][j - 1] = __builtin_popcount(adj[i]);
                else L[i - 1][j - 1] = -(ll)((adj[i] >> j) & 1);
            }
            ll t = (m == 1) ? 1 : det_bareiss(L, m - 1);
            if (t <= 0 || t >= MAXV) { fprintf(stderr, "bad tau %lld\n", t); return 1; }
            seenCur[t] = 1;
            if (!seenPrev[t]) isnew++;
        }
        int nvals = 0, newvals = 0, maxv = 0;
        for (int t = 1; t < MAXV; t++) if (seenCur[t]) { nvals++; maxv = t; if (!seenPrev[t]) newvals++; }
        /* smallest positive integer not in S_tau(m) */
        int gapcount = 0; for (int t = 1; t <= maxv; t++) if (!seenCur[t]) gapcount++;
        printf("m=%d connected labelled=%lld |S_tau(m)|=%d max=%d new values=%d  labelled share NEW=%.4f  missing below max=%d\n",
               m, connected, nvals, maxv, newvals, (double)isnew / connected, gapcount);
        fflush(stdout);
        if (getenv("ALPHA")) { printf("first-at-%d:", m); for (int t = 1; t <= 120; t++) if (seenCur[t] && !seenPrev[t]) printf(" %d", t); printf("\n"); }
        for (int t = 1; t < MAXV; t++) seenPrev[t] |= seenCur[t];
    }
    return 0;
}
