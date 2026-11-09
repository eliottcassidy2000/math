/* aw_dump.c - AUDIT D1: dump every labelled graph on NMIN..NMAX (<= 9) vertices with non-increasing degree sequence
 * (so every isomorphism class appears at least once) whose number of dominating sets lies in a target set.
 * Used for an independent census of the small d-values in Prop 3.2 (3, 7, 13), Remark 7.1 (11, 23, 47) and
 * HYP-9625 (5, 37).  Completeness: a connected graph on n >= 2 vertices has complementary dominating sets A, B (Ore),
 * so d >= 2^|A| + 2^|B| - 1 >= 2^(n/2 + 1) - 1; hence every connected realiser of d <= 62 has <= 9 vertices.
 * usage: aw_dump NMAX NMIN d1 d2 ... > out      (output: "n d u-v u-v ...")
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
int main(int argc, char **argv) {
    int NMAX = atoi(argv[1]), NMIN = atoi(argv[2]);
    int want[1024]; memset(want, 0, sizeof want);
    for (int i = 3; i < argc; i++) want[atoi(argv[i])] = 1;
    for (int n = NMIN; n <= NMAX; n++) {
        int m = n * (n - 1) / 2, eu[64], ev[64], t = 0;
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) { eu[t] = u; ev[t] = v; t++; }
        uint32_t adj[16]; int deg[16]; memset(adj, 0, sizeof adj); memset(deg, 0, sizeof deg);
        uint64_t total = (uint64_t)1 << m, hits = 0;
        for (uint64_t i = 0; i < total; i++) {
            if (i > 0) {
                int b = __builtin_ctzll(i), u = eu[b], v = ev[b];
                if (adj[u] >> v & 1) { adj[u] &= ~(1u << v); adj[v] &= ~(1u << u); deg[u]--; deg[v]--; }
                else { adj[u] |= 1u << v; adj[v] |= 1u << u; deg[u]++; deg[v]++; }
            }
            int ok = 1; for (int x = 0; x + 1 < n; x++) if (deg[x] < deg[x + 1]) { ok = 0; break; }
            if (!ok) continue;
            uint32_t full = (1u << n) - 1, cl[16], dm[512]; int d = 0;
            for (int x = 0; x < n; x++) cl[x] = adj[x] | (1u << x);
            dm[0] = 0;
            for (uint32_t S = 1; S <= full; S++) { dm[S] = dm[S & (S - 1)] | cl[__builtin_ctz(S)]; if (dm[S] == full) d++; }
            if (d < 1024 && want[d]) {
                hits++;
                printf("%d %d", n, d);
                for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (adj[u] >> v & 1) printf(" %d-%d", u, v);
                printf("\n");
            }
        }
        fprintf(stderr, "n=%d hits=%llu\n", n, (unsigned long long)hits);
    }
    return 0;
}
