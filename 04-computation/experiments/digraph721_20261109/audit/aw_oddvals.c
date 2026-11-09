/* aw_oddvals.c - AUDIT D1 (HYP-9625 evidence): for n <= 9, which odd d < 512 are realised by an ODD graph
 * (von Bromssen) on exactly n vertices, by exhaustive enumeration of labelled graphs with non-increasing degree
 * sequence (every isomorphism class appears).  Parity by aw_verify.c's stabiliser chain; one parity call per graph
 * whose d is not yet known to be odd-realised.  Prints the odd numbers < 2^n NOT realised by odd graphs on n vertices.
 * usage: aw_oddvals NMAX */
#define AW_NO_MAIN
#include "aw_verify.c"
static unsigned char odd_real[10][512];
int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 8; if (NMAX > 9) NMAX = 9;
    for (int nn = 1; nn <= NMAX; nn++) {
        n = nn; int m = n * (n - 1) / 2, eu[64], ev[64], t = 0;
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) { eu[t] = u; ev[t] = v; t++; }
        memset(adj, 0, sizeof adj); int deg[MAXN]; memset(deg, 0, sizeof deg);
        uint64_t total = (uint64_t)1 << m, calls = 0;
        for (uint64_t i = 0; i < total; i++) {
            if (i > 0) {
                int b = __builtin_ctzll(i), u = eu[b], v = ev[b];
                if (adj[u] >> v & 1) { adj[u] &= ~(1u << v); adj[v] &= ~(1u << u); deg[u]--; deg[v]--; }
                else { adj[u] |= 1u << v; adj[v] |= 1u << u; deg[u]++; deg[v]++; }
            }
            int okd = 1; for (int x = 0; x + 1 < n; x++) if (deg[x] < deg[x + 1]) { okd = 0; break; }
            if (!okd) continue;
            uint32_t full = (1u << n) - 1, cl[MAXN], dm[512]; int d = 0;
            for (int x = 0; x < n; x++) cl[x] = adj[x] | (1u << x);
            dm[0] = 0;
            for (uint32_t S = 1; S <= full; S++) { dm[S] = dm[S & (S - 1)] | cl[__builtin_ctz(S)]; if (dm[S] == full) d++; }
            if (odd_real[n][d]) continue;
            calls++;
            if (parity_odd()) odd_real[n][d] = 1;
        }
        printf("n=%d (parity calls %llu): odd k < 2^n not realised by an odd graph on n vertices:", n, (unsigned long long)calls);
        for (int k = 1; k < (1 << n); k += 2) if (!odd_real[n][k]) printf(" %d", k);
        printf("\n"); fflush(stdout);
    }
    return 0;
}
