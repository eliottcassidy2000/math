/* aw_enum.c - AUDIT D1 (THM-5055): exhaustive enumeration of ALL labelled graphs on n <= NMAX (default 8) vertices
 * whose degree sequence is non-increasing in the labels (every isomorphism class has such a labelling), by Gray code.
 * For each graph: d(G) (subset enumeration), flags (U, #pendants, #pendants with non-universal neighbour, #isolated,
 * has a degree-2 vertex), and - unless its flag tuple is already known to be realised by an even graph - its
 * von Bromssen parity (aw_verify.c stabiliser chain).
 * Output: for each n, every tuple "n d U min(pend,2) min(pnu,2) min(iso,2) r3" realised by an EVEN graph, in the
 * format of mkbase.py (base9f.txt), plus a summary of the even d-values versus F(n).
 * Shares no code with dom_even.c / mkbase.py.
 * usage: aw_enum NMAX [NMIN] [rep] > out   (NMAX <= 9; rep: also print a representative edge list per tuple)
 */
#define AW_NO_MAIN
#include "aw_verify.c"

#define MAXT (1 << 22)
static unsigned char seen_even[10][1 << 9][2][3][3][3][2];   /* [n][d][U][pend][pnu][iso][r3] even-realised */
static unsigned char seen_any[10][1 << 9][2][3][3][3][2];    /* realised by some graph */
static uint64_t rep[10][1 << 9][2][3][3][3][2];              /* edge mask (Gray-code bit order) of a representative */

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 8;
    int NMIN = argc > 2 ? atoi(argv[2]) : 1;
    int REP = argc > 3 && !strcmp(argv[3], "rep");
    if (NMAX > 9) NMAX = 9;
    for (int nn = NMIN; nn <= NMAX; nn++) {
        n = nn;
        int m = n * (n - 1) / 2;
        int eu[64], ev[64], t = 0;
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) { eu[t] = u; ev[t] = v; t++; }
        memset(adj, 0, sizeof adj);
        int deg[MAXN]; memset(deg, 0, sizeof deg);
        uint64_t total = (uint64_t)1 << m, filtered = 0, parity_calls = 0;
        for (uint64_t i = 0; i < total; i++) {
            if (i > 0) {   /* Gray code: flip edge at bit ctz(i) */
                int b = __builtin_ctzll(i);
                int u = eu[b], v = ev[b];
                if (adj[u] >> v & 1) { adj[u] &= ~(1u << v); adj[v] &= ~(1u << u); deg[u]--; deg[v]--; }
                else { adj[u] |= 1u << v; adj[v] |= 1u << u; deg[u]++; deg[v]++; }
            }
            int okd = 1; for (int x = 0; x + 1 < n; x++) if (deg[x] < deg[x + 1]) { okd = 0; break; }
            if (!okd) continue;
            filtered++;
            /* d by subset enumeration */
            uint32_t full = (1u << n) - 1, cl[MAXN]; for (int x = 0; x < n; x++) cl[x] = adj[x] | (1u << x);
            uint32_t dm[1 << 9]; dm[0] = 0; int d = 0;
            for (uint32_t S = 1; S <= full; S++) { dm[S] = dm[S & (S - 1)] | cl[__builtin_ctz(S)]; if (dm[S] == full) d++; }
            int U = 0, pend = 0, pnu = 0, iso = 0, r3 = 0;
            for (int x = 0; x < n; x++) if (deg[x] == n - 1) U = 1;
            for (int x = 0; x < n; x++) {
                if (deg[x] == 1) { pend++; int q = __builtin_ctz(adj[x]); if (deg[q] != n - 1) pnu++; }
                if (deg[x] == 0) iso++;
                if (deg[x] == 2) r3 = 1;
            }
            if (pend > 2) pend = 2; if (pnu > 2) pnu = 2; if (iso > 2) iso = 2;
            seen_any[n][d][U][pend][pnu][iso][r3] = 1;
            if (seen_even[n][d][U][pend][pnu][iso][r3]) continue;
            parity_calls++;
            if (!parity_odd()) { seen_even[n][d][U][pend][pnu][iso][r3] = 1; rep[n][d][U][pend][pnu][iso][r3] = i ^ (i >> 1); }
        }
        fprintf(stderr, "n=%d labelled graphs %llu, degree-sorted %llu, parity calls %llu\n", n,
                (unsigned long long)total, (unsigned long long)filtered, (unsigned long long)parity_calls);
        int maxd = 0; int has[1 << 9]; memset(has, 0, sizeof has);
        for (int d = 1; d < (1 << n); d++)
            for (int U = 0; U < 2; U++) for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) for (int c = 0; c < 3; c++) for (int r = 0; r < 2; r++)
                if (seen_even[n][d][U][a][b][c][r]) {
                    printf("%d %d %d %d %d %d %d", n, d, U, a, b, c, r);
                    if (REP) { uint64_t g = rep[n][d][U][a][b][c][r]; printf(" |"); for (int e = 0; e < m; e++) if (g >> e & 1) printf(" %d-%d", eu[e], ev[e]); }
                    printf("\n"); has[d] = 1; if (d > maxd) maxd = d; }
        int Mn = (n & 1) ? (1 << n) - n : (1 << n) - n - 1;
        int Xn = ((n & 1) && n >= 5) ? (1 << n) - n - 2 : -1;
        int miss = 0, extra = 0;
        for (int k = 1; k < (1 << n); k += 2) {
            int want = k <= Mn && k != 3 && k != 7 && k != 13 && k != Xn;
            if (want && !has[k]) miss++;
            if (!want && has[k]) extra++;
        }
        fprintf(stderr, "n=%d even d-values vs F(n): missing %d extra %d; max even d = %d, M(n) = %d; X(n) %s; 3/7/13 even: %d %d %d\n",
                n, miss, extra, maxd, Mn, Xn < 0 ? "none" : (has[Xn] ? "REALISED" : "not realised"), has[3], has[7], has[13]);
    }
    return 0;
}
