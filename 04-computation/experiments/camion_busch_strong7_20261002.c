/* Companion to 05-knowledge/results/camion_busch_gaps_polyhedra_collatz_20261002.md.
   All 2^21 labelled 7-tournaments: strong ones, H (Hamiltonian paths), odd-cycle counts c3, c5, c7, disjoint
   3+3 pairs d33; checks H = 1 + 2(c3+c5+c7) + 4 d33 (OCF at n = 7); reports the floor, the alpha_1 minimum, and the
   isomorphism classes of the floor tournaments (canonical form = least arc mask over all 5040 relabellings). */
#include <stdio.h>
#include <string.h>
#define N 7
static int adj[N][N];
static long long hp_count(int n, const int *vs) {           /* Hamiltonian paths of the induced subtournament */
    static long long dp[1 << N][N]; int full = (1 << n) - 1;
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) if (dp[S][v])
        for (int u = 0; u < n; u++) if (!(S >> u & 1) && adj[vs[v]][vs[u]]) dp[S | 1 << u][u] += dp[S][v];
    long long h = 0; for (int v = 0; v < n; v++) h += dp[full][v]; return h;
}
static int perms[5040][N], np = 0;
static void gen(int *a, int k) {
    if (k == N) { for (int i = 0; i < N; i++) perms[np][i] = a[i]; np++; return; }
    for (int i = k; i < N; i++) { int t = a[k]; a[k] = a[i]; a[i] = t; gen(a, k + 1); t = a[k]; a[k] = a[i]; a[i] = t; }
}
static long long hc_count(int n, const int *vs) {           /* directed Hamiltonian cycles of the induced subtournament */
    static long long dp[1 << N][N]; int full = (1 << n) - 1;
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    dp[1][0] = 1;                                            /* cycles through vs[0] as start */
    for (int S = 1; S <= full; S++) if (S & 1) for (int v = 0; v < n; v++) if (dp[S][v])
        for (int u = 1; u < n; u++) if (!(S >> u & 1) && adj[vs[v]][vs[u]]) dp[S | 1 << u][u] += dp[S][v];
    long long c = 0; for (int v = 1; v < n; v++) if (adj[vs[v]][vs[0]]) c += dp[full][v]; return c;
}
int main(void) {
    int pr[21][2], t = 0;
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { pr[t][0] = i; pr[t][1] = j; t++; }
    long long minH = 1 << 30, minA1 = 1 << 30, cntStrong = 0, hist[400]; memset(hist, 0, sizeof hist);
    long long minA1H = 0; int seen21 = 0; long long bad = 0;
    static long long low[32][40]; /* low[H][alpha_1] for strong H <= 31 */
    static long floorMasks[6000]; int nFloor = 0;
    int perm0[N]; for (int i = 0; i < N; i++) perm0[i] = i; gen(perm0, 0);
    for (long mask = 0; mask < (1L << 21); mask++) {
        for (int k = 0; k < 21; k++) { int a = pr[k][0], b = pr[k][1]; int f = (mask >> k) & 1; adj[a][b] = f; adj[b][a] = !f; }
        /* strong connectivity via reachability closure from vertex 0 both ways */
        int fw = 1, bw = 1, ch = 1;
        while (ch) { ch = 0; for (int v = 0; v < N; v++) for (int u = 0; u < N; u++) {
            if ((fw >> v & 1) && adj[v][u] && !(fw >> u & 1)) { fw |= 1 << u; ch = 1; }
            if ((bw >> v & 1) && adj[u][v] && !(bw >> u & 1)) { bw |= 1 << u; ch = 1; } } }
        int all[N] = {0, 1, 2, 3, 4, 5, 6};
        long long H = hp_count(N, all);
        if (H == 21) seen21 = 1;
        if (fw != 127 || bw != 127) continue;
        cntStrong++;
        int c3 = 0; int tri[35][3]; 
        for (int a = 0; a < N; a++) for (int b = a + 1; b < N; b++) for (int c = b + 1; c < N; c++)
            if ((adj[a][b] && adj[b][c] && adj[c][a]) || (adj[a][c] && adj[c][b] && adj[b][a])) { tri[c3][0] = a; tri[c3][1] = b; tri[c3][2] = c; c3++; }
        long long c5 = 0;
        for (int x = 0; x < N; x++) for (int y = x + 1; y < N; y++) {      /* 5-subsets = complement of a pair */
            int vs[5], m = 0; for (int v = 0; v < N; v++) if (v != x && v != y) vs[m++] = v;
            c5 += hc_count(5, vs); }
        long long c7 = hc_count(N, all);
        int d33 = 0;
        for (int i = 0; i < c3; i++) for (int j = i + 1; j < c3; j++) {
            int dis = 1; for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) if (tri[i][p] == tri[j][q]) dis = 0;
            d33 += dis; }
        long long A1 = c3 + c5 + c7;
        if (H != 1 + 2 * A1 + 4 * d33) bad++;
        if (H < minH) minH = H;
        if (A1 < minA1) { minA1 = A1; minA1H = H; }
        if (H < 400) hist[H]++;
        if (H <= 31) low[H][A1]++;
        if (H == 25 && nFloor < 6000) floorMasks[nFloor++] = mask;
    }
    printf("strong 7-tournaments (labelled): %lld; OCF failures: %lld; H = 21 anywhere at n = 7: %d\n", cntStrong, bad, seen21);
    printf("floor f(7) = min H over strong = %lld; min alpha_1 = c3+c5+c7 over strong = %lld (attained with H = %lld)\n", minH, minA1, minA1H);
    for (int h = 25; h <= 31; h += 2) { printf("strong H = %d: (alpha_1, alpha_2 = d33, labelled count):", h);
        for (int a = 0; a < 40; a++) if (low[h][a]) printf(" (%d, %d, %lld)", a, (int)((h - 1 - 2 * a) / 4), low[h][a]);
        printf("\n"); }
    static long canon[6000]; int nc = 0;
    for (int f = 0; f < nFloor; f++) {
        long mask = floorMasks[f];
        for (int k = 0; k < 21; k++) { int a = pr[k][0], b = pr[k][1]; int g = (mask >> k) & 1; adj[a][b] = g; adj[b][a] = !g; }
        long best = -1;
        for (int p = 0; p < np; p++) {
            long m2 = 0;
            for (int k = 0; k < 21; k++) if (adj[perms[p][pr[k][0]]][perms[p][pr[k][1]]]) m2 |= 1L << k;
            if (best < 0 || m2 < best) best = m2;
        }
        int seen = 0; for (int c = 0; c < nc; c++) if (canon[c] == best) seen = 1;
        if (!seen) canon[nc++] = best;
    }
    printf("floor tournaments (strong, H = 25): %d labelled, %d isomorphism class(es)", nFloor, nc);
    if (nc) {
        long m = canon[0]; int sc[N] = {0};
        for (int k = 0; k < 21; k++) { if ((m >> k) & 1) sc[pr[k][0]]++; else sc[pr[k][1]]++; }
        printf("; score vector of the canonical form:"); for (int v = 0; v < N; v++) printf(" %d", sc[v]);
    }
    printf("\n");
    printf("strong H values below 60:"); for (int h = 1; h < 60; h++) if (hist[h]) printf(" %d", h); printf("\n");
    return 0;
}
