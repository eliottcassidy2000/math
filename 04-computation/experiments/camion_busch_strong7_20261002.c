/* Companion to 05-knowledge/results/camion_busch_gaps_polyhedra_collatz_20261002.md.
   All 2^21 labelled 7-tournaments: strong ones, H (Hamiltonian paths), odd-cycle counts c3, c5, c7, disjoint
   3+3 pairs d33; checks H = 1 + 2(c3+c5+c7) + 4 d33 (OCF at n = 7). Reports the floor f(7), the minimum number of
   odd cycles alpha_1 and the alpha_2 values that occur with it, and the isomorphism classes of the floor tournaments
   and of the alpha_1-minimal ones (canonical form = least arc mask over all 5040 relabellings), compared with Moon's
   T'_7 (p_i -> p_j iff i = j - 1 or i >= j + 2). */
#include <stdio.h>
#include <string.h>
#define N 7
static int adj[N][N], pr[21][2];
static int perms[5040][N], np = 0;
static void gen(int *a, int k) {
    if (k == N) { for (int i = 0; i < N; i++) perms[np][i] = a[i]; np++; return; }
    for (int i = k; i < N; i++) { int t = a[k]; a[k] = a[i]; a[i] = t; gen(a, k + 1); t = a[k]; a[k] = a[i]; a[i] = t; }
}
static void load(long mask) {
    for (int k = 0; k < 21; k++) { int a = pr[k][0], b = pr[k][1]; int f = (mask >> k) & 1; adj[a][b] = f; adj[b][a] = !f; }
}
static long canon_of(long mask) {
    load(mask);
    long best = -1;
    for (int p = 0; p < np; p++) {
        long m2 = 0;
        for (int k = 0; k < 21; k++) if (adj[perms[p][pr[k][0]]][perms[p][pr[k][1]]]) m2 |= 1L << k;
        if (best < 0 || m2 < best) best = m2;
    }
    return best;
}
static int classes(const long *masks, int n, long *out) {
    int nc = 0;
    for (int f = 0; f < n; f++) {
        long c = canon_of(masks[f]); int seen = 0;
        for (int i = 0; i < nc; i++) if (out[i] == c) seen = 1;
        if (!seen) out[nc++] = c;
    }
    return nc;
}
static long long hp_count(int n, const int *vs) {           /* Hamiltonian paths of the induced subtournament */
    static long long dp[1 << N][N]; int full = (1 << n) - 1;
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) if (dp[S][v])
        for (int u = 0; u < n; u++) if (!(S >> u & 1) && adj[vs[v]][vs[u]]) dp[S | 1 << u][u] += dp[S][v];
    long long h = 0; for (int v = 0; v < n; v++) h += dp[full][v]; return h;
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
    int t = 0;
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { pr[t][0] = i; pr[t][1] = j; t++; }
    int perm0[N]; for (int i = 0; i < N; i++) perm0[i] = i; gen(perm0, 0);
    long long minH = 1 << 30, minA1 = 1 << 30, cntStrong = 0, hist[400], bad = 0; memset(hist, 0, sizeof hist);
    int seen21 = 0;
    static long long low[32][40];        /* low[H][alpha_1] for strong H <= 31 */
    static long long a1d[40][20];        /* a1d[alpha_1][alpha_2] over all strong */
    static long floorMasks[6000], a1Masks[6000]; int nFloor = 0, nA1 = 0;
    for (long mask = 0; mask < (1L << 21); mask++) {
        load(mask);
        int fw = 1, bw = 1, ch = 1;      /* strong connectivity via reachability from vertex 0 both ways */
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
        if (A1 < minA1) minA1 = A1;
        if (H < 400) hist[H]++;
        if (H <= 31) low[H][A1]++;
        if (A1 < 40 && d33 < 20) a1d[A1][d33]++;
        if (H == 25 && nFloor < 6000) floorMasks[nFloor++] = mask;
        if (A1 == 9 && nA1 < 6000) a1Masks[nA1++] = mask;
    }
    printf("strong 7-tournaments (labelled): %lld; OCF failures: %lld; H = 21 anywhere at n = 7: %d\n", cntStrong, bad, seen21);
    printf("floor f(7) = min H over strong = %lld; min alpha_1 = c3+c5+c7 over strong = %lld\n", minH, minA1);
    printf("alpha_2 values occurring with alpha_1 = %lld:", minA1);
    for (int d = 0; d < 20; d++) if (a1d[minA1][d]) printf(" alpha_2 = %d (H = %lld, %lld labelled)", d, 1 + 2 * minA1 + 4 * d, a1d[minA1][d]);
    printf("\n");
    for (int h = 25; h <= 31; h += 2) { printf("strong H = %d: (alpha_1, alpha_2 = d33, labelled count):", h);
        for (int a = 0; a < 40; a++) if (low[h][a]) printf(" (%d, %d, %lld)", a, (int)((h - 1 - 2 * a) / 4), low[h][a]);
        printf("\n"); }
    static long canonF[6000], canonA[6000];
    int ncF = classes(floorMasks, nFloor, canonF);
    printf("floor tournaments (strong, H = 25): %d labelled, %d isomorphism class(es)", nFloor, ncF);
    if (ncF) {
        long m = canonF[0]; int sc[N] = {0};
        for (int k = 0; k < 21; k++) { if ((m >> k) & 1) sc[pr[k][0]]++; else sc[pr[k][1]]++; }
        printf("; score vector of the canonical form:"); for (int v = 0; v < N; v++) printf(" %d", sc[v]);
    }
    printf("\n");
    int ncA = (minA1 == 9) ? classes(a1Masks, nA1, canonA) : -1;
    long tprime = 0;                     /* Moon's T'_7: a -> b (a < b) iff b = a + 1 */
    for (int k = 0; k < 21; k++) if (pr[k][1] == pr[k][0] + 1) tprime |= 1L << k;
    long ct = canon_of(tprime);
    printf("alpha_1-minimal tournaments (alpha_1 = 9): %d labelled, %d isomorphism class(es); Moon's T'_7 is that class: %s\n",
           nA1, ncA, (ncA == 1 && canonA[0] == ct) ? "yes" : "no");
    printf("strong H values below 60:"); for (int h = 1; h < 60; h++) if (hist[h]) printf(" %d", h); printf("\n");
    return 0;
}
