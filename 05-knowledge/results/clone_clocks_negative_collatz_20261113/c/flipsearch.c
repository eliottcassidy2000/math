/* flipsearch.c -- local search for the C2-line two-clock (sigma, sig2) = (1, 0).
   stdin: lines "n u v row0 .. row(n-1)".  For each skeleton, every set of at most F reversed arcs (all pairs,
   including arcs at u and v) is applied and the pair (u, v) is tested with a complete two-module table
   (rigorous in-neighbour caps).  Prints "HIT idx f flips... A sigma al sig2" for every two-clock and
   "C2LINE ..." for (sigma, sig2) = (1, 0).   usage: flipsearch F part nparts                          */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int n; static int adj[24][24];
static ll *dp; static size_t dpcap = 0;
static int table(int u, int v, ll t[32][32]) {
    int inu = 0, inv = 0; for (int x = 0; x < n; x++) { if (adj[x][u]) inu++; if (adj[x][v]) inv++; }
    int CU, CV; if (adj[v][u]) { CV = inv + 1; CU = inu + inv + 1; } else { CU = inu + 1; CV = inv + inu + 1; }
    if (CU > 30) CU = 30; if (CV > 30) CV = 30;
    int others[24], m = 0; for (int x = 0; x < n; x++) if (x != u && x != v) others[m++] = x;
    int pos[24]; for (int i = 0; i < m; i++) pos[others[i]] = i;
    size_t S = (size_t)1 << m; size_t W = (size_t)(CU + 1) * (CV + 1) * n;
    if (S * W > dpcap) { free(dp); dpcap = S * W; dp = malloc(dpcap * sizeof(ll)); }
    memset(dp, 0, S * W * sizeof(ll));
    #define I(mask,cu,cv,last) ((size_t)(mask) * W + ((size_t)(cu) * (CV + 1) + (cv)) * n + (last))
    dp[I(0,1,0,u)] = 1; dp[I(0,0,1,v)] = 1; for (int i = 0; i < m; i++) dp[I((size_t)1 << i,0,0,others[i])] = 1;
    memset(t, 0, sizeof(ll) * 32 * 32);
    for (size_t mask = 0; mask < S; mask++) for (int cu = 0; cu <= CU; cu++) for (int cv = 0; cv <= CV; cv++) for (int last = 0; last < n; last++) {
        ll val = dp[I(mask,cu,cv,last)]; if (!val) continue;
        if (mask == S - 1) t[cu][cv] += val;
        for (int w = 0; w < n; w++) { if (!adj[last][w]) continue;
            if (w == u) { if (last == u || cu == CU) continue; dp[I(mask,cu+1,cv,u)] += val; }
            else if (w == v) { if (last == v || cv == CV) continue; dp[I(mask,cu,cv+1,v)] += val; }
            else { size_t b = (size_t)1 << pos[w]; if (mask & b) continue; dp[I(mask|b,cu,cv,w)] += val; } }
    }
    for (int i = 1; i < 32; i++) for (int ip = 1; ip < 32; ip++) if ((i >= 4 || ip >= 3) && t[i][ip]) return 0;
    if (t[3][1] <= 0 || t[3][2] <= 0) return 0;
    if (t[2][1] != 3 * t[3][1] || t[2][2] != 3 * t[3][2]) return 0;
    return 1;
}
static int u0, v0, idxg; static long hits = 0, c2 = 0;
static int pa[300], pb[300], np;
static void flip(int k) { int a = pa[k], b = pb[k]; int t_ = adj[a][b]; adj[a][b] = adj[b][a]; adj[b][a] = t_; }
static void test(int f, int *ks) {
    ll t[32][32];
    if (!table(u0, v0, t)) return;
    ll A = t[3][1], al = t[3][2], sig = 3*A - t[1][1], sig2 = 3*al - t[1][2];
    hits++;
    printf("HIT %d %d", idxg, f); for (int i = 0; i < f; i++) printf(" %d-%d", pa[ks[i]], pb[ks[i]]);
    printf(" A=%lld sigma=%lld al=%lld sig2=%lld\n", A, sig, al, sig2);
    if (sig == 1 && sig2 == 0) { c2++; printf("C2LINE %d\n", idxg); }
}
int main(int argc, char **argv) {
    int F = atoi(argv[1]), part = argc > 2 ? atoi(argv[2]) : 0, nparts = argc > 3 ? atoi(argv[3]) : 1;
    long long rows[24]; int idx = 0;
    while (scanf("%d %d %d", &n, &u0, &v0) == 3) {
        for (int i = 0; i < n; i++) scanf("%lld", &rows[i]);
        idxg = idx++; if (idxg % nparts != part) continue;
        memset(adj, 0, sizeof adj);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj[i][j] = (rows[i] >> j) & 1;
        np = 0; for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) { pa[np] = a; pb[np] = b; np++; }
        int ks[4];
        for (int i = 0; i < np; i++) { flip(i); ks[0] = i; test(1, ks);
            if (F >= 2) for (int j = i + 1; j < np; j++) { flip(j); ks[1] = j; test(2, ks);
                if (F >= 3) for (int k = j + 1; k < np; k++) { flip(k); ks[2] = k; test(3, ks); flip(k); }
                flip(j); }
            flip(i); }
        fflush(stdout);
    }
    fprintf(stderr, "part %d: hits %ld c2line %ld\n", part, hits, c2);
    return 0;
}
