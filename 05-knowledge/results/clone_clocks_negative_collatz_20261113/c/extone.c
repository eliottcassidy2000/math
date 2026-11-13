/* extone.c -- one-vertex extensions of certified two-clock skeletons.
   stdin: lines "n u v row0 .. row(n-1)".  For every attachment pattern of a new vertex z (2^n patterns) test the
   pairs (u,v) and (u,z) for the two-clock shape (complete table with rigorous caps) and print the digits.
   Output: "HIT idx pat which A sigma al sig2 c11" for every two-clock found; c11 = 2 sig2 - sigma.            */
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
int main(int argc, char **argv) {
    int part = argc > 1 ? atoi(argv[1]) : 0, np = argc > 2 ? atoi(argv[2]) : 1;
    int n0, u, v, idx = 0; long long rows[24];
    while (scanf("%d %d %d", &n0, &u, &v) == 3) {
        for (int i = 0; i < n0; i++) scanf("%lld", &rows[i]);
        if (idx++ % np != part) continue;
        n = n0 + 1; int z = n0;
        for (long pat = 0; pat < (1L << n0); pat++) {
            memset(adj, 0, sizeof adj);
            for (int i = 0; i < n0; i++) for (int j = 0; j < n0; j++) adj[i][j] = (rows[i] >> j) & 1;
            for (int j = 0; j < n0; j++) { if ((pat >> j) & 1) adj[z][j] = 1; else adj[j][z] = 1; }
            ll t[32][32];
            for (int which = 0; which < 2; which++) {
                int vv = which ? z : v;
                if (!table(u, vv, t)) continue;
                ll A = t[3][1], al = t[3][2], sig = 3*A - t[1][1], sig2 = 3*al - t[1][2];
                printf("HIT %d %ld %d %lld %lld %lld %lld %lld\n", idx - 1, pat, which, A, sig, al, sig2, 2*sig2 - sig);
                if (2*sig2 - sig == -1) { printf("C11MINUS1 %d %ld %d\n", idx - 1, pat, which); }
            }
        }
        fflush(stdout);
    }
    return 0;
}
