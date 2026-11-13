/* certify.c -- complete two-module visit table with rigorous caps.
   stdin: lines "n u v row0 row1 ... row(n-1)" (rows = out-neighbour bitmasks).
   Caps: if v->u then i' <= |N^-(v)|+1 and i <= |N^-(u)|+|N^-(v)|+1 (each internal excursion of a symbol ends at a
   distinct in-neighbour vertex or at a symbol of the other module); symmetric if u->v.  The DP never needs to stop
   early below these caps, so the table is complete.  Output: the nonzero entries a[i][i'] with i,i' >= 1.        */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
int main(void) {
    int n, u, v;
    while (scanf("%d %d %d", &n, &u, &v) == 3) {
        int adj[24][24] = {{0}}; long long rows[24];
        for (int i = 0; i < n; i++) { scanf("%lld", &rows[i]); for (int j = 0; j < n; j++) adj[i][j] = (rows[i] >> j) & 1; }
        int inu = 0, inv = 0; for (int x = 0; x < n; x++) { if (adj[x][u]) inu++; if (adj[x][v]) inv++; }
        int CU, CV;
        if (adj[v][u]) { CV = inv + 1; CU = inu + inv + 1; } else { CU = inu + 1; CV = inv + inu + 1; }
        int others[24], m = 0; for (int x = 0; x < n; x++) if (x != u && x != v) others[m++] = x;
        int pos[24]; for (int i = 0; i < m; i++) pos[others[i]] = i;
        size_t S = (size_t)1 << m; size_t W = (size_t)(CU + 1) * (CV + 1) * n;
        ll *dp = calloc(S * W, sizeof(ll));
        #define I(mask,cu,cv,last) ((size_t)(mask) * W + ((size_t)(cu) * (CV + 1) + (cv)) * n + (last))
        dp[I(0,1,0,u)] = 1; dp[I(0,0,1,v)] = 1; for (int i = 0; i < m; i++) dp[I((size_t)1 << i,0,0,others[i])] = 1;
        ll t[32][32]; memset(t, 0, sizeof t);
        for (size_t mask = 0; mask < S; mask++) for (int cu = 0; cu <= CU; cu++) for (int cv = 0; cv <= CV; cv++) for (int last = 0; last < n; last++) {
            ll val = dp[I(mask,cu,cv,last)]; if (!val) continue;
            if (mask == S - 1) t[cu][cv] += val;
            for (int w = 0; w < n; w++) { if (!adj[last][w]) continue;
                if (w == u) { if (last == u) continue; if (cu == CU) { fprintf(stderr, "CAP EXCEEDED\n"); return 2; } dp[I(mask,cu+1,cv,u)] += val; }
                else if (w == v) { if (last == v) continue; if (cv == CV) { fprintf(stderr, "CAP EXCEEDED\n"); return 2; } dp[I(mask,cu,cv+1,v)] += val; }
                else { size_t b = (size_t)1 << pos[w]; if (mask & b) continue; dp[I(mask|b,cu,cv,w)] += val; } }
        }
        free(dp);
        printf("TABLE");
        for (int i = 1; i <= CU; i++) for (int ip = 1; ip <= CV; ip++) if (t[i][ip]) printf(" %d,%d:%lld", i, ip, t[i][ip]);
        printf("\n"); fflush(stdout);
    }
    return 0;
}
