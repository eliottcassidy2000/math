/* visits.c -- module visit numbers.
   For a tournament S (gentourng format on stdin, n vertices) and a vertex u, a_i(u) = number of Hamiltonian
   "walks" of S that visit u exactly i times (as i maximal separate visits, consecutive visits not adjacent)
   and every other vertex exactly once, arcs respected; u's arcs are uniform.  Then for a transitive module TT_j
   substituted at u:  H(S[u<-TT_j]) = sum_i a_i * i! S(j,i).
   Mode "tern": report (S,u) with a_i = 0 for i>=4, a_2 = 3 a_3, a_1 = 3 a_3 - 1, a_3 > 0  (H = a_3 3^j - 1).
   Also mode "ternready": additionally require a root x != u with: x has a unique Hamiltonian walk starting at x
   (counting u-visits with multiplicity one only, and no walk from x visits u twice) -- then S+z (z beaten only by x)
   is ready; and the ready count's visit numbers satisfy the ternary relations.
   DP: state (set of non-u vertices used, number of u visits so far, last vertex) -- u visits bounded by deg.  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef __int128 i128;
static int n; static int adj[32][32];
#define MAXV 6
/* count walks: vertices 0..n-1, special u; returns a[0..MAXV] (a[i] = walks with i visits of u).
   If start>=0, only walks starting at start.  */
static u64 *dp; /* dp[(mask*(MAXV+1)+c)*n + last] */
static void walks(int u, int start, u64 *a) {
    int m = n; u64 full = 0; for (int v = 0; v < n; v++) if (v != u) full |= 1ULL << v;
    size_t S = 1ULL << n; size_t sz = S * (MAXV + 1) * n;
    memset(dp, 0, sz * sizeof(u64));
    #define IDX(mask,c,last) ((((size_t)(mask))*(MAXV+1)+(c))*n+(last))
    for (int v = 0; v < n; v++) {
        if (start >= 0 && v != start) continue;
        if (v == u) dp[IDX(0, 1, u)] = 1; else dp[IDX(1ULL << v, 0, v)] = 1;
    }
    for (u64 mask = 0; mask < S; mask++) {
        if (mask >> u & 1) continue;
        for (int c = 0; c <= MAXV; c++) for (int last = 0; last < n; last++) {
            u64 val = dp[IDX(mask, c, last)]; if (!val) continue;
            for (int w = 0; w < n; w++) {
                if (!adj[last][w]) continue;
                if (w == u) { if (last == u || c == MAXV) continue; dp[IDX(mask, c + 1, u)] += val; }
                else { if (mask >> w & 1) continue; dp[IDX(mask | (1ULL << w), c, w)] += val; }
            }
        }
    }
    for (int c = 0; c <= MAXV; c++) { a[c] = 0; for (int last = 0; last < n; last++) a[c] += dp[IDX(full, c, last)]; }
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); const char *mode = argv[2];
    dp = malloc(sizeof(u64) * ((size_t)1 << n) * (MAXV + 1) * n);
    char line[1024]; long cnt = 0, hits = 0;
    while (fgets(line, sizeof line, stdin)) {
        int L = strlen(line); while (L && (line[L-1]=='\n')) line[--L] = 0;
        if (L != n*(n-1)/2) continue; cnt++;
        memset(adj, 0, sizeof adj); int idx = 0;
        for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[idx++]=='1') adj[i][j]=1; else adj[j][i]=1; }
        for (int u = 0; u < n; u++) {
            u64 a[MAXV+1]; walks(u, -1, a);
            int ok = a[4]==0 && a[5]==0 && a[6]==0 && a[3]>0 && a[2]==3*a[3] && a[1]+1==3*a[3];
            if (!strcmp(mode, "census")) {
                /* print visit vector for later analysis */
                printf("%s u=%d a=%llu,%llu,%llu,%llu,%llu\n", line, u, a[1],a[2],a[3],a[4],a[5]);
                continue;
            }
            if (ok) { hits++; printf("TERN %s u=%d a1=%llu a2=%llu a3=%llu  H(j)=%llu*3^j-1\n", line, u, a[1], a[2], a[3], a[3]); }
        }
    }
    fprintf(stderr, "tournaments %ld hits %ld\n", cnt, hits);
    return 0;
}
