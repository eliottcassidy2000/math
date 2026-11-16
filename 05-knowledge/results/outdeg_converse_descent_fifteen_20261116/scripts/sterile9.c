/* sterile9.c <n>: reads gentourng -q output for order n on stdin, keeps the STRONG classes, and for each computes h = H and
   the number of digit-0 vertices (root response B_t = h; B_t = sum over masks m containing t of (paths from t covering m)
   times H(complement)).  Reports the strong floor (least h), the number of sterile classes (no digit-0 vertex), the largest
   number of digit-0 vertices in one class, the counts h < 200 carried by a sterile strong class, and the counts h < 200
   carried only by fertile strong classes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n; static u64 adj[12], radj[12]; static i64 Hm[1<<10], dp[1<<10][10], dpr[1<<10][10];
#define MAXH 100000
static unsigned char ster[MAXH], fert[MAXH];
static u64 reachfrom(u64 s, u64 *A) { u64 r = s, p = 0; while (r != p) { p = r; for (int v = 0; v < n; v++) if ((r>>v)&1) r |= A[v]; } return r; }
int main(int argc, char **argv) {
    n = atoi(argv[1]); char line[256]; long classes = 0, strongc = 0, sterile = 0; i64 floor_h = -1; int maxd0 = 0; i64 maxd0_h = 0;
    while (fgets(line, sizeof line, stdin)) {
        int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue; classes++;
        for (int i = 0; i < n; i++) adj[i] = 0; int k = 0;
        for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; }
        for (int v = 0; v < n; v++) { radj[v] = 0; for (int u = 0; u < n; u++) if ((adj[u]>>v)&1) radj[v] |= 1ULL<<u; }
        u64 full = (1ULL<<n)-1; if (reachfrom(1, adj) != full || reachfrom(1, radj) != full) continue; strongc++;
        int N = 1 << n;
        for (int m = 0; m < N; m++) for (int l = 0; l < n; l++) { dp[m][l] = 0; dpr[m][l] = 0; }
        for (int v = 0; v < n; v++) { dp[1<<v][v] = 1; dpr[1<<v][v] = 1; }
        for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) {
            i64 a = dp[m][l]; if (a) { u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dp[m|(1<<v)][v] += a; } }
            i64 b = dpr[m][l]; if (b) { u64 nx = radj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dpr[m|(1<<v)][v] += b; } }
        }
        for (int m = 1; m < N; m++) { i64 s = 0; for (int l = 0; l < n; l++) s += dp[m][l]; Hm[m] = s; } Hm[0] = 1;
        i64 h = Hm[N-1]; int d0 = 0;
        for (int t = 0; t < n; t++) { i64 B = 0; for (int m = 1; m < N; m++) if ((m>>t)&1) B += dpr[m][t] * Hm[(N-1) & ~m]; if (B == h) d0++; }
        if (floor_h < 0 || h < floor_h) floor_h = h;
        if (d0 > maxd0) { maxd0 = d0; maxd0_h = h; }
        if (h < MAXH) { if (d0) fert[h] = 1; else ster[h] = 1; }
        if (!d0) sterile++;
    }
    printf("order %d: %ld classes, %ld strong, %ld sterile; strong floor %lld; max digit-0 vertices in one strong class %d (at h = %lld)\n", n, classes, strongc, sterile, floor_h, maxd0, maxd0_h);
    printf("  counts < 200 carried by a sterile strong class:"); for (int h = 1; h < 200; h++) if (ster[h]) printf(" %d", h); printf("\n");
    printf("  counts < 200 carried only by fertile strong classes:"); for (int h = 1; h < 200; h++) if (fert[h] && !ster[h]) printf(" %d", h); printf("\n");
    return 0;
}
