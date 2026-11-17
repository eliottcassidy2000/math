/* creators9.c (R-5, R-6): first creators (classes whose count first appears at their order) for orders 3..9.  Reads
   counts.txt (order h per class, gentourng order) to build the new-count sets per order, then re-reads the gentourng
   streams: for every first creator F of order n and every vertex v, H(F - v) and the out-degree of v; a genealogy edge is
   (F, v) with H(F - v) new at order n - 1 (so F - v is a first creator of its count, THM-4922).  Reports per order: new
   counts, first creators, those with a genealogy edge, with a sole-exit/sole-entry genealogy edge (d+(v) in {1, n-2}), with
   a digit-0 vertex (B_t = h), the out-degree distribution of genealogy vertices, and the orphans. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXH 4096
static unsigned char seen[MAXH], newc[10][MAXH];
static int n; static u64 adj[12], radj[12]; static i64 Hm[1<<10], dp[1<<10][10], dpr[1<<10][10];
static void dps(u64 *A, i64 d[][10], int nn) {
    int N = 1 << nn; for (int m = 0; m < N; m++) for (int l = 0; l < nn; l++) d[m][l] = 0;
    for (int v = 0; v < nn; v++) d[1<<v][v] = 1;
    for (int m = 1; m < N; m++) for (int l = 0; l < nn; l++) { i64 a = d[m][l]; if (!a) continue; u64 nx = A[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; d[m|(1<<v)][v] += a; } }
}
int main(void) {
    FILE *c = fopen("counts.txt", "r"); int o; long h;
    while (fscanf(c, "%d %ld", &o, &h) == 2) { if (h < MAXH && !seen[h]) { seen[h] = 1; newc[o][h] = 1; } }
    fclose(c);
    printf("order: new counts, first creators, with a genealogy edge, with a sole-exit/sole-entry genealogy edge, with a digit-0 vertex | out-degree distribution of genealogy vertices | orphans (count)\n");
    char line[256];
    for (n = 3; n <= 9; n++) {
        long nnew = 0; for (int x = 0; x < MAXH; x++) nnew += newc[n][x];
        long fc = 0, gen = 0, sole = 0, d0 = 0; long deg[12] = {0}; long orph[64]; int norph = 0;
        char cmd[64]; snprintf(cmd, 64, "gentourng -q %d", n); FILE *g = popen(cmd, "r");
        while (fgets(line, sizeof line, g)) {
            int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue;
            for (int i = 0; i < n; i++) adj[i] = 0; int k = 0;
            for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; }
            int N = 1 << n; dps(adj, dp, n);
            for (int m = 1; m < N; m++) { i64 s = 0; for (int l = 0; l < n; l++) s += dp[m][l]; Hm[m] = s; } Hm[0] = 1;
            i64 hh = Hm[N-1]; if (hh >= MAXH || !newc[n][hh]) continue;
            fc++;
            for (int v = 0; v < n; v++) { radj[v] = 0; for (int u = 0; u < n; u++) if ((adj[u]>>v)&1) radj[v] |= 1ULL<<u; }
            dps(radj, dpr, n);
            int isd0 = 0; for (int t = 0; t < n && !isd0; t++) { i64 B = 0; for (int m = 1; m < N; m++) if ((m>>t)&1) B += dpr[m][t] * Hm[(N-1) & ~m]; if (B == hh) isd0 = 1; }
            d0 += isd0;
            int any = 0, anysole = 0;
            for (int v = 0; v < n; v++) {
                i64 hv = Hm[(N-1) & ~(1 << v)];            /* count of F - v: paths covering V minus v */
                if (hv < MAXH && newc[n-1][hv]) { any = 1; int d = __builtin_popcountll(adj[v]); deg[d]++; if (d == 1 || d == n - 2) anysole = 1; }
            }
            gen += any; sole += anysole; if (!any && norph < 64) orph[norph++] = hh;
        }
        pclose(g);
        printf("  %d: %ld, %ld, %ld, %ld, %ld |", n, nnew, fc, gen, sole, d0);
        for (int d = 0; d < n; d++) if (deg[d]) printf(" %d:%ld", d, deg[d]);
        printf(" | orphans:"); for (int i = 0; i < norph; i++) printf(" %ld", orph[i]); printf("\n"); fflush(stdout);
    }
    return 0;
}
