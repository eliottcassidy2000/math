/* nj10.c <L>: near-join values h + 2^k B (k >= 1) from every rooted class of orders 1..10 (gentourng -q n on stdin per
   order via popen), bitmap up to L; prints the first odd miss above 67 and the misses per dyadic window; writes nj10.bits. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n; static u64 adj[12], radj[12]; static i64 Hm[1<<11], dp[1<<11][11], dpr[1<<11][11];
static unsigned char *nj; static long L;
static void setbit(long v) { nj[v >> 3] |= 1 << (v & 7); }
static void process(void) {
    for (int v = 0; v < n; v++) { radj[v] = 0; for (int u = 0; u < n; u++) if ((adj[u]>>v)&1) radj[v] |= 1ULL<<u; }
    int N = 1 << n;
    for (int m = 0; m < N; m++) for (int l = 0; l < n; l++) { dp[m][l] = 0; dpr[m][l] = 0; }
    for (int v = 0; v < n; v++) { dp[1<<v][v] = 1; dpr[1<<v][v] = 1; }
    for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) {
        i64 a = dp[m][l]; if (a) { u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dp[m|(1<<v)][v] += a; } }
        i64 b = dpr[m][l]; if (b) { u64 nx = radj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dpr[m|(1<<v)][v] += b; } }
    }
    for (int m = 1; m < N; m++) { i64 s = 0; for (int l = 0; l < n; l++) s += dp[m][l]; Hm[m] = s; } Hm[0] = 1;
    i64 h = Hm[N-1];
    for (int t = 0; t < n; t++) {
        i64 B = 0; for (int m = 1; m < N; m++) if ((m>>t)&1) B += dpr[m][t] * Hm[(N-1) & ~m];
        for (int k = 1; k < 40; k++) { i64 v = h + (B << k); if (v > L) break; setbit(v); }
    }
}
int main(int argc, char **argv) {
    L = atol(argv[1]); nj = calloc(L/8 + 2, 1); char line[256]; long classes = 0;
    n = 1; adj[0] = 0; process(); n = 2; adj[0] = 2; adj[1] = 0; process();
    for (n = 3; n <= 10; n++) {
        char cmd[64]; snprintf(cmd, 64, "gentourng -q %d", n); FILE *g = popen(cmd, "r");
        while (fgets(line, sizeof line, g)) {
            int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue; classes++;
            for (int i = 0; i < n; i++) adj[i] = 0; int k = 0;
            for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; }
            process();
            if (n == 10 && classes % 1000000 == 0) fprintf(stderr, "order 10: %ld classes so far\n", classes);
        }
        pclose(g); fprintf(stderr, "order %d done (%ld classes cumulative)\n", n, classes);
    }
    FILE *o = fopen("nj10.bits", "wb"); fwrite(nj, 1, L/8 + 2, o); fclose(o);
    long first = -1, miss = 0; printf("classes of orders 3-10: %ld; L = %ld\n", classes, L);
    for (long x = 69; x <= L; x += 2) if (!(nj[x>>3] & (1 << (x&7)))) { miss++; if (first < 0) first = x; }
    printf("near-join (k >= 1), bases <= 10: every odd <= %ld except the seven is realised; first miss %ld; odd misses in (67, %ld]: %ld\n", first - 2, first, L, miss);
    printf("misses by dyadic window [2^j, 2^(j+1)):"); for (int j = 18; (1L << j) < L; j++) { long lo = 1L << j, hi = (1L << (j+1)) < L ? (1L << (j+1)) : L, m = 0, tot = 0; for (long x = lo | 1; x < hi; x += 2) { tot++; if (!(nj[x>>3] & (1 << (x&7)))) m++; } printf(" j=%d: %ld/%ld", j, m, tot); } printf("\n");
    return 0;
}
