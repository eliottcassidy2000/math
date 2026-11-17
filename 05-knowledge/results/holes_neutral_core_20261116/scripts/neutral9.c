/* neutral9.c: for every class of order 9 (gentourng -q 9, ascii), the pair matrix M and the balanced arcs a -> b
   (M(a, b) = M(b, a), THM-5267 (A)); writes balanced9.txt (class index, a, b, h) and flips9.d6 (the flipped tournament
   in digraph6) for canonical labelling by labelg. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n = 9; static u64 adj[12], radj[12]; static i64 Hm[1<<10], dp[1<<10][10], dpr[1<<10][10], M[10][10];
static void dps(u64 *A, i64 d[][10]) { int N = 1 << n; for (int m = 0; m < N; m++) for (int l = 0; l < n; l++) d[m][l] = 0; for (int v = 0; v < n; v++) d[1<<v][v] = 1;
    for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) { i64 a = d[m][l]; if (!a) continue; u64 nx = A[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; d[m|(1<<v)][v] += a; } } }
static void encode9(u64 *A, char *out) { out[0] = '&'; out[1] = 9 + 63; int p = 2, acc = 0, nb = 0;
    for (int i = 0; i < 9; i++) for (int j = 0; j < 9; j++) { acc = (acc << 1) | (int)((A[i] >> j) & 1); nb++; if (nb == 6) { out[p++] = acc + 63; acc = 0; nb = 0; } }
    if (nb) { acc <<= (6 - nb); out[p++] = acc + 63; } out[p] = 0; }
int main(void) {
    char line[256]; FILE *g = popen("gentourng -q 9", "r"); FILE *fb = fopen("balanced9.txt", "w"), *fd = fopen("flips9.d6", "w");
    int N = 1 << n, full = N - 1; long idx = 0, nbal = 0;
    while (fgets(line, sizeof line, g)) {
        int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue;
        for (int i = 0; i < n; i++) adj[i] = 0; int k = 0;
        for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; }
        for (int v = 0; v < n; v++) { radj[v] = 0; for (int u = 0; u < n; u++) if ((adj[u]>>v)&1) radj[v] |= 1ULL<<u; }
        dps(adj, dp); dps(radj, dpr);
        for (int m = 1; m < N; m++) { i64 s = 0; for (int l = 0; l < n; l++) s += dp[m][l]; Hm[m] = s; } Hm[0] = 1;
        i64 h = Hm[full];
        for (int y = 0; y < n; y++) for (int x = 0; x < n; x++) M[y][x] = 0;
        for (int m = 1; m < full; m++) { int c = full & ~m; for (int y = 0; y < n; y++) { i64 a = dp[m][y]; if (!a) continue; for (int x = 0; x < n; x++) { i64 b = dpr[c][x]; if (b) M[y][x] += a * b; } } }
        for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (((adj[a] >> b) & 1) && M[a][b] == M[b][a]) {
            u64 F[12]; memcpy(F, adj, sizeof F); F[a] &= ~(1ULL << b); F[b] |= 1ULL << a; char s[20]; encode9(F, s);
            fprintf(fb, "%ld %d %d %lld\n", idx, a, b, h); fprintf(fd, "%s\n", s); nbal++;
        }
        idx++;
    }
    pclose(g); fclose(fb); fclose(fd); printf("order 9: %ld classes, %ld balanced arcs\n", idx, nbal); return 0;
}
