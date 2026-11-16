/* strongB.c <n>: for every STRONG tournament class of order n (gentourng, stdin) and every vertex t: h, e_t, B_t.  Per count
   h: flags digit0 (B = h), digit1any (2^k B = h + 1 for some k >= 0), digit2 (B = h + 2), charge0any (a_t = B - e - h + 1 = 0,
   or 2^k B = h - 1 for some k >= 1), and the sorted set of B values (<= 40 listed). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n; static u64 adj[12]; static i64 Hm[1<<10]; static i64 dp[1<<10][10], dps[1<<10][10];
#define MAXH 4000
static int att[MAXH], d0[MAXH], d1[MAXH], d2[MAXH], c0[MAXH]; static i64 *Bvals[MAXH]; static int nB[MAXH];
static u64 reachfrom(u64 s, u64 *A) { u64 r = s, p = 0; while (r != p) { p = r; for (int v = 0; v < n; v++) if ((r>>v)&1) r |= A[v]; } return r; }
int main(int argc, char **argv) {
    int nn = atoi(argv[1]); char line[256]; long classes = 0, strongc = 0;
    while (fgets(line, sizeof line, stdin)) {
        int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue; classes++;
        n = nn; for (int i = 0; i < n; i++) adj[i] = 0; int k = 0;
        for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; }
        u64 radj[12]; for (int v = 0; v < n; v++) { radj[v] = 0; for (int u = 0; u < n; u++) if ((adj[u]>>v)&1) radj[v] |= 1ULL<<u; }
        u64 full = (1ULL<<n)-1; if (reachfrom(1, adj) != full || reachfrom(1, radj) != full) continue; strongc++;
        int N = 1 << n;
        for (int m = 0; m < N; m++) { Hm[m] = 0; for (int l = 0; l < n; l++) dp[m][l] = 0; }
        for (int v = 0; v < n; v++) dp[1<<v][v] = 1;
        for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) { i64 a = dp[m][l]; if (!a) continue; u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dp[m|(1<<v)][v] += a; } }
        for (int m = 1; m < N; m++) { i64 s = 0; for (int l = 0; l < n; l++) s += dp[m][l]; Hm[m] = s; } Hm[0] = 1;
        i64 h = Hm[N-1]; if (h >= MAXH) continue; att[h] = 1;
        for (int t = 0; t < n; t++) {
            for (int m = 0; m < N; m++) for (int l = 0; l < n; l++) dps[m][l] = 0;
            dps[1<<t][t] = 1;
            for (int m = 1; m < N; m++) { if (!((m>>t)&1)) continue; for (int l = 0; l < n; l++) { i64 a = dps[m][l]; if (!a) continue; u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dps[m|(1<<v)][v] += a; } } }
            i64 e = dp[N-1][t], B = 0;
            for (int m = 1; m < N; m++) if ((m>>t)&1) { i64 st = 0; for (int l = 0; l < n; l++) st += dps[m][l]; B += st * Hm[(N-1) & ~m]; }
            if (B == h) d0[h] = 1; if (B == h + 2) d2[h] = 1;
            for (int kk = 0; kk < 12; kk++) { if ((B << kk) == h + 1) d1[h] = 1; if (kk >= 1 && (B << kk) == h - 1) c0[h] = 1; }
            if (B - e - h + 1 == 0) c0[h] = 1;
            /* record B value */
            int found = 0; for (int i = 0; i < nB[h]; i++) if (Bvals[h][i] == B) { found = 1; break; }
            if (!found) { Bvals[h] = realloc(Bvals[h], sizeof(i64) * (nB[h] + 1)); Bvals[h][nB[h]++] = B; }
        }
    }
    printf("order %d: %ld classes, %ld strong\n", nn, classes, strongc);
    printf("h: digit0 digit1any digit2 charge0any | B values of strong realisers' vertices\n");
    for (int h = 1; h < MAXH; h++) if (att[h]) {
        printf("  %d: %d %d %d %d |", h, d0[h], d1[h], d2[h], c0[h]);
        /* sort */ for (int i = 0; i < nB[h]; i++) for (int j = i+1; j < nB[h]; j++) if (Bvals[h][j] < Bvals[h][i]) { i64 tmp = Bvals[h][i]; Bvals[h][i] = Bvals[h][j]; Bvals[h][j] = tmp; }
        for (int i = 0; i < nB[h] && i < 40; i++) printf(" %lld", Bvals[h][i]); if (nB[h] > 40) printf(" ...(%d)", nB[h]);
        printf("\n");
    }
    return 0;
}
