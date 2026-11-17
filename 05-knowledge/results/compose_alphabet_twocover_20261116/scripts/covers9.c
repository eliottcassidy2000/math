/* covers9.c: h, P(2), P(3) for every class of order 9 (gentourng -q 9), P(3) by submask enumeration. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n = 9; static u64 adj[12]; static i64 Hm[1<<10], dp[1<<10][10];
int main(void) {
    char line[256]; FILE *g = popen("gentourng -q 9", "r"); int N = 1 << n, full = N - 1;
    while (fgets(line, sizeof line, g)) {
        int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue;
        for (int i = 0; i < n; i++) adj[i] = 0; int k = 0;
        for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; }
        for (int m = 0; m < N; m++) for (int l = 0; l < n; l++) dp[m][l] = 0;
        for (int v = 0; v < n; v++) dp[1<<v][v] = 1;
        for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) { i64 a = dp[m][l]; if (!a) continue; u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dp[m|(1<<v)][v] += a; } }
        for (int m = 1; m < N; m++) { i64 s = 0; for (int l = 0; l < n; l++) s += dp[m][l]; Hm[m] = s; } Hm[0] = 1;
        i64 P2 = 0, P3 = 0;
        for (int W = 0; W < N; W++) P2 += Hm[W] * Hm[full & ~W];
        for (int W1 = 0; W1 < N; W1++) { int c = full & ~W1; i64 h1 = Hm[W1]; int sub = c; while (1) { P3 += h1 * Hm[sub] * Hm[c & ~sub]; if (!sub) break; sub = (sub - 1) & c; } }
        printf("%lld %lld %lld\n", Hm[full], P2, P3);
    }
    pclose(g); return 0;
}
