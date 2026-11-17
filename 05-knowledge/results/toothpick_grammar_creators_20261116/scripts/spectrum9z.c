/* spectrum9z.c: counts h for every class of orders 3..9 read from gentourng -q -z (digraph6), to align with the
   canonical-parent files; orders 1, 2 built in.  Output "order h" per class (counts_z.txt), to be compared with counts.txt. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n; static u64 adj[12]; static i64 dp[1<<10][10];
static i64 hcount(void) {
    int N = 1 << n; for (int m = 0; m < N; m++) for (int l = 0; l < n; l++) dp[m][l] = 0;
    for (int v = 0; v < n; v++) dp[1<<v][v] = 1;
    for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) { i64 a = dp[m][l]; if (!a) continue; u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dp[m|(1<<v)][v] += a; } }
    i64 h = 0; for (int l = 0; l < n; l++) h += dp[N-1][l]; return h;
}
int main(void) {
    char line[256];
    n = 1; adj[0] = 0; printf("1 %lld\n", hcount());
    n = 2; adj[0] = 2; adj[1] = 0; printf("2 %lld\n", hcount());
    for (n = 3; n <= 9; n++) {
        char cmd[64]; snprintf(cmd, 64, "gentourng -q -z %d", n); FILE *g = popen(cmd, "r");
        while (fgets(line, sizeof line, g)) {
            int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (len < 2 || line[0] != '&') continue;
            int nn = line[1] - 63; if (nn != n) { fprintf(stderr, "order mismatch\n"); return 1; }
            for (int i = 0; i < n; i++) adj[i] = 0;
            int bit = 0;
            for (int p = 2; p < len; p++) { int c = line[p] - 63; for (int b = 5; b >= 0; b--, bit++) { if (bit >= n * n) break; if ((c >> b) & 1) { int i = bit / n, j = bit % n; adj[i] |= 1ULL << j; } } }
            printf("%d %lld\n", n, hcount());
        }
        pclose(g);
    }
    return 0;
}
