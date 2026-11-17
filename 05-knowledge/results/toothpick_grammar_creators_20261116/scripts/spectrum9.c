/* spectrum9.c: for every tournament class of orders 1..9 (gentourng -q on stdin per order), prints one line
   "order index h strong nparents" is NOT done here; this program only prints, per class, "order h" so that the spectrum by
   order, the first creators and their multiplicities can be tabulated (counts.txt).  Orders 1, 2 are built in. */
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
        char cmd[64]; snprintf(cmd, 64, "gentourng -q %d", n); FILE *g = popen(cmd, "r");
        while (fgets(line, sizeof line, g)) {
            int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue;
            for (int i = 0; i < n; i++) adj[i] = 0; int k = 0;
            for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; }
            printf("%d %lld\n", n, hcount());
        }
        pclose(g);
    }
    return 0;
}
