/* hval.c (klein-2026-11-05): for each "hc digraph6" line on stdin print "hc h" where h = number of Hamiltonian paths (subset DP) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int parse_d6(const char *s, u64 *a) { if (s[0] != '&') return -1; int n = s[1] - 63; int bits = n * n; const char *p = s + 2;
    for (int i = 0; i < n; i++) a[i] = 0;
    for (int k = 0; k < bits; k++) { int ch = p[k / 6] - 63; int b = (ch >> (5 - k % 6)) & 1; if (b) a[k / n] |= 1ULL << (k % n); } return n; }
int main(void) { char line[512]; u64 adj[32]; u64 *dp = malloc(sizeof(u64) * (1ULL << 17) * 17);
    while (fgets(line, sizeof line, stdin)) { int h; char d6[256]; if (sscanf(line, "%d %255s", &h, d6) != 2) continue;
        int n = parse_d6(d6, adj); u64 N = 1ULL << n; memset(dp, 0, sizeof(u64) * N * n);
        for (int v = 0; v < n; v++) dp[(1ULL << v) * n + v] = 1;
        for (u64 S = 1; S < N; S++) for (int v = 0; v < n; v++) { u64 c = dp[S * n + v]; if (!c) continue; u64 nx = adj[v] & ~S;
            while (nx) { int u = __builtin_ctzll(nx); nx &= nx - 1; dp[(S | (1ULL << u)) * n + u] += c; } }
        u64 tot = 0; for (int v = 0; v < n; v++) tot += dp[(N - 1) * n + v];
        printf("%d %llu\n", h, tot); }
    return 0; }
