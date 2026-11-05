/* minhc2conn.c (klein-2026-11-05): over all tournaments of order n (gentourng format on stdin) that are 2-connected
   (every one-vertex deletion strong), the histogram of the number of Hamiltonian cycles and its minimum. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static u64 adj[16]; static u64 *dp;
static int strong_sub(u64 verts) {
    int s = __builtin_ctzll(verts); u64 seen = 1ULL << s, fr = seen;
    while (fr) { u64 nf = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nf |= adj[v] & verts; nf &= ~seen; seen |= nf; fr = nf; }
    if (seen != verts) return 0;
    seen = 1ULL << s; fr = seen;
    while (fr) { u64 nf = 0; for (int v = 0; v < n; v++) if ((fr >> v & 1)) for (int u = 0; u < n; u++) if ((verts >> u & 1) && (adj[u] >> v & 1)) nf |= 1ULL << u; nf &= ~seen; seen |= nf; fr = nf; }
    return seen == verts;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); u64 N = 1ULL << n; dp = malloc(sizeof(u64) * N * n); char line[256];
    u64 total = 0, twoconn = 0; static u64 hist[4096]; u64 minhc = ~0ULL;
    while (fgets(line, sizeof line, stdin)) { int L = strlen(line); while (L && (line[L-1]=='\n'||line[L-1]=='\r')) line[--L] = 0; if (L != n*(n-1)/2) continue;
        memset(adj, 0, sizeof adj); int idx = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
        total++;
        u64 all = N - 1; int ok = 1;
        for (int v = 0; v < n && ok; v++) if (!strong_sub(all & ~(1ULL << v))) ok = 0;
        if (!ok) continue;
        twoconn++;
        memset(dp, 0, sizeof(u64) * N * n); dp[1 * n + 0] = 1;
        for (u64 S = 1; S < N; S += 2) for (int v = 0; v < n; v++) { u64 c = dp[S * n + v]; if (!c) continue; u64 nx = adj[v] & ~S; while (nx) { int u = __builtin_ctzll(nx); nx &= nx - 1; dp[(S | (1ULL << u)) * n + u] += c; } }
        u64 hc = 0; for (int v = 1; v < n; v++) if (adj[v] & 1ULL) hc += dp[(N - 1) * n + v];
        if (hc < 4096) hist[hc]++; if (hc < minhc) minhc = hc;
    }
    printf("order %d: %llu tournaments, %llu 2-connected, minimum hc among 2-connected = %llu; histogram of hc (2-connected, hc < 40):", n, total, twoconn, minhc);
    for (int k = 0; k < 40; k++) if (hist[k]) printf(" %d:%llu", k, hist[k]); printf("\n"); return 0; }
