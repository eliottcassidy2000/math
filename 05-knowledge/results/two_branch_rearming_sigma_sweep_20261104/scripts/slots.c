/* slots.c (klein-2026-11-04) -- for each tournament on stdin (gentourng default format): h = [P], K = [P^2], K3 = [P^3]
   in the square-zero vertex ring, i.e. h(S) for every subset S, then K = sum_S h(S) h(V\S), K3 = sum over ordered
   triples of disjoint subsets covering V of h(S1) h(S2) h(S3) (empty subsets allowed, h(empty) = 1). Prints "h K K3 code". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static u64 adj[16]; static u64 *f; static u64 *hs;
int main(int argc, char **argv) {
    n = atoi(argv[1]); u64 N = 1ULL << n; f = malloc(sizeof(u64) * N * n); hs = malloc(sizeof(u64) * N);
    char line[256];
    while (fgets(line, sizeof line, stdin)) {
        int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
        if (L != n * (n - 1) / 2) continue;
        memset(adj, 0, sizeof adj); int idx = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
        memset(f, 0, sizeof(u64) * N * n);
        for (int v = 0; v < n; v++) f[(1ULL << v) * n + v] = 1;
        for (u64 S = 1; S < N; S++) { u64 tot = 0; for (int v = 0; v < n; v++) { u64 c = f[S * n + v]; if (!c) continue; tot += c; u64 out = adj[v] & ~S; while (out) { int w = __builtin_ctzll(out); out &= out - 1; f[(S | (1ULL << w)) * n + w] += c; } } hs[S] = tot; }
        hs[0] = 1;
        u64 full = N - 1, K = 0, K3 = 0;
        for (u64 S = 0; S < N; S++) K += hs[S] * hs[full ^ S];
        for (u64 S1 = 0; S1 < N; S1++) { u64 rest = full ^ S1; u64 acc = 0; for (u64 S2 = rest; ; S2 = (S2 - 1) & rest) { acc += hs[S2] * hs[rest ^ S2]; if (S2 == 0) break; } K3 += hs[S1] * acc; }
        printf("%llu %llu %llu %s\n", hs[full], K, K3, line);
    }
    return 0;
}
