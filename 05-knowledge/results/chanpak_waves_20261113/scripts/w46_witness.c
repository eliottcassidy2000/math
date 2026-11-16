/* print one recursive tree (parent list) realizing each target H, searching n = 1..NM */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long ull;
static ull *dp;
static ull Hcount(int n, const unsigned *out) { int N = 1 << n; memset(dp, 0, sizeof(ull) * (size_t)N * n);
    for (int v = 0; v < n; v++) dp[(size_t)(1u << v) * n + v] = 1;
    for (int S = 1; S < N; S++) { ull *row = dp + (size_t)S * n; for (int v = 0; v < n; v++) { ull c = row[v]; if (!c) continue;
        unsigned t = out[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dp[(size_t)(S | 1u << u) * n + u] += c; } } }
    ull s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; return s; }
static int p[32]; static unsigned out[32]; static ull *targets; static int nt, *found; static int fp[64][32], fn[64];
static void rec(int j, int n) { if (j == n) { for (int i = 0; i < n; i++) out[i] = 0;
        for (int b = 1; b < n; b++) for (int i = 0; i < b; i++) { if (p[b] == i) out[b] |= 1u << i; else out[i] |= 1u << b; }
        ull h = Hcount(n, out); for (int k = 0; k < nt; k++) if (!found[k] && h == targets[k]) { found[k] = 1; fn[k] = n; memcpy(fp[k], p, sizeof p); } return; }
    for (int q = 0; q < j; q++) { p[j] = q; rec(j + 1, n); } }
int main(int argc, char **argv) { int NM = atoi(argv[1]); nt = argc - 2; targets = malloc(sizeof(ull) * nt); found = calloc(nt, sizeof(int));
    for (int k = 0; k < nt; k++) targets[k] = strtoull(argv[k + 2], 0, 10);
    dp = malloc(sizeof(ull) * ((size_t)1 << NM) * NM);
    for (int n = 1; n <= NM; n++) rec(1, n);
    for (int k = 0; k < nt; k++) { printf("%llu:", targets[k]); if (found[k]) { printf(" n=%d parents", fn[k]); for (int j = 1; j < fn[k]; j++) printf(" %d", fp[k][j]); } else printf(" none"); printf("\n"); } }
