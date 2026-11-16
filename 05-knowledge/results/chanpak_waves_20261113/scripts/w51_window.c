/* Matched windows as binary port pairs: for a base count h, look for a tournament T with H(T) = h (one-parent trees,
   n <= NB) and two one-vertex extensions with counts 2h-1 and 2h+1 (the window centred at 2h).  Usage: window NB H */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long ull;
static ull *dp;
static ull Hc(int n, const unsigned *out) { int N = 1 << n; memset(dp, 0, sizeof(ull) * (size_t)N * n);
    for (int v = 0; v < n; v++) dp[(size_t)(1u << v) * n + v] = 1;
    for (int S = 1; S < N; S++) { ull *row = dp + (size_t)S * n; for (int v = 0; v < n; v++) { ull c = row[v]; if (!c) continue;
        unsigned t = out[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dp[(size_t)(S | 1u << u) * n + u] += c; } } }
    ull s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; return s; }
static int p[32]; static unsigned out[32]; static ull TARGET; static long bases, lowonly, highonly, both; static int printed;
static void build(int n) { for (int i = 0; i < n; i++) out[i] = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) { if (p[j] == i) out[j] |= 1u << i; else out[i] |= 1u << j; } }
static void try_tree(int n) {
    build(n); if (Hc(n, out) != TARGET) return;
    bases++; unsigned base[32]; memcpy(base, out, sizeof(unsigned) * n);
    int lo = 0, hi = 0; unsigned Slo = 0, Shi = 0;
    for (unsigned S = 0; S < (1u << n); S++) {
        for (int i = 0; i < n; i++) out[i] = base[i] | ((S >> i & 1) ? 0 : (1u << n));
        out[n] = S; ull h2 = Hc(n + 1, out);
        if (h2 == 2 * TARGET - 1 && !lo) { lo = 1; Slo = S; }
        if (h2 == 2 * TARGET + 1 && !hi) { hi = 1; Shi = S; }
    }
    if (lo && hi) { both++; if (printed < 3) { printed++; printf("  base n=%d parents", n); for (int j = 1; j < n; j++) printf(" %d", p[j]);
        printf(" | new vertex beats mask %u -> %llu, mask %u -> %llu\n", Slo, 2 * TARGET - 1, Shi, 2 * TARGET + 1); } }
    else if (lo) lowonly++; else if (hi) { highonly++; if (printed < 2) { printed++; printf("  127-type witness: base n=%d parents", n); for (int j = 1; j < n; j++) printf(" %d", p[j]); printf(" | new vertex beats mask %u -> %llu\n", Shi, 2 * TARGET + 1); } }
}
static void rec(int j, int n) { if (j == n) { try_tree(n); return; } for (int q = 0; q < j; q++) { p[j] = q; rec(j + 1, n); } }
int main(int argc, char **argv) { int NB = atoi(argv[1]); TARGET = strtoull(argv[2], 0, 10);
    dp = malloc(sizeof(ull) * ((size_t)1 << (NB + 1)) * (NB + 1));
    for (int n = 1; n <= NB; n++) rec(1, n);
    printf("base count %llu (window %llu/%llu, centre %llu): one-parent bases n<=%d: %ld; both repairs: %ld; only %llu: %ld; only %llu: %ld\n",
           TARGET, 2 * TARGET - 1, 2 * TARGET + 1, 2 * TARGET, NB, bases, both, 2 * TARGET - 1, lowonly, 2 * TARGET + 1, highonly); }
