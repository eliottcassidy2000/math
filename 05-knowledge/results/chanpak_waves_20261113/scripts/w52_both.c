/* One base, two repairs: a tournament T (one-parent tree, n <= NB) and two one-vertex extensions with counts A and B.
   Usage: both NB A B HLO HHI   (try bases with HLO <= H(T) <= HHI) */
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
static int p[32]; static unsigned out[32]; static ull A, B, LO, HI; static long bases, hitA, hitB, both; static int printed;
static long byh[4096][3];
static void build(int n) { for (int i = 0; i < n; i++) out[i] = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) { if (p[j] == i) out[j] |= 1u << i; else out[i] |= 1u << j; } }
static void try_tree(int n) {
    build(n); ull h = Hc(n, out); if (h < LO || h > HI) return;
    bases++; unsigned base[32]; memcpy(base, out, sizeof(unsigned) * n);
    int a = 0, b = 0; unsigned Sa = 0, Sb = 0;
    for (unsigned S = 0; S < (1u << n); S++) {
        for (int i = 0; i < n; i++) out[i] = base[i] | ((S >> i & 1) ? 0 : (1u << n));
        out[n] = S; ull h2 = Hc(n + 1, out);
        if (h2 == A && !a) { a = 1; Sa = S; }
        if (h2 == B && !b) { b = 1; Sb = S; }
        if (a && b) break;
    }
    byh[h][0]++; if (a) { hitA++; byh[h][1]++; } if (b) { hitB++; byh[h][2]++; }
    if (a && b) { both++; if (printed < 4) { printed++; printf("  BOTH: base H=%llu n=%d parents", h, n); for (int j = 1; j < n; j++) printf(" %d", p[j]);
        printf(" | x beats mask %u -> %llu ; mask %u -> %llu\n", Sa, A, Sb, B); } }
}
static void rec(int j, int n) { if (j == n) { try_tree(n); return; } for (int q = 0; q < j; q++) { p[j] = q; rec(j + 1, n); } }
int main(int argc, char **argv) { int NB = atoi(argv[1]); A = strtoull(argv[2], 0, 10); B = strtoull(argv[3], 0, 10);
    LO = strtoull(argv[4], 0, 10); HI = strtoull(argv[5], 0, 10);
    dp = malloc(sizeof(ull) * ((size_t)1 << (NB + 1)) * (NB + 1));
    for (int n = 1; n <= NB; n++) rec(1, n);
    printf("bases with %llu <= H <= %llu (one-parent, n <= %d): %ld; give %llu: %ld; give %llu: %ld; give both: %ld\n", LO, HI, NB, bases, A, hitA, B, hitB, both);
    printf("  per base count H: (bases, ->%llu, ->%llu):", A, B);
    for (ull h = LO; h <= HI && h < 4096; h++) if (byh[h][0]) printf(" %llu:(%ld,%ld,%ld)", h, byh[h][0], byh[h][1], byh[h][2]);
    printf("\n"); }
