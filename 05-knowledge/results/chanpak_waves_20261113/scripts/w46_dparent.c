/* d-parent tournaments: vertex j beats exactly min(d, j) earlier vertices (any choice) and loses to the rest.
   d = 1 is the one-parent (tree) language.  Report the spectrum completeness. */
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
#define HMAX 2000000
static unsigned char seen[HMAX]; static unsigned back[32]; static int D; static ull hmax;
static void rec(int j, int n) {
    if (j == n) { unsigned out[32] = {0};
        for (int b = 0; b < n; b++) for (int i = 0; i < b; i++) { if (back[b] >> i & 1) out[b] |= 1u << i; else out[i] |= 1u << b; }
        ull h = Hcount(n, out); if (h < HMAX) seen[h] = 1; if (h > hmax) hmax = h; return; }
    int k = j < D ? j : D;
    for (unsigned m = 0; m < (1u << j); m++) if (__builtin_popcount(m) == k) { back[j] = m; rec(j + 1, n); }
}
int main(int argc, char **argv) { D = atoi(argv[1]); int NM = atoi(argv[2]); dp = malloc(sizeof(ull) * ((size_t)1 << NM) * NM);
    for (int n = 1; n <= NM; n++) { hmax = 0; back[0] = 0; rec(1, n); ull bound = 0; int nm = 0;
        for (ull h = 1; h < HMAX; h += 2) if (!seen[h] && h != 7 && h != 21) { bound = h - 2; break; }
        printf("d=%d n<=%d: max H %llu, complete except 7,21 up to %llu; misses:", D, n, hmax, bound);
        for (ull h = 1; h <= hmax && h < HMAX && nm < 10; h += 2) if (!seen[h]) { printf(" %llu", h); nm++; } printf("\n"); fflush(stdout); } }
