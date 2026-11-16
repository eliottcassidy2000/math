/* Charges of tree-reversal tournaments: a_v = H(T + z_v) - 2H(T) + 1, where z_v beats only v.
   For every recursive tree T with n vertices (n <= NM-1) collect the charge set; report statistics:
   how often 0 and 2 are charges, the distribution of min/max charge, and binary closure:
   T is "binary" if some vertex has charge 0 and some vertex has charge 2. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long ull; typedef long long ll;
static ull *dp;
static ull Hcount(int n, const unsigned *out) {
    int N = 1 << n; memset(dp, 0, sizeof(ull) * (size_t)N * n);
    for (int v = 0; v < n; v++) dp[(size_t)(1u << v) * n + v] = 1;
    for (int S = 1; S < N; S++) { ull *row = dp + (size_t)S * n;
        for (int v = 0; v < n; v++) { ull c = row[v]; if (!c) continue; unsigned t = out[v] & ~(unsigned)S;
            while (t) { int u = __builtin_ctz(t); t &= t - 1; dp[(size_t)(S | 1u << u) * n + u] += c; } } }
    ull s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; return s; }
static void build(int n, const int *p, unsigned *out) {
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) { if (p[j] == i) out[j] |= 1u << i; else out[i] |= 1u << j; } }
static int p[32]; static unsigned out[32];
static long ntrees, has0, has2, has02, both_children_02, kids02;
static ll hist[64];
static int NM;
/* returns 1 if tree has charges 0 and 2 */
static int charges(int n, ll *a, ull *H0) {
    build(n, p, out); ull H = Hcount(n, out); *H0 = H;
    for (int v = 0; v < n; v++) { p[n] = v; build(n + 1, p, out); a[v] = (ll)Hcount(n + 1, out) - 2 * (ll)H + 1; }
    int z = 0, t = 0; for (int v = 0; v < n; v++) { if (a[v] == 0) z = 1; if (a[v] == 2) t = 1; }
    return z && t;
}
static void rec(int j, int n) {
    if (j == n) {
        ll a[32]; ull H; int b = charges(n, a, &H); ntrees++;
        int z = 0, t = 0; for (int v = 0; v < n; v++) { if (a[v] == 0) z = 1; if (a[v] == 2) t = 1;
            ll c = a[v]; int bin = c < -30 ? 0 : (c > 30 ? 62 : (int)(c + 31)); hist[bin]++; }
        has0 += z; has2 += t; has02 += b;
        if (b && n + 1 < NM) {   /* do the children at a charge-0 and a charge-2 vertex both stay binary? */
            int v0 = -1, v2 = -1; for (int v = 0; v < n; v++) { if (a[v] == 0 && v0 < 0) v0 = v; if (a[v] == 2 && v2 < 0) v2 = v; }
            int ok = 1; ll a2[32]; ull H2;
            for (int k = 0; k < 2; k++) { p[n] = k ? v2 : v0; if (!charges(n + 1, a2, &H2)) ok = 0; }
            kids02++; both_children_02 += ok;
        }
        return; }
    for (int q = 0; q < j; q++) { p[j] = q; rec(j + 1, n); }
}
int main(int argc, char **argv) {
    NM = atoi(argv[1]); dp = malloc(sizeof(ull) * ((size_t)1 << (NM + 1)) * (NM + 1));
    for (int n = 2; n < NM; n++) {
        ntrees = has0 = has2 = has02 = both_children_02 = kids02 = 0; memset(hist, 0, sizeof hist);
        rec(1, n);
        printf("n=%d trees %ld: charge 0 present %ld, charge 2 present %ld, both %ld; binary trees whose two children stay binary: %ld of %ld\n",
               n, ntrees, has0, has2, has02, both_children_02, kids02);
        printf("   charge histogram (-30..30):"); for (int c = -8; c <= 8; c += 2) printf(" %lld:%lld", (ll)c, hist[c + 31]); printf("  below -30: %lld, above 30: %lld\n", hist[0], hist[62]);
        fflush(stdout);
    }
}
