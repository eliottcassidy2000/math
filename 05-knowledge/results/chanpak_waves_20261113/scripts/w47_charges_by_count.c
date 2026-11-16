/* (H, charge) pairs of one-parent tournaments, n <= NM: for each H <= 45 print the set of charges seen. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long ull; typedef long long ll;
static ull *dp;
static ull Hc(int n, const unsigned *out) { int N = 1 << n; memset(dp, 0, sizeof(ull) * (size_t)N * n);
    for (int v = 0; v < n; v++) dp[(size_t)(1u << v) * n + v] = 1;
    for (int S = 1; S < N; S++) { ull *row = dp + (size_t)S * n; for (int v = 0; v < n; v++) { ull c = row[v]; if (!c) continue;
        unsigned t = out[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dp[(size_t)(S | 1u << u) * n + u] += c; } } }
    ull s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; return s; }
static void build(int n, const int *p, unsigned *out) { for (int i = 0; i < n; i++) out[i] = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) { if (p[j] == i) out[j] |= 1u << i; else out[i] |= 1u << j; } }
static int p[32]; static unsigned out[32];
#define HM 46
#define AOFF 200
static unsigned char seen[HM][401];
static void rec(int j, int n) {
    if (j == n) { build(n, p, out); ull H = Hc(n, out); if (H >= HM) return;
        for (int v = 0; v < n; v++) { p[n] = v; build(n + 1, p, out); ll a = (ll)Hc(n + 1, out) - 2 * (ll)H + 1;
            if (a >= -AOFF && a <= AOFF) seen[H][a + AOFF] = 1; } return; }
    for (int q = 0; q < j; q++) { p[j] = q; rec(j + 1, n); } }
int main(int argc, char **argv) { int NM = atoi(argv[1]); dp = malloc(sizeof(ull) * ((size_t)1 << (NM + 1)) * (NM + 1));
    for (int n = 1; n <= NM; n++) rec(1, n);
    for (int H = 1; H < HM; H += 2) { int any = 0; for (int a = 0; a <= 400; a++) any |= seen[H][a]; if (!any) continue;
        printf("H=%2d charges:", H); for (int a = 0; a <= 400; a++) if (seen[H][a]) printf(" %d", a - AOFF);
        printf("   | would give 7 with a=%d, 21 with a=%d\n", 8 - 2 * H, 22 - 2 * H); } }
