/* r2resid.c (opus-2026-11-06-S82, exploratory): F2 multiplicity r2(N) by residue class and its lower tail.
 * Usage: r2resid LOGX bankfile
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll;
static uint32_t *gtab; static ll X; static uint16_t *r2;
int main(int argc, char **argv) {
    int LOGX = atoi(argv[1]); X = 1LL << LOGX;
    FILE *f = fopen(argv[2], "r"); int nb = 0; static ll yh[4096], ye[4096], yR[4096];
    while (fscanf(f, "%lld %lld %lld", &yh[nb], &ye[nb], &yR[nb]) == 3) nb++;
    fclose(f);
    gtab = malloc(sizeof(uint32_t) * (size_t)(X + 2)); gtab[0] = 1;
    for (ll v = 1; v <= X + 1; v++) { if (v & 1) gtab[v] = gtab[v >> 1]; else { ll u = v >> 1; gtab[v] = (uint32_t)((ll)gtab[u] + u + (u & 1)); } }
    r2 = calloc((size_t)(X / 2 + 1), 2);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nb; i++)
        for (ll v = 0; (1 + v) * yR[i] <= X; v++) {
            ll ev = (v == 0 || !(v & 1)); ll N = (ll)gtab[v] * yh[i] - ev * ye[i] + (1 + v) * yR[i];
            if (N >= 1 && N <= X && (N & 1)) __atomic_fetch_add(&r2[N >> 1], 1, __ATOMIC_RELAXED);
        }
    int K = LOGX - 1; ll lo = 1LL << K, hi = 1LL << (K + 1);
    printf("# r2resid LOGX=%d, block [2^%d, 2^%d)\n", LOGX, K, K + 1);
    int mods[] = {3, 9, 27, 8, 16, 5, 7};
    for (int mi = 0; mi < 7; mi++) {
        int m = mods[mi]; double s[64] = {0}; ll c[64] = {0}; ll mn[64]; for (int r = 0; r < 64; r++) mn[r] = 1 << 30;
        for (ll N = lo + 1; N < hi; N += 2) { int r = (int)(N % m); s[r] += r2[N >> 1]; c[r]++; if (r2[N >> 1] < mn[r]) mn[r] = r2[N >> 1]; }
        printf("mod %2d (mean/min):", m);
        for (int r = 0; r < m; r++) if (c[r]) printf(" %d:%.2f/%lld", r, s[r] / c[r], mn[r]);
        printf("\n");
    }
    /* lower tail histogram */
    ll hist[64] = {0}; for (ll N = lo + 1; N < hi; N += 2) { int t = r2[N >> 1]; if (t < 64) hist[t]++; }
    printf("lower tail: #N with r2 = t, t = 0..30:"); for (int t = 0; t <= 30; t++) printf(" %lld", hist[t]); printf("\n");
    printf("hardest N (r2 <= 7):"); int sh = 0; for (ll N = lo + 1; N < hi && sh < 40; N += 2) if (r2[N >> 1] <= 7) { printf(" %lld(%d)", N, r2[N >> 1]); sh++; } printf("\n");
    return 0;
}
