/* repcount.c (opus-2026-11-06-S82, exploratory): representation counts of odd N <= 2^LOGX in the one-apex join families
 *   F1: (1+u)(1+v) + g(u)g(v) - e(u)e(v), 0 <= u <= v            (unordered pairs)
 *   F2: g(v) y_h - e(v) y_e + (1+v) y_R, v >= 0, y in the order-<=7 bank (2,732 triples)
 * r1(N), r2(N) = number of representations (16-bit counters, atomic). Reports per dyadic block the minimum of r1+r2,
 * the number of N with small total counts, the mean, and the hardest N.
 * Usage: repcount LOGX bankfile
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll;
static uint32_t *gtab; static ll X;
static uint16_t *r1, *r2;
static inline void inc(uint16_t *r, ll N) { if (N >= 1 && N <= X && (N & 1)) __atomic_fetch_add(&r[N >> 1], 1, __ATOMIC_RELAXED); }
int main(int argc, char **argv) {
    int LOGX = atoi(argv[1]); X = 1LL << LOGX;
    FILE *f = fopen(argv[2], "r"); int nb = 0; ll yh[4096], ye[4096], yR[4096];
    while (fscanf(f, "%lld %lld %lld", &yh[nb], &ye[nb], &yR[nb]) == 3) nb++;
    fclose(f);
    gtab = malloc(sizeof(uint32_t) * (size_t)(X + 2)); gtab[0] = 1;
    for (ll v = 1; v <= X + 1; v++) { if (v & 1) gtab[v] = gtab[v >> 1]; else { ll u = v >> 1; gtab[v] = (uint32_t)((ll)gtab[u] + u + (u & 1)); } }
    r1 = calloc((size_t)(X / 2 + 1), 2); r2 = calloc((size_t)(X / 2 + 1), 2);
    ll umax = 0; while ((umax + 2) * (umax + 2) <= X + 1) umax++;
    #pragma omp parallel for schedule(dynamic, 1)
    for (ll u = 0; u <= umax; u++) {
        ll gu = gtab[u], eu = (u == 0 || !(u & 1));
        for (ll v = u; (1 + u) * (1 + v) <= X; v++) {
            ll ev = (v == 0 || !(v & 1));
            inc(r1, (1 + u) * (1 + v) + gu * (ll)gtab[v] - eu * ev);
        }
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nb; i++)
        for (ll v = 0; (1 + v) * yR[i] <= X; v++) {
            ll ev = (v == 0 || !(v & 1));
            inc(r2, (ll)gtab[v] * yh[i] - ev * ye[i] + (1 + v) * yR[i]);
        }
    printf("# repcount LOGX=%d bank=%d triples\n", LOGX, nb);
    printf("# k: min(r1+r2) argmin | #N with total<=1, <=2, <=4, <=8 | mean r1, mean r2 | min r2 over F1-misses | #F1-miss\n");
    for (int k = 4; k < LOGX; k++) {
        ll lo = 1LL << k, hi = 1LL << (k + 1), mn = 1LL << 40, arg = -1, c1 = 0, c2 = 0, c4 = 0, c8 = 0, s1 = 0, s2 = 0, cnt = 0, mn2 = 1LL << 40, fm = 0;
        for (ll N = lo + 1; N < hi; N += 2) {
            if (N == 7 || N == 21) continue;
            ll a = r1[N >> 1], b = r2[N >> 1], t = a + b; cnt++; s1 += a; s2 += b;
            if (t < mn) { mn = t; arg = N; }
            if (t <= 1) c1++; if (t <= 2) c2++; if (t <= 4) c4++; if (t <= 8) c8++;
            if (a == 0) { fm++; if (b < mn2) mn2 = b; }
        }
        printf("k=%2d: %lld @%lld | %lld %lld %lld %lld | %.2f %.2f | %lld | %lld\n", k, mn, arg, c1, c2, c4, c8, (double)s1 / cnt, (double)s2 / cnt, fm ? mn2 : -1, fm);
    }
    /* the hardest N overall above 2^12 */
    printf("# N > 4096 with r1+r2 <= 3:\n");
    ll shown = 0;
    for (ll N = 4097; N <= X && shown < 60; N += 2) { ll t = r1[N >> 1] + r2[N >> 1]; if (t <= 3) { printf(" %lld(%d,%d)", N, r1[N >> 1], r2[N >> 1]); shown++; } }
    printf("\n");
    return 0;
}
