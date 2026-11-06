/*
 * chanpak_join_multiplicity_20261106.c  (opus-2026-11-06-S82)
 *
 * Representation counts of odd N <= 2^LOGX in the one-arc join families of THM-4875 (D):
 *   r1(N) = #{0 <= u <= v : (1+u)(1+v) + g(u)g(v) - e(u)e(v) = N}
 *   r2(N) = #{(y, v) : y in the bank, v >= 0, y_h g(v) - y_e e(v) + y_R (1+v) = N}
 * (16-bit atomic counters; the largest counts are far below 65,535 at LOGX <= 30).
 * Per dyadic block [2^k, 2^(k+1)): mean r1, increment, mean r2, min r2, # of r2 = 0, lower-tail fractions,
 * mean r2 by residue mod 3; and the mass constant M_2 = sum_y (2/y_h) ln(1 + y_h/y_R) of the bank.
 * Usage: multiplicity LOGX bankfile
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
typedef long long ll;
static uint32_t *gtab; static ll X; static uint16_t *r1, *r2;
static inline void inc(uint16_t *r, ll N) { if (N >= 1 && N <= X && (N & 1)) __atomic_fetch_add(&r[N >> 1], 1, __ATOMIC_RELAXED); }
int main(int argc, char **argv) {
    int LOGX = atoi(argv[1]); X = 1LL << LOGX;
    FILE *f = fopen(argv[2], "r"); static ll yh[4096], ye[4096], yR[4096]; int nb = 0;
    while (nb < 4096 && fscanf(f, "%lld %lld %lld", &yh[nb], &ye[nb], &yR[nb]) == 3) nb++;
    fclose(f);
    double M2 = 0; for (int i = 0; i < nb; i++) M2 += (2.0 / yh[i]) * log(1.0 + (double)yh[i] / yR[i]);
    printf("# chanpak_join_multiplicity_20261106 LOGX=%d bank=%d triples; M_2 = %.6f; (pi^2/12) ln 2 = %.6f\n", LOGX, nb, M2, M_PI * M_PI / 12 * log(2.0));
    gtab = malloc(sizeof(uint32_t) * (size_t)(X + 2)); gtab[0] = 1;   /* g(v) <= v + 1 < 2^32 for v <= 2^30 */
    for (ll v = 1; v <= X + 1; v++) { if (v & 1) gtab[v] = gtab[v >> 1]; else { ll u = v >> 1; gtab[v] = (uint32_t)((ll)gtab[u] + u + (u & 1)); } }
    r1 = calloc((size_t)(X / 2 + 1), 2); r2 = calloc((size_t)(X / 2 + 1), 2);
    ll umax = 0; while ((umax + 2) * (umax + 2) <= X + 1) umax++;
    #pragma omp parallel for schedule(dynamic, 1)
    for (ll u = 0; u <= umax; u++) {
        ll gu = gtab[u], eu = (u == 0 || !(u & 1));
        for (ll v = u; (1 + u) * (1 + v) <= X; v++) { ll ev = (v == 0 || !(v & 1)); inc(r1, (1 + u) * (1 + v) + gu * (ll)gtab[v] - eu * ev); }
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nb; i++)
        for (ll v = 0; (1 + v) * yR[i] <= X; v++) { ll ev = (v == 0 || !(v & 1)); inc(r2, (ll)gtab[v] * yh[i] - ev * ye[i] + (1 + v) * yR[i]); }
    printf("# k | mean r1 | incr | mean r2 | min r2 | #r2=0 | frac r2<=10, <=12, <=20 | mean r2 by N mod 3 (0,1,2)\n");
    double prev = -1;
    for (int k = 10; k < LOGX; k++) {
        ll lo = 1LL << k, hi = 1LL << (k + 1), n = 0, z = 0, mn = 1 << 30, c10 = 0, c12 = 0, c20 = 0; double s1 = 0, s2 = 0, sm[3] = {0, 0, 0}; ll cm[3] = {0, 0, 0};
        for (ll N = lo + 1; N < hi; N += 2) {
            ll a = r1[N >> 1], b = r2[N >> 1]; n++; s1 += a; s2 += b;
            if (b == 0) z++; if (b < mn) mn = b; if (b <= 10) c10++; if (b <= 12) c12++; if (b <= 20) c20++;
            sm[N % 3] += b; cm[N % 3]++;
        }
        double m1 = s1 / n;
        printf("k=%2d | %.4f | %+.4f | %.4f | %lld | %lld | %.3e %.3e %.3e | %.3f %.3f %.3f\n", k, m1, prev < 0 ? 0.0 : m1 - prev, s2 / n, mn, z,
               (double)c10 / n, (double)c12 / n, (double)c20 / n, sm[0] / cm[0], sm[1] / cm[1], sm[2] / cm[2]);
        prev = m1;
    }
    ll miss = 0; printf("odd N <= 2^%d with r2 = 0:", LOGX);
    for (ll N = 1; N <= X; N += 2) if (!r2[N >> 1]) { miss++; if (miss <= 50) printf(" %lld", N); }
    printf("  [total %lld]\n", miss);
    return 0;
}
