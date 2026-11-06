/*
 * aud82_kappa.c  (independent audit of opus-2026-11-06-S82)
 *
 * Theorem B constant.  Count(X) = #{0 <= u <= v : F1(u,v) <= X} = (pi^2/24) X ln X + kappa X + o(X) with
 *   kappa = kappa1 - c',
 *   kappa1 = lim_U [ sum_{u=0}^{U-1} phi(t_u)/(1+u) - (pi^2/12) ln U ],  phi(t) = ln(1+t)/t,  t_u = g(u)/(1+u),
 *   c'     = int_0^1 dlam int_0^lam dmu P(lam mu (1 + r s) <= 1)  (r, s iid uniform)   [the pairs v < u],
 * so the block mean of r1 is (pi^2/12)(k+2) ln 2 + 2 kappa.  This program prints the partial sums
 *   D(m) = sum_{u < 2^m} phi(t_u)/(1+u) - (pi^2/12) ln 2^m
 * and the block deviations B(m) = sum_{2^m <= u < 2^(m+1)} phi(t_u)/(1+u) - (pi^2/12) ln 2 for m <= MMAX.
 * c' is computed separately (aud82_kappa.py).
 * Usage: aud82_kappa MMAX
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
typedef long long ll;
#define SB 16
#define SN (1 << SB)
static ll gs[SN], Ct[SN], Et[SN], Kt[SN];
static ll g_closed(ll v) {
    if (v == 0) return 1;
    int top = 63 - __builtin_clzll((unsigned long long)v);
    ll s = 1;
    for (int i = 0; i < top; i++) if (!((v >> i) & 1)) s += (v >> (i + 1)) + ((v >> (i + 1)) & 1);
    return s;
}
static inline ll G(ll v) {
    if (v < SN) return gs[v];
    ll k = v >> SB, r = v & (SN - 1);
    return G(k) + k * Ct[r] + (k & 1) * Et[r] + Kt[r];
}
int main(int argc, char **argv) {
    int MMAX = argc > 1 ? atoi(argv[1]) : 30;
    for (ll v = 0; v < SN; v++) gs[v] = g_closed(v);
    for (int r = 0; r < SN; r++) {
        ll C = 0, K = 0;
        for (int i = 0; i < SB; i++) if (!((r >> i) & 1)) { C += 1LL << (SB - 1 - i); K += r >> (i + 1); }
        for (int i = 0; i + 1 < SB; i++) if (((r >> (i + 1)) & 1) && !((r >> i) & 1)) K++;
        Ct[r] = C; Et[r] = !((r >> (SB - 1)) & 1); Kt[r] = K;
    }
    { ll bad = 0; for (ll v = 0; v < (1LL << 22); v++) if (G(v) != g_closed(v)) bad++;
      for (ll v = (1LL << 40) - 5000; v < (1LL << 40) + 5000; v++) if (G(v) != g_closed(v)) bad++;
      printf("# split law check: %lld failures\n", bad); if (bad) return 1; }
    const long double c0 = (long double)M_PI * M_PI / 12.0L;
    long double S = 1.0L * logl(2.0L);   /* u = 0: t = g(0)/1 = 1, phi(1) = ln 2 */
    printf("# m  D(m) = sum_{u<2^m} phi/(1+u) - (pi^2/12) m ln2      B(m-1) block deviation\n");
    for (int m = 1; m <= MMAX; m++) {
        /* block u in [2^(m-1), 2^m) */
        ll lo = 1LL << (m - 1), hi = 1LL << m;
        long double bs = 0;
        #pragma omp parallel for reduction(+:bs) schedule(static)
        for (ll u = lo; u < hi; u++) {
            long double t = (long double)G(u) / (long double)(1 + u);
            bs += log1pl(t) / t / (long double)(1 + u);
        }
        S += bs;
        printf("m=%2d  D=%.12Lf   B(%d)=%+.3Le\n", m, S - c0 * m * logl(2.0L), m - 1, bs - c0 * logl(2.0L));
        fflush(stdout);
    }
    return 0;
}
