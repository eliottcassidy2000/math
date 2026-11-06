/*
 * aud82_sieve.c  (independent audit of opus-2026-11-06-S82; no code shared with the session)
 *
 * Representation counts of the odd N <= 2^LOGX:
 *   r1(N) = #{0 <= u <= v : (1+u)(1+v) + g(u) g(v) - e(u) e(v) = N}
 *   r2(N) = #{(y, v) : y in the bank, v >= 0 : y_h g(v) - y_e e(v) + y_R (1+v) = N}
 * g from the auditor's closed form g(v) = 1 + sum_{zero bits i < top} (floor(v/2^(i+1)) + bit_(i+1)(v)) on v < 2^15,
 * and beyond by the auditor's split law (derived from the closed form), for v = k 2^15 + r, k >= 1:
 *   g(v) = g(k) + k C(r) + [k odd] eps(r) + K(r),
 *   C(r) = sum_{zero bits i < 15 of r} 2^(14-i),  eps(r) = [bit 14 of r = 0],
 *   K(r) = #{i <= 13 : bit_(i+1)(r) = 1, bit_i(r) = 0} + sum_{zero bits i < 15 of r} floor(r / 2^(i+1)).
 * The split law is checked against the recursion on all v < 2^24 and on 10^7 random v <= 2^LOGX plus the endpoints.
 * 16-bit counters (atomic); a saturation counter reports any overflow.  Every produced value is checked to be odd.
 * Usage: aud82_sieve LOGX bankfile [doF1=1] [statk=25]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
typedef long long ll;
#define SB 15
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
static ll grec_f(ll v) { if (v == 0) return 1; if (v & 1) return grec_f(v >> 1); ll u = v >> 1; return grec_f(u) + u + (u & 1); }
static uint16_t *r1, *r2; static ll X; static ll nsat = 0, neven = 0;
static inline void bump(uint16_t *r, ll N) {
    if (N < 1 || N > X) return;
    if (!(N & 1)) { __atomic_fetch_add(&neven, 1, __ATOMIC_RELAXED); return; }
    uint16_t old = __atomic_fetch_add(&r[N >> 1], 1, __ATOMIC_RELAXED);
    if (old == 65535) __atomic_fetch_add(&nsat, 1, __ATOMIC_RELAXED);
}
int main(int argc, char **argv) {
    int LOGX = atoi(argv[1]); X = 1LL << LOGX;
    int doF1 = argc > 3 ? atoi(argv[3]) : 1; int statk = argc > 4 ? atoi(argv[4]) : 25;
    FILE *f = fopen(argv[2], "r"); if (!f) { printf("no bank\n"); return 1; }
    int cap = 1 << 17; ll *yh = malloc(sizeof(ll) * cap), *ye = malloc(sizeof(ll) * cap), *yR = malloc(sizeof(ll) * cap); int nb = 0;
    while (nb < cap && fscanf(f, "%lld %lld %lld", &yh[nb], &ye[nb], &yR[nb]) == 3) nb++;
    fclose(f);
    long double M = 0; for (int i = 0; i < nb; i++) M += (2.0L / yh[i]) * logl(1.0L + (long double)yh[i] / yR[i]);
    printf("# aud82_sieve LOGX=%d bank=%s (%d triples), M = %.6Lf, threads=%d\n", LOGX, argv[2], nb, M, omp_get_max_threads());
    for (ll v = 0; v < SN; v++) gs[v] = g_closed(v);
    for (int r = 0; r < SN; r++) {
        ll C = 0, K = 0;
        for (int i = 0; i < SB; i++) if (!((r >> i) & 1)) { C += 1LL << (SB - 1 - i); K += r >> (i + 1); }
        for (int i = 0; i + 1 < SB; i++) if (((r >> (i + 1)) & 1) && !((r >> i) & 1)) K++;
        Ct[r] = C; Et[r] = !((r >> (SB - 1)) & 1); Kt[r] = K;
    }
    {   ll bad = 0, VM = 1LL << 24;
        ll *t = malloc(sizeof(ll) * (VM + 1)); t[0] = 1;
        for (ll v = 1; v <= VM; v++) t[v] = (v & 1) ? t[v >> 1] : t[v >> 1] + (v >> 1) + ((v >> 1) & 1);
        for (ll v = 0; v <= VM; v++) if (t[v] != G(v)) bad++;
        free(t);
        unsigned long long s = 0x2545F4914F6CDD1DULL;
        for (int i = 0; i < 10000000; i++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; ll v = (ll)(s % (unsigned long long)(X + 2)); if (grec_f(v) != G(v)) bad++; }
        for (ll v = X - 3; v <= X + 1; v++) if (grec_f(v) != G(v)) bad++;
        printf("split law vs recursion (all v <= 2^24, 10^7 random v <= 2^LOGX + 1, endpoints): %lld failures\n", bad);
        if (bad) return 2;
    }
    r1 = calloc((size_t)(X / 2 + 1), 2); r2 = calloc((size_t)(X / 2 + 1), 2);
    if (!r1 || !r2) { printf("alloc failed\n"); return 3; }
    double t0 = omp_get_wtime();
    ll npairs = 0;
    if (doF1) {
        /* work items: (u, v-chunk) */
        ll umax = 0; while ((umax + 1) * (umax + 1) <= X) umax++; umax--;   /* (1+u)^2 <= X */
        const ll CH = 1 << 18;
        ll nitems = 0; ll *iu = NULL, *iv0 = NULL; ll capi = 0;
        for (ll u = 0; u <= umax; u++) {
            ll vmax = X / (1 + u) - 1;   /* (1+u)(1+v) <= X */
            for (ll v0 = u; v0 <= vmax; v0 += CH) {
                if (nitems == capi) { capi = capi ? 2 * capi : 1 << 16; iu = realloc(iu, sizeof(ll) * capi); iv0 = realloc(iv0, sizeof(ll) * capi); }
                iu[nitems] = u; iv0[nitems] = v0; nitems++;
            }
        }
        #pragma omp parallel for schedule(dynamic, 1) reduction(+:npairs)
        for (ll it = 0; it < nitems; it++) {
            ll u = iu[it], v0 = iv0[it];
            ll vmax = X / (1 + u) - 1; ll v1 = v0 + CH - 1; if (v1 > vmax) v1 = vmax;
            ll gu = G(u), eu = (u == 0 || !(u & 1));
            for (ll v = v0; v <= v1; v++) {
                ll ev = (v == 0 || !(v & 1));
                bump(r1, (1 + u) * (1 + v) + gu * G(v) - eu * ev); npairs++;
            }
        }
        free(iu); free(iv0);
        printf("F1: %lld pairs with (1+u)(1+v) <= X, %.1f s\n", npairs, omp_get_wtime() - t0); fflush(stdout);
    }
    ll nev = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:nev)
    for (int i = 0; i < nb; i++) {
        ll a = yh[i], b = ye[i], c = yR[i];
        for (ll v = 0; (1 + v) * c <= X; v++) {
            ll ev = (v == 0 || !(v & 1));
            bump(r2, a * G(v) - b * ev + c * (1 + v)); nev++;
        }
    }
    printf("F2: %lld evaluations with (1+v) y_R <= X, %.1f s; even values produced: %lld; counter saturations: %lld\n", nev, omp_get_wtime() - t0, neven, nsat);
    fflush(stdout);
    /* per-block statistics */
    printf("# k | mean r1 | incr | 2kappa=mean r1-(pi^2/12)(k+2)ln2 | #r1=0 | mean r2 | min r2 (argmin) | #r2=0 | frac r2<=10 <=12 <=20 | mean r2 by N mod 3 | max r1 max r2 | cum #pairs<=2^(k+1)\n");
    double prev = -1; const double c0 = M_PI * M_PI / 12.0 * log(2.0);
    ll cum1 = 0;
    for (ll N = 1; N < 2; N += 2) cum1 += r1[N >> 1];
    for (int k = 1; k < LOGX; k++) {
        ll lo = 1LL << k, hi = 1LL << (k + 1), n = 0, z1 = 0, z2 = 0, mn = 1 << 30, arg = -1, c10 = 0, c12 = 0, c20 = 0, mx1 = 0, mx2 = 0;
        double s1 = 0, s2 = 0, sm[3] = {0}; ll cm[3] = {0};
        for (ll N = lo + 1; N < hi; N += 2) {
            ll a = r1[N >> 1], b = r2[N >> 1]; n++; s1 += a; s2 += b;
            if (!a) z1++; if (!b) z2++;
            if (b < mn) { mn = b; arg = N; }
            if (b <= 10) c10++; if (b <= 12) c12++; if (b <= 20) c20++;
            if (a > mx1) mx1 = a; if (b > mx2) mx2 = b;
            sm[N % 3] += b; cm[N % 3]++;
        }
        cum1 += (ll)s1;
        double m1 = s1 / n;
        printf("k=%2d | %.4f | %+.4f | %+.4f | %lld | %.4f | %lld (%lld) | %lld | %.3e %.3e %.3e | %.3f %.3f %.3f | %lld %lld | %lld\n",
               k, m1, prev < 0 ? 0.0 : m1 - prev, m1 - c0 * (k + 2), z1, s2 / n, mn, arg, z2, (double)c10 / n, (double)c12 / n, (double)c20 / n,
               cm[0] ? sm[0] / cm[0] : 0, cm[1] ? sm[1] / cm[1] : 0, cm[2] ? sm[2] / cm[2] : 0, mx1, mx2, cum1);
        prev = m1;
    }
    ll miss = 0; printf("odd N <= 2^%d with r2 = 0:", LOGX);
    for (ll N = 1; N <= X; N += 2) if (!r2[N >> 1]) { miss++; if (miss <= 60) printf(" %lld", N); }
    printf("  [total %lld]\n", miss);
    if (doF1) {
        ll mu = 0; printf("odd N <= 2^%d with r1 + r2 = 0:", LOGX);
        for (ll N = 1; N <= X; N += 2) if (!r2[N >> 1] && !r1[N >> 1]) { mu++; if (mu <= 60) printf(" %lld", N); }
        printf("  [total %lld]\n", mu);
    }
    /* residue statistics in block statk */
    for (int kk = statk; kk <= statk + 2 && kk < LOGX; kk += 2) {
        ll lo = 1LL << kk, hi = 1LL << (kk + 1);
#ifdef MODS2
        int mods[] = {3, 5, 7, 25, 49, 21};   /* build with -DMODS2 for the 5- and 7-adic class means */
#else
        int mods[] = {3, 8, 16, 9, 27, 81};
#endif
        for (int mi = 0; mi < 6; mi++) {
            int m = mods[mi]; double s1[81] = {0}, s2[81] = {0}; ll c[81] = {0}, z1[81] = {0};
            for (ll N = lo + 1; N < hi; N += 2) { int r = (int)(N % m); s1[r] += r1[N >> 1]; s2[r] += r2[N >> 1]; c[r]++; if (!r1[N >> 1]) z1[r]++; }
            double mn1 = 1e9, mx1 = 0, mn2 = 1e9, mx2 = 0; int amn1 = -1;
            for (int r = 0; r < m; r++) if (c[r]) { double a = s1[r] / c[r], b = s2[r] / c[r]; if (a < mn1) { mn1 = a; amn1 = r; } if (a > mx1) mx1 = a; if (b < mn2) mn2 = b; if (b > mx2) mx2 = b; }
            printf("k=%d mod %2d: r2 class means in [%.3f, %.3f]; r1 class means in [%.3f, %.3f], min at r=%d (F1-miss frac there %.5f)", kk, m, mn2, mx2, mn1, mx1, amn1, (double)z1[amn1] / c[amn1]);
            if (m == 3) printf("; r2 = %.3f %.3f %.3f; r1 = %.3f %.3f %.3f; F1-miss frac = %.5f %.5f %.5f", s2[0] / c[0], s2[1] / c[1], s2[2] / c[2], s1[0] / c[0], s1[1] / c[1], s1[2] / c[2], (double)z1[0] / c[0], (double)z1[1] / c[1], (double)z1[2] / c[2]);
#ifdef MODS2
            if (m == 5 || m == 7) { printf("\n    r2 by class:"); for (int r = 0; r < m; r++) printf(" %d:%.3f", r, s2[r] / c[r]); printf("\n    r1 by class:"); for (int r = 0; r < m; r++) printf(" %d:%.3f", r, s1[r] / c[r]); }
#endif
            printf("\n");
        }
    }
    printf("elapsed %.1f s\n", omp_get_wtime() - t0);
    return 0;
}
