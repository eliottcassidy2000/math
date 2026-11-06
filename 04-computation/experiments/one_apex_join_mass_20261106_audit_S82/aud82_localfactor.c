/*
 * aud82_localfactor.c  (independent audit of opus-2026-11-06-S82)
 *
 * Does the local factor of r1, r2 stay bounded below when N avoids every small prime?  Per dyadic block k, the mean of
 * r1 and r2 over odd N with gcd(N, P_j) = 1 for P_j = 3, 15, 105, 1155, 15015, 255255, and over prime N, plus the
 * block minima over primes.  (The note, section 5: "the 3-adic deficit is a constant-factor local density, like a
 * singular series, not an unbounded depletion"; section 7: "S is 2-adically trivial and 3-adically bounded".)
 * Same g / split law / sieve as aud82_sieve.c.
 * Usage: aud82_localfactor LOGX bankfile
 */
#include <stdio.h>
#include <stdlib.h>
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
static uint16_t *r1, *r2; static ll X;
static inline void bump(uint16_t *r, ll N) { if (N >= 1 && N <= X && (N & 1)) __atomic_fetch_add(&r[N >> 1], 1, __ATOMIC_RELAXED); }
static ll gcdll(ll a, ll b) { while (b) { ll t = a % b; a = b; b = t; } return a; }
int main(int argc, char **argv) {
    int LOGX = atoi(argv[1]); X = 1LL << LOGX;
    FILE *f = fopen(argv[2], "r"); int cap = 1 << 17; ll *yh = malloc(8 * cap), *ye = malloc(8 * cap), *yR = malloc(8 * cap); int nb = 0;
    while (nb < cap && fscanf(f, "%lld %lld %lld", &yh[nb], &ye[nb], &yR[nb]) == 3) nb++;
    fclose(f);
    for (ll v = 0; v < SN; v++) gs[v] = g_closed(v);
    for (int r = 0; r < SN; r++) {
        ll C = 0, K = 0;
        for (int i = 0; i < SB; i++) if (!((r >> i) & 1)) { C += 1LL << (SB - 1 - i); K += r >> (i + 1); }
        for (int i = 0; i + 1 < SB; i++) if (((r >> (i + 1)) & 1) && !((r >> i) & 1)) K++;
        Ct[r] = C; Et[r] = !((r >> (SB - 1)) & 1); Kt[r] = K;
    }
    r1 = calloc((size_t)(X / 2 + 1), 2); r2 = calloc((size_t)(X / 2 + 1), 2);
    ll umax = 0; while ((umax + 1) * (umax + 1) <= X) umax++; umax--;
    #pragma omp parallel for schedule(dynamic, 1)
    for (ll u = 0; u <= umax; u++) {
        ll gu = G(u), eu = (u == 0 || !(u & 1)), vmax = X / (1 + u) - 1;
        for (ll v = u; v <= vmax; v++) bump(r1, (1 + u) * (1 + v) + gu * G(v) - eu * (v == 0 || !(v & 1)));
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < nb; i++)
        for (ll v = 0; (1 + v) * yR[i] <= X; v++) bump(r2, yh[i] * G(v) - ye[i] * (v == 0 || !(v & 1)) + yR[i] * (1 + v));
    /* primes up to X: odd sieve */
    uint8_t *comp = calloc((size_t)(X / 2 + 1), 1);
    for (ll p = 3; p * p <= X; p += 2) if (!comp[p >> 1]) for (ll q = p * p; q <= X; q += 2 * p) comp[q >> 1] = 1;
    ll P[] = {3, 15, 105, 1155, 15015, 255255};
    printf("# aud82_localfactor LOGX=%d bank=%s (%d triples)\n", LOGX, argv[2], nb);
    printf("# k | mean r1, r2 over: all odd | gcd(N,3)=1 | gcd(N,15)=1 | gcd(N,105)=1 | gcd(N,1155)=1 | gcd(N,15015)=1 | gcd(N,255255)=1 | primes || min r1, min r2 over primes | #primes\n");
    for (int k = 14; k < LOGX; k++) {
        ll lo = 1LL << k, hi = 1LL << (k + 1);
        double s1[8] = {0}, s2[8] = {0}; ll c[8] = {0}; ll mn1 = 1 << 30, mn2 = 1 << 30;
        for (ll N = lo + 1; N < hi; N += 2) {
            ll a = r1[N >> 1], b = r2[N >> 1];
            s1[0] += a; s2[0] += b; c[0]++;
            for (int j = 0; j < 6; j++) if (gcdll(N, P[j]) == 1) { s1[j + 1] += a; s2[j + 1] += b; c[j + 1]++; }
            if (!comp[N >> 1]) { s1[7] += a; s2[7] += b; c[7]++; if (a < mn1) mn1 = a; if (b < mn2) mn2 = b; }
        }
        printf("k=%2d |", k);
        for (int j = 0; j < 8; j++) printf(" %.3f,%.3f |", s1[j] / c[j], s2[j] / c[j]);
        printf("| %lld %lld | %lld\n", mn1, mn2, c[7]);
        fflush(stdout);
    }
    return 0;
}
