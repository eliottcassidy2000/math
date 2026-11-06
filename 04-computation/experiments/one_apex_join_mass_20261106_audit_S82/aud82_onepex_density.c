/*
 * aud82_onepex_density.c  (independent audit of opus-2026-11-06-S82)
 *
 * The one-apex value set S_1 = {g(v)} (THM-4876 (6), HYP-9469 (3), THM-4878 (5)): since g(2w+1) = g(w) and
 * g(2w) = g(w) + w + (w mod 2) >= w + 1, every value <= X is g(2w) for some w < X.  This program marks g(2w) <= X for
 * all w < X and prints, per dyadic block of odd N, the fraction attained and the mean multiplicity per odd number
 * (#{w : g(2w) in block} / #odd in block; Theorem A with (a,b,c) = (1,1,1) predicts 2 ln 2 = 1.386294).
 * Usage: aud82_onepex_density LOGX
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
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
    int LOGX = atoi(argv[1]); ll X = 1LL << LOGX;
    for (ll v = 0; v < SN; v++) gs[v] = g_closed(v);
    for (int r = 0; r < SN; r++) {
        ll C = 0, K = 0;
        for (int i = 0; i < SB; i++) if (!((r >> i) & 1)) { C += 1LL << (SB - 1 - i); K += r >> (i + 1); }
        for (int i = 0; i + 1 < SB; i++) if (((r >> (i + 1)) & 1) && !((r >> i) & 1)) K++;
        Ct[r] = C; Et[r] = !((r >> (SB - 1)) & 1); Kt[r] = K;
    }
    uint8_t *cnt = calloc((size_t)(X / 2 + 1), 1);   /* saturating multiplicity of g(2w) = N */
    ll even = 0;
    for (ll w = 0; w < X; w++) {
        ll N = G(2 * w);
        if (N > X) continue;
        if (!(N & 1)) { even++; continue; }
        if (cnt[N >> 1] < 255) cnt[N >> 1]++;
    }
    printf("# one-apex values g(2w) <= 2^%d; even values: %lld\n", LOGX, even);
    printf("# k | attained fraction of odd N in [2^k, 2^(k+1)) | missed | mean multiplicity per odd N (2 ln 2 = %.6f) | max mult\n", 2 * log(2.0));
    ll totatt = 0, totodd = 0;
    for (int k = 1; k < LOGX; k++) {
        ll lo = 1LL << k, hi = 1LL << (k + 1), n = 0, a = 0, s = 0, mx = 0;
        for (ll N = lo + 1; N < hi; N += 2) { n++; int c = cnt[N >> 1]; if (c) a++; s += c; if (c > mx) mx = c; }
        totatt += a; totodd += n;
        printf("k=%2d | %.6f | %lld | %.5f | %lld\n", k, (double)a / n, n - a, (double)s / n, mx);
    }
    /* below 2^21 total, for comparison with THM-4878 (5): 521,915 of the 1,048,576 odd numbers below 2^21 */
    if (LOGX >= 21) { ll a = 0; for (ll N = 1; N < (1LL << 21); N += 2) if (cnt[N >> 1]) a++; printf("attained odd N < 2^21: %lld of %lld\n", a, 1LL << 20); }
    { ll a = 0; for (ll N = 1; N < 200000; N += 2) if (cnt[N >> 1]) a++; printf("attained odd N < 200,000: %lld of 100,000\n", a); }
    printf("first missed odd N:"); int sh = 0; for (ll N = 1; N < 200 && sh < 30; N += 2) if (!cnt[N >> 1]) { printf(" %lld", N); sh++; } printf("\n");
    return 0;
}
