/*
 * aud82_onepex.c  (independent audit of opus-2026-11-06-S82)
 *
 * (1) One-apex triple: T(v) = transitive chain c_0..c_{L-1} (c_i -> c_j iff i > j) plus an apex x with x -> c_i iff
 *     bit i of v is 1.  For every L <= LMAX and every v < 2^L, the triple (h, e_x, R_x) is computed by subset DP and
 *     compared with (g(v), [v = 0 or v even], 1 + v), g from the auditor's closed form
 *         g(v) = 1 + sum over zero bits i < top(v) of ( floor(v / 2^(i+1)) + bit_(i+1)(v) ),   g(0) = 1,
 *     and the closed form is compared with the recursion g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2) for v < 2^VMAXB.
 * (2) Linearisation: delta(v) = g(v) - v rho(v), rho(v) = sum over zero bits i < top(v) of 2^-(i+1); per top bit n the
 *     max and min of delta (exact rational arithmetic: 2^n delta is an integer... we use long double, exact here).
 * (3) One-arc join law: A = T(v) rooted at its apex a, B a random tournament of order <= 7 rooted at b; the join
 *     A => B with b -> a reversed; H counted directly and compared with x_h y_h - x_e y_e + x_R y_R, where
 *     x = triple of (A, a) and y = triple of (B^op, b).  Also F1: B = T(u)^op rooted at its apex.
 * Usage: aud82_onepex LMAX VMAXB NJOIN
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef long long ll;

static ll g_closed(ll v) {
    if (v == 0) return 1;
    int top = 63 - __builtin_clzll((unsigned long long)v);
    ll s = 1;
    for (int i = 0; i < top; i++) if (!((v >> i) & 1)) s += (v >> (i + 1)) + ((v >> (i + 1)) & 1);
    return s;
}
static ll *grec; /* recursion table */

#define MAXV 20
static ll endp[1 << MAXV][MAXV];
static ll count_paths(int n, const unsigned *out, ll *e_at, int t) {   /* H(T); also #paths ending at t */
    int full = (1 << n) - 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) endp[S][v] = 0;
    for (int v = 0; v < n; v++) endp[1 << v][v] = 1;
    for (int S = 1; S <= full; S++)
        for (int v = 0; v < n; v++) {
            ll c = endp[S][v]; if (!c) continue;
            unsigned nb = out[v] & ~(unsigned)S;
            while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; endp[S | (1 << w)][w] += c; }
        }
    ll h = 0; for (int v = 0; v < n; v++) h += endp[full][v];
    if (e_at) *e_at = endp[full][t];
    return h;
}
/* triple at t (small n only) */
static ll hsub[1 << 12], stsub[1 << 12];
static void triple(int n, const unsigned *out, int t, ll *h, ll *e, ll *R) {
    static ll ep[1 << 12][12], sp[1 << 12][12];
    int full = (1 << n) - 1;
    memset(ep, 0, sizeof ep); memset(sp, 0, sizeof sp);
    for (int v = 0; v < n; v++) ep[1 << v][v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { ll c = ep[S][v]; if (!c) continue; unsigned nb = out[v] & ~(unsigned)S; while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; ep[S | (1 << w)][w] += c; } }
    sp[1 << t][t] = 1;
    for (int S = 1; S <= full; S++) { if (!((S >> t) & 1)) continue; for (int v = 0; v < n; v++) { ll c = sp[S][v]; if (!c) continue; unsigned nb = out[v] & ~(unsigned)S; while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; sp[S | (1 << w)][w] += c; } } }
    for (int S = 0; S <= full; S++) { hsub[S] = 0; stsub[S] = 0; for (int v = 0; v < n; v++) { hsub[S] += ep[S][v]; stsub[S] += sp[S][v]; } }
    hsub[0] = 1;
    *h = hsub[full]; *e = ep[full][t]; *R = 0;
    for (int Y = 1; Y <= full; Y++) if ((Y >> t) & 1) *R += hsub[full & ~Y] * stsub[Y];
}
static void one_apex(int L, ll v, unsigned *out) {   /* vertices 0..L-1 chain, L = apex */
    for (int i = 0; i <= L; i++) out[i] = 0;
    for (int i = 0; i < L; i++) for (int j = 0; j < i; j++) out[i] |= 1u << j;   /* c_i -> c_j iff i > j */
    for (int i = 0; i < L; i++) { if ((v >> i) & 1) out[L] |= 1u << i; else out[i] |= 1u << L; }
}
static unsigned long long rs = 0x9E3779B97F4A7C15ULL;
static unsigned long long rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

int main(int argc, char **argv) {
    int LMAX = argc > 1 ? atoi(argv[1]) : 12, VB = argc > 2 ? atoi(argv[2]) : 24, NJ = argc > 3 ? atoi(argv[3]) : 400;
    ll VM = 1LL << VB;
    grec = malloc(sizeof(ll) * (VM + 1)); grec[0] = 1;
    for (ll v = 1; v <= VM; v++) grec[v] = (v & 1) ? grec[v >> 1] : grec[v >> 1] + (v >> 1) + ((v >> 1) & 1);
    ll bad = 0; for (ll v = 0; v <= VM; v++) if (grec[v] != g_closed(v)) bad++;
    printf("(0) closed form vs recursion, v <= 2^%d: %lld failures\n", VB, bad);
    ll bad2 = 0; for (ll v = 0; v <= VM; v++) if (grec[v] > v + 1 || grec[v] < 1) bad2++;
    printf("    1 <= g(v) <= v + 1 for v <= 2^%d: %lld failures\n", VB, bad2);
    /* (1) one-apex triple */
    ll checked = 0, fail = 0;
    for (int L = 0; L <= LMAX; L++)
        for (ll v = 0; v < (1LL << L); v++) {
            unsigned out[32]; one_apex(L, v, out);
            ll h, e, R;
            if (L + 1 <= 9) triple(L + 1, out, L, &h, &e, &R);
            else {  /* larger L: h and e only by the big DP; R by the decomposition sum (direct) */
                ll ee; h = count_paths(L + 1, out, &ee, L); e = ee; R = -1;
            }
            checked++;
            ll ev = (v == 0 || !(v & 1));
            if (h != g_closed(v) || e != ev || (R >= 0 && R != 1 + v)) { fail++; if (fail < 10) printf("   FAIL L=%d v=%lld: (%lld,%lld,%lld)\n", L, v, h, e, R); }
        }
    printf("(1) one-apex triple (h, e_x, R_x) = (g(v), [v=0 or even], 1+v) for all L <= %d, v < 2^L (R checked for L <= 8): %lld cases, %lld failures\n", LMAX, checked, fail);
    /* (2) delta range per top bit */
    printf("(2) delta(v) = g(v) - v rho(v) per top bit n: [min, max]   (note's bound |delta| <= 2n + 1)\n");
    for (int n = 1; n <= VB - 1; n++) {
        long double mn = 1e30L, mx = -1e30L;
        for (ll v = 1LL << n; v < (2LL << n); v++) {
            long double rho = 0; for (int i = 0; i < n; i++) if (!((v >> i) & 1)) rho += ldexpl(1.0L, -(i + 1));
            long double d = (long double)grec[v] - (long double)v * rho;
            if (d < mn) mn = d; if (d > mx) mx = d;
        }
        if (n <= 6 || n % 4 == 0 || n == VB - 1) printf("    n=%2d: [%.4Lf, %.4Lf]   bound %d\n", n, mn, mx, 2 * n + 1);
    }
    /* (3) join law */
    ll jf = 0, jc = 0, nodd = 0;
    for (int it = 0; it < NJ; it++) {
        int L = 1 + rnd() % 9; ll v = rnd() % (1LL << L);
        int m = 1 + rnd() % 7;
        unsigned A[32], B[16], Bop[16], T[32];
        one_apex(L, v, A);
        for (int i = 0; i < m; i++) B[i] = 0;
        for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (rnd() & 1) B[i] |= 1u << j; else B[j] |= 1u << i; }
        int b = rnd() % m;
        for (int i = 0; i < m; i++) { Bop[i] = 0; for (int j = 0; j < m; j++) if ((B[j] >> i) & 1) Bop[i] |= 1u << j; }
        ll xh, xe, xR, yh, ye, yR;
        triple(L + 1, A, L, &xh, &xe, &xR);
        triple(m, Bop, b, &yh, &ye, &yR);
        int n = L + 1 + m;   /* A: 0..L (apex L), B: L+1 .. L+m */
        for (int i = 0; i <= L; i++) T[i] = A[i] | (((1u << m) - 1) << (L + 1));
        for (int i = 0; i < m; i++) T[L + 1 + i] = B[i] << (L + 1);
        int a = L, bb = L + 1 + b;
        T[a] &= ~(1u << bb); T[bb] |= 1u << a;   /* reverse a -> b into b -> a */
        ll H = count_paths(n, T, NULL, 0);
        ll pred = xh * yh - xe * ye + xR * yR;
        ll pred2 = yh * g_closed(v) - ye * (v == 0 || !(v & 1)) + yR * (1 + v);
        jc++; if (H != pred || H != pred2) { jf++; if (jf < 10) printf("   JOIN FAIL L=%d v=%lld m=%d: H=%lld pred=%lld pred2=%lld\n", L, v, m, H, pred, pred2); }
        if (H & 1) nodd++;
    }
    printf("(3a) one-arc join law, one-apex x random rooted tournament of order <= 7: %lld random cases, %lld failures, %lld odd\n", jc, jf, nodd);
    jf = 0; jc = 0;
    for (int it = 0; it < NJ; it++) {
        int L1 = rnd() % 9, L2 = rnd() % 9; ll u = rnd() % (1LL << L1), v = rnd() % (1LL << L2);
        unsigned A[32], C[32], Cop[32], T[32];
        one_apex(L1, u, A); one_apex(L2, v, C);
        int m = L2 + 1;
        for (int i = 0; i < m; i++) { Cop[i] = 0; for (int j = 0; j < m; j++) if ((C[j] >> i) & 1) Cop[i] |= 1u << j; }
        /* B = T(v)^op rooted at its apex (index L2); its converse is T(v), whose apex triple is (g(v), e(v), 1+v) */
        int n = L1 + 1 + m;
        for (int i = 0; i <= L1; i++) T[i] = A[i] | (((1u << m) - 1) << (L1 + 1));
        for (int i = 0; i < m; i++) T[L1 + 1 + i] = Cop[i] << (L1 + 1);
        int a = L1, bb = L1 + 1 + L2;
        T[a] &= ~(1u << bb); T[bb] |= 1u << a;
        ll H = count_paths(n, T, NULL, 0);
        ll eu = (u == 0 || !(u & 1)), ev = (v == 0 || !(v & 1));
        ll pred = (1 + u) * (1 + v) + g_closed(u) * g_closed(v) - eu * ev;
        jc++; if (H != pred) { jf++; if (jf < 10) printf("   F1 FAIL u=%lld v=%lld: H=%lld pred=%lld\n", u, v, H, pred); }
    }
    printf("(3b) F1 law, T(u) joined with T(v)^op at the apexes: %lld random cases, %lld failures\n", jc, jf);
    return 0;
}
