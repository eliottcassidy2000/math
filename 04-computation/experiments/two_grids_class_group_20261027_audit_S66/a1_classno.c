/* Independent class numbers for the audit of THM-4825 (3) and THM-4826 (2).
 *
 * Method 1 (complete count, no early exit): h(D) = number of reduced primitive forms (a, b, c), b^2 - 4ac = D < 0,
 *   |b| <= a <= c, b >= 0 if |b| = a or a = c.
 * Method 2 (analytic, for fundamental D < -4): h(D) = (1 / (2 - chi(2))) * sum_{1 <= a <= |D|/2} chi(a),
 *   chi = Kronecker symbol (D / .), computed with my own Jacobi-symbol routine.
 *
 * Modes:
 *   sweep DMAX AMAX : all squarefree d = 1 (mod 4), 1 <= d < DMAX: h(-4d) by method 1; list d with h = 2; for d < AMAX
 *                     also method 2 and compare; histogram of small class numbers.
 *   cn2 DABS        : all fundamental discriminants D with -DABS <= D < 0 and h(D) = 2 (method 1, cross-checked by
 *                     method 2); compare with OEIS A014603 by eye.
 *   ratio DMAX      : h(-4d)/h(-d) for squarefree d = 3 (mod 4), 3 <= d < DMAX (method 1 for both; -4d is a
 *                     non-fundamental discriminant of conductor 2).
 *   one D           : print h(D) and the reduced forms.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;

static ll gcdll(ll a, ll b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { ll t = a % b; a = b; b = t; } return a; }

static int squarefree(ll n) { for (ll p = 2; p * p <= n; p++) if (n % (p * p) == 0) return 0; return 1; }

static ll classno(ll D, int print) {        /* method 1 */
    ll h = 0, N = -D;
    for (ll a = 1; 3 * a * a <= N; a++) {
        for (ll b = -a + 1; b <= a; b++) {
            ll num = b * b - D;            /* = 4ac */
            if (num % (4 * a)) continue;
            ll c = num / (4 * a);
            if (c < a) continue;
            if (b < 0 && a == c) continue;
            if (gcdll(gcdll(a, b), c) != 1) continue;
            h++;
            if (print) printf("  (%lld, %lld, %lld)\n", a, b, c);
        }
    }
    return h;
}

static int jacobi(ll a, ll n) {             /* n odd positive */
    int s = 1;
    a %= n; if (a < 0) a += n;
    while (a) {
        while ((a & 1) == 0) { a >>= 1; ll r = n & 7; if (r == 3 || r == 5) s = -s; }
        ll t = a; a = n; n = t;
        if ((a & 3) == 3 && (n & 3) == 3) s = -s;
        a %= n;
    }
    return n == 1 ? s : 0;
}

static int kronecker(ll D, ll m) {          /* (D / m), m >= 1, D = 0 or 1 mod 4 */
    int s = 1;
    while ((m & 1) == 0) {
        m >>= 1;
        if ((D & 1) == 0) return 0;
        ll r = ((D % 8) + 8) % 8;
        if (r == 3 || r == 5) s = -s;
    }
    if (m == 1) return s;
    return s * jacobi(D, m);
}

static ll classno_analytic(ll D) {          /* fundamental D < -4 */
    ll N = -D, S = 0;
    for (ll a = 1; a <= N / 2; a++) S += kronecker(D, a);
    int c2 = kronecker(D, 2);
    ll den = 2 - c2;
    if (S % den) { fprintf(stderr, "non-integral analytic h at D = %lld\n", D); exit(3); }
    return S / den;
}

static int fundamental(ll D) {              /* D < 0 */
    ll N = -D;
    if ((N % 4) == 3) return squarefree(N);            /* D = 1 mod 4 */
    if (N % 4 == 0) { ll m = N / 4; if ((m % 4 == 1 || m % 4 == 2) && squarefree(m)) return 1; }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) return 1;
    if (!strcmp(argv[1], "sweep")) {
        ll DMAX = atoll(argv[2]), AMAX = atoll(argv[3]);
        ll count = 0, cmp = 0, hist[11] = {0};
        printf("squarefree d = 1 mod 4, d < %lld: d with h(-4d) = 2:", DMAX);
        for (ll d = 1; d < DMAX; d += 4) {
            if (!squarefree(d)) continue;
            count++;
            ll h = classno(-4 * d, 0);
            if (h <= 10) hist[h]++;
            if (h == 2) printf(" %lld", d);
            if (d > 1 && d < AMAX) {
                ll h2 = classno_analytic(-4 * d);
                if (h2 != h) { printf("\nMISMATCH d = %lld: %lld vs %lld\n", d, h, h2); return 4; }
                cmp++;
            }
        }
        printf("\ncount of squarefree d = 1 mod 4 below %lld (including d = 1): %lld\n", DMAX, count);
        printf("analytic class number formula agrees for all %lld values 1 < d < %lld\n", cmp, AMAX);
        printf("histogram h(-4d) = 1..10:");
        for (int k = 1; k <= 10; k++) printf(" %d:%lld", k, hist[k]);
        printf("\n");
        return 0;
    }
    if (!strcmp(argv[1], "cn2")) {
        ll DA = atoll(argv[2]);
        printf("fundamental D, -%lld <= D < 0, with h(D) = 2 (method 1; method 2 cross-check):", DA);
        ll n = 0;
        for (ll N = 3; N <= DA; N++) {
            ll D = -N;
            if (!fundamental(D)) continue;
            ll h = classno(D, 0);
            if (N > 4) { ll h2 = classno_analytic(D); if (h2 != h) { printf("\nMISMATCH D = %lld\n", D); return 4; } }
            if (h == 2) { printf(" %lld", D); n++; }
        }
        printf("\n  total %lld\n", n);
        return 0;
    }
    if (!strcmp(argv[1], "ratio")) {
        ll DMAX = atoll(argv[2]);
        int bad = 0, nvals = 0;
        for (ll d = 3; d < DMAX; d += 4) {
            if (!squarefree(d)) continue;
            ll h4 = classno(-4 * d, 0), h1 = classno(-d, 0);
            /* expected: 1 at d = 3; 3 if d = 3 mod 8, d > 3; 1 if d = 7 mod 8 */
            ll num = h4, den = h1;
            int expect = (d == 3) ? 1 : ((d % 8 == 3) ? 3 : 1);
            nvals++;
            if (num != expect * den) { printf("  ratio fails at d = %lld: h(-4d) = %lld, h(-d) = %lld\n", d, h4, h1); bad++; }
        }
        printf("h(-4d)/h(-d) for squarefree d = 3 mod 4 in [3, %lld): %d values, %d failures of the rule (1 at 3; 3 for d = 3 mod 8; 1 for d = 7 mod 8)\n", DMAX, nvals, bad);
        return 0;
    }
    if (!strcmp(argv[1], "one")) {
        ll D = atoll(argv[2]);
        ll h = classno(D, 1);
        printf("h(%lld) = %lld\n", D, h);
        return 0;
    }
    return 1;
}
