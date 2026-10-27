/* Sums of four fourth powers of Eisenstein integers equal to a rational integer n (HYP-9429; opus-2026-10-27-S66).
 *
 * Elements are pairs (x, y) for x + y*omega, omega^2 + omega + 1 = 0:
 *   (a, b)(c, d) = (ac - bd, ad + bc - bd),  conj(x, y) = (x - y, -y).
 *
 * mode "all" B lo hi:   bases a + b*omega with |a|, |b| <= B; every pair sum alpha^4 + beta^4 (alpha, beta up to sign)
 *                       sorted by (y, x); n = P + Q is rational iff y_Q = -y_P; for each n in [lo, hi] report a witness
 *                       or "none".
 * mode "pat" B n e1 e2: the class-restricted search used for n = 6 (mod 16) [e1 = 1, e2 = 9] and n = 14 (mod 16):
 *                       by THM-4826 (6) a representation of n = 6 (mod 8) has two bases = omega and two = omega^2
 *                       (mod 2); for n = 6 (mod 16) each pair has fourth powers = omega*1 and omega*9 (mod 16).
 *                       Conjugation maps the omega^2 pair to an omega pair, so n = P + conj(Q) with P, Q in
 *                       S = {alpha^4 + beta^4 : alpha^4 = e1*omega, beta^4 = e2*omega (mod 16)}, i.e. y_P = y_Q and
 *                       x_P + x_Q = n + y_P.  Every hit is re-verified by recomputing the four fourth powers.
 *                       Note: this covers omega^2 bases in the conjugate of the box, not the box itself.
 * mode "box" B n:       the exact box version for n = 6 (mod 16): omega-class bases (a even, b odd) and omega^2-class
 *                       bases (a, b odd), all with |a|, |b| <= B; S1 = {alpha^4 + beta^4 : alpha^4 = omega, beta^4 =
 *                       9 omega (mod 16)}, S2 = {gamma^4 + delta^4 : gamma^4 = omega^2, delta^4 = 9 omega^2 (mod 16)};
 *                       n = P + Q with P in S1, Q in S2 iff y_Q = -y_P and x_P + x_Q = n.
 * All arithmetic is exact in int64 (|coordinates| < 2^40 for B <= 150, checked).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef long long ll;
typedef struct { ll x, y; int i, j; } pr;           /* pair sum and the indices of its two terms */
typedef struct { ll x, y; int a, b; } pw;           /* fourth power and its base */

static void mul(ll a, ll b, ll c, ll d, ll *x, ll *y) { *x = a * c - b * d; *y = a * d + b * c - b * d; }
static void p4(ll a, ll b, ll *x, ll *y) { ll s, t; mul(a, b, a, b, &s, &t); mul(s, t, s, t, x, y); }
static ll md(ll v, ll m) { ll r = v % m; return r < 0 ? r + m : r; }

static int cmp_pr(const void *p, const void *q) {
    const pr *u = p, *v = q;
    if (u->y != v->y) return u->y < v->y ? -1 : 1;
    if (u->x != v->x) return u->x < v->x ? -1 : 1;
    return 0;
}
static int cmp_pw(const void *p, const void *q) {
    const pw *u = p, *v = q;
    if (u->x != v->x) return u->x < v->x ? -1 : 1;
    if (u->y != v->y) return u->y < v->y ? -1 : 1;
    return 0;
}

/* all fourth powers of bases with |a|,|b| <= B, deduplicated (alpha and -alpha give the same value).
   pattern 1: alpha = omega (mod 2) [a even, b odd] and alpha^4 = e*omega (mod 16), i.e. (0, e);
   pattern 2: alpha = omega^2 (mod 2) [a, b odd] and alpha^4 = e*omega^2 = (-e, -e) (mod 16). */
static pw *powers(int B, int *count, int pattern, ll e) {
    int cap = (2 * B + 1) * (2 * B + 1), n = 0;
    pw *P = malloc(sizeof(pw) * cap);
    for (int a = -B; a <= B; a++)
        for (int b = -B; b <= B; b++) {
            if (a < 0 || (a == 0 && b < 0)) continue;          /* one of +-alpha (the box is symmetric) */
            ll x, y; p4(a, b, &x, &y);
            if (pattern == 1) {
                if (md(a, 2) != 0 || md(b, 2) != 1) continue;
                if (md(x, 16) != 0 || md(y, 16) != e) continue;
            }
            if (pattern == 2) {
                if (md(a, 2) != 1 || md(b, 2) != 1) continue;
                if (md(x, 16) != md(-e, 16) || md(y, 16) != md(-e, 16)) continue;
            }
            if (x > (1LL << 40) || x < -(1LL << 40) || y > (1LL << 40) || y < -(1LL << 40)) { fprintf(stderr, "overflow\n"); exit(2); }
            P[n].x = x; P[n].y = y; P[n].a = a; P[n].b = b; n++;
        }
    qsort(P, n, sizeof(pw), cmp_pw);
    int m = 0;
    for (int k = 0; k < n; k++)
        if (m == 0 || P[k].x != P[m - 1].x || P[k].y != P[m - 1].y) P[m++] = P[k];
    *count = m;
    return P;
}

static int verify(ll n, int a[4], int b[4]) {
    ll sx = 0, sy = 0;
    for (int k = 0; k < 4; k++) { ll x, y; p4(a[k], b[k], &x, &y); sx += x; sy += y; }
    return sx == n && sy == 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 1; }
    if (!strcmp(argv[1], "all")) {
        int B = atoi(argv[2]); ll lo = atoll(argv[3]), hi = atoll(argv[4]);
        int np; pw *P = powers(B, &np, 0, 0);
        size_t npr = (size_t)np * (np + 1) / 2, k = 0;
        pr *S = malloc(sizeof(pr) * npr);
        if (!S) { fprintf(stderr, "alloc\n"); return 2; }
        for (int i = 0; i < np; i++)
            for (int j = i; j < np; j++) { S[k].x = P[i].x + P[j].x; S[k].y = P[i].y + P[j].y; S[k].i = i; S[k].j = j; k++; }
        qsort(S, npr, sizeof(pr), cmp_pr);
        int width = (int)(hi - lo + 1);
        long *wit = malloc(sizeof(long) * 2 * width);
        for (int t = 0; t < 2 * width; t++) wit[t] = -1;
        /* for each P (index p), look in the bucket y = -y_P for x in [lo - x_P, hi - x_P] */
        size_t p = 0;
        while (p < npr) {
            size_t q0 = p; while (p < npr && S[p].y == S[q0].y) p++;      /* bucket [q0, p) with y = Y */
            ll Y = S[q0].y;
            if (Y > 0) continue;                                          /* handle each {Y, -Y} once (Y <= 0) */
            /* find bucket -Y */
            size_t L = 0, R = npr;
            while (L < R) { size_t M = (L + R) / 2; if (S[M].y < -Y) L = M + 1; else R = M; }
            size_t b0 = L; while (R < npr && S[R].y == -Y) R++;
            size_t b1 = R;
            if (b0 == b1 || S[b0].y != -Y) continue;
            for (size_t u = q0; u < p; u++) {
                ll xlo = lo - S[u].x, xhi = hi - S[u].x;
                size_t l2 = b0, r2 = b1;
                while (l2 < r2) { size_t M = (l2 + r2) / 2; if (S[M].x < xlo) l2 = M + 1; else r2 = M; }
                for (size_t v = l2; v < b1 && S[v].x <= xhi; v++) {
                    ll n = S[u].x + S[v].x;
                    int t = (int)(n - lo);
                    if (wit[2 * t] < 0) { wit[2 * t] = (long)u; wit[2 * t + 1] = (long)v; }
                }
            }
        }
        int none = 0;
        for (int t = 0; t < width; t++) {
            ll n = lo + t;
            if (wit[2 * t] < 0) { printf("n = %lld: none\n", n); none++; continue; }
            pr U = S[wit[2 * t]], V = S[wit[2 * t + 1]];
            int a[4] = {P[U.i].a, P[U.j].a, P[V.i].a, P[V.j].a}, b[4] = {P[U.i].b, P[U.j].b, P[V.i].b, P[V.j].b};
            if (!verify(n, a, b)) { printf("VERIFY FAILED at %lld\n", n); return 3; }
            printf("n = %lld: (%d,%d) (%d,%d) (%d,%d) (%d,%d)\n", n, a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3]);
        }
        printf("B = %d: %d fourth powers, %zu pair sums; %d of %d targets without a representation\n", B, np, npr, none, width);
        return 0;
    }
    if (!strcmp(argv[1], "pat")) {
        int B = atoi(argv[2]); ll n = atoll(argv[3]); ll e1 = atoll(argv[4]), e2 = atoll(argv[5]);
        int n1, n2; pw *P1 = powers(B, &n1, 1, e1), *P2 = powers(B, &n2, 1, e2);
        size_t npr = (size_t)n1 * n2, k = 0;
        if (e1 == e2) npr = (size_t)n1 * (n1 + 1) / 2;
        pr *S = malloc(sizeof(pr) * npr);
        if (!S) { fprintf(stderr, "alloc\n"); return 2; }
        if (e1 == e2) {
            for (int i = 0; i < n1; i++) for (int j = i; j < n1; j++) { S[k].x = P1[i].x + P1[j].x; S[k].y = P1[i].y + P1[j].y; S[k].i = i; S[k].j = j; k++; }
        } else {
            for (int i = 0; i < n1; i++) for (int j = 0; j < n2; j++) { S[k].x = P1[i].x + P2[j].x; S[k].y = P1[i].y + P2[j].y; S[k].i = i; S[k].j = j; k++; }
        }
        qsort(S, npr, sizeof(pr), cmp_pr);
        long hits = 0; int shown = 0;
        size_t p = 0;
        while (p < npr) {
            size_t q0 = p; while (p < npr && S[p].y == S[q0].y) p++;
            ll Y = S[q0].y, target = n + Y;                 /* x_P + x_Q = n + Y within the bucket */
            size_t l = q0, r = p - 1;
            while (l <= r && r < npr && l < p) {
                ll s = S[l].x + S[r].x;
                if (s == target) {
                    hits++;
                    if (shown < 3) {
                        pw A = P1[S[l].i], Bq = (e1 == e2 ? P1[S[l].j] : P2[S[l].j]);
                        pw C = P1[S[r].i], D = (e1 == e2 ? P1[S[r].j] : P2[S[r].j]);
                        /* the omega^2 pair is the conjugate of (C, D): conj(a + b omega) = (a - b) - b omega */
                        int a[4] = {A.a, Bq.a, C.a - C.b, D.a - D.b}, b[4] = {A.b, Bq.b, -C.b, -D.b};
                        if (!verify(n, a, b)) { printf("VERIFY FAILED\n"); return 3; }
                        printf("  hit: (%d,%d) (%d,%d) (%d,%d) (%d,%d)\n", a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3]);
                        shown++;
                    }
                    if (l == r) break;
                    l++;
                } else if (s < target) l++;
                else { if (r == 0) break; r--; }
            }
        }
        printf("pattern search n = %lld, e = (%lld, %lld), B = %d: class sizes %d, %d; %zu pair sums; hits %ld\n", n, e1, e2, B, n1, n2, npr, hits);
        return 0;
    }
    if (!strcmp(argv[1], "box")) {
        int B = atoi(argv[2]); ll n = atoll(argv[3]);
        int c1, c9, d1, d9;
        pw *A1 = powers(B, &c1, 1, 1), *A9 = powers(B, &c9, 1, 9), *G1 = powers(B, &d1, 2, 1), *G9 = powers(B, &d9, 2, 9);
        size_t n1 = (size_t)c1 * c9, n2 = (size_t)d1 * d9, k = 0;
        pr *S1 = malloc(sizeof(pr) * n1), *S2 = malloc(sizeof(pr) * n2);
        if (!S1 || !S2) { fprintf(stderr, "alloc\n"); return 2; }
        for (int i = 0; i < c1; i++) for (int j = 0; j < c9; j++) { S1[k].x = A1[i].x + A9[j].x; S1[k].y = A1[i].y + A9[j].y; S1[k].i = i; S1[k].j = j; k++; }
        k = 0;
        for (int i = 0; i < d1; i++) for (int j = 0; j < d9; j++) { S2[k].x = G1[i].x + G9[j].x; S2[k].y = G1[i].y + G9[j].y; S2[k].i = i; S2[k].j = j; k++; }
        qsort(S1, n1, sizeof(pr), cmp_pr);
        qsort(S2, n2, sizeof(pr), cmp_pr);
        long hits = 0; int shown = 0;
        size_t p = 0;
        while (p < n1) {
            size_t q0 = p; while (p < n1 && S1[p].y == S1[q0].y) p++;
            ll Y = S1[q0].y;
            size_t L = 0, R = n2;                                   /* bucket y = -Y in S2 */
            while (L < R) { size_t M = (L + R) / 2; if (S2[M].y < -Y) L = M + 1; else R = M; }
            if (L == n2 || S2[L].y != -Y) continue;
            size_t b0 = L, b1 = L; while (b1 < n2 && S2[b1].y == -Y) b1++;
            for (size_t u = q0; u < p; u++) {
                ll want = n - S1[u].x;
                size_t l2 = b0, r2 = b1;
                while (l2 < r2) { size_t M = (l2 + r2) / 2; if (S2[M].x < want) l2 = M + 1; else r2 = M; }
                for (size_t v = l2; v < b1 && S2[v].x == want; v++) {
                    hits++;
                    if (shown < 3) {
                        int a[4] = {A1[S1[u].i].a, A9[S1[u].j].a, G1[S2[v].i].a, G9[S2[v].j].a};
                        int b[4] = {A1[S1[u].i].b, A9[S1[u].j].b, G1[S2[v].i].b, G9[S2[v].j].b};
                        if (!verify(n, a, b)) { printf("VERIFY FAILED\n"); return 3; }
                        printf("  hit: (%d,%d) (%d,%d) (%d,%d) (%d,%d)\n", a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3]);
                        shown++;
                    }
                }
            }
        }
        printf("box search n = %lld (pattern omega{1,9} + omega^2{1,9}), |a|,|b| <= %d: class sizes %d, %d, %d, %d; pair sums %zu, %zu; hits %ld\n",
               n, B, c1, c9, d1, d9, n1, n2, hits);
        return 0;
    }
    return 1;
}
