/* procgen_repunit_20260924_collatz_repunits.c   (needs GMP)
 *
 * Collatz statistics of repunit families, as functions of the length k.
 *   family M  : n = 2^k - 1          (base-2 repunits, Mersenne numbers), k = 2..K
 *   family W  : n = (2^k + 1)/3      (base -2 repunits of odd length, Wagstaff / minus-trunk numbers), odd k = 3..K
 *   family R3 : n = (3^k - 1)/2      (base-3 repunits), k = 2..K
 *   family R10: n = (10^k - 1)/9     (base-10 repunits), k = 2..K
 * For every n the shortcut orbit T(x) = x/2, (3x+1)/2 is followed to 1 (Syracuse-compressed) and one line
 *   FAM k tau a sigma log2h Dword Dbridge extra...
 * is printed:
 *   tau     total stopping time (T-steps to reach 1); a = number of odd steps;
 *   sigma   stopping time: least j >= 1 with T^j(n) < n (exact, also inside halving runs);
 *   log2h   log2 of the height max_j T^j(n);
 *   Dword   max_j |a_j - j a/tau|  (a_j = odd steps among the first j): discrepancy of the parity word
 *           from the balanced word of the same density;
 *   Dbridge max_j |ln T^j(n) - ln n (1 - j/tau)|: sup of the log-orbit bridge.
 * Family M also prints: [T^k(n) == 3^k - 1] [T^(k+1)(n) == (3^k-1)/2] [height == 3^k - 1]
 *   v2(3^k - 1), the stopping time of (3^k-1)/2 measured inside the same orbit, and the statistics of the
 *   tail = the orbit of (3^k-1)/2 = T^(k+1)(n): tau_tail a_tail Dword_tail Dbridge_tail log2(height of tail).
 * Family W prints [T(n) == 2^(k-1)+1] [T^k(n) == 3^((k-1)/2)+1].
 * usage: collatz_repunits FAM K
 * Session collatz-procgen-20260922, lane procgen_repunit_20260924.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

static double lnz(const mpz_t x) {           /* natural log of a positive mpz */
    long e; double d = mpz_get_d_2exp(&e, x);
    return log(d) + (double)e * M_LN2;
}

typedef struct { long j; double ln; long a; } Ck;
static Ck *ck; static long nck, cap;
static void push(long j, double ln, long a) {
    if (nck == cap) { cap = cap ? 2 * cap : 1 << 16; ck = realloc(ck, cap * sizeof(Ck)); }
    ck[nck].j = j; ck[nck].ln = ln; ck[nck].a = a; nck++;
}

/* follow the orbit of n0; ref2 (optional, may be NULL) is a second threshold whose first passage
 * after T-index jref is also recorded (used for the stopping time of (3^k-1)/2 inside the M orbit). */
static double tail_tau, tail_a, tail_dw, tail_db, tail_log2h;   /* statistics of the orbit from T-index jref on */
static void orbit(const mpz_t n0, long *tau, long *a_out, long *sigma, double *log2h, mpz_t hmax,
                  double *Dword, double *Dbridge, const mpz_t ref2, long jref, long *sigma2) {
    mpz_t hmax2; mpz_init(hmax2); if (ref2) mpz_set(hmax2, ref2);
    mpz_t x, y, t; mpz_inits(x, y, t, NULL);
    mpz_set(x, n0); mpz_set(hmax, n0);
    long j = 0, a = 0; *sigma = -1; if (sigma2) *sigma2 = -1;
    nck = 0; push(0, lnz(x), 0);
    if (mpz_even_p(x)) {
        long v = mpz_scan1(x, 0);
        for (long i = 1; i <= v && *sigma < 0; i++) { mpz_tdiv_q_2exp(t, x, i); if (mpz_cmp(t, n0) < 0) *sigma = i; }
        j += v; mpz_tdiv_q_2exp(x, x, v); push(j, lnz(x), a);
    }
    while (mpz_cmp_ui(x, 1) != 0) {
        mpz_mul_ui(y, x, 3); mpz_add_ui(y, y, 1); a++;
        long v = mpz_scan1(y, 0);
        mpz_tdiv_q_2exp(t, y, 1);
        if (mpz_cmp(t, hmax) > 0) mpz_set(hmax, t);
        if (ref2 && j + 1 > jref && mpz_cmp(t, hmax2) > 0) mpz_set(hmax2, t);
        push(j + 1, lnz(y) - M_LN2, a);
        if (*sigma < 0) {
            mpz_tdiv_q_2exp(t, y, v);
            if (mpz_cmp(t, n0) < 0)
                for (long i = 1; i <= v; i++) { mpz_tdiv_q_2exp(t, y, i); if (mpz_cmp(t, n0) < 0) { *sigma = j + i; break; } }
        }
        if (sigma2 && *sigma2 < 0 && j + v > jref) {        /* run indices j+1..j+v; only those > jref count */
            long i0 = jref - j + 1; if (i0 < 1) i0 = 1;
            mpz_tdiv_q_2exp(t, y, v);
            if (mpz_cmp(t, ref2) < 0)
                for (long i = i0; i <= v; i++) { mpz_tdiv_q_2exp(t, y, i); if (mpz_cmp(t, ref2) < 0) { *sigma2 = j + i - jref; break; } }
        }
        j += v; mpz_tdiv_q_2exp(x, y, v);
        push(j, lnz(x), a);
    }
    *tau = j; *a_out = a;
    { long e; double d = mpz_get_d_2exp(&e, hmax); *log2h = log2(d) + (double)e; }
    double rho = (double)a / (double)j, L0 = ck[0].ln, dw = 0, db = 0;
    for (long i = 0; i < nck; i++) {
        double w = fabs((double)ck[i].a - ck[i].j * rho); if (w > dw) dw = w;
        double b = fabs(ck[i].ln - L0 * (1.0 - (double)ck[i].j / (double)j)); if (b > db) db = b;
    }
    *Dword = dw; *Dbridge = db;
    if (ref2) {                                   /* tail = orbit of ref2 = T^jref(n0), started at index jref */
        long ta = -1; double lnr = lnz(ref2);
        for (long i = 0; i < nck; i++) if (ck[i].j <= jref) ta = ck[i].a;   /* odd steps before jref */
        double tt = (double)(j - jref), at = (double)(a - ta), r2 = at / tt, w0 = 0, b0 = 0;
        for (long i = 0; i < nck; i++) {
            if (ck[i].j < jref) continue;
            double jj = (double)(ck[i].j - jref), aa = (double)(ck[i].a - ta);
            double w = fabs(aa - jj * r2); if (w > w0) w0 = w;
            double lv = (ck[i].j == jref) ? lnr : ck[i].ln;
            double b = fabs(lv - lnr * (1.0 - jj / tt)); if (b > b0) b0 = b;
        }
        /* the start point itself (index jref, inside a halving run) */
        { double b = 0; if (b > b0) b0 = b; }
        tail_tau = tt; tail_a = at; tail_dw = w0; tail_db = b0;
        long e; double d = mpz_get_d_2exp(&e, hmax2); tail_log2h = log2(d) + (double)e;
    }
    mpz_clear(hmax2);
    mpz_clears(x, y, t, NULL);
}

/* T^m(n) by direct iteration (for identity checks) */
static void Tpow(mpz_t out, const mpz_t n, long m) {
    mpz_set(out, n);
    for (long i = 0; i < m; i++) {
        if (mpz_odd_p(out)) { mpz_mul_ui(out, out, 3); mpz_add_ui(out, out, 1); }
        mpz_tdiv_q_2exp(out, out, 1);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s FAM K\n", argv[0]); return 1; }
    const char *fam = argv[1]; long K = atol(argv[2]);
    mpz_t n, hmax, p3, r3, t, u; mpz_inits(n, hmax, p3, r3, t, u, NULL);
    long tau, a, sigma, sigma2; double log2h, dw, db;
    for (long k = 2; k <= K; k++) {
        if (!strcmp(fam, "M")) {
            mpz_ui_pow_ui(n, 2, k); mpz_sub_ui(n, n, 1);
            mpz_ui_pow_ui(p3, 3, k); mpz_sub_ui(p3, p3, 1);          /* 3^k - 1 */
            mpz_tdiv_q_2exp(r3, p3, 1);                               /* (3^k-1)/2 */
            orbit(n, &tau, &a, &sigma, &log2h, hmax, &dw, &db, r3, k + 1, &sigma2);
            int idk = -1, idk1 = -1;
            if (k <= 4000) {                                          /* direct identity checks */
                Tpow(t, n, k); idk = (mpz_cmp(t, p3) == 0);
                Tpow(u, t, 1); idk1 = (mpz_cmp(u, r3) == 0);
            }
            long v2 = mpz_scan1(p3, 0);
            printf("M %ld %ld %ld %ld %.4f %.4f %.4f %d %d %d %ld %ld %.0f %.0f %.4f %.4f %.4f\n", k, tau, a, sigma, log2h,
                   dw, db, idk, idk1, mpz_cmp(hmax, p3) == 0, v2, sigma2, tail_tau, tail_a, tail_dw, tail_db, tail_log2h);
        } else if (!strcmp(fam, "W")) {
            if (k % 2 == 0) continue;
            mpz_ui_pow_ui(n, 2, k); mpz_add_ui(n, n, 1); mpz_divexact_ui(n, n, 3);
            orbit(n, &tau, &a, &sigma, &log2h, hmax, &dw, &db, NULL, 0, NULL);
            int id1 = -1, idk = -1;
            if (k <= 4000) {
                Tpow(t, n, 1); mpz_ui_pow_ui(u, 2, k - 1); mpz_add_ui(u, u, 1); id1 = (mpz_cmp(t, u) == 0);
                Tpow(t, n, k); mpz_ui_pow_ui(u, 3, (k - 1) / 2); mpz_add_ui(u, u, 1); idk = (mpz_cmp(t, u) == 0);
            }
            printf("W %ld %ld %ld %ld %.4f %.4f %.4f %d %d\n", k, tau, a, sigma, log2h, dw, db, id1, idk);
        } else if (!strcmp(fam, "R3")) {
            mpz_ui_pow_ui(n, 3, k); mpz_sub_ui(n, n, 1); mpz_tdiv_q_2exp(n, n, 1);
            orbit(n, &tau, &a, &sigma, &log2h, hmax, &dw, &db, NULL, 0, NULL);
            printf("R3 %ld %ld %ld %ld %.4f %.4f %.4f\n", k, tau, a, sigma, log2h, dw, db);
        } else if (!strcmp(fam, "R10")) {
            mpz_ui_pow_ui(n, 10, k); mpz_sub_ui(n, n, 1); mpz_divexact_ui(n, n, 9);
            orbit(n, &tau, &a, &sigma, &log2h, hmax, &dw, &db, NULL, 0, NULL);
            printf("R10 %ld %ld %ld %ld %.4f %.4f %.4f\n", k, tau, a, sigma, log2h, dw, db);
        }
        fflush(stdout);
    }
    return 0;
}
