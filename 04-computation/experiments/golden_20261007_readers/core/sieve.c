/* sieve.c -- exact Collatz (Terras map T) residue sieves, lane "core" (golden session 2026-10-07).

   A class is a parity word w of length t (<-> n mod 2^t, Terras bijection), optionally refined by n mod 3^J.
   On the class, x_s = T^s(n) = (3^{a_s} n + c_s)/2^s for s <= t (affine in n).

   Certificates (each is "for all large n in the class, n is not the minimum of its grand orbit"):
     DESC   : 3^{a_s} < 2^s for some 1 <= s <= t                         (Terras descent)
     MERGE  : T^s(n-d) = T^s(n) uniformly on the class, some s <= t, d>=1  (pair chain from (0,-d) absorbed;
              THM-4581 chain; only d < 2^{s-a_s} can occur, so d <= DMAX is exhaustive when DMAX >= 2^{max(s-a_s)})
     BRANCH : a grand-orbit element m = (2^i x_s - c(u))/3^b reached by a backward word u (i steps, b odd)
              from some x_s, s <= t, with b <= a_s + J (decided by the class) and slope ratio 3^{a_s-b} 2^{i-s} < 1.
              (Includes 3-adic predecessor certificates such as n = 2 mod 3 and n = 4 mod 9 at s = 0.)
   A backward word with ratio exactly 1 is a translation (i = s, b = a_s), i.e. a MERGE; ratio < 1 is BRANCH.
   So DESC+MERGE(all d)+BRANCH is the maximal sieve whose certificates are decided by (n mod 2^t, n mod 3^J).

   usage: sieve KMAX DMAX BRANCH J
     DMAX   = number of translation chains d = 1..DMAX (0 = no merges)
     BRANCH = 0/1
     J      = 3-adic refinement depth (0..3)
   output per depth t: number of uncertified classes (out of 2^t 3^J) and the fraction. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef __int128 i128;
typedef unsigned __int128 u128;

static int KMAX, DMAX, BRANCH, J, IMAX = 0, MERGE_AFTER = 0;
static uint64_t c_mergeonly[128];
#define B 39                     /* 3-adic precision digits; 3^39 < 2^63 */
static uint64_t P3[48];
static i128 P3L[100];
static uint64_t INV2;            /* 1/2 mod 3^B */
static uint64_t INV2POW[128];    /* 2^{-t} mod 3^B */
static int NSUB;                 /* 3^J */

static uint64_t alive_cnt[128], c_desc[128], c_merge[128], c_branch[128], nodes[128];
static uint64_t bnodes = 0;
static int overflow = 0;

typedef struct { int k; i128 N; } chain_t;
static chain_t *CH;              /* CH[t*DMAX + (d-1)] */

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) { return (uint64_t)(((u128)a * b) % m); }

/* one Terras step of the pair chain u = 3^k v + e, e = N / 3^max(0,-k); beta = parity of v */
static inline chain_t cstep(chain_t s, int beta) {
    chain_t r; int k = s.k; i128 N = s.N;
    int Nodd = (int)(N & 1);
    if (k >= 0) {
        if (!Nodd) {
            if (!beta) { r.k = k; r.N = N / 2; }
            else { r.k = k; r.N = (3 * N + 1 - P3L[k]) / 2; }
        } else {
            if (!beta) { r.k = k + 1; r.N = (3 * N + 1) / 2; }
            else if (k >= 1) { r.k = k - 1; r.N = (N - P3L[k - 1]) / 2; }
            else { r.k = -1; r.N = (3 * N - 1) / 2; }
        }
    } else {
        int a = -k;
        if (!Nodd) {
            if (!beta) { r.k = k; r.N = N / 2; }
            else { r.k = k; r.N = (3 * N + P3L[a] - 1) / 2; }
        } else {
            if (!beta) { r.k = k + 1; r.N = (N + P3L[a - 1]) / 2; }
            else { r.k = k - 1; r.N = (3 * N - 1) / 2; }
        }
    }
    i128 lim = ((i128)1) << 120;
    if (r.N > lim || r.N < -lim) overflow = 1;
    return r;
}

/* is 2^e2 3^e3 < 1 ? exact */
static int lt1(int e2, int e3) {
    if (e2 <= 0 && e3 <= 0) return !(e2 == 0 && e3 == 0);
    if (e2 >= 0 && e3 >= 0) return 0;
    if (e2 < 0) { /* 3^e3 < 2^-e2 */
        if (-e2 >= 127) return 1;
        u128 p2 = ((u128)1) << (-e2);
        u128 p3 = 1; for (int i = 0; i < e3; i++) { p3 *= 3; if (p3 >= p2) return 0; }
        return p3 < p2;
    } else { /* 2^e2 < 3^-e3 */
        if (e2 >= 127) return 0;
        u128 p2 = ((u128)1) << e2;
        u128 p3 = 1; for (int i = 0; i < -e3; i++) { p3 *= 3; if (p3 > p2) return 1; }
        return p2 < p3;
    }
}

/* backward search from a node whose value is z mod 3^p; ratio (node / n) = 2^e2 3^e3.
   returns 1 if some backward descendant (including the node itself) has ratio < 1. */
static int bsrch(uint64_t z, int p, int e2, int e3, int dep) {
    bnodes++;
    if (lt1(e2, e3)) return 1;
    if (IMAX && dep >= IMAX) return 0;
    if (!lt1(e2 + p, e3 - p)) return 0;          /* even all-odd continuation cannot get below 1 */
    if (p == 0) return 0;                        /* no 3-adic info: only even steps possible, ratio grows */
    int r = (int)(z % 3);
    if (r == 0) return 0;                        /* multiples of 3 have no odd preimage, ever */
    if (r == 2) {                                /* odd backward step: (2z-1)/3 */
        uint64_t w = (2 * z - 1) / 3;            /* exact since z = 2 mod 3; z < 3^p <= 3^39 so 2z fits */
        w %= P3[p - 1];
        if (bsrch(w, p - 1, e2 + 1, e3 - 1, dep + 1)) return 1;
    }
    /* even backward step: 2z */
    uint64_t w = (2 * z) % P3[p];
    return bsrch(w, p, e2 + 1, e3, dep + 1);
}

static int descent(int t, int a) { /* 3^a < 2^t */
    if (t >= 127) return 1;
    u128 p2 = ((u128)1) << t, p3 = 1;
    for (int i = 0; i < a; i++) { p3 *= 3; if (p3 >= p2) return 0; }
    return p3 < p2;
}

/* node at depth t: word fixed, a = #odd, Y = c_t 2^{-t} mod 3^B (x_t = 3^a n 2^{-t} + Y mod 3^B),
   chains at CH[t], mask of alive 3-adic subclasses (bit rho for n = rho mod 3^J) */
static inline int pc128(u128 m) { return __builtin_popcountll((uint64_t)m) + __builtin_popcountll((uint64_t)(m >> 64)); }
static void dfs(int t, int a, uint64_t Y, u128 mask) {
    nodes[t]++;
    int pc = pc128(mask);
    if (t >= 1 && descent(t, a)) { c_desc[t] += pc; return; }
    if (DMAX > 0 && t >= 1 && !MERGE_AFTER) {
        chain_t *cur = CH + (size_t)t * DMAX;
        for (int d = 0; d < DMAX; d++) if (cur[d].k == 0 && cur[d].N == 0) { c_merge[t] += pc; return; }
    }
    if (BRANCH) {
        int p = a + J; if (p > B) p = B;
        uint64_t mod = P3[p];
        for (int rho = 0; rho < NSUB; rho++) {
            if (!((mask >> rho) & 1)) continue;
            /* x_t mod 3^p = Y + 3^a rho 2^{-t} */
            uint64_t z = Y % mod;
            if (J > 0 && a < B) {
                uint64_t term = mulmod(mulmod(P3[a] % mod, (uint64_t)rho, mod), INV2POW[t] % mod, mod);
                z = (z + term) % mod;
            }
            if (bsrch(z, p, -t, a, 0)) { mask &= ~(((u128)1) << rho); c_branch[t]++; }
        }
        if (!mask) return;
        pc = pc128(mask);
    }
    if (DMAX > 0 && t >= 1 && MERGE_AFTER) {
        chain_t *cur = CH + (size_t)t * DMAX;
        for (int d = 0; d < DMAX; d++) if (cur[d].k == 0 && cur[d].N == 0) { c_mergeonly[t] += pc; return; }
    }
    alive_cnt[t] += pc;
    if (t == KMAX) return;
    for (int beta = 0; beta <= 1; beta++) {
        if (DMAX > 0) {
            chain_t *cur = CH + (size_t)t * DMAX, *nxt = CH + (size_t)(t + 1) * DMAX;
            for (int d = 0; d < DMAX; d++) nxt[d] = cstep(cur[d], beta);
        }
        uint64_t Y2 = beta ? mulmod((3 * Y + 1) % P3[B], INV2, P3[B]) : mulmod(Y, INV2, P3[B]);
        dfs(t + 1, a + beta, Y2, mask);
    }
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: sieve KMAX DMAX BRANCH J [IMAX] [MERGE_AFTER]\n"); return 1; }
    KMAX = atoi(argv[1]); DMAX = atoi(argv[2]); BRANCH = atoi(argv[3]); J = atoi(argv[4]);
    if (argc > 5) IMAX = atoi(argv[5]);
    if (argc > 6) MERGE_AFTER = atoi(argv[6]);
    P3[0] = 1; for (int i = 1; i < 48; i++) P3[i] = (i <= 40) ? P3[i - 1] * 3 : 0;
    P3L[0] = 1; for (int i = 1; i < 80; i++) P3L[i] = P3L[i - 1] * 3;
    INV2 = (P3[B] + 1) / 2;
    INV2POW[0] = 1; for (int i = 1; i < 128; i++) INV2POW[i] = mulmod(INV2POW[i - 1], INV2, P3[B]);
    NSUB = 1; for (int i = 0; i < J; i++) NSUB *= 3;
    if (NSUB > 81) { fprintf(stderr, "J <= 4\n"); return 1; }
    u128 mask = (((u128)1) << NSUB) - 1;
    if (DMAX > 0) {
        CH = (chain_t *)malloc(sizeof(chain_t) * (size_t)(KMAX + 2) * DMAX);
        for (int d = 0; d < DMAX; d++) { CH[d].k = 0; CH[d].N = -(i128)(d + 1); }
    }
    clock_t t0 = clock();
    dfs(0, 0, 0, mask);
    double sec = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("# sieve KMAX=%d DMAX=%d BRANCH=%d J=%d IMAX=%d MERGE_AFTER=%d (classes per depth t: 2^t * 3^J)  time %.1fs  bnodes %llu  overflow %d\n",
           KMAX, DMAX, BRANCH, J, IMAX, MERGE_AFTER, sec, (unsigned long long)bnodes, overflow);
    printf("# t  uncertified  fraction  newly_cert_desc  newly_cert_merge  newly_cert_branch  dfs_nodes  merge_only_after_branch\n");
    for (int t = 0; t <= KMAX; t++) {
        long double tot = ldexpl((long double)NSUB, t);
        printf("%2d %14llu %.10Lf %12llu %12llu %12llu %12llu %8llu\n", t, (unsigned long long)alive_cnt[t],
               (long double)alive_cnt[t] / tot, (unsigned long long)c_desc[t], (unsigned long long)c_merge[t],
               (unsigned long long)c_branch[t], (unsigned long long)nodes[t], (unsigned long long)c_mergeonly[t]);
    }
    return 0;
}
