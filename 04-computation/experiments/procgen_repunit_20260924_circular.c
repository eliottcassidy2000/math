/* procgen_repunit_20260924_circular.c
 *
 * Circular and permutable primes (base 10).
 *
 *   circular sieve LIMIT
 *       odd-only bitset sieve to LIMIT (LIMIT = 10^9: 62.5 MB).  For every prime p < LIMIT:
 *         - "eligible" = at least two digits, all in {1,3,7,9}; counted per length;
 *         - circular = every rotation prime (printed "CIRC p");
 *         - rotation orbits with exactly one composite rotation ("near misses") are printed once,
 *           from their smallest prime member ("NEAR k p composite");
 *         - permutable = every permutation of the digits prime (checked for the circular ones; "PERM p").
 *       Summary lines "LEN k eligible circular orbits near" per length.
 *   circular necklace KMIN KMAX
 *       all rotation orbits (necklaces, FKM algorithm) of strings over {1,3,7,9} of length KMIN..KMAX
 *       (KMAX <= 19, so every number is < 2^64); orbits with digit sum = 0 mod 3 are skipped (all
 *       rotations divisible by 3); rotations are filtered by trial division and a base-2 strong
 *       probable-prime test until two composites are found; survivors are certified by a
 *       deterministic Miller--Rabin (12 prime bases, valid below 3.3e24).
 *       Prints circular orbits ("NCIRC"), near-miss orbits ("NNEAR", printed for k <= 12) and a
 *       summary "NLEN k necklaces tested circular near".
 * Session collatz-procgen-20260922, lane procgen_repunit_20260924.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned __int128 u128;

/* ---------------------------------------------------------------- deterministic Miller-Rabin */
static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) { return (uint64_t)((u128)a * b % m); }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) {
    uint64_t r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
static int is_prime_u64(uint64_t n) {
    static const uint64_t sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    if (n < 2) return 0;
    for (int i = 0; i < 25; i++) { if (n == sp[i]) return 1; if (n % sp[i] == 0) return 0; }
    uint64_t d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    static const uint64_t bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int i = 0; i < 12; i++) {
        uint64_t x = powmod(bases[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { ok = 1; break; } }
        if (!ok) return 0;
    }
    return 1;
}

/* quick filter: trial division, then one strong probable-prime test to base 2.  A number failing it is
 * certainly composite; numbers passing it are re-checked with is_prime_u64 before being counted prime. */
static int quick_prime(uint64_t n) {
    static const uint64_t sp[] = {7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    for (int i = 0; i < 22; i++) { if (n == sp[i]) return 1; if (n % sp[i] == 0) return 0; }
    uint64_t d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    uint64_t x = powmod(2, d, n);
    if (x == 1 || x == n - 1) return 1;
    for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) return 1; }
    return 0;
}

/* ---------------------------------------------------------------- sieve mode */
static uint64_t *bits;          /* bit i <-> odd number 2i+1 is composite */
static uint64_t LIMIT;
static inline int isp(uint64_t n) {
    if (n < 2) return 0;
    if (n == 2) return 1;
    if (!(n & 1)) return 0;
    if (n >= LIMIT) return is_prime_u64(n);
    uint64_t i = n >> 1;
    return !((bits[i >> 6] >> (i & 63)) & 1);
}
static int ndig(uint64_t n) { int k = 0; do { k++; n /= 10; } while (n); return k; }
static uint64_t p10[20];

static int all_perms_prime(char *d, int k) {
    /* next_permutation over sorted digits */
    char s[24]; memcpy(s, d, k); s[k] = 0;
    for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) if (s[j] < s[i]) { char t = s[i]; s[i] = s[j]; s[j] = t; }
    for (;;) {
        uint64_t v = 0; for (int i = 0; i < k; i++) v = v * 10 + (s[i] - '0');
        if (!isp(v)) return 0;
        int i = k - 2; while (i >= 0 && s[i] >= s[i + 1]) i--;
        if (i < 0) return 1;
        int j = k - 1; while (s[j] <= s[i]) j--;
        char t = s[i]; s[i] = s[j]; s[j] = t;
        for (int a = i + 1, b = k - 1; a < b; a++, b--) { t = s[a]; s[a] = s[b]; s[b] = t; }
    }
}

static int mode_sieve(uint64_t lim) {
    LIMIT = lim;
    uint64_t nb = lim / 2 + 1;
    bits = calloc(nb / 64 + 1, 8);
    for (uint64_t i = 1; ; i++) {                     /* odd p = 2i+1 */
        uint64_t p = 2 * i + 1;
        if (p * p >= lim) break;
        if ((bits[i >> 6] >> (i & 63)) & 1) continue;
        for (uint64_t j = (p * p) >> 1; j < nb; j += p) bits[j >> 6] |= 1ULL << (j & 63);
    }
    bits[0] |= 1;                                     /* 1 is not prime */
    long elig[20] = {0}, circ[20] = {0}, orb[20] = {0}, nearc[20] = {0};
    uint64_t primes_seen = 0;
    for (uint64_t n = 2; n < lim; n = (n == 2) ? 3 : n + 2) {
        if (!isp(n)) continue;
        primes_seen++;
        int k = ndig(n);
        if (k == 1) { printf("CIRC %llu\n", (unsigned long long)n); circ[1]++; orb[1]++; printf("PERM %llu\n", (unsigned long long)n); continue; }
        uint64_t t = n; int ok = 1;
        while (t) { int d = t % 10; if (d != 1 && d != 3 && d != 7 && d != 9) { ok = 0; break; } t /= 10; }
        if (!ok) continue;
        elig[k]++;
        uint64_t r = n, minprime = n, composite = 0; int ncomp = 0;
        for (int j = 1; j < k; j++) {
            uint64_t lead = r / p10[k - 1];
            r = (r - lead * p10[k - 1]) * 10 + lead;
            if (isp(r)) { if (r < minprime) minprime = r; }
            else { ncomp++; composite = r; }
        }
        if (ncomp == 0) {
            circ[k]++;
            if (minprime == n) orb[k]++;
            printf("CIRC %llu\n", (unsigned long long)n);
            char d[24]; sprintf(d, "%llu", (unsigned long long)n);
            if (all_perms_prime(d, k)) printf("PERM %llu\n", (unsigned long long)n);
        } else if (ncomp == 1 && minprime == n) {
            nearc[k]++;
            printf("NEAR %d %llu %llu\n", k, (unsigned long long)n, (unsigned long long)composite);
        }
    }
    for (int k = 1; k <= ndig(lim - 1); k++)
        printf("LEN %d %ld %ld %ld %ld\n", k, elig[k], circ[k], orb[k], nearc[k]);
    printf("PRIMES_BELOW %llu %llu\n", (unsigned long long)lim, (unsigned long long)primes_seen);
    return 0;
}

/* ---------------------------------------------------------------- necklace mode */
static int K;                         /* current length */
static int a[24];
static const int DIG[4] = {1, 3, 7, 9};
static long nneck, ntested, ncirc, nnear;

static void visit(void) {
    nneck++;
    int s3 = 0; uint64_t v = 0;
    for (int i = 1; i <= K; i++) { s3 += DIG[a[i]]; v = v * 10 + DIG[a[i]]; }
    if (s3 % 3 == 0) return;
    ntested++;
    int ncomp = 0, npass = 0; uint64_t r = v, comp = 0, pass[24];
    for (int j = 0; j < K; j++) {
        if (j) { uint64_t lead = r / p10[K - 1]; r = (r - lead * p10[K - 1]) * 10 + lead; }
        if (!quick_prime(r)) { ncomp++; comp = r; if (ncomp >= 2) return; }
        else pass[npass++] = r;
    }
    for (int j = 0; j < npass; j++)                       /* certify the survivors */
        if (!is_prime_u64(pass[j])) { ncomp++; comp = pass[j]; if (ncomp >= 2) return; }
    if (ncomp == 0) { ncirc++; printf("NCIRC %d %llu\n", K, (unsigned long long)v); }
    else { nnear++; if (K <= 12) printf("NNEAR %d %llu %llu\n", K, (unsigned long long)v, (unsigned long long)comp); }
}
static void gen(int t, int p) {
    if (t > K) { if (K % p == 0) visit(); return; }
    a[t] = a[t - p]; gen(t + 1, p);
    for (int j = a[t - p] + 1; j < 4; j++) { a[t] = j; gen(t + 1, t); }
}
static int mode_necklace(int kmin, int kmax) {
    for (K = kmin; K <= kmax; K++) {
        nneck = ntested = ncirc = nnear = 0;
        a[0] = 0; gen(1, 1);
        printf("NLEN %d %ld %ld %ld %ld\n", K, nneck, ntested, ncirc, nnear);
        fflush(stdout);
    }
    return 0;
}

int main(int argc, char **argv) {
    p10[0] = 1; for (int i = 1; i < 20; i++) p10[i] = p10[i - 1] * 10;
    if (argc >= 3 && !strcmp(argv[1], "sieve")) return mode_sieve(strtoull(argv[2], 0, 10));
    if (argc >= 4 && !strcmp(argv[1], "necklace")) return mode_necklace(atoi(argv[2]), atoi(argv[3]));
    fprintf(stderr, "usage\n");
    return 1;
}
