/* procgen_repunit_20260924_primes.c
 *
 * Consecutive-prime residue transitions (Lemke Oliver--Soundararajan reproduction).
 * Segmented, odd-only sieve of Eratosthenes with a 3*5*7*11*13 presieve pattern.
 * Every consecutive pair (p_n, p_{n+1}) is recorded in a 120 x 120 table indexed by
 * (p_n mod 120, p_{n+1} mod 120).  Every modulus q | 120 (q = 3, 4, 5, 8, 10, 12, ...) is
 * an aggregation of this table; a pair belongs to the mod-q count iff both residues are
 * reduced mod q, which automatically excludes the pairs whose first or second prime divides q.
 *
 * Snapshots (cumulative tables) are written
 *   - by x:     pairs with p_n <= x          (the LO-S definition of pi(x; q, (a,b)))
 *   - by index: pairs with n <= N            (for "the first N primes" tables)
 *
 * usage: primes XMAX OUTFILE
 *   x-snapshots at 10^(k/4) for k = 16..4*log10(XMAX), plus XMAX itself;
 *   index snapshots at 10^j, 10^j+2, 10^j+3, 10^j+4 (j = 6..9; 10^7 without +2).
 * Each SNAP line also records the last pair counted (p_n, p_{n+1}), so that the variant
 * "both primes <= x" (= this count minus that boundary pair) can be formed exactly.
 * Memory: ~ 0.5 MB (segment 256 KiB + base primes + 120x120 table).
 * Session collatz-procgen-20260922, lane procgen_repunit_20260924.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define SEG (1u << 18)          /* odd numbers per segment */
#define PAT 15015u              /* 3*5*7*11*13 */

static uint64_t tab[120][120];
static FILE *fo;

static void dump(const char *kind, double xval, uint64_t pix, uint64_t npairs, uint64_t plast, uint64_t pnext) {
    /* plast = p_n (first prime of the last pair counted), pnext = p_{n+1} */
    fprintf(fo, "SNAP %s %.0f %llu %llu %llu %llu\n", kind, xval, (unsigned long long)pix, (unsigned long long)npairs,
            (unsigned long long)plast, (unsigned long long)pnext);
    for (int a = 0; a < 120; a++)
        for (int b = 0; b < 120; b++)
            if (tab[a][b]) fprintf(fo, "T %d %d %llu\n", a, b, (unsigned long long)tab[a][b]);
    fprintf(fo, "END\n");
    fflush(fo);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s XMAX OUTFILE\n", argv[0]); return 1; }
    uint64_t XMAX = strtoull(argv[1], NULL, 10);
    fo = fopen(argv[2], "w");
    if (!fo) { perror("fopen"); return 1; }
    clock_t t0 = clock();

    /* x snapshots */
    double xs[64]; int nxs = 0;
    for (int k = 16; ; k++) {
        double x = pow(10.0, k / 4.0);
        if (x > (double)XMAX * 1.0000001) break;
        xs[nxs++] = floor(x + 0.5);
    }
    if (nxs == 0 || xs[nxs - 1] < (double)XMAX) xs[nxs++] = (double)XMAX;
    int ix = 0;
    /* index snapshots */
    uint64_t Ns[] = {1000000ULL, 1000002ULL, 1000003ULL, 1000004ULL, 10000000ULL, 10000003ULL, 10000004ULL,
                     100000000ULL, 100000002ULL, 100000003ULL, 100000004ULL,
                     1000000000ULL, 1000000002ULL, 1000000003ULL, 1000000004ULL, 0};
    int iN = 0;

    uint64_t LIM = XMAX + (1u << 22);            /* sieve a little beyond XMAX to find p_{n+1} */
    uint64_t R = (uint64_t)sqrt((double)LIM) + 2;
    /* base primes up to R (simple sieve) */
    char *small = calloc(R + 1, 1);
    uint32_t *bp = malloc(sizeof(uint32_t) * (R / 2 + 16));
    uint32_t nbp = 0;
    for (uint64_t i = 2; i <= R; i++) {
        if (!small[i]) {
            if (i >= 17) bp[nbp++] = (uint32_t)i;
            for (uint64_t j = i * i; j <= R; j += i) small[j] = 1;
        }
    }
    free(small);
    uint64_t *off = malloc(sizeof(uint64_t) * (nbp + 1));
    for (uint32_t t = 0; t < nbp; t++) off[t] = ((uint64_t)bp[t] * bp[t] - 1) / 2;   /* index of p^2 */

    /* presieve pattern over odd indices i <-> n = 2i+1 */
    unsigned char *pat = malloc(PAT + SEG);
    for (uint32_t i = 0; i < PAT + SEG; i++) {
        uint64_t n = 2ULL * i + 1;
        pat[i] = (n % 3 == 0 || n % 5 == 0 || n % 7 == 0 || n % 11 == 0 || n % 13 == 0);
    }
    unsigned char *seg = malloc(SEG);

    uint64_t npairs = 0, pix = 0;   /* pix = number of primes processed so far */
    int rp = -1;                    /* residue of previous prime mod 120 */
    uint64_t prevp = 0;             /* previous prime */
    int done = 0;

    /* process one prime p (in increasing order) */
#define PROCESS(pval) do {                                                         \
        uint64_t _p = (pval);                                                      \
        int _r = (int)(_p % 120u);                                                 \
        if (rp >= 0) {                                                             \
            tab[rp][_r]++; npairs++;                                               \
            while (ix < nxs && (double)_p > xs[ix]) {                              \
                dump("X", xs[ix], pix, npairs, prevp, _p); ix++;                  \
            }                                                                      \
            while (Ns[iN] && npairs == Ns[iN]) { dump("N", (double)Ns[iN], pix, npairs, prevp, _p); iN++; } \
            if (_p > XMAX) { done = 1; }                                           \
        }                                                                          \
        rp = _r; prevp = _p; pix++;                                                \
    } while (0)

    uint64_t smallp[] = {2, 3, 5, 7, 11, 13};
    for (int t = 0; t < 6; t++) PROCESS(smallp[t]);

    for (uint64_t L0 = 0; !done; L0 += SEG) {          /* L0 = first odd index of the segment */
        uint64_t phase = L0 % PAT;
        memcpy(seg, pat + phase, SEG);
        for (uint32_t t = 0; t < nbp; t++) {
            uint64_t j = off[t];
            uint32_t p = bp[t];
            for (; j < SEG; j += p) seg[j] = 1;
            off[t] = j - SEG;
        }
        uint32_t jstart = (L0 == 0) ? 7 : 0;           /* skip n = 1..13 (handled above) */
        uint64_t lo = 2 * L0 + 1;
        for (uint32_t j = jstart; j < SEG && !done; j++) {
            if (!seg[j]) PROCESS(lo + 2ULL * j);
        }
        if (2 * (L0 + SEG) + 1 > LIM && !done) { fprintf(stderr, "ran past LIM\n"); return 2; }
    }
    fprintf(fo, "DONE xmax=%llu primes_processed=%llu pairs=%llu cpu_s=%.1f\n",
            (unsigned long long)XMAX, (unsigned long long)pix, (unsigned long long)npairs,
            (double)(clock() - t0) / CLOCKS_PER_SEC);
    fclose(fo);
    return 0;
}
