/* procgen_repunit_20260924_collatz.c
 *
 * Residue-transition tallies for Collatz orbits.
 *   U(n) = odd part of 3n+1   (Syracuse map on odd n; "consecutive odd terms")
 *   T(n) = n/2 (n even), (3n+1)/2 (n odd)   (shortcut map)
 * One-step transitions are tallied mod 360 = 8*9*5 (every modulus dividing 360 -- 3, 5, 8, 9,
 * 10, 40, 72, 360 -- is an aggregation), two-step transitions mod 40 (Markov checks).
 *
 * usage:
 *   collatz bulk NSTARTS SEED STEPS_U STEPS_T OUTFILE
 *       random starts n in [2^100, 2^101) (splitmix64; odd for U); at most STEPS_U Syracuse steps
 *       (U tables) and STEPS_T shortcut steps (T tables) per start.  "terras" tables hold the
 *       transitions determined by the 100 random low bits (exactly Haar-distributed), "post_high" /
 *       "post_low" the later ones (deterministic Collatz dynamics) with current value >= 2^40 / < 2^40.
 *       Orbits stop below 10^6; overflow guard 2^124.
 *   collatz vhist NSTARTS SEED
 *       histogram of v = v_2(3n+1) over the Terras-regime Syracuse steps of NSTARTS random starts
 *       (an exact i.i.d. Geom(1/2) sample; P(v even) = 1/3 is the mod-3 law of consecutive odd terms).
 *   collatz full NMAX THRESH OUTFILE
 *       every start n <= NMAX followed to 1 (U: odd n; T: all n); the fixed/loop transitions at 1
 *       are not counted.  Tables "all" count every transition, tables "big" only those with the
 *       current value > THRESH.
 * Memory: a few MB.  Session collatz-procgen-20260922, lane procgen_repunit_20260924.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned __int128 u128;
#define M1 360
#define M2 40

static uint64_t sm_state;
static uint64_t splitmix64(void) {
    uint64_t z = (sm_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static inline int ctz128(u128 m) {
    uint64_t lo = (uint64_t)m;
    if (lo) return __builtin_ctzll(lo);
    return 64 + __builtin_ctzll((uint64_t)(m >> 64));
}
static inline unsigned mod128(u128 n, unsigned m) { return (unsigned)(n % m); }

typedef struct { uint64_t one[M1][M1]; uint64_t two[M2][M2][M2]; uint64_t steps, cut_low, cut_high; } Tally;

static void dump(FILE *fo, const char *name, Tally *t) {
    fprintf(fo, "TABLE %s steps=%llu cut_low=%llu cut_high=%llu\n", name,
            (unsigned long long)t->steps, (unsigned long long)t->cut_low, (unsigned long long)t->cut_high);
    for (int a = 0; a < M1; a++)
        for (int b = 0; b < M1; b++)
            if (t->one[a][b]) fprintf(fo, "O %d %d %llu\n", a, b, (unsigned long long)t->one[a][b]);
    for (int a = 0; a < M2; a++)
        for (int b = 0; b < M2; b++)
            for (int c = 0; c < M2; c++)
                if (t->two[a][b][c]) fprintf(fo, "W %d %d %d %llu\n", a, b, c, (unsigned long long)t->two[a][b][c]);
    fprintf(fo, "END\n");
}

static inline u128 Ustep(u128 n) { u128 m = 3 * n + 1; return m >> ctz128(m); }
static inline u128 Tstep(u128 n) { return (n & 1) ? (3 * n + 1) >> 1 : n >> 1; }

int main(int argc, char **argv) {
    if (argc < 2) return 1;
    if (!strcmp(argv[1], "bulk")) {
        /* Starts are uniform in [2^100, 2^101): their low 100 bits are uniform, so by the Terras bijection
         * the first 100 T-parity bits are i.i.d. fair.  A U-transition is recorded in the "terras" table while
         * the halvings used so far are <= 60 (the transition then needs at most ~64 bits except with
         * probability 2^-36), a T-transition while t <= 90; later transitions go to the "post" tables
         * (actual Collatz dynamics beyond the random bits), split at 2^40: "post_high" (current value >= 2^40)
         * and "post_low" (10^6 <= value < 2^40, where the 10^7 sampled orbits funnel through few integers and
         * merge, so transitions are counted with multiplicity).  Orbits stop below 10^6. */
        uint64_t NST = strtoull(argv[2], 0, 10);
        sm_state = strtoull(argv[3], 0, 10);
        int SU = atoi(argv[4]), ST = atoi(argv[5]);
        FILE *fo = fopen(argv[6], "w");
        Tally *ut = calloc(1, sizeof(Tally)), *uh = calloc(1, sizeof(Tally)), *ul = calloc(1, sizeof(Tally));
        Tally *tt = calloc(1, sizeof(Tally)), *th = calloc(1, sizeof(Tally)), *tl0 = calloc(1, sizeof(Tally));
        const u128 LOW = 1000000, HIGH = ((u128)1) << 124, MID = ((u128)1) << 40;
        /* merging diagnostic: U-transitions from odd values in the band [10^6, 10^8), every visit (ubA) and
         * only the first visit of each value (ubD, a bitset over the odd values of the band: 6.25 MB) */
        const uint64_t B0 = 1000000ULL, B1 = 100000000ULL;
        uint64_t *seen = calloc((B1 - B0) / 2 / 64 + 1, 8);
        Tally *ubA = calloc(1, sizeof(Tally)), *ubD = calloc(1, sizeof(Tally));
        for (uint64_t s = 0; s < NST; s++) {
            u128 n = (((u128)(splitmix64() >> 28)) << 64) | (u128)splitmix64();   /* 100 random bits */
            n |= ((u128)1) << 100; n |= 1;
            int H = 0, have_prev = 0; unsigned qp = 0;
            for (int k = 0; k < SU; k++) {
                if (n < LOW) { ut->cut_low++; break; }
                if (n > HIGH) { ut->cut_high++; break; }
                u128 m = 3 * n + 1; int v = ctz128(m); u128 n2 = m >> v;
                Tally *tl = (H <= 60) ? ut : (n >= MID ? uh : ul);
                unsigned a = mod128(n, M1), b = mod128(n2, M1), qa = mod128(n, M2), qb = mod128(n2, M2);
                tl->one[a][b]++; tl->steps++;
                if (n < B1) {                                    /* n >= LOW = B0 here */
                    uint64_t i = ((uint64_t)n - B0) >> 1;
                    ubA->one[a][b]++; ubA->steps++;
                    if (!((seen[i >> 6] >> (i & 63)) & 1)) { seen[i >> 6] |= 1ULL << (i & 63); ubD->one[a][b]++; ubD->steps++; }
                }
                if (have_prev) tl->two[qp][qa][qb]++;
                qp = qa; have_prev = 1; H += v; n = n2;
            }
            n = (((u128)(splitmix64() >> 28)) << 64) | (u128)splitmix64();
            n |= ((u128)1) << 100;
            have_prev = 0;
            for (int k = 0; k < ST; k++) {
                if (n < LOW) { tt->cut_low++; break; }
                if (n > HIGH) { tt->cut_high++; break; }
                u128 n2 = Tstep(n);
                Tally *tl = (k <= 90) ? tt : (n >= MID ? th : tl0);
                unsigned a = mod128(n, M1), b = mod128(n2, M1), qa = mod128(n, M2), qb = mod128(n2, M2);
                tl->one[a][b]++; tl->steps++;
                if (have_prev) tl->two[qp][qa][qb]++;
                qp = qa; have_prev = 1; n = n2;
            }
        }
        dump(fo, "U_terras", ut); dump(fo, "U_post_high", uh); dump(fo, "U_post_low", ul);
        dump(fo, "U_band_all", ubA); dump(fo, "U_band_distinct", ubD);
        dump(fo, "T_terras", tt); dump(fo, "T_post_high", th); dump(fo, "T_post_low", tl0);
        fclose(fo);
        return 0;
    }
    if (!strcmp(argv[1], "vhist")) {
        /* v = v_2(3n+1) along Terras-regime Syracuse steps (halvings so far <= 60) of random starts in [2^100, 2^101) */
        uint64_t NST = strtoull(argv[2], 0, 10);
        sm_state = strtoull(argv[3], 0, 10);
        uint64_t cnt[128] = {0}, N = 0;
        for (uint64_t s = 0; s < NST; s++) {
            u128 n = (((u128)(splitmix64() >> 28)) << 64) | (u128)splitmix64();
            n |= ((u128)1) << 100; n |= 1;
            int H = 0;
            while (H <= 60) {
                u128 m = 3 * n + 1; int v = ctz128(m);
                cnt[v < 127 ? v : 127]++; N++;
                H += v; n = m >> v;
            }
        }
        printf("VHIST N=%llu", (unsigned long long)N);
        for (int j = 1; j <= 40; j++) printf(" %llu", (unsigned long long)cnt[j]);
        printf("\n");
        return 0;
    }
    if (!strcmp(argv[1], "full")) {
        uint64_t NMAX = strtoull(argv[2], 0, 10), TH = strtoull(argv[3], 0, 10);
        FILE *fo = fopen(argv[4], "w");
        Tally *ua = calloc(1, sizeof(Tally)), *ub = calloc(1, sizeof(Tally));
        Tally *ta = calloc(1, sizeof(Tally)), *tb = calloc(1, sizeof(Tally));
        for (uint64_t s = 1; s <= NMAX; s++) {
            if (s & 1) {                                  /* U orbit of odd s, until 1 */
                u128 n = s; int have_prev = 0; unsigned qp = 0;
                while (n != 1) {
                    u128 m = Ustep(n);
                    unsigned a = mod128(n, M1), b = mod128(m, M1), qa = mod128(n, M2), qb = mod128(m, M2);
                    ua->one[a][b]++; ua->steps++;
                    if (have_prev) ua->two[qp][qa][qb]++;
                    if (n > TH) { ub->one[a][b]++; ub->steps++; }
                    qp = qa; have_prev = 1; n = m;
                }
            }
            {                                             /* T orbit of s, until 1 */
                u128 n = s; int have_prev = 0; unsigned qp = 0;
                while (n != 1) {
                    u128 m = Tstep(n);
                    unsigned a = mod128(n, M1), b = mod128(m, M1), qa = mod128(n, M2), qb = mod128(m, M2);
                    ta->one[a][b]++; ta->steps++;
                    if (have_prev) ta->two[qp][qa][qb]++;
                    if (n > TH) { tb->one[a][b]++; tb->steps++; }
                    qp = qa; have_prev = 1; n = m;
                }
            }
        }
        dump(fo, "U_full_all", ua); dump(fo, "U_full_big", ub);
        dump(fo, "T_full_all", ta); dump(fo, "T_full_big", tb);
        fclose(fo);
        return 0;
    }
    return 1;
}
