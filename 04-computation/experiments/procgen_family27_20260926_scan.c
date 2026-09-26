/*
 * procgen_family27_20260926_scan.c
 *
 * Exhaustive scan of 1 <= n <= X for the shortcut Collatz map
 *     T(x) = x/2 (x even),  (3x+1)/2 (x odd),
 * trajectories stopped at the first arrival at 1 (the K-L / OEIS convention).
 *
 * Per n it computes
 *   dT(n)    total stopping time in T-steps            (K-L sigma_infinity)
 *   ones(n)  number of odd steps before reaching 1      (standard delay = dT + ones, A006878)
 *   glT(n)   glide (dropping time) in T-steps           (A060413), glStd = glT + odd steps in it (A217934)
 *   tT(n)    max of the T-trajectory                    (K-L t(n));  tStd = max(n, 2*maxOdd)  (A006885)
 *   win(n)   max_{0 <= j <= floor(log2 n)} T^j(n)       (Terras window maximum)
 *   lab(n)   index of the first "root" (from a fixed root list) met by the trajectory
 * and accumulates histograms per dyadic block b = floor(log2 n).
 *
 * Method: memo arrays for all n < N0 = 2^24 (filled in increasing order); for odd m >= N0 the
 * trajectory is iterated until it falls below N0 (and at least floor(log2 m) steps, for the
 * window), then the memo supplies the rest.  Even n = 2^a m (m odd) are handled in closed form:
 *   dT(2^a m) = a + dT(m), ones same, tT = max(n, tT(m)), win = max(n, win(m)), glide = 1,
 *   lab = first root among 2^a m, ..., 2m, else lab(m).
 * Records are detected in increasing n order; even delay-record candidates are kept in a small
 * heap (an even n = 2k can be a delay record only if d(k) >= running max; see the note).
 *
 * Output: a plain-text report on stdout, parsed by procgen_family27_20260926_run.py.
 * Usage:  scan X   (X <= 2^32)
 * Memory: about 220 MB for N0 = 2^24.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;

#ifndef LOGN0
#define LOGN0 24
#endif
#define N0 (1u << LOGN0)
#define ROOTLIM (1u << 24)
#define NB 34          /* dyadic blocks 0..33 */
#define NRB 512        /* ratio bins (width 1/4) */
#define NRI 130        /* riser index bins: W = 2^(i/2) */
#define MAXROOTS 250
#define NLC 6          /* glide classes for the excursion histogram */
#define NUB 250        /* u bins of width 0.02 on [0,5) */
#define LMAXEX 1200

static uint16_t *dTm, *onesm;
static uint64_t *maxOddm;
static uint8_t *labm;

static uint64_t roots[MAXROOTS];
static int nroots = 0;
static int16_t *rootidx;      /* rootidx[v] for v < ROOTLIM, -1 if not a root */
static uint8_t reach[MAXROOTS][MAXROOTS]; /* reach[i][k] = root k lies on orbit of root i */

/* histograms */
static uint64_t H_all[NB];
static uint64_t H_full[NB][NRI], H_win[NB][NRI];
static uint64_t H_dT[NB][NRB], H_dStd[NB][NRB], H_gl[NB][NRB];
static uint64_t H_lab[NB][MAXROOTS];
static uint64_t H_B27cls[NB][3];
static uint64_t H_glen[NB][1024];     /* exact glide length (T-steps), capped 1023 */
static uint64_t H_dTB[NB][NRB];       /* delay-ratio histogram restricted to 27's branch */
static uint64_t C27_rise[NB], C27_winrise[NB], C27_delay[NB], C27_glide[NB], C27_rho[NB];
static uint64_t EX_hist[NLC][NUB + 1];
static uint64_t EX_cnt[LMAXEX + 1];
static double EX_su[LMAXEX + 1], EX_su2[LMAXEX + 1];
static int inB27[MAXROOTS];

static const double SIGMA = 0.53013833; /* sqrt(p(1-p)) ln 3, p = log_3 2; recomputed in main */
static double sigma_tilt;
static double LN27, LN4616;

static inline uint64_t Tstep(uint64_t x) {
    return (x & 1) ? x + (x >> 1) + 1 : (x >> 1);
}

static void overflow_die(uint64_t n) {
    fprintf(stderr, "overflow risk at n=%llu\n", (unsigned long long)n);
    exit(3);
}

static inline int ilog2_u64(uint64_t n) { return 63 - __builtin_clzll(n); }

/* exact: t >= 2^(i/2) * n ? */
static inline int riser_check(uint64_t t, uint64_t n, int i) {
    if (i <= 0) return t >= n;
    if ((i & 1) == 0) {
        int e = i / 2;
        if (e >= 64) return 0;
        return (t >> e) >= n;
    } else {
        u128 tt = (u128)t * t;
        u128 nn = (u128)n * n;
        if (i >= 128) return 0;
        return (tt >> i) >= nn;
    }
}

static inline int riser_index(uint64_t t, uint64_t n) {
    double q = (double)t / (double)n;
    int e;
    double mant = frexp(q, &e); /* q = mant*2^e, mant in [0.5,1) */
    int i0 = 2 * e - 2 + (mant >= 0.70710678118654752 ? 1 : 0);
    if (i0 < 0) i0 = 0;
    while (i0 + 1 < NRI && riser_check(t, n, i0 + 1)) i0++;
    while (i0 > 0 && !riser_check(t, n, i0)) i0--;
    if (i0 >= NRI) i0 = NRI - 1;
    return i0;
}

/* ---------------- record bookkeeping (increasing n order) ---------------- */
typedef struct { uint64_t pos; uint32_t val; } cand_t;
typedef struct { cand_t a[4096]; int n; } heap_t;
static void hpush(heap_t *h, uint64_t pos, uint32_t val) {
    if (h->n >= 4096) { fprintf(stderr, "heap full\n"); exit(4); }
    int i = h->n++;
    h->a[i].pos = pos; h->a[i].val = val;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (h->a[p].pos <= h->a[i].pos) break;
        cand_t tmp = h->a[p]; h->a[p] = h->a[i]; h->a[i] = tmp; i = p;
    }
}
static cand_t hpop(heap_t *h) {
    cand_t top = h->a[0];
    h->a[0] = h->a[--h->n];
    int i = 0;
    for (;;) {
        int l = 2 * i + 1, r = l + 1, s = i;
        if (l < h->n && h->a[l].pos < h->a[s].pos) s = l;
        if (r < h->n && h->a[r].pos < h->a[s].pos) s = r;
        if (s == i) break;
        cand_t tmp = h->a[s]; h->a[s] = h->a[i]; h->a[i] = tmp; i = s;
    }
    return top;
}

static heap_t hT, hS;
static uint32_t DmaxT = 0, DmaxS = 0;
static uint64_t X;

static void delay_rec_pending(heap_t *h, uint32_t *Dmax, uint64_t upto, const char *tag) {
    while (h->n > 0 && h->a[0].pos < upto) {
        cand_t c = hpop(h);
        if (c.val > *Dmax) {
            *Dmax = c.val;
            printf("REC %s %llu %u\n", tag, (unsigned long long)c.pos, c.val);
        }
        if (c.val >= *Dmax && 2 * c.pos <= X) hpush(h, 2 * c.pos, c.val + 1);
    }
}
static void delay_rec_odd(heap_t *h, uint32_t *Dmax, uint64_t m, uint32_t val, const char *tag) {
    if (val > *Dmax) {
        *Dmax = val;
        printf("REC %s %llu %u\n", tag, (unsigned long long)m, val);
    }
    if (val >= *Dmax && 2 * m <= X) hpush(h, 2 * m, val + 1);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: scan X\n"); return 1; }
    X = strtoull(argv[1], NULL, 10);
    if (X < N0 || X > (1ull << 32)) { fprintf(stderr, "need 2^24 <= X <= 2^32\n"); return 1; }
    {
        double p = log(2.0) / log(3.0);
        sigma_tilt = sqrt(p * (1 - p)) * log(3.0);
        (void)SIGMA;
    }
    LN27 = log(27.0); LN4616 = log(4616.0);

    dTm = (uint16_t *)calloc(N0, sizeof(uint16_t));
    onesm = (uint16_t *)calloc(N0, sizeof(uint16_t));
    maxOddm = (uint64_t *)calloc(N0, sizeof(uint64_t));
    labm = (uint8_t *)calloc(N0, sizeof(uint8_t));
    rootidx = (int16_t *)malloc(ROOTLIM * sizeof(int16_t));
    if (!dTm || !onesm || !maxOddm || !labm || !rootidx) { fprintf(stderr, "alloc\n"); return 2; }
    for (uint32_t v = 0; v < ROOTLIM; v++) rootidx[v] = -1;

    /* ---- roots: 27's T-orbit, the ladder of 2051 = E(3077), the trunk (4^j-1)/3, the ladder of 27 ---- */
    {
        uint64_t x = 27;
        for (;;) { if (rootidx[x] < 0) { rootidx[x] = nroots; roots[nroots++] = x; } if (x == 1) break; x = Tstep(x); }
        uint64_t p = 2051;
        for (int j = 1; j <= 6; j++) { p = 4 * p + 1; if (p < ROOTLIM && rootidx[p] < 0) { rootidx[p] = nroots; roots[nroots++] = p; } }
        uint64_t tr = 1;
        for (int j = 1; j <= 12; j++) { if (tr < ROOTLIM && rootidx[tr] < 0) { rootidx[tr] = nroots; roots[nroots++] = tr; } tr = 4 * tr + 1; }
        p = 27;
        for (int j = 1; j <= 9; j++) { p = 4 * p + 1; if (p < ROOTLIM && rootidx[p] < 0) { rootidx[p] = nroots; roots[nroots++] = p; } }
    }
    for (int i = 0; i < nroots; i++) {
        uint64_t x = roots[i];
        for (;;) {
            if (x < ROOTLIM && rootidx[x] >= 0) reach[i][rootidx[x]] = 1;
            if (x == 1) break;
            x = Tstep(x);
        }
    }
    int idx3077 = rootidx[3077];
    for (int i = 0; i < nroots; i++) inB27[i] = reach[i][idx3077];
    printf("ROOTS %d\n", nroots);
    for (int i = 0; i < nroots; i++) printf("ROOT %d %llu %d\n", i, (unsigned long long)roots[i], inB27[i]);

    /* ---- memo phase: all 1 <= n < N0 in increasing order ---- */
    dTm[1] = 0; onesm[1] = 0; maxOddm[1] = 0; labm[1] = (uint8_t)rootidx[1];
    for (uint64_t n = 2; n < N0; n++) {
        uint64_t x = n, segmax = 0; uint32_t j = 0, c = 0; int r = rootidx[n];
        while (x >= n) {
            if (x & 1) { x = x + (x >> 1) + 1; c++; if (x > segmax) segmax = x; }
            else x >>= 1;
            j++;
            if (r < 0 && x >= n && x < ROOTLIM && rootidx[x] >= 0) r = rootidx[x];
            if (x > 0x5555555555555554ull) overflow_die(n);
        }
        uint32_t d = j + dTm[x], o = c + onesm[x];
        if (d > 65535 || o > 65535) { fprintf(stderr, "u16 overflow\n"); return 5; }
        dTm[n] = (uint16_t)d; onesm[n] = (uint16_t)o;
        maxOddm[n] = segmax > maxOddm[x] ? segmax : maxOddm[x];
        labm[n] = (uint8_t)(r >= 0 ? r : labm[x]);
    }
    fprintf(stderr, "memo done\n");

    /* ---- stream over odd m <= X (all n = 2^a m handled here) ---- */
    uint64_t pathmaxT = 0, pathmaxS = 0; uint32_t glmaxT = 0, glmaxS = 0;
    double gamrec = -1, glratrec = -1, rhorec = -1;
    /* n = 1 */
    printf("REC dT 1 0\nREC dS 1 0\nREC tT 1 1\nREC tS 1 1\n");
    pathmaxT = 1; pathmaxS = 1;
    hpush(&hT, 2, 1); hpush(&hS, 2, 1);
    for (uint64_t m = 1; m <= X; m += 2) {
        uint64_t tTm, winm; uint32_t dT, ones, glT = 0, glC = 0; uint64_t glmax = 0; int lab;
        int w = ilog2_u64(m);
        if (m < N0) {
            dT = dTm[m]; ones = onesm[m]; lab = labm[m];
            tTm = maxOddm[m] > m ? maxOddm[m] : m;
            /* glide + window by direct iteration */
            uint64_t x = m, segmax = 0; uint32_t j = 0, c = 0; winm = m;
            if (m > 1) {
                for (;;) {
                    if (x & 1) { x = x + (x >> 1) + 1; c++; if (x > segmax) segmax = x; if ((int)j < w && x > winm) winm = x; }
                    else x >>= 1;
                    j++;
                    if (!glT && x < m) { glT = j; glC = c; glmax = segmax; }
                    if (glT && (int)j >= w) break;
                }
            }
        } else {
            uint64_t x = m, segmax = 0, x1 = 0; uint32_t j = 0, c = 0; winm = m;
            int r = (m < ROOTLIM) ? rootidx[m] : -1;   /* never a root when N0 = ROOTLIM */
            for (;;) {
                if (x & 1) { x = x + (x >> 1) + 1; c++; if (x > segmax) segmax = x; if ((int)j < w && x > winm) winm = x; }
                else x >>= 1;
                j++;
                if (!glT && x < m) { glT = j; glC = c; glmax = segmax; }
                if (!x1 && r < 0 && x >= N0 && x < ROOTLIM && rootidx[x] >= 0) r = rootidx[x];
                if (!x1 && x < N0) x1 = x;   /* first value below N0: labm[x1] labels the rest */
                if (x < N0 && (int)j >= w) break;
                if (x > 0x5555555555555554ull) overflow_die(m);
            }
            dT = j + dTm[x]; ones = c + onesm[x]; lab = (r >= 0) ? r : labm[x1];
            uint64_t mo = segmax > maxOddm[x] ? segmax : maxOddm[x];
            tTm = mo > m ? mo : m;
            if (tTm > 0x4000000000000000ull) overflow_die(m);
        }
        uint64_t tSm = (tTm > m) ? 2 * tTm : m; /* standard map maximum */

        /* ---- records (increasing order): pending evens < m first ---- */
        delay_rec_pending(&hT, &DmaxT, m, "dT");
        delay_rec_pending(&hS, &DmaxS, m, "dS");
        if (m > 1) {
            if (m == 3) { /* n = 2 precedes m = 3: t(2) = 2, glide(2) = 1 */
                printf("REC tT 2 2\nREC tS 2 2\nREC gT 2 1\nREC gS 2 1\n");
                pathmaxT = 2; pathmaxS = 2; glmaxT = 1; glmaxS = 1;
            }
            delay_rec_odd(&hT, &DmaxT, m, dT, "dT");
            delay_rec_odd(&hS, &DmaxS, m, dT + ones, "dS");
            if (tTm > pathmaxT) { pathmaxT = tTm; printf("REC tT %llu %llu\n", (unsigned long long)m, (unsigned long long)tTm); }
            if (tSm > pathmaxS) { pathmaxS = tSm; printf("REC tS %llu %llu\n", (unsigned long long)m, (unsigned long long)tSm); }
            if (glT > glmaxT) { glmaxT = glT; printf("REC gT %llu %u\n", (unsigned long long)m, glT); }
            if (glT + glC > glmaxS) { glmaxS = glT + glC; printf("REC gS %llu %u\n", (unsigned long long)m, glT + glC); }
            double lnm = log((double)m);
            double gam = dT / lnm;
            if (gam > gamrec) { gamrec = gam; printf("REC gam %llu %.9f\n", (unsigned long long)m, gam); }
            double glr = glT / (lnm / log(2.0));
            if (glr > glratrec) { glratrec = glr; printf("REC glr %llu %.9f\n", (unsigned long long)m, glr); }
            double rho = log((double)tTm) / lnm;
            if (rho > rhorec) { rhorec = rho; printf("REC rho %llu %.9f\n", (unsigned long long)m, rho); }
            /* excursion statistics of the glide */
            if (glT >= 20) {
                double u = log((double)glmax / (double)m) / (sigma_tilt * sqrt((double)glT));
                int lc = glT < 30 ? 0 : glT < 40 ? 1 : glT < 60 ? 2 : glT < 80 ? 3 : glT < 120 ? 4 : 5;
                int ub = (int)(u / 0.02); if (ub > NUB) ub = NUB; if (ub < 0) ub = 0;
                EX_hist[lc][ub]++;
                int L = glT > LMAXEX ? LMAXEX : glT;
                EX_cnt[L]++; EX_su[L] += u; EX_su2[L] += u * u;
            }
        }

        /* ---- histogram contributions of n = 2^a m <= X ---- */
        double lnm = (m > 1) ? log((double)m) : 0.0;
        uint64_t n = m;
        int labn = lab;
        for (int a = 0; ; a++) {
            if (a > 0) {
                if (n > (X >> 1)) break;
                n <<= 1;
            }
            int b = w + a;
            H_all[b]++;
            /* label: first root among 2^a m ... 2m, then lab(m) */
            if (a > 0) {
                labn = lab;
                uint64_t y = n;
                for (int i = a; i >= 1; i--, y >>= 1) {
                    if (y < ROOTLIM && rootidx[y] >= 0) { labn = rootidx[y]; break; }
                }
            }
            H_lab[b][labn]++;
            if (inB27[labn]) H_B27cls[b][n % 3]++;
            uint64_t tT = tTm > n ? tTm : n;
            uint64_t wn = winm > n ? winm : n;
            H_full[b][riser_index(tT, n)]++;
            H_win[b][riser_index(wn, n)]++;
            if ((u128)tT * 27 >= (u128)4616 * n) C27_rise[b]++;
            if ((u128)wn * 27 >= (u128)4616 * n) C27_winrise[b]++;
            if (n >= 2) {
                double lnn = lnm + a * 0.69314718055994531;
                uint32_t dTn = dT + a, dSn = dT + a + ones;
                int k1 = (int)(4.0 * dTn / lnn); if (k1 >= NRB) k1 = NRB - 1;
                int k2 = (int)(4.0 * dSn / lnn); if (k2 >= NRB) k2 = NRB - 1;
                H_dT[b][k1]++; H_dStd[b][k2]++;
                if (inB27[labn]) H_dTB[b][k1]++;
                uint32_t gl = (a > 0) ? 1 : glT;
                H_glen[b][gl > 1023 ? 1023 : gl]++;
                double l2 = lnn / 0.69314718055994531;
                int k3 = (int)(4.0 * gl / l2); if (k3 >= NRB) k3 = NRB - 1;
                H_gl[b][k3]++;
                if (dTn * LN27 >= 70.0 * lnn - 1e-9) C27_delay[b]++;
                if (gl * LN27 >= 59.0 * lnn - 1e-9) C27_glide[b]++;
                if (log((double)tT) * LN27 >= LN4616 * lnn - 1e-12) C27_rho[b]++;
                /* rho > 2 list: tT > n^2 */
                if ((u128)tT > (u128)n * n)
                    printf("RHO2 %llu %llu %d\n", (unsigned long long)n, (unsigned long long)tT, labn);
            }
        }
    }
    /* flush remaining even candidates */
    delay_rec_pending(&hT, &DmaxT, X + 1, "dT");
    delay_rec_pending(&hS, &DmaxS, X + 1, "dS");

    /* ---- dump histograms ---- */
    printf("X %llu\n", (unsigned long long)X);
    printf("SIGMA %.10f\n", sigma_tilt);
    for (int b = 0; b < NB; b++) {
        if (!H_all[b]) continue;
        printf("BLK %d %llu\n", b, (unsigned long long)H_all[b]);
        printf("HFULL %d", b); for (int i = 0; i < NRI; i++) printf(" %llu", (unsigned long long)H_full[b][i]); printf("\n");
        printf("HWIN %d", b); for (int i = 0; i < NRI; i++) printf(" %llu", (unsigned long long)H_win[b][i]); printf("\n");
        printf("HDT %d", b); for (int i = 0; i < NRB; i++) printf(" %llu", (unsigned long long)H_dT[b][i]); printf("\n");
        printf("HDS %d", b); for (int i = 0; i < NRB; i++) printf(" %llu", (unsigned long long)H_dStd[b][i]); printf("\n");
        printf("HGL %d", b); for (int i = 0; i < NRB; i++) printf(" %llu", (unsigned long long)H_gl[b][i]); printf("\n");
        printf("HLAB %d", b); for (int i = 0; i < nroots; i++) printf(" %llu", (unsigned long long)H_lab[b][i]); printf("\n");
        printf("HDTB %d", b); for (int i = 0; i < NRB; i++) printf(" %llu", (unsigned long long)H_dTB[b][i]); printf("\n");
        for (int i = 0; i < 1024; i++) if (H_glen[b][i]) printf("HGLEN %d %d %llu\n", b, i, (unsigned long long)H_glen[b][i]);
        printf("HB27 %d %llu %llu %llu\n", b, (unsigned long long)H_B27cls[b][0], (unsigned long long)H_B27cls[b][1], (unsigned long long)H_B27cls[b][2]);
        printf("C27 %d %llu %llu %llu %llu %llu\n", b, (unsigned long long)C27_rise[b], (unsigned long long)C27_winrise[b],
               (unsigned long long)C27_delay[b], (unsigned long long)C27_glide[b], (unsigned long long)C27_rho[b]);
    }
    for (int lc = 0; lc < NLC; lc++) {
        printf("EXH %d", lc); for (int i = 0; i <= NUB; i++) printf(" %llu", (unsigned long long)EX_hist[lc][i]); printf("\n");
    }
    for (int L = 20; L <= LMAXEX; L++) if (EX_cnt[L]) printf("EXL %d %llu %.10e %.10e\n", L, (unsigned long long)EX_cnt[L], EX_su[L], EX_su2[L]);
    printf("END\n");
    return 0;
}
