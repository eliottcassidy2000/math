/* evenclos.c - THM-5055 base computations (mac-mini-2026-11-09-digraph721).
 *
 * Value-level closure of even-graph realisations by order n. Sets of odd values per order n:
 *   E  d(G) for an even graph G on n vertices (isolated vertices allowed)
 *   W  ... G has no universal vertex
 *   P  ... G has a pendant vertex
 *   S  ... G has a pendant vertex whose neighbour is not universal
 *   S2 ... G has >= 2 pendant vertices whose neighbours are not universal          (Seed2)
 *   S3 ... G has a vertex u with |N[u]| = 3                                       (Seed3)
 *   G  ... both S2 and S3 ("good")
 * Bases: every even graph on <= NB vertices with its flags (file from mkbase.py; exhaustive for NB <= 9).
 * Rules (each is a lemma of PROOF.md section 2):
 *   padding: X_{n-1} -> X_n for every set; E_{n-1} -> W_n; P_{n-1} -> S_n;
 *   cone (L2):           W_{n-1} + 2^{n-1}          -> E_n
 *   E_2-join (L3):       E_{n-2} + 3*2^{n-2} - 2    -> E_n ;  W_{n-2} + 3*2^{n-2} - 2 -> W_n
 *   cone over G u K_1:   E_{n-2} + 2^{n-1}          -> P_n
 *   R = PH(pendant,1):   S_{n-2} + 3*2^{n-3}        -> S_n ;  S2_{n-2} + 3*2^{n-3} -> G_n     (L4, C5)
 *   PH(r3 vertex,2):     S3_{n-3} + 7*2^{n-4}       -> G_n                                      (L4, C5)
 *   seed A:  K_1 v (X u E_2) u K_1, |X| = n-4:  E_{n-4} + 2^{n-2}        -> S2_n
 *   seed B:  E_2 v (Z u K_1),       |Z| = n-3:  E_{n-3} + 3*2^{n-2} - 2  -> S3_n
 *   G_n -> S2_n, S3_n, S_n, P_n, W_n, E_n (a good graph has no universal vertex).
 * Mode "noR" drops the R move (pre-registered R3).
 * Checks printed:
 *   (1) E_n = formula F(n) = odd [1, M(n)] minus {3,7,13} minus X(n) for every n <= NM (missing / extra counts);
 *   (2) Gamma(c) = {n : 2^{n-1} + c in G_n} for c in {1,3,5,7,11,13};
 *   (3) for every odd j with L(j) <= JMAX: the least n0 with n0, n0+1 in Gamma(-j); must satisfy n0 <= L(j) + 8.
 * usage: evenclos NM basefile NB JMAX [noR] [dumpC]   (dumpC: print Gamma(c) for odd |c| <= dumpC instead of (3))
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define NX 40
enum {SE,SW,SP,SS,S2,S3,SG,NS};
static uint8_t *A[NS][NX];
static int NB, NM;
static inline void setb(uint8_t *a, uint64_t k) { a[k >> 4] |= 1 << ((k >> 1) & 7); }
static inline int getb(uint8_t *a, uint64_t k) { return (a[k >> 4] >> ((k >> 1) & 7)) & 1; }
static uint64_t sz(int n) { return (((uint64_t)1 << n) >> 4) + 1; }
static void orset(uint8_t *dst, uint8_t *src, int nsrc) { uint64_t m = sz(nsrc); for (uint64_t i = 0; i < m; i++) dst[i] |= src[i]; }
static void shiftadd(uint8_t *dst, int ndst, uint8_t *src, int nsrc, int64_t off) {
    uint64_t lim = (uint64_t)1 << nsrc, dl = (uint64_t)1 << ndst;
    for (uint64_t k = 1; k < lim; k += 2) if (getb(src, k)) { int64_t v = (int64_t)k + off; if (v > 0 && (uint64_t)v < dl) setb(dst, (uint64_t)v); }
}
static int L(uint64_t k) { int l = 0; while (k) { l++; k >>= 1; } return l; }
static uint64_t M(int n) { return (n & 1) ? ((uint64_t)1 << n) - n : ((uint64_t)1 << n) - n - 1; }
static int inF(int n, uint64_t k) {
    if (k < 1 || !(k & 1) || k > M(n)) return 0;
    if (k == 3 || k == 7 || k == 13) return 0;
    if ((n & 1) && n >= 5 && k == ((uint64_t)1 << n) - n - 2) return 0;
    return 1;
}
int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: evenclos NM basefile NB JMAX [noR] [dumpC]\n"); return 1; }
    NM = atoi(argv[1]); NB = atoi(argv[3]); int JMAX = atoi(argv[4]);
    int noR = (argc > 5 && !strcmp(argv[5], "noR")); int dumpC = (argc > 6) ? atoi(argv[6]) : 0;
    for (int s = 0; s < NS; s++) for (int n = 1; n <= NM; n++) A[s][n] = calloc(sz(n), 1);
    FILE *f = fopen(argv[2], "r"); int n, U, pany, pnu, iso, r3; unsigned long long d;
    while (fscanf(f, "%d %llu %d %d %d %d %d", &n, &d, &U, &pany, &pnu, &iso, &r3) == 7) {
        if (n > NB || n > NM) continue;
        setb(A[SE][n], d); if (!U) setb(A[SW][n], d); if (pany) setb(A[SP][n], d); if (pnu) setb(A[SS][n], d);
        if (pnu >= 2) setb(A[S2][n], d); if (r3) setb(A[S3][n], d); if (pnu >= 2 && r3) setb(A[SG][n], d);
    }
    fclose(f);
    for (n = 2; n <= NM; n++) {
        for (int s = 0; s < NS; s++) orset(A[s][n], A[s][n - 1], n - 1);
        orset(A[SW][n], A[SE][n - 1], n - 1); orset(A[SS][n], A[SP][n - 1], n - 1);
        if (n > NB) {
            if (!noR) shiftadd(A[SS][n], n, A[SS][n - 2], n - 2, (int64_t)3 << (n - 3));
            orset(A[SP][n], A[SS][n], n); shiftadd(A[SP][n], n, A[SE][n - 2], n - 2, (int64_t)1 << (n - 1));
            shiftadd(A[SW][n], n, A[SW][n - 2], n - 2, ((int64_t)3 << (n - 2)) - 2); orset(A[SW][n], A[SS][n], n);
            orset(A[SE][n], A[SW][n], n); shiftadd(A[SE][n], n, A[SW][n - 1], n - 1, (int64_t)1 << (n - 1));
            shiftadd(A[SE][n], n, A[SE][n - 2], n - 2, ((int64_t)3 << (n - 2)) - 2);
        }
        if (!noR) {
            if (n >= 3) shiftadd(A[SG][n], n, A[S2][n - 2], n - 2, (int64_t)3 << (n - 3));
            if (n >= 4) shiftadd(A[SG][n], n, A[S3][n - 3], n - 3, (int64_t)7 << (n - 4));
            if (n >= 6) shiftadd(A[S2][n], n, A[SE][n - 4], n - 4, (int64_t)1 << (n - 2));
            if (n >= 5) shiftadd(A[S3][n], n, A[SE][n - 3], n - 3, ((int64_t)3 << (n - 2)) - 2);
            orset(A[S2][n], A[SG][n], n); orset(A[S3][n], A[SG][n], n);
            /* good graphs feed the order-n sets too (they have pendants, no universal vertex) */
            orset(A[SS][n], A[SG][n], n); orset(A[SP][n], A[SG][n], n); orset(A[SW][n], A[SG][n], n); orset(A[SE][n], A[SG][n], n);
        }
    }
    long totmiss = 0, totextra = 0;
    for (n = 1; n <= NM; n++) {
        long miss = 0, extra = 0; uint64_t first = 0;
        for (uint64_t k = 1; k < ((uint64_t)1 << n); k += 2) {
            int want = inF(n, k), has = getb(A[SE][n], k);
            if (want && !has) { miss++; if (!first) first = k; }
            if (has && !want) extra++;
        }
        totmiss += miss; totextra += extra;
        printf("E(%d): missing %ld, extra %ld%s", n, miss, extra, miss ? "" : "\n");
        if (miss) printf(", first missing %llu\n", (unsigned long long)first);
    }
    printf("TOTAL E(n) vs formula for n <= %d: missing %ld, extra %ld\n", NM, totmiss, totextra);
    if (noR) return 0;
    if (dumpC) {
        for (int c = -dumpC; c <= dumpC; c += 2) { printf("G %d", c);
            for (n = 3; n <= NM; n++) { int64_t v = ((int64_t)1 << (n - 1)) + c; if (v > 0 && v < ((int64_t)1 << n) && getb(A[SG][n], (uint64_t)v)) printf(" %d", n); }
            printf("\n"); }
        return 0;
    }
    int cs[6] = {1, 3, 5, 7, 11, 13};
    for (int t = 0; t < 6; t++) { printf("Gamma(+%d):", cs[t]);
        for (n = 3; n <= NM; n++) { uint64_t v = ((uint64_t)1 << (n - 1)) + cs[t]; if (v < ((uint64_t)1 << n) && getb(A[SG][n], v)) printf(" %d", n); }
        printf("\n"); }
    long bad = 0; int worst = -99;
    for (uint64_t j = 1; j < ((uint64_t)1 << JMAX); j += 2) {
        int n0 = -1;
        for (n = 3; n + 1 <= NM; n++) {
            if (((uint64_t)1 << (n - 1)) <= j) continue;
            if (getb(A[SG][n], ((uint64_t)1 << (n - 1)) - j) && getb(A[SG][n + 1], ((uint64_t)1 << n) - j)) { n0 = n; break; }
        }
        if (n0 < 0 || n0 > L(j) + 8) { bad++; if (bad <= 20) printf("seed FAIL j=%llu L=%d n0=%d\n", (unsigned long long)j, L(j), n0); }
        else if (n0 - L(j) > worst) worst = n0 - L(j);
    }
    printf("seed base (odd j with L(j) <= %d): failures %ld; max over j of n0 - L(j) = %d (need <= 8)\n", JMAX, bad, worst);
    return 0;
}
