/* audit_rand.c  --  independent audit program (not derived from the author's code).
 *
 * Conventions: vertices 0..n-1, out[a] = bit mask of out-neighbours of a, A[a][b] = 1 iff a -> b.
 * Cycles are counted DIRECTLY: c_3 over all triples, c_5 over all 5-subsets and all 24 cyclic arrangements.
 *
 * Modes
 *   p1 n N seed        N random labelled tournaments of order n: c_3, c_5 (direct) against the score formulas of Theorem P (1)
 *   p3 m N seed        N random labelled bases of order m (m <= 10): for all 2^m out-sets X the tournament R + v_X is built and
 *                      c_3, c_5 are counted directly; sum_X (-1)^(c_3) against 2^ceil(m/2) (Theorem T (2)) and
 *                      sum_X (-1)^(c_3+c_5) against the table of Theorem P (3), K_0, K_1 found by trying all u in F_2^m.
 *                      Every third base is forced to have a source, a sink, or a pair of twins-like structure to exercise
 *                      singular adjacency matrices; statistics of |K_0|, K_1 are printed.
 *   p3class m          the same for every class of order m read from stdin (gentourng format, '1' at (i,j), i<j, means i -> j)
 *                      and for its converse.
 *   tower n N seed     tr(A^3+A^5+A^7) mod 2 against C(|A^3 1|,2) + floor(n/4), and the two true analogues
 *   zero m N seed      fraction of random bases with K_1 non-empty and all of K_1 of weight 0 mod 4 (the criterion of
 *                      Theorem P (3) for G(R) = 0 when m = 1 mod 4); option: direct check on the first D bases (argv[5])
 *   kernel m           exhaustive: sum of |ker A| over all labelled tournaments of order m (m <= 7)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define MAXN 40
static int n;
static uint64_t out[MAXN + 1];
static int popc(uint64_t x) { return __builtin_popcountll(x); }
static int C2(int k) { return k * (k - 1) / 2; }

static uint64_t rs[4];
static uint64_t splitmix(uint64_t *x) { uint64_t z = (*x += 0x9E3779B97F4A7C15ull); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t rnd(void) { uint64_t r = rotl(rs[1] * 5, 7) * 9, t = rs[1] << 17; rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3]; rs[2] ^= t; rs[3] = rotl(rs[3], 45); return r; }
static void seed_rng(uint64_t seed) { uint64_t sx = seed; for (int i = 0; i < 4; i++) rs[i] = splitmix(&sx); }
static void random_tournament(int nn, uint64_t *o) {
    for (int i = 0; i < nn; i++) o[i] = 0;
    uint64_t bits = 0; int nb = 0;
    for (int i = 0; i < nn; i++) for (int j = i + 1; j < nn; j++) {
        if (!nb) { bits = rnd(); nb = 64; }
        int b = (int)(bits & 1); bits >>= 1; nb--;
        if (b) o[i] |= 1ull << j; else o[j] |= 1ull << i;
    }
}
#define ARC(o, a, b) (((o)[a] >> (b)) & 1)

/* direct counts */
static long count_c3(const uint64_t *o, int nn) {
    long c = 0;
    for (int a = 0; a < nn; a++) for (int b = a + 1; b < nn; b++) for (int d = b + 1; d < nn; d++)
        if ((ARC(o, a, b) && ARC(o, b, d) && ARC(o, d, a)) || (ARC(o, a, d) && ARC(o, d, b) && ARC(o, b, a))) c++;
    return c;
}
static const int P4[24][4] = {{0,1,2,3},{0,1,3,2},{0,2,1,3},{0,2,3,1},{0,3,1,2},{0,3,2,1},{1,0,2,3},{1,0,3,2},{1,2,0,3},{1,2,3,0},{1,3,0,2},{1,3,2,0},
    {2,0,1,3},{2,0,3,1},{2,1,0,3},{2,1,3,0},{2,3,0,1},{2,3,1,0},{3,0,1,2},{3,0,2,1},{3,1,0,2},{3,1,2,0},{3,2,0,1},{3,2,1,0}};
static long count_c5(const uint64_t *o, int nn) {
    long c = 0; int v[5];
    for (v[0] = 0; v[0] < nn; v[0]++) for (v[1] = v[0] + 1; v[1] < nn; v[1]++) for (v[2] = v[1] + 1; v[2] < nn; v[2]++)
    for (v[3] = v[2] + 1; v[3] < nn; v[3]++) for (v[4] = v[3] + 1; v[4] < nn; v[4]++)
        for (int k = 0; k < 24; k++) {
            int w1 = v[1 + P4[k][0]], w2 = v[1 + P4[k][1]], w3 = v[1 + P4[k][2]], w4 = v[1 + P4[k][3]];
            if (ARC(o, v[0], w1) && ARC(o, w1, w2) && ARC(o, w2, w3) && ARC(o, w3, w4) && ARC(o, w4, v[0])) c++;
        }
    return c;
}
/* number of 5-cycles through the vertex z (direct: all 4-subsets of the others, all 24 orders) */
static long count_c5_through(const uint64_t *o, int nn, int z) {
    long c = 0; int v[4];
    for (v[0] = 0; v[0] < nn; v[0]++) { if (v[0] == z) continue;
    for (v[1] = v[0] + 1; v[1] < nn; v[1]++) { if (v[1] == z) continue;
    for (v[2] = v[1] + 1; v[2] < nn; v[2]++) { if (v[2] == z) continue;
    for (v[3] = v[2] + 1; v[3] < nn; v[3]++) { if (v[3] == z) continue;
        for (int k = 0; k < 24; k++) {
            int w1 = v[P4[k][0]], w2 = v[P4[k][1]], w3 = v[P4[k][2]], w4 = v[P4[k][3]];
            if (ARC(o, z, w1) && ARC(o, w1, w2) && ARC(o, w2, w3) && ARC(o, w3, w4) && ARC(o, w4, z)) c++;
        } } } } }
    return c;
}
static long count_c3_through(const uint64_t *o, int nn, int z) {
    long c = 0;
    for (int a = 0; a < nn; a++) { if (a == z) continue; for (int b = 0; b < nn; b++) { if (b == z || b == a) continue;
        if (ARC(o, z, a) && ARC(o, a, b) && ARC(o, b, z)) c++; } }
    return c;
}

static long fail_T2, fail_P3, fail_wt, fail_through, nsing, nK1, nG0, nbases;
static long Ghist_val[256], Ghist_cnt[256]; static int nGh;
static void check_base(const uint64_t *ob, int m) {
    /* direct cube sums */
    uint64_t o[MAXN + 1]; long sT = 0, sP = 0;
    long c3R = count_c3(ob, m), c5R = count_c5(ob, m);
    for (uint64_t X = 0; X < (1ull << m); X++) {
        for (int a = 0; a < m; a++) o[a] = ob[a] | (((X >> a) & 1) ? 0 : (1ull << m));
        o[m] = X;
        long c3, c5;
        if (m <= 8) { c3 = count_c3(o, m + 1); c5 = count_c5(o, m + 1);
            if (c3 != c3R + count_c3_through(o, m + 1, m) || c5 != c5R + count_c5_through(o, m + 1, m)) fail_through++; }
        else { c3 = c3R + count_c3_through(o, m + 1, m); c5 = c5R + count_c5_through(o, m + 1, m); }
        sT += (c3 & 1) ? -1 : 1; sP += ((c3 + c5) & 1) ? -1 : 1;
    }
    if (sT != (1L << ((m + 1) / 2))) fail_T2++;
    long n00 = 0, n01 = 0, n10 = 0, n13 = 0, bad = 0;
    uint64_t full = (1ull << m) - 1;
    uint64_t inb[MAXN]; for (int c = 0; c < m; c++) { inb[c] = 0; for (int a = 0; a < m; a++) if (ARC(ob, a, c)) inb[c] |= 1ull << a; }
    for (uint64_t u = 0; u <= full; u++) {
        uint64_t im = 0; for (int c = 0; c < m; c++) if (popc(u & inb[c]) & 1) im |= 1ull << c;     /* A^T u */
        int w = popc(u) & 3;
        if (im == 0) { if (w == 0) n00++; else if (w == 1) n01++; else bad++; }
        if (im == full) { if (w == 0) n10++; else if (w == 3) n13++; else bad++; }
    }
    if (bad) fail_wt++;
    long K0 = n00 + n01, G;
    switch (m & 3) {
        case 0: G = (1L << (m / 2)) * K0; break;
        case 1: G = (1L << ((m - 1) / 2)) * (K0 + n13 - n10); break;
        case 2: G = (1L << (m / 2)) * (n00 + n13); break;
        default: G = (1L << ((m + 1) / 2)) * (n00 + n10); break;
    }
    if (G != sP) fail_P3++;
    if (K0 > 1) nsing++;
    if (n10 + n13 > 0) nK1++;
    if (sP == 0) nG0++;
    nbases++;
    int k; for (k = 0; k < nGh; k++) if (Ghist_val[k] == sP) { Ghist_cnt[k]++; break; }
    if (k == nGh && nGh < 256) { Ghist_val[nGh] = sP; Ghist_cnt[nGh] = 1; nGh++; }
}
static void report_base(const char *tag, int m) {
    printf("%s m=%d bases %ld : failures T2 %ld, P3 %ld, kernel weights %ld, through-vertex consistency %ld ; singular A: %ld ; K_1 non-empty: %ld ; G = 0: %ld\n",
           tag, m, nbases, fail_T2, fail_P3, fail_wt, fail_through, nsing, nK1, nG0);
    printf("   values of G (count):");
    for (int k = 0; k < nGh; k++) printf(" %ld(%ld)", Ghist_val[k], Ghist_cnt[k]);
    printf("\n");
}

/* F_2 matrix helpers (rows as bit masks) */
static void matmul(const uint64_t *A, const uint64_t *B, uint64_t *Cm, int nn) {      /* C = A B over F_2 */
    for (int i = 0; i < nn; i++) { uint64_t r = 0, a = A[i]; while (a) { int k = __builtin_ctzll(a); a &= a - 1; r ^= B[k]; } Cm[i] = r; }
}
static int trace2(const uint64_t *A, int nn) { int t = 0; for (int i = 0; i < nn; i++) t ^= (int)((A[i] >> i) & 1); return t; }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: see source\n"); return 1; }
    const char *mode = argv[1];
    if (!strcmp(mode, "p1")) {
        n = atoi(argv[2]); long N = atol(argv[3]); seed_rng(strtoull(argv[4], NULL, 10));
        long f1 = 0, f2 = 0, f3 = 0, odd3 = 0, odd5 = 0;
        for (long it = 0; it < N; it++) {
            random_tournament(n, out);
            long c3 = count_c3(out, n), c5 = count_c5(out, n);
            uint64_t s = 0, t = 0;
            for (int a = 0; a < n; a++) if (popc(out[a]) & 1) s |= 1ull << a;
            for (int a = 0; a < n; a++) if (popc(out[a] & s) & 1) t |= 1ull << a;
            int q1 = popc(s), q2 = popc(t);
            if ((c3 + C2(q1) + n / 4) & 1) f1++;
            if ((c3 + c5 + C2(q2) + n / 4) & 1) f2++;
            if ((c5 + C2(q1) + C2(q2)) & 1) f3++;
            /* classical formula c3 = C(n,3) - sum C(s_i,2) as a control of the direct count */
            long cl = (long)n * (n - 1) * (n - 2) / 6; for (int a = 0; a < n; a++) cl -= C2(popc(out[a]));
            if (cl != c3) f1 += 1000000;
            odd3 += c3 & 1; odd5 += c5 & 1;
        }
        printf("p1 n=%d N=%ld : failures c3 %ld, c3+c5 %ld, c5 %ld ; c3 odd in %ld, c5 odd in %ld samples\n", n, N, f1, f2, f3, odd3, odd5);
        return 0;
    }
    if (!strcmp(mode, "p3")) {
        int m = atoi(argv[2]); long N = atol(argv[3]); seed_rng(strtoull(argv[4], NULL, 10));
        uint64_t ob[MAXN];
        for (long it = 0; it < N; it++) {
            random_tournament(m, ob);
            int kind = (int)(it % 6);
            if (kind == 1) { /* vertex 0 becomes a source */ for (int a = 1; a < m; a++) { ob[a] &= ~1ull; ob[0] |= 1ull << a; } }
            if (kind == 3) { /* vertex 0 becomes a sink */ ob[0] = 0; for (int a = 1; a < m; a++) ob[a] |= 1ull; }
            if (kind == 5 && m >= 3) { /* vertices 0 and 1 get the same out-neighbours among 2..m-1 */
                for (int a = 2; a < m; a++) { if (ARC(ob, 0, a)) { ob[1] |= 1ull << a; ob[a] &= ~2ull; } else { ob[1] &= ~(1ull << a); ob[a] |= 2ull; } } }
            check_base(ob, m);
        }
        report_base("p3 random", m);
        return 0;
    }
    if (!strcmp(mode, "p3class")) {
        int m = atoi(argv[2]); int E = C2(m); char line[256]; uint64_t ob[MAXN], oc[MAXN];
        while (fgets(line, sizeof line, stdin)) {
            int len = 0; while (line[len] == '0' || line[len] == '1') len++;
            if (len != E) continue;
            int pos = 0; for (int i = 0; i < m; i++) ob[i] = oc[i] = 0;
            for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { int b = line[pos++] == '1';
                if (b) { ob[i] |= 1ull << j; oc[j] |= 1ull << i; } else { ob[j] |= 1ull << i; oc[i] |= 1ull << j; } }
            check_base(ob, m); check_base(oc, m);
        }
        report_base("p3class (each class and its converse)", m);
        return 0;
    }
    if (!strcmp(mode, "tower")) {
        n = atoi(argv[2]); long N = atol(argv[3]); seed_rng(strtoull(argv[4], NULL, 10));
        long f3 = 0, f5 = 0, f7 = 0;
        for (long it = 0; it < N; it++) {
            random_tournament(n, out);
            uint64_t A2[MAXN], A3[MAXN], A5[MAXN], A7[MAXN];
            matmul(out, out, A2, n); matmul(A2, out, A3, n); matmul(A3, A2, A5, n); matmul(A5, A2, A7, n);
            uint64_t full = (n == 64) ? ~0ull : ((1ull << n) - 1);
            int w1 = 0, w2 = 0, w3 = 0;             /* |A 1|, |A^2 1|, |A^3 1| over F_2 */
            for (int a = 0; a < n; a++) { w1 += popc(out[a] & full) & 1; w2 += popc(A2[a] & full) & 1; w3 += popc(A3[a] & full) & 1; }
            int t3 = trace2(A3, n), t5 = trace2(A5, n), t7 = trace2(A7, n);
            if (t3 != ((C2(w1) + n / 4) & 1)) f3++;
            if ((t3 ^ t5) != ((C2(w2) + n / 4) & 1)) f5++;
            if ((t3 ^ t5 ^ t7) != ((C2(w3) + n / 4) & 1)) f7++;
        }
        printf("tower n=%d N=%ld : tr A^3 = C(|A1|,2)+floor(n/4) fails %ld ; tr(A^3+A^5) = C(|A^2 1|,2)+floor(n/4) fails %ld ; tr(A^3+A^5+A^7) = C(|A^3 1|,2)+floor(n/4) fails %ld (%.4f)\n",
               n, N, f3, f5, f7, (double)f7 / N);
        return 0;
    }
    if (!strcmp(mode, "zero")) {
        int m = atoi(argv[2]); long N = atol(argv[3]); seed_rng(strtoull(argv[4], NULL, 10));
        long D = argc > 5 ? atol(argv[5]) : 0;
        long nz = 0, nk1 = 0, fdirect = 0, ndirect = 0; uint64_t ob[MAXN];
        uint64_t full = (1ull << m) - 1;
        for (long it = 0; it < N; it++) {
            random_tournament(m, ob);
            /* solve A^T u = 1 and find ker A^T by elimination on the rows r_c = (in-mask of c | rhs bit) */
            uint64_t rows[MAXN]; int rhs[MAXN];
            for (int c = 0; c < m; c++) { rows[c] = 0; for (int a = 0; a < m; a++) if (ARC(ob, a, c)) rows[c] |= 1ull << a; rhs[c] = 1; }
            int piv[MAXN], npiv = 0, where[MAXN]; for (int a = 0; a < m; a++) where[a] = -1;
            int r = 0;
            for (int col = 0; col < m && r < m; col++) {
                int p = -1; for (int i = r; i < m; i++) if ((rows[i] >> col) & 1) { p = i; break; }
                if (p < 0) continue;
                uint64_t tr = rows[p]; rows[p] = rows[r]; rows[r] = tr; int tb = rhs[p]; rhs[p] = rhs[r]; rhs[r] = tb;
                for (int i = 0; i < m; i++) if (i != r && ((rows[i] >> col) & 1)) { rows[i] ^= rows[r]; rhs[i] ^= rhs[r]; }
                where[col] = r; piv[npiv++] = col; r++;
            }
            int solvable = 1; for (int i = r; i < m; i++) if (rhs[i]) solvable = 0;
            int isz = 0;
            if (solvable) {
                nk1++;
                /* particular solution and kernel basis */
                uint64_t u0 = 0; for (int col = 0; col < m; col++) if (where[col] >= 0 && rhs[where[col]]) u0 |= 1ull << col;
                uint64_t basis[MAXN]; int nb = 0;
                for (int col = 0; col < m; col++) if (where[col] < 0) { uint64_t v = 1ull << col;
                    for (int c2 = 0; c2 < m; c2++) if (where[c2] >= 0 && ((rows[where[c2]] >> col) & 1)) v |= 1ull << c2; basis[nb++] = v; }
                isz = 1;
                for (uint64_t k = 0; k < (1ull << nb) && isz; k++) { uint64_t u = u0; for (int j = 0; j < nb; j++) if ((k >> j) & 1) u ^= basis[j];
                    if (popc(u) & 3) isz = 0; }
                (void)piv;
            }
            if (isz) nz++;
            if (it < D) {
                /* direct: cube sum of (-1)^(increment of c3+c5), increment = number of walks of length 1 and 3 from X to V\X mod 2 */
                uint64_t A2[MAXN], A3[MAXN]; matmul(ob, ob, A2, m); matmul(A2, ob, A3, m);
                long s = 0;
                for (uint64_t X = 0; X <= full; X++) { int par = 0; uint64_t rem = X, nX = full & ~X;
                    while (rem) { int a = __builtin_ctzll(rem); rem &= rem - 1; par ^= (popc(ob[a] & nX) ^ popc(A3[a] & nX)) & 1; }
                    s += par ? -1 : 1; }
                ndirect++; if ((s == 0) != (isz == 1)) fdirect++;
            }
        }
        printf("zero m=%d N=%ld : K_1 non-empty %ld (%.5f) ; criterion for G = 0 holds in %ld (%.5f +- %.5f) ; direct cube sums on %ld bases, disagreements with the criterion %ld\n",
               m, N, nk1, (double)nk1 / N, nz, (double)nz / N, sqrt((double)nz / N * (1 - (double)nz / N) / N), ndirect, fdirect);
        return 0;
    }
    if (!strcmp(mode, "kernel")) {
        int m = atoi(argv[2]); int E = C2(m); unsigned long long tot = 0;
        for (uint64_t t = 0; t < (1ull << E); t++) {
            uint64_t o[8]; for (int i = 0; i < m; i++) o[i] = 0; int pos = 0;
            for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if ((t >> pos) & 1) o[i] |= 1ull << j; else o[j] |= 1ull << i; pos++; }
            for (uint64_t x = 0; x < (1ull << m); x++) { int ok = 1; for (int a = 0; a < m && ok; a++) if (popc(o[a] & x) & 1) ok = 0; tot += ok; }
        }
        printf("kernel m=%d : sum over all 2^%d labelled tournaments of |ker A| = %llu\n", m, E, tot);
        return 0;
    }
    return 1;
}
