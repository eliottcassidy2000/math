/* audit_small.c  --  independent audit program (not derived from the author's code).
 *
 * Exhaustive enumeration of ALL labelled tournaments of order n = 1..NMAX (NMAX <= 7), everything by brute force:
 *   - directed cycles of every length by depth-first enumeration (each cycle once, from its least vertex);
 *   - Hamiltonian paths by depth-first enumeration;
 *   - self-converse test and non-trivial-automorphism test by trying all n! permutations.
 * Conventions: vertices 0..n-1; the pair (i,j), i<j, has bit position C(j,2)+i in the tournament index;
 *   bit = 1 means i -> j.  A[a][b] = 1 iff a -> b.  The base R of T is T minus the last vertex v = n-1 (index = low bits);
 *   X = out-neighbourhood of v in T.
 * Checks, per tournament T of order n:
 *   (H)   H = 1 + 2 c_odd (mod 4)
 *   (P1)  c_3 = C(q1,2) + floor(n/4), c_3 + c_5 = C(q2,2) + floor(n/4) (mod 2), s = A1, t = As, q1 = |s|, q2 = |t|
 *   (T1)  c_3(T) - c_3(R) = e_R(X, V\X) = sum_{a in X} s_R(a) - C(|X|,2)            (exact integers)
 *   (P2)  (c_3+c_5)(T) - (c_3+c_5)(R) = x^T (A + A^3)(1+x) = C(|y|,2) + u_sigma . y   (mod 2)
 * Per base R of order m = n-1 (after the loop):
 *   (T2)  sum_X (-1)^(c_3(R+v_X)) = 2^ceil(m/2)
 *   (P3)  sum_X (-1)^((c_3+c_5)(R+v_X)) = table of Theorem P (3), kernel sets K_0, K_1 found by trying all u in F_2^m
 * Per order: joint law of (s,t) (P4), sums of signs for every set K of odd lengths, D_n, H mod 8 histogram,
 *   labelled self-converse counts and the three strata.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NMAXMAX 7
static int n;
static int A[NMAXMAX][NMAXMAX];
static unsigned outm[NMAXMAX];
static long cyc[NMAXMAX + 2];
static long Hc;

static void dfs_cycle(int start, int last, int len, unsigned vis) {
    if (len >= 3 && A[last][start]) cyc[len]++;
    for (int w = start + 1; w < n; w++)
        if (!((vis >> w) & 1) && A[last][w]) dfs_cycle(start, w, len + 1, vis | (1u << w));
}
static void dfs_ham(int last, unsigned vis, int cnt) {
    if (cnt == n) { Hc++; return; }
    for (int w = 0; w < n; w++)
        if (!((vis >> w) & 1) && A[last][w]) dfs_ham(w, vis | (1u << w), cnt + 1);
}
static int popc(unsigned x) { return __builtin_popcount(x); }
static int C2(int k) { return k * (k - 1) / 2; }

/* permutations */
static int nperm; static unsigned char *perms;
static void gen_perms(int nn) {
    int f = 1; for (int i = 2; i <= nn; i++) f *= i;
    perms = realloc(perms, (size_t)f * nn + 8); nperm = 0;
    int p[NMAXMAX]; for (int i = 0; i < nn; i++) p[i] = i;
    for (;;) {
        for (int i = 0; i < nn; i++) perms[(size_t)nperm * nn + i] = (unsigned char)p[i];
        nperm++;
        int i = nn - 2; while (i >= 0 && p[i] > p[i + 1]) i--;
        if (i < 0) break;
        int j = nn - 1; while (p[j] < p[i]) j--;
        int t = p[i]; p[i] = p[j]; p[j] = t;
        for (int a = i + 1, b = nn - 1; a < b; a++, b--) { t = p[a]; p[a] = p[b]; p[b] = t; }
    }
}
static int is_self_converse(void) {
    for (int k = 0; k < nperm; k++) {
        const unsigned char *p = perms + (size_t)k * n; int ok = 1;
        for (int j = 1; j < n && ok; j++) for (int i = 0; i < j; i++)
            if (A[i][j] != A[p[j]][p[i]]) { ok = 0; break; }
        if (ok) return 1;
    }
    return 0;
}
static int has_nontrivial_aut(void) {
    for (int k = 1; k < nperm; k++) {           /* k = 0 is the identity (lexicographically first) */
        const unsigned char *p = perms + (size_t)k * n; int ok = 1;
        for (int j = 1; j < n && ok; j++) for (int i = 0; i < j; i++)
            if (A[i][j] != A[p[i]][p[j]]) { ok = 0; break; }
        if (ok) return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 6;
    int do_sym = argc > 2 ? atoi(argv[2]) : 1;   /* 1: brute-force self-converse / automorphism tests */
    if (NMAX > NMAXMAX) NMAX = NMAXMAX;
    /* parities of c3, c5 of every tournament of the previous order (the bases) */
    unsigned char *prev3 = NULL, *prev5 = NULL; long *prevc3 = NULL;
    long fails_H = 0, fails_P1a = 0, fails_P1b = 0, fails_T1 = 0, fails_P2a = 0, fails_P2b = 0, fails_T2 = 0, fails_P3 = 0,
         fails_P4 = 0, fails_wt = 0;
    for (n = 1; n <= NMAX; n++) {
        int E = C2(n), m = n - 1, Eb = C2(m);
        long NT = 1L << E, NB = 1L << Eb;
        unsigned char *cur3 = malloc(NT), *cur5 = malloc(NT); long *curc3 = malloc(NT * sizeof(long));
        int *sumT = calloc(NB, sizeof(int)), *sumP = calloc(NB, sizeof(int));
        long *joint = calloc((size_t)1 << (2 * n), sizeof(long));
        long sgn[8]; memset(sgn, 0, sizeof sgn);        /* index: bit0 = 3, bit1 = 5, bit2 = 7 */
        long Dsum = 0, h8[8]; memset(h8, 0, sizeof h8);
        long SC = 0, SCsigned = 0, symN = 0, symS = 0, rigN = 0, rigS = 0, sc_h8[8]; memset(sc_h8, 0, sizeof sc_h8);
        long maxG = 0, minG = 1L << 40; long Gvals[64]; int nG = 0;
        if (do_sym) gen_perms(n);
        for (long t = 0; t < NT; t++) {
            for (int j = 0; j < n; j++) { outm[j] = 0; }
            for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
                int b = (int)((t >> (C2(j) + i)) & 1);
                A[i][j] = b; A[j][i] = 1 - b;
                if (b) outm[i] |= 1u << j; else outm[j] |= 1u << i;
            }
            for (int i = 0; i < n; i++) A[i][i] = 0;
            memset(cyc, 0, sizeof cyc);
            for (int v = 0; v < n; v++) dfs_cycle(v, v, 1, 1u << v);
            Hc = 0; for (int v = 0; v < n; v++) dfs_ham(v, 1u << v, 1);
            long c3 = cyc[3], c5 = cyc[5], c7 = cyc[7];
            long codd = c3 + c5 + c7;
            if (((Hc - 1 - 2 * codd) % 4 + 4) % 4 != 0 || (Hc & 1) == 0) fails_H++;
            /* P1 */
            unsigned s = 0, tt = 0;
            for (int a = 0; a < n; a++) if (popc(outm[a]) & 1) s |= 1u << a;
            for (int a = 0; a < n; a++) if (popc(outm[a] & s) & 1) tt |= 1u << a;
            int q1 = popc(s), q2 = popc(tt);
            if (((c3 - C2(q1) - n / 4) & 1) != 0) fails_P1a++;
            if (((c3 + c5 - C2(q2) - n / 4) & 1) != 0) fails_P1b++;
            joint[((long)s << n) | tt]++;
            cur3[t] = (unsigned char)(c3 & 1); cur5[t] = (unsigned char)(c5 & 1); curc3[t] = c3;
            int b3 = (int)(c3 & 1), b5 = (int)(c5 & 1), b7 = (int)(c7 & 1);
            for (int K = 0; K < 8; K++) {
                int e = ((K & 1) ? b3 : 0) ^ ((K & 2) ? b5 : 0) ^ ((K & 4) ? b7 : 0);
                sgn[K] += e ? -1 : 1;
            }
            int chi = (Hc % 4 == 1) ? 1 : -1;
            Dsum += chi; h8[Hc % 8]++;
            if (do_sym) {
                int sc = is_self_converse();
                if (sc) { SC++; SCsigned += chi; sc_h8[Hc % 8]++; }
                else if (has_nontrivial_aut()) { symN++; symS += chi; }
                else { rigN++; rigS += chi; }
            }
            if (n >= 2) {
                /* base R = T - v, v = n-1; X = out-set of v */
                long b = t & (NB - 1);
                unsigned X = outm[m] & ((1u << m) - 1);
                unsigned full = (1u << m) - 1;
                unsigned ob[NMAXMAX]; int sR[NMAXMAX];
                for (int a = 0; a < m; a++) { ob[a] = outm[a] & full; sR[a] = popc(ob[a]); }
                /* T1 exact */
                long e = 0, ss = 0;
                for (int a = 0; a < m; a++) if ((X >> a) & 1) { e += popc(ob[a] & ~X & full); ss += sR[a]; }
                if (c3 - prevc3[b] != e || e != ss - C2(popc(X))) fails_T1++;
                /* P2 mod 2 */
                int sigma = popc(X) & 1;
                unsigned sv = 0, Asv = 0, y = 0;
                for (int a = 0; a < m; a++) if (sR[a] & 1) sv |= 1u << a;
                for (int a = 0; a < m; a++) { if (popc(ob[a] & sv) & 1) Asv |= 1u << a; if ((popc(ob[a] & X) + sigma) & 1) y |= 1u << a; }
                unsigned u = sigma ? Asv : (Asv ^ sv ^ full);
                int rhs = (C2(popc(y)) + popc(u & y)) & 1;
                /* x^T (A + A^3)(1+x): walks of length 1 and 3 from X to V\X, by matrix-vector products over the integers */
                long w1 = 0, w3 = 0; long vec[NMAXMAX], v2[NMAXMAX], v3[NMAXMAX];
                for (int a = 0; a < m; a++) vec[a] = ((X >> a) & 1) ? 0 : 1;        /* indicator of V\X */
                for (int a = 0; a < m; a++) { v2[a] = 0; for (int c = 0; c < m; c++) if ((ob[a] >> c) & 1) v2[a] += vec[c]; }  /* A(1+x) */
                for (int a = 0; a < m; a++) if ((X >> a) & 1) w1 += v2[a];
                for (int a = 0; a < m; a++) { v3[a] = 0; for (int c = 0; c < m; c++) if ((ob[a] >> c) & 1) v3[a] += v2[c]; }  /* A^2(1+x) */
                for (int a = 0; a < m; a++) { v2[a] = 0; for (int c = 0; c < m; c++) if ((ob[a] >> c) & 1) v2[a] += v3[c]; }  /* A^3(1+x) */
                for (int a = 0; a < m; a++) if ((X >> a) & 1) w3 += v2[a];
                int lhs = (int)((c3 + c5) & 1) ^ prev3[b] ^ prev5[b];
                if (lhs != (int)((w1 + w3) & 1)) fails_P2a++;
                if (lhs != rhs) fails_P2b++;
                sumT[b] += (c3 & 1) ? -1 : 1;
                sumP[b] += ((c3 + c5) & 1) ? -1 : 1;
            }
        }
        /* per-base checks */
        if (n >= 2) {
            for (long b = 0; b < NB; b++) {
                unsigned ob[NMAXMAX]; memset(ob, 0, sizeof ob);
                for (int j = 1; j < m; j++) for (int i = 0; i < j; i++) {
                    int bit = (int)((b >> (C2(j) + i)) & 1);
                    if (bit) ob[i] |= 1u << j; else ob[j] |= 1u << i;
                }
                if (sumT[b] != (1 << ((m + 1) / 2))) fails_T2++;
                /* kernel sets by trying all u: (A^T u)_c = parity of #{a in u : a -> c} */
                long n00 = 0, n01 = 0, n10 = 0, n13 = 0, bad = 0;
                for (unsigned u = 0; u < (1u << m); u++) {
                    unsigned img = 0;
                    for (int c = 0; c < m; c++) { int par = 0; for (int a = 0; a < m; a++) if (((u >> a) & 1) && ((ob[a] >> c) & 1)) par ^= 1; if (par) img |= 1u << c; }
                    int w = popc(u) & 3;
                    if (img == 0) { if (w == 0) n00++; else if (w == 1) n01++; else bad++; }
                    if (m > 0 && img == (1u << m) - 1) { if (w == 0) n10++; else if (w == 3) n13++; else bad++; }
                }
                if (bad) fails_wt++;
                long K0 = n00 + n01, G;
                switch (m & 3) {
                    case 0: G = (1L << (m / 2)) * K0; break;
                    case 1: G = (1L << ((m - 1) / 2)) * (K0 + n13 - n10); break;
                    case 2: G = (1L << (m / 2)) * (n00 + n13); break;
                    default: G = (1L << ((m + 1) / 2)) * (n00 + n10); break;
                }
                if (G != sumP[b]) fails_P3++;
                if (sumP[b] > maxG) maxG = sumP[b];
                if (sumP[b] < minG) minG = sumP[b];
                int found = 0; for (int k = 0; k < nG; k++) if (Gvals[k] == sumP[b]) found = 1;
                if (!found && nG < 64) Gvals[nG++] = sumP[b];
            }
        }
        /* P4: joint law */
        {
            long expect_s = (n >= 1) ? (1L << (E - (n - 1))) : 1;
            for (unsigned s = 0; s < (1u << n); s++) {
                long tot = 0; for (unsigned tt = 0; tt < (1u << n); tt++) tot += joint[((long)s << n) | tt];
                int adm_s = ((popc(s) - E) & 1) == 0;
                if (tot != (adm_s ? expect_s : 0)) fails_P4++;
                if (!adm_s) continue;
                int q = popc(s);
                for (unsigned tt = 0; tt < (1u << n); tt++) {
                    long c = joint[((long)s << n) | tt], want;
                    if (q == 0) want = (tt == 0) ? expect_s : 0;
                    else if (q == n) want = (tt == (1u << n) - 1) ? expect_s : 0;
                    else {
                        int in = popc(tt & s), out = popc(tt & ~s);
                        int ok = (((in - C2(q)) & 1) == 0) && (((out - C2(n - q)) & 1) == 0);
                        want = ok ? (1L << (E - (2 * n - 3))) : 0;
                    }
                    if (c != want) fails_P4++;
                }
            }
        }
        printf("n=%d  N=2^%d\n", n, E);
        printf("  signs (numerators over 2^%d): E[{}]=%ld E[{3}]=%ld E[{5}]=%ld E[{3,5}]=%ld E[{7}]=%ld E[{3,7}]=%ld E[{5,7}]=%ld E[{3,5,7}]=%ld\n",
               E, sgn[0], sgn[1], sgn[2], sgn[3], sgn[4], sgn[5], sgn[6], sgn[7]);
        printf("  D_n numerator (sum of chi) = %ld ; H mod 8 counts: 1:%ld 3:%ld 5:%ld 7:%ld\n", Dsum, h8[1], h8[3], h8[5], h8[7]);
        if (do_sym)
            printf("  labelled self-converse SC_n = %ld (signed %ld; H mod 8: 1:%ld 3:%ld 5:%ld 7:%ld); sym non-sc: %ld (signed %ld); rigid non-sc: %ld (signed %ld)\n",
                   SC, SCsigned, sc_h8[1], sc_h8[3], sc_h8[5], sc_h8[7], symN, symS, rigN, rigS);
        if (n >= 2) {
            printf("  bases of order m=%d: %ld ; cube sums of c3+c5: min %ld max %ld ; distinct values:", m, NB, minG, maxG);
            for (int k = 0; k < nG; k++) printf(" %ld", Gvals[k]);
            printf("\n");
        }
        fflush(stdout);
        free(prev3); free(prev5); free(prevc3); prev3 = cur3; prev5 = cur5; prevc3 = curc3;
        free(sumT); free(sumP); free(joint);
    }
    printf("FAILURES: H=1+2codd mod 4: %ld ; P1 c3: %ld ; P1 c3+c5: %ld ; T1: %ld ; P2 (walk form): %ld ; P2 (weight form): %ld ; T2: %ld ; P3: %ld ; kernel weights: %ld ; P4: %ld\n",
           fails_H, fails_P1a, fails_P1b, fails_T1, fails_P2a, fails_P2b, fails_T2, fails_P3, fails_wt, fails_P4);
    return 0;
}
