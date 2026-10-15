/* Audit E (2026-10-14), THM-4655 statement 1, independent brute force.
   For every word w in {0,1}^m (ALL words, strong or not), m = 1..MMAX, build the ACTUAL tournament
   (u_1..u_m transitive u_i -> u_j for i<j; z -> u_i iff w_i = 1, else u_i -> z), and:
     (a) decide strong connectivity by transitive closure (Warshall), compare with [w_1 = 1 and w_m = 0];
     (b) enumerate ALL Hamiltonian paths by depth-first search on the adjacency matrix (no use of the A-bijection),
         map each path to A = {u's before z}, check injectivity, and check that the image is exactly the set of A
         satisfying (A empty or w_maxA = 0) and (B empty or w_minB = 1);
     (c) compare H with the closed form |K| + sum_{beta<alpha, w_beta=1, w_alpha=0} 2^(alpha-beta-1)
         (general K = {k : (k=0 or w_k=0) and (k=m or w_(k+1)=1)}; for strong words also K = {0,m} u {01-boundaries});
     (d) for strong words, m >= 3: exact P(u_2 before z) from the path tallies; check 1/3 < P <= 2/3, equality 2/3 iff
         w = 1 0 1^(m-3) 0; report min P and its argmin;
     (e) exact delta(T) over all pairs from the path tallies; family minimum per m with argmin; also delta at the
         equality word, and whether some pair OTHER than (z,u_2) is balanced there.
   Usage: ./tpo_dfs_bruteforce MMAX   (MMAX <= 13)  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n, m;
static int adjm[16][16];          /* adjm[x][y] = 1 iff x -> y */
static int path[16];
static u64 bef[16][16];           /* bef[x][y] = #paths with x before y */
static u64 H;
static unsigned char *seenA;      /* size 2^m */
static long long injfail = 0, orderfail = 0;
static int pos[16];

static void record(void) {
    /* path complete: path[0..n-1] */
    H++;
    int zpos = -1;
    for (int i = 0; i < n; i++) { pos[path[i]] = i; if (path[i] == m) zpos = i; }
    unsigned A = 0;
    for (int i = 0; i < zpos; i++) A |= 1u << path[i];
    /* the u-runs must be increasing (consequence of the tournament, checked anyway) */
    for (int i = 0; i + 1 < zpos; i++) if (path[i] > path[i + 1]) orderfail++;
    for (int i = zpos + 1; i + 1 < n; i++) if (path[i] > path[i + 1]) orderfail++;
    if (seenA[A]) injfail++;
    seenA[A] = 1;
    for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) if (x != y && pos[x] < pos[y]) bef[x][y]++;
}

static void dfs(int depth, unsigned used) {
    if (depth == n) { record(); return; }
    int last = path[depth - 1];
    for (int v = 0; v < n; v++) if (!(used >> v & 1) && adjm[last][v]) {
        path[depth] = v; dfs(depth + 1, used | 1u << v);
    }
}

static int is_strong(void) {
    int R[16][16];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) R[i][j] = (i == j) || adjm[i][j];
    for (int k = 0; k < n; k++) for (int i = 0; i < n; i++) if (R[i][k]) for (int j = 0; j < n; j++) if (R[k][j]) R[i][j] = 1;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (!R[i][j]) return 0;
    return 1;
}

/* compare fractions a/b < c/d */
static int lt(u64 a, u64 b, u64 c, u64 d) { return (unsigned __int128)a * d < (unsigned __int128)c * b; }

int main(int argc, char **argv) {
    int MMAX = atoi(argv[1]);
    long long tot_strongfail = 0, tot_bijfail = 0, tot_Hfail = 0, tot_Kfail = 0, tot_Pfail = 0, tot_eqfail = 0, tot_words = 0, tot_strong = 0;
    for (m = 1; m <= MMAX; m++) {
        n = m + 1;
        seenA = malloc((size_t)1 << m);
        long long words = 0, strongc = 0, strongfail = 0, bijfail = 0, Hfail = 0, Kfail = 0, Pfail = 0, eqfail = 0;
        u64 pmin_num = 1, pmin_den = 1, pmax_num = 0, pmax_den = 1; unsigned pmin_w = 0;
        u64 dmin_num = 1, dmin_den = 1; unsigned dmin_w = 0; int dmin_mult = 0;
        u64 eq_delta_num = 0, eq_delta_den = 1; int eq_other_balanced = -1; u64 eq_other_num = 0, eq_other_den = 1;
        u64 Hmax = 0;
        long long vw_ok = 0, vw_tight = 0; u64 vwmin_num = 1, vwmin_den = 1; unsigned vwmin_w = 0; int eq_z_tight = -1;
        for (unsigned w = 0; w < (1u << m); w++) {    /* bit i-1 = w_i */
            words++;
            memset(adjm, 0, sizeof adjm);
            for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) adjm[i][j] = 1;
            for (int i = 0; i < m; i++) { if (w >> i & 1) adjm[m][i] = 1; else adjm[i][m] = 1; }
            int strong = is_strong();
            int crit = (w & 1) && !(w >> (m - 1) & 1);
            if (strong != crit) strongfail++;
            /* enumerate paths */
            memset(seenA, 0, (size_t)1 << m); memset(bef, 0, sizeof bef); H = 0; injfail = 0; orderfail = 0;
            for (int s = 0; s < n; s++) { path[0] = s; dfs(1, 1u << s); }
            if (H > Hmax) Hmax = H;
            /* image = valid A ? */
            int bad = (injfail || orderfail);
            unsigned full = (1u << m) - 1;
            for (unsigned A = 0; A <= full; A++) {
                unsigned B = full ^ A;
                int valid = 1;
                if (A) { int mx = 31 - __builtin_clz(A); if (w >> mx & 1) valid = 0; }       /* need w_maxA = 0 */
                if (B) { int mn = __builtin_ctz(B); if (!(w >> mn & 1)) valid = 0; }          /* need w_minB = 1 */
                if (valid != seenA[A]) bad = 1;
            }
            if (bad) bijfail++;
            /* closed form for H (1-based indices) */
            u64 K = 0;
            for (int k = 0; k <= m; k++) {
                int ok = (k == 0 || !(w >> (k - 1) & 1)) && (k == m || (w >> k & 1));
                K += ok;
            }
            u64 cross = 0;
            for (int be = 1; be <= m; be++) if (w >> (be - 1) & 1)
                for (int al = be + 1; al <= m; al++) if (!(w >> (al - 1) & 1)) cross += 1ULL << (al - be - 1);
            if (K + cross != H) Hfail++;
            if (strong) {
                strongc++;
                u64 K2 = 2; for (int k = 1; k < m; k++) if (!(w >> (k - 1) & 1) && (w >> k & 1)) K2++;
                if (K2 != K) Kfail++;
            }
            if (strong && m >= 3) {
                u64 in2 = bef[1][m];                 /* u_2 is vertex 1, z is vertex m */
                if (!(3 * in2 > H && 3 * in2 <= 2 * H)) Pfail++;
                unsigned eqw = 1u; for (int i = 3; i <= m - 1; i++) eqw |= 1u << (i - 1);   /* 1 0 1^(m-3) 0 */
                int is_eq = (3 * in2 == 2 * H);
                int pred = (m >= 4 && w == eqw);
                if (is_eq != pred) eqfail++;
                if (lt(in2, H, pmin_num, pmin_den)) { pmin_num = in2; pmin_den = H; pmin_w = w; }
                if (lt(pmax_num, pmax_den, in2, H)) { pmax_num = in2; pmax_den = H; }
                /* delta over all pairs */
                u64 best = 0; u64 bestother = 0;
                for (int x = 0; x < n; x++) for (int y = x + 1; y < n; y++) {
                    u64 c = bef[x][y]; u64 mn = c < H - c ? c : H - c;
                    if (mn > best) best = mn;
                    if (!((x == 1 && y == m) || (x == m && y == 1)) && mn > bestother) bestother = mn;
                }
                if (lt(best, H, dmin_num, dmin_den)) { dmin_num = best; dmin_den = H; dmin_w = w; dmin_mult = 1; }
                else if ((unsigned __int128)best * dmin_den == (unsigned __int128)dmin_num * H) dmin_mult++;
                /* HYP-9279(B) vertex-wise balance */
                {
                    u64 vmn = 1, vmd = 1; int allv = 1;
                    for (int v = 0; v < n; v++) {
                        u64 bv = 0;
                        for (int x = 0; x < n; x++) if (x != v) { u64 c = bef[v][x], mn = c < H - c ? c : H - c; if (mn > bv) bv = mn; }
                        if (3 * bv < H) allv = 0;
                        if (lt(bv, H, vmn, vmd)) { vmn = bv; vmd = H; }
                        if (pred && v == m) eq_z_tight = (3 * bv == H);
                    }
                    if (allv) vw_ok++;
                    if (3 * vmn == vmd) vw_tight++;
                    if (lt(vmn, vmd, vwmin_num, vwmin_den)) { vwmin_num = vmn; vwmin_den = vmd; vwmin_w = w; }
                }
                if (pred) { eq_delta_num = best; eq_delta_den = H; eq_other_num = bestother; eq_other_den = H;
                            eq_other_balanced = (3 * bestother >= H); }
            }
        }
        printf("m=%2d n=%2d: words %lld (strong %lld), Hmax %llu | strong-criterion mismatches %lld | bijection failures %lld | "
               "H-formula failures %lld | K-form failures %lld\n", m, n, words, strongc, Hmax, strongfail, bijfail, Hfail, Kfail);
        if (m >= 3) {
            char wb[40]; for (int i = 0; i < m; i++) wb[i] = '0' + (pmin_w >> i & 1); wb[m] = 0;
            char wd[40]; for (int i = 0; i < m; i++) wd[i] = '0' + (dmin_w >> i & 1); wd[m] = 0;
            printf("        P(u2<z): violations of (1/3,2/3] %lld, equality-iff failures %lld, min %llu/%llu=%.6f at w=%s, max %llu/%llu=%.6f\n",
                   Pfail, eqfail, pmin_num, pmin_den, (double)pmin_num / pmin_den, wb, pmax_num, pmax_den, (double)pmax_num / pmax_den);
            printf("        family min delta %llu/%llu = %.6f at w=%s (#words attaining %d)", dmin_num, dmin_den, (double)dmin_num / dmin_den, wd, dmin_mult);
            if (m >= 4) printf(" | equality word: delta %llu/%llu=%.6f, best pair other than (z,u2) %llu/%llu=%.6f balanced=%d",
                               eq_delta_num, eq_delta_den, (double)eq_delta_num / eq_delta_den, eq_other_num, eq_other_den,
                               (double)eq_other_num / eq_other_den, eq_other_balanced);
            printf("\n");
            char wv[40]; for (int i = 0; i < m; i++) wv[i] = '0' + (vwmin_w >> i & 1); wv[m] = 0;
            printf("        HYP-9279(B): strong words with every vertex balanced %lld/%lld; min vertex balance %llu/%llu=%.6f at w=%s; "
                   "words where some vertex balance is exactly 1/3: %lld; at the equality word z's best partner is exactly 1/3: %d\n",
                   vw_ok, strongc, vwmin_num, vwmin_den, (double)vwmin_num / vwmin_den, wv, vw_tight, eq_z_tight);
        }
        fflush(stdout);
        tot_strongfail += strongfail; tot_bijfail += bijfail; tot_Hfail += Hfail; tot_Kfail += Kfail; tot_Pfail += Pfail; tot_eqfail += eqfail;
        tot_words += words; tot_strong += strongc;
        free(seenA);
    }
    printf("TOTAL: %lld words (%lld strong); mismatches: strong %lld, bijection %lld, H %lld, K %lld, P-range %lld, equality %lld\n",
           tot_words, tot_strong, tot_strongfail, tot_bijfail, tot_Hfail, tot_Kfail, tot_Pfail, tot_eqfail);
    return 0;
}
