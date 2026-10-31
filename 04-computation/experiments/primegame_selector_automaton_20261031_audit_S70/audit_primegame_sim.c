/* Independent audit (S70 audit) of THM-4852 (1) / HYP-9448: PRIMEGAME phase law and the agreement of the two
   published fraction lists.  Written from scratch; shares no code with primegame_variant_phase_law_20261031.c.

   Design: every denominator is squarefree, so applicability of a fraction depends only on the SUPPORT of the
   configuration (which of the ten registers 2..29 are nonzero).  A 1024-entry table gives the first applicable
   fraction for every support mask; a step is one table lookup plus at most six register updates.

   Lists (labels checked against the sources, see AUDIT.md):
     C87 = Conway 1987, Theorem 1 (= Wikipedia / OEIS A007542, A202138/A203363; the owner's list):
           17/91 78/85 19/51 23/38 29/33 77/29 95/23 77/19 1/17 11/13 13/11 15/2 1/7 55/1
     G83 = Guy 1983 (OEIS A202138 comment; = the "fast version" of OEIS A007546):
           17/91 78/85 19/51 23/38 29/33 77/29 95/23 77/19 1/17 11/13 13/11 15/14 15/2 55/1

   Checks, for both lists, up to phase NMAX (k-th phase entry = k-th time the orbit ENTERS the state-free set,
   counting the start 2 as entry 1):
     (a) the k-th entry is exactly 2^k 7^(D(k)-1), D(k) = k/spf(k), D(1) = 1  (r2 = k is checked, not inferred);
     (b) the rule word of every state-free stretch is the predicted one:
           C87: (15/2)^a (1/7)^b (55/1)        G83: (15/14)^b (15/2)^(a-b) (55/1)        with (a,b) = (k, D(k)-1);
         and the configuration right after 55/1 is 3^a 5^(a+1) 11;
     (c) every stateful configuration has exactly one state prime, with exponent 1, and the fired rule is one of
         the first eleven (A..K), which are common to both lists;
     (d) the two lists have identical phase-entry and post-55/1 configurations (stored and compared);
     (e) the powers of 2 reached are exactly 2^p, p prime <= NMAX, in increasing order;
     (f) the number of steps from entry k-1 to entry k equals the closed form
           G83: f(k) = k-1 + (6k+2)(k-D(k)) + 2 sum_{d=D(k)}^{k-1} floor(k/d)    (A. P. Heinz's formula, OEIS A007546)
           C87: f(k) + D(k-1) - 1                                             (the extra 1/7 steps);
     (g) the cumulative step counts at 2^p agree with the 34 terms of OEIS A007546 (fast) and A007547 (slow) shown on
         the OEIS pages (fetched 2026-10-31).
   usage: audit_primegame_sim NMAX [control]
          ("control" falsifies the expected law at k = 97 only, to show that the checker fires)
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define NP 10
#define STATEMASK 0x3F0 /* bits 4..9 = registers 11,13,17,19,23,29 */
static const int PR[NP] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};

static const int C87_NUM[14] = {17, 78, 19, 23, 29, 77, 95, 77, 1, 11, 13, 15, 1, 55};
static const int C87_DEN[14] = {91, 85, 51, 38, 33, 29, 23, 19, 17, 13, 11, 2, 7, 1};
static const int G83_NUM[14] = {17, 78, 19, 23, 29, 77, 95, 77, 1, 11, 13, 15, 15, 55};
static const int G83_DEN[14] = {91, 85, 51, 38, 33, 29, 23, 19, 17, 13, 11, 14, 2, 1};

/* OEIS A007546 (fast) and A007547 (slow), the terms displayed on the OEIS pages */
static const long long A007546[34] = {19, 69, 280, 707, 2363, 3876, 8068, 11319, 19201, 36866, 45551, 75224, 101112,
    117831, 152025, 215384, 293375, 327020, 428553, 507519, 555694, 700063, 808331, 989526, 1273490, 1434366, 1530213,
    1710923, 1818254, 2019962, 2833089, 3104685, 3546320, 3720785};
static const long long A007547[34] = {19, 69, 281, 710, 2375, 3893, 8102, 11361, 19268, 36981, 45680, 75417, 101354,
    118093, 152344, 215797, 293897, 327571, 429229, 508284, 556494, 701008, 809381, 990746, 1274952, 1435957, 1531854,
    1712701, 1820085, 2021938, 2835628, 3107393, 3549288, 3723821};

typedef struct {
    const char *name;
    int dmask[14];
    int nd[14], di[14][NP], dv[14][NP];
    unsigned char first[1 << NP];
} Prog;

static void factor(int x, int e[NP]) {
    for (int i = 0; i < NP; i++) { e[i] = 0; while (x % PR[i] == 0) { x /= PR[i]; e[i]++; } }
    if (x != 1) { fprintf(stderr, "not 29-smooth\n"); exit(2); }
}

static void build(Prog *P, const char *name, const int *num, const int *den) {
    P->name = name;
    for (int k = 0; k < 14; k++) {
        int en[NP], ed[NP];
        factor(num[k], en); factor(den[k], ed);
        P->dmask[k] = 0; P->nd[k] = 0;
        for (int i = 0; i < NP; i++) {
            if (ed[i] > 1) { fprintf(stderr, "denominator not squarefree\n"); exit(2); }
            if (ed[i]) P->dmask[k] |= 1 << i;
            if (en[i] && ed[i]) { fprintf(stderr, "fraction not in lowest terms\n"); exit(2); }
            if (en[i] - ed[i]) { P->di[k][P->nd[k]] = i; P->dv[k][P->nd[k]] = en[i] - ed[i]; P->nd[k]++; }
        }
    }
    for (int m = 0; m < (1 << NP); m++) {
        P->first[m] = 255;
        for (int k = 0; k < 14; k++) if ((P->dmask[k] & ~m) == 0) { P->first[m] = (unsigned char)k; break; }
    }
}

static int spf(int n) { for (int p = 2; p * p <= n; p++) if (n % p == 0) return p; return n; }
static int CONTROL = 0; /* negative control: if set, the EXPECTED law is deliberately wrong at k = 97 */
static int Dof(int n) { return n == 1 ? 1 : n / spf(n); }
static int Dexp(int n) { return Dof(n) + (CONTROL && n == 97); }

typedef struct { int r2, r3, r5, r7; } C4;

static long long closed_form_g83(int k) { /* steps from entry k-1 to entry k, Guy's list, k >= 2 */
    int b = Dof(k);
    long long s = 0;
    for (int d = b; d <= k - 1; d++) s += k / d;
    return (long long)(k - 1) + (long long)(6 * (long long)k + 2) * (k - b) + 2 * s;
}

static int fail(const char *msg, const char *name, long long t, int k) {
    printf("FAIL [%s] %s at step %lld (phase %d)\n", name, msg, t, k);
    return 0;
}

/* returns 1 on success; fills entry[1..nmax], post[1..nmax-1], entry_step[1..nmax], pows/powsteps */
static int run(const Prog *P, int is_guy, int nmax, C4 *entry, C4 *post, long long *entry_step,
               int *pows, long long *powsteps, int *npows) {
    int x[NP] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int mask = 1, inphase = 0, k = 0;
    int a = 0, b = 0, stage = 0, cnt[3] = {0, 0, 0};
    *npows = 0;
    for (long long t = 0;; t++) {
        if ((mask & STATEMASK) == 0) {
            if (!inphase) {
                k++;
                int D = Dexp(k);
                if (x[0] != k || x[1] != 0 || x[2] != 0 || x[3] != D - 1)
                    return fail("phase law (k-th entry != 2^k 7^(D(k)-1))", P->name, t, k);
                entry[k].r2 = x[0]; entry[k].r3 = x[1]; entry[k].r5 = x[2]; entry[k].r7 = x[3];
                entry_step[k] = t;
                if (k == nmax) return 1;
                inphase = 1; a = k; b = D - 1; stage = 0; cnt[0] = cnt[1] = cnt[2] = 0;
            }
        } else {
            inphase = 0;
            int s = mask & STATEMASK;
            if (s & (s - 1)) return fail("two state primes present", P->name, t, k);
            int idx = __builtin_ctz(s);
            if (x[idx] != 1) return fail("state exponent != 1", P->name, t, k);
        }
        int r = P->first[mask];
        if (r == 255) return fail("halt", P->name, t, k);
        if (!inphase) {
            if (r > 10) return fail("a last-three rule fired in a stateful configuration", P->name, t, k);
        } else {
            if (r < 11) return fail("an A..K rule fired in a state-free configuration", P->name, t, k);
            int st = r - 11; /* 0,1,2 = the three last rules in list order */
            if (st < stage) return fail("state-free rule word out of order", P->name, t, k);
            stage = st; cnt[st]++;
        }
        for (int j = 0; j < P->nd[r]; j++) {
            int i = P->di[r][j];
            x[i] += P->dv[r][j];
            if (x[i] < 0) return fail("negative exponent", P->name, t, k);
            if (x[i]) mask |= 1 << i; else mask &= ~(1 << i);
        }
        if (inphase && r == 13) { /* 55/1 fired: end of the state-free stretch */
            int e0 = is_guy ? b : a, e1 = is_guy ? a - b : b;
            if (cnt[0] != e0 || cnt[1] != e1 || cnt[2] != 1)
                return fail("state-free rule word counts", P->name, t, k);
            if (x[0] != 0 || x[1] != a || x[2] != a + 1 || x[3] != 0 || x[4] != 1 || (mask & ~0x10 & STATEMASK))
                return fail("post-55/1 configuration != 3^a 5^(a+1) 11", P->name, t, k);
            post[k].r2 = x[0]; post[k].r3 = x[1]; post[k].r5 = x[2]; post[k].r7 = x[3];
        }
        if (mask == 1) { pows[*npows] = x[0]; powsteps[*npows] = t + 1; (*npows)++; }
    }
}

int main(int argc, char **argv) {
    int nmax = argc > 1 ? atoi(argv[1]) : 1000;
    if (argc > 2 && !strcmp(argv[2], "control")) { CONTROL = 1; printf("NEGATIVE CONTROL: expected law falsified at k = 97\n"); }
    static Prog C87, G83;
    build(&C87, "Conway1987/owner", C87_NUM, C87_DEN);
    build(&G83, "Guy1983", G83_NUM, G83_DEN);
    C4 *eC = calloc(nmax + 2, sizeof(C4)), *eG = calloc(nmax + 2, sizeof(C4));
    C4 *pC = calloc(nmax + 2, sizeof(C4)), *pG = calloc(nmax + 2, sizeof(C4));
    long long *sC = calloc(nmax + 2, sizeof(long long)), *sG = calloc(nmax + 2, sizeof(long long));
    int *powC = calloc(nmax + 2, sizeof(int)), *powG = calloc(nmax + 2, sizeof(int));
    long long *psC = calloc(nmax + 2, sizeof(long long)), *psG = calloc(nmax + 2, sizeof(long long));
    int nC, nG;
    clock_t t0 = clock();
    int okC = run(&C87, 0, nmax, eC, pC, sC, powC, psC, &nC);
    double tC = (double)(clock() - t0) / CLOCKS_PER_SEC;
    t0 = clock();
    int okG = run(&G83, 1, nmax, eG, pG, sG, powG, psG, &nG);
    double tG = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("NMAX = %d\n", nmax);
    printf("(a)(b)(c) Conway1987/owner list: %s  (%lld steps to entry %d, %.1f s)\n", okC ? "OK" : "FAILED", sC[nmax], nmax, tC);
    printf("(a)(b)(c) Guy1983 list:          %s  (%lld steps to entry %d, %.1f s)\n", okG ? "OK" : "FAILED", sG[nmax], nmax, tG);
    if (!okC || !okG) return 1;
    /* (d) */
    int same = 1;
    for (int k = 1; k <= nmax; k++)
        if (memcmp(&eC[k], &eG[k], sizeof(C4))) same = 0;
    for (int k = 1; k < nmax; k++)
        if (memcmp(&pC[k], &pG[k], sizeof(C4))) same = 0;
    printf("(d) identical phase-entry configurations (k = 1..%d) and post-55/1 configurations: %s\n", nmax, same ? "yes" : "NO");
    /* (e) */
    int np = 0, primes_ok = 1;
    for (int p = 2; p <= nmax; p++) if (spf(p) == p) { if (np >= nC || np >= nG || powC[np] != p || powG[np] != p) primes_ok = 0; np++; }
    if (nC != np || nG != np) primes_ok = 0;
    printf("(e) powers of 2 = 2^p for exactly the %d primes p <= %d, in order, both lists: %s\n", np, nmax, primes_ok ? "yes" : "NO");
    /* (f) */
    int cf_ok = 1; long long maxdiff = 0, extra = 0;
    for (int k = 2; k <= nmax; k++) {
        long long g = closed_form_g83(k);
        long long dg = sG[k] - sG[k - 1], dc = sC[k] - sC[k - 1];
        if (dg != g) cf_ok = 0;
        if (dc != g + Dof(k - 1) - 1) cf_ok = 0;
        if (dc - dg > maxdiff) maxdiff = dc - dg;
        extra += dc - dg;
    }
    printf("(f) per-phase step counts equal Heinz's closed form (Guy) and closed form + D(k-1)-1 (Conway), k = 2..%d: %s\n", nmax, cf_ok ? "yes" : "NO");
    printf("    total extra steps of the 1/7 list = sum_{k<%d} (D(k)-1) = %lld\n", nmax, extra);
    /* (g) */
    int m = np < 34 ? np : 34, oeis_ok = 1;
    for (int i = 0; i < m; i++) {
        if (psG[i] != A007546[i]) oeis_ok = 0;
        if (psC[i] != A007547[i]) oeis_ok = 0;
    }
    printf("(g) cumulative steps at 2^p for the first %d primes: Guy = OEIS A007546 (fast), Conway/owner = OEIS A007547 (slow): %s\n", m, oeis_ok ? "yes" : "NO");
    printf("    first values Guy: %lld %lld %lld %lld %lld   Conway/owner: %lld %lld %lld %lld %lld\n",
           psG[0], psG[1], psG[2], psG[3], psG[4], psC[0], psC[1], psC[2], psC[3], psC[4]);
    printf("    first phase entries (r2,r3,r5,r7):");
    for (int k = 1; k <= 12 && k <= nmax; k++) printf(" (%d,%d,%d,%d)", eC[k].r2, eC[k].r3, eC[k].r5, eC[k].r7);
    printf("\n");
    int all = okC && okG && same && primes_ok && cf_ok && oeis_ok;
    printf("ALL CHECKS %s\n", all ? "PASSED" : "FAILED");
    return all ? 0 : 1;
}
