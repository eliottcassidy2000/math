/* lrc14_k11_band_mu.c -- HYP-5775: the k=11 hfloor band, C port (mac-mini S65 cont.7)
 *
 * Certifies mu(E) >= BAR + MARGIN for every primitive 11-shape E = {0=e1<...<e11=d},
 * reflection-reduced, for a given diameter d and first-middle-element range [e2lo, e2hi].
 *
 * mu(E) = meas{x in [0,1) : maxCircGap{frac(e_i x)} > 1/7}.
 * Exact structure: the verdict is constant between candidate breakpoints
 *   x = m/delta  (adjacency)  and  x = (7m+-1)/(7 delta)  (gap = 1/7 crossings),
 * over all pairwise differences delta. Symmetry mu-region x <-> 1-x: scan [0, 1/2] only,
 * mu = 2 * acc. Early exit once acc >= (BAR+MARGIN)/2. Shapes that fail to certify are
 * PRINTED for exact re-evaluation in Python (expected: none; near-degenerate midpoints
 * within 1e-11 of 1/7 are also flagged for safety).
 *
 * Modes:  certify d e2lo e2hi        -> "CERT d e2lo e2hi <#shapes> <#flagged>" + FLAG lines
 *         evalone e1 e2 ... e11      -> full mu (no early exit), for cross-validation
 *
 * Validated against lrc14_k11_band_mu_macmini_S65cont6.py (float=exact=MC) on d=19..24.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double BAR = 0.3312;
static const double MARGIN = 0.02;
#define K 11
#define MAXD 36

static int cmpd(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* returns: 1 = certified (mu >= BAR+MARGIN, or exact full mu if fullmu!=NULL),
   0 = not certified / degenerate  */
static int eval_shape(const int *E, double *fullmu) {
    static _Thread_local double cand[4200];
    int ncand = 0;
    unsigned char seen[MAXD] = {0};
    for (int i = 0; i < K; i++)
        for (int j = i + 1; j < K; j++) {
            int dl = E[j] - E[i];
            if (dl > 0 && dl < MAXD) seen[dl] = 1;
        }
    for (int dl = 1; dl < MAXD; dl++) {
        if (!seen[dl]) continue;
        for (int m = 0; ; m++) {
            double x0 = (double)m / dl;
            double x1 = (7.0 * m + 1.0) / (7.0 * dl);
            double x2 = (7.0 * m + 6.0) / (7.0 * dl);
            int any = 0;
            if (x0 <= 0.5) { cand[ncand++] = x0; any = 1; }
            if (x1 <= 0.5) { cand[ncand++] = x1; any = 1; }
            if (x2 <= 0.5) { cand[ncand++] = x2; any = 1; }
            if (!any) break;
        }
    }
    cand[ncand++] = 0.5;
    qsort(cand, ncand, sizeof(double), cmpd);
    double need = fullmu ? 2.0 : (BAR + MARGIN) / 2.0;
    double acc = 0.0, prev = 0.0;
    for (int c = 0; c < ncand; c++) {
        double cur = cand[c];
        if (cur <= prev) continue;
        double mid = 0.5 * (prev + cur);
        /* phases */
        double ph[K];
        for (int i = 0; i < K; i++) {
            double v = E[i] * mid;
            ph[i] = v - (long)v;                     /* E[i]*mid < 18: (long) floor ok, v>=0 */
        }
        /* insertion sort */
        for (int i = 1; i < K; i++) {
            double key = ph[i]; int j = i - 1;
            while (j >= 0 && ph[j] > key) { ph[j + 1] = ph[j]; j--; }
            ph[j + 1] = key;
        }
        double mg = 1.0 - ph[K - 1] + ph[0];
        for (int i = 1; i < K; i++) {
            double g = ph[i] - ph[i - 1];
            if (g > mg) mg = g;
        }
        double diff = mg - 1.0 / 7.0;
        if (diff > -1e-11 && diff < 1e-11) return 0;   /* degenerate midpoint: flag */
        if (diff > 0) {
            acc += (cur - prev);
            if (!fullmu && acc >= need) return 1;      /* certified */
        }
        prev = cur;
    }
    if (fullmu) { *fullmu = 2.0 * acc; return 1; }
    return 0;                                          /* completed scan below threshold */
}

static int gcd_i(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

int main(int argc, char **argv) {
    if (argc >= 13 && !strcmp(argv[1], "evalone")) {
        int E[K];
        for (int i = 0; i < K; i++) E[i] = atoi(argv[2 + i]);
        double mu;
        eval_shape(E, &mu);
        printf("MU %.15f\n", mu);
        return 0;
    }
    if (argc != 5 || strcmp(argv[1], "certify")) {
        fprintf(stderr, "usage: %s certify d e2lo e2hi | evalone e1..e11\n", argv[0]);
        return 1;
    }
    int d = atoi(argv[2]), lo = atoi(argv[3]), hi = atoi(argv[4]);
    long nshapes = 0, nflag = 0;
    int mid9[9];                                       /* the 9 middle elements */
    for (int e2 = lo; e2 <= hi && e2 <= d - 9; e2++) {
        /* odometer over e3..e10 in (e2, d) choosing 8 more */
        mid9[0] = e2;
        int idx[8];
        for (int i = 0; i < 8; i++) idx[i] = e2 + 1 + i;
        while (1) {
            int E[K];
            E[0] = 0;
            E[1] = e2;
            for (int i = 0; i < 8; i++) E[2 + i] = idx[i];
            E[10] = d;
            /* primitivity: gcd of all nonzero elements */
            int g = E[1];
            for (int i = 2; i < K; i++) g = gcd_i(g, E[i]);
            if (g == 1) {
                /* reflection reduction: R = sorted{d - e} ; skip if R < E lexicographically */
                int R[K];
                for (int i = 0; i < K; i++) R[i] = d - E[K - 1 - i];
                int cmp = 0;
                for (int i = 0; i < K; i++) {
                    if (R[i] != E[i]) { cmp = (R[i] < E[i]) ? -1 : 1; break; }
                }
                if (cmp >= 0) {
                    nshapes++;
                    if (!eval_shape(E, NULL)) {
                        nflag++;
                        printf("FLAG %d :", d);
                        for (int i = 0; i < K; i++) printf(" %d", E[i]);
                        printf("\n");
                    }
                }
            }
            /* next combination idx[0..7] within (e2, d-1], ending before d */
            int p = 7;
            while (p >= 0 && idx[p] == d - 1 - (7 - p)) p--;
            if (p < 0) break;
            idx[p]++;
            for (int q = p + 1; q < 8; q++) idx[q] = idx[q - 1] + 1;
        }
    }
    printf("CERT %d %d %d %ld %ld\n", d, lo, hi, nshapes, nflag);
    return 0;
}
