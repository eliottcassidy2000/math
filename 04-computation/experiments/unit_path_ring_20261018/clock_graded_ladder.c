/* C port of merge_18_19_11_22_20261018/clock_graded_ladder.py (mac-mini-2026-10-18-unitpath): THM-4667 synchronous cross
   terms E_d - E_(d-1) resolved by literal clock lag in the ideal band model (clock_duality R7/3c).
   F(e, r, b): signed clock profile of all backward words of length e from a node of class r mod 3^(e+1) in band b.
   Children: a_j = 2 if B^j(t) = 1 mod 3, 1 if = 2 mod 3, none if 0 (B(t) = 4t + 1); clock a_j + 2j, sign (-1)^a_j,
   class (2^a_j B^j(t) - 1)/3 mod 3^e, band b + 1 - a_j - 2j. F(0, r, b) = h(r mod 3, 3 * 2^-b).
   Node weight per class mod 3^(d+1) and band b: (3/8)/3^(d+1) * M_d/2^(b+1), M_d = 4 (3/2)^d.
   Usage: ./clock_graded_ladder DMIN DMAX BMAX   (prints the same fields as the Python version, plus lags 4..8) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef struct { int n; int *c; long long *v; } Prof;     /* sorted clocks, integer values */
static int BMAX;
static long long P3[40];
static Prof **memo[40];          /* memo[e][r * BMAX + b] */
static const int kap[3] = {0, -3, 3}, lam[3] = {-1, -1, 2};
static long long hval(int r, int b) {   /* h(r, y) with y = 3 * 2^-b */
    long long tot = 0; int k = 0;
    /* 4^k * y <= 2  <=>  3 * 2^(2k - b) <= 2 */
    while (1) { double y = 3.0 * ldexp(1.0, 2 * k - b); if (y <= 2.0) { tot += kap[(r + k) % 3]; k++; } else break; }
    { double y = 3.0 * ldexp(1.0, 2 * k - b); if (y <= 4.0) tot += lam[(r + k) % 3]; }
    return tot;
}
/* children of (e, r, b): arrays of sign, clock, child class, child band; returns count */
static int children(int e, long long r, int b, int *sg, int *cl, long long *rc, int *bc) {
    long long mod = P3[e + 1], rho = r % mod; int j = 0, n = 0;
    while (b - 2 * j >= 0) {
        int a = (rho % 3 == 1) ? 2 : (rho % 3 == 2) ? 1 : 0;
        if (a) {
            int band = b + 1 - a - 2 * j;
            if (band >= 0) {
                sg[n] = (a == 2) ? 1 : -1; cl[n] = a + 2 * j;
                long long t = ((a == 2 ? 4 : 2) * rho - 1) % mod; if (t < 0) t += mod;
                rc[n] = t / 3; bc[n] = band; n++;
            }
        }
        rho = (4 * rho + 1) % mod; j++;
    }
    return n;
}
static int cmpint(const void *a, const void *b) { return (*(const int *)a) - (*(const int *)b); }
static Prof *F(int e, long long r, int b) {
    if (b < 0) return NULL;
    Prof **slot = &memo[e][r * BMAX + b];
    if (*slot) return *slot;
    Prof *P = calloc(1, sizeof(Prof));
    if (e == 0) {
        long long v = hval((int)(r % 3), b);
        if (v) { P->n = 1; P->c = malloc(sizeof(int)); P->v = malloc(sizeof(long long)); P->c[0] = 0; P->v[0] = v; }
        *slot = P; return P;
    }
    int sg[64], cl[64], bc[64]; long long rc[64];
    int nk = children(e, r, b, sg, cl, rc, bc);
    /* accumulate into a dense buffer over clocks */
    int maxc = 0; Prof *kid[64];
    for (int i = 0; i < nk; i++) { kid[i] = F(e - 1, rc[i], bc[i]); if (kid[i] && kid[i]->n) { int m = kid[i]->c[kid[i]->n - 1] + cl[i]; if (m > maxc) maxc = m; } }
    long long *acc = calloc(maxc + 1, sizeof(long long));
    for (int i = 0; i < nk; i++) if (kid[i]) for (int t = 0; t < kid[i]->n; t++) acc[kid[i]->c[t] + cl[i]] += sg[i] * kid[i]->v[t];
    int cnt = 0; for (int c = 0; c <= maxc; c++) if (acc[c]) cnt++;
    P->n = cnt; P->c = malloc(sizeof(int) * (cnt ? cnt : 1)); P->v = malloc(sizeof(long long) * (cnt ? cnt : 1));
    cnt = 0; for (int c = 0; c <= maxc; c++) if (acc[c]) { P->c[cnt] = c; P->v[cnt] = acc[c]; cnt++; }
    free(acc); *slot = P; return P;
}
int main(int argc, char **argv) {
    int dmin = atoi(argv[1]), dmax = atoi(argv[2]); BMAX = atoi(argv[3]);
    P3[0] = 1; for (int i = 1; i < 40; i++) P3[i] = P3[i - 1] * 3;
    for (int e = 0; e < dmax; e++) memo[e] = calloc((size_t)P3[e + 1] * BMAX, sizeof(Prof *));
    for (int d = dmin; d <= dmax; d++) {
        double lag[4096]; memset(lag, 0, sizeof(lag)); double Ed = 0, M = 4.0 * pow(1.5, d);
        for (int b = 0; b < BMAX; b++) {
            double w = (3.0 / 8.0) / (double)P3[d + 1] * M / ldexp(1.0, b + 1);
            for (long long r = 0; r < P3[d + 1]; r++) {
                int sg[64], cl[64], bc[64]; long long rc[64];
                int nk = children(d, r, b, sg, cl, rc, bc);
                Prof *kid[64]; double tot = 0;
                for (int i = 0; i < nk; i++) { kid[i] = F(d - 1, rc[i], bc[i]); if (kid[i]) for (int t = 0; t < kid[i]->n; t++) tot += sg[i] * (double)kid[i]->v[t]; }
                Ed += w * tot * tot;
                for (int i = 0; i < nk; i++) { if (!kid[i] || !kid[i]->n) continue;
                    for (int k = i + 1; k < nk; k++) { if (!kid[k] || !kid[k]->n) continue;
                        for (int s = 0; s < kid[i]->n; s++) { int c1 = kid[i]->c[s] + cl[i]; double v1 = sg[i] * (double)kid[i]->v[s];
                            for (int t = 0; t < kid[k]->n; t++) { int c2 = kid[k]->c[t] + cl[k]; int L = c1 > c2 ? c1 - c2 : c2 - c1;
                                lag[L] += w * v1 * sg[k] * (double)kid[k]->v[t]; } } } }
            }
        }
        double cross = 0; for (int L = 0; L < 4096; L++) cross += 2 * lag[L];
        printf("d = %d: E_d = %.12f ; cross term E_d - E_(d-1) = %+.12f ; lag 0: %+.8f lag 1: %+.8f lag 2: %+.8f lag 3: %+.8f lag 4: %+.8f lag 5: %+.8f lag 6: %+.8f (x2 = THM-4667 units)\n",
               d, Ed, cross, 2 * lag[0], 2 * lag[1], 2 * lag[2], 2 * lag[3], 2 * lag[4], 2 * lag[5], 2 * lag[6]);
        fflush(stdout);
    }
    return 0;
}
