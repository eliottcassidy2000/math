// Independent audit of THM-4808 (f) / HYP-9418 (2): branch-character height sums by depth-first traversal of the
// Berggren tree (own code; chi is propagated along the tree, no continued fractions, no lane code).
//   mode 1:  height_sum 1 XMAX            -> S(X), N(X) and E(X) = S - 1 - N_U - N_D at 1-2-5 checkpoints <= XMAX
//   mode 2:  height_sum 2 YMAX WINDOW [X0] -> full resolution on every integer X in [1, YMAX]:
//            sup |E(X)|/sqrt X over X >= 100, >= 1000, >= 10^4 (with argmax), number of sign changes of E,
//            the log-mean (1/log Y) int_{X0}^{Y} E(X) X^{-3/2} dX, and the same for S(X)/sqrt X.
// Nodes: (m, n) with c = m^2 + n^2; children U = (2m-n, m), A = (2m+n, m), D = (m+2n, n); chi(A-child) = -chi.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

typedef struct { int64_t m, n; int8_t s; } node;

static node *stk;
static int64_t stkcap;

static void push(int64_t *top, int64_t m, int64_t n, int8_t s) {
    if (*top >= stkcap) { stkcap *= 2; stk = realloc(stk, stkcap * sizeof(node)); if (!stk) { fprintf(stderr, "oom\n"); exit(1); } }
    stk[*top].m = m; stk[*top].n = n; stk[*top].s = s; (*top)++;
}

// N_U(X) = #{k >= 1 : (k+2)^2 + (k+1)^2 <= X},  N_D(X) = #{k >= 1 : (2k+2)^2 + 1 <= X}
static int64_t NU(int64_t X) { int64_t k = (int64_t)floor((sqrt(2.0 * (double)X - 1.0) - 3.0) / 2.0); if (k < 0) k = 0;
    while (k > 0 && (k + 2) * (k + 2) + (k + 1) * (k + 1) > X) k--; while ((k + 3) * (k + 3) + (k + 2) * (k + 2) <= X) k++; return k; }
static int64_t ND(int64_t X) { int64_t k = (int64_t)floor((sqrt((double)X - 1.0) - 2.0) / 2.0); if (k < 0) k = 0;
    while (k > 0 && (2 * k + 2) * (2 * k + 2) + 1 > X) k--; while ((2 * k + 4) * (2 * k + 4) + 1 <= X) k++; return k; }

int main(int argc, char **argv) {
    int mode = atoi(argv[1]);
    stkcap = 1 << 20; stk = malloc(stkcap * sizeof(node));
    if (mode == 1) {
        int64_t XMAX = (int64_t)atof(argv[2]);
        int64_t cps[64]; int ncp = 0;
        for (int64_t p = 10; p <= XMAX; p *= 10) { int64_t a[3] = {p, 2 * p, 5 * p}; for (int i = 0; i < 3; i++) if (a[i] <= XMAX) cps[ncp++] = a[i]; }
        // sort checkpoints
        for (int i = 0; i < ncp; i++) for (int j = i + 1; j < ncp; j++) if (cps[j] < cps[i]) { int64_t t = cps[i]; cps[i] = cps[j]; cps[j] = t; }
        int64_t Sb[64] = {0}, Nb[64] = {0};
        int64_t top = 0; push(&top, 2, 1, 1);
        while (top > 0) {
            node v = stk[--top];
            int64_t c = v.m * v.m + v.n * v.n;
            if (c > XMAX) continue;
            int lo = 0, hi = ncp - 1;                 // first checkpoint >= c
            while (lo < hi) { int mid = (lo + hi) / 2; if (cps[mid] >= c) hi = mid; else lo = mid + 1; }
            Sb[lo] += v.s; Nb[lo] += 1;
            push(&top, 2 * v.m - v.n, v.m, v.s);
            push(&top, 2 * v.m + v.n, v.m, (int8_t)-v.s);
            push(&top, v.m + 2 * v.n, v.n, v.s);
        }
        printf("X N(X) S(X) S/sqrtX E(X) E/sqrtX N_U N_D\n");
        int64_t S = 0, N = 0;
        for (int i = 0; i < ncp; i++) {
            S += Sb[i]; N += Nb[i];
            int64_t E = S - 1 - NU(cps[i]) - ND(cps[i]);
            printf("%lld %lld %lld %.5f %lld %.5f %lld %lld\n", (long long)cps[i], (long long)N, (long long)S,
                   S / sqrt((double)cps[i]), (long long)E, E / sqrt((double)cps[i]), (long long)NU(cps[i]), (long long)ND(cps[i]));
        }
        return 0;
    }
    // mode 2: windows of the per-c array
    int64_t YMAX = (int64_t)atof(argv[2]), W = (int64_t)atof(argv[3]);
    int8_t *arr = malloc(W);
    if (!arr) { fprintf(stderr, "oom\n"); return 1; }
    int64_t S = 0, nu = 0, nd = 0;            // running S(X), N_U(X), N_D(X)
    int64_t knu = 1, knd = 1;                 // next k to be added to N_U, N_D
    double best[3] = {0, 0, 0}; int64_t arg[3] = {0, 0, 0}; int64_t lim[3] = {100, 1000, 10000};
    int64_t signchanges = 0; int lastsign = 0;
    double logint = 0.0, logintS = 0.0;       // int_{X0}^{Y} E X^{-3/2} dX, int S X^{-3/2} dX
    int64_t X0 = argc > 4 ? (int64_t)atof(argv[4]) : 5;   // start of the log integral (default: root height)
    // sampled heights of the lane's table
    int64_t samp[] = {100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000,
                      5000000, 10000000, 20000000, 50000000, 100000000, 200000000, 500000000, 1000000000};
    int nsamp = sizeof(samp) / sizeof(samp[0]); int si = 0;
    // sup over the windows [10^k, 10^(k+1)) too
    double decsup[20] = {0}; int64_t decarg[20] = {0};
    for (int64_t L = 1; L <= YMAX; L += W) {
        int64_t R = L + W - 1; if (R > YMAX) R = YMAX;
        memset(arr, 0, W);
        int64_t top = 0; push(&top, 2, 1, 1);
        while (top > 0) {
            node v = stk[--top];
            int64_t c = v.m * v.m + v.n * v.n;
            if (c > R) continue;
            if (c >= L) arr[c - L] += v.s;
            push(&top, 2 * v.m - v.n, v.m, v.s);
            push(&top, 2 * v.m + v.n, v.m, (int8_t)-v.s);
            push(&top, v.m + 2 * v.n, v.n, v.s);
        }
        for (int64_t X = L; X <= R; X++) {
            S += arr[X - L];
            while ((knu + 2) * (knu + 2) + (knu + 1) * (knu + 1) <= X) { nu++; knu++; }
            while ((2 * knd + 2) * (2 * knd + 2) + 1 <= X) { nd++; knd++; }
            int64_t E = (X >= 5) ? S - 1 - nu - nd : 0;
            double r = fabs((double)E) / sqrt((double)X);
            for (int i = 0; i < 3; i++) if (X >= lim[i] && r > best[i]) { best[i] = r; arg[i] = X; }
            int dec = (int)floor(log10((double)X) + 1e-12);
            if (r > decsup[dec]) { decsup[dec] = r; decarg[dec] = X; }
            int sg = (E > 0) - (E < 0);
            if (X >= 100 && sg != 0) { if (lastsign != 0 && sg != lastsign) signchanges++; lastsign = sg; }
            if (X >= X0) {   // E and S are constant on [X, X+1)
                double w = 2.0 * (1.0 / sqrt((double)X) - 1.0 / sqrt((double)X + 1.0));
                logint += (double)E * w; logintS += (double)S * w;
            }
            if (si < nsamp && X == samp[si]) {
                printf("sample X=%lld S=%lld E=%lld E/sqrtX=%.5f\n", (long long)X, (long long)S, (long long)E, E / sqrt((double)X));
                si++;
            }
        }
        fprintf(stderr, "window [%lld, %lld] done\n", (long long)L, (long long)R);
    }
    for (int i = 0; i < 3; i++) printf("sup_{%lld <= X <= %lld} |E(X)|/sqrt X = %.5f at X = %lld\n", (long long)lim[i], (long long)YMAX, best[i], (long long)arg[i]);
    for (int d = 2; d < 20 && pow(10, d) <= YMAX; d++) printf("decade [1e%d, 1e%d): sup |E|/sqrtX = %.5f at X = %lld\n", d, d + 1, decsup[d], (long long)decarg[d]);
    printf("sign changes of E(X) over integers 100 <= X <= %lld (zeros skipped): %lld\n", (long long)YMAX, (long long)signchanges);
    double lY = log((double)YMAX + 1.0) - log((double)X0);
    printf("log-mean of E/sqrtX on [%lld, %lld]: %.5f ; log-mean of S/sqrtX: %.5f  (1/2 + 1/sqrt2 = %.5f)\n",
           (long long)X0, (long long)YMAX, logint / lY, logintS / lY, 0.5 + 1.0 / sqrt(2.0));
    return 0;
}
