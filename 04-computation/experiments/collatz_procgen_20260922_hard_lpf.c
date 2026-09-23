/*
 * collatz_procgen_20260922_hard_lpf.c  (hard-class lane, collatz-procgen-20260922)
 *
 * Repetition profile of a finite binary word via the Longest Previous Factor array.
 *
 *   LPF[j] = max_{0 <= a < j} LCE(a, j)      (LCE = longest common extension, overlaps allowed)
 *
 * For the periodic approximant with U = x[0,a), V = x[a,j) the common prefix of x and U V^inf
 * has length lambda = j + LCE(a,j), and |UV| = j.  Hence (Adamczewski-Bugeaud Diophantine exponent)
 *
 *   Dio(x) = 1 + limsup_j LPF[j]/j,
 *
 * and Bugeaud-Kim's function r(n) (length of the shortest prefix containing two occurrences of
 * some length-n factor) is r(n) = n + min{ j : LPF[j] >= n }.
 *
 * Algorithm: suffix array by prefix doubling with counting sort (O(N log N)), Kasai LCP,
 * previous/next-smaller-position scans with chained range minima (O(N)).
 *
 * Usage: hard_lpf WORDFILE [mu] [w0 w1] [topK] [jmin]
 *   WORDFILE : ASCII '0'/'1' (other bytes ignored)
 *   mu       : height exponent for the linear gain estimate G(j) = LPF[j] - (mu-1) j  (default 1.585)
 *   w0 w1    : letter weights (e.g. left Perron-Frobenius vector of a substitution; default 1 1)
 *   topK     : number of top linear-gain candidates printed (default 12)
 *   jmin     : minimal |UV| for candidates (default 64)
 * Output (stdout, one record per line):
 *   N <n>
 *   SCALE s jlo jhi maxratio j a lpf  wmaxratio wj wa wlpf  maxgain gj ga glpf   (uncapped j only)
 *   CAND j a lpf gain                                                          (top linear gains)
 *   RN n r(n)            for n = 1..min(N/4, 2^16) at geometric-ish sampling
 *   REP n1 n2 min_{n1<=n<=n2} r(n)/n
 *   DIOTAIL frac  max_{j >= frac*N/2, uncapped} (j+LPF)/j   (tail estimate)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static unsigned char *S;
static int N;

static void die(const char *m) { fprintf(stderr, "%s\n", m); exit(1); }

static void build_sa(int *sa, int *rk, int *tmp) {
    int n = N, i, j, k, m, r;
    int *cnt = (int *)calloc((size_t)(n > 3 ? n : 3) + 2, sizeof(int));
    if (!cnt) die("oom cnt");
    /* initial ranks 1,2 ; 0 reserved for positions beyond the end */
    for (i = 0; i < n; i++) rk[i] = S[i] + 1;
    m = 2;
    memset(cnt, 0, (size_t)(m + 1) * sizeof(int));
    for (i = 0; i < n; i++) cnt[rk[i]]++;
    for (i = 1; i <= m; i++) cnt[i] += cnt[i - 1];
    for (i = n - 1; i >= 0; i--) sa[--cnt[rk[i]]] = i;
    for (k = 1;; k <<= 1) {
        int p = 0;
        for (i = n - k; i < n; i++) if (i >= 0) tmp[p++] = i;
        for (j = 0; j < n; j++) if (sa[j] >= k) tmp[p++] = sa[j] - k;
        memset(cnt, 0, (size_t)(m + 1) * sizeof(int));
        for (i = 0; i < n; i++) cnt[rk[i]]++;
        for (i = 1; i <= m; i++) cnt[i] += cnt[i - 1];
        for (j = n - 1; j >= 0; j--) sa[--cnt[rk[tmp[j]]]] = tmp[j];
        tmp[sa[0]] = 1; r = 1;
        for (j = 1; j < n; j++) {
            int a = sa[j - 1], b = sa[j];
            int ra2 = (a + k < n) ? rk[a + k] : 0, rb2 = (b + k < n) ? rk[b + k] : 0;
            if (rk[a] != rk[b] || ra2 != rb2) r++;
            tmp[b] = r;
        }
        memcpy(rk, tmp, (size_t)n * sizeof(int));
        m = r;
        if (r == n) break;
        if (k > n) break;
    }
    free(cnt);
    /* rk becomes the inverse suffix array (0-based) */
    for (i = 0; i < n; i++) rk[sa[i]] = i;
}

static void kasai(const int *sa, const int *isa, int *lcp) {
    int n = N, h = 0, i;
    lcp[0] = 0;
    for (i = 0; i < n; i++) {
        int r = isa[i];
        if (r > 0) {
            int j = sa[r - 1];
            while (i + h < n && j + h < n && S[i + h] == S[j + h]) h++;
            lcp[r] = h;
            if (h > 0) h--;
        } else h = 0;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) die("usage: hard_lpf WORDFILE [mu] [w0 w1] [topK] [jmin]");
    double mu = argc > 2 ? atof(argv[2]) : 1.585;
    double w0 = argc > 4 ? atof(argv[3]) : 1.0, w1 = argc > 4 ? atof(argv[4]) : 1.0;
    int topK = argc > 5 ? atoi(argv[5]) : 12;
    int jmin = argc > 6 ? atoi(argv[6]) : 64;
    FILE *f = fopen(argv[1], "rb");
    if (!f) die("cannot open word file");
    fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
    unsigned char *raw = (unsigned char *)malloc((size_t)sz + 1);
    if (!raw || fread(raw, 1, (size_t)sz, f) != (size_t)sz) die("read");
    fclose(f);
    S = (unsigned char *)malloc((size_t)sz + 1);
    N = 0;
    for (long i = 0; i < sz; i++) if (raw[i] == '0' || raw[i] == '1') S[N++] = (unsigned char)(raw[i] - '0');
    free(raw);
    if (N < 16) die("word too short");
    int *sa = (int *)malloc((size_t)N * sizeof(int)), *rk = (int *)malloc((size_t)N * sizeof(int));
    int *tmp = (int *)malloc((size_t)N * sizeof(int)), *lcp = (int *)malloc((size_t)(N + 1) * sizeof(int));
    if (!sa || !rk || !tmp || !lcp) die("oom");
    build_sa(sa, rk, tmp);
    kasai(sa, rk, lcp);
    /* previous / next smaller text position in SA order, with lcp to it */
    int *lpf = (int *)calloc((size_t)N, sizeof(int)), *occ = (int *)malloc((size_t)N * sizeof(int));
    int *st = tmp; /* reuse */
    int *dd = (int *)malloc((size_t)N * sizeof(int));
    if (!lpf || !occ || !dd) die("oom2");
    for (int i = 0; i < N; i++) occ[i] = -1;
    const int INF = 0x3fffffff;
    int top = -1, cur = INF;
    for (int k = 0; k < N; k++) {
        if (k > 0 && lcp[k] < cur) cur = lcp[k];
        while (top >= 0 && sa[st[top]] > sa[k]) {
            if (top >= 1 && dd[top - 1] < cur) cur = dd[top - 1];
            top--;
        }
        if (top >= 0) {
            int v = cur == INF ? 0 : cur;
            if (v > lpf[sa[k]] || occ[sa[k]] < 0) { lpf[sa[k]] = v; occ[sa[k]] = sa[st[top]]; }
            dd[top] = cur;
        }
        st[++top] = k; cur = INF;
    }
    top = -1; cur = INF;
    for (int k = N - 1; k >= 0; k--) {
        if (k + 1 < N && lcp[k + 1] < cur) cur = lcp[k + 1];
        while (top >= 0 && sa[st[top]] > sa[k]) {
            if (top >= 1 && dd[top - 1] < cur) cur = dd[top - 1];
            top--;
        }
        if (top >= 0) {
            int v = cur == INF ? 0 : cur;
            if (v > lpf[sa[k]] || occ[sa[k]] < 0) { lpf[sa[k]] = v; occ[sa[k]] = sa[st[top]]; }
            dd[top] = cur;
        }
        st[++top] = k; cur = INF;
    }
    /* prefix weights */
    double *W = (double *)malloc((size_t)(N + 1) * sizeof(double));
    W[0] = 0.0;
    for (int i = 0; i < N; i++) W[i + 1] = W[i] + (S[i] ? w1 : w0);
    printf("N %d\n", N);
    /* per-scale profile */
    for (int s = 0; (1 << s) < N; s++) {
        int jlo = 1 << s, jhi = (1 << (s + 1)) < N ? (1 << (s + 1)) : N;
        double best = 0, wbest = 0, gbest = -1e300; int bj = -1, bwj = -1, bgj = -1, nun = 0;
        for (int j = jlo; j < jhi; j++) {
            if (j + lpf[j] >= N) continue; /* capped: repetition reaches the end of the window */
            nun++;
            double r = (double)(j + lpf[j]) / j;
            double wr = W[j + lpf[j]] / W[j];
            double g = lpf[j] - (mu - 1.0) * j;
            if (r > best) { best = r; bj = j; }
            if (wr > wbest) { wbest = wr; bwj = j; }
            if (g > gbest) { gbest = g; bgj = j; }
        }
        if (nun == 0) { printf("SCALE %d %d %d capped\n", s, jlo, jhi); continue; }
        printf("SCALE %d %d %d %.6f %d %d %d %.6f %d %d %d %.1f %d %d %d\n", s, jlo, jhi,
               best, bj, occ[bj], lpf[bj], wbest, bwj, occ[bwj], lpf[bwj], gbest, bgj, occ[bgj], lpf[bgj]);
    }
    /* top-K linear gains, at most one per dyadic scale first, then fill */
    {
        int *used = (int *)calloc(64, sizeof(int));
        int printed = 0;
        for (int round = 0; round < 3 && printed < topK; round++) {
            for (int s = 31; s >= 0 && printed < topK; s--) {
                int jlo = 1 << s; if (jlo >= N) continue;
                int jhi = (s < 30 && (1 << (s + 1)) < N) ? (1 << (s + 1)) : N;
                if (jhi <= jmin) continue;
                if (used[s] > round) continue;
                double gb = -1e300; int bj = -1;
                for (int j = (jlo > jmin ? jlo : jmin); j < jhi; j++) {
                    if (j + lpf[j] >= N) continue;
                    double g = lpf[j] - (mu - 1.0) * j;
                    if (round > 0 && g <= 0) continue;
                    if (g > gb) { gb = g; bj = j; }
                }
                if (bj < 0) continue;
                if (round > 0) { /* second/third best in the scale: exclude a neighbourhood */
                    double gb2 = -1e300; int bj2 = -1;
                    for (int j = (jlo > jmin ? jlo : jmin); j < jhi; j++) {
                        if (j + lpf[j] >= N || (j > bj - 64 && j < bj + 64)) continue;
                        double g = lpf[j] - (mu - 1.0) * j;
                        if (g > gb2) { gb2 = g; bj2 = j; }
                    }
                    if (bj2 < 0) continue; bj = bj2; gb = gb2;
                }
                used[s]++;
                printf("CAND %d %d %d %.1f\n", bj, occ[bj], lpf[bj], gb);
                printed++;
            }
        }
        free(used);
    }
    /* Bugeaud-Kim r(n): r(n) = n + min{ j : LPF[j] >= n } */
    {
        int nmax = N / 4; if (nmax > (1 << 20)) nmax = 1 << 20;
        int *mj = (int *)malloc((size_t)(nmax + 1) * sizeof(int));
        for (int n = 0; n <= nmax; n++) mj[n] = -1;
        int run = 0;
        for (int j = 1; j < N && run < nmax; j++) {
            int L = lpf[j]; if (L > nmax) L = nmax;
            if (L > run) { for (int n = run + 1; n <= L; n++) mj[n] = j; run = L; }
        }
        for (int n = 1; n <= nmax; n++) {
            if (mj[n] < 0) break;
            if (n <= 64 || (n & (n - 1)) == 0 || n % 997 == 0)
                printf("RN %d %d\n", n, n + mj[n]);
        }
        /* rep estimate on windows [n1, n2] */
        int wins[][2] = {{16, 256}, {256, 4096}, {4096, 65536}, {65536, 1 << 20}};
        for (int w = 0; w < 4; w++) {
            int n1 = wins[w][0], n2 = wins[w][1]; if (n2 > nmax) n2 = nmax;
            if (n1 >= n2) continue;
            double best = 1e300; int bn = -1, ok = 1;
            for (int n = n1; n <= n2; n++) {
                if (mj[n] < 0) { ok = 0; break; }
                /* r(n) must be an uncapped value: the occurrence pair lies inside the window */
                double r = (double)(n + mj[n]) / n;
                if (r < best) { best = r; bn = n; }
            }
            if (ok) printf("REP %d %d %.6f %d\n", n1, n2, best, bn);
            else printf("REP %d %d incomplete\n", n1, n2);
        }
        free(mj);
    }
    /* tail estimates of Dio over the second half of the usable range */
    {
        double fr[] = {0.125, 0.25, 0.5};
        for (int q = 0; q < 3; q++) {
            int j0 = (int)(fr[q] * N / 2);
            double best = 0, wbest = 0; int bj = -1;
            for (int j = j0; j < N / 2; j++) {
                if (j + lpf[j] >= N) continue;
                double r = (double)(j + lpf[j]) / j, wr = W[j + lpf[j]] / W[j];
                if (r > best) { best = r; bj = j; }
                if (wr > wbest) wbest = wr;
            }
            printf("DIOTAIL %.3f %d %.6f %.6f %d\n", fr[q], j0, best, wbest, bj);
        }
    }
    free(sa); free(rk); free(tmp); free(lcp); free(lpf); free(occ); free(dd); free(W); free(S);
    return 0;
}
