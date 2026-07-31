/*
 * sens_search.c — exact decision: does there exist f:{0,1}^n -> {0,1},
 * multilinear real degree <= d, with f(0)=0 and f(e_i)=1 for ALL i?
 *
 * Motivation (Epoch FrontierMath "degree vs sensitivity"):
 *   restricting a Boolean P (on any number of variables) to the coordinates
 *   sensitive at 0 preserves degree<=d, Booleanness, and the values at 0/e_i,
 *   so WLOG n = s0.  Hence:
 *     SAT(n,d)   => explicit P with s0(P) = n = deg^{log_d n}
 *     UNSAT(n,d) => theorem: no degree-<=d Boolean function has n
 *                   sensitive coordinates at an input (max 0-sensitivity bound).
 *   Known record: Kushilevitz (d,s)=(3,6), a = log 6/log 3 ~ 1.6309.
 *   Target: (3,7) a~1.7712, (4,10) a~1.6610, (5,14) a~1.6397.
 *
 * Method: masks processed in increasing numeric (colex-compatible) order.
 *   |m|<=1: fixed.  2<=|m|<=d: free bit (branch).  |m|>d: value forced by
 *   vanishing of the Moebius coefficient: c_m = sum_{T<=m} (-1)^{|m|-|T|} f(T) = 0,
 *   and every proper submask of m is numerically smaller => already assigned.
 *   Prune when forced value is not in {0,1}.  Exhaustive => UNSAT is a proof.
 *
 * Usage: ./sens_search n d [max_solutions_to_print] [stop_after_first]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, d, NPOW;
static signed char f[1 << 20];
static unsigned long long nodes = 0;
static unsigned long long sols = 0;
static int max_print = 10, stop_first = 0;
static int done = 0;

static inline int pc(int x) { return __builtin_popcount(x); }

static void print_solution(void) {
    sols++;
    if (sols <= (unsigned)max_print) {
        printf("SOLUTION #%llu\n", sols);
        /* truth table */
        printf("  truth table (mask:val): ");
        for (int m = 0; m < NPOW; m++) if (f[m]) printf("%d ", m);
        printf("\n  polynomial: ");
        int first = 1, maxdeg = 0;
        for (int S = 0; S < NPOW; S++) {
            if (pc(S) > d) continue;
            /* Moebius coefficient c_S = sum_{T subset S} (-1)^{|S|-|T|} f(T) */
            long c = 0;
            int T = S;
            for (;;) {
                c += ((pc(S) - pc(T)) & 1) ? -f[T] : f[T];
                if (T == 0) break;
                T = (T - 1) & S;
            }
            if (c != 0) {
                if (pc(S) > maxdeg) maxdeg = pc(S);
                if (!first) printf(" + ");
                first = 0;
                if (S == 0) { printf("%ld", c); }
                else {
                    if (c != 1) printf("%ld*", c);
                    int firstv = 1;
                    for (int i = 0; i < n; i++) if (S & (1 << i)) {
                        if (!firstv) printf("*");
                        firstv = 0;
                        printf("x%d", i + 1);
                    }
                }
            }
        }
        printf("\n  actual degree: %d\n", maxdeg);
        fflush(stdout);
    }
}

static void solve(int m) {
    if (done) return;
    if (m == NPOW) {
        print_solution();
        if (stop_first) done = 1;
        return;
    }
    nodes++;
    if ((nodes & ((1ULL << 31) - 1)) == 0) {
        fprintf(stderr, "  ... %llu nodes, %llu sols, at mask %d\n", nodes, sols, m);
    }
    int k = pc(m);
    if (k <= 1) {           /* fixed: f(0)=0, f(e_i)=1 */
        f[m] = (k == 1);
        solve(m + 1);
        return;
    }
    if (k <= d) {           /* free bit */
        f[m] = 0; solve(m + 1);
        if (done) return;
        f[m] = 1; solve(m + 1);
        return;
    }
    /* forced: f(m) = -sum_{T proper subset m} (-1)^{|m|-|T|} f(T) */
    long s = 0;
    int T = (m - 1) & m;
    for (;;) {
        s += ((k - pc(T)) & 1) ? -f[T] : f[T];
        if (T == 0) break;
        T = (T - 1) & m;
    }
    long v = -s;
    if (v == 0 || v == 1) {
        f[m] = (signed char)v;
        solve(m + 1);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s n d [max_print] [stop_first]\n", argv[0]); return 1; }
    n = atoi(argv[1]); d = atoi(argv[2]);
    if (argc > 3) max_print = atoi(argv[3]);
    if (argc > 4) stop_first = atoi(argv[4]);
    NPOW = 1 << n;
    fprintf(stderr, "sens_search: n=%d d=%d (free bits: ", n, d);
    long fb = 0;
    for (int k = 2; k <= d && k <= n; k++) {
        /* C(n,k) */
        long c = 1;
        for (int i = 0; i < k; i++) c = c * (n - i) / (i + 1);
        fb += c;
    }
    fprintf(stderr, "%ld)\n", fb);
    solve(0);
    printf("RESULT n=%d d=%d: %s  (solutions=%llu, nodes=%llu)%s\n",
           n, d, sols ? "SAT" : "UNSAT", sols, nodes,
           (done ? " [stopped at first solution]" : " [exhaustive]"));
    return 0;
}
