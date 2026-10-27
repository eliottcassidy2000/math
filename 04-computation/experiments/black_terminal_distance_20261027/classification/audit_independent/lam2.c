/* lam2.c -- exhaustive search of the lambda = 2 < delta regime at small orders (Theorem D control), auditor's code.
 * Every exception is self-converse (Part I) with an involutive anti-automorphism sigma (F11) reversing the chain of
 * minimum cuts (III.3), and |L_0| = |L_t| >= 2*delta+1 >= 7 (III.1).
 * mode t1  (n = 2a):   bottom S (a vertices), top sigma(S) with S^op; forward arcs x -> sigma(u) for (x,u) in a
 *                      sigma-symmetric relation F on S with |F| = 2; all other cross pairs point down.
 *                      At n = 14 (a = 7) this is every configuration of the regime (t = 1 forced).
 * mode mid (n = 2a+1): bottom S, sigma-fixed middle c, top sigma(S).  Forward arcs: x -> c and c -> sigma(x) for x in A,
 *                      and x -> sigma(u) for (x,u) in symmetric F, with |A| + |F| = 2.  At n = 15 (a = 7) this is every
 *                      configuration of the regime (|L_0| = 7, one middle vertex, fixed by sigma).
 * Keeps T with 2 <= lambda < delta, strong; tests for a black descent candidate (exact enumeration).
 * Usage: gentourng -cq a | ./lam2 a t1|mid
 */
#include "at_nauty.h"
#include "cand.h"

static long long built = 0, regime = 0, fails = 0, lamDelta[8][8];
static void test(const tour *T, const int *sg) {
    built++;
    if (!t_valid(T)) { printf("invalid\n"); exit(1); }
    if (!t_is_anti(T, sg)) { printf("sigma not anti\n"); exit(1); }
    int lam; static uint32_t cuts[1 << 16]; int nc;
    int del = t_delta(T);
    if (del < 3) return;                             /* need delta > lambda >= 2 */
    lam = t_mincuts(T, cuts, 1 << 16, &nc);
    if (!(lam >= 2 && lam < del)) return;
    regime++; lamDelta[lam][del]++;
    int Q[MAXV], l2; long long ncand, nblue;
    int found = black_descent(T, 0, Q, &l2, &ncand, &nblue);
    if (!found) { fails++; char code[512]; t_code(T, code); printf("NOBLACK lambda=%d delta=%d ncuts=%d cand=%lld %s\n", lam, del, nc, ncand, code); }
}

int main(int argc, char **argv) {
    int a = atoi(argv[1]); int mid = strcmp(argv[2], "mid") == 0; int n = 2 * a + mid;
    nauty_check(WORDSIZE, 1, n, NAUTYVERSIONID);
    char line[512]; tour S;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        t_parse(line, a, &S);
        int top0 = a + mid;
        /* base tournament: S, top = S^op, middle c = a (if mid), all cross pairs down */
        tour B; t_clear(&B, n);
        int layer[MAXV]; for (int u = 0; u < a; u++) { layer[u] = 0; layer[top0 + u] = 2; } if (mid) layer[a] = 1;
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (layer[u] != layer[v]) { if (layer[u] > layer[v]) t_arc(&B, u, v); else t_arc(&B, v, u); }
        for (int u = 0; u < a; u++) for (int v = 0; v < a; v++) if (u != v && t_has(&S, u, v)) { t_arc(&B, u, v); t_arc(&B, top0 + v, top0 + u); }
        int sg[MAXV]; for (int u = 0; u < a; u++) { sg[u] = top0 + u; sg[top0 + u] = u; } if (mid) sg[a] = a;
        /* symmetric relations F with |F| = k, as lists of unordered "orbits": diagonal (x,x) or pair {(x,u),(u,x)} */
        int ndiag = a, npair = a * (a - 1) / 2;
        int px[64], pu[64]; int k = 0; for (int x = 0; x < a; x++) for (int u = x + 1; u < a; u++) { px[k] = x; pu[k] = u; k++; }
        /* F of size 2: two diagonals, or one pair */
        for (int fsz = 0; fsz <= 2; fsz++) {
            int asz = 2 - fsz;
            if (!mid && asz) continue;
            /* enumerate A subsets of size asz */
            for (uint32_t Am = 0; Am < (1u << a); Am++) {
                if (pc(Am) != asz) continue;
                if (fsz == 0) {
                    tour T = B; for (int x = 0; x < a; x++) if (Am >> x & 1u) { t_arc(&T, x, a); t_arc(&T, a, top0 + x); }
                    t_fin(&T); test(&T, sg);
                } else if (fsz == 1) {
                    for (int y = 0; y < ndiag; y++) {
                        tour T = B; for (int x = 0; x < a; x++) if (Am >> x & 1u) { t_arc(&T, x, a); t_arc(&T, a, top0 + x); }
                        t_arc(&T, y, top0 + y); t_fin(&T); test(&T, sg);
                    }
                } else {
                    for (int y = 0; y < ndiag; y++) for (int z = y + 1; z < ndiag; z++) {
                        tour T = B; t_arc(&T, y, top0 + y); t_arc(&T, z, top0 + z); t_fin(&T); test(&T, sg);
                    }
                    for (int p = 0; p < npair; p++) {
                        tour T = B; t_arc(&T, px[p], top0 + pu[p]); t_arc(&T, pu[p], top0 + px[p]); t_fin(&T); test(&T, sg);
                    }
                }
            }
        }
    }
    printf("a=%d n=%d mode=%s built=%lld in-regime(2<=lambda<delta)=%lld no-black-descent=%lld   (lambda,delta) counts:", a, n, argv[2], built, regime, fails);
    for (int l = 0; l < 8; l++) for (int d = 0; d < 8; d++) if (lamDelta[l][d]) printf(" (%d,%d):%lld", l, d, lamDelta[l][d]);
    printf("\n");
    return 0;
}
