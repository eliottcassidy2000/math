/* lam2gen.c -- general exhaustive search of self-converse configurations in the regime 2 <= lambda < delta with
 * |L_0| = a and k middle vertices (auditor's code; superset of every exception of the regime with |X_1| = a).
 * Structure: bottom S = T[L_0] (strong, a vertices), top sigma(L_0) carrying S^op, middles M (k vertices) with sigma
 * acting as one fixed point (k odd) plus swapped pairs; every orientation among middles compatible with sigma; up-arcs
 * U from L_0 to M and F from L_0 to sigma(L_0) (x -> sigma(u), F symmetric) with |U| + |F| = 2 (so f(L_0) = 2); all
 * images under sigma; every other pair points down (top -> M -> bottom, top -> bottom).
 * Keeps T strong with 2 <= lambda < delta and tests for a black descent candidate.
 * Usage: gentourng -cq a | ./lam2gen a k
 */
#include "at_nauty.h"
#include "cand.h"

static long long built = 0, regime = 0, fails = 0, lamDelta[8][8];
static void test(const tour *T, const int *sg) {
    built++;
    if (!t_valid(T)) { printf("invalid\n"); exit(1); }
    if (!t_is_anti(T, sg)) { printf("sigma not anti\n"); exit(1); }
    int del = t_delta(T);
    if (del < 3) return;
    static uint32_t cuts[1 << 17]; int nc;
    int lam = t_mincuts(T, cuts, 1 << 17, &nc);
    if (!(lam >= 2 && lam < del)) return;
    regime++; lamDelta[lam][del]++;
    int Q[MAXV], l2; long long ncand, nblue;
    if (!black_descent(T, 0, Q, &l2, &ncand, &nblue)) { fails++; char code[512]; t_code(T, code); printf("NOBLACK lambda=%d delta=%d ncuts=%d cand=%lld %s\n", lam, del, nc, ncand, code); }
}

int main(int argc, char **argv) {
    int a = atoi(argv[1]), k = atoi(argv[2]); int n = 2 * a + k;
    nauty_check(WORDSIZE, 1, n, NAUTYVERSIONID);
    int top0 = a + k;
    int sg[MAXV];
    for (int u = 0; u < a; u++) { sg[u] = top0 + u; sg[top0 + u] = u; }
    int j = a;
    if (k % 2) { sg[a] = a; j = a + 1; }
    for (; j + 1 < a + k; j += 2) { sg[j] = j + 1; sg[j + 1] = j; }
    /* orbits of unordered middle pairs under sigma: choose a representative orientation bit per orbit */
    int po_u[64], po_v[64], npo = 0; int seen[MAXV][MAXV]; memset(seen, 0, sizeof seen);
    for (int u = a; u < a + k; u++) for (int v = u + 1; v < a + k; v++) {
        if (seen[u][v]) continue;
        int su = sg[u], sv = sg[v];
        seen[u][v] = seen[v][u] = 1; seen[su][sv] = seen[sv][su] = 1;
        po_u[npo] = u; po_v[npo] = v; npo++;
    }
    /* possible up-arcs from L_0 to M: list */
    int ux[256], um[256], nu = 0; for (int x = 0; x < a; x++) for (int m = a; m < a + k; m++) { ux[nu] = x; um[nu] = m; nu++; }
    char line[512]; tour S;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        t_parse(line, a, &S);
        for (int mb = 0; mb < (1 << npo); mb++) {
            tour B; t_clear(&B, n);
            int layer[MAXV]; for (int u = 0; u < a; u++) { layer[u] = 0; layer[top0 + u] = 2; } for (int m = a; m < a + k; m++) layer[m] = 1;
            for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (layer[u] != layer[v]) { if (layer[u] > layer[v]) t_arc(&B, u, v); else t_arc(&B, v, u); }
            for (int u = 0; u < a; u++) for (int v = 0; v < a; v++) if (u != v && t_has(&S, u, v)) { t_arc(&B, u, v); t_arc(&B, top0 + v, top0 + u); }
            for (int p = 0; p < npo; p++) {
                int u = po_u[p], v = po_v[p];
                if (mb >> p & 1) { t_arc(&B, u, v); t_arc(&B, sg[v], sg[u]); } else { t_arc(&B, v, u); t_arc(&B, sg[u], sg[v]); }
            }
            /* up-arc choices: |U| + |F| = 2 */
            /* |U| = 2 */
            for (int i1 = 0; i1 < nu; i1++) for (int i2 = i1 + 1; i2 < nu; i2++) {
                tour T = B;
                int L[2] = {i1, i2};
                for (int q = 0; q < 2; q++) { int x = ux[L[q]], m = um[L[q]]; t_arc(&T, x, m); t_arc(&T, sg[m], sg[x]); }
                t_fin(&T); test(&T, sg);
            }
            /* |U| = 1, |F| = 1 (diagonal) */
            for (int i1 = 0; i1 < nu; i1++) for (int y = 0; y < a; y++) {
                tour T = B; int x = ux[i1], m = um[i1]; t_arc(&T, x, m); t_arc(&T, sg[m], sg[x]);
                t_arc(&T, y, top0 + y); t_fin(&T); test(&T, sg);
            }
            /* |F| = 2 */
            for (int y = 0; y < a; y++) for (int z = y + 1; z < a; z++) {
                tour T = B; t_arc(&T, y, top0 + y); t_arc(&T, z, top0 + z); t_fin(&T); test(&T, sg);
                tour T2 = B; t_arc(&T2, y, top0 + z); t_arc(&T2, z, top0 + y); t_fin(&T2); test(&T2, sg);
            }
        }
    }
    printf("a=%d k=%d n=%d built=%lld in-regime(2<=lambda<delta)=%lld no-black-descent=%lld  (lambda,delta) counts:", a, k, n, built, regime, fails);
    for (int l = 0; l < 8; l++) for (int d = 0; d < 8; d++) if (lamDelta[l][d]) printf(" (%d,%d):%lld", l, d, lamDelta[l][d]);
    printf("\n");
    return 0;
}
