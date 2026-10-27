/* lam1gen.c -- EXHAUSTIVE search of the regime lambda = 1 < delta for a given (a, k): every self-converse T whose
 * smallest minimum cut L_0 has a vertices, with k further vertices between L_0 and L_t = sigma(L_0).
 * Justification (PROOF.md Part I, F11, III.1-III.3): an exception is self-converse with an involutive anti-automorphism
 * sigma; the minimum cuts form a chain reversed by sigma, so sigma(L_0) = L_t, and the k = n - 2a remaining vertices
 * ("middles") are permuted by sigma with at most one fixed point.  Every arc from L_0 to the middles or to L_t counts in
 * f(L_0) = 1, and by sigma every arc from the middles or L_0 into L_t counts in f(V \ L_t) = 1.  So T is determined by:
 * S = T[L_0] (strong), a sigma-compatible orientation of the middles, and ONE up-arc out of L_0 (x -> m, with its image
 * sigma(m) -> sigma(x), or x -> sigma(x)); all other pairs point down.  We enumerate all of these (superset; T need not
 * have L_0 as its smallest cut), keep strong T with lambda = 1 < delta, test for a black descent candidate, and compare
 * every exception with the auditor's own E list (nauty canonical forms).
 * Usage: gentourng -cq a | ./lam1gen a k
 */
#include "at_nauty.h"
#include "cand.h"
#include "fam.h"

static char Ecanon[4000][420]; static int nE = 0;
static long long built = 0, regime = 0, fails = 0, failsE = 0, failsNew = 0;
static void test(const tour *T, const int *sg) {
    built++;
    if (!t_valid(T) || !t_is_anti(T, sg)) { printf("construction error\n"); exit(1); }
    int del = t_delta(T); if (del < 2) return;
    if (!t_strong_mask(T, fullmask(T->n))) return;
    static uint32_t cuts[1 << 17]; int nc; int lam = t_mincuts(T, cuts, 1 << 17, &nc);
    if (!(lam == 1 && lam < del)) return;
    regime++;
    int Q[MAXV], l2; long long ncand, nblue;
    if (!black_descent(T, 0, Q, &l2, &ncand, &nblue)) {
        fails++; char cs[420]; nt_canon_str(T, cs);
        int in = 0; for (int i = 0; i < nE; i++) if (!strcmp(Ecanon[i], cs)) in = 1;
        if (in) failsE++; else { failsNew++; char code[512]; t_code(T, code); printf("NEW-EXCEPTION n=%d delta=%d ncuts=%d %s\n", T->n, del, nc, code); }
    }
}
int main(int argc, char **argv) {
    int a = atoi(argv[1]), k = atoi(argv[2]); int n = 2 * a + k;
    nauty_check(WORDSIZE, 1, n, NAUTYVERSIONID);
    /* own E list at order n */
    int parts[16]; for (int aa = 5; aa <= n / 2; aa++) gen_comp(aa - 4, parts, 0, aa);
    for (int bi = 0; bi < nBL; bi++) { tour S0; int bs[17]; build_S(&BL[bi], &S0, bs); int t = n - 2 * S0.n + 1; if (t >= 1) { tour T0; int s0[MAXV]; build_E(&S0, t, &T0, s0); nt_canon_str(&T0, Ecanon[nE++]); } }
    int top0 = a + k; int sg[MAXV];
    for (int u = 0; u < a; u++) { sg[u] = top0 + u; sg[top0 + u] = u; }
    int j = a; if (k % 2) { sg[a] = a; j = a + 1; } for (; j + 1 < a + k; j += 2) { sg[j] = j + 1; sg[j + 1] = j; }
    int po_u[64], po_v[64], npo = 0; static int seen[MAXV][MAXV];
    for (int u = a; u < a + k; u++) for (int v = u + 1; v < a + k; v++) {
        if (seen[u][v]) continue; seen[u][v] = seen[v][u] = 1; seen[sg[u]][sg[v]] = seen[sg[v]][sg[u]] = 1;
        po_u[npo] = u; po_v[npo] = v; npo++;
    }
    char line[512]; tour S;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        t_parse(line, a, &S);
        for (long mb = 0; mb < (1L << npo); mb++) {
            tour B; t_clear(&B, n);
            int layer[MAXV]; for (int u = 0; u < a; u++) { layer[u] = 0; layer[top0 + u] = 2; } for (int m = a; m < a + k; m++) layer[m] = 1;
            for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (layer[u] != layer[v]) { if (layer[u] > layer[v]) t_arc(&B, u, v); else t_arc(&B, v, u); }
            for (int u = 0; u < a; u++) for (int v = 0; v < a; v++) if (u != v && t_has(&S, u, v)) { t_arc(&B, u, v); t_arc(&B, top0 + v, top0 + u); }
            for (int p = 0; p < npo; p++) { int u = po_u[p], v = po_v[p]; if (mb >> p & 1) { t_arc(&B, u, v); t_arc(&B, sg[v], sg[u]); } else { t_arc(&B, v, u); t_arc(&B, sg[u], sg[v]); } }
            for (int x = 0; x < a; x++) {
                for (int m = a; m < a + k; m++) { tour T = B; t_arc(&T, x, m); t_arc(&T, sg[m], sg[x]); t_fin(&T); test(&T, sg); }
                { tour T = B; t_arc(&T, x, top0 + x); t_fin(&T); test(&T, sg); }
            }
        }
    }
    printf("a=%d k=%d n=%d E-members-at-n=%d built=%lld in-regime(lambda=1<delta)=%lld no-black-descent=%lld (in E %lld, NEW %lld)\n",
           a, k, n, nE, built, regime, fails, failsE, failsNew);
    return 0;
}
