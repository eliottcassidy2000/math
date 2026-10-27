/* fshape.c -- exhaustive search over the self-converse shapes that contain every lambda = 1 exception (PROOF.md V.3):
 *   even: F(S, x0): bottom S, top sigma(S) carrying S^op, single upward arc x0 -> sigma(x0), all other pairs down.
 *   odd : G(S, x0): bottom S, one sigma-fixed middle vertex c, top sigma(S); upward arcs only x0 -> c -> sigma(x0).
 * For every strong S on a vertices (gentourng -cq a) and every x0: build T, test strong, lambda, delta, and whether
 * a black descent candidate exists (exact candidate enumeration, cand.h).  Prints every T without black descent
 * (with its nauty canonical string, to be matched against the auditor's own E list).
 * Usage: gentourng -cq a | ./fshape a even|odd
 */
#include "at_nauty.h"
#include "cand.h"

int main(int argc, char **argv) {
    int a = atoi(argv[1]); int odd = strcmp(argv[2], "odd") == 0;
    int n = 2 * a + odd;
    nauty_check(WORDSIZE, 1, n, NAUTYVERSIONID);
    char line[512]; tour S;
    long long built = 0, strong = 0, fails = 0, lamHist[8] = {0}, regime = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        line[strcspn(line, "\r\n")] = 0;
        t_parse(line, a, &S);
        for (int x0 = 0; x0 < a; x0++) {
            tour T; t_clear(&T, n);
            int top0 = a + odd;            /* sigma(u) = top0 + u */
            int layer[MAXV];
            for (int u = 0; u < a; u++) { layer[u] = 0; layer[top0 + u] = 2; }
            if (odd) layer[a] = 1;
            for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (layer[u] != layer[v]) { if (layer[u] > layer[v]) t_arc(&T, u, v); else t_arc(&T, v, u); }
            for (int u = 0; u < a; u++) for (int v = 0; v < a; v++) if (u != v && t_has(&S, u, v)) { t_arc(&T, u, v); t_arc(&T, top0 + v, top0 + u); }
            if (odd) { t_arc(&T, x0, a); t_arc(&T, a, top0 + x0); } else t_arc(&T, x0, top0 + x0);
            t_fin(&T); built++;
            if (!t_valid(&T)) { printf("invalid\n"); return 1; }
            /* self-converse by construction: check the involution */
            int sg[MAXV]; for (int u = 0; u < a; u++) { sg[u] = top0 + u; sg[top0 + u] = u; } if (odd) sg[a] = a;
            if (!t_is_anti(&T, sg)) { printf("sigma not anti\n"); return 1; }
            if (!t_strong_mask(&T, fullmask(n))) continue;
            strong++;
            int Q[MAXV], lam; long long nc, nb;
            int found = black_descent(&T, 0, Q, &lam, &nc, &nb);
            int del = t_delta(&T);
            lamHist[lam]++;
            if (lam < del) regime++;
            if (!found) {
                fails++;
                char code[512], cs[400]; t_code(&T, code); nt_canon_str(&T, cs);
                printf("NOBLACK n=%d S=%s x0=%d lambda=%d delta=%d cand=%lld canon=%s\n", n, line, x0, lam, del, nc, cs);
            }
        }
    }
    printf("a=%d n=%d shape=%s built=%lld strong=%lld lambda<delta=%lld no-black-descent=%lld lambda-hist:", a, n, argv[2], built, strong, regime, fails);
    for (int l = 0; l < 8; l++) if (lamHist[l]) printf(" %d:%lld", l, lamHist[l]);
    printf("\n");
    return 0;
}
