/* thmCrand.c -- random control of Theorem C beyond the census: random self-converse tournaments (random orientation
 * of every orbit of vertex pairs under a fixed involution sigma with at most one fixed point), kept if strong and
 * nonregular; for lambda = delta we require a black descent of singleton type (s(q_1) = lambda or d^-(q_n) = lambda);
 * for lambda < delta we record whether a black descent exists (Theorem A) and, if not, whether T is in E (nauty).
 * Usage: ./thmCrand n samples seed
 */
#include "at_nauty.h"
#include "cand.h"
#include "fam.h"
static uint64_t rs = 0x9E3779B97F4A7C15ull;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int main(int argc, char **argv) {
    int n = atoi(argv[1]); long N = atol(argv[2]); rs ^= (uint64_t)atoll(argv[3]) * 0x2545F4914F6CDD1Dull;
    nauty_check(WORDSIZE, 1, n, NAUTYVERSIONID);
    static char Ecanon[400][420]; int nE = 0; int parts[16];
    for (int aa = 5; aa <= n / 2; aa++) gen_comp(aa - 4, parts, 0, aa);
    for (int bi = 0; bi < nBL; bi++) { tour S0; int bs[17]; build_S(&BL[bi], &S0, bs); int t = n - 2 * S0.n + 1; if (t >= 1) { tour T0; int s0[MAXV]; build_E(&S0, t, &T0, s0); nt_canon_str(&T0, Ecanon[nE++]); } }
    int k = n / 2, sg[MAXV];
    for (int i = 0; i < k; i++) { sg[i] = k + i; sg[k + i] = i; } if (n % 2) sg[n - 1] = n - 1;
    long kept = 0, eq = 0, failC = 0, lt = 0, ltNoBlack = 0, ltNoBlackNotE = 0;
    for (long it = 0; it < N; it++) {
        tour T; t_clear(&T, n);
        /* bias: with probability 1/2 use a random 'layered' bias towards lambda < delta structures */
        static int done[MAXV][MAXV]; memset(done, 0, sizeof done);
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) {
            if (done[u][v]) continue;
            int a = sg[u], b = sg[v];
            done[u][v] = done[v][u] = done[a][b] = done[b][a] = 1;
            int dir = rnd() & 1;
            if (dir) { t_arc(&T, u, v); t_arc(&T, b, a); } else { t_arc(&T, v, u); t_arc(&T, a, b); }
        }
        t_fin(&T);
        if (!t_is_anti(&T, sg)) { printf("construction error\n"); return 1; }
        if (!t_strong_mask(&T, fullmask(n)) || t_regular(&T)) continue;
        kept++;
        int lam, Q[MAXV]; long long nc, nb;
        static uint32_t cuts[1 << 17]; int ncut; lam = t_mincuts(&T, cuts, 1 << 17, &ncut);
        int del = t_delta(&T);
        if (lam == del) {
            eq++;
            if (!black_descent(&T, 1, Q, &lam, &nc, &nb)) { failC++; char code[512]; t_code(&T, code); printf("THEOREM-C-FAIL n=%d %s\n", n, code); }
        } else {
            lt++;
            if (!black_descent(&T, 0, Q, &lam, &nc, &nb)) {
                ltNoBlack++; char cs[420]; nt_canon_str(&T, cs); int in = 0; for (int i = 0; i < nE; i++) if (!strcmp(cs, Ecanon[i])) in = 1;
                if (!in) { ltNoBlackNotE++; char code[512]; t_code(&T, code); printf("NEW-EXCEPTION n=%d %s\n", n, code); }
            }
        }
    }
    printf("n=%d samples=%ld kept(strong nonregular self-converse)=%ld lambda=delta:%ld TheoremC-fail:%ld  lambda<delta:%ld no-black:%ld not-in-E:%ld\n",
           n, N, kept, eq, failC, lt, ltNoBlack, ltNoBlackNotE);
    return 0;
}
