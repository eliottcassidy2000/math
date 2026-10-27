/* twoorbit.c -- targeted test of the l' = 2 branch of the last-block lemma (PROOF.md Part VI (c)):
 * m = 2p, Hamiltonian cycle C = (u_0, w_0, u_1, w_1, ...), rho^2 an automorphism, B = one rho^2-orbit, A = the other.
 * All configurations: skew connection sets S_U, S_W of Z_p, cross set S_X (u_i -> w_(i+d) iff d in S_X) with
 * 0 in S_X, -1 notin S_X.  For B = U-orbit and for B = W-orbit, test (H1) (G = Aut(R;A) transitive on Ham(R -> B),
 * via |G| = |H'| -- G acts freely, LEM-003 -- and also directly by the cycle argument (a): does R have a Ham cycle
 * with two consecutive B vertices?) and (H2).
 * Usage: ./twoorbit p
 */
#include "at_nauty.h"

int main(int argc, char **argv) {
    int p = atoi(argv[1]); int m = 2 * p;
    nauty_check(WORDSIZE, 1, m, NAUTYVERSIONID);
    int half = (p - 1) / 2;
    long long configs = 0, h1 = 0, h1h2 = 0, consecCycle = 0, h1butConsec = 0;
    /* skew set chosen by sign vector over d = 1..half: d in S or -d in S */
    for (int su = 0; su < (1 << half); su++) for (int sw = 0; sw < (1 << half); sw++) for (int sx = 0; sx < (1 << p); sx++) {
        if (!(sx & 1)) continue;                 /* 0 in S_X */
        if (sx >> (p - 1) & 1) continue;         /* -1 notin S_X */
        tour R; t_clear(&R, m);
        /* u_i = 2i, w_i = 2i+1 */
        for (int i = 0; i < p; i++) for (int d = 1; d <= half; d++) {
            int j = (i + d) % p;
            if (su >> (d - 1) & 1) t_arc(&R, 2 * i, 2 * j); else t_arc(&R, 2 * j, 2 * i);
            if (sw >> (d - 1) & 1) t_arc(&R, 2 * i + 1, 2 * j + 1); else t_arc(&R, 2 * j + 1, 2 * i + 1);
        }
        for (int i = 0; i < p; i++) for (int d = 0; d < p; d++) {
            int j = (i + d) % p;
            if (sx >> d & 1) t_arc(&R, 2 * i, 2 * j + 1); else t_arc(&R, 2 * j + 1, 2 * i);
        }
        t_fin(&R);
        if (!t_valid(&R)) { printf("invalid\n"); return 1; }
        /* sanity: C is a Ham cycle, rho^2 automorphism */
        int rho2[MAXV]; for (int v = 0; v < m; v++) rho2[v] = (v + 2) % m;
        if (!t_is_aut(&R, rho2)) { printf("rho2 not aut\n"); return 1; }
        for (int v = 0; v < m; v++) if (!t_has(&R, v, (v + 1) % m)) { printf("C not a cycle\n"); return 1; }
        configs++;
        for (int which = 0; which < 2; which++) {
            uint32_t U = 0; for (int i = 0; i < p; i++) U |= 1u << (2 * i);
            uint32_t B = which ? (fullmask(m) & ~U) : U, A = fullmask(m) & ~B;
            int col[MAXV]; for (int v = 0; v < m; v++) col[v] = (A >> v) & 1u;
            double g = nt_autsize_col(&R, col);
            long long hB = t_count_hp_end(&R, fullmask(m), B);
            int H1 = ((double)hB == g);
            /* (H2) */
            int H2 = 1; for (uint32_t W = 1; W <= fullmask(m) && H2; W++) if (t_f(&R, W) + pc(W & B) < 2) H2 = 0;
            /* exists Ham cycle with two consecutive B vertices?  = Ham path from b to b' (both in B) ... closing b' -> b?
               consecutive b -> b2 on a Ham cycle: a Ham path starting at b2 and ending at b with b -> b2 */
            int consec = 0;
            for (int b = 0; b < m && !consec; b++) if (B >> b & 1u) for (int b2 = 0; b2 < m && !consec; b2++) if ((B >> b2 & 1u) && b2 != b && t_has(&R, b, b2)) {
                /* Ham path from b2 to b: count paths of R ending at b among those starting at b2 -- use DFS */
                int P[MAXV]; uint32_t endset = 1u << b;
                if (t_hp_in(&R, fullmask(m), b2, endset, P)) consec = 1;
            }
            if (H1) h1++;
            if (H1 && H2) h1h2++;
            if (consec) consecCycle++;
            if (H1 && consec) h1butConsec++;
            if (H1) {
                char code[256]; t_code(&R, code);
                printf("p=%d S_U=%x S_W=%x S_X=%x B=%s: (H1) holds |G|=%g |H'|=%lld (H2)=%d consecutive-B Ham cycle=%d  %s\n", p, su, sw, sx, which ? "W" : "U", g, hB, H2, consec, code);
            }
        }
    }
    printf("p=%d m=%d configurations=%lld (x2 choices of B): (H1)=%lld (H1)+(H2)=%lld  consecutive-B cycle exists=%lld  (H1) with consecutive-B cycle=%lld\n",
           p, m, configs, h1, h1h2, consecCycle, h1butConsec);
    return 0;
}
