// lastblock_check.cpp -- finite check for Lemma BL(c) (the last block) used in PROOF.md, Part VI:
// for every strong tournament R on m vertices (gentourng -cq m) and every subset A of V(R) (B = V(R) \ A nonempty),
// test  (H1) the group G = Aut(R; A) (automorphisms preserving A) acts transitively on the Hamiltonian paths of R
//            ending in B  (equivalently h_B = |G|, the action being free), and
//       (H2) every nonempty W subset of V(R) has e_R(W -> V(R)\W) + |W cap B| >= 2.
// The lemma asserts that (H1)+(H2) force A = empty and R = C3.  The proof covers all m except the sub-case
// m = 6 (p = 3, l' = 2), which this program settles; it also re-checks every m <= 7 directly.
// Usage: gentourng -cq m | ./lastblock_check m
#include "tour.hpp"
struct AutA {
    const Tour *T; uint32_t A; int sig[20]; bool used[20]; long long cnt;
    void rec(int k) {
        int n = T->n;
        if (k == n) { cnt++; return; }
        for (int c = 0; c < n; c++) {
            if (used[c] || T->s[c] != T->s[k]) continue;
            if (((A >> k) & 1) != ((A >> c) & 1)) continue;
            bool ok = true;
            for (int j = 0; j < k && ok; j++) if (T->arc(k, j) != T->arc(c, sig[j])) ok = false;
            if (!ok) continue;
            sig[k] = c; used[c] = true; rec(k + 1); used[c] = false;
        }
    }
    long long run(const Tour &TT, uint32_t AA) { T = &TT; A = AA; memset(used, 0, sizeof used); cnt = 0; rec(0); return cnt; }
};
int main(int argc, char **argv) {
    int m = atoi(argv[1]); char line[256]; Tour R; long long pairs = 0, h2 = 0, both = 0; uint32_t full = (1u << m) - 1;
    while (fgets(line, sizeof line, stdin)) {
        if (!parseLine(line, m, R)) continue;
        for (uint32_t A = 0; A < full; A++) {           // B = full \ A nonempty
            uint32_t B = full & ~A; pairs++;
            bool ok = true;
            for (uint32_t W = 1; W <= full && ok; W++) if (fOut(R, W) + __builtin_popcount(W & B) < 2) ok = false;
            if (!ok) continue;
            h2++;
            long long hB = 0;
            for (int s0 = 0; s0 < m; s0++) for (int e = 0; e < m; e++) if (e != s0 && (B >> e & 1)) pathsAB(R, s0, e, [&](const int *) { hB++; return false; });
            long long g = AutA().run(R, A);
            if (hB == g) { both++; printf("  (H1)+(H2): R=%s A=%x |G|=%lld h_B=%lld %s\n", encode(R).c_str(), A, g, hB, (A == 0 && m == 3) ? "(C3, A empty: allowed)" : "COUNTEREXAMPLE"); }
        }
    }
    printf("m=%d: (R,A) pairs=%lld  satisfying (H2)=%lld  satisfying (H1)+(H2)=%lld\n", m, pairs, h2, both);
    return 0;
}
