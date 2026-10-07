/* Audit B (independent) of THM-4602 (2): exhaustive equivalence (a) <=> (b) <=> (c) for n <= 9.
 *
 *  B_n : all labelled tournaments with every 4-Pfaffian = +-1, enumerated by backtracking over arcs
 *        (vertex by vertex), pruning a branch as soon as a fully decided 4-set has |Pf| = 3.
 *        Each member is checked for (c) local transitivity BY DEFINITION and for (a) by a constructive
 *        integer-vector certificate (switch so that 0 is a source; transitive => u_i = eps_i (1, r_i)).
 *  C_n : all labelled locally transitive tournaments, enumerated independently by vertex extension
 *        (LT is hereditary), checking LT by definition (out/in-neighbourhood scores distinct).
 *  Expect |B_n| = |C_n| = (n-1)! 2^(n-1), no certificate failure.
 *  Build: cc -O2 -o lt_enum lt_enum.c ; run: ./lt_enum 9
 */
#include <stdio.h>
#include <stdlib.h>

static int n;
static int b[16][16];
static long long cntB = 0, failLT = 0, failReal = 0, cntC = 0;
static int PI[128], PL[128], NP;

static int pf(int i, int j, int k, int l) { return b[i][j]*b[k][l] - b[i][k]*b[j][l] + b[i][l]*b[j][k]; }

static int lt_masks(int nn, const unsigned *out) {
    unsigned all = (1u << nn) - 1;
    for (int v = 0; v < nn; v++) {
        for (int side = 0; side < 2; side++) {
            unsigned S = side ? (all & ~out[v] & ~(1u << v)) : out[v];
            unsigned seen = 0, SS = S;
            while (SS) {
                int u = __builtin_ctz(SS); SS &= SS - 1;
                int sc = __builtin_popcount(out[u] & S);
                if (seen & (1u << sc)) return 0;
                seen |= 1u << sc;
            }
        }
    }
    return 1;
}

static int realisable_certificate(void) {
    int eps[16], score[16]; unsigned seen = 0;
    eps[0] = 1;
    for (int j = 1; j < n; j++) eps[j] = b[0][j];
    for (int i = 0; i < n; i++) {
        int s = 0;
        for (int j = 0; j < n; j++) if (j != i && eps[i]*eps[j]*b[i][j] == 1) s++;
        score[i] = s;
        if (seen & (1u << s)) return 0;
        seen |= 1u << s;
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) {
        long long x1 = eps[i], y1 = (long long)eps[i]*(n - 1 - score[i]);
        long long x2 = eps[j], y2 = (long long)eps[j]*(n - 1 - score[j]);
        long long d = x1*y2 - y1*x2;
        int sg = (d > 0) - (d < 0);
        if (sg != b[i][j]) return 0;
    }
    return 1;
}

static void leafB(void) {
    unsigned out[16];
    cntB++;
    for (int i = 0; i < n; i++) { out[i] = 0; for (int j = 0; j < n; j++) if (j != i && b[i][j] == 1) out[i] |= 1u << j; }
    if (!lt_masks(n, out)) failLT++;
    if (!realisable_certificate()) failReal++;
}

static void recB(int idx) {
    if (idx == NP) { leafB(); return; }
    int i = PI[idx], l = PL[idx];
    for (int s = -1; s <= 1; s += 2) {
        b[i][l] = s; b[l][i] = -s;
        int ok = 1;
        for (int a = 0; a < i && ok; a++)
            for (int c = a + 1; c < i; c++)
                if (abs(pf(a, c, i, l)) == 3) { ok = 0; break; }
        if (ok) recB(idx + 1);
    }
    b[i][l] = b[l][i] = 0;
}

static void recC(int nn, const unsigned *out) {
    if (nn == n) { cntC++; return; }
    for (unsigned m = 0; m < (1u << nn); m++) {          /* m = old vertices beaten by the new vertex nn */
        unsigned o2[16];
        for (int v = 0; v < nn; v++) o2[v] = out[v] | (((m >> v) & 1) ? 0u : (1u << nn));
        o2[nn] = m;
        if (lt_masks(nn + 1, o2)) recC(nn + 1, o2);
    }
}

int main(int argc, char **argv) {
    int nmax = argc > 1 ? atoi(argv[1]) : 9;
    for (n = 4; n <= nmax; n++) {
        NP = 0;
        for (int l = 1; l < n; l++) for (int i = 0; i < l; i++) { PI[NP] = i; PL[NP] = l; NP++; }
        for (int i = 0; i < 16; i++) for (int j = 0; j < 16; j++) b[i][j] = 0;
        cntB = failLT = failReal = cntC = 0;
        recB(0);
        unsigned out0[16] = {0};
        recC(1, out0);
        long long expect = 1; for (int k = 2; k <= n - 1; k++) expect *= k; expect <<= (n - 1);
        printf("n=%d: |B| (all 4-Pfaffians +-1) = %lld, |C| (LT, independent enumeration) = %lld, "
               "(n-1)! 2^(n-1) = %lld; B-members failing LT: %lld, failing integer realisation: %lld -> %s\n",
               n, cntB, cntC, expect, failLT, failReal,
               (cntB == cntC && cntB == expect && failLT == 0 && failReal == 0) ? "OK" : "FAIL");
        fflush(stdout);
    }
    return 0;
}
