/* epred.c -- what exactly does Theorem C's case (c2) prove?  For every nonregular strong SELF-CONVERSE class:
 *   C'  : some Hamiltonian path starting in m (min-score set) or ending in M (max-score set) is BLACK
 *         (Steps 1-3 of Theorem C never use lambda = delta; (c1) via F10 also gives such a path)
 *   E   : some Hamiltonian path starts in m and ends outside M   (the literal open-path form of the
 *         extremal-predecessor statement of cycles.md Section 7)
 *   Ecl : some CLOSED Hamiltonian path starts in m and ends outside M (the OPEN closed form, cycles.md Section 7)
 * Usage: gentourng -cq n | ./epred n
 */
#include "at_nauty.h"

typedef struct { const tour *T; uint32_t mset, Mset; int Cp, E, Ecl; } ectx;
static int ecb(const int *Q, void *vc) {
    ectx *c = (ectx *)vc; int n = c->T->n; int a = Q[0], b = Q[n - 1];
    int startm = c->mset >> a & 1u, endM = c->Mset >> b & 1u;
    if (startm && !endM) { c->E = 1; if (t_has(c->T, b, a)) c->Ecl = 1; }
    if ((startm || endM) && !c->Cp) { if (!t_blue(c->T, Q)) c->Cp = 1; }
    return c->Cp && c->E && c->Ecl;
}
int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    nauty_check(WORDSIZE, 1, n, NAUTYVERSIONID);
    char line[1024]; tour T; long long sc = 0, fCp = 0, fE = 0, fEcl = 0, lamLt = 0, fElamLt = 0;
    uint32_t cuts[4096]; int nc;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        t_parse(line, n, &T);
        if (t_regular(&T)) continue;
        int ss[MAXV]; for (int v = 0; v < n; v++) ss[v] = T.s[v];
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (ss[j] < ss[i]) { int x = ss[i]; ss[i] = ss[j]; ss[j] = x; }
        int sym = 1; for (int i = 0; i < n; i++) if (ss[i] + ss[n - 1 - i] != n - 1) sym = 0;
        if (!sym || !nt_selfconverse(&T)) continue;
        sc++;
        int smin = ss[0], smax = ss[n - 1];
        ectx c; memset(&c, 0, sizeof c); c.T = &T;
        for (int v = 0; v < n; v++) { if (T.s[v] == smin) c.mset |= 1u << v; if (T.s[v] == smax) c.Mset |= 1u << v; }
        t_hpaths(&T, ecb, &c);
        int lam = t_mincuts(&T, cuts, 4096, &nc); int del = t_delta(&T);
        if (lam < del) lamLt++;
        char code[512]; line[strcspn(line, "\r\n")] = 0; strcpy(code, line);
        if (!c.Cp) { fCp++; printf("C'-FAIL %s lambda=%d delta=%d\n", code, lam, del); }
        if (!c.E) { fE++; if (lam < del) fElamLt++; printf("E-FAIL %s lambda=%d delta=%d |m|=%d\n", code, lam, del, pc(c.mset)); }
        if (!c.Ecl) { fEcl++; printf("Ecl-FAIL %s lambda=%d delta=%d\n", code, lam, del); }
    }
    printf("n=%d nonregular strong self-converse=%lld (lambda<delta: %lld): C' fails=%lld  E fails=%lld (with lambda<delta %lld)  Ecl fails=%lld\n",
           n, sc, lamLt, fCp, fE, fElamLt, fEcl);
    return 0;
}
