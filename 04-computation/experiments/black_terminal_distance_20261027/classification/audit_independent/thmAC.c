/* thmAC.c -- auditor's independent check of Theorem A and Theorem C (PROOF.md) on gentourng -c output.
 *
 * Usage:  gentourng -cq n | ./thmAC n all     (every strong class; n <= 9)
 *         gentourng -cq n | ./thmAC n sc      (only nonregular self-converse classes, self-converseness tested
 *                                               by nauty canonical forms of T and its converse; the others are
 *                                               only tallied for lambda / delta statistics)
 * For every examined class: all minimum out-cuts, all Hamiltonian paths (plain DFS), and for each path whether it is
 * a descent candidate, whether it is black, whether it is of singleton type.  Reports
 *   - classes with NO black descent candidate (Theorem A exceptions), with codes;
 *   - classes with lambda = delta and no singleton-type black descent (Theorem C failures);
 *   - unconditional controls: minimum cuts never cross (III.2); lambda = delta for n <= 9 (III.1);
 *     F10 (self-converse or not: a min-score vertex v with T - v non-strong => every Ham-cycle cut from v is black);
 *     F6 (H = |Aut| only for C3) via exact path counts and nauty group orders.
 */
#include "at_nauty.h"

typedef struct {
    const tour *T; uint32_t sep[MAXV]; int lambda;
    int foundBlackCand, foundSingleton;
    long long nPaths, nCand, nCandBlue;
    uint32_t f10start; long long f10viol, f10checked;
    int needAll;
} ctx_t;

static int cb(const int *Q, void *vc) {
    ctx_t *c = (ctx_t *)vc; const tour *T = c->T; int n = T->n;
    c->nPaths++;
    int a = Q[0], b = Q[n - 1];
    int closed = t_has(T, b, a);
    int blueKnown = 0, blue = 0;
    if ((c->f10start >> a & 1u) && closed) {
        blue = t_blue(T, Q); blueKnown = 1; c->f10checked++;
        if (blue) c->f10viol++;
    }
    if (c->sep[a] >> b & 1u) {
        c->nCand++;
        if (!blueKnown) { blue = t_blue(T, Q); blueKnown = 1; }
        if (blue) c->nCandBlue++;
        else {
            c->foundBlackCand = 1;
            int single = (T->s[a] == c->lambda) || (n - 1 - T->s[b] == c->lambda);
            if (single) c->foundSingleton = 1;
        }
    }
    if (!c->needAll && c->foundBlackCand && c->foundSingleton) return 1;
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: thmAC n all|sc\n"); return 2; }
    int n = atoi(argv[1]); int scMode = strcmp(argv[2], "sc") == 0;
    int fullScan = (argc > 3 && strcmp(argv[3], "full") == 0);
    nauty_check(WORDSIZE, 1, n, NAUTYVERSIONID);
    char line[1024]; tour T;
    long long classes = 0, examined = 0, lamHist[16] = {0}, lamLtDelta = 0, cross = 0, excA = 0, failC = 0, lamEqDelta = 0;
    long long f10viol = 0, f10checked = 0, f6hits = 0, scCount = 0, regCount = 0, scDisagree = 0;
    uint32_t cuts[4096];
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        if (!t_parse(line, n, &T)) continue;
        if (!t_valid(&T) || !t_strong_mask(&T, fullmask(n))) { fprintf(stderr, "bad/nonstrong input\n"); return 1; }
        classes++;
        int nc; int lam = t_mincuts(&T, cuts, 4096, &nc);
        if (nc > 4096) { fprintf(stderr, "too many cuts\n"); return 1; }
        int del = t_delta(&T);
        lamHist[lam]++;
        if (lam < del) lamLtDelta++;
        if (lam > del) { printf("ERROR lambda > delta\n"); return 1; }
        int reg = t_regular(&T);
        if (reg) regCount++;
        if (scMode) {
            if (reg) continue;
            int ss[MAXV]; for (int v = 0; v < n; v++) ss[v] = T.s[v];
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (ss[j] < ss[i]) { int x = ss[i]; ss[i] = ss[j]; ss[j] = x; }
            int sym = 1; for (int i = 0; i < n; i++) if (ss[i] + ss[n - 1 - i] != n - 1) sym = 0;
            if (!sym) continue;
            int scN = nt_selfconverse(&T); int sgx[MAXV]; int scB = t_find_anti(&T, sgx);
            if (scB && !t_is_anti(&T, sgx)) { printf("ERROR anti search\n"); return 1; }
            if (scN != scB) { scDisagree++; t_code(&T, line); printf("SC-DISAGREE %s nauty=%d backtrack=%d\n", line, scN, scB); }
            if (!scN) continue;
            scCount++;
        }
        examined++;
        /* III.2: no two minimum cuts cross (all four Venn cells nonempty) */
        uint32_t F = fullmask(n);
        for (int i = 0; i < nc; i++) for (int j = i + 1; j < nc; j++) {
            uint32_t X = cuts[i], Z = cuts[j];
            if ((X & Z) && (X & ~Z) && (Z & ~X) && (F & ~(X | Z))) cross++;
        }
        ctx_t c; memset(&c, 0, sizeof c); c.T = &T; c.lambda = lam; c.needAll = fullScan;
        for (int i = 0; i < nc; i++) for (int a = 0; a < n; a++) if (cuts[i] >> a & 1u) c.sep[a] |= F & ~cuts[i];
        /* F10 starts: min-score vertices whose deletion is non-strong */
        int smin = 99; for (int v = 0; v < n; v++) if (T.s[v] < smin) smin = T.s[v];
        for (int v = 0; v < n; v++) if (T.s[v] == smin && !t_strong_mask(&T, F & ~(1u << v))) c.f10start |= 1u << v;
        if (c.f10start) c.needAll = 1;   /* scan all paths so that F10 is fully tested */
        t_hpaths(&T, cb, &c);
        f10viol += c.f10viol; f10checked += c.f10checked;
        char code[512];
        if (!c.foundBlackCand) {
            excA++; t_code(&T, code);
            printf("NO-BLACK-DESCENT %s lambda=%d delta=%d ncuts=%d regular=%d paths=%lld cand=%lld candBlue=%lld\n",
                   code, lam, del, nc, reg, c.nPaths, c.nCand, c.nCandBlue);
        }
        if (lam == del) {
            lamEqDelta++;
            if (!c.foundSingleton) {
                failC++; t_code(&T, code);
                printf("NO-SINGLETON-BLACK-DESCENT %s lambda=%d delta=%d regular=%d\n", code, lam, del, reg);
            }
        }
        /* F6: H = |Aut| */
        if (!scMode || 1) {
            double aut = nt_autsize(&T);
            if (aut > 1.5) {
                long long H = t_count_hp_end(&T, F, F);
                if ((double)H == aut) { f6hits++; t_code(&T, code); printf("F6-HIT H=|Aut|=%lld %s\n", H, code); }
            }
        }
    }
    printf("n=%d mode=%s strong classes read=%lld examined=%lld (self-converse nonregular=%lld, regular=%lld)\n",
           n, argv[2], classes, examined, scCount, regCount);
    printf("  lambda histogram over all classes read:");
    for (int l = 0; l < 16; l++) if (lamHist[l]) printf(" %d:%lld", l, lamHist[l]);
    printf("\n  classes with lambda < delta (all read): %lld\n", lamLtDelta);
    printf("  examined: crossing min-cut pairs=%lld  no-black-descent=%lld  lambda=delta=%lld  TheoremC-failures=%lld\n",
           cross, excA, lamEqDelta, failC);
    printf("  F10 control: Ham-cycle cuts from min-score v with T-v nonstrong: checked=%lld blue(violations)=%lld\n", f10checked, f10viol);
    printf("  F6 control: classes (examined) with H=|Aut|>1: %lld\n", f6hits);
    printf("  self-converse test disagreements (nauty vs own backtracking, sc mode only): %lld\n", scDisagree);
    return 0;
}
