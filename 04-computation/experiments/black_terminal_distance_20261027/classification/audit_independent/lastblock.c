/* lastblock.c -- auditor's independent finite check of the last-block lemma (PROOF.md Part VI).
 * For every strong R on m vertices (gentourng -cq m) and every A subset V(R) with B = V \ A nonempty:
 *   (H2) every nonempty W subset of V(R): e_R(W -> V\W) + |W cap B| >= 2   (W = V(R) included: |B| >= 2)
 *   (H1) G = Aut(R; A) (automorphisms with g(A) = A) acts transitively on H' = Ham paths of R ending in B.
 * (H1) is tested DIRECTLY: G is enumerated by the auditor's own backtracking, the orbit of one path of H' is
 * computed as a set and compared with H' (no appeal to freeness).  |G| is cross-checked against nauty.
 * Usage: gentourng -cq m | ./lastblock m [skiprigid]
 */
#include "at_nauty.h"

#define MAXP 400000
static int paths[MAXP][10]; static int npaths;
static uint32_t Bset; static const tour *gT;
static int hp_collect(const int *Q, void *ctx) {
    (void)ctx; int n = gT->n;
    if (!(Bset >> Q[n - 1] & 1u)) return 0;
    if (npaths >= MAXP) { fprintf(stderr, "MAXP\n"); exit(1); }
    memcpy(paths[npaths++], Q, sizeof(int) * n); return 0;
}
static int perms[100000][10]; static int nperm;
static void aut_collect(const int *p, void *ctx) { (void)ctx; if (nperm < 100000) memcpy(perms[nperm], p, sizeof(int) * gT->n); nperm++; }

static int cmpPath(const void *x, const void *y) { return memcmp(x, y, sizeof(int) * 10); }

int main(int argc, char **argv) {
    int m = atoi(argv[1]); int skiprigid = argc > 2;
    nauty_check(WORDSIZE, 1, m, NAUTYVERSIONID);
    char line[512]; tour R; uint32_t F = fullmask(m);
    long long pairs = 0, h2 = 0, both = 0, autDisagree = 0, rigidSkipped = 0, classes = 0;
    static int fW[1 << 10];
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        t_parse(line, m, &R); gT = &R; classes++;
        long long autR = t_auts(&R, NULL, NULL, NULL);
        if (skiprigid && autR == 1) { rigidSkipped++; continue; }
        for (uint32_t W = 1; W <= F; W++) fW[W] = t_f(&R, W);
        /* all Ham paths once */
        for (uint32_t A = 0; A < F; A++) {
            uint32_t B = F & ~A; pairs++;
            int ok = 1;
            for (uint32_t W = 1; W <= F && ok; W++) if (fW[W] + pc(W & B) < 2) ok = 0;
            if (!ok) continue;
            h2++;
            /* G = Aut(R; A) */
            int col[MAXV]; for (int v = 0; v < m; v++) col[v] = (A >> v) & 1u;
            nperm = 0; long long g = t_auts(&R, col, aut_collect, NULL);
            double gn = nt_autsize_col(&R, col);
            if ((double)g != gn) autDisagree++;
            /* H' */
            Bset = B; npaths = 0; t_hpaths(&R, hp_collect, NULL);
            if (npaths == 0) continue;
            /* orbit of paths[0] under G, as a set */
            static int orbit[100000][10]; int no = 0;
            for (int k = 0; k < nperm && k < 100000; k++) { int img[10] = {0}; for (int i = 0; i < m; i++) img[i] = perms[k][paths[0][i]]; memcpy(orbit[no++], img, sizeof img); }
            qsort(orbit, no, sizeof orbit[0], cmpPath);
            int distinct = 0; for (int k = 0; k < no; k++) if (k == 0 || memcmp(orbit[k], orbit[k - 1], sizeof orbit[0])) distinct++;
            int transitive = (distinct == npaths);
            if (transitive) {
                both++;
                char code[64]; t_code(&R, code);
                printf("(H1)+(H2): R=%s A=%x |G|=%lld |H'|=%d %s\n", code, A, g, npaths, (A == 0 && m == 3) ? "(C3, A empty: the lemma's conclusion)" : "COUNTEREXAMPLE");
            }
        }
    }
    printf("m=%d classes=%lld rigid-skipped=%lld pairs=%lld (H2)=%lld (H1)+(H2)=%lld  |G| own-vs-nauty disagreements=%lld\n",
           m, classes, rigidSkipped, pairs, h2, both, autDisagree);
    return 0;
}
