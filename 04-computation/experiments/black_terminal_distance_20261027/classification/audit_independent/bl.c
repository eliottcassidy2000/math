/* bl.c -- auditor's independent check of Lemma BL and of the bottom-layer census (PROOF.md Sections 3, 9).
 * For every strong S on a vertices (gentourng -cq a) and every z in S:
 *   trans  : Stab_{Aut S}(z) acts transitively on Ham(S -> z)   (tested directly: orbit of one path = all paths)
 *   condA  : every nonempty W subset of S - z has e_S(W -> S\W) >= 2                     (Condition A)
 *   H2     : with R_q the last strong component of S - z and B_q = R_q cap N^-(z):
 *            every nonempty W subset of R_q has e_{R_q}(W -> R_q\W) + |W cap B_q| >= 2, and |R_q| >= 2
 *   BLform : S - z = R_1 => ... => R_q, q >= 2, every block K1 or C3, z uniform on every block, z => R_1, R_q = C3, R_q => z
 * Lemma BL claims trans & H2 => BLform.  The bottom-layer census claims #{trans & condA} = b(a) = 0,0,1,2,5,12,28 (a=3..9),
 * all of BLform.  Pairs are counted per canonical class representative and per vertex z (z in one Aut(S)-orbit are
 * counted separately; we also report the number of Aut(S)-orbits).
 * Usage: gentourng -cq a | ./bl a
 */
#include "at_nauty.h"

static int stabperm[2000][MAXV]; static int nstab; static int gz;
static void stab_cb(const int *p, void *ctx) { (void)ctx; if (p[gz] == gz && nstab < 2000) memcpy(stabperm[nstab++], p, sizeof(int) * MAXV); }
static int pathsZ[200000][MAXV]; static int npz; static int gzEnd; static const tour *gS;
static int pz_cb(const int *Q, void *ctx) { (void)ctx; if (Q[gS->n - 1] == gzEnd) { memcpy(pathsZ[npz++], Q, sizeof(int) * gS->n); } return 0; }
static int cmpP(const void *x, const void *y) { return memcmp(x, y, sizeof(int) * MAXV); }

/* strong components of T[mask] in condensation order (first = source component) */
static int condense(const tour *T, uint32_t mask, uint32_t *comp) {
    int k = 0; uint32_t left = mask;
    while (left) {
        /* source component: the component of a vertex that reaches all of left */
        uint32_t found = 0;
        for (uint32_t c = left; c && !found; c &= c - 1) {
            int v = lowbit(c);
            uint32_t R = 1u << v, fr = R;
            while (fr) { int x = lowbit(fr); fr &= fr - 1; uint32_t nb = T->o[x] & left & ~R; R |= nb; fr |= nb; }
            if (R == left) {
                /* vertices of left that reach v */
                uint32_t Rb = 1u << v; fr = Rb;
                while (fr) { int x = lowbit(fr); fr &= fr - 1; uint32_t nb = t_in(T, x) & left & ~Rb; Rb |= nb; fr |= nb; }
                found = Rb;
            }
        }
        comp[k++] = found; left &= ~found;
    }
    return k;
}

int main(int argc, char **argv) {
    int a = atoi(argv[1]);
    nauty_check(WORDSIZE, 1, a, NAUTYVERSIONID);
    char line[512]; tour S; uint32_t F = fullmask(a);
    long long pairs = 0, ntrans = 0, ntransH2 = 0, ntransH2BL = 0, ntransCondA = 0, ntransCondABL = 0, transBL = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        t_parse(line, a, &S); gS = &S;
        for (int z = 0; z < a; z++) {
            pairs++;
            /* Stab(z) by enumeration with colour z distinct */
            int col[MAXV]; for (int v = 0; v < a; v++) col[v] = (v == z);
            nstab = 0; gz = z; t_auts(&S, col, stab_cb, NULL);
            npz = 0; gzEnd = z; t_hpaths(&S, pz_cb, NULL);
            if (npz == 0) { printf("no path ending at z?\n"); return 1; }
            /* orbit of pathsZ[0] */
            static int orb[2000][MAXV]; int no = 0;
            for (int k = 0; k < nstab; k++) { memset(orb[no], 0, sizeof orb[no]); for (int i = 0; i < a; i++) orb[no][i] = stabperm[k][pathsZ[0][i]]; no++; }
            qsort(orb, no, sizeof orb[0], cmpP);
            int distinct = 0; for (int k = 0; k < no; k++) if (k == 0 || memcmp(orb[k], orb[k - 1], sizeof orb[0])) distinct++;
            int trans = (distinct == npz);
            uint32_t rest = F & ~(1u << z);
            int condA = 1; for (uint32_t W = rest; W && condA; W = (W - 1) & rest) if (t_f(&S, W) < 2) condA = 0;
            uint32_t comp[MAXV]; int q = condense(&S, rest, comp);
            uint32_t Rq = comp[q - 1]; uint32_t Bq = Rq & ~S.o[z];
            int H2 = pc(Rq) >= 2;
            for (uint32_t W = Rq; W && H2; W = (W - 1) & Rq) {
                int e = 0; for (uint32_t c = W; c; c &= c - 1) e += pc(S.o[lowbit(c)] & Rq & ~W);
                if (e + pc(W & Bq) < 2) H2 = 0;
            }
            int BL = (q >= 2);
            for (int i = 0; i < q && BL; i++) {
                int sz = pc(comp[i]);
                if (sz != 1 && sz != 3) BL = 0;
                uint32_t A = comp[i] & S.o[z];
                if (A != 0 && A != comp[i]) BL = 0;
            }
            if (BL) { if ((comp[0] & S.o[z]) != comp[0]) BL = 0; if (pc(Rq) != 3 || (Rq & S.o[z])) BL = 0; }
            if (trans) ntrans++;
            if (trans && BL) transBL++;
            if (trans && H2) { ntransH2++; if (BL) ntransH2BL++; else { char code[256]; t_code(&S, code); printf("LEMMA-BL-FAIL S=%s z=%d\n", code, z); } }
            if (trans && condA) {
                ntransCondA++; if (BL) ntransCondABL++;
                char code[256]; t_code(&S, code);
                printf("BOTTOM S=%s z=%d |Stab|=%d h=%d q=%d BLform=%d blocks:", code, z, nstab, npz, q, BL);
                for (int i = 0; i < q; i++) printf(" %d%c", pc(comp[i]), (comp[i] & S.o[z]) == comp[i] ? '+' : ((comp[i] & S.o[z]) ? '?' : '-'));
                printf("\n");
            }
        }
    }
    printf("a=%d pairs(S,z)=%lld transitive=%lld (of which BL-form %lld)  transitive&H2=%lld (BL-form %lld)  transitive&CondA=%lld (BL-form %lld)\n",
           a, pairs, ntrans, transBL, ntransH2, ntransH2BL, ntransCondA, ntransCondABL);
    return 0;
}
