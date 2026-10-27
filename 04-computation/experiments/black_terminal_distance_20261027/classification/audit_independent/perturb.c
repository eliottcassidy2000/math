/* perturb.c -- try to break Theorem A near the family E: sigma-symmetric arc flips of members of E (n <= nmax).
 * A flip of a sigma-orbit of vertex pairs keeps sigma an anti-automorphism, so the result stays self-converse.
 * Single-orbit flips, all double-orbit flips (n <= dmax), and random multi-orbit flips; every strong result without a
 * black descent candidate is compared (nauty canonical form) with the auditor's own list of E members and R5.
 * Usage: ./perturb nmax dmax randomSamplesPerMember seed
 */
#include "at_nauty.h"
#include "cand.h"
#include "fam.h"

static char Ecanon[400][400]; static int nE = 0;
static int inE(const char *cs) { for (int i = 0; i < nE; i++) if (!strcmp(Ecanon[i], cs)) return 1; return 0; }
static long long tested = 0, strongc = 0, exc = 0, excE = 0, excNew = 0, regimeLt = 0;
static void test(const tour *T, const char *tag) {
    tested++;
    if (!t_strong_mask(T, fullmask(T->n))) return;
    strongc++;
    int Q[MAXV], lam; long long nc, nb;
    int found = black_descent(T, 0, Q, &lam, &nc, &nb);
    if (lam < t_delta(T)) regimeLt++;
    if (!found) {
        exc++;
        char cs[400]; nt_canon_str(T, cs);
        if (inE(cs)) excE++;
        else { excNew++; char code[512]; t_code(T, code); printf("NEW-EXCEPTION %s n=%d lambda=%d delta=%d %s\n", tag, T->n, lam, t_delta(T), code); }
    }
}
static uint64_t rs = 88172645463325252ull;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

int main(int argc, char **argv) {
    int nmax = atoi(argv[1]), dmax = atoi(argv[2]); int R = atoi(argv[3]); rs ^= (uint64_t)atoll(argv[4]) * 2654435761ull;
    nauty_check(WORDSIZE, 1, nmax, NAUTYVERSIONID);
    int parts[16];
    for (int a = 5; a <= nmax / 2; a++) gen_comp(a - 4, parts, 0, a);
    /* E list up to nmax */
    for (int bi = 0; bi < nBL; bi++) {
        tour S; int bs[17]; build_S(&BL[bi], &S, bs);
        for (int t = 1; 2 * S.n + t - 1 <= nmax; t++) { tour T; int sg[MAXV]; build_E(&S, t, &T, sg); nt_canon_str(&T, Ecanon[nE++]); }
    }
    for (int bi = 0; bi < nBL; bi++) {
        tour S; int bs[17]; build_S(&BL[bi], &S, bs);
        for (int t = 1; 2 * S.n + t - 1 <= nmax; t++) {
            tour T; int sg[MAXV]; build_E(&S, t, &T, sg); int n = T.n;
            /* sigma-orbits of unordered pairs */
            int ou[400], ov[400], no = 0; static int seen[MAXV][MAXV]; memset(seen, 0, sizeof seen);
            for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) {
                if (seen[u][v]) continue;
                seen[u][v] = seen[v][u] = 1; seen[sg[u]][sg[v]] = seen[sg[v]][sg[u]] = 1;
                ou[no] = u; ov[no] = v; no++;
            }
            #define FLIP(TT, k) do { int _u = ou[k], _v = ov[k]; int _uv = t_has(&TT, _u, _v); \
                if (_uv) t_arc(&TT, _v, _u); else t_arc(&TT, _u, _v); \
                int _a = sg[_u], _b = sg[_v]; if (!((_a == _u && _b == _v) || (_a == _v && _b == _u))) { \
                    int _ab = t_has(&TT, _a, _b); if (_ab) t_arc(&TT, _b, _a); else t_arc(&TT, _a, _b); } } while (0)
            for (int k1 = 0; k1 < no; k1++) {
                tour U = T; FLIP(U, k1); t_fin(&U);
                if (!t_is_anti(&U, sg)) { printf("flip broke sigma\n"); return 1; }
                test(&U, "single");
                if (n <= dmax) for (int k2 = k1 + 1; k2 < no; k2++) { tour W = U; FLIP(W, k2); t_fin(&W); test(&W, "double"); }
            }
            for (int r = 0; r < R; r++) {
                tour U = T; int nf = 2 + (int)(rnd() % 9);
                for (int f = 0; f < nf; f++) { int k = (int)(rnd() % no); FLIP(U, k); }
                t_fin(&U); test(&U, "random");
            }
        }
    }
    printf("nmax=%d: E members=%d tested=%lld strong=%lld (lambda<delta among strong %lld) no-black-descent=%lld (in E: %lld, NEW: %lld)\n",
           nmax, nE, tested, strongc, regimeLt, exc, excE, excNew);
    return 0;
}
