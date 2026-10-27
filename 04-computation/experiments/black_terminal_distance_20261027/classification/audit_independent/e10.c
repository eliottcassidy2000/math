/* e10.c -- check the explicit claims of PROOF.md Section 1 about the order-10 exception E10:
 * scores, lambda, delta, unique minimum cut, 9 candidates all blue and closed, |Aut| = 9, 9 anti-automorphisms of which
 * 3 involutions, 657 Ham paths of which 558 black, labelled black BFS: level 1 = 558 moves (targets by lambda), level 2
 * reaches how many distinct labelled non-strong tournaments.
 */
#include "at_nauty.h"
#include "cand.h"

static const tour *gT; static long long nh = 0, nblack = 0;
static int Qs[20000][MAXV]; static int nq = 0;
static int cb(const int *Q, void *c) { (void)c; nh++; if (!t_blue(gT, Q)) { nblack++; memcpy(Qs[nq++], Q, sizeof(int) * MAXV); } return 0; }
static int antis = 0, invol = 0;
static void permrec(const tour *T, int *p, uint32_t used, int k) {
    int n = T->n;
    if (k == n) { if (t_is_anti(T, p)) { antis++; int inv = 1; for (int v = 0; v < n; v++) if (p[p[v]] != v) inv = 0; invol += inv; } return; }
    for (int c = 0; c < n; c++) if (!(used >> c & 1u) && T->s[c] == n - 1 - T->s[k]) {
        int ok = 1; for (int j = 0; j < k && ok; j++) if (t_has(T, k, j) != t_has(T, p[j], c)) ok = 0;
        if (ok) { p[k] = c; permrec(T, p, used | (1u << c), k + 1); }
    }
}
/* hash set of labelled tournaments (45-bit codes) */
#define HS (1 << 22)
static uint64_t *hs;
static int hs_insert(uint64_t key) { uint64_t h = (key * 0x9E3779B97F4A7C15ull) >> 42; while (hs[h]) { if (hs[h] == key + 1) return 0; h = (h + 1) & (HS - 1); } hs[h] = key + 1; return 1; }
static uint64_t lcode(const tour *T) { uint64_t k = 0; int b = 0; for (int i = 0; i < T->n; i++) for (int j = i + 1; j < T->n; j++) { if (t_has(T, i, j)) k |= 1ull << b; b++; } return k; }
static const tour *gK; static long long lvl2moves = 0, lvl2nonstrong = 0;
static int cb2(const int *Q, void *c) { (void)c; if (t_blue(gK, Q)) return 0; tour K2; t_move(gK, Q, &K2); lvl2moves++; if (!t_strong_mask(&K2, fullmask(K2.n))) { lvl2nonstrong++; hs_insert(lcode(&K2)); } return 0; }

int main(void) {
    nauty_check(WORDSIZE, 1, 10, NAUTYVERSIONID);
    tour T; t_parse("111001111101111111111111111111111111000111101", 10, &T); gT = &T;
    printf("scores:"); for (int v = 0; v < 10; v++) printf(" %d", T.s[v]); printf("\n");
    static uint32_t cuts[1024]; int nc; int lam = t_mincuts(&T, cuts, 1024, &nc);
    printf("lambda=%d delta=%d ncuts=%d cut0=", lam, t_delta(&T), nc); for (int v = 0; v < 10; v++) if (cuts[0] >> v & 1u) printf("%d", v); printf("\n");
    printf("|Aut| own=%lld nauty=%g\n", t_auts(&T, NULL, NULL, NULL), nt_autsize(&T));
    int p[MAXV]; permrec(&T, p, 0, 0); printf("anti-automorphisms=%d involutions=%d\n", antis, invol);
    t_hpaths(&T, cb, NULL); printf("Ham paths=%lld black=%lld\n", nh, nblack);
    long long tgt[8] = {0}; hs = (uint64_t *)calloc(HS, 8);
    for (int i = 0; i < nq; i++) { tour K; t_move(&T, Qs[i], &K); int ncK; static uint32_t ck[1024]; int lk = t_mincuts(&K, ck, 1024, &ncK); tgt[lk]++;
        gK = &K; t_hpaths(&K, cb2, NULL); }
    printf("level-1 targets by lambda: 0:%lld 1:%lld 2:%lld 3:%lld\n", tgt[0], tgt[1], tgt[2], tgt[3]);
    long long distinct = 0; for (long long h = 0; h < HS; h++) if (hs[h]) distinct++;
    printf("level-2: black moves from level-1 targets=%lld, landing non-strong=%lld, distinct labelled non-strong tournaments=%lld\n", lvl2moves, lvl2nonstrong, distinct);
    return 0;
}
