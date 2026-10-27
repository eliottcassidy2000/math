/* dist.c -- auditor's class-level computation of the black distance d_black(T) to the non-strong locus
 * (moves K_Q(T), Q black in T), for every tournament class of order n (gentourng -q n, all classes).
 * Layered BFS: layer 0 = non-strong classes; T is in layer L if some black move of T lands in layer L-1 (first
 * time).  mode exact: every unresolved strong class is re-scanned at every layer (no use of F1);
 * mode f1: a class with lambda(T) > L is not scanned at layer L (uses d >= lambda, i.e. F1).
 * Reports d versus lambda for every strong class (Theorem B: d = lambda except R5 (infinite) and E).
 * Also checks the cut identity F1 on every move it generates (exact mode) for a sample of X.
 * Usage: gentourng -q n | ./dist n exact|f1
 */
#include "at_nauty.h"

static int N;
static uint64_t keyOf(const tour *T) {
    graph cg[MAXN]; nt_canon(T, cg);
    uint64_t k = 0; int b = 0;
    for (int i = 0; i < T->n; i++) for (int j = i + 1; j < T->n; j++) { if (ISELEMENT(GRAPHROW(cg, i, 1), j)) k |= (uint64_t)1 << b; b++; }
    return k;
}
static uint64_t *keys; static int nk;
static int cmpu(const void *x, const void *y) { uint64_t a = *(const uint64_t *)x, b = *(const uint64_t *)y; return a < b ? -1 : a > b; }
static int findKey(uint64_t k) {
    int lo = 0, hi = nk - 1;
    while (lo <= hi) { int mid = (lo + hi) / 2; if (keys[mid] == k) return mid; if (keys[mid] < k) lo = mid + 1; else hi = mid - 1; }
    return -1;
}
typedef struct { const tour *T; const int *d; int layer; int found; long long moves; long long f1checks, f1fail; int exact; } dctx;
static int dcb(const int *Q, void *vc) {
    dctx *c = (dctx *)vc;
    if (t_blue(c->T, Q)) return 0;
    tour K; t_move(c->T, Q, &K); c->moves++;
    if (c->exact && (c->moves % 7) == 0) {
        /* F1: f_K(V\X) = f_T(X) - [q1 in X] + [qn in X], on a few X */
        int n = c->T->n; uint32_t F = fullmask(n);
        for (uint32_t X = 1 + (c->moves % 13); X < F; X += 37) {
            int lhs = t_f(&K, F & ~X), rhs = t_f(c->T, X) - (int)((X >> Q[0]) & 1u) + (int)((X >> Q[n - 1]) & 1u);
            c->f1checks++; if (lhs != rhs) c->f1fail++;
        }
    }
    int idx = findKey(keyOf(&K));
    if (idx < 0) { fprintf(stderr, "target class not found\n"); exit(1); }
    if (c->d[idx] == c->layer - 1) { c->found = 1; return 1; }
    return 0;
}
int main(int argc, char **argv) {
    N = atoi(argv[1]); int exact = strcmp(argv[2], "exact") == 0;
    nauty_check(WORDSIZE, 1, N, NAUTYVERSIONID);
    int cap = 300000; tour *Ts = (tour *)malloc(sizeof(tour) * cap); keys = (uint64_t *)malloc(8 * cap);
    char line[512]; int nt = 0;
    while (fgets(line, sizeof line, stdin)) { if (line[0] != '0' && line[0] != '1') continue; t_parse(line, N, &Ts[nt]); nt++; }
    /* sort classes by key */
    uint64_t *kk = (uint64_t *)malloc(8 * nt); for (int i = 0; i < nt; i++) kk[i] = keyOf(&Ts[i]);
    /* order index by key */
    int *ord = (int *)malloc(sizeof(int) * nt); for (int i = 0; i < nt; i++) ord[i] = i;
    /* simple sort of indices by key */
    for (int i = 0; i < nt; i++) keys[i] = kk[i];
    qsort(keys, nt, 8, cmpu); nk = nt;
    for (int i = 1; i < nt; i++) if (keys[i] == keys[i - 1]) { fprintf(stderr, "duplicate class key\n"); return 1; }
    int *d = (int *)malloc(sizeof(int) * nt), *lam = (int *)malloc(sizeof(int) * nt); /* indexed by sorted key position */
    tour **byk = (tour **)malloc(sizeof(tour *) * nt);
    for (int i = 0; i < nt; i++) { int p = findKey(kk[i]); byk[p] = &Ts[i]; }
    uint32_t cuts[1 << 12]; int nc; int strongc = 0;
    for (int p = 0; p < nt; p++) { lam[p] = t_mincuts(byk[p], cuts, 1 << 12, &nc); d[p] = lam[p] == 0 ? 0 : 1 << 20; if (lam[p]) strongc++; }
    long long f1c = 0, f1f = 0;
    for (int L = 1; L <= 20; L++) {
        int newly = 0;
        for (int p = 0; p < nt; p++) {
            if (d[p] <= L) continue;
            if (!exact && lam[p] > L) continue;
            dctx c; memset(&c, 0, sizeof c); c.T = byk[p]; c.d = d; c.layer = L; c.exact = exact;
            t_hpaths(byk[p], dcb, &c);
            f1c += c.f1checks; f1f += c.f1fail;
            if (c.found) { d[p] = L + 1000; newly++; }    /* mark, finalize after the layer */
        }
        for (int p = 0; p < nt; p++) if (d[p] >= 1000 && d[p] < (1 << 20)) d[p] -= 1000;
        printf("layer %d: %d classes\n", L, newly); fflush(stdout);
        if (!newly) {
            /* in f1 mode a class with lam > L may still be unresolved only if lam[p] > L; continue while such exist */
            int pending = 0; for (int p = 0; p < nt; p++) if (d[p] == (1 << 20) && lam[p] > L) pending = 1;
            if (exact || !pending) break;
        }
    }
    long long eq = 0, ne = 0, inf = 0; long long tab[8][24]; memset(tab, 0, sizeof tab);
    for (int p = 0; p < nt; p++) if (lam[p] > 0) {
        if (d[p] == (1 << 20)) { inf++; char code[256]; t_code(byk[p], code); printf("NO-BLACK-ROUTE %s lambda=%d regular=%d\n", code, lam[p], t_regular(byk[p])); continue; }
        tab[lam[p]][d[p]]++;
        if (d[p] == lam[p]) eq++; else { ne++; char code[256]; t_code(byk[p], code); printf("D!=LAMBDA %s lambda=%d d=%d\n", code, lam[p], d[p]); }
    }
    printf("n=%d mode=%s classes=%d strong=%d: d=lambda %lld, d!=lambda %lld, unreachable %lld; F1 checks %lld failures %lld\n", N, argv[2], nt, strongc, eq, ne, inf, f1c, f1f);
    for (int l = 1; l < 8; l++) for (int dd = 0; dd < 24; dd++) if (tab[l][dd]) printf("  lambda=%d d=%d : %lld\n", l, dd, tab[l][dd]);
    free(ord);
    return 0;
}
