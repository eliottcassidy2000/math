/* audit_bfs_onfly.c -- class-level native-move BFS without storing the class graph (for order 9).
 * Written from scratch for the audit.  Each class is expanded once per BFS: all Hamiltonian paths are enumerated,
 * the move tau_Q(T) = T^op xor E(Q) is canonised with nauty (DEFAULTOPTIONS_DIGRAPH, libnauty.a), and unvisited
 * neighbour classes are put on the next level.  While a class is expanded we also record whether it has a neighbour
 * whose potential is exactly one lower (local-minimum test) -- valid for every class, since every class is expanded.
 * Targets (argv[2..]): REG (odd n: regular, R), REGB (same, black moves only), TT (transitive, L), SING (source or
 * sink, delta0), NS (non-strong, lambda).
 * Reports: classes, unreachable, d < potential (must be 0), failures (d > potential) with symmetric count, local
 * minima, strong failures with kappa >= 2, max excess, and the distance histogram by (potential, excess).
 * Build: cc -O3 -march=native -I/opt/homebrew/opt/nauty/include audit_bfs_onfly.c /opt/homebrew/opt/nauty/lib/libnauty.a -o audit_bfs_onfly
 * Run:   gentourng -q 9 | ./audit_bfs_onfly 9 REG REGB TT SING NS
 */
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, m_words; static uint32_t full;
typedef struct { uint32_t o[12]; } Tour;
static uint64_t canon(const Tour *T, double *grp) {
    DYNALLSTAT(graph, g, g_sz); DYNALLSTAT(graph, cg, cg_sz);
    DYNALLSTAT(int, lab, lab_sz); DYNALLSTAT(int, ptn, ptn_sz); DYNALLSTAT(int, orbits, orbits_sz);
    static DEFAULTOPTIONS_DIGRAPH(options); statsblk stats; int m = m_words;
    DYNALLOC2(graph, g, g_sz, m, n, "c"); DYNALLOC2(graph, cg, cg_sz, m, n, "c");
    DYNALLOC1(int, lab, lab_sz, n, "c"); DYNALLOC1(int, ptn, ptn_sz, n, "c"); DYNALLOC1(int, orbits, orbits_sz, n, "c");
    options.getcanon = TRUE; EMPTYGRAPH(g, m, n);
    for (int u = 0; u < n; u++) for (uint32_t x = T->o[u]; x; x &= x - 1) ADDONEARC(g, u, __builtin_ctz(x), m);
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, cg);
    uint64_t code = 0; int b = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (ISELEMENT(GRAPHROW(cg, i, m), j)) code |= (uint64_t)1 << b; b++; }
    if (grp) { double s = stats.grpsize1; for (int e = 0; e < stats.grpsize2; e++) s *= 10.0; *grp = s; }
    return code;
}
static Tour converse(const Tour *T) { Tour R; memset(&R, 0, sizeof R); for (int u = 0; u < n; u++) R.o[u] = ~T->o[u] & full & ~(1u << u); return R; }
static int isBlue(const Tour *T, const int *Q) { int r[16]; for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) for (uint32_t x = T->o[u]; x; x &= x - 1) { int v = __builtin_ctz(x); if (!(T->o[r[v]] >> r[u] & 1)) return 0; } return 1; }
static uint64_t *hk; static int *hv; static size_t hcap;
static size_t hh(uint64_t k) { k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33; return (size_t)k; }
static void hput(uint64_t k, int v) { size_t i = hh(k) & (hcap - 1); while (hv[i] >= 0) { if (hk[i] == k) { hv[i] = v; return; } i = (i + 1) & (hcap - 1); } hk[i] = k; hv[i] = v; }
static int hget(uint64_t k) { size_t i = hh(k) & (hcap - 1); while (hv[i] >= 0) { if (hk[i] == k) return hv[i]; i = (i + 1) & (hcap - 1); } return -1; }
static uint32_t reach_in(const Tour *T, uint32_t vs, int v) { uint32_t seen = 1u << v, fr = seen; while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= T->o[__builtin_ctz(f)] & vs; nx &= ~seen; seen |= nx; fr = nx; } return seen; }
static int strong_in(const Tour *T, uint32_t vs) { int pc = __builtin_popcount(vs); if (pc <= 1) return 1; if (pc == 2) return 0; for (uint32_t w = vs; w; w &= w - 1) if (reach_in(T, vs, __builtin_ctz(w)) != vs) return 0; return 1; }
static int kappa_of(const Tour *T) { if (!strong_in(T, full)) return 0; for (int k = 1; k <= n - 2; k++) for (uint32_t S = 1; S < full; S++) if (__builtin_popcount(S) == k && !strong_in(T, full & ~S)) return k; return n - 1; }
static int lambda_of(const Tour *T) { int best = 1 << 20; for (uint32_t X = 1; X < full; X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T->o[x] & ~X & full); if (f < best) best = f; } return best; }
static int cmpint(const void *a, const void *b) { return *(const int *)a - *(const int *)b; }

static int C; static Tour *reps; static int8_t *dist; static int *pot, *tgt; static char *hasdesc;
static int *cur, *nxt; static int ncur, nnxt; static int blackOnly;
static void expand(int c) {
    const Tour *T = &reps[c]; int Q[16]; uint32_t cand[16];
    for (int st = 0; st < n; st++) {
        int depth = 0; Q[0] = st; uint32_t used = 1u << st; cand[0] = T->o[st] & ~used;
        while (depth >= 0) {
            if (depth == n - 1) {
                if (!blackOnly || !isBlue(T, Q)) {
                    Tour U = converse(T);
                    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; U.o[u] |= 1u << v; U.o[v] &= ~(1u << u); }
                    int id = hget(canon(&U, NULL));
                    if (id < 0) { fprintf(stderr, "missing class\n"); exit(1); }
                    if (pot[id] == pot[c] - 1) hasdesc[c] = 1;
                    if (dist[id] < 0) { dist[id] = dist[c] + 1; nxt[nnxt++] = id; }
                }
                used &= ~(1u << Q[depth]); depth--; continue;
            }
            if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
            depth++; Q[depth] = v; used |= 1u << v; cand[depth] = T->o[v] & ~used;
        }
    }
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); full = (1u << n) - 1; m_words = SETWORDSNEEDED(n); nauty_check(WORDSIZE, m_words, n, NAUTYVERSIONID);
    int cap = 1 << 18; reps = malloc(sizeof(Tour) * cap); C = 0; hcap = 1 << 20; hk = malloc(8 * hcap); hv = malloc(sizeof(int) * hcap); for (size_t i = 0; i < hcap; i++) hv[i] = -1;
    char line[256];
    while (fgets(line, sizeof line, stdin)) { if (line[0] != '0' && line[0] != '1') continue; if (C == cap) { cap *= 2; reps = realloc(reps, sizeof(Tour) * cap); }
        Tour T; memset(&T, 0, sizeof T); int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.o[i] |= 1u << j; else T.o[j] |= 1u << i; k++; }
        hput(canon(&T, NULL), C); reps[C++] = T; }
    dist = malloc(C); pot = malloc(sizeof(int) * C); tgt = malloc(sizeof(int) * C); hasdesc = malloc(C); cur = malloc(sizeof(int) * C); nxt = malloc(sizeof(int) * C);
    int *aut = malloc(sizeof(int) * C); char *sc = malloc(C);
    for (int c = 0; c < C; c++) { double g; uint64_t code = canon(&reps[c], &g); aut[c] = (int)(g + 0.5); Tour Tc = converse(&reps[c]); sc[c] = canon(&Tc, NULL) == code; }
    printf("n=%d classes=%d\n", n, C); fflush(stdout);
    for (int a = 2; a < argc; a++) {
        const char *job = argv[a]; blackOnly = !strcmp(job, "REGB") || !strcmp(job, "TTB") || !strcmp(job, "SINGB");
        for (int c = 0; c < C; c++) {
            const Tour *T = &reps[c]; int s[16], src = 0; for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(T->o[v]); if (s[v] == 0 || s[v] == n - 1) src = 1; }
            if (!strncmp(job, "REG", 3)) { int r = 0; for (int v = 0; v < n; v++) r += abs(2 * s[v] - (n - 1)); pot[c] = r / 4; tgt[c] = (r == 0); }
            else if (!strncmp(job, "TT", 2)) { int t[16]; memcpy(t, s, sizeof(int) * n); qsort(t, n, sizeof(int), cmpint); int l = 0; for (int i = 0; i < n; i++) l += abs(t[i] - i); pot[c] = l / 2; tgt[c] = (l == 0); }
            else if (!strncmp(job, "SING", 4)) { int d = 99; for (int v = 0; v < n; v++) { int x = s[v] < n - 1 - s[v] ? s[v] : n - 1 - s[v]; if (x < d) d = x; } pot[c] = d; tgt[c] = src; }
            else if (!strcmp(job, "NS")) { pot[c] = lambda_of(T); tgt[c] = (pot[c] == 0); }
            else { fprintf(stderr, "unknown job %s\n", job); return 1; }
        }
        memset(dist, -1, C); memset(hasdesc, 0, C); ncur = 0;
        for (int c = 0; c < C; c++) if (tgt[c]) { dist[c] = 0; cur[ncur++] = c; }
        int level = 0;
        while (ncur) { nnxt = 0; for (int i = 0; i < ncur; i++) expand(cur[i]); int *t = cur; cur = nxt; nxt = t; ncur = nnxt; level++; }
        /* classes never reached were never expanded: expand them now for the local-minimum test (no BFS effect) */
        long unreach = 0; for (int c = 0; c < C; c++) if (dist[c] < 0) { unreach++; }
        long below = 0, fail = 0, failSym = 0, lmin = 0, lminSym = 0, failStrongK2 = 0; int maxex = 0; static long hist[64][16];
        memset(hist, 0, sizeof hist);
        for (int c = 0; c < C; c++) {
            if (dist[c] < 0) continue;
            if (dist[c] < pot[c]) below++;
            int ex = dist[c] - pot[c]; if (ex >= 0 && ex < 16 && pot[c] < 64) hist[pot[c]][ex]++;
            if (ex > 0) { fail++; if (aut[c] > 1 || sc[c]) failSym++; if (ex > maxex) maxex = ex;
                if (strong_in(&reps[c], full) && kappa_of(&reps[c]) >= 2) failStrongK2++; }
            if (!tgt[c] && !hasdesc[c]) { lmin++; if (aut[c] > 1 || sc[c]) lminSym++; }
        }
        printf("%-5s levels %d unreachable %ld  d<pot %ld  failures %ld (symmetric %ld; strong 2-connected %ld; max excess %d)  local minima %ld (symmetric %ld)\n   (pot,+excess):count", job, level, unreach, below, fail, failSym, failStrongK2, maxex, lmin, lminSym);
        for (int p = 0; p < 64; p++) for (int e = 0; e < 16; e++) if (hist[p][e]) printf(" (%d,+%d):%ld", p, e, hist[p][e]);
        printf("\n"); fflush(stdout);
    }
    return 0;
}
