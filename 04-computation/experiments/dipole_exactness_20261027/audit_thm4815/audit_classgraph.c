/* audit_classgraph.c -- independent audit of THM-4815 section D (class-level native-move graph, orders 4..8).
 * Written from scratch for the audit.  Input: gentourng ascii on stdin.
 * Move: tau_Q(T) = T^op xor E(Q) along a Hamiltonian path Q = (q_0..q_{n-1}) (keep Q's arcs, reverse all others).
 * Colour: blue iff the positional reversal q_i -> q_{n-1-i} is an anti-automorphism of T (u->v implies r(v)->r(u)),
 * else black.  Canonical forms: nauty densenauty, DEFAULTOPTIONS_DIGRAPH, libnauty.a (dynamic MAXN).
 * Targets / potentials:
 *   TT   transitive class             L      = (1/2) sum |s_(i) - i|
 *   SING a source or a sink           delta0 = min_v min(s_v, n-1-s_v)
 *   NS   non-strong                   lambda = arc-strong connectivity
 *   REG  odd n: regular               R      = (1/2) sum |s_v - (n-1)/2|
 *        even n: almost regular       ceil(R'/2), R' = sum_v (|2 s_v - n + 1| - 1)/2
 * For each (target, colour in {any, black}): multi-source BFS distances; unreachable; d < bound (must be 0);
 * failures (d > bound) and symmetric failures (|Aut| > 1 or self-converse); local minima (reachable, off target,
 * no neighbour with potential exactly one lower); (lambda, kappa) of local minima; colour obstructions (black local
 * minima with an any-colour descending neighbour) and how many are self-converse.  Also: unit-tariff check on every
 * edge, undirectedness of both class graphs, colour symmetry (Q blue in T iff blue in tau_Q T), black components,
 * pure-blue classes, kappa of every strong failure, rigid failures, and a few named classes.
 * Build: cc -O3 -march=native -I/opt/homebrew/opt/nauty/include audit_classgraph.c /opt/homebrew/opt/nauty/lib/libnauty.a -o audit_classgraph
 * Run:   gentourng -q n | ./audit_classgraph n
 */
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, m_words;
static uint32_t full;
typedef struct { uint32_t o[16]; } Tour;

static uint64_t canon(const Tour *T, double *grp) {
    DYNALLSTAT(graph, g, g_sz); DYNALLSTAT(graph, cg, cg_sz);
    DYNALLSTAT(int, lab, lab_sz); DYNALLSTAT(int, ptn, ptn_sz); DYNALLSTAT(int, orbits, orbits_sz);
    static DEFAULTOPTIONS_DIGRAPH(options); statsblk stats;
    int m = m_words;
    DYNALLOC2(graph, g, g_sz, m, n, "canon"); DYNALLOC2(graph, cg, cg_sz, m, n, "canon");
    DYNALLOC1(int, lab, lab_sz, n, "canon"); DYNALLOC1(int, ptn, ptn_sz, n, "canon"); DYNALLOC1(int, orbits, orbits_sz, n, "canon");
    options.getcanon = TRUE;
    EMPTYGRAPH(g, m, n);
    for (int u = 0; u < n; u++) for (uint32_t x = T->o[u]; x; x &= x - 1) ADDONEARC(g, u, __builtin_ctz(x), m);
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, cg);
    uint64_t code = 0; int b = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (ISELEMENT(GRAPHROW(cg, i, m), j)) code |= (uint64_t)1 << b; b++; }
    if (grp) { double s = stats.grpsize1; for (int e = 0; e < stats.grpsize2; e++) s *= 10.0; *grp = s; }
    return code;
}
static Tour converse(const Tour *T) { Tour R; for (int u = 0; u < n; u++) R.o[u] = ~T->o[u] & full & ~(1u << u); for (int u = n; u < 16; u++) R.o[u] = 0; return R; }
static Tour tau(const Tour *T, const int *Q) {
    Tour R = converse(T);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; R.o[u] |= 1u << v; R.o[v] &= ~(1u << u); }
    return R;
}
static int isBlue(const Tour *T, const int *Q) {
    int r[16]; for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) for (uint32_t x = T->o[u]; x; x &= x - 1) { int v = __builtin_ctz(x); if (!(T->o[r[v]] >> r[u] & 1)) return 0; }
    return 1;
}
/* --- hash map uint64 -> int --- */
static uint64_t *hk; static int *hv; static size_t hcap;
static size_t hh(uint64_t k) { k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33; return (size_t)k; }
static void hput(uint64_t k, int v) { size_t i = hh(k) & (hcap - 1); while (hv[i] >= 0) { if (hk[i] == k) { hv[i] = v; return; } i = (i + 1) & (hcap - 1); } hk[i] = k; hv[i] = v; }
static int hget(uint64_t k) { size_t i = hh(k) & (hcap - 1); while (hv[i] >= 0) { if (hk[i] == k) return hv[i]; i = (i + 1) & (hcap - 1); } return -1; }
/* --- structure --- */
static uint32_t reach_in(const Tour *T, uint32_t vs, int v) { uint32_t seen = 1u << v, fr = seen; while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= T->o[__builtin_ctz(f)] & vs; nx &= ~seen; seen |= nx; fr = nx; } return seen; }
static int strong_in(const Tour *T, uint32_t vs) { int pc = __builtin_popcount(vs); if (pc <= 1) return 1; if (pc == 2) return 0; for (uint32_t w = vs; w; w &= w - 1) if (reach_in(T, vs, __builtin_ctz(w)) != vs) return 0; return 1; }
static int kappa_of(const Tour *T) { if (!strong_in(T, full)) return 0; for (int k = 1; k <= n - 2; k++) for (uint32_t S = 1; S < full; S++) if (__builtin_popcount(S) == k && !strong_in(T, full & ~S)) return k; return n - 1; }
static int lambda_of(const Tour *T) { int best = 1 << 20; for (uint32_t X = 1; X < full; X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T->o[x] & ~X & full); if (f < best) best = f; } return best; }
static int cmpint(const void *a, const void *b) { return *(const int *)a - *(const int *)b; }

typedef struct { int *a; int len, cap; } Vec;
static void vpush(Vec *v, int x) { if (v->len == v->cap) { v->cap = v->cap ? 2 * v->cap : 8; v->a = realloc(v->a, sizeof(int) * v->cap); } v->a[v->len++] = x; }
static int cmpi(const void *a, const void *b) { int x = *(const int *)a, y = *(const int *)b; return (x > y) - (x < y); }
static void vuniq(Vec *v) { if (!v->len) return; qsort(v->a, v->len, sizeof(int), cmpi); int k = 1; for (int i = 1; i < v->len; i++) if (v->a[i] != v->a[k - 1]) v->a[k++] = v->a[i]; v->len = k; }
static int vhas(const Vec *v, int x) { int lo = 0, hi = v->len - 1; while (lo <= hi) { int md = (lo + hi) / 2; if (v->a[md] == x) return 1; if (v->a[md] < x) lo = md + 1; else hi = md - 1; } return 0; }

static int C; static Tour *reps; static Vec *adjA, *adjB;
static int *aut, *sc, *lam, *kap, *L2, *d0, *R4, *Rev, *strongv, *nHP, *nBlackHP, *selfloopBlack;
static char **arcstr;
static int *ufp; static int uf(int x) { while (ufp[x] != x) { ufp[x] = ufp[ufp[x]]; x = ufp[x]; } return x; }

int main(int argc, char **argv) {
    n = atoi(argv[1]); full = (1u << n) - 1; m_words = SETWORDSNEEDED(n); nauty_check(WORDSIZE, m_words, n, NAUTYVERSIONID);
    int capC = 1 << 16; reps = malloc(sizeof(Tour) * capC); arcstr = malloc(sizeof(char *) * capC); C = 0;
    hcap = 1 << 18; hk = malloc(sizeof(uint64_t) * hcap); hv = malloc(sizeof(int) * hcap); for (size_t i = 0; i < hcap; i++) hv[i] = -1;
    char line[256];
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        if (C == capC) { capC *= 2; reps = realloc(reps, sizeof(Tour) * capC); arcstr = realloc(arcstr, sizeof(char *) * capC); }
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.o[i] |= 1u << j; else T.o[j] |= 1u << i; k++; }
        uint64_t code = canon(&T, NULL);
        if (hget(code) >= 0) { fprintf(stderr, "duplicate class in input\n"); return 1; }
        hput(code, C); reps[C] = T; line[strcspn(line, "\n")] = 0; arcstr[C] = strdup(line); C++;
    }
    adjA = calloc(C, sizeof(Vec)); adjB = calloc(C, sizeof(Vec));
    aut = calloc(C, sizeof(int)); sc = calloc(C, sizeof(int)); lam = calloc(C, sizeof(int)); kap = calloc(C, sizeof(int));
    L2 = calloc(C, sizeof(int)); d0 = calloc(C, sizeof(int)); R4 = calloc(C, sizeof(int)); Rev = calloc(C, sizeof(int));
    strongv = calloc(C, sizeof(int)); nHP = calloc(C, sizeof(int)); nBlackHP = calloc(C, sizeof(int)); selfloopBlack = calloc(C, sizeof(int));
    long colourAsym = 0, involFail = 0, missing = 0;
    for (int c = 0; c < C; c++) {
        Tour *T = &reps[c]; double g; uint64_t code = canon(T, &g); aut[c] = (int)(g + 0.5);
        Tour Tc = converse(T); sc[c] = canon(&Tc, NULL) == code;
        lam[c] = lambda_of(T); kap[c] = kappa_of(T); strongv[c] = lam[c] > 0;
        int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(T->o[v]);
        int t[16]; memcpy(t, s, sizeof(int) * n); qsort(t, n, sizeof(int), cmpint);
        int l2 = 0; for (int i = 0; i < n; i++) l2 += abs(t[i] - i); L2[c] = l2;
        int dd = 1 << 20, r4 = 0, rev = 0;
        for (int v = 0; v < n; v++) { int x = s[v] < n - 1 - s[v] ? s[v] : n - 1 - s[v]; if (x < dd) dd = x; r4 += abs(2 * s[v] - (n - 1)); rev += (abs(2 * s[v] - n + 1) - 1) / 2; }
        d0[c] = dd; R4[c] = r4; Rev[c] = rev;
        /* enumerate Hamiltonian paths by DFS */
        int Q[16]; uint32_t cand[16]; int depth;
        for (int st = 0; st < n; st++) {
            depth = 0; Q[0] = st; uint32_t used = 1u << st; cand[0] = T->o[st] & ~used;
            if (n == 1) { break; }
            while (depth >= 0) {
                if (depth == n - 1) {
                    nHP[c]++;
                    Tour U = tau(T, Q);
                    /* A1 checks: Q is a Hamiltonian path of U, tau_Q(U) = T, colour symmetric */
                    Tour W = tau(&U, Q); if (memcmp(&W, T, sizeof(Tour))) involFail++;
                    int blue = isBlue(T, Q); if (blue != isBlue(&U, Q)) colourAsym++;
                    int id = hget(canon(&U, NULL)); if (id < 0) { missing++; }
                    else { vpush(&adjA[c], id); if (!blue) { nBlackHP[c]++; vpush(&adjB[c], id); if (id == c) selfloopBlack[c] = 1; } }
                    used &= ~(1u << Q[depth]); depth--; continue;
                }
                if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
                int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
                depth++; Q[depth] = v; used |= 1u << v; cand[depth] = T->o[v] & ~used;
            }
        }
        vuniq(&adjA[c]); vuniq(&adjB[c]);
    }
    long asymA = 0, asymB = 0; for (int c = 0; c < C; c++) { for (int i = 0; i < adjA[c].len; i++) if (!vhas(&adjA[adjA[c].a[i]], c)) asymA++; for (int i = 0; i < adjB[c].len; i++) if (!vhas(&adjB[adjB[c].a[i]], c)) asymB++; }
    printf("n=%d classes=%d  moves: involution failures %ld, colour asymmetries %ld, missing classes %ld; class graph asymmetric edges any %ld black %ld\n", n, C, involFail, colourAsym, missing, asymA, asymB);
    /* black components */
    ufp = malloc(sizeof(int) * C); for (int c = 0; c < C; c++) ufp[c] = c;
    for (int c = 0; c < C; c++) for (int i = 0; i < adjB[c].len; i++) { int a = uf(c), b = uf(adjB[c].a[i]); if (a != b) ufp[a] = b; }
    int *csize = calloc(C, sizeof(int)); for (int c = 0; c < C; c++) csize[uf(c)]++;
    int nbig = 0, nsingle = 0, pureblue = 0, pbSingle = 0, singleNotPB = 0;
    for (int c = 0; c < C; c++) { if (uf(c) == c) { if (csize[c] > 1) nbig++; else nsingle++; } }
    for (int c = 0; c < C; c++) { int pb = (nBlackHP[c] == 0); pureblue += pb; if (csize[uf(c)] == 1) { if (pb) pbSingle++; else singleNotPB++; } }
    printf("black class graph: components of size > 1: %d; isolated classes: %d (pure blue %d, isolated but not pure blue %d); pure-blue classes total %d\n", nbig, nsingle, pbSingle, singleNotPB, pureblue);
    for (int c = 0; c < C; c++) if (nBlackHP[c] == 0) { int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(reps[c].o[v]); qsort(s, n, sizeof(int), cmpint); printf("   pure blue: %s scores", arcstr[c]); for (int v = 0; v < n; v++) printf(" %d", s[v]); printf("  |Aut| %d SC %d HP %d\n", aut[c], sc[c], nHP[c]); }
    /* targets & potentials */
    const char *tname[4] = { "TT/L", "SING/delta0", "NS/lambda", (n & 1) ? "REG/R" : "REG(almost)/ceil(R'/2)" };
    int *pot = malloc(sizeof(int) * C), *tgt = malloc(sizeof(int) * C), *dist = malloc(sizeof(int) * C), *queue = malloc(sizeof(int) * C);
    for (int job = 0; job < 4; job++) {
        for (int c = 0; c < C; c++) {
            int s[16]; int src = 0; for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(reps[c].o[v]); if (s[v] == 0 || s[v] == n - 1) src = 1; }
            switch (job) {
                case 0: pot[c] = L2[c] / 2; tgt[c] = (L2[c] == 0); break;
                case 1: pot[c] = d0[c]; tgt[c] = src; break;
                case 2: pot[c] = lam[c]; tgt[c] = (lam[c] == 0); break;
                case 3: if (n & 1) { pot[c] = R4[c] / 4; tgt[c] = (R4[c] == 0); } else { pot[c] = (Rev[c] + 1) / 2; tgt[c] = (Rev[c] == 0); } break;
            }
            if (tgt[c] != (pot[c] == 0)) { printf("   WARNING %s: target/potential-zero mismatch at class %d\n", tname[job], c); }
        }
        long tariff = 0; for (int c = 0; c < C; c++) for (int i = 0; i < adjA[c].len; i++) if (abs(pot[adjA[c].a[i]] - pot[c]) > 1) tariff++;
        for (int col = 0; col < 2; col++) {
            Vec *adj = col ? adjB : adjA;
            int qh = 0, qt = 0; for (int c = 0; c < C; c++) { dist[c] = -1; if (tgt[c]) { dist[c] = 0; queue[qt++] = c; } }
            while (qh < qt) { int c = queue[qh++]; for (int i = 0; i < adj[c].len; i++) { int t = adj[c].a[i]; if (dist[t] < 0) { dist[t] = dist[c] + 1; queue[qt++] = t; } } }
            int unr = 0, below = 0, fail = 0, failSym = 0, lmin = 0, lminSym = 0, colour = 0, colourSC = 0, failStrongK2 = 0, failStrong = 0, failRigid = 0, failRigidStrong = 0;
            int hist[20][20]; memset(hist, 0, sizeof hist); int maxexcess = 0;
            for (int c = 0; c < C; c++) {
                if (dist[c] < 0) { unr++; continue; }
                if (dist[c] < pot[c]) below++;
                if (dist[c] > pot[c]) { fail++; if (aut[c] > 1 || sc[c]) failSym++; else { failRigid++; if (strongv[c]) failRigidStrong++; }
                    if (dist[c] - pot[c] > maxexcess) maxexcess = dist[c] - pot[c];
                    if (strongv[c]) { failStrong++; if (kap[c] >= 2) { failStrongK2++;
                        if (failStrongK2 <= 12) { int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(reps[c].o[v]); qsort(s, n, sizeof(int), cmpint);
                            int dsc = 0; for (int i = 0; i < adj[c].len; i++) if (pot[adj[c].a[i]] == pot[c] - 1) dsc = 1;
                            printf("   2-connected strong failure (%s, %s): %s scores", tname[job], col ? "black" : "any", arcstr[c]); for (int v = 0; v < n; v++) printf(" %d", s[v]);
                            printf(" |Aut| %d SC %d lambda %d kappa %d pot %d dist %d has-descending-neighbour %d\n", aut[c], sc[c], lam[c], kap[c], pot[c], dist[c], dsc); } } } }
                if (!tgt[c]) {
                    int desc = 0; for (int i = 0; i < adj[c].len; i++) if (pot[adj[c].a[i]] == pot[c] - 1) { desc = 1; break; }
                    if (!desc) { lmin++; if (aut[c] > 1 || sc[c]) lminSym++; hist[lam[c]][kap[c]]++;
                        if (col) { int ad = 0; for (int i = 0; i < adjA[c].len; i++) if (pot[adjA[c].a[i]] == pot[c] - 1) { ad = 1; break; } if (ad) { colour++; colourSC += sc[c];
                            int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(reps[c].o[v]); qsort(s, n, sizeof(int), cmpint);
                            printf("   colour obstruction (%s): %s scores", tname[job], arcstr[c]); for (int v = 0; v < n; v++) printf(" %d", s[v]); printf(" |Aut| %d SC %d lambda %d kappa %d pot %d dist %d\n", aut[c], sc[c], lam[c], kap[c], pot[c], dist[c]); } } }
                }
            }
            printf("%-24s %-5s unreachable %d  d<bound %d  tariff>1 edges %ld  failures %d (symmetric %d, rigid %d of which strong %d; strong %d, strong kappa>=2 %d; max excess %d)  local minima %d (symmetric %d)",
                   tname[job], col ? "black" : "any", unr, below, tariff, fail, failSym, failRigid, failRigidStrong, failStrong, failStrongK2, maxexcess, lmin, lminSym);
            printf("  LM by (lambda,kappa):"); for (int a = 0; a < 20; a++) for (int b = 0; b < 20; b++) if (hist[a][b]) printf(" (%d,%d):%d", a, b, hist[a][b]);
            if (col) printf("  colour obstructions %d (self-converse %d)", colour, colourSC);
            printf("\n");
            if (col == 0) {   /* breakdown of local minima and failures by (|Aut| = 1 ?, self-converse ?) */
                int bLM[2][2] = {{0, 0}, {0, 0}}, bF[2][2] = {{0, 0}, {0, 0}};
                for (int c = 0; c < C; c++) { if (dist[c] < 0) continue; int r = aut[c] == 1, s2 = sc[c] ? 1 : 0;
                    if (dist[c] > pot[c]) bF[r][s2]++;
                    if (!tgt[c]) { int desc = 0; for (int i = 0; i < adj[c].len; i++) if (pot[adj[c].a[i]] == pot[c] - 1) { desc = 1; break; } if (!desc) bLM[r][s2]++; } }
                printf("   breakdown %s any: local minima [|Aut|=1 & not SC %d, |Aut|=1 & SC %d, |Aut|>1 & not SC %d, |Aut|>1 & SC %d]; failures [%d, %d, %d, %d]\n", tname[job],
                       bLM[1][0], bLM[1][1], bLM[0][0], bLM[0][1], bF[1][0], bF[1][1], bF[0][0], bF[0][1]);
            }
            if (job == 0 && col == 0 && n <= 6) {   /* list rigid L-failures */
                for (int c = 0; c < C; c++) if (dist[c] > pot[c] && !(aut[c] > 1 || sc[c])) { int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(reps[c].o[v]); qsort(s, n, sizeof(int), cmpint);
                    printf("   rigid L-failure: %s scores", arcstr[c]); for (int v = 0; v < n; v++) printf(" %d", s[v]); printf(" strong %d lambda %d kappa %d L %d dist %d\n", strongv[c], lam[c], kap[c], pot[c], dist[c]); }
            }
            if (job == 0 && col == 0) { /* named class */
                for (int c = 0; c < C; c++) if (!strcmp(arcstr[c], "111011111111101") && n == 6) printf("   named 111011111111101: |Aut| %d SC %d strong %d lambda %d kappa %d L %d dist(TT) %d HP %d\n", aut[c], sc[c], strongv[c], lam[c], kap[c], pot[c], dist[c], nHP[c]);
            }
            if (job == 3 && (n & 1)) { /* scores 1,3^5,5 */
                for (int c = 0; c < C; c++) { int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(reps[c].o[v]); qsort(s, n, sizeof(int), cmpint);
                    if (n == 7 && s[0] == 1 && s[6] == 5 && s[1] == 3 && s[5] == 3) printf("   scores 1,3^5,5: %s |Aut| %d SC %d strong %d kappa %d R %d dist(%s) %d\n", arcstr[c], aut[c], sc[c], strongv[c], kap[c], pot[c], col ? "black" : "any", dist[c]); }
            }
        }
    }
    long hp = 0; for (int c = 0; c < C; c++) hp += nHP[c];
    printf("total Hamiltonian paths over class representatives: %ld\n", hp);
    return 0;
}
