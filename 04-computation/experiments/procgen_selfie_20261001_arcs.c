/* procgen_selfie_20261001_arcs.c  (procgen selfie lane, 2026-10-01)
 *
 * Arc-Hamiltonian-path counts c(e) = #{Hamiltonian paths through arc e} for tournaments,
 * plus the "shaving" invariants of a tournament.  Prints to stdout only.
 *
 * Modes
 *   stats N            read gentourng ASCII lines (upper triangle, row by row; bit 1 means i->j)
 *                      from stdin; aggregate per-class statistics, also weighted by N!/|Aut|
 *   lab N              all 2^C(N,2) labeled tournaments (N <= 7), same statistics, weight 1
 *   rnd N COUNT SEED   COUNT uniform random labeled tournaments
 *   shave N            read gentourng lines (N <= 7); per class: HP-blocking number beta,
 *                      two-source/two-sink bound, parity-break number rho, Redei-shavability
 *   anf N              truth table of c(0->1) mod 2 over all labeled tournaments containing 0->1
 *                      (other C(N,2)-1 arcs free; bit k of the index = arc k in the standard order
 *                      after removing arc (0,1)), printed as a 0/1 string
 *   search N ITERS SEED  random-restart local search for an all-odd tournament (prints any found)
 *   beta N             read gentourng lines (N <= 11); branch-and-bound test of
 *                      beta(T) = min(sigma2-(T), sigma2+(T)) (HP-blocking number = two-source/two-sink bound);
 *                      prints a summary line
 *   list N             read gentourng lines; print 'T H #odd #zero strong |Aut|' per class
 *   minodd N ITERS SEED  local search (from near-transitive starts) minimising the number of odd arcs;
 *                      prints the smallest number found and a witness (tests Conjecture C4 for even N)
 *   oddfilter N K      read gentourng lines; mod-2 bitmask DP (independent of hp_counts) for the arc
 *                      parities; print every class with at most K odd arcs; final #odd histogram
 *   hlab N             H(T) for every labeled tournament, one per line, in index order
 *                      (bit b of the index = arc b of the standard order; 1 means i->j, i<j)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 20
typedef unsigned long long u64;

static int N;
static int A[MAXN][MAXN];          /* A[i][j]=1 iff i->j */
static u64 *fbuf, *gbuf;           /* (1<<N) x N arrays, allocated in main */
#define f(S, v) fbuf[(size_t)(S) * N + (v)]   /* paths covering S ending at v   */
#define g(S, v) gbuf[(size_t)(S) * N + (v)]   /* paths covering S starting at v */
static u64 C[MAXN][MAXN];          /* C[u][v] = #HPs through arc u->v */
static u64 Hval, startc[MAXN], endc[MAXN];
static int only_H = 0;   /* hlab mode: skip the arc counts */

static void hp_counts(void) {
    int full = (1 << N) - 1;
    for (int S = 0; S <= full; S++) for (int v = 0; v < N; v++) { f(S, v) = 0; g(S, v) = 0; }
    for (int v = 0; v < N; v++) { f(1 << v, v) = 1; g(1 << v, v) = 1; }
    for (int S = 1; S <= full; S++) {
        for (int v = 0; v < N; v++) {
            if (!((S >> v) & 1)) continue;
            u64 fv = f(S, v), gv = g(S, v);
            if (!fv && !gv) continue;
            for (int w = 0; w < N; w++) {
                if ((S >> w) & 1) continue;
                if (fv && A[v][w]) f(S | (1 << w), w) += fv;
                if (gv && A[w][v]) g(S | (1 << w), w) += gv;
            }
        }
    }
    Hval = 0;
    for (int v = 0; v < N; v++) { Hval += f(full, v); endc[v] = f(full, v); startc[v] = g(full, v); }
    if (only_H) return;
    for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) C[u][v] = 0;
    for (int S = 1; S < full; S++) {
        int Cm = full ^ S;
        for (int u = 0; u < N; u++) {
            if (!((S >> u) & 1)) continue;
            u64 fu = f(S, u);
            if (!fu) continue;
            for (int v = 0; v < N; v++) {
                if (!((Cm >> v) & 1) || !A[u][v]) continue;
                C[u][v] += fu * g(Cm, v);
            }
        }
    }
}

/* automorphism count by backtracking */
static int perm[MAXN], usedv[MAXN], outdeg[MAXN], inv2[MAXN];
static u64 autcount;
static void aut_bt(int i) {
    if (i == N) { autcount++; return; }
    for (int w = 0; w < N; w++) {
        if (usedv[w] || outdeg[w] != outdeg[i] || inv2[w] != inv2[i]) continue;
        int ok = 1;
        for (int j = 0; j < i && ok; j++) if (A[i][j] != A[w][perm[j]]) ok = 0;
        if (!ok) continue;
        usedv[w] = 1; perm[i] = w;
        aut_bt(i + 1);
        usedv[w] = 0;
    }
}
static u64 aut_size(void) {
    for (int v = 0; v < N; v++) { outdeg[v] = 0; for (int w = 0; w < N; w++) outdeg[v] += A[v][w]; }
    for (int v = 0; v < N; v++) { inv2[v] = 0; for (int w = 0; w < N; w++) if (A[v][w]) inv2[v] += outdeg[w]; }
    for (int v = 0; v < N; v++) usedv[v] = 0;
    autcount = 0; aut_bt(0); return autcount;
}

static int is_strong(void) {
    /* reachability from 0 forwards and backwards */
    int seen = 1, stack[MAXN], sp = 0; stack[sp++] = 0;
    while (sp) { int v = stack[--sp]; for (int w = 0; w < N; w++) if (A[v][w] && !((seen >> w) & 1)) { seen |= 1 << w; stack[sp++] = w; } }
    if (seen != (1 << N) - 1) return 0;
    seen = 1; sp = 0; stack[sp++] = 0;
    while (sp) { int v = stack[--sp]; for (int w = 0; w < N; w++) if (A[w][v] && !((seen >> w) & 1)) { seen |= 1 << w; stack[sp++] = w; } }
    return seen == (1 << N) - 1;
}
static int n_components(void) { /* number of strong components = number of distinct reach-classes */
    int reach[MAXN];
    for (int s = 0; s < N; s++) {
        int seen = 1 << s, stack[MAXN], sp = 0; stack[sp++] = s;
        while (sp) { int v = stack[--sp]; for (int w = 0; w < N; w++) if (A[v][w] && !((seen >> w) & 1)) { seen |= 1 << w; stack[sp++] = w; } }
        reach[s] = seen;
    }
    int cnt = 0, done = 0;
    for (int s = 0; s < N; s++) {
        if ((done >> s) & 1) continue;
        int comp = 0;
        for (int t = 0; t < N; t++) if (((reach[s] >> t) & 1) && ((reach[t] >> s) & 1)) comp |= 1 << t;
        done |= comp; cnt++;
    }
    return cnt;
}

/* ---------------- statistics ---------------- */
#define NE_MAX 190
typedef struct {
    double cls, lab;
} acc;
static acc tot, allpos, allodd, alleven, allequal, strong, strong_allpos, some_on_all, comp_le2;
static acc zhist[NE_MAX + 1], oddhist[NE_MAX + 1], ocomp[MAXN + 1], oiso[MAXN + 1];
static long id_fail = 0, lemma_fail = 0;
static int printed_odd = 0, printed_eq = 0, printed_strong_bad = 0;

static void addacc(acc *a, double w) { a->cls += 1; a->lab += w; }
static u64 maxH = 0; static int nmax = 0, maxodd[64];

static void process(const char *label, double w, int verbose) {
    hp_counts();
    int ne = N * (N - 1) / 2, z = 0, odd = 0, eq = 1, on_all = 0;
    u64 c0 = 0; int first = 1;
    for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v]) {
        u64 c = C[u][v];
        if (c == 0) z++;
        if (c & 1) odd++;
        if (first) { c0 = c; first = 0; } else if (c != c0) eq = 0;
        if (c == Hval) on_all = 1;
    }
    /* identities: sum over out-arcs = H - end(v); sum over in-arcs = H - start(v) */
    u64 se = 0, ss = 0;
    for (int v = 0; v < N; v++) {
        u64 so = 0, si = 0;
        for (int w = 0; w < N; w++) { if (A[v][w]) so += C[v][w]; if (A[w][v]) si += C[w][v]; }
        if (so != Hval - endc[v] || si != Hval - startc[v]) id_fail++;
        se += endc[v]; ss += startc[v];
    }
    if (se != Hval || ss != Hval || !(Hval & 1)) id_fail++;
    int st = is_strong();
    int k = n_components();
    /* lemma: all arcs on some HP  <=>  k <= 2 and every strong component's arcs are on HPs of the component.
       Here we check the necessary direction cheaply: z == 0 implies k <= 2. */
    if (z == 0 && k > 2) lemma_fail++;
    if (Hval > maxH) { maxH = Hval; nmax = 0; }
    if (Hval == maxH && nmax < 64) maxodd[nmax++] = odd;
    addacc(&tot, w);
    if (z == 0) addacc(&allpos, w);
    if (odd == ne) addacc(&allodd, w);
    if (odd == 0) addacc(&alleven, w);
    if (eq) addacc(&allequal, w);
    if (st) { addacc(&strong, w); if (z == 0) addacc(&strong_allpos, w); }
    if (on_all) addacc(&some_on_all, w);
    if (k <= 2) addacc(&comp_le2, w);
    addacc(&zhist[z], w);
    addacc(&oddhist[odd], w);
    { /* odd-arc graph O: number of connected components and isolated vertices */
        int comp[MAXN], nc = 0, iso = 0;
        for (int v = 0; v < N; v++) comp[v] = -1;
        for (int s0 = 0; s0 < N; s0++) {
            if (comp[s0] >= 0) continue;
            int stack[MAXN], sp = 0; stack[sp++] = s0; comp[s0] = nc;
            while (sp) { int v = stack[--sp];
                for (int x = 0; x < N; x++) if (comp[x] < 0 && ((A[v][x] && (C[v][x] & 1)) || (A[x][v] && (C[x][v] & 1)))) { comp[x] = nc; stack[sp++] = x; } }
            nc++;
        }
        for (int v = 0; v < N; v++) { int d = 0; for (int x = 0; x < N; x++) if ((A[v][x] && (C[v][x] & 1)) || (A[x][v] && (C[x][v] & 1))) d++; if (!d) iso++; }
        addacc(&ocomp[nc], w); addacc(&oiso[iso], w);
    }
    if (verbose && odd == ne && printed_odd < 40) {
        printed_odd++;
        int sc[MAXN]; for (int v = 0; v < N; v++) { sc[v] = 0; for (int x = 0; x < N; x++) sc[v] += A[v][x]; }
        printf("ALLODD N=%d T=%s H=%llu scores=", N, label, Hval);
        for (int v = 0; v < N; v++) printf("%d%s", sc[v], v < N - 1 ? "," : "");
        printf(" weight=%.0f cvals=", w);
        for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v]) printf("%llu.", C[u][v]);
        printf("\n");
    }
    if (verbose && eq && printed_eq < 20) {
        printed_eq++;
        printf("ALLEQUAL N=%d T=%s H=%llu c=%llu weight=%.0f\n", N, label, Hval, c0, w);
    }
    if (verbose && st && z > 0 && printed_strong_bad < 3) {
        printed_strong_bad++;
        printf("STRONG_WITH_DEAD_ARC N=%d T=%s H=%llu dead=", N, label, Hval);
        for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v] && C[u][v] == 0) printf("%d>%d ", u, v);
        printf("\n");
    }
}

static void report(void) {
    int ne = N * (N - 1) / 2;
    printf("SUMMARY N=%d classes=%.0f labeled=%.0f\n", N, tot.cls, tot.lab);
#define P(name, a) printf("  %-22s classes=%.0f labeled=%.0f\n", name, a.cls, a.lab)
    P("all_arcs_on_some_HP", allpos);
    P("all_arcs_odd", allodd);
    P("all_arcs_even", alleven);
    P("all_arcs_equal", allequal);
    P("strong", strong);
    P("strong_and_all_on_HP", strong_allpos);
    P("<=2_strong_components", comp_le2);
    P("some_arc_on_every_HP", some_on_all);
    printf("  zero-arc histogram (z: classes/labeled):");
    for (int z = 0; z <= ne; z++) if (zhist[z].cls > 0) printf(" %d:%.0f/%.0f", z, zhist[z].cls, zhist[z].lab);
    printf("\n  odd-arc histogram (#odd: classes/labeled):");
    for (int z = 0; z <= ne; z++) if (oddhist[z].cls > 0) printf(" %d:%.0f/%.0f", z, oddhist[z].cls, oddhist[z].lab);
    printf("\n  odd-arc-graph components (k: classes/labeled):");
    for (int z = 0; z <= N; z++) if (ocomp[z].cls > 0) printf(" %d:%.0f/%.0f", z, ocomp[z].cls, ocomp[z].lab);
    printf("\n  odd-arc-graph isolated vertices (k: classes/labeled):");
    for (int z = 0; z <= N; z++) if (oiso[z].cls > 0) printf(" %d:%.0f/%.0f", z, oiso[z].cls, oiso[z].lab);
    printf("\n  max_H=%llu maximizing_classes=%d odd_arcs_of_maximizers:", maxH, nmax);
    for (int i = 0; i < nmax; i++) printf(" %d", maxodd[i]);
    printf("\n  identity_failures=%ld lemma_failures=%ld\n", id_fail, lemma_fail);
}

static double factorial(int n) { double r = 1; for (int i = 2; i <= n; i++) r *= i; return r; }

static void set_from_string(const char *s) {
    int k = 0;
    for (int i = 0; i < N; i++) { A[i][i] = 0; for (int j = i + 1; j < N; j++) { int b = s[k++] == '1'; A[i][j] = b; A[j][i] = !b; } }
}

/* Hall-deficiency bound: min over X (|X|>=2) and Y (|Y| <= |X|-2) of #arcs entering X from outside Y
   (in-version) and the out-version.  Deleting those arcs leaves the predecessor (successor) bipartite
   graph with deficiency >= 2, hence no Hamiltonian path. */
static int match_bound(void) {
    int best = 1 << 30;
    for (int dir = 0; dir < 2; dir++) {
        for (int X = 0; X < (1 << N); X++) {
            int k = __builtin_popcount(X);
            if (k < 2) continue;
            int a[MAXN], tot = 0;
            for (int y = 0; y < N; y++) {
                a[y] = 0;
                for (int x = 0; x < N; x++) if ((X >> x) & 1) { int arc = dir == 0 ? A[y][x] : A[x][y]; a[y] += arc; }
                tot += a[y];
            }
            /* subtract the k-2 largest a[y] */
            for (int r = 0; r < k - 2; r++) { int bi = -1; for (int y = 0; y < N; y++) if (a[y] >= 0 && (bi < 0 || a[y] > a[bi])) bi = y; tot -= a[bi]; a[bi] = -1; }
            if (tot < best) best = tot;
        }
    }
    return best;
}

/* ---------------- shaving (N <= 7) ---------------- */
static uint32_t *tab; /* tab[mask] = #HPs whose arc set is contained in mask (after zeta) */
static uint8_t *good;
static int arcidx[MAXN][MAXN];

static void enum_hps(int *path, int len, int used, int ne) {
    if (len == N) {
        uint32_t m = 0;
        for (int i = 0; i + 1 < N; i++) m |= 1u << arcidx[path[i]][path[i + 1]];
        tab[m] += 1; return;
    }
    int last = path[len - 1];
    for (int w = 0; w < N; w++) if (!((used >> w) & 1) && A[last][w]) { path[len] = w; enum_hps(path, len + 1, used | (1 << w), ne); }
}

static void shave_one(const char *label) {
    int ne = N * (N - 1) / 2, k = 0;
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { arcidx[i][j] = k; arcidx[j][i] = k; k++; }
    uint32_t M = 1u << ne;
    memset(tab, 0, sizeof(uint32_t) * M);
    int path[MAXN];
    for (int s = 0; s < N; s++) { path[0] = s; enum_hps(path, 1, 1 << s, ne); }
    for (int b = 0; b < ne; b++) for (uint32_t m = 0; m < M; m++) if ((m >> b) & 1) tab[m] += tab[m ^ (1u << b)];
    uint32_t full = M - 1, H = tab[full];
    /* beta: min |F| with H(T - F) = 0 ; rho: min |F| with H(T-F) even */
    int beta = 99, rho = 99;
    for (uint32_t m = 0; m < M; m++) {
        int del = ne - __builtin_popcount(m);
        if (tab[m] == 0 && del < beta) beta = del;
        if (!(tab[m] & 1) && del < rho) rho = del;
    }
    /* two-source / two-sink bound */
    int indeg[MAXN], od[MAXN];
    for (int v = 0; v < N; v++) { indeg[v] = 0; od[v] = 0; for (int w = 0; w < N; w++) { indeg[v] += A[w][v]; od[v] += A[v][w]; } }
    int loc = 99;
    for (int u = 0; u < N; u++) for (int v = u + 1; v < N; v++) {
        int a = indeg[u] + indeg[v], b = od[u] + od[v];
        if (a < loc) loc = a; if (b < loc) loc = b;
    }
    /* Redei shaving: good[m] = H(m) odd and (H(m)==1 or exists arc e in m with good[m-e]) */
    for (uint32_t m = 0; m < M; m++) {
        good[m] = 0;
        if (!(tab[m] & 1)) continue;
        if (tab[m] == 1) { good[m] = 1; continue; }
        for (int b = 0; b < ne; b++) if (((m >> b) & 1) && good[m ^ (1u << b)]) { good[m] = 1; break; }
    }
    /* number of arcs with c odd, and core size */
    int oddc = 0, core = 0;
    for (int b = 0; b < ne; b++) { uint32_t c = H - tab[full ^ (1u << b)]; if (c & 1) oddc++; if (c) core++; }
    printf("SHAVE N=%d T=%s H=%u core=%d oddarcs=%d beta=%d twosrc=%d hall=%d rho=%d redei_shavable=%d\n",
           N, label, H, core, oddc, beta, loc, match_bound(), rho, good[full]);
}

/* ---------------- mod-2 arc parities by bitmask DP ---------------- */
static uint32_t *f2, *g2;   /* f2[S] bit v = parity of #paths covering S ending at v; g2 = starting at v */
static uint32_t outmask[MAXN], inmask[MAXN];
static int odd_arcs_mod2(int *cnt_out) {
    uint32_t M = 1u << N, full = M - 1;
    for (int v = 0; v < N; v++) { outmask[v] = inmask[v] = 0; }
    for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v]) { outmask[u] |= 1u << v; inmask[v] |= 1u << u; }
    memset(f2, 0, sizeof(uint32_t) * M); memset(g2, 0, sizeof(uint32_t) * M);
    for (int v = 0; v < N; v++) { f2[1u << v] = 1u << v; g2[1u << v] = 1u << v; }
    for (uint32_t S = 1; S < M; S++) {
        uint32_t F = f2[S], G = g2[S], rest = full & ~S;
        while (rest) {
            int w = __builtin_ctz(rest); rest &= rest - 1;
            if (__builtin_popcount(F & inmask[w]) & 1) f2[S | (1u << w)] ^= 1u << w;
            if (__builtin_popcount(G & outmask[w]) & 1) g2[S | (1u << w)] ^= 1u << w;
        }
    }
    static uint8_t par[MAXN][MAXN];
    memset(par, 0, sizeof par);
    for (uint32_t S = 1; S < full; S++) {
        uint32_t F = f2[S]; if (!F) continue;
        uint32_t G = g2[full ^ S]; if (!G) continue;
        while (F) { int u = __builtin_ctz(F); F &= F - 1; uint32_t t = G & outmask[u];
            while (t) { int v = __builtin_ctz(t); t &= t - 1; par[u][v] ^= 1; } }
    }
    int odd = 0;
    for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v]) odd += par[u][v];
    *cnt_out = __builtin_popcount(f2[full]) & 1;   /* H parity */
    return odd;
}

/* ---------------- HP-blocking number by branch and bound ---------------- */
static u64 *hpm; static long nhp, hpcap;
static long long bb_nodes;
static void collect_hp(int *path, int len, int used) {
    if (len == N) {
        u64 m = 0;
        for (int i = 0; i + 1 < N; i++) m |= 1ULL << arcidx[path[i]][path[i + 1]];
        if (nhp == hpcap) { hpcap = hpcap ? 2 * hpcap : 1024; hpm = realloc(hpm, sizeof(u64) * hpcap); }
        hpm[nhp++] = m; return;
    }
    int last = path[len - 1];
    for (int w = 0; w < N; w++) if (!((used >> w) & 1) && A[last][w]) { path[len] = w; collect_hp(path, len + 1, used | (1 << w)); }
}
static int hits(u64 F, int depth) {
    bb_nodes++;
    long i;
    for (i = 0; i < nhp; i++) if (!(hpm[i] & F)) break;
    if (i == nhp) return 1;
    if (depth == 0) return 0;
    u64 m = hpm[i];
    while (m) { u64 b = m & (~m + 1); m ^= b; if (hits(F | b, depth - 1)) return 1; }
    return 0;
}

/* ---------------- main ---------------- */
static uint64_t rng_state;
static uint64_t rng(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage\n"); return 1; }
    const char *mode = argv[1]; N = atoi(argv[2]);
    char line[1024];
    fbuf = malloc(sizeof(u64) * ((size_t)1 << N) * N); gbuf = malloc(sizeof(u64) * ((size_t)1 << N) * N);
    if (!fbuf || !gbuf) { fprintf(stderr, "alloc\n"); return 1; }
    if (!strcmp(mode, "stats")) {
        int noaut = argc > 3 && !strcmp(argv[3], "noaut");
        while (fgets(line, sizeof line, stdin)) {
            if (line[0] != '0' && line[0] != '1') continue;
            line[strcspn(line, "\r\n")] = 0;
            set_from_string(line);
            double w = noaut ? 1.0 : factorial(N) / (double)aut_size();
            process(line, w, 1);
        }
        report();
    } else if (!strcmp(mode, "lab")) {
        int ne = N * (N - 1) / 2; char s[256];
        for (uint64_t x = 0; x < (1ULL << ne); x++) {
            for (int b = 0; b < ne; b++) s[b] = ((x >> b) & 1) ? '1' : '0';
            s[ne] = 0; set_from_string(s); process(s, 1.0, 0);
        }
        report();
    } else if (!strcmp(mode, "rnd")) {
        long cnt = atol(argv[3]); rng_state = strtoull(argv[4], 0, 10) * 2654435761ULL + 88172645463325252ULL;
        int ne = N * (N - 1) / 2; char s[256];
        for (long t = 0; t < cnt; t++) {
            for (int b = 0; b < ne; b++) s[b] = (rng() >> 33) & 1 ? '1' : '0';
            s[ne] = 0; set_from_string(s); process(s, 1.0, 1);
        }
        report();
    } else if (!strcmp(mode, "shave")) {
        int ne = N * (N - 1) / 2;
        tab = malloc(sizeof(uint32_t) << ne); good = malloc((size_t)1 << ne);
        while (fgets(line, sizeof line, stdin)) {
            if (line[0] != '0' && line[0] != '1') continue;
            line[strcspn(line, "\r\n")] = 0;
            set_from_string(line); shave_one(line);
        }
    } else if (!strcmp(mode, "anf")) {
        int ne = N * (N - 1) / 2; char s[256];
        for (uint64_t x = 0; x < (1ULL << (ne - 1)); x++) {
            s[0] = '1'; /* arc (0,1) oriented 0->1 */
            for (int b = 1; b < ne; b++) s[b] = ((x >> (b - 1)) & 1) ? '1' : '0';
            s[ne] = 0; set_from_string(s); hp_counts();
            putchar((C[0][1] & 1) ? '1' : '0');
        }
        putchar('\n');
    } else if (!strcmp(mode, "beta")) {
        int k = 0; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { arcidx[i][j] = k; arcidx[j][i] = k; k++; }
        long cls = 0, bad = 0, nhall = 0; long long hist[64] = {0};
        while (fgets(line, sizeof line, stdin)) {
            if (line[0] != '0' && line[0] != '1') continue;
            line[strcspn(line, "\r\n")] = 0;
            set_from_string(line);
            nhp = 0; int path[MAXN];
            for (int s0 = 0; s0 < N; s0++) { path[0] = s0; collect_hp(path, 1, 1 << s0); }
            int indeg[MAXN], od[MAXN], loc = 99;
            for (int v = 0; v < N; v++) { indeg[v] = 0; od[v] = 0; for (int w = 0; w < N; w++) { indeg[v] += A[w][v]; od[v] += A[v][w]; } }
            for (int u = 0; u < N; u++) for (int v = u + 1; v < N; v++) {
                int a = indeg[u] + indeg[v], b = od[u] + od[v]; if (a < loc) loc = a; if (b < loc) loc = b; }
            int hb = match_bound(); if (hb < loc) { nhall++; } if (hb < loc) loc = hb;
            cls++; hist[loc]++;
            if (loc >= 1 && hits(0, loc - 1)) { bad++; printf("BETA_BELOW_BOUND T=%s bound=%d\n", line, loc); }
        }
        printf("BETA N=%d classes=%ld hall_below_two_source=%ld beta_below_hall_bound=%ld bb_nodes=%lld bound_hist:", N, cls, nhall, bad, bb_nodes);
        for (int i = 0; i < 64; i++) if (hist[i]) printf(" %d:%lld", i, hist[i]);
        printf("\n");
    } else if (!strcmp(mode, "oddfilter")) {
        int K = atoi(argv[3]); long long hist[200] = {0}; long cls = 0, hbad = 0;
        f2 = malloc(sizeof(uint32_t) << N); g2 = malloc(sizeof(uint32_t) << N);
        while (fgets(line, sizeof line, stdin)) {
            if (line[0] != '0' && line[0] != '1') continue;
            line[strcspn(line, "\r\n")] = 0;
            set_from_string(line);
            int hp, odd = odd_arcs_mod2(&hp);
            if (!hp) hbad++;
            hist[odd]++; cls++;
            if (odd <= K) printf("FEWODD N=%d T=%s odd=%d\n", N, line, odd);
        }
        printf("ODDFILTER N=%d classes=%ld H_even_failures=%ld odd_hist:", N, cls, hbad);
        for (int i = 0; i < 200; i++) if (hist[i]) printf(" %d:%lld", i, hist[i]);
        printf("\n");
    } else if (!strcmp(mode, "list")) {
        while (fgets(line, sizeof line, stdin)) {
            if (line[0] != '0' && line[0] != '1') continue;
            line[strcspn(line, "\r\n")] = 0;
            set_from_string(line); hp_counts();
            int odd = 0, z = 0;
            for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v]) { if (C[u][v] & 1) odd++; if (!C[u][v]) z++; }
            printf("%s %llu %d %d %d %llu\n", line, Hval, odd, z, is_strong(), aut_size());
        }
    } else if (!strcmp(mode, "hlab")) {
        int ne = N * (N - 1) / 2; char s[256]; only_H = 1;
        for (uint64_t x = 0; x < (1ULL << ne); x++) {
            for (int b = 0; b < ne; b++) s[b] = ((x >> b) & 1) ? '1' : '0';
            s[ne] = 0; set_from_string(s); hp_counts();
            printf("%llu\n", Hval);
        }
    } else if (!strcmp(mode, "minodd")) {
        long iters = atol(argv[3]); rng_state = strtoull(argv[4], 0, 10) * 2654435761ULL + 7046029254386353131ULL;
        int ne = N * (N - 1) / 2, best = 1 << 30; char bests[256] = {0};
        for (long r = 0; r < iters; r++) {
            /* near-transitive start: transitive order with each arc reversed w.p. 1/8 */
            for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { int b = ((rng() >> 20) & 7) == 0; A[i][j] = !b; A[j][i] = b; }
            hp_counts();
            int cur = 0; for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v] && (C[u][v] & 1)) cur++;
            for (int step = 0; step < 4 * ne; step++) {
                int b = (int)(rng() % ne), i = 0, j = 0, kk = 0;
                for (int a = 0; a < N; a++) for (int c = a + 1; c < N; c++) { if (kk == b) { i = a; j = c; } kk++; }
                A[i][j] ^= 1; A[j][i] ^= 1;
                hp_counts();
                int nw = 0; for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v] && (C[u][v] & 1)) nw++;
                if (nw <= cur || (rng() % 16) == 0) cur = nw; else { A[i][j] ^= 1; A[j][i] ^= 1; }
                if (cur < best) {
                    best = cur; int k2 = 0;
                    for (int a = 0; a < N; a++) for (int c = a + 1; c < N; c++) bests[k2++] = A[a][c] ? '1' : '0';
                    bests[k2] = 0;
                }
            }
        }
        printf("MINODD N=%d restarts=%ld min_odd_found=%d witness=%s\n", N, iters, best, bests);
    } else if (!strcmp(mode, "search")) {
        long iters = atol(argv[3]); rng_state = strtoull(argv[4], 0, 10) * 2654435761ULL + 1442695040888963407ULL;
        int ne = N * (N - 1) / 2; char s[256]; int found = 0;
        for (long r = 0; r < iters && found < 3; r++) {
            for (int b = 0; b < ne; b++) s[b] = (rng() >> 33) & 1 ? '1' : '0';
            s[ne] = 0; set_from_string(s);
            /* greedy: minimise #even arcs by single arc flips, 4*ne steps */
            for (int step = 0; step < 6 * ne; step++) {
                hp_counts();
                int ev = 0; for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v] && !(C[u][v] & 1)) ev++;
                if (ev == 0) {
                    found++;
                    printf("FOUND_ALLODD N=%d T=", N);
                    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) putchar(A[i][j] ? '1' : '0');
                    printf(" H=%llu restart=%ld\n", Hval, r);
                    break;
                }
                /* try a random flip; accept if not worse */
                int b = (int)(rng() % ne), i = 0, j = 0, kk = 0;
                for (int a = 0; a < N; a++) for (int c = a + 1; c < N; c++) { if (kk == b) { i = a; j = c; } kk++; }
                A[i][j] ^= 1; A[j][i] ^= 1;
                hp_counts();
                int ev2 = 0; for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v] && !(C[u][v] & 1)) ev2++;
                if (ev2 > ev && (rng() % 8)) { A[i][j] ^= 1; A[j][i] ^= 1; }
            }
        }
        printf("SEARCH_DONE N=%d found=%d\n", N, found);
    }
    return 0;
}
