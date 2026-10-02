// Shaved tournaments, n = 9 and the audit of n <= 8 (thread session, 2026-10-01).
// Helper for shaved_tournaments_n9_20261001.py.  Build: gcc -O2 -fopenmp -o shv9 shaved_tournaments_n9_20261001.c
//
// A shaving of order n is a spanning oriented graph S on n vertices contained in every n-tournament.  Shavings are
// acyclic (they lie in the transitive tournament) and hereditary (spanning subgraphs of shavings are shavings).
//
// Modes (tournaments are read as gentourng ascii: the upper triangle row by row, '1' = i->j for i<j):
//   orderly  n host classes e out [killfrom] [seedk]
//       all e-arc acyclic spanning subgraphs of the fixed host tournament, up to Aut(host), that embed in every
//       tournament of the class list.  Every n-shaving embeds in the host, so the output decides u(n) >= e.
//       Orderly generation: a set is canonical if it is the lexicographically least image of itself under
//       Aut(host) (arcs of the host indexed 0..NP-1); deleting the largest arc of a canonical set leaves a canonical
//       set, and acyclicity and containment in a fixed tournament pass to subsets, so a DFS that only extends
//       canonical, acyclic, killer-surviving sets by larger arcs reaches every canonical solution.  Killers are
//       tournaments that avoided an earlier candidate (move-to-front); they only prune, the final test is a full
//       pass over the class list.  Output: one arc mask per line (bit b = b-th arc of the host).
//   count    n classes candidates
//       exhaustive embedding counts (plain DFS, no look-ahead): for each candidate (line "u1 v1 u2 v2 ..."),
//       the minimum number of embeddings over all classes and the number of classes with an even count.
//   bruteext n classes candidates
//       brute force over all n! bijections (no pruning): prints every candidate contained in every class.
//   blockers n classes candidates [auts]
//       which classes avoid which candidates: the most frequent blockers and a greedy blocking set.  With auts
//       (one |Aut| per class, in class order), also the candidates that no symmetric class (|Aut| > 1) avoids
//       (if there is one, every blocking set contains a rigid class), the candidates that no rigid class avoids
//       (if there is one, every blocking set contains a symmetric class), and the candidate with fewest avoiders.
//   thma     n            (stdin: class list)
//       counts tournaments with no Hamiltonian path whose first vertex beats its last one (no copy of H_n).
//   lrctight k B
//       primitive tight lonely-runner instances {n_1 < ... < n_k}, gcd 1, n_i <= B, lon = 1/(k+1) exactly, and
//       whether each is the connection set of a circulant tournament on Z_(2k+1).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXN 12
typedef struct { uint16_t out[MAXN]; uint16_t in[MAXN]; } tour_t;
typedef struct { uint16_t pred[MAXN]; uint8_t nsucc[MAXN]; } dag_t;   // forward-labelled: pred[k] = {j<k : j->k}

static int NV, NP, PI[66], PJ[66];
static void init_pairs(int n) {
    NV = n; NP = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { PI[NP] = i; PJ[NP] = j; NP++; }
}
static int load_tours(const char *fn, tour_t **out) {
    FILE *f = fopen(fn, "r"); if (!f) { perror(fn); exit(1); }
    int cap = 1 << 16, cnt = 0; tour_t *T = malloc(sizeof(tour_t) * cap); char buf[256];
    while (fgets(buf, sizeof buf, f)) {
        int L = strlen(buf); while (L && (buf[L-1] == '\n' || buf[L-1] == '\r')) buf[--L] = 0;
        if (L != NP) continue;
        if (cnt == cap) { cap *= 2; T = realloc(T, sizeof(tour_t) * cap); }
        memset(&T[cnt], 0, sizeof(tour_t));
        for (int b = 0; b < NP; b++) {
            int i = PI[b], j = PJ[b];
            if (buf[b] == '1') { T[cnt].out[i] |= 1u << j; T[cnt].in[j] |= 1u << i; }
            else { T[cnt].out[j] |= 1u << i; T[cnt].in[i] |= 1u << j; }
        }
        cnt++;
    }
    fclose(f); *out = T; return cnt;
}
// does the forward-labelled S embed in T?  look-ahead: the image of k needs >= nsucc[k] unused out-neighbours
static int embed(const dag_t *d, const tour_t *t) {
    int n = NV; uint8_t img[MAXN]; uint16_t cand[MAXN]; uint16_t used = 0; int k = 0;
    uint16_t full = (uint16_t)((1u << n) - 1), c2 = 0;
    for (int y = 0; y < n; y++) if (__builtin_popcount(t->out[y]) >= d->nsucc[0]) c2 |= 1u << y;
    cand[0] = c2;
    for (;;) {
        if (!cand[k]) { if (k == 0) return 0; k--; used &= ~(1u << img[k]); continue; }
        int y = __builtin_ctz(cand[k]); cand[k] &= cand[k] - 1;
        img[k] = y; used |= 1u << y;
        if (k == n - 1) return 1;
        k++;
        uint16_t c = full & ~used, p = d->pred[k], fr = full & ~used;
        while (p) { int j = __builtin_ctz(p); p &= p - 1; c &= t->out[img[j]]; }
        uint16_t c3 = 0;
        while (c) { int z = __builtin_ctz(c); c &= c - 1;
            if (__builtin_popcount(t->out[z] & fr & ~(1u << z)) >= d->nsucc[k]) c3 |= 1u << z; }
        cand[k] = c3;
    }
}
// exhaustive count of embeddings (no look-ahead)
static long count_emb(const dag_t *d, const tour_t *t) {
    int n = NV; uint8_t img[MAXN]; uint16_t cand[MAXN]; uint16_t used = 0; int k = 0; long cnt = 0;
    uint16_t full = (uint16_t)((1u << n) - 1);
    cand[0] = full;
    for (;;) {
        if (!cand[k]) { if (k == 0) return cnt; k--; used &= ~(1u << img[k]); continue; }
        int y = __builtin_ctz(cand[k]); cand[k] &= cand[k] - 1;
        img[k] = y; used |= 1u << y;
        if (k == n - 1) { cnt++; used &= ~(1u << y); continue; }
        k++;
        uint16_t c = full & ~used, p = d->pred[k];
        while (p) { int j = __builtin_ctz(p); p &= p - 1; c &= t->out[img[j]]; }
        cand[k] = c;
    }
}
// arc list -> forward-labelled dag via a topological order; returns 0 if cyclic
static int dag_from_arcs(int e, const int *U, const int *V, dag_t *d) {
    int n = NV; uint16_t inm[MAXN] = {0};
    for (int t = 0; t < e; t++) inm[V[t]] |= 1u << U[t];
    int pos[MAXN]; uint16_t placed = 0;
    for (int k = 0; k < n; k++) {
        int v; for (v = 0; v < n; v++) if (!(placed >> v & 1) && (inm[v] & ~placed) == 0) break;
        if (v == n) return 0;
        pos[v] = k; placed |= 1u << v;
    }
    memset(d, 0, sizeof *d);
    for (int t = 0; t < e; t++) { d->pred[pos[V[t]]] |= 1u << pos[U[t]]; d->nsucc[pos[U[t]]]++; }
    return 1;
}
static int read_arcs(char *s, int *U, int *V) {
    int e = 0, a, b, k;
    while (sscanf(s, "%d %d%n", &a, &b, &k) == 2) { U[e] = a; V[e] = b; e++; s += k; }
    return e;
}
static int next_perm(int *a, int n) {
    int i = n - 2; while (i >= 0 && a[i] >= a[i + 1]) i--; if (i < 0) return 0;
    int j = n - 1; while (a[j] <= a[i]) j--; int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

// ------------------------------------------------------------------------------------------------ orderly
static tour_t HOST; static int HA_u[66], HA_v[66]; static int NG; static uint8_t *GP; static int TGT, KILL_FROM;
static tour_t *CL; static int NC;
static inline uint64_t img_mask(uint64_t m, const uint8_t *g) {
    uint64_t r = 0; while (m) { int b = __builtin_ctzll(m); m &= m - 1; r |= 1ull << g[b]; } return r;
}
static inline int lex_less(uint64_t a, uint64_t b) { uint64_t x = a ^ b; return x && (a & x & -x); }
static int canonical(uint64_t m) {
    for (int g = 1; g < NG; g++) if (lex_less(img_mask(m, GP + (size_t)g * NP), m)) return 0;
    return 1;
}
static int mask_to_dag(uint64_t m, dag_t *d) {
    int U[66], V[66], e = 0;
    for (uint64_t x = m; x; x &= x - 1) { int b = __builtin_ctzll(x); U[e] = HA_u[b]; V[e] = HA_v[b]; e++; }
    return dag_from_arcs(e, U, V, d);
}
#define MAXK 4096
typedef struct { int nk; int k[MAXK]; int *ord; long full_checks, sols; FILE *fo; } ctx_t;
static int alive_killers(const dag_t *d, ctx_t *c) {
    for (int q = 0; q < c->nk; q++)
        if (!embed(d, &CL[c->k[q]])) {
            if (q) { int t = c->k[q]; memmove(c->k + 1, c->k, q * sizeof(int)); c->k[0] = t; }
            return 0;
        }
    return 1;
}
static int full_check(const dag_t *d, ctx_t *c) {
    c->full_checks++;
    for (int q = 0; q < NC; q++) {
        int ci = c->ord[q];
        if (!embed(d, &CL[ci])) {
            if (q) { memmove(c->ord + 1, c->ord, q * sizeof(int)); c->ord[0] = ci; }
            if (c->nk < MAXK) { memmove(c->k + 1, c->k, c->nk * sizeof(int)); c->k[0] = ci; c->nk++; }
            return 0;
        }
    }
    return 1;
}
static void emit(uint64_t m, ctx_t *c) {
    c->sols++;
    if (c->fo) {
        #pragma omp critical
        { fprintf(c->fo, "%llu\n", (unsigned long long)m); fflush(c->fo); }
    }
}
static void dfs(uint64_t m, int size, int maxb, ctx_t *c) {
    for (int b = maxb + 1; b < NP; b++) {
        if (NP - b < TGT - size) break;
        uint64_t m2 = m | (1ull << b); dag_t d;
        if (!mask_to_dag(m2, &d)) continue;
        if (!canonical(m2)) continue;
        if (size + 1 >= KILL_FROM && !alive_killers(&d, c)) continue;
        if (size + 1 == TGT) { if (full_check(&d, c)) emit(m2, c); }
        else dfs(m2, size + 1, b, c);
    }
}
static uint64_t *SEEDS; static int NS, SEEDCAP, SEEDK;
static void seeds_rec(uint64_t m, int size, int maxb) {
    if (size == SEEDK) {
        if (NS == SEEDCAP) { SEEDCAP *= 2; SEEDS = realloc(SEEDS, sizeof(uint64_t) * SEEDCAP); }
        SEEDS[NS++] = m; return;
    }
    for (int b = maxb + 1; b < NP; b++) {
        if (NP - b < TGT - size) break;
        uint64_t m2 = m | (1ull << b); dag_t d;
        if (!mask_to_dag(m2, &d) || !canonical(m2)) continue;
        seeds_rec(m2, size + 1, b);
    }
}
static int mode_orderly(int argc, char **argv) {
    int n = atoi(argv[2]); init_pairs(n);
    tour_t *H; if (load_tours(argv[3], &H) != 1) { fprintf(stderr, "host file must hold one tournament\n"); return 1; }
    HOST = H[0]; NC = load_tours(argv[4], &CL); TGT = atoi(argv[5]);
    FILE *fo = fopen(argv[6], "w");
    KILL_FROM = argc > 7 ? atoi(argv[7]) : TGT - 4;
    SEEDK = argc > 8 ? atoi(argv[8]) : 4; if (SEEDK > TGT) SEEDK = TGT;
    int na = 0, arcid[MAXN][MAXN];
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (HOST.out[u] >> v & 1) { arcid[u][v] = na; HA_u[na] = u; HA_v[na] = v; na++; }
    int p[MAXN]; for (int i = 0; i < n; i++) p[i] = i;
    int cap = 1024; GP = malloc((size_t)cap * NP); NG = 0;
    do {   // Aut(host) by brute force; the identity comes first
        int ok = 1;
        for (int a = 0; a < na && ok; a++) if (!(HOST.out[p[HA_u[a]]] >> p[HA_v[a]] & 1)) ok = 0;
        if (ok) { if (NG == cap) { cap *= 2; GP = realloc(GP, (size_t)cap * NP); }
            for (int a = 0; a < na; a++) GP[(size_t)NG * NP + a] = arcid[p[HA_u[a]]][p[HA_v[a]]];
            NG++; }
    } while (next_perm(p, n));
    SEEDCAP = 1 << 16; SEEDS = malloc(sizeof(uint64_t) * SEEDCAP); NS = 0;
    seeds_rec(0, 0, -1);
    long full = 0, sols = 0;
    #pragma omp parallel reduction(+:full,sols)
    {
        ctx_t *c = malloc(sizeof(ctx_t)); c->nk = 0; c->full_checks = 0; c->sols = 0; c->fo = fo;
        c->ord = malloc(sizeof(int) * NC); for (int i = 0; i < NC; i++) c->ord[i] = i;
        #pragma omp for schedule(dynamic, 1)
        for (int s = 0; s < NS; s++) {
            uint64_t m = SEEDS[s]; int sz = __builtin_popcountll(m), maxb = m ? 63 - __builtin_clzll(m) : -1;
            if (sz == TGT) { dag_t d; mask_to_dag(m, &d);
                if ((sz < KILL_FROM || alive_killers(&d, c)) && full_check(&d, c)) emit(m, c);
                continue; }
            dfs(m, sz, maxb, c);
        }
        full += c->full_checks; sols += c->sols;
    }
    printf("orderly n=%d e=%d host|Aut|=%d classes=%d: canonical e-arc shavings inside the host (Aut-orbits) = %ld (full checks %ld)\n",
           n, TGT, NG, NC, sols, full);
    fclose(fo);
    return 0;
}

// ------------------------------------------------------------------------------------------------ count
static int load_candidates(const char *fn, dag_t **out) {
    FILE *g = fopen(fn, "r"); if (!g) { perror(fn); exit(1); }
    char buf[2048]; int cap = 4096, nc = 0; dag_t *D = malloc(sizeof(dag_t) * cap);
    while (fgets(buf, sizeof buf, g)) {
        int U[66], V[66], e = read_arcs(buf, U, V);
        if (!e) continue;
        if (nc == cap) { cap *= 2; D = realloc(D, sizeof(dag_t) * cap); }
        if (!dag_from_arcs(e, U, V, &D[nc])) { fprintf(stderr, "cyclic candidate\n"); exit(1); }
        nc++;
    }
    fclose(g); *out = D; return nc;
}
static int mode_count(char **argv) {
    int n = atoi(argv[2]); init_pairs(n);
    tour_t *T; int nt = load_tours(argv[3], &T);
    dag_t *D; int nc = load_candidates(argv[4], &D);
    long *even = calloc(nc, sizeof(long)), *zero = calloc(nc, sizeof(long)), *mn = malloc(sizeof(long) * nc);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int c = 0; c < nc; c++) {
        long e = 0, z = 0, m = 1L << 60;
        for (int i = 0; i < nt; i++) { long x = count_emb(&D[c], &T[i]); if (!(x & 1)) e++; if (!x) z++; if (x < m) m = x; }
        even[c] = e; zero[c] = z; mn[c] = m;
    }
    int pos = 0, allodd = 0; long gmin = 1L << 60;
    for (int c = 0; c < nc; c++) { if (mn[c] > 0) pos++; if (!even[c]) allodd++; if (mn[c] < gmin) gmin = mn[c]; }
    printf("count n=%d classes=%d candidates=%d: embedded in every class = %d (least minimum %ld); odd count in every class = %d\n",
           n, nt, nc, pos, gmin, allodd);
    printf("  histogram of the minimum embedding count over the classes:");
    for (long v = 0; v <= 64; v++) { int h = 0; for (int c = 0; c < nc; c++) if (mn[c] == v) h++; if (h) printf(" %ld:%d", v, h); }
    printf("\n");
    for (int c = 0; c < nc; c++) if (!even[c]) printf("  odd in every class: candidate %d\n", c);
    if (nc <= 5) for (int c = 0; c < nc; c++)
        printf("  candidate %d: minimum %ld, classes with an even count %ld, classes avoiding it %ld\n", c, mn[c], even[c], zero[c]);
    return 0;
}

// ------------------------------------------------------------------------------------------------ bruteext
static int mode_bruteext(char **argv) {
    int n = atoi(argv[2]); init_pairs(n);
    tour_t *T; int nt = load_tours(argv[3], &T);
    int *ord = malloc(sizeof(int) * nt); for (int i = 0; i < nt; i++) ord[i] = i;
    FILE *g = fopen(argv[4], "r"); if (!g) { perror(argv[4]); return 1; }
    char buf[2048]; long cands = 0, all = 0;
    while (fgets(buf, sizeof buf, g)) {
        int U[66], V[66], e = read_arcs(buf, U, V);
        if (!e) continue;
        cands++;
        int killed = 0;
        for (int q = 0; q < nt && !killed; q++) {
            int ci = ord[q], p[MAXN], found = 0;
            for (int i = 0; i < n; i++) p[i] = i;
            do { int good = 1;
                for (int t = 0; t < e; t++) if (!(T[ci].out[p[U[t]]] >> p[V[t]] & 1)) { good = 0; break; }
                if (good) { found = 1; break; }
            } while (next_perm(p, n));
            if (!found) { killed = 1; if (q) { memmove(ord + 1, ord, q * sizeof(int)); ord[0] = ci; } }
        }
        if (!killed) { all++; printf("  contained in every class: %s", buf); }
    }
    printf("bruteext n=%d classes=%d candidates=%ld: contained in every class = %ld\n", n, nt, cands, all);
    return 0;
}

// ------------------------------------------------------------------------------------------------ blockers
// which classes avoid which candidates: the most frequent blockers and a greedy blocking set (a certificate that
// none of the candidates is a shaving, if every candidate is avoided by some class)
static int mode_blockers(int argc, char **argv) {
    int n = atoi(argv[2]); init_pairs(n);
    tour_t *T; int nt = load_tours(argv[3], &T);
    dag_t *D; int nc = load_candidates(argv[4], &D);
    int *aut = NULL;
    if (argc == 6) {
        FILE *f = fopen(argv[5], "r"); if (!f) { perror(argv[5]); exit(1); }
        aut = malloc(nt * sizeof(int));
        for (int i = 0; i < nt; i++) if (fscanf(f, "%d", &aut[i]) != 1) { fprintf(stderr, "short aut file\n"); exit(1); }
        fclose(f);
    }
    int words = (nc + 63) / 64;
    uint64_t *M = calloc((size_t)nt * words, sizeof(uint64_t));     // M[i] = candidates avoided by class i
    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < nt; i++)
        for (int c = 0; c < nc; c++) if (!embed(&D[c], &T[i])) M[(size_t)i * words + c / 64] |= 1ull << (c % 64);
    int *tot = calloc(nt, sizeof(int)); long unblocked = 0;
    for (int i = 0; i < nt; i++) for (int w = 0; w < words; w++) tot[i] += __builtin_popcountll(M[(size_t)i * words + w]);
    for (int c = 0; c < nc; c++) { int any = 0; for (int i = 0; i < nt && !any; i++) if (M[(size_t)i * words + c / 64] >> (c % 64) & 1) any = 1; if (!any) unblocked++; }
    printf("blockers n=%d classes=%d candidates=%d: candidates avoided by no class = %ld\n", n, nt, nc, unblocked);
    if (aut) {   // per candidate: how many classes avoid it, and whether a symmetric or a rigid one does
        long onlyrigid = 0, onlysym = 0; int fewest = -1, fewc = -1;
        for (int c = 0; c < nc; c++) {
            int k = 0, sym = 0, rig = 0;
            for (int i = 0; i < nt; i++)
                if (M[(size_t)i * words + c / 64] >> (c % 64) & 1) { k++; if (aut[i] > 1) sym = 1; else rig = 1; }
            if (k > 0 && !sym) onlyrigid++;
            if (k > 0 && !rig) onlysym++;
            if (k > 0 && (fewest < 0 || k < fewest)) { fewest = k; fewc = c; }
        }
        printf("  candidates avoided only by rigid classes = %ld\n", onlyrigid);
        printf("  candidates avoided only by symmetric classes = %ld\n", onlysym);
        printf("  fewest classes avoiding one candidate = %d (candidate %d)", fewest, fewc);
        if (fewest > 0 && fewest <= 5) {
            printf(", avoided by class");
            for (int i = 0; i < nt; i++) if (M[(size_t)i * words + fewc / 64] >> (fewc % 64) & 1) printf(" %d", i);
        }
        printf("\n");
    }
    for (int r = 0; r < 15; r++) {   // the 15 most frequent blockers (class index, number of candidates avoided)
        int best = -1; for (int i = 0; i < nt; i++) if (tot[i] >= 0 && (best < 0 || tot[i] > tot[best])) best = i;
        if (best < 0 || tot[best] <= 0) break;
        printf("  top blocker %d: class %d avoids %d\n", r + 1, best, tot[best]); tot[best] = -1 - tot[best];
    }
    uint64_t *left = malloc(words * sizeof(uint64_t)); int remaining = nc - (int)unblocked;
    for (int w = 0; w < words; w++) left[w] = 0;
    for (int c = 0; c < nc; c++) { int any = 0; for (int i = 0; i < nt && !any; i++) if (M[(size_t)i * words + c / 64] >> (c % 64) & 1) any = 1; if (any) left[c / 64] |= 1ull << (c % 64); }
    int steps = 0;
    while (remaining > 0) {
        int best = -1, bc = 0;
        for (int i = 0; i < nt; i++) { int k = 0; for (int w = 0; w < words; w++) k += __builtin_popcountll(M[(size_t)i * words + w] & left[w]); if (k > bc) { bc = k; best = i; } }
        for (int w = 0; w < words; w++) left[w] &= ~M[(size_t)best * words + w];
        remaining -= bc; steps++;
        printf("  greedy %d: class %d blocks %d more (left %d)\n", steps, best, bc, remaining);
    }
    printf("greedy blocking set size = %d\n", steps);
    return 0;
}

// ------------------------------------------------------------------------------------------------ thma
static uint16_t ST[1 << MAXN][MAXN];
static int mode_thma(char **argv) {
    int n = atoi(argv[2]); int full = (1 << n) - 1; char buf[256]; long cnt = 0, bad = 0;
    while (fgets(buf, sizeof buf, stdin)) {
        int L = strlen(buf); while (L && (buf[L-1] == '\n' || buf[L-1] == '\r')) buf[--L] = 0;
        if (L != n * (n - 1) / 2) continue;
        uint16_t out[MAXN] = {0}, in[MAXN] = {0}; int b = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, b++) {
            if (buf[b] == '1') { out[i] |= 1 << j; in[j] |= 1 << i; } else { out[j] |= 1 << i; in[i] |= 1 << j; }
        }
        cnt++;
        // ST[m][e] = set of start vertices of Hamiltonian paths of T[m] ending at e
        for (int m = 1; m <= full; m++) for (int e = 0; e < n; e++) ST[m][e] = 0;
        for (int v = 0; v < n; v++) ST[1 << v][v] = 1 << v;
        int found = 0;
        for (int m = 1; m <= full && !found; m++)
            for (int e = 0; e < n; e++) {
                uint16_t s = ST[m][e]; if (!s) continue;
                if (m == full) { if (s & in[e]) { found = 1; break; } continue; }   // some start beats the end
                uint16_t o = out[e] & ~m;
                while (o) { int w = __builtin_ctz(o); o &= o - 1; ST[m | 1 << w][w] |= s; }
            }
        if (!found) { bad++; printf("  no copy of H_%d: %s\n", n, buf); }
    }
    printf("thma n=%d tournaments=%ld without a copy of H_n (path + first->last arc) = %ld\n", n, cnt, bad);
    return 0;
}

// ------------------------------------------------------------------------------------------------ lrctight
static int LK, LB, LS[16]; static long LT, LC, LX, LN;
static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static int fval(int p, int q) {   // q * min_i ||(p/q) n_i||
    int best = q;
    for (int i = 0; i < LK; i++) { long r = ((long)p * LS[i]) % q; int d = r < q - r ? (int)r : (int)(q - r); if (d < best) best = d; }
    return best;
}
static void lon(int *num, int *den) {   // max over the critical times a/(n_i + n_j), a/|n_i - n_j|, (2a+1)/(2 n_i)
    int bn = 0, bd = 1;
    #define CONSIDER(P, Q) do { int v_ = fval((P), (Q)); if ((long)v_ * bd > (long)bn * (Q)) { bn = v_; bd = (Q); } } while (0)
    for (int i = 0; i < LK; i++) {
        for (int a = 0; a < LS[i]; a++) CONSIDER(2 * a + 1, 2 * LS[i]);
        for (int j = i + 1; j < LK; j++) {
            int s = LS[i] + LS[j], d = LS[j] - LS[i];
            for (int a = 1; a < s; a++) CONSIDER(a, s);
            for (int a = 1; a < d; a++) CONSIDER(a, d);
        }
    }
    *num = bn; *den = bd;
}
static int circulant_connection_set(void) {
    int m = 2 * LK + 1, seen[64] = {0};
    for (int i = 0; i < LK; i++) {
        int r = LS[i] % m; if (r == 0) return 0;
        int s = (m - r) % m; if (seen[r] || seen[s]) return 0;
        seen[r] = seen[s] = 1;
    }
    return 1;
}
static void lrc_rec(int idx, int start) {
    if (idx == LK) {
        int g = 0; for (int i = 0; i < LK; i++) g = gcd(g, LS[i]);
        if (g != 1) return;
        LN++;
        int num, den; lon(&num, &den);
        long lhs = (long)num * (LK + 1);
        if (lhs < den) { LX++; printf("  lon < 1/(k+1):"); for (int i = 0; i < LK; i++) printf(" %d", LS[i]); printf("\n"); return; }
        if (lhs == den) {
            int c = circulant_connection_set(); LT++; LC += c;
            printf("  tight:"); for (int i = 0; i < LK; i++) printf(" %d", LS[i]);
            printf("   circulant tournament connection set mod %d: %s\n", 2 * LK + 1, c ? "yes" : "no");
        }
        return;
    }
    for (int v = start; v <= LB; v++) { LS[idx] = v; lrc_rec(idx + 1, v + 1); }
}
static int mode_lrctight(char **argv) {
    LK = atoi(argv[2]); LB = atoi(argv[3]); LT = LC = LX = LN = 0;
    lrc_rec(0, 1);
    printf("lrctight k=%d B=%d primitive sets=%ld tight=%ld (connection sets mod %d: %ld) below 1/(k+1)=%ld\n",
           LK, LB, LN, LT, 2 * LK + 1, LC, LX);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 7 && !strcmp(argv[1], "orderly")) return mode_orderly(argc, argv);
    if (argc == 5 && !strcmp(argv[1], "count")) return mode_count(argv);
    if (argc == 5 && !strcmp(argv[1], "bruteext")) return mode_bruteext(argv);
    if ((argc == 5 || argc == 6) && !strcmp(argv[1], "blockers")) return mode_blockers(argc, argv);
    if (argc == 3 && !strcmp(argv[1], "thma")) return mode_thma(argv);
    if (argc == 4 && !strcmp(argv[1], "lrctight")) return mode_lrctight(argv);
    fprintf(stderr, "usage: see the header comment\n");
    return 2;
}
