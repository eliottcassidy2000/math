/* Tour-blocking number of the 6x6 TORUS knight graph (chessboard-weave, 2026-10-06).
 *
 * G: vertices (i,j) in Z6 x Z6, edges to (i+-1,j+-2),(i+-2,j+-1): 8-regular, 144 edges.
 * Claim checked: for EVERY set D of K edges, G - D still has a Hamiltonian cycle
 * (or, for K = 7, list the D for which it does not).
 * By edge-transitivity (verified by the Python driver) it suffices to take D = {e0} u D'
 * with D' a (K-1)-subset of the other 143 edges: C(143,5) = 464,306,843 sets for K = 6.
 *
 * Method: a pool of Hamiltonian cycles (144-bit edge masks).  Nested loops over
 * a<b<c<d (prefix with e0) keep the sub-pool of cycles avoiding the prefix; for the
 * last edge f>d, D is certified iff some sub-pool cycle avoids f, i.e. f is not in
 * the AND of the sub-pool masks.  Uncertified D are sent to an exact DFS Hamiltonian
 * cycle search on G - D; a found cycle is verified and added to the pool; a D with
 * no cycle (or a search hitting the node limit) is printed as UNRESOLVED/BLOCKING.
 *
 * usage: ./torusblock K [pool_size] [seed] [pooldump]   (K = deletion-set size, e.g. 6 or 7;
 *        K = 7 also lists every BLOCKING 7-set containing e0 -- negative control)
 * build: cc -O2 -o torusblock chessboard_weave_20261006_knight_torusblock.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define NV 36
#define NE 144
typedef struct { uint64_t w[3]; } emask;

static int eu[NE], ev[NE];
static int eid[NV][NV];
static uint64_t nbr[NV];          /* vertex adjacency in G */
static emask *pool; static int npool, cap;

static inline int pc(uint64_t x) { return __builtin_popcountll(x); }
static inline int has(const emask *m, int e) { return (m->w[e >> 6] >> (e & 63)) & 1; }
static inline void setb(emask *m, int e) { m->w[e >> 6] |= 1ULL << (e & 63); }

/* ---- random number generator */
static uint64_t rs = 88172645463325252ULL;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

/* ---- exact / randomized Hamiltonian cycle DFS on graph with adjacency adj[] */
static uint64_t adj[NV];
static int path[NV], plen;
static long long nodes, nodelimit;
static int randomize;
static int found;

static void hc_dfs(int head, uint64_t U) {
    if (found || nodes > nodelimit) return;
    nodes++;
    if (U == 0) { if (adj[head] & 1ULL) found = 1; return; }  /* close to vertex 0 */
    /* pruning: each unvisited v needs >= 2 available among U u {head} u {0};
       vertex 0 (start) still needs one more neighbour: must be in U u {head} */
    uint64_t hb = 1ULL << head;
    int forced = -1;
    uint64_t avail0 = adj[0] & (U | hb);
    if (!avail0) return;
    for (uint64_t W = U; W; W &= W - 1) {
        int v = __builtin_ctzll(W);
        uint64_t av = adj[v] & (U | hb | 1ULL);
        int eff = pc(av);
        if (eff < 2) return;
        if (eff == 2 && (adj[v] & hb)) {
            /* v's two path-neighbours are head and one other: v must follow head,
               unless the other is vertex 0 and v is the last vertex -- still follows head */
            if (forced >= 0 && forced != v) return;
            forced = v;
        }
    }
    uint64_t cand = adj[head] & U;
    if (forced >= 0) cand = 1ULL << forced;
    int list[NV], n = 0;
    for (uint64_t W = cand; W; W &= W - 1) list[n++] = __builtin_ctzll(W);
    /* Warnsdorff order (fewest onward moves first), random tie-break */
    int key[NV];
    for (int k = 0; k < n; k++) key[k] = pc(adj[list[k]] & U) * 64 + (randomize ? (int)(rnd() & 63) : 0);
    for (int a = 1; a < n; a++) { int x = list[a], y = key[a], b = a - 1;
        while (b >= 0 && key[b] > y) { list[b + 1] = list[b]; key[b + 1] = key[b]; b--; }
        list[b + 1] = x; key[b + 1] = y; }
    for (int k = 0; k < n && !found; k++) {
        int v = list[k];
        path[plen++] = v;
        hc_dfs(v, U & ~(1ULL << v));
        if (found) return;
        plen--;
    }
}

/* returns 1 and fills *cyc if a Hamiltonian cycle of G - D exists; 0 if none; -1 if aborted */
static int find_hc(const emask *D, emask *cyc, long long limit, int rnd_order) {
    for (int v = 0; v < NV; v++) adj[v] = nbr[v];
    for (int e = 0; e < NE; e++) if (has(D, e)) { adj[eu[e]] &= ~(1ULL << ev[e]); adj[ev[e]] &= ~(1ULL << eu[e]); }
    for (int v = 0; v < NV; v++) if (pc(adj[v]) < 2) return 0;
    found = 0; nodes = 0; nodelimit = limit; randomize = rnd_order;
    plen = 0; path[plen++] = 0;
    uint64_t all = (1ULL << NV) - 1;
    hc_dfs(0, all & ~1ULL);
    if (!found) return nodes > limit ? -1 : 0;
    memset(cyc, 0, sizeof *cyc);
    for (int k = 0; k < NV; k++) {
        int a = path[k], b = path[(k + 1) % NV];
        int e = eid[a][b];
        if (e < 0 || has(D, e)) { fprintf(stderr, "BUG: bad cycle\n"); exit(2); }
        setb(cyc, e);
    }
    /* verify Hamiltonicity */
    uint64_t seen = 0; for (int k = 0; k < NV; k++) seen |= 1ULL << path[k];
    if (seen != all) { fprintf(stderr, "BUG: not spanning\n"); exit(2); }
    return 1;
}

static void pool_add(const emask *c) {
    if (npool == cap) { cap = cap ? 2 * cap : 1024; pool = realloc(pool, cap * sizeof(emask)); }
    pool[npool++] = *c;
}

/* ---- recursive enumeration of D = {e0, chosen[1..K-1]} with sub-pools */
#define LMAX 8000000
static int K = 6;
static int *L[8], n[8];
static int chosen[8];
static long long checked = 0, finder_calls = 0, unresolved = 0, blocking = 0, blocking_iso = 0;

static void handle_uncertified(int depth_full) {
    emask D; memset(&D, 0, sizeof D);
    for (int t = 0; t < K; t++) setb(&D, chosen[t]);
    emask cyc;
    finder_calls++;
    int r = find_hc(&D, &cyc, 200000000LL, 0);
    if (r == -1) r = find_hc(&D, &cyc, 200000000LL, 1);   /* retry only if aborted */
    if (r == 1) {
        pool_add(&cyc);
        int id = npool - 1;
        for (int l = 0; l < K - 1; l++) L[l][n[l]++] = id;   /* avoids every chosen edge */
        return;
    }
    /* classify: does D contain 7 edges at a single vertex (vertex left with degree 1)? */
    int iso = 0;
    for (int v = 0; v < NV; v++) { int cnt = 0; for (int t = 0; t < K; t++) if (eu[chosen[t]] == v || ev[chosen[t]] == v) cnt++; if (cnt >= 7) iso = 1; }
    if (r == 0) { blocking++; blocking_iso += iso; } else unresolved++;
    printf("%s%s D =", r == 0 ? "BLOCKING" : "UNRESOLVED", iso ? " (vertex-isolating)" : "");
    for (int t = 0; t < K; t++) printf(" {%d,%d}", eu[chosen[t]], ev[chosen[t]]);
    printf("\n");
    (void)depth_full;
}

static void rec(int level, int from) {
    /* chosen[0..level-1] fixed; L[level-1] = pool cycles avoiding them */
    if (level == K - 1) {
        /* last edge f: D certified iff some sub-pool cycle avoids f, i.e. f not in AND */
        emask I; I.w[0] = I.w[1] = I.w[2] = ~0ULL;
        for (int k = 0; k < n[level - 1]; k++) { const emask *c = &pool[L[level - 1][k]]; I.w[0] &= c->w[0]; I.w[1] &= c->w[1]; I.w[2] &= c->w[2]; }
        for (int f = from; f < NE; f++) {
            checked++;
            if (n[level - 1] > 0 && !has(&I, f)) continue;
            chosen[level] = f;
            int before = npool;
            handle_uncertified(level);
            if (npool > before) { const emask *c = &pool[npool - 1]; I.w[0] &= c->w[0]; I.w[1] &= c->w[1]; I.w[2] &= c->w[2]; }
        }
        return;
    }
    for (int e = from; e < NE; e++) {
        chosen[level] = e;
        n[level] = 0;
        for (int k = 0; k < n[level - 1]; k++) if (!has(&pool[L[level - 1][k]], e)) L[level][n[level]++] = L[level - 1][k];
        rec(level + 1, e + 1);
        if (level == 1) {
            printf("  first free edge %3d done: sets checked %lld, finder calls %lld, pool %d, blocking %lld (iso %lld), unresolved %lld\n",
                   e, checked, finder_calls, npool, blocking, blocking_iso, unresolved);
            fflush(stdout);
        }
    }
}

int main(int argc, char **argv) {
    K = argc > 1 ? atoi(argv[1]) : 6;
    int P = argc > 2 ? atoi(argv[2]) : 1000;
    if (argc > 3) rs ^= (uint64_t)atoll(argv[3]) * 0x9E3779B97F4A7C15ULL;
    const char *dumpf = argc > 4 ? argv[4] : NULL;
    if (K < 2 || K > 7) { fprintf(stderr, "K must be 2..7\n"); return 1; }
    memset(eid, -1, sizeof eid);
    int ne = 0;
    int d[8][2] = {{1,2},{2,1},{-1,2},{-2,1},{1,-2},{2,-1},{-1,-2},{-2,-1}};
    for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) {
        int u = 6 * i + j;
        for (int k = 0; k < 8; k++) {
            int v = 6 * ((i + d[k][0] + 6) % 6) + (j + d[k][1] + 6) % 6;
            nbr[u] |= 1ULL << v;
            if (u < v && eid[u][v] < 0) { eid[u][v] = eid[v][u] = ne; eu[ne] = u; ev[ne] = v; ne++; }
        }
    }
    for (int v = 0; v < NV; v++) if (pc(nbr[v]) != 8) { fprintf(stderr, "not 8-regular\n"); return 1; }
    printf("6x6 torus knight graph: %d vertices, %d edges (8-regular). e0 = {%d,%d}\n", NV, ne, eu[0], ev[0]);
    /* initial random pool */
    emask none; memset(&none, 0, sizeof none);
    for (int t = 0; t < P; t++) {
        emask c;
        int r = find_hc(&none, &c, 100000000LL, 1);
        if (r != 1) { fprintf(stderr, "pool generation failed\n"); return 1; }
        pool_add(&c);
    }
    printf("initial pool: %d random Hamiltonian cycles\n", npool);
    fflush(stdout);
    if (dumpf) {
        FILE *fp = fopen(dumpf, "w");
        for (int k = 0; k < npool; k++) fprintf(fp, "%016llx %016llx %016llx\n",
            (unsigned long long)pool[k].w[0], (unsigned long long)pool[k].w[1], (unsigned long long)pool[k].w[2]);
        fclose(fp);
    }
    /* enumerate D = {e0} u (K-1)-subset of the other edges, recursively */
    for (int l = 0; l <= K; l++) L[l] = malloc(sizeof(int) * LMAX);
    n[0] = 0; for (int k = 0; k < npool; k++) if (!has(&pool[k], 0)) L[0][n[0]++] = k;
    chosen[0] = 0;
    rec(1, 1);
    printf("TOTAL %d-sets containing e0 checked: %lld\n", K, checked);
    printf("finder calls %lld, final pool %d, BLOCKING sets %lld (vertex-isolating type %lld), UNRESOLVED %lld\n",
           finder_calls, npool, blocking, blocking_iso, unresolved);
    if (dumpf) {
        FILE *fp = fopen(dumpf, "w");
        for (int k = 0; k < npool; k++) fprintf(fp, "%016llx %016llx %016llx\n",
            (unsigned long long)pool[k].w[0], (unsigned long long)pool[k].w[1], (unsigned long long)pool[k].w[2]);
        fclose(fp);
    }
    return 0;
}
