/*
 * chessboard_weave_20261006_hp_blocking_engine.c
 *
 * HYP-9168 (beta(T) = hall(T)) engine.
 *
 *   beta(T) = fewest arcs whose deletion leaves no Hamiltonian path (HP).
 *   hall(T) = min over X (|X|>=2), Y (|Y|<=|X|-2) of #arcs u->x, x in X, u not in Y
 *             (or the dual with out-arcs).
 *   KONIG form (checked here against the brute force):
 *   hall(T) = min over A,B subsets of V with |A|+|B| = N+2 of e(A,B),
 *             e(A,B) = #{arcs u->v of T : u in A, v in B}.
 *   Fast evaluation: for each A, the best B of size N+2-|A| takes the vertices of
 *   smallest in-degree from A.
 *
 * Input: one tournament per line, upper-triangle string (gentourng format:
 *   char k for pair (i,j), i<j in row order, '1' means i->j), optionally followed by
 *   automorphisms "A:p0,p1,...,p_{N-1}" (used only for root symmetry breaking in bnb).
 *
 * Modes (argv[1]):
 *   hall      : hall via (A,B) formula, via brute (X,Y) in- and out-versions, sigma
 *   brute     : beta by exhaustive subset enumeration (small N) + hall
 *   bnb       : hall, then exact test "is there a blocking set of size hall-1?" by
 *               branch and bound (HP-guided branching with protection, packing LB);
 *               if yes, exact beta by decreasing k. Optional argv[2] = node limit.
 *   minsets   : beta by brute force, then enumerate ALL minimum blocking sets and
 *               classify: cluster-type (deleted pairs = disjoint union of cliques),
 *               has 1-path-cycle factor (bipartite matching of size >= N-1).
 *   both      : brute and bnb, compare.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define MAXN 16
typedef uint32_t mask;

static int N;
static mask T_out[MAXN], T_in[MAXN];
static mask G_out[MAXN], G_in[MAXN];
static mask PR[MAXN];               /* protected out-arcs: PR[u] bit v <=> arc u->v protected */
static mask BASE_PR[MAXN];          /* permanently protected arcs (kpath mode: arcs at universal vertices) */
static uint16_t ends_buf[1 << MAXN];
static int naut;
static int aut[256][MAXN];

static long long nodes, node_limit = -1;
static long long q_count, q_nodes, q_maxnodes, q_hist[64], q_hall_lt_sigma, q_bad, q_abort, q_counter;
static char q_maxtour[512];
static int aborted;

static inline int pc(mask x) { return __builtin_popcount(x); }
static inline int ctz(mask x) { return __builtin_ctz(x); }

/* ---------------- parsing ---------------- */
static int parse_line(char *line) {
    char *tok = strtok(line, " \t\r\n");
    if (!tok) return -1;
    int len = (int)strlen(tok);
    int n = 1;
    while (n * (n - 1) / 2 < len) n++;
    if (n * (n - 1) / 2 != len || n > MAXN) return -1;
    N = n;
    memset(T_out, 0, sizeof T_out);
    memset(T_in, 0, sizeof T_in);
    int k = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++, k++) {
            if (tok[k] == '1') { T_out[i] |= 1u << j; T_in[j] |= 1u << i; }
            else if (tok[k] == '2') {   /* semicomplete extension: both arcs (2-cycle) */
                T_out[i] |= 1u << j; T_in[j] |= 1u << i;
                T_out[j] |= 1u << i; T_in[i] |= 1u << j;
            }
            else { T_out[j] |= 1u << i; T_in[i] |= 1u << j; }
        }
    naut = 0;
    while ((tok = strtok(NULL, " \t\r\n"))) {
        if (tok[0] == 'D' && tok[1] == ':') {
            /* pre-deleted arcs "D:u>v,u>v" -> T becomes a general oriented graph (validation only) */
            char *p = tok + 2;
            while (*p) {
                int u = (int)strtol(p, &p, 10);
                if (*p == '>') p++;
                int v = (int)strtol(p, &p, 10);
                T_out[u] &= ~(1u << v); T_in[v] &= ~(1u << u);
                if (*p == ',') p++; else break;
            }
            continue;
        }
        if (tok[0] == 'A' && tok[1] == ':' && naut < 256) {
            char *p = tok + 2;
            for (int i = 0; i < N; i++) {
                aut[naut][i] = (int)strtol(p, &p, 10);
                if (*p == ',') p++;
            }
            naut++;
        }
    }
    return N;
}

/* ---------------- HP existence: exact subset DP ---------------- */
/* ends[S] = set of v in S such that G[S] has an HP ending at v */
static int hp_exact(const mask *out, const mask *in, int *path) {
    (void)out;
    uint32_t full = (N == 32) ? 0xffffffffu : ((1u << N) - 1);
    uint16_t *ends = ends_buf;
    memset(ends, 0, sizeof(uint16_t) << N);
    for (int v = 0; v < N; v++) ends[1u << v] = (uint16_t)(1u << v);
    for (uint32_t S = 1; S < full; S++) {
        uint16_t e = ends[S];
        if (!e) continue;
        mask cand = full & ~S;
        while (cand) {
            int w = ctz(cand);
            cand &= cand - 1;
            if (e & in[w]) ends[S | (1u << w)] |= (uint16_t)(1u << w);
        }
    }
    if (!ends[full]) return 0;
    if (path) {
        uint32_t S = full;
        int v = ctz(ends[full]);
        path[N - 1] = v;
        for (int pos = N - 2; pos >= 0; pos--) {
            uint32_t S2 = S & ~(1u << v);
            mask cand = ends[S2] & in[v];
            int u = ctz(cand);
            path[pos] = u;
            S = S2;
            v = u;
        }
    }
    return 1;
}

/* verify a path is an HP of (out) */
static int check_path(const mask *out, const int *path) {
    mask seen = 0;
    for (int i = 0; i < N; i++) {
        if (path[i] < 0 || path[i] >= N || (seen >> path[i] & 1)) return 0;
        seen |= 1u << path[i];
        if (i + 1 < N && !((out[path[i]] >> path[i + 1]) & 1)) return 0;
    }
    return 1;
}

/* ---------------- HP heuristic: Redei insertion with protected preference ---------------- */
static uint64_t rng_state = 88172645463325252ULL;
static inline uint64_t rng(void) {
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return rng_state;
}

static int try_insert(const mask *out, int *path, int *len, int w, const mask *prot) {
    int L = *len;
    if (L == 0) { path[0] = w; *len = 1; return 1; }
    int best = -100, bpos = -1;
    /* prepend */
    if ((out[w] >> path[0]) & 1) {
        int sc = prot ? ((prot[w] >> path[0]) & 1) : 0;
        if (sc > best) { best = sc; bpos = 0; }
    }
    if ((out[path[L - 1]] >> w) & 1) {
        int sc = prot ? ((prot[path[L - 1]] >> w) & 1) : 0;
        if (sc > best) { best = sc; bpos = L; }
    }
    for (int i = 1; i < L; i++) {
        int a = path[i - 1], b = path[i];
        if (((out[a] >> w) & 1) && ((out[w] >> b) & 1)) {
            int sc = 0;
            if (prot) sc = (int)((prot[a] >> w) & 1) + (int)((prot[w] >> b) & 1) - (int)((prot[a] >> b) & 1);
            if (sc > best) { best = sc; bpos = i; }
        }
    }
    if (bpos < 0) return 0;
    for (int i = L; i > bpos; i--) path[i] = path[i - 1];
    path[bpos] = w;
    *len = L + 1;
    return 1;
}

static int hp_heur(const mask *out, int *path, const mask *prot, int tries) {
    int order[MAXN];
    for (int t = 0; t < tries; t++) {
        for (int i = 0; i < N; i++) order[i] = i;
        for (int i = N - 1; i > 0; i--) { int j = (int)(rng() % (uint64_t)(i + 1)); int x = order[i]; order[i] = order[j]; order[j] = x; }
        int len = 0, pend[MAXN], np = 0;
        for (int i = 0; i < N; i++) if (!try_insert(out, path, &len, order[i], prot)) pend[np++] = order[i];
        int progress = 1;
        while (np && progress) {
            progress = 0;
            for (int i = 0; i < np; i++) {
                if (try_insert(out, path, &len, pend[i], prot)) { pend[i] = pend[--np]; i--; progress = 1; }
            }
        }
        if (np == 0) return 1;
    }
    return 0;
}

static int find_hp(const mask *out, const mask *in, int *path, const mask *prot, int exact) {
    if (hp_heur(out, path, prot, 3)) return 1;
    if (!exact) return 0;
    return hp_exact(out, in, path);
}

/* ---------------- hall / sigma ---------------- */
static int hall_extra = 2; /* hall_k uses |A|+|B| = N + k + 1; default k = 1 */
static int hall_formula(mask *bestA, mask *bestB) {
    uint32_t full = (1u << N) - 1;
    int best = 1 << 30;
    for (uint32_t A = 0; A <= full; A++) {
        int a = pc(A);
        int b = N + hall_extra - a; /* |A| + |B| = N + hall_extra */
        if (b < 0 || b > N) continue;
        int d[MAXN];
        int cnt[MAXN + 1];
        memset(cnt, 0, sizeof cnt);
        for (int v = 0; v < N; v++) { d[v] = pc(T_in[v] & A); cnt[d[v]]++; }
        int s = 0, need = b;
        for (int x = 0; x <= N && need > 0; x++) {
            int take = cnt[x] < need ? cnt[x] : need;
            s += take * x;
            need -= take;
        }
        if (s < best) {
            best = s;
            if (bestA) {
                *bestA = A;
                /* reconstruct B */
                int idx[MAXN];
                for (int v = 0; v < N; v++) idx[v] = v;
                for (int i = 0; i < N; i++)
                    for (int j = i + 1; j < N; j++)
                        if (d[idx[j]] < d[idx[i]]) { int t = idx[i]; idx[i] = idx[j]; idx[j] = t; }
                mask B = 0;
                for (int i = 0; i < b; i++) B |= 1u << idx[i];
                *bestB = B;
            }
        }
    }
    return best;
}

static void hall_brute(int *hin, int *hout) {
    uint32_t full = (1u << N) - 1;
    int bi = 1 << 30, bo = 1 << 30;
    for (uint32_t X = 0; X <= full; X++) {
        int x = pc(X);
        if (x < 2) continue;
        for (uint32_t Y = 0; Y <= full; Y++) {
            if (pc(Y) > x - 2) continue;
            int si = 0, so = 0;
            for (int v = 0; v < N; v++)
                if ((X >> v) & 1) { si += pc(T_in[v] & ~Y & full); so += pc(T_out[v] & ~Y & full); }
            if (si < bi) bi = si;
            if (so < bo) bo = so;
        }
    }
    *hin = bi; *hout = bo;
}

static int sigma_bound(void) {
    int best = 1 << 30;
    for (int u = 0; u < N; u++)
        for (int v = u + 1; v < N; v++) {
            int a = pc(T_in[u]) + pc(T_in[v]);
            int b = pc(T_out[u]) + pc(T_out[v]);
            if (a < best) best = a;
            if (b < best) best = b;
        }
    return best;
}

/* ---------------- matching (1-path-cycle factor) ---------------- */
static int mt_right[MAXN];
static mask mt_vis;
static int kuhn(const mask *out, int u) {
    mask c = out[u];
    while (c) {
        int v = ctz(c); c &= c - 1;
        if ((mt_vis >> v) & 1) continue;
        mt_vis |= 1u << v;
        if (mt_right[v] < 0 || kuhn(out, mt_right[v])) { mt_right[v] = u; return 1; }
    }
    return 0;
}
static int max_matching(const mask *out) {
    for (int v = 0; v < N; v++) mt_right[v] = -1;
    int m = 0;
    for (int u = 0; u < N; u++) { mt_vis = 0; m += kuhn(out, u); }
    return m;
}

/* ---------------- arcs list ---------------- */
static int na;
static int arc_u[MAXN * MAXN], arc_v[MAXN * MAXN];
static void build_arcs(void) {
    na = 0;
    for (int u = 0; u < N; u++)
        for (int v = 0; v < N; v++)
            if ((T_out[u] >> v) & 1) { arc_u[na] = u; arc_v[na] = v; na++; }
}

/* ---------------- brute-force beta ---------------- */
/* returns beta; enumerates subsets of arcs of increasing size (combinations) */
static int comb_idx[64];
static int next_comb(int k, int n) {
    int i = k - 1;
    while (i >= 0 && comb_idx[i] == n - k + i) i--;
    if (i < 0) return 0;
    comb_idx[i]++;
    for (int j = i + 1; j < k; j++) comb_idx[j] = comb_idx[j - 1] + 1;
    return 1;
}
static void apply_set(const int *idx, int k) {
    memcpy(G_out, T_out, sizeof(mask) * N);
    memcpy(G_in, T_in, sizeof(mask) * N);
    for (int i = 0; i < k; i++) {
        int u = arc_u[idx[i]], v = arc_v[idx[i]];
        G_out[u] &= ~(1u << v); G_in[v] &= ~(1u << u);
    }
}
static int beta_brute(int maxk) {
    build_arcs();
    for (int k = 0; k <= maxk; k++) {
        for (int i = 0; i < k; i++) comb_idx[i] = i;
        do {
            apply_set(comb_idx, k);
            if (!hp_exact(G_out, G_in, NULL)) return k;
        } while (k > 0 && next_comb(k, na));
    }
    return -1;
}

/* cluster test for a set of deleted pairs given as adjacency masks (undirected) */
static int is_cluster(const mask *und) {
    for (int v = 0; v < N; v++) {
        mask cl = und[v] | (1u << v);
        mask nb = und[v];
        while (nb) {
            int u = ctz(nb); nb &= nb - 1;
            if ((und[u] | (1u << u)) != cl) return 0;
        }
    }
    return 1;
}

/* ---------------- branch and bound ---------------- */
static int best_witness_k;
static int wit_u[MAXN * MAXN], wit_v[MAXN * MAXN], wit_n;
static int cur_u[MAXN * MAXN], cur_v[MAXN * MAXN], cur_n;

/* packing lower bound; also returns an HP of G with fewest unprotected arcs found.
   returns -1 if G has no HP at all (blocking achieved), else t = #HPs found with
   pairwise disjoint unprotected-arc sets (capped at cap); t = 1<<20 if some HP is
   fully protected. */
static int packing(int cap, int *bpath, int *bcnt) {
    mask o[MAXN], in[MAXN];
    memcpy(o, G_out, sizeof(mask) * N);
    memcpy(in, G_in, sizeof(mask) * N);
    int t = 0;
    *bcnt = 1 << 20;
    int path[MAXN];
    while (t < cap) {
        int ok = find_hp(o, in, path, PR, 1);
        if (!ok) {
            if (t == 0) return -1;
            break;
        }
        if (!check_path(o, path)) { fprintf(stderr, "FATAL: invalid HP returned\n"); exit(2); }
        int c = 0;
        for (int i = 0; i + 1 < N; i++) if (!((PR[path[i]] >> path[i + 1]) & 1)) c++;
        if (c == 0) { memcpy(bpath, path, sizeof(int) * N); *bcnt = 0; return 1 << 20; }
        if (c < *bcnt) { *bcnt = c; memcpy(bpath, path, sizeof(int) * N); }
        t++;
        for (int i = 0; i + 1 < N; i++) {
            int a = path[i], b = path[i + 1];
            if (!((PR[a] >> b) & 1)) { o[a] &= ~(1u << b); in[b] &= ~(1u << a); }
        }
    }
    return t;
}

static int bnb(int b) {
    nodes++;
    if (node_limit > 0 && nodes > node_limit) { aborted = 1; return 0; }
    int path[MAXN], cnt;
    int t = packing(b + 1, path, &cnt);
    if (t < 0) {
        /* G is HP-free */
        wit_n = cur_n;
        memcpy(wit_u, cur_u, sizeof(int) * cur_n);
        memcpy(wit_v, cur_v, sizeof(int) * cur_n);
        return 1;
    }
    if (b == 0) return 0;
    if (t > b) return 0;
    int ua[MAXN], ub[MAXN], m = 0;
    for (int i = 0; i + 1 < N; i++) {
        int a = path[i], c = path[i + 1];
        if (!((PR[a] >> c) & 1)) { ua[m] = a; ub[m] = c; m++; }
    }
    int res = 0;
    int i;
    for (i = 0; i < m; i++) {
        int a = ua[i], c = ub[i];
        G_out[a] &= ~(1u << c); G_in[c] &= ~(1u << a);
        cur_u[cur_n] = a; cur_v[cur_n] = c; cur_n++;
        int r = bnb(b - 1);
        cur_n--;
        G_out[a] |= 1u << c; G_in[c] |= 1u << a;
        if (r) { res = 1; i++; break; }
        if (aborted) { i++; break; }
        PR[a] |= 1u << c;
    }
    for (int j = 0; j < i && j < m; j++) PR[ua[j]] &= ~(1u << ub[j]);
    /* note: arcs j < i were protected (the last one only if loop completed); clearing all is safe
       because they were all unprotected on entry */
    return res;
}

/* enumerate ALL blocking sets of size exactly all_target (each exactly once, by the
   protection scheme); for each, record whether T - D still has a 1-path-cycle factor
   and whether D is the star of a vertex */
static int all_target;
static long long all_count, all_exotic, all_exotic_star, all_star;
static int max_matching(const mask *out);
static void bnb_all(int b) {
    nodes++;
    int path[MAXN], cnt;
    int t = packing(b + 1, path, &cnt);
    if (t < 0) {
        if (cur_n == all_target) {
            all_count++;
            int mm = max_matching(G_out);
            int deg[MAXN];
            memset(deg, 0, sizeof deg);
            for (int i = 0; i < cur_n; i++) { deg[cur_u[i]]++; deg[cur_v[i]]++; }
            int star = 0;  /* D = all N-1 arcs at one vertex (isolates it) */
            for (int v = 0; v < N; v++) if (deg[v] == cur_n && cur_n == N - 1) star = 1;
            all_star += star;
            if (mm >= N - 1) {
                all_exotic++;
                all_exotic_star += star;
                if (!star) {
                    printf("EXOTIC_NONSTAR D=");
                    for (int i = 0; i < cur_n; i++) printf("%d>%d%s", cur_u[i], cur_v[i], i + 1 < cur_n ? "," : "");
                    printf(" ");
                }
            }
        }
        return;
    }
    if (b == 0 || t > b) return;
    int ua[MAXN], ub[MAXN], m = 0;
    for (int i = 0; i + 1 < N; i++) {
        int a = path[i], c = path[i + 1];
        if (!((PR[a] >> c) & 1)) { ua[m] = a; ub[m] = c; m++; }
    }
    for (int i = 0; i < m; i++) {
        int a = ua[i], c = ub[i];
        G_out[a] &= ~(1u << c); G_in[c] &= ~(1u << a);
        cur_u[cur_n] = a; cur_v[cur_n] = c; cur_n++;
        bnb_all(b - 1);
        cur_n--;
        G_out[a] |= 1u << c; G_in[c] |= 1u << a;
        PR[a] |= 1u << c;
    }
    for (int j = 0; j < m; j++) PR[ua[j]] &= ~(1u << ub[j]);
}

/* root with symmetry breaking via arc orbits under the supplied automorphisms */
static int arc_orbit[MAXN][MAXN];
static int feasible(int k) {
    /* is there a blocking set of size <= k ? */
    memcpy(G_out, T_out, sizeof(mask) * N);
    memcpy(G_in, T_in, sizeof(mask) * N);
    memcpy(PR, BASE_PR, sizeof PR);
    cur_n = 0;
    aborted = 0;
    if (!hp_exact(G_out, G_in, NULL)) { wit_n = 0; return 1; }
    if (k == 0) return 0;
    if (naut <= 1) return bnb(k);
    /* orbits of arcs */
    int norb = 0;
    for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) arc_orbit[u][v] = -1;
    for (int u = 0; u < N; u++)
        for (int v = 0; v < N; v++) {
            if (!((T_out[u] >> v) & 1) || arc_orbit[u][v] >= 0) continue;
            /* BFS closure under all automorphisms supplied (we assume the full group is supplied) */
            for (int g = 0; g < naut; g++) arc_orbit[aut[g][u]][aut[g][v]] = norb;
            norb++;
        }
    for (int o = 0; o < norb; o++) {
        /* representative */
        int ru = -1, rv = -1;
        for (int u = 0; u < N && ru < 0; u++) for (int v = 0; v < N; v++) if (arc_orbit[u][v] == o) { ru = u; rv = v; break; }
        memcpy(PR, BASE_PR, sizeof PR);
        for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (arc_orbit[u][v] >= 0 && arc_orbit[u][v] < o) PR[u] |= 1u << v;
        G_out[ru] &= ~(1u << rv); G_in[rv] &= ~(1u << ru);
        cur_u[0] = ru; cur_v[0] = rv; cur_n = 1;
        int r = bnb(k - 1);
        cur_n = 0;
        G_out[ru] |= 1u << rv; G_in[rv] |= 1u << ru;
        if (r) return 1;
        if (aborted) return 0;
    }
    return 0;
}

static int verify_witness(void) {
    memcpy(G_out, T_out, sizeof(mask) * N);
    memcpy(G_in, T_in, sizeof(mask) * N);
    for (int i = 0; i < wit_n; i++) {
        if (!((T_out[wit_u[i]] >> wit_v[i]) & 1)) return 0;
        G_out[wit_u[i]] &= ~(1u << wit_v[i]); G_in[wit_v[i]] &= ~(1u << wit_u[i]);
    }
    return !hp_exact(G_out, G_in, NULL);
}

/* check that hall obstruction (A,B) really kills all HPs (sanity) */
static int check_hall_kills(mask A, mask B) {
    memcpy(G_out, T_out, sizeof(mask) * N);
    memcpy(G_in, T_in, sizeof(mask) * N);
    int cnt = 0;
    for (int u = 0; u < N; u++) {
        if (!((A >> u) & 1)) continue;
        mask c = T_out[u] & B;
        cnt += pc(c);
        G_out[u] &= ~c;
        while (c) { int v = ctz(c); c &= c - 1; G_in[v] &= ~(1u << u); }
    }
    int has = hp_exact(G_out, G_in, NULL);
    int mm = max_matching(G_out);
    return (!has && mm <= N - 2) ? cnt : -1;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s mode [node_limit]\n", argv[0]); return 1; }
    const char *mode = argv[1];
    if (argc >= 3) node_limit = atoll(argv[2]);
    char line[8192];
    long long lineno = 0;
    while (fgets(line, sizeof line, stdin)) {
        char keep[8192];
        strcpy(keep, line);
        char *nl = strchr(keep, '\n'); if (nl) *nl = 0;
        char *sp = strchr(keep, ' '); if (sp) *sp = 0;
        if (line[0] == '>' || line[0] == '#') continue;
        if (parse_line(line) < 0) continue;
        lineno++;
        mask bA = 0, bB = 0;
        int h = hall_formula(&bA, &bB);
        int sg = sigma_bound();
        if (!strcmp(mode, "hall")) {
            int hi, ho;
            hall_brute(&hi, &ho);
            int hk = check_hall_kills(bA, bB);
            printf("%s N=%d hall_formula=%d hall_in=%d hall_out=%d sigma=%d obstruction_ok=%d %s\n", keep, N, h, hi, ho, sg,
                   hk == h, (h == hi && h == ho) ? "AGREE" : "DISAGREE");
        } else if (!strcmp(mode, "brute")) {
            int bt = beta_brute(N);
            printf("%s N=%d hall=%d sigma=%d beta_brute=%d %s\n", keep, N, h, sg, bt, bt == h ? "EQ" : "NEQ");
        } else if (!strcmp(mode, "bnbq")) {
            /* quiet exhaustive mode: print only non-EQ lines; summary at the end */
            nodes = 0;
            int hk = check_hall_kills(bA, bB);
            int f = feasible(h - 1);
            q_count++;
            q_nodes += nodes;
            if (nodes > q_maxnodes) { q_maxnodes = nodes; strcpy(q_maxtour, keep); }
            if (h <= MAXN * 2) q_hist[h]++;
            if (h < sg) q_hall_lt_sigma++;
            if (hk != h) { q_bad++; printf("OBSTRUCTION_FAIL %s hall=%d\n", keep, h); }
            if (aborted) { q_abort++; printf("ABORTED %s hall=%d nodes=%lld\n", keep, h, nodes); }
            else if (f) {
                q_counter++;
                printf("COUNTEREXAMPLE %s hall=%d sigma=%d witness_size=%d witness=", keep, h, sg, wit_n);
                for (int i = 0; i < wit_n; i++) printf("%d>%d%s", wit_u[i], wit_v[i], i + 1 < wit_n ? "," : "");
                printf(" verified=%d\n", verify_witness());
            }
            fflush(stdout);
        } else if (!strcmp(mode, "bnb") || !strcmp(mode, "both")) {
            clock_t t0 = clock();
            nodes = 0;
            int hk = check_hall_kills(bA, bB);
            int f = feasible(h - 1);
            int beta = h;
            int status_ok = 1;
            if (aborted) { beta = -1; status_ok = 0; }
            else if (f) {
                /* counterexample: find exact beta */
                if (!verify_witness()) status_ok = 0;
                beta = wit_n;
                while (beta > 0) {
                    if (feasible(beta - 1)) { beta = wit_n; if (!verify_witness()) status_ok = 0; }
                    else break;
                }
            }
            double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
            printf("%s N=%d hall=%d sigma=%d beta=%d nodes=%lld time=%.3f obstruction_ok=%d %s", keep, N, h, sg, beta, nodes, secs,
                   hk == h, aborted ? "ABORTED" : (beta == h ? "EQ" : "COUNTEREXAMPLE"));
            if (!status_ok && !aborted) printf(" WITNESS_FAIL");
            if (!strcmp(mode, "both")) {
                int bt = beta_brute(N);
                printf(" beta_brute=%d %s", bt, bt == beta ? "MATCH" : "MISMATCH");
            }
            if (f && !aborted) {
                printf(" witness=");
                for (int i = 0; i < wit_n; i++) printf("%d>%d%s", wit_u[i], wit_v[i], i + 1 < wit_n ? "," : "");
            }
            printf("\n");
            fflush(stdout);
        } else if (!strcmp(mode, "deep") || !strcmp(mode, "deepsym")) {
            /* beta by iterative deepening with the branch and bound, vs brute force.
               "deep": no symmetry breaking; "deepsym": root symmetry breaking with the supplied
               automorphism group (validates the orbit code). */
            nodes = 0;
            int saved_naut = naut;
            if (!strcmp(mode, "deep")) naut = 0;
            int k = 0;
            while (!feasible(k)) k++;
            int okw = verify_witness() && wit_n == k;
            naut = saved_naut;
            int bt = beta_brute(N);
            printf("%s N=%d hall=%d beta_deep=%d beta_brute=%d witness_ok=%d %s\n", keep, N, h, k, bt, okw, (k == bt && okw) ? "MATCH" : "MISMATCH");
            fflush(stdout);
        } else if (!strcmp(mode, "exotic")) {
            /* list every minimum blocking set D such that T - D still has a 1-path-cycle factor */
            int bt = beta_brute(N);
            int k = bt;
            for (int i = 0; i < k; i++) comb_idx[i] = i;
            do {
                apply_set(comb_idx, k);
                if (hp_exact(G_out, G_in, NULL)) continue;
                int mm = max_matching(G_out);
                if (mm < N - 1) continue;
                printf("%s N=%d beta=%d hall=%d D=", keep, N, bt, h);
                for (int i = 0; i < k; i++) printf("%d>%d%s", arc_u[comb_idx[i]], arc_v[comb_idx[i]], i + 1 < k ? "," : "");
                printf(" matching=");
                for (int v = 0; v < N; v++) if (mt_right[v] >= 0) printf("%d>%d,", mt_right[v], v);
                printf("\n");
            } while (k > 0 && next_comb(k, na));
            fflush(stdout);
        } else if (!strcmp(mode, "gutin")) {
            /* all set partitions of V: D = arcs inside parts; check HP <=> 1-path-cycle factor
               (Gutin 1993, Corollary 1, for the multipartite tournament T - D) */
            int rgs[MAXN];
            for (int i = 0; i < N; i++) rgs[i] = 0;
            long long parts = 0, mism = 0;
            int minblock = 1 << 30;
            int bt = beta_brute(N);
            while (1) {
                memcpy(G_out, T_out, sizeof(mask) * N);
                memcpy(G_in, T_in, sizeof(mask) * N);
                int dcount = 0;
                for (int u = 0; u < N; u++)
                    for (int v = 0; v < N; v++)
                        if (u != v && rgs[u] == rgs[v] && ((T_out[u] >> v) & 1)) {
                            G_out[u] &= ~(1u << v); G_in[v] &= ~(1u << u); dcount++;
                        }
                int hp = hp_exact(G_out, G_in, NULL);
                int pcf = max_matching(G_out) >= N - 1;
                parts++;
                if (hp != pcf) mism++;
                if (!hp && dcount < minblock) minblock = dcount;
                /* next restricted growth string */
                int i = N - 1;
                while (i > 0) {
                    int mx = 0;
                    for (int j = 0; j < i; j++) if (rgs[j] > mx) mx = rgs[j];
                    if (rgs[i] <= mx) { rgs[i]++; for (int j = i + 1; j < N; j++) rgs[j] = 0; break; }
                    i--;
                }
                if (i == 0) break;
            }
            printf("%s N=%d hall=%d beta=%d partitions=%lld hp_vs_1pcf_mismatch=%lld min_cluster_blocking=%d\n",
                   keep, N, h, bt, parts, mism, minblock);
            fflush(stdout);
        } else if (!strcmp(mode, "kpath")) {
            /* beta_K = fewest deletions leaving no cover by <= K paths, vs hall_K = min_{|A|+|B|=N+K+1} e(A,B).
               K-1 universal vertices (2-cycles to everything, permanently protected) turn
               "cover by <= K paths" into "Hamiltonian path". argv[2] = K. */
            int K = (argc >= 4) ? atoi(argv[3]) : 2;
            hall_extra = K + 1;
            int hK = hall_formula(NULL, NULL);
            hall_extra = 2;
            build_arcs();          /* original arcs only */
            int N0 = N;
            int N1 = N0 + K - 1;
            for (int u = N0; u < N1; u++) { T_out[u] = 0; T_in[u] = 0; }
            memset(BASE_PR, 0, sizeof BASE_PR);
            for (int u = N0; u < N1; u++)
                for (int v = 0; v < N1; v++)
                    if (v != u) {
                        T_out[u] |= 1u << v; T_in[v] |= 1u << u;
                        T_out[v] |= 1u << u; T_in[u] |= 1u << v;
                        BASE_PR[u] |= 1u << v; BASE_PR[v] |= 1u << u;
                    }
            N = N1;
            int kk = 0;
            nodes = 0;
            while (!feasible(kk)) kk++;
            int bbrute = -1;
            if (argc >= 5 && atoi(argv[4]) == 1) {
                /* brute force over subsets of original arcs */
                for (int k2 = 0; k2 <= na && bbrute < 0; k2++) {
                    for (int i = 0; i < k2; i++) comb_idx[i] = i;
                    do {
                        apply_set(comb_idx, k2);
                        if (!hp_exact(G_out, G_in, NULL)) { bbrute = k2; break; }
                    } while (k2 > 0 && next_comb(k2, na));
                }
            }
            memset(BASE_PR, 0, sizeof BASE_PR);
            printf("%s N=%d K=%d hall_K=%d beta_K=%d beta_K_brute=%d %s\n", keep, N0, K, hK, kk, bbrute,
                   kk == hK ? "EQ" : "NEQ");
            fflush(stdout);
        } else if (!strcmp(mode, "allmin")) {
            /* all minimum blocking sets via the enumerating branch and bound (beta = hall assumed
               only as the target size; verified by the absence of smaller sets in bnbq runs) */
            memcpy(G_out, T_out, sizeof(mask) * N);
            memcpy(G_in, T_in, sizeof(mask) * N);
            memset(PR, 0, sizeof PR);
            cur_n = 0; nodes = 0;
            all_target = h; all_count = all_exotic = all_exotic_star = all_star = 0;
            printf("%s ", keep);
            bnb_all(h);
            printf("N=%d hall=%d min_sets=%lld isolating_stars=%lld with_1pcf=%lld with_1pcf_isolating_star=%lld nodes=%lld\n",
                   N, h, all_count, all_star, all_exotic, all_exotic_star, nodes);
            fflush(stdout);
        } else if (!strcmp(mode, "minsets")) {
            int bt = beta_brute(N);
            /* enumerate all blocking sets of size bt */
            long long nmin = 0, ncl = 0, npcf = 0, ncl_pcf = 0, ncf = 0;
            int k = bt;
            for (int i = 0; i < k; i++) comb_idx[i] = i;
            do {
                apply_set(comb_idx, k);
                if (hp_exact(G_out, G_in, NULL)) continue;
                nmin++;
                mask und[MAXN];
                memset(und, 0, sizeof und);
                for (int i = 0; i < k; i++) { int u = arc_u[comb_idx[i]], v = arc_v[comb_idx[i]]; und[u] |= 1u << v; und[v] |= 1u << u; }
                int cl = is_cluster(und);
                int mm = max_matching(G_out);
                int pcf = (mm >= N - 1);
                if (mm >= N) ncf++;
                ncl += cl; npcf += pcf; ncl_pcf += (cl && pcf);
            } while (k > 0 && next_comb(k, na));
            printf("%s N=%d hall=%d sigma=%d beta=%d nmin=%lld cluster=%lld with_1pcf=%lld cluster_and_1pcf=%lld with_cycle_factor=%lld\n",
                   keep, N, h, sg, bt, nmin, ncl, npcf, ncl_pcf, ncf);
            fflush(stdout);
        }
    }
    if (!strcmp(mode, "bnbq")) {
        printf("SUMMARY tournaments=%lld total_nodes=%lld max_nodes=%lld (%s) hall_lt_sigma=%lld obstruction_fail=%lld aborted=%lld counterexamples=%lld\n",
               q_count, q_nodes, q_maxnodes, q_maxtour, q_hall_lt_sigma, q_bad, q_abort, q_counter);
        printf("HALL_HISTOGRAM");
        for (int i = 0; i < 64; i++) if (q_hist[i]) printf(" %d:%lld", i, q_hist[i]);
        printf("\n");
    }
    return 0;
}
