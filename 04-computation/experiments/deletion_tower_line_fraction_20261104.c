/*
 * deletion_tower_line_fraction_20261104.c  (opus-2026-11-04-S72)
 *
 * Companion to deletion_tower_census_20261104.c (same canonical form and class generation, copied verbatim).
 *
 * (1) The line fraction of a merged node.  For a fixed-path tiling t of order n let d(t) delete a path endpoint and
 *     kappa complement all free (non-path) arcs.  THM-796 (2.5): d kappa = kappa d.  For a merged node u put
 *
 *         Phi(u) = sum over tilings t in the fibre of u of ( e_[kappa d t] - e_[d t] )   in Z^(M_(n-1)),
 *
 *     i.e. the fraction prod_t p_[kappa d t] / p_[d t] when the order n-1 nodes are read as primes.  Computed from all
 *     2^binom(n-1,2) tilings for n <= 8, through both endpoint faces (terminal and initial vertex deleted).
 *     Checks: sum_u Phi(u) = 0, and the two faces agree.
 *
 * (2) Every collision of the unweighted converse-merged deck (multiset of [T - x] over all vertices x) for n <= NMAX,
 *     printed with the decks and the path-weighted decks W of the colliding nodes.
 *
 * Usage: deletion_tower_line_fraction NMAX   (NMAX <= 9; Phi is computed for n <= min(NMAX, 8))
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 9
typedef unsigned long long u64;

/* ---------------- canonical form ---------------- */

static int g_n;
static const uint16_t *g_out;
static int g_cellof[MAXN];      /* cell index of position */
static int g_cellstart[MAXN + 1], g_ncell;
static int g_cellverts[MAXN][MAXN], g_cellsize[MAXN];
static int g_used[MAXN], g_pos[MAXN];
static u64 g_best;
static int g_have, g_cnt, g_bestperm[MAXN];

static void canon_rec(int q, u64 code) {
    int n = g_n;
    if (q == n) {
        if (!g_have || code < g_best) { g_best = code; g_have = 1; g_cnt = 1; memcpy(g_bestperm, g_pos, sizeof(int) * n); }
        else if (code == g_best) g_cnt++;
        return;
    }
    int c = g_cellof[q];
    for (int k = 0; k < g_cellsize[c]; k++) {
        int v = g_cellverts[c][k];
        if (g_used[v]) continue;
        u64 nc = code;
        for (int p = 0; p < q; p++) nc = (nc << 1) | ((g_out[g_pos[p]] >> v) & 1u);
        if (g_have) {
            /* compare prefix: bits so far = q(q+1)/2 */
            int nb = (q + 1) * q / 2, tot = n * (n - 1) / 2;
            u64 bestpre = g_best >> (tot - nb);
            if (nc > bestpre) continue;
        }
        g_used[v] = 1; g_pos[q] = v;
        canon_rec(q + 1, nc);
        g_used[v] = 0;
    }
}

static int popc16(unsigned x) { return __builtin_popcount(x); }

/* returns canonical code; *aut = |Aut|; perm (optional) receives best order */
static u64 canon(int n, const uint16_t *out, int *aut, int *perm) {
    int col[MAXN], ncol;
    /* initial colours: (score, c3(v)) */
    int key[MAXN];
    for (int v = 0; v < n; v++) {
        int s = popc16(out[v]);
        unsigned inv = (~out[v]) & ((1u << n) - 1) & ~(1u << v);
        int c3 = 0;
        for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) c3 += popc16(out[u] & inv);
        key[v] = s * 1000 + c3;
    }
    /* rank keys */
    {
        int tmp[MAXN]; memcpy(tmp, key, sizeof(int) * n);
        for (int i = 1; i < n; i++) { int x = tmp[i], j = i - 1; while (j >= 0 && tmp[j] > x) { tmp[j + 1] = tmp[j]; j--; } tmp[j + 1] = x; }
        int d[MAXN], nd = 0;
        for (int i = 0; i < n; i++) if (i == 0 || tmp[i] != tmp[i - 1]) d[nd++] = tmp[i];
        for (int v = 0; v < n; v++) for (int i = 0; i < nd; i++) if (d[i] == key[v]) { col[v] = i; break; }
        ncol = nd;
    }
    /* refinement */
    for (;;) {
        /* signature: col[v], then counts of out-neighbours per colour */
        int sig[MAXN][MAXN + 1];
        for (int v = 0; v < n; v++) {
            sig[v][0] = col[v];
            for (int c = 0; c < ncol; c++) sig[v][c + 1] = 0;
            for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) sig[v][col[u] + 1]++;
        }
        int ord[MAXN]; for (int i = 0; i < n; i++) ord[i] = i;
        for (int i = 1; i < n; i++) {
            int x = ord[i], j = i - 1;
            while (j >= 0 && memcmp(sig[ord[j]], sig[x], sizeof(int) * (ncol + 1)) > 0) { ord[j + 1] = ord[j]; j--; }
            ord[j + 1] = x;
        }
        int ncol2 = 0, newcol[MAXN];
        for (int i = 0; i < n; i++) {
            if (i > 0 && memcmp(sig[ord[i]], sig[ord[i - 1]], sizeof(int) * (ncol + 1)) != 0) ncol2++;
            newcol[ord[i]] = ncol2;
        }
        ncol2++;
        memcpy(col, newcol, sizeof(int) * n);
        if (ncol2 == ncol) break;
        ncol = ncol2;
    }
    /* cells */
    g_ncell = ncol;
    for (int c = 0; c < ncol; c++) g_cellsize[c] = 0;
    for (int v = 0; v < n; v++) { int c = col[v]; g_cellverts[c][g_cellsize[c]++] = v; }
    int p = 0;
    for (int c = 0; c < ncol; c++) { g_cellstart[c] = p; for (int k = 0; k < g_cellsize[c]; k++) g_cellof[p++] = c; }
    g_n = n; g_out = out; g_have = 0; g_cnt = 0;
    memset(g_used, 0, sizeof(g_used));
    canon_rec(0, 0);
    if (aut) *aut = g_cnt;
    if (perm) memcpy(perm, g_bestperm, sizeof(int) * n);
    return g_best;
}

static void decode(int n, u64 code, uint16_t *out) {
    int tot = n * (n - 1) / 2, b = tot - 1;
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int q = 1; q < n; q++) for (int p = 0; p < q; p++) {
        int bit = (code >> b) & 1; b--;
        if (bit) out[p] |= 1u << q; else out[q] |= 1u << p;
    }
}

static void converse(int n, const uint16_t *out, uint16_t *op) {
    for (int i = 0; i < n; i++) op[i] = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if ((out[i] >> j) & 1) op[j] |= 1u << i;
}

static void delete_vertex(int n, const uint16_t *out, int x, uint16_t *res) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (i == x) continue;
        uint16_t r = 0; int kk = 0;
        for (int j = 0; j < n; j++) { if (j == x) continue; if ((out[i] >> j) & 1) r |= 1u << kk; kk++; }
        res[k++] = r;
    }
}

/* Hamiltonian path counts: hs[v] (paths starting at v), he[v] (ending at v) */
static u64 dp[1 << MAXN][MAXN];
static void ham_counts(int n, const uint16_t *out, u64 *hs, u64 *he) {
    int full = (1 << n) - 1;
    /* ending counts: dp[S][v] = paths covering S ending at v */
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
        u64 c = dp[S][v]; if (!c) continue;
        unsigned nb = out[v] & ~S;
        while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[S | (1 << w)][w] += c; }
    }
    for (int v = 0; v < n; v++) he[v] = dp[full][v];
    /* starting counts: run on converse */
    uint16_t op[MAXN]; converse(n, out, op);
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
        u64 c = dp[S][v]; if (!c) continue;
        unsigned nb = op[v] & ~S;
        while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[S | (1 << w)][w] += c; }
    }
    for (int v = 0; v < n; v++) hs[v] = dp[full][v];
}

/* ---------------- hash set of codes ---------------- */
typedef struct { u64 *key; int *val; size_t cap; size_t cnt; } HMap;
static void hm_init(HMap *h, size_t cap) { h->cap = cap; h->cnt = 0; h->key = malloc(sizeof(u64) * cap); h->val = malloc(sizeof(int) * cap); for (size_t i = 0; i < cap; i++) h->key[i] = ~0ULL; }
static size_t hm_slot(HMap *h, u64 k) { u64 x = k * 0x9E3779B97F4A7C15ULL; size_t i = (size_t)(x >> 17) % h->cap; while (h->key[i] != ~0ULL && h->key[i] != k) { i++; if (i == h->cap) i = 0; } return i; }
static int hm_get(HMap *h, u64 k) { size_t i = hm_slot(h, k); return h->key[i] == k ? h->val[i] : -1; }
static int hm_put(HMap *h, u64 k, int v) { size_t i = hm_slot(h, k); if (h->key[i] == k) return 0; h->key[i] = k; h->val[i] = v; h->cnt++; return 1; }

/* ---------------- per-order data ---------------- */
typedef struct {
    int nclass;            /* ordinary classes */
    u64 *code;             /* canonical codes, sorted */
    u64 *conv;             /* canonical code of converse */
    int *aut;
    HMap cls;              /* code -> class index */
    int nmerged;
    int *mid;              /* class -> merged id */
    int *rep;              /* merged id -> representative class (least code) */
} Order;

static Order ord_[MAXN + 1];

static int cmp_u64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }

static void build_order(int n) {
    Order *O = &ord_[n];
    size_t cap = 1;
    { size_t want = 4; if (n >= 8) want = 32768; if (n >= 9) want = 1u << 20; if (n <= 7) want = 2048; cap = want; }
    HMap seen; hm_init(&seen, cap);
    u64 *list = malloc(sizeof(u64) * cap); int cnt = 0;
    if (n == 1) { list[cnt++] = 0; hm_put(&seen, 0, 0); }
    else {
        Order *P = &ord_[n - 1];
        uint16_t base[MAXN], t[MAXN];
        for (int c = 0; c < P->nclass; c++) {
            decode(n - 1, P->code[c], base);
            for (unsigned mask = 0; mask < (1u << (n - 1)); mask++) {
                for (int i = 0; i < n - 1; i++) t[i] = base[i] | ((mask >> i) & 1 ? 0 : (1u << (n - 1)));
                /* new vertex n-1 beats i iff bit i of mask */
                t[n - 1] = (uint16_t)mask;
                int aut; u64 cc = canon(n, t, &aut, NULL);
                if (hm_put(&seen, cc, cnt)) list[cnt++] = cc;
            }
        }
    }
    qsort(list, cnt, sizeof(u64), cmp_u64);
    O->nclass = cnt; O->code = list;
    O->conv = malloc(sizeof(u64) * cnt); O->aut = malloc(sizeof(int) * cnt);
    hm_init(&O->cls, cap);
    for (int i = 0; i < cnt; i++) hm_put(&O->cls, list[i], i);
    uint16_t t[MAXN], op[MAXN];
    for (int i = 0; i < cnt; i++) {
        decode(n, list[i], t);
        int aut; u64 c2 = canon(n, t, &aut, NULL);
        if (c2 != list[i]) { fprintf(stderr, "canonical code not idempotent at n=%d\n", n); exit(1); }
        O->aut[i] = aut;
        converse(n, t, op);
        O->conv[i] = canon(n, op, NULL, NULL);
        if (hm_get(&O->cls, O->conv[i]) < 0) { fprintf(stderr, "converse class missing\n"); exit(1); }
    }
    O->mid = malloc(sizeof(int) * cnt); O->rep = malloc(sizeof(int) * cnt);
    int m = 0;
    for (int i = 0; i < cnt; i++) {
        int j = hm_get(&O->cls, O->conv[i]);
        if (j < i) O->mid[i] = O->mid[j];
        else { O->mid[i] = m; O->rep[m] = i; m++; }
    }
    O->nmerged = m;
    free(seen.key); free(seen.val);
}


/* ---------------- node data: deck and path-weighted deck ---------------- */
typedef struct { int len; int idx[MAXN]; long long w[MAXN]; } SV;
static void sv_add(SV *s, int id, long long w) {
    if (!w) return;
    for (int k = 0; k < s->len; k++) if (s->idx[k] == id) { s->w[k] += w; if (!s->w[k]) { for (int j = k; j + 1 < s->len; j++) { s->idx[j] = s->idx[j + 1]; s->w[j] = s->w[j + 1]; } s->len--; } return; }
    int k = s->len++;
    while (k > 0 && s->idx[k - 1] > id) { s->idx[k] = s->idx[k - 1]; s->w[k] = s->w[k - 1]; k--; }
    s->idx[k] = id; s->w[k] = w;
}
static int sv_cmp(const SV *a, const SV *b) {
    if (a->len != b->len) return a->len < b->len ? -1 : 1;
    for (int k = 0; k < a->len; k++) {
        if (a->idx[k] != b->idx[k]) return a->idx[k] < b->idx[k] ? -1 : 1;
        if (a->w[k] != b->w[k]) return a->w[k] < b->w[k] ? -1 : 1;
    }
    return 0;
}
typedef struct { int sc, aut; u64 H; SV deck, W; } NodeL;
static NodeL *nl_[MAXN + 1];

static void build_nodes_l(int n) {
    Order *O = &ord_[n], *P = &ord_[n - 1];
    int M = O->nmerged;
    nl_[n] = calloc(M, sizeof(NodeL));
    uint16_t t[MAXN], sub[MAXN];
    u64 hs[MAXN], he[MAXN];
    for (int m = 0; m < M; m++) {
        NodeL *N = &nl_[n][m];
        int c = O->rep[m];
        decode(n, O->code[c], t);
        N->sc = (O->conv[c] == O->code[c]);
        N->aut = O->aut[c];
        ham_counts(n, t, hs, he);
        u64 H = 0; for (int v = 0; v < n; v++) H += he[v];
        N->H = H;
        for (int x = 0; x < n; x++) {
            delete_vertex(n, t, x, sub);
            int v = P->mid[hm_get(&P->cls, canon(n - 1, sub, NULL, NULL))];
            sv_add(&N->deck, v, 1);
            sv_add(&N->W, v, (long long)(hs[x] + he[x]));
        }
    }
}

static int g_M; static NodeL *g_nodes;
static int cmp_deck(const void *a, const void *b) {
    int i = *(const int *)a, j = *(const int *)b;
    int c = sv_cmp(&g_nodes[i].deck, &g_nodes[j].deck);
    if (c) return c;
    return i < j ? -1 : i > j;
}
static void print_sv(const SV *s) { printf("{"); for (int k = 0; k < s->len; k++) printf("%s%d:%lld", k ? "," : "", s->idx[k], s->w[k]); printf("}"); }

static void deck_collisions(int n) {
    int M = ord_[n].nmerged;
    int *ix = malloc(sizeof(int) * M); for (int i = 0; i < M; i++) ix[i] = i;
    g_M = M; g_nodes = nl_[n];
    qsort(ix, M, sizeof(int), cmp_deck);
    int pairs = 0;
    for (int i = 0; i < M;) {
        int j = i + 1;
        while (j < M && !sv_cmp(&nl_[n][ix[i]].deck, &nl_[n][ix[j]].deck)) j++;
        if (j - i > 1) {
            pairs++;
            printf("  n=%d deck collision #%d (size %d), deck ", n, pairs, j - i); print_sv(&nl_[n][ix[i]].deck); printf("\n");
            for (int k = i; k < j; k++) {
                NodeL *z = &nl_[n][ix[k]];
                printf("     node %d: H=%llu |Aut|=%d %s, code %llu, path-weighted deck W=", ix[k], z->H, z->aut, z->sc ? "self-converse" : "converse pair", ord_[n].code[ord_[n].rep[ix[k]]]);
                print_sv(&z->W); printf("\n");
            }
        }
        i = j;
    }
    printf("  n=%d: %d unweighted-deck collision classes\n", n, pairs);
    free(ix);
}

/* line fraction through tilings: face 0 deletes the terminal vertex n-1, face 1 the initial vertex 0 */
static void line_fraction(int n) {
    Order *O = &ord_[n], *P = &ord_[n - 1];
    int M = O->nmerged, Mp = P->nmerged;
    long long *Phi0 = calloc((size_t)M * Mp, sizeof(long long)), *Phi1 = calloc((size_t)M * Mp, sizeof(long long));
    int pi[64], pj[64], m = 0;
    for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { pi[m] = i; pj[m] = j; m++; }
    uint16_t t[MAXN], tk[MAXN], sub[MAXN], subk[MAXN];
    int refl[64];
    for (int k = 0; k < m; k++) { int ri = n - 1 - pj[k], rj = n - 1 - pi[k]; for (int l = 0; l < m; l++) if (pi[l] == ri && pj[l] == rj) refl[k] = l; }
    long long *blue = calloc(M, sizeof(long long)), *fib = calloc(M, sizeof(long long));
    for (u64 mask = 0; mask < (1ULL << m); mask++) {
        for (int i = 0; i < n; i++) { t[i] = 0; tk[i] = 0; }
        for (int i = 0; i + 1 < n; i++) { t[i] |= 1u << (i + 1); tk[i] |= 1u << (i + 1); }
        for (int k = 0; k < m; k++) {
            if ((mask >> k) & 1) { t[pj[k]] |= 1u << pi[k]; tk[pi[k]] |= 1u << pj[k]; }
            else { t[pi[k]] |= 1u << pj[k]; tk[pj[k]] |= 1u << pi[k]; }
        }
        int u = O->mid[hm_get(&O->cls, canon(n, t, NULL, NULL))];
        int isblue = 1; for (int k = 0; k < m; k++) if (((mask >> k) & 1) != ((mask >> refl[k]) & 1)) { isblue = 0; break; }
        blue[u] += isblue; fib[u]++;
        /* face 0 */
        delete_vertex(n, t, n - 1, sub); delete_vertex(n, tk, n - 1, subk);
        int v = P->mid[hm_get(&P->cls, canon(n - 1, sub, NULL, NULL))];
        int vk = P->mid[hm_get(&P->cls, canon(n - 1, subk, NULL, NULL))];
        Phi0[(size_t)u * Mp + vk]++; Phi0[(size_t)u * Mp + v]--;
        /* face 1 */
        delete_vertex(n, t, 0, sub); delete_vertex(n, tk, 0, subk);
        v = P->mid[hm_get(&P->cls, canon(n - 1, sub, NULL, NULL))];
        vk = P->mid[hm_get(&P->cls, canon(n - 1, subk, NULL, NULL))];
        Phi1[(size_t)u * Mp + vk]++; Phi1[(size_t)u * Mp + v]--;
    }
    int faces_agree = 1, zero_sum = 1, nzero = 0, nsingle = 0;
    for (int u = 0; u < M; u++) for (int v = 0; v < Mp; v++) if (Phi0[(size_t)u * Mp + v] != Phi1[(size_t)u * Mp + v]) faces_agree = 0;
    for (int v = 0; v < Mp; v++) { long long s = 0; for (int u = 0; u < M; u++) s += Phi0[(size_t)u * Mp + v]; if (s) zero_sum = 0; }
    /* Phi type by colour class: [pure-blue, mixed, pure-black] x [zero, single ratio, other] */
    int tab[3][3] = {{0}};
    /* distinctness of Phi */
    int distinct = 0;
    for (int u = 0; u < M; u++) {
        int dup = 0;
        for (int w = 0; w < u && !dup; w++) {
            int eq = 1;
            for (int v = 0; v < Mp; v++) if (Phi0[(size_t)u * Mp + v] != Phi0[(size_t)w * Mp + v]) { eq = 0; break; }
            if (eq) dup = 1;
        }
        if (!dup) distinct++;
        int nz = 0, pos = 0, neg = 0; long long ap = 0, an = 0;
        for (int v = 0; v < Mp; v++) { long long x = Phi0[(size_t)u * Mp + v]; if (x) { nz++; if (x > 0) { pos++; ap = x; } else { neg++; an = -x; } } }
        if (!nz) nzero++;
        int single = (nz == 2 && pos == 1 && neg == 1 && ap == 1 && an == 1);
        if (single) nsingle++;
        int colour = blue[u] == fib[u] ? 0 : blue[u] == 0 ? 2 : 1;
        tab[colour][!nz ? 0 : single ? 1 : 2]++;
    }
    printf("n=%d Phi type by colour: pure-blue (zero/single/other) = %d/%d/%d; mixed = %d/%d/%d; pure-black = %d/%d/%d\n", n,
           tab[0][0], tab[0][1], tab[0][2], tab[1][0], tab[1][1], tab[1][2], tab[2][0], tab[2][1], tab[2][2]);
    printf("n=%d line fraction Phi: faces agree=%s, sum over nodes = 0: %s, Phi = 0 at %d nodes, Phi a single ratio p_a/p_b at %d nodes, %d distinct values over %d nodes\n",
           n, faces_agree ? "yes" : "NO", zero_sum ? "yes" : "NO", nzero, nsingle, distinct, M);
    {
        for (int u = 0; u < M; u++) {
            NodeL *z = &nl_[n][u];
            int nzc = 0; for (int v = 0; v < Mp; v++) nzc += Phi0[(size_t)u * Mp + v] != 0;
            if (n > 6 && !(blue[u] == fib[u] || nzc <= 2)) continue;
            printf("   node %d (H=%llu |Aut|=%d %s, %s, tilings %lld): Phi = ", u, z->H, z->aut, z->sc ? "sc" : "pair", blue[u] == fib[u] ? "pure-blue" : blue[u] ? "mixed" : "pure-black", fib[u]);
            int first = 1;
            for (int v = 0; v < Mp; v++) { long long x = Phi0[(size_t)u * Mp + v]; if (x) { printf("%s%+lld*[%d]", first ? "" : " ", x, v); first = 0; } }
            if (first) printf("0");
            printf("\n");
        }
    }
    free(Phi0); free(Phi1); free(blue); free(fib);
}

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 8;
    if (NMAX > MAXN) NMAX = MAXN;
    printf("deletion_tower_line_fraction NMAX=%d\n", NMAX);
    for (int n = 1; n <= NMAX; n++) {
        build_order(n);
        if (n >= 2) build_nodes_l(n);
        if (n >= 3) deck_collisions(n);
        if (n >= 4 && n <= 8) line_fraction(n);
        fflush(stdout);
    }
    return 0;
}
