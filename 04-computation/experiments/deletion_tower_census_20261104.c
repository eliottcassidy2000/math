/*
 * deletion_tower_census_20261104.c  (opus-2026-11-04-S72)
 *
 * Exact census of converse-merged tournament nodes through order NMAX (<= 9),
 * with the path-endpoint-weighted vertex deck of every node:
 *
 *   W(T) = sum over vertices x of (h_s(x) + h_e(x)) * e_[T - x],
 *
 * where h_s(x) / h_e(x) count Hamiltonian paths of T starting / ending at x and
 * [T - x] is the converse-merged class of the vertex-deleted subtournament.
 * Lemma W (note, section 1): THM-796's face matrix is R_n(u,v) = W(T)_v / (2|Aut T|)
 * for a self-converse node u = {T} and W(T)_v / |Aut T| for a converse pair.
 *
 * Usage:  deletion_tower_census NMAX [DATA_DIR]
 *   prints, for every order n <= NMAX: class counts, merged counts, injectivity of
 *   primitive / raw face profiles, of the unweighted merged deck (multiset of
 *   [T - x]), of the support set; and for n <= TILECHECK it recomputes R_n from all
 *   2^binom(n-1,2) fixed-path tilings (both endpoint faces) and compares with Lemma W.
 *   If DATA_DIR is given and n <= DATAMAX, writes DATA_DIR/tower_n.txt with one line
 *   per merged node: id sc aut H tilings score W(sparse) R(sparse).
 *
 * Canonical form: lexicographically least adjacency code over all vertex orders that
 * respect an isomorphism-invariant ordered partition (score and 3-cycle count, then
 * 1-dimensional colour refinement on out-neighbourhoods).  Code bits are ordered by
 * (q = 1..n-1, p = 0..q-1), bit = 1 iff the vertex at position p beats the one at q;
 * the number of orders attaining the least code is |Aut T|.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 9
#define TILECHECK 7
#define DATAMAX 9

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

/* ---------------- sparse vectors ---------------- */
typedef struct { int len; int idx[MAXN]; u64 w[MAXN]; } SVec;

static void sv_add(SVec *s, int id, u64 w) {
    if (!w) return;
    for (int k = 0; k < s->len; k++) if (s->idx[k] == id) { s->w[k] += w; return; }
    int k = s->len++;
    while (k > 0 && s->idx[k - 1] > id) { s->idx[k] = s->idx[k - 1]; s->w[k] = s->w[k - 1]; k--; }
    s->idx[k] = id; s->w[k] = w;
}
static u64 gcdu(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }
static void sv_prim(const SVec *s, SVec *p) {
    u64 g = 0; for (int k = 0; k < s->len; k++) g = gcdu(g, s->w[k]);
    *p = *s; if (g) for (int k = 0; k < p->len; k++) p->w[k] /= g;
}
static int sv_cmp(const SVec *a, const SVec *b) {
    if (a->len != b->len) return a->len < b->len ? -1 : 1;
    for (int k = 0; k < a->len; k++) {
        if (a->idx[k] != b->idx[k]) return a->idx[k] < b->idx[k] ? -1 : 1;
        if (a->w[k] != b->w[k]) return a->w[k] < b->w[k] ? -1 : 1;
    }
    return 0;
}

typedef struct {
    int sc, aut; u64 H, tilings; int score[MAXN];
    SVec W, R, Wp, deck, deckset, Wsup;   /* deck: multiplicities of [T-x]; deckset: its support; Wsup: support of W */
    SVec We, Ws;                     /* ordinary (unmerged) end / start weighted decks over ordinary classes of T(rep) */
} Node;

static Node *nodes_[MAXN + 1];

static int g_sortn; static Node *g_sortnodes; static int g_sortfield;
static int cmp_nodeidx(const void *a, const void *b) {
    int i = *(const int *)a, j = *(const int *)b;
    const Node *x = &g_sortnodes[i], *y = &g_sortnodes[j];
    int c;
    switch (g_sortfield) {
        case 0: c = sv_cmp(&x->Wp, &y->Wp); break;
        case 1: c = sv_cmp(&x->R, &y->R); break;
        case 2: c = sv_cmp(&x->deck, &y->deck); break;
        case 3: c = sv_cmp(&x->deckset, &y->deckset); break;
        case 4: c = sv_cmp(&x->Wsup, &y->Wsup); break;
        default: c = 0;
    }
    if (c) return c;
    return i < j ? -1 : i > j;
}
static int field_cmp(const Node *x, const Node *y, int field) {
    switch (field) {
        case 0: return sv_cmp(&x->Wp, &y->Wp);
        case 1: return sv_cmp(&x->R, &y->R);
        case 2: return sv_cmp(&x->deck, &y->deck);
        case 3: return sv_cmp(&x->deckset, &y->deckset);
        default: return sv_cmp(&x->Wsup, &y->Wsup);
    }
}

/* rank of the face matrix R_n over Q, certified from below by the rank of the Gram matrix R^T R
 * modulo the prime 2^31-1 (rank_Q(R) = rank_Q(R^T R) >= rank_p(R^T R)); returns that modular rank. */
static long long modrank_gram(int n) {
    const u64 P = 2147483647ULL;
    int M = ord_[n].nmerged, C = ord_[n - 1].nmerged;
    u64 *G = calloc((size_t)C * C, sizeof(u64));
    for (int m = 0; m < M; m++) {
        const SVec *r = &nodes_[n][m].R;
        for (int a = 0; a < r->len; a++) for (int b = 0; b < r->len; b++) {
            size_t k = (size_t)r->idx[a] * C + r->idx[b];
            G[k] = (G[k] + (r->w[a] % P) * (r->w[b] % P)) % P;
        }
    }
    long long rk = 0;
    for (int col = 0; col < C && rk < C; col++) {
        int piv = -1;
        for (int i = (int)rk; i < C; i++) if (G[(size_t)i * C + col]) { piv = i; break; }
        if (piv < 0) continue;
        if (piv != rk) for (int j = 0; j < C; j++) { u64 t = G[(size_t)piv * C + j]; G[(size_t)piv * C + j] = G[(size_t)rk * C + j]; G[(size_t)rk * C + j] = t; }
        /* inverse of pivot */
        u64 a = G[(size_t)rk * C + col], e = P - 2, inv = 1;
        while (e) { if (e & 1) inv = inv * a % P; a = a * a % P; e >>= 1; }
        for (int j = col; j < C; j++) G[(size_t)rk * C + j] = G[(size_t)rk * C + j] * inv % P;
        for (int i = 0; i < C; i++) {
            if (i == rk) continue;
            u64 f = G[(size_t)i * C + col];
            if (!f) continue;
            u64 *ri = &G[(size_t)i * C], *rp = &G[(size_t)rk * C];
            for (int j = col; j < C; j++) if (rp[j]) ri[j] = (ri[j] + (P - f) * rp[j]) % P;
        }
        rk++;
    }
    free(G);
    return rk;
}

/* number of nodes lying in a collision class and number of classes, for a field */
static void collisions(int n, int field, int *ncls, int *nbad, int *maxsize, int printlimit, const char *name) {
    int M = ord_[n].nmerged;
    int *ix = malloc(sizeof(int) * M); for (int i = 0; i < M; i++) ix[i] = i;
    g_sortnodes = nodes_[n]; g_sortfield = field; g_sortn = n;
    qsort(ix, M, sizeof(int), cmp_nodeidx);
    int classes = 0, bad = 0, mx = 0, printed = 0;
    for (int i = 0; i < M;) {
        int j = i + 1;
        while (j < M) {
            const Node *x = &nodes_[n][ix[i]], *y = &nodes_[n][ix[j]];
            if (field_cmp(x, y, field)) break;
            j++;
        }
        classes++;
        if (j - i > 1) {
            bad += j - i;
            if (printed < printlimit) {
                printf("    %s collision class of size %d:", name, j - i);
                for (int k = i; k < j; k++) { const Node *z = &nodes_[n][ix[k]]; printf(" [node %d H=%llu aut=%d sc=%d]", ix[k], z->H, z->aut, z->sc); }
                printf("\n");
                printed++;
            }
        }
        if (j - i > mx) mx = j - i;
        i = j;
    }
    *ncls = classes; *nbad = bad; *maxsize = mx;
    free(ix);
}

static void build_nodes(int n) {
    Order *O = &ord_[n], *P = &ord_[n - 1];
    int M = O->nmerged;
    nodes_[n] = calloc(M, sizeof(Node));
    uint16_t t[MAXN], sub[MAXN];
    u64 hs[MAXN], he[MAXN];
    for (int m = 0; m < M; m++) {
        Node *N = &nodes_[n][m];
        int c = O->rep[m];
        decode(n, O->code[c], t);
        N->sc = (O->conv[c] == O->code[c]);
        N->aut = O->aut[c];
        ham_counts(n, t, hs, he);
        u64 H = 0; for (int v = 0; v < n; v++) H += he[v];
        u64 H2 = 0; for (int v = 0; v < n; v++) H2 += hs[v];
        if (H != H2) { fprintf(stderr, "path count mismatch\n"); exit(1); }
        N->H = H;
        N->tilings = N->sc ? H / N->aut : 2 * H / N->aut;
        if (H % N->aut) { fprintf(stderr, "H not divisible by aut\n"); exit(1); }
        int sc[MAXN]; for (int v = 0; v < n; v++) sc[v] = popc16(t[v]);
        for (int i = 1; i < n; i++) { int x = sc[i], j = i - 1; while (j >= 0 && sc[j] > x) { sc[j + 1] = sc[j]; j--; } sc[j + 1] = x; }
        memcpy(N->score, sc, sizeof(sc));
        for (int x = 0; x < n; x++) {
            delete_vertex(n, t, x, sub);
            u64 cc = canon(n - 1, sub, NULL, NULL);
            int cls = hm_get(&P->cls, cc);
            if (cls < 0) { fprintf(stderr, "deleted class missing\n"); exit(1); }
            int v = P->mid[cls];
            sv_add(&N->W, v, hs[x] + he[x]);
            sv_add(&N->deck, v, 1);
            sv_add(&N->We, cls, he[x]);
            sv_add(&N->Ws, cls, hs[x]);
        }
        for (int k = 0; k < N->deck.len; k++) { N->deckset.idx[k] = N->deck.idx[k]; N->deckset.w[k] = 1; }
        N->deckset.len = N->deck.len;
        for (int k = 0; k < N->W.len; k++) { N->Wsup.idx[k] = N->W.idx[k]; N->Wsup.w[k] = 1; }
        N->Wsup.len = N->W.len;
        sv_prim(&N->W, &N->Wp);
        /* raw face row R = W / (2|Aut|) (self-converse) or W / |Aut| (pair) */
        N->R = N->W;
        u64 den = N->sc ? 2ULL * N->aut : (u64)N->aut;
        for (int k = 0; k < N->R.len; k++) {
            if (N->R.w[k] % den) { fprintf(stderr, "Lemma W integrality failure n=%d node=%d\n", n, m); exit(1); }
            N->R.w[k] /= den;
        }
    }
}

/* direct face matrix from all fixed-path tilings (path 0->1->...->n-1; free pair (i,j), j>=i+2;
 * bit set => j -> i).  face 0 deletes the terminal vertex n-1, face 1 deletes the initial vertex 0. */
static int tile_check(int n) {
    Order *O = &ord_[n], *P = &ord_[n - 1];
    int M = O->nmerged, Mp = P->nmerged;
    u64 *R0 = calloc((size_t)M * Mp, sizeof(u64)), *R1 = calloc((size_t)M * Mp, sizeof(u64));
    int pi[64], pj[64], m = 0;
    for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { pi[m] = i; pj[m] = j; m++; }
    uint16_t t[MAXN], sub[MAXN];
    u64 *fib = calloc(M, sizeof(u64));
    for (u64 mask = 0; mask < (1ULL << m); mask++) {
        for (int i = 0; i < n; i++) t[i] = 0;
        for (int i = 0; i + 1 < n; i++) t[i] |= 1u << (i + 1);
        for (int k = 0; k < m; k++) { if ((mask >> k) & 1) t[pj[k]] |= 1u << pi[k]; else t[pi[k]] |= 1u << pj[k]; }
        int u = O->mid[hm_get(&O->cls, canon(n, t, NULL, NULL))];
        fib[u]++;
        delete_vertex(n, t, n - 1, sub);
        int v0 = P->mid[hm_get(&P->cls, canon(n - 1, sub, NULL, NULL))];
        delete_vertex(n, t, 0, sub);
        int v1 = P->mid[hm_get(&P->cls, canon(n - 1, sub, NULL, NULL))];
        R0[(size_t)u * Mp + v0]++; R1[(size_t)u * Mp + v1]++;
    }
    int bad = 0;
    for (int u = 0; u < M; u++) {
        Node *N = &nodes_[n][u];
        if (fib[u] != N->tilings) bad++;
        for (int v = 0; v < Mp; v++) {
            u64 pred = 0;
            for (int k = 0; k < N->R.len; k++) if (N->R.idx[k] == v) pred = N->R.w[k];
            if (R0[(size_t)u * Mp + v] != pred || R1[(size_t)u * Mp + v] != pred) bad++;
        }
    }
    free(R0); free(R1); free(fib);
    return bad;
}

static void write_data(int n, const char *dir) {
    char path[512]; snprintf(path, sizeof path, "%s/tower_%d.txt", dir, n);
    FILE *f = fopen(path, "w"); if (!f) { perror(path); return; }
    Order *O = &ord_[n];
    fprintf(f, "# n=%d merged nodes=%d ordinary classes=%d; fields: id sc aut H tilings score | W | R | canonical_code\n", n, O->nmerged, O->nclass);
    for (int m = 0; m < O->nmerged; m++) {
        Node *N = &nodes_[n][m];
        fprintf(f, "%d %d %d %llu %llu ", m, N->sc, N->aut, N->H, N->tilings);
        for (int v = 0; v < n; v++) fprintf(f, "%s%d", v ? "," : "", N->score[v]);
        fprintf(f, " |");
        for (int k = 0; k < N->W.len; k++) fprintf(f, " %d:%llu", N->W.idx[k], N->W.w[k]);
        fprintf(f, " |");
        for (int k = 0; k < N->R.len; k++) fprintf(f, " %d:%llu", N->R.idx[k], N->R.w[k]);
        fprintf(f, " | %llu\n", O->code[O->rep[m]]);
    }
    fclose(f);
}

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 8;
    const char *dir = argc > 2 ? argv[2] : NULL;
    if (NMAX > MAXN) NMAX = MAXN;
    printf("deletion_tower_census NMAX=%d\n", NMAX);
    for (int n = 1; n <= NMAX; n++) {
        build_order(n);
        Order *O = &ord_[n];
        int sc = 0; for (int i = 0; i < O->nclass; i++) sc += (O->conv[i] == O->code[i]);
        printf("n=%d ordinary_classes=%d self_converse=%d merged_nodes=%d (Burnside (a+s)/2=%d)\n", n, O->nclass, sc, O->nmerged, (O->nclass + sc) / 2);
        if (n >= 2) {
            build_nodes(n);
            u64 tot = 0; for (int m = 0; m < O->nmerged; m++) tot += nodes_[n][m].tilings;
            printf("  sum of node tilings = %llu (2^binom(n-1,2) = %llu)\n", tot, 1ULL << ((n - 1) * (n - 2) / 2));
            if (n >= 3) {
                int ncls, nbad, mx;
                collisions(n, 0, &ncls, &nbad, &mx, 5, "primitive-W");
                printf("  primitive face profile: %d distinct of %d nodes; nodes in collisions=%d; max class=%d\n", ncls, O->nmerged, nbad, mx);
                collisions(n, 1, &ncls, &nbad, &mx, 0, "raw-R");
                printf("  raw face row R: %d distinct\n", ncls);
                collisions(n, 2, &ncls, &nbad, &mx, 3, "merged-deck");
                printf("  unweighted merged deck (multiset of [T-x]): %d distinct; nodes in collisions=%d; max class=%d\n", ncls, nbad, mx);
                collisions(n, 3, &ncls, &nbad, &mx, 0, "deck-support");
                printf("  deck support set: %d distinct\n", ncls);
                collisions(n, 4, &ncls, &nbad, &mx, 0, "W-support");
                printf("  support of the face row (endpoint cards only): %d distinct\n", ncls);
                int axis = 0; for (int m = 0; m < O->nmerged; m++) axis += (nodes_[n][m].W.len == 1);
                printf("  axis-profile nodes (face row supported on one lower node): %d\n", axis);
                long long rk = modrank_gram(n);
                printf("  rank of R_n: modular Gram rank = %lld of %d columns (%s)\n", rk, ord_[n - 1].nmerged, rk == ord_[n - 1].nmerged ? "full column rank over Q certified" : "not certified full");
            }
            if (n >= 3 && n <= TILECHECK) {
                int bad = tile_check(n);
                printf("  Lemma W check against all fixed-path tilings (both endpoint faces, fibre sizes): %s (%d mismatches)\n", bad ? "FAIL" : "PASS", bad);
            }
            if (dir && n <= DATAMAX) write_data(n, dir);
        }
        fflush(stdout);
    }
    return 0;
}
