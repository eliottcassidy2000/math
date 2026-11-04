/*
 * audit_tower.c  -- independent audit code for THM-4863 / HYP-9455 (S72 audit).
 * Written from scratch; shares no code with the session's programs.
 *
 * Canonical form: individualization-refinement (equitable refinement on out-neighbour colour counts, starting from
 * (score, 3-cycles through v, sum of out-neighbour scores)), first non-singleton cell individualized, full search tree
 * (no pruning), code = bits (j = 1..n-1, i = 0..j-1, MSB first), bit = 1 iff the vertex at position j BEATS the one at
 * position i; the canonical code is the MAXIMUM over leaves; |Aut| = number of leaves attaining it.
 *
 * Usage:
 *   audit_tower NMAX OUTDIR          (NMAX <= 10)  census, decks, W, collisions, rank (n <= 9), dumps for n <= 9
 *   audit_tower NMAX OUTDIR tiles K  also checks Lemma W / Phi / line classes from all fixed-path tilings for n <= K (<= 8)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

typedef unsigned long long u64;
typedef uint16_t u16;
typedef uint32_t u32;
#define MAXN 10

/* ------------------------------------------------------------------ canonical form */
static int sigcmp(const int *a, const int *b, int len) {
    for (int i = 0; i < len; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}
static void refine(int n, const u16 *adj, int *col, int *ncol) {
    for (;;) {
        int k = *ncol;
        int sig[MAXN][MAXN + 1];
        for (int v = 0; v < n; v++) {
            sig[v][0] = col[v];
            for (int c = 0; c < k; c++) sig[v][c + 1] = 0;
            unsigned m = adj[v];
            while (m) { int u = __builtin_ctz(m); m &= m - 1; sig[v][col[u] + 1]++; }
        }
        int ord[MAXN];
        for (int i = 0; i < n; i++) ord[i] = i;
        for (int i = 1; i < n; i++) {
            int x = ord[i], j = i - 1;
            while (j >= 0 && sigcmp(sig[ord[j]], sig[x], k + 1) > 0) { ord[j + 1] = ord[j]; j--; }
            ord[j + 1] = x;
        }
        int nc = 0, newcol[MAXN];
        for (int i = 0; i < n; i++) {
            if (i > 0 && sigcmp(sig[ord[i]], sig[ord[i - 1]], k + 1)) nc++;
            newcol[ord[i]] = nc;
        }
        nc++;
        memcpy(col, newcol, sizeof(int) * n);
        if (nc == k) return;
        *ncol = nc;
    }
}
typedef struct { int n; const u16 *adj; u64 best; long long nbest; int have; } CC;
static void search(CC *C, const int *col, int ncol) {
    int n = C->n;
    if (ncol == n) {
        int at[MAXN];
        for (int v = 0; v < n; v++) at[col[v]] = v;
        u64 code = 0;
        for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) code = (code << 1) | ((C->adj[at[j]] >> at[i]) & 1u);
        if (!C->have || code > C->best) { C->best = code; C->nbest = 1; C->have = 1; }
        else if (code == C->best) C->nbest++;
        return;
    }
    int cnt[MAXN] = {0};
    for (int v = 0; v < n; v++) cnt[col[v]]++;
    int c = 0;
    while (cnt[c] == 1) c++;
    for (int v = 0; v < n; v++) {
        if (col[v] != c) continue;
        int col2[MAXN];
        for (int u = 0; u < n; u++) col2[u] = col[u] > c ? col[u] + 1 : (col[u] == c ? c + 1 : col[u]);
        col2[v] = c;
        int nc2 = ncol + 1;
        refine(n, C->adj, col2, &nc2);
        search(C, col2, nc2);
    }
}
static u64 canon(int n, const u16 *adj, long long *aut) {
    int key[MAXN], score[MAXN];
    for (int v = 0; v < n; v++) score[v] = __builtin_popcount(adj[v]);
    for (int v = 0; v < n; v++) {
        int c3 = 0, so = 0;
        unsigned in = 0;
        for (int u = 0; u < n; u++) if (u != v && !((adj[v] >> u) & 1)) in |= 1u << u;
        unsigned m = adj[v];
        while (m) { int u = __builtin_ctz(m); m &= m - 1; c3 += __builtin_popcount(adj[u] & in); so += score[u]; }
        key[v] = score[v] * 10000 + c3 * 100 + so;
    }
    /* initial colours = rank of key */
    int col[MAXN], ncol = 0;
    {
        int s[MAXN];
        memcpy(s, key, sizeof(int) * n);
        for (int i = 1; i < n; i++) { int x = s[i], j = i - 1; while (j >= 0 && s[j] > x) { s[j + 1] = s[j]; j--; } s[j + 1] = x; }
        int d[MAXN], nd = 0;
        for (int i = 0; i < n; i++) if (!i || s[i] != s[i - 1]) d[nd++] = s[i];
        for (int v = 0; v < n; v++) for (int i = 0; i < nd; i++) if (d[i] == key[v]) col[v] = i;
        ncol = nd;
    }
    refine(n, adj, col, &ncol);
    CC C; C.n = n; C.adj = adj; C.have = 0; C.nbest = 0; C.best = 0;
    search(&C, col, ncol);
    if (aut) *aut = C.nbest;
    return C.best;
}
static void decode(int n, u64 code, u16 *adj) {
    int b = n * (n - 1) / 2 - 1;
    for (int i = 0; i < n; i++) adj[i] = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
        if ((code >> b) & 1) adj[j] |= 1u << i; else adj[i] |= 1u << j;
        b--;
    }
}
static void converse(int n, const u16 *adj, u16 *op) {
    for (int i = 0; i < n; i++) op[i] = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if ((adj[i] >> j) & 1) op[j] |= 1u << i;
}
static void delv(int n, const u16 *adj, int x, u16 *r) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (i == x) continue;
        unsigned lo = adj[i] & ((1u << x) - 1), hi = (adj[i] >> (x + 1)) << x;
        r[k++] = (u16)(lo | hi);
    }
}
/* paths ending at v (he) and starting at v (hs) */
static void hampaths(int n, const u16 *adj, u64 *hs, u64 *he) {
    static __thread u64 dp[1 << MAXN][MAXN];
    int full = (1 << n) - 1;
    for (int pass = 0; pass < 2; pass++) {
        u16 g[MAXN];
        if (pass == 0) memcpy(g, adj, sizeof(u16) * n); else converse(n, adj, g);
        for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
        for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
        for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
            u64 c = dp[S][v];
            if (!c) continue;
            unsigned m = g[v] & ~(unsigned)S;
            while (m) { int w = __builtin_ctz(m); m &= m - 1; dp[S | (1 << w)][w] += c; }
        }
        for (int v = 0; v < n; v++) (pass == 0 ? he : hs)[v] = dp[full][v];
    }
}
static u64 hamcycles(int n, const u16 *adj) {
    static __thread u64 dp[1 << MAXN][MAXN];
    if (n < 3) return 0;
    int full = (1 << n) - 1;
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    dp[1][0] = 1;
    for (int S = 1; S <= full; S += 2) for (int v = 0; v < n; v++) {
        u64 c = dp[S][v];
        if (!c) continue;
        unsigned m = adj[v] & ~(unsigned)S;
        while (m) { int w = __builtin_ctz(m); m &= m - 1; dp[S | (1 << w)][w] += c; }
    }
    u64 h = 0;
    for (int v = 1; v < n; v++) if (adj[v] & 1u) h += dp[full][v];
    return h;
}

/* ------------------------------------------------------------------ per-order data */
typedef struct {
    int n; long long ncls; u64 *code;
    int *conv; int *aut; int *mid; long long nm; int *rep; long long nsc;
    u64 *H, *hc;
    u32 *odeck;     /* ncls x n sorted lower class ids */
    u32 *mdeck;     /* nm x n sorted lower merged ids */
    u32 *Widx, *Ww; /* nm x n sparse W (merged ids, weights); length in Wlen */
    unsigned char *Wlen;
    u64 *Wg;        /* gcd of W */
} Ord;
static Ord O[MAXN + 1];

static long long find(const Ord *P, u64 code) {
    long long lo = 0, hi = P->ncls - 1;
    while (lo <= hi) { long long mid = (lo + hi) / 2; if (P->code[mid] == code) return mid; if (P->code[mid] < code) lo = mid + 1; else hi = mid - 1; }
    return -1;
}
static void radix_sort(u64 *a, u64 *tmp, size_t N, int bits) {
    for (int sh = 0; sh < bits; sh += 16) {
        size_t *cnt = calloc(65537, sizeof(size_t));
        for (size_t i = 0; i < N; i++) cnt[((a[i] >> sh) & 0xFFFF) + 1]++;
        for (int i = 0; i < 65536; i++) cnt[i + 1] += cnt[i];
        for (size_t i = 0; i < N; i++) tmp[cnt[(a[i] >> sh) & 0xFFFF]++] = a[i];
        memcpy(a, tmp, N * sizeof(u64));
        free(cnt);
    }
}
static void generate(int n) {
    Ord *Q = &O[n];
    Q->n = n;
    if (n == 1) { Q->ncls = 1; Q->code = malloc(8); Q->code[0] = 0; return; }
    Ord *P = &O[n - 1];
    size_t per = (size_t)1 << (n - 1), tot = (size_t)P->ncls * per;
    u64 *all = malloc(tot * sizeof(u64));
#pragma omp parallel for schedule(dynamic, 64)
    for (long long c = 0; c < P->ncls; c++) {
        u16 base[MAXN], t[MAXN];
        decode(n - 1, P->code[c], base);
        for (size_t mask = 0; mask < per; mask++) {
            for (int i = 0; i < n - 1; i++) t[i] = base[i] | (((mask >> i) & 1) ? 0 : (1u << (n - 1)));
            t[n - 1] = (u16)mask;
            all[(size_t)c * per + mask] = canon(n, t, NULL);
        }
    }
    u64 *tmp = malloc(tot * sizeof(u64));
    radix_sort(all, tmp, tot, 48);
    free(tmp);
    size_t m = 0;
    for (size_t i = 0; i < tot; i++) if (!i || all[i] != all[i - 1]) all[m++] = all[i];
    Q->ncls = (long long)m;
    Q->code = realloc(all, m * sizeof(u64));
}
static int cmp_u32(const void *a, const void *b) { u32 x = *(const u32 *)a, y = *(const u32 *)b; return x < y ? -1 : x > y; }
static u64 gcd64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

static void analyse(int n) {
    Ord *Q = &O[n], *P = n > 1 ? &O[n - 1] : NULL;
    long long N = Q->ncls;
    Q->conv = malloc(sizeof(int) * N); Q->aut = malloc(sizeof(int) * N); Q->mid = malloc(sizeof(int) * N);
    Q->H = malloc(sizeof(u64) * N); Q->hc = malloc(sizeof(u64) * N);
    Q->odeck = n > 1 ? malloc(sizeof(u32) * N * n) : NULL;
    int bad = 0;
#pragma omp parallel for schedule(dynamic, 256) reduction(+:bad)
    for (long long i = 0; i < N; i++) {
        u16 t[MAXN], op[MAXN], sub[MAXN];
        decode(n, Q->code[i], t);
        long long aut;
        u64 c2 = canon(n, t, &aut);
        if (c2 != Q->code[i]) bad++;
        Q->aut[i] = (int)aut;
        converse(n, t, op);
        long long j = find(Q, canon(n, op, NULL));
        if (j < 0) bad++;
        Q->conv[i] = (int)j;
        u64 hs[MAXN], he[MAXN];
        hampaths(n, t, hs, he);
        u64 H = 0, H2 = 0;
        for (int v = 0; v < n; v++) { H += he[v]; H2 += hs[v]; }
        if (H != H2) bad++;
        Q->H[i] = H;
        Q->hc[i] = hamcycles(n, t);
        if (n > 1) {
            u32 *d = &Q->odeck[(size_t)i * n];
            for (int x = 0; x < n; x++) {
                delv(n, t, x, sub);
                long long k = find(P, canon(n - 1, sub, NULL));
                if (k < 0) bad++;
                d[x] = (u32)k;
            }
            qsort(d, n, sizeof(u32), cmp_u32);
        }
    }
    if (bad) { fprintf(stderr, "consistency failure at n=%d (%d)\n", n, bad); exit(1); }
    long long m = 0, sc = 0;
    Q->rep = malloc(sizeof(int) * N);
    for (long long i = 0; i < N; i++) {
        if (Q->conv[i] == i) sc++;
        if (Q->conv[i] < i) Q->mid[i] = Q->mid[Q->conv[i]];
        else { Q->mid[i] = (int)m; Q->rep[m] = (int)i; m++; }
    }
    Q->nm = m; Q->nsc = sc;
    if (n < 2) return;
    /* merged nodes: W and merged deck */
    Q->mdeck = malloc(sizeof(u32) * m * n);
    Q->Widx = malloc(sizeof(u32) * m * n); Q->Ww = malloc(sizeof(u32) * m * n);
    Q->Wlen = malloc(m); Q->Wg = malloc(sizeof(u64) * m);
#pragma omp parallel for schedule(dynamic, 256)
    for (long long u = 0; u < m; u++) {
        long long i = Q->rep[u];
        u16 t[MAXN], sub[MAXN];
        decode(n, Q->code[i], t);
        u64 hs[MAXN], he[MAXN];
        hampaths(n, t, hs, he);
        u32 card[MAXN]; u64 w[MAXN];
        for (int x = 0; x < n; x++) {
            delv(n, t, x, sub);
            long long k = find(P, canon(n - 1, sub, NULL));
            card[x] = (u32)P->mid[k];
            w[x] = hs[x] + he[x];
        }
        /* merged deck */
        u32 *md = &Q->mdeck[(size_t)u * n];
        memcpy(md, card, sizeof(u32) * n);
        qsort(md, n, sizeof(u32), cmp_u32);
        /* W sparse: aggregate by card */
        u32 idx[MAXN]; u64 ww[MAXN]; int len = 0;
        for (int x = 0; x < n; x++) {
            if (!w[x]) continue;
            int f = -1;
            for (int k = 0; k < len; k++) if (idx[k] == card[x]) f = k;
            if (f < 0) { idx[len] = card[x]; ww[len] = w[x]; len++; } else ww[f] += w[x];
        }
        for (int a = 1; a < len; a++) { u32 xi = idx[a]; u64 xw = ww[a]; int b = a - 1; while (b >= 0 && idx[b] > xi) { idx[b + 1] = idx[b]; ww[b + 1] = ww[b]; b--; } idx[b + 1] = xi; ww[b + 1] = xw; }
        u64 g = 0;
        for (int k = 0; k < len; k++) g = gcd64(g, ww[k]);
        for (int k = 0; k < len; k++) { Q->Widx[(size_t)u * n + k] = idx[k]; Q->Ww[(size_t)u * n + k] = (u32)ww[k]; }
        Q->Wlen[u] = (unsigned char)len; Q->Wg[u] = g;
    }
}

/* comparators over merged nodes / classes */
static Ord *gQ;
static int cmp_primW(const void *a, const void *b) {
    long long x = *(const long long *)a, y = *(const long long *)b;
    int n = gQ->n;
    if (gQ->Wlen[x] != gQ->Wlen[y]) return gQ->Wlen[x] < gQ->Wlen[y] ? -1 : 1;
    for (int k = 0; k < gQ->Wlen[x]; k++) {
        u32 ix = gQ->Widx[(size_t)x * n + k], iy = gQ->Widx[(size_t)y * n + k];
        if (ix != iy) return ix < iy ? -1 : 1;
        u64 wx = gQ->Ww[(size_t)x * n + k] / gQ->Wg[x], wy = gQ->Ww[(size_t)y * n + k] / gQ->Wg[y];
        if (wx != wy) return wx < wy ? -1 : 1;
    }
    return x < y ? -1 : x > y;
}
static int eq_primW(long long x, long long y) {
    int n = gQ->n;
    if (gQ->Wlen[x] != gQ->Wlen[y]) return 0;
    for (int k = 0; k < gQ->Wlen[x]; k++) {
        if (gQ->Widx[(size_t)x * n + k] != gQ->Widx[(size_t)y * n + k]) return 0;
        if (gQ->Ww[(size_t)x * n + k] / gQ->Wg[x] != gQ->Ww[(size_t)y * n + k] / gQ->Wg[y]) return 0;
    }
    return 1;
}
static int cmp_mdeck(const void *a, const void *b) {
    long long x = *(const long long *)a, y = *(const long long *)b;
    int n = gQ->n;
    for (int k = 0; k < n; k++) { u32 p = gQ->mdeck[(size_t)x * n + k], q = gQ->mdeck[(size_t)y * n + k]; if (p != q) return p < q ? -1 : 1; }
    return x < y ? -1 : x > y;
}
static int eq_mdeck(long long x, long long y) {
    int n = gQ->n;
    for (int k = 0; k < n; k++) if (gQ->mdeck[(size_t)x * n + k] != gQ->mdeck[(size_t)y * n + k]) return 0;
    return 1;
}
static int cmp_odeck(const void *a, const void *b) {
    long long x = *(const long long *)a, y = *(const long long *)b;
    int n = gQ->n;
    for (int k = 0; k < n; k++) { u32 p = gQ->odeck[(size_t)x * n + k], q = gQ->odeck[(size_t)y * n + k]; if (p != q) return p < q ? -1 : 1; }
    return x < y ? -1 : x > y;
}
static int eq_odeck(long long x, long long y) {
    int n = gQ->n;
    for (int k = 0; k < n; k++) if (gQ->odeck[(size_t)x * n + k] != gQ->odeck[(size_t)y * n + k]) return 0;
    return 1;
}

/* rank of R_n (R = W / (2|Aut|) or W/|Aut|) over GF(p): full rank of R^T D R mod p certifies full column rank over Q */
static long long rank_RtDR(int n, u64 seed) {
    const u64 p = 2147483629ULL; /* prime < 2^31 */
    Ord *Q = &O[n];
    long long M = Q->nm, C = O[n - 1].nm;
    u64 *G = calloc((size_t)C * C, sizeof(u64));
    u64 s = seed;
    for (long long u = 0; u < M; u++) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        u64 d = (s >> 33) % p; if (!d) d = 1;
        long long i = Q->rep[u];
        u64 den = (Q->conv[i] == i) ? 2ULL * Q->aut[i] : (u64)Q->aut[i];
        int len = Q->Wlen[u];
        u64 r[MAXN]; u32 ix[MAXN];
        for (int k = 0; k < len; k++) {
            u64 w = Q->Ww[(size_t)u * n + k];
            if (w % den) { fprintf(stderr, "Lemma W integrality failure n=%d node %lld\n", n, u); exit(1); }
            r[k] = (w / den) % p; ix[k] = Q->Widx[(size_t)u * n + k];
        }
        for (int a = 0; a < len; a++) for (int b = 0; b < len; b++) {
            size_t pos = (size_t)ix[a] * C + ix[b];
            G[pos] = (G[pos] + d * r[a] % p * r[b]) % p;
        }
    }
    long long rk = 0;
    for (long long col = 0; col < C; col++) {
        long long piv = -1;
        for (long long i = rk; i < C; i++) if (G[(size_t)i * C + col]) { piv = i; break; }
        if (piv < 0) continue;
        if (piv != rk) for (long long j = 0; j < C; j++) { u64 t = G[(size_t)piv * C + j]; G[(size_t)piv * C + j] = G[(size_t)rk * C + j]; G[(size_t)rk * C + j] = t; }
        u64 a = G[(size_t)rk * C + col], e = p - 2, inv = 1;
        while (e) { if (e & 1) inv = inv * a % p; a = a * a % p; e >>= 1; }
        for (long long j = col; j < C; j++) G[(size_t)rk * C + j] = G[(size_t)rk * C + j] * inv % p;
        u64 *prow = &G[(size_t)rk * C];
#pragma omp parallel for schedule(static)
        for (long long i = 0; i < C; i++) {
            if (i == rk) continue;
            u64 f = G[(size_t)i * C + col];
            if (!f) continue;
            u64 *ri = &G[(size_t)i * C];
            u64 nf = p - f;
            for (long long j = col; j < C; j++) if (prow[j]) ri[j] = (ri[j] + nf * prow[j]) % p;
        }
        rk++;
    }
    free(G);
    return rk;
}

static void print_class(int n, long long i) {
    Ord *Q = &O[n];
    printf("code=%llu H=%llu hc=%llu aut=%d %s", Q->code[i], Q->H[i], Q->hc[i], Q->aut[i], Q->conv[i] == i ? "self-converse" : "non-self-converse");
}

static void collisions_report(int n) {
    Ord *Q = &O[n];
    gQ = Q;
    long long M = Q->nm;
    long long *ix = malloc(sizeof(long long) * M);
    /* primitive W */
    for (long long u = 0; u < M; u++) ix[u] = u;
    qsort(ix, M, sizeof(long long), cmp_primW);
    long long distinct = 0, cls = 0;
    for (long long a = 0; a < M;) { long long b = a + 1; while (b < M && eq_primW(ix[a], ix[b])) b++; distinct++; if (b - a > 1) cls++; a = b; }
    printf("n=%d primitive W: %lld distinct of %lld merged nodes (%lld collision classes)\n", n, distinct, M, cls);
    /* raw R distinct */
    /* merged deck */
    for (long long u = 0; u < M; u++) ix[u] = u;
    qsort(ix, M, sizeof(long long), cmp_mdeck);
    distinct = 0; cls = 0;
    for (long long a = 0; a < M;) {
        long long b = a + 1;
        while (b < M && eq_mdeck(ix[a], ix[b])) b++;
        distinct++;
        if (b - a > 1) {
            cls++;
            printf("  n=%d merged-deck collision class %lld (size %lld):\n", n, cls, b - a);
            for (long long k = a; k < b; k++) {
                long long i = Q->rep[ix[k]];
                printf("     node %lld: ", ix[k]); print_class(n, i);
                printf(" | W=");
                for (int q = 0; q < Q->Wlen[ix[k]]; q++) printf("%s%u:%u", q ? "," : "", Q->Widx[(size_t)ix[k] * n + q], Q->Ww[(size_t)ix[k] * n + q]);
                printf("\n");
            }
            /* Theorem H check on all pairs of the class */
            for (long long k = a; k < b; k++) for (long long l = k + 1; l < b; l++) {
                long long i = Q->rep[ix[k]], j = Q->rep[ix[l]];
                long long dH = (long long)Q->H[i] - (long long)Q->H[j];
                long long dc = (long long)Q->hc[i] - (long long)Q->hc[j];
                int ok = (n % 2 == 0) ? (dH == 0) : (dH == 2 * dc);
                printf("     Theorem H check: dH=%lld, 2*dhc=%lld, n %s -> %s\n", dH, 2 * dc, n % 2 ? "odd" : "even", ok ? "OK" : "VIOLATED");
            }
        }
        a = b;
    }
    printf("n=%d unweighted merged deck: %lld distinct of %lld (%lld collision classes)\n", n, distinct, M, cls);
    free(ix);
    /* ordinary decks over ordinary classes */
    long long N = Q->ncls;
    long long *jx = malloc(sizeof(long long) * N);
    for (long long i = 0; i < N; i++) jx[i] = i;
    qsort(jx, N, sizeof(long long), cmp_odeck);
    distinct = 0; cls = 0;
    long long convpairs = 0;
    for (long long a = 0; a < N;) {
        long long b = a + 1;
        while (b < N && eq_odeck(jx[a], jx[b])) b++;
        distinct++;
        if (b - a > 1) {
            cls++;
            int isconv = (b - a == 2 && Q->conv[jx[a]] == jx[a + 1]);
            convpairs += isconv;
            printf("  n=%d same ordinary deck (size %lld):", n, b - a);
            for (long long k = a; k < b; k++) { printf(" ["); print_class(n, jx[k]); printf("]"); }
            printf(" -> %s\n", isconv ? "converse pair" : "NOT converse pair");
        }
        a = b;
    }
    printf("n=%d ordinary deck: %lld distinct of %lld classes; %lld same-deck classes (%lld converse pairs)\n", n, distinct, N, cls, convpairs);
    free(jx);
}

static void dump(int n, const char *dir) {
    char path[600];
    Ord *Q = &O[n];
    snprintf(path, sizeof path, "%s/audit_classes_%d.txt", dir, n);
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return; }
    fprintf(f, "# idx code conv aut mid H hc | odeck\n");
    for (long long i = 0; i < Q->ncls; i++) {
        fprintf(f, "%lld %llu %d %d %d %llu %llu |", i, Q->code[i], Q->conv[i], Q->aut[i], Q->mid[i], Q->H[i], Q->hc[i]);
        if (n > 1) for (int k = 0; k < n; k++) fprintf(f, " %u", Q->odeck[(size_t)i * n + k]);
        fprintf(f, "\n");
    }
    fclose(f);
    if (n < 2) return;
    snprintf(path, sizeof path, "%s/audit_nodes_%d.txt", dir, n);
    f = fopen(path, "w");
    if (!f) { perror(path); return; }
    fprintf(f, "# mid repidx code sc aut H hc | W (merged id:weight) | mdeck\n");
    for (long long u = 0; u < Q->nm; u++) {
        long long i = Q->rep[u];
        fprintf(f, "%lld %lld %llu %d %d %llu %llu |", u, i, Q->code[i], Q->conv[i] == i, Q->aut[i], Q->H[i], Q->hc[i]);
        for (int k = 0; k < Q->Wlen[u]; k++) fprintf(f, " %u:%u", Q->Widx[(size_t)u * n + k], Q->Ww[(size_t)u * n + k]);
        fprintf(f, " |");
        for (int k = 0; k < n; k++) fprintf(f, " %u", Q->mdeck[(size_t)u * n + k]);
        fprintf(f, "\n");
    }
    fclose(f);
}

/* ------------------------------------------------------------------ tilings (THM-796 conventions, 0-indexed)
 * vertices 0..n-1, base path (n-1) -> (n-2) -> ... -> 0 (vertex i+1 beats i); free pairs (x,y), x - y >= 2;
 * bit 0: x beats y, bit 1: y beats x.  d_L deletes vertex 0 (path endpoint "1"), d_H deletes vertex n-1.
 * kappa flips every free bit; sigma(t)(x,y) = t(n-1-y, n-1-x). */
static void tiling_checks(int n) {
    Ord *Q = &O[n], *P = &O[n - 1];
    int px[64], py[64], m = 0;
    for (int y = 0; y < n; y++) for (int x = y + 2; x < n; x++) { px[m] = x; py[m] = y; m++; }
    int refl[64];
    for (int k = 0; k < m; k++) { int rx = n - 1 - py[k], ry = n - 1 - px[k]; for (int l = 0; l < m; l++) if (px[l] == rx && py[l] == ry) refl[k] = l; }
    long long M = Q->nm, Mp = P->nm;
    long long *RL = calloc((size_t)M * Mp, sizeof(long long)), *RH = calloc((size_t)M * Mp, sizeof(long long));
    long long *PhL = calloc((size_t)M * Mp, sizeof(long long)), *PhH = calloc((size_t)M * Mp, sizeof(long long));
    long long *fib = calloc(M, sizeof(long long)), *blue = calloc(M, sizeof(long long));
    u64 total = 1ULL << m;
#pragma omp parallel
    {
        long long *rl = calloc((size_t)M * Mp, sizeof(long long)), *rh = calloc((size_t)M * Mp, sizeof(long long));
        long long *pl = calloc((size_t)M * Mp, sizeof(long long)), *ph = calloc((size_t)M * Mp, sizeof(long long));
        long long *fb = calloc(M, sizeof(long long)), *bl = calloc(M, sizeof(long long));
#pragma omp for schedule(dynamic, 1024)
        for (long long mask = 0; mask < (long long)total; mask++) {
            u16 t[MAXN], tk[MAXN], s1[MAXN], s2[MAXN];
            for (int i = 0; i < n; i++) { t[i] = 0; tk[i] = 0; }
            for (int i = 1; i < n; i++) { t[i] |= 1u << (i - 1); tk[i] |= 1u << (i - 1); }
            for (int k = 0; k < m; k++) {
                int b = (mask >> k) & 1;
                if (!b) { t[px[k]] |= 1u << py[k]; tk[py[k]] |= 1u << px[k]; }
                else { t[py[k]] |= 1u << px[k]; tk[px[k]] |= 1u << py[k]; }
            }
            long long u = Q->mid[find(Q, canon(n, t, NULL))];
            int isblue = 1;
            for (int k = 0; k < m; k++) if (((mask >> k) & 1) != ((mask >> refl[k]) & 1)) { isblue = 0; break; }
            fb[u]++; bl[u] += isblue;
            /* d_L: delete vertex 0 */
            delv(n, t, 0, s1); delv(n, tk, 0, s2);
            long long v = P->mid[find(P, canon(n - 1, s1, NULL))], vk = P->mid[find(P, canon(n - 1, s2, NULL))];
            rl[u * Mp + v]++; pl[u * Mp + vk]++; pl[u * Mp + v]--;
            /* d_H: delete vertex n-1 */
            delv(n, t, n - 1, s1); delv(n, tk, n - 1, s2);
            v = P->mid[find(P, canon(n - 1, s1, NULL))]; vk = P->mid[find(P, canon(n - 1, s2, NULL))];
            rh[u * Mp + v]++; ph[u * Mp + vk]++; ph[u * Mp + v]--;
        }
#pragma omp critical
        {
            for (size_t z = 0; z < (size_t)M * Mp; z++) { RL[z] += rl[z]; RH[z] += rh[z]; PhL[z] += pl[z]; PhH[z] += ph[z]; }
            for (long long z = 0; z < M; z++) { fib[z] += fb[z]; blue[z] += bl[z]; }
        }
        free(rl); free(rh); free(pl); free(ph); free(fb); free(bl);
    }
    /* Lemma W comparison */
    long long mism = 0, fibbad = 0;
    for (long long u = 0; u < M; u++) {
        long long i = Q->rep[u];
        int sc = Q->conv[i] == i;
        u64 den = sc ? 2ULL * Q->aut[i] : (u64)Q->aut[i];
        u64 fsz = sc ? Q->H[i] / Q->aut[i] : 2 * Q->H[i] / Q->aut[i];
        if ((u64)fib[u] != fsz) fibbad++;
        for (long long v = 0; v < Mp; v++) {
            u64 pred = 0;
            for (int k = 0; k < Q->Wlen[u]; k++) if (Q->Widx[(size_t)u * n + k] == v) pred = Q->Ww[(size_t)u * n + k];
            if (pred % den) mism++;
            pred /= den;
            if ((u64)RL[u * Mp + v] != pred || (u64)RH[u * Mp + v] != pred) mism++;
        }
    }
    printf("n=%d tilings=%llu: Lemma W vs all tilings (both faces): %lld mismatches; fibre-size mismatches %lld\n", n, total, mism, fibbad);
    /* Phi */
    long long facesdiff = 0, sumnz = 0;
    for (size_t z = 0; z < (size_t)M * Mp; z++) if (PhL[z] != PhH[z]) facesdiff++;
    for (long long v = 0; v < Mp; v++) { long long s = 0; for (long long u = 0; u < M; u++) s += PhL[u * Mp + v]; if (s) sumnz++; }
    long long distinct = 0;
    long long zero = 0, single = 0;
    for (long long u = 0; u < M; u++) {
        int dup = 0;
        for (long long w = 0; w < u && !dup; w++) { int eq = 1; for (long long v = 0; v < Mp; v++) if (PhL[u * Mp + v] != PhL[w * Mp + v]) { eq = 0; break; } if (eq) dup = 1; }
        if (!dup) distinct++;
        int nz = 0, pos = 0, neg = 0; long long ap = 0, an = 0;
        for (long long v = 0; v < Mp; v++) { long long x = PhL[u * Mp + v]; if (x) { nz++; if (x > 0) { pos++; ap = x; } else { neg++; an = -x; } } }
        if (!nz) zero++;
        int sgl = nz == 2 && pos == 1 && neg == 1 && ap == 1 && an == 1;
        if (sgl) single++;
        if (n <= 5 || sgl || blue[u] == fib[u]) {
            long long i = Q->rep[u];
            printf("   n=%d node %lld (H=%llu aut=%d %s, tilings %lld, blue %lld = %s): Phi =", n, u, Q->H[i], Q->aut[i], Q->conv[i] == i ? "sc" : "pair", fib[u], blue[u], blue[u] == fib[u] ? "pure-blue" : blue[u] ? "mixed" : "pure-black");
            int first = 1;
            for (long long v = 0; v < Mp; v++) { long long x = PhL[u * Mp + v]; if (x) { printf(" %+lld[%lld]", x, v); first = 0; } }
            if (first) printf(" 0");
            printf("%s\n", sgl ? "  (single ratio)" : "");
        }
    }
    printf("n=%d Phi: faces differ at %lld entries; nonzero column sums %lld; %lld distinct values over %lld nodes; Phi=0 at %lld; single ratio at %lld\n", n, facesdiff, sumnz, distinct, M, zero, single);
    /* line classes: orbits of <kappa, sigma>; and the adjacency map */
    {
        unsigned char *seen = calloc(total, 1);
        long long orbits = 0;
        u64 full = total - 1;
        /* adjacency set: unordered pair of nodes (and colour) per orbit */
        long long cap = 1 << 22; long long *adjkey = malloc(sizeof(long long) * cap); long long na = 0;
        for (u64 t0 = 0; t0 < total; t0++) {
            if (seen[t0]) continue;
            u64 s = 0;
            for (int k = 0; k < m; k++) if ((t0 >> k) & 1) s |= 1ULL << refl[k];
            u64 orb[4] = {t0, t0 ^ full, s, s ^ full};
            for (int z = 0; z < 4; z++) seen[orb[z]] = 1;
            orbits++;
            if (n <= 7) {
                u16 a[MAXN], b[MAXN];
                for (int i = 0; i < n; i++) { a[i] = 0; b[i] = 0; }
                for (int i = 1; i < n; i++) { a[i] |= 1u << (i - 1); b[i] |= 1u << (i - 1); }
                for (int k = 0; k < m; k++) { if (!((t0 >> k) & 1)) { a[px[k]] |= 1u << py[k]; b[py[k]] |= 1u << px[k]; } else { a[py[k]] |= 1u << px[k]; b[px[k]] |= 1u << py[k]; } }
                long long u1 = Q->mid[find(Q, canon(n, a, NULL))], u2 = Q->mid[find(Q, canon(n, b, NULL))];
                if (u1 > u2) { long long tt = u1; u1 = u2; u2 = tt; }
                if (na < cap) adjkey[na++] = u1 * 1000000 + u2;
            }
        }
        long long nadj = 0;
        if (n <= 7) {
            /* count distinct unordered node pairs */
            for (long long i = 1; i < na; i++) { long long x = adjkey[i], j = i - 1; while (j >= 0 && adjkey[j] > x) { adjkey[j + 1] = adjkey[j]; j--; } adjkey[j + 1] = x; }
            for (long long i = 0; i < na; i++) if (!i || adjkey[i] != adjkey[i - 1]) nadj++;
        }
        printf("n=%d line classes (<kappa,sigma> orbits) = %lld; distinct node adjacencies (incl. loops) among them = %lld\n", n, orbits, nadj);
        free(seen); free(adjkey);
    }
    free(RL); free(RH); free(PhL); free(PhH); free(fib); free(blue);
}

/* canonical-form self test: invariance under random relabelling (all n), |Aut| by brute force (n <= 7) */
static void selftest(int n) {
    Ord *Q = &O[n];
    long long bad = 0;
    unsigned seed = 12345u + n;
    long long lim = Q->ncls < 20000 ? Q->ncls : 20000;
    for (long long s = 0; s < lim; s++) {
        long long i = Q->ncls <= 20000 ? s : (long long)(((u64)s * 2654435761ULL) % (u64)Q->ncls);
        u16 t[MAXN], r[MAXN];
        decode(n, Q->code[i], t);
        for (int rep = 0; rep < 3; rep++) {
            int perm[MAXN];
            for (int k = 0; k < n; k++) perm[k] = k;
            for (int k = n - 1; k > 0; k--) { seed = seed * 1103515245u + 12345u; int j = (seed >> 8) % (k + 1); int tt = perm[k]; perm[k] = perm[j]; perm[j] = tt; }
            for (int k = 0; k < n; k++) r[k] = 0;
            for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if ((t[a] >> b) & 1) r[perm[a]] |= 1u << perm[b];
            if (canon(n, r, NULL) != Q->code[i]) bad++;
        }
        if (n <= 7) {
            /* brute-force automorphism count */
            int perm[MAXN]; for (int k = 0; k < n; k++) perm[k] = k;
            long long cnt = 0;
            for (;;) {
                int ok = 1;
                for (int a = 0; a < n && ok; a++) for (int b = 0; b < n; b++) if (((t[a] >> b) & 1) != ((t[perm[a]] >> perm[b]) & 1)) { ok = 0; break; }
                cnt += ok;
                int k = n - 2; while (k >= 0 && perm[k] > perm[k + 1]) k--;
                if (k < 0) break;
                int l = n - 1; while (perm[l] < perm[k]) l--;
                int tt = perm[k]; perm[k] = perm[l]; perm[l] = tt;
                for (int a = k + 1, b = n - 1; a < b; a++, b--) { tt = perm[a]; perm[a] = perm[b]; perm[b] = tt; }
            }
            if (cnt != Q->aut[i]) bad++;
        }
    }
    printf("n=%d canonical-form self-test on %lld classes (3 random relabellings each%s): %lld failures\n", n, lim, n <= 7 ? ", brute-force |Aut|" : "", bad);
}

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 8;
    const char *dir = argc > 2 ? argv[2] : ".";
    int tiles = (argc > 4 && !strcmp(argv[3], "tiles")) ? atoi(argv[4]) : 0;
    double t0 = omp_get_wtime();
    printf("audit_tower NMAX=%d threads=%d\n", NMAX, omp_get_max_threads());
    for (int n = 1; n <= NMAX; n++) {
        generate(n);
        analyse(n);
        Ord *Q = &O[n];
        printf("n=%d classes=%lld self-converse=%lld merged=%lld (%.1f s)\n", n, Q->ncls, Q->nsc, Q->nm, omp_get_wtime() - t0);
        if (n >= 2) selftest(n);
        if (n >= 3) {
            collisions_report(n);
            if (n <= 9) {
                long long rk = rank_RtDR(n, 0x5EED0000ULL + n);
                printf("n=%d rank_p(R^T D R) = %lld of %lld columns -> %s\n", n, rk, O[n - 1].nm, rk == O[n - 1].nm ? "R_n has full column rank over Q" : "NOT certified");
            }
            if (n <= tiles) tiling_checks(n);
        }
        if (n <= 9) dump(n, dir);
        fflush(stdout);
    }
    printf("total %.1f s\n", omp_get_wtime() - t0);
    return 0;
}
