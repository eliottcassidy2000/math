/* at.h -- auditor's own tournament toolkit (independent of classification/tour.hpp).
 * Plain C99.  Vertices 0..n-1, n <= 24.  o[v] = out-neighbourhood bitmask.
 * Conventions follow PROOF.md: f(X) = e(X -> V\X); lambda = min f over proper nonempty X;
 * delta = min(min outdegree, min indegree); Q = (q_1..q_n) Hamiltonian path; Q blue iff the positional
 * reversal r_Q (q_i -> q_(n+1-i)) is an anti-automorphism; Q is a descent candidate iff some minimum
 * out-cut X has q_1 in X and q_n not in X.
 */
#ifndef AT_H
#define AT_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXV 24
typedef struct { int n; uint32_t o[MAXV]; int s[MAXV]; } tour;

static inline int pc(uint32_t x) { return __builtin_popcount(x); }
static inline int lowbit(uint32_t x) { return __builtin_ctz(x); }
static inline uint32_t fullmask(int n) { return (n >= 32) ? 0xffffffffu : ((1u << n) - 1u); }

static void t_clear(tour *T, int n) { memset(T, 0, sizeof *T); T->n = n; }
static void t_arc(tour *T, int u, int v) { T->o[u] |= 1u << v; T->o[v] &= ~(1u << u); } /* force u -> v */
static void t_fin(tour *T) { for (int v = 0; v < T->n; v++) T->s[v] = pc(T->o[v]); }
static inline int t_has(const tour *T, int u, int v) { return (T->o[u] >> v) & 1u; }
static inline uint32_t t_in(const tour *T, int v) { return fullmask(T->n) & ~T->o[v] & ~(1u << v); }

/* sanity: every pair has exactly one arc */
static int t_valid(const tour *T) {
    for (int u = 0; u < T->n; u++) {
        if (T->o[u] >> u & 1u) return 0;
        for (int v = u + 1; v < T->n; v++) if (t_has(T, u, v) == t_has(T, v, u)) return 0;
    }
    return 1;
}

/* gentourng default output: upper triangle row by row, '1' at (i,j), i<j, means i -> j */
static int t_parse(const char *line, int n, tour *T) {
    t_clear(T, n);
    int k = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        char c = line[k++];
        if (c == '1') T->o[i] |= 1u << j; else if (c == '0') T->o[j] |= 1u << i; else return 0;
    }
    t_fin(T);
    return 1;
}
static void t_code(const tour *T, char *buf) {
    int k = 0;
    for (int i = 0; i < T->n; i++) for (int j = i + 1; j < T->n; j++) buf[k++] = t_has(T, i, j) ? '1' : '0';
    buf[k] = 0;
}

static void t_converse(const tour *T, tour *C) {
    t_clear(C, T->n);
    for (int u = 0; u < T->n; u++) C->o[u] = t_in(T, u);
    t_fin(C);
}

static int t_delta(const tour *T) {
    int d = 99;
    for (int v = 0; v < T->n; v++) { int a = T->s[v], b = T->n - 1 - T->s[v]; if (a < d) d = a; if (b < d) d = b; }
    return d;
}
static int t_regular(const tour *T) { for (int v = 0; v < T->n; v++) if (2 * T->s[v] != T->n - 1) return 0; return 1; }

/* f(X) for one X */
static inline int t_f(const tour *T, uint32_t X) {
    int f = 0; uint32_t c = ~X & fullmask(T->n);
    for (uint32_t r = X; r; r &= r - 1) f += pc(T->o[lowbit(r)] & c);
    return f;
}

/* strong connectivity of T[mask] by forward and backward reachability */
static int t_strong_mask(const tour *T, uint32_t mask) {
    if (!mask) return 0;
    int st = lowbit(mask);
    uint32_t R = 1u << st, fr = R;
    while (fr) { int x = lowbit(fr); fr &= fr - 1; uint32_t nb = T->o[x] & mask & ~R; R |= nb; fr |= nb; }
    if (R != mask) return 0;
    R = 1u << st; fr = R;
    while (fr) { int x = lowbit(fr); fr &= fr - 1; uint32_t nb = t_in(T, x) & mask & ~R; R |= nb; fr |= nb; }
    return R == mask;
}

/* all minimum out-cuts; returns lambda; cuts written to arr (cap) and count to *nc */
static int t_mincuts(const tour *T, uint32_t *arr, int cap, int *nc) {
    uint32_t F = fullmask(T->n); int best = 1 << 30; int k = 0;
    for (uint32_t X = 1; X < F; X++) {
        int f = t_f(T, X);
        if (f < best) { best = f; k = 0; }
        if (f == best) { if (k < cap) arr[k] = X; k++; }
    }
    *nc = k; return best;
}

/* positional reversal of Q is an anti-automorphism? */
static int t_blue(const tour *T, const int *Q) {
    int n = T->n, r[MAXV];
    for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) {
        /* image of out(u) under r must equal in(r(u)) */
        uint32_t img = 0;
        for (uint32_t c = T->o[u]; c; c &= c - 1) img |= 1u << r[lowbit(c)];
        if (img != t_in(T, r[u])) return 0;
    }
    return 1;
}

/* is sigma (a permutation) an anti-automorphism / automorphism? */
static int t_is_anti(const tour *T, const int *sg) {
    for (int u = 0; u < T->n; u++) for (int v = 0; v < T->n; v++) if (u != v && t_has(T, u, v) && !t_has(T, sg[v], sg[u])) return 0;
    return 1;
}
static int t_is_aut(const tour *T, const int *sg) {
    for (int u = 0; u < T->n; u++) for (int v = 0; v < T->n; v++) if (u != v && t_has(T, u, v) && !t_has(T, sg[u], sg[v])) return 0;
    return 1;
}

/* the native complement K_Q(T): keep the arcs of Q, reverse every other pair */
static void t_move(const tour *T, const int *Q, tour *K) {
    int n = T->n; uint32_t F = fullmask(n);
    t_clear(K, n);
    for (int u = 0; u < n; u++) K->o[u] = F & ~T->o[u] & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) t_arc(K, Q[i], Q[i + 1]);
    t_fin(K);
}

/* Is Q a Hamiltonian path of T? */
static int t_is_hpath(const tour *T, const int *Q) {
    uint32_t seen = 0;
    for (int i = 0; i < T->n; i++) { if (Q[i] < 0 || Q[i] >= T->n || (seen >> Q[i] & 1u)) return 0; seen |= 1u << Q[i]; }
    for (int i = 0; i + 1 < T->n; i++) if (!t_has(T, Q[i], Q[i + 1])) return 0;
    return 1;
}

/* Enumerate all Hamiltonian paths (plain DFS from every start).  cb returns nonzero to stop.  Returns 1 if stopped. */
typedef int (*hp_cb)(const int *Q, void *ctx);
static int t_hpaths_rec(const tour *T, int *Q, int d, uint32_t used, hp_cb cb, void *ctx) {
    if (d == T->n) return cb(Q, ctx);
    for (uint32_t c = T->o[Q[d - 1]] & ~used; c; c &= c - 1) {
        int w = lowbit(c); Q[d] = w;
        if (t_hpaths_rec(T, Q, d + 1, used | (1u << w), cb, ctx)) return 1;
    }
    return 0;
}
static int t_hpaths(const tour *T, hp_cb cb, void *ctx) {
    int Q[MAXV];
    for (int a = 0; a < T->n; a++) { Q[0] = a; if (t_hpaths_rec(T, Q, 1, 1u << a, cb, ctx)) return 1; }
    return 0;
}

/* Hamiltonian path of T[mask] with prescribed start (start<0: any) and optional end set; plain DFS, returns length or 0 */
static int t_hp_in_rec(const tour *T, uint32_t mask, uint32_t endset, int *P, int d, int k, uint32_t used) {
    if (d == k) return (endset >> P[d - 1]) & 1u;
    for (uint32_t c = T->o[P[d - 1]] & mask & ~used; c; c &= c - 1) {
        int w = lowbit(c); P[d] = w;
        if (t_hp_in_rec(T, mask, endset, P, d + 1, k, used | (1u << w))) return 1;
    }
    return 0;
}
static int t_hp_in(const tour *T, uint32_t mask, int start, uint32_t endset, int *P) {
    int k = pc(mask);
    for (uint32_t c = mask; c; c &= c - 1) {
        int a = lowbit(c); if (start >= 0 && a != start) continue;
        P[0] = a;
        if (t_hp_in_rec(T, mask, endset, P, 1, k, 1u << a)) return k;
    }
    return 0;
}

/* Count Hamiltonian paths of T[mask] ending in endset (bitmask DP, exact). */
static long long t_count_hp_end(const tour *T, uint32_t mask, uint32_t endset) {
    int n = T->n; size_t N = (size_t)1 << n;
    long long *dp = (long long *)calloc(N * n, sizeof(long long));
    for (int v = 0; v < n; v++) if (mask >> v & 1u) dp[((size_t)(1u << v)) * n + v] = 1;
    long long tot = 0;
    for (uint32_t S = 1; S < N; S++) {
        if ((S & mask) != S) continue;
        for (int v = 0; v < n; v++) {
            long long c = dp[(size_t)S * n + v]; if (!c) continue;
            if (S == mask) { if (endset >> v & 1u) tot += c; continue; }
            for (uint32_t w = T->o[v] & mask & ~S; w; w &= w - 1) { int x = lowbit(w); dp[(size_t)(S | (1u << x)) * n + x] += c; }
        }
    }
    free(dp);
    return tot;
}

/* Own backtracking enumeration of all automorphisms of T preserving the colouring col (col may be NULL).
 * cb(perm, ctx) is called for each automorphism (may be NULL).  Returns the number of automorphisms. */
typedef void (*aut_cb)(const int *perm, void *ctx);
static long long t_auts_rec(const tour *T, const int *col, int *img, uint32_t usedimg, int k, aut_cb cb, void *ctx) {
    int n = T->n;
    if (k == n) { if (cb) cb(img, ctx); return 1; }
    long long cnt = 0;
    for (int c = 0; c < n; c++) {
        if (usedimg >> c & 1u) continue;
        if (T->s[c] != T->s[k]) continue;
        if (col && col[c] != col[k]) continue;
        int ok = 1;
        for (int j = 0; j < k && ok; j++) if (t_has(T, k, j) != t_has(T, c, img[j])) ok = 0;
        if (!ok) continue;
        img[k] = c;
        cnt += t_auts_rec(T, col, img, usedimg | (1u << c), k + 1, cb, ctx);
    }
    return cnt;
}
static long long t_auts(const tour *T, const int *col, aut_cb cb, void *ctx) { int img[MAXV]; return t_auts_rec(T, col, img, 0, 0, cb, ctx); }

/* Own backtracking search for an anti-automorphism; returns 1 and fills sg if found */
static int t_anti_rec(const tour *T, int *sg, uint32_t used, int k) {
    int n = T->n;
    if (k == n) return 1;
    for (int c = 0; c < n; c++) {
        if (used >> c & 1u) continue;
        if (T->s[c] != n - 1 - T->s[k]) continue;
        int ok = 1;
        for (int j = 0; j < k && ok; j++) if (t_has(T, k, j) != t_has(T, sg[j], c)) ok = 0;
        if (!ok) continue;
        sg[k] = c;
        if (t_anti_rec(T, sg, used | (1u << c), k + 1)) return 1;
    }
    return 0;
}
static int t_find_anti(const tour *T, int *sg) { return t_anti_rec(T, sg, 0, 0); }
#endif
