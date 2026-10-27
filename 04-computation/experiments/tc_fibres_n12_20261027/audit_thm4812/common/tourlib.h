// tourlib.h -- independent audit helpers for THM-4812 (written from scratch, no code shared with the repository).
// Tournaments on n <= 32 vertices as out-neighbour bitmasks out[v] (bit w set iff v -> w).
#pragma once
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static inline int pc32(uint32_t x) { return __builtin_popcount(x); }
static inline int ctz32(uint32_t x) { return __builtin_ctz(x); }

static void make_in(int n, const uint32_t *out, uint32_t *in) {
    for (int v = 0; v < n; v++) in[v] = 0;
    for (int u = 0; u < n; u++) for (uint32_t o = out[u]; o; o &= o - 1) in[ctz32(o)] |= 1u << u;
}

// sanity: exactly one arc per pair, no loops
static int is_tournament(int n, const uint32_t *out) {
    for (int u = 0; u < n; u++) {
        if (out[u] >> u & 1u) return 0;
        for (int v = u + 1; v < n; v++) { int a = (out[u] >> v) & 1u, b = (out[v] >> u) & 1u; if (a + b != 1) return 0; }
    }
    return 1;
}

// Strong connectivity by Landau's criterion (independent of reachability):
// T is strong iff the ascending score sequence s_1 <= ... <= s_n has s_1 + ... + s_k > C(k,2) for every 1 <= k < n.
static int strong_landau(int n, const uint32_t *out) {
    int cnt[33] = {0};
    for (int v = 0; v < n; v++) cnt[pc32(out[v])]++;
    long sum = 0; int k = 0;
    for (int s = 0; s < n; s++) {
        for (int j = 0; j < cnt[s]; j++) {
            sum += s; k++;
            if (k < n && sum <= (long)k * (k - 1) / 2) return 0;
        }
    }
    return 1;
}

// Hamiltonian paths, pull-style subset DP: f[mask][v] = #paths with vertex set mask ending at v.
// Uses a caller-provided buffer of size (1<<n)*n uint64.
static uint64_t ham_pull(int n, const uint32_t *in, uint64_t *f) {
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    for (int v = 0; v < n; v++) f[((size_t)1 << v) * n + v] = 1;
    for (uint32_t mask = 1; mask <= full && mask != 0; mask++) {
        if (!(mask & (mask - 1))) continue;
        uint64_t *row = f + (size_t)mask * n;
        for (uint32_t vs = mask; vs; vs &= vs - 1) {
            int v = ctz32(vs); uint32_t pm = mask ^ (1u << v);
            const uint64_t *prow = f + (size_t)pm * n;
            uint64_t s = 0;
            for (uint32_t us = pm & in[v]; us; us &= us - 1) s += prow[ctz32(us)];
            row[v] = s;
        }
        if (mask == full) break;
    }
    uint64_t H = 0; for (int v = 0; v < n; v++) H += f[(size_t)full * n + v];
    return H;
}

// ---------------- automorphism group order (orbit-stabilizer over a backtracking existence test) ----------------
static uint64_t mix64(uint64_t x) { x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL; x ^= x >> 27; x *= 0x94d049bb133111ebULL; x ^= x >> 31; return x; }

// isomorphism-invariant colour refinement (hash of colour + multiset of out-neighbour colours)
static void refine_colours(int n, const uint32_t *out, uint64_t *col) {
    for (int v = 0; v < n; v++) col[v] = mix64(0x1234 + pc32(out[v]));
    int ncol = -1;
    for (int round = 0; round < n + 2; round++) {
        uint64_t nc[32];
        for (int v = 0; v < n; v++) {
            uint64_t s = 0, t = 0;
            for (uint32_t o = out[v]; o; o &= o - 1) s += mix64(col[ctz32(o)] ^ 0x9e3779b97f4a7c15ULL);
            for (int w = 0; w < n; w++) if (w != v && !((out[v] >> w) & 1u)) t += mix64(col[w] ^ 0x7f4a7c159e3779b9ULL);
            nc[v] = mix64(col[v] * 31 + s * 7 + t * 3 + 1);
        }
        // count distinct
        int d = 0; uint64_t seen[32];
        for (int v = 0; v < n; v++) { int f = 0; for (int j = 0; j < d; j++) if (seen[j] == nc[v]) { f = 1; break; } if (!f) seen[d++] = nc[v]; }
        memcpy(col, nc, sizeof(uint64_t) * n);
        if (d == ncol) break;
        ncol = d;
    }
}

static int ext_exists(int n, const uint32_t *out, const uint32_t *in, const uint32_t *cand_in, uint32_t unassigned) {
    if (!unassigned) return 1;
    int v = -1, bc = 99;
    for (uint32_t us = unassigned; us; us &= us - 1) { int w = ctz32(us); int c = pc32(cand_in[w]); if (c == 0) return 0; if (c < bc) { bc = c; v = w; } }
    uint32_t nu = unassigned & ~(1u << v);
    for (uint32_t xs = cand_in[v]; xs; xs &= xs - 1) {
        int x = ctz32(xs);
        uint32_t cand[32]; int ok = 1;
        for (uint32_t ws = nu; ws; ws &= ws - 1) {
            int w = ctz32(ws);
            uint32_t m = ((out[v] >> w) & 1u) ? out[x] : in[x];
            cand[w] = cand_in[w] & m & ~(1u << x);
            if (!cand[w]) { ok = 0; break; }
        }
        if (ok && ext_exists(n, out, in, cand, nu)) return 1;
    }
    return 0;
}

// propagate assignment v -> x into a candidate table; returns 0 on wipe-out
static int assign_prop(int n, const uint32_t *out, const uint32_t *in, const uint32_t *cand_in, uint32_t unassigned, int v, int x, uint32_t *cand_outp) {
    uint32_t nu = unassigned & ~(1u << v);
    for (int w = 0; w < n; w++) cand_outp[w] = cand_in[w];
    for (uint32_t ws = nu; ws; ws &= ws - 1) {
        int w = ctz32(ws);
        uint32_t m = ((out[v] >> w) & 1u) ? out[x] : in[x];
        cand_outp[w] = cand_in[w] & m & ~(1u << x);
        if (!cand_outp[w]) return 0;
    }
    cand_outp[v] = 1u << x;
    return 1;
}

static long long aut_order(int n, const uint32_t *out) {
    uint32_t in[32]; make_in(n, out, in);
    uint64_t col[32]; refine_colours(n, out, col);
    uint32_t cand[32];
    for (int v = 0; v < n; v++) { cand[v] = 0; for (int x = 0; x < n; x++) if (col[x] == col[v]) cand[v] |= 1u << x; }
    uint32_t unassigned = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    long long order = 1;
    while (unassigned) {
        // stop if every unassigned candidate set is a singleton (then only the identity extends)
        int b = -1, bc = 99, allsingle = 1;
        for (uint32_t us = unassigned; us; us &= us - 1) { int w = ctz32(us); int c = pc32(cand[w]); if (c > 1) allsingle = 0; if (c > 1 && c < bc) { bc = c; b = w; } }
        if (allsingle) break;
        long long cnt = 0;
        for (uint32_t xs = cand[b]; xs; xs &= xs - 1) {
            int x = ctz32(xs); uint32_t c2[32];
            if (assign_prop(n, out, in, cand, unassigned, b, x, c2) && ext_exists(n, out, in, c2, unassigned & ~(1u << b))) cnt++;
        }
        order *= cnt;
        uint32_t c3[32];
        if (!assign_prop(n, out, in, cand, unassigned, b, b, c3)) { fprintf(stderr, "aut_order: identity propagation failed\n"); exit(2); }
        memcpy(cand, c3, sizeof cand);
        unassigned &= ~(1u << b);
    }
    return order;
}

// brute-force automorphism count by full backtracking (for testing aut_order on small n)
static long long aut_brute_rec(int n, const uint32_t *out, int *img, uint32_t used, int i) {
    if (i == n) return 1;
    long long c = 0;
    for (int x = 0; x < n; x++) {
        if (used >> x & 1u) continue;
        int ok = 1;
        for (int j = 0; j < i && ok; j++) { int a = (out[j] >> i) & 1u, b = (out[img[j]] >> x) & 1u; if (a != b) ok = 0; }
        if (!ok) continue;
        img[i] = x; c += aut_brute_rec(n, out, img, used | (1u << x), i + 1);
    }
    return c;
}
static long long aut_brute(int n, const uint32_t *out) { int img[32]; return aut_brute_rec(n, out, img, 0, 0); }

// strong connectivity by BFS (second method, for cross-checking)
static int strong_bfs(int n, const uint32_t *out) {
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1), in[32]; make_in(n, out, in);
    uint32_t seen = 1, fr = 1;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= out[ctz32(f)]; nx &= ~seen; seen |= nx; fr = nx; }
    if (seen != full) return 0;
    seen = 1; fr = 1;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= in[ctz32(f)]; nx &= ~seen; seen |= nx; fr = nx; }
    return seen == full;
}
