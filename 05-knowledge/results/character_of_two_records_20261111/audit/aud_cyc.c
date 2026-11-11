/* AUDIT program (independent): directed cycles by length, cycles through each vertex by length, and
 * directed paths by number of vertices, of a tournament (or digraph) on n <= 24 vertices.
 * All counts are modulo 2^32 (exact when below 2^32), except mode dfs (exact, 64 bits).
 *
 * usage: aud_cyc MODE < input          input: n, then n rows of n characters 0/1 (row i, column j: 1 iff i -> j)
 *   MODE = global : c_l for l = 3..n; a cycle is counted at its smallest vertex r: paths r -> ... -> v inside
 *                   {r+1, ..., n-1} by dynamic programming over (set, last vertex), closed by an arc v -> r.
 *   MODE = vertex : for every vertex v, the number of cycles of each length through v (the same programme
 *                   rooted at v on all the other vertices). Then l * c_l = sum over v, an independent check.
 *   MODE = paths  : p_j = number of directed paths on j vertices, j = 1..n (n <= 23).
 *   MODE = dfs    : every cycle enumerated by depth-first search from its smallest vertex (n <= 16).
 * Output: one line per length: "c l value" / "v vertex l value" / "p j value" / "d l value".
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int n;
static uint32_t outm[32];
static int adj[32][32];

/* cycles through the vertex `root`, using only the vertices listed in vs[0..m-1] besides the root */
static void rooted(int root, const int *vs, int m, uint32_t *dp, const uint32_t *off, uint64_t *cyc) {
    uint32_t o[32], start = 0, back = 0, full = (m == 32) ? 0xffffffffu : ((1u << m) - 1), mask;
    int i, j;
    for (i = 0; i < m; i++) {
        o[i] = 0;
        for (j = 0; j < m; j++) if (adj[vs[i]][vs[j]]) o[i] |= 1u << j;
        if (adj[root][vs[i]]) start |= 1u << i;
        if (adj[vs[i]][root]) back |= 1u << i;
    }
    memset(dp, 0, (size_t)off[full + 1] * sizeof(uint32_t));
    for (i = 0; i < m; i++) if (start >> i & 1) dp[off[1u << i]] = 1;
    for (mask = 1; mask <= full && mask != 0; mask++) {
        uint32_t rem = ~mask & full, mm, idx = 0;
        const uint32_t *rowp = dp + off[mask];
        int k = __builtin_popcount(mask);
        for (mm = mask; mm; mm &= mm - 1, idx++) {
            uint32_t x = rowp[idx];
            if (!x) continue;
            int u = __builtin_ctz(mm);
            if (k >= 2 && (back >> u & 1)) cyc[k + 1] += x;
            uint32_t t;
            for (t = o[u] & rem; t; t &= t - 1) {
                int w = __builtin_ctz(t);
                uint32_t nm = mask | (1u << w);
                dp[off[nm] + (uint32_t)__builtin_popcount(nm & ((1u << w) - 1))] += x;
            }
        }
        if (mask == full) break;
    }
}

static uint64_t dcyc[40];
static int path[32], onp[32], droot;
static void dfs(int len) {
    int u = path[len - 1], w;
    if (len >= 3 && adj[u][droot]) dcyc[len]++;
    for (w = droot + 1; w < n; w++)
        if (!onp[w] && adj[u][w]) { onp[w] = 1; path[len] = w; dfs(len + 1); onp[w] = 0; }
}

int main(int argc, char **argv) {
    const char *mode = argc > 1 ? argv[1] : "global";
    int i, j, l;
    char row[64];
    if (scanf("%d", &n) != 1 || n < 3 || n > 24) { fprintf(stderr, "bad n\n"); return 1; }
    for (i = 0; i < n; i++) {
        if (scanf("%63s", row) != 1 || (int)strlen(row) != n) { fprintf(stderr, "bad row\n"); return 1; }
        outm[i] = 0;
        for (j = 0; j < n; j++) { adj[i][j] = row[j] == '1'; if (adj[i][j]) outm[i] |= 1u << j; }
    }
    if (!strcmp(mode, "dfs")) {
        if (n > 16) { fprintf(stderr, "dfs: n <= 16\n"); return 1; }
        for (droot = 0; droot < n; droot++) { memset(onp, 0, sizeof onp); path[0] = droot; onp[droot] = 1; dfs(1); }
        for (l = 3; l <= n; l++) printf("d %d %llu\n", l, (unsigned long long)dcyc[l]);
        return 0;
    }
    int mbits = !strcmp(mode, "paths") ? n : n - 1;
    if (mbits > 23) { fprintf(stderr, "too large for this mode\n"); return 1; }
    uint32_t NM = 1u << mbits, mask;
    uint32_t *off = malloc((size_t)(NM + 1) * sizeof(uint32_t));
    if (!off) { fprintf(stderr, "no memory\n"); return 1; }
    off[0] = 0;
    for (mask = 0; mask < NM; mask++) off[mask + 1] = off[mask] + (uint32_t)__builtin_popcount(mask);
    uint32_t *dp = malloc((size_t)off[NM] * sizeof(uint32_t));
    if (!dp) { fprintf(stderr, "no memory\n"); return 1; }
    if (!strcmp(mode, "global")) {
        uint64_t cyc[40];
        memset(cyc, 0, sizeof cyc);
        int r, vs[32];
        for (r = 0; r <= n - 3; r++) {
            int m = n - 1 - r;
            for (i = 0; i < m; i++) vs[i] = r + 1 + i;
            rooted(r, vs, m, dp, off, cyc);
        }
        for (l = 3; l <= n; l++) printf("c %d %u\n", l, (uint32_t)cyc[l]);
    } else if (!strcmp(mode, "vertex")) {
        int v, vs[32];
        for (v = 0; v < n; v++) {
            uint64_t cyc[40];
            memset(cyc, 0, sizeof cyc);
            int m = 0;
            for (i = 0; i < n; i++) if (i != v) vs[m++] = i;
            rooted(v, vs, m, dp, off, cyc);
            for (l = 3; l <= n; l++) printf("v %d %d %u\n", v, l, (uint32_t)cyc[l]);
        }
    } else if (!strcmp(mode, "paths")) {
        uint64_t p[40];
        memset(p, 0, sizeof p);
        memset(dp, 0, (size_t)off[NM] * sizeof(uint32_t));
        for (i = 0; i < n; i++) dp[off[1u << i]] = 1;
        for (mask = 1; mask < NM; mask++) {
            uint32_t rem = ~mask & (NM - 1), mm, idx = 0;
            const uint32_t *rowp = dp + off[mask];
            int k = __builtin_popcount(mask);
            for (mm = mask; mm; mm &= mm - 1, idx++) {
                uint32_t x = rowp[idx];
                if (!x) continue;
                p[k] += x;
                int u = __builtin_ctz(mm);
                uint32_t t;
                for (t = outm[u] & rem; t; t &= t - 1) {
                    int w = __builtin_ctz(t);
                    uint32_t nm = mask | (1u << w);
                    dp[off[nm] + (uint32_t)__builtin_popcount(nm & ((1u << w) - 1))] += x;
                }
            }
        }
        for (l = 1; l <= n; l++) printf("p %d %u\n", l, (uint32_t)p[l]);
    } else { fprintf(stderr, "unknown mode\n"); return 1; }
    return 0;
}
