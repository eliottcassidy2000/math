/* hampath.c (klein-2026-11-02) -- Hamiltonian path counts H(T) of tournaments.
   modes:  census n      : read gentourng's default format (upper triangle row by row, one tournament per line) on stdin,
                           print the spectrum (value -> number of tournaments) and the first tournament attaining each value
           paley q       : H of the Paley tournament P_q (q = 3 mod 4 prime), via vertex-transitivity: H = q * paths from 0
           climb n target seed : random arc flips on an n-tournament, hill-climbing |H - target| to 0; prints a witness
   The count uses the subset DP f[S][v] = number of paths covering S and ending at v (64-bit). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
static int n; static u64 adj[64]; /* adj[u] bit v set iff u -> v */
static u64 *f;                    /* f[S * n + v] */
static u64 count_paths(void) {   /* all Hamiltonian paths */
    u64 N = 1ULL << n; memset(f, 0, sizeof(u64) * N * n);
    for (int v = 0; v < n; v++) f[(1ULL << v) * n + v] = 1;
    for (u64 S = 1; S < N; S++) for (int v = 0; v < n; v++) {
        u64 c = f[S * n + v]; if (!c) continue;
        u64 out = adj[v] & ~S;
        while (out) { int w = __builtin_ctzll(out); out &= out - 1; f[(S | (1ULL << w)) * n + w] += c; }
    }
    u64 h = 0; for (int v = 0; v < n; v++) h += f[(N - 1) * n + v]; return h;
}
static u64 count_from(int s) {   /* Hamiltonian paths starting at s, only subsets containing s (memory n * 2^(n-1)) */
    /* relabel so that s is vertex n-1 and enumerate subsets of the other n-1 vertices plus s fixed */
    int m = n - 1; u64 N = 1ULL << m;                /* S ranges over subsets of the m other vertices */
    int map[64], k = 0; for (int v = 0; v < n; v++) if (v != s) map[k++] = v;
    memset(f, 0, sizeof(u64) * N * n);              /* f[S * n + v]: v in {0..m-1} index into map, or v = m for s itself */
    f[0 * n + m] = 1;                                /* empty set of others, ending at s */
    for (u64 S = 0; S < N; S++) {
        for (int vi = 0; vi <= m; vi++) {
            u64 c = f[S * n + vi]; if (!c) continue;
            int v = (vi == m) ? s : map[vi];
            for (int wi = 0; wi < m; wi++) if (!(S >> wi & 1) && (adj[v] >> map[wi] & 1)) f[(S | (1ULL << wi)) * n + wi] += c;
        }
    }
    u64 h = 0; for (int vi = 0; vi < m; vi++) h += f[(N - 1) * n + vi]; return h;
}

static u64 count_cycles(void) {   /* directed Hamiltonian cycles = paths from 0 covering all, ending at v with v -> 0 */
    int m = n - 1; u64 N = 1ULL << m;
    int map[64], k = 0; for (int v = 0; v < n; v++) if (v != 0) map[k++] = v;
    memset(f, 0, sizeof(u64) * N * n);
    f[0 * n + m] = 1;
    for (u64 S = 0; S < N; S++) for (int vi = 0; vi <= m; vi++) {
        u64 c = f[S * n + vi]; if (!c) continue;
        int v = (vi == m) ? 0 : map[vi];
        for (int wi = 0; wi < m; wi++) if (!(S >> wi & 1) && (adj[v] >> map[wi] & 1)) f[(S | (1ULL << wi)) * n + wi] += c;
    }
    u64 cyc = 0; for (int vi = 0; vi < m; vi++) if (adj[map[vi]] >> 0 & 1) cyc += f[(N - 1) * n + vi]; return cyc;
}
static int is_strong(void) {
    u64 seen = 1, frontier = 1;
    while (frontier) { u64 nf = 0; for (int v = 0; v < n; v++) if (frontier >> v & 1) nf |= adj[v]; nf &= ~seen; seen |= nf; frontier = nf; }
    if (seen != (1ULL << n) - 1) return 0;
    u64 radj[64]; for (int v = 0; v < n; v++) radj[v] = 0;
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (adj[u] >> v & 1) radj[v] |= 1ULL << u;
    seen = 1; frontier = 1;
    while (frontier) { u64 nf = 0; for (int v = 0; v < n; v++) if (frontier >> v & 1) nf |= radj[v]; nf &= ~seen; seen |= nf; frontier = nf; }
    return seen == (1ULL << n) - 1;
}
static int isqr(int a, int q) { for (int x = 1; x < q; x++) if ((x * x) % q == a) return 1; return 0; }
static u64 xorshift(u64 *s) { u64 x = *s; x ^= x << 13; x ^= x >> 7; x ^= x << 17; *s = x; return x; }
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: hampath census n | paley q | climb n target seed\n"); return 1; }
    if (!strcmp(argv[1], "census")) {
        n = atoi(argv[2]); u64 N = 1ULL << n; f = malloc(sizeof(u64) * N * n);
        char line[4096]; u64 *spec = calloc(1 << 20, sizeof(u64)); char **first = calloc(1 << 20, sizeof(char *));
        u64 total = 0, maxh = 0;
        while (fgets(line, sizeof line, stdin)) {
            int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
            if (L != n * (n - 1) / 2) continue;
            memset(adj, 0, sizeof adj); int idx = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
            u64 h = count_paths(); total++;
            if (h < (1 << 20)) { if (!spec[h]) first[h] = strdup(line); spec[h]++; } if (h > maxh) maxh = h;
        }
        printf("order %d: %llu tournaments, max H = %llu\n", n, total, maxh);
        printf("spectrum (H: count, first witness):\n");
        for (u64 h = 0; h < (1 << 20); h++) if (spec[h]) printf("  %llu: %llu  %s\n", h, spec[h], first[h]);
        return 0;
    }
    if (!strcmp(argv[1], "census2")) {   /* per line: H C strong code */
        n = atoi(argv[2]); u64 N = 1ULL << n; f = malloc(sizeof(u64) * N * n);
        char line[4096];
        while (fgets(line, sizeof line, stdin)) {
            int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
            if (L != n * (n - 1) / 2) continue;
            memset(adj, 0, sizeof adj); int idx = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
            u64 h = count_paths(); u64 c = count_cycles(); int st = is_strong();
            printf("%llu %llu %d %s\n", h, c, st, line);
        }
        return 0;
    }
    if (!strcmp(argv[1], "starts")) {   /* per line: H then h_start(v) for each vertex v */
        n = atoi(argv[2]); u64 N = 1ULL << (n - 1); f = malloc(sizeof(u64) * (N << 1) * n);
        char line[4096];
        while (fgets(line, sizeof line, stdin)) {
            int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
            if (L != n * (n - 1) / 2) continue;
            memset(adj, 0, sizeof adj); int idx = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
            u64 h = count_paths(); printf("%llu", h);
            for (int v = 0; v < n; v++) printf(" %llu", count_from(v));
            printf(" %s\n", line);
        }
        return 0;
    }
    if (!strcmp(argv[1], "paley")) {
        int q = atoi(argv[2]); n = q; memset(adj, 0, sizeof adj);
        for (int i = 0; i < q; i++) for (int j = 0; j < q; j++) if (i != j && isqr(((j - i) % q + q) % q, q)) adj[i] |= 1ULL << j;
        u64 N = 1ULL << (n - 1); f = malloc(sizeof(u64) * N * n); if (!f) { fprintf(stderr, "alloc failed\n"); return 2; }
        u64 h0 = count_from(0); u64 H = h0 * (u64)q; u64 aut = (u64)q * (q - 1) / 2;
        printf("Paley q = %d: paths from 0 = %llu, H = q * that = %llu, |Aut| = q(q-1)/2 = %llu, H / |Aut| = %llu (remainder %llu)\n", q, h0, H, aut, H / aut, H % aut);
        return 0;
    }
    if (!strcmp(argv[1], "climb")) {
        n = atoi(argv[2]); long long target = atoll(argv[3]); u64 seed = atoll(argv[4]) * 2654435761ULL + 1;
        u64 N = 1ULL << n; f = malloc(sizeof(u64) * N * n);
        memset(adj, 0, sizeof adj);
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (xorshift(&seed) & 1) adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
        long long h = (long long)count_paths(); long long best = llabs(h - target); int it = 0;
        while (h != target && it < 200000) {
            int i = xorshift(&seed) % n, j = xorshift(&seed) % n; if (i == j) continue; it++;
            /* flip arc i-j */
            if (adj[i] >> j & 1) { adj[i] &= ~(1ULL << j); adj[j] |= 1ULL << i; } else { adj[j] &= ~(1ULL << i); adj[i] |= 1ULL << j; }
            long long h2 = (long long)count_paths(); long long d = llabs(h2 - target);
            if (d <= best || (xorshift(&seed) % 100) < 3) { h = h2; best = d; }
            else { if (adj[i] >> j & 1) { adj[i] &= ~(1ULL << j); adj[j] |= 1ULL << i; } else { adj[j] &= ~(1ULL << i); adj[i] |= 1ULL << j; } }
        }
        printf("order %d target %lld: %s after %d flips, H = %lld\n", n, target, h == target ? "FOUND" : "not found", it, h);
        if (h == target) { printf("witness (upper triangle, row by row): "); for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) putchar(adj[i] >> j & 1 ? '1' : '0'); putchar('\n'); }
        return h == target ? 0 : 1;
    }
    return 1;
}
