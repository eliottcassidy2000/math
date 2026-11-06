/* hflip.c (klein-2026-11-06) -- the flip gradient of the Hamiltonian-path count on the tournament metagraph.
   usage: hflip n < gentourng-output
   For every tournament T read (gentourng default format) and every arc u -> v, let T' be T with that arc reversed.
   Prints, per class, one line:  "F h c3 strong  d1 d2 ... dm"  where m = n(n-1)/2 and d_i = h(T') - h(T) for the
   i-th arc (upper-triangle order), followed by the c3 change in the same order as "G e1 ... em"  (e_i = c3(T') - c3(T)).
   Also prints aggregate lines "# ...": the distribution of the gradient sign, the maximum |d| by level, and the
   union-find components of the attained values under the relation "some flip joins h and h'". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef unsigned int u32; typedef u32 cell;
#define HMAX (1u << 18)
static int n; static u64 adj[32], radj[32]; static cell *f;
static u64 count_paths(void) {
    u64 N = 1ULL << n;
    for (int v = 0; v < n; v++) f[(1ULL << v) * n + v] = 1;
    for (u64 S = 1; S < N; S++) {
        if (!(S & (S - 1))) continue;
        u64 vs = S; cell *dst = f + S * n;
        while (vs) { int v = __builtin_ctzll(vs); vs &= vs - 1; u64 R = S & ~(1ULL << v); u64 in = radj[v] & R; u64 c = 0; const cell *row = f + R * n;
            while (in) { int u = __builtin_ctzll(in); in &= in - 1; c += row[u]; } dst[v] = c; }
    }
    u64 h = 0; for (int v = 0; v < n; v++) h += f[(N - 1) * n + v]; return h;
}
static int c3_of(void) { int c = n * (n - 1) * (n - 2) / 6; for (int v = 0; v < n; v++) { int s = __builtin_popcountll(adj[v]); c -= s * (s - 1) / 2; } return c; }
static int is_strong(void) {
    u64 all = (1ULL << n) - 1, seen = 1, fr = 1;
    while (fr) { u64 nf = 0, x = fr; while (x) { int v = __builtin_ctzll(x); x &= x - 1; nf |= adj[v]; } nf &= ~seen; seen |= nf; fr = nf; }
    if (seen != all) return 0; seen = 1; fr = 1;
    while (fr) { u64 nf = 0, x = fr; while (x) { int v = __builtin_ctzll(x); x &= x - 1; nf |= radj[v]; } nf &= ~seen; seen |= nf; fr = nf; }
    return seen == all;
}
static int *uf; static int find(int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }
static void unite(int a, int b) { a = find(a); b = find(b); if (a != b) uf[a] = b; }
int main(int argc, char **argv) {
    if (argc < 2) return 1; n = atoi(argv[1]); u64 N = 1ULL << n; f = malloc(sizeof(cell) * N * n);
    uf = malloc(sizeof(int) * HMAX); for (u32 i = 0; i < HMAX; i++) uf[i] = i;
    unsigned char *seen = calloc(HMAX, 1);
    char line[1024]; int L0 = n * (n - 1) / 2; int m = L0;
    long long up = 0, down = 0, zero = 0, flips = 0; long long maxabs[64] = {0}; long long maxup[64] = {0}, maxdown[64] = {0};
    long long levelflips[64] = {0};
    long long *d = malloc(sizeof(long long) * m); int *e = malloc(sizeof(int) * m);
    while (fgets(line, sizeof line, stdin)) {
        int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0; if (L != L0) continue;
        memset(adj, 0, sizeof adj); memset(radj, 0, sizeof radj); int idx = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') { adj[i] |= 1ULL << j; radj[j] |= 1ULL << i; } else { adj[j] |= 1ULL << i; radj[i] |= 1ULL << j; } }
        u64 h = count_paths(); int c3 = c3_of(); int st = is_strong(); seen[h] = 1;
        idx = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            int u, v; if (adj[i] >> j & 1) { u = i; v = j; } else { u = j; v = i; }
            adj[u] &= ~(1ULL << v); radj[v] &= ~(1ULL << u); adj[v] |= 1ULL << u; radj[u] |= 1ULL << v;
            u64 h2 = count_paths(); int c32 = c3_of();
            adj[v] &= ~(1ULL << u); radj[u] &= ~(1ULL << v); adj[u] |= 1ULL << v; radj[v] |= 1ULL << u;
            d[idx] = (long long)h2 - (long long)h; e[idx] = c32 - c3; idx++;
            flips++; if (h2 > h) up++; else if (h2 < h) down++; else zero++;
            long long a = llabs((long long)h2 - (long long)h); if (a > maxabs[c3]) maxabs[c3] = a;
            if ((long long)h2 - (long long)h > maxup[c3]) maxup[c3] = (long long)h2 - (long long)h;
            if ((long long)h - (long long)h2 > maxdown[c3]) maxdown[c3] = (long long)h - (long long)h2;
            if (c32 == c3) levelflips[c3]++;
            unite((int)h, (int)h2); seen[h2] = 1;
        }
        printf("F %llu %d %d ", h, c3, st); for (int i = 0; i < m; i++) printf("%lld ", d[i]); printf("\nG"); for (int i = 0; i < m; i++) printf(" %d", e[i]); printf("\n");
    }
    printf("# order %d: %lld flips (class, arc): h up %lld, down %lld, unchanged %lld\n", n, flips, up, down, zero);
    printf("# by level k: max |dh|, max up, max down, level-preserving flips:\n");
    for (int k = 0; k < 64; k++) if (maxabs[k] || levelflips[k]) printf("#   k=%d: %lld %lld %lld %lld\n", k, maxabs[k], maxup[k], maxdown[k], levelflips[k]);
    int comps = 0; for (u32 h = 0; h < HMAX; h++) if (seen[h] && find(h) == (int)h) comps++;
    printf("# values joined by flips: components of the attained-value graph = %d\n", comps);
    return 0;
}
