/* Audit E (2026-10-14), HYP-9279 (B_loc) and (Z) (added in commit 93231731e6), independent of tbal6/tbal7.
   Input: digraph6 lines (gentourng -q -z -c n).  DFS enumeration of all Hamiltonian paths (exact tallies), then
   (B_loc) every vertex lies on an arc x -> y with 1/2 < P(x before y) <= 2/3 that belongs to a directed 3-cycle;
   (Z)     for every z with T - z not strong, components D_1 => ... => D_k of T - z (topological order):
           (Z1) k >= 3: 1/3 <= P(v before z) <= 2/3 for all v in D_2 u ... u D_(k-1);  minimum of min(P,1-P) there;
           (Z2) z has a balanced partner; and one lies in D_1 u D_2 u D_(k-1) u D_k;  for k = 2 the minimum over instances of
                z's best-partner balance.                                                                                 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
static int n; static unsigned out[16]; static int path[16]; static u64 bef[16][16]; static u64 H;
static void dfs(int d, unsigned used) {
    if (d == n) { H++; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) bef[path[i]][path[j]]++; return; }
    unsigned c = out[path[d - 1]] & ~used;
    while (c) { int v = __builtin_ctz(c); c &= c - 1; path[d] = v; dfs(d + 1, used | 1u << v); }
}
static unsigned reach(unsigned S, int s) {          /* vertices reachable from s inside S */
    unsigned r = 1u << s, fr = r;
    while (fr) { int v = __builtin_ctz(fr); fr &= fr - 1; unsigned nw = out[v] & S & ~r; r |= nw; fr |= nw; }
    return r;
}
int main(void) {
    char line[256];
    long long tot = 0, bloc_ok = 0, inst = 0, inst_k3 = 0, z1_fail = 0, z2_fail = 0, z2_outer_fail = 0;
    u64 z1min_num = 1, z1min_den = 1, k2min_num = 1, k2min_den = 1; char k2arg[128] = ""; int k2z = -1; char k2sizes[64] = "";
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '&') continue;
        line[strcspn(line, "\r\n")] = 0;
        n = line[1] - 63; memset(out, 0, sizeof out);
        for (int k = 0; k < n * n; k++) { int b = line[2 + k / 6] - 63; if (b >> (5 - k % 6) & 1) out[k / n] |= 1u << (k % n); }
        unsigned full = (1u << n) - 1;
        if (reach(full, 0) != full) continue;
        { int sg = 1; for (int v = 0; v < n; v++) if (reach(full, v) != full) sg = 0; if (!sg) continue; }
        tot++;
        memset(bef, 0, sizeof bef); H = 0;
        for (int s = 0; s < n; s++) { path[0] = s; dfs(1, 1u << s); }
        /* (B_loc) */
        unsigned covered = 0;
        for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) if (out[x] >> y & 1) {
            u64 c = bef[x][y];
            if (!(2 * c > H && 3 * c <= 2 * H)) continue;
            int tri = 0; for (int w = 0; w < n; w++) if (w != x && w != y && (out[y] >> w & 1) && (out[w] >> x & 1)) tri = 1;
            if (tri) covered |= (1u << x) | (1u << y);
        }
        if (covered == full) bloc_ok++;
        /* (Z) */
        for (int z = 0; z < n; z++) {
            unsigned R = full ^ (1u << z);
            int strongz = 1;
            for (int v = 0; v < n; v++) if (R >> v & 1) if (reach(R, v) != R) { strongz = 0; break; }
            if (strongz) continue;
            inst++;
            /* components and their order: sort by size of reach set (descending) */
            int comp[16], nc = 0; unsigned cm[16]; unsigned rs[16];
            for (int v = 0; v < n; v++) comp[v] = -1;
            for (int v = 0; v < n; v++) if ((R >> v & 1) && comp[v] < 0) {
                unsigned rv = reach(R, v), C = 0;
                for (int u = 0; u < n; u++) if ((rv >> u & 1) && (reach(R, u) >> v & 1)) C |= 1u << u;
                for (int u = 0; u < n; u++) if (C >> u & 1) comp[u] = nc;
                cm[nc] = C; rs[nc] = rv; nc++;
            }
            /* order components: D_1 reaches everything (largest reach) */
            int ord[16]; for (int i = 0; i < nc; i++) ord[i] = i;
            for (int i = 0; i < nc; i++) for (int j = i + 1; j < nc; j++) if (__builtin_popcount(rs[ord[j]]) > __builtin_popcount(rs[ord[i]])) { int t = ord[i]; ord[i] = ord[j]; ord[j] = t; }
            int k = nc;
            if (k >= 3) {
                inst_k3++;
                for (int b = 1; b < k - 1; b++) { unsigned C = cm[ord[b]];
                    for (int v = 0; v < n; v++) if (C >> v & 1) {
                        u64 c = bef[v][z], m = c < H - c ? c : H - c;
                        if (3 * m < H) z1_fail++;
                        if ((u128)m * z1min_den < (u128)z1min_num * H) { z1min_num = m; z1min_den = H; }
                    } }
            }
            u64 bz = 0; int outer = 0;
            unsigned O = cm[ord[0]] | cm[ord[k - 1]] | (k >= 2 ? cm[ord[1]] | cm[ord[k - 2]] : 0);
            for (int v = 0; v < n; v++) if (v != z) {
                u64 c = bef[z][v], m = c < H - c ? c : H - c; if (m > bz) bz = m;
                if ((O >> v & 1) && 3 * m >= H) outer = 1;
            }
            if (3 * bz < H) z2_fail++;
            if (!outer) z2_outer_fail++;
            if (k == 2 && (u128)bz * k2min_den < (u128)k2min_num * H) { k2min_num = bz; k2min_den = H; strncpy(k2arg, line, 127); k2z = z;
                snprintf(k2sizes, 64, "|D1|=%d |D2|=%d", __builtin_popcount(cm[ord[0]]), __builtin_popcount(cm[ord[1]])); }
        }
    }
    u64 g1 = z1min_num, g2 = z1min_den; while (g2) { u64 t = g1 % g2; g1 = g2; g2 = t; }
    u64 h1 = k2min_num, h2 = k2min_den; while (h2) { u64 t = h1 % h2; h1 = h2; h2 = t; }
    printf("n=%d: strong %lld | (B_loc) holds in %lld | cut-vertex instances %lld (k >= 3: %lld) | (Z1) failures %lld, min balance of "
           "interior-block vertices vs z %llu/%llu | (Z2) failures %lld, partner not in D1uD2uD(k-1)uDk %lld | k = 2: min best-partner "
           "balance of z %llu/%llu\n", n, tot, bloc_ok, inst, inst_k3, z1_fail, z1min_num / g1, z1min_den / g1, z2_fail, z2_outer_fail,
           k2min_num / h1, k2min_den / h1);
    printf("   k = 2 argmin: %s z=%d %s\n", k2arg, k2z, k2sizes);
    return 0;
}
