/* Audit E (2026-10-14), HYP-9279 (B), (B+), sharpness and refuted variant (H); independent of tbal2-tbal6.
   Input: digraph6 lines (gentourng -q -z -c n).  For each strong tournament: all Hamiltonian paths by DFS (exact
   tallies bef[x][y]), then
     - vertex balance b(v) = max_w min(P, 1-P); tight tournaments (min_v b(v) = 1/3 exactly) and, for each tight vertex v,
       whether v lies on a cyclic triangle {v,a,b} such that V minus the triangle is transitive and exactly one arc leaves
       the triangle towards the rest ("C_3 hung at the end of the order by one back arc");
     - the balanced graph (edges = pairs with 1/3 <= P <= 2/3): connected? any vertex of degree 1?
     - Hamiltonian cycles (paths from vertex 0 that close): number, number with no balanced arc, minimum number of
       balanced arcs per cycle and how many tournaments attain it;
     - if n <= 8 (argument "short"): all cycles of length 3..n-1 (DFS from their least vertex) with no balanced arc.
   Usage: gentourng -q -z -c n | ./strong_dfs_extra [short]                                                       */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
static int n; static unsigned out[16]; static int path[16]; static u64 bef[16][16]; static u64 H;
static u64 hamcyc, hamcyc_nobal; static int hc_minbal; static int balmat[16][16];
static void dfs(int d, unsigned used) {
    if (d == n) {
        H++;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) bef[path[i]][path[j]]++;
        return;
    }
    unsigned c = out[path[d - 1]] & ~used;
    while (c) { int v = __builtin_ctz(c); c &= c - 1; path[d] = v; dfs(d + 1, used | 1u << v); }
}
/* Hamiltonian cycles: paths from 0 that close, evaluated after bef is known */
static void dfs_hc(int d, unsigned used) {
    if (d == n) {
        if (out[path[n - 1]] & 1u) {
            hamcyc++; int nb = 0;
            for (int i = 0; i < n; i++) { int a = path[i], b = path[(i + 1) % n]; nb += balmat[a][b]; }
            if (!nb) hamcyc_nobal++;
            if (nb < hc_minbal) hc_minbal = nb;
        }
        return;
    }
    unsigned c = out[path[d - 1]] & ~used;
    while (c) { int v = __builtin_ctz(c); c &= c - 1; path[d] = v; dfs_hc(d + 1, used | 1u << v); }
}
/* short cycles: start s = least vertex, other vertices > s */
static u64 shortcyc[16], shortcyc_nobal[16];
static void dfs_sc(int s, int d, unsigned used) {
    int last = path[d - 1];
    if (d >= 3 && (out[last] >> s & 1) && d < n) {
        shortcyc[d]++; int nb = 0;
        for (int i = 0; i < d; i++) nb += balmat[path[i]][path[(i + 1) % d]];
        if (!nb) shortcyc_nobal[d]++;
    }
    if (d == n) return;
    unsigned c = out[last] & ~used & ~((2u << s) - 1);
    while (c) { int v = __builtin_ctz(c); c &= c - 1; path[d] = v; dfs_sc(s, d + 1, used | 1u << v); }
}
static int strong(void) {
    unsigned R[16];
    for (int i = 0; i < n; i++) R[i] = out[i] | 1u << i;
    for (int k = 0; k < n; k++) for (int i = 0; i < n; i++) if (R[i] >> k & 1) R[i] |= R[k];
    for (int i = 0; i < n; i++) if (R[i] != (1u << n) - 1) return 0;
    return 1;
}
static int transitive_on(unsigned S) {
    unsigned seen = 0;
    for (int v = 0; v < n; v++) if (S >> v & 1) { int d = __builtin_popcount(out[v] & S); if (seen >> d & 1) return 0; seen |= 1u << d; }
    return 1;
}
int main(int argc, char **argv) {
    int doshort = argc > 1 && !strcmp(argv[1], "short");
    char line[256];
    long long tot = 0, tight = 0, tight_struct = 0, conn = 0, deg1 = 0, hc_tour_nobal = 0, minbal_count = 0;
    u64 hc_total = 0, hc_nobal_total = 0; int global_minbal = 99;
    char tight_ex[4][128]; int ntex = 0; char hcex[4][128]; int nhcex = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '&') continue;
        line[strcspn(line, "\r\n")] = 0;
        n = line[1] - 63; memset(out, 0, sizeof out);
        for (int k = 0; k < n * n; k++) { int byte = line[2 + k / 6] - 63; if (byte >> (5 - k % 6) & 1) out[k / n] |= 1u << (k % n); }
        if (!strong()) continue;
        tot++;
        memset(bef, 0, sizeof bef); H = 0;
        for (int s = 0; s < n; s++) { path[0] = s; dfs(1, 1u << s); }
        for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) balmat[x][y] = (x != y) && 3 * bef[x][y] >= H && 3 * bef[x][y] <= 2 * H;
        /* vertex balance */
        int is_tight = 0, struct_ok = 0;
        u64 minb_num = 1, minb_den = 1;
        for (int v = 0; v < n; v++) {
            u64 bv = 0; for (int w = 0; w < n; w++) if (w != v) { u64 c = bef[v][w], m = c < H - c ? c : H - c; if (m > bv) bv = m; }
            if ((u128)bv * minb_den < (u128)minb_num * H) { minb_num = bv; minb_den = H; }
        }
        if (3 * minb_num == minb_den) {
            is_tight = 1; tight++;
            /* structural test for each tight vertex */
            int all_struct = 1;
            for (int v = 0; v < n; v++) {
                u64 bv = 0; for (int w = 0; w < n; w++) if (w != v) { u64 c = bef[v][w], m = c < H - c ? c : H - c; if (m > bv) bv = m; }
                if (3 * bv != H) continue;
                int found = 0;
                for (int a = 0; a < n && !found; a++) for (int b = 0; b < n && !found; b++) {
                    if (a == v || b == v || a == b) continue;
                    if (!((out[v] >> a & 1) && (out[a] >> b & 1) && (out[b] >> v & 1))) continue;   /* v->a->b->v */
                    unsigned Tri = (1u << v) | (1u << a) | (1u << b), Rest = ((1u << n) - 1) ^ Tri;
                    if (!transitive_on(Rest)) continue;
                    int back = 0; for (int t = 0; t < n; t++) if (Tri >> t & 1) back += __builtin_popcount(out[t] & Rest);
                    if (back == 1) found = 1;
                }
                if (!found) all_struct = 0;
            }
            if (all_struct) { struct_ok = 1; tight_struct++; }
            if (ntex < 4) { snprintf(tight_ex[ntex], 128, "%s(struct=%d,H=%llu)", line, struct_ok, H); ntex++; }
        }
        (void)is_tight;
        /* balanced graph */
        unsigned seen = 1, frontier = 1;
        while (frontier) { int v = __builtin_ctz(frontier); frontier &= frontier - 1;
            for (int w = 0; w < n; w++) if (balmat[v][w] && !(seen >> w & 1)) { seen |= 1u << w; frontier |= 1u << w; } }
        if (seen == (1u << n) - 1) conn++;
        int hasdeg1 = 0; for (int v = 0; v < n; v++) { int d = 0; for (int w = 0; w < n; w++) d += balmat[v][w]; if (d == 1) hasdeg1 = 1; }
        deg1 += hasdeg1;
        /* Hamiltonian cycles */
        hamcyc = 0; hamcyc_nobal = 0; hc_minbal = 99; path[0] = 0; dfs_hc(1, 1u);
        hc_total += hamcyc; hc_nobal_total += hamcyc_nobal;
        if (hamcyc_nobal) { hc_tour_nobal++; if (nhcex < 4) { strncpy(hcex[nhcex], line, 127); hcex[nhcex][127] = 0; nhcex++; } }
        if (hc_minbal < global_minbal) { global_minbal = hc_minbal; minbal_count = 0; }
        if (hc_minbal == global_minbal) minbal_count++;
        if (doshort) for (int s = 0; s < n; s++) { path[0] = s; dfs_sc(s, 1, 1u << s); }
    }
    printf("n=%d: strong %lld | vertex-wise tight (min_v b(v) = 1/3): %lld, of which every tight vertex lies on a C_3 hung by one back arc "
           "over a transitive rest: %lld | balanced graph connected: %lld | some vertex of balanced degree 1: %lld\n",
           n, tot, tight, tight_struct, conn, deg1);
    for (int i = 0; i < ntex; i++) printf("   tight example: %s\n", tight_ex[i]);
    printf("   Hamiltonian cycles: %llu in total; without a balanced arc: %llu (in %lld tournaments); min #balanced arcs on a "
           "Hamiltonian cycle = %d, attained in %lld tournaments\n", hc_total, hc_nobal_total, hc_tour_nobal, global_minbal, minbal_count);
    for (int i = 0; i < nhcex; i++) printf("   tournament with a balanced-arc-free Hamiltonian cycle: %s\n", hcex[i]);
    if (doshort) { printf("   cycles of length k with no balanced arc (k: total/none):"); for (int k = 3; k < n; k++) printf(" %d: %llu/%llu", k, shortcyc[k], shortcyc_nobal[k]); printf("\n"); }
    return 0;
}
