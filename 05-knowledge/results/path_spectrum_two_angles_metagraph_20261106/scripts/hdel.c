/* hdel.c (klein-2026-11-06) -- the vertical gradient of the Hamiltonian-path count: h(T) against h(T - v).
   usage: hdel n < gentourng-output
   For every tournament T read and every vertex v, computes h(T - v) (subset DP on the deletion). Prints per class
     "D h c3 strong  hmin hmax"   (the least and greatest h(T - v) over v)
   and aggregate lines: the distribution of the integer part of the ratios h / hmax and h / hmin, the extreme ratios,
   and for the values new at order n (not attained at order n - 1: read from the census file census_<n-1>.txt given as
   argv[2]) the distribution of h / hmax.  By the insertion lemma h(T - v) <= h(T) <= n * h(T - v). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef unsigned int u32; typedef u32 cell;
#define HMAX (1u << 18)
static int n; static u64 adj[32], radj[32]; static cell *f;
static u64 count_paths_sub(u64 V) {   /* Hamiltonian paths of the subtournament induced on vertex set V (bitmask) */
    int m = __builtin_popcountll(V); int map[32], k = 0; u64 x = V; while (x) { map[k++] = __builtin_ctzll(x); x &= x - 1; }
    u64 ra[32]; for (int i = 0; i < m; i++) { ra[i] = 0; for (int j = 0; j < m; j++) if (radj[map[i]] >> map[j] & 1) ra[i] |= 1ULL << j; }
    u64 N = 1ULL << m;
    for (int v = 0; v < m; v++) f[(1ULL << v) * m + v] = 1;
    for (u64 S = 1; S < N; S++) {
        if (!(S & (S - 1))) continue;
        u64 vs = S; cell *dst = f + S * m;
        while (vs) { int v = __builtin_ctzll(vs); vs &= vs - 1; u64 R = S & ~(1ULL << v); u64 in = ra[v] & R; u64 c = 0; const cell *row = f + R * m;
            while (in) { int u = __builtin_ctzll(in); in &= in - 1; c += row[u]; } dst[v] = c; }
    }
    u64 h = 0; for (int v = 0; v < m; v++) h += f[(N - 1) * m + v]; return h;
}
static int is_strong(void) {
    u64 all = (1ULL << n) - 1, seen = 1, fr = 1;
    while (fr) { u64 nf = 0, x = fr; while (x) { int v = __builtin_ctzll(x); x &= x - 1; nf |= adj[v]; } nf &= ~seen; seen |= nf; fr = nf; }
    if (seen != all) return 0; seen = 1; fr = 1;
    while (fr) { u64 nf = 0, x = fr; while (x) { int v = __builtin_ctzll(x); x &= x - 1; nf |= radj[v]; } nf &= ~seen; seen |= nf; fr = nf; }
    return seen == all;
}
int main(int argc, char **argv) {
    if (argc < 2) return 1; n = atoi(argv[1]); u64 N = 1ULL << n; f = malloc(sizeof(cell) * N * n);
    unsigned char *prev = calloc(HMAX, 1);
    if (argc > 2) { FILE *fp = fopen(argv[2], "r"); char ln[1024]; while (fp && fgets(ln, sizeof ln, fp)) { if (ln[0] == 'T') { u64 h = strtoull(ln + 2, 0, 10); if (h < HMAX) prev[h] = 1; } } if (fp) fclose(fp); }
    char line[1024]; int L0 = n * (n - 1) / 2; long long cls = 0;
    long long hist_max[64] = {0}, hist_min[64] = {0}, hist_new[64] = {0}; double rmax_max = 0, rmin_max = 0, rmax_min = 1e9; long long newcls = 0;
    double rnew_min = 1e9; u64 rnew_min_h = 0;
    while (fgets(line, sizeof line, stdin)) {
        int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0; if (L != L0) continue;
        memset(adj, 0, sizeof adj); memset(radj, 0, sizeof radj); int idx = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') { adj[i] |= 1ULL << j; radj[j] |= 1ULL << i; } else { adj[j] |= 1ULL << i; radj[i] |= 1ULL << j; } }
        u64 h = count_paths_sub(N - 1); int c3 = n * (n - 1) * (n - 2) / 6; for (int v = 0; v < n; v++) { int s = __builtin_popcountll(adj[v]); c3 -= s * (s - 1) / 2; }
        int st = is_strong(); u64 hmin = ~0ULL, hmax = 0;
        for (int v = 0; v < n; v++) { u64 hv = count_paths_sub((N - 1) & ~(1ULL << v)); if (hv < hmin) hmin = hv; if (hv > hmax) hmax = hv; }
        cls++;
        double rmax = (double)h / hmax, rmin = (double)h / hmin;   /* rmax <= rmin */
        { int ia = (int)(rmax * 4); if (ia > 63) ia = 63; hist_max[ia]++; int ib = (int)(rmin * 4); if (ib > 63) ib = 63; hist_min[ib]++; }
        if (rmax > rmax_max) rmax_max = rmax; if (rmin > rmin_max) rmin_max = rmin; if (rmax < rmax_min) rmax_min = rmax;
        if (!prev[h]) { newcls++; { int ic = (int)(rmax * 4); if (ic > 63) ic = 63; hist_new[ic]++; } if (rmax < rnew_min) { rnew_min = rmax; rnew_min_h = h; } }
        printf("D %llu %d %d %llu %llu\n", h, c3, st, hmin, hmax);
    }
    printf("# order %d: %lld classes; h / max_v h(T-v): min %.4f max %.4f; h / min_v h(T-v): max %.4f\n", n, cls, rmax_min, rmax_max, rmin_max);
    printf("# histogram of h / max_v h(T-v) in quarter-steps [x/4, (x+1)/4):"); for (int i = 0; i < 64; i++) if (hist_max[i]) printf(" %d:%lld", i, hist_max[i]); printf("\n");
    printf("# histogram of h / min_v h(T-v) (bin 63 = ratio >= 15.75):"); for (int i = 0; i < 64; i++) if (hist_min[i]) printf(" %d:%lld", i, hist_min[i]); printf("\n");
    if (argc > 2) { printf("# classes with a value new at order %d: %lld; histogram of h / max_v h(T-v):", n, newcls); for (int i = 0; i < 64; i++) if (hist_new[i]) printf(" %d:%lld", i, hist_new[i]); printf("; least ratio %.4f at h = %llu\n", rnew_min, rnew_min_h); }
    return 0;
}
