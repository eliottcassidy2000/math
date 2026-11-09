/* lines_census.c -- the complement lines of the order-n metagraph (01-canon/definitions.md, THM-643/4729).
   Base path n -> n-1 -> ... -> 1; tiles (a,b) with a >= b+2; tiling bit 0: a -> b, bit 1: b -> a; line {t, t^c}.
   Every line has exactly one CLOSED end (tile (n,1) = 1, i.e. 1 -> n), so lines are enumerated by their closed
   tiling t.  For each line: the cut set of the open end (position i, between i+1 and i, is a cut of T_{t^c} iff every
   tile crossing i is 1 in t), its number of strong components, the blue flag (t grid-symmetric), the cyclic-triangle
   counts, and (with the H option) the Hamiltonian-path counts of both ends.
   usage: lines_census n [part nparts] [H]      -> lines_n{n}_p{part}.txt                                          */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n, m, tile_n1;
static int ta[64], tb[64], cross[64], sig[64];
static uint64_t hamcount(const uint32_t *adj) {
    static uint64_t dp[1 << 10][10];
    int full = (1 << n) - 1;
    for (int mask = 1; mask <= full; mask++) for (int v = 0; v < n; v++) dp[mask][v] = 0;
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int mask = 1; mask <= full; mask++) {
        for (int v = 0; v < n; v++) {
            uint64_t c = dp[mask][v]; if (!c) continue;
            uint32_t out = adj[v] & ~(uint32_t)mask;
            while (out) { int u = __builtin_ctz(out); out &= out - 1; dp[mask | (1 << u)][u] += c; }
        }
    }
    uint64_t H = 0; for (int v = 0; v < n; v++) H += dp[full][v];
    return H;
}
static int c3count(const uint32_t *adj) {
    int tot = n * (n - 1) * (n - 2) / 6;
    for (int v = 0; v < n; v++) { int d = __builtin_popcount(adj[v]); tot -= d * (d - 1) / 2; }
    return tot;
}
#define HMAX 65536
static uint64_t hist_c[HMAX], hist_o[HMAX];
#define PCAP (1u << 22)
static uint64_t *pkey; static uint32_t *pcnt;
static void pair_add(uint64_t key) {
    uint64_t h = (key * 0x9E3779B97F4A7C15ull) >> 42; h &= PCAP - 1;
    while (pcnt[h] && pkey[h] != key) h = (h + 1) & (PCAP - 1);
    pkey[h] = key; pcnt[h]++;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); int part = argc > 3 ? atoi(argv[2]) : 0, nparts = argc > 3 ? atoi(argv[3]) : 1;
    int withH = (argc > 4 && argv[4][0] == 'H') || (argc == 3 && argv[2][0] == 'H');
    if (argc == 3 && argv[2][0] == 'H') { part = 0; nparts = 1; }
    m = 0;
    for (int a = 3; a <= n; a++) for (int b = 1; b + 2 <= a; b++) { if (a == n && b == 1) continue; ta[m] = a; tb[m] = b; m++; }
    ta[m] = n; tb[m] = 1; tile_n1 = m; m++;   /* the closing tile is the top bit */
    for (int k = 0; k < m; k++) {
        cross[k] = 0; for (int i = tb[k]; i < ta[k]; i++) cross[k] |= 1 << (i - 1);
        int sa = n + 1 - tb[k], sb = n + 1 - ta[k]; sig[k] = -1;
        for (int l = 0; l < m; l++) if (ta[l] == sa && tb[l] == sb) sig[k] = l;
        if (sig[k] < 0) { fprintf(stderr, "sigma failure\n"); return 1; }
    }
    pkey = calloc(PCAP, 8); pcnt = calloc(PCAP, 4);
    uint64_t total = 1ull << (m - 1);
    uint64_t lo = total / nparts * part, hi = (part == nparts - 1) ? total : total / nparts * (part + 1);
    static uint64_t cutall[1 << 9], cutblue[1 << 9], kall[16], kblue[16], c3sum[200], c3joint[90][90];
    uint64_t nblue = 0, nlines = 0, strong_strong = 0, blue_ss = 0;
    uint32_t adj[16], adjc[16];
    for (uint64_t u = lo; u < hi; u++) {
        uint64_t t = u | (1ull << tile_n1);
        for (int v = 0; v < n; v++) { adj[v] = 0; adjc[v] = 0; }
        for (int a = 2; a <= n; a++) { adj[a - 1] |= 1u << (a - 2); adjc[a - 1] |= 1u << (a - 2); }
        int cuts = (1 << (n - 1)) - 1;   /* cuts of the open end: positions where every crossing tile is 1 in t */
        int blue = 1;
        for (int k = 0; k < m; k++) {
            int bit = (t >> k) & 1;
            if (bit) { adj[tb[k] - 1] |= 1u << (ta[k] - 1); adjc[ta[k] - 1] |= 1u << (tb[k] - 1); }
            else { adj[ta[k] - 1] |= 1u << (tb[k] - 1); adjc[tb[k] - 1] |= 1u << (ta[k] - 1); cuts &= ~cross[k]; }
            if (((t >> sig[k]) & 1) != (unsigned)bit) blue = 0;
        }
        int kc = __builtin_popcount(cuts) + 1;
        nlines++; nblue += blue; cutall[cuts]++; kall[kc]++; if (blue) { cutblue[cuts]++; kblue[kc]++; }
        if (kc == 1) { strong_strong++; if (blue) blue_ss++; }
        int c1 = c3count(adj), c2 = c3count(adjc); c3sum[c1 + c2]++; if (c1 < 90 && c2 < 90) c3joint[c1][c2]++;
        if (withH) {
            uint64_t Hc = hamcount(adj), Ho = hamcount(adjc);
            if (Hc < HMAX) hist_c[Hc]++; if (Ho < HMAX) hist_o[Ho]++;
            pair_add((Hc << 20) | Ho | ((uint64_t)blue << 40) | ((uint64_t)(kc == 1) << 41));
        }
    }
    char name[256]; snprintf(name, sizeof name, "lines_n%d_p%d.txt", n, part);
    FILE *o = fopen(name, "w");
    fprintf(o, "# n=%d m=%d part %d/%d: lines %llu blue %llu both-ends-strong %llu (blue %llu)\n", n, m, part, nparts, (unsigned long long)nlines, (unsigned long long)nblue, (unsigned long long)strong_strong, (unsigned long long)blue_ss);
    fprintf(o, "#k\n"); for (int k = 1; k <= n; k++) fprintf(o, "%d %llu %llu\n", k, (unsigned long long)kall[k], (unsigned long long)kblue[k]);
    fprintf(o, "#cuts\n"); for (int c = 0; c < (1 << (n - 1)); c++) if (cutall[c]) fprintf(o, "%d %llu %llu\n", c, (unsigned long long)cutall[c], (unsigned long long)cutblue[c]);
    fprintf(o, "#c3sum\n"); for (int s = 0; s < 200; s++) if (c3sum[s]) fprintf(o, "%d %llu\n", s, (unsigned long long)c3sum[s]);
    fprintf(o, "#c3joint\n"); for (int a = 0; a < 90; a++) for (int b = 0; b < 90; b++) if (c3joint[a][b]) fprintf(o, "%d %d %llu\n", a, b, (unsigned long long)c3joint[a][b]);
    if (withH) {
        fprintf(o, "#Hclosed\n"); for (int h = 0; h < HMAX; h++) if (hist_c[h]) fprintf(o, "%d %llu\n", h, (unsigned long long)hist_c[h]);
        fprintf(o, "#Hopen\n"); for (int h = 0; h < HMAX; h++) if (hist_o[h]) fprintf(o, "%d %llu\n", h, (unsigned long long)hist_o[h]);
        fprintf(o, "#pairs Hc Ho blue openstrong count\n");
        for (uint64_t i = 0; i < PCAP; i++) if (pcnt[i]) fprintf(o, "%llu %llu %d %d %u\n", (unsigned long long)(pkey[i] >> 20) & 0xFFFFF, (unsigned long long)(pkey[i] & 0xFFFFF), (int)((pkey[i] >> 40) & 1), (int)((pkey[i] >> 41) & 1), pcnt[i]);
    }
    fclose(o);
    fprintf(stderr, "n=%d part %d done: lines %llu blue %llu\n", n, part, (unsigned long long)nlines, (unsigned long long)nblue);
    return 0;
}
