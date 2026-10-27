// tc_strong_fibres.cpp -- mac-mini-2026-10-27-potentials (strong-only variant of tc_spectrum_sweep.cpp)
// Strong presentation fibres: for each tc value <= VMAX, the number of STRONG classes of order n with that tc.
// Input: `gentourng -cq n` (strong classes only).
// Reads `gentourng -q n` (upper triangle row by row) for one order and prints every attained tc value <= VMAX with
// the number of classes attaining it.  Merge the per-order files to get the spectrum of orders <= 10.
// Build: g++ -O3 -march=native -std=c++17 -I/opt/homebrew/opt/nauty/include -o tc_spectrum_sweep tc_spectrum_sweep.cpp /opt/homebrew/opt/nauty/lib/libnauty1.a
// Run:   gentourng -q n | ./tc_spectrum_sweep n VMAX > tc_spectrum_n.txt
#define MAXN WORDSIZE
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <map>
static int n;
static uint32_t outm[32];
static uint64_t dp[1 << 10][10];
static uint64_t ham() {
    int full = (1 << n) - 1;
    for (int m = 0; m <= full; m++) for (int v = 0; v < n; v++) dp[m][v] = 0;
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int m = 1; m < full; m++) for (uint32_t vs = (uint32_t)m; vs; vs &= vs - 1) {
        int v = __builtin_ctz(vs); uint64_t c = dp[m][v]; if (!c) continue;
        for (uint32_t o = outm[v] & ~(uint32_t)m; o; o &= o - 1) { int w = __builtin_ctz(o); dp[m | (1 << w)][w] += c; }
    }
    uint64_t H = 0; for (int v = 0; v < n; v++) H += dp[full][v]; return H;
}
static long long autOrder() {
    graph g[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    static DEFAULTOPTIONS_DIGRAPH(options); statsblk stats;
    EMPTYGRAPH(g, 1, n);
    for (int u = 0; u < n; u++) for (uint32_t o = outm[u]; o; o &= o - 1) ADDONEARC(g, u, __builtin_ctz(o), 1);
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, NULL);
    return (long long)(stats.grpsize1 * pow(10.0, stats.grpsize2) + 0.5);
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); long long vmax = atoll(argv[2]);
    char line[256]; std::map<long long, long long> spec; long long classes = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (n > 1 && line[0] != '0' && line[0] != '1') continue;
        memset(outm, 0, sizeof outm); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') outm[i] |= 1u << j; else outm[j] |= 1u << i; k++; }
        classes++;
        uint64_t H = ham(); long long A = autOrder();
        if (H % (uint64_t)A) { printf("ERROR divisibility\n"); return 1; }
        long long tc = (long long)(H / (uint64_t)A);
        if (tc <= vmax) spec[tc]++;
    }
    printf("# order %d: %lld strong classes; attained tc values <= %lld (value count)\n", n, classes, vmax);
    for (auto &p : spec) printf("%lld %lld\n", p.first, p.second);
    return 0;
}
