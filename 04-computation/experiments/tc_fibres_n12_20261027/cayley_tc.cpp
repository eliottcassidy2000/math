// cayley_tc.cpp -- mac-mini-2026-10-27-potentials
// tc = H/|Aut| for every Cayley tournament on an abelian group Z_m1 x Z_m2 (odd order m = m1*m2 <= 19):
// connection sets S with S disjoint from -S and S u -S = G \ 0 (2^((m-1)/2) of them).  Reports the minimum tc and
// every Cayley tournament with tc <= 21 (evidence for the vertex-transitive lemma of THM-4812).
// Build: g++ -O3 -march=native -std=c++17 -I/opt/homebrew/opt/nauty/include -o cayley_tc cayley_tc.cpp /opt/homebrew/opt/nauty/lib/libnauty1.a
// Run:   ./cayley_tc m1 m2
#define MAXN WORDSIZE
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <map>
static int n;
static uint32_t outm[32];
static std::vector<uint64_t> dp;
static uint64_t hamCount() {
    int full = (1 << n) - 1;
    std::fill(dp.begin(), dp.end(), 0);
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
    for (int mask = 1; mask < full; mask++) {
        uint64_t *row = &dp[(size_t)mask * n];
        for (uint32_t vs = (uint32_t)mask; vs; vs &= vs - 1) {
            int v = __builtin_ctz(vs); uint64_t c = row[v]; if (!c) continue;
            for (uint32_t o = outm[v] & ~(uint32_t)mask; o; o &= o - 1) { int w = __builtin_ctz(o); dp[(size_t)(mask | (1 << w)) * n + w] += c; }
        }
    }
    uint64_t H = 0; for (int v = 0; v < n; v++) H += dp[(size_t)full * n + v];
    return H;
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
    int m1 = atoi(argv[1]), m2 = argc > 2 ? atoi(argv[2]) : 1; n = m1 * m2;
    auto id = [&](int a, int b) { return ((a % m1 + m1) % m1) * m2 + ((b % m2 + m2) % m2); };
    // pairs {g, -g}
    std::vector<std::pair<int,int>> reps; std::vector<char> seen(n, 0); seen[0] = 1;
    for (int a = 0; a < m1; a++) for (int b = 0; b < m2; b++) { int g = id(a, b); if (seen[g]) continue; int h = id(-a, -b); seen[g] = seen[h] = 1; reps.push_back({a, b}); }
    int K = (int)reps.size();
    dp.assign((size_t)(1 << n) * n, 0);
    uint64_t minTc = ~0ull; long long count = 0; std::map<uint64_t, long long> small;
    for (long long mask = 0; mask < (1ll << K); mask++) {
        std::vector<std::pair<int,int>> S;
        for (int k = 0; k < K; k++) { auto r = reps[k]; if (mask >> k & 1) S.push_back({-r.first, -r.second}); else S.push_back(r); }
        memset(outm, 0, sizeof outm);
        for (int a = 0; a < m1; a++) for (int b = 0; b < m2; b++) { int u = id(a, b); for (auto &s : S) outm[u] |= 1u << id(a + s.first, b + s.second); }
        uint64_t H = hamCount(); long long A = autOrder(); count++;
        uint64_t tc = H / (uint64_t)A;
        if (H % (uint64_t)A) { printf("ERROR: H not divisible by |Aut|\n"); return 1; }
        if (tc < minTc) minTc = tc;
        if (tc <= 21) small[tc]++;
    }
    printf("Cayley tournaments on Z_%d x Z_%d (order %d): %lld connection sets, min tc %llu, sets with tc <= 21:", m1, m2, n, count, (unsigned long long)minTc);
    for (auto &p : small) printf(" tc=%llu:%lld", (unsigned long long)p.first, p.second);
    printf("\n");
    return 0;
}
