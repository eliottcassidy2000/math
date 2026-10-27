// tc_invariant_scan.cpp -- mac-mini-2026-10-27-potentials
//
// Extends THM-4790 (presentation fibres tc = H/|Aut| never 7 or 21 through n = 11) by exhaustive scans of the
// tournaments invariant under a given permutation sigma of odd prime-power cycle structure.
//
// For n = 12: a strong S with tc(S) in {7,21} has |Aut(S)| >= f(12)/21 = 375/21, hence |Aut(S)| >= 19, and the
// type lemma of the accompanying note shows Aut(S) then contains an element of one of the types
//   (11), (7), (5)(5), (5)(3), (5)(3)(3), (3)(3)(3), (3)(3)(3)(3).
// This program enumerates all 2^(#pair orbits) sigma-invariant tournaments (Gray code over pair orbits), keeps the
// strong ones, computes H exactly (Held-Karp over 2^n subsets), and for every H of the form 7a or 21a with a odd,
// a >= amin and ord(sigma) | a, computes |Aut| with nauty and tests tc = H/|Aut|.
//
// Build: g++ -O3 -march=native -std=c++17 -I/opt/homebrew/opt/nauty/include -o tcscan tc_invariant_scan.cpp \
//            /opt/homebrew/opt/nauty/lib/libnauty1.a
// Run:   ./tcscan n amin cyc1 cyc2 ... [slice r m]      e.g. ./tcscan 12 19 3 3 3 ;  ./tcscan 12 19 3 3 3 slice 0 4
#define MAXN WORDSIZE
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <map>
#include <chrono>
using std::vector; using std::map; using std::pair; using std::min; using std::max;
namespace chrono = std::chrono;

static int n;
static uint16_t outm[16];

static bool strong() {
    // forward and backward reachability from 0
    uint16_t full = (uint16_t)((1u << n) - 1);
    uint16_t seen = 1, frontier = 1;
    while (frontier) { uint16_t nxt = 0; for (uint16_t f = frontier; f; f &= f - 1) nxt |= outm[__builtin_ctz(f)]; nxt &= ~seen; seen |= nxt; frontier = nxt; }
    if (seen != full) return false;
    uint16_t inm[16] = {0};
    for (int u = 0; u < n; u++) for (uint16_t o = outm[u]; o; o &= o - 1) inm[__builtin_ctz(o)] |= (uint16_t)(1u << u);
    seen = 1; frontier = 1;
    while (frontier) { uint16_t nxt = 0; for (uint16_t f = frontier; f; f &= f - 1) nxt |= inm[__builtin_ctz(f)]; nxt &= ~seen; seen |= nxt; frontier = nxt; }
    return seen == full;
}

static vector<uint64_t> dp;
static uint64_t hamCount() {
    int full = (1 << n) - 1;
    // dp[mask*n + v] = number of directed paths covering mask ending at v
    std::fill(dp.begin(), dp.end(), 0);
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
    for (int mask = 1; mask < full; mask++) {
        uint64_t *row = &dp[(size_t)mask * n];
        for (uint16_t vs = (uint16_t)mask; vs; vs &= vs - 1) {
            int v = __builtin_ctz(vs); uint64_t c = row[v]; if (!c) continue;
            for (uint16_t o = outm[v] & ~mask; o; o &= o - 1) { int w = __builtin_ctz(o); dp[(size_t)(mask | (1 << w)) * n + w] += c; }
        }
    }
    uint64_t H = 0; for (int v = 0; v < n; v++) H += dp[(size_t)full * n + v];
    return H;
}

static long long autOrder() {
    graph g[MAXN * 1]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    static DEFAULTOPTIONS_DIGRAPH(options); statsblk stats;
    EMPTYGRAPH(g, 1, n);
    for (int u = 0; u < n; u++) for (uint16_t o = outm[u]; o; o &= o - 1) ADDONEARC(g, u, __builtin_ctz(o), 1);
    options.getcanon = FALSE;
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, NULL);
    return (long long)(stats.grpsize1 * pow(10.0, stats.grpsize2) + 0.5);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: tcscan n amin cyc... [slice r m]\n"); return 1; }
    n = atoi(argv[1]); long long amin = atoll(argv[2]);
    vector<int> cyc; int ai = 3; long long sliceR = 0, sliceM = 1;
    for (; ai < argc; ai++) { if (!strcmp(argv[ai], "slice")) { sliceR = atoll(argv[ai + 1]); sliceM = atoll(argv[ai + 2]); break; } cyc.push_back(atoi(argv[ai])); }
    // sigma: consecutive cycles, remaining points fixed
    int sigma[16]; for (int i = 0; i < n; i++) sigma[i] = i;
    int pos = 0, ord = 1;
    for (int c : cyc) { for (int i = 0; i < c; i++) sigma[pos + i] = pos + (i + 1) % c; pos += c; long long g = ord; long long a = c; while (a) { long long t = g % a; g = a; a = t; } ord = ord / (int)g * c; }
    if (pos > n) { fprintf(stderr, "cycles exceed n\n"); return 1; }
    // pair orbits under <sigma>: each orbit stored as list of ordered pairs (u,v) obtained from a representative u<v
    vector<vector<pair<int,int>>> orbitsP; int seenP[16][16]; memset(seenP, 0, sizeof seenP);
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) {
        if (seenP[u][v]) continue;
        vector<pair<int,int>> orb; int a = u, b = v;
        while (!seenP[min(a,b)][max(a,b)]) { seenP[min(a,b)][max(a,b)] = 1; orb.push_back({a, b}); int na = sigma[a], nb = sigma[b]; a = na; b = nb; }
        // odd order: no pair is mapped to its own reversal
        orbitsP.push_back(orb);
    }
    int K = (int)orbitsP.size();
    printf("n=%d type", n); for (int c : cyc) printf(" (%d)", c); printf(" order %d, pair orbits %d, invariant tournaments 2^%d, amin %lld, slice %lld/%lld\n", ord, K, K, amin, sliceR, sliceM);
    fflush(stdout);
    dp.assign((size_t)(1 << n) * n, 0);
    // initial tournament: bit 0 for every orbit (orientation as stored), Gray code thereafter
    memset(outm, 0, sizeof outm);
    for (auto &orb : orbitsP) for (auto &p : orb) outm[p.first] |= (uint16_t)(1u << p.second);
    auto flipOrbit = [&](int k) { for (auto &p : orbitsP[k]) { int a = p.first, b = p.second; outm[a] ^= (uint16_t)(1u << b); outm[b] ^= (uint16_t)(1u << a); } };
    uint64_t total = 1ull << K; long long strongCount = 0, candidates = 0, hits = 0; uint64_t minH = ~0ull;
    map<long long, long long> autHist;
    auto t0 = chrono::steady_clock::now();
    for (uint64_t i = 0; i < total; i++) {
        if (i) { int k = __builtin_ctzll(i); flipOrbit(k); }
        if ((long long)(i % sliceM) != sliceR) continue;
        if (!strong()) continue;
        strongCount++;
        uint64_t H = hamCount();
        if (H < minH) minH = H;
        for (int tcv : {7, 21}) {
            if (H % tcv) continue; uint64_t a = H / tcv;
            if (!(a & 1) || (long long)a < amin || a % ord) continue;
            candidates++;
            long long A = autOrder(); autHist[A]++;
            if ((uint64_t)A == a) { hits++; printf("HIT tc=%d H=%llu |Aut|=%lld code", tcv, (unsigned long long)H, A); for (int u = 0; u < n; u++) printf(" %04x", outm[u]); printf("\n"); fflush(stdout); }
        }
    }
    double secs = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    printf("done: strong invariant %lld, min H %llu, |Aut| computations %lld, hits %lld, %.1f s\n", strongCount, (unsigned long long)minH, candidates, hits, secs);
    printf("|Aut| histogram of candidates:"); for (auto &p : autHist) printf(" %lld:%lld", p.first, p.second); printf("\n");
    return 0;
}
