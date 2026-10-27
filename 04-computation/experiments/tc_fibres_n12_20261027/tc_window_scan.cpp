// tc_window_scan.cpp -- mac-mini-2026-10-27-potentials
// Exhaustive scan of one window job written by window_configs.py: every G-invariant tournament with the prescribed
// orbit structure (within-orbit arcs fixed, one orientation per pair of mutually inverse between-orbit orbitals).
// For each strong one, H is computed exactly -- by the block-run (substitution) DP when every between-orbit orbital is a
// full rectangle (all orbits are modules), otherwise by Held-Karp -- and whenever H = 21a or 7a with |G| dividing a,
// |Aut| is computed with nauty and tc = H/|Aut| is tested.  With "check" the two H methods are compared on every
// module tournament with n <= 16.
// Build: g++ -O3 -march=native -std=c++17 -I/opt/homebrew/opt/nauty/include -o tc_window_scan tc_window_scan.cpp /opt/homebrew/opt/nauty/lib/libnauty1.a
// Run:   ./tc_window_scan job [check]
#define MAXN WORDSIZE
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
using std::vector; using std::string;
static int n;
static uint32_t outm[32];
static long long autOrder() {
    graph g[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    static DEFAULTOPTIONS_DIGRAPH(options); statsblk stats;
    EMPTYGRAPH(g, 1, n);
    for (int u = 0; u < n; u++) for (uint32_t o = outm[u]; o; o &= o - 1) ADDONEARC(g, u, __builtin_ctz(o), 1);
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, NULL);
    return (long long)(stats.grpsize1 * pow(10.0, stats.grpsize2) + 0.5);
}
static bool strong() {
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    uint32_t seen = 1, fr = 1;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= outm[__builtin_ctz(f)]; nx &= ~seen; seen |= nx; fr = nx; }
    if (seen != full) return false;
    uint32_t inm[32] = {0};
    for (int u = 0; u < n; u++) for (uint32_t o = outm[u]; o; o &= o - 1) inm[__builtin_ctz(o)] |= 1u << u;
    seen = 1; fr = 1;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= inm[__builtin_ctz(f)]; nx &= ~seen; seen |= nx; fr = nx; }
    return seen == full;
}
static vector<uint64_t> hk;
static uint64_t heldKarp() {
    int full = (1 << n) - 1;
    if (hk.size() != (size_t)(1 << n) * n) hk.assign((size_t)(1 << n) * n, 0); else std::fill(hk.begin(), hk.end(), 0);
    for (int v = 0; v < n; v++) hk[(size_t)(1 << v) * n + v] = 1;
    for (int mask = 1; mask < full; mask++) {
        uint64_t *row = &hk[(size_t)mask * n];
        for (uint32_t vs = (uint32_t)mask; vs; vs &= vs - 1) {
            int v = __builtin_ctz(vs); uint64_t c = row[v]; if (!c) continue;
            for (uint32_t o = outm[v] & ~(uint32_t)mask; o; o &= o - 1) { int w = __builtin_ctz(o); hk[(size_t)(mask | (1 << w)) * n + w] += c; }
        }
    }
    uint64_t H = 0; for (int v = 0; v < n; v++) H += hk[(size_t)full * n + v];
    return H;
}
// ---- block-run DP for substitutions ----
static int nb;                                   // number of blocks (orbits)
static vector<vector<int>> blocks;
static vector<vector<uint64_t>> pcov;            // pcov[i][k] ordered k-path covers of block i
static vector<int> blockOf;
static void pathCovers() {
    pcov.assign(nb, {});
    for (int i = 0; i < nb; i++) {
        int m = (int)blocks[i].size(); int full = (1 << m) - 1;
        // W[mask] = directed paths with vertex set exactly mask
        vector<vector<uint64_t>> e(1 << m, vector<uint64_t>(m, 0));
        for (int v = 0; v < m; v++) e[1 << v][v] = 1;
        for (int mask = 1; mask <= full; mask++) for (int v = 0; v < m; v++) if (e[mask][v])
            for (int w = 0; w < m; w++) if (!(mask >> w & 1) && (outm[blocks[i][v]] >> blocks[i][w] & 1)) e[mask | 1 << w][w] += e[mask][v];
        vector<uint64_t> W(1 << m, 0); for (int mask = 1; mask <= full; mask++) for (int v = 0; v < m; v++) W[mask] += e[mask][v];
        vector<vector<uint64_t>> F(m + 1, vector<uint64_t>(1 << m, 0)); F[0][0] = 1;
        for (int k = 1; k <= m; k++) for (int mask = 1; mask <= full; mask++) {
            uint64_t s = 0; for (int sub = mask; sub; sub = (sub - 1) & mask) if (F[k - 1][mask ^ sub]) s += W[sub] * F[k - 1][mask ^ sub];
            F[k][mask] = s;
        }
        pcov[i].assign(m + 1, 0); for (int k = 1; k <= m; k++) pcov[i][k] = F[k][full];
    }
}
static uint64_t moduleDP() {
    // quotient arcs
    vector<vector<int>> q(nb, vector<int>(nb, 0));
    for (int i = 0; i < nb; i++) for (int j = 0; j < nb; j++) if (i != j) q[i][j] = (outm[blocks[i][0]] >> blocks[j][0]) & 1;
    // state: counts vector (mixed radix) x last block
    vector<int> rad(nb), mul(nb); long long S = 1;
    for (int i = 0; i < nb; i++) { rad[i] = (int)blocks[i].size() + 1; mul[i] = (int)S; S *= rad[i]; }
    vector<uint64_t> cur((size_t)S * nb, 0);
    for (int i = 0; i < nb; i++) cur[(size_t)mul[i] * nb + i] = 1;
    // process states in increasing total count order: iterate over code; counts only increase so code increases
    uint64_t H = 0;
    vector<int> cnt(nb);
    for (long long code = 0; code < S; code++) {
        long long c = code; int tot = 0; bool ok = true;
        for (int i = 0; i < nb; i++) { cnt[i] = (int)(c % rad[i]); c /= rad[i]; tot += cnt[i]; }
        for (int last = 0; last < nb; last++) {
            uint64_t v = cur[(size_t)code * nb + last]; if (!v) continue;
            bool all = true; uint64_t prod = 1;
            for (int i = 0; i < nb; i++) { if (cnt[i] == 0) { all = false; break; } prod *= pcov[i][cnt[i]]; }
            if (all) H += v * prod;
            for (int j = 0; j < nb; j++) if (j != last && q[last][j] && cnt[j] + 1 < rad[j]) cur[(size_t)(code + mul[j]) * nb + j] += v;
        }
        (void)ok; (void)tot;
    }
    return H;
}
int main(int argc, char **argv) {
    std::ifstream in(argv[1]); bool check = argc > 2 && !strcmp(argv[2], "check");
    string line; std::getline(in, line); int K; long long Gsize; { std::istringstream ss(line); ss >> n >> K >> Gsize; }
    std::getline(in, line); { std::istringstream ss(line); string tok; while (ss >> tok) { vector<int> b; std::istringstream bs(tok); string x; while (std::getline(bs, x, '.')) b.push_back(atoi(x.c_str())); blocks.push_back(b); } }
    nb = (int)blocks.size(); blockOf.assign(n, -1); for (int i = 0; i < nb; i++) for (int v : blocks[i]) blockOf[v] = i;
    auto parseArcs = [](const string &l) { vector<std::pair<int,int>> a; std::istringstream ss(l); string tok; while (ss >> tok) { size_t c = tok.find(','); a.push_back({atoi(tok.substr(0, c).c_str()), atoi(tok.substr(c + 1).c_str())}); } return a; };
    std::getline(in, line); auto fixedArcs = parseArcs(line);
    vector<vector<std::pair<int,int>>> orb(K);
    bool rect = true;
    for (int k = 0; k < K; k++) { std::getline(in, line); orb[k] = parseArcs(line); int a = blockOf[orb[k][0].first], b = blockOf[orb[k][0].second]; if ((int)orb[k].size() != (int)(blocks[a].size() * blocks[b].size())) rect = false; }
    memset(outm, 0, sizeof outm);
    for (auto &e : fixedArcs) outm[e.first] |= 1u << e.second;
    for (int k = 0; k < K; k++) for (auto &e : orb[k]) outm[e.first] |= 1u << e.second;
    // sanity: tournament
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) { int a = (outm[u] >> v) & 1, b = (outm[v] >> u) & 1; if (a + b != 1) { fprintf(stderr, "not a tournament at %d %d\n", u, v); return 1; } }
    if (rect) pathCovers();
    if (!rect && n > 18) { fprintf(stderr, "non-module job at n=%d too large for Held-Karp\n", n); return 1; }
    long long strongCount = 0, cand = 0, hits = 0, checks = 0, checkFail = 0; uint64_t minH = ~0ull;
    for (long long i = 0; i < (1ll << K); i++) {
        if (i) { int k = __builtin_ctzll(i); for (auto &e : orb[k]) { outm[e.first] ^= 1u << e.second; outm[e.second] ^= 1u << e.first; } }
        if (!strong()) continue;
        strongCount++;
        uint64_t H;
        if (rect) { H = moduleDP(); if (check && n <= 16) { checks++; if (heldKarp() != H) checkFail++; } }
        else H = heldKarp();
        if (H < minH) minH = H;
        for (int t : {7, 21}) {
            if (H % t) continue; uint64_t a = H / t; if (a % (uint64_t)Gsize) continue;
            cand++; long long A = autOrder();
            if ((uint64_t)A == a) { hits++; printf("HIT tc=%d H=%llu |Aut|=%lld\n", t, (unsigned long long)H, A); }
        }
    }
    printf("%s: n=%d blocks=%d K=%d |G|=%lld modules=%s strong=%lld minH=%llu autcalls=%lld hits=%lld", argv[1], n, nb, K, Gsize, rect ? "yes" : "no", strongCount, (unsigned long long)minH, cand, hits);
    if (check) printf(" moduleDP-vs-HeldKarp checks=%lld failures=%lld", checks, checkFail);
    printf("\n");
    return 0;
}
