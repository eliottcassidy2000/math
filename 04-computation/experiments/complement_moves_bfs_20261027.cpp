// complement_moves_bfs_20261027.cpp -- opus-2026-10-27-S63
// The native complement move on labelled tournaments: for a Hamiltonian path P of T, T -> T^op xor E(P), i.e. every arc
// not on P is reversed and P is kept (the "legal path complement" of the weighted-renewal / clock-even / session73 grammar;
// re-marking the path is free, so the move graph has an edge T -- T' for every Hamiltonian path P of T, and the same P
// applied twice is the identity).  This program runs BFS in that graph on all 2^(n(n-1)/2) labelled tournaments, n <= 7:
// from the transitive tournament (distance = the number of complement moves needed = the exception budget of the native
// grammar relative to the transitive baseline), and from the Paley clock P_n when n is a prime 3 mod 4.  For every
// tournament it also computes f(T) = the minimum feedback-arc count (subset DP), and compares the distance with the
// exact lower bound ceil(f(T)/(n-1)) (one move toggles n-1 ... exactly the M-(n-1) pairs off the path; relative to the
// reversed order of a best order of T, the backward set changes by the n-1 path pairs, so f changes by at most n-1 per move).
// With the flag "diam" (n <= 6) it also computes the exact diameter of the move graph by BFS from every tournament.
// Build: g++ -O2 -std=c++17 -o cmb complement_moves_bfs_20261027.cpp ; run: ./cmb n [diam]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;
typedef unsigned int u32; typedef unsigned long long u64;

static int n, M; static int pidx[8][8]; static vector<pair<int,int>> pairs; static u32 fullV, fullP;

static inline void decode(u32 code, u32* out) { for (int v = 0; v < n; v++) out[v] = 0;
    for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; if (code >> k & 1) out[i] |= 1u << j; else out[j] |= 1u << i; } }
static void hpaths(const u32* out, int v, u32 used, u32 pm, vector<u32>& res) {
    if (used == fullV) { res.push_back(pm); return; }
    u32 ext = out[v] & ~used; while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; hpaths(out, w, used | (1u << w), pm | (1u << pidx[v][w]), res); } }
static inline int fas(const u32* out) { static uint8_t f[128]; f[0] = 0;
    for (u32 S = 1; S <= fullV; S++) { int best = 255; u32 t = S; while (t) { int v = __builtin_ctz(t); t &= t - 1; u32 R = S ^ (1u << v); int c = f[R] + __builtin_popcount(out[v] & R); if (c < best) best = c; } f[S] = (uint8_t)best; }
    return f[fullV]; }
static u32 paleyCode() { u32 code = 0; vector<int> sq(n, 0); for (int x = 1; x < n; x++) sq[(x * x) % n] = 1;
    for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; if (sq[(j - i) % n]) code |= 1u << k; } return code; }

static int bfs(u32 src, vector<uint8_t>& dist, vector<u32>& queue) {
    fill(dist.begin(), dist.end(), 255); dist[src] = 0; queue.clear(); queue.push_back(src); size_t head = 0; int ecc = 0; u32 out[8]; vector<u32> res;
    while (head < queue.size()) { u32 c = queue[head++]; int d = dist[c]; if (d > ecc) ecc = d; decode(c, out); res.clear();
        for (int v = 0; v < n; v++) hpaths(out, v, 1u << v, 0, res);
        for (u32 pm : res) { u32 nb = c ^ (fullP ^ pm); if (dist[nb] == 255) { dist[nb] = (uint8_t)(d + 1); queue.push_back(nb); } } }
    return ecc; }

int main(int argc, char** argv) {
    n = argc > 1 ? atoi(argv[1]) : 5; bool diam = argc > 2 && string(argv[2]) == "diam"; if (n < 3 || n > 7) { fprintf(stderr, "3 <= n <= 7\n"); return 1; }
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pidx[i][j] = pidx[j][i] = (int)pairs.size(); pairs.push_back({i, j}); }
    M = (int)pairs.size(); u32 N = 1u << M; fullV = (1u << n) - 1; fullP = N - 1;
    vector<uint8_t> dist(N); vector<u32> queue; queue.reserve(N);
    // f(T) for every tournament
    vector<uint8_t> F(N); { u32 out[8]; for (u32 c = 0; c < N; c++) { decode(c, out); F[c] = (uint8_t)fas(out); } }
    u32 trans = fullP;                                  // i -> j for all i < j
    int ecc = bfs(trans, dist, queue); u32 reached = (u32)queue.size();
    printf("n = %d: %u labelled tournaments; BFS from the transitive tournament: reached %u (%s), eccentricity %d\n", n, N, reached, reached == N ? "connected" : "NOT connected", ecc);
    vector<u64> level(ecc + 2, 0); for (u32 c = 0; c < N; c++) if (dist[c] != 255) level[dist[c]]++;
    printf("level counts:"); for (int d = 0; d <= ecc; d++) printf(" %d:%llu", d, level[d]); printf("\n");
    { u64 par[2] = {0, 0}, parAll[2] = {0, 0}; for (u32 c = 0; c < N; c++) { int p = __builtin_popcount(c) & 1; parAll[p]++; if (dist[c] != 255) par[p]++; }
      printf("forward-arc parity of reached states: even %llu, odd %llu (all tournaments: even %llu, odd %llu); one move toggles (n-1)(n-2)/2 = %d pairs (%s)\n", par[0], par[1], parAll[0], parAll[1], (n - 1) * (n - 2) / 2, ((n - 1) * (n - 2) / 2) % 2 ? "odd: no parity obstruction" : "even: parity is invariant"); }
    // comparison with the lower bound ceil(f/(n-1))
    vector<u64> slack(16, 0); u64 tight = 0; int maxF = 0; vector<int> minFlev(ecc + 1, 999), maxFlev(ecc + 1, -1);
    for (u32 c = 0; c < N; c++) { if (dist[c] == 255) continue; int f = F[c], d = dist[c]; int lb = (f + n - 2) / (n - 1); int s = d - lb; if (s < 0) { printf("LOWER BOUND VIOLATED at code %u: f = %d, dist = %d\n", c, f, d); return 2; }
        slack[min(s, 15)]++; if (s == 0) tight++; if (f > maxF) maxF = f; minFlev[d] = min(minFlev[d], f); maxFlev[d] = max(maxFlev[d], f); }
    printf("max f(T) = %d; dist - ceil(f/(n-1)) histogram:", maxF); for (int s = 0; s < 16; s++) if (slack[s]) printf(" %d:%llu", s, slack[s]); printf("  (tight: %llu)\n", tight);
    printf("per level [min f, max f]:"); for (int d = 0; d <= ecc; d++) printf(" %d:[%d,%d]", d, minFlev[d], maxFlev[d]); printf("\n");
    if (n % 4 == 3 && (n == 3 || n == 7)) { u32 pc = paleyCode(); printf("Paley P_%d: code %u, f = %d, lower bound %d, forest one-step ceiling 2n-4 = %d, distance from transitive = %d\n", n, pc, F[pc], (F[pc] + n - 2) / (n - 1), 2 * n - 4, dist[pc]);
        vector<uint8_t> d2(N); int e2 = bfs(pc, d2, queue); vector<u64> lv(e2 + 1, 0); for (u32 c = 0; c < N; c++) if (d2[c] != 255) lv[d2[c]]++;
        printf("BFS from P_%d: reached %zu, eccentricity %d, level counts:", n, queue.size(), e2); for (int d = 0; d <= e2; d++) printf(" %d:%llu", d, lv[d]); printf("\n"); }
    if (n == 5) { u32 pc = paleyCode(); printf("regular tournament on 5 vertices (Z_5 with {1,4}): f = %d, distance from transitive = %d\n", F[pc], dist[pc]); }
    if (diam && n <= 6) { int D = 0; u64 pairsAt[32] = {0}; vector<uint8_t> dd(N);
        for (u32 s = 0; s < N; s++) { int e = bfs(s, dd, queue); if (e > D) D = e; for (u32 c = 0; c < N; c++) pairsAt[dd[c] == 255 ? 31 : dd[c]]++; }
        printf("diameter of the move graph: %d; ordered-pair distance counts:", D); for (int d = 0; d < 32; d++) if (pairsAt[d]) printf(" %d:%llu", d == 31 ? -1 : d, pairsAt[d]); printf("\n"); }
    return 0;
}
