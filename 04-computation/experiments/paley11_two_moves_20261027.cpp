// paley11_two_moves_20261027.cpp -- opus-2026-10-27-S63
// Exact depth-two census of native complement moves from the Paley clock P11 (x -> y iff y - x is a nonzero square mod 11).
// One move: T -> T^op xor E(P) for a Hamiltonian path P of T (every arc off P reversed).  Two moves compose to
// T -> T xor (E(P1) xor E(P2)) with P1 a Hamiltonian path of T and P2 one of the intermediate tournament, so the
// depth-two states are enumerated exactly by two nested Hamiltonian-path enumerations (H(P11) * H(T') path pairs).
// For every depth-two state U the program computes the score-deviation bound L(U) = (1/2) sum_i |s_(i) - i| <= f(U)
// (s_(i) the sorted out-degrees; every backward arc of an order moves two scores by one), then for the states with
// L <= 10 the exact minimum feedback-arc count f (subset DP), and for every state with f <= 10 an exact forest test
// (is some total order's backward-pair graph acyclic?  depth-first placement with union-find pruning).  It reports H(P11),
// the number of depth-two path pairs, the exact minimum f at depth one and two, the counts passing each filter, and
// whether any depth-two state is forest-framed (cycle debt 0), which decides session73's two-complement forest question.
// Build: g++ -O2 -std=c++17 -fopenmp -o p11 paley11_two_moves_20261027.cpp ; run: ./p11
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace std;
typedef unsigned int u32; typedef unsigned long long u64;
static const int n = 11, M = 55;
static int pidx[11][11]; static vector<pair<int,int>> pairs; static const u32 fullV = (1u << n) - 1;

struct Tour { u32 out[11]; };
static inline u64 code(const Tour& T) { u64 c = 0; for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; if (T.out[i] >> j & 1) c |= 1ULL << k; } return c; }
static void hpaths(const Tour& T, int v, u32 used, u64 pm, vector<u64>& res) {
    if (used == fullV) { res.push_back(pm); return; }
    u32 ext = T.out[v] & ~used; while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; hpaths(T, w, used | (1u << w), pm | (1ULL << pidx[v][w]), res); } }
static inline Tour applyToggle(const Tour& T, u64 toggle) {      // reverse the pairs in toggle
    Tour R = T; while (toggle) { int k = __builtin_ctzll(toggle); toggle &= toggle - 1; int i = pairs[k].first, j = pairs[k].second;
        if (R.out[i] >> j & 1) { R.out[i] &= ~(1u << j); R.out[j] |= 1u << i; } else { R.out[j] &= ~(1u << i); R.out[i] |= 1u << j; } } return R; }
static inline int scoreBound2(const Tour& T) {                      // 2 L(T) = sum_i |s_(i) - i|
    int s[11]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(T.out[v]); sort(s, s + n); int d = 0; for (int i = 0; i < n; i++) d += abs(s[i] - i); return d; }
static int fasDP(const Tour& T) { static thread_local uint8_t f[1 << n]; f[0] = 0;
    for (u32 S = 1; S <= fullV; S++) { int best = 255; u32 t = S; while (t) { int v = __builtin_ctz(t); t &= t - 1; u32 R = S ^ (1u << v); int c = f[R] + __builtin_popcount(T.out[v] & R); if (c < best) best = c; } f[S] = (uint8_t)best; }
    return f[fullV]; }
static bool forestDFS(const Tour& T, u32 placed, int* parent, int edges) {
    if (placed == fullV) return true;
    u32 cand = fullV & ~placed;
    while (cand) { int v = __builtin_ctz(cand); cand &= cand - 1; u32 back = T.out[v] & placed; int nb = __builtin_popcount(back);
        if (edges + nb > n - 1) continue;
        int save[11]; memcpy(save, parent, sizeof(save)); bool ok = true; u32 b = back;
        while (b && ok) { int u = __builtin_ctz(b); b &= b - 1; int a = v, c = u; while (parent[a] != a) a = parent[a]; while (parent[c] != c) c = parent[c]; if (a == c) ok = false; else parent[a] = c; }
        if (ok && forestDFS(T, placed | (1u << v), parent, edges + nb)) return true;
        memcpy(parent, save, sizeof(save)); }
    return false; }
static bool forestFramed(const Tour& T) { int parent[11]; iota(parent, parent + n, 0); return forestDFS(T, 0, parent, 0); }

int main() {
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pidx[i][j] = pidx[j][i] = (int)pairs.size(); pairs.push_back({i, j}); }
    Tour P; vector<int> sq(n, 0); for (int x = 1; x < n; x++) sq[(x * x) % n] = 1;
    for (int i = 0; i < n; i++) { P.out[i] = 0; for (int j = 0; j < n; j++) if (i != j && sq[((j - i) % n + n) % n]) P.out[i] |= 1u << j; }
    const u64 fullP = (1ULL << M) - 1;
    vector<u64> paths1; for (int v = 0; v < n; v++) hpaths(P, v, 1u << v, 0, paths1);
    printf("P11: f = %d, score bound L = %d, H(P11) = %zu Hamiltonian paths\n", fasDP(P), scoreBound2(P) / 2, paths1.size()); fflush(stdout);
    int minF1 = 999; vector<Tour> T1(paths1.size()); for (size_t a = 0; a < paths1.size(); a++) T1[a] = applyToggle(P, fullP ^ paths1[a]);
    #pragma omp parallel for reduction(min:minF1) schedule(dynamic, 64)
    for (long a = 0; a < (long)paths1.size(); a++) { int f = fasDP(T1[a]); if (f < minF1) minF1 = f; }
    printf("depth one: %zu states (one per path, with repetitions), exact minimum f = %d\n", paths1.size(), minF1); fflush(stdout);
    u64 pairsCount = 0, passL = 0, passF = 0; int minF2 = 999, minL2 = 999; vector<u64> candidates; u64 fHist[16] = {0}; u64 LHist[64] = {0}; long done = 0;
    #pragma omp parallel
    {   vector<u64> loc; u64 pc = 0, lp = 0, fp = 0, fh[16] = {0}, lh[64] = {0}; int mf = 999, ml = 999; vector<u64> paths2;
        #pragma omp for schedule(dynamic, 8)
        for (long a = 0; a < (long)paths1.size(); a++) { const Tour& T = T1[a]; paths2.clear(); for (int v = 0; v < n; v++) hpaths(T, v, 1u << v, 0, paths2);
            for (u64 pm : paths2) { pc++; Tour U = applyToggle(T, fullP ^ pm); int L2 = scoreBound2(U); if (L2 < ml) ml = L2; lh[min(L2, 63)]++; if (L2 > 20) continue; lp++;
                int f = fasDP(U); if (f < mf) mf = f; if (f < 16) fh[f]++; if (f <= 10) { fp++; loc.push_back(code(U)); } }
            #pragma omp critical
            { done++; if (done % 5000 == 0) { printf("  progress: %ld of %zu depth-one states expanded\n", done, paths1.size()); fflush(stdout); } } }
        #pragma omp critical
        { pairsCount += pc; passL += lp; passF += fp; if (mf < minF2) minF2 = mf; if (ml < minL2) minL2 = ml; for (int i = 0; i < 16; i++) fHist[i] += fh[i]; for (int i = 0; i < 64; i++) LHist[i] += lh[i]; candidates.insert(candidates.end(), loc.begin(), loc.end()); }
    }
    sort(candidates.begin(), candidates.end()); candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
    printf("depth two: %llu path pairs; minimum score bound L = %d; states with L <= 10: %llu; exact minimum f among them = %d (every other state has f >= 11)\n", pairsCount, (minL2 + 1) / 2, passL, minF2);
    printf("score-bound histogram 2L over all depth-two states:"); for (int i = 0; i < 64; i++) if (LHist[i]) printf(" %d:%llu", i, LHist[i]); printf("\n");
    printf("f histogram over the states with L <= 10 (with repetitions):"); for (int i = 0; i < 16; i++) if (fHist[i]) printf(" %d:%llu", i, fHist[i]); printf("\n");
    printf("states with f <= 10: %llu (%zu distinct)\n", passF, candidates.size()); fflush(stdout);
    int forests = 0;
    for (u64 c : candidates) { Tour U; for (int v = 0; v < n; v++) U.out[v] = 0; for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; if (c >> k & 1) U.out[i] |= 1u << j; else U.out[j] |= 1u << i; }
        bool fr = forestFramed(U); if (fr) { forests++; if (forests <= 5) printf("  forest-framed depth-two state: code %llu, f = %d\n", c, fasDP(U)); } }
    printf("forest-framed (cycle debt 0) states at depth two from P11: %d  => P11 %s two complement moves from a forest frame\n", forests, forests ? "IS within" : "is NOT within");
    return 0;
}
