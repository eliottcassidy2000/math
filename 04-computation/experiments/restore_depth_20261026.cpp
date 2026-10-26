// restore_depth_20261026.cpp -- opus-2026-10-26-S62
// Restore depth on the marked carrier X_n: for each state x, the least k such that k of its own Hamiltonian paths have
// free pairs covering all M = (n-1)(n-2)/2 free pairs (then the point indicator of x is the k-fold product of named
// guards g_p, "test p, restore, test q, ...").  k-fold products are subcube indicators fixing at most k(n-1) free pairs,
// so k >= ceil(M/(n-1)) always (degree bound).  Exhaustive for n <= 7 (all 2^M states); random sample for n = 8.
// Build: g++ -O2 -std=c++17 -o rdepth restore_depth_20261026.cpp ; run: ./rdepth [maxExhaustive=7] [samples8=10000]
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace std;
typedef unsigned long long u64; typedef unsigned int u32;

static int n, M; static int adj[9][9]; static int pairIndex[9][9]; static int order[9];
static vector<u32> masks;                       // free-pair masks of the Hamiltonian paths of the current state

static void dfs(int pos, int used, u32 mask) {
    if (pos == n) { masks.push_back(mask); return; }
    int last = order[pos - 1];
    for (int v = 0; v < n; v++) if (!(used >> v & 1) && adj[last][v]) {
        int k = pairIndex[last][v]; order[pos] = v; dfs(pos + 1, used | 1 << v, k >= 0 ? (mask | (1u << k)) : mask);
    }
}
static int depthOf() {                           // least k in {1,2,3,4} covering all M free pairs, or 5 if none
    u32 full = (1u << M) - 1;
    sort(masks.begin(), masks.end()); masks.erase(unique(masks.begin(), masks.end()), masks.end());
    for (u32 m : masks) if (m == full) return 1;
    vector<char> cover(1u << M, 0);              // cover[s] = 1 iff some path mask contains s (superset closure)
    for (u32 m : masks) cover[m] = 1;
    for (int b = 0; b < M; b++) for (u32 s = 0; s < (1u << M); s++) if ((s >> b & 1) && cover[s]) cover[s ^ (1u << b)] = 1;
    // depth 2: exists a with cover[full & ~a]
    for (u32 a : masks) if (cover[full & ~a]) return 2;
    // depth 3: exists a, b with cover[full & ~(a|b)]
    vector<u32> un; un.reserve(masks.size() * (masks.size() + 1) / 2);
    for (size_t i = 0; i < masks.size(); i++) for (size_t j = i; j < masks.size(); j++) { u32 ab = masks[i] | masks[j]; if (cover[full & ~ab]) return 3; un.push_back(ab); }
    sort(un.begin(), un.end()); un.erase(unique(un.begin(), un.end()), un.end());
    for (u32 ab : un) for (u32 c : masks) if (cover[full & ~(ab | c)]) return 4;
    return 5;
}
int main(int argc, char** argv) {
    int maxEx = argc > 1 ? atoi(argv[1]) : 7; int samples8 = argc > 2 ? atoi(argv[2]) : 10000;
    u64 rng = 0x9E3779B97F4A7C15ULL; auto rnd = [&]() { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng * 0x2545F4914F6CDD1DULL; };
    for (n = 4; n <= 8; n++) {
        vector<pair<int,int>> pairs; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) pairIndex[i][j] = -1;
        for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { pairIndex[i][j] = pairIndex[j][i] = (int)pairs.size(); pairs.push_back({i, j}); }
        M = (int)pairs.size(); u64 N = 1ULL << M; bool exhaustive = (n <= maxEx);
        u64 total = exhaustive ? N : (u64)samples8; vector<u64> hist(6, 0); int maxd = 0;
        for (u64 t = 0; t < total; t++) {
            u64 code = exhaustive ? t : (rnd() & (N - 1));
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj[i][j] = 0;
            for (int i = 0; i + 1 < n; i++) adj[i][i + 1] = 1;
            for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; if (code >> k & 1) adj[i][j] = 1; else adj[j][i] = 1; }
            masks.clear(); for (int v = 0; v < n; v++) { order[0] = v; dfs(1, 1 << v, 0); }
            int d = depthOf(); hist[d]++; if (d > maxd) maxd = d;
        }
        int bound = (M + n - 2) / (n - 1);
        printf("n = %d: M = %d, degree bound ceil(M/(n-1)) = %d, %s %llu states: depth histogram 1:%llu 2:%llu 3:%llu 4:%llu >=5:%llu; restore depth d(%d) = %d%s\n",
               n, M, bound, exhaustive ? "exhaustive" : "sampled", (unsigned long long)total, (unsigned long long)hist[1], (unsigned long long)hist[2], (unsigned long long)hist[3], (unsigned long long)hist[4], (unsigned long long)hist[5], n, maxd, (maxd == bound ? " (= bound)" : ""));
        fflush(stdout);
    }
    return 0;
}
