// metagraph_pureblue_n13_threaded_kps_S134.cpp
// (kind-pasteur-2026-07-26-S134; threaded pruned pure-blue detector, n <= 13)
//
// Same sound design as metagraph_pureblue_n12_pruned_kps_S133.cpp
// (sampled alternative-Hamiltonian-path impurity witness; only false
// SURVIVORS possible, eliminated by exhaustive path enumeration), with
// std::thread range partitioning for the 2^36 blue cube at n = 13.
// Controls: n = 11, 12 must reproduce the known inventories (9; 6).
// Prediction under test (THM-2454 alphabet law, pre-registered):
//   pure-blue(13) = 13 = 10 rigid [1,3,9,9,9,9,27,27,27,81]
//                      + (15,5,3) + 2 x (135,45,3).
// A 14th class or any inventory deviation = a NEW pure-blue atom.
//
// Build: g++ -O2 -std=c++17 -pthread -o pb13 04-computation/metagraph_pureblue_n13_threaded_kps_S134.cpp
// Run:   ./pb13 [n ...]        (default 11 12 13)
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
using namespace std;

static int N;
typedef array<uint16_t, 13> Tour;
struct Tile { int x, y; };
static vector<Tile> tiles;
static vector<array<int, 2>> orbits;

static void build_tiles() {
    tiles.clear(); orbits.clear();
    for (int y = 1; y <= N - 2; ++y)
        for (int x = N; x >= y + 2; --x)
            tiles.push_back({x, y});
    int m = (int)tiles.size();
    vector<int> where(64 * 64, -1);
    for (int i = 0; i < m; ++i) where[tiles[i].x * 64 + tiles[i].y] = i;
    vector<char> seen(m, 0);
    for (int i = 0; i < m; ++i) {
        if (seen[i]) continue;
        int x = tiles[i].x, y = tiles[i].y;
        int j = where[(N - y + 1) * 64 + (N - x + 1)];
        seen[i] = 1; seen[j] = 1;
        orbits.push_back({min(i, j), max(i, j)});
    }
}

static inline void tour_from_code(uint64_t code, Tour& A) {
    for (int i = 0; i < N; ++i) A[i] = 0;
    for (int k = 1; k < N; ++k) A[k] |= (uint16_t)(1u << (k - 1));
    int e = (int)orbits.size();
    for (int oi = 0; oi < e; ++oi) {
        int v = (int)((code >> oi) & 1u);
        for (int t = 0; t < 2; ++t) {
            int idx = orbits[oi][t];
            int x = tiles[idx].x - 1, y = tiles[idx].y - 1;
            if (v == 0) A[x] |= (uint16_t)(1u << y);
            else        A[y] |= (uint16_t)(1u << x);
        }
    }
}

static inline bool tiling_gridsym(const Tour& A, const int pos[13]) {
    for (int a = 0; a < N; ++a)
        for (int b = a + 1; b < N; ++b) {
            int u = pos[a], v = pos[b];
            int ru = pos[N - 1 - b], rv = pos[N - 1 - a];
            if ((((A[u] >> v) & 1)) != ((A[ru] >> rv) & 1)) return false;
        }
    return true;
}

static inline void ham_path_insertion(const Tour& A, const int* order, int pos[13]) {
    int len = 0;
    for (int oi = 0; oi < N; ++oi) {
        int v = order[oi];
        int at = len;
        for (int k = 0; k < len; ++k)
            if ((A[v] >> pos[k]) & 1) { at = k; break; }
        for (int k = len; k > at; --k) pos[k] = pos[k - 1];
        pos[at] = v;
        ++len;
    }
}

static const int RAW_ORDERS[6][13] = {
    {0,1,2,3,4,5,6,7,8,9,10,11,12},
    {12,11,10,9,8,7,6,5,4,3,2,1,0},
    {5,8,2,11,0,7,3,10,1,6,4,9,12},
    {6,1,9,4,0,11,2,8,5,10,3,7,12},
    {3,7,11,1,5,9,0,4,8,2,6,10,12},
    {9,2,6,0,10,4,8,1,11,5,3,7,12},
};
static int ORDERS[6][13];
static void build_orders() {
    for (int s = 0; s < 6; ++s) {
        bool seen[13] = {false};
        int nb = 0;
        for (int k = 0; k < 13 && nb < N; ++k) {
            int v = RAW_ORDERS[s][k] % N;
            if (!seen[v]) { seen[v] = true; ORDERS[s][nb++] = v; }
        }
        for (int v = 0; v < N && nb < N; ++v)
            if (!seen[v]) { seen[v] = true; ORDERS[s][nb++] = v; }
    }
}

// ---- exact-phase machinery (single-threaded post-pass) ----
static bool isomorphic(const Tour& A, const Tour& B) {
    int sa[13], sb[13];
    for (int i = 0; i < N; ++i) sa[i] = __builtin_popcount(A[i]);
    for (int i = 0; i < N; ++i) sb[i] = __builtin_popcount(B[i]);
    int ta[13], tb[13];
    memcpy(ta, sa, sizeof(int) * N); memcpy(tb, sb, sizeof(int) * N);
    sort(ta, ta + N); sort(tb, tb + N);
    for (int i = 0; i < N; ++i) if (ta[i] != tb[i]) return false;
    int mapping[13]; bool used[13];
    memset(used, 0, sizeof(used));
    struct BT {
        const Tour *A, *B; int n; const int *sa, *sb; int* mapping; bool* used;
        bool go(int k) {
            if (k == n) return true;
            for (int j = 0; j < n; ++j) {
                if (used[j] || sb[j] != sa[k]) continue;
                bool ok = true;
                for (int k2 = 0; k2 < k && ok; ++k2) {
                    if ((((*A)[k] >> k2) & 1) != (((*B)[j] >> mapping[k2]) & 1)) ok = false;
                    else if ((((*A)[k2] >> k) & 1) != (((*B)[mapping[k2]] >> j) & 1)) ok = false;
                }
                if (!ok) continue;
                mapping[k] = j; used[j] = true;
                if (go(k + 1)) return true;
                used[j] = false;
            }
            return false;
        }
    } bt{&A, &B, N, sa, sb, mapping, used};
    return bt.go(0);
}

static long long aut_size(const Tour& A) {
    int sa[13];
    for (int i = 0; i < N; ++i) sa[i] = __builtin_popcount(A[i]);
    int mapping[13]; bool used[13];
    memset(used, 0, sizeof(used));
    long long cnt = 0;
    struct BT {
        const Tour* A; int n; const int* sa; int* mapping; bool* used; long long* cnt;
        void go(int k) {
            if (k == n) { ++*cnt; return; }
            for (int j = 0; j < n; ++j) {
                if (used[j] || sa[j] != sa[k]) continue;
                bool ok = true;
                for (int k2 = 0; k2 < k && ok; ++k2) {
                    if ((((*A)[k] >> k2) & 1) != (((*A)[j] >> mapping[k2]) & 1)) ok = false;
                    else if ((((*A)[k2] >> k) & 1) != (((*A)[mapping[k2]] >> j) & 1)) ok = false;
                }
                if (!ok) continue;
                mapping[k] = j; used[j] = true;
                go(k + 1);
                used[j] = false;
            }
        }
    } bt{&A, N, sa, mapping, used, &cnt};
    bt.go(0);
    return cnt;
}

static long long path_count;
static bool all_pure;
static int gpos[13];
static uint16_t gused;
static const Tour* gA;
static void path_bt(int depth) {
    if (path_count > 4000000) { all_pure = false; return; }
    if (depth == N) {
        ++path_count;
        if (all_pure && !tiling_gridsym(*gA, gpos)) all_pure = false;
        return;
    }
    for (int v = 0; v < N; ++v) {
        if ((gused >> v) & 1) continue;
        if (depth > 0 && !(((*gA)[gpos[depth - 1]] >> v) & 1)) continue;
        gpos[depth] = v; gused |= (uint16_t)(1u << v);
        path_bt(depth + 1);
        gused &= (uint16_t)~(1u << v);
    }
}

int main(int argc, char** argv) {
    vector<int> ns;
    for (int i = 1; i < argc; ++i) ns.push_back(atoi(argv[i]));
    if (ns.empty()) ns = {11, 12, 13};
    for (int n : ns) {
        N = n;
        build_tiles();
        build_orders();
        int e = (int)orbits.size();
        if (e != (N - 1) * (N - 1) / 4) { printf("ORBIT COUNT WRONG\n"); return 1; }
        uint64_t total = 1ULL << e;
        unsigned nthreads = max(2u, thread::hardware_concurrency()) - 1;
        if (total < (1ULL << 22)) nthreads = 1;
        vector<vector<Tour>> tl_surv(nthreads);
        vector<uint64_t> tl_raw(nthreads, 0);
        atomic<uint64_t> progress{0};
        auto worker = [&](unsigned tid) {
            uint64_t lo = total / nthreads * tid;
            uint64_t hi = (tid == nthreads - 1) ? total : total / nthreads * (tid + 1);
            int pos[13];
            for (uint64_t code = lo; code < hi; ++code) {
                Tour A; tour_from_code(code, A);
                bool alive = true;
                for (int s = 0; s < 6 && alive; ++s) {
                    ham_path_insertion(A, ORDERS[s], pos);
                    if (!tiling_gridsym(A, pos)) alive = false;
                }
                if (!alive) continue;
                ++tl_raw[tid];
                bool dup = false;
                for (auto& S : tl_surv[tid])
                    if (isomorphic(A, S)) { dup = true; break; }
                if (!dup) tl_surv[tid].push_back(A);
                if ((code & 0xFFFFFFF) == 0xFFFFFFF) {
                    uint64_t done = progress.fetch_add(0x10000000) + 0x10000000;
                    fprintf(stderr, "n=%d: ~%llu/%llu done, tid%u raw %llu\n",
                            N, (unsigned long long)done, (unsigned long long)total,
                            tid, (unsigned long long)tl_raw[tid]);
                }
            }
        };
        vector<thread> pool;
        for (unsigned t = 0; t < nthreads; ++t) pool.emplace_back(worker, t);
        for (auto& th : pool) th.join();
        // merge
        vector<Tour> survivors;
        uint64_t surv_raw = 0;
        for (unsigned t = 0; t < nthreads; ++t) {
            surv_raw += tl_raw[t];
            for (auto& S : tl_surv[t]) {
                bool dup = false;
                for (auto& G : survivors)
                    if (isomorphic(S, G)) { dup = true; break; }
                if (!dup) survivors.push_back(S);
            }
        }
        int pure = 0;
        vector<array<long long, 3>> inv;
        for (auto& S : survivors) {
            gA = &S; path_count = 0; all_pure = true; gused = 0;
            path_bt(0);
            if (!all_pure) continue;
            long long H = path_count;
            long long au = aut_size(S);
            if (H % au != 0) { printf("TC NOT INTEGER\n"); return 1; }
            inv.push_back({H, au, H / au});
            ++pure;
        }
        sort(inv.begin(), inv.end());
        printf("n=%d: blue cube 2^%d; threads %u; raw survivors %llu; "
               "survivor classes %zu; PURE-BLUE=%d\n", N, e, nthreads,
               (unsigned long long)surv_raw, survivors.size(), pure);
        for (auto& t : inv)
            printf("    pure-blue class: H=%lld |Aut|=%lld tc=%lld\n",
                   t[0], t[1], t[2]);
        fflush(stdout);
    }
    printf("ALL RUNS COMPLETE\n");
    return 0;
}
