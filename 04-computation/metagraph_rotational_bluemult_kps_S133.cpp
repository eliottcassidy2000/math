// metagraph_rotational_bluemult_kps_S133.cpp
// (kind-pasteur-2026-07-26-S133)
//
// Exact blue multiplicities of the rotational (circulant) tournaments
// T_p inside the blue sub-cube, p in {5, 7, 11}, plus their H, |Aut|,
// tc.  Data so far: bm(T3)=1=tc, bm(T5)=3=tc, bm(T7)=3<9=tc.
// Question (THM-2454 open): does bm(T_p) stall at 3 for all p >= 5,
// making T5 the unique pure rotational atom automatically?
//
// Method: iterate the 2^e blue cube at n = p; filter to regular score
// vectors (all (p-1)/2) -- T_p is regular, so only regular tilings can
// be isomorphic to it; exact iso test on the survivors.
//
// Build: g++ -O2 -std=c++17 -o rotbm 04-computation/metagraph_rotational_bluemult_kps_S133.cpp
// Run:   ./rotbm
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
using namespace std;

static int N;
typedef array<uint16_t, 12> Tour;
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

// iso via backtracking with score classes (all equal here, so plain);
// prune with out-neighborhood score multisets first (cheap 2nd invariant)
static bool isomorphic(const Tour& A, const Tour& B) {
    // order: fix image of vertex 0 among all, then extend
    int mapping[12]; bool used[12];
    memset(used, 0, sizeof(used));
    struct BT {
        const Tour *A, *B; int n; int* mapping; bool* used;
        bool go(int k) {
            if (k == n) return true;
            for (int j = 0; j < n; ++j) {
                if (used[j]) continue;
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
    } bt{&A, &B, N, mapping, used};
    return bt.go(0);
}

static long long aut_size(const Tour& A) {
    int mapping[12]; bool used[12];
    memset(used, 0, sizeof(used));
    long long cnt = 0;
    struct BT {
        const Tour* A; int n; int* mapping; bool* used; long long* cnt;
        void go(int k) {
            if (k == n) { ++*cnt; return; }
            for (int j = 0; j < n; ++j) {
                if (used[j]) continue;
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
    } bt{&A, N, mapping, used, &cnt};
    bt.go(0);
    return cnt;
}

static long long ham(const Tour& A) {
    int full = (1 << N) - 1;
    vector<long long> dp((size_t)(1 << N) * N, 0);
    for (int i = 0; i < N; ++i) dp[(size_t)(1 << i) * N + i] = 1;
    for (int mask = 1; mask <= full; ++mask)
        for (int last = 0; last < N; ++last) {
            long long c = dp[(size_t)mask * N + last];
            if (!c) continue;
            uint16_t nx = A[last] & (uint16_t)(~mask);
            while (nx) {
                int v = __builtin_ctz(nx); nx &= (uint16_t)(nx - 1);
                dp[(size_t)(mask | (1 << v)) * N + v] += c;
            }
        }
    long long t = 0;
    for (int i = 0; i < N; ++i) t += dp[(size_t)full * N + i];
    return t;
}

int main() {
    int ps[3] = {5, 7, 11};
    // residue sets: QR-type circulants used in canon: T5 {1,2}, T7 {1,2,4},
    // T11 {1,3,4,5,9} (QRs mod 11)
    for (int pi = 0; pi < 3; ++pi) {
        int p = ps[pi];
        N = p;
        build_tiles();
        int e = (int)orbits.size();
        Tour Tp;
        for (int i = 0; i < p; ++i) Tp[i] = 0;
        for (int i = 0; i < p; ++i)
            for (int j = 0; j < p; ++j) {
                if (i == j) continue;
                int d = ((j - i) % p + p) % p;
                bool in;
                if (p == 5) in = (d == 1 || d == 2);
                else if (p == 7) in = (d == 1 || d == 2 || d == 4);
                else in = (d == 1 || d == 3 || d == 4 || d == 5 || d == 9);
                if (in) Tp[i] |= (uint16_t)(1u << j);
            }
        long long H = ham(Tp), au = aut_size(Tp);
        int reg = (p - 1) / 2;
        uint64_t total = 1ULL << e;
        long long bm = 0, regcount = 0;
        for (uint64_t code = 0; code < total; ++code) {
            Tour A; tour_from_code(code, A);
            bool isreg = true;
            for (int i = 0; i < N && isreg; ++i)
                if (__builtin_popcount(A[i]) != reg) isreg = false;
            if (!isreg) continue;
            ++regcount;
            if (isomorphic(A, Tp)) ++bm;
        }
        printf("p=%d: H=%lld |Aut|=%lld tc=%lld  regular blue tilings=%lld  "
               "blue-mult(T_p)=%lld  %s\n",
               p, H, au, H / au, regcount, bm,
               (bm == H / au) ? "PURE" : "IMPURE");
        fflush(stdout);
    }
    printf("DONE\n");
    return 0;
}
