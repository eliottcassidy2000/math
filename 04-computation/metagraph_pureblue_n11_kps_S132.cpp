// metagraph_pureblue_n11_kps_S132.cpp
// (kind-pasteur-2026-07-26-S132; pure-blue census engine, n <= 11)
//
// Decide pure-blue(n) for n = 9, 10 (controls: 6 and 4, THM-2444) and
// n = 11 (new point).  Blue sub-cube method (kps-S66): enumerate the
// 2^e grid-symmetric tilings, e = floor((n-1)^2/4); group by tournament
// isomorphism (WL fingerprint buckets + exact colour-pruned iso);
// a touched class is pure-blue iff blue-mult == tc = H/|Aut| (LEM-003).
// Repaired conjecture under test (THM-2444/2450):
//   pure-blue(11) = rigid-SC-towers(11) + 1.
//
// Build:  g++ -O2 -std=c++17 -o metagraph_pureblue_n11_kps_S132 \
//              04-computation/metagraph_pureblue_n11_kps_S132.cpp
// Run:    ./metagraph_pureblue_n11_kps_S132 [n ...]   (default: 9 10 11)
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>
using namespace std;

static int N;                       // current n <= 11
typedef array<uint16_t, 11> Tour;   // adj[i] bit j set: i beats j

struct Tile { int x, y; };          // 1-indexed, y < x, x-y >= 2

static vector<Tile> tiles;
static vector<array<int, 2>> orbits;   // tile indices per grid orbit

static void build_tiles() {
    tiles.clear();
    for (int y = 1; y <= N - 2; ++y)
        for (int x = N; x >= y + 2; --x)
            tiles.push_back({x, y});
    // grid orbits under (x,y) -> (N-y+1, N-x+1)
    orbits.clear();
    int m = (int)tiles.size();
    vector<int> pos(1 << 10, -1);
    auto key = [](int x, int y) { return x * 32 + y; };
    vector<int> where(32 * 32, -1);
    for (int i = 0; i < m; ++i) where[key(tiles[i].x, tiles[i].y)] = i;
    vector<char> seen(m, 0);
    for (int i = 0; i < m; ++i) {
        if (seen[i]) continue;
        int x = tiles[i].x, y = tiles[i].y;
        int j = where[key(N - y + 1, N - x + 1)];
        seen[i] = 1; seen[j] = 1;
        orbits.push_back({min(i, j), max(i, j)});
    }
}

static inline void tour_from_code(uint32_t code, Tour& A) {
    for (int i = 0; i < N; ++i) A[i] = 0;
    for (int k = 1; k < N; ++k) A[k] |= (uint16_t)(1u << (k - 1));
    int e = (int)orbits.size();
    for (int oi = 0; oi < e; ++oi) {
        int v = (code >> oi) & 1u;
        for (int t = 0; t < 2; ++t) {
            int idx = orbits[oi][t];
            int x = tiles[idx].x - 1, y = tiles[idx].y - 1;
            if (v == 0) A[x] |= (uint16_t)(1u << y);
            else        A[y] |= (uint16_t)(1u << x);
        }
    }
}

// ---- WL colours (3 rounds), fingerprint hash ----
static inline uint64_t mix(uint64_t h, uint64_t v) {
    h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

static void wl_colors(const Tour& A, int col[11]) {
    int nc[11];
    for (int i = 0; i < N; ++i) col[i] = __builtin_popcount(A[i]);
    for (int r = 0; r < 3; ++r) {
        for (int i = 0; i < N; ++i) {
            int outs[11], ins[11]; int no = 0, ni = 0;
            for (int j = 0; j < N; ++j) {
                if (j == i) continue;
                if ((A[i] >> j) & 1) outs[no++] = col[j];
                else ins[ni++] = col[j];
            }
            sort(outs, outs + no); sort(ins, ins + ni);
            uint64_t h = (uint64_t)col[i] * 1315423911ULL;
            for (int k = 0; k < no; ++k) h = mix(h, outs[k] * 2 + 0);
            for (int k = 0; k < ni; ++k) h = mix(h, ins[k] * 2 + 1);
            nc[i] = (int)(h & 0x3fffffff);
        }
        // compress
        int sorted_vals[11];
        memcpy(sorted_vals, nc, sizeof(int) * N);
        sort(sorted_vals, sorted_vals + N);
        int uniq[11]; int nu = 0;
        for (int i = 0; i < N; ++i)
            if (i == 0 || sorted_vals[i] != sorted_vals[i - 1]) uniq[nu++] = sorted_vals[i];
        for (int i = 0; i < N; ++i)
            col[i] = (int)(lower_bound(uniq, uniq + nu, nc[i]) - uniq);
    }
}

static uint64_t fingerprint(const Tour& A) {
    int col[11];
    wl_colors(A, col);
    pair<int, int> per[11];
    for (int i = 0; i < N; ++i) per[i] = {col[i], __builtin_popcount(A[i])};
    sort(per, per + N);
    uint64_t h = 1469598103934665603ULL ^ (uint64_t)N;
    for (int i = 0; i < N; ++i) { h = mix(h, per[i].first); h = mix(h, per[i].second); }
    return h;
}

// ---- exact iso / aut via colour-pruned backtracking ----
struct IsoCtx {
    const Tour *A, *B;
    int colA[11], colB[11];
    int order[11];
    int mapping[11];
    bool used[11];
    bool count_mode;
    long long count;
};

static bool iso_bt(IsoCtx& c, int k) {
    if (k == N) {
        if (c.count_mode) { c.count++; return false; }
        return true;
    }
    int i = c.order[k];
    for (int j = 0; j < N; ++j) {
        if (c.used[j] || c.colB[j] != c.colA[i]) continue;
        bool ok = true;
        for (int k2 = 0; k2 < k && ok; ++k2) {
            int i2 = c.order[k2], j2 = c.mapping[i2];
            if ((((*c.A)[i] >> i2) & 1) != (((*c.B)[j] >> j2) & 1)) ok = false;
            else if ((((*c.A)[i2] >> i) & 1) != (((*c.B)[j2] >> j) & 1)) ok = false;
        }
        if (!ok) continue;
        c.mapping[i] = j; c.used[j] = true;
        if (iso_bt(c, k + 1)) return true;
        c.used[j] = false;
    }
    return false;
}

static void iso_prepare(IsoCtx& c, const Tour& A, const Tour& B) {
    c.A = &A; c.B = &B;
    wl_colors(A, c.colA);
    wl_colors(B, c.colB);
    // order by ascending colour-class size
    int csize[11] = {0};
    for (int i = 0; i < N; ++i) csize[c.colA[i]]++;
    for (int i = 0; i < N; ++i) c.order[i] = i;
    sort(c.order, c.order + N, [&](int a, int b) {
        if (csize[c.colA[a]] != csize[c.colA[b]]) return csize[c.colA[a]] < csize[c.colA[b]];
        return a < b;
    });
    memset(c.used, 0, sizeof(c.used));
    c.count = 0;
}

static bool isomorphic(const Tour& A, const Tour& B) {
    int sa[11], sb[11];
    for (int i = 0; i < N; ++i) sa[i] = __builtin_popcount(A[i]);
    for (int i = 0; i < N; ++i) sb[i] = __builtin_popcount(B[i]);
    int ta[11], tb[11];
    memcpy(ta, sa, sizeof(int) * N); memcpy(tb, sb, sizeof(int) * N);
    sort(ta, ta + N); sort(tb, tb + N);
    for (int i = 0; i < N; ++i) if (ta[i] != tb[i]) return false;
    IsoCtx c; c.count_mode = false;
    iso_prepare(c, A, B);
    // colour multisets must agree
    int ca[11] = {0}, cb[11] = {0};
    for (int i = 0; i < N; ++i) { ca[c.colA[i]]++; cb[c.colB[i]]++; }
    for (int i = 0; i < N; ++i) if (ca[i] != cb[i]) return false;
    return iso_bt(c, 0);
}

static long long aut_size(const Tour& A) {
    IsoCtx c; c.count_mode = true;
    iso_prepare(c, A, A);
    iso_bt(c, 0);
    return c.count;
}

// ---- Hamiltonian path count DP ----
static long long ham(const Tour& A) {
    int full = (1 << N) - 1;
    static vector<long long> dp;
    dp.assign((size_t)(1 << N) * N, 0);
    for (int i = 0; i < N; ++i) dp[(size_t)(1 << i) * N + i] = 1;
    for (int mask = 1; mask <= full; ++mask) {
        for (int last = 0; last < N; ++last) {
            long long cval = dp[(size_t)mask * N + last];
            if (!cval) continue;
            uint16_t nxts = A[last] & (uint16_t)(~mask);
            while (nxts) {
                int nx = __builtin_ctz(nxts);
                nxts &= (uint16_t)(nxts - 1);
                dp[(size_t)(mask | (1 << nx)) * N + nx] += cval;
            }
        }
    }
    long long tot = 0;
    for (int i = 0; i < N; ++i) tot += dp[(size_t)full * N + i];
    return tot;
}

struct Cls { Tour rep; uint64_t bm; };

int main(int argc, char** argv) {
    vector<int> ns;
    for (int i = 1; i < argc; ++i) ns.push_back(atoi(argv[i]));
    if (ns.empty()) ns = {9, 10, 11};
    for (int n : ns) {
        N = n;
        build_tiles();
        int e = (int)orbits.size();
        if (e != (N - 1) * (N - 1) / 4) { printf("ORBIT COUNT WRONG\n"); return 1; }
        unordered_map<uint64_t, vector<int>> buckets;
        vector<Cls> classes;
        buckets.reserve(1 << 20);
        uint64_t total = 1ULL << e;
        for (uint64_t code = 0; code < total; ++code) {
            Tour A;
            tour_from_code((uint32_t)code, A);
            uint64_t fp = fingerprint(A);
            auto& lst = buckets[fp];
            int hit = -1;
            for (int ci : lst)
                if (isomorphic(A, classes[ci].rep)) { hit = ci; break; }
            if (hit < 0) {
                classes.push_back({A, 0});
                lst.push_back((int)classes.size() - 1);
                hit = (int)classes.size() - 1;
            }
            classes[hit].bm++;
            if ((code & 0xFFFFFF) == 0xFFFFFF)
                fprintf(stderr, "n=%d: %llu / %llu, classes %zu\n",
                        N, (unsigned long long)(code + 1),
                        (unsigned long long)total, classes.size());
        }
        uint64_t chk = 0;
        int pure = 0;
        vector<array<long long, 3>> inventory;
        for (auto& c : classes) {
            chk += c.bm;
            long long H = ham(c.rep);
            long long au = aut_size(c.rep);
            if (H % au != 0) { printf("TC NOT INTEGER\n"); return 1; }
            long long tc = H / au;
            if ((long long)c.bm == tc) {
                pure++;
                inventory.push_back({H, au, tc});
            }
        }
        if (chk != total) { printf("BLUE MASS MISMATCH\n"); return 1; }
        sort(inventory.begin(), inventory.end());
        printf("n=%d: blue cube 2^%d; classes touched=%zu; PURE-BLUE=%d\n",
               N, e, classes.size(), pure);
        for (auto& t : inventory)
            printf("    pure-blue class: H=%lld |Aut|=%lld tc=%lld\n", t[0], t[1], t[2]);
        fflush(stdout);
    }
    printf("ALL CHECKS PASSED\n");
    return 0;
}
