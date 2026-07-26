// metagraph_pureblue_n12_pruned_kps_S133.cpp
// (kind-pasteur-2026-07-26-S133; pruned pure-blue detector, n <= 12)
//
// Finds ALL pure-blue classes without a full class census.  For each
// grid-symmetric tiling (2^e, e = floor((n-1)^2/4); 2^30 at n = 12):
// sample alternative Hamiltonian paths of its tournament by insertion
// sort under several fixed vertex orders; if any sampled path induces
// a non-grid-symmetric tiling, the class is impure -> skip (cheap).
// Survivors are deduplicated by exact isomorphism and settled by
// EXHAUSTIVE Hamiltonian-path enumeration: pure-blue iff every path's
// tiling is grid-symmetric.  Exactness: the sampler can only create
// false survivors, never false rejections (a gridsym tiling's class
// containing ANY non-gridsym tiling is impure by definition), and
// false survivors are eliminated in the exact phase.
//
// Grid-symmetry test at tournament level (THM-644(a)): with vertices
// labeled by path position (label n = path start ... label 1 = end),
// the tiling is grid-symmetric iff  arc(i -> j) <=> arc(n+1-j -> n+1-i)
// for all i, j (rho(i) = n+1-i is an anti-automorphism).
//
// Controls: n = 9, 10, 11 must reproduce the exact known inventories
// (6, 4, 9 classes; THM-2444/2454).  Prediction under test at n = 12
// (THM-2453/2454): exactly 6 pure-blue classes, all rigid,
// H-multiset [1, 9, 9, 9, 9, 81].
//
// Build: g++ -O2 -std=c++17 -o pb12 04-computation/metagraph_pureblue_n12_pruned_kps_S133.cpp
// Run:   ./pb12 [n ...]      (default 9 10 11 12)
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
using namespace std;

static int N;
typedef array<uint16_t, 12> Tour;   // adj[i] bit j: i beats j (0-indexed)

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

// gridsym test for the tiling induced by path perm: pos[k] = vertex at
// path position k (position 0 = path start = label n ... position n-1 = label 1).
// label of vertex v: L(v) = n - pos-index.  Condition: for all label pairs
// (i, j): arc(i->j) <=> arc(n+1-j -> n+1-i).  In pos indices: vertex at
// position a has label n-a; n+1-(n-a) = a+1 -> label n-a maps to position
// index (n-1)-a.  So: arc(pos[a] -> pos[b]) <=> arc(pos[n-1-b] -> pos[n-1-a]).
static inline bool tiling_gridsym(const Tour& A, const int pos[12]) {
    for (int a = 0; a < N; ++a) {
        for (int b = a + 1; b < N; ++b) {
            int u = pos[a], v = pos[b];
            int ru = pos[N - 1 - b], rv = pos[N - 1 - a];
            int arc1 = (A[u] >> v) & 1;
            int arc2 = (A[ru] >> rv) & 1;
            if (arc1 != arc2) return false;
        }
    }
    return true;
}

// insertion-sort Hamiltonian path under a given vertex order:
// maintain path as list; insert next vertex at first position where it
// beats the occupant (standard tournament Ham path construction).
static inline void ham_path_insertion(const Tour& A, const int* order, int pos[12]) {
    int len = 0;
    for (int oi = 0; oi < N; ++oi) {
        int v = order[oi];
        int at = len;               // default: append at end
        for (int k = 0; k < len; ++k) {
            if ((A[v] >> pos[k]) & 1) { at = k; break; }  // v beats pos[k]
        }
        for (int k = len; k > at; --k) pos[k] = pos[k - 1];
        pos[at] = v;
        ++len;
    }
}

// several fixed vertex orders for the sampler
static const int ORDERS12[6][12] = {
    {0,1,2,3,4,5,6,7,8,9,10,11},
    {11,10,9,8,7,6,5,4,3,2,1,0},
    {5,8,2,11,0,7,3,10,1,6,4,9},
    {6,1,9,4,0,11,2,8,5,10,3,7},
    {3,7,11,1,5,9,0,4,8,2,6,10},
    {9,2,6,0,10,4,8,1,11,5,3,7},
};

// ---- exact phase helpers (WL, iso, aut, exhaustive paths) ----
static inline uint64_t mixh(uint64_t h, uint64_t v) {
    h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}
static void wl_colors(const Tour& A, int col[12]) {
    int nc[12];
    for (int i = 0; i < N; ++i) col[i] = __builtin_popcount(A[i]);
    for (int r = 0; r < 3; ++r) {
        for (int i = 0; i < N; ++i) {
            int outs[12], ins[12]; int no = 0, ni = 0;
            for (int j = 0; j < N; ++j) {
                if (j == i) continue;
                if ((A[i] >> j) & 1) outs[no++] = col[j];
                else ins[ni++] = col[j];
            }
            sort(outs, outs + no); sort(ins, ins + ni);
            uint64_t h = (uint64_t)col[i] * 1315423911ULL;
            for (int k = 0; k < no; ++k) h = mixh(h, outs[k] * 2);
            for (int k = 0; k < ni; ++k) h = mixh(h, ins[k] * 2 + 1);
            nc[i] = (int)(h & 0x3fffffff);
        }
        int sv[12]; memcpy(sv, nc, sizeof(int) * N); sort(sv, sv + N);
        int uq[12]; int nu = 0;
        for (int i = 0; i < N; ++i)
            if (i == 0 || sv[i] != sv[i - 1]) uq[nu++] = sv[i];
        for (int i = 0; i < N; ++i)
            col[i] = (int)(lower_bound(uq, uq + nu, nc[i]) - uq);
    }
}
struct IsoCtx {
    const Tour *A, *B;
    int colA[12], colB[12], order[12], mapping[12];
    bool used[12]; bool count_mode; long long count;
};
static bool iso_bt(IsoCtx& c, int k) {
    if (k == N) { if (c.count_mode) { c.count++; return false; } return true; }
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
static void iso_prep(IsoCtx& c, const Tour& A, const Tour& B) {
    c.A = &A; c.B = &B;
    wl_colors(A, c.colA); wl_colors(B, c.colB);
    int cs[12] = {0};
    for (int i = 0; i < N; ++i) cs[c.colA[i]]++;
    for (int i = 0; i < N; ++i) c.order[i] = i;
    sort(c.order, c.order + N, [&](int a, int b) {
        if (cs[c.colA[a]] != cs[c.colA[b]]) return cs[c.colA[a]] < cs[c.colA[b]];
        return a < b;
    });
    memset(c.used, 0, sizeof(c.used)); c.count = 0;
}
static bool isomorphic(const Tour& A, const Tour& B) {
    IsoCtx c; c.count_mode = false; iso_prep(c, A, B);
    int ca[12] = {0}, cb[12] = {0};
    for (int i = 0; i < N; ++i) { ca[c.colA[i]]++; cb[c.colB[i]]++; }
    for (int i = 0; i < N; ++i) if (ca[i] != cb[i]) return false;
    return iso_bt(c, 0);
}
static long long aut_size(const Tour& A) {
    IsoCtx c; c.count_mode = true; iso_prep(c, A, A);
    iso_bt(c, 0); return c.count;
}

// exhaustive path enumeration with purity check; returns H, sets pure
static long long path_count;
static bool all_pure;
static int gpos[12];
static uint16_t gused;
static const Tour* gA;
static void path_bt(int depth) {
    if (path_count > 2000000) { all_pure = false; return; }  // safety cap
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
        if (!all_pure && path_count > 0) {
            // keep counting H even after impurity found (H needed for report)
        }
    }
}

int main(int argc, char** argv) {
    vector<int> ns;
    for (int i = 1; i < argc; ++i) ns.push_back(atoi(argv[i]));
    if (ns.empty()) ns = {9, 10, 11, 12};
    for (int n : ns) {
        N = n;
        build_tiles();
        int e = (int)orbits.size();
        if (e != (N - 1) * (N - 1) / 4) { printf("ORBIT COUNT WRONG\n"); return 1; }
        uint64_t total = 1ULL << e;
        vector<Tour> survivors;
        uint64_t surv_raw = 0;
        int order_buf[12];
        for (uint64_t code = 0; code < total; ++code) {
            Tour A; tour_from_code(code, A);
            bool alive = true;
            for (int s = 0; s < 6 && alive; ++s) {
                for (int k = 0; k < N; ++k) order_buf[k] = ORDERS12[s][k] % N;
                // dedupe order entries for n < 12: use modulo then stable unique
                bool seen[12] = {false};
                int ob[12], nb = 0;
                for (int k = 0; k < 12 && nb < N; ++k) {
                    int v = ORDERS12[s][k] % N;
                    if (!seen[v]) { seen[v] = true; ob[nb++] = v; }
                }
                int pos[12];
                ham_path_insertion(A, ob, pos);
                if (!tiling_gridsym(A, pos)) alive = false;
            }
            if (!alive) continue;
            ++surv_raw;
            bool dup = false;
            for (auto& S : survivors)
                if (isomorphic(A, S)) { dup = true; break; }
            if (!dup) survivors.push_back(A);
            if ((code & 0x3FFFFFF) == 0x3FFFFFF)
                fprintf(stderr, "n=%d: %llu/%llu, raw survivors %llu, classes %zu\n",
                        N, (unsigned long long)(code + 1), (unsigned long long)total,
                        (unsigned long long)surv_raw, survivors.size());
        }
        // exact phase
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
        printf("n=%d: blue cube 2^%d; raw survivors %llu; survivor classes %zu; "
               "PURE-BLUE=%d\n", N, e, (unsigned long long)surv_raw,
               survivors.size(), pure);
        for (auto& t : inv)
            printf("    pure-blue class: H=%lld |Aut|=%lld tc=%lld\n",
                   t[0], t[1], t[2]);
        fflush(stdout);
    }
    printf("ALL RUNS COMPLETE\n");
    return 0;
}
