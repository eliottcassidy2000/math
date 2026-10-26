// legality_language_reconstruction_20261026.cpp -- opus-2026-10-26-S61
// The named Hamilton-path legality language L(T) = {permutations p : p is a directed Hamilton path of T}.
// Theorem (note, section on re-marking guards): for every tournament T and every pair {u, v}, the arc of T between
// u and v lies on a Hamiltonian path of T, or the reversed arc lies on a Hamiltonian path of the tournament T^(uv)
// obtained by flipping it; hence T -> L(T) is injective on labelled tournaments (legal re-marking guards reconstruct
// every tournament at every size).  This program checks the statement exhaustively:
//   mode A (n <= 7): all 2^C(n,2) labelled tournaments; for every tournament and every pair, test the pair directly
//                    (bright in T, or bright after the flip); count dark arcs (arcs on no Hamiltonian path) and the
//                    tournaments with dark arcs, split by strong / non-strong; also hash L(T) and look for collisions.
//   mode B (n = 8):  the marked-path carrier X_8 (tournaments containing 0->1->...->7): hash L(T) over all 2^21 members,
//                    look for collisions (injectivity on X_8), and count dark arcs.
// Build: g++ -O2 -std=c++17 -o llr legality_language_reconstruction_20261026.cpp ; run: ./llr [maxA=7] [doB=1]
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace std;
typedef unsigned long long u64;

static int n; static int adj[9][9];                 // adj[i][j] = 1 iff i -> j
static u64 arcsUsed;                                 // bit i*n+j set iff arc i->j lies on some Hamiltonian path
static u64 langHash; static u64 pathCount;
static int order[9];

static inline u64 mix(u64 h, u64 x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return h ^ (x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2)); }

static void dfs(int pos, int used, u64 code) {
    if (pos == n) {
        pathCount++;
        for (int i = 0; i + 1 < n; i++) arcsUsed |= 1ULL << (order[i] * n + order[i + 1]);
        langHash += mix(0x1234567ULL, code);          // commutative accumulation: the multiset of path codes
        return;
    }
    int last = order[pos - 1];
    for (int v = 0; v < n; v++) if (!(used >> v & 1) && adj[last][v]) { order[pos] = v; dfs(pos + 1, used | 1 << v, code * 9 + v + 1); }
}
static void language() { arcsUsed = 0; langHash = 0; pathCount = 0; for (int v = 0; v < n; v++) { order[0] = v; dfs(1, 1 << v, v + 1); } }
static bool strongT() {                              // reachability from 0 and to 0
    int r = 1, prev = 0; while (r != prev) { prev = r; for (int i = 0; i < n; i++) if (r >> i & 1) for (int j = 0; j < n; j++) if (adj[i][j]) r |= 1 << j; }
    if (r != (1 << n) - 1) return false;
    r = 1; prev = 0; while (r != prev) { prev = r; for (int i = 0; i < n; i++) if (r >> i & 1) for (int j = 0; j < n; j++) if (adj[j][i]) r |= 1 << j; }
    return r == (1 << n) - 1;
}

int main(int argc, char** argv) {
    int maxA = argc > 1 ? atoi(argv[1]) : 7; int doB = argc > 2 ? atoi(argv[2]) : 1;
    vector<pair<int,int>> pairs;
    for (n = 3; n <= maxA; n++) {
        pairs.clear(); for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) pairs.push_back({i, j});
        int M = (int)pairs.size(); u64 total = 1ULL << M;
        u64 withDark = 0, strongWithDark = 0, strongCnt = 0, darkArcs = 0, maxDark = 0, violations = 0;
        vector<u64> hashes; hashes.reserve(total);
        for (u64 code = 0; code < total; code++) {
            for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; adj[i][j] = (code >> k) & 1; adj[j][i] = 1 - adj[i][j]; }
            language(); u64 used = arcsUsed, h = langHash; hashes.push_back(h);
            bool st = strongT(); if (st) strongCnt++;
            int dark = 0;
            for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; int a = adj[i][j] ? i : j, b = adj[i][j] ? j : i;
                if (!(used >> (a * n + b) & 1)) {
                    dark++;
                    // direct test: flip the pair and check that the reversed arc is on a Hamiltonian path
                    adj[i][j] ^= 1; adj[j][i] ^= 1; language();
                    if (!(arcsUsed >> (b * n + a) & 1)) violations++;
                    adj[i][j] ^= 1; adj[j][i] ^= 1;
                }
            }
            if (dark) { withDark++; if (st) strongWithDark++; darkArcs += dark; if ((u64)dark > maxDark) maxDark = dark; }
        }
        sort(hashes.begin(), hashes.end()); u64 coll = 0; for (size_t i = 1; i < hashes.size(); i++) if (hashes[i] == hashes[i - 1]) coll++;
        printf("n = %d: %llu tournaments (%llu strong); with a dark arc %llu (strong ones %llu); dark arcs total %llu, max per tournament %llu; pair-theorem violations %llu; language-hash collisions %llu\n",
               n, total, strongCnt, withDark, strongWithDark, darkArcs, maxDark, violations, coll);
    }
    if (doB) {
        n = 8; pairs.clear(); for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) pairs.push_back({i, j});   // free pairs of X_8
        int M = (int)pairs.size(); u64 total = 1ULL << M;
        vector<u64> hashes; hashes.reserve(total); u64 withDark = 0, darkArcs = 0, strongCnt = 0, strongWithDark = 0;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj[i][j] = 0;
        for (int i = 0; i + 1 < n; i++) { adj[i][i + 1] = 1; adj[i + 1][i] = 0; }
        for (u64 code = 0; code < total; code++) {
            for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; adj[i][j] = (code >> k) & 1; adj[j][i] = 1 - adj[i][j]; }
            language(); hashes.push_back(langHash); bool st = strongT(); if (st) strongCnt++;
            int dark = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (adj[i][j] && !(arcsUsed >> (i * n + j) & 1)) dark++;
            if (dark) { withDark++; darkArcs += dark; if (st) strongWithDark++; }
        }
        sort(hashes.begin(), hashes.end()); u64 coll = 0; for (size_t i = 1; i < hashes.size(); i++) if (hashes[i] == hashes[i - 1]) coll++;
        printf("X_8: %llu tournaments containing the standard path (%llu strong); with a dark arc %llu (strong %llu); dark arcs total %llu; language-hash collisions %llu\n", total, strongCnt, withDark, strongWithDark, darkArcs, coll);
    }
    return 0;
}
