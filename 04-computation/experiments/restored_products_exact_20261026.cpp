// restored_products_exact_20261026.cpp -- opus-2026-10-26-S62
// Exact spans of products of at most k named guards on the marked carrier X_n (tournaments on {0, ..., n-1} that contain
// the path 0 -> 1 -> ... -> n-1; the M = (n-1)(n-2)/2 free pairs {i, j}, j >= i+2, carry one bit each, bit 1 = arc i -> j).
// The named guard g_p of a permutation p is the indicator of "p is a Hamiltonian path".  On X_n it is the indicator of a
// subcube: every free pair stepped by p must carry the arc p[i] -> p[i+1], and a step against a marked arc makes the
// guard identically zero.  A product of guards is the indicator of the intersection of subcubes, i.e. the subcube of the
// union of the constraint patterns when they are consistent and zero otherwise; so the distinct nonzero products of at
// most k guards are exactly the consistent unions of at most k constraint patterns.  The program enumerates them, builds
// each indicator once, inserts it into a GF(2) basis (bitsets over the 2^M states) and prints the GF(2) rank.
// A full GF(2) rank 2^M certifies rational rank 2^M (a 2^M x 2^M minor that is odd is nonzero over Z), i.e. that k
// restored named tests span every function on X_n; with the Boolean-degree bound (THM-4809 B, session72) that gives
// d_n = k whenever k = ceil((n-2)/2).
// With a third argument a certificate is written, one line per basis vector: "mask ori q1 q2 q3" where mask / ori are the
// fixed free pairs and their orientations (pair t = the t-th pair (i, j) in lexicographic order, bit 1 = arc i -> j) and
// q1..q3 are lexicographic indices (0-based, among all n! permutations) of guards whose product is that subcube (-1 =
// unused).  An independent program can recompute the products from q1..q3 and check that the 2^M vectors are independent.
// Build: g++ -O2 -std=c++17 -o rpe restored_products_exact_20261026.cpp ; run: ./rpe n k [certificate_path]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <unordered_map>
using namespace std;
typedef unsigned long long u64; typedef unsigned int u32;

int main(int argc, char** argv) {
    int n = argc > 1 ? atoi(argv[1]) : 7, k = argc > 2 ? atoi(argv[2]) : 3; const char* cert = argc > 3 ? argv[3] : nullptr;
    int M = (n - 1) * (n - 2) / 2; if (n < 3 || M > 30) { fprintf(stderr, "need 3 <= n <= 9\n"); return 1; }
    u32 N = 1u << M; int W = (N + 63) / 64;
    vector<pair<int,int>> pairs; vector<vector<int>> idx(n, vector<int>(n, -1));
    for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { idx[i][j] = idx[j][i] = (int)pairs.size(); pairs.push_back({i, j}); }
    vector<vector<int>> perms; { vector<int> p(n); iota(p.begin(), p.end(), 0); do perms.push_back(p); while (next_permutation(p.begin(), p.end())); }
    // constraint pattern of each permutation: key = (mask << 32) | ori; witness = up to three permutation indices
    unordered_map<u64, array<int,3>> unions; vector<u64> P; int zero_guards = 0;
    for (size_t q = 0; q < perms.size(); q++) { const vector<int>& p = perms[q]; u32 mask = 0, ori = 0; bool ok = true;
        for (int i = 0; i + 1 < n && ok; i++) { int a = p[i], b = p[i + 1]; int t = idx[a][b];
            if (t < 0) { if (b != a + 1) ok = false; } else { mask |= 1u << t; if (a < b) ori |= 1u << t; } }
        if (!ok) { zero_guards++; continue; }
        u64 key = ((u64)mask << 32) | ori;
        if (!unions.count(key)) { unions[key] = {(int)q, -1, -1}; P.push_back(key); } }
    printf("n = %d, M = %d free pairs, %u states: %zu permutations, %d zero guards, %zu distinct nonzero single guards\n", n, M, N, perms.size(), zero_guards, P.size());
    vector<u64> level = P;                       // unions first found at the previous level
    for (int j = 2; j <= k; j++) {
        vector<u64> next;
        for (u64 u : level) { array<int,3> wu = unions[u]; u32 mu = u >> 32, ou = (u32)u;
            for (u64 q : P) { u32 mq = q >> 32, oq = (u32)q; u32 c = mu & mq; if (((ou ^ oq) & c) != 0) continue;
                u64 z = ((u64)(mu | mq) << 32) | (ou | oq); if (unions.count(z)) continue;
                array<int,3> w = wu; w[j - 1] = unions[q][0]; unions[z] = w; next.push_back(z); } }
        printf("  distinct nonzero products of <= %d guards: %zu (%zu first appearing with %d factors)\n", j, unions.size(), next.size(), j); fflush(stdout);
        level.swap(next);
    }
    vector<u64> all; all.reserve(unions.size()); for (auto& kv : unions) all.push_back(kv.first);
    sort(all.begin(), all.end(), [](u64 a, u64 b) { int pa = __builtin_popcount((u32)(a >> 32)), pb = __builtin_popcount((u32)(b >> 32)); return pa != pb ? pa > pb : a < b; });
    vector<vector<u64>> basis; vector<int> pivot; vector<u64> bkey; u32 rank = 0; size_t done = 0;
    for (u64 u : all) {
        u32 mask = u >> 32, ori = (u32)u, fr = (N - 1) & ~mask; vector<u64> v(W, 0);
        for (u32 t = fr; ; t = (t - 1) & fr) { u32 s = ori | t; v[s >> 6] |= 1ULL << (s & 63); if (t == 0) break; }
        for (size_t i = 0; i < basis.size(); i++) { int pv = pivot[i]; if (v[pv >> 6] >> (pv & 63) & 1) { const u64* b = basis[i].data(); for (int w = 0; w < W; w++) v[w] ^= b[w]; } }
        int pv = -1; for (int w = 0; w < W && pv < 0; w++) if (v[w]) pv = w * 64 + __builtin_ctzll(v[w]);
        done++; if (done % 500000 == 0) { printf("  inserted %zu of %zu, GF(2) rank %u\n", done, all.size(), rank); fflush(stdout); }
        if (pv < 0) continue; basis.push_back(v); pivot.push_back(pv); bkey.push_back(u); rank++;
        if (rank == N) break;
    }
    printf("n = %d, k = %d: GF(2) rank of all products of <= %d named guards = %u of %u%s\n", n, k, k, rank, N, rank == N ? " (FULL: rational rank full as well; k restored tests span the presentation law)" : " (not full over GF(2))");
    if (cert) { FILE* f = fopen(cert, "w"); if (!f) { perror(cert); return 1; }
        fprintf(f, "# restored_products_exact_20261026 certificate: n = %d, k = %d, %u basis vectors; columns: mask ori q1 q2 q3\n", n, k, rank);
        fprintf(f, "# pair t = t-th pair (i, j) with j >= i+2 in lexicographic order, bit 1 = arc i -> j; q = 0-based lexicographic index of a permutation of 0..n-1; -1 = unused\n");
        for (size_t i = 0; i < bkey.size(); i++) { array<int,3> w = unions[bkey[i]]; fprintf(f, "%u %u %d %d %d\n", (u32)(bkey[i] >> 32), (u32)bkey[i], w[0], w[1], w[2]); }
        fclose(f); printf("certificate written to %s\n", cert); }
    return 0;
}
