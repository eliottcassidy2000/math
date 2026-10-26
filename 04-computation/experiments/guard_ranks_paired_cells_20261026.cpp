// guard_ranks_paired_cells_20261026.cpp -- opus-2026-10-26-S62
// On the marked-path carrier X_n (tournaments on 0..n-1 containing 0 -> 1 -> ... -> n-1; 2^M states, M = (n-1)(n-2)/2):
//  (1) named Hamilton-path guards g_p(x) = [p is a directed Hamiltonian path of x], p in S_n:
//      rank of the single guards modulo a large prime (= rational rank unless a modular accident; two primes reported)
//      and the GF(2) rank of the pairwise products g_p g_q ("test, then restore, then test"), a lower bound for the
//      rational rank that certifies full rank 2^M when it reaches it (for n = 7 the products are sampled);
//  (2) the anonymous reset protocol's paired cells ([x], [tau x]) with tau = complement of every non-path pair, [.] the
//      isomorphism class: the bipartite graph with the source classes and the target classes as vertices and one edge
//      per distinct pair; its edges E (distinguishable paired states), vertices V, components c, rank V - c
//      (independent terminal responses), and cycle rank E - V + c (lost joint-law directions).
// Build: g++ -O2 -std=c++17 -o grp guard_ranks_paired_cells_20261026.cpp ; run: ./grp [maxn=7] [pairSamples7=200000]
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <numeric>
using namespace std;
typedef unsigned long long u64; typedef long long i64;

static int n; static int adj[8][8]; static int order[8];
static vector<int> perm; static vector<vector<int>> perms;

static void genperms() { perms.clear(); vector<int> p(n); iota(p.begin(), p.end(), 0); do perms.push_back(p); while (next_permutation(p.begin(), p.end())); }
static bool isPath(const vector<int>& p) { for (int i = 0; i + 1 < n; i++) if (!adj[p[i]][p[i + 1]]) return false; return true; }
static u64 canon() {                                      // canonical code: min over relabelings of the adjacency bit string
    u64 best = ~0ULL; vector<int> p(n); iota(p.begin(), p.end(), 0);
    do { u64 code = 0; int b = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (adj[p[i]][p[j]]) code |= 1ULL << b; b++; }
         if (code < best) best = code; } while (next_permutation(p.begin(), p.end()));
    return best;
}
static int rank_mod(vector<vector<i64>> rows, i64 mod) {
    int r = 0, ncol = rows.empty() ? 0 : (int)rows[0].size();
    for (int c = 0; c < ncol && r < (int)rows.size(); c++) {
        int piv = -1; for (int i = r; i < (int)rows.size(); i++) if (rows[i][c] % mod) { piv = i; break; }
        if (piv < 0) continue; swap(rows[r], rows[piv]);
        i64 inv = 1, base = ((rows[r][c] % mod) + mod) % mod, e = mod - 2; while (e) { if (e & 1) inv = (__int128)inv * base % mod; base = (__int128)base * base % mod; e >>= 1; }
        for (auto& a : rows[r]) a = (__int128)(((a % mod) + mod) % mod) * inv % mod;
        for (int i = 0; i < (int)rows.size(); i++) if (i != r && rows[i][c] % mod) { i64 f = ((rows[i][c] % mod) + mod) % mod; for (int k = c; k < ncol; k++) rows[i][k] = ((rows[i][k] - (__int128)f * rows[r][k]) % mod + mod) % mod; }
        r++;
    }
    return r;
}
// GF(2) rank with bitsets, incremental with early stop at full rank
struct GF2Basis { int ncol, words; vector<vector<u64>> basis; vector<int> pivot; int rank = 0;
    GF2Basis(int c) : ncol(c), words((c + 63) / 64) {}
    bool add(vector<u64> v) { for (size_t i = 0; i < basis.size(); i++) { int p = pivot[i]; if (v[p >> 6] >> (p & 63) & 1) for (int w = 0; w < words; w++) v[w] ^= basis[i][w]; }
        int p = -1; for (int w = 0; w < words && p < 0; w++) if (v[w]) { p = w * 64 + __builtin_ctzll(v[w]); }
        if (p < 0) return false; basis.push_back(v); pivot.push_back(p); rank++; return true; } };

int main(int argc, char** argv) {
    int maxn = argc > 1 ? atoi(argv[1]) : 7; int samples7 = argc > 2 ? atoi(argv[2]) : 200000;
    u64 rng = 0x9E3779B97F4A7C15ULL; auto rnd = [&]() { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng * 0x2545F4914F6CDD1DULL; };
    for (n = 4; n <= maxn; n++) {
        vector<pair<int,int>> pairs; for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) pairs.push_back({i, j});
        int M = (int)pairs.size(); int N = 1 << M; genperms();
        // guard matrix: rows = permutations, columns = states
        vector<vector<i64>> G(perms.size(), vector<i64>(N, 0));
        vector<u64> srcClass(N), tgtClass(N);
        for (int code = 0; code < N; code++) {
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj[i][j] = 0;
            for (int i = 0; i + 1 < n; i++) adj[i][i + 1] = 1;
            for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; if (code >> k & 1) adj[i][j] = 1; else adj[j][i] = 1; }
            for (size_t q = 0; q < perms.size(); q++) G[q][code] = isPath(perms[q]) ? 1 : 0;
            srcClass[code] = canon();
            // tau: complement every non-path pair
            for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; adj[i][j] ^= 1; adj[j][i] ^= 1; }
            tgtClass[code] = canon();
        }
        // distinct single guard vectors
        set<vector<i64>> distinct(G.begin(), G.end()); vector<vector<i64>> singles(distinct.begin(), distinct.end());
        int r1a = rank_mod(singles, 1000000007LL), r1b = rank_mod(singles, 998244353LL);
        // pair products over GF(2)
        GF2Basis B(N); auto tovec = [&](const vector<i64>& v) { vector<u64> w((N + 63) / 64, 0); for (int c = 0; c < N; c++) if (v[c]) w[c >> 6] |= 1ULL << (c & 63); return w; };
        for (auto& s : singles) B.add(tovec(s));
        int rSinglesGF2 = B.rank; u64 tested = 0;
        if (n <= 6) {
            for (size_t a = 0; a < singles.size() && B.rank < N; a++) for (size_t b = a + 1; b < singles.size() && B.rank < N; b++) {
                vector<i64> v(N); bool nz = false; for (int c = 0; c < N; c++) { v[c] = singles[a][c] & singles[b][c]; nz |= v[c]; }
                if (nz) { B.add(tovec(v)); tested++; } }
        } else {
            for (int t = 0; t < samples7 && B.rank < N; t++) { size_t a = rnd() % singles.size(), b = rnd() % singles.size(); if (a == b) continue;
                vector<i64> v(N); bool nz = false; for (int c = 0; c < N; c++) { v[c] = singles[a][c] & singles[b][c]; nz |= v[c]; }
                if (nz) { B.add(tovec(v)); tested++; } }
        }
        // paired cells bipartite graph
        set<pair<u64,u64>> edges; set<u64> S, T;
        for (int code = 0; code < N; code++) { edges.insert({srcClass[code], tgtClass[code]}); S.insert(srcClass[code]); T.insert(tgtClass[code]); }
        // components by union-find over vertices (source classes and target classes kept distinct)
        map<u64,int> idS, idT; int V = 0; for (u64 s : S) idS[s] = V++; for (u64 t : T) idT[t] = V++;
        vector<int> par(V); iota(par.begin(), par.end(), 0);
        auto find = [&](int x) { while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; } return x; };
        for (auto& e : edges) { int a = find(idS[e.first]), b = find(idT[e.second]); if (a != b) par[a] = b; }
        set<int> roots; for (int v = 0; v < V; v++) roots.insert(find(v)); int c = (int)roots.size(); int E = (int)edges.size();
        printf("n = %d: states %d; single guards: %zu distinct of %zu permutations, rank %d / %d (mod two primes), GF(2) rank %d; "
               "pair products: %llu tested, singles+pairs GF(2) rank %d%s; paired cells: E = %d distinct ([x],[tau x]) pairs, source classes %zu, target classes %zu, V = %d, components %d, rank V - c = %d, cycle rank E - V + c = %d\n",
               n, N, singles.size(), perms.size(), r1a, r1b, rSinglesGF2, (unsigned long long)tested, B.rank, (B.rank == N ? " (full)" : ""), E, S.size(), T.size(), V, c, V - c, E - V + c);
        fflush(stdout);
    }
    return 0;
}
