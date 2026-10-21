// blue_line_deep_detours_20261020.cpp -- opus-2026-10-20-S47
//
// Shortest literal black detour of one closed blue line (T, P) by bidirectional breadth-first search in
// the labelled black K-move graph: vertices = labelled tournaments on [n], T ~ K_Q(T) for a black
// directed Hamiltonian path Q of T (K is an involution on marked pairs, and Q is black in T iff it is
// black in K_Q(T), so the relation is symmetric).  Reports the distance from T to K_P(T), or that it
// exceeds the depth bound.  Conventions as in blue_line_black_detours_20261020.cpp.
// Build: g++ -O2 -std=c++17 -o bldd blue_line_deep_detours_20261020.cpp ; run: ./bldd n tiling [maxdepth]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <unordered_set>
#include <algorithm>
using namespace std;
static int n, m; static int cellI[64], cellJ[64];
struct Tour { uint16_t adj[16]; };
static inline bool arc(const Tour &T, int u, int v) { return (T.adj[u] >> v) & 1; }
static inline bool isBlack(const Tour &T, const int *Q) { int r[16]; for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) { uint16_t out = T.adj[u]; while (out) { int v = __builtin_ctz(out); out &= out - 1; if (!arc(T, r[v], r[u])) return true; } } return false; }
static inline Tour kmove(const Tour &T, const int *Q) { Tour R; uint16_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) R.adj[u] = (~T.adj[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; R.adj[u] |= 1u << v; R.adj[v] &= ~(1u << u); } return R; }
static inline uint64_t enc(const Tour &T) { uint64_t c = 0; int b = 0; for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) { if (arc(T, u, v)) c |= 1ull << b; b++; } return c; }
static inline Tour dec(uint64_t c) { Tour T; for (int u = 0; u < n; u++) T.adj[u] = 0; int b = 0; for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) { if ((c >> b) & 1) T.adj[u] |= 1u << v; else T.adj[v] |= 1u << u; b++; } return T; }
template <class F> static void hamPaths(const Tour &T, F &&f) { int Q[16];
    for (int s = 0; s < n; s++) { int depth = 0; Q[0] = s; uint16_t used = 1u << s; uint16_t cand[16]; cand[0] = T.adj[s] & ~used;
        while (depth >= 0) { if (depth == n - 1) { f(Q); depth--; if (depth >= 0) used &= ~(1u << Q[depth + 1]); continue; }
            if (!cand[depth]) { depth--; if (depth >= 0) used &= ~(1u << Q[depth + 1]); continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1; Q[depth + 1] = v; used |= 1u << v; depth++; cand[depth] = T.adj[v] & ~used; } } }
static vector<uint64_t> blackNeighbours(uint64_t c) { Tour T = dec(c); vector<uint64_t> out; hamPaths(T, [&](const int *Q) { if (isBlack(T, Q)) out.push_back(enc(kmove(T, Q))); }); sort(out.begin(), out.end()); out.erase(unique(out.begin(), out.end()), out.end()); return out; }
int main(int argc, char **argv) {
    n = atoi(argv[1]); uint64_t t = strtoull(argv[2], 0, 10); if (n > 12) { fprintf(stderr, "n <= 12 (tiling in 64 bits)\n"); return 1; } int maxd = argc > 3 ? atoi(argv[3]) : 6;
    m = 0; for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { cellI[m] = i; cellJ[m] = j; m++; }
    Tour T; for (int u = 0; u < n; u++) T.adj[u] = 0; for (int i = 0; i + 1 < n; i++) T.adj[i] |= 1u << (i + 1);
    for (int k = 0; k < m; k++) { int i = cellI[k], j = cellJ[k]; if ((t >> k) & 1) T.adj[i] |= 1u << j; else T.adj[j] |= 1u << i; }
    int P[16]; for (int i = 0; i < n; i++) P[i] = i; Tour KP = kmove(T, P);
    uint64_t s = enc(T), g = enc(KP);
    printf("n=%d tiling=%llu scores:", n, (unsigned long long)t); for (int u = 0; u < n; u++) printf(" %d", __builtin_popcount(T.adj[u])); printf("; partner scores:"); for (int u = 0; u < n; u++) printf(" %d", __builtin_popcount(KP.adj[u])); printf("\n");
    // bidirectional BFS: layers from s and from g
    vector<unordered_set<uint64_t>> A(1), B(1); A[0].insert(s); B[0].insert(g);
    unordered_set<uint64_t> seenA = {s}, seenB = {g};
    for (int d = 0; d <= maxd; d++) {
        // check meeting: distance = da + db
        for (int da = 0; da < (int)A.size(); da++) for (int db = 0; db < (int)B.size(); db++) if (da + db == d) {
            const auto &X = A[da].size() < B[db].size() ? A[da] : B[db]; const auto &Y = A[da].size() < B[db].size() ? B[db] : A[da];
            for (uint64_t x : X) if (Y.count(x)) { printf("black distance %d (layers %d + %d; frontier sizes", d, da, db); for (auto &L : A) printf(" %zu", L.size()); printf(" |"); for (auto &L : B) printf(" %zu", L.size()); printf(")\n"); return 0; } }
        // expand the smaller frontier
        bool expandA = A.back().size() <= B.back().size();
        auto &F = expandA ? A : B; auto &seen = expandA ? seenA : seenB;
        unordered_set<uint64_t> next;
        for (uint64_t c : F.back()) for (uint64_t y : blackNeighbours(c)) if (!seen.count(y)) { seen.insert(y); next.insert(y); }
        F.push_back(move(next));
        if (F.back().empty()) break;
    }
    printf("black distance > %d (frontier sizes", maxd); for (auto &L : A) printf(" %zu", L.size()); printf(" |"); for (auto &L : B) printf(" %zu", L.size()); printf(")\n");
    return 0;
}
