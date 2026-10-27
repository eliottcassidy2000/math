// forward_closure.cpp -- mac-mini-2026-10-27-slab
// The good slab classes near the transitive end: B_k = classes at level L = k reachable from TT by R-descending native
// moves (paths from a vertex above the mean score to one below it).  By THM-4819 a slab class at level k is on a
// geodesic between the transitive and the regular class (and is exactly k moves from TT) iff it lies in B_k.
// Computes B_0, ..., B_K by BFS over canonical forms and writes them (one canonical key per line, level-tagged) so that a
// top-down search from a regular tournament can stop as soon as it lands in B_K (meet in the middle).
// Build: g++ -O3 -march=native -std=c++17 -I/opt/homebrew/opt/nauty/include -o forward_closure forward_closure.cpp /opt/homebrew/opt/nauty/lib/libnauty1.a
// Run:   ./forward_closure n K outfile
#define MAXN WORDSIZE
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <unordered_set>
#include <algorithm>
static int n, m;
struct Tour { uint32_t out[32]; };
struct Key { uint64_t w[4]; bool operator==(const Key &o) const { return !memcmp(w, o.w, sizeof w); } };
struct KH { size_t operator()(const Key &k) const { return (k.w[0] * 0x9E3779B97F4A7C15ULL) ^ (k.w[1] * 0xC2B2AE3D27D4EB4FULL) ^ (k.w[2] * 0x165667B19E3779F9ULL) ^ (k.w[3] >> 7); } };
static Key canon(const Tour &T, Tour *canonT) {
    graph g[MAXN], cg[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    static DEFAULTOPTIONS_DIGRAPH(options); statsblk stats; options.getcanon = TRUE;
    EMPTYGRAPH(g, 1, n);
    for (int u = 0; u < n; u++) for (uint32_t o = T.out[u]; o; o &= o - 1) ADDONEARC(g, u, __builtin_ctz(o), 1);
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, cg);
    Key k; memset(&k, 0, sizeof k); int pos = 0;
    if (canonT) memset(canonT, 0, sizeof(Tour));
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        bool a = ISELEMENT(GRAPHROW(cg, i, 1), j);
        if (a) k.w[pos >> 6] |= 1ULL << (pos & 63);
        if (canonT) { if (a) canonT->out[i] |= 1u << j; else canonT->out[j] |= 1u << i; }
        pos++;
    }
    return k;
}
static Tour kmove(const Tour &T, const int *Q) {
    Tour R; uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    for (int u = 0; u < n; u++) R.out[u] = (~T.out[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; R.out[u] |= 1u << v; R.out[v] &= ~(1u << u); }
    return R;
}
template <class F> static void paths_from(const Tour &T, int s, uint32_t endmask, F &&f) {
    int Q[32]; uint32_t cand[32]; int depth = 0; Q[0] = s; uint32_t used = 1u << s; cand[0] = T.out[s] & ~used;
    while (depth >= 0) {
        if (depth == n - 1) { if (endmask >> Q[depth] & 1) f(Q); used &= ~(1u << Q[depth]); depth--; continue; }
        if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
        int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
        depth++; Q[depth] = v; used |= 1u << v; cand[depth] = T.out[v] & ~used;
    }
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); m = (n - 1) / 2; int K = atoi(argv[2]); FILE *fo = argc > 3 ? fopen(argv[3], "w") : nullptr;
    Tour TT; memset(&TT, 0, sizeof TT); for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) TT.out[i] |= 1u << j;
    Tour c0; Key k0 = canon(TT, &c0);
    std::vector<Tour> level{c0}; std::unordered_set<Key, KH> seen{k0};
    if (fo) { fprintf(fo, "0"); for (int q = 0; q < 4; q++) fprintf(fo, " %016llx", (unsigned long long)k0.w[q]); fprintf(fo, "\n"); }
    printf("n=%d level 0: 1 class\n", n); fflush(stdout);
    for (int L = 1; L <= K; L++) {
        std::vector<Tour> nxt; std::unordered_set<Key, KH> lev; long moves = 0;
        for (const Tour &T : level) {
            int s[32]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(T.out[v]);
            uint32_t low = 0; for (int v = 0; v < n; v++) if (s[v] < m) low |= 1u << v;
            for (int a = 0; a < n; a++) {
                if (s[a] <= m) continue;
                paths_from(T, a, low, [&](const int *Q) {
                    moves++; Tour R = kmove(T, Q); Tour cR; Key k = canon(R, &cR);
                    if (lev.insert(k).second) nxt.push_back(cR);
                });
            }
        }
        if (fo) for (const Key &k : lev) { fprintf(fo, "%d", L); for (int q = 0; q < 4; q++) fprintf(fo, " %016llx", (unsigned long long)k.w[q]); fprintf(fo, "\n"); }
        printf("n=%d level %d: %zu classes (from %ld R-descending moves)\n", n, L, nxt.size(), moves); fflush(stdout);
        level.swap(nxt);
        if (level.size() > 20000000) { printf("stopping: level too large\n"); break; }
    }
    if (fo) fclose(fo);
    return 0;
}
