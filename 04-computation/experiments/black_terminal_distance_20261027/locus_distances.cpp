// locus_distances.cpp -- mac-mini-2026-10-27-potentials
// Exact native-move distances on isomorphism classes of order n (moves of any colour, and black moves only) from every
// class to three target loci: NS (non-strong), SING (some vertex is a source or a sink), CORE (a source or sink whose
// deletion is strong: "a singleton plus a strong core").  Compared with the potentials lambda (arc-strong connectivity)
// and delta0 = min_v min(s_v, n-1-s_v).  The move graph on classes is undirected (moves are involutions and blackness is
// preserved), so a multi-source BFS from the targets gives the distances.
// Build: g++ -O3 -march=native -std=c++17 -I/opt/homebrew/opt/nauty/include -o locus_distances locus_distances.cpp /opt/homebrew/opt/nauty/lib/libnauty1.a
// Run:   gentourng -q n | ./locus_distances n
#define MAXN WORDSIZE
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <map>
#include <deque>
#include <algorithm>
static int n;
struct Tour { uint32_t out[16]; };
static uint64_t canonCode(const Tour &T) {
    graph g[MAXN], cg[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    static DEFAULTOPTIONS_DIGRAPH(options); statsblk stats; options.getcanon = TRUE;
    EMPTYGRAPH(g, 1, n);
    for (int u = 0; u < n; u++) for (uint32_t o = T.out[u]; o; o &= o - 1) ADDONEARC(g, u, __builtin_ctz(o), 1);
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, cg);
    uint64_t code = 0; int b = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (ISELEMENT(GRAPHROW(cg, i, 1), j)) code |= 1ull << b; b++; }
    return code;
}
static bool strongMask(const Tour &T, uint32_t vs) {
    if (!vs) return true; int start = __builtin_ctz(vs);
    uint32_t seen = 1u << start, fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= T.out[__builtin_ctz(f)] & vs; nx &= ~seen; seen |= nx; fr = nx; }
    if (seen != vs) return false;
    seen = 1u << start; fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) { int v = __builtin_ctz(f); for (uint32_t w = vs; w; w &= w - 1) { int u = __builtin_ctz(w); if (T.out[u] >> v & 1) nx |= 1u << u; } } nx &= ~seen; seen |= nx; fr = nx; }
    return seen == vs;
}
static inline bool isBlack(const Tour &T, const int *Q) {
    int r[16]; for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) for (uint32_t o = T.out[u]; o; o &= o - 1) { int v = __builtin_ctz(o); if (!((T.out[r[v]] >> r[u]) & 1)) return true; }
    return false;
}
static Tour kmove(const Tour &T, const int *Q) {
    Tour R; uint32_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) R.out[u] = (~T.out[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; R.out[u] |= 1u << v; R.out[v] &= ~(1u << u); }
    return R;
}
template <class F> static void paths(const Tour &T, F &&f) {
    for (int s = 0; s < n; s++) {
        int Q[16]; uint32_t cand[16]; int depth = 0; Q[0] = s; uint32_t used = 1u << s; cand[0] = T.out[s] & ~used;
        while (depth >= 0) {
            if (depth == n - 1) { f(Q); used &= ~(1u << Q[depth]); depth--; continue; }
            if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
            depth++; Q[depth] = v; used |= 1u << v; cand[depth] = T.out[v] & ~used;
        }
    }
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); uint32_t full = (1u << n) - 1;
    std::vector<Tour> reps; std::unordered_map<uint64_t, int> id;
    char line[256];
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.out[i] |= 1u << j; else T.out[j] |= 1u << i; k++; }
        id[canonCode(T)] = (int)reps.size(); reps.push_back(T);
    }
    int C = (int)reps.size();
    std::vector<std::vector<int>> adjAny(C), adjBlack(C);
    std::vector<int> lam(C), d0(C); std::vector<char> tNS(C), tSING(C), tCORE(C);
    for (int c = 0; c < C; c++) {
        const Tour &T = reps[c];
        int best = 99; for (uint32_t X = 1; X < full; X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T.out[x] & ~X & full); best = std::min(best, f); }
        lam[c] = best; int dd = 99; bool sing = false, core = false;
        for (int v = 0; v < n; v++) { int s = __builtin_popcount(T.out[v]); dd = std::min(dd, std::min(s, n - 1 - s)); if (s == 0 || s == n - 1) { sing = true; if (strongMask(T, full & ~(1u << v))) core = true; } }
        d0[c] = dd; tNS[c] = best == 0; tSING[c] = sing; tCORE[c] = core;
        paths(T, [&](const int *Q) {
            Tour R = kmove(T, Q); int t = id.at(canonCode(R));
            adjAny[c].push_back(t); if (isBlack(T, Q)) adjBlack[c].push_back(t);
        });
        std::sort(adjAny[c].begin(), adjAny[c].end()); adjAny[c].erase(std::unique(adjAny[c].begin(), adjAny[c].end()), adjAny[c].end());
        std::sort(adjBlack[c].begin(), adjBlack[c].end()); adjBlack[c].erase(std::unique(adjBlack[c].begin(), adjBlack[c].end()), adjBlack[c].end());
    }
    auto bfs = [&](const std::vector<char> &tgt, const std::vector<std::vector<int>> &adj) {
        std::vector<int> d(C, -1); std::deque<int> q; for (int c = 0; c < C; c++) if (tgt[c]) { d[c] = 0; q.push_back(c); }
        while (!q.empty()) { int c = q.front(); q.pop_front(); for (int t : adj[c]) if (d[t] < 0) { d[t] = d[c] + 1; q.push_back(t); } }
        return d;
    };
    const char *names[3] = {"NS", "SING", "CORE"}; std::vector<char> *tg[3] = {&tNS, &tSING, &tCORE};
    printf("n=%d classes=%d\n", n, C);
    for (int k = 0; k < 3; k++) {
        for (int col = 0; col < 2; col++) {
            auto d = bfs(*tg[k], col ? adjBlack : adjAny);
            const std::vector<int> &pot = (k == 0) ? lam : d0;
            std::map<std::pair<int,int>, int> hist; int unreachable = 0, below = 0;
            for (int c = 0; c < C; c++) { if (d[c] < 0) { unreachable++; continue; } hist[{pot[c], d[c] - pot[c]}]++; if (d[c] < pot[c]) below++; }
            printf("  target %-4s %-5s: potential %s; unreachable %d; d < potential %d; (potential, d - potential): ", names[k], col ? "black" : "any", k == 0 ? "lambda" : "delta0", unreachable, below);
            for (auto &p : hist) printf("(%d,+%d):%d ", p.first.first, p.first.second, p.second);
            printf("\n");
        }
    }
    return 0;
}
