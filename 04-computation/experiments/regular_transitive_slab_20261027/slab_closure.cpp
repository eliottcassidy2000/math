// slab_closure.cpp -- mac-mini-2026-10-27-slab
// The transitive-regular slab at odd order n: classes with L + R = (n^2-1)/8 (sorted scores between i and m).
// A slab class lies on a geodesic between the transitive and the regular class iff it is reachable from TT by
// R-descending native moves (THM-4819), and then d(T, TT) = L(T).  This program
//   (1) reads every class from gentourng, keeps the slab classes (canonical keys), counts them per level L;
//   (2) computes the closure of TT under R-descending moves (paths from a vertex above the mean score to one below),
//       level by level, inside the slab;
//   (3) prints per level: slab classes, reached (= on a TT-regular geodesic, = L-exact), unreached, and for the
//       unreached ones: with a source or sink / strong with arc-connectivity 1 / other.
// Build: g++ -O3 -march=native -std=c++17 -I/opt/homebrew/opt/nauty/include -o slab_closure slab_closure.cpp /opt/homebrew/opt/nauty/lib/libnauty1.a
// Run:   gentourng -q n | ./slab_closure n
#define MAXN WORDSIZE
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <algorithm>
static int n, m;
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
static Tour kmove(const Tour &T, const int *Q) {
    Tour R; uint32_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) R.out[u] = (~T.out[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; R.out[u] |= 1u << v; R.out[v] &= ~(1u << u); }
    return R;
}
template <class F> static void paths_from(const Tour &T, int s, F &&f) {
    int Q[16]; uint32_t cand[16]; int depth = 0; Q[0] = s; uint32_t used = 1u << s; cand[0] = T.out[s] & ~used;
    while (depth >= 0) {
        if (depth == n - 1) { f(Q); used &= ~(1u << Q[depth]); depth--; continue; }
        if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
        int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
        depth++; Q[depth] = v; used |= 1u << v; cand[depth] = T.out[v] & ~used;
    }
}
static void scores(const Tour &T, int *s) { for (int v = 0; v < n; v++) s[v] = __builtin_popcount(T.out[v]); }
static int Lof(const int *s0) { int s[16]; memcpy(s, s0, n * sizeof(int)); std::sort(s, s + n); int l = 0; for (int i = 0; i < n; i++) l += abs(s[i] - i); return l / 2; }
static int Rof(const int *s) { int r = 0; for (int v = 0; v < n; v++) r += abs(s[v] - m); return r / 2; }
static int arcConn(const Tour &T) {
    uint32_t full = (1u << n) - 1; int best = 99;
    for (uint32_t X = 1; X < full; X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T.out[x] & ~X & full); if (f < best) best = f; }
    return best;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); m = (n - 1) / 2; int c = (n * n - 1) / 8;
    std::vector<Tour> reps; std::vector<int> lev; std::unordered_map<uint64_t, int> id;
    char line[256]; long all = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        all++;
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.out[i] |= 1u << j; else T.out[j] |= 1u << i; k++; }
        int s[16]; scores(T, s); int L = Lof(s), R = Rof(s);
        if (L + R != c) continue;
        id[canonCode(T)] = (int)reps.size(); reps.push_back(T); lev.push_back(L);
    }
    int C = (int)reps.size();
    std::vector<char> reached(C, 0);
    std::vector<int> frontier;
    for (int i = 0; i < C; i++) if (lev[i] == 0) { reached[i] = 1; frontier.push_back(i); }
    long moves = 0;
    for (int L = 0; L < c; L++) {
        std::vector<int> nxt;
        for (int u : frontier) {
            const Tour &T = reps[u]; int s[16]; scores(T, s);
            for (int a = 0; a < n; a++) {
                if (s[a] <= m) continue;    // R-descending: start above the mean ...
                paths_from(T, a, [&](const int *Q) {
                    if (s[Q[n - 1]] >= m) return;   // ... and end below it
                    moves++;
                    Tour Rr = kmove(T, Q);
                    auto it = id.find(canonCode(Rr));
                    if (it == id.end()) { fprintf(stderr, "image outside the slab\n"); exit(2); }
                    int t = it->second;
                    if (lev[t] != L + 1) { fprintf(stderr, "level error\n"); exit(3); }
                    if (!reached[t]) { reached[t] = 1; nxt.push_back(t); }
                });
            }
        }
        frontier.swap(nxt);
    }
    printf("n=%d: classes %ld, slab classes %d, R-descending moves explored %ld\n", n, all, C, moves);
    printf("  level L: slab / reached (on a TT-regular geodesic, L-exact) / unreached (source or sink, strong lambda = 1, other)\n");
    for (int L = 0; L <= c; L++) {
        int tot = 0, re = 0, ss = 0, l1 = 0, ot = 0;
        for (int i = 0; i < C; i++) if (lev[i] == L) {
            tot++; if (reached[i]) { re++; continue; }
            int s[16]; scores(reps[i], s); bool srcsink = false; for (int v = 0; v < n; v++) if (s[v] == 0 || s[v] == n - 1) srcsink = true;
            if (srcsink) ss++; else if (arcConn(reps[i]) == 1) l1++; else ot++;
        }
        printf("   L=%2d R=%2d: %7d / %7d / %7d  (%d, %d, %d)\n", L, c - L, tot, re, tot - re, ss, l1, ot);
    }
    return 0;
}
