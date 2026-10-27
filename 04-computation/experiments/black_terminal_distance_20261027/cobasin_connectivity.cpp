// cobasin_connectivity.cpp -- mac-mini-2026-10-27-potentials
// A sufficient test for black-core connectivity (HYP-9348 (b)) built on THM-4811.
// Converse-merged classes (the metagraph's nodes).  Bipartite graph B_n: lambda = 1 strong classes S joined to every
// non-strong class reached from S by a BLACK descent word (start in a minimum out-cut, end outside; since lambda = 1 the
// target is non-strong).  If B_n is connected and covers every non-strong class that has a black line, then the black
// core is connected at order n: every non-pure-blue strong class black-descends to lambda = 1 (THM-4811; the E-family
// takes one extra black move inside lambda = 1), and every black line at a non-strong class ends at a lambda = 1 class.
// Also counts non-strong classes with no black line at all (pure blue) for comparison.
// Build: g++ -O3 -march=native -std=c++17 -I/opt/homebrew/opt/nauty/include -o cobasin cobasin_connectivity.cpp /opt/homebrew/opt/nauty/lib/libnauty1.a
// Run:   gentourng -q n | ./cobasin n
#define MAXN WORDSIZE
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <unordered_map>
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
static Tour converse(const Tour &T) { Tour R; memset(&R, 0, sizeof R); for (int u = 0; u < n; u++) for (uint32_t o = T.out[u]; o; o &= o - 1) R.out[__builtin_ctz(o)] |= 1u << u; return R; }
static uint64_t mergedCode(const Tour &T) { uint64_t a = canonCode(T), b = canonCode(converse(T)); return a < b ? a : b; }
static bool strongT(const Tour &T) {
    uint32_t full = (1u << n) - 1, seen = 1, fr = 1;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= T.out[__builtin_ctz(f)]; nx &= ~seen; seen |= nx; fr = nx; }
    if (seen != full) return false;
    Tour R = converse(T); seen = 1; fr = 1;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= R.out[__builtin_ctz(f)]; nx &= ~seen; seen |= nx; fr = nx; }
    return seen == full;
}
static int lambdaCuts(const Tour &T, std::vector<uint32_t> &cuts) {
    int full = (1 << n) - 1, best = 1 << 20; cuts.clear();
    for (int X = 1; X < full; X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T.out[x] & ~X & full); if (f < best) { best = f; cuts.clear(); } if (f == best) cuts.push_back(X); }
    return best;
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
static std::vector<int> par;
static int findp(int x) { while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; } return x; }
int main(int argc, char **argv) {
    n = atoi(argv[1]);
    std::unordered_map<uint64_t, int> id; std::vector<char> isStrong, hasBlackLine;
    auto getId = [&](uint64_t c, bool strong) { auto it = id.find(c); if (it != id.end()) return it->second; int k = (int)par.size(); id[c] = k; par.push_back(k); isStrong.push_back(strong); hasBlackLine.push_back(0); return k; };
    char line[256]; long long nsClasses = 0, lam1 = 0, lam1WithBlack = 0, edges = 0;
    std::vector<Tour> nonstrong;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.out[i] |= 1u << j; else T.out[j] |= 1u << i; k++; }
        if (!strongT(T)) { nonstrong.push_back(T); continue; }
        std::vector<uint32_t> cuts; int lam = lambdaCuts(T, cuts);
        if (lam != 1) continue;
        lam1++;
        bool sep[16][16]; memset(sep, 0, sizeof sep);
        for (uint32_t X : cuts) for (int a = 0; a < n; a++) if (X >> a & 1) for (int b = 0; b < n; b++) if (!(X >> b & 1)) sep[a][b] = true;
        int sId = -1; bool any = false;
        paths(T, [&](const int *Q) {
            if (!sep[Q[0]][Q[n - 1]] || !isBlack(T, Q)) return;
            Tour R = kmove(T, Q);
            uint64_t c = mergedCode(R);
            if (sId < 0) sId = getId(mergedCode(T), true);
            int t = getId(c, false); hasBlackLine[t] = 1; any = true; edges++;
            int a = findp(sId), b = findp(t); if (a != b) par[a] = b;
        });
        if (any) lam1WithBlack++;
    }
    // non-strong merged classes and which have a black line at all (to any class)
    long long nsMerged = 0, nsBlack = 0, nsCovered = 0; std::unordered_map<uint64_t, char> seenNS;
    for (auto &T : nonstrong) {
        uint64_t c = mergedCode(T); if (seenNS.count(c)) continue; seenNS[c] = 1; nsMerged++;
        bool black = false;
        paths(T, [&](const int *Q) { if (!black && isBlack(T, Q)) black = true; });
        if (black) { nsBlack++; if (id.count(c)) nsCovered++; }
    }
    // components among covered classes
    std::unordered_map<int, int> comps; for (auto &p : id) comps[findp(p.second)]++;
    printf("n=%d: non-strong merged classes %lld, with a black line %lld, covered by B_n %lld; lambda=1 strong classes %lld, with a black descent %lld; black descent lines %lld; components of B_n: %zu\n",
           n, nsMerged, nsBlack, nsCovered, lam1, lam1WithBlack, edges, comps.size());
    return 0;
}
