// singleton_distance.cpp -- mac-mini-2026-10-27-potentials
// For the strong classes where lambda < delta0 = min_v min(s_v, n-1-s_v) (first possible at order 10), decide whether
// the native distance to the singleton locus SING (some source or sink) equals delta0, by a depth-delta0 search in which
// every move must raise max_v |d_v| by one (the only moves compatible with distance delta0, by the singleton potential).
// Also reports whether the final state can be a singleton plus a STRONG core (CORE), and whether all moves can be black.
// Depth-first search with memoised failures (hashed labelled states).
// Input: one class per line (gentourng upper-triangle code).  Build: g++ -O3 -std=c++17 -o singleton_distance singleton_distance.cpp
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <set>
#include <algorithm>
#include <unordered_set>
static int n;
struct Tour { uint32_t out[16]; bool operator<(const Tour &o) const { return memcmp(out, o.out, sizeof out) < 0; } };
static int score(const Tour &T, int v) { return __builtin_popcount(T.out[v]); }
static int maxdev2(const Tour &T) { int m = 0; for (int v = 0; v < n; v++) m = std::max(m, std::abs(2 * score(T, v) - (n - 1))); return m; }
static bool strongMask(const Tour &T, uint32_t vs) {
    if (!vs) return true; int st = __builtin_ctz(vs); uint32_t seen = 1u << st, fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= T.out[__builtin_ctz(f)] & vs; nx &= ~seen; seen |= nx; fr = nx; }
    if (seen != vs) return false;
    seen = 1u << st; fr = seen;
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
// enumerate Hamiltonian paths with a prescribed start (s >= 0) or end (e >= 0)
template <class F> static void pathsStart(const Tour &T, int s, F &&f) {
    int Q[16]; uint32_t cand[16]; int depth = 0; Q[0] = s; uint32_t used = 1u << s; cand[0] = T.out[s] & ~used;
    while (depth >= 0) {
        if (depth == n - 1) { f(Q); used &= ~(1u << Q[depth]); depth--; continue; }
        if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
        int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
        depth++; Q[depth] = v; used |= 1u << v; cand[depth] = T.out[v] & ~used;
    }
}
// raising moves: some v with |d_v| maximal and (d_v <= 0, v = start) or (d_v >= 0, v = end).  Ends handled on the converse.
static void raisingMoves(const Tour &T, std::vector<std::pair<Tour,bool>> &res) {
    int M = maxdev2(T); res.clear();
    Tour Rv; memset(&Rv, 0, sizeof Rv); for (int u = 0; u < n; u++) for (uint32_t o = T.out[u]; o; o &= o - 1) Rv.out[__builtin_ctz(o)] |= 1u << u;
    for (int v = 0; v < n; v++) {
        int d2 = 2 * score(T, v) - (n - 1); if (std::abs(d2) != M) continue;
        if (d2 <= 0) pathsStart(T, v, [&](const int *Q) { res.push_back({kmove(T, Q), isBlack(T, Q)}); });
        if (d2 >= 0) pathsStart(Rv, v, [&](const int *Qr) { int Q[16]; for (int i = 0; i < n; i++) Q[i] = Qr[n - 1 - i]; res.push_back({kmove(T, Q), isBlack(T, Q)}); });
    }
}
struct H64 { size_t operator()(const std::pair<uint64_t,int> &p) const { return p.first * 1000003ull ^ (uint64_t)p.second; } };
static uint64_t hashT(const Tour &T) { uint64_t h = 1469598103934665603ull; for (int u = 0; u < n; u++) { h ^= T.out[u]; h *= 1099511628211ull; } return h; }
static std::unordered_set<std::pair<uint64_t,int>, H64> failed[2];
static bool isTarget(const Tour &S, bool core) {
    uint32_t full = (1u << n) - 1;
    for (int v = 0; v < n; v++) { int s = score(S, v); if (s == 0 || s == n - 1) { if (!core || strongMask(S, full & ~(1u << v))) return true; } }
    return false;
}
static bool dfs(const Tour &S, int rem, bool black, bool core) {
    if (rem == 0) return isTarget(S, core);
    auto key = std::make_pair(hashT(S), rem * 2 + (core ? 1 : 0));
    if (failed[black].count(key)) return false;
    std::vector<std::pair<Tour,bool>> mv; raisingMoves(S, mv);
    for (auto &m : mv) { if (black && !m.second) continue; if (dfs(m.first, rem - 1, black, core)) return true; }
    failed[black].insert(key);
    return false;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); char line[256]; long long cls = 0, okAny = 0, okCore = 0, okBlack = 0, okBlackCore = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.out[i] |= 1u << j; else T.out[j] |= 1u << i; k++; }
        cls++;
        int d0 = (n - 1 - maxdev2(T)) / 2;
        failed[0].clear(); failed[1].clear();
        bool a = dfs(T, d0, false, false), ac = dfs(T, d0, false, true), b = dfs(T, d0, true, false), bc = dfs(T, d0, true, true);
        okAny += a; okCore += ac; okBlack += b; okBlackCore += bc;
        if (!a || !ac || !bc) printf("class %lld delta0 %d: SING %d CORE %d (black %d %d) %s", cls, d0, a, ac, b, bc, line);
    }
    printf("n=%d classes %lld: distance delta0 to SING %lld, to CORE %lld; with black moves: SING %lld, CORE %lld\n", n, cls, okAny, okCore, okBlack, okBlackCore);
    return 0;
}
