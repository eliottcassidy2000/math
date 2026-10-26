// carousel_line_partners_20261026.cpp -- opus-2026-10-26-S59
//
// The complement-line neighbourhood of the order-12 carousel fibre of cyclictriple_composition_20261026:
// T_w = Q_3[A_{w_0}, A_{w_1}, A_{w_2}] with A_0 = A+ = (C3 -> 1), A_1 = A- = (1 -> C3), Q_3 = C3 on the blocks.
// For every Hamiltonian path P of T_w the line partner is K_P(T_w) = T_w^op xor E(P) (the arcs of P kept, every
// other pair reversed; the free-arc involution tau of the packets).  We record, for every presentation:
//   * whether the partner is prime (no module of size 2..n-1; pair-closure test), else its smallest module size;
//   * that none of the six proper modules of T_w (three blocks, three triangles) is a module of the partner
//     (the module-shattering theorem, checked literally);
//   * whether the presentation is blue (positional reversal an anti-automorphism) or black;
//   * the partner's isomorphism class (colour-refinement canonical form with a 128-bit word) -> number of
//     distinct partner classes, their multiplicities, and whether any partner class is itself in the fibre.
// Build: g++ -O2 -std=c++17 -o clp carousel_line_partners_20261026.cpp ; run: ./clp [m=3]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <chrono>
using namespace std;
typedef unsigned __int128 u128;
static int n;
struct Tour { uint16_t adj[16]; };
static inline bool arc(const Tour &T, int u, int v) { return (T.adj[u] >> v) & 1; }

// ---- canonical form (refinement + individualisation; word = adjacency bits in canonical order) ----
struct Canon { u128 best; int aut; };
struct Canoniser {
    uint16_t adj[16]; int nn; u128 best; int autc;
    void refine(int col[16]) {
        for (int iter = 0; iter < 16; iter++) {
            uint64_t key[16];
            for (int v = 0; v < nn; v++) { int c[16] = {0}; for (int w = 0; w < nn; w++) if ((adj[v] >> w) & 1) c[col[w]]++;
                uint64_t k = (uint64_t)col[v]; for (int x = 0; x < nn; x++) k = k * 17 + c[x]; key[v] = k; }
            uint64_t sorted[16]; for (int v = 0; v < nn; v++) sorted[v] = key[v]; sort(sorted, sorted + nn); int nd = (int)(unique(sorted, sorted + nn) - sorted);
            int s2[16]; for (int v = 0; v < nn; v++) s2[v] = col[v]; sort(s2, s2 + nn); int oldnd = (int)(unique(s2, s2 + nn) - s2);
            int newcol[16]; for (int v = 0; v < nn; v++) newcol[v] = (int)(lower_bound(sorted, sorted + nd, key[v]) - sorted);
            bool same = true; for (int v = 0; v < nn; v++) if (newcol[v] != col[v]) { same = false; break; }
            for (int v = 0; v < nn; v++) col[v] = newcol[v];
            if (same || nd == oldnd) { if (nd == oldnd && !same) continue; break; }
            if (nd == nn) break;
        }
    }
    void leaf(const int col[16]) { int perm[16]; for (int v = 0; v < nn; v++) perm[col[v]] = v; u128 w = 0;
        for (int a = 0; a < nn; a++) for (int b = a + 1; b < nn; b++) w = (w << 1) | (u128)((adj[perm[a]] >> perm[b]) & 1u);
        if (w < best) { best = w; autc = 1; } else if (w == best) autc++; }
    void search(int col[16]) { refine(col); int s2[16]; for (int v = 0; v < nn; v++) s2[v] = col[v]; sort(s2, s2 + nn); int nd = (int)(unique(s2, s2 + nn) - s2);
        if (nd == nn) { leaf(col); return; } int cnt[16] = {0}; for (int v = 0; v < nn; v++) cnt[col[v]]++; int cell = 0; while (cnt[cell] < 2) cell++;
        for (int v = 0; v < nn; v++) if (col[v] == cell) { int c2[16]; for (int u = 0; u < nn; u++) c2[u] = col[u] + (col[u] > cell ? 1 : 0); c2[v] = cell; for (int u = 0; u < nn; u++) if (u != v && col[u] == cell) c2[u] = cell + 1; search(c2); } }
    Canon run(const Tour &T) { nn = n; for (int v = 0; v < nn; v++) adj[v] = T.adj[v]; best = ~(u128)0; autc = 0; int col[16], sc[16];
        for (int v = 0; v < nn; v++) sc[v] = __builtin_popcount(adj[v]); int s2[16]; for (int v = 0; v < nn; v++) s2[v] = sc[v]; sort(s2, s2 + nn); int nd = (int)(unique(s2, s2 + nn) - s2);
        for (int v = 0; v < nn; v++) col[v] = (int)(lower_bound(s2, s2 + nd, sc[v]) - s2); search(col); return {best, autc}; }
};
struct H128 { size_t operator()(const u128 &x) const { uint64_t lo = (uint64_t)x, hi = (uint64_t)(x >> 64); return std::hash<uint64_t>()(lo * 0x9E3779B97F4A7C15ull ^ (hi + 0x7F4A7C15ull)); } };

// ---- modules ----
static inline bool isModule(const Tour &T, uint16_t M) {
    uint16_t full = (1u << n) - 1, out = full & ~M; if (M == 0 || out == 0) return true;
    int m0 = __builtin_ctz(M);
    for (int x = 0; x < n; x++) if ((out >> x) & 1) { bool o = arc(T, x, m0); uint16_t MM = M & ~(1u << m0);
        while (MM) { int y = __builtin_ctz(MM); MM &= MM - 1; if (arc(T, x, y) != o) return false; } }
    return true;
}
// smallest module containing u and v: closure by adding non-uniform outside vertices
static inline uint16_t closureModule(const Tour &T, int u, int v) {
    uint16_t M = (1u << u) | (1u << v), full = (1u << n) - 1;
    bool changed = true;
    while (changed) { changed = false; uint16_t out = full & ~M;
        for (int x = 0; x < n; x++) if ((out >> x) & 1) { int m0 = __builtin_ctz(M); bool o = arc(T, x, m0); uint16_t MM = M & ~(1u << m0); bool bad = false;
            while (MM) { int y = __builtin_ctz(MM); MM &= MM - 1; if (arc(T, x, y) != o) { bad = true; break; } }
            if (bad) { M |= 1u << x; changed = true; } } }
    return M;
}
static inline int smallestModuleSize(const Tour &T) { int best = n + 1; uint16_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) { uint16_t M = closureModule(T, u, v); if (M != full) best = min(best, __builtin_popcount(M)); }
    return best; }   // n+1 means prime

static inline bool isBlue(const Tour &T, const int *P) { int r[16]; for (int i = 0; i < n; i++) r[P[i]] = P[n - 1 - i];
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (u != v && arc(T, u, v) && !arc(T, r[v], r[u])) return false; return true; }
static inline Tour kmove(const Tour &T, const int *Q) { Tour R; uint16_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) R.adj[u] = (~T.adj[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; R.adj[u] |= 1u << v; R.adj[v] &= ~(1u << u); } return R; }

template <class F> static void hamPaths(const Tour &T, F &&f) { int Q[16];
    for (int s = 0; s < n; s++) { int depth = 0; Q[0] = s; uint16_t used = 1u << s; uint16_t cand[16]; cand[0] = T.adj[s] & ~used;
        while (depth >= 0) { if (depth == n - 1) { f(Q); depth--; if (depth >= 0) used &= ~(1u << Q[depth + 1]); continue; }
            if (!cand[depth]) { depth--; if (depth >= 0) used &= ~(1u << Q[depth + 1]); continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1; Q[depth + 1] = v; used |= 1u << v; depth++; cand[depth] = T.adj[v] & ~used; } } }

int main(int argc, char **argv) {
    int m = argc > 1 ? atoi(argv[1]) : 3; n = 4 * m; if (n > 16) { fprintf(stderr, "n <= 16\n"); return 1; }
    auto t0 = chrono::steady_clock::now();
    // canonical forms of the fibre members themselves
    vector<Tour> fibre; vector<u128> fibreCanon; Canoniser C;
    for (int w = 0; w < (1 << m); w++) {
        Tour T; for (int u = 0; u < n; u++) T.adj[u] = 0;
        for (int b = 0; b < m; b++) { int o = 4 * b; T.adj[o] |= 1u << (o + 1); T.adj[o + 1] |= 1u << (o + 2); T.adj[o + 2] |= 1u << o;   // C3
            bool minus = (w >> b) & 1; for (int k = 0; k < 3; k++) { if (minus) T.adj[o + 3] |= 1u << (o + k); else T.adj[o + k] |= 1u << (o + 3); } }
        for (int b = 0; b < m; b++) for (int c = 0; c < m; c++) if (b != c) { int d = ((c - b) % m + m) % m; if (d >= 1 && d <= (m - 1) / 2) for (int k = 0; k < 4; k++) for (int l = 0; l < 4; l++) T.adj[4 * b + k] |= 1u << (4 * c + l); }
        fibre.push_back(T); fibreCanon.push_back(C.run(T).best);
    }
    // distinct classes of the fibre (necklaces)
    vector<int> rep; for (int w = 0; w < (1 << m); w++) { bool seen = false; for (int r : rep) if (fibreCanon[r] == fibreCanon[w]) { seen = true; break; } if (!seen) rep.push_back(w); }
    printf("# carousel fibre m = %d (n = %d): %zu words, %zu classes:", m, n, fibre.size(), rep.size());
    for (int r : rep) { Canon c = C.run(fibre[r]); printf(" w=%d(|Aut|=%d)", r, c.aut); } printf("\n");
    unordered_map<u128, long long, H128> partnerClasses; partnerClasses.reserve(1 << 22);
    unordered_map<u128, int, H128> partnerSmallestModule;
    for (int r : rep) {
        const Tour &T = fibre[r]; long long H = 0, blue = 0, prime = 0, nonprime = 0, oldModuleSurvived = 0, inFibre = 0; map<int, long long> smallMod;
        // the six proper modules of T_w: blocks and triangles
        vector<uint16_t> oldMods; for (int b = 0; b < m; b++) { oldMods.push_back(0xFu << (4 * b)); oldMods.push_back(0x7u << (4 * b)); }
        for (auto M : oldMods) if (!isModule(T, M)) { fprintf(stderr, "internal: block/triangle not a module of T\n"); return 2; }
        unordered_map<u128, long long, H128> local; local.reserve(1 << 20);
        hamPaths(T, [&](const int *P) {
            H++; if (isBlue(T, P)) blue++;
            Tour K = kmove(T, P);
            for (auto M : oldMods) if (isModule(K, M)) oldModuleSurvived++;
            int sm = smallestModuleSize(K); if (sm > n) prime++; else { nonprime++; smallMod[sm]++; }
            u128 cw = C.run(K).best; local[cw]++; partnerClasses[cw]++; partnerSmallestModule[cw] = sm;
            for (size_t i = 0; i < fibreCanon.size(); i++) if (fibreCanon[i] == cw) { inFibre++; break; }
        });
        printf("w=%d: H=%lld blue=%lld black=%lld; partners prime=%lld nonprime=%lld (smallest-module sizes:", r, H, blue, H - blue, prime, nonprime);
        for (auto &kv : smallMod) printf(" %d:%lld", kv.first, kv.second); printf("); old modules surviving=%lld; partners in the fibre=%lld; distinct partner classes=%zu\n", oldModuleSurvived, inFibre, local.size()); fflush(stdout);
    }
    // multiplicity histogram over all partner classes
    map<long long, long long> hist; long long primeClasses = 0; for (auto &kv : partnerClasses) { hist[kv.second]++; if (partnerSmallestModule[kv.first] > n) primeClasses++; }
    printf("all representatives: distinct partner classes %zu (prime %lld); multiplicity histogram (presentations per class: classes):", partnerClasses.size(), primeClasses);
    int shown = 0; for (auto &kv : hist) { if (shown++ < 12) printf(" %lld:%lld", kv.first, kv.second); } printf("%s\n", hist.size() > 12 ? " ..." : "");
    printf("done in %.0f s\n", chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    return 0;
}
