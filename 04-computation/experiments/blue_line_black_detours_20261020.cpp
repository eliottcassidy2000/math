// blue_line_black_detours_20261020.cpp -- opus-2026-10-20-S47
//
// Literal black detours for closed blue lines, exhaustively over all closed blue presentations at a
// given order n (every sigma-fixed tiling with the apex tile backward), grouped by the endpoint layer
// k = d_T(z) (score pattern (n-1-k, k) of the endpoint_three_repairs_20261020h packet).
//
// Conventions (endpoint_two/three packets): P = (0,1,...,n-1) is a Hamiltonian path of T; K_Q(T) =
// T^op xor E(Q) keeps the arcs of Q and reverses every other pair; a word Q is blue when its positional
// reversal is an anti-automorphism of the tournament it lives in, black otherwise.  A detour of length L
// is a sequence of black words Q_1,...,Q_L, Q_i a directed Hamiltonian path of T_(i-1), T_i = K_(Q_i)(T_(i-1)),
// with T_L = K_P(T) as a LABELLED tournament (equivalently E(Q_1) xor ... xor E(Q_L) = E(P) for odd L).
// One move is impossible (Q_1 = P is blue).  We search L = 2 (needs E(P) xor E(Q_1) to be a black
// directed Hamiltonian path of T_1) and L = 3 (all black Q_1, all black Q_2 of T_1; Q_3 is then forced:
// E_3 = E(P) xor E_1 xor E_2 must be a black directed Hamiltonian path of T_2), stopping at the first
// certificate of each line and reporting the lines with none (exhaustive for L <= 3).  For an even number of
// moves the path-edge sets XOR to the COMPLEMENT of E(P), so L = 2 needs C(n-1,2) <= 2(n-1), i.e. n <= 6.
//
// Build: g++ -O2 -std=c++17 -o bld blue_line_black_detours_20261020.cpp ; run: ./bld n [k | all] [codeLo codeHi]
// (the optional code range slices the 2^e sigma-orbit codes for parallel runs; n <= 11 since pair masks use 64 bits;
//  a trailing "sub" restricts to the closed-middle sub-case with y, c in I, the frontier of endpoint_four_repairs_20261021)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <chrono>
#include <algorithm>
using namespace std;

static int n, m;
static int cellI[64], cellJ[64], cellIdx[16][16], sigmaCell[64];
static uint64_t pairBit[16][16];   // bit index of the unordered pair {u,v} among all C(n,2) pairs
static uint64_t EP;                // pair mask of the base path
static uint64_t ALLP;              // mask of all C(n,2) pairs

struct Tour { uint16_t adj[16]; };

static inline bool arc(const Tour &T, int u, int v) { return (T.adj[u] >> v) & 1; }
static inline bool isBlack(const Tour &T, const int *Q) {
    int r[16]; for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) { uint16_t out = T.adj[u]; while (out) { int v = __builtin_ctz(out); out &= out - 1; if (!arc(T, r[v], r[u])) return true; } }
    return false;
}
static inline uint64_t pathMask(const int *Q) { uint64_t e = 0; for (int i = 0; i + 1 < n; i++) e |= pairBit[Q[i]][Q[i + 1]]; return e; }
static inline Tour kmove(const Tour &T, const int *Q) {
    Tour R; uint16_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) R.adj[u] = (~T.adj[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; R.adj[u] |= 1u << v; R.adj[v] &= ~(1u << u); }   // Q is directed u->v in T
    return R;
}
// Is the pair mask E the edge set of a directed Hamiltonian path of T?  If so write it into Q.
static bool directedPath(uint64_t E, const Tour &T, int *Q) {
    if (__builtin_popcountll(E) != n - 1) return false;
    int deg[16] = {0}; int nb[16][2];
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if ((E >> __builtin_ctzll(pairBit[u][v])) & 1) {
        if (deg[u] == 2 || deg[v] == 2) return false; nb[u][deg[u]++] = v; nb[v][deg[v]++] = u; }
    int ends[2], ne = 0; for (int u = 0; u < n; u++) { if (deg[u] == 0) return false; if (deg[u] == 1) { if (ne == 2) return false; ends[ne++] = u; } }
    if (ne != 2) return false;
    for (int s = 0; s < 2; s++) {
        int prev = -1, cur = ends[s]; bool ok = true; Q[0] = cur;
        for (int i = 1; i < n; i++) { int nxt = (nb[cur][0] != prev) ? nb[cur][0] : nb[cur][1]; if (deg[cur] == 1 && i > 1) { ok = false; break; } if (!arc(T, cur, nxt)) { ok = false; break; } prev = cur; cur = nxt; Q[i] = cur; }
        if (ok) { // connectivity: all n vertices visited
            uint16_t seen = 0; for (int i = 0; i < n; i++) seen |= 1u << Q[i]; if (seen == (uint16_t)((1u << n) - 1)) return true; return false; }
    }
    return false;
}

// Enumerate directed Hamiltonian paths of T by DFS; callback returns true to stop.
template <class F> static bool hamPaths(const Tour &T, F &&f) {
    int Q[16];
    struct Frame { int v; uint16_t rest; };
    for (int s = 0; s < n; s++) {
        // iterative DFS
        int depth = 0; Q[0] = s; uint16_t used = 1u << s; uint16_t cand[16]; cand[0] = T.adj[s] & ~used;
        while (depth >= 0) {
            if (depth == n - 1) { if (f(Q)) return true; depth--; if (depth >= 0) used &= ~(1u << Q[depth + 1]); continue; }
            if (!cand[depth]) { depth--; if (depth >= 0) used &= ~(1u << Q[depth + 1]); continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
            Q[depth + 1] = v; used |= 1u << v; depth++; cand[depth] = T.adj[v] & ~used;
        }
    }
    return false;
}

int main(int argc, char **argv) {
    n = argc > 1 ? atoi(argv[1]) : 8; int layer = -1; if (argc > 2 && strcmp(argv[2], "all") != 0) layer = atoi(argv[2]);
    if (n < 4 || n > 11) { fprintf(stderr, "n in 4..11 (pair masks in 64 bits)\n"); return 1; }
    unsigned long long codeLo = 0, codeHi = 0; if (argc > 4) { codeLo = strtoull(argv[3], 0, 10); codeHi = strtoull(argv[4], 0, 10); }
    bool subcase = argc > 5 && strcmp(argv[5], "sub") == 0;   // closed middle (c -> b) with y, c in I: the sub-case left open by endpoint_four_repairs_20261021
    m = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cellIdx[i][j] = -1;
    for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { cellI[m] = i; cellJ[m] = j; cellIdx[i][j] = m; m++; }
    for (int k = 0; k < m; k++) sigmaCell[k] = cellIdx[n - 1 - cellJ[k]][n - 1 - cellI[k]];
    { int b = 0; for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) { pairBit[u][v] = pairBit[v][u] = 1ull << b; b++; } }
    EP = 0; for (int i = 0; i + 1 < n; i++) EP |= pairBit[i][i + 1];
    ALLP = 0; for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) ALLP |= pairBit[u][v];
    // sigma-orbit representatives of tiles (fixed tiles and the lower member of each 2-orbit)
    vector<int> reps; for (int k = 0; k < m; k++) if (sigmaCell[k] >= k) reps.push_back(k);
    int apex = cellIdx[0][n - 1]; int e = (int)reps.size();
    printf("# black detours of closed blue lines, n = %d, m = %d, e = %d, closed blue presentations = %d, layer = %s\n", n, m, e, 1 << (e - 1), layer < 0 ? "all" : argv[2]); fflush(stdout);
    auto t0 = chrono::steady_clock::now();
    const int K = 16; long long lines[K] = {0}, len2[K] = {0}, len3[K] = {0}, none[K] = {0}, pbinv[K] = {0}, genuine[K] = {0}; double maxSec[K] = {0}; long long maxH[K] = {0};
    long long done = 0;
    if (codeHi == 0) codeHi = 1ull << e;
    printf("# code range [%llu, %llu)%s\n", codeLo, codeHi, subcase ? ", sub-case: closed middle with y, c in I" : ""); fflush(stdout);
    for (unsigned long long code = codeLo; code < codeHi; code++) {
        uint64_t t = 0; for (int i = 0; i < e; i++) if ((code >> i) & 1) { t |= 1ull << reps[i]; t |= 1ull << sigmaCell[reps[i]]; }
        if ((t >> apex) & 1) continue;   // closed: apex tile backward (z -> a)
        Tour T; for (int u = 0; u < n; u++) T.adj[u] = 0;
        for (int i = 0; i + 1 < n; i++) T.adj[i] |= 1u << (i + 1);
        for (int k = 0; k < m; k++) { int i = cellI[k], j = cellJ[k]; if ((t >> k) & 1) T.adj[i] |= 1u << j; else T.adj[j] |= 1u << i; }
        int kk = __builtin_popcount(T.adj[n - 1]);
        if (layer >= 0 && kk != layer) continue;
        if (subcase && !(arc(T, n - 2, 1) && arc(T, n - 3, 0) && arc(T, n - 2, 0))) continue;
        int P[16]; for (int i = 0; i < n; i++) P[i] = i;
        if (isBlack(T, P)) { fprintf(stderr, "internal: base word not blue\n"); return 2; }
        lines[kk]++;
        auto c0 = chrono::steady_clock::now();
        // pure-blue flags of the two classes: every Hamiltonian path blue
        Tour KP = kmove(T, P);
        bool pbT = !hamPaths(T, [&](const int *Q) { return isBlack(T, Q); });
        bool pbK = !hamPaths(KP, [&](const int *Q) { return isBlack(KP, Q); });
        if (pbT || pbK) { pbinv[kk]++; none[kk]++; continue; }
        // --- length 2 ---
        bool found2 = false; long long H = 0;
        hamPaths(T, [&](const int *Q1) {
            H++; if (!isBlack(T, Q1)) return false;
            // two moves: E1 xor E2 must be the complement of E(P) (even number of K moves), possible only for n <= 6
            Tour T1 = kmove(T, Q1); uint64_t E2 = ALLP ^ EP ^ pathMask(Q1); int Q2[16];
            if (directedPath(E2, T1, Q2) && isBlack(T1, Q2)) {
                Tour T2 = kmove(T1, Q2); for (int u = 0; u < n; u++) if (T2.adj[u] != KP.adj[u]) { fprintf(stderr, "internal: 2-move target mismatch\n"); exit(4); }
                found2 = true; return true; }
            return false; });
        if (H > maxH[kk]) maxH[kk] = H;
        if (found2) { len2[kk]++; }
        else {
            bool found3 = false;
            hamPaths(T, [&](const int *Q1) {
                if (!isBlack(T, Q1)) return false;
                Tour T1 = kmove(T, Q1); uint64_t E1 = pathMask(Q1);
                return hamPaths(T1, [&](const int *Q2) {
                    if (!isBlack(T1, Q2)) return false;
                    uint64_t E3 = EP ^ E1 ^ pathMask(Q2); if (__builtin_popcountll(E3) != n - 1) return false;
                    Tour T2 = kmove(T1, Q2); int Q3[16];
                    if (directedPath(E3, T2, Q3) && isBlack(T2, Q3)) {
                        // verify the literal target
                        Tour T3 = kmove(T2, Q3), KP = kmove(T, P);
                        for (int u = 0; u < n; u++) if (T3.adj[u] != KP.adj[u]) { fprintf(stderr, "internal: target mismatch\n"); exit(3); }
                        found3 = true; return true; }
                    return false; });
            });
            if (found3) len3[kk]++; else { none[kk]++; genuine[kk]++; printf("NO 3-detour (both classes non-pure-blue): n=%d layer=%d tiling=%llu scores:", n, kk, (unsigned long long)t); for (int u = 0; u < n; u++) printf(" %d", __builtin_popcount(T.adj[u])); printf("\n"); fflush(stdout); }
        }
        double cs = chrono::duration<double>(chrono::steady_clock::now() - c0).count(); if (cs > maxSec[kk]) maxSec[kk] = cs;
        done++; if ((done & 0x3FFF) == 0) { fprintf(stderr, "  %lld lines done, %.0f s\n", done, chrono::duration<double>(chrono::steady_clock::now() - t0).count()); }
    }
    printf("layer k: closed blue lines / shortest detour 2 / shortest detour 3 / none within 3 = pure-blue partner + genuine   (max H(T), max seconds per line)\n");
    for (int k = 0; k < K; k++) if (lines[k]) printf("   k=%d: %lld / %lld / %lld / %lld = %lld + %lld   (%lld, %.3f)\n", k, lines[k], len2[k], len3[k], none[k], pbinv[k], genuine[k], maxH[k], maxSec[k]);
    printf("done in %.0f s\n", chrono::duration<double>(chrono::steady_clock::now() - t0).count());
    return 0;
}
