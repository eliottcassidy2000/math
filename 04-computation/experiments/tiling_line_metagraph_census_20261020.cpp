// tiling_line_metagraph_census_20261020.cpp -- opus-2026-10-20-S46
//
// Exhaustive census of the owner's tiling line metagraph for 3 <= n <= 8.
//
// Objects (THM-643 frame, base path 0 -> 1 -> ... -> n-1):
//   tiles = non-path pairs {i < j}, j >= i+2, m = C(n-1,2); bit 1 = forward arc i -> j.
//   kappa = flip all tiles (the owner's LINE pairing {t, kappa t}).
//   sigma = staircase reflection (i,j) -> (n-1-j, n-1-i) (converse along the reversed path).
//   blue tiling = sigma-fixed; a line is blue iff its tilings are (kappa commutes with sigma).
//   depth(t) = least k >= 0 such that t restricted to the concentric sub-staircase
//              Delta_k = {(i,j): i >= k, j <= n-1-k} is sigma-symmetric (depth 0 = blue);
//              kappa preserves depth, so every line has a depth.
// For every tiling we compute the isomorphism class of its tournament (canonical form by
// score-block pruned brute force, exact |Aut|), the class of sigma t, and the class of kappa t.
// Output: per-class data (H, |Aut|, tc = H/|Aut|, bm = blue multiplicity, SC flag, depth
// profile), the line multigraph by colour / depth / node type, self-lines (ordinary and
// converse-merged), sigma-twin and Klein-orbit checks, the handshake parity constraint, and
// the controls of THM-643 / THM-793 / THM-849.
//
// Build: g++ -O2 -std=c++17 -o tlmc tiling_line_metagraph_census_20261020.cpp
// Run:   ./tlmc n            (n = 3..8; n = 8 enumerates 2^21 tilings, ~1 min)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <array>
#include <numeric>
#include <cstring>
using namespace std;

static int n, m;
static int cellI[64], cellJ[64];      // cell k -> (i,j)
static int cellIdx[16][16];           // (i,j) -> k, -1 if not a cell
static int sigmaCell[64];             // sigma permutation of cells
static uint64_t deltaMask[16];        // cells of Delta_k

struct ClassInfo { uint64_t canon; int aut; long long H; int tc; int bm; int sigmaClass; int depthProfile[8]; int score[16]; uint32_t rep; };

static inline void buildAdj(uint32_t t, int adj[16][16]) {
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) adj[u][v] = 0;
    for (int i = 0; i + 1 < n; i++) adj[i][i + 1] = 1;
    for (int k = 0; k < m; k++) {
        int i = cellI[k], j = cellJ[k];
        if ((t >> k) & 1) adj[i][j] = 1; else adj[j][i] = 1;
    }
}

// canonical form: min over score-respecting labellings of the upper-triangular bit word
static uint64_t canonForm(const int adj[16][16], int &autCount, int scoreOut[16]) {
    int sc[16];
    for (int v = 0; v < n; v++) { sc[v] = 0; for (int w = 0; w < n; w++) sc[v] += adj[v][w]; }
    int order[16]; iota(order, order + n, 0);
    sort(order, order + n, [&](int a, int b) { return sc[a] < sc[b]; });
    for (int v = 0; v < n; v++) scoreOut[v] = sc[order[v]];
    // groups of equal score
    vector<pair<int,int>> groups; // [start, end)
    for (int s = 0; s < n;) { int e = s; while (e < n && sc[order[e]] == sc[order[s]]) e++; groups.push_back({s, e}); s = e; }
    int perm[16]; for (int v = 0; v < n; v++) perm[v] = order[v];
    uint64_t best = ~0ULL; autCount = 0;
    // enumerate permutations within each group via next_permutation on segments (odometer)
    for (auto &g : groups) sort(perm + g.first, perm + g.second);
    while (true) {
        uint64_t w = 0;
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) w = (w << 1) | (uint64_t)adj[perm[a]][perm[b]];
        if (w < best) { best = w; autCount = 1; } else if (w == best) autCount++;
        // odometer increment
        int gi = (int)groups.size() - 1;
        while (gi >= 0) {
            if (next_permutation(perm + groups[gi].first, perm + groups[gi].second)) break;
            gi--;
        }
        if (gi < 0) break;
    }
    return best;
}

static long long hamPaths(const int adj[16][16]) {
    int N = 1 << n; vector<long long> dp((size_t)N * n, 0);
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
    for (int mask = 1; mask < N; mask++) for (int v = 0; v < n; v++) {
        long long c = dp[(size_t)mask * n + v]; if (!c) continue;
        for (int w = 0; w < n; w++) if (!((mask >> w) & 1) && adj[v][w]) dp[(size_t)(mask | (1 << w)) * n + w] += c;
    }
    long long tot = 0; for (int v = 0; v < n; v++) tot += dp[(size_t)(N - 1) * n + v];
    return tot;
}

static inline uint32_t sigmaOf(uint32_t t) { uint32_t s = 0; for (int k = 0; k < m; k++) if ((t >> k) & 1) s |= 1u << sigmaCell[k]; return s; }
static inline int depthOf(uint32_t t) { uint32_t s = sigmaOf(t); for (int k = 0; ; k++) { if (((t ^ s) & deltaMask[k]) == 0) return k; } }

int main(int argc, char **argv) {
    n = argc > 1 ? atoi(argv[1]) : 6;
    if (n < 3 || n > 8) { fprintf(stderr, "n in 3..8\n"); return 1; }
    m = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cellIdx[i][j] = -1;
    for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { cellI[m] = i; cellJ[m] = j; cellIdx[i][j] = m; m++; }
    for (int k = 0; k < m; k++) sigmaCell[k] = cellIdx[n - 1 - cellJ[k]][n - 1 - cellI[k]];
    int maxDepth = 0;
    for (int k = 0; k < 16; k++) { deltaMask[k] = 0; for (int c = 0; c < m; c++) if (cellI[c] >= k && cellJ[c] <= n - 1 - k) deltaMask[k] |= 1ull << c; if (deltaMask[k] == 0 && maxDepth == 0) maxDepth = k; }
    uint32_t NT = 1u << m, ALL = NT - 1;
    printf("# tiling line metagraph census, n = %d, m = %d, tilings = %u, max depth = %d\n", n, m, NT, maxDepth);

    vector<int> cls(NT, -1); vector<uint8_t> depth(NT, 0);
    vector<ClassInfo> classes; unordered_map<uint64_t, int> byCanon; byCanon.reserve(20000);
    int adj[16][16];
    for (uint32_t t = 0; t < NT; t++) {
        buildAdj(t, adj);
        int aut, score[16]; uint64_t c = canonForm(adj, aut, score);
        auto it = byCanon.find(c); int id;
        if (it == byCanon.end()) {
            id = (int)classes.size(); byCanon[c] = id;
            ClassInfo ci; ci.canon = c; ci.aut = aut; ci.H = hamPaths(adj); ci.tc = 0; ci.bm = 0; ci.sigmaClass = -1; ci.rep = t;
            for (int d = 0; d < 8; d++) ci.depthProfile[d] = 0; for (int v = 0; v < n; v++) ci.score[v] = score[v];
            classes.push_back(ci);
        } else id = it->second;
        cls[t] = id; classes[id].tc++;
        int d = depthOf(t); depth[t] = (uint8_t)d; classes[id].depthProfile[d]++;
        if (d == 0) classes[id].bm++;
    }
    int C = (int)classes.size();
    // sigma class of each class (class of sigma(rep))
    int SC = 0;
    for (int id = 0; id < C; id++) { classes[id].sigmaClass = cls[sigmaOf(classes[id].rep)]; if (classes[id].sigmaClass == id) SC++; }
    // fibre checks
    long long sumTc = 0; int badFibre = 0, badParity = 0, badBlue = 0;
    for (auto &ci : classes) {
        sumTc += ci.tc;
        if ((long long)ci.tc * ci.aut != ci.H) badFibre++;
        if (ci.tc % 2 == 0) badParity++;
        bool sc = ci.sigmaClass >= 0 && classes[ci.sigmaClass].canon == ci.canon;
        if (sc && (ci.bm % 2 == 0)) badBlue++;
        if (!sc && ci.bm != 0) badBlue++;
    }
    int merged = (C + SC) / 2;
    printf("classes %d (SC %d, merged %d); sum tc = %lld (2^m = %u); fibre law violations %d; even tc %d; blue-parity violations %d\n", C, SC, merged, sumTc, NT, badFibre, badParity, badBlue);
    // node types
    int PB = 0, MX = 0, PK = 0;
    auto nodeType = [&](int id) { auto &ci = classes[id]; bool sc = ci.sigmaClass == id; if (!sc) return 2; return ci.bm == ci.tc ? 0 : 1; };
    for (int id = 0; id < C; id++) { int ty = nodeType(id); if (ty == 0) PB++; else if (ty == 1) MX++; else PK++; }
    printf("node types: pure-blue %d, mixed %d, pure-black %d\n", PB, MX, PK);
    // handshake parity constraint on the fundamental black graph B (sigma-twins identified):
    int oddSC = 0; for (int id = 0; id < C; id++) if (classes[id].sigmaClass == id && (((classes[id].tc - classes[id].bm) / 2) % 2 == 1)) oddSC++;
    int nsMerged = (C - SC) / 2;
    printf("handshake: #NS merged nodes = %d, #SC with (tc-bm)/2 odd = %d, sum parity = %d (must be 0)\n", nsMerged, oddSC, (nsMerged + oddSC) % 2);
    // lines
    long long blueLines = 0, blackLines = 0, selfB = 0, selfK = 0, mselfB = 0, mselfK = 0, pathRevSelf = 0;
    long long linesByDepth[8] = {0}; long long byType[3][3][2]; memset(byType, 0, sizeof byType);
    map<pair<int,int>, pair<int,int>> edgeMult; // (c1<=c2) -> (blue, black) multiplicities (class level)
    long long twinsOK = 0, twinsBad = 0;
    for (uint32_t t = 0; t < NT; t++) {
        uint32_t u = t ^ ALL; if (u < t) continue;
        int c1 = cls[t], c2 = cls[u]; bool blue = depth[t] == 0;
        if (blue) blueLines++; else blackLines++;
        linesByDepth[depth[t]]++;
        int a = min(c1, c2), b = max(c1, c2);
        auto &e = edgeMult[{a, b}]; if (blue) e.first++; else e.second++;
        int ta = nodeType(a), tb = nodeType(b); if (ta > tb) swap(ta, tb); byType[ta][tb][blue ? 0 : 1]++;
        if (c1 == c2) { if (blue) selfB++; else selfK++; }
        int s1 = classes[c1].sigmaClass, s2 = classes[c2].sigmaClass;
        bool mergedSelf = (c1 == c2) || (c2 == s1);
        if (mergedSelf) { if (blue) mselfB++; else mselfK++; }
        if (c2 == s1 && c1 != c2) pathRevSelf++;
        if (!blue) { // sigma twin: line {sigma t, sigma u} joins classes (s1, s2): same merged pair
            uint32_t st = sigmaOf(t); int d1 = cls[st], d2 = cls[st ^ ALL];
            bool same = (min(d1, d2) == min(s1, s2) && max(d1, d2) == max(s1, s2));
            if (same) twinsOK++; else twinsBad++;
        }
    }
    printf("lines: blue %lld (2^(e-1) check), black %lld, total %lld = 2^(m-1)\n", blueLines, blackLines, blueLines + blackLines);
    printf("self-lines (class level): blue %lld, black %lld; converse-merged self-loops: blue %lld, black %lld; path-reversal self-lines (kappa t in the converse class, c1 != c2) %lld\n", selfB, selfK, mselfB, mselfK, pathRevSelf);
    printf("sigma-twin check on black lines (twin joins the sigma-image classes): ok %lld, bad %lld\n", twinsOK, twinsBad);
    printf("lines by depth:"); for (int d = 0; d <= maxDepth; d++) printf(" d%d=%lld", d, linesByDepth[d]); printf("\n");
    const char *tn[3] = {"PB", "MX", "PK"};
    printf("lines by node types (blue, black):"); for (int a = 0; a < 3; a++) for (int b = a; b < 3; b++) printf(" (%s,%s)=(%lld,%lld)", tn[a], tn[b], byType[a][b][0], byType[a][b][1]); printf("\n");
    // multiplicity structure of the class-level line multigraph
    map<int,int> blueMultHist, blackMultHist; int distinctPairs = 0;
    for (auto &kv : edgeMult) { distinctPairs++; if (kv.second.first) blueMultHist[kv.second.first]++; if (kv.second.second) blackMultHist[kv.second.second]++; }
    printf("distinct class pairs joined: %d; blue multiplicity histogram:", distinctPairs); for (auto &kv : blueMultHist) printf(" %d:%d", kv.first, kv.second);
    printf("; black multiplicity histogram:"); for (auto &kv : blackMultHist) printf(" %d:%d", kv.first, kv.second); printf("\n");
    // merged-level black multiplicities (must all be even by the sigma-twin theorem)
    map<pair<int,int>, int> mergedBlack; auto mid = [&](int c) { return min(c, classes[c].sigmaClass); };
    for (auto &kv : edgeMult) { int a = mid(kv.first.first), b = mid(kv.first.second); if (a > b) swap(a, b); mergedBlack[{a, b}] += kv.second.second; }
    int oddMergedBlack = 0; for (auto &kv : mergedBlack) if (kv.second % 2) oddMergedBlack++;
    printf("merged-level black multiplicities: %zu node pairs, odd multiplicities %d (must be 0)\n", mergedBlack.size(), oddMergedBlack);
    // depth profile summary: number of tilings per depth, and per-class depth profiles aggregated by node type
    long long tilesByDepth[8] = {0}; for (uint32_t t = 0; t < NT; t++) tilesByDepth[depth[t]]++;
    printf("tilings by depth:"); for (int d = 0; d <= maxDepth; d++) printf(" d%d=%lld", d, tilesByDepth[d]); printf("\n");
    // Klein orbits on black self-line endpoints: Q^K divisible by 4 (THM-849)
    long long QK = 2 * selfK, QB = 2 * selfB; printf("Q^B = %lld, Q^K = %lld (Q^K mod 4 = %lld)\n", QB, QK, QK % 4);
    // per-class dump
    char fname[128]; snprintf(fname, sizeof fname, "tiling_line_metagraph_census_20261020_classes_n%d.tsv", n);
    FILE *f = fopen(fname, "w");
    fprintf(f, "class\tcanon\tH\tAut\ttc\tbm\tSC\ttype\tsigmaClass\tscores\tdepthProfile\n");
    for (int id = 0; id < C; id++) { auto &ci = classes[id]; fprintf(f, "%d\t%llu\t%lld\t%d\t%d\t%d\t%d\t%s\t%d\t", id, (unsigned long long)ci.canon, ci.H, ci.aut, ci.tc, ci.bm, ci.sigmaClass == id ? 1 : 0, tn[nodeType(id)], ci.sigmaClass);
        for (int v = 0; v < n; v++) fprintf(f, "%d%s", ci.score[v], v + 1 < n ? "," : "\t"); for (int d = 0; d <= maxDepth; d++) fprintf(f, "%d%s", ci.depthProfile[d], d < maxDepth ? "," : "\n"); }
    fclose(f);
    snprintf(fname, sizeof fname, "tiling_line_metagraph_census_20261020_lines_n%d.tsv", n);
    f = fopen(fname, "w"); fprintf(f, "classA\tclassB\tblue\tblack\n");
    for (auto &kv : edgeMult) fprintf(f, "%d\t%d\t%d\t%d\n", kv.first.first, kv.first.second, kv.second.first, kv.second.second);
    fclose(f);
    printf("wrote per-class and line files for n = %d\n", n);
    return 0;
}
