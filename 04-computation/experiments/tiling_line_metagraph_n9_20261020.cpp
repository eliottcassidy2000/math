// tiling_line_metagraph_n9_20261020.cpp -- opus-2026-10-20-S47 (S48: H-descent columns minHblack / maxHblack and sink statistics)
//
// Line metagraph census for n <= 9 with a colour-refinement / individualisation canonical form
// (exact |Aut| from the number of leaves attaining the canonical word).  Frame as in
// tiling_line_metagraph_census_20261020.cpp (THM-4731): base path 0 -> ... -> n-1, tiles = non-path
// pairs, bit 1 = forward arc; kappa = flip all tiles (the line pairing); sigma = staircase reflection
// (converse along the reversed path); blue = sigma-fixed.
//
// Computes per class: H, |Aut|, tc = H/|Aut|, bm (blue multiplicity), self-converse flag, node type
// (pure-blue / mixed / pure-black); the line multigraph's union-find components for all lines, for
// blue lines (on self-converse classes) and for black lines; blue/black degrees; self-lines and
// converse-merged self-loops; lines by node type and by symmetric-core depth; and, for every blue
// line, the closed member's endpoint layer k = s_(n-1) (the layer (n-1-k, k) of the detour
// program), with the test whether its two classes lie in the same black component.
//
// Memory: a uint32 class id per tiling (2^28 * 4 bytes = 1 GB at n = 9).
// Build: g++ -O2 -std=c++17 -o tlm9 tiling_line_metagraph_n9_20261020.cpp ; run: ./tlm9 n
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <numeric>
#include <cstring>
#include <chrono>
#include <climits>
using namespace std;

static int n, m;
static int cellI[64], cellJ[64], cellIdx[16][16], sigmaCell[64];
static uint64_t deltaMask[16]; static int maxDepth;
static uint32_t ALL;

struct Canon { uint64_t best; int aut; };

static inline void buildAdj(uint32_t t, uint16_t adj[16]) {
    for (int u = 0; u < n; u++) adj[u] = 0;
    for (int i = 0; i + 1 < n; i++) adj[i] |= 1u << (i + 1);
    for (int k = 0; k < m; k++) { int i = cellI[k], j = cellJ[k]; if ((t >> k) & 1) adj[i] |= 1u << j; else adj[j] |= 1u << i; }
}

// ---- canonical form by colour refinement + individualisation ----
struct Canoniser {
    uint16_t adj[16]; int nn;
    uint64_t best; int autc;
    // refine colours (array of ints, values 0..nn-1 meaning colour classes ordered) until stable
    void refine(int col[16]) {
        for (int iter = 0; iter < 16; iter++) {
            // signature: (col[v], sorted multiset of colours of out-neighbours) encoded as a sortable key
            uint64_t key[16]; int cnt[16][16];
            for (int v = 0; v < nn; v++) {
                int c[16] = {0};
                for (int w = 0; w < nn; w++) if ((adj[v] >> w) & 1) c[col[w]]++;
                uint64_t k = (uint64_t)col[v];
                for (int x = 0; x < nn; x++) k = k * 16 + c[x];   // nn <= 9: counts < 16, 9 digits of 4 bits fits in 64 bits with the colour prefix
                key[v] = k; (void)cnt;
            }
            // new colours = rank of key among distinct keys
            uint64_t sorted[16]; for (int v = 0; v < nn; v++) sorted[v] = key[v];
            sort(sorted, sorted + nn); int nd = (int)(unique(sorted, sorted + nn) - sorted);
            int newcol[16]; int oldnd = 0; { int s2[16]; for (int v = 0; v < nn; v++) s2[v] = col[v]; sort(s2, s2 + nn); oldnd = (int)(unique(s2, s2 + nn) - s2); }
            for (int v = 0; v < nn; v++) newcol[v] = (int)(lower_bound(sorted, sorted + nd, key[v]) - sorted);
            bool same = true; for (int v = 0; v < nn; v++) if (newcol[v] != col[v]) { same = false; break; }
            for (int v = 0; v < nn; v++) col[v] = newcol[v];
            if (same || nd == oldnd) { if (nd == oldnd && !same) continue; break; }
            if (nd == nn) break;
        }
    }
    void leaf(const int col[16]) {
        int perm[16]; for (int v = 0; v < nn; v++) perm[col[v]] = v;   // position = colour (discrete)
        uint64_t w = 0;
        for (int a = 0; a < nn; a++) for (int b = a + 1; b < nn; b++) w = (w << 1) | ((adj[perm[a]] >> perm[b]) & 1u);
        if (w < best) { best = w; autc = 1; } else if (w == best) autc++;
    }
    void search(int col[16]) {
        refine(col);
        int nd = 0; { int s2[16]; for (int v = 0; v < nn; v++) s2[v] = col[v]; sort(s2, s2 + nn); nd = (int)(unique(s2, s2 + nn) - s2); }
        if (nd == nn) { leaf(col); return; }
        // first non-singleton cell (smallest colour with >= 2 members)
        int cnt[16] = {0}; for (int v = 0; v < nn; v++) cnt[col[v]]++;
        int cell = 0; while (cnt[cell] < 2) cell++;
        for (int v = 0; v < nn; v++) if (col[v] == cell) {
            int c2[16];
            for (int u = 0; u < nn; u++) c2[u] = col[u] + (col[u] > cell ? 1 : 0);   // open a slot after `cell`
            c2[v] = cell;                // v keeps the colour, the others of the cell move up
            for (int u = 0; u < nn; u++) if (u != v && col[u] == cell) c2[u] = cell + 1;
            search(c2);
        }
    }
    Canon run(const uint16_t a[16], int n_) {
        nn = n_; for (int v = 0; v < nn; v++) adj[v] = a[v];
        best = ~0ULL; autc = 0;
        int col[16]; int sc[16];
        for (int v = 0; v < nn; v++) sc[v] = __builtin_popcount(adj[v]);
        int s2[16]; for (int v = 0; v < nn; v++) s2[v] = sc[v]; sort(s2, s2 + nn); int nd = (int)(unique(s2, s2 + nn) - s2);
        for (int v = 0; v < nn; v++) col[v] = (int)(lower_bound(s2, s2 + nd, sc[v]) - s2);
        search(col);
        return {best, autc};
    }
};

static long long hamPaths(const uint16_t adj[16]) {
    int N = 1 << n; vector<long long> dp((size_t)N * n, 0);
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
    for (int mask = 1; mask < N; mask++) for (int v = 0; v < n; v++) {
        long long c = dp[(size_t)mask * n + v]; if (!c) continue;
        uint16_t nb = adj[v] & ~mask;
        while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(mask | (1 << w)) * n + w] += c; }
    }
    long long tot = 0; for (int v = 0; v < n; v++) tot += dp[(size_t)(N - 1) * n + v];
    return tot;
}
static inline uint32_t sigmaOf(uint32_t t) { uint32_t s = 0; for (int k = 0; k < m; k++) if ((t >> k) & 1) s |= 1u << sigmaCell[k]; return s; }
static inline int depthOf(uint32_t t) { uint32_t s = sigmaOf(t); for (int k = 0;; k++) if (((t ^ s) & deltaMask[k]) == 0) return k; }

struct UF { vector<int> p; UF(int N) : p(N) { iota(p.begin(), p.end(), 0); } int f(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; } void u(int a, int b) { a = f(a); b = f(b); if (a != b) p[a] = b; } };

int main(int argc, char **argv) {
    n = argc > 1 ? atoi(argv[1]) : 7;
    if (n < 4 || n > 9) { fprintf(stderr, "n in 4..9\n"); return 1; }
    m = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cellIdx[i][j] = -1;
    for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { cellI[m] = i; cellJ[m] = j; cellIdx[i][j] = m; m++; }
    for (int k = 0; k < m; k++) sigmaCell[k] = cellIdx[n - 1 - cellJ[k]][n - 1 - cellI[k]];
    maxDepth = 0;
    for (int k = 0; k < 16; k++) { deltaMask[k] = 0; for (int c = 0; c < m; c++) if (cellI[c] >= k && cellJ[c] <= n - 1 - k) deltaMask[k] |= 1ull << c; if (deltaMask[k] == 0 && maxDepth == 0) maxDepth = k; }
    uint32_t NT = 1u << m; ALL = NT - 1;
    int apex = cellIdx[0][n - 1];
    auto t0 = chrono::steady_clock::now();
    printf("# line metagraph census (refinement canoniser), n = %d, m = %d, tilings = %u\n", n, m, NT); fflush(stdout);

    vector<uint32_t> cls(NT, 0xFFFFFFFFu);
    unordered_map<uint64_t, uint32_t> byCanon; byCanon.reserve(1 << 19);
    vector<int> aut; vector<long long> H; vector<uint32_t> rep; vector<int> tc; vector<int> bm; vector<int> closedCnt;   // closed presentations (apex backward): > 0 iff strong (THM-4729 L2)
    Canoniser C; uint16_t adj[16];
    for (uint32_t t = 0; t < NT; t++) {
        buildAdj(t, adj);
        Canon c = C.run(adj, n);
        auto it = byCanon.find(c.best); uint32_t id;
        if (it == byCanon.end()) { id = (uint32_t)aut.size(); byCanon.emplace(c.best, id); aut.push_back(c.aut); H.push_back(hamPaths(adj)); rep.push_back(t); tc.push_back(0); bm.push_back(0); closedCnt.push_back(0); }
        else id = it->second;
        cls[t] = id; tc[id]++; if (sigmaOf(t) == t) bm[id]++; if (!((t >> apex) & 1)) closedCnt[id]++;
        if ((t & 0xFFFFFF) == 0 && t) { double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count(); fprintf(stderr, "  %u / %u tilings, %zu classes, %.0f s\n", t, NT, aut.size(), el); }
    }
    uint32_t Cn = (uint32_t)aut.size();
    vector<uint32_t> sigmaClass(Cn); int SC = 0;
    for (uint32_t c = 0; c < Cn; c++) { sigmaClass[c] = cls[sigmaOf(rep[c])]; if (sigmaClass[c] == c) SC++; }
    long long sumTc = 0; int badFibre = 0, badBlue = 0;
    for (uint32_t c = 0; c < Cn; c++) { sumTc += tc[c]; if ((long long)tc[c] * aut[c] != H[c]) badFibre++; bool sc = sigmaClass[c] == c; if ((sc && bm[c] % 2 == 0) || (!sc && bm[c] != 0)) badBlue++; }
    auto type = [&](uint32_t c) { if (sigmaClass[c] != c) return 2; return bm[c] == tc[c] ? 0 : 1; };
    int PB = 0, MX = 0, PK = 0; for (uint32_t c = 0; c < Cn; c++) { int ty = type(c); if (ty == 0) PB++; else if (ty == 1) MX++; else PK++; }
    double el1 = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    printf("classes %u (SC %d, merged %u); PB/MX/PK %d/%d/%d; sum tc %lld (= 2^m: %s); fibre violations %d; blue-parity violations %d; pass 1 %.0f s\n", Cn, SC, (Cn + SC) / 2, PB, MX, PK, sumTc, sumTc == (long long)NT ? "yes" : "NO", badFibre, badBlue, el1); fflush(stdout);

    // ---- pass 2: lines ----
    UF ufAll(Cn), ufBlue(Cn), ufBlack(Cn);
    vector<int> blueDeg(Cn, 0), blackDeg(Cn, 0), selfB(Cn, 0), selfK(Cn, 0);
    vector<long long> minHb(Cn, LLONG_MAX), maxHb(Cn, 0);   // H-descent: extreme H over black neighbours (other classes)
    long long blueLines = 0, blackLines = 0, mselfB = 0, mselfK = 0, pathRev = 0; long long byType[3][3][2]; memset(byType, 0, sizeof byType);
    long long linesByDepth[16] = {0};
    for (uint32_t t = 0; t < NT; t++) {
        uint32_t u = t ^ ALL; if (u < t) continue;
        uint32_t c1 = cls[t], c2 = cls[u]; bool blue = sigmaOf(t) == t; int d = depthOf(t); linesByDepth[d]++;
        int ta = type(c1), tb = type(c2); if (ta > tb) swap(ta, tb); byType[ta][tb][blue ? 0 : 1]++;
        ufAll.u(c1, c2);
        if (blue) { blueLines++; ufBlue.u(c1, c2); if (c1 == c2) { selfB[c1]++; blueDeg[c1] += 2; } else { blueDeg[c1]++; blueDeg[c2]++; } }
        else { blackLines++; ufBlack.u(c1, c2); if (c1 == c2) { selfK[c1]++; blackDeg[c1] += 2; } else { blackDeg[c1]++; blackDeg[c2]++; minHb[c1] = min(minHb[c1], H[c2]); minHb[c2] = min(minHb[c2], H[c1]); maxHb[c1] = max(maxHb[c1], H[c2]); maxHb[c2] = max(maxHb[c2], H[c1]); } }
        bool merged = (c1 == c2) || (c2 == sigmaClass[c1]);
        if (merged) { if (blue) mselfB++; else mselfK++; }
        if (c2 == sigmaClass[c1] && c1 != c2) pathRev++;
    }
    int badDeg = 0; for (uint32_t c = 0; c < Cn; c++) if (blueDeg[c] != bm[c] || blackDeg[c] != tc[c] - bm[c]) badDeg++;
    // components
    auto comps = [&](UF &uf, bool onlySC, bool skipPB) { map<int,int> sz; for (uint32_t c = 0; c < Cn; c++) { if (onlySC && sigmaClass[c] != c) continue; if (skipPB && type(c) == 0) continue; sz[uf.f(c)]++; } int big = 0; for (auto &kv : sz) big = max(big, kv.second); return make_pair((int)sz.size(), big); };
    auto cAll = comps(ufAll, false, false), cBlue = comps(ufBlue, true, false), cBlack = comps(ufBlack, false, false), cBlackNoPB = comps(ufBlack, false, true);
    long long totalSB = 0, totalSK = 0; for (uint32_t c = 0; c < Cn; c++) { totalSB += selfB[c]; totalSK += selfK[c]; }
    printf("lines: blue %lld, black %lld; degree law violations %d; self-lines blue %lld black %lld; merged self-loops blue %lld black %lld; reversal-neutral lines %lld\n", blueLines, blackLines, badDeg, totalSB, totalSK, mselfB, mselfK, pathRev);
    printf("components: all lines %d (largest %d); blue on SC classes %d (largest %d); black %d (largest %d); black on non-pure-blue classes %d (largest %d of %d)\n", cAll.first, cAll.second, cBlue.first, cBlue.second, cBlack.first, cBlack.second, cBlackNoPB.first, cBlackNoPB.second, (int)Cn - PB);
    // ---- H-descent statistics: a non-pure-blue class is a sink when no black neighbour has strictly smaller H ----
    { long long sinks = 0, sinksSC = 0, sinkDegSum = 0; map<long long,int> sinkH; long long upOnly = 0;
      for (uint32_t c = 0; c < Cn; c++) { if (type(c) == 0) continue; bool hasNb = minHb[c] != LLONG_MAX; if (!hasNb || minHb[c] >= H[c]) { sinks++; if (sigmaClass[c] == c) sinksSC++; sinkDegSum += blackDeg[c]; sinkH[H[c]]++; } if (hasNb && minHb[c] > H[c]) upOnly++; }
      printf("H-descent: non-pure-blue sinks (no black neighbour with smaller H) %lld (self-converse %lld; strictly uphill only %lld); sink H histogram:", sinks, sinksSC, upOnly); int shown = 0; for (auto &kv : sinkH) { if (shown++ < 14) printf(" %lld:%d", kv.first, kv.second); } printf("%s\n", (int)sinkH.size() > 14 ? " ..." : ""); }
    const char *tn[3] = {"PB", "MX", "PK"}; printf("lines by node types (blue, black):"); for (int a = 0; a < 3; a++) for (int b = a; b < 3; b++) printf(" (%s,%s)=(%lld,%lld)", tn[a], tn[b], byType[a][b][0], byType[a][b][1]); printf("\n");
    printf("lines by depth:"); for (int d = 0; d <= maxDepth; d++) printf(" d%d=%lld", d, linesByDepth[d]); printf("\n");
    // ---- blue lines by closed endpoint layer k = s_(n-1) of the closed member; black-connectivity of their classes ----
    map<int, array<long long, 4>> layer;   // k -> [lines, cross-class lines, cross-class with both non-PB, of those black-connected]
    for (uint32_t t = 0; t < NT; t++) {
        uint32_t u = t ^ ALL; if (u < t) continue;
        if (sigmaOf(t) != t) continue;
        uint32_t closed = ((t >> apex) & 1) ? u : t;   // closed member: apex tile backward (arc n-1 -> 0)
        buildAdj(closed, adj); int k = __builtin_popcount(adj[n - 1]);
        uint32_t c1 = cls[t], c2 = cls[u]; auto &L = layer[k]; L[0]++;
        if (c1 != c2) { L[1]++; if (type(c1) != 0 && type(c2) != 0) { L[2]++; if (ufBlack.f(c1) == ufBlack.f(c2)) L[3]++; } }
    }
    printf("blue lines by closed-member endpoint layer k = s_(n-1) (pattern (n-1-k, k)): k: lines / cross-class / both non-PB / black-connected\n");
    for (auto &kv : layer) printf("   k=%d: %lld / %lld / %lld / %lld\n", kv.first, kv.second[0], kv.second[1], kv.second[2], kv.second[3]);
    // ---- per-class dump ----
    char fname[160]; snprintf(fname, sizeof fname, "tiling_line_metagraph_n9_20261020_classes_n%d.tsv", n);
    FILE *f = fopen(fname, "w"); fprintf(f, "class\tH\tAut\ttc\tbm\tSC\ttype\tselfBlue\tselfBlack\tblackComponent\tminHblack\tmaxHblack\tclosed\n");
    for (uint32_t c = 0; c < Cn; c++) fprintf(f, "%u\t%lld\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%lld\t%lld\t%d\n", c, H[c], aut[c], tc[c], bm[c], sigmaClass[c] == c ? 1 : 0, tn[type(c)], selfB[c], selfK[c], ufBlack.f(c), minHb[c] == LLONG_MAX ? -1 : minHb[c], maxHb[c], closedCnt[c]);
    fclose(f);
    double el2 = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    printf("done in %.0f s; per-class file %s\n", el2, fname);
    return 0;
}
