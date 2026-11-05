// class_census_audit.cpp -- independent audit code.  Isomorphism classes of tournaments of order <= NMAX (<= 8)
// by vertex extension + canonical form (colour refinement, then permutations within colour classes, lexicographic
// minimum of the adjacency matrix).  For each class: H (subset DP), strong, self-converse, |Aut|.
// Since Busch's f(9) = 75, every strong class with H < 75 has order <= 8, so m(h) is exact for h < 75.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
using namespace std;
typedef unsigned long long u64;
int NMAX = 8;
struct Tour { int n; int adj[8]; };
static u64 countH(const Tour& T) {
    int n = T.n; static u64 dp[256][8]; int full = (1 << n) - 1;
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { u64 c = dp[S][v]; if (!c) continue;
        int nb = T.adj[v] & ~S; while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[S | (1 << w)][w] += c; } }
    u64 h = 0; for (int v = 0; v < n; v++) h += dp[full][v]; return h; }
static bool isStrong(const Tour& T) {
    int n = T.n; if (n == 1) return true; int full = (1 << n) - 1;
    int reach = 1, frontier = 1; while (frontier) { int nf = 0; for (int v = 0; v < n; v++) if ((frontier >> v) & 1) nf |= T.adj[v]; nf &= ~reach; reach |= nf; frontier = nf; }
    if (reach != full) return false;
    int radj[8] = {0}; for (int v = 0; v < n; v++) for (int w = 0; w < n; w++) if ((T.adj[v] >> w) & 1) radj[w] |= 1 << v;
    reach = 1; frontier = 1; while (frontier) { int nf = 0; for (int v = 0; v < n; v++) if ((frontier >> v) & 1) nf |= radj[v]; nf &= ~reach; reach |= nf; frontier = nf; }
    return reach == full; }
// canonical form: returns min matrix code; aut = number of colour-respecting permutations attaining it
static int cn; static const Tour* cT; static vector<vector<int>> classes; static int pos[8]; static int usedMask; static u64 best; static u64 autCount; static int perm[8];
static void rec(int ci, int k) {
    if (ci == (int)classes.size()) { // perm[old] = new position
        u64 code = 0; for (int a = 0; a < cn; a++) for (int b = 0; b < cn; b++) if ((cT->adj[a] >> b) & 1) code |= 1ULL << (perm[a] * cn + perm[b]);
        if (code < best) { best = code; autCount = 1; } else if (code == best) autCount++; return; }
    const vector<int>& cl = classes[ci]; if (k == (int)cl.size()) { rec(ci + 1, 0); return; }
    int base = 0; for (int i = 0; i < ci; i++) base += classes[i].size();
    for (int v : cl) if (!((usedMask >> v) & 1)) { usedMask |= 1 << v; perm[v] = base + k; rec(ci, k + 1); usedMask &= ~(1 << v); } }
static u64 canon(const Tour& T, u64* aut) {
    int n = T.n; vector<int> col(n); for (int v = 0; v < n; v++) col[v] = __builtin_popcount(T.adj[v]);
    // refinement
    while (true) { vector<pair<vector<int>, int>> sig(n); for (int v = 0; v < n; v++) { vector<int> s; s.push_back(col[v]); vector<int> outc, inc; for (int w = 0; w < n; w++) if ((T.adj[v] >> w) & 1) outc.push_back(col[w]); else if (w != v) inc.push_back(col[w]); sort(outc.begin(), outc.end()); sort(inc.begin(), inc.end()); s.push_back(-1); for (int x : outc) s.push_back(x); s.push_back(-2); for (int x : inc) s.push_back(x); sig[v] = {s, v}; }
        vector<pair<vector<int>, int>> ss = sig; sort(ss.begin(), ss.end()); vector<int> ncol(n); int c = 0; for (int i = 0; i < n; i++) { if (i > 0 && ss[i].first != ss[i - 1].first) c++; ncol[ss[i].second] = c; }
        int oldc = *max_element(col.begin(), col.end()), newc = c; col = ncol; if (newc == oldc) break; }
    int ncl = *max_element(col.begin(), col.end()) + 1; classes.assign(ncl, vector<int>()); for (int v = 0; v < n; v++) classes[col[v]].push_back(v);
    cn = n; cT = &T; usedMask = 0; best = ~0ULL; autCount = 0; rec(0, 0); if (aut) *aut = autCount; return best; }
static Tour fromCode(int n, u64 code) { Tour T; T.n = n; for (int a = 0; a < n; a++) { T.adj[a] = 0; for (int b = 0; b < n; b++) if ((code >> (a * n + b)) & 1) T.adj[a] |= 1 << b; } return T; }
int main(int argc, char** argv) {
    if (argc > 1) NMAX = atoi(argv[1]);
    vector<u64> cur; { Tour T; T.n = 1; T.adj[0] = 0; cur.push_back(canon(T, 0)); }
    map<u64, vector<int>> mByH; // h -> list of orders of strong classes (h < 75 complete when NMAX >= 8)
    map<u64, int> mSC; // h -> number of self-converse strong classes
    long long expectedClasses[9] = {0, 1, 1, 2, 4, 12, 56, 456, 6880}, expectedStrong[9] = {0, 1, 0, 1, 1, 6, 35, 353, 6008};
    for (int n = 1; n <= NMAX; n++) {
        if (n > 1) { set<u64> next; for (u64 code : cur) { Tour T = fromCode(n - 1, code); for (int pat = 0; pat < (1 << (n - 1)); pat++) { Tour U; U.n = n; for (int a = 0; a < n - 1; a++) { U.adj[a] = T.adj[a]; if (!((pat >> a) & 1)) U.adj[a] |= 1 << (n - 1); } U.adj[n - 1] = pat; next.insert(canon(U, 0)); } } cur.assign(next.begin(), next.end()); }
        long long nstrong = 0; printf("order %d: %zu classes (expected %lld)", n, cur.size(), expectedClasses[n]);
        for (u64 code : cur) { Tour T = fromCode(n, code); bool st = isStrong(T); if (!st) continue; nstrong++; u64 h = countH(T);
            Tour C; C.n = n; for (int a = 0; a < n; a++) { C.adj[a] = 0; for (int b = 0; b < n; b++) if (a != b && !((T.adj[a] >> b) & 1)) C.adj[a] |= 1 << b; }
            u64 aut; u64 cc = canon(C, 0); u64 own = canon(T, &aut); bool sc = (cc == own); mByH[h].push_back(n); if (sc) mSC[h]++;
            if (n <= 5 || (h <= 25)) printf("\n   strong class order %d: H=%llu self-converse=%d |Aut|=%llu adj=", n, h, (int)sc, aut), ({ for (int a = 0; a < n; a++) printf("%x", T.adj[a]); 0; });
            if (aut && h % aut != 0) printf("\n   WARNING |Aut| does not divide H"); }
        printf("\n   strong classes: %lld (expected %lld)\n", nstrong, expectedStrong[n]); }
    printf("\nm(h) over strong classes of orders <= %d (exact for h < f(%d)):\n", NMAX, NMAX + 1);
    for (auto& kv : mByH) { if (kv.first >= 200) break; printf("  h=%llu m=%zu orders=", kv.first, kv.second.size()); for (int o : kv.second) printf("%d,", o); printf(" selfconverse=%d\n", mSC[kv.first]); }
    printf("odd h < 75 with m(h) = 0:"); for (u64 h = 1; h < 75; h += 2) if (!mByH.count(h)) printf(" %llu", h); printf("\n");
    return 0; }
