// hankel_block_audit.cpp -- independent audit code.  Lower bound on the Hankel rank of the sub-language series of
// G(0,2) (modes 0..3) WITHOUT the author's automaton: words are decoded directly into chain tournaments (vertex i -> j
// for i < j except the terminated backward arcs), H(T_uv) is computed by subset DP, and the rank of the block
// [H(T_uv)]_{|u| <= Lp, |v| <= Ls} is computed modulo two primes.  Letter semantics follow the engine's lettersFor(c):
// for T = 0..2^c-1, b = 0..W-(c-|T|), filtered by the mode; a letter index out of range for the current c gives 0.
// Usage: hankel mode Lp Ls
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
using namespace std;
typedef unsigned long long u64; typedef long long i64;
static int MODE, W = 2;
static bool allowed(int c, int T, int b) { if (MODE == 0) return true; if (b > 1 || __builtin_popcount(T) > 1) return false;
    if (MODE == 2 && c == 2 && T == 2) return false; if (MODE == 3 && c == 2 && T == 1) return false; return true; }
static vector<pair<int,int>> lettersFor(int c) { vector<pair<int,int>> L; for (int T = 0; T < (1 << c); T++) { int t = __builtin_popcount(T); for (int b = 0; c - t + b <= W; b++) if (allowed(c, T, b)) L.push_back({T, b}); } return L; }
struct Word { vector<int> li; };
// decode: returns false if invalid; fills backward arcs (tail, head) and final c
static bool decode(const vector<int>& w, vector<pair<int,int>>& arcs, int& cOut) {
    vector<int> heads; arcs.clear(); // pending slots in creation order, each = head vertex
    for (size_t v = 0; v < w.size(); v++) { int c = heads.size(); vector<pair<int,int>> L = lettersFor(c); if (w[v] < 0 || w[v] >= (int)L.size()) return false;
        int T = L[w[v]].first, b = L[w[v]].second; vector<int> nh; vector<int> termHeads;
        for (int i = 0; i < c; i++) { if ((T >> i) & 1) termHeads.push_back(heads[i]); else nh.push_back(heads[i]); }
        sort(termHeads.begin(), termHeads.end()); for (size_t i = 1; i < termHeads.size(); i++) if (termHeads[i] == termHeads[i - 1]) return false; // double arc
        for (int h : termHeads) arcs.push_back({(int)v, h});
        for (int j = 0; j < b; j++) nh.push_back((int)v); heads = nh; }
    cOut = heads.size(); return true; }
static u64 countH(int n, const vector<pair<int,int>>& arcs) {
    if (n == 0) return 1; vector<int> adj(n, 0); for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) adj[i] |= 1 << j;
    for (auto& a : arcs) { adj[a.second] &= ~(1 << a.first); adj[a.first] |= 1 << a.second; }
    int full = (1 << n) - 1; vector<u64> dp((size_t)(1 << n) * n, 0); for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { u64 c = dp[(size_t)S * n + v]; if (!c) continue; int nb = adj[v] & ~S; while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * n + w] += c; } }
    u64 h = 0; for (int v = 0; v < n; v++) h += dp[(size_t)full * n + v]; return h; }
static void genWords(int maxLen, vector<vector<int>>& out) { // all valid words (prefixes from c = 0) of length <= maxLen
    out.push_back({}); vector<vector<int>> frontier = {{}};
    for (int len = 1; len <= maxLen; len++) { vector<vector<int>> nxt; for (auto& w : frontier) { vector<pair<int,int>> arcs; int c; decode(w, arcs, c); int nl = lettersFor(c).size();
            for (int li = 0; li < nl; li++) { vector<int> w2 = w; w2.push_back(li); vector<pair<int,int>> a2; int c2; if (decode(w2, a2, c2)) { nxt.push_back(w2); out.push_back(w2); } } } frontier = nxt; } }
static int rankMod(vector<vector<i64>> M, i64 p) { int r = 0; int rows = M.size(); if (!rows) return 0; int cols = M[0].size();
    for (int c = 0; c < cols && r < rows; c++) { int piv = -1; for (int i = r; i < rows; i++) if (M[i][c] % p) { piv = i; break; } if (piv < 0) continue; swap(M[r], M[piv]);
        i64 inv = 1, base = ((M[r][c] % p) + p) % p, e = p - 2; while (e) { if (e & 1) inv = (__int128)inv * base % p; base = (__int128)base * base % p; e >>= 1; }
        for (int j = c; j < cols; j++) M[r][j] = (__int128)(((M[r][j] % p) + p) % p) * inv % p;
        for (int i = 0; i < rows; i++) if (i != r && M[i][c] % p) { i64 f = ((M[i][c] % p) + p) % p; for (int j = c; j < cols; j++) M[i][j] = ((M[i][j] - (__int128)f * M[r][j]) % p + p) % p; } r++; }
    return r; }
int main(int argc, char** argv) {
    MODE = atoi(argv[1]); int Lp = atoi(argv[2]), Ls = atoi(argv[3]);
    vector<vector<int>> pre; genWords(Lp, pre);
    // suffixes: all letter sequences of length <= Ls (validity judged on the concatenation); we take the sequences over
    // the maximal letter index range, i.e. every li in 0..maxLetters-1
    int maxL = 0; for (int c = 0; c <= W; c++) maxL = max(maxL, (int)lettersFor(c).size());
    vector<vector<int>> suf; suf.push_back({}); { vector<vector<int>> frontier = {{}}; for (int len = 1; len <= Ls; len++) { vector<vector<int>> nxt; for (auto& w : frontier) for (int li = 0; li < maxL; li++) { vector<int> w2 = w; w2.push_back(li); nxt.push_back(w2); suf.push_back(w2); } frontier = nxt; } }
    printf("mode %d: %zu prefixes (len <= %d), %zu suffix letter-strings (len <= %d, alphabet %d)\n", MODE, pre.size(), Lp, suf.size(), Ls, maxL);
    vector<vector<i64>> M(pre.size(), vector<i64>(suf.size(), 0)); long long nz = 0;
    for (size_t i = 0; i < pre.size(); i++) for (size_t j = 0; j < suf.size(); j++) { vector<int> w = pre[i]; w.insert(w.end(), suf[j].begin(), suf[j].end()); vector<pair<int,int>> arcs; int c; if (!decode(w, arcs, c)) continue; u64 h = countH(w.size(), arcs); M[i][j] = h; if (h) nz++; }
    int r1 = rankMod(M, 2147483647LL), r2 = rankMod(M, 2147483629LL);
    printf("nonzero entries %lld; block rank mod 2^31-1 = %d, mod 2^31-19 = %d (lower bound on the Hankel rank; words with pending slots included, empty word = 1)\n", nz, r1, r2);
    // also the rank restricted to closed words (uv ends with c = 0): the series on complete chain tournaments
    vector<vector<i64>> M0 = M; for (size_t i = 0; i < pre.size(); i++) for (size_t j = 0; j < suf.size(); j++) { vector<int> w = pre[i]; w.insert(w.end(), suf[j].begin(), suf[j].end()); vector<pair<int,int>> arcs; int c; if (!decode(w, arcs, c) || c != 0) M0[i][j] = 0; }
    printf("rank of the block restricted to words ending with c = 0: %d\n", rankMod(M0, 2147483647LL));
    return 0; }
