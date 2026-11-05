// hp.cpp -- independent auditor tool for the S79 packet (written from the definitions in the note, not from the engine).
// Hamiltonian-path counts by a brute-force subset dynamic programme; graded spine grammars G(k,w) decoded from the
// letter semantics (T = set of pending slots terminated here, b = slots opened here, A = apexes beating this vertex;
// letter index = position in the list of all (T,b,A), T outermost over all 2^c subsets, then b = 0..w-(c-|T|), then A).
// Modes:
//   batch                     : stdin lines "n <n*n adjacency chars>" (row i col j = 1 iff i -> j); prints H per line
//   word k w idx,idx,...      : decode a word of G(k,w) and print the tournament and H (subset DP)
//   spine n t>h t>h ...       : chain 0..n-1 with the listed reversed arcs; prints H
//   chain r                   : all-parts-1 upset chain on r+1 vertices (arcs i+1 -> i for all i); prints H
//   census k w X L [outfile]  : DFS over all words of G(k,w) with chain length <= L, prefix pruning H(prefix) <= X by the
//                               insertion lemma, incremental subset DP; prints attained odd values <= X (direct, with
//                               ordered-join products), misses, and writes one witness per value to outfile
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <random>
using namespace std;
typedef unsigned long long u64;

static u64 countHP(int n, const vector<u64>& in) {
    if (n == 0) return 1;
    vector<u64> dp((size_t)(1ULL << n) * n, 0);
    for (int v = 0; v < n; v++) dp[((size_t)(1ULL << v)) * n + v] = 1;
    u64 full = (1ULL << n) - 1;
    for (u64 S = 1; S <= full; S++) for (int v = 0; v < n; v++) if (S >> v & 1) {
        u64 T = S & ~(1ULL << v); if (!T) continue;
        u64 P = T & in[v]; u64 acc = 0;
        while (P) { int u = __builtin_ctzll(P); P &= P - 1; acc += dp[(size_t)T * n + u]; }
        dp[(size_t)S * n + v] = acc; }
    u64 h = 0; for (int v = 0; v < n; v++) h += dp[(size_t)full * n + v]; return h; }

struct Letter { int T, b, A; };
// my own enumeration from the note: T outermost over all subsets of the c pending slots, then b, then A
static vector<Letter> enumerateLetters(int c, int w, int k) {
    vector<Letter> L;
    for (int T = 0; T < (1 << c); T++) { int t = __builtin_popcount(T);
        for (int b = 0; c - t + b <= w; b++) for (int A = 0; A < (1 << k); A++) L.push_back({T, b, A}); }
    return L; }

struct Slot { int head; int group; };
// decode a word; returns false if invalid (two slots of one head group terminated together, index out of range) or incomplete
static bool decodeWord(int k, int w, const vector<int>& word, int& n, vector<u64>& in, vector<pair<int,int>>& arcs, vector<int>& apexbits, bool requireComplete = true) {
    vector<Slot> pending; int groups = 0; arcs.clear(); apexbits.clear();
    int m = (int)word.size(); n = m + k; in.assign(n, 0);
    // apexes are vertices 0..k-1 (x_1 = 0 -> x_2 = 1); chain vertex t is vertex k + t
    if (k == 2) in[1] |= 1ULL << 0;
    for (int t = 0; t < m; t++) {
        int c = (int)pending.size(); vector<Letter> L = enumerateLetters(c, w, k);
        int li = word[t]; if (li < 0 || li >= (int)L.size()) return false;
        Letter l = L[li];
        // group check
        { vector<int> seen; for (int i = 0; i < c; i++) if (l.T >> i & 1) { for (int g : seen) if (g == pending[i].group) return false; seen.push_back(pending[i].group); } }
        int d = k + t; u64 inmask = 0;
        for (int s = 0; s < t; s++) inmask |= 1ULL << (k + s);            // all earlier chain vertices beat d (forward arcs)
        vector<Slot> np;
        for (int i = 0; i < c; i++) {
            if (l.T >> i & 1) { int h = pending[i].head; inmask &= ~(1ULL << h); in[h] |= 1ULL << d; arcs.push_back({d - k, h - k}); }   // reversed arc d -> head
            else np.push_back(pending[i]); }
        for (int j = 0; j < l.b; j++) np.push_back({d, groups});
        if (l.b) groups++;
        for (int a = 0; a < k; a++) { if (l.A >> a & 1) inmask |= 1ULL << a; else in[a] |= 1ULL << d; }
        in[d] = inmask; apexbits.push_back(l.A); pending = np; }
    if (requireComplete && !pending.empty()) return false;
    return true; }

// ---------------- census with incremental DP ----------------
static int K, W, LMAX; static u64 X; static int NV;   // NV = max vertices = K + LMAX
static vector<u64> DP;                                 // DP[S * NV + v]
static vector<u64> IN;                                 // in-neighbour masks (stack-managed)
static vector<unsigned char> att; static vector<string> wit; static vector<int> witFa, witM;
static u64 nodes = 0, completed = 0; static vector<int> curWord; static int curFa = 0;
static string wordStr() { string s; for (size_t i = 0; i < curWord.size(); i++) { if (i) s += ","; s += to_string(curWord[i]); } return s; }
static void addVertex(int d) {   // IN[d] set; IN[v] for v<d updated with bit d where d -> v
    u64 lim = 1ULL << d;
    for (u64 S = 0; S < lim; S++) {
        u64 Sp = S | lim;
        // v = d
        u64 acc = 0; if (S == 0) acc = 1; else { u64 P = S & IN[d]; while (P) { int u = __builtin_ctzll(P); P &= P - 1; acc += DP[(size_t)S * NV + u]; } }
        DP[(size_t)Sp * NV + d] = acc;
        // v in S
        u64 Q = S; while (Q) { int v = __builtin_ctzll(Q); Q &= Q - 1; u64 R = Sp & ~(1ULL << v); u64 P = R & IN[v]; u64 a2 = 0;
            while (P) { int u = __builtin_ctzll(P); P &= P - 1; a2 += DP[(size_t)R * NV + u]; }
            DP[(size_t)Sp * NV + v] = a2; } } }
static u64 prefixH(int d) { u64 full = (1ULL << (d + 1)) - 1; u64 h = 0; for (int v = 0; v <= d; v++) h += DP[(size_t)full * NV + v]; return h; }
static void recordValue(u64 h) {
    if (h > X) return; int m = (int)curWord.size();
    if (!att[h] || curFa < witFa[h] || (curFa == witFa[h] && m < witM[h])) { att[h] = 1; witFa[h] = curFa; witM[h] = m; wit[h] = wordStr(); } }
static void dfs(int t, vector<Slot>& pending, int groups) {
    nodes++;
    int c = (int)pending.size(); vector<Letter> L = enumerateLetters(c, W, K);
    int d = K + t;
    for (size_t li = 0; li < L.size(); li++) {
        Letter l = L[li];
        bool ok = true; { vector<int> seen; for (int i = 0; i < c && ok; i++) if (l.T >> i & 1) { for (int g : seen) if (g == pending[i].group) ok = false; seen.push_back(pending[i].group); } }
        if (!ok) continue;
        bool slotFree = (l.T == 0 && l.b == 0);
        if (t == 0 && (K == 0 ? l.b == 0 : (slotFree && l.A == 0))) continue;   // padding (1-vertex word / leading source)
        int cNext = c - __builtin_popcount(l.T) + l.b;
        // arcs of d
        u64 inmask = 0; for (int s = 0; s < t; s++) inmask |= 1ULL << (K + s);
        vector<int> headsOut;
        for (int i = 0; i < c; i++) if (l.T >> i & 1) { int h = pending[i].head; inmask &= ~(1ULL << h); headsOut.push_back(h); }
        vector<int> apexOut;
        for (int a = 0; a < K; a++) { if (l.A >> a & 1) inmask |= 1ULL << a; else apexOut.push_back(a); }
        IN[d] = inmask; for (int h : headsOut) IN[h] |= 1ULL << d; for (int a : apexOut) IN[a] |= 1ULL << d;
        addVertex(d);
        u64 h = prefixH(d);
        if (h <= X) {
            curWord.push_back((int)li); curFa += __builtin_popcount(l.T);
            if (cNext == 0) { completed++; recordValue(h); }
            bool cont = (t + 1 < LMAX) && (K >= 1 || cNext > 0);
            if (cont) {
                vector<Slot> np; for (int i = 0; i < c; i++) if (!(l.T >> i & 1)) np.push_back(pending[i]);
                for (int j = 0; j < l.b; j++) np.push_back({d, groups});
                dfs(t + 1, np, groups + (l.b ? 1 : 0)); }
            curWord.pop_back(); curFa -= __builtin_popcount(l.T); }
        for (int hh : headsOut) IN[hh] &= ~(1ULL << d); for (int a : apexOut) IN[a] &= ~(1ULL << d); IN[d] = 0; } }

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "modes: batch | word k w idx,... | spine n t>h ... | chain r | census k w X L [outfile]\n"); return 1; }
    string mode = argv[1];
    if (mode == "batch") {
        int n; char buf[1 << 16];
        while (scanf("%d %65535s", &n, buf) == 2) { vector<u64> in(n, 0);
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (buf[i * n + j] == '1') in[j] |= 1ULL << i;
            printf("%llu\n", countHP(n, in)); }
        return 0; }
    if (mode == "word") { int k = atoi(argv[2]), w = atoi(argv[3]); vector<int> word; string s = argv[4]; size_t p = 0;
        while (p < s.size()) { size_t q = s.find(',', p); string tok = s.substr(p, q == string::npos ? string::npos : q - p); p = (q == string::npos) ? s.size() : q + 1; if (!tok.empty()) word.push_back(atoi(tok.c_str())); }
        int n; vector<u64> in; vector<pair<int,int>> arcs; vector<int> ab;
        if (!decodeWord(k, w, word, n, in, arcs, ab)) { printf("INVALID\n"); return 0; }
        printf("m=%d n=%d arcs(tail>head):", (int)word.size(), n); for (auto& a : arcs) printf(" %d>%d", a.first, a.second); printf(" apex:"); for (int x : ab) printf("%d", x);
        printf("\nH = %llu\n", countHP(n, in)); return 0; }
    if (mode == "spine") { int n = atoi(argv[2]); vector<u64> in(n, 0); for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) in[j] |= 1ULL << i;
        for (int a = 3; a < argc; a++) { int t, h; if (sscanf(argv[a], "%d>%d", &t, &h) == 2) { in[t] &= ~(1ULL << h); in[h] |= 1ULL << t; } }
        printf("H = %llu\n", countHP(n, in)); return 0; }
    if (mode == "chain") { int r = atoi(argv[2]); int n = r + 1; vector<u64> in(n, 0); for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) in[j] |= 1ULL << i;
        for (int i = 0; i + 1 < n; i++) { in[i + 1] &= ~(1ULL << i); in[i] |= 1ULL << (i + 1); }
        printf("chain r=%d (n=%d): H = %llu\n", r, n, countHP(n, in)); return 0; }
    if (mode == "randsearch") {   // random valid words (my own generator), H by the subset DP; which odd targets <= X get a witness
        int k = atoi(argv[2]), w = atoi(argv[3]); u64 Xr = strtoull(argv[4], 0, 10); int Lmin = atoi(argv[5]), Lmax = atoi(argv[6]); long trials = atol(argv[7]); u64 seed = strtoull(argv[8], 0, 10); const char* outf = argc > 9 ? argv[9] : 0;
        mt19937_64 rng(seed); vector<string> found(Xr + 1, ""); vector<unsigned char> have(Xr + 1, 0); have[1] = 1; long valid = 0;
        for (long t = 0; t < trials; t++) {
            int m = Lmin + (int)(rng() % (Lmax - Lmin + 1)); vector<int> word; vector<Slot> pending; int groups = 0; bool ok = true;
            for (int v = 0; v < m && ok; v++) { int c = (int)pending.size(); vector<Letter> L = enumerateLetters(c, w, k); bool last = (v == m - 1);
                vector<int> cands; for (size_t li = 0; li < L.size(); li++) { Letter l = L[li]; bool g = true; vector<int> seen; for (int i = 0; i < c && g; i++) if (l.T >> i & 1) { for (int x : seen) if (x == pending[i].group) g = false; seen.push_back(pending[i].group); } if (!g) continue;
                    int cn = c - __builtin_popcount(l.T) + l.b; if (last && cn != 0) continue; if (k == 0 && !last && cn == 0) continue; if (v == 0 && (k == 0 ? l.b == 0 : (l.T == 0 && l.b == 0 && l.A == 0))) continue; cands.push_back((int)li); }
                if (cands.empty()) { ok = false; break; }
                int li = cands[rng() % cands.size()]; Letter l = L[li]; word.push_back(li);
                vector<Slot> np; for (int i = 0; i < c; i++) if (!(l.T >> i & 1)) np.push_back(pending[i]); for (int j = 0; j < l.b; j++) np.push_back({v, groups}); if (l.b) groups++; pending = np; }
            if (!ok || !pending.empty()) continue;
            int n; vector<u64> in; vector<pair<int,int>> arcs; vector<int> ab; if (!decodeWord(k, w, word, n, in, arcs, ab)) continue; valid++;
            u64 h = countHP(n, in); if (h <= Xr && !have[h]) { have[h] = 1; string s; for (size_t i = 0; i < word.size(); i++) { if (i) s += ","; s += to_string(word[i]); } found[h] = s; } }
        u64 cnt = 0; vector<u64> miss; for (u64 x = 1; x <= Xr; x += 2) { if (have[x]) cnt++; else miss.push_back(x); }
        printf("randsearch G(%d,%d) X=%llu lengths %d..%d trials %ld (valid %ld): odd values found %llu of %llu; not found:", k, w, Xr, Lmin, Lmax, trials, valid, cnt, (Xr + 1) / 2); for (u64 x : miss) printf(" %llu", x); printf("\n");
        if (outf) { FILE* f = fopen(outf, "w"); for (u64 x = 1; x <= Xr; x += 2) if (have[x]) fprintf(f, "%llu = %s\n", x, found[x].c_str()); fclose(f); }
        return 0; }
    if (mode == "census") { K = atoi(argv[2]); W = atoi(argv[3]); X = strtoull(argv[4], 0, 10); LMAX = atoi(argv[5]); NV = K + LMAX;
        if (NV > 22) { fprintf(stderr, "too many vertices\n"); return 1; }
        DP.assign((size_t)(1ULL << NV) * NV, 0); IN.assign(NV, 0);
        att.assign(X + 1, 0); wit.assign(X + 1, ""); witFa.assign(X + 1, 1 << 30); witM.assign(X + 1, 1 << 30);
        att[1] = 1; wit[1] = "(transitive)"; witFa[1] = 0; witM[1] = 0;
        // apexes first
        if (K == 2) IN[1] |= 1ULL << 0;
        for (int a = 0; a < K; a++) addVertex(a);
        vector<Slot> pend; dfs(0, pend, 0);
        u64 cnt = 0; for (u64 x = 1; x <= X; x += 2) if (att[x]) cnt++;
        printf("census G(%d,%d) X=%llu L<=%d: %llu DFS nodes, %llu completed words; odd values attained directly: %llu of %llu\n", K, W, X, LMAX, nodes, completed, cnt, (X + 1) / 2);
        vector<unsigned char> cl(att.begin(), att.end()); for (u64 a = 3; a <= X; a += 2) if (cl[a]) for (u64 b = 3; a * b <= X; b += 2) if (att[b]) cl[a * b] = 1;
        u64 cc = 0; vector<u64> miss; for (u64 x = 1; x <= X; x += 2) { if (cl[x]) cc++; else miss.push_back(x); }
        printf("with products: %llu of %llu; %zu misses:", cc, (X + 1) / 2, miss.size()); for (size_t i = 0; i < miss.size(); i++) printf(" %llu", miss[i]); printf("\n");
        printf("misses per dyadic scale:"); for (int kk = 2; (1ULL << kk) <= X; kk++) { u64 lo = 1ULL << kk, hi = min(X, (1ULL << (kk + 1)) - 1), c1 = 0, tot = 0; for (u64 x = lo | 1; x <= hi; x += 2) { tot++; if (!cl[x]) c1++; } printf(" [2^%d:%llu/%llu]", kk, c1, tot); } printf("\n");
        map<int, int> byFa; for (u64 x = 1; x <= X; x += 2) if (att[x]) byFa[witFa[x]]++; printf("direct values by minimal number of reversed arcs (exact within L):"); for (auto& kv : byFa) printf(" fa%d:%d", kv.first, kv.second); printf("\n");
        map<int, int> byM; for (u64 x = 1; x <= X; x += 2) if (att[x]) byM[witM[x]]++; printf("direct values by shortest chain length (at minimal fa):"); for (auto& kv : byM) printf(" m%d:%d", kv.first, kv.second); printf("\n");
        if (argc > 6) { FILE* f = fopen(argv[6], "w"); for (u64 x = 1; x <= X; x += 2) if (att[x]) fprintf(f, "%llu = %s | fa=%d m=%d\n", x, wit[x].c_str(), witFa[x], witM[x]); fclose(f); }
        return 0; }
    return 0; }
