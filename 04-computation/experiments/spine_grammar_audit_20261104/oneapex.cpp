// oneapex.cpp -- auditor's independent one-apex spectrum census via the odd-cycle formula.
// One-apex tournament: transitive chain v_1 < ... < v_m, apex x with x -> v_i iff b_i = 1 (else v_i -> x).
// Every cycle passes through x (the chain is transitive) and no two odd cycles are vertex-disjoint, so by the odd-cycle
// collection formula H = 1 + 2 * #(odd cycles).  An odd cycle is x -> v_i -> (even subset of v_(i+1..j-1), increasing)
// -> v_j -> x with b_i = 1, b_j = 0, i < j: w(j - i) of them, w(1) = 1, w(d) = 2^(d-2) for d >= 2.
// Modes:  verify L : compare the formula with a subset DP for every word of length <= L
//         census X L : all words of length <= L, attained odd values <= X (direct, with products), dyadic misses, residues
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
using namespace std;
typedef unsigned long long u64;
static u64 wgt(int d) { return d == 1 ? 1 : (1ULL << (d - 2)); }
static u64 formulaH(const vector<int>& b) { int m = (int)b.size(); u64 cyc = 0;
    for (int i = 0; i < m; i++) if (b[i] == 1) for (int j = i + 1; j < m; j++) if (b[j] == 0) cyc += wgt(j - i);
    return 1 + 2 * cyc; }
static u64 countHP(int n, const vector<u64>& in) {
    vector<u64> dp((size_t)(1ULL << n) * n, 0); for (int v = 0; v < n; v++) dp[((size_t)(1ULL << v)) * n + v] = 1;
    u64 full = (1ULL << n) - 1;
    for (u64 S = 1; S <= full; S++) for (int v = 0; v < n; v++) if (S >> v & 1) { u64 T = S & ~(1ULL << v); if (!T) continue;
        u64 P = T & in[v]; u64 acc = 0; while (P) { int u = __builtin_ctzll(P); P &= P - 1; acc += dp[(size_t)T * n + u]; } dp[(size_t)S * n + v] = acc; }
    u64 h = 0; for (int v = 0; v < n; v++) h += dp[(size_t)full * n + v]; return h; }
static u64 dpH(const vector<int>& b) { int m = (int)b.size(), n = m + 1; vector<u64> in(n, 0);   // apex = vertex m
    for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) in[j] |= 1ULL << i;
    for (int i = 0; i < m; i++) { if (b[i]) in[i] |= 1ULL << m; else in[m] |= 1ULL << i; }
    return countHP(n, in); }
int main(int argc, char** argv) {
    string mode = argv[1];
    if (mode == "verify") { int L = atoi(argv[2]); long bad = 0, tot = 0;
        for (int m = 1; m <= L; m++) for (u64 code = 0; code < (1ULL << m); code++) { vector<int> b(m); for (int i = 0; i < m; i++) b[i] = (code >> i) & 1;
            u64 f = formulaH(b), d = dpH(b); tot++; if (f != d) { bad++; if (bad < 5) printf("MISMATCH m=%d code=%llu formula=%llu dp=%llu\n", m, code, f, d); } }
        printf("one-apex formula vs subset DP: %ld words of length <= %d, %ld mismatches\n", tot, L, bad); return 0; }
    if (mode == "census") { u64 X = strtoull(argv[2], 0, 10); int L = atoi(argv[3]); vector<unsigned char> att(X + 1, 0); att[1] = 1;
        // DFS over bits with incremental cycle count: adding bit 0 at position j adds sum_{i<j, b_i = 1} w(j - i)
        vector<int> b; u64 nodes = 0;
        struct R { static void go(vector<int>& b, int L, u64 X, u64 cyc, vector<unsigned char>& att, u64& nodes) {
            nodes++; u64 h = 1 + 2 * cyc; if (h > X) return; if (h <= X) att[h] = 1;
            if ((int)b.size() == L) return;
            int j = (int)b.size();
            // bit 1: no new cycles
            b.push_back(1); go(b, L, X, cyc, att, nodes); b.pop_back();
            // bit 0: new cycles ending at j
            u64 add = 0; for (int i = 0; i < j; i++) if (b[i] == 1) add += wgt(j - i);
            b.push_back(0); go(b, L, X, cyc + add, att, nodes); b.pop_back(); } };
        R::go(b, L, X, 0, att, nodes);
        u64 cnt = 0; for (u64 x = 1; x <= X; x += 2) if (att[x]) cnt++;
        printf("one-apex census X=%llu L<=%d: %llu nodes; odd values attained directly: %llu of %llu\n", X, L, nodes, cnt, (X + 1) / 2);
        vector<unsigned char> cl(att.begin(), att.end()); for (u64 a = 3; a <= X; a += 2) if (cl[a]) for (u64 c = 3; a * c <= X; c += 2) if (att[c]) cl[a * c] = 1;
        u64 cc = 0; for (u64 x = 1; x <= X; x += 2) if (cl[x]) cc++; printf("with products: %llu of %llu\n", cc, (X + 1) / 2);
        printf("direct misses per dyadic scale [2^j, 2^(j+1)):"); for (int j = 2; (1ULL << j) <= X; j++) { u64 lo = 1ULL << j, hi = min(X, (1ULL << (j + 1)) - 1), c1 = 0; for (u64 x = lo | 1; x <= hi; x += 2) if (!att[x]) c1++; printf(" %llu", c1); } printf("\n");
        printf("first direct misses:"); int shown = 0; for (u64 x = 1; x <= X && shown < 12; x += 2) if (!att[x]) { printf(" %llu", x); shown++; } printf("\n");
        int mods[4] = {16, 3, 5, 7}; for (int t = 0; t < 4; t++) { int q = mods[t]; vector<u64> miss(q, 0), tot(q, 0); for (u64 x = 1; x <= X; x += 2) { tot[x % q]++; if (!att[x]) miss[x % q]++; }
            printf("direct miss fraction by residue mod %d:", q); for (int r = 0; r < q; r++) if (tot[r]) printf(" %d:%.4f", r, (double)miss[r] / tot[r]); printf("\n"); }
        return 0; }
    return 0; }
