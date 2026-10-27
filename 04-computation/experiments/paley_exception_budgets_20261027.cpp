// paley_exception_budgets_20261027.cpp -- opus-2026-10-27-S63
// Exact exception budgets of Paley clocks P_q (vertices Z/q, i -> j iff j - i is a nonzero square, q prime, q = 3 mod 4)
// and of the first substitutions at the {7, 21} orders:
//   f(T) = minimum feedback-arc count over total orders (subset DP, 2^q states; q <= 23 by default, q = 31 needs 2 GB),
//   b(T) = minimum cycle debt beta(D) = |E(D)| - |V| + components(D) of the backward-pair graph D over all q! orders
//          (brute force, q <= 11), together with the number of orders attaining b and the edge counts they use,
//   the lower bound ceil(f/(q-1)) on the number of native complement moves from the transitive tournament,
//   session73's one-step forest ceiling 2q - 4,
//   and session73's weighted substitution law G1 checked on P7[C3] and C3[P7] (21 vertices, subset DP).
// Build: g++ -O2 -std=c++17 -o peb paley_exception_budgets_20261027.cpp ; run: ./peb [q ...]   (default 7 11 19 23)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace std;
typedef unsigned int u32; typedef unsigned long long u64;

static int fasDP(int n, const vector<u32>& out) {            // minimum feedback-arc count
    u64 N = 1ULL << n; vector<uint8_t> f(N); f[0] = 0;
    for (u64 S = 1; S < N; S++) { int best = 255; u64 t = S; while (t) { int v = __builtin_ctzll(t); t &= t - 1; u64 R = S ^ (1ULL << v); int c = f[R] + __builtin_popcount(out[v] & (u32)R); if (c < best) best = c; } f[S] = (uint8_t)best; }
    return f[N - 1]; }
static vector<u32> paley(int q) { vector<u32> out(q, 0); vector<int> sq(q, 0); for (int x = 1; x < q; x++) sq[(x * x) % q] = 1;
    for (int i = 0; i < q; i++) for (int j = 0; j < q; j++) if (i != j && sq[((j - i) % q + q) % q]) out[i] |= 1u << j; return out; }
static vector<u32> c3() { vector<u32> out(3, 0); out[0] = 1u << 1; out[1] = 1u << 2; out[2] = 1u << 0; return out; }
static vector<u32> subst(const vector<u32>& Q, const vector<u32>& T) {   // Q[T,...,T], vertex (a, b) -> a*|T| + b
    int q = (int)Q.size(), m = (int)T.size(), n = q * m; vector<u32> out(n, 0);
    for (int a = 0; a < q; a++) for (int b = 0; b < m; b++) for (int a2 = 0; a2 < q; a2++) for (int b2 = 0; b2 < m; b2++) { if (a == a2 && b == b2) continue;
        bool arc = (a != a2) ? (Q[a] >> a2 & 1) : (T[b] >> b2 & 1); if (arc) out[a * m + b] |= 1u << (a2 * m + b2); }
    return out; }
struct UF { vector<int> p; UF(int n) : p(n) { iota(p.begin(), p.end(), 0); } int f(int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; } bool u(int a, int b) { a = f(a); b = f(b); if (a == b) return false; p[a] = b; return true; } };
static void cycleDebt(int n, const vector<u32>& out, int& bmin, u64& nb, int& eAtB, int& eMin) {
    vector<int> perm(n); iota(perm.begin(), perm.end(), 0); bmin = 1 << 30; nb = 0; eAtB = 1 << 30; eMin = 1 << 30;
    do { int E = 0, comp = n; UF uf(n);
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { int u = perm[i], w = perm[j]; if (out[w] >> u & 1) { E++; if (uf.u(u, w)) comp--; } }
        int beta = E - n + comp; if (E < eMin) eMin = E;
        if (beta < bmin) { bmin = beta; nb = 1; eAtB = E; } else if (beta == bmin) { nb++; if (E < eAtB) eAtB = E; }
    } while (next_permutation(perm.begin(), perm.end())); }

int main(int argc, char** argv) {
    vector<int> qs; for (int i = 1; i < argc; i++) qs.push_back(atoi(argv[i])); if (qs.empty()) qs = {7, 11, 19, 23};
    for (int q : qs) { vector<u32> P = paley(q); int f = fasDP(q, P);
        printf("P_%d: f = %d (minimum feedback arcs); lower bound on complement moves ceil(f/(q-1)) = %d; forest one-step ceiling 2q-4 = %d; %s\n", q, f, (f + q - 2) / (q - 1), 2 * q - 4, f > 2 * q - 4 ? "beyond one step from any forest frame" : "within the one-step ceiling");
        if (q <= 11) { int b, eAtB, eMin; u64 nb; cycleDebt(q, P, b, nb, eAtB, eMin); printf("P_%d: cycle debt b = %d attained by %llu of %d! orders (fewest backward edges among them %d); edge minimum over all orders %d (= f: %s)\n", q, b, nb, q, eAtB, eMin, eMin == f ? "yes" : "NO"); }
        fflush(stdout); }
    // G1 at the {7, 21} orders
    vector<u32> P7 = paley(7), C3 = c3(); int f7 = fasDP(7, P7), f3 = fasDP(3, C3);
    vector<u32> A = subst(P7, C3), B = subst(C3, P7); int fA = fasDP(21, A), fB = fasDP(21, B);
    printf("G1 check: f(C3) = %d, f(P7) = %d; f(P7[C3]) = %d vs G1 7*f(C3) + 9*f(P7) = %d; f(C3[P7]) = %d vs G1 3*f(P7) + 49*f(C3) = %d\n", f3, f7, fA, 7 * f3 + 9 * f7, fB, 3 * f7 + 49 * f3);
    return 0;
}
