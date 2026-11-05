// spectrum_polynomial_congruence_20261105.cpp -- opus-2026-11-05-S81
// The spectrum polynomial P_n(x) = sum over all labelled tournaments T on {0..n-1} of x^H(T), where H = number of
// Hamiltonian paths, and the circulant polynomial C_n(x) = sum over circulant tournaments on Z_n (arc i -> j iff
// j - i in S, S a half-system of Z_n^*, n odd) of x^H.  Orbit counting under the rotation group Z_n (which acts on
// labelled tournaments by relabelling; H is an isomorphism invariant) gives, for every prime p dividing n, the
// congruence P_n(x) = F_{n,p}(x) (mod p), where F_{n,p} sums over the tournaments fixed by the rotation i -> i + n/p
// (for n = p prime: the circulants).  This program computes P_n exactly for n <= 7 (2^21 labelled tournaments at
// n = 7) and the fixed-tournament polynomials, and prints the coefficients, their residues, and the check.
// Usage: congruence n   (n <= 7)
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>
#include <cstdint>
using namespace std;
typedef unsigned long long u64;
static int n; static int pairIndex[8][8];
static u64 hamilton(u64 mask) {   // mask bit pairIndex[i][j] (i<j) = 1 means i -> j
    int full = (1 << n) - 1; vector<u64> dp((size_t)(1 << n) * n, 0);
    auto arc = [&](int a, int b) { if (a < b) return (mask >> pairIndex[a][b]) & 1; return ((mask >> pairIndex[b][a]) & 1) ^ 1ULL; };
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { u64 c = dp[(size_t)S * n + v]; if (!c) continue;
        for (int w = 0; w < n; w++) if (!(S >> w & 1) && arc(v, w)) dp[(size_t)(S | (1 << w)) * n + w] += c; }
    u64 h = 0; for (int v = 0; v < n; v++) h += dp[(size_t)full * n + v]; return h; }
int main(int argc, char** argv) {
    n = atoi(argv[1]); int P = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) pairIndex[i][j] = P++;
    map<u64, u64> Pn; map<int, map<u64, u64>> fixedPoly;    // per prime p | n
    vector<int> primes; for (int p = 2; p <= n; p++) { bool pr = true; for (int d = 2; d * d <= p; d++) if (p % d == 0) pr = false; if (pr && n % p == 0) primes.push_back(p); }
    for (u64 mask = 0; mask < (1ULL << P); mask++) {
        u64 h = hamilton(mask); Pn[h]++;
        for (int p : primes) { int r = n / p; bool fixed = true;   // invariant under i -> i + r (mod n)
            for (int i = 0; i < n && fixed; i++) for (int j = i + 1; j < n && fixed; j++) {
                int a = (i + r) % n, b = (j + r) % n; u64 ij = (mask >> pairIndex[i][j]) & 1;   // i -> j ?
                u64 ab; if (a < b) ab = (mask >> pairIndex[a][b]) & 1; else ab = ((mask >> pairIndex[b][a]) & 1) ^ 1ULL;
                if (ij != ab) fixed = false; }
            if (fixed) fixedPoly[p][h]++; } }
    u64 total = 0; for (auto& kv : Pn) total += kv.second;
    printf("n = %d: %llu labelled tournaments; spectrum polynomial P_%d(x) coefficients (h : count):\n", n, total, n);
    for (auto& kv : Pn) printf("  %llu : %llu\n", kv.first, kv.second);
    for (int p : primes) { printf("prime p = %d (rotation by %d): fixed tournaments %llu; coefficients mod %d of P vs fixed polynomial F:\n", p, n / p, [&]{ u64 s = 0; for (auto& kv : fixedPoly[p]) s += kv.second; return s; }(), p);
        bool ok = true; for (auto& kv : Pn) { u64 f = fixedPoly[p].count(kv.first) ? fixedPoly[p][kv.first] : 0; if (kv.second % p != f % p) ok = false; printf("  h=%llu: P=%llu (%llu mod %d), F=%llu (%llu mod %d)%s\n", kv.first, kv.second, kv.second % p, p, f, f % p, p, (kv.second % p == f % p) ? "" : "  MISMATCH"); }
        for (auto& kv : fixedPoly[p]) if (!Pn.count(kv.first)) { ok = false; printf("  h=%llu only in F  MISMATCH\n", kv.first); }
        printf("  congruence P_%d = F_{%d,%d} (mod %d): %s\n", n, n, p, p, ok ? "holds" : "FAILS"); }
    return 0; }
