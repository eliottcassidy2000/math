// spectrum_poly_audit.cpp -- independent audit code (not the author's).
// For n <= 7: enumerate all labelled tournaments on {0..n-1}, count Hamiltonian paths by subset DP,
// build P_n(h) = #labelled tournaments with H = h, and for every prime p | n the polynomial F_{n,p}(h)
// over tournaments invariant under i -> i + n/p.  Also: circulant counts for odd prime n (every half-system S),
// and for n = 6, p = 3 the H distribution of the invariant tournaments.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <map>
#include <vector>
using namespace std;
typedef unsigned long long u64;
int n; int pid[8][8];
static inline int arcbit(u64 mask, int a, int b) { // 1 iff a -> b
    if (a < b) return (mask >> pid[a][b]) & 1; return 1 - ((mask >> pid[b][a]) & 1); }
u64 countH(const int adj[8]) { // adj[v] = bitmask of out-neighbours
    static u64 dp[256][8]; int full = (1 << n) - 1;
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { u64 c = dp[S][v]; if (!c) continue; if (!((S >> v) & 1)) continue;
        int nb = adj[v] & ~S; while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[S | (1 << w)][w] += c; } }
    u64 h = 0; for (int v = 0; v < n; v++) h += dp[full][v]; return h; }
int main(int argc, char** argv) {
    n = atoi(argv[1]); int P = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) pid[i][j] = P++;
    vector<int> primes; for (int p = 2; p <= n; p++) { bool pr = true; for (int d = 2; d * d <= p; d++) if (p % d == 0) pr = false; if (pr && n % p == 0) primes.push_back(p); }
    map<u64, u64> Pn; map<int, map<u64, u64>> F; map<u64, u64> circ; vector<pair<u64,u64>> circList;
    for (u64 mask = 0; mask < (1ULL << P); mask++) {
        int adj[8] = {0}; for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (a != b && arcbit(mask, a, b)) adj[a] |= 1 << b;
        u64 h = countH(adj); Pn[h]++;
        for (int p : primes) { int r = n / p; bool fixed = true;
            for (int a = 0; a < n && fixed; a++) for (int b = 0; b < n && fixed; b++) if (a != b && arcbit(mask, a, b) != arcbit(mask, (a + r) % n, (b + r) % n)) fixed = false;
            if (fixed) { F[p][h]++; if (p == n) { u64 S = 0; for (int d = 1; d < n; d++) if (arcbit(mask, 0, d)) S |= 1ULL << d; circList.push_back({S, h}); } } } }
    u64 tot = 0; for (auto& kv : Pn) tot += kv.second;
    printf("n = %d: %llu labelled tournaments. P_%d coefficients h:count", n, tot, n); for (auto& kv : Pn) printf(" %llu:%llu", kv.first, kv.second); printf("\n");
    for (int p : primes) { u64 nf = 0; for (auto& kv : F[p]) nf += kv.second; printf("p = %d (rotation by %d): %llu invariant tournaments; F coefficients h:count", p, n / p, nf); for (auto& kv : F[p]) printf(" %llu:%llu", kv.first, kv.second); printf("\n");
        bool ok = true; printf("  residues mod %d where P != 0 mod %d:", p, p); for (auto& kv : Pn) { u64 f = F[p].count(kv.first) ? F[p][kv.first] : 0; if (kv.second % p != f % p) ok = false; if (kv.second % p) printf(" h=%llu:P=%llu(%llu mod %d),F=%llu", kv.first, kv.second, kv.second % p, p, f); }
        for (auto& kv : F[p]) if (!Pn.count(kv.first)) ok = false; printf("\n  congruence P_%d == F_{%d,%d} mod %d: %s\n", n, n, p, p, ok ? "HOLDS" : "FAILS"); }
    if (!circList.empty()) { printf("circulants on Z_%d (S = out-neighbours of 0): ", n); for (auto& c : circList) { printf(" S={"); for (int d = 1; d < n; d++) if ((c.first >> d) & 1) printf("%d", d); printf("}:H=%llu", c.second); } printf("\n"); }
    return 0; }
