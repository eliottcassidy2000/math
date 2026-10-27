// score_sum_descent.cpp -- mac-mini-2026-10-27-potentials
// For every strong class of order n: does some descent candidate Q (start in a minimum out-cut X, end outside X)
// have s(start)+s(end) != n-1 (hence black, since a blue word reverses outdegrees)?  Also: restricted to closed
// candidates (Hamiltonian-cycle cuts).  Failures are listed with their score sequences and regularity.
// Build: g++ -O2 -std=c++17 -o score_sum_descent score_sum_descent.cpp ; run: gentourng -cq n | ./score_sum_descent n
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
static int n;
struct Tour { uint16_t out[16]; };
int main(int argc, char **argv) {
    n = argc > 1 ? atoi(argv[1]) : 7;
    char line[256]; long long classes = 0, ok = 0, okClosed = 0, regular = 0, failNonreg = 0; vector<string> bad;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.out[i] |= (uint16_t)(1u << j); else T.out[j] |= (uint16_t)(1u << i); k++; }
        classes++;
        int s[16]; bool reg = true; for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(T.out[v]); if (2 * s[v] != n - 1) reg = false; }
        if (reg) regular++;
        int full = (1 << n) - 1, best = 1 << 20; vector<int> cuts;
        for (int X = 1; X < full; X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T.out[x] & ~X & full); if (f < best) { best = f; cuts.clear(); } if (f == best) cuts.push_back(X); }
        bool sep[16][16]; memset(sep, 0, sizeof sep);
        for (int X : cuts) for (int a = 0; a < n; a++) if (X >> a & 1) for (int b = 0; b < n; b++) if (!(X >> b & 1)) sep[a][b] = true;
        // reachable (start,end) pairs and closed pairs: DP over Hamiltonian paths by endpoint pairs
        // dp[mask][v] as bitset of starts
        vector<uint16_t> dp((size_t)1 << n << 4, 0);
        auto D = [&](int mask, int v) -> uint16_t & { return dp[((size_t)mask << 4) | v]; };
        for (int v = 0; v < n; v++) D(1 << v, v) = (uint16_t)(1u << v);
        for (int mask = 1; mask <= full; mask++) for (int v = 0; v < n; v++) { uint16_t st = D(mask, v); if (!st) continue; uint16_t o = T.out[v] & ~mask; while (o) { int w = __builtin_ctz(o); o &= o - 1; D(mask | (1 << w), w) |= st; } }
        bool found = false, foundClosed = false;
        for (int e = 0; e < n; e++) { uint16_t st = D(full, e); while (st) { int a = __builtin_ctz(st); st &= st - 1; if (sep[a][e] && s[a] + s[e] != n - 1) { found = true; if ((T.out[e] >> a) & 1) foundClosed = true; } } }
        if (found) ok++; if (foundClosed) okClosed++;
        if (!found && !reg) { failNonreg++; if (bad.size() < 40) { vector<int> sc(s, s + n); sort(sc.begin(), sc.end()); string ss; for (int x : sc) ss += to_string(x); bad.push_back(string(line, strlen(line) - 1) + " scores " + ss + " lambda " + to_string(best)); } }
    }
    printf("n=%d strong=%lld regular=%lld  score-sum certified descent=%lld (closed=%lld)  nonregular failures=%lld\n", n, classes, regular, ok, okClosed, failNonreg);
    for (auto &b : bad) printf("    %s\n", b.c_str());
    return 0;
}
