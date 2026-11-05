// hamcount.cpp -- independent subset-DP profile of a one-apex tournament given by its chain word (MSB first = v_1).
// Apex x = vertex 0; chain v_1..v_m = vertices 1..m, v_i -> v_j for i < j; x -> v_i iff bit = '1'.
// Prints H, E (Hamiltonian paths ending at x), s (starting at x), U (ordered covers (P1, x P2), P1 possibly empty),
// d = H(T - x).  Input lines: a b r n0 h D E U d word   (from recipe_test.py); checks the recorded (D,E,U,d) and the response.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <cstdint>
using namespace std;
typedef unsigned long long ull;
int main() {
    char line[4096]; int lines = 0, ok = 0, skipped = 0;
    while (fgets(line, sizeof line, stdin)) {
        long long a, b, r, n0, h, D, E, U, d; char word[2048];
        if (sscanf(line, "%lld %lld %lld %lld %lld %lld %lld %lld %lld %s", &a, &b, &r, &n0, &h, &D, &E, &U, &d, word) != 10) continue;
        lines++;
        int m = strlen(word); int n = m + 1;
        if (n > 20) { skipped++; continue; }
        vector<uint32_t> adj(n, 0);  // adj[u] bit v: u -> v
        for (int i = 1; i <= m; i++) { for (int j = i + 1; j <= m; j++) adj[i] |= 1u << j; if (word[i - 1] == '1') adj[0] |= 1u << i; else adj[i] |= 1u << 0; }
        size_t N = (size_t)1 << n;
        vector<ull> endc(N * n, 0), startc(N * n, 0);
        for (int v = 0; v < n; v++) { endc[((size_t)1 << v) * n + v] = 1; startc[((size_t)1 << v) * n + v] = 1; }
        for (size_t S = 1; S < N; S++) for (int v = 0; v < n; v++) {
            if (!((S >> v) & 1)) continue; ull e = endc[S * n + v], s = startc[S * n + v]; if (!e && !s) continue;
            uint32_t rest = (uint32_t)((N - 1) ^ S);
            if (e) { uint32_t out = adj[v] & rest; while (out) { int w = __builtin_ctz(out); out &= out - 1; endc[(S | ((size_t)1 << w)) * n + w] += e; } }
            if (s) { for (int w = 0; w < n; w++) if (((rest >> w) & 1) && ((adj[w] >> v) & 1)) startc[(S | ((size_t)1 << w)) * n + w] += s; }
        }
        size_t F = N - 1; ull H = 0; for (int v = 0; v < n; v++) H += endc[F * n + v];
        ull Ex = endc[F * n + 0], sx = startc[F * n + 0];
        // d = H(T - x): paths within the set F ^ 1
        ull dd = 0; for (int v = 1; v < n; v++) dd += endc[(F ^ 1) * n + v];
        // U = sum over X subset of chain (X excludes x) of H(T[X]) * start_x(T[V - X])
        ull UU = 0;
        for (size_t X = 0; X < N; X += 2) {  // even X: bit 0 (x) not in X
            ull hx = 0; if (X == 0) hx = 1; else for (int v = 1; v < n; v++) hx += endc[X * n + v];
            ull st = startc[(F ^ X) * n + 0];
            UU += hx * st;
        }
        bool good = (H == (ull)D) && (Ex == (ull)E) && (UU == (ull)U) && (dd == (ull)d) && (sx == 1) && (a * (long long)UU - b * (long long)H == r);
        if (good) ok++;
        printf("slope %lld/%lld seed n0=%lld h=%lld order %d: DP (H,E,U,d,s) = (%llu,%llu,%llu,%llu,%llu) recorded (%lld,%lld,%lld,%lld) response %lld expected %lld -> %s\n",
               a, b, n0, h, n, H, Ex, UU, dd, sx, D, E, U, d, a * (long long)UU - b * (long long)H, r, good ? "OK" : "MISMATCH");
    }
    printf("%d lines, %d verified physically, %d skipped (order > 20)\n", lines, ok, skipped);
    return 0;
}
