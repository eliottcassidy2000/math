// Independent audit engine for THM-4864 (written from scratch by the referee).
//
// Two modes, selected by the first token of each job on stdin:
//
//   H  N  out[0..N-1]  K  (a c)*K
//      Hamilton-path count of the augmented tournament on N+2 vertices:
//      old names 0..N-1 with out-neighbour bitmasks out[i];
//      x = N, y = N+1, arc x->y; v->x iff bit v of a, else x->v;
//      v->y iff bit v of c, else y->v.
//      Prints K exact unsigned 64-bit counts (one line).
//
//   C  n  out[0..n-1]  wmask  K  (a c)*K
//      Parity of the number of directed cycles of the same augmented
//      tournament that contain at least one port and have odd total
//      weight; old vertex v weighs bit v of wmask, each port weighs 1.
//      Brute force by subset DP from each port (no use of any path formula).
//      Prints K bits (one line).
//
// Input is terminated by EOF.
#include <cstdio>
#include <cstdint>
#include <vector>
#include <cstring>
using namespace std;
typedef unsigned long long u64;

static vector<u64> dp;

static u64 hamcount(int V, const vector<unsigned>& out) {
    size_t full = (1u << V);
    dp.assign(full * V, 0ULL);
    for (int v = 0; v < V; ++v) dp[((size_t)1 << v) * V + v] = 1;
    for (size_t S = 1; S < full; ++S) {
        for (int v = 0; v < V; ++v) {
            if (!(S >> v & 1)) continue;
            u64 val = dp[S * V + v];
            if (!val) continue;
            unsigned nb = out[v] & ~(unsigned)S;
            while (nb) {
                int w = __builtin_ctz(nb);
                nb &= nb - 1;
                dp[(S | (1u << w)) * V + w] += val;
            }
        }
    }
    u64 tot = 0;
    for (int v = 0; v < V; ++v) tot += dp[(full - 1) * V + v];
    return tot;
}

// number (mod 2) of directed cycles through `start` using only vertices in
// `allowed` (bitmask, contains start) whose vertex-weight parity is odd.
static int cycpar(int V, const vector<unsigned>& out, unsigned wmask_aug, int start, unsigned allowed) {
    size_t full = (1u << V);
    // f[S][v] parity of paths start -> ... -> v with vertex set S
    vector<unsigned char> f(full * V, 0);
    f[((size_t)1 << start) * V + start] = 1;
    int total = 0;
    for (size_t S = 1; S < full; ++S) {
        if (!(S >> start & 1)) continue;
        if ((S & ~(size_t)allowed) != 0) continue;
        int wpar = __builtin_popcount((unsigned)S & wmask_aug) & 1;
        for (int v = 0; v < V; ++v) {
            if (!(S >> v & 1)) continue;
            unsigned char val = f[S * V + v];
            if (!val) continue;
            if (v != start && (out[v] >> start & 1) && wpar) total ^= 1;  // close cycle v->start
            unsigned nb = out[v] & ~(unsigned)S & allowed;
            while (nb) {
                int w = __builtin_ctz(nb);
                nb &= nb - 1;
                f[(S | (1u << w)) * V + w] ^= val;
            }
        }
    }
    return total;
}

int main() {
    char mode[4];
    while (scanf("%3s", mode) == 1) {
        int n;
        if (scanf("%d", &n) != 1) return 1;
        vector<unsigned> base(n);
        for (int i = 0; i < n; ++i) scanf("%u", &base[i]);
        unsigned wmask = 0;
        if (mode[0] == 'C') scanf("%u", &wmask);
        int K;
        scanf("%d", &K);
        int V = n + 2, x = n, y = n + 1;
        for (int q = 0; q < K; ++q) {
            unsigned a, c;
            scanf("%u %u", &a, &c);
            vector<unsigned> out(V, 0);
            for (int v = 0; v < n; ++v) {
                out[v] = base[v];
                if (a >> v & 1) out[v] |= 1u << x;
                if (c >> v & 1) out[v] |= 1u << y;
            }
            out[x] = 1u << y;
            for (int v = 0; v < n; ++v) if (!(a >> v & 1)) out[x] |= 1u << v;
            out[y] = 0;
            for (int v = 0; v < n; ++v) if (!(c >> v & 1)) out[y] |= 1u << v;
            if (mode[0] == 'H') {
                printf("%llu%c", hamcount(V, out), q + 1 == K ? '\n' : ' ');
            } else {
                unsigned waug = wmask | (1u << x) | (1u << y);
                unsigned all = (1u << V) - 1;
                int g = cycpar(V, out, waug, x, all) ^ cycpar(V, out, waug, y, all & ~(1u << x));
                printf("%d%c", g, q + 1 == K ? '\n' : ' ');
            }
        }
        if (K == 0) printf("\n");
        fflush(stdout);
    }
    return 0;
}
