// Strong-tournament H spectrum for n = 3..7: enumerate all labeled
// tournaments, filter strongly connected, collect the SET of H values.
#include <cstdio>
#include <cstdint>
#include <set>
#include <vector>
using namespace std;
int N;
long long ham(const uint8_t adj[8]) {
    int full = (1 << N) - 1;
    static long long dp[128][7];
    for (int m = 0; m <= full; ++m) for (int i = 0; i < N; ++i) dp[m][i] = 0;
    for (int i = 0; i < N; ++i) dp[1 << i][i] = 1;
    for (int m = 1; m <= full; ++m)
        for (int l = 0; l < N; ++l) {
            long long c = dp[m][l];
            if (!c) continue;
            int nx = adj[l] & ~m;
            while (nx) { int v = __builtin_ctz(nx); nx &= nx - 1; dp[m | (1 << v)][v] += c; }
        }
    long long t = 0;
    for (int i = 0; i < N; ++i) t += dp[full][i];
    return t;
}
bool strongly_connected(const uint8_t adj[8]) {
    // reachability from 0 and to 0
    int reach = 1, prev = 0;
    while (reach != prev) { prev = reach;
        for (int i = 0; i < N; ++i) if ((reach >> i) & 1) reach |= adj[i]; }
    if (reach != (1 << N) - 1) return false;
    // reverse
    uint8_t radj[8] = {0};
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j)
        if (i != j && ((adj[i] >> j) & 1)) radj[j] |= 1 << i;
    reach = 1; prev = 0;
    while (reach != prev) { prev = reach;
        for (int i = 0; i < N; ++i) if ((reach >> i) & 1) reach |= radj[i]; }
    return reach == (1 << N) - 1;
}
int main() {
    for (N = 3; N <= 7; ++N) {
        int m = N * (N - 1) / 2;
        set<long long> vals;
        for (long long code = 0; code < (1LL << m); ++code) {
            uint8_t adj[8] = {0};
            int bit = 0;
            for (int i = 0; i < N; ++i)
                for (int j = i + 1; j < N; ++j, ++bit) {
                    if ((code >> bit) & 1) adj[i] |= 1 << j;
                    else adj[j] |= 1 << i;
                }
            if (!strongly_connected(adj)) continue;
            vals.insert(ham(adj));
        }
        printf("n=%d strong-H values (%zu):", N, vals.size());
        for (long long v : vals) printf(" %lld", v);
        printf("\n");
        fflush(stdout);
    }
    return 0;
}
