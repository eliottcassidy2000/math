// Strong-H spectrum at n=8: 2^28 labeled tournaments.
#include <cstdio>
#include <cstdint>
#include <set>
using namespace std;
const int N = 8;
long long ham(const uint8_t adj[8]) {
    static long long dp[256][8];
    int full = 255;
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
bool strong(const uint8_t adj[8]) {
    int reach = 1, prev = 0;
    while (reach != prev) { prev = reach;
        for (int i = 0; i < N; ++i) if ((reach >> i) & 1) reach |= adj[i]; }
    if (reach != 255) return false;
    uint8_t r[8] = {0};
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j)
        if (i != j && ((adj[i] >> j) & 1)) r[j] |= 1 << i;
    reach = 1; prev = 0;
    while (reach != prev) { prev = reach;
        for (int i = 0; i < N; ++i) if ((reach >> i) & 1) reach |= r[i]; }
    return reach == 255;
}
int main() {
    set<long long> vals;
    for (long long code = 0; code < (1LL << 28); ++code) {
        uint8_t adj[8] = {0};
        int bit = 0;
        for (int i = 0; i < N; ++i)
            for (int j = i + 1; j < N; ++j, ++bit) {
                if ((code >> bit) & 1) adj[i] |= 1 << j;
                else adj[j] |= 1 << i;
            }
        if (!strong(adj)) continue;
        vals.insert(ham(adj));
        if ((code & 0xFFFFFF) == 0xFFFFFF) fprintf(stderr, "%lld/268435456 vals=%zu\n", code+1, vals.size());
    }
    printf("n=8 strong-H: %zu values, min=%lld max=%lld\n", vals.size(),
           *vals.begin(), *vals.rbegin());
    printf("values:");
    for (long long v : vals) printf(" %lld", v);
    printf("\n");
    return 0;
}
