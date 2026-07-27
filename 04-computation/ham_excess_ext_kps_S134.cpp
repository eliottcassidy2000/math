// Extend HYP-9028 excess data: R_19, R_21, R_23 and QR_23.
// Bit-DP over subsets; memory 2^n * n * 8B (1.5 GB at n = 23).
#include <cstdio>
#include <cstdint>
#include <vector>
#include <cmath>
using namespace std;

static long double ham_excess(int n, const vector<uint32_t>& adj, unsigned long long* Hout) {
    size_t states = (size_t)1 << n;
    vector<unsigned long long> dp(states * (size_t)n, 0ULL);
    for (int i = 0; i < n; ++i) dp[((size_t)1 << i) * n + i] = 1;
    for (size_t mask = 1; mask < states; ++mask) {
        for (int last = 0; last < n; ++last) {
            unsigned long long c = dp[mask * (size_t)n + last];
            if (!c) continue;
            uint32_t nx = adj[last] & ~(uint32_t)mask;
            while (nx) {
                int v = __builtin_ctz(nx); nx &= nx - 1;
                dp[(mask | ((size_t)1 << v)) * (size_t)n + v] += c;
            }
        }
    }
    unsigned long long H = 0;
    size_t full = states - 1;
    for (int i = 0; i < n; ++i) H += dp[full * (size_t)n + i];
    *Hout = H;
    long double mean = 1.0L;
    for (int k = 2; k <= n; ++k) mean *= k;
    for (int k = 1; k < n; ++k) mean /= 2.0L;
    return (long double)H / mean;
}

int main() {
    // rotational
    for (int n : {19, 21, 23}) {
        vector<uint32_t> adj(n, 0);
        for (int i = 0; i < n; ++i)
            for (int d = 1; d <= (n - 1) / 2; ++d)
                adj[i] |= (uint32_t)1 << ((i + d) % n);
        unsigned long long H;
        long double e = ham_excess(n, adj, &H);
        printf("R_%d: H=%llu excess=%.6Lf n|H:%s\n", n, H, e,
               (H % (unsigned long long)n == 0) ? "yes" : "NO");
        fflush(stdout);
    }
    // Paley 23
    {
        int p = 23;
        bool qr[23] = {false};
        for (int i = 1; i < p; ++i) qr[(i * i) % p] = true;
        vector<uint32_t> adj(p, 0);
        for (int i = 0; i < p; ++i)
            for (int d = 1; d < p; ++d)
                if (qr[d]) adj[i] |= (uint32_t)1 << ((i + d) % p);
        unsigned long long H;
        long double e = ham_excess(p, adj, &H);
        unsigned long long aut = (unsigned long long)p * (p - 1) / 2;
        printf("QR_%d: H=%llu excess=%.6Lf tc=%llu p|H:%s\n", p, H, e,
               H % aut == 0 ? H / aut : 0ULL,
               (H % (unsigned long long)p == 0) ? "yes" : "NO");
    }
    printf("DONE\n");
    return 0;
}
