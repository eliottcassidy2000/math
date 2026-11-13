#include <cstdint>
#include <iostream>
#include <vector>

struct Result { uint64_t h, q, end; };

Result count(const std::vector<unsigned>& rows, int root) {
    const int n = rows.size();
    const unsigned full = (1u << n) - 1;
    std::vector<uint64_t> dp((full + 1) * n, 0);
    for (int v = 0; v < n; ++v) dp[(1u << v) * n + v] = 1;
    for (unsigned mask = 1; mask <= full; ++mask) {
        unsigned used = mask;
        while (used) {
            const int last = __builtin_ctz(used); used &= used - 1;
            const uint64_t value = dp[mask * n + last];
            if (!value) continue;
            unsigned choices = rows[last] & (full ^ mask);
            while (choices) {
                const unsigned bit = choices & -choices; choices ^= bit;
                const int next = __builtin_ctz(bit);
                dp[(mask | bit) * n + next] += value;
            }
        }
    }
    uint64_t h = 0, q = 0;
    for (int v = 0; v < n; ++v) {
        h += dp[full * n + v];
        if (root >= 0) q += dp[(full ^ (1u << root)) * n + v];
    }
    return {h, q, root < 0 ? 0 : dp[full * n + root]};
}

int main() {
    int n, root;
    while (std::cin >> n >> root) {
        std::vector<unsigned> rows(n);
        for (auto& x : rows) std::cin >> x;
        auto first = count(rows, root);
        uint64_t s = 0;
        if (root >= 0) {
            const unsigned full = (1u << n) - 1;
            for (int i = 0; i < n; ++i) rows[i] = full ^ (1u << i) ^ rows[i];
            s = count(rows, root).end;
        }
        std::cout << first.h << ' ' << first.q << ' ' << s << '\n';
    }
}
