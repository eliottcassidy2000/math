// Literal central-layer Paley extensions. Only orders <=20 are used, so
// every intermediate path count is bounded by 20! and fits in uint64_t.
#include <cstdint>
#include <iostream>
#include <vector>

std::uint64_t central_count(int prime) {
    const int n = prime + 1;
    const int size = 1 << n;
    std::vector<int> incoming(n);
    for (int u = 0; u < prime; ++u) {
        for (int j = 1; j <= (prime - 1) / 2; ++j) {
            const int v = (u + j*j) % prime;
            incoming[v] |= 1 << u;
        }
        if (u < (prime - 1) / 2) {
            incoming[u] |= 1 << prime;
        } else {
            incoming[prime] |= 1 << u;
        }
    }
    std::vector<std::uint64_t> paths(static_cast<std::size_t>(n) * size);
    for (int i = 0; i < n; ++i) {
        paths[n * (1 << i) + i] = 1;
    }
    for (int mask = 1; mask < size; ++mask) {
        for (int ends = mask; ends; ends &= ends - 1) {
            const int end = __builtin_ctz(ends);
            const int rest = mask ^ (1 << end);
            for (int pred = rest & incoming[end]; pred; pred &= pred - 1) {
                paths[n * mask + end] += paths[n * rest + __builtin_ctz(pred)];
            }
        }
    }
    std::uint64_t total = 0;
    for (int i = 0; i < n; ++i) {
        total += paths[n * (size - 1) + i];
    }
    return total;
}

int main() {
    for (int p : {3, 7, 11, 19}) {
        std::cout << p << " " << central_count(p) << "\n";
    }
}
