// Exact all-apex extension counts of a small parent tournament; standard C++17.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using Rows = std::array<int, 14>;
using Counts = std::array<std::int64_t, 14>;

struct Generator {
    std::uint64_t state = 20261122;
    std::uint64_t next() {
        auto z = (state += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }
};

Rows parent(int n, int percent, Generator &rng) {
    Rows rows{};
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            if (int(rng.next() % 100) < percent) rows[j] |= 1 << i;
            else rows[i] |= 1 << j;
        }
    return rows;
}

std::vector<Counts> endpoint_counts(int n, const Rows &incoming) {
    std::vector<Counts> dp(1 << n);
    for (int mask = 1; mask < (1 << n); ++mask)
        for (int bits = mask; bits; bits &= bits - 1) {
            int v = __builtin_ctz(bits), rest = mask ^ (1 << v);
            std::int64_t ways = (rest == 0);
            for (int prev = incoming[v] & rest; prev; prev &= prev - 1)
                ways += dp[rest][__builtin_ctz(prev)];
            dp[mask][v] = ways;
        }
    return dp;
}

std::vector<std::int64_t> all_children(int n, const Rows &out) {
    int full = (1 << n) - 1;
    Rows in{};
    for (int i = 0; i < n; ++i) in[i] = full ^ out[i] ^ (1 << i);
    auto last = endpoint_counts(n, in), first = endpoint_counts(n, out);
    std::array<Counts, 14> cross{};
    for (int mask = 1; mask < full; ++mask) {
        int rest = full ^ mask;
        for (int p = mask; p; p &= p - 1) {
            int i = __builtin_ctz(p);
            auto left = last[mask][i];
            if (left)
                for (int q = rest; q; q &= q - 1) {
                    int j = __builtin_ctz(q);
                    cross[i][j] += left * first[rest][j];
                }
        }
    }
    Counts linear{};
    for (int j = 0; j < n; ++j) {
        linear[j] = first[full][j] - last[full][j];
        for (int i = 0; i < n; ++i) linear[j] += cross[i][j];
    }
    std::vector<std::int64_t> values(1 << n);
    values[0] = std::accumulate(last[full].begin(), last[full].end(), std::int64_t{0});
    for (int mask = 1; mask <= full; ++mask) {
        int j = __builtin_ctz(mask), rest = mask ^ (1 << j);
        auto value = values[rest] + linear[j];
        for (int q = rest; q; q &= q - 1) {
            int i = __builtin_ctz(q);
            value -= cross[i][j] + cross[j][i];
        }
        values[mask] = value;
    }
    return values;
}

Rows extend(int n, Rows out, int mask) {
    for (int i = 0; i < n; ++i)
        if (!(mask >> i & 1)) out[i] |= 1 << n;
    out[n] = mask;
    return out;
}

void check_child(int n, const Rows &out, int mask, std::int64_t expected) {
    // Ending DP on the converse counts all paths of the original by reversal.
    auto dp = endpoint_counts(n + 1, extend(n, out, mask));
    auto actual = std::accumulate(dp.back().begin(), dp.back().end(), std::int64_t{0});
    if (actual != expected) throw std::logic_error("child count failed its direct subset DP");
}

int main(int argc, char **argv) {
    if (argc != 3) return 2;
    int shell = std::stoi(argv[1]), witness = std::stoi(argv[2]);
    if (shell != 14 && shell != 15) return 2;
    std::vector<bool> needed;
    int target = 0, left = 0;
    while (std::cin >> target) {
        if (target < 1) return 2;
        if (int(needed.size()) <= target) needed.resize(target + 1);
        if (!needed[target]) { needed[target] = true; ++left; }
    }
    for (int stage = 0; stage < (shell == 14 ? 1 : 2) && left; ++stage) {
        Generator rng;
        int budget = shell == 14 ? 950 : (stage == 0 ? 20000 : 1057);
        int trial = 0;
        for (; trial < budget && left; ++trial) {
            int n = shell == 14 ? 11 : (stage == 1 ? 13 : 12 + (trial % 3 == 0));
            int percent = 15 + 5 * (trial % 8);
            if (shell == 15) percent = stage == 1 ? 18 + trial % 10 :
                                      (n == 12 ? 10 + 5 * (trial % 8) : 8 + 2 * (trial % 11));
            auto out = parent(n, percent, rng);
            auto values = all_children(n, out);
            if (trial == 0)
                for (int mask : {0, 1, (1 << n) / 3, (1 << n) - 1})
                    check_child(n, out, mask, values[mask]);
            for (int mask = 0; mask < (1 << n); ++mask) {
                auto h = values[mask];
                if (h < 1 || h % 2 == 0) throw std::logic_error("invalid Hamiltonian count");
                if (h >= int(needed.size()) || !needed[h]) continue;
                needed[h] = false;
                --left;
                if (h == witness) {
                    check_child(n, out, mask, h);
                    std::cout << "Witness " << h << " on " << n + 1 << " vertices, outgoing masks:";
                    auto rows = extend(n, out, mask);
                    for (int v = 0; v <= n; ++v) std::cout << ' ' << rows[v];
                    std::cout << '\n';
                }
            }
        }
        std::cout << "shell " << shell << ", stage " << stage + 1 << ", parents " << trial
                  << ", targets remaining " << left << '\n';
    }
    return left ? 1 : 0;
}
