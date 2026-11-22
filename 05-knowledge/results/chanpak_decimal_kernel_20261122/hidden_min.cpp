// Compare complete slot profiles of single apex-contact flips through order nine.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

using Count = std::int64_t;
using Slots = std::array<Count, 4>;
using Endpoints = std::array<std::array<int, 8>, 8>;

static int first(int mask) { return __builtin_ctz(static_cast<unsigned>(mask)); }

std::vector<Slots> extensions(const std::vector<int>& rows) {
    const int n = rows.size(), full = (1 << n) - 1;
    std::vector<Endpoints> end(full + 1);
    std::vector<int> cycles(full + 1);
    std::vector<Slots> slot(full + 1);
    slot[0] = {1,1,1,1};
    for (int mask = 1; mask <= full; ++mask) {
        int root = first(mask), rest = mask ^ (1 << root);
        if (rest == 0) end[mask][root][root] = 1;
        for (int a = 0; a < n; ++a) {
            for (int b = 0; b < n; ++b) {
                int ways = end[mask][a][b];
                if (!ways) continue;
                for (int options = rows[b] & ~mask; options; options &= options-1) {
                    int c = first(options);
                    end[mask | (1 << c)][a][c] += ways;
                }
            }
        }
        if (rest && __builtin_popcount(static_cast<unsigned>(mask)) % 2)
            for (int b = 0; b < n; ++b)
                if (rows[b] & (1 << root)) cycles[mask] += end[mask][root][b];
        // The least vertex is uncovered, or lies in its unique packed odd cycle.
        slot[mask] = slot[rest];
        for (int sub = rest; sub; sub = (sub - 1) & rest) {
            int number = cycles[sub | (1 << root)];
            if (number)
                for (int j = 0; j < 4; ++j) slot[mask][j] += 2*number*slot[rest ^ sub][j];
        }
        for (int j = 0; j < 4; ++j) slot[mask][j] *= j+1;
    }

    std::array<std::array<Slots,8>,8> weight{};
    for (int mask = 1; mask <= full; ++mask) {
        if (__builtin_popcount(static_cast<unsigned>(mask)) % 2) continue;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                for (int j = 0; j < 4; ++j)
                    weight[a][b][j] += 2*(j+1)*end[mask][a][b]*slot[full ^ mask][j];
    }
    std::vector<Slots> child(full + 1);
    child[0] = slot[full];
    for (int j = 0; j < 4; ++j) child[0][j] *= j+1;
    for (int s = 1; s <= full; ++s) {
        int i = first(s), base = s ^ (1 << i);
        child[s] = child[base];
        for (int u = 0; u < n; ++u) {
            if (u == i) continue;
            for (int j = 0; j < 4; ++j)
                child[s][j] += (base & (1 << u)) ? -weight[u][i][j] : weight[i][u][j];
        }
    }
    return child;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    int n = std::stoi(argv[1]);
    if (n < 1 || n > 8) return 2;
    std::vector<int> row(n);
    int parents = 0;
    Count flags = 0, delayed = 0;
    while (std::cin >> row[0]) {
        for (int i = 1; i < n; ++i) if (!(std::cin >> row[i])) return 2;
        ++parents;
        auto child = extensions(row);
        for (int s = 0; s < (1 << n); ++s)
            for (int i = 0; i < n; ++i) {
                if (s & (1 << i)) continue;
                ++flags;
                const auto& a = child[s];
                const auto& b = child[s | (1 << i)];
                if (a[0] == b[0] && a[1] == b[1] && (a[2] != b[2] || a[3] != b[3])) ++delayed;
            }
    }
    constexpr int parent_counts[] = {0,1,1,2,4,12,56,456,6880};
    std::cout << "order " << n+1 << " parent-mask-contact choices " << flags
              << " delayed flips " << delayed << '\n';
    return parents == parent_counts[n] && delayed == 0 ? 0 : 1;
}
