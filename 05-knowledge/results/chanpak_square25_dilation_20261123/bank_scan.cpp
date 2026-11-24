// Uses the inherited all-apex height table from birth_cut_cover.cpp.
// bank.py --census builds this reader and streams the complete gentourng bank.
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using Rows = std::array<int, 14>;
std::vector<std::int64_t> all_children(int n, const Rows &rows);
Rows extend(int n, Rows rows, int mask);

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    int n = std::stoi(argv[1]), height = std::stoi(argv[2]);
    if (!((n == 5 && height == 23) || (n == 9 && height == 12767))) return 2;
    int parents = 0;
    std::string record;
    while (std::cin >> record) {
        if (record.size() != static_cast<std::size_t>(n*(n-1)/2)
            || record.find_first_not_of("01") != std::string::npos) return 1;
        Rows rows{};
        int offset = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                if (record[offset++] == '1') rows[i] |= 1 << j;
                else rows[j] |= 1 << i;
            }
        }
        ++parents;
        auto heights = all_children(n, rows);
        for (int mask = 0; mask < (1 << n); ++mask) {
            if (heights[mask] != height) continue;
            for (int i = 0; i < n; ++i) {
                if (heights[mask ^ (1 << i)] != height+2) continue;
                auto lower = extend(n, rows, mask);
                for (int v = 0; v <= n; ++v) std::cout << lower[v] << ' ';
                std::cout << i << ' ' << n << '\n';
            }
        }
    }
    return parents == (n == 5 ? 12 : 191536) ? 0 : 1;
}
