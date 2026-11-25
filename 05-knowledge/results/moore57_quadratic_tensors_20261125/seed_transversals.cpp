// Usage: seed_transversals [--words] < seed19.txt
// Input rows are the strict upper triangle of a cyclic 56-sheet chart.
#include <array>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Mask = std::uint64_t;
constexpr int sheets = 56;
constexpr Mask full = (Mask(1) << sheets) - 1;

Mask rotate(Mask values, int offset) {
    return ((values << offset) | (values >> (sheets - offset))) & full;
}

struct Search {
    int n;
    bool print_words;
    std::vector<std::vector<Mask>> forbidden, covered;
    std::vector<int> word;
    std::uint64_t solutions = 0;

    Search(const std::vector<std::vector<int>>& shifts, bool output)
        : n(static_cast<int>(shifts.size())), print_words(output),
          forbidden(n, std::vector<Mask>(n)), covered(forbidden), word(n) {
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            Mask used = Mask(1) << shifts[i][j];
            for (int k = 0; k < n; ++k) {
                if (k == i || k == j) continue;
                int diff = (shifts[i][k] + shifts[k][j]) % sheets;
                used |= Mask(1) << diff;
            }
            if (__builtin_popcountll(used) != n - 1)
                throw std::runtime_error("collision in the seed");
            forbidden[i][j] = used;
        }
    }

    void record() {
        ++solutions;
        for (int i = 0; i < n; ++i) for (int j = i+1; j < n; ++j)
            covered[i][j] |= Mask(1) << ((word[j] - word[i] + sheets) % sheets);
        if (print_words) {
            for (int i = 0; i < n; ++i) std::cout << (i ? " " : "") << word[i];
            std::cout << '\n';
        }
    }

    void extend(const std::array<Mask, sheets>& domain, Mask remaining) {
        if (!remaining) { record(); return; }
        int next = -1, size = sheets + 1;
        for (Mask bits = remaining; bits; bits &= bits - 1) {
            int i = __builtin_ctzll(bits);
            int count = __builtin_popcountll(domain[i]);
            if (count < size) { next = i; size = count; }
        }
        remaining ^= Mask(1) << next;
        for (Mask values = domain[next]; values; values &= values - 1) {
            int x = __builtin_ctzll(values);
            auto child = domain;
            bool possible = true;
            for (Mask bits = remaining; bits; bits &= bits - 1) {
                int i = __builtin_ctzll(bits);
                child[i] &= ~rotate(forbidden[next][i], x);
                if (!child[i]) { possible = false; break; }
            }
            if (possible) { word[next] = x; extend(child, remaining); }
        }
    }
};

int main(int argc, char** argv) {
    if (argc > 2 || (argc == 2 && std::string(argv[1]) != "--words")) return 2;
    std::vector<std::vector<int>> rows;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream input(line);
        std::vector<int> row;
        int value;
        while (input >> value) {
            if (value < 0 || value >= sheets) return 2;
            row.push_back(value);
        }
        if (!input.eof() || row.empty()) return 2;
        rows.push_back(row);
    }
    int n = static_cast<int>(rows.size()) + 1;
    if (n < 2 || n > sheets) return 2;
    std::vector<std::vector<int>> shifts(n, std::vector<int>(n));
    for (int i = 0; i+1 < n; ++i) {
        if (static_cast<int>(rows[i].size()) != n-i-1) return 2;
        for (int j = i+1; j < n; ++j) {
            shifts[i][j] = rows[i][j-i-1];
            shifts[j][i] = (sheets - shifts[i][j]) % sheets;
        }
    }
    Search search(shifts, argc == 2);
    std::array<Mask, sheets> domain;
    domain.fill(full);
    domain[0] = 1;
    search.extend(domain, (Mask(1) << n) - 1);
    int missed = 0;
    for (int i = 0; i < n; ++i) for (int j = i+1; j < n; ++j)
        missed += __builtin_popcountll(full & ~(search.covered[i][j] | search.forbidden[i][j]));
    std::cerr << n << " branches, " << search.solutions << " normalized words, "
              << missed << " unavailable differences\n";
}
