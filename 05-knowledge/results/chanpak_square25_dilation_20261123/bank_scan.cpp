// Input is a complete gentourng parent bank, or marked lower cores for --lift.
// all_children and extend are linked from birth_cut_cover.cpp.
#include <array>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using Rows = std::array<int, 14>;
std::vector<std::int64_t> all_children(int n, const Rows &rows);
Rows extend(int n, Rows rows, int mask);

int scan(int n, int height) {
    constexpr int parent_counts[] = {0, 1, 1, 2, 4, 12, 56, 456, 6880, 191536};
    if (n < 5 || n > 9 || height < 1 || height % 2 == 0) return 2;
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
    return parents == parent_counts[n] ? 0 : 1;
}

struct LiftCounts {
    std::array<std::int64_t,14> expanded{};
    std::array<std::int64_t,14> hits{};
};

void lift(const Rows &low, const Rows &high, int n, int cap, int target, LiftCounts &counts) {
    ++counts.expanded[n];
    auto a = all_children(n, low);
    auto b = all_children(n, high);
    // The new vertex has identical contacts in both completions of the mark.
    for (int mask = 0; mask < (1 << n); ++mask) {
        if (a[mask] > target || b[mask] > target+2) continue;
        if (a[mask] == target && b[mask] == target+2) ++counts.hits[n+1];
        if (n+1 < cap)
            lift(extend(n,low,mask), extend(n,high,mask), n+1, cap, target, counts);
    }
}

int lift_cores(int cap, int target) {
    if (cap < 3 || cap > 13 || target < 1 || target % 2 == 0) return 2;
    std::array<LiftCounts,2> counts{};
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream input(line);
        std::vector<int> values;
        int value;
        while (input >> value) values.push_back(value);
        int n = int(values.size())-2;
        if (n < 2 || n >= cap) return 2;
        int u = values[n], v = values[n+1];
        if (u < 0 || u >= n || v < 0 || v >= n || u == v) return 2;
        Rows first{};
        for (int i = 0; i < n; ++i) first[i] = values[i];
        Rows second = first;
        second[u] ^= 1 << v;
        second[v] ^= 1 << u;
        lift(first,second,n,cap,target,counts[0]);
        lift(second,first,n,cap,target,counts[1]);
    }
    for (const auto &reading : counts) {
        for (int n = 0; n < cap; ++n) std::cout << reading.expanded[n] << ' ';
        for (int n = 0; n <= cap; ++n) std::cout << reading.hits[n] << ' ';
        std::cout << '\n';
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 3) return scan(std::stoi(argv[1]), std::stoi(argv[2]));
    if (argc == 4 && std::string(argv[1]) == "--lift")
        return lift_cores(std::stoi(argv[2]), std::stoi(argv[3]));
    return 2;
}
