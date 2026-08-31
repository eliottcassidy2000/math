#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Pair = std::pair<int, int>;
using Rows = std::set<Pair>;

Rows read_pairs(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);
    Rows out;
    std::string line;
    Pair prior{-1, -1};
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line == "q,r") continue;
        const auto comma = line.find(',');
        if (comma == std::string::npos) throw std::runtime_error("bad row");
        Pair row{std::stoi(line.substr(0, comma)),
                 std::stoi(line.substr(comma + 1))};
        if (!(prior < row)) throw std::runtime_error("not strictly sorted");
        if (!(0 < row.first && row.first < row.second))
            throw std::runtime_error("bad ordered pair");
        prior = row;
        out.insert(row);
    }
    return out;
}

Rows unite(const Rows& a, const Rows& b) {
    Rows out;
    std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                   std::inserter(out, out.end()));
    return out;
}

Rows intersect(const Rows& a, const Rows& b) {
    Rows out;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                          std::inserter(out, out.end()));
    return out;
}

Rows difference(const Rows& a, const Rows& b) {
    Rows out;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(out, out.end()));
    return out;
}

bool proper_subset(const Rows& a, const Rows& b) {
    return a.size() < b.size() &&
           std::includes(b.begin(), b.end(), a.begin(), a.end());
}

std::uint64_t fnv(const Rows& rows) {
    std::uint64_t state = 0xcbf29ce484222325ULL;
    for (auto [q, r] : rows) {
        for (std::uint64_t value : {static_cast<std::uint64_t>(q),
                                    static_cast<std::uint64_t>(r)}) {
            for (unsigned shift = 0; shift < 64; shift += 8) {
                state ^= (value >> shift) & 0xffULL;
                state *= 0x100000001b3ULL;
            }
        }
    }
    return state;
}

void describe(const std::string& label, const Rows& rows) {
    std::cout << label << " rows=" << std::dec << rows.size()
              << " fnv=" << std::hex << std::setfill('0') << std::setw(16)
              << fnv(rows) << '\n';
}

void write_rows(const std::string& path, const Rows& rows) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + path);
    for (auto [q, r] : rows) out << q << ',' << r << '\n';
}

int main(int argc, char** argv) {
    if (argc != 12) {
        std::cerr << "usage: audit LIVE S J E C H19 H294 H372 UNION RESIDUAL NEW\n";
        return 2;
    }
    const Rows live = read_pairs(argv[1]);
    std::vector<std::pair<std::string, Rows>> nodes;
    const char* names[] = {"S", "J", "E", "C", "H19", "H294", "H372"};
    for (int i = 0; i < 7; ++i)
        nodes.push_back({names[i], read_pairs(argv[i + 2])});

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        for (std::size_t j = i + 1; j < nodes.size(); ++j) {
            describe(nodes[i].first + "&" + nodes[j].first,
                     intersect(nodes[i].second, nodes[j].second));
        }
    }

    Rows ours;
    for (int i = 0; i < 3; ++i) ours = unite(ours, nodes[i].second);
    Rows incoming;
    for (int i = 3; i < 7; ++i) incoming = unite(incoming, nodes[i].second);
    const Rows covered = intersect(ours, incoming);
    const Rows added = difference(incoming, ours);
    const Rows all = unite(ours, incoming);
    const Rows complement = difference(live, all);

    describe("OURS", ours);
    describe("INCOMING", incoming);
    describe("COVERED", covered);
    describe("ADDED", added);
    describe("UNION", all);
    describe("COMPLEMENT", complement);

    const Rows expected_added{{147, 294}, {147, 590}, {372, 619}};
    const Rows expected_top{{100, 626}, {256, 626}};
    Rows top;
    int maximum = -1;
    for (auto [q, r] : complement) maximum = std::max(maximum, r);
    for (auto row : complement)
        if (row.second == maximum) top.insert(row);
    if (!proper_subset(nodes[3].second, nodes[2].second))
        throw std::runtime_error("C is not a proper subset of E");
    if (nodes[4].second != nodes[1].second)
        throw std::runtime_error("H19 differs from J");
    if (added != expected_added || all.size() != 1324 ||
        fnv(all) != 0xf55ee025df29bb65ULL || complement.size() != 21323 ||
        fnv(complement) != 0x09a0dfc4515d556bULL || maximum != 626 ||
        top != expected_top || unite(all, complement) != live ||
        !intersect(all, complement).empty()) {
        throw std::runtime_error("cross-packet identity check failed");
    }

    write_rows(argv[9], all);
    write_rows(argv[10], complement);
    write_rows(argv[11], added);
    std::cout << "MAX " << std::dec << maximum
              << " TOP (100,626) (256,626)\n";
    std::cout << "TYPE_GUARD typed_row_set_join_only_no_common_deck_no_LRC14\n";
}
