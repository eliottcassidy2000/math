// Exact native all-parity head for THM-4437.
// Reuse the audited six-sheet literal engine from THM-4434, but enumerate
// every parity mask and classify the three residual signed circuit types.
#define main canonical_odd_main
#include "lrc14_universal_literal_empty_core_sep06.cpp"
#undef main

#include <map>
#include <set>
#include <tuple>
#include <vector>

struct MixedLeader {
    I numerator = -1;
    I denominator = 1;
    V w{0, 0, 0};
    V projections{0, 0, 0};

    void update(I p, const Native& row, const V& speeds) {
        if (numerator < 0 || Wide(p) * denominator > Wide(numerator) * row.denominator) {
            numerator = p;
            denominator = row.denominator;
            w = speeds;
            projections = row.projections;
        }
    }
};

int main(int argc, char** argv) {
    const I H = argc > 1 ? std::stoll(argv[1]) : 199;
    std::map<int, I> rows, failures, equalities;
    std::map<int, I> low_signature_rows, low_signature_failures;
    std::map<int, MixedLeader> leaders;
    auto has_pattern = [](const V& w, V pattern) {
        std::sort(pattern.begin(), pattern.end());
        do {
            for (int signs = 0; signs < 8; ++signs) {
                I dot = 0;
                for (int i = 0; i < 3; ++i)
                    dot += ((signs >> i) & 1 ? -pattern[i] : pattern[i]) * w[i];
                if (dot == 0) return true;
            }
        } while (std::next_permutation(pattern.begin(), pattern.end()));
        return false;
    };
    I all_rows = 0, all_failures = 0, failures_outside_low = 0;
    I low_rows = 0, outside_low_any_gt = 0, outside_low_projection_equalities = 0;
    MixedLeader outside_low_leader, outside_low_individual_leader;
    std::vector<std::pair<V, int>> outside_low_equality_locus;
    std::uint64_t checksum = 14695981039346656037ULL;
    auto hash_integer = [&](I x) {
        auto value = static_cast<std::uint64_t>(x);
        for (int i = 0; i < 8; ++i) {
            checksum ^= value & 255U;
            checksum *= 1099511628211ULL;
            value >>= 8;
        }
    };
    V first_failure{0, 0, 0};
    for (I c = 3; c <= H; ++c) {
        if (c % 3 == 0) continue;
        for (I b = 2; b < c; ++b) {
            if (b % 3 == 0) continue;
            const I bc_gcd = std::gcd(b, c);
            for (I a = 1; a < b; ++a) {
                if (a % 3 == 0 || std::gcd(a, bc_gcd) != 1) continue;
                const V w{a, b, c};
                const int parity_mask = int((a & 1) | ((b & 1) << 1) | ((c & 1) << 2));
                const Native row = literal(w);
                const I minimum = *std::min_element(row.projections.begin(), row.projections.end());
                const bool p111 = has_pattern(w, V{1, 1, 1});
                const bool p112 = has_pattern(w, V{1, 1, 2});
                const bool p122 = has_pattern(w, V{1, 2, 2});
                const int low_bits = int(p111) | (int(p112) << 1) | (int(p122) << 2);
                const bool low = low_bits != 0;
                ++rows[parity_mask];
                ++all_rows;
                leaders[parity_mask].update(minimum, row, w);
                if (low) {
                    ++low_rows;
                    ++low_signature_rows[low_bits];
                } else {
                    outside_low_leader.update(minimum, row, w);
                    for (int coordinate = 0; coordinate < 3; ++coordinate) {
                        const I projection = row.projections[coordinate];
                        outside_low_individual_leader.update(projection, row, w);
                        const Wide projection_comparison =
                            Wide(77) * projection - Wide(6) * row.denominator;
                        if (projection_comparison > 0) ++outside_low_any_gt;
                        if (projection_comparison == 0) {
                            ++outside_low_projection_equalities;
                            outside_low_equality_locus.push_back({w, coordinate});
                        }
                    }
                }
                for (I value : w) hash_integer(value);
                hash_integer(row.denominator);
                for (I value : row.projections) hash_integer(value);
                hash_integer(low_bits);
                const Wide comparison = Wide(77) * minimum - Wide(6) * row.denominator;
                if (comparison > 0) {
                    ++failures[parity_mask];
                    ++all_failures;
                    if (low) ++low_signature_failures[low_bits];
                    else ++failures_outside_low;
                    if (first_failure[2] == 0) {
                        first_failure = w;
                        std::cout << "FIRST_FAILURE w=" << triple(w)
                                  << " E=" << fraction(row.projections[0], row.denominator) << ","
                                  << fraction(row.projections[1], row.denominator) << ","
                                  << fraction(row.projections[2], row.denominator)
                                  << " mass=" << fraction(row.mass, row.denominator)
                                  << " contacts=" << row.contacts << "\n";
                    }
                } else if (comparison == 0) {
                    ++equalities[parity_mask];
                }
            }
        }
    }
    std::cout << "MIXED_PARITY_HEAD H=" << H << " rows=" << all_rows
              << " min_failures=" << all_failures << " outside_low=" << failures_outside_low
              << " low_rows=" << low_rows << " outside_low_any_E_gt=" << outside_low_any_gt
              << " outside_low_projection_equalities=" << outside_low_projection_equalities
              << " first_failure=" << triple(first_failure) << "\n";
    std::cout << "OUTSIDE_LOW max_min=" << fraction(outside_low_leader.numerator, outside_low_leader.denominator)
              << " at=" << triple(outside_low_leader.w)
              << " max_individual=" << fraction(outside_low_individual_leader.numerator, outside_low_individual_leader.denominator)
              << " at=" << triple(outside_low_individual_leader.w) << "\n";
    for (const auto& entry : outside_low_equality_locus)
        std::cout << "OUTSIDE_LOW_EQUALITY w=" << triple(entry.first)
                  << " coordinate=" << entry.second << "\n";
    for (const auto& item : low_signature_rows)
        std::cout << "low_bits=" << item.first << " rows=" << item.second
                  << " min_failures=" << low_signature_failures[item.first] << "\n";
    for (const auto& item : rows) {
        const int mask = item.first;
        const auto& lead = leaders[mask];
        std::cout << "mask=" << mask << " rows=" << item.second
                  << " failures=" << failures[mask] << " equalities=" << equalities[mask]
                  << " max_min=" << fraction(lead.numerator, lead.denominator)
                  << " at=" << triple(lead.w)
                  << " E=" << fraction(lead.projections[0], lead.denominator) << ","
                  << fraction(lead.projections[1], lead.denominator) << ","
                  << fraction(lead.projections[2], lead.denominator) << "\n";
    }
    std::cout << "SEMANTIC_FNV1A64 " << checksum << "\n";
    need(failures_outside_low == 0, "nonlow minimum failure");
    need(outside_low_any_gt == 0, "nonlow individual projection failure");
    need(first_failure == V{2, 5, 7}, "first low-circuit hostile");
    if (H == 199) need(all_rows == 333960, "H199 row count");
    if (H == 611) {
        need(all_rows == 9720930, "H611 row count");
        need(all_failures == 14220, "H611 low-circuit failure count");
        const std::vector<std::pair<V, int>> expected{
            {{7, 16, 22}, 1}, {{14, 17, 22}, 2}, {{4, 19, 22}, 2}
        };
        need(outside_low_equality_locus == expected, "individual equality locus");
    }
    std::cout << "PASS: outside the three low signed circuits every projection is at most 6/77 and the minimum is strict.\n";
    return 0;
}
