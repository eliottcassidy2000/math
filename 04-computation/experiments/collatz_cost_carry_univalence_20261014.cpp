// Exact finite certificate and independent frontier census; standard C++17 + GNU uint128.
// Build: g++ -O2 -std=c++17 -Wall -Wextra this_file.cpp -o /tmp/cost-carry
// Run:   /tmp/cost-carry > result.json
// Resource guard: length <= 20, at most 14,000,000 records in any frontier.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = std::uint64_t;
struct Cell { u64 key, B; };
constexpr u64 MASK = (u64(1) << 40) - 1;
u64 checks = 0;
void require(bool condition, const char* message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
unsigned log_floor(u64 n) { return 63 - __builtin_clzll(n); }
u64 pow3(unsigned L) { u64 q = 1; while (L--) q *= 3; return q; }
u64 multiply_mod(u64 a, u64 b, u64 q) {
    return static_cast<u64>((static_cast<__uint128_t>(a) * b) % q);
}
std::vector<unsigned> decode(unsigned L, unsigned A, u64 B) {
    std::vector<unsigned> word;
    for (u64 p = pow3(L - 1); L > 1; --L, p /= 3) {
        require(B > p, "decoder positive remainder");
        u64 remainder = B - p;
        unsigned a = __builtin_ctzll(remainder);
        require(a >= 1 && A > a, "decoder positive valuations");
        word.push_back(a);
        A -= a;
        B = remainder >> a;
    }
    require(B == 1 && A >= 1, "decoder terminal valuation");
    word.push_back(A);
    return word;
}
std::pair<unsigned, u64> affine(const std::vector<unsigned>& word, bool check_prefix) {
    unsigned A = 0;
    u64 B = 0, q = 1;
    for (unsigned a : word) {
        B = 3 * B + (u64(1) << A);
        A += a;
        q *= 3;
        if (check_prefix) require((u64(1) << A) < q, "native slope prefix");
    }
    return {A, B};
}
struct Bounds { u64 low, high; };
Bounds bounds(unsigned L, unsigned A) {
    u64 low = pow3(L) - (u64(1) << L), high = 0;
    for (unsigned i = 0; i < L; ++i) {
        unsigned prefix = std::min(log_floor(pow3(i)), A - L + i);
        high += pow3(L - 1 - i) * (u64(1) << prefix);
    }
    return {low, high};
}

int main() {
    constexpr unsigned maximum_length = 20;
    std::array<u64, 21> expected = {1,1,2,3,7,12,30,85,173,476,961,2652,
        8045,17637,51033,108950,312455,663535,1900470,5936673,13472296};
    // Tiny exact range certificate: every cost at every level, including extremizers.
    std::cout << "{\n  \"maximum_length\":20,\n  \"range_rows\":[\n";
    u64 previous_width = 0;
    for (unsigned L = 1; L <= maximum_length + 1; ++L) {
        u64 q = pow3(L);
        unsigned budget = log_floor(q);
        u64 largest_width = 0;
        for (unsigned A = L; A <= budget; ++A) {
            auto [low, high] = bounds(L, A);
            require(high >= low, "ordered range");
            std::vector<unsigned> minimum(L, 1), maximum;
            minimum.back() += A - L;
            unsigned previous = 0;
            for (unsigned i = 1; i <= L; ++i) {
                unsigned cumulative = std::min(log_floor(pow3(i)), A - L + i);
                require(cumulative > previous, "extremal positive increments");
                maximum.push_back(cumulative - previous);
                previous = cumulative;
            }
            require(affine(minimum, true) == std::make_pair(A, low), "minimum attained");
            require(affine(maximum, true) == std::make_pair(A, high), "maximum attained");
            largest_width = std::max(largest_width, high - low);
            if (L <= maximum_length) require(high - low < 4 * q, "mod-four range certificate");
        }
        require(largest_width == 3 * previous_width + (u64(1) << log_floor(pow3(L - 1)))
                - (u64(1) << (L - 1)), "normalized range monotonicity recurrence");
        previous_width = largest_width;
        if (L > 1) std::cout << ",\n";
        std::cout << "    {\"length\":" << L << ",\"budget\":" << budget
                  << ",\"maximum_width\":" << largest_width << ",\"four_times_modulus\":"
                  << 4 * q << ",\"range_proves_univalence\":"
                  << (largest_width < 4 * q ? "true" : "false") << "}";
    }
    std::cout << "\n  ],\n  \"census_rows\":[\n";
    std::vector<Cell> frontier{{0, 0}};
    std::array<u64, 64> cost_counts{};
    cost_counts[0] = 1;
    for (unsigned L = 1; L <= maximum_length; ++L) {
        u64 q = pow3(L);
        unsigned budget = log_floor(q);
        std::array<u64, 64> inverse{}, next_counts{}, observed_counts{};
        inverse[0] = 1;
        for (unsigned a = 1; a <= budget; ++a)
            inverse[a] = multiply_mod(inverse[a - 1], (q + 1) / 2, q);
        u64 count = 0;
        for (unsigned A = 0; A < budget; ++A)
            for (unsigned a = 1; a <= budget - A; ++a)
                next_counts[A + a] += cost_counts[A];
        for (u64 c : next_counts) count += c;
        require(count <= 14000000, "explicit resource guard");
        require(count == expected[L], "independent composition count");
        std::vector<Cell> next;
        next.reserve(count);
        for (const auto& cell : frontier) {
            unsigned A = cell.key >> 40;
            u64 B = 3 * cell.B + (u64(1) << A);
            u64 numerator = 3 * (cell.key & MASK) + 1;
            for (unsigned a = 1; a <= budget - A; ++a) {
                unsigned total = A + a;
                u64 carry = multiply_mod(numerator, inverse[a], q);
                require(carry == multiply_mod(B, inverse[total], q), "direct affine carry identity");
                require(B % 4 == (pow3(L - 1) + 2) % 4 || L == 1,
                        "common residue modulo four");
                next.push_back({(u64(total) << 40) | carry, B});
                ++observed_counts[total];
            }
        }
        require(next.size() == count && observed_counts == next_counts, "frontier cost counts");
        std::sort(next.begin(), next.end(), [](const Cell& a, const Cell& b) {
            return a.key < b.key || (a.key == b.key && a.B < b.B);
        });
        for (std::size_t i = 1; i < next.size(); ++i)
            require(next[i].key != next[i - 1].key, "same-cost carry collision");
        // Bounded inverse-decoder samples, independent of the carry transition.
        for (std::size_t i = 0; i < std::min<std::size_t>(16, next.size()); ++i) {
            unsigned A = next[i].key >> 40;
            auto word = decode(L, A, next[i].B);
            require(affine(word, true) == std::make_pair(A, next[i].B), "decoder round trip");
        }
        if (L > 1) std::cout << ",\n";
        std::cout << "    {\"length\":" << L << ",\"budget\":" << budget
                  << ",\"word_count\":" << count << ",\"distinct_cost_carry_pairs\":"
                  << next.size() << ",\"collision\":false}";
        frontier.swap(next);
        cost_counts = next_counts;
    }
    // A counterexample if prefix admissibility is replaced by terminal slope alone.
    std::vector<unsigned> w{1,3,1,1,1,1,3}, v{2,3,1,1,1,2,1};
    auto [Aw, Bw] = affine(w, false);
    auto [Av, Bv] = affine(v, false);
    require(Aw == 11 && Av == Aw && Bw == 4591 && Bv == 8965, "terminal-only witness affine data");
    require((u64(1) << Aw) < pow3(7), "terminal-only witness slope");
    require(Bv - Bw == 2 * pow3(7), "terminal-only witness equal carry");
    require(multiply_mod(911, u64(1) << Aw, pow3(7)) == Bw % pow3(7), "counterexample carry 911");
    require((u64(1) << (w[0] + w[1])) > pow3(2), "first witness fails prefix two");
    require((u64(1) << v[0]) > 3, "second witness fails prefix one");
    require(decode(7, Aw, Bw) == w && decode(7, Av, Bv) == v, "terminal-only witness decoder");
    std::cout << "\n  ],\n  \"terminal_only_counterexample\":{\"word1\":[1,3,1,1,1,1,3],"
              << "\"word2\":[2,3,1,1,1,2,1],\"length\":7,\"cost\":11,\"B1\":4591,"
              << "\"B2\":8965,\"carry\":911,\"modulus\":2187},\n"
              << "  \"checks\":" << checks << ",\n"
              << "  \"conclusion\":\"Cost-carry univalence is proved through length 20 by the exact mod-four range certificate and independently confirmed by exhaustive frontiers.\",\n"
              << "  \"limitations\":[\"No all-length univalence theorem.\","
              << "\"Failure of the length-21 range test is not a collision.\","
              << "\"No implication for convergence of a fixed positive integer source.\"]\n}\n";
}
