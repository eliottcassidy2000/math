#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Primary exact audit for THM-4190.
//
// THM-4188 transports the q=50 repair layers outside its exact 23-resonance
// set.  At those 23 labels this program audits the theorem's consequence
// directly.  It integrates each newcomer comb inside the fixed pool-wall
// cells, projects labelled FAILURE masks to the 25 optional coordinates of
// each primitive triple anchor, performs a dense subset zeta, and evaluates
// every one of the C(25,7)=480,700 labelled complements.

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr i64 COMMON = INT64_C(18241159416480);
constexpr std::array<int, 23> RESONANCES = {
    6, 22, 24, 25, 48, 70, 72, 96, 100, 105, 110, 128,
    130, 140, 186, 192, 206, 210, 220, 256, 260, 294, 366};
constexpr u64 FNV_OFFSET = UINT64_C(0xcbf29ce484222325);
constexpr u64 FNV_PRIME = UINT64_C(0x100000001b3);
constexpr u64 EXPECTED_LEDGER = UINT64_C(0x72db7e3090d5626e);

constexpr u64 binomial(int n, int k) {
    u64 result = 1;
    for (int index = 1; index <= k; ++index) {
        result = result * static_cast<u64>(n - k + index) /
                 static_cast<u64>(index);
    }
    return result;
}

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

i64 checked_lcm(i64 left, i64 right) {
    const i64 divisor = std::gcd(left, right);
    const i128 value = static_cast<i128>(left / divisor) * right;
    require(value <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(value);
}

u32 next_combination(u32 mask) {
    const u32 low = mask & (~mask + 1u);
    const u32 ripple = mask + low;
    if (ripple == 0) return 0;
    return ripple | (((mask ^ ripple) >> 2) / low);
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    if (negative) value = -value;
    std::string result;
    while (value != 0) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

std::string hex64(u64 value) {
    std::ostringstream out;
    out << std::hex << std::nouppercase << std::setw(16) << std::setfill('0')
        << value;
    return out.str();
}

class Ledger {
  public:
    void add(u64 word) {
        for (int shift = 0; shift < 64; shift += 8) {
            value_ ^= static_cast<std::uint8_t>((word >> shift) & 0xffu);
            value_ *= FNV_PRIME;
        }
    }
    u64 value() const { return value_; }

  private:
    u64 value_ = FNV_OFFSET;
};

u32 label_mask(std::initializer_list<int> labels) {
    u32 result = 0;
    for (int label : labels) {
        const auto found = std::find(POOL.begin(), POOL.end(), label);
        require(found != POOL.end(), "label outside pool");
        result |= u32{1} << static_cast<unsigned>(found - POOL.begin());
    }
    return result;
}

std::string labels(u32 mask) {
    std::string result;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((mask & (u32{1} << vertex)) == 0) continue;
        if (!result.empty()) result += ',';
        result += std::to_string(POOL[vertex]);
    }
    return result.empty() ? "-" : result;
}

bool midpoint_safe(int speed, i64 left, i64 right) {
    const i64 residue = static_cast<i64>(
        (static_cast<i128>(speed) * (left + right)) %
        (static_cast<i128>(2) * COMMON));
    return static_cast<i128>(7) * residue >= COMMON &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * COMMON;
}

struct Cell {
    i64 left;
    i64 right;
    u32 failed_pool;
};

std::vector<Cell> build_pool_cells() {
    i64 denominator = 1;
    for (int speed : POOL) denominator = checked_lcm(denominator, 14LL * speed);
    require(denominator == COMMON, "fixed-pool denominator changed");

    std::vector<i64> walls = {0, COMMON};
    for (int speed : POOL) {
        const i64 unit = COMMON / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.size() == 7134, "fixed-pool wall count changed");

    std::vector<Cell> cells;
    cells.reserve(walls.size() - 1);
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        u32 failure = 0;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if (!midpoint_safe(POOL[vertex], walls[index], walls[index + 1])) {
                failure |= u32{1} << vertex;
            }
        }
        cells.push_back({walls[index], walls[index + 1], failure});
    }
    return cells;
}

// Integral of the q-safe comb on [0,tick/COMMON], normalized so that a full
// circle has denominator 14*q*COMMON.
i64 safe_prefix(i64 tick, int q) {
    const i128 product = static_cast<i128>(q) * tick;
    const i64 whole = static_cast<i64>(product / COMMON);
    const i64 remainder = static_cast<i64>(product % COMMON);
    const i64 scaled = 14 * remainder;
    i64 partial = 0;
    if (scaled <= COMMON) {
        partial = 0;
    } else if (scaled >= 13 * COMMON) {
        partial = 12 * COMMON;
    } else {
        partial = scaled - COMMON;
    }
    return static_cast<i64>(static_cast<i128>(12) * whole * COMMON + partial);
}

struct BodyAudit {
    u64 presentations = 0;
    u64 below = 0;
    u64 equalities = 0;
    i64 minimum_numerator = 0;
    i128 minimum_delta = 0;
    u32 minimum_k = 0;
};

BodyAudit audit_bodies(const std::vector<Cell>& cells,
                       int q,
                       u32 anchor,
                       u32 allowed) {
    require(std::popcount(anchor) == 3, "anchor is not a triple");
    require(std::popcount(allowed) == 25, "allowed ground is not 25");
    std::array<int, 30> local_index{};
    local_index.fill(-1);
    std::array<int, 25> local_to_pool{};
    int local_count = 0;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((allowed & (u32{1} << vertex)) != 0) {
            local_index[vertex] = local_count;
            local_to_pool[local_count++] = vertex;
        }
    }
    require(local_count == 25, "local coordinate count changed");

    constexpr u32 local_size = u32{1} << 25;
    std::vector<i64> subset_mass(local_size, 0);
    i64 previous = safe_prefix(0, q);
    i128 total_q_safe = 0;
    for (const Cell& cell : cells) {
        const i64 current = safe_prefix(cell.right, q);
        const i64 contribution = current - previous;
        previous = current;
        require(contribution >= 0, "negative q-safe cell contribution");
        total_q_safe += contribution;
        if ((cell.failed_pool & anchor) != 0) continue;
        u32 local_failure = 0;
        u32 failure = cell.failed_pool & allowed;
        while (failure != 0) {
            const u32 bit = failure & (~failure + 1u);
            local_failure |= u32{1} << local_index[std::countr_zero(bit)];
            failure ^= bit;
        }
        subset_mass[local_failure] += contribution;
    }
    require(total_q_safe == static_cast<i128>(12) * q * COMMON,
            "newcomer safe mass is not 6/7");

    // After the subset zeta, subset_mass[S] is the weight of all anchor-safe
    // cells whose optional failure mask is contained in S.
    for (int bit = 0; bit < 25; ++bit) {
        const u32 half = u32{1} << bit;
        const u32 block = half << 1;
        for (u32 base = 0; base < local_size; base += block) {
            for (u32 offset = 0; offset < half; ++offset) {
                subset_mass[base + half + offset] += subset_mass[base + offset];
            }
        }
    }

    BodyAudit result;
    bool first = true;
    const u32 all_local = local_size - 1;
    u32 local_k = (u32{1} << 7) - 1;
    while (local_k != 0 && local_k < local_size) {
        const i64 numerator = subset_mass[all_local ^ local_k];
        const i128 delta = static_cast<i128>(9) * numerator -
                           static_cast<i128>(8) * q * COMMON;
        u32 global_k = 0;
        u32 bits = local_k;
        while (bits != 0) {
            const u32 bit = bits & (~bits + 1u);
            global_k |= u32{1} << local_to_pool[std::countr_zero(bit)];
            bits ^= bit;
        }
        ++result.presentations;
        if (delta < 0) ++result.below;
        if (delta == 0) ++result.equalities;
        if (first || delta < result.minimum_delta) {
            result.minimum_numerator = numerator;
            result.minimum_delta = delta;
            result.minimum_k = global_k;
            first = false;
        }
        const u32 next = next_combination(local_k);
        if (next <= local_k) break;
        local_k = next;
    }
    require(result.presentations == 480700, "presentation count changed");
    return result;
}

struct AnchorRow {
    const char* name;
    u32 anchor;
};

}  // namespace

int main() {
    try {
        std::cout.setf(std::ios::unitbuf);
        const std::vector<Cell> cells = build_pool_cells();
        const u32 full_pool = (u32{1} << 30) - 1;
        const u32 excluded_originals = label_mask({120, 126});
        const std::array<AnchorRow, 3> anchors = {{
            {"A40", label_mask({40, 143, 252})},
            {"A80", label_mask({80, 143, 252})},
            {"A240", label_mask({143, 240, 252})},
        }};
        require(3 * binomial(25, 7) == 1442100,
                "anchor-presentation formula changed");
        require(binomial(26, 8) - binomial(23, 8) == 1071961,
                "unique-body inclusion-exclusion formula changed");
        require(binomial(25, 7) + binomial(24, 7) + binomial(23, 7) ==
                    1071961,
                "priority partition formula changed");

        Ledger ledger;
        u64 total_presentations = 0;
        u64 total_below = 0;
        u64 total_equalities = 0;
        int global_q = 0;
        const char* global_anchor = nullptr;
        i64 global_num = 0;
        i64 global_den = 1;
        i128 global_delta = 0;
        u32 global_k = 0;

        std::cout << "AUDIT all_newcomer_one_original_pool_wall_failure_subset_zeta\n";
        std::cout << "OUTSIDE_Q7_CONTRACT A40_Q50_DEPTH 6 A40_Q50_EDGES 44562"
                  << " A40_TAU 8 A80_Q50_DEPTH 7 A80_Q50_EDGES 298279"
                  << " A80_TAU 8 A240_Q50_DEPTH 7 A240_Q50_EDGES 286291"
                  << " A240_TAU 8\n";
        std::cout << "CHOICE_RULE 40_FIRST_THEN_80_THEN_240"
                  << " PRIORITY_COUNTS 480700,346104,245157"
                  << " UNIQUE_BODIES_PER_Q 1071961"
                  << " UNIQUE_FORMULA C26_8-C23_8"
                  << " ANCHOR_PRESENTATIONS_PER_Q 1442100"
                  << " PRESENTATION_FORMULA 3*C25_7\n";
        for (int q : RESONANCES) {
            for (const AnchorRow& row : anchors) {
                const u32 allowed = full_pool & ~(row.anchor | excluded_originals);
                const BodyAudit audit = audit_bodies(cells, q, row.anchor, allowed);
                const i64 denominator = static_cast<i64>(
                    static_cast<i128>(14) * q * COMMON);
                total_presentations += audit.presentations;
                total_below += audit.below;
                total_equalities += audit.equalities;
                if (global_anchor == nullptr ||
                    static_cast<i128>(audit.minimum_numerator) * global_den <
                        static_cast<i128>(global_num) * denominator) {
                    global_q = q;
                    global_anchor = row.name;
                    global_num = audit.minimum_numerator;
                    global_den = denominator;
                    global_delta = audit.minimum_delta;
                    global_k = audit.minimum_k;
                }
                ledger.add(static_cast<u64>(q));
                ledger.add(static_cast<u64>(row.anchor));
                ledger.add(static_cast<u64>(audit.presentations));
                ledger.add(static_cast<u64>(audit.below));
                ledger.add(static_cast<u64>(audit.equalities));
                ledger.add(static_cast<u64>(audit.minimum_numerator));
                ledger.add(static_cast<u64>(
                    static_cast<i64>(audit.minimum_delta)));
                ledger.add(static_cast<u64>(audit.minimum_k));
                std::cout << "BODY_ROW q " << q
                          << " " << row.name
                          << " ANCHOR " << labels(row.anchor)
                          << " PRESENTATIONS " << audit.presentations
                          << " BELOW " << audit.below
                          << " EQUALITIES " << audit.equalities
                          << " MIN_MASS_NUM " << audit.minimum_numerator
                          << " MASS_DEN " << denominator
                          << " MIN_DELTA " << decimal(audit.minimum_delta)
                          << " MIN_K " << labels(audit.minimum_k)
                          << " MIN_BODY " << labels(row.anchor | audit.minimum_k)
                          << '\n';
            }
        }
        require(total_presentations == 33168300,
                "23x3 presentation total changed");
        require(total_below == 0 && total_equalities == 0,
                "unsafe or equality body appeared");
        require(global_q == 140 && std::string(global_anchor) == "A80" &&
                    global_num == INT64_C(5459802690013840) &&
                    global_den == INT64_C(35752672456300800) &&
                    global_delta == static_cast<i128>(INT64_C(28708125663666960)) &&
                    global_k == label_mask({10, 42, 60, 85, 95, 145, 264}),
                "global resonance minimum changed");
        require(ledger.value() == EXPECTED_LEDGER,
                "semantic ledger changed");
        std::cout << "SUMMARY ROWS 69 PRESENTATIONS " << total_presentations
                  << " UNIQUE_BODIES_PER_Q 1071961 UNIQUE_BODY_Q_PAIRS 24655103"
                  << " BELOW " << total_below
                  << " EQUALITIES " << total_equalities << '\n';
        std::cout << "GLOBAL_MIN q " << global_q
                  << " ANCHOR " << global_anchor
                  << " NUM " << global_num
                  << " DEN " << global_den
                  << " DELTA " << decimal(global_delta)
                  << " K " << labels(global_k)
                  << " BODY " << labels(label_mask({80, 143, 252}) | global_k)
                  << '\n';
        std::cout << "SEMANTIC_LEDGER_FNV1A64_LE " << hex64(ledger.value()) << '\n';
        std::cout << "VERDICT PASS\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
