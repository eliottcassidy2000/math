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

// Independent exact audit for THM-4190's all-newcomer exactly-one-original slice.
//
// The primary path keeps the fixed P-wall cells, integrates the newcomer comb
// by a prefix primitive, projects FAILURE masks, and performs a subset zeta.
// This path instead constructs the complete joint P-union-{q} wall refinement,
// reclassifies every pool runner on every joint atom, projects SAFE masks, and
// performs the dual superset zeta.  Every reported minimizing body is then
// measured once more by a literal runner-by-runner joint-wall sweep.

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr i64 POOL_DENOMINATOR = INT64_C(18241159416480);
constexpr std::array<int, 23> RESONANCES = {
    6, 22, 24, 25, 48, 70, 72, 96, 100, 105, 110, 128,
    130, 140, 186, 192, 206, 210, 220, 256, 260, 294, 366};
constexpr u64 FNV_OFFSET = UINT64_C(0xcbf29ce484222325);
constexpr u64 FNV_PRIME = UINT64_C(0x100000001b3);
constexpr u64 EXPECTED_LEDGER = UINT64_C(0x72db7e3090d5626e);

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

bool safe_at_midpoint(int speed, i64 left, i64 right, i64 denominator) {
    const i128 raw = static_cast<i128>(speed) * (left + right);
    const i64 residue = static_cast<i64>(raw % (static_cast<i128>(2) * denominator));
    return static_cast<i128>(7) * residue >= denominator &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * denominator;
}

struct JointAtom {
    i64 length;
    u32 safe_pool;
};

struct JointGeometry {
    i64 denominator;
    std::vector<i64> walls;
    std::vector<JointAtom> q_safe_atoms;
    i64 q_safe_mass;
};

JointGeometry build_joint_geometry(int q) {
    JointGeometry geometry{1, {0}, {}, 0};
    for (int speed : POOL) {
        geometry.denominator = checked_lcm(geometry.denominator, 14LL * speed);
    }
    require(geometry.denominator == POOL_DENOMINATOR,
            "fixed pool denominator changed");
    geometry.denominator = checked_lcm(geometry.denominator, 14LL * q);
    geometry.walls.push_back(geometry.denominator);

    const auto add_walls = [&](int speed) {
        const i64 unit = geometry.denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            geometry.walls.push_back((14LL * tooth + 1) * unit);
            geometry.walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : POOL) add_walls(speed);
    add_walls(q);
    std::sort(geometry.walls.begin(), geometry.walls.end());
    geometry.walls.erase(
        std::unique(geometry.walls.begin(), geometry.walls.end()),
        geometry.walls.end());

    for (std::size_t index = 0; index + 1 < geometry.walls.size(); ++index) {
        const i64 left = geometry.walls[index];
        const i64 right = geometry.walls[index + 1];
        if (!safe_at_midpoint(q, left, right, geometry.denominator)) continue;
        u32 safe_pool = 0;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if (safe_at_midpoint(POOL[vertex], left, right,
                                 geometry.denominator)) {
                safe_pool |= u32{1} << vertex;
            }
        }
        geometry.q_safe_atoms.push_back({right - left, safe_pool});
        geometry.q_safe_mass += right - left;
    }
    require(static_cast<i128>(7) * geometry.q_safe_mass ==
                static_cast<i128>(6) * geometry.denominator,
            "newcomer safe mass is not 6/7");
    return geometry;
}

struct BodyAudit {
    u64 presentations = 0;
    u64 below = 0;
    u64 equalities = 0;
    i64 minimum_joint_mass = 0;
    i64 minimum_primary_numerator = 0;
    i128 minimum_primary_delta = 0;
    u32 minimum_k = 0;
};

BodyAudit audit_bodies(const JointGeometry& geometry,
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
    std::vector<i64> superset_mass(local_size, 0);
    for (const JointAtom& atom : geometry.q_safe_atoms) {
        if ((atom.safe_pool & anchor) != anchor) continue;
        u32 local_safe = 0;
        u32 safe = atom.safe_pool & allowed;
        while (safe != 0) {
            const u32 bit = safe & (~safe + 1u);
            local_safe |= u32{1} << local_index[std::countr_zero(bit)];
            safe ^= bit;
        }
        superset_mass[local_safe] += atom.length;
    }

    // Dual to the primary failure-mask subset zeta: after this transform,
    // superset_mass[K] is the weight of all joint atoms safe on every label K.
    for (int bit = 0; bit < 25; ++bit) {
        const u32 half = u32{1} << bit;
        const u32 block = half << 1;
        for (u32 base = 0; base < local_size; base += block) {
            for (u32 offset = 0; offset < half; ++offset) {
                superset_mass[base + offset] +=
                    superset_mass[base + half + offset];
            }
        }
    }

    const i64 primary_denominator = static_cast<i64>(
        static_cast<i128>(14) * q * POOL_DENOMINATOR);
    require(primary_denominator % geometry.denominator == 0,
            "joint denominator does not divide primary normalization");
    const i64 scale = primary_denominator / geometry.denominator;

    BodyAudit result;
    bool first = true;
    u32 local_k = (u32{1} << 7) - 1;
    while (local_k != 0 && local_k < local_size) {
        const i64 joint_mass = superset_mass[local_k];
        const i64 primary_numerator = static_cast<i64>(
            static_cast<i128>(joint_mass) * scale);
        const i128 primary_delta = static_cast<i128>(9) * primary_numerator -
                                   static_cast<i128>(8) * q * POOL_DENOMINATOR;
        u32 global_k = 0;
        u32 bits = local_k;
        while (bits != 0) {
            const u32 bit = bits & (~bits + 1u);
            global_k |= u32{1} << local_to_pool[std::countr_zero(bit)];
            bits ^= bit;
        }
        ++result.presentations;
        if (primary_delta < 0) ++result.below;
        if (primary_delta == 0) ++result.equalities;
        if (first || primary_delta < result.minimum_primary_delta) {
            result.minimum_joint_mass = joint_mass;
            result.minimum_primary_numerator = primary_numerator;
            result.minimum_primary_delta = primary_delta;
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

i64 direct_body_mass(const JointGeometry& geometry,
                     int q,
                     u32 body) {
    i64 result = 0;
    for (std::size_t index = 0; index + 1 < geometry.walls.size(); ++index) {
        const i64 left = geometry.walls[index];
        const i64 right = geometry.walls[index + 1];
        if (!safe_at_midpoint(q, left, right, geometry.denominator)) continue;
        bool safe = true;
        u32 runners = body;
        while (runners != 0) {
            const u32 bit = runners & (~runners + 1u);
            if (!safe_at_midpoint(POOL[std::countr_zero(bit)], left, right,
                                  geometry.denominator)) {
                safe = false;
                break;
            }
            runners ^= bit;
        }
        if (safe) result += right - left;
    }
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
        const u32 full_pool = (u32{1} << 30) - 1;
        const u32 excluded_originals = label_mask({120, 126});
        const std::array<AnchorRow, 3> anchors = {{
            {"A40", label_mask({40, 143, 252})},
            {"A80", label_mask({80, 143, 252})},
            {"A240", label_mask({143, 240, 252})},
        }};

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

        std::cout << "AUDIT all_newcomer_one_original_joint_wall_safe_superset_zeta\n";
        for (int q : RESONANCES) {
            const JointGeometry geometry = build_joint_geometry(q);
            std::cout << "JOINT q " << q
                      << " DENOMINATOR " << geometry.denominator
                      << " WALLS " << geometry.walls.size()
                      << " Q_SAFE_ATOMS " << geometry.q_safe_atoms.size()
                      << '\n';
            for (const AnchorRow& row : anchors) {
                const u32 allowed = full_pool & ~(row.anchor | excluded_originals);
                const BodyAudit audit = audit_bodies(geometry, q, row.anchor, allowed);
                const i64 primary_denominator = static_cast<i64>(
                    static_cast<i128>(14) * q * POOL_DENOMINATOR);
                const u32 body = row.anchor | audit.minimum_k;
                const i64 direct = direct_body_mass(geometry, q, body);
                require(direct == audit.minimum_joint_mass,
                        "direct minimizer sweep disagrees with superset zeta");
                total_presentations += audit.presentations;
                total_below += audit.below;
                total_equalities += audit.equalities;
                if (global_anchor == nullptr ||
                    static_cast<i128>(audit.minimum_primary_numerator) * global_den <
                        static_cast<i128>(global_num) * primary_denominator) {
                    global_q = q;
                    global_anchor = row.name;
                    global_num = audit.minimum_primary_numerator;
                    global_den = primary_denominator;
                    global_delta = audit.minimum_primary_delta;
                    global_k = audit.minimum_k;
                }
                ledger.add(static_cast<u64>(q));
                ledger.add(static_cast<u64>(row.anchor));
                ledger.add(static_cast<u64>(audit.presentations));
                ledger.add(static_cast<u64>(audit.below));
                ledger.add(static_cast<u64>(audit.equalities));
                ledger.add(static_cast<u64>(audit.minimum_primary_numerator));
                ledger.add(static_cast<u64>(
                    static_cast<i64>(audit.minimum_primary_delta)));
                ledger.add(static_cast<u64>(audit.minimum_k));
                std::cout << "BODY_ROW " << row.name
                          << " ANCHOR " << labels(row.anchor)
                          << " PRESENTATIONS " << audit.presentations
                          << " BELOW " << audit.below
                          << " EQUALITIES " << audit.equalities
                          << " MIN_MASS_NUM " << audit.minimum_primary_numerator
                          << " MASS_DEN " << primary_denominator
                          << " MIN_DELTA " << decimal(audit.minimum_primary_delta)
                          << " MIN_K " << labels(audit.minimum_k)
                          << " MIN_BODY " << labels(body)
                          << " DIRECT_JOINT_MASS " << direct
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
