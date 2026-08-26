// Primary exact audit for THM-4228.
//
// This path inherits the canonical THM-4188 fixed-pool cell geometry, groups
// failure atoms, performs an ordinary-colex superset zeta transform, and uses
// Gosper enumeration for both repair and body universes.  The independent
// audit uses a separately implemented wall geometry, reverse-colex storage,
// recursive subset enumeration, and ungrouped incidence contributions.

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main thm4188_primary_original_main
#include "lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp"
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace {

constexpr u64 EXPECTED_REPAIRS = UINT64_C(5852925);
constexpr u64 EXPECTED_BODIES = UINT64_C(14307150);
constexpr u64 EXPECTED_STRICT_EDGES = UINT64_C(5052990);
constexpr u64 EXPECTED_EDGES_3466 = UINT64_C(821120);
constexpr u64 EXPECTED_EDGES_3467 = UINT64_C(821923);
constexpr u64 EXPECTED_BODY_CHECKS = UINT64_C(44319858349);
constexpr u64 EXPECTED_MAX_CHECKS = UINT64_C(821129);
constexpr u64 EXPECTED_EDGE_LEDGER = UINT64_C(0xb270307777887d42);

std::array<std::array<u64, 9>, 31> choose8{};

void initialize_choose8() {
    for (int n = 0; n <= 30; ++n) {
        choose8[n][0] = 1;
        for (int k = 1; k <= 8; ++k) {
            choose8[n][k] = n == 0 ? 0 :
                choose8[n - 1][k] + choose8[n - 1][k - 1];
        }
    }
    require(choose8[30][8] == EXPECTED_REPAIRS,
            "rank-eight universe changed");
}

u64 colex_rank8(u32 mask) {
    u64 rank = 0;
    int ordinal = 1;
    for (int bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        rank += choose8[bit][ordinal];
        ++ordinal;
    }
    require(ordinal == 9 && rank < choose8[30][8],
            "invalid ordinary-colex rank");
    return rank;
}

u32 mask_for_labels(std::initializer_list<int> values) {
    u32 result = 0;
    for (int value : values) {
        const auto found = std::find(POOL.begin(), POOL.end(), value);
        require(found != POOL.end(), "control label absent from pool");
        result |= u32{1} << static_cast<int>(found - POOL.begin());
    }
    return result;
}

i128 gcd128(i128 left, i128 right) {
    if (left < 0) left = -left;
    if (right < 0) right = -right;
    while (right != 0) {
        const i128 remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

struct Fraction128 {
    i128 numerator;
    i128 denominator;
};

Fraction128 reduce_fraction(i128 numerator, i128 denominator) {
    require(denominator > 0, "nonpositive fraction denominator");
    const i128 divisor = gcd128(numerator, denominator);
    return {numerator / divisor, denominator / divisor};
}

// For coprime 1<=a<b, the exact danger-overlap mass is
// overlap_numerator(a,b)/(14ab).  The theorem proves the formula for all
// pairs; this routine supplies hostile finite and equality controls.
i128 overlap_numerator(int a, int b) {
    const int divisor = std::gcd(a, b);
    a /= divisor;
    b /= divisor;
    if (a > b) std::swap(a, b);
    require(a > 0 && a < b, "overlap control requires distinct positives");
    const int sum = a + b;
    i128 result = 2 * a;
    for (int k = 1; 14 * k < sum; ++k) {
        result += 2 * std::min(2 * a, sum - 14 * k);
    }
    return result;
}

struct AnalyticControls {
    u64 coprime_pairs = 0;
    u64 equality_pairs = 0;
};

AnalyticControls audit_overlap_formula() {
    AnalyticControls result;
    require(overlap_numerator(1, 13) == 2,
            "sharp 1:13 overlap numerator changed");
    require(overlap_numerator(1, 2) == 2,
            "1:2 overlap control changed");
    require(overlap_numerator(50, 51) == 730,
            "50:51 overlap control changed");
    require(static_cast<i128>(1650) * 819 * 91 ==
                static_cast<i128>(7425) * 8281 * 2,
            "activation-constant reduction changed");

    for (int b = 2; b <= 512; ++b) {
        for (int a = 1; a < b; ++a) {
            if (std::gcd(a, b) != 1) continue;
            ++result.coprime_pairs;
            const i128 numerator = overlap_numerator(a, b);
            const i128 denominator = static_cast<i128>(14) * a * b;
            require(91 * numerator >= denominator,
                    "finite pair-overlap hostile control failed");
            if (91 * numerator == denominator) {
                ++result.equality_pairs;
                require(a == 1 && b == 13,
                        "unexpected finite overlap equality");
            }
        }
    }
    require(result.equality_pairs == 1,
            "sharp overlap equality count changed");
    return result;
}

struct AtomStats {
    i64 length = 0;
    std::uint16_t cells = 0;
    std::uint16_t adjacencies = 0;
};

void add_supersets(u32 atom, int need, int start, u32 extra,
                   const AtomStats& value,
                   std::vector<i64>& lengths,
                   std::vector<std::uint16_t>& cells,
                   std::vector<std::uint16_t>& adjacencies,
                   u64& operations) {
    if (need == 0) {
        const u64 rank = colex_rank8(atom | extra);
        lengths[rank] += value.length;
        const unsigned new_cells =
            static_cast<unsigned>(cells[rank]) + value.cells;
        const unsigned new_adjacencies =
            static_cast<unsigned>(adjacencies[rank]) + value.adjacencies;
        require(new_cells <= std::numeric_limits<std::uint16_t>::max(),
                "safe-cell count overflow");
        require(new_adjacencies <=
                    std::numeric_limits<std::uint16_t>::max(),
                "safe-adjacency count overflow");
        cells[rank] = static_cast<std::uint16_t>(new_cells);
        adjacencies[rank] =
            static_cast<std::uint16_t>(new_adjacencies);
        ++operations;
        return;
    }
    for (int bit = start; bit <= 30 - need; ++bit) {
        const u32 flag = u32{1} << bit;
        if ((atom & flag) != 0) continue;
        add_supersets(atom, need - 1, bit + 1, extra | flag, value,
                      lengths, cells, adjacencies, operations);
    }
}

struct RayEdge {
    u64 cutoff = 0;
    u32 mask = 0;
    i64 mass = 0;
    std::uint16_t components = 0;
    i128 surplus = 0;
};

struct DirectStats {
    i64 mass = 0;
    i64 components = 0;
};

DirectStats direct_stats(const std::vector<Cell>& cells, u32 repair) {
    DirectStats result;
    std::vector<unsigned char> safe(cells.size(), 0);
    for (std::size_t index = 0; index < cells.size(); ++index) {
        safe[index] = static_cast<unsigned char>(
            (cells[index].failed_pool & ~repair) == 0);
        if (safe[index]) result.mass += cells[index].right - cells[index].left;
    }
    for (std::size_t index = 0; index < safe.size(); ++index) {
        const std::size_t previous =
            (index + safe.size() - 1) % safe.size();
        if (safe[index] && !safe[previous]) ++result.components;
    }
    require(result.mass > 0 && result.components > 0,
            "invalid direct witness geometry");
    return result;
}

u64 edge_ledger(const std::vector<RayEdge>& edges) {
    Fnv1a64 ledger;
    for (const RayEdge& edge : edges) {
        require(edge.surplus <= std::numeric_limits<u64>::max(),
                "semantic-ledger surplus overflow");
        ledger.add_u64_le(edge.cutoff);
        ledger.add_u64_le(edge.mask);
        ledger.add_u64_le(static_cast<u64>(edge.mass));
        ledger.add_u64_le(edge.components);
        ledger.add_u64_le(static_cast<u64>(edge.surplus));
    }
    return ledger.value();
}

}  // namespace

int main() {
    try {
        const AnalyticControls analytic = audit_overlap_formula();
        initialize_choose8();
        const std::vector<Cell> pool_cells = build_pool_cells();
        require(pool_cells.size() == 7133, "fixed-pool cell count changed");

        std::unordered_map<u32, AtomStats> atoms;
        atoms.reserve(8192);
        for (const Cell& cell : pool_cells) {
            if (std::popcount(cell.failed_pool) > 8) continue;
            AtomStats& value = atoms[cell.failed_pool];
            value.length += cell.right - cell.left;
            ++value.cells;
        }
        for (std::size_t index = 0; index < pool_cells.size(); ++index) {
            const Cell& previous = pool_cells[
                (index + pool_cells.size() - 1) % pool_cells.size()];
            const Cell& current = pool_cells[index];
            const u32 joined = previous.failed_pool | current.failed_pool;
            if (std::popcount(joined) <= 8) ++atoms[joined].adjacencies;
        }

        const u64 universe = choose8[30][8];
        std::vector<i64> lengths(universe, 0);
        std::vector<std::uint16_t> cells(universe, 0);
        std::vector<std::uint16_t> adjacencies(universe, 0);
        u64 zeta_operations = 0;
        for (const auto& [atom, value] : atoms) {
            add_supersets(atom, 8 - std::popcount(atom), 0, 0, value,
                          lengths, cells, adjacencies, zeta_operations);
        }

        std::vector<RayEdge> edges;
        edges.reserve(universe);
        u64 equalities = 0;
        u32 deletion = (u32{1} << 8) - 1;
        const u32 limit = u32{1} << 30;
        u64 rank = 0;
        while (deletion != 0 && deletion < limit) {
            const i64 mass = lengths[rank];
            const i64 components = static_cast<i64>(cells[rank]) -
                                   static_cast<i64>(adjacencies[rank]);
            require(mass > 0 && components > 0 && components <= 65535,
                    "invalid base safe-set geometry");
            const i128 surplus = static_cast<i128>(297) * mass -
                                 static_cast<i128>(26) * COMMON;
            if (surplus == 0) ++equalities;
            if (surplus > 0) {
                const i128 numerator = static_cast<i128>(7425) *
                                       components * COMMON;
                const i128 denominator = static_cast<i128>(91) * surplus;
                const i128 cutoff =
                    (numerator + denominator - 1) / denominator;
                require(cutoff > 0 &&
                            cutoff <= std::numeric_limits<u64>::max(),
                        "activation cutoff overflow");
                edges.push_back({static_cast<u64>(cutoff), deletion, mass,
                                 static_cast<std::uint16_t>(components),
                                 surplus});
            }
            ++rank;
            const u32 next = next_combination(deletion);
            if (next <= deletion) break;
            deletion = next;
        }
        require(rank == universe, "rank-eight repair universe incomplete");
        require(edges.size() == EXPECTED_STRICT_EDGES && equalities == 0,
                "strict ray deck changed");
        std::sort(edges.begin(), edges.end(),
                  [](const RayEdge& left, const RayEdge& right) {
                      if (left.cutoff != right.cutoff) {
                          return left.cutoff < right.cutoff;
                      }
                      return left.mask < right.mask;
                  });
        const u64 semantic_ledger = edge_ledger(edges);
        require(semantic_ledger == EXPECTED_EDGE_LEDGER,
                "primary semantic edge ledger changed");

        u64 activation = 0;
        u32 hostile_body = 0;
        RayEdge witness{};
        u64 bodies = 0;
        u64 checks = 0;
        u64 maximum_checks = 0;
        u32 body = (u32{1} << 9) - 1;
        while (body < limit) {
            u64 local = 0;
            const RayEdge* first = nullptr;
            for (const RayEdge& edge : edges) {
                ++local;
                if ((body & edge.mask) == 0) {
                    first = &edge;
                    break;
                }
            }
            require(first != nullptr,
                    "nine-body covers every strict ray edge");
            checks += local;
            maximum_checks = std::max(maximum_checks, local);
            ++bodies;
            if (first->cutoff > activation ||
                (first->cutoff == activation && body < hostile_body)) {
                activation = first->cutoff;
                hostile_body = body;
                witness = *first;
            }
            const u32 next = next_combination(body);
            if (next <= body) break;
            body = next;
        }
        require(bodies == EXPECTED_BODIES && checks == EXPECTED_BODY_CHECKS &&
                    maximum_checks == EXPECTED_MAX_CHECKS,
                "nine-body scan changed");
        require(activation == 3467,
                "sharp sufficient-certificate activation changed");

        const auto active_count = [&](u64 g) {
            return static_cast<u64>(std::upper_bound(
                edges.begin(), edges.end(), g,
                [](u64 value, const RayEdge& edge) {
                    return value < edge.cutoff;
                }) - edges.begin());
        };
        const u64 edges_3466 = active_count(3466);
        const u64 edges_3467 = active_count(3467);
        require(edges_3466 == EXPECTED_EDGES_3466 &&
                    edges_3467 == EXPECTED_EDGES_3467,
                "activation filtration counts changed");
        for (u64 index = 0; index < edges_3466; ++index) {
            require((hostile_body & edges[index].mask) != 0,
                    "predecessor body is not an E_3466 cover");
        }

        const u32 expected_body = mask_for_labels(
            {88,95,170,193,240,252,264,286,290});
        const u32 expected_repair = mask_for_labels(
            {8,16,85,132,143,145,168,176});
        require(hostile_body == expected_body && witness.mask == expected_repair &&
                    witness.cutoff == 3467 &&
                    (hostile_body & witness.mask) == 0,
                "sharp predecessor/activating witness changed");
        require(witness.mass == INT64_C(1920678988176) &&
                    witness.components == 224 &&
                    witness.surplus == static_cast<i128>(96171514659792LL),
                "activating witness constants changed");
        const DirectStats direct = direct_stats(pool_cells, witness.mask);
        require(direct.mass == witness.mass &&
                    direct.components == witness.components,
                "direct witness geometry disagrees with zeta geometry");

        const Fraction128 ratio = reduce_fraction(
            static_cast<i128>(7425) * witness.components * COMMON,
            static_cast<i128>(91) * witness.surplus);
        require(ratio.numerator == static_cast<i128>(210474916344000LL) &&
                    ratio.denominator == static_cast<i128>(60714340063LL),
                "activating witness ratio changed");
        const auto residual = [&](i64 g) {
            const i128 numerator =
                static_cast<i128>(2) * witness.surplus * 8281 * g -
                static_cast<i128>(1650) * witness.components * 819 * COMMON;
            const i128 denominator =
                static_cast<i128>(819) * COMMON * 8281 * g;
            return reduce_fraction(numerator, denominator);
        };
        const Fraction128 residual_3466 = residual(3466);
        const Fraction128 residual_3467 = residual(3467);
        require(residual_3466.numerator ==
                    -static_cast<i128>(19506842821LL) &&
                    residual_3466.denominator ==
                    static_cast<i128>(8172402168912345LL) &&
                    residual_3467.numerator ==
                    static_cast<i128>(21700654421LL) &&
                    residual_3467.denominator ==
                    static_cast<i128>(16349520092105655LL),
                "boundary residuals changed");

        std::cout << "THM-4228 primary pool-cell activation audit\n";
        std::cout << "pair_overlap_controls=" << analytic.coprime_pairs
                  << " equalities=" << analytic.equality_pairs
                  << " beta_floor=66/91 osc_ceiling=1650/8281"
                  << " equality_ratio=1:13\n";
        std::cout << "fixed_cells=" << pool_cells.size()
                  << " atom_support=" << atoms.size()
                  << " zeta_operations=" << zeta_operations << '\n';
        std::cout << "repairs=" << universe
                  << " strict_edges=" << edges.size()
                  << " equalities=" << equalities
                  << " semantic_ledger=" << hex64(semantic_ledger) << '\n';
        std::cout << "bodies=" << bodies
                  << " checks=" << checks
                  << " max_checks=" << maximum_checks
                  << " E3466=" << edges_3466
                  << " E3467=" << edges_3467
                  << " sharp_certificate_g=" << activation << '\n';
        std::cout << "predecessor_cover={" << labels(hostile_body) << "}"
                  << " activating_repair={" << labels(witness.mask) << "}\n";
        std::cout << "repair_mass=" << witness.mass << '/' << COMMON
                  << " repair_components=" << witness.components
                  << " surplus297=" << decimal(witness.surplus) << '\n';
        std::cout << "activation_ratio=" << decimal(ratio.numerator) << '/'
                  << decimal(ratio.denominator)
                  << " residual3466=" << decimal(residual_3466.numerator)
                  << '/' << decimal(residual_3466.denominator)
                  << " residual3467=" << decimal(residual_3467.numerator)
                  << '/' << decimal(residual_3467.denominator) << '\n';
        std::cout << "checks=PASS analytic_constants,overlap_hostiles,"
                     "pool_cells,all_repairs,activation_filtration,all_bodies,"
                     "predecessor_cover,direct_witness\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
