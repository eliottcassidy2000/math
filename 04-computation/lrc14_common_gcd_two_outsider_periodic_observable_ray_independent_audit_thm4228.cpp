// Independent exact audit for THM-4228.
//
// Unlike the primary grouped-atom/colex path, this program inherits the
// independently implemented THM-4188 fixed-wall geometry, recursively
// enumerates subsets, stores repairs in reverse colex, and scatters every
// fixed cell and adjacency separately.  It also constructs literal joint
// walls for a hostile bank of primitive pairs rather than using the closed
// overlap formula implemented by the primary audit.

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main thm4188_independent_original_main
#include "lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.cpp"
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
#include <unordered_set>
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

std::array<std::array<u64, 9>, 31> binomial{};

void initialize_binomial() {
    for (int n = 0; n <= 30; ++n) {
        binomial[n][0] = 1;
        for (int k = 1; k <= 8; ++k) {
            binomial[n][k] = n == 0 ? 0 :
                binomial[n - 1][k] + binomial[n - 1][k - 1];
        }
    }
    require(binomial[30][8] == EXPECTED_REPAIRS,
            "rank-eight universe changed");
}

// This reverse-colex bijection differs from the primary ordinary-colex
// indexing even though both arrays have the same cardinality.
u64 reverse_colex_rank8(u32 mask) {
    u64 rank = 0;
    int ordinal = 1;
    for (int reflected = 0; reflected < 30; ++reflected) {
        const int original = 29 - reflected;
        if ((mask & (u32{1} << original)) == 0) continue;
        rank += binomial[reflected][ordinal];
        ++ordinal;
    }
    require(ordinal == 9 && rank < binomial[30][8],
            "invalid reverse-colex rank");
    return rank;
}

u32 mask_for_labels_independent(std::initializer_list<int> values) {
    u32 result = 0;
    for (int value : values) {
        const auto found = std::find(POOL.begin(), POOL.end(), value);
        require(found != POOL.end(), "control label absent from pool");
        result |= u32{1} << static_cast<int>(found - POOL.begin());
    }
    return result;
}

i128 gcd128_independent(i128 left, i128 right) {
    if (left < 0) left = -left;
    if (right < 0) right = -right;
    while (right != 0) {
        const i128 remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

struct Fraction128Independent {
    i128 numerator;
    i128 denominator;
};

Fraction128Independent reduce_independent(i128 numerator,
                                          i128 denominator) {
    require(denominator > 0, "nonpositive fraction denominator");
    const i128 divisor = gcd128_independent(numerator, denominator);
    return {numerator / divisor, denominator / divisor};
}

struct PairWallStats {
    i64 denominator = 0;
    i64 safe_length = 0;
    i128 primitive_oscillation_numerator = 0;
    u64 cells = 0;
};

PairWallStats literal_pair_stats(int first, int second) {
    const int divisor = std::gcd(first, second);
    first /= divisor;
    second /= divisor;
    require(first > 0 && second > 0 && first != second,
            "literal pair control requires distinct positives");
    const i64 denominator = checked_lcm(14LL * first, 14LL * second);
    std::vector<i64> walls = {0, denominator};
    for (int speed : {first, second}) {
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    std::vector<unsigned char> safe(walls.size() - 1, 0);
    i64 safe_length = 0;
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        safe[index] = static_cast<unsigned char>(
            safe_at_midpoint(first, walls[index], walls[index + 1], denominator) &&
            safe_at_midpoint(second, walls[index], walls[index + 1], denominator));
        if (safe[index]) safe_length += walls[index + 1] - walls[index];
    }

    i128 primitive = 0;
    i128 minimum = 0;
    i128 maximum = 0;
    for (std::size_t index = 0; index < safe.size(); ++index) {
        const i64 length = walls[index + 1] - walls[index];
        primitive += (static_cast<i128>(safe[index] ? denominator : 0) -
                      safe_length) * length;
        minimum = std::min(minimum, primitive);
        maximum = std::max(maximum, primitive);
    }
    require(primitive == 0, "centered pair primitive did not close");
    return {denominator, safe_length, maximum - minimum,
            static_cast<u64>(safe.size())};
}

struct PairControls {
    u64 coprime_pairs = 0;
    u64 beta_equalities = 0;
    i128 maximum_oscillation_numerator = 0;
    i128 maximum_oscillation_denominator = 1;
    int maximum_first = 0;
    int maximum_second = 0;
};

PairControls audit_literal_pairs() {
    PairControls result;
    for (int second = 2; second <= 128; ++second) {
        for (int first = 1; first < second; ++first) {
            if (std::gcd(first, second) != 1) continue;
            ++result.coprime_pairs;
            const PairWallStats stats = literal_pair_stats(first, second);
            require(static_cast<i128>(91) * stats.safe_length >=
                        static_cast<i128>(66) * stats.denominator,
                    "literal pair violates beta>=66/91");
            if (static_cast<i128>(91) * stats.safe_length ==
                static_cast<i128>(66) * stats.denominator) {
                ++result.beta_equalities;
                require(first == 1 && second == 13,
                        "unexpected literal beta equality");
            }
            const i128 square =
                static_cast<i128>(stats.denominator) * stats.denominator;
            require(static_cast<i128>(8281) *
                        stats.primitive_oscillation_numerator <=
                    static_cast<i128>(1650) * square,
                    "literal pair violates primitive oscillation ceiling");
            if (result.maximum_first == 0 ||
                stats.primitive_oscillation_numerator *
                        result.maximum_oscillation_denominator >
                    result.maximum_oscillation_numerator * square) {
                result.maximum_oscillation_numerator =
                    stats.primitive_oscillation_numerator;
                result.maximum_oscillation_denominator = square;
                result.maximum_first = first;
                result.maximum_second = second;
            }
        }
    }
    require(result.beta_equalities == 1,
            "literal beta equality count changed");
    const Fraction128Independent maximum = reduce_independent(
        result.maximum_oscillation_numerator,
        result.maximum_oscillation_denominator);
    require(result.maximum_first == 1 && result.maximum_second == 13 &&
                maximum.numerator == 990 && maximum.denominator == 8281,
            "finite literal primitive maximum changed");
    return result;
}

struct DenseStats {
    std::vector<i64> length;
    std::vector<std::uint16_t> safe_cells;
    std::vector<std::uint16_t> safe_adjacencies;
    u64 contribution_count = 0;
    u64 support_count = 0;
};

void contribute_supersets(const std::vector<int>& available, int need,
                          int next, u32 mask, i64 length,
                          unsigned cells, unsigned adjacencies,
                          DenseStats& result) {
    if (need == 0) {
        const u64 rank = reverse_colex_rank8(mask);
        result.length[rank] += length;
        const unsigned new_cells = result.safe_cells[rank] + cells;
        const unsigned new_adjacencies =
            result.safe_adjacencies[rank] + adjacencies;
        require(new_cells <= std::numeric_limits<std::uint16_t>::max(),
                "safe-cell count overflow");
        require(new_adjacencies <=
                    std::numeric_limits<std::uint16_t>::max(),
                "safe-adjacency count overflow");
        result.safe_cells[rank] =
            static_cast<std::uint16_t>(new_cells);
        result.safe_adjacencies[rank] =
            static_cast<std::uint16_t>(new_adjacencies);
        ++result.contribution_count;
        return;
    }
    for (int index = next;
         index <= static_cast<int>(available.size()) - need; ++index) {
        contribute_supersets(available, need - 1, index + 1,
                             mask | (u32{1} << available[index]),
                             length, cells, adjacencies, result);
    }
}

void contribute(u32 forced, i64 length, unsigned cells,
                unsigned adjacencies, DenseStats& result) {
    const int fixed = std::popcount(forced);
    if (fixed > 8) return;
    std::vector<int> available;
    available.reserve(30 - fixed);
    for (int bit = 0; bit < 30; ++bit) {
        if ((forced & (u32{1} << bit)) == 0) available.push_back(bit);
    }
    contribute_supersets(available, 8 - fixed, 0, forced, length,
                         cells, adjacencies, result);
}

DenseStats build_dense_stats_independently(const FixedGeometry& fixed) {
    DenseStats result;
    const u64 universe = binomial[30][8];
    result.length.assign(universe, 0);
    result.safe_cells.assign(universe, 0);
    result.safe_adjacencies.assign(universe, 0);

    std::unordered_set<u32> support;
    support.reserve(8192);
    for (const FixedCell& cell : fixed.cells) {
        if (std::popcount(cell.failed) <= 8) support.insert(cell.failed);
        contribute(cell.failed, cell.right - cell.left, 1, 0, result);
    }
    for (std::size_t index = 0; index < fixed.cells.size(); ++index) {
        const u32 joined = fixed.cells[index].failed |
            fixed.cells[(index + 1) % fixed.cells.size()].failed;
        if (std::popcount(joined) <= 8) support.insert(joined);
        contribute(joined, 0, 0, 1, result);
    }
    result.support_count = support.size();
    return result;
}

struct ActivatedEdge {
    u64 cutoff;
    u32 repair;
    i64 mass;
    std::uint16_t components;
    i128 surplus;
};

struct DirectStats {
    i64 mass = 0;
    i64 components = 0;
};

DirectStats direct_stats(const FixedGeometry& fixed, u32 repair) {
    DirectStats result;
    std::vector<unsigned char> safe(fixed.cells.size(), 0);
    for (std::size_t index = 0; index < fixed.cells.size(); ++index) {
        safe[index] = static_cast<unsigned char>(
            (fixed.cells[index].failed & ~repair) == 0);
        if (safe[index]) {
            result.mass += fixed.cells[index].right - fixed.cells[index].left;
        }
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

u64 edge_ledger_independent(const std::vector<ActivatedEdge>& edges) {
    Ledger ledger;
    for (const ActivatedEdge& edge : edges) {
        require(edge.surplus <= std::numeric_limits<u64>::max(),
                "semantic-ledger surplus overflow");
        ledger.add(edge.cutoff);
        ledger.add(edge.repair);
        ledger.add(static_cast<u64>(edge.mass));
        ledger.add(edge.components);
        ledger.add(static_cast<u64>(edge.surplus));
    }
    return ledger.value();
}

}  // namespace

int main() {
    try {
        const PairControls pair_controls = audit_literal_pairs();
        initialize_binomial();
        const FixedGeometry fixed = build_fixed_geometry();
        require(fixed.walls.size() == 7134 && fixed.cells.size() == 7133,
                "independent fixed geometry changed");
        const DenseStats dense = build_dense_stats_independently(fixed);
        require(dense.support_count == 2721 &&
                    dense.contribution_count == UINT64_C(2493648782),
                "independent incidence scatter changed");

        std::vector<ActivatedEdge> edges;
        edges.reserve(binomial[30][8]);
        u64 repairs = 0;
        u64 equalities = 0;
        for_each_k_subset(30, 8, [&](u32 repair) {
            const u64 rank = reverse_colex_rank8(repair);
            const i64 mass = dense.length[rank];
            const i64 components =
                static_cast<i64>(dense.safe_cells[rank]) -
                static_cast<i64>(dense.safe_adjacencies[rank]);
            require(mass > 0 && components > 0 && components <= 65535,
                    "invalid independently aggregated geometry");
            const i128 surplus = static_cast<i128>(297) * mass -
                                 static_cast<i128>(26) * POOL_DENOMINATOR;
            if (surplus == 0) ++equalities;
            if (surplus > 0) {
                const i128 numerator = static_cast<i128>(7425) *
                                       components * POOL_DENOMINATOR;
                const i128 denominator = static_cast<i128>(91) * surplus;
                const i128 cutoff =
                    (numerator + denominator - 1) / denominator;
                require(cutoff > 0 &&
                            cutoff <= std::numeric_limits<u64>::max(),
                        "independent activation cutoff overflow");
                edges.push_back({static_cast<u64>(cutoff), repair, mass,
                                 static_cast<std::uint16_t>(components),
                                 surplus});
            }
            ++repairs;
        });
        require(repairs == EXPECTED_REPAIRS &&
                    edges.size() == EXPECTED_STRICT_EDGES && equalities == 0,
                "independent strict ray deck changed");
        std::sort(edges.begin(), edges.end(),
                  [](const ActivatedEdge& left, const ActivatedEdge& right) {
                      if (left.cutoff != right.cutoff) {
                          return left.cutoff < right.cutoff;
                      }
                      return left.repair < right.repair;
                  });
        const u64 semantic_ledger = edge_ledger_independent(edges);
        require(semantic_ledger == EXPECTED_EDGE_LEDGER,
                "independent semantic edge ledger changed");

        u64 activation = 0;
        u32 hostile_body = 0;
        ActivatedEdge witness{};
        u64 bodies = 0;
        u64 checks = 0;
        u64 maximum_checks = 0;
        for_each_k_subset(30, 9, [&](u32 body) {
            u64 local = 0;
            const ActivatedEdge* first = nullptr;
            for (const ActivatedEdge& edge : edges) {
                ++local;
                if ((body & edge.repair) == 0) {
                    first = &edge;
                    break;
                }
            }
            require(first != nullptr,
                    "nine-body covers every independent strict ray edge");
            checks += local;
            maximum_checks = std::max(maximum_checks, local);
            ++bodies;
            if (first->cutoff > activation ||
                (first->cutoff == activation && body < hostile_body)) {
                activation = first->cutoff;
                hostile_body = body;
                witness = *first;
            }
        });
        require(bodies == EXPECTED_BODIES && checks == EXPECTED_BODY_CHECKS &&
                    maximum_checks == EXPECTED_MAX_CHECKS &&
                    activation == 3467,
                "independent nine-body activation scan changed");

        const auto active_count = [&](u64 g) {
            return static_cast<u64>(std::upper_bound(
                edges.begin(), edges.end(), g,
                [](u64 value, const ActivatedEdge& edge) {
                    return value < edge.cutoff;
                }) - edges.begin());
        };
        const u64 edges_3466 = active_count(3466);
        const u64 edges_3467 = active_count(3467);
        require(edges_3466 == EXPECTED_EDGES_3466 &&
                    edges_3467 == EXPECTED_EDGES_3467,
                "independent activation filtration changed");
        for (u64 index = 0; index < edges_3466; ++index) {
            require((hostile_body & edges[index].repair) != 0,
                    "independent predecessor is not an E_3466 cover");
        }

        const u32 expected_body = mask_for_labels_independent(
            {88,95,170,193,240,252,264,286,290});
        const u32 expected_repair = mask_for_labels_independent(
            {8,16,85,132,143,145,168,176});
        require(hostile_body == expected_body &&
                    witness.repair == expected_repair &&
                    witness.cutoff == 3467 &&
                    (hostile_body & witness.repair) == 0,
                "independent sharp witness changed");
        require(witness.mass == INT64_C(1920678988176) &&
                    witness.components == 224 &&
                    witness.surplus == static_cast<i128>(96171514659792LL),
                "independent witness constants changed");
        const DirectStats direct = direct_stats(fixed, witness.repair);
        require(direct.mass == witness.mass &&
                    direct.components == witness.components,
                "independent direct witness disagrees with scatter geometry");

        const Fraction128Independent ratio = reduce_independent(
            static_cast<i128>(7425) * witness.components * POOL_DENOMINATOR,
            static_cast<i128>(91) * witness.surplus);
        const auto residual = [&](i64 g) {
            const i128 numerator =
                static_cast<i128>(2) * witness.surplus * 8281 * g -
                static_cast<i128>(1650) * witness.components * 819 *
                    POOL_DENOMINATOR;
            const i128 denominator = static_cast<i128>(819) *
                POOL_DENOMINATOR * 8281 * g;
            return reduce_independent(numerator, denominator);
        };
        const Fraction128Independent residual_3466 = residual(3466);
        const Fraction128Independent residual_3467 = residual(3467);
        require(ratio.numerator == static_cast<i128>(210474916344000LL) &&
                    ratio.denominator == static_cast<i128>(60714340063LL) &&
                    residual_3466.numerator ==
                        -static_cast<i128>(19506842821LL) &&
                    residual_3466.denominator ==
                        static_cast<i128>(8172402168912345LL) &&
                    residual_3467.numerator ==
                        static_cast<i128>(21700654421LL) &&
                    residual_3467.denominator ==
                        static_cast<i128>(16349520092105655LL),
                "independent boundary arithmetic changed");

        const Fraction128Independent maximum_pair_oscillation =
            reduce_independent(pair_controls.maximum_oscillation_numerator,
                               pair_controls.maximum_oscillation_denominator);
        std::cout << "THM-4228 independent fixed-wall activation audit\n";
        std::cout << "literal_pair_controls=" << pair_controls.coprime_pairs
                  << " beta_equalities=" << pair_controls.beta_equalities
                  << " beta_floor=66/91 osc_ceiling=1650/8281\n";
        std::cout << "finite_max_pair_osc="
                  << decimal(maximum_pair_oscillation.numerator) << '/'
                  << decimal(maximum_pair_oscillation.denominator)
                  << " at=" << pair_controls.maximum_first << ':'
                  << pair_controls.maximum_second << '\n';
        std::cout << "fixed_walls=" << fixed.walls.size()
                  << " fixed_cells=" << fixed.cells.size()
                  << " support=" << dense.support_count
                  << " raw_contributions=" << dense.contribution_count << '\n';
        std::cout << "repairs=" << repairs
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
                  << " activating_repair={" << labels(witness.repair) << "}\n";
        std::cout << "repair_mass=" << witness.mass << '/'
                  << POOL_DENOMINATOR
                  << " repair_components=" << witness.components
                  << " surplus297=" << decimal(witness.surplus) << '\n';
        std::cout << "activation_ratio=" << decimal(ratio.numerator) << '/'
                  << decimal(ratio.denominator)
                  << " residual3466=" << decimal(residual_3466.numerator)
                  << '/' << decimal(residual_3466.denominator)
                  << " residual3467=" << decimal(residual_3467.numerator)
                  << '/' << decimal(residual_3467.denominator) << '\n';
        std::cout << "checks=PASS literal_pair_walls,independent_fixed_walls,"
                     "reverse_colex,raw_incidence_scatter,all_repairs,"
                     "semantic_ledger,activation_filtration,all_bodies,"
                     "predecessor_cover,direct_witness\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
