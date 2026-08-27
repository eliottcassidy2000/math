#include "independent_common.hpp"

#include <unordered_map>

namespace {
using namespace ray_audit;

struct Atom {
    i128 width = 0;
    i128 cells = 0;
};

struct ChooseTable {
    std::array<std::array<u64, 9>, 31> c{};
    ChooseTable() {
        for (int n = 0; n <= 30; ++n) {
            c[n][0] = 1;
            for (int k = 1; k <= 8; ++k) {
                c[n][k] = n == 0 ? 0 : c[n - 1][k] + c[n - 1][k - 1];
            }
        }
        ensure(c[30][8] == REPAIR_COUNT, "rank-eight binomial mismatch");
    }
    u64 rank8(u32 mask) const {
        ensure(std::popcount(mask) == 8, "rank called outside arity eight");
        u64 result = 0;
        int ordinal = 1;
        for (int bit = 0; bit < 30; ++bit) {
            if ((mask & (u32{1} << bit)) == 0) continue;
            result += c[bit][ordinal];
            ++ordinal;
        }
        ensure(result < REPAIR_COUNT, "repair rank overflow");
        return result;
    }
};

template <class Callback>
void enumerate_completions(u32 fixed, int need, int start, u32 additions,
                           const Callback& callback) {
    if (need == 0) {
        callback(fixed | additions);
        return;
    }
    for (int bit = start; bit < 30; ++bit) {
        const u32 flag = u32{1} << bit;
        if ((fixed & flag) != 0) continue;
        int available = 0;
        for (int probe = bit; probe < 30; ++probe) {
            if ((fixed & (u32{1} << probe)) == 0) ++available;
        }
        if (available < need) break;
        enumerate_completions(fixed, need - 1, bit + 1, additions | flag,
                              callback);
    }
}

struct GeometryTransform {
    std::vector<i128> width;
    std::vector<i128> cells;
    std::vector<i128> adjacency;
    u64 width_ops = 0;
    u64 cell_ops = 0;
    u64 adjacency_ops = 0;
};

GeometryTransform transform_geometry(const std::vector<Cell>& arrangement,
                                     const ChooseTable& choose) {
    std::map<u32, Atom> atoms;
    std::map<u32, i128> boundary_atoms;
    for (std::size_t i = 0; i < arrangement.size(); ++i) {
        const Cell& cell = arrangement[i];
        if (std::popcount(cell.failed) <= 8) {
            atoms[cell.failed].width += cell.right - cell.left;
            atoms[cell.failed].cells += 1;
        }
        const u32 both = cell.failed |
                         arrangement[(i + 1) % arrangement.size()].failed;
        if (std::popcount(both) <= 8) boundary_atoms[both] += 1;
    }

    GeometryTransform result;
    result.width.assign(REPAIR_COUNT, 0);
    result.cells.assign(REPAIR_COUNT, 0);
    result.adjacency.assign(REPAIR_COUNT, 0);
    for (const auto& [fixed, atom] : atoms) {
        const int need = 8 - std::popcount(fixed);
        enumerate_completions(fixed, need, 0, 0, [&](u32 repair) {
            const u64 rank = choose.rank8(repair);
            result.width[rank] += atom.width;
            result.cells[rank] += atom.cells;
            ++result.width_ops;
            ++result.cell_ops;
        });
    }
    for (const auto& [fixed, count] : boundary_atoms) {
        const int need = 8 - std::popcount(fixed);
        enumerate_completions(fixed, need, 0, 0, [&](u32 repair) {
            result.adjacency[choose.rank8(repair)] += count;
            ++result.adjacency_ops;
        });
    }
    return result;
}

struct ThresholdRepair {
    u64 order = 0;
    u32 mask = 0;
    i64 threshold = 0;
};

struct ThresholdAtlas {
    std::vector<ThresholdRepair> positive;
    u64 geometry_controls = 0;
};

ThresholdAtlas build_threshold_atlas(const std::vector<Cell>& arrangement,
                                     const ChooseTable& choose,
                                     const GeometryTransform& transform) {
    ThresholdAtlas result;
    result.positive.reserve(REPAIR_COUNT);
    u32 repair = (u32{1} << 8) - 1;
    const u32 limit = u32{1} << 30;
    u64 rank = 0;
    while (repair != 0 && repair < limit) {
        ensure(choose.rank8(repair) == rank, "colex traversal/rank mismatch");
        const i128 width = transform.width[rank];
        const i128 components = transform.cells[rank] - transform.adjacency[rank];
        ensure(width > 0 && components > 0, "nonpositive repair geometry");

        // Periodic direct-integration bound, written on the theorem's
        // D*N^2 denominator: g*A - 63*c*D*DeltaC >= 0.
        const i128 slope = static_cast<i128>(63) * N * S * width
                         - static_cast<i128>(4) * D * N * N;
        const i128 loss = static_cast<i128>(63) * components * D * DELTA_C;
        if (slope > 0) {
            const i128 threshold = (loss + slope - 1) / slope;
            ensure(threshold > 0 && threshold <= std::numeric_limits<i64>::max(),
                   "robust threshold overflow");
            result.positive.push_back({
                splitmix64(static_cast<u64>(repair) ^ ORDER_SEED), repair,
                static_cast<i64>(threshold)});
        }

        // A direct cell-block reconstruction is a hostile control on the
        // incidence transform.  The deterministic congruence gives 143
        // controls spread across the full colex universe.
        if (rank % 41001 == 17 || rank + 1 == REPAIR_COUNT) {
            const Geometry direct = direct_geometry(repair, arrangement);
            ensure(static_cast<i128>(direct.width) == width,
                   "direct/transformed width disagreement");
            ensure(static_cast<i128>(direct.components.size()) == components,
                   "direct/transformed component disagreement");
            ++result.geometry_controls;
        }

        ++rank;
        const u32 following = next_mask(repair);
        if (following <= repair) break;
        repair = following;
    }
    ensure(rank == REPAIR_COUNT, "threshold repair universe incomplete");
    std::sort(result.positive.begin(), result.positive.end(),
              [](const ThresholdRepair& a, const ThresholdRepair& b) {
        return a.order != b.order ? a.order < b.order : a.mask < b.mask;
    });
    return result;
}

std::vector<u32> deck_at(const std::vector<ThresholdRepair>& atlas, i64 scale) {
    std::vector<u32> answer;
    for (const ThresholdRepair& item : atlas) {
        if (item.threshold <= scale) answer.push_back(item.mask);
    }
    return answer;
}

void smoke_mode(const std::vector<Cell>& arrangement) {
    Fnv64 ledger;
    u32 repair = (u32{1} << 8) - 1;
    for (int sample = 0; sample < 256; ++sample) {
        const Geometry geometry = direct_geometry(repair, arrangement);
        ledger.word(repair);
        ledger.word(static_cast<u64>(geometry.width));
        ledger.word(geometry.components.size());
        ledger.word(static_cast<u64>(exact_margin(geometry, 1015)));
        ledger.word(static_cast<u64>(exact_margin(geometry, 1016)));
        for (int jump = 0; jump < 97; ++jump) repair = next_mask(repair);
        ensure(repair != 0, "smoke repair traversal ended early");
    }
    std::cout << "LRC_23_RAY_INDEPENDENT_ROBUST_SMOKE_V1\n"
              << "POOL_CELLS " << arrangement.size() << " SAMPLE_REPAIRS 256"
              << " SAMPLE_LEDGER " << hex64(ledger.value()) << "\n"
              << "VERDICT COMPONENT_AND_PRIMITIVE_SMOKE_PASS\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        using namespace ray_audit;
        const bool smoke = argc == 2 && std::string(argv[1]) == "--smoke";
        ensure(argc == 1 || smoke, "usage: independent_robust_tail [--smoke]");
        check_primitive_constants();
        const std::vector<Cell> arrangement = build_pool_arrangement();
        if (smoke) {
            smoke_mode(arrangement);
            return 0;
        }

        const ChooseTable choose;
        GeometryTransform transform = transform_geometry(arrangement, choose);
        ensure(transform.width_ops == UINT64_C(152170690) &&
               transform.cell_ops == UINT64_C(152170690) &&
               transform.adjacency_ops == UINT64_C(112487480),
               "incidence transform operation ledger mismatch");
        ThresholdAtlas atlas = build_threshold_atlas(arrangement, choose, transform);
        ensure(atlas.positive.size() == UINT64_C(5654814),
               "positive-slope census mismatch");
        transform.width.clear();
        transform.cells.clear();
        transform.adjacency.clear();
        transform.width.shrink_to_fit();
        transform.cells.shrink_to_fit();
        transform.adjacency.shrink_to_fit();

        const std::vector<u32> current = deck_at(atlas.positive, 1016);
        const std::vector<u32> prior = deck_at(atlas.positive, 1015);
        ensure(current.size() == 968835, "g=1016 robust deck census mismatch");
        ensure(prior.size() == 965228, "g=1015 robust deck census mismatch");
        Fnv64 deck_ledger;
        for (u32 repair : current) deck_ledger.word(repair);
        ensure(deck_ledger.value() == UINT64_C(0xd809676561e7da8e),
               "robust deck ledger mismatch");

        const std::vector<u32> bodies = all_bodies();
        const BodyScan current_scan = scan_bodies(current, bodies);
        const BodyScan prior_scan = scan_bodies(prior, bodies);
        ensure(current_scan.failures == 0, "robust deck fails at 1016");
        ensure(prior_scan.failures == 1, "robust predecessor failure census mismatch");
        ensure(labels(prior_scan.first_failure) ==
                   "88,95,170,193,240,252,264,286,290",
               "robust predecessor hostile body changed");

        Fnv64 prefix_ledger;
        ensure(current_scan.max_checks <= current.size(), "prefix length overflow");
        for (u64 i = 0; i < current_scan.max_checks; ++i) {
            prefix_ledger.word(current[static_cast<std::size_t>(i)]);
        }
        ensure(prefix_ledger.value() == UINT64_C(0x810f1ca450105a52),
               "robust prefix ledger mismatch");

        std::cout << "LRC_23_RAY_INDEPENDENT_ROBUST_TAIL_V1\n";
        std::cout << "PRIMITIVE_GRID " << N << " SAFE_TICKS " << S
                  << " CENTERED_RANGE " << DELTA_C << " POOL_DENOMINATOR "
                  << D << '\n';
        std::cout << "TRANSFORM_OPS WIDTH " << transform.width_ops
                  << " CELLS " << transform.cell_ops << " ADJACENCY "
                  << transform.adjacency_ops << " DIRECT_GEOMETRY_CONTROLS "
                  << atlas.geometry_controls << '\n';
        std::cout << "POSITIVE_SLOPE_REPAIRS " << atlas.positive.size()
                  << " G1016_REPAIRS " << current.size() << " DECK_FNV "
                  << hex64(deck_ledger.value()) << '\n';
        std::cout << "G1016_SCAN BODIES " << current_scan.bodies
                  << " FAILURES " << current_scan.failures << " CHECKS "
                  << current_scan.checks << " MAX_CHECKS "
                  << current_scan.max_checks << " WORST_BODY {"
                  << labels(current_scan.worst_body) << "}\n";
        std::cout << "G1015_SHARPNESS_ONLY REPAIRS " << prior.size()
                  << " FAILURES " << prior_scan.failures << " CHECKS "
                  << prior_scan.checks << " FIRST_FAILURE {"
                  << labels(prior_scan.first_failure) << "}\n";
        std::cout << "PREFIX_FNV " << hex64(prefix_ledger.value())
                  << " ORDER_SEED " << hex64(ORDER_SEED) << '\n';
        std::cout << "VERDICT EVERY_BODY_UNIFORMLY_CLOSED_FOR_ALL_G_GE_1016\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
