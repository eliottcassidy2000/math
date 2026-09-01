// Exact audit for the THM-4302 C596 carrier on the protected
// consumer K596 union the complete residual endpoint-595 layer.  This file
// is intentionally not a canonical artifact.  Pair identities and literal
// grids remain attached to every sign computation (MISTAKE-532).

#define ENDPOINT617_RAW_VERIFY_MAIN protected391_hidden_main
#include "04-computation/lrc14_size_preserving_response_staircase_thm4300/endpoint617_exchanged_carrier_raw_verify.cpp"
#undef ENDPOINT617_RAW_VERIFY_MAIN

#include <array>
#include <fstream>
#include <map>

namespace {

constexpr u64 kRepairFnv391 = UINT64_C(0x64ce5f9d1ec8c4c2);
constexpr u64 kAdditionFnv391 = UINT64_C(0xdc0eebaebf688c65);
constexpr u64 kDeleteFnv391 = UINT64_C(0x9240b264ab65aa62);
constexpr u64 kAugmentedFnv391 = UINT64_C(0x55e8588798885ae5);
constexpr u64 kCarrierFnv391 = UINT64_C(0x892fef44a9e6b37e);
constexpr u64 kOldUnionFnv391 = UINT64_C(0x11414a33ab91fef6);
constexpr u64 kThm4302UnionFnv391 = UINT64_C(0xb1c8ecf1dd4a71c5);
constexpr u64 kThm4302ResidualFnv391 = UINT64_C(0x7da11cd038486887);
constexpr u64 kK596Fnv391 = UINT64_C(0xfbf4e6bd7a593649);
constexpr u64 kTop28Fnv391 = UINT64_C(0x47981ce64825ef2a);

using PairKey391 = std::pair<int, int>;

std::set<PairKey391> read_pairs391(const std::filesystem::path& path,
                                   std::size_t expected_count,
                                   u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pair ledger");
    std::set<PairKey391> rows;
    Fnv ledger;
    std::string line;
    PairKey391 previous{0, 0};
    bool have_previous = false;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank pair row");
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed pair row");
        const PairKey391 pair{std::stoi(line.substr(0, comma)),
                              std::stoi(line.substr(comma + 1))};
        require(pair.first > 0 && pair.first < pair.second,
                "invalid pair row");
        require(!have_previous || previous < pair,
                "pair ledger order changed");
        require(rows.insert(pair).second, "duplicate pair row");
        previous = pair;
        have_previous = true;
        ledger.add(pair.first);
        ledger.add(pair.second);
    }
    require(input.eof() && rows.size() == expected_count &&
                ledger.state == expected_fnv,
            "pair ledger identity changed");
    return rows;
}

u64 pair_fnv391(const auto& rows) {
    Fnv ledger;
    for (const auto& pair : rows) {
        ledger.add(pair.first);
        ledger.add(pair.second);
    }
    return ledger.state;
}

std::set<PairKey391> set_union391(const std::set<PairKey391>& left,
                                  const std::set<PairKey391>& right) {
    std::set<PairKey391> answer = left;
    answer.insert(right.begin(), right.end());
    return answer;
}

std::set<PairKey391> set_difference391(const std::set<PairKey391>& left,
                                       const std::set<PairKey391>& right) {
    std::set<PairKey391> answer;
    std::set_difference(left.begin(), left.end(), right.begin(), right.end(),
                        std::inserter(answer, answer.end()));
    return answer;
}

u64 mask_fnv391(const auto& masks) {
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

// Retain every literal-wall failure class.  The carrier itself has ranks
// eight and nine, but the same audit also evaluates a higher-rank hostile
// control, so the rank-at-most-nine truncation is not sufficient here.
Geometry build_full_geometry391(int pair_q, int pair_r) {
    i64 grid = fixed_grid();
    grid = checked_lcm(grid, 14LL * pair_q);
    grid = checked_lcm(grid, 14LL * pair_r);
    std::vector<i64> walls = {0, grid};
    auto add_walls = [&](int speed) {
        const i64 divisor = 14LL * speed;
        require(grid % divisor == 0, "nonintegral wall unit");
        const i64 unit = grid / divisor;
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(pair_q);
    add_walls(pair_r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::map<u32, i64> by_failure;
    Geometry geometry;
    geometry.grid = grid;
    geometry.cells = walls.size() - 1;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(pair_q, grid, left, right) ||
            !safe_midpoint(pair_r, grid, left, right))
            continue;
        u32 failure = 0;
        for (unsigned label = 0; label < kPool.size(); ++label)
            if (!safe_midpoint(kPool[label], grid, left, right))
                failure |= u32{1} << label;
        const i64 width = right - left;
        geometry.pair_ticks += width;
        by_failure[failure] += width;
    }
    for (const auto& entry : by_failure) geometry.classes.push_back(entry);
    return geometry;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 13,
                "usage: protected391 JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION1624 REPAIRS76 ADDITIONS4 DELETE73 ROW_CSV "
                "MASK_CSV CANDIDATE_CSV");

        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> repairs =
            read_mixed617(argv[7], 76, kRepairFnv391);
        const std::vector<u32> additions =
            read_mixed617(argv[8], 4, kAdditionFnv391);
        const std::vector<u32> deletions =
            read_mixed617(argv[9], 73, kDeleteFnv391);

        std::vector<u32> augmented =
            build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(augmented.begin(), augmented.end());
        for (u32 repair : repairs) {
            require(distinct.insert(repair).second, "repair overlap");
            augmented.push_back(repair);
        }
        require(augmented.size() == 9088 &&
                    masks_fnv_agent(augmented) == kAugmentedFnv391,
                "augmented carrier changed");
        const std::set<u32> deletion_set(deletions.begin(), deletions.end());
        std::vector<u32> carrier;
        for (u32 mask : augmented)
            if (!deletion_set.contains(mask)) carrier.push_back(mask);
        for (u32 addition : additions) {
            require(distinct.insert(addition).second, "addition overlap");
            carrier.push_back(addition);
        }
        require(carrier.size() == 9019 &&
                    masks_fnv_agent(carrier) == kCarrierFnv391,
                "C596 identity changed");
        const u64 carrier_rank8 = std::count_if(
            carrier.begin(), carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(carrier_rank8 == 8961, "C596 ranks changed");
        for (u32 mask : joint)
            require(std::find(carrier.begin(), carrier.end(), mask) !=
                        carrier.end(),
                    "joint mask absent from C596");

        const std::set<PairKey391> universe =
            read_pairs391(argv[5], 22647, kResidualFnvAgent);
        const std::set<PairKey391> old_union =
            read_pairs391(argv[6], 1624, kOldUnionFnv391);
        require(std::includes(universe.begin(), universe.end(),
                              old_union.begin(), old_union.end()),
                "old union escaped universe");
        std::set<PairKey391> k596;
        for (PairKey391 pair : universe)
            if (pair.second >= 596) k596.insert(pair);
        require(k596.size() == 363 && pair_fnv391(k596) == kK596Fnv391,
                "K596 identity changed");
        const std::set<PairKey391> thm4302_union =
            set_union391(old_union, k596);
        const std::set<PairKey391> thm4302_residual =
            set_difference391(universe, thm4302_union);
        require(thm4302_union.size() == 1633 &&
                    pair_fnv391(thm4302_union) == kThm4302UnionFnv391 &&
                    thm4302_residual.size() == 21014 &&
                    pair_fnv391(thm4302_residual) == kThm4302ResidualFnv391,
                "THM4302 typed partition changed");
        int maximum_endpoint = 0;
        for (PairKey391 pair : thm4302_residual)
            maximum_endpoint = std::max(maximum_endpoint, pair.second);
        std::set<PairKey391> top28;
        for (PairKey391 pair : thm4302_residual)
            if (pair.second == maximum_endpoint) top28.insert(pair);
        require(maximum_endpoint == 595 && top28.size() == 28 &&
                    pair_fnv391(top28) == kTop28Fnv391,
                "top-595 identity changed");
        const std::set<PairKey391> protected_set =
            set_union391(k596, top28);
        require(protected_set.size() == 391, "protected union size changed");

        // Evaluation order is endpoint descending, then q ascending, matching
        // the maintained endpoint scans.  Canonical pair-set FNV is also
        // recorded separately in q,r order.
        std::vector<PairKey391> protected_rows(protected_set.begin(),
                                               protected_set.end());
        std::sort(protected_rows.begin(), protected_rows.end(),
                  [](PairKey391 left, PairKey391 right) {
                      if (left.second != right.second)
                          return left.second > right.second;
                      return left.first < right.first;
                  });
        Fnv protected_order_ledger;
        for (PairKey391 pair : protected_rows) {
            protected_order_ledger.add(pair.first);
            protected_order_ledger.add(pair.second);
        }

        std::vector<Geometry> geometries;
        geometries.reserve(protected_rows.size());
        std::ofstream row_out(argv[10]);
        require(static_cast<bool>(row_out), "cannot create row grid ledger");
        row_out << "q,r,grid,cells,pair_ticks,classes\n";
        Fnv grid_ledger;
        i64 minimum_grid = std::numeric_limits<i64>::max();
        i64 maximum_grid = 0;
        u64 total_cells = 0;
        u64 total_classes = 0;
        for (PairKey391 pair : protected_rows) {
            Geometry geometry =
                build_full_geometry391(pair.first, pair.second);
            minimum_grid = std::min(minimum_grid, geometry.grid);
            maximum_grid = std::max(maximum_grid, geometry.grid);
            total_cells += geometry.cells;
            total_classes += geometry.classes.size();
            grid_ledger.add(pair.first);
            grid_ledger.add(pair.second);
            grid_ledger.add(static_cast<u64>(geometry.grid));
            grid_ledger.add(geometry.cells);
            grid_ledger.add(static_cast<u64>(geometry.pair_ticks));
            grid_ledger.add(geometry.classes.size());
            row_out << pair.first << ',' << pair.second << ','
                    << geometry.grid << ',' << geometry.cells << ','
                    << geometry.pair_ticks << ',' << geometry.classes.size()
                    << '\n';
            geometries.push_back(std::move(geometry));
        }
        require(row_out.good(), "row grid ledger write failed");

        constexpr u32 candidate = UINT32_C(0x00910bf6);
        require(std::popcount(candidate) == 12 &&
                    std::find(carrier.begin(), carrier.end(), candidate) ==
                        carrier.end(),
                "rank-twelve candidate identity changed");
        std::ofstream candidate_out(argv[12]);
        require(static_cast<bool>(candidate_out),
                "cannot create candidate activity ledger");
        candidate_out << "q,r,grid,mass,ticks,active\n";
        Fnv candidate_ledger;
        u64 candidate_active_rows = 0;
        u64 candidate_equalities = 0;
        unsigned candidate_failure_signature = 0;
        std::map<PairKey391, Margin> candidate_failure_margins;
        for (std::size_t index = 0; index < protected_rows.size(); ++index) {
            const PairKey391 pair = protected_rows[index];
            const Geometry& geometry = geometries[index];
            const Margin exact = margin(geometry, candidate);
            const bool active = exact.ticks >= 0;
            candidate_active_rows += active;
            candidate_equalities += exact.ticks == 0;
            if (pair == PairKey391{96, 595} && active)
                candidate_failure_signature |= 1u;
            if (pair == PairKey391{100, 595} && active)
                candidate_failure_signature |= 2u;
            if (pair == PairKey391{210, 595} && active)
                candidate_failure_signature |= 4u;
            if (pair == PairKey391{96, 595} ||
                pair == PairKey391{100, 595} ||
                pair == PairKey391{210, 595})
                candidate_failure_margins.emplace(pair, exact);
            candidate_ledger.add(candidate);
            candidate_ledger.add(pair.first);
            candidate_ledger.add(pair.second);
            candidate_ledger.add(static_cast<u64>(geometry.grid));
            candidate_ledger.add(static_cast<u64>(exact.mass));
            add_i128(candidate_ledger, exact.ticks);
            candidate_ledger.add(active);
            candidate_out << pair.first << ',' << pair.second << ','
                          << geometry.grid << ',' << exact.mass << ','
                          << decimal(exact.ticks) << ',' << active << '\n';
        }
        require(candidate_out.good(), "candidate activity ledger write failed");

        std::ofstream mask_out(argv[11]);
        require(static_cast<bool>(mask_out), "cannot create mask ledger");
        mask_out << "mask_hex,rank,joint,active_rows,equality_rows,"
                    "failure_pair_signature\n";
        Fnv bool_sign_ledger;
        Fnv exact_margin_ledger;
        std::vector<u32> common_inactive;
        std::vector<u32> inactive_on_k596;
        std::array<u64, 392> active_row_histogram{};
        std::array<u64, 8> failure_signature_histogram{};
        std::map<PairKey391, u64> active_per_row;
        std::map<PairKey391, u64> equal_per_row;
        u64 equality_cells = 0;
        for (u32 mask : carrier) {
            u64 active_rows = 0;
            u64 equality_rows = 0;
            bool active_k596 = false;
            unsigned failure_signature = 0;
            for (std::size_t index = 0; index < protected_rows.size(); ++index) {
                const PairKey391 pair = protected_rows[index];
                const Margin exact = margin(geometries[index], mask);
                const bool active = exact.ticks >= 0;
                const bool equality = exact.ticks == 0;
                active_rows += active;
                equality_rows += equality;
                equality_cells += equality;
                if (pair.second >= 596) active_k596 |= active;
                if (pair == PairKey391{96, 595} && active)
                    failure_signature |= 1u;
                if (pair == PairKey391{100, 595} && active)
                    failure_signature |= 2u;
                if (pair == PairKey391{210, 595} && active)
                    failure_signature |= 4u;
                active_per_row[pair] += active;
                equal_per_row[pair] += equality;
                bool_sign_ledger.add(mask);
                bool_sign_ledger.add(pair.first);
                bool_sign_ledger.add(pair.second);
                bool_sign_ledger.add(static_cast<u64>(geometries[index].grid));
                bool_sign_ledger.add(active);
                exact_margin_ledger.add(mask);
                exact_margin_ledger.add(pair.first);
                exact_margin_ledger.add(pair.second);
                exact_margin_ledger.add(static_cast<u64>(geometries[index].grid));
                add_i128(exact_margin_ledger, exact.ticks);
            }
            ++active_row_histogram[active_rows];
            ++failure_signature_histogram[failure_signature];
            if (active_rows == 0) common_inactive.push_back(mask);
            if (!active_k596) inactive_on_k596.push_back(mask);
            mask_out << hex8(mask) << ',' << std::popcount(mask) << ','
                     << joint_set.contains(mask) << ',' << active_rows << ','
                     << equality_rows << ',' << failure_signature << '\n';
        }
        require(mask_out.good(), "mask activity ledger write failed");
        std::sort(common_inactive.begin(), common_inactive.end());
        std::sort(inactive_on_k596.begin(), inactive_on_k596.end());
        const u64 inactive_rank8 = std::count_if(
            common_inactive.begin(), common_inactive.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        const u64 inactive_rank9 = common_inactive.size() - inactive_rank8;
        const u64 inactive_joint = std::count_if(
            common_inactive.begin(), common_inactive.end(),
            [&](u32 mask) { return joint_set.contains(mask); });

        Fnv row_activity_ledger;
        for (PairKey391 pair : protected_rows) {
            row_activity_ledger.add(pair.first);
            row_activity_ledger.add(pair.second);
            row_activity_ledger.add(active_per_row.at(pair));
            row_activity_ledger.add(equal_per_row.at(pair));
        }

        std::cout << "LRC14_C596_PROTECTED391_ACTIVITY_AUDIT_V1\n"
                  << "CARRIER 9019 FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " RANK8 "
                  << carrier_rank8 << " RANK9 "
                  << carrier.size() - carrier_rank8
                  << " JOINT_RETAINED 421\n"
                  << "ROW_CONSUMER K596 363 FNV " << std::hex
                  << pair_fnv391(k596) << " TOP595 28 FNV "
                  << pair_fnv391(top28) << " PROTECTED 391 CANONICAL_FNV "
                  << pair_fnv391(protected_set) << " EVALUATION_ORDER_FNV "
                  << protected_order_ledger.state << std::dec << '\n'
                  << "GRID_LEDGER FNV " << std::hex << grid_ledger.state
                  << std::dec << " MIN " << minimum_grid << " MAX "
                  << maximum_grid << " TOTAL_CELLS " << total_cells
                  << " TOTAL_CLASSES " << total_classes << '\n'
                  << "SIGN_CELLS " << carrier.size() * protected_rows.size()
                  << " BOOL_SIGN_FNV " << std::hex << bool_sign_ledger.state
                  << " EXACT_MARGIN_FNV " << exact_margin_ledger.state
                  << " ROW_ACTIVITY_FNV " << row_activity_ledger.state
                  << std::dec << " EQUALITIES " << equality_cells << '\n'
                  << "COMMON_INACTIVE " << common_inactive.size() << " FNV "
                  << std::hex << mask_fnv391(common_inactive) << std::dec
                  << " RANK8 " << inactive_rank8 << " RANK9 "
                  << inactive_rank9 << " JOINT " << inactive_joint
                  << " NONJOINT " << common_inactive.size() - inactive_joint
                  << '\n'
                  << "COMMON_INACTIVE_MASKS";
        for (u32 mask : common_inactive)
            std::cout << ' ' << hex8(mask);
        std::cout << '\n'
                  << "K596_INACTIVE_RETAINED " << inactive_on_k596.size()
                  << " FNV " << std::hex << mask_fnv391(inactive_on_k596)
                  << std::dec << " MASKS";
        for (u32 mask : inactive_on_k596)
            std::cout << ' ' << hex8(mask);
        std::cout << '\n' << "ACTIVE_ROW_HISTOGRAM";
        for (std::size_t rows = 0; rows < active_row_histogram.size(); ++rows)
            if (active_row_histogram[rows] != 0)
                std::cout << ' ' << rows << ':' << active_row_histogram[rows];
        std::cout << '\n' << "FAILURE_PAIR_SIGNATURE_HISTOGRAM";
        for (unsigned signature = 0; signature < 8; ++signature)
            if (failure_signature_histogram[signature] != 0)
                std::cout << ' ' << signature << ':'
                          << failure_signature_histogram[signature];
        std::cout << '\n' << "FAILURE_PAIR_ACTIVE_COUNTS"
                  << " 96,595:" << active_per_row.at({96, 595})
                  << " 100,595:" << active_per_row.at({100, 595})
                  << " 210,595:" << active_per_row.at({210, 595}) << '\n'
                  << "RANK12_CANDIDATE 00910bf6 ACTIVE_ROWS "
                  << candidate_active_rows << " EQUALITIES "
                  << candidate_equalities << " FAILURE_SIGNATURE "
                  << candidate_failure_signature << " EXACT_LEDGER_FNV "
                  << std::hex << candidate_ledger.state << std::dec << '\n'
                  << "RANK12_FAILURE_MARGINS";
        for (PairKey391 pair : std::array<PairKey391, 3>{
                 PairKey391{96,595}, PairKey391{100,595},
                 PairKey391{210,595}}) {
            const Margin exact = candidate_failure_margins.at(pair);
            const Geometry geometry = build_geometry(pair.first, pair.second);
            std::cout << ' ' << pair.first << ',' << pair.second
                      << ":grid=" << geometry.grid << ":mass=" << exact.mass
                      << ":ticks=" << decimal(exact.ticks);
        }
        std::cout << '\n'
                  << "DELETION_CAPACITY_STRICT_COMMON_INACTIVE "
                  << common_inactive.size() << '\n'
                  << "SCOPE FINITE_EXACT_FIXED_C596_CARRIER_PROTECTED_"
                     "391_ROW_CONSUMER_ONLY_PAIR_GRIDS_RETAINED_NO_CROSS_"
                     "GRID_MARGIN_ORDER_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS COMPLETE_CARRIER_ACTIVITY_AUDIT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "PROTECTED391_ACTIVITY_AUDIT_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
