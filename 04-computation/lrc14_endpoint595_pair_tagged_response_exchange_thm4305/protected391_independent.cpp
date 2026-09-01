// Structurally independent scratch replay of the protected-391 activity
// audit.  It consumes THM-4303's self-contained carrier reconstruction and
// literal-wall geometry rather than the maintained endpoint primitives.

#define main protected391_independent_hidden_main
#include "04-computation/lrc14_endpoint595_twentyfive_closure_thm4303/endpoint595_independent_replay.cpp"
#undef main

#include <array>

namespace {

using u128_local = __uint128_t;

void add_i128_local(Fnv& ledger, i128 value) {
    const u128_local bits = static_cast<u128_local>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

u64 mask_fnv_local(const std::vector<u32>& masks) {
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

Geometry build_full_geometry_local(Pair pair) {
    i64 grid = fixed_grid();
    grid = checked_lcm(grid, 14LL * pair.q);
    grid = checked_lcm(grid, 14LL * pair.r);
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
    add_walls(pair.q);
    add_walls(pair.r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::map<u32, i64> by_failure;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(pair.q, grid, left, right) ||
            !safe_midpoint(pair.r, grid, left, right))
            continue;
        u32 failure = 0;
        for (unsigned label = 0; label < kPool.size(); ++label)
            if (!safe_midpoint(kPool[label], grid, left, right))
                failure |= u32{1} << label;
        by_failure[failure] += right - left;
    }
    Geometry geometry;
    geometry.grid = grid;
    for (const auto& entry : by_failure) geometry.classes.push_back(entry);
    return geometry;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 10,
                "usage: independent391 JOINT BASE8951 ADD45 SUFFIX9 "
                "REPAIRS76 DELETE73 ADDITIONS4 UNIVERSE OLD_UNION1624");
        const std::vector<u32> joint = read_masks(
            argv[1], 421, UINT64_C(0x20d63dd42fe8150e), 8, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> carrier = reconstruct_carrier(
            argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);

        std::set<Pair> universe;
        std::set<Pair> typed_union;
        const std::vector<Pair> top_vector = derive_top_rows(
            argv[8], argv[9], universe, typed_union);
        const std::set<Pair> top(top_vector.begin(), top_vector.end());
        std::set<Pair> k596;
        for (Pair pair : universe)
            if (pair.r >= 596) k596.insert(pair);
        require(k596.size() == 363 &&
                    pair_fingerprint(k596) == UINT64_C(0xfbf4e6bd7a593649),
                "K596 identity changed");
        std::set<Pair> protected_set = k596;
        protected_set.insert(top.begin(), top.end());
        require(protected_set.size() == 391 &&
                    pair_fingerprint(protected_set) ==
                        UINT64_C(0xc732a1532c12b9f6),
                "protected set identity changed");
        std::vector<Pair> rows(protected_set.begin(), protected_set.end());
        std::sort(rows.begin(), rows.end(), [](Pair left, Pair right) {
            if (left.r != right.r) return left.r > right.r;
            return left.q < right.q;
        });
        Fnv order_ledger;
        for (Pair pair : rows) {
            order_ledger.add(pair.q);
            order_ledger.add(pair.r);
        }
        require(order_ledger.state == UINT64_C(0x38e567eb72a62c3e),
                "evaluation order changed");

        std::vector<Geometry> geometries;
        geometries.reserve(rows.size());
        Fnv grid_class_ledger;
        for (Pair pair : rows) {
            Geometry geometry = build_full_geometry_local(pair);
            grid_class_ledger.add(pair.q);
            grid_class_ledger.add(pair.r);
            grid_class_ledger.add(static_cast<u64>(geometry.grid));
            grid_class_ledger.add(geometry.classes.size());
            for (const auto& [failure, width] : geometry.classes) {
                grid_class_ledger.add(failure);
                grid_class_ledger.add(static_cast<u64>(width));
            }
            geometries.push_back(std::move(geometry));
        }

        constexpr u32 candidate = UINT32_C(0x00910bf6);
        require(std::popcount(candidate) == 12 &&
                    std::find(carrier.begin(), carrier.end(), candidate) ==
                        carrier.end(),
                "rank-twelve candidate identity changed");
        Fnv candidate_ledger;
        u64 candidate_active_rows = 0;
        u64 candidate_equalities = 0;
        std::map<Pair, i128> candidate_failure_ticks;
        constexpr u32 remainder_cap = UINT32_C(0x00000022);
        std::map<Pair, i128> remainder_cap_ticks;
        for (std::size_t index = 0; index < rows.size(); ++index) {
            const Pair pair = rows[index];
            const Geometry& geometry = geometries[index];
            i64 mass = 0;
            for (const auto& [failure, width] : geometry.classes)
                if ((failure & ~candidate) == 0) mass += width;
            const i128 ticks = static_cast<i128>(63) * mass -
                               static_cast<i128>(4) * geometry.grid;
            const bool active = ticks >= 0;
            candidate_active_rows += active;
            candidate_equalities += ticks == 0;
            if (pair == Pair{96,595} || pair == Pair{100,595} ||
                pair == Pair{210,595}) {
                candidate_failure_ticks.emplace(pair, ticks);
                i64 cap_mass = 0;
                for (const auto& [failure, width] : geometry.classes)
                    if ((failure & ~remainder_cap) == 0) cap_mass += width;
                remainder_cap_ticks.emplace(
                    pair, static_cast<i128>(63) * cap_mass -
                              static_cast<i128>(4) * geometry.grid);
            }
            candidate_ledger.add(candidate);
            candidate_ledger.add(pair.q);
            candidate_ledger.add(pair.r);
            candidate_ledger.add(static_cast<u64>(geometry.grid));
            candidate_ledger.add(static_cast<u64>(mass));
            add_i128_local(candidate_ledger, ticks);
            candidate_ledger.add(active);
        }
        require(candidate_active_rows == 391 && candidate_equalities == 0 &&
                    candidate_ledger.state == UINT64_C(0xfe3bbad41d0c9ce3) &&
                    candidate_failure_ticks.at({96,595}) ==
                        static_cast<i128>(3303345791970LL) &&
                    candidate_failure_ticks.at({100,595}) ==
                        static_cast<i128>(31968167963256LL) &&
                    candidate_failure_ticks.at({210,595}) ==
                        static_cast<i128>(8209878640866LL) &&
                    remainder_cap_ticks.at({96,595}) ==
                        static_cast<i128>(-58664013767604LL) &&
                    remainder_cap_ticks.at({100,595}) ==
                        static_cast<i128>(-142736095694304LL) &&
                    remainder_cap_ticks.at({210,595}) ==
                        static_cast<i128>(-28403370492336LL),
                "rank-twelve candidate activity changed");

        Fnv bool_sign_ledger;
        Fnv exact_margin_ledger;
        std::vector<u32> common_inactive;
        std::vector<u32> k596_inactive;
        std::array<u64, 8> signature_histogram{};
        std::map<Pair, u64> active_per_row;
        std::map<Pair, u64> equality_per_row;
        u64 equality_cells = 0;
        for (u32 mask : carrier) {
            u64 active_rows = 0;
            bool ever_active_k596 = false;
            unsigned signature = 0;
            for (std::size_t index = 0; index < rows.size(); ++index) {
                const Pair pair = rows[index];
                const Geometry& geometry = geometries[index];
                i64 mass = 0;
                for (const auto& [failure, width] : geometry.classes)
                    if ((failure & ~mask) == 0) mass += width;
                const i128 ticks = static_cast<i128>(63) * mass -
                                   static_cast<i128>(4) * geometry.grid;
                const bool active = ticks >= 0;
                const bool equality = ticks == 0;
                active_rows += active;
                equality_cells += equality;
                ever_active_k596 |= pair.r >= 596 && active;
                if (pair == Pair{96, 595} && active) signature |= 1u;
                if (pair == Pair{100, 595} && active) signature |= 2u;
                if (pair == Pair{210, 595} && active) signature |= 4u;
                active_per_row[pair] += active;
                equality_per_row[pair] += equality;
                bool_sign_ledger.add(mask);
                bool_sign_ledger.add(pair.q);
                bool_sign_ledger.add(pair.r);
                bool_sign_ledger.add(static_cast<u64>(geometry.grid));
                bool_sign_ledger.add(active);
                exact_margin_ledger.add(mask);
                exact_margin_ledger.add(pair.q);
                exact_margin_ledger.add(pair.r);
                exact_margin_ledger.add(static_cast<u64>(geometry.grid));
                add_i128_local(exact_margin_ledger, ticks);
            }
            ++signature_histogram[signature];
            if (active_rows == 0) common_inactive.push_back(mask);
            if (!ever_active_k596) k596_inactive.push_back(mask);
        }
        std::sort(common_inactive.begin(), common_inactive.end());
        std::sort(k596_inactive.begin(), k596_inactive.end());
        const std::vector<u32> expected_inactive = {
            UINT32_C(0x380086a0), UINT32_C(0x388088c0)};
        require(common_inactive == expected_inactive &&
                    k596_inactive == expected_inactive &&
                    mask_fnv_local(common_inactive) ==
                        UINT64_C(0x3b5ca775eedae38b),
                "common inactive pool changed");
        require(equality_cells == 0 &&
                    bool_sign_ledger.state == UINT64_C(0x855c6f82deef952f) &&
                    exact_margin_ledger.state == UINT64_C(0xeee2c2c38d5b2e55),
                "protected sign ledger changed");
        require(signature_histogram ==
                    std::array<u64, 8>{4405,303,288,296,891,253,783,1800},
                "failure signature histogram changed");
        require(active_per_row.at({96,595}) == 2652 &&
                    active_per_row.at({100,595}) == 3167 &&
                    active_per_row.at({210,595}) == 3727,
                "failure-row activity counts changed");

        Fnv row_activity_ledger;
        for (Pair pair : rows) {
            row_activity_ledger.add(pair.q);
            row_activity_ledger.add(pair.r);
            row_activity_ledger.add(active_per_row.at(pair));
            row_activity_ledger.add(equality_per_row.at(pair));
        }
        require(row_activity_ledger.state ==
                    UINT64_C(0x61549efd32b7ed3b),
                "row activity ledger changed");

        std::cout << "LRC14_C596_PROTECTED391_INDEPENDENT_AUDIT_V1\n"
                  << "CARRIER " << carrier.size() << " FNV "
                  << hex16(fingerprint(carrier)) << " RANK8 8961 RANK9 58\n"
                  << "PROTECTED_ROWS " << rows.size() << " CANONICAL_FNV "
                  << hex16(pair_fingerprint(protected_set))
                  << " ORDER_FNV " << hex16(order_ledger.state) << '\n'
                  << "GRID_CLASS_LEDGER_FNV "
                  << hex16(grid_class_ledger.state) << '\n'
                  << "SIGN_CELLS " << carrier.size() * rows.size()
                  << " BOOL_SIGN_FNV " << hex16(bool_sign_ledger.state)
                  << " EXACT_MARGIN_FNV "
                  << hex16(exact_margin_ledger.state)
                  << " ROW_ACTIVITY_FNV "
                  << hex16(row_activity_ledger.state)
                  << " EQUALITIES " << equality_cells << '\n'
                  << "COMMON_INACTIVE 2 FNV "
                  << hex16(mask_fnv_local(common_inactive))
                  << " RANK8 2 RANK9 0 JOINT 0 NONJOINT 2 MASKS"
                  << " 380086a0 388088c0\n"
                  << "FAILURE_PAIR_ACTIVITY 96,595:2652 100,595:3167 "
                     "210,595:3727\n"
                  << "RANK12_CANDIDATE 00910bf6 ACTIVE_ROWS "
                  << candidate_active_rows << " EQUALITIES "
                  << candidate_equalities << " EXACT_LEDGER_FNV "
                  << hex16(candidate_ledger.state) << '\n'
                  << "REMAINDER_CAP 00000022 RANK 2 TICKS "
                     "96,595:-58664013767604 "
                     "100,595:-142736095694304 "
                     "210,595:-28403370492336 ALL_INACTIVE\n"
                  << "VERDICT PASS STRUCTURALLY_INDEPENDENT_PROTECTED391_"
                     "ACTIVITY_REPLAY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "PROTECTED391_INDEPENDENT_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
