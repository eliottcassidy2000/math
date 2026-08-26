// FINITE-EXACT scout for the first post-full-pool lane: two fixed newcomers
// and nine pool labels. Geometry is rebuilt on the literal joint wall lattice.

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main thm4188_canonical_main
#include "lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp"
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace {

i64 pair_lcm(i64 left, i64 right) {
    const i64 common = std::gcd(left, right);
    const i128 product = static_cast<i128>(left / common) * right;
    require(product <= std::numeric_limits<i64>::max(), "pair lattice overflow");
    return static_cast<i64>(product);
}

bool safe_mid(int speed, i64 left, i64 right, i64 denominator) {
    const i128 twice = static_cast<i128>(2) * denominator;
    i128 raw = static_cast<i128>(speed) * (left + right);
    i128 residue = raw % twice;
    if (residue < 0) residue += twice;
    return static_cast<i128>(7) * residue >= denominator &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * denominator;
}

struct PairMass {
    i64 denominator = 0;
    AtomMass atoms;
    u64 cells = 0;
};

PairMass build_pair_mass(int q1, int q2, int max_arity) {
    require(q1 > 0 && q2 > 0 && q1 != q2, "newcomers must be positive/distinct");
    require(std::find(POOL.begin(), POOL.end(), q1) == POOL.end(),
            "q1 lies in pool");
    require(std::find(POOL.begin(), POOL.end(), q2) == POOL.end(),
            "q2 lies in pool");
    PairMass result;
    result.denominator = pair_lcm(pair_lcm(COMMON, 14LL * q1), 14LL * q2);
    std::vector<i64> walls = {0, result.denominator};
    auto add_walls = [&](int speed) {
        const i64 unit = result.denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : POOL) add_walls(speed);
    add_walls(q1);
    add_walls(q2);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    result.cells = walls.size() - 1;
    result.atoms.reserve(4096);
    i128 total_pair_safe = 0;
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        const i64 left = walls[index];
        const i64 right = walls[index + 1];
        if (!safe_mid(q1, left, right, result.denominator) ||
            !safe_mid(q2, left, right, result.denominator)) continue;
        total_pair_safe += right - left;
        u32 failed = 0;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if (!safe_mid(POOL[vertex], left, right, result.denominator)) {
                failed |= u32{1} << vertex;
            }
        }
        if (std::popcount(failed) <= max_arity) {
            result.atoms[failed] += right - left;
        }
    }
    require(total_pair_safe > 0 && total_pair_safe < result.denominator,
            "pair-safe control is degenerate");
    return result;
}

EdgeLayer build_pair_layer(const PairMass& mass, int arity) {
    EdgeLayer layer{arity, {}, 0};
    const u32 limit = u32{1} << 30;
    u32 deletion = (u32{1} << arity) - 1;
    while (deletion != 0 && deletion < limit) {
        const i64 good_mass = deletion_mass(deletion, mass.atoms);
        const i128 comparison = static_cast<i128>(63) * good_mass -
                                static_cast<i128>(4) * mass.denominator;
        if (comparison == 0) ++layer.equalities;
        if (comparison >= 0) layer.edges.push_back(deletion);
        const u32 next = next_combination(deletion);
        if (next <= deletion) break;
        deletion = next;
    }
    return layer;
}

u64 mix64(u64 value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

struct CoverScan {
    u64 bodies = 0;
    u64 checks = 0;
    u64 maximum = 0;
    u32 closest = 0;
    u32 missed = 0;
    u32 cover = 0;
    u64 cover_count = 0;
};

CoverScan scan_nine(std::vector<u32> edges) {
    std::sort(edges.begin(), edges.end(), [](u32 left, u32 right) {
        const u64 lh = mix64(left), rh = mix64(right);
        return lh == rh ? left < right : lh < rh;
    });
    CoverScan result;
    u32 body = (u32{1} << 9) - 1;
    const u32 limit = u32{1} << 30;
    while (body < limit) {
        u64 checked = 0;
        u32 missed = 0;
        for (u32 edge : edges) {
            ++checked;
            if ((body & edge) == 0) {
                missed = edge;
                break;
            }
        }
        ++result.bodies;
        result.checks += checked;
        if (checked > result.maximum) {
            result.maximum = checked;
            result.closest = body;
            result.missed = missed;
        }
        if (missed == 0) {
            if (result.cover == 0) result.cover = body;
            ++result.cover_count;
        }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    return result;
}

void probe(int q1, int q2, int arity) {
    // Retain every observed pool-failure atom.  The shallow layer only reads
    // submasks of its deletion, while the full support lets us measure a
    // discovered nine-cover body against its 21-label complement exactly.
    const PairMass mass = build_pair_mass(q1, q2, 30);
    const EdgeLayer layer = build_pair_layer(mass, arity);
    const CoverScan scan = scan_nine(layer.edges);
    require(scan.bodies == UINT64_C(14307150),
            "nine-body universe changed");
    i128 cover_delta = 0;
    if (scan.cover != 0) {
        const u32 all_pool = (u32{1} << 30) - 1;
        const u32 complement = all_pool ^ scan.cover;
        const i64 target_mass = deletion_mass(complement, mass.atoms);
        cover_delta = static_cast<i128>(63) * target_mass -
                      static_cast<i128>(4) * mass.denominator;
    }
    if (q1 == 49 && q2 == 50 && arity == 7) {
        require(layer.edges.size() == 145462 && layer.equalities == 0 &&
                    scan.cover_count == 120,
                "49,50 depth-seven control changed");
    } else if (q1 == 49 && q2 == 50 && arity == 8) {
        require(layer.edges.size() == 1536023 && layer.equalities == 0 &&
                    scan.cover_count == 0,
                "49,50 depth-eight control changed");
    } else if (q1 == 6 && q2 == 50 && arity == 8) {
        require(layer.edges.size() == 13497 && layer.equalities == 0 &&
                    scan.cover_count == 472050 &&
                    cover_delta == static_cast<i128>(601044495065784LL),
                "6,50 depth-eight hostile changed");
    }
    std::cout << "PAIR " << q1 << ',' << q2
              << " D " << arity
              << " DENOM " << mass.denominator
              << " CELLS " << mass.cells
              << " EDGES " << layer.edges.size()
              << " EQUALITIES " << layer.equalities
              << " BODIES " << scan.bodies
              << " CHECKS " << scan.checks
              << " MAX " << scan.maximum
              << " CLOSEST " << labels(scan.closest)
              << " MISSED " << labels(scan.missed)
              << " COVER " << labels(scan.cover)
              << " COVER_COUNT " << scan.cover_count
              << " COVER_BODY_DELTA " << decimal(cover_delta) << '\n';
}

void probe_double_limit() {
    const auto cells = build_pool_cells();
    AtomMass fixed;
    fixed.reserve(4096);
    for (const Cell& cell : cells) {
        if (std::popcount(cell.failed_pool) <= 8) {
            fixed[cell.failed_pool] += cell.right - cell.left;
        }
    }
    EdgeLayer layer{8, {}, 0};
    i128 minimum_strict_numerator = 0;
    u32 minimum_strict_edge = 0;
    const u32 limit = u32{1} << 30;
    u32 deletion = (u32{1} << 8) - 1;
    while (deletion != 0 && deletion < limit) {
        const i64 good_mass = deletion_mass(deletion, fixed);
        // Two asymptotically independent safe combs retain density (6/7)^2.
        // (36/49) M > 4/63 is exactly 81 M > 7.
        const i128 comparison = static_cast<i128>(81) * good_mass -
                                static_cast<i128>(7) * COMMON;
        if (comparison == 0) ++layer.equalities;
        if (comparison > 0) {
            layer.edges.push_back(deletion);
            if (minimum_strict_edge == 0 || comparison < minimum_strict_numerator) {
                minimum_strict_numerator = comparison;
                minimum_strict_edge = deletion;
            }
        }
        const u32 next = next_combination(deletion);
        if (next <= deletion) break;
        deletion = next;
    }
    const CoverScan scan = scan_nine(layer.edges);
    require(layer.edges.size() == 5267460 && layer.equalities == 0 &&
                minimum_strict_numerator == static_cast<i128>(944928) &&
                scan.bodies == UINT64_C(14307150) && scan.cover_count == 0,
            "strict double-limit control changed");
    std::cout << "DOUBLE_LIMIT D 8 EDGES " << layer.edges.size()
              << " EQUALITIES " << layer.equalities
              << " MIN_STRICT_NUM " << decimal(minimum_strict_numerator)
              << " MIN_STRICT_EDGE " << labels(minimum_strict_edge)
              << " BODIES " << scan.bodies
              << " CHECKS " << scan.checks
              << " MAX " << scan.maximum
              << " CLOSEST " << labels(scan.closest)
              << " MISSED " << labels(scan.missed)
              << " COVER " << labels(scan.cover)
              << " COVER_COUNT " << scan.cover_count << '\n';
    const i128 component_cap = 7133;
    const i128 q1_numerator = static_cast<i128>(162) * component_cap * COMMON;
    const i128 q1_denominator = static_cast<i128>(7) * minimum_strict_numerator;
    const i128 q1_cutoff =
        (q1_numerator + q1_denominator - 1) / q1_denominator;
    const i128 q2_numerator = static_cast<i128>(27) * COMMON;
    const i128 q2_multiplier =
        (q2_numerator + minimum_strict_numerator - 1) /
        minimum_strict_numerator;
    require(q1_cutoff == static_cast<i128>(3186712759230LL) &&
                q2_multiplier == static_cast<i128>(521215695),
            "scale-separated cutoff changed");
    std::cout << "SCALE_SEPARATED_WEDGE Q1_AT_LEAST "
              << decimal(q1_cutoff)
              << " Q2_AT_LEAST CEIL(" << decimal(q2_multiplier)
              << "*(7133+Q1)) OR SWAP" << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "all") {
        probe(49, 50, 7);
        probe(49, 50, 8);
        probe(6, 50, 8);
        probe_double_limit();
        return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "limit") {
        probe_double_limit();
        return 0;
    }
    require(argc == 4,
            "usage: probe_two_newcomers all|limit|q1 q2 deletion_arity");
    probe(std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]));
}
