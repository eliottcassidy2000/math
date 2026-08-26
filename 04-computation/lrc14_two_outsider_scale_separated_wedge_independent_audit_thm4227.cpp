// Independent exact audit for THM-4227.
//
// The primary scout inherits the THM-4188 pool-cell implementation and uses
// Gosper masks.  This audit instead inherits the explicit fixed-wall geometry
// from the independent THM-4188 joint-wall referee, rebuilds the base masses,
// recursively enumerates both subset universes, and uses a different edge
// ordering in the no-nine-cover scan.

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main thm4188_independent_main
#include "lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.cpp"
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace {

struct BaseSupport {
    std::vector<u32> masks;
    std::unordered_map<u32, u32> position;
    std::vector<i64> atom;
    std::unordered_map<u32, i64> cell_count;
    std::unordered_map<u32, i64> adjacency_count;
};

BaseSupport build_base_support(const FixedGeometry& fixed) {
    BaseSupport result;
    for (const FixedCell& cell : fixed.cells) {
        if (std::popcount(cell.failed) <= 8) result.masks.push_back(cell.failed);
    }
    std::sort(result.masks.begin(), result.masks.end());
    result.masks.erase(std::unique(result.masks.begin(), result.masks.end()),
                       result.masks.end());
    result.position.reserve(2 * result.masks.size());
    result.atom.assign(result.masks.size(), 0);
    for (u32 i = 0; i < result.masks.size(); ++i) {
        result.position.emplace(result.masks[i], i);
    }
    for (const FixedCell& cell : fixed.cells) {
        const auto found = result.position.find(cell.failed);
        if (found != result.position.end()) {
            result.atom[found->second] += cell.right - cell.left;
            ++result.cell_count[cell.failed];
        }
    }
    for (std::size_t i = 0; i < fixed.cells.size(); ++i) {
        const u32 joined = fixed.cells[i].failed |
            fixed.cells[(i + 1) % fixed.cells.size()].failed;
        if (std::popcount(joined) <= 8) ++result.adjacency_count[joined];
    }
    return result;
}

i64 base_mass(u32 deletion, const BaseSupport& support) {
    i64 result = 0;
    u32 subset = deletion;
    while (true) {
        const auto found = support.position.find(subset);
        if (found != support.position.end()) result += support.atom[found->second];
        if (subset == 0) break;
        subset = (subset - 1) & deletion;
    }
    return result;
}

i64 sparse_subset_sum(u32 deletion,
                      const std::unordered_map<u32, i64>& values) {
    i64 result = 0;
    u32 subset = deletion;
    while (true) {
        const auto found = values.find(subset);
        if (found != values.end()) result += found->second;
        if (subset == 0) break;
        subset = (subset - 1) & deletion;
    }
    return result;
}

i64 base_components(u32 deletion, const BaseSupport& support) {
    const i64 cells = sparse_subset_sum(deletion, support.cell_count);
    const i64 adjacencies = sparse_subset_sum(deletion, support.adjacency_count);
    require(cells > adjacencies, "invalid base component count");
    return cells - adjacencies;
}

u64 independent_mix(u64 value) {
    value ^= UINT64_C(0xd1b54a32d192ed03);
    value *= UINT64_C(0x9e3779b185ebca87);
    value ^= value >> 29;
    value *= UINT64_C(0xc2b2ae3d27d4eb4f);
    return value ^ (value >> 32);
}

struct BodyScan {
    u64 bodies = 0;
    u64 checks = 0;
    u64 maximum = 0;
    u64 covers = 0;
    u32 closest = 0;
    u32 missed = 0;
};

BodyScan scan_bodies(std::vector<u32> edges) {
    std::sort(edges.begin(), edges.end(), [](u32 left, u32 right) {
        const u64 lhs = independent_mix(left);
        const u64 rhs = independent_mix(right);
        return lhs == rhs ? left < right : lhs < rhs;
    });
    BodyScan result;
    for_each_k_subset(30, 9, [&](u32 body) {
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
        if (missed == 0) ++result.covers;
    });
    return result;
}

u64 primary_mix(u64 value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

BodyScan scan_ordered_prefix(const std::vector<u32>& edges) {
    BodyScan result;
    for_each_k_subset(30, 9, [&](u32 body) {
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
        if (missed == 0) ++result.covers;
    });
    return result;
}

i64 nonnull_components(const FixedGeometry& fixed, u32 deletion) {
    std::vector<bool> safe;
    safe.reserve(fixed.cells.size());
    for (const FixedCell& cell : fixed.cells) {
        safe.push_back((cell.failed & ~deletion) == 0);
    }
    const u64 safe_count = static_cast<u64>(
        std::count(safe.begin(), safe.end(), true));
    require(safe_count > 0, "base safe set is empty");
    if (safe_count == safe.size()) return 1;
    i64 components = 0;
    for (std::size_t i = 0; i < safe.size(); ++i) {
        const std::size_t previous = (i + safe.size() - 1) % safe.size();
        if (safe[i] && !safe[previous]) ++components;
    }
    require(components > 0, "component census failed");
    return components;
}

struct QualityEdge {
    i64 delta;
    u64 first_q;
    u32 mask;
    u32 components;
};

BodyScan scan_quality_edges(const std::vector<QualityEdge>& edges) {
    BodyScan result;
    for_each_k_subset(30, 9, [&](u32 body) {
        u64 checked = 0;
        u32 missed = 0;
        for (const QualityEdge& edge : edges) {
            ++checked;
            if ((body & edge.mask) == 0) {
                missed = edge.mask;
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
        if (missed == 0) ++result.covers;
    });
    return result;
}

i128 ceil_div(i128 numerator, i128 denominator) {
    require(numerator >= 0 && denominator > 0, "invalid ceiling division");
    return (numerator + denominator - 1) / denominator;
}

}  // namespace

int main() {
    try {
        const FixedGeometry fixed = build_fixed_geometry();
        const BaseSupport support = build_base_support(fixed);
        std::vector<u32> edges;
        edges.reserve(UINT64_C(5267460));
        u64 deletions = 0;
        u64 equalities = 0;
        i128 minimum = 0;
        u32 minimum_edge = 0;
        Ledger ledger;

        for_each_k_subset(30, 8, [&](u32 deletion) {
            ++deletions;
            const i64 mass = base_mass(deletion, support);
            const i128 delta = static_cast<i128>(81) * mass -
                               static_cast<i128>(7) * POOL_DENOMINATOR;
            ledger.add(deletion);
            ledger.add(static_cast<u64>(delta));
            if (delta == 0) ++equalities;
            if (delta > 0) {
                edges.push_back(deletion);
                if (minimum_edge == 0 || delta < minimum) {
                    minimum = delta;
                    minimum_edge = deletion;
                }
            }
        });

        require(fixed.walls.size() == UINT64_C(7134) &&
                    fixed.cells.size() == UINT64_C(7133) &&
                    support.masks.size() == UINT64_C(2718),
                "fixed geometry changed");
        require(deletions == UINT64_C(5852925),
                "eight-deletion universe changed");
        require(edges.size() == UINT64_C(5267460) && equalities == 0,
                "double-limit layer changed");
        require(ledger.value() == UINT64_C(0xed44e240af3d3a59),
                "double-limit semantic ledger changed");
        require(minimum == static_cast<i128>(944928) &&
                    minimum_edge == label_mask(
                        {15,16,20,40,170,190,193,240}),
                "minimum strict double-limit edge changed");

        const BodyScan scan = scan_bodies(edges);
        require(scan.bodies == UINT64_C(14307150) && scan.covers == 0,
                "a nine-body covers the double-limit layer");
        require(scan.checks == UINT64_C(413986509) && scan.maximum == 501 &&
                    scan.closest == label_mask(
                        {16,42,60,88,95,132,143,170,176}) &&
                    scan.missed == label_mask(
                        {8,15,20,63,84,85,120,240}),
                "independent no-cover traversal changed");

        std::vector<u32> prefix = edges;
        std::sort(prefix.begin(), prefix.end(), [](u32 left, u32 right) {
            const u64 lhs = primary_mix(left);
            const u64 rhs = primary_mix(right);
            return lhs == rhs ? left < right : lhs < rhs;
        });
        prefix.resize(549);
        const BodyScan prefix_scan = scan_ordered_prefix(prefix);
        require(prefix_scan.bodies == UINT64_C(14307150) &&
                    prefix_scan.checks == UINT64_C(407411762) &&
                    prefix_scan.maximum == 549 && prefix_scan.covers == 0 &&
                    prefix_scan.closest == label_mask(
                        {16,42,60,85,88,145,168,193,252}) &&
                    prefix_scan.missed == label_mask(
                        {8,20,40,80,132,143,176,286}),
                "primary 549-edge prefix no-cover certificate changed");

        i128 prefix_minimum = 0;
        u32 prefix_minimum_edge = 0;
        i64 prefix_minimum_components = 0;
        i128 prefix_q = 0;
        u32 prefix_q_edge = 0;
        Ledger prefix_ledger;
        for (u32 edge : prefix) {
            const i64 mass = base_mass(edge, support);
            const i128 delta = static_cast<i128>(81) * mass -
                               static_cast<i128>(7) * POOL_DENOMINATOR;
            const i64 components = nonnull_components(fixed, edge);
            require(delta > 0, "prefix contains a non-strict edge");
            prefix_ledger.add(edge);
            prefix_ledger.add(static_cast<u64>(delta));
            prefix_ledger.add(static_cast<u64>(components));
            if (prefix_minimum_edge == 0 || delta < prefix_minimum) {
                prefix_minimum = delta;
                prefix_minimum_edge = edge;
                prefix_minimum_components = components;
            }
            const i128 numerator =
                static_cast<i128>(162) * components * POOL_DENOMINATOR;
            const i128 denominator = static_cast<i128>(7) * delta;
            const i128 threshold =
                (numerator + denominator - 1) / denominator;
            if (threshold > prefix_q) {
                prefix_q = threshold;
                prefix_q_edge = edge;
            }
        }
        require(prefix_minimum == static_cast<i128>(11314922220LL) &&
                    prefix_minimum_edge == label_mask(
                        {30,40,42,126,170,190,193,264}) &&
                    prefix_minimum_components == 170,
                "prefix minimum changed");
        require(prefix_q == static_cast<i128>(6342592),
                "prefix first-stage threshold changed");

        i128 prefix_r = 0;
        u32 prefix_r_edge = 0;
        for (u32 edge : prefix) {
            const i64 mass = base_mass(edge, support);
            const i128 delta = static_cast<i128>(81) * mass -
                               static_cast<i128>(7) * POOL_DENOMINATOR;
            const i64 components = nonnull_components(fixed, edge);
            const i128 numerator = static_cast<i128>(27) * POOL_DENOMINATOR *
                                   (prefix_q + components);
            const i128 threshold = (numerator + delta - 1) / delta;
            if (threshold > prefix_r) {
                prefix_r = threshold;
                prefix_r_edge = edge;
            }
            require(prefix_minimum * (prefix_q + components) <=
                        delta * (prefix_q + prefix_minimum_components),
                    "prefix minimum edge does not dominate the tail envelope");
        }
        require(prefix_r == static_cast<i128>(276085148833LL) &&
                    prefix_r_edge == prefix_minimum_edge,
                "prefix second-stage threshold changed");
        const i64 prefix_ratio_gcd = std::gcd(
            static_cast<i64>(27 * POOL_DENOMINATOR),
            static_cast<i64>(prefix_minimum));
        const i64 prefix_ratio_numerator =
            static_cast<i64>(27 * POOL_DENOMINATOR) / prefix_ratio_gcd;
        const i64 prefix_ratio_denominator =
            static_cast<i64>(prefix_minimum) / prefix_ratio_gcd;
        require(prefix_ratio_numerator == INT64_C(17019288) &&
                    prefix_ratio_denominator == 391,
                "prefix tail ratio changed");

        std::vector<QualityEdge> quality;
        quality.reserve(edges.size());
        for (u32 edge : edges) {
            const i64 mass = base_mass(edge, support);
            const i128 delta128 = static_cast<i128>(81) * mass -
                                  static_cast<i128>(7) * POOL_DENOMINATOR;
            const i64 components = base_components(edge, support);
            require(delta128 > 0 &&
                        delta128 <= std::numeric_limits<i64>::max() &&
                        components > 0,
                    "invalid quality edge");
            const i128 first_q128 = ceil_div(
                static_cast<i128>(162) * components * POOL_DENOMINATOR,
                static_cast<i128>(7) * delta128);
            require(first_q128 <= std::numeric_limits<u64>::max(),
                    "quality threshold overflow");
            quality.push_back({static_cast<i64>(delta128),
                               static_cast<u64>(first_q128), edge,
                               static_cast<u32>(components)});
        }
        std::sort(quality.begin(), quality.end(),
                  [](const QualityEdge& left, const QualityEdge& right) {
            if (left.first_q != right.first_q) return left.first_q < right.first_q;
            if (left.delta != right.delta) return left.delta > right.delta;
            return left.mask < right.mask;
        });
        const BodyScan quality_scan = scan_quality_edges(quality);
        require(quality_scan.bodies == UINT64_C(14307150) &&
                    quality_scan.checks == UINT64_C(45855551044) &&
                    quality_scan.maximum == UINT64_C(859307) &&
                    quality_scan.covers == 0 &&
                    quality_scan.closest == label_mask(
                        {88,95,170,193,240,252,264,286,290}) &&
                    quality_scan.missed == label_mask(
                        {8,16,85,132,143,145,168,176}),
                "quality-prefix no-cover certificate changed");
        quality.resize(static_cast<std::size_t>(quality_scan.maximum));
        require(quality.back().first_q == UINT64_C(3391),
                "quality-prefix activation changed");
        Ledger quality_ledger;
        quality_ledger.add(quality.size());
        for (const QualityEdge& edge : quality) {
            quality_ledger.add(edge.mask);
            quality_ledger.add(static_cast<u64>(edge.delta));
            quality_ledger.add(edge.components);
            quality_ledger.add(edge.first_q);
        }
        require(quality_ledger.value() == UINT64_C(0x232b38cfc255acfb),
                "quality-prefix semantic ledger changed");

        const u32 quality_dominant = label_mask(
            {60,84,170,240,252,264,286,290});
        const u32 quality_runner = label_mask(
            {20,60,168,170,240,252,264,290});
        const QualityEdge* dominant = nullptr;
        const QualityEdge* runner = nullptr;
        for (const QualityEdge& edge : quality) {
            if (edge.mask == quality_dominant) dominant = &edge;
            if (edge.mask == quality_runner) runner = &edge;
        }
        require(dominant != nullptr && runner != nullptr &&
                    dominant->delta == INT64_C(16269324968430) &&
                    dominant->components == 130 && dominant->first_q == 3374 &&
                    runner->delta == INT64_C(16415601867918) &&
                    runner->components == 130 && runner->first_q == 3344,
                "quality envelope controls changed");
        const i128 quality_q = 3391;
        for (const QualityEdge& edge : quality) {
            require(edge.delta >= dominant->delta &&
                    static_cast<i128>(dominant->delta) *
                        (quality_q + edge.components) <=
                    static_cast<i128>(edge.delta) *
                        (quality_q + dominant->components),
                    "dominant quality line failed");
        }
        const i128 runner_gap =
            (quality_q + dominant->components) *
            (static_cast<i128>(runner->delta) - dominant->delta);
        require(runner_gap == static_cast<i128>(515040963097248LL),
                "quality runner gap changed");
        const i64 quality_ratio_gcd = std::gcd(
            static_cast<i64>(27 * POOL_DENOMINATOR), dominant->delta);
        const i64 quality_ratio_numerator =
            static_cast<i64>(27 * POOL_DENOMINATOR) / quality_ratio_gcd;
        const i64 quality_ratio_denominator = dominant->delta / quality_ratio_gcd;
        require(quality_ratio_numerator == INT64_C(321902813232) &&
                    quality_ratio_denominator == INT64_C(10633545731),
                "quality tail ratio changed");
        const i128 quality_half_r = ceil_div(
            static_cast<i128>(27) * POOL_DENOMINATOR *
                (quality_q + dominant->components),
            dominant->delta);
        require(quality_half_r == static_cast<i128>(106590),
                "quality half-slack corner changed");
        i128 quality_optimized_r = 0;
        u32 quality_optimized_edge = 0;
        for (const QualityEdge& edge : quality) {
            const i128 denominator = static_cast<i128>(2) *
                (static_cast<i128>(7) * edge.delta * quality_q -
                 static_cast<i128>(81) * edge.components * POOL_DENOMINATOR);
            const i128 numerator = static_cast<i128>(189) * POOL_DENOMINATOR *
                quality_q * (quality_q + edge.components);
            const i128 threshold = ceil_div(numerator, denominator);
            if (threshold > quality_optimized_r) {
                quality_optimized_r = threshold;
                quality_optimized_edge = edge.mask;
            }
        }
        require(quality_optimized_r == static_cast<i128>(106033) &&
                    quality_optimized_edge == quality_dominant,
                "quality optimized corner changed");

        const i64 epsilon_numerator = 4 * static_cast<i64>(minimum);
        const i64 epsilon_denominator = 441 * POOL_DENOMINATOR;
        const i64 epsilon_gcd = std::gcd(epsilon_numerator, epsilon_denominator);
        const i128 curve_a_numerator =
            static_cast<i128>(81) * 7133 * POOL_DENOMINATOR;
        const i128 curve_a_denominator = static_cast<i128>(7) * minimum;
        const i128 twice_b_numerator =
            static_cast<i128>(27) * POOL_DENOMINATOR;
        require(curve_a_numerator % curve_a_denominator == 0 &&
                    twice_b_numerator % minimum == 0,
                "analytic wedge constants are nonintegral");
        const i128 curve_a = curve_a_numerator / curve_a_denominator;
        const i128 twice_b = twice_b_numerator / minimum;
        require(curve_a == static_cast<i128>(1593356379615LL) &&
                    twice_b == static_cast<i128>(521215695LL),
                "analytic wedge constants changed");

        std::cout << "THM-4227 independent fixed-geometry audit\n";
        std::cout << "fixed_walls=" << fixed.walls.size()
                  << " fixed_cells=" << fixed.cells.size()
                  << " low_support=" << support.masks.size() << '\n';
        std::cout << "eight_deletions=" << deletions
                  << " strict_edges=" << edges.size()
                  << " equalities=" << equalities
                  << " ledger=" << hex64(ledger.value()) << '\n';
        std::cout << "min_strict_num=" << decimal(minimum)
                  << " min_edge=" << labels(minimum_edge) << '\n';
        std::cout << "nine_bodies=" << scan.bodies
                  << " checks=" << scan.checks
                  << " max=" << scan.maximum
                  << " closest=" << labels(scan.closest)
                  << " missed=" << labels(scan.missed)
                  << " covers=" << scan.covers << '\n';
        std::cout << "prefix_edges=" << prefix.size()
                  << " ledger=" << hex64(prefix_ledger.value())
                  << " checks=" << prefix_scan.checks
                  << " max=" << prefix_scan.maximum
                  << " covers=" << prefix_scan.covers << '\n';
        std::cout << "prefix_min_num=" << decimal(prefix_minimum)
                  << " prefix_min_edge=" << labels(prefix_minimum_edge)
                  << " components=" << prefix_minimum_components
                  << " first_q=" << decimal(prefix_q)
                  << " first_q_edge=" << labels(prefix_q_edge) << '\n';
        std::cout << "prefix_tail=" << prefix_ratio_numerator << '/'
                  << prefix_ratio_denominator << "*(q+"
                  << prefix_minimum_components << ")"
                  << " at_first_q=" << decimal(prefix_r) << '\n';
        std::cout << "quality_prefix_edges=" << quality.size()
                  << " ledger=" << hex64(quality_ledger.value())
                  << " checks=" << quality_scan.checks
                  << " max=" << quality_scan.maximum
                  << " covers=" << quality_scan.covers << '\n';
        std::cout << "quality_first_q=" << decimal(quality_q)
                  << " tail=" << quality_ratio_numerator << '/'
                  << quality_ratio_denominator << "*(q+"
                  << dominant->components << ")"
                  << " half_corner=" << decimal(quality_half_r)
                  << " optimized_corner=" << decimal(quality_optimized_r)
                  << " optimized_edge=" << labels(quality_optimized_edge) << '\n';
        std::cout << "epsilon="
                  << epsilon_numerator / epsilon_gcd << '/'
                  << epsilon_denominator / epsilon_gcd
                  << " curve_A=" << decimal(curve_a)
                  << " curve_twice_B=" << decimal(twice_b) << '\n';
        std::cout << "simple_q=" << decimal(2 * curve_a)
                  << " simple_multiplier=" << decimal(twice_b) << '\n';
        std::cout << "checks=PASS independent_geometry,all_eight_deletions,"
                     "strict_minimum,all_nine_bodies,prefix_certificate,"
                     "component_envelope,quality_prefix,analytic_constants\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
