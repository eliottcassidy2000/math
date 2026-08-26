// Primary finite certificate for the sharp THM-4227 scale-separated wedge.
// It reuses the pool-cell path underlying the earlier two-newcomer scout, while the
// maintained independent audit uses a separately implemented fixed geometry.

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

u64 prefix_mix(u64 value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

u32 prefix_label_mask(std::initializer_list<int> values) {
    u32 result = 0;
    for (int value : values) {
        const auto found = std::find(POOL.begin(), POOL.end(), value);
        require(found != POOL.end(), "prefix label absent from pool");
        result |= u32{1} << static_cast<int>(found - POOL.begin());
    }
    return result;
}

struct PrefixCoverScan {
    u64 bodies = 0;
    u64 checks = 0;
    u64 maximum = 0;
    u64 covers = 0;
};

PrefixCoverScan prefix_scan_nine(const std::vector<u32>& edges) {
    PrefixCoverScan result;
    u32 body = (u32{1} << 9) - 1;
    const u32 limit = u32{1} << 30;
    while (body < limit) {
        u64 checked = 0;
        bool missed = false;
        for (u32 edge : edges) {
            ++checked;
            if ((body & edge) == 0) {
                missed = true;
                break;
            }
        }
        ++result.bodies;
        result.checks += checked;
        result.maximum = std::max(result.maximum, checked);
        if (!missed) ++result.covers;
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    return result;
}

struct PrefixEnvelope {
    i128 minimum = 0;
    u32 minimum_edge = 0;
    i64 minimum_components = 0;
    i128 first_q = 0;
    u32 first_q_edge = 0;
    i128 first_r = 0;
    u32 first_r_edge = 0;
    u64 ledger = 0;
};

PrefixEnvelope prefix_envelope(const std::vector<u32>& prefix,
                               const FixedSafeStats& stats) {
    PrefixEnvelope result;
    Fnv1a64 ledger;
    for (u32 edge : prefix) {
        const i64 mass = subset_sum(edge, stats.length);
        const i64 cells = subset_sum(edge, stats.cell_count);
        const i64 adjacencies = subset_sum(edge, stats.adjacency_union_count);
        const i64 components = cells - adjacencies;
        const i128 delta = static_cast<i128>(81) * mass -
                           static_cast<i128>(7) * COMMON;
        require(delta > 0 && components > 0, "invalid prefix edge");
        ledger.add_u64_le(edge);
        ledger.add_u64_le(static_cast<u64>(delta));
        ledger.add_u64_le(static_cast<u64>(components));
        if (result.minimum_edge == 0 || delta < result.minimum) {
            result.minimum = delta;
            result.minimum_edge = edge;
            result.minimum_components = components;
        }
        const i128 numerator = static_cast<i128>(162) * components * COMMON;
        const i128 denominator = static_cast<i128>(7) * delta;
        const i128 threshold = (numerator + denominator - 1) / denominator;
        if (threshold > result.first_q) {
            result.first_q = threshold;
            result.first_q_edge = edge;
        }
    }
    for (u32 edge : prefix) {
        const i64 mass = subset_sum(edge, stats.length);
        const i64 cells = subset_sum(edge, stats.cell_count);
        const i64 adjacencies = subset_sum(edge, stats.adjacency_union_count);
        const i64 components = cells - adjacencies;
        const i128 delta = static_cast<i128>(81) * mass -
                           static_cast<i128>(7) * COMMON;
        const i128 numerator =
            static_cast<i128>(27) * COMMON * (result.first_q + components);
        const i128 threshold = (numerator + delta - 1) / delta;
        if (threshold > result.first_r) {
            result.first_r = threshold;
            result.first_r_edge = edge;
        }
        require(result.minimum * (result.first_q + components) <=
                    delta * (result.first_q + result.minimum_components),
                "minimum edge does not dominate the cofinal envelope");
    }
    result.ledger = ledger.value();
    return result;
}

struct PrimaryQualityEdge {
    i64 delta;
    u64 first_q;
    u32 mask;
    u32 components;
};

struct PrimaryQualityScan {
    u64 bodies = 0;
    u64 checks = 0;
    u64 maximum = 0;
    u64 covers = 0;
    u32 closest = 0;
    u32 missed = 0;
};

PrimaryQualityScan primary_quality_scan(
        const std::vector<PrimaryQualityEdge>& edges) {
    PrimaryQualityScan result;
    u32 body = (u32{1} << 9) - 1;
    const u32 limit = u32{1} << 30;
    while (body < limit) {
        u64 checked = 0;
        u32 missed = 0;
        for (const PrimaryQualityEdge& edge : edges) {
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
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    return result;
}

i128 primary_ceil_div(i128 numerator, i128 denominator) {
    require(numerator >= 0 && denominator > 0, "invalid ceiling division");
    return (numerator + denominator - 1) / denominator;
}

}  // namespace

int main() {
    try {
        const std::vector<Cell> cells = build_pool_cells();
        const FixedSafeStats stats = build_fixed_safe_stats(cells, 8);
        AtomMass atoms;
        atoms.reserve(4096);
        for (const Cell& cell : cells) {
            if (std::popcount(cell.failed_pool) <= 8) {
                atoms[cell.failed_pool] += cell.right - cell.left;
            }
        }

        std::vector<u32> edges;
        edges.reserve(UINT64_C(5267460));
        u64 deletions = 0;
        u64 equalities = 0;
        u32 deletion = (u32{1} << 8) - 1;
        const u32 limit = u32{1} << 30;
        while (deletion != 0 && deletion < limit) {
            ++deletions;
            const i64 mass = subset_sum(deletion, atoms);
            const i128 delta = static_cast<i128>(81) * mass -
                               static_cast<i128>(7) * COMMON;
            if (delta == 0) ++equalities;
            if (delta > 0) edges.push_back(deletion);
            const u32 next = next_combination(deletion);
            if (next <= deletion) break;
            deletion = next;
        }
        require(deletions == UINT64_C(5852925) &&
                    edges.size() == UINT64_C(5267460) && equalities == 0,
                "double-limit layer changed");

        std::vector<PrimaryQualityEdge> quality;
        quality.reserve(edges.size());
        for (u32 edge : edges) {
            const i64 mass = subset_sum(edge, stats.length);
            const i64 cell_count = subset_sum(edge, stats.cell_count);
            const i64 adjacency_count =
                subset_sum(edge, stats.adjacency_union_count);
            const i64 components = cell_count - adjacency_count;
            const i128 delta128 = static_cast<i128>(81) * mass -
                                  static_cast<i128>(7) * COMMON;
            require(delta128 > 0 &&
                        delta128 <= std::numeric_limits<i64>::max() &&
                        components > 0,
                    "invalid quality edge");
            const i128 first_q128 = primary_ceil_div(
                static_cast<i128>(162) * components * COMMON,
                static_cast<i128>(7) * delta128);
            require(first_q128 <= std::numeric_limits<u64>::max(),
                    "quality threshold overflow");
            quality.push_back({static_cast<i64>(delta128),
                               static_cast<u64>(first_q128), edge,
                               static_cast<u32>(components)});
        }
        std::sort(quality.begin(), quality.end(),
                  [](const PrimaryQualityEdge& left,
                     const PrimaryQualityEdge& right) {
            if (left.first_q != right.first_q) return left.first_q < right.first_q;
            if (left.delta != right.delta) return left.delta > right.delta;
            return left.mask < right.mask;
        });
        const PrimaryQualityScan quality_scan = primary_quality_scan(quality);
        require(quality_scan.bodies == UINT64_C(14307150) &&
                    quality_scan.checks == UINT64_C(45855551044) &&
                    quality_scan.maximum == UINT64_C(859307) &&
                    quality_scan.covers == 0 &&
                    quality_scan.closest == prefix_label_mask(
                        {88,95,170,193,240,252,264,286,290}) &&
                    quality_scan.missed == prefix_label_mask(
                        {8,16,85,132,143,145,168,176}),
                "quality-prefix no-cover certificate changed");
        quality.resize(static_cast<std::size_t>(quality_scan.maximum));
        require(quality.back().first_q == UINT64_C(3391),
                "quality-prefix activation changed");
        Fnv1a64 quality_ledger;
        quality_ledger.add_u64_le(quality.size());
        for (const PrimaryQualityEdge& edge : quality) {
            quality_ledger.add_u64_le(edge.mask);
            quality_ledger.add_u64_le(static_cast<u64>(edge.delta));
            quality_ledger.add_u64_le(edge.components);
            quality_ledger.add_u64_le(edge.first_q);
        }
        require(quality_ledger.value() == UINT64_C(0x232b38cfc255acfb),
                "quality-prefix semantic ledger changed");
        const u32 dominant_mask = prefix_label_mask(
            {60,84,170,240,252,264,286,290});
        const PrimaryQualityEdge* dominant = nullptr;
        for (const PrimaryQualityEdge& edge : quality) {
            if (edge.mask == dominant_mask) dominant = &edge;
        }
        require(dominant != nullptr &&
                    dominant->delta == INT64_C(16269324968430) &&
                    dominant->components == 130 && dominant->first_q == 3374,
                "quality dominant edge changed");
        const i128 quality_q = 3391;
        for (const PrimaryQualityEdge& edge : quality) {
            require(edge.delta >= dominant->delta &&
                    static_cast<i128>(dominant->delta) *
                        (quality_q + edge.components) <=
                    static_cast<i128>(edge.delta) *
                        (quality_q + dominant->components),
                    "quality linear envelope changed");
        }
        const i64 quality_divisor = std::gcd(
            static_cast<i64>(27 * COMMON), dominant->delta);
        const i64 quality_numerator = 27 * COMMON / quality_divisor;
        const i64 quality_denominator = dominant->delta / quality_divisor;
        require(quality_numerator == INT64_C(321902813232) &&
                    quality_denominator == INT64_C(10633545731),
                "quality tail ratio changed");
        const i128 quality_half_r = primary_ceil_div(
            static_cast<i128>(27) * COMMON *
                (quality_q + dominant->components), dominant->delta);
        i128 quality_optimized_r = 0;
        u32 quality_optimized_edge = 0;
        for (const PrimaryQualityEdge& edge : quality) {
            const i128 denominator = static_cast<i128>(2) *
                (static_cast<i128>(7) * edge.delta * quality_q -
                 static_cast<i128>(81) * edge.components * COMMON);
            const i128 numerator = static_cast<i128>(189) * COMMON *
                quality_q * (quality_q + edge.components);
            const i128 threshold = primary_ceil_div(numerator, denominator);
            if (threshold > quality_optimized_r) {
                quality_optimized_r = threshold;
                quality_optimized_edge = edge.mask;
            }
        }
        require(quality_half_r == static_cast<i128>(106590) &&
                    quality_optimized_r == static_cast<i128>(106033) &&
                    quality_optimized_edge == dominant_mask,
                "quality corner changed");

        std::sort(edges.begin(), edges.end(), [](u32 left, u32 right) {
            const u64 lhs = prefix_mix(left), rhs = prefix_mix(right);
            return lhs == rhs ? left < right : lhs < rhs;
        });
        edges.resize(549);
        const PrefixCoverScan scan = prefix_scan_nine(edges);
        require(scan.bodies == UINT64_C(14307150) &&
                    scan.checks == UINT64_C(407411762) &&
                    scan.maximum == 549 && scan.covers == 0,
                "prefix no-cover certificate changed");

        const PrefixEnvelope envelope = prefix_envelope(edges, stats);
        const u32 expected_edge = prefix_label_mask(
            {30,40,42,126,170,190,193,264});
        require(envelope.minimum == static_cast<i128>(11314922220LL) &&
                    envelope.minimum_edge == expected_edge &&
                    envelope.minimum_components == 170 &&
                    envelope.first_q == static_cast<i128>(6342592) &&
                    envelope.first_q_edge == expected_edge &&
                    envelope.first_r == static_cast<i128>(276085148833LL) &&
                    envelope.first_r_edge == expected_edge &&
                    envelope.ledger == UINT64_C(0x4257f9e48be9706d),
                "prefix envelope changed");

        const i64 divisor = std::gcd(
            static_cast<i64>(27 * COMMON),
            static_cast<i64>(envelope.minimum));
        require(27 * COMMON / divisor == INT64_C(17019288) &&
                    static_cast<i64>(envelope.minimum) / divisor == 391,
                "reduced tail ratio changed");

        std::cout << "THM-4227 primary prefix certificate\n";
        std::cout << "deletions=" << deletions
                  << " strict_edges=" << UINT64_C(5267460)
                  << " equalities=" << equalities << '\n';
        std::cout << "prefix_edges=" << edges.size()
                  << " ledger=" << hex64(envelope.ledger)
                  << " bodies=" << scan.bodies
                  << " checks=" << scan.checks
                  << " max=" << scan.maximum
                  << " covers=" << scan.covers << '\n';
        std::cout << "min_num=" << decimal(envelope.minimum)
                  << " min_edge=" << labels(envelope.minimum_edge)
                  << " components=" << envelope.minimum_components << '\n';
        std::cout << "first_q=" << decimal(envelope.first_q)
                  << " first_r=" << decimal(envelope.first_r)
                  << " tail=17019288/391*(q+170)\n";
        std::cout << "quality_prefix_edges=" << quality.size()
                  << " ledger=" << hex64(quality_ledger.value())
                  << " bodies=" << quality_scan.bodies
                  << " checks=" << quality_scan.checks
                  << " max=" << quality_scan.maximum
                  << " closest=" << labels(quality_scan.closest)
                  << " missed=" << labels(quality_scan.missed)
                  << " covers=" << quality_scan.covers << '\n';
        std::cout << "quality_first_q=" << decimal(quality_q)
                  << " dominant=" << labels(dominant->mask)
                  << " delta=" << dominant->delta
                  << " components=" << dominant->components
                  << " tail=" << quality_numerator << '/'
                  << quality_denominator << "*(q+"
                  << dominant->components << ")"
                  << " half_corner=" << decimal(quality_half_r)
                  << " optimized_corner=" << decimal(quality_optimized_r)
                  << '\n';
        std::cout << "checks=PASS pool_cells,all_eight_deletions,"
                     "prefix_no_cover,component_slack_envelope,"
                     "quality_prefix\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}
