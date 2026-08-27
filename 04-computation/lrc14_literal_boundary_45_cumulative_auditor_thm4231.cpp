#define main all_q_referee_main
#if __has_include("lrc14_all_fixed_outsider_ray_census_independent_audit_thm4231.cpp")
#include "lrc14_all_fixed_outsider_ray_census_independent_audit_thm4231.cpp"
#else
#include "../../04-computation/lrc14_all_fixed_outsider_ray_census_independent_audit_thm4231.cpp"
#endif
#undef main

namespace {

struct JointCoefficients {
    i64 denominator;
    std::vector<i64> mass;
    std::vector<i64> components;
    std::size_t walls;
    i64 raw_mass;
    i64 raw_components;
};

JointCoefficients build_joint_coefficients(
    int q,
    int r,
    const BaseGeometry& base,
    const std::vector<u32>& potential_keys,
    const std::unordered_map<u32, u16>& key_index
) {
    const i64 denominator = std::lcm(
        std::lcm(base.denominator, i64{14} * q), i64{14} * r
    );
    const i64 scale = denominator / base.denominator;
    std::vector<ScaledEvent> events;
    events.reserve(base.events.size() + 2 * q + 2 * r);
    std::vector<i64> walls = {0, denominator};
    walls.reserve(events.capacity() + 2);
    for (const RawEvent& event : base.events) {
        events.push_back({event.position * scale, event.vertex, event.enter_safe});
        walls.push_back(event.position * scale);
    }
    const auto append_outsider = [&](int speed, int vertex) {
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            const i64 enter = (14LL * tooth + 1) * unit;
            const i64 leave = (14LL * tooth + 13) * unit;
            events.push_back({enter, vertex, true});
            events.push_back({leave, vertex, false});
            walls.push_back(enter);
            walls.push_back(leave);
        }
    };
    append_outsider(q, 30);
    append_outsider(r, 31);
    std::sort(events.begin(), events.end(), [](const ScaledEvent& a, const ScaledEvent& b) {
        return std::tie(a.position, a.vertex, a.enter_safe) <
               std::tie(b.position, b.vertex, b.enter_safe);
    });
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    std::vector<i64> mass(potential_keys.size(), 0);
    std::vector<i64> components(potential_keys.size(), 0);
    u32 failure = (u32{1} << 30) - 1;
    bool q_safe = false;
    bool r_safe = false;
    bool previous_joint_safe = false;
    u32 previous_failure = failure;
    i64 raw_mass = 0;
    i64 raw_components = 0;
    std::size_t event_index = 0;
    for (std::size_t wall_index = 0; wall_index + 1 < walls.size(); ++wall_index) {
        const i64 left = walls[wall_index];
        while (event_index < events.size() && events[event_index].position == left) {
            const ScaledEvent& event = events[event_index++];
            if (event.vertex == 30) q_safe = event.enter_safe;
            else if (event.vertex == 31) r_safe = event.enter_safe;
            else if (event.enter_safe) failure &= ~(u32{1} << event.vertex);
            else failure |= u32{1} << event.vertex;
        }
        const bool joint_safe = q_safe && r_safe;
        if (joint_safe) {
            const i64 width = walls[wall_index + 1] - left;
            raw_mass += width;
            ++raw_components;
            if (std::popcount(failure) <= 30 - BODY_WEIGHT) {
                const auto it = key_index.find(failure);
                require(it != key_index.end(), "joint mass key absent");
                mass[it->second] += width;
                components[it->second] += 1;
            }
        }
        if (joint_safe && previous_joint_safe) {
            --raw_components;
            const u32 joined = failure | previous_failure;
            if (std::popcount(joined) <= 30 - BODY_WEIGHT) {
                const auto it = key_index.find(joined);
                require(it != key_index.end(), "joint component key absent");
                components[it->second] -= 1;
            }
        }
        previous_joint_safe = joint_safe;
        previous_failure = failure;
    }
    require(event_index == events.size(), "joint event exhaustion");
    require(raw_mass > 0 && raw_components > 0, "raw joint geometry");
    return {
        denominator, std::move(mass), std::move(components), walls.size(),
        raw_mass, raw_components
    };
}

struct PairResult {
    int q;
    int r;
    i64 denominator;
    std::size_t walls;
    i64 raw_mass;
    i64 raw_components;
    u64 positive;
    u64 zero;
    u64 negative;
    i64 minimum_margin;
    u64 minimum_count;
    u32 witness;
    i64 witness_mass;
    i64 witness_components;
    u64 digest_xor;
    u64 digest_sum;
    u64 ordered_fingerprint;
};

PairResult evaluate_pair(
    int q,
    int r,
    const BaseGeometry& base,
    const MaskSpace& space,
    const std::vector<u32>& potential_keys,
    const std::unordered_map<u32, u16>& key_index,
    const Incidence& incidence,
    const std::array<std::vector<PackedEdge>, 30>& edges,
    Fnv1a64& combined_fingerprint
) {
    const JointCoefficients coefficients = build_joint_coefficients(
        q, r, base, potential_keys, key_index
    );
    std::vector<i64> mass(space.total_masks);
    std::vector<i64> components(space.total_masks);
    for (int weight = 0; weight <= BODY_WEIGHT; ++weight) {
        const u32 begin = space.offsets[weight];
        const u32 end = space.offsets[weight + 1];
        const i64 sign = (weight & 1) ? -1 : 1;
#pragma omp parallel for schedule(static)
        for (i64 raw_index = begin; raw_index < end; ++raw_index) {
            const u32 index = static_cast<u32>(raw_index);
            i64 m = 0;
            i64 c = 0;
            for (u32 position = incidence.offsets[index];
                 position < incidence.offsets[index + 1]; ++position) {
                const u16 key = incidence.keys[position];
                m += coefficients.mass[key];
                c += coefficients.components[key];
            }
            mass[index] = sign * m;
            components[index] = sign * c;
        }
    }
    for (int bit = 0; bit < 30; ++bit) {
        const std::vector<PackedEdge>& list = edges[bit];
#pragma omp parallel for schedule(static)
        for (i64 raw_edge = 0; raw_edge < static_cast<i64>(list.size()); ++raw_edge) {
            const PackedEdge edge = list[raw_edge];
            const u32 source = static_cast<u32>(edge);
            const u32 destination = static_cast<u32>(edge >> 32);
            mass[destination] += mass[source];
            components[destination] += components[source];
        }
    }

    struct Summary {
        u64 positive = 0;
        u64 zero = 0;
        u64 negative = 0;
        i64 minimum_margin = std::numeric_limits<i64>::max();
        u64 minimum_count = 0;
        u32 witness = 0;
        i64 witness_mass = 0;
        i64 witness_components = 0;
        u64 digest_xor = 0;
        u64 digest_sum = 0;
    } total;
    const u32 body_offset = space.offsets[BODY_WEIGHT];
#pragma omp parallel
    {
        Summary local;
#pragma omp for schedule(static)
        for (i64 raw_body = 0; raw_body < static_cast<i64>(space.body_masks.size()); ++raw_body) {
            const u32 body_index = static_cast<u32>(raw_body);
            const u32 body = space.body_masks[body_index];
            const i64 m = mass[body_offset + body_index];
            const i64 c = components[body_offset + body_index];
            require(m > 0 && c > 0, "literal body mass/components");
            const i64 margin = 63 * m - 4 * coefficients.denominator;
            if (margin > 0) ++local.positive;
            else if (margin == 0) ++local.zero;
            else ++local.negative;
            if (margin < local.minimum_margin) {
                local.minimum_margin = margin;
                local.minimum_count = 1;
                local.witness = body;
                local.witness_mass = m;
                local.witness_components = c;
            } else if (margin == local.minimum_margin) {
                ++local.minimum_count;
                if (body < local.witness) {
                    local.witness = body;
                    local.witness_mass = m;
                    local.witness_components = c;
                }
            }
            const u64 word = cheap_mix(
                (u64{static_cast<u32>(q)} << 32) ^
                (u64{static_cast<u32>(r)} << 48) ^ body ^
                std::rotl(static_cast<u64>(m), 13) ^
                std::rotl(static_cast<u64>(c), 39) ^
                std::rotl(static_cast<u64>(margin), 21)
            );
            local.digest_xor ^= word;
            local.digest_sum += word;
        }
#pragma omp critical
        {
            total.positive += local.positive;
            total.zero += local.zero;
            total.negative += local.negative;
            if (local.minimum_margin < total.minimum_margin) {
                total.minimum_margin = local.minimum_margin;
                total.minimum_count = local.minimum_count;
                total.witness = local.witness;
                total.witness_mass = local.witness_mass;
                total.witness_components = local.witness_components;
            } else if (local.minimum_margin == total.minimum_margin) {
                total.minimum_count += local.minimum_count;
                if (local.witness < total.witness) {
                    total.witness = local.witness;
                    total.witness_mass = local.witness_mass;
                    total.witness_components = local.witness_components;
                }
            }
            total.digest_xor ^= local.digest_xor;
            total.digest_sum += local.digest_sum;
        }
    }
    require(total.positive + total.zero + total.negative == BODY_COUNT, "pair universe");

    Fnv1a64 ordered;
    combined_fingerprint.add(q);
    combined_fingerprint.add(r);
    for (u32 body_index = 0; body_index < space.body_masks.size(); ++body_index) {
        const u32 body = space.body_masks[body_index];
        const i64 m = mass[body_offset + body_index];
        const i64 c = components[body_offset + body_index];
        const i64 margin = 63 * m - 4 * coefficients.denominator;
        ordered.add(body);
        ordered.add(static_cast<u64>(m));
        ordered.add(static_cast<u64>(c));
        ordered.add(static_cast<u64>(margin));
        combined_fingerprint.add(body);
        combined_fingerprint.add(static_cast<u64>(m));
        combined_fingerprint.add(static_cast<u64>(c));
        combined_fingerprint.add(static_cast<u64>(margin));
    }
    return {
        q, r, coefficients.denominator, coefficients.walls,
        coefficients.raw_mass, coefficients.raw_components,
        total.positive, total.zero, total.negative,
        total.minimum_margin, total.minimum_count, total.witness,
        total.witness_mass, total.witness_components,
        total.digest_xor, total.digest_sum, ordered.value
    };
}

}  // namespace

int main() {
    const BaseGeometry base = build_base_geometry();
    const MaskSpace space;
    const std::vector<u32> potential_keys = build_potential_keys(base);
    std::unordered_map<u32, u16> key_index;
    key_index.reserve(2 * potential_keys.size());
    for (u32 index = 0; index < potential_keys.size(); ++index) {
        key_index.emplace(potential_keys[index], static_cast<u16>(index));
    }
    const Incidence incidence = build_incidence(potential_keys, space);
    const auto edges = build_transform_edges(space);
    u64 edge_count = 0;
    for (const auto& list : edges) edge_count += list.size();

    Fnv1a64 combined_fingerprint;
    std::vector<PairResult> results;
    const std::array<std::pair<int, int>, 45> pairs = {{
        {744, 824}, {744, 822}, {744, 821}, {744, 820}, {744, 818},
        {744, 817}, {744, 815}, {744, 814}, {744, 813}, {744, 812},
        {744, 811}, {744, 810}, {744, 809}, {744, 805}, {744, 803},
        {744, 800}, {766, 800}, {744, 798}, {744, 794}, {744, 793},
        {744, 791}, {744, 790}, {744, 789}, {744, 787}, {744, 780},
        {765, 780}, {766, 780}, {768, 780}, {616, 777}, {616, 776},
        {616, 775}, {744, 775}, {616, 774}, {616, 773}, {744, 773},
        {616, 772}, {616, 771}, {721, 771}, {616, 770}, {721, 770},
        {744, 770}, {750, 770}, {765, 770}, {766, 770}, {768, 770}
    }};
    for (const auto& [q, r] : pairs) {
        results.push_back(evaluate_pair(
            q, r, base, space, potential_keys, key_index,
            incidence, edges, combined_fingerprint
        ));
    }

    u64 combined_xor = 0;
    u64 combined_sum = 0;
    u64 all_positive = 0;
    for (const PairResult& result : results) {
        combined_xor ^= result.digest_xor;
        combined_sum += result.digest_sum;
        if (result.zero == 0 && result.negative == 0) ++all_positive;
    }
    std::cout << "LRC14_LITERAL_BOUNDARY_45_CUMULATIVE_AUDITOR\n";
    std::cout << "UNIVERSE PAIRS " << results.size()
              << " BODIES_PER_PAIR " << BODY_COUNT
              << " CASES " << results.size() * BODY_COUNT
              << " TRUNCATED_MASKS " << space.total_masks
              << " POTENTIAL_KEYS " << potential_keys.size()
              << " INCIDENCES " << incidence.keys.size()
              << " TRANSFORM_EDGES " << edge_count << '\n';
    for (const PairResult& result : results) {
        const i64 divisor = std::gcd(result.minimum_margin, result.denominator);
        std::cout << "PAIR " << result.q << ',' << result.r
                  << " D " << result.denominator
                  << " WALLS " << result.walls
                  << " RAW_MASS " << result.raw_mass
                  << " RAW_COMPONENTS " << result.raw_components
                  << " POSITIVE " << result.positive
                  << " ZERO " << result.zero
                  << " NEGATIVE " << result.negative << '\n';
        std::cout << "PAIR_MIN " << result.q << ',' << result.r
                  << " MARGIN63_TICKS " << result.minimum_margin
                  << " COUNT " << result.minimum_count
                  << " WITNESS " << labels(result.witness)
                  << " MASS_TICKS " << result.witness_mass
                  << " COMPONENTS " << result.witness_components
                  << " REDUCED " << result.minimum_margin / divisor
                  << '/' << result.denominator / divisor << '\n';
        std::cout << "PAIR_DIGEST " << result.q << ',' << result.r
                  << " XOR " << std::hex << std::setw(16) << std::setfill('0')
                  << result.digest_xor
                  << " SUM " << std::setw(16) << result.digest_sum
                  << " ORDERED_FNV " << std::setw(16) << result.ordered_fingerprint
                  << std::dec << std::setfill(' ') << '\n';
    }
    std::cout << "COMBINED_DIGEST XOR " << std::hex << std::setw(16)
              << std::setfill('0') << combined_xor
              << " SUM " << std::setw(16) << combined_sum
              << " CONCATENATED_FNV " << std::setw(16) << combined_fingerprint.value
              << std::dec << std::setfill(' ') << '\n';
    std::cout << "VERDICT ALL_POSITIVE_PAIRS " << all_positive
              << " OF " << results.size() << " LRC14_OPEN\n";
    return 0;
}
