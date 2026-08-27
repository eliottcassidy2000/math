// Clean-room exact referee for all fixed-q one-tail rays, 1 <= q <= 1306.
// No primary range table or source is imported.
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <omp.h>

using i64 = std::int64_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

namespace {

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
};
constexpr int BODY_WEIGHT = 9;
constexpr int Q_MAXIMUM = 1306;
constexpr u64 BODY_COUNT = 14'307'150;
constexpr u64 EXPECTED_CASES = 18'255'923'400ULL;

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "FAIL " << message << '\n';
    std::exit(1);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

std::string decimal(u128 value) {
    if (value == 0) return "0";
    std::string answer;
    while (value != 0) {
        answer.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    std::reverse(answer.begin(), answer.end());
    return answer;
}

u64 cheap_mix(u64 x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

u64 ceil_div(u128 numerator, u128 denominator) {
    require(denominator != 0, "zero ceiling denominator");
    const u128 answer = (numerator + denominator - 1) / denominator;
    require(answer <= std::numeric_limits<u64>::max(), "ceiling overflow");
    return static_cast<u64>(answer);
}

struct Fnv1a64 {
    u64 value = 0xcbf29ce484222325ULL;
    void add(u64 word) {
        for (int shift = 0; shift < 64; shift += 8) {
            value ^= (word >> shift) & 0xffULL;
            value *= 0x100000001b3ULL;
        }
    }
};

std::string labels(u32 body) {
    std::ostringstream out;
    out << '{';
    bool first = true;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((body & (u32{1} << vertex)) == 0) continue;
        if (!first) out << ',';
        out << POOL[vertex];
        first = false;
    }
    out << '}';
    return out.str();
}

bool is_pool_label(int q) {
    return std::find(POOL.begin(), POOL.end(), q) != POOL.end();
}

struct RawEvent {
    i64 position;
    int vertex;
    bool enter_safe;
};

struct BaseCell {
    i64 left;
    i64 right;
    u32 failure;
};

struct BaseGeometry {
    i64 denominator;
    std::vector<RawEvent> events;
    std::vector<BaseCell> cells;
};

BaseGeometry build_base_geometry() {
    i64 denominator = 1;
    for (int speed : POOL) denominator = std::lcm(denominator, i64{14} * speed);

    std::vector<RawEvent> events;
    std::vector<i64> walls = {0, denominator};
    for (int vertex = 0; vertex < 30; ++vertex) {
        const int speed = POOL[vertex];
        const i64 unit = denominator / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            const i64 enter = (14LL * tooth + 1) * unit;
            const i64 leave = (14LL * tooth + 13) * unit;
            events.push_back({enter, vertex, true});
            events.push_back({leave, vertex, false});
            walls.push_back(enter);
            walls.push_back(leave);
        }
    }
    std::sort(events.begin(), events.end(), [](const RawEvent& a, const RawEvent& b) {
        return std::tie(a.position, a.vertex, a.enter_safe) <
               std::tie(b.position, b.vertex, b.enter_safe);
    });
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    u32 failure = (u32{1} << 30) - 1;
    std::size_t event_index = 0;
    std::vector<BaseCell> cells;
    cells.reserve(walls.size() - 1);
    for (std::size_t i = 0; i + 1 < walls.size(); ++i) {
        const i64 left = walls[i];
        while (event_index < events.size() && events[event_index].position == left) {
            const RawEvent& event = events[event_index++];
            if (event.enter_safe) {
                failure &= ~(u32{1} << event.vertex);
            } else {
                failure |= u32{1} << event.vertex;
            }
        }
        cells.push_back({left, walls[i + 1], failure});
    }
    require(event_index == events.size(), "base event exhaustion");
    require(walls.size() == 7'134 && cells.size() == 7'133, "base arrangement");
    return {denominator, std::move(events), std::move(cells)};
}

struct MaskSpace {
    std::array<std::array<u32, 11>, 31> choose{};
    std::array<u32, BODY_WEIGHT + 2> offsets{};
    std::vector<u32> body_masks;
    u32 total_masks = 0;

    MaskSpace() {
        for (int n = 0; n <= 30; ++n) {
            choose[n][0] = 1;
            for (int k = 1; k <= 10; ++k) {
                if (k > n) choose[n][k] = 0;
                else if (k == n) choose[n][k] = 1;
                else choose[n][k] = choose[n - 1][k - 1] + choose[n - 1][k];
            }
        }
        offsets[0] = 0;
        for (int weight = 0; weight <= BODY_WEIGHT; ++weight) {
            offsets[weight + 1] = offsets[weight] + choose[30][weight];
        }
        total_masks = offsets[BODY_WEIGHT + 1];
        require(total_masks == 22'964'087, "truncated mask universe");

        body_masks.reserve(choose[30][BODY_WEIGHT]);
        u32 mask = (u32{1} << BODY_WEIGHT) - 1;
        const u32 limit = u32{1} << 30;
        while (mask < limit) {
            body_masks.push_back(mask);
            const u32 low = mask & (~mask + 1);
            const u32 high = mask + low;
            mask = high | (((high ^ mask) >> 2) / low);
        }
        require(body_masks.size() == BODY_COUNT, "body mask universe");
    }

    u32 rank(u32 mask, int weight) const {
        u32 answer = 0;
        int ordinal = 1;
        while (mask != 0) {
            const int position = std::countr_zero(mask);
            answer += choose[position][ordinal];
            ++ordinal;
            mask &= mask - 1;
        }
        require(ordinal == weight + 1, "rank weight mismatch");
        return answer;
    }

    u32 index(u32 mask) const {
        const int weight = std::popcount(mask);
        require(weight <= BODY_WEIGHT, "index outside truncated space");
        return offsets[weight] + rank(mask, weight);
    }
};

std::vector<u32> build_potential_keys(const BaseGeometry& base) {
    std::vector<u32> keys;
    keys.reserve(2 * base.cells.size());
    for (std::size_t i = 0; i < base.cells.size(); ++i) {
        const u32 current = base.cells[i].failure;
        const u32 previous = base.cells[(i + base.cells.size() - 1) % base.cells.size()].failure;
        if (std::popcount(current) <= 30 - BODY_WEIGHT) keys.push_back(current);
        const u32 joined = current | previous;
        if (std::popcount(joined) <= 30 - BODY_WEIGHT) keys.push_back(joined);
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    require(keys.size() < std::numeric_limits<u16>::max(), "too many coefficient keys");
    return keys;
}

struct Incidence {
    std::vector<u32> offsets;
    std::vector<u16> keys;
};

template <class Visitor>
void visit_truncated_submasks(u32 mask, const MaskSpace& space, Visitor&& visitor) {
    u32 submask = mask;
    while (true) {
        if (std::popcount(submask) <= BODY_WEIGHT) visitor(space.index(submask));
        if (submask == 0) break;
        submask = (submask - 1) & mask;
    }
}

Incidence build_incidence(const std::vector<u32>& potential_keys, const MaskSpace& space) {
    std::vector<u32> counts(space.total_masks, 0);
    u64 incidence_count = 0;
    for (u32 key : potential_keys) {
        visit_truncated_submasks(key, space, [&](u32 index) {
            require(counts[index] != std::numeric_limits<u32>::max(), "incidence count overflow");
            ++counts[index];
            ++incidence_count;
        });
    }
    require(incidence_count < std::numeric_limits<u32>::max(), "incidence universe overflow");

    std::vector<u32> offsets(space.total_masks + 1, 0);
    for (u32 i = 0; i < space.total_masks; ++i) offsets[i + 1] = offsets[i] + counts[i];
    require(offsets.back() == incidence_count, "incidence prefix");
    std::vector<u16> entries(incidence_count);
    std::vector<u32> cursor = offsets;
    for (u32 key_index = 0; key_index < potential_keys.size(); ++key_index) {
        visit_truncated_submasks(potential_keys[key_index], space, [&](u32 index) {
            entries[cursor[index]++] = static_cast<u16>(key_index);
        });
    }
    return {std::move(offsets), std::move(entries)};
}

using PackedEdge = u64;

std::array<std::vector<PackedEdge>, 30> build_transform_edges(const MaskSpace& space) {
    std::array<std::vector<PackedEdge>, 30> edges;
    u64 per_bit = 0;
    for (int weight = 1; weight <= BODY_WEIGHT; ++weight) {
        per_bit += space.choose[29][weight - 1];
    }
    for (auto& list : edges) list.reserve(per_bit);

    for (int weight = 1; weight <= BODY_WEIGHT; ++weight) {
        u32 mask = (u32{1} << weight) - 1;
        const u32 limit = u32{1} << 30;
        u32 rank = 0;
        while (mask < limit) {
            const u32 destination = space.offsets[weight] + rank;
            u32 bits = mask;
            while (bits != 0) {
                const int bit = std::countr_zero(bits);
                const u32 source_mask = mask ^ (u32{1} << bit);
                const u32 source = space.offsets[weight - 1] + space.rank(source_mask, weight - 1);
                edges[bit].push_back((u64{destination} << 32) | source);
                bits &= bits - 1;
            }
            ++rank;
            const u32 low = mask & (~mask + 1);
            const u32 high = mask + low;
            mask = high | (((high ^ mask) >> 2) / low);
        }
        require(rank == space.choose[30][weight], "edge layer enumeration");
    }
    for (const auto& list : edges) require(list.size() == per_bit, "edge symmetry");
    return edges;
}

struct ScaledEvent {
    i64 position;
    int vertex;
    bool enter_safe;
};

struct Coefficients {
    i64 denominator;
    std::vector<i64> mass;
    std::vector<i64> components;
    std::size_t walls;
};

Coefficients build_q_coefficients(
    int q,
    const BaseGeometry& base,
    const std::vector<u32>& potential_keys,
    const std::unordered_map<u32, u16>& key_index
) {
    const i64 denominator = std::lcm(base.denominator, i64{14} * q);
    const i64 scale = denominator / base.denominator;
    std::vector<ScaledEvent> events;
    events.reserve(base.events.size() + 2 * q);
    std::vector<i64> walls = {0, denominator};
    walls.reserve(base.events.size() + 2 * q + 2);
    for (const RawEvent& event : base.events) {
        events.push_back({event.position * scale, event.vertex, event.enter_safe});
        walls.push_back(event.position * scale);
    }
    const i64 unit = denominator / (14LL * q);
    for (int tooth = 0; tooth < q; ++tooth) {
        const i64 enter = (14LL * tooth + 1) * unit;
        const i64 leave = (14LL * tooth + 13) * unit;
        events.push_back({enter, 30, true});
        events.push_back({leave, 30, false});
        walls.push_back(enter);
        walls.push_back(leave);
    }
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
    bool previous_q_safe = false;
    u32 previous_failure = failure;
    i64 raw_q_mass = 0;
    i64 raw_q_components = 0;
    std::size_t event_index = 0;
    for (std::size_t wall_index = 0; wall_index + 1 < walls.size(); ++wall_index) {
        const i64 left = walls[wall_index];
        while (event_index < events.size() && events[event_index].position == left) {
            const ScaledEvent& event = events[event_index++];
            if (event.vertex == 30) {
                q_safe = event.enter_safe;
            } else if (event.enter_safe) {
                failure &= ~(u32{1} << event.vertex);
            } else {
                failure |= u32{1} << event.vertex;
            }
        }
        if (q_safe) {
            const i64 width = walls[wall_index + 1] - left;
            raw_q_mass += width;
            ++raw_q_components;
            if (std::popcount(failure) <= 30 - BODY_WEIGHT) {
                auto it = key_index.find(failure);
                require(it != key_index.end(), "mass key absent");
                mass[it->second] += width;
                components[it->second] += 1;
            }
        }
        if (q_safe && previous_q_safe) {
            --raw_q_components;
            const u32 joined = failure | previous_failure;
            if (std::popcount(joined) <= 30 - BODY_WEIGHT) {
                auto it = key_index.find(joined);
                require(it != key_index.end(), "component key absent");
                components[it->second] -= 1;
            }
        }
        previous_q_safe = q_safe;
        previous_failure = failure;
    }
    require(event_index == events.size(), "q event exhaustion");
    require(7 * raw_q_mass == 6 * denominator, "empty-body q mass");
    require(raw_q_components == q, "empty-body q components");
    return {denominator, std::move(mass), std::move(components), walls.size()};
}

struct RaySummary {
    u64 cases = 0;
    u64 nonpositive_mass_or_components = 0;
    u64 nonpositive_surplus = 0;
    u64 violations_931 = 0;
    u64 exact_931 = 0;
    int exact_931_q = 0;
    u32 exact_931_body = 0;
    i64 exact_931_denominator = 0;
    i64 exact_931_mass = 0;
    i64 exact_931_components = 0;
    i64 exact_931_surplus = 0;
    i64 minimum_surplus = std::numeric_limits<i64>::max();
    u32 minimum_surplus_body = 0;
    i64 minimum_surplus_mass = 0;
    i64 minimum_surplus_components = 0;
    u64 maximum_threshold = 0;
    u32 maximum_threshold_witness = 0;
    i64 maximum_threshold_mass = 0;
    i64 maximum_threshold_components = 0;
    i64 maximum_threshold_surplus = 0;
    u64 digest_xor = 0;
    u64 digest_sum = 0;
};

struct GlobalSummary {
    u64 cases = 0;
    u64 nonpositive_mass_or_components = 0;
    u64 nonpositive_surplus = 0;
    u64 violations_931 = 0;
    u64 exact_931 = 0;
    int exact_931_q = 0;
    u32 exact_931_body = 0;
    i64 exact_931_denominator = 0;
    i64 exact_931_mass = 0;
    i64 exact_931_components = 0;
    i64 exact_931_surplus = 0;
    i64 normalized_surplus = 0;
    i64 normalized_denominator = 1;
    int normalized_q = 0;
    u32 normalized_body = 0;
    i64 normalized_mass = 0;
    i64 normalized_components = 0;
    u64 digest_xor = 0;
    u64 digest_sum = 0;
};

}  // namespace

int main() {
    const auto started = std::chrono::steady_clock::now();
    const BaseGeometry base = build_base_geometry();
    const MaskSpace space;
    const std::vector<u32> potential_keys = build_potential_keys(base);
    std::unordered_map<u32, u16> key_index;
    key_index.reserve(2 * potential_keys.size());
    for (u32 i = 0; i < potential_keys.size(); ++i) {
        key_index.emplace(potential_keys[i], static_cast<u16>(i));
    }
    const Incidence incidence = build_incidence(potential_keys, space);
    const auto edges = build_transform_edges(space);

    u64 edge_count = 0;
    for (const auto& list : edges) edge_count += list.size();
    std::vector<i64> mass(space.total_masks);
    std::vector<i64> components(space.total_masks);
    std::array<u64, Q_MAXIMUM + 1> ray_threshold{};
    std::array<u32, Q_MAXIMUM + 1> ray_witness{};
    GlobalSummary global;
    int q_count = 0;
    std::size_t wall_minimum = std::numeric_limits<std::size_t>::max();
    std::size_t wall_maximum = 0;

    for (int q = 1; q <= Q_MAXIMUM; ++q) {
        if (is_pool_label(q)) continue;
        ++q_count;
        const Coefficients coefficients = build_q_coefficients(q, base, potential_keys, key_index);
        wall_minimum = std::min(wall_minimum, coefficients.walls);
        wall_maximum = std::max(wall_maximum, coefficients.walls);

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

        RaySummary ray;
        const u32 body_offset = space.offsets[BODY_WEIGHT];
#pragma omp parallel
        {
            RaySummary local;
#pragma omp for schedule(static)
            for (i64 raw_body = 0; raw_body < static_cast<i64>(space.body_masks.size()); ++raw_body) {
                const u32 body_index = static_cast<u32>(raw_body);
                const u32 body = space.body_masks[body_index];
                const i64 m = mass[body_offset + body_index];
                const i64 c = components[body_offset + body_index];
                ++local.cases;
                if (m <= 0 || c <= 0) ++local.nonpositive_mass_or_components;
                const i64 surplus = 54 * m - 4 * coefficients.denominator;
                if (surplus <= 0) ++local.nonpositive_surplus;
                if (
                    surplus < local.minimum_surplus ||
                    (surplus == local.minimum_surplus && body < local.minimum_surplus_body)
                ) {
                    local.minimum_surplus = surplus;
                    local.minimum_surplus_body = body;
                    local.minimum_surplus_mass = m;
                    local.minimum_surplus_components = c;
                }
                if (surplus > 0 && c > 0) {
                    const u128 numerator = u128{54} * static_cast<u64>(c) *
                                           static_cast<u64>(coefficients.denominator);
                    const u128 denominator = u128{7} * static_cast<u64>(surplus);
                    if (numerator > u128{local.maximum_threshold} * denominator) {
                        const u64 threshold = ceil_div(numerator, denominator);
                        if (
                            threshold > local.maximum_threshold ||
                            (threshold == local.maximum_threshold &&
                             (local.maximum_threshold_witness == 0 ||
                              body < local.maximum_threshold_witness))
                        ) {
                            local.maximum_threshold = threshold;
                            local.maximum_threshold_witness = body;
                            local.maximum_threshold_mass = m;
                            local.maximum_threshold_components = c;
                            local.maximum_threshold_surplus = surplus;
                        }
                    }
                    const u128 bound_931 = u128{6517} * static_cast<u64>(surplus);
                    const u128 bound_930 = u128{6510} * static_cast<u64>(surplus);
                    if (numerator > bound_931) ++local.violations_931;
                    if (numerator > bound_930 && numerator <= bound_931) {
                        ++local.exact_931;
                        local.exact_931_q = q;
                        if (local.exact_931_body == 0 || body < local.exact_931_body) {
                            local.exact_931_body = body;
                            local.exact_931_denominator = coefficients.denominator;
                            local.exact_931_mass = m;
                            local.exact_931_components = c;
                            local.exact_931_surplus = surplus;
                        }
                    }
                }
                const u64 word = cheap_mix(
                    (u64{static_cast<u32>(q)} << 32) ^ body ^
                    std::rotl(static_cast<u64>(m), 11) ^
                    std::rotl(static_cast<u64>(c), 37) ^
                    std::rotl(static_cast<u64>(surplus), 23)
                );
                local.digest_xor ^= word;
                local.digest_sum += word;
            }
#pragma omp critical
            {
                ray.cases += local.cases;
                ray.nonpositive_mass_or_components += local.nonpositive_mass_or_components;
                ray.nonpositive_surplus += local.nonpositive_surplus;
                ray.violations_931 += local.violations_931;
                if (local.exact_931 != 0) {
                    ray.exact_931 += local.exact_931;
                    ray.exact_931_q = local.exact_931_q;
                    if (ray.exact_931_body == 0 || local.exact_931_body < ray.exact_931_body) {
                        ray.exact_931_body = local.exact_931_body;
                        ray.exact_931_denominator = local.exact_931_denominator;
                        ray.exact_931_mass = local.exact_931_mass;
                        ray.exact_931_components = local.exact_931_components;
                        ray.exact_931_surplus = local.exact_931_surplus;
                    }
                }
                if (
                    local.minimum_surplus < ray.minimum_surplus ||
                    (local.minimum_surplus == ray.minimum_surplus &&
                     local.minimum_surplus_body < ray.minimum_surplus_body)
                ) {
                    ray.minimum_surplus = local.minimum_surplus;
                    ray.minimum_surplus_body = local.minimum_surplus_body;
                    ray.minimum_surplus_mass = local.minimum_surplus_mass;
                    ray.minimum_surplus_components = local.minimum_surplus_components;
                }
                if (
                    local.maximum_threshold > ray.maximum_threshold ||
                    (local.maximum_threshold == ray.maximum_threshold &&
                     local.maximum_threshold_witness != 0 &&
                     (ray.maximum_threshold_witness == 0 ||
                      local.maximum_threshold_witness < ray.maximum_threshold_witness))
                ) {
                    ray.maximum_threshold = local.maximum_threshold;
                    ray.maximum_threshold_witness = local.maximum_threshold_witness;
                    ray.maximum_threshold_mass = local.maximum_threshold_mass;
                    ray.maximum_threshold_components = local.maximum_threshold_components;
                    ray.maximum_threshold_surplus = local.maximum_threshold_surplus;
                }
                ray.digest_xor ^= local.digest_xor;
                ray.digest_sum += local.digest_sum;
            }
        }
        require(ray.cases == BODY_COUNT, "ray body count");
        require(ray.maximum_threshold > 0, "ray maximum threshold");
        ray_threshold[q] = ray.maximum_threshold;
        ray_witness[q] = ray.maximum_threshold_witness;
        global.cases += ray.cases;
        global.nonpositive_mass_or_components += ray.nonpositive_mass_or_components;
        global.nonpositive_surplus += ray.nonpositive_surplus;
        global.violations_931 += ray.violations_931;
        if (ray.exact_931 != 0) {
            global.exact_931 += ray.exact_931;
            global.exact_931_q = ray.exact_931_q;
            global.exact_931_body = ray.exact_931_body;
            global.exact_931_denominator = ray.exact_931_denominator;
            global.exact_931_mass = ray.exact_931_mass;
            global.exact_931_components = ray.exact_931_components;
            global.exact_931_surplus = ray.exact_931_surplus;
        }
        if (
            global.normalized_q == 0 ||
            i128{ray.minimum_surplus} * global.normalized_denominator <
                i128{global.normalized_surplus} * coefficients.denominator ||
            (
                i128{ray.minimum_surplus} * global.normalized_denominator ==
                    i128{global.normalized_surplus} * coefficients.denominator &&
                (q < global.normalized_q ||
                 (q == global.normalized_q &&
                  ray.minimum_surplus_body < global.normalized_body))
            )
        ) {
            global.normalized_surplus = ray.minimum_surplus;
            global.normalized_denominator = coefficients.denominator;
            global.normalized_q = q;
            global.normalized_body = ray.minimum_surplus_body;
            global.normalized_mass = ray.minimum_surplus_mass;
            global.normalized_components = ray.minimum_surplus_components;
        }
        global.digest_xor ^= ray.digest_xor;
        global.digest_sum += ray.digest_sum;

        if (q == 1 || q == 543 || q == 930 || q == 1305 || q == 1306) {
            std::cerr << "PROGRESS Q " << q
                      << " MIN_SURPLUS " << ray.minimum_surplus
                      << " EXACT931 " << ray.exact_931
                      << " ELAPSED "
                      << std::chrono::duration_cast<std::chrono::seconds>(
                             std::chrono::steady_clock::now() - started).count()
                      << "s\n";
        }
    }

    require(q_count == 1'276, "q universe");
    require(global.cases == EXPECTED_CASES, "case universe");
    require(global.nonpositive_mass_or_components == 0, "nonpositive mass/components");
    require(global.nonpositive_surplus == 0, "nonpositive limit surplus");
    require(global.violations_931 == 0, "threshold above 931");

    Fnv1a64 ray_table_fingerprint;
    std::map<u64, u64> ray_threshold_histogram;
    u64 ray_maximum = 0;
    int ray_maximum_q = 0;
    u32 ray_maximum_body = 0;
    for (int q = 1; q <= Q_MAXIMUM; ++q) {
        if (is_pool_label(q)) continue;
        require(ray_threshold[q] != 0, "missing ray threshold");
        ray_table_fingerprint.add(q);
        ray_table_fingerprint.add(ray_threshold[q]);
        ray_table_fingerprint.add(ray_witness[q]);
        ++ray_threshold_histogram[ray_threshold[q]];
        if (ray_threshold[q] > ray_maximum) {
            ray_maximum = ray_threshold[q];
            ray_maximum_q = q;
            ray_maximum_body = ray_witness[q];
        }
    }
    require(ray_maximum == 931, "ray-table maximum");

    int last_diagonal_bad = 0;
    u64 last_diagonal_bad_threshold = 0;
    u64 diagonal_bad_count = 0;
    Fnv1a64 diagonal_bad_fingerprint;
    for (int q = 1; q <= 1'289; ++q) {
        if (is_pool_label(q)) continue;
        if (ray_threshold[q] > static_cast<u64>(q)) {
            ++diagonal_bad_count;
            last_diagonal_bad = q;
            last_diagonal_bad_threshold = ray_threshold[q];
            diagonal_bad_fingerprint.add(q);
            diagonal_bad_fingerprint.add(ray_threshold[q]);
        }
    }

    int conservative_cutoff = 0;
    int cutoff_below_argmax = 0;
    u64 cutoff_below_maximum = 0;
    for (int candidate = 1; candidate <= Q_MAXIMUM + 1; ++candidate) {
        u64 below_maximum = 0;
        int below_argmax = 0;
        for (int q = 1; q < candidate && q <= Q_MAXIMUM; ++q) {
            if (is_pool_label(q)) continue;
            if (ray_threshold[q] > below_maximum) {
                below_maximum = ray_threshold[q];
                below_argmax = q;
            }
        }
        bool diagonal_ok = true;
        for (int q = candidate; q <= 1'289; ++q) {
            if (is_pool_label(q)) continue;
            if (ray_threshold[q] > static_cast<u64>(q)) {
                diagonal_ok = false;
                break;
            }
        }
        if (below_maximum <= static_cast<u64>(candidate) && diagonal_ok) {
            conservative_cutoff = candidate;
            cutoff_below_argmax = below_argmax;
            cutoff_below_maximum = below_maximum;
            break;
        }
    }
    require(conservative_cutoff != 0, "conservative cutoff absent");

    const auto first_outsider_at_or_above = [](int value) {
        while (is_pool_label(value)) ++value;
        return value;
    };
    const auto next_outsider = [&](int value) {
        return first_outsider_at_or_above(value + 1);
    };
    int smaller_ray_cutoff = 0;
    for (int candidate = 1; candidate <= Q_MAXIMUM + 1; ++candidate) {
        const int first = first_outsider_at_or_above(candidate);
        bool good = true;
        for (int q = 1; q < first; ++q) {
            if (!is_pool_label(q) && ray_threshold[q] > static_cast<u64>(first)) {
                good = false;
                break;
            }
        }
        if (!good) continue;
        for (int q = first; q <= Q_MAXIMUM; ++q) {
            if (!is_pool_label(q) &&
                ray_threshold[q] > static_cast<u64>(next_outsider(q))) {
                good = false;
                break;
            }
        }
        if (good) {
            smaller_ray_cutoff = candidate;
            break;
        }
    }
    require(smaller_ray_cutoff != 0, "smaller-ray cutoff absent");

    u64 residual_uncertified_pairs = 0;
    Fnv1a64 residual_pair_fingerprint;
    int residual_last_a = 0;
    int residual_last_b = 0;
    for (int a = 1; a < conservative_cutoff; ++a) {
        if (is_pool_label(a)) continue;
        for (int b = a + 1; b < conservative_cutoff; ++b) {
            if (is_pool_label(b)) continue;
            if (static_cast<u64>(b) >= ray_threshold[a]) continue;
            ++residual_uncertified_pairs;
            residual_pair_fingerprint.add(a);
            residual_pair_fingerprint.add(b);
            residual_last_a = a;
            residual_last_b = b;
        }
    }

    u64 two_orientation_residual_pairs = 0;
    Fnv1a64 two_orientation_fingerprint;
    int two_orientation_lex_last_a = 0;
    int two_orientation_lex_last_b = 0;
    int two_orientation_max_endpoint_a = 0;
    int two_orientation_max_endpoint_b = 0;
    for (int a = 1; a <= Q_MAXIMUM; ++a) {
        if (is_pool_label(a)) continue;
        for (int b = a + 1; b <= Q_MAXIMUM; ++b) {
            if (is_pool_label(b)) continue;
            if (static_cast<u64>(b) >= ray_threshold[a]) continue;
            if (static_cast<u64>(a) >= ray_threshold[b]) continue;
            ++two_orientation_residual_pairs;
            two_orientation_fingerprint.add(a);
            two_orientation_fingerprint.add(b);
            two_orientation_lex_last_a = a;
            two_orientation_lex_last_b = b;
            if (b > two_orientation_max_endpoint_b) {
                two_orientation_max_endpoint_a = a;
                two_orientation_max_endpoint_b = b;
            }
        }
    }
    const int strongest_cutoff = two_orientation_max_endpoint_b + 1;

    std::cout << "LRC14_ALL_FIXED_Q_RAYS_INDEPENDENT_REFEREE\n";
    std::cout << "UNIVERSES Q " << q_count
              << " BODIES_PER_Q " << BODY_COUNT
              << " CASES " << global.cases
              << " TRUNCATED_MASKS " << space.total_masks
              << " POTENTIAL_KEYS " << potential_keys.size()
              << " INCIDENCES " << incidence.keys.size()
              << " TRANSFORM_EDGES " << edge_count << '\n';
    std::cout << "GEOMETRY BASE_D " << base.denominator
              << " BASE_WALLS " << base.cells.size() + 1
              << " Q_WALL_RANGE " << wall_minimum << ".." << wall_maximum << '\n';
    std::cout << "POSITIVITY MASS_OR_COMPONENT_NONPOSITIVE "
              << global.nonpositive_mass_or_components
              << " LIMIT_SURPLUS_NONPOSITIVE " << global.nonpositive_surplus << '\n';
    std::cout << "GLOBAL_CERTIFICATE MAX_K 931 VIOLATIONS " << global.violations_931
              << " EXACT_COUNT " << global.exact_931
              << " WITNESS_Q " << global.exact_931_q
              << " WITNESS_B " << labels(global.exact_931_body) << '\n';
    const u128 maximum_numerator =
        u128{54} * static_cast<u64>(global.exact_931_components) *
        static_cast<u64>(global.exact_931_denominator);
    const u128 maximum_bound_930 = u128{6510} * static_cast<u64>(global.exact_931_surplus);
    const u128 maximum_bound_931 = u128{6517} * static_cast<u64>(global.exact_931_surplus);
    std::cout << "MAX_WITNESS_GEOMETRY D " << global.exact_931_denominator
              << " MASS_TICKS " << global.exact_931_mass
              << " COMPONENTS " << global.exact_931_components
              << " SURPLUS54_TICKS " << global.exact_931_surplus
              << " NUMERATOR54CD " << decimal(maximum_numerator)
              << " BOUND930 " << decimal(maximum_bound_930)
              << " BOUND931 " << decimal(maximum_bound_931) << '\n';
    std::cout << "NORMALIZED_SURPLUS_MIN " << global.normalized_surplus
              << '/' << global.normalized_denominator
              << " Q " << global.normalized_q
              << " B " << labels(global.normalized_body)
              << " MASS_TICKS " << global.normalized_mass
              << " COMPONENTS " << global.normalized_components << '\n';
    std::cout << "COMMUTATIVE_LEDGER_DIGEST XOR " << std::hex << std::setw(16)
              << std::setfill('0') << global.digest_xor
              << " SUM " << std::setw(16) << global.digest_sum
              << std::dec << std::setfill(' ') << '\n';
    std::cout << "RAY_TABLE MAX_K " << ray_maximum
              << " Q " << ray_maximum_q
              << " B " << labels(ray_maximum_body)
              << " FNV_Q_K_BODY " << std::hex << std::setw(16)
              << std::setfill('0') << ray_table_fingerprint.value
              << std::dec << std::setfill(' ')
              << " DISTINCT_K " << ray_threshold_histogram.size() << '\n';
    std::cout << "CONSERVATIVE_Q_LEQ_SELF_CUTOFF N " << conservative_cutoff
              << " BELOW_MAX_K " << cutoff_below_maximum
              << " BELOW_ARGMAX_Q " << cutoff_below_argmax
              << " DIAGONAL_BAD_COUNT " << diagonal_bad_count
              << " LAST_DIAGONAL_BAD_Q " << last_diagonal_bad
              << " LAST_DIAGONAL_BAD_K " << last_diagonal_bad_threshold
              << " BAD_FNV " << std::hex << std::setw(16) << std::setfill('0')
              << diagonal_bad_fingerprint.value << std::dec << std::setfill(' ') << '\n';
    std::cout << "SMALLER_RAY_EXACT_CUTOFF N " << smaller_ray_cutoff
              << " FIRST_OUTSIDER " << first_outsider_at_or_above(smaller_ray_cutoff)
              << '\n';
    std::cout << "ONE_ORIENTATION_RESIDUAL_PAIRS COUNT " << residual_uncertified_pairs
              << " LAST " << residual_last_a << ',' << residual_last_b
              << " FNV " << std::hex << std::setw(16) << std::setfill('0')
              << residual_pair_fingerprint.value << std::dec << std::setfill(' ') << '\n';
    std::cout << "TWO_ORIENTATION_RESIDUAL_PAIRS COUNT "
              << two_orientation_residual_pairs
              << " STRONGEST_CUTOFF " << strongest_cutoff
              << " MAX_ENDPOINT_EDGE " << two_orientation_max_endpoint_a << ','
              << two_orientation_max_endpoint_b
              << " LEX_LAST " << two_orientation_lex_last_a << ','
              << two_orientation_lex_last_b
              << " FNV " << std::hex << std::setw(16) << std::setfill('0')
              << two_orientation_fingerprint.value << std::dec << std::setfill(' ') << '\n';
    std::cout << "VERDICT EXACT_ALL_Q_COFINAL_TAIL_MAX_931 LRC14_OPEN\n";
    return 0;
}
