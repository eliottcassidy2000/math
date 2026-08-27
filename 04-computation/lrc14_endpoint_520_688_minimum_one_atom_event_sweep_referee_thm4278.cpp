// Independent finite-exact event-sweep referee for THM-4278 at the former
// LRC14 fixed-pool top edge.
//
// This implementation deliberately does not include or call the endpoint-zeta,
// cocycle-atom, superset-transform, or primary minimum-augmentation code.  It
// reconstructs the inherited carrier only from frozen canonical transcripts,
// then obtains every mass from a literal wall-event sweep on one exact common
// denominator.  The sweep changes safety states at the displayed wall events;
// it does not sample cell midpoints.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};

constexpr std::array<std::pair<int, int>, 59> BASE_PAIRS = {{
    {616,755}, {616,756}, {616,757}, {616,758}, {616,759},
    {616,760}, {616,761}, {616,762}, {616,763}, {616,764},
    {616,765}, {616,766}, {616,767}, {616,768}, {698,755},
    {698,757}, {704,755}, {704,757}, {704,758}, {704,759},
    {704,761}, {704,762}, {704,763}, {704,764}, {704,765},
    {721,755}, {721,757}, {721,758}, {721,759}, {721,761},
    {721,762}, {721,763}, {721,764}, {721,765}, {721,766},
    {721,767}, {721,768}, {726,755}, {726,757}, {726,758},
    {726,761}, {732,755}, {732,757}, {732,761}, {732,762},
    {732,763}, {744,762}, {744,763}, {744,765}, {744,766},
    {744,768}, {750,762}, {750,763}, {750,765}, {750,766},
    {750,768}, {765,766}, {765,768}, {766,768}
}};

constexpr u64 BODY_COUNT = UINT64_C(14307150);
constexpr u32 EXPECTED_FIRST_FAILURE = UINT32_C(0x07187008);
constexpr u32 EXPECTED_SECOND_FAILURE = UINT32_C(0x27503008);
constexpr u32 EXPECTED_REPAIR = UINT32_C(0x00048ec1);

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    if (negative) value = -value;
    std::string answer;
    while (value != 0) {
        answer.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

i128 gcd128(i128 left, i128 right) {
    if (left < 0) left = -left;
    if (right < 0) right = -right;
    while (right != 0) {
        const i128 next = left % right;
        left = right;
        right = next;
    }
    return left;
}

std::string fraction(i128 numerator, i128 denominator) {
    require(denominator > 0, "nonpositive denominator");
    const i128 divisor = gcd128(numerator, denominator);
    return decimal(numerator / divisor) + "/" + decimal(denominator / divisor);
}

i64 exact_lcm(i64 left, i64 right) {
    require(left > 0 && right > 0, "nonpositive lcm input");
    const i128 wide = static_cast<i128>(left / std::gcd(left, right)) * right;
    require(wide <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(wide);
}

struct Fnv {
    u64 value = UINT64_C(0xcbf29ce484222325);

    void add(u64 word) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            value ^= (word >> (8 * byte)) & UINT64_C(0xff);
            value *= UINT64_C(0x100000001b3);
        }
    }
};

u64 parse_hex(const std::string& text) {
    require(!text.empty(), "empty hexadecimal word");
    std::size_t used = 0;
    const u64 value = std::stoull(text, &used, 16);
    require(used == text.size(), "trailing hexadecimal text");
    return value;
}

std::string labels(u32 mask) {
    std::string answer;
    for (unsigned bit = 0; bit < POOL.size(); ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        if (!answer.empty()) answer.push_back(',');
        answer += std::to_string(POOL[bit]);
    }
    return answer;
}

struct Transcript {
    int q = -1;
    int r = -1;
    u64 declared_count = 0;
    u64 declared_fnv = 0;
    std::vector<u32> masks;
};

Transcript read_transcript(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open " + path.string());
    Transcript result;
    bool saw_pair = false;
    bool saw_cert = false;
    bool saw_masks = false;
    bool saw_closed = false;
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("PAIR ", 0) == 0) {
            std::istringstream row(line.substr(5));
            char comma = 0;
            row >> result.q >> comma >> result.r;
            require(static_cast<bool>(row) && comma == ',',
                    "bad PAIR row in " + path.string());
            saw_pair = true;
        } else if (line.rfind("PREFIX_CERT ", 0) == 0) {
            std::istringstream row(line);
            std::string prefix, size_word, fnv_word, hex_word;
            row >> prefix >> size_word >> result.declared_count
                >> fnv_word >> hex_word;
            require(static_cast<bool>(row) && prefix == "PREFIX_CERT" &&
                        size_word == "SIZE" && fnv_word == "FNV",
                    "bad PREFIX_CERT row in " + path.string());
            result.declared_fnv = parse_hex(hex_word);
            saw_cert = true;
        } else if (line.rfind("PREFIX_MASKS_HEX ", 0) == 0) {
            const std::string words = line.substr(std::string("PREFIX_MASKS_HEX ").size());
            std::size_t begin = 0;
            while (begin < words.size()) {
                const std::size_t comma = words.find(',', begin);
                const std::string token = words.substr(
                    begin, comma == std::string::npos ? std::string::npos
                                                      : comma - begin);
                const u64 wide = parse_hex(token);
                require(wide < (UINT64_C(1) << 30),
                        "prefix mask outside pool");
                result.masks.push_back(static_cast<u32>(wide));
                if (comma == std::string::npos) break;
                begin = comma + 1;
            }
            saw_masks = true;
        } else if (line.rfind("VERDICT EVERY_BODY_CLOSED", 0) == 0) {
            saw_closed = true;
        }
    }
    require(saw_pair && saw_cert && saw_masks && saw_closed,
            "incomplete transcript " + path.string());
    require(result.masks.size() == result.declared_count,
            "declared prefix count mismatch");
    Fnv ledger;
    std::set<u32> unique;
    for (u32 mask : result.masks) {
        require(std::popcount(mask) == 8, "non-rank-eight prefix mask");
        require(unique.insert(mask).second, "duplicate mask inside prefix");
        ledger.add(mask);
    }
    require(ledger.value == result.declared_fnv,
            "declared prefix FNV mismatch");
    return result;
}

struct Carrier {
    std::vector<u32> masks;
    std::array<u64, 5> stage_count{};
    std::array<u64, 5> stage_fnv{};
};

Carrier reconstruct_carrier(const std::filesystem::path& replay,
                            const std::filesystem::path& round_results) {
    Carrier carrier;
    std::set<u32> seen;
    auto append = [&](const Transcript& transcript) {
        for (u32 mask : transcript.masks) {
            if (seen.insert(mask).second) carrier.masks.push_back(mask);
        }
    };
    for (const auto& [q, r] : BASE_PAIRS) {
        const Transcript transcript = read_transcript(
            replay / (std::to_string(q) + "_" + std::to_string(r) + ".out"));
        require(transcript.q == q && transcript.r == r,
                "base transcript pair mismatch");
        append(transcript);
    }
    Fnv ledger;
    for (u32 mask : carrier.masks) ledger.add(mask);
    carrier.stage_count[0] = carrier.masks.size();
    carrier.stage_fnv[0] = ledger.value;
    require(carrier.masks.size() == 4675 &&
                ledger.value == UINT64_C(0xce4e76ec11df057c),
            "base carrier reconstruction mismatch");

    struct Source {
        const char* name;
        int q;
        int r;
        u64 prefix_count;
        u64 prefix_fnv;
        u64 union_count;
        u64 union_fnv;
    };
    constexpr std::array<Source, 4> sources = {{
        {"full_discovery_416_704_O3.semantic.out", 416, 704, 2608,
         UINT64_C(0x18ff663a123e684e), 4733,
         UINT64_C(0xa7b046289655c733)},
        {"full_discovery_416_700_O3.semantic.out", 416, 700, 5894,
         UINT64_C(0x701c0233c8f8abeb), 7703, 0},
        {"full_discovery_520_700_O3.semantic.out", 520, 700, 4557,
         UINT64_C(0xd8466558bcd8ef9d), 7986,
         UINT64_C(0xbaef1d2f49444638)},
        {"full_discovery_384_694_O3.semantic.out", 384, 694, 5307,
         UINT64_C(0x7b9a46c60e8514f8), 8319,
         UINT64_C(0xe08b227730f6793c)},
    }};
    for (std::size_t index = 0; index < sources.size(); ++index) {
        const Source& source = sources[index];
        const Transcript transcript = read_transcript(round_results / source.name);
        require(transcript.q == source.q && transcript.r == source.r &&
                    transcript.masks.size() == source.prefix_count &&
                    transcript.declared_fnv == source.prefix_fnv,
                "round source transcript mismatch");
        append(transcript);
        ledger = Fnv{};
        for (u32 mask : carrier.masks) ledger.add(mask);
        carrier.stage_count[index + 1] = carrier.masks.size();
        carrier.stage_fnv[index + 1] = ledger.value;
        require(carrier.masks.size() == source.union_count,
                "round carrier union count mismatch");
        if (source.union_fnv != 0) {
            require(ledger.value == source.union_fnv,
                    "round carrier union FNV mismatch");
        }
    }
    require(carrier.masks.size() == 8319 &&
                carrier.stage_fnv[4] == UINT64_C(0xe08b227730f6793c),
            "final carrier reconstruction mismatch");
    return carrier;
}

struct Event {
    i64 position = 0;
    unsigned speed_index = 0;
    bool entering = false;
};

struct LiteralWall {
    i64 grid = 0;
    i64 pair_ticks = 0;
    u64 cells = 0;
    u64 pair_safe_cells = 0;
    std::vector<std::pair<u32, i64>> weighted_patterns;
};

LiteralWall build_event_sweep(int q, int r) {
    std::array<int, 32> speeds{};
    std::copy(POOL.begin(), POOL.end(), speeds.begin());
    speeds[30] = q;
    speeds[31] = r;

    i64 grid = 1;
    for (int speed : speeds) grid = exact_lcm(grid, 14LL * speed);

    std::vector<Event> events;
    for (unsigned index = 0; index < speeds.size(); ++index) {
        const i64 speed = speeds[index];
        const i64 unit = grid / (14 * speed);
        require(unit * 14 * speed == grid, "nonintegral literal wall unit");
        for (i64 tooth = 0; tooth < speed; ++tooth) {
            events.push_back({(14 * tooth + 1) * unit, index, true});
            events.push_back({(14 * tooth + 13) * unit, index, false});
        }
    }
    std::sort(events.begin(), events.end(), [](const Event& left,
                                               const Event& right) {
        if (left.position != right.position) return left.position < right.position;
        if (left.speed_index != right.speed_index) {
            return left.speed_index < right.speed_index;
        }
        return left.entering < right.entering;
    });

    std::array<bool, 32> safe{};
    std::map<u32, i64> weights;
    LiteralWall result;
    result.grid = grid;
    i64 previous = 0;
    std::size_t cursor = 0;
    auto consume_cell = [&](i64 right) {
        require(right > previous, "nonpositive literal cell");
        ++result.cells;
        const i64 width = right - previous;
        if (safe[30] && safe[31]) {
            ++result.pair_safe_cells;
            u32 failure = 0;
            for (unsigned bit = 0; bit < POOL.size(); ++bit) {
                if (!safe[bit]) failure |= u32{1} << bit;
            }
            weights[failure] += width;
            result.pair_ticks += width;
        }
        previous = right;
    };
    while (cursor < events.size()) {
        const i64 position = events[cursor].position;
        require(position > 0 && position < grid, "wall outside open circle");
        if (position > previous) consume_cell(position);
        while (cursor < events.size() && events[cursor].position == position) {
            const Event event = events[cursor++];
            require(safe[event.speed_index] != event.entering,
                    "literal wall state did not toggle");
            safe[event.speed_index] = event.entering;
        }
    }
    if (previous < grid) consume_cell(grid);
    require(std::none_of(safe.begin(), safe.end(), [](bool value) { return value; }),
            "literal circle failed to return to unsafe origin state");
    result.weighted_patterns.assign(weights.begin(), weights.end());
    i128 total = 0;
    for (const auto& [mask, width] : result.weighted_patterns) {
        (void)mask;
        total += width;
    }
    require(total == result.pair_ticks, "pattern compression changed mass");
    return result;
}

i64 mass_ticks(const LiteralWall& wall, u32 repair) {
    i64 total = 0;
    for (const auto& [failed, width] : wall.weighted_patterns) {
        if ((failed & ~repair) == 0) total += width;
    }
    return total;
}

std::vector<u32> active_subdeck(const std::vector<u32>& carrier,
                                const LiteralWall& wall,
                                u64& equalities) {
    std::vector<u32> active;
    equalities = 0;
    for (u32 repair : carrier) {
        const i128 margin = static_cast<i128>(63) * mass_ticks(wall, repair) -
                            static_cast<i128>(4) * wall.grid;
        if (margin == 0) ++equalities;
        if (margin >= 0) active.push_back(repair);
    }
    return active;
}

void generate(const std::vector<int>& bits, int begin, int need,
              u32 mask, std::vector<u32>& output) {
    if (need == 0) {
        output.push_back(mask);
        return;
    }
    const int last = static_cast<int>(bits.size()) - need;
    for (int cursor = begin; cursor <= last; ++cursor) {
        generate(bits, cursor + 1, need - 1,
                 mask | (u32{1} << bits[cursor]), output);
    }
}

struct BodyScan {
    u64 bodies = 0;
    u64 checks = 0;
    u64 max_checks = 0;
    std::vector<u32> failures;
};

BodyScan scan_bodies(const std::vector<u32>& bodies,
                     const std::vector<u32>& active) {
    const unsigned detected = std::thread::hardware_concurrency();
    const unsigned threads_count = std::max(1u, std::min(8u, detected ? detected : 1u));
    std::vector<BodyScan> local(threads_count);
    std::vector<std::thread> threads;
    for (unsigned thread = 0; thread < threads_count; ++thread) {
        threads.emplace_back([&, thread]() {
            const std::size_t begin = bodies.size() * thread / threads_count;
            const std::size_t end = bodies.size() * (thread + 1) / threads_count;
            BodyScan report;
            for (std::size_t index = begin; index < end; ++index) {
                const u32 body = bodies[index];
                u64 used = 0;
                bool covered = false;
                for (u32 repair : active) {
                    ++used;
                    if ((body & repair) == 0) {
                        covered = true;
                        break;
                    }
                }
                ++report.bodies;
                report.checks += used;
                report.max_checks = std::max(report.max_checks, used);
                if (!covered) report.failures.push_back(body);
            }
            local[thread] = std::move(report);
        });
    }
    for (std::thread& thread : threads) thread.join();
    BodyScan answer;
    for (BodyScan& report : local) {
        answer.bodies += report.bodies;
        answer.checks += report.checks;
        answer.max_checks = std::max(answer.max_checks, report.max_checks);
        answer.failures.insert(answer.failures.end(),
                               report.failures.begin(), report.failures.end());
    }
    std::sort(answer.failures.begin(), answer.failures.end());
    return answer;
}

u64 choose(unsigned n, unsigned k) {
    if (k > n) return 0;
    k = std::min(k, n - k);
    u64 result = 1;
    for (unsigned i = 1; i <= k; ++i) result = result * (n - k + i) / i;
    return result;
}

u64 colex_rank0(u32 mask) {
    u64 rank = 0;
    unsigned order = 1;
    for (unsigned bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        rank += choose(bit, order);
        ++order;
    }
    require(order == 9, "colex rank applied outside rank eight");
    return rank;
}

Fnv fnv_of(const std::vector<u32>& masks) {
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3,
                "usage: top-pair-event-referee REPLAY_BAND THM4266_RESULTS");
        const Carrier carrier = reconstruct_carrier(argv[1], argv[2]);
        const LiteralWall wall = build_event_sweep(520, 688);
        require(wall.grid == INT64_C(784369854908640),
                "literal common denominator changed");
        require(static_cast<i128>(wall.pair_ticks) * 19565 ==
                    static_cast<i128>(wall.grid) * 14374,
                "pair-safe mass disagrees with primitive quotient");

        u64 carrier_equalities = 0;
        std::vector<u32> active = active_subdeck(
            carrier.masks, wall, carrier_equalities);
        require(active.size() == 5934 && carrier_equalities == 0,
                "active inherited subdeck changed");
        const Fnv active_ordered_fnv = fnv_of(active);
        std::vector<u32> active_sorted = active;
        std::sort(active_sorted.begin(), active_sorted.end());
        const Fnv active_sorted_fnv = fnv_of(active_sorted);

        std::vector<int> all_bits(30);
        std::iota(all_bits.begin(), all_bits.end(), 0);
        std::vector<u32> bodies;
        bodies.reserve(BODY_COUNT);
        generate(all_bits, 0, 9, 0, bodies);
        require(bodies.size() == BODY_COUNT, "body universe count changed");
        const BodyScan inherited = scan_bodies(bodies, active);
        require(inherited.bodies == BODY_COUNT &&
                    inherited.failures == std::vector<u32>(
                        {EXPECTED_FIRST_FAILURE, EXPECTED_SECOND_FAILURE}),
                "inherited missed-body set changed");
        const Fnv failure_fnv = fnv_of(inherited.failures);

        const u32 common_union = inherited.failures[0] | inherited.failures[1];
        const u32 common_intersection = inherited.failures[0] & inherited.failures[1];
        std::vector<int> permitted_bits;
        for (int bit = 0; bit < 30; ++bit) {
            if ((common_union & (u32{1} << bit)) == 0) {
                permitted_bits.push_back(bit);
            }
        }
        require(permitted_bits.size() == 19 && std::popcount(common_union) == 11 &&
                    std::popcount(common_intersection) == 7,
                "missed-body overlap geometry changed");
        std::vector<u32> common_candidates;
        generate(permitted_bits, 0, 8, 0, common_candidates);
        require(common_candidates.size() == 75582, "common candidate count changed");

        std::vector<u32> common_active;
        u64 common_equalities = 0;
        i128 least_margin = 0;
        i128 greatest_margin = 0;
        u32 least_margin_mask = 0;
        u32 greatest_margin_mask = 0;
        bool first_active = true;
        for (u32 repair : common_candidates) {
            const i128 margin = static_cast<i128>(63) * mass_ticks(wall, repair) -
                                static_cast<i128>(4) * wall.grid;
            if (margin == 0) ++common_equalities;
            if (margin < 0) continue;
            common_active.push_back(repair);
            if (first_active || margin < least_margin ||
                (margin == least_margin && repair < least_margin_mask)) {
                least_margin = margin;
                least_margin_mask = repair;
            }
            if (first_active || margin > greatest_margin ||
                (margin == greatest_margin && repair < greatest_margin_mask)) {
                greatest_margin = margin;
                greatest_margin_mask = repair;
            }
            first_active = false;
        }
        std::sort(common_active.begin(), common_active.end());
        require(common_active.size() == 72 && common_equalities == 0 &&
                    common_active.front() == EXPECTED_REPAIR,
                "common-active repair classification changed");
        const Fnv common_fnv = fnv_of(common_active);
        require(common_fnv.value == UINT64_C(0xed1bfbaf6eaa06a3),
                "common-active repair FNV changed");
        std::set<u32> carrier_set(carrier.masks.begin(), carrier.masks.end());
        u64 carrier_overlap = 0;
        for (u32 repair : common_active) {
            if (carrier_set.contains(repair)) ++carrier_overlap;
        }
        require(carrier_overlap == 0, "common-active repair already in carrier");

        const i64 repair_mass = mass_ticks(wall, EXPECTED_REPAIR);
        const i128 repair_margin = static_cast<i128>(63) * repair_mass -
                                   static_cast<i128>(4) * wall.grid;
        require(repair_margin > 0 &&
                    fraction(repair_mass, wall.grid) ==
                        "1559831620541/24511557965895" &&
                    colex_rank0(EXPECTED_REPAIR) == 51083,
                "least common repair exact mass or rank changed");

        std::vector<u32> augmented = active;
        augmented.push_back(EXPECTED_REPAIR);
        const BodyScan repaired = scan_bodies(bodies, augmented);
        require(repaired.bodies == BODY_COUNT && repaired.failures.empty(),
                "one-atom augmented deck did not close all bodies");

        std::cout << "LRC14_ENDPOINT_520_688_MINIMUM_ONE_ATOM_EVENT_SWEEP_REFEREE_THM4278\n";
        std::cout << "PAIR 520,688 GRID " << wall.grid
                  << " CELLS " << wall.cells
                  << " PAIR_SAFE_CELLS " << wall.pair_safe_cells
                  << " PATTERNS " << wall.weighted_patterns.size()
                  << " PAIR_TICKS " << wall.pair_ticks
                  << " PAIR_MASS " << fraction(wall.pair_ticks, wall.grid)
                  << "\n";
        std::cout << "CARRIER_STAGES";
        for (std::size_t index = 0; index < carrier.stage_count.size(); ++index) {
            std::cout << ' ' << carrier.stage_count[index] << '/'
                      << std::hex << carrier.stage_fnv[index] << std::dec;
        }
        std::cout << "\n";
        std::cout << "CARRIER " << carrier.masks.size()
                  << " FNV " << std::hex << carrier.stage_fnv[4] << std::dec
                  << " ACTIVE " << active.size()
                  << " EQUALITIES " << carrier_equalities
                  << " ACTIVE_ORDER_FNV " << std::hex << active_ordered_fnv.value
                  << " ACTIVE_SORTED_FNV " << active_sorted_fnv.value << std::dec
                  << "\n";
        std::cout << "INHERITED_BODY_SCAN BODIES " << inherited.bodies
                  << " FAILURES " << inherited.failures.size()
                  << " CHECKS " << inherited.checks
                  << " MAX_CHECKS " << inherited.max_checks
                  << " FAILURE_FNV " << std::hex << failure_fnv.value << std::dec
                  << "\n";
        for (std::size_t index = 0; index < inherited.failures.size(); ++index) {
            std::cout << "FAILURE_" << index << " MASK " << std::hex
                      << inherited.failures[index] << std::dec << " LABELS {"
                      << labels(inherited.failures[index]) << "}\n";
        }
        std::cout << "FAILURE_INTERSECTION SIZE "
                  << std::popcount(common_intersection) << " MASK " << std::hex
                  << common_intersection << std::dec << " LABELS {"
                  << labels(common_intersection) << "} UNION SIZE "
                  << std::popcount(common_union) << " MASK " << std::hex
                  << common_union << std::dec << " LABELS {"
                  << labels(common_union) << "}\n";
        std::cout << "COMMON_DISJOINT PERMITTED " << permitted_bits.size()
                  << " CANDIDATES " << common_candidates.size()
                  << " ACTIVE " << common_active.size()
                  << " EQUALITIES " << common_equalities
                  << " CARRIER_OVERLAP " << carrier_overlap
                  << " ACTIVE_SORTED_FNV " << std::hex << common_fnv.value
                  << std::dec << "\n";
        std::cout << "LEAST_REPAIR MASK " << std::hex << EXPECTED_REPAIR
                  << std::dec << " LABELS {" << labels(EXPECTED_REPAIR)
                  << "} COLEX_RANK0 " << colex_rank0(EXPECTED_REPAIR)
                  << " MASS_TICKS " << repair_mass
                  << " MASS " << fraction(repair_mass, wall.grid)
                  << " TARGET 4/63 MARGIN_TICKS63 " << decimal(repair_margin)
                  << " MARGIN_REDUCED "
                  << fraction(repair_margin, static_cast<i128>(63) * wall.grid)
                  << "\n";
        std::cout << "ACTIVE_MARGIN_RANGE MIN " << decimal(least_margin)
                  << " MASK " << std::hex << least_margin_mask << std::dec
                  << " LABELS {" << labels(least_margin_mask) << "} MAX "
                  << decimal(greatest_margin) << " MASK " << std::hex
                  << greatest_margin_mask << std::dec << " LABELS {"
                  << labels(greatest_margin_mask) << "}\n";
        std::cout << "AUGMENTED_BODY_SCAN ACTIVE " << augmented.size()
                  << " BODIES " << repaired.bodies
                  << " FAILURES " << repaired.failures.size()
                  << " CHECKS " << repaired.checks
                  << " MAX_CHECKS " << repaired.max_checks << "\n";
        std::cout << "COMMON_ACTIVE_MASKS_HEX";
        for (u32 repair : common_active) {
            std::cout << ' ' << std::hex << repair;
        }
        std::cout << std::dec << "\n";
        std::cout << "STATUS FINITE_EXACT INDEPENDENT_EVENT_SWEEP PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "REFEREE_ERROR " << error.what() << '\n';
        return 1;
    }
}
