// Import-free exact literal-wall certificate for the first obstruction after
// the endpoint-636 carrier exchange.  Two independent enumerations prove that
// body 1d106401 has no active rank-eight responder at (256,632).  A ten-body
// activity-incompatible packing and an explicit ten-mask cover prove that the
// other 65 failed bodies have exact repair-cover number ten.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr std::array<int, 30> kPool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr i64 kFixedGrid = INT64_C(18241159416480);
constexpr u32 kHostile = UINT32_C(0x1d106401);
constexpr u64 kFailuresFnv = UINT64_C(0xbde16ce2f1f9043e);
constexpr u64 kEmptyFnv = UINT64_C(0xcbf29ce484222325);
constexpr u64 kRank8Count = UINT64_C(5852925);
constexpr u64 kHostileCandidates = UINT64_C(203490);
constexpr std::array<std::size_t, 10> kPacking = {
    6, 7, 9, 14, 24, 28, 29, 41, 52, 58};
constexpr std::array<u32, 10> kWitness = {
    0x201e8088, 0x03448124, 0x12a01206, 0x282c1048, 0x22285204,
    0x02201a85, 0x2098012c, 0x0024b209, 0x31481280, 0x016490c0};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

struct Fnv {
    u64 state = kEmptyFnv;
    void add(u64 word) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (word >> (8 * byte)) & UINT64_C(0xff);
            state *= UINT64_C(0x100000001b3);
        }
    }
};

void add_i128(Fnv& ledger, i128 value) {
    const u128 bits = static_cast<u128>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    u128 magnitude = static_cast<u128>(value);
    if (negative) magnitude = u128{0} - magnitude;
    std::string answer;
    while (magnitude != 0) {
        answer.push_back(static_cast<char>('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

std::string hex8(u32 value) {
    std::ostringstream out;
    out << std::hex << std::setw(8) << std::setfill('0') << value;
    return out.str();
}

std::string labels(u32 mask) {
    std::ostringstream out;
    bool first = true;
    for (unsigned bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        if (!first) out << ',';
        out << bit;
        first = false;
    }
    return out.str();
}

i128 abs128(i128 value) {
    return value < 0 ? -value : value;
}

i128 gcd128(i128 left, i128 right) {
    left = abs128(left);
    right = abs128(right);
    while (right != 0) {
        const i128 remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

i64 checked_lcm(i64 left, i64 right) {
    require(left > 0 && right > 0, "nonpositive LCM input");
    const i128 wide = static_cast<i128>(left / std::gcd(left, right)) * right;
    require(wide <= std::numeric_limits<i64>::max(), "grid overflow");
    return static_cast<i64>(wide);
}

i64 fixed_grid() {
    i64 grid = 1;
    for (int speed : kPool) grid = checked_lcm(grid, 14LL * speed);
    require(grid == kFixedGrid, "fixed grid changed");
    return grid;
}

bool safe_midpoint(int speed, i64 grid, i64 left, i64 right) {
    i128 residue = static_cast<i128>(speed) *
                   (static_cast<i128>(left) + right);
    residue %= static_cast<i128>(2) * grid;
    if (residue < 0) residue += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * residue >= grid &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * grid;
}

struct Geometry {
    i64 grid = 0;
    u64 cells = 0;
    i64 pair_ticks = 0;
    std::vector<std::pair<u32, i64>> classes;
};

Geometry build_geometry(int pair_q = 256, int pair_r = 632) {
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
        for (unsigned vertex = 0; vertex < kPool.size(); ++vertex)
            if (!safe_midpoint(kPool[vertex], grid, left, right))
                failure |= u32{1} << vertex;
        const i64 width = right - left;
        geometry.pair_ticks += width;
        by_failure[failure] += width;
    }
    for (const auto& entry : by_failure)
        if (std::popcount(entry.first) <= 9) geometry.classes.push_back(entry);
    return geometry;
}

struct Margin {
    i64 mass = 0;
    i128 ticks = 0;
};

Margin margin(const Geometry& geometry, u32 mask) {
    i64 mass = 0;
    for (const auto& [failure, width] : geometry.classes)
        if ((failure & ~mask) == 0) mass += width;
    return {mass, static_cast<i128>(63) * mass -
                      static_cast<i128>(4) * geometry.grid};
}

std::vector<u32> read_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,body_hex",
            "failure header changed");
    std::vector<u32> bodies;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        int q = 0;
        int r = 0;
        std::string token;
        row >> q >> r >> token;
        require(row && r == 632, "malformed failure row");
        if (q != 256) continue;
        const u32 body = static_cast<u32>(std::stoul(token, nullptr, 16));
        require(std::popcount(body) == 9, "body rank changed");
        bodies.push_back(body);
    }
    Fnv ledger;
    for (u32 body : bodies) ledger.add(body);
    require(bodies.size() == 66 && ledger.state == kFailuresFnv,
            "(256,632) failure ledger changed");
    require(bodies[48] == kHostile, "hostile ordinal changed");
    return bodies;
}

u32 next_combination(u32 value) {
    const u32 low = value & (u32{0} - value);
    require(low != 0, "zero combination low bit");
    const u32 ripple = value + low;
    return ripple | (((ripple ^ value) >> 2) / low);
}

template <class Callback>
void choose_rec(const std::vector<unsigned char>& positions,
                std::size_t start, unsigned chosen, u32 mask,
                Callback& callback) {
    if (chosen == 8) {
        callback(mask);
        return;
    }
    const std::size_t needed = 8 - chosen;
    for (std::size_t index = start;
         index + needed <= positions.size(); ++index)
        choose_rec(positions, index + 1, chosen + 1,
                   mask | (u32{1} << positions[index]), callback);
}

template <class Callback>
u64 enumerate_complement(u32 forbidden, Callback callback) {
    std::vector<unsigned char> positions;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((forbidden & (u32{1} << bit)) == 0)
            positions.push_back(static_cast<unsigned char>(bit));
    if (positions.size() < 8) return 0;
    u64 count = 0;
    auto counted = [&](u32 mask) {
        require(std::popcount(mask) == 8 && (mask & forbidden) == 0,
                "complement enumeration escaped");
        callback(mask);
        ++count;
    };
    choose_rec(positions, 0, 0, 0, counted);
    return count;
}

struct Audit {
    u64 candidates = 0;
    u64 active = 0;
    bool have = false;
    Margin maximum;
    u32 least_maximum = 0;
    u64 maximum_ties = 0;
    Margin minimum;
    u32 least_minimum = 0;
    Fnv ledger;
    std::vector<u32> masks;
};

void consume(Audit& audit, const Geometry& geometry, u32 mask) {
    const Margin value = margin(geometry, mask);
    ++audit.candidates;
    audit.active += value.ticks >= 0;
    audit.ledger.add(mask);
    audit.ledger.add(static_cast<u64>(value.mass));
    add_i128(audit.ledger, value.ticks);
    audit.masks.push_back(mask);
    if (!audit.have || value.ticks > audit.maximum.ticks) {
        audit.maximum = value;
        audit.least_maximum = mask;
        audit.maximum_ties = 1;
    } else if (value.ticks == audit.maximum.ticks) {
        ++audit.maximum_ties;
        audit.least_maximum = std::min(audit.least_maximum, mask);
    }
    if (!audit.have || value.ticks < audit.minimum.ticks ||
        (value.ticks == audit.minimum.ticks && mask < audit.least_minimum)) {
        audit.minimum = value;
        audit.least_minimum = mask;
    }
    audit.have = true;
}

u64 choose(unsigned n, unsigned k) {
    if (k > n) return 0;
    if (k > n - k) k = n - k;
    u64 answer = 1;
    for (unsigned i = 1; i <= k; ++i)
        answer = answer * (n - k + i) / i;
    return answer;
}

}  // namespace


int main(int argc, char** argv) {
    try {
        require(argc == 3,
                "usage: detached R632_FAILURES RESPONSE_CSV");
        const std::vector<u32> failures = read_failures(argv[1]);
        const Geometry geometry = build_geometry();
        std::cout << "R632_DETACHED_HOSTILE_SURVIVOR_V1\n"
                  << "PAIR 256,632 GRID " << geometry.grid << " CELLS "
                  << geometry.cells << " PAIR_TICKS " << geometry.pair_ticks
                  << " FAILURE_CLASSES_RANK_LE9 " << geometry.classes.size()
                  << '\n'
                  << "FAILURES 66 FNV " << std::hex << kFailuresFnv
                  << std::dec << " HOSTILE_ORDINAL 48 HOSTILE "
                  << hex8(kHostile) << " LABELS {" << labels(kHostile)
                  << "}\n";

        Audit complement;
        const u64 complement_count = enumerate_complement(
            kHostile, [&](u32 mask) { consume(complement, geometry, mask); });
        require(complement_count == kHostileCandidates &&
                    complement.candidates == kHostileCandidates,
                "complement candidate count changed");

        Audit full_scan;
        u64 universe = 0;
        const u32 limit = u32{1} << 30;
        for (u32 mask = (u32{1} << 8) - 1; mask < limit;
             mask = next_combination(mask)) {
            ++universe;
            if ((mask & kHostile) == 0) consume(full_scan, geometry, mask);
        }
        require(universe == kRank8Count &&
                    full_scan.candidates == kHostileCandidates,
                "full rank-eight scan changed");
        std::sort(complement.masks.begin(), complement.masks.end());
        std::sort(full_scan.masks.begin(), full_scan.masks.end());
        require(complement.masks == full_scan.masks,
                "two hostile candidate sets disagree");
        require(complement.active == 0 && full_scan.active == 0 &&
                    complement.maximum.ticks < 0 &&
                    complement.maximum.ticks == full_scan.maximum.ticks &&
                    complement.maximum.mass == full_scan.maximum.mass &&
                    complement.least_maximum == full_scan.least_maximum &&
                    complement.maximum_ties == 1 &&
                    full_scan.maximum_ties == 1,
                "hostile maximum changed");
        const i128 surplus_denominator =
            static_cast<i128>(63) * geometry.grid;
        const i128 divisor =
            gcd128(complement.maximum.ticks, surplus_denominator);
        std::cout << "HOSTILE_COMPLEMENT CANDIDATES "
                  << complement.candidates << " ACTIVE " << complement.active
                  << " LEDGER_FNV " << std::hex << complement.ledger.state
                  << std::dec << '\n'
                  << "HOSTILE_FULL_SCAN UNIVERSE " << universe
                  << " CANDIDATES " << full_scan.candidates << " ACTIVE "
                  << full_scan.active << " LEDGER_FNV " << std::hex
                  << full_scan.ledger.state << std::dec << '\n'
                  << "MAX_MASS " << complement.maximum.mass
                  << " MAX_MARGIN_TICKS63 "
                  << decimal(complement.maximum.ticks)
                  << " SURPLUS_DENOMINATOR_63GRID "
                  << decimal(surplus_denominator)
                  << " REDUCED_SURPLUS "
                  << decimal(complement.maximum.ticks / divisor) << '/'
                  << decimal(surplus_denominator / divisor)
                  << " UNIQUE_MAX_MASK " << hex8(complement.least_maximum)
                  << " LABELS {" << labels(complement.least_maximum) << "}"
                  << " MAX_TIES " << complement.maximum_ties << '\n';

        // Rank-sharpness of the obstruction: the best disjoint rank-eight
        // mask has thirteen one-label extensions that remain disjoint from
        // the hostile body.  Check the complete extension fibre literally.
        u64 extension_count = 0;
        u64 strictly_active_extensions = 0;
        u32 least_extension = 0;
        u32 least_margin_mask = 0;
        Margin least_extension_margin{};
        bool have_extension = false;
        Fnv extension_ledger;
        for (unsigned bit = 0; bit < 30; ++bit) {
            const u32 label = u32{1} << bit;
            if (((complement.least_maximum | kHostile) & label) != 0)
                continue;
            const u32 extension = complement.least_maximum | label;
            require(std::popcount(extension) == 9 &&
                        (extension & kHostile) == 0,
                    "rank-nine extension escaped hostile complement");
            const Margin value = margin(geometry, extension);
            ++extension_count;
            strictly_active_extensions += value.ticks > 0;
            if (!have_extension || extension < least_extension)
                least_extension = extension;
            if (!have_extension || value.ticks < least_extension_margin.ticks ||
                (value.ticks == least_extension_margin.ticks &&
                 extension < least_margin_mask)) {
                least_extension_margin = value;
                least_margin_mask = extension;
            }
            have_extension = true;
            extension_ledger.add(extension);
            extension_ledger.add(static_cast<u64>(value.mass));
            add_i128(extension_ledger, value.ticks);
        }
        require(have_extension && extension_count == 13 &&
                    strictly_active_extensions == extension_count &&
                    least_extension_margin.ticks > 0,
                "rank-nine extension fibre is not strictly active");
        std::cout << "RANK9_ONE_LABEL_EXTENSIONS COUNT " << extension_count
                  << " STRICTLY_ACTIVE " << strictly_active_extensions
                  << " LEAST_MASK " << hex8(least_extension)
                  << " MIN_MARGIN_TICKS63 "
                  << decimal(least_extension_margin.ticks)
                  << " MIN_MARGIN_MASK " << hex8(least_margin_mask)
                  << " LEDGER_FNV " << std::hex << extension_ledger.state
                  << std::dec << '\n';

        Fnv packing_ledger;
        u64 packing_candidate_checks = 0;
        i128 closest_pair_margin = 0;
        u32 closest_pair_mask = 0;
        std::size_t closest_left = 0;
        std::size_t closest_right = 0;
        bool have_pair = false;
        for (std::size_t left = 0; left < kPacking.size(); ++left) {
            for (std::size_t right = left + 1; right < kPacking.size(); ++right) {
                const u32 united = failures[kPacking[left]] |
                                   failures[kPacking[right]];
                Audit pair_audit;
                const u64 count = enumerate_complement(
                    united,
                    [&](u32 mask) { consume(pair_audit, geometry, mask); });
                require(count == choose(30 - std::popcount(united), 8) &&
                            pair_audit.active == 0 &&
                            pair_audit.maximum.ticks < 0,
                        "packing pair has an active common responder");
                packing_candidate_checks += count;
                packing_ledger.add(kPacking[left]);
                packing_ledger.add(kPacking[right]);
                packing_ledger.add(failures[kPacking[left]]);
                packing_ledger.add(failures[kPacking[right]]);
                packing_ledger.add(count);
                add_i128(packing_ledger, pair_audit.maximum.ticks);
                packing_ledger.add(pair_audit.least_maximum);
                if (!have_pair || pair_audit.maximum.ticks > closest_pair_margin) {
                    have_pair = true;
                    closest_pair_margin = pair_audit.maximum.ticks;
                    closest_pair_mask = pair_audit.least_maximum;
                    closest_left = left;
                    closest_right = right;
                }
            }
        }
        std::cout << "PACKING SIZE " << kPacking.size()
                  << " PAIRS 45 CANDIDATE_CHECKS " << packing_candidate_checks
                  << " LEDGER_FNV " << std::hex << packing_ledger.state
                  << std::dec << " CLOSEST_PAIR " << kPacking[closest_left]
                  << ',' << kPacking[closest_right] << " BODIES "
                  << hex8(failures[kPacking[closest_left]]) << ','
                  << hex8(failures[kPacking[closest_right]])
                  << " MAX_MARGIN_TICKS63 " << decimal(closest_pair_margin)
                  << " MASK " << hex8(closest_pair_mask) << '\n';

        std::ofstream responses(argv[2]);
        require(static_cast<bool>(responses), "cannot create response ledger");
        responses << "mask_hex,body_ordinal,body_hex,margin_ticks63,grid\n";
        std::vector<unsigned char> covered(failures.size(), 0);
        Fnv witness_ledger;
        Fnv response_ledger;
        std::set<u32> distinct_witness;
        u64 incidences = 0;
        for (u32 mask : kWitness) {
            require(std::popcount(mask) == 8 &&
                        distinct_witness.insert(mask).second,
                    "invalid witness mask");
            const Margin value = margin(geometry, mask);
            require(value.ticks > 0, "witness mask is not strictly active");
            witness_ledger.add(mask);
            witness_ledger.add(static_cast<u64>(value.mass));
            add_i128(witness_ledger, value.ticks);
            u64 hits = 0;
            u64 packing_hits = 0;
            for (std::size_t index = 0; index < failures.size(); ++index) {
                if ((mask & failures[index]) != 0) continue;
                require(index != 48, "witness hits hostile body");
                covered[index] = 1;
                ++hits;
                ++incidences;
                packing_hits +=
                    std::find(kPacking.begin(), kPacking.end(), index) !=
                    kPacking.end();
                response_ledger.add(mask);
                response_ledger.add(index);
                response_ledger.add(failures[index]);
                responses << hex8(mask) << ',' << index << ','
                          << hex8(failures[index]) << ','
                          << decimal(value.ticks) << ',' << geometry.grid
                          << '\n';
            }
            require(packing_hits == 1,
                    "witness does not pay exactly one packing obligation");
            std::cout << "WITNESS MASK " << hex8(mask) << " LABELS {"
                      << labels(mask) << "} MASS " << value.mass
                      << " MARGIN_TICKS63 " << decimal(value.ticks)
                      << " HITS " << hits << " PACKING_HITS "
                      << packing_hits << '\n';
        }
        require(responses.good(), "response-ledger write failed");
        u64 covered_count = 0;
        for (std::size_t index = 0; index < covered.size(); ++index) {
            if (index == 48) require(!covered[index], "hostile marked covered");
            else require(covered[index], "reachable survivor not covered");
            covered_count += covered[index];
        }
        std::cout << "SURVIVOR_WITNESS MASKS " << kWitness.size()
                  << " COVERED " << covered_count << " MISSED 1 HOSTILE "
                  << hex8(kHostile) << " INCIDENCES " << incidences
                  << " WITNESS_FNV " << std::hex << witness_ledger.state
                  << " RESPONSE_FNV " << response_ledger.state << std::dec
                  << '\n'
                  << "SURVIVOR_MINIMUM 10 LOWER activity-incompatible-packing"
                     " UPPER explicit-ten-mask-cover\n"
                  << "SCOPE IMPORT_FREE_LITERAL_WALL_FIXED_PAIR_LABELLED_"
                     "RANK8_NO_HIGHER_RANK_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS UNIQUE_UNREPAIRABLE_HOSTILE_AND_EXACT_"
                     "TEN_MASK_SURVIVOR_COVER\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "R632_DETACHED_ERROR " << error.what() << '\n';
        return 1;
    }
}
