// Finite-exact direct literal audit for every body missed by the unchanged
// THM-4318 carrier on a runtime-selected endpoint.  The truncated mass uses
// only pair-safe
// wall cells whose pool-failure set has rank at most nine, hence is a lower
// bound on the complete literal safe mass.

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
int gEndpoint = 0;
constexpr u64 kEmptyFnv = UINT64_C(0xcbf29ce484222325);

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
    std::ostringstream output;
    output << std::hex << std::setw(8) << std::setfill('0') << value;
    return output.str();
}

i64 checked_lcm(i64 left, i64 right) {
    const i64 divisor = std::gcd(left, right);
    const i128 result = static_cast<i128>(left / divisor) * right;
    require(result <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(result);
}

i64 base_grid() {
    i64 grid = 1;
    for (int speed : kPool) grid = checked_lcm(grid, 14LL * speed);
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
    u64 pair_safe_cells = 0;
    i64 pair_ticks = 0;
    std::map<u32, i64> all_classes;
    std::map<u32, i64> low_classes;
    u64 all_class_fnv = 0;
    u64 low_class_fnv = 0;
};

Geometry build_geometry(int q) {
    Geometry geometry;
    geometry.grid = checked_lcm(base_grid(), 14LL * q);
    geometry.grid = checked_lcm(geometry.grid, 14LL * gEndpoint);
    std::vector<i64> walls = {0, geometry.grid};
    auto add_walls = [&](int speed) {
        const i64 divisor = 14LL * speed;
        require(geometry.grid % divisor == 0, "nonintegral wall unit");
        const i64 unit = geometry.grid / divisor;
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(q);
    add_walls(gEndpoint);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    geometry.cells = walls.size() - 1;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(q, geometry.grid, left, right) ||
            !safe_midpoint(gEndpoint, geometry.grid, left, right))
            continue;
        ++geometry.pair_safe_cells;
        const i64 width = right - left;
        geometry.pair_ticks += width;
        u32 failure = 0;
        for (unsigned bit = 0; bit < kPool.size(); ++bit)
            if (!safe_midpoint(kPool[bit], geometry.grid, left, right))
                failure |= UINT32_C(1) << bit;
        geometry.all_classes[failure] += width;
        if (std::popcount(failure) <= 9)
            geometry.low_classes[failure] += width;
    }
    Fnv all_ledger, low_ledger;
    for (const auto& [mask, width] : geometry.all_classes) {
        all_ledger.add(mask); all_ledger.add(width);
    }
    for (const auto& [mask, width] : geometry.low_classes) {
        low_ledger.add(mask); low_ledger.add(width);
    }
    geometry.all_class_fnv = all_ledger.state;
    geometry.low_class_fnv = low_ledger.state;
    return geometry;
}

struct Failures {
    std::map<int, std::vector<u32>> rows;
    std::map<int, u64> body_fnv;
    u64 global_fnv = 0;
    u64 count = 0;
};

Failures read_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "failure header changed");
    Failures result;
    std::map<int, std::set<u32>> distinct;
    std::map<int, Fnv> row_ledgers;
    Fnv global;
    int previous_q = -1;
    u32 previous_body = 0;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q = 0, r = 0;
        std::string token, extra;
        fields >> q >> r >> token;
        require(fields && !(fields >> extra) && r == gEndpoint,
                "malformed or wrong-endpoint failure row");
        const u64 wide = std::stoull(token, nullptr, 16);
        require(wide < (UINT64_C(1) << 30), "body escaped label universe");
        const u32 body = static_cast<u32>(wide);
        require(std::popcount(body) == 9 && distinct[q].insert(body).second,
                "body rank/distinctness changed");
        require(q > previous_q || (q == previous_q && body > previous_body),
                "failure ledger ordering changed");
        previous_q = q; previous_body = body;
        result.rows[q].push_back(body);
        row_ledgers[q].add(body);
        global.add(q); global.add(r); global.add(body);
        ++result.count;
    }
    require(input.eof() && !result.rows.empty(), "empty/incomplete failure ledger");
    for (const auto& [q, ledger] : row_ledgers)
        result.body_fnv[q] = ledger.state;
    result.global_fnv = global.state;
    return result;
}

struct Summary {
    i128 low_min = 0;
    i128 low_max = 0;
    i128 full_min = 0;
    i128 full_max = 0;
    u32 low_min_body = 0;
    u32 low_max_body = 0;
    u32 full_min_body = 0;
    u32 full_max_body = 0;
    u64 low_positive = 0;
    u64 full_positive = 0;
    u64 equality = 0;
    u64 ledger_fnv = 0;
};
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 4, "usage: generic_endpoint_direct ENDPOINT FAILURES DETAIL_CSV");
        gEndpoint = std::stoi(argv[1]);
        require(gEndpoint > 0, "invalid endpoint");
        const Failures failures = read_failures(argv[2]);
        std::ofstream detail(argv[3], std::ios::binary);
        require(static_cast<bool>(detail), "cannot create detail ledger");
        detail << "q,r,ordinal,body_hex,truncated_mass,full_mass,truncated_scaled_ticks,full_scaled_ticks\n";
        std::map<int, Geometry> geometries;
        std::map<int, Summary> summaries;
        for (const auto& [q, bodies] : failures.rows) {
            Geometry geometry = build_geometry(q);
            Summary summary;
            Fnv ledger;
            for (std::size_t ordinal = 0; ordinal < bodies.size(); ++ordinal) {
                const u32 body = bodies[ordinal];
                i64 low_mass = 0;
                i64 full_mass = 0;
                for (const auto& [failure, width] : geometry.low_classes)
                    if ((failure & body) == 0) low_mass += width;
                for (const auto& [failure, width] : geometry.all_classes)
                    if ((failure & body) == 0) full_mass += width;
                require(low_mass <= full_mass, "lower-bound direction failed");
                const i128 low_ticks = static_cast<i128>(63) * low_mass -
                                       static_cast<i128>(4) * geometry.grid;
                const i128 full_ticks = static_cast<i128>(63) * full_mass -
                                        static_cast<i128>(4) * geometry.grid;
                detail << q << ',' << gEndpoint << ',' << ordinal << ','
                       << hex8(body) << ',' << low_mass << ',' << full_mass
                       << ',' << decimal(low_ticks) << ',' << decimal(full_ticks)
                       << '\n';
                ledger.add(q); ledger.add(gEndpoint); ledger.add(ordinal);
                ledger.add(body); ledger.add(low_mass); ledger.add(full_mass);
                add_i128(ledger, low_ticks); add_i128(ledger, full_ticks);
                if (ordinal == 0 || low_ticks < summary.low_min ||
                    (low_ticks == summary.low_min && body < summary.low_min_body)) {
                    summary.low_min = low_ticks; summary.low_min_body = body;
                }
                if (ordinal == 0 || low_ticks > summary.low_max ||
                    (low_ticks == summary.low_max && body < summary.low_max_body)) {
                    summary.low_max = low_ticks; summary.low_max_body = body;
                }
                if (ordinal == 0 || full_ticks < summary.full_min ||
                    (full_ticks == summary.full_min && body < summary.full_min_body)) {
                    summary.full_min = full_ticks; summary.full_min_body = body;
                }
                if (ordinal == 0 || full_ticks > summary.full_max ||
                    (full_ticks == summary.full_max && body < summary.full_max_body)) {
                    summary.full_max = full_ticks; summary.full_max_body = body;
                }
                summary.low_positive += low_ticks > 0;
                summary.full_positive += full_ticks > 0;
                summary.equality += low_mass == full_mass;
            }
            summary.ledger_fnv = ledger.state;
            geometries.emplace(q, std::move(geometry));
            summaries.emplace(q, summary);
        }
        require(detail.good(), "detail ledger write failed");

        u64 low_positive = 0, full_positive = 0, equality = 0;
        Fnv summary_ledger;
        std::cout << "LRC14_GENERIC_ENDPOINT_DIRECT_LITERAL_PRIMARY_V1\n"
                  << "ENDPOINT " << gEndpoint << '\n'
                  << "FAILURE_ROWS " << failures.rows.size() << " BODIES "
                  << failures.count << " FAILURE_FNV " << std::hex
                  << failures.global_fnv << std::dec << '\n';
        for (const auto& [q, bodies] : failures.rows) {
            const Geometry& geometry = geometries.at(q);
            const Summary& summary = summaries.at(q);
            low_positive += summary.low_positive;
            full_positive += summary.full_positive;
            equality += summary.equality;
            summary_ledger.add(q); summary_ledger.add(bodies.size());
            summary_ledger.add(geometry.grid); summary_ledger.add(geometry.cells);
            summary_ledger.add(geometry.pair_safe_cells);
            summary_ledger.add(geometry.pair_ticks);
            summary_ledger.add(geometry.all_classes.size());
            summary_ledger.add(geometry.low_classes.size());
            summary_ledger.add(geometry.all_class_fnv);
            summary_ledger.add(geometry.low_class_fnv);
            add_i128(summary_ledger, summary.low_min);
            summary_ledger.add(summary.low_min_body);
            add_i128(summary_ledger, summary.full_min);
            summary_ledger.add(summary.full_min_body);
            summary_ledger.add(summary.low_positive);
            summary_ledger.add(summary.full_positive);
            summary_ledger.add(summary.equality);
            summary_ledger.add(summary.ledger_fnv);
            std::cout << "ROW " << q << ',' << gEndpoint << " BODIES "
                      << bodies.size() << " BODY_FNV " << std::hex
                      << failures.body_fnv.at(q) << std::dec << " GRID "
                      << geometry.grid << " CELLS " << geometry.cells
                      << " PAIR_SAFE_CELLS " << geometry.pair_safe_cells
                      << " PAIR_TICKS " << geometry.pair_ticks
                      << " ALL_CLASSES " << geometry.all_classes.size()
                      << " LOW_CLASSES " << geometry.low_classes.size()
                      << " ALL_CLASS_FNV " << std::hex
                      << geometry.all_class_fnv << " LOW_CLASS_FNV "
                      << geometry.low_class_fnv << std::dec
                      << " LOW_POSITIVE " << summary.low_positive
                      << " FULL_POSITIVE " << summary.full_positive
                      << " EQUALITY " << summary.equality
                      << " LOW_MIN_TICKS " << decimal(summary.low_min)
                      << " LOW_MIN_BODY " << hex8(summary.low_min_body)
                      << " FULL_MIN_TICKS " << decimal(summary.full_min)
                      << " FULL_MIN_BODY " << hex8(summary.full_min_body)
                      << " LOW_MAX_TICKS " << decimal(summary.low_max)
                      << " LOW_MAX_BODY " << hex8(summary.low_max_body)
                      << " FULL_MAX_TICKS " << decimal(summary.full_max)
                      << " FULL_MAX_BODY " << hex8(summary.full_max_body)
                      << " DETAIL_FNV " << std::hex << summary.ledger_fnv
                      << std::dec << '\n';
        }
        std::cout << "SUMMARY BODIES " << failures.count
                  << " LOW_POSITIVE " << low_positive
                  << " FULL_POSITIVE " << full_positive
                  << " EQUALITY " << equality << " SUMMARY_FNV "
                  << std::hex << summary_ledger.state << std::dec << '\n'
                  << "INEQUALITY TRUNCATED_LOW_RANK_CELL_MASS_LE_LITERAL_SAFE_MASS\n"
                  << "SCOPE FINITE_EXACT_FIXED_POOL_SINGLE_ENDPOINT_FAILURE_BODIES_"
                     "NO_CARRIER_EXCHANGE_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT "
                  << (low_positive == failures.count ? "PASS" : "HOSTILE_FAIL")
                  << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "GENERIC_ENDPOINT_DIRECT_LITERAL_PRIMARY_ERROR "
                  << error.what() << '\n';
        return 1;
    }
}
