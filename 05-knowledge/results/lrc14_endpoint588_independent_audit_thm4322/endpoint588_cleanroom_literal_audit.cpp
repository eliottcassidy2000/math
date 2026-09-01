// Clean-room exact wall-cell audit for every endpoint-588 carrier failure.
// It does not import any maintained or scout implementation.

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
constexpr u64 kFnvBasis = UINT64_C(0xcbf29ce484222325);
constexpr u64 kFnvPrime = UINT64_C(0x100000001b3);

[[noreturn]] void die(const std::string& message) {
    throw std::runtime_error(message);
}
void check(bool condition, const std::string& message) {
    if (!condition) die(message);
}
struct Fnv {
    u64 state = kFnvBasis;
    void add(u64 value) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (value >> (8 * byte)) & UINT64_C(255);
            state *= kFnvPrime;
        }
    }
};
void add_signed128(Fnv& hash, i128 value) {
    const u128 bits = static_cast<u128>(value);
    hash.add(static_cast<u64>(bits));
    hash.add(static_cast<u64>(bits >> 64));
}
std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    u128 magnitude = static_cast<u128>(value);
    if (negative) magnitude = u128{0} - magnitude;
    std::string result;
    while (magnitude) {
        result.push_back(static_cast<char>('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}
std::string hex8(u32 value) {
    std::ostringstream result;
    result << std::hex << std::setw(8) << std::setfill('0') << value;
    return result.str();
}
u32 parse_hex(const std::string& token) {
    const u64 wide = std::stoull(token, nullptr, 16);
    check(wide < (UINT64_C(1) << 30), "body escaped label universe");
    return static_cast<u32>(wide);
}
i64 checked_lcm(i64 a, i64 b) {
    const i128 wide = static_cast<i128>(a / std::gcd(a, b)) * b;
    check(wide <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(wide);
}
bool safe_on_cell(int speed, i64 grid, i64 left, i64 right) {
    i128 phase2 = static_cast<i128>(speed) * (left + static_cast<i128>(right));
    phase2 %= static_cast<i128>(2) * grid;
    if (phase2 < 0) phase2 += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * phase2 >= grid &&
           static_cast<i128>(7) * phase2 <= static_cast<i128>(13) * grid;
}
struct Geometry {
    i64 grid = 0;
    u64 cells = 0, pair_safe_cells = 0;
    i64 pair_ticks = 0;
    std::vector<std::pair<u32, i64>> all_classes;
    std::vector<std::pair<u32, i64>> low_classes;
    u64 all_fnv = 0, low_fnv = 0;
};
Geometry geometry_for(int q) {
    Geometry geometry;
    geometry.grid = 1;
    for (int speed : kPool)
        geometry.grid = checked_lcm(geometry.grid, 14LL * speed);
    geometry.grid = checked_lcm(geometry.grid, 14LL * q);
    geometry.grid = checked_lcm(geometry.grid, 14LL * 588);
    std::vector<i64> walls{0, geometry.grid};
    const auto add_walls = [&](int speed) {
        check(geometry.grid % (14LL * speed) == 0, "nonintegral wall unit");
        const i64 unit = geometry.grid / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(q); add_walls(588);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    geometry.cells = walls.size() - 1;
    std::map<u32, i64> classes;
    for (std::size_t i = 1; i < walls.size(); ++i) {
        const i64 left = walls[i - 1], right = walls[i];
        if (!safe_on_cell(q, geometry.grid, left, right) ||
            !safe_on_cell(588, geometry.grid, left, right)) continue;
        ++geometry.pair_safe_cells;
        const i64 width = right - left;
        geometry.pair_ticks += width;
        u32 failure = 0;
        for (unsigned bit = 0; bit < kPool.size(); ++bit)
            if (!safe_on_cell(kPool[bit], geometry.grid, left, right))
                failure |= UINT32_C(1) << bit;
        classes[failure] += width;
    }
    Fnv all_hash, low_hash;
    for (const auto& item : classes) {
        geometry.all_classes.push_back(item);
        all_hash.add(item.first); all_hash.add(item.second);
        if (std::popcount(item.first) <= 9) {
            geometry.low_classes.push_back(item);
            low_hash.add(item.first); low_hash.add(item.second);
        }
    }
    geometry.all_fnv = all_hash.state;
    geometry.low_fnv = low_hash.state;
    return geometry;
}
struct FailureLedger {
    std::map<int, std::vector<u32>> by_q;
    std::map<int, u64> row_fnv;
    u64 count = 0, global_fnv = 0;
};
FailureLedger read_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    check(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    check(std::getline(input, line) && line == "q,r,body_hex",
          "failure header changed");
    FailureLedger result;
    std::map<int, Fnv> row_hashes;
    std::set<std::pair<int, u32>> distinct;
    Fnv global;
    int prior_q = -1; u32 prior_body = 0;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto first = line.find(','), second = line.find(',', first + 1);
        check(first != std::string::npos && second != std::string::npos,
              "malformed failure row");
        const int q = std::stoi(line.substr(0, first));
        const int r = std::stoi(line.substr(first + 1, second - first - 1));
        const u32 body = parse_hex(line.substr(second + 1));
        check(r == 588 && std::popcount(body) == 9 &&
                  distinct.insert({q, body}).second,
              "failure row endpoint/rank/distinctness changed");
        check(q > prior_q || (q == prior_q && body > prior_body),
              "failure row order changed");
        prior_q = q; prior_body = body;
        result.by_q[q].push_back(body);
        row_hashes[q].add(body);
        global.add(q); global.add(r); global.add(body);
        ++result.count;
    }
    check(input.eof() && result.count == 144867 && result.by_q.size() == 10,
          "failure census changed");
    for (const auto& [q, hash] : row_hashes) result.row_fnv[q] = hash.state;
    result.global_fnv = global.state;
    check(result.global_fnv == UINT64_C(0x6f51fa88f3b09cdc),
          "failure ledger FNV changed");
    return result;
}
struct RowSummary {
    i128 low_min = 0, low_max = 0, full_min = 0, full_max = 0;
    u32 low_min_body = 0, low_max_body = 0;
    u32 full_min_body = 0, full_max_body = 0;
    u64 low_positive = 0, full_positive = 0, equality = 0, detail_fnv = 0;
};
}

int main(int argc, char** argv) {
    try {
        check(argc == 3, "usage: literal_audit FAILURES DETAIL_OUT");
        const FailureLedger failures = read_failures(argv[1]);
        std::ofstream detail(argv[2], std::ios::binary);
        check(static_cast<bool>(detail), "cannot create detail ledger");
        detail << "q,r,ordinal,body_hex,truncated_mass,full_mass,"
                  "truncated_scaled_ticks,full_scaled_ticks\n";
        std::map<int, Geometry> geometries;
        std::map<int, RowSummary> summaries;
        for (const auto& [q, bodies] : failures.by_q) {
            const Geometry geometry = geometry_for(q);
            RowSummary summary;
            Fnv detail_hash;
            for (std::size_t ordinal = 0; ordinal < bodies.size(); ++ordinal) {
                const u32 body = bodies[ordinal];
                i64 low_mass = 0, full_mass = 0;
                for (const auto& [failure, width] : geometry.low_classes)
                    if ((failure & body) == 0) low_mass += width;
                for (const auto& [failure, width] : geometry.all_classes)
                    if ((failure & body) == 0) full_mass += width;
                check(low_mass <= full_mass, "truncated mass exceeded full mass");
                const i128 low_ticks = static_cast<i128>(63) * low_mass -
                                       static_cast<i128>(4) * geometry.grid;
                const i128 full_ticks = static_cast<i128>(63) * full_mass -
                                        static_cast<i128>(4) * geometry.grid;
                detail << q << ",588," << ordinal << ',' << hex8(body) << ','
                       << low_mass << ',' << full_mass << ',' << decimal(low_ticks)
                       << ',' << decimal(full_ticks) << '\n';
                detail_hash.add(q); detail_hash.add(588); detail_hash.add(ordinal);
                detail_hash.add(body); detail_hash.add(low_mass); detail_hash.add(full_mass);
                add_signed128(detail_hash, low_ticks);
                add_signed128(detail_hash, full_ticks);
                const bool first = ordinal == 0;
                if (first || low_ticks < summary.low_min ||
                    (low_ticks == summary.low_min && body < summary.low_min_body)) {
                    summary.low_min = low_ticks; summary.low_min_body = body;
                }
                if (first || low_ticks > summary.low_max ||
                    (low_ticks == summary.low_max && body < summary.low_max_body)) {
                    summary.low_max = low_ticks; summary.low_max_body = body;
                }
                if (first || full_ticks < summary.full_min ||
                    (full_ticks == summary.full_min && body < summary.full_min_body)) {
                    summary.full_min = full_ticks; summary.full_min_body = body;
                }
                if (first || full_ticks > summary.full_max ||
                    (full_ticks == summary.full_max && body < summary.full_max_body)) {
                    summary.full_max = full_ticks; summary.full_max_body = body;
                }
                summary.low_positive += low_ticks > 0;
                summary.full_positive += full_ticks > 0;
                summary.equality += low_mass == full_mass;
            }
            summary.detail_fnv = detail_hash.state;
            geometries[q] = geometry;
            summaries[q] = summary;
        }
        check(detail.good(), "detail write failed");
        u64 total_low = 0, total_full = 0, total_equality = 0;
        Fnv summary_hash;
        std::cout << "LRC14_ENDPOINT588_CLEANROOM_LITERAL_AUDIT_V1\n"
                  << "FAILURE_ROWS " << failures.by_q.size() << " BODIES "
                  << failures.count << " FAILURE_FNV " << std::hex
                  << failures.global_fnv << std::dec << '\n';
        for (const auto& [q, bodies] : failures.by_q) {
            const auto& geometry = geometries.at(q);
            const auto& summary = summaries.at(q);
            total_low += summary.low_positive; total_full += summary.full_positive;
            total_equality += summary.equality;
            summary_hash.add(q); summary_hash.add(bodies.size());
            summary_hash.add(geometry.grid); summary_hash.add(geometry.cells);
            summary_hash.add(geometry.pair_safe_cells);
            summary_hash.add(geometry.pair_ticks);
            summary_hash.add(geometry.all_classes.size());
            summary_hash.add(geometry.low_classes.size());
            summary_hash.add(geometry.all_fnv); summary_hash.add(geometry.low_fnv);
            add_signed128(summary_hash, summary.low_min);
            summary_hash.add(summary.low_min_body);
            add_signed128(summary_hash, summary.full_min);
            summary_hash.add(summary.full_min_body);
            summary_hash.add(summary.low_positive); summary_hash.add(summary.full_positive);
            summary_hash.add(summary.equality); summary_hash.add(summary.detail_fnv);
            std::cout << "ROW " << q << ",588 BODIES " << bodies.size()
                      << " BODY_FNV " << std::hex << failures.row_fnv.at(q)
                      << std::dec << " GRID " << geometry.grid << " CELLS "
                      << geometry.cells << " PAIR_SAFE_CELLS "
                      << geometry.pair_safe_cells << " PAIR_TICKS "
                      << geometry.pair_ticks << " ALL_CLASSES "
                      << geometry.all_classes.size() << " LOW_CLASSES "
                      << geometry.low_classes.size() << " ALL_CLASS_FNV "
                      << std::hex << geometry.all_fnv << " LOW_CLASS_FNV "
                      << geometry.low_fnv << std::dec << " LOW_POSITIVE "
                      << summary.low_positive << " FULL_POSITIVE "
                      << summary.full_positive << " EQUALITY " << summary.equality
                      << " LOW_MIN_TICKS " << decimal(summary.low_min)
                      << " LOW_MIN_BODY " << hex8(summary.low_min_body)
                      << " FULL_MIN_TICKS " << decimal(summary.full_min)
                      << " FULL_MIN_BODY " << hex8(summary.full_min_body)
                      << " LOW_MAX_TICKS " << decimal(summary.low_max)
                      << " LOW_MAX_BODY " << hex8(summary.low_max_body)
                      << " FULL_MAX_TICKS " << decimal(summary.full_max)
                      << " FULL_MAX_BODY " << hex8(summary.full_max_body)
                      << " DETAIL_FNV " << std::hex << summary.detail_fnv
                      << std::dec << '\n';
        }
        std::cout << "SUMMARY BODIES " << failures.count << " LOW_POSITIVE "
                  << total_low << " FULL_POSITIVE " << total_full << " EQUALITY "
                  << total_equality << " SUMMARY_FNV " << std::hex
                  << summary_hash.state << std::dec << '\n'
                  << "INEQUALITY TRUNCATED_RANK_LE9_MASS_LE_FULL_LITERAL_MASS\n"
                  << "SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT588_FAILURE_BODIES_"
                     "NO_CARRIER_EXCHANGE_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT "
                  << (total_low == failures.count && total_full == failures.count
                          ? "PASS" : "HOSTILE_FAIL") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CLEANROOM_LITERAL_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
