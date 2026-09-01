// Independent finite-exact lower-bound audit for the residual bodies of the
// endpoint-589 exchanged-carrier scan.  Only pair-safe cells whose pool
// failure mask has rank at most nine are counted.  Consequently the recorded
// mass is a lower bound for the literal safe-set mass of each nine-body.

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
constexpr int kEndpoint = 589;
constexpr u64 kFailureLedgerFnv = UINT64_C(0x2228c33b4ef3b01d);
constexpr u64 kQ50BodyFnv = UINT64_C(0xff421454f02d9099);
constexpr u64 kQ96BodyFnv = UINT64_C(0x6f70a2d28ff6a957);
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
    i64 pair_ticks = 0;
    std::map<u32, i64> low_classes;
    u64 class_fnv = 0;
};

Geometry build_geometry(int q) {
    Geometry geometry;
    geometry.grid = checked_lcm(base_grid(), 14LL * q);
    geometry.grid = checked_lcm(geometry.grid, 14LL * kEndpoint);
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
    add_walls(kEndpoint);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    geometry.cells = walls.size() - 1;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(q, geometry.grid, left, right) ||
            !safe_midpoint(kEndpoint, geometry.grid, left, right))
            continue;
        const i64 width = right - left;
        geometry.pair_ticks += width;
        u32 failure = 0;
        for (unsigned bit = 0; bit < kPool.size(); ++bit)
            if (!safe_midpoint(kPool[bit], geometry.grid, left, right))
                failure |= UINT32_C(1) << bit;
        if (std::popcount(failure) <= 9)
            geometry.low_classes[failure] += width;
    }
    Fnv class_ledger;
    for (const auto& [failure, width] : geometry.low_classes) {
        class_ledger.add(failure);
        class_ledger.add(width);
    }
    geometry.class_fnv = class_ledger.state;
    return geometry;
}

struct FailureRows {
    std::array<std::vector<u32>, 2> bodies;
};

FailureRows read_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "failure header changed");
    FailureRows result;
    std::array<std::set<u32>, 2> distinct;
    std::array<Fnv, 2> body_ledgers;
    Fnv global;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q = 0, r = 0;
        std::string token;
        fields >> q >> r >> token;
        require(fields && r == kEndpoint && (q == 50 || q == 96),
                "failure row escaped hostile endpoint-589 rows");
        const u64 wide = std::stoull(token, nullptr, 16);
        require(wide < (UINT64_C(1) << 30), "body escaped label universe");
        const u32 body = static_cast<u32>(wide);
        const std::size_t slot = q == 50 ? 0 : 1;
        require(std::popcount(body) == 9 && distinct[slot].insert(body).second,
                "body rank/distinctness changed");
        result.bodies[slot].push_back(body);
        body_ledgers[slot].add(body);
        global.add(q);
        global.add(r);
        global.add(body);
    }
    require(input.eof() && result.bodies[0].size() == 20025 &&
                result.bodies[1].size() == 11 &&
                body_ledgers[0].state == kQ50BodyFnv &&
                body_ledgers[1].state == kQ96BodyFnv &&
                global.state == kFailureLedgerFnv,
            "failure ledger identity changed");
    return result;
}

struct RowSummary {
    i128 minimum = 0;
    u32 minimum_body = 0;
    i128 maximum = 0;
    u32 maximum_body = 0;
    u64 positive = 0;
    u64 ledger_fnv = 0;
};
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3, "usage: primary FAILURES DETAIL_CSV");
        const FailureRows failures = read_failures(argv[1]);
        std::ofstream detail(argv[2], std::ios::binary);
        require(static_cast<bool>(detail), "cannot create detail ledger");
        detail << "q,r,ordinal,body_hex,truncated_mass,scaled_ticks\n";
        std::array<Geometry, 2> geometries{
            build_geometry(50), build_geometry(96)};
        std::array<RowSummary, 2> summaries{};
        for (std::size_t slot = 0; slot < 2; ++slot) {
            const int q = slot == 0 ? 50 : 96;
            const Geometry& geometry = geometries[slot];
            Fnv ledger;
            for (std::size_t ordinal = 0;
                 ordinal < failures.bodies[slot].size(); ++ordinal) {
                const u32 body = failures.bodies[slot][ordinal];
                i64 mass = 0;
                for (const auto& [failure, width] : geometry.low_classes)
                    if ((failure & body) == 0) mass += width;
                const i128 ticks = static_cast<i128>(63) * mass -
                                   static_cast<i128>(4) * geometry.grid;
                detail << q << ',' << kEndpoint << ',' << ordinal << ','
                       << hex8(body) << ',' << mass << ',' << decimal(ticks)
                       << '\n';
                ledger.add(q);
                ledger.add(kEndpoint);
                ledger.add(ordinal);
                ledger.add(body);
                ledger.add(mass);
                add_i128(ledger, ticks);
                RowSummary& summary = summaries[slot];
                if (ordinal == 0 || ticks < summary.minimum ||
                    (ticks == summary.minimum && body < summary.minimum_body)) {
                    summary.minimum = ticks;
                    summary.minimum_body = body;
                }
                if (ordinal == 0 || ticks > summary.maximum ||
                    (ticks == summary.maximum && body < summary.maximum_body)) {
                    summary.maximum = ticks;
                    summary.maximum_body = body;
                }
                summary.positive += ticks > 0;
            }
            summaries[slot].ledger_fnv = ledger.state;
        }
        require(detail.good(), "detail ledger write failed");

        require(geometries[0].grid == INT64_C(2827379709554400) &&
                    geometries[0].cells == 8389 &&
                    geometries[0].pair_ticks == INT64_C(2077256602813392) &&
                    geometries[0].low_classes.size() == 2383 &&
                    geometries[0].class_fnv ==
                        UINT64_C(0x88d3eb2d7a477232),
                "q50 geometry changed");
        require(geometries[1].grid == INT64_C(1130951883821760) &&
                    geometries[1].cells == 8501 &&
                    geometries[1].pair_ticks == INT64_C(830901383902590) &&
                    geometries[1].low_classes.size() == 2352 &&
                    geometries[1].class_fnv ==
                        UINT64_C(0xcb74e72091a68363),
                "q96 geometry changed");
        require(summaries[0].positive == failures.bodies[0].size() &&
                    summaries[0].minimum ==
                        static_cast<i128>(14566818763788984LL) &&
                    summaries[0].minimum_body == UINT32_C(0x013c6401) &&
                    summaries[0].maximum ==
                        static_cast<i128>(26685137010259728LL) &&
                    summaries[0].maximum_body == UINT32_C(0x21884443) &&
                    summaries[0].ledger_fnv ==
                        UINT64_C(0x9a0cd88b508499a2) &&
                    summaries[1].positive == failures.bodies[1].size() &&
                    summaries[1].minimum ==
                        static_cast<i128>(7172391058639758LL) &&
                    summaries[1].minimum_body == UINT32_C(0x0d0c6401) &&
                    summaries[1].maximum ==
                        static_cast<i128>(8558758749944214LL) &&
                    summaries[1].maximum_body == UINT32_C(0x35126400) &&
                    summaries[1].ledger_fnv ==
                        UINT64_C(0x313a0cdba0e5ac5c),
                "direct literal positivity/minimum changed");

        std::cout << "LRC14_ENDPOINT589_DIRECT_LITERAL_PRIMARY_V1\n";
        for (std::size_t slot = 0; slot < 2; ++slot) {
            const int q = slot == 0 ? 50 : 96;
            const Geometry& geometry = geometries[slot];
            const RowSummary& summary = summaries[slot];
            std::cout << "ROW " << q << ',' << kEndpoint << " BODIES "
                      << failures.bodies[slot].size() << " GRID "
                      << geometry.grid << " CELLS " << geometry.cells
                      << " PAIR_TICKS " << geometry.pair_ticks
                      << " LOW_CLASSES " << geometry.low_classes.size()
                      << " CLASS_FNV " << std::hex << geometry.class_fnv
                      << std::dec << " POSITIVE " << summary.positive
                      << " MIN_TICKS " << decimal(summary.minimum)
                      << " MIN_BODY " << hex8(summary.minimum_body)
                      << " MAX_TICKS " << decimal(summary.maximum)
                      << " MAX_BODY " << hex8(summary.maximum_body)
                      << " DETAIL_FNV " << std::hex << summary.ledger_fnv
                      << std::dec << '\n';
        }
        std::cout << "TOTAL_BODIES 20036 ALL_STRICTLY_POSITIVE 1\n"
                  << "INEQUALITY TRUNCATED_LOW_RANK_CELL_MASS_LE_LITERAL_SAFE_MASS\n"
                  << "SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT589_FAILURE_BODIES_"
                     "NO_CARRIER_EXCHANGE_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT589_DIRECT_LITERAL_PRIMARY_ERROR "
                  << error.what() << '\n';
        return 1;
    }
}
