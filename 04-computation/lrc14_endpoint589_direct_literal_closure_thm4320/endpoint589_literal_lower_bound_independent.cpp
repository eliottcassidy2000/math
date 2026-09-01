// Scratch-only, self-contained literal-wall audit for the two hostile rows
// left by the endpoint-589 exchanged carrier.  This source intentionally does
// not include or read the pre-existing q50 geometry dump.

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

constexpr std::array<int, 30> kPool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr i64 kFixedGrid = INT64_C(18241159416480);
constexpr u32 kPoolMask = (UINT32_C(1) << 30) - 1;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
    if (negative) value = -value;
    std::string result;
    while (value != 0) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

std::string hex8(u32 value) {
    std::ostringstream output;
    output << std::hex << std::setfill('0') << std::setw(8) << value;
    return output.str();
}

u32 parse_hex(const std::string& token) {
    std::size_t used = 0;
    const unsigned long value = std::stoul(token, &used, 16);
    require(used == token.size() && value <= kPoolMask,
            "malformed/out-of-pool body mask");
    return static_cast<u32>(value);
}

struct Fnv {
    u64 state = UINT64_C(0xcbf29ce484222325);
    void add(u64 word) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (word >> (8 * byte)) & UINT64_C(0xff);
            state *= UINT64_C(0x100000001b3);
        }
    }
};

i64 checked_lcm(i64 left, i64 right) {
    require(left > 0 && right > 0, "nonpositive LCM input");
    const i128 value = static_cast<i128>(left / std::gcd(left, right)) * right;
    require(value <= std::numeric_limits<i64>::max(), "grid overflow");
    return static_cast<i64>(value);
}

bool safe_midpoint(int speed, i64 grid, i64 left, i64 right) {
    // (left+right)/(2*grid) is an interior point of one literal wall cell.
    i128 residue = static_cast<i128>(speed) *
                   (static_cast<i128>(left) + right);
    residue %= static_cast<i128>(2) * grid;
    if (residue < 0) residue += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * residue >= grid &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * grid;
}

struct Cell {
    u32 failure = 0;
    i64 width = 0;
};

struct Geometry {
    int q = 0;
    int r = 0;
    i64 grid = 0;
    u64 wall_intervals = 0;
    u64 pair_safe_cells = 0;
    i64 pair_ticks = 0;
    std::vector<Cell> cells;
    std::map<u32, i64> classes;
    std::map<u32, i64> low_classes;
    Fnv class_ledger;
    Fnv low_class_ledger;
};

Geometry build_geometry(int q, int r) {
    i64 grid = 1;
    for (int speed : kPool) grid = checked_lcm(grid, 14LL * speed);
    require(grid == kFixedGrid, "fixed grid changed");
    grid = checked_lcm(grid, 14LL * q);
    grid = checked_lcm(grid, 14LL * r);

    std::vector<i64> walls = {0, grid};
    const auto add_walls = [&](int speed) {
        const i64 divisor = 14LL * speed;
        require(grid % divisor == 0, "nonintegral wall unit");
        const i64 unit = grid / divisor;
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(q);
    add_walls(r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.front() == 0 && walls.back() == grid,
            "literal wall boundary lost");

    Geometry geometry;
    geometry.q = q;
    geometry.r = r;
    geometry.grid = grid;
    geometry.wall_intervals = walls.size() - 1;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        require(left < right, "nonpositive wall cell");
        if (!safe_midpoint(q, grid, left, right) ||
            !safe_midpoint(r, grid, left, right)) {
            continue;
        }
        u32 failure = 0;
        for (unsigned vertex = 0; vertex < kPool.size(); ++vertex) {
            if (!safe_midpoint(kPool[vertex], grid, left, right))
                failure |= UINT32_C(1) << vertex;
        }
        const i64 width = right - left;
        ++geometry.pair_safe_cells;
        geometry.pair_ticks += width;
        geometry.cells.push_back({failure, width});
        geometry.classes[failure] += width;
        if (std::popcount(failure) <= 9)
            geometry.low_classes[failure] += width;
    }
    i64 class_sum = 0;
    for (const auto& [failure, width] : geometry.classes) {
        require(width > 0, "nonpositive aggregated class");
        class_sum += width;
        geometry.class_ledger.add(failure);
        geometry.class_ledger.add(static_cast<u64>(width));
    }
    i64 low_sum = 0;
    for (const auto& [failure, width] : geometry.low_classes) {
        require(std::popcount(failure) <= 9 && width > 0,
                "bad low-rank class");
        low_sum += width;
        geometry.low_class_ledger.add(failure);
        geometry.low_class_ledger.add(static_cast<u64>(width));
    }
    require(class_sum == geometry.pair_ticks && low_sum <= class_sum,
            "class aggregation changed mass");
    return geometry;
}

struct FailureRows {
    std::vector<u32> q50;
    std::vector<u32> q96;
    Fnv total_ledger;
    Fnv q50_ledger;
    Fnv q96_ledger;
};

FailureRows read_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open failure ledger");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "failure header changed");
    FailureRows rows;
    std::set<std::pair<int, u32>> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q = 0, r = 0;
        std::string token, extra;
        require(static_cast<bool>(fields >> q >> r >> token) &&
                    !(fields >> extra),
                "malformed failure row");
        require((q == 50 || q == 96) && r == 589,
                "unexpected failed pair");
        const u32 body = parse_hex(token);
        require(std::popcount(body) == 9 && distinct.emplace(q, body).second,
                "failure rank/distinctness changed");
        rows.total_ledger.add(q);
        rows.total_ledger.add(r);
        rows.total_ledger.add(body);
        if (q == 50) {
            rows.q50.push_back(body);
            rows.q50_ledger.add(body);
        } else {
            rows.q96.push_back(body);
            rows.q96_ledger.add(body);
        }
    }
    require(input.eof() && rows.q50.size() == 20025 &&
                rows.q96.size() == 11 &&
                rows.total_ledger.state == UINT64_C(0x2228c33b4ef3b01d) &&
                rows.q50_ledger.state == UINT64_C(0xff421454f02d9099) &&
                rows.q96_ledger.state == UINT64_C(0x6f70a2d28ff6a957),
            "failure identity changed");
    return rows;
}

struct BodyMass {
    i64 low = 0;
    i64 full = 0;
};

BodyMass mass_from_raw_cells(const Geometry& geometry, u32 body) {
    BodyMass result;
    for (const Cell& cell : geometry.cells) {
        if ((cell.failure & body) != 0) continue;
        result.full += cell.width;
        if (std::popcount(cell.failure) <= 9) result.low += cell.width;
    }
    return result;
}

BodyMass mass_from_classes(const Geometry& geometry, u32 body) {
    BodyMass result;
    for (const auto& [failure, width] : geometry.classes)
        if ((failure & body) == 0) result.full += width;
    for (const auto& [failure, width] : geometry.low_classes)
        if ((failure & body) == 0) result.low += width;
    return result;
}

struct RowAudit {
    i128 minimum_low_surplus = 0;
    i128 minimum_full_surplus = 0;
    i64 minimum_low_mass = 0;
    i64 minimum_full_mass = 0;
    u32 minimum_low_body = 0;
    u32 minimum_full_body = 0;
    i64 minimum_omitted = 0;
    i64 maximum_omitted = 0;
    u64 positive_low = 0;
    u64 equal_full_truncated = 0;
    u64 strict_full_dominance = 0;
    Fnv ledger;
};

RowAudit audit_row(const Geometry& geometry, const std::vector<u32>& bodies,
                   std::ofstream& output) {
    require(!bodies.empty(), "empty hostile row");
    RowAudit audit;
    audit.minimum_low_surplus = std::numeric_limits<i64>::max();
    audit.minimum_full_surplus = std::numeric_limits<i64>::max();
    audit.minimum_omitted = std::numeric_limits<i64>::max();
    for (u32 body : bodies) {
        const BodyMass raw = mass_from_raw_cells(geometry, body);
        const BodyMass aggregate = mass_from_classes(geometry, body);
        require(raw.low == aggregate.low && raw.full == aggregate.full,
                "raw-cell/class path mismatch");
        require(raw.low <= raw.full, "truncated mass exceeds full mass");
        const i128 low_surplus = static_cast<i128>(63) * raw.low -
                                 static_cast<i128>(4) * geometry.grid;
        const i128 full_surplus = static_cast<i128>(63) * raw.full -
                                  static_cast<i128>(4) * geometry.grid;
        require(full_surplus - low_surplus ==
                    static_cast<i128>(63) * (raw.full - raw.low),
                "saturation/lower-bound identity failed");
        audit.positive_low += low_surplus > 0;
        audit.equal_full_truncated += raw.full == raw.low;
        audit.strict_full_dominance += raw.full > raw.low;
        if (low_surplus < audit.minimum_low_surplus) {
            audit.minimum_low_surplus = low_surplus;
            audit.minimum_low_mass = raw.low;
            audit.minimum_low_body = body;
        }
        if (full_surplus < audit.minimum_full_surplus) {
            audit.minimum_full_surplus = full_surplus;
            audit.minimum_full_mass = raw.full;
            audit.minimum_full_body = body;
        }
        audit.minimum_omitted =
            std::min(audit.minimum_omitted, raw.full - raw.low);
        audit.maximum_omitted =
            std::max(audit.maximum_omitted, raw.full - raw.low);
        audit.ledger.add(body);
        audit.ledger.add(static_cast<u64>(raw.low));
        audit.ledger.add(static_cast<u64>(low_surplus));
        audit.ledger.add(static_cast<u64>(raw.full));
        audit.ledger.add(static_cast<u64>(full_surplus));
        output << geometry.q << ',' << geometry.r << ',' << hex8(body) << ','
               << raw.low << ',' << decimal(low_surplus) << ',' << raw.full
               << ',' << decimal(full_surplus) << ',' << raw.full - raw.low
               << '\n';
    }
    require(audit.positive_low == bodies.size(),
            "nonpositive truncated lower bound survived");
    require(audit.equal_full_truncated + audit.strict_full_dominance ==
                bodies.size(),
            "full/truncated comparison partition failed");
    return audit;
}

void print_geometry(const Geometry& geometry) {
    std::cout << "GEOMETRY " << geometry.q << ',' << geometry.r
              << " GRID " << geometry.grid
              << " GRID_OVER_FIXED " << geometry.grid / kFixedGrid
              << " WALL_INTERVALS " << geometry.wall_intervals
              << " PAIR_SAFE_CELLS " << geometry.pair_safe_cells
              << " PAIR_TICKS " << geometry.pair_ticks
              << " CLASSES_ALL " << geometry.classes.size()
              << " CLASSES_RANK_LE_9 " << geometry.low_classes.size()
              << " CLASS_FNV " << std::hex << geometry.class_ledger.state
              << " LOW_CLASS_FNV " << geometry.low_class_ledger.state
              << std::dec << '\n';
}

void print_row(const Geometry& geometry, const std::vector<u32>& bodies,
               const RowAudit& audit) {
    std::cout << "ROW " << geometry.q << ',' << geometry.r
              << " BODIES " << bodies.size()
              << " POSITIVE_TRUNCATED " << audit.positive_low
              << " MIN_TRUNCATED_MASS " << audit.minimum_low_mass
              << " MIN_TRUNCATED_SURPLUS "
              << decimal(audit.minimum_low_surplus)
              << " MIN_TRUNCATED_BODY " << hex8(audit.minimum_low_body)
              << " MIN_FULL_MASS " << audit.minimum_full_mass
              << " MIN_FULL_SURPLUS " << decimal(audit.minimum_full_surplus)
              << " MIN_FULL_BODY " << hex8(audit.minimum_full_body)
              << " MIN_OMITTED_TICKS " << audit.minimum_omitted
              << " MAX_OMITTED_TICKS " << audit.maximum_omitted
              << " FULL_EQUALS_TRUNCATED " << audit.equal_full_truncated
              << " FULL_STRICTLY_DOMINATES " << audit.strict_full_dominance
              << " LEDGER_FNV " << std::hex << audit.ledger.state << std::dec
              << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3,
                "usage: endpoint589_literal_lower_bound_audit "
                "FAILURES_CSV OUTPUT_LEDGER");
        const FailureRows rows = read_failures(argv[1]);
        const Geometry q50 = build_geometry(50, 589);
        const Geometry q96 = build_geometry(96, 589);
        require(q50.grid == INT64_C(2827379709554400) &&
                    q50.wall_intervals == 8389 &&
                    q50.pair_safe_cells == 5798 &&
                    q50.pair_ticks == INT64_C(2077256602813392) &&
                    q50.classes.size() == 2434 &&
                    q50.low_classes.size() == 2383 &&
                    q50.class_ledger.state == UINT64_C(0x3dfcb37c5aa053f9) &&
                    q50.low_class_ledger.state ==
                        UINT64_C(0x88d3eb2d7a477232),
                "q50 geometry identity changed");
        require(q96.grid == INT64_C(1130951883821760) &&
                    q96.wall_intervals == 8501 &&
                    q96.pair_safe_cells == 5780 &&
                    q96.pair_ticks == INT64_C(830901383902590) &&
                    q96.classes.size() == 2423 &&
                    q96.low_classes.size() == 2352 &&
                    q96.class_ledger.state == UINT64_C(0x3f384eb372969245) &&
                    q96.low_class_ledger.state ==
                        UINT64_C(0xcb74e72091a68363),
                "q96 geometry identity changed");

        std::ofstream output(argv[2]);
        require(static_cast<bool>(output), "cannot create body-mass ledger");
        output << "q,r,body_hex,truncated_mass_ticks,truncated_surplus,"
                  "full_mass_ticks,full_surplus,omitted_high_rank_ticks\n";
        const RowAudit audit50 = audit_row(q50, rows.q50, output);
        const RowAudit audit96 = audit_row(q96, rows.q96, output);
        require(output.good(), "body-mass ledger write failed");

        require(audit50.minimum_low_surplus ==
                    static_cast<i128>(INT64_C(14566818763788984)) &&
                    audit50.minimum_low_body == UINT32_C(0x013c6401) &&
                    audit50.minimum_low_mass == INT64_C(410735517492168) &&
                    audit50.minimum_full_surplus == audit50.minimum_low_surplus &&
                    audit50.minimum_full_body == audit50.minimum_low_body &&
                    audit50.minimum_full_mass == audit50.minimum_low_mass &&
                    audit50.equal_full_truncated == 16788 &&
                    audit50.strict_full_dominance == 3237 &&
                    audit50.maximum_omitted == INT64_C(8146986595540) &&
                    audit50.ledger.state == UINT64_C(0x1ae6cd55bb7f7c32),
                "q50 minimum changed");
        require(audit96.minimum_low_surplus ==
                    static_cast<i128>(INT64_C(7172391058639758)) &&
                    audit96.minimum_low_body == UINT32_C(0x0d0c6401) &&
                    audit96.minimum_low_mass == INT64_C(185653945935346) &&
                    audit96.minimum_full_surplus == audit96.minimum_low_surplus &&
                    audit96.minimum_full_body == audit96.minimum_low_body &&
                    audit96.minimum_full_mass == audit96.minimum_low_mass &&
                    audit96.equal_full_truncated == 4 &&
                    audit96.strict_full_dominance == 7 &&
                    audit96.maximum_omitted == INT64_C(544959807964) &&
                    audit96.ledger.state == UINT64_C(0x8c0f0a16ba384407),
                "q96 minimum changed");

        std::cout << "LRC14_ENDPOINT589_LITERAL_LOWER_BOUND_AUDIT_V1\n";
        print_geometry(q50);
        print_geometry(q96);
        print_row(q50, rows.q50, audit50);
        print_row(q96, rows.q96, audit96);
        std::cout << "NORMALIZATION TRUNCATED_SURPLUS=63*TRUNCATED_MASS-4*GRID "
                     "FULL_MINUS_TRUNCATED=63*OMITTED_HIGH_RANK_MASS\n"
                  << "QUANTIFIER ALL_20025_Q50_AND_ALL_11_Q96_"
                     "EXACT_CARRIER_FAILURE_BODIES\n"
                  << "SCOPE FINITE_EXACT_LITERAL_TARGET_BODY_LOWER_BOUND_"
                     "FIXED_POOL_ENDPOINT589_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT589_LITERAL_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
