// Maintained import-free exact audit of the mixed rank-8/rank-9 optimum at
// (100,629).  The response atlas is treated as claimed data: this consumer
// reconstructs literal wall geometry, all representative responses/margins,
// the exact 7/2 dual, and the stated four-mask primal certificate.

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
constexpr i64 kFixedGrid = INT64_C(18241159416480);
constexpr i64 kPairGrid = INT64_C(3374614492048800);
constexpr i64 kPairTicks = INT64_C(2479303131865560);
constexpr u64 kEmptyFnv = UINT64_C(0xcbf29ce484222325);
constexpr u64 kFailureFnv = UINT64_C(0x61ee3b61e7de8594);
constexpr u64 kClassFnv = UINT64_C(0x288fc0fd18371855);
constexpr std::array<u32, 4> kWitness = {
    UINT32_C(0x002ac4c0), UINT32_C(0x3882a082),
    UINT32_C(0x0041c325), UINT32_C(0x08c28e40)};
constexpr std::array<u32, 4> kResponse = {
    UINT32_C(0x0c00000), UINT32_C(0x0000c27),
    UINT32_C(0x0687198), UINT32_C(0xf178240)};
constexpr std::array<unsigned, 4> kRank = {8, 9, 9, 9};
constexpr std::array<i64, 4> kMass = {
    INT64_C(214760471161330), INT64_C(215227333274140),
    INT64_C(214424024475422), INT64_C(215167289205628)};
constexpr std::array<i128, 4> kTicks = {
    static_cast<i128>(INT64_C(31451714968590)),
    static_cast<i128>(INT64_C(60864028075620)),
    static_cast<i128>(INT64_C(10255573756386)),
    static_cast<i128>(INT64_C(57081251759364))};
constexpr std::array<unsigned char, 28> kDualQuarterUnits = {
    1, 1, 0, 1, 1, 2, 0, 1, 1, 2, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0};

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

i128 parse_decimal(const std::string& text) {
    require(!text.empty(), "empty integer field");
    std::size_t index = 0;
    bool negative = false;
    if (text[0] == '-') {
        negative = true;
        index = 1;
    }
    require(index < text.size(), "sign-only integer field");
    i128 value = 0;
    for (; index < text.size(); ++index) {
        require(text[index] >= '0' && text[index] <= '9',
                "nondigit integer field");
        value = 10 * value + (text[index] - '0');
    }
    return negative ? -value : value;
}

std::string hex8(u32 value) {
    std::ostringstream output;
    output << std::hex << std::setw(8) << std::setfill('0') << value;
    return output.str();
}

std::vector<std::string> split(const std::string& line, char separator) {
    std::vector<std::string> fields;
    std::stringstream input(line);
    std::string field;
    while (std::getline(input, field, separator)) fields.push_back(field);
    return fields;
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

Geometry build_geometry() {
    i64 grid = fixed_grid();
    grid = checked_lcm(grid, 14LL * 100);
    grid = checked_lcm(grid, 14LL * 629);
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
    add_walls(100);
    add_walls(629);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::map<u32, i64> by_failure;
    Geometry geometry;
    geometry.grid = grid;
    geometry.cells = walls.size() - 1;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(100, grid, left, right) ||
            !safe_midpoint(629, grid, left, right))
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
        if (std::popcount(entry.first) <= 9)
            geometry.classes.push_back(entry);
    require(geometry.grid == kPairGrid && geometry.cells == 8557 &&
                geometry.pair_ticks == kPairTicks &&
                geometry.classes.size() == 2397,
            "literal geometry identity changed");
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
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto fields = split(line, ',');
        require(fields.size() == 3 && fields[0] == "100" &&
                    fields[1] == "629",
                "malformed/non-r629 failure row");
        const u32 body = static_cast<u32>(std::stoul(fields[2], nullptr, 16));
        require(std::popcount(body) == 9, "failure body rank changed");
        bodies.push_back(body);
    }
    Fnv ledger;
    for (u32 body : bodies) ledger.add(body);
    require(input.eof() && bodies.size() == 28 &&
                std::is_sorted(bodies.begin(), bodies.end()) &&
                std::adjacent_find(bodies.begin(), bodies.end()) == bodies.end() &&
                ledger.state == kFailureFnv,
            "failure ledger identity changed");
    return bodies;
}

struct ClassRow {
    u32 pattern = 0;
    u64 count = 0;
    u32 mask = 0;
    unsigned rank = 0;
    u64 count8 = 0;
    u64 count9 = 0;
    i64 mass = 0;
    i128 ticks = 0;
};

std::vector<ClassRow> read_atlas(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mixed atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "a_hex\tb_hex\tcover\tcount\tleast_mask\tleast_rank"
                        "\tcount8\tcount9\tmass\tmargin_ticks63",
            "atlas header changed");
    std::vector<ClassRow> rows;
    u32 previous = 0;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto f = split(line, '\t');
        require(f.size() == 10 && std::stoull(f[1], nullptr, 16) == 0,
                "malformed atlas row");
        ClassRow row;
        row.pattern = static_cast<u32>(std::stoull(f[0], nullptr, 16));
        const unsigned stated_cover = std::stoul(f[2]);
        row.count = std::stoull(f[3]);
        row.mask = static_cast<u32>(std::stoul(f[4], nullptr, 16));
        row.rank = std::stoul(f[5]);
        row.count8 = std::stoull(f[6]);
        row.count9 = std::stoull(f[7]);
        row.mass = std::stoll(f[8]);
        row.ticks = parse_decimal(f[9]);
        require(row.pattern != 0 && row.pattern < (u32{1} << 28) &&
                    (rows.empty() || previous < row.pattern) &&
                    stated_cover == std::popcount(row.pattern) &&
                    row.count == row.count8 + row.count9 && row.count > 0 &&
                    (row.rank == 8 || row.rank == 9) &&
                    std::popcount(row.mask) == row.rank &&
                    ((row.count8 > 0 && row.rank == 8) ||
                     (row.count8 == 0 && row.count9 > 0 && row.rank == 9)),
                "atlas row invariant changed");
        previous = row.pattern;
        rows.push_back(row);
    }
    require(input.eof() && rows.size() == 419,
            "atlas class count changed");
    return rows;
}

std::array<std::array<u64, 10>, 31> choose{};

void init_choose() {
    for (unsigned n = 0; n <= 30; ++n) {
        choose[n][0] = 1;
        for (unsigned k = 1; k <= 9; ++k)
            choose[n][k] = k > n ? 0 :
                (k == n ? 1 : choose[n - 1][k - 1] + choose[n - 1][k]);
    }
    require(choose[30][8] == 5852925 && choose[30][9] == 14307150,
            "combination universe changed");
}

u64 colex_rank(u32 mask) {
    const unsigned rank = std::popcount(mask);
    require(rank == 8 || rank == 9, "colex rank outside mixed universe");
    u64 answer = 0;
    unsigned ordinal = 1;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((mask >> bit) & 1u) answer += choose[bit][ordinal++];
    require(ordinal == rank + 1 && answer < choose[30][rank],
            "colex escaped universe");
    return answer;
}

u32 response(const std::vector<u32>& bodies, u32 mask) {
    u32 answer = 0;
    for (unsigned index = 0; index < bodies.size(); ++index)
        if ((mask & bodies[index]) == 0) answer |= u32{1} << index;
    return answer;
}

} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3, "usage: detached-r629 FAILURE_CSV ATLAS_TSV");
        init_choose();
        const std::vector<u32> bodies = read_failures(argv[1]);
        const std::vector<ClassRow> rows = read_atlas(argv[2]);
        const Geometry geometry = build_geometry();

        const unsigned dual_units = std::accumulate(
            kDualQuarterUnits.begin(), kDualQuarterUnits.end(), 0u);
        require(dual_units == 14, "dual total is not 7/2");
        std::array<u64, 5> load_histogram{};
        unsigned maximum_load = 0;
        u32 full = 0;
        Fnv class_ledger;
        Fnv dual_ledger;
        std::map<u32, ClassRow> by_mask;
        for (const ClassRow& row : rows) {
            const u32 direct_response = response(bodies, row.mask);
            const Margin direct_margin = margin(geometry, row.mask);
            require(direct_response == row.pattern &&
                        direct_margin.mass == row.mass &&
                        direct_margin.ticks == row.ticks && row.ticks >= 0,
                    "representative response/activity changed");
            unsigned load = 0;
            for (unsigned bit = 0; bit < 28; ++bit)
                if ((row.pattern >> bit) & 1u)
                    load += kDualQuarterUnits[bit];
            require(load <= 4, "exact dual constraint violated");
            ++load_histogram[load];
            maximum_load = std::max(maximum_load, load);
            full |= row.pattern;
            require(by_mask.emplace(row.mask, row).second,
                    "repeated representative mask");
            class_ledger.add(row.pattern);
            class_ledger.add(row.count8);
            class_ledger.add(row.count9);
            class_ledger.add(row.mask);
            class_ledger.add(row.rank);
            class_ledger.add(static_cast<u64>(row.mass));
            add_i128(class_ledger, row.ticks);
            dual_ledger.add(row.pattern);
            dual_ledger.add(load);
        }
        require(full == (u32{1} << 28) - 1 && maximum_load == 4 &&
                    class_ledger.state == kClassFnv,
                "class/dual universe identity changed");
        for (unsigned bit = 0; bit < 28; ++bit) {
            dual_ledger.add(bit);
            dual_ledger.add(kDualQuarterUnits[bit]);
        }

        u32 joined = 0;
        std::array<unsigned, 28> multiplicity{};
        Fnv witness_ledger;
        Fnv mask_only_ledger;
        Fnv response_ledger;
        std::cout << "R629_MIXED_OPTIMUM_DETACHED_AUDIT_V1\n"
                  << "FAILURES 28 FNV " << std::hex << kFailureFnv
                  << " CLASSES 419 CLASS_FNV " << kClassFnv << std::dec
                  << '\n'
                  << "GEOMETRY GRID " << geometry.grid << " CELLS "
                  << geometry.cells << " PAIR_TICKS " << geometry.pair_ticks
                  << " FAILURE_CLASSES_RANK_LE9 " << geometry.classes.size()
                  << '\n'
                  << "DUAL UNITS_DENOMINATOR 4 TOTAL_UNITS " << dual_units
                  << " TOTAL 7/2 MAXIMUM_LOAD_UNITS " << maximum_load
                  << " LOAD_HISTOGRAM";
        for (unsigned load = 0; load <= 4; ++load)
            std::cout << ' ' << load << ':' << load_histogram[load];
        std::cout << " FNV " << std::hex << dual_ledger.state << std::dec
                  << '\n';
        for (unsigned ordinal = 0; ordinal < kWitness.size(); ++ordinal) {
            const u32 mask = kWitness[ordinal];
            const auto found = by_mask.find(mask);
            require(found != by_mask.end(), "witness absent from atlas");
            const ClassRow& row = found->second;
            const Margin value = margin(geometry, mask);
            const u32 direct_response = response(bodies, mask);
            const u64 colex = colex_rank(mask);
            require(row.rank == kRank[ordinal] &&
                        direct_response == kResponse[ordinal] &&
                        row.pattern == kResponse[ordinal] &&
                        value.mass == kMass[ordinal] &&
                        value.ticks == kTicks[ordinal] && value.ticks > 0,
                    "witness rank/response/activity changed");
            joined |= direct_response;
            for (unsigned bit = 0; bit < 28; ++bit)
                multiplicity[bit] += (direct_response >> bit) & 1u;
            witness_ledger.add(mask);
            witness_ledger.add(row.rank);
            witness_ledger.add(colex);
            witness_ledger.add(direct_response);
            witness_ledger.add(static_cast<u64>(value.mass));
            add_i128(witness_ledger, value.ticks);
            mask_only_ledger.add(mask);
            response_ledger.add(mask);
            response_ledger.add(row.rank);
            response_ledger.add(direct_response);
            response_ledger.add(static_cast<u64>(value.mass));
            add_i128(response_ledger, value.ticks);
            std::cout << "WITNESS " << ordinal << " MASK " << hex8(mask)
                      << " RANK " << row.rank << " COLEX " << colex
                      << " RESPONSE " << std::hex << std::setw(7)
                      << std::setfill('0') << direct_response << std::dec
                      << std::setfill(' ') << " COVER "
                      << std::popcount(direct_response) << " MASS "
                      << value.mass << " MARGIN_TICKS63 "
                      << decimal(value.ticks) << '\n';
        }
        require(joined == (u32{1} << 28) - 1,
                "four witnesses do not cover the universe");
        const u64 once = std::count(multiplicity.begin(), multiplicity.end(), 1u);
        const u64 twice = std::count(multiplicity.begin(), multiplicity.end(), 2u);
        require(once == 27 && twice == 1,
                "witness multiplicity changed");
        require(mask_only_ledger.state == UINT64_C(0x5601ef0ed7c3ecaa) &&
                    response_ledger.state == UINT64_C(0xe83f7f763d141b63),
                "witness/response identity changed");
        std::cout << "PRIMAL SIZE 4 MULTIPLICITY_ONE " << once
                  << " MULTIPLICITY_TWO " << twice << " FNV " << std::hex
                  << witness_ledger.state << " MASK_FNV "
                  << mask_only_ledger.state << " RESPONSE_FNV "
                  << response_ledger.state << std::dec << '\n'
                  << "OPTIMUM LOWER_CEIL_7_OVER_2 4 UPPER_WITNESS 4 EXACT 4\n"
                  << "SCOPE FIXED_PAIR_LITERAL_MIXED_RANK8_OR_RANK9_"
                     "RESPONSE_QUOTIENT_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_R629_MIXED_MINIMUM_FOUR\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "DETACHED_R629_ERROR " << error.what() << '\n';
        return 1;
    }
}
