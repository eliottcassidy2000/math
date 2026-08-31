// Detached literal-wall and exhaustive-body audit of the endpoint-636
// signature-{107,173,247} fibre surgery.  This translation unit imports no
// project source.  It independently reads the frozen joint deck, exact
// 21-row fibre, and four-mask witness; rebuilds the deck; exhausts every
// labelled rank-nine body; and reconstructs literal walls for every pair.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
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
constexpr std::size_t kJointCount = 421;
constexpr u64 kJointFnv = UINT64_C(0x20d63dd42fe8150e);
constexpr std::size_t kFibreCount = 21;
constexpr u64 kFibreFnv = UINT64_C(0xeadefa2fae582ca7);
constexpr std::size_t kWitnessCount = 4;
constexpr u64 kWitnessFnv = UINT64_C(0x1796af626bcea63b);
constexpr std::array<std::size_t, 3> kDeletedIndices{107, 173, 247};
constexpr std::array<u32, 3> kDeletedMasks{
    UINT32_C(0x128c8900), UINT32_C(0x0280b206), UINT32_C(0x2380c108)};
constexpr std::array<u32, 4> kExpectedWitness{
    UINT32_C(0x12848902), UINT32_C(0x0380a520),
    UINT32_C(0x0380e012), UINT32_C(0x1288b200)};
constexpr std::size_t kRebuiltCount = 422;
constexpr u64 kRebuiltFnv = UINT64_C(0x1de1e54f8262ac49);
constexpr u64 kBodyCount = UINT64_C(14307150);
constexpr u64 kBodyChecks = UINT64_C(405392382);
constexpr u64 kBodyRowFnv = UINT64_C(0x15a176ecbab025c3);
constexpr i64 kFixedGrid = INT64_C(18241159416480);
constexpr u64 kEmptyFnv = UINT64_C(0xcbf29ce484222325);
constexpr u64 kTotalCells = UINT64_C(178475);
constexpr u64 kTotalFailureClasses = UINT64_C(52275);
constexpr u64 kActivityAggregateFnv = UINT64_C(0x75c1967679e26e88);
constexpr int kWeakestQ = 294;
constexpr int kWeakestR = 526;
constexpr u32 kWeakestMask = UINT32_C(0x01e08c40);
constexpr i128 kWeakestMargin = static_cast<i128>(UINT64_C(43848255939894));
constexpr i64 kWeakestGrid = INT64_C(33581974485739680);

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

u32 parse_mask(const std::string& token) {
    std::size_t used = 0;
    const u64 wide = std::stoull(token, &used, 16);
    require(used == token.size() && wide < (UINT64_C(1) << 30),
            "bad mask token");
    const u32 mask = static_cast<u32>(wide);
    require(std::popcount(mask) == 8, "mask is not rank eight");
    return mask;
}

std::vector<u32> read_masks(const std::string& path,
                            std::size_t expected_count,
                            u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mask ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    Fnv ledger;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask(token);
        require(distinct.insert(mask).second, "duplicate mask in ledger");
        masks.push_back(mask);
        ledger.add(mask);
    }
    require(input.eof(), "mask-ledger read failed");
    require(masks.size() == expected_count && ledger.state == expected_fnv,
            "mask-ledger identity changed");
    return masks;
}

u64 mask_fnv(const std::vector<u32>& masks) {
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

struct Pair {
    int q = 0;
    int r = 0;
    friend bool operator<(Pair left, Pair right) {
        return std::tie(left.q, left.r) < std::tie(right.q, right.r);
    }
    friend bool operator==(Pair, Pair) = default;
};

std::vector<Pair> read_fibre(const std::string& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open fibre ledger");
    std::vector<Pair> pairs;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank fibre row");
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed fibre row");
        std::size_t used_q = 0;
        std::size_t used_r = 0;
        Pair pair;
        pair.q = std::stoi(line.substr(0, comma), &used_q);
        pair.r = std::stoi(line.substr(comma + 1), &used_r);
        require(used_q == comma && used_r == line.size() - comma - 1 &&
                    pair.q > 0 && pair.q < pair.r &&
                    (pairs.empty() || pairs.back() < pair),
                "invalid or out-of-order fibre pair");
        pairs.push_back(pair);
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(input.eof(), "fibre-ledger read failed");
    require(pairs.size() == kFibreCount && ledger.state == kFibreFnv &&
                std::find(pairs.begin(), pairs.end(), Pair{294, 636}) !=
                    pairs.end(),
            "fibre identity changed");
    return pairs;
}

i64 checked_lcm(i64 left, i64 right) {
    require(left > 0 && right > 0, "nonpositive LCM input");
    const i128 wide = static_cast<i128>(left / std::gcd(left, right)) * right;
    require(wide <= std::numeric_limits<i64>::max(),
            "literal grid overflows i64");
    return static_cast<i64>(wide);
}

i64 fixed_grid() {
    i64 grid = 1;
    for (int speed : kPool) grid = checked_lcm(grid, 14LL * speed);
    require(grid == kFixedGrid, "fixed-pool literal grid changed");
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
    std::map<u32, i64> safe_width_by_failure;
};

Geometry build_geometry(Pair pair) {
    i64 grid = fixed_grid();
    grid = checked_lcm(grid, 14LL * pair.q);
    grid = checked_lcm(grid, 14LL * pair.r);
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
    add_walls(pair.q);
    add_walls(pair.r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.front() == 0 && walls.back() == grid,
            "literal wall boundary changed");

    Geometry geometry;
    geometry.grid = grid;
    geometry.cells = walls.size() - 1;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(pair.q, grid, left, right) ||
            !safe_midpoint(pair.r, grid, left, right))
            continue;
        u32 failure = 0;
        for (unsigned vertex = 0; vertex < kPool.size(); ++vertex)
            if (!safe_midpoint(kPool[vertex], grid, left, right))
                failure |= u32{1} << vertex;
        const i64 width = right - left;
        geometry.pair_ticks += width;
        geometry.safe_width_by_failure[failure] += width;
    }
    return geometry;
}

u32 next_combination(u32 value) {
    const u32 low = value & (u32{0} - value);
    require(low != 0, "zero low bit in combination successor");
    const u32 ripple = value + low;
    return ripple | (((ripple ^ value) >> 2) / low);
}

struct BodyAudit {
    u64 bodies = 0;
    u64 checks = 0;
    u64 maximum_prefix = 0;
    u64 failures = 0;
    u64 row_fnv = kEmptyFnv;
};

BodyAudit audit_bodies(const std::vector<u32>& deck) {
    BodyAudit audit;
    Fnv ledger;
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++audit.bodies;
        u64 prefix = 0;
        bool covered = false;
        for (u32 mask : deck) {
            ++prefix;
            ++audit.checks;
            if ((body & mask) == 0) {
                covered = true;
                break;
            }
        }
        if (!covered) {
            ++audit.failures;
            prefix = 0;
        } else {
            audit.maximum_prefix = std::max(audit.maximum_prefix, prefix);
        }
        ledger.add(body);
        ledger.add(prefix);
    }
    audit.row_fnv = ledger.state;
    require(audit.bodies == kBodyCount && audit.checks == kBodyChecks &&
                audit.maximum_prefix == kRebuiltCount &&
                audit.failures == 0 && audit.row_fnv == kBodyRowFnv,
            "rebuilt exhaustive body scan changed");
    return audit;
}

struct PairAudit {
    Pair pair;
    u64 cells = 0;
    u64 failure_classes = 0;
    i64 grid = 0;
    i64 pair_ticks = 0;
    i128 minimum_margin = 0;
    i128 maximum_margin = 0;
    u32 minimum_mask = 0;
    u64 equalities = 0;
    u64 ledger_fnv = kEmptyFnv;
};

PairAudit audit_pair(Pair pair, const std::vector<u32>& deck) {
    const Geometry geometry = build_geometry(pair);
    PairAudit audit;
    audit.pair = pair;
    audit.cells = geometry.cells;
    audit.failure_classes = geometry.safe_width_by_failure.size();
    audit.grid = geometry.grid;
    audit.pair_ticks = geometry.pair_ticks;
    bool first = true;
    Fnv ledger;
    for (u32 mask : deck) {
        i64 mass = 0;
        for (const auto& [failure, width] : geometry.safe_width_by_failure)
            if ((failure & ~mask) == 0) mass += width;
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * geometry.grid;
        audit.equalities += margin == 0;
        require(margin > 0, "rebuilt mask is not strictly active");
        if (first || margin < audit.minimum_margin) {
            audit.minimum_margin = margin;
            audit.minimum_mask = mask;
        }
        if (first || margin > audit.maximum_margin)
            audit.maximum_margin = margin;
        first = false;
        ledger.add(pair.q);
        ledger.add(pair.r);
        ledger.add(mask);
        add_i128(ledger, margin);
    }
    audit.ledger_fnv = ledger.state;
    return audit;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 4,
                "usage: detached-signature294 JOINT421 FIBRE21 WITNESS4");
        const std::vector<u32> joint =
            read_masks(argv[1], kJointCount, kJointFnv);
        const std::vector<Pair> fibre = read_fibre(argv[2]);
        const std::vector<u32> witness =
            read_masks(argv[3], kWitnessCount, kWitnessFnv);
        require(std::equal(witness.begin(), witness.end(),
                           kExpectedWitness.begin()),
                "two-mask witness identity changed");
        for (std::size_t i = 0; i < kDeletedIndices.size(); ++i)
            require(joint[kDeletedIndices[i]] == kDeletedMasks[i],
                    "deleted-mask indices changed");

        std::vector<u32> rebuilt;
        rebuilt.reserve(kRebuiltCount);
        std::set<u32> distinct;
        for (std::size_t index = 0; index < joint.size(); ++index) {
            if (std::find(kDeletedIndices.begin(), kDeletedIndices.end(),
                          index) != kDeletedIndices.end())
                continue;
            require(distinct.insert(joint[index]).second,
                    "retained deck duplicate");
            rebuilt.push_back(joint[index]);
        }
        for (u32 mask : witness) {
            require(distinct.insert(mask).second,
                    "witness overlaps retained deck");
            rebuilt.push_back(mask);
        }
        require(rebuilt.size() == kRebuiltCount &&
                    mask_fnv(rebuilt) == kRebuiltFnv,
                "rebuilt deck identity changed");

        const BodyAudit bodies = audit_bodies(rebuilt);
        Fnv aggregate;
        u64 total_cells = 0;
        u64 total_classes = 0;
        u64 total_equalities = 0;
        bool weakest_set = false;
        Pair weakest_pair;
        u32 weakest_mask = 0;
        i128 weakest_margin = 0;
        i64 weakest_grid = 1;

        std::cout << "LRC14_SIGNATURE294_DETACHED_LITERAL_AUDIT_V1\n"
                  << "JOINT " << joint.size() << " FNV " << std::hex
                  << kJointFnv << " FIBRE " << std::dec << fibre.size()
                  << " FNV " << std::hex << kFibreFnv << " WITNESS "
                  << std::dec << witness.size() << " FNV " << std::hex
                  << kWitnessFnv << std::dec << '\n'
                  << "DELETE_INDICES 107 173 247 DELETE_MASKS " << std::hex
                  << std::setw(8) << std::setfill('0') << kDeletedMasks[0]
                  << ' ' << std::setw(8) << kDeletedMasks[1] << ' '
                  << std::setw(8) << kDeletedMasks[2] << " APPEND";
        for (u32 mask : witness)
            std::cout << ' ' << std::setw(8) << mask;
        std::cout << std::dec << std::setfill(' ') << '\n'
                  << "REBUILT " << rebuilt.size() << " FNV " << std::hex
                  << mask_fnv(rebuilt) << std::dec << " BODY_UNIVERSE "
                  << bodies.bodies << " CHECKS " << bodies.checks
                  << " MAX_PREFIX " << bodies.maximum_prefix
                  << " FAILURES " << bodies.failures << " ROW_FNV "
                  << std::hex << bodies.row_fnv << std::dec << '\n';

        for (Pair pair : fibre) {
            const PairAudit audit = audit_pair(pair, rebuilt);
            total_cells += audit.cells;
            total_classes += audit.failure_classes;
            total_equalities += audit.equalities;
            aggregate.add(pair.q);
            aggregate.add(pair.r);
            aggregate.add(audit.grid);
            aggregate.add(audit.cells);
            aggregate.add(audit.failure_classes);
            aggregate.add(audit.pair_ticks);
            add_i128(aggregate, audit.minimum_margin);
            add_i128(aggregate, audit.maximum_margin);
            aggregate.add(audit.minimum_mask);
            aggregate.add(audit.equalities);
            aggregate.add(audit.ledger_fnv);
            if (!weakest_set ||
                audit.minimum_margin * weakest_grid <
                    weakest_margin * audit.grid) {
                weakest_set = true;
                weakest_pair = pair;
                weakest_mask = audit.minimum_mask;
                weakest_margin = audit.minimum_margin;
                weakest_grid = audit.grid;
            }
            std::cout << "PAIR " << pair.q << ',' << pair.r << " GRID "
                      << audit.grid << " CELLS " << audit.cells
                      << " FAILURE_CLASSES " << audit.failure_classes
                      << " PAIR_TICKS " << audit.pair_ticks
                      << " ACTIVITY_TESTS " << rebuilt.size()
                      << " EQUALITIES " << audit.equalities
                      << " MIN_MARGIN " << decimal(audit.minimum_margin)
                      << " MIN_MASK " << std::hex << std::setw(8)
                      << std::setfill('0') << audit.minimum_mask << std::dec
                      << std::setfill(' ') << " MAX_MARGIN "
                      << decimal(audit.maximum_margin) << " LEDGER_FNV "
                      << std::hex << audit.ledger_fnv << std::dec << '\n';
        }
        require(total_cells == kTotalCells &&
                    total_classes == kTotalFailureClasses &&
                    total_equalities == 0 &&
                    aggregate.state == kActivityAggregateFnv &&
                    weakest_pair == Pair{kWeakestQ, kWeakestR} &&
                    weakest_mask == kWeakestMask &&
                    weakest_margin == kWeakestMargin &&
                    weakest_grid == kWeakestGrid,
                "detached activity aggregate changed");
        std::cout << "ACTIVITY PAIRS " << fibre.size() << " MASKS_PER_PAIR "
                  << rebuilt.size() << " TESTS "
                  << fibre.size() * rebuilt.size() << " TOTAL_CELLS "
                  << total_cells << " TOTAL_FAILURE_CLASSES "
                  << total_classes << " EQUALITIES " << total_equalities
                  << " AGGREGATE_FNV " << std::hex << aggregate.state
                  << std::dec << '\n'
                  << "WEAKEST PAIR " << weakest_pair.q << ','
                  << weakest_pair.r << " MASK " << std::hex << std::setw(8)
                  << std::setfill('0') << weakest_mask << std::dec
                  << std::setfill(' ') << " MARGIN_NUM "
                  << decimal(weakest_margin) << " MARGIN_DEN "
                  << weakest_grid << '\n'
                  << "VERDICT PASS DETACHED_EXACT_21_ROW_COMMON_DECK_SURGERY\n"
                  << "SCOPE LITERAL_ACTIVITY_AND_EXHAUSTIVE_BODY_SCAN_"
                     "FIXED_POOL_NO_PHYSICAL_ENTRY_NO_LRC14\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "LRC14_SIGNATURE294_DETACHED_AUDIT_ERROR "
                  << error.what() << '\n';
        return 1;
    }
}
