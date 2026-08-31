// Import-free literal-wall audit of the recursive exact inactive-signature
// ideal {19} and its minimum-two common-deck surgery.

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
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
constexpr std::size_t kDeckCount = 421;
constexpr std::size_t kTargetIndex = 19;
constexpr u32 kOldMask = UINT32_C(0x1804aa01);
constexpr u64 kDeckFnv = UINT64_C(0x20d63dd42fe8150e);
constexpr std::size_t kAtlasRows = 24223;
constexpr std::size_t kLiveCount = 22647;
constexpr u64 kLiveFnv = UINT64_C(0xdf5374d4aca67677);
constexpr std::size_t kIdealCount = 36;
constexpr u64 kIdealFnv = UINT64_C(0x5c8af37cf2f002e7);
constexpr u64 kRepairCount = UINT64_C(5852925);
constexpr u64 kBodyCount = UINT64_C(14307150);
constexpr i64 kFixedGrid = INT64_C(18241159416480);
constexpr std::array<u32, 2> kWitness = {
    UINT32_C(0x0000e649), UINT32_C(0x0184a205)};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
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
    std::string out;
    while (magnitude != 0) {
        out.push_back(static_cast<char>('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}

struct U256 {
    std::array<u64, 4> limb{};
};

U256 multiply_u128(u128 left, u128 right) {
    const std::array<u64, 2> a = {
        static_cast<u64>(left), static_cast<u64>(left >> 64)};
    const std::array<u64, 2> b = {
        static_cast<u64>(right), static_cast<u64>(right >> 64)};
    U256 result;
    for (unsigned i = 0; i < 2; ++i) {
        u128 carry = 0;
        for (unsigned j = 0; j < 2; ++j) {
            const u128 cell = static_cast<u128>(a[i]) * b[j] +
                              result.limb[i + j] + carry;
            result.limb[i + j] = static_cast<u64>(cell);
            carry = cell >> 64;
        }
        result.limb[i + 2] = static_cast<u64>(carry);
    }
    return result;
}

bool fraction_less(i128 left_num, i128 left_den,
                   i128 right_num, i128 right_den) {
    require(left_num >= 0 && left_den > 0 && right_num >= 0 && right_den > 0,
            "invalid fraction comparison");
    const U256 left = multiply_u128(static_cast<u128>(left_num),
                                    static_cast<u128>(right_den));
    const U256 right = multiply_u128(static_cast<u128>(right_num),
                                     static_cast<u128>(left_den));
    for (int index = 3; index >= 0; --index) {
        if (left.limb[index] != right.limb[index])
            return left.limb[index] < right.limb[index];
    }
    return false;
}

i128 gcd128(i128 left, i128 right) {
    if (left < 0) left = -left;
    if (right < 0) right = -right;
    while (right != 0) {
        const i128 remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

std::vector<std::string> split(const std::string& line, char delimiter) {
    std::vector<std::string> out;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, delimiter)) out.push_back(field);
    return out;
}

struct Pair {
    int q = 0;
    int r = 0;
    friend bool operator<(Pair left, Pair right) {
        return std::tie(left.q, left.r) < std::tie(right.q, right.r);
    }
    friend bool operator==(Pair left, Pair right) {
        return left.q == right.q && left.r == right.r;
    }
};

std::vector<u32> read_deck(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open joint deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    Fnv ledger;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "invalid deck token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "invalid or repeated deck mask");
        deck.push_back(mask);
        ledger.add(mask);
    }
    require(input.eof() && deck.size() == kDeckCount &&
                ledger.state == kDeckFnv && deck[kTargetIndex] == kOldMask,
            "joint deck identity changed");
    return deck;
}

std::set<Pair> read_live(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open live residual");
    std::set<Pair> live;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank live row");
        const std::vector<std::string> fields = split(line, ',');
        require(fields.size() == 2, "malformed live row");
        const Pair pair{std::stoi(fields[0]), std::stoi(fields[1])};
        require(pair.q > 0 && pair.q < pair.r && live.insert(pair).second,
                "invalid or repeated live row");
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(input.eof() && live.size() == kLiveCount &&
                ledger.state == kLiveFnv,
            "THM4287 residual identity changed");
    return live;
}

std::vector<Pair> read_ideal(const std::filesystem::path& path,
                             const std::set<Pair>& live) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open signature atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::vector<Pair> ideal;
    Fnv ledger;
    std::size_t rows = 0;
    while (std::getline(input, line)) {
        ++rows;
        const std::vector<std::string> fields = split(line, ',');
        require(fields.size() == 10, "malformed signature row");
        const Pair pair{std::stoi(fields[0]), std::stoi(fields[1])};
        if (!live.contains(pair) || std::stoul(fields[2]) != 1) continue;
        bool exact = true;
        for (std::size_t word = 0; word < 7; ++word) {
            const u64 value = std::stoull(fields[3 + word], nullptr, 16);
            const u64 expected = word == 0
                ? (UINT64_C(1) << kTargetIndex) : 0;
            exact &= value == expected;
        }
        if (!exact) continue;
        ideal.push_back(pair);
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(input.eof() && rows == kAtlasRows && ideal.size() == kIdealCount &&
                ledger.state == kIdealFnv &&
                std::find(ideal.begin(), ideal.end(), Pair{338, 636}) !=
                    ideal.end(),
            "exact signature ideal {19} changed");
    return ideal;
}

i64 checked_lcm(i64 left, i64 right) {
    require(left > 0 && right > 0, "nonpositive LCM input");
    const i128 wide = static_cast<i128>(left / std::gcd(left, right)) * right;
    require(wide <= std::numeric_limits<i64>::max(), "literal grid overflow");
    return static_cast<i64>(wide);
}

i64 fixed_grid() {
    i64 grid = 1;
    for (int speed : kPool) grid = checked_lcm(grid, 14LL * speed);
    require(grid == kFixedGrid, "fixed pool grid changed");
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
    std::map<u32, i64> widths;
};

Geometry build_geometry(Pair pair) {
    i64 grid = fixed_grid();
    grid = checked_lcm(grid, 14LL * pair.q);
    grid = checked_lcm(grid, 14LL * pair.r);
    std::vector<i64> walls = {0, grid};
    auto add_walls = [&](int speed) {
        const i64 divisor = 14LL * speed;
        require(grid % divisor == 0, "nonintegral literal wall unit");
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
            "literal boundary changed");
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
        geometry.widths[failure] += right - left;
    }
    return geometry;
}

std::array<std::array<u64, 9>, 31> choose{};

void init_choose() {
    for (int n = 0; n <= 30; ++n) {
        choose[n][0] = 1;
        for (int k = 1; k <= 8; ++k)
            choose[n][k] = n == 0 ? 0 :
                choose[n - 1][k] + choose[n - 1][k - 1];
    }
    require(choose[30][8] == kRepairCount, "rank-eight universe changed");
}

u64 colex_rank(u32 mask) {
    u64 rank = 0;
    int ordinal = 1;
    for (int bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        rank += choose[bit][ordinal++];
    }
    require(ordinal == 9 && rank < kRepairCount, "invalid rank-eight mask");
    return rank;
}

u32 next_combination(u32 value) {
    const u32 low = value & (u32{0} - value);
    require(low != 0, "zero combination low bit");
    const u32 ripple = value + low;
    return ripple | (((ripple ^ value) >> 2) / low);
}

void add_supersets(u32 atom, int need, int start, u32 extra, i64 width,
                   std::vector<i64>& masses, u64& operations) {
    if (need == 0) {
        masses[colex_rank(atom | extra)] += width;
        ++operations;
        return;
    }
    for (int bit = start; bit <= 30 - need; ++bit) {
        const u32 flag = u32{1} << bit;
        if ((atom & flag) != 0) continue;
        add_supersets(atom, need - 1, bit + 1, extra | flag, width,
                      masses, operations);
    }
}

struct RowActivity {
    Pair pair;
    i64 grid = 0;
    u64 cells = 0;
    u64 failure_classes = 0;
    u64 operations = 0;
    u64 active_count = 0;
    std::vector<u64> active;
    i128 old_margin = 0;
    std::vector<i128> required_margins;
};

RowActivity compute_activity(Pair pair,
                             const std::vector<u32>& required_masks) {
    const Geometry geometry = build_geometry(pair);
    std::vector<i64> masses(kRepairCount, 0);
    u64 operations = 0;
    for (const auto& [failure, width] : geometry.widths) {
        const int weight = std::popcount(failure);
        if (weight <= 8)
            add_supersets(failure, 8 - weight, 0, 0, width,
                          masses, operations);
    }
    RowActivity out;
    out.pair = pair;
    out.grid = geometry.grid;
    out.cells = geometry.cells;
    out.failure_classes = geometry.widths.size();
    out.operations = operations;
    out.active.assign((kRepairCount + 63) / 64, 0);
    for (u64 rank = 0; rank < kRepairCount; ++rank) {
        const i128 margin = static_cast<i128>(63) * masses[rank] -
                            static_cast<i128>(4) * geometry.grid;
        if (margin >= 0) {
            out.active[rank / 64] |= UINT64_C(1) << (rank % 64);
            ++out.active_count;
        }
    }
    out.old_margin = static_cast<i128>(63) * masses[colex_rank(kOldMask)] -
                     static_cast<i128>(4) * geometry.grid;
    out.required_margins.reserve(required_masks.size());
    for (u32 mask : required_masks)
        out.required_margins.push_back(
            static_cast<i128>(63) * masses[colex_rank(mask)] -
            static_cast<i128>(4) * geometry.grid);
    return out;
}

struct Obligations {
    std::vector<u32> bodies;
    u64 checks = 0;
    u64 fnv = 0;
    u32 body_union = 0;
};

Obligations private_obligations(const std::vector<u32>& deck) {
    Obligations out;
    Fnv ledger;
    u32 body = (u32{1} << 9) - 1;
    for (u64 ordinal = 0; ordinal < kBodyCount; ++ordinal) {
        if ((body & kOldMask) == 0) {
            bool retained_hit = false;
            for (std::size_t index = 0; index < deck.size(); ++index) {
                if (index == kTargetIndex) continue;
                ++out.checks;
                if ((body & deck[index]) == 0) {
                    retained_hit = true;
                    break;
                }
            }
            if (!retained_hit) {
                out.bodies.push_back(body);
                out.body_union |= body;
                ledger.add(body);
            }
        }
        if (ordinal + 1 < kBodyCount) body = next_combination(body);
    }
    require(body == UINT32_C(0x3fe00000) && out.bodies.size() == 8 &&
                out.body_union == UINT32_C(0x27fb548a),
            "old-mask private obligations changed");
    out.fnv = ledger.state;
    return out;
}

struct BodyAudit {
    u64 checks = 0;
    u64 failures = 0;
    u64 fnv = 0;
};

BodyAudit scan_body_cover(const std::vector<u32>& deck) {
    BodyAudit out;
    Fnv ledger;
    u32 body = (u32{1} << 9) - 1;
    for (u64 ordinal = 0; ordinal < kBodyCount; ++ordinal) {
        bool hit = false;
        for (u32 mask : deck) {
            ++out.checks;
            if ((body & mask) == 0) {
                hit = true;
                break;
            }
        }
        out.failures += !hit;
        ledger.add(body);
        ledger.add(hit ? 1 : 0);
        if (ordinal + 1 < kBodyCount) body = next_combination(body);
    }
    require(body == UINT32_C(0x3fe00000), "body scan endpoint changed");
    out.fnv = ledger.state;
    return out;
}

}  // namespace


int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 5,
                "usage: detached-j19 DECK SIGNATURES LIVE THREADS");
        init_choose();
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::set<Pair> live = read_live(argv[3]);
        const std::vector<Pair> ideal = read_ideal(argv[2], live);
        const Obligations obligations = private_obligations(deck);

        std::vector<u32> rebuilt;
        rebuilt.reserve(422);
        for (std::size_t index = 0; index < deck.size(); ++index)
            if (index != kTargetIndex) rebuilt.push_back(deck[index]);
        rebuilt.insert(rebuilt.end(), kWitness.begin(), kWitness.end());
        require(rebuilt.size() == 422 &&
                    std::set<u32>(rebuilt.begin(), rebuilt.end()).size() == 422,
                "rebuilt deck rank/distinctness changed");

        const unsigned thread_count =
            static_cast<unsigned>(std::stoul(argv[4]));
        require(thread_count >= 1 && thread_count <= 8,
                "thread count outside 1..8");
        std::vector<RowActivity> rows(ideal.size());
        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        for (unsigned thread = 0; thread < thread_count; ++thread) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t index = next.fetch_add(1);
                    if (index >= ideal.size()) break;
                    rows[index] = compute_activity(ideal[index], rebuilt);
                }
            });
        }
        for (std::thread& worker : workers) worker.join();

        std::vector<u64> common((kRepairCount + 63) / 64, ~UINT64_C(0));
        if (kRepairCount % 64 != 0)
            common.back() = (UINT64_C(1) << (kRepairCount % 64)) - 1;
        Fnv activity_ledger;
        Fnv strict_ledger;
        u64 strict_cells = 0;
        u64 equalities = 0;
        u64 operation_total = 0;
        for (std::size_t row = 0; row < rows.size(); ++row) {
            require(rows[row].pair == ideal[row] && rows[row].old_margin < 0,
                    "old index 19 is not strictly inactive on ideal row");
            require(rows[row].required_margins.size() == rebuilt.size(),
                    "required margin vector changed");
            for (std::size_t index = 0; index < rebuilt.size(); ++index) {
                const i128 margin = rows[row].required_margins[index];
                require(margin > 0,
                        "rebuilt deck mask is not strictly active on ideal");
                ++strict_cells;
                equalities += margin == 0;
                strict_ledger.add(ideal[row].q);
                strict_ledger.add(ideal[row].r);
                strict_ledger.add(rebuilt[index]);
                add_i128(strict_ledger, margin);
            }
            for (std::size_t word = 0; word < common.size(); ++word)
                common[word] &= rows[row].active[word];
            operation_total += rows[row].operations;
            activity_ledger.add(ideal[row].q);
            activity_ledger.add(ideal[row].r);
            activity_ledger.add(rows[row].grid);
            activity_ledger.add(rows[row].cells);
            activity_ledger.add(rows[row].failure_classes);
            activity_ledger.add(rows[row].operations);
            activity_ledger.add(rows[row].active_count);
            add_i128(activity_ledger, rows[row].old_margin);
        }

        struct ResponseClass { u64 count = 0; u32 least = 0; };
        std::map<unsigned, ResponseClass> classes;
        Fnv common_ledger;
        Fnv response_ledger;
        u64 common_count = 0;
        u32 repair = (u32{1} << 8) - 1;
        for (u64 rank = 0; rank < kRepairCount; ++rank) {
            if ((common[rank / 64] & (UINT64_C(1) << (rank % 64))) != 0) {
                ++common_count;
                common_ledger.add(repair);
                unsigned response = 0;
                for (unsigned index = 0; index < obligations.bodies.size();
                     ++index)
                    if ((repair & obligations.bodies[index]) == 0)
                        response |= 1u << index;
                response_ledger.add(repair);
                response_ledger.add(response);
                ResponseClass& cls = classes[response];
                ++cls.count;
                if (cls.count == 1 || repair < cls.least) cls.least = repair;
            }
            if (rank + 1 < kRepairCount) repair = next_combination(repair);
        }
        require(repair == UINT32_C(0x3fc00000) &&
                    common_count == UINT64_C(2212775) &&
                    common_ledger.state == UINT64_C(0x8133b80b94da565c) &&
                    response_ledger.state == UINT64_C(0x0fc1d9bd1c626c8f) &&
                    classes.size() == 41 && !classes.contains(0xff),
                "common response quotient changed");

        unsigned witness_union = 0;
        Fnv witness_margin_ledger;
        for (u32 witness : kWitness) {
            unsigned response = 0;
            for (unsigned index = 0; index < obligations.bodies.size(); ++index)
                if ((witness & obligations.bodies[index]) == 0)
                    response |= 1u << index;
            require(common[colex_rank(witness) / 64] &
                        (UINT64_C(1) << (colex_rank(witness) % 64)),
                    "stated witness is not common-active");
            witness_union |= response;
            bool first = true;
            i128 min_num = 0;
            i128 min_den = 1;
            Pair weakest{};
            for (std::size_t row = 0; row < rows.size(); ++row) {
                const auto found = std::find(rebuilt.begin(), rebuilt.end(), witness);
                require(found != rebuilt.end(), "witness absent from rebuilt deck");
                const std::size_t rebuilt_index = found - rebuilt.begin();
                const i128 margin = rows[row].required_margins[rebuilt_index];
                const i128 denominator = static_cast<i128>(63) * rows[row].grid;
                // The row grids approach 64-bit scale, so the cross-products
                // can exceed signed 128-bit range even though each individual
                // margin and denominator fits.  Compare in an unbounded ring.
                if (first || fraction_less(margin, denominator,
                                           min_num, min_den)) {
                    min_num = margin;
                    min_den = denominator;
                    weakest = ideal[row];
                    first = false;
                }
                witness_margin_ledger.add(witness);
                witness_margin_ledger.add(ideal[row].q);
                witness_margin_ledger.add(ideal[row].r);
                add_i128(witness_margin_ledger, margin);
                add_i128(witness_margin_ledger, denominator);
            }
            const i128 divisor = gcd128(min_num, min_den);
            std::cout << "WITNESS " << std::hex << std::setw(8)
                      << std::setfill('0') << witness << " RESPONSE "
                      << std::setw(2) << response << std::dec
                      << std::setfill(' ') << " WEAKEST " << weakest.q << ','
                      << weakest.r << " GAP " << decimal(min_num / divisor)
                      << '/' << decimal(min_den / divisor) << '\n';
        }
        require(witness_union == 0xff,
                "two witnesses do not cover all private obligations");

        Fnv rebuilt_ledger;
        for (u32 mask : rebuilt) rebuilt_ledger.add(mask);
        require(rebuilt_ledger.state == UINT64_C(0xdc50478119bc6c12),
                "rebuilt deck identity changed");
        const BodyAudit body_audit = scan_body_cover(rebuilt);
        require(body_audit.failures == 0 &&
                    body_audit.fnv == UINT64_C(0x10c4c3ed46d44bf1),
                "rebuilt body cover changed");

        std::cout << "DETACHED_J19_COMMON_DECK_AUDIT_V1\n"
                  << "IDEAL {19} ROWS " << ideal.size() << " FNV "
                  << std::hex << kIdealFnv << std::dec << " ROWS";
        for (Pair pair : ideal) std::cout << ' ' << pair.q << ',' << pair.r;
        std::cout << '\n'
                  << "OLD INDEX 19 MASK " << std::hex << std::setw(8)
                  << std::setfill('0') << kOldMask << std::dec
                  << std::setfill(' ') << " OBLIGATIONS "
                  << obligations.bodies.size() << " FNV " << std::hex
                  << obligations.fnv << " UNION " << obligations.body_union
                  << std::dec << " RETAINED_CHECKS " << obligations.checks
                  << " BODIES";
        for (u32 body : obligations.bodies)
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << body << std::dec
                      << std::setfill(' ');
        std::cout << '\n'
                  << "ACTIVITY ROWS " << rows.size() << " REPAIRS_PER_ROW "
                  << kRepairCount << " ZETA_OPERATIONS " << operation_total
                  << " ACTIVITY_FNV " << std::hex << activity_ledger.state
                  << std::dec << '\n'
                  << "COMMON_ACTIVE " << common_count << " FNV " << std::hex
                  << common_ledger.state << " RESPONSE_FNV "
                  << response_ledger.state << std::dec << " CLASSES "
                  << classes.size() << " FULL_RESPONDERS 0 EXACT_MINIMUM 2\n"
                  << "STRICT_REBUILT_CELLS " << strict_cells
                  << " EQUALITIES " << equalities << " FNV " << std::hex
                  << strict_ledger.state << " WITNESS_MARGIN_FNV "
                  << witness_margin_ledger.state << std::dec << '\n'
                  << "REBUILT_DECK " << rebuilt.size() << " FNV " << std::hex
                  << rebuilt_ledger.state << std::dec << " BODY_SCAN "
                  << kBodyCount << " CHECKS " << body_audit.checks
                  << " FAILURES " << body_audit.failures << " BODY_FNV "
                  << std::hex << body_audit.fnv << std::dec << '\n'
                  << "SCOPE DETACHED_LITERAL_FIXED_POOL_COMMON_DECK_NO_"
                     "PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_J19_MINIMUM_TWO_COMMON_DECK\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "DETACHED_J19_ERROR " << error.what() << '\n';
        return 1;
    }
}
