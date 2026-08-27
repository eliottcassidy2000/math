// Detached exact audit of the explicit THM-4281 delete-seven/add-eight
// surgery branch.  This translation unit imports no project source and uses
// neither primitive/cocycle nor atom/zeta activity code.

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
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};

// These are the seven deleted masks in the response-pattern bit order.
constexpr std::array<u32, 7> DELETED = {
    UINT32_C(0x21948006), UINT32_C(0x128c8900),
    UINT32_C(0x10259240), UINT32_C(0x20027016),
    UINT32_C(0x10a41016), UINT32_C(0x2051d200),
    UINT32_C(0x01188550)};

// The final eight lines of surgery422_masks.txt, treated only as an explicit
// witness family.  This audit makes no optimality or minimum-cardinality claim.
constexpr std::array<u32, 8> APPENDED = {
    UINT32_C(0x0110a550), UINT32_C(0x04871108),
    UINT32_C(0x10241207), UINT32_C(0x1042d088),
    UINT32_C(0x12848902), UINT32_C(0x21141284),
    UINT32_C(0x30249140), UINT32_C(0x31202206)};

constexpr std::size_t RETAINED_COUNT = 414;
constexpr std::size_t DECK_COUNT = 422;
constexpr std::size_t COMMON_COUNT = 188;
constexpr std::size_t RESIDUAL_COUNT = 24223;
constexpr std::size_t COMPLEMENT_COUNT = 24035;
constexpr std::size_t OBLIGATION_COUNT = 71;
constexpr u64 BODY_COUNT = UINT64_C(14307150);
constexpr i64 FIXED_GRID = INT64_C(18241159416480);

constexpr u64 EXPECTED_DECK_FNV = UINT64_C(0x1c97b54ece61b351);
constexpr u64 EXPECTED_APPENDED_FNV = UINT64_C(0x3a69f2f733d2309c);
constexpr u64 EXPECTED_COMMON_FNV = UINT64_C(0x6588121dbec57bcb);
constexpr u64 EXPECTED_RESIDUAL_FNV = UINT64_C(0x80ec0687d8c7dba7);
constexpr u64 EXPECTED_COMPLEMENT_FNV = UINT64_C(0x5d900a28f7d280c9);
constexpr u64 EXPECTED_BODY_FNV = UINT64_C(0x414d30a143d2e22d);
constexpr u64 EXPECTED_OBLIGATION_FNV = UINT64_C(0x287281de2900cca7);

constexpr u64 EXPECTED_ACTIVITY_FNV = UINT64_C(0xf68329205d16d97d);
constexpr u64 EXPECTED_GEOMETRY_FNV = UINT64_C(0xaf23d48e2d17ca24);
constexpr u64 EXPECTED_REPLACEMENT_RESPONSE_FNV =
    UINT64_C(0xb9c73d1da397988d);

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

struct Fnv {
    u64 state = UINT64_C(0xcbf29ce484222325);
    void add(u64 word) {
        for (unsigned shift = 0; shift < 64; shift += 8) {
            state ^= (word >> shift) & UINT64_C(0xff);
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

u32 parse_hex32(const std::string& token) {
    std::size_t used = 0;
    const u64 wide = std::stoull(token, &used, 16);
    require(used == token.size() && wide < (UINT64_C(1) << 32),
            "bad hexadecimal token: " + token);
    return static_cast<u32>(wide);
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream input(line);
    std::string field;
    while (std::getline(input, field, ',')) fields.push_back(field);
    return fields;
}

struct Pair {
    int q = 0;
    int r = 0;
    friend bool operator<(const Pair& left, const Pair& right) {
        return std::tie(left.q, left.r) < std::tie(right.q, right.r);
    }
    friend bool operator==(const Pair& left, const Pair& right) {
        return left.q == right.q && left.r == right.r;
    }
};

struct PairRows {
    std::vector<Pair> rows;
    u64 fnv = 0;
};

PairRows read_pair_rows(const std::filesystem::path& path,
                        std::size_t expected_count,
                        u64 expected_fnv,
                        const std::string& label) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open " + label);
    PairRows out;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty() && line.back() != '\r',
                label + " has blank/CRLF row");
        const std::vector<std::string> fields = split_csv(line);
        require(fields.size() == 2, "malformed " + label + " row");
        Pair row{std::stoi(fields[0]), std::stoi(fields[1])};
        require(row.q > 0 && row.q < row.r,
                "invalid ordered pair in " + label);
        require(out.rows.empty() || out.rows.back() < row,
                label + " is not strictly lexicographically ordered");
        out.rows.push_back(row);
        ledger.add(static_cast<u64>(row.q));
        ledger.add(static_cast<u64>(row.r));
    }
    require(out.rows.size() == expected_count,
            label + " count changed");
    require(ledger.state == expected_fnv, label + " FNV changed");
    out.fnv = ledger.state;
    return out;
}

std::vector<u32> read_deck(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open surgery deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    Fnv ledger;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_hex32(token);
        require(mask < (u32{1} << 30) && std::popcount(mask) == 8,
                "deck mask has wrong range/arity");
        require(distinct.insert(mask).second, "duplicate deck mask");
        deck.push_back(mask);
        ledger.add(mask);
    }
    require(deck.size() == DECK_COUNT, "surgery deck count changed");
    require(ledger.state == EXPECTED_DECK_FNV,
            "surgery deck FNV changed");
    Fnv appended_ledger;
    for (std::size_t i = 0; i < APPENDED.size(); ++i) {
        require(deck[RETAINED_COUNT + i] == APPENDED[i],
                "appended witness order changed");
        appended_ledger.add(APPENDED[i]);
    }
    require(appended_ledger.state == EXPECTED_APPENDED_FNV,
            "appended witness FNV changed");
    for (u32 removed : DELETED)
        require(distinct.count(removed) == 0,
                "deleted mask remains in surgery deck");
    return deck;
}

struct Obligation {
    u32 body = 0;
    unsigned response = 0;
};

std::vector<Obligation> read_obligations(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open vulnerable obligations");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "body,response_pattern,full_carrier_subset",
            "vulnerable-obligation header changed");
    std::vector<Obligation> rows;
    Fnv body_ledger;
    Fnv row_ledger;
    while (std::getline(input, line)) {
        require(!line.empty() && line.back() != '\r',
                "bad vulnerable-obligation line ending");
        const std::vector<std::string> fields = split_csv(line);
        require(fields.size() == 3 && fields[2] == "0",
                "malformed vulnerable-obligation row");
        const u32 body = parse_hex32(fields[0]);
        const unsigned response = parse_hex32(fields[1]);
        require(body < (u32{1} << 30) && std::popcount(body) == 9,
                "obligation body has wrong range/arity");
        require(response > 0 && response < 128,
                "obligation response has wrong range");
        require(rows.empty() || rows.back().body < body,
                "obligation bodies are not strictly ordered");
        rows.push_back({body, response});
        body_ledger.add(body);
        row_ledger.add(body);
        row_ledger.add(response);
    }
    require(rows.size() == OBLIGATION_COUNT,
            "vulnerable-obligation count changed");
    require(body_ledger.state == EXPECTED_BODY_FNV,
            "vulnerable body FNV changed");
    require(row_ledger.state == EXPECTED_OBLIGATION_FNV,
            "typed vulnerable-obligation FNV changed");
    return rows;
}

i64 checked_lcm(i64 left, i64 right) {
    require(left > 0 && right > 0, "nonpositive LCM input");
    const i128 wide = static_cast<i128>(left / std::gcd(left, right)) * right;
    require(wide <= std::numeric_limits<i64>::max(),
            "literal joint grid overflows i64");
    return static_cast<i64>(wide);
}

i64 fixed_grid() {
    i64 grid = 1;
    for (int speed : POOL) grid = checked_lcm(grid, 14LL * speed);
    require(grid == FIXED_GRID, "fixed-pool literal grid changed");
    return grid;
}

bool literal_safe_at_midpoint(int speed, i64 grid,
                              i64 left, i64 right) {
    i128 residue = static_cast<i128>(speed) *
                   (static_cast<i128>(left) + right);
    residue %= static_cast<i128>(2) * grid;
    if (residue < 0) residue += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * residue >= grid &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * grid;
}

constexpr std::size_t WORDS = (DECK_COUNT + 63) / 64;
using Bits = std::array<u64, WORDS>;

struct MaskIndex {
    std::array<Bits, 30> contains{};
    Bits all{};
};

MaskIndex build_mask_index(const std::vector<u32>& deck) {
    MaskIndex out;
    for (std::size_t i = 0; i < deck.size(); ++i) {
        out.all[i / 64] |= UINT64_C(1) << (i % 64);
        for (unsigned vertex = 0; vertex < 30; ++vertex)
            if ((deck[i] >> vertex) & 1U)
                out.contains[vertex][i / 64] |= UINT64_C(1) << (i % 64);
    }
    return out;
}

struct Geometry {
    Pair pair;
    i64 grid = 0;
    u64 cells = 0;
    i64 pair_ticks = 0;
    std::map<u32, i64> safe_width_by_failure;
};

Geometry build_geometry(Pair pair) {
    require(pair.q > 0 && pair.q < pair.r,
            "invalid pair for literal geometry");
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
    for (int speed : POOL) add_walls(speed);
    add_walls(pair.q);
    add_walls(pair.r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(!walls.empty() && walls.front() == 0 && walls.back() == grid,
            "literal wall boundary changed");

    Geometry out;
    out.pair = pair;
    out.grid = grid;
    out.cells = walls.size() - 1;
    for (std::size_t i = 1; i < walls.size(); ++i) {
        const i64 left = walls[i - 1];
        const i64 right = walls[i];
        require(left < right, "nonpositive literal cell");
        if (!literal_safe_at_midpoint(pair.q, grid, left, right) ||
            !literal_safe_at_midpoint(pair.r, grid, left, right))
            continue;
        u32 failed = 0;
        for (unsigned vertex = 0; vertex < POOL.size(); ++vertex)
            if (!literal_safe_at_midpoint(POOL[vertex], grid, left, right))
                failed |= u32{1} << vertex;
        const i64 width = right - left;
        out.pair_ticks += width;
        out.safe_width_by_failure[failed] += width;
    }
    return out;
}

struct Classification {
    Geometry geometry;
    std::vector<i64> mass;
    std::vector<i128> margin;
};

Classification classify_all(Pair pair, const std::vector<u32>& deck,
                            const MaskIndex& index) {
    Classification out;
    out.geometry = build_geometry(pair);
    out.mass.assign(deck.size(), 0);
    for (const auto& [failure, width] :
         out.geometry.safe_width_by_failure) {
        if (std::popcount(failure) > 8) continue;
        Bits candidates = index.all;
        u32 remaining = failure;
        while (remaining != 0) {
            const unsigned vertex = std::countr_zero(remaining);
            remaining &= remaining - 1;
            for (std::size_t word = 0; word < WORDS; ++word)
                candidates[word] &= index.contains[vertex][word];
        }
        for (std::size_t word = 0; word < WORDS; ++word) {
            u64 bits = candidates[word];
            while (bits != 0) {
                const unsigned offset = std::countr_zero(bits);
                const std::size_t mask_index = 64 * word + offset;
                require(mask_index < out.mass.size(),
                        "mask bitset tail escaped deck");
                out.mass[mask_index] += width;
                bits &= bits - 1;
            }
        }
    }
    out.margin.reserve(deck.size());
    for (i64 mass : out.mass)
        out.margin.push_back(static_cast<i128>(63) * mass -
                             static_cast<i128>(4) * out.geometry.grid);
    return out;
}

struct LiteralAudit {
    u64 checks = 0;
    u64 equalities = 0;
    Fnv activity;
    Fnv geometry;
    i128 global_min_margin = 0;
    Pair global_min_pair;
    u32 global_min_mask = 0;
    i64 global_min_mass = 0;
    bool minimum_set = false;
};

LiteralAudit audit_common_literal(
    const PairRows& common,
    const std::vector<u32>& deck,
    const MaskIndex& index,
    const std::filesystem::path& output_path,
    const std::filesystem::path& verified_common_path) {
    std::ofstream output(output_path);
    require(static_cast<bool>(output), "cannot create common literal CSV");
    std::ofstream verified_common(verified_common_path);
    require(static_cast<bool>(verified_common),
            "cannot create verified common CSV");
    output << "q,r,grid,cells,pair_ticks,failure_patterns,min_mask,"
              "min_mass,min_margin\n";
    LiteralAudit audit;
    for (const Pair pair : common.rows) {
        const Classification row = classify_all(pair, deck, index);
        audit.geometry.add(static_cast<u64>(pair.q));
        audit.geometry.add(static_cast<u64>(pair.r));
        audit.geometry.add(static_cast<u64>(row.geometry.grid));
        audit.geometry.add(row.geometry.cells);
        audit.geometry.add(static_cast<u64>(row.geometry.pair_ticks));
        audit.geometry.add(row.geometry.safe_width_by_failure.size());

        i128 row_min_margin = 0;
        std::size_t row_min_index = 0;
        bool row_min_set = false;
        for (std::size_t i = 0; i < deck.size(); ++i) {
            ++audit.checks;
            audit.equalities += row.margin[i] == 0;
            require(row.margin[i] >= 0,
                    "claimed common row has literal-inactive mask " +
                    std::to_string(pair.q) + "," +
                    std::to_string(pair.r) + ":" + hex8(deck[i]) +
                    " margin=" + decimal(row.margin[i]));
            audit.activity.add(static_cast<u64>(pair.q));
            audit.activity.add(static_cast<u64>(pair.r));
            audit.activity.add(deck[i]);
            audit.activity.add(static_cast<u64>(row.mass[i]));
            add_i128(audit.activity, row.margin[i]);
            if (!row_min_set || row.margin[i] < row_min_margin ||
                (row.margin[i] == row_min_margin && deck[i] < deck[row_min_index])) {
                row_min_set = true;
                row_min_margin = row.margin[i];
                row_min_index = i;
            }
            if (!audit.minimum_set ||
                row.margin[i] < audit.global_min_margin ||
                (row.margin[i] == audit.global_min_margin &&
                 std::tie(pair.q, pair.r, deck[i]) <
                 std::tie(audit.global_min_pair.q,
                          audit.global_min_pair.r,
                          audit.global_min_mask))) {
                audit.minimum_set = true;
                audit.global_min_margin = row.margin[i];
                audit.global_min_pair = pair;
                audit.global_min_mask = deck[i];
                audit.global_min_mass = row.mass[i];
            }
        }
        output << pair.q << ',' << pair.r << ',' << row.geometry.grid << ','
               << row.geometry.cells << ',' << row.geometry.pair_ticks << ','
               << row.geometry.safe_width_by_failure.size() << ','
               << hex8(deck[row_min_index]) << ',' << row.mass[row_min_index]
               << ',' << decimal(row_min_margin) << '\n';
        verified_common << pair.q << ',' << pair.r << '\n';
    }
    require(output.good(), "failed writing common literal CSV");
    require(verified_common.good(), "failed writing verified common CSV");
    require(audit.checks == COMMON_COUNT * DECK_COUNT,
            "common literal check count changed");
    require(audit.equalities == 0,
            "common literal audit encountered equality");
    require(audit.activity.state == EXPECTED_ACTIVITY_FNV,
            "common literal activity FNV changed");
    require(audit.geometry.state == EXPECTED_GEOMETRY_FNV,
            "common literal geometry FNV changed");
    require(audit.global_min_pair == Pair{51, 420} &&
                audit.global_min_mask == UINT32_C(0x08c0a980) &&
                audit.global_min_mass == INT64_C(1159430348218) &&
                audit.global_min_margin == INT64_C(79474271814),
            "common literal global positive margin changed");
    return audit;
}

u32 next_combination(u32 value) {
    const u32 low = value & (u32{0} - value);
    require(low != 0, "zero low bit in combination successor");
    const u32 ripple = value + low;
    return ripple | (((ripple ^ value) >> 2) / low);
}

struct BodyAudit {
    u64 bodies = 0;
    u64 retained_checks = 0;
    u64 retained_covered = 0;
    u64 exposed = 0;
    u64 replacement_incidences = 0;
    u64 replacement_minimum = std::numeric_limits<u64>::max();
    u64 replacement_maximum = 0;
    Fnv replacement_responses;
};

BodyAudit audit_body_decomposition(
    const std::vector<u32>& deck,
    const std::vector<Obligation>& obligations,
    const std::filesystem::path& output_path) {
    std::ofstream output(output_path);
    require(static_cast<bool>(output),
            "cannot create exposed-body response CSV");
    output << "body,deleted_response,appended_response,appended_hits\n";
    BodyAudit audit;
    std::size_t obligation_index = 0;
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++audit.bodies;
        bool retained_hit = false;
        for (std::size_t i = 0; i < RETAINED_COUNT; ++i) {
            ++audit.retained_checks;
            if ((deck[i] & body) == 0) {
                retained_hit = true;
                break;
            }
        }
        if (retained_hit) {
            ++audit.retained_covered;
            continue;
        }

        ++audit.exposed;
        require(obligation_index < obligations.size(),
                "more exposed bodies than supplied obligations");
        const Obligation& obligation = obligations[obligation_index++];
        require(body == obligation.body,
                "exposed body differs from supplied obligation at index " +
                std::to_string(obligation_index - 1));

        unsigned deleted_response = 0;
        for (std::size_t i = 0; i < DELETED.size(); ++i)
            if ((DELETED[i] & body) == 0)
                deleted_response |= 1U << i;
        require(deleted_response == obligation.response,
                "deleted-mask response pattern changed for body " +
                hex8(body));

        unsigned appended_response = 0;
        u64 hits = 0;
        for (std::size_t i = 0; i < APPENDED.size(); ++i) {
            const bool disjoint = (APPENDED[i] & body) == 0;
            if (disjoint) {
                appended_response |= 1U << i;
                ++hits;
                ++audit.replacement_incidences;
                require((APPENDED[i] & body) == 0,
                        "nondisjoint appended incidence");
            }
        }
        require(appended_response != 0,
                "explicit appended family misses exposed body " +
                hex8(body));
        audit.replacement_minimum =
            std::min(audit.replacement_minimum, hits);
        audit.replacement_maximum =
            std::max(audit.replacement_maximum, hits);
        audit.replacement_responses.add(body);
        audit.replacement_responses.add(deleted_response);
        audit.replacement_responses.add(appended_response);
        audit.replacement_responses.add(hits);
        output << hex8(body) << ',' << std::hex << std::setw(2)
               << std::setfill('0') << deleted_response << ','
               << std::setw(2) << appended_response << std::dec << ','
               << hits << '\n';
    }
    require(output.good(), "failed writing exposed-body response CSV");
    require(audit.bodies == BODY_COUNT, "nine-body universe count changed");
    require(audit.exposed == obligations.size() &&
                obligation_index == obligations.size(),
            "exposed-body set is not exactly the supplied obligations");
    require(audit.retained_covered + audit.exposed == audit.bodies,
            "body decomposition is not exhaustive/disjoint");
    require(audit.replacement_incidences == 82 &&
                audit.replacement_minimum == 1 &&
                audit.replacement_maximum == 2,
            "appended response census changed");
    require(audit.retained_checks == UINT64_C(405775415),
            "retained-deck body search ledger changed");
    require(audit.replacement_responses.state ==
                EXPECTED_REPLACEMENT_RESPONSE_FNV,
            "appended response FNV changed");
    return audit;
}

struct PartitionAudit {
    u64 common = 0;
    u64 complement = 0;
    Fnv common_ledger;
    Fnv complement_ledger;
};

PartitionAudit audit_partition(
    const PairRows& residual,
    const PairRows& common,
    const std::filesystem::path& output_path) {
    std::ofstream output(output_path);
    require(static_cast<bool>(output), "cannot create residual complement");
    PartitionAudit audit;
    std::size_t common_index = 0;
    for (const Pair row : residual.rows) {
        while (common_index < common.rows.size() &&
               common.rows[common_index] < row)
            fail("claimed common row is absent from residual: " +
                 std::to_string(common.rows[common_index].q) + "," +
                 std::to_string(common.rows[common_index].r));
        if (common_index < common.rows.size() &&
            common.rows[common_index] == row) {
            ++audit.common;
            audit.common_ledger.add(static_cast<u64>(row.q));
            audit.common_ledger.add(static_cast<u64>(row.r));
            ++common_index;
        } else {
            ++audit.complement;
            audit.complement_ledger.add(static_cast<u64>(row.q));
            audit.complement_ledger.add(static_cast<u64>(row.r));
            output << row.q << ',' << row.r << '\n';
        }
    }
    require(common_index == common.rows.size(),
            "trailing claimed common row is absent from residual");
    require(output.good(), "failed writing residual complement");
    require(audit.common == COMMON_COUNT &&
                audit.complement == COMPLEMENT_COUNT &&
                audit.common + audit.complement == residual.rows.size(),
            "residual partition counts changed");
    require(audit.common_ledger.state == EXPECTED_COMMON_FNV,
            "partition common FNV changed");
    require(audit.complement_ledger.state == EXPECTED_COMPLEMENT_FNV,
            "partition complement FNV changed");
    return audit;
}

struct DirectStatus {
    i64 mass = 0;
    i128 margin = 0;
};

DirectStatus classify_one(Pair pair, u32 mask) {
    require(mask < (u32{1} << 30) && std::popcount(mask) == 8,
            "single literal control mask has wrong range/arity");
    const Geometry geometry = build_geometry(pair);
    i64 mass = 0;
    for (const auto& [failure, width] : geometry.safe_width_by_failure)
        if ((failure & ~mask) == 0) mass += width;
    return {mass, static_cast<i128>(63) * mass -
                  static_cast<i128>(4) * geometry.grid};
}

void print_hostile_controls() {
    const DirectStatus target_appended =
        classify_one({256, 663}, UINT32_C(0x21141284));
    const DirectStatus target_deleted =
        classify_one({256, 663}, UINT32_C(0x21948006));
    require(target_appended.margin == INT64_C(530687385396),
            "positive literal control changed");
    require(target_deleted.margin == -INT64_C(17044098303582),
            "inactive deleted-mask control changed");
    std::cout << "CONTROL 256,663 APPENDED 21141284 MASS "
              << target_appended.mass << " MARGIN "
              << decimal(target_appended.margin) << " ACTIVE 1\n"
              << "CONTROL 256,663 DELETED 21948006 MASS "
              << target_deleted.mass << " MARGIN "
              << decimal(target_deleted.margin) << " ACTIVE 0\n";
    for (Pair pair : {Pair{256, 287}, Pair{256, 558}, Pair{256, 575}}) {
        const DirectStatus hostile =
            classify_one(pair, UINT32_C(0x04871108));
        require(hostile.margin < 0,
                "hostile complement control became active");
        std::cout << "CONTROL " << pair.q << ',' << pair.r
                  << " APPENDED 04871108 MASS " << hostile.mass
                  << " MARGIN " << decimal(hostile.margin)
                  << " ACTIVE 0\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 6,
                "usage: surgery188-audit SURGERY_MASKS COMMON_CSV "
                "RESIDUAL_CSV OBLIGATIONS_CSV OUTPUT_DIR");
        const std::filesystem::path output_dir = argv[5];
        require(std::filesystem::is_directory(output_dir),
                "output directory does not exist");

        const std::vector<u32> deck = read_deck(argv[1]);
        const PairRows common = read_pair_rows(
            argv[2], COMMON_COUNT, EXPECTED_COMMON_FNV, "common rows");
        const PairRows residual = read_pair_rows(
            argv[3], RESIDUAL_COUNT, EXPECTED_RESIDUAL_FNV,
            "post-THM-4281 residual");
        const std::vector<Obligation> obligations =
            read_obligations(argv[4]);

        const MaskIndex index = build_mask_index(deck);
        const PartitionAudit partition = audit_partition(
            residual, common, output_dir / "residual_complement.csv");
        const BodyAudit body = audit_body_decomposition(
            deck, obligations,
            output_dir / "exposed_body_replacements.csv");
        const LiteralAudit literal = audit_common_literal(
            common, deck, index,
            output_dir / "common_literal_rows.csv",
            output_dir / "verified_common.csv");

        std::cout << "SURGERY188_DETACHED_LITERAL_AUDIT_V1\n"
                  << "INPUT DECK " << deck.size() << " FNV " << std::hex
                  << EXPECTED_DECK_FNV << " COMMON " << std::dec
                  << common.rows.size() << " FNV " << std::hex
                  << common.fnv << " RESIDUAL " << std::dec
                  << residual.rows.size() << " FNV " << std::hex
                  << residual.fnv << std::dec << '\n'
                  << "PARTITION COMMON " << partition.common
                  << " COMPLEMENT " << partition.complement
                  << " COMPLEMENT_FNV " << std::hex
                  << partition.complement_ledger.state << std::dec << '\n'
                  << "BODY_DECOMPOSITION TOTAL " << body.bodies
                  << " RETAINED_COVERED " << body.retained_covered
                  << " EXPOSED " << body.exposed
                  << " RETAINED_CHECKS " << body.retained_checks << '\n'
                  << "APPENDED_WITNESSES COUNT " << APPENDED.size()
                  << " INCIDENCES " << body.replacement_incidences
                  << " MIN_HITS " << body.replacement_minimum
                  << " MAX_HITS " << body.replacement_maximum
                  << " RESPONSE_FNV " << std::hex
                  << body.replacement_responses.state << std::dec << '\n'
                  << "COMMON_LITERAL ROWS " << common.rows.size()
                  << " MASKS " << deck.size() << " CHECKS "
                  << literal.checks << " EQUALITIES " << literal.equalities
                  << " GEOMETRY_FNV " << std::hex << literal.geometry.state
                  << " ACTIVITY_FNV " << literal.activity.state << std::dec
                  << '\n'
                  << "GLOBAL_MINIMUM ROW " << literal.global_min_pair.q << ','
                  << literal.global_min_pair.r << " MASK "
                  << hex8(literal.global_min_mask) << " MASS "
                  << literal.global_min_mass << " MARGIN "
                  << decimal(literal.global_min_margin) << '\n';
        print_hostile_controls();
        std::cout << "SCOPE EXPLICIT_APPEND_WITNESS_ONLY_NO_MINIMALITY_CLAIM\n"
                  << "VERDICT PASS EXACT_BODY_DECOMPOSITION_"
                     "COMMON_LITERAL_ACTIVITY_AND_RESIDUAL_PARTITION\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SURGERY188_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
