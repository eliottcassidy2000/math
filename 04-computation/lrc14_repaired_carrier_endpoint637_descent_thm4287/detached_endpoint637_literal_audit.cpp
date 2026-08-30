// Detached literal-wall and exhaustive-body audit of the THM-4287 repaired
// carrier boundary.  This translation unit imports no project source.  It
// independently reconstructs literal walls, exact activity, and every
// labelled rank-nine body for the complete endpoint-637 and endpoint-636
// layers of the post-THM-4283 residual.

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
#include <set>
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
constexpr std::size_t kInheritedCount = 8951;
constexpr std::size_t kAdditionCount = 45;
constexpr std::size_t kBaseCount = 8996;
constexpr std::size_t kNestedCount = 8997;
constexpr std::size_t kWitnessCount = 9;
constexpr std::size_t kRepairedCount = 9006;
constexpr u64 kInheritedFnv = UINT64_C(0x188f82ab9dd1695a);
constexpr u64 kAdditionsFnv = UINT64_C(0xec083b65cc8c34e3);
constexpr u64 kBaseFnv = UINT64_C(0xfd899660f14b311c);
constexpr u32 kPriorRepair = UINT32_C(0x014c9084);
constexpr u64 kNestedFnv = UINT64_C(0x8e1860a25d0fcf87);
constexpr u64 kWitnessFnv = UINT64_C(0x02b936529030e4bc);
constexpr u64 kRepairSuffixFnv = UINT64_C(0x0aeb35a950785473);
constexpr u64 kRepairedFnv = UINT64_C(0xfdc1c57ae4dc1bb6);
constexpr std::size_t kResidualCount = 22682;
constexpr u64 kResidualFnv = UINT64_C(0xf7563445f15efebf);
constexpr std::size_t kBandCount = 9;
constexpr u64 kBandFnv = UINT64_C(0x53c3da1b8069ca12);
constexpr u64 kLayer637Fnv = UINT64_C(0xf4b48b16a28826ad);
constexpr u64 kLayer636Fnv = UINT64_C(0xb3a1869a95dec7aa);
constexpr u64 kBodyCount = UINT64_C(14307150);
constexpr i64 kFixedGrid = INT64_C(18241159416480);
constexpr u64 kEmptyFnv = UINT64_C(0xcbf29ce484222325);
constexpr u64 kFailureLedgerFnv = UINT64_C(0x603a550bc39265fa);
constexpr u64 kTotalBodyChecks = UINT64_C(4443963867);
constexpr u64 kPairLedgerFnv = UINT64_C(0xeb9f876da9477a31);

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
            "bad carrier mask token");
    const u32 mask = static_cast<u32>(wide);
    require(std::popcount(mask) == 8, "carrier mask is not rank eight");
    return mask;
}

std::vector<u32> read_masks(const std::filesystem::path& path,
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
    friend bool operator<(const Pair& left, const Pair& right) {
        return std::tie(left.q, left.r) < std::tie(right.q, right.r);
    }
    friend bool operator==(const Pair& left, const Pair& right) {
        return left.q == right.q && left.r == right.r;
    }
};

constexpr std::array<Pair, kBandCount> kExpectedBand = {{
    {100, 637}, {294, 637}, {520, 637},
    {100, 636}, {256, 636}, {294, 636},
    {338, 636}, {372, 636}, {384, 636}}};

u64 pair_fnv(const std::vector<Pair>& pairs) {
    Fnv ledger;
    for (Pair pair : pairs) {
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    return ledger.state;
}

std::vector<Pair> read_boundary_band(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open current residual");
    std::vector<Pair> residual;
    std::vector<Pair> band;
    Fnv residual_ledger;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank current-residual row");
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed current-residual row");
        std::size_t used_q = 0;
        std::size_t used_r = 0;
        Pair pair;
        pair.q = std::stoi(line.substr(0, comma), &used_q);
        pair.r = std::stoi(line.substr(comma + 1), &used_r);
        require(used_q == comma && used_r == line.size() - comma - 1 &&
                    pair.q > 0 && pair.q < pair.r &&
                    (residual.empty() || residual.back() < pair),
                "invalid or out-of-order current-residual pair");
        residual.push_back(pair);
        residual_ledger.add(pair.q);
        residual_ledger.add(pair.r);
        if (pair.r >= 636) band.push_back(pair);
    }
    require(input.eof(), "current-residual read failed");
    require(residual.size() == kResidualCount &&
                residual_ledger.state == kResidualFnv,
            "post-THM4283 residual identity changed");
    std::sort(band.begin(), band.end(), [](Pair left, Pair right) {
        if (left.r != right.r) return left.r > right.r;
        return left.q < right.q;
    });
    require(band.size() == kExpectedBand.size() &&
                pair_fnv(band) == kBandFnv,
            "endpoint-637/636 boundary identity changed");
    for (std::size_t index = 0; index < band.size(); ++index)
        require(band[index] == kExpectedBand[index],
                "endpoint-637/636 boundary row changed");
    return band;
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

struct CarrierIndex {
    std::vector<std::vector<u64>> contains;
    std::vector<u64> all;
};

CarrierIndex build_index(const std::vector<u32>& carrier) {
    const std::size_t words = (carrier.size() + 63) / 64;
    CarrierIndex index;
    index.contains.assign(30, std::vector<u64>(words, 0));
    index.all.assign(words, 0);
    for (std::size_t mask_index = 0; mask_index < carrier.size(); ++mask_index) {
        index.all[mask_index / 64] |= UINT64_C(1) << (mask_index % 64);
        for (unsigned vertex = 0; vertex < 30; ++vertex)
            if ((carrier[mask_index] >> vertex) & 1)
                index.contains[vertex][mask_index / 64] |=
                    UINT64_C(1) << (mask_index % 64);
    }
    return index;
}

std::vector<i64> classify(const Geometry& geometry,
                          const std::vector<u32>& carrier,
                          const CarrierIndex& index) {
    std::vector<i64> mass(carrier.size(), 0);
    for (const auto& [failure, width] : geometry.safe_width_by_failure) {
        if (std::popcount(failure) > 8) continue;
        std::vector<u64> candidates = index.all;
        u32 remaining = failure;
        while (remaining != 0) {
            const unsigned vertex = std::countr_zero(remaining);
            remaining &= remaining - 1;
            for (std::size_t word = 0; word < candidates.size(); ++word)
                candidates[word] &= index.contains[vertex][word];
        }
        for (std::size_t word = 0; word < candidates.size(); ++word) {
            u64 bits = candidates[word];
            while (bits != 0) {
                const unsigned offset = std::countr_zero(bits);
                const std::size_t mask_index = 64 * word + offset;
                require(mask_index < mass.size(), "carrier bitset tail escaped");
                mass[mask_index] += width;
                bits &= bits - 1;
            }
        }
    }
    return mass;
}

u32 next_combination(u32 value) {
    const u32 low = value & (u32{0} - value);
    require(low != 0, "zero low bit in combination successor");
    const u32 ripple = value + low;
    return ripple | (((ripple ^ value) >> 2) / low);
}

struct PairAudit {
    Pair pair;
    u64 cells = 0;
    u64 failure_classes = 0;
    i64 grid = 0;
    i64 pair_ticks = 0;
    u64 active_carrier = 0;
    u64 active_carrier_fnv = 0;
    u64 active_base = 0;
    u64 active_suffix = 0;
    u64 active_suffix_fnv = 0;
    i128 minimum_suffix_margin = 0;
    i128 maximum_suffix_margin = 0;
    bool have_suffix = false;
    u64 equalities = 0;
    i128 minimum_active_margin = 0;
    i128 maximum_inactive_margin = 0;
    bool have_active = false;
    bool have_inactive = false;
    u64 bodies = 0;
    u64 body_checks = 0;
    std::vector<u32> failures;
    std::vector<u32> base_failures;
};

PairAudit audit_pair(Pair pair, const std::vector<u32>& carrier,
                     const CarrierIndex& index) {
    PairAudit audit;
    audit.pair = pair;
    const Geometry geometry = build_geometry(pair);
    audit.cells = geometry.cells;
    audit.failure_classes = geometry.safe_width_by_failure.size();
    audit.grid = geometry.grid;
    audit.pair_ticks = geometry.pair_ticks;
    const std::vector<i64> mass = classify(geometry, carrier, index);

    std::vector<u32> active;
    active.reserve(carrier.size());
    Fnv active_ledger;
    Fnv suffix_ledger;
    for (std::size_t mask_index = 0; mask_index < carrier.size(); ++mask_index) {
        const i128 margin = static_cast<i128>(63) * mass[mask_index] -
                            static_cast<i128>(4) * geometry.grid;
        audit.equalities += margin == 0;
        if (margin >= 0) {
            if (!audit.have_active || margin < audit.minimum_active_margin)
                audit.minimum_active_margin = margin;
            audit.have_active = true;
            active.push_back(carrier[mask_index]);
            active_ledger.add(carrier[mask_index]);
            if (mask_index < kBaseCount) ++audit.active_base;
        } else {
            if (!audit.have_inactive || margin > audit.maximum_inactive_margin)
                audit.maximum_inactive_margin = margin;
            audit.have_inactive = true;
        }
        if (mask_index >= kBaseCount) {
            if (!audit.have_suffix || margin < audit.minimum_suffix_margin)
                audit.minimum_suffix_margin = margin;
            if (!audit.have_suffix || margin > audit.maximum_suffix_margin)
                audit.maximum_suffix_margin = margin;
            audit.have_suffix = true;
            if (margin >= 0) {
                ++audit.active_suffix;
                suffix_ledger.add(carrier[mask_index]);
            }
        }
    }
    require(audit.have_active && audit.have_inactive && audit.have_suffix,
            "degenerate detached activity partition");
    audit.active_carrier = active.size();
    audit.active_carrier_fnv = active_ledger.state;
    audit.active_suffix_fnv = suffix_ledger.state;

    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++audit.bodies;
        bool covered = false;
        bool base_covered = false;
        for (std::size_t active_index = 0; active_index < active.size();
             ++active_index) {
            ++audit.body_checks;
            if ((active[active_index] & body) == 0) {
                covered = true;
                base_covered = active_index < audit.active_base;
                break;
            }
        }
        if (!covered) audit.failures.push_back(body);
        if (!base_covered) audit.base_failures.push_back(body);
    }
    require(audit.bodies == kBodyCount, "nine-body universe count changed");
    return audit;
}

u64 body_fnv(const std::vector<u32>& bodies) {
    Fnv ledger;
    for (u32 body : bodies) ledger.add(body);
    return ledger.state;
}

struct ExpectedFailures {
    Pair pair;
    std::size_t count = 0;
    u64 fnv = kEmptyFnv;
};

constexpr std::array<ExpectedFailures, kBandCount> kExpectedFailures = {{
    {{100, 637}, 0, kEmptyFnv},
    {{294, 637}, 0, kEmptyFnv},
    {{520, 637}, 0, kEmptyFnv},
    {{100, 636}, 64, UINT64_C(0xd611500ea833ff83)},
    {{256, 636}, 37, UINT64_C(0xee7792a8a2fd51c9)},
    {{294, 636}, 0, kEmptyFnv},
    {{338, 636}, 0, kEmptyFnv},
    {{372, 636}, 0, kEmptyFnv},
    {{384, 636}, 0, kEmptyFnv}}};

const ExpectedFailures& expected_failures(Pair pair) {
    const auto found = std::find_if(
        kExpectedFailures.begin(), kExpectedFailures.end(),
        [&](const ExpectedFailures& expected) { return expected.pair == pair; });
    require(found != kExpectedFailures.end(), "unexpected boundary pair");
    return *found;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 7,
                "usage: detached-endpoint637 BASE8951 ADDITIONS45 WITNESS9 "
                "CURRENT22682 FAILURE_LEDGER THREADS");

        std::vector<u32> carrier =
            read_masks(argv[1], kInheritedCount, kInheritedFnv);
        const std::vector<u32> additions =
            read_masks(argv[2], kAdditionCount, kAdditionsFnv);
        const std::vector<u32> witness =
            read_masks(argv[3], kWitnessCount, kWitnessFnv);
        std::set<u32> distinct(carrier.begin(), carrier.end());
        for (u32 mask : additions) {
            require(distinct.insert(mask).second,
                    "addition overlaps inherited carrier");
            carrier.push_back(mask);
        }
        require(carrier.size() == kBaseCount &&
                    mask_fnv(carrier) == kBaseFnv,
                "base8996 identity changed");
        require(std::popcount(kPriorRepair) == 8 &&
                    distinct.insert(kPriorRepair).second,
                "prior repair is invalid or overlaps base8996");
        carrier.push_back(kPriorRepair);
        require(carrier.size() == kNestedCount &&
                    mask_fnv(carrier) == kNestedFnv,
                "nested8997 identity changed");
        std::vector<u32> repair_suffix = {kPriorRepair};
        for (u32 mask : witness) {
            require(distinct.insert(mask).second,
                    "endpoint-638 witness overlaps nested8997");
            repair_suffix.push_back(mask);
            carrier.push_back(mask);
        }
        require(mask_fnv(repair_suffix) == kRepairSuffixFnv,
                "ten-mask repair suffix identity changed");
        require(carrier.size() == kRepairedCount &&
                    mask_fnv(carrier) == kRepairedFnv,
                "repaired9006 identity changed");

        const std::vector<Pair> band = read_boundary_band(argv[4]);
        const CarrierIndex carrier_index = build_index(carrier);
        const unsigned requested_threads =
            static_cast<unsigned>(std::stoul(argv[6]));
        require(requested_threads >= 1 && requested_threads <= 8,
                "thread count outside 1..8");

        std::vector<PairAudit> audits(band.size());
        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        for (unsigned thread = 0; thread < requested_threads; ++thread) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t index = next.fetch_add(1);
                    if (index >= band.size()) break;
                    audits[index] = audit_pair(
                        band[index], carrier, carrier_index);
                }
            });
        }
        for (std::thread& worker : workers) worker.join();

        std::ofstream failures_out(argv[5]);
        require(static_cast<bool>(failures_out),
                "cannot create detached failure ledger");
        failures_out << "q,r,body_hex\n";

        std::cout << "THM4287_ENDPOINT637_636_DETACHED_LITERAL_AUDIT_V1\n"
                  << "CARRIER INHERITED " << kInheritedCount << " FNV "
                  << std::hex << kInheritedFnv << " ADDITIONS " << std::dec
                  << kAdditionCount << " FNV " << std::hex << kAdditionsFnv
                  << " BASE8996_FNV " << kBaseFnv << " PRIOR_REPAIR "
                  << std::setw(8) << std::setfill('0') << kPriorRepair
                  << std::setfill(' ') << " NESTED8997_FNV " << kNestedFnv
                  << " WITNESS9_FNV " << kWitnessFnv
                  << " REPAIRED9006_FNV " << kRepairedFnv << std::dec
                  << '\n'
                  << "CURRENT_RESIDUAL " << kResidualCount << " FNV "
                  << std::hex << kResidualFnv << std::dec << " BAND "
                  << band.size() << " FNV " << std::hex << pair_fnv(band)
                  << std::dec << " BODY_UNIVERSE " << kBodyCount << '\n';

        Fnv pair_ledger;
        Fnv complete_failure_ledger;
        u64 total_bodies = 0;
        u64 total_checks = 0;
        u64 total_equalities = 0;
        u64 total_failures = 0;
        std::size_t layer_begin = 0;
        while (layer_begin < audits.size()) {
            const int endpoint = audits[layer_begin].pair.r;
            std::size_t layer_end = layer_begin;
            std::vector<Pair> layer_pairs;
            u64 layer_failures = 0;
            while (layer_end < audits.size() &&
                   audits[layer_end].pair.r == endpoint) {
                const PairAudit& audit = audits[layer_end];
                const ExpectedFailures& expected =
                    expected_failures(audit.pair);
                require(audit.equalities == 0,
                        "detached carrier activity equality encountered");
                require(audit.failures.size() == expected.count &&
                            body_fnv(audit.failures) == expected.fnv,
                        "detached boundary failure family changed");
                if (audit.pair.r == 637)
                    require(audit.base_failures.empty(),
                            "base8996 does not close endpoint 637");
                for (u32 body : audit.failures) {
                    failures_out << audit.pair.q << ',' << audit.pair.r << ','
                                 << std::hex << std::setw(8)
                                 << std::setfill('0') << body << std::dec
                                 << std::setfill(' ') << '\n';
                    complete_failure_ledger.add(audit.pair.q);
                    complete_failure_ledger.add(audit.pair.r);
                    complete_failure_ledger.add(body);
                }
                layer_pairs.push_back(audit.pair);
                layer_failures += audit.failures.size();
                total_bodies += audit.bodies;
                total_checks += audit.body_checks;
                total_equalities += audit.equalities;
                total_failures += audit.failures.size();
                pair_ledger.add(audit.pair.q);
                pair_ledger.add(audit.pair.r);
                pair_ledger.add(audit.grid);
                pair_ledger.add(audit.cells);
                pair_ledger.add(audit.failure_classes);
                pair_ledger.add(audit.pair_ticks);
                pair_ledger.add(audit.active_carrier);
                pair_ledger.add(audit.active_carrier_fnv);
                pair_ledger.add(audit.active_suffix);
                pair_ledger.add(audit.active_suffix_fnv);
                add_i128(pair_ledger, audit.minimum_suffix_margin);
                add_i128(pair_ledger, audit.maximum_suffix_margin);
                pair_ledger.add(audit.equalities);
                add_i128(pair_ledger, audit.minimum_active_margin);
                add_i128(pair_ledger, audit.maximum_inactive_margin);
                pair_ledger.add(audit.bodies);
                pair_ledger.add(audit.body_checks);
                pair_ledger.add(audit.failures.size());
                pair_ledger.add(body_fnv(audit.failures));

                std::cout << "PAIR " << audit.pair.q << ',' << audit.pair.r
                          << " GRID " << audit.grid << " CELLS "
                          << audit.cells << " FAILURE_CLASSES "
                          << audit.failure_classes << " PAIR_TICKS "
                          << audit.pair_ticks << " ACTIVE_CARRIER "
                          << audit.active_carrier << " ACTIVE_FNV "
                          << std::hex << audit.active_carrier_fnv << std::dec
                          << " ACTIVE_BASE8996 " << audit.active_base
                          << " BASE8996_FAILURES "
                          << audit.base_failures.size()
                          << " BASE8996_FAILURE_FNV " << std::hex
                          << body_fnv(audit.base_failures) << std::dec
                          << " ACTIVE_SUFFIX " << audit.active_suffix
                          << " ACTIVE_SUFFIX_FNV " << std::hex
                          << audit.active_suffix_fnv << std::dec
                          << " SUFFIX_MARGIN_RANGE "
                          << decimal(audit.minimum_suffix_margin) << ".."
                          << decimal(audit.maximum_suffix_margin)
                          << " EQUALITIES " << audit.equalities
                          << " ACTIVE_MARGIN_MIN "
                          << decimal(audit.minimum_active_margin)
                          << " INACTIVE_MARGIN_MAX "
                          << decimal(audit.maximum_inactive_margin)
                          << " BODIES " << audit.bodies << " BODY_CHECKS "
                          << audit.body_checks << " FAILURES "
                          << audit.failures.size() << " FAILURE_FNV "
                          << std::hex << body_fnv(audit.failures) << std::dec
                          << '\n';
                ++layer_end;
            }
            const u64 layer_hash = pair_fnv(layer_pairs);
            if (endpoint == 637)
                require(layer_pairs.size() == 3 &&
                            layer_hash == kLayer637Fnv &&
                            layer_failures == 0,
                        "endpoint-637 detached layer changed");
            else if (endpoint == 636)
                require(layer_pairs.size() == 6 &&
                            layer_hash == kLayer636Fnv &&
                            layer_failures == 101,
                        "endpoint-636 detached layer changed");
            else
                fail("unexpected detached endpoint layer");
            std::cout << "LAYER " << endpoint << " ROWS "
                      << layer_pairs.size() << " FNV " << std::hex
                      << layer_hash << std::dec << " FAILURES "
                      << layer_failures << '\n';
            layer_begin = layer_end;
        }
        require(failures_out.good(), "detached failure-ledger write failed");
        require(total_bodies == kBandCount * kBodyCount &&
                    total_checks == kTotalBodyChecks &&
                    total_equalities == 0 && total_failures == 101 &&
                    complete_failure_ledger.state == kFailureLedgerFnv &&
                    pair_ledger.state == kPairLedgerFnv,
                "detached aggregate boundary changed");

        std::cout << "FAILURE_LEDGER ROWS " << total_failures << " FNV "
                  << std::hex << complete_failure_ledger.state << std::dec
                  << " TOTAL_BODIES " << total_bodies
                  << " TOTAL_BODY_CHECKS " << total_checks
                  << " EQUALITIES " << total_equalities
                  << " PAIR_LEDGER_FNV " << std::hex << pair_ledger.state
                  << std::dec << '\n'
                  << "SCOPE LITERAL_ACTIVITY_AND_EXHAUSTIVE_BODY_SCAN_"
                     "FIXED_POOL_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS DETACHED_ENDPOINT637_CLOSURE_AND_"
                     "ENDPOINT636_BOUNDARY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4287_DETACHED_LITERAL_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
