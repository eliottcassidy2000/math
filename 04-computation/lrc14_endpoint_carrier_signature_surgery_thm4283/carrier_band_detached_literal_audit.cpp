// Detached literal-wall and exhaustive-body audit of the THM-4283 carrier
// band.  This translation unit imports no project source and does not use the
// primitive cocycle, atom quotient, colex transform, or primary scanner.

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
constexpr std::size_t kBaseCount = 8996;
constexpr std::size_t kNestedCount = 8997;
constexpr u64 kInherited8951Fnv = UINT64_C(0x188f82ab9dd1695a);
constexpr u64 kAdditions45Fnv = UINT64_C(0xec083b65cc8c34e3);
constexpr u64 kBaseFnv = UINT64_C(0xfd899660f14b311c);
constexpr u64 kNestedFnv = UINT64_C(0x8e1860a25d0fcf87);
constexpr u32 kNestedRepair = UINT32_C(0x014c9084);
constexpr u64 kWitness9Fnv = UINT64_C(0x02b936529030e4bc);
constexpr u64 kRepaired9006Fnv = UINT64_C(0xfdc1c57ae4dc1bb6);
constexpr std::size_t kResidualCount = 23373;
constexpr u64 kResidualFnv = UINT64_C(0xc6ab0ae49ee32273);
constexpr u64 kBodyCount = UINT64_C(14307150);
constexpr i64 kFixedGrid = INT64_C(18241159416480);

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
    require(masks.size() == expected_count && ledger.state == expected_fnv,
            "mask-ledger identity changed");
    return masks;
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

std::vector<Pair> read_band(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open final residual");
    std::vector<Pair> residual;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank residual row");
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed residual row");
        Pair pair{std::stoi(line.substr(0, comma)),
                  std::stoi(line.substr(comma + 1))};
        require(pair.q > 0 && pair.q < pair.r &&
                    (residual.empty() || residual.back() < pair),
                "residual order changed");
        residual.push_back(pair);
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(residual.size() == kResidualCount && ledger.state == kResidualFnv,
            "final residual identity changed");
    std::vector<Pair> band;
    for (Pair pair : residual)
        if (pair.r >= 638) band.push_back(pair);
    std::sort(band.begin(), band.end(), [](Pair left, Pair right) {
        if (left.r != right.r) return left.r > right.r;
        return left.q < right.q;
    });
    require(band.size() == 64 && band.front() == Pair{220, 644} &&
                band.back() == Pair{520, 638},
            "638+ band changed");
    return band;
}

struct FailureRow {
    Pair pair;
    u32 body = 0;
    friend bool operator<(const FailureRow& left, const FailureRow& right) {
        return std::tie(left.pair.q, left.pair.r, left.body) <
               std::tie(right.pair.q, right.pair.r, right.body);
    }
    friend bool operator==(const FailureRow& left, const FailureRow& right) {
        return left.pair == right.pair && left.body == right.body;
    }
};

std::vector<FailureRow> read_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open primary failure ledger");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,body_hex",
            "failure-ledger header changed");
    std::vector<FailureRow> rows;
    while (std::getline(input, line)) {
        const std::size_t first = line.find(',');
        const std::size_t second = line.find(',', first + 1);
        require(first != std::string::npos && second != std::string::npos &&
                    line.find(',', second + 1) == std::string::npos,
                "malformed failure row");
        FailureRow row{{std::stoi(line.substr(0, first)),
                        std::stoi(line.substr(first + 1,
                                              second - first - 1))},
                       static_cast<u32>(
                           std::stoul(line.substr(second + 1), nullptr, 16))};
        require(std::popcount(row.body) == 9,
                "failure body is not rank nine");
        rows.push_back(row);
    }
    std::sort(rows.begin(), rows.end());
    return rows;
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
    i64 grid = 0;
    i64 pair_ticks = 0;
    u64 active_base = 0;
    u64 active_base_fnv = 0;
    bool repair_active = false;
    i128 repair_margin = 0;
    u64 equalities = 0;
    i128 minimum_active_margin = 0;
    i128 maximum_inactive_margin = 0;
    bool have_active = false;
    bool have_inactive = false;
    u64 body_checks = 0;
    std::vector<u32> base_failures;
    std::vector<u32> nested_failures;
};

PairAudit audit_pair(Pair pair, const std::vector<u32>& carrier,
                     const CarrierIndex& index) {
    PairAudit audit;
    audit.pair = pair;
    const Geometry geometry = build_geometry(pair);
    audit.cells = geometry.cells;
    audit.grid = geometry.grid;
    audit.pair_ticks = geometry.pair_ticks;
    const std::vector<i64> mass = classify(geometry, carrier, index);
    std::vector<u32> active_base;
    Fnv active_ledger;
    for (std::size_t mask_index = 0; mask_index < carrier.size(); ++mask_index) {
        const i128 margin = static_cast<i128>(63) * mass[mask_index] -
                            static_cast<i128>(4) * geometry.grid;
        audit.equalities += margin == 0;
        if (margin >= 0) {
            if (!audit.have_active || margin < audit.minimum_active_margin)
                audit.minimum_active_margin = margin;
            audit.have_active = true;
            if (mask_index < kBaseCount) {
                active_base.push_back(carrier[mask_index]);
                active_ledger.add(carrier[mask_index]);
            }
        } else {
            if (!audit.have_inactive || margin > audit.maximum_inactive_margin)
                audit.maximum_inactive_margin = margin;
            audit.have_inactive = true;
        }
        if (mask_index + 1 == carrier.size()) {
            audit.repair_active = margin >= 0;
            audit.repair_margin = margin;
        }
    }
    audit.active_base = active_base.size();
    audit.active_base_fnv = active_ledger.state;
    const u32 limit = u32{1} << 30;
    u64 bodies = 0;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++bodies;
        bool covered = false;
        for (u32 mask : active_base) {
            ++audit.body_checks;
            if ((mask & body) == 0) {
                covered = true;
                break;
            }
        }
        if (!covered) audit.base_failures.push_back(body);
    }
    require(bodies == kBodyCount, "nine-body universe count changed");
    for (u32 body : audit.base_failures)
        if (!audit.repair_active || (body & kNestedRepair) != 0)
            audit.nested_failures.push_back(body);
    return audit;
}

u64 body_fnv(const std::vector<u32>& bodies) {
    Fnv ledger;
    for (u32 body : bodies) ledger.add(body);
    return ledger.state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 9,
                "usage: carrier-detached BASE8951 ADDITIONS45 FINAL23373 "
                "PRIMARY_BASE_FAILURES PRIMARY_NESTED_FAILURES WITNESS9 "
                "OUT_CSV THREADS");
        std::vector<u32> carrier = read_masks(
            argv[1], 8951, kInherited8951Fnv);
        const std::vector<u32> additions = read_masks(
            argv[2], 45, kAdditions45Fnv);
        std::set<u32> distinct(carrier.begin(), carrier.end());
        for (u32 mask : additions) {
            require(distinct.insert(mask).second,
                    "addition overlaps inherited carrier");
            carrier.push_back(mask);
        }
        Fnv base_ledger;
        for (u32 mask : carrier) base_ledger.add(mask);
        require(carrier.size() == kBaseCount &&
                    base_ledger.state == kBaseFnv,
                "base8996 identity changed");
        require(distinct.insert(kNestedRepair).second,
                "nested repair overlaps base carrier");
        carrier.push_back(kNestedRepair);
        Fnv nested_ledger;
        for (u32 mask : carrier) nested_ledger.add(mask);
        require(carrier.size() == kNestedCount &&
                    nested_ledger.state == kNestedFnv,
                "nested8997 identity changed");
        const std::vector<Pair> band = read_band(argv[3]);
        const std::vector<FailureRow> expected_base = read_failures(argv[4]);
        const std::vector<FailureRow> expected_nested = read_failures(argv[5]);
        const std::vector<u32> witness9 =
            read_masks(argv[6], 9, kWitness9Fnv);
        require(expected_base.size() == 42 && expected_nested.size() == 40,
                "primary failure counts changed");
        const CarrierIndex carrier_index = build_index(carrier);

        const unsigned requested_threads =
            static_cast<unsigned>(std::stoul(argv[8]));
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

        std::vector<FailureRow> actual_base;
        std::vector<FailureRow> actual_nested;
        std::ofstream output(argv[7]);
        require(static_cast<bool>(output), "cannot create detached pair CSV");
        output << "q,r,grid,cells,pair_ticks,active_base,active_base_fnv,"
                  "repair_active,repair_margin,equalities,body_checks,"
                  "base_failures,base_failure_fnv,nested_failures,"
                  "nested_failure_fnv,min_active_margin,max_inactive_margin\n";
        Fnv pair_ledger;
        u64 total_checks = 0;
        u64 total_equalities = 0;
        for (const PairAudit& audit : audits) {
            for (u32 body : audit.base_failures)
                actual_base.push_back({audit.pair, body});
            for (u32 body : audit.nested_failures)
                actual_nested.push_back({audit.pair, body});
            if (audit.pair.r >= 639) {
                require(audit.nested_failures.empty(),
                        "nested carrier failed inside 639..644 band");
            } else if (audit.pair == Pair{256, 638}) {
                require(audit.base_failures.size() == 40 &&
                            audit.nested_failures == audit.base_failures,
                        "hostile (256,638) failure family changed");
            } else {
                require(audit.base_failures.empty() &&
                            audit.nested_failures.empty(),
                        "unexpected endpoint-638 failure");
            }
            require(audit.equalities == 0,
                    "detached carrier activity equality encountered");
            pair_ledger.add(audit.pair.q);
            pair_ledger.add(audit.pair.r);
            pair_ledger.add(audit.active_base);
            pair_ledger.add(audit.active_base_fnv);
            pair_ledger.add(audit.repair_active);
            add_i128(pair_ledger, audit.repair_margin);
            pair_ledger.add(audit.body_checks);
            pair_ledger.add(audit.base_failures.size());
            pair_ledger.add(body_fnv(audit.base_failures));
            pair_ledger.add(audit.nested_failures.size());
            pair_ledger.add(body_fnv(audit.nested_failures));
            total_checks += audit.body_checks;
            total_equalities += audit.equalities;
            output << audit.pair.q << ',' << audit.pair.r << ','
                   << audit.grid << ',' << audit.cells << ','
                   << audit.pair_ticks << ',' << audit.active_base << ','
                   << std::hex << audit.active_base_fnv << std::dec << ','
                   << audit.repair_active << ','
                   << decimal(audit.repair_margin) << ',' << audit.equalities
                   << ',' << audit.body_checks << ','
                   << audit.base_failures.size() << ',' << std::hex
                   << body_fnv(audit.base_failures) << std::dec << ','
                   << audit.nested_failures.size() << ',' << std::hex
                   << body_fnv(audit.nested_failures) << std::dec << ','
                   << decimal(audit.minimum_active_margin) << ','
                   << decimal(audit.maximum_inactive_margin) << '\n';
        }
        require(output.good(), "detached pair CSV write failed");
        std::sort(actual_base.begin(), actual_base.end());
        std::sort(actual_nested.begin(), actual_nested.end());
        require(actual_base == expected_base && actual_nested == expected_nested,
                "detached and primary failure ledgers disagree");

        const auto boundary_found = std::find_if(
            audits.begin(), audits.end(), [](const PairAudit& audit) {
                return audit.pair == Pair{256, 638};
            });
        require(boundary_found != audits.end() &&
                    boundary_found->nested_failures.size() == 40,
                "detached boundary row missing");
        const Geometry boundary_geometry = build_geometry({256, 638});
        i128 witness_minimum_margin = 0;
        bool witness_minimum_set = false;
        u64 witness_incidences = 0;
        for (u32 witness : witness9) {
            i64 mass = 0;
            for (const auto& [failure, width] :
                 boundary_geometry.safe_width_by_failure)
                if ((failure & ~witness) == 0) mass += width;
            const i128 margin = static_cast<i128>(63) * mass -
                                static_cast<i128>(4) *
                                    boundary_geometry.grid;
            require(margin > 0,
                    "boundary witness is not strict literal-active");
            if (!witness_minimum_set || margin < witness_minimum_margin)
                witness_minimum_margin = margin;
            witness_minimum_set = true;
        }
        for (u32 body : boundary_found->nested_failures) {
            bool covered = false;
            for (u32 witness : witness9)
                if ((body & witness) == 0) {
                    covered = true;
                    ++witness_incidences;
                }
            require(covered, "boundary witness family misses failed body");
        }
        Fnv repaired_ledger;
        for (u32 mask : carrier) repaired_ledger.add(mask);
        for (u32 witness : witness9) repaired_ledger.add(witness);
        require(repaired_ledger.state == kRepaired9006Fnv,
                "repaired9006 identity changed");

        std::cout << "THM4283_CARRIER_BAND_DETACHED_LITERAL_AUDIT_V1\n"
                  << "CARRIER BASE " << kBaseCount << " FNV " << std::hex
                  << kBaseFnv << " NESTED " << std::dec << kNestedCount
                  << " FNV " << std::hex << kNestedFnv << std::dec << '\n';
        for (const PairAudit& audit : audits)
            std::cout << "PAIR " << audit.pair.q << ',' << audit.pair.r
                      << " ACTIVE_BASE " << audit.active_base
                      << " ACTIVE_FNV " << std::hex
                      << audit.active_base_fnv << std::dec
                      << " REPAIR_ACTIVE " << audit.repair_active
                      << " REPAIR_MARGIN " << decimal(audit.repair_margin)
                      << " BODY_CHECKS " << audit.body_checks
                      << " BASE_FAILURES " << audit.base_failures.size()
                      << " BASE_FAILURE_FNV " << std::hex
                      << body_fnv(audit.base_failures) << std::dec
                      << " NESTED_FAILURES "
                      << audit.nested_failures.size()
                      << " NESTED_FAILURE_FNV " << std::hex
                      << body_fnv(audit.nested_failures) << std::dec << '\n';
        std::cout << "BAND PAIRS " << audits.size()
                  << " CLOSED_639_644 55 HOSTILE_638_ROWS 9"
                  << " TOTAL_BODY_CHECKS " << total_checks
                  << " EQUALITIES " << total_equalities
                  << " PAIR_LEDGER_FNV " << std::hex << pair_ledger.state
                  << std::dec << '\n'
                  << "FAILURES BASE " << actual_base.size()
                  << " NESTED " << actual_nested.size() << '\n'
                  << "BOUNDARY_WITNESS COUNT " << witness9.size()
                  << " FNV " << std::hex << kWitness9Fnv
                  << " REPAIRED9006_FNV " << kRepaired9006Fnv << std::dec
                  << " INCIDENCES " << witness_incidences
                  << " MIN_MARGIN " << decimal(witness_minimum_margin)
                  << " REPLAY_FAILURES 0\n"
                  << "SCOPE LITERAL_ACTIVITY_AND_EXHAUSTIVE_BODY_SCAN\n"
                  << "VERDICT PASS DETACHED_CARRIER_BAND_AND_SHARP_BOUNDARY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4283_CARRIER_DETACHED_ERROR "
                  << error.what() << '\n';
        return 1;
    }
}
