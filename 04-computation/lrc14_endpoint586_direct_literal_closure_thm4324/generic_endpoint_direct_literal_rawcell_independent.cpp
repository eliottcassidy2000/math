// Independent raw-cell replay for a runtime-selected endpoint literal audit.
// Unlike the primary implementation, this retains every pair-safe open wall
// cell separately and parallelizes over bodies; no failure-class aggregation
// is used when computing either the truncated or complete mass.

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
#include <vector>

namespace {
using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr std::array<int, 30> kPool{
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
int gEndpoint = 0;
u64 gExpectedFailureFnv = 0;
u64 gExpectedBodyCount = 0;
constexpr u64 kFnvBasis = UINT64_C(0xcbf29ce484222325);

void demand(bool truth, const std::string& message) {
    if (!truth) throw std::runtime_error(message);
}

struct Fnv64 {
    u64 value = kFnvBasis;
    void word(u64 input) {
        for (unsigned shift = 0; shift != 64; shift += 8) {
            value ^= (input >> shift) & UINT64_C(255);
            value *= UINT64_C(1099511628211);
        }
    }
};

std::string to_decimal(i128 number) {
    if (number == 0) return "0";
    const bool minus = number < 0;
    u128 value = static_cast<u128>(number);
    if (minus) value = u128{0} - value;
    std::string text;
    do {
        text.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    } while (value != 0);
    if (minus) text.push_back('-');
    std::reverse(text.begin(), text.end());
    return text;
}

std::string to_hex(u32 mask) {
    std::ostringstream text;
    text << std::hex << std::setfill('0') << std::setw(8) << mask;
    return text.str();
}

i64 lcm_exact(i64 a, i64 b) {
    const i128 candidate = static_cast<i128>(a / std::gcd(a, b)) * b;
    demand(candidate <= std::numeric_limits<i64>::max(), "grid overflow");
    return static_cast<i64>(candidate);
}

bool is_safe(int velocity, i64 modulus, i64 a, i64 b) {
    // The test point is (a+b)/(2*modulus).  Multiplying by 14*velocity
    // avoids floating point and reduces modulo one turn.
    i128 doubled_residue = static_cast<i128>(velocity) * (a + i128{b});
    doubled_residue %= i128{2} * modulus;
    if (doubled_residue < 0) doubled_residue += i128{2} * modulus;
    return i128{7} * doubled_residue >= modulus &&
           i128{7} * doubled_residue <= i128{13} * modulus;
}

struct Cell {
    u32 failed;
    i64 width;
    bool low;
};

struct WallModel {
    i64 grid = 1;
    u64 open_cells = 0;
    std::vector<Cell> pair_safe;
};

WallModel make_wall_model(int q) {
    WallModel model;
    for (int velocity : kPool)
        model.grid = lcm_exact(model.grid, i64{14} * velocity);
    model.grid = lcm_exact(model.grid, i64{14} * q);
    model.grid = lcm_exact(model.grid, i64{14} * gEndpoint);
    std::vector<i64> cut{0, model.grid};
    auto cut_for = [&](int velocity) {
        demand(model.grid % (i64{14} * velocity) == 0,
               "nonintegral tooth step");
        const i64 step = model.grid / (i64{14} * velocity);
        for (int period = 0; period < velocity; ++period) {
            cut.push_back((i64{14} * period + 1) * step);
            cut.push_back((i64{14} * period + 13) * step);
        }
    };
    for (int velocity : kPool) cut_for(velocity);
    cut_for(q);
    cut_for(gEndpoint);
    std::sort(cut.begin(), cut.end());
    cut.erase(std::unique(cut.begin(), cut.end()), cut.end());
    model.open_cells = cut.size() - 1;
    for (std::size_t i = 1; i < cut.size(); ++i) {
        const i64 a = cut[i - 1], b = cut[i];
        if (!is_safe(q, model.grid, a, b) ||
            !is_safe(gEndpoint, model.grid, a, b)) continue;
        u32 failed = 0;
        for (unsigned label = 0; label < kPool.size(); ++label)
            if (!is_safe(kPool[label], model.grid, a, b))
                failed |= u32{1} << label;
        model.pair_safe.push_back(
            Cell{failed, b - a, std::popcount(failed) <= 9});
    }
    return model;
}

struct FailureInput {
    std::map<int, std::vector<u32>> row;
    u64 count = 0;
    u64 fnv = 0;
};

FailureInput load_failures(const std::filesystem::path& name) {
    std::ifstream file(name);
    demand(bool(file), "cannot open failures");
    std::string record;
    demand(std::getline(file, record) && record == "q,r,body_hex",
           "unexpected header");
    FailureInput input;
    std::map<int, std::set<u32>> unique;
    Fnv64 fingerprint;
    int last_q = -1;
    u32 last_body = 0;
    while (std::getline(file, record)) {
        if (record.empty()) continue;
        const std::size_t first = record.find(',');
        const std::size_t second = record.find(',', first + 1);
        demand(first != std::string::npos && second != std::string::npos &&
                   record.find(',', second + 1) == std::string::npos,
               "bad CSV row");
        const int q = std::stoi(record.substr(0, first));
        const int r = std::stoi(record.substr(first + 1, second - first - 1));
        const u64 wide = std::stoull(record.substr(second + 1), nullptr, 16);
        demand(r == gEndpoint && wide < (u64{1} << 30),
               "failure escaped universe");
        const u32 body = static_cast<u32>(wide);
        demand(std::popcount(body) == 9 && unique[q].insert(body).second,
               "invalid or duplicate body");
        demand(q > last_q || (q == last_q && body > last_body),
               "noncanonical failure ordering");
        last_q = q; last_body = body;
        input.row[q].push_back(body);
        fingerprint.word(q); fingerprint.word(r); fingerprint.word(body);
        ++input.count;
    }
    demand(file.eof() && input.count == gExpectedBodyCount &&
               fingerprint.value == gExpectedFailureFnv,
           "failure census/fingerprint changed");
    input.fnv = fingerprint.value;
    return input;
}

struct Mass {
    i64 low = 0;
    i64 full = 0;
};

std::vector<Mass> evaluate(const WallModel& model,
                           const std::vector<u32>& bodies,
                           unsigned threads) {
    std::vector<Mass> answer(bodies.size());
    std::atomic<std::size_t> cursor{0};
    auto work = [&]() {
        while (true) {
            const std::size_t i = cursor.fetch_add(1, std::memory_order_relaxed);
            if (i >= bodies.size()) return;
            Mass mass;
            for (const Cell& cell : model.pair_safe) {
                if ((cell.failed & bodies[i]) != 0) continue;
                mass.full += cell.width;
                if (cell.low) mass.low += cell.width;
            }
            demand(mass.low <= mass.full, "truncation exceeded full mass");
            answer[i] = mass;
        }
    };
    std::vector<std::thread> pool;
    for (unsigned i = 0; i < threads; ++i) pool.emplace_back(work);
    for (auto& thread : pool) thread.join();
    return answer;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        demand(argc == 7,
               "usage: independent ENDPOINT COUNT FAILURE_FNV_HEX FAILURES DETAIL_CSV THREADS");
        gEndpoint = std::stoi(argv[1]);
        gExpectedBodyCount = std::stoull(argv[2]);
        gExpectedFailureFnv = std::stoull(argv[3], nullptr, 16);
        demand(gEndpoint > 0 && gExpectedBodyCount > 0,
               "invalid endpoint/count");
        const unsigned threads = std::stoul(argv[6]);
        demand(threads >= 1 && threads <= 64, "thread count escaped 1..64");
        const FailureInput failures = load_failures(argv[4]);
        std::ofstream detail(argv[5], std::ios::binary);
        demand(bool(detail), "cannot create independent detail");
        detail << "q,r,ordinal,body_hex,truncated_mass,full_mass,truncated_scaled_ticks,full_scaled_ticks\n";
        u64 low_positive = 0, full_positive = 0, equality = 0;
        Fnv64 global;
        std::cout << "LRC14_GENERIC_ENDPOINT_DIRECT_LITERAL_RAWCELL_INDEPENDENT_V1\n"
                  << "ENDPOINT " << gEndpoint << '\n'
                  << "ROWS " << failures.row.size() << " BODIES "
                  << failures.count << " FAILURE_FNV " << std::hex
                  << failures.fnv << std::dec << " THREADS " << threads << '\n';
        for (const auto& [q, bodies] : failures.row) {
            const WallModel model = make_wall_model(q);
            const std::vector<Mass> masses = evaluate(model, bodies, threads);
            i128 min_low = 0, min_full = 0;
            u32 min_low_body = 0, min_full_body = 0;
            u64 row_low = 0, row_full = 0, row_equal = 0;
            Fnv64 ledger;
            for (std::size_t i = 0; i < bodies.size(); ++i) {
                const i128 low_ticks = i128{63} * masses[i].low -
                                       i128{4} * model.grid;
                const i128 full_ticks = i128{63} * masses[i].full -
                                        i128{4} * model.grid;
                detail << q << ',' << gEndpoint << ',' << i << ',' << to_hex(bodies[i])
                       << ',' << masses[i].low << ',' << masses[i].full << ','
                       << to_decimal(low_ticks) << ',' << to_decimal(full_ticks)
                       << '\n';
                if (i == 0 || low_ticks < min_low ||
                    (low_ticks == min_low && bodies[i] < min_low_body)) {
                    min_low = low_ticks; min_low_body = bodies[i];
                }
                if (i == 0 || full_ticks < min_full ||
                    (full_ticks == min_full && bodies[i] < min_full_body)) {
                    min_full = full_ticks; min_full_body = bodies[i];
                }
                row_low += low_ticks > 0;
                row_full += full_ticks > 0;
                row_equal += masses[i].low == masses[i].full;
                ledger.word(q); ledger.word(gEndpoint); ledger.word(i);
                ledger.word(bodies[i]); ledger.word(masses[i].low);
                ledger.word(masses[i].full);
                const u128 low_bits = static_cast<u128>(low_ticks);
                const u128 full_bits = static_cast<u128>(full_ticks);
                ledger.word(static_cast<u64>(low_bits));
                ledger.word(static_cast<u64>(low_bits >> 64));
                ledger.word(static_cast<u64>(full_bits));
                ledger.word(static_cast<u64>(full_bits >> 64));
            }
            low_positive += row_low; full_positive += row_full;
            equality += row_equal;
            global.word(q); global.word(bodies.size()); global.word(model.grid);
            global.word(model.open_cells); global.word(model.pair_safe.size());
            global.word(row_low); global.word(row_full); global.word(row_equal);
            global.word(ledger.value);
            std::cout << "ROW " << q << ',' << gEndpoint << " BODIES "
                      << bodies.size() << " GRID " << model.grid
                      << " CELLS " << model.open_cells << " PAIR_SAFE_CELLS "
                      << model.pair_safe.size() << " LOW_POSITIVE " << row_low
                      << " FULL_POSITIVE " << row_full << " EQUALITY "
                      << row_equal << " LOW_MIN_TICKS " << to_decimal(min_low)
                      << " LOW_MIN_BODY " << to_hex(min_low_body)
                      << " FULL_MIN_TICKS " << to_decimal(min_full)
                      << " FULL_MIN_BODY " << to_hex(min_full_body)
                      << " DETAIL_FNV " << std::hex << ledger.value << std::dec
                      << '\n';
        }
        demand(detail.good(), "independent detail write failed");
        std::cout << "SUMMARY BODIES " << failures.count
                  << " LOW_POSITIVE " << low_positive
                  << " FULL_POSITIVE " << full_positive
                  << " EQUALITY " << equality << " RAW_SUMMARY_FNV "
                  << std::hex << global.value << std::dec << '\n'
                  << "REPRESENTATION RAW_PAIR_SAFE_WALL_CELLS_NO_CLASS_AGGREGATION\n"
                  << "SCOPE FINITE_EXACT_FIXED_POOL_SINGLE_ENDPOINT_FAILURE_BODIES_"
                     "NO_CARRIER_EXCHANGE_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT "
                  << (low_positive == failures.count ? "PASS" : "HOSTILE_FAIL")
                  << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "GENERIC_ENDPOINT_RAWCELL_INDEPENDENT_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
