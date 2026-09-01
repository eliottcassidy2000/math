// Import-free endpoint-587 replay of the serialized THM-4318 carrier.
// Activity is reconstructed directly from unaggregated pair-safe wall cells.

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
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <utility>
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
constexpr u64 kFnvBasis = UINT64_C(0xcbf29ce484222325);
constexpr u64 kFnvPrime = UINT64_C(0x100000001b3);
constexpr u64 kBodyCount = UINT64_C(14307150);

[[noreturn]] void die(const std::string& message) {
    throw std::runtime_error(message);
}
void check(bool condition, const std::string& message) {
    if (!condition) die(message);
}

struct Fnv {
    u64 value = kFnvBasis;
    void add(u64 word) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            value ^= (word >> (8 * byte)) & UINT64_C(0xff);
            value *= kFnvPrime;
        }
    }
};

u32 parse_hex(const std::string& token) {
    const u64 wide = std::stoull(token, nullptr, 16);
    check(wide < (UINT64_C(1) << 30), "mask escaped label universe");
    return static_cast<u32>(wide);
}

std::string hex8(u32 value) {
    std::ostringstream out;
    out << std::hex << std::setw(8) << std::setfill('0') << value;
    return out.str();
}

i64 checked_lcm(i64 left, i64 right) {
    const i128 wide = static_cast<i128>(left / std::gcd(left, right)) * right;
    check(wide <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(wide);
}

bool safe_midpoint(int speed, i64 grid, i64 left, i64 right) {
    i128 twice_phase = static_cast<i128>(speed) *
                       (static_cast<i128>(left) + right);
    twice_phase %= static_cast<i128>(2) * grid;
    if (twice_phase < 0) twice_phase += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * twice_phase >= grid &&
           static_cast<i128>(7) * twice_phase <= static_cast<i128>(13) * grid;
}

struct Geometry {
    i64 grid = 0;
    std::vector<std::pair<u32, i64>> low_cells;
};

Geometry make_geometry(int q, int r) {
    Geometry answer;
    answer.grid = 1;
    for (int speed : kPool)
        answer.grid = checked_lcm(answer.grid, 14LL * speed);
    answer.grid = checked_lcm(answer.grid, 14LL * q);
    answer.grid = checked_lcm(answer.grid, 14LL * r);

    std::vector<i64> walls{0, answer.grid};
    const auto append_walls = [&](int speed) {
        const i64 divisor = 14LL * speed;
        check(answer.grid % divisor == 0, "nonintegral wall unit");
        const i64 unit = answer.grid / divisor;
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) append_walls(speed);
    append_walls(q);
    append_walls(r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());

    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(q, answer.grid, left, right) ||
            !safe_midpoint(r, answer.grid, left, right)) continue;
        u32 failure = 0;
        for (unsigned bit = 0; bit < kPool.size(); ++bit)
            if (!safe_midpoint(kPool[bit], answer.grid, left, right))
                failure |= UINT32_C(1) << bit;
        if (std::popcount(failure) <= 9)
            answer.low_cells.emplace_back(failure, right - left);
    }
    return answer;
}

i128 activity_ticks(const Geometry& geometry, u32 mask) {
    i64 mass = 0;
    for (const auto& [failure, width] : geometry.low_cells)
        if ((failure & ~mask) == 0) mass += width;
    return static_cast<i128>(63) * mass -
           static_cast<i128>(4) * geometry.grid;
}

std::vector<u32> read_masks(const std::filesystem::path& path,
                            std::size_t expected, int rank_a, int rank_b,
                            u64 expected_fnv) {
    std::ifstream input(path);
    check(static_cast<bool>(input), "cannot open mask ledger");
    std::vector<u32> result;
    std::set<u32> distinct;
    Fnv hash;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_hex(token);
        const int rank = std::popcount(mask);
        check((rank == rank_a || rank == rank_b) && distinct.insert(mask).second,
              "mask rank/distinctness changed");
        result.push_back(mask);
        hash.add(mask);
    }
    check(input.eof() && result.size() == expected && hash.value == expected_fnv,
          "mask ledger identity changed");
    return result;
}

struct Pair { int q = 0; int r = 0; };

std::vector<Pair> read_rows(const std::filesystem::path& path) {
    std::ifstream input(path);
    check(static_cast<bool>(input), "cannot open endpoint-587 rows");
    std::vector<Pair> result;
    Fnv hash;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto comma = line.find(',');
        check(comma != std::string::npos, "malformed pair row");
        Pair pair{std::stoi(line.substr(0, comma)),
                  std::stoi(line.substr(comma + 1))};
        check(pair.q > 0 && pair.q < pair.r && pair.r == 587,
              "row escaped endpoint 587");
        if (!result.empty())
            check(std::tie(result.back().q, result.back().r) <
                      std::tie(pair.q, pair.r), "row order changed");
        result.push_back(pair);
        hash.add(pair.q); hash.add(pair.r);
    }
    check(input.eof() && result.size() == 10 &&
              hash.value == UINT64_C(0xf48ca5f1904d6f52),
          "endpoint-587 row identity changed");
    return result;
}

u32 next_combination(u32 value) {
    const u32 least = value & (u32{0} - value);
    const u32 ripple = value + least;
    return ripple | (((ripple ^ value) >> 2) / least);
}

struct Audit {
    u64 active = 0;
    u64 active_fnv = 0;
    u64 active_joint = 0;
    u64 active_nonjoint = 0;
    u64 exposed = 0;
    u64 exposed_fnv = 0;
    u64 hit_incidences = 0;
    u64 minimum_hits = std::numeric_limits<u64>::max();
    u64 maximum_hits = 0;
    u64 failures = 0;
    u64 failure_fnv = 0;
    std::vector<u32> failure_bodies;
};

Audit audit_pair(Pair pair, const std::vector<u32>& joint,
                 const std::vector<u32>& carrier,
                 const std::unordered_set<u32>& joint_set) {
    const Geometry geometry = make_geometry(pair.q, pair.r);
    std::vector<u32> active_joint;
    std::vector<u32> active_nonjoint;
    Audit result;
    Fnv active_hash;
    for (u32 mask : joint)
        if (activity_ticks(geometry, mask) >= 0) active_joint.push_back(mask);
    for (u32 mask : carrier) {
        if (activity_ticks(geometry, mask) < 0) continue;
        ++result.active;
        active_hash.add(mask);
        if (!joint_set.contains(mask)) active_nonjoint.push_back(mask);
    }
    result.active_fnv = active_hash.value;
    result.active_joint = active_joint.size();
    result.active_nonjoint = active_nonjoint.size();
    check(result.active == result.active_joint + result.active_nonjoint,
          "active joint partition failed");

    Fnv exposed_hash;
    Fnv failure_hash;
    u64 body_count = 0;
    const u32 limit = UINT32_C(1) << 30;
    for (u32 body = (UINT32_C(1) << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++body_count;
        bool joint_hit = false;
        for (u32 mask : active_joint) {
            if ((mask & body) == 0) {
                joint_hit = true;
                break;
            }
        }
        if (joint_hit) continue;

        ++result.exposed;
        exposed_hash.add(body);
        u64 hits = 0;
        for (u32 mask : active_nonjoint) hits += (mask & body) == 0;
        result.hit_incidences += hits;
        if (hits == 0) {
            ++result.failures;
            failure_hash.add(body);
            result.failure_bodies.push_back(body);
        } else {
            result.minimum_hits = std::min(result.minimum_hits, hits);
            result.maximum_hits = std::max(result.maximum_hits, hits);
        }
    }
    check(body_count == kBodyCount, "body universe changed");
    result.exposed_fnv = exposed_hash.value;
    result.failure_fnv = failure_hash.value;
    if (result.exposed == 0) result.minimum_hits = 0;
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        check(argc == 7,
              "usage: audit JOINT CARRIER ROWS PAIR_OUT FAILURE_OUT THREADS");
        const auto joint = read_masks(argv[1], 421, 8, 8,
                                      UINT64_C(0x20d63dd42fe8150e));
        const auto carrier = read_masks(argv[2], 3925, 8, 9,
                                        UINT64_C(0xeeae5518d84ccac5));
        const auto rows = read_rows(argv[3]);
        const unsigned thread_count = std::stoul(argv[6]);
        check(thread_count >= 1 && thread_count <= 64,
              "thread count escaped 1..64");
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        for (u32 mask : joint)
            check(std::find(carrier.begin(), carrier.end(), mask) != carrier.end(),
                  "protected joint mask absent from carrier");

        std::vector<Audit> audits(rows.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        const auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= rows.size()) return;
                    audits[index] = audit_pair(rows[index], joint, carrier,
                                               joint_set);
                }
            } catch (...) {
                std::lock_guard<std::mutex> guard(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
        };
        std::vector<std::thread> workers;
        for (unsigned index = 0; index < thread_count; ++index)
            workers.emplace_back(worker);
        for (auto& thread : workers) thread.join();
        if (worker_error) std::rethrow_exception(worker_error);

        std::ofstream pair_out(argv[4], std::ios::binary);
        std::ofstream failure_out(argv[5], std::ios::binary);
        check(pair_out && failure_out, "cannot create audit ledgers");
        pair_out << "q,r,active,active_fnv,active_joint,active_nonjoint,exposed,"
                    "exposed_fnv,minimum_hits,maximum_hits,failures,failure_fnv\n";
        failure_out << "q,r,body_hex\n";

        Fnv pair_hash;
        Fnv failure_hash;
        Fnv row_hash;
        u64 exposed_total = 0;
        u64 hit_total = 0;
        u64 failure_total = 0;
        u64 failed_rows = 0;
        for (std::size_t index = 0; index < rows.size(); ++index) {
            const Pair pair = rows[index];
            const Audit& audit = audits[index];
            row_hash.add(pair.q); row_hash.add(pair.r);
            exposed_total += audit.exposed;
            hit_total += audit.hit_incidences;
            failure_total += audit.failures;
            failed_rows += audit.failures != 0;
            pair_hash.add(pair.q); pair_hash.add(pair.r);
            pair_hash.add(audit.active); pair_hash.add(audit.active_fnv);
            pair_hash.add(audit.active_joint); pair_hash.add(audit.active_nonjoint);
            pair_hash.add(audit.exposed); pair_hash.add(audit.exposed_fnv);
            pair_hash.add(audit.minimum_hits); pair_hash.add(audit.maximum_hits);
            pair_hash.add(audit.failures); pair_hash.add(audit.failure_fnv);
            pair_out << pair.q << ',' << pair.r << ',' << audit.active << ','
                     << std::hex << audit.active_fnv << std::dec << ','
                     << audit.active_joint << ',' << audit.active_nonjoint << ','
                     << audit.exposed << ',' << std::hex << audit.exposed_fnv
                     << std::dec << ',' << audit.minimum_hits << ','
                     << audit.maximum_hits << ',' << audit.failures << ','
                     << std::hex << audit.failure_fnv << std::dec << '\n';
            for (u32 body : audit.failure_bodies) {
                failure_out << pair.q << ',' << pair.r << ',' << hex8(body)
                            << '\n';
                failure_hash.add(pair.q); failure_hash.add(pair.r);
                failure_hash.add(body);
            }
        }
        check(pair_out.good() && failure_out.good(), "audit ledger write failed");

        const u64 rank8 = std::count_if(carrier.begin(), carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        std::cout << "LRC14_ENDPOINT587_CLEANROOM_CARRIER_AUDIT_V1\n"
                  << "CARRIER " << carrier.size() << " FNV " << std::hex
                  << UINT64_C(0xeeae5518d84ccac5) << std::dec << " RANK8 "
                  << rank8 << " RANK9 " << carrier.size() - rank8
                  << " JOINT_RETAINED " << joint.size() << '\n'
                  << "ROWS " << rows.size() << " ENDPOINT 587 ROW_FNV "
                  << std::hex << row_hash.value << std::dec << " WORKERS "
                  << thread_count << " BODY_TESTS "
                  << rows.size() * kBodyCount << '\n';
        for (std::size_t index = 0; index < rows.size(); ++index) {
            const Pair pair = rows[index];
            const Audit& audit = audits[index];
            std::cout << "PAIR " << pair.q << ',' << pair.r << " ACTIVE "
                      << audit.active << " EXPOSED " << audit.exposed
                      << " HIT_RANGE " << audit.minimum_hits << ".."
                      << audit.maximum_hits << " FAILURES " << audit.failures
                      << " FAILURE_FNV " << std::hex << audit.failure_fnv
                      << std::dec << '\n';
        }
        std::cout << "SUMMARY EXPOSED " << exposed_total
                  << " HIT_INCIDENCES " << hit_total
                  << " FAILURES " << failure_total
                  << " FAILED_ROWS " << failed_rows
                  << " FAILURE_FNV " << std::hex << failure_hash.value
                  << " PAIR_LEDGER_FNV " << pair_hash.value << std::dec << '\n'
                  << "REPRESENTATION RAW_LOW_RANK_PAIR_SAFE_WALL_CELLS_"
                     "SERIALIZED_CARRIER_NO_SOURCE_IMPORT\n"
                  << "SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT587_LAYER_"
                     "UNCHANGED_THM4318_CARRIER_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT "
                  << (failure_total == 0 ? "PASS" : "HOSTILE_FAIL") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT587_CLEANROOM_AUDIT_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}
