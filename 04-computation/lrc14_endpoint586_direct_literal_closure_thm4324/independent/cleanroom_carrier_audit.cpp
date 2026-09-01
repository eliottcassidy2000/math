// Import-free replay of the serialized THM-4318 carrier at endpoint 586.
// Activity is reconstructed from exact unaggregated pair-safe wall cells.

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

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
struct Fnv {
    u64 state = kFnvBasis;
    void add(u64 word) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (word >> (8 * byte)) & UINT64_C(255);
            state *= kFnvPrime;
        }
    }
};
u32 parse_hex(const std::string& token) {
    const u64 wide = std::stoull(token, nullptr, 16);
    check(wide < (UINT64_C(1) << 30), "mask escaped 30 labels");
    return static_cast<u32>(wide);
}
std::string hex8(u32 value) {
    std::ostringstream out;
    out << std::hex << std::setw(8) << std::setfill('0') << value;
    return out.str();
}
i64 checked_lcm(i64 a, i64 b) {
    const i128 wide = static_cast<i128>(a / std::gcd(a, b)) * b;
    check(wide <= std::numeric_limits<i64>::max(), "wall grid overflow");
    return static_cast<i64>(wide);
}
bool safe_midpoint(int speed, i64 grid, i64 left, i64 right) {
    i128 phase2 = static_cast<i128>(speed) *
                  (static_cast<i128>(left) + right);
    phase2 %= static_cast<i128>(2) * grid;
    if (phase2 < 0) phase2 += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * phase2 >= grid &&
           static_cast<i128>(7) * phase2 <= static_cast<i128>(13) * grid;
}
struct Geometry {
    i64 grid = 1;
    std::vector<std::pair<u32, i64>> low_cells;
};
Geometry make_geometry(int q) {
    Geometry geometry;
    for (int speed : kPool)
        geometry.grid = checked_lcm(geometry.grid, 14LL * speed);
    geometry.grid = checked_lcm(geometry.grid, 14LL * q);
    geometry.grid = checked_lcm(geometry.grid, 14LL * 586);
    std::vector<i64> walls{0, geometry.grid};
    const auto add_walls = [&](int speed) {
        const i64 divisor = 14LL * speed;
        check(geometry.grid % divisor == 0, "nonintegral wall unit");
        const i64 unit = geometry.grid / divisor;
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(q);
    add_walls(586);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1], right = walls[index];
        if (!safe_midpoint(q, geometry.grid, left, right) ||
            !safe_midpoint(586, geometry.grid, left, right)) continue;
        u32 failed = 0;
        for (unsigned bit = 0; bit < kPool.size(); ++bit)
            if (!safe_midpoint(kPool[bit], geometry.grid, left, right))
                failed |= UINT32_C(1) << bit;
        if (std::popcount(failed) <= 9)
            geometry.low_cells.emplace_back(failed, right - left);
    }
    return geometry;
}
i128 activity_ticks(const Geometry& geometry, u32 repair) {
    i64 mass = 0;
    for (const auto& [failed, width] : geometry.low_cells)
        if ((failed & ~repair) == 0) mass += width;
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
        result.push_back(mask); hash.add(mask);
    }
    check(input.eof() && result.size() == expected && hash.state == expected_fnv,
          "mask ledger identity changed");
    return result;
}
struct Pair { int q = 0, r = 0; };
std::vector<Pair> read_rows(const std::filesystem::path& path) {
    std::ifstream input(path);
    check(static_cast<bool>(input), "cannot open reconstructed frontier");
    std::vector<Pair> rows;
    Fnv hash;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto comma = line.find(',');
        check(comma != std::string::npos &&
              line.find(',', comma + 1) == std::string::npos,
              "malformed frontier row");
        const Pair pair{std::stoi(line.substr(0, comma)),
                        std::stoi(line.substr(comma + 1))};
        check(pair.q > 0 && pair.q < pair.r && pair.r == 586,
              "frontier row escaped endpoint 586");
        if (!rows.empty())
            check(std::tie(rows.back().q, rows.back().r) <
                  std::tie(pair.q, pair.r), "frontier order changed");
        rows.push_back(pair); hash.add(pair.q); hash.add(pair.r);
    }
    check(input.eof() && rows.size() == 12 &&
          hash.state == UINT64_C(0xa1b617faa2e7f63f),
          "endpoint-586 frontier identity changed");
    return rows;
}
u32 next_combination(u32 value) {
    const u32 least = value & (u32{0} - value);
    const u32 ripple = value + least;
    return ripple | (((ripple ^ value) >> 2) / least);
}
struct Audit {
    u64 active = 0, active_fnv = 0, active_joint = 0, active_nonjoint = 0;
    u64 exposed = 0, exposed_fnv = 0, hit_incidences = 0;
    u64 minimum_hits = std::numeric_limits<u64>::max(), maximum_hits = 0;
    u64 failures = 0, failure_fnv = 0;
    std::vector<u32> failure_bodies;
};
Audit audit_pair(Pair pair, const std::vector<u32>& joint,
                 const std::vector<u32>& carrier,
                 const std::unordered_set<u32>& joint_set) {
    const Geometry geometry = make_geometry(pair.q);
    std::vector<u32> active_joint, active_nonjoint;
    Audit result;
    Fnv active_hash;
    for (u32 mask : joint)
        if (activity_ticks(geometry, mask) >= 0) active_joint.push_back(mask);
    for (u32 mask : carrier) {
        if (activity_ticks(geometry, mask) < 0) continue;
        ++result.active; active_hash.add(mask);
        if (!joint_set.contains(mask)) active_nonjoint.push_back(mask);
    }
    result.active_fnv = active_hash.state;
    result.active_joint = active_joint.size();
    result.active_nonjoint = active_nonjoint.size();
    check(result.active == result.active_joint + result.active_nonjoint,
          "active partition failed");
    Fnv exposed_hash, failure_hash;
    u64 body_count = 0;
    const u32 limit = UINT32_C(1) << 30;
    for (u32 body = (UINT32_C(1) << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++body_count;
        bool protected_hit = false;
        for (u32 mask : active_joint)
            if ((mask & body) == 0) { protected_hit = true; break; }
        if (protected_hit) continue;
        ++result.exposed; exposed_hash.add(body);
        u64 hits = 0;
        for (u32 mask : active_nonjoint) hits += (mask & body) == 0;
        result.hit_incidences += hits;
        if (hits == 0) {
            ++result.failures; failure_hash.add(body);
            result.failure_bodies.push_back(body);
        } else {
            result.minimum_hits = std::min(result.minimum_hits, hits);
            result.maximum_hits = std::max(result.maximum_hits, hits);
        }
    }
    check(body_count == kBodyCount, "body universe changed");
    result.exposed_fnv = exposed_hash.state;
    result.failure_fnv = failure_hash.state;
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
        const unsigned threads = std::stoul(argv[6]);
        check(threads >= 1 && threads <= 64, "thread count escaped 1..64");
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        for (u32 mask : joint)
            check(std::find(carrier.begin(), carrier.end(), mask) != carrier.end(),
                  "protected mask absent from carrier");
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
        for (unsigned i = 0; i < threads; ++i) workers.emplace_back(worker);
        for (auto& thread : workers) thread.join();
        if (worker_error) std::rethrow_exception(worker_error);

        std::ofstream pair_out(argv[4]);
        std::ofstream failure_out(argv[5]);
        check(pair_out && failure_out, "cannot create audit ledgers");
        pair_out << "q,r,active,active_fnv,active_joint,active_nonjoint,exposed,"
                    "exposed_fnv,minimum_hits,maximum_hits,failures,failure_fnv\n";
        failure_out << "q,r,body_hex\n";
        Fnv pair_hash, failure_hash, row_hash;
        u64 exposed_total = 0, hit_total = 0, failure_total = 0, failed_rows = 0;
        for (std::size_t index = 0; index < rows.size(); ++index) {
            const Pair pair = rows[index]; const Audit& audit = audits[index];
            row_hash.add(pair.q); row_hash.add(pair.r);
            exposed_total += audit.exposed; hit_total += audit.hit_incidences;
            failure_total += audit.failures; failed_rows += audit.failures != 0;
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
                failure_out << pair.q << ',' << pair.r << ',' << hex8(body) << '\n';
                failure_hash.add(pair.q); failure_hash.add(pair.r);
                failure_hash.add(body);
            }
        }
        check(pair_out.good() && failure_out.good(), "audit ledger write failed");
        const u64 rank8 = std::count_if(carrier.begin(), carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        std::cout << "LRC14_ENDPOINT586_CLEANROOM_CARRIER_AUDIT_V1\n"
                  << "CARRIER " << carrier.size() << " FNV " << std::hex
                  << UINT64_C(0xeeae5518d84ccac5) << std::dec << " RANK8 "
                  << rank8 << " RANK9 " << carrier.size() - rank8
                  << " JOINT_RETAINED " << joint.size() << '\n'
                  << "ROWS " << rows.size() << " ENDPOINT 586 ROW_FNV "
                  << std::hex << row_hash.state << std::dec << " WORKERS "
                  << threads << " BODY_TESTS " << rows.size() * kBodyCount << '\n';
        for (std::size_t index = 0; index < rows.size(); ++index)
            std::cout << "PAIR " << rows[index].q << ",586 ACTIVE "
                      << audits[index].active << " EXPOSED " << audits[index].exposed
                      << " HIT_RANGE " << audits[index].minimum_hits << ".."
                      << audits[index].maximum_hits << " FAILURES "
                      << audits[index].failures << " FAILURE_FNV " << std::hex
                      << audits[index].failure_fnv << std::dec << '\n';
        std::cout << "SUMMARY EXPOSED " << exposed_total << " HIT_INCIDENCES "
                  << hit_total << " FAILURES " << failure_total << " FAILED_ROWS "
                  << failed_rows << " FAILURE_FNV " << std::hex
                  << failure_hash.state << " PAIR_LEDGER_FNV " << pair_hash.state
                  << std::dec << '\n'
                  << "REPRESENTATION RAW_LOW_RANK_PAIR_SAFE_WALL_CELLS_"
                     "SERIALIZED_CARRIER_NO_SOURCE_IMPORT\n"
                  << "SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT586_LAYER_"
                     "UNCHANGED_THM4318_CARRIER_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT " << (failure_total == 0 ? "PASS" : "HOSTILE_FAIL")
                  << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT586_CLEANROOM_CARRIER_ERROR " << error.what() << '\n';
        return 1;
    }
}
