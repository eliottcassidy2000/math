// Clean-room endpoint-588 audit.  It imports no maintained audit source and
// consumes only the inherited carrier/joint mask ledgers and the typed rows.

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
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_set>
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
    u64 state = kFnvBasis;
    void add(u64 word) {
        for (unsigned byte = 0; byte != 8; ++byte) {
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
std::string hex8(u32 word) {
    std::ostringstream out;
    out << std::hex << std::setw(8) << std::setfill('0') << word;
    return out.str();
}
i64 checked_lcm(i64 a, i64 b) {
    const i128 value = static_cast<i128>(a / std::gcd(a, b)) * b;
    check(value <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(value);
}
bool is_safe(int speed, i64 grid, i64 left, i64 right) {
    i128 twice_phase = static_cast<i128>(speed) * (left + static_cast<i128>(right));
    twice_phase %= static_cast<i128>(2) * grid;
    if (twice_phase < 0) twice_phase += static_cast<i128>(2) * grid;
    return static_cast<i128>(7) * twice_phase >= grid &&
           static_cast<i128>(7) * twice_phase <= static_cast<i128>(13) * grid;
}
struct Geometry {
    i64 grid = 0;
    std::vector<std::pair<u32, i64>> low_classes;
};
Geometry make_geometry(int q, int r) {
    Geometry geometry;
    geometry.grid = 1;
    for (int speed : kPool)
        geometry.grid = checked_lcm(geometry.grid, 14LL * speed);
    geometry.grid = checked_lcm(geometry.grid, 14LL * q);
    geometry.grid = checked_lcm(geometry.grid, 14LL * r);
    std::vector<i64> walls{0, geometry.grid};
    const auto append = [&](int speed) {
        check(geometry.grid % (14LL * speed) == 0, "nonintegral wall unit");
        const i64 unit = geometry.grid / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) append(speed);
    append(q);
    append(r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::map<u32, i64> aggregated;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1], right = walls[index];
        if (!is_safe(q, geometry.grid, left, right) ||
            !is_safe(r, geometry.grid, left, right)) continue;
        u32 failure = 0;
        for (unsigned bit = 0; bit < kPool.size(); ++bit)
            if (!is_safe(kPool[bit], geometry.grid, left, right))
                failure |= UINT32_C(1) << bit;
        if (std::popcount(failure) <= 9)
            aggregated[failure] += right - left;
    }
    geometry.low_classes.assign(aggregated.begin(), aggregated.end());
    return geometry;
}
i128 activity_ticks(const Geometry& geometry, u32 mask) {
    i64 mass = 0;
    for (const auto& [failure, width] : geometry.low_classes)
        if ((failure & ~mask) == 0) mass += width;
    return static_cast<i128>(63) * mass - static_cast<i128>(4) * geometry.grid;
}
std::vector<u32> read_masks(const std::filesystem::path& path,
                            std::size_t expected, int allowed_rank_a,
                            int allowed_rank_b, u64 expected_fnv) {
    std::ifstream input(path);
    check(static_cast<bool>(input), "cannot open mask ledger");
    std::vector<u32> result;
    std::set<u32> distinct;
    Fnv ledger;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_hex(token);
        const int rank = std::popcount(mask);
        check((rank == allowed_rank_a || rank == allowed_rank_b) &&
                  distinct.insert(mask).second,
              "mask rank/distinctness changed");
        result.push_back(mask);
        ledger.add(mask);
    }
    check(input.eof() && result.size() == expected &&
              ledger.state == expected_fnv,
          "mask ledger identity changed");
    return result;
}
struct Pair { int q = 0; int r = 0; };
std::vector<Pair> read_rows(const std::filesystem::path& path) {
    std::ifstream input(path);
    check(static_cast<bool>(input), "cannot open row ledger");
    std::vector<Pair> rows;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto comma = line.find(',');
        check(comma != std::string::npos, "malformed row");
        Pair pair{std::stoi(line.substr(0, comma)),
                  std::stoi(line.substr(comma + 1))};
        check(pair.q > 0 && pair.q < pair.r && pair.r == 588,
              "row escaped endpoint 588");
        if (!rows.empty())
            check(std::tie(rows.back().q, rows.back().r) <
                      std::tie(pair.q, pair.r), "row order changed");
        rows.push_back(pair);
        ledger.add(pair.q); ledger.add(pair.r);
    }
    check(input.eof() && rows.size() == 66 &&
              ledger.state == UINT64_C(0x18cf9a572cf9a5be),
          "typed endpoint-588 row identity changed");
    return rows;
}
u32 next_combination(u32 value) {
    const u32 smallest = value & (u32{0} - value);
    const u32 ripple = value + smallest;
    return ripple | (((value ^ ripple) >> 2) / smallest);
}
struct Audit {
    u64 active = 0, active_fnv = 0, active_joint = 0, active_nonjoint = 0;
    u64 exposed = 0, exposed_fnv = 0, hit_incidences = 0;
    u64 minimum_hits = std::numeric_limits<u64>::max(), maximum_hits = 0;
    u64 failures = 0, failure_fnv = 0;
    std::vector<u32> bodies;
};
Audit audit_row(Pair pair, const std::vector<u32>& joint,
                const std::vector<u32>& carrier,
                const std::unordered_set<u32>& joint_set) {
    const Geometry geometry = make_geometry(pair.q, pair.r);
    std::vector<u32> active_joint, active_nonjoint;
    Audit answer;
    Fnv active_hash;
    for (u32 mask : joint)
        if (activity_ticks(geometry, mask) >= 0) active_joint.push_back(mask);
    for (u32 mask : carrier) {
        if (activity_ticks(geometry, mask) < 0) continue;
        ++answer.active; active_hash.add(mask);
        if (!joint_set.contains(mask)) active_nonjoint.push_back(mask);
    }
    answer.active_fnv = active_hash.state;
    answer.active_joint = active_joint.size();
    answer.active_nonjoint = active_nonjoint.size();
    check(answer.active == answer.active_joint + answer.active_nonjoint,
          "active partition failed");
    Fnv exposed_hash, failure_hash;
    u64 tested = 0;
    const u32 limit = UINT32_C(1) << 30;
    for (u32 body = (UINT32_C(1) << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++tested;
        bool shielded = false;
        for (u32 mask : active_joint)
            if ((mask & body) == 0) { shielded = true; break; }
        if (shielded) continue;
        ++answer.exposed; exposed_hash.add(body);
        u64 hits = 0;
        for (u32 mask : active_nonjoint) hits += (mask & body) == 0;
        answer.hit_incidences += hits;
        if (hits == 0) {
            ++answer.failures; failure_hash.add(body); answer.bodies.push_back(body);
        } else {
            answer.minimum_hits = std::min(answer.minimum_hits, hits);
            answer.maximum_hits = std::max(answer.maximum_hits, hits);
        }
    }
    check(tested == kBodyCount, "body universe changed");
    answer.exposed_fnv = exposed_hash.state;
    answer.failure_fnv = failure_hash.state;
    if (answer.exposed == 0) answer.minimum_hits = 0;
    return answer;
}
}

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
        check(thread_count >= 1 && thread_count <= 64, "thread count escaped");
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        for (u32 mask : joint)
            check(std::find(carrier.begin(), carrier.end(), mask) != carrier.end(),
                  "joint mask absent from carrier");
        std::vector<Audit> audits(rows.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr failure;
        std::mutex failure_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const auto index = cursor.fetch_add(1);
                    if (index >= rows.size()) return;
                    audits[index] = audit_row(rows[index], joint, carrier, joint_set);
                }
            } catch (...) {
                std::lock_guard<std::mutex> guard(failure_mutex);
                if (!failure) failure = std::current_exception();
            }
        };
        std::vector<std::thread> threads;
        for (unsigned i = 0; i < thread_count; ++i) threads.emplace_back(worker);
        for (auto& thread : threads) thread.join();
        if (failure) std::rethrow_exception(failure);

        std::ofstream pairs(argv[4], std::ios::binary);
        std::ofstream failures(argv[5], std::ios::binary);
        check(pairs && failures, "cannot create outputs");
        pairs << "q,r,active,active_fnv,active_joint,active_nonjoint,exposed,"
                 "exposed_fnv,minimum_hits,maximum_hits,failures,failure_fnv\n";
        failures << "q,r,body_hex\n";
        Fnv pair_hash, failure_hash;
        u64 total_exposed = 0, total_hits = 0, total_failures = 0, failed_rows = 0;
        std::cout << "LRC14_ENDPOINT588_CLEANROOM_CARRIER_AUDIT_V1\n"
                  << "CARRIER 3925 FNV eeae5518d84ccac5 JOINT 421\n"
                  << "ROWS 66 ENDPOINT 588 ROW_FNV 18cf9a572cf9a5be THREADS "
                  << thread_count << " BODY_TESTS " << 66 * kBodyCount << '\n';
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const auto pair = rows[i]; const auto& audit = audits[i];
            total_exposed += audit.exposed; total_hits += audit.hit_incidences;
            total_failures += audit.failures; failed_rows += audit.failures != 0;
            pair_hash.add(pair.q); pair_hash.add(pair.r); pair_hash.add(audit.active);
            pair_hash.add(audit.active_fnv); pair_hash.add(audit.active_joint);
            pair_hash.add(audit.active_nonjoint); pair_hash.add(audit.exposed);
            pair_hash.add(audit.exposed_fnv); pair_hash.add(audit.minimum_hits);
            pair_hash.add(audit.maximum_hits); pair_hash.add(audit.failures);
            pair_hash.add(audit.failure_fnv);
            pairs << pair.q << ',' << pair.r << ',' << audit.active << ','
                  << std::hex << audit.active_fnv << std::dec << ','
                  << audit.active_joint << ',' << audit.active_nonjoint << ','
                  << audit.exposed << ',' << std::hex << audit.exposed_fnv
                  << std::dec << ',' << audit.minimum_hits << ','
                  << audit.maximum_hits << ',' << audit.failures << ','
                  << std::hex << audit.failure_fnv << std::dec << '\n';
            for (u32 body : audit.bodies) {
                failures << pair.q << ',' << pair.r << ',' << hex8(body) << '\n';
                failure_hash.add(pair.q); failure_hash.add(pair.r); failure_hash.add(body);
            }
            std::cout << "PAIR " << pair.q << ",588 ACTIVE " << audit.active
                      << " EXPOSED " << audit.exposed << " HIT_RANGE "
                      << audit.minimum_hits << ".." << audit.maximum_hits
                      << " FAILURES " << audit.failures << " FAILURE_FNV "
                      << std::hex << audit.failure_fnv << std::dec << '\n';
        }
        check(pairs.good() && failures.good(), "output write failed");
        std::cout << "SUMMARY EXPOSED " << total_exposed << " HIT_INCIDENCES "
                  << total_hits << " FAILURES " << total_failures
                  << " FAILED_ROWS " << failed_rows << " FAILURE_FNV "
                  << std::hex << failure_hash.state << " PAIR_LEDGER_FNV "
                  << pair_hash.state << std::dec << '\n'
                  << "SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT588_LAYER_"
                     "INHERITED_THM4318_CARRIER_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT " << (total_failures ? "HOSTILE_FAIL" : "PASS") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CLEANROOM_CARRIER_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}
