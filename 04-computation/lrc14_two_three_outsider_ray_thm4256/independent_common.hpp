#ifndef LRC_23_RAY_INDEPENDENT_COMMON_HPP
#define LRC_23_RAY_INDEPENDENT_COMMON_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace ray_audit {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr std::array<int, 30> POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
constexpr i64 D = INT64_C(18241159416480);
constexpr i64 N = 84;
constexpr i64 S = 64;
constexpr i64 DELTA_C = 536;
constexpr u64 REPAIR_COUNT = UINT64_C(5852925);
constexpr u64 BODY_COUNT = UINT64_C(14307150);
constexpr u64 ORDER_SEED = UINT64_C(0x4245422842334245);

inline void ensure(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

inline std::string decimal(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    if (negative) value = -value;
    std::string result;
    while (value != 0) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

inline std::string hex64(u64 value) {
    std::ostringstream out;
    out << std::hex << std::nouppercase << std::setw(16) << std::setfill('0')
        << value;
    return out.str();
}

inline i64 exact_lcm(i64 a, i64 b) {
    const i64 divisor = std::gcd(a, b);
    const i128 result = static_cast<i128>(a / divisor) * b;
    ensure(result <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(result);
}

inline u64 splitmix64(u64 x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

class Fnv64 {
  public:
    void word(u64 value) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state_ ^= static_cast<std::uint8_t>((value >> (8 * byte)) & 0xffu);
            state_ *= UINT64_C(0x100000001b3);
        }
    }
    u64 value() const { return state_; }

  private:
    u64 state_ = UINT64_C(0xcbf29ce484222325);
};

inline u32 next_mask(u32 mask) {
    const u32 least = mask & (~mask + 1u);
    const u32 ripple = mask + least;
    if (ripple == 0) return 0;
    return ripple | (((mask ^ ripple) >> 2) / least);
}

inline std::string labels(u32 mask) {
    std::string answer;
    for (int bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        if (!answer.empty()) answer += ',';
        answer += std::to_string(POOL[bit]);
    }
    return answer;
}

struct Cell {
    i64 left = 0;
    i64 right = 0;
    u32 failed = 0;
};

inline bool safe_at_midpoint(int speed, i64 left, i64 right) {
    // The midpoint is (left+right)/(2D).  The residue below is twice its
    // fractional speed phase, on denominator 2D.
    const i64 modulus = 2 * D;
    const i64 residue = static_cast<i64>(
        (static_cast<i128>(speed) * (left + right)) % modulus);
    return static_cast<i128>(7) * residue >= D &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * D;
}

inline std::vector<Cell> build_pool_arrangement() {
    i64 recomputed = 1;
    for (int speed : POOL) recomputed = exact_lcm(recomputed, 14LL * speed);
    ensure(recomputed == D, "pool denominator mismatch");

    std::vector<i64> wall;
    wall.reserve(2 + 2 * std::accumulate(POOL.begin(), POOL.end(), 0));
    wall.push_back(0);
    wall.push_back(D);
    for (int speed : POOL) {
        const i64 quantum = D / (14LL * speed);
        for (int tooth = 0; tooth < speed; ++tooth) {
            wall.push_back((14LL * tooth + 1) * quantum);
            wall.push_back((14LL * tooth + 13) * quantum);
        }
    }
    std::sort(wall.begin(), wall.end());
    wall.erase(std::unique(wall.begin(), wall.end()), wall.end());
    ensure(wall.size() == 7134, "pool wall census mismatch");

    std::vector<Cell> result;
    result.reserve(wall.size() - 1);
    for (std::size_t i = 0; i + 1 < wall.size(); ++i) {
        u32 failed = 0;
        for (int bit = 0; bit < 30; ++bit) {
            if (!safe_at_midpoint(POOL[bit], wall[i], wall[i + 1])) {
                failed |= u32{1} << bit;
            }
        }
        result.push_back({wall[i], wall[i + 1], failed});
    }
    ensure(!result.empty() && result.front().failed == (u32{1} << 30) - 1,
           "origin-neighbourhood control failed");
    ensure(result.back().failed == (u32{1} << 30) - 1,
           "cyclic origin-neighbourhood control failed");
    return result;
}

// For A=G_2 intersect G_3, the unsafe arcs on the 84-grid merge to
// [0,3], [26,30], [39,45], [54,58], [81,84].  Thus A is the following
// four-interval union.  safe_prefix_numerator(r) is 84D times the safe
// length in [0,r/D].
constexpr std::array<std::pair<i64, i64>, 4> PRIMITIVE_ARCS = {
    std::pair<i64, i64>{3, 26}, {30, 39}, {45, 54}, {58, 81}};

inline i128 primitive_safe_prefix_numerator(i64 remainder) {
    ensure(0 <= remainder && remainder < D, "primitive phase out of range");
    const i128 y = static_cast<i128>(N) * remainder;
    i128 answer = 0;
    for (const auto& [a, b] : PRIMITIVE_ARCS) {
        const i128 left = static_cast<i128>(a) * D;
        const i128 right = static_cast<i128>(b) * D;
        if (y <= left) continue;
        answer += std::min(y, right) - left;
    }
    return answer;
}

// 84D times integral_0^(z/D) 1_A(t)dt, where z may span many periods.
inline i128 primitive_integral_numerator(i128 z) {
    ensure(z >= 0, "negative primitive coordinate");
    const i128 whole = z / D;
    const i64 remainder = static_cast<i64>(z % D);
    return whole * S * D + primitive_safe_prefix_numerator(remainder);
}

inline void check_primitive_constants() {
    i64 safe_ticks = 0;
    for (const auto& [a, b] : PRIMITIVE_ARCS) safe_ticks += b - a;
    ensure(safe_ticks == S, "primitive safe length mismatch");

    // C(t)=84*safe_length([0,t])-64*t at primitive wall ticks.
    constexpr std::array<i64, 12> points = {0, 2, 3, 26, 30, 39,
                                            45, 54, 58, 81, 82, 84};
    i64 minimum = 0;
    i64 maximum = 0;
    for (i64 point : points) {
        i64 prefix = 0;
        for (const auto& [a, b] : PRIMITIVE_ARCS) {
            prefix += std::max<i64>(0, std::min(point, b) - a);
        }
        const i64 c = N * prefix - S * point;
        minimum = std::min(minimum, c);
        maximum = std::max(maximum, c);
    }
    ensure(maximum - minimum == DELTA_C && minimum == -268 && maximum == 268,
           "primitive centered range mismatch");
    ensure(primitive_safe_prefix_numerator(0) == 0,
           "primitive prefix origin mismatch");
}

struct Geometry {
    u32 repair = 0;
    i64 width = 0;
    std::vector<std::pair<i64, i64>> components;
};

inline Geometry direct_geometry(u32 repair, const std::vector<Cell>& cells) {
    ensure(std::popcount(repair) == 8, "repair arity mismatch");
    Geometry result;
    result.repair = repair;
    std::size_t i = 0;
    while (i < cells.size()) {
        const bool safe = (cells[i].failed & ~repair) == 0;
        if (!safe) {
            ++i;
            continue;
        }
        const i64 left = cells[i].left;
        i64 right = cells[i].right;
        ++i;
        while (i < cells.size() && (cells[i].failed & ~repair) == 0) {
            right = cells[i].right;
            ++i;
        }
        result.components.push_back({left, right});
        result.width += right - left;
    }
    ensure(result.width > 0 && !result.components.empty(),
           "empty positive-length repair geometry");
    return result;
}

// Numerator on denominator 84*D*g of the exact mass of
// U_R intersect G_(2g) intersect G_(3g).
inline i128 exact_mass_numerator(const Geometry& geometry, i64 g) {
    ensure(g > 0, "nonpositive scale");
    i128 numerator = 0;
    for (const auto& [left, right] : geometry.components) {
        numerator += primitive_integral_numerator(static_cast<i128>(g) * right)
                   - primitive_integral_numerator(static_cast<i128>(g) * left);
    }
    return numerator;
}

inline i128 exact_margin(const Geometry& geometry, i64 g) {
    return static_cast<i128>(63) * exact_mass_numerator(geometry, g)
         - static_cast<i128>(4) * N * D * g;
}

struct BodyScan {
    u64 bodies = 0;
    u64 failures = 0;
    u64 checks = 0;
    u64 max_checks = 0;
    u32 worst_body = 0;
    u32 worst_repair = 0;
    u32 first_failure = std::numeric_limits<u32>::max();
};

inline std::vector<u32> all_bodies() {
    std::vector<u32> result;
    result.reserve(BODY_COUNT);
    u32 body = (u32{1} << 9) - 1;
    const u32 limit = u32{1} << 30;
    while (body != 0 && body < limit) {
        result.push_back(body);
        const u32 following = next_mask(body);
        if (following <= body) break;
        body = following;
    }
    ensure(result.size() == BODY_COUNT, "body census mismatch");
    return result;
}

inline BodyScan scan_bodies(const std::vector<u32>& deck,
                            const std::vector<u32>& bodies) {
    ensure(bodies.size() == BODY_COUNT, "incomplete body universe");
    const unsigned available = std::thread::hardware_concurrency();
    const unsigned threads = std::max(1u, std::min(16u, available ? available : 1u));
    std::vector<BodyScan> local(threads);
    std::vector<std::thread> worker;
    worker.reserve(threads);
    for (unsigned lane = 0; lane < threads; ++lane) {
        worker.emplace_back([&, lane]() {
            BodyScan report;
            const std::size_t begin = bodies.size() * lane / threads;
            const std::size_t end = bodies.size() * (lane + 1) / threads;
            for (std::size_t index = begin; index < end; ++index) {
                const u32 body = bodies[index];
                u64 used = 0;
                u32 witness = 0;
                for (u32 repair : deck) {
                    ++used;
                    if ((repair & body) == 0) {
                        witness = repair;
                        break;
                    }
                }
                ++report.bodies;
                report.checks += used;
                if (witness == 0) {
                    ++report.failures;
                    report.first_failure = std::min(report.first_failure, body);
                }
                if (used > report.max_checks ||
                    (used == report.max_checks && body < report.worst_body)) {
                    report.max_checks = used;
                    report.worst_body = body;
                    report.worst_repair = witness;
                }
            }
            local[lane] = report;
        });
    }
    for (std::thread& item : worker) item.join();

    BodyScan result;
    for (const BodyScan& report : local) {
        result.bodies += report.bodies;
        result.failures += report.failures;
        result.checks += report.checks;
        result.first_failure = std::min(result.first_failure, report.first_failure);
        if (report.max_checks > result.max_checks ||
            (report.max_checks == result.max_checks &&
             report.worst_body < result.worst_body)) {
            result.max_checks = report.max_checks;
            result.worst_body = report.worst_body;
            result.worst_repair = report.worst_repair;
        }
    }
    return result;
}

}  // namespace ray_audit

#endif
