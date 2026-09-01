#ifndef ENDPOINT589_Q50_RESPONSE_COMMON_HPP
#define ENDPOINT589_Q50_RESPONSE_COMMON_HPP

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
#include <ranges>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace q50_common {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr std::array<int, 30> kPool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290};

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
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

std::string hex8(u32 value) {
    std::ostringstream output;
    output << std::hex << std::setfill('0') << std::setw(8) << value;
    return output.str();
}

u32 parse_mask_agent(const std::string& text) {
    std::size_t consumed = 0;
    const unsigned long value = std::stoul(text, &consumed, 16);
    require(consumed == text.size() && value < (UINT64_C(1) << 30),
            "malformed 30-bit mask");
    return static_cast<u32>(value);
}

u32 next_combination(u32 value) {
    const u32 low = value & (0U - value);
    const u32 ripple = value + low;
    if (ripple == 0) return 0;
    return ripple | (((value ^ ripple) >> 2) / low);
}

i64 checked_lcm(i64 left, i64 right) {
    require(left > 0 && right > 0, "nonpositive LCM input");
    const i128 value = static_cast<i128>(left / std::gcd(left, right)) * right;
    require(value <= std::numeric_limits<i64>::max(), "LCM overflow");
    return static_cast<i64>(value);
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
    std::vector<std::pair<u32, i64>> low_classes;
};

Geometry build_geometry(int q, int r) {
    Geometry geometry;
    geometry.grid = 1;
    for (int speed : kPool)
        geometry.grid = checked_lcm(geometry.grid, 14LL * speed);
    geometry.grid = checked_lcm(geometry.grid, 14LL * q);
    geometry.grid = checked_lcm(geometry.grid, 14LL * r);

    std::vector<i64> walls = {0, geometry.grid};
    const auto add_walls = [&](int speed) {
        const i64 divisor = 14LL * speed;
        require(geometry.grid % divisor == 0, "nonintegral wall unit");
        const i64 unit = geometry.grid / divisor;
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(q);
    add_walls(r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.front() == 0 && walls.back() == geometry.grid,
            "wall boundary changed");

    std::map<u32, i64> by_failure;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(q, geometry.grid, left, right) ||
            !safe_midpoint(r, geometry.grid, left, right))
            continue;
        u32 failure = 0;
        for (unsigned bit = 0; bit < kPool.size(); ++bit)
            if (!safe_midpoint(kPool[bit], geometry.grid, left, right))
                failure |= UINT32_C(1) << bit;
        if (std::popcount(failure) <= 9)
            by_failure[failure] += right - left;
    }

    geometry.low_classes.assign(by_failure.begin(), by_failure.end());

    Fnv ledger;
    for (const auto& [failure, width] : geometry.low_classes) {
        ledger.add(failure);
        ledger.add(width);
    }
    require(q == 50 && r == 589 &&
                geometry.grid == INT64_C(2827379709554400) &&
                geometry.low_classes.size() == 2383 &&
                ledger.state == UINT64_C(0x88d3eb2d7a477232),
            "q50 geometry changed");
    return geometry;
}

struct Margin {
    i64 mass = 0;
    i128 ticks = 0;
};

Margin margin(const Geometry& geometry, u32 deletion) {
    Margin result;
    for (const auto& [failure, width] : geometry.low_classes)
        if ((failure & ~deletion) == 0) result.mass += width;
    result.ticks = static_cast<i128>(63) * result.mass -
                   static_cast<i128>(4) * geometry.grid;
    return result;
}

}  // namespace q50_common

#endif

