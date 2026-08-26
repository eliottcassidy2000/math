#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static std::string show(u128 x) {
  if (x == 0) return "0";
  std::string s;
  while (x != 0) {
    s.push_back(static_cast<char>('0' + x % 10));
    x /= 10;
  }
  std::reverse(s.begin(), s.end());
  return s;
}

static u128 checked_mul(u128 x, u128 y) {
  const u128 top = ~u128{0};
  if (x != 0 && y > top / x) {
    throw std::overflow_error("unsigned-128 cross-product overflow");
  }
  return x * y;
}

struct Stats {
  u64 H{};
  u64 W{};
  std::vector<u64> v;
  u128 m{};
  u128 rhs{};
};

static std::size_t order_of(const std::string& bits) {
  std::size_t n = 0;
  while (n * (n - 1) / 2 < bits.size()) ++n;
  if (n * (n - 1) / 2 != bits.size()) {
    throw std::runtime_error("bad upper-triangle tournament word length");
  }
  return n;
}

static Stats evaluate(const std::string& bits) {
  const std::size_t n = order_of(bits);
  if (n == 0 || n > 9) {
    throw std::runtime_error("this finite audit is intentionally restricted to 1<=n<=9");
  }
  std::vector<u64> out(n, 0);
  std::size_t p = 0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j, ++p) {
      if (bits[p] == '1') out[i] |= u64{1} << j;
      else if (bits[p] == '0') out[j] |= u64{1} << i;
      else throw std::runtime_error("nonbinary tournament word");
    }
  }

  const u64 count = u64{1} << n;
  const u64 full = count - 1;
  std::vector<u64> finish(static_cast<std::size_t>(count) * n, 0);
  const auto ix = [n](u64 mask, std::size_t i) {
    return static_cast<std::size_t>(mask) * n + i;
  };
  for (std::size_t i = 0; i < n; ++i) finish[ix(u64{1} << i, i)] = 1;
  std::vector<u64> h(count, 0);
  h[0] = 1;
  for (u64 mask = 1; mask < count; ++mask) {
    for (std::size_t last = 0; last < n; ++last) {
      if ((mask & (u64{1} << last)) == 0) continue;
      const u64 rest = mask ^ (u64{1} << last);
      if (rest == 0) continue;
      u64 value = 0;
      for (std::size_t prev = 0; prev < n; ++prev) {
        if ((rest & (u64{1} << prev)) != 0 &&
            (out[prev] & (u64{1} << last)) != 0) {
          value += finish[ix(rest, prev)];
        }
      }
      finish[ix(mask, last)] = value;
    }
    for (std::size_t last = 0; last < n; ++last) h[mask] += finish[ix(mask, last)];
  }

  Stats s;
  s.H = h[full];
  s.v.assign(n, 0);
  for (u64 mask = 1; mask < count; ++mask) {
    const u64 complement_h = h[full ^ mask];
    for (std::size_t i = 0; i < n; ++i) {
      if ((mask & (u64{1} << i)) != 0) {
        s.v[i] += finish[ix(mask, i)] * complement_h;
      }
    }
  }
  u64 a = 0;
  for (u64 value : s.v) {
    a += value;
    s.m += u128{value} * value;
  }
  if (a < s.H) throw std::runtime_error("negative W");
  s.W = a - s.H;
  s.rhs = u128{s.W + s.H} * u128{s.W + 3 * s.H};
  return s;
}

int main() {
  // The caller supplies one canonical representative of every strong class
  // via `gentourng -q -c n`.  This evaluator deliberately does not duplicate
  // nauty's strong filter; transitive hostile controls can therefore be piped
  // to it separately.
  std::string bits;
  std::size_t rows = 0;
  std::size_t bad = 0;
  std::size_t equal = 0;
  std::size_t n = 0;
  u128 min_gap = ~u128{0};
  u128 ratio_num = 0;
  u128 ratio_den = 1;
  std::string first_bad_word;
  std::string min_gap_word;
  std::string min_ratio_word;
  Stats first_bad_stats;
  Stats min_gap_stats;
  Stats min_ratio_stats;

  while (std::cin >> bits) {
    if (n == 0) n = order_of(bits);
    if (order_of(bits) != n) throw std::runtime_error("mixed tournament orders");
    const Stats s = evaluate(bits);
    ++rows;
    const u128 five_m = 5 * s.m;
    if (s.rhs < five_m) {
      ++bad;
      if (first_bad_word.empty()) {
        first_bad_word = bits;
        first_bad_stats = s;
      }
      continue;
    }
    const u128 gap = s.rhs - five_m;
    if (gap == 0) ++equal;
    if (gap < min_gap) {
      min_gap = gap;
      min_gap_word = bits;
      min_gap_stats = s;
    }
    if (min_ratio_word.empty() ||
        checked_mul(s.rhs, ratio_den) < checked_mul(ratio_num, s.m)) {
      ratio_num = s.rhs;
      ratio_den = s.m;
      min_ratio_word = bits;
      min_ratio_stats = s;
    }
  }

  std::cout << "n=" << n << " rows=" << rows << " bad=" << bad << " equal=" << equal << "\n";
  if (!first_bad_word.empty()) {
    std::cout << "first_bad=" << first_bad_word << " H=" << first_bad_stats.H
              << " W=" << first_bad_stats.W << " m=" << show(first_bad_stats.m)
              << " rhs=" << show(first_bad_stats.rhs) << "\n";
  }
  if (!min_gap_word.empty()) {
    std::cout << "min_gap=" << show(min_gap) << " word=" << min_gap_word
              << " H=" << min_gap_stats.H << " W=" << min_gap_stats.W
              << " m=" << show(min_gap_stats.m) << "\n";
    const long double decimal = static_cast<long double>(ratio_num) /
                                static_cast<long double>(ratio_den);
    std::cout << std::setprecision(18) << "min_ratio=" << decimal
              << " fraction=" << show(ratio_num) << "/" << show(ratio_den)
              << " word=" << min_ratio_word << " H=" << min_ratio_stats.H
              << " W=" << min_ratio_stats.W << " m=" << show(min_ratio_stats.m)
              << " gap=" << show(min_ratio_stats.rhs - 5 * min_ratio_stats.m) << "\n";
    std::cout << "v_min_ratio=";
    for (std::size_t i = 0; i < min_ratio_stats.v.size(); ++i) {
      if (i != 0) std::cout << ',';
      std::cout << min_ratio_stats.v[i];
    }
    std::cout << "\n";
  }
}
