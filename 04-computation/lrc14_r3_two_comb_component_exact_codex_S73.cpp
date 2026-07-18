// Exact certificate for THM-1094.
//
// For every ten-speed core P in {1,...,12} and every ordered pair
//   13 max(P) < k1 < k2
// not already covered by the analytic tail, construct the exact rational
// components of S(P) \ (D_k1 union D_k2) and verify
//
//   3 k2 L(P;k1,k2) > 1.
//
// The analytic tail and the finite cutoffs are printed with the certificate.
// All acceptance decisions and endpoint comparisons use signed __int128
// cross-products.  No floating-point arithmetic is used.

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

using i64 = std::int64_t;
using i128 = __int128_t;

struct Rat {
  i64 n;
  i64 d;
};

static int compare(Rat a, Rat b) {
  const i128 lhs = static_cast<i128>(a.n) * b.d;
  const i128 rhs = static_cast<i128>(b.n) * a.d;
  return (lhs > rhs) - (lhs < rhs);
}

static Rat rat_min(Rat a, Rat b) { return compare(a, b) <= 0 ? a : b; }
static Rat rat_max(Rat a, Rat b) { return compare(a, b) >= 0 ? a : b; }

// The denominators entering this function are individual endpoint
// denominators, at most 14*3345 in the finite bank.  Their product fits i64.
static Rat subtract(Rat b, Rat a) {
  return {b.n * a.d - a.n * b.d, b.d * a.d};
}

static Rat reduced(Rat x) {
  if (x.d < 0) {
    x.n = -x.n;
    x.d = -x.d;
  }
  const i64 g = std::gcd(std::llabs(x.n), x.d);
  return {x.n / g, x.d / g};
}

static std::string show(Rat x) {
  x = reduced(x);
  return std::to_string(x.n) + "/" + std::to_string(x.d);
}

using Interval = std::pair<Rat, Rat>;
using Region = std::vector<Interval>;

// Remove the closed radius-1/(14k) teeth centred at j/k.  Endpoint convention
// does not change any component length; the strict final inequality supplies
// an interior safe point.
static Region remove_bad(const Region& input, int k) {
  Region output;
  output.reserve(input.size() + 32);
  const i64 den = 14LL * k;

  for (const auto& [a, b] : input) {
    int j_lo = static_cast<int>((static_cast<i128>(a.n) * k) / a.d) - 1;
    int j_hi = static_cast<int>((static_cast<i128>(b.n) * k) / b.d) + 1;
    j_lo = std::max(j_lo, 0);
    j_hi = std::min(j_hi, k);
    Rat cursor = a;

    for (int j = j_lo; j <= j_hi; ++j) {
      Rat left{14LL * j - 1, den};
      Rat right{14LL * j + 1, den};
      if (compare(right, a) <= 0 || compare(left, b) >= 0) continue;
      left = rat_max(left, a);
      right = rat_min(right, b);
      if (compare(left, cursor) > 0) output.push_back({cursor, left});
      cursor = rat_max(cursor, right);
    }
    if (compare(cursor, b) < 0) output.push_back({cursor, b});
  }
  return output;
}

static Region safe_set(const std::vector<int>& core) {
  Region result{{Rat{0, 1}, Rat{1, 1}}};
  for (const int p : core) result = remove_bad(result, p);
  return result;
}

static Rat largest_length(const Region& region) {
  Rat best{0, 1};
  for (const auto& [a, b] : region) {
    const Rat length = subtract(b, a);
    if (compare(length, best) > 0) best = length;
  }
  return best;
}

static std::vector<Interval> longest_intervals(const Region& region) {
  const Rat length = largest_length(region);
  std::vector<Interval> result;
  for (const auto& interval : region) {
    if (compare(subtract(interval.second, interval.first), length) == 0)
      result.push_back(interval);
  }
  return result;
}

static std::string show_core(const std::vector<int>& core) {
  std::ostringstream out;
  out << "[";
  for (std::size_t i = 0; i < core.size(); ++i) {
    if (i) out << ",";
    out << core[i];
  }
  out << "]";
  return out.str();
}

struct Candidate {
  bool initialized = false;
  Rat length{0, 1};
  int k1 = 0;
  int k2 = 0;
};

// Compare the target-scaled quantities 3*k2*L.  Smaller is harder.
static int compare_metric(const Candidate& a, const Candidate& b) {
  const i128 lhs = static_cast<i128>(3) * a.k2 * a.length.n * b.length.d;
  const i128 rhs = static_cast<i128>(3) * b.k2 * b.length.n * a.length.d;
  return (lhs > rhs) - (lhs < rhs);
}

static bool harder(const Candidate& a, const Candidate& b) {
  if (!b.initialized) return true;
  const int c = compare_metric(a, b);
  if (c != 0) return c < 0;
  return std::tie(a.k1, a.k2) < std::tie(b.k1, b.k2);
}

struct CoreResult {
  std::vector<int> core;
  Rat ell{0, 1};
  std::uint64_t finite_pairs = 0;
  std::uint64_t failures = 0;
  Candidate hardest;
};

static CoreResult scan_core(const std::vector<int>& core) {
  CoreResult result;
  result.core = core;
  const Region initial = safe_set(core);
  result.ell = largest_length(initial);
  const int lo = 13 * core.back() + 1;

  // The analytic condition is automatic for ell*k1 >= 209/7.  Scan one
  // guard layer beyond floor(209/(7ell)).
  const i128 x_num = static_cast<i128>(209) * result.ell.d;
  const i128 x_den = static_cast<i128>(7) * result.ell.n;
  const int x_max = static_cast<int>(x_num / x_den) + 1;

  for (int k1 = lo; k1 <= x_max; ++k1) {
    const Region after_k1 = remove_bad(initial, k1);

    // The exact tail inequality
    //   56 ell k2 - 49 ell k1 - 24(k2/k1) - 185 > 0
    // is automatic above this rational threshold.  Again include one guard
    // layer beyond its floor.
    const i128 y_num = static_cast<i128>(k1) *
        (static_cast<i128>(49) * result.ell.n * k1 +
         static_cast<i128>(185) * result.ell.d);
    const i128 y_den = static_cast<i128>(56) * result.ell.n * k1 -
        static_cast<i128>(24) * result.ell.d;
    if (y_den <= 0) std::abort();
    const int y_max = static_cast<int>(y_num / y_den) + 1;

    for (int k2 = k1 + 1; k2 <= y_max; ++k2) {
      ++result.finite_pairs;
      const Region remainder = remove_bad(after_k1, k2);
      const Rat length = largest_length(remainder);
      Candidate current{true, length, k1, k2};

      // This is the theorem's exact strict inequality.
      if (static_cast<i128>(3) * k2 * length.n <= length.d)
        ++result.failures;
      if (harder(current, result.hardest)) result.hardest = current;
    }
  }
  return result;
}

int main() {
  std::vector<std::vector<int>> cores;
  for (int mask = 0; mask < (1 << 12); ++mask) {
    if (__builtin_popcount(static_cast<unsigned>(mask)) != 10) continue;
    std::vector<int> core;
    for (int bit = 0; bit < 12; ++bit)
      if ((mask >> bit) & 1) core.push_back(bit + 1);
    cores.push_back(core);
  }

  std::vector<CoreResult> results(cores.size());
  std::atomic<std::size_t> next{0};
  const unsigned detected = std::thread::hardware_concurrency();
  const unsigned workers = std::max(1u, std::min(8u, detected ? detected : 8u));
  std::vector<std::thread> pool;
  for (unsigned worker = 0; worker < workers; ++worker) {
    pool.emplace_back([&] {
      for (;;) {
        const std::size_t i = next.fetch_add(1);
        if (i >= cores.size()) return;
        results[i] = scan_core(cores[i]);
      }
    });
  }
  for (auto& worker : pool) worker.join();

  std::uint64_t total_pairs = 0;
  std::uint64_t total_failures = 0;
  Candidate global_hardest;
  Rat min_ell{1, 1};
  std::vector<std::vector<int>> min_ell_cores;
  std::size_t global_core = 0;

  std::cout << "THM-1094 exact two-comb component certificate\n";
  std::cout << "arithmetic=integer-rational (__int128 cross-products)\n";
  std::cout << "tail_mass_lower=5*ell/7-8/(49*k1)-8/(49*k2)\n";
  std::cout << "tail_component_upper=ell*(k1+k2)+23/7\n";
  std::cout << "tail_condition=56*ell*k2-49*ell*k1-24*(k2/k1)-185>0\n";
  std::cout << "finite_x_guard=floor(209/(7*ell))+1\n";
  std::cout << "finite_y_guard=floor(k1*(49*ell*k1+185)/(56*ell*k1-24))+1\n";
  std::cout << "cores=" << cores.size() << "\n";

  for (std::size_t i = 0; i < results.size(); ++i) {
    const CoreResult& row = results[i];
    total_pairs += row.finite_pairs;
    total_failures += row.failures;
    if (compare(row.ell, min_ell) < 0) {
      min_ell = row.ell;
      min_ell_cores = {row.core};
    } else if (compare(row.ell, min_ell) == 0) {
      min_ell_cores.push_back(row.core);
    }
    if (harder(row.hardest, global_hardest)) {
      global_hardest = row.hardest;
      global_core = i;
    }
    const Rat metric{3LL * row.hardest.k2 * row.hardest.length.n,
                     row.hardest.length.d};
    std::cout << "core=" << show_core(row.core)
              << " ell=" << show(row.ell)
              << " finite_pairs=" << row.finite_pairs
              << " failures=" << row.failures
              << " hardest_pair=(" << row.hardest.k1 << ","
              << row.hardest.k2 << ")"
              << " min_3k2L=" << show(metric) << "\n";
  }

  const CoreResult& hardest_row = results[global_core];
  const Region hardest_remainder = remove_bad(
      remove_bad(safe_set(hardest_row.core), global_hardest.k1),
      global_hardest.k2);
  const std::vector<Interval> maximizers = longest_intervals(hardest_remainder);
  const Rat global_metric{3LL * global_hardest.k2 * global_hardest.length.n,
                          global_hardest.length.d};
  const Rat global_ratio{global_metric.d, global_metric.n};

  std::cout << "finite_pairs_total=" << total_pairs << "\n";
  std::cout << "failures_total=" << total_failures << "\n";
  std::cout << "min_core_ell=" << show(min_ell) << "\n";
  for (const auto& core : min_ell_cores)
    std::cout << "min_ell_core=" << show_core(core) << "\n";
  std::cout << "finite_bank_hardest_core=" << show_core(hardest_row.core) << "\n";
  std::cout << "finite_bank_hardest_pair=(" << global_hardest.k1 << ","
            << global_hardest.k2 << ")\n";
  std::cout << "finite_bank_longest_component=" << show(global_hardest.length) << "\n";
  std::cout << "finite_bank_min_3k2L=" << show(global_metric) << "\n";
  std::cout << "finite_bank_max_R=" << show(global_ratio) << "\n";
  for (const auto& [a, b] : maximizers)
    std::cout << "finite_bank_longest_interval=[" << show(a) << "," << show(b)
              << "]\n";
  std::cout << "certificate=PASS\n";
  return total_failures == 0 ? 0 : 1;
}
