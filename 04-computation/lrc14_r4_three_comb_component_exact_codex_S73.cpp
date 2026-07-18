// Exact endpoint certificate for THM-1097.
//
// For every nine-speed core P in {1,...,12} and every ordered triple
//
//   13 max(P) < k1 < k2 < k3,
//
// not already covered by the analytic tail, construct the exact rational
// components of
//
//   S(P) \ (D_k1 union D_k2 union D_k3)
//
// and verify 7*k3*L > 1.  The sharp periodic-comb discrepancy then prevents
// any fourth killer k4>k3 from covering that component.
//
// All acceptance decisions, endpoint comparisons, and cutoff floors use
// integer arithmetic.  No floating-point arithmetic is used.

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

// Region endpoints have denominator at most 14*642 in the guarded bank, so
// endpoint-difference products fit comfortably in signed i64.
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

// Remove the closed radius-1/(14k) teeth centred at j/k.  Open versus closed
// endpoint convention changes no component length, and the final strict
// inequality supplies an interior safe point.
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
  int k3 = 0;
};

// Compare 7*k3*L; smaller is harder.
static int compare_metric(const Candidate& a, const Candidate& b) {
  const i128 lhs = static_cast<i128>(7) * a.k3 * a.length.n * b.length.d;
  const i128 rhs = static_cast<i128>(7) * b.k3 * b.length.n * a.length.d;
  return (lhs > rhs) - (lhs < rhs);
}

static bool harder(const Candidate& a, const Candidate& b) {
  if (!a.initialized) return false;
  if (!b.initialized) return true;
  const int c = compare_metric(a, b);
  if (c != 0) return c < 0;
  return std::tie(a.k1, a.k2, a.k3) <
         std::tie(b.k1, b.k2, b.k3);
}

struct CoreResult {
  std::vector<int> core;
  Rat ell{0, 1};
  std::uint64_t finite_triples = 0;
  std::uint64_t failures = 0;
  Candidate hardest;
};

static CoreResult scan_core(const std::vector<int>& core) {
  CoreResult result;
  result.core = core;
  const Region initial = safe_set(core);
  result.ell = largest_length(initial);
  const int lo = 13 * core.back() + 1;

  // Q>0 is automatic for ell*k1>=7.  Include one guard layer beyond the
  // floor of 7/ell.
  const i128 a_num = static_cast<i128>(7) * result.ell.d;
  const i128 a_den = result.ell.n;
  const int a_max = static_cast<int>(a_num / a_den) + 1;

  for (int k1 = lo; k1 <= a_max; ++k1) {
    const Region after_k1 = remove_bad(initial, k1);

    // At k3=k2, Q>0 becomes
    //   k2(14ell-6/k1) > 7ell*k1+43.
    // Q is increasing in k3, so larger k2 are wholly analytic.  Include one
    // guard layer beyond the exact floor.
    const i128 b_num = static_cast<i128>(k1) *
        (static_cast<i128>(7) * result.ell.n * k1 +
         static_cast<i128>(43) * result.ell.d);
    const i128 b_den = static_cast<i128>(14) * result.ell.n * k1 -
        static_cast<i128>(6) * result.ell.d;
    if (b_den <= 0) std::abort();
    const int b_max = static_cast<int>(b_num / b_den) + 1;

    for (int k2 = k1 + 1; k2 <= b_max; ++k2) {
      const Region after_k2 = remove_bad(after_k1, k2);

      // Solve Q>0 exactly for k3 and include one redundant guard layer:
      //
      // k3 > (7ell(k1+k2)+37)/(21ell-6/k1-6/k2).
      const i128 c_num =
          (static_cast<i128>(7) * result.ell.n * (k1 + k2) +
           static_cast<i128>(37) * result.ell.d) * k1 * k2;
      const i128 c_den = static_cast<i128>(21) * result.ell.n * k1 * k2 -
          static_cast<i128>(6) * result.ell.d * (k1 + k2);
      if (c_den <= 0) std::abort();
      const int c_max = static_cast<int>(c_num / c_den) + 1;

      for (int k3 = k2 + 1; k3 <= c_max; ++k3) {
        ++result.finite_triples;
        const Region remainder = remove_bad(after_k2, k3);
        const Rat length = largest_length(remainder);
        Candidate current{true, length, k1, k2, k3};

        if (static_cast<i128>(7) * k3 * length.n <= length.d)
          ++result.failures;
        if (harder(current, result.hardest)) result.hardest = current;
      }
    }
  }
  return result;
}

int main() {
  std::vector<std::vector<int>> cores;
  for (int mask = 0; mask < (1 << 12); ++mask) {
    if (__builtin_popcount(static_cast<unsigned>(mask)) != 9) continue;
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

  std::uint64_t total_triples = 0;
  std::uint64_t total_failures = 0;
  Candidate finite_bank_hardest;
  Rat min_ell{1, 1};
  std::vector<std::vector<int>> min_ell_cores;
  std::size_t finite_bank_core = 0;

  std::cout << "THM-1097 exact three-comb component certificate\n";
  std::cout << "arithmetic=integer-rational (__int128 cross-products)\n";
  std::cout << "target=7*k3*L>1\n";
  std::cout << "tail_mass_lower=4*ell/7-(6/49)*(1/k1+1/k2+1/k3)\n";
  std::cout << "tail_component_upper=ell*(k1+k2+k3)+31/7\n";
  std::cout << "tail_condition=21*ell*k3-7*ell*(k1+k2)-6*(k3/k1+k3/k2)-37>0\n";
  std::cout << "finite_k1_guard=floor(7/ell)+1\n";
  std::cout << "finite_k2_guard=floor(k1*(7*ell*k1+43)/(14*ell*k1-6))+1\n";
  std::cout << "finite_k3_guard=floor((7*ell*(k1+k2)+37)/(21*ell-6/k1-6/k2))+1\n";
  std::cout << "cores=" << cores.size() << "\n";

  for (std::size_t i = 0; i < results.size(); ++i) {
    const CoreResult& row = results[i];
    total_triples += row.finite_triples;
    total_failures += row.failures;
    if (compare(row.ell, min_ell) < 0) {
      min_ell = row.ell;
      min_ell_cores = {row.core};
    } else if (compare(row.ell, min_ell) == 0) {
      min_ell_cores.push_back(row.core);
    }
    if (harder(row.hardest, finite_bank_hardest)) {
      finite_bank_hardest = row.hardest;
      finite_bank_core = i;
    }
    std::cout << "core=" << show_core(row.core)
              << " ell=" << show(row.ell)
              << " finite_triples=" << row.finite_triples
              << " failures=" << row.failures;
    if (row.hardest.initialized) {
      const Rat metric{7LL * row.hardest.k3 * row.hardest.length.n,
                       row.hardest.length.d};
      std::cout << " hardest_triple=(" << row.hardest.k1 << ","
                << row.hardest.k2 << "," << row.hardest.k3 << ")"
                << " min_7k3L=" << show(metric);
    } else {
      std::cout << " hardest_triple=none min_7k3L=analytic";
    }
    std::cout << "\n";
  }

  if (!finite_bank_hardest.initialized) std::abort();

  const CoreResult& hardest_row = results[finite_bank_core];
  const Region hardest_remainder = remove_bad(
      remove_bad(remove_bad(safe_set(hardest_row.core), finite_bank_hardest.k1),
                 finite_bank_hardest.k2),
      finite_bank_hardest.k3);
  const std::vector<Interval> maximizers = longest_intervals(hardest_remainder);
  const Rat finite_bank_metric{
      7LL * finite_bank_hardest.k3 * finite_bank_hardest.length.n,
      finite_bank_hardest.length.d};
  const Rat finite_bank_ratio{finite_bank_metric.d, finite_bank_metric.n};

  std::cout << "finite_triples_total=" << total_triples << "\n";
  std::cout << "failures_total=" << total_failures << "\n";
  std::cout << "min_core_ell=" << show(min_ell) << "\n";
  for (const auto& core : min_ell_cores)
    std::cout << "min_ell_core=" << show_core(core) << "\n";
  std::cout << "finite_bank_hardest_core=" << show_core(hardest_row.core) << "\n";
  std::cout << "finite_bank_hardest_triple=(" << finite_bank_hardest.k1 << ","
            << finite_bank_hardest.k2 << "," << finite_bank_hardest.k3 << ")\n";
  std::cout << "finite_bank_longest_component="
            << show(finite_bank_hardest.length) << "\n";
  std::cout << "finite_bank_min_7k3L=" << show(finite_bank_metric) << "\n";
  std::cout << "finite_bank_max_R=" << show(finite_bank_ratio) << "\n";
  for (const auto& [a, b] : maximizers)
    std::cout << "finite_bank_longest_interval=[" << show(a) << "," << show(b)
              << "]\n";
  std::cout << "certificate=" << (total_failures == 0 ? "PASS" : "FAIL") << "\n";
  return total_failures == 0 ? 0 : 1;
}
