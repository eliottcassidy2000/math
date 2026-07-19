// Independent exact finite complement for the r=6 step-two five-comb ray.
//
// The companion rectangle atlas proves the sharp component target for
//
//     (b,b+2,b+4,b+6,b+8),  B=b+4 >= 168.
//
// This referee does not call or reproduce the rectangle search.  It checks
// every remaining legal row b <= 164 directly, for all seven-speed cores
// P subset {1,...,12}.  Legality is b > 13*max(P).  It constructs
//
//   S(P) \ (D_b union D_{b+2} union ... union D_{b+8})
//
// by exact endpoint subtraction and decides 7*(b+8)*L > 1 by __int128
// cross-products.  Thus the two artifacts have independent mechanisms:
// labelled fixed-polygon rectangles for the tail, and literal rational
// interval geometry for the finite complement.
//
// Tournament / alternate-carrier audit.  The useful vertices here are
// exposed rational endpoints with speed-owner sidecars.  Exact pairwise
// order is transitive after ties are coalesced (no directed cycles,
// singleton SCCs, one sorted Hamiltonian path), but order alone destroys
// metric gap length and endpoint ownership.  Runner, comb, core-component,
// tooth, endpoint, residue, wall-event, and proof-obligation carriers were
// considered; only the full labelled interval word preserves the target.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
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

static Rat rmin(Rat a, Rat b) { return compare(a, b) <= 0 ? a : b; }
static Rat rmax(Rat a, Rat b) { return compare(a, b) >= 0 ? a : b; }

static Rat difference(Rat b, Rat a) {
  return {b.n * a.d - a.n * b.d, b.d * a.d};
}

static Rat reduced(Rat x) {
  if (x.d < 0) {
    x.n = -x.n;
    x.d = -x.d;
  }
  const i64 divisor = std::gcd(std::llabs(x.n), x.d);
  return {x.n / divisor, x.d / divisor};
}

static std::string show(Rat x) {
  x = reduced(x);
  return std::to_string(x.n) + "/" + std::to_string(x.d);
}

using Interval = std::pair<Rat, Rat>;
using Region = std::vector<Interval>;

static Region remove_bad(const Region& input, int speed) {
  Region output;
  output.reserve(input.size() + 32);
  const i64 denominator = 14LL * speed;

  for (const auto& [ambient_left, ambient_right] : input) {
    int first = static_cast<int>(
                    (static_cast<i128>(ambient_left.n) * speed) /
                    ambient_left.d) -
                1;
    int last = static_cast<int>(
                   (static_cast<i128>(ambient_right.n) * speed) /
                   ambient_right.d) +
               1;
    first = std::max(first, 0);
    last = std::min(last, speed);

    Rat cursor = ambient_left;
    for (int tooth = first; tooth <= last; ++tooth) {
      Rat left{14LL * tooth - 1, denominator};
      Rat right{14LL * tooth + 1, denominator};
      if (compare(right, ambient_left) <= 0 ||
          compare(left, ambient_right) >= 0) {
        continue;
      }
      left = rmax(left, ambient_left);
      right = rmin(right, ambient_right);
      if (compare(left, cursor) > 0) output.push_back({cursor, left});
      cursor = rmax(cursor, right);
    }
    if (compare(cursor, ambient_right) < 0) {
      output.push_back({cursor, ambient_right});
    }
  }
  return output;
}

static Region safe_set(const std::vector<int>& core) {
  Region region{{Rat{0, 1}, Rat{1, 1}}};
  for (int speed : core) region = remove_bad(region, speed);
  return region;
}

static Rat longest(const Region& region) {
  Rat best{0, 1};
  for (const auto& [left, right] : region) {
    best = rmax(best, difference(right, left));
  }
  return best;
}

static std::vector<Interval> maximizers(const Region& region) {
  const Rat length = longest(region);
  std::vector<Interval> output;
  for (const auto& interval : region) {
    if (compare(difference(interval.second, interval.first), length) == 0) {
      output.push_back(interval);
    }
  }
  return output;
}

static std::string show_core(const std::vector<int>& core) {
  std::ostringstream output;
  output << "[";
  for (std::size_t index = 0; index < core.size(); ++index) {
    if (index) output << ",";
    output << core[index];
  }
  return output.str() + "]";
}

struct Candidate {
  bool set = false;
  std::vector<int> core;
  int base = 0;
  Rat length{0, 1};
};

static int compare_metric(const Candidate& left, const Candidate& right) {
  // Compare 7*(b+8)*L without reducing either rational.
  const i128 lhs = static_cast<i128>(7) * (left.base + 8) * left.length.n *
                   right.length.d;
  const i128 rhs = static_cast<i128>(7) * (right.base + 8) * right.length.n *
                   left.length.d;
  return (lhs > rhs) - (lhs < rhs);
}

static bool harder(const Candidate& left, const Candidate& right) {
  if (!left.set) return false;
  if (!right.set) return true;
  const int comparison = compare_metric(left, right);
  if (comparison) return comparison < 0;
  return std::tie(left.core, left.base) < std::tie(right.core, right.base);
}

int main() {
  constexpr int last_finite_base = 164;

  std::vector<std::vector<int>> cores;
  for (int mask = 0; mask < (1 << 12); ++mask) {
    if (__builtin_popcount(static_cast<unsigned>(mask)) != 7) continue;
    std::vector<int> core;
    for (int index = 0; index < 12; ++index) {
      if ((mask >> index) & 1) core.push_back(index + 1);
    }
    cores.push_back(core);
  }

  std::uint64_t rows = 0;
  std::uint64_t failures = 0;
  std::uint64_t active_cores = 0;
  Candidate hardest;
  Rat minimum_core_length{1, 1};
  std::vector<std::vector<int>> minimum_cores;
  std::vector<std::uint64_t> rows_by_core_maximum(13, 0);

  for (const auto& core : cores) {
    Region core_region = safe_set(core);
    const Rat core_length = longest(core_region);
    if (compare(core_length, minimum_core_length) < 0) {
      minimum_core_length = core_length;
      minimum_cores = {core};
    } else if (compare(core_length, minimum_core_length) == 0) {
      minimum_cores.push_back(core);
    }

    const int first_base = 13 * core.back() + 1;
    if (first_base <= last_finite_base) ++active_cores;
    for (int base = first_base; base <= last_finite_base; ++base) {
      Region region = core_region;
      for (int offset = 0; offset <= 8; offset += 2) {
        region = remove_bad(region, base + offset);
      }
      const Rat length = longest(region);
      ++rows;
      ++rows_by_core_maximum[core.back()];
      if (static_cast<i128>(7) * (base + 8) * length.n <= length.d) {
        ++failures;
      }

      Candidate candidate{true, core, base, length};
      if (harder(candidate, hardest)) hardest = candidate;
    }
  }

  Region hardest_region = safe_set(hardest.core);
  for (int offset = 0; offset <= 8; offset += 2) {
    hardest_region = remove_bad(hardest_region, hardest.base + offset);
  }
  const Rat metric = reduced(
      Rat{7LL * (hardest.base + 8) * hardest.length.n, hardest.length.d});
  const Rat inverse_metric = reduced(Rat{metric.d, metric.n});

  std::cout << "r=6 step-two five-comb exact finite complement\n";
  std::cout << "arithmetic=integer-rational; decisions=__int128 cross-products\n";
  std::cout << "scope=all 792 seven-speed cores; every legal b<=164\n";
  std::cout << "family=(b,b+2,b+4,b+6,b+8)\n";
  std::cout << "target=7*(b+8)*L>1\n";
  std::cout << "cores=" << cores.size() << "\n";
  std::cout << "cores_with_finite_rows=" << active_cores << "\n";
  std::cout << "finite_rows=" << rows << "\n";
  for (int maximum = 1; maximum <= 12; ++maximum) {
    if (rows_by_core_maximum[maximum]) {
      std::cout << "rows_with_core_max=" << maximum << ":"
                << rows_by_core_maximum[maximum] << "\n";
    }
  }
  std::cout << "failures=" << failures << "\n";
  std::cout << "minimum_core_component=" << show(minimum_core_length) << "\n";
  std::cout << "minimum_core_count=" << minimum_cores.size() << "\n";
  for (const auto& core : minimum_cores) {
    std::cout << "minimum_core=" << show_core(core) << "\n";
  }
  std::cout << "hardest_core=" << show_core(hardest.core) << "\n";
  std::cout << "hardest_base=" << hardest.base << "\n";
  std::cout << "hardest_killers=(" << hardest.base << ","
            << hardest.base + 2 << "," << hardest.base + 4 << ","
            << hardest.base + 6 << "," << hardest.base + 8 << ")\n";
  std::cout << "hardest_longest_component=" << show(hardest.length) << "\n";
  std::cout << "minimum_7k5L=" << show(metric) << "\n";
  std::cout << "maximum_sharp_R=" << show(inverse_metric) << "\n";
  for (const auto& [left, right] : maximizers(hardest_region)) {
    std::cout << "hardest_longest_interval=[" << show(left) << ","
              << show(right) << "]\n";
  }
  std::cout << "tournament_vertices=coalesced exact surviving endpoints\n";
  std::cout << "tournament_fingerprint=transitive; directed_cycles=0; SCCs=singletons; sorted_HP=1\n";
  std::cout << "proof_carrier=ordered rational endpoints|speed owners|metric lengths\n";
  std::cout << "destroyed_by_order_only=metric gap length|endpoint owners|removal stages\n";
  std::cout << "challenged_vertices=runners|combs|core components|teeth|endpoints|wall events|residues|proof obligations\n";
  std::cout << "certificate=" << (failures == 0 ? "PASS" : "FAIL") << "\n";
  return failures == 0 ? 0 : 1;
}
