// Exact finite-complement referee for THM-1133.
//
// THM-1129 assigns every pair
//
//   P subset {1,...,12}, |P|=8,  A={0<a<b<c<=30}
//
// an exact one-sided 1/7-rectangle tail.  THM-1123 owns the initial rows
// with K+c <= 13*max(P)+40.  This program reconstructs the rectangle tail
// from rational core endpoints and offset-center collisions, then checks
// every one of the exactly 3,539,936 rows lying strictly between those two
// certificates by exact endpoint subtraction.
//
// No floating point enters a decision.  Rational comparisons and the sharp
// test 7*(K+c)*L>1 use signed __int128 cross-products.
//
// Tournament audit: exact candidate-boundary order is transitive.  It loses
// the metric cyclic gap, adjacent safe side, wall slopes, and surviving
// endpoint distances.  The faithful carrier is the labelled boundary word
// with those sidecars; runner or naked wall-order tournaments are not enough.

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
using u64 = std::uint64_t;

struct Rat {
  i64 n;
  i64 d;
};

static Rat reduced(Rat x) {
  if (x.d < 0) {
    x.n = -x.n;
    x.d = -x.d;
  }
  const i64 g = std::gcd(std::llabs(x.n), x.d);
  return {x.n / g, x.d / g};
}

static int compare(Rat a, Rat b) {
  const i128 x = static_cast<i128>(a.n) * b.d;
  const i128 y = static_cast<i128>(b.n) * a.d;
  return (x > y) - (x < y);
}

static Rat rmin(Rat a, Rat b) { return compare(a, b) <= 0 ? a : b; }
static Rat rmax(Rat a, Rat b) { return compare(a, b) >= 0 ? a : b; }

static Rat add(Rat a, Rat b) {
  return reduced({a.n * b.d + b.n * a.d, a.d * b.d});
}

static Rat subtract(Rat a, Rat b) {
  return reduced({a.n * b.d - b.n * a.d, a.d * b.d});
}

static Rat divide_by(Rat a, int k) {
  return reduced({a.n, a.d * static_cast<i64>(k)});
}

static std::string show(Rat x) {
  x = reduced(x);
  return std::to_string(x.n) + "/" + std::to_string(x.d);
}

static i64 positive_mod(i128 x, i64 d) {
  i64 r = static_cast<i64>(x % d);
  if (r < 0) r += d;
  return r;
}

using Interval = std::pair<Rat, Rat>;
using Region = std::vector<Interval>;

static Region remove_bad(const Region& input, int k) {
  Region output;
  output.reserve(input.size() + 32);
  const i64 den = 14LL * k;
  for (const auto& [a, b] : input) {
    int lo = static_cast<int>((static_cast<i128>(a.n) * k) / a.d) - 1;
    int hi = static_cast<int>((static_cast<i128>(b.n) * k) / b.d) + 1;
    lo = std::max(lo, 0);
    hi = std::min(hi, k);
    Rat cursor = a;
    for (int tooth = lo; tooth <= hi; ++tooth) {
      Rat left{14LL * tooth - 1, den};
      Rat right{14LL * tooth + 1, den};
      if (compare(right, a) <= 0 || compare(left, b) >= 0) continue;
      left = rmax(left, a);
      right = rmin(right, b);
      if (compare(left, cursor) > 0) output.push_back({cursor, left});
      cursor = rmax(cursor, right);
    }
    if (compare(cursor, b) < 0) output.push_back({cursor, b});
  }
  return output;
}

static Region safe_set(const std::vector<int>& speeds) {
  Region region{{Rat{0, 1}, Rat{1, 1}}};
  for (int speed : speeds) region = remove_bad(region, speed);
  return region;
}

static Rat longest(const Region& region) {
  Rat best{0, 1};
  for (const auto& [left, right] : region) {
    const Rat length = subtract(right, left);
    if (compare(length, best) > 0) best = length;
  }
  return best;
}

static std::vector<Interval> maximizers(const Region& region) {
  const Rat best = longest(region);
  std::vector<Interval> output;
  for (const auto& interval : region) {
    if (compare(subtract(interval.second, interval.first), best) == 0) {
      output.push_back(interval);
    }
  }
  return output;
}

static std::string show_core(const std::vector<int>& core) {
  std::ostringstream out;
  out << "[";
  for (std::size_t i = 0; i < core.size(); ++i) {
    if (i) out << ",";
    out << core[i];
  }
  return out.str() + "]";
}

struct Shape {
  int a;
  int b;
  int c;
};

static std::vector<Rat> core_endpoints() {
  std::vector<Rat> points{{0, 1}, {1, 1}};
  for (int p = 1; p <= 12; ++p) {
    for (int tooth = -1; tooth <= p + 1; ++tooth) {
      for (int sign : {-1, 1}) {
        Rat t = reduced({14LL * tooth + sign, 14LL * p});
        if (compare(t, Rat{0, 1}) >= 0 && compare(t, Rat{1, 1}) <= 0) {
          points.push_back(t);
        }
      }
    }
  }
  std::sort(points.begin(), points.end(),
            [](Rat x, Rat y) { return compare(x, y) < 0; });
  points.erase(std::unique(points.begin(), points.end(),
                           [](Rat x, Rat y) { return compare(x, y) == 0; }),
               points.end());
  return points;
}

static std::vector<Rat> candidate_points(const Shape& shape,
                                         const std::vector<Rat>& base) {
  const int offsets[4] = {0, shape.a, shape.b, shape.c};
  std::vector<Rat> points = base;
  for (int i = 0; i < 4; ++i) {
    for (int j = i + 1; j < 4; ++j) {
      const int difference = offsets[j] - offsets[i];
      for (int numerator = 0; numerator <= difference; ++numerator) {
        points.push_back(reduced({numerator, difference}));
      }
    }
  }
  std::sort(points.begin(), points.end(),
            [](Rat x, Rat y) { return compare(x, y) < 0; });
  points.erase(std::unique(points.begin(), points.end(),
                           [](Rat x, Rat y) { return compare(x, y) == 0; }),
               points.end());
  return points;
}

static unsigned safe_speed_mask(Rat t) {
  unsigned mask = 0;
  for (int p = 1; p <= 12; ++p) {
    const i64 residue = positive_mod(static_cast<i128>(p) * t.n, t.d);
    const i64 distance = std::min(residue, t.d - residue);
    if (14 * distance >= t.d) mask |= 1U << (p - 1);
  }
  return mask;
}

static Rat midpoint(Rat left, Rat right) {
  return divide_by(add(left, right), 2);
}

static Rat vertical_safe_width(const Shape& shape, Rat t) {
  const int offsets[4] = {0, shape.a, shape.b, shape.c};
  i64 centers[4];
  for (int i = 0; i < 4; ++i) {
    centers[i] = positive_mod(-static_cast<i128>(offsets[i]) * t.n, t.d);
  }
  std::sort(centers, centers + 4);
  i64 maximum_gap = t.d - centers[3] + centers[0];
  for (int i = 0; i < 3; ++i) {
    maximum_gap = std::max(maximum_gap, centers[i + 1] - centers[i]);
  }
  return reduced({7 * maximum_gap - t.d, 7 * t.d});
}

static int tail_start(Rat time_width) {
  // ceil((8/7)/time_width)
  const i128 numerator = static_cast<i128>(8) * time_width.d;
  const i128 denominator = static_cast<i128>(7) * time_width.n;
  return static_cast<int>((numerator + denominator - 1) / denominator);
}

static u64 fnv_mix(u64 hash, u64 value) {
  constexpr u64 prime = 1099511628211ULL;
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (8 * byte)) & 0xffU;
    hash *= prime;
  }
  return hash;
}

struct Candidate {
  bool set = false;
  Rat length{0, 1};
  int core_index = 0;
  Shape shape{0, 0, 0};
  int K = 0;
};

static int compare_metric(const Candidate& x, const Candidate& y) {
  const i128 lhs = static_cast<i128>(7) * (x.K + x.shape.c) * x.length.n *
                   y.length.d;
  const i128 rhs = static_cast<i128>(7) * (y.K + y.shape.c) * y.length.n *
                   x.length.d;
  return (lhs > rhs) - (lhs < rhs);
}

static bool harder(const Candidate& x, const Candidate& y) {
  if (!x.set) return false;
  if (!y.set) return true;
  const int metric = compare_metric(x, y);
  if (metric) return metric < 0;
  return std::tie(x.core_index, x.shape.a, x.shape.b, x.shape.c, x.K) <
         std::tie(y.core_index, y.shape.a, y.shape.b, y.shape.c, y.K);
}

struct ShapeResult {
  u64 pairs = 0;
  u64 closed_at_legal_start = 0;
  u64 raw_finite = 0;
  u64 bottom = 0;
  u64 residual = 0;
  u64 tested = 0;
  u64 failures = 0;
  int maximum_tail = 0;
  Candidate hardest;
  Candidate first_failure;
  u64 digest = 1469598103934665603ULL;
};

static ShapeResult scan_shape(
    const Shape& shape, const std::vector<Rat>& base,
    const std::vector<std::vector<int>>& eligible,
    const std::vector<std::vector<int>>& cores,
    const std::vector<Region>& core_regions) {
  ShapeResult result;
  const std::vector<Rat> points = candidate_points(shape, base);
  std::vector<Rat> widths;
  widths.reserve(points.size());
  for (Rat t : points) widths.push_back(vertical_safe_width(shape, t));

  std::vector<Rat> best(cores.size(), Rat{0, 1});
  for (std::size_t cell = 0; cell + 1 < points.size(); ++cell) {
    const Rat left = points[cell];
    const Rat right = points[cell + 1];
    const auto& core_indices = eligible[safe_speed_mask(midpoint(left, right))];
    if (core_indices.empty()) continue;
    const Rat cell_width = subtract(right, left);
    for (std::size_t endpoint : {cell, cell + 1}) {
      const Rat width = widths[endpoint];
      if (compare(width, Rat{2, 5}) < 0) continue;
      const Rat drift_time = divide_by(subtract(width, Rat{1, 7}), shape.c);
      const Rat certificate = rmin(cell_width, drift_time);
      for (int core_index : core_indices) {
        if (compare(certificate, best[core_index]) > 0) {
          best[core_index] = certificate;
        }
      }
    }
  }

  for (std::size_t core_index = 0; core_index < cores.size(); ++core_index) {
    ++result.pairs;
    if (compare(best[core_index], Rat{0, 1}) <= 0) {
      std::cerr << "missing rectangle certificate at core=" << core_index
                << " shape=(0," << shape.a << "," << shape.b << ","
                << shape.c << ")\n";
      std::abort();
    }
    const int tail = tail_start(best[core_index]);
    result.maximum_tail = std::max(result.maximum_tail, tail);
    const int legal = 13 * cores[core_index].back() + 1;
    if (tail <= legal) ++result.closed_at_legal_start;
    const int finite = std::max(0, tail - legal);
    const int bottom = std::min(finite, 40 - shape.c);
    const int residual = finite - bottom;
    result.raw_finite += finite;
    result.bottom += bottom;
    result.residual += residual;

    const int first_K = legal + (40 - shape.c);
    for (int K = first_K; K < tail; ++K) {
      ++result.tested;
      Region region = core_regions[core_index];
      region = remove_bad(region, K);
      region = remove_bad(region, K + shape.a);
      region = remove_bad(region, K + shape.b);
      region = remove_bad(region, K + shape.c);
      const Rat length = longest(region);
      Candidate candidate{true, length, static_cast<int>(core_index), shape, K};
      if (static_cast<i128>(7) * (K + shape.c) * length.n <= length.d) {
        ++result.failures;
        if (!result.first_failure.set) result.first_failure = candidate;
      }
      if (harder(candidate, result.hardest)) result.hardest = candidate;

      const Rat normalized = reduced(length);
      result.digest = fnv_mix(result.digest, core_index);
      result.digest = fnv_mix(result.digest, shape.a);
      result.digest = fnv_mix(result.digest, shape.b);
      result.digest = fnv_mix(result.digest, shape.c);
      result.digest = fnv_mix(result.digest, K);
      result.digest = fnv_mix(result.digest, static_cast<u64>(normalized.n));
      result.digest = fnv_mix(result.digest, static_cast<u64>(normalized.d));
    }
  }
  return result;
}

int main() {
  std::vector<std::vector<int>> cores;
  for (int mask = 0; mask < (1 << 12); ++mask) {
    if (__builtin_popcount(static_cast<unsigned>(mask)) != 8) continue;
    std::vector<int> core;
    for (int bit = 0; bit < 12; ++bit) {
      if ((mask >> bit) & 1) core.push_back(bit + 1);
    }
    cores.push_back(core);
  }
  if (cores.size() != 495) return 2;

  std::vector<unsigned> core_masks;
  std::vector<Region> core_regions;
  for (const auto& core : cores) {
    unsigned mask = 0;
    for (int p : core) mask |= 1U << (p - 1);
    core_masks.push_back(mask);
    core_regions.push_back(safe_set(core));
  }
  std::vector<std::vector<int>> eligible(1 << 12);
  for (unsigned safe = 0; safe < (1U << 12); ++safe) {
    for (std::size_t core = 0; core < core_masks.size(); ++core) {
      if ((core_masks[core] & ~safe) == 0) {
        eligible[safe].push_back(static_cast<int>(core));
      }
    }
  }

  std::vector<Shape> shapes;
  for (int a = 1; a <= 28; ++a) {
    for (int b = a + 1; b <= 29; ++b) {
      for (int c = b + 1; c <= 30; ++c) shapes.push_back({a, b, c});
    }
  }
  if (shapes.size() != 4060) return 3;

  const std::vector<Rat> base = core_endpoints();
  std::vector<ShapeResult> results(shapes.size());
  std::atomic<std::size_t> next{0};
  const unsigned detected = std::thread::hardware_concurrency();
  const unsigned workers = std::max(1U, std::min(8U, detected ? detected : 8U));
  std::vector<std::thread> pool;
  for (unsigned worker = 0; worker < workers; ++worker) {
    pool.emplace_back([&] {
      for (;;) {
        const std::size_t index = next.fetch_add(1);
        if (index >= shapes.size()) return;
        results[index] =
            scan_shape(shapes[index], base, eligible, cores, core_regions);
      }
    });
  }
  for (auto& thread : pool) thread.join();

  ShapeResult total;
  total.digest = 1469598103934665603ULL;
  for (const auto& result : results) {
    total.pairs += result.pairs;
    total.closed_at_legal_start += result.closed_at_legal_start;
    total.raw_finite += result.raw_finite;
    total.bottom += result.bottom;
    total.residual += result.residual;
    total.tested += result.tested;
    total.failures += result.failures;
    total.maximum_tail = std::max(total.maximum_tail, result.maximum_tail);
    if (harder(result.hardest, total.hardest)) total.hardest = result.hardest;
    if (!total.first_failure.set && result.first_failure.set) {
      total.first_failure = result.first_failure;
    }
    total.digest = fnv_mix(total.digest, result.digest);
    total.digest = fnv_mix(total.digest, result.tested);
  }

  if (total.pairs != 2'009'700ULL ||
      total.closed_at_legal_start != 1'802'872ULL ||
      total.raw_finite != 6'040'056ULL || total.bottom != 2'500'120ULL ||
      total.residual != 3'539'936ULL || total.tested != total.residual ||
      total.maximum_tail != 832) {
    std::cerr << "THM-1129 accounting mismatch\n";
    return 4;
  }

  // Independent shared-row guard against THM-1123's stored extremal.
  const std::vector<int> shared_core{1, 2, 4, 5, 7, 9, 11, 12};
  Region shared = safe_set(shared_core);
  for (int speed : {158, 160, 162, 164}) shared = remove_bad(shared, speed);
  const Rat shared_length = reduced(longest(shared));
  if (compare(shared_length, Rat{41, 25920}) != 0) return 5;

  const Candidate& hardest = total.hardest;
  Region hard_region = core_regions[hardest.core_index];
  for (int speed : {hardest.K, hardest.K + hardest.shape.a,
                    hardest.K + hardest.shape.b,
                    hardest.K + hardest.shape.c}) {
    hard_region = remove_bad(hard_region, speed);
  }
  const Rat hard_metric = reduced(
      {7LL * (hardest.K + hardest.shape.c) * hardest.length.n,
       hardest.length.d});

  std::cout << "THM-1133 bounded-offset all-scale finite complement\n";
  std::cout << "arithmetic=integer-rational; decisions=__int128 cross-products\n";
  std::cout << "cores=" << cores.size() << "\n";
  std::cout << "offset_shapes=" << shapes.size() << "\n";
  std::cout << "core_shape_pairs=" << total.pairs << "\n";
  std::cout << "maximum_individual_tail_K=" << total.maximum_tail << "\n";
  std::cout << "pairs_tail_closed_at_legal_start="
            << total.closed_at_legal_start << "\n";
  std::cout << "per_pair_raw_finite_rows=" << total.raw_finite << "\n";
  std::cout << "rows_already_in_THM1123_bottom_bank=" << total.bottom << "\n";
  std::cout << "residual_rows_expected=" << total.residual << "\n";
  std::cout << "residual_rows_tested=" << total.tested << "\n";
  std::cout << "failures_total=" << total.failures << "\n";
  std::cout << "row_digest_fnv64=" << total.digest << "\n";
  std::cout << "THM1123_shared_row_L=" << show(shared_length) << "\n";
  std::cout << "THM1123_shared_row_7k4L="
            << show(reduced({7LL * 164 * shared_length.n, shared_length.d}))
            << "\n";
  std::cout << "hardest_core=" << show_core(cores[hardest.core_index]) << "\n";
  std::cout << "hardest_offsets=(0," << hardest.shape.a << ","
            << hardest.shape.b << "," << hardest.shape.c << ")\n";
  std::cout << "hardest_K=" << hardest.K << "\n";
  std::cout << "hardest_killers=(" << hardest.K << ","
            << hardest.K + hardest.shape.a << ","
            << hardest.K + hardest.shape.b << ","
            << hardest.K + hardest.shape.c << ")\n";
  std::cout << "hardest_longest_component=" << show(hardest.length) << "\n";
  std::cout << "min_7k4L=" << show(hard_metric) << "\n";
  for (const auto& [left, right] : maximizers(hard_region)) {
    std::cout << "hardest_longest_interval=[" << show(left) << ","
              << show(right) << "]\n";
  }
  if (total.first_failure.set) {
    const auto& failure = total.first_failure;
    std::cout << "first_failure_core=" << show_core(cores[failure.core_index])
              << "\n";
    std::cout << "first_failure_offsets=(0," << failure.shape.a << ","
              << failure.shape.b << "," << failure.shape.c << ")\n";
    std::cout << "first_failure_K=" << failure.K << "\n";
  }
  std::cout << "tournament_vertices=labelled candidate boundaries and surviving endpoints\n";
  std::cout << "pairwise_observable=exact rational order; tournament=transitive\n";
  std::cout << "order_only_destroyed=metric gaps|owners|wall slopes|safe side|stage\n";
  std::cout << "challenged_vertices=runners|combs|core cells|boundaries|wall events|residues|proof obligations\n";
  std::cout << "scope=all legal K on every normalized offset shape of span <=30\n";
  std::cout << "certificate=" << (total.failures == 0 ? "PASS" : "FAIL")
            << "\n";
  return total.failures == 0 ? 0 : 1;
}
