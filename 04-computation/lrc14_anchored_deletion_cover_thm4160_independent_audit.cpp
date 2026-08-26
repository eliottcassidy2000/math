// Independent exact integer-lattice audit for THM-4160.
//
// This implementation uses only reduced int64 rational endpoints and
// boost::multiprecision::cpp_int measure numerators.  Every int64 product is
// checked before multiplication; the observed and theoretical bounds are
// printed in the frozen certificate.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace {

constexpr int kMaxQ = 27816;
constexpr std::int64_t kEndpointBound = 14LL * kMaxQ;
constexpr std::int64_t kEndpointCrossBound = kEndpointBound * kEndpointBound;

std::int64_t observed_max_raw_numerator = 0;
std::int64_t observed_max_raw_denominator = 0;
std::int64_t observed_max_reduced_numerator = 0;
std::int64_t observed_max_reduced_denominator = 0;
std::int64_t observed_max_checked_product = 0;

void require(bool condition, const std::string& label) {
  if (!condition) throw std::runtime_error("requirement failed: " + label);
}

std::int64_t checked_product(std::int64_t left, std::int64_t right) {
  require(left >= 0 && right >= 0, "checked product is nonnegative");
  if (left != 0 && right > std::numeric_limits<std::int64_t>::max() / left) {
    throw std::overflow_error("int64 product overflow");
  }
  const std::int64_t answer = left * right;
  observed_max_checked_product = std::max(observed_max_checked_product, answer);
  return answer;
}

struct Rat {
  std::int64_t n;
  std::int64_t d;

  Rat(std::int64_t numerator = 0, std::int64_t denominator = 1) {
    require(numerator >= 0 && denominator > 0, "nonnegative rational input");
    observed_max_raw_numerator = std::max(observed_max_raw_numerator, numerator);
    observed_max_raw_denominator = std::max(observed_max_raw_denominator, denominator);
    const std::int64_t divisor = std::gcd(numerator, denominator);
    n = numerator / divisor;
    d = denominator / divisor;
    observed_max_reduced_numerator = std::max(observed_max_reduced_numerator, n);
    observed_max_reduced_denominator = std::max(observed_max_reduced_denominator, d);
  }
};

bool less_rat(const Rat& left, const Rat& right) {
  return checked_product(left.n, right.d) < checked_product(right.n, left.d);
}

bool equal_rat(const Rat& left, const Rat& right) {
  return checked_product(left.n, right.d) == checked_product(right.n, left.d);
}

Rat maximum(const Rat& left, const Rat& right) {
  return less_rat(left, right) ? right : left;
}

Rat minimum(const Rat& left, const Rat& right) {
  return less_rat(left, right) ? left : right;
}

using Interval = std::pair<Rat, Rat>;
using Bank = std::vector<Interval>;

Bank intersect_speed(const Bank& input, int speed) {
  Bank output;
  for (const auto& [left, right] : input) {
    const std::int64_t scaled_left = checked_product(speed, left.n);
    const std::int64_t scaled_right = checked_product(speed, right.n);
    std::int64_t low = scaled_left / left.d - 1;
    std::int64_t high = scaled_right / right.d + 1;
    low = std::max<std::int64_t>(0, low);
    high = std::min<std::int64_t>(speed - 1, high);
    for (std::int64_t tooth = low; tooth <= high; ++tooth) {
      const Rat tooth_left(14 * tooth + 1, 14LL * speed);
      const Rat tooth_right(14 * tooth + 13, 14LL * speed);
      const Rat a = maximum(left, tooth_left);
      const Rat b = minimum(right, tooth_right);
      if (less_rat(b, a)) continue;
      if (!output.empty() && equal_rat(output.back().second, a)) {
        output.back().second = b;
      } else {
        output.emplace_back(a, b);
      }
    }
  }
  return output;
}

Bank safe_components(const std::vector<int>& speeds) {
  Bank result{{Rat(0), Rat(1)}};
  for (const int speed : speeds) result = intersect_speed(result, speed);
  return result;
}

cpp_int gcd_small(const cpp_int& value, std::uint64_t small) {
  const std::uint64_t remainder = (value % small).convert_to<std::uint64_t>();
  return cpp_int(std::gcd(remainder, small));
}

cpp_int gcd_big(cpp_int left, cpp_int right) {
  while (right != 0) {
    cpp_int remainder = left % right;
    left = right;
    right = remainder;
  }
  return left;
}

cpp_int lcm_denominators(const std::vector<Bank>& banks) {
  cpp_int result = 1;
  for (const Bank& bank : banks) {
    for (const auto& [left, right] : bank) {
      for (const std::int64_t denominator : {left.d, right.d}) {
        result = (result / gcd_small(result, static_cast<std::uint64_t>(denominator)))
               * denominator;
      }
    }
  }
  return result;
}

cpp_int prefix_numerator(
    const Rat& value, int speed, const cpp_int& common_denominator) {
  const std::int64_t product = checked_product(speed, value.n);
  const std::int64_t whole = product / value.d;
  const std::int64_t remainder = product % value.d;
  std::int64_t partial = 0;
  if (14 * remainder <= value.d) {
    partial = 0;
  } else if (14 * remainder >= 13 * value.d) {
    partial = 12 * value.d;
  } else {
    partial = 14 * remainder - value.d;
  }
  return cpp_int(12 * whole) * common_denominator
       + cpp_int(partial) * (common_denominator / value.d);
}

int singleton_comparison(
    const Bank& bank, int speed, const cpp_int& common_denominator) {
  cpp_int safe_numerator = 0;
  for (const auto& [left, right] : bank) {
    safe_numerator += prefix_numerator(right, speed, common_denominator)
                    - prefix_numerator(left, speed, common_denominator);
  }
  // Denominator is 14*speed*common.  Comparison with 4/63 is 9*N vs 8*q*C.
  const cpp_int difference = cpp_int(9) * safe_numerator
                           - cpp_int(8) * speed * common_denominator;
  return difference == 0 ? 0 : (difference > 0 ? 1 : -1);
}

struct MeasureResult {
  int comparison;
  cpp_int numerator;
  cpp_int denominator;
};

MeasureResult exact_measure(const Bank& bank) {
  if (bank.empty()) return MeasureResult{-1, 0, 1};
  cpp_int denominator = 1;
  for (const auto& [left, right] : bank) {
    denominator = (denominator / gcd_small(denominator, left.d)) * left.d;
    denominator = (denominator / gcd_small(denominator, right.d)) * right.d;
  }
  cpp_int numerator = 0;
  for (const auto& [left, right] : bank) {
    numerator += cpp_int(right.n) * (denominator / right.d)
               - cpp_int(left.n) * (denominator / left.d);
  }
  const cpp_int difference = cpp_int(63) * numerator - cpp_int(4) * denominator;
  const cpp_int divisor = gcd_big(numerator, denominator);
  return MeasureResult{
      difference == 0 ? 0 : (difference > 0 ? 1 : -1),
      numerator / divisor,
      denominator / divisor,
  };
}

struct Candidate {
  int q;
  std::uint32_t mask;
};

struct MaximumSummary {
  std::array<std::uint64_t, 5> tried{};
  std::array<std::uint64_t, 5> survived{};
  std::uint64_t equalities = 0;
  int maximum = 0;
  std::vector<int> witness;
  std::vector<std::pair<int, int>> surviving_pairs;
};

void maximum_dfs(
    const std::vector<int>& candidates,
    std::size_t start,
    const Bank& bank,
    std::vector<int>& selected,
    MaximumSummary& summary) {
  for (std::size_t index = start; index < candidates.size(); ++index) {
    const int depth = static_cast<int>(selected.size()) + 1;
    require(depth < static_cast<int>(summary.tried.size()), "maximum depth storage");
    ++summary.tried[depth];
    Bank next_bank = intersect_speed(bank, candidates[index]);
    const MeasureResult result = exact_measure(next_bank);
    if (result.comparison < 0) continue;
    if (result.comparison == 0) ++summary.equalities;
    ++summary.survived[depth];
    selected.push_back(candidates[index]);
    if (depth > summary.maximum) {
      summary.maximum = depth;
      summary.witness = selected;
    }
    if (depth == 2) {
      summary.surviving_pairs.emplace_back(selected[0], selected[1]);
    }
    maximum_dfs(candidates, index + 1, next_bank, selected, summary);
    selected.pop_back();
  }
}

std::string vector_text(const std::vector<int>& values) {
  std::string answer = "(";
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) answer += ',';
    answer += std::to_string(values[index]);
  }
  answer += ')';
  return answer;
}

}  // namespace

int main() {
  static_assert(kEndpointCrossBound < std::numeric_limits<std::int64_t>::max());
  const std::vector<int> anchors{120, 126, 143};
  const std::vector<int> pool{
      8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
      80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
      168, 170, 176, 190, 193, 240, 252, 264, 286, 290};
  const std::set<int> anchor_set(anchors.begin(), anchors.end());
  const std::set<int> pool_set(pool.begin(), pool.end());
  std::vector<int> optional;
  for (const int value : pool) {
    if (!anchor_set.contains(value)) optional.push_back(value);
  }
  require(optional.size() == 27, "optional size");

  std::vector<Bank> banks;
  for (const int removed : optional) {
    std::vector<int> speeds;
    for (const int value : pool) {
      if (value != removed) speeds.push_back(value);
    }
    banks.push_back(safe_components(speeds));
  }
  const cpp_int common_denominator = lcm_denominators(banks);

  std::vector<Candidate> positive_rows;
  std::map<int, std::uint64_t> histogram;
  std::uint64_t singleton_equalities = 0;
  for (int q = 1; q <= kMaxQ; ++q) {
    if (pool_set.contains(q)) continue;
    std::uint32_t mask = 0;
    for (std::size_t repair = 0; repair < banks.size(); ++repair) {
      const int comparison = singleton_comparison(banks[repair], q, common_denominator);
      if (comparison >= 0) mask |= (1U << repair);
      if (comparison == 0) ++singleton_equalities;
    }
    ++histogram[std::popcount(mask)];
    if (mask != 0) positive_rows.push_back(Candidate{q, mask});
  }

  const std::map<int, std::uint64_t> expected_histogram{
      {0, 27213}, {1, 541}, {2, 8}, {3, 4}, {4, 5}, {5, 3}, {6, 1},
      {7, 1}, {8, 1}, {9, 2}, {11, 3}, {12, 1}, {15, 1}, {16, 1}, {18, 1}};
  require(histogram == expected_histogram, "complete singleton histogram");
  require(positive_rows.size() == 573, "positive singleton rows");
  require(singleton_equalities == 0, "singleton equality rows");

  std::array<std::vector<int>, 9> threshold_candidates;
  for (int size = 1; size <= 8; ++size) {
    const int needed = 9 - size;
    for (const Candidate& candidate : positive_rows) {
      if (std::popcount(candidate.mask) >= needed) {
        threshold_candidates[size].push_back(candidate.q);
      }
    }
  }
  require(threshold_candidates[1] ==
          std::vector<int>({5, 66, 182, 298, 336, 340, 380, 386, 528, 572}),
          "singleton deletion-cover candidates");

  std::vector<MaximumSummary> maxima(optional.size());
  std::uint64_t all_multi_equalities = 0;
  for (std::size_t repair = 0; repair < optional.size(); ++repair) {
    std::vector<int> candidates;
    for (const Candidate& candidate : positive_rows) {
      if ((candidate.mask >> repair) & 1U) candidates.push_back(candidate.q);
    }
    std::vector<int> selected;
    maximum_dfs(candidates, 0, banks[repair], selected, maxima[repair]);
    require(maxima[repair].survived[1] == candidates.size(),
            "singleton mask recomputation");
    all_multi_equalities += maxima[repair].equalities;
  }

  int maximum_over_repairs = 0;
  int maximum_owner = 0;
  for (std::size_t repair = 0; repair < optional.size(); ++repair) {
    if (maxima[repair].maximum > maximum_over_repairs) {
      maximum_over_repairs = maxima[repair].maximum;
      maximum_owner = optional[repair];
    }
    if (optional[repair] == 252) {
      require(maxima[repair].maximum == 2, "r252 maximum");
      require(maxima[repair].survived[2] == 45, "r252 pair count");
      require(maxima[repair].survived[3] == 0, "r252 triple count");
    } else {
      require(maxima[repair].maximum <= 1, "non-252 repair maximum");
    }
  }
  require(maximum_over_repairs == 2 && maximum_owner == 252,
          "global joint-newcomer maximum");
  require(all_multi_equalities == 0, "all recomputed equality rows");

  const std::size_t repair252_index = static_cast<std::size_t>(
      std::find(optional.begin(), optional.end(), 252) - optional.begin());
  Bank witness_bank = intersect_speed(banks[repair252_index], 5);
  witness_bank = intersect_speed(witness_bank, 66);
  const MeasureResult witness_measure = exact_measure(witness_bank);
  require(witness_measure.comparison > 0, "pair witness is strict");
  require(witness_measure.numerator == cpp_int(347277745)
          && witness_measure.denominator == cpp_int(5333672344),
          "pair witness exact measure");

  require(observed_max_raw_numerator <= kEndpointBound,
          "raw endpoint numerator bound");
  require(observed_max_raw_denominator <= kEndpointBound,
          "raw endpoint denominator bound");
  require(observed_max_reduced_numerator <= kEndpointBound,
          "reduced endpoint numerator bound");
  require(observed_max_reduced_denominator <= kEndpointBound,
          "reduced endpoint denominator bound");
  require(observed_max_checked_product <= kEndpointCrossBound,
          "checked cross-product bound");

  std::cout << "LRC14_ANCHORED_DELETION_COVER_THM4160_CPP_AUDIT_20260826\n";
  std::cout << "universe=positive q notin P;finite_scan=1..27816;"
               "infinite_tail_closed_by_B252<27817\n";
  std::cout << "pool_size=30;optional_size=27;positive_rows="
            << positive_rows.size() << ";singleton_equalities="
            << singleton_equalities << "\n";
  std::cout << "singleton_histogram=";
  for (const auto& [count, frequency] : histogram) {
    std::cout << count << ':' << frequency << ',';
  }
  std::cout << "\n";
  for (int size = 1; size <= 8; ++size) {
    std::cout << "member_threshold_for_size=" << size
              << ";needed_singleton_repairs=" << 9 - size
              << ";candidate_count=" << threshold_candidates[size].size();
    if (size <= 3) {
      std::cout << ";candidates=" << vector_text(threshold_candidates[size]);
    }
    std::cout << "\n";
  }
  for (std::size_t repair = 0; repair < optional.size(); ++repair) {
    const MaximumSummary& summary = maxima[repair];
    std::cout << "repair=" << optional[repair]
              << ";singleton_candidates=" << summary.survived[1]
              << ";pair_tests=" << summary.tried[2]
              << ";pair_survivors=" << summary.survived[2]
              << ";triple_tests=" << summary.tried[3]
              << ";triple_survivors=" << summary.survived[3]
              << ";maximum_joint_size=" << summary.maximum
              << ";witness=" << vector_text(summary.witness)
              << ";equalities=" << summary.equalities << "\n";
  }
  const MaximumSummary& summary252 = maxima[repair252_index];
  std::cout << "repair252_pair_rows=(";
  for (std::size_t index = 0; index < summary252.surviving_pairs.size(); ++index) {
    if (index != 0) std::cout << ',';
    std::cout << '(' << summary252.surviving_pairs[index].first << ','
              << summary252.surviving_pairs[index].second << ')';
  }
  std::cout << ")\n";
  std::cout << "pair_witness=(5,66);repair=252;measure="
            << witness_measure.numerator << '/' << witness_measure.denominator
            << ";strict=True\n";
  std::cout << "hierarchy=max_D_size1:18;max_D_size2:1;"
               "max_D_sizes3to8:0;required:8,7,6,5,4,3,2,1\n";
  std::cout << "integer_width_audit=max_q=" << kMaxQ
            << ";theoretical_endpoint_bound=" << kEndpointBound
            << ";theoretical_cross_product_bound=" << kEndpointCrossBound
            << ";observed_raw_numerator=" << observed_max_raw_numerator
            << ";observed_raw_denominator=" << observed_max_raw_denominator
            << ";observed_reduced_numerator=" << observed_max_reduced_numerator
            << ";observed_reduced_denominator=" << observed_max_reduced_denominator
            << ";observed_checked_product=" << observed_max_checked_product << "\n";
  std::cout << "all_exact_threshold_comparisons_use_cpp_int=True;PASS\n";
}
