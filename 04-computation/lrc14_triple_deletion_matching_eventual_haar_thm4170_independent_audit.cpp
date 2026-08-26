#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using i128 = __int128_t;
using u128 = __uint128_t;
constexpr std::int64_t kCommon = 18241159416480LL;

void require(bool condition, const std::string& label) {
  if (!condition) throw std::runtime_error("requirement failed: " + label);
}

struct Rat {
  std::int64_t n;
  std::int64_t d;
  Rat(std::int64_t numerator = 0, std::int64_t denominator = 1) {
    const std::int64_t divisor = std::gcd(numerator, denominator);
    n = numerator / divisor;
    d = denominator / divisor;
  }
};

struct RatLess {
  bool operator()(const Rat& left, const Rat& right) const {
    return i128(left.n) * right.d < i128(right.n) * left.d;
  }
};

Rat midpoint(const Rat& left, const Rat& right) {
  return Rat(left.n * right.d + right.n * left.d, 2 * left.d * right.d);
}

bool safe_at(const Rat& point, int speed) {
  const std::int64_t residue = static_cast<std::int64_t>(
      (i128(speed) * point.n) % point.d);
  return 14 * residue >= point.d && 14 * residue <= 13 * point.d;
}

i128 prefix_numerator(const Rat& point, int speed) {
  const std::int64_t product = speed * point.n;
  const std::int64_t whole = product / point.d;
  const std::int64_t remainder = product % point.d;
  std::int64_t partial = 0;
  if (14 * remainder <= point.d) {
    partial = 0;
  } else if (14 * remainder >= 13 * point.d) {
    partial = 12 * point.d;
  } else {
    partial = 14 * remainder - point.d;
  }
  require(kCommon % point.d == 0, "common wall lattice");
  return i128(12) * whole * kCommon
       + i128(partial) * (kCommon / point.d);
}

std::string i128_text(i128 value) {
  if (value == 0) return "0";
  bool negative = value < 0;
  u128 magnitude = negative ? u128(-value) : u128(value);
  std::string answer;
  while (magnitude != 0) {
    answer.push_back(static_cast<char>('0' + magnitude % 10));
    magnitude /= 10;
  }
  if (negative) answer.push_back('-');
  std::reverse(answer.begin(), answer.end());
  return answer;
}

i128 gcd128(i128 first, i128 second) {
  if (first < 0) first = -first;
  if (second < 0) second = -second;
  while (second != 0) {
    const i128 remainder = first % second;
    first = second;
    second = remainder;
  }
  return first;
}

std::string fraction_text(i128 numerator, i128 denominator) {
  const i128 divisor = gcd128(numerator, denominator);
  return i128_text(numerator / divisor) + "/" + i128_text(denominator / divisor);
}

bool cover_search(const std::vector<std::uint32_t>& edges,
                  std::uint32_t chosen, int remaining,
                  std::uint32_t& witness) {
  std::uint32_t uncovered = 0;
  int greedy_matching = 0;
  std::uint32_t used = 0;
  for (const std::uint32_t edge : edges) {
    if ((edge & chosen) != 0) continue;
    if (uncovered == 0) uncovered = edge;
    if ((edge & used) == 0) {
      used |= edge;
      ++greedy_matching;
      if (greedy_matching > remaining) return false;
    }
  }
  if (uncovered == 0) {
    witness = chosen;
    return true;
  }
  if (remaining == 0) return false;
  while (uncovered != 0) {
    const std::uint32_t bit = uncovered & (0U - uncovered);
    if (cover_search(edges, chosen | bit, remaining - 1, witness)) return true;
    uncovered &= uncovered - 1;
  }
  return false;
}

std::uint32_t greedy_matching_mask(const std::vector<std::uint32_t>& edges,
                                   std::vector<std::uint32_t>& selected) {
  std::uint32_t used = 0;
  for (const std::uint32_t edge : edges) {
    if ((edge & used) != 0) continue;
    selected.push_back(edge);
    used |= edge;
  }
  return used;
}

std::string labels_text(std::uint32_t mask, const std::vector<int>& optional) {
  std::string answer = "(";
  bool first = true;
  for (int index = 0; index < 27; ++index) {
    if (((mask >> index) & 1U) == 0) continue;
    if (!first) answer += ',';
    first = false;
    answer += std::to_string(optional[index]);
  }
  answer += ')';
  return answer;
}

struct BankSummary {
  std::uint32_t mask;
  std::int64_t mass_numerator;
  int components;
  bool strict_limit;
  bool limit_equality;
  std::int64_t eventual_cutoff;
};

struct ExactRow {
  int q;
  int hyperedges;
  std::vector<std::uint32_t> hyperedge_masks;
  bool cover7;
  std::uint32_t cover;
  std::vector<std::uint32_t> matching;
  std::uint64_t equalities;
  i128 smallest_positive;
  i128 smallest_negative;
};

}  // namespace

int main() {
  const std::vector<int> anchors{120, 126, 143};
  const std::vector<int> pool{
      8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
      120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
      264, 286, 290};
  const std::set<int> anchor_set(anchors.begin(), anchors.end());
  std::vector<int> optional;
  for (const int value : pool) {
    if (!anchor_set.contains(value)) optional.push_back(value);
  }
  require(optional.size() == 27, "optional size");
  std::array<int, 291> optional_index{};
  optional_index.fill(-1);
  for (int index = 0; index < 27; ++index) optional_index[optional[index]] = index;

  std::set<Rat, RatLess> wall_set{Rat(0), Rat(1)};
  for (const int speed : pool) {
    for (int tooth = 0; tooth < speed; ++tooth) {
      wall_set.insert(Rat(14 * tooth + 1, 14 * speed));
      wall_set.insert(Rat(14 * tooth + 13, 14 * speed));
    }
  }
  std::vector<Rat> walls(wall_set.begin(), wall_set.end());
  require(walls.size() == 7134, "wall count");

  std::vector<std::uint32_t> failure_masks;
  std::array<std::uint64_t, 5> cell_histogram{};
  for (std::size_t cell = 0; cell + 1 < walls.size(); ++cell) {
    const Rat point = midpoint(walls[cell], walls[cell + 1]);
    std::uint32_t mask = 0;
    bool invalid = false;
    for (const int speed : pool) {
      if (safe_at(point, speed)) continue;
      if (anchor_set.contains(speed)) {
        invalid = true;
        break;
      }
      mask |= 1U << optional_index[speed];
      if (std::popcount(mask) > 3) {
        invalid = true;
        break;
      }
    }
    if (invalid) {
      failure_masks.push_back(1U << 31);
      ++cell_histogram[4];
    } else {
      failure_masks.push_back(mask);
      ++cell_histogram[std::popcount(mask)];
    }
  }

  std::vector<BankSummary> banks;
  for (int first = 0; first < 27; ++first) {
    for (int second = first + 1; second < 27; ++second) {
      for (int third = second + 1; third < 27; ++third) {
        const std::uint32_t triple = (1U << first) | (1U << second) | (1U << third);
        std::int64_t numerator = 0;
        int components = 0;
        bool previous = false;
        for (std::size_t cell = 0; cell < failure_masks.size(); ++cell) {
          const bool included = (failure_masks[cell] & ~triple) == 0;
          if (included) {
            const Rat left = walls[cell];
            const Rat right = walls[cell + 1];
            numerator += right.n * (kCommon / right.d)
                       - left.n * (kCommon / left.d);
          }
          if (included && !previous) ++components;
          previous = included;
        }
        const i128 limit_difference = i128(54) * numerator - i128(4) * kCommon;
        const bool strict_limit = limit_difference > 0;
        const bool limit_equality = limit_difference == 0;
        std::int64_t eventual_cutoff = -1;
        if (limit_difference < 0) {
        } else if (limit_difference > 0) {
          const i128 bound_numerator = i128(54) * components * kCommon;
          const i128 bound_denominator = i128(7) * limit_difference;
          eventual_cutoff = static_cast<std::int64_t>(
              (bound_numerator + bound_denominator - 1) / bound_denominator);
        }
        banks.push_back(BankSummary{triple, numerator, components, strict_limit,
                                    limit_equality, eventual_cutoff});
      }
    }
  }
  require(banks.size() == 2925, "triple bank count");

  auto exact_row = [&](int q) {
    std::vector<i128> cell_contributions;
    cell_contributions.reserve(failure_masks.size());
    i128 previous = prefix_numerator(walls[0], q);
    for (std::size_t cell = 0; cell < failure_masks.size(); ++cell) {
      const i128 current = prefix_numerator(walls[cell + 1], q);
      cell_contributions.push_back(current - previous);
      previous = current;
    }
    std::vector<std::uint32_t> hyperedges;
    std::uint64_t equalities = 0;
    i128 smallest_positive = 0;
    i128 smallest_negative = 0;
    for (const BankSummary& bank : banks) {
      i128 numerator = 0;
      for (std::size_t cell = 0; cell < failure_masks.size(); ++cell) {
        if ((failure_masks[cell] & ~bank.mask) == 0) {
          numerator += cell_contributions[cell];
        }
      }
      const i128 difference = i128(9) * numerator - i128(8) * q * kCommon;
      if (difference >= 0) {
        hyperedges.push_back(bank.mask);
        if (difference == 0) ++equalities;
        if (difference > 0
            && (smallest_positive == 0 || difference < smallest_positive)) {
          smallest_positive = difference;
        }
      } else if (smallest_negative == 0 || -difference < smallest_negative) {
        smallest_negative = -difference;
      }
    }
    std::uint32_t cover = 0;
    const bool cover7 = cover_search(hyperedges, 0, 7, cover);
    std::vector<std::uint32_t> matching;
    greedy_matching_mask(hyperedges, matching);
    return ExactRow{q, static_cast<int>(hyperedges.size()), hyperedges,
                    cover7, cover, matching, equalities,
                    smallest_positive, smallest_negative};
  };

  const std::vector<int> probes{9699, 9700};
  std::vector<ExactRow> exact_rows;
  for (const int q : probes) exact_rows.push_back(exact_row(q));

  std::vector<std::uint32_t> stable_edges;
  std::uint64_t limit_equalities = 0;
  for (const BankSummary& bank : banks) {
    if (bank.strict_limit) stable_edges.push_back(bank.mask);
    if (bank.limit_equality) ++limit_equalities;
  }
  std::uint32_t stable_cover = 0;
  const bool stable_cover7 = cover_search(stable_edges, 0, 7, stable_cover);
  const std::array<std::array<int, 3>, 8> matching_labels{{
      {{8, 84, 252}}, {{15, 88, 95}}, {{40, 85, 170}},
      {{63, 176, 193}}, {{10, 145, 290}}, {{42, 132, 264}},
      {{16, 168, 286}}, {{80, 190, 240}},
  }};
  std::vector<std::uint32_t> matching8;
  std::uint32_t matching_used = 0;
  std::int64_t matching_cutoff = -1;
  for (const auto& labels : matching_labels) {
    std::uint32_t edge = 0;
    for (const int label : labels) {
      require(optional_index[label] >= 0, "matching label is optional");
      edge |= 1U << optional_index[label];
    }
    require((edge & matching_used) == 0, "matching edges are disjoint");
    matching_used |= edge;
    const auto bank_it = std::find_if(
        banks.begin(), banks.end(),
        [edge](const BankSummary& bank) { return bank.mask == edge; });
    require(bank_it != banks.end() && bank_it->strict_limit,
            "matching edge is strict-limit");
    matching_cutoff = std::max(matching_cutoff, bank_it->eventual_cutoff);
    matching8.push_back(edge);
  }
  require(std::popcount(matching_used) == 24, "matching label count");
  require(matching_cutoff == 9700, "chosen matching cutoff");
  for (const ExactRow& row : exact_rows) {
    if (row.q != 9700) continue;
    for (const std::uint32_t edge : matching8) {
      require(std::find(row.hyperedge_masks.begin(), row.hyperedge_masks.end(),
                        edge) != row.hyperedge_masks.end(),
              "chosen matching edge is active at q=9700");
    }
  }
  int activated_at_9700 = 0;
  int activated_at_9699 = 0;
  std::int64_t previous_event = -1;
  for (const BankSummary& bank : banks) {
    if (!bank.strict_limit) continue;
    if (bank.eventual_cutoff <= 9700) ++activated_at_9700;
    if (bank.eventual_cutoff <= 9699) ++activated_at_9699;
    if (bank.eventual_cutoff < 9700) {
      previous_event = std::max(previous_event, bank.eventual_cutoff);
    }
  }
  require(activated_at_9700 == 666 && activated_at_9699 == 665,
          "activation counts");
  require(previous_event == 9687, "previous activation event");

  constexpr std::int64_t kEndpointBound = 14 * 290;
  constexpr std::int64_t kMidpointProductBound =
      2 * kEndpointBound * kEndpointBound;
  require(kMidpointProductBound < std::numeric_limits<std::int64_t>::max(),
          "midpoint products fit int64");
  for (const Rat& wall : walls) {
    require(wall.n >= 0 && wall.n <= kEndpointBound
                && wall.d > 0 && wall.d <= kEndpointBound,
            "reduced endpoint width");
  }

  std::cout << "LRC14_TRIPLE_DELETION_MATCHING_EVENTUAL_HAAR_THM4170_"
               "INDEPENDENT_CPP_20260826\n";
  std::cout << "walls=" << walls.size() << ";cells=" << failure_masks.size()
            << ";cell_hist=";
  for (const auto value : cell_histogram) std::cout << value << ',';
  std::cout << "\n";
  std::cout << "triple_banks=" << banks.size() << ";stable_edges="
            << stable_edges.size() << ";stable_cover7=" << stable_cover7
            << ";stable_cover=" << labels_text(stable_cover, optional)
            << ";limit_equalities=" << limit_equalities << "\n";
  std::cout << "chosen_matching_cutoff=" << matching_cutoff
            << ";previous_activation_event=" << previous_event
            << ";activated_at_9699=" << activated_at_9699
            << ";activated_at_9700=" << activated_at_9700 << "\n";
  std::cout << "eventual_matching8=(";
  for (std::size_t index = 0; index < matching8.size(); ++index) {
    if (index != 0) std::cout << ',';
    std::cout << labels_text(matching8[index], optional);
  }
  std::cout << ")\n";
  std::cout << "integer_width_endpoint_bound=" << kEndpointBound
            << ";midpoint_product_bound=" << kMidpointProductBound
            << ";int64_max=" << std::numeric_limits<std::int64_t>::max()
            << ";comparisons_use_int128=1\n";
  for (const std::uint32_t edge : matching8) {
    const auto bank_it = std::find_if(
        banks.begin(), banks.end(),
        [edge](const BankSummary& bank) { return bank.mask == edge; });
    require(bank_it != banks.end(), "matching bank lookup");
    const BankSummary& bank = *bank_it;
    const i128 limit_difference =
        i128(54) * bank.mass_numerator - i128(4) * kCommon;
    const i128 bound_numerator = i128(54) * bank.components * kCommon;
    const i128 bound_denominator = i128(7) * limit_difference;
    std::cout << "matching_edge=" << labels_text(edge, optional)
              << ";components=" << bank.components
              << ";bank_mass="
              << fraction_text(bank.mass_numerator, kCommon)
              << ";limit_surplus="
              << fraction_text(limit_difference, i128(63) * kCommon)
              << ";discrepancy_bound_B="
              << fraction_text(bound_numerator, bound_denominator)
              << ";eventual_cutoff=" << bank.eventual_cutoff << "\n";
  }
  for (const ExactRow& row : exact_rows) {
    std::cout << "q=" << row.q << ";hyperedges=" << row.hyperedges
              << ";cover7=" << row.cover7
              << ";cover=" << labels_text(row.cover, optional)
              << ";greedy_matching=" << row.matching.size()
              << ";matching_first8=(";
    for (std::size_t index = 0; index < std::min<std::size_t>(8, row.matching.size()); ++index) {
      if (index != 0) std::cout << ',';
      std::cout << labels_text(row.matching[index], optional);
    }
    std::cout << ");equalities=" << row.equalities;
    if (row.smallest_positive != 0) {
      std::cout << ";min_positive="
                << fraction_text(row.smallest_positive, i128(126) * row.q * kCommon);
    }
    if (row.smallest_negative != 0) {
      std::cout << ";min_negative="
                << fraction_text(row.smallest_negative, i128(126) * row.q * kCommon);
    }
    std::cout << "\n";
  }
  std::cout << "PASS\n";
}
