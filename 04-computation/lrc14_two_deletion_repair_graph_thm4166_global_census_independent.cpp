#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
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

constexpr int kLimit = 49493;
constexpr std::int64_t kCommon = 18241159416480LL;

void require(bool condition, const std::string& label) {
  if (!condition) throw std::runtime_error("requirement failed: " + label);
}

struct Rat {
  std::int64_t n;
  std::int64_t d;
  Rat(std::int64_t numerator = 0, std::int64_t denominator = 1) {
    require(denominator > 0, "positive denominator");
    const std::int64_t divisor = std::gcd(
        numerator < 0 ? -numerator : numerator, denominator);
    n = numerator / divisor;
    d = denominator / divisor;
  }
};

struct RatLess {
  bool operator()(const Rat& left, const Rat& right) const {
    return i128(left.n) * right.d < i128(right.n) * left.d;
  }
};

bool safe_at(const Rat& point, int speed) {
  const std::int64_t residue = static_cast<std::int64_t>(
      (i128(speed) * point.n) % point.d);
  return 14 * residue >= point.d && 14 * residue <= 13 * point.d;
}

Rat midpoint(const Rat& left, const Rat& right) {
  return Rat(
      left.n * right.d + right.n * left.d,
      2 * left.d * right.d);
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
  require(kCommon % point.d == 0, "wall denominator divides common lattice");
  return i128(12) * whole * kCommon
       + i128(partial) * (kCommon / point.d);
}

std::string i128_text(i128 value) {
  if (value == 0) return "0";
  bool negative = value < 0;
  u128 magnitude = negative ? u128(-value) : u128(value);
  std::string result;
  while (magnitude != 0) {
    result.push_back(static_cast<char>('0' + magnitude % 10));
    magnitude /= 10;
  }
  if (negative) result.push_back('-');
  std::reverse(result.begin(), result.end());
  return result;
}

i128 gcd128(i128 left, i128 right) {
  if (left < 0) left = -left;
  if (right < 0) right = -right;
  while (right != 0) {
    const i128 remainder = left % right;
    left = right;
    right = remainder;
  }
  return left;
}

std::string fraction_text(i128 numerator, i128 denominator) {
  const i128 divisor = gcd128(numerator, denominator);
  numerator /= divisor;
  denominator /= divisor;
  return i128_text(numerator) + "/" + i128_text(denominator);
}

struct CellClass {
  // kind -1 ignored, 0 base, 1 single, 2 pair.
  int kind = -1;
  int first = -1;
  int second = -1;
};

struct GraphRow {
  int q;
  int edges;
  int alpha;
  int tau;
  std::uint32_t witness;
  std::array<std::uint32_t, 27> adjacency;
};

std::pair<int, std::uint32_t> maximum_independent_set(
    const std::array<std::uint32_t, 27>& adjacency) {
  constexpr std::uint32_t universe = (1U << 27) - 1;
  std::array<std::uint32_t, 27> complement{};
  for (int vertex = 0; vertex < 27; ++vertex) {
    complement[vertex] = universe & ~(1U << vertex) & ~adjacency[vertex];
  }
  int best_size = 0;
  std::uint32_t best_mask = 0;
  auto search = [&](auto&& self, std::uint32_t chosen, std::uint32_t candidates,
                    std::uint32_t excluded) -> void {
    const int chosen_size = std::popcount(chosen);
    if (chosen_size + std::popcount(candidates) <= best_size) return;
    if (candidates == 0 && excluded == 0) {
      if (chosen_size > best_size) {
        best_size = chosen_size;
        best_mask = chosen;
      }
      return;
    }
    const std::uint32_t pivot_pool = candidates | excluded;
    std::uint32_t extensions = candidates;
    if (pivot_pool != 0) {
      int pivot = std::countr_zero(pivot_pool);
      int pivot_score = -1;
      std::uint32_t remaining = pivot_pool;
      while (remaining != 0) {
        const int vertex = std::countr_zero(remaining);
        const int score = std::popcount(candidates & complement[vertex]);
        if (score > pivot_score) {
          pivot = vertex;
          pivot_score = score;
        }
        remaining &= remaining - 1;
      }
      extensions = candidates & ~complement[pivot];
    }
    while (extensions != 0) {
      const std::uint32_t bit = extensions & (0U - extensions);
      const int vertex = std::countr_zero(bit);
      self(self, chosen | bit, candidates & complement[vertex],
           excluded & complement[vertex]);
      candidates &= ~bit;
      excluded |= bit;
      extensions &= ~bit;
      if (chosen_size + std::popcount(candidates) <= best_size) return;
    }
  };
  search(search, 0, universe, 0);
  return {best_size, best_mask};
}

bool has_vertex_cover_at_most(
    const std::array<std::uint32_t, 27>& adjacency,
    std::uint32_t remaining,
    int budget) {
  for (int first = 0; first < 27; ++first) {
    if (((remaining >> first) & 1U) == 0U) continue;
    const std::uint32_t neighbours = adjacency[first] & remaining;
    if (neighbours == 0U) continue;
    if (budget == 0) return false;
    const int second = std::countr_zero(neighbours);
    return has_vertex_cover_at_most(
               adjacency, remaining & ~(1U << first), budget - 1)
        || has_vertex_cover_at_most(
               adjacency, remaining & ~(1U << second), budget - 1);
  }
  return true;
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

struct MarginRecord {
  bool initialized = false;
  i128 difference = 0;
  int q = 1;
  int first = 0;
  int second = 0;
};

void update_margin(MarginRecord& record, i128 difference, int q,
                   int first, int second) {
  require(difference > 0, "positive margin magnitude");
  if (!record.initialized
      || difference * record.q < record.difference * q
      || (difference * record.q == record.difference * q
          && std::tuple(q, first, second)
             < std::tuple(record.q, record.first, record.second))) {
    record = MarginRecord{true, difference, q, first, second};
  }
}

}  // namespace

int main() {
  const std::vector<int> anchors{120, 126, 143};
  const std::vector<int> pool{
      8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
      120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
      264, 286, 290};
  const std::set<int> anchor_set(anchors.begin(), anchors.end());
  const std::set<int> pool_set(pool.begin(), pool.end());
  std::vector<int> optional;
  for (const int value : pool) {
    if (!anchor_set.contains(value)) optional.push_back(value);
  }
  require(optional.size() == 27, "optional size");
  std::array<int, 291> optional_index{};
  optional_index.fill(-1);
  for (int index = 0; index < 27; ++index) optional_index[optional[index]] = index;

  std::set<Rat, RatLess> wall_set;
  wall_set.insert(Rat(0));
  wall_set.insert(Rat(1));
  for (const int speed : pool) {
    for (int tooth = 0; tooth < speed; ++tooth) {
      wall_set.insert(Rat(14 * tooth + 1, 14 * speed));
      wall_set.insert(Rat(14 * tooth + 13, 14 * speed));
    }
  }
  std::vector<Rat> walls(wall_set.begin(), wall_set.end());
  require(walls.size() == 7134, "global wall count");
  std::vector<CellClass> classes;
  std::array<std::array<int, 27>, 27> pair_index{};
  for (auto& row : pair_index) row.fill(-1);
  int next_pair = 0;
  for (int first = 0; first < 27; ++first) {
    for (int second = first + 1; second < 27; ++second) {
      pair_index[first][second] = pair_index[second][first] = next_pair++;
    }
  }
  require(next_pair == 351, "pair index count");

  std::array<std::uint64_t, 4> class_histogram{};
  for (std::size_t cell = 0; cell + 1 < walls.size(); ++cell) {
    const Rat point = midpoint(walls[cell], walls[cell + 1]);
    std::vector<int> failures;
    bool anchor_failure = false;
    for (const int speed : pool) {
      if (safe_at(point, speed)) continue;
      if (anchor_set.contains(speed)) {
        anchor_failure = true;
        break;
      }
      failures.push_back(optional_index[speed]);
      if (failures.size() > 2) break;
    }
    CellClass category;
    if (!anchor_failure && failures.size() <= 2) {
      category.kind = static_cast<int>(failures.size());
      if (!failures.empty()) category.first = failures[0];
      if (failures.size() == 2) category.second = failures[1];
      ++class_histogram[category.kind];
    } else {
      ++class_histogram[3];
    }
    classes.push_back(category);
  }
  require(classes.size() == 7133, "cell class count");

  std::map<int, int> tau_histogram;
  std::map<int, int> tau_histogram_200;
  std::vector<GraphRow> admitted;
  std::vector<int> admitted_200;
  std::uint64_t equality_edges = 0;
  MarginRecord positive_margin;
  MarginRecord negative_margin;
  std::uint64_t semantic_hash = 1469598103934665603ULL;
  auto fnv = [&](std::uint64_t value) {
    semantic_hash ^= value;
    semantic_hash *= 1099511628211ULL;
  };
  GraphRow q27{};
  GraphRow q235{};
  GraphRow q8265{};
  GraphRow q8266{};
  GraphRow q49493{};
  GraphRow q49494{};

  auto evaluate = [&](int q, bool track_margins) {
    i128 base = 0;
    std::array<i128, 27> singles{};
    std::array<i128, 351> pairs{};
    i128 previous = prefix_numerator(walls[0], q);
    for (std::size_t cell = 0; cell < classes.size(); ++cell) {
      const i128 current = prefix_numerator(walls[cell + 1], q);
      const i128 contribution = current - previous;
      require(contribution >= 0, "nonnegative safe cell contribution");
      previous = current;
      const CellClass category = classes[cell];
      if (category.kind == 0) {
        base += contribution;
      } else if (category.kind == 1) {
        singles[category.first] += contribution;
      } else if (category.kind == 2) {
        pairs[pair_index[category.first][category.second]] += contribution;
      }
    }
    std::array<std::uint32_t, 27> adjacency{};
    int edge_count = 0;
    for (int first = 0; first < 27; ++first) {
      for (int second = first + 1; second < 27; ++second) {
        const i128 numerator = base + singles[first] + singles[second]
                             + pairs[pair_index[first][second]];
        // Measure denominator is 14*q*kCommon; threshold difference numerator
        // over 126*q*kCommon is 9*numerator - 8*q*kCommon.
        const i128 difference = i128(9) * numerator - i128(8) * q * kCommon;
        if (difference >= 0) {
          adjacency[first] |= 1U << second;
          adjacency[second] |= 1U << first;
          ++edge_count;
          if (difference == 0 && track_margins) {
            ++equality_edges;
          } else if (difference > 0 && track_margins) {
            update_margin(positive_margin, difference, q,
                          optional[first], optional[second]);
          }
        } else if (track_margins) {
          update_margin(negative_margin, -difference, q,
                        optional[first], optional[second]);
        }
      }
    }
    const auto [alpha, witness] = maximum_independent_set(adjacency);
    return GraphRow{q, edge_count, alpha, 27 - alpha, witness, adjacency};
  };

  std::uint64_t universe_count = 0;
  std::uint64_t bounded_cover_checks = 0;
  for (int q = 1; q <= kLimit; ++q) {
    if (pool_set.contains(q)) continue;
    ++universe_count;
    const GraphRow row = evaluate(q, true);
    const bool cover_at_most_seven =
        has_vertex_cover_at_most(row.adjacency, (1U << 27) - 1U, 7);
    require(cover_at_most_seven == (row.tau <= 7),
            "bounded vertex-cover recursion disagrees with alpha/tau");
    ++bounded_cover_checks;
    ++tau_histogram[row.tau];
    if (q <= 200) ++tau_histogram_200[row.tau];
    fnv(q);
    fnv(row.edges);
    fnv(row.alpha);
    for (const std::uint32_t mask : row.adjacency) fnv(mask);
    if (row.tau > 7) {
      admitted.push_back(row);
      if (q <= 200) admitted_200.push_back(q);
    }
    if (q == 27) q27 = row;
    if (q == 235) q235 = row;
    if (q == 8265) q8265 = row;
    if (q == 49493) q49493 = row;
  }
  q8266 = evaluate(8266, false);
  q49494 = evaluate(49494, false);

  require(universe_count == 49463, "global q universe count");
  require(bounded_cover_checks == universe_count,
          "bounded vertex-cover check count");
  const std::map<int, int> expected_tau_histogram{
      {0,45},{1,127},{2,124},{3,596},{4,793},{5,6241},{6,38003},
      {7,2502},{8,377},{9,435},{10,81},{11,45},{12,32},{13,19},
      {14,11},{15,13},{16,4},{17,5},{18,3},{19,5},{20,2}};
  require(tau_histogram == expected_tau_histogram, "global tau histogram");
  const std::map<int, int> expected_tau_histogram_200{
      {0,25},{1,25},{2,25},{3,7},{4,10},{5,11},{6,11},{7,8},{8,6},
      {9,12},{10,15},{11,3},{12,8},{13,3},{15,2},{16,1},{18,2},{19,1}};
  require(tau_histogram_200 == expected_tau_histogram_200,
          "q<=200 Fraction histogram control");
  const std::vector<int> expected_admitted_200{
      4,5,9,13,18,21,27,32,33,34,39,41,43,49,54,56,66,68,77,78,
      81,82,86,91,97,98,101,107,115,117,118,119,121,131,133,138,
      139,142,149,152,156,164,167,171,174,178,181,182,187,191,194,
      196,199};
  require(admitted_200 == expected_admitted_200, "q<=200 Fraction census control");
  require(admitted.size() == 1032 && admitted.back().q == 8265,
          "global admitted count and last label");
  require((q27.edges == 217 && q27.tau == 16), "q27 control");
  require((q235.edges == 238 && q235.tau == 17), "q235 control");
  require((q8265.edges == 52 && q8265.tau == 8), "last admitted control");
  require((q8266.edges == 45 && q8266.tau == 6), "immediate hostile control");
  require(q49494.tau <= 7, "post-cutoff hostile");
  const std::vector<int> expected_one_deletion{
      5, 66, 182, 298, 336, 340, 380, 386, 528, 572};
  for (const int q : expected_one_deletion) {
    require(std::any_of(admitted.begin(), admitted.end(),
                        [q](const GraphRow& row) { return row.q == q; }),
            "one-deletion family is subsumed");
  }
  require(semantic_hash == 0xcbc7947c3ed0e41eULL, "semantic FNV64");

  std::cout << "THM4166_TWO_DELETION_GLOBAL_CPP_CENSUS_20260826\n";
  std::cout << "q_universe=1..49493 outside P;count=" << universe_count << "\n";
  std::cout << "walls=" << walls.size() << ";cells=" << classes.size()
            << ";common_denominator=" << kCommon << ";cell_class_hist=";
  for (const std::uint64_t value : class_histogram) std::cout << value << ',';
  std::cout << "\n";
  std::cout << "threshold_equalities=" << equality_edges << "\n";
  std::cout << "bounded_cover7_checks=" << bounded_cover_checks << "\n";
  std::cout << "tau_histogram=";
  for (const auto& [tau, count] : tau_histogram) {
    std::cout << tau << ':' << count << ',';
  }
  std::cout << "\n";
  std::cout << "q_le_200_tau_histogram=";
  for (const auto& [tau, count] : tau_histogram_200) {
    std::cout << tau << ':' << count << ',';
  }
  std::cout << ";q_le_200_admitted_count=" << admitted_200.size() << "\n";
  std::cout << "admitted_count=" << admitted.size() << "\n";
  std::cout << "admitted_q=(";
  for (std::size_t index = 0; index < admitted.size(); ++index) {
    if (index != 0) std::cout << ',';
    std::cout << admitted[index].q;
  }
  std::cout << ")\n";
  std::cout << "admitted_rows=(";
  for (std::size_t index = 0; index < admitted.size(); ++index) {
    if (index != 0) std::cout << ',';
    const GraphRow& row = admitted[index];
    std::cout << '(' << row.q << ',' << row.edges << ',' << row.alpha
              << ',' << row.tau << ')';
  }
  std::cout << ")\n";
  const i128 margin_denominator_factor = i128(126) * kCommon;
  std::cout << "smallest_positive_margin="
            << fraction_text(positive_margin.difference,
                             margin_denominator_factor * positive_margin.q)
            << ";q_r_s=" << positive_margin.q << ',' << positive_margin.first
            << ',' << positive_margin.second << "\n";
  std::cout << "smallest_negative_deficit="
            << fraction_text(negative_margin.difference,
                             margin_denominator_factor * negative_margin.q)
            << ";q_r_s=" << negative_margin.q << ',' << negative_margin.first
            << ',' << negative_margin.second << "\n";
  const auto row_text = [&](const GraphRow& row) {
    std::cout << "q=" << row.q << ";edges=" << row.edges
              << ";alpha=" << row.alpha << ";tau=" << row.tau
              << ";independent_witness=" << labels_text(row.witness, optional)
              << "\n";
  };
  row_text(q27);
  row_text(q235);
  row_text(q8265);
  row_text(q8266);
  row_text(q49493);
  row_text(q49494);
  std::cout << "one_deletion_ten_subsumed=True\n";
  constexpr std::uint64_t bodies_per_q = 888030;
  const std::uint64_t extension_bodies = admitted.size() * bodies_per_q;
  std::cout << "bodies_per_q=" << bodies_per_q
            << ";extension_bodies=" << extension_bodies
            << ";with_old_thm4156=" << extension_bodies + 2220075 << "\n";
  std::cout << "semantic_fnv64=" << std::hex << semantic_hash << std::dec << "\n";
  std::cout << "PASS\n";
}
