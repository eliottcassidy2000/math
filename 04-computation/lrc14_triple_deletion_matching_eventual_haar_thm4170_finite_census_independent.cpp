#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {
using i128 = __int128_t;
using u128 = __uint128_t;
constexpr std::int64_t kCommon = 18241159416480LL;

void require(bool value, const std::string& label) {
  if (!value) throw std::runtime_error("requirement failed: " + label);
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
  return Rat(left.n * right.d + right.n * left.d,
             2 * left.d * right.d);
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
  const bool negative = value < 0;
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

struct CellClass {
  int kind = -1;
  int first = -1;
  int second = -1;
  int third = -1;
};

struct Triple {
  int first;
  int second;
  int third;
  std::uint32_t mask;
};

int greedy_matching(const std::vector<std::uint32_t>& edges,
                    std::vector<std::uint32_t>* witness = nullptr) {
  std::uint32_t used = 0;
  int count = 0;
  for (const std::uint32_t edge : edges) {
    if ((edge & used) != 0) continue;
    used |= edge;
    ++count;
    if (witness != nullptr) witness->push_back(edge);
  }
  return count;
}

bool cover_search(const std::vector<std::uint32_t>& edges,
                  std::uint32_t chosen, int remaining,
                  std::uint32_t& witness) {
  std::uint32_t uncovered = 0;
  std::uint32_t matching_used = 0;
  int matching = 0;
  for (const std::uint32_t edge : edges) {
    if ((edge & chosen) != 0) continue;
    if (uncovered == 0) uncovered = edge;
    if ((edge & matching_used) == 0) {
      matching_used |= edge;
      if (++matching > remaining) return false;
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

struct Control {
  int q = 0;
  int edges = 0;
  int greedy = 0;
  bool cover7 = false;
  std::uint32_t cover = 0;
  std::uint64_t equalities = 0;
  std::vector<std::uint32_t> matching;
  std::vector<i128> matching_differences;
};
}  // namespace

int main() {
#ifdef _WIN32
  _setmode(_fileno(stdout), _O_BINARY);
#endif
  const std::vector<int> anchors{120, 126, 143};
  const std::vector<int> pool{
      8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
      120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
      264, 286, 290};
  const std::set<int> anchor_set(anchors.begin(), anchors.end());
  const std::set<int> pool_set(pool.begin(), pool.end());
  std::vector<int> optional;
  for (const int speed : pool) {
    if (!anchor_set.contains(speed)) optional.push_back(speed);
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
  const std::vector<Rat> walls(wall_set.begin(), wall_set.end());
  require(walls.size() == 7134, "wall count");

  std::array<std::array<int, 27>, 27> pair_index{};
  for (auto& row : pair_index) row.fill(-1);
  int pair_count = 0;
  for (int first = 0; first < 27; ++first) {
    for (int second = first + 1; second < 27; ++second) {
      pair_index[first][second] = pair_index[second][first] = pair_count++;
    }
  }
  std::array<std::array<std::array<int, 27>, 27>, 27> triple_index{};
  for (auto& plane : triple_index) {
    for (auto& row : plane) row.fill(-1);
  }
  std::vector<Triple> triples;
  for (int first = 0; first < 27; ++first) {
    for (int second = first + 1; second < 27; ++second) {
      for (int third = second + 1; third < 27; ++third) {
        const int index = static_cast<int>(triples.size());
        triple_index[first][second][third] = index;
        triples.push_back(Triple{first, second, third,
            (1U << first) | (1U << second) | (1U << third)});
      }
    }
  }
  require(pair_count == 351 && triples.size() == 2925, "bank index counts");

  std::vector<CellClass> classes;
  std::array<std::uint64_t, 5> cell_hist{};
  for (std::size_t cell = 0; cell + 1 < walls.size(); ++cell) {
    const Rat point = midpoint(walls[cell], walls[cell + 1]);
    std::vector<int> failures;
    bool invalid = false;
    for (const int speed : pool) {
      if (safe_at(point, speed)) continue;
      if (anchor_set.contains(speed)) {
        invalid = true;
        break;
      }
      failures.push_back(optional_index[speed]);
      if (failures.size() > 3) {
        invalid = true;
        break;
      }
    }
    CellClass row;
    if (!invalid) {
      row.kind = static_cast<int>(failures.size());
      if (row.kind >= 1) row.first = failures[0];
      if (row.kind >= 2) row.second = failures[1];
      if (row.kind >= 3) row.third = failures[2];
      ++cell_hist[row.kind];
    } else {
      ++cell_hist[4];
    }
    classes.push_back(row);
  }
  require(cell_hist == std::array<std::uint64_t, 5>{150,328,518,678,5459},
          "cell histogram");

  std::uint64_t equality_total = 0;
  int qualifiers = 0;
  int hostile = 0;
  int last_hostile = -1;
  int first_qualifier = -1;
  int greedy_qualifiers = 0;
  int branch_qualifiers = 0;
  std::map<int, int> hostile_tau_hist;
  std::vector<std::pair<int, std::uint32_t>> hostile_rows;
  std::vector<Control> controls;
  const std::set<int> control_set{1, 7, 200, 924, 925, 4959, 4960, 8266, 9699};
  std::uint64_t qualifier_hash = 14695981039346656037ULL;

  for (int q = 1; q <= 9699; ++q) {
    if (pool_set.contains(q)) continue;
    i128 base = 0;
    std::array<i128, 27> singles{};
    std::array<i128, 351> pairs{};
    std::array<i128, 2925> threes{};
    i128 previous = prefix_numerator(walls.front(), q);
    for (std::size_t cell = 0; cell < classes.size(); ++cell) {
      const i128 current = prefix_numerator(walls[cell + 1], q);
      const i128 contribution = current - previous;
      previous = current;
      require(contribution >= 0, "cell contribution");
      const CellClass& row = classes[cell];
      if (row.kind == 0) base += contribution;
      else if (row.kind == 1) singles[row.first] += contribution;
      else if (row.kind == 2) pairs[pair_index[row.first][row.second]] += contribution;
      else if (row.kind == 3) {
        threes[triple_index[row.first][row.second][row.third]] += contribution;
      }
    }

    std::vector<std::uint32_t> edges;
    std::vector<i128> edge_differences;
    std::uint64_t equalities = 0;
    for (std::size_t index = 0; index < triples.size(); ++index) {
      const Triple& edge = triples[index];
      const i128 numerator = base
          + singles[edge.first] + singles[edge.second] + singles[edge.third]
          + pairs[pair_index[edge.first][edge.second]]
          + pairs[pair_index[edge.first][edge.third]]
          + pairs[pair_index[edge.second][edge.third]] + threes[index];
      const i128 difference = i128(9) * numerator - i128(8) * q * kCommon;
      if (difference >= 0) {
        edges.push_back(edge.mask);
        edge_differences.push_back(difference);
      }
      if (difference == 0) ++equalities;
    }
    equality_total += equalities;
    std::vector<std::uint32_t> greedy_witness;
    const int greedy = greedy_matching(edges, &greedy_witness);
    std::uint32_t cover = 0;
    bool cover7 = false;
    if (greedy < 8) cover7 = cover_search(edges, 0, 7, cover);
    if (!cover7) {
      ++qualifiers;
      if (first_qualifier < 0) first_qualifier = q;
      if (greedy >= 8) ++greedy_qualifiers;
      else ++branch_qualifiers;
      qualifier_hash ^= static_cast<std::uint64_t>(q);
      qualifier_hash *= 1099511628211ULL;
    } else {
      ++hostile;
      last_hostile = q;
      int tau = -1;
      for (int size = 0; size <= 7; ++size) {
        std::uint32_t minimum_cover = 0;
        if (cover_search(edges, 0, size, minimum_cover)) {
          tau = size;
          cover = minimum_cover;
          break;
        }
      }
      require(tau >= 0, "hostile minimum cover");
      ++hostile_tau_hist[tau];
      hostile_rows.emplace_back(q, cover);
    }
    if (control_set.contains(q)) {
      std::vector<i128> matching_differences;
      for (const std::uint32_t edge : greedy_witness) {
        const auto position = std::find(edges.begin(), edges.end(), edge);
        require(position != edges.end(), "control matching edge lookup");
        matching_differences.push_back(
            edge_differences[static_cast<std::size_t>(position - edges.begin())]);
      }
      controls.push_back(Control{q, static_cast<int>(edges.size()), greedy,
                                 cover7, cover, equalities, greedy_witness,
                                 matching_differences});
    }
  }

  const std::vector<int> expected_hostile{
      3,6,22,24,25,46,48,50,55,64,70,72,75,83,93,96,100,103,105,
      110,122,127,128,140,147,153,158,166,172,173,183,186,192,206,
      210,220,256,260,270,282,294,306,313,320,325,332,346,366,372,
      384,416,440,462,512,516,520,550,567,744,768,924};
  std::vector<int> hostile_q;
  for (const auto& row : hostile_rows) hostile_q.push_back(row.first);
  require(hostile_q == expected_hostile, "hostile q ledger");
  require(hostile_tau_hist == std::map<int, int>{{0,4},{1,6},{2,5},{3,8},
          {4,5},{5,9},{6,10},{7,14}}, "hostile tau histogram");
  require(qualifiers == 9608 && hostile == 61 && last_hostile == 924,
          "qualifier census");
  require(greedy_qualifiers == 8837 && branch_qualifiers == 771,
          "qualifier proof split");
  require(equality_total == 0, "finite equality audit");
  require(qualifier_hash == 0x02784121a66537acULL,
          "qualifier wordwise XOR/multiply hash");
  const auto control_at = [&](int q) -> const Control& {
    const auto row = std::find_if(controls.begin(), controls.end(),
        [q](const Control& control) { return control.q == q; });
    require(row != controls.end(), "control lookup");
    return *row;
  };
  require(control_at(924).edges == 258 && control_at(924).cover7
              && labels_text(control_at(924).cover, optional)
                 == "(16,85,88,145,168,252)",
          "q924 hostile control");
  require(control_at(925).edges == 1705 && !control_at(925).cover7
              && control_at(925).greedy == 9,
          "q925 qualifier control");

  std::cout << "LRC14_TRIPLE_DELETION_THM4170_FINITE_CENSUS_"
               "INDEPENDENT_CPP_20260826\n";
  std::cout << "universe=1<=q<=9699,q_notin_P;count="
            << 9699 - static_cast<int>(pool.size()) << "\n";
  std::cout << "walls=" << walls.size() << ";cells=" << classes.size()
            << ";cell_hist=";
  for (const auto value : cell_hist) std::cout << value << ',';
  std::cout << "\n";
  std::cout << "qualifiers=" << qualifiers << ";hostile=" << hostile
            << ";first_qualifier=" << first_qualifier
            << ";last_hostile=" << last_hostile
            << ";greedy_matching8_qualifiers=" << greedy_qualifiers
            << ";branch_only_qualifiers=" << branch_qualifiers
            << ";equality_total=" << equality_total << "\n";
  std::cout << "hostile_tau_hist=";
  for (const auto& [tau, count] : hostile_tau_hist) {
    std::cout << tau << ':' << count << ',';
  }
  std::cout << "\n";
  std::cout << "hostile_rows=";
  for (const auto& [q, cover] : hostile_rows) {
    std::cout << q << ':' << labels_text(cover, optional) << ';';
  }
  std::cout << "\n";
  for (const Control& row : controls) {
    std::cout << "control_q=" << row.q << ";edges=" << row.edges
              << ";greedy=" << row.greedy << ";cover7=" << row.cover7
              << ";cover=" << labels_text(row.cover, optional)
              << ";equalities=" << row.equalities
              << ";matching_first8=(";
    for (std::size_t index = 0;
         index < std::min<std::size_t>(8, row.matching.size()); ++index) {
      if (index != 0) std::cout << ',';
      std::cout << labels_text(row.matching[index], optional) << ':'
                << fraction_text(row.matching_differences[index],
                                 i128(126) * row.q * kCommon);
    }
    std::cout << ")\n";
  }
  std::cout << "qualifier_q_word_xor_mul64=" << std::hex << std::setfill('0')
            << std::setw(16) << qualifier_hash << std::dec << "\n";
  std::cout << "certificate=all_q>=925_by_finite_925_to_9699_plus_"
               "eventual_q>=9700\n";
  std::cout << "PASS\n";
}
