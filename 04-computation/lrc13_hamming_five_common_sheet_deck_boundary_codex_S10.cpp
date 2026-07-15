#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Deterministic exact replay for THM-823.
//
// Every scalar comparison is made after multiplication by SCALE=360360,
// which is divisible by 13 and by every finite order in {2,...,12}.  The
// formal order INF has constant capacity 2/13.  Sheet decisions use literal
// bit sets on the common lcm, and maximin decisions use rational branch
// intersections only.

namespace {

constexpr int P = 13;
constexpr int INF = 13;  // display/order code; finite scan orders are 2..12
constexpr int SCALE = 360360;
constexpr std::array<int, 4> H = {1, 5, 8, 12};

using Labels = std::array<int, 5>;
using Orders = std::array<int, 5>;
using Pattern = std::array<int, 5>;

struct Row {
  Labels labels{};
  Orders orders{};
};

int inv13(int x) {
  for (int y = 1; y < P; ++y) {
    if (x * y % P == 1) return y;
  }
  assert(false);
  return 0;
}

int signed_mod(long long x, int modulus) {
  int r = static_cast<int>(x % modulus);
  if (r < 0) r += modulus;
  return r <= modulus / 2 ? r : r - modulus;
}

int base_count(int d, int replacement, int owner) {
  assert(1 <= d && d <= 12);
  const int target = d * replacement * inv13(owner) % P;
  int answer = 0;
  for (int z = -d + 1; z <= d; ++z) {
    int residue = (z - target) % P;
    if (residue < 0) residue += P;
    answer += residue == 0;
  }
  return answer;
}

int count_capacity(int D, int replacement, int owner) {
  assert(D > 0 && D % P != 0);
  const int q = D / P;
  const int d = D % P;
  return 2 * q + base_count(d, replacement, owner);
}

long long gcdll(long long a, long long b) {
  while (b) {
    const long long r = a % b;
    a = b;
    b = r;
  }
  return a < 0 ? -a : a;
}

struct Fraction {
  long long n = 0;
  long long d = 1;

  Fraction() = default;
  Fraction(long long numerator, long long denominator) : n(numerator), d(denominator) {
    assert(d != 0);
    if (d < 0) {
      n = -n;
      d = -d;
    }
    const long long g = gcdll(n, d);
    n /= g;
    d /= g;
  }
};

bool operator<(const Fraction& a, const Fraction& b) {
  return static_cast<__int128>(a.n) * b.d < static_cast<__int128>(b.n) * a.d;
}
bool operator==(const Fraction& a, const Fraction& b) {
  return a.n == b.n && a.d == b.d;
}
bool operator>(const Fraction& a, const Fraction& b) { return b < a; }

std::ostream& operator<<(std::ostream& out, const Fraction& x) {
  if (x.d == 1) return out << x.n;
  return out << x.n << '/' << x.d;
}

std::string order_string(int D) {
  return D == INF ? "inf" : std::to_string(D);
}

std::string pattern_string(const Pattern& pattern) {
  std::ostringstream out;
  out << '(';
  for (int i = 0; i < 5; ++i) {
    if (i) out << ',';
    out << order_string(pattern[i]);
  }
  out << ')';
  return out.str();
}

std::string decorated_string(const Row& row) {
  std::ostringstream out;
  out << '(';
  for (int i = 0; i < 5; ++i) {
    if (i) out << ',';
    out << row.labels[i] << ':' << order_string(row.orders[i]);
  }
  out << ')';
  return out.str();
}

Pattern order_pattern(const Orders& orders) {
  Pattern answer = orders;
  std::sort(answer.begin(), answer.end());
  return answer;
}

using Canonical = std::array<unsigned char, 12>;

Canonical canonical_decorated(const Row& row) {
  Canonical best{};
  best.fill(255);
  for (int multiplier = 1; multiplier < P; ++multiplier) {
    Canonical candidate{};
    candidate.fill(0);
    for (int i = 0; i < 5; ++i) {
      const int label = multiplier * row.labels[i] % P;
      candidate[label - 1] = static_cast<unsigned char>(row.orders[i]);
    }
    if (candidate < best) best = candidate;
  }
  return best;
}

uint64_t fnv_append(uint64_t hash, uint64_t value) {
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (8 * byte)) & 0xffU;
    hash *= 1099511628211ULL;
  }
  return hash;
}

std::string hex64(uint64_t value) {
  std::ostringstream out;
  out << std::hex << std::setfill('0') << std::setw(16) << value;
  return out.str();
}

std::array<std::array<std::array<int, 13>, 13>, 14> scaled_capacity{};

void initialize_capacities() {
  for (int D = 2; D <= 12; ++D) {
    for (int r = 1; r < P; ++r) {
      for (int o = 1; o < P; ++o) {
        scaled_capacity[D][r][o] = count_capacity(D, r, o) * (SCALE / D);
      }
    }
  }
  for (int r = 1; r < P; ++r) {
    for (int o = 1; o < P; ++o) {
      scaled_capacity[INF][r][o] = 2 * (SCALE / P);
    }
  }
}

struct ScalarCensus {
  long long total = 0;
  std::array<long long, 6> by_infinity{};
  std::set<Pattern> patterns;
  std::vector<Row> finite_rows;
  std::map<Pattern, long long> infinity_pattern_rows;
  std::map<Pattern, std::pair<int, int>> infinity_margin_range;
  std::map<Pattern, std::set<Canonical>> infinity_orbits;
  std::map<Pattern, Row> infinity_first;
  uint64_t row_hash = 1469598103934665603ULL;
};

void scan_orders_for_labels(const Labels& labels, ScalarCensus& census) {
  constexpr std::array<int, 12> alphabet = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, INF};
  std::array<std::array<int, 5>, 6> optimistic{};
  for (int owner = 0; owner < 5; ++owner) optimistic[5][owner] = 0;
  for (int index = 4; index >= 0; --index) {
    for (int owner = 0; owner < 5; ++owner) {
      int best = 0;
      for (int D : alphabet) {
        best = std::max(best, scaled_capacity[D][labels[index]][labels[owner]]);
      }
      optimistic[index][owner] = optimistic[index + 1][owner] + best;
    }
  }

  Orders orders{};
  std::array<int, 5> sums{};
  auto dfs = [&](auto&& self, int index) -> void {
    if (index == 5) {
      for (int owner = 0; owner < 5; ++owner) {
        if (sums[owner] < SCALE) return;
      }
      Row row{labels, orders};
      ++census.total;
      int infinity_count = 0;
      for (int D : orders) infinity_count += D == INF;
      ++census.by_infinity[infinity_count];
      const Pattern pattern = order_pattern(orders);
      census.patterns.insert(pattern);
      for (int x : labels) census.row_hash = fnv_append(census.row_hash, x);
      for (int x : orders) census.row_hash = fnv_append(census.row_hash, x);
      for (int x : sums) census.row_hash = fnv_append(census.row_hash, x);
      if (infinity_count == 0) {
        census.finite_rows.push_back(row);
      } else {
        ++census.infinity_pattern_rows[pattern];
        const int margin = *std::min_element(sums.begin(), sums.end()) - SCALE;
        auto it = census.infinity_margin_range.find(pattern);
        if (it == census.infinity_margin_range.end()) {
          census.infinity_margin_range[pattern] = {margin, margin};
          census.infinity_first[pattern] = row;
        } else {
          it->second.first = std::min(it->second.first, margin);
          it->second.second = std::max(it->second.second, margin);
          if (decorated_string(row) < decorated_string(census.infinity_first[pattern])) {
            census.infinity_first[pattern] = row;
          }
        }
        census.infinity_orbits[pattern].insert(canonical_decorated(row));
      }
      return;
    }

    for (int D : alphabet) {
      orders[index] = D;
      bool possible = true;
      for (int owner = 0; owner < 5; ++owner) {
        sums[owner] += scaled_capacity[D][labels[index]][labels[owner]];
        if (sums[owner] + optimistic[index + 1][owner] < SCALE) possible = false;
      }
      if (possible) self(self, index + 1);
      for (int owner = 0; owner < 5; ++owner) {
        sums[owner] -= scaled_capacity[D][labels[index]][labels[owner]];
      }
    }
  };
  dfs(dfs, 0);
}

ScalarCensus scalar_census() {
  ScalarCensus census;
  for (int a = 2; a <= 9; ++a) {
    for (int b = a + 1; b <= 10; ++b) {
      for (int c = b + 1; c <= 11; ++c) {
        for (int d = c + 1; d <= 12; ++d) {
          scan_orders_for_labels(Labels{1, a, b, c, d}, census);
        }
      }
    }
  }
  return census;
}

int crt_u(int D, int label, int e) {
  assert(std::gcd(D, e) == 1);
  for (int u = 1; u < P * D; ++u) {
    if (u % D == e % D && u % P == D * label % P) return u;
  }
  assert(false);
  return 0;
}

std::vector<int> units(int D) {
  std::vector<int> answer;
  for (int e = 1; e <= D; ++e) {
    if (std::gcd(D, e) == 1) answer.push_back(e);
  }
  return answer;
}

bool eligible_base_sheet(int D, int label, int e, int owner, int ell) {
  const int u = crt_u(D, label, e);
  const int residue = signed_mod(1LL * u * (inv13(owner) + P * ell), P * D);
  return -D < residue && residue <= D;
}

using Bits = std::vector<uint64_t>;

struct SheetResult {
  long long compatible_choices = 0;
  std::vector<Orders> unit_words;
  int best_min_covered = -1;
  long long best_choice_count = 0;
  std::map<std::array<int, 5>, long long> best_profiles;
};

SheetResult sheet_choices(const Row& row, bool collect_all, bool collect_best) {
  int common = 1;
  for (int D : row.orders) common = std::lcm(common, D);
  const int words = (common + 63) / 64;
  std::array<std::vector<int>, 5> choices;
  for (int i = 0; i < 5; ++i) choices[i] = units(row.orders[i]);

  // bits[colour][choice][owner][word]
  std::array<std::vector<std::array<Bits, 5>>, 5> bits;
  for (int colour = 0; colour < 5; ++colour) {
    const int D = row.orders[colour];
    bits[colour].resize(choices[colour].size());
    for (int choice = 0; choice < static_cast<int>(choices[colour].size()); ++choice) {
      const int e = choices[colour][choice];
      for (int owner = 0; owner < 5; ++owner) {
        bits[colour][choice][owner].assign(words, 0);
        for (int ell0 = 0; ell0 < D; ++ell0) {
          if (!eligible_base_sheet(D, row.labels[colour], e, row.labels[owner], ell0)) continue;
          for (int ell = ell0; ell < common; ell += D) {
            bits[colour][choice][owner][ell / 64] |= 1ULL << (ell % 64);
          }
        }
      }
    }
  }

  SheetResult result;
  std::array<int, 5> selected{};
  std::array<int, 5> selected_units{};
  auto dfs = [&](auto&& self, int colour) -> bool {
    if (colour < 5) {
      bool found = false;
      for (int choice = 0; choice < static_cast<int>(choices[colour].size()); ++choice) {
        selected[colour] = choice;
        selected_units[colour] = choices[colour][choice];
        found = self(self, colour + 1) || found;
        if (found && !collect_all && !collect_best) return true;
      }
      return found;
    }

    bool covers_all = true;
    std::array<int, 5> profile{};
    for (int owner = 0; owner < 5; ++owner) {
      int covered = 0;
      for (int word = 0; word < words; ++word) {
        uint64_t union_word = 0;
        for (int i = 0; i < 5; ++i) union_word |= bits[i][selected[i]][owner][word];
        if (word + 1 == words && common % 64) {
          union_word &= (1ULL << (common % 64)) - 1;
        }
        covered += __builtin_popcountll(union_word);
      }
      profile[owner] = covered;
      covers_all = covers_all && covered == common;
    }
    const int minimum = *std::min_element(profile.begin(), profile.end());
    if (collect_best) {
      if (minimum > result.best_min_covered) {
        result.best_min_covered = minimum;
        result.best_choice_count = 1;
        result.best_profiles.clear();
        result.best_profiles[profile] = 1;
      } else if (minimum == result.best_min_covered) {
        ++result.best_choice_count;
        ++result.best_profiles[profile];
      }
    }
    if (covers_all) {
      ++result.compatible_choices;
      result.unit_words.push_back(selected_units);
    }
    return covers_all;
  };
  dfs(dfs, 0);
  return result;
}

struct ExactMaximum {
  Fraction maximum{};
  Fraction first_witness{};
};

ExactMaximum exact_maximum(const std::vector<int>& speeds) {
  ExactMaximum answer;
  bool initialized = false;

  auto test = [&](long long k, long long denominator) {
    if (denominator <= 0 || k < 0 || k > denominator) return;
    long long minimum = denominator;
    for (int speed : speeds) {
      long long residue = (static_cast<__int128>(speed) * k) % denominator;
      minimum = std::min(minimum, std::min(residue, denominator - residue));
    }
    const Fraction value(minimum, denominator);
    const Fraction witness(k, denominator);
    if (!initialized || value > answer.maximum ||
        (value == answer.maximum && witness < answer.first_witness)) {
      answer.maximum = value;
      answer.first_witness = witness;
      initialized = true;
    }
  };

  for (int speed : speeds) {
    const int denominator = 2 * speed;
    for (int k = 0; k <= denominator; ++k) test(k, denominator);
  }
  for (int i = 0; i < static_cast<int>(speeds.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(speeds.size()); ++j) {
      for (int denominator : {speeds[i] + speeds[j], std::abs(speeds[i] - speeds[j])}) {
        if (denominator == 0) continue;
        for (int k = 0; k <= denominator; ++k) test(k, denominator);
      }
    }
  }
  return answer;
}

std::vector<int> packet_speeds(const Row& row, const Orders& unit_word) {
  assert(std::all_of(row.orders.begin(), row.orders.end(), [](int D) { return D == 3; }));
  std::set<int> speeds;
  for (int label = 1; label < P; ++label) {
    if (std::find(row.labels.begin(), row.labels.end(), label) == row.labels.end()) {
      speeds.insert(3 * label);
    }
  }
  for (int i = 0; i < 5; ++i) speeds.insert(crt_u(3, row.labels[i], unit_word[i]));
  assert(speeds.size() == 12);
  return std::vector<int>(speeds.begin(), speeds.end());
}

struct SheetCensus {
  long long scalar_rows = 0;
  std::set<Pattern> scalar_patterns;
  long long sheet_rows = 0;
  std::set<Pattern> sheet_patterns;
  long long compatible_unit_words = 0;
  std::vector<std::pair<Row, Orders>> least_packets;
  uint64_t sheet_hash = 1469598103934665603ULL;
};

SheetCensus finite_sheet_census(const std::vector<Row>& rows) {
  SheetCensus census;
  census.scalar_rows = static_cast<long long>(rows.size());
  for (const Row& row : rows) {
    census.scalar_patterns.insert(order_pattern(row.orders));
    const bool all_three = std::all_of(row.orders.begin(), row.orders.end(), [](int D) { return D == 3; });
    const SheetResult result = sheet_choices(row, all_three, false);
    if (result.compatible_choices == 0) continue;
    ++census.sheet_rows;
    census.sheet_patterns.insert(order_pattern(row.orders));
    census.compatible_unit_words += result.compatible_choices;
    for (int x : row.labels) census.sheet_hash = fnv_append(census.sheet_hash, x);
    for (int x : row.orders) census.sheet_hash = fnv_append(census.sheet_hash, x);
    census.sheet_hash = fnv_append(census.sheet_hash, result.compatible_choices);
    assert(all_three);
    for (const Orders& word : result.unit_words) census.least_packets.push_back({row, word});
  }
  return census;
}

std::array<int, 4> coset(int multiplier) {
  std::array<int, 4> answer{};
  for (int i = 0; i < 4; ++i) answer[i] = multiplier * H[i] % P;
  std::sort(answer.begin(), answer.end());
  return answer;
}

bool in_coset(int value, const std::array<int, 4>& C) {
  return std::find(C.begin(), C.end(), value) != C.end();
}

struct ClockAudit {
  int quartic_parity_rows = 0;
  int all_d3_flag_rows = 0;
  int mixed_d1_rows = 0;
  int flag_clock_min = 100;
  int flag_clock_max = -1;
  std::array<std::array<int, 3>, 3> tournament{};
};

ClockAudit clock_and_coset_audit() {
  ClockAudit audit;
  const std::array<int, 3> multipliers = {1, 2, 4};
  const std::array<std::array<int, 8>, 3> expected_clocks = {{
      {{2, 10, 11, 16, 23, 28, 29, 37}},
      {{1, 5, 8, 14, 25, 31, 34, 38}},
      {{4, 7, 17, 19, 20, 22, 32, 35}},
  }};

  std::array<std::array<int, 4>, 3> cosets{};
  for (int i = 0; i < 3; ++i) cosets[i] = coset(multipliers[i]);

  for (int ci = 0; ci < 3; ++ci) {
    const auto C = cosets[ci];
    const auto next = cosets[(ci + 1) % 3];
    std::vector<std::pair<int, int>> opposite_pairs;
    for (int r : C) {
      if (r < P - r && in_coset(P - r, C)) opposite_pairs.push_back({r, P - r});
    }
    assert(opposite_pairs.size() == 2);

    for (int bit0 = 1; bit0 <= 2; ++bit0) {
      for (int bit1 = 1; bit1 <= 2; ++bit1) {
        std::map<int, int> parity;
        parity[opposite_pairs[0].first] = parity[opposite_pairs[0].second] = bit0;
        parity[opposite_pairs[1].first] = parity[opposite_pairs[1].second] = bit1;
        std::vector<int> quartic;
        for (int label = 1; label < P; ++label) {
          if (!in_coset(label, C)) quartic.push_back(3 * label);
        }
        for (int r : C) quartic.push_back(crt_u(3, r, parity[r]));
        assert(quartic.size() == 12);

        std::vector<int> clock;
        for (int a = 1; a < 39; ++a) {
          if (std::gcd(a, 39) != 1) continue;
          int margin = 39;
          for (int speed : quartic) margin = std::min(margin, std::abs(signed_mod(1LL * a * speed, 39)));
          if (margin >= 3) {
            assert(margin == 3);
            clock.push_back(a);
          }
        }
        assert(clock.size() == 8);
        assert(std::equal(clock.begin(), clock.end(), expected_clocks[ci].begin()));
        ++audit.quartic_parity_rows;

        for (int b : next) {
          for (int e = 1; e <= 2; ++e) {
            std::vector<int> five = quartic;
            const auto it = std::find(five.begin(), five.end(), 3 * b);
            assert(it != five.end());
            for (int a : clock) {
              // The removed core coordinate is never one of the quartic
              // equality constraints on the forward flag.
              assert(std::abs(signed_mod(1LL * a * (3 * b), 39)) > 3);
            }
            five.erase(it);
            five.push_back(crt_u(3, b, e));
            int safe_clock = 0;
            for (int a : clock) {
              int margin = 39;
              for (int speed : five) margin = std::min(margin, std::abs(signed_mod(1LL * a * speed, 39)));
              if (margin >= 3) {
                assert(margin == 3);
                ++safe_clock;
              }
            }
            assert(safe_clock == 6);
            audit.flag_clock_min = std::min(audit.flag_clock_min, safe_clock);
            audit.flag_clock_max = std::max(audit.flag_clock_max, safe_clock);
            ++audit.all_d3_flag_rows;
          }
        }

        // In the D=1 + D=3^4 branch the actual D=1 replacement at common
        // scale 3 is 3u_b with u_b=b (mod 13), hence is 3b (mod 39).  It
        // exactly replaces the removed quartic core coordinate on the clock.
        for (int b = 1; b < P; ++b) {
          if (in_coset(b, C)) continue;
          std::vector<int> mixed = quartic;
          const auto it = std::find(mixed.begin(), mixed.end(), 3 * b);
          assert(it != mixed.end());
          mixed.erase(it);
          mixed.push_back(3 * (b + 13 * 7));  // arbitrary-height congruence check
          for (int a : clock) {
            int margin = 39;
            for (int speed : mixed) margin = std::min(margin, std::abs(signed_mod(1LL * a * speed, 39)));
            assert(margin == 3);
          }
          ++audit.mixed_d1_rows;
        }
      }
    }
  }

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const int owner = cosets[j][0];
      int hits = 0;
      for (int r : cosets[i]) hits += count_capacity(3, r, owner);
      audit.tournament[i][j] = hits;
    }
  }
  const std::array<std::array<int, 3>, 3> expected = {{{{3, 2, 0}}, {{0, 3, 2}}, {{2, 0, 3}}}};
  assert(audit.tournament == expected);
  return audit;
}

struct RunnerTournament {
  std::array<int, 5> scores{};
  int ties = 0;
  int triangles = 0;
  std::vector<int> scc_sizes;
  int hamiltonian_paths = 0;
};

RunnerTournament runner_tournament() {
  const Labels labels = {1, 2, 5, 8, 12};
  bool edge[5][5]{};
  RunnerTournament result;
  for (int i = 0; i < 5; ++i) {
    for (int j = i + 1; j < 5; ++j) {
      const int delta = count_capacity(3, labels[i], labels[j]) - count_capacity(3, labels[j], labels[i]);
      result.ties += delta == 0;
      if (delta > 0 || (delta == 0 && labels[i] < labels[j])) edge[i][j] = true;
      else edge[j][i] = true;
    }
  }
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) result.scores[i] += edge[i][j];
  }
  for (int i = 0; i < 5; ++i) {
    for (int j = i + 1; j < 5; ++j) {
      for (int k = j + 1; k < 5; ++k) {
        result.triangles += (edge[i][j] && edge[j][k] && edge[k][i]) ||
                            (edge[j][i] && edge[k][j] && edge[i][k]);
      }
    }
  }
  std::array<int, 5> permutation = {0, 1, 2, 3, 4};
  do {
    bool path = true;
    for (int i = 0; i < 4; ++i) path = path && edge[permutation[i]][permutation[i + 1]];
    result.hamiltonian_paths += path;
  } while (std::next_permutation(permutation.begin(), permutation.end()));

  bool reach[5][5]{};
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) reach[i][j] = i == j || edge[i][j];
  }
  for (int k = 0; k < 5; ++k) {
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j]);
    }
  }
  std::set<int> unseen = {0, 1, 2, 3, 4};
  while (!unseen.empty()) {
    const int i = *unseen.begin();
    std::vector<int> component;
    for (int j : unseen) {
      if (reach[i][j] && reach[j][i]) component.push_back(j);
    }
    result.scc_sizes.push_back(static_cast<int>(component.size()));
    for (int j : component) unseen.erase(j);
  }
  std::sort(result.scc_sizes.rbegin(), result.scc_sizes.rend());
  return result;
}

void verify_attenuation_and_counterfamily() {
  for (int D = 1; D <= 999; ++D) {
    if (D % P == 0) continue;
    std::array<int, 13> direct{};
    for (int z = -D + 1; z <= D; ++z) {
      int r = z % P;
      if (r < 0) r += P;
      ++direct[r];
    }
    for (int r = 1; r < P; ++r) {
      for (int o = 1; o < P; ++o) {
        const int target = D % P * r * inv13(o) % P;
        assert(count_capacity(D, r, o) == direct[target]);
        assert(std::abs(P * count_capacity(D, r, o) - 2 * D) < P);
      }
    }
  }

  const Labels labels = {1, 2, 3, 5, 10};
  for (int q = 3; q <= 1000; ++q) {
    const int D = P * q + 1;
    const Orders orders = {2, 5, D, 2, D};
    for (int owner : labels) {
      Fraction sum(0, 1);
      long long denominator = 1;
      for (int order : orders) denominator = std::lcm(denominator, static_cast<long long>(order));
      long long numerator = 0;
      for (int i = 0; i < 5; ++i) {
        numerator += count_capacity(orders[i], labels[i], owner) * (denominator / orders[i]);
      }
      assert(numerator >= denominator);
    }
  }
}

void verify_mixed_sheet_boundary() {
  for (int multiplier : {1, 2, 4}) {
    const auto C = coset(multiplier);
    std::vector<std::pair<int, int>> pairs;
    for (int r : C) {
      if (r < P - r) pairs.push_back({r, P - r});
    }
    assert(pairs.size() == 2);
    for (int b = 1; b < P; ++b) {
      if (in_coset(b, C)) continue;
      Labels labels = {C[0], C[1], C[2], C[3], b};
      Orders orders = {3, 3, 3, 3, 1};
      int compatible = 0;
      for (int mask = 0; mask < 16; ++mask) {
        Orders word{};
        for (int i = 0; i < 4; ++i) word[i] = 1 + ((mask >> i) & 1);
        word[4] = 1;
        bool okay = true;
        for (int owner = 0; owner < 5; ++owner) {
          std::array<bool, 3> covered{};
          for (int colour = 0; colour < 5; ++colour) {
            const int D = orders[colour];
            for (int ell0 = 0; ell0 < D; ++ell0) {
              if (!eligible_base_sheet(D, labels[colour], word[colour], labels[owner], ell0)) continue;
              for (int ell = ell0; ell < 3; ell += D) covered[ell] = true;
            }
          }
          okay = okay && std::all_of(covered.begin(), covered.end(), [](bool x) { return x; });
        }
        const auto index_of = [&](int label) {
          return static_cast<int>(std::find(labels.begin(), labels.end(), label) - labels.begin());
        };
        const bool pair_rule = word[index_of(pairs[0].first)] == word[index_of(pairs[0].second)] &&
                               word[index_of(pairs[1].first)] == word[index_of(pairs[1].second)];
        assert(okay == pair_rule);
        compatible += okay;
      }
      assert(compatible == 4);
    }
  }
}

}  // namespace

int main() {
  initialize_capacities();
  verify_attenuation_and_counterfamily();

  const ScalarCensus scalar = scalar_census();
  assert(scalar.total == 2410);
  assert(scalar.by_infinity[0] == 2190);
  assert(scalar.by_infinity[1] == 140);
  assert(scalar.by_infinity[2] == 80);
  assert(scalar.by_infinity[3] == 0 && scalar.by_infinity[4] == 0 && scalar.by_infinity[5] == 0);
  assert(scalar.patterns.size() == 86);
  assert(scalar.infinity_pattern_rows.size() == 8);

  const std::map<Pattern, std::tuple<long long, long long, Fraction>> expected_infinity = {
      {{{2, 2, 5, INF, INF}}, {50, 10, Fraction(1, 130)}},
      {{{2, 2, 10, INF, INF}}, {30, 6, Fraction(1, 130)}},
      {{{2, 2, 8, 9, INF}}, {60, 12, Fraction(1, 936)}},
      {{{2, 3, 3, 10, INF}}, {20, 4, Fraction(4, 195)}},
      {{{2, 3, 3, 11, INF}}, {20, 4, Fraction(1, 429)}},
      {{{3, 3, 3, 5, INF}}, {20, 4, Fraction(4, 195)}},
      {{{3, 3, 3, 10, INF}}, {10, 2, Fraction(4, 195)}},
      {{{3, 3, 3, 11, INF}}, {10, 2, Fraction(1, 429)}},
  };
  for (const auto& [pattern, expected] : expected_infinity) {
    const auto [rows, orbits, margin] = expected;
    assert(scalar.infinity_pattern_rows.at(pattern) == rows);
    assert(static_cast<long long>(scalar.infinity_orbits.at(pattern).size()) == orbits);
    const auto range = scalar.infinity_margin_range.at(pattern);
    assert(range.first == range.second);
    assert(Fraction(range.first, SCALE) == margin);
  }

  const SheetCensus sheets = finite_sheet_census(scalar.finite_rows);
  assert(sheets.scalar_rows == 2190);
  assert(sheets.scalar_patterns.size() == 78);
  assert(sheets.sheet_rows == 5);
  const std::set<Pattern> expected_sheet_patterns = {{{3, 3, 3, 3, 3}}};
  assert(sheets.sheet_patterns == expected_sheet_patterns);
  assert(sheets.compatible_unit_words == 40);
  assert(sheets.least_packets.size() == 40);

  std::map<std::pair<long long, long long>, int> maximum_histogram;
  Fraction metric_minimum(1, 1);
  int metric_minimum_count = 0;
  uint64_t metric_hash = 1469598103934665603ULL;
  for (const auto& [row, word] : sheets.least_packets) {
    const std::vector<int> speeds = packet_speeds(row, word);
    const ExactMaximum maximum = exact_maximum(speeds);
    assert(maximum.maximum > Fraction(1, 13));
    ++maximum_histogram[{maximum.maximum.n, maximum.maximum.d}];
    if (maximum.maximum < metric_minimum) {
      metric_minimum = maximum.maximum;
      metric_minimum_count = 1;
    } else if (maximum.maximum == metric_minimum) {
      ++metric_minimum_count;
    }
    for (int x : row.labels) metric_hash = fnv_append(metric_hash, x);
    for (int x : word) metric_hash = fnv_append(metric_hash, x);
    for (int x : speeds) metric_hash = fnv_append(metric_hash, x);
    metric_hash = fnv_append(metric_hash, maximum.maximum.n);
    metric_hash = fnv_append(metric_hash, maximum.maximum.d);
    metric_hash = fnv_append(metric_hash, maximum.first_witness.n);
    metric_hash = fnv_append(metric_hash, maximum.first_witness.d);
  }
  assert(metric_minimum == Fraction(2, 17));
  assert(metric_minimum_count == 4);
  const std::map<std::pair<long long, long long>, int> expected_histogram = {
      {{2, 17}, 4}, {{1, 8}, 6}, {{6, 43}, 8}, {{1, 7}, 4}, {{9, 46}, 4},
      {{6, 41}, 1}, {{3, 20}, 2}, {{5, 31}, 4}, {{6, 37}, 2}, {{11, 64}, 1}, {{1, 5}, 4},
  };
  assert(maximum_histogram == expected_histogram);

  const Row liar{{1, 2, 3, 5, 10}, {2, 5, 40, 2, 40}};
  const SheetResult liar_sheet = sheet_choices(liar, true, true);
  assert(liar_sheet.compatible_choices == 0);
  assert(liar_sheet.best_min_covered == 29);
  assert(liar_sheet.best_choice_count == 96);

  verify_mixed_sheet_boundary();
  const ClockAudit clock = clock_and_coset_audit();
  assert(clock.quartic_parity_rows == 12);
  assert(clock.all_d3_flag_rows == 96);
  assert(clock.mixed_d1_rows == 96);
  assert(clock.flag_clock_min == 6 && clock.flag_clock_max == 6);

  const RunnerTournament runners = runner_tournament();
  const std::array<int, 5> expected_runner_scores = {4, 1, 3, 2, 0};
  assert(runners.scores == expected_runner_scores);
  assert(runners.ties == 8 && runners.triangles == 0);
  assert(runners.scc_sizes == std::vector<int>({1, 1, 1, 1, 1}));
  assert(runners.hamiltonian_paths == 1);

  std::cout << "THM-823 exact replay\n";
  std::cout << "threshold.guardrail=1/13 arbitrary_lift_strictness=OPEN full_H5_closure=NO\n";
  std::cout << "attenuation.direct_D<=999=PASS error=|C_D-2/13|<1/D reciprocal=sum(1/D_i)>3/13 min_order<=21\n";
  std::cout << "scalar.counterfamily.labels=(1,2,3,5,10) orders=(2,5,13q+1,2,13q+1) q>=3 small_mass=7/10 large_floor=4q/(13q+1)\n";
  std::cout << "compact.rows=" << scalar.total << " by_inf=(" << scalar.by_infinity[0] << ','
            << scalar.by_infinity[1] << ',' << scalar.by_infinity[2] << ",0,0,0) patterns="
            << scalar.patterns.size() << " row_hash=" << hex64(scalar.row_hash) << "\n";
  for (const auto& [pattern, expected] : expected_infinity) {
    const auto [rows, orbits, margin] = expected;
    std::cout << "compact.pattern=" << pattern_string(pattern) << " rows=" << rows
              << " orbits=" << orbits << " margin=" << margin
              << " first=" << decorated_string(scalar.infinity_first.at(pattern)) << "\n";
  }
  std::cout << "compact.robustness=replace_inf_by_Dj_if_sum(1/Dj)<=pattern_margin\n";
  std::cout << "finite.no_D1.orders<=12.scalar_rows=" << sheets.scalar_rows
            << " scalar_patterns=" << sheets.scalar_patterns.size() << " sheet_rows="
            << sheets.sheet_rows << " sheet_patterns=" << sheets.sheet_patterns.size()
            << " unit_words=" << sheets.compatible_unit_words
            << " sheet_hash=" << hex64(sheets.sheet_hash) << "\n";
  std::cout << "finite.sheet_structure=C_union_{b}, C=a{1,5,8,12}, b_in_2C, quartet_opposites_equal, b_parity_free\n";
  std::cout << "least_CRT.packets=" << sheets.least_packets.size() << " min_M=" << metric_minimum
            << " min_count=" << metric_minimum_count << " histogram=";
  bool first = true;
  for (const auto& [key, count] : maximum_histogram) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << key.first << '/' << key.second << ':' << count;
  }
  std::cout << " metric_hash=" << hex64(metric_hash) << "\n";
  std::cout << "scalar_liar.D40.sheet_choices=1024 covers=0 best_min=29/40 best_count=96\n";
  std::cout << "D1.classification=all_D1_or_one_D1_plus_D3_coset common_sheet=quartet_pair_rule\n";
  std::cout << "clock.mod39.quartic_rows=" << clock.quartic_parity_rows
            << " D3_flag_tests=" << clock.all_d3_flag_rows << " safe_points_each="
            << clock.flag_clock_min << " mixed_D1_tests=" << clock.mixed_d1_rows
            << " conclusion=M>=1/13_only\n";
  std::cout << "tournament.coset.observable=A(Ci,Cj)=#D3_quartet_hits matrix=((3,2,0),(0,3,2),(2,0,3)) scores=(1,1,1) triangles=1 scc=(3) H=3\n";
  std::cout << "tournament.runner.scores=(4,1,3,2,0) ties=" << runners.ties
            << " triangles=" << runners.triangles << " scc=(1,1,1,1,1) H="
            << runners.hamiltonian_paths << "\n";
  std::cout << "carrier.audit=coset_cycle_preserves_D3_scalar_flags_but_loses_unit_parity_sheet_overlap_and_metric_lifts\n";
  std::cout << "scope.no_claim=unbounded_common_sheet_classification,arbitrary_lift_strictness,full_H5_closure\n";
}
