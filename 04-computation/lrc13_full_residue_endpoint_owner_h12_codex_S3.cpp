// Exact endpoint-owner census for full nonzero residue packets modulo 13.
//
// The theorem-facing bank is
//
//   W(k) = {r + 13*k_r : 1 <= r <= 12},   0 <= k_r <= 12.
//
// It contains 13^12 = 23,298,085,122,481 labelled packets.  A literal loop
// over that bank is unnecessary.  All possible endpoints of the open danger
// teeth ||w t|| < 1/13 cut the circle into finitely many open atomic cells.
// Every candidate runner is either dangerous on a whole cell or on none of
// it.  Moreover a strict witness is an open condition, so it exists exactly
// when one atomic cell is uncovered.  This turns zero endpoint defect into an
// exact grouped set-cover problem: choose one speed in each nonzero residue
// class and cover every cell.
//
// The DFS uses an exact owner rule.  At a partial assignment, a cell which no
// remaining residue group except r can possibly cover must be covered by the
// eventual choice in group r.  Options failing any such uniquely owned cell
// are impossible.  Branching on the group with the fewest surviving options
// is lossless; it is not a heuristic certificate.
//
// Build and run (about two minutes on the reference Apple workstation):
//
//   c++ -O3 -std=c++17 \
//     04-computation/lrc13_full_residue_endpoint_owner_h12_codex_S3.cpp \
//     -o /tmp/lrc13_full_residue_endpoint_owner_h12
//   /tmp/lrc13_full_residue_endpoint_owner_h12
//
// All comparisons use integers (__int128 cross-products where needed).  No
// floating-point arithmetic, numerical maximization, or unproved q cutoff is
// used.  The built-in height-two literal endpoint sweep independently checks
// the cell-cover reduction on all 3^12 = 531,441 packets.
//
// Tournament Analysis.  The diagnostic vertices are the twelve residue-choice
// obligations, not runners.  For a pair r,s, erase the possible coverage of
// the other ten groups and compare the numbers of remaining cells exclusively
// ownable by r and by s.  The sign is the switch/gauge; residue order
// 1->2->...->12 is the tie Hamiltonian path.  The resulting score histogram,
// directed triangles, SCCs, edge flips, and Hamiltonian-path count are printed.
// This tournament preserves pairwise owner pressure but destroys simultaneous
// option compatibility.  The exact proof carrier is the cell/choice incidence
// hypergraph and its unique-owner propagation, not the tournament quotient.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

constexpr int Q = 13;
constexpr int RESIDUES = 12;
constexpr int HEIGHT = 12;

struct Rational {
  std::int64_t numerator;
  std::int64_t denominator;
};

bool rational_less(const Rational& a, const Rational& b) {
  return static_cast<__int128>(a.numerator) * b.denominator <
         static_cast<__int128>(b.numerator) * a.denominator;
}

bool rational_equal(const Rational& a, const Rational& b) {
  return static_cast<__int128>(a.numerator) * b.denominator ==
         static_cast<__int128>(b.numerator) * a.denominator;
}

Rational midpoint_on_universal_cover(Rational a, Rational b) {
  const std::int64_t common_gcd = std::gcd(a.denominator, b.denominator);
  const __int128 common_denominator =
      static_cast<__int128>(a.denominator / common_gcd) * b.denominator;
  const __int128 numerator =
      static_cast<__int128>(a.numerator) *
          (common_denominator / a.denominator) +
      static_cast<__int128>(b.numerator) *
          (common_denominator / b.denominator);
  const __int128 denominator = 2 * common_denominator;

  assert(numerator >= 0);
  assert(numerator <= std::numeric_limits<std::int64_t>::max());
  assert(denominator <= std::numeric_limits<std::int64_t>::max());
  const std::int64_t n = static_cast<std::int64_t>(numerator);
  const std::int64_t d = static_cast<std::int64_t>(denominator);
  const std::int64_t g = std::gcd(n, d);
  return {n / g, d / g};
}

using Bits = std::vector<std::uint64_t>;

void set_bit(Bits& bits, std::size_t index) {
  bits[index >> 6] |= std::uint64_t{1} << (index & 63);
}

bool bits_empty(const Bits& bits) {
  for (std::uint64_t word : bits) {
    if (word != 0) return false;
  }
  return true;
}

bool bits_subset(const Bits& a, const Bits& b) {
  assert(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    if ((a[i] & ~b[i]) != 0) return false;
  }
  return true;
}

struct EndpointEvent {
  Rational time;
  int starts;
  int ends;
};

struct EndpointStats {
  int kappa;
  int protected_splices;
  int defect;
};

EndpointStats literal_endpoint_stats(
    const std::array<int, RESIDUES>& speeds) {
  std::vector<EndpointEvent> events;
  int sum_speeds = 0;
  for (int w : speeds) sum_speeds += w;
  events.reserve(2 * sum_speeds);

  for (int w : speeds) {
    const int denominator = Q * w;
    for (int m = 0; m < w; ++m) {
      int start = Q * m - 1;
      if (start < 0) start += denominator;
      events.push_back({{start, denominator}, 1, 0});
      events.push_back({{Q * m + 1, denominator}, 0, 1});
    }
  }

  std::sort(events.begin(), events.end(), [](const EndpointEvent& a,
                                              const EndpointEvent& b) {
    return rational_less(a.time, b.time);
  });

  // At t=0 the zero-centred tooth of every speed is active.
  int depth = RESIDUES;
  int kappa = 0;
  int protected_splices = 0;
  int strict_safe_components = 0;

  for (std::size_t i = 0; i < events.size();) {
    std::size_t j = i + 1;
    int starts = events[i].starts;
    int ends = events[i].ends;
    while (j < events.size() &&
           rational_equal(events[i].time, events[j].time)) {
      starts += events[j].starts;
      ends += events[j].ends;
      ++j;
    }

    const int before = depth;
    const int continuing = before - ends;
    const int after = continuing + starts;
    if (continuing == 0 && after > 0) {
      ++kappa;
      if (before > 0) ++protected_splices;
    }
    if (before > 0 && after == 0) ++strict_safe_components;
    depth = after;
    i = j;
  }

  assert(depth == RESIDUES);
  assert(kappa - protected_splices == strict_safe_components);
  return {kappa, protected_splices, strict_safe_components};
}

int common_gcd(const std::array<int, RESIDUES>& speeds) {
  int g = 0;
  for (int w : speeds) g = std::gcd(g, w);
  return g;
}

bool dilated_ap_scale(const std::array<int, RESIDUES>& labelled,
                      int* scale) {
  std::array<int, RESIDUES> sorted = labelled;
  std::sort(sorted.begin(), sorted.end());
  const int c = sorted.front();
  for (int i = 0; i < RESIDUES; ++i) {
    if (sorted[i] != c * (i + 1)) return false;
  }
  if (scale != nullptr) *scale = c;
  return true;
}

std::uint64_t integer_power(std::uint64_t base, int exponent) {
  std::uint64_t result = 1;
  for (int i = 0; i < exponent; ++i) result *= base;
  return result;
}

struct SearchStats {
  std::uint64_t nodes = 0;
  std::uint64_t prunes = 0;
  std::uint64_t leaves = 0;
  std::uint64_t uncovered_leaves = 0;
};

class EndpointOwnerSolver {
 public:
  explicit EndpointOwnerSolver(int height) : height_(height) {
    build_choices();
    build_atomic_cells();
    build_coverage();
  }

  void solve() {
    Bits uncovered(word_count_, ~std::uint64_t{0});
    if (cell_count_ % 64 != 0) {
      uncovered.back() =
          (std::uint64_t{1} << (cell_count_ % 64)) - 1;
    }
    dfs((std::uint16_t{1} << RESIDUES) - 1, uncovered);
  }

  int height() const { return height_; }
  int endpoint_count() const { return endpoint_count_; }
  int cell_count() const { return cell_count_; }
  int word_count() const { return word_count_; }
  const SearchStats& stats() const { return stats_; }
  const std::vector<std::array<int, RESIDUES>>& solutions() const {
    return solutions_;
  }

  void print_tournament_analysis() const {
    bool adjacency[RESIDUES][RESIDUES] = {};
    int edge_flips = 0;
    int burden_ties = 0;
    std::uint64_t forced_owner_checksum = 0;

    for (int r = 0; r < RESIDUES; ++r) {
      for (int s = r + 1; s < RESIDUES; ++s) {
        Bits other(word_count_, 0);
        for (int h = 0; h < RESIDUES; ++h) {
          if (h == r || h == s) continue;
          for (int z = 0; z < word_count_; ++z) {
            other[z] |= group_union_[h][z];
          }
        }

        std::uint64_t forced_r = 0;
        std::uint64_t forced_s = 0;
        for (int z = 0; z < word_count_; ++z) {
          const std::uint64_t residual = ~other[z];
          forced_r += __builtin_popcountll(
              residual & group_union_[r][z] & ~group_union_[s][z]);
          forced_s += __builtin_popcountll(
              residual & group_union_[s][z] & ~group_union_[r][z]);
        }
        forced_owner_checksum +=
            static_cast<std::uint64_t>(r + 1) * 1000003ULL * forced_r +
            static_cast<std::uint64_t>(s + 1) * 1000033ULL * forced_s;

        int winner;
        int loser;
        if (forced_r > forced_s || forced_r == forced_s) {
          winner = r;
          loser = s;
          if (forced_r == forced_s) ++burden_ties;
        } else {
          winner = s;
          loser = r;
          ++edge_flips;
        }
        adjacency[winner][loser] = true;
      }
    }

    std::array<int, RESIDUES> scores{};
    std::map<int, int> score_histogram;
    for (int r = 0; r < RESIDUES; ++r) {
      for (int s = 0; s < RESIDUES; ++s) {
        if (adjacency[r][s]) ++scores[r];
      }
      ++score_histogram[scores[r]];
    }

    int directed_triangles = 0;
    for (int a = 0; a < RESIDUES; ++a) {
      for (int b = a + 1; b < RESIDUES; ++b) {
        for (int c = b + 1; c < RESIDUES; ++c) {
          if ((adjacency[a][b] && adjacency[b][c] && adjacency[c][a]) ||
              (adjacency[a][c] && adjacency[c][b] && adjacency[b][a])) {
            ++directed_triangles;
          }
        }
      }
    }

    bool reach[RESIDUES][RESIDUES] = {};
    for (int r = 0; r < RESIDUES; ++r) {
      for (int s = 0; s < RESIDUES; ++s) {
        reach[r][s] = r == s || adjacency[r][s];
      }
    }
    for (int k = 0; k < RESIDUES; ++k) {
      for (int r = 0; r < RESIDUES; ++r) {
        for (int s = 0; s < RESIDUES; ++s) {
          reach[r][s] = reach[r][s] || (reach[r][k] && reach[k][s]);
        }
      }
    }
    std::vector<int> scc_sizes;
    std::uint16_t unseen = (std::uint16_t{1} << RESIDUES) - 1;
    while (unseen != 0) {
      const int r = __builtin_ctz(unseen);
      std::uint16_t component = 0;
      for (int s = 0; s < RESIDUES; ++s) {
        if ((unseen & (std::uint16_t{1} << s)) && reach[r][s] &&
            reach[s][r]) {
          component |= std::uint16_t{1} << s;
        }
      }
      scc_sizes.push_back(__builtin_popcount(component));
      unseen &= ~component;
    }
    std::sort(scc_sizes.begin(), scc_sizes.end(), std::greater<int>());

    // Exact number of directed Hamiltonian paths by subset DP.
    std::vector<std::array<std::uint64_t, RESIDUES>> dp(
        std::size_t{1} << RESIDUES);
    for (int r = 0; r < RESIDUES; ++r) {
      dp[std::size_t{1} << r][r] = 1;
    }
    for (std::size_t mask = 1; mask < dp.size(); ++mask) {
      for (int last = 0; last < RESIDUES; ++last) {
        if (dp[mask][last] == 0) continue;
        for (int next = 0; next < RESIDUES; ++next) {
          if ((mask & (std::size_t{1} << next)) == 0 &&
              adjacency[last][next]) {
            dp[mask | (std::size_t{1} << next)][next] += dp[mask][last];
          }
        }
      }
    }
    std::uint64_t hamiltonian_paths = 0;
    for (int last = 0; last < RESIDUES; ++last) {
      hamiltonian_paths += dp.back()[last];
    }

    std::cout << "TOURNAMENT ANALYSIS (DIAGNOSTIC OWNER-PRESSURE QUOTIENT)\n";
    std::cout << "vertices=residue-choice obligations 1..12\n";
    std::cout << "observable=forced owner cells after deleting the other ten "
                 "groups\n";
    std::cout << "switch/gauge=sign(forced_r-forced_s)\n";
    std::cout << "tie_hamiltonian_path=1->2->...->12\n";
    std::cout << "score_histogram={";
    bool first = true;
    for (const auto& [score, count] : score_histogram) {
      if (!first) std::cout << ", ";
      std::cout << score << ":" << count;
      first = false;
    }
    std::cout << "}\n";
    std::cout << "directed_3cycles=" << directed_triangles
              << " scc_sizes=[";
    for (std::size_t i = 0; i < scc_sizes.size(); ++i) {
      if (i != 0) std::cout << ",";
      std::cout << scc_sizes[i];
    }
    std::cout << "] edge_flips=" << edge_flips
              << " burden_ties=" << burden_ties
              << " hamiltonian_paths=" << hamiltonian_paths << "\n";
    std::cout << "forced_owner_checksum=" << forced_owner_checksum << "\n";
    std::cout << "guardrail=pairwise owner pressure destroys simultaneous "
                 "choice compatibility\n";
  }

 private:
  void build_choices() {
    for (int r = 1; r <= RESIDUES; ++r) {
      for (int k = 0; k <= height_; ++k) {
        choices_[r - 1].push_back(r + Q * k);
      }
    }
  }

  void build_atomic_cells() {
    std::vector<Rational> endpoints;
    for (int r = 0; r < RESIDUES; ++r) {
      for (int w : choices_[r]) {
        const int denominator = Q * w;
        for (int m = 0; m < w; ++m) {
          int start = Q * m - 1;
          if (start < 0) start += denominator;
          endpoints.push_back({start, denominator});
          endpoints.push_back({Q * m + 1, denominator});
        }
      }
    }
    std::sort(endpoints.begin(), endpoints.end(), rational_less);
    endpoints.erase(
        std::unique(endpoints.begin(), endpoints.end(), rational_equal),
        endpoints.end());
    endpoint_count_ = static_cast<int>(endpoints.size());

    cell_samples_.resize(endpoints.size());
    for (std::size_t i = 0; i < endpoints.size(); ++i) {
      Rational a = endpoints[i];
      Rational b = endpoints[(i + 1) % endpoints.size()];
      if (i + 1 == endpoints.size()) b.numerator += b.denominator;
      cell_samples_[i] = midpoint_on_universal_cover(a, b);
    }
    cell_count_ = static_cast<int>(cell_samples_.size());
    word_count_ = (cell_count_ + 63) / 64;
  }

  void build_coverage() {
    for (int r = 0; r < RESIDUES; ++r) {
      group_union_[r] = Bits(word_count_, 0);
      for (int w : choices_[r]) {
        Bits covered(word_count_, 0);
        for (int cell = 0; cell < cell_count_; ++cell) {
          const Rational t = cell_samples_[cell];
          std::int64_t residue = static_cast<std::int64_t>(
              (static_cast<__int128>(w) * t.numerator) % t.denominator);
          if (residue < 0) residue += t.denominator;
          const std::int64_t distance =
              std::min(residue, t.denominator - residue);
          if (static_cast<__int128>(Q) * distance < t.denominator) {
            set_bit(covered, cell);
          }
        }
        for (int z = 0; z < word_count_; ++z) {
          group_union_[r][z] |= covered[z];
        }
        coverage_[r].push_back(std::move(covered));
      }
    }
  }

  void dfs(std::uint16_t remaining_groups, const Bits& uncovered) {
    ++stats_.nodes;
    if (remaining_groups == 0) {
      ++stats_.leaves;
      if (bits_empty(uncovered)) {
        solutions_.push_back(picked_);
      } else {
        ++stats_.uncovered_leaves;
      }
      return;
    }

    int best_group = -1;
    std::vector<int> best_options;

    for (int r = 0; r < RESIDUES; ++r) {
      if ((remaining_groups & (std::uint16_t{1} << r)) == 0) continue;

      Bits coverage_by_others(word_count_, 0);
      for (int h = 0; h < RESIDUES; ++h) {
        if (h == r ||
            (remaining_groups & (std::uint16_t{1} << h)) == 0) {
          continue;
        }
        for (int z = 0; z < word_count_; ++z) {
          coverage_by_others[z] |= group_union_[h][z];
        }
      }

      Bits uniquely_owned(word_count_, 0);
      for (int z = 0; z < word_count_; ++z) {
        uniquely_owned[z] = uncovered[z] & ~coverage_by_others[z];
      }

      std::vector<int> viable_options;
      for (int option = 0; option <= height_; ++option) {
        if (bits_subset(uniquely_owned, coverage_[r][option])) {
          viable_options.push_back(option);
        }
      }
      if (viable_options.empty()) {
        ++stats_.prunes;
        return;
      }
      if (best_group < 0 || viable_options.size() < best_options.size()) {
        best_group = r;
        best_options = std::move(viable_options);
        if (best_options.size() == 1) break;
      }
    }

    assert(best_group >= 0);
    for (int option : best_options) {
      Bits next_uncovered(word_count_, 0);
      for (int z = 0; z < word_count_; ++z) {
        next_uncovered[z] =
            uncovered[z] & ~coverage_[best_group][option][z];
      }
      picked_[best_group] = choices_[best_group][option];
      dfs(remaining_groups ^ (std::uint16_t{1} << best_group),
          next_uncovered);
    }
  }

  int height_;
  int endpoint_count_ = 0;
  int cell_count_ = 0;
  int word_count_ = 0;
  std::array<std::vector<int>, RESIDUES> choices_;
  std::vector<Rational> cell_samples_;
  std::array<std::vector<Bits>, RESIDUES> coverage_;
  std::array<Bits, RESIDUES> group_union_;
  std::array<int, RESIDUES> picked_{};
  SearchStats stats_;
  std::vector<std::array<int, RESIDUES>> solutions_;
};

std::vector<std::array<int, RESIDUES>> literal_height_two_zeros(
    std::uint64_t* primitive_non_ap_count,
    int* primitive_non_ap_min_defect) {
  constexpr int small_height = 2;
  const std::uint64_t patterns =
      integer_power(small_height + 1, RESIDUES);
  std::vector<std::array<int, RESIDUES>> zeros;
  *primitive_non_ap_count = 0;
  *primitive_non_ap_min_defect = std::numeric_limits<int>::max();

  for (std::uint64_t code = 0; code < patterns; ++code) {
    std::uint64_t digits = code;
    std::array<int, RESIDUES> speeds{};
    for (int r = 1; r <= RESIDUES; ++r) {
      const int k = digits % (small_height + 1);
      digits /= small_height + 1;
      speeds[r - 1] = r + Q * k;
    }
    const EndpointStats stats = literal_endpoint_stats(speeds);
    if (stats.defect == 0) zeros.push_back(speeds);

    int scale = 0;
    const bool is_ap = dilated_ap_scale(speeds, &scale);
    if (common_gcd(speeds) == 1 && !is_ap) {
      ++*primitive_non_ap_count;
      *primitive_non_ap_min_defect =
          std::min(*primitive_non_ap_min_defect, stats.defect);
    }
  }
  return zeros;
}

void sort_by_ap_scale(std::vector<std::array<int, RESIDUES>>* rows) {
  std::sort(rows->begin(), rows->end(), [](const auto& a, const auto& b) {
    int ca = 0;
    int cb = 0;
    const bool apa = dilated_ap_scale(a, &ca);
    const bool apb = dilated_ap_scale(b, &cb);
    if (apa != apb) return apa > apb;
    if (apa && ca != cb) return ca < cb;
    return a < b;
  });
}

void print_labelled_row(const std::array<int, RESIDUES>& row) {
  std::cout << "[";
  for (int i = 0; i < RESIDUES; ++i) {
    if (i != 0) std::cout << ",";
    std::cout << row[i];
  }
  std::cout << "]";
}

}  // namespace

int main() {
  std::cout << std::string(88, '=') << "\n";
  std::cout << "LRC(13) FULL-RESIDUE ENDPOINT-OWNER H=12 CENSUS "
               "(codex-2026-07-14-S3)\n";
  std::cout << std::string(88, '=') << "\n";
  std::cout << "arithmetic=integer exact (__int128 rational comparisons); "
               "no floating point\n";
  std::cout << "predicate=zero strict 1/13-safe components "
               "(equivalently endpoint defect chi_13=0)\n";
  std::cout << "finite_carrier=open atomic endpoint cells x residue-choice "
               "options\n";
  std::cout << "owner_rule=a cell coverable by only one remaining residue "
               "group forces that owner\n\n";

  std::cout << "INDEPENDENT LITERAL ENDPOINT-SWEEP AUDIT\n";
  std::uint64_t primitive_non_ap_count = 0;
  int primitive_non_ap_min_defect = 0;
  auto literal_zeros = literal_height_two_zeros(
      &primitive_non_ap_count, &primitive_non_ap_min_defect);
  sort_by_ap_scale(&literal_zeros);

  EndpointOwnerSolver height_two_solver(2);
  height_two_solver.solve();
  auto cell_zeros = height_two_solver.solutions();
  sort_by_ap_scale(&cell_zeros);
  assert(literal_zeros == cell_zeros);
  assert(literal_zeros.size() == 3);
  assert(primitive_non_ap_min_defect == 2);
  std::cout << "height=2 patterns=531441 literal_zero_rows="
            << literal_zeros.size()
            << " cell_CSP_zero_rows=" << cell_zeros.size()
            << " classifications_agree=true\n";
  std::cout << "primitive_non_AP_rows=" << primitive_non_ap_count
            << " primitive_non_AP_min_defect="
            << primitive_non_ap_min_defect << "\n\n";

  std::cout << "HEIGHT-TWELVE EXACT OWNER-CSP\n";
  EndpointOwnerSolver solver(HEIGHT);
  std::cout << "height=" << HEIGHT
            << " conceptual_patterns="
            << integer_power(HEIGHT + 1, RESIDUES)
            << " distinct_endpoints=" << solver.endpoint_count()
            << " atomic_cells=" << solver.cell_count()
            << " bit_words=" << solver.word_count() << "\n";
  solver.solve();
  auto solutions = solver.solutions();
  sort_by_ap_scale(&solutions);

  const SearchStats& search = solver.stats();
  std::cout << "search_nodes=" << search.nodes
            << " prunes=" << search.prunes
            << " complete_leaves=" << search.leaves
            << " uncovered_complete_leaves=" << search.uncovered_leaves
            << " zero_rows=" << solutions.size() << "\n";

  assert(solutions.size() == 13);
  assert(search.uncovered_leaves == 0);
  int primitive_zero_rows = 0;
  std::vector<int> scales;
  for (const auto& row : solutions) {
    int scale = 0;
    const bool is_ap = dilated_ap_scale(row, &scale);
    assert(is_ap);
    scales.push_back(scale);
    if (common_gcd(row) == 1) ++primitive_zero_rows;

    const EndpointStats endpoint = literal_endpoint_stats(row);
    assert(endpoint.defect == 0);
    assert(endpoint.kappa == RESIDUES * scale);
    assert(endpoint.protected_splices == RESIDUES * scale);
    std::cout << "  scale=" << scale << " gcd=" << common_gcd(row)
              << " kappa=" << endpoint.kappa
              << " protected=" << endpoint.protected_splices
              << " defect=" << endpoint.defect << " labelled=";
    print_labelled_row(row);
    std::cout << "\n";
  }
  const std::vector<int> expected_scales =
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14};
  assert(scales == expected_scales);
  assert(primitive_zero_rows == 1);
  std::cout << "classification=zero defect iff W=c*{1,...,12} in the box\n";
  std::cout << "scales=[1,2,3,4,5,6,7,8,9,10,11,12,14]\n";
  std::cout << "primitive_zero_rows=" << primitive_zero_rows
            << " primitive_row={1,...,12}\n\n";

  solver.print_tournament_analysis();
  std::cout << "\nPRESERVE / DESTROY\n";
  std::cout << "atomic cell incidence: preserves the exact finite strict-witness "
               "predicate\n";
  std::cout << "unique-owner propagation: preserves all completions; destroys "
               "only impossible options\n";
  std::cout << "owner-pressure tournament: preserves pairwise burden; destroys "
               "joint option compatibility\n\n";
  std::cout << "SCOPE\n";
  std::cout << "This is a finite-exact theorem for 0<=k_r<=12, not the unbounded "
               "primitive rigidity theorem.\n";
  std::cout << "Dilation/gcd descent reduces the unbounded problem to primitive "
               "zero-defect packets, but does not classify them.\n";
  return 0;
}
