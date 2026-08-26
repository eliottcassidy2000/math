#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

// Independent exact C++ audit for THM-4166.
//
// Unlike the primary sequential interval-intersection implementation, this
// program constructs one global arrangement of every safe-comb wall needed
// for the fixed pool and 1 <= q <= 200.  It classifies each open atom by its
// unsafe optional labels.  Two-deletion masses are then sums of atom lengths.
// It also enumerates all C(27,20)=888030 twenty-vertex subsets for every q,
// independently checking the tau(Gamma_q)>7 admission boundary.

namespace {

using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_rational;
using Rat = cpp_rational;

constexpr std::array<int, 3> kAnchors = {120, 126, 143};
constexpr std::array<int, 30> kPool = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
};
constexpr std::array<int, 27> kOptional = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    132, 145, 168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
};
constexpr int kSearchLimit = 200;
constexpr std::uint32_t kAllVertices = (std::uint32_t{1} << 27U) - 1U;

struct Atom {
  Rat midpoint;
  Rat length;
  bool anchors_safe;
  std::uint32_t unsafe_optional;
};

struct GraphRecord {
  int q;
  std::array<std::uint32_t, 27> adjacency;
  int edge_count;
  int equality_count;
  int alpha;
  int tau;
  std::uint64_t subsets20_tested;
  std::uint64_t independent20_count;
};

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
  if (!condition) {
    fail(message);
  }
}

Rat make_rat(std::int64_t numerator_value, std::int64_t denominator_value) {
  return Rat(numerator_value) / Rat(denominator_value);
}

std::string rat_text(const Rat& value) {
  std::ostringstream output;
  output << numerator(value) << '/' << denominator(value);
  return output.str();
}

bool in_array(int value, const std::array<int, 30>& values) {
  return std::find(values.begin(), values.end(), value) != values.end();
}

bool safe_at(const Rat& x, int speed) {
  const Rat scaled = x * speed;
  const cpp_int den = denominator(scaled);
  cpp_int rem = numerator(scaled) % den;
  if (rem < 0) {
    rem += den;
  }
  return 14 * rem >= den && 14 * rem <= 13 * den;
}

std::vector<Rat> build_walls() {
  std::set<int> speeds;
  for (int speed = 1; speed <= kSearchLimit; ++speed) {
    speeds.insert(speed);
  }
  speeds.insert(kPool.begin(), kPool.end());

  std::vector<Rat> walls;
  walls.emplace_back(0);
  walls.emplace_back(1);
  for (const int speed : speeds) {
    for (int tooth = 0; tooth < speed; ++tooth) {
      walls.push_back(make_rat(14LL * tooth + 1, 14LL * speed));
      walls.push_back(make_rat(14LL * tooth + 13, 14LL * speed));
    }
  }
  std::sort(walls.begin(), walls.end());
  walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
  require(walls.front() == 0 && walls.back() == 1,
          "global wall endpoints changed");
  return walls;
}

std::vector<Atom> build_atoms(const std::vector<Rat>& walls) {
  std::vector<Atom> atoms;
  atoms.reserve(walls.size() - 1U);
  for (std::size_t index = 0; index + 1U < walls.size(); ++index) {
    const Rat left = walls[index];
    const Rat right = walls[index + 1U];
    require(left < right, "nonpositive global atom");
    const Rat midpoint = (left + right) / 2;
    bool anchors_safe = true;
    for (const int anchor : kAnchors) {
      anchors_safe = anchors_safe && safe_at(midpoint, anchor);
    }
    std::uint32_t unsafe = 0;
    for (std::size_t vertex = 0; vertex < kOptional.size(); ++vertex) {
      if (!safe_at(midpoint, kOptional[vertex])) {
        unsafe |= std::uint32_t{1} << vertex;
      }
    }
    atoms.push_back(Atom{midpoint, right - left, anchors_safe, unsafe});
  }
  return atoms;
}

bool selected_by_deletion_pair(const Atom& atom, int first, int second) {
  if (!atom.anchors_safe) {
    return false;
  }
  const std::uint32_t allowed =
      (std::uint32_t{1} << static_cast<unsigned>(first)) |
      (std::uint32_t{1} << static_cast<unsigned>(second));
  return (atom.unsafe_optional & ~allowed) == 0U;
}

std::array<std::array<Rat, 27>, 27> base_pair_masses(
    const std::vector<Atom>& atoms) {
  std::array<std::array<Rat, 27>, 27> result{};
  for (int first = 0; first < 27; ++first) {
    for (int second = first + 1; second < 27; ++second) {
      Rat total = 0;
      for (const Atom& atom : atoms) {
        if (selected_by_deletion_pair(atom, first, second)) {
          total += atom.length;
        }
      }
      result[static_cast<std::size_t>(first)][static_cast<std::size_t>(second)] = total;
      result[static_cast<std::size_t>(second)][static_cast<std::size_t>(first)] = total;
    }
  }
  return result;
}

std::array<std::array<int, 27>, 27> base_component_counts(
    const std::vector<Atom>& atoms) {
  std::array<std::array<int, 27>, 27> result{};
  for (int first = 0; first < 27; ++first) {
    for (int second = first + 1; second < 27; ++second) {
      int components = 0;
      bool previous = false;
      for (const Atom& atom : atoms) {
        const bool selected = selected_by_deletion_pair(atom, first, second);
        if (selected && !previous) {
          ++components;
        }
        previous = selected;
      }
      result[static_cast<std::size_t>(first)][static_cast<std::size_t>(second)] =
          components;
      result[static_cast<std::size_t>(second)][static_cast<std::size_t>(first)] =
          components;
    }
  }
  return result;
}

std::array<std::array<Rat, 27>, 27> q_pair_masses(
    const std::vector<Atom>& atoms, int q) {
  Rat mass0 = 0;
  std::array<Rat, 27> mass1{};
  std::array<std::array<Rat, 27>, 27> mass2{};
  for (const Atom& atom : atoms) {
    if (!atom.anchors_safe || !safe_at(atom.midpoint, q)) {
      continue;
    }
    const int unsafe_count = std::popcount(atom.unsafe_optional);
    if (unsafe_count == 0) {
      mass0 += atom.length;
    } else if (unsafe_count == 1) {
      const unsigned vertex = std::countr_zero(atom.unsafe_optional);
      mass1[vertex] += atom.length;
    } else if (unsafe_count == 2) {
      std::uint32_t rest = atom.unsafe_optional;
      const unsigned first = std::countr_zero(rest);
      rest &= rest - 1U;
      const unsigned second = std::countr_zero(rest);
      mass2[first][second] += atom.length;
      mass2[second][first] += atom.length;
    }
  }

  std::array<std::array<Rat, 27>, 27> result{};
  for (int first = 0; first < 27; ++first) {
    for (int second = first + 1; second < 27; ++second) {
      const Rat value = mass0 + mass1[static_cast<std::size_t>(first)] +
                        mass1[static_cast<std::size_t>(second)] +
                        mass2[static_cast<std::size_t>(first)]
                             [static_cast<std::size_t>(second)];
      result[static_cast<std::size_t>(first)][static_cast<std::size_t>(second)] = value;
      result[static_cast<std::size_t>(second)][static_cast<std::size_t>(first)] = value;
    }
  }
  return result;
}

class MaximumIndependentSet {
 public:
  explicit MaximumIndependentSet(const std::array<std::uint32_t, 27>& adjacency)
      : adjacency_(adjacency), best_(0) {}

  int run() {
    search(kAllVertices, 0);
    return best_;
  }

 private:
  void search(std::uint32_t candidates, int chosen) {
    if (chosen + std::popcount(candidates) <= best_) {
      return;
    }
    if (candidates == 0U) {
      best_ = std::max(best_, chosen);
      return;
    }
    unsigned selected = 0;
    int largest_degree = -1;
    std::uint32_t scan = candidates;
    while (scan != 0U) {
      const unsigned vertex = std::countr_zero(scan);
      scan &= scan - 1U;
      const int degree = std::popcount(adjacency_[vertex] & candidates);
      if (degree > largest_degree) {
        largest_degree = degree;
        selected = vertex;
      }
    }
    const std::uint32_t bit = std::uint32_t{1} << selected;
    search(candidates & ~bit & ~adjacency_[selected], chosen + 1);
    search(candidates & ~bit, chosen);
  }

  const std::array<std::uint32_t, 27>& adjacency_;
  int best_;
};

bool is_independent20(
    std::uint32_t subset,
    const std::array<std::uint32_t, 27>& adjacency) {
  std::uint32_t rest = subset;
  while (rest != 0U) {
    const unsigned vertex = std::countr_zero(rest);
    rest &= rest - 1U;
    if ((adjacency[vertex] & rest) != 0U) {
      return false;
    }
  }
  return true;
}

std::vector<std::uint32_t> all_twenty_subsets() {
  std::vector<std::uint32_t> answer;
  answer.reserve(888030U);
  std::uint64_t mask = (std::uint64_t{1} << 20U) - 1U;
  const std::uint64_t limit = std::uint64_t{1} << 27U;
  while (mask < limit) {
    answer.push_back(static_cast<std::uint32_t>(mask));
    const std::uint64_t smallest = mask & (~mask + 1U);
    const std::uint64_t ripple = mask + smallest;
    mask = ripple | (((mask ^ ripple) >> 2U) / smallest);
  }
  require(answer.size() == 888030U, "C(27,20) enumeration changed");
  for (const std::uint32_t subset : answer) {
    require(std::popcount(subset) == 20, "non-twenty-set in exhaustive bank");
  }
  return answer;
}

std::string semantic_row(const GraphRecord& row) {
  std::ostringstream output;
  output << "q=" << std::setw(3) << std::setfill('0') << row.q
         << ";edges=" << std::setw(3) << row.edge_count
         << ";eq=" << row.equality_count
         << ";alpha=" << std::setw(2) << row.alpha
         << ";tau=" << std::setw(2) << row.tau << ";adj=";
  output << std::hex << std::nouppercase;
  for (std::size_t vertex = 0; vertex < row.adjacency.size(); ++vertex) {
    if (vertex != 0U) {
      output << ',';
    }
    output << std::setw(7) << row.adjacency[vertex];
  }
  return output.str();
}

}  // namespace

int main() {
#ifdef _WIN32
  _setmode(_fileno(stdout), _O_BINARY);
#endif
  try {
    const Rat threshold = make_rat(4, 63);
    const Rat density = make_rat(6, 7);
    const Rat oscillation = make_rat(6, 49);
    const std::vector<Rat> walls = build_walls();
    const std::vector<Atom> atoms = build_atoms(walls);
    require(atoms.size() + 1U == walls.size(), "atom/wall count mismatch");

    const auto base_masses = base_pair_masses(atoms);
    const auto component_counts = base_component_counts(atoms);
    std::array<std::uint32_t, 27> limiting_adjacency{};
    int strict_limiting_edges = 0;
    int limiting_equalities = 0;
    bool have_worst_transient = false;
    Rat worst_transient = 0;
    int worst_first = -1;
    int worst_second = -1;
    for (int first = 0; first < 27; ++first) {
      for (int second = first + 1; second < 27; ++second) {
        const Rat limit_mass =
            density * base_masses[static_cast<std::size_t>(first)]
                                      [static_cast<std::size_t>(second)];
        if (limit_mass > threshold) {
          limiting_adjacency[static_cast<std::size_t>(first)] |=
              std::uint32_t{1} << static_cast<unsigned>(second);
          limiting_adjacency[static_cast<std::size_t>(second)] |=
              std::uint32_t{1} << static_cast<unsigned>(first);
          ++strict_limiting_edges;
        } else if (limit_mass == threshold) {
          ++limiting_equalities;
        } else {
          const Rat bound =
              oscillation *
              component_counts[static_cast<std::size_t>(first)]
                              [static_cast<std::size_t>(second)] /
              (threshold - limit_mass);
          if (!have_worst_transient || bound > worst_transient) {
            have_worst_transient = true;
            worst_transient = bound;
            worst_first = first;
            worst_second = second;
          }
        }
      }
    }
    require(have_worst_transient, "missing transient discrepancy bound");
    const int limiting_alpha = MaximumIndependentSet(limiting_adjacency).run();
    const int limiting_tau = 27 - limiting_alpha;
    require(strict_limiting_edges == 39, "strict limiting edge count changed");
    require(limiting_equalities == 0, "limiting equality appeared");
    require(limiting_alpha == 21 && limiting_tau == 6,
            "limiting graph alpha/tau changed");
    require(worst_first >= 0 && worst_second >= 0,
            "missing worst transient pair");

    const std::vector<std::uint32_t> subsets20 = all_twenty_subsets();
    std::vector<GraphRecord> rows;
    Rat smallest_positive = 0;
    Rat smallest_negative = 0;
    int positive_q = -1;
    int positive_first = -1;
    int positive_second = -1;
    int negative_q = -1;
    int negative_first = -1;
    int negative_second = -1;
    bool have_positive = false;
    bool have_negative = false;

    for (int q = 1; q <= kSearchLimit; ++q) {
      if (in_array(q, kPool)) {
        continue;
      }
      const auto masses = q_pair_masses(atoms, q);
      GraphRecord row{q, {}, 0, 0, 0, 0, 0, 0};
      for (int first = 0; first < 27; ++first) {
        for (int second = first + 1; second < 27; ++second) {
          const Rat difference =
              masses[static_cast<std::size_t>(first)]
                    [static_cast<std::size_t>(second)] -
              threshold;
          if (difference >= 0) {
            row.adjacency[static_cast<std::size_t>(first)] |=
                std::uint32_t{1} << static_cast<unsigned>(second);
            row.adjacency[static_cast<std::size_t>(second)] |=
                std::uint32_t{1} << static_cast<unsigned>(first);
            ++row.edge_count;
            if (difference == 0) {
              ++row.equality_count;
            } else if (!have_positive || difference < smallest_positive) {
              have_positive = true;
              smallest_positive = difference;
              positive_q = q;
              positive_first = first;
              positive_second = second;
            }
          } else if (!have_negative || -difference < smallest_negative) {
            have_negative = true;
            smallest_negative = -difference;
            negative_q = q;
            negative_first = first;
            negative_second = second;
          }
        }
      }
      row.alpha = MaximumIndependentSet(row.adjacency).run();
      row.tau = 27 - row.alpha;
      for (const std::uint32_t subset : subsets20) {
        ++row.subsets20_tested;
        if (is_independent20(subset, row.adjacency)) {
          ++row.independent20_count;
        }
      }
      require(row.subsets20_tested == 888030U,
              "not every twenty-set was tested");
      require((row.independent20_count == 0U) == (row.tau > 7),
              "twenty-set exhaustion disagrees with exact alpha/tau");
      rows.push_back(row);
    }

    require(rows.size() == 175U, "q census size changed");
    const int equality_total = std::accumulate(
        rows.begin(), rows.end(), 0,
        [](int total, const GraphRecord& row) {
          return total + row.equality_count;
        });
    require(equality_total == 0, "finite threshold equality appeared");

    std::map<int, int> tau_histogram;
    std::vector<int> admitted;
    std::uint64_t total_twenty_tests = 0;
    for (const GraphRecord& row : rows) {
      ++tau_histogram[row.tau];
      total_twenty_tests += row.subsets20_tested;
      if (row.tau > 7) {
        admitted.push_back(row.q);
      }
    }
    require(admitted.size() == 53U, "admitted q count changed");
    require(total_twenty_tests == 155405250U,
            "total twenty-set test count changed");

    std::cout << "THM4166_INDEPENDENT_GLOBAL_WALL_CPP_AUDIT\n";
    std::cout << "global_walls " << walls.size() << " atoms " << atoms.size()
              << "\n";
    std::cout << "q_universe 1..200 outside_P " << rows.size() << "\n";
    std::cout << "threshold " << rat_text(threshold) << " comparison >=\n";
    std::cout << "threshold_equalities " << equality_total << "\n";
    std::cout << "strict_limit_edges " << strict_limiting_edges
              << " limit_equalities " << limiting_equalities
              << " limit_alpha " << limiting_alpha
              << " limit_tau " << limiting_tau << "\n";
    std::cout << "worst_transient_bound " << rat_text(worst_transient)
              << " floor "
              << (numerator(worst_transient) / denominator(worst_transient))
              << " pair " << kOptional[static_cast<std::size_t>(worst_first)]
              << ',' << kOptional[static_cast<std::size_t>(worst_second)]
              << "\n";
    std::cout << "twenty_subsets_per_q " << subsets20.size()
              << " total_twenty_subset_tests " << total_twenty_tests << "\n";
    std::cout << "tau_histogram";
    for (const auto& [tau, count] : tau_histogram) {
      std::cout << ' ' << tau << ':' << count;
    }
    std::cout << "\nadmitted_count " << admitted.size() << "\nadmitted_q";
    for (const int q : admitted) {
      std::cout << ' ' << q;
    }
    std::cout << "\nsmallest_positive_margin " << rat_text(smallest_positive)
              << " q_r_s " << positive_q << ','
              << kOptional[static_cast<std::size_t>(positive_first)] << ','
              << kOptional[static_cast<std::size_t>(positive_second)] << "\n";
    std::cout << "smallest_negative_deficit " << rat_text(smallest_negative)
              << " q_r_s " << negative_q << ','
              << kOptional[static_cast<std::size_t>(negative_first)] << ','
              << kOptional[static_cast<std::size_t>(negative_second)] << "\n";
    std::cout << "EXHAUSTIVE_TWENTY_SET_ROWS_BEGIN\n";
    for (const GraphRecord& row : rows) {
      std::cout << "q=" << row.q << ";tested=" << row.subsets20_tested
                << ";independent20=" << row.independent20_count << "\n";
    }
    std::cout << "EXHAUSTIVE_TWENTY_SET_ROWS_END\n";
    std::cout << "SEMANTIC_LEDGER_BEGIN\n";
    for (const GraphRecord& row : rows) {
      std::cout << semantic_row(row) << '\n';
    }
    std::cout << "SEMANTIC_LEDGER_END\n";
    std::cout << "INDEPENDENT_CPP_AUDIT_PASS\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "AUDIT_FAILURE: " << error.what() << '\n';
    return 1;
  }
}
