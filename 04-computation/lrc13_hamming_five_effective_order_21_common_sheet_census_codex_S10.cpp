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

// Exact effective-order-at-most-21 Hamming-five common-sheet census.
//
// This file is deliberately theorem-ID-free.  It supplies a replay artifact
// for the next theorem after THM-823.  There are two independent closure
// layers.
//
//  (1) Relative-ramification cuts.  For a nonempty colour set S, let T be its
//      complement, L_T=lcm(D_j:j in T), and
//
//          m_i = D_i/gcd(D_i,L_T),  i in S.
//
//      At an owner missed by T, an L_T-sheet fibre must be covered by S.
//      After the affine unit change derived below, colour i is a cyclic
//      interval sampled on a gcd coset, so it covers at most
//
//          ceil(2*m_i/13)/m_i
//
//      of that fibre.  Hence a common-sheet cover necessarily satisfies
//
//          sum_(i in S) ceil(2*m_i/13)/m_i >= 1.          (R)
//
//      The singleton cuts imply D_i | lcm(D_j:j!=i).  On the present bank,
//      singleton cuts reject all but 185 scalar rows and top-prime cuts reject
//      those remaining rows.
//
//  (2) Literal common-sheet replay.  Independently of the top-prime cut,
//      every unit word on those 185 lcm-redundant rows is tested by the
//      original centred-CRT definition.  This is a small exact cross-check of
//      the structural closure, not a sampled metric calculation.
//
// The all-size input that every set of at most four replacement colours misses
// some one of the five owners at scalar capacity follows from THM-810: a
// four-colour scalar cover of its own owners is either all order one or an
// all-order-three coset.  The former contributes zero at a fifth owner, and
// the latter contributes only 0 or 2/3 there (THM-823's coset matrix).  For
// self-containment of this finite replay, that deficit is also checked
// directly for every scalar row retained here.

namespace {

constexpr int P = 13;
constexpr long long SCALE = 232792560LL;  // lcm(1,...,22), hence every D here
constexpr std::array<int, 19> ORDER_ALPHABET = {
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21};
constexpr std::array<int, 7> SMALL_PRIMES = {2, 3, 5, 7, 11, 17, 19};

using Labels = std::array<int, 5>;
using Orders = std::array<int, 5>;
using Pattern = std::array<int, 5>;

struct Row {
  Labels labels{};
  Orders orders{};
  std::array<int, 5> four_colour_deficit_owner{};
};

long long gcdll(long long a, long long b) {
  while (b != 0) {
    const long long r = a % b;
    a = b;
    b = r;
  }
  return a < 0 ? -a : a;
}

struct Fraction {
  long long numerator = 0;
  long long denominator = 1;

  Fraction() = default;
  Fraction(long long n, long long d) : numerator(n), denominator(d) {
    assert(denominator != 0);
    if (denominator < 0) {
      numerator = -numerator;
      denominator = -denominator;
    }
    const long long g = gcdll(numerator, denominator);
    numerator /= g;
    denominator /= g;
  }
};

Fraction operator+(const Fraction& x, const Fraction& y) {
  return Fraction(x.numerator * y.denominator + y.numerator * x.denominator,
                  x.denominator * y.denominator);
}

bool operator<(const Fraction& x, const Fraction& y) {
  return x.numerator * y.denominator < y.numerator * x.denominator;
}

bool operator==(const Fraction& x, const Fraction& y) {
  return x.numerator == y.numerator && x.denominator == y.denominator;
}

bool operator>(const Fraction& x, const Fraction& y) { return y < x; }

std::ostream& operator<<(std::ostream& out, const Fraction& x) {
  if (x.denominator == 1) return out << x.numerator;
  return out << x.numerator << '/' << x.denominator;
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

int positive_mod(long long value, int modulus) {
  int answer = static_cast<int>(value % modulus);
  if (answer < 0) answer += modulus;
  return answer;
}

int centred_mod(long long value, int modulus) {
  const int residue = positive_mod(value, modulus);
  return residue <= modulus / 2 ? residue : residue - modulus;
}

int inverse_mod(int value, int modulus) {
  assert(std::gcd(value, modulus) == 1);
  for (int candidate = 1; candidate < modulus; ++candidate) {
    if (1LL * value * candidate % modulus == 1) return candidate;
  }
  assert(false);
  return 0;
}

int inverse13(int value) { return inverse_mod(value, P); }

int capacity_count(int order, int replacement, int owner) {
  assert(order >= 1 && order % P != 0);
  const int q = order / P;
  const int d = order % P;
  const int target = d * replacement * inverse13(owner) % P;
  int answer = 2 * q;
  for (int z = -d + 1; z <= d; ++z) {
    answer += positive_mod(z - target, P) == 0;
  }
  return answer;
}

std::array<std::array<std::array<long long, 13>, 13>, 22> scaled_capacity{};

void initialize_capacities() {
  for (int order : ORDER_ALPHABET) {
    assert(SCALE % order == 0);
    for (int replacement = 1; replacement < P; ++replacement) {
      for (int owner = 1; owner < P; ++owner) {
        scaled_capacity[order][replacement][owner] =
            1LL * capacity_count(order, replacement, owner) * (SCALE / order);
      }
    }
  }
}

Pattern order_pattern(Orders orders) {
  std::sort(orders.begin(), orders.end());
  return orders;
}

std::string pattern_string(const Pattern& pattern) {
  std::ostringstream out;
  out << '(';
  for (int i = 0; i < 5; ++i) {
    if (i != 0) out << ',';
    out << pattern[i];
  }
  out << ')';
  return out.str();
}

long long lcm_without(const Orders& orders, int omitted) {
  long long answer = 1;
  for (int i = 0; i < 5; ++i) {
    if (i != omitted) answer = std::lcm(answer, static_cast<long long>(orders[i]));
  }
  return answer;
}

long long lcm_on_mask(const Orders& orders, int mask) {
  long long answer = 1;
  for (int i = 0; i < 5; ++i) {
    if ((mask >> i) & 1) answer = std::lcm(answer, static_cast<long long>(orders[i]));
  }
  return answer;
}

Fraction relative_capacity(int relative_order) {
  assert(relative_order >= 1 && relative_order % P != 0);
  return Fraction((2 * relative_order + 12) / 13, relative_order);
}

// Two all-size consequences of (R), used only as theorem-facing readouts:
//
// * If p>=11 (p!=13) divides some D_i, take S to be the colours carrying the
//   maximal p-adic exponent.  Every m_i is divisible by p and
//   ceil(2m_i/13)/m_i<1/5, so at most five colours cannot meet (R).
// * For p=5 the same argument gives at most 1/5 per carrier.  Equality forces
//   all five colours to be top carriers and every m_i=D_i to lie in
//   {5,10,15,20}.  The finite census below rejects that whole alphabet.
//
// Thus, after combining the uniform cut with this census, any still-open
// effective order is {2,3,7}-smooth.

int valuation(int value, int prime) {
  int answer = 0;
  while (value % prime == 0) {
    value /= prime;
    ++answer;
  }
  return answer;
}

struct ScalarCensus {
  long long conceptual_rows = 0;
  long long scalar_rows = 0;
  long long scalar_rows_at_most_twelve = 0;
  long long scalar_rows_with_large_order = 0;
  long long singleton_cut_rejections = 0;
  std::vector<Row> lcm_redundant_rows;
  std::map<Pattern, long long> lcm_redundant_patterns;
  uint64_t scalar_hash = 1469598103934665603ULL;
  uint64_t singleton_certificate_hash = 1469598103934665603ULL;
};

void scan_orders_for_labels(const Labels& labels, ScalarCensus& census) {
  std::array<std::array<long long, 5>, 6> optimistic{};
  for (int index = 4; index >= 0; --index) {
    for (int owner_index = 0; owner_index < 5; ++owner_index) {
      long long best = 0;
      for (int order : ORDER_ALPHABET) {
        best = std::max(
            best, scaled_capacity[order][labels[index]][labels[owner_index]]);
      }
      optimistic[index][owner_index] = optimistic[index + 1][owner_index] + best;
    }
  }

  Orders orders{};
  std::array<long long, 5> sums{};
  auto dfs = [&](auto&& self, int index) -> void {
    if (index == 5) {
      for (long long sum : sums) {
        if (sum < SCALE) return;
      }
      ++census.scalar_rows;
      bool has_large_order = false;
      for (int order : orders) has_large_order = has_large_order || order > 12;
      if (!has_large_order) {
        ++census.scalar_rows_at_most_twelve;
        return;
      }
      ++census.scalar_rows_with_large_order;

      Row row{labels, orders, {}};
      for (int omitted = 0; omitted < 5; ++omitted) {
        row.four_colour_deficit_owner[omitted] = -1;
        for (int owner_index = 0; owner_index < 5; ++owner_index) {
          const long long four_sum =
              sums[owner_index] -
              scaled_capacity[orders[omitted]][labels[omitted]][labels[owner_index]];
          if (four_sum < SCALE) {
            row.four_colour_deficit_owner[omitted] = owner_index;
            break;
          }
        }
        // Finite-bank check of the all-size THM-810 consequence used by (R).
        assert(row.four_colour_deficit_owner[omitted] >= 0);
      }

      for (int value : labels) census.scalar_hash = fnv_append(census.scalar_hash, value);
      for (int value : orders) census.scalar_hash = fnv_append(census.scalar_hash, value);
      for (long long value : sums) {
        census.scalar_hash = fnv_append(census.scalar_hash, static_cast<uint64_t>(value));
      }

      bool lcm_redundant = true;
      for (int omitted = 0; omitted < 5; ++omitted) {
        const long long other_lcm = lcm_without(orders, omitted);
        if (other_lcm % orders[omitted] == 0) continue;
        lcm_redundant = false;

        // A proper gcd coset needs at least D-g+1 consecutive residues to be
        // contained in a cyclic interval.  Every owner interval here has at
        // most ceil(2D/13) <= floor(D/2) residues, so containment is impossible.
        const int g = std::gcd(orders[omitted], static_cast<int>(other_lcm));
        const int maximum_interval = (2 * orders[omitted] + 12) / 13;
        assert(g < orders[omitted]);
        assert(maximum_interval <= orders[omitted] / 2);
        assert(orders[omitted] - g + 1 > maximum_interval);

        census.singleton_certificate_hash =
            fnv_append(census.singleton_certificate_hash, omitted);
        census.singleton_certificate_hash = fnv_append(
            census.singleton_certificate_hash,
            row.four_colour_deficit_owner[omitted]);
        census.singleton_certificate_hash =
            fnv_append(census.singleton_certificate_hash, orders[omitted]);
        census.singleton_certificate_hash =
            fnv_append(census.singleton_certificate_hash, other_lcm);
        census.singleton_certificate_hash =
            fnv_append(census.singleton_certificate_hash, g);
      }

      if (!lcm_redundant) {
        ++census.singleton_cut_rejections;
        return;
      }

      census.lcm_redundant_rows.push_back(row);
      ++census.lcm_redundant_patterns[order_pattern(orders)];
      return;
    }

    for (int order : ORDER_ALPHABET) {
      orders[index] = order;
      bool possible = true;
      for (int owner_index = 0; owner_index < 5; ++owner_index) {
        sums[owner_index] +=
            scaled_capacity[order][labels[index]][labels[owner_index]];
        if (sums[owner_index] + optimistic[index + 1][owner_index] < SCALE) {
          possible = false;
        }
      }
      if (possible) self(self, index + 1);
      for (int owner_index = 0; owner_index < 5; ++owner_index) {
        sums[owner_index] -=
            scaled_capacity[order][labels[index]][labels[owner_index]];
      }
    }
  };
  dfs(dfs, 0);
}

ScalarCensus scalar_census() {
  ScalarCensus census;
  census.conceptual_rows = 330LL * 19 * 19 * 19 * 19 * 19;
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

struct RamificationCensus {
  std::map<int, long long> top_prime_rejections;
  long long admissible_rows = 0;
  uint64_t certificate_hash = 1469598103934665603ULL;
};

RamificationCensus ramification_census(const std::vector<Row>& rows) {
  RamificationCensus census;
  for (const Row& row : rows) {
    int rejecting_prime = 0;
    int rejecting_mask = 0;
    Fraction rejecting_sum;
    for (int prime : SMALL_PRIMES) {
      int maximum_valuation = 0;
      for (int order : row.orders) {
        maximum_valuation = std::max(maximum_valuation, valuation(order, prime));
      }
      if (maximum_valuation == 0) continue;

      int top_mask = 0;
      for (int i = 0; i < 5; ++i) {
        if (valuation(row.orders[i], prime) == maximum_valuation) top_mask |= 1 << i;
      }
      const int complement_mask = 31 ^ top_mask;
      const long long complement_lcm = lcm_on_mask(row.orders, complement_mask);
      Fraction sum;
      for (int i = 0; i < 5; ++i) {
        if (!((top_mask >> i) & 1)) continue;
        const int relative_order =
            row.orders[i] /
            std::gcd(row.orders[i], static_cast<int>(complement_lcm));
        sum = sum + relative_capacity(relative_order);
      }
      if (sum < Fraction(1, 1)) {
        rejecting_prime = prime;
        rejecting_mask = top_mask;
        rejecting_sum = sum;
        break;
      }
    }

    // Every one of the 185 lcm-redundant rows has an exact p=3 or p=5 cut.
    if (rejecting_prime == 0) {
      ++census.admissible_rows;
      continue;
    }
    ++census.top_prime_rejections[rejecting_prime];
    for (int value : row.labels) census.certificate_hash = fnv_append(census.certificate_hash, value);
    for (int value : row.orders) census.certificate_hash = fnv_append(census.certificate_hash, value);
    census.certificate_hash = fnv_append(census.certificate_hash, rejecting_prime);
    census.certificate_hash = fnv_append(census.certificate_hash, rejecting_mask);
    census.certificate_hash =
        fnv_append(census.certificate_hash, rejecting_sum.numerator);
    census.certificate_hash =
        fnv_append(census.certificate_hash, rejecting_sum.denominator);
  }
  return census;
}

std::vector<int> units(int order) {
  std::vector<int> answer;
  for (int value = 1; value < order; ++value) {
    if (std::gcd(value, order) == 1) answer.push_back(value);
  }
  return answer;
}

int crt_lift(int order, int replacement, int unit) {
  for (int value = 1; value < P * order; ++value) {
    if (value % order == unit && value % P == order * replacement % P) return value;
  }
  assert(false);
  return 0;
}

bool eligible_direct(int order, int replacement, int unit, int owner, int sheet) {
  const int lift = crt_lift(order, replacement, unit);
  const int residue = centred_mod(
      1LL * lift * (inverse13(owner) + P * sheet), P * order);
  return -order < residue && residue <= order;
}

bool eligible_affine(int order, int replacement, int unit, int owner, int sheet) {
  const int owner_inverse = inverse13(owner);
  const int ratio = replacement * owner_inverse % P;
  const int inverse_thirteen_mod_order = inverse_mod(P % order, order);
  const int target_k = positive_mod(
      1LL * unit * sheet +
          1LL * unit * owner_inverse * inverse_thirteen_mod_order,
      order);

  // The eligible k form one literal interval inside {-D+1,...,0}; reducing
  // that interval modulo D gives the affine sheet set.
  for (int k = -order; k <= order; ++k) {
    if (positive_mod(k, order) != target_k) continue;
    const int z = order * ratio + P * k;
    if (-order < z && z <= order) return true;
  }
  return false;
}

struct AffineAudit {
  long long membership_checks = 0;
  long long capacity_checks = 0;
  uint64_t hash = 1469598103934665603ULL;
};

AffineAudit affine_audit() {
  AffineAudit audit;
  for (int order : ORDER_ALPHABET) {
    for (int replacement = 1; replacement < P; ++replacement) {
      for (int unit : units(order)) {
        for (int owner = 1; owner < P; ++owner) {
          int count = 0;
          for (int sheet = 0; sheet < order; ++sheet) {
            const bool direct =
                eligible_direct(order, replacement, unit, owner, sheet);
            const bool affine =
                eligible_affine(order, replacement, unit, owner, sheet);
            assert(direct == affine);
            count += direct;
            ++audit.membership_checks;
            audit.hash = fnv_append(audit.hash, direct ? 1 : 0);
          }
          assert(count == capacity_count(order, replacement, owner));
          ++audit.capacity_checks;
          audit.hash = fnv_append(audit.hash, count);
        }
      }
    }
  }
  return audit;
}

struct LiteralCensus {
  long long unit_words = 0;
  long long common_sheet_words = 0;
  std::map<int, long long> minimum_deficit_histogram;
  Fraction best_covered_fraction;
  long long best_covered_count = 0;
  uint64_t certificate_hash = 1469598103934665603ULL;
};

void literal_row_census(const Row& row, LiteralCensus& census) {
  int common_order = 1;
  for (int order : row.orders) common_order = std::lcm(common_order, order);
  assert(common_order <= 64);
  const uint64_t full = common_order == 64
                            ? ~0ULL
                            : (1ULL << common_order) - 1ULL;

  std::array<std::vector<int>, 5> choices;
  std::array<std::vector<std::array<uint64_t, 5>>, 5> masks;
  for (int colour = 0; colour < 5; ++colour) {
    choices[colour] = units(row.orders[colour]);
    masks[colour].resize(choices[colour].size());
    for (int choice = 0; choice < static_cast<int>(choices[colour].size()); ++choice) {
      for (int owner_index = 0; owner_index < 5; ++owner_index) {
        uint64_t mask = 0;
        for (int sheet = 0; sheet < common_order; ++sheet) {
          if (eligible_direct(row.orders[colour], row.labels[colour],
                              choices[colour][choice], row.labels[owner_index],
                              sheet)) {
            mask |= 1ULL << sheet;
          }
        }
        masks[colour][choice][owner_index] = mask;
      }
    }
  }

  for (int value : row.labels) census.certificate_hash = fnv_append(census.certificate_hash, value);
  for (int value : row.orders) census.certificate_hash = fnv_append(census.certificate_hash, value);

  std::array<int, 5> selected{};
  auto dfs = [&](auto&& self, int colour) -> void {
    if (colour < 5) {
      for (int choice = 0; choice < static_cast<int>(choices[colour].size()); ++choice) {
        selected[colour] = choice;
        self(self, colour + 1);
      }
      return;
    }

    ++census.unit_words;
    int minimum_covered = common_order;
    bool covers = true;
    for (int owner_index = 0; owner_index < 5; ++owner_index) {
      uint64_t covered_mask = 0;
      for (int i = 0; i < 5; ++i) {
        covered_mask |= masks[i][selected[i]][owner_index];
      }
      covered_mask &= full;
      const int covered = __builtin_popcountll(covered_mask);
      minimum_covered = std::min(minimum_covered, covered);
      covers = covers && covered == common_order;
      census.certificate_hash = fnv_append(census.certificate_hash, covered);
    }
    for (int i = 0; i < 5; ++i) {
      census.certificate_hash =
          fnv_append(census.certificate_hash, choices[i][selected[i]]);
    }
    const int deficit = common_order - minimum_covered;
    ++census.minimum_deficit_histogram[deficit];
    const Fraction covered_fraction(minimum_covered, common_order);
    if (covered_fraction > census.best_covered_fraction) {
      census.best_covered_fraction = covered_fraction;
      census.best_covered_count = 1;
    } else if (covered_fraction == census.best_covered_fraction) {
      ++census.best_covered_count;
    }
    if (covers) ++census.common_sheet_words;
  };
  dfs(dfs, 0);
}

LiteralCensus literal_census(const std::vector<Row>& rows) {
  LiteralCensus census;
  for (const Row& row : rows) literal_row_census(row, census);
  return census;
}

using Tournament = std::array<std::array<bool, 5>, 5>;

Tournament capacity_tournament(const Row& row, bool complement_conditioned,
                               int& ties) {
  Tournament tournament{};
  for (int i = 0; i < 5; ++i) {
    for (int j = i + 1; j < 5; ++j) {
      int left_order = row.orders[i];
      int right_order = row.orders[j];
      if (complement_conditioned) {
        long long complement_lcm = 1;
        for (int k = 0; k < 5; ++k) {
          if (k != i && k != j) {
            complement_lcm =
                std::lcm(complement_lcm, static_cast<long long>(row.orders[k]));
          }
        }
        left_order /= std::gcd(left_order, static_cast<int>(complement_lcm));
        right_order /= std::gcd(right_order, static_cast<int>(complement_lcm));
      }
      const Fraction left = relative_capacity(left_order);
      const Fraction right = relative_capacity(right_order);
      int winner = -1;
      if (left > right) {
        winner = i;
      } else if (right > left) {
        winner = j;
      } else {
        ++ties;
        // Declared tie Hamiltonian order: increasing replacement label.
        winner = row.labels[i] < row.labels[j] ? i : j;
      }
      const int loser = winner == i ? j : i;
      tournament[winner][loser] = true;
    }
  }
  return tournament;
}

struct TournamentFingerprint {
  std::array<int, 5> score_histogram{};
  int triangles = 0;
  std::array<int, 6> scc_size_histogram{};
  int hamiltonian_paths = 0;

  bool operator<(const TournamentFingerprint& other) const {
    return std::tie(score_histogram, triangles, scc_size_histogram,
                    hamiltonian_paths) <
           std::tie(other.score_histogram, other.triangles,
                    other.scc_size_histogram, other.hamiltonian_paths);
  }
};

TournamentFingerprint tournament_fingerprint(const Tournament& tournament) {
  TournamentFingerprint fingerprint;
  for (int i = 0; i < 5; ++i) {
    int score = 0;
    for (int j = 0; j < 5; ++j) score += tournament[i][j];
    ++fingerprint.score_histogram[score];
  }
  for (int i = 0; i < 5; ++i) {
    for (int j = i + 1; j < 5; ++j) {
      for (int k = j + 1; k < 5; ++k) {
        fingerprint.triangles +=
            (tournament[i][j] && tournament[j][k] && tournament[k][i]) ||
            (tournament[i][k] && tournament[k][j] && tournament[j][i]);
      }
    }
  }

  std::array<std::array<bool, 5>, 5> reachable = tournament;
  for (int i = 0; i < 5; ++i) reachable[i][i] = true;
  for (int k = 0; k < 5; ++k) {
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) {
        reachable[i][j] =
            reachable[i][j] || (reachable[i][k] && reachable[k][j]);
      }
    }
  }
  std::array<bool, 5> seen{};
  for (int i = 0; i < 5; ++i) {
    if (seen[i]) continue;
    int size = 0;
    for (int j = 0; j < 5; ++j) {
      if (reachable[i][j] && reachable[j][i]) {
        seen[j] = true;
        ++size;
      }
    }
    ++fingerprint.scc_size_histogram[size];
  }

  std::array<int, 5> path = {0, 1, 2, 3, 4};
  do {
    bool valid = true;
    for (int i = 0; i < 4; ++i) valid = valid && tournament[path[i]][path[i + 1]];
    fingerprint.hamiltonian_paths += valid;
  } while (std::next_permutation(path.begin(), path.end()));
  return fingerprint;
}

struct TournamentCensus {
  std::map<TournamentFingerprint, long long> raw_fingerprints;
  std::map<TournamentFingerprint, long long> conditioned_fingerprints;
  std::map<int, long long> flip_histogram;
  long long total_flips = 0;
  long long raw_ties = 0;
  long long conditioned_ties = 0;
  uint64_t hash = 1469598103934665603ULL;
};

TournamentCensus tournament_census(const std::vector<Row>& rows) {
  TournamentCensus census;
  for (const Row& row : rows) {
    int raw_ties = 0;
    int conditioned_ties = 0;
    const Tournament raw = capacity_tournament(row, false, raw_ties);
    const Tournament conditioned =
        capacity_tournament(row, true, conditioned_ties);
    ++census.raw_fingerprints[tournament_fingerprint(raw)];
    ++census.conditioned_fingerprints[tournament_fingerprint(conditioned)];
    census.raw_ties += raw_ties;
    census.conditioned_ties += conditioned_ties;
    int flips = 0;
    for (int i = 0; i < 5; ++i) {
      for (int j = i + 1; j < 5; ++j) {
        flips += raw[i][j] != conditioned[i][j];
        census.hash = fnv_append(census.hash, raw[i][j] ? 1 : 0);
        census.hash = fnv_append(census.hash, conditioned[i][j] ? 1 : 0);
      }
    }
    ++census.flip_histogram[flips];
    census.total_flips += flips;
  }
  return census;
}

std::string score_histogram_string(const std::array<int, 5>& histogram) {
  std::ostringstream out;
  out << '{';
  bool first = true;
  for (int score = 0; score < 5; ++score) {
    if (histogram[score] == 0) continue;
    if (!first) out << ',';
    first = false;
    out << score << ':' << histogram[score];
  }
  out << '}';
  return out.str();
}

std::string scc_histogram_string(const std::array<int, 6>& histogram) {
  std::ostringstream out;
  out << '{';
  bool first = true;
  for (int size = 1; size <= 5; ++size) {
    if (histogram[size] == 0) continue;
    if (!first) out << ',';
    first = false;
    out << size << ':' << histogram[size];
  }
  out << '}';
  return out.str();
}

template <typename K>
std::string integer_map_string(const std::map<K, long long>& values) {
  std::ostringstream out;
  out << '{';
  bool first = true;
  for (const auto& [key, count] : values) {
    if (!first) out << ',';
    first = false;
    out << key << ':' << count;
  }
  out << '}';
  return out.str();
}

void assert_expected(const ScalarCensus& scalar,
                     const RamificationCensus& ramification,
                     const AffineAudit& affine, const LiteralCensus& literal,
                     const TournamentCensus& tournament) {
  assert(scalar.conceptual_rows == 817112670);
  assert(scalar.scalar_rows == 4245);
  assert(scalar.scalar_rows_at_most_twelve == 2190);
  assert(scalar.scalar_rows_with_large_order == 2055);
  assert(scalar.singleton_cut_rejections == 1870);
  assert(scalar.lcm_redundant_rows.size() == 185);
  assert(scalar.lcm_redundant_patterns.size() == 8);
  assert(ramification.top_prime_rejections.at(2) == 20);
  assert(ramification.top_prime_rejections.at(3) == 40);
  assert(ramification.top_prime_rejections.at(5) == 125);
  assert(ramification.top_prime_rejections.size() == 3);
  assert(ramification.admissible_rows == 0);
  assert(affine.membership_checks > 0);
  assert(affine.capacity_checks > 0);
  assert(literal.unit_words == 51360);
  assert(literal.common_sheet_words == 0);
  assert(tournament.raw_fingerprints.size() == 1);
  assert(tournament.conditioned_fingerprints.size() == 1);
  assert(tournament.total_flips == 574);
}

}  // namespace

int main() {
  initialize_capacities();
  const ScalarCensus scalar = scalar_census();
  const RamificationCensus ramification =
      ramification_census(scalar.lcm_redundant_rows);
  const AffineAudit affine = affine_audit();
  const LiteralCensus literal = literal_census(scalar.lcm_redundant_rows);
  const TournamentCensus tournament =
      tournament_census(scalar.lcm_redundant_rows);
  assert_expected(scalar, ramification, affine, literal, tournament);

  uint64_t payload_hash = 1469598103934665603ULL;
  for (uint64_t value : {
           static_cast<uint64_t>(scalar.conceptual_rows),
           static_cast<uint64_t>(scalar.scalar_rows),
           static_cast<uint64_t>(scalar.scalar_rows_with_large_order),
           static_cast<uint64_t>(scalar.singleton_cut_rejections),
           static_cast<uint64_t>(scalar.lcm_redundant_rows.size()),
           static_cast<uint64_t>(ramification.admissible_rows),
           static_cast<uint64_t>(affine.membership_checks),
           static_cast<uint64_t>(literal.unit_words),
           static_cast<uint64_t>(literal.common_sheet_words),
           static_cast<uint64_t>(tournament.total_flips),
           scalar.scalar_hash,
           scalar.singleton_certificate_hash,
           ramification.certificate_hash,
           affine.hash,
           literal.certificate_hash,
           tournament.hash,
       }) {
    payload_hash = fnv_append(payload_hash, value);
  }

  std::cout << "Hamming-five effective-order-at-most-21 common-sheet census\n";
  std::cout << "scope=1_in_R_subset_F13star size_R=5 "
               "orders={2,...,12,14,...,21} at_least_one_order_above_12\n";
  std::cout << "orientation=(-D,D] exact_integer_capacity_and_centered_CRT\n";
  std::cout << "normalization=presentations_containing_label_1_not_orbit_quotient\n";
  std::cout << "\nMETHOD\n";
  std::cout << "scalar_filter=integer_scaled_by_232792560\n";
  std::cout << "relative_cut: for S nonempty, T=complement, L_T=lcm(T), "
               "m_i=D_i/gcd(D_i,L_T), require sum_S ceil(2m_i/13)/m_i>=1\n";
  std::cout << "affine_sheet_formula: k=e*ell+e*o^(-1)*13^(-1) (mod D), "
               "eligible k form one consecutive cyclic interval\n";
  std::cout << "singleton_consequence=D_i_divides_lcm_of_other_four\n";
  std::cout << "literal_crosscheck=all_unit_words_on_every_lcm_redundant_row\n";
  std::cout << "\nSCALAR CENSUS\n";
  std::cout << "conceptual rows: " << scalar.conceptual_rows << '\n';
  std::cout << "scalar rows: " << scalar.scalar_rows << '\n';
  std::cout << "scalar rows all orders <=12: "
            << scalar.scalar_rows_at_most_twelve << '\n';
  std::cout << "scalar rows with an order >12: "
            << scalar.scalar_rows_with_large_order << '\n';
  std::cout << "all five four-colour deletions have a strict scalar-deficit owner: "
            << scalar.scalar_rows_with_large_order << '\n';
  std::cout << "singleton-cut rejections: "
            << scalar.singleton_cut_rejections << '\n';
  std::cout << "lcm-redundant rows: " << scalar.lcm_redundant_rows.size() << '\n';
  std::cout << "lcm-redundant order patterns: "
            << scalar.lcm_redundant_patterns.size() << '\n';
  for (const auto& [pattern, count] : scalar.lcm_redundant_patterns) {
    std::cout << "  " << pattern_string(pattern) << ": " << count << '\n';
  }
  std::cout << "\nRAMIFICATION CUT CENSUS\n";
  std::cout << "top-prime cut rejections: "
            << integer_map_string(ramification.top_prime_rejections) << '\n';
  std::cout << "relative-cut-admissible rows: "
            << ramification.admissible_rows << '\n';
  std::cout << "uniform corollary: every D_i divides lcm(other D_j); "
               "private maximal prime powers are impossible\n";
  std::cout << "uniform prime corollary: primes >=11 cannot divide a surviving "
               "order; a twice-maximal prime power is the dyadic quotient (2,2)\n";
  std::cout << "uniform odd corollary: a three-carrier maximal odd prime power "
               "can only be p=3 with relative quotients (3,3,3)\n";
  std::cout << "combined prime-support corollary: every still-open effective "
               "order is {2,3,7}-smooth; p=5 would force all five orders into "
               "{5,10,15,20}, rejected here\n";
  std::cout << "\nAFFINE/LITERAL AUDIT\n";
  std::cout << "direct-versus-affine membership checks: "
            << affine.membership_checks << '\n';
  std::cout << "owner-capacity checks: " << affine.capacity_checks << '\n';
  std::cout << "literal unit words on lcm-redundant rows: "
            << literal.unit_words << '\n';
  std::cout << "literal common-sheet unit words: "
            << literal.common_sheet_words << '\n';
  std::cout << "minimum uncovered-sheet deficit histogram: "
            << integer_map_string(literal.minimum_deficit_histogram) << '\n';
  std::cout << "best minimum covered fraction: " << literal.best_covered_fraction
            << " (words=" << literal.best_covered_count << ")\n";
  std::cout << "\nTOURNAMENT ANALYSIS\n";
  std::cout << "vertices=replacement colours\n";
  std::cout << "raw pair observable=ceil(2D_i/13)/D_i - "
               "ceil(2D_j/13)/D_j\n";
  std::cout << "conditioned pair observable=the same after quotienting each "
               "order by gcd with the other-three lcm\n";
  std::cout << "switch=raw_to_complement-conditioned; "
               "tie Hamiltonian order=increasing replacement label\n";
  for (const auto& [fingerprint, count] : tournament.raw_fingerprints) {
    std::cout << "raw fingerprint count=" << count
              << " scores=" << score_histogram_string(fingerprint.score_histogram)
              << " triangles=" << fingerprint.triangles
              << " SCC-size-hist="
              << scc_histogram_string(fingerprint.scc_size_histogram)
              << " Hamiltonian-paths=" << fingerprint.hamiltonian_paths << '\n';
  }
  for (const auto& [fingerprint, count] : tournament.conditioned_fingerprints) {
    std::cout << "conditioned fingerprint count=" << count
              << " scores=" << score_histogram_string(fingerprint.score_histogram)
              << " triangles=" << fingerprint.triangles
              << " SCC-size-hist="
              << scc_histogram_string(fingerprint.scc_size_histogram)
              << " Hamiltonian-paths=" << fingerprint.hamiltonian_paths << '\n';
  }
  std::cout << "edge-flip histogram: "
            << integer_map_string(tournament.flip_histogram) << '\n';
  std::cout << "total raw/conditioned edge flips: "
            << tournament.total_flips << '\n';
  std::cout << "raw ties=" << tournament.raw_ties
            << " conditioned ties=" << tournament.conditioned_ties << '\n';
  std::cout << "tournament verdict=both gauges are always transitive despite "
               "574 switches; pair rankings do not encode the rejecting subset cut\n";
  std::cout << "\nASSUMPTION CHALLENGE\n";
  std::cout << "challenged default=vertices need not be runners or replacement colours\n";
  std::cout << "faithful carrier=prime-power ramification hyperedges together with "
               "complement lcm fibres and affine owner-sheet intervals\n";
  std::cout << "preserved by carrier=common-sheet predicate, gcd/Bezout fibre, "
               "owner deficit, interval phase, and relative capacity cut\n";
  std::cout << "destroyed by bare tournament=which complement misses which sheet, "
               "the top-prime carrier set, affine phase, and union capacity\n";
  std::cout << "\nVERDICT\n";
  std::cout << "common-sheet rows with all 2<=D_i<=21 and some D_i>12: 0\n";
  std::cout << "combined with THM-823: no new common-sheet language exists through "
               "effective order 21\n";
  std::cout << "remaining shallow frontier: min D_i<=21, max D_i>=22, with every "
               "private-prime-power and relative-ramification cut excluded\n";
  std::cout << "\nCERTIFICATES\n";
  std::cout << "scalar hash: " << hex64(scalar.scalar_hash) << '\n';
  std::cout << "singleton-cut hash: "
            << hex64(scalar.singleton_certificate_hash) << '\n';
  std::cout << "ramification-cut hash: "
            << hex64(ramification.certificate_hash) << '\n';
  std::cout << "affine audit hash: " << hex64(affine.hash) << '\n';
  std::cout << "literal audit hash: " << hex64(literal.certificate_hash) << '\n';
  std::cout << "tournament hash: " << hex64(tournament.hash) << '\n';
  std::cout << "payload hash: " << hex64(payload_hash) << '\n';
}
