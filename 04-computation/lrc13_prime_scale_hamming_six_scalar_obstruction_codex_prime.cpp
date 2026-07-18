#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Exact primary certificate for the uniform prime-scale scalar obstruction in
// the primitive proper AP-centred Hamming-six bank.  For D=p, unit-independent
// sheet cardinality is the residue count
//
//   a_p(r) = #{ x : -p < x <= p and x = p*r (mod 13) }.
//
// Enlarging p by thirteen adds two complete blocks of thirteen consecutive
// integers, hence a_{p+13}(r)=a_p(r)+2.  The program reconstructs the exact
// base-scale order statistics B6 and B5, checks the two exceptional empty
// scans p=23,29, anchors the p=17 QR/NQR pair, and audits all 139 primes from
// 19 through 839 occurring in the THM-860 finite range.

using Labels = std::array<uint8_t, 6>;

static constexpr int MODULUS = 13;
static constexpr std::array<int, 12> EXPECTED_B6{
    1, 3, 5, 6, 6, 6, 7, 9, 11, 12, 12, 12};
static constexpr std::array<int, 12> EXPECTED_B5{
    1, 3, 5, 5, 5, 5, 6, 8, 10, 10, 10, 10};

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

void fnv_byte(uint64_t &digest, uint8_t value) {
  digest ^= value;
  digest *= 1'099'511'628'211ULL;
}

void fnv_u32(uint64_t &digest, uint32_t value) {
  for (int shift = 0; shift < 32; shift += 8)
    fnv_byte(digest, static_cast<uint8_t>(value >> shift));
}

int positive_mod(int value, int modulus) {
  int result = value % modulus;
  if (result < 0) result += modulus;
  return result;
}

int inverse_mod_13(int value) {
  for (int candidate = 1; candidate < MODULUS; ++candidate)
    if (value * candidate % MODULUS == 1) return candidate;
  fail("nonunit modulo thirteen");
}

int residue_count_literal(int scale, int ratio) {
  require(scale >= 1, "nonpositive prime scale");
  require(1 <= ratio && ratio < MODULUS, "bad nonzero ratio");
  const int target = positive_mod(scale * ratio, MODULUS);
  int count = 0;
  for (int value = -scale + 1; value <= scale; ++value)
    count += positive_mod(value, MODULUS) == target;
  return count;
}

int d1_cardinality_literal(int scale, int ratio) {
  int centered_ratio = ratio;
  if (2 * centered_ratio > MODULUS) centered_ratio -= MODULUS;
  int count = 0;
  for (int sheet = 0; sheet < scale; ++sheet) {
    static_cast<void>(sheet);
    count += -1 < centered_ratio && centered_ratio <= 1;
  }
  return count;
}

std::array<int, 12> residue_row(int scale) {
  std::array<int, 12> result{};
  for (int ratio = 1; ratio < MODULUS; ++ratio)
    result[ratio - 1] = residue_count_literal(scale, ratio);
  return result;
}

void maximize_subset(const std::array<int, 12> &weights, int next_ratio,
                     int remaining, int sum, int &maximum,
                     int &maximizer_count) {
  if (remaining == 0) {
    if (sum > maximum) {
      maximum = sum;
      maximizer_count = 1;
    } else if (sum == maximum) {
      ++maximizer_count;
    }
    return;
  }
  if (13 - next_ratio < remaining) return;
  for (int ratio = next_ratio; ratio < MODULUS; ++ratio)
    maximize_subset(weights, ratio + 1, remaining - 1,
                    sum + weights[ratio - 1], maximum, maximizer_count);
}

std::pair<int, int> b6_bound(const std::array<int, 12> &weights) {
  // Owner-normalized support ratios are a six-subset containing ratio one.
  int maximum = -1;
  int maximizers = 0;
  maximize_subset(weights, 2, 5, weights[0], maximum, maximizers);
  return {maximum, maximizers};
}

std::pair<int, int> b5_bound(const std::array<int, 12> &weights) {
  // With a D=1 coordinate, choose a different (D=p) owner.  At most five
  // D=p providers remain, bounded by an arbitrary five-subset of ratios.
  int maximum = -1;
  int maximizers = 0;
  maximize_subset(weights, 1, 5, 0, maximum, maximizers);
  return {maximum, maximizers};
}

bool is_prime(int value) {
  if (value < 2) return false;
  for (int divisor = 2; divisor * divisor <= value; ++divisor)
    if (value % divisor == 0) return false;
  return true;
}

int all_dp_capacity(int prime, const Labels &support, int owner) {
  const int inverse = inverse_mod_13(owner);
  int capacity = 0;
  for (int provider : support) {
    const int ratio = provider * inverse % MODULUS;
    capacity += residue_count_literal(prime, ratio);
  }
  return capacity;
}

std::set<Labels> scan_all_dp_supports(int prime, uint64_t &digest,
                                     uint64_t &supports_scanned) {
  const std::array<int, 12> counts = residue_row(prime);
  std::set<Labels> survivors;
  fnv_u32(digest, static_cast<uint32_t>(prime));
  for (int a = 1; a <= 7; ++a)
    for (int b = a + 1; b <= 8; ++b)
      for (int c = b + 1; c <= 9; ++c)
        for (int d = c + 1; d <= 10; ++d)
          for (int e = d + 1; e <= 11; ++e)
            for (int f = e + 1; f <= 12; ++f) {
              const Labels labels{static_cast<uint8_t>(a),
                                  static_cast<uint8_t>(b),
                                  static_cast<uint8_t>(c),
                                  static_cast<uint8_t>(d),
                                  static_cast<uint8_t>(e),
                                  static_cast<uint8_t>(f)};
              ++supports_scanned;
              bool scalar = true;
              for (int owner : labels) {
                const int inverse = inverse_mod_13(owner);
                int capacity = 0;
                for (int provider : labels)
                  capacity += counts[provider * inverse % MODULUS - 1];
                if (capacity < prime) {
                  scalar = false;
                  break;
                }
              }
              for (uint8_t label : labels) fnv_byte(digest, label);
              fnv_byte(digest, static_cast<uint8_t>(scalar));
              if (scalar) survivors.insert(labels);
            }
  return survivors;
}

Labels multiply_support(const Labels &labels, int multiplier) {
  Labels result{};
  for (int coordinate = 0; coordinate < 6; ++coordinate)
    result[coordinate] = static_cast<uint8_t>(
        multiplier * labels[coordinate] % MODULUS);
  std::sort(result.begin(), result.end());
  return result;
}

struct RatioTournament {
  int ties = 0;
  int flips = 0;
  int directed_triangles = 0;
  uint64_t hamiltonian_paths = 0;
  std::array<uint8_t, 12> scores{};
};

RatioTournament ratio_tournament(const std::array<int, 12> &weights) {
  // Pair observable is the exact ordered cardinality pair (a_s(r),a_s(u)).
  // Orient toward the larger value; ties follow 1->2->...->12.
  std::array<uint16_t, 12> out{};
  RatioTournament result;
  for (int left = 0; left < 12; ++left)
    for (int right = left + 1; right < 12; ++right) {
      int winner = left;
      if (weights[left] == weights[right]) {
        ++result.ties;
      } else if (weights[right] > weights[left]) {
        winner = right;
        ++result.flips;
      }
      const int loser = left + right - winner;
      out[winner] |= static_cast<uint16_t>(1U << loser);
      ++result.scores[winner];
    }
  for (int i = 0; i < 12; ++i)
    for (int j = i + 1; j < 12; ++j)
      for (int k = j + 1; k < 12; ++k) {
        const bool forward = ((out[i] >> j) & 1U) &&
                             ((out[j] >> k) & 1U) &&
                             ((out[k] >> i) & 1U);
        const bool reverse = ((out[i] >> k) & 1U) &&
                             ((out[k] >> j) & 1U) &&
                             ((out[j] >> i) & 1U);
        result.directed_triangles += forward || reverse;
      }
  std::array<std::array<uint64_t, 12>, 1U << 12> paths{};
  for (int last = 0; last < 12; ++last) paths[1U << last][last] = 1;
  for (int mask = 1; mask < (1 << 12); ++mask)
    for (int last = 0; last < 12; ++last)
      if ((mask >> last) & 1) {
        const int previous_mask = mask ^ (1 << last);
        for (int previous = 0; previous < 12; ++previous)
          if (((previous_mask >> previous) & 1) &&
              ((out[previous] >> last) & 1U))
            paths[mask][last] += paths[previous_mask][previous];
      }
  for (int last = 0; last < 12; ++last)
    result.hamiltonian_paths += paths.back()[last];
  return result;
}

std::string array_string(const std::array<int, 12> &values) {
  std::ostringstream out;
  for (int index = 0; index < 12; ++index) {
    if (index) out << ',';
    out << values[index];
  }
  return out.str();
}

std::string support_string(const Labels &labels) {
  std::ostringstream out;
  for (int index = 0; index < 6; ++index) {
    if (index) out << ',';
    out << static_cast<int>(labels[index]);
  }
  return out.str();
}

int main() {
  // Every block of thirteen consecutive integers contains every residue once;
  // checking one full period of starting residues is the finite kernel of the
  // uniform block lemma used by the recurrence proof.
  for (int first = 0; first < MODULUS; ++first)
    for (int residue = 0; residue < MODULUS; ++residue) {
      int count = 0;
      for (int offset = 0; offset < MODULUS; ++offset)
        count += positive_mod(first + offset, MODULUS) == residue;
      require(count == 1, "complete residue-block lemma failed");
    }

  // Literal recurrence audit.  The proof is uniform: the difference between
  // (-p-13,p+13] and (-p,p] is two disjoint length-thirteen blocks, each of
  // which contains exactly one representative of every residue modulo 13.
  for (int scale = 1; scale <= 839; ++scale)
    for (int ratio = 1; ratio < MODULUS; ++ratio)
      require(residue_count_literal(scale + 13, ratio) ==
                  residue_count_literal(scale, ratio) + 2,
              "a_(p+13)=a_p+2 recurrence mismatch");

  std::array<int, 12> b6{};
  std::array<int, 12> b5{};
  std::array<int, 12> b6_maximizers{};
  std::array<int, 12> b5_maximizers{};
  std::array<int, 12> tournament_ties{};
  std::array<int, 12> tournament_flips{};
  uint64_t base_digest = 14'695'981'039'346'656'037ULL;
  uint64_t tournament_digest = 14'695'981'039'346'656'037ULL;
  for (int residue = 1; residue < MODULUS; ++residue) {
    const std::array<int, 12> weights = residue_row(residue);
    const auto [six_bound, six_maximizers] = b6_bound(weights);
    const auto [five_bound, five_maximizers] = b5_bound(weights);
    b6[residue - 1] = six_bound;
    b5[residue - 1] = five_bound;
    b6_maximizers[residue - 1] = six_maximizers;
    b5_maximizers[residue - 1] = five_maximizers;
    for (int value : weights) fnv_u32(base_digest, value);
    fnv_u32(base_digest, six_bound);
    fnv_u32(base_digest, five_bound);

    const RatioTournament tournament = ratio_tournament(weights);
    require(tournament.directed_triangles == 0,
            "ratio tournament has a directed triangle");
    require(tournament.hamiltonian_paths == 1,
            "ratio tournament path count mismatch");
    std::array<uint8_t, 12> sorted_scores = tournament.scores;
    std::sort(sorted_scores.begin(), sorted_scores.end());
    for (int score = 0; score < 12; ++score)
      require(sorted_scores[score] == score,
              "ratio tournament score histogram mismatch");
    tournament_ties[residue - 1] = tournament.ties;
    tournament_flips[residue - 1] = tournament.flips;
    fnv_u32(tournament_digest, tournament.ties);
    fnv_u32(tournament_digest, tournament.flips);
  }
  require(b6 == EXPECTED_B6, "B6 vector mismatch");
  require(b5 == EXPECTED_B5, "B5 vector mismatch");
  require(b6_maximizers ==
              std::array<int, 12>{462, 84, 7, 6, 56, 252,
                                  462, 84, 7, 6, 56, 252},
          "B6 maximizer-count vector mismatch");
  require(b5_maximizers ==
              std::array<int, 12>{330, 36, 1, 21, 126, 462,
                                  330, 36, 1, 21, 126, 462},
          "B5 maximizer-count vector mismatch");
  require(base_digest == 0x8a078ddec8902cc1ULL,
          "base cardinality/bound digest mismatch");
  for (int residue = 1; residue < MODULUS; ++residue) {
    require(b6[residue - 1] <= residue + 2,
            "B6 exceeds uniform residue bound");
    require(b5[residue - 1] <= residue + 2,
            "B5 exceeds uniform residue bound");
  }

  // D=1 cardinality is p at ratio one and zero at every other nonzero ratio.
  // If any D=1 provider occurs, choose a distinct D=p owner (heredity supplies
  // at least two).  Its scalar capacity is at most 10q+B5(s), strictly below
  // p=13q+s because B5(s)<=s+2 and q>=1.
  for (int quotient = 1; quotient <= 64; ++quotient)
    for (int residue = 1; residue < MODULUS; ++residue) {
      const int scale = 13 * quotient + residue;
      require(10 * quotient + b5[residue - 1] < scale,
              "D=1 exclusion inequality failed");
      for (int ratio = 1; ratio < MODULUS; ++ratio) {
        const int d1_cardinality = d1_cardinality_literal(scale, ratio);
        require(d1_cardinality == (ratio == 1 ? scale : 0),
                "D=1 cardinality mismatch");
      }
    }

  // With all providers at D=p, capacity is at most 12q+B6(s).  For q>=3,
  // B6(s)<=s+2 makes the deficit at least q-2>=1.
  for (int quotient = 3; quotient <= 64; ++quotient)
    for (int residue = 1; residue < MODULUS; ++residue)
      require(12 * quotient + b6[residue - 1] <
                  13 * quotient + residue,
              "uniform all-Dp inequality failed");

  std::vector<int> audit_primes;
  for (int value = 19; value <= 839; ++value)
    if (is_prime(value)) audit_primes.push_back(value);
  require(audit_primes.size() == 139 && audit_primes.front() == 19 &&
              audit_primes.back() == 839,
          "THM-860 prime audit list mismatch");

  const Labels qr{1, 3, 4, 9, 10, 12};
  const Labels nqr{2, 5, 6, 7, 8, 11};
  const std::set<Labels> expected_anchor{qr, nqr};
  uint64_t support_digest = 14'695'981'039'346'656'037ULL;
  uint64_t scanned_supports = 0;
  const std::set<Labels> anchor =
      scan_all_dp_supports(17, support_digest, scanned_supports);
  require(anchor == expected_anchor, "p=17 QR/NQR anchor mismatch");

  std::set<Labels> anchor_orbit;
  for (int multiplier = 1; multiplier < MODULUS; ++multiplier)
    anchor_orbit.insert(multiply_support(qr, multiplier));
  require(anchor_orbit == expected_anchor,
          "p=17 multiplication orbit mismatch");

  int uniform_q3 = 0;
  int low_bound = 0;
  int exceptional_scans = 0;
  uint64_t prime_digest = 14'695'981'039'346'656'037ULL;
  for (int prime : audit_primes) {
    const int quotient = prime / MODULUS;
    const int residue = prime % MODULUS;
    require(residue != 0, "thirteen entered prime audit unexpectedly");
    fnv_u32(prime_digest, static_cast<uint32_t>(prime));
    fnv_u32(prime_digest, static_cast<uint32_t>(quotient));
    fnv_u32(prime_digest, static_cast<uint32_t>(residue));
    require(10 * quotient + b5[residue - 1] < prime,
            "prime D=1 exclusion failed");
    if (quotient >= 3) {
      ++uniform_q3;
      require(12 * quotient + b6[residue - 1] < prime,
              "q>=3 prime scalar bound failed");
    } else if (12 * quotient + b6[residue - 1] < prime) {
      ++low_bound;
    } else {
      ++exceptional_scans;
      require(prime == 23 || prime == 29,
              "unexpected low-prime exceptional case");
    }
    const std::set<Labels> survivors =
        scan_all_dp_supports(prime, support_digest, scanned_supports);
    require(survivors.empty(), "prime audit has a scalar support survivor");
  }
  require(uniform_q3 == 134 && low_bound == 3 && exceptional_scans == 2,
          "prime proof-route census mismatch");
  require(scanned_supports == 140ULL * 924ULL,
          "support scan census mismatch");
  require(prime_digest == 0x1e008620107180afULL,
          "prime-list digest mismatch");
  require(support_digest == 0x339ca4032c7ab2ddULL,
          "support-scan digest mismatch");
  require(tournament_ties ==
              std::array<int, 12>{55, 39, 31, 31, 39, 55,
                                  55, 39, 31, 31, 39, 55},
          "ratio tournament tie vector mismatch");
  require(tournament_flips ==
              std::array<int, 12>{0, 8, 12, 12, 8, 0,
                                  0, 8, 12, 12, 8, 0},
          "ratio tournament flip vector mismatch");
  require(tournament_digest == 0x9191074f3eae0065ULL,
          "ratio tournament digest mismatch");

  std::cout << "uniform prime-scale AP-centred Hamming-six scalar obstruction\n";
  std::cout << "literal cardinality a_p(r)=#{-p<x<=p:x congruent p*r mod13}\n";
  std::cout << "recurrence a_(p+13)(r)=a_p(r)+2: two added length-thirteen residue blocks; literal audit p=1..839, r=1..12\n";
  std::cout << "B6 " << array_string(b6) << '\n';
  std::cout << "B5 " << array_string(b5) << '\n';
  std::cout << "B6 maximizer counts " << array_string(b6_maximizers) << '\n';
  std::cout << "B5 maximizer counts " << array_string(b5_maximizers) << '\n';
  std::cout << "base cardinality/bound FNV64 " << std::hex << base_digest
            << std::dec << '\n';
  std::cout << "D1 exclusion 10q+B5(s)<13q+s for q>=1; all-Dp obstruction 12q+B6(s)<13q+s for q>=3\n";
  std::cout << "p17 anchor QR [" << support_string(qr) << "] NQR ["
            << support_string(nqr) << "]; one multiplication orbit of size "
            << anchor_orbit.size() << '\n';
  std::cout << "exceptional exact scans p23 924 supports empty; p29 924 supports empty\n";
  std::cout << "THM860 audit primes 19..839 count " << audit_primes.size()
            << "; routes q>=3:" << uniform_q3 << " low-bound:" << low_bound
            << " exceptional-scan:" << exceptional_scans << '\n';
  std::cout << "full all-Dp scans primes19..839 plus p17 supports "
            << scanned_supports << "; survivors beyond p17 0\n";
  std::cout << "prime-list FNV64 " << std::hex << prime_digest
            << "; support-scan FNV64 " << support_digest << std::dec << '\n';
  std::cout << "ratio tournament pair observable exact (a_s(r),a_s(u)); orient larger, ties along 1->2->...->12\n";
  std::cout << "ratio tournament ties by s " << array_string(tournament_ties)
            << "; flips by s " << array_string(tournament_flips) << '\n';
  std::cout << "ratio tournament fingerprints all twelve transitive: scores 0..11, cycles 0, SCCs 12, Hamiltonian paths 1; FNV64 "
            << std::hex << tournament_digest << std::dec << '\n';
  std::cout << "challenged vertices ratio weights preserve the uniform capacity bounds, support/owner incidence is restored for exceptional scans, and the oriented tournament alone loses weight magnitudes and the mandatory ratio-one constraint; sheet vertices arise only after scalar survival\n";
}
