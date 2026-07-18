#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
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

// Independent exact referee for the primitive proper AP-centred common-scale-
// twenty-five Hamming-six face.  This implementation does not import the
// Python primary's algebraic CRT, scalar classifier, or immutable-set DP.  It
// finds CRT representatives by literal search, enumerates the leave-one-out
// grammar, discovers the scalar survivors from local cardinalities, and uses
// sorted uint32_t vectors for the full owner-local union banks.  Forward and
// reverse provider traversals must produce the same bank byte for byte.

namespace {

constexpr int P = 13;
constexpr int C = 25;
constexpr uint32_t FULL = (1U << C) - 1U;
constexpr std::array<uint8_t, 3> ORDERS{1, 5, 25};
constexpr std::array<uint8_t, 1> U1{0};
constexpr std::array<uint8_t, 4> U5{1, 2, 3, 4};
constexpr std::array<uint8_t, 20> U25{
    1, 2, 3, 4, 6, 7, 8, 9, 11, 12,
    13, 14, 16, 17, 18, 19, 21, 22, 23, 24};

using Support = std::array<uint8_t, 6>;
using Word = std::array<uint8_t, 6>;
using Capacities = std::array<uint8_t, 6>;
using MaskFibre = std::array<uint32_t, 20>;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

// Compact SHA-256 used only to compare canonical byte streams with the
// independently frozen primary.  No C++ object layout enters a digest.
class Sha256 {
 public:
  void byte(uint8_t value) {
    buffer_[used_++] = value;
    bit_count_ += 8;
    if (used_ == 64) {
      compress(buffer_.data());
      used_ = 0;
    }
  }

  void u16le(uint16_t value) {
    byte(static_cast<uint8_t>(value));
    byte(static_cast<uint8_t>(value >> 8));
  }

  void u32le(uint32_t value) {
    for (int shift = 0; shift < 32; shift += 8)
      byte(static_cast<uint8_t>(value >> shift));
  }

  void u64le(uint64_t value) {
    for (int shift = 0; shift < 64; shift += 8)
      byte(static_cast<uint8_t>(value >> shift));
  }

  std::string finish() {
    const uint64_t original_bits = bit_count_;
    byte(0x80);
    while (used_ != 56) byte(0);
    for (int shift = 56; shift >= 0; shift -= 8)
      byte(static_cast<uint8_t>(original_bits >> shift));
    require(used_ == 0, "SHA-256 final block did not close");
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (uint32_t word : state_)
      for (int shift = 24; shift >= 0; shift -= 8)
        out << std::setw(2) << static_cast<int>(
            static_cast<uint8_t>(word >> shift));
    return out.str();
  }

 private:
  static constexpr std::array<uint32_t, 64> K{
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U,
      0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
      0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
      0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
      0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
      0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
      0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
      0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
      0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
      0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
      0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U,
      0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
      0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U,
      0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
      0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
      0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};

  static uint32_t rotate(uint32_t value, int count) {
    return (value >> count) | (value << (32 - count));
  }

  static uint32_t load_be(const uint8_t *p) {
    return (static_cast<uint32_t>(p[0]) << 24) |
           (static_cast<uint32_t>(p[1]) << 16) |
           (static_cast<uint32_t>(p[2]) << 8) |
           static_cast<uint32_t>(p[3]);
  }

  void compress(const uint8_t *block) {
    std::array<uint32_t, 64> words{};
    for (int i = 0; i < 16; ++i) words[i] = load_be(block + 4 * i);
    for (int i = 16; i < 64; ++i) {
      const uint32_t s0 = rotate(words[i - 15], 7) ^
                          rotate(words[i - 15], 18) ^
                          (words[i - 15] >> 3);
      const uint32_t s1 = rotate(words[i - 2], 17) ^
                          rotate(words[i - 2], 19) ^
                          (words[i - 2] >> 10);
      words[i] = words[i - 16] + s0 + words[i - 7] + s1;
    }
    uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
    uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
    for (int i = 0; i < 64; ++i) {
      const uint32_t s1 = rotate(e, 6) ^ rotate(e, 11) ^ rotate(e, 25);
      const uint32_t choose = (e & f) ^ (~e & g);
      const uint32_t t1 = h + s1 + choose + K[i] + words[i];
      const uint32_t s0 = rotate(a, 2) ^ rotate(a, 13) ^ rotate(a, 22);
      const uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const uint32_t t2 = s0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }
    state_[0] += a;
    state_[1] += b;
    state_[2] += c;
    state_[3] += d;
    state_[4] += e;
    state_[5] += f;
    state_[6] += g;
    state_[7] += h;
  }

  std::array<uint32_t, 8> state_{0x6a09e667U, 0xbb67ae85U,
                                 0x3c6ef372U, 0xa54ff53aU,
                                 0x510e527fU, 0x9b05688cU,
                                 0x1f83d9abU, 0x5be0cd19U};
  std::array<uint8_t, 64> buffer_{};
  uint64_t bit_count_ = 0;
  std::size_t used_ = 0;
};

int order_index(int order) {
  if (order == 1) return 0;
  if (order == 5) return 1;
  if (order == 25) return 2;
  fail("unknown effective order");
}

int unit_count(int order) {
  if (order == 1) return 1;
  if (order == 5) return 4;
  if (order == 25) return 20;
  fail("unknown unit fibre");
}

int unit_at(int order, int index) {
  if (order == 1) return U1[index];
  if (order == 5) return U5[index];
  if (order == 25) return U25[index];
  fail("unknown unit fibre");
}

int inverse13(int value) {
  for (int candidate = 1; candidate < P; ++candidate)
    if (value * candidate % P == 1) return candidate;
  fail("nonunit modulo thirteen");
}

int centered(int value, int modulus) {
  int residue = value % modulus;
  if (residue < 0) residue += modulus;
  return 2 * residue > modulus ? residue - modulus : residue;
}

int literal_crt_base(int label, int order, int unit) {
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P && value % order == unit % order)
      return value;
  fail("literal CRT search found no representative");
}

uint32_t literal_mask(int label, int order, int unit, int owner) {
  const int base = literal_crt_base(label, order, unit);
  const int inverse = inverse13(owner);
  uint32_t mask = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    if (-order < value && value <= order) mask |= 1U << sheet;
  }
  return mask;
}

struct Tables {
  // [label][order-index][owner][unit-index]
  std::array<std::array<std::array<MaskFibre, P>, 3>, P> masks{};
  std::array<std::array<std::array<uint8_t, P>, 3>, P> cards{};
  std::string mask_digest;
  std::string base_digest;
};

Tables build_tables() {
  Tables table;
  Sha256 mask_hash;
  Sha256 base_hash;
  for (int label = 1; label < P; ++label)
    for (int oi = 0; oi < 3; ++oi) {
      const int order = ORDERS[oi];
      for (int ui = 0; ui < unit_count(order); ++ui) {
        const int unit = unit_at(order, ui);
        const int base = literal_crt_base(label, order, unit);
        base_hash.byte(label);
        base_hash.byte(order);
        base_hash.byte(unit);
        base_hash.u16le(base);
        for (int owner = 1; owner < P; ++owner) {
          const uint32_t mask = literal_mask(label, order, unit, owner);
          table.masks[label][oi][owner][ui] = mask;
          mask_hash.byte(label);
          mask_hash.byte(order);
          mask_hash.byte(unit);
          mask_hash.byte(owner);
          mask_hash.u32le(mask);
          const int card = std::popcount(mask);
          if (ui == 0)
            table.cards[label][oi][owner] = static_cast<uint8_t>(card);
          else
            require(table.cards[label][oi][owner] == card,
                    "local cardinality depends on unit");
        }
      }
    }
  table.mask_digest = mask_hash.finish();
  table.base_digest = base_hash.finish();
  return table;
}

bool hereditary_lcm(const Word &word) {
  for (int omitted = 0; omitted < 6; ++omitted) {
    int value = 1;
    for (int i = 0; i < 6; ++i)
      if (i != omitted) value = std::lcm(value, static_cast<int>(word[i]));
    if (value != C) return false;
  }
  return true;
}

void enumerate_words(std::vector<Word> &words, Word &word, int coordinate,
                     uint64_t &state_words, Sha256 &grammar_hash) {
  if (coordinate == 6) {
    const bool hereditary = hereditary_lcm(word);
    const bool prime_square =
        std::count(word.begin(), word.end(), static_cast<uint8_t>(25)) >= 2;
    require(hereditary == prime_square,
            "leave-one-out and prime-square grammars disagree");
    if (!hereditary) return;
    words.push_back(word);
    uint64_t fibre = 1;
    for (uint8_t order : word) {
      grammar_hash.byte(order);
      fibre *= unit_count(order);
    }
    grammar_hash.u64le(fibre);
    state_words += fibre;
    return;
  }
  for (uint8_t order : ORDERS) {
    word[coordinate] = order;
    enumerate_words(words, word, coordinate + 1, state_words, grammar_hash);
  }
}

std::vector<Support> all_supports() {
  std::vector<Support> result;
  for (int a = 1; a <= 7; ++a)
    for (int b = a + 1; b <= 8; ++b)
      for (int c = b + 1; c <= 9; ++c)
        for (int d = c + 1; d <= 10; ++d)
          for (int e = d + 1; e <= 11; ++e)
            for (int f = e + 1; f <= 12; ++f)
              result.push_back(Support{static_cast<uint8_t>(a),
                                       static_cast<uint8_t>(b),
                                       static_cast<uint8_t>(c),
                                       static_cast<uint8_t>(d),
                                       static_cast<uint8_t>(e),
                                       static_cast<uint8_t>(f)});
  return result;
}

struct Survivor {
  Support support;
  Word word;
  Capacities capacities;
};

std::string format(const std::map<uint64_t, uint64_t> &histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, count] : histogram) {
    if (!first) out << ' ';
    first = false;
    out << key << ':' << count;
  }
  return out.str();
}

std::vector<Survivor> scalar_census(
    const std::vector<Support> &supports, const std::vector<Word> &words,
    const Tables &table, std::map<uint64_t, uint64_t> &feasible_histogram,
    std::map<uint64_t, uint64_t> &support_histogram,
    std::string &scalar_digest) {
  std::vector<Survivor> survivors;
  Sha256 scalar_hash;
  for (const Support &support : supports) {
    int on_support = 0;
    for (const Word &word : words) {
      Capacities capacities{};
      int feasible = 0;
      for (int owner_index = 0; owner_index < 6; ++owner_index) {
        const int owner = support[owner_index];
        int capacity = 0;
        for (int i = 0; i < 6; ++i)
          capacity += table.cards[support[i]][order_index(word[i])][owner];
        capacities[owner_index] = static_cast<uint8_t>(capacity);
        feasible += capacity >= C;
      }
      ++feasible_histogram[feasible];
      if (feasible != 6) continue;
      survivors.push_back(Survivor{support, word, capacities});
      ++on_support;
      for (uint8_t value : support) scalar_hash.byte(value);
      for (uint8_t value : word) scalar_hash.byte(value);
      for (uint8_t value : capacities) scalar_hash.byte(value);
    }
    ++support_histogram[on_support];
  }
  scalar_digest = scalar_hash.finish();
  return survivors;
}

std::set<int> quadratic_class() {
  std::set<int> result;
  for (int value = 1; value < P; ++value) result.insert(value * value % P);
  return result;
}

void structural_scalar_audit(const Tables &table,
                             const std::vector<Survivor> &survivors) {
  std::set<int> forbidden;
  for (int ratio = 1; ratio < P; ++ratio)
    if (table.cards[ratio][order_index(5)][1] == 0) forbidden.insert(ratio);
  require(forbidden == std::set<int>({4, 9, 12}),
          "literal order-five forbidden relation");
  std::set<int> inverse_forbidden;
  for (int value : forbidden) inverse_forbidden.insert(inverse13(value));
  require(inverse_forbidden == std::set<int>({3, 10, 12}),
          "inverse forbidden relation");
  const std::set<int> quadratic = quadratic_class();
  std::set<std::pair<int, int>> edges;
  int arcs = 0;
  int reciprocal = 0;
  for (int owner = 1; owner < P; ++owner)
    for (int provider = 1; provider < P; ++provider) {
      if (owner == provider) continue;
      const int ratio = provider * inverse13(owner) % P;
      if (forbidden.contains(ratio)) ++arcs;
      if (forbidden.contains(ratio) || inverse_forbidden.contains(ratio))
        edges.insert(std::minmax(owner, provider));
      if (owner < provider && forbidden.contains(ratio) &&
          forbidden.contains(inverse13(ratio)))
        ++reciprocal;
    }
  require(arcs == 36 && edges.size() == 30 && reciprocal == 6,
          "Cayley K6+K6 fingerprint");
  for (const auto &[left, right] : edges)
    require(quadratic.contains(left) == quadratic.contains(right),
            "Cayley edge crosses square classes");

  for (const Survivor &row : survivors) {
    std::vector<int> five;
    for (int i = 0; i < 6; ++i)
      if (row.word[i] == 5) five.push_back(row.support[i]);
    require(five.size() == 2 &&
                std::count(row.word.begin(), row.word.end(), 25) == 4,
            "survivor order multiplicity");
    require(quadratic.contains(five[0]) != quadratic.contains(five[1]),
            "order-five providers are not in opposite square classes");
    std::set<int> predicted;
    for (int value = 1; value < P; ++value) predicted.insert(value);
    for (int provider : five)
      for (int multiplier : inverse_forbidden)
        predicted.erase(multiplier * provider % P);
    require(std::set<int>(row.support.begin(), row.support.end()) == predicted,
            "survivor is not the forbidden-complement support");
  }
}

std::vector<uint32_t> options(const Tables &table, int label, int order,
                              int owner) {
  std::vector<uint32_t> result;
  for (int ui = 0; ui < unit_count(order); ++ui)
    result.push_back(table.masks[label][order_index(order)][owner][ui]);
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

std::vector<uint32_t> advance(const std::vector<uint32_t> &bank,
                              const std::vector<uint32_t> &choices) {
  std::vector<uint32_t> result;
  result.reserve(bank.size() * choices.size());
  for (uint32_t prior : bank)
    for (uint32_t choice : choices) result.push_back(prior | choice);
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

std::vector<uint32_t> owner_bank(const Survivor &row, int owner,
                                 const Tables &table, bool reverse) {
  std::vector<uint32_t> bank{0};
  for (int step = 0; step < 6; ++step) {
    const int index = reverse ? 5 - step : step;
    bank = advance(bank, options(table, row.support[index], row.word[index], owner));
  }
  return bank;
}

int structural_owner_bound(const Survivor &row, int owner,
                           const Tables &table) {
  std::vector<int> five_indices;
  std::vector<int> twenty_five_indices;
  for (int i = 0; i < 6; ++i)
    (row.word[i] == 5 ? five_indices : twenty_five_indices).push_back(i);
  require(five_indices.size() == 2 && twenty_five_indices.size() == 4,
          "owner incidence grammar");
  int worst_upper = 0;
  for (int left_unit = 0; left_unit < 4; ++left_unit)
    for (int right_unit = 0; right_unit < 4; ++right_unit) {
      uint32_t thick =
          table.masks[row.support[five_indices[0]]][order_index(5)][owner]
                     [left_unit] |
          table.masks[row.support[five_indices[1]]][order_index(5)][owner]
                     [right_unit];
      int forced_intersections = 0;
      int needle_mass = 0;
      for (int index : twenty_five_indices) {
        int minimum = C;
        for (int ui = 0; ui < 20; ++ui) {
          const uint32_t needle =
              table.masks[row.support[index]][order_index(25)][owner][ui];
          minimum = std::min(minimum, std::popcount(needle & thick));
        }
        forced_intersections += minimum;
        needle_mass += table.cards[row.support[index]][order_index(25)][owner];
      }
      // With two distinct thick fibres every one of the four needles is
      // forced to meet them.  If the D5 fibres coincide, their starting mass
      // is only five and the unconditional upper bound is already 20.
      require(std::popcount(thick) == 5 || forced_intersections >= 4,
              "distinct five-fibres did not force four incidences");
      worst_upper = std::max(
          worst_upper,
          std::popcount(thick) + needle_mass - forced_intersections);
    }
  const auto owner_it = std::find(row.support.begin(), row.support.end(), owner);
  require(owner_it != row.support.end(), "owner outside support");
  const int owner_order = row.word[owner_it - row.support.begin()];
  require(worst_upper <= (owner_order == 5 ? 22 : 21),
          "literal five-fibre structural bound");
  return worst_upper;
}

struct TournamentResult {
  int ties;
  int flips;
  std::array<int, 6> scores;
  int triangles;
  int hamiltonian_paths;
};

TournamentResult tournament(const std::array<bool, 6> &feasible,
                            const std::array<int, 6> &maxima,
                            const Capacities &capacities) {
  bool edge[6][6]{};
  int ties = 0;
  int flips = 0;
  for (int left = 0; left < 6; ++left)
    for (int right = left + 1; right < 6; ++right) {
      const auto left_key =
          std::tuple(feasible[left], maxima[left], capacities[left]);
      const auto right_key =
          std::tuple(feasible[right], maxima[right], capacities[right]);
      int winner = left, loser = right;
      if (left_key == right_key) {
        ++ties;
      } else if (right_key > left_key) {
        winner = right;
        loser = left;
        ++flips;
      }
      edge[winner][loser] = true;
    }
  std::array<int, 6> scores{};
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) scores[i] += edge[i][j];
  std::sort(scores.begin(), scores.end());
  int triangles = 0;
  for (int a = 0; a < 6; ++a)
    for (int b = a + 1; b < 6; ++b)
      for (int c = b + 1; c < 6; ++c)
        triangles += (edge[a][b] && edge[b][c] && edge[c][a]) ||
                     (edge[a][c] && edge[c][b] && edge[b][a]);
  int paths[1 << 6][6]{};
  for (int last = 0; last < 6; ++last) paths[1 << last][last] = 1;
  for (int mask = 1; mask < (1 << 6); ++mask)
    for (int last = 0; last < 6; ++last) {
      if (!(mask & (1 << last))) continue;
      const int previous_mask = mask ^ (1 << last);
      for (int previous = 0; previous < 6; ++previous)
        if ((previous_mask & (1 << previous)) && edge[previous][last])
          paths[mask][last] += paths[previous_mask][previous];
    }
  int hamiltonian_paths = 0;
  for (int last = 0; last < 6; ++last)
    hamiltonian_paths += paths[(1 << 6) - 1][last];
  return {ties, flips, scores, triangles, hamiltonian_paths};
}

void check_histogram(const std::map<uint64_t, uint64_t> &actual,
                     const std::map<uint64_t, uint64_t> &expected,
                     const std::string &name) {
  require(actual == expected, name + " mismatch: " + format(actual));
}

}  // namespace

int main() {
  const Tables table = build_tables();
  Word scratch{};
  std::vector<Word> words;
  uint64_t state_words = 0;
  Sha256 grammar_hash;
  enumerate_words(words, scratch, 0, state_words, grammar_hash);
  const std::string grammar_digest = grammar_hash.finish();
  const std::vector<Support> supports = all_supports();

  std::map<uint64_t, uint64_t> scalar_feasible;
  std::map<uint64_t, uint64_t> support_histogram;
  std::string scalar_digest;
  const std::vector<Survivor> survivors = scalar_census(
      supports, words, table, scalar_feasible, support_histogram, scalar_digest);
  structural_scalar_audit(table, survivors);

  require(words.size() == 473, "hereditary word count");
  require(state_words == 243750000ULL, "state words per support");
  require(supports.size() == 924, "support count");
  require(survivors.size() == 36, "scalar survivor count");
  check_histogram(scalar_feasible,
                  {{0, 1156}, {1, 149868}, {2, 171636}, {3, 90884},
                   {4, 21864}, {5, 1608}, {6, 36}},
                  "scalar feasible-owner histogram");
  check_histogram(support_histogram, {{0, 888}, {1, 36}},
                  "support survivor histogram");

  require(table.mask_digest ==
              "741748b977fd90f0b506a15780e33d87bb04fdd60ae4f844ea1e40349ff8c47d",
          "primary mask digest");
  require(grammar_digest ==
              "7ae50439ddbd7e09d37516d067fe20f35d1f36b7830c25dc7156c900c6fde62f",
          "primary grammar digest");
  require(scalar_digest ==
              "ad266f55f820615eb2f7b4e323b6599024842c56352ff65c1b6dcdd117c250f9",
          "primary scalar digest");

  std::map<uint64_t, uint64_t> feasible_contexts;
  std::map<uint64_t, uint64_t> maximum_histogram;
  std::map<uint64_t, uint64_t> bank_size_histogram;
  std::map<uint64_t, uint64_t> maximizing_histogram;
  std::map<uint64_t, uint64_t> tie_histogram;
  std::map<uint64_t, uint64_t> flip_histogram;
  uint64_t total_masks = 0;
  uint64_t maximum_bank = 0;
  Sha256 owner_hash;
  for (const Survivor &row : survivors) {
    std::array<bool, 6> feasible{};
    std::array<int, 6> maxima{};
    for (int oi = 0; oi < 6; ++oi) {
      const int owner = row.support[oi];
      const std::vector<uint32_t> forward = owner_bank(row, owner, table, false);
      const std::vector<uint32_t> reverse = owner_bank(row, owner, table, true);
      require(forward == reverse, "forward/reverse reachable banks differ");
      feasible[oi] = std::binary_search(forward.begin(), forward.end(), FULL);
      for (uint32_t mask : forward)
        maxima[oi] = std::max(maxima[oi], std::popcount(mask));
      int maximizing = 0;
      for (uint32_t mask : forward)
        maximizing += std::popcount(mask) == maxima[oi];
      const int proof_bound = structural_owner_bound(row, owner, table);
      require(maxima[oi] <= proof_bound, "DP exceeds structural bound");
      ++maximum_histogram[maxima[oi]];
      ++bank_size_histogram[forward.size()];
      ++maximizing_histogram[maximizing];
      total_masks += forward.size();
      maximum_bank = std::max<uint64_t>(maximum_bank, forward.size());
      owner_hash.byte(owner);
      owner_hash.byte(feasible[oi]);
      owner_hash.byte(maxima[oi]);
      owner_hash.u64le(forward.size());
      owner_hash.u64le(maximizing);
      for (uint32_t mask : forward) owner_hash.u32le(mask);
    }
    ++feasible_contexts[std::count(feasible.begin(), feasible.end(), true)];
    const TournamentResult fingerprint =
        tournament(feasible, maxima, row.capacities);
    require(fingerprint.scores == std::array<int, 6>{0, 1, 2, 3, 4, 5} &&
                fingerprint.triangles == 0 &&
                fingerprint.hamiltonian_paths == 1,
            "nontransitive tournament fingerprint");
    ++tie_histogram[fingerprint.ties];
    ++flip_histogram[fingerprint.flips];
  }
  const std::string owner_digest = owner_hash.finish();
  require(owner_digest ==
              "a9064e057fcd7169395a28a3c411476f7b8178ec6137fbf74235aaee6db1f85c",
          "primary owner-bank digest");
  check_histogram(feasible_contexts, {{0, 36}},
                  "owner feasible-context histogram");
  check_histogram(maximum_histogram, {{21, 144}, {22, 72}},
                  "owner maximum histogram");
  check_histogram(bank_size_histogram,
                  {{45200, 24}, {48380, 24}, {48540, 24},
                   {133390, 48}, {140330, 48}, {141430, 48}},
                  "reachable-bank-size histogram");
  check_histogram(maximizing_histogram,
                  {{80, 24}, {90, 144}, {100, 24}, {140, 24}},
                  "maximizing-mask histogram");
  check_histogram(tie_histogram, {{7, 36}}, "tournament tie histogram");
  check_histogram(flip_histogram,
                  {{0, 1}, {1, 2}, {2, 5}, {3, 2}, {4, 16},
                   {5, 2}, {6, 5}, {7, 2}, {8, 1}},
                  "tournament flip histogram");
  require(total_masks == 23338080ULL && maximum_bank == 141430,
          "reachable-bank totals");

  std::cout << "scale-twenty-five independent literal-CRT C++ referee\n";
  std::cout << "method literal CRT; leave-one-out lcm grammar; discovered scalar rows; sorted-vector union DP; forward/reverse equality\n";
  std::cout << "orders 1,5,25 unit-counts 1,4,20\n";
  std::cout << "hereditary words " << words.size()
            << " state words/support " << state_words
            << " raw labelled states " << supports.size() * state_words << '\n';
  std::cout << "literal CRT-base SHA256 " << table.base_digest << '\n';
  std::cout << "primary mask SHA256 " << table.mask_digest << '\n';
  std::cout << "primary grammar SHA256 " << grammar_digest << '\n';
  std::cout << "scalar feasible-owner histogram " << format(scalar_feasible) << '\n';
  std::cout << "scalar support histogram " << format(support_histogram) << '\n';
  std::cout << "scalar survivors " << survivors.size() << '\n';
  std::cout << "primary scalar SHA256 " << scalar_digest << '\n';
  std::cout << "structural scalar audit forbidden Cayley relation 36 arcs; 30 edges; 6 reciprocal pairs; K6+K6 on QR/NQR\n";
  std::cout << "structural owner audit distinct D5 fibres force four needle incidences; coincident D5 fibres have total union at most 20\n";
  std::cout << "owner feasible-context histogram " << format(feasible_contexts) << '\n';
  std::cout << "owner maximum-union histogram " << format(maximum_histogram) << '\n';
  std::cout << "owner reachable-bank-size histogram " << format(bank_size_histogram) << '\n';
  std::cout << "owner maximizing-mask histogram " << format(maximizing_histogram) << '\n';
  std::cout << "owner rows " << 6 * survivors.size()
            << " total reachable masks " << total_masks
            << " maximum bank " << maximum_bank << '\n';
  std::cout << "primary owner SHA256 " << owner_digest << '\n';
  std::cout << "tournament observable (feasible,max-union,capacity), coordinate tie switch\n";
  std::cout << "tournament all transitive; scores 0,1,2,3,4,5; cycles 0; Hamiltonian paths 1; tie histogram "
            << format(tie_histogram) << " flip histogram " << format(flip_histogram) << '\n';
  std::cout << "verdict every one of the 36 scalar contexts has zero feasible owners; scale 25 empty\n";
  return 0;
}
