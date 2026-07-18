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
#include <unordered_map>
#include <utility>
#include <vector>

// Independent literal-CRT referee for the primitive proper AP-centred
// common-scale-twenty-four Hamming-six face (THM-990).
//
// Design was frozen before the Python primary or its output was inspected.
// This implementation deliberately uses a different carrier:
//   * orders and exact-order residues are generated from gcd/lcm definitions;
//   * every CRT representative is found by bounded literal search;
//   * the six owner-local unit fibres are joined by sorted-vector union DP;
//   * multiplication covariance is audited but never used to remove rows.
//
// Only the final, commutative owner-local DP is memoized.  Its key is the
// exact multiset of (provider/owner ratio, provider order) pairs.  The table
// audit below proves that this is an exact change of sheet gauge, not a
// quotient of the owner predicate.

namespace {

constexpr int P = 13;
constexpr int C = 24;
constexpr uint32_t FULL = (uint32_t{1} << C) - 1;
constexpr std::array<uint8_t, 8> ORDERS{1, 2, 3, 4, 6, 8, 12, 24};

using Support = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;  // indices into ORDERS
using Capacities = std::array<uint8_t, 6>;
using Multiplicity = std::array<uint8_t, 8>;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

// Deterministic FNV-1a byte serialization.  File-level SHA-256 hashes are
// recorded separately; these internal digests make the exact enumerated banks
// cheap to compare across compiler modes.
struct Fnv64 {
  uint64_t state = 14695981039346656037ULL;

  void byte(uint8_t value) {
    state ^= value;
    state *= 1099511628211ULL;
  }
  void u32(uint32_t value) {
    for (int shift = 0; shift < 32; shift += 8)
      byte(static_cast<uint8_t>(value >> shift));
  }
  void u64(uint64_t value) {
    for (int shift = 0; shift < 64; shift += 8)
      byte(static_cast<uint8_t>(value >> shift));
  }
};

std::string hex64(uint64_t value) {
  std::ostringstream out;
  out << std::hex << std::setfill('0') << std::setw(16) << value;
  return out.str();
}

int inverse_mod_13(int value) {
  for (int candidate = 1; candidate < P; ++candidate)
    if (value * candidate % P == 1) return candidate;
  fail("attempted to invert a nonunit modulo thirteen");
}

int centered(int value, int modulus) {
  int residue = value % modulus;
  if (residue < 0) residue += modulus;
  return 2 * residue > modulus ? residue - modulus : residue;
}

int literal_crt_base(int label, int order, int unit) {
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P &&
        value % order == unit % order)
      return value;
  fail("bounded CRT search found no representative");
}

uint32_t literal_mask(int label, int order, int unit, int owner) {
  const int base = literal_crt_base(label, order, unit);
  const int inverse = inverse_mod_13(owner);
  uint32_t mask = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    if (-order < value && value <= order) mask |= uint32_t{1} << sheet;
  }
  return mask;
}

uint32_t rotate24(uint32_t mask, int amount) {
  amount %= C;
  if (amount < 0) amount += C;
  if (amount == 0) return mask & FULL;
  return ((mask << amount) | (mask >> (C - amount))) & FULL;
}

uint32_t residue_class_mask(int residue) {
  uint32_t result = 0;
  for (int sheet = residue; sheet < C; sheet += 3)
    result |= uint32_t{1} << sheet;
  return result;
}

std::vector<uint8_t> units_mod_order(int order) {
  if (order == 1) return {0};
  std::vector<uint8_t> result;
  for (int unit = 1; unit < order; ++unit)
    if (std::gcd(unit, order) == 1)
      result.push_back(static_cast<uint8_t>(unit));
  return result;
}

struct Tables {
  std::array<std::vector<uint8_t>, 8> units;
  // Normalized owner-one fibres [ratio][order index].
  std::array<std::array<std::vector<uint32_t>, 8>, P> fibres;
  // Unit-independent cardinalities [provider][order index][owner].
  std::array<std::array<std::array<uint8_t, P>, 8>, P> cards{};
  uint64_t base_digest = 0;
  uint64_t mask_digest = 0;
};

Tables build_tables() {
  Tables table;
  Fnv64 bases;
  Fnv64 masks;
  for (int oi = 0; oi < 8; ++oi) {
    table.units[oi] = units_mod_order(ORDERS[oi]);
  }

  // Direct gcd census, checked against the explicit Euler-phi values.
  const std::array<std::size_t, 8> expected_phi{1, 1, 2, 2, 2, 4, 4, 8};
  for (int oi = 0; oi < 8; ++oi)
    require(table.units[oi].size() == expected_phi[oi],
            "exact-order residue census mismatch");

  for (int ratio = 1; ratio < P; ++ratio)
    for (int oi = 0; oi < 8; ++oi) {
      const int order = ORDERS[oi];
      for (uint8_t unit : table.units[oi])
        table.fibres[ratio][oi].push_back(
            literal_mask(ratio, order, unit, 1));
      std::sort(table.fibres[ratio][oi].begin(),
                table.fibres[ratio][oi].end());
      table.fibres[ratio][oi].erase(
          std::unique(table.fibres[ratio][oi].begin(),
                      table.fibres[ratio][oi].end()),
          table.fibres[ratio][oi].end());
      require(!table.fibres[ratio][oi].empty(),
              "an exact-order residue fibre has empty mask image");
    }

  for (int label = 1; label < P; ++label)
    for (int oi = 0; oi < 8; ++oi) {
      const int order = ORDERS[oi];
      for (uint8_t unit : table.units[oi]) {
        bases.byte(static_cast<uint8_t>(label));
        bases.byte(static_cast<uint8_t>(oi));
        bases.byte(unit);
        bases.u32(static_cast<uint32_t>(literal_crt_base(label, order, unit)));
      }
      for (int owner = 1; owner < P; ++owner) {
        const int ratio = label * inverse_mod_13(owner) % P;
        int shift = -1;
        for (int candidate = 0; candidate < C; ++candidate)
          if ((inverse_mod_13(owner) + P * candidate) % C == 1) {
            shift = candidate;
            break;
          }
        require(shift >= 0, "owner sheet-gauge shift does not exist");
        int common_card = -1;
        for (std::size_t ui = 0; ui < table.units[oi].size(); ++ui) {
          const uint8_t unit = table.units[oi][ui];
          const uint32_t actual = literal_mask(label, order, unit, owner);
          const uint32_t normalized = literal_mask(ratio, order, unit, 1);
          require(actual == rotate24(normalized, shift),
                  "literal owner mask violates the common cyclic gauge");
          const int card = std::popcount(actual);
          if (common_card < 0) common_card = card;
          require(card == common_card,
                  "local sheet cardinality depends on exact-order residue");

          // A direct one-period scan is independent of the 24-sheet count.
          const int base = literal_crt_base(ratio, order, unit);
          int period_card = 0;
          for (int sheet = 0; sheet < order; ++sheet) {
            const int value = centered(base * (1 + P * sheet), P * order);
            period_card += -order < value && value <= order;
          }
          require(card == (C / order) * period_card,
                  "full sheet scan disagrees with one-period count");

          masks.byte(static_cast<uint8_t>(label));
          masks.byte(static_cast<uint8_t>(oi));
          masks.byte(static_cast<uint8_t>(owner));
          masks.byte(unit);
          masks.u32(actual);
        }
        table.cards[label][oi][owner] =
            static_cast<uint8_t>(common_card);
      }
    }
  table.base_digest = bases.state;
  table.mask_digest = masks.state;
  return table;
}

uint64_t audit_cubic_sheet_nerve(const Tables &table) {
  const uint32_t class0 = residue_class_mask(0);
  const uint32_t class1 = residue_class_mask(1);
  const uint32_t class2 = residue_class_mask(2);
  const std::vector<uint32_t> classes01{class0, class1};
  const std::vector<uint32_t> empty{0};
  require(table.fibres[1][2] == std::vector<uint32_t>{class2},
          "D3 ratio-one fibre is not sheet class two");
  require(table.fibres[5][2] == classes01 &&
              table.fibres[8][2] == classes01,
          "D3 active C0 fibres are not sheet classes zero and one");
  require(table.fibres[12][2] == empty,
          "D3 inactive C0 fibre is not empty");
  require(table.fibres[4][2] == classes01 &&
              table.fibres[9][2] == classes01 &&
              table.fibres[6][2] == empty &&
              table.fibres[7][2] == empty,
          "D3 C2 fibres do not have the claimed two-active/two-empty split");

  Fnv64 digest;
  for (const auto &[oi, partner] :
       std::vector<std::pair<int, int>>{{5, 5}, {5, 8}, {5, 12},
                                        {7, 5}, {7, 8}}) {
    const int order = ORDERS[oi];
    const int expected_left_card = order == 8 ? 6 : 4;
    const int expected_right_card = order == 8 ? 3 : 4;
    for (uint32_t mask : table.fibres[1][oi])
      require(std::popcount(mask) == expected_left_card,
              "high-order ratio-one mask cardinality mismatch");
    for (uint32_t mask : table.fibres[partner][oi])
      require(std::popcount(mask) == expected_right_card,
              "high-order partner mask cardinality mismatch");
    for (int residue = 0; residue < 3; ++residue) {
      const uint32_t sheet_class = residue_class_mask(residue);
      int maximum = 0;
      int witnesses = 0;
      for (uint32_t left : table.fibres[1][oi])
        for (uint32_t right : table.fibres[partner][oi]) {
          const int card = std::popcount((left | right) & sheet_class);
          if (card > maximum) {
            maximum = card;
            witnesses = 1;
          } else if (card == maximum) {
            ++witnesses;
          }
        }
      require(maximum == 3,
              "high-order pair adds more or less than three points to a "
              "mod-three sheet class");
      digest.byte(static_cast<uint8_t>(order));
      digest.byte(static_cast<uint8_t>(partner));
      digest.byte(static_cast<uint8_t>(residue));
      digest.byte(static_cast<uint8_t>(maximum));
      digest.u32(static_cast<uint32_t>(witnesses));
    }
  }
  return digest.state;
}

bool hereditary_prime_power(const OrderWord &word) {
  int carries_eight = 0;
  int carries_three = 0;
  for (uint8_t oi : word) {
    const int order = ORDERS[oi];
    carries_eight += order % 8 == 0;
    carries_three += order % 3 == 0;
  }
  return carries_eight >= 2 && carries_three >= 2;
}

bool hereditary_lcm(const OrderWord &word) {
  for (int omitted = 0; omitted < 6; ++omitted) {
    int residual = 1;
    for (int coordinate = 0; coordinate < 6; ++coordinate)
      if (coordinate != omitted)
        residual = std::lcm(residual,
                            static_cast<int>(ORDERS[word[coordinate]]));
    if (residual != C) return false;
  }
  return true;
}

void enumerate_words(int coordinate, OrderWord &word,
                     const Tables &table, std::vector<OrderWord> &accepted,
                     uint64_t &state_words, Fnv64 &grammar_digest) {
  if (coordinate == 6) {
    const bool prime_power = hereditary_prime_power(word);
    require(prime_power == hereditary_lcm(word),
            "prime-power and leave-one-out-lcm grammars disagree");
    if (!prime_power) return;
    uint64_t fibre = 1;
    for (uint8_t oi : word) {
      fibre *= table.units[oi].size();
      grammar_digest.byte(oi);
    }
    grammar_digest.u64(fibre);
    state_words += fibre;
    accepted.push_back(word);
    return;
  }
  for (uint8_t oi = 0; oi < 8; ++oi) {
    word[coordinate] = oi;
    enumerate_words(coordinate + 1, word, table, accepted, state_words,
                    grammar_digest);
  }
}

std::vector<Support> all_supports() {
  std::vector<Support> result;
  Support support{};
  const auto visit = [&](auto &&self, int coordinate, int next) -> void {
    if (coordinate == 6) {
      result.push_back(support);
      return;
    }
    for (int label = next; label <= 12 - (5 - coordinate); ++label) {
      support[coordinate] = static_cast<uint8_t>(label);
      self(self, coordinate + 1, label + 1);
    }
  };
  visit(visit, 0, 1);
  return result;
}

Multiplicity multiplicity(const OrderWord &word) {
  Multiplicity result{};
  for (uint8_t oi : word) ++result[oi];
  return result;
}

struct ScalarRow {
  Support support{};
  OrderWord word{};
  Capacities capacities{};
};

struct FourFeasibleRow {
  Support support{};
  OrderWord word{};
  Capacities capacities{};
  std::array<uint8_t, 6> maxima{};
  uint8_t feasible_mask = 0;
};

uint64_t scalar_row_key(const Support &support, const OrderWord &word) {
  uint64_t result = 0;
  for (int i = 0; i < 6; ++i) {
    result <<= 7;
    result |= static_cast<uint64_t>((support[i] << 3) | word[i]);
  }
  return result;
}

uint64_t multiply_row_key(const Support &support, const OrderWord &word,
                          int multiplier) {
  std::array<std::pair<uint8_t, uint8_t>, 6> pairs{};
  for (int i = 0; i < 6; ++i)
    pairs[i] = {static_cast<uint8_t>(support[i] * multiplier % P), word[i]};
  std::sort(pairs.begin(), pairs.end());
  Support transformed_support{};
  OrderWord transformed_word{};
  for (int i = 0; i < 6; ++i) {
    transformed_support[i] = pairs[i].first;
    transformed_word[i] = pairs[i].second;
  }
  return scalar_row_key(transformed_support, transformed_word);
}

struct LocalSummary {
  bool feasible = false;
  uint8_t maximum = 0;
  uint32_t reachable_count = 0;
  uint32_t maximum_count = 0;
  uint32_t largest_layer = 0;
  uint64_t reachable_digest = 0;
};

std::vector<uint32_t> union_dp(
    const std::array<uint8_t, 6> &encoded_pairs, const Tables &table,
    bool reverse, uint32_t &largest_layer) {
  std::vector<uint32_t> reachable{0};
  largest_layer = 1;
  for (int step = 0; step < 6; ++step) {
    const int provider = reverse ? 5 - step : step;
    const int ratio = encoded_pairs[provider] >> 3;
    const int oi = encoded_pairs[provider] & 7;
    const std::vector<uint32_t> &fibre = table.fibres[ratio][oi];
    std::vector<uint32_t> next;
    next.reserve(reachable.size() * fibre.size());
    for (uint32_t partial : reachable)
      for (uint32_t mask : fibre) next.push_back(partial | mask);
    std::sort(next.begin(), next.end());
    next.erase(std::unique(next.begin(), next.end()), next.end());
    largest_layer =
        std::max(largest_layer, static_cast<uint32_t>(next.size()));
    reachable = std::move(next);
  }
  return reachable;
}

uint64_t local_key(const Support &support, const OrderWord &word, int owner) {
  std::array<uint8_t, 6> pairs{};
  const int inverse = inverse_mod_13(owner);
  for (int provider = 0; provider < 6; ++provider) {
    const int ratio = support[provider] * inverse % P;
    pairs[provider] = static_cast<uint8_t>((ratio << 3) | word[provider]);
  }
  std::sort(pairs.begin(), pairs.end());
  uint64_t result = 0;
  for (uint8_t pair : pairs) result = (result << 7) | pair;
  return result;
}

std::array<uint8_t, 6> decode_local_key(uint64_t key) {
  std::array<uint8_t, 6> result{};
  for (int i = 5; i >= 0; --i) {
    result[i] = static_cast<uint8_t>(key & 0x7fU);
    key >>= 7;
  }
  require(key == 0, "local key overflow");
  return result;
}

LocalSummary compute_local(uint64_t key, const Tables &table) {
  const std::array<uint8_t, 6> pairs = decode_local_key(key);
  uint32_t forward_largest = 0;
  uint32_t reverse_largest = 0;
  const std::vector<uint32_t> forward =
      union_dp(pairs, table, false, forward_largest);
  const std::vector<uint32_t> reverse =
      union_dp(pairs, table, true, reverse_largest);
  require(forward == reverse,
          "forward and reverse sorted-vector joins disagree");

  LocalSummary result;
  result.reachable_count = static_cast<uint32_t>(forward.size());
  result.largest_layer = std::max(forward_largest, reverse_largest);
  Fnv64 digest;
  for (uint32_t mask : forward) {
    digest.u32(mask);
    const int card = std::popcount(mask);
    if (card > result.maximum) {
      result.maximum = static_cast<uint8_t>(card);
      result.maximum_count = 1;
    } else if (card == result.maximum) {
      ++result.maximum_count;
    }
    result.feasible |= mask == FULL;
  }
  result.reachable_digest = digest.state;
  require(result.feasible == (result.maximum == C),
          "FULL-mask and maximum-union feasibility disagree");
  return result;
}

struct TournamentSummary {
  int ties = 0;
  int flips = 0;
  int triangles = 0;
  int sccs = 0;
  uint64_t paths = 0;
  std::array<uint8_t, 6> scores{};
};

TournamentSummary tournament(const Capacities &capacities,
                             const std::array<LocalSummary, 6> &locals) {
  std::array<uint8_t, 6> out{};
  TournamentSummary result;
  for (int left = 0; left < 6; ++left)
    for (int right = left + 1; right < 6; ++right) {
      const auto left_key = std::tuple{
          locals[left].feasible, locals[left].maximum, capacities[left],
          locals[left].reachable_count, locals[left].maximum_count};
      const auto right_key = std::tuple{
          locals[right].feasible, locals[right].maximum, capacities[right],
          locals[right].reachable_count, locals[right].maximum_count};
      int winner = left;
      if (left_key == right_key) {
        ++result.ties;  // coordinate order is the tie Hamiltonian path
      } else if (right_key > left_key) {
        winner = right;
        ++result.flips;
      }
      const int loser = left + right - winner;
      out[winner] |= static_cast<uint8_t>(1U << loser);
    }

  for (int vertex = 0; vertex < 6; ++vertex)
    result.scores[vertex] =
        static_cast<uint8_t>(std::popcount(out[vertex]));
  for (int a = 0; a < 6; ++a)
    for (int b = a + 1; b < 6; ++b)
      for (int c = b + 1; c < 6; ++c) {
        const bool forward = ((out[a] >> b) & 1U) &&
                             ((out[b] >> c) & 1U) &&
                             ((out[c] >> a) & 1U);
        const bool reverse = ((out[a] >> c) & 1U) &&
                             ((out[c] >> b) & 1U) &&
                             ((out[b] >> a) & 1U);
        result.triangles += forward || reverse;
      }

  std::array<uint8_t, 6> reach = out;
  for (int vertex = 0; vertex < 6; ++vertex)
    reach[vertex] |= static_cast<uint8_t>(1U << vertex);
  for (int middle = 0; middle < 6; ++middle)
    for (int source = 0; source < 6; ++source)
      if ((reach[source] >> middle) & 1U) reach[source] |= reach[middle];
  uint8_t assigned = 0;
  for (int root = 0; root < 6; ++root)
    if (!((assigned >> root) & 1U)) {
      ++result.sccs;
      for (int vertex = 0; vertex < 6; ++vertex)
        if (((reach[root] >> vertex) & 1U) &&
            ((reach[vertex] >> root) & 1U))
          assigned |= static_cast<uint8_t>(1U << vertex);
    }

  std::array<std::array<uint64_t, 6>, 1U << 6> paths{};
  for (int last = 0; last < 6; ++last) paths[1U << last][last] = 1;
  for (int mask = 1; mask < (1 << 6); ++mask)
    for (int last = 0; last < 6; ++last)
      if ((mask >> last) & 1) {
        const int previous_mask = mask ^ (1 << last);
        for (int previous = 0; previous < 6; ++previous)
          if (((previous_mask >> previous) & 1) &&
              ((out[previous] >> last) & 1U))
            paths[mask][last] += paths[previous_mask][previous];
      }
  for (int last = 0; last < 6; ++last) result.paths += paths.back()[last];
  return result;
}

template <typename Key>
std::string histogram(const std::map<Key, uint64_t> &values) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, count] : values) {
    if (!first) out << ' ';
    first = false;
    out << key << ':' << count;
  }
  return out.str();
}

std::string feasible_mask_histogram(const std::map<int, uint64_t> &values) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, count] : values) {
    if (!first) out << ' ';
    first = false;
    out << "0x" << std::hex << key << std::dec << ':' << count;
  }
  return out.str();
}

template <typename Values>
std::string comma_values(const Values &values) {
  std::ostringstream out;
  bool first = true;
  for (const auto value : values) {
    if (!first) out << ',';
    first = false;
    out << static_cast<int>(value);
  }
  return out.str();
}

int multiplicative_exponent(int label) {
  int value = 1;
  for (int exponent = 0; exponent < 12; ++exponent) {
    if (value == label) return exponent;
    value = value * 2 % P;
  }
  fail("nonzero label is absent from the primitive-root table");
}

int unordered_ratio_distance(int left, int right) {
  int difference =
      (multiplicative_exponent(right) - multiplicative_exponent(left) + 12) %
      12;
  return std::min(difference, 12 - difference);
}

}  // namespace

int main(int argc, char **argv) {
  const Tables table = build_tables();
  const uint64_t cubic_nerve_digest = audit_cubic_sheet_nerve(table);

  std::vector<OrderWord> words;
  OrderWord scratch{};
  uint64_t state_words_per_support = 0;
  Fnv64 grammar_digest;
  enumerate_words(0, scratch, table, words, state_words_per_support,
                  grammar_digest);
  require(words.size() == 108'813,
          "hereditary divisor-word census mismatch");
  require(state_words_per_support == 167'165'952ULL,
          "literal state-word census mismatch");

  if (argc == 2 && std::string(argv[1]) == "--preflight") {
    std::cout << "scale-24 referee preflight green: literal CRT/mask gauge, "
                 "one-period cardinalities, hereditary grammar "
              << words.size() << ", weighted state words/support "
              << state_words_per_support << '\n';
    return 0;
  }
  require(argc == 1, "usage: referee [--preflight]");

  const std::vector<Support> supports = all_supports();
  require(supports.size() == 924, "six-support census mismatch");

  std::vector<ScalarRow> scalar_bank;
  scalar_bank.reserve(70'000);
  std::set<Support> scalar_supports;
  std::map<Support, int> contexts_per_support;
  std::map<int, uint64_t> all_support_context_histogram;
  std::map<Multiplicity, uint64_t> multiplicity_histogram;
  std::set<Capacities> capacity_vectors;
  Fnv64 scalar_digest;
  Fnv64 capacity_digest;

  for (const Support &support : supports) {
    // Six independent eight-bit lanes fit without carry: each is a sum of
    // six mask cardinalities and is at most 6*C.
    std::array<std::array<uint64_t, 8>, 6> contribution{};
    for (int provider = 0; provider < 6; ++provider)
      for (int oi = 0; oi < 8; ++oi)
        for (int owner_index = 0; owner_index < 6; ++owner_index)
          contribution[provider][oi] |=
              static_cast<uint64_t>(table.cards[support[provider]][oi]
                                               [support[owner_index]])
              << (8 * owner_index);

    int support_contexts = 0;
    for (const OrderWord &word : words) {
      uint64_t packed = 0;
      for (int provider = 0; provider < 6; ++provider)
        packed += contribution[provider][word[provider]];
      Capacities capacities{};
      bool scalar = true;
      for (int owner_index = 0; owner_index < 6; ++owner_index) {
        capacities[owner_index] =
            static_cast<uint8_t>(packed >> (8 * owner_index));
        scalar &= capacities[owner_index] >= C;
      }
      if (!scalar) continue;

      scalar_bank.push_back(ScalarRow{support, word, capacities});
      scalar_supports.insert(support);
      ++contexts_per_support[support];
      ++support_contexts;
      ++multiplicity_histogram[multiplicity(word)];
      capacity_vectors.insert(capacities);
      for (uint8_t label : support) {
        scalar_digest.byte(label);
        capacity_digest.byte(label);
      }
      for (uint8_t oi : word) {
        scalar_digest.byte(oi);
        capacity_digest.byte(oi);
      }
      for (uint8_t capacity : capacities) capacity_digest.byte(capacity);
    }
    ++all_support_context_histogram[support_contexts];
  }

  require(scalar_bank.size() == 66'984,
          "scalar survivor census disagrees with THM-990 primary");
  require(scalar_supports.size() == 854,
          "scalar-support census disagrees with THM-990 primary");
  require(multiplicity_histogram.size() == 202,
          "scalar multiplicity-profile census disagrees with THM-990");

  // Exact multiplication covariance of every labelled (support,order) row.
  std::set<uint64_t> scalar_row_keys;
  for (const ScalarRow &row : scalar_bank)
    scalar_row_keys.insert(scalar_row_key(row.support, row.word));
  require(scalar_row_keys.size() == scalar_bank.size(),
          "duplicate scalar row encoding");
  for (const ScalarRow &row : scalar_bank)
    for (int multiplier = 1; multiplier < P; ++multiplier)
      require(scalar_row_keys.contains(
                  multiply_row_key(row.support, row.word, multiplier)),
              "scalar row bank violates multiplication covariance");

  std::unordered_map<uint64_t, LocalSummary> local_cache;
  local_cache.reserve(50'000);
  std::map<int, uint64_t> feasible_contexts;
  std::map<int, uint64_t> feasible_masks;
  std::map<int, uint64_t> maximum_union;
  std::map<int, uint64_t> reachable_counts;
  std::map<int, uint64_t> maximum_mask_counts;
  std::map<int, uint64_t> largest_layers;
  std::map<int, uint64_t> tournament_ties;
  std::map<int, uint64_t> tournament_flips;
  std::map<int, uint64_t> tournament_scores;
  std::set<std::array<uint8_t, 6>> owner_maximum_vectors;
  Fnv64 owner_bank_digest;
  Fnv64 tournament_digest;
  uint64_t feasible_owner_rows = 0;
  uint64_t total_reachable_masks = 0;
  uint32_t largest_reachable_bank = 0;
  uint32_t largest_intermediate_layer = 0;
  std::vector<FourFeasibleRow> four_feasible_rows;

  for (const ScalarRow &row : scalar_bank) {
    std::array<LocalSummary, 6> locals{};
    std::array<uint8_t, 6> owner_vector{};
    int feasible_count = 0;
    int feasible_mask = 0;
    for (int owner_index = 0; owner_index < 6; ++owner_index) {
      const int owner = row.support[owner_index];
      const uint64_t key = local_key(row.support, row.word, owner);
      auto found = local_cache.find(key);
      if (found == local_cache.end())
        found = local_cache.emplace(key, compute_local(key, table)).first;
      locals[owner_index] = found->second;
      const LocalSummary &local = locals[owner_index];
      feasible_count += local.feasible;
      feasible_mask |= static_cast<int>(local.feasible) << owner_index;
      feasible_owner_rows += local.feasible;
      total_reachable_masks += local.reachable_count;
      owner_vector[owner_index] = local.maximum;
      ++maximum_union[local.maximum];
      ++reachable_counts[local.reachable_count];
      ++maximum_mask_counts[local.maximum_count];
      ++largest_layers[local.largest_layer];
      largest_reachable_bank =
          std::max(largest_reachable_bank, local.reachable_count);
      largest_intermediate_layer =
          std::max(largest_intermediate_layer, local.largest_layer);

      for (uint8_t label : row.support) owner_bank_digest.byte(label);
      for (uint8_t oi : row.word) owner_bank_digest.byte(oi);
      owner_bank_digest.byte(static_cast<uint8_t>(owner_index));
      owner_bank_digest.byte(static_cast<uint8_t>(owner));
      owner_bank_digest.byte(static_cast<uint8_t>(local.feasible));
      owner_bank_digest.byte(local.maximum);
      owner_bank_digest.byte(row.capacities[owner_index]);
      owner_bank_digest.u32(local.reachable_count);
      owner_bank_digest.u32(local.maximum_count);
      owner_bank_digest.u64(local.reachable_digest);
    }
    ++feasible_contexts[feasible_count];
    ++feasible_masks[feasible_mask];
    owner_maximum_vectors.insert(owner_vector);
    owner_bank_digest.byte(static_cast<uint8_t>(feasible_mask));
    for (uint8_t value : owner_vector) owner_bank_digest.byte(value);
    if (feasible_count == 4)
      four_feasible_rows.push_back(FourFeasibleRow{
          row.support, row.word, row.capacities, owner_vector,
          static_cast<uint8_t>(feasible_mask)});

    const TournamentSummary audit = tournament(row.capacities, locals);
    ++tournament_ties[audit.ties];
    ++tournament_flips[audit.flips];
    for (uint8_t score : audit.scores) ++tournament_scores[score];
    require(audit.triangles == 0 && audit.sccs == 6 && audit.paths == 1,
            "owner-obligation tournament is not transitive");
    std::array<uint8_t, 6> sorted_scores = audit.scores;
    std::sort(sorted_scores.begin(), sorted_scores.end());
    require(sorted_scores == std::array<uint8_t, 6>{0, 1, 2, 3, 4, 5},
            "owner-obligation tournament score word mismatch");
    tournament_digest.byte(static_cast<uint8_t>(audit.ties));
    tournament_digest.byte(static_cast<uint8_t>(audit.flips));
    for (uint8_t score : audit.scores) tournament_digest.byte(score);
  }

  const std::map<int, uint64_t> expected_feasible{
      {0, 64'962}, {1, 1'800}, {2, 192}, {4, 30}};
  const std::map<int, uint64_t> expected_maximum{
      {12, 72},    {14, 2'136}, {15, 1'644},  {16, 15'876},
      {17, 24'420}, {18, 76'296}, {19, 94'872}, {20, 104'592},
      {21, 53'040}, {22, 24'948}, {23, 1'704},  {24, 2'304}};
  require(feasible_contexts == expected_feasible,
          "feasible-owner histogram disagrees with THM-990 primary");
  require(maximum_union == expected_maximum,
          "maximum-union histogram disagrees with THM-990 primary");
  require(owner_maximum_vectors.size() == 20'302,
          "owner maximum-vector census disagrees with THM-990 primary");
  require(largest_reachable_bank == 7'728,
          "largest reachable-mask bank disagrees with THM-990 primary");
  require(total_reachable_masks == 101'961'528ULL &&
              reachable_counts.size() == 674,
          "reachable-mask total/bin census disagrees with THM-990 primary");
  require(feasible_owner_rows == 2'304,
          "feasible owner-row total does not follow the context histogram");
  require(four_feasible_rows.size() == 30,
          "four-feasible exceptional-row census mismatch");

  std::map<Multiplicity, uint64_t> four_order_profiles;
  std::map<int, uint64_t> missing_ratio_distances;
  std::map<std::pair<int, int>, uint64_t> missing_edges;
  std::map<std::pair<int, int>, uint64_t> four_types;
  std::array<uint64_t, P> missing_degrees{};
  std::set<uint64_t> four_row_keys;
  Fnv64 four_row_digest;
  for (const FourFeasibleRow &row : four_feasible_rows) {
    ++four_order_profiles[multiplicity(row.word)];
    four_row_keys.insert(scalar_row_key(row.support, row.word));
    std::array<int, 2> missing{};
    int next_missing = 0;
    for (int owner_index = 0; owner_index < 6; ++owner_index)
      if (!((row.feasible_mask >> owner_index) & 1U))
        missing[next_missing++] = row.support[owner_index];
    require(next_missing == 2, "four-feasible row does not miss two owners");
    if (missing[1] < missing[0]) std::swap(missing[0], missing[1]);
    ++missing_edges[{missing[0], missing[1]}];
    ++missing_degrees[missing[0]];
    ++missing_degrees[missing[1]];
    const int distance = unordered_ratio_distance(missing[0], missing[1]);
    ++missing_ratio_distances[distance];

    int high_order = -1;
    std::set<int> low_labels;
    std::set<int> high_labels;
    for (int owner_index = 0; owner_index < 6; ++owner_index) {
      const int order = ORDERS[row.word[owner_index]];
      if ((row.feasible_mask >> owner_index) & 1U) {
        require(order == 3 && row.maxima[owner_index] == 24,
                "four-feasible row has a non-D3 or non-full low owner");
        low_labels.insert(row.support[owner_index]);
      } else {
        if (high_order < 0) high_order = order;
        require(order == high_order && row.maxima[owner_index] == 19,
                "four-feasible row high owners do not share order/max 19");
        high_labels.insert(row.support[owner_index]);
      }
    }
    require(high_order == 8 || high_order == 24,
            "four-feasible high order is neither eight nor twenty-four");
    for (int owner_index = 0; owner_index < 6; ++owner_index) {
      const bool low = (row.feasible_mask >> owner_index) & 1U;
      const int expected_capacity =
          high_order == 8 ? (low ? 30 : 25) : (low ? 32 : 24);
      require(row.capacities[owner_index] == expected_capacity,
              "four-feasible capacity vector violates the cubic type law");
    }
    require((distance == 3 && (high_order == 8 || high_order == 24)) ||
                (distance == 6 && high_order == 8),
            "four-feasible high pair has an unclassified ratio/order type");
    ++four_types[{high_order, distance}];

    // Taking either high label as a, the support is exactly
    //     a*C2 union E,  E subset a*C0,
    // while a*C1 is the absent cubic coset.  Here the primitive root is 2.
    const int a = missing[0];
    const std::array<int, 4> C0{1, 5, 8, 12};
    const std::array<int, 4> C1{2, 3, 10, 11};
    const std::array<int, 4> C2{4, 6, 7, 9};
    std::set<int> expected_high_coset;
    std::set<int> expected_absent_coset;
    std::set<int> expected_low_coset;
    for (int value : C0) expected_high_coset.insert(a * value % P);
    for (int value : C1) expected_absent_coset.insert(a * value % P);
    for (int value : C2) expected_low_coset.insert(a * value % P);
    require(low_labels == expected_low_coset,
            "four-feasible low owners are not the oriented cubic coset a*C2");
    require(std::includes(expected_high_coset.begin(),
                          expected_high_coset.end(), high_labels.begin(),
                          high_labels.end()),
            "four-feasible high pair is not contained in a*C0");
    for (uint8_t label : row.support)
      require(!expected_absent_coset.contains(label),
              "four-feasible support meets the absent cubic coset a*C1");

    for (uint8_t label : row.support) four_row_digest.byte(label);
    for (uint8_t oi : row.word) four_row_digest.byte(oi);
    for (uint8_t capacity : row.capacities) four_row_digest.byte(capacity);
    for (uint8_t maximum : row.maxima) four_row_digest.byte(maximum);
    four_row_digest.byte(row.feasible_mask);
    four_row_digest.byte(static_cast<uint8_t>(distance));
  }
  require(four_types ==
              std::map<std::pair<int, int>, uint64_t>{{{8, 3}, 12},
                                                       {{8, 6}, 6},
                                                       {{24, 3}, 12}},
          "four-feasible cubic type census mismatch");

  // Multiplication must preserve not merely scalar survival but this extremal
  // four-feasible stratum.  Orbit sizes are computed on the labelled rows.
  std::set<uint64_t> four_remaining = four_row_keys;
  std::map<int, uint64_t> four_orbit_sizes;
  while (!four_remaining.empty()) {
    const uint64_t seed = *four_remaining.begin();
    const auto found = std::find_if(
        four_feasible_rows.begin(), four_feasible_rows.end(),
        [&](const FourFeasibleRow &row) {
          return scalar_row_key(row.support, row.word) == seed;
        });
    require(found != four_feasible_rows.end(),
            "four-feasible orbit seed has no row");
    std::set<uint64_t> orbit;
    for (int multiplier = 1; multiplier < P; ++multiplier)
      orbit.insert(multiply_row_key(found->support, found->word, multiplier));
    require(std::includes(four_row_keys.begin(), four_row_keys.end(),
                          orbit.begin(), orbit.end()),
            "four-feasible stratum violates multiplication covariance");
    require(std::includes(four_remaining.begin(), four_remaining.end(),
                          orbit.begin(), orbit.end()),
            "four-feasible multiplication orbits overlap");
    ++four_orbit_sizes[orbit.size()];
    for (uint64_t key : orbit) four_remaining.erase(key);
  }

  std::cout << "THM-990 scale-24 literal-CRT sorted-vector referee\n";
  std::cout << "scope: primitive proper AP-centred common-scale Hamming-six "
               "owner-local gate only\n";
  std::cout << "supports " << supports.size() << "; hereditary divisor words "
            << words.size() << "; labelled order rows "
            << supports.size() * words.size() << '\n';
  std::cout << "literal state words/support " << state_words_per_support
            << "; raw labelled state contexts "
            << supports.size() * state_words_per_support << '\n';
  std::cout << "grammar: at least two orders carrying 8 and at least two "
               "carrying 3; all six leave-one-out lcms audited\n";
  std::cout << "literal CRT-base FNV64 " << hex64(table.base_digest) << '\n';
  std::cout << "literal owner-mask FNV64 " << hex64(table.mask_digest) << '\n';
  std::cout << "cubic three-class nerve FNV64 "
            << hex64(cubic_nerve_digest) << '\n';
  std::cout << "hereditary weighted-grammar FNV64 "
            << hex64(grammar_digest.state) << '\n';
  std::cout << "scalar contexts " << scalar_bank.size() << " on "
            << scalar_supports.size() << " supports across "
            << multiplicity_histogram.size() << " multiplicity profiles\n";
  std::cout << "scalar-bank FNV64 " << hex64(scalar_digest.state) << '\n';
  std::cout << "capacity-bank FNV64 " << hex64(capacity_digest.state) << '\n';
  std::cout << "distinct scalar capacity vectors " << capacity_vectors.size()
            << '\n';
  std::cout << "all-support scalar-context histogram "
            << histogram(all_support_context_histogram) << '\n';
  std::cout << "owner-local rows " << 6 * scalar_bank.size()
            << "; feasible rows " << feasible_owner_rows << '\n';
  std::cout << "feasible-owner/context histogram "
            << histogram(feasible_contexts) << '\n';
  std::cout << "feasible-owner-mask histogram "
            << feasible_mask_histogram(feasible_masks) << '\n';
  std::cout << "maximum reachable sheet-union histogram "
            << histogram(maximum_union) << '\n';
  std::cout << "distinct owner maximum vectors "
            << owner_maximum_vectors.size() << '\n';
  std::cout << "unique exact normalized owner problems "
            << local_cache.size() << '\n';
  std::cout << "largest reachable bank " << largest_reachable_bank
            << "; largest intermediate layer " << largest_intermediate_layer
            << "; total reachable masks " << total_reachable_masks
            << "; bank-size bins " << reachable_counts.size() << '\n';
  std::cout << "canonical owner-bank FNV64 "
            << hex64(owner_bank_digest.state) << '\n';
  std::cout << "reachable-count histogram " << histogram(reachable_counts)
            << '\n';
  std::cout << "maximum-mask-count histogram "
            << histogram(maximum_mask_counts) << '\n';
  std::cout << "largest-layer histogram " << histogram(largest_layers) << '\n';
  std::cout << "tournament carrier owner obligations; pair observable "
               "(feasible,max-union,capacity,reachable-count,maximum-mask-"
               "count); lexicographic switch; coordinate tie path\n";
  std::cout << "tournament fingerprints: all " << scalar_bank.size()
            << " transitive; score word 0,1,2,3,4,5; triangles 0; SCCs 6; "
               "Hamiltonian paths 1\n";
  std::cout << "tournament tie-edge histogram " << histogram(tournament_ties)
            << '\n';
  std::cout << "tournament edge-flip histogram "
            << histogram(tournament_flips) << '\n';
  std::cout << "tournament aggregate scores "
            << histogram(tournament_scores) << '\n';
  std::cout << "tournament-bank FNV64 " << hex64(tournament_digest.state)
            << '\n';
  std::cout << "four-feasible exceptional rows " << four_feasible_rows.size()
            << "; multiplication-orbit sizes "
            << histogram(four_orbit_sizes) << '\n';
  std::cout << "four-feasible cubic type counts (high-order,ratio-distance)";
  for (const auto &[type, count] : four_types)
    std::cout << " (" << type.first << ',' << type.second << "):" << count;
  std::cout << "; exact-row FNV64 " << hex64(four_row_digest.state) << '\n';
  std::cout << "four-feasible order profiles n1,n2,n3,n4,n6,n8,n12,n24";
  for (const auto &[profile, count] : four_order_profiles)
    std::cout << ' ' << comma_values(profile) << ':' << count;
  std::cout << '\n';
  std::cout << "missing-owner unordered ratio-exponent distances "
            << histogram(missing_ratio_distances) << '\n';
  std::cout << "missing-owner graph edges";
  for (const auto &[edge, count] : missing_edges)
    std::cout << ' ' << edge.first << '-' << edge.second << ':' << count;
  std::cout << '\n';
  std::cout << "missing-owner graph degrees";
  for (int label = 1; label < P; ++label)
    std::cout << ' ' << label << ':' << missing_degrees[label];
  std::cout << '\n';
  std::cout << "cubic classification: for primitive root 2, C0={1,5,8,12}, "
               "C1={2,3,10,11}, C2={4,6,7,9}; every exceptional support "
               "is a*C2 union E with E a two-subset of a*C0 and a*C1 "
               "absent; the feasible owners are exactly a*C2\n";
  std::cout << "missing-owner Cayley multigraph: exponent steps +/-3 have "
               "multiplicity two (high orders 8 and 24), while step 6 has "
               "multiplicity one (high order 8); every label has degree 5\n";
  std::cout << "three-class nerve proof: at a low owner, D3 ratios in C0 "
               "cover sheet classes 2,0,1 mod 3.  At a high owner, D3 ratios "
               "in C2 have two empty fibres and two fibres confined to one "
               "of classes 0,1; if those classes differ they cover 16, and "
               "the two high-order fibres add at most 3 points in the "
               "remaining class.  Therefore the high-owner maximum is 19 "
               "(attained), versus low-owner maximum 24\n";
  for (std::size_t index = 0; index < four_feasible_rows.size(); ++index) {
    const FourFeasibleRow &row = four_feasible_rows[index];
    std::array<uint8_t, 6> actual_orders{};
    std::array<uint8_t, 2> missing{};
    int next_missing = 0;
    for (int coordinate = 0; coordinate < 6; ++coordinate) {
      actual_orders[coordinate] = ORDERS[row.word[coordinate]];
      if (!((row.feasible_mask >> coordinate) & 1U))
        missing[next_missing++] = row.support[coordinate];
    }
    std::cout << "four-row " << index << " support="
              << comma_values(row.support) << " orders="
              << comma_values(actual_orders) << " missing="
              << comma_values(missing) << " capacities="
              << comma_values(row.capacities) << " maxima="
              << comma_values(row.maxima) << " feasible-mask=0x" << std::hex
              << static_cast<int>(row.feasible_mask) << std::dec << '\n';
  }
  std::cout << "preserved: exact FULL-mask feasibility, maximum deficit, "
               "absolute 24-sheet masks, all local unit choices, and every "
               "labelled owner obligation\n";
  std::cout << "lost: an owner projection forgets simultaneous cross-owner "
               "unit gluing; the tournament also forgets exact masks, sheet "
               "identities, threshold magnitude, and witness incidence\n";
  std::cout << "challenged vertices: providers, sheets, wall events, residue "
               "states, gaps, boundaries, cover arcs, Fourier modes, and "
               "matroid circuits all need incidence sidecars; owner proof "
               "obligations are terminal-faithful here\n";
  std::cout << "verdict: independent literal replay agrees; every scalar row "
               "has at least two impossible owners, so the scale-24 common "
               "H6 face is empty\n";
}
