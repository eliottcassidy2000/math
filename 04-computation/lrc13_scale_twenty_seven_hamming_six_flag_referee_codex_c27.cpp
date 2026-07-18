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
#include <unordered_set>
#include <utility>
#include <vector>

// Independent literal-CRT and 3-adic flag referee for the primitive proper
// AP-centred common-scale-twenty-seven Hamming-six face.
//
// The program was designed without reading the preliminary scale-27 primary.
// It reconstructs the divisor grammar, all scalar rows, and every owner-local
// union bank from definitions.  Its proof-facing quotient is smaller: write
// the 27 sheets as nine residue fibres modulo 9, each containing three sheets.
// Orders 1, 3, and 9 fill whole fibres, while every order-27 mask meets a fibre
// in at most one point.  A saturated three-bit hit count is therefore a
// necessary cover observable.  On the complete scalar survivor bank this flag
// quotient is in fact equivalent to exact owner-local cover feasibility.
//
// The decisive structural fact is that every order-27 owner has flag score at
// most 26.  Hereditary lcm supplies at least two such owners in every row, so a
// global unit word is impossible.  The exact union DP is an independent audit,
// not the logical carrier of this obstruction.

namespace {

constexpr int P = 13;
constexpr int C = 27;
constexpr uint32_t FULL = (uint32_t{1} << C) - 1;
constexpr std::array<int, 4> ORDERS{1, 3, 9, 27};

using Support = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;  // indices into ORDERS
using Capacities = std::array<uint8_t, 6>;
using Multiplicity = std::array<uint8_t, 4>;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

struct Fnv64 {
  uint64_t state = 14695981039346656037ULL;

  void byte(uint8_t value) {
    state ^= value;
    state *= 1099511628211ULL;
  }
  void u16(uint16_t value) {
    byte(static_cast<uint8_t>(value));
    byte(static_cast<uint8_t>(value >> 8));
  }
  void u32(uint32_t value) {
    for (int shift = 0; shift < 32; shift += 8)
      byte(static_cast<uint8_t>(value >> shift));
  }
  void u64(uint64_t value) {
    for (int shift = 0; shift < 64; shift += 8)
      byte(static_cast<uint8_t>(value >> shift));
  }
  void text(const std::string &value) {
    u64(value.size());
    for (unsigned char character : value) byte(character);
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

std::vector<uint8_t> units_mod_order(int order) {
  if (order == 1) return {0};
  std::vector<uint8_t> result;
  for (int unit = 1; unit < order; ++unit)
    if (std::gcd(unit, order) == 1)
      result.push_back(static_cast<uint8_t>(unit));
  return result;
}

int literal_crt_base(int label, int order, int unit) {
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P && value % order == unit % order)
      return value;
  fail("bounded CRT search found no representative");
}

uint32_t literal_mask(int label, int order, int unit, int owner) {
  const int base = literal_crt_base(label, order, unit);
  const int owner_inverse = inverse_mod_13(owner);
  uint32_t mask = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value =
        centered(base * (owner_inverse + P * sheet), P * order);
    if (-order < value && value <= order)
      mask |= uint32_t{1} << sheet;
  }
  return mask;
}

uint32_t rotate27(uint32_t mask, int amount) {
  amount %= C;
  if (amount < 0) amount += C;
  if (amount == 0) return mask & FULL;
  return ((mask << amount) | (mask >> (C - amount))) & FULL;
}

uint16_t occupied_fibres(uint32_t mask) {
  uint16_t signature = 0;
  for (int residue = 0; residue < 9; ++residue)
    for (int sheet = residue; sheet < C; sheet += 9)
      if ((mask >> sheet) & 1) {
        signature |= uint16_t{1} << residue;
        break;
      }
  return signature;
}

uint16_t complete_fibres(uint32_t mask) {
  uint16_t signature = 0;
  for (int residue = 0; residue < 9; ++residue) {
    int count = 0;
    for (int sheet = residue; sheet < C; sheet += 9)
      count += static_cast<int>((mask >> sheet) & 1);
    if (count == 3) signature |= uint16_t{1} << residue;
  }
  return signature;
}

uint32_t expand_complete_fibres(uint16_t signature) {
  uint32_t mask = 0;
  for (int residue = 0; residue < 9; ++residue)
    if ((signature >> residue) & 1)
      for (int sheet = residue; sheet < C; sheet += 9)
        mask |= uint32_t{1} << sheet;
  return mask;
}

struct Tables {
  std::array<std::vector<uint8_t>, 4> units;
  // Owner-one normalized fibres [provider/owner ratio][order index].
  std::array<std::array<std::vector<uint32_t>, 4>, P> fibres;
  std::array<std::array<uint8_t, 4>, P> cards{};
  std::array<std::array<uint8_t, 4>, P> signature_counts{};
  uint64_t base_digest = 0;
  uint64_t mask_digest = 0;
  uint64_t signature_digest = 0;
};

Tables build_tables() {
  Tables table;
  Fnv64 bases;
  Fnv64 masks;
  Fnv64 signatures;
  const std::array<std::size_t, 4> expected_phi{1, 2, 6, 18};

  for (int order_index = 0; order_index < 4; ++order_index) {
    table.units[order_index] = units_mod_order(ORDERS[order_index]);
    require(table.units[order_index].size() == expected_phi[order_index],
            "exact-order residue census mismatch");
  }

  for (int ratio = 1; ratio < P; ++ratio)
    for (int order_index = 0; order_index < 4; ++order_index) {
      const int order = ORDERS[order_index];
      std::set<uint16_t> signature_bank;
      int common_card = -1;
      for (uint8_t unit : table.units[order_index]) {
        const uint32_t mask = literal_mask(ratio, order, unit, 1);
        table.fibres[ratio][order_index].push_back(mask);
        const int card = std::popcount(mask);
        if (common_card < 0) common_card = card;
        require(card == common_card,
                "normalized mask cardinality depends on exact-order residue");

        if (order < C) {
          const uint16_t signature = complete_fibres(mask);
          require(expand_complete_fibres(signature) == mask,
                  "lower-order mask is not a union of mod-nine fibres");
          signature_bank.insert(signature);
        } else {
          const uint16_t signature = occupied_fibres(mask);
          for (int residue = 0; residue < 9; ++residue) {
            int count = 0;
            for (int sheet = residue; sheet < C; sheet += 9)
              count += static_cast<int>((mask >> sheet) & 1);
            require(count <= 1,
                    "order-twenty-seven mask repeats a mod-nine fibre");
          }
          signature_bank.insert(signature);
        }
      }
      auto &fibre = table.fibres[ratio][order_index];
      std::sort(fibre.begin(), fibre.end());
      fibre.erase(std::unique(fibre.begin(), fibre.end()), fibre.end());
      table.cards[ratio][order_index] = static_cast<uint8_t>(common_card);
      table.signature_counts[ratio][order_index] =
          static_cast<uint8_t>(signature_bank.size());
      signatures.byte(static_cast<uint8_t>(ratio));
      signatures.byte(static_cast<uint8_t>(order));
      for (uint16_t signature : signature_bank) signatures.u16(signature);
    }

  // Audit the ratio reduction and common cyclic sheet gauge on all literal
  // label/order/unit/owner masks.  No multiplication orbit is skipped.
  for (int label = 1; label < P; ++label)
    for (int order_index = 0; order_index < 4; ++order_index) {
      const int order = ORDERS[order_index];
      for (uint8_t unit : table.units[order_index]) {
        bases.byte(static_cast<uint8_t>(label));
        bases.byte(static_cast<uint8_t>(order));
        bases.byte(unit);
        bases.u16(static_cast<uint16_t>(literal_crt_base(label, order, unit)));
      }
      for (int owner = 1; owner < P; ++owner) {
        const int ratio = label * inverse_mod_13(owner) % P;
        int shift = -1;
        for (int candidate = 0; candidate < C; ++candidate)
          if ((inverse_mod_13(owner) + P * candidate) % C == 1) {
            shift = candidate;
            break;
          }
        require(shift >= 0, "owner sheet-gauge rotation does not exist");

        int common_card = -1;
        for (uint8_t unit : table.units[order_index]) {
          const uint32_t actual = literal_mask(label, order, unit, owner);
          const uint32_t normalized = literal_mask(ratio, order, unit, 1);
          require(actual == rotate27(normalized, shift),
                  "literal owner mask violates the cyclic ratio gauge");
          const int card = std::popcount(actual);
          if (common_card < 0) common_card = card;
          require(card == common_card,
                  "owner mask cardinality depends on exact-order residue");

          const int base = literal_crt_base(ratio, order, unit);
          int one_period = 0;
          for (int sheet = 0; sheet < order; ++sheet) {
            const int value = centered(base * (1 + P * sheet), P * order);
            one_period += static_cast<int>(-order < value && value <= order);
          }
          require(card == (C / order) * one_period,
                  "literal full scan disagrees with one-period cardinality");

          masks.byte(static_cast<uint8_t>(label));
          masks.byte(static_cast<uint8_t>(order));
          masks.byte(unit);
          masks.byte(static_cast<uint8_t>(owner));
          masks.u32(actual);
        }
        require(common_card == table.cards[ratio][order_index],
                "ratio cardinality table disagrees with literal owner mask");
      }
    }

  const std::array<std::array<uint8_t, 12>, 4> expected_cards{{
      {27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {9, 0, 0, 9, 9, 0, 0, 9, 9, 0, 0, 0},
      {6, 6, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3},
      {5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4},
  }};
  const std::array<std::array<uint8_t, 12>, 4> expected_signatures{{
      {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      {1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1},
      {6, 6, 6, 2, 6, 6, 6, 6, 2, 6, 6, 6},
      {6, 3, 6, 6, 6, 6, 6, 6, 6, 6, 3, 6},
  }};
  for (int order_index = 0; order_index < 4; ++order_index)
    for (int ratio = 1; ratio < P; ++ratio) {
      require(table.cards[ratio][order_index] ==
                  expected_cards[order_index][ratio - 1],
              "ratio-cardinality row changed");
      require(table.signature_counts[ratio][order_index] ==
                  expected_signatures[order_index][ratio - 1],
              "mod-nine signature count changed");
    }

  table.base_digest = bases.state;
  table.mask_digest = masks.state;
  table.signature_digest = signatures.state;
  return table;
}

bool hereditary(const OrderWord &word) {
  for (int omitted = 0; omitted < 6; ++omitted) {
    int value = 1;
    for (int index = 0; index < 6; ++index)
      if (index != omitted) value = std::lcm(value, ORDERS[word[index]]);
    if (value != C) return false;
  }
  return true;
}

struct Grammar {
  std::vector<OrderWord> words;
  uint64_t weighted_states = 0;
  uint64_t digest = 0;
};

Grammar build_grammar(const Tables &table) {
  Grammar grammar;
  Fnv64 digest;
  for (int encoded = 0; encoded < 4096; ++encoded) {
    int remainder = encoded;
    OrderWord word{};
    int count_27 = 0;
    uint64_t fibre_size = 1;
    for (int index = 0; index < 6; ++index) {
      word[index] = static_cast<uint8_t>(remainder % 4);
      remainder /= 4;
      count_27 += static_cast<int>(word[index] == 3);
      fibre_size *= table.units[word[index]].size();
    }
    const bool by_lcm = hereditary(word);
    const bool by_top_power = count_27 >= 2;
    require(by_lcm == by_top_power,
            "hereditary lcm is not equivalent to two top-power providers");
    if (!by_lcm) continue;
    grammar.words.push_back(word);
    grammar.weighted_states += fibre_size;
    for (uint8_t order_index : word)
      digest.byte(static_cast<uint8_t>(ORDERS[order_index]));
    digest.u64(fibre_size);
  }
  require(grammar.words.size() == 1909, "hereditary word census mismatch");
  require(grammar.weighted_states == 380511756ULL,
          "weighted exact-order word census mismatch");
  grammar.digest = digest.state;
  return grammar;
}

struct ScalarRow {
  Support support{};
  OrderWord word{};
  Capacities capacities{};
};

Multiplicity multiplicity(const OrderWord &word) {
  Multiplicity result{};
  for (uint8_t order_index : word) ++result[order_index];
  return result;
}

std::string row_key(const Support &support, const OrderWord &word) {
  std::string result;
  result.reserve(12);
  for (int index = 0; index < 6; ++index) {
    result.push_back(static_cast<char>(support[index]));
    result.push_back(static_cast<char>(word[index]));
  }
  return result;
}

struct ScalarBank {
  std::vector<ScalarRow> rows;
  std::map<Multiplicity, int> profile_histogram;
  std::map<int, int> capacity_histogram;
  int support_count = 0;
  std::map<int, int> orbit_size_histogram;
  uint64_t digest = 0;
};

ScalarBank build_scalar_bank(const Tables &table, const Grammar &grammar) {
  ScalarBank bank;
  Fnv64 digest;
  std::set<Support> surviving_supports;

  for (int a = 1; a <= 7; ++a)
    for (int b = a + 1; b <= 8; ++b)
      for (int c = b + 1; c <= 9; ++c)
        for (int d = c + 1; d <= 10; ++d)
          for (int e = d + 1; e <= 11; ++e)
            for (int f = e + 1; f <= 12; ++f) {
              const Support support{static_cast<uint8_t>(a),
                                    static_cast<uint8_t>(b),
                                    static_cast<uint8_t>(c),
                                    static_cast<uint8_t>(d),
                                    static_cast<uint8_t>(e),
                                    static_cast<uint8_t>(f)};
              for (const OrderWord &word : grammar.words) {
                ScalarRow row{support, word, {}};
                bool survives = true;
                for (int owner_index = 0; owner_index < 6; ++owner_index) {
                  const int inverse = inverse_mod_13(support[owner_index]);
                  int capacity = 0;
                  for (int provider = 0; provider < 6; ++provider) {
                    const int ratio = support[provider] * inverse % P;
                    capacity += table.cards[ratio][word[provider]];
                  }
                  row.capacities[owner_index] =
                      static_cast<uint8_t>(capacity);
                  survives = survives && capacity >= C;
                }
                if (!survives) continue;
                bank.rows.push_back(row);
                surviving_supports.insert(support);
                ++bank.profile_histogram[multiplicity(word)];
                for (uint8_t capacity : row.capacities)
                  ++bank.capacity_histogram[capacity];
                for (uint8_t label : support) digest.byte(label);
                for (uint8_t order_index : word)
                  digest.byte(static_cast<uint8_t>(ORDERS[order_index]));
                for (uint8_t capacity : row.capacities) digest.byte(capacity);
              }
            }

  require(bank.rows.size() == 450, "scalar survivor census mismatch");
  require(surviving_supports.size() == 84,
          "scalar surviving-support census mismatch");
  const std::map<Multiplicity, int> expected_profiles{
      {{0, 0, 3, 3}, 12}, {{0, 0, 4, 2}, 18},
      {{0, 2, 1, 3}, 60}, {{0, 2, 2, 2}, 294},
      {{0, 3, 0, 3}, 12}, {{0, 3, 1, 2}, 36},
      {{0, 4, 0, 2}, 18},
  };
  require(bank.profile_histogram == expected_profiles,
          "scalar multiplicity classification changed");

  // Multiplication covariance is audited after the full labelled scan; it is
  // never used to remove a support/order row.
  std::set<std::string> row_keys;
  for (const ScalarRow &row : bank.rows)
    row_keys.insert(row_key(row.support, row.word));
  std::set<std::string> remaining = row_keys;
  while (!remaining.empty()) {
    const std::string seed = *remaining.begin();
    std::vector<std::pair<int, int>> seed_pairs;
    for (int index = 0; index < 6; ++index)
      seed_pairs.push_back({static_cast<uint8_t>(seed[2 * index]),
                            static_cast<uint8_t>(seed[2 * index + 1])});
    std::set<std::string> orbit;
    for (int multiplier = 1; multiplier < P; ++multiplier) {
      auto transformed = seed_pairs;
      for (auto &[label, order_index] : transformed)
        label = multiplier * label % P;
      std::sort(transformed.begin(), transformed.end());
      std::string key;
      for (const auto &[label, order_index] : transformed) {
        key.push_back(static_cast<char>(label));
        key.push_back(static_cast<char>(order_index));
      }
      require(row_keys.contains(key),
              "multiplication orbit left the scalar survivor bank");
      orbit.insert(key);
    }
    ++bank.orbit_size_histogram[static_cast<int>(orbit.size())];
    for (const std::string &key : orbit) remaining.erase(key);
  }
  require(bank.orbit_size_histogram == std::map<int, int>{{6, 3}, {12, 36}},
          "multiplication-orbit census mismatch");

  bank.support_count = static_cast<int>(surviving_supports.size());
  bank.digest = digest.state;
  return bank;
}

std::string normalized_owner_key(const ScalarRow &row, int owner_index) {
  const int inverse = inverse_mod_13(row.support[owner_index]);
  std::vector<std::pair<int, int>> providers;
  for (int index = 0; index < 6; ++index)
    providers.push_back(
        {row.support[index] * inverse % P, row.word[index]});
  std::sort(providers.begin(), providers.end());
  std::string key;
  key.reserve(12);
  for (const auto &[ratio, order_index] : providers) {
    key.push_back(static_cast<char>(ratio));
    key.push_back(static_cast<char>(order_index));
  }
  return key;
}

struct ExactSummary {
  int maximum = 0;
  bool feasible = false;
  uint64_t bank_size = 0;
  uint64_t bank_digest = 0;
};

ExactSummary exact_owner_bank(const std::string &key, const Tables &table) {
  std::unordered_set<uint32_t> reachable{0};
  for (int offset = 0; offset < 12; offset += 2) {
    const int ratio = static_cast<uint8_t>(key[offset]);
    const int order_index = static_cast<uint8_t>(key[offset + 1]);
    std::unordered_set<uint32_t> next;
    next.reserve(reachable.size() * table.fibres[ratio][order_index].size());
    for (uint32_t partial : reachable)
      for (uint32_t option : table.fibres[ratio][order_index])
        next.insert(partial | option);
    reachable = std::move(next);
  }

  std::vector<uint32_t> sorted(reachable.begin(), reachable.end());
  std::sort(sorted.begin(), sorted.end());
  ExactSummary summary;
  summary.bank_size = sorted.size();
  Fnv64 digest;
  for (uint32_t mask : sorted) {
    summary.maximum = std::max(summary.maximum, std::popcount(mask));
    summary.feasible = summary.feasible || mask == FULL;
    digest.u32(mask);
  }
  summary.bank_digest = digest.state;
  return summary;
}

struct FlagState {
  uint16_t complete = 0;  // fibres filled by an order below 27
  uint32_t hits = 0;      // nine saturated base-four counters for D27 hits

  bool operator<(const FlagState &other) const {
    return std::tie(complete, hits) < std::tie(other.complete, other.hits);
  }
};

uint32_t add_saturated_hits(uint32_t packed, uint16_t signature) {
  for (int fibre = 0; fibre < 9; ++fibre)
    if ((signature >> fibre) & 1) {
      const int count = static_cast<int>((packed >> (2 * fibre)) & 3);
      if (count < 3) packed += uint32_t{1} << (2 * fibre);
    }
  return packed;
}

struct FlagSummary {
  int maximum = 0;
  bool feasible = false;
  uint64_t state_count = 0;
  uint64_t state_digest = 0;
};

FlagSummary flag_owner_bank(const std::string &key, const Tables &table) {
  std::set<FlagState> states{{}};
  for (int offset = 0; offset < 12; offset += 2) {
    const int ratio = static_cast<uint8_t>(key[offset]);
    const int order_index = static_cast<uint8_t>(key[offset + 1]);
    std::set<uint16_t> options;
    for (uint32_t mask : table.fibres[ratio][order_index])
      options.insert(order_index < 3 ? complete_fibres(mask)
                                     : occupied_fibres(mask));

    std::set<FlagState> next;
    for (const FlagState &state : states)
      for (uint16_t option : options) {
        FlagState successor = state;
        if (order_index < 3)
          successor.complete |= option;
        else
          successor.hits = add_saturated_hits(successor.hits, option);
        next.insert(successor);
      }
    states = std::move(next);
  }

  FlagSummary summary;
  summary.state_count = states.size();
  Fnv64 digest;
  for (const FlagState &state : states) {
    int score = 0;
    for (int fibre = 0; fibre < 9; ++fibre) {
      const int hits = static_cast<int>((state.hits >> (2 * fibre)) & 3);
      score += ((state.complete >> fibre) & 1) ? 3 : hits;
    }
    summary.maximum = std::max(summary.maximum, score);
    summary.feasible = summary.feasible || score == C;
    digest.u16(state.complete);
    digest.u32(state.hits);
  }
  summary.state_digest = digest.state;
  return summary;
}

template <typename Key>
int winner_from_keys(const Key &left_key, const Key &right_key, int left,
                     int right, int *ties) {
  if (left_key == right_key) {
    ++*ties;
    return left;  // coordinate-order tie Hamiltonian path
  }
  return right_key > left_key ? right : left;
}

struct TournamentFingerprint {
  int ties = 0;
  int natural_flips = 0;
  std::array<uint8_t, 6> out{};
};

template <typename Key>
TournamentFingerprint completed_tournament(const std::array<Key, 6> &keys) {
  TournamentFingerprint result;
  for (int left = 0; left < 6; ++left)
    for (int right = left + 1; right < 6; ++right) {
      const int winner = winner_from_keys(keys[left], keys[right], left, right,
                                          &result.ties);
      const int loser = left + right - winner;
      result.out[winner] |= uint8_t{1} << loser;
      result.natural_flips += static_cast<int>(winner == right);
    }
  return result;
}

void require_transitive_fingerprint(const TournamentFingerprint &fingerprint) {
  std::array<int, 6> scores{};
  for (int vertex = 0; vertex < 6; ++vertex)
    scores[vertex] = std::popcount(fingerprint.out[vertex]);
  std::sort(scores.begin(), scores.end());
  require(scores == std::array<int, 6>{0, 1, 2, 3, 4, 5},
          "completed owner tournament has nontransitive score sequence");

  int triangles = 0;
  for (int a = 0; a < 4; ++a)
    for (int b = a + 1; b < 5; ++b)
      for (int c = b + 1; c < 6; ++c)
        triangles += static_cast<int>(
            (((fingerprint.out[a] >> b) & 1) &&
             ((fingerprint.out[b] >> c) & 1) &&
             ((fingerprint.out[c] >> a) & 1)) ||
            (((fingerprint.out[a] >> c) & 1) &&
             ((fingerprint.out[c] >> b) & 1) &&
             ((fingerprint.out[b] >> a) & 1)));
  require(triangles == 0, "completed owner tournament has a directed triangle");

  std::array<std::array<uint64_t, 6>, 64> paths{};
  for (int vertex = 0; vertex < 6; ++vertex)
    paths[1 << vertex][vertex] = 1;
  for (int subset = 1; subset < 64; ++subset)
    for (int last = 0; last < 6; ++last)
      if ((subset >> last) & 1) {
        const int previous_subset = subset ^ (1 << last);
        for (int previous = 0; previous < 6; ++previous)
          if (((previous_subset >> previous) & 1) &&
              ((fingerprint.out[previous] >> last) & 1))
            paths[subset][last] += paths[previous_subset][previous];
      }
  uint64_t hamiltonian_paths = 0;
  for (int last = 0; last < 6; ++last)
    hamiltonian_paths += paths[63][last];
  require(hamiltonian_paths == 1,
          "completed owner tournament has more than one Hamiltonian path");
  // A transitive tournament has six singleton SCCs; the score/triangle/path
  // checks above are intentionally redundant fingerprints of that fact.
}

struct OwnerAudit {
  std::map<int, ExactSummary> exact_by_position;
  std::map<int, FlagSummary> flag_by_position;
};

template <typename Key>
std::string format_counter(const std::map<Key, int> &counter) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, count] : counter) {
    if (!first) out << ' ';
    first = false;
    out << key << ':' << count;
  }
  return out.str();
}

std::string format_profile(const Multiplicity &profile) {
  std::ostringstream out;
  out << '(' << static_cast<int>(profile[0]) << ','
      << static_cast<int>(profile[1]) << ','
      << static_cast<int>(profile[2]) << ','
      << static_cast<int>(profile[3]) << ')';
  return out.str();
}

struct Audit {
  std::map<std::string, ExactSummary> exact_memo;
  std::map<std::string, FlagSummary> flag_memo;
  std::map<int, int> feasible_owner_histogram;
  std::map<std::pair<int, int>, int> exact_maximum_by_order;
  std::map<std::pair<int, int>, int> flag_maximum_by_order;
  std::map<int, int> exact_tie_histogram;
  std::map<int, int> flag_tie_histogram;
  std::map<int, int> gauge_flip_histogram;
  std::map<int, int> exact_natural_flip_histogram;
  std::map<int, int> flag_natural_flip_histogram;
  std::map<int, std::set<std::string>> keys_by_owner_order;
  std::map<Multiplicity, std::map<int, int>> profile_feasible_histogram;
  std::map<Multiplicity, int> d27_profile_flag_bound;
  uint64_t reachable_total = 0;
  uint64_t exact_digest = 0;
  uint64_t flag_digest = 0;
};

Audit audit_owners(const Tables &table, const ScalarBank &scalar) {
  Audit audit;
  for (const ScalarRow &row : scalar.rows) {
    int feasible_owners = 0;
    std::array<std::tuple<int, int, int, uint64_t>, 6> exact_keys{};
    std::array<std::tuple<int, int, int, uint64_t>, 6> flag_keys{};
    for (int owner_index = 0; owner_index < 6; ++owner_index) {
      const std::string key = normalized_owner_key(row, owner_index);
      const int owner_order = ORDERS[row.word[owner_index]];
      audit.keys_by_owner_order[owner_order].insert(key);
      if (!audit.exact_memo.contains(key))
        audit.exact_memo.emplace(key, exact_owner_bank(key, table));
      if (!audit.flag_memo.contains(key))
        audit.flag_memo.emplace(key, flag_owner_bank(key, table));
      const ExactSummary &exact = audit.exact_memo.at(key);
      const FlagSummary &flag = audit.flag_memo.at(key);
      require(exact.feasible == flag.feasible,
              "mod-nine flag changed exact owner-cover feasibility");
      feasible_owners += static_cast<int>(exact.feasible);
      ++audit.exact_maximum_by_order[{owner_order, exact.maximum}];
      ++audit.flag_maximum_by_order[{owner_order, flag.maximum}];
      audit.reachable_total += exact.bank_size;
      if (owner_order == 27) {
        const Multiplicity profile = multiplicity(row.word);
        audit.d27_profile_flag_bound[profile] =
            std::max(audit.d27_profile_flag_bound[profile], flag.maximum);
        require(!flag.feasible,
                "an order-twenty-seven owner passed the flag obstruction");
      }
      exact_keys[owner_index] =
          {static_cast<int>(exact.feasible), exact.maximum,
           row.capacities[owner_index], exact.bank_size};
      flag_keys[owner_index] =
          {static_cast<int>(flag.feasible), flag.maximum,
           row.capacities[owner_index], flag.state_count};
    }

    ++audit.feasible_owner_histogram[feasible_owners];
    ++audit.profile_feasible_histogram[multiplicity(row.word)][feasible_owners];
    require(feasible_owners <= 4,
            "a scalar survivor has fewer than two impossible owners");

    const TournamentFingerprint exact_tournament =
        completed_tournament(exact_keys);
    const TournamentFingerprint flag_tournament =
        completed_tournament(flag_keys);
    require_transitive_fingerprint(exact_tournament);
    require_transitive_fingerprint(flag_tournament);
    ++audit.exact_tie_histogram[exact_tournament.ties];
    ++audit.flag_tie_histogram[flag_tournament.ties];
    ++audit.exact_natural_flip_histogram[exact_tournament.natural_flips];
    ++audit.flag_natural_flip_histogram[flag_tournament.natural_flips];
    int gauge_flips = 0;
    for (int left = 0; left < 6; ++left)
      for (int right = left + 1; right < 6; ++right) {
        const bool exact_arc = (exact_tournament.out[left] >> right) & 1;
        const bool flag_arc = (flag_tournament.out[left] >> right) & 1;
        gauge_flips += static_cast<int>(exact_arc != flag_arc);
      }
    ++audit.gauge_flip_histogram[gauge_flips];
  }

  require(audit.exact_memo.size() == 225,
          "normalized exact owner-key census mismatch");
  require(audit.flag_memo.size() == 225,
          "normalized flag owner-key census mismatch");
  require(audit.keys_by_owner_order[3].size() == 77 &&
              audit.keys_by_owner_order[9].size() == 66 &&
              audit.keys_by_owner_order[27].size() == 82,
          "owner-key order classification changed");
  require(audit.feasible_owner_histogram ==
              std::map<int, int>{{0, 336}, {1, 96}, {4, 18}},
          "feasible-owner histogram mismatch");

  const std::map<Multiplicity, std::map<int, int>> expected_profile_feasible{
      {{0, 0, 3, 3}, {{0, 12}}},
      {{0, 0, 4, 2}, {{0, 18}}},
      {{0, 2, 1, 3}, {{0, 12}, {1, 48}}},
      {{0, 2, 2, 2}, {{0, 294}}},
      {{0, 3, 0, 3}, {{1, 12}}},
      {{0, 3, 1, 2}, {{1, 36}}},
      {{0, 4, 0, 2}, {{4, 18}}},
  };
  require(audit.profile_feasible_histogram == expected_profile_feasible,
          "profile/feasible-owner classification changed");

  const std::map<Multiplicity, int> expected_d27_bounds{
      {{0, 0, 3, 3}, 26}, {{0, 0, 4, 2}, 22},
      {{0, 2, 1, 3}, 24}, {{0, 2, 2, 2}, 24},
      {{0, 3, 0, 3}, 24}, {{0, 3, 1, 2}, 22},
      {{0, 4, 0, 2}, 22},
  };
  require(audit.d27_profile_flag_bound == expected_d27_bounds,
          "order-twenty-seven profile flag bound changed");

  // Hash every normalized exact bank and every reduced flag state bank in key
  // order, not in unordered-set iteration order.
  Fnv64 exact_digest;
  for (const auto &[key, summary] : audit.exact_memo) {
    exact_digest.text(key);
    exact_digest.byte(static_cast<uint8_t>(summary.maximum));
    exact_digest.byte(static_cast<uint8_t>(summary.feasible));
    exact_digest.u64(summary.bank_size);
    exact_digest.u64(summary.bank_digest);
  }
  Fnv64 flag_digest;
  for (const auto &[key, summary] : audit.flag_memo) {
    flag_digest.text(key);
    flag_digest.byte(static_cast<uint8_t>(summary.maximum));
    flag_digest.byte(static_cast<uint8_t>(summary.feasible));
    flag_digest.u64(summary.state_count);
    flag_digest.u64(summary.state_digest);
  }
  audit.exact_digest = exact_digest.state;
  audit.flag_digest = flag_digest.state;
  return audit;
}

std::string format_order_rows(
    const std::map<std::pair<int, int>, int> &histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, count] : histogram) {
    if (!first) out << ' ';
    first = false;
    out << 'D' << key.first << '/' << key.second << ':' << count;
  }
  return out.str();
}

std::string format_profiles(const std::map<Multiplicity, int> &histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[profile, count] : histogram) {
    if (!first) out << ' ';
    first = false;
    out << format_profile(profile) << ':' << count;
  }
  return out.str();
}

std::string format_profile_feasible(
    const std::map<Multiplicity, std::map<int, int>> &histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[profile, counts] : histogram) {
    if (!first) out << ' ';
    first = false;
    out << format_profile(profile) << "={" << format_counter(counts) << '}';
  }
  return out.str();
}

std::string format_profile_bounds(
    const std::map<Multiplicity, int> &bounds) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[profile, bound] : bounds) {
    if (!first) out << ' ';
    first = false;
    out << format_profile(profile) << "<=" << bound;
  }
  return out.str();
}

template <std::size_t N>
std::string format_array(const std::array<uint8_t, N> &values) {
  std::ostringstream out;
  out << '(';
  for (std::size_t index = 0; index < N; ++index) {
    if (index) out << ',';
    out << static_cast<int>(values[index]);
  }
  out << ')';
  return out.str();
}

}  // namespace

int main() {
  const Tables table = build_tables();
  const Grammar grammar = build_grammar(table);
  const ScalarBank scalar = build_scalar_bank(table, grammar);
  const Audit audit = audit_owners(table, scalar);

  require(grammar.words.size() * 924ULL == 1763916ULL,
          "labelled support/order context census mismatch");
  require(grammar.weighted_states * 924ULL == 351592862544ULL,
          "unquotiented exact-order context census mismatch");

  std::cout << "scale=27 p=13 hamming=6 referee=literal-crt-plus-mod9-flag\n";
  std::cout << "orders=(1,3,9,27) unit-counts=(1,2,6,18)\n";
  std::cout << "hereditary-words=" << grammar.words.size()
            << " weighted-states/support=" << grammar.weighted_states
            << " labelled-support-order-contexts="
            << grammar.words.size() * 924ULL
            << " unquotiented-exact-order-contexts="
            << grammar.weighted_states * 924ULL << '\n';
  for (int order_index = 0; order_index < 4; ++order_index) {
    std::array<uint8_t, 12> cards{};
    std::array<uint8_t, 12> signatures{};
    for (int ratio = 1; ratio < P; ++ratio) {
      cards[ratio - 1] = table.cards[ratio][order_index];
      signatures[ratio - 1] = table.signature_counts[ratio][order_index];
    }
    std::cout << "D" << ORDERS[order_index]
              << " ratio-cardinalities=" << format_array(cards)
              << " mod9-signature-counts=" << format_array(signatures) << '\n';
  }
  std::cout << "scalar-survivors=" << scalar.rows.size()
            << " supports=" << scalar.support_count
            << " multiplicities=" << format_profiles(scalar.profile_histogram)
            << '\n';
  std::cout << "scalar-owner-capacities="
            << format_counter(scalar.capacity_histogram) << '\n';
  std::cout << "multiplication-orbits="
            << format_counter(scalar.orbit_size_histogram) << '\n';
  std::cout << "normalized-owner-keys=225 D3:"
            << audit.keys_by_owner_order.at(3).size() << " D9:"
            << audit.keys_by_owner_order.at(9).size() << " D27:"
            << audit.keys_by_owner_order.at(27).size() << '\n';
  std::cout << "exact-feasible-owners-per-row="
            << format_counter(audit.feasible_owner_histogram) << '\n';
  std::cout << "profile-feasible-owners="
            << format_profile_feasible(audit.profile_feasible_histogram) << '\n';
  std::cout << "exact-maximum-by-owner-order="
            << format_order_rows(audit.exact_maximum_by_order) << '\n';
  std::cout << "flag-maximum-by-owner-order="
            << format_order_rows(audit.flag_maximum_by_order) << '\n';
  std::cout << "D27-profile-flag-bounds="
            << format_profile_bounds(audit.d27_profile_flag_bound) << '\n';
  std::cout << "D27-owner-obstruction=984/984 flag-infeasible; "
               "every-row-has-at-least-two-D27-owners; all-six-owner-row=0\n";
  std::cout << "exact-reachable-mask-incidences=" << audit.reachable_total
            << " flag-equivalence-checks=" << 6 * scalar.rows.size() << '\n';
  std::cout << "owner-tournament exact-observable="
               "(feasible,max-union,scalar-capacity,bank-size) "
               "flag-observable=(feasible,flag-score,scalar-capacity,state-count)\n";
  std::cout << "owner-tournament tie-switch=coordinate-order "
               "scores=(0,1,2,3,4,5) triangles=0 SCCs=(1,1,1,1,1,1) "
               "Hamiltonian-paths=1 rows=450\n";
  std::cout << "owner-tournament exact-ties="
            << format_counter(audit.exact_tie_histogram)
            << " flag-ties=" << format_counter(audit.flag_tie_histogram)
            << " gauge-edge-flips="
            << format_counter(audit.gauge_flip_histogram) << '\n';
  std::cout << "owner-tournament natural-order-flips exact="
            << format_counter(audit.exact_natural_flip_histogram)
            << " flag=" << format_counter(audit.flag_natural_flip_histogram)
            << '\n';
  std::cout << "alternate-carrier=9 mod9 sheet-fibres; "
               "lower-order masks are full-fibre hyperedges and D27 masks "
               "are transversal point-hyperedges; not naturally a tournament\n";
  std::cout << "flag-preserves=owner-cover-feasibility-on-all-2700-survivor-"
               "owner-rows; destroys=within-fibre-offsets,unit-multiplicity,"
               "exact-maximum-union\n";
  std::cout << "challenged-assumption=tournament vertices need not be runners "
               "or owners; the proof carrier is a colored 3-by-9 flag incidence "
               "system\n";
  std::cout << "FNV64 crt-bases=" << hex64(table.base_digest)
            << " literal-masks=" << hex64(table.mask_digest)
            << " flag-signatures=" << hex64(table.signature_digest) << '\n';
  std::cout << "FNV64 grammar=" << hex64(grammar.digest)
            << " scalar=" << hex64(scalar.digest)
            << " exact-owner-banks=" << hex64(audit.exact_digest)
            << " flag-state-banks=" << hex64(audit.flag_digest) << '\n';
}
