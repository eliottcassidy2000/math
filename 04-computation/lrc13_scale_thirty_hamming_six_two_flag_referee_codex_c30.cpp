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

// Independent standard-library referee for the primitive proper AP-centred
// common-scale-30 Hamming-six face (THM-1090).
//
// The implementation was developed from the theorem specification without
// consulting the primary source.  It begins with bounded CRT search, checks
// every literal label/owner gauge, enumerates the squarefree hereditary
// divisor grammar, and scans every labelled support/order context.  Its first
// proof carrier retains the masks whose effective orders divide six as whole
// fibres of Z/30 -> Z/6.  The rows surviving all six owner bounds are then
// tested with the independent anchor family whose orders divide ten, i.e. the
// fibres of Z/30 -> Z/10.  Both tests use the sound anchor/nonanchor upper
// relaxation: nonanchor providers may independently choose their best unit
// outside an exact anchor union.  This only enlarges literal unions.
//
// An immutable-union DP, evaluated in both provider orders, supplies an exact
// sidecar on the complete two-flag residual.  No multiplicative orbit removes
// a support, order word, row, or owner obligation.

namespace {

constexpr int P = 13;
constexpr int C = 30;
constexpr int WIDTH = 6;
constexpr uint32_t FULL = (uint32_t{1} << C) - 1;
constexpr std::array<int, 8> ORDERS{1, 2, 3, 5, 6, 10, 15, 30};
constexpr std::array<int, 3> PRIMES{2, 3, 5};

using Support = std::array<uint8_t, WIDTH>;
using Word = std::array<uint8_t, WIDTH>;  // indices into ORDERS
using Capacities = std::array<uint8_t, WIDTH>;
using Multiplicity = std::array<uint8_t, 8>;
using OwnerKey = std::array<uint8_t, 2 * WIDTH>;  // sorted ratio/order pairs

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
};

std::string hex64(uint64_t value) {
  std::ostringstream out;
  out << std::hex << std::setfill('0') << std::setw(16) << value;
  return out.str();
}

int positive_mod(int value, int modulus) {
  const int residue = value % modulus;
  return residue < 0 ? residue + modulus : residue;
}

int centered(int value, int modulus) {
  const int residue = positive_mod(value, modulus);
  return 2 * residue > modulus ? residue - modulus : residue;
}

int inverse_mod_13(int value) {
  for (int candidate = 1; candidate < P; ++candidate)
    if (value * candidate % P == 1) return candidate;
  fail("inverse modulo thirteen requested for a nonunit");
}

std::vector<uint8_t> exact_units(int order) {
  if (order == 1) return {0};
  std::vector<uint8_t> result;
  for (int unit = 1; unit < order; ++unit)
    if (std::gcd(unit, order) == 1)
      result.push_back(static_cast<uint8_t>(unit));
  return result;
}

int literal_crt_base(int label, int order, int unit) {
  int answer = -1;
  for (int candidate = 0; candidate < P * order; ++candidate) {
    if (candidate % P != order * label % P) continue;
    if (candidate % order != unit % order) continue;
    require(answer < 0, "bounded CRT search found two representatives");
    answer = candidate;
  }
  require(answer >= 0, "bounded CRT search found no representative");
  return answer;
}

uint32_t literal_mask(int label, int order, int unit, int owner) {
  const int base = literal_crt_base(label, order, unit);
  const int inverse = inverse_mod_13(owner);
  uint32_t mask = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    if (-order < value && value <= order)
      mask |= uint32_t{1} << sheet;
  }
  return mask;
}

uint32_t rotate30(uint32_t mask, int amount) {
  amount = positive_mod(amount, C);
  if (amount == 0) return mask & FULL;
  return ((mask << amount) | (mask >> (C - amount))) & FULL;
}

bool periodic_by(uint32_t mask, int period) {
  for (int sheet = 0; sheet < C; ++sheet) {
    const int next = (sheet + period) % C;
    if (((mask >> sheet) & 1U) != ((mask >> next) & 1U)) return false;
  }
  return true;
}

uint16_t whole_fibre_signature(uint32_t mask, int modulus) {
  require(C % modulus == 0, "fibre modulus does not divide the scale");
  uint16_t signature = 0;
  for (int residue = 0; residue < modulus; ++residue) {
    const bool first = (mask >> residue) & 1U;
    for (int sheet = residue; sheet < C; sheet += modulus)
      require(static_cast<bool>((mask >> sheet) & 1U) == first,
              "putative anchor is not a union of whole quotient fibres");
    if (first) signature |= uint16_t{1} << residue;
  }
  return signature;
}

uint32_t expand_fibre_signature(uint16_t signature, int modulus) {
  uint32_t mask = 0;
  for (int residue = 0; residue < modulus; ++residue)
    if ((signature >> residue) & 1U)
      for (int sheet = residue; sheet < C; sheet += modulus)
        mask |= uint32_t{1} << sheet;
  return mask;
}

struct Tables {
  std::array<std::vector<uint8_t>, 8> units;
  // Owner-one, ratio-normalized literal masks.  Unit multiplicity is retained.
  std::array<std::array<std::vector<uint32_t>, 8>, P> fibres;
  std::array<std::array<uint8_t, 8>, P> cards{};
  uint64_t base_digest = 0;
  uint64_t mask_digest = 0;
  uint64_t card_digest = 0;
};

Tables build_tables() {
  Tables table;
  Fnv64 base_digest;
  Fnv64 mask_digest;
  Fnv64 card_digest;
  const std::array<std::size_t, 8> expected_phi{1, 1, 2, 4, 2, 4, 8, 8};

  for (int oi = 0; oi < 8; ++oi) {
    table.units[oi] = exact_units(ORDERS[oi]);
    require(table.units[oi].size() == expected_phi[oi],
            "exact-order unit census mismatch");
  }

  // Build normalized masks by literal bounded CRT search and independently
  // recover each cardinality from one effective-order period.
  for (int ratio = 1; ratio < P; ++ratio)
    for (int oi = 0; oi < 8; ++oi) {
      const int order = ORDERS[oi];
      int common_card = -1;
      for (uint8_t unit : table.units[oi]) {
        const uint32_t mask = literal_mask(ratio, order, unit, 1);
        table.fibres[ratio][oi].push_back(mask);
        const int card = std::popcount(mask);
        if (common_card < 0) common_card = card;
        require(card == common_card,
                "mask cardinality depends on the exact-order unit");
        require(periodic_by(mask, order),
                "normalized mask violates effective-order periodicity");
        if (6 % order == 0) {
          require(periodic_by(mask, 6),
                  "six-anchor mask violates quotient periodicity");
          require(expand_fibre_signature(whole_fibre_signature(mask, 6), 6) ==
                      mask,
                  "six-fibre signature does not reconstruct anchor mask");
        }
        if (10 % order == 0) {
          require(periodic_by(mask, 10),
                  "ten-anchor mask violates quotient periodicity");
          require(expand_fibre_signature(whole_fibre_signature(mask, 10), 10) ==
                      mask,
                  "ten-fibre signature does not reconstruct anchor mask");
        }
      }
      const int target = order * ratio % P;
      int one_period = 0;
      for (int value = -order + 1; value <= order; ++value)
        one_period += static_cast<int>(positive_mod(value, P) == target);
      require(common_card == (C / order) * one_period,
              "single-period interval count disagrees with literal mask");
      table.cards[ratio][oi] = static_cast<uint8_t>(common_card);
      card_digest.byte(static_cast<uint8_t>(ratio));
      card_digest.byte(static_cast<uint8_t>(order));
      card_digest.byte(static_cast<uint8_t>(common_card));
    }

  // Exhaustively establish the common ratio/rotation gauge against literal
  // masks for all labels, owners, orders, and exact-order units.
  for (int label = 1; label < P; ++label)
    for (int oi = 0; oi < 8; ++oi) {
      const int order = ORDERS[oi];
      for (uint8_t unit : table.units[oi]) {
        base_digest.byte(static_cast<uint8_t>(label));
        base_digest.byte(static_cast<uint8_t>(order));
        base_digest.byte(unit);
        base_digest.u16(
            static_cast<uint16_t>(literal_crt_base(label, order, unit)));
      }
      for (int owner = 1; owner < P; ++owner) {
        const int ratio = label * inverse_mod_13(owner) % P;
        int shift = -1;
        for (int candidate = 0; candidate < C; ++candidate)
          if ((inverse_mod_13(owner) + P * candidate) % C == 1) {
            shift = candidate;
            break;
          }
        require(shift >= 0, "common owner sheet rotation was not found");
        int common_card = -1;
        for (std::size_t ui = 0; ui < table.units[oi].size(); ++ui) {
          const uint8_t unit = table.units[oi][ui];
          const uint32_t actual = literal_mask(label, order, unit, owner);
          const uint32_t normalized = table.fibres[ratio][oi][ui];
          require(actual == rotate30(normalized, shift),
                  "literal owner mask violates the ratio gauge");
          const int card = std::popcount(actual);
          if (common_card < 0) common_card = card;
          require(card == common_card,
                  "owner cardinality depends on exact-order unit");
          require(periodic_by(actual, order),
                  "actual mask violates effective-order periodicity");
          if (6 % order == 0)
            require(expand_fibre_signature(
                        whole_fibre_signature(actual, 6), 6) == actual,
                    "actual six-anchor mask is not whole-fibre");
          if (10 % order == 0)
            require(expand_fibre_signature(
                        whole_fibre_signature(actual, 10), 10) == actual,
                    "actual ten-anchor mask is not whole-fibre");
          mask_digest.byte(static_cast<uint8_t>(label));
          mask_digest.byte(static_cast<uint8_t>(order));
          mask_digest.byte(unit);
          mask_digest.byte(static_cast<uint8_t>(owner));
          mask_digest.u32(actual);
        }
        require(common_card == table.cards[ratio][oi],
                "ratio card table disagrees with actual masks");
      }
    }

  table.base_digest = base_digest.state;
  table.mask_digest = mask_digest.state;
  table.card_digest = card_digest.state;
  return table;
}

bool leave_one_out_lcm(const Word &word) {
  for (int omitted = 0; omitted < WIDTH; ++omitted) {
    int value = 1;
    for (int index = 0; index < WIDTH; ++index)
      if (index != omitted) value = std::lcm(value, ORDERS[word[index]]);
    if (value != C) return false;
  }
  return true;
}

uint64_t intersection_count(int prime_subset, const Tables &table,
                            bool weighted) {
  // For every selected prime, count order letters carrying it at most once.
  // A bit in `state` means that prime has already appeared once.
  std::array<uint64_t, 8> states{};
  states[0] = 1;
  for (int coordinate = 0; coordinate < WIDTH; ++coordinate) {
    std::array<uint64_t, 8> next{};
    for (int state = 0; state < 8; ++state) {
      if (states[state] == 0) continue;
      for (int oi = 0; oi < 8; ++oi) {
        int carriers = 0;
        for (int pi = 0; pi < 3; ++pi)
          if ((prime_subset >> pi) & 1)
            carriers |= static_cast<int>(ORDERS[oi] % PRIMES[pi] == 0) << pi;
        if (state & carriers) continue;
        const uint64_t letter_weight =
            weighted ? table.units[oi].size() : uint64_t{1};
        next[state | carriers] += states[state] * letter_weight;
      }
    }
    states = next;
  }
  return std::accumulate(states.begin(), states.end(), uint64_t{0});
}

struct WeightedWord {
  Word word{};
  uint32_t code = 0;  // base-eight code, provider zero most significant
  uint64_t weight = 0;
};

struct Grammar {
  std::vector<WeightedWord> words;
  std::array<uint64_t, 8> intersections{};
  std::array<uint64_t, 8> weighted_intersections{};
  uint64_t weighted_states = 0;
  uint64_t digest = 0;
};

Grammar build_grammar(const Tables &table) {
  Grammar grammar;
  Fnv64 digest;
  uint32_t total_codes = 1;
  for (int i = 0; i < WIDTH; ++i) total_codes *= 8;

  for (uint32_t code = 0; code < total_codes; ++code) {
    uint32_t remainder = code;
    Word word{};
    std::array<int, 3> providers{};
    uint64_t weight = 1;
    for (int index = WIDTH - 1; index >= 0; --index) {
      word[index] = static_cast<uint8_t>(remainder % 8);
      remainder /= 8;
      const int order = ORDERS[word[index]];
      for (int pi = 0; pi < 3; ++pi)
        providers[pi] += static_cast<int>(order % PRIMES[pi] == 0);
      weight *= table.units[word[index]].size();
    }
    const bool valuation_grammar =
        std::all_of(providers.begin(), providers.end(),
                    [](int count) { return count >= 2; });
    require(valuation_grammar == leave_one_out_lcm(word),
            "squarefree valuation grammar disagrees with leave-one-out lcm");
    if (!valuation_grammar) continue;
    grammar.words.push_back({word, code, weight});
    grammar.weighted_states += weight;
    for (uint8_t oi : word)
      digest.byte(static_cast<uint8_t>(ORDERS[oi]));
    digest.u64(weight);
  }

  uint64_t total_unweighted = 1;
  uint64_t total_weighted = 1;
  uint64_t alphabet_weight = 0;
  for (const auto &units : table.units) alphabet_weight += units.size();
  for (int index = 0; index < WIDTH; ++index) {
    total_unweighted *= ORDERS.size();
    total_weighted *= alphabet_weight;
  }
  require(total_codes == total_unweighted,
          "base-eight code census does not cover the order alphabet");
  require(alphabet_weight == C,
          "exact-order unit alphabet does not sum to the common scale");

  int64_t inclusion_unweighted = static_cast<int64_t>(total_unweighted);
  int64_t inclusion_weighted = static_cast<int64_t>(total_weighted);
  for (int subset = 1; subset < 8; ++subset) {
    grammar.intersections[subset] =
        intersection_count(subset, table, false);
    grammar.weighted_intersections[subset] =
        intersection_count(subset, table, true);
    const int sign = std::popcount(static_cast<unsigned>(subset)) % 2 ? -1 : 1;
    inclusion_unweighted +=
        sign * static_cast<int64_t>(grammar.intersections[subset]);
    inclusion_weighted +=
        sign * static_cast<int64_t>(grammar.weighted_intersections[subset]);
  }
  require(inclusion_unweighted == static_cast<int64_t>(grammar.words.size()),
          "unweighted inclusion-exclusion disagrees with grammar");
  require(inclusion_weighted == static_cast<int64_t>(grammar.weighted_states),
          "weighted inclusion-exclusion disagrees with grammar");
  grammar.digest = digest.state;
  return grammar;
}

Multiplicity multiplicity(const Word &word) {
  Multiplicity result{};
  for (uint8_t oi : word) ++result[oi];
  return result;
}

struct ScalarRow {
  Support support{};
  Word word{};
  Capacities capacities{};
};

struct ScalarBank {
  std::vector<ScalarRow> rows;
  std::map<int, uint64_t> feasible_owner_histogram;
  std::map<int, uint64_t> support_survivor_histogram;
  std::map<Multiplicity, uint64_t> multiplicity_histogram;
  std::set<Capacities> capacity_vectors;
  uint64_t literal_survivor_words = 0;
  std::size_t support_count = 0;
  uint64_t digest = 0;
};

uint64_t pack_owner_cards(const std::array<uint8_t, WIDTH> &cards) {
  uint64_t packed = 0;
  for (int owner = 0; owner < WIDTH; ++owner)
    packed |= uint64_t{cards[owner]} << (8 * owner);
  return packed;
}

ScalarBank build_scalar_bank(const Tables &table, const Grammar &grammar) {
  ScalarBank bank;
  Fnv64 digest;
  std::set<Support> surviving_supports;
  uint64_t support_count = 0;

  for (int a = 1; a <= 7; ++a)
    for (int b = a + 1; b <= 8; ++b)
      for (int c = b + 1; c <= 9; ++c)
        for (int d = c + 1; d <= 10; ++d)
          for (int e = d + 1; e <= 11; ++e)
            for (int f = e + 1; f <= 12; ++f) {
              ++support_count;
              const Support support{static_cast<uint8_t>(a),
                                    static_cast<uint8_t>(b),
                                    static_cast<uint8_t>(c),
                                    static_cast<uint8_t>(d),
                                    static_cast<uint8_t>(e),
                                    static_cast<uint8_t>(f)};
              std::array<std::array<uint64_t, 8>, WIDTH> contribution{};
              for (int provider = 0; provider < WIDTH; ++provider)
                for (int oi = 0; oi < 8; ++oi) {
                  std::array<uint8_t, WIDTH> cards{};
                  for (int owner = 0; owner < WIDTH; ++owner) {
                    const int ratio = support[provider] *
                        inverse_mod_13(support[owner]) % P;
                    cards[owner] = table.cards[ratio][oi];
                  }
                  contribution[provider][oi] = pack_owner_cards(cards);
                }

              // Prefix expansion indexes all 8^6 order words in the same
              // provider-major base-eight order used by Grammar::code.  Each
              // owner occupies one byte.  No byte can carry into the next:
              // six literal cards total at most 6*30=180 < 256.
              std::vector<uint64_t> packed_capacities{0};
              for (int provider = 0; provider < WIDTH; ++provider) {
                std::vector<uint64_t> next;
                next.reserve(packed_capacities.size() * 8);
                for (uint64_t partial : packed_capacities)
                  for (int oi = 0; oi < 8; ++oi)
                    next.push_back(partial + contribution[provider][oi]);
                packed_capacities.swap(next);
              }
              require(packed_capacities.size() == 262144,
                      "packed scalar table has the wrong base-eight size");

              int support_survivors = 0;
              for (const WeightedWord &weighted : grammar.words) {
                const uint64_t packed = packed_capacities[weighted.code];
                ScalarRow row{support, weighted.word, {}};
                int feasible = 0;
                for (int owner = 0; owner < WIDTH; ++owner) {
                  const uint8_t capacity =
                      static_cast<uint8_t>(packed >> (8 * owner));
                  row.capacities[owner] = capacity;
                  feasible += static_cast<int>(capacity >= C);
                }
                require((packed >> 48) == 0,
                        "packed scalar capacity overflowed six owner bytes");
                ++bank.feasible_owner_histogram[feasible];
                if (feasible != WIDTH) continue;

                ++support_survivors;
                bank.rows.push_back(row);
                surviving_supports.insert(support);
                ++bank.multiplicity_histogram[multiplicity(weighted.word)];
                bank.capacity_vectors.insert(row.capacities);
                bank.literal_survivor_words += weighted.weight;
                for (uint8_t label : support) digest.byte(label);
                for (uint8_t oi : weighted.word)
                  digest.byte(static_cast<uint8_t>(ORDERS[oi]));
                for (uint8_t capacity : row.capacities) digest.byte(capacity);
              }
              ++bank.support_survivor_histogram[support_survivors];
            }

  require(support_count == 924, "labelled support census is not binom(12,6)");
  uint64_t visited_contexts = 0;
  for (const auto &[feasible, count] : bank.feasible_owner_histogram) {
    require(0 <= feasible && feasible <= WIDTH,
            "invalid scalar feasible-owner count");
    visited_contexts += count;
  }
  require(visited_contexts == support_count * grammar.words.size(),
          "scalar scan did not visit every labelled support/order context");
  for (const ScalarRow &row : bank.rows)
    require(leave_one_out_lcm(row.word),
            "scalar bank contains a nonhereditary order word");
  bank.support_count = surviving_supports.size();
  bank.digest = digest.state;
  return bank;
}

OwnerKey normalized_owner_key(const ScalarRow &row, int owner) {
  std::array<std::pair<uint8_t, uint8_t>, WIDTH> pairs{};
  const int inverse = inverse_mod_13(row.support[owner]);
  for (int provider = 0; provider < WIDTH; ++provider)
    pairs[provider] = {
        static_cast<uint8_t>(row.support[provider] * inverse % P),
        row.word[provider]};
  std::sort(pairs.begin(), pairs.end());
  OwnerKey key{};
  for (int index = 0; index < WIDTH; ++index) {
    key[2 * index] = pairs[index].first;
    key[2 * index + 1] = pairs[index].second;
  }
  return key;
}

std::vector<uint32_t> distinct_options(const Tables &table, int ratio,
                                       int oi) {
  std::vector<uint32_t> result = table.fibres[ratio][oi];
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

std::vector<uint32_t> advance_bank(const std::vector<uint32_t> &bank,
                                   const std::vector<uint32_t> &options) {
  std::vector<uint32_t> next;
  next.reserve(bank.size() * options.size());
  for (uint32_t partial : bank)
    for (uint32_t option : options) next.push_back(partial | option);
  std::sort(next.begin(), next.end());
  next.erase(std::unique(next.begin(), next.end()), next.end());
  return next;
}

struct RelaxedSummary {
  int bound = 0;
  uint64_t anchor_bank_size = 0;
  uint64_t anchor_bank_digest = 0;
};

RelaxedSummary relaxed_owner_bound(const OwnerKey &key, const Tables &table,
                                   int modulus) {
  require(C % modulus == 0, "anchor modulus does not divide scale thirty");
  std::array<std::vector<uint32_t>, WIDTH> options;
  std::array<bool, WIDTH> is_anchor{};
  std::vector<uint32_t> anchors{0};
  for (int layer = 0; layer < WIDTH; ++layer) {
    const int ratio = key[2 * layer];
    const int oi = key[2 * layer + 1];
    const int order = ORDERS[oi];
    options[layer] = distinct_options(table, ratio, oi);
    is_anchor[layer] = modulus % order == 0;
    if (!is_anchor[layer]) continue;
    for (uint32_t mask : options[layer]) {
      require(periodic_by(mask, modulus),
              "relaxation anchor is not quotient-periodic");
      require(expand_fibre_signature(
                  whole_fibre_signature(mask, modulus), modulus) == mask,
              "relaxation anchor is not a whole-fibre union");
    }
    anchors = advance_bank(anchors, options[layer]);
  }

  RelaxedSummary summary;
  summary.anchor_bank_size = anchors.size();
  Fnv64 digest;
  for (uint32_t anchor : anchors) {
    digest.u32(anchor);
    int score = std::popcount(anchor);
    for (int layer = 0; layer < WIDTH; ++layer) {
      if (is_anchor[layer]) continue;
      int best = 0;
      for (uint32_t option : options[layer])
        best = std::max(best, std::popcount(option & (FULL ^ anchor)));
      for (uint32_t option : options[layer])
        require(std::popcount(option & (FULL ^ anchor)) <= best,
                "pointwise nonanchor maximum is not an upper bound");
      score += best;
    }
    summary.bound = std::max(summary.bound, score);
  }
  summary.anchor_bank_digest = digest.state;
  return summary;
}

std::vector<uint32_t> exact_bank(const OwnerKey &key, const Tables &table,
                                 bool reverse) {
  std::vector<uint32_t> bank{0};
  for (int layer = 0; layer < WIDTH; ++layer) {
    const int index = reverse ? WIDTH - 1 - layer : layer;
    bank = advance_bank(
        bank, distinct_options(table, key[2 * index], key[2 * index + 1]));
  }
  return bank;
}

struct ExactSummary {
  bool feasible = false;
  int maximum = 0;
  uint64_t bank_size = 0;
  uint64_t bank_digest = 0;
};

ExactSummary exact_owner_summary(const OwnerKey &key, const Tables &table) {
  const std::vector<uint32_t> forward = exact_bank(key, table, false);
  const std::vector<uint32_t> reverse = exact_bank(key, table, true);
  require(forward == reverse,
          "forward/reverse immutable-union banks disagree");
  ExactSummary summary;
  summary.bank_size = forward.size();
  Fnv64 digest;
  for (uint32_t mask : forward) {
    const int card = std::popcount(mask);
    summary.maximum = std::max(summary.maximum, card);
    summary.feasible = summary.feasible || mask == FULL;
    digest.u32(mask);
  }
  summary.bank_digest = digest.state;
  require(summary.feasible == (summary.maximum == C),
          "exact feasibility disagrees with maximum union size");
  return summary;
}

struct FirstFlagAudit {
  std::map<OwnerKey, RelaxedSummary> memo;
  std::map<int, uint64_t> bound_histogram;
  std::map<uint64_t, uint64_t> anchor_bank_histogram;
  std::map<int, uint64_t> live_owner_histogram;
  std::vector<ScalarRow> residual_rows;
  uint64_t owner_checks = 0;
  uint64_t digest = 0;
};

FirstFlagAudit audit_six_flag(const Tables &table, const ScalarBank &scalar) {
  FirstFlagAudit audit;
  Fnv64 digest;
  for (const ScalarRow &row : scalar.rows) {
    int live = 0;
    for (int owner = 0; owner < WIDTH; ++owner) {
      const OwnerKey key = normalized_owner_key(row, owner);
      auto [iterator, inserted] = audit.memo.try_emplace(key);
      if (inserted)
        iterator->second = relaxed_owner_bound(key, table, 6);
      const RelaxedSummary &local = iterator->second;
      require(local.bound <= row.capacities[owner],
              "six-flag relaxation exceeds scalar capacity");
      live += static_cast<int>(local.bound >= C);
      ++audit.bound_histogram[local.bound];
      ++audit.anchor_bank_histogram[local.anchor_bank_size];
      ++audit.owner_checks;
      for (uint8_t value : key) digest.byte(value);
      digest.byte(static_cast<uint8_t>(local.bound));
      digest.u64(local.anchor_bank_size);
      digest.u64(local.anchor_bank_digest);
    }
    ++audit.live_owner_histogram[live];
    if (live == WIDTH) audit.residual_rows.push_back(row);
  }
  require(audit.owner_checks == WIDTH * scalar.rows.size(),
          "six-flag audit missed a labelled owner obligation");
  require(audit.residual_rows.size() == 120,
          "six-flag all-owner residual is not the specified 120 rows");
  audit.digest = digest.state;
  return audit;
}

struct TournamentFingerprint {
  int ties = 0;
  int flips = 0;
  int triangles = 0;
  int hamiltonian_paths = 0;
  std::array<uint8_t, WIDTH> out{};
};

using TournamentKey = std::tuple<int, int, int, int, uint64_t>;

TournamentFingerprint tournament(const std::array<TournamentKey, WIDTH> &keys) {
  TournamentFingerprint result;
  for (int left = 0; left < WIDTH; ++left)
    for (int right = left + 1; right < WIDTH; ++right) {
      int winner = left;
      if (keys[left] == keys[right]) {
        ++result.ties;
      } else if (keys[right] > keys[left]) {
        winner = right;
        ++result.flips;
      }
      const int loser = left + right - winner;
      result.out[winner] |= uint8_t{1} << loser;
    }

  std::array<int, WIDTH> scores{};
  for (int vertex = 0; vertex < WIDTH; ++vertex)
    scores[vertex] = std::popcount(result.out[vertex]);
  std::sort(scores.begin(), scores.end());
  require(scores == std::array<int, WIDTH>{0, 1, 2, 3, 4, 5},
          "owner tournament has nontransitive score word");

  for (int a = 0; a < 4; ++a)
    for (int b = a + 1; b < 5; ++b)
      for (int c = b + 1; c < WIDTH; ++c)
        result.triangles += static_cast<int>(
            (((result.out[a] >> b) & 1U) &&
             ((result.out[b] >> c) & 1U) &&
             ((result.out[c] >> a) & 1U)) ||
            (((result.out[a] >> c) & 1U) &&
             ((result.out[c] >> b) & 1U) &&
             ((result.out[b] >> a) & 1U)));
  require(result.triangles == 0,
          "owner tournament contains a directed triangle");

  uint64_t paths[64][WIDTH]{};
  for (int vertex = 0; vertex < WIDTH; ++vertex)
    paths[1 << vertex][vertex] = 1;
  for (int subset = 1; subset < 64; ++subset)
    for (int last = 0; last < WIDTH; ++last) {
      if (!((subset >> last) & 1)) continue;
      const int previous_subset = subset ^ (1 << last);
      for (int previous = 0; previous < WIDTH; ++previous)
        if (((previous_subset >> previous) & 1) &&
            ((result.out[previous] >> last) & 1U))
          paths[subset][last] += paths[previous_subset][previous];
    }
  uint64_t path_count = 0;
  for (int last = 0; last < WIDTH; ++last) path_count += paths[63][last];
  result.hamiltonian_paths = static_cast<int>(path_count);
  require(result.hamiltonian_paths == 1,
          "owner tournament Hamiltonian-path count changed");
  return result;
}

struct SecondFlagAudit {
  std::map<OwnerKey, RelaxedSummary> ten_memo;
  std::map<OwnerKey, ExactSummary> exact_memo;
  std::map<int, uint64_t> ten_bound_histogram;
  std::map<uint64_t, uint64_t> ten_anchor_bank_histogram;
  std::map<int, uint64_t> exact_maximum_histogram;
  std::map<uint64_t, uint64_t> exact_bank_histogram;
  std::map<int, uint64_t> ten_minus_exact_histogram;
  std::map<int, uint64_t> tie_histogram;
  std::map<int, uint64_t> flip_histogram;
  uint64_t owner_checks = 0;
  uint64_t terminal_checks = 0;
  uint64_t exact_reachable_total = 0;
  uint64_t largest_exact_bank = 0;
  uint64_t digest = 0;
};

SecondFlagAudit audit_ten_flag_and_exact(
    const Tables &table, const FirstFlagAudit &six_audit) {
  SecondFlagAudit audit;
  Fnv64 digest;
  for (const ScalarRow &row : six_audit.residual_rows) {
    std::array<TournamentKey, WIDTH> keys{};
    for (int owner = 0; owner < WIDTH; ++owner) {
      const OwnerKey key = normalized_owner_key(row, owner);
      auto [ten_iterator, ten_inserted] = audit.ten_memo.try_emplace(key);
      if (ten_inserted)
        ten_iterator->second = relaxed_owner_bound(key, table, 10);
      auto [exact_iterator, exact_inserted] = audit.exact_memo.try_emplace(key);
      if (exact_inserted)
        exact_iterator->second = exact_owner_summary(key, table);
      const RelaxedSummary &ten = ten_iterator->second;
      const ExactSummary &exact = exact_iterator->second;
      const RelaxedSummary &six = six_audit.memo.at(key);

      require(ten.bound <= row.capacities[owner],
              "ten-flag relaxation exceeds scalar capacity");
      require(exact.maximum <= ten.bound,
              "exact union maximum exceeds ten-flag upper bound");
      require(exact.maximum <= six.bound,
              "exact union maximum exceeds six-flag upper bound");
      require(!exact.feasible,
              "two-flag residual contains an exact owner cover");
      ++audit.owner_checks;
      if (ten.bound < C) ++audit.terminal_checks;
      ++audit.ten_bound_histogram[ten.bound];
      ++audit.ten_anchor_bank_histogram[ten.anchor_bank_size];
      ++audit.exact_maximum_histogram[exact.maximum];
      ++audit.exact_bank_histogram[exact.bank_size];
      ++audit.ten_minus_exact_histogram[ten.bound - exact.maximum];
      audit.exact_reachable_total += exact.bank_size;
      audit.largest_exact_bank =
          std::max(audit.largest_exact_bank, exact.bank_size);
      keys[owner] = {ten.bound, six.bound, exact.maximum,
                     row.capacities[owner], exact.bank_size};

      for (uint8_t value : key) digest.byte(value);
      digest.byte(static_cast<uint8_t>(six.bound));
      digest.byte(static_cast<uint8_t>(ten.bound));
      digest.byte(static_cast<uint8_t>(exact.maximum));
      digest.u64(ten.anchor_bank_size);
      digest.u64(exact.bank_size);
      digest.u64(exact.bank_digest);
    }
    const TournamentFingerprint fingerprint = tournament(keys);
    ++audit.tie_histogram[fingerprint.ties];
    ++audit.flip_histogram[fingerprint.flips];
  }
  require(audit.owner_checks == 720,
          "ten-flag audit did not visit all 120*6 labelled obligations");
  require(audit.terminal_checks == audit.owner_checks,
          "a ten-flag residual owner bound reaches scale thirty");
  audit.digest = digest.state;
  return audit;
}

template <typename Key, typename Value>
std::string format_map(const std::map<Key, Value> &histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, value] : histogram) {
    if (!first) out << ' ';
    first = false;
    out << key << ':' << value;
  }
  return out.str();
}

std::string format_array(const std::array<uint64_t, 8> &values) {
  std::ostringstream out;
  for (int subset = 1; subset < 8; ++subset) {
    if (subset > 1) out << ' ';
    out << subset << ':' << values[subset];
  }
  return out.str();
}

std::string format_cards(const Tables &table, int oi) {
  std::ostringstream out;
  out << '(';
  for (int ratio = 1; ratio < P; ++ratio) {
    if (ratio > 1) out << ',';
    out << static_cast<int>(table.cards[ratio][oi]);
  }
  out << ')';
  return out.str();
}

std::string format_multiplicities(
    const std::map<Multiplicity, uint64_t> &histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[profile, count] : histogram) {
    if (!first) out << ' ';
    first = false;
    out << '(';
    for (int index = 0; index < 8; ++index) {
      if (index) out << ',';
      out << static_cast<int>(profile[index]);
    }
    out << "):" << count;
  }
  return out.str();
}

}  // namespace

int main() {
  const Tables table = build_tables();
  const Grammar grammar = build_grammar(table);
  const ScalarBank scalar = build_scalar_bank(table, grammar);
  const FirstFlagAudit six = audit_six_flag(table, scalar);
  const SecondFlagAudit ten = audit_ten_flag_and_exact(table, six);

  std::cout << "scale=30 p=13 hamming=6 referee=literal-crt-two-flag-upper-relaxation\n";
  std::cout << "orders=(1,2,3,5,6,10,15,30) unit-counts=(1,1,2,4,2,4,8,8)\n";
  std::cout << "hereditary-grammar=at-least-two-orders-divisible-by-each-of-2,3,5\n";
  std::cout << "hereditary-words=" << grammar.words.size()
            << " weighted-states/support=" << grammar.weighted_states
            << " labelled-support-order-contexts="
            << 924ULL * grammar.words.size()
            << " raw-labelled-states=" << 924ULL * grammar.weighted_states
            << '\n';
  std::cout << "bad-prime-intersections-unweighted bits(2,3,5)="
            << format_array(grammar.intersections) << '\n';
  std::cout << "bad-prime-intersections-weighted bits(2,3,5)="
            << format_array(grammar.weighted_intersections) << '\n';
  for (int oi = 0; oi < 8; ++oi)
    std::cout << "D" << ORDERS[oi]
              << " ratio-cardinalities=" << format_cards(table, oi) << '\n';
  std::cout << "scalar-feasible-owners-per-context="
            << format_map(scalar.feasible_owner_histogram) << '\n';
  std::cout << "scalar-supports-by-survivor-count="
            << format_map(scalar.support_survivor_histogram) << '\n';
  std::cout << "scalar-survivors=" << scalar.rows.size()
            << " supports=" << scalar.support_count
            << " literal-unit-words=" << scalar.literal_survivor_words
            << " capacity-vectors=" << scalar.capacity_vectors.size() << '\n';
  std::cout << "scalar-multiplicities="
            << format_multiplicities(scalar.multiplicity_histogram) << '\n';
  std::cout << "normalized-six-owner-keys=" << six.memo.size()
            << " labelled-six-owner-checks=" << six.owner_checks << '\n';
  std::cout << "mod6-anchor-bank-size="
            << format_map(six.anchor_bank_histogram) << '\n';
  std::cout << "mod6-relaxed-bound=" << format_map(six.bound_histogram) << '\n';
  std::cout << "mod6-live-owners-per-context="
            << format_map(six.live_owner_histogram) << '\n';
  std::cout << "mod6-all-owner-residual-rows=" << six.residual_rows.size()
            << '\n';
  std::cout << "normalized-residual-owner-keys ten=" << ten.ten_memo.size()
            << " exact=" << ten.exact_memo.size()
            << " labelled-residual-owner-checks=" << ten.owner_checks << '\n';
  std::cout << "mod10-anchor-bank-size="
            << format_map(ten.ten_anchor_bank_histogram) << '\n';
  std::cout << "mod10-relaxed-bound="
            << format_map(ten.ten_bound_histogram) << '\n';
  std::cout << "mod10-terminal-owner-checks=" << ten.terminal_checks
            << " of " << ten.owner_checks << '\n';
  std::cout << "residual-exact-maximum-union="
            << format_map(ten.exact_maximum_histogram) << '\n';
  std::cout << "mod10-minus-exact="
            << format_map(ten.ten_minus_exact_histogram) << '\n';
  std::cout << "residual-exact-bank-size="
            << format_map(ten.exact_bank_histogram) << '\n';
  std::cout << "residual-exact-reachable-total=" << ten.exact_reachable_total
            << " largest-bank=" << ten.largest_exact_bank << '\n';
  std::cout << "literal-cover-implies-flag-cover=checked-on-all-"
            << six.owner_checks << "-mod6-obligations-and-all-"
            << ten.owner_checks
            << "-residual-mod10-obligations; all-residual-mod10-bounds-below-30\n";
  std::cout << "owner-tournament observable=(U10,U6,exact-max,scalar-capacity,bank-size); "
               "switch=harder-key-wins; tie-path=coordinate-order\n";
  std::cout << "owner-tournament fingerprints=all " << six.residual_rows.size()
            << " transitive scores=(0,1,2,3,4,5) triangles=0 "
               "SCCs=(1,1,1,1,1,1) Hamiltonian-paths=1\n";
  std::cout << "owner-tournament ties=" << format_map(ten.tie_histogram)
            << " natural-order-flips=" << format_map(ten.flip_histogram)
            << '\n';
  std::cout << "proof-carriers=Z/6 fibres anchor D1,D2,D3,D6 then Z/10 fibres "
               "anchor D1,D2,D5,D10\n";
  std::cout << "carrier-preserves=absolute owner upper bounds and literal-cover "
               "implication; destroys=nonanchor overlaps, transverse positions, "
               "and shared-unit compatibility\n";
  std::cout << "alternate-carrier-audit=vertices considered runners,gaps,sections,"
               "boundaries,wall-events,residues,cover-arcs,Fourier-modes,matroid-"
               "circuits,and-owner-obligations; only quotient fibres retain the "
               "monotone absolute threshold used here\n";
  std::cout << "challenged-assumption=tournament vertices need not be runners; "
               "owner obligations rank difficulty while the proof vertices are "
               "six-then-ten quotient fibres\n";
  std::cout << "FNV64 crt-bases=" << hex64(table.base_digest)
            << " literal-masks=" << hex64(table.mask_digest)
            << " scalar-cards=" << hex64(table.card_digest) << '\n';
  std::cout << "FNV64 grammar=" << hex64(grammar.digest)
            << " scalar=" << hex64(scalar.digest)
            << " mod6-audit=" << hex64(six.digest)
            << " residual-mod10-exact=" << hex64(ten.digest) << '\n';
}
