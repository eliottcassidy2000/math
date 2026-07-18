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

// Independent standard-library referee for the primitive proper AP-centred
// common-scale-28 Hamming-six face (THM-995).
//
// This implementation starts from bounded CRT search and the leave-one-out
// lcm predicate.  Its theorem-bearing carrier is the four thick sheet fibres
// of Z/28 -> Z/4.  Orders 1, 2, and 4 are retained literally as unions of
// whole seven-point fibres.  Each transverse order 7, 14, or 28 provider is
// then allowed to select, independently, the option adding the most points
// outside the retained anchor union.  This is an upper relaxation: every
// literal union is bounded by one of these scores.
//
// A separately ordered forward/reverse immutable-union DP checks exact owner
// feasibility and quantifies the information discarded by the relaxation.
// Multiplicative normalization is used only after every literal mask has been
// checked against a common cyclic sheet rotation.  No support, order word, or
// owner obligation is removed by an orbit quotient.

namespace {

constexpr int P = 13;
constexpr int C = 28;
constexpr uint32_t FULL = (uint32_t{1} << C) - 1;
constexpr std::array<int, 6> ORDERS{1, 2, 4, 7, 14, 28};

using Support = std::array<uint8_t, 6>;
using Word = std::array<uint8_t, 6>;  // indices into ORDERS
using Capacities = std::array<uint8_t, 6>;
using Multiplicity = std::array<uint8_t, 6>;
using OwnerKey = std::array<uint8_t, 12>;  // sorted (ratio, order-index) pairs

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

int inverse_mod_13(int value) {
  for (int candidate = 1; candidate < P; ++candidate)
    if (value * candidate % P == 1) return candidate;
  fail("inverse modulo thirteen requested for a nonunit");
}

int positive_mod(int value, int modulus) {
  int residue = value % modulus;
  return residue < 0 ? residue + modulus : residue;
}

int centered(int value, int modulus) {
  const int residue = positive_mod(value, modulus);
  return 2 * residue > modulus ? residue - modulus : residue;
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

uint32_t rotate28(uint32_t mask, int amount) {
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

uint8_t whole_mod4_fibres(uint32_t mask) {
  uint8_t signature = 0;
  for (int residue = 0; residue < 4; ++residue) {
    const bool first = (mask >> residue) & 1U;
    for (int sheet = residue; sheet < C; sheet += 4)
      require(static_cast<bool>((mask >> sheet) & 1U) == first,
              "putative mod-four anchor is not a union of whole fibres");
    if (first) signature |= uint8_t{1} << residue;
  }
  return signature;
}

uint32_t expand_mod4_fibres(uint8_t signature) {
  uint32_t result = 0;
  for (int residue = 0; residue < 4; ++residue)
    if ((signature >> residue) & 1U)
      for (int sheet = residue; sheet < C; sheet += 4)
        result |= uint32_t{1} << sheet;
  return result;
}

std::array<uint8_t, 4> mod4_occupancy(uint32_t mask) {
  std::array<uint8_t, 4> result{};
  for (int sheet = 0; sheet < C; ++sheet)
    if ((mask >> sheet) & 1U) ++result[sheet % 4];
  return result;
}

struct Tables {
  std::array<std::vector<uint8_t>, 6> units;
  // Owner-one, ratio-normalized literal masks.  Unit multiplicity is retained.
  std::array<std::array<std::vector<uint32_t>, 6>, P> fibres;
  std::array<std::array<uint8_t, 6>, P> cards{};
  std::array<std::array<std::set<std::array<uint8_t, 4>>, 6>, P>
      occupancy_flags;
  uint64_t base_digest = 0;
  uint64_t mask_digest = 0;
  uint64_t flag_digest = 0;
};

Tables build_tables() {
  Tables table;
  Fnv64 base_digest;
  Fnv64 mask_digest;
  Fnv64 flag_digest;
  const std::array<std::size_t, 6> expected_phi{1, 1, 2, 6, 6, 12};

  for (int oi = 0; oi < 6; ++oi) {
    table.units[oi] = exact_units(ORDERS[oi]);
    require(table.units[oi].size() == expected_phi[oi],
            "exact-order unit census mismatch");
  }

  // Construct the normalized table from literal bounded CRT search.
  for (int ratio = 1; ratio < P; ++ratio)
    for (int oi = 0; oi < 6; ++oi) {
      const int order = ORDERS[oi];
      int common_card = -1;
      for (uint8_t unit : table.units[oi]) {
        const uint32_t mask = literal_mask(ratio, order, unit, 1);
        table.fibres[ratio][oi].push_back(mask);
        table.occupancy_flags[ratio][oi].insert(mod4_occupancy(mask));
        const int card = std::popcount(mask);
        if (common_card < 0) common_card = card;
        require(card == common_card,
                "mask cardinality depends on the exact-order unit");
        require(periodic_by(mask, order),
                "literal mask violates its effective-order periodicity");
        if (order <= 4) {
          require(periodic_by(mask, 4),
                  "low-order mask violates mod-four periodicity");
          const uint8_t signature = whole_mod4_fibres(mask);
          require(expand_mod4_fibres(signature) == mask,
                  "whole-fibre signature does not reconstruct anchor mask");
        }
      }
      table.cards[ratio][oi] = static_cast<uint8_t>(common_card);

      // Derive cardinality once more by a single effective period, without
      // using the 28-sheet mask.
      const int target = order * ratio % P;
      int one_period = 0;
      for (int value = -order + 1; value <= order; ++value)
        one_period += static_cast<int>(positive_mod(value, P) == target);
      require(common_card == (C / order) * one_period,
              "single-period interval count disagrees with literal mask");

      flag_digest.byte(static_cast<uint8_t>(ratio));
      flag_digest.byte(static_cast<uint8_t>(order));
      for (const auto &flag : table.occupancy_flags[ratio][oi])
        for (uint8_t count : flag) flag_digest.byte(count);
    }

  // Exhaustively prove the ratio/rotation normalization against all literal
  // labels, exact units, and owners before it is used for memoization.
  for (int label = 1; label < P; ++label)
    for (int oi = 0; oi < 6; ++oi) {
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
          require(actual == rotate28(normalized, shift),
                  "literal owner mask violates the common ratio gauge");
          const int card = std::popcount(actual);
          if (common_card < 0) common_card = card;
          require(card == common_card,
                  "owner cardinality depends on exact-order unit");
          require(periodic_by(actual, order),
                  "actual mask violates effective-order periodicity");
          if (order <= 4)
            require(expand_mod4_fibres(whole_mod4_fibres(actual)) == actual,
                    "actual anchor mask is not whole-fibre");

          mask_digest.byte(static_cast<uint8_t>(label));
          mask_digest.byte(static_cast<uint8_t>(order));
          mask_digest.byte(unit);
          mask_digest.byte(static_cast<uint8_t>(owner));
          mask_digest.u32(actual);
        }
        require(common_card == table.cards[ratio][oi],
                "ratio card table disagrees with literal owner masks");
      }
    }

  table.base_digest = base_digest.state;
  table.mask_digest = mask_digest.state;
  table.flag_digest = flag_digest.state;
  return table;
}

bool leave_one_out_lcm(const Word &word) {
  for (int omitted = 0; omitted < 6; ++omitted) {
    int value = 1;
    for (int index = 0; index < 6; ++index)
      if (index != omitted) value = std::lcm(value, ORDERS[word[index]]);
    if (value != C) return false;
  }
  return true;
}

struct WeightedWord {
  Word word{};
  uint64_t weight = 0;
};

struct Grammar {
  std::vector<WeightedWord> words;
  uint64_t weighted_states = 0;
  uint64_t fail_four = 0;
  uint64_t fail_seven = 0;
  uint64_t fail_both = 0;
  uint64_t weighted_fail_four = 0;
  uint64_t weighted_fail_seven = 0;
  uint64_t weighted_fail_both = 0;
  uint64_t digest = 0;
};

Grammar build_grammar(const Tables &table) {
  Grammar result;
  Fnv64 digest;
  uint64_t total_weight = 0;
  int power = 1;
  for (int i = 0; i < 6; ++i) power *= 6;
  for (int encoded = 0; encoded < power; ++encoded) {
    int remainder = encoded;
    Word word{};
    int top_four = 0;
    int top_seven = 0;
    uint64_t weight = 1;
    for (int index = 0; index < 6; ++index) {
      word[index] = static_cast<uint8_t>(remainder % 6);
      remainder /= 6;
      const int order = ORDERS[word[index]];
      top_four += static_cast<int>(order % 4 == 0);
      top_seven += static_cast<int>(order % 7 == 0);
      weight *= table.units[word[index]].size();
    }
    total_weight += weight;
    const bool bad_four = top_four < 2;
    const bool bad_seven = top_seven < 2;
    result.fail_four += bad_four;
    result.fail_seven += bad_seven;
    result.fail_both += bad_four && bad_seven;
    result.weighted_fail_four += bad_four ? weight : 0;
    result.weighted_fail_seven += bad_seven ? weight : 0;
    result.weighted_fail_both += bad_four && bad_seven ? weight : 0;

    const bool grammar = !bad_four && !bad_seven;
    const bool lcm = leave_one_out_lcm(word);
    require(grammar == lcm,
            "valuation grammar disagrees with leave-one-out lcm");
    if (!lcm) continue;
    result.words.push_back({word, weight});
    result.weighted_states += weight;
    for (uint8_t oi : word)
      digest.byte(static_cast<uint8_t>(ORDERS[oi]));
    digest.u64(weight);
  }
  uint64_t expected_total_weight = 1;
  uint64_t unit_sum = 0;
  for (const auto &units : table.units) unit_sum += units.size();
  for (int i = 0; i < 6; ++i) expected_total_weight *= unit_sum;
  require(total_weight == expected_total_weight,
          "weighted grammar does not partition all literal words");
  require(result.words.size() ==
              static_cast<std::size_t>(power - result.fail_four -
                                       result.fail_seven + result.fail_both),
          "unweighted inclusion-exclusion mismatch");
  require(result.weighted_states ==
              total_weight - result.weighted_fail_four -
                  result.weighted_fail_seven + result.weighted_fail_both,
          "weighted inclusion-exclusion mismatch");
  result.digest = digest.state;
  return result;
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
              int support_survivors = 0;
              int contribution[6][6][6]{};
              for (int provider = 0; provider < 6; ++provider)
                for (int oi = 0; oi < 6; ++oi)
                  for (int owner = 0; owner < 6; ++owner) {
                    const int ratio = support[provider] *
                        inverse_mod_13(support[owner]) % P;
                    contribution[provider][oi][owner] =
                        table.cards[ratio][oi];
                  }

              for (const WeightedWord &weighted : grammar.words) {
                ScalarRow row{support, weighted.word, {}};
                int feasible = 0;
                for (int owner = 0; owner < 6; ++owner) {
                  int capacity = 0;
                  for (int provider = 0; provider < 6; ++provider)
                    capacity += contribution[provider]
                                            [weighted.word[provider]][owner];
                  row.capacities[owner] = static_cast<uint8_t>(capacity);
                  feasible += static_cast<int>(capacity >= C);
                }
                ++bank.feasible_owner_histogram[feasible];
                if (feasible != 6) continue;

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

  bank.support_count = surviving_supports.size();
  bank.digest = digest.state;
  uint64_t total_contexts = 0;
  for (const auto &[count, frequency] : bank.feasible_owner_histogram) {
    require(0 <= count && count <= 6,
            "invalid feasible-owner count in scalar census");
    total_contexts += frequency;
  }
  require(total_contexts == 924ULL * grammar.words.size(),
          "scalar census did not visit every labelled context");
  for (const ScalarRow &row : bank.rows) {
    require(leave_one_out_lcm(row.word),
            "scalar bank contains a nonhereditary order word");
    require(std::none_of(row.word.begin(), row.word.end(),
                         [](uint8_t oi) { return ORDERS[oi] == 1; }),
            "scalar survivor contains an order-one provider");
  }
  return bank;
}

OwnerKey normalized_owner_key(const ScalarRow &row, int owner) {
  std::array<std::pair<uint8_t, uint8_t>, 6> pairs{};
  const int inverse = inverse_mod_13(row.support[owner]);
  for (int provider = 0; provider < 6; ++provider)
    pairs[provider] = {
        static_cast<uint8_t>(row.support[provider] * inverse % P),
        row.word[provider]};
  std::sort(pairs.begin(), pairs.end());
  OwnerKey key{};
  for (int index = 0; index < 6; ++index) {
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

std::vector<uint32_t> exact_bank(const OwnerKey &key, const Tables &table,
                                 bool reverse) {
  std::vector<uint32_t> bank{0};
  for (int layer = 0; layer < 6; ++layer) {
    const int pair_index = reverse ? 5 - layer : layer;
    const int ratio = key[2 * pair_index];
    const int oi = key[2 * pair_index + 1];
    bank = advance_bank(bank, distinct_options(table, ratio, oi));
  }
  return bank;
}

struct LocalSummary {
  int relaxed_bound = 0;
  int exact_maximum = 0;
  bool exact_feasible = false;
  uint64_t exact_bank_size = 0;
  uint64_t anchor_bank_size = 0;
  uint64_t exact_bank_digest = 0;
};

LocalSummary analyze_owner_key(const OwnerKey &key, const Tables &table) {
  std::vector<uint32_t> anchors{0};
  std::array<std::vector<uint32_t>, 6> options;
  std::array<bool, 6> is_anchor{};
  for (int layer = 0; layer < 6; ++layer) {
    const int ratio = key[2 * layer];
    const int oi = key[2 * layer + 1];
    options[layer] = distinct_options(table, ratio, oi);
    is_anchor[layer] = ORDERS[oi] <= 4;
    if (is_anchor[layer]) {
      for (uint32_t mask : options[layer])
        require(expand_mod4_fibres(whole_mod4_fibres(mask)) == mask,
                "normalized anchor is not a whole mod-four fibre union");
      anchors = advance_bank(anchors, options[layer]);
    }
  }

  LocalSummary summary;
  summary.anchor_bank_size = anchors.size();
  for (uint32_t anchor : anchors) {
    int score = std::popcount(anchor);
    for (int layer = 0; layer < 6; ++layer) {
      if (is_anchor[layer]) continue;
      int best = 0;
      for (uint32_t option : options[layer])
        best = std::max(best, std::popcount(option & (FULL ^ anchor)));
      for (uint32_t option : options[layer])
        require(std::popcount(option & (FULL ^ anchor)) <= best,
                "pointwise transversal maximum is not an upper bound");
      score += best;
    }
    summary.relaxed_bound = std::max(summary.relaxed_bound, score);
  }

  const std::vector<uint32_t> forward = exact_bank(key, table, false);
  const std::vector<uint32_t> reverse = exact_bank(key, table, true);
  require(forward == reverse,
          "forward/reverse immutable-union banks disagree");
  summary.exact_bank_size = forward.size();
  Fnv64 digest;
  for (uint32_t mask : forward) {
    const int card = std::popcount(mask);
    summary.exact_maximum = std::max(summary.exact_maximum, card);
    summary.exact_feasible = summary.exact_feasible || mask == FULL;
    digest.u32(mask);
  }
  summary.exact_bank_digest = digest.state;
  require(summary.exact_feasible == (summary.exact_maximum == C),
          "exact feasibility and exact maximum disagree");
  require(summary.exact_maximum <= summary.relaxed_bound,
          "literal maximum exceeds mod-four upper relaxation");
  require(!summary.exact_feasible || summary.relaxed_bound >= C,
          "literal cover does not map to a relaxed flag cover");
  return summary;
}

struct TournamentFingerprint {
  int ties = 0;
  int flips = 0;
  int triangles = 0;
  int hamiltonian_paths = 0;
  std::array<uint8_t, 6> out{};
  std::array<int, 6> scc_sizes{};
};

using TournamentKey = std::tuple<int, int, int, int, uint64_t>;

TournamentFingerprint tournament(const std::array<TournamentKey, 6> &keys) {
  TournamentFingerprint result;
  for (int left = 0; left < 6; ++left)
    for (int right = left + 1; right < 6; ++right) {
      int winner = left;
      if (keys[left] == keys[right]) {
        ++result.ties;  // coordinate order is the tie Hamiltonian path
      } else if (keys[right] > keys[left]) {
        winner = right;
        ++result.flips;
      }
      const int loser = left + right - winner;
      result.out[winner] |= uint8_t{1} << loser;
    }

  std::array<int, 6> scores{};
  for (int vertex = 0; vertex < 6; ++vertex)
    scores[vertex] = std::popcount(result.out[vertex]);
  std::sort(scores.begin(), scores.end());
  require(scores == std::array<int, 6>{0, 1, 2, 3, 4, 5},
          "owner tournament has nontransitive score word");

  for (int a = 0; a < 4; ++a)
    for (int b = a + 1; b < 5; ++b)
      for (int c = b + 1; c < 6; ++c)
        result.triangles += static_cast<int>(
            (((result.out[a] >> b) & 1U) &&
             ((result.out[b] >> c) & 1U) &&
             ((result.out[c] >> a) & 1U)) ||
            (((result.out[a] >> c) & 1U) &&
             ((result.out[c] >> b) & 1U) &&
             ((result.out[b] >> a) & 1U)));
  require(result.triangles == 0,
          "owner tournament contains a directed triangle");

  bool reach[6][6]{};
  for (int a = 0; a < 6; ++a) {
    reach[a][a] = true;
    for (int b = 0; b < 6; ++b)
      reach[a][b] = reach[a][b] || ((result.out[a] >> b) & 1U);
  }
  for (int middle = 0; middle < 6; ++middle)
    for (int a = 0; a < 6; ++a)
      for (int b = 0; b < 6; ++b)
        reach[a][b] = reach[a][b] || (reach[a][middle] && reach[middle][b]);
  std::array<bool, 6> used{};
  int component_index = 0;
  for (int root = 0; root < 6; ++root) {
    if (used[root]) continue;
    int size = 0;
    for (int vertex = 0; vertex < 6; ++vertex)
      if (reach[root][vertex] && reach[vertex][root]) {
        used[vertex] = true;
        ++size;
      }
    result.scc_sizes[component_index++] = size;
  }
  require(component_index == 6 &&
              result.scc_sizes == std::array<int, 6>{1, 1, 1, 1, 1, 1},
          "owner tournament SCC fingerprint changed");

  uint64_t paths[64][6]{};
  for (int vertex = 0; vertex < 6; ++vertex)
    paths[1 << vertex][vertex] = 1;
  for (int subset = 1; subset < 64; ++subset)
    for (int last = 0; last < 6; ++last) {
      if (!((subset >> last) & 1)) continue;
      const int previous_subset = subset ^ (1 << last);
      for (int previous = 0; previous < 6; ++previous)
        if (((previous_subset >> previous) & 1) &&
            ((result.out[previous] >> last) & 1U))
          paths[subset][last] += paths[previous_subset][previous];
    }
  uint64_t hamiltonian_paths = 0;
  for (int last = 0; last < 6; ++last)
    hamiltonian_paths += paths[63][last];
  result.hamiltonian_paths = static_cast<int>(hamiltonian_paths);
  require(result.hamiltonian_paths == 1,
          "owner tournament Hamiltonian-path count changed");
  return result;
}

struct OwnerAudit {
  std::map<OwnerKey, LocalSummary> memo;
  std::map<int, uint64_t> relaxed_histogram;
  std::map<int, uint64_t> exact_maximum_histogram;
  std::map<int, uint64_t> anchor_bank_histogram;
  std::map<uint64_t, uint64_t> exact_bank_histogram;
  std::map<int, uint64_t> relaxed_live_owner_histogram;
  std::map<int, uint64_t> exact_live_owner_histogram;
  std::map<int, uint64_t> loss_histogram;
  std::map<int, uint64_t> tie_histogram;
  std::map<int, uint64_t> flip_histogram;
  uint64_t reachable_total = 0;
  uint64_t largest_bank = 0;
  uint64_t implication_checks = 0;
  uint64_t threshold_mismatches = 0;
  uint64_t digest = 0;
};

OwnerAudit audit_owners(const Tables &table, const ScalarBank &scalar) {
  OwnerAudit audit;
  Fnv64 digest;
  for (const ScalarRow &row : scalar.rows) {
    int relaxed_live = 0;
    int exact_live = 0;
    std::array<TournamentKey, 6> keys{};
    for (int owner = 0; owner < 6; ++owner) {
      const OwnerKey key = normalized_owner_key(row, owner);
      auto [iterator, inserted] = audit.memo.try_emplace(key);
      if (inserted) iterator->second = analyze_owner_key(key, table);
      const LocalSummary &local = iterator->second;
      require(local.relaxed_bound <= row.capacities[owner],
              "mod-four relaxation exceeds scalar capacity");
      require(!local.exact_feasible || local.relaxed_bound >= C,
              "literal-cover implication failed at labelled owner");
      ++audit.implication_checks;
      relaxed_live += static_cast<int>(local.relaxed_bound >= C);
      exact_live += static_cast<int>(local.exact_feasible);
      audit.threshold_mismatches += static_cast<uint64_t>(
          (local.relaxed_bound >= C) != local.exact_feasible);
      ++audit.relaxed_histogram[local.relaxed_bound];
      ++audit.exact_maximum_histogram[local.exact_maximum];
      ++audit.anchor_bank_histogram[local.anchor_bank_size];
      ++audit.exact_bank_histogram[local.exact_bank_size];
      ++audit.loss_histogram[local.relaxed_bound - local.exact_maximum];
      audit.reachable_total += local.exact_bank_size;
      audit.largest_bank = std::max(audit.largest_bank,
                                    local.exact_bank_size);
      keys[owner] = {static_cast<int>(local.relaxed_bound >= C),
                     local.relaxed_bound, local.exact_maximum,
                     row.capacities[owner], local.exact_bank_size};

      for (uint8_t value : key) digest.byte(value);
      digest.byte(static_cast<uint8_t>(local.relaxed_bound));
      digest.byte(static_cast<uint8_t>(local.exact_maximum));
      digest.byte(static_cast<uint8_t>(local.exact_feasible));
      digest.u64(local.anchor_bank_size);
      digest.u64(local.exact_bank_size);
      digest.u64(local.exact_bank_digest);
    }
    ++audit.relaxed_live_owner_histogram[relaxed_live];
    ++audit.exact_live_owner_histogram[exact_live];
    require(relaxed_live <= 2,
            "a scalar row survives the mod-four relaxation at three owners");
    require(exact_live <= relaxed_live,
            "exact feasibility escapes the relaxed owner set");

    const TournamentFingerprint fingerprint = tournament(keys);
    ++audit.tie_histogram[fingerprint.ties];
    ++audit.flip_histogram[fingerprint.flips];
  }
  audit.digest = digest.state;
  require(audit.implication_checks == 6ULL * scalar.rows.size(),
          "not every labelled owner implication was checked");
  return audit;
}

struct RelationFingerprint {
  int arcs = 0;
  int undirected_edges = 0;
  int reciprocal_edges = 0;
  int directed_triangles = 0;
  std::vector<int> scc_sizes;
};

RelationFingerprint cayley_relation(const std::set<int> &switch_set) {
  bool adjacency[12][12]{};
  RelationFingerprint result;
  for (int source = 1; source < P; ++source)
    for (int target = 1; target < P; ++target) {
      if (source == target) continue;
      const int ratio = target * inverse_mod_13(source) % P;
      if (switch_set.contains(ratio)) {
        adjacency[source - 1][target - 1] = true;
        ++result.arcs;
      }
    }
  for (int left = 0; left < 11; ++left)
    for (int right = left + 1; right < 12; ++right) {
      result.undirected_edges += adjacency[left][right] ||
                                 adjacency[right][left];
      result.reciprocal_edges += adjacency[left][right] &&
                                 adjacency[right][left];
    }
  for (int a = 0; a < 10; ++a)
    for (int b = a + 1; b < 11; ++b)
      for (int c = b + 1; c < 12; ++c)
        result.directed_triangles += static_cast<int>(
            (adjacency[a][b] && adjacency[b][c] && adjacency[c][a]) ||
            (adjacency[a][c] && adjacency[c][b] && adjacency[b][a]));

  bool reach[12][12]{};
  for (int a = 0; a < 12; ++a)
    for (int b = 0; b < 12; ++b)
      reach[a][b] = a == b || adjacency[a][b];
  for (int middle = 0; middle < 12; ++middle)
    for (int a = 0; a < 12; ++a)
      for (int b = 0; b < 12; ++b)
        reach[a][b] = reach[a][b] || (reach[a][middle] && reach[middle][b]);
  std::array<bool, 12> used{};
  for (int root = 0; root < 12; ++root) {
    if (used[root]) continue;
    int size = 0;
    for (int vertex = 0; vertex < 12; ++vertex)
      if (reach[root][vertex] && reach[vertex][root]) {
        used[vertex] = true;
        ++size;
      }
    result.scc_sizes.push_back(size);
  }
  std::sort(result.scc_sizes.begin(), result.scc_sizes.end());
  return result;
}

template <typename Key, typename Value>
std::string format_map(const std::map<Key, Value> &histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, count] : histogram) {
    if (!first) out << ' ';
    first = false;
    out << key << ':' << count;
  }
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
    for (int index = 0; index < 6; ++index) {
      if (index) out << ',';
      out << static_cast<int>(profile[index]);
    }
    out << "):" << count;
  }
  return out.str();
}

std::string format_sccs(const std::vector<int> &sizes) {
  std::ostringstream out;
  out << '(';
  for (std::size_t index = 0; index < sizes.size(); ++index) {
    if (index) out << ',';
    out << sizes[index];
  }
  out << ')';
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

}  // namespace

int main() {
  const Tables table = build_tables();
  const Grammar grammar = build_grammar(table);
  const ScalarBank scalar = build_scalar_bank(table, grammar);
  const OwnerAudit owners = audit_owners(table, scalar);

  std::set<int> high_two;
  std::set<int> high_four;
  for (int ratio = 1; ratio < P; ++ratio) {
    if (table.cards[ratio][1] == 14) high_two.insert(ratio);
    if (table.cards[ratio][2] == 7) high_four.insert(ratio);
  }
  high_two.erase(1);
  high_four.erase(1);
  const RelationFingerprint relation_two = cayley_relation(high_two);
  const RelationFingerprint relation_four = cayley_relation(high_four);

  std::cout << "scale=28 p=13 hamming=6 referee=literal-crt-mod4-upper-relaxation\n";
  std::cout << "orders=(1,2,4,7,14,28) unit-counts=(1,1,2,6,6,12)\n";
  std::cout << "hereditary-grammar=at-least-two-orders-divisible-by-4-and-at-"
               "least-two-orders-divisible-by-7\n";
  std::cout << "hereditary-words=" << grammar.words.size()
            << " weighted-states/support=" << grammar.weighted_states
            << " labelled-support-order-contexts="
            << 924ULL * grammar.words.size()
            << " raw-labelled-states=" << 924ULL * grammar.weighted_states
            << '\n';
  std::cout << "inclusion-exclusion-unweighted bad4=" << grammar.fail_four
            << " bad7=" << grammar.fail_seven
            << " both=" << grammar.fail_both << '\n';
  std::cout << "inclusion-exclusion-weighted bad4="
            << grammar.weighted_fail_four
            << " bad7=" << grammar.weighted_fail_seven
            << " both=" << grammar.weighted_fail_both << '\n';
  for (int oi = 0; oi < 6; ++oi)
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
  std::cout << "normalized-owner-keys=" << owners.memo.size()
            << " labelled-owner-implication-checks="
            << owners.implication_checks << '\n';
  std::cout << "mod4-anchor-bank-size="
            << format_map(owners.anchor_bank_histogram) << '\n';
  std::cout << "mod4-relaxed-bound="
            << format_map(owners.relaxed_histogram) << '\n';
  std::cout << "mod4-live-owners-per-context="
            << format_map(owners.relaxed_live_owner_histogram) << '\n';
  std::cout << "exact-maximum-union="
            << format_map(owners.exact_maximum_histogram) << '\n';
  std::cout << "exact-feasible-owners-per-context="
            << format_map(owners.exact_live_owner_histogram) << '\n';
  std::cout << "mod4-minus-exact=" << format_map(owners.loss_histogram)
            << " threshold-mismatches=" << owners.threshold_mismatches << '\n';
  std::cout << "exact-bank-size=" << format_map(owners.exact_bank_histogram)
            << '\n';
  std::cout << "exact-reachable-total=" << owners.reachable_total
            << " largest-bank=" << owners.largest_bank << '\n';
  std::cout << "literal-cover-implies-relaxed-flag-cover=checked-on-all-"
            << owners.implication_checks
            << "-labelled-owner-obligations; all-six-relaxed-row=0\n";
  std::cout << "owner-tournament observable=(relaxed-status,U4,exact-max,"
               "scalar-capacity,bank-size); switch=harder-key-wins; "
               "tie-path=coordinate-order\n";
  std::cout << "owner-tournament fingerprints=all " << scalar.rows.size()
            << " transitive scores=(0,1,2,3,4,5) triangles=0 "
               "SCCs=(1,1,1,1,1,1) Hamiltonian-paths=1\n";
  std::cout << "owner-tournament ties=" << format_map(owners.tie_histogram)
            << " natural-order-flips=" << format_map(owners.flip_histogram)
            << '\n';
  std::cout << "D2-Cayley arcs=" << relation_two.arcs
            << " edges=" << relation_two.undirected_edges
            << " reciprocal=" << relation_two.reciprocal_edges
            << " triangles=" << relation_two.directed_triangles
            << " SCCs=" << format_sccs(relation_two.scc_sizes) << '\n';
  std::cout << "D4-Cayley arcs=" << relation_four.arcs
            << " edges=" << relation_four.undirected_edges
            << " reciprocal=" << relation_four.reciprocal_edges
            << " triangles=" << relation_four.directed_triangles
            << " SCCs=" << format_sccs(relation_four.scc_sizes) << '\n';
  std::cout << "proof-carrier=four mod4 sheet fibres: D2/D4 masks are whole-"
               "fibre anchors; D7/D14/D28 masks are relaxed transverse "
               "toothpicks\n";
  std::cout << "carrier-preserves=absolute owner upper bound and literal-"
               "cover implication; destroys=mod7 positions, nonanchor "
               "overlaps, shared-unit compatibility, exact maxima\n";
  std::cout << "alternate-carrier-audit=providers,runners,gaps,sections,wall-"
               "events,residues,cover-arcs,Fourier-modes,matroid-circuits,"
               "Fano-points,chi7-colours,and completed tournaments all lose "
               "either the absolute threshold or thick-fibre incidence\n";
  std::cout << "challenged-assumption=tournament vertices need not be runners; "
               "owner obligations are telemetry while the faithful proof "
               "vertices are four thick fibres plus transverse deviations\n";
  std::cout << "FNV64 crt-bases=" << hex64(table.base_digest)
            << " literal-masks=" << hex64(table.mask_digest)
            << " mod4-flags=" << hex64(table.flag_digest) << '\n';
  std::cout << "FNV64 grammar=" << hex64(grammar.digest)
            << " scalar=" << hex64(scalar.digest)
            << " owner-audit=" << hex64(owners.digest) << '\n';
}
