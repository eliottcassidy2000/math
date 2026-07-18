#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Independent standard-library referee for the primitive proper AP-centred
// common-scale-32 Hamming-six face reserved as THM-1096.
//
// The implementation begins with literal bounded CRT search.  It enumerates
// every hereditary labelled support/order context, applies only the necessary
// scalar-capacity filter, and then retains every mask whose effective order
// divides eight as a literal thick-fibre anchor.  For a fixed anchor union Q,
// each order-16/order-32 provider is independently allowed to choose the unit
// that maximizes |M\Q|.  Hence
//
//   |union_i M_i| <= |Q| + sum_(nonanchors i) max_e |M_i(e)\Q|.
//
// Maximizing the right side over the exact anchor bank is a sound upper
// relaxation.  A separately ordered exact immutable-union DP is telemetry:
// it checks the implication, measures the relaxation loss, and never replaces
// the displayed inequality as the proof certificate.  Multiplicative owner
// normalization is used only after exhaustive literal label/owner masks have
// been matched to one common cyclic sheet rotation.
//
// Reproduce from the repository root with:
//
//   clang++ -std=c++20 -O3 -Wall -Wextra -pedantic \
//     04-computation/lrc13_scale_thirty_two_hamming_six_eight_fibre_referee_codex_c32.cpp \
//     -o /tmp/thm1096-c32
//   /tmp/thm1096-c32 | cmp - \
//     05-knowledge/results/lrc13_scale_thirty_two_hamming_six_eight_fibre_referee_codex_c32.out

namespace {

constexpr int PRIME = 13;
constexpr int SCALE = 32;
constexpr int WIDTH = 6;
constexpr uint32_t FULL = std::numeric_limits<uint32_t>::max();
constexpr std::array<int, WIDTH> ORDERS{1, 2, 4, 8, 16, 32};

using Support = std::array<uint8_t, WIDTH>;
using OrderWord = std::array<uint8_t, WIDTH>;  // indices in ORDERS
using Capacities = std::array<uint8_t, WIDTH>;
using Multiplicity = std::array<uint8_t, WIDTH>;
using OwnerKey = std::array<uint8_t, 2 * WIDTH>;

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

int mod(int value, int modulus) {
  const int residue = value % modulus;
  return residue < 0 ? residue + modulus : residue;
}

int inverse_mod_13(int value) {
  for (int candidate = 1; candidate < PRIME; ++candidate)
    if (value * candidate % PRIME == 1) return candidate;
  fail("inverse modulo thirteen requested for a nonunit");
}

int centred_residue(int value, int modulus) {
  const int residue = mod(value, modulus);
  return 2 * residue > modulus ? residue - modulus : residue;
}

std::vector<uint8_t> units_of_exact_order(int order) {
  if (order == 1) return {0};
  std::vector<uint8_t> units;
  for (int value = 1; value < order; ++value)
    if (std::gcd(value, order) == 1)
      units.push_back(static_cast<uint8_t>(value));
  return units;
}

int bounded_crt_base(int label, int order, int unit) {
  int answer = -1;
  for (int candidate = 0; candidate < PRIME * order; ++candidate) {
    if (candidate % PRIME != order * label % PRIME) continue;
    if (candidate % order != unit % order) continue;
    require(answer < 0, "bounded CRT search returned two bases");
    answer = candidate;
  }
  require(answer >= 0, "bounded CRT search returned no base");
  return answer;
}

uint32_t sheet_mask(int label, int order, int unit, int owner) {
  const int base = bounded_crt_base(label, order, unit);
  const int owner_inverse = inverse_mod_13(owner);
  uint32_t mask = 0;
  for (int sheet = 0; sheet < SCALE; ++sheet) {
    const int residue = centred_residue(
        base * (owner_inverse + PRIME * sheet), PRIME * order);
    if (-order < residue && residue <= order)
      mask |= uint32_t{1} << sheet;
  }
  return mask;
}

uint32_t rotate32(uint32_t mask, int amount) {
  amount = mod(amount, SCALE);
  if (amount == 0) return mask;
  return (mask << amount) | (mask >> (SCALE - amount));
}

bool has_period(uint32_t mask, int period) {
  for (int sheet = 0; sheet < SCALE; ++sheet) {
    const int translated = (sheet + period) % SCALE;
    if (((mask >> sheet) & 1U) != ((mask >> translated) & 1U))
      return false;
  }
  return true;
}

uint8_t mod8_signature(uint32_t mask) {
  uint8_t signature = 0;
  for (int residue = 0; residue < 8; ++residue) {
    const bool occupied = (mask >> residue) & 1U;
    for (int sheet = residue; sheet < SCALE; sheet += 8)
      require(static_cast<bool>((mask >> sheet) & 1U) == occupied,
              "putative Z/8 anchor is not a union of thick fibres");
    if (occupied) signature |= uint8_t{1} << residue;
  }
  return signature;
}

uint32_t expand_mod8(uint8_t signature) {
  uint32_t mask = 0;
  for (int residue = 0; residue < 8; ++residue)
    if ((signature >> residue) & 1U)
      for (int sheet = residue; sheet < SCALE; sheet += 8)
        mask |= uint32_t{1} << sheet;
  return mask;
}

struct LiteralTables {
  std::array<std::vector<uint8_t>, WIDTH> units;
  std::array<std::array<std::vector<uint32_t>, WIDTH>, PRIME> masks;
  std::array<std::array<uint8_t, WIDTH>, PRIME> cardinalities{};
  uint64_t crt_digest = 0;
  uint64_t mask_digest = 0;
  uint64_t anchor_digest = 0;
};

LiteralTables build_literal_tables() {
  LiteralTables result;
  const std::array<std::size_t, WIDTH> expected_units{1, 1, 2, 4, 8, 16};
  for (int oi = 0; oi < WIDTH; ++oi) {
    result.units[oi] = units_of_exact_order(ORDERS[oi]);
    require(result.units[oi].size() == expected_units[oi],
            "Euler-unit count mismatch");
  }

  Fnv64 anchor_digest;
  for (int ratio = 1; ratio < PRIME; ++ratio) {
    for (int oi = 0; oi < WIDTH; ++oi) {
      const int order = ORDERS[oi];
      int common_cardinality = -1;
      for (uint8_t unit : result.units[oi]) {
        const uint32_t mask = sheet_mask(ratio, order, unit, 1);
        result.masks[ratio][oi].push_back(mask);
        const int cardinality = std::popcount(mask);
        if (common_cardinality < 0) common_cardinality = cardinality;
        require(cardinality == common_cardinality,
                "local cardinality depends on the unit");
        require(has_period(mask, order),
                "literal mask violates effective-order periodicity");
        if (order <= 8) {
          require(has_period(mask, 8),
                  "D|8 mask violates the Z/8 anchor period");
          const uint8_t signature = mod8_signature(mask);
          require(expand_mod8(signature) == mask,
                  "Z/8 signature does not reconstruct its anchor mask");
          anchor_digest.byte(static_cast<uint8_t>(ratio));
          anchor_digest.byte(static_cast<uint8_t>(order));
          anchor_digest.byte(unit);
          anchor_digest.byte(signature);
        }
      }
      result.cardinalities[ratio][oi] =
          static_cast<uint8_t>(common_cardinality);

      int hits_in_one_period = 0;
      const int target = order * ratio % PRIME;
      for (int value = -order + 1; value <= order; ++value)
        hits_in_one_period += static_cast<int>(mod(value, PRIME) == target);
      require(common_cardinality == (SCALE / order) * hits_in_one_period,
              "literal mask and interval-cardinality formula disagree");
    }
  }

  // Prove the ratio gauge before it is used for owner-key memoization.
  Fnv64 crt_digest;
  Fnv64 mask_digest;
  for (int label = 1; label < PRIME; ++label) {
    for (int oi = 0; oi < WIDTH; ++oi) {
      const int order = ORDERS[oi];
      for (uint8_t unit : result.units[oi]) {
        crt_digest.byte(static_cast<uint8_t>(label));
        crt_digest.byte(static_cast<uint8_t>(order));
        crt_digest.byte(unit);
        crt_digest.u16(
            static_cast<uint16_t>(bounded_crt_base(label, order, unit)));
      }
      for (int owner = 1; owner < PRIME; ++owner) {
        const int ratio = label * inverse_mod_13(owner) % PRIME;
        int common_shift = -1;
        for (int candidate = 0; candidate < SCALE; ++candidate)
          if (mod(inverse_mod_13(owner) + PRIME * candidate, SCALE) == 1) {
            common_shift = candidate;
            break;
          }
        require(common_shift >= 0,
                "owner normalization has no common sheet rotation");

        for (std::size_t ui = 0; ui < result.units[oi].size(); ++ui) {
          const uint8_t unit = result.units[oi][ui];
          const uint32_t actual = sheet_mask(label, order, unit, owner);
          const uint32_t normalized = result.masks[ratio][oi][ui];
          require(actual == rotate32(normalized, common_shift),
                  "literal mask violates ratio/rotation normalization");
          require(std::popcount(actual) ==
                      result.cardinalities[ratio][oi],
                  "literal owner cardinality violates the ratio table");
          require(has_period(actual, order),
                  "literal owner mask violates its exact period");
          if (order <= 8)
            require(expand_mod8(mod8_signature(actual)) == actual,
                    "literal owner anchor is not a whole-fibre union");

          mask_digest.byte(static_cast<uint8_t>(label));
          mask_digest.byte(static_cast<uint8_t>(order));
          mask_digest.byte(unit);
          mask_digest.byte(static_cast<uint8_t>(owner));
          mask_digest.u32(actual);
        }
      }
    }
  }

  result.crt_digest = crt_digest.state;
  result.mask_digest = mask_digest.state;
  result.anchor_digest = anchor_digest.state;
  return result;
}

bool hereditary_lcm(const OrderWord &word) {
  for (int omitted = 0; omitted < WIDTH; ++omitted) {
    int lcm = 1;
    for (int coordinate = 0; coordinate < WIDTH; ++coordinate)
      if (coordinate != omitted)
        lcm = std::lcm(lcm, ORDERS[word[coordinate]]);
    if (lcm != SCALE) return false;
  }
  return true;
}

struct WeightedWord {
  OrderWord word{};
  uint64_t unit_weight = 0;
};

struct Grammar {
  std::vector<WeightedWord> words;
  uint64_t weighted_states = 0;
  uint64_t no_top_words = 0;
  uint64_t one_top_words = 0;
  uint64_t weighted_no_top = 0;
  uint64_t weighted_one_top = 0;
  uint64_t digest = 0;
};

Grammar enumerate_grammar(const LiteralTables &tables) {
  Grammar grammar;
  Fnv64 digest;
  int word_count = 1;
  for (int coordinate = 0; coordinate < WIDTH; ++coordinate) word_count *= 6;

  uint64_t all_weight = 0;
  for (int encoded = 0; encoded < word_count; ++encoded) {
    int value = encoded;
    OrderWord word{};
    int top_count = 0;
    uint64_t weight = 1;
    for (int coordinate = 0; coordinate < WIDTH; ++coordinate) {
      const int oi = value % 6;
      value /= 6;
      word[coordinate] = static_cast<uint8_t>(oi);
      top_count += static_cast<int>(ORDERS[oi] == SCALE);
      weight *= tables.units[oi].size();
    }
    all_weight += weight;
    if (top_count == 0) {
      ++grammar.no_top_words;
      grammar.weighted_no_top += weight;
    } else if (top_count == 1) {
      ++grammar.one_top_words;
      grammar.weighted_one_top += weight;
    }

    const bool valuation_grammar = top_count >= 2;
    require(valuation_grammar == hereditary_lcm(word),
            "top-order grammar disagrees with leave-one-out lcm");
    if (!valuation_grammar) continue;

    grammar.words.push_back({word, weight});
    grammar.weighted_states += weight;
    for (uint8_t oi : word)
      digest.byte(static_cast<uint8_t>(ORDERS[oi]));
    digest.u64(weight);
  }

  uint64_t unit_sum = 0;
  for (const auto &units : tables.units) unit_sum += units.size();
  uint64_t expected_all_weight = 1;
  for (int coordinate = 0; coordinate < WIDTH; ++coordinate)
    expected_all_weight *= unit_sum;
  require(all_weight == expected_all_weight,
          "weighted order words do not partition the literal state space");
  require(grammar.words.size() ==
              static_cast<std::size_t>(word_count - grammar.no_top_words -
                                       grammar.one_top_words),
          "unweighted hereditary grammar count mismatch");
  require(grammar.weighted_states ==
              all_weight - grammar.weighted_no_top -
                  grammar.weighted_one_top,
          "weighted hereditary grammar count mismatch");
  grammar.digest = digest.state;
  return grammar;
}

std::vector<Support> all_supports() {
  std::vector<Support> supports;
  for (int a = 1; a <= 7; ++a)
    for (int b = a + 1; b <= 8; ++b)
      for (int c = b + 1; c <= 9; ++c)
        for (int d = c + 1; d <= 10; ++d)
          for (int e = d + 1; e <= 11; ++e)
            for (int f = e + 1; f <= 12; ++f)
              supports.push_back(
                  {static_cast<uint8_t>(a), static_cast<uint8_t>(b),
                   static_cast<uint8_t>(c), static_cast<uint8_t>(d),
                   static_cast<uint8_t>(e), static_cast<uint8_t>(f)});
  require(supports.size() == 924, "six-support census mismatch");
  return supports;
}

Multiplicity order_multiplicity(const OrderWord &word) {
  Multiplicity result{};
  for (uint8_t oi : word) ++result[oi];
  return result;
}

struct ScalarRow {
  Support support{};
  OrderWord word{};
  Capacities capacities{};
};

struct ScalarCensus {
  std::vector<ScalarRow> survivors;
  std::map<int, uint64_t> feasible_owner_histogram;
  std::map<int, uint64_t> support_survivor_histogram;
  std::map<Multiplicity, uint64_t> multiplicity_histogram;
  std::set<Capacities> capacity_vectors;
  uint64_t literal_unit_words = 0;
  uint64_t order_one_survivors = 0;
  uint64_t digest = 0;
};

ScalarCensus scalar_census(const LiteralTables &tables,
                           const Grammar &grammar) {
  ScalarCensus census;
  Fnv64 digest;

  for (const Support &support : all_supports()) {
    uint8_t local[WIDTH][WIDTH][WIDTH]{};
    for (int provider = 0; provider < WIDTH; ++provider)
      for (int oi = 0; oi < WIDTH; ++oi)
        for (int owner = 0; owner < WIDTH; ++owner) {
          const int ratio = support[provider] *
                            inverse_mod_13(support[owner]) % PRIME;
          local[provider][oi][owner] = tables.cardinalities[ratio][oi];
        }

    int support_survivors = 0;
    for (const WeightedWord &weighted : grammar.words) {
      ScalarRow row{support, weighted.word, {}};
      int feasible_owners = 0;
      for (int owner = 0; owner < WIDTH; ++owner) {
        int capacity = 0;
        for (int provider = 0; provider < WIDTH; ++provider)
          capacity += local[provider][weighted.word[provider]][owner];
        require(capacity <= std::numeric_limits<uint8_t>::max(),
                "scalar capacity does not fit its storage type");
        row.capacities[owner] = static_cast<uint8_t>(capacity);
        feasible_owners += static_cast<int>(capacity >= SCALE);
      }
      ++census.feasible_owner_histogram[feasible_owners];
      if (feasible_owners != WIDTH) continue;

      ++support_survivors;
      census.survivors.push_back(row);
      census.literal_unit_words += weighted.unit_weight;
      ++census.multiplicity_histogram[order_multiplicity(weighted.word)];
      census.capacity_vectors.insert(row.capacities);
      census.order_one_survivors += static_cast<uint64_t>(
          std::any_of(weighted.word.begin(), weighted.word.end(),
                      [](uint8_t oi) { return ORDERS[oi] == 1; }));

      for (uint8_t label : row.support) digest.byte(label);
      for (uint8_t oi : row.word)
        digest.byte(static_cast<uint8_t>(ORDERS[oi]));
      for (uint8_t capacity : row.capacities) digest.byte(capacity);
    }
    ++census.support_survivor_histogram[support_survivors];
  }

  uint64_t contexts = 0;
  for (const auto &[owners, count] : census.feasible_owner_histogram) {
    require(0 <= owners && owners <= WIDTH,
            "invalid scalar feasible-owner count");
    contexts += count;
  }
  require(contexts == 924ULL * grammar.words.size(),
          "scalar scan missed labelled support/order contexts");
  for (const ScalarRow &row : census.survivors)
    require(hereditary_lcm(row.word),
            "scalar bank contains a nonhereditary word");
  census.digest = digest.state;
  return census;
}

OwnerKey owner_key(const ScalarRow &row, int owner_coordinate) {
  std::array<std::pair<uint8_t, uint8_t>, WIDTH> decorated{};
  const int inverse = inverse_mod_13(row.support[owner_coordinate]);
  for (int provider = 0; provider < WIDTH; ++provider)
    decorated[provider] = {
        static_cast<uint8_t>(row.support[provider] * inverse % PRIME),
        row.word[provider]};
  std::sort(decorated.begin(), decorated.end());

  OwnerKey key{};
  for (int coordinate = 0; coordinate < WIDTH; ++coordinate) {
    key[2 * coordinate] = decorated[coordinate].first;
    key[2 * coordinate + 1] = decorated[coordinate].second;
  }
  return key;
}

std::vector<uint32_t> distinct_masks(const LiteralTables &tables, int ratio,
                                     int oi) {
  std::vector<uint32_t> options = tables.masks[ratio][oi];
  std::sort(options.begin(), options.end());
  options.erase(std::unique(options.begin(), options.end()), options.end());
  return options;
}

std::vector<uint32_t> add_provider(const std::vector<uint32_t> &bank,
                                   const std::vector<uint32_t> &options) {
  std::vector<uint32_t> next;
  next.reserve(bank.size() * options.size());
  for (uint32_t partial : bank)
    for (uint32_t option : options) next.push_back(partial | option);
  std::sort(next.begin(), next.end());
  next.erase(std::unique(next.begin(), next.end()), next.end());
  return next;
}

std::vector<uint32_t> immutable_union_bank(
    const std::array<std::vector<uint32_t>, WIDTH> &options, bool reverse) {
  std::vector<uint32_t> bank{0};
  for (int step = 0; step < WIDTH; ++step) {
    const int provider = reverse ? WIDTH - 1 - step : step;
    bank = add_provider(bank, options[provider]);
  }
  return bank;
}

struct LocalOwnerAudit {
  int upper_bound = 0;
  int exact_maximum = 0;
  bool exact_cover = false;
  uint64_t anchor_bank_size = 0;
  uint64_t exact_bank_size = 0;
  uint64_t exact_bank_digest = 0;
};

LocalOwnerAudit analyze_owner(const OwnerKey &key,
                              const LiteralTables &tables) {
  std::array<std::vector<uint32_t>, WIDTH> options;
  std::array<bool, WIDTH> anchor{};
  std::vector<uint32_t> anchor_bank{0};

  for (int provider = 0; provider < WIDTH; ++provider) {
    const int ratio = key[2 * provider];
    const int oi = key[2 * provider + 1];
    options[provider] = distinct_masks(tables, ratio, oi);
    anchor[provider] = ORDERS[oi] <= 8;
    if (!anchor[provider]) continue;
    for (uint32_t mask : options[provider])
      require(expand_mod8(mod8_signature(mask)) == mask,
              "owner-key anchor is not a whole Z/8-fibre union");
    anchor_bank = add_provider(anchor_bank, options[provider]);
  }

  LocalOwnerAudit audit;
  audit.anchor_bank_size = anchor_bank.size();
  for (uint32_t thick_union : anchor_bank) {
    int candidate_bound = std::popcount(thick_union);
    for (int provider = 0; provider < WIDTH; ++provider) {
      if (anchor[provider]) continue;
      int best_outside = 0;
      for (uint32_t option : options[provider])
        best_outside = std::max(
            best_outside, std::popcount(option & ~thick_union));
      candidate_bound += best_outside;
    }
    audit.upper_bound = std::max(audit.upper_bound, candidate_bound);
  }

  const std::vector<uint32_t> forward = immutable_union_bank(options, false);
  const std::vector<uint32_t> reverse = immutable_union_bank(options, true);
  require(forward == reverse,
          "forward/reverse immutable-union banks disagree");
  audit.exact_bank_size = forward.size();
  Fnv64 bank_digest;
  for (uint32_t mask : forward) {
    audit.exact_maximum = std::max(audit.exact_maximum, std::popcount(mask));
    audit.exact_cover = audit.exact_cover || mask == FULL;
    bank_digest.u32(mask);
  }
  audit.exact_bank_digest = bank_digest.state;

  require(audit.exact_cover == (audit.exact_maximum == SCALE),
          "exact cover status disagrees with exact maximum");
  require(audit.exact_maximum <= audit.upper_bound,
          "exact union exceeds the Z/8 upper relaxation");
  require(!audit.exact_cover || audit.upper_bound >= SCALE,
          "literal cover does not imply relaxed cover");
  return audit;
}

template <typename Key>
struct Tournament {
  std::array<uint8_t, WIDTH> out{};
  int ties = 0;
  int natural_flips = 0;
  int triangles = 0;
  int scc_count = 0;
  uint64_t hamiltonian_paths = 0;
};

template <typename Key>
Tournament<Key> make_tournament(const std::array<Key, WIDTH> &keys) {
  Tournament<Key> tournament;
  for (int left = 0; left < WIDTH; ++left) {
    for (int right = left + 1; right < WIDTH; ++right) {
      int winner = left;
      if (keys[left] == keys[right]) {
        ++tournament.ties;
      } else if (keys[right] > keys[left]) {
        winner = right;
        ++tournament.natural_flips;
      }
      const int loser = left + right - winner;
      tournament.out[winner] |= uint8_t{1} << loser;
    }
  }

  std::array<int, WIDTH> scores{};
  for (int vertex = 0; vertex < WIDTH; ++vertex)
    scores[vertex] = std::popcount(tournament.out[vertex]);
  std::sort(scores.begin(), scores.end());
  require(scores == std::array<int, WIDTH>{0, 1, 2, 3, 4, 5},
          "proof-obligation tournament is not transitive");

  for (int a = 0; a < WIDTH - 2; ++a)
    for (int b = a + 1; b < WIDTH - 1; ++b)
      for (int c = b + 1; c < WIDTH; ++c)
        tournament.triangles += static_cast<int>(
            (((tournament.out[a] >> b) & 1U) &&
             ((tournament.out[b] >> c) & 1U) &&
             ((tournament.out[c] >> a) & 1U)) ||
            (((tournament.out[a] >> c) & 1U) &&
             ((tournament.out[c] >> b) & 1U) &&
             ((tournament.out[b] >> a) & 1U)));
  require(tournament.triangles == 0,
          "proof-obligation tournament has a directed triangle");

  bool reach[WIDTH][WIDTH]{};
  for (int source = 0; source < WIDTH; ++source) {
    reach[source][source] = true;
    for (int target = 0; target < WIDTH; ++target)
      reach[source][target] = reach[source][target] ||
                              ((tournament.out[source] >> target) & 1U);
  }
  for (int middle = 0; middle < WIDTH; ++middle)
    for (int source = 0; source < WIDTH; ++source)
      for (int target = 0; target < WIDTH; ++target)
        reach[source][target] =
            reach[source][target] ||
            (reach[source][middle] && reach[middle][target]);
  std::array<bool, WIDTH> used{};
  for (int root = 0; root < WIDTH; ++root) {
    if (used[root]) continue;
    ++tournament.scc_count;
    for (int vertex = 0; vertex < WIDTH; ++vertex)
      if (reach[root][vertex] && reach[vertex][root]) used[vertex] = true;
  }
  require(tournament.scc_count == WIDTH,
          "proof-obligation tournament has a nonsingleton SCC");

  uint64_t paths[1 << WIDTH][WIDTH]{};
  for (int vertex = 0; vertex < WIDTH; ++vertex)
    paths[1 << vertex][vertex] = 1;
  for (int subset = 1; subset < (1 << WIDTH); ++subset)
    for (int last = 0; last < WIDTH; ++last) {
      if (!((subset >> last) & 1)) continue;
      const int prior = subset ^ (1 << last);
      for (int previous = 0; previous < WIDTH; ++previous)
        if (((prior >> previous) & 1) &&
            ((tournament.out[previous] >> last) & 1U))
          paths[subset][last] += paths[prior][previous];
    }
  for (int last = 0; last < WIDTH; ++last)
    tournament.hamiltonian_paths += paths[(1 << WIDTH) - 1][last];
  require(tournament.hamiltonian_paths == 1,
          "proof-obligation tournament has multiple Hamiltonian paths");
  return tournament;
}

using RelaxKey = std::tuple<int, int, int, uint64_t>;
using ExactKey = std::tuple<int, int, int, uint64_t>;

struct GlobalOwnerAudit {
  std::map<OwnerKey, LocalOwnerAudit> memo;
  std::map<uint64_t, uint64_t> anchor_bank_histogram;
  std::map<int, uint64_t> upper_bound_histogram;
  std::map<int, uint64_t> exact_maximum_histogram;
  std::map<uint64_t, uint64_t> exact_bank_histogram;
  std::map<int, uint64_t> relaxation_loss_histogram;
  std::map<int, uint64_t> relaxed_live_histogram;
  std::map<int, uint64_t> exact_live_histogram;
  std::map<int, uint64_t> relax_tie_histogram;
  std::map<int, uint64_t> exact_tie_histogram;
  std::map<int, uint64_t> relax_natural_flip_histogram;
  std::map<int, uint64_t> exact_natural_flip_histogram;
  std::map<int, uint64_t> gauge_flip_histogram;
  uint64_t labelled_implications = 0;
  uint64_t threshold_mismatches = 0;
  uint64_t exact_reachable_total = 0;
  uint64_t largest_exact_bank = 0;
  uint64_t digest = 0;
};

GlobalOwnerAudit audit_all_owners(const LiteralTables &tables,
                                  const ScalarCensus &scalar) {
  GlobalOwnerAudit audit;
  Fnv64 digest;

  for (const ScalarRow &row : scalar.survivors) {
    int relaxed_live = 0;
    int exact_live = 0;
    std::array<RelaxKey, WIDTH> relax_keys{};
    std::array<ExactKey, WIDTH> exact_keys{};

    for (int owner = 0; owner < WIDTH; ++owner) {
      const OwnerKey key = owner_key(row, owner);
      auto [position, inserted] = audit.memo.try_emplace(key);
      if (inserted) position->second = analyze_owner(key, tables);
      const LocalOwnerAudit &local = position->second;

      require(local.upper_bound <= row.capacities[owner],
              "Z/8 upper relaxation exceeds scalar capacity");
      require(!local.exact_cover || local.upper_bound >= SCALE,
              "labelled literal cover escapes the Z/8 relaxation");
      ++audit.labelled_implications;

      const bool relaxed = local.upper_bound >= SCALE;
      relaxed_live += static_cast<int>(relaxed);
      exact_live += static_cast<int>(local.exact_cover);
      audit.threshold_mismatches +=
          static_cast<uint64_t>(relaxed != local.exact_cover);
      ++audit.anchor_bank_histogram[local.anchor_bank_size];
      ++audit.upper_bound_histogram[local.upper_bound];
      ++audit.exact_maximum_histogram[local.exact_maximum];
      ++audit.exact_bank_histogram[local.exact_bank_size];
      ++audit.relaxation_loss_histogram[local.upper_bound -
                                        local.exact_maximum];
      audit.exact_reachable_total += local.exact_bank_size;
      audit.largest_exact_bank =
          std::max(audit.largest_exact_bank, local.exact_bank_size);

      relax_keys[owner] =
          {static_cast<int>(relaxed), local.upper_bound,
           row.capacities[owner], local.anchor_bank_size};
      exact_keys[owner] =
          {static_cast<int>(local.exact_cover), local.exact_maximum,
           row.capacities[owner], local.exact_bank_size};

      for (uint8_t value : key) digest.byte(value);
      digest.u16(static_cast<uint16_t>(local.upper_bound));
      digest.byte(static_cast<uint8_t>(local.exact_maximum));
      digest.byte(static_cast<uint8_t>(local.exact_cover));
      digest.u64(local.anchor_bank_size);
      digest.u64(local.exact_bank_size);
      digest.u64(local.exact_bank_digest);
    }

    ++audit.relaxed_live_histogram[relaxed_live];
    ++audit.exact_live_histogram[exact_live];
    require(relaxed_live <= 2,
            "a scalar survivor is live at three Z/8-relaxed owners");
    require(exact_live <= relaxed_live,
            "exact owner feasibility escapes relaxed owner feasibility");

    const Tournament<RelaxKey> relaxed_tournament =
        make_tournament(relax_keys);
    const Tournament<ExactKey> exact_tournament =
        make_tournament(exact_keys);
    ++audit.relax_tie_histogram[relaxed_tournament.ties];
    ++audit.exact_tie_histogram[exact_tournament.ties];
    ++audit.relax_natural_flip_histogram[
        relaxed_tournament.natural_flips];
    ++audit.exact_natural_flip_histogram[exact_tournament.natural_flips];
    int gauge_flips = 0;
    for (int left = 0; left < WIDTH; ++left)
      for (int right = left + 1; right < WIDTH; ++right)
        gauge_flips += static_cast<int>(
            ((relaxed_tournament.out[left] >> right) & 1U) !=
            ((exact_tournament.out[left] >> right) & 1U));
    ++audit.gauge_flip_histogram[gauge_flips];
  }

  require(audit.labelled_implications ==
              WIDTH * scalar.survivors.size(),
          "not every labelled owner implication was checked");
  audit.digest = digest.state;
  return audit;
}

template <typename Key, typename Count>
std::string format_map(const std::map<Key, Count> &histogram) {
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
    for (int oi = 0; oi < WIDTH; ++oi) {
      if (oi) out << ',';
      out << static_cast<int>(profile[oi]);
    }
    out << "):" << count;
  }
  return out.str();
}

std::string format_cardinalities(const LiteralTables &tables, int oi) {
  std::ostringstream out;
  out << '(';
  for (int ratio = 1; ratio < PRIME; ++ratio) {
    if (ratio > 1) out << ',';
    out << static_cast<int>(tables.cardinalities[ratio][oi]);
  }
  out << ')';
  return out.str();
}

}  // namespace

int main() {
  const LiteralTables tables = build_literal_tables();
  const Grammar grammar = enumerate_grammar(tables);
  const ScalarCensus scalar = scalar_census(tables, grammar);
  const GlobalOwnerAudit owners = audit_all_owners(tables, scalar);

  std::size_t surviving_supports = 0;
  for (const auto &[survivor_count, frequency] :
       scalar.support_survivor_histogram)
    if (survivor_count > 0) surviving_supports += frequency;

  std::cout << "scale=32 p=13 hamming=6 referee=literal-crt-Z8-upper-relaxation\n";
  std::cout << "orders=(1,2,4,8,16,32) unit-counts=(1,1,2,4,8,16)\n";
  std::cout << "hereditary-grammar=at-least-two-coordinates-have-order-32\n";
  std::cout << "hereditary-words=" << grammar.words.size()
            << " weighted-states/support=" << grammar.weighted_states
            << " labelled-support-order-contexts="
            << 924ULL * grammar.words.size()
            << " raw-labelled-states="
            << 924ULL * grammar.weighted_states << '\n';
  std::cout << "grammar-complement-unweighted no-D32="
            << grammar.no_top_words << " exactly-one-D32="
            << grammar.one_top_words << '\n';
  std::cout << "grammar-complement-weighted no-D32="
            << grammar.weighted_no_top << " exactly-one-D32="
            << grammar.weighted_one_top << '\n';
  for (int oi = 0; oi < WIDTH; ++oi)
    std::cout << 'D' << ORDERS[oi] << " ratio-cardinalities="
              << format_cardinalities(tables, oi) << '\n';
  std::cout << "scalar-feasible-owners-per-context="
            << format_map(scalar.feasible_owner_histogram) << '\n';
  std::cout << "scalar-supports-by-survivor-count="
            << format_map(scalar.support_survivor_histogram) << '\n';
  std::cout << "scalar-survivors=" << scalar.survivors.size()
            << " supports=" << surviving_supports
            << " literal-unit-words=" << scalar.literal_unit_words
            << " capacity-vectors=" << scalar.capacity_vectors.size()
            << " order-one-survivors=" << scalar.order_one_survivors << '\n';
  std::cout << "scalar-multiplicities="
            << format_multiplicities(scalar.multiplicity_histogram) << '\n';
  std::cout << "normalized-owner-keys=" << owners.memo.size()
            << " labelled-owner-implication-checks="
            << owners.labelled_implications << '\n';
  std::cout << "Z8-anchor-bank-size="
            << format_map(owners.anchor_bank_histogram) << '\n';
  std::cout << "Z8-upper-bound="
            << format_map(owners.upper_bound_histogram) << '\n';
  std::cout << "Z8-live-owners-per-context="
            << format_map(owners.relaxed_live_histogram) << '\n';
  std::cout << "exact-maximum-union="
            << format_map(owners.exact_maximum_histogram) << '\n';
  std::cout << "exact-feasible-owners-per-context="
            << format_map(owners.exact_live_histogram) << '\n';
  std::cout << "Z8-minus-exact="
            << format_map(owners.relaxation_loss_histogram)
            << " threshold-mismatches=" << owners.threshold_mismatches
            << '\n';
  std::cout << "exact-bank-size="
            << format_map(owners.exact_bank_histogram) << '\n';
  std::cout << "exact-reachable-total=" << owners.exact_reachable_total
            << " largest-bank=" << owners.largest_exact_bank << '\n';
  std::cout << "literal-cover-implies-Z8-cover=checked-on-all-"
            << owners.labelled_implications
            << "-labelled-owner-obligations; all-six-Z8-live-row=0\n";
  std::cout << "owner-tournament vertices=owner-proof-obligations; "
               "observable=lexicographic-sign(K_i-K_j); "
               "K8=(U8>=32,U8,scalar-capacity,anchor-bank-size); "
               "Kexact=(cover,exact-max,scalar-capacity,exact-bank-size)\n";
  std::cout << "owner-tournament switch=K8-to-Kexact; "
               "tie-Hamiltonian-path=increasing-owner-label-coordinate\n";
  std::cout << "owner-tournament fingerprints=both gauges all "
            << scalar.survivors.size()
            << " transitive scores=(0,1,2,3,4,5) triangles=0 "
               "SCCs=(1,1,1,1,1,1) Hamiltonian-paths=1\n";
  std::cout << "owner-tournament K8-ties="
            << format_map(owners.relax_tie_histogram)
            << " Kexact-ties=" << format_map(owners.exact_tie_histogram)
            << '\n';
  std::cout << "owner-tournament K8-natural-flips="
            << format_map(owners.relax_natural_flip_histogram)
            << " Kexact-natural-flips="
            << format_map(owners.exact_natural_flip_histogram) << '\n';
  std::cout << "owner-tournament gauge-edge-flips="
            << format_map(owners.gauge_flip_histogram) << '\n';
  std::cout << "proof-carrier=eight Z8 sheet fibres; D|8 masks are exact "
               "whole-fibre anchors and D16/D32 masks are independently "
               "maximized nonanchor toothpicks\n";
  std::cout << "carrier-preserves=absolute owner upper bound and literal-"
               "cover implication; destroys=within-fibre offsets, "
               "nonanchor overlaps, and shared-unit compatibility\n";
  std::cout << "assumption-challenge=runner tournaments, gaps, fixed "
               "sections, boundaries, wall-crossing events, residues, "
               "cover arcs, Fourier modes, matroid circuits, and bare "
               "proof rankings all lose the absolute threshold or anchor "
               "incidence; owner obligations are the least honest "
               "tournament sidecar, not the proof carrier\n";
  std::cout << "FNV64 crt-bases=" << hex64(tables.crt_digest)
            << " literal-masks=" << hex64(tables.mask_digest)
            << " Z8-anchors=" << hex64(tables.anchor_digest) << '\n';
  std::cout << "FNV64 grammar=" << hex64(grammar.digest)
            << " scalar=" << hex64(scalar.digest)
            << " owner-audit=" << hex64(owners.digest) << '\n';
}
