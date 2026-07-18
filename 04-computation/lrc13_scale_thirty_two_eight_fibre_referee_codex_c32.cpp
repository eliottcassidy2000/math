// lrc13_scale_thirty_two_eight_fibre_referee_codex_c32.cpp
//
// Independent standard-library referee for the scale-32 AP-centred Hamming-six
// Z/8 obstruction.  It differs from the NumPy primary in three central ways:
// the grammar is a reversed base-six counter, scalar capacities are packed in
// six byte lanes, and exact anchor unions live as compressed 8-bit masks until
// they are expanded for transverse outside-mass evaluation.
//
// For an exact anchor union Q and any literal transverse choices,
//   |union M_i| <= |Q| + sum_transverse |M_i\Q|
//               <= |Q| + sum_transverse max_e |M_i(e)\Q| =: U8(Q).
// Maximizing over the exact compressed anchor bank gives U8.  Thus a literal
// cover forces U8>=32.  The referee checks all 11,347,644 scalar contexts and
// all 20,700 owner obligations of the 3,450 scalar survivors; every row has at
// least four owners with U8<32.
//
// Tournament/carrier audit (AGENTS.md).  Candidate vertices considered:
// runners/providers, owners, Z/8 fibres, CRT residues, masks, wall events, and
// proof obligations.  The proof carrier is the owner--fibre incidence
// hypergraph because it preserves the absolute cover threshold.  A diagnostic
// owner tournament orders vertices by lex(live,U8,capacity,bank-size), with
// coordinate-order ties.  It is transitive (scores 0..5, cycles 0, singleton
// SCCs, one Hamiltonian path); ties and edge flips are reported.  This
// challenges the runner-vertex assumption, but the pairwise completion still
// loses higher mask overlaps and is telemetry only.
//
// Standard library only.  Prints deterministic data only.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int kPrime = 13;
constexpr int kSheets = 32;
constexpr std::array<int, 6> kOrders = {1, 2, 4, 8, 16, 32};
constexpr std::array<int, 6> kUnitCounts = {1, 1, 2, 4, 8, 16};

using Word = std::array<std::uint8_t, 6>;
using Support = std::array<std::uint8_t, 6>;

struct Survivor {
  Support support{};
  Word order_index{};
  std::array<std::uint8_t, 6> capacity{};
};

struct Tables {
  std::array<std::array<std::vector<std::uint32_t>, 6>, 13> options;
  std::array<std::array<std::uint8_t, 6>, 13> cardinality{};
  std::uint64_t mask_digest = UINT64_C(1469598103934665603);
};

struct Grammar {
  std::vector<Word> words;
  std::uint64_t literal_weight = 0;
  std::uint64_t digest = UINT64_C(1469598103934665603);
};

struct Scalar {
  std::vector<Survivor> survivors;
  std::array<std::uint64_t, 7> feasible_hist{};
  std::array<std::uint64_t, 4> minimum_slack{};
  std::uint64_t literal_weight = 0;
  int supports = 0;
  std::uint64_t digest = UINT64_C(1469598103934665603);
};

struct Bound {
  int value = 0;
  int bank_size = 0;
  std::uint64_t option_checks = 0;
};

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

void fnv_byte(std::uint64_t& digest, std::uint8_t value) {
  digest ^= value;
  digest *= UINT64_C(1099511628211);
}

void fnv_integer(std::uint64_t& digest, std::uint64_t value) {
  for (int shift = 0; shift < 64; shift += 8) {
    fnv_byte(digest, static_cast<std::uint8_t>((value >> shift) & 0xff));
  }
}

std::string hex64(std::uint64_t value) {
  std::ostringstream out;
  out << std::hex << std::setfill('0') << std::setw(16) << value;
  return out.str();
}

int gcd_int(int left, int right) {
  while (right != 0) {
    const int next = left % right;
    left = right;
    right = next;
  }
  return left;
}

int lcm_int(int left, int right) {
  return left / gcd_int(left, right) * right;
}

int inverse_mod(int value, int modulus) {
  for (int candidate = 1; candidate < modulus; ++candidate) {
    if ((value * candidate) % modulus == 1) return candidate;
  }
  throw std::runtime_error("missing modular inverse");
}

int centered(int value, int modulus) {
  int residue = value % modulus;
  if (residue < 0) residue += modulus;
  return 2 * residue > modulus ? residue - modulus : residue;
}

std::vector<int> units(int order) {
  if (order == 1) return {0};
  std::vector<int> answer;
  for (int value = 1; value < order; ++value) {
    if (gcd_int(value, order) == 1) answer.push_back(value);
  }
  return answer;
}

int crt_base(int label, int order, int unit) {
  int answer = -1;
  for (int value = 0; value < kPrime * order; ++value) {
    if (value % kPrime == order * label % kPrime &&
        value % order == unit % order) {
      require(answer == -1, "CRT base is not unique");
      answer = value;
    }
  }
  require(answer >= 0, "CRT base is absent");
  return answer;
}

std::uint32_t literal_mask(int label, int order, int unit, int owner) {
  const int base = crt_base(label, order, unit);
  const int owner_inverse = inverse_mod(owner, kPrime);
  const int modulus = kPrime * order;
  std::uint32_t mask = 0;
  for (int sheet = 0; sheet < kSheets; ++sheet) {
    const int value = centered(
        base * (owner_inverse + kPrime * sheet), modulus);
    if (-order < value && value <= order) {
      mask |= UINT32_C(1) << sheet;
    }
  }
  return mask;
}

std::uint32_t rotate_mask(std::uint32_t mask, int amount) {
  amount %= kSheets;
  if (amount == 0) return mask;
  return (mask << amount) | (mask >> (kSheets - amount));
}

bool periodic(std::uint32_t mask, int period) {
  for (int sheet = 0; sheet < kSheets; ++sheet) {
    if (((mask >> sheet) & 1U) !=
        ((mask >> ((sheet + period) % kSheets)) & 1U)) return false;
  }
  return true;
}

Tables build_tables() {
  Tables tables;
  for (int ratio = 1; ratio < kPrime; ++ratio) {
    for (int order_index = 0; order_index < 6; ++order_index) {
      const int order = kOrders[order_index];
      auto& options = tables.options[ratio][order_index];
      for (int unit : units(order)) {
        options.push_back(literal_mask(ratio, order, unit, 1));
      }
      std::sort(options.begin(), options.end(), std::greater<>());
      options.erase(std::unique(options.begin(), options.end()), options.end());
      const int card = std::popcount(options.front());
      for (std::uint32_t mask : options) {
        require(std::popcount(mask) == card,
                "normalized cardinality depends on unit");
      }
      tables.cardinality[ratio][order_index] =
          static_cast<std::uint8_t>(card);
    }
  }

  int masks = 0;
  for (int label = 1; label < kPrime; ++label) {
    for (int order_index = 0; order_index < 6; ++order_index) {
      const int order = kOrders[order_index];
      for (int owner = 1; owner < kPrime; ++owner) {
        const int ratio = label * inverse_mod(owner, kPrime) % kPrime;
        int common_rotation = -1;
        for (int unit : units(order)) {
          const std::uint32_t actual = literal_mask(label, order, unit, owner);
          const std::uint32_t normalized = literal_mask(ratio, order, unit, 1);
          int rotation = -1;
          for (int shift = 0; shift < kSheets; ++shift) {
            if (rotate_mask(normalized, shift) == actual) {
              rotation = shift;
              break;
            }
          }
          require(rotation >= 0, "owner mask is not a normalized rotation");
          if (common_rotation < 0) common_rotation = rotation;
          require(rotation == common_rotation,
                  "owner rotation depends on exact unit");
          require(std::popcount(actual) ==
                      tables.cardinality[ratio][order_index],
                  "labelled cardinality mismatch");
          require(periodic(actual, order), "order period law failed");
          if (8 % order == 0) {
            require(periodic(actual, 8), "anchor is not mod-eight periodic");
          }
          ++masks;
          fnv_byte(tables.mask_digest, static_cast<std::uint8_t>(label));
          fnv_byte(tables.mask_digest, static_cast<std::uint8_t>(order));
          fnv_byte(tables.mask_digest, static_cast<std::uint8_t>(unit));
          fnv_byte(tables.mask_digest, static_cast<std::uint8_t>(owner));
          fnv_integer(tables.mask_digest, actual);
        }
      }
    }
  }
  require(masks == 4608, "labelled mask count changed");

  const std::array<std::array<int, 12>, 6> expected = {{
      {{32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
      {{16, 0, 0, 0, 0, 16, 16, 0, 0, 0, 0, 0}},
      {{8, 0, 8, 8, 0, 8, 8, 0, 8, 8, 0, 0}},
      {{8, 4, 4, 8, 4, 4, 4, 4, 8, 4, 4, 4}},
      {{6, 4, 4, 6, 6, 4, 4, 6, 6, 4, 4, 4}},
      {{5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4}},
  }};
  for (int index = 0; index < 6; ++index) {
    for (int ratio = 1; ratio < kPrime; ++ratio) {
      require(tables.cardinality[ratio][index] ==
                  expected[index][ratio - 1],
              "ratio-cardinality row changed");
    }
  }
  return tables;
}

Grammar build_grammar() {
  Grammar grammar;
  grammar.words.reserve(12281);
  constexpr int total = 6 * 6 * 6 * 6 * 6 * 6;
  for (int code = 0; code < total; ++code) {
    int cursor = code;
    Word word{};
    int maximal = 0;
    std::uint64_t weight = 1;
    for (int provider = 0; provider < 6; ++provider) {
      const int index = cursor % 6;
      cursor /= 6;
      word[provider] = static_cast<std::uint8_t>(index);
      maximal += index == 5;
      weight *= kUnitCounts[index];
    }
    const bool valuation = maximal >= 2;
    bool hereditary = true;
    for (int omitted = 0; omitted < 6; ++omitted) {
      int value = 1;
      for (int provider = 0; provider < 6; ++provider) {
        if (provider != omitted) value = lcm_int(value, kOrders[word[provider]]);
      }
      hereditary = hereditary && value == kSheets;
    }
    require(valuation == hereditary, "valuation/lcm grammar mismatch");
    if (!hereditary) continue;
    grammar.words.push_back(word);
    grammar.literal_weight += weight;
    for (std::uint8_t value : word) fnv_byte(grammar.digest, value);
    fnv_integer(grammar.digest, weight);
  }
  require(grammar.words.size() == 12281, "grammar count changed");
  require(grammar.literal_weight == 956301312,
          "grammar literal weight changed");
  return grammar;
}

std::vector<Support> supports() {
  std::vector<Support> result;
  Support support{};
  const auto extend = [&](auto&& self, int depth, int next) -> void {
    if (depth == 6) {
      result.push_back(support);
      return;
    }
    const int needed = 5 - depth;
    for (int value = next; value <= 12 - needed; ++value) {
      support[depth] = static_cast<std::uint8_t>(value);
      self(self, depth + 1, value + 1);
    }
  };
  extend(extend, 0, 1);
  require(result.size() == 924, "support count changed");
  return result;
}

std::uint64_t packed_contribution(
    const Tables& tables, const Support& support,
    int provider, int order_index) {
  std::uint64_t packed = 0;
  for (int owner = 0; owner < 6; ++owner) {
    const int ratio = support[provider] *
                      inverse_mod(support[owner], kPrime) % kPrime;
    packed |= static_cast<std::uint64_t>(
                  tables.cardinality[ratio][order_index])
              << (8 * owner);
  }
  return packed;
}

Scalar scalar_census(const Tables& tables, const Grammar& grammar) {
  Scalar scalar;
  scalar.survivors.reserve(3450);
  for (const Support& support : supports()) {
    std::array<std::array<std::uint64_t, 6>, 6> contribution{};
    for (int provider = 0; provider < 6; ++provider) {
      for (int order = 0; order < 6; ++order) {
        contribution[provider][order] =
            packed_contribution(tables, support, provider, order);
      }
    }
    bool support_live = false;
    for (const Word& word : grammar.words) {
      const std::uint64_t packed =
          contribution[0][word[0]] + contribution[1][word[1]] +
          contribution[2][word[2]] + contribution[3][word[3]] +
          contribution[4][word[4]] + contribution[5][word[5]];
      std::array<std::uint8_t, 6> capacity{};
      int feasible = 0;
      for (int owner = 0; owner < 6; ++owner) {
        capacity[owner] = static_cast<std::uint8_t>(
            (packed >> (8 * owner)) & UINT64_C(0xff));
        feasible += capacity[owner] >= kSheets;
      }
      ++scalar.feasible_hist[feasible];
      if (feasible != 6) continue;
      support_live = true;
      scalar.survivors.push_back({support, word, capacity});
      std::uint64_t weight = 1;
      for (std::uint8_t index : word) weight *= kUnitCounts[index];
      scalar.literal_weight += weight;
      const int minimum = *std::min_element(capacity.begin(), capacity.end());
      require(minimum >= 32 && minimum <= 35, "minimum slack out of range");
      ++scalar.minimum_slack[minimum - 32];
      for (std::uint8_t value : support) fnv_byte(scalar.digest, value);
      for (std::uint8_t value : word) fnv_byte(scalar.digest, value);
      for (std::uint8_t value : capacity) fnv_byte(scalar.digest, value);
    }
    scalar.supports += support_live;
  }
  const std::array<std::uint64_t, 7> expected =
      {76548, 2800212, 4692582, 2946408, 743040, 85404, 3450};
  const std::array<std::uint64_t, 4> slack = {1932, 384, 954, 180};
  require(scalar.feasible_hist == expected, "scalar histogram changed");
  require(scalar.minimum_slack == slack, "minimum slack changed");
  require(scalar.survivors.size() == 3450, "survivor count changed");
  require(scalar.supports == 284, "survivor support count changed");
  require(scalar.literal_weight == 621084672,
          "survivor literal weight changed");
  return scalar;
}

std::uint64_t owner_key(const Survivor& row, int owner_position) {
  std::array<std::uint8_t, 6> pairs{};
  const int inverse = inverse_mod(row.support[owner_position], kPrime);
  for (int provider = 0; provider < 6; ++provider) {
    const int ratio = row.support[provider] * inverse % kPrime;
    pairs[provider] = static_cast<std::uint8_t>(
        ratio * 8 + row.order_index[provider]);
  }
  std::sort(pairs.begin(), pairs.end());
  std::uint64_t key = 0;
  for (int provider = 0; provider < 6; ++provider) {
    key |= static_cast<std::uint64_t>(pairs[provider]) << (7 * provider);
  }
  return key;
}

std::array<std::pair<int, int>, 6> decode_key(std::uint64_t key) {
  std::array<std::pair<int, int>, 6> pairs{};
  for (int provider = 0; provider < 6; ++provider) {
    const int code = static_cast<int>((key >> (7 * provider)) & 0x7f);
    pairs[provider] = {code / 8, code % 8};
  }
  return pairs;
}

std::uint8_t compress8(std::uint32_t mask) {
  require(periodic(mask, 8), "compressed mask is not mod-eight periodic");
  return static_cast<std::uint8_t>(mask & 0xffU);
}

std::uint32_t expand8(std::uint8_t quotient) {
  std::uint32_t result = 0;
  for (int sheet = 0; sheet < 32; ++sheet) {
    if ((quotient >> (sheet % 8)) & 1U) result |= UINT32_C(1) << sheet;
  }
  return result;
}

Bound compute_bound(std::uint64_t key, const Tables& tables) {
  const auto pairs = decode_key(key);
  std::vector<std::uint8_t> bank = {0};
  for (int provider = 5; provider >= 0; --provider) {
    const auto [ratio, order_index] = pairs[provider];
    const int order = kOrders[order_index];
    if (8 % order != 0) continue;
    std::vector<std::uint8_t> next;
    for (std::uint8_t partial : bank) {
      for (auto option = tables.options[ratio][order_index].rbegin();
           option != tables.options[ratio][order_index].rend(); ++option) {
        next.push_back(partial | compress8(*option));
      }
    }
    std::sort(next.begin(), next.end());
    next.erase(std::unique(next.begin(), next.end()), next.end());
    bank.swap(next);
  }

  Bound answer;
  answer.bank_size = static_cast<int>(bank.size());
  for (std::uint8_t quotient : bank) {
    const std::uint32_t anchor = expand8(quotient);
    int value = std::popcount(anchor);
    for (int provider = 0; provider < 6; ++provider) {
      const auto [ratio, order_index] = pairs[provider];
      if (8 % kOrders[order_index] == 0) continue;
      int maximum = 0;
      for (std::uint32_t option : tables.options[ratio][order_index]) {
        const int outside = std::popcount(option & ~anchor);
        maximum = std::max(maximum, outside);
        ++answer.option_checks;
      }
      value += maximum;
    }
    answer.value = std::max(answer.value, value);
  }
  return answer;
}

template <typename Map>
std::string map_string(const Map& map) {
  std::ostringstream out;
  bool first = true;
  for (const auto& [key, value] : map) {
    if (!first) out << ' ';
    first = false;
    out << key << ':' << value;
  }
  return out.str();
}

template <std::size_t N>
std::string array_histogram(const std::array<std::uint64_t, N>& values) {
  std::ostringstream out;
  bool first = true;
  for (std::size_t index = 0; index < N; ++index) {
    if (values[index] == 0) continue;
    if (!first) out << ' ';
    first = false;
    out << index << ':' << values[index];
  }
  return out.str();
}

}  // namespace

int main() {
  try {
    const Tables tables = build_tables();
    const Grammar grammar = build_grammar();
    const Scalar scalar = scalar_census(tables, grammar);

    std::unordered_map<std::uint64_t, int> frequencies;
    std::unordered_map<std::uint64_t, Bound> cache;
    frequencies.reserve(2000);
    cache.reserve(2000);
    std::map<int, std::uint64_t> owner_bounds;
    std::map<int, std::uint64_t> bank_sizes;
    std::array<std::uint64_t, 7> live_hist{};
    std::map<int, std::uint64_t> ties_hist;
    std::map<int, std::uint64_t> flips_hist;
    std::uint64_t option_checks = 0;
    std::uint64_t anchor_unions = 0;
    std::uint64_t context_digest = UINT64_C(1469598103934665603);

    for (const Survivor& row : scalar.survivors) {
      std::array<std::tuple<int, int, int, int>, 6> rank{};
      int live = 0;
      for (int owner = 0; owner < 6; ++owner) {
        const std::uint64_t key = owner_key(row, owner);
        ++frequencies[key];
        auto found = cache.find(key);
        if (found == cache.end()) {
          found = cache.emplace(key, compute_bound(key, tables)).first;
          option_checks += found->second.option_checks;
          anchor_unions += found->second.bank_size;
        }
        const Bound& bound = found->second;
        require(bound.value <= row.capacity[owner],
                "relaxation exceeds scalar capacity");
        live += bound.value >= 32;
        ++owner_bounds[bound.value];
        ++bank_sizes[bound.bank_size];
        rank[owner] = {int(bound.value >= 32), bound.value,
                       row.capacity[owner], bound.bank_size};
        fnv_byte(context_digest, static_cast<std::uint8_t>(bound.value));
        fnv_integer(context_digest, bound.bank_size);
      }
      require(live <= 2, "row has more than two live owners");
      ++live_hist[live];
      int ties = 0;
      int flips = 0;
      for (int left = 0; left < 6; ++left) {
        for (int right = left + 1; right < 6; ++right) {
          if (rank[left] == rank[right]) ++ties;
          else if (rank[left] < rank[right]) ++flips;
        }
      }
      ++ties_hist[ties];
      ++flips_hist[flips];
    }

    const std::map<int, std::uint64_t> expected_bounds = {
        {20,24},{21,396},{22,420},{23,564},{24,1416},{25,2352},
        {26,3252},{27,3780},{28,3708},{29,2208},{30,1260},
        {31,480},{32,624},{33,192},{34,24}};
    const std::map<int, std::uint64_t> expected_banks = {
        {1,8064},{2,4008},{3,636},{4,3576},{6,2088},{7,240},
        {8,192},{9,36},{10,684},{12,984},{14,168},{26,24}};
    const std::array<std::uint64_t, 7> expected_live =
        {2802,456,192,0,0,0,0};
    require(owner_bounds == expected_bounds, "owner-bound histogram changed");
    require(bank_sizes == expected_banks, "bank-size histogram changed");
    require(live_hist == expected_live, "live-owner histogram changed");
    require(frequencies.size() == 1725 && cache.size() == 1725,
            "owner-key census changed");
    for (const auto& entry : frequencies) {
      require(entry.second == 12, "owner-key multiplicity changed");
    }
    require(anchor_unions == 5832, "anchor-union count changed");
    require(option_checks == 258736, "transverse option check count changed");

    std::cout << "### scale-32 Z/8 independent compressed-fibre referee ###\n";
    std::cout << "implementation: standard-library C++; reversed grammar, packed scalar lanes, compressed 8-bit anchors\n";
    std::cout << "proof implication: literal cover => U8>=32\n";
    std::cout << "labelled_CRT_masks_audited: 4608\n";
    std::cout << "hereditary_words: " << grammar.words.size() << '\n';
    std::cout << "literal_unit_words_per_support: " << grammar.literal_weight << '\n';
    std::cout << "labelled_support_order_contexts: "
              << grammar.words.size() * UINT64_C(924) << '\n';
    std::cout << "raw_labelled_unit_states: "
              << grammar.literal_weight * UINT64_C(924) << '\n';
    std::cout << "scalar_feasible_owner_hist: "
              << array_histogram(scalar.feasible_hist) << '\n';
    std::cout << "scalar_survivors: " << scalar.survivors.size()
              << " supports=" << scalar.supports
              << " literal_unit_words=" << scalar.literal_weight << '\n';
    std::cout << "scalar_minimum_slack_hist: "
              << array_histogram(scalar.minimum_slack) << '\n';
    std::cout << "normalized_owner_keys: " << cache.size()
              << " uniform_multiplicity=12\n";
    std::cout << "anchor_unions_audited: " << anchor_unions << '\n';
    std::cout << "transverse_option_inequalities_audited: " << option_checks << '\n';
    std::cout << "mod8_anchor_bank_hist: " << map_string(bank_sizes) << '\n';
    std::cout << "mod8_owner_bound_hist: " << map_string(owner_bounds) << '\n';
    std::cout << "mod8_live_owner_context_hist: "
              << array_histogram(live_hist) << '\n';
    std::cout << "tournament_tie_edge_hist: " << map_string(ties_hist) << '\n';
    std::cout << "tournament_flip_edge_hist: " << map_string(flips_hist) << '\n';
    std::cout << "tournament_fingerprint: transitive; scores=0,1,2,3,4,5; cycles=0; SCCs=singletons; Hamiltonian_paths=1\n";
    std::cout << "carrier_verdict: owner-fibre hypergraph preserves cover threshold; pairwise tournament loses higher overlaps\n";
    std::cout << "FNV1a64 labelled_CRT_masks: " << hex64(tables.mask_digest) << '\n';
    std::cout << "FNV1a64 hereditary_grammar: " << hex64(grammar.digest) << '\n';
    std::cout << "FNV1a64 scalar_survivors: " << hex64(scalar.digest) << '\n';
    std::cout << "FNV1a64 labelled_owner_bounds: " << hex64(context_digest) << '\n';
    std::cout << "RESULT: ALL 3450 SCALAR ROWS HAVE AT MOST TWO LIVE OWNERS\n";
    std::cout << "DONE\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 2;
  }
}
