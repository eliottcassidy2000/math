// lrc13_scale_thirty_complementary_fibre_referee_codex_c30.cpp
//
// Independent standard-library referee for the AP-centred common-scale-30
// Hamming-six complementary-fibre obstruction.
//
// This implementation is deliberately organized differently from the NumPy
// primary.  It:
//   (1) generates every labelled CRT mask by bounded integer search;
//   (2) enumerates the hereditary grammar with a base-eight counter;
//   (3) packs six scalar owner capacities into independent byte lanes;
//   (4) forms anchor unions in compressed Z/6 or Z/10 masks, expanding them to
//       Z/30 only when transverse outside mass is evaluated; and
//   (5) exhaustively enumerates all literal unit words at every one of the 720
//       owner obligations in the 120-row mod-six residual.
//
// Proof implication checked here.  For an anchor-unit choice let Q be its
// literal union.  For every literal choice of the nonanchor units,
//
//   |union_i M_i| <= |Q| + sum_nonanchor |M_i \ Q|
//                 <= |Q| + sum_nonanchor max_e |M_i(e) \ Q|.
//
// Maximizing the right side over the exact compressed anchor bank gives U_6
// or U_10.  Therefore a literal thirty-sheet cover forces BOTH U_6 >= 30 and
// U_10 >= 30.  Besides constructing those exact upper bounds, this referee
// checks the stronger pointwise inequalities |literal union| <= U_6,U_10 for
// every residual literal unit word.  Every residual owner has U_10 < 30.
//
// Tournament/carrier audit (AGENTS.md).  Candidate vertices considered were
// providers/runners, owners, quotient fibres, CRT residues, masks, and cover
// obligations.  The proof carrier is the owner--fibre incidence hypergraph:
// it preserves the absolute thirty-sheet threshold and the paired quotient
// predicates.  A diagnostic owner tournament uses the pairwise observable
// lex(min(U6,U10), U10, U6, scalar capacity), with larger key winning and the
// support-coordinate order as tie gauge.  This completion is transitive, so
// its score histogram is 0,...,5, its SCCs are singletons, it has no directed
// cycle, and its tie Hamiltonian path is unique.  Edge flips and ties are
// reported.  The challenged assumption is that tournament vertices must be
// runners: owner obligations are more faithful here, but even their pairwise
// tournament destroys the higher-order mask overlaps, so it is telemetry and
// not used in the certificate.
//
// Standard library only.  Prints deterministic data only.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <functional>
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
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int kPrime = 13;
constexpr int kSheets = 30;
constexpr std::uint32_t kFullMask = (UINT32_C(1) << kSheets) - 1;
constexpr std::array<int, 8> kOrders = {1, 2, 3, 5, 6, 10, 15, 30};
constexpr std::array<int, 8> kUnitCounts = {1, 1, 2, 4, 2, 4, 8, 8};

using Word = std::array<std::uint8_t, 6>;
using Support = std::array<std::uint8_t, 6>;

struct Survivor {
  Support support{};
  Word order_index{};
  std::array<std::uint8_t, 6> capacity{};
};

struct Bounds {
  int u6 = 0;
  int u10 = 0;
  int bank6 = 0;
  int bank10 = 0;
};

struct Residual {
  Survivor row;
  std::array<std::uint8_t, 6> u6{};
  std::array<std::uint8_t, 6> u10{};
  std::array<std::uint8_t, 6> literal_max{};
};

struct Tables {
  std::array<std::array<std::vector<std::uint32_t>, 8>, 13> normalized;
  std::array<
      std::array<std::array<std::vector<std::uint32_t>, 13>, 8>, 13>
      labelled;
  std::array<std::array<std::uint8_t, 8>, 13> cardinality{};
  std::uint64_t mask_digest = UINT64_C(1469598103934665603);
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

int gcd_int(int a, int b) {
  while (b != 0) {
    const int next = a % b;
    a = b;
    b = next;
  }
  return a;
}

int lcm_int(int a, int b) { return a / gcd_int(a, b) * b; }

int inverse_mod(int value, int modulus) {
  value %= modulus;
  for (int candidate = 1; candidate < modulus; ++candidate) {
    if ((value * candidate) % modulus == 1) return candidate;
  }
  throw std::runtime_error("modular inverse does not exist");
}

int centered(int value, int modulus) {
  int residue = value % modulus;
  if (residue < 0) residue += modulus;
  if (2 * residue > modulus) residue -= modulus;
  return residue;
}

std::vector<int> exact_units(int order) {
  if (order == 1) return {0};
  std::vector<int> units;
  for (int unit = 1; unit < order; ++unit) {
    if (gcd_int(unit, order) == 1) units.push_back(unit);
  }
  return units;
}

int crt_base_literal(int label, int order, int unit) {
  const int modulus = kPrime * order;
  const int target13 = (order * label) % kPrime;
  int answer = -1;
  for (int value = 0; value < modulus; ++value) {
    if (value % kPrime == target13 && value % order == unit % order) {
      require(answer == -1, "CRT representative is not unique");
      answer = value;
    }
  }
  require(answer >= 0, "CRT representative is missing");
  return answer;
}

std::uint32_t literal_mask(int label, int order, int unit, int owner) {
  const int base = crt_base_literal(label, order, unit);
  const int owner_inverse = inverse_mod(owner, kPrime);
  const int modulus = kPrime * order;
  std::uint32_t mask = 0;
  for (int sheet = 0; sheet < kSheets; ++sheet) {
    const int value = centered(base * (owner_inverse + kPrime * sheet), modulus);
    if (-order < value && value <= order) {
      mask |= UINT32_C(1) << sheet;
    }
  }
  return mask;
}

std::uint32_t rotate_mask(std::uint32_t mask, int amount) {
  amount %= kSheets;
  if (amount < 0) amount += kSheets;
  if (amount == 0) return mask;
  return ((mask << amount) | (mask >> (kSheets - amount))) & kFullMask;
}

bool periodic(std::uint32_t mask, int period) {
  for (int sheet = 0; sheet < kSheets; ++sheet) {
    const bool left = ((mask >> sheet) & 1U) != 0;
    const bool right = ((mask >> ((sheet + period) % kSheets)) & 1U) != 0;
    if (left != right) return false;
  }
  return true;
}

Tables build_tables() {
  Tables tables;
  for (int ratio = 1; ratio < kPrime; ++ratio) {
    for (int order_index = 0; order_index < 8; ++order_index) {
      const int order = kOrders[order_index];
      auto& options = tables.normalized[ratio][order_index];
      for (int unit : exact_units(order)) {
        options.push_back(literal_mask(ratio, order, unit, 1));
      }
      std::sort(options.begin(), options.end(), std::greater<>());
      options.erase(std::unique(options.begin(), options.end()), options.end());
      require(!options.empty(), "normalized option bank is empty");
      const int card = std::popcount(options.front());
      for (std::uint32_t mask : options) {
        require(std::popcount(mask) == card,
                "unit-dependent normalized cardinality");
      }
      tables.cardinality[ratio][order_index] =
          static_cast<std::uint8_t>(card);
    }
  }

  int labelled_masks = 0;
  for (int label = 1; label < kPrime; ++label) {
    for (int order_index = 0; order_index < 8; ++order_index) {
      const int order = kOrders[order_index];
      const std::vector<int> units = exact_units(order);
      for (int owner = 1; owner < kPrime; ++owner) {
        const int ratio = label * inverse_mod(owner, kPrime) % kPrime;
        auto& masks = tables.labelled[label][order_index][owner];
        int common_rotation = -1;
        for (int unit : units) {
          const std::uint32_t actual = literal_mask(label, order, unit, owner);
          const std::uint32_t normalized = literal_mask(ratio, order, unit, 1);
          int rotation = -1;
          for (int shift = 0; shift < kSheets; ++shift) {
            if (rotate_mask(normalized, shift) == actual) {
              rotation = shift;
              break;
            }
          }
          require(rotation >= 0, "owner mask is not a rotation of ratio mask");
          if (common_rotation < 0) common_rotation = rotation;
          require(rotation == common_rotation,
                  "owner gauge rotation depends on exact unit");
          require(std::popcount(actual) ==
                      tables.cardinality[ratio][order_index],
                  "labelled and normalized cardinalities disagree");
          if (6 % order == 0) {
            require(periodic(actual, 6), "mod-six anchor is not periodic");
          }
          if (10 % order == 0) {
            require(periodic(actual, 10), "mod-ten anchor is not periodic");
          }
          masks.push_back(actual);
          ++labelled_masks;
          fnv_byte(tables.mask_digest, static_cast<std::uint8_t>(label));
          fnv_byte(tables.mask_digest, static_cast<std::uint8_t>(order));
          fnv_byte(tables.mask_digest, static_cast<std::uint8_t>(unit));
          fnv_byte(tables.mask_digest, static_cast<std::uint8_t>(owner));
          fnv_integer(tables.mask_digest, actual);
        }
      }
    }
  }
  require(labelled_masks == 4320, "labelled CRT mask count changed");

  const std::array<std::array<int, 12>, 8> expected = {{
      {{30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
      {{15, 0, 0, 0, 0, 15, 15, 0, 0, 0, 0, 0}},
      {{10, 0, 0, 10, 10, 0, 0, 10, 10, 0, 0, 0}},
      {{6, 6, 6, 0, 6, 6, 6, 6, 0, 6, 6, 0}},
      {{5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0}},
      {{6, 6, 6, 3, 3, 6, 6, 3, 3, 6, 6, 3}},
      {{6, 4, 4, 4, 4, 6, 6, 4, 4, 4, 4, 4}},
      {{5, 4, 5, 5, 4, 5, 5, 4, 5, 5, 4, 4}},
  }};
  for (int order_index = 0; order_index < 8; ++order_index) {
    for (int ratio = 1; ratio < kPrime; ++ratio) {
      require(tables.cardinality[ratio][order_index] ==
                  expected[order_index][ratio - 1],
              "ratio-cardinality ledger changed");
    }
  }
  return tables;
}

struct Grammar {
  std::vector<Word> words;
  std::uint64_t literal_words_per_support = 0;
  std::uint64_t digest = UINT64_C(1469598103934665603);
};

Grammar build_grammar() {
  Grammar grammar;
  grammar.words.reserve(185193);
  constexpr int total_words = 8 * 8 * 8 * 8 * 8 * 8;
  for (int code = 0; code < total_words; ++code) {
    int cursor = code;
    Word word{};
    std::array<int, 3> flags{};
    std::uint64_t weight = 1;
    for (int provider = 0; provider < 6; ++provider) {
      const int order_index = cursor % 8;
      cursor /= 8;
      word[provider] = static_cast<std::uint8_t>(order_index);
      const int order = kOrders[order_index];
      flags[0] += order % 2 == 0;
      flags[1] += order % 3 == 0;
      flags[2] += order % 5 == 0;
      weight *= kUnitCounts[order_index];
    }
    const bool via_flags = flags[0] >= 2 && flags[1] >= 2 && flags[2] >= 2;
    bool via_lcm = true;
    for (int omitted = 0; omitted < 6; ++omitted) {
      int value = 1;
      for (int provider = 0; provider < 6; ++provider) {
        if (provider != omitted) {
          value = lcm_int(value, kOrders[word[provider]]);
        }
      }
      via_lcm = via_lcm && value == kSheets;
    }
    require(via_flags == via_lcm, "prime flags and hereditary lcm disagree");
    if (!via_flags) continue;
    grammar.words.push_back(word);
    grammar.literal_words_per_support += weight;
    for (std::uint8_t value : word) fnv_byte(grammar.digest, value);
    fnv_integer(grammar.digest, weight);
  }
  require(grammar.words.size() == 185193, "hereditary grammar count changed");
  require(grammar.literal_words_per_support == 636667200,
          "literal grammar weight changed");
  return grammar;
}

std::vector<Support> make_supports() {
  std::vector<Support> supports;
  Support support{};
  const auto extend = [&](auto&& self, int depth, int next) -> void {
    if (depth == 6) {
      supports.push_back(support);
      return;
    }
    const int needed_after = 5 - depth;
    for (int value = next; value <= 12 - needed_after; ++value) {
      support[depth] = static_cast<std::uint8_t>(value);
      self(self, depth + 1, value + 1);
    }
  };
  extend(extend, 0, 1);
  require(supports.size() == 924, "support count changed");
  return supports;
}

std::uint64_t pack_capacities(
    const Tables& tables, const Support& support, int provider, int order_index) {
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

struct ScalarCensus {
  std::vector<Survivor> survivors;
  std::array<std::uint64_t, 7> feasible_hist{};
  std::array<std::uint64_t, 6> minimum_slack_hist{};
  std::uint64_t literal_survivor_words = 0;
  int survivor_supports = 0;
  std::uint64_t digest = UINT64_C(1469598103934665603);
};

ScalarCensus scalar_census(const Tables& tables, const Grammar& grammar) {
  ScalarCensus result;
  result.survivors.reserve(54050);
  const std::vector<Support> supports = make_supports();
  for (const Support& support : supports) {
    std::array<std::array<std::uint64_t, 8>, 6> contribution{};
    for (int provider = 0; provider < 6; ++provider) {
      for (int order_index = 0; order_index < 8; ++order_index) {
        contribution[provider][order_index] =
            pack_capacities(tables, support, provider, order_index);
      }
    }

    bool support_survives = false;
    for (const Word& word : grammar.words) {
      const std::uint64_t packed =
          contribution[0][word[0]] + contribution[1][word[1]] +
          contribution[2][word[2]] + contribution[3][word[3]] +
          contribution[4][word[4]] + contribution[5][word[5]];
      int feasible = 0;
      std::array<std::uint8_t, 6> capacity{};
      for (int owner = 0; owner < 6; ++owner) {
        capacity[owner] = static_cast<std::uint8_t>(
            (packed >> (8 * owner)) & UINT64_C(0xff));
        feasible += capacity[owner] >= kSheets;
      }
      ++result.feasible_hist[feasible];
      if (feasible != 6) continue;

      support_survives = true;
      result.survivors.push_back({support, word, capacity});
      std::uint64_t literal_weight = 1;
      for (std::uint8_t index : word) literal_weight *= kUnitCounts[index];
      result.literal_survivor_words += literal_weight;
      int minimum = 255;
      for (std::uint8_t value : capacity) minimum = std::min(minimum, int(value));
      require(minimum >= kSheets && minimum - kSheets < 6,
              "minimum scalar slack is out of range");
      ++result.minimum_slack_hist[minimum - kSheets];
      for (std::uint8_t label : support) fnv_byte(result.digest, label);
      for (std::uint8_t index : word) fnv_byte(result.digest, index);
      for (std::uint8_t value : capacity) fnv_byte(result.digest, value);
    }
    result.survivor_supports += support_survives;
  }

  const std::array<std::uint64_t, 7> expected_hist = {
      1401966, 36143640, 66874158, 49478260,
      15326622, 1839636, 54050};
  const std::array<std::uint64_t, 6> expected_slack = {
      39152, 9672, 4248, 744, 198, 36};
  require(result.feasible_hist == expected_hist,
          "scalar feasible-owner histogram changed");
  require(result.minimum_slack_hist == expected_slack,
          "scalar minimum-slack histogram changed");
  require(result.survivors.size() == 54050, "scalar survivor count changed");
  require(result.survivor_supports == 772, "scalar survivor support count changed");
  require(result.literal_survivor_words == 64678912,
          "scalar survivor literal weight changed");
  return result;
}

std::uint64_t owner_key(const Survivor& row, int owner_position) {
  std::array<std::uint8_t, 6> pair_code{};
  const int owner_inverse = inverse_mod(row.support[owner_position], kPrime);
  for (int provider = 0; provider < 6; ++provider) {
    const int ratio = row.support[provider] * owner_inverse % kPrime;
    pair_code[provider] = static_cast<std::uint8_t>(
        ratio * 8 + row.order_index[provider]);
  }
  std::sort(pair_code.begin(), pair_code.end());
  std::uint64_t key = 0;
  for (int provider = 0; provider < 6; ++provider) {
    key |= static_cast<std::uint64_t>(pair_code[provider]) << (7 * provider);
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

std::uint16_t compress_periodic(std::uint32_t mask, int modulus) {
  require(periodic(mask, modulus), "attempt to compress a nonperiodic mask");
  std::uint16_t quotient = 0;
  for (int residue = 0; residue < modulus; ++residue) {
    if ((mask >> residue) & 1U) quotient |= UINT16_C(1) << residue;
  }
  return quotient;
}

std::uint32_t expand_periodic(std::uint16_t quotient, int modulus) {
  std::uint32_t mask = 0;
  for (int sheet = 0; sheet < kSheets; ++sheet) {
    if ((quotient >> (sheet % modulus)) & 1U) {
      mask |= UINT32_C(1) << sheet;
    }
  }
  return mask;
}

std::pair<int, int> quotient_bound(
    std::uint64_t key, int modulus, const Tables& tables) {
  const auto pairs = decode_key(key);
  std::vector<std::uint16_t> bank = {0};
  // Reverse traversal is intentional: it differs from the primary's anchor
  // provider order while producing the same exact union bank.
  for (int provider = 5; provider >= 0; --provider) {
    const auto [ratio, order_index] = pairs[provider];
    const int order = kOrders[order_index];
    if (modulus % order != 0) continue;
    std::vector<std::uint16_t> next;
    for (std::uint16_t partial : bank) {
      for (auto option_it = tables.normalized[ratio][order_index].rbegin();
           option_it != tables.normalized[ratio][order_index].rend(); ++option_it) {
        next.push_back(partial | compress_periodic(*option_it, modulus));
      }
    }
    std::sort(next.begin(), next.end());
    next.erase(std::unique(next.begin(), next.end()), next.end());
    bank.swap(next);
  }

  int best = 0;
  for (std::uint16_t quotient_union : bank) {
    const std::uint32_t anchor_union =
        expand_periodic(quotient_union, modulus);
    int bound = std::popcount(anchor_union);
    for (int provider = 0; provider < 6; ++provider) {
      const auto [ratio, order_index] = pairs[provider];
      const int order = kOrders[order_index];
      if (modulus % order == 0) continue;
      int outside = 0;
      for (std::uint32_t option : tables.normalized[ratio][order_index]) {
        outside = std::max(
            outside, std::popcount(option & (~anchor_union & kFullMask)));
      }
      bound += outside;
    }
    best = std::max(best, bound);
  }
  return {best, static_cast<int>(bank.size())};
}

Bounds compute_bounds(std::uint64_t key, const Tables& tables) {
  const auto six = quotient_bound(key, 6, tables);
  const auto ten = quotient_bound(key, 10, tables);
  return {six.first, ten.first, six.second, ten.second};
}

std::string array_string(const Support& values) {
  std::ostringstream out;
  out << '[';
  for (int i = 0; i < 6; ++i) {
    if (i != 0) out << ',';
    out << int(values[i]);
  }
  out << ']';
  return out.str();
}

std::string order_string(const Word& indices) {
  std::ostringstream out;
  out << '[';
  for (int i = 0; i < 6; ++i) {
    if (i != 0) out << ',';
    out << kOrders[indices[i]];
  }
  out << ']';
  return out.str();
}

std::string byte_array_string(const std::array<std::uint8_t, 6>& values) {
  return array_string(values);
}

struct FibreCensus {
  std::array<std::uint64_t, 7> live6{};
  std::array<std::uint64_t, 7> live10{};
  std::array<std::uint64_t, 7> common{};
  std::map<int, std::uint64_t> residual_u10_hist;
  std::map<int, std::uint64_t> residual_literal_max_hist;
  std::map<std::array<int, 8>, std::uint64_t> residual_order_profiles;
  std::map<int, std::uint64_t> tournament_ties;
  std::map<int, std::uint64_t> tournament_flips;
  std::vector<Residual> residuals;
  std::unordered_map<std::uint64_t, int> key_frequency;
  std::unordered_map<std::uint64_t, Bounds> bound_cache;
};

FibreCensus fibre_census(const ScalarCensus& scalar, const Tables& tables) {
  FibreCensus result;
  result.key_frequency.reserve(30000);
  result.bound_cache.reserve(30000);

  for (const Survivor& row : scalar.survivors) {
    std::array<std::uint8_t, 6> u6{};
    std::array<std::uint8_t, 6> u10{};
    int count6 = 0;
    int count10 = 0;
    int count_common = 0;
    std::array<std::tuple<int, int, int, int>, 6> rank{};
    for (int owner = 0; owner < 6; ++owner) {
      const std::uint64_t key = owner_key(row, owner);
      ++result.key_frequency[key];
      auto found = result.bound_cache.find(key);
      if (found == result.bound_cache.end()) {
        found = result.bound_cache.emplace(key, compute_bounds(key, tables)).first;
      }
      const Bounds& bounds = found->second;
      u6[owner] = static_cast<std::uint8_t>(bounds.u6);
      u10[owner] = static_cast<std::uint8_t>(bounds.u10);
      count6 += bounds.u6 >= kSheets;
      count10 += bounds.u10 >= kSheets;
      count_common += std::min(bounds.u6, bounds.u10) >= kSheets;
      rank[owner] = {
          std::min(bounds.u6, bounds.u10), bounds.u10,
          bounds.u6, row.capacity[owner]};
    }
    ++result.live6[count6];
    ++result.live10[count10];
    ++result.common[count_common];

    int ties = 0;
    int flips = 0;
    for (int left = 0; left < 6; ++left) {
      for (int right = left + 1; right < 6; ++right) {
        if (rank[left] == rank[right]) {
          ++ties;
        } else if (rank[left] < rank[right]) {
          ++flips;
        }
      }
    }
    ++result.tournament_ties[ties];
    ++result.tournament_flips[flips];

    if (count6 == 6) {
      Residual residual;
      residual.row = row;
      residual.u6 = u6;
      residual.u10 = u10;
      for (std::uint8_t score : u10) {
        ++result.residual_u10_hist[score];
        require(score < kSheets,
                "a mod-six residual has a mod-ten-live owner");
      }
      std::array<int, 8> profile{};
      for (std::uint8_t index : row.order_index) ++profile[index];
      ++result.residual_order_profiles[profile];
      result.residuals.push_back(residual);
    }
  }

  const std::array<std::uint64_t, 7> expected6 =
      {45110, 7536, 1284, 0, 0, 0, 120};
  const std::array<std::uint64_t, 7> expected10 =
      {33944, 10344, 6288, 2232, 744, 132, 366};
  const std::array<std::uint64_t, 7> expected_common =
      {45998, 6852, 1200, 0, 0, 0, 0};
  const std::map<int, std::uint64_t> expected_residual_scores =
      {{23, 48}, {24, 120}, {25, 192}, {26, 240}, {27, 72}, {28, 48}};
  const std::map<std::array<int, 8>, std::uint64_t> expected_profiles = {
      {{{0, 0, 0, 1, 0, 3, 0, 2}}, 48},
      {{{0, 0, 0, 0, 0, 4, 0, 2}}, 48},
      {{{0, 0, 0, 0, 0, 2, 3, 1}}, 24},
  };
  require(result.live6 == expected6, "mod-six live histogram changed");
  require(result.live10 == expected10, "mod-ten live histogram changed");
  require(result.common == expected_common, "common-live histogram changed");
  require(result.residuals.size() == 120, "mod-six residual count changed");
  require(result.residual_u10_hist == expected_residual_scores,
          "residual mod-ten score histogram changed");
  require(result.residual_order_profiles == expected_profiles,
          "residual order-profile census changed");
  require(result.key_frequency.size() == 27025,
          "normalized owner key count changed");
  require(result.bound_cache.size() == result.key_frequency.size(),
          "bound cache did not cover every owner key");
  for (const auto& entry : result.key_frequency) {
    require(entry.second == 12, "owner key multiplicity is not twelve");
  }
  return result;
}

std::uint64_t literal_residual_audit(
    FibreCensus& fibre, const Tables& tables, std::uint64_t& residual_digest) {
  std::sort(fibre.residuals.begin(), fibre.residuals.end(),
            [](const Residual& left, const Residual& right) {
              if (left.row.support != right.row.support) {
                return left.row.support < right.row.support;
              }
              return left.row.order_index < right.row.order_index;
            });

  std::uint64_t literal_words = 0;
  for (Residual& residual : fibre.residuals) {
    for (int owner = 0; owner < 6; ++owner) {
      std::array<const std::vector<std::uint32_t>*, 6> options{};
      for (int provider = 0; provider < 6; ++provider) {
        options[provider] = &tables.labelled
            [residual.row.support[provider]]
            [residual.row.order_index[provider]]
            [residual.row.support[owner]];
      }

      int literal_max = 0;
      const auto extend = [&](auto&& self, int provider,
                              std::uint32_t union_mask) -> void {
        if (provider == 6) {
          ++literal_words;
          const int size = std::popcount(union_mask);
          literal_max = std::max(literal_max, size);
          // This pointwise check independently validates the two upper
          // relaxations on every residual literal unit word.
          require(size <= residual.u6[owner],
                  "literal union exceeds mod-six relaxation");
          require(size <= residual.u10[owner],
                  "literal union exceeds mod-ten relaxation");
          if (union_mask == kFullMask) {
            require(residual.u6[owner] >= kSheets &&
                        residual.u10[owner] >= kSheets,
                    "literal cover violates paired implication");
          }
          return;
        }
        for (std::uint32_t option : *options[provider]) {
          self(self, provider + 1, union_mask | option);
        }
      };
      extend(extend, 0, 0);
      residual.literal_max[owner] = static_cast<std::uint8_t>(literal_max);
      ++fibre.residual_literal_max_hist[literal_max];
      require(literal_max < kSheets,
              "a residual literal owner obligation covers all sheets");
    }

    for (std::uint8_t value : residual.row.support) fnv_byte(residual_digest, value);
    for (std::uint8_t value : residual.row.order_index) fnv_byte(residual_digest, value);
    for (std::uint8_t value : residual.u6) fnv_byte(residual_digest, value);
    for (std::uint8_t value : residual.u10) fnv_byte(residual_digest, value);
    for (std::uint8_t value : residual.literal_max) fnv_byte(residual_digest, value);
  }
  require(literal_words == 18874368,
          "residual literal owner-word count changed");
  const std::map<int, std::uint64_t> expected_literal_max =
      {{23, 72}, {24, 144}, {25, 192}, {26, 264}, {27, 48}};
  require(fibre.residual_literal_max_hist == expected_literal_max,
          "residual literal-maximum histogram changed");
  return literal_words;
}

template <std::size_t N>
std::string histogram_string(const std::array<std::uint64_t, N>& histogram) {
  std::ostringstream out;
  bool first = true;
  for (std::size_t index = 0; index < N; ++index) {
    if (histogram[index] == 0) continue;
    if (!first) out << ' ';
    first = false;
    out << index << ':' << histogram[index];
  }
  return out.str();
}

std::string map_string(const std::map<int, std::uint64_t>& histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto& [key, value] : histogram) {
    if (!first) out << ' ';
    first = false;
    out << key << ':' << value;
  }
  return out.str();
}

std::string profile_map_string(
    const std::map<std::array<int, 8>, std::uint64_t>& profiles) {
  std::ostringstream out;
  bool first_profile = true;
  for (const auto& [profile, count] : profiles) {
    if (!first_profile) out << ' ';
    first_profile = false;
    out << '[';
    for (int index = 0; index < 8; ++index) {
      if (index != 0) out << ',';
      out << profile[index];
    }
    out << "]:" << count;
  }
  return out.str();
}

}  // namespace

int main() {
  try {
    const Tables tables = build_tables();
    const Grammar grammar = build_grammar();
    const ScalarCensus scalar = scalar_census(tables, grammar);
    FibreCensus fibre = fibre_census(scalar, tables);
    std::uint64_t residual_digest = UINT64_C(1469598103934665603);
    const std::uint64_t literal_words =
        literal_residual_audit(fibre, tables, residual_digest);

    std::cout << "### scale-30 complementary-fibre independent referee ###\n";
    std::cout << "scope: primitive proper AP-centred common-scale Hamming-six owner gate\n";
    std::cout << "implementation: standard-library C++; compressed Z/6 and Z/10 anchor banks\n";
    std::cout << "proof implication: literal cover => U6>=30 and U10>=30\n";
    std::cout << "CRT masks: labelled=4320 independently generated; owner rotations and periods exact\n";
    std::cout << "orders: [1,2,3,5,6,10,15,30] unit_counts: [1,1,2,4,2,4,8,8]\n";
    std::cout << "hereditary_words: " << grammar.words.size() << '\n';
    std::cout << "literal_unit_words_per_support: "
              << grammar.literal_words_per_support << '\n';
    std::cout << "labelled_support_order_contexts: "
              << grammar.words.size() * UINT64_C(924) << '\n';
    std::cout << "raw_labelled_unit_states: "
              << grammar.literal_words_per_support * UINT64_C(924) << '\n';
    std::cout << "scalar_feasible_owner_hist: "
              << histogram_string(scalar.feasible_hist) << '\n';
    std::cout << "scalar_survivors: " << scalar.survivors.size()
              << " supports=" << scalar.survivor_supports
              << " literal_unit_words=" << scalar.literal_survivor_words << '\n';
    std::cout << "scalar_minimum_slack_hist: "
              << histogram_string(scalar.minimum_slack_hist) << '\n';
    std::cout << "normalized_owner_keys: " << fibre.key_frequency.size()
              << " uniform_multiplicity=12\n";
    std::cout << "mod6_live_owner_context_hist: "
              << histogram_string(fibre.live6) << '\n';
    std::cout << "mod10_live_owner_context_hist: "
              << histogram_string(fibre.live10) << '\n';
    std::cout << "common_live_owner_context_hist: "
              << histogram_string(fibre.common) << '\n';
    std::cout << "mod6_all_six_residual_rows: " << fibre.residuals.size() << '\n';
    std::cout << "residual_owner_obligations: " << fibre.residuals.size() * 6 << '\n';
    std::cout << "residual_mod10_score_hist: "
              << map_string(fibre.residual_u10_hist) << '\n';
    std::cout << "residual_order_profiles_[D1,D2,D3,D5,D6,D10,D15,D30]: "
              << profile_map_string(fibre.residual_order_profiles) << '\n';
    std::cout << "residual_literal_owner_unit_words_checked: "
              << literal_words << '\n';
    std::cout << "residual_literal_maximum_hist: "
              << map_string(fibre.residual_literal_max_hist) << '\n';
    std::cout << "tournament_tie_edge_hist: "
              << map_string(fibre.tournament_ties) << '\n';
    std::cout << "tournament_flip_edge_hist: "
              << map_string(fibre.tournament_flips) << '\n';
    std::cout << "tournament_fingerprint: transitive; scores=0,1,2,3,4,5; "
                 "cycles=0; SCCs=singletons; Hamiltonian_paths=1\n";
    std::cout << "carrier_verdict: paired owner-fibre hypergraph preserves cover; "
                 "pairwise tournament loses higher overlaps\n";
    std::cout << "FNV1a64 labelled_CRT_masks: "
              << hex64(tables.mask_digest) << '\n';
    std::cout << "FNV1a64 hereditary_grammar: "
              << hex64(grammar.digest) << '\n';
    std::cout << "FNV1a64 scalar_survivors: "
              << hex64(scalar.digest) << '\n';
    std::cout << "FNV1a64 residual_certificate: "
              << hex64(residual_digest) << "\n\n";

    for (std::size_t index = 0; index < fibre.residuals.size(); ++index) {
      const Residual& residual = fibre.residuals[index];
      std::cout << "residual=" << index
                << " support=" << array_string(residual.row.support)
                << " orders=" << order_string(residual.row.order_index)
                << " U6=" << byte_array_string(residual.u6)
                << " U10=" << byte_array_string(residual.u10)
                << " literal_max=" << byte_array_string(residual.literal_max)
                << '\n';
    }

    std::cout << "\npointwise_relaxation_checks: " << literal_words
              << "/" << literal_words << " passed for both U6 and U10\n";
    std::cout << "residual_mod10_max: 28 < 30\n";
    std::cout << "literal_covers_in_720_residual_obligations: 0\n";
    std::cout << "RESULT: ALL 120 MOD6 RESIDUAL ROWS REJECTED BY MOD10 AT ALL 720 OWNERS\n";
    std::cout << "DONE\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 2;
  }
}
