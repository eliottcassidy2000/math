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
#include <utility>
#include <vector>

// Independent standard-library C++20 referee for the primitive proper
// AP-centred common-scale-34 Hamming-six sheet bank (THM-1125).
//
// The construction starts at the mathematical definition.  It finds the CRT
// base by bounded search, builds all 34 sheet bits literally, enumerates all
// 924 labelled supports and all leave-one-out-lcm order words, applies the
// scalar union bound, and only then evaluates the proof-facing Z/2 anchor
// relaxation.  No orbit quotient or precomputed mask/capacity table is used.
//
// For an owner and a fixed exact anchor union Q, every literal choice obeys
//
//   |union_i M_i| <= |Q| + sum_nonanchor max_u |M_i(u) \ Q|.
//
// Orders 1 and 2 are retained in Q.  The independent maxima for orders 17
// and 34 forget their shared global unit word and all mutual overlaps; this
// only enlarges the bound, so a value below 34 is a sound obstruction.

namespace {

constexpr int P = 13;
constexpr int SCALE = 34;
constexpr int WIDTH = 6;
constexpr std::array<int, 4> ORDERS{1, 2, 17, 34};
constexpr uint64_t FULL = (uint64_t{1} << SCALE) - 1;

using Word = std::array<uint8_t, WIDTH>;
using Support = std::array<uint8_t, WIDTH>;
using Capacities = std::array<uint8_t, WIDTH>;
using Multiplicity = std::array<uint8_t, ORDERS.size()>;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

int positive_mod(int value, int modulus) {
  const int residue = value % modulus;
  return residue < 0 ? residue + modulus : residue;
}

int inverse_mod_13(int value) {
  for (int candidate = 1; candidate < P; ++candidate)
    if (value * candidate % P == 1) return candidate;
  fail("inverse modulo thirteen requested for a nonunit");
}

int centred_residue(int value, int modulus) {
  const int residue = positive_mod(value, modulus);
  return 2 * residue > modulus ? residue - modulus : residue;
}

std::vector<uint8_t> exact_units(int order) {
  if (order == 1) return {0};
  std::vector<uint8_t> answer;
  for (int value = 1; value < order; ++value)
    if (std::gcd(value, order) == 1)
      answer.push_back(static_cast<uint8_t>(value));
  return answer;
}

int literal_crt_base(int label, int order, int unit) {
  int answer = -1;
  int solutions = 0;
  for (int candidate = 0; candidate < P * order; ++candidate) {
    if (candidate % P != order * label % P) continue;
    if (candidate % order != unit % order) continue;
    answer = candidate;
    ++solutions;
  }
  require(solutions == 1, "bounded CRT search was not unique");
  return answer;
}

uint64_t literal_mask(int label, int order, int unit, int owner) {
  const int base = literal_crt_base(label, order, unit);
  const int owner_inverse = inverse_mod_13(owner);
  uint64_t result = 0;
  for (int sheet = 0; sheet < SCALE; ++sheet) {
    const int residue = centred_residue(
        base * (owner_inverse + P * sheet), P * order);
    if (-order < residue && residue <= order)
      result |= uint64_t{1} << sheet;
  }
  return result;
}

bool periodic_by(uint64_t mask, int step) {
  for (int sheet = 0; sheet < SCALE; ++sheet) {
    const int translated = (sheet + step) % SCALE;
    if (((mask >> sheet) & 1U) != ((mask >> translated) & 1U))
      return false;
  }
  return true;
}

struct Fnv64 {
  uint64_t value = 14695981039346656037ULL;

  void byte(uint8_t item) {
    value ^= item;
    value *= 1099511628211ULL;
  }
  void u16(uint16_t item) {
    byte(static_cast<uint8_t>(item));
    byte(static_cast<uint8_t>(item >> 8));
  }
  void u32(uint32_t item) {
    for (int shift = 0; shift < 32; shift += 8)
      byte(static_cast<uint8_t>(item >> shift));
  }
  void u64(uint64_t item) {
    u32(static_cast<uint32_t>(item));
    u32(static_cast<uint32_t>(item >> 32));
  }
};

std::string hex64(uint64_t value) {
  std::ostringstream out;
  out << std::hex << std::setfill('0') << std::setw(16) << value;
  return out.str();
}

struct LiteralTables {
  std::array<std::vector<uint8_t>, ORDERS.size()> units;
  std::array<
      std::array<std::array<std::vector<uint64_t>, P>, ORDERS.size()>, P>
      masks;
  std::array<std::array<std::array<uint8_t, P>, ORDERS.size()>, P>
      cardinalities{};
  uint64_t base_digest = 0;
  uint64_t mask_digest = 0;
};

LiteralTables construct_literal_tables() {
  LiteralTables tables;
  const std::array<std::size_t, 4> expected_phi{1, 1, 16, 16};
  for (std::size_t oi = 0; oi < ORDERS.size(); ++oi) {
    tables.units[oi] = exact_units(ORDERS[oi]);
    require(tables.units[oi].size() == expected_phi[oi],
            "Euler-unit count mismatch");
  }

  Fnv64 bases;
  Fnv64 masks;
  for (int label = 1; label < P; ++label) {
    for (std::size_t oi = 0; oi < ORDERS.size(); ++oi) {
      const int order = ORDERS[oi];
      for (uint8_t unit : tables.units[oi]) {
        bases.byte(static_cast<uint8_t>(label));
        bases.byte(static_cast<uint8_t>(order));
        bases.byte(unit);
        bases.u16(static_cast<uint16_t>(literal_crt_base(label, order, unit)));
      }

      for (int owner = 1; owner < P; ++owner) {
        int common_cardinality = -1;
        auto &choices = tables.masks[label][oi][owner];
        for (uint8_t unit : tables.units[oi]) {
          const uint64_t mask = literal_mask(label, order, unit, owner);
          const int card = std::popcount(mask);
          if (common_cardinality < 0) common_cardinality = card;
          require(card == common_cardinality,
                  "mask cardinality depends on the exact unit");
          require(periodic_by(mask, order),
                  "literal mask violates effective-order periodicity");
          if (order == 1 || order == 2)
            require(periodic_by(mask, 2),
                    "putative Z/2 anchor is not a thick-fibre union");
          choices.push_back(mask);

          masks.byte(static_cast<uint8_t>(label));
          masks.byte(static_cast<uint8_t>(order));
          masks.byte(unit);
          masks.byte(static_cast<uint8_t>(owner));
          masks.u64(mask);
        }
        tables.cardinalities[label][oi][owner] =
            static_cast<uint8_t>(common_cardinality);

        // A second derivation counts admissible centred residues in one
        // effective-order period, independently of the 34-bit construction.
        const int ratio = label * inverse_mod_13(owner) % P;
        const int target = order * ratio % P;
        int interval_hits = 0;
        for (int value = -order + 1; value <= order; ++value)
          interval_hits += positive_mod(value, P) == target;
        require(common_cardinality == (SCALE / order) * interval_hits,
                "literal mask and interval-cardinality formula disagree");
      }
    }
  }

  // The scalar layer depends only on the provider/owner ratio.  Check that
  // covariance on every literal labelled entry; it is not used to skip work.
  for (int label = 1; label < P; ++label)
    for (std::size_t oi = 0; oi < ORDERS.size(); ++oi)
      for (int owner = 1; owner < P; ++owner) {
        const int ratio = label * inverse_mod_13(owner) % P;
        require(tables.cardinalities[label][oi][owner] ==
                    tables.cardinalities[ratio][oi][1],
                "scalar ratio covariance failed");
      }

  tables.base_digest = bases.value;
  tables.mask_digest = masks.value;
  return tables;
}

bool hereditary_word(const Word &word) {
  bool literal = true;
  for (int omitted = 0; omitted < WIDTH; ++omitted) {
    int order = 1;
    for (int coordinate = 0; coordinate < WIDTH; ++coordinate)
      if (coordinate != omitted)
        order = std::lcm(order, ORDERS[word[coordinate]]);
    literal &= order == SCALE;
  }

  int carriers_two = 0;
  int carriers_seventeen = 0;
  for (uint8_t oi : word) {
    carriers_two += ORDERS[oi] % 2 == 0;
    carriers_seventeen += ORDERS[oi] % 17 == 0;
  }
  const bool carrier_grammar =
      carriers_two >= 2 && carriers_seventeen >= 2;
  require(literal == carrier_grammar,
          "leave-one-out lcm and carrier grammar disagree");
  return literal;
}

struct WeightedWord {
  Word order_indices{};
  uint64_t unit_weight = 0;
};

std::vector<WeightedWord> enumerate_hereditary_words(
    const LiteralTables &tables, uint64_t &weighted_total) {
  std::vector<WeightedWord> result;
  weighted_total = 0;
  constexpr int ALL_WORDS = 4 * 4 * 4 * 4 * 4 * 4;
  for (int code = 0; code < ALL_WORDS; ++code) {
    int remaining = code;
    WeightedWord candidate;
    candidate.unit_weight = 1;
    for (int coordinate = 0; coordinate < WIDTH; ++coordinate) {
      const int oi = remaining % static_cast<int>(ORDERS.size());
      remaining /= static_cast<int>(ORDERS.size());
      candidate.order_indices[coordinate] = static_cast<uint8_t>(oi);
      candidate.unit_weight *= tables.units[oi].size();
    }
    if (!hereditary_word(candidate.order_indices)) continue;
    result.push_back(candidate);
    weighted_total += candidate.unit_weight;
  }
  return result;
}

std::vector<Support> enumerate_supports() {
  std::vector<Support> result;
  for (uint32_t bits = 0; bits < (uint32_t{1} << 12); ++bits) {
    if (std::popcount(bits) != WIDTH) continue;
    Support support{};
    int position = 0;
    for (int label = 1; label < P; ++label)
      if ((bits >> (label - 1)) & 1U)
        support[position++] = static_cast<uint8_t>(label);
    require(position == WIDTH, "support popcount/fill mismatch");
    result.push_back(support);
  }
  return result;
}

struct Survivor {
  Support support{};
  Word word{};
  Capacities capacities{};
  uint64_t unit_weight = 0;
};

struct ScalarCensus {
  std::array<uint64_t, WIDTH + 1> owner_capacity_histogram{};
  std::map<int, int> support_survivor_histogram;
  std::set<Capacities> capacity_vectors;
  std::map<Multiplicity, int> multiplicity_profile_histogram;
  std::vector<Survivor> survivors;
  uint64_t contexts = 0;
  uint64_t survivor_unit_states = 0;
  uint64_t survivor_digest = 0;
  int survivor_supports = 0;
};

ScalarCensus scalar_census(const LiteralTables &tables,
                           const std::vector<Support> &supports,
                           const std::vector<WeightedWord> &words) {
  ScalarCensus census;
  Fnv64 digest;
  for (const Support &support : supports) {
    int on_support = 0;
    for (const WeightedWord &weighted_word : words) {
      ++census.contexts;
      Capacities capacities{};
      int live_owners = 0;
      for (int owner_index = 0; owner_index < WIDTH; ++owner_index) {
        int capacity = 0;
        for (int provider = 0; provider < WIDTH; ++provider)
          capacity += tables.cardinalities[support[provider]]
                                          [weighted_word.order_indices[provider]]
                                          [support[owner_index]];
        require(capacity <= 6 * SCALE,
                "scalar capacity does not fit the census type");
        capacities[owner_index] = static_cast<uint8_t>(capacity);
        live_owners += capacity >= SCALE;
      }
      ++census.owner_capacity_histogram[live_owners];
      if (live_owners != WIDTH) continue;

      ++on_support;
      census.survivor_unit_states += weighted_word.unit_weight;
      Multiplicity multiplicity{};
      for (uint8_t oi : weighted_word.order_indices) ++multiplicity[oi];
      ++census.multiplicity_profile_histogram[multiplicity];
      census.capacity_vectors.insert(capacities);
      census.survivors.push_back(
          {support, weighted_word.order_indices, capacities,
           weighted_word.unit_weight});

      for (uint8_t label : support) digest.byte(label);
      for (uint8_t oi : weighted_word.order_indices)
        digest.byte(static_cast<uint8_t>(ORDERS[oi]));
      for (uint8_t capacity : capacities) digest.byte(capacity);
      digest.u64(weighted_word.unit_weight);
    }
    if (on_support != 0) ++census.survivor_supports;
    ++census.support_survivor_histogram[on_support];
  }
  census.survivor_digest = digest.value;
  return census;
}

int z2_relaxed_bound(const LiteralTables &tables,
                     const Survivor &survivor, int owner_index,
                     int &anchor_bank_size) {
  const int owner = survivor.support[owner_index];
  std::set<uint64_t> anchor_bank{0};
  for (int provider = 0; provider < WIDTH; ++provider) {
    const int oi = survivor.word[provider];
    const int order = ORDERS[oi];
    if (order != 1 && order != 2) continue;
    std::set<uint64_t> next;
    for (uint64_t partial : anchor_bank)
      for (uint64_t option :
           tables.masks[survivor.support[provider]][oi][owner])
        next.insert(partial | option);
    anchor_bank = std::move(next);
  }
  require(!anchor_bank.empty(), "empty Z/2 anchor bank");
  anchor_bank_size = static_cast<int>(anchor_bank.size());

  int best = 0;
  for (uint64_t anchor_union : anchor_bank) {
    const uint64_t outside_anchor = FULL ^ anchor_union;
    int bound = std::popcount(anchor_union);
    for (int provider = 0; provider < WIDTH; ++provider) {
      const int oi = survivor.word[provider];
      const int order = ORDERS[oi];
      if (order == 1 || order == 2) continue;
      int best_increment = 0;
      for (uint64_t option :
           tables.masks[survivor.support[provider]][oi][owner])
        best_increment = std::max(
            best_increment, std::popcount(option & outside_anchor));
      bound += best_increment;
    }
    best = std::max(best, bound);
  }
  return best;
}

struct RelaxationCensus {
  std::map<int, int> bound_histogram;
  std::map<int, int> anchor_bank_histogram;
  std::map<std::pair<int, int>, int> own_order_bound_histogram;
  std::map<int, int> live_owner_histogram;
  uint64_t obligations = 0;
  uint64_t digest = 0;
  int global_ceiling = 0;
  int seventeen_carrier_ceiling = 0;
};

RelaxationCensus relaxation_census(const LiteralTables &tables,
                                   const ScalarCensus &scalar) {
  RelaxationCensus census;
  Fnv64 digest;
  for (const Survivor &survivor : scalar.survivors) {
    int live_owners = 0;
    for (int owner_index = 0; owner_index < WIDTH; ++owner_index) {
      int bank_size = 0;
      const int bound =
          z2_relaxed_bound(tables, survivor, owner_index, bank_size);
      const int own_order = ORDERS[survivor.word[owner_index]];
      ++census.obligations;
      ++census.bound_histogram[bound];
      ++census.anchor_bank_histogram[bank_size];
      ++census.own_order_bound_histogram[{own_order, bound}];
      census.global_ceiling = std::max(census.global_ceiling, bound);
      if (own_order % 17 == 0)
        census.seventeen_carrier_ceiling =
            std::max(census.seventeen_carrier_ceiling, bound);
      live_owners += bound >= SCALE;

      for (uint8_t label : survivor.support) digest.byte(label);
      for (uint8_t oi : survivor.word)
        digest.byte(static_cast<uint8_t>(ORDERS[oi]));
      digest.byte(static_cast<uint8_t>(owner_index));
      digest.byte(static_cast<uint8_t>(bank_size));
      digest.byte(static_cast<uint8_t>(bound));
    }
    ++census.live_owner_histogram[live_owners];
  }
  census.digest = digest.value;
  return census;
}

template <class Map>
void print_integer_histogram(const Map &histogram) {
  bool first = true;
  for (const auto &[key, count] : histogram) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << key << ':' << count;
  }
  std::cout << '\n';
}

void print_scalar_table(const LiteralTables &tables) {
  std::cout << "scalar-cardinalities columns=1,2,17,34\n";
  for (int ratio = 1; ratio < P; ++ratio) {
    std::cout << "r=" << ratio << ':';
    for (std::size_t oi = 0; oi < ORDERS.size(); ++oi)
      std::cout << (oi == 0 ? "" : ",")
                << static_cast<int>(tables.cardinalities[ratio][oi][1]);
    std::cout << '\n';
  }
}

void print_own_order_bounds(
    const std::map<std::pair<int, int>, int> &histogram) {
  for (int order : ORDERS) {
    std::cout << "own-order-" << order << "-bounds=";
    bool first = true;
    for (const auto &[key, count] : histogram) {
      if (key.first != order) continue;
      if (!first) std::cout << ',';
      first = false;
      std::cout << key.second << ':' << count;
    }
    if (first) std::cout << "none";
    std::cout << '\n';
  }
}

void print_multiplicity_profiles(
    const std::map<Multiplicity, int> &histogram) {
  bool first_profile = true;
  for (const auto &[profile, count] : histogram) {
    if (!first_profile) std::cout << ';';
    first_profile = false;
    for (std::size_t oi = 0; oi < ORDERS.size(); ++oi) {
      if (oi != 0) std::cout << ',';
      std::cout << static_cast<int>(profile[oi]);
    }
    std::cout << ':' << count;
  }
  std::cout << '\n';
}

}  // namespace

int main() {
  const LiteralTables tables = construct_literal_tables();
  uint64_t weighted_words = 0;
  const std::vector<WeightedWord> words =
      enumerate_hereditary_words(tables, weighted_words);
  const std::vector<Support> supports = enumerate_supports();
  const ScalarCensus scalar = scalar_census(tables, supports, words);
  const RelaxationCensus relaxation = relaxation_census(tables, scalar);

  // Closed-form squarefree carrier checks, independent of enumeration.
  const uint64_t carrier_subsets = (uint64_t{1} << WIDTH) - 1 - WIDTH;
  const uint64_t expected_words = carrier_subsets * carrier_subsets;
  const uint64_t weighted_two =
      uint64_t{2} * 2 * 2 * 2 * 2 * 2 - 1 - WIDTH;
  uint64_t seventeen_to_six = 1;
  for (int i = 0; i < WIDTH; ++i) seventeen_to_six *= 17;
  const uint64_t weighted_seventeen =
      seventeen_to_six - 1 - WIDTH * uint64_t{16};
  const uint64_t expected_weighted_words =
      weighted_two * weighted_seventeen;
  const uint64_t expected_contexts = supports.size() * expected_words;
  const uint64_t raw_labelled_unit_states =
      supports.size() * expected_weighted_words;

  require(supports.size() == 924, "labelled support count mismatch");
  require(words.size() == expected_words,
          "hereditary order-word count mismatch");
  require(weighted_words == expected_weighted_words,
          "weighted hereditary word count mismatch");
  require(scalar.contexts == expected_contexts,
          "labelled scalar context count mismatch");
  require(scalar.owner_capacity_histogram ==
              std::array<uint64_t, 7>{33'112, 387'180, 1'283'004,
                                      1'099'712, 180'156, 18'360, 552},
          "frozen scalar live-owner histogram mismatch");
  require(scalar.survivors.size() == 552,
          "frozen scalar survivor count mismatch");
  require(scalar.survivor_supports == 36,
          "frozen scalar survivor-support count mismatch");
  require(scalar.survivor_unit_states == 36'175'872,
          "frozen scalar survivor unit-weight mismatch");
  require(scalar.support_survivor_histogram ==
              std::map<int, int>{{0, 888}, {15, 24}, {16, 12}},
          "frozen support survivor histogram mismatch");
  require(scalar.multiplicity_profile_histogram ==
              std::map<Multiplicity, int>{
                  {Multiplicity{0, 2, 0, 4}, 36},
                  {Multiplicity{0, 2, 1, 3}, 144},
                  {Multiplicity{0, 2, 2, 2}, 216},
                  {Multiplicity{0, 2, 3, 1}, 144},
                  {Multiplicity{0, 2, 4, 0}, 12}},
          "frozen order-multiplicity histogram mismatch");
  require(scalar.capacity_vectors.size() == 356,
          "frozen capacity-vector count mismatch");
  require(relaxation.obligations == scalar.survivors.size() * WIDTH,
          "owner-obligation count mismatch");
  require(relaxation.anchor_bank_histogram ==
              std::map<int, int>{{1, 3'312}},
          "frozen anchor-bank histogram mismatch");
  require(relaxation.bound_histogram ==
              std::map<int, int>{{25, 48}, {26, 1'056}, {27, 408},
                                 {28, 1'728}, {29, 72}},
          "frozen Z/2 bound histogram mismatch");
  require(relaxation.own_order_bound_histogram ==
              std::map<std::pair<int, int>, int>{
                  {{2, 25}, 48}, {{2, 26}, 288}, {{2, 27}, 408},
                  {{2, 28}, 288}, {{2, 29}, 72}, {{17, 26}, 384},
                  {{17, 28}, 672}, {{34, 26}, 384}, {{34, 28}, 768}},
          "frozen own-order Z/2 histogram mismatch");
  require(relaxation.global_ceiling == 29,
          "frozen global Z/2 ceiling mismatch");
  require(relaxation.seventeen_carrier_ceiling == 28,
          "frozen own-order-17/34 Z/2 ceiling mismatch");
  require(relaxation.live_owner_histogram ==
              std::map<int, int>{{0, 552}},
          "a scalar row survived the Z/2 owner relaxation");
  require(tables.base_digest == 0x8837d6f407684800ULL &&
              tables.mask_digest == 0xbf343e2c845dcc21ULL &&
              scalar.survivor_digest == 0xf835d05ee73ee9e9ULL &&
              relaxation.digest == 0x0ecdd794c695c4c1ULL,
          "frozen internal digest mismatch");

  std::cout << "scale-34 literal-CRT Z/2 referee\n";
  std::cout << "scope=primitive-proper-AP-centred-common-scale-34-H6\n";
  std::cout << "orders=1,2,17,34 supports=" << supports.size()
            << " hereditary-words=" << words.size()
            << " labelled-contexts=" << scalar.contexts << '\n';
  std::cout << "unit-words-per-support=" << weighted_words
            << " raw-labelled-unit-states=" << raw_labelled_unit_states
            << '\n';
  print_scalar_table(tables);
  std::cout << "scalar-live-owner-histogram=";
  for (int owners = 0; owners <= WIDTH; ++owners) {
    if (owners != 0) std::cout << ',';
    std::cout << owners << ':' << scalar.owner_capacity_histogram[owners];
  }
  std::cout << '\n';
  std::cout << "scalar-survivors=" << scalar.survivors.size()
            << " survivor-supports=" << scalar.survivor_supports
            << " literal-survivor-unit-states="
            << scalar.survivor_unit_states
            << " multiplicity-profiles="
            << scalar.multiplicity_profile_histogram.size()
            << " capacity-vectors=" << scalar.capacity_vectors.size() << '\n';
  std::cout << "multiplicity-profile-histogram=";
  print_multiplicity_profiles(scalar.multiplicity_profile_histogram);
  std::cout << "support-survivor-histogram=";
  print_integer_histogram(scalar.support_survivor_histogram);
  std::cout << "anchor-bank-size-histogram=";
  print_integer_histogram(relaxation.anchor_bank_histogram);
  std::cout << "z2-bound-histogram=";
  print_integer_histogram(relaxation.bound_histogram);
  print_own_order_bounds(relaxation.own_order_bound_histogram);
  std::cout << "z2-live-owner-histogram=";
  print_integer_histogram(relaxation.live_owner_histogram);
  std::cout << "owner-obligations=" << relaxation.obligations
            << " global-ceiling=" << relaxation.global_ceiling
            << " carrier-17-ceiling="
            << relaxation.seventeen_carrier_ceiling << " threshold=" << SCALE
            << '\n';
  std::cout << "digests base=" << hex64(tables.base_digest)
            << " masks=" << hex64(tables.mask_digest)
            << " scalar-survivors=" << hex64(scalar.survivor_digest)
            << " relaxation=" << hex64(relaxation.digest) << '\n';
  std::cout << "VERDICT=EMPTY (every scalar survivor has Z/2 bound <=29<34)\n";
  return 0;
}
