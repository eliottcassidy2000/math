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

// Exact frontier scout for the primitive proper AP-centred common-scale-
// twenty Hamming-six sheet bank.  The program reconstructs the literal CRT
// masks, traverses the full hereditary divisor/unit grammar, scans all 924
// labelled supports, applies the unit-independent scalar owner-capacity gate,
// and performs exact owner-local set-union reachability on every survivor.

using Labels = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;
using StateWord = std::array<uint8_t, 6>;
using Multiplicity = std::array<uint8_t, 6>;
using CapacityVector = std::array<uint8_t, 6>;
using OwnerVector = std::array<uint8_t, 6>;

static constexpr int P = 13;
static constexpr int C = 20;
static constexpr std::array<int, 6> DIVISORS{1, 2, 4, 5, 10, 20};
static constexpr std::array<int, 20> STATE_ORDER{
    1, 2, 4, 4, 5, 5, 5, 5, 10, 10,
    10, 10, 20, 20, 20, 20, 20, 20, 20, 20};
static constexpr std::array<int, 20> STATE_UNIT{
    0, 1, 1, 3, 1, 2, 3, 4, 1, 3,
    7, 9, 1, 3, 7, 9, 11, 13, 17, 19};
static constexpr std::array<int, 7> STATE_BEGIN{0, 1, 2, 4, 8, 12, 20};
static constexpr std::array<int, 6> ORDER_REP{0, 1, 2, 4, 8, 12};

static std::array<std::array<std::array<uint32_t, 13>, 20>, 13> MASK{};
static std::array<std::array<std::array<uint8_t, 13>, 6>, 13> CARD{};
static std::array<uint32_t, 1U << C> MARK{};
static uint32_t MARK_EPOCH = 0;

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

void fnv_u64(uint64_t &digest, uint64_t value) {
  for (int shift = 0; shift < 64; shift += 8)
    fnv_byte(digest, static_cast<uint8_t>(value >> shift));
}

int centered(int value, int modulus) {
  int residue = value % modulus;
  if (residue < 0) residue += modulus;
  return 2 * residue > modulus ? residue - modulus : residue;
}

int inverse_mod_13(int value) {
  for (int candidate = 1; candidate < P; ++candidate)
    if (value * candidate % P == 1) return candidate;
  fail("nonunit modulo thirteen");
}

int crt_base(int label, int state) {
  const int order = STATE_ORDER[state];
  const int unit = STATE_UNIT[state];
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P &&
        value % order == unit % order)
      return value;
  fail("literal CRT search found no base");
}

uint32_t local_mask(int label, int state, int owner) {
  const int order = STATE_ORDER[state];
  const int base = crt_base(label, state);
  const int inverse = inverse_mod_13(owner);
  uint32_t result = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    if (-order < value && value <= order) result |= 1U << sheet;
  }
  return result;
}

int analytic_cardinality(int label, int order_index, int owner) {
  // As a sheet runs through one period modulo D, multiplication by the CRT
  // base permutes exactly the D centred integers congruent to D*(label/owner)
  // modulo 13.  The 20-sheet bank repeats that period 20/D times.
  const int order = DIVISORS[static_cast<std::size_t>(order_index)];
  const int ratio = label * inverse_mod_13(owner) % P;
  const int target = order * ratio % P;
  int period_count = 0;
  for (int value = -order + 1; value <= order; ++value) {
    int residue = value % P;
    if (residue < 0) residue += P;
    period_count += residue == target;
  }
  return (C / order) * period_count;
}

std::pair<int, int> state_range(int order_index) {
  require(0 <= order_index && order_index < 6, "bad order index");
  return {STATE_BEGIN[order_index], STATE_BEGIN[order_index + 1]};
}

bool hereditary(const OrderWord &word) {
  for (int omitted = 0; omitted < 6; ++omitted) {
    int residual = 1;
    for (int coordinate = 0; coordinate < 6; ++coordinate)
      if (coordinate != omitted)
        residual = std::lcm(residual, DIVISORS[word[coordinate]]);
    if (residual != C) return false;
  }
  return true;
}

bool hereditary_prime_power_audit(const OrderWord &word) {
  int four_providers = 0;
  int five_providers = 0;
  for (uint8_t order_index : word) {
    const int order = DIVISORS[order_index];
    four_providers += order % 4 == 0;
    five_providers += order % 5 == 0;
  }
  return four_providers >= 2 && five_providers >= 2;
}

void enumerate_orders(std::vector<OrderWord> &answer, OrderWord &word,
                      int coordinate) {
  if (coordinate == 6) {
    require(hereditary(word) == hereditary_prime_power_audit(word),
            "lcm and prime-power hereditary audits disagree");
    if (hereditary(word)) answer.push_back(word);
    return;
  }
  for (int order_index = 0; order_index < 6; ++order_index) {
    word[coordinate] = static_cast<uint8_t>(order_index);
    enumerate_orders(answer, word, coordinate + 1);
  }
}

uint64_t fibre_size(const OrderWord &word) {
  uint64_t result = 1;
  for (int order_index : word) {
    const auto [begin, end] = state_range(order_index);
    result *= static_cast<uint64_t>(end - begin);
  }
  return result;
}

void enumerate_states(const OrderWord &word, StateWord &states, int coordinate,
                      uint64_t &count, uint64_t &digest) {
  if (coordinate == 6) {
    ++count;
    for (uint8_t state : states) fnv_byte(digest, state + 1U);
    fnv_byte(digest, 0xffU);
    return;
  }
  const auto [begin, end] = state_range(word[coordinate]);
  for (int state = begin; state < end; ++state) {
    states[coordinate] = static_cast<uint8_t>(state);
    enumerate_states(word, states, coordinate + 1, count, digest);
  }
}

Multiplicity multiplicity(const OrderWord &word) {
  Multiplicity result{};
  for (int order_index : word) ++result[order_index];
  return result;
}

CapacityVector scalar_capacities(const Labels &labels,
                                 const OrderWord &word) {
  CapacityVector result{};
  for (int owner_index = 0; owner_index < 6; ++owner_index)
    for (int provider = 0; provider < 6; ++provider)
      result[owner_index] = static_cast<uint8_t>(
          result[owner_index] +
          CARD[labels[provider]][word[provider]][labels[owner_index]]);
  return result;
}

struct LocalAudit {
  bool feasible = false;
  uint8_t maximum_union = 0;
  uint32_t reachable_count = 0;
};

LocalAudit owner_local_audit(const Labels &labels, const OrderWord &word,
                             int owner) {
  std::vector<uint32_t> reachable{0};
  for (int provider = 0; provider < 6; ++provider) {
    ++MARK_EPOCH;
    require(MARK_EPOCH != 0, "owner-local marker epoch overflow");
    std::vector<uint32_t> next;
    const auto [begin, end] = state_range(word[provider]);
    for (uint32_t partial : reachable)
      for (int state = begin; state < end; ++state) {
        const uint32_t joined =
            partial | MASK[labels[provider]][state][owner];
        if (MARK[joined] != MARK_EPOCH) {
          MARK[joined] = MARK_EPOCH;
          next.push_back(joined);
        }
      }
    reachable = std::move(next);
  }
  LocalAudit result;
  result.reachable_count = static_cast<uint32_t>(reachable.size());
  for (uint32_t mask : reachable) {
    result.maximum_union = static_cast<uint8_t>(
        std::max<int>(result.maximum_union, std::popcount(mask)));
    result.feasible |= std::popcount(mask) == C;
  }
  require(result.feasible == (result.maximum_union == C),
          "owner feasibility and maximum union disagree");
  return result;
}

struct TournamentAudit {
  int ties = 0;
  int flips = 0;
  int directed_triangles = 0;
  int strongly_connected_components = 0;
  uint64_t hamiltonian_paths = 0;
  std::array<uint8_t, 6> scores{};
};

TournamentAudit tournament_audit(
    const CapacityVector &capacities,
    const std::array<LocalAudit, 6> &locals) {
  // Owner obligations are vertices.  The exact pair observable compares the
  // terminal-faithful local summary first, then scalar capacity.  Orient
  // lexicographically; an exact tie uses the support-coordinate path.
  std::array<uint8_t, 6> out{};
  TournamentAudit result;
  for (int left = 0; left < 6; ++left)
    for (int right = left + 1; right < 6; ++right) {
      int winner = left;
      const std::array<int, 3> left_key{
          locals[left].feasible, locals[left].maximum_union, capacities[left]};
      const std::array<int, 3> right_key{locals[right].feasible,
                                         locals[right].maximum_union,
                                         capacities[right]};
      if (left_key == right_key) {
        ++result.ties;
      } else if (right_key > left_key) {
        winner = right;
        ++result.flips;
      }
      const int loser = left + right - winner;
      out[winner] |= static_cast<uint8_t>(1U << loser);
      ++result.scores[winner];
    }

  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j)
      for (int k = j + 1; k < 6; ++k) {
        const bool forward = ((out[i] >> j) & 1U) &&
                             ((out[j] >> k) & 1U) &&
                             ((out[k] >> i) & 1U);
        const bool reverse = ((out[i] >> k) & 1U) &&
                             ((out[k] >> j) & 1U) &&
                             ((out[j] >> i) & 1U);
        result.directed_triangles += forward || reverse;
      }

  std::array<std::array<bool, 6>, 6> reach{};
  for (int i = 0; i < 6; ++i) {
    reach[i][i] = true;
    for (int j = 0; j < 6; ++j)
      if ((out[i] >> j) & 1U) reach[i][j] = true;
  }
  for (int k = 0; k < 6; ++k)
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j)
        reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j]);
  std::array<bool, 6> assigned{};
  for (int i = 0; i < 6; ++i)
    if (!assigned[i]) {
      ++result.strongly_connected_components;
      for (int j = 0; j < 6; ++j)
        if (reach[i][j] && reach[j][i]) assigned[j] = true;
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
  for (int last = 0; last < 6; ++last)
    result.hamiltonian_paths += paths.back()[last];
  return result;
}

Labels multiply_support(const Labels &labels, int multiplier) {
  Labels result{};
  for (int coordinate = 0; coordinate < 6; ++coordinate)
    result[coordinate] =
        static_cast<uint8_t>(multiplier * labels[coordinate] % P);
  std::sort(result.begin(), result.end());
  return result;
}

std::string histogram_string(const std::map<int, uint64_t> &histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, value] : histogram) {
    if (!first) out << ' ';
    first = false;
    out << key << ':' << value;
  }
  return out.str();
}

std::string multiplicity_string(
    const std::map<Multiplicity, uint64_t> &histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[profile, count] : histogram) {
    if (!first) out << ' ';
    first = false;
    for (int index = 0; index < 6; ++index) {
      if (index) out << ',';
      out << static_cast<int>(profile[static_cast<std::size_t>(index)]);
    }
    out << ':' << count;
  }
  return out.str();
}

int main() {
  require(STATE_ORDER.size() == STATE_UNIT.size(), "state arrays disagree");
  for (int order_index = 0; order_index < 6; ++order_index) {
    const int divisor = DIVISORS[order_index];
    const auto [begin, end] = state_range(order_index);
    require(begin == ORDER_REP[order_index], "wrong order representative");
    std::set<int> actual_units;
    for (int state = begin; state < end; ++state) {
      require(STATE_ORDER[state] == divisor, "state/order mismatch");
      require(divisor == 1 || std::gcd(STATE_UNIT[state], divisor) == 1,
              "listed state is not a unit");
      actual_units.insert(STATE_UNIT[state]);
    }
    std::set<int> expected_units;
    if (divisor == 1) {
      expected_units.insert(0);
    } else {
      for (int unit = 1; unit < divisor; ++unit)
        if (std::gcd(unit, divisor) == 1) expected_units.insert(unit);
    }
    require(actual_units == expected_units, "incomplete unit grammar");
  }

  uint64_t mask_digest = 14'695'981'039'346'656'037ULL;
  for (int label = 1; label < P; ++label)
    for (int state = 0; state < 20; ++state)
      for (int owner = 1; owner < P; ++owner) {
        MASK[label][state][owner] = local_mask(label, state, owner);
        fnv_u32(mask_digest, MASK[label][state][owner]);
      }
  for (int label = 1; label < P; ++label)
    for (int order_index = 0; order_index < 6; ++order_index)
      for (int owner = 1; owner < P; ++owner) {
        const auto [begin, end] = state_range(order_index);
        std::set<int> cardinalities;
        for (int state = begin; state < end; ++state)
          cardinalities.insert(std::popcount(MASK[label][state][owner]));
        require(cardinalities.size() == 1,
                "unit-dependent scalar cardinality");
        CARD[label][order_index][owner] =
            static_cast<uint8_t>(*cardinalities.begin());
        require(CARD[label][order_index][owner] ==
                    analytic_cardinality(label, order_index, owner),
                "literal and analytic scalar cardinalities disagree");
      }

  std::vector<OrderWord> order_words;
  OrderWord order_scratch{};
  enumerate_orders(order_words, order_scratch, 0);
  uint64_t order_digest = 14'695'981'039'346'656'037ULL;
  uint64_t weighted_states = 0;
  for (const OrderWord &word : order_words) {
    weighted_states += fibre_size(word);
    for (uint8_t order : word) fnv_byte(order_digest, order);
    fnv_byte(order_digest, 0xffU);
  }
  uint64_t literal_states = 0;
  uint64_t grammar_digest = 14'695'981'039'346'656'037ULL;
  StateWord state_scratch{};
  for (const OrderWord &word : order_words)
    enumerate_states(word, state_scratch, 0, literal_states, grammar_digest);
  require(literal_states == weighted_states,
          "literal and weighted state censuses disagree");
  require(order_words.size() == 26'961,
          "hereditary order-word census mismatch");
  require(weighted_states == 56'908'800,
          "weighted state-word census mismatch");

  uint64_t supports = 0;
  uint64_t scalar_contexts = 0;
  uint64_t scalar_digest = 14'695'981'039'346'656'037ULL;
  uint64_t capacity_digest = 14'695'981'039'346'656'037ULL;
  uint64_t owner_digest = 14'695'981'039'346'656'037ULL;
  std::set<Labels> scalar_supports;
  std::map<Labels, uint64_t> contexts_per_support;
  std::map<Multiplicity, uint64_t> scalar_patterns;
  std::map<CapacityVector, uint64_t> capacity_vectors;
  std::map<int, uint64_t> minimum_slack_histogram;
  std::map<int, uint64_t> maximum_slack_histogram;
  std::map<int, uint64_t> tight_owner_histogram;
  std::map<int, uint64_t> tournament_ties;
  std::map<int, uint64_t> tournament_flips;
  std::map<int, uint64_t> feasible_owner_histogram;
  std::map<int, uint64_t> maximum_union_histogram;
  std::map<int, uint64_t> minimum_owner_histogram;
  std::map<OwnerVector, uint64_t> owner_vectors;
  uint64_t feasible_owner_rows = 0;
  uint64_t all_owner_contexts = 0;

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
              ++supports;
              std::array<std::array<uint64_t, 6>, 6> packed{};
              for (int provider = 0; provider < 6; ++provider)
                for (int order_index = 0; order_index < 6; ++order_index)
                  for (int owner_index = 0; owner_index < 6; ++owner_index)
                    packed[provider][order_index] |=
                        static_cast<uint64_t>(CARD[labels[provider]][order_index]
                                                   [labels[owner_index]])
                        << (7 * owner_index);

              for (const OrderWord &word : order_words) {
                uint64_t capacity = 0;
                for (int provider = 0; provider < 6; ++provider)
                  capacity += packed[provider][word[provider]];
                bool scalar = true;
                for (int owner_index = 0; owner_index < 6; ++owner_index)
                  if (((capacity >> (7 * owner_index)) & 127U) < C) {
                    scalar = false;
                    break;
                  }
                if (!scalar) continue;

                const CapacityVector capacities = scalar_capacities(labels, word);
                for (int owner_index = 0; owner_index < 6; ++owner_index)
                  require(capacities[owner_index] ==
                              ((capacity >> (7 * owner_index)) & 127U),
                          "packed and direct capacities disagree");

                ++scalar_contexts;
                scalar_supports.insert(labels);
                ++contexts_per_support[labels];
                ++scalar_patterns[multiplicity(word)];
                ++capacity_vectors[capacities];
                int minimum_slack = 127;
                int maximum_slack = 0;
                int tight_owners = 0;
                for (uint8_t value : capacities) {
                  minimum_slack = std::min(minimum_slack, value - C);
                  maximum_slack = std::max(maximum_slack, value - C);
                  tight_owners += value == C;
                }
                ++minimum_slack_histogram[minimum_slack];
                ++maximum_slack_histogram[maximum_slack];
                ++tight_owner_histogram[tight_owners];

                for (uint8_t label : labels) fnv_byte(scalar_digest, label);
                for (uint8_t order : word) fnv_byte(scalar_digest, order);
                for (uint8_t value : capacities)
                  fnv_byte(capacity_digest, value);

                std::array<LocalAudit, 6> locals{};
                OwnerVector maxima{};
                int feasible_owners = 0;
                int minimum_owner = C;
                uint8_t feasible_mask = 0;
                for (int owner_index = 0; owner_index < 6; ++owner_index) {
                  locals[owner_index] = owner_local_audit(
                      labels, word, labels[owner_index]);
                  maxima[owner_index] = locals[owner_index].maximum_union;
                  feasible_owners += locals[owner_index].feasible;
                  feasible_owner_rows += locals[owner_index].feasible;
                  if (locals[owner_index].feasible)
                    feasible_mask |= static_cast<uint8_t>(1U << owner_index);
                  minimum_owner = std::min<int>(
                      minimum_owner, locals[owner_index].maximum_union);
                  ++maximum_union_histogram[locals[owner_index].maximum_union];
                }
                ++feasible_owner_histogram[feasible_owners];
                ++minimum_owner_histogram[minimum_owner];
                ++owner_vectors[maxima];
                all_owner_contexts += feasible_owners == 6;
                for (uint8_t label : labels) fnv_byte(owner_digest, label);
                for (uint8_t order : word) fnv_byte(owner_digest, order);
                fnv_byte(owner_digest, feasible_mask);
                for (const LocalAudit &audit : locals) {
                  fnv_byte(owner_digest, audit.maximum_union);
                  fnv_u32(owner_digest, audit.reachable_count);
                }

                const TournamentAudit tournament =
                    tournament_audit(capacities, locals);
                ++tournament_ties[tournament.ties];
                ++tournament_flips[tournament.flips];
                require(tournament.directed_triangles == 0,
                        "capacity tournament has a directed triangle");
                require(tournament.strongly_connected_components == 6,
                        "capacity tournament is not acyclic");
                require(tournament.hamiltonian_paths == 1,
                        "capacity tournament path count mismatch");
                std::array<uint8_t, 6> sorted_scores = tournament.scores;
                std::sort(sorted_scores.begin(), sorted_scores.end());
                require(sorted_scores ==
                            std::array<uint8_t, 6>{0, 1, 2, 3, 4, 5},
                        "capacity tournament score histogram mismatch");
              }
            }

  require(supports == 924, "support census mismatch");

  std::map<int, uint64_t> contexts_per_support_histogram;
  for (const auto &[labels, count] : contexts_per_support) {
    static_cast<void>(labels);
    ++contexts_per_support_histogram[static_cast<int>(count)];
  }

  std::set<Labels> remaining = scalar_supports;
  std::map<int, uint64_t> orbit_size_histogram;
  while (!remaining.empty()) {
    const Labels representative = *remaining.begin();
    std::set<Labels> orbit;
    for (int multiplier = 1; multiplier < P; ++multiplier)
      orbit.insert(multiply_support(representative, multiplier));
    require(std::includes(scalar_supports.begin(), scalar_supports.end(),
                          orbit.begin(), orbit.end()),
            "multiplication orbit leaves scalar support bank");
    int removed = 0;
    for (const Labels &labels : orbit) removed += remaining.erase(labels);
    require(removed == static_cast<int>(orbit.size()),
            "multiplication orbits overlap");
    ++orbit_size_histogram[static_cast<int>(orbit.size())];
  }

  uint64_t multiplicity_digest = 14'695'981'039'346'656'037ULL;
  for (const auto &[profile, count] : scalar_patterns) {
    for (uint8_t value : profile) fnv_byte(multiplicity_digest, value);
    fnv_u64(multiplicity_digest, count);
  }
  uint64_t capacity_vector_digest = 14'695'981'039'346'656'037ULL;
  for (const auto &[profile, count] : capacity_vectors) {
    for (uint8_t value : profile) fnv_byte(capacity_vector_digest, value);
    fnv_u64(capacity_vector_digest, count);
  }

  require(supports * order_words.size() == 24'911'964ULL,
          "labelled order-context census mismatch");
  require(supports * weighted_states == 52'583'731'200ULL,
          "raw labelled state-context census mismatch");
  require(mask_digest == 0x2f21825c4e17e321ULL,
          "mask-table digest mismatch");
  require(order_digest == 0x0d945578b6015e84ULL,
          "order-grammar digest mismatch");
  require(grammar_digest == 0x853c1ff8f3d58fa5ULL,
          "literal state-grammar digest mismatch");
  require(scalar_contexts == 12'584 && scalar_supports.size() == 830 &&
              scalar_patterns.size() == 65,
          "scalar classification mismatch");
  require(multiplicity_digest == 0xb60273a4a59e5d86ULL &&
              scalar_digest == 0x6141390944657ae9ULL,
          "scalar bank digest mismatch");
  require(contexts_per_support_histogram ==
              std::map<int, uint64_t>{{1, 12},  {2, 96},  {3, 24},
                                      {4, 96},  {6, 48},  {7, 96},
                                      {13, 216}, {15, 6}, {17, 24},
                                      {18, 24}, {19, 24}, {25, 24},
                                      {29, 24}, {30, 6}, {38, 48},
                                      {44, 24}, {61, 24}, {65, 12},
                                      {85, 2}},
          "contexts-per-support histogram mismatch");
  require(orbit_size_histogram ==
              std::map<int, uint64_t>{{2, 1}, {6, 2}, {12, 68}},
          "multiplication orbit-size histogram mismatch");
  require(capacity_vectors.size() == 3'941 &&
              capacity_vector_digest == 0xdeeb0167bf8b89e7ULL &&
              capacity_digest == 0x990a5b2b2ba5c94dULL,
          "capacity-vector audit mismatch");
  require(minimum_slack_histogram ==
              std::map<int, uint64_t>{{0, 11'276}, {1, 1'152}, {2, 156}},
          "minimum-slack histogram mismatch");
  require(maximum_slack_histogram ==
              std::map<int, uint64_t>{{0, 96},   {1, 1'266}, {2, 2'358},
                                      {3, 1'330}, {4, 1'336}, {5, 1'692},
                                      {6, 1'968}, {7, 474},   {8, 336},
                                      {9, 24},    {10, 1'464}, {11, 72},
                                      {16, 108},  {18, 60}},
          "maximum-slack histogram mismatch");
  require(tight_owner_histogram ==
              std::map<int, uint64_t>{{0, 1'308}, {1, 1'560}, {2, 3'384},
                                      {3, 4'052}, {4, 1'860}, {5, 324},
                                      {6, 96}},
          "tight-owner histogram mismatch");
  require(feasible_owner_rows == 3'720,
          "feasible owner-row census mismatch");
  require(feasible_owner_histogram ==
              std::map<int, uint64_t>{{0, 9'632}, {1, 2'184}, {2, 768}},
          "feasible-owner histogram mismatch");
  require(maximum_union_histogram ==
              std::map<int, uint64_t>{{12, 1'320}, {13, 2'904},
                                      {14, 9'720}, {15, 15'924},
                                      {16, 22'800}, {17, 12'876},
                                      {18, 5'700}, {19, 540}, {20, 3'720}},
          "maximum-union histogram mismatch");
  require(minimum_owner_histogram ==
              std::map<int, uint64_t>{{12, 1'194}, {13, 1'422},
                                      {14, 4'282}, {15, 3'748},
                                      {16, 1'668}, {17, 222}, {18, 48}},
          "minimum-owner histogram mismatch");
  require(owner_vectors.size() == 4'061 &&
              owner_digest == 0xdb14aac97cf1f251ULL,
          "owner-profile audit mismatch");
  require(all_owner_contexts == 0,
          "a context is locally feasible at all six owners");
  require(tournament_ties ==
              std::map<int, uint64_t>{{0, 1'248}, {1, 4'272},
                                      {2, 3'936}, {3, 1'602},
                                      {4, 1'044}, {6, 212},
                                      {7, 198}, {10, 72}},
          "tournament tie-edge histogram mismatch");
  require(tournament_flips ==
              std::map<int, uint64_t>{{0, 72},   {1, 170}, {2, 357},
                                      {3, 655},  {4, 1'208}, {5, 1'689},
                                      {6, 2'124}, {7, 2'025}, {8, 1'712},
                                      {9, 1'234}, {10, 719}, {11, 413},
                                      {12, 143}, {13, 50}, {14, 13}},
          "tournament edge-flip histogram mismatch");

  std::cout << "scale-twenty AP-centred Hamming-six scalar frontier scout\n";
  std::cout << "divisor grammar 1,2,4,5,10,20; literal states 20\n";
  std::cout << "supports " << supports << "; hereditary order words "
            << order_words.size() << "; labelled order contexts "
            << supports * order_words.size() << '\n';
  std::cout << "state words/support " << weighted_states
            << "; raw labelled states " << supports * weighted_states << '\n';
  std::cout << "mask FNV64 " << std::hex << mask_digest
            << "; order FNV64 " << order_digest
            << "; state FNV64 " << grammar_digest << std::dec << '\n';
  std::cout << "scalar contexts " << scalar_contexts << " on "
            << scalar_supports.size() << " supports; multiplicity patterns "
            << scalar_patterns.size() << "; multiplicity FNV64 " << std::hex
            << multiplicity_digest << "; scalar-bank FNV64 " << scalar_digest
            << std::dec << '\n';
  std::cout << "scalar multiplicities n1,n2,n4,n5,n10,n20 "
            << multiplicity_string(scalar_patterns) << '\n';
  std::cout << "contexts-per-support histogram "
            << histogram_string(contexts_per_support_histogram) << '\n';
  std::cout << "multiplication orbit-size histogram "
            << histogram_string(orbit_size_histogram)
            << " (telemetry; no quotient)\n";
  std::cout << "capacity vectors " << capacity_vectors.size()
            << "; ordered-vector FNV64 " << std::hex
            << capacity_vector_digest << "; row-capacity FNV64 "
            << capacity_digest << std::dec << '\n';
  std::cout << "minimum scalar-slack histogram "
            << histogram_string(minimum_slack_histogram) << '\n';
  std::cout << "maximum scalar-slack histogram "
            << histogram_string(maximum_slack_histogram) << '\n';
  std::cout << "tight-owner/context histogram "
            << histogram_string(tight_owner_histogram) << '\n';
  std::cout << "owner-local rows " << scalar_contexts * 6
            << "; feasible rows " << feasible_owner_rows << '\n';
  std::cout << "feasible-owner/context histogram "
            << histogram_string(feasible_owner_histogram) << '\n';
  std::cout << "maximum reachable sheet-union histogram "
            << histogram_string(maximum_union_histogram) << '\n';
  std::cout << "minimum owner maximum/context histogram "
            << histogram_string(minimum_owner_histogram) << '\n';
  std::cout << "distinct owner max-union vectors " << owner_vectors.size()
            << "; owner-profile FNV64 " << std::hex << owner_digest
            << std::dec << '\n';
  std::cout << "owner-local all-six contexts " << all_owner_contexts << '\n';
  std::cout << "tournament pair observable exact ordered (feasible,max-union,capacity) owner summaries; lexicographic switch and coordinate tie Hamiltonian path\n";
  std::cout << "tournament fingerprints all " << scalar_contexts
            << " transitive: scores 0,1,2,3,4,5; cycles 0; SCCs 6; Hamiltonian paths 1\n";
  std::cout << "tournament tie-edge histogram "
            << histogram_string(tournament_ties) << '\n';
  std::cout << "tournament edge-flip histogram "
            << histogram_string(tournament_flips) << '\n';
  std::cout << "challenged vertices owner obligations preserve the terminal deficit through their feasibility/max-union vector and the scalar gate through capacities; the tournament loses thresholds and magnitudes, while provider, divisor, residue, isolated sheet, and wall-event vertices lose shared-unit incidence\n";
  std::cout << "frontier verdict scalar-empty "
            << (scalar_contexts == 0 ? "yes" : "no")
            << "; owner-local all-six empty "
            << (all_owner_contexts == 0 ? "yes" : "no")
            << "; next exact obstruction "
            << (all_owner_contexts == 0
                    ? "none inside the common-scale Hamming-six face"
                    : "globally covariant unit replay across owners")
            << '\n';
  std::cout << "local D20 mask table at owner one (units 1,3,7,9,11,13,17,19; ratios 1..12 in hex)\n";
  for (int state = 12; state < 20; ++state) {
    std::cout << "  e=" << STATE_UNIT[state] << ':';
    for (int ratio = 1; ratio < P; ++ratio)
      std::cout << ' ' << std::hex << MASK[ratio][state][1] << std::dec;
    std::cout << '\n';
  }
}
