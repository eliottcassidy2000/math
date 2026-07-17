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

// Exact primary frontier scout for the primitive proper AP-centred common-
// scale-sixteen Hamming-six sheet bank.  It reconstructs literal CRT masks,
// traverses the hereditary divisor/unit grammar, scans all 924 supports, and
// applies exact owner-local union reachability.  Any context surviving all six
// owner-local gates is retained for a literal globally covariant unit replay.

using Labels = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;
using StateWord = std::array<uint8_t, 6>;
using Multiplicity = std::array<uint8_t, 5>;
using OwnerVector = std::array<uint8_t, 6>;

static constexpr int P = 13;
static constexpr int C = 16;
static constexpr uint16_t FULL = UINT16_MAX;

static constexpr std::array<int, 16> STATE_ORDER{
    1, 2, 4, 4, 8, 8, 8, 8, 16, 16, 16, 16, 16, 16, 16, 16};
static constexpr std::array<int, 16> STATE_UNIT{
    0, 1, 1, 3, 1, 3, 5, 7, 1, 3, 5, 7, 9, 11, 13, 15};
static constexpr std::array<int, 5> DIVISORS{1, 2, 4, 8, 16};
static constexpr std::array<int, 5> ORDER_REP{0, 1, 2, 4, 8};
static constexpr std::array<int, 6> STATE_BEGIN{0, 1, 2, 4, 8, 16};

static std::array<std::array<std::array<uint16_t, 13>, 16>, 13> MASK{};

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
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
  fail("literal CRT search found no solution");
}

uint16_t local_mask(int label, int state, int owner) {
  const int order = STATE_ORDER[state];
  const int base = crt_base(label, state);
  const int inverse = inverse_mod_13(owner);
  uint16_t result = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    if (-order < value && value <= order)
      result |= static_cast<uint16_t>(1U << sheet);
  }
  return result;
}

std::pair<int, int> state_range(int order_index) {
  require(0 <= order_index && order_index < 5, "bad order index");
  return {STATE_BEGIN[order_index], STATE_BEGIN[order_index + 1]};
}

bool hereditary(const OrderWord &orders) {
  for (int omitted = 0; omitted < 6; ++omitted) {
    int residual_lcm = 1;
    for (int coordinate = 0; coordinate < 6; ++coordinate)
      if (coordinate != omitted)
        residual_lcm =
            std::lcm(residual_lcm, DIVISORS[orders[coordinate]]);
    if (residual_lcm != C) return false;
  }
  return true;
}

void enumerate_orders(std::vector<OrderWord> &answer, OrderWord &word,
                      int coordinate) {
  if (coordinate == 6) {
    if (hereditary(word)) answer.push_back(word);
    return;
  }
  for (int order_index = 0; order_index < 5; ++order_index) {
    word[coordinate] = static_cast<uint8_t>(order_index);
    enumerate_orders(answer, word, coordinate + 1);
  }
}

uint64_t weighted_state_fibre(const OrderWord &orders) {
  uint64_t result = 1;
  for (int order_index : orders) {
    const auto [begin, end] = state_range(order_index);
    result *= static_cast<uint64_t>(end - begin);
  }
  return result;
}

void enumerate_state_grammar(const OrderWord &orders, StateWord &states,
                             int coordinate, uint64_t &count,
                             uint64_t &digest) {
  if (coordinate == 6) {
    ++count;
    for (uint8_t state : states) {
      digest ^= static_cast<uint64_t>(state + 1);
      digest *= 1'099'511'628'211ULL;
    }
    digest ^= 0xffULL;
    digest *= 1'099'511'628'211ULL;
    return;
  }
  const auto [begin, end] = state_range(orders[coordinate]);
  for (int state = begin; state < end; ++state) {
    states[coordinate] = static_cast<uint8_t>(state);
    enumerate_state_grammar(orders, states, coordinate + 1, count, digest);
  }
}

Multiplicity multiplicity(const OrderWord &orders) {
  Multiplicity result{};
  for (int order_index : orders) ++result[order_index];
  return result;
}

bool scalar_capacity(const Labels &labels, const OrderWord &orders, int owner) {
  int capacity = 0;
  for (int provider = 0; provider < 6; ++provider)
    capacity += std::popcount(
        MASK[labels[provider]][ORDER_REP[orders[provider]]][owner]);
  return capacity >= C;
}

struct LocalAudit {
  bool feasible = false;
  uint8_t maximum_union = 0;
  uint32_t reachable_count = 0;
};

LocalAudit owner_local_audit(const Labels &labels, const OrderWord &orders,
                             int owner) {
  std::vector<uint16_t> reachable{0};
  std::array<uint8_t, 1U << C> seen{};
  for (int provider = 0; provider < 6; ++provider) {
    std::vector<uint16_t> next;
    const auto [begin, end] = state_range(orders[provider]);
    for (uint16_t partial : reachable)
      for (int state = begin; state < end; ++state) {
        const uint16_t joined =
            partial | MASK[labels[provider]][state][owner];
        if (!seen[joined]) {
          seen[joined] = 1;
          next.push_back(joined);
        }
      }
    for (uint16_t value : next) seen[value] = 0;
    reachable = std::move(next);
  }
  LocalAudit result;
  result.reachable_count = static_cast<uint32_t>(reachable.size());
  for (uint16_t mask : reachable) {
    result.maximum_union = static_cast<uint8_t>(
        std::max<int>(result.maximum_union, std::popcount(mask)));
    result.feasible |= mask == FULL;
  }
  require(result.feasible == (result.maximum_union == C),
          "feasibility and maximum union disagree");
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

TournamentAudit tournament_audit(const std::array<LocalAudit, 6> &locals) {
  // Pair observable: compare (feasible,max-union).  The switch is the sign of
  // this lexicographic comparison; exact ties use coordinate 0->...->5.
  std::array<uint8_t, 6> out{};
  TournamentAudit result;
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j) {
      const std::pair<bool, uint8_t> left{locals[i].feasible,
                                          locals[i].maximum_union};
      const std::pair<bool, uint8_t> right{locals[j].feasible,
                                           locals[j].maximum_union};
      int winner = i;
      if (left == right) {
        ++result.ties;
      } else if (right > left) {
        winner = j;
        ++result.flips;
      }
      const int loser = i + j - winner;
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
    result.hamiltonian_paths += paths[(1U << 6) - 1U][last];
  return result;
}

struct Context {
  Labels labels{};
  OrderWord orders{};
  uint64_t fibre_size = 0;
};

void replay_context(const Context &context, StateWord &states, int coordinate,
                    uint64_t &visited, uint64_t &survivors,
                    uint64_t &survivor_digest) {
  if (coordinate == 6) {
    ++visited;
    for (int owner : context.labels) {
      uint16_t covered = 0;
      for (int provider = 0; provider < 6; ++provider)
        covered |= MASK[context.labels[provider]][states[provider]][owner];
      if (covered != FULL) return;
    }
    ++survivors;
    for (uint8_t state : states) {
      survivor_digest ^= static_cast<uint64_t>(state + 1);
      survivor_digest *= 1'099'511'628'211ULL;
    }
    return;
  }
  const auto [begin, end] = state_range(context.orders[coordinate]);
  for (int state = begin; state < end; ++state) {
    states[coordinate] = static_cast<uint8_t>(state);
    replay_context(context, states, coordinate + 1, visited, survivors,
                   survivor_digest);
  }
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
    for (int index = 0; index < 5; ++index) {
      if (index) out << ',';
      out << static_cast<int>(profile[index]);
    }
    out << ':' << count;
  }
  return out.str();
}

int main() {
  require(STATE_ORDER.size() == STATE_UNIT.size(), "state arrays disagree");
  for (int order_index = 0; order_index < 5; ++order_index) {
    const int divisor = DIVISORS[order_index];
    require(C % divisor == 0, "listed order does not divide sixteen");
    const auto [begin, end] = state_range(order_index);
    require(begin == ORDER_REP[order_index], "wrong order representative");
    std::set<int> actual_units;
    for (int state = begin; state < end; ++state) {
      require(STATE_ORDER[state] == divisor,
              "state assigned to wrong effective order");
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
    require(actual_units == expected_units, "incomplete divisor/unit grammar");
  }

  for (int label = 1; label < P; ++label)
    for (int state = 0; state < 16; ++state)
      for (int owner = 1; owner < P; ++owner)
        MASK[label][state][owner] = local_mask(label, state, owner);

  for (int order_index = 0; order_index < 5; ++order_index)
    for (int provider = 1; provider < P; ++provider)
      for (int owner = 1; owner < P; ++owner) {
        std::set<int> cardinalities;
        const auto [begin, end] = state_range(order_index);
        for (int state = begin; state < end; ++state)
          cardinalities.insert(
              std::popcount(MASK[provider][state][owner]));
        require(cardinalities.size() == 1,
                "mask cardinality depends on unit");
      }

  std::vector<OrderWord> order_words;
  OrderWord order_scratch{};
  enumerate_orders(order_words, order_scratch, 0);
  uint64_t weighted_states_per_support = 0;
  for (const OrderWord &orders : order_words)
    weighted_states_per_support += weighted_state_fibre(orders);
  uint64_t literal_states_per_support = 0;
  uint64_t grammar_digest = 14'695'981'039'346'656'037ULL;
  StateWord state_scratch{};
  for (const OrderWord &orders : order_words)
    enumerate_state_grammar(orders, state_scratch, 0,
                            literal_states_per_support, grammar_digest);
  require(order_words.size() == 5'385, "hereditary order-word mismatch");
  require(weighted_states_per_support == 14'942'208,
          "weighted state-word mismatch");
  require(literal_states_per_support == weighted_states_per_support,
          "literal and weighted state counts disagree");
  require(grammar_digest == 0xa25ff24dc14911c5ULL,
          "literal state-grammar digest mismatch");

  uint64_t supports = 0;
  uint64_t scalar_contexts = 0;
  uint64_t feasible_owner_rows = 0;
  uint64_t all_owner_local_contexts = 0;
  uint64_t owner_vector_digest = 14'695'981'039'346'656'037ULL;
  std::set<Labels> scalar_supports;
  std::set<Labels> owner_local_supports;
  std::map<Multiplicity, uint64_t> scalar_patterns;
  std::map<Multiplicity, uint64_t> owner_local_patterns;
  std::map<Labels, uint64_t> contexts_per_support;
  std::map<int, uint64_t> feasible_owner_histogram;
  std::map<int, uint64_t> maximum_union_histogram;
  std::map<int, uint64_t> minimum_owner_histogram;
  std::map<int, uint64_t> tournament_tie_histogram;
  std::map<int, uint64_t> tournament_flip_histogram;
  std::map<OwnerVector, uint64_t> owner_vectors;
  std::vector<Context> local_survivors;

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
              for (const OrderWord &orders : order_words) {
                bool scalar = true;
                for (int owner : labels)
                  if (!scalar_capacity(labels, orders, owner)) {
                    scalar = false;
                    break;
                  }
                if (!scalar) continue;
                ++scalar_contexts;
                scalar_supports.insert(labels);
                ++contexts_per_support[labels];
                ++scalar_patterns[multiplicity(orders)];

                int feasible_owners = 0;
                int minimum_owner = C;
                OwnerVector owner_vector{};
                std::array<LocalAudit, 6> locals{};
                for (int owner_index = 0; owner_index < 6; ++owner_index) {
                  locals[owner_index] = owner_local_audit(
                      labels, orders, labels[owner_index]);
                  const LocalAudit &audit = locals[owner_index];
                  owner_vector[owner_index] = audit.maximum_union;
                  feasible_owners += audit.feasible;
                  feasible_owner_rows += audit.feasible;
                  minimum_owner =
                      std::min<int>(minimum_owner, audit.maximum_union);
                  ++maximum_union_histogram[audit.maximum_union];
                }
                const TournamentAudit tournament = tournament_audit(locals);
                ++tournament_tie_histogram[tournament.ties];
                ++tournament_flip_histogram[tournament.flips];
                require(tournament.directed_triangles == 0,
                        "owner-summary tournament has a directed triangle");
                require(tournament.strongly_connected_components == 6,
                        "owner-summary tournament is not acyclic");
                require(tournament.hamiltonian_paths == 1,
                        "owner-summary tournament lacks a unique path");
                std::array<uint8_t, 6> sorted_scores = tournament.scores;
                std::sort(sorted_scores.begin(), sorted_scores.end());
                require(sorted_scores ==
                            std::array<uint8_t, 6>{0, 1, 2, 3, 4, 5},
                        "owner-summary score histogram changed");
                ++feasible_owner_histogram[feasible_owners];
                ++minimum_owner_histogram[minimum_owner];
                ++owner_vectors[owner_vector];
                for (uint8_t value : owner_vector) {
                  owner_vector_digest ^= static_cast<uint64_t>(value + 1);
                  owner_vector_digest *= 1'099'511'628'211ULL;
                }
                if (feasible_owners == 6) {
                  ++all_owner_local_contexts;
                  owner_local_supports.insert(labels);
                  ++owner_local_patterns[multiplicity(orders)];
                  local_survivors.push_back(
                      {labels, orders, weighted_state_fibre(orders)});
                }
              }
            }

  require(supports == 924, "support census mismatch");
  require(supports * order_words.size() == 4'975'740,
          "labelled order-context census mismatch");
  require(supports * weighted_states_per_support == 13'806'600'192ULL,
          "raw labelled state-context census mismatch");

  uint64_t local_fibre_words = 0;
  for (const Context &context : local_survivors)
    local_fibre_words += context.fibre_size;

  constexpr uint64_t REPLAY_LIMIT = 500'000'000ULL;
  uint64_t replayed_words = 0;
  uint64_t global_survivors = 0;
  uint64_t global_survivor_digest = 14'695'981'039'346'656'037ULL;
  bool replay_complete = local_fibre_words <= REPLAY_LIMIT;
  if (replay_complete)
    for (const Context &context : local_survivors)
      replay_context(context, state_scratch, 0, replayed_words,
                     global_survivors, global_survivor_digest);
  require(!replay_complete || replayed_words == local_fibre_words,
          "literal global replay did not traverse its full bank");

  std::set<Labels> remaining_supports = scalar_supports;
  std::map<int, uint64_t> scalar_orbit_histogram;
  while (!remaining_supports.empty()) {
    const Labels representative = *remaining_supports.begin();
    std::set<Labels> orbit;
    for (int multiplier = 1; multiplier < P; ++multiplier)
      orbit.insert(multiply_support(representative, multiplier));
    require(std::includes(scalar_supports.begin(), scalar_supports.end(),
                          orbit.begin(), orbit.end()),
            "multiplication orbit leaves scalar support bank");
    int removed = 0;
    for (const Labels &labels : orbit) removed += remaining_supports.erase(labels);
    require(removed == static_cast<int>(orbit.size()),
            "scalar multiplication orbits overlap");
    ++scalar_orbit_histogram[static_cast<int>(orbit.size())];
  }

  std::map<int, uint64_t> contexts_per_support_histogram;
  for (const auto &[labels, count] : contexts_per_support) {
    static_cast<void>(labels);
    ++contexts_per_support_histogram[static_cast<int>(count)];
  }

  const std::map<Multiplicity, uint64_t> expected_scalar_patterns{
      {{{0, 0, 0, 3, 3}}, 28},  {{{0, 0, 0, 4, 2}}, 54},
      {{{0, 0, 1, 0, 5}}, 48},  {{{0, 0, 1, 1, 4}}, 216},
      {{{0, 0, 1, 2, 3}}, 360}, {{{0, 0, 1, 3, 2}}, 444},
      {{{0, 0, 2, 2, 2}}, 132}, {{{0, 0, 3, 0, 3}}, 16},
      {{{0, 0, 3, 1, 2}}, 192}, {{{0, 0, 4, 0, 2}}, 30},
      {{{0, 1, 1, 2, 2}}, 48},  {{{0, 1, 2, 0, 3}}, 120},
      {{{0, 1, 2, 1, 2}}, 120}, {{{0, 1, 3, 0, 2}}, 96},
      {{{0, 2, 0, 0, 4}}, 36},  {{{0, 2, 0, 1, 3}}, 144},
      {{{0, 2, 0, 2, 2}}, 216}, {{{0, 2, 1, 0, 3}}, 48},
      {{{0, 2, 1, 1, 2}}, 48},  {{{0, 3, 1, 0, 2}}, 144}};
  require(scalar_contexts == 2'540 && scalar_supports.size() == 404,
          "scalar classification mismatch");
  require(scalar_patterns == expected_scalar_patterns,
          "scalar multiplicity classification mismatch");
  require(contexts_per_support_histogram ==
              std::map<int, uint64_t>{{1, 102}, {2, 96}, {3, 24},
                                      {4, 24},  {7, 24}, {8, 24},
                                      {11, 36}, {13, 24}, {16, 24},
                                      {17, 24}, {109, 2}},
          "contexts-per-support histogram mismatch");
  require(scalar_orbit_histogram ==
              std::map<int, uint64_t>{{2, 1}, {6, 1}, {12, 33}},
          "scalar support orbit histogram mismatch");
  require(feasible_owner_rows == 636,
          "feasible owner-row census mismatch");
  require(feasible_owner_histogram ==
              std::map<int, uint64_t>{{0, 2'006}, {1, 432}, {2, 102}},
          "feasible-owner/context histogram mismatch");
  require(maximum_union_histogram ==
              std::map<int, uint64_t>{{9, 144},  {10, 468}, {11, 876},
                                      {12, 2'316}, {13, 4'068},
                                      {14, 4'740}, {15, 1'992}, {16, 636}},
          "maximum-union histogram mismatch");
  require(minimum_owner_histogram ==
              std::map<int, uint64_t>{{9, 120}, {10, 462}, {11, 336},
                                      {12, 720}, {13, 886}, {14, 16}},
          "minimum-owner histogram mismatch");
  require(owner_vectors.size() == 1'210 &&
              owner_vector_digest == 0x1be49d1ad3c21515ULL,
          "owner max-union vector audit mismatch");
  require(all_owner_local_contexts == 0 && owner_local_supports.empty() &&
              owner_local_patterns.empty() && local_survivors.empty() &&
              local_fibre_words == 0,
          "owner-local bank unexpectedly survives");
  require(replay_complete && replayed_words == 0 && global_survivors == 0 &&
              global_survivor_digest == 0xcbf29ce484222325ULL,
          "empty global replay audit mismatch");
  require(tournament_tie_histogram ==
              std::map<int, uint64_t>{{1, 120}, {2, 648}, {3, 390},
                                      {4, 744}, {6, 344}, {7, 270},
                                      {10, 24}},
          "tournament tie-edge histogram mismatch");
  require(tournament_flip_histogram ==
              std::map<int, uint64_t>{{0, 47},  {1, 43},  {2, 115},
                                      {3, 209}, {4, 387}, {5, 384},
                                      {6, 510}, {7, 360}, {8, 256},
                                      {9, 151}, {10, 48}, {11, 25},
                                      {12, 4}, {14, 1}},
          "tournament edge-flip histogram mismatch");

  std::cout << "scale-sixteen AP-centred Hamming-six frontier scout\n";
  std::cout << "divisor grammar 1,2,4,8,16; literal states 16\n";
  std::cout << "supports " << supports << '\n';
  std::cout << "hereditary order words " << order_words.size() << '\n';
  std::cout << "hereditary labelled order contexts "
            << supports * order_words.size() << '\n';
  std::cout << "hereditary state words per support "
            << weighted_states_per_support << '\n';
  std::cout << "literal state-grammar FNV64 " << std::hex << grammar_digest
            << std::dec << '\n';
  std::cout << "hereditary labelled state contexts "
            << supports * weighted_states_per_support << '\n';
  std::cout << "scalar contexts " << scalar_contexts << " on "
            << scalar_supports.size() << " supports; patterns "
            << scalar_patterns.size() << '\n';
  std::cout << "scalar multiplicities n1,n2,n4,n8,n16 "
            << multiplicity_string(scalar_patterns) << '\n';
  std::cout << "contexts-per-support histogram "
            << histogram_string(contexts_per_support_histogram) << '\n';
  std::cout << "scalar multiplication orbit-size histogram "
            << histogram_string(scalar_orbit_histogram)
            << " (telemetry; no quotient)\n";
  std::cout << "owner-local rows " << scalar_contexts * 6
            << "; feasible rows " << feasible_owner_rows << '\n';
  std::cout << "feasible-owner/context histogram "
            << histogram_string(feasible_owner_histogram) << '\n';
  std::cout << "maximum reachable sheet-union histogram "
            << histogram_string(maximum_union_histogram) << '\n';
  std::cout << "minimum owner maximum/context histogram "
            << histogram_string(minimum_owner_histogram) << '\n';
  std::cout << "distinct owner max-union vectors " << owner_vectors.size()
            << "; vector FNV64 " << std::hex << owner_vector_digest
            << std::dec << '\n';
  std::cout << "owner-local contexts " << all_owner_local_contexts << " on "
            << owner_local_supports.size() << " supports; patterns "
            << owner_local_patterns.size() << '\n';
  if (owner_local_patterns.empty()) {
    std::cout << "owner-local multiplicities none\n";
  } else {
    std::cout << "owner-local multiplicities "
              << multiplicity_string(owner_local_patterns) << '\n';
  }
  std::cout << "owner-local literal unit words " << local_fibre_words << '\n';
  if (replay_complete)
    std::cout << "global replay complete " << replayed_words
              << "; common-sheet survivors " << global_survivors
              << "; survivor FNV64 " << std::hex << global_survivor_digest
              << std::dec << '\n';
  else
    std::cout << "global replay deferred: " << local_fibre_words
              << " words exceed scout limit " << REPLAY_LIMIT << '\n';
  std::cout << "tournament pair observable exact ordered owner summaries ((feasible,max_i),(feasible,max_j)); all pair data reconstruct the terminal-faithful vector; switch is lexicographic sign with fixed support-coordinate tie path\n";
  std::cout << "tournament fingerprints all " << scalar_contexts
            << " transitive: score histogram 0,1,2,3,4,5; directed triangles 0; SCCs 6; Hamiltonian paths 1\n";
  std::cout << "tournament tie-edge histogram "
            << histogram_string(tournament_tie_histogram) << '\n';
  std::cout << "tournament edge-flip histogram against tie path "
            << histogram_string(tournament_flip_histogram) << '\n';
  std::cout << "challenged vertices owner obligations preserve the terminal obstruction; sheet positions are faithful only inside one owner DP, while providers, residues, and divisor words destroy shared-unit incidence; the oriented tournament alone also loses deficit magnitudes\n";
  std::cout << "local D16 mask table at owner one (units 1,3,5,7,9,11,13,15; ratios 1..12 in hex)\n";
  for (int state = 8; state < 16; ++state) {
    std::cout << "  e=" << STATE_UNIT[state] << ':';
    for (int ratio = 1; ratio < P; ++ratio)
      std::cout << ' ' << std::hex << MASK[ratio][state][1] << std::dec;
    std::cout << '\n';
  }
}
