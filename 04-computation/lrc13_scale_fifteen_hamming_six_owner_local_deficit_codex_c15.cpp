#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Exact primary certificate for the primitive proper AP-centred common-scale-
// fifteen Hamming-six sheet bank.  The program reconstructs every literal CRT
// mask, enumerates the complete hereditary divisor/unit grammar, and scans all
// 924 labelled supports.  Scalar capacity leaves 2,184 divisor contexts.
// Exact owner-local union reachability shows that no context is feasible at
// all six owners, so no global unit fibre remains to inspect.
//
// The terminal-faithful object is the six-entry vector of exact owner-local
// maximum unions (equivalently, together with its entries equal to fifteen,
// the feasible-owner subset).  A tournament on owner obligations is reported
// only as structural telemetry: compare (feasible,max-union), then orient
// exact ties along the fixed support-coordinate Hamiltonian path.  This
// quotient is transitive and deliberately is not used as the proof step.

using Labels = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;
using StateWord = std::array<uint8_t, 6>;
using Multiplicity = std::array<uint8_t, 4>;
using OwnerVector = std::array<uint8_t, 6>;

static constexpr int P = 13;
static constexpr int C = 15;
static constexpr uint16_t FULL = (1U << C) - 1U;

// Each effective order D carries precisely the units modulo D (with the
// unique zero representative at D=1).
static constexpr std::array<int, 15> STATE_ORDER{
    1, 3, 3, 5, 5, 5, 5, 15, 15, 15, 15, 15, 15, 15, 15};
static constexpr std::array<int, 15> STATE_UNIT{
    0, 1, 2, 1, 2, 3, 4, 1, 2, 4, 7, 8, 11, 13, 14};
static constexpr std::array<int, 4> DIVISORS{1, 3, 5, 15};
static constexpr std::array<int, 4> ORDER_REP{0, 1, 3, 7};
static constexpr std::array<int, 5> STATE_BEGIN{0, 1, 3, 7, 15};

static std::array<std::array<std::array<uint16_t, 13>, 15>, 13> MASK{};

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

// Literal CRT search, intentionally independent of any closed mask formula.
int crt_base(int label, int state) {
  const int order = STATE_ORDER[state];
  const int unit = STATE_UNIT[state];
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P &&
        value % order == unit % order)
      return value;
  fail("CRT base missing");
}

uint16_t local_mask(int label, int state, int owner) {
  const int order = STATE_ORDER[state];
  const int base = crt_base(label, state);
  const int inverse = inverse_mod_13(owner);
  uint16_t answer = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    // Strict-safe interval (-D,D] on the fifteen AP sheets.
    if (-order < value && value <= order)
      answer |= static_cast<uint16_t>(1U << sheet);
  }
  return answer;
}

std::pair<int, int> state_range(int order_index) {
  require(0 <= order_index && order_index < 4, "bad order index");
  return {STATE_BEGIN[order_index], STATE_BEGIN[order_index + 1]};
}

bool hereditary(const OrderWord &orders) {
  // Removing any one provider must leave lcm fifteen.
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
  for (int order_index = 0; order_index < 4; ++order_index) {
    word[coordinate] = static_cast<uint8_t>(order_index);
    enumerate_orders(answer, word, coordinate + 1);
  }
}

uint64_t weighted_state_fibre(const OrderWord &orders) {
  uint64_t answer = 1;
  for (int order_index : orders) {
    const auto [begin, end] = state_range(order_index);
    answer *= static_cast<uint64_t>(end - begin);
  }
  return answer;
}

void enumerate_states_for_order(const OrderWord &orders, StateWord &states,
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
    enumerate_states_for_order(orders, states, coordinate + 1, count, digest);
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
  uint16_t reachable_count = 0;
};

LocalAudit owner_local_audit(const Labels &labels, const OrderWord &orders,
                             int owner) {
  // Exact union DP: each provider chooses one unit independently for this
  // owner.  No multiplication quotient or covariance shortcut is used.
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
  require(reachable.size() <= UINT16_MAX,
          "owner reachable-mask count overflows telemetry field");
  result.reachable_count = static_cast<uint16_t>(reachable.size());
  for (uint16_t mask : reachable) {
    const int size = std::popcount(mask);
    result.maximum_union =
        static_cast<uint8_t>(std::max<int>(result.maximum_union, size));
    result.feasible |= mask == FULL;
  }
  require(result.feasible == (result.maximum_union == C),
          "full-mask feasibility and maximum union disagree");
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
  // Owner i beats owner j when (feasible,max-union) is lexicographically
  // larger.  Exact ties follow the fixed coordinate path 0->1->...->5.
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
        const bool cycle_forward = ((out[i] >> j) & 1U) &&
                                   ((out[j] >> k) & 1U) &&
                                   ((out[k] >> i) & 1U);
        const bool cycle_reverse = ((out[i] >> k) & 1U) &&
                                   ((out[k] >> j) & 1U) &&
                                   ((out[j] >> i) & 1U);
        result.directed_triangles += cycle_forward || cycle_reverse;
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
        if (previous_mask == 0) continue;
        for (int previous = 0; previous < 6; ++previous)
          if (((previous_mask >> previous) & 1) &&
              ((out[previous] >> last) & 1U))
            paths[mask][last] += paths[previous_mask][previous];
      }
  for (int last = 0; last < 6; ++last)
    result.hamiltonian_paths += paths[(1U << 6) - 1U][last];
  return result;
}

Labels multiply_support(const Labels &labels, int multiplier) {
  Labels answer{};
  for (int coordinate = 0; coordinate < 6; ++coordinate)
    answer[coordinate] =
        static_cast<uint8_t>(multiplier * labels[coordinate] % P);
  std::sort(answer.begin(), answer.end());
  return answer;
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
    out << static_cast<int>(profile[0]) << ','
        << static_cast<int>(profile[1]) << ','
        << static_cast<int>(profile[2]) << ','
        << static_cast<int>(profile[3]) << ':' << count;
  }
  return out.str();
}

int main() {
  require(STATE_ORDER.size() == STATE_UNIT.size(), "state arrays disagree");
  for (int order_index = 0; order_index < 4; ++order_index) {
    const int divisor = DIVISORS[order_index];
    require(C % divisor == 0, "listed order does not divide fifteen");
    const auto [begin, end] = state_range(order_index);
    require(begin == ORDER_REP[order_index], "wrong order representative");
    std::set<int> actual_units;
    for (int state = begin; state < end; ++state) {
      require(STATE_ORDER[state] == divisor,
              "state assigned to the wrong order");
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
    for (int state = 0; state < 15; ++state)
      for (int owner = 1; owner < P; ++owner)
        MASK[label][state][owner] = local_mask(label, state, owner);

  // Scalar capacity uses an order representative only after checking that
  // literal mask cardinality is independent of the unit at every ratio.
  for (int order_index = 0; order_index < 4; ++order_index)
    for (int provider = 1; provider < P; ++provider)
      for (int owner = 1; owner < P; ++owner) {
        const auto [begin, end] = state_range(order_index);
        std::set<int> sizes;
        for (int state = begin; state < end; ++state)
          sizes.insert(std::popcount(MASK[provider][state][owner]));
        require(sizes.size() == 1,
                "mask cardinality unexpectedly depends on the unit");
      }

  std::vector<OrderWord> order_words;
  OrderWord order_scratch{};
  enumerate_orders(order_words, order_scratch, 0);
  uint64_t weighted_states_per_support = 0;
  for (const OrderWord &orders : order_words)
    weighted_states_per_support += weighted_state_fibre(orders);

  // Traverse all 11,169,600 literal state words, separately from the product
  // of totients used in the weighted count.
  uint64_t literal_states_per_support = 0;
  uint64_t state_grammar_digest = 14'695'981'039'346'656'037ULL;
  StateWord state_scratch{};
  for (const OrderWord &orders : order_words)
    enumerate_states_for_order(orders, state_scratch, 0,
                               literal_states_per_support,
                               state_grammar_digest);

  require(order_words.size() == 3'249,
          "hereditary divisor-word census mismatch");
  require(weighted_states_per_support == 11'169'600,
          "weighted hereditary state-word census mismatch");
  require(literal_states_per_support == weighted_states_per_support,
          "literal and weighted state-word censuses disagree");
  require(state_grammar_digest == 0x555244fbea49f335ULL,
          "literal state-grammar digest mismatch");

  uint64_t supports = 0;
  uint64_t scalar_contexts = 0;
  uint64_t feasible_owner_rows = 0;
  uint64_t all_owner_contexts = 0;
  uint64_t owner_vector_digest = 14'695'981'039'346'656'037ULL;
  std::set<Labels> scalar_supports;
  std::map<Labels, uint64_t> contexts_per_support;
  std::map<Multiplicity, uint64_t> scalar_patterns;
  std::map<int, uint64_t> feasible_owner_histogram;
  std::map<int, uint64_t> maximum_union_histogram;
  std::map<int, uint64_t> context_minimum_histogram;
  std::map<int, uint64_t> tournament_tie_histogram;
  std::map<int, uint64_t> tournament_flip_histogram;
  std::map<OwnerVector, uint64_t> owner_vectors;

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

                std::array<LocalAudit, 6> locals{};
                OwnerVector owner_vector{};
                int feasible_owners = 0;
                int context_minimum = C;
                for (int owner_index = 0; owner_index < 6; ++owner_index) {
                  locals[owner_index] = owner_local_audit(
                      labels, orders, labels[owner_index]);
                  owner_vector[owner_index] =
                      locals[owner_index].maximum_union;
                  feasible_owners += locals[owner_index].feasible;
                  feasible_owner_rows += locals[owner_index].feasible;
                  context_minimum = std::min<int>(
                      context_minimum, locals[owner_index].maximum_union);
                  ++maximum_union_histogram[locals[owner_index].maximum_union];
                }
                ++feasible_owner_histogram[feasible_owners];
                ++context_minimum_histogram[context_minimum];
                ++owner_vectors[owner_vector];
                all_owner_contexts += feasible_owners == 6;
                for (uint8_t value : owner_vector) {
                  owner_vector_digest ^= static_cast<uint64_t>(value + 1);
                  owner_vector_digest *= 1'099'511'628'211ULL;
                }
                owner_vector_digest ^= static_cast<uint64_t>(feasible_owners);
                owner_vector_digest *= 1'099'511'628'211ULL;

                const TournamentAudit tournament = tournament_audit(locals);
                ++tournament_tie_histogram[tournament.ties];
                ++tournament_flip_histogram[tournament.flips];
                require(tournament.directed_triangles == 0,
                        "owner-summary tournament has a directed triangle");
                require(tournament.strongly_connected_components == 6,
                        "owner-summary tournament is not acyclic");
                require(tournament.hamiltonian_paths == 1,
                        "owner-summary tournament lacks its unique path");
                std::array<uint8_t, 6> sorted_scores = tournament.scores;
                std::sort(sorted_scores.begin(), sorted_scores.end());
                require(sorted_scores == std::array<uint8_t, 6>{0, 1, 2, 3, 4, 5},
                        "owner-summary tournament score histogram changed");
              }
            }

  require(supports == 924, "support census mismatch");
  require(supports * order_words.size() == 3'002'076,
          "labelled divisor-context census mismatch");
  require(supports * weighted_states_per_support == 10'320'710'400ULL,
          "labelled raw state-context census mismatch");
  require(scalar_contexts == 2'184 && scalar_supports.size() == 462,
          "scalar-capacity classification mismatch");
  const std::map<Multiplicity, uint64_t> expected_scalar_patterns{
      {{{0, 0, 2, 4}}, 12},  {{{0, 1, 3, 2}}, 48},
      {{{0, 1, 4, 1}}, 84},  {{{0, 2, 1, 3}}, 144},
      {{{0, 2, 2, 2}}, 276}, {{{0, 2, 3, 1}}, 24},
      {{{0, 2, 4, 0}}, 6},   {{{0, 3, 0, 3}}, 12},
      {{{0, 3, 1, 2}}, 36},  {{{0, 3, 2, 1}}, 36},
      {{{0, 3, 3, 0}}, 12},  {{{0, 4, 0, 2}}, 18},
      {{{0, 4, 1, 1}}, 792}, {{{0, 4, 2, 0}}, 180},
      {{{1, 3, 1, 1}}, 336}, {{{1, 3, 2, 0}}, 168}};
  require(scalar_patterns == expected_scalar_patterns,
          "scalar multiplicity-pattern classification mismatch");
  require(feasible_owner_histogram ==
              std::map<int, uint64_t>{{0, 750}, {1, 456}, {2, 912}, {4, 66}},
          "feasible-owner histogram mismatch");
  require(maximum_union_histogram ==
              std::map<int, uint64_t>{{11, 804}, {12, 7'512}, {13, 1'812},
                                      {14, 432}, {15, 2'544}},
          "maximum-union histogram mismatch");
  require(feasible_owner_rows == 2'544,
          "feasible owner-row census mismatch");
  require(all_owner_contexts == 0,
          "a scalar context is locally feasible at all six owners");
  require(context_minimum_histogram ==
              std::map<int, uint64_t>{{11, 474}, {12, 1'692}, {13, 18}},
          "minimum-owner maximum histogram mismatch");
  require(owner_vectors.size() == 333,
          "distinct owner max-union vector census mismatch");
  require(owner_vector_digest == 0x7ef173744757e5d1ULL,
          "owner max-union vector digest mismatch");
  require(tournament_tie_histogram ==
              std::map<int, uint64_t>{{2, 96},   {3, 96},  {4, 360},
                                      {6, 144},  {7, 1'098}, {10, 288},
                                      {15, 102}},
          "owner-summary tournament tie histogram mismatch");
  require(tournament_flip_histogram ==
              std::map<int, uint64_t>{{0, 231}, {1, 132}, {2, 198},
                                      {3, 246}, {4, 432}, {5, 319},
                                      {6, 272}, {7, 172}, {8, 145},
                                      {9, 22},  {10, 9},  {11, 5},
                                      {12, 1}},
          "owner-summary tournament edge-flip histogram mismatch");

  // Multiplication by F_13^* is telemetry only; the full scan above did not
  // quotient by it.  Verify that it really preserves this support bank and
  // census its exact orbit sizes.
  std::set<Labels> remaining_supports = scalar_supports;
  std::map<int, uint64_t> orbit_size_histogram;
  while (!remaining_supports.empty()) {
    const Labels representative = *remaining_supports.begin();
    std::set<Labels> orbit;
    for (int multiplier = 1; multiplier < P; ++multiplier)
      orbit.insert(multiply_support(representative, multiplier));
    require(std::includes(scalar_supports.begin(), scalar_supports.end(),
                          orbit.begin(), orbit.end()),
            "multiplication orbit leaves the scalar support bank");
    int removed = 0;
    for (const Labels &labels : orbit) removed += remaining_supports.erase(labels);
    require(removed == static_cast<int>(orbit.size()),
            "multiplication support orbits overlap");
    ++orbit_size_histogram[static_cast<int>(orbit.size())];
  }

  std::map<int, uint64_t> contexts_per_support_histogram;
  for (const auto &[labels, count] : contexts_per_support) {
    static_cast<void>(labels);
    ++contexts_per_support_histogram[static_cast<int>(count)];
  }
  require(contexts_per_support_histogram ==
              std::map<int, uint64_t>{{1, 192}, {2, 120}, {3, 24},
                                      {4, 24},  {6, 24},  {8, 24},
                                      {11, 24}, {18, 12}, {37, 12},
                                      {54, 6}},
          "contexts-per-support histogram mismatch");
  require(orbit_size_histogram ==
              std::map<int, uint64_t>{{6, 1}, {12, 38}},
          "multiplication orbit-size histogram mismatch");

  std::cout << "scale-fifteen AP-centred Hamming-six owner-local deficit\n";
  std::cout << "divisor grammar 1,3,5,15; literal states 15\n";
  std::cout << "states (D,e) (1,0) (3,1) (3,2) (5,1) (5,2) (5,3) (5,4) (15,1) (15,2) (15,4) (15,7) (15,8) (15,11) (15,13) (15,14)\n";
  std::cout << "supports " << supports << '\n';
  std::cout << "hereditary order words " << order_words.size() << '\n';
  std::cout << "hereditary labelled order contexts "
            << supports * order_words.size() << '\n';
  std::cout << "hereditary state words per support "
            << weighted_states_per_support << '\n';
  std::cout << "literal state-grammar FNV64 " << std::hex
            << state_grammar_digest << std::dec << '\n';
  std::cout << "hereditary labelled state contexts "
            << supports * weighted_states_per_support << '\n';
  std::cout << "scalar-capacity contexts " << scalar_contexts << " on "
            << scalar_supports.size() << " supports; multiplicity patterns "
            << scalar_patterns.size() << '\n';
  std::cout << "scalar multiplicities n1,n3,n5,n15 "
            << multiplicity_string(scalar_patterns) << '\n';
  std::cout << "contexts-per-support histogram "
            << histogram_string(contexts_per_support_histogram) << '\n';
  std::cout << "multiplication orbit-size histogram "
            << histogram_string(orbit_size_histogram)
            << " (telemetry; no quotient)\n";
  std::cout << "owner-local rows " << scalar_contexts * 6
            << "; feasible rows " << feasible_owner_rows << '\n';
  std::cout << "feasible-owner/context histogram "
            << histogram_string(feasible_owner_histogram) << '\n';
  std::cout << "maximum reachable sheet-union histogram "
            << histogram_string(maximum_union_histogram) << '\n';
  std::cout << "minimum owner maximum/context histogram "
            << histogram_string(context_minimum_histogram) << '\n';
  std::cout << "distinct exact owner max-union vectors "
            << owner_vectors.size() << "; vector FNV64 " << std::hex
            << owner_vector_digest << std::dec << '\n';
  std::cout << "owner-locally feasible at all six contexts "
            << all_owner_contexts << '\n';
  std::cout << "global literal unit fibres 0; global common-sheet survivors 0\n";
  std::cout << "tournament observable compare owner (feasible,max-union); tie gauge follows support-coordinate path 0->1->2->3->4->5\n";
  std::cout << "tournament fingerprints all 2184 transitive: score histogram 0,1,2,3,4,5; directed triangles 0; SCCs 6; Hamiltonian paths 1\n";
  std::cout << "tournament tie-edge histogram "
            << histogram_string(tournament_tie_histogram) << '\n';
  std::cout << "tournament edge-flip histogram against tie path "
            << histogram_string(tournament_flip_histogram) << '\n';
  std::cout << "challenged vertices exact feasible-owner subset/max-union vector is faithful for the terminal local obstruction; its owner tournament loses deficit magnitudes, provider vertices lose shared-unit choices, and sheet vertices are faithful only inside one owner DP\n";
  std::cout << "local D15 mask table at owner one (units 1,2,4,7,8,11,13,14; ratios 1..12 in hex)\n";
  for (int state = 7; state < 15; ++state) {
    std::cout << "  e=" << STATE_UNIT[state] << ':';
    for (int ratio = 1; ratio < P; ++ratio)
      std::cout << ' ' << std::hex << MASK[ratio][state][1] << std::dec;
    std::cout << '\n';
  }
}
