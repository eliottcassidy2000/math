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
// fourteen Hamming-six sheet bank.  The program reconstructs every literal
// CRT mask, enumerates the complete hereditary divisor/unit grammar, and scans
// every one of the 924 labelled supports.  Scalar capacity leaves 576 divisor
// contexts.  Exact owner-by-owner mask-union reachability then kills every
// context before there is a global unit fibre to replay.
//
// The terminal obstruction is an owner-local cardinality deficit, so the
// faithful vertices at this stage are the fourteen sheet positions for one
// fixed owner.  Runner/provider vertices and pairwise tournaments forget the
// shared choice of a provider unit.  Since no context survives all six local
// gates, there is no nonempty owner-obligation fibre on which an intersection
// graph or a tournament completion could honestly be formed.

using Labels = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;
using StateWord = std::array<uint8_t, 6>;
using Multiplicity = std::array<uint8_t, 4>;

static constexpr int P = 13;
static constexpr int C = 14;
static constexpr uint16_t FULL = (1U << C) - 1;

// All effective-order/unit states, listed literally.  D=1 has its unique
// residue; the other rows contain every and only unit modulo D.
static constexpr std::array<int, 14> STATE_ORDER{
    1, 2, 7, 7, 7, 7, 7, 7, 14, 14, 14, 14, 14, 14};
static constexpr std::array<int, 14> STATE_UNIT{
    0, 1, 1, 2, 3, 4, 5, 6, 1, 3, 5, 9, 11, 13};
static constexpr std::array<int, 4> DIVISORS{1, 2, 7, 14};
static constexpr std::array<int, 4> ORDER_REP{0, 1, 2, 8};
static constexpr std::array<int, 5> STATE_BEGIN{0, 1, 2, 8, 14};

static std::array<std::array<std::array<uint16_t, 13>, 14>, 13> MASK{};

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

// Deliberately use a literal CRT search rather than the closed formula in the
// exploratory Python scout.  This makes the two implementations independent
// at the most error-prone point of the sheet construction.
int crt_base(int label, int state) {
  const int order = STATE_ORDER[state], unit = STATE_UNIT[state];
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P &&
        value % order == unit % order)
      return value;
  fail("CRT base missing");
}

uint16_t local_mask(int label, int state, int owner) {
  const int order = STATE_ORDER[state], base = crt_base(label, state);
  const int inverse = inverse_mod_13(owner);
  uint16_t answer = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    // Literal strict-safe interval (-D,D] in the AP sheet model.
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
  // Deleting any provider must leave lcm fourteen.
  for (int omitted = 0; omitted < 6; ++omitted) {
    int value = 1;
    for (int coordinate = 0; coordinate < 6; ++coordinate)
      if (coordinate != omitted)
        value = std::lcm(value, DIVISORS[orders[coordinate]]);
    if (value != C) return false;
  }
  return true;
}

void enumerate_orders(std::vector<OrderWord> &answer, OrderWord &word,
                      int coordinate) {
  if (coordinate == 6) {
    if (hereditary(word)) answer.push_back(word);
    return;
  }
  for (int order = 0; order < 4; ++order) {
    word[coordinate] = static_cast<uint8_t>(order);
    enumerate_orders(answer, word, coordinate + 1);
  }
}

uint64_t state_fibre_size(const OrderWord &orders) {
  uint64_t answer = 1;
  for (int order : orders) {
    const auto [begin, end] = state_range(order);
    answer *= static_cast<uint64_t>(end - begin);
  }
  return answer;
}

void enumerate_states_for_order(const OrderWord &orders, StateWord &states,
                                int coordinate, uint64_t &count,
                                uint64_t &digest) {
  if (coordinate == 6) {
    ++count;
    // A deterministic FNV-1a audit word over the literal state tuples.  It is
    // telemetry, not a cryptographic proof component; the exact count and all
    // later mask computations use the enumerated states directly.
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
  Multiplicity answer{};
  for (int order : orders) ++answer[order];
  return answer;
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
  int maximum_union = 0;
  std::size_t reachable_masks = 0;
};

LocalAudit owner_local_audit(const Labels &labels, const OrderWord &orders,
                             int owner) {
  // Exact union DP over the independently selectable units for this owner.
  // No covariance or multiplication quotient is used.
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

  LocalAudit answer;
  answer.feasible =
      std::find(reachable.begin(), reachable.end(), FULL) != reachable.end();
  answer.reachable_masks = reachable.size();
  for (uint16_t mask : reachable)
    answer.maximum_union =
        std::max(answer.maximum_union, std::popcount(mask));
  return answer;
}

Labels multiply_support(const Labels &labels, int multiplier) {
  Labels answer{};
  for (int coordinate = 0; coordinate < 6; ++coordinate)
    answer[coordinate] =
        static_cast<uint8_t>(multiplier * labels[coordinate] % P);
  std::sort(answer.begin(), answer.end());
  return answer;
}

std::string support_string(const Labels &labels) {
  std::ostringstream out;
  for (int coordinate = 0; coordinate < 6; ++coordinate) {
    if (coordinate) out << ',';
    out << static_cast<int>(labels[coordinate]);
  }
  return out.str();
}

std::string orbit_string(
    const std::vector<std::pair<Labels, int>> &orbits) {
  std::ostringstream out;
  for (std::size_t orbit = 0; orbit < orbits.size(); ++orbit) {
    if (orbit) out << "; ";
    out << '[' << support_string(orbits[orbit].first) << "]:"
        << orbits[orbit].second;
  }
  return out.str();
}

int main() {
  require(STATE_ORDER.size() == STATE_UNIT.size(), "state arrays disagree");
  for (int order_index = 0; order_index < 4; ++order_index) {
    const int divisor = DIVISORS[order_index];
    require(C % divisor == 0, "listed order does not divide fourteen");
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
    for (int state = 0; state < 14; ++state)
      for (int owner = 1; owner < P; ++owner)
        MASK[label][state][owner] = local_mask(label, state, owner);

  // The scalar relaxation may use one representative state at each order
  // only after literal unit independence of cardinality has been checked.
  for (int order = 0; order < 4; ++order)
    for (int provider = 1; provider < P; ++provider)
      for (int owner = 1; owner < P; ++owner) {
        const auto [begin, end] = state_range(order);
        std::set<int> sizes;
        for (int state = begin; state < end; ++state)
          sizes.insert(std::popcount(MASK[provider][state][owner]));
        require(sizes.size() == 1,
                "mask cardinality unexpectedly depends on the unit");
      }

  std::vector<OrderWord> order_words;
  OrderWord scratch_order{};
  enumerate_orders(order_words, scratch_order, 0);
  uint64_t state_words_per_support = 0;
  for (const OrderWord &orders : order_words)
    state_words_per_support += state_fibre_size(orders);

  // Literally traverse the complete state grammar once, independently of the
  // product-of-totients count used above.
  uint64_t literal_state_words_per_support = 0;
  uint64_t state_grammar_digest = 14'695'981'039'346'656'037ULL;
  StateWord scratch_state{};
  for (const OrderWord &orders : order_words)
    enumerate_states_for_order(orders, scratch_state, 0,
                               literal_state_words_per_support,
                               state_grammar_digest);

  require(order_words.size() == 3'249,
          "hereditary divisor-word census mismatch");
  require(state_words_per_support == 6'703'884,
          "hereditary state-word census mismatch");
  require(literal_state_words_per_support == state_words_per_support,
          "literal and weighted state-word censuses disagree");
  require(state_grammar_digest == 0x3f93d84053a6bdd3ULL,
          "literal state-grammar digest mismatch");

  uint64_t supports = 0, scalar_contexts = 0, owner_feasible_rows = 0;
  uint64_t owner_local_contexts = 0;
  std::map<Multiplicity, uint64_t> scalar_patterns;
  std::map<int, uint64_t> maximum_union_histogram;
  std::set<Labels> scalar_supports, owner_local_supports;
  std::map<Labels, uint64_t> contexts_per_support;
  std::map<Labels, std::set<uint8_t>> d2_positions_per_support;
  std::map<Labels, std::set<uint8_t>> d14_positions_per_support;

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
                ++scalar_patterns[multiplicity(orders)];
                scalar_supports.insert(labels);
                ++contexts_per_support[labels];

                uint8_t d2_positions = 0, d14_positions = 0;
                for (int coordinate = 0; coordinate < 6; ++coordinate) {
                  if (orders[coordinate] == 1)
                    d2_positions |= static_cast<uint8_t>(1U << coordinate);
                  if (orders[coordinate] == 3)
                    d14_positions |= static_cast<uint8_t>(1U << coordinate);
                }
                d2_positions_per_support[labels].insert(d2_positions);
                d14_positions_per_support[labels].insert(d14_positions);

                bool every_owner_feasible = true;
                for (int owner : labels) {
                  const LocalAudit audit =
                      owner_local_audit(labels, orders, owner);
                  owner_feasible_rows += audit.feasible;
                  every_owner_feasible &= audit.feasible;
                  ++maximum_union_histogram[audit.maximum_union];
                }
                if (every_owner_feasible) {
                  ++owner_local_contexts;
                  owner_local_supports.insert(labels);
                }
              }
            }

  require(supports == 924, "support census mismatch");
  require(supports * order_words.size() == 3'002'076,
          "labelled divisor-context census mismatch");
  require(supports * state_words_per_support == 6'194'388'816ULL,
          "labelled raw state-context census mismatch");

  const std::map<Multiplicity, uint64_t> expected_patterns{
      {{{0, 2, 0, 4}}, 36}, {{{0, 2, 1, 3}}, 144},
      {{{0, 2, 2, 2}}, 216}, {{{0, 2, 3, 1}}, 144},
      {{{0, 2, 4, 0}}, 36}};
  require(scalar_contexts == 576 && scalar_supports.size() == 36,
          "scalar-capacity classification mismatch");
  require(scalar_patterns == expected_patterns,
          "scalar multiplicity-pattern classification mismatch");

  // On each surviving support one fixed provider pair has order two; all
  // sixteen independent D7/D14 choices on the other four providers occur.
  for (const Labels &labels : scalar_supports) {
    require(contexts_per_support[labels] == 16,
            "scalar support does not carry sixteen contexts");
    require(d2_positions_per_support[labels].size() == 1,
            "D2 provider pair varies on one scalar support");
    const uint8_t d2_positions = *d2_positions_per_support[labels].begin();
    require(std::popcount(static_cast<unsigned>(d2_positions)) == 2,
            "scalar context does not have exactly two D2 providers");
    require(d14_positions_per_support[labels].size() == 16,
            "D7/D14 assignments are not the complete four-cube");
    for (uint8_t d14_positions : d14_positions_per_support[labels])
      require((d14_positions & d2_positions) == 0,
              "D14 position overlaps a D2 provider");
  }

  const std::map<int, uint64_t> expected_maxima{
      {10, 96}, {11, 1'056}, {12, 2'304}};
  require(owner_feasible_rows == 0,
          "some scalar context is feasible even for one owner");
  require(maximum_union_histogram == expected_maxima,
          "owner-local maximum-union histogram mismatch");
  require(owner_local_contexts == 0 && owner_local_supports.empty(),
          "owner-local bank unexpectedly survives");

  // Multiplication by F_13^* is only structural telemetry: no orbit quotient
  // was used in either the scalar scan or the owner-local DP.
  std::set<Labels> remaining = scalar_supports;
  std::vector<std::pair<Labels, int>> orbit_data;
  while (!remaining.empty()) {
    const Labels representative = *remaining.begin();
    std::set<Labels> orbit;
    for (int multiplier = 1; multiplier < P; ++multiplier)
      orbit.insert(multiply_support(representative, multiplier));
    require(std::includes(scalar_supports.begin(), scalar_supports.end(),
                          orbit.begin(), orbit.end()),
            "multiplication orbit leaves the scalar support bank");
    int removed = 0;
    for (const Labels &labels : orbit) removed += remaining.erase(labels);
    require(removed == static_cast<int>(orbit.size()),
            "multiplication orbit decomposition overlaps");
    orbit_data.push_back({representative, removed});
  }
  require(orbit_data.size() == 3 &&
              std::all_of(orbit_data.begin(), orbit_data.end(),
                          [](const auto &entry) { return entry.second == 12; }),
          "unexpected scalar-support multiplication orbits");
  const std::vector<std::pair<Labels, int>> expected_orbits{
      {{{1, 2, 3, 5, 10, 11}}, 12},
      {{{1, 2, 3, 6, 7, 11}}, 12},
      {{{1, 2, 4, 5, 8, 11}}, 12}};
  require(orbit_data == expected_orbits,
          "scalar-support orbit representatives changed");

  std::cout << "scale-fourteen AP-centred Hamming-six owner-local deficit\n";
  std::cout << "divisor grammar 1,2,7,14; literal states 14\n";
  std::cout << "states (D,e) (1,0) (2,1) (7,1) (7,2) (7,3) (7,4) (7,5) (7,6) (14,1) (14,3) (14,5) (14,9) (14,11) (14,13)\n";
  std::cout << "supports " << supports << '\n';
  std::cout << "hereditary order words " << order_words.size() << '\n';
  std::cout << "hereditary labelled order contexts "
            << supports * order_words.size() << '\n';
  std::cout << "hereditary state words per support "
            << state_words_per_support << '\n';
  std::cout << "literal state-grammar FNV64 " << std::hex
            << state_grammar_digest << std::dec << '\n';
  std::cout << "hereditary labelled state contexts "
            << supports * state_words_per_support << '\n';
  std::cout << "scalar-capacity contexts " << scalar_contexts
            << " on " << scalar_supports.size() << " supports\n";
  std::cout << "scalar pattern two D2 providers plus four D7/D14 providers; sixteen contexts per support\n";
  std::cout << "multiplication orbits " << orbit_string(orbit_data)
            << " (telemetry; no quotient)\n";
  std::cout << "owner-local rows " << scalar_contexts * 6
            << "; feasible rows " << owner_feasible_rows << '\n';
  std::cout << "maximum reachable sheet-union histogram 10:96 11:1056 12:2304\n";
  std::cout << "owner-local contexts " << owner_local_contexts << '\n';
  std::cout << "global literal unit fibres 0; global common-sheet survivors 0\n";
  std::cout << "tournament N/A: the owner-local gate is empty, so there is no owner-obligation global fibre or pair-intersection observable to orient\n";
  std::cout << "challenged vertices fourteen sheet positions preserve the exact owner-local cover deficit; runner/provider or divisor vertices destroy shared-unit incidence, while owner-obligation vertices arise only after a global fibre survives\n";
  std::cout << "local D14 mask table at owner one (units 1,3,5,9,11,13; ratios 1..12 in hex)\n";
  for (int state = 8; state < 14; ++state) {
    std::cout << "  e=" << STATE_UNIT[state] << ':';
    for (int ratio = 1; ratio < P; ++ratio)
      std::cout << ' ' << std::hex << MASK[ratio][state][1] << std::dec;
    std::cout << '\n';
  }
}
