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
#include <tuple>
#include <utility>
#include <vector>

// Exact primary certificate for the primitive proper AP-centred common-scale-
// twelve Hamming-six sheet bank.  The program reconstructs the complete
// divisor/unit grammar, scans every hereditary divisor word on every support,
// applies only proved scalar-capacity and owner-local necessary conditions,
// and literally replays every unit word in every remaining context.
// Multiplication orbits and tournaments are telemetry only: neither is used
// to quotient the exact census or the terminal replay.

using Labels = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;
using StateWord = std::array<uint8_t, 6>;
using Profile = std::array<uint64_t, 64>;

static constexpr int P = 13;
static constexpr int C = 12;
static constexpr uint16_t FULL = (1U << C) - 1;

// The twelve literal (effective order, unit) states.  State zero is the
// unique residue modulo one.  All other units are listed rather than inferred
// from a cardinality formula, so the CRT generator below audits every state.
static constexpr std::array<int, 12> STATE_ORDER{
    1, 2, 3, 3, 4, 4, 6, 6, 12, 12, 12, 12};
static constexpr std::array<int, 12> STATE_UNIT{
    0, 1, 1, 2, 1, 3, 1, 5, 1, 5, 7, 11};
static constexpr std::array<int, 6> DIVISORS{1, 2, 3, 4, 6, 12};
static constexpr std::array<int, 6> ORDER_REP{0, 1, 2, 4, 6, 8};
static constexpr std::array<int, 7> STATE_BEGIN{0, 1, 2, 4, 6, 8, 12};

static std::array<std::array<std::array<uint16_t, 13>, 12>, 13> MASK{};

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
    // This is the literal strict-safe interval (-D,D] from the sheet model.
    if (-order < value && value <= order)
      answer |= static_cast<uint16_t>(1U << sheet);
  }
  return answer;
}

int divisor_for_index(int order_index) { return DIVISORS[order_index]; }

std::pair<int, int> state_range(int order_index) {
  require(0 <= order_index && order_index < 6, "bad order index");
  return {STATE_BEGIN[order_index], STATE_BEGIN[order_index + 1]};
}

bool hereditary(const OrderWord &orders) {
  // Hereditary leave-one-out lcm: deletion of any provider still sees the
  // full common scale twelve.
  for (int omitted = 0; omitted < 6; ++omitted) {
    int value = 1;
    for (int coordinate = 0; coordinate < 6; ++coordinate)
      if (coordinate != omitted)
        value = std::lcm(value, divisor_for_index(orders[coordinate]));
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
  for (int order = 0; order < 6; ++order) {
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

using Multiplicity = std::array<uint8_t, 6>;

Multiplicity multiplicity(const OrderWord &orders) {
  Multiplicity answer{};
  for (int order : orders) ++answer[order];
  return answer;
}

bool scalar_capacity(const Labels &labels, const OrderWord &orders, int owner) {
  int capacity = 0;
  for (int provider = 0; provider < 6; ++provider) {
    capacity += std::popcount(
        MASK[labels[provider]][ORDER_REP[orders[provider]]][owner]);
  }
  return capacity >= C;
}

bool owner_locally_feasible(const Labels &labels, const OrderWord &orders,
                            int owner) {
  // Exact six-fold mask-union dynamic program for one owner.  This is a
  // relaxation only in that every owner chooses its units independently.
  std::vector<uint16_t> reachable{0};
  std::array<uint8_t, 4096> seen{};
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
  return std::find(reachable.begin(), reachable.end(), FULL) != reachable.end();
}

int projective_class(int label) { return std::min(label, P - label); }

bool sign_transversal(const Labels &labels) {
  std::set<int> classes;
  for (int label : labels) classes.insert(projective_class(label));
  return classes == std::set<int>{1, 2, 3, 4, 5, 6};
}

Labels multiply_support(const Labels &labels, int multiplier) {
  Labels answer{};
  for (int i = 0; i < 6; ++i)
    answer[i] = static_cast<uint8_t>(multiplier * labels[i] % P);
  std::sort(answer.begin(), answer.end());
  return answer;
}

uint8_t owner_satisfaction_mask(const Labels &labels,
                                const StateWord &states) {
  uint8_t answer = 0;
  for (int owner_index = 0; owner_index < 6; ++owner_index) {
    uint16_t cover = 0;
    for (int provider = 0; provider < 6; ++provider)
      cover |= MASK[labels[provider]][states[provider]][labels[owner_index]];
    if (cover == FULL)
      answer |= static_cast<uint8_t>(1U << owner_index);
  }
  return answer;
}

Profile replay_all_d12_fibre(const Labels &labels) {
  Profile profile{};
  // Four D12 units at every one of the six labelled coordinates: all 4^6
  // words are visited literally, with no multiplication or sign quotient.
  for (int code = 0; code < 4096; ++code) {
    int quotient = code;
    StateWord states{};
    for (int coordinate = 0; coordinate < 6; ++coordinate) {
      states[coordinate] = static_cast<uint8_t>(8 + quotient % 4);
      quotient /= 4;
    }
    ++profile[owner_satisfaction_mask(labels, states)];
  }
  return profile;
}

std::array<uint64_t, 7> satisfaction_histogram(const Profile &profile) {
  std::array<uint64_t, 7> answer{};
  for (int subset = 0; subset < 64; ++subset)
    answer[std::popcount(static_cast<unsigned>(subset))] += profile[subset];
  return answer;
}

std::array<uint64_t, 6> owner_sizes(const Profile &profile) {
  std::array<uint64_t, 6> answer{};
  for (int owner = 0; owner < 6; ++owner)
    for (int subset = 0; subset < 64; ++subset)
      if (subset & (1 << owner)) answer[owner] += profile[subset];
  return answer;
}

std::array<std::array<uint64_t, 6>, 6>
pair_intersections(const Profile &profile) {
  std::array<std::array<uint64_t, 6>, 6> answer{};
  for (int first = 0; first < 6; ++first)
    for (int second = 0; second < 6; ++second)
      for (int subset = 0; subset < 64; ++subset)
        if ((subset & (1 << first)) && (subset & (1 << second)))
          answer[first][second] += profile[subset];
  return answer;
}

struct TournamentFingerprint {
  std::array<int, 6> scores{};
  int directed_triangles = 0;
  std::array<int, 6> scc_sizes{};
  int gauge_edge_flips = 0;
  int hamiltonian_paths = 0;
};

bool tournament_edge(int source, int target) {
  // The faithful relation |O_source intersect O_target|>0 is empty.  After
  // the projective sign gauge, every pair is therefore a tie, broken by the
  // Hamiltonian path on projective classes 1<2<...<6.
  return projective_class(source) < projective_class(target);
}

bool tournament_reachable(int source, int target, const Labels &labels) {
  std::vector<int> stack{source};
  std::set<int> seen{source};
  while (!stack.empty()) {
    const int vertex = stack.back();
    stack.pop_back();
    for (int neighbor : labels) {
      if (neighbor == vertex || !tournament_edge(vertex, neighbor) ||
          seen.count(neighbor))
        continue;
      seen.insert(neighbor);
      stack.push_back(neighbor);
    }
  }
  return seen.count(target);
}

TournamentFingerprint tournament_fingerprint(const Labels &labels) {
  require(sign_transversal(labels), "tournament support is not a sign transversal");
  TournamentFingerprint answer;

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j)
      if (i != j) answer.scores[i] += tournament_edge(labels[i], labels[j]);
    for (int j = i + 1; j < 6; ++j)
      answer.gauge_edge_flips += tournament_edge(labels[j], labels[i]);
  }
  std::sort(answer.scores.begin(), answer.scores.end());

  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j)
      for (int k = j + 1; k < 6; ++k) {
        const std::array<int, 3> triple{labels[i], labels[j], labels[k]};
        bool cycle = true;
        for (int vertex : triple) {
          int outdegree = 0;
          for (int other : triple)
            if (other != vertex) outdegree += tournament_edge(vertex, other);
          cycle &= outdegree == 1;
        }
        answer.directed_triangles += cycle;
      }

  std::set<int> unused(labels.begin(), labels.end());
  int component = 0;
  while (!unused.empty()) {
    const int seed = *unused.begin();
    std::vector<int> members;
    for (int vertex : labels)
      if (tournament_reachable(seed, vertex, labels) &&
          tournament_reachable(vertex, seed, labels))
        members.push_back(vertex);
    require(!members.empty(), "empty tournament SCC");
    answer.scc_sizes[component++] = static_cast<int>(members.size());
    for (int vertex : members) unused.erase(vertex);
  }
  std::sort(answer.scc_sizes.begin(), answer.scc_sizes.end(), std::greater<>());

  Labels path = labels;
  do {
    bool directed = true;
    for (int i = 0; i < 5; ++i)
      directed &= tournament_edge(path[i], path[i + 1]);
    answer.hamiltonian_paths += directed;
  } while (std::next_permutation(path.begin(), path.end()));
  return answer;
}

std::string orbit_sizes_string(
    const std::vector<std::pair<Labels, int>> &orbits) {
  std::ostringstream out;
  for (std::size_t i = 0; i < orbits.size(); ++i) {
    if (i) out << ',';
    out << orbits[i].second;
  }
  return out.str();
}

int main() {
  require(STATE_ORDER.size() == STATE_UNIT.size(), "state arrays disagree");
  for (int order_index = 0; order_index < 6; ++order_index) {
    const int divisor = DIVISORS[order_index];
    require(C % divisor == 0, "listed order does not divide twelve");
    const auto [begin, end] = state_range(order_index);
    require(begin == ORDER_REP[order_index], "wrong order representative");
    std::set<int> actual_units;
    for (int state = begin; state < end; ++state) {
      require(STATE_ORDER[state] == divisor, "state assigned to wrong order");
      require(std::gcd(STATE_UNIT[state], divisor) == 1 || divisor == 1,
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
    for (int state = 0; state < 12; ++state)
      for (int owner = 1; owner < P; ++owner)
        MASK[label][state][owner] = local_mask(label, state, owner);

  // Unit-independent cardinality is the scalar-capacity relaxation.  Check it
  // from the literal masks at every divisor and every provider/owner pair;
  // the census below uses these direct masks and trusts no covariance quotient.
  for (int order = 0; order < 6; ++order)
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

  require(order_words.size() == 26'961,
          "hereditary divisor-word census mismatch");
  require(state_words_per_support == 2'611'968,
          "hereditary state-word census mismatch");

  uint64_t supports = 0, scalar_contexts = 0, local_contexts = 0;
  std::set<Multiplicity> scalar_patterns;
  std::set<Labels> scalar_supports, local_supports;
  std::vector<std::pair<Labels, OrderWord>> local_bank;

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
                  scalar &= scalar_capacity(labels, orders, owner);
                if (!scalar) continue;
                ++scalar_contexts;
                scalar_patterns.insert(multiplicity(orders));
                scalar_supports.insert(labels);

                bool local = true;
                for (int owner : labels) {
                  if (!owner_locally_feasible(labels, orders, owner)) {
                    local = false;
                    break;
                  }
                }
                if (!local) continue;
                ++local_contexts;
                local_supports.insert(labels);
                local_bank.push_back({labels, orders});
              }
            }

  require(supports == 924, "support census mismatch");
  require(supports * order_words.size() == 24'911'964,
          "labelled divisor-context census mismatch");
  require(supports * state_words_per_support == 2'413'458'432ULL,
          "labelled raw state-context census mismatch");
  require(scalar_contexts == 36'830 && scalar_patterns.size() == 85,
          "scalar-capacity classification mismatch");
  require(scalar_supports.size() == 912,
          "scalar-capacity support projection mismatch");
  require(local_contexts == 64 && local_supports.size() == 64 &&
              local_bank.size() == 64,
          "owner-local bank size mismatch");

  std::set<Labels> expected_transversals;
  for (int signs = 0; signs < 64; ++signs) {
    Labels labels{};
    for (int projective = 1; projective <= 6; ++projective)
      labels[projective - 1] = static_cast<uint8_t>(
          signs & (1 << (projective - 1)) ? P - projective : projective);
    std::sort(labels.begin(), labels.end());
    expected_transversals.insert(labels);
  }
  require(local_supports == expected_transversals,
          "owner-local supports are not exactly the sign transversals");
  for (const auto &[labels, orders] : local_bank) {
    require(sign_transversal(labels), "local support fails sign-transversal test");
    require(std::all_of(orders.begin(), orders.end(),
                        [](uint8_t order) { return order == 5; }),
            "owner-local context is not all-D12");
  }

  // Record multiplicative structure, but replay every support below.  The
  // expected decomposition is asserted only as an independent checksum.
  std::set<Labels> remaining = local_supports;
  std::vector<std::pair<Labels, int>> orbit_data;
  while (!remaining.empty()) {
    const Labels representative = *remaining.begin();
    std::set<Labels> orbit;
    for (int multiplier = 1; multiplier < P; ++multiplier)
      orbit.insert(multiply_support(representative, multiplier));
    require(std::includes(local_supports.begin(), local_supports.end(),
                          orbit.begin(), orbit.end()),
            "multiplication orbit leaves the local support bank");
    int removed = 0;
    for (const Labels &labels : orbit) removed += remaining.erase(labels);
    require(removed == static_cast<int>(orbit.size()),
            "multiplication orbit decomposition overlaps");
    orbit_data.push_back({representative, removed});
  }
  const std::vector<std::pair<Labels, int>> expected_orbits{
      {{{1, 2, 3, 4, 5, 6}}, 12}, {{{1, 2, 3, 4, 5, 7}}, 12},
      {{{1, 2, 3, 4, 6, 8}}, 12}, {{{1, 2, 3, 5, 6, 9}}, 4},
      {{{1, 2, 3, 5, 7, 9}}, 12}, {{{1, 2, 6, 8, 9, 10}}, 12}};
  require(orbit_data == expected_orbits,
          "unexpected multiplication-orbit decomposition");

  const std::array<uint64_t, 7> expected_histogram{3'808, 288, 0, 0, 0, 0, 0};
  const std::array<uint64_t, 6> expected_owner_sizes{48, 48, 48, 48, 48, 48};
  std::map<int, int> edge_flip_histogram;
  uint64_t replayed_words = 0, global_survivors = 0;
  for (const Labels &labels : local_supports) {
    const Profile profile = replay_all_d12_fibre(labels);
    replayed_words += 4096;
    global_survivors += profile[63];
    require(satisfaction_histogram(profile) == expected_histogram,
            "owner-satisfaction histogram mismatch");
    require(owner_sizes(profile) == expected_owner_sizes,
            "owner-obligation size mismatch");
    require(profile[0] == 3'808, "empty owner subset has wrong multiplicity");
    for (int owner = 0; owner < 6; ++owner)
      require(profile[1 << owner] == 48,
              "singleton owner subset has wrong multiplicity");
    for (int subset = 1; subset < 64; ++subset)
      if (std::popcount(static_cast<unsigned>(subset)) >= 2)
        require(profile[subset] == 0,
                "two or more owners are simultaneously satisfiable");
    const auto pairs = pair_intersections(profile);
    for (int first = 0; first < 6; ++first)
      for (int second = first + 1; second < 6; ++second)
        require(pairs[first][second] == 0,
                "distinct owner obligations intersect");

    const TournamentFingerprint fingerprint =
        tournament_fingerprint(labels);
    require(fingerprint.scores == std::array<int, 6>{0, 1, 2, 3, 4, 5},
            "tournament score sequence mismatch");
    require(fingerprint.directed_triangles == 0,
            "tournament unexpectedly has a directed triangle");
    require(fingerprint.scc_sizes == std::array<int, 6>{1, 1, 1, 1, 1, 1},
            "tournament SCC fingerprint mismatch");
    require(fingerprint.hamiltonian_paths == 1,
            "tournament Hamiltonian-path count mismatch");
    ++edge_flip_histogram[fingerprint.gauge_edge_flips];
  }
  const std::map<int, int> expected_edge_flip_histogram{
      {0, 2}, {1, 2}, {2, 2}, {3, 4}, {4, 4}, {5, 6},
      {6, 6}, {7, 6}, {8, 6}, {9, 6}, {10, 6}, {11, 4},
      {12, 4}, {13, 2}, {14, 2}, {15, 2}};
  require(edge_flip_histogram == expected_edge_flip_histogram,
          "projective-gauge edge-flip histogram mismatch");
  require(replayed_words == 262'144 && global_survivors == 0,
          "global literal replay mismatch");

  std::cout << "scale-twelve AP-centred Hamming-six owner orthogonality\n";
  std::cout << "divisor grammar 1,2,3,4,6,12; literal states 12\n";
  std::cout << "states (D,e) (1,0) (2,1) (3,1) (3,2) (4,1) (4,3) (6,1) (6,5) (12,1) (12,5) (12,7) (12,11)\n";
  std::cout << "supports " << supports << '\n';
  std::cout << "hereditary order words " << order_words.size() << '\n';
  std::cout << "hereditary labelled order contexts "
            << supports * order_words.size() << '\n';
  std::cout << "hereditary state words per support "
            << state_words_per_support << '\n';
  std::cout << "hereditary labelled state contexts "
            << supports * state_words_per_support << '\n';
  std::cout << "scalar-capacity contexts " << scalar_contexts
            << " across multiplicity patterns " << scalar_patterns.size() << '\n';
  std::cout << "owner-local contexts " << local_contexts
            << "; all D12 on the 64 projective sign transversals\n";
  std::cout << "multiplication orbits " << orbit_data.size() << " sizes "
            << orbit_sizes_string(orbit_data) << " (telemetry; no quotient)\n";
  std::cout << "replayed literal unit words " << replayed_words << '\n';
  std::cout << "uniform profile 0:3808 1:288; each owner size 48; distinct pairs disjoint\n";
  std::cout << "global common-sheet survivors " << global_survivors << '\n';
  std::cout << "tournament pair observable |O_i intersect O_j|; positive-intersection relation empty, all 15 pairs tied\n";
  std::cout << "tournament projective sign gauge owner +/-r -> r; tie Hamiltonian path 1<2<3<4<5<6\n";
  std::cout << "tournament fingerprints scores 0,1,2,3,4,5; directed triangles 0; SCCs 1,1,1,1,1,1; Hamiltonian paths 1\n";
  std::cout << "tournament gauge edge-flip histogram 0:2 1:2 2:2 3:4 4:4 5:6 6:6 7:6 8:6 9:6 10:6 11:4 12:4 13:2 14:2 15:2\n";
  std::cout << "challenged vertices owner obligations preserve common-sheet intersections; sign classes preserve the local-bank support classification but destroy unit-mask incidence; runner/provider vertices destroy the owner nerve\n";
  std::cout << "local D12 mask table at owner one (units 1,5,7,11; ratios 1..12 in hex)\n";
  for (int state = 8; state < 12; ++state) {
    std::cout << "  e=" << STATE_UNIT[state] << ':';
    for (int ratio = 1; ratio < P; ++ratio)
      std::cout << ' ' << std::hex << MASK[ratio][state][1] << std::dec;
    std::cout << '\n';
  }
}
