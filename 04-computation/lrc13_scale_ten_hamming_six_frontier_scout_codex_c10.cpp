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
#include <tuple>
#include <vector>

// Exact frontier scout for the primitive proper AP-centred common-scale-ten
// Hamming-six sheet bank.  This deliberately separates the hereditary-lcm
// divisor grammar, unit-independent scalar capacity, owner-local feasibility,
// and the globally shared unit word.

using Labels = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;
using StateWord = std::array<uint8_t, 6>;

static constexpr int P = 13;
static constexpr int C = 10;
static constexpr std::array<int, 10> ORDERS{1, 2, 5, 5, 5, 5,
                                             10, 10, 10, 10};
static constexpr std::array<int, 10> UNITS{0, 1, 1, 2, 3, 4,
                                            1, 3, 7, 9};
static constexpr std::array<int, 4> DIVISORS{1, 2, 5, 10};
static constexpr std::array<int, 4> ORDER_REP{0, 1, 2, 6};
static std::array<std::array<std::array<uint16_t, 13>, 10>, 13> MASK{};

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
  const int order = ORDERS[state], unit = UNITS[state];
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P && value % order == unit % order)
      return value;
  fail("CRT base missing");
}

uint16_t local_mask(int label, int state, int owner) {
  const int order = ORDERS[state], base = crt_base(label, state);
  const int inverse = inverse_mod_13(owner);
  uint16_t answer = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    if (-order < value && value <= order)
      answer |= static_cast<uint16_t>(1U << sheet);
  }
  return answer;
}

const std::vector<int> &states_for_order(int order_index) {
  static const std::vector<int> d1{0};
  static const std::vector<int> d2{1};
  static const std::vector<int> d5{2, 3, 4, 5};
  static const std::vector<int> d10{6, 7, 8, 9};
  if (order_index == 0) return d1;
  if (order_index == 1) return d2;
  if (order_index == 2) return d5;
  if (order_index == 3) return d10;
  fail("bad order index");
}

int divisor_for_index(int order_index) { return DIVISORS[order_index]; }

bool hereditary(const OrderWord &orders) {
  for (int omitted = 0; omitted < 6; ++omitted) {
    int value = 1;
    for (int i = 0; i < 6; ++i)
      if (i != omitted)
        value = std::lcm(value, divisor_for_index(orders[i]));
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
  for (int order : orders) answer *= states_for_order(order).size();
  return answer;
}

std::string pattern(const OrderWord &orders) {
  std::array<int, 4> counts{};
  for (int order : orders) ++counts[order];
  std::ostringstream out;
  out << "1^" << counts[0] << " 2^" << counts[1]
      << " 5^" << counts[2] << " 10^" << counts[3];
  return out.str();
}

bool scalar_capacity(const Labels &labels, const OrderWord &orders, int owner) {
  int capacity = 0;
  for (int provider = 0; provider < 6; ++provider)
    capacity += std::popcount(
        MASK[labels[provider]][ORDER_REP[orders[provider]]][owner]);
  return capacity >= C;
}

bool owner_locally_feasible(const Labels &labels, const OrderWord &orders,
                            int owner) {
  std::array<bool, 1024> reachable{};
  reachable[0] = true;
  for (int provider = 0; provider < 6; ++provider) {
    std::array<bool, 1024> next{};
    for (int partial = 0; partial < 1024; ++partial) {
      if (!reachable[partial]) continue;
      for (int state : states_for_order(orders[provider]))
        next[partial | MASK[labels[provider]][state][owner]] = true;
    }
    reachable = next;
  }
  return reachable[1023];
}

uint8_t owner_satisfaction_mask(const Labels &labels, const StateWord &states) {
  uint8_t answer = 0;
  for (int owner = 0; owner < 6; ++owner) {
    uint16_t cover = 0;
    for (int provider = 0; provider < 6; ++provider)
      cover |= MASK[labels[provider]][states[provider]][labels[owner]];
    if (cover == 1023) answer |= static_cast<uint8_t>(1U << owner);
  }
  return answer;
}

void enumerate_global_fibre(const Labels &labels, const OrderWord &orders,
                            StateWord &states, int coordinate,
                            std::array<uint64_t, 64> &profile) {
  if (coordinate == 6) {
    ++profile[owner_satisfaction_mask(labels, states)];
    return;
  }
  for (int state : states_for_order(orders[coordinate])) {
    states[coordinate] = static_cast<uint8_t>(state);
    enumerate_global_fibre(labels, orders, states, coordinate + 1, profile);
  }
}

std::array<uint64_t, 7> satisfaction_histogram(
    const std::array<uint64_t, 64> &profile) {
  std::array<uint64_t, 7> answer{};
  for (int subset = 0; subset < 64; ++subset)
    answer[std::popcount(static_cast<unsigned>(subset))] += profile[subset];
  return answer;
}

std::array<uint64_t, 6> owner_intersections(
    const std::array<uint64_t, 64> &profile) {
  std::array<uint64_t, 6> answer{};
  for (int i = 0; i < 6; ++i)
    for (int subset = 0; subset < 64; ++subset)
      if (subset & (1 << i)) answer[i] += profile[subset];
  return answer;
}

std::array<std::array<uint64_t, 6>, 6> pair_intersections(
    const std::array<uint64_t, 64> &profile) {
  std::array<std::array<uint64_t, 6>, 6> answer{};
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      for (int subset = 0; subset < 64; ++subset)
        if ((subset & (1 << i)) && (subset & (1 << j)))
          answer[i][j] += profile[subset];
  return answer;
}

std::string support_string(const Labels &labels) {
  std::ostringstream out;
  for (int i = 0; i < 6; ++i) {
    if (i) out << ',';
    out << static_cast<int>(labels[i]);
  }
  return out.str();
}

int projective_class(int label) { return std::min(label, P - label); }

bool sign_transversal(const Labels &labels) {
  std::set<int> classes;
  for (int label : labels) classes.insert(projective_class(label));
  return classes == std::set<int>{1, 2, 3, 4, 5, 6};
}

bool zero_obligation_pair(int first, int second) {
  const int quotient = second * inverse_mod_13(first) % P;
  return quotient == 2 || quotient == 6 || quotient == 7 || quotient == 11;
}

std::set<std::pair<int, int>> projective_cycle_edges(const Labels &labels) {
  // Switching every selected residue to its positive representative gauges the
  // zero-pair C6 as 1 -> 2 -> 4 -> 5 -> 3 -> 6 -> 1.
  static constexpr std::array<int, 6> cycle{1, 2, 4, 5, 3, 6};
  std::array<int, 7> successor{};
  for (int i = 0; i < 6; ++i) successor[cycle[i]] = cycle[(i + 1) % 6];
  std::set<std::pair<int, int>> edges;
  for (int source : labels)
    for (int target : labels)
      if (source != target &&
          successor[projective_class(source)] == projective_class(target))
        edges.insert({source, target});
  return edges;
}

bool reachable(int source, int target, const Labels &labels,
               const std::set<std::pair<int, int>> &edges) {
  std::vector<int> stack{source};
  std::set<int> seen{source};
  while (!stack.empty()) {
    const int vertex = stack.back();
    stack.pop_back();
    for (int neighbor : labels) {
      if (!edges.count({vertex, neighbor}) || seen.count(neighbor)) continue;
      seen.insert(neighbor);
      stack.push_back(neighbor);
    }
  }
  return seen.count(target);
}

using TournamentFingerprint =
    std::tuple<std::array<int, 6>, int, std::array<int, 6>, int, int>;

TournamentFingerprint tournament_fingerprint(const Labels &labels) {
  const auto sparse = projective_cycle_edges(labels);
  require(sparse.size() == 6, "projective sparse relation is not C6");
  Labels tie_path = labels;
  bool found = false;
  do {
    bool path = true;
    for (int i = 0; i < 5; ++i)
      path &= !sparse.count({tie_path[i], tie_path[i + 1]}) &&
              !sparse.count({tie_path[i + 1], tie_path[i]});
    if (path) {
      found = true;
      break;
    }
  } while (std::next_permutation(tie_path.begin(), tie_path.end()));
  require(found, "tie graph has no Hamiltonian path");

  std::array<int, 13> position{};
  for (int i = 0; i < 6; ++i) position[tie_path[i]] = i;
  std::set<std::pair<int, int>> tournament;
  int sparse_flips = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j) {
      const int first = labels[i], second = labels[j];
      std::pair<int, int> edge;
      if (sparse.count({first, second})) edge = {first, second};
      else if (sparse.count({second, first})) edge = {second, first};
      else edge = position[first] < position[second]
                      ? std::pair{first, second}
                      : std::pair{second, first};
      tournament.insert(edge);
      sparse_flips += sparse.count(edge) &&
                      position[edge.first] > position[edge.second];
    }

  std::array<int, 6> scores{};
  for (int i = 0; i < 6; ++i)
    for (int other : labels)
      scores[i] += tournament.count({labels[i], other});
  std::sort(scores.begin(), scores.end());

  int triangles = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j)
      for (int k = j + 1; k < 6; ++k) {
        const std::array<int, 3> triple{labels[i], labels[j], labels[k]};
        bool directed_cycle = true;
        for (int vertex : triple) {
          int outdegree = 0;
          for (int other : triple)
            outdegree += tournament.count({vertex, other});
          directed_cycle &= outdegree == 1;
        }
        triangles += directed_cycle;
      }

  std::array<int, 6> scc_sizes{};
  std::set<int> unused(labels.begin(), labels.end());
  int component_index = 0;
  while (!unused.empty()) {
    const int seed = *unused.begin();
    std::set<int> component;
    for (int vertex : labels)
      if (reachable(seed, vertex, labels, tournament) &&
          reachable(vertex, seed, labels, tournament))
        component.insert(vertex);
    require(!component.empty(), "empty tournament SCC");
    for (int vertex : component) unused.erase(vertex);
    scc_sizes[component_index++] = static_cast<int>(component.size());
  }
  std::sort(scc_sizes.begin(), scc_sizes.end(), std::greater<>());

  int hamiltonian_paths = 0;
  Labels path = labels;
  do {
    bool directed = true;
    for (int i = 0; i < 5; ++i)
      directed &= tournament.count({path[i], path[i + 1]});
    hamiltonian_paths += directed;
  } while (std::next_permutation(path.begin(), path.end()));
  return {scores, triangles, scc_sizes, sparse_flips, hamiltonian_paths};
}

int main() {
  for (int label = 1; label < P; ++label)
    for (int state = 0; state < 10; ++state)
      for (int owner = 1; owner < P; ++owner)
        MASK[label][state][owner] = local_mask(label, state, owner);

  for (int order = 0; order < 4; ++order)
    for (int ratio = 1; ratio < P; ++ratio) {
      std::set<int> sizes;
      for (int state : states_for_order(order))
        sizes.insert(std::popcount(MASK[ratio][state][1]));
      require(sizes.size() == 1,
              "mask cardinality unexpectedly depends on the unit");
    }

  std::vector<OrderWord> order_words;
  OrderWord scratch_order{};
  enumerate_orders(order_words, scratch_order, 0);
  uint64_t hereditary_state_words = 0;
  for (const auto &orders : order_words)
    hereditary_state_words += state_fibre_size(orders);

  uint64_t supports = 0, scalar_contexts = 0, local_contexts = 0;
  uint64_t local_state_words = 0, global_survivors = 0;
  std::set<std::string> scalar_supports, local_supports, surviving_supports;
  std::map<std::string, uint64_t> scalar_patterns, local_patterns,
      local_pattern_words;
  std::vector<std::pair<Labels, OrderWord>> local_bank;
  std::set<Labels> local_label_bank, sign_transversal_bank;

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
              if (sign_transversal(labels)) sign_transversal_bank.insert(labels);
              for (const auto &orders : order_words) {
                bool scalar = true;
                for (int owner : labels)
                  scalar &= scalar_capacity(labels, orders, owner);
                if (!scalar) continue;
                ++scalar_contexts;
                ++scalar_patterns[pattern(orders)];
                scalar_supports.insert(support_string(labels));
                bool local = true;
                for (int owner : labels)
                  local &= owner_locally_feasible(labels, orders, owner);
                if (!local) continue;
                ++local_contexts;
                ++local_patterns[pattern(orders)];
                local_pattern_words[pattern(orders)] += state_fibre_size(orders);
                local_state_words += state_fibre_size(orders);
                local_supports.insert(support_string(labels));
                local_bank.push_back({labels, orders});
                local_label_bank.insert(labels);
              }
            }

  std::map<TournamentFingerprint, int> tournament_fingerprints;
  for (const auto &[labels, orders] : local_bank) {
    StateWord states{};
    std::array<uint64_t, 64> profile{};
    enumerate_global_fibre(labels, orders, states, 0, profile);
    global_survivors += profile[63];
    if (profile[63]) surviving_supports.insert(support_string(labels));
    const auto histogram = satisfaction_histogram(profile);
    const auto owners = owner_intersections(profile);
    const auto pairs = pair_intersections(profile);
    require(std::all_of(orders.begin(), orders.end(),
                        [](int order) { return order == 3; }),
            "owner-local bank contains a non-D10 word");
    require(sign_transversal(labels),
            "owner-local support is not a sign transversal");
    require(histogram ==
                std::array<uint64_t, 7>{3'940, 120, 36, 0, 0, 0, 0},
            "owner-obligation profile mismatch");
    require(owners == std::array<uint64_t, 6>{32, 32, 32, 32, 32, 32},
            "owner-obligation cardinality mismatch");
    for (int i = 0; i < 6; ++i)
      for (int j = i + 1; j < 6; ++j)
        require(pairs[i][j] ==
                    (zero_obligation_pair(labels[i], labels[j]) ? 0 : 4),
                "owner-obligation pair nerve mismatch");
    ++tournament_fingerprints[tournament_fingerprint(labels)];
  }

  const std::map<std::string, uint64_t> expected_scalar_patterns{
      {"1^0 2^0 5^0 10^6", 64}, {"1^0 2^0 5^3 10^3", 120},
      {"1^0 2^2 5^0 10^4", 36}, {"1^0 2^2 5^2 10^2", 48},
      {"1^0 2^2 5^3 10^1", 48}, {"1^0 2^2 5^4 10^0", 12},
      {"1^0 2^3 5^0 10^3", 344}, {"1^0 2^3 5^1 10^2", 336},
      {"1^0 2^3 5^2 10^1", 144}, {"1^0 2^3 5^3 10^0", 48}};
  require(supports == 924 && order_words.size() == 3'249 &&
              hereditary_state_words == 889'200,
          "raw hereditary census mismatch");
  require(scalar_contexts == 1'200 && scalar_supports.size() == 388 &&
              scalar_patterns == expected_scalar_patterns,
          "scalar-capacity census mismatch");
  require(local_contexts == 64 && local_supports.size() == 64 &&
              local_label_bank == sign_transversal_bank &&
              sign_transversal_bank.size() == 64 &&
              local_patterns ==
                  std::map<std::string, uint64_t>{{"1^0 2^0 5^0 10^6", 64}} &&
              local_state_words == 262'144,
          "owner-local sign-transversal classification mismatch");
  require(global_survivors == 0 && surviving_supports.empty(),
          "global common-sheet bank is nonempty");

  std::map<int, int> flip_histogram, path_histogram;
  for (const auto &[fingerprint, multiplicity] : tournament_fingerprints) {
    const auto &[scores, triangles, scc, flips, paths] = fingerprint;
    require(scores == std::array<int, 6>{1, 2, 2, 3, 3, 4} &&
                triangles == 6 &&
                scc == std::array<int, 6>{6, 0, 0, 0, 0, 0},
            "unexpected completed-tournament fingerprint");
    flip_histogram[flips] += multiplicity;
    path_histogram[paths] += multiplicity;
  }
  require(tournament_fingerprints.size() == 5 &&
              flip_histogram == std::map<int, int>{{2, 8}, {3, 52}, {4, 4}} &&
              path_histogram ==
                  std::map<int, int>{{29, 32}, {31, 20}, {37, 12}},
          "completed-tournament histogram mismatch");

  std::cout << "scale-ten AP-centred Hamming-six frontier scout\n";
  std::cout << "supports " << supports << '\n';
  std::cout << "hereditary order words " << order_words.size() << '\n';
  std::cout << "hereditary order contexts " << supports * order_words.size() << '\n';
  std::cout << "hereditary state words per support " << hereditary_state_words << '\n';
  std::cout << "hereditary labelled state contexts "
            << supports * hereditary_state_words << '\n';
  std::cout << "scalar-capacity contexts " << scalar_contexts << '\n';
  std::cout << "scalar-capacity supports " << scalar_supports.size() << '\n';
  std::cout << "owner-local contexts " << local_contexts << '\n';
  std::cout << "owner-local supports " << local_supports.size() << '\n';
  std::cout << "owner-local state words " << local_state_words << '\n';
  std::cout << "global common-sheet survivors " << global_survivors << '\n';
  std::cout << "global surviving supports " << surviving_supports.size() << '\n';
  std::cout << "owner-local structure 64 sign transversals of F13*/{+-1}, all D10\n";
  std::cout << "obligation profile 0:3940 1:120 2:36; every owner size 32\n";
  std::cout << "obligation nerve K6\\C6 (triangular-prism graph), edge intersections 4, no higher faces\n";
  std::cout << "zero-pair C6 quotient ratios {+-2,+-6}; projective cycle 1-2-4-5-3-6-1\n";
  std::cout << "tournament observable zero-obligation C6; switch each residue to positive projective representative\n";
  std::cout << "tournament tie path lexicographically first Hamiltonian path in K6\\C6\n";
  std::cout << "tournament joint fingerprints 5; scores 1,2,2,3,3,4; triangles 6; SCC 6\n";
  std::cout << "tournament flips {2:8,3:52,4:4}; Hamiltonian paths {29:32,31:20,37:12}\n";
  std::cout << "challenged vertices projective owner-obligations preserve the nerve; runners and completed tournaments destroy intersection data\n";
  std::cout << "scalar pattern histogram\n";
  for (const auto &[key, value] : scalar_patterns)
    std::cout << "  " << key << ": " << value << '\n';
  std::cout << "owner-local pattern histogram\n";
  for (const auto &[key, value] : local_patterns)
    std::cout << "  " << key << ": " << value
              << " contexts, " << local_pattern_words[key] << " state words\n";
  std::cout << "local mask table at owner one (state: ratios 1..12 in hex)\n";
  for (int state = 0; state < 10; ++state) {
    std::cout << "  (" << ORDERS[state] << ',' << UNITS[state] << "):";
    for (int ratio = 1; ratio < P; ++ratio)
      std::cout << ' ' << std::hex << MASK[ratio][state][1] << std::dec;
    std::cout << '\n';
  }
}
