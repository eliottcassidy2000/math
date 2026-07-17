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

// Exact frontier scout for the primitive proper AP-centred common-scale-nine
// Hamming-six sheet bank.  It separates hereditary-lcm order grammar, scalar
// capacity, owner-local feasibility (units may vary by owner), and globally
// compatible literal sheet coverage.

using Labels = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;
using StateWord = std::array<uint8_t, 6>;

static constexpr int P = 13;
static constexpr int C = 9;
static constexpr std::array<int, 9> ORDERS{1, 3, 3, 9, 9, 9, 9, 9, 9};
static constexpr std::array<int, 9> UNITS{0, 1, 2, 1, 2, 4, 5, 7, 8};
static constexpr std::array<int, 3> ORDER_REP{0, 1, 3};
static std::array<std::array<std::array<uint16_t, 13>, 9>, 13> MASK{};

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
  static const std::vector<int> d3{1, 2};
  static const std::vector<int> d9{3, 4, 5, 6, 7, 8};
  if (order_index == 0) return d1;
  if (order_index == 1) return d3;
  if (order_index == 2) return d9;
  fail("bad order index");
}

bool owner_locally_feasible(const Labels &labels, const OrderWord &orders,
                            int owner) {
  std::array<bool, 512> reachable{};
  reachable[0] = true;
  for (int provider = 0; provider < 6; ++provider) {
    std::array<bool, 512> next{};
    for (int partial = 0; partial < 512; ++partial) {
      if (!reachable[partial]) continue;
      for (int state : states_for_order(orders[provider]))
        next[partial | MASK[labels[provider]][state][owner]] = true;
    }
    reachable = next;
  }
  return reachable[511];
}

bool scalar_capacity(const Labels &labels, const OrderWord &orders, int owner) {
  int capacity = 0;
  for (int provider = 0; provider < 6; ++provider)
    capacity += std::popcount(
        MASK[labels[provider]][ORDER_REP[orders[provider]]][owner]);
  return capacity >= C;
}

uint64_t state_fibre_size(const OrderWord &orders) {
  uint64_t answer = 1;
  for (int order : orders) answer *= states_for_order(order).size();
  return answer;
}

std::string pattern(const OrderWord &orders) {
  std::array<int, 3> counts{};
  for (int order : orders) ++counts[order];
  std::ostringstream out;
  out << "1^" << counts[0] << " 3^" << counts[1] << " 9^" << counts[2];
  return out.str();
}

uint8_t owner_satisfaction_mask(const Labels &labels, const StateWord &states) {
  uint8_t answer = 0;
  for (int owner = 0; owner < 6; ++owner) {
    uint16_t cover = 0;
    for (int provider = 0; provider < 6; ++provider)
      cover |= MASK[labels[provider]][states[provider]][labels[owner]];
    if (cover == 511) answer |= static_cast<uint8_t>(1U << owner);
  }
  return answer;
}

uint64_t enumerate_global_fibre(const Labels &labels, const OrderWord &orders,
                                StateWord &states, int coordinate,
                                uint64_t &tested,
                                std::array<uint64_t, 64> &profile) {
  if (coordinate == 6) {
    ++tested;
    const uint8_t owners = owner_satisfaction_mask(labels, states);
    ++profile[owners];
    return owners == 63 ? 1 : 0;
  }
  uint64_t survivors = 0;
  for (int state : states_for_order(orders[coordinate])) {
    states[coordinate] = static_cast<uint8_t>(state);
    survivors += enumerate_global_fibre(labels, orders, states, coordinate + 1,
                                        tested, profile);
  }
  return survivors;
}

int ratio(int provider, int owner) {
  return provider * inverse_mod_13(owner) % P;
}

std::string support_string(const Labels &labels) {
  std::ostringstream out;
  for (int i = 0; i < 6; ++i) {
    if (i) out << ',';
    out << static_cast<int>(labels[i]);
  }
  return out.str();
}

std::set<std::pair<int, int>> doubling_edges(const Labels &labels) {
  std::set<std::pair<int, int>> edges;
  for (int provider : labels)
    for (int owner : labels)
      if (provider != owner) {
        const int relative = owner * inverse_mod_13(provider) % P;
        if (relative == 2 || relative == 11) edges.insert({provider, owner});
      }
  return edges;
}

bool signed_doubling_cycle(const Labels &labels) {
  const auto edges = doubling_edges(labels);
  if (edges.size() != 6) return false;
  for (int vertex : labels) {
    int indegree = 0, outdegree = 0;
    for (const auto &[source, target] : edges) {
      indegree += target == vertex;
      outdegree += source == vertex;
    }
    if (indegree != 1 || outdegree != 1) return false;
  }
  int current = labels[0], visited = 0;
  do {
    ++visited;
    current = std::find_if(edges.begin(), edges.end(),
                           [current](const auto &edge) {
                             return edge.first == current;
                           })->second;
  } while (current != labels[0] && visited <= 6);
  return visited == 6 && current == labels[0];
}

std::array<int, 6> doubling_cycle_order(const Labels &labels) {
  require(signed_doubling_cycle(labels), "support is not a signed doubling cycle");
  const auto edges = doubling_edges(labels);
  std::array<int, 6> cycle{};
  cycle[0] = labels[0];
  for (int i = 1; i < 6; ++i)
    cycle[i] = std::find_if(edges.begin(), edges.end(),
                            [previous = cycle[i - 1]](const auto &edge) {
                              return edge.first == previous;
                            })->second;
  return cycle;
}

bool mixed_orbit_context(const Labels &labels, const OrderWord &orders) {
  std::set<int> d3, d9;
  for (int i = 0; i < 6; ++i)
    (orders[i] == 1 ? d3 : d9).insert(labels[i]);
  if (d3.size() != 2 || d9.size() != 4) return false;
  for (int anchor : d3) {
    const std::set<int> expected_d3{anchor, 5 * anchor % P};
    const std::set<int> expected_d9{2 * anchor % P, 11 * anchor % P,
                                    3 * anchor % P, 10 * anchor % P};
    if (d3 == expected_d3 && d9 == expected_d9) return true;
  }
  return false;
}

std::array<uint64_t, 7> satisfaction_histogram(
    const std::array<uint64_t, 64> &profile) {
  std::array<uint64_t, 7> answer{};
  for (int owners = 0; owners < 64; ++owners)
    answer[std::popcount(static_cast<unsigned>(owners))] += profile[owners];
  return answer;
}

std::array<uint64_t, 6> owner_intersections(
    const std::array<uint64_t, 64> &profile) {
  std::array<uint64_t, 6> answer{};
  for (int i = 0; i < 6; ++i)
    for (int owners = 0; owners < 64; ++owners)
      if (owners & (1 << i)) answer[i] += profile[owners];
  return answer;
}

std::array<std::array<uint64_t, 6>, 6> pair_intersections(
    const std::array<uint64_t, 64> &profile) {
  std::array<std::array<uint64_t, 6>, 6> answer{};
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      for (int owners = 0; owners < 64; ++owners)
        if ((owners & (1 << i)) && (owners & (1 << j)))
          answer[i][j] += profile[owners];
  return answer;
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
      if (neighbor == target) return true;
      seen.insert(neighbor);
      stack.push_back(neighbor);
    }
  }
  return source == target;
}

using TournamentFingerprint =
    std::tuple<std::array<int, 6>, int, std::array<int, 6>, int, int>;

TournamentFingerprint tournament_fingerprint(const Labels &labels) {
  const auto sparse = doubling_edges(labels);
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
  require(found, "absence graph lacks a tie Hamiltonian path");

  std::array<int, 13> position{};
  for (int i = 0; i < 6; ++i) position[tie_path[i]] = i;
  std::set<std::pair<int, int>> tournament;
  int sparse_flips = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j) {
      const int a = labels[i], b = labels[j];
      std::pair<int, int> edge;
      if (sparse.count({a, b})) edge = {a, b};
      else if (sparse.count({b, a})) edge = {b, a};
      else edge = position[a] < position[b] ? std::pair{a, b} : std::pair{b, a};
      tournament.insert(edge);
      sparse_flips += sparse.count(edge) && position[edge.first] > position[edge.second];
    }

  std::array<int, 6> scores{};
  for (int i = 0; i < 6; ++i)
    for (int other : labels) scores[i] += tournament.count({labels[i], other});
  std::sort(scores.begin(), scores.end());

  int triangles = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j)
      for (int k = j + 1; k < 6; ++k) {
        const std::array<int, 3> triple{labels[i], labels[j], labels[k]};
        bool directed_cycle = true;
        for (int vertex : triple) {
          int outdegree = 0;
          for (int other : triple) outdegree += tournament.count({vertex, other});
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
    for (int state = 0; state < 9; ++state)
      for (int owner = 1; owner < P; ++owner)
        MASK[label][state][owner] = local_mask(label, state, owner);

  for (int order = 0; order < 3; ++order)
    for (int ratio_value = 1; ratio_value < P; ++ratio_value) {
      const int cardinality = std::popcount(MASK[ratio_value][ORDER_REP[order]][1]);
      for (int state : states_for_order(order))
        require(std::popcount(MASK[ratio_value][state][1]) == cardinality,
                "mask cardinality depends on unit");
    }

  uint64_t supports = 0, order_words = 0, scalar_contexts = 0;
  uint64_t local_contexts = 0, global_state_words_tested = 0;
  uint64_t global_survivors = 0;
  std::map<std::string, uint64_t> scalar_patterns, local_patterns;
  std::map<std::string, uint64_t> local_fibre_words;
  std::set<std::string> scalar_supports, local_supports, global_supports;
  std::map<int, uint64_t> all_d9_special_degree_hist;
  std::vector<std::pair<Labels, OrderWord>> local_bank;
  std::set<Labels> local_all_d9_supports;

  Labels labels{};
  for (int a = 1; a <= 7; ++a)
    for (int b = a + 1; b <= 8; ++b)
      for (int c = b + 1; c <= 9; ++c)
        for (int d = c + 1; d <= 10; ++d)
          for (int e = d + 1; e <= 11; ++e)
            for (int f = e + 1; f <= 12; ++f) {
              labels = {static_cast<uint8_t>(a), static_cast<uint8_t>(b),
                        static_cast<uint8_t>(c), static_cast<uint8_t>(d),
                        static_cast<uint8_t>(e), static_cast<uint8_t>(f)};
              ++supports;
              for (int code = 0; code < 729; ++code) {
                int work = code, d9_count = 0;
                OrderWord orders{};
                for (int i = 0; i < 6; ++i) {
                  orders[i] = static_cast<uint8_t>(work % 3);
                  d9_count += orders[i] == 2;
                  work /= 3;
                }
                if (d9_count < 2) continue;
                ++order_words;
                bool scalar = true;
                for (int owner : labels)
                  scalar &= scalar_capacity(labels, orders, owner);
                if (!scalar) continue;
                ++scalar_contexts;
                scalar_patterns[pattern(orders)]++;
                scalar_supports.insert(support_string(labels));
                bool local = true;
                for (int owner : labels)
                  local &= owner_locally_feasible(labels, orders, owner);
                if (!local) continue;
                ++local_contexts;
                local_patterns[pattern(orders)]++;
                local_fibre_words[pattern(orders)] += state_fibre_size(orders);
                local_supports.insert(support_string(labels));
                local_bank.push_back({labels, orders});
                if (pattern(orders) == "1^0 3^0 9^6")
                  local_all_d9_supports.insert(labels);
              }

              OrderWord all_d9{2, 2, 2, 2, 2, 2};
              bool local_all_d9 = true;
              int minimum_special_degree = 100;
              for (int owner : labels) {
                local_all_d9 &= owner_locally_feasible(labels, all_d9, owner);
                int degree = 0;
                for (int provider : labels)
                  degree += ratio(provider, owner) == 2 || ratio(provider, owner) == 5 ||
                            ratio(provider, owner) == 8 || ratio(provider, owner) == 11;
                minimum_special_degree = std::min(minimum_special_degree, degree);
              }
              if (local_all_d9) ++all_d9_special_degree_hist[minimum_special_degree];
            }

  StateWord states{};
  int all_d9_nerve_contexts = 0, mixed_nerve_contexts = 0;
  std::map<uint64_t, uint64_t> mixed_d3_d9_pair_histogram;
  std::map<TournamentFingerprint, int> all_d9_tournaments, mixed_tournaments;
  for (const auto &[support, orders] : local_bank) {
    uint64_t tested = 0;
    std::array<uint64_t, 64> profile{};
    const uint64_t survivors_here =
        enumerate_global_fibre(support, orders, states, 0, tested, profile);
    require(tested == state_fibre_size(orders), "state fibre count mismatch");
    global_state_words_tested += tested;
    global_survivors += survivors_here;
    if (survivors_here) global_supports.insert(support_string(support));

    const auto histogram = satisfaction_histogram(profile);
    const auto owners = owner_intersections(profile);
    const auto pairs = pair_intersections(profile);
    const bool all_d9 = std::all_of(orders.begin(), orders.end(),
                                    [](int order) { return order == 2; });
    if (all_d9) {
      ++all_d9_nerve_contexts;
      ++all_d9_tournaments[tournament_fingerprint(support)];
      require(signed_doubling_cycle(support),
              "an all-D9 owner-local support is not a signed doubling cycle");
      require(histogram == std::array<uint64_t, 7>{44'100, 2'520, 36, 0, 0, 0, 0},
              "unexpected all-D9 owner-satisfaction profile");
      require(owners == std::array<uint64_t, 6>{432, 432, 432, 432, 432, 432},
              "unexpected all-D9 owner-obligation size");
      const auto cycle = doubling_cycle_order(support);
      std::array<int, 13> position{};
      for (int i = 0; i < 6; ++i) position[cycle[i]] = i;
      for (int i = 0; i < 6; ++i)
        for (int j = i + 1; j < 6; ++j) {
          int distance = std::abs(position[support[i]] - position[support[j]]);
          distance = std::min(distance, 6 - distance);
          require(pairs[i][j] == (distance == 3 ? 12 : 0),
                  "all-D9 pair nerve is not the antipodal matching");
        }
    } else {
      ++mixed_nerve_contexts;
      ++mixed_tournaments[tournament_fingerprint(support)];
      require(mixed_orbit_context(support, orders),
              "mixed owner-local context lies outside the claimed orbit");
      require(histogram == std::array<uint64_t, 7>{2'928, 1'776, 336, 144, 0, 0, 0},
              "unexpected mixed owner-satisfaction profile");
      for (int i = 0; i < 6; ++i)
        require(owners[i] == (orders[i] == 1 ? 1'152 : 144),
                "unexpected mixed owner-obligation size");
      for (int i = 0; i < 6; ++i)
        for (int j = i + 1; j < 6; ++j) {
          if (orders[i] == 1 && orders[j] == 1) {
            require(pairs[i][j] == 192,
                    "mixed D3-D3 pair intersection mismatch");
          } else if (orders[i] == 1 || orders[j] == 1) {
            require(pairs[i][j] == 24 || pairs[i][j] == 48,
                    "mixed D3-D9 pair intersection mismatch");
            ++mixed_d3_d9_pair_histogram[pairs[i][j]];
          } else {
            const bool antipodal = support[i] + support[j] == P;
            require(pairs[i][j] == (antipodal ? 144 : 0),
                    "mixed D9 pair nerve is not two antipodal edges");
          }
        }
    }
  }

  std::set<Labels> signed_cycle_supports;
  Labels candidate{};
  for (candidate[0] = 1; candidate[0] <= 7; ++candidate[0])
  for (candidate[1] = candidate[0] + 1; candidate[1] <= 8; ++candidate[1])
  for (candidate[2] = candidate[1] + 1; candidate[2] <= 9; ++candidate[2])
  for (candidate[3] = candidate[2] + 1; candidate[3] <= 10; ++candidate[3])
  for (candidate[4] = candidate[3] + 1; candidate[4] <= 11; ++candidate[4])
  for (candidate[5] = candidate[4] + 1; candidate[5] <= 12; ++candidate[5])
    if (signed_doubling_cycle(candidate)) signed_cycle_supports.insert(candidate);
  require(local_all_d9_supports == signed_cycle_supports &&
              signed_cycle_supports.size() == 64,
          "the all-D9 local bank is not exactly the 64 signed cycles");
  require(all_d9_nerve_contexts == 64 && mixed_nerve_contexts == 12,
          "unexpected structural context split");
  require(mixed_d3_d9_pair_histogram ==
              std::map<uint64_t, uint64_t>{{24, 48}, {48, 48}},
          "unexpected mixed cross-order pair histogram");

  std::map<int, int> all_d9_flip_histogram, all_d9_path_histogram;
  for (const auto &[fingerprint, multiplicity] : all_d9_tournaments) {
    const auto &[scores, triangles, scc, flips, paths] = fingerprint;
    require(scores == std::array<int, 6>{1, 2, 2, 3, 3, 4} && triangles == 6 &&
                scc == std::array<int, 6>{6, 0, 0, 0, 0, 0},
            "unexpected all-D9 tournament fingerprint");
    all_d9_flip_histogram[flips] += multiplicity;
    all_d9_path_histogram[paths] += multiplicity;
  }
  require(all_d9_tournaments.size() == 5 &&
              all_d9_flip_histogram == std::map<int, int>{{2, 8}, {3, 52}, {4, 4}} &&
              all_d9_path_histogram == std::map<int, int>{{29, 32}, {31, 20}, {37, 12}},
          "all-D9 tournament histogram mismatch");

  std::map<int, int> mixed_triangle_histogram, mixed_flip_histogram;
  std::map<int, int> mixed_path_histogram;
  std::map<std::array<int, 6>, int> mixed_scc_histogram;
  for (const auto &[fingerprint, multiplicity] : mixed_tournaments) {
    const auto &[scores, triangles, scc, flips, paths] = fingerprint;
    (void)scores;
    mixed_triangle_histogram[triangles] += multiplicity;
    mixed_scc_histogram[scc] += multiplicity;
    mixed_flip_histogram[flips] += multiplicity;
    mixed_path_histogram[paths] += multiplicity;
  }
  require(mixed_tournaments.size() == 6 &&
              mixed_triangle_histogram ==
                  std::map<int, int>{{0, 1}, {1, 4}, {3, 2}, {4, 2}, {6, 3}} &&
              mixed_scc_histogram ==
                  std::map<std::array<int, 6>, int>{
                      {{1, 1, 1, 1, 1, 1}, 1},
                      {{3, 1, 1, 1, 0, 0}, 4},
                      {{5, 1, 0, 0, 0, 0}, 2},
                      {{6, 0, 0, 0, 0, 0}, 5}} &&
              mixed_flip_histogram ==
                  std::map<int, int>{{0, 1}, {1, 4}, {2, 2}, {3, 2}, {4, 3}} &&
              mixed_path_histogram ==
                  std::map<int, int>{{1, 1}, {3, 4}, {9, 2}, {17, 2}, {33, 3}},
          "mixed tournament histogram mismatch");

  require(supports == 924 && order_words == 437'052,
          "support/order-context census mismatch");
  require(scalar_contexts == 1'186 && local_contexts == 76,
          "capacity/owner-local census mismatch");
  require(global_state_words_tested == 3'048'192 && global_survivors == 0,
          "reduced global census mismatch");

  std::cout << "scale-nine AP-centred Hamming-six frontier scout\n";
  std::cout << "supports " << supports << '\n';
  std::cout << "hereditary state words per support 521964\n";
  std::cout << "hereditary labelled state contexts 482294736\n";
  std::cout << "hereditary order contexts " << order_words << '\n';
  std::cout << "scalar-capacity contexts " << scalar_contexts << '\n';
  std::cout << "scalar-capacity supports " << scalar_supports.size() << '\n';
  std::cout << "owner-local contexts " << local_contexts << '\n';
  std::cout << "owner-local supports " << local_supports.size() << '\n';
  std::cout << "owner-local patterns\n";
  for (const auto &[key, value] : local_patterns)
    std::cout << "  " << key << ": " << value << " contexts, "
              << local_fibre_words[key] << " state words\n";
  std::cout << "global state words tested after owner-local reduction "
            << global_state_words_tested << '\n';
  std::cout << "global common-sheet survivors " << global_survivors << '\n';
  std::cout << "global surviving supports " << global_supports.size() << '\n';
  std::cout << "all-D9 structure 64 signed-doubling C6 supports\n";
  std::cout << "all-D9 obligation profile 0:44100 1:2520 2:36\n";
  std::cout << "all-D9 owner size 432; pair nerve 3K2 with antipodal intersection 12\n";
  std::cout << "mixed structure orbit {D3:a,5a; D9:+-2a,+-3a}, size 12\n";
  std::cout << "mixed obligation profile 0:2928 1:1776 2:336 3:144\n";
  std::cout << "mixed owner sizes D3:1152 D9:144; D9 pair nerve 2K2\n";
  std::cout << "mixed pair intersections D3-D3:192 D3-D9:{24:48,48:48} D9 antipodes:144\n";
  std::cout << "tournament observable o/p in {+2,-2}; switch signed edge; absent pairs tied\n";
  std::cout << "tournament tie path lexicographically first absence-graph Hamiltonian path\n";
  std::cout << "all-D9 tournament joint fingerprints 5; scores 1,2,2,3,3,4; triangles 6; SCC 6\n";
  std::cout << "all-D9 tournament flips {2:8,3:52,4:4}; Hamiltonian paths {29:32,31:20,37:12}\n";
  std::cout << "mixed tournament joint fingerprints 6; triangle histogram {0:1,1:4,3:2,4:2,6:3}\n";
  std::cout << "mixed tournament SCC {1+1+1+1+1+1:1,3+1+1+1:4,5+1:2,6:5}\n";
  std::cout << "challenged vertices owner obligations preserve 3K2/2K2 nerves; completed tournaments destroy intersection data\n";
  std::cout << "all-D9 owner-local minimum special-degree histogram\n";
  for (const auto &[key, value] : all_d9_special_degree_hist)
    std::cout << "  " << key << ": " << value << '\n';
  std::cout << "scalar pattern histogram\n";
  for (const auto &[key, value] : scalar_patterns)
    std::cout << "  " << key << ": " << value << '\n';
  std::cout << "local mask table at owner one (state: ratios 1..12 in hex)\n";
  for (int state = 0; state < 9; ++state) {
    std::cout << "  (" << ORDERS[state] << ',' << UNITS[state] << "):";
    for (int label = 1; label < P; ++label) {
      std::ostringstream hex;
      hex << std::hex << MASK[label][state][1];
      std::cout << ' ' << hex.str();
    }
    std::cout << '\n';
  }
  return 0;
}
