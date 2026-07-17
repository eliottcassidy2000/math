#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// Exact exploratory census for the AP-centred, common-scale c=8,
// Hamming-six sheet bank.  This is deliberately a frontier scout rather than
// a theorem-numbered certificate: it separates the hereditary-lcm order
// grammar, scalar owner capacity, owner-local unit feasibility, and the
// globally compatible literal sheet predicate.

using StateWord = std::array<uint8_t, 6>;
using Labels = std::array<uint8_t, 6>;

static constexpr int P = 13;
static constexpr int C = 8;
static constexpr std::array<int, 8> ORDERS{1, 2, 4, 4, 8, 8, 8, 8};
static constexpr std::array<int, 8> UNITS{0, 1, 1, 3, 1, 3, 5, 7};
static constexpr std::array<int, 4> ORDER_REP{0, 1, 2, 4};
static std::array<std::array<std::array<uint8_t, 13>, 8>, 13> MASK{};

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
    if ((value * candidate) % P == 1) return candidate;
  fail("nonunit modulo 13");
}

int crt_base(int label, int state) {
  const int order = ORDERS[state], unit = UNITS[state];
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P && value % order == unit % order)
      return value;
  fail("CRT base missing");
}

uint8_t local_mask(int label, int state, int owner) {
  const int order = ORDERS[state], base = crt_base(label, state);
  const int inverse = inverse_mod_13(owner);
  uint8_t answer = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    if (-order < value && value <= order)
      answer |= static_cast<uint8_t>(1U << sheet);
  }
  return answer;
}

int order_pattern_key(const StateWord &word) {
  std::array<int, 4> counts{};
  for (uint8_t state : word) {
    const int order = ORDERS[state];
    ++counts[order == 1 ? 0 : order == 2 ? 1 : order == 4 ? 2 : 3];
  }
  return counts[0] + 7 * counts[1] + 49 * counts[2] + 343 * counts[3];
}

std::string order_pattern_string(int key) {
  std::array<int, 4> counts{};
  for (int i = 0; i < 4; ++i) {
    counts[i] = key % 7;
    key /= 7;
  }
  static constexpr std::array<int, 4> values{1, 2, 4, 8};
  std::ostringstream out;
  bool first = true;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < counts[i]; ++j) {
      if (!first) out << ',';
      first = false;
      out << values[i];
    }
  }
  return out.str();
}

const std::vector<int> &states_for_order_index(int order_index) {
  static const std::vector<int> d1{0};
  static const std::vector<int> d2{1};
  static const std::vector<int> d4{2, 3};
  static const std::vector<int> d8{4, 5, 6, 7};
  switch (order_index) {
    case 0: return d1;
    case 1: return d2;
    case 2: return d4;
    case 3: return d8;
    default: fail("bad order index");
  }
}

bool owner_locally_feasible(const Labels &labels,
                            const std::array<uint8_t, 6> &order_word,
                            int owner) {
  // The owner-local unit word is allowed to depend on the owner.  Dynamic
  // programming in the 256-element mask semilattice avoids a unit-product
  // loop and is an exact necessary relaxation of global compatibility.
  std::array<bool, 256> reachable{};
  reachable[0] = true;
  for (int provider = 0; provider < 6; ++provider) {
    std::array<bool, 256> next{};
    for (int partial = 0; partial < 256; ++partial) {
      if (!reachable[partial]) continue;
      for (int state : states_for_order_index(order_word[provider]))
        next[partial | MASK[labels[provider]][state][owner]] = true;
    }
    reachable = next;
  }
  return reachable[255];
}

std::set<std::pair<int, int>> doubling_edges(const Labels &labels) {
  std::set<std::pair<int, int>> edges;
  for (int provider : labels)
    for (int owner : labels)
      if (provider != owner &&
          (owner * inverse_mod_13(provider) % P == 2 ||
           owner * inverse_mod_13(provider) % P == 11))
        edges.insert({provider, owner});
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

std::string cyclic_gap_shape(uint8_t owner_subset,
                             const Labels &labels,
                             const std::array<int, 6> &cycle) {
  std::vector<int> positions;
  for (int position = 0; position < 6; ++position) {
    const auto found = std::find(labels.begin(), labels.end(), cycle[position]);
    const int label_index = static_cast<int>(found - labels.begin());
    if (owner_subset & (1U << label_index)) positions.push_back(position);
  }
  if (positions.empty()) return "empty";
  std::vector<int> gaps;
  for (size_t i = 0; i < positions.size(); ++i)
    gaps.push_back((positions[(i + 1) % positions.size()] - positions[i] + 6) % 6);
  std::vector<int> best;
  for (int reflection = 0; reflection < 2; ++reflection) {
    std::vector<int> candidate = gaps;
    if (reflection) std::reverse(candidate.begin(), candidate.end());
    for (size_t shift = 0; shift < candidate.size(); ++shift) {
      std::rotate(candidate.begin(), candidate.begin() + 1, candidate.end());
      if (best.empty() || candidate < best) best = candidate;
    }
  }
  std::ostringstream out;
  for (size_t i = 0; i < best.size(); ++i) {
    if (i) out << ',';
    out << best[i];
  }
  return out.str();
}

int mod4(int value) {
  value %= 4;
  return value < 0 ? value + 4 : value;
}

bool symbolic_all_d8_owner_cover(int owner_position,
                                 const std::array<int, 6> &unit_symbols,
                                 const std::array<int, 6> &edge_signs) {
  // Write the odd unit as e=2x+1 with x in Z/4.  Relative to a signed
  // doubling cycle, an owner receives: fixed sheet 7 from distance 5; fixed
  // sheet 3 plus one even sheet from itself; one even singleton from distance
  // 1; one odd+even pair from distance 2; one even singleton from distance 3;
  // and one odd singleton from distance 4.  Thus literal coverage is exactly
  // the conjunction below: the two odd sheets differ, and the four even-sheet
  // symbols are all distinct in Z/4.
  const auto at = [](const auto &values, int index) {
    return values[(index % 6 + 6) % 6];
  };
  const int i = owner_position;
  const int n1 = at(edge_signs, i);
  const int n2 = n1 ^ at(edge_signs, i + 1);
  const int n3_power = n2 ^ at(edge_signs, i + 2);
  const int n4 = n3_power ^ at(edge_signs, i + 3);

  const int odd_from_distance_2 = (at(unit_symbols, i + 2) & 1) ^ n2;
  const int odd_from_distance_4 = (at(unit_symbols, i + 4) & 1) ^ n4;
  if (odd_from_distance_2 == odd_from_distance_4) return false;

  std::array<int, 4> even_symbols{
      mod4(1 - at(unit_symbols, i)),
      n1 ? mod4(at(unit_symbols, i + 1) + 2)
         : mod4(1 - at(unit_symbols, i + 1)),
      n2 ? mod4(3 - at(unit_symbols, i + 2))
         : at(unit_symbols, i + 2),
      n3_power ? at(unit_symbols, i + 3)
               : mod4(3 - at(unit_symbols, i + 3)),
  };
  std::sort(even_symbols.begin(), even_symbols.end());
  return even_symbols == std::array<int, 4>{0, 1, 2, 3};
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
  return false;
}

using TournamentFingerprint =
    std::tuple<std::array<int, 6>, int, bool, int, int>;

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
  require(found, "signed-cycle absence graph lacks a tie Hamiltonian path");

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
          int out = 0;
          for (int other : triple) out += tournament.count({vertex, other});
          directed_cycle &= out == 1;
        }
        triangles += directed_cycle;
      }
  bool strong = true;
  for (int a : labels)
    for (int b : labels)
      if (a != b) strong &= reachable(a, b, labels, tournament);

  int hamiltonian_paths = 0;
  Labels path = labels;
  do {
    bool directed = true;
    for (int i = 0; i < 5; ++i)
      directed &= tournament.count({path[i], path[i + 1]});
    hamiltonian_paths += directed;
  } while (std::next_permutation(path.begin(), path.end()));
  return {scores, triangles, strong, sparse_flips, hamiltonian_paths};
}

int main() {
  for (int label = 1; label < P; ++label)
    for (int state = 0; state < 8; ++state)
      for (int owner = 1; owner < P; ++owner)
        MASK[label][state][owner] = local_mask(label, state, owner);
  for (int label = 1; label < P; ++label)
    for (int owner = 1; owner < P; ++owner) {
      require(std::popcount(MASK[label][2][owner]) ==
                  std::popcount(MASK[label][3][owner]),
              "D=4 mask cardinality depends on the unit");
      for (int state = 5; state < 8; ++state)
        require(std::popcount(MASK[label][4][owner]) ==
                    std::popcount(MASK[label][state][owner]),
                "D=8 mask cardinality depends on the unit");
    }

  // Freeze the exact D|8 state grammar and local owner-one masks.
  std::map<std::pair<int, int>, std::map<int, std::vector<int>>> local_table;
  for (int state = 0; state < 8; ++state)
    for (int label = 1; label < P; ++label)
      local_table[{ORDERS[state], UNITS[state]}][MASK[label][state][1]]
          .push_back(label);

  // Hereditary lcm(D_j:j!=i)=8 is equivalent to at least two D=8 entries.
  struct StateEntry {
    StateWord word;
    int pattern;
  };
  std::vector<StateEntry> state_words;
  std::map<int, uint64_t> state_words_by_d8;
  for (int code = 0; code < (1 << 18); ++code) {  // 8^6
    int work = code, d8_count = 0;
    StateWord word{};
    for (int i = 0; i < 6; ++i) {
      word[i] = work & 7;
      work >>= 3;
      d8_count += ORDERS[word[i]] == 8;
    }
    if (d8_count < 2) continue;
    state_words.push_back({word, order_pattern_key(word)});
    ++state_words_by_d8[d8_count];
  }
  require(state_words.size() == 233472, "unexpected admissible state-word count");

  std::vector<std::array<uint8_t, 6>> order_words;
  std::map<int, uint64_t> order_words_by_d8;
  for (int code = 0; code < 4096; ++code) {  // 4^6
    int work = code, d8_count = 0;
    std::array<uint8_t, 6> word{};
    for (int i = 0; i < 6; ++i) {
      word[i] = work & 3;
      work >>= 2;
      d8_count += word[i] == 3;
    }
    if (d8_count < 2) continue;
    order_words.push_back(word);
    ++order_words_by_d8[d8_count];
  }
  require(order_words.size() == 1909, "unexpected admissible order-word count");

  uint64_t order_contexts = 0, capacity_contexts = 0, owner_local_contexts = 0;
  uint64_t literal_contexts = 0, literal_survivors = 0;
  std::map<int, uint64_t> order_by_pattern, capacity_by_pattern;
  std::map<int, uint64_t> owner_local_by_pattern;
  std::map<int, uint64_t> tested_by_pattern, survived_by_pattern;
  std::map<Labels, uint64_t> survivor_support_multiplicity;
  std::map<StateWord, uint64_t> survivor_state_words;
  std::set<Labels> capacity_all_d8_supports;
  std::set<Labels> owner_local_all_d8_supports;

  Labels labels{};
  for (labels[0] = 1; labels[0] <= 7; ++labels[0])
  for (labels[1] = labels[0] + 1; labels[1] <= 8; ++labels[1])
  for (labels[2] = labels[1] + 1; labels[2] <= 9; ++labels[2])
  for (labels[3] = labels[2] + 1; labels[3] <= 10; ++labels[3])
  for (labels[4] = labels[3] + 1; labels[4] <= 11; ++labels[4])
  for (labels[5] = labels[4] + 1; labels[5] <= 12; ++labels[5]) {
    std::array<std::array<uint64_t, 8>, 6> packed{};
    for (int provider = 0; provider < 6; ++provider)
      for (int state = 0; state < 8; ++state)
        for (int owner = 0; owner < 6; ++owner)
          packed[provider][state] |=
              static_cast<uint64_t>(MASK[labels[provider]][state][labels[owner]])
              << (C * owner);

    for (const auto &order_word : order_words) {
      StateWord representatives{};
      for (int i = 0; i < 6; ++i) representatives[i] = ORDER_REP[order_word[i]];
      const int pattern = order_pattern_key(representatives);
      ++order_contexts;
      ++order_by_pattern[pattern];

      bool capacity = true;
      for (int owner = 0; owner < 6 && capacity; ++owner) {
        int total = 0;
        for (int provider = 0; provider < 6; ++provider)
          total += std::popcount(
              MASK[labels[provider]][representatives[provider]][labels[owner]]);
        capacity = total >= C;
      }
      if (!capacity) continue;
      ++capacity_contexts;
      ++capacity_by_pattern[pattern];
      if (pattern == 6 * 343) capacity_all_d8_supports.insert(labels);

      bool owner_local = true;
      for (int owner = 0; owner < 6 && owner_local; ++owner)
        owner_local = owner_locally_feasible(labels, order_word, labels[owner]);
      if (!owner_local) continue;
      ++owner_local_contexts;
      ++owner_local_by_pattern[pattern];
      if (pattern == 6 * 343) owner_local_all_d8_supports.insert(labels);
    }

    constexpr uint64_t full = (uint64_t{1} << 48) - 1;
    for (const StateEntry &entry : state_words) {
      const StateWord &word = entry.word;
      const int pattern = entry.pattern;
      ++literal_contexts;
      ++tested_by_pattern[pattern];
      uint64_t cover = 0;
      for (int provider = 0; provider < 6; ++provider)
        cover |= packed[provider][word[provider]];
      if (cover != full) continue;
      ++literal_survivors;
      ++survived_by_pattern[pattern];
      ++survivor_support_multiplicity[labels];
      ++survivor_state_words[word];
    }
  }

  require(order_contexts == 924ULL * order_words.size(),
          "order-context count mismatch");
  require(literal_contexts == 924ULL * state_words.size(),
          "literal-context count mismatch");
  require(capacity_contexts == 3166,
          "unexpected scalar-capacity context count");
  require(owner_local_contexts == 64,
          "unexpected owner-local context count");
  require(literal_survivors == 0 && survivor_support_multiplicity.empty() &&
              survivor_state_words.empty(),
          "a literal scale-eight common-sheet context survived");

  std::set<Labels> signed_cycle_supports;
  Labels candidate{};
  for (candidate[0] = 1; candidate[0] <= 7; ++candidate[0])
  for (candidate[1] = candidate[0] + 1; candidate[1] <= 8; ++candidate[1])
  for (candidate[2] = candidate[1] + 1; candidate[2] <= 9; ++candidate[2])
  for (candidate[3] = candidate[2] + 1; candidate[3] <= 10; ++candidate[3])
  for (candidate[4] = candidate[3] + 1; candidate[4] <= 11; ++candidate[4])
  for (candidate[5] = candidate[4] + 1; candidate[5] <= 12; ++candidate[5])
    if (signed_doubling_cycle(candidate)) signed_cycle_supports.insert(candidate);
  require(owner_local_all_d8_supports.size() == 64,
          "unexpected all-D8 owner-local support count");
  require(owner_local_all_d8_supports == signed_cycle_supports,
          "all-D8 owner-local supports differ from signed doubling C6 bank");
  require(capacity_all_d8_supports.size() == 66,
          "unexpected all-D8 scalar-capacity support count");
  std::set<Labels> capacity_only_supports;
  std::set_difference(capacity_all_d8_supports.begin(), capacity_all_d8_supports.end(),
                      owner_local_all_d8_supports.begin(),
                      owner_local_all_d8_supports.end(),
                      std::inserter(capacity_only_supports,
                                    capacity_only_supports.end()));
  const std::set<Labels> quadratic_cosets{
      Labels{1, 3, 4, 9, 10, 12}, Labels{2, 5, 6, 7, 8, 11}};
  require(capacity_only_supports == quadratic_cosets,
          "the two all-D8 capacity-only supports are not the quadratic cosets");

  // Resolve the 64-support global unit CSP completely and retain its nerve.
  std::map<std::array<int, 7>, int> satisfaction_profile_variants;
  std::map<std::array<int, 7>, int> minimal_obstruction_size_variants;
  std::map<std::array<int, 4>, int> exact_pair_distance_variants;
  std::map<std::string, uint64_t> minimal_obstruction_gap_shapes;
  int maximum_satisfied_owners = 0;
  std::map<TournamentFingerprint, int> tournament_fingerprints;
  for (const Labels &support : owner_local_all_d8_supports) {
    ++tournament_fingerprints[tournament_fingerprint(support)];
    std::array<int, 7> satisfaction_profile{};
    std::array<int, 4> exact_pairs_by_cycle_distance{};
    std::array<int, 6> owner_witnesses{};
    std::array<std::array<int, 6>, 6> pair_witnesses{};
    std::array<bool, 64> compatible{};
    compatible[0] = true;
    const auto cycle = doubling_cycle_order(support);
    std::array<int, 6> cycle_position_by_label_index{};
    std::array<int, 6> edge_signs{};
    for (int position = 0; position < 6; ++position) {
      const int label_index = static_cast<int>(
          std::find(support.begin(), support.end(), cycle[position]) - support.begin());
      cycle_position_by_label_index[label_index] = position;
      const int next = cycle[(position + 1) % 6];
      edge_signs[position] = next * inverse_mod_13(cycle[position]) % P == 11;
    }
    require(std::accumulate(edge_signs.begin(), edge_signs.end(), 0,
                            std::bit_xor<int>()) == 1,
            "signed doubling cycle does not have odd sign parity");
    for (int code = 0; code < 4096; ++code) {
      int work = code;
      StateWord word{};
      for (int provider = 0; provider < 6; ++provider) {
        word[provider] = static_cast<uint8_t>(4 + (work & 3));
        work >>= 2;
      }
      std::array<int, 6> unit_symbols_by_cycle_position{};
      for (int position = 0; position < 6; ++position) {
        const int label_index = static_cast<int>(
            std::find(support.begin(), support.end(), cycle[position]) - support.begin());
        unit_symbols_by_cycle_position[position] = word[label_index] - 4;
      }
      uint8_t owner_mask = 0;
      for (int owner = 0; owner < 6; ++owner) {
        uint8_t cover = 0;
        for (int provider = 0; provider < 6; ++provider)
          cover |= MASK[support[provider]][word[provider]][support[owner]];
        const bool literal = cover == 255;
        const bool symbolic = symbolic_all_d8_owner_cover(
            cycle_position_by_label_index[owner], unit_symbols_by_cycle_position,
            edge_signs);
        require(literal == symbolic,
                "symbolic all-D8 odd/all-different grammar mismatch");
        if (literal) owner_mask |= static_cast<uint8_t>(1U << owner);
      }
      const int satisfied = std::popcount(owner_mask);
      ++satisfaction_profile[satisfied];
      maximum_satisfied_owners = std::max(maximum_satisfied_owners, satisfied);
      for (int owner = 0; owner < 6; ++owner)
        owner_witnesses[owner] += (owner_mask >> owner) & 1U;
      if (satisfied == 2) {
        int first = -1, second = -1;
        for (int owner = 0; owner < 6; ++owner)
          if (owner_mask & (1U << owner)) {
            if (first < 0) first = owner;
            else second = owner;
          }
        int distance = std::abs(cycle_position_by_label_index[first] -
                                cycle_position_by_label_index[second]);
        distance = std::min(distance, 6 - distance);
        ++exact_pairs_by_cycle_distance[distance];
        ++pair_witnesses[first][second];
        ++pair_witnesses[second][first];
      }
      for (uint8_t subset = owner_mask;; subset = (subset - 1) & owner_mask) {
        compatible[subset] = true;
        if (subset == 0) break;
      }
    }
    ++satisfaction_profile_variants[satisfaction_profile];
    ++exact_pair_distance_variants[exact_pairs_by_cycle_distance];
    require(std::all_of(owner_witnesses.begin(), owner_witnesses.end(),
                        [](int count) { return count == 192; }),
            "all-D8 owner obligations do not each have 192 unit witnesses");
    for (int first = 0; first < 6; ++first)
      for (int second = first + 1; second < 6; ++second) {
        int distance = std::abs(cycle_position_by_label_index[first] -
                                cycle_position_by_label_index[second]);
        distance = std::min(distance, 6 - distance);
        const int expected = distance == 1 ? 8 : distance == 2 ? 0 : 16;
        require(pair_witnesses[first][second] == expected,
                "all-D8 owner-pair intersection size mismatch");
      }

    std::array<int, 7> minimal_by_size{};
    for (uint8_t subset = 1; subset < 64; ++subset) {
      if (compatible[subset]) continue;
      bool minimal = true;
      for (int owner = 0; owner < 6; ++owner)
        if ((subset & (1U << owner)) &&
            !compatible[subset ^ static_cast<uint8_t>(1U << owner)])
          minimal = false;
      if (!minimal) continue;
      ++minimal_by_size[std::popcount(subset)];
      ++minimal_obstruction_gap_shapes[cyclic_gap_shape(subset, support, cycle)];
    }
    ++minimal_obstruction_size_variants[minimal_by_size];
  }
  std::array<int, 7> expected_satisfaction_profile{};
  expected_satisfaction_profile[0] = 3040;
  expected_satisfaction_profile[1] = 960;
  expected_satisfaction_profile[2] = 96;
  require(maximum_satisfied_owners == 2,
          "an all-D8 unit word satisfies more than two owners");
  require(satisfaction_profile_variants ==
              std::map<std::array<int, 7>, int>{{expected_satisfaction_profile, 64}},
          "all-D8 owner-satisfaction profile mismatch");
  std::array<int, 4> expected_pair_distances{};
  expected_pair_distances[1] = 48;
  expected_pair_distances[3] = 48;
  require(exact_pair_distance_variants ==
              std::map<std::array<int, 4>, int>{{expected_pair_distances, 64}},
          "all-D8 exact-two-owner distance profile mismatch");
  std::array<int, 7> expected_minimal_obstructions{};
  expected_minimal_obstructions[2] = 6;
  require(minimal_obstruction_size_variants ==
              std::map<std::array<int, 7>, int>{{expected_minimal_obstructions, 64}},
          "all-D8 minimal-obstruction size profile mismatch");
  require(minimal_obstruction_gap_shapes == std::map<std::string, uint64_t>{{"2,4", 384}},
          "all-D8 minimal-obstruction cycle shape mismatch");

  std::cout << "SCALE-EIGHT HAMMING-SIX COMMON-SHEET FRONTIER SCOUT\n";
  std::cout << "states=(1,0),(2,1),(4,1),(4,3),(8,1),(8,3),(8,5),(8,7)\n";
  std::cout << "hereditary_lcm_criterion=D8_count_at_least_2\n";
  std::cout << "admissible_order_words=" << order_words.size() << '\n';
  std::cout << "admissible_state_words=" << state_words.size() << '\n';
  std::cout << "order_words_by_D8_count";
  for (const auto &[count, number] : order_words_by_d8)
    std::cout << ' ' << count << ':' << number;
  std::cout << '\n';
  std::cout << "state_words_by_D8_count";
  for (const auto &[count, number] : state_words_by_d8)
    std::cout << ' ' << count << ':' << number;
  std::cout << '\n';
  std::cout << "local_owner_one_table\n";
  for (const auto &[state, groups] : local_table) {
    std::cout << "D=" << state.first << ",e=" << state.second;
    for (const auto &[mask, values] : groups) {
      std::cout << " " << std::hex << mask << std::dec << '=';
      for (size_t i = 0; i < values.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << values[i];
      }
    }
    std::cout << '\n';
  }
  std::cout << "order_contexts=" << order_contexts << '\n';
  std::cout << "scalar_capacity_contexts=" << capacity_contexts << '\n';
  std::cout << "owner_locally_feasible_order_contexts=" << owner_local_contexts << '\n';
  std::cout << "all_D8_scalar_capacity_supports=66=64_signed_cycles_plus_2_quadratic_cosets\n";
  std::cout << "all_D8_quadratic_cosets_fail_owner_local_units=1\n";
  std::cout << "owner_local_supports_equal_signed_doubling_C6="
            << (owner_local_all_d8_supports == signed_cycle_supports) << '\n';
  std::cout << "literal_contexts=" << literal_contexts << '\n';
  std::cout << "literal_common_sheet_survivors=" << literal_survivors << '\n';
  std::cout << "surviving_supports=" << survivor_support_multiplicity.size() << '\n';
  std::cout << "surviving_state_words=" << survivor_state_words.size() << '\n';
  std::cout << "all_D8_maximum_simultaneously_satisfied_owners="
            << maximum_satisfied_owners << '\n';
  std::cout << "all_D8_each_owner_unit_witnesses=192\n";
  std::cout << "all_D8_owner_grammar=odd_complement_plus_four_even_symbols_all_different\n";
  std::cout << "all_D8_pair_nerve=K3,3_cycle_distance_1_or_3\n";
  std::cout << "all_D8_pair_intersections=distance1:8,distance2:0,distance3:16\n";
  std::cout << "all_D8_satisfaction_profile_variants\n";
  for (const auto &[profile, support_count] : satisfaction_profile_variants) {
    std::cout << "supports=" << support_count;
    for (int satisfied = 0; satisfied <= 6; ++satisfied)
      if (profile[satisfied])
        std::cout << " sat" << satisfied << '=' << profile[satisfied];
    std::cout << '\n';
  }
  std::cout << "all_D8_exact_two_owner_cycle_distance_variants\n";
  for (const auto &[profile, support_count] : exact_pair_distance_variants) {
    std::cout << "supports=" << support_count;
    for (int distance = 1; distance <= 3; ++distance)
      std::cout << " distance" << distance << '=' << profile[distance];
    std::cout << '\n';
  }
  std::cout << "all_D8_minimal_obstruction_size_variants\n";
  for (const auto &[profile, support_count] : minimal_obstruction_size_variants) {
    std::cout << "supports=" << support_count;
    for (int size = 1; size <= 6; ++size)
      if (profile[size]) std::cout << " size" << size << '=' << profile[size];
    std::cout << '\n';
  }
  std::cout << "all_D8_minimal_obstruction_cyclic_gap_shapes";
  for (const auto &[shape, count] : minimal_obstruction_gap_shapes)
    std::cout << " (" << shape << "):" << count;
  std::cout << '\n';
  std::map<int, int> tournament_flip_histogram, tournament_path_histogram;
  for (const auto &[fingerprint, count] : tournament_fingerprints) {
    const auto &[scores, triangles, strong, flips, paths] = fingerprint;
    require(scores == std::array<int, 6>{1, 2, 2, 3, 3, 4},
            "tournament score fingerprint mismatch");
    require(triangles == 6 && strong,
            "tournament triangle/SCC fingerprint mismatch");
    tournament_flip_histogram[flips] += count;
    tournament_path_histogram[paths] += count;
  }
  require(tournament_flip_histogram == std::map<int, int>{{2, 8}, {3, 52}, {4, 4}},
          "tournament sparse-flip histogram mismatch");
  require(tournament_path_histogram == std::map<int, int>{{29, 32}, {31, 20}, {37, 12}},
          "tournament Hamiltonian-path histogram mismatch");
  std::cout << "tournament_pair_observable=signed_doubling_ratio_plus_or_minus_2\n";
  std::cout << "tournament_switch=edge_sign_with_absent_pairs_tied\n";
  std::cout << "tournament_tie_path=lexicographically_first_absence_graph_Hamiltonian_path\n";
  std::cout << "tournament_joint_fingerprints=" << tournament_fingerprints.size() << '\n';
  std::cout << "tournament_scores=1,2,2,3,3,4 triangles=6 SCC=6\n";
  std::cout << "tournament_sparse_flip_histogram";
  for (const auto &[value, count] : tournament_flip_histogram)
    std::cout << ' ' << value << ':' << count;
  std::cout << '\n';
  std::cout << "tournament_Hamiltonian_path_histogram";
  for (const auto &[value, count] : tournament_path_histogram)
    std::cout << ' ' << value << ':' << count;
  std::cout << '\n';
  std::cout << "challenged_vertices=owner_obligations_preserve_K3,3_nerve_bare_tournament_destroys_it\n";
  std::cout << "patterns_with_owner_local_or_literal_survivors\n";
  std::set<int> patterns;
  for (const auto &[pattern, count] : owner_local_by_pattern) patterns.insert(pattern);
  for (const auto &[pattern, count] : survived_by_pattern) patterns.insert(pattern);
  for (int pattern : patterns)
    std::cout << order_pattern_string(pattern) << " order=" << order_by_pattern[pattern]
              << " capacity=" << capacity_by_pattern[pattern]
              << " owner_local=" << owner_local_by_pattern[pattern]
              << " literal_tested=" << tested_by_pattern[pattern]
              << " literal_survived=" << survived_by_pattern[pattern] << '\n';
  std::cout << "patterns_with_scalar_capacity\n";
  for (const auto &[pattern, count] : capacity_by_pattern)
    std::cout << order_pattern_string(pattern) << " order=" << order_by_pattern[pattern]
              << " capacity=" << count
              << " owner_local=" << owner_local_by_pattern[pattern] << '\n';
}
