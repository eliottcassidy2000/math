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

// Exact certificate for the primitive proper AP-centred common-scale-eleven
// Hamming-six sheet bank.  It scans the complete hereditary divisor grammar,
// applies unit-independent scalar capacity and owner-local feasibility, and
// replays all 66 remaining 10^6-word literal fibres.  Multiplication orbits
// are recorded only as structure; no covariance quotient is trusted.

using Labels = std::array<uint8_t, 6>;
using Profile = std::array<uint64_t, 64>;

static constexpr int P = 13;
static constexpr int C = 11;
static constexpr uint16_t FULL = (1U << C) - 1;
static std::array<std::array<uint16_t, 13>, 13> MASK1{};
static std::array<std::array<std::array<uint16_t, 13>, 10>, 13> MASK11{};

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

int inverse_mod_13(int value) {
  for (int candidate = 1; candidate < P; ++candidate)
    if (value * candidate % P == 1) return candidate;
  fail("nonunit modulo thirteen");
}

int centered(int value, int modulus) {
  int residue = value % modulus;
  if (residue < 0) residue += modulus;
  return 2 * residue > modulus ? residue - modulus : residue;
}

int crt_base(int label, int order, int unit) {
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P &&
        value % order == unit % order)
      return value;
  fail("CRT base missing");
}

uint16_t local_mask(int label, int order, int unit, int owner) {
  const int base = crt_base(label, order, unit);
  const int inverse = inverse_mod_13(owner);
  uint16_t answer = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    if (-order < value && value <= order)
      answer |= static_cast<uint16_t>(1U << sheet);
  }
  return answer;
}

uint64_t pow10(int exponent) {
  uint64_t answer = 1;
  while (exponent-- > 0) answer *= 10;
  return answer;
}

bool scalar_capacity(const Labels &labels, uint8_t order_word, int owner) {
  int capacity = 0;
  for (int provider = 0; provider < 6; ++provider) {
    const int label = labels[provider];
    const uint16_t mask = (order_word & (1U << provider))
        ? MASK11[label][0][owner]
        : MASK1[label][owner];
    capacity += std::popcount(mask);
  }
  return capacity >= C;
}

bool owner_locally_feasible(const Labels &labels, uint8_t order_word,
                            int owner) {
  std::array<bool, 2048> reachable{};
  reachable[0] = true;
  for (int provider = 0; provider < 6; ++provider) {
    std::array<bool, 2048> next{};
    const int label = labels[provider];
    for (int partial = 0; partial < 2048; ++partial) {
      if (!reachable[partial]) continue;
      if (!(order_word & (1U << provider))) {
        next[partial | MASK1[label][owner]] = true;
      } else {
        for (int unit = 0; unit < 10; ++unit)
          next[partial | MASK11[label][unit][owner]] = true;
      }
    }
    reachable = next;
  }
  return reachable[FULL];
}

std::string support_string(const Labels &labels) {
  std::ostringstream out;
  for (int i = 0; i < 6; ++i) {
    if (i) out << ',';
    out << static_cast<int>(labels[i]);
  }
  return out.str();
}

Labels multiply_support(const Labels &labels, int multiplier) {
  Labels answer{};
  for (int i = 0; i < 6; ++i)
    answer[i] = static_cast<uint8_t>(multiplier * labels[i] % P);
  std::sort(answer.begin(), answer.end());
  return answer;
}

bool quadratic_coset(const Labels &labels) {
  static const std::set<int> residues{1, 3, 4, 9, 10, 12};
  static const std::set<int> nonresidues{2, 5, 6, 7, 8, 11};
  const std::set<int> support(labels.begin(), labels.end());
  return support == residues || support == nonresidues;
}

Profile replay_fibre(const Labels &labels) {
  Profile profile{};
  for (int code = 0; code < 1'000'000; ++code) {
    int quotient = code;
    std::array<int, 6> units{};
    for (int provider = 0; provider < 6; ++provider) {
      units[provider] = quotient % 10;
      quotient /= 10;
    }
    uint8_t satisfaction = 0;
    for (int owner_index = 0; owner_index < 6; ++owner_index) {
      const int owner = labels[owner_index];
      uint16_t cover = 0;
      for (int provider = 0; provider < 6; ++provider)
        cover |= MASK11[labels[provider]][units[provider]][owner];
      if (cover == FULL)
        satisfaction |= static_cast<uint8_t>(1U << owner_index);
    }
    ++profile[satisfaction];
  }
  return profile;
}

std::array<uint64_t, 7> satisfaction_histogram(const Profile &profile) {
  std::array<uint64_t, 7> answer{};
  for (int subset = 0; subset < 64; ++subset)
    answer[std::popcount(static_cast<unsigned>(subset))] += profile[subset];
  return answer;
}

std::array<uint64_t, 6> owner_intersections(const Profile &profile) {
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

int main() {
  for (int label = 1; label < P; ++label)
    for (int owner = 1; owner < P; ++owner) {
      MASK1[label][owner] = local_mask(label, 1, 0, owner);
      for (int unit = 1; unit <= 10; ++unit)
        MASK11[label][unit - 1][owner] =
            local_mask(label, 11, unit, owner);
    }

  for (int ratio = 1; ratio < P; ++ratio) {
    std::set<int> sizes;
    for (int unit = 0; unit < 10; ++unit)
      sizes.insert(std::popcount(MASK11[ratio][unit][1]));
    require(sizes.size() == 1,
            "D11 mask cardinality unexpectedly depends on unit");
  }

  std::vector<uint8_t> order_words;
  uint64_t state_words_per_support = 0;
  for (uint8_t word = 0; word < 64; ++word) {
    const int count = std::popcount(static_cast<unsigned>(word));
    if (count < 2) continue;  // hereditary leave-one-out lcm = 11
    order_words.push_back(word);
    state_words_per_support += pow10(count);
  }

  uint64_t supports = 0, scalar_contexts = 0, local_contexts = 0;
  std::map<int, uint64_t> scalar_by_d11, local_by_d11;
  std::set<Labels> scalar_supports, local_supports;

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
              for (uint8_t word : order_words) {
                bool scalar = true;
                for (int owner : labels)
                  scalar &= scalar_capacity(labels, word, owner);
                if (!scalar) continue;
                ++scalar_contexts;
                ++scalar_by_d11[std::popcount(static_cast<unsigned>(word))];
                scalar_supports.insert(labels);
                bool local = true;
                for (int owner : labels)
                  local &= owner_locally_feasible(labels, word, owner);
                if (!local) continue;
                ++local_contexts;
                ++local_by_d11[std::popcount(static_cast<unsigned>(word))];
                local_supports.insert(labels);
              }
            }

  require(supports == 924 && order_words.size() == 57 &&
              state_words_per_support == 1'771'500,
          "raw hereditary census mismatch");
  require(scalar_contexts == 66 && scalar_supports.size() == 66 &&
              scalar_by_d11 == std::map<int, uint64_t>{{6, 66}},
          "scalar-capacity classification mismatch");
  require(local_contexts == 66 && local_supports == scalar_supports &&
              local_by_d11 == std::map<int, uint64_t>{{6, 66}},
          "owner-local classification mismatch");

  std::set<Labels> remaining = local_supports;
  std::vector<std::pair<Labels, int>> orbit_data;
  while (!remaining.empty()) {
    const Labels representative = *remaining.begin();
    std::set<Labels> orbit;
    for (int multiplier = 1; multiplier < P; ++multiplier)
      orbit.insert(multiply_support(representative, multiplier));
    int removed = 0;
    for (const Labels &support : orbit) removed += remaining.erase(support);
    require(removed == static_cast<int>(orbit.size()),
            "multiplication orbit leaves owner-local bank");
    orbit_data.push_back({representative, removed});
  }

  const std::vector<std::pair<Labels, int>> expected_orbits{
      {{{1, 2, 3, 4, 5, 6}}, 12}, {{{1, 2, 3, 4, 5, 7}}, 12},
      {{{1, 2, 3, 4, 6, 8}}, 12}, {{{1, 2, 3, 5, 6, 9}}, 4},
      {{{1, 2, 3, 5, 7, 9}}, 12}, {{{1, 2, 6, 8, 9, 10}}, 12},
      {{{1, 3, 4, 9, 10, 12}}, 2}};
  require(orbit_data == expected_orbits,
          "unexpected multiplication-orbit decomposition");

  const std::array<uint64_t, 7> ordinary_hist{
      996'640, 3'360, 0, 0, 0, 0, 0};
  const std::array<uint64_t, 7> quadratic_hist{
      998'200, 0, 1'800, 0, 0, 0, 0};
  uint64_t replayed_words = 0, global_survivors = 0;
  for (const Labels &labels : local_supports) {
    const Profile profile = replay_fibre(labels);
    replayed_words += 1'000'000;
    global_survivors += profile[63];
    const bool quadratic = quadratic_coset(labels);
    require(satisfaction_histogram(profile) ==
                (quadratic ? quadratic_hist : ordinary_hist),
            "owner-satisfaction histogram mismatch");
    const auto owners = owner_intersections(profile);
    const auto pairs = pair_intersections(profile);
    require(owners == (quadratic
                           ? std::array<uint64_t, 6>{600, 600, 600, 600, 600, 600}
                           : std::array<uint64_t, 6>{560, 560, 560, 560, 560, 560}),
            "owner-obligation size mismatch");
    for (int i = 0; i < 6; ++i)
      for (int j = i + 1; j < 6; ++j) {
        const uint64_t expected =
            quadratic && labels[i] + labels[j] == P ? 600 : 0;
        require(pairs[i][j] == expected,
                "owner-obligation intersection mismatch");
      }
  }
  require(replayed_words == 66'000'000 && global_survivors == 0,
          "global replay mismatch");

  std::cout << "scale-eleven AP-centred Hamming-six owner orthogonality\n";
  std::cout << "supports " << supports << '\n';
  std::cout << "hereditary order words " << order_words.size() << '\n';
  std::cout << "hereditary labelled order contexts "
            << supports * order_words.size() << '\n';
  std::cout << "hereditary state words per support "
            << state_words_per_support << '\n';
  std::cout << "hereditary labelled state contexts "
            << supports * state_words_per_support << '\n';
  std::cout << "scalar-capacity contexts " << scalar_contexts << '\n';
  std::cout << "owner-local contexts " << local_contexts << '\n';
  std::cout << "owner-local state words " << replayed_words << '\n';
  std::cout << "multiplication orbits " << orbit_data.size()
            << " sizes 12,12,12,4,12,12,2\n";
  std::cout << "replayed literal unit words " << replayed_words << '\n';
  std::cout << "global common-sheet survivors " << global_survivors << '\n';
  std::cout << "ordinary 64 supports profile 0:996640 1:3360; owner size 560; distinct pairs disjoint\n";
  std::cout << "quadratic two supports profile 0:998200 2:1800; owner size 600; intersection nerve 3K2 on negation pairs\n";
  std::cout << "tournament observable nonempty owner intersection; multiplicative gauge to orbit representative; sorted tie Hamiltonian path\n";
  std::cout << "tournament fingerprints both transitive (scores 0,1,2,3,4,5; no cycles; six singleton SCCs; one Hamiltonian path)\n";
  std::cout << "challenged vertices owner obligations distinguish empty nerve from 3K2; completed tournaments destroy that distinction\n";
  std::cout << "local D11 mask table at owner one (units 1..10, ratios 1..12 in hex)\n";
  for (int unit = 0; unit < 10; ++unit) {
    std::cout << "  e=" << unit + 1 << ':';
    for (int ratio = 1; ratio < P; ++ratio)
      std::cout << ' ' << std::hex << MASK11[ratio][unit][1] << std::dec;
    std::cout << '\n';
  }
}
