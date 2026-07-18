#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

// Independent standard-library referee for THM-1124.  This program starts
// from bounded literal CRT search, enumerates every labelled support/order
// context, and independently checks the proof-facing Z/3 upper relaxation.
// It deliberately does not share code, NumPy tables, digests, or assertions
// with the Python primary.
//
// Reproduce from the repository root with:
//
//   clang++ -std=c++20 -O3 -Wall -Wextra -pedantic \
//     04-computation/lrc13_scale_thirty_three_hamming_six_three_fibre_referee_codex_c33.cpp \
//     -o /tmp/thm1124-c33
//   /tmp/thm1124-c33 | cmp - \
//     05-knowledge/results/lrc13_scale_thirty_three_hamming_six_three_fibre_referee_codex_c33.out

namespace {

constexpr int P = 13;
constexpr int C = 33;
constexpr int W = 6;
constexpr std::array<int, 4> ORDERS{1, 3, 11, 33};
constexpr uint64_t FULL = (uint64_t{1} << C) - 1;

using Word = std::array<uint8_t, W>;
using Support = std::array<uint8_t, W>;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

int mod(int value, int modulus) {
  int residue = value % modulus;
  return residue < 0 ? residue + modulus : residue;
}

int inverse_mod_13(int value) {
  for (int candidate = 1; candidate < P; ++candidate)
    if (value * candidate % P == 1) return candidate;
  fail("nonunit requested modulo thirteen");
}

int centred(int value, int modulus) {
  int residue = mod(value, modulus);
  return 2 * residue > modulus ? residue - modulus : residue;
}

std::vector<uint8_t> units(int order) {
  if (order == 1) return {0};
  std::vector<uint8_t> answer;
  for (int value = 1; value < order; ++value)
    if (std::gcd(value, order) == 1)
      answer.push_back(static_cast<uint8_t>(value));
  return answer;
}

int bounded_crt(int label, int order, int unit) {
  int answer = -1;
  for (int candidate = 0; candidate < P * order; ++candidate) {
    if (candidate % P != order * label % P) continue;
    if (candidate % order != unit % order) continue;
    require(answer < 0, "CRT search returned two representatives");
    answer = candidate;
  }
  require(answer >= 0, "CRT search returned no representative");
  return answer;
}

uint64_t literal_mask(int label, int order, int unit, int owner) {
  const int base = bounded_crt(label, order, unit);
  const int owner_inverse = inverse_mod_13(owner);
  uint64_t mask = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centred(
        base * (owner_inverse + P * sheet), P * order);
    if (-order < value && value <= order)
      mask |= uint64_t{1} << sheet;
  }
  return mask;
}

bool has_period(uint64_t mask, int period) {
  for (int sheet = 0; sheet < C; ++sheet) {
    const int translated = (sheet + period) % C;
    if (((mask >> sheet) & 1U) != ((mask >> translated) & 1U))
      return false;
  }
  return true;
}

struct Tables {
  std::array<std::vector<uint8_t>, 4> unit_bank;
  std::array<std::array<std::array<std::vector<uint64_t>, P>, 4>, P> masks;
  std::array<std::array<std::array<uint8_t, P>, 4>, P> cards{};
};

Tables build_tables() {
  Tables tables;
  const std::array<std::size_t, 4> expected{1, 2, 10, 20};
  for (int oi = 0; oi < 4; ++oi) {
    tables.unit_bank[oi] = units(ORDERS[oi]);
    require(tables.unit_bank[oi].size() == expected[oi],
            "Euler-unit count mismatch");
  }

  for (int label = 1; label < P; ++label) {
    for (int oi = 0; oi < 4; ++oi) {
      const int order = ORDERS[oi];
      for (int owner = 1; owner < P; ++owner) {
        std::vector<uint64_t> choices;
        int common_cardinality = -1;
        for (uint8_t unit : tables.unit_bank[oi]) {
          const uint64_t mask = literal_mask(label, order, unit, owner);
          const int cardinality = std::popcount(mask);
          if (common_cardinality < 0) common_cardinality = cardinality;
          require(cardinality == common_cardinality,
                  "literal cardinality depends on unit");
          require(has_period(mask, order),
                  "literal mask violates its effective-order period");
          if (order == 1 || order == 3)
            require(has_period(mask, 3),
                    "putative Z/3 anchor is not a thick-fibre union");
          choices.push_back(mask);
        }
        std::sort(choices.begin(), choices.end());
        choices.erase(std::unique(choices.begin(), choices.end()), choices.end());
        tables.masks[label][oi][owner] = choices;
        tables.cards[label][oi][owner] =
            static_cast<uint8_t>(common_cardinality);

        const int ratio = label * inverse_mod_13(owner) % P;
        const int target = order * ratio % P;
        int interval_hits = 0;
        for (int value = -order + 1; value <= order; ++value)
          interval_hits += static_cast<int>(mod(value, P) == target);
        require(common_cardinality == (C / order) * interval_hits,
                "literal mask and interval count disagree");
      }
    }
  }
  for (int label = 1; label < P; ++label)
    for (int oi = 0; oi < 4; ++oi)
      for (int owner = 1; owner < P; ++owner) {
        const int ratio = label * inverse_mod_13(owner) % P;
        require(tables.cards[label][oi][owner]
                    == tables.cards[ratio][oi][1],
                "ratio covariance failed");
      }
  return tables;
}

bool hereditary(const Word &word) {
  bool literal = true;
  for (int omitted = 0; omitted < W; ++omitted) {
    int value = 1;
    for (int coordinate = 0; coordinate < W; ++coordinate)
      if (coordinate != omitted)
        value = std::lcm(value, ORDERS[word[coordinate]]);
    literal &= value == C;
  }
  int carriers3 = 0;
  int carriers11 = 0;
  for (uint8_t oi : word) {
    carriers3 += ORDERS[oi] % 3 == 0;
    carriers11 += ORDERS[oi] % 11 == 0;
  }
  const bool carrier_form = carriers3 >= 2 && carriers11 >= 2;
  require(literal == carrier_form, "hereditary/carrier grammar mismatch");
  return literal;
}

std::vector<Word> enumerate_words(const Tables &tables, uint64_t &weighted) {
  std::vector<Word> words;
  weighted = 0;
  for (int code = 0; code < 4 * 4 * 4 * 4 * 4 * 4; ++code) {
    int remaining = code;
    Word word{};
    uint64_t unit_weight = 1;
    for (int coordinate = 0; coordinate < W; ++coordinate) {
      word[coordinate] = static_cast<uint8_t>(remaining % 4);
      remaining /= 4;
      unit_weight *= tables.unit_bank[word[coordinate]].size();
    }
    if (!hereditary(word)) continue;
    words.push_back(word);
    weighted += unit_weight;
  }
  require(words.size() == 3249, "hereditary word count mismatch");
  require(weighted == 1'268'394'000ULL,
          "weighted hereditary grammar mismatch");
  return words;
}

int z3_upper_bound(const Tables &tables, const Support &support,
                   const Word &word, int owner_index, int &bank_size) {
  const int owner = support[owner_index];
  std::set<uint64_t> bank{0};
  for (int provider = 0; provider < W; ++provider) {
    const int order = ORDERS[word[provider]];
    if (order != 1 && order != 3) continue;
    std::set<uint64_t> next;
    for (uint64_t partial : bank)
      for (uint64_t option : tables.masks[support[provider]][word[provider]][owner])
        next.insert(partial | option);
    bank = std::move(next);
  }
  bank_size = static_cast<int>(bank.size());
  int best = 0;
  for (uint64_t anchor_union : bank) {
    const uint64_t outside = FULL ^ anchor_union;
    int bound = std::popcount(anchor_union);
    for (int provider = 0; provider < W; ++provider) {
      const int order = ORDERS[word[provider]];
      if (order == 1 || order == 3) continue;
      int contribution = 0;
      for (uint64_t option :
           tables.masks[support[provider]][word[provider]][owner])
        contribution = std::max(
            contribution, std::popcount(option & outside));
      bound += contribution;
    }
    best = std::max(best, bound);
  }
  return best;
}

struct Census {
  std::array<uint64_t, 7> scalar_histogram{};
  std::map<int, int> support_histogram;
  std::map<int, int> bound_histogram;
  std::map<int, int> bank_histogram;
  std::map<int, int> live_histogram;
  std::map<std::pair<int, int>, int> own_order_bound;
  uint64_t contexts = 0;
  uint64_t survivors = 0;
  uint64_t literal_survivors = 0;
  int survivor_supports = 0;
  int carrier_ceiling = 0;
};

Census enumerate_contexts(const Tables &tables, const std::vector<Word> &words) {
  Census census;
  Support support{};
  for (int a = 1; a <= 7; ++a)
  for (int b = a + 1; b <= 8; ++b)
  for (int c = b + 1; c <= 9; ++c)
  for (int d = c + 1; d <= 10; ++d)
  for (int e = d + 1; e <= 11; ++e)
  for (int f = e + 1; f <= 12; ++f) {
    support = {static_cast<uint8_t>(a), static_cast<uint8_t>(b),
               static_cast<uint8_t>(c), static_cast<uint8_t>(d),
               static_cast<uint8_t>(e), static_cast<uint8_t>(f)};
    int support_survivors = 0;
    for (const Word &word : words) {
      ++census.contexts;
      int feasible_owners = 0;
      for (int owner_index = 0; owner_index < W; ++owner_index) {
        int capacity = 0;
        for (int provider = 0; provider < W; ++provider)
          capacity += tables.cards[support[provider]][word[provider]]
                                  [support[owner_index]];
        feasible_owners += capacity >= C;
      }
      ++census.scalar_histogram[feasible_owners];
      if (feasible_owners != W) continue;

      ++census.survivors;
      ++support_survivors;
      uint64_t unit_weight = 1;
      int carrier_owners = 0;
      int live_owners = 0;
      for (uint8_t oi : word)
        unit_weight *= tables.unit_bank[oi].size();
      census.literal_survivors += unit_weight;
      for (int owner_index = 0; owner_index < W; ++owner_index) {
        int bank_size = 0;
        const int bound = z3_upper_bound(
            tables, support, word, owner_index, bank_size);
        ++census.bound_histogram[bound];
        ++census.bank_histogram[bank_size];
        ++census.own_order_bound[{ORDERS[word[owner_index]], bound}];
        live_owners += bound >= C;
        if (ORDERS[word[owner_index]] % 11 == 0) {
          ++carrier_owners;
          census.carrier_ceiling = std::max(census.carrier_ceiling, bound);
          require(bound < C, "eleven-carrier owner survived Z/3 relaxation");
        }
      }
      require(carrier_owners >= 2,
              "hereditary survivor lost the eleven-carrier bridge");
      ++census.live_histogram[live_owners];
    }
    ++census.support_histogram[support_survivors];
    census.survivor_supports += support_survivors != 0;
  }
  return census;
}

template <class Map>
void print_map(const Map &values) {
  bool first = true;
  for (const auto &[key, count] : values) {
    if (!first) std::cout << ' ';
    std::cout << key << ':' << count;
    first = false;
  }
  std::cout << '\n';
}

}  // namespace

int main() {
  const Tables tables = build_tables();
  uint64_t weighted = 0;
  const std::vector<Word> words = enumerate_words(tables, weighted);
  const Census census = enumerate_contexts(tables, words);

  require(census.contexts == 3'002'076ULL, "labelled context count mismatch");
  require(census.scalar_histogram
              == std::array<uint64_t, 7>{16'086, 446'568, 1'446'252,
                  891'684, 171'594, 28'548, 1'344},
          "scalar histogram mismatch");
  require(census.support_histogram
              == std::map<int, int>{{0, 774}, {2, 48}, {3, 24}, {4, 24},
                  {9, 24}, {18, 12}, {36, 18}},
          "support histogram mismatch");
  require(census.survivors == 1'344, "scalar survivor count mismatch");
  require(census.survivor_supports == 150,
          "scalar survivor support count mismatch");
  require(census.literal_survivors == 40'022'400ULL,
          "literal survivor weight mismatch");
  require(census.bank_histogram
              == std::map<int, int>{{1, 444}, {2, 3'336}, {3, 4'284}},
          "Z/3 bank histogram mismatch");
  require(census.bound_histogram
              == std::map<int, int>{{25, 996}, {26, 4'248}, {27, 384},
                  {28, 288}, {29, 228}, {30, 96}, {33, 1'824}},
          "Z/3 bound histogram mismatch");
  require(census.live_histogram
              == std::map<int, int>{{0, 306}, {1, 360}, {2, 624}, {4, 54}},
          "Z/3 live-owner histogram mismatch");
  require(census.carrier_ceiling == 29,
          "eleven-carrier ceiling mismatch");

  std::cout << "scale=33 p=13 hamming=6 referee=literal-crt-Z3-upper-relaxation\n";
  std::cout << "orders=(1,3,11,33) unit-counts=(1,2,10,20)\n";
  std::cout << "hereditary-grammar=at-least-two-3-carriers-and-two-11-carriers\n";
  std::cout << "hereditary-words=" << words.size()
            << " weighted-states/support=" << weighted
            << " labelled-support-order-contexts=" << census.contexts
            << " raw-labelled-states=" << 924ULL * weighted << '\n';
  std::cout << "scalar-feasible-owners-per-context=";
  for (int index = 0; index <= W; ++index) {
    if (index) std::cout << ' ';
    std::cout << index << ':' << census.scalar_histogram[index];
  }
  std::cout << '\n';
  std::cout << "scalar-supports-by-survivor-count=";
  print_map(census.support_histogram);
  std::cout << "scalar-survivors=" << census.survivors
            << " supports=" << census.survivor_supports
            << " literal-unit-words=" << census.literal_survivors << '\n';
  std::cout << "Z3-anchor-bank-size=";
  print_map(census.bank_histogram);
  std::cout << "Z3-upper-bound=";
  print_map(census.bound_histogram);
  std::cout << "Z3-live-owners-per-context=";
  print_map(census.live_histogram);
  std::cout << "11-carrier-owner-Z3-ceiling=" << census.carrier_ceiling
            << "<33\n";
  std::cout << "literal-cover-implies-Z3-cover=checked-on-all-8064-owner-obligations; "
               "hereditary-supplies-at-least-two-terminal-11-carrier-owners\n";
  std::cout << "proof-carrier=owner-coloured-Z3-thick-fibre-incidence; "
               "destroyed=shared-units-and-nonanchor-overlaps\n";
  std::cout << "VERDICT=primitive-proper-AP-centred-common-scale-33-H6-face-EMPTY\n";
}
