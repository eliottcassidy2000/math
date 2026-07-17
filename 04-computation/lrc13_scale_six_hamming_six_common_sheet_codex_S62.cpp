#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

using Word = std::array<uint8_t, 6>;

static constexpr int P = 13;
static constexpr int C = 6;
static constexpr std::array<int, 6> ORDERS{1, 2, 3, 3, 6, 6};
static constexpr std::array<int, 6> UNITS{0, 1, 1, 2, 1, 5};

int centered(int value, int modulus) {
  int residue = value % modulus;
  if (residue < 0) residue += modulus;
  return 2 * residue > modulus ? residue - modulus : residue;
}

int inverse_mod_13(int value) {
  for (int candidate = 1; candidate < P; ++candidate)
    if ((value * candidate) % P == 1) return candidate;
  std::abort();
}

int crt_base(int label, int state) {
  int order = ORDERS[state], unit = UNITS[state];
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P && value % order == unit % order)
      return value;
  std::abort();
}

uint64_t local_mask(int label, int state, int owner) {
  int order = ORDERS[state], base = crt_base(label, state);
  int inverse = inverse_mod_13(owner);
  uint64_t answer = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    int value = centered(base * (inverse + P * sheet), P * order);
    if (-order < value && value <= order) answer |= uint64_t{1} << sheet;
  }
  return answer;
}

std::string order_pattern(const Word &word) {
  std::array<int, 6> values{};
  for (int i = 0; i < 6; ++i) values[i] = ORDERS[word[i]];
  std::sort(values.begin(), values.end());
  std::string answer;
  for (int value : values) {
    if (!answer.empty()) answer += ',';
    answer += std::to_string(value);
  }
  return answer;
}

int main() {
  std::vector<Word> words;
  for (int code = 0; code < 46656; ++code) {
    int work = code, even_count = 0, three_count = 0;
    Word word{};
    for (int i = 0; i < 6; ++i) {
      word[i] = work % 6;
      work /= 6;
      even_count += ORDERS[word[i]] % 2 == 0;
      three_count += ORDERS[word[i]] % 3 == 0;
    }
    if (even_count >= 2 && three_count >= 2) words.push_back(word);
  }

  std::map<std::string, uint64_t> tested_by_pattern;
  std::map<std::string, uint64_t> survived_by_pattern;
  std::map<std::string, uint64_t> order_tested_by_pattern;
  std::map<std::string, uint64_t> capacity_by_pattern;
  std::map<std::string, uint64_t> local_feasible_by_pattern;
  std::map<std::array<int, 6>, uint64_t> survived_by_labels;
  uint64_t tested = 0, survived = 0;
  uint64_t order_tested = 0, capacity_survived = 0, local_feasible = 0;
  std::array<int, 6> labels{};
  for (labels[0] = 1; labels[0] <= 7; ++labels[0])
  for (labels[1] = labels[0] + 1; labels[1] <= 8; ++labels[1])
  for (labels[2] = labels[1] + 1; labels[2] <= 9; ++labels[2])
  for (labels[3] = labels[2] + 1; labels[3] <= 10; ++labels[3])
  for (labels[4] = labels[3] + 1; labels[4] <= 11; ++labels[4])
  for (labels[5] = labels[4] + 1; labels[5] <= 12; ++labels[5]) {
    std::array<std::array<uint64_t, 6>, 6> packed{};
    for (int provider = 0; provider < 6; ++provider)
      for (int state = 0; state < 6; ++state)
        for (int owner = 0; owner < 6; ++owner)
          packed[provider][state] |=
              local_mask(labels[provider], state, labels[owner]) << (C * owner);
    constexpr uint64_t full = (uint64_t{1} << 36) - 1;

    // Unit-free order layer.  A word passes `capacity` if the sums of the
    // literal mask cardinalities can reach six at each owner.  It passes
    // `local feasible` if every owner can be covered by some unit assignment,
    // allowing a different assignment at each owner.  Both are necessary
    // relaxations of the global common-sheet CSP.
    for (int order_code = 0; order_code < 4096; ++order_code) {
      int work = order_code, even_count = 0, three_count = 0;
      std::array<int, 6> order_word{};
      Word representatives{};
      for (int i = 0; i < 6; ++i) {
        order_word[i] = work % 4;
        work /= 4;
        int order = std::array<int, 4>{1, 2, 3, 6}[order_word[i]];
        even_count += order % 2 == 0;
        three_count += order % 3 == 0;
        representatives[i] = std::array<int, 4>{0, 1, 2, 4}[order_word[i]];
      }
      if (even_count < 2 || three_count < 2) continue;
      std::string pattern = order_pattern(representatives);
      ++order_tested;
      ++order_tested_by_pattern[pattern];
      bool capacity_ok = true;
      for (int owner = 0; owner < 6; ++owner) {
        int cardinality = 0;
        for (int provider = 0; provider < 6; ++provider) {
          uint64_t mask = local_mask(labels[provider], representatives[provider], labels[owner]);
          cardinality += __builtin_popcountll(mask);
        }
        capacity_ok &= cardinality >= C;
      }
      if (!capacity_ok) continue;
      ++capacity_survived;
      ++capacity_by_pattern[pattern];

      bool all_owners_locally_feasible = true;
      for (int owner = 0; owner < 6 && all_owners_locally_feasible; ++owner) {
        bool owner_feasible = false;
        for (int unit_bits = 0; unit_bits < 64 && !owner_feasible; ++unit_bits) {
          uint64_t cover = 0;
          for (int provider = 0; provider < 6; ++provider) {
            int state = representatives[provider];
            if (order_word[provider] >= 2 && ((unit_bits >> provider) & 1)) ++state;
            cover |= local_mask(labels[provider], state, labels[owner]);
          }
          owner_feasible = cover == ((uint64_t{1} << C) - 1);
        }
        all_owners_locally_feasible &= owner_feasible;
      }
      if (!all_owners_locally_feasible) continue;
      ++local_feasible;
      ++local_feasible_by_pattern[pattern];
    }

    for (const Word &word : words) {
      auto pattern = order_pattern(word);
      ++tested_by_pattern[pattern];
      ++tested;
      uint64_t cover = 0;
      for (int provider = 0; provider < 6; ++provider)
        cover |= packed[provider][word[provider]];
      if (cover != full) continue;
      ++survived;
      ++survived_by_pattern[pattern];
      ++survived_by_labels[labels];
    }
  }

  std::cout << "THM-960 SCALE-SIX HAMMING-SIX COMMON-SHEET EXHAUSTION\n";
  std::cout << "admissible state words " << words.size() << '\n';
  std::cout << "admissible order words 3249\n";
  std::cout << "order contexts " << order_tested << '\n';
  std::cout << "capacity order contexts " << capacity_survived << '\n';
  std::cout << "owner-locally-feasible order contexts " << local_feasible << '\n';
  std::cout << "tested " << tested << '\n';
  std::cout << "survived " << survived << '\n';
  std::cout << "surviving label sets " << survived_by_labels.size() << '\n';
  std::cout << "patterns\n";
  for (const auto &[pattern, count] : tested_by_pattern)
    std::cout << pattern << " tested=" << count
              << " survived=" << survived_by_pattern[pattern]
              << " order_contexts=" << order_tested_by_pattern[pattern]
              << " capacity=" << capacity_by_pattern[pattern]
              << " owner_local=" << local_feasible_by_pattern[pattern] << '\n';
  std::cout << "label multiplicity histogram\n";
  std::map<uint64_t, int> histogram;
  for (const auto &[labels_key, count] : survived_by_labels) ++histogram[count];
  for (const auto &[count, frequency] : histogram)
    std::cout << count << ':' << frequency << '\n';
}
