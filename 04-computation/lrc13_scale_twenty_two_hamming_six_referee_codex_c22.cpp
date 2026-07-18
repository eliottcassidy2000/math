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
#include <tuple>
#include <vector>

// Independent exact referee for the primitive proper AP-centred common-scale-
// twenty-two Hamming-six face.  In contrast to the Python primary, this replay
// finds every CRT representative by literal search, represents reachable
// owner-local unions as sorted integer vectors, and derives the hereditary
// grammar directly from the two prime-power coverage conditions.

namespace {

constexpr int P = 13;
constexpr int C = 22;
constexpr uint32_t FULL = (1U << C) - 1U;
constexpr std::array<uint8_t, 4> ORDERS{1, 2, 11, 22};
constexpr std::array<std::array<uint8_t, 10>, 4> UNITS{{
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}},
    {{1, 3, 5, 7, 9, 13, 15, 17, 19, 21}},
}};
constexpr std::array<int, 4> UNIT_COUNTS{1, 1, 10, 10};

using Support = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;  // actual orders, not indices
using Capacities = std::array<uint8_t, 6>;
using Multiplicity = std::array<uint8_t, 4>;
using MaskFibre = std::array<uint32_t, 10>;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

// Compact self-contained SHA-256.  All certificate serializations below are
// explicit byte streams; they do not depend on C++ object layout or hashing.
class Sha256 {
 public:
  void byte(uint8_t value) {
    buffer_[used_++] = value;
    bit_count_ += 8;
    if (used_ == 64) {
      compress(buffer_.data());
      used_ = 0;
    }
  }

  void u32le(uint32_t value) {
    for (int shift = 0; shift < 32; shift += 8)
      byte(static_cast<uint8_t>(value >> shift));
  }

  std::string finish() {
    const uint64_t original_bits = bit_count_;
    byte(0x80);
    while (used_ != 56) byte(0);
    // byte() counts the padding, so serialize the saved pre-padding length.
    for (int shift = 56; shift >= 0; shift -= 8)
      byte(static_cast<uint8_t>(original_bits >> shift));
    require(used_ == 0, "SHA-256 final block did not close");
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (uint32_t word : state_) out << std::setw(8) << word;
    return out.str();
  }

 private:
  static constexpr std::array<uint32_t, 64> K{
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U,
      0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
      0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
      0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
      0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
      0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
      0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
      0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
      0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
      0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
      0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U,
      0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
      0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U,
      0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
      0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
      0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};

  static uint32_t load_be(const uint8_t *p) {
    return (static_cast<uint32_t>(p[0]) << 24) |
           (static_cast<uint32_t>(p[1]) << 16) |
           (static_cast<uint32_t>(p[2]) << 8) |
           static_cast<uint32_t>(p[3]);
  }

  static uint32_t rotate(uint32_t x, int count) {
    return (x >> count) | (x << (32 - count));
  }

  void compress(const uint8_t *block) {
    std::array<uint32_t, 64> words{};
    for (int i = 0; i < 16; ++i) words[i] = load_be(block + 4 * i);
    for (int i = 16; i < 64; ++i) {
      const uint32_t s0 = rotate(words[i - 15], 7) ^
                          rotate(words[i - 15], 18) ^ (words[i - 15] >> 3);
      const uint32_t s1 = rotate(words[i - 2], 17) ^
                          rotate(words[i - 2], 19) ^ (words[i - 2] >> 10);
      words[i] = words[i - 16] + s0 + words[i - 7] + s1;
    }
    uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
    uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
    for (int i = 0; i < 64; ++i) {
      const uint32_t s1 = rotate(e, 6) ^ rotate(e, 11) ^ rotate(e, 25);
      const uint32_t choose = (e & f) ^ (~e & g);
      const uint32_t temp1 = h + s1 + choose + K[i] + words[i];
      const uint32_t s0 = rotate(a, 2) ^ rotate(a, 13) ^ rotate(a, 22);
      const uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const uint32_t temp2 = s0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + temp1;
      d = c;
      c = b;
      b = a;
      a = temp1 + temp2;
    }
    state_[0] += a;
    state_[1] += b;
    state_[2] += c;
    state_[3] += d;
    state_[4] += e;
    state_[5] += f;
    state_[6] += g;
    state_[7] += h;
  }

  std::array<uint32_t, 8> state_{0x6a09e667U, 0xbb67ae85U,
                                 0x3c6ef372U, 0xa54ff53aU,
                                 0x510e527fU, 0x9b05688cU,
                                 0x1f83d9abU, 0x5be0cd19U};
  std::array<uint8_t, 64> buffer_{};
  uint64_t bit_count_ = 0;
  std::size_t used_ = 0;
};

int order_index(int order) {
  for (int i = 0; i < 4; ++i)
    if (ORDERS[i] == order) return i;
  fail("unknown effective order");
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

int literal_crt_base(int label, int order, int unit) {
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P &&
        value % order == unit % order)
      return value;
  fail("literal CRT search found no representative");
}

uint32_t literal_local_mask(int label, int order, int unit, int owner) {
  const int base = literal_crt_base(label, order, unit);
  const int inverse = inverse_mod_13(owner);
  uint32_t mask = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(base * (inverse + P * sheet), P * order);
    if (-order < value && value <= order) mask |= 1U << sheet;
  }
  return mask;
}

struct Tables {
  // [label][order-index][owner][unit-index]
  std::array<std::array<std::array<MaskFibre, P>, 4>, P> masks{};
  std::array<std::array<std::array<uint8_t, P>, 4>, P> cards{};
  std::string mask_digest;
};

Tables build_tables() {
  Tables table;
  Sha256 digest;
  for (int label = 1; label < P; ++label)
    for (int oi = 0; oi < 4; ++oi) {
      const int order = ORDERS[oi];
      for (int owner = 1; owner < P; ++owner) {
        int common_card = -1;
        for (int ui = 0; ui < UNIT_COUNTS[oi]; ++ui) {
          const int unit = UNITS[oi][ui];
          const uint32_t mask =
              literal_local_mask(label, order, unit, owner);
          table.masks[label][oi][owner][ui] = mask;
          const int card = std::popcount(mask);
          if (common_card < 0) common_card = card;
          require(card == common_card,
                  "sheet cardinality depends on the unit choice");
          digest.byte(static_cast<uint8_t>(label));
          digest.byte(static_cast<uint8_t>(order));
          digest.byte(static_cast<uint8_t>(unit));
          digest.byte(static_cast<uint8_t>(owner));
          digest.u32le(mask);
        }
        table.cards[label][oi][owner] =
            static_cast<uint8_t>(common_card);
        const int ratio = label * inverse_mod_13(owner) % P;
        require(common_card ==
                    std::popcount(literal_local_mask(ratio, order,
                                                     UNITS[oi][0], 1)),
                "provider/owner ratio reduction failed");
      }
    }
  table.mask_digest = digest.finish();
  return table;
}

bool hereditary_prime_power(const OrderWord &word) {
  int even = 0;
  int eleven = 0;
  for (int order : word) {
    even += order % 2 == 0;
    eleven += order % 11 == 0;
  }
  return even >= 2 && eleven >= 2;
}

bool hereditary_lcm(const OrderWord &word) {
  for (int omitted = 0; omitted < 6; ++omitted) {
    int residual = 1;
    for (int coordinate = 0; coordinate < 6; ++coordinate)
      if (coordinate != omitted)
        residual = std::lcm(residual, static_cast<int>(word[coordinate]));
    if (residual != C) return false;
  }
  return true;
}

void enumerate_words(std::vector<OrderWord> &words, OrderWord &word,
                     int coordinate, Sha256 &digest) {
  if (coordinate == 6) {
    require(hereditary_prime_power(word) == hereditary_lcm(word),
            "prime-power and leave-one-out-lcm grammars disagree");
    if (hereditary_prime_power(word)) {
      words.push_back(word);
      for (uint8_t order : word) digest.byte(order);
    }
    return;
  }
  for (uint8_t order : ORDERS) {
    word[coordinate] = order;
    enumerate_words(words, word, coordinate + 1, digest);
  }
}

std::vector<Support> all_supports() {
  std::vector<Support> result;
  for (int a = 1; a <= 7; ++a)
    for (int b = a + 1; b <= 8; ++b)
      for (int c = b + 1; c <= 9; ++c)
        for (int d = c + 1; d <= 10; ++d)
          for (int e = d + 1; e <= 11; ++e)
            for (int f = e + 1; f <= 12; ++f)
              result.push_back(Support{static_cast<uint8_t>(a),
                                       static_cast<uint8_t>(b),
                                       static_cast<uint8_t>(c),
                                       static_cast<uint8_t>(d),
                                       static_cast<uint8_t>(e),
                                       static_cast<uint8_t>(f)});
  return result;
}

Capacities capacity_vector(const Support &support, const OrderWord &word,
                           const Tables &table) {
  Capacities result{};
  for (int owner_index = 0; owner_index < 6; ++owner_index)
    for (int provider = 0; provider < 6; ++provider)
      result[owner_index] +=
          table.cards[support[provider]][order_index(word[provider])]
                     [support[owner_index]];
  return result;
}

Multiplicity multiplicity(const OrderWord &word) {
  Multiplicity result{};
  for (int order : word) ++result[order_index(order)];
  return result;
}

struct ScalarRow {
  Support support{};
  OrderWord word{};
  Capacities capacities{};
};

struct LocalRow {
  bool feasible = false;
  uint8_t maximum = 0;
  uint32_t reachable_count = 0;
  std::vector<uint32_t> reachable;
};

LocalRow owner_local(const Support &support, const OrderWord &word, int owner,
                     const Tables &table) {
  std::vector<uint32_t> reachable{0};
  for (int provider = 0; provider < 6; ++provider) {
    const int oi = order_index(word[provider]);
    std::vector<uint32_t> next;
    next.reserve(reachable.size() * UNIT_COUNTS[oi]);
    for (uint32_t partial : reachable)
      for (int ui = 0; ui < UNIT_COUNTS[oi]; ++ui)
        next.push_back(partial |
                       table.masks[support[provider]][oi][owner][ui]);
    std::sort(next.begin(), next.end());
    next.erase(std::unique(next.begin(), next.end()), next.end());
    reachable = std::move(next);
  }
  LocalRow result;
  result.reachable_count = static_cast<uint32_t>(reachable.size());
  result.reachable = std::move(reachable);
  for (uint32_t mask : result.reachable) {
    result.maximum = static_cast<uint8_t>(
        std::max<int>(result.maximum, std::popcount(mask)));
    result.feasible |= mask == FULL;
  }
  require(result.feasible == (result.maximum == C),
          "full-mask feasibility and maximum-union test disagree");
  return result;
}

struct Tournament {
  int ties = 0;
  int flips = 0;
  int triangles = 0;
  int sccs = 0;
  uint64_t paths = 0;
  std::array<uint8_t, 6> sorted_scores{};
};

Tournament tournament(const Capacities &capacities,
                      const std::array<LocalRow, 6> &locals) {
  std::array<uint8_t, 6> out{};
  Tournament result;
  for (int left = 0; left < 6; ++left)
    for (int right = left + 1; right < 6; ++right) {
      const std::array<int, 3> left_key{
          locals[left].feasible, locals[left].maximum, capacities[left]};
      const std::array<int, 3> right_key{
          locals[right].feasible, locals[right].maximum, capacities[right]};
      int winner = left;
      if (left_key == right_key)
        ++result.ties;  // coordinate tie path: the earlier vertex wins
      else if (right_key > left_key) {
        winner = right;
        ++result.flips;
      }
      out[winner] |= static_cast<uint8_t>(1U << (left + right - winner));
    }

  for (int vertex = 0; vertex < 6; ++vertex)
    result.sorted_scores[vertex] =
        static_cast<uint8_t>(std::popcount(out[vertex]));
  std::sort(result.sorted_scores.begin(), result.sorted_scores.end());

  for (int a = 0; a < 6; ++a)
    for (int b = a + 1; b < 6; ++b)
      for (int c = b + 1; c < 6; ++c) {
        const bool forward = ((out[a] >> b) & 1U) &&
                             ((out[b] >> c) & 1U) &&
                             ((out[c] >> a) & 1U);
        const bool reverse = ((out[a] >> c) & 1U) &&
                             ((out[c] >> b) & 1U) &&
                             ((out[b] >> a) & 1U);
        result.triangles += forward || reverse;
      }

  std::array<uint8_t, 6> reach = out;
  for (int v = 0; v < 6; ++v) reach[v] |= static_cast<uint8_t>(1U << v);
  for (int middle = 0; middle < 6; ++middle)
    for (int source = 0; source < 6; ++source)
      if ((reach[source] >> middle) & 1U) reach[source] |= reach[middle];
  uint8_t assigned = 0;
  for (int root = 0; root < 6; ++root)
    if (!((assigned >> root) & 1U)) {
      ++result.sccs;
      for (int v = 0; v < 6; ++v)
        if (((reach[root] >> v) & 1U) && ((reach[v] >> root) & 1U))
          assigned |= static_cast<uint8_t>(1U << v);
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
  for (int last = 0; last < 6; ++last) result.paths += paths.back()[last];
  return result;
}

template <typename Key>
std::string histogram(const std::map<Key, uint64_t> &values) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, count] : values) {
    if (!first) out << ' ';
    first = false;
    out << key << ':' << count;
  }
  return out.str();
}

std::string multiplicity_histogram(
    const std::map<Multiplicity, uint64_t> &values) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, count] : values) {
    if (!first) out << ' ';
    first = false;
    out << static_cast<int>(key[0]) << ',' << static_cast<int>(key[1]) << ','
        << static_cast<int>(key[2]) << ',' << static_cast<int>(key[3]) << ':'
        << count;
  }
  return out.str();
}

Support multiply_support(const Support &support, int multiplier) {
  Support result{};
  for (int i = 0; i < 6; ++i)
    result[i] = static_cast<uint8_t>(support[i] * multiplier % P);
  std::sort(result.begin(), result.end());
  return result;
}

}  // namespace

int main() {
  const Tables table = build_tables();
  require(table.mask_digest ==
              "a77964e49c10fc3731f7948059b315dc9f5d94b98ba611ebf6f1c1f9fa5fb26b",
          "literal CRT mask-table digest differs from the frozen primary");

  std::vector<OrderWord> words;
  OrderWord scratch{};
  Sha256 order_digest_state;
  enumerate_words(words, scratch, 0, order_digest_state);
  const std::string order_digest = order_digest_state.finish();
  require(words.size() == 3249, "hereditary order-word census mismatch");
  require(order_digest ==
              "9c9ea6b5101659f4a3e958bb81bd859b73b05b9f1f04cfbfb65b358352a31f11",
          "hereditary order-word digest differs from the frozen primary");

  uint64_t state_words_per_support = 0;
  for (const OrderWord &word : words) {
    uint64_t fibre = 1;
    for (int order : word) fibre *= UNIT_COUNTS[order_index(order)];
    state_words_per_support += fibre;
  }
  require(state_words_per_support == 100'975'500ULL,
          "literal state-word count mismatch");

  const std::vector<Support> supports = all_supports();
  require(supports.size() == 924, "support census mismatch");
  std::vector<ScalarRow> bank;
  std::set<Support> scalar_supports;
  std::map<int, uint64_t> contexts_per_support;
  std::map<Multiplicity, uint64_t> multiplicities;
  std::set<Capacities> distinct_capacities;
  std::map<int, uint64_t> minimum_slack;
  std::map<int, uint64_t> maximum_slack;
  std::map<int, uint64_t> tight_owners;
  Sha256 scalar_digest_state;

  for (const Support &support : supports) {
    int support_contexts = 0;
    for (const OrderWord &word : words) {
      const Capacities capacities = capacity_vector(support, word, table);
      if (*std::min_element(capacities.begin(), capacities.end()) < C)
        continue;
      bank.push_back(ScalarRow{support, word, capacities});
      scalar_supports.insert(support);
      ++support_contexts;
      ++multiplicities[multiplicity(word)];
      distinct_capacities.insert(capacities);
      ++minimum_slack[*std::min_element(capacities.begin(), capacities.end()) -
                      C];
      ++maximum_slack[*std::max_element(capacities.begin(), capacities.end()) -
                      C];
      ++tight_owners[std::count(capacities.begin(), capacities.end(), C)];
      for (uint8_t value : support) scalar_digest_state.byte(value);
      for (uint8_t value : word) scalar_digest_state.byte(value);
      for (uint8_t value : capacities) scalar_digest_state.byte(value);
    }
    if (support_contexts) ++contexts_per_support[support_contexts];
  }
  const std::string scalar_digest = scalar_digest_state.finish();
  require(bank.size() == 984 && scalar_supports.size() == 180,
          "scalar survivor census mismatch");
  require(scalar_digest ==
              "fb618d8e443ddfa5f118dbcff16c5e196d8693240ede4844e8266aa6b16980a1",
          "scalar-bank digest differs from the frozen primary");

  const std::map<Multiplicity, uint64_t> expected_multiplicities{
      {{0, 2, 0, 4}, 36},  {{0, 2, 1, 3}, 144},
      {{0, 2, 2, 2}, 216}, {{0, 2, 3, 1}, 144},
      {{0, 2, 4, 0}, 36},  {{0, 3, 1, 2}, 288},
      {{0, 3, 2, 1}, 96},  {{0, 3, 3, 0}, 24}};
  require(multiplicities == expected_multiplicities,
          "scalar multiplicity ledger mismatch");
  require(contexts_per_support ==
              std::map<int, uint64_t>{{2, 96}, {3, 24}, {6, 24}, {16, 36}},
          "contexts-per-support histogram mismatch");

  std::map<int, uint64_t> feasible_contexts;
  std::map<int, uint64_t> maximum_union;
  std::map<int, uint64_t> minimum_owner_maximum;
  std::map<int, uint64_t> reachable_count;
  std::set<std::array<uint8_t, 6>> owner_vectors;
  std::map<int, uint64_t> tie_histogram;
  std::map<int, uint64_t> flip_histogram;
  uint64_t feasible_rows = 0;
  Sha256 reachable_digest_state;

  for (const ScalarRow &row : bank) {
    std::array<LocalRow, 6> locals;
    std::array<uint8_t, 6> owner_vector{};
    int feasible_count = 0;
    for (int owner_index = 0; owner_index < 6; ++owner_index) {
      const int owner = row.support[owner_index];
      locals[owner_index] = owner_local(row.support, row.word, owner, table);
      feasible_count += locals[owner_index].feasible;
      feasible_rows += locals[owner_index].feasible;
      ++maximum_union[locals[owner_index].maximum];
      ++reachable_count[locals[owner_index].reachable_count];
      owner_vector[owner_index] = locals[owner_index].maximum;
      for (uint8_t value : row.support) reachable_digest_state.byte(value);
      for (uint8_t value : row.word) reachable_digest_state.byte(value);
      reachable_digest_state.byte(static_cast<uint8_t>(owner));
      for (uint32_t mask : locals[owner_index].reachable)
        reachable_digest_state.u32le(mask);
    }
    ++feasible_contexts[feasible_count];
    ++minimum_owner_maximum
         [*std::min_element(owner_vector.begin(), owner_vector.end())];
    owner_vectors.insert(owner_vector);

    const Tournament audit = tournament(row.capacities, locals);
    ++tie_histogram[audit.ties];
    ++flip_histogram[audit.flips];
    require(audit.sorted_scores ==
                std::array<uint8_t, 6>{0, 1, 2, 3, 4, 5} &&
                audit.triangles == 0 && audit.sccs == 6 && audit.paths == 1,
            "owner-obligation tournament fingerprint mismatch");
  }
  const std::string reachable_digest = reachable_digest_state.finish();
  require(reachable_digest ==
              "baf8aa9ee67d7686b25e8665bec8f94514d7abb6e7be780bebc2f98039675f1b",
          "reachable-union-bank digest differs from the frozen primary");
  require(feasible_contexts ==
              std::map<int, uint64_t>{{0, 792}, {1, 192}},
          "feasible-owner deficit mismatch");
  require(maximum_union == std::map<int, uint64_t>{{16, 864},
                                                   {17, 1584},
                                                   {18, 2784},
                                                   {19, 480},
                                                   {22, 192}},
          "maximum reachable-union histogram mismatch");

  std::set<Support> remaining = scalar_supports;
  std::map<int, uint64_t> orbit_sizes;
  while (!remaining.empty()) {
    const Support seed = *remaining.begin();
    std::set<Support> orbit;
    for (int multiplier = 1; multiplier < P; ++multiplier)
      orbit.insert(multiply_support(seed, multiplier));
    require(std::includes(scalar_supports.begin(), scalar_supports.end(),
                          orbit.begin(), orbit.end()),
            "scalar support bank is not multiplication-invariant");
    for (const Support &support : orbit) remaining.erase(support);
    ++orbit_sizes[orbit.size()];
  }

  std::cout << "scale-twenty-two H6 literal-CRT sorted-union referee\n";
  std::cout << "hereditary grammar: at least two even orders and at least "
               "two eleven-divisible orders; leave-one-out lcm audited\n";
  std::cout << "supports " << supports.size() << "; hereditary words "
            << words.size() << "; labelled order contexts "
            << supports.size() * words.size() << '\n';
  std::cout << "literal state words/support " << state_words_per_support
            << "; raw labelled states "
            << supports.size() * state_words_per_support << '\n';
  std::cout << "mask SHA256 " << table.mask_digest << '\n';
  std::cout << "order SHA256 " << order_digest << '\n';
  std::cout << "scalar contexts " << bank.size() << " on "
            << scalar_supports.size() << " supports; scalar-bank SHA256 "
            << scalar_digest << '\n';
  std::cout << "scalar multiplicities n1,n2,n11,n22 "
            << multiplicity_histogram(multiplicities) << '\n';
  std::cout << "contexts-per-support histogram "
            << histogram(contexts_per_support) << '\n';
  std::cout << "multiplication orbit-size histogram "
            << histogram(orbit_sizes) << " (telemetry; no quotient)\n";
  std::cout << "capacity vectors " << distinct_capacities.size() << '\n';
  std::cout << "minimum scalar-slack histogram " << histogram(minimum_slack)
            << '\n';
  std::cout << "maximum scalar-slack histogram " << histogram(maximum_slack)
            << '\n';
  std::cout << "tight-owner/context histogram " << histogram(tight_owners)
            << '\n';
  std::cout << "owner-local rows " << 6 * bank.size() << "; feasible rows "
            << feasible_rows << '\n';
  std::cout << "feasible-owner/context histogram "
            << histogram(feasible_contexts) << '\n';
  std::cout << "maximum reachable sheet-union histogram "
            << histogram(maximum_union) << '\n';
  std::cout << "minimum owner maximum/context histogram "
            << histogram(minimum_owner_maximum) << '\n';
  std::cout << "distinct owner max-union vectors " << owner_vectors.size()
            << '\n';
  std::cout << "reachable-union-bank SHA256 " << reachable_digest << '\n';
  std::cout << "reachable-count histogram " << histogram(reachable_count)
            << '\n';
  std::cout << "owner-local all-six contexts " << feasible_contexts[6]
            << '\n';
  std::cout << "tournament carrier owner obligations; observable "
               "(feasible,max-union,capacity); lexicographic switch; "
               "coordinate tie path\n";
  std::cout << "tournament fingerprints all 984 transitive: scores "
               "0,1,2,3,4,5; directed triangles 0; SCCs 6; Hamiltonian "
               "paths 1\n";
  std::cout << "tournament tie-edge histogram " << histogram(tie_histogram)
            << '\n';
  std::cout << "tournament edge-flip histogram " << histogram(flip_histogram)
            << '\n';
  std::cout << "alternate-carrier audit: owner obligations preserve the "
               "terminal projection deficit; runner/provider, divisor, "
               "residue, sheet, and wall-event vertices each destroy "
               "shared-unit incidence, while the tournament itself forgets "
               "the absolute 22-sheet threshold and empty-owner count\n";
  std::cout << "verdict: scalar gate survives, but every context has at most "
               "one feasible owner; the scale-22 common H6 face is empty\n";
}
