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
#include <utility>
#include <vector>

// Independent exact referee for the primitive proper AP-centred common-scale-
// twenty-four Hamming-six face.  Unlike the NumPy/Python primary, this replay
// finds every CRT representative by literal congruence search, enumerates the
// hereditary grammar from prime-power coverage and leave-one-out lcms, and
// represents every owner-local reachable bank as a sorted vector of masks.
// Forward and reverse provider traversals must give the same vector exactly.

namespace {

constexpr int P = 13;
constexpr int C = 24;
constexpr uint32_t FULL = (1U << C) - 1U;
constexpr std::array<uint8_t, 8> ORDERS{1, 2, 3, 4, 6, 8, 12, 24};
constexpr std::array<std::array<uint8_t, 8>, 8> UNITS{{
    {{0, 0, 0, 0, 0, 0, 0, 0}},
    {{1, 0, 0, 0, 0, 0, 0, 0}},
    {{1, 2, 0, 0, 0, 0, 0, 0}},
    {{1, 3, 0, 0, 0, 0, 0, 0}},
    {{1, 5, 0, 0, 0, 0, 0, 0}},
    {{1, 3, 5, 7, 0, 0, 0, 0}},
    {{1, 5, 7, 11, 0, 0, 0, 0}},
    {{1, 5, 7, 11, 13, 17, 19, 23}},
}};
constexpr std::array<int, 8> UNIT_COUNTS{1, 1, 2, 2, 2, 4, 4, 8};

using Support = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;  // actual orders
using Capacities = std::array<uint16_t, 6>;
using Multiplicity = std::array<uint8_t, 8>;
using MaskFibre = std::array<uint32_t, 8>;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

// Compact self-contained SHA-256.  Certificate streams below are serialized
// byte by byte and never depend on C++ object layout or library hashes.
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

  void u16le(uint16_t value) {
    byte(static_cast<uint8_t>(value));
    byte(static_cast<uint8_t>(value >> 8));
  }

  void u24le(uint32_t value) {
    byte(static_cast<uint8_t>(value));
    byte(static_cast<uint8_t>(value >> 8));
    byte(static_cast<uint8_t>(value >> 16));
  }

  void u32le(uint32_t value) {
    for (int shift = 0; shift < 32; shift += 8)
      byte(static_cast<uint8_t>(value >> shift));
  }

  void u64le(uint64_t value) {
    for (int shift = 0; shift < 64; shift += 8)
      byte(static_cast<uint8_t>(value >> shift));
  }

  std::array<uint8_t, 32> finish_bytes() {
    const uint64_t original_bits = bit_count_;
    byte(0x80);
    while (used_ != 56) byte(0);
    for (int shift = 56; shift >= 0; shift -= 8)
      byte(static_cast<uint8_t>(original_bits >> shift));
    require(used_ == 0, "SHA-256 final block did not close");
    std::array<uint8_t, 32> answer{};
    for (int i = 0; i < 8; ++i)
      for (int offset = 0; offset < 4; ++offset)
        answer[4 * i + offset] =
            static_cast<uint8_t>(state_[i] >> (24 - 8 * offset));
    return answer;
  }

  std::string finish() {
    const auto digest = finish_bytes();
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (uint8_t value : digest)
      out << std::setw(2) << static_cast<int>(value);
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

  static uint32_t rotate(uint32_t value, int count) {
    return (value >> count) | (value << (32 - count));
  }

  void compress(const uint8_t *block) {
    std::array<uint32_t, 64> words{};
    for (int i = 0; i < 16; ++i) words[i] = load_be(block + 4 * i);
    for (int i = 16; i < 64; ++i) {
      const uint32_t s0 = rotate(words[i - 15], 7) ^
                          rotate(words[i - 15], 18) ^
                          (words[i - 15] >> 3);
      const uint32_t s1 = rotate(words[i - 2], 17) ^
                          rotate(words[i - 2], 19) ^
                          (words[i - 2] >> 10);
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
  for (int i = 0; i < static_cast<int>(ORDERS.size()); ++i)
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

int analytic_cardinality(int label, int order, int owner) {
  const int ratio = label * inverse_mod_13(owner) % P;
  const int target = order * ratio % P;
  int period = 0;
  for (int value = -order + 1; value <= order; ++value) {
    int residue = value % P;
    if (residue < 0) residue += P;
    period += residue == target;
  }
  return (C / order) * period;
}

struct Tables {
  // [label][order-index][owner][unit-index]
  std::array<std::array<std::array<MaskFibre, P>, 8>, P> masks{};
  std::array<std::array<std::array<uint8_t, P>, 8>, P> cards{};
  std::string primary_mask_digest;
  std::string literal_base_digest;
};

Tables build_tables() {
  Tables table;
  Sha256 primary_mask_digest;
  Sha256 literal_base_digest;
  for (int label = 1; label < P; ++label)
    for (int oi = 0; oi < static_cast<int>(ORDERS.size()); ++oi) {
      const int order = ORDERS[oi];
      for (int ui = 0; ui < UNIT_COUNTS[oi]; ++ui) {
        const int unit = UNITS[oi][ui];
        const int base = literal_crt_base(label, order, unit);
        literal_base_digest.byte(static_cast<uint8_t>(label));
        literal_base_digest.byte(static_cast<uint8_t>(order));
        literal_base_digest.byte(static_cast<uint8_t>(unit));
        literal_base_digest.u16le(static_cast<uint16_t>(base));
      }
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
                  "sheet cardinality depends on unit choice");
          require(card == analytic_cardinality(label, order, owner),
                  "literal mask and independent period formula disagree");
          // This serialization intentionally matches the Python primary.
          primary_mask_digest.byte(static_cast<uint8_t>(label));
          primary_mask_digest.byte(static_cast<uint8_t>(order));
          primary_mask_digest.byte(static_cast<uint8_t>(unit));
          primary_mask_digest.byte(static_cast<uint8_t>(owner));
          primary_mask_digest.u32le(mask);
        }
        table.cards[label][oi][owner] =
            static_cast<uint8_t>(common_card);
        const int ratio = label * inverse_mod_13(owner) % P;
        require(common_card ==
                    std::popcount(literal_local_mask(
                        ratio, order, UNITS[oi][0], 1)),
                "provider/owner ratio reduction failed");
      }
    }
  table.primary_mask_digest = primary_mask_digest.finish();
  table.literal_base_digest = literal_base_digest.finish();
  return table;
}

bool hereditary_prime_power(const OrderWord &word) {
  int eight = 0;
  int three = 0;
  for (int order : word) {
    eight += order % 8 == 0;
    three += order % 3 == 0;
  }
  return eight >= 2 && three >= 2;
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
                     int coordinate, Sha256 &order_digest,
                     Sha256 &weighted_digest, uint64_t &state_words) {
  if (coordinate == 6) {
    const bool prime_power = hereditary_prime_power(word);
    require(prime_power == hereditary_lcm(word),
            "prime-power and leave-one-out-lcm grammars disagree");
    if (prime_power) {
      words.push_back(word);
      uint64_t fibre = 1;
      for (uint8_t order : word) {
        const uint8_t index = static_cast<uint8_t>(order_index(order));
        order_digest.byte(index);  // matches the Python primary
        weighted_digest.byte(index);
        fibre *= UNIT_COUNTS[index];
      }
      weighted_digest.u64le(fibre);
      state_words += fibre;
    }
    return;
  }
  for (uint8_t order : ORDERS) {
    word[coordinate] = order;
    enumerate_words(words, word, coordinate + 1, order_digest,
                    weighted_digest, state_words);
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
  for (int provider = 0; provider < 6; ++provider) {
    const int oi = order_index(word[provider]);
    for (int owner = 0; owner < 6; ++owner)
      result[owner] += table.cards[support[provider]][oi][support[owner]];
  }
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

uint64_t mix64(uint64_t value) {
  value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ULL;
  value = (value ^ (value >> 27)) * 0x94D049BB133111EBULL;
  return value ^ (value >> 31);
}

struct LocalRow {
  bool feasible = false;
  uint8_t maximum = 0;
  uint32_t maximum_count = 0;
  uint32_t reachable_count = 0;
  uint64_t mask_sum = 0;
  uint64_t mask_xor = 0;
  std::array<uint32_t, 6> forward_layers{};
  std::array<uint32_t, 6> reverse_layers{};
  std::vector<uint32_t> reachable;
};

LocalRow owner_local(const Support &support, const OrderWord &word, int owner,
                     const Tables &table) {
  const auto traverse = [&](bool reverse) {
    std::vector<uint32_t> reachable{0};
    std::array<uint32_t, 6> layers{};
    for (int step = 0; step < 6; ++step) {
      const int provider = reverse ? 5 - step : step;
      const int oi = order_index(word[provider]);
      std::vector<uint32_t> next;
      next.reserve(reachable.size() * UNIT_COUNTS[oi]);
      for (uint32_t partial : reachable)
        for (int ui = 0; ui < UNIT_COUNTS[oi]; ++ui)
          next.push_back(partial |
                         table.masks[support[provider]][oi][owner][ui]);
      std::sort(next.begin(), next.end());
      next.erase(std::unique(next.begin(), next.end()), next.end());
      layers[step] = static_cast<uint32_t>(next.size());
      reachable = std::move(next);
    }
    return std::pair{reachable, layers};
  };

  auto [reachable, forward_layers] = traverse(false);
  auto [reverse_reachable, reverse_layers] = traverse(true);
  require(reachable == reverse_reachable,
          "forward/reverse provider reachability mismatch");
  LocalRow result;
  result.reachable = std::move(reachable);
  result.forward_layers = forward_layers;
  result.reverse_layers = reverse_layers;
  result.reachable_count = static_cast<uint32_t>(result.reachable.size());
  for (uint32_t mask : result.reachable) {
    result.maximum = static_cast<uint8_t>(
        std::max<int>(result.maximum, std::popcount(mask)));
    result.feasible |= mask == FULL;
    result.mask_sum += mask;
    result.mask_xor ^= mix64(mask);
  }
  result.maximum_count = static_cast<uint32_t>(std::count_if(
      result.reachable.begin(), result.reachable.end(), [&](uint32_t mask) {
        return std::popcount(mask) == result.maximum;
      }));
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
  std::array<uint8_t, 6> scores{};
};

Tournament tournament(const Capacities &capacities,
                      const std::array<LocalRow, 6> &locals) {
  std::array<uint8_t, 6> out{};
  Tournament result;
  for (int left = 0; left < 6; ++left)
    for (int right = left + 1; right < 6; ++right) {
      // This is exactly the primary's faithful owner-obligation observable.
      const std::array<int, 3> left_key{
          locals[left].feasible, locals[left].maximum, capacities[left]};
      const std::array<int, 3> right_key{
          locals[right].feasible, locals[right].maximum, capacities[right]};
      int winner = left;
      if (left_key == right_key)
        ++result.ties;  // coordinate tie path: earlier vertex wins
      else if (right_key > left_key) {
        winner = right;
        ++result.flips;
      }
      const int loser = left + right - winner;
      out[winner] |= static_cast<uint8_t>(1U << loser);
    }

  for (int vertex = 0; vertex < 6; ++vertex)
    result.scores[vertex] =
        static_cast<uint8_t>(std::popcount(out[vertex]));

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
  for (int vertex = 0; vertex < 6; ++vertex)
    reach[vertex] |= static_cast<uint8_t>(1U << vertex);
  for (int middle = 0; middle < 6; ++middle)
    for (int source = 0; source < 6; ++source)
      if ((reach[source] >> middle) & 1U) reach[source] |= reach[middle];
  uint8_t assigned = 0;
  for (int root = 0; root < 6; ++root)
    if (!((assigned >> root) & 1U)) {
      ++result.sccs;
      for (int vertex = 0; vertex < 6; ++vertex)
        if (((reach[root] >> vertex) & 1U) &&
            ((reach[vertex] >> root) & 1U))
          assigned |= static_cast<uint8_t>(1U << vertex);
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

std::string orbit_context_histogram(
    const std::map<std::pair<int, int>, uint64_t> &values) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, count] : values) {
    if (!first) out << ' ';
    first = false;
    out << '(' << key.first << ',' << key.second << "):" << count;
  }
  return out.str();
}

template <std::size_t N>
std::string tuple_counter_digest(
    const std::map<std::array<uint8_t, N>, uint64_t> &values) {
  Sha256 digest;
  for (const auto &[key, count] : values) {
    for (uint8_t value : key) digest.byte(value);
    digest.u64le(count);
  }
  return digest.finish();
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
  require(table.literal_base_digest ==
              "0366120df932f618cfabe5e5a5590c40eedfb1674dbb5de0d02bb711392986a7",
          "literal CRT-base digest mismatch");
  require(table.primary_mask_digest ==
              "490f80bf0f0b982e6451e01d04a229210af90c7572eae318c7bdf41a842b4e6f",
          "literal CRT mask table differs from the frozen primary");

  std::vector<OrderWord> words;
  OrderWord scratch{};
  Sha256 order_digest_state;
  Sha256 weighted_digest_state;
  uint64_t state_words_per_support = 0;
  enumerate_words(words, scratch, 0, order_digest_state,
                  weighted_digest_state, state_words_per_support);
  const std::string order_digest = order_digest_state.finish();
  const std::string weighted_digest = weighted_digest_state.finish();
  require(words.size() == 108'813, "hereditary order-word census mismatch");
  require(order_digest ==
              "3fced9d39f778694e85547f7660b3b99c0ee5e67b81c5da9957a30fc128fef24",
          "hereditary order-word digest differs from the frozen primary");
  require(weighted_digest ==
              "b1d7eb335a2656e75c6db90c251b574f4dbb5d263b47f69e0dba58b0c5748a56",
          "weighted literal-state grammar digest mismatch");
  require(state_words_per_support == 167'165'952ULL,
          "literal state-word count mismatch");

  const std::vector<Support> supports = all_supports();
  require(supports.size() == 924, "support census mismatch");
  std::vector<ScalarRow> bank;
  std::set<Support> scalar_supports;
  std::map<Support, int> support_context_count;
  std::map<int, uint64_t> all_support_contexts;
  std::map<int, uint64_t> scalar_feasible_owner;
  std::map<Multiplicity, uint64_t> multiplicities;
  std::set<Capacities> distinct_capacities;
  std::map<int, uint64_t> minimum_slack;
  std::map<int, uint64_t> maximum_slack;
  std::map<int, uint64_t> tight_owners;
  Sha256 scalar_digest_state;
  Sha256 independent_capacity_digest_state;

  for (const Support &support : supports) {
    int support_contexts = 0;
    for (const OrderWord &word : words) {
      const Capacities capacities = capacity_vector(support, word, table);
      const int feasible_count = static_cast<int>(std::count_if(
          capacities.begin(), capacities.end(),
          [](uint16_t value) { return value >= C; }));
      ++scalar_feasible_owner[feasible_count];
      if (feasible_count != 6) continue;
      bank.push_back(ScalarRow{support, word, capacities});
      scalar_supports.insert(support);
      ++support_contexts;
      ++support_context_count[support];
      ++multiplicities[multiplicity(word)];
      distinct_capacities.insert(capacities);
      ++minimum_slack[*std::min_element(capacities.begin(), capacities.end()) -
                      C];
      ++maximum_slack[*std::max_element(capacities.begin(), capacities.end()) -
                      C];
      ++tight_owners[std::count(capacities.begin(), capacities.end(), C)];
      // Exactly the Python primary's scalar-bank serialization.
      for (uint8_t value : support) scalar_digest_state.byte(value);
      for (uint8_t order : word) scalar_digest_state.byte(order);
      for (uint16_t value : capacities) scalar_digest_state.u16le(value);
      // A differently serialized audit digest records order indices too.
      for (uint8_t value : support) independent_capacity_digest_state.byte(value);
      for (uint8_t order : word)
        independent_capacity_digest_state.byte(
            static_cast<uint8_t>(order_index(order)));
      for (uint16_t value : capacities)
        independent_capacity_digest_state.u16le(value);
    }
    ++all_support_contexts[support_contexts];
  }
  const std::string scalar_digest = scalar_digest_state.finish();
  const std::string capacity_digest = independent_capacity_digest_state.finish();
  require(bank.size() == 66'984 && scalar_supports.size() == 854,
          "scalar survivor census mismatch");
  require(scalar_digest ==
              "5bf236e92aedb4f226f75cdc2b5218cf5acdd780c6c1fca34d770d513033993f",
          "scalar bank differs from the frozen primary");
  require(capacity_digest ==
              "55ba3cc44355d47642cee180a763031451d7bf585bfbba8bcda30e9c94297951",
          "independent capacity-bank digest mismatch");
  require(multiplicities.size() == 202,
          "scalar multiplicity-profile census mismatch");
  require(scalar_feasible_owner ==
              std::map<int, uint64_t>{{0, 544'572},
                                      {1, 18'881'520},
                                      {2, 40'784'532},
                                      {3, 29'283'168},
                                      {4, 9'613'632},
                                      {5, 1'368'804},
                                      {6, 66'984}},
          "scalar capacity feasible-owner histogram mismatch");
  require(all_support_contexts ==
              std::map<int, uint64_t>{{0, 70},   {12, 24}, {14, 24},
                                      {22, 24},  {26, 24}, {31, 96},
                                      {32, 30},  {46, 96}, {52, 96},
                                      {53, 24},  {54, 96}, {55, 24},
                                      {63, 96},  {74, 24}, {102, 12},
                                      {132, 24}, {142, 24}, {152, 12},
                                      {173, 24}, {185, 6}, {224, 24},
                                      {226, 24}, {374, 24}, {801, 2}},
          "all-support contexts histogram mismatch");
  require(distinct_capacities.size() == 18'432,
          "distinct capacity-vector census mismatch");
  require(minimum_slack ==
              std::map<int, uint64_t>{{0, 56'286}, {1, 7'884},
                                      {2, 2'682}, {3, 132}},
          "minimum scalar-slack histogram mismatch");
  require(maximum_slack ==
              std::map<int, uint64_t>{{0, 24},      {1, 4'372},
                                      {2, 5'872},   {3, 14'622},
                                      {4, 10'266},  {5, 12'922},
                                      {6, 8'694},   {7, 5'808},
                                      {8, 2'436},   {9, 816},
                                      {10, 480},    {11, 540},
                                      {12, 60},     {16, 72}},
          "maximum scalar-slack histogram mismatch");
  require(tight_owners ==
              std::map<int, uint64_t>{{0, 10'698}, {1, 12'504},
                                      {2, 20'862}, {3, 12'036},
                                      {4, 6'252},  {5, 4'608},
                                      {6, 24}},
          "tight-owner histogram mismatch");

  std::map<int, uint64_t> feasible_contexts;
  std::map<int, uint64_t> feasible_masks;
  std::map<int, uint64_t> maximum_union;
  std::map<int, uint64_t> minimum_owner_maximum;
  std::map<int, uint64_t> reachable_count;
  std::map<int, uint64_t> maximum_mask_count;
  std::map<std::array<uint8_t, 6>, uint64_t> owner_vectors;
  std::map<int, uint64_t> tie_histogram;
  std::map<int, uint64_t> flip_histogram;
  std::map<int, uint64_t> aggregate_scores;
  std::map<std::array<uint8_t, 6>, uint64_t> score_vectors;
  uint64_t feasible_rows = 0;
  uint64_t reachable_total = 0;
  uint32_t reachable_maximum = 0;
  Sha256 primary_owner_digest_state;
  Sha256 reachable_bank_digest_state;
  Sha256 layer_digest_state;

  for (const ScalarRow &row : bank) {
    std::array<LocalRow, 6> locals;
    std::array<uint8_t, 6> owner_vector{};
    int feasible_count = 0;
    uint8_t feasible_mask = 0;
    for (int owner_index = 0; owner_index < 6; ++owner_index) {
      const int owner = row.support[owner_index];
      locals[owner_index] = owner_local(row.support, row.word, owner, table);
      const LocalRow &local = locals[owner_index];
      feasible_count += local.feasible;
      feasible_mask |= static_cast<uint8_t>(local.feasible << owner_index);
      feasible_rows += local.feasible;
      ++maximum_union[local.maximum];
      ++reachable_count[local.reachable_count];
      ++maximum_mask_count[local.maximum_count];
      reachable_total += local.reachable_count;
      reachable_maximum = std::max(reachable_maximum, local.reachable_count);
      owner_vector[owner_index] = local.maximum;

      // Exactly the Python primary's owner-summary stream.
      for (uint8_t value : row.support) primary_owner_digest_state.byte(value);
      for (uint8_t order : row.word) primary_owner_digest_state.byte(order);
      primary_owner_digest_state.byte(static_cast<uint8_t>(owner));
      primary_owner_digest_state.byte(static_cast<uint8_t>(local.feasible));
      primary_owner_digest_state.byte(local.maximum);
      primary_owner_digest_state.u32le(local.reachable_count);
      primary_owner_digest_state.u64le(local.mask_sum);
      primary_owner_digest_state.u64le(local.mask_xor);

      Sha256 local_reachable_digest;
      for (uint32_t mask : local.reachable) local_reachable_digest.u24le(mask);
      const auto local_digest = local_reachable_digest.finish_bytes();
      for (uint8_t value : row.support) reachable_bank_digest_state.byte(value);
      for (uint8_t order : row.word)
        reachable_bank_digest_state.byte(
            static_cast<uint8_t>(order_index(order)));
      reachable_bank_digest_state.byte(static_cast<uint8_t>(owner_index));
      reachable_bank_digest_state.byte(static_cast<uint8_t>(owner));
      reachable_bank_digest_state.u16le(row.capacities[owner_index]);
      reachable_bank_digest_state.u32le(local.reachable_count);
      reachable_bank_digest_state.u32le(local.maximum_count);
      for (uint8_t value : local_digest) reachable_bank_digest_state.byte(value);

      for (uint8_t value : row.support) layer_digest_state.byte(value);
      for (uint8_t order : row.word)
        layer_digest_state.byte(static_cast<uint8_t>(order_index(order)));
      layer_digest_state.byte(static_cast<uint8_t>(owner_index));
      for (uint32_t size : local.forward_layers) layer_digest_state.u32le(size);
      for (uint32_t size : local.reverse_layers) layer_digest_state.u32le(size);
    }
    ++feasible_contexts[feasible_count];
    ++feasible_masks[feasible_mask];
    ++minimum_owner_maximum
         [*std::min_element(owner_vector.begin(), owner_vector.end())];
    ++owner_vectors[owner_vector];

    const Tournament audit = tournament(row.capacities, locals);
    ++tie_histogram[audit.ties];
    ++flip_histogram[audit.flips];
    for (uint8_t score : audit.scores) ++aggregate_scores[score];
    ++score_vectors[audit.scores];
    auto sorted_scores = audit.scores;
    std::sort(sorted_scores.begin(), sorted_scores.end());
    require(sorted_scores == std::array<uint8_t, 6>{0, 1, 2, 3, 4, 5} &&
                audit.triangles == 0 && audit.sccs == 6 && audit.paths == 1,
            "owner-obligation tournament fingerprint mismatch");
  }

  const std::string primary_owner_digest = primary_owner_digest_state.finish();
  const std::string reachable_bank_digest = reachable_bank_digest_state.finish();
  const std::string layer_digest = layer_digest_state.finish();
  require(primary_owner_digest ==
              "4ee53bdb8f3967a45b1fef546135454a5e85419f0ab6382ad22ae91d674618c8",
          "owner summary differs from the frozen primary");
  require(feasible_contexts ==
              std::map<int, uint64_t>{{0, 64'962}, {1, 1'800},
                                      {2, 192}, {4, 30}},
          "feasible-owner deficit mismatch");
  require(feasible_masks ==
              std::map<int, uint64_t>{{0, 64'962}, {1, 349}, {2, 250},
                                      {3, 15},     {4, 301}, {5, 13},
                                      {6, 15},     {8, 301}, {9, 14},
                                      {10, 13},    {12, 10}, {16, 250},
                                      {17, 12},    {18, 10}, {20, 13},
                                      {23, 2},     {24, 15}, {29, 2},
                                      {30, 1},     {32, 349}, {33, 8},
                                      {34, 12},    {36, 14}, {39, 2},
                                      {40, 13},    {45, 7},  {46, 2},
                                      {48, 15},    {51, 10}, {57, 2},
                                      {58, 2}},
          "feasible-owner mask histogram mismatch");
  require(maximum_union ==
              std::map<int, uint64_t>{{12, 72},      {14, 2'136},
                                      {15, 1'644},   {16, 15'876},
                                      {17, 24'420},  {18, 76'296},
                                      {19, 94'872},  {20, 104'592},
                                      {21, 53'040},  {22, 24'948},
                                      {23, 1'704},   {24, 2'304}},
          "maximum reachable-union histogram mismatch");
  require(minimum_owner_maximum ==
              std::map<int, uint64_t>{{12, 72},     {14, 2'016},
                                      {15, 1'050},  {16, 10'668},
                                      {17, 13'716}, {18, 23'932},
                                      {19, 12'540}, {20, 2'906},
                                      {21, 84}},
          "minimum owner-maximum histogram mismatch");
  require(owner_vectors.size() == 20'302,
          "owner maximum-vector census mismatch");
  require(feasible_rows == 2'304 && reachable_total == 101'961'528ULL &&
              reachable_count.size() == 674 && reachable_maximum == 7'728,
          "reachable-bank headline ledger mismatch");
  require(tie_histogram ==
              std::map<int, uint64_t>{{0, 14'112}, {1, 28'488},
                                      {2, 14'496}, {3, 5'700},
                                      {4, 3'408},  {6, 504},
                                      {7, 276}},
          "tournament tie histogram mismatch");
  require(flip_histogram ==
              std::map<int, uint64_t>{{0, 323},    {1, 753},
                                      {2, 1'857},  {3, 3'485},
                                      {4, 5'881},  {5, 8'067},
                                      {6, 9'954},  {7, 10'563},
                                      {8, 9'264},  {9, 7'075},
                                      {10, 4'776}, {11, 2'898},
                                      {12, 1'390}, {13, 543},
                                      {14, 139},   {15, 16}},
          "tournament edge-flip histogram mismatch");
  require(aggregate_scores ==
              std::map<int, uint64_t>{{0, 66'984}, {1, 66'984},
                                      {2, 66'984}, {3, 66'984},
                                      {4, 66'984}, {5, 66'984}},
          "tournament aggregate score histogram mismatch");

  std::set<Support> remaining = scalar_supports;
  std::map<int, uint64_t> orbit_sizes;
  std::map<std::pair<int, int>, uint64_t> orbit_contexts;
  while (!remaining.empty()) {
    const Support seed = *remaining.begin();
    std::set<Support> orbit;
    for (int multiplier = 1; multiplier < P; ++multiplier)
      orbit.insert(multiply_support(seed, multiplier));
    require(std::includes(scalar_supports.begin(), scalar_supports.end(),
                          orbit.begin(), orbit.end()),
            "scalar support bank is not multiplication-invariant");
    require(std::includes(remaining.begin(), remaining.end(),
                          orbit.begin(), orbit.end()),
            "multiplication support orbits overlap");
    std::set<int> context_counts;
    for (const Support &support : orbit)
      context_counts.insert(support_context_count[support]);
    require(context_counts.size() == 1,
            "context count changes within a support orbit");
    for (const Support &support : orbit) remaining.erase(support);
    ++orbit_sizes[orbit.size()];
    ++orbit_contexts[{static_cast<int>(orbit.size()),
                      *context_counts.begin()}];
  }
  require(orbit_sizes ==
              std::map<int, uint64_t>{{2, 1}, {6, 2}, {12, 70}},
          "multiplication orbit-size histogram mismatch");
  require(orbit_contexts ==
              std::map<std::pair<int, int>, uint64_t>{
                  {{2, 801}, 1}, {{6, 32}, 1},  {{6, 185}, 1},
                  {{12, 12}, 2}, {{12, 14}, 2}, {{12, 22}, 2},
                  {{12, 26}, 2}, {{12, 31}, 8}, {{12, 32}, 2},
                  {{12, 46}, 8}, {{12, 52}, 8}, {{12, 53}, 2},
                  {{12, 54}, 8}, {{12, 55}, 2}, {{12, 63}, 8},
                  {{12, 74}, 2}, {{12, 102}, 1}, {{12, 132}, 2},
                  {{12, 142}, 2}, {{12, 152}, 1}, {{12, 173}, 2},
                  {{12, 224}, 2}, {{12, 226}, 2}, {{12, 374}, 2}},
          "support-orbit context histogram mismatch");

  const std::string multiplicity_digest = tuple_counter_digest(multiplicities);
  const std::string owner_vector_digest = tuple_counter_digest(owner_vectors);
  const std::string score_vector_digest = tuple_counter_digest(score_vectors);
  require(multiplicity_digest ==
              "a746f2d8eca4f35bb8b377252a5ff7b3f7ffdbce8499b03990ddefe5566ed70d",
          "scalar multiplicity-ledger digest mismatch");
  require(owner_vector_digest ==
              "844840639474a8afcf382a523ea644be2c2954d5bfb799f5b0df0aa12d00a3e2",
          "owner maximum-vector ledger digest mismatch");
  require(reachable_bank_digest ==
              "daafa6b27d1c56c32defccf49167a472ce6f9318d9f1ac70424a8fd574cbe7d7",
          "full sorted reachable-bank digest mismatch");
  require(layer_digest ==
              "d89d963f9ed57959c9f12c0e047f8bac29ac0632f63b8fcd01f3aba30a052e58",
          "forward/reverse provider-layer digest mismatch");
  require(score_vectors.size() == 720 &&
              score_vector_digest ==
                  "675f74f85d1bc747f92027b63126df0cba1b30fca8dfbef2fa33637de55d25e5",
          "tournament score-vector ledger mismatch");

  std::cout << "scale-twenty-four H6 literal-CRT sorted-union referee\n";
  std::cout << "hereditary grammar: at least two 8-divisible orders and at "
               "least two 3-divisible orders; every leave-one-out lcm audited\n";
  std::cout << "supports " << supports.size() << "; hereditary words "
            << words.size() << "; labelled order contexts "
            << supports.size() * words.size() << '\n';
  std::cout << "literal state words/support " << state_words_per_support
            << "; raw labelled states "
            << supports.size() * state_words_per_support << '\n';
  std::cout << "literal CRT-base SHA256 " << table.literal_base_digest << '\n';
  std::cout << "primary-compatible mask-table SHA256 "
            << table.primary_mask_digest << '\n';
  std::cout << "order-grammar SHA256 " << order_digest << '\n';
  std::cout << "weighted-literal-grammar SHA256 " << weighted_digest << '\n';
  std::cout << "scalar capacity feasible-owner histogram "
            << histogram(scalar_feasible_owner) << '\n';
  std::cout << "scalar contexts " << bank.size() << " on "
            << scalar_supports.size() << " supports; multiplicity patterns "
            << multiplicities.size() << '\n';
  std::cout << "all-support contexts histogram "
            << histogram(all_support_contexts) << '\n';
  std::cout << "multiplication orbit-size histogram "
            << histogram(orbit_sizes) << " (telemetry; no quotient)\n";
  std::cout << "multiplication orbit (size,contexts/support) histogram "
            << orbit_context_histogram(orbit_contexts) << '\n';
  std::cout << "primary-compatible scalar-bank SHA256 " << scalar_digest
            << '\n';
  std::cout << "independent capacity-bank SHA256 " << capacity_digest << '\n';
  std::cout << "scalar-multiplicity SHA256 " << multiplicity_digest << '\n';
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
  std::cout << "feasible-owner-mask histogram " << histogram(feasible_masks)
            << '\n';
  std::cout << "maximum reachable sheet-union histogram "
            << histogram(maximum_union) << '\n';
  std::cout << "minimum owner maximum/context histogram "
            << histogram(minimum_owner_maximum) << '\n';
  std::cout << "distinct owner max-union vectors " << owner_vectors.size()
            << "; owner-vector SHA256 " << owner_vector_digest << '\n';
  std::cout << "reachable banks " << 6 * bank.size() << "; total masks "
            << reachable_total << "; distinct bank-size bins "
            << reachable_count.size() << "; maximum bank "
            << reachable_maximum << '\n';
  std::cout << "primary-compatible owner-summary SHA256 "
            << primary_owner_digest << '\n';
  std::cout << "full reachable-bank SHA256 " << reachable_bank_digest << '\n';
  std::cout << "forward/reverse layer-bank SHA256 " << layer_digest << '\n';
  std::cout << "maximum-mask-count bins " << maximum_mask_count.size()
            << "; maximum-mask-count SHA256 ";
  {
    Sha256 digest;
    for (const auto &[key, value] : maximum_mask_count) {
      digest.u32le(static_cast<uint32_t>(key));
      digest.u64le(value);
    }
    const std::string maximum_mask_digest = digest.finish();
    require(maximum_mask_count.size() == 52 &&
                maximum_mask_digest ==
                    "f84d2990631d8b76fe6c80ce7d80da65058f61fcd7b8a291418009324c02729c",
            "maximum-mask-count ledger mismatch");
    std::cout << maximum_mask_digest << '\n';
  }
  const auto all_six = feasible_contexts.find(6);
  std::cout << "owner-local all-six contexts "
            << (all_six == feasible_contexts.end() ? 0 : all_six->second)
            << '\n';
  std::cout << "tournament carrier owner obligations; observable "
               "(feasible,max-union,capacity); lexicographic switch; "
               "coordinate tie Hamiltonian path\n";
  std::cout << "tournament fingerprints all " << bank.size()
            << " transitive: scores 0,1,2,3,4,5; directed triangles 0; "
               "SCCs 6; Hamiltonian paths 1\n";
  std::cout << "tournament tie-edge histogram " << histogram(tie_histogram)
            << '\n';
  std::cout << "tournament edge-flip histogram " << histogram(flip_histogram)
            << '\n';
  std::cout << "tournament aggregate score histogram "
            << histogram(aggregate_scores) << '\n';
  std::cout << "tournament score vectors " << score_vectors.size()
            << "; score-vector SHA256 " << score_vector_digest << '\n';
  std::cout << "preserved audit: exact owner banks retain FULL-mask "
               "feasibility, maximum deficit, sheet identities, unit "
               "incidence, capacities, and both provider traversal ledgers\n";
  std::cout << "lost audit: owner projection forgets simultaneous cross-owner "
               "unit gluing; the tournament further forgets exact masks, "
               "absolute thresholds, magnitudes, and witness incidence\n";
  std::cout << "challenged vertex assumption: owner proof obligations retain "
               "the terminal predicate; providers, gaps, fixed sections, "
               "section boundaries, wall events, residues, cover arcs, "
               "Fourier modes, and matroid circuits need extra incidence data\n";
  std::cout << "Kakeya/toothpick audit: isolated sheet needles and recursive "
               "mask silhouettes do not preserve shared-unit compatibility; "
               "the sorted union bank is the faithful discrete carrier\n";
  std::cout << "verdict: scalar gate survives, but every context has at most "
               "four feasible owners; the scale-24 common H6 face is empty\n";
}
