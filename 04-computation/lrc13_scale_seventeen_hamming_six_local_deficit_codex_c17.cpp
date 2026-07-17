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
#include <utility>
#include <vector>

// Exact primary certificate for the primitive proper AP-centred common-scale-
// seventeen Hamming-six sheet bank.  Literal CRT masks and the full hereditary
// divisor/unit grammar are reconstructed without a quotient.  Scalar capacity
// isolates the quadratic-residue and nonresidue supports with all providers at
// order seventeen.  Exact owner-local union DP then misses at least one of the
// seventeen sheets in every one of the twelve owner rows.

using Labels = std::array<uint8_t, 6>;
using OrderWord = std::array<uint8_t, 6>;
using StateWord = std::array<uint8_t, 6>;

static constexpr int P = 13;
static constexpr int C = 17;
static constexpr uint32_t FULL = (1U << C) - 1U;
static constexpr std::array<int, 2> DIVISORS{1, 17};
static constexpr std::array<int, 17> STATE_ORDER{
    1, 17, 17, 17, 17, 17, 17, 17, 17,
    17, 17, 17, 17, 17, 17, 17, 17};
static constexpr std::array<int, 17> STATE_UNIT{
    0, 1, 2, 3, 4, 5, 6, 7, 8,
    9, 10, 11, 12, 13, 14, 15, 16};
static constexpr std::array<int, 3> STATE_BEGIN{0, 1, 17};
static constexpr std::array<int, 2> ORDER_REP{0, 1};

static std::array<std::array<std::array<uint32_t, 13>, 17>, 13> MASK{};

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition) fail(message);
}

void fnv_byte(uint64_t &digest, uint8_t value) {
  digest ^= value;
  digest *= 1'099'511'628'211ULL;
}

void fnv_u32(uint64_t &digest, uint32_t value) {
  for (int shift = 0; shift < 32; shift += 8)
    fnv_byte(digest, static_cast<uint8_t>(value >> shift));
}

int centered(int value, int modulus) {
  int residue = value % modulus;
  if (residue < 0) residue += modulus;
  return 2 * residue > modulus ? residue - modulus : residue;
}

int inverse_mod(int value, int modulus) {
  for (int candidate = 1; candidate < modulus; ++candidate)
    if (value * candidate % modulus == 1) return candidate;
  fail("modular inverse missing");
}

// Literal CRT search is deliberately used in the primary implementation.
int crt_base(int label, int state) {
  const int order = STATE_ORDER[state];
  const int unit = STATE_UNIT[state];
  for (int value = 0; value < P * order; ++value)
    if (value % P == order * label % P &&
        value % order == unit % order)
      return value;
  fail("CRT base missing");
}

uint32_t local_mask(int label, int state, int owner) {
  const int order = STATE_ORDER[state];
  const int base = crt_base(label, state);
  const int owner_inverse = inverse_mod(owner, P);
  uint32_t result = 0;
  for (int sheet = 0; sheet < C; ++sheet) {
    const int value = centered(
        base * (owner_inverse + P * sheet), P * order);
    if (-order < value && value <= order)
      result |= 1U << sheet;
  }
  return result;
}

std::pair<int, int> state_range(int order_index) {
  require(0 <= order_index && order_index < 2, "bad order index");
  return {STATE_BEGIN[order_index], STATE_BEGIN[order_index + 1]};
}

bool hereditary(const OrderWord &word) {
  for (int omitted = 0; omitted < 6; ++omitted) {
    int residual = 1;
    for (int coordinate = 0; coordinate < 6; ++coordinate)
      if (coordinate != omitted)
        residual = std::lcm(residual, DIVISORS[word[coordinate]]);
    if (residual != C) return false;
  }
  return true;
}

void enumerate_orders(std::vector<OrderWord> &answer, OrderWord &word,
                      int coordinate) {
  if (coordinate == 6) {
    if (hereditary(word)) answer.push_back(word);
    return;
  }
  for (int order_index = 0; order_index < 2; ++order_index) {
    word[coordinate] = static_cast<uint8_t>(order_index);
    enumerate_orders(answer, word, coordinate + 1);
  }
}

uint64_t fibre_size(const OrderWord &word) {
  uint64_t result = 1;
  for (int order_index : word) {
    const auto [begin, end] = state_range(order_index);
    result *= static_cast<uint64_t>(end - begin);
  }
  return result;
}

void enumerate_states(const OrderWord &word, StateWord &states, int coordinate,
                      uint64_t &count, uint64_t &digest) {
  if (coordinate == 6) {
    ++count;
    for (uint8_t state : states) fnv_byte(digest, state + 1U);
    fnv_byte(digest, 0xffU);
    return;
  }
  const auto [begin, end] = state_range(word[coordinate]);
  for (int state = begin; state < end; ++state) {
    states[coordinate] = static_cast<uint8_t>(state);
    enumerate_states(word, states, coordinate + 1, count, digest);
  }
}

bool scalar_capacity(const Labels &labels, const OrderWord &word, int owner) {
  int capacity = 0;
  for (int provider = 0; provider < 6; ++provider)
    capacity += std::popcount(
        MASK[labels[provider]][ORDER_REP[word[provider]]][owner]);
  return capacity >= C;
}

struct LocalAudit {
  bool feasible = false;
  int maximum_union = 0;
  std::size_t reachable_count = 0;
  std::size_t maximum_count = 0;
  std::array<uint64_t, C> missing_at_maximum{};
  std::vector<uint32_t> sorted_reachable;
};

LocalAudit owner_local_audit(const Labels &labels, const OrderWord &word,
                             int owner) {
  std::vector<uint32_t> reachable{0};
  std::array<uint8_t, 1U << C> seen{};
  for (int provider = 0; provider < 6; ++provider) {
    std::vector<uint32_t> next;
    const auto [begin, end] = state_range(word[provider]);
    for (uint32_t partial : reachable)
      for (int state = begin; state < end; ++state) {
        const uint32_t joined =
            partial | MASK[labels[provider]][state][owner];
        if (!seen[joined]) {
          seen[joined] = 1;
          next.push_back(joined);
        }
      }
    for (uint32_t mask : next) seen[mask] = 0;
    reachable = std::move(next);
  }

  LocalAudit result;
  result.reachable_count = reachable.size();
  for (uint32_t mask : reachable) {
    result.maximum_union =
        std::max(result.maximum_union, std::popcount(mask));
    result.feasible |= mask == FULL;
  }
  for (uint32_t mask : reachable)
    if (std::popcount(mask) == result.maximum_union) {
      ++result.maximum_count;
      if (result.maximum_union == C - 1) {
        const uint32_t missing = FULL ^ mask;
        require(std::has_single_bit(missing),
                "maximum mask does not miss exactly one sheet");
        ++result.missing_at_maximum[std::countr_zero(missing)];
      }
    }
  std::sort(reachable.begin(), reachable.end());
  result.sorted_reachable = std::move(reachable);
  return result;
}

Labels multiply_support(const Labels &labels, int multiplier) {
  Labels result{};
  for (int coordinate = 0; coordinate < 6; ++coordinate)
    result[coordinate] =
        static_cast<uint8_t>(multiplier * labels[coordinate] % P);
  std::sort(result.begin(), result.end());
  return result;
}

std::string histogram_string(const std::map<int, uint64_t> &histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto &[key, value] : histogram) {
    if (!first) out << ' ';
    first = false;
    out << key << ':' << value;
  }
  return out.str();
}

std::string support_string(const Labels &labels) {
  std::ostringstream out;
  for (int coordinate = 0; coordinate < 6; ++coordinate) {
    if (coordinate) out << ',';
    out << static_cast<int>(labels[coordinate]);
  }
  return out.str();
}

int main() {
  require(STATE_ORDER.size() == STATE_UNIT.size(), "state arrays disagree");
  require(STATE_ORDER[0] == 1 && STATE_UNIT[0] == 0,
          "order-one state mismatch");
  std::set<int> units17;
  for (int state = 1; state < 17; ++state) {
    require(STATE_ORDER[state] == 17, "D17 state order mismatch");
    require(std::gcd(STATE_UNIT[state], 17) == 1, "D17 nonunit listed");
    units17.insert(STATE_UNIT[state]);
  }
  std::set<int> expected_units17;
  for (int unit = 1; unit < 17; ++unit) expected_units17.insert(unit);
  require(units17 == expected_units17, "D17 unit grammar incomplete");

  uint64_t mask_digest = 14'695'981'039'346'656'037ULL;
  for (int label = 1; label < P; ++label)
    for (int state = 0; state < 17; ++state)
      for (int owner = 1; owner < P; ++owner) {
        MASK[label][state][owner] = local_mask(label, state, owner);
        fnv_u32(mask_digest, MASK[label][state][owner]);
      }

  for (int order_index = 0; order_index < 2; ++order_index)
    for (int label = 1; label < P; ++label)
      for (int owner = 1; owner < P; ++owner) {
        const auto [begin, end] = state_range(order_index);
        std::set<int> cardinalities;
        for (int state = begin; state < end; ++state)
          cardinalities.insert(std::popcount(MASK[label][state][owner]));
        require(cardinalities.size() == 1,
                "unit-dependent scalar mask cardinality");
      }

  std::vector<OrderWord> order_words;
  OrderWord order_scratch{};
  enumerate_orders(order_words, order_scratch, 0);
  uint64_t order_digest = 14'695'981'039'346'656'037ULL;
  for (const OrderWord &word : order_words) {
    for (uint8_t order : word) fnv_byte(order_digest, order);
    fnv_byte(order_digest, 0xffU);
  }
  uint64_t weighted_states = 0;
  for (const OrderWord &word : order_words) weighted_states += fibre_size(word);
  uint64_t literal_states = 0;
  uint64_t grammar_digest = 14'695'981'039'346'656'037ULL;
  StateWord state_scratch{};
  for (const OrderWord &word : order_words)
    enumerate_states(word, state_scratch, 0, literal_states, grammar_digest);
  require(order_words.size() == 57, "hereditary order-word census mismatch");
  require(weighted_states == 24'137'472,
          "weighted state-word census mismatch");
  require(literal_states == weighted_states,
          "literal and weighted state censuses disagree");
  require(mask_digest == 0xd7751aac7824587dULL,
          "literal mask-table digest mismatch");
  require(order_digest == 0xe5de68761e24680eULL,
          "order-grammar digest mismatch");
  require(grammar_digest == 0xa30487f2a4edede5ULL,
          "literal state-grammar digest mismatch");

  uint64_t supports = 0;
  uint64_t scalar_contexts = 0;
  uint64_t scalar_digest = 14'695'981'039'346'656'037ULL;
  uint64_t owner_digest = 14'695'981'039'346'656'037ULL;
  std::set<Labels> scalar_supports;
  std::vector<std::pair<Labels, OrderWord>> scalar_bank;
  std::map<int, uint64_t> d17_count_histogram;

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
              for (const OrderWord &word : order_words) {
                bool scalar = true;
                for (int owner : labels)
                  if (!scalar_capacity(labels, word, owner)) {
                    scalar = false;
                    break;
                  }
                if (!scalar) continue;
                ++scalar_contexts;
                scalar_supports.insert(labels);
                scalar_bank.push_back({labels, word});
                const int d17_count = std::accumulate(word.begin(), word.end(), 0);
                ++d17_count_histogram[d17_count];
                for (uint8_t label : labels) fnv_byte(scalar_digest, label);
                for (uint8_t order : word) fnv_byte(scalar_digest, order);
              }
            }

  const Labels qr{1, 3, 4, 9, 10, 12};
  const Labels nqr{2, 5, 6, 7, 8, 11};
  const std::set<Labels> expected_supports{qr, nqr};
  const OrderWord all_d17{1, 1, 1, 1, 1, 1};
  require(supports == 924, "support census mismatch");
  require(supports * order_words.size() == 52'668,
          "labelled order-context census mismatch");
  // Arithmetic correction: 924*24,137,472 is 22,303,024,128.  The earlier
  // 22,304,224,128 transcription is not divisible by 924 and is false.
  require(supports * weighted_states == 22'303'024'128ULL,
          "raw labelled state-context census mismatch");
  require(scalar_contexts == 2 && scalar_supports == expected_supports,
          "scalar support classification mismatch");
  require(d17_count_histogram == std::map<int, uint64_t>{{6, 2}},
          "scalar order classification mismatch");
  require(scalar_digest == 0xf70fa5a77389e62bULL,
          "scalar-bank digest mismatch");
  for (const auto &[labels, word] : scalar_bank) {
    require(word == all_d17, "scalar context is not all-D17");
    require(labels == qr || labels == nqr, "unexpected scalar support");
  }

  std::map<int, uint64_t> reachable_count_histogram;
  std::map<int, uint64_t> maximum_count_histogram;
  std::array<uint64_t, C> aggregate_missing{};
  uint64_t owner_rows = 0;
  uint64_t feasible_rows = 0;
  std::map<int, uint64_t> tournament_ties;
  std::map<int, uint64_t> tournament_flips;
  for (const auto &[labels, word] : scalar_bank) {
    std::array<std::pair<bool, int>, 6> owner_summaries{};
    int owner_index = 0;
    for (int owner : labels) {
      const LocalAudit audit = owner_local_audit(labels, word, owner);
      ++owner_rows;
      feasible_rows += audit.feasible;
      require(audit.maximum_union == 16,
              "owner row does not have exact sixteen-sheet maximum");
      ++reachable_count_histogram[static_cast<int>(audit.reachable_count)];
      ++maximum_count_histogram[static_cast<int>(audit.maximum_count)];
      for (int sheet = 0; sheet < C; ++sheet)
        aggregate_missing[sheet] += audit.missing_at_maximum[sheet];
      for (uint8_t label : labels) fnv_byte(owner_digest, label);
      fnv_byte(owner_digest, static_cast<uint8_t>(owner));
      fnv_byte(owner_digest, static_cast<uint8_t>(audit.maximum_union));
      for (uint32_t mask : audit.sorted_reachable) fnv_u32(owner_digest, mask);
      owner_summaries[owner_index++] = {audit.feasible, audit.maximum_union};
    }
    int ties = 0;
    int flips = 0;
    std::array<int, 6> scores{};
    for (int left = 0; left < 6; ++left)
      for (int right = left + 1; right < 6; ++right) {
        int winner = left;
        if (owner_summaries[left] == owner_summaries[right]) {
          ++ties;
        } else if (owner_summaries[right] > owner_summaries[left]) {
          winner = right;
          ++flips;
        }
        ++scores[winner];
      }
    std::sort(scores.begin(), scores.end());
    require(scores == std::array<int, 6>{0, 1, 2, 3, 4, 5},
            "owner tournament score histogram mismatch");
    ++tournament_ties[ties];
    ++tournament_flips[flips];
  }
  require(owner_rows == 12 && feasible_rows == 0,
          "owner-local deficit census mismatch");
  require(reachable_count_histogram ==
              std::map<int, uint64_t>{{38'540, 12}},
          "reachable-mask count histogram mismatch");
  require(maximum_count_histogram == std::map<int, uint64_t>{{16, 12}},
          "maximum-mask count histogram mismatch");
  const std::array<uint64_t, C> expected_missing{
      12, 11, 11, 11, 12, 11, 11, 11, 12,
      11, 11, 11, 12, 11, 11, 11, 12};
  require(aggregate_missing == expected_missing,
          "maximum-mask missing-sheet profile mismatch");
  require(owner_digest == 0x020e891b7358b525ULL,
          "owner-reachable-bank digest mismatch");

  // Multiplication by F_13^* swaps/preserves the QR/NQR color classes.  This
  // is telemetry only: neither scalar nor local scans used the quotient.
  std::set<Labels> orbit;
  for (int multiplier = 1; multiplier < P; ++multiplier)
    orbit.insert(multiply_support(qr, multiplier));
  require(orbit == expected_supports,
          "QR/NQR multiplication orbit mismatch");

  // Pair observable: exact owner summaries (feasible,max-union).  Every pair
  // ties, so the fixed coordinate Hamiltonian path orients all fifteen edges.
  // The tournament is transitive with scores 5,4,3,2,1,0, no cycle, six SCCs,
  // and one Hamiltonian path; it carries no deficit information beyond the
  // absolute summary retained at each owner vertex.
  require(tournament_ties == std::map<int, uint64_t>{{15, 2}} &&
              tournament_flips == std::map<int, uint64_t>{{0, 2}},
          "owner tournament tie/flip census mismatch");

  std::cout << "scale-seventeen AP-centred Hamming-six local deficit\n";
  std::cout << "divisor grammar 1,17; literal states 17\n";
  std::cout << "supports " << supports << "; hereditary order words "
            << order_words.size() << "; labelled order contexts "
            << supports * order_words.size() << '\n';
  std::cout << "state words/support " << weighted_states
            << "; raw labelled states " << supports * weighted_states << '\n';
  std::cout << "arithmetic audit earlier value 22304224128 is false: not divisible by 924; exact product is 22303024128\n";
  std::cout << "mask-table FNV64 " << std::hex << mask_digest
            << "; order-grammar FNV64 " << order_digest
            << "; literal-state-grammar FNV64 " << grammar_digest
            << std::dec << '\n';
  std::cout << "scalar contexts " << scalar_contexts
            << "; all-D17; supports QR [" << support_string(qr)
            << "] and NQR [" << support_string(nqr) << "]\n";
  std::cout << "scalar-bank FNV64 " << std::hex << scalar_digest
            << std::dec << "; multiplication orbit size " << orbit.size()
            << " (telemetry; no quotient)\n";
  std::cout << "owner-local rows " << owner_rows << "; feasible "
            << feasible_rows << "; maximum-union histogram 16:12\n";
  std::cout << "reachable-mask count histogram "
            << histogram_string(reachable_count_histogram) << '\n';
  std::cout << "maximum-mask count histogram "
            << histogram_string(maximum_count_histogram) << '\n';
  std::cout << "missing-sheet counts among maximum masks";
  for (int sheet = 0; sheet < C; ++sheet)
    std::cout << ' ' << sheet << ':' << aggregate_missing[sheet];
  std::cout << '\n';
  std::cout << "owner-reachable-bank FNV64 " << std::hex << owner_digest
            << std::dec << '\n';
  std::cout << "global literal unit fibres 0; common-sheet survivors 0\n";
  std::cout << "tournament pair observable exact owner (feasible,max-union) summaries; all ties use coordinate path 0->1->2->3->4->5\n";
  std::cout << "tournament fingerprints two transitive: scores 0,1,2,3,4,5; cycles 0; SCCs 6; Hamiltonian paths 1; ties "
            << histogram_string(tournament_ties) << "; flips "
            << histogram_string(tournament_flips) << '\n';
  std::cout << "challenged vertices owner obligations with absolute maxima preserve the deficit; the oriented tournament forgets it, while providers, residues, divisor words, and isolated sheets destroy shared-unit incidence\n";
  std::cout << "local D17 mask table at owner one (units 1..16; ratios 1..12 in hex)\n";
  for (int state = 1; state < 17; ++state) {
    std::cout << "  e=" << STATE_UNIT[state] << ':';
    for (int ratio = 1; ratio < P; ++ratio)
      std::cout << ' ' << std::hex << MASK[ratio][state][1] << std::dec;
    std::cout << '\n';
  }
}
