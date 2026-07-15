// Independent closed-danger-union replay for the scale-one Hamming-six tree.
//
// The primary recursion propagates intersections of open safe bands.  This
// replay instead rebuilds the CLOSED danger union of every expanded prefix
// from scratch, takes its OPEN complementary gaps, and drives the same exact
// THM-815 next-speed recursion.  Candidate shortcuts are verified by directly
// subtracting one closed danger comb from the reconstructed parent gaps.
//
// Each root receives a canonical SHA-256 certificate over expanded states,
// full-safe-gap witnesses, early-cap witnesses, and terminal verdicts.  The
// root hashes are deterministic across optimization levels and can later be
// combined in root-index order without storing a multi-gigabyte text trace.
//
// This source deliberately does not include or call the primary verifier.

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

class Sha256 {
  array<uint32_t, 8> h_ = {0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U,
                           0xa54ff53aU, 0x510e527fU, 0x9b05688cU,
                           0x1f83d9abU, 0x5be0cd19U};
  array<uint8_t, 64> buffer_{};
  size_t used_ = 0;
  uint64_t bytes_ = 0;
  bool finished_ = false;

  static uint32_t rotr(uint32_t x, unsigned n) {
    return (x >> n) | (x << (32U - n));
  }

  void transform(const uint8_t *block) {
    static constexpr array<uint32_t, 64> K = {
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
    array<uint32_t, 64> w{};
    for (size_t i = 0; i < 16; ++i)
      w[i] = (uint32_t(block[4 * i]) << 24U) |
             (uint32_t(block[4 * i + 1]) << 16U) |
             (uint32_t(block[4 * i + 2]) << 8U) |
             uint32_t(block[4 * i + 3]);
    for (size_t i = 16; i < 64; ++i) {
      uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^
                    (w[i - 15] >> 3U);
      uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^
                    (w[i - 2] >> 10U);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3];
    uint32_t e = h_[4], f = h_[5], g = h_[6], hh = h_[7];
    for (size_t i = 0; i < 64; ++i) {
      uint32_t s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
      uint32_t ch = (e & f) ^ ((~e) & g);
      uint32_t t1 = hh + s1 + ch + K[i] + w[i];
      uint32_t s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
      uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
      uint32_t t2 = s0 + maj;
      hh = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }
    h_[0] += a;
    h_[1] += b;
    h_[2] += c;
    h_[3] += d;
    h_[4] += e;
    h_[5] += f;
    h_[6] += g;
    h_[7] += hh;
  }

public:
  void update(const string &s) {
    if (finished_)
      throw runtime_error("SHA-256 update after finish");
    if (numeric_limits<uint64_t>::max() - bytes_ < s.size())
      throw runtime_error("SHA-256 byte count overflow");
    bytes_ += static_cast<uint64_t>(s.size());
    for (char raw : s) {
      buffer_[used_++] =
          static_cast<uint8_t>(static_cast<unsigned char>(raw));
      if (used_ == 64) {
        transform(buffer_.data());
        used_ = 0;
      }
    }
  }

  string finish() {
    if (finished_)
      throw runtime_error("SHA-256 finish called twice");
    finished_ = true;
    uint64_t bits = bytes_ * 8U;
    buffer_[used_++] = 0x80U;
    if (used_ > 56) {
      while (used_ < 64)
        buffer_[used_++] = 0;
      transform(buffer_.data());
      used_ = 0;
    }
    while (used_ < 56)
      buffer_[used_++] = 0;
    for (int shift = 56; shift >= 0; shift -= 8)
      buffer_[used_++] = uint8_t(bits >> unsigned(shift));
    transform(buffer_.data());
    ostringstream out;
    out.imbue(locale::classic());
    out << hex << setfill('0');
    for (uint32_t x : h_)
      out << setw(8) << x;
    return out.str();
  }
};

long long narrow_i128(__int128 x) {
  if (x < numeric_limits<long long>::min() ||
      x > numeric_limits<long long>::max())
    throw overflow_error("integer intermediate exceeds int64");
  return static_cast<long long>(x);
}

struct Rat {
  long long n = 0;
  long long d = 1;

  Rat(long long a = 0, long long b = 1) {
    if (b <= 0)
      throw invalid_argument("Rat denominator must be positive");
    long long g = gcd(a < 0 ? -a : a, b);
    n = a / g;
    d = b / g;
  }

  bool operator<(Rat const &other) const {
    return (__int128)n * other.d < (__int128)other.n * d;
  }
  bool operator<=(Rat const &other) const { return !(other < *this); }
  bool operator==(Rat const &other) const {
    return n == other.n && d == other.d;
  }
};

Rat operator-(Rat a, Rat b) {
  long long g = gcd(a.d, b.d);
  __int128 nn = (__int128)a.n * (b.d / g) -
                (__int128)b.n * (a.d / g);
  __int128 dd = (__int128)(a.d / g) * b.d;
  return Rat(narrow_i128(nn), narrow_i128(dd));
}

ostream &operator<<(ostream &out, Rat value) {
  return out << value.n << '/' << value.d;
}

struct Interval {
  Rat lo;
  Rat hi;
};

bool same_intervals(const vector<Interval> &a, const vector<Interval> &b) {
  if (a.size() != b.size())
    return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (!(a[i].lo == b[i].lo) || !(a[i].hi == b[i].hi))
      return false;
  return true;
}

// Reconstruct all closed danger teeth from scratch and return their open
// complementary gaps.  Closed teeth that merely touch are merged.
vector<Interval> closed_danger_complement(const vector<int> &speeds) {
  if (speeds.empty())
    throw invalid_argument("empty speed prefix");
  set<int> distinct;
  size_t tooth_count = 0;
  for (int u : speeds) {
    if (u <= 0 || !distinct.insert(u).second)
      throw runtime_error("prefix speeds must be positive and distinct");
    tooth_count += static_cast<size_t>(u) + 1U;
  }

  vector<Interval> danger;
  danger.reserve(tooth_count);
  for (int u : speeds) {
    danger.push_back({Rat(0), Rat(1, 13LL * u)});
    for (int k = 1; k < u; ++k)
      danger.push_back(
          {Rat(13LL * k - 1, 13LL * u),
           Rat(13LL * k + 1, 13LL * u)});
    danger.push_back({Rat(13LL * u - 1, 13LL * u), Rat(1)});
  }
  sort(danger.begin(), danger.end(), [](Interval const &a, Interval const &b) {
    if (a.lo == b.lo)
      return a.hi < b.hi;
    return a.lo < b.lo;
  });

  vector<Interval> merged;
  for (Interval tooth : danger) {
    if (merged.empty() || merged.back().hi < tooth.lo)
      merged.push_back(tooth);
    else if (merged.back().hi < tooth.hi)
      merged.back().hi = tooth.hi;
  }
  if (merged.empty() || !(merged.front().lo == Rat(0)) ||
      !(merged.back().hi == Rat(1)))
    throw runtime_error("closed danger union misses a circle boundary");

  vector<Interval> gaps;
  Rat cursor(0);
  for (Interval covered : merged) {
    if (cursor < covered.lo)
      gaps.push_back({cursor, covered.lo});
    if (cursor < covered.hi)
      cursor = covered.hi;
  }
  if (cursor < Rat(1))
    gaps.push_back({cursor, Rat(1)});
  return gaps;
}

Interval longest_gap(const vector<Interval> &gaps) {
  if (gaps.empty())
    throw runtime_error("longest gap requested from empty complement");
  Interval best = gaps.front();
  Rat best_length = best.hi - best.lo;
  for (Interval gap : gaps) {
    Rat length = gap.hi - gap.lo;
    if (best_length < length) {
      best = gap;
      best_length = length;
    }
  }
  return best;
}

long long discrepancy_cap(int remaining, Rat length) {
  if (remaining < 1 || remaining > 6 || length.n <= 0 ||
      13 - 2 * remaining <= 0)
    throw invalid_argument("invalid discrepancy-cap input");
  __int128 numerator = (__int128)22 * remaining * length.d;
  __int128 denominator =
      (__int128)13 * (13 - 2 * remaining) * length.n;
  __int128 quotient = numerator / denominator;
  return narrow_i128(quotient);
}

bool cap_below_candidate(int remaining, Rat lo, Rat hi, int candidate) {
  __int128 length_numerator =
      (__int128)hi.n * lo.d - (__int128)lo.n * hi.d;
  if (length_numerator <= 0)
    throw runtime_error("nonpositive gap in cap witness");
  __int128 lhs = (__int128)22 * remaining * lo.d * hi.d;
  __int128 rhs = (__int128)candidate * 13 * (13 - 2 * remaining) *
                 length_numerator;
  return lhs < rhs;
}

long long first_overlapping_safe_gap(Rat left, int u) {
  // (13k+12)/(13u) > left.
  __int128 numerator = (__int128)13 * u * left.n -
                       (__int128)12 * left.d;
  __int128 denominator = (__int128)13 * left.d;
  if (numerator < 0)
    return 0;
  return narrow_i128(numerator / denominator + 1);
}

long long last_overlapping_safe_gap(Rat right, int u) {
  // (13k+1)/(13u) < right.  The -1 implements strictness exactly.
  __int128 numerator = (__int128)13 * u * right.n - right.d;
  __int128 denominator = (__int128)13 * right.d;
  if (numerator <= 0)
    return -1;
  return narrow_i128((numerator - 1) / denominator);
}

struct GapWitness {
  size_t parent_component = 0;
  int k = -1;
  Interval gap;
};

bool find_full_safe_gap(const vector<Interval> &parent, int u,
                        GapWitness &witness) {
  for (size_t index = 0; index < parent.size(); ++index) {
    Interval component = parent[index];
    long long first = max<long long>(
        0, first_overlapping_safe_gap(component.lo, u));
    long long last = min<long long>(
        u - 1, last_overlapping_safe_gap(component.hi, u));
    for (long long k = first; k <= last; ++k) {
      Interval safe = {Rat(13LL * k + 1, 13LL * u),
                       Rat(13LL * (k + 1) - 1, 13LL * u)};
      if (component.lo <= safe.lo && safe.hi <= component.hi) {
        witness = {index, static_cast<int>(k), safe};
        return true;
      }
    }
  }
  return false;
}

enum class SubtractVerdict { materialized, dead_completion, loose_terminal };

// Subtract one closed danger comb from the reconstructed open parent gaps.
// The gaps between consecutive closed teeth are the open intervals used below.
SubtractVerdict subtract_closed_comb(const vector<Interval> &parent, int u,
                                     int remaining_after, int least_future,
                                     vector<Interval> &child,
                                     GapWitness &witness) {
  child.clear();
  child.reserve(parent.size() + 4);
  for (size_t index = 0; index < parent.size(); ++index) {
    Interval component = parent[index];
    long long first = max<long long>(
        0, first_overlapping_safe_gap(component.lo, u));
    long long last = min<long long>(
        u - 1, last_overlapping_safe_gap(component.hi, u));
    for (long long k = first; k <= last; ++k) {
      Interval safe = {Rat(13LL * k + 1, 13LL * u),
                       Rat(13LL * (k + 1) - 1, 13LL * u)};
      Rat lo = component.lo < safe.lo ? safe.lo : component.lo;
      Rat hi = component.hi < safe.hi ? component.hi : safe.hi;
      if (!(lo < hi))
        continue;
      Interval piece = {lo, hi};
      if (remaining_after == 0) {
        witness = {index, static_cast<int>(k), piece};
        return SubtractVerdict::loose_terminal;
      }
      if (cap_below_candidate(remaining_after, lo, hi, least_future)) {
        witness = {index, static_cast<int>(k), piece};
        return SubtractVerdict::dead_completion;
      }
      child.push_back(piece);
    }
  }
  return SubtractVerdict::materialized;
}

struct RootResult {
  int root_index = -1;
  array<int, 6> missing{};
  array<unsigned long long, 7> nodes{};
  array<unsigned long long, 7> dead{};
  array<unsigned long long, 7> full_prunes{};
  array<unsigned long long, 7> early_prunes{};
  unsigned long long candidate_edges = 0;
  unsigned long long covers = 0;
  unsigned long long loose = 0;
  vector<string> cover_rows;
  string certificate_sha256;
  double runtime_seconds = 0;
};

class RootReplay {
  int root_index_;
  array<int, 6> missing_;
  vector<int> core_;
  array<int, 6> labels_{};
  array<int, 6> speeds_{};
  RootResult result_;
  Sha256 certificate_;

  string prefix_word(int depth) const {
    ostringstream out;
    out.imbue(locale::classic());
    for (int i = 0; i < depth; ++i)
      out << labels_[static_cast<size_t>(i)] << ':'
          << speeds_[static_cast<size_t>(i)] << ',';
    return out.str();
  }

  string missing_word() const {
    ostringstream out;
    out.imbue(locale::classic());
    for (int q : missing_)
      out << q << ',';
    return out.str();
  }

  void hash_line(const string &line) { certificate_.update(line); }

  void hash_node(int depth, const vector<Interval> &gaps, long long cap,
                 size_t candidates) {
    ostringstream line;
    line.imbue(locale::classic());
    line << "N|d=" << depth << "|R=" << missing_word()
         << "|p=" << prefix_word(depth) << "|g=" << gaps.size();
    for (Interval gap : gaps)
      line << '|' << gap.lo << ',' << gap.hi;
    line << "|D=" << cap << "|C=" << candidates << '\n';
    hash_line(line.str());
  }

  void hash_witness(char kind, int child_depth, int u, int label,
                    int remaining_after, int least_future,
                    const GapWitness &witness) {
    Rat length = witness.gap.hi - witness.gap.lo;
    long long cap = remaining_after ? discrepancy_cap(remaining_after, length)
                                    : -1;
    ostringstream line;
    line.imbue(locale::classic());
    line << kind << "|d=" << child_depth << "|R=" << missing_word()
         << "|p=" << prefix_word(child_depth) << "|u=" << label << ':' << u
         << "|pc=" << witness.parent_component << "|k=" << witness.k
         << "|J=" << witness.gap.lo << ',' << witness.gap.hi
         << "|r=" << remaining_after << "|y=" << least_future
         << "|D=" << cap << '\n';
    hash_line(line.str());
  }

  vector<int> prefix_speeds(int depth) const {
    vector<int> out = core_;
    for (int i = 0; i < depth; ++i)
      out.push_back(speeds_[static_cast<size_t>(i)]);
    return out;
  }

  long long least_future_speed(int depth_after, int last) const {
    long long best = numeric_limits<long long>::max();
    for (int label : missing_) {
      bool used = false;
      for (int i = 0; i < depth_after; ++i)
        used |= labels_[static_cast<size_t>(i)] == label;
      if (used)
        continue;
      long long speed = label + 13;
      if (speed <= last)
        speed += 13 * ((last - speed) / 13 + 1);
      best = min(best, speed);
    }
    if (best == numeric_limits<long long>::max() ||
        best > numeric_limits<int>::max())
      throw runtime_error("missing or overflowing future speed");
    return best;
  }

  string cover_word() const {
    ostringstream out;
    out.imbue(locale::classic());
    out << "missing=" << missing_word() << " ordered=" << prefix_word(6);
    return out.str();
  }

  void recurse(int depth, int last,
               const vector<Interval> *expected_gaps = nullptr) {
    result_.nodes.at(static_cast<size_t>(depth))++;
    vector<Interval> gaps = closed_danger_complement(prefix_speeds(depth));
    if (expected_gaps && !same_intervals(gaps, *expected_gaps))
      throw runtime_error("closed-danger reconstruction disagrees with direct subtraction");

    if (gaps.empty()) {
      if (depth != 6)
        throw runtime_error("nonterminal closed-danger cover contradicts LRC(<=13)");
      result_.covers++;
      result_.cover_rows.push_back(cover_word());
      hash_node(depth, gaps, -2, 0);
      return;
    }
    if (depth == 6) {
      result_.loose++;
      hash_node(depth, gaps, -1, 0);
      return;
    }

    Interval longest = longest_gap(gaps);
    long long cap = discrepancy_cap(6 - depth, longest.hi - longest.lo);
    if (cap > 10000000)
      throw runtime_error("independent cap exceeded 10^7 guard");
    vector<pair<int, int>> candidates;
    for (int label : missing_) {
      bool used = false;
      for (int i = 0; i < depth; ++i)
        used |= labels_[static_cast<size_t>(i)] == label;
      if (used)
        continue;
      long long speed = label + 13;
      if (speed <= last)
        speed += 13 * ((last - speed) / 13 + 1);
      for (; speed <= cap; speed += 13)
        candidates.emplace_back(static_cast<int>(speed), label);
    }
    sort(candidates.begin(), candidates.end());
    if (candidates.empty())
      result_.dead.at(static_cast<size_t>(depth))++;
    result_.candidate_edges += candidates.size();
    hash_node(depth, gaps, cap, candidates.size());

    for (auto [u, label] : candidates) {
      labels_.at(static_cast<size_t>(depth)) = label;
      speeds_.at(static_cast<size_t>(depth)) = u;
      int child_depth = depth + 1;
      int remaining_after = 5 - depth;

      GapWitness witness;
      if (depth >= 2 && find_full_safe_gap(gaps, u, witness)) {
        result_.nodes.at(static_cast<size_t>(child_depth))++;
        result_.full_prunes.at(static_cast<size_t>(child_depth))++;
        if (depth == 5)
          result_.loose++;
        else
          result_.dead.at(static_cast<size_t>(child_depth))++;
        int least_future = remaining_after
                               ? static_cast<int>(least_future_speed(child_depth, u))
                               : -1;
        if (remaining_after &&
            !cap_below_candidate(remaining_after, witness.gap.lo,
                                 witness.gap.hi, least_future))
          throw runtime_error("full safe gap does not imply the claimed cap contradiction");
        hash_witness('F', child_depth, u, label, remaining_after,
                     least_future, witness);
        continue;
      }

      int least_future = remaining_after
                             ? static_cast<int>(least_future_speed(child_depth, u))
                             : -1;
      vector<Interval> child;
      SubtractVerdict verdict = subtract_closed_comb(
          gaps, u, remaining_after, least_future, child, witness);
      if (verdict != SubtractVerdict::materialized) {
        result_.nodes.at(static_cast<size_t>(child_depth))++;
        result_.early_prunes.at(static_cast<size_t>(child_depth))++;
        if (verdict == SubtractVerdict::dead_completion) {
          result_.dead.at(static_cast<size_t>(child_depth))++;
          if (!cap_below_candidate(remaining_after, witness.gap.lo,
                                   witness.gap.hi, least_future))
            throw runtime_error("invalid independent early-cap witness");
          hash_witness('E', child_depth, u, label, remaining_after,
                       least_future, witness);
        } else {
          result_.loose++;
          hash_witness('L', child_depth, u, label, 0, -1, witness);
        }
        continue;
      }
      recurse(child_depth, u, &child);
    }
  }

public:
  RootReplay(int root_index, array<int, 6> missing)
      : root_index_(root_index), missing_(missing) {
    bool removed[13] = {};
    for (int q : missing_)
      removed[q] = true;
    for (int q = 1; q <= 12; ++q)
      if (!removed[q])
        core_.push_back(q);
    if (core_.size() != 6)
      throw runtime_error("invalid H6 deletion core");
  }

  RootResult run() {
    ostringstream header;
    header.imbue(locale::classic());
    header << "LRC-H6-CLOSED-DANGER|v=1|delta=1/13|root=" << root_index_
           << "|R=" << missing_word() << '\n';
    hash_line(header.str());
    auto started = chrono::steady_clock::now();
    recurse(0, 0);
    result_.runtime_seconds =
        chrono::duration<double>(chrono::steady_clock::now() - started).count();
    result_.root_index = root_index_;
    result_.missing = missing_;
    result_.certificate_sha256 = certificate_.finish();
    return result_;
  }
};

string counts_word(const array<unsigned long long, 7> &counts) {
  ostringstream out;
  out.imbue(locale::classic());
  for (auto value : counts)
    out << value << ',';
  return out.str();
}

void check_diagnostic(const RootResult &result) {
  if (result.root_index == 0) {
    const array<unsigned long long, 7> nodes =
        {1, 66, 4841, 199427, 9146, 14, 0};
    const array<unsigned long long, 7> dead =
        {0, 0, 0, 196130, 9135, 14, 0};
    const array<unsigned long long, 7> full =
        {0, 0, 0, 162063, 0, 0, 0};
    const array<unsigned long long, 7> early =
        {0, 0, 0, 34067, 9135, 14, 0};
    if (result.nodes != nodes || result.dead != dead ||
        result.full_prunes != full || result.early_prunes != early ||
        result.candidate_edges != 213494 || result.covers != 0 ||
        result.loose != 0)
      throw runtime_error("root 0 diagnostic census mismatch");
  }
  if (result.root_index == 286) {
    const array<unsigned long long, 7> nodes =
        {1, 116, 14520, 1021458, 11289, 94, 3};
    const array<unsigned long long, 7> dead =
        {0, 0, 0, 1018566, 11248, 92, 0};
    const array<unsigned long long, 7> full =
        {0, 0, 0, 924432, 0, 0, 0};
    const array<unsigned long long, 7> early =
        {0, 0, 0, 94134, 11248, 92, 2};
    const string cover =
        "missing=1,3,5,7,9,11, ordered=1:14,3:16,5:18,7:20,9:22,11:24,";
    if (result.nodes != nodes || result.dead != dead ||
        result.full_prunes != full || result.early_prunes != early ||
        result.candidate_edges != 1047480 || result.covers != 1 ||
        result.loose != 2 || result.cover_rows != vector<string>{cover})
      throw runtime_error("root 286 diagnostic census mismatch");
  }
  if (result.root_index == 462) {
    const array<unsigned long long, 7> nodes =
        {1, 156, 25406, 2361302, 89674, 570, 0};
    const array<unsigned long long, 7> dead =
        {0, 0, 0, 2337077, 89438, 570, 0};
    const array<unsigned long long, 7> full =
        {0, 0, 0, 2130175, 0, 0, 0};
    const array<unsigned long long, 7> early =
        {0, 0, 0, 206902, 89438, 570, 0};
    if (result.nodes != nodes || result.dead != dead ||
        result.full_prunes != full || result.early_prunes != early ||
        result.candidate_edges != 2477108 || result.covers != 0 ||
        result.loose != 0)
      throw runtime_error("root 462 diagnostic census mismatch");
  }
}

void check_all_root_aggregate(
    const array<unsigned long long, 7> &nodes,
    const array<unsigned long long, 7> &dead,
    const array<unsigned long long, 7> &full,
    const array<unsigned long long, 7> &early,
    unsigned long long edges, unsigned long long covers,
    unsigned long long loose,
    const vector<pair<int, string>> &cover_rows) {
  const array<unsigned long long, 7> expected_nodes =
      {924, 83881, 8906315, 559202706, 12671505, 53812, 21};
  const array<unsigned long long, 7> expected_dead =
      {0, 0, 0, 555565824, 12638291, 53792, 0};
  const array<unsigned long long, 7> expected_full =
      {0, 0, 0, 495797163, 0, 0, 0};
  const array<unsigned long long, 7> expected_early =
      {0, 0, 0, 59768661, 12638291, 53792, 20};
  const vector<pair<int, string>> expected_cover_rows = {{
      286,
      "missing=1,3,5,7,9,11, ordered=1:14,3:16,5:18,7:20,9:22,11:24,"}};
  if (nodes != expected_nodes || dead != expected_dead ||
      full != expected_full || early != expected_early ||
      edges != 580918240 || covers != 1 || loose != 20 ||
      cover_rows != expected_cover_rows)
    throw runtime_error("all-root independent replay census mismatch");
}

int main(int argc, char **argv) {
  int root_start = 0;
  int root_limit = 924;
  for (int i = 1; i < argc; ++i) {
    string argument = argv[i];
    if (argument == "--root-start" && i + 1 < argc)
      root_start = stoi(argv[++i]);
    else if (argument == "--root-limit" && i + 1 < argc)
      root_limit = stoi(argv[++i]);
    else
      throw runtime_error("unknown or incomplete argument: " + argument);
  }
  if (root_start < 0 || root_limit < 0 || root_start + root_limit > 924)
    throw runtime_error("invalid root range");

  Sha256 sha_test;
  sha_test.update("abc");
  if (sha_test.finish() !=
      "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
    throw runtime_error("local SHA-256 self-test failed");

  cout << "H6_CLOSED_DANGER_UNION_REPLAY\n";
  cout << "method=rebuild_closed_danger_union_at_expanded_prefixes"
          "+subtract_closed_candidate_comb\n";
  cout << "root_start=" << root_start << " root_limit=" << root_limit << '\n';

  Sha256 manifest;
  array<unsigned long long, 7> total_nodes{}, total_dead{}, total_full{},
      total_early{};
  unsigned long long total_edges = 0, total_covers = 0, total_loose = 0;
  vector<pair<int, string>> all_cover_rows;
  int root_index = 0;
  int completed = 0;
  auto all_started = chrono::steady_clock::now();
  for (int a = 1; a <= 7; ++a)
    for (int b = a + 1; b <= 8; ++b)
      for (int c = b + 1; c <= 9; ++c)
        for (int d = c + 1; d <= 10; ++d)
          for (int e = d + 1; e <= 11; ++e)
            for (int f = e + 1; f <= 12; ++f) {
              int this_root = root_index++;
              if (this_root < root_start ||
                  this_root >= root_start + root_limit)
                continue;
              RootResult result =
                  RootReplay(this_root, {a, b, c, d, e, f}).run();
              check_diagnostic(result);
              for (size_t i = 0; i < 7; ++i) {
                total_nodes[i] += result.nodes[i];
                total_dead[i] += result.dead[i];
                total_full[i] += result.full_prunes[i];
                total_early[i] += result.early_prunes[i];
              }
              total_edges += result.candidate_edges;
              total_covers += result.covers;
              total_loose += result.loose;

              ostringstream line;
              line.imbue(locale::classic());
              line << "root=" << result.root_index << " missing=";
              for (int q : result.missing)
                line << q << ',';
              line << " nodes=" << counts_word(result.nodes)
                   << " dead=" << counts_word(result.dead)
                   << " full_prunes=" << counts_word(result.full_prunes)
                   << " early_prunes=" << counts_word(result.early_prunes)
                   << " edges=" << result.candidate_edges
                   << " covers=" << result.covers
                   << " loose=" << result.loose
                   << " certificate_sha256=" << result.certificate_sha256;
              manifest.update(line.str() + "\n");
              cout << line.str() << '\n';
              cerr << "root_runtime root=" << result.root_index
                   << " seconds=" << fixed << setprecision(6)
                   << result.runtime_seconds << '\n';
              for (string const &cover : result.cover_rows)
                all_cover_rows.emplace_back(result.root_index, cover);
              for (string const &cover : result.cover_rows)
                cout << "COVER root=" << result.root_index << ' ' << cover
                     << '\n';
              ++completed;
            }

  if (root_start == 0 && root_limit == 924)
    check_all_root_aggregate(total_nodes, total_dead, total_full, total_early,
                             total_edges, total_covers, total_loose,
                             all_cover_rows);

  double seconds = chrono::duration<double>(chrono::steady_clock::now() -
                                             all_started)
                       .count();
  cout << "TOTAL roots=" << completed << " nodes=" << counts_word(total_nodes)
       << " dead=" << counts_word(total_dead)
       << " full_prunes=" << counts_word(total_full)
       << " early_prunes=" << counts_word(total_early)
       << " edges=" << total_edges << " covers=" << total_covers
       << " loose=" << total_loose
       << " manifest_sha256=" << manifest.finish() << '\n';
  cerr << "total_runtime roots=" << completed << " seconds=" << fixed
       << setprecision(6) << seconds << '\n';
  cout << "PASS: independent closed-danger-union H6 replay completed\n";
}
