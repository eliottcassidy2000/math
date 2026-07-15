// Independent closed-danger-union replay for THM-842.
//
// Unlike the primary verifier, this program does not propagate intersections
// of open safe bands.  At every prefix it reconstructs all CLOSED danger teeth
// from scratch, unions them exactly, and takes their OPEN complementary gaps.
// It then independently regenerates the same exhaustive THM-815/820 tree and
// byte-level state certificate.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
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
    for (size_t i = 0; i < 16; i++)
      w[i] = (uint32_t(block[4 * i]) << 24U) |
             (uint32_t(block[4 * i + 1]) << 16U) |
             (uint32_t(block[4 * i + 2]) << 8U) |
             uint32_t(block[4 * i + 3]);
    for (size_t i = 16; i < 64; i++) {
      uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^
                    (w[i - 15] >> 3U);
      uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^
                    (w[i - 2] >> 10U);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3];
    uint32_t e = h_[4], f = h_[5], g = h_[6], hh = h_[7];
    for (size_t i = 0; i < 64; i++) {
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
      uint8_t ch = static_cast<uint8_t>(static_cast<unsigned char>(raw));
      buffer_[used_++] = ch;
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
    throw overflow_error("rational intermediate exceeds int64");
  return static_cast<long long>(x);
}

struct Rat {
  long long n = 0, d = 1;
  Rat(long long a = 0, long long b = 1) {
    if (b == 0)
      throw invalid_argument("zero rational denominator");
    if (b < 0) {
      a = -a;
      b = -b;
    }
    long long g = gcd(a < 0 ? -a : a, b);
    n = a / g;
    d = b / g;
  }
  bool operator<(Rat const &o) const {
    return (__int128)n * o.d < (__int128)o.n * d;
  }
  bool operator<=(Rat const &o) const { return !(o < *this); }
  bool operator==(Rat const &o) const { return n == o.n && d == o.d; }
};
Rat operator-(Rat a, Rat b) {
  return Rat(narrow_i128((__int128)a.n * b.d - (__int128)b.n * a.d),
             narrow_i128((__int128)a.d * b.d));
}
ostream &operator<<(ostream &out, Rat q) { return out << q.n << '/' << q.d; }
struct I {
  Rat lo, hi;
};
struct Summary {
  vector<I> gaps;
};

// Form every closed tooth of D_u={t:||ut||<=1/13} on [0,1], merge the
// closed union, and return its open complementary components.  The two pieces
// at 0 and 1 are the split copy of the same circular tooth.
Summary summarize(const vector<int> &speeds) {
  if (speeds.empty())
    throw invalid_argument("empty speed prefix");
  set<int> distinct;
  size_t count = 0;
  for (int u : speeds) {
    if (u <= 0 || !distinct.insert(u).second)
      throw runtime_error("speeds must be positive and distinct");
    count += static_cast<size_t>(u) + 1U;
  }
  vector<I> danger;
  danger.reserve(count);
  for (int u : speeds) {
    danger.push_back({Rat(0), Rat(1, 13LL * u)});
    for (int k = 1; k < u; k++)
      danger.push_back(
          {Rat(13LL * k - 1, 13LL * u), Rat(13LL * k + 1, 13LL * u)});
    danger.push_back({Rat(13LL * u - 1, 13LL * u), Rat(1)});
  }
  sort(danger.begin(), danger.end(), [](I const &a, I const &b) {
    if (a.lo == b.lo)
      return a.hi < b.hi;
    return a.lo < b.lo;
  });
  vector<I> merged;
  for (I x : danger) {
    if (!(Rat(0) <= x.lo && x.lo <= x.hi && x.hi <= Rat(1)))
      throw runtime_error("danger tooth outside [0,1]");
    if (merged.empty() || merged.back().hi < x.lo)
      merged.push_back(x);
    else if (merged.back().hi < x.hi)
      merged.back().hi = x.hi;
  }
  if (merged.empty() || !(merged.front().lo == Rat(0)) ||
      !(merged.back().hi == Rat(1)))
    throw runtime_error("closed danger union misses a circle boundary");
  Summary summary;
  Rat previous(0);
  for (I x : merged) {
    if (previous < x.lo)
      summary.gaps.push_back({previous, x.lo});
    if (previous < x.hi)
      previous = x.hi;
  }
  if (previous < Rat(1))
    summary.gaps.push_back({previous, Rat(1)});
  for (size_t i = 0; i < summary.gaps.size(); i++) {
    I x = summary.gaps[i];
    if (!(Rat(0) <= x.lo && x.lo < x.hi && x.hi <= Rat(1)))
      throw runtime_error("invalid open complementary gap");
    if (i && !(summary.gaps[i - 1].hi <= x.lo))
      throw runtime_error("complementary gaps are not ordered");
  }
  return summary;
}

I longest(const vector<I> &gaps) {
  if (gaps.empty())
    throw runtime_error("longest gap requested from empty safe set");
  I best = gaps.front();
  Rat length = best.hi - best.lo;
  for (I x : gaps) {
    Rat candidate = x.hi - x.lo;
    if (length < candidate) {
      length = candidate;
      best = x;
    }
  }
  return best;
}

long long disc_cap(int remaining, Rat length) {
  if (remaining < 1 || remaining > 5 || !(Rat(0) < length) ||
      13 - 2 * remaining <= 0)
    throw invalid_argument("invalid THM-815 discrepancy-cap input");
  __int128 numerator = (__int128)22 * remaining * length.d;
  __int128 denominator =
      (__int128)13 * (13 - 2 * remaining) * length.n;
  if (denominator <= 0)
    throw runtime_error("nonpositive discrepancy denominator");
  return narrow_i128(numerator / denominator);
}

ofstream certificate;
Sha256 certificate_sha;
unsigned long long certificate_rows = 0;
array<int, 5> current_speed{}, current_label{};
array<unsigned long long, 6> nodes{}, branch_b_nodes{};
unsigned long long loose = 0, tight = 0;
set<array<int, 4>> cycle_sets;
vector<int> core;

void emit_certificate(const string &line) {
  certificate.write(line.data(), static_cast<streamsize>(line.size()));
  if (!certificate)
    throw runtime_error("failed while writing state certificate");
  certificate_sha.update(line);
}

void certificate_line(char kind, int depth, char branch,
                      const array<int, 5> &labels, const Summary &summary,
                      long long discrepancy = -1, long long effective = -1) {
  ostringstream line;
  line.imbue(locale::classic());
  line << kind << "|d=" << depth << "|b=" << branch << "|R=";
  for (int q : labels)
    line << q << ',';
  line << "|p=";
  for (size_t i = 0; i < static_cast<size_t>(depth); i++)
    line << current_label[i] << ':' << current_speed[i] << ',';
  line << "|c=" << summary.gaps.size();
  if (!summary.gaps.empty()) {
    I interval = longest(summary.gaps);
    line << "|J=" << interval.lo << ',' << interval.hi
         << "|L=" << (interval.hi - interval.lo);
  } else {
    line << "|J=-|L=0/1";
  }
  line << "|D=" << discrepancy << "|K=" << effective << '\n';
  emit_certificate(line.str());
  certificate_rows++;
}

Summary prefix_summary(int depth) {
  vector<int> speeds = core;
  for (size_t i = 0; i < static_cast<size_t>(depth); i++)
    speeds.push_back(current_speed[i]);
  return summarize(speeds);
}

void recurse(const array<int, 5> &labels, int depth, bool branch_b,
             int v_anchor) {
  if (depth < 2 || depth > 5)
    throw runtime_error("invalid recursive depth");
  nodes.at(static_cast<size_t>(depth))++;
  if (branch_b)
    branch_b_nodes.at(static_cast<size_t>(depth))++;
  Summary summary = prefix_summary(depth);
  if (depth == 5) {
    certificate_line('N', depth, branch_b ? 'B' : 'A', labels, summary);
    if (summary.gaps.empty())
      tight++;
    else
      loose++;
    return;
  }
  if (summary.gaps.empty())
    throw runtime_error("LRC(<=13) nonterminal nonemptiness invariant failed");
  if (branch_b != (current_speed[1] > 2 * current_speed[0]))
    throw runtime_error("THM-820 branch flag mismatch");
  if (branch_b && v_anchor != current_speed[1])
    throw runtime_error("exceptional-branch anchor mismatch");
  I interval = longest(summary.gaps);
  Rat length = interval.hi - interval.lo;
  long long discrepancy = disc_cap(5 - depth, length);
  long long effective = discrepancy;
  int previous = current_speed.at(static_cast<size_t>(depth - 1));
  if (branch_b)
    effective = min<long long>(effective, 4LL * v_anchor);
  else
    effective = min<long long>(effective, 2LL * previous);
  certificate_line('N', depth, branch_b ? 'B' : 'A', labels, summary,
                   discrepancy, effective);
  for (size_t j = 0; j < labels.size(); j++) {
    bool used = false;
    for (size_t k = 0; k < static_cast<size_t>(depth); k++)
      used |= current_label[k] == labels[j];
    if (used)
      continue;
    for (int u = labels[j] + 13; u <= effective; u += 13) {
      if (u <= previous)
        continue;
      current_label.at(static_cast<size_t>(depth)) = labels[j];
      current_speed.at(static_cast<size_t>(depth)) = u;
      recurse(labels, depth + 1, branch_b, v_anchor);
    }
  }
}

int main() {
  const string certificate_path =
      "/tmp/lrc_h5_closed_danger_union_certificate.tsv";
  const string expected_certificate_sha =
      "6524ac6dd2d1f8c59256816c86b95d9ee52cc94766d4d3f993425e7071434a29";
  certificate.open(certificate_path, ios::binary | ios::trunc);
  if (!certificate)
    throw runtime_error("cannot open independent state certificate");
  emit_certificate("LRC-H5-CERT|v=1|delta=1/13|proper=1|sorted=1\n");
  certificate_rows = 1;

  for (int a = 1; a <= 12; a++) {
    array<int, 4> translated{};
    size_t index = 0;
    for (int q : {1, 2, 4, 8})
      translated[index++] = (a * q) % 13;
    sort(translated.begin(), translated.end());
    cycle_sets.insert(translated);
  }
  if (cycle_sets.size() != 12)
    throw runtime_error("exceptional label bank is incomplete");

  auto start = chrono::steady_clock::now();
  unsigned long long core_rows = 0, first_rows = 0;
  unsigned long long exceptional_first = 0, exceptional_killed = 0,
                     exceptional_v_choices = 0;
  long long maximum_root_cap = 0;
  unsigned long long maximum_root_cap_count = 0;

  for (int a = 1; a <= 8; a++)
    for (int b = a + 1; b <= 9; b++)
      for (int c = b + 1; c <= 10; c++)
        for (int d = c + 1; d <= 11; d++)
          for (int e = d + 1; e <= 12; e++) {
            array<int, 5> labels = {a, b, c, d, e};
            bool missing[13] = {};
            for (int q : labels)
              missing[q] = true;
            core.clear();
            int core_maximum = 0;
            for (int q = 1; q <= 12; q++)
              if (!missing[q]) {
                core.push_back(q);
                core_maximum = q;
              }
            if (core.size() != 7 || core_maximum < 7 || core_maximum > 12)
              throw runtime_error("invalid seven-speed core");
            Summary root = summarize(core);
            if (root.gaps.empty())
              throw runtime_error("seven-speed core has no strict-safe gap");
            I root_interval = longest(root.gaps);
            long long root_cap = disc_cap(5, root_interval.hi - root_interval.lo);
            if (root_cap > maximum_root_cap) {
              maximum_root_cap = root_cap;
              maximum_root_cap_count = 1;
            } else if (root_cap == maximum_root_cap) {
              maximum_root_cap_count++;
            }
            certificate_line('C', 0, '-', labels, root, root_cap, 146);
            core_rows++;

            for (size_t xi = 0; xi < labels.size(); xi++)
              for (int x = labels[xi] + 13; x <= 146; x += 13) {
                current_label[0] = labels[xi];
                current_speed[0] = x;
                Summary one = prefix_summary(1);
                if (one.gaps.empty())
                  throw runtime_error("eight-speed prefix has no strict-safe gap");
                I interval = longest(one.gaps);
                long long cap = disc_cap(4, interval.hi - interval.lo);
                certificate_line('N', 1, 'U', labels, one, cap, cap);
                nodes[1]++;
                first_rows++;

                array<int, 4> top{};
                size_t top_index = 0;
                for (size_t j = 0; j < labels.size(); j++)
                  if (j != xi)
                    top[top_index++] = labels[j];
                sort(top.begin(), top.end());
                bool exceptional_labels = cycle_sets.count(top);
                if (exceptional_labels)
                  exceptional_first++;

                for (size_t vi = 0; vi < labels.size(); vi++) {
                  if (vi == xi)
                    continue;
                  for (int v = labels[vi] + 13; v <= cap; v += 13) {
                    if (v <= x)
                      continue;
                    bool branch_b = v > 2 * x;
                    if (branch_b) {
                      if (!exceptional_labels)
                        continue;
                      long long exceptional_cap = 819LL * x / 40;
                      Rat rhs = Rat(15, 104LL * core_maximum) - Rat(1, x);
                      if (rhs.n > 0)
                        exceptional_cap = min<long long>(
                            exceptional_cap,
                            narrow_i128((__int128)7 * rhs.d /
                                        ((__int128)2 * rhs.n)));
                      if (v > exceptional_cap)
                        continue;
                      exceptional_v_choices++;
                    }
                    current_label[1] = labels[vi];
                    current_speed[1] = v;
                    recurse(labels, 2, branch_b, v);
                  }
                }
                if (exceptional_labels && cap <= 2 * x)
                  exceptional_killed++;
              }
          }

  certificate.close();
  if (!certificate)
    throw runtime_error("failed to close independent state certificate");
  string observed_certificate_sha = certificate_sha.finish();

  const array<unsigned long long, 6> expected_nodes = {
      0, 40590, 612221, 111675, 7255, 9};
  const array<unsigned long long, 6> expected_branch_b = {0, 0, 415, 178, 1,
                                                          0};
  if (core_rows != 792 || first_rows != 40590 || nodes != expected_nodes ||
      branch_b_nodes != expected_branch_b || exceptional_first != 984 ||
      exceptional_killed != 897 || exceptional_v_choices != 415 ||
      maximum_root_cap != 146 || maximum_root_cap_count != 1 || loose != 9 ||
      tight != 0 || certificate_rows != 772543 ||
      observed_certificate_sha != expected_certificate_sha)
    throw runtime_error("frozen independent replay mismatch");

  double seconds =
      chrono::duration<double>(chrono::steady_clock::now() - start).count();
  cerr << "runtime_s=" << seconds
       << " method=closed_danger_union_recomputed_at_every_prefix\n";
  cout << "THM842_HAMMING_FIVE_CLOSED_DANGER_UNION_REPLAY\n";
  cout << "method=union_closed_danger_teeth_then_take_open_complementary_gaps\n";
  cout << "core_rows=" << core_rows << " root_global_cap=" << maximum_root_cap
       << " root_cap_maximizers=" << maximum_root_cap_count << '\n';
  cout << "first_prefixes=" << first_rows
       << " exceptional_label_prefixes=" << exceptional_first
       << " exceptional_killed_before_v=" << exceptional_killed
       << " exceptional_v_choices=" << exceptional_v_choices << '\n';
  cout << "nodes depth0..5=";
  for (auto count : nodes)
    cout << count << ',';
  cout << " full_loose=" << loose << " full_tight=" << tight << '\n';
  cout << "exceptional_branch_nodes depth0..5=";
  for (auto count : branch_b_nodes)
    cout << count << ',';
  cout << '\n';
  cout << "certificate=" << certificate_path << " rows=" << certificate_rows
       << '\n';
  cout << "certificate_sha256=" << observed_certificate_sha << '\n';
  cout << "PASS: independent closed-danger-union replay byte-matches the "
          "primary THM-842 certificate\n";
}
