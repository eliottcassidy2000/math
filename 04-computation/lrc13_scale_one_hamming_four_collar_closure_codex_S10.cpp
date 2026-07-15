// Exact replay for THM-815, scale-one Hamming-four collar closure.
//
// Build and run:
//   c++ -O3 -std=c++17 \
//     04-computation/lrc13_scale_one_hamming_four_collar_closure_codex_S10.cpp \
//     -o /tmp/lrc13_h4 && /tmp/lrc13_h4
//
// Every coverage comparison uses integers.  A strict-safe component endpoint
// is always one of the input fractions (13k+/-1)/(13u), so unreduced pairs of
// 64-bit integers suffice; cross-products use __int128.  SHA-256 is implemented
// locally so the certificate has no nonstandard library dependency.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

class Sha256 {
  array<uint32_t, 8> state_ = {0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U,
                               0xa54ff53aU, 0x510e527fU, 0x9b05688cU,
                               0x1f83d9abU, 0x5be0cd19U};
  array<uint8_t, 64> data_{};
  size_t used_ = 0;
  uint64_t transformed_bits_ = 0;

  static uint32_t rotr(uint32_t x, unsigned n) {
    return (x >> n) | (x << (32 - n));
  }

  void transform() {
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
    for (int i = 0; i < 16; ++i) {
      w[i] = (uint32_t(data_[4 * i]) << 24) |
             (uint32_t(data_[4 * i + 1]) << 16) |
             (uint32_t(data_[4 * i + 2]) << 8) | uint32_t(data_[4 * i + 3]);
    }
    for (int i = 16; i < 64; ++i) {
      uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
      uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    auto [a, b, c, d, e, f, g, h] = state_;
    for (int i = 0; i < 64; ++i) {
      uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
      uint32_t ch = (e & f) ^ (~e & g);
      uint32_t t1 = h + S1 + ch + K[i] + w[i];
      uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
      uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
      uint32_t t2 = S0 + maj;
      h = g; g = f; f = e; e = d + t1;
      d = c; c = b; b = a; a = t1 + t2;
    }
    state_[0] += a; state_[1] += b; state_[2] += c; state_[3] += d;
    state_[4] += e; state_[5] += f; state_[6] += g; state_[7] += h;
  }

 public:
  void update(const uint8_t* bytes, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      data_[used_++] = bytes[i];
      if (used_ == 64) {
        transform();
        transformed_bits_ += 512;
        used_ = 0;
      }
    }
  }

  string hex_digest() const {
    Sha256 copy = *this;
    uint64_t total_bits = copy.transformed_bits_ + 8 * copy.used_;
    copy.data_[copy.used_++] = 0x80;
    if (copy.used_ > 56) {
      while (copy.used_ < 64) copy.data_[copy.used_++] = 0;
      copy.transform();
      copy.used_ = 0;
    }
    while (copy.used_ < 56) copy.data_[copy.used_++] = 0;
    for (int i = 7; i >= 0; --i)
      copy.data_[copy.used_++] = uint8_t(total_bits >> (8 * i));
    copy.transform();
    ostringstream out;
    out << hex << setfill('0');
    for (uint32_t word : copy.state_)
      out << setw(8) << word;
    return out.str();
  }
};

struct Rat { long long n, d; };
struct Interval { Rat lo, hi; };

static bool rat_lt(Rat a, Rat b) {
  return (__int128)a.n * b.d < (__int128)b.n * a.d;
}
static bool rat_le(Rat a, Rat b) { return !rat_lt(b, a); }
static Rat rat_max(Rat a, Rat b) { return rat_lt(a, b) ? b : a; }
static Rat rat_min(Rat a, Rat b) { return rat_lt(a, b) ? a : b; }

static vector<Interval> safe_bands(int speed) {
  vector<Interval> out;
  out.reserve(speed);
  for (int k = 0; k < speed; ++k)
    out.push_back({{13LL * k + 1, 13LL * speed},
                   {13LL * (k + 1) - 1, 13LL * speed}});
  return out;
}

static vector<Interval> intersect_unions(const vector<Interval>& a,
                                         const vector<Interval>& b) {
  vector<Interval> out;
  size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    Rat lo = rat_max(a[i].lo, b[j].lo);
    Rat hi = rat_min(a[i].hi, b[j].hi);
    if (rat_lt(lo, hi)) out.push_back({lo, hi});
    if (rat_le(a[i].hi, b[j].hi)) ++i;
    else ++j;
  }
  return out;
}

static vector<Interval> strict_safe_components(vector<int> speeds) {
  sort(speeds.begin(), speeds.end());
  vector<Interval> current = {{{0, 1}, {1, 1}}};
  for (int speed : speeds) {
    current = intersect_unions(current, safe_bands(speed));
    if (current.empty()) break;
  }
  return current;
}

struct Failure {
  bool failed;
  int index;
  long long surplus, denominator, lhs;
  Rat lo, hi;
};

static Failure first_containment_failure(const vector<Interval>& components,
                                          int speed) {
  for (int i = 0; i < int(components.size()); ++i) {
    Rat lo = components[i].lo, hi = components[i].hi;
    long long denominator = 2 * lo.d * hi.d;
    long long centre_num = lo.n * hi.d + hi.n * lo.d;
    long long halfwidth_num = hi.n * lo.d - lo.n * hi.d;
    long long residue = (speed * centre_num) % denominator;
    if (residue < 0) residue += denominator;
    long long distance_num = min(residue, denominator - residue);
    long long lhs = 13 * (distance_num + speed * halfwidth_num);
    if (lhs > denominator)
      return {true, i, lhs - denominator, denominator, lhs, lo, hi};
  }
  return {false, -1, 0, 1, 0, {0, 1}, {0, 1}};
}

static void digest_u64(Sha256& digest, uint64_t value) {
  uint8_t bytes[8];
  for (int i = 7; i >= 0; --i) {
    bytes[i] = uint8_t(value);
    value >>= 8;
  }
  digest.update(bytes, 8);
}

static void digest_row(Sha256& digest, const array<int, 4>& labels,
                       const array<int, 4>& speeds, const Failure& failure) {
  for (int x : labels) digest_u64(digest, x);
  for (int x : speeds) digest_u64(digest, x);
  digest_u64(digest, failure.index);
  digest_u64(digest, failure.surplus);
  digest_u64(digest, failure.denominator);
  digest_u64(digest, failure.lhs);
  digest_u64(digest, failure.lo.n);
  digest_u64(digest, failure.lo.d);
  digest_u64(digest, failure.hi.n);
  digest_u64(digest, failure.hi.d);
}

static int mod13(int x) {
  x %= 13;
  return x < 0 ? x + 13 : x;
}

static vector<int> residue_values(int label, int low, int high) {
  int first_height = max(1, (low - label + 12) / 13);
  int last_height = (high - label) / 13;
  vector<int> out;
  for (int h = first_height; h <= last_height; ++h)
    out.push_back(label + 13 * h);
  return out;
}

static long long gcdll(long long a, long long b) {
  while (b) { long long r = a % b; a = b; b = r; }
  return a < 0 ? -a : a;
}

static pair<long long, long long> exact_maximin(const vector<int>& speeds) {
  set<int> denominators;
  for (int u : speeds) denominators.insert(2 * u);
  for (size_t i = 0; i < speeds.size(); ++i)
    for (size_t j = i + 1; j < speeds.size(); ++j) {
      denominators.insert(speeds[i] + speeds[j]);
      denominators.insert(abs(speeds[i] - speeds[j]));
    }
  denominators.erase(0);
  long long best_n = 0, best_d = 1;
  for (int d : denominators)
    for (int n = 1; n < d; ++n) {
      int value_n = d;
      for (int u : speeds) {
        int r = int((1LL * u * n) % d);
        value_n = min(value_n, min(r, d - r));
      }
      if ((__int128)value_n * best_d > (__int128)best_n * d) {
        best_n = value_n;
        best_d = d;
      }
    }
  long long g = gcdll(best_n, best_d);
  return {best_n / g, best_d / g};
}

static void finite_band_lemma() {
  const vector<int> centres = {2, 3, 4, 5, 6, 7, 8,
                               9, 10, 11, 12, 15, 16, 17};
  int rows = 0, survivors = 0;
  array<int, 4> survivor{};
  for (size_t i = 0; i < centres.size(); ++i)
    for (size_t j = i; j < centres.size(); ++j)
      for (size_t k = j; k < centres.size(); ++k)
        for (size_t l = k; l < centres.size(); ++l) {
          array<int, 4> a = {centres[i], centres[j], centres[k], centres[l]};
          ++rows;
          long long lower = 1, upper = 1;
          int residue_product = 1;
          for (int x : a) {
            lower *= x - 1;
            upper *= x + 1;
            residue_product = residue_product * (x % 13) % 13;
          }
          if (lower <= 16 && 16 < upper && residue_product == 1) {
            ++survivors;
            survivor = a;
          }
        }
  if (rows != 2380 || survivors != 1 || survivor != array<int, 4>{2, 2, 2, 5})
    throw runtime_error("finite band lemma failed");
  cout << "FOUR_CYCLE_BAND_LEMMA\n";
  cout << "arrow=provider_to_owner z=provider_label/owner_label "
          "lambda=provider_speed/owner_speed\n";
  cout << "allowed_integer_centres={2..12,15,16,17} sorted_multisets=" << rows
       << " survivors=" << survivors << " survivor=(2,2,2,5)\n";
  cout << "normal_form_labels=a*{1,2,4,8} "
          "cycle=a->8a->4a->2a->a residue_word=(5,2,2,2)\n\n";
}

static bool left_handoff(int provider, int owner) {
  int z = mod13((provider % 13) * [] (int x) {
    for (int y = 1; y < 13; ++y) if (x * y % 13 == 1) return y;
    return 0;
  }(owner % 13));
  for (int m = -20; m <= 20; ++m) {
    long long scaled = 1LL * z * owner - 2LL * provider - 13LL * m * owner;
    if (-owner < scaled && scaled <= owner) return true;
  }
  return false;
}

struct TournamentFingerprint {
  map<int, int> score_histogram;
  int directed_triangles = 0;
  vector<int> scc_sizes;
  int hamiltonian_paths = 0;
};

static TournamentFingerprint tournament_fingerprint(
    const array<array<bool, 4>, 4>& edge) {
  TournamentFingerprint out;
  for (int i = 0; i < 4; ++i) {
    int score = 0;
    for (int j = 0; j < 4; ++j) score += edge[i][j];
    ++out.score_histogram[score];
  }
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j)
      for (int k = j + 1; k < 4; ++k)
        if ((edge[i][j] && edge[j][k] && edge[k][i]) ||
            (edge[j][i] && edge[i][k] && edge[k][j]))
          ++out.directed_triangles;
  array<array<bool, 4>, 4> reach = edge;
  for (int i = 0; i < 4; ++i) reach[i][i] = true;
  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j]);
  array<bool, 4> used{};
  for (int i = 0; i < 4; ++i) if (!used[i]) {
    int size = 0;
    for (int j = 0; j < 4; ++j)
      if (!used[j] && reach[i][j] && reach[j][i]) { used[j] = true; ++size; }
    out.scc_sizes.push_back(size);
  }
  sort(out.scc_sizes.rbegin(), out.scc_sizes.rend());
  array<int, 4> p = {0, 1, 2, 3};
  do {
    if (edge[p[0]][p[1]] && edge[p[1]][p[2]] && edge[p[2]][p[3]])
      ++out.hamiltonian_paths;
  } while (next_permutation(p.begin(), p.end()));
  return out;
}

template <class T>
static string map_string(const map<int, T>& values) {
  ostringstream out;
  out << '{';
  bool first = true;
  for (auto [k, v] : values) {
    if (!first) out << ',';
    first = false;
    out << k << ':' << v;
  }
  out << '}';
  return out.str();
}

static void tournament_audit() {
  const array<int, 4> speeds = {79, 54, 30, 34};
  array<array<bool, 4>, 4> live{};
  int live_edges = 0;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      if (i != j && left_handoff(speeds[i], speeds[j])) {
        live[i][j] = true;
        ++live_edges;
      }
  if (live_edges != 4 || !live[0][3] || !live[3][2] ||
      !live[2][1] || !live[1][0])
    throw runtime_error("method-witness handoff cycle failed");

  array<array<bool, 4>, 4> forward{}, reverse{};
  int silent_pairs = 0, edge_flips = 0;
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j) {
      if (live[i][j] || live[j][i]) {
        forward[i][j] = reverse[i][j] = live[i][j];
        forward[j][i] = reverse[j][i] = live[j][i];
      } else {
        ++silent_pairs;
        forward[i][j] = true;
        reverse[j][i] = true;
        ++edge_flips;
      }
    }
  auto f = tournament_fingerprint(forward);
  auto r = tournament_fingerprint(reverse);
  if (silent_pairs != 2 || edge_flips != 2 ||
      f.score_histogram != map<int, int>{{1, 2}, {2, 2}} ||
      r.score_histogram != f.score_histogram ||
      f.directed_triangles != 2 || r.directed_triangles != 2 ||
      f.scc_sizes != vector<int>{4} || r.scc_sizes != vector<int>{4} ||
      f.hamiltonian_paths != 5 || r.hamiltonian_paths != 5)
    throw runtime_error("tournament fingerprint failed");

  vector<int> packet = {3, 5, 6, 7, 9, 10, 11, 12, 30, 34, 54, 79};
  auto maximum = exact_maximin(packet);
  if (maximum != pair<long long, long long>{3, 19})
    throw runtime_error("method-witness maximin failed");

  cout << "TOURNAMENT_TELEMETRY\n";
  cout << "method_witness_labels=(1,2,4,8) replacements=(79,54,30,34) "
          "live_cycle=1->8->4->2->1 silent_diagonals={1,4},{2,8}\n";
  cout << "pair_observable=left_handoff_difference switch=reverse_two_silent_ties "
       << "edge_flips=" << edge_flips << " shared_hamiltonian_path=(1,8,4,2)\n";
  cout << "both_score_histogram=" << map_string(f.score_histogram)
       << " directed_triangles=" << f.directed_triangles
       << " scc_sizes=[4] hamiltonian_paths=" << f.hamiltonian_paths << '\n';
  cout << "method_witness_M=3/19 interpretation=bare_tournament_is_telemetry\n\n";
}

int main() {
  cout << "LRC13 SCALE-ONE HAMMING-FOUR COLLAR CLOSURE - EXACT REPLAY\n";
  cout << "arithmetic=integer endpoint_fractions=unreduced cross_products=__int128\n";
  cout << "scope=all four distinct labels and all four positive lift heights\n\n";

  finite_band_lemma();

  const array<int, 4> multipliers = {7, 4, 2, 1};
  long long exceptional_rows = 0, exceptional_tight = 0, exceptional_cores = 0;
  map<array<int, 4>, long long> exceptional_label_histogram;
  long long exceptional_best_surplus = 1, exceptional_best_denominator = 0;
  array<int, 8> exceptional_best_row{};
  Failure exceptional_best_failure{};
  Sha256 exceptional_digest;
  vector<array<int, 4>> seen_label_sets;

  for (int scalar = 1; scalar <= 12; ++scalar) {
    array<int, 4> cycle_labels;
    for (int i = 0; i < 4; ++i)
      cycle_labels[i] = mod13(scalar * multipliers[i]);
    auto labels = cycle_labels;
    sort(labels.begin(), labels.end());
    if (find(seen_label_sets.begin(), seen_label_sets.end(), labels) !=
        seen_label_sets.end())
      continue;
    seen_label_sets.push_back(labels);

    vector<int> core;
    for (int r = 1; r <= 12; ++r)
      if (!binary_search(labels.begin(), labels.end(), r)) core.push_back(r);
    int core_max = *max_element(core.begin(), core.end());
    int minimum_bound = 819 * core_max / 40;
    array<vector<int>, 4> bank;
    for (int i = 0; i < 4; ++i)
      bank[i] = residue_values(cycle_labels[i], 25, 4 * minimum_bound);

    unordered_map<long long, vector<Interval>> cache;
    cache.reserve(70000);
    long long local_rows = 0;
    for (int B : bank[1]) for (int A : bank[0]) {
      if (!(2 * B <= A && A < 3 * B)) continue;
      for (int C : bank[2]) {
        if (!(C <= 2 * B && 2 * B < 3 * C)) continue;
        for (int D : bank[3]) {
          if (!(D <= 2 * C && 2 * C < 3 * D &&
                A <= 2 * D && 2 * D < 3 * A))
            continue;
          array<int, 4> speeds = {A, B, C, D};
          sort(speeds.begin(), speeds.end());
          if (speeds[0] > minimum_bound || speeds[3] > 4 * speeds[0]) continue;
          ++exceptional_rows;
          ++local_rows;
          long long key = (static_cast<long long>(speeds[0]) << 40) |
                          (static_cast<long long>(speeds[1]) << 20) | speeds[2];
          auto it = cache.find(key);
          if (it == cache.end()) {
            vector<int> eleven = core;
            eleven.insert(eleven.end(), speeds.begin(), speeds.begin() + 3);
            it = cache.emplace(key, strict_safe_components(eleven)).first;
            ++exceptional_cores;
          }
          Failure failure = first_containment_failure(it->second, speeds[3]);
          digest_row(exceptional_digest, labels, speeds, failure);
          if (!failure.failed) ++exceptional_tight;
          else if (exceptional_best_denominator == 0 ||
                   (__int128)failure.surplus * exceptional_best_denominator <
                       (__int128)exceptional_best_surplus * failure.denominator) {
            exceptional_best_surplus = failure.surplus;
            exceptional_best_denominator = failure.denominator;
            exceptional_best_row = {labels[0], labels[1], labels[2], labels[3],
                                    speeds[0], speeds[1], speeds[2], speeds[3]};
            exceptional_best_failure = failure;
          }
        }
      }
    }
    exceptional_label_histogram[labels] = local_rows;
  }

  const string expected_exceptional_digest =
      "27c45d31f19370b8b3c30e79f378b5b3ed9b1f9538062ac2f80e7dd056a6a64e";
  if (seen_label_sets.size() != 12 || exceptional_rows != 626962 ||
      exceptional_cores != 64404 || exceptional_tight != 0 ||
      exceptional_best_surplus * 550 != exceptional_best_denominator ||
      exceptional_digest.hex_digest() != expected_exceptional_digest)
    throw runtime_error("exceptional all-large audit failed");

  cout << "ALL_LARGE_EXCEPTIONAL_CHART\n";
  cout << "analytic_bound=min_replacement<=floor(819*max(retained_core)/40)<=245 "
          "max_replacement<=4*min_replacement\n";
  cout << "label_sets=" << seen_label_sets.size()
       << " exact_rows=" << exceptional_rows
       << " distinct_eleven_speed_cores=" << exceptional_cores
       << " tight_rows=" << exceptional_tight << '\n';
  cout << "label_row_histogram=";
  for (auto [labels, count] : exceptional_label_histogram) {
    cout << '(' << labels[0] << ',' << labels[1] << ',' << labels[2] << ','
         << labels[3] << "):" << count << ';';
  }
  cout << '\n';
  cout << "closest_scaled_containment_surplus=1/550 row_labels_speeds=";
  for (int x : exceptional_best_row) cout << x << ',';
  cout << " first_component=" << exceptional_best_failure.lo.n << '/'
       << exceptional_best_failure.lo.d << ',' << exceptional_best_failure.hi.n
       << '/' << exceptional_best_failure.hi.d
       << " actual_margin_surplus=1/7150\n";
  cout << "certificate_sha256=" << exceptional_digest.hex_digest() << "\n\n";

  long long anchored_rows = 0, anchored_tight = 0, anchored_cores = 0;
  map<int, long long> anchor_histogram;
  long long anchored_best_surplus = 1, anchored_best_denominator = 0;
  array<int, 8> anchored_best_row{};
  Failure anchored_best_failure{};
  Sha256 anchored_digest;

  for (int r0 = 1; r0 <= 9; ++r0)
    for (int r1 = r0 + 1; r1 <= 10; ++r1)
      for (int r2 = r1 + 1; r2 <= 11; ++r2)
        for (int r3 = r2 + 1; r3 <= 12; ++r3) {
          array<int, 4> labels = {r0, r1, r2, r3};
          vector<int> core;
          for (int r = 1; r <= 12; ++r)
            if (!binary_search(labels.begin(), labels.end(), r)) core.push_back(r);
          for (int anchor_index = 0; anchor_index < 4; ++anchor_index) {
            int x = labels[anchor_index] + 13;
            if (x > 24) continue;
            vector<int> other_labels;
            for (int i = 0; i < 4; ++i)
              if (i != anchor_index) other_labels.push_back(labels[i]);
            array<vector<int>, 3> bank;
            for (int i = 0; i < 3; ++i)
              bank[i] = residue_values(other_labels[i], x + 1, 8 * x);
            unordered_map<long long, vector<Interval>> cache;
            cache.reserve(2000);
            for (int a : bank[0]) for (int b : bank[1]) for (int c : bank[2]) {
              array<int, 3> tail = {a, b, c};
              sort(tail.begin(), tail.end());
              int v = tail[0], w = tail[1], z = tail[2];
              if (v > 2 * x || w > 2 * v || z > 2 * w) continue;
              ++anchored_rows;
              ++anchor_histogram[x];
              long long key = (static_cast<long long>(x) << 40) |
                              (static_cast<long long>(v) << 20) | w;
              auto it = cache.find(key);
              if (it == cache.end()) {
                vector<int> eleven = core;
                eleven.insert(eleven.end(), {x, v, w});
                it = cache.emplace(key, strict_safe_components(eleven)).first;
                ++anchored_cores;
              }
              array<int, 4> speeds = {x, v, w, z};
              Failure failure = first_containment_failure(it->second, z);
              digest_row(anchored_digest, labels, speeds, failure);
              if (!failure.failed) ++anchored_tight;
              else if (anchored_best_denominator == 0 ||
                       (__int128)failure.surplus * anchored_best_denominator <
                           (__int128)anchored_best_surplus * failure.denominator) {
                anchored_best_surplus = failure.surplus;
                anchored_best_denominator = failure.denominator;
                anchored_best_row = {labels[0], labels[1], labels[2], labels[3],
                                     x, v, w, z};
                anchored_best_failure = failure;
              }
            }
          }
        }

  const map<int, long long> expected_anchor_histogram = {
      {14, 4515}, {15, 5725}, {16, 7144}, {17, 8693}, {18, 10399},
      {19, 12387}, {20, 13342}, {21, 15747}, {22, 18353},
      {23, 21182}, {24, 24286}};
  const string expected_anchored_digest =
      "07594ab0e69196583fdf667b4d54c8a048a1b4d2b2a87924d26a7da4d8bc7542";
  if (anchored_rows != 141773 || anchored_cores != 38196 || anchored_tight != 0 ||
      anchor_histogram != expected_anchor_histogram ||
      anchored_best_surplus * 41 != anchored_best_denominator ||
      anchored_digest.hex_digest() != expected_anchored_digest)
    throw runtime_error("anchored doubling audit failed");

  cout << "ANCHORED_DOUBLING_CHART\n";
  cout << "symbolic_box=14<=x<=24, x<v<=2x, v<w<=2v, w<z<=2w\n";
  cout << "exact_rows=" << anchored_rows
       << " distinct_eleven_speed_cores=" << anchored_cores
       << " tight_rows=" << anchored_tight << '\n';
  cout << "anchor_row_histogram=" << map_string(anchor_histogram) << '\n';
  cout << "closest_scaled_containment_surplus=1/41 row_labels_speeds=";
  for (int x : anchored_best_row) cout << x << ',';
  cout << " first_component=" << anchored_best_failure.lo.n << '/'
       << anchored_best_failure.lo.d << ',' << anchored_best_failure.hi.n
       << '/' << anchored_best_failure.hi.d
       << " actual_margin_surplus=1/533\n";
  cout << "certificate_sha256=" << anchored_digest.hex_digest() << "\n\n";

  int height_one_rows = 0;
  pair<long long, long long> height_one_minimum = {1, 1};
  vector<array<int, 4>> height_one_minimizers;
  for (int r0 = 1; r0 <= 9; ++r0)
    for (int r1 = r0 + 1; r1 <= 10; ++r1)
      for (int r2 = r1 + 1; r2 <= 11; ++r2)
        for (int r3 = r2 + 1; r3 <= 12; ++r3) {
          array<int, 4> labels = {r0, r1, r2, r3};
          vector<int> speeds;
          for (int r = 1; r <= 12; ++r) {
            bool missing = binary_search(labels.begin(), labels.end(), r);
            speeds.push_back(missing ? r + 13 : r);
          }
          auto value = exact_maximin(speeds);
          ++height_one_rows;
          if ((__int128)value.first * height_one_minimum.second <
              (__int128)height_one_minimum.first * value.second) {
            height_one_minimum = value;
            height_one_minimizers = {labels};
          } else if ((__int128)value.first * height_one_minimum.second ==
                     (__int128)height_one_minimum.first * value.second) {
            height_one_minimizers.push_back(labels);
          }
        }
  if (height_one_rows != 495 ||
      height_one_minimum != pair<long long, long long>{1, 11} ||
      height_one_minimizers != vector<array<int, 4>>{{1, 3, 5, 7}})
    throw runtime_error("height-one audit failed");
  cout << "HEIGHT_ONE_CROSSCHECK\n";
  cout << "rows=495 exact_maximin_minimum=1/11 unique_missing_labels=(1,3,5,7)\n\n";

  tournament_audit();

  cout << "TOTAL_CERTIFIED_ROWS=" << exceptional_rows + anchored_rows
       << " TOTAL_TIGHT_ROWS=" << exceptional_tight + anchored_tight << '\n';
  cout << "PASS: every proper labelled scale-one Hamming-four lift has M>1/13.\n";
}
