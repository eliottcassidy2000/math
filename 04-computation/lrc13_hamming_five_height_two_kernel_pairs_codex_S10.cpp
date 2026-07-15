// Exact replay for THM-822: the height-{1,2} Hamming-five kernel boundary.
//
// Build and run:
//   c++ -O3 -std=c++17 \
//     04-computation/lrc13_hamming_five_height_two_kernel_pairs_codex_S10.cpp \
//     -o /tmp/lrc13_h5_h2_kernel && /tmp/lrc13_h5_h2_kernel
//
// All maximin comparisons and safe-component endpoints are exact.  The
// largest speed is 38, so signed 64-bit rational fields with __int128 cross
// products are ample.  The three faces are 1234, 1235, and 1245 in the
// increasing missing-label order.  Relation rows are ordered kernel pairs.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
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
        transform(); transformed_bits_ += 512; used_ = 0;
      }
    }
  }
  void update(const string& text) {
    update(reinterpret_cast<const uint8_t*>(text.data()), text.size());
  }
  string hex_digest() const {
    Sha256 copy = *this;
    uint64_t total_bits = copy.transformed_bits_ + 8 * copy.used_;
    copy.data_[copy.used_++] = 0x80;
    if (copy.used_ > 56) {
      while (copy.used_ < 64) copy.data_[copy.used_++] = 0;
      copy.transform(); copy.used_ = 0;
    }
    while (copy.used_ < 56) copy.data_[copy.used_++] = 0;
    for (int i = 7; i >= 0; --i)
      copy.data_[copy.used_++] = uint8_t(total_bits >> (8 * i));
    copy.transform();
    ostringstream out; out << hex << setfill('0');
    for (uint32_t word : copy.state_) out << setw(8) << word;
    return out.str();
  }
};

static void digest_u64(Sha256& digest, uint64_t value) {
  uint8_t bytes[8];
  for (int i = 7; i >= 0; --i) { bytes[i] = uint8_t(value); value >>= 8; }
  digest.update(bytes, 8);
}

struct Rat {
  long long n = 0, d = 1;
  Rat() = default;
  Rat(long long numerator, long long denominator = 1) {
    if (denominator < 0) { numerator = -numerator; denominator = -denominator; }
    long long g = gcd(numerator < 0 ? -numerator : numerator, denominator);
    n = numerator / g; d = denominator / g;
  }
  bool operator<(const Rat& other) const {
    __int128 lhs = (__int128)n * other.d;
    __int128 rhs = (__int128)other.n * d;
    return lhs < rhs;
  }
  bool operator<=(const Rat& other) const { return !(other < *this); }
  bool operator==(const Rat& other) const { return n == other.n && d == other.d; }
  bool operator!=(const Rat& other) const { return !(*this == other); }
};

struct Interval { Rat lo, hi; };

static vector<Interval> safe_bands(int speed) {
  vector<Interval> out;
  out.reserve(speed);
  for (int k = 0; k < speed; ++k)
    out.push_back({Rat(13LL * k + 1, 13LL * speed),
                   Rat(13LL * (k + 1) - 1, 13LL * speed)});
  return out;
}

static vector<Interval> intersect_unions(const vector<Interval>& a,
                                         const vector<Interval>& b) {
  vector<Interval> out;
  size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    Rat lo = a[i].lo < b[j].lo ? b[j].lo : a[i].lo;
    Rat hi = a[i].hi < b[j].hi ? a[i].hi : b[j].hi;
    if (lo < hi) out.push_back({lo, hi});
    if (a[i].hi <= b[j].hi) ++i; else ++j;
  }
  return out;
}

static vector<Interval> residual_components(vector<int> speeds,
                                            const array<vector<Interval>, 39>& bands) {
  sort(speeds.begin(), speeds.end());
  vector<Interval> current = {{Rat(0), Rat(1)}};
  for (int speed : speeds) {
    current = intersect_unions(current, bands[speed]);
    if (current.empty()) break;
  }
  return current;
}

static vector<long long> component_signature(const vector<Interval>& components) {
  vector<long long> out;
  out.reserve(4 * components.size());
  for (const auto& interval : components) {
    out.insert(out.end(), {interval.lo.n, interval.lo.d,
                           interval.hi.n, interval.hi.d});
  }
  return out;
}

struct Maximin {
  Rat value;
  vector<Rat> witnesses;
};

static Maximin exact_maximin(const vector<int>& speeds) {
  set<int> denominators;
  for (int u : speeds) denominators.insert(2 * u);
  for (size_t i = 0; i < speeds.size(); ++i)
    for (size_t j = i + 1; j < speeds.size(); ++j) {
      denominators.insert(speeds[i] + speeds[j]);
      denominators.insert(abs(speeds[i] - speeds[j]));
    }
  denominators.erase(0);
  long long best_n = 0, best_d = 1;
  set<Rat> witnesses;
  for (int d : denominators) {
    for (int n = 1; n < d; ++n) {
      int value_n = d;
      for (int u : speeds) {
        int residue = (u * n) % d;
        value_n = min(value_n, min(residue, d - residue));
      }
      __int128 comparison = (__int128)value_n * best_d - (__int128)best_n * d;
      if (comparison > 0) {
        best_n = value_n; best_d = d;
        witnesses.clear(); witnesses.insert(Rat(n, d));
      } else if (comparison == 0 && value_n != 0) {
        witnesses.insert(Rat(n, d));
      }
    }
  }
  return {Rat(best_n, best_d), vector<Rat>(witnesses.begin(), witnesses.end())};
}

static int inverse_mod13(int x) {
  x %= 13; if (x < 0) x += 13;
  for (int y = 1; y < 13; ++y) if (x * y % 13 == 1) return y;
  throw runtime_error("nonunit mod 13");
}

struct LiveData {
  uint32_t mask = 0;
  vector<int> centres;
  vector<vector<bool>> edge;
};

static LiveData live_data(const vector<int>& speeds) {
  int n = int(speeds.size());
  LiveData out;
  out.edge.assign(n, vector<bool>(n, false));
  int bit = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) if (i != j) {
      int z = (speeds[i] % 13) * inverse_mod13(speeds[j] % 13) % 13;
      bool found = false; int centre = 0;
      for (int m = -10; m <= 10; ++m) {
        long long scaled = 1LL * z * speeds[j] - 2LL * speeds[i]
                         - 13LL * m * speeds[j];
        if (-speeds[j] < scaled && scaled <= speeds[j]) {
          if (found) throw runtime_error("two handoff centres");
          found = true; centre = z - 13 * m;
        }
      }
      if (found) {
        out.mask |= uint32_t(1) << bit;
        out.edge[i][j] = true;
      }
      out.centres.push_back(found ? centre : 0);
      ++bit;
    }
  }
  return out;
}

struct Fingerprint {
  map<int, int> score;
  int triangles = 0;
  vector<int> scc;
  int hamiltonian_paths = 0;
  vector<int> first_path;
};

static vector<vector<bool>> complete_tournament(const vector<vector<bool>>& live,
                                                const vector<int>& labels,
                                                bool reverse_ties,
                                                int* silent_pairs = nullptr) {
  vector<vector<bool>> edge = live;
  int silent = 0;
  for (int i = 0; i < int(labels.size()); ++i)
    for (int j = i + 1; j < int(labels.size()); ++j) {
      if (edge[i][j] && edge[j][i]) throw runtime_error("two-cycle");
      if (!edge[i][j] && !edge[j][i]) {
        ++silent;
        bool forward = (labels[i] < labels[j]) ^ reverse_ties;
        edge[i][j] = forward; edge[j][i] = !forward;
      }
    }
  if (silent_pairs) *silent_pairs = silent;
  return edge;
}

static Fingerprint fingerprint(const vector<vector<bool>>& edge) {
  Fingerprint out; int n = int(edge.size());
  for (int i = 0; i < n; ++i) {
    int score = 0; for (int j = 0; j < n; ++j) score += edge[i][j];
    ++out.score[score];
  }
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      for (int k = j + 1; k < n; ++k)
        if ((edge[i][j] && edge[j][k] && edge[k][i]) ||
            (edge[j][i] && edge[i][k] && edge[k][j])) ++out.triangles;
  vector<vector<bool>> reach = edge;
  for (int i = 0; i < n; ++i) reach[i][i] = true;
  for (int k = 0; k < n; ++k)
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j]);
  vector<bool> used(n, false);
  for (int i = 0; i < n; ++i) if (!used[i]) {
    int size = 0;
    for (int j = 0; j < n; ++j)
      if (!used[j] && reach[i][j] && reach[j][i]) {
        used[j] = true; ++size;
      }
    out.scc.push_back(size);
  }
  sort(out.scc.rbegin(), out.scc.rend());
  vector<int> path(n); iota(path.begin(), path.end(), 0);
  do {
    bool ok = true;
    for (int i = 0; i + 1 < n; ++i) ok = ok && edge[path[i]][path[i + 1]];
    if (ok) {
      ++out.hamiltonian_paths;
      if (out.first_path.empty()) out.first_path = path;
    }
  } while (next_permutation(path.begin(), path.end()));
  return out;
}

template <size_t N>
static vector<int> heights_from_code(int code) {
  vector<int> out(N);
  for (size_t i = 0; i < N; ++i) out[i] = 1 + ((code >> (N - 1 - i)) & 1);
  return out;
}

static vector<array<int, 5>> label_sets() {
  vector<array<int, 5>> out;
  for (int a = 1; a <= 8; ++a)
    for (int b = a + 1; b <= 9; ++b)
      for (int c = b + 1; c <= 10; ++c)
        for (int d = c + 1; d <= 11; ++d)
          for (int e = d + 1; e <= 12; ++e)
            out.push_back({a, b, c, d, e});
  return out;
}

static vector<long long> prefix(const array<int, 5>& labels) {
  return vector<long long>(labels.begin(), labels.end());
}

static void append(vector<long long>& key, long long x) { key.push_back(x); }
static void append(vector<long long>& key, const vector<int>& xs) {
  key.insert(key.end(), xs.begin(), xs.end());
}
static void append(vector<long long>& key, const vector<long long>& xs) {
  key.insert(key.end(), xs.begin(), xs.end());
}

struct Group {
  long long n = 0;
  set<Rat> values;
  set<bool> tight;
};

using Groups = map<vector<long long>, Group>;

static void add_group(Groups& groups, const vector<long long>& key,
                      const Rat& value, bool tight) {
  Group& group = groups[key];
  ++group.n; group.values.insert(value); group.tight.insert(tight);
}

struct Stats {
  long long support = 0, relation_rows = 0, off_diagonal = 0;
  long long unordered_pairs = 0, collision_fibres = 0, max_fibre = 0;
  long long mixed_M_fibres = 0, rows_in_mixed_M = 0;
  long long mixed_tight_fibres = 0, rows_in_mixed_tight = 0;
};

static Stats stats(const Groups& groups, long long total) {
  Stats out; out.support = groups.size();
  for (const auto& [key, group] : groups) {
    (void)key;
    out.relation_rows += group.n * group.n;
    out.unordered_pairs += group.n * (group.n - 1) / 2;
    out.collision_fibres += group.n > 1;
    out.max_fibre = max(out.max_fibre, group.n);
    if (group.values.size() > 1) {
      ++out.mixed_M_fibres; out.rows_in_mixed_M += group.n;
    }
    if (group.tight.size() > 1) {
      ++out.mixed_tight_fibres; out.rows_in_mixed_tight += group.n;
    }
  }
  out.off_diagonal = out.relation_rows - total;
  return out;
}

static string tuple_repr(const vector<long long>& xs) {
  ostringstream out; out << '(';
  for (size_t i = 0; i < xs.size(); ++i) {
    if (i) out << ", ";
    out << xs[i];
  }
  if (xs.size() == 1) out << ',';
  out << ')'; return out.str();
}

static string row_repr(const array<int, 5>& labels, const vector<int>& heights,
                       const vector<int>& replacements) {
  return "(" + tuple_repr(vector<long long>(labels.begin(), labels.end())) + ", " +
         tuple_repr(vector<long long>(heights.begin(), heights.end())) + ", " +
         tuple_repr(vector<long long>(replacements.begin(), replacements.end())) + ")";
}

static string fraction_repr(const Rat& value) {
  return "Fraction(" + to_string(value.n) + ", " + to_string(value.d) + ")";
}

static string h2_face_key_repr(const array<int, 5>& labels, int role,
                               const vector<long long>& signature) {
  return "(" + tuple_repr(vector<long long>(labels.begin(), labels.end())) + ", " +
         to_string(role) + ", " + tuple_repr(signature) + ")";
}

static string full_h2_repr_record(const array<int, 5>& labels,
                                  const vector<int>& heights,
                                  const vector<int>& replacements,
                                  const Rat& value,
                                  const vector<long long>& signature) {
  string key = "(" + tuple_repr(vector<long long>(labels.begin(), labels.end())) +
               ", " + tuple_repr(signature) + ")";
  return "(" + row_repr(labels, heights, replacements) + ", " + key + ", " +
         fraction_repr(value) + ")";
}

static string composite_h2_repr_record(const array<int, 5>& labels,
                                       const vector<int>& heights,
                                       const vector<int>& replacements,
                                       const Rat& value,
                                       const array<vector<long long>, 3>& signatures) {
  string faces = "(";
  for (int role = 0; role < 3; ++role) {
    if (role) faces += ", ";
    faces += h2_face_key_repr(labels, role, signatures[role]);
  }
  faces += ")";
  string key = "(" + tuple_repr(vector<long long>(labels.begin(), labels.end())) +
               ", " + faces + ")";
  return "(" + row_repr(labels, heights, replacements) + ", " + key + ", " +
         fraction_repr(value) + ")";
}

struct FaceState {
  uint32_t mask = 0;
  vector<int> centres;
  vector<Interval> components;
  vector<long long> signature;
  string h0_repr;
  string h1_repr;
};

static string python_score_tuple_repr(const map<int, int>& score) {
  ostringstream out; out << '('; size_t index = 0;
  for (auto [k, v] : score) {
    if (index++) out << ", ";
    out << '(' << k << ", " << v << ')';
  }
  if (score.size() == 1) out << ',';
  out << ')'; return out.str();
}

static string python_fingerprint_repr(const Fingerprint& fp) {
  return "(" + python_score_tuple_repr(fp.score) + ", " +
         to_string(fp.triangles) + ", " +
         tuple_repr(vector<long long>(fp.scc.begin(), fp.scc.end())) + ", " +
         to_string(fp.hamiltonian_paths) + ")";
}

static string h0_face_repr(const array<int, 5>& full_labels, int role,
                           const vector<int>& face_labels, const LiveData& live) {
  auto forward = complete_tournament(live.edge, face_labels, false);
  auto reverse = complete_tournament(live.edge, face_labels, true);
  string bare = "(" + tuple_repr(vector<long long>(face_labels.begin(), face_labels.end())) +
                ", " + to_string(live.mask) + ", " +
                python_fingerprint_repr(fingerprint(forward)) + ", " +
                python_fingerprint_repr(fingerprint(reverse)) + ")";
  return "(" + tuple_repr(vector<long long>(full_labels.begin(), full_labels.end())) +
         ", " + to_string(role) + ", " + bare + ")";
}

static string composite_repr_record(const array<int, 5>& labels,
                                    const vector<int>& heights,
                                    const vector<int>& replacements,
                                    const Rat& value,
                                    const array<string, 3>& face_keys) {
  string faces = "(" + face_keys[0] + ", " + face_keys[1] + ", " + face_keys[2] + ")";
  string key = "(" + tuple_repr(vector<long long>(labels.begin(), labels.end())) +
               ", " + faces + ")";
  return "(" + row_repr(labels, heights, replacements) + ", " + key + ", " +
         fraction_repr(value) + ")";
}

static string score_string(const map<int, int>& score) {
  ostringstream out; out << '{'; bool first = true;
  for (auto [k, v] : score) { if (!first) out << ','; first = false; out << k << ':' << v; }
  out << '}'; return out.str();
}

static string vector_string(const vector<int>& xs, const vector<int>* labels = nullptr) {
  ostringstream out; out << '(';
  for (size_t i = 0; i < xs.size(); ++i) {
    if (i) out << ',';
    out << (labels ? (*labels)[xs[i]] : xs[i]);
  }
  out << ')'; return out.str();
}

int main() {
  const array<array<int, 4>, 3> faces = {{{0, 1, 2, 3},
                                          {0, 1, 2, 4},
                                          {0, 1, 3, 4}}};
  const array<string, 3> face_names = {"A1234", "B1235", "C1245"};
  array<vector<Interval>, 39> bands;
  for (int speed = 1; speed <= 38; ++speed) bands[speed] = safe_bands(speed);
  vector<array<int, 5>> labels_bank = label_sets();
  if (labels_bank.size() != 792) throw runtime_error("label-set count");

  // Face-state atlas: [label-set index][role][four-height code].
  vector<array<array<FaceState, 16>, 3>> atlas(labels_bank.size());
  array<Groups, 3> face_H0, face_H1, face_H2;
  Sha256 face_digest;
  for (size_t li = 0; li < labels_bank.size(); ++li) {
    const auto& labels = labels_bank[li];
    vector<int> core;
    for (int r = 1; r <= 12; ++r)
      if (find(labels.begin(), labels.end(), r) == labels.end()) core.push_back(r);
    for (int role = 0; role < 3; ++role) {
      vector<int> face_labels;
      for (int index : faces[role]) face_labels.push_back(labels[index]);
      for (int code = 0; code < 16; ++code) {
        vector<int> heights = heights_from_code<4>(code), replacements;
        for (int i = 0; i < 4; ++i)
          replacements.push_back(face_labels[i] + 13 * heights[i]);
        LiveData live = live_data(replacements);
        vector<int> partial = core;
        partial.insert(partial.end(), replacements.begin(), replacements.end());
        FaceState state;
        state.mask = live.mask;
        state.centres = live.centres;
        state.components = residual_components(partial, bands);
        state.signature = component_signature(state.components);
        state.h0_repr = h0_face_repr(labels, role, face_labels, live);
        state.h1_repr = "(" + state.h0_repr + ", " +
                        tuple_repr(vector<long long>(state.centres.begin(),
                                                      state.centres.end())) + ")";
        atlas[li][role][code] = state;

        vector<long long> key0 = prefix(labels);
        append(key0, role); append(key0, state.mask);
        vector<long long> key1 = key0; append(key1, state.centres);
        vector<long long> key2 = prefix(labels);
        append(key2, role); append(key2, state.signature);
        add_group(face_H0[role], key0, Rat(0), false);
        add_group(face_H1[role], key1, Rat(0), false);
        add_group(face_H2[role], key2, Rat(0), false);

        for (int x : labels) digest_u64(face_digest, x);
        digest_u64(face_digest, role);
        for (int x : heights) digest_u64(face_digest, x);
        digest_u64(face_digest, state.mask);
        for (int x : state.centres) digest_u64(face_digest, x);
        for (long long x : state.signature) digest_u64(face_digest, x);
      }
    }
  }

  Groups full_H0, full_H1, full_H2;
  array<Groups, 3> AB, ABC;
  Sha256 oracle_digest, full_h2_repr_digest, composite_h0_repr_digest,
         composite_h1_repr_digest, composite_h2_repr_digest;
  Rat global_minimum(1); vector<string> minimizers;
  long long total = 0, tight_rows = 0;

  for (size_t li = 0; li < labels_bank.size(); ++li) {
    const auto& labels = labels_bank[li];
    vector<int> core;
    for (int r = 1; r <= 12; ++r)
      if (find(labels.begin(), labels.end(), r) == labels.end()) core.push_back(r);
    for (int code = 0; code < 32; ++code) {
      vector<int> heights = heights_from_code<5>(code), replacements;
      for (int i = 0; i < 5; ++i)
        replacements.push_back(labels[i] + 13 * heights[i]);
      vector<int> packet = core;
      packet.insert(packet.end(), replacements.begin(), replacements.end());
      Maximin maximum = exact_maximin(packet);
      bool tight = maximum.value == Rat(1, 13);
      ++total; tight_rows += tight;
      string row_text = row_repr(labels, heights, replacements);
      if (maximum.value < global_minimum) {
        global_minimum = maximum.value; minimizers = {row_text};
      } else if (maximum.value == global_minimum) minimizers.push_back(row_text);

      for (int x : labels) digest_u64(oracle_digest, x);
      for (int x : heights) digest_u64(oracle_digest, x);
      digest_u64(oracle_digest, maximum.value.n);
      digest_u64(oracle_digest, maximum.value.d);
      digest_u64(oracle_digest, maximum.witnesses.size());
      for (const Rat& witness : maximum.witnesses) {
        digest_u64(oracle_digest, witness.n); digest_u64(oracle_digest, witness.d);
      }

      LiveData full_live = live_data(replacements);
      vector<Interval> full_components = residual_components(packet, bands);
      vector<long long> full_signature = component_signature(full_components);
      vector<long long> full0 = prefix(labels); append(full0, full_live.mask);
      vector<long long> full1 = full0; append(full1, full_live.centres);
      vector<long long> full2 = prefix(labels); append(full2, full_signature);
      add_group(full_H0, full0, maximum.value, tight);
      add_group(full_H1, full1, maximum.value, tight);
      add_group(full_H2, full2, maximum.value, tight);
      full_h2_repr_digest.update(full_h2_repr_record(
          labels, heights, replacements, maximum.value, full_signature));

      array<FaceState, 3> states;
      array<vector<long long>, 3> sigs;
      array<string, 3> h0_reprs, h1_reprs;
      for (int role = 0; role < 3; ++role) {
        int face_code = 0;
        for (int index : faces[role]) face_code = 2 * face_code + (heights[index] - 1);
        states[role] = atlas[li][role][face_code];
        sigs[role] = states[role].signature;
        h0_reprs[role] = states[role].h0_repr;
        h1_reprs[role] = states[role].h1_repr;
      }

      // H0 composite and AB intermediate.
      vector<long long> c0 = prefix(labels), b0 = prefix(labels);
      for (int role = 0; role < 3; ++role) append(c0, states[role].mask);
      for (int role = 0; role < 2; ++role) append(b0, states[role].mask);
      add_group(ABC[0], c0, maximum.value, tight);
      add_group(AB[0], b0, maximum.value, tight);

      // H1: retain a length marker before each ordered centre word.
      vector<long long> c1 = prefix(labels), b1 = prefix(labels);
      for (int role = 0; role < 3; ++role) {
        append(c1, states[role].mask);
        append(c1, static_cast<long long>(states[role].centres.size()));
        append(c1, states[role].centres);
      }
      for (int role = 0; role < 2; ++role) {
        append(b1, states[role].mask);
        append(b1, static_cast<long long>(states[role].centres.size()));
        append(b1, states[role].centres);
      }
      add_group(ABC[1], c1, maximum.value, tight);
      add_group(AB[1], b1, maximum.value, tight);

      // H2: literal ordered strict-safe endpoint words.
      vector<long long> c2 = prefix(labels), b2 = prefix(labels);
      for (int role = 0; role < 3; ++role) {
        append(c2, static_cast<long long>(states[role].signature.size()));
        append(c2, states[role].signature);
      }
      for (int role = 0; role < 2; ++role) {
        append(b2, static_cast<long long>(states[role].signature.size()));
        append(b2, states[role].signature);
      }
      add_group(ABC[2], c2, maximum.value, tight);
      add_group(AB[2], b2, maximum.value, tight);
      composite_h0_repr_digest.update(composite_repr_record(
          labels, heights, replacements, maximum.value, h0_reprs));
      composite_h1_repr_digest.update(composite_repr_record(
          labels, heights, replacements, maximum.value, h1_reprs));
      composite_h2_repr_digest.update(composite_h2_repr_record(
          labels, heights, replacements, maximum.value, sigs));
    }
  }

  const string expected_oracle =
      "a704bcad9cf023838f17e77d8853b03c9c98a011a92b861960f165ed08f816bd";
  const string expected_face =
      "ccee941334229f071d4b118e8cf1852c39cdc5decccc8ba4c561b0924e31a297";
  const string expected_full_h2 =
      "a750c6c1d8332278f9a71adf790f7f079a0158da117bb17a52148d30b3621dce";
  const string expected_composite_h2 =
      "9b5b2757f76e596c6ac5be91e2cc00a4c804830cd178d643e3e9bd4e722e36ec";
  const string expected_composite_h0 =
      "06186f17284850ca02ca4f032db113991609a5095f5579e8dad88aba13332ca9";
  const string expected_composite_h1 =
      "606ea18f5e8768ae38140ceebd9d2ca3acb087b8dff4d1eb102ad6b334767b2d";

  array<Stats, 3> f0, f1, f2;
  for (int role = 0; role < 3; ++role) {
    f0[role] = stats(face_H0[role], 12672);
    f1[role] = stats(face_H1[role], 12672);
    f2[role] = stats(face_H2[role], 12672);
  }
  Stats g0 = stats(full_H0, total), g1 = stats(full_H1, total),
        g2 = stats(full_H2, total);
  array<Stats, 3> ab, abc;
  for (int level = 0; level < 3; ++level) {
    ab[level] = stats(AB[level], total);
    abc[level] = stats(ABC[level], total);
  }

  if (total != 25344 || tight_rows != 0 || global_minimum != Rat(1, 12) ||
      minimizers.size() != 1 || oracle_digest.hex_digest() != expected_oracle ||
      face_digest.hex_digest() != expected_face ||
      full_h2_repr_digest.hex_digest() != expected_full_h2 ||
      composite_h0_repr_digest.hex_digest() != expected_composite_h0 ||
      composite_h1_repr_digest.hex_digest() != expected_composite_h1 ||
      composite_h2_repr_digest.hex_digest() != expected_composite_h2)
    throw runtime_error("global census or digest mismatch");
  const array<long long, 3> expected_h0_rows = {54614, 51686, 57582};
  const array<long long, 3> expected_h2_rows = {12680, 12678, 12676};
  for (int role = 0; role < 3; ++role) {
    if (f0[role].relation_rows != expected_h0_rows[role] ||
        f1[role].relation_rows != f0[role].relation_rows ||
        f1[role].support != f0[role].support ||
        f2[role].relation_rows != expected_h2_rows[role] ||
        f2[role].max_fibre != 2)
      throw runtime_error("face relation mismatch");
  }
  if (abc[0].relation_rows != 111006 || abc[0].off_diagonal != 85662 ||
      abc[0].collision_fibres != 5855 || abc[0].mixed_M_fibres != 3810 ||
      abc[0].rows_in_mixed_M != 15354 || abc[0].max_fibre != 20 ||
      abc[0].mixed_tight_fibres != 0 ||
      abc[1].relation_rows != abc[0].relation_rows ||
      abc[1].support != abc[0].support ||
      abc[2].relation_rows != 25344 || abc[2].off_diagonal != 0 ||
      g0.relation_rows != 111006 || g1.relation_rows != g0.relation_rows ||
      g2.relation_rows != 25372 || g2.off_diagonal != 28 ||
      g2.collision_fibres != 14 || g2.mixed_M_fibres != 0)
    throw runtime_error("upper relation mismatch");

  // Explicit quantitative liar pair and Tournament Analysis.
  vector<int> liar_labels = {1, 2, 3, 4, 5};
  vector<int> liar_low = {14, 15, 16, 17, 18};
  vector<int> liar_high = {27, 28, 29, 30, 31};
  LiveData low_live = live_data(liar_low), high_live = live_data(liar_high);
  if (low_live.mask != high_live.mask || low_live.centres != high_live.centres)
    throw runtime_error("liar relation mismatch");
  int silent_forward = 0, silent_reverse = 0;
  auto forward = complete_tournament(low_live.edge, liar_labels, false, &silent_forward);
  auto reverse = complete_tournament(low_live.edge, liar_labels, true, &silent_reverse);
  Fingerprint fp_forward = fingerprint(forward), fp_reverse = fingerprint(reverse);
  Maximin liar_M_low = exact_maximin(vector<int>{6,7,8,9,10,11,12,14,15,16,17,18});
  Maximin liar_M_high = exact_maximin(vector<int>{6,7,8,9,10,11,12,27,28,29,30,31});
  if (liar_M_low.value != Rat(1, 4) || liar_M_high.value != Rat(12, 37) ||
      silent_forward != 7 || silent_reverse != 7)
    throw runtime_error("liar telemetry mismatch");

  long long h0_face_total = 0, h2_face_total = 0;
  for (int role = 0; role < 3; ++role) {
    h0_face_total += f0[role].relation_rows;
    h2_face_total += f2[role].relation_rows;
  }

  cout << "LRC13 HAMMING-FIVE HEIGHT-{1,2} KERNEL PAIRS - EXACT REPLAY\n";
  cout << "theorem=THM-822 rows=C(12,5)*2^5=25344 arithmetic=integer+reduced_rational\n";
  cout << "scope=bounded_height_bank_not_arbitrary_height_closure\n\n";

  cout << "EXACT MAXIMIN CENSUS\n";
  cout << "rows=" << total << " tight_rows=" << tight_rows
       << " loose_rows=" << total << " global_minimum=1/12 unique_minimizers=1\n";
  cout << "minimum_missing=(1,3,5,7,9) heights=(1,1,1,1,1) "
          "replacements=(14,16,18,20,22) packet=2*[11]+{11}\n";
  cout << "minimum_witnesses={1,5,7,17,19,23}/24\n";
  cout << "oracle_digest=" << oracle_digest.hex_digest() << "\n\n";

  cout << "CODECS\n";
  cout << "H0=labelled_live_provider_to_owner_relation+two_silent_tie_fingerprints\n";
  cout << "H1=H0+ordered_integer_centre_on_each_live_edge\n";
  cout << "H2=context+(ordered_reduced_endpoints_of_every_strict_safe_component)\n";
  cout << "H0_equals_H1_reason=live_ratio_in_[1/2,38/14]_subset_[1/2,3) "
          "forces_k_in_{2,3,4,5,6}; labelled_z_fixes_k\n\n";

  cout << "FACE KERNEL RELATIONS unique_states_per_face=12672\n";
  for (int role = 0; role < 3; ++role)
    cout << "H0/H1 " << face_names[role] << " support=" << f0[role].support
         << " relation_rows=" << f0[role].relation_rows
         << " off_diagonal=" << f0[role].off_diagonal
         << " collision_fibres=" << f0[role].collision_fibres
         << " max_fibre=" << f0[role].max_fibre << '\n';
  for (int role = 0; role < 3; ++role)
    cout << "H2 " << face_names[role] << " support=" << f2[role].support
         << " relation_rows=" << f2[role].relation_rows
         << " off_diagonal=" << f2[role].off_diagonal
         << " collision_fibres=" << f2[role].collision_fibres
         << " max_fibre=" << f2[role].max_fibre << '\n';
  cout << "face_state_digest=" << face_digest.hex_digest() << "\n\n";

  cout << "ORDERED THREE-FACE JOIN faces=(1234,1235,1245)\n";
  cout << "literal_overlap_gluing=unique every_replacement_pair_seen_by_at_least_one_face=true\n";
  cout << "H0_AB relation_rows=" << ab[0].relation_rows
       << " off_diagonal=" << ab[0].off_diagonal
       << " collision_fibres=" << ab[0].collision_fibres
       << " mixed_M_fibres=" << ab[0].mixed_M_fibres
       << " rows_in_mixed_M=" << ab[0].rows_in_mixed_M
       << " max_fibre=" << ab[0].max_fibre << '\n';
  cout << "H0_ABC relation_rows=" << abc[0].relation_rows
       << " off_diagonal=" << abc[0].off_diagonal
       << " support=" << abc[0].support
       << " collision_fibres=" << abc[0].collision_fibres
       << " mixed_M_fibres=" << abc[0].mixed_M_fibres
       << " rows_in_mixed_M=" << abc[0].rows_in_mixed_M
       << " max_fibre=" << abc[0].max_fibre
       << " mixed_tight_fibres=" << abc[0].mixed_tight_fibres << '\n';
  cout << "H1_ABC_relation_rows=" << abc[1].relation_rows
       << " identical_partition_to_H0=true\n";
  cout << "H0_ABC_digest=" << composite_h0_repr_digest.hex_digest() << '\n';
  cout << "H1_ABC_digest=" << composite_h1_repr_digest.hex_digest() << '\n';
  cout << "H2_AB relation_rows=" << ab[2].relation_rows
       << " off_diagonal=" << ab[2].off_diagonal
       << " collision_fibres=" << ab[2].collision_fibres << '\n';
  cout << "H2_ABC relation_rows=" << abc[2].relation_rows
       << " off_diagonal=" << abc[2].off_diagonal
       << " support=" << abc[2].support
       << " injective=true\n";
  cout << "H2_ABC_digest=" << composite_h2_repr_digest.hex_digest() << "\n\n";

  cout << "FULL-ROW KERNELS\n";
  cout << "H0/H1 relation_rows=" << g0.relation_rows
       << " off_diagonal=" << g0.off_diagonal
       << " collision_fibres=" << g0.collision_fibres
       << " mixed_M_fibres=" << g0.mixed_M_fibres
       << " rows_in_mixed_M=" << g0.rows_in_mixed_M
       << " max_fibre=" << g0.max_fibre << '\n';
  cout << "H2_full support=" << g2.support
       << " relation_rows=" << g2.relation_rows
       << " off_diagonal=" << g2.off_diagonal
       << " unordered_collision_pairs=" << g2.unordered_pairs
       << " collision_fibres=" << g2.collision_fibres
       << " max_fibre=" << g2.max_fibre
       << " mixed_M_fibres=" << g2.mixed_M_fibres << '\n';
  cout << "H2_full_digest=" << full_h2_repr_digest.hex_digest() << "\n\n";

  cout << "EXPLICIT H0=H1 QUANTITATIVE LIAR\n";
  cout << "missing=(1,2,3,4,5) low_heights=(1,1,1,1,1) "
          "low_replacements=(14,15,16,17,18) M=1/4\n";
  cout << "high_heights=(2,2,2,2,2) high_replacements=(27,28,29,30,31) M=12/37\n";
  cout << "shared_live_edges_with_centres={(2->1,k=2),(3->1,k=3),(4->2,k=2)}\n";
  cout << "both_loose=true interpretation=mixed_exact_M_not_mixed_tightness\n\n";

  cout << "TOURNAMENT ANALYSIS\n";
  cout << "vertices=replacement_labels pair_observable=antisymmetric_left_handoff "
          "switch=increasing_vs_decreasing_silent_ties edge_flips=" << silent_forward << '\n';
  cout << "increasing score=" << score_string(fp_forward.score)
       << " triangles=" << fp_forward.triangles
       << " SCC=" << vector_string(fp_forward.scc)
       << " Hamiltonian_paths=" << fp_forward.hamiltonian_paths
       << " tie_Hamiltonian_path=" << vector_string(fp_forward.first_path, &liar_labels) << '\n';
  cout << "decreasing score=" << score_string(fp_reverse.score)
       << " triangles=" << fp_reverse.triangles
       << " SCC=" << vector_string(fp_reverse.scc)
       << " Hamiltonian_paths=" << fp_reverse.hamiltonian_paths
       << " tie_Hamiltonian_path=" << vector_string(fp_reverse.first_path, &liar_labels) << '\n';
  cout << "predicate_carrier=ordered_strict_safe_components_with_labelled_remaining_tooth_incidence\n\n";

  cout << "MEMORY PREFLIGHT\n";
  cout << "H0_three_face_relation_rows=" << h0_face_total
       << " packed_two_uint32_bytes=" << 8 * h0_face_total
       << " conservative_56byte_rows_plus_indexes=" << 56 * h0_face_total << '\n';
  cout << "H0_upper_relation_rows=" << abc[0].relation_rows
       << " packed_two_uint64_bytes=" << 16 * abc[0].relation_rows << '\n';
  cout << "H2_three_face_relation_rows=" << h2_face_total
       << " conservative_56byte_rows_plus_indexes=" << 56 * h2_face_total << '\n';
  cout << "verdict=canonical_join_is_small;H0_is_quantitatively_impure;H2_is_finitely_static_injective\n\n";

  cout << "SCOPE GUARDRAILS\n";
  cout << "all_codecs_tight_pure_only_because_bank_has_zero_tight_rows\n";
  cout << "H2_face_injectivity_is_finite_for_heights_{1,2}_not_an_all_height_Markov_theorem\n";
  cout << "PASS: all 25344 rows loose; exact relation and digest audits matched.\n";
}
