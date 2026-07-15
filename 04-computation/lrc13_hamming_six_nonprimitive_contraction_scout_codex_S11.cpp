// Exact scout for the nonprimitive scale-one Hamming-six branch at q=13.
//
// A scale-one Hamming-six packet has the form
//
//   W = ([12] \ R) union {r + 13 h_r : r in R},  |R|=6, h_r>=1.
//
// If gcd(W)>1, the retained six labels force gcd(W)=2 and
// R={1,3,5,7,9,11}.  Writing h_{2i-1}=2k_i+1 and dividing by two gives
//
//   {1,...,6} union {6+i+13k_i : i=1,...,6}.
//
// The cases with zero, one through five, and six positive k_i are respectively
// [12], proper scale-one Hamming radii <=5, and the one top-half H6 chamber.
// Existing exact results close radii <=5.  This file closes the top-half
// chamber by the same longest-safe-component discrepancy used at H4/H5.
// Thus it certifies the new, deliberately narrow conclusion:
//
//   among nonprimitive scale-one H6 packets, only 2[12] can be tight.
//
// It does NOT close primitive H6 packets or the global n=12 sporadic branch.
//
// Kakeya/assumption-challenge sidecar.  The useful "needles" are the six
// periodic danger combs D_7,...,D_12 acting on the twelve components of the
// strict-safe set E({1,...,6}), not six runner vertices in a tournament.  The
// exact component--comb incidence graph preserves union coverage.  Pairwise
// comb overlap is symmetric, so switching it into a tournament would require
// an arbitrary gauge and would destroy the higher-order covering predicate.
// The replay records the four zero-overlap-debt components of the AP equality.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

struct Rat {
  long long n = 0;
  long long d = 1;

  Rat(long long a = 0, long long b = 1) {
    if (b == 0)
      throw runtime_error("zero rational denominator");
    if (b < 0) {
      a = -a;
      b = -b;
    }
    const long long g = gcd(a < 0 ? -a : a, b);
    n = a / g;
    d = b / g;
  }
};

bool operator<(Rat a, Rat b) {
  return (__int128)a.n * b.d < (__int128)b.n * a.d;
}
bool operator<=(Rat a, Rat b) { return !(b < a); }
bool operator==(Rat a, Rat b) { return a.n == b.n && a.d == b.d; }
bool operator!=(Rat a, Rat b) { return !(a == b); }

Rat add_signed(Rat a, Rat b, int sign) {
  const long long g = gcd(a.d, b.d);
  const __int128 nn = (__int128)a.n * (b.d / g) +
                      (__int128)sign * b.n * (a.d / g);
  const __int128 dd = (__int128)(a.d / g) * b.d;
  if (nn > INT64_MAX || nn < INT64_MIN || dd > INT64_MAX)
    throw runtime_error("rational arithmetic overflow");
  return Rat((long long)nn, (long long)dd);
}
Rat operator+(Rat a, Rat b) { return add_signed(a, b, +1); }
Rat operator-(Rat a, Rat b) { return add_signed(a, b, -1); }

ostream &operator<<(ostream &out, Rat q) { return out << q.n << "/" << q.d; }

struct Interval {
  Rat lo;
  Rat hi;
};

vector<Interval> safe_bands(int u) {
  vector<Interval> out;
  out.reserve(u);
  for (int k = 0; k < u; ++k) {
    out.push_back({Rat(13LL * k + 1, 13LL * u),
                   Rat(13LL * (k + 1) - 1, 13LL * u)});
  }
  return out;
}

vector<Interval> danger_bands(int u) {
  vector<Interval> out;
  out.reserve(u + 1);
  for (int k = 0; k <= u; ++k) {
    Rat lo(13LL * k - 1, 13LL * u);
    Rat hi(13LL * k + 1, 13LL * u);
    if (lo < Rat(0))
      lo = Rat(0);
    if (Rat(1) < hi)
      hi = Rat(1);
    if (lo < hi)
      out.push_back({lo, hi});
  }
  return out;
}

vector<Interval> meet(vector<Interval> const &a, vector<Interval> const &b) {
  vector<Interval> out;
  out.reserve(a.size() + 8);
  size_t i = 0;
  size_t j = 0;
  while (i < a.size() && j < b.size()) {
    const Rat lo = a[i].lo < b[j].lo ? b[j].lo : a[i].lo;
    const Rat hi = a[i].hi < b[j].hi ? a[i].hi : b[j].hi;
    if (lo < hi)
      out.push_back({lo, hi});
    if (a[i].hi <= b[j].hi)
      ++i;
    else
      ++j;
  }
  return out;
}

vector<Interval> merged(vector<Interval> v) {
  sort(v.begin(), v.end(), [](Interval const &a, Interval const &b) {
    if (a.lo == b.lo)
      return a.hi < b.hi;
    return a.lo < b.lo;
  });
  vector<Interval> out;
  for (Interval x : v) {
    if (!out.empty() && x.lo <= out.back().hi) {
      if (out.back().hi < x.hi)
        out.back().hi = x.hi;
    } else {
      out.push_back(x);
    }
  }
  return out;
}

Rat measure(vector<Interval> const &v) {
  Rat out(0);
  for (Interval x : v)
    out = out + (x.hi - x.lo);
  return out;
}

Rat longest_length(vector<Interval> const &v) {
  if (v.empty())
    throw runtime_error("longest component requested from empty safe set");
  Rat best = v.front().hi - v.front().lo;
  for (Interval x : v) {
    const Rat length = x.hi - x.lo;
    if (best < length)
      best = length;
  }
  return best;
}

vector<Interval> safe_by_closed_danger_union(vector<int> const &speeds) {
  vector<Interval> danger;
  for (int u : speeds) {
    const vector<Interval> teeth = danger_bands(u);
    danger.insert(danger.end(), teeth.begin(), teeth.end());
  }
  const vector<Interval> covered = merged(danger);
  vector<Interval> gaps;
  Rat cursor(0);
  for (Interval x : covered) {
    if (cursor < x.lo)
      gaps.push_back({cursor, x.lo});
    if (cursor < x.hi)
      cursor = x.hi;
  }
  if (cursor < Rat(1))
    gaps.push_back({cursor, Rat(1)});
  return gaps;
}

bool same_intervals(vector<Interval> const &a, vector<Interval> const &b) {
  if (a.size() != b.size())
    return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i].lo != b[i].lo || a[i].hi != b[i].hi)
      return false;
  return true;
}

// If m periodic danger combs cover a strict-safe component of length L, the
// sharp one-period discrepancy gives
//
//   min remaining speed <= floor(22 m / (13 (13-2m) L)).
//
// The denominator stays positive exactly through m=6; m=7 is the known wall.
long long discrepancy_cap(int m, Rat length) {
  if (m < 1 || m > 6 || length.n <= 0)
    throw runtime_error("invalid discrepancy-cap input");
  return (long long)((__int128)22 * m * length.d /
                     ((__int128)13 * (13 - 2 * m) * length.n));
}

unordered_map<int, vector<Interval>> band_cache;

vector<Interval> const &bands_for(int u) {
  const auto found = band_cache.find(u);
  if (found != band_cache.end())
    return found->second;
  return band_cache.emplace(u, safe_bands(u)).first->second;
}

struct Trace128 {
  uint64_t a = 0xcbf29ce484222325ULL;
  uint64_t b = 0x9e3779b97f4a7c15ULL;

  void mix(uint64_t x) {
    a ^= x;
    a *= 0x100000001b3ULL;
    b ^= x + 0x9e3779b97f4a7c15ULL + (b << 6) + (b >> 2);
  }
};

array<unsigned long long, 7> node_count{};
array<unsigned long long, 6> edge_count{};
array<unsigned long long, 6> dead_count{};
array<unsigned long long, 7> covering_count{};
array<unsigned long long, 7> nonempty_count{};
array<long long, 6> max_cap{};
array<Rat, 7> minimum_longest{};
array<int, 6> owner{};
array<int, 6> speed{};
Trace128 trace_digest;

struct DeadLeaf {
  array<int, 6> owners{};
  array<int, 6> speeds{};
  vector<Interval> components;
  long long cap = 0;
};
vector<DeadLeaf> deepest_dead_leaves;

void trace_state(int depth, vector<Interval> const &components, long long cap) {
  trace_digest.mix((uint64_t)depth);
  trace_digest.mix((uint64_t)cap);
  trace_digest.mix((uint64_t)components.size());
  for (int i = 0; i < depth; ++i) {
    trace_digest.mix((uint64_t)owner[i]);
    trace_digest.mix((uint64_t)speed[i]);
  }
  for (Interval x : components) {
    trace_digest.mix((uint64_t)x.lo.n);
    trace_digest.mix((uint64_t)x.lo.d);
    trace_digest.mix((uint64_t)x.hi.n);
    trace_digest.mix((uint64_t)x.hi.d);
  }
}

constexpr array<int, 6> TOP_LABELS = {7, 8, 9, 10, 11, 12};

void recurse_top_half(vector<Interval> const &components, int depth) {
  ++node_count[depth];
  if (components.empty()) {
    ++covering_count[depth];
    trace_state(depth, components, -2);
    return;
  }

  ++nonempty_count[depth];
  const Rat length = longest_length(components);
  if (nonempty_count[depth] == 1 || length < minimum_longest[depth])
    minimum_longest[depth] = length;

  if (depth == 6) {
    trace_state(depth, components, -1);
    return;
  }

  const int remaining = 6 - depth;
  const long long cap = discrepancy_cap(remaining, length);
  max_cap[depth] = max(max_cap[depth], cap);
  trace_state(depth, components, cap);

  const int previous = depth == 0 ? 0 : speed[depth - 1];
  const unsigned long long edges_before = edge_count[depth];
  for (int label : TOP_LABELS) {
    bool used = false;
    for (int i = 0; i < depth; ++i)
      if (owner[i] == label)
        used = true;
    if (used)
      continue;

    int u = label + 13;
    if (u <= previous)
      u += 13 * ((previous - u) / 13 + 1);
    for (; u <= cap; u += 13) {
      ++edge_count[depth];
      owner[depth] = label;
      speed[depth] = u;
      recurse_top_half(meet(components, bands_for(u)), depth + 1);
    }
  }
  if (edge_count[depth] == edges_before) {
    ++dead_count[depth];
    if (depth == 5)
      deepest_dead_leaves.push_back({owner, speed, components, cap});
  }
}

string label_set(vector<int> const &labels) {
  string out = "{";
  for (size_t i = 0; i < labels.size(); ++i) {
    if (i)
      out += ",";
    out += to_string(labels[i]);
  }
  return out + "}";
}

void gcd_core_census() {
  int rows = 0;
  vector<vector<int>> candidates;
  vector<int> candidate_gcd;
  for (int mask = 0; mask < (1 << 12); ++mask) {
    if (__builtin_popcount((unsigned)mask) != 6)
      continue;
    ++rows;
    int retained_gcd = 0;
    vector<int> missing;
    for (int r = 1; r <= 12; ++r) {
      if (mask & (1 << (r - 1)))
        missing.push_back(r);
      else
        retained_gcd = gcd(retained_gcd, r);
    }
    if (retained_gcd > 1) {
      candidates.push_back(missing);
      candidate_gcd.push_back(retained_gcd);
    }
  }
  const vector<int> odds = {1, 3, 5, 7, 9, 11};
  if (rows != 924 || candidates.size() != 1 || candidates[0] != odds ||
      candidate_gcd[0] != 2)
    throw runtime_error("gcd-core census mismatch");

  cout << "gcd_core_census rows=" << rows
       << " candidate_rows=" << candidates.size()
       << " unique_missing=" << label_set(candidates[0])
       << " retained_gcd=" << candidate_gcd[0] << "\n";
  cout << "nonprimitive_criterion=missing_odds_and_all_six_heights_odd"
       << " gcd_exactly=2\n";
  cout << "divide_by_two={1,...,6}_union_{6+i+13k_i:i=1,...,6}"
       << " branch_by_positive_k=0:AP,1..5:H1..H5,6:top_half_H6\n";
}

void kakeya_needle_ledger(vector<Interval> const &root) {
  array<vector<Interval>, 13> danger;
  for (int u = 7; u <= 12; ++u)
    danger[u] = danger_bands(u);

  Rat total_measure = measure(root);
  Rat longest = longest_length(root);
  int covered = 0;
  int zero_debt = 0;
  int incidence_edges = 0;
  array<int, 7> component_degree_hist{};
  array<int, 13> comb_degrees{};
  vector<string> zero_records;

  for (Interval component : root) {
    vector<Interval> all_pieces;
    Rat sum_piece_measure(0);
    vector<int> touching;
    vector<int> full_owners;
    for (int u = 7; u <= 12; ++u) {
      const vector<Interval> pieces =
          meet(vector<Interval>{component}, danger[u]);
      const Rat piece_measure = measure(pieces);
      if (Rat(0) < piece_measure) {
        touching.push_back(u);
        ++comb_degrees[u];
        ++incidence_edges;
        all_pieces.insert(all_pieces.end(), pieces.begin(), pieces.end());
      }
      if (piece_measure == component.hi - component.lo)
        full_owners.push_back(u);
      sum_piece_measure = sum_piece_measure + piece_measure;
    }

    const vector<Interval> union_pieces = merged(all_pieces);
    const Rat union_measure = measure(union_pieces);
    const Rat component_measure = component.hi - component.lo;
    if (union_measure == component_measure)
      ++covered;
    else
      throw runtime_error("AP needle union failed to cover a core component");
    const Rat debt = sum_piece_measure - union_measure;
    ++component_degree_hist[touching.size()];
    if (debt == Rat(0)) {
      ++zero_debt;
      if (touching.size() != 1 || full_owners != touching)
        throw runtime_error("zero-debt component lacks a unique full owner");
      zero_records.push_back("[" + to_string(component.lo.n) + "/" +
                             to_string(component.lo.d) + "," +
                             to_string(component.hi.n) + "/" +
                             to_string(component.hi.d) + "]:" +
                             to_string(touching[0]));
    }
  }

  if (root.size() != 12 || total_measure != Rat(27, 65) ||
      longest != Rat(1, 13) || covered != 12 || zero_debt != 4 ||
      incidence_edges != 34)
    throw runtime_error("Kakeya needle ledger mismatch");

  cout << "kakeya_core={1,...,6} needles={7,...,12} components="
       << root.size() << " measure=" << total_measure
       << " longest=" << longest << " covered=" << covered
       << " zero_overlap_debt=" << zero_debt
       << " incidence_edges=" << incidence_edges << "\n";
  cout << "component_degree_hist=1:4,2:2,3:2,4:2,6:2"
       << " comb_degrees=7:6,8:4,9:6,10:4,11:10,12:4\n";
  cout << "zero_debt_unique_owners=";
  for (size_t i = 0; i < zero_records.size(); ++i) {
    if (i)
      cout << ",";
    cout << zero_records[i];
  }
  cout << "\n";
  cout << "tournament_guardrail=pairwise_comb_overlap_is_symmetric;"
       << " faithful_carrier=component_comb_incidence_with_lengths\n";
}

void closed_danger_leaf_crosscheck() {
  if (deepest_dead_leaves.size() != 2)
    throw runtime_error("unexpected deepest dead-leaf count");
  constexpr int EXPECTED_OWNERS[2][5] = {
      {7, 10, 9, 8, 11}, {7, 10, 9, 8, 12}};
  constexpr int EXPECTED_SPEEDS[2][5] = {
      {33, 36, 48, 60, 63}, {33, 36, 48, 60, 64}};
  constexpr int EXPECTED_UNUSED[2] = {12, 11};
  constexpr int EXPECTED_LEAST_NEXT[2] = {64, 76};
  constexpr size_t EXPECTED_COMPONENTS[2] = {58, 54};
  const Rat expected_longest[2] = {Rat(17, 3120), Rat(47, 8580)};
  int matches = 0;
  for (size_t index = 0; index < deepest_dead_leaves.size(); ++index) {
    const DeadLeaf &leaf = deepest_dead_leaves[index];
    for (int i = 0; i < 5; ++i) {
      if (leaf.owners[i] != EXPECTED_OWNERS[index][i] ||
          leaf.speeds[i] != EXPECTED_SPEEDS[index][i])
        throw runtime_error("deepest dead-leaf path mismatch");
    }
    if (leaf.cap != 28 || leaf.components.size() != EXPECTED_COMPONENTS[index] ||
        longest_length(leaf.components) != expected_longest[index])
      throw runtime_error("deepest dead-leaf metric mismatch");
    vector<int> packet = {1, 2, 3, 4, 5, 6};
    for (int i = 0; i < 5; ++i)
      packet.push_back(leaf.speeds[i]);
    const vector<Interval> independent = safe_by_closed_danger_union(packet);
    if (!same_intervals(independent, leaf.components))
      throw runtime_error("closed-danger leaf reconstruction mismatch");
    ++matches;

    int unused = -1;
    for (int label : TOP_LABELS) {
      bool used = false;
      for (int i = 0; i < 5; ++i)
        if (leaf.owners[i] == label)
          used = true;
      if (!used)
        unused = label;
    }
    int least_next = unused + 13;
    while (least_next <= leaf.speeds[4])
      least_next += 13;
    if (unused != EXPECTED_UNUSED[index] ||
        least_next != EXPECTED_LEAST_NEXT[index] || least_next <= leaf.cap)
      throw runtime_error("deepest leaf is not arithmetically dead");

    cout << "deepest_dead_leaf=" << index + 1 << " path=";
    for (int i = 0; i < 5; ++i) {
      if (i)
        cout << ",";
      cout << leaf.owners[i] << ":" << leaf.speeds[i];
    }
    cout << " unused=" << unused << " least_next=" << least_next
         << " cap=" << leaf.cap
         << " components=" << leaf.components.size()
         << " longest=" << longest_length(leaf.components) << "\n";
  }
  cout << "closed_danger_deepest_leaf_checks=" << deepest_dead_leaves.size()
       << " exact_component_matches=" << matches << "\n";
}

void assert_frozen_recursion() {
  constexpr array<unsigned long long, 7> EXPECTED_NODES = {
      1, 54, 3612, 130515, 2104, 2, 0};
  constexpr array<unsigned long long, 6> EXPECTED_EDGES = {
      54, 3612, 130515, 2104, 2, 0};
  constexpr array<unsigned long long, 6> EXPECTED_DEAD = {
      0, 0, 0, 129772, 2103, 2};
  constexpr array<unsigned long long, 7> EXPECTED_COVERING = {
      0, 0, 0, 0, 0, 0, 0};
  constexpr array<long long, 6> EXPECTED_MAX_CAP = {
      132, 430, 683, 608, 315, 28};
  constexpr uint64_t EXPECTED_TRACE_A = 0x919c6848d4e1187aULL;
  constexpr uint64_t EXPECTED_TRACE_B = 0x2cef093e58982ae6ULL;

  if (node_count != EXPECTED_NODES || edge_count != EXPECTED_EDGES ||
      dead_count != EXPECTED_DEAD || covering_count != EXPECTED_COVERING ||
      max_cap != EXPECTED_MAX_CAP || trace_digest.a != EXPECTED_TRACE_A ||
      trace_digest.b != EXPECTED_TRACE_B || band_cache.size() != 313)
    throw runtime_error("frozen top-half recursion certificate mismatch");
}

int main() {
  cout << "LRC13 SCALE-ONE HAMMING-SIX NONPRIMITIVE CONTRACTION SCOUT\n";
  cout << "arithmetic=exact_int64_rationals comparisons=int128"
       << " endpoint_convention=strict_safe_open\n";

  gcd_core_census();

  vector<Interval> root = {{Rat(0), Rat(1)}};
  for (int u = 1; u <= 6; ++u)
    root = meet(root, bands_for(u));
  kakeya_needle_ledger(root);

  for (Rat &q : minimum_longest)
    q = Rat(1);
  recurse_top_half(root, 0);

  cout << "top_half_packet={1,...,6}_union_{r+13h_r:r=7,...,12,h_r>=1}\n";
  for (int depth = 0; depth <= 6; ++depth) {
    cout << "depth=" << depth << " nodes=" << node_count[depth];
    if (depth < 6) {
      cout << " edges=" << edge_count[depth]
           << " dead=" << dead_count[depth]
           << " max_cap=" << max_cap[depth];
    }
    cout << " covering=" << covering_count[depth] << " min_longest=";
    if (nonempty_count[depth])
      cout << minimum_longest[depth];
    else
      cout << "-";
    cout << "\n";
  }

  closed_danger_leaf_crosscheck();
  assert_frozen_recursion();

  cout << "trace128=" << hex << setfill('0') << setw(16) << trace_digest.a
       << setw(16) << trace_digest.b << dec
       << " cached_speeds=" << band_cache.size() << "\n";

  cout << "PASS: top-half H6 chamber has no covering completion\n";
  cout << "CONCLUSION: among nonprimitive scale-one H6 packets, only 2[12]"
       << " can be tight\n";
  cout << "SCOPE: primitive scale-one H6 and global n=12 sporadic emptiness"
       << " remain open\n";
  return 0;
}
