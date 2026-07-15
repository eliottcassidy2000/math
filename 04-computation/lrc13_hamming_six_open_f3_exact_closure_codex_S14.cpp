// Exact closure of the six product-1/16 primitive-core Hamming-six rows.
//
// THM-815 C.2 partitions the twenty scale-one H6 rows made from three full
// antipodal pairs into six germ-cycle closures, eight forced-coordinate slice
// closures, and the six rows below.  Their one-step handoff graphs contain a
// contracting cycle, so that quotient cannot close them.  The proof-bearing
// carrier here is the exact residual strict-safe interval union together with
// the labelled bank of unplaced danger combs.
//
// At every prefix, order the six replacement speeds numerically.  If k<=6
// unplaced danger combs cover a strict-safe component of length L, the sharp
// one-component discrepancy gives
//
//   next speed <= floor(22 k / (13 (13-2k) L)).
//
// Enumerating every unused residue class and every proper congruent lift up to
// this cap includes every hypothetical tight completion exactly once.  All
// arithmetic and endpoint comparisons are exact.  Frozen row traces and depth
// counts certify the trees; every deepest leaf is independently reconstructed
// as the complement of the full closed-danger union.
//
// Tournament analysis: numerical speed comparison is a transitive enumeration
// gauge (score sequence 0,...,5, singleton SCCs, one Hamiltonian path), not a
// quotient of union coverage.  The earlier germ carrier is a weighted digraph,
// not a tournament.  Vertices must therefore be allowed to be residual safe
// components / proof obligations rather than runners alone.

#include <algorithm>
#include <array>
#include <cstdint>
#include <future>
#include <iomanip>
#include <iostream>
#include <memory>
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
  for (int k = 0; k < u; ++k)
    out.push_back({Rat(13LL * k + 1, 13LL * u),
                   Rat(13LL * (k + 1) - 1, 13LL * u)});
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

vector<Interval> safe_by_closed_danger_union(vector<int> const &speeds) {
  vector<Interval> danger;
  for (int u : speeds) {
    vector<Interval> teeth = danger_bands(u);
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

long long discrepancy_cap(int remaining, Rat length) {
  if (remaining < 1 || remaining > 6 || length.n <= 0)
    throw runtime_error("invalid discrepancy-cap input");
  return (long long)((__int128)22 * remaining * length.d /
                     ((__int128)13 * (13 - 2 * remaining) * length.n));
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

struct RowSpec {
  array<int, 6> missing;
  array<int, 6> core;
  array<unsigned long long, 7> nodes;
  array<unsigned long long, 6> edges;
  array<unsigned long long, 6> dead;
  array<long long, 6> max_cap;
  array<Rat, 6> min_longest;
  uint64_t trace_a;
  uint64_t trace_b;
  size_t cached_speeds;
};

vector<RowSpec> row_specs() {
  return {
      {{1, 2, 4, 9, 11, 12},
       {3, 5, 6, 7, 8, 10},
       {1, 64, 4769, 195705, 7340, 50, 0},
       {64, 4769, 195705, 7340, 50, 0},
       {0, 0, 0, 193060, 7305, 50},
       {152, 506, 809, 692, 337, 45},
       {Rat(1, 15), Rat(11, 1976), Rat(1, 598), Rat(11, 10504),
        Rat(11, 9880), Rat(4, 1189)},
       0x063c1b3f17d4520fULL,
       0xeaf265b370a1c09dULL,
       371},
      {{1, 2, 6, 7, 11, 12},
       {3, 4, 5, 8, 9, 10},
       {1, 64, 4761, 195502, 5875, 12, 0},
       {64, 4761, 195502, 5875, 12, 0},
       {0, 0, 0, 193257, 5865, 12},
       {152, 500, 800, 681, 293, 42},
       {Rat(1, 15), Rat(11, 1950), Rat(11, 6500), Rat(11, 10335),
        Rat(1, 780), Rat(41, 11284)},
       0x3b468adfcff92e6bULL,
       0xafc5cf50ea6b5000ULL,
       368},
      {{1, 3, 6, 7, 10, 12},
       {2, 4, 5, 8, 9, 11},
       {1, 97, 10343, 620068, 17195, 39, 0},
       {97, 10343, 620068, 17195, 39, 0},
       {0, 0, 0, 613727, 17170, 39},
       {223, 740, 1184, 1011, 488, 33},
       {Rat(1, 22), Rat(11, 2886), Rat(11, 9620), Rat(11, 15340),
        Rat(1, 1300), Rat(23, 4940)},
       0x4844a2fbb06e7c9aULL,
       0xad37d047eefeecabULL,
       544},
      {{2, 4, 5, 8, 9, 11},
       {1, 3, 6, 7, 10, 12},
       {1, 93, 9560, 550797, 10885, 21, 0},
       {93, 9560, 550797, 10885, 21, 0},
       {0, 0, 0, 547550, 10870, 21},
       {213, 710, 1136, 972, 444, 35},
       {Rat(1, 21), Rat(11, 2769), Rat(11, 9230), Rat(11, 14755),
        Rat(11, 13000), Rat(1, 228)},
       0xa60c4579e7f6cff7ULL,
       0x404360e8db9bc6f4ULL,
       524},
      {{3, 4, 5, 8, 9, 10},
       {1, 2, 6, 7, 11, 12},
       {1, 79, 7032, 349248, 10523, 10, 0},
       {79, 7032, 349248, 10523, 10, 0},
       {0, 0, 0, 345927, 10515, 10},
       {185, 616, 985, 844, 382, 53},
       {Rat(47, 858), Rat(11, 2405), Rat(1, 728), Rat(11, 12805),
        Rat(11, 11180), Rat(3, 1034)},
       0x144d8b3b67d164cdULL,
       0xb04afbf0bafffa3cULL,
       453},
      {{3, 5, 6, 7, 8, 10},
       {1, 2, 4, 9, 11, 12},
       {1, 97, 10348, 620350, 22597, 70, 0},
       {97, 10348, 620350, 22597, 70, 0},
       {0, 0, 0, 611134, 22560, 70},
       {225, 746, 1193, 1022, 521, 111},
       {Rat(58, 1287), Rat(11, 2912), Rat(11, 9698), Rat(11, 15509),
        Rat(4222, 5853081), Rat(28, 20301)},
       0x242306835e1bb0c1ULL,
       0x289ce8d030c29a7aULL,
       548},
  };
}

struct Replay {
  RowSpec const &spec;
  array<unsigned long long, 7> nodes{};
  array<unsigned long long, 6> edges{};
  array<unsigned long long, 6> dead{};
  array<unsigned long long, 7> covering{};
  array<unsigned long long, 7> nonempty{};
  array<long long, 6> max_cap{};
  array<Rat, 7> min_longest{};
  array<int, 6> owner{};
  array<int, 6> speed{};
  unordered_map<int, vector<Interval>> cache;
  Trace128 trace;
  unsigned long long root_crosschecks = 0;
  unsigned long long deepest_crosschecks = 0;

  explicit Replay(RowSpec const &row) : spec(row) {
    for (Rat &q : min_longest)
      q = Rat(1);
  }

  vector<Interval> const &bands_for(int u) {
    const auto found = cache.find(u);
    if (found != cache.end())
      return found->second;
    return cache.emplace(u, safe_bands(u)).first->second;
  }

  void trace_state(int depth, vector<Interval> const &components, long long cap) {
    trace.mix((uint64_t)depth);
    trace.mix((uint64_t)cap);
    trace.mix((uint64_t)components.size());
    for (int i = 0; i < depth; ++i) {
      trace.mix((uint64_t)owner[i]);
      trace.mix((uint64_t)speed[i]);
    }
    for (Interval x : components) {
      trace.mix((uint64_t)x.lo.n);
      trace.mix((uint64_t)x.lo.d);
      trace.mix((uint64_t)x.hi.n);
      trace.mix((uint64_t)x.hi.d);
    }
  }

  void crosscheck(vector<Interval> const &components, int depth) {
    vector<int> packet(spec.core.begin(), spec.core.end());
    for (int i = 0; i < depth; ++i)
      packet.push_back(speed[i]);
    if (!same_intervals(components, safe_by_closed_danger_union(packet)))
      throw runtime_error("independent closed-danger reconstruction mismatch");
  }

  void recurse(vector<Interval> const &components, int depth) {
    ++nodes[depth];
    if (components.empty()) {
      ++covering[depth];
      trace_state(depth, components, -2);
      return;
    }

    ++nonempty[depth];
    const Rat length = longest_length(components);
    if (nonempty[depth] == 1 || length < min_longest[depth])
      min_longest[depth] = length;

    if (depth == 6) {
      trace_state(depth, components, -1);
      return;
    }

    const int remaining = 6 - depth;
    const long long cap = discrepancy_cap(remaining, length);
    max_cap[depth] = max(max_cap[depth], cap);
    trace_state(depth, components, cap);

    const int previous = depth == 0 ? 0 : speed[depth - 1];
    const unsigned long long before = edges[depth];
    for (int label : spec.missing) {
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
        ++edges[depth];
        owner[depth] = label;
        speed[depth] = u;
        recurse(meet(components, bands_for(u)), depth + 1);
      }
    }
    if (edges[depth] == before) {
      ++dead[depth];
      if (depth == 5) {
        crosscheck(components, depth);
        ++deepest_crosschecks;
      }
    }
  }

  void run() {
    vector<Interval> root = {{Rat(0), Rat(1)}};
    for (int u : spec.core)
      root = meet(root, bands_for(u));
    crosscheck(root, 0);
    ++root_crosschecks;
    recurse(root, 0);
  }

  void assert_frozen() const {
    const array<unsigned long long, 7> zero_cover{};
    if (nodes != spec.nodes || edges != spec.edges || dead != spec.dead ||
        covering != zero_cover || max_cap != spec.max_cap ||
        trace.a != spec.trace_a || trace.b != spec.trace_b ||
        cache.size() != spec.cached_speeds || root_crosschecks != 1 ||
        deepest_crosschecks != spec.nodes[5])
      throw runtime_error("frozen row certificate mismatch");
    for (int d = 0; d < 6; ++d)
      if (min_longest[d] != spec.min_longest[d])
        throw runtime_error("frozen minimum-longest mismatch");
  }
};

string set_text(array<int, 6> const &values) {
  string out = "{";
  for (int i = 0; i < 6; ++i) {
    if (i)
      out += ",";
    out += to_string(values[i]);
  }
  return out + "}";
}

string counts_text(array<unsigned long long, 7> const &values) {
  string out;
  for (int i = 0; i < 7; ++i) {
    if (i)
      out += ",";
    out += to_string(values[i]);
  }
  return out;
}

string counts_text(array<unsigned long long, 6> const &values) {
  string out;
  for (int i = 0; i < 6; ++i) {
    if (i)
      out += ",";
    out += to_string(values[i]);
  }
  return out;
}

int main() {
  cout << "LRC13 SCALE-ONE HAMMING-SIX SIX OPEN F=3 ROW CLOSURE\n";
  cout << "arithmetic=exact_int64_rationals comparisons=int128"
       << " endpoint_convention=strict_safe_open\n";
  cout << "method=numerically_ordered_rowwise_longest_component_recursion\n";

  const vector<RowSpec> specs = row_specs();
  array<unsigned long long, 7> aggregate_nodes{};
  array<unsigned long long, 6> aggregate_edges{};
  unsigned long long total_crosschecks = 0;

  vector<future<unique_ptr<Replay>>> jobs;
  for (RowSpec const &spec : specs) {
    RowSpec const *row = &spec;
    jobs.push_back(async(launch::async, [row]() {
      auto replay = make_unique<Replay>(*row);
      replay->run();
      replay->assert_frozen();
      return replay;
    }));
  }

  for (future<unique_ptr<Replay>> &job : jobs) {
    unique_ptr<Replay> replay_holder = job.get();
    Replay const &replay = *replay_holder;
    RowSpec const &spec = replay.spec;
    for (int d = 0; d < 7; ++d)
      aggregate_nodes[d] += replay.nodes[d];
    for (int d = 0; d < 6; ++d)
      aggregate_edges[d] += replay.edges[d];
    total_crosschecks += replay.root_crosschecks + replay.deepest_crosschecks;

    cout << "row=" << set_text(spec.missing)
         << " core=" << set_text(spec.core)
         << " nodes=" << counts_text(replay.nodes)
         << " edges=" << counts_text(replay.edges)
         << " deepest_crosschecks=" << replay.deepest_crosschecks
         << " trace128=" << hex << setfill('0') << setw(16) << replay.trace.a
         << setw(16) << replay.trace.b << dec
         << " cached_speeds=" << replay.cache.size() << "\n";
  }

  const array<unsigned long long, 7> expected_nodes = {
      6, 494, 46813, 2531670, 74415, 202, 0};
  const array<unsigned long long, 6> expected_edges = {
      494, 46813, 2531670, 74415, 202, 0};
  if (aggregate_nodes != expected_nodes || aggregate_edges != expected_edges ||
      total_crosschecks != 208)
    throw runtime_error("aggregate certificate mismatch");

  const unsigned long long total_states =
      accumulate(aggregate_nodes.begin(), aggregate_nodes.end(), 0ULL);
  if (total_states != 2653600)
    throw runtime_error("aggregate state total mismatch");

  cout << "aggregate_nodes=" << counts_text(aggregate_nodes)
       << " total_states=" << total_states << "\n";
  cout << "aggregate_edges=" << counts_text(aggregate_edges)
       << " independent_closed_danger_crosschecks=" << total_crosschecks
       << "\n";
  cout << "tournament_guardrail=speed_order_is_transitive_enumeration_only;"
       << " score_histogram=0,1,2,3,4,5; cycles=0; scc=1,1,1,1,1,1;"
       << " Hamiltonian_paths=1\n";
  cout << "faithful_carrier=residual_endpoint_union_x_labelled_remaining_comb_bank\n";
  cout << "PASS: all six product-1/16 primitive-core f=3 rows are loose\n";
  cout << "CONCLUSION: all twenty primitive-core f=3 missing-label rows are closed\n";
  cout << "FRONTIER: primitive-core label rows 923->903\n";
  cout << "SCOPE: 903 primitive-core f<=2 rows and the exceptional odd-row"
       << " primitive mixed-parity branch remain open\n";
  return 0;
}
