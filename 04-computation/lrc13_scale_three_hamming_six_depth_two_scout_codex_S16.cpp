// Exact scale-three Hamming-six depth-two metric scout.
//
// Reconstruct all 1,504 primitive common-sheet contexts independently from
// the literal CRT sheet predicate, enumerate their step-39 replacement rays
// in numerical order, and apply THM-815's longest-component cap at depths zero
// and one.  Every cap-admissible depth-two child is intersected exactly and
// classified by its exact four-comb continuation cap.  There is no height
// cutoff and no floating point.
//
// Geometry is cached at (missing-label set, x1) and (missing-label set,x1,x2).
// The cache never identifies proof states: every context and every remaining
// labelled ray language is counted separately in the logical census.
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

struct Rat {
  long long n = 0;
  long long d = 1;

  Rat(long long a = 0, long long b = 1) {
    if (b <= 0)
      throw runtime_error("Rat denominator must be positive");
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
  long long g = gcd(a.d, b.d);
  long long as = b.d / g;
  long long bs = a.d / g;
  __int128 nn = (__int128)a.n * as - (__int128)b.n * bs;
  __int128 dd = (__int128)a.d * as;
  if (nn < numeric_limits<long long>::min() ||
      nn > numeric_limits<long long>::max() ||
      dd > numeric_limits<long long>::max())
    throw runtime_error("Rat subtraction overflow");
  return Rat((long long)nn, (long long)dd);
}

ostream &operator<<(ostream &out, Rat q) {
  return out << q.n << "/" << q.d;
}

struct Interval {
  Rat lo;
  Rat hi;
};

unordered_map<int, vector<Interval>> safe_band_cache;

vector<Interval> const &safe_bands(int speed) {
  auto found = safe_band_cache.find(speed);
  if (found != safe_band_cache.end())
    return found->second;
  vector<Interval> answer;
  answer.reserve(speed);
  for (int k = 0; k < speed; ++k)
    answer.push_back({Rat(13LL * k + 1, 13LL * speed),
                      Rat(13LL * (k + 1) - 1, 13LL * speed)});
  return safe_band_cache.emplace(speed, std::move(answer)).first->second;
}

vector<Interval> meet_speed(vector<Interval> const &current, int speed) {
  vector<Interval> const &bands = safe_bands(speed);
  vector<Interval> answer;
  answer.reserve(current.size() + 8);
  size_t i = 0, j = 0;
  while (i < current.size() && j < bands.size()) {
    Rat lo = current[i].lo < bands[j].lo ? bands[j].lo : current[i].lo;
    Rat hi = current[i].hi < bands[j].hi ? current[i].hi : bands[j].hi;
    if (lo < hi)
      answer.push_back({lo, hi});
    if (current[i].hi <= bands[j].hi)
      ++i;
    else
      ++j;
  }
  return answer;
}

long long floor_mul(Rat q, int speed) {
  if (q.n < 0)
    throw runtime_error("floor_mul expects nonnegative rational");
  return (long long)((__int128)q.n * speed / q.d);
}

bool contains_full_safe_band(vector<Interval> const &components, int speed) {
  for (Interval const &component : components) {
    long long centre = floor_mul(component.lo, speed);
    long long first = max<long long>(0, centre - 1);
    long long last = min<long long>(speed - 1, centre + 2);
    for (long long k = first; k <= last; ++k) {
      bool left_inside =
          (__int128)component.lo.n * 13 * speed <=
          (__int128)(13 * k + 1) * component.lo.d;
      bool right_inside =
          (__int128)(13 * k + 12) * component.hi.d <=
          (__int128)component.hi.n * 13 * speed;
      if (left_inside && right_inside)
        return true;
    }
  }
  return false;
}

bool cap_strictly_below(int remaining, Rat lo, Rat hi, int candidate) {
  __int128 length_numerator =
      (__int128)hi.n * lo.d - (__int128)lo.n * hi.d;
  if (length_numerator <= 0)
    throw runtime_error("nonpositive interval in streaming cap gate");
  __int128 lhs = (__int128)22 * remaining * lo.d * hi.d;
  __int128 rhs = (__int128)candidate * 13 * (13 - 2 * remaining) *
                 length_numerator;
  return lhs < rhs;
}

enum class MeetVerdict { materialized, dead_completion, loose_terminal };

MeetVerdict meet_speed_early(vector<Interval> const &current, int speed,
                             int remaining_after, int least_future,
                             vector<Interval> &answer) {
  vector<Interval> const &bands = safe_bands(speed);
  answer.clear();
  answer.reserve(current.size() + 8);
  size_t i = 0, j = 0;
  while (i < current.size() && j < bands.size()) {
    Rat lo = current[i].lo < bands[j].lo ? bands[j].lo : current[i].lo;
    Rat hi = current[i].hi < bands[j].hi ? current[i].hi : bands[j].hi;
    if (lo < hi) {
      if (remaining_after == 0)
        return MeetVerdict::loose_terminal;
      if (cap_strictly_below(remaining_after, lo, hi, least_future))
        return MeetVerdict::dead_completion;
      answer.push_back({lo, hi});
    }
    if (current[i].hi <= bands[j].hi)
      ++i;
    else
      ++j;
  }
  return MeetVerdict::materialized;
}

Interval longest_component(vector<Interval> const &components) {
  if (components.empty())
    throw runtime_error("longest_component on empty residual");
  Interval best = components.front();
  Rat length = best.hi - best.lo;
  for (Interval const &component : components) {
    Rat candidate = component.hi - component.lo;
    if (length < candidate) {
      best = component;
      length = candidate;
    }
  }
  return best;
}

long long discrepancy_cap(int remaining, Rat length) {
  if (remaining < 1 || remaining > 6 || length.n <= 0)
    throw runtime_error("invalid discrepancy cap input");
  __int128 numerator = (__int128)22 * remaining * length.d;
  __int128 denominator =
      (__int128)13 * (13 - 2 * remaining) * length.n;
  __int128 quotient = numerator / denominator;
  if (quotient > numeric_limits<long long>::max())
    throw runtime_error("discrepancy cap overflow");
  return (long long)quotient;
}

int centred(int value, int modulus) {
  int residue = value % modulus;
  if (residue < 0)
    residue += modulus;
  return 2 * residue > modulus ? residue - modulus : residue;
}

int inverse_mod_13(int value) {
  for (int inverse = 1; inverse < 13; ++inverse)
    if (value * inverse % 13 == 1)
      return inverse;
  throw runtime_error("nonunit modulo 13");
}

int crt_base(int label, int order, int unit) {
  for (int value = 0; value < 13 * order; ++value)
    if (value % 13 == order * label % 13 &&
        value % order == unit % order)
      return value;
  throw runtime_error("CRT base does not exist");
}

bool direct_sheet_cover(unsigned order_three_mask, unsigned unit_mask) {
  for (int owner = 1; owner <= 12; ++owner) {
    if (!(order_three_mask & (1u << (owner - 1))))
      continue;
    int sheet_mask = 0;
    int inverse = inverse_mod_13(owner);
    int unit_index = 0;
    for (int provider = 1; provider <= 12; ++provider) {
      if (!(order_three_mask & (1u << (provider - 1))))
        continue;
      int unit = (unit_mask & (1u << unit_index)) ? 2 : 1;
      int base = crt_base(provider, 3, unit);
      for (int sheet = 0; sheet < 3; ++sheet) {
        int z = centred(base * (inverse + 13 * sheet), 39);
        if (-3 < z && z <= 3)
          sheet_mask |= 1 << sheet;
      }
      ++unit_index;
    }
    if (sheet_mask != 7)
      return false;
  }
  return true;
}

int edge_sign(int ratio) {
  if (ratio == 4 || ratio == 5)
    return -1;
  if (ratio == 8 || ratio == 9)
    return 1;
  return 0;
}

bool symbolic_nae_cover(unsigned order_three_mask, unsigned unit_mask) {
  array<int, 13> sign{};
  int unit_index = 0;
  for (int label = 1; label <= 12; ++label)
    if (order_three_mask & (1u << (label - 1))) {
      sign[label] = (unit_mask & (1u << unit_index)) ? -1 : 1;
      ++unit_index;
    }
  for (int owner = 1; owner <= 12; ++owner) {
    if (!(order_three_mask & (1u << (owner - 1))))
      continue;
    bool plus = false, minus = false;
    for (int provider = 1; provider <= 12; ++provider) {
      if (provider == owner ||
          !(order_three_mask & (1u << (provider - 1))))
        continue;
      int ratio = provider * inverse_mod_13(owner) % 13;
      int relation = edge_sign(ratio);
      if (relation * sign[provider] == 1)
        plus = true;
      if (relation * sign[provider] == -1)
        minus = true;
    }
    if (!plus || !minus)
      return false;
  }
  return true;
}

struct Context {
  unsigned root_mask = 0;
  array<int, 6> labels{};
  array<int, 6> orders{};
  array<int, 6> units{};
};

bool operator<(Context const &a, Context const &b) {
  return tie(a.labels, a.orders, a.units) < tie(b.labels, b.orders, b.units);
}

vector<Context> reconstruct_contexts(unsigned long long &direct_checks) {
  vector<pair<unsigned, unsigned>> valid;
  direct_checks = 0;
  for (unsigned c_mask = 1; c_mask < (1u << 12); ++c_mask) {
    int k = popcount(c_mask);
    if (k > 6)
      continue;
    for (unsigned units = 0; units < (1u << k); ++units) {
      bool direct = direct_sheet_cover(c_mask, units);
      bool symbolic = symbolic_nae_cover(c_mask, units);
      if (direct != symbolic)
        throw runtime_error("literal sheet / signed NAE mismatch");
      ++direct_checks;
      if (direct)
        valid.emplace_back(c_mask, units);
    }
  }

  vector<Context> contexts;
  constexpr unsigned all = (1u << 12) - 1;
  for (auto [c_mask, unit_mask] : valid) {
    int k = popcount(c_mask);
    unsigned outside = all ^ c_mask;
    for (unsigned d1_mask = outside;; d1_mask = (d1_mask - 1) & outside) {
      if (popcount(d1_mask) == 6 - k) {
        Context row;
        row.root_mask = c_mask | d1_mask;
        int position = 0, unit_position = 0;
        for (int label = 1; label <= 12; ++label) {
          unsigned bit = 1u << (label - 1);
          if (!(row.root_mask & bit))
            continue;
          row.labels[position] = label;
          if (c_mask & bit) {
            row.orders[position] = 3;
            row.units[position] =
                (unit_mask & (1u << unit_position)) ? 2 : 1;
            ++unit_position;
          } else {
            row.orders[position] = 1;
            row.units[position] = 0;
          }
          ++position;
        }
        if (position != 6 || unit_position != k)
          throw runtime_error("context decoration mismatch");
        contexts.push_back(row);
      }
      if (d1_mask == 0)
        break;
    }
  }
  sort(contexts.begin(), contexts.end());
  if (adjacent_find(contexts.begin(), contexts.end(),
                    [](Context const &a, Context const &b) {
                      return a.labels == b.labels && a.orders == b.orders &&
                             a.units == b.units;
                    }) != contexts.end())
    throw runtime_error("duplicate reconstructed context");
  return contexts;
}

int stratum(Context const &context) {
  int order_one = count(context.orders.begin(), context.orders.end(), 1);
  if (order_one < 0 || order_one > 2)
    throw runtime_error("unexpected scale-three stratum");
  return 2 - order_one; // 0=(2,4), 1=(1,5), 2=(0,6)
}

string type_word(int type) {
  if (type == 0)
    return "2,4";
  if (type == 1)
    return "1,5";
  if (type == 2)
    return "0,6";
  throw runtime_error("bad stratum");
}

int ray_base(Context const &context, int position) {
  if (context.orders[position] == 1)
    return 3 * context.labels[position] + 39;
  return crt_base(context.labels[position], 3, context.units[position]);
}

int next_ray_speed(Context const &context, int position, int last_speed) {
  int speed = ray_base(context, position);
  if (speed <= last_speed)
    speed += 39 * ((last_speed - speed) / 39 + 1);
  return speed;
}

struct Candidate {
  int speed = 0;
  int position = 0;
};

vector<Candidate> candidates(Context const &context, unsigned used_mask,
                             int last_speed, long long cap) {
  vector<Candidate> answer;
  for (int position = 0; position < 6; ++position) {
    if (used_mask & (1u << position))
      continue;
    int speed = next_ray_speed(context, position, last_speed);
    for (; speed <= cap; speed += 39)
      answer.push_back({speed, position});
  }
  sort(answer.begin(), answer.end(), [](Candidate a, Candidate b) {
    return tie(a.speed, a.position) < tie(b.speed, b.position);
  });
  return answer;
}

struct GeometryOneKey {
  unsigned root_mask = 0;
  int x1 = 0;
  bool operator==(GeometryOneKey const &) const = default;
};

struct GeometryOneHash {
  size_t operator()(GeometryOneKey const &key) const {
    return ((size_t)key.root_mask << 32) ^ (unsigned)key.x1;
  }
};

struct GeometryTwoKey {
  unsigned root_mask = 0;
  int x1 = 0;
  int x2 = 0;
  bool operator==(GeometryTwoKey const &) const = default;
};

struct GeometryTwoHash {
  size_t operator()(GeometryTwoKey const &key) const {
    size_t value = key.root_mask;
    value = value * 1000003u + (unsigned)key.x1;
    return value * 1000003u + (unsigned)key.x2;
  }
};

struct GeometryOne {
  vector<Interval> components;
  Rat longest;
  long long cap = 0;
  unsigned long long multiplicity = 0;
};

struct GeometryTwo {
  Interval longest_interval;
  Rat longest;
  long long cap = 0;
};

struct Extremum {
  Rat minimum = Rat(1);
  unsigned long long minimum_count = 0;
  long long maximum_cap = 0;
  bool initialized = false;

  void note(Rat length, long long cap) {
    if (!initialized || length < minimum) {
      minimum = length;
      minimum_count = 1;
      initialized = true;
    } else if (length == minimum) {
      ++minimum_count;
    }
    maximum_cap = max(maximum_cap, cap);
  }
};

struct RowStats {
  unsigned long long depth1 = 0, depth2 = 0;
  unsigned long long first_d1 = 0, first_d3 = 0;
  unsigned long long second_d1 = 0, second_d3 = 0;
  unsigned long long dead0 = 0, dead1 = 0, dead2 = 0;
  unsigned long long stream2 = 0;
  unsigned long long cover1 = 0, cover2 = 0, frontier_live = 0;
  unsigned long long conditioned_flips = 0;
  array<unsigned long long, 11> conditioned_flip_hist{};
  array<Extremum, 3> extrema{};
  int signed_edges = 0, signed_positive = 0, signed_negative = 0;
  int signed_triangles = 0, signed_hamiltonian_paths = 0;
  string signed_scc;
};

vector<int> d3_positions(Context const &context) {
  vector<int> answer;
  for (int i = 0; i < 6; ++i)
    if (context.orders[i] == 3)
      answer.push_back(i);
  return answer;
}

void signed_graph_audit(Context const &context, RowStats &stats) {
  vector<int> vertices = d3_positions(context);
  int k = vertices.size();
  vector<vector<bool>> edge(k, vector<bool>(k));
  for (int i = 0; i < k; ++i)
    for (int j = 0; j < k; ++j) {
      if (i == j)
        continue;
      int ratio = context.labels[vertices[i]] *
                  inverse_mod_13(context.labels[vertices[j]]) % 13;
      int sign = edge_sign(ratio);
      if (sign) {
        edge[i][j] = true;
        ++stats.signed_edges;
        stats.signed_positive += sign == 1;
        stats.signed_negative += sign == -1;
      }
    }
  for (int a = 0; a < k; ++a)
    for (int b = a + 1; b < k; ++b)
      for (int c = b + 1; c < k; ++c)
        stats.signed_triangles +=
            (edge[a][b] && edge[b][c] && edge[c][a]) ||
            (edge[a][c] && edge[c][b] && edge[b][a]);

  vector<vector<bool>> reach = edge;
  for (int i = 0; i < k; ++i)
    reach[i][i] = true;
  for (int mid = 0; mid < k; ++mid)
    for (int i = 0; i < k; ++i)
      for (int j = 0; j < k; ++j)
        reach[i][j] = reach[i][j] || (reach[i][mid] && reach[mid][j]);
  vector<bool> used(k);
  vector<int> sizes;
  for (int i = 0; i < k; ++i) {
    if (used[i])
      continue;
    int size = 0;
    for (int j = 0; j < k; ++j)
      if (!used[j] && reach[i][j] && reach[j][i]) {
        used[j] = true;
        ++size;
      }
    sizes.push_back(size);
  }
  sort(sizes.begin(), sizes.end(), greater<int>());
  ostringstream scc;
  for (size_t i = 0; i < sizes.size(); ++i) {
    if (i)
      scc << ",";
    scc << sizes[i];
  }
  stats.signed_scc = scc.str();

  vector<int> permutation(k);
  iota(permutation.begin(), permutation.end(), 0);
  do {
    bool path = true;
    for (int i = 0; i + 1 < k; ++i)
      path &= edge[permutation[i]][permutation[i + 1]];
    stats.signed_hamiltonian_paths += path;
  } while (next_permutation(permutation.begin(), permutation.end()));
}

int conditioned_flip_count(Context const &context, int chosen_position,
                           int last_speed) {
  int flips = 0;
  for (int a = 0; a < 6; ++a) {
    if (a == chosen_position)
      continue;
    for (int b = a + 1; b < 6; ++b) {
      if (b == chosen_position)
        continue;
      bool raw = ray_base(context, a) < ray_base(context, b);
      bool conditioned = next_ray_speed(context, a, last_speed) <
                         next_ray_speed(context, b, last_speed);
      flips += raw != conditioned;
    }
  }
  return flips;
}

string list_word(array<int, 6> const &values, char separator) {
  ostringstream out;
  for (int i = 0; i < 6; ++i) {
    if (i)
      out << separator;
    out << values[i];
  }
  return out.str();
}

string flip_hist_word(array<unsigned long long, 11> const &histogram) {
  ostringstream out;
  bool first = true;
  for (int i = 0; i <= 10; ++i) {
    if (!histogram[i])
      continue;
    if (!first)
      out << ",";
    first = false;
    out << i << ":" << histogram[i];
  }
  return out.str();
}

string row_line(int index, Context const &context, RowStats const &stats) {
  ostringstream out;
  out << "ROW|" << index << "|" << type_word(stratum(context)) << "|"
      << context.root_mask << "|" << list_word(context.labels, ',') << "|"
      << list_word(context.orders, ',') << "|"
      << list_word(context.units, ',') << "|" << stats.depth1 << "|"
      << stats.depth2 << "|" << stats.first_d1 << "|" << stats.first_d3
      << "|" << stats.second_d1 << "|" << stats.second_d3 << "|"
      << stats.dead0 << "|" << stats.dead1 << "|" << stats.dead2 << "|"
      << stats.stream2 << "|" << stats.cover1 << "|" << stats.cover2 << "|"
      << stats.frontier_live;
  for (int depth = 0; depth <= 2; ++depth)
    out << "|" << stats.extrema[depth].minimum << "|"
        << stats.extrema[depth].minimum_count << "|"
        << stats.extrema[depth].maximum_cap;
  out << "|" << stats.conditioned_flips << "|"
      << flip_hist_word(stats.conditioned_flip_hist) << "|"
      << stats.signed_edges << "|" << stats.signed_positive << "|"
      << stats.signed_negative << "|" << stats.signed_triangles << "|"
      << stats.signed_scc << "|" << stats.signed_hamiltonian_paths;
  return out.str();
}

string lane_language_key(Context const &context, Candidate first) {
  ostringstream out;
  out << context.root_mask << ":" << first.speed << ":";
  for (int position = 0; position < 6; ++position) {
    if (position == first.position)
      continue;
    out << context.labels[position] << "," << context.orders[position] << ","
        << context.units[position] << ","
        << next_ray_speed(context, position, first.speed) << ";";
  }
  return out.str();
}

struct GenericStats {
  array<unsigned long long, 7> nodes{};
  array<unsigned long long, 7> dead{};
  array<unsigned long long, 7> full_tooth{};
  array<unsigned long long, 7> streaming_cap{};
  unsigned long long candidate_edges = 0;
  unsigned long long covers = 0;
  unsigned long long loose = 0;
};

void generic_recurse(Context const &context, vector<Interval> const &components,
                     int depth, int depth_limit, unsigned used_mask,
                     int last_speed, bool early_gates, GenericStats &stats) {
  ++stats.nodes[depth];
  if (components.empty()) {
    if (depth != 6)
      throw runtime_error("nonterminal residual is empty");
    ++stats.covers;
    return;
  }
  if (depth == 6) {
    ++stats.loose;
    return;
  }
  if (depth == depth_limit)
    return;

  Interval longest = longest_component(components);
  long long cap = discrepancy_cap(6 - depth, longest.hi - longest.lo);
  vector<Candidate> next = candidates(context, used_mask, last_speed, cap);
  if (next.empty())
    ++stats.dead[depth];
  stats.candidate_edges += next.size();

  for (Candidate candidate : next) {
    unsigned child_used = used_mask | (1u << candidate.position);
    vector<Interval> child;
    if (early_gates && depth_limit == 6) {
      if (depth >= 2 && contains_full_safe_band(components, candidate.speed)) {
        ++stats.nodes[depth + 1];
        ++stats.full_tooth[depth + 1];
        if (depth == 5)
          ++stats.loose;
        else
          ++stats.dead[depth + 1];
        continue;
      }
      int remaining_after = 5 - depth;
      int least_future = numeric_limits<int>::max();
      if (remaining_after) {
        for (int position = 0; position < 6; ++position)
          if (!(child_used & (1u << position)))
            least_future = min(
                least_future,
                next_ray_speed(context, position, candidate.speed));
        if (least_future == numeric_limits<int>::max())
          throw runtime_error("missing future ray in streaming cap gate");
      }
      MeetVerdict verdict =
          meet_speed_early(components, candidate.speed, remaining_after,
                           least_future, child);
      if (verdict != MeetVerdict::materialized) {
        ++stats.nodes[depth + 1];
        ++stats.streaming_cap[depth + 1];
        if (verdict == MeetVerdict::dead_completion)
          ++stats.dead[depth + 1];
        else
          ++stats.loose;
        continue;
      }
    } else {
      child = meet_speed(components, candidate.speed);
    }
    generic_recurse(context, child, depth + 1, depth_limit, child_used,
                    candidate.speed, early_gates, stats);
  }
}

string array_word(array<unsigned long long, 7> const &values) {
  ostringstream out;
  for (int i = 0; i <= 6; ++i) {
    if (i)
      out << ",";
    out << values[i];
  }
  return out.str();
}

int run_generic(vector<Context> const &contexts, int context_start,
                int context_end, int depth_limit, bool early_gates) {
  unordered_map<unsigned, vector<Interval>> root_components;
  cout << "SCALE_THREE_HAMMING_SIX_GENERIC_RECURSION_SHARD\n";
  cout << "arithmetic=integer+rational floating_point=none height_cutoff=none "
          "depth_limit="
       << depth_limit << " early_gates=" << (early_gates ? 1 : 0) << "\n";
  cout << "context_start=" << context_start << " context_end=" << context_end
       << " context_count=" << context_end - context_start
       << " all_contexts=1504\n";
  GenericStats aggregate;
  auto started = chrono::steady_clock::now();
  for (int index = context_start; index < context_end; ++index) {
    Context const &context = contexts[index];
    auto root = root_components.find(context.root_mask);
    if (root == root_components.end()) {
      vector<Interval> components = {{Rat(0), Rat(1)}};
      for (int label = 1; label <= 12; ++label)
        if (!(context.root_mask & (1u << (label - 1))))
          components = meet_speed(components, 3 * label);
      root = root_components.emplace(context.root_mask,
                                     std::move(components)).first;
    }
    GenericStats stats;
    generic_recurse(context, root->second, 0, depth_limit, 0, 0,
                    early_gates, stats);
    for (int depth = 0; depth <= 6; ++depth) {
      aggregate.nodes[depth] += stats.nodes[depth];
      aggregate.dead[depth] += stats.dead[depth];
      aggregate.full_tooth[depth] += stats.full_tooth[depth];
      aggregate.streaming_cap[depth] += stats.streaming_cap[depth];
    }
    aggregate.candidate_edges += stats.candidate_edges;
    aggregate.covers += stats.covers;
    aggregate.loose += stats.loose;
    cout << "GENERIC_ROW|" << index << "|" << type_word(stratum(context))
         << "|" << array_word(stats.nodes) << "|" << array_word(stats.dead)
         << "|" << array_word(stats.full_tooth) << "|"
         << array_word(stats.streaming_cap) << "|" << stats.candidate_edges
         << "|" << stats.covers << "|" << stats.loose << "\n";
    int completed = index - context_start + 1;
    if (completed <= 3 || completed % 25 == 0 || index + 1 == context_end) {
      double seconds = chrono::duration<double>(chrono::steady_clock::now() -
                                                started)
                           .count();
      cerr << "generic_contexts=" << completed << "/"
           << context_end - context_start << " index=" << index
           << " nodes=" << array_word(aggregate.nodes)
           << " seconds=" << seconds << "\n";
    }
  }
  cout << "GENERIC_SHARD_SUMMARY|nodes=" << array_word(aggregate.nodes)
       << "|dead=" << array_word(aggregate.dead)
       << "|full_tooth=" << array_word(aggregate.full_tooth)
       << "|streaming_cap=" << array_word(aggregate.streaming_cap)
       << "|candidate_edges=" << aggregate.candidate_edges
       << "|covers=" << aggregate.covers << "|loose=" << aggregate.loose
       << "\nGENERIC_SHARD_DONE\n";
  return 0;
}

int main(int argc, char **argv) {
  int context_start = 0;
  int context_limit = numeric_limits<int>::max();
  int depth_limit = 2;
  bool emit_geometry_one = true;
  bool emit_language_lines = true;
  bool early_gates = true;
  for (int i = 1; i < argc; ++i) {
    string argument = argv[i];
    if (argument == "--context-start" && i + 1 < argc)
      context_start = stoi(argv[++i]);
    else if (argument == "--context-limit" && i + 1 < argc)
      context_limit = stoi(argv[++i]);
    else if (argument == "--depth" && i + 1 < argc)
      depth_limit = stoi(argv[++i]);
    else if (argument == "--no-geometry-one-lines")
      emit_geometry_one = false;
    else if (argument == "--no-language-lines")
      emit_language_lines = false;
    else if (argument == "--no-early-cap-gate")
      early_gates = false;
    else
      throw runtime_error("unknown or incomplete argument: " + argument);
  }

  unsigned long long direct_checks = 0;
  vector<Context> contexts = reconstruct_contexts(direct_checks);
  if (direct_checks != 94448 || contexts.size() != 1504)
    throw runtime_error("scale-three context reconstruction mismatch");
  array<int, 3> expected_strata = {336, 672, 496};
  array<int, 3> found_strata{};
  for (Context const &context : contexts)
    ++found_strata[stratum(context)];
  if (found_strata != expected_strata)
    throw runtime_error("scale-three context stratum mismatch");
  if (context_start < 0 || context_limit < 0 || depth_limit < 0 ||
      depth_limit > 6 ||
      context_start > (int)contexts.size())
    throw runtime_error("invalid context range");
  int context_end = min<int>(contexts.size(), context_start + context_limit);

  if (depth_limit != 2)
    return run_generic(contexts, context_start, context_end, depth_limit,
                       early_gates);

  unordered_map<unsigned, vector<Interval>> root_components;
  unordered_map<GeometryOneKey, GeometryOne, GeometryOneHash> geometry_one;
  unordered_map<GeometryTwoKey, GeometryTwo, GeometryTwoHash> geometry_two;
  unordered_set<string> lane_languages;
  vector<string> rows;
  vector<string> dead_depth_two_certificates;
  rows.reserve(context_end - context_start);
  auto started = chrono::steady_clock::now();

  unsigned long long logical_depth1 = 0, logical_depth2 = 0;
  for (int context_index = context_start; context_index < context_end;
       ++context_index) {
    Context const &context = contexts[context_index];
    auto root_found = root_components.find(context.root_mask);
    if (root_found == root_components.end()) {
      vector<Interval> components = {{Rat(0), Rat(1)}};
      for (int label = 1; label <= 12; ++label)
        if (!(context.root_mask & (1u << (label - 1))))
          components = meet_speed(components, 3 * label);
      if (components.empty())
        throw runtime_error("six-speed retained core is covering");
      root_found =
          root_components.emplace(context.root_mask, std::move(components)).first;
    }

    RowStats stats;
    signed_graph_audit(context, stats);
    Rat root_longest = longest_component(root_found->second).hi -
                       longest_component(root_found->second).lo;
    long long root_cap = discrepancy_cap(6, root_longest);
    stats.extrema[0].note(root_longest, root_cap);
    vector<Candidate> first = candidates(context, 0, 0, root_cap);
    if (first.empty())
      ++stats.dead0;

    for (Candidate candidate_one : first) {
      ++stats.depth1;
      ++logical_depth1;
      if (context.orders[candidate_one.position] == 1)
        ++stats.first_d1;
      else
        ++stats.first_d3;
      string language = lane_language_key(context, candidate_one);
      if (!lane_languages.insert(language).second)
        throw runtime_error("distinct logical first-edge languages collided");

      GeometryOneKey key_one{context.root_mask, candidate_one.speed};
      auto one_found = geometry_one.find(key_one);
      if (one_found == geometry_one.end()) {
        GeometryOne geometry;
        geometry.components = meet_speed(root_found->second, candidate_one.speed);
        if (geometry.components.empty())
          throw runtime_error("seven-speed depth-one prefix is covering");
        Interval longest = longest_component(geometry.components);
        geometry.longest = longest.hi - longest.lo;
        geometry.cap = discrepancy_cap(5, geometry.longest);
        one_found = geometry_one.emplace(key_one, std::move(geometry)).first;
      }
      ++one_found->second.multiplicity;
      stats.extrema[1].note(one_found->second.longest, one_found->second.cap);

      int flips = conditioned_flip_count(context, candidate_one.position,
                                         candidate_one.speed);
      stats.conditioned_flips += flips;
      ++stats.conditioned_flip_hist[flips];

      unsigned used_one = 1u << candidate_one.position;
      vector<Candidate> second = candidates(
          context, used_one, candidate_one.speed, one_found->second.cap);
      if (second.empty())
        ++stats.dead1;

      for (Candidate candidate_two : second) {
        ++stats.depth2;
        ++logical_depth2;
        if (context.orders[candidate_two.position] == 1)
          ++stats.second_d1;
        else
          ++stats.second_d3;

        GeometryTwoKey key_two{context.root_mask, candidate_one.speed,
                               candidate_two.speed};
        auto two_found = geometry_two.find(key_two);
        if (two_found == geometry_two.end()) {
          vector<Interval> components =
              meet_speed(one_found->second.components, candidate_two.speed);
          if (components.empty()) {
            ++stats.cover2;
            throw runtime_error("eight-speed depth-two prefix is covering");
          }
          Interval longest = longest_component(components);
          GeometryTwo geometry;
          geometry.longest_interval = longest;
          geometry.longest = longest.hi - longest.lo;
          geometry.cap = discrepancy_cap(4, geometry.longest);
          two_found = geometry_two.emplace(key_two, geometry).first;
        }
        stats.extrema[2].note(two_found->second.longest, two_found->second.cap);

        unsigned used_two = used_one | (1u << candidate_two.position);
        int least_third = numeric_limits<int>::max();
        for (int position = 0; position < 6; ++position)
          if (!(used_two & (1u << position)))
            least_third = min(
                least_third,
                next_ray_speed(context, position, candidate_two.speed));
        if (least_third == numeric_limits<int>::max())
          throw runtime_error("missing third ray at depth-two frontier");
        if (least_third > two_found->second.cap) {
          ++stats.dead2;
          ++stats.stream2;
          ostringstream certificate;
          certificate << "D2DEAD|" << context_index << "|"
                      << type_word(stratum(context)) << "|"
                      << context.root_mask << "|"
                      << list_word(context.labels, ',') << "|"
                      << list_word(context.orders, ',') << "|"
                      << list_word(context.units, ',') << "|"
                      << candidate_one.speed << "|"
                      << context.labels[candidate_one.position] << "|"
                      << candidate_two.speed << "|"
                      << context.labels[candidate_two.position] << "|"
                      << two_found->second.longest_interval.lo << ","
                      << two_found->second.longest_interval.hi << "|"
                      << two_found->second.longest << "|"
                      << two_found->second.cap << "|" << least_third << "|";
          bool first_remaining = true;
          for (int position = 0; position < 6; ++position) {
            if (used_two & (1u << position))
              continue;
            if (!first_remaining)
              certificate << ",";
            first_remaining = false;
            certificate << context.labels[position] << ":"
                        << context.orders[position] << ":"
                        << context.units[position] << ":"
                        << next_ray_speed(context, position,
                                          candidate_two.speed);
          }
          dead_depth_two_certificates.push_back(certificate.str());
        } else {
          ++stats.frontier_live;
        }
      }
    }
    if (stats.depth1 != stats.first_d1 + stats.first_d3 ||
        stats.depth2 != stats.second_d1 + stats.second_d3 ||
        stats.depth2 != stats.dead2 + stats.frontier_live ||
        stats.dead2 != stats.stream2 || stats.cover1 || stats.cover2)
      throw runtime_error("context logical accounting mismatch");
    if (context_index == 1448) {
      const string expected_dead =
          "D2DEAD|1448|1,5|2536|4,6,7,8,9,12|3,3,3,3,3,1|"
          "2,1,1,1,2,0|14|9|38|4|183/494,56/143|115/5434|63|70|"
          "6:3:1:70,7:3:1:73,8:3:1:76,12:1:0:75";
      if (stats.depth1 != 66 || stats.depth2 != 4609 ||
          stats.first_d1 != 10 || stats.first_d3 != 56 ||
          stats.second_d1 != 765 || stats.second_d3 != 3844 ||
          stats.dead2 != 1 || stats.frontier_live != 4608 ||
          dead_depth_two_certificates.empty() ||
          dead_depth_two_certificates.back() != expected_dead)
        throw runtime_error("unique depth-two dead certificate mismatch");
    } else if (stats.dead2 != 0) {
      throw runtime_error("unexpected second depth-two dead context");
    }
    rows.push_back(row_line(context_index, context, stats));

    int completed = context_index - context_start + 1;
    if (completed <= 3 || completed % 25 == 0 || context_index + 1 == context_end) {
      double seconds = chrono::duration<double>(chrono::steady_clock::now() -
                                                started)
                           .count();
      cerr << "contexts=" << completed << "/" << context_end - context_start
           << " index=" << context_index << " depth1=" << logical_depth1
           << " depth2=" << logical_depth2
           << " g1=" << geometry_one.size() << " g2=" << geometry_two.size()
           << " seconds=" << seconds << "\n";
    }
  }

  cout << "SCALE_THREE_HAMMING_SIX_DEPTH_TWO_SHARD\n";
  cout << "arithmetic=integer+rational floating_point=none height_cutoff=none "
          "depth_limit=2\n";
  cout << "context_start=" << context_start
       << " context_end=" << context_end
       << " context_count=" << context_end - context_start
       << " all_contexts=1504 direct_NAE_checks=" << direct_checks << "\n";
  for (string const &row : rows)
    cout << row << "\n";
  for (string const &certificate : dead_depth_two_certificates)
    cout << certificate << "\n";

  if (emit_geometry_one) {
    vector<pair<GeometryOneKey, unsigned long long>> keys;
    keys.reserve(geometry_one.size());
    for (auto const &[key, geometry] : geometry_one)
      keys.push_back({key, geometry.multiplicity});
    sort(keys.begin(), keys.end(), [](auto const &a, auto const &b) {
      return tie(a.first.root_mask, a.first.x1) <
             tie(b.first.root_mask, b.first.x1);
    });
    for (auto const &[key, multiplicity] : keys)
      cout << "G1|" << key.root_mask << "|" << key.x1 << "|"
           << multiplicity << "\n";
  }
  if (emit_language_lines) {
    vector<string> languages(lane_languages.begin(), lane_languages.end());
    sort(languages.begin(), languages.end());
    for (string const &language : languages)
      cout << "L1|" << language << "\n";
  }
  cout << "SHARD_SUMMARY|contexts=" << context_end - context_start
       << "|depth1=" << logical_depth1 << "|depth2=" << logical_depth2
       << "|geometry1=" << geometry_one.size()
       << "|geometry2=" << geometry_two.size()
       << "|language_keys=" << lane_languages.size() << "\n";
  cout << "SHARD_DONE\n";
}
