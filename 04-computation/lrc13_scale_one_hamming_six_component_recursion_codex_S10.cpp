// Exploratory exact scale-one Hamming-six component recursion.
//
// For each six-label deletion core P subset [12], enumerate the six proper
// residue-preserving lifts in increasing numerical order.  At every prefix,
// THM-815's longest-component discrepancy inequality gives a complete bound
// for the next speed.  There is no height cutoff and no floating point.
//
// This file is intentionally exploratory until the full census, independent
// replay, and theorem scope are frozen.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
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
  __int128 nn = (__int128)a.n * b.d - (__int128)b.n * a.d;
  __int128 dd = (__int128)a.d * b.d;
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

long long floor_mul(Rat q, int u) {
  if (q.n < 0)
    throw runtime_error("floor_mul expects nonnegative q");
  return (long long)((__int128)q.n * u / q.d);
}

unordered_map<int, vector<Interval>> safe_band_cache;

vector<Interval> const &safe_bands(int u) {
  auto found = safe_band_cache.find(u);
  if (found != safe_band_cache.end())
    return found->second;
  vector<Interval> bands;
  bands.reserve(u);
  for (int k = 0; k < u; ++k)
    bands.push_back(
        {Rat(13LL * k + 1, 13LL * u),
         Rat(13LL * (k + 1) - 1, 13LL * u)});
  return safe_band_cache.emplace(u, std::move(bands)).first->second;
}

// Intersect two sorted interval unions.  Caching each comb's safe bands is
// substantially faster on the H6 tree, where the same candidate speeds recur
// against many different prefixes.
vector<Interval> meet_speed(vector<Interval> const &current, int u) {
  vector<Interval> const &bands = safe_bands(u);
  vector<Interval> out;
  out.reserve(current.size() + 4);
  size_t i = 0, j = 0;
  while (i < current.size() && j < bands.size()) {
    Rat lo = current[i].lo < bands[j].lo ? bands[j].lo : current[i].lo;
    Rat hi = current[i].hi < bands[j].hi ? current[i].hi : bands[j].hi;
    if (lo < hi)
      out.push_back({lo, hi});
    if (current[i].hi <= bands[j].hi)
      ++i;
    else
      ++j;
  }
  return out;
}

Interval longest_component(vector<Interval> const &components) {
  if (components.empty())
    throw runtime_error("longest_component called on empty residual");
  Interval best = components.front();
  Rat best_length = best.hi - best.lo;
  for (Interval const &component : components) {
    Rat length = component.hi - component.lo;
    if (best_length < length) {
      best = component;
      best_length = length;
    }
  }
  return best;
}

long long discrepancy_cap(int remaining, Rat longest_length) {
  if (remaining < 1 || remaining > 6 || longest_length.n <= 0)
    throw runtime_error("invalid discrepancy-cap inputs");
  __int128 numerator = (__int128)22 * remaining * longest_length.d;
  __int128 denominator =
      (__int128)13 * (13 - 2 * remaining) * longest_length.n;
  return (long long)(numerator / denominator);
}

array<unsigned long long, 7> nodes{};
array<unsigned long long, 7> dead_no_candidate{};
array<long long, 6> max_cap{};
array<Rat, 6> min_longest{};
array<unsigned long long, 6> min_longest_count{};
unsigned long long covering_terminals = 0;
unsigned long long loose_terminals = 0;
unsigned long long root_count = 0;
unsigned long long global_cap_first_prefix_states = 0;
unsigned long long candidate_edges = 0;
array<int, 6> chosen_label{};
array<int, 6> chosen_speed{};
vector<string> covering_rows;
int exploratory_depth_limit = 6;

void note_prefix(int depth, vector<Interval> const &components) {
  Interval longest = longest_component(components);
  Rat length = longest.hi - longest.lo;
  if (length < min_longest[depth]) {
    min_longest[depth] = length;
    min_longest_count[depth] = 1;
  } else if (length == min_longest[depth]) {
    min_longest_count[depth]++;
  }
  int remaining = 6 - depth;
  if (remaining)
    max_cap[depth] = max(max_cap[depth], discrepancy_cap(remaining, length));
}

string row_word(array<int, 6> const &missing) {
  string out = "missing=";
  for (int label : missing)
    out += to_string(label) + ",";
  out += " ordered=";
  for (int i = 0; i < 6; ++i)
    out += to_string(chosen_label[i]) + ":" + to_string(chosen_speed[i]) + ",";
  return out;
}

void recurse(array<int, 6> const &missing,
             vector<Interval> const &components, int depth, int last_speed) {
  nodes[depth]++;

  if (components.empty()) {
    if (depth != 6)
      throw runtime_error("nonterminal residual is empty, contradicting LRC(<=13)");
    covering_terminals++;
    covering_rows.push_back(row_word(missing));
    return;
  }

  if (depth == 6) {
    loose_terminals++;
    return;
  }

  if (depth == exploratory_depth_limit)
    return;

  note_prefix(depth, components);
  int remaining = 6 - depth;
  Interval longest = longest_component(components);
  long long cap = discrepancy_cap(remaining, longest.hi - longest.lo);
  if (cap > 10000000)
    throw runtime_error("exploratory cap exceeded 10^7 guard");

  vector<pair<int, int>> candidates;
  for (int label : missing) {
    bool used = false;
    for (int i = 0; i < depth; ++i)
      used |= chosen_label[i] == label;
    if (used)
      continue;
    long long speed = label + 13;
    if (speed <= last_speed)
      speed += 13 * ((last_speed - speed) / 13 + 1);
    for (; speed <= cap; speed += 13)
      candidates.emplace_back((int)speed, label);
  }
  sort(candidates.begin(), candidates.end());
  if (candidates.empty())
    dead_no_candidate[depth]++;
  candidate_edges += candidates.size();

  for (auto [speed, label] : candidates) {
    chosen_label[depth] = label;
    chosen_speed[depth] = speed;
    vector<Interval> child = meet_speed(components, speed);
    recurse(missing, child, depth + 1, speed);
  }
}

int main(int argc, char **argv) {
  int root_start = 0;
  int root_limit = 924;
  for (int i = 1; i < argc; ++i) {
    string argument = argv[i];
    if (argument == "--depth" && i + 1 < argc)
      exploratory_depth_limit = stoi(argv[++i]);
    else if (argument == "--root-start" && i + 1 < argc)
      root_start = stoi(argv[++i]);
    else if (argument == "--root-limit" && i + 1 < argc)
      root_limit = stoi(argv[++i]);
    else
      throw runtime_error("unknown or incomplete argument: " + argument);
  }
  if (exploratory_depth_limit < 0 || exploratory_depth_limit > 6 ||
      root_start < 0 || root_limit < 0 || root_start + root_limit > 924)
    throw runtime_error("invalid exploratory range");

  for (Rat &value : min_longest)
    value = Rat(1);

  auto started = chrono::steady_clock::now();
  int root_index = 0;
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
              array<int, 6> missing = {a, b, c, d, e, f};
              bool removed[13] = {};
              for (int label : missing)
                removed[label] = true;
              // THM-815 records the coarser census obtained by applying the
              // worst-core cap 468 uniformly.  The recursion below uses each
              // core's own longest component, so its nodes[1] count is a
              // strictly sharper quantity.
              for (int label : missing)
                for (int speed = label + 13; speed <= 468; speed += 13)
                  global_cap_first_prefix_states++;
              vector<Interval> components = {{Rat(0), Rat(1)}};
              for (int label = 1; label <= 12; ++label)
                if (!removed[label])
                  components = meet_speed(components, label);
              if (components.empty())
                throw runtime_error("six-speed deletion core has empty residual");
              root_count++;
              recurse(missing, components, 0, 0);
              if (root_count <= 10 || root_count % 50 == 0) {
                double seconds = chrono::duration<double>(
                                     chrono::steady_clock::now() - started)
                                     .count();
                cerr << "roots=" << root_count << "/" << root_limit
                     << " index=" << this_root << " depth=";
                for (auto value : nodes)
                  cerr << value << ",";
                cerr << " covers=" << covering_terminals
                     << " seconds=" << seconds << "\n";
              }
            }

  double seconds =
      chrono::duration<double>(chrono::steady_clock::now() - started).count();
  if (root_start == 0 && root_limit == 924 && exploratory_depth_limit >= 1 &&
      (root_count != 924 || global_cap_first_prefix_states != 194040 ||
       nodes[1] > global_cap_first_prefix_states))
    throw runtime_error("THM-815 root/first-prefix census mismatch");

  cout << "EXPLORATORY_SCALE_ONE_HAMMING_SIX_COMPONENT_RECURSION\n";
  cout << "arithmetic=integer+rational floating_point=none height_cutoff=none\n";
  cout << "root_start=" << root_start << " root_limit=" << root_limit
       << " depth_limit=" << exploratory_depth_limit << "\n";
  cout << "roots=" << root_count << " candidate_edges=" << candidate_edges
       << " runtime_seconds=" << seconds << "\n";
  cout << "global_cap_first_prefix_states="
       << global_cap_first_prefix_states
       << " local_cap_first_prefix_states=" << nodes[1] << "\n";
  cout << "nodes_depth0..6=";
  for (auto value : nodes)
    cout << value << ",";
  cout << "\n";
  cout << "dead_no_candidate_depth0..6=";
  for (auto value : dead_no_candidate)
    cout << value << ",";
  cout << "\n";
  for (int depth = 0; depth < 6; ++depth)
    cout << "depth=" << depth << " min_longest=" << min_longest[depth]
         << " min_count=" << min_longest_count[depth]
         << " max_next_cap=" << max_cap[depth] << "\n";
  cout << "covering_terminals=" << covering_terminals
       << " loose_terminals=" << loose_terminals << "\n";
  for (string const &row : covering_rows)
    cout << "COVER " << row << "\n";
  cout << "EXPLORATORY_DONE\n";
}
