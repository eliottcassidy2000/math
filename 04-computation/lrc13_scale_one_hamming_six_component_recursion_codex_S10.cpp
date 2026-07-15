// Exact scale-one Hamming-six component recursion.
//
// For each six-label deletion core P subset [12], enumerate the six proper
// residue-preserving lifts in increasing numerical order.  At every prefix,
// THM-815's longest-component discrepancy inequality gives a complete bound
// for the next speed.  There is no height cutoff and no floating point.
//
// The full 924-root census is frozen below.  A companion implementation
// independently reconstructs every expanded prefix from its closed danger
// union and certifies the shortcut witnesses.
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
  // Form the difference over the lcm denominator.  Cancelling the common
  // denominator factor before multiplying is important on the deep H6 tree:
  // the exact endpoints are small, but their unreduced denominator product
  // need not be.
  long long g = gcd(a.d, b.d);
  long long a_scale = b.d / g;
  long long b_scale = a.d / g;
  __int128 nn = (__int128)a.n * a_scale - (__int128)b.n * b_scale;
  __int128 dd = (__int128)a.d * a_scale;
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

enum class MeetVerdict { materialized, dead_completion, loose_terminal };

bool contains_full_safe_band(vector<Interval> const &components, int u) {
  for (Interval const &component : components) {
    long long centre = floor_mul(component.lo, u);
    long long first = max<long long>(0, centre - 1);
    long long last = min<long long>(u - 1, centre + 2);
    for (long long k = first; k <= last; ++k) {
      bool left_inside =
          (__int128)component.lo.n * 13 * u <=
          (__int128)(13 * k + 1) * component.lo.d;
      bool right_inside =
          (__int128)(13 * k + 12) * component.hi.d <=
          (__int128)component.hi.n * 13 * u;
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
    throw runtime_error("nonpositive interval in early cap gate");
  __int128 lhs = (__int128)22 * remaining * lo.d * hi.d;
  __int128 rhs = (__int128)candidate * 13 * (13 - 2 * remaining) *
                 length_numerator;
  return lhs < rhs;
}

// Fuse intersection with the child's next-cap test.  If one emerging child
// component already makes THM-815's cap smaller than the least legal future
// lift, the actual longest child component can only make that cap smaller.
// At a terminal, the first emerging component directly certifies looseness.
MeetVerdict meet_speed_early(vector<Interval> const &current, int u,
                             int remaining_after, int least_future,
                             vector<Interval> &out) {
  vector<Interval> const &bands = safe_bands(u);
  out.clear();
  out.reserve(current.size() + 4);
  size_t i = 0, j = 0;
  while (i < current.size() && j < bands.size()) {
    Rat lo = current[i].lo < bands[j].lo ? bands[j].lo : current[i].lo;
    Rat hi = current[i].hi < bands[j].hi ? current[i].hi : bands[j].hi;
    if (lo < hi) {
      if (remaining_after == 0)
        return MeetVerdict::loose_terminal;
      if (cap_strictly_below(remaining_after, lo, hi, least_future))
        return MeetVerdict::dead_completion;
      out.push_back({lo, hi});
    }
    if (current[i].hi <= bands[j].hi)
      ++i;
    else
      ++j;
  }
  return MeetVerdict::materialized;
}

array<unsigned long long, 7> nodes{};
array<unsigned long long, 7> dead_no_candidate{};
array<unsigned long long, 7> full_safe_band_prunes{};
array<unsigned long long, 7> early_cap_prunes{};
array<unsigned long long, 6> materialized_prefixes{};
array<unsigned long long, 16> tournament_flip_histogram{};
array<long long, 6> max_cap{};
array<Rat, 6> min_longest{};
array<unsigned long long, 6> min_longest_count{};
unsigned long long covering_terminals = 0;
unsigned long long loose_terminals = 0;
unsigned long long root_count = 0;
unsigned long long global_cap_first_prefix_states = 0;
unsigned long long candidate_edges = 0;
unsigned long long tournament_total_flips = 0;
unsigned long long tournament_conditioned_ties = 0;
array<int, 6> chosen_label{};
array<int, 6> chosen_speed{};
vector<string> covering_rows;
int exploratory_depth_limit = 6;
bool use_early_cap_gate = true;

// Root-level Tournament Analysis.  The raw gauge orders the six missing
// labels by their least proper lifts.  The conditioned gauge orders them by
// the exact five-comb cap left after adjoining that least lift to the whole
// deletion core; ties use increasing label.  Each gauge is therefore
// transitive by construction.  Edge flips measure how much the literal core
// changes this pairwise planning telemetry, not the cover predicate.
void audit_root_tournament(array<int, 6> const &missing,
                           vector<Interval> const &components) {
  array<long long, 6> conditioned_cap{};
  for (int i = 0; i < 6; ++i) {
    vector<Interval> child = meet_speed(components, missing[i] + 13);
    if (child.empty())
      throw runtime_error("root tournament child unexpectedly covers");
    Interval longest = longest_component(child);
    conditioned_cap[i] =
        discrepancy_cap(5, longest.hi - longest.lo);
  }
  int flips = 0;
  for (int i = 0; i < 6; ++i)
    for (int j = i + 1; j < 6; ++j) {
      // The raw winner is i because missing[] is strictly increasing.
      int conditioned_winner = i;
      if (conditioned_cap[j] < conditioned_cap[i])
        conditioned_winner = j;
      else if (conditioned_cap[i] == conditioned_cap[j])
        tournament_conditioned_ties++;
      flips += conditioned_winner != i;
    }
  tournament_flip_histogram[flips]++;
  tournament_total_flips += flips;
}

void note_prefix(int depth, vector<Interval> const &components) {
  materialized_prefixes[depth]++;
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
    vector<Interval> child;
    if (use_early_cap_gate && exploratory_depth_limit == 6) {
      if (depth >= 2 && contains_full_safe_band(components, speed)) {
        nodes[depth + 1]++;
        full_safe_band_prunes[depth + 1]++;
        if (depth == 5)
          loose_terminals++;
        else
          dead_no_candidate[depth + 1]++;
        continue;
      }
      int remaining_after = 5 - depth;
      int least_future = numeric_limits<int>::max();
      if (remaining_after) {
        for (int future_label : missing) {
          bool used = future_label == label;
          for (int i = 0; i < depth; ++i)
            used |= chosen_label[i] == future_label;
          if (used)
            continue;
          int future_speed = future_label + 13;
          if (future_speed <= speed)
            future_speed += 13 * ((speed - future_speed) / 13 + 1);
          least_future = min(least_future, future_speed);
        }
        if (least_future == numeric_limits<int>::max())
          throw runtime_error("missing least future lift");
      }
      MeetVerdict verdict = meet_speed_early(
          components, speed, remaining_after, least_future, child);
      if (verdict != MeetVerdict::materialized) {
        nodes[depth + 1]++;
        if (verdict == MeetVerdict::dead_completion) {
          dead_no_candidate[depth + 1]++;
          early_cap_prunes[depth + 1]++;
        } else {
          loose_terminals++;
          early_cap_prunes[depth + 1]++;
        }
        continue;
      }
    } else {
      child = meet_speed(components, speed);
    }
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
    else if (argument == "--no-early-cap-gate")
      use_early_cap_gate = false;
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
              audit_root_tournament(missing, components);
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

  if (root_start == 0 && root_limit == 924 && exploratory_depth_limit == 6 &&
      use_early_cap_gate) {
    const array<unsigned long long, 7> expected_nodes = {
        924, 83881, 8906315, 559202706, 12671505, 53812, 21};
    const array<unsigned long long, 7> expected_dead = {
        0, 0, 0, 555565824, 12638291, 53792, 0};
    const array<unsigned long long, 7> expected_full = {
        0, 0, 0, 495797163, 0, 0, 0};
    const array<unsigned long long, 7> expected_early = {
        0, 0, 0, 59768661, 12638291, 53792, 20};
    const array<unsigned long long, 6> expected_materialized = {
        924, 83881, 8906315, 3636882, 33214, 20};
    const array<Rat, 6> expected_min_longest = {
        Rat(31, 1430), Rat(11, 6071), Rat(11, 20215),
        Rat(2, 5941), Rat(2251, 2834949), Rat(1, 325)};
    const array<long long, 6> expected_max_cap = {
        468, 1556, 2488, 2154, 473, 50};
    const array<unsigned long long, 6> expected_min_count = {
        1, 1, 1, 1, 1, 1};
    const array<unsigned long long, 16> expected_flip_histogram = {
        47, 61, 93, 115, 146, 108, 110, 78,
        49, 46, 33, 18, 16, 4, 0, 0};
    const vector<string> expected_covers = {
        "missing=1,3,5,7,9,11, ordered="
        "1:14,3:16,5:18,7:20,9:22,11:24,"};
    if (nodes != expected_nodes || dead_no_candidate != expected_dead ||
        full_safe_band_prunes != expected_full ||
        early_cap_prunes != expected_early ||
        materialized_prefixes != expected_materialized ||
        min_longest != expected_min_longest ||
        min_longest_count != expected_min_count || max_cap != expected_max_cap ||
        candidate_edges != 580918240ULL || covering_terminals != 1 ||
        loose_terminals != 20 || covering_rows != expected_covers ||
        tournament_flip_histogram != expected_flip_histogram ||
        tournament_conditioned_ties != 406 || tournament_total_flips != 4500)
      throw runtime_error("frozen full H6 census mismatch");
  }

  cout << "SCALE_ONE_HAMMING_SIX_COMPONENT_RECURSION\n";
  cout << "arithmetic=integer+rational floating_point=none height_cutoff=none\n";
  cout << "root_start=" << root_start << " root_limit=" << root_limit
       << " depth_limit=" << exploratory_depth_limit << "\n";
  cout << "roots=" << root_count << " candidate_edges=" << candidate_edges
       << "\n";
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
  cout << "full_safe_band_prunes_depth0..6=";
  for (auto value : full_safe_band_prunes)
    cout << value << ",";
  cout << "\n";
  cout << "early_cap_prunes_depth0..6=";
  for (auto value : early_cap_prunes)
    cout << value << ",";
  cout << "\n";
  for (int depth = 0; depth < 6; ++depth)
    cout << "materialized_depth=" << depth
         << " prefixes=" << materialized_prefixes[depth]
         << " min_longest=" << min_longest[depth]
         << " min_count=" << min_longest_count[depth]
         << " max_next_cap=" << max_cap[depth] << "\n";
  cout << "covering_terminals=" << covering_terminals
       << " loose_terminals=" << loose_terminals << "\n";
  for (string const &row : covering_rows)
    cout << "COVER " << row << "\n";
  cout << "TOURNAMENT_ANALYSIS\n";
  cout << "vertices=missing_labels raw=least_proper_lift_order "
          "conditioned=five_comb_cap_after_least_lift\n";
  cout << "tie_hamiltonian_order=increasing_label raw_and_conditioned_"
          "fingerprints=transitive_scores_0_to_5_triangles_0_"
          "singleton_SCCs_6_Hamiltonian_paths_1\n";
  cout << "conditioned_ties=" << tournament_conditioned_ties
       << " total_edge_flips=" << tournament_total_flips << "\n";
  cout << "edge_flip_histogram=";
  for (int flips = 0; flips <= 15; ++flips)
    if (tournament_flip_histogram[flips])
      cout << flips << ":" << tournament_flip_histogram[flips] << ",";
  cout << "\n";
  cout << "tournament_verdict=pairwise_planning_telemetry_only;_cover_"
          "predicate_requires_literal_components_and_remaining_combs\n";
  cout << "CERTIFIED_DONE\n";
  cerr << "runtime_seconds=" << seconds << "\n";
}
