// Exact generic terminal recursion for the THM-952 scale-four H6 bank.
//
// Reconstruct the 256 all-order-four common-sheet contexts directly from the
// literal CRT masks.  Their replacement rays are u+52k.  Enumerate the rays
// in numerical order and apply THM-815's exact longest-component cap at every
// prefix.  At depth six, an empty residual is a covering packet and a
// nonempty residual is a certified loose packet.  There is no height cutoff
// and no floating-point arithmetic.
//
// The optional full-safe-tooth and streaming-cap exits are the same
// proof-preserving gates audited for THM-862.  --no-early-cap-gate disables
// them for independent small-shard comparison.
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace std;

struct Rat {
  long long n = 0;
  long long d = 1;

  Rat(long long numerator = 0, long long denominator = 1) {
    if (denominator <= 0)
      throw runtime_error("Rat denominator must be positive");
    long long divisor = gcd(numerator < 0 ? -numerator : numerator,
                            denominator);
    n = numerator / divisor;
    d = denominator / divisor;
  }
  bool operator<(Rat const &other) const {
    return (__int128)n * other.d < (__int128)other.n * d;
  }
  bool operator<=(Rat const &other) const { return !(other < *this); }
  bool operator==(Rat const &other) const {
    return n == other.n && d == other.d;
  }
};

Rat operator-(Rat left, Rat right) {
  long long common = gcd(left.d, right.d);
  long long left_scale = right.d / common;
  long long right_scale = left.d / common;
  __int128 numerator =
      (__int128)left.n * left_scale - (__int128)right.n * right_scale;
  __int128 denominator = (__int128)left.d * left_scale;
  if (numerator < numeric_limits<long long>::min() ||
      numerator > numeric_limits<long long>::max() ||
      denominator > numeric_limits<long long>::max())
    throw runtime_error("Rat subtraction overflow");
  return Rat((long long)numerator, (long long)denominator);
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
  for (int integer = 0; integer < speed; ++integer)
    answer.push_back(
        {Rat(13LL * integer + 1, 13LL * speed),
         Rat(13LL * (integer + 1) - 1, 13LL * speed)});
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

long long floor_mul(Rat value, int speed) {
  if (value.n < 0)
    throw runtime_error("floor_mul expects a nonnegative rational");
  return (long long)((__int128)value.n * speed / value.d);
}

bool contains_full_safe_band(vector<Interval> const &components, int speed) {
  for (Interval const &component : components) {
    long long centre = floor_mul(component.lo, speed);
    long long first = max<long long>(0, centre - 1);
    long long last = min<long long>(speed - 1, centre + 2);
    for (long long integer = first; integer <= last; ++integer) {
      bool left_inside =
          (__int128)component.lo.n * 13 * speed <=
          (__int128)(13 * integer + 1) * component.lo.d;
      bool right_inside =
          (__int128)(13 * integer + 12) * component.hi.d <=
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
  __int128 left = (__int128)22 * remaining * lo.d * hi.d;
  __int128 right = (__int128)candidate * 13 * (13 - 2 * remaining) *
                   length_numerator;
  return left < right;
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

int crt_base(int label, int unit) {
  for (int value = 0; value < 52; ++value)
    if (value % 13 == 4 * label % 13 && value % 4 == unit)
      return value;
  throw runtime_error("scale-four CRT base does not exist");
}

struct Context {
  unsigned root_mask = 0;
  array<int, 6> labels{};
  array<int, 6> units{};
};

bool operator<(Context const &left, Context const &right) {
  return tie(left.labels, left.units) < tie(right.labels, right.units);
}

bool direct_sheet_cover(unsigned root_mask, unsigned unit_mask) {
  for (int owner = 1; owner <= 12; ++owner) {
    if (!(root_mask & (1u << (owner - 1))))
      continue;
    int union_mask = 0;
    int position = 0;
    for (int provider = 1; provider <= 12; ++provider) {
      if (!(root_mask & (1u << (provider - 1))))
        continue;
      int unit = (unit_mask & (1u << position)) ? 3 : 1;
      int base = crt_base(provider, unit);
      int inverse = inverse_mod_13(owner);
      for (int sheet = 0; sheet < 4; ++sheet) {
        int z = centred(base * (inverse + 13 * sheet), 52);
        if (-4 < z && z <= 4)
          union_mask |= 1 << sheet;
      }
      ++position;
    }
    if (union_mask != 15)
      return false;
  }
  return true;
}

vector<Context> reconstruct_contexts(unsigned long long &direct_checks) {
  vector<Context> contexts;
  direct_checks = 0;
  for (unsigned root_mask = 0; root_mask < (1u << 12); ++root_mask) {
    if (popcount(root_mask) != 6)
      continue;
    for (unsigned unit_mask = 0; unit_mask < (1u << 6); ++unit_mask) {
      ++direct_checks;
      if (!direct_sheet_cover(root_mask, unit_mask))
        continue;
      Context context;
      context.root_mask = root_mask;
      int position = 0;
      for (int label = 1; label <= 12; ++label) {
        if (!(root_mask & (1u << (label - 1))))
          continue;
        context.labels[position] = label;
        context.units[position] =
            (unit_mask & (1u << position)) ? 3 : 1;
        ++position;
      }
      if (position != 6)
        throw runtime_error("context label count mismatch");
      contexts.push_back(context);
    }
  }
  sort(contexts.begin(), contexts.end());
  if (adjacent_find(contexts.begin(), contexts.end(),
                    [](Context const &left, Context const &right) {
                      return left.labels == right.labels &&
                             left.units == right.units;
                    }) != contexts.end())
    throw runtime_error("duplicate scale-four context");
  return contexts;
}

int ray_base(Context const &context, int position) {
  return crt_base(context.labels[position], context.units[position]);
}

int next_ray_speed(Context const &context, int position, int last_speed) {
  long long speed = ray_base(context, position);
  if (speed <= last_speed)
    speed += 52LL * ((last_speed - speed) / 52 + 1);
  if (speed > numeric_limits<int>::max())
    throw runtime_error("next ray speed exceeds integer carrier");
  return (int)speed;
}

struct Candidate {
  int speed = 0;
  int position = 0;
};

vector<Candidate> candidates(Context const &context, unsigned used_mask,
                             int last_speed, long long cap) {
  if (cap > numeric_limits<int>::max())
    throw runtime_error("candidate cap exceeds integer carrier");
  vector<Candidate> answer;
  for (int position = 0; position < 6; ++position) {
    if (used_mask & (1u << position))
      continue;
    int speed = next_ray_speed(context, position, last_speed);
    for (; speed <= cap;) {
      answer.push_back({speed, position});
      if (speed > cap - 52)
        break;
      speed += 52;
    }
  }
  sort(answer.begin(), answer.end(), [](Candidate left, Candidate right) {
    return tie(left.speed, left.position) < tie(right.speed, right.position);
  });
  return answer;
}

struct Stats {
  array<unsigned long long, 7> nodes{};
  array<unsigned long long, 7> dead{};
  array<unsigned long long, 7> full_tooth{};
  array<unsigned long long, 7> streaming_cap{};
  unsigned long long candidate_edges = 0;
  unsigned long long covers = 0;
  unsigned long long loose = 0;
  long long maximum_cap = 0;
  int maximum_candidate_speed = 0;
};

void recurse(Context const &context, vector<Interval> const &components,
             int depth, int depth_limit, unsigned used_mask, int last_speed,
             bool early_gates, Stats &stats) {
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
  stats.maximum_cap = max(stats.maximum_cap, cap);
  vector<Candidate> next = candidates(context, used_mask, last_speed, cap);
  if (next.empty())
    ++stats.dead[depth];
  stats.candidate_edges += next.size();

  for (Candidate candidate : next) {
    stats.maximum_candidate_speed =
        max(stats.maximum_candidate_speed, candidate.speed);
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
    recurse(context, child, depth + 1, depth_limit, child_used,
            candidate.speed, early_gates, stats);
  }
}

string array_word(array<unsigned long long, 7> const &values) {
  string answer;
  for (int depth = 0; depth <= 6; ++depth) {
    if (depth)
      answer += ",";
    answer += to_string(values[depth]);
  }
  return answer;
}

string list_word(array<int, 6> const &values) {
  string answer;
  for (int position = 0; position < 6; ++position) {
    if (position)
      answer += ",";
    answer += to_string(values[position]);
  }
  return answer;
}

int main(int argc, char **argv) {
  int context_start = 0;
  int context_limit = numeric_limits<int>::max();
  int depth_limit = 2;
  bool early_gates = true;
  for (int index = 1; index < argc; ++index) {
    string argument = argv[index];
    if (argument == "--context-start" && index + 1 < argc)
      context_start = stoi(argv[++index]);
    else if (argument == "--context-limit" && index + 1 < argc)
      context_limit = stoi(argv[++index]);
    else if (argument == "--depth" && index + 1 < argc)
      depth_limit = stoi(argv[++index]);
    else if (argument == "--no-early-cap-gate")
      early_gates = false;
    else
      throw runtime_error("unknown or incomplete argument: " + argument);
  }

  unsigned long long direct_checks = 0;
  vector<Context> contexts = reconstruct_contexts(direct_checks);
  if (direct_checks != 59'136 || contexts.size() != 256)
    throw runtime_error("scale-four context reconstruction mismatch");
  if (context_start < 0 || context_limit < 0 || depth_limit < 0 ||
      depth_limit > 6 || context_start > (int)contexts.size())
    throw runtime_error("invalid context range");
  int available = (int)contexts.size() - context_start;
  int context_end = context_start + min(context_limit, available);

  cout << "SCALE_FOUR_HAMMING_SIX_GENERIC_RECURSION_SHARD\n";
  cout << "arithmetic=integer+rational floating_point=none height_cutoff=none "
          "depth_limit="
       << depth_limit << " early_gates=" << (early_gates ? 1 : 0) << "\n";
  cout << "context_start=" << context_start << " context_end=" << context_end
       << " context_count=" << context_end - context_start
       << " all_contexts=256 direct_sheet_checks=" << direct_checks << "\n";

  unordered_map<unsigned, vector<Interval>> root_components;
  Stats aggregate;
  auto started = chrono::steady_clock::now();
  for (int context_index = context_start; context_index < context_end;
       ++context_index) {
    Context const &context = contexts[context_index];
    auto root = root_components.find(context.root_mask);
    if (root == root_components.end()) {
      vector<Interval> components = {{Rat(0), Rat(1)}};
      for (int label = 1; label <= 12; ++label)
        if (!(context.root_mask & (1u << (label - 1))))
          components = meet_speed(components, 4 * label);
      if (components.empty())
        throw runtime_error("six-speed retained root is covering");
      root = root_components.emplace(context.root_mask,
                                     std::move(components)).first;
    }

    Stats stats;
    recurse(context, root->second, 0, depth_limit, 0, 0, early_gates, stats);
    for (int depth = 0; depth <= 6; ++depth) {
      aggregate.nodes[depth] += stats.nodes[depth];
      aggregate.dead[depth] += stats.dead[depth];
      aggregate.full_tooth[depth] += stats.full_tooth[depth];
      aggregate.streaming_cap[depth] += stats.streaming_cap[depth];
    }
    aggregate.candidate_edges += stats.candidate_edges;
    aggregate.covers += stats.covers;
    aggregate.loose += stats.loose;
    aggregate.maximum_cap = max(aggregate.maximum_cap, stats.maximum_cap);
    aggregate.maximum_candidate_speed =
        max(aggregate.maximum_candidate_speed, stats.maximum_candidate_speed);

    cout << "GENERIC_ROW|" << context_index << "|" << context.root_mask << "|"
         << list_word(context.labels) << "|" << list_word(context.units) << "|"
         << array_word(stats.nodes) << "|" << array_word(stats.dead) << "|"
         << array_word(stats.full_tooth) << "|"
         << array_word(stats.streaming_cap) << "|" << stats.candidate_edges
         << "|" << stats.covers << "|" << stats.loose << "|"
         << stats.maximum_cap << "|" << stats.maximum_candidate_speed << "\n";

    int completed = context_index - context_start + 1;
    if (completed <= 3 || completed % 16 == 0 ||
        context_index + 1 == context_end) {
      double seconds = chrono::duration<double>(chrono::steady_clock::now() -
                                                started)
                           .count();
      cerr << "contexts=" << completed << "/" << context_end - context_start
           << " index=" << context_index
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
       << "|maximum_cap=" << aggregate.maximum_cap
       << "|maximum_candidate_speed=" << aggregate.maximum_candidate_speed
       << "\nGENERIC_SHARD_DONE\n";
}
