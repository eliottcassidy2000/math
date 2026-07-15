// Geometry-batched depth-three audit for the c=3 Hamming-six branch.
//
// This deliberately imports the independently checked depth-two engine so the
// CRT context reconstruction and exact interval arithmetic have one source of
// truth.  The cache is an execution cache only: every labelled lane is still
// enumerated and retains its context, used-label mask, least future ray, and
// shortcut verdict.  Only literal metric geometry is shared.
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main lrc13_depth_two_imported_main
#include "lrc13_scale_three_hamming_six_depth_two_scout_codex_S16.cpp"
#undef main
#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#include <ctime>
#include <optional>

enum class ThirdVerdict : uint64_t {
  materialized = 1,
  full_tooth = 2,
  streaming_dead = 3,
  covering = 4,
};

struct LaneCount {
  unsigned long long d1 = 0, d2 = 0, d3 = 0;
  unsigned long long tooth3 = 0, stream3 = 0, materialized3 = 0, cover3 = 0;
  unsigned long long flips2 = 0;
  array<unsigned long long, 7> flip_hist2{};
  bool operator==(LaneCount const &) const = default;
};

struct Audit {
  vector<LaneCount> rows;
  unsigned long long roots = 0, meet1 = 0, meet2 = 0, meet3 = 0;
  unsigned long long unique1 = 0, unique2 = 0, unique3 = 0;
  uint64_t fingerprint_xor = 0, fingerprint_sum = 0;
  double seconds = 0, cpu_seconds = 0;
};

bool literal_crosscheck = false;

void mix64(uint64_t &h, uint64_t value) {
  h ^= value;
  h *= 1099511628211ULL;
}

void mix_rat(uint64_t &h, Rat q) {
  mix64(h, static_cast<uint64_t>(q.n));
  mix64(h, static_cast<uint64_t>(q.d));
}

void mix_components(uint64_t &h, vector<Interval> const &components) {
  mix64(h, components.size());
  for (Interval const &component : components) {
    mix_rat(h, component.lo);
    mix_rat(h, component.hi);
  }
}

bool same_components(vector<Interval> const &a, vector<Interval> const &b) {
  if (a.size() != b.size())
    return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (!(a[i].lo == b[i].lo) || !(a[i].hi == b[i].hi))
      return false;
  return true;
}

// Intersect directly with the speed-s comb without materializing and retaining
// all s safe bands.  A component meets only the bands around floor(s*t), so
// this is both exact and memory-bounded.
vector<Interval> metric_meet(vector<Interval> const &current, int speed) {
  vector<Interval> answer;
  answer.reserve(current.size() + 8);
  for (Interval const &component : current) {
    long long first = max<long long>(0, floor_mul(component.lo, speed) - 1);
    for (long long k = first; k < speed; ++k) {
      Interval band{Rat(13 * k + 1, 13LL * speed),
                    Rat(13 * (k + 1) - 1, 13LL * speed)};
      if (!(band.lo < component.hi))
        break;
      if (band.hi <= component.lo)
        continue;
      Rat lo = component.lo < band.lo ? band.lo : component.lo;
      Rat hi = component.hi < band.hi ? component.hi : band.hi;
      if (lo < hi)
        answer.push_back({lo, hi});
    }
  }
  if (literal_crosscheck) {
    vector<Interval> literal = meet_speed(current, speed);
    if (!same_components(answer, literal))
      throw runtime_error("streaming comb / literal safe-band mismatch");
  }
  return answer;
}

int remaining_ray_flips(Context const &context, unsigned used_mask,
                        int last_speed) {
  vector<int> remaining;
  for (int p = 0; p < 6; ++p)
    if (!(used_mask & (1u << p)))
      remaining.push_back(p);
  int flips = 0;
  for (size_t i = 0; i < remaining.size(); ++i)
    for (size_t j = i + 1; j < remaining.size(); ++j) {
      int a = remaining[i], b = remaining[j];
      bool raw = make_pair(ray_base(context, a), context.labels[a]) <
                 make_pair(ray_base(context, b), context.labels[b]);
      bool conditioned =
          make_pair(next_ray_speed(context, a, last_speed), context.labels[a]) <
          make_pair(next_ray_speed(context, b, last_speed), context.labels[b]);
      flips += raw != conditioned;
    }
  return flips;
}

int least_future(Context const &context, unsigned used_mask, int last_speed) {
  int answer = numeric_limits<int>::max();
  for (int p = 0; p < 6; ++p)
    if (!(used_mask & (1u << p)))
      answer = min(answer, next_ray_speed(context, p, last_speed));
  return answer;
}

bool streaming_dead(vector<Interval> const &child, int least) {
  if (least == numeric_limits<int>::max())
    throw runtime_error("depth-three lane lost its three future rays");
  for (Interval const &component : child)
    if (cap_strictly_below(3, component.lo, component.hi, least))
      return true;
  return false;
}

vector<Interval> root_geometry(Context const &context) {
  vector<Interval> components = {{Rat(0), Rat(1)}};
  for (int label = 1; label <= 12; ++label)
    if (!(context.root_mask & (1u << (label - 1))))
      components = metric_meet(components, 3 * label);
  if (components.empty())
    throw runtime_error("covering retained core");
  return components;
}

struct OneValue {
  vector<Interval> components;
  long long cap = 0;
};

struct TwoValue {
  struct Lane {
    int context_index = 0;
    Candidate one, two;
    unsigned used2 = 0;
  };
  vector<Lane> lanes;
};

void note_lane(Audit &audit, int local_index, int context_index,
               Candidate one, Candidate two, Candidate three, unsigned used_two,
               vector<Interval> const &child, ThirdVerdict verdict,
               int least) {
  LaneCount &row = audit.rows[local_index];
  ++row.d3;
  if (verdict == ThirdVerdict::full_tooth)
    ++row.tooth3;
  else if (verdict == ThirdVerdict::streaming_dead)
    ++row.stream3;
  else if (verdict == ThirdVerdict::covering)
    ++row.cover3;
  else
    ++row.materialized3;

  unsigned used_three = used_two | (1u << three.position);
  uint64_t lane_fingerprint = 1469598103934665603ULL;
  mix64(lane_fingerprint, static_cast<uint64_t>(context_index));
  mix64(lane_fingerprint, used_three);
  mix64(lane_fingerprint, static_cast<uint64_t>(one.speed));
  mix64(lane_fingerprint, static_cast<uint64_t>(two.speed));
  mix64(lane_fingerprint, static_cast<uint64_t>(three.speed));
  mix64(lane_fingerprint, static_cast<uint64_t>(verdict));
  mix64(lane_fingerprint, static_cast<uint64_t>(least));
  mix_components(lane_fingerprint, child);
  audit.fingerprint_xor ^= lane_fingerprint;
  audit.fingerprint_sum += lane_fingerprint;
}

Audit run_raw(vector<Context> const &contexts, int start, int end,
              bool gate_audit) {
  Audit audit;
  audit.rows.resize(end - start);
  unordered_map<unsigned, vector<Interval>> roots;
  auto begun = chrono::steady_clock::now();
  clock_t cpu_begun = clock();
  for (int index = start; index < end; ++index) {
    Context const &context = contexts[index];
    LaneCount &row = audit.rows[index - start];
    auto root = roots.find(context.root_mask);
    if (root == roots.end()) {
      root = roots.emplace(context.root_mask, root_geometry(context)).first;
      ++audit.roots;
    }
    Interval root_longest = longest_component(root->second);
    vector<Candidate> first = candidates(
        context, 0, 0, discrepancy_cap(6, root_longest.hi - root_longest.lo));
    for (Candidate one : first) {
      ++row.d1;
      ++audit.meet1;
      vector<Interval> child1 = metric_meet(root->second, one.speed);
      Interval longest1 = longest_component(child1);
      vector<Candidate> second = candidates(
          context, 1u << one.position, one.speed,
          discrepancy_cap(5, longest1.hi - longest1.lo));
      for (Candidate two : second) {
        ++row.d2;
        ++audit.meet2;
        vector<Interval> child2 = metric_meet(child1, two.speed);
        Interval longest2 = longest_component(child2);
        unsigned used2 = (1u << one.position) | (1u << two.position);
        int flips = remaining_ray_flips(context, used2, two.speed);
        row.flips2 += flips;
        ++row.flip_hist2[flips];
        vector<Candidate> third = candidates(
            context, used2, two.speed,
            discrepancy_cap(4, longest2.hi - longest2.lo));
        for (Candidate three : third) {
          int least = least_future(context, used2 | (1u << three.position),
                                   three.speed);
          vector<Interval> child3;
          ThirdVerdict verdict;
          if (gate_audit && contains_full_safe_band(child2, three.speed)) {
            verdict = ThirdVerdict::full_tooth;
          } else {
            ++audit.meet3;
            child3 = metric_meet(child2, three.speed);
            if (child3.empty())
              verdict = ThirdVerdict::covering;
            else if (gate_audit && streaming_dead(child3, least))
              verdict = ThirdVerdict::streaming_dead;
            else
              verdict = ThirdVerdict::materialized;
          }
          note_lane(audit, index - start, index, one, two, three,
                    used2, child3, verdict, least);
        }
      }
    }
  }
  audit.seconds =
      chrono::duration<double>(chrono::steady_clock::now() - begun).count();
  audit.cpu_seconds = double(clock() - cpu_begun) / CLOCKS_PER_SEC;
  return audit;
}

Audit run_batched(vector<Context> const &contexts, int start, int end,
                  bool gate_audit) {
  Audit audit;
  audit.rows.resize(end - start);
  auto begun = chrono::steady_clock::now();
  clock_t cpu_begun = clock();
  // Contexts with one root mask are contiguous.  No metric key can cross root
  // masks, so flushing after each root block bounds memory without losing a
  // single possible reuse.
  for (int block_start = start; block_start < end;) {
    int block_end = block_start + 1;
    while (block_end < end &&
           contexts[block_end].root_mask == contexts[block_start].root_mask)
      ++block_end;
    ++audit.roots;
    vector<Interval> root = root_geometry(contexts[block_start]);
    Interval root_longest = longest_component(root);
    long long root_cap = discrepancy_cap(6, root_longest.hi - root_longest.lo);
    unordered_map<GeometryOneKey, OneValue, GeometryOneHash> one_cache;
    unordered_map<GeometryTwoKey, TwoValue, GeometryTwoHash> two_cache;

    // Phase one keeps each proof lane but attaches it to its metric parent.
    for (int index = block_start; index < block_end; ++index) {
      Context const &context = contexts[index];
      LaneCount &row = audit.rows[index - start];
      vector<Candidate> first = candidates(context, 0, 0, root_cap);
      for (Candidate one : first) {
        ++row.d1;
        GeometryOneKey key1{context.root_mask, one.speed};
        auto found1 = one_cache.find(key1);
        if (found1 == one_cache.end()) {
          OneValue value;
          ++audit.meet1;
          value.components = metric_meet(root, one.speed);
          Interval longest = longest_component(value.components);
          value.cap = discrepancy_cap(5, longest.hi - longest.lo);
          found1 = one_cache.emplace(key1, std::move(value)).first;
        }
        vector<Candidate> second = candidates(
            context, 1u << one.position, one.speed, found1->second.cap);
        for (Candidate two : second) {
          ++row.d2;
          GeometryTwoKey key2{context.root_mask, one.speed, two.speed};
          auto found2 = two_cache.find(key2);
          if (found2 == two_cache.end())
            found2 = two_cache.emplace(key2, TwoValue{}).first;
          unsigned used2 = (1u << one.position) | (1u << two.position);
          int flips = remaining_ray_flips(context, used2, two.speed);
          row.flips2 += flips;
          ++row.flip_hist2[flips];
          found2->second.lanes.push_back({index, one, two, used2});
        }
      }
    }
    audit.unique1 += one_cache.size();
    audit.unique2 += two_cache.size();

    // Phase two transposes lanes by proposed x3 inside each metric parent.
    // The child geometry is computed once, consumed by all corresponding
    // labelled lanes, and immediately discarded.
    struct ThirdLane {
      TwoValue::Lane lane;
      Candidate three;
      int least = 0;
    };
    for (auto &[key2, value] : two_cache) {
      auto found1 = one_cache.find({key2.root_mask, key2.x1});
      if (found1 == one_cache.end())
        throw runtime_error("depth-two batch lost its depth-one parent");
      ++audit.meet2;
      vector<Interval> parent =
          metric_meet(found1->second.components, key2.x2);
      Interval longest = longest_component(parent);
      long long parent_cap = discrepancy_cap(4, longest.hi - longest.lo);
      map<int, vector<ThirdLane>> by_speed;
      for (TwoValue::Lane const &lane : value.lanes) {
        Context const &context = contexts[lane.context_index];
        vector<Candidate> third =
            candidates(context, lane.used2, lane.two.speed, parent_cap);
        for (Candidate candidate : third) {
          int least = least_future(
              context, lane.used2 | (1u << candidate.position), candidate.speed);
          by_speed[candidate.speed].push_back({lane, candidate, least});
        }
      }
      for (auto &[speed, lanes] : by_speed) {
        ++audit.unique3;
        bool tooth = gate_audit && contains_full_safe_band(parent, speed);
        vector<Interval> child;
        if (!tooth) {
          ++audit.meet3;
          child = metric_meet(parent, speed);
        }
        for (ThirdLane const &third_lane : lanes) {
          ThirdVerdict verdict;
          if (tooth)
            verdict = ThirdVerdict::full_tooth;
          else if (child.empty())
            verdict = ThirdVerdict::covering;
          else if (gate_audit && streaming_dead(child, third_lane.least))
            verdict = ThirdVerdict::streaming_dead;
          else
            verdict = ThirdVerdict::materialized;
          TwoValue::Lane const &lane = third_lane.lane;
          note_lane(audit, lane.context_index - start, lane.context_index,
                    lane.one, lane.two, third_lane.three, lane.used2, child,
                    verdict, third_lane.least);
        }
      }
    }
    block_start = block_end;
  }
  audit.seconds =
      chrono::duration<double>(chrono::steady_clock::now() - begun).count();
  audit.cpu_seconds = double(clock() - cpu_begun) / CLOCKS_PER_SEC;
  return audit;
}

string flip_hist(vector<LaneCount> const &rows) {
  array<unsigned long long, 7> total{};
  for (LaneCount const &row : rows)
    for (int i = 0; i <= 6; ++i)
      total[i] += row.flip_hist2[i];
  ostringstream out;
  for (int i = 0; i <= 6; ++i) {
    if (i)
      out << ",";
    out << i << ":" << total[i];
  }
  return out.str();
}

void print_audit(string const &name, Audit const &audit) {
  LaneCount total;
  for (LaneCount const &row : audit.rows) {
    total.d1 += row.d1;
    total.d2 += row.d2;
    total.d3 += row.d3;
    total.tooth3 += row.tooth3;
    total.stream3 += row.stream3;
    total.materialized3 += row.materialized3;
    total.cover3 += row.cover3;
    total.flips2 += row.flips2;
  }
  cout << name << "|d1=" << total.d1 << "|d2=" << total.d2
       << "|d3=" << total.d3 << "|tooth3=" << total.tooth3
       << "|stream3=" << total.stream3
       << "|materialized3=" << total.materialized3
       << "|cover3=" << total.cover3 << "|roots=" << audit.roots
       << "|meet1=" << audit.meet1 << "|meet2=" << audit.meet2
       << "|meet3=" << audit.meet3 << "|unique1=" << audit.unique1
       << "|unique2=" << audit.unique2 << "|unique3=" << audit.unique3
       << "|seconds=" << audit.seconds << "|cpu_seconds=" << audit.cpu_seconds
       << "|fingerprint_xor=" << hex
       << audit.fingerprint_xor << "|fingerprint_sum=" << audit.fingerprint_sum
       << dec << "|ray_flips2=" << total.flips2
       << "|ray_flip_hist2=" << flip_hist(audit.rows) << "\n";
}

int main(int argc, char **argv) {
  int start = 0, limit = 1;
  bool gate_audit = false, list = false;
  string mode = "compare";
  for (int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if (arg == "--context-start" && i + 1 < argc)
      start = stoi(argv[++i]);
    else if (arg == "--context-limit" && i + 1 < argc)
      limit = stoi(argv[++i]);
    else if (arg == "--mode" && i + 1 < argc)
      mode = argv[++i];
    else if (arg == "--gate-audit")
      gate_audit = true;
    else if (arg == "--literal-crosscheck")
      literal_crosscheck = true;
    else if (arg == "--list-contexts")
      list = true;
    else
      throw runtime_error("unknown or incomplete argument: " + arg);
  }
  unsigned long long checks = 0;
  vector<Context> contexts = reconstruct_contexts(checks);
  if (checks != 94448 || contexts.size() != 1504)
    throw runtime_error("context reconstruction mismatch");
  int end = min<int>(contexts.size(), start + limit);
  if (start < 0 || limit < 0 || start > end)
    throw runtime_error("invalid context range");
  if (list) {
    for (int i = start; i < end; ++i)
      cout << "CONTEXT|" << i << "|" << type_word(stratum(contexts[i]))
           << "|" << contexts[i].root_mask << "|"
           << list_word(contexts[i].labels, ',') << "|"
           << list_word(contexts[i].orders, ',') << "|"
           << list_word(contexts[i].units, ',') << "\n";
    return 0;
  }
  cout << "SCALE_THREE_HAMMING_SIX_DEPTH_THREE_GEOMETRY_BATCH\n"
       << "context_start=" << start << " context_end=" << end
       << " contexts=" << end - start << " gate_audit=" << gate_audit
       << " arithmetic=integer+rational floating_point=none\n";
  optional<Audit> raw, batch;
  if (mode == "raw" || mode == "compare") {
    safe_band_cache.clear();
    raw = run_raw(contexts, start, end, gate_audit);
    print_audit("RAW", *raw);
  }
  if (mode == "batch" || mode == "compare") {
    safe_band_cache.clear();
    batch = run_batched(contexts, start, end, gate_audit);
    print_audit("BATCH", *batch);
  }
  if (mode != "raw" && mode != "batch" && mode != "compare")
    throw runtime_error("mode must be raw, batch, or compare");
  if (raw && batch) {
    if (raw->rows != batch->rows ||
        raw->fingerprint_xor != batch->fingerprint_xor ||
        raw->fingerprint_sum != batch->fingerprint_sum)
      throw runtime_error("raw/batch logical-lane mismatch");
    cout << "COMPARE|rows_equal=1|fingerprints_equal=1|meet3_reduction="
         << raw->meet3 << "/" << batch->meet3 << "|wall_speedup="
         << raw->seconds / batch->seconds << "|cpu_speedup="
         << raw->cpu_seconds / batch->cpu_seconds << "\n";
  }
  cout << "TOURNAMENT|vertices=remaining_labelled_rays"
          "|pair_observable=least_legal_speed_difference"
          "|switch=numerical_order|ties=increasing_label"
          "|triangles=0|SCCs=singletons|Hamiltonian_paths=1\n"
       << "ASSUMPTION_AUDIT|cache_vertices=metric_prefix_geometries"
          "|proof_vertices=labelled_future_ray_obligations"
          "|preserves=literal_residual_components"
          "|destroys=label_availability,least_future,ray_tournament,shortcut_ancestry\n"
       << "DONE\n";
}
