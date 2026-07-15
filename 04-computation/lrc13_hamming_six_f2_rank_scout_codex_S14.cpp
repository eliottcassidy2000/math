// Exact scout for the 360 primitive-core scale-one H6 rows with two full
// antipodal missing pairs.
//
// This deliberately reuses the hash-frozen interval arithmetic and recursive
// discrepancy carrier from the completed six-row f=3 certificate.  It does
// not alter that certificate.  `--rank` only ranks root states by their exact
// first-speed cap.  `--row i` exhausts one lexicographically indexed f=2 row.
// A row is proved loose exactly when the exhaustive tree has no empty residual
// state; nonempty depth-six leaves are loose completions, not failures.
//
// Tournament guardrail: increasing speed order is merely a transitive
// enumeration tournament.  The proof state is the literal residual endpoint
// union together with the labelled bank of unused residue classes.  The f=2
// label is a routing stratum and forgets both endpoint geometry and heights.

#define main lrc13_h6_f3_frozen_certificate_main
#include "lrc13_hamming_six_open_f3_exact_closure_codex_S14.cpp"
#undef main

#include <cstdlib>
#include <tuple>

int full_pair_count(array<int, 6> const &missing) {
  int count = 0;
  for (int r = 1; r <= 6; ++r) {
    const bool left = find(missing.begin(), missing.end(), r) != missing.end();
    const bool right =
        find(missing.begin(), missing.end(), 13 - r) != missing.end();
    if (left && right)
      ++count;
  }
  return count;
}

vector<RowSpec> f2_specs() {
  vector<RowSpec> specs;
  for (int a = 1; a <= 7; ++a)
    for (int b = a + 1; b <= 8; ++b)
      for (int c = b + 1; c <= 9; ++c)
        for (int d = c + 1; d <= 10; ++d)
          for (int e = d + 1; e <= 11; ++e)
            for (int f = e + 1; f <= 12; ++f) {
              array<int, 6> missing = {a, b, c, d, e, f};
              if (full_pair_count(missing) != 2)
                continue;
              array<int, 6> core{};
              int at = 0;
              for (int r = 1; r <= 12; ++r)
                if (find(missing.begin(), missing.end(), r) == missing.end())
                  core[at++] = r;
              int g = 0;
              for (int r : core)
                g = gcd(g, r);
              if (g != 1)
                continue;
              RowSpec spec{};
              spec.missing = missing;
              spec.core = core;
              specs.push_back(spec);
            }
  if (specs.size() != 360)
    throw runtime_error("f=2 atlas size mismatch");
  return specs;
}

struct RootRank {
  size_t index;
  array<int, 6> missing;
  array<int, 6> core;
  size_t components;
  Rat longest;
  long long cap;
};

RootRank rank_root(size_t index, RowSpec const &spec) {
  vector<Interval> root = {{Rat(0), Rat(1)}};
  for (int u : spec.core)
    root = meet(root, safe_bands(u));
  if (root.empty())
    throw runtime_error("primitive f=2 root unexpectedly empty");
  const Rat longest = longest_length(root);
  return {index, spec.missing, spec.core, root.size(), longest,
          discrepancy_cap(6, longest)};
}

void print_rank(RootRank const &rank) {
  cout << "index=" << rank.index << " missing=" << set_text(rank.missing)
       << " core=" << set_text(rank.core)
       << " components=" << rank.components << " longest=" << rank.longest
       << " first_cap=" << rank.cap << "\n";
}

struct FrozenF2 {
  size_t index;
  array<unsigned long long, 7> nodes;
  uint64_t trace_a;
  uint64_t trace_b;
  size_t cached_speeds;
  unsigned long long crosschecks;
};

vector<FrozenF2> cap156_frozen() {
  return {
      {77, {1, 45, 2557, 77709, 3702, 2, 0},
       0x85a7697025e0e5aeULL, 0x6f3287bbc440ffbdULL, 275, 3},
      {5, {1, 56, 3674, 132834, 5252, 3, 0},
       0xea40c77864c66fbfULL, 0x5cb78e611942473aULL, 324, 4},
      {47, {1, 56, 3738, 136081, 6709, 36, 0},
       0x5c9b96273e3f43eaULL, 0xd997469908694974ULL, 325, 37},
      {62, {1, 56, 3752, 136711, 4005, 6, 0},
       0xd3e41a35c09d3057ULL, 0x49fc1e008c635594ULL, 325, 7},
      {71, {1, 56, 3756, 137387, 8745, 47, 0},
       0x9027fd6e6ad82d10ULL, 0xd66512e9019adf18ULL, 325, 48},
      {74, {1, 56, 3759, 137400, 4436, 20, 0},
       0xa42fc5f303eac307ULL, 0x9e615106e826fab6ULL, 325, 21},
      {76, {1, 56, 3762, 137747, 5283, 17, 0},
       0xf778c269845c4c6fULL, 0x08004eaedd6a6215ULL, 324, 18},
      {128, {1, 55, 3657, 132503, 5398, 20, 0},
       0x0c469ed636149c19ULL, 0x6f063ff0e3559e12ULL, 318, 21},
      {158, {1, 55, 3660, 132877, 8694, 74, 0},
       0xbd89169f8e179c8aULL, 0x92a078ef520bcbb0ULL, 318, 75},
      {331, {1, 54, 3546, 126277, 3445, 2, 0},
       0xe4cfc62613369d5dULL, 0xd69c82a9836ad2a4ULL, 313, 3},
      {346, {1, 54, 3554, 126862, 6774, 2, 0},
       0x906297eda7bb524dULL, 0xa6d6f73cad95d006ULL, 313, 3},
      {11, {1, 58, 3939, 146701, 6052, 18, 0},
       0x9dfd2ad13c83ba30ULL, 0x919bd28ef8c9c60fULL, 328, 19},
      {17, {1, 58, 3940, 147191, 5652, 7, 0},
       0x0c1326e1ea529a78ULL, 0xafdb1d35524a5d1eULL, 332, 8},
      {20, {1, 58, 3940, 147437, 4263, 8, 0},
       0xfc71271ea0f1f86eULL, 0xade43866800e064cULL, 334, 9},
      {23, {1, 58, 3949, 148058, 5702, 17, 0},
       0x0f2a3345d11aa555ULL, 0x02f78834446a9fd3ULL, 337, 18},
      {26, {1, 57, 3837, 141812, 5259, 17, 0},
       0x783b6f2c0dbb3a6fULL, 0xe97a7d3453ca485aULL, 326, 18},
      {359, {1, 58, 4055, 154487, 4731, 6, 0},
       0x149302a47c26ba62ULL, 0x87dd49a94552dad2ULL, 338, 7},
      {257, {1, 59, 4143, 159295, 5450, 26, 0},
       0x8f5cdd713bf2b66cULL, 0x7dee095ebd180dc0ULL, 346, 27},
      {326, {1, 60, 4314, 169208, 8247, 11, 0},
       0xef3fb9706f7de64cULL, 0xa293bb2d2b1d4f81ULL, 347, 12},
      {319, {1, 62, 4537, 181094, 8579, 3, 0},
       0xb79c7408b53fa573ULL, 0x6269c1b262ea81cbULL, 361, 4},
      {2, {1, 64, 4681, 190234, 8541, 45, 0},
       0x2461b4d61c03ae01ULL, 0x0d8ce890da83ae25ULL, 363, 46},
      {4, {1, 64, 4679, 190329, 9644, 69, 0},
       0xdd0a6615ee9a6657ULL, 0x119d43cd4defdac0ULL, 358, 70},
      {174, {1, 64, 4807, 197729, 3803, 4, 0},
       0xf5cec418dc344162ULL, 0xd2a3068d62ab79e6ULL, 369, 5},
      {350, {1, 64, 4822, 199796, 5908, 13, 0},
       0xd224adaae473b4d9ULL, 0xb43b437e109aa1abULL, 369, 14},
      {356, {1, 64, 4835, 200477, 5978, 26, 0},
       0xb9ebd772ba89f77eULL, 0x5f346dc68324f4a4ULL, 369, 27},
      {347, {1, 63, 4723, 193483, 5915, 11, 0},
       0x7d4435af92a8eddfULL, 0x7b3141583f557372ULL, 369, 12},
      {32, {1, 64, 4741, 192868, 8145, 51, 0},
       0x54733cb19949231fULL, 0xaf8ad59641a147efULL, 363, 52},
      {38, {1, 64, 4748, 193673, 7157, 38, 0},
       0x6ec12d791a27d906ULL, 0x491e26998739811eULL, 364, 39},
      {41, {1, 64, 4748, 193824, 5543, 9, 0},
       0xf294e5a1c201c8b1ULL, 0xdaaad8e2111d108eULL, 366, 10},
      {44, {1, 64, 4750, 194566, 7047, 25, 0},
       0x21dee79e6b77313dULL, 0x47653f0f42fcad06ULL, 366, 26},
      {53, {1, 64, 4758, 194399, 9721, 49, 0},
       0x66f53cb60fe8b263ULL, 0x49a6580dab0f0a26ULL, 366, 50},
      {56, {1, 64, 4761, 194706, 5285, 13, 0},
       0x68f56be8f35e09eeULL, 0xf6bab662dbbdc9cfULL, 368, 14},
      {61, {1, 64, 4779, 196178, 6571, 40, 0},
       0xa27bd23536f2a67cULL, 0x0a1d3f5d3fff5951ULL, 372, 41},
      {69, {1, 64, 4770, 195861, 7678, 35, 0},
       0xbeac26a5449d2b91ULL, 0xf25a3b412f1b6946ULL, 368, 36},
      {70, {1, 64, 4788, 197184, 9139, 35, 0},
       0x768981bd81722bd3ULL, 0xcb5121b0352a9428ULL, 373, 36},
      {72, {1, 64, 4770, 196065, 5117, 18, 0},
       0xb4da2a785d7997b5ULL, 0x4ea20ec565059ea9ULL, 368, 19},
      {73, {1, 64, 4776, 196804, 4816, 8, 0},
       0x9afbbffb8b5cc166ULL, 0xb1efe2a301fd439fULL, 373, 9},
      {75, {1, 64, 4797, 197940, 7663, 25, 0},
       0x0bbe09e0be342e03ULL, 0xedf9ae9d768206a3ULL, 373, 26},
      {280, {1, 66, 5105, 217388, 7170, 33, 0},
       0xfd04e1349e2d92cfULL, 0xe9c582827548a950ULL, 380, 34},
      {169, {1, 66, 5067, 214895, 6402, 9, 0},
       0x971380b5f27286a9ULL, 0x40a5476e19569373ULL, 379, 10},
      {323, {1, 66, 5085, 216479, 10688, 49, 0},
       0xc40003e1d77cbb5fULL, 0x28ab07641cc5bc57ULL, 379, 50},
      {358, {1, 66, 5111, 218092, 8461, 20, 0},
       0xe3e723168dca39b9ULL, 0xf61075a1bf29f5e9ULL, 377, 21},
      {301, {1, 66, 5056, 213427, 6486, 8, 0},
       0x36b7724faf172d25ULL, 0xa9506193a69d64d4ULL, 377, 9},
      {8, {1, 66, 4972, 207731, 8961, 60, 0},
       0xa6e11d3f37be1315ULL, 0x04990eb8a65b609bULL, 380, 61},
      {10, {1, 66, 4986, 208561, 6260, 18, 0},
       0xa29396a7079685b7ULL, 0x29804641f102683bULL, 377, 19},
      {14, {1, 66, 4971, 207086, 5831, 5, 0},
       0x8c8cb1841443be4bULL, 0x1bcd228bd20099fbULL, 379, 6},
      {16, {1, 66, 4993, 209280, 12632, 31, 0},
       0x3267a59299de2fc8ULL, 0x35a603f5c84b621fULL, 378, 32},
      {19, {1, 66, 5002, 209839, 6914, 24, 0},
       0xf7b5ad54ddbd435fULL, 0xc08d4f29560a8c5dULL, 380, 25},
      {22, {1, 66, 5001, 210461, 9064, 68, 0},
       0x6330492aa915060cULL, 0x0580784a29b4448fULL, 380, 69},
      {25, {1, 66, 5015, 210978, 9470, 65, 0},
       0xcd04172987ef0582ULL, 0x3f440d60cc60dc45ULL, 380, 66},
      {176, {1, 66, 5090, 215933, 4217, 5, 0},
       0xd337b2978bf7917bULL, 0x82aed333ceb1088dULL, 379, 6},
      {278, {1, 66, 5095, 216460, 6281, 20, 0},
       0xb0e1d3849b12d9e9ULL, 0xb9b5467dd49d3882ULL, 379, 21},
  };
}

unique_ptr<Replay> run_frozen_f2(RowSpec const &spec, FrozenF2 const &frozen) {
  auto replay = make_unique<Replay>(spec);
  replay->run();
  const array<unsigned long long, 7> zero_cover{};
  for (int depth = 0; depth < 6; ++depth)
    if (replay->edges[depth] != replay->nodes[depth + 1])
      throw runtime_error("f=2 tree edge/node mismatch");
  if (replay->nodes != frozen.nodes || replay->covering != zero_cover ||
      replay->trace.a != frozen.trace_a || replay->trace.b != frozen.trace_b ||
      replay->cache.size() != frozen.cached_speeds ||
      replay->root_crosschecks + replay->deepest_crosschecks !=
          frozen.crosschecks)
    throw runtime_error("f=2 frozen row certificate mismatch");
  return replay;
}

int scout_main(int argc, char **argv) {
  if (argc != 2 && argc != 3)
    throw runtime_error("usage: --rank [N] | --row INDEX | --cap156");
  const vector<RowSpec> specs = f2_specs();

  cout << "LRC13 PRIMITIVE-CORE H6 F=2 EXACT RANK/ROW SCOUT\n";
  cout << "dependency=lrc13_hamming_six_open_f3_exact_closure_codex_S14.cpp\n";
  cout << "dependency_sha256=af8fd1cc6b4eaff076c2cc25d20d7613efde5802ce4a81e7b545df19d8314d87\n";
  cout << "atlas_rows=" << specs.size()
       << " faithful_carrier=residual_endpoint_union_x_labelled_comb_bank\n";

  const string mode = argv[1];
  if (mode == "--rank") {
    size_t limit = specs.size();
    if (argc == 3)
      limit = min(limit, (size_t)stoull(argv[2]));
    vector<RootRank> ranks;
    for (size_t i = 0; i < specs.size(); ++i)
      ranks.push_back(rank_root(i, specs[i]));
    sort(ranks.begin(), ranks.end(), [](RootRank const &x, RootRank const &y) {
      return tie(x.cap, x.longest.d, x.longest.n, x.missing) <
             tie(y.cap, y.longest.d, y.longest.n, y.missing);
    });
    cout << "root_first_cap_range=" << ranks.front().cap << ","
         << ranks.back().cap << " requested=" << limit << "\n";
    for (size_t i = 0; i < limit; ++i)
      print_rank(ranks[i]);
    cout << "RANK ONLY: no row closure claimed\n";
    return 0;
  }

  if (mode == "--cap156" && argc == 2) {
    const vector<FrozenF2> frozen = cap156_frozen();
    vector<size_t> cap156_indices;
    for (size_t i = 0; i < specs.size(); ++i)
      if (rank_root(i, specs[i]).cap <= 156)
        cap156_indices.push_back(i);
    vector<size_t> expected_indices;
    for (FrozenF2 const &row : frozen)
      expected_indices.push_back(row.index);
    sort(expected_indices.begin(), expected_indices.end());
    vector<size_t> sorted_actual = cap156_indices;
    sort(sorted_actual.begin(), sorted_actual.end());
    if (sorted_actual != expected_indices)
      throw runtime_error("cap<=156 sub-stratum mismatch");

    vector<future<unique_ptr<Replay>>> jobs;
    for (FrozenF2 const &row : frozen) {
      RowSpec const *spec = &specs[row.index];
      FrozenF2 const *certificate = &row;
      jobs.push_back(async(launch::async, [spec, certificate]() {
        return run_frozen_f2(*spec, *certificate);
      }));
    }

    array<unsigned long long, 7> aggregate_nodes{};
    unsigned long long aggregate_crosschecks = 0;
    for (size_t i = 0; i < frozen.size(); ++i) {
      unique_ptr<Replay> replay = jobs[i].get();
      FrozenF2 const &row = frozen[i];
      RootRank root = rank_root(row.index, replay->spec);
      for (int depth = 0; depth < 7; ++depth)
        aggregate_nodes[depth] += replay->nodes[depth];
      aggregate_crosschecks +=
          replay->root_crosschecks + replay->deepest_crosschecks;
      cout << "row_index=" << row.index
           << " missing=" << set_text(replay->spec.missing)
           << " first_cap=" << root.cap
           << " nodes=" << counts_text(replay->nodes)
           << " trace128=" << hex << setfill('0') << setw(16)
           << replay->trace.a << setw(16) << replay->trace.b << dec
           << " crosschecks="
           << replay->root_crosschecks + replay->deepest_crosschecks << "\n";
    }
    const array<unsigned long long, 7> expected_nodes = {
        52, 3202, 232351, 9302397, 348886, 1271, 0};
    if (aggregate_nodes != expected_nodes || aggregate_crosschecks != 1323)
      throw runtime_error("cap<=156 aggregate mismatch");
    const unsigned long long total_states =
        accumulate(aggregate_nodes.begin(), aggregate_nodes.end(), 0ULL);
    if (total_states != 9888159)
      throw runtime_error("cap<=156 state total mismatch");
    cout << "aggregate_nodes=" << counts_text(aggregate_nodes)
         << " total_states=" << total_states
         << " independent_crosschecks=" << aggregate_crosschecks << "\n";
    cout << "PROVED LOOSE: all 52 primitive-core f=2 rows with root first_cap<=156\n";
    cout << "INDEPENDENT SUBCERTIFICATE: 52 of 360 f=2 roots; THM-857 subsequently closes the full proper scale-one H6 chart\n";
    return 0;
  }

  if (mode != "--row" || argc != 3)
    throw runtime_error("usage: --rank [N] | --row INDEX | --cap156");
  const size_t index = stoull(argv[2]);
  if (index >= specs.size())
    throw runtime_error("row index out of range");
  RowSpec const &spec = specs[index];
  print_rank(rank_root(index, spec));
  Replay replay(spec);
  replay.run();
  const unsigned long long total_states =
      accumulate(replay.nodes.begin(), replay.nodes.end(), 0ULL);
  const unsigned long long total_covering =
      accumulate(replay.covering.begin(), replay.covering.end(), 0ULL);
  cout << "nodes=" << counts_text(replay.nodes)
       << " edges=" << counts_text(replay.edges)
       << " dead=" << counts_text(replay.dead)
       << " covering=" << counts_text(replay.covering)
       << " total_states=" << total_states << "\n";
  cout << "trace128=" << hex << setfill('0') << setw(16) << replay.trace.a
       << setw(16) << replay.trace.b << dec
       << " cached_speeds=" << replay.cache.size()
       << " independent_crosschecks="
       << replay.root_crosschecks + replay.deepest_crosschecks << "\n";
  cout << "tournament=speed_order_transitive_enumeration_only"
       << " score_histogram=0,1,2,3,4,5 cycles=0 scc=1,1,1,1,1,1 hp=1\n";
  if (total_covering == 0)
    cout << "PROVED LOOSE ROW index=" << index
         << " missing=" << set_text(spec.missing) << "\n";
  else
    cout << "TIGHT COMPLETIONS FOUND=" << total_covering
         << " row_not_closed index=" << index << "\n";
  return 0;
}

int main(int argc, char **argv) {
  try {
    return scout_main(argc, argv);
  } catch (exception const &error) {
    cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }
}
