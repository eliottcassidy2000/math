#define main lrc_two_mask_concepts_hidden_main
#include "../lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_detached_hostile_survivor.cpp"
#undef main

#include <map>

namespace {

constexpr u32 kFull = (u32{1} << 30) - 1;
constexpr std::array<int, 3> kQ = {96, 100, 210};

// Unlike the imported rank-8/9 hostile evaluator, retain every literal-wall
// failure class.  This makes margin() exact for masks of arbitrary rank.
Geometry build_full_geometry(int pair_q, int pair_r) {
  i64 grid = fixed_grid();
  grid = checked_lcm(grid, 14LL * pair_q);
  grid = checked_lcm(grid, 14LL * pair_r);
  std::vector<i64> walls = {0, grid};
  auto add_walls = [&](int speed) {
    const i64 divisor = 14LL * speed;
    require(grid % divisor == 0, "nonintegral wall unit");
    const i64 unit = grid / divisor;
    for (int tooth = 0; tooth < speed; ++tooth) {
      walls.push_back((14LL * tooth + 1) * unit);
      walls.push_back((14LL * tooth + 13) * unit);
    }
  };
  for (int speed : kPool) add_walls(speed);
  add_walls(pair_q);
  add_walls(pair_r);
  std::sort(walls.begin(), walls.end());
  walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
  std::map<u32, i64> by_failure;
  Geometry geometry;
  geometry.grid = grid;
  geometry.cells = walls.size() - 1;
  for (std::size_t index = 1; index < walls.size(); ++index) {
    const i64 left = walls[index - 1];
    const i64 right = walls[index];
    if (!safe_midpoint(pair_q, grid, left, right) ||
        !safe_midpoint(pair_r, grid, left, right))
      continue;
    u32 failure = 0;
    for (unsigned vertex = 0; vertex < kPool.size(); ++vertex)
      if (!safe_midpoint(kPool[vertex], grid, left, right))
        failure |= u32{1} << vertex;
    const i64 width = right - left;
    geometry.pair_ticks += width;
    by_failure[failure] += width;
  }
  for (const auto& entry : by_failure) geometry.classes.push_back(entry);
  return geometry;
}

std::array<std::vector<u32>, 3> read_tagged_failures(
    const std::filesystem::path& path) {
  std::ifstream input(path);
  require(static_cast<bool>(input), "cannot open failure ledger");
  std::string line;
  require(std::getline(input, line) && line == "q,r,body_hex",
          "failure header changed");
  std::array<std::vector<u32>, 3> out;
  while (std::getline(input, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    int q = 0, r = 0;
    std::string token;
    row >> q >> r >> token;
    require(row && r == 595, "bad failure row");
    const int p = q == 96 ? 0 : q == 100 ? 1 : q == 210 ? 2 : -1;
    require(p >= 0, "unexpected pair");
    const u32 body = static_cast<u32>(std::stoul(token, nullptr, 16));
    require(std::popcount(body) == 9, "body rank changed");
    out[p].push_back(body);
  }
  require(out[0].size() == 116 && out[1].size() == 13 &&
              out[2].size() == 16,
          "failure counts changed");
  return out;
}

u32 family_union(const std::vector<u32>& family) {
  u32 answer = 0;
  for (u32 body : family) answer |= body;
  return answer;
}

std::array<bool, 3> active_signature(
    u32 mask, const std::array<Geometry, 3>& geometry) {
  std::array<bool, 3> answer{};
  for (int p = 0; p < 3; ++p)
    answer[p] = margin(geometry[p], mask).ticks >= 0;
  return answer;
}

std::array<unsigned, 3> coverage(
    u32 mask, const std::array<std::vector<u32>, 3>& failures,
    const std::array<Geometry, 3>& geometry) {
  const auto active = active_signature(mask, geometry);
  std::array<unsigned, 3> answer{};
  for (int p = 0; p < 3; ++p)
    if (active[p])
      for (u32 body : failures[p]) answer[p] += (mask & body) == 0;
  return answer;
}

bool covers_all(u32 a, u32 b,
                const std::array<std::vector<u32>, 3>& failures,
                const std::array<Geometry, 3>& geometry) {
  const auto aa = active_signature(a, geometry);
  const auto ab = active_signature(b, geometry);
  for (int p = 0; p < 3; ++p)
    for (u32 body : failures[p])
      if (!((aa[p] && !(a & body)) || (ab[p] && !(b & body)))) return false;
  return true;
}

// derivation(X) is the common-neighbor set in the no-cooccurrence relation.
u32 derivation(u32 x, const std::array<u32, 30>& relation) {
  u32 answer = kFull;
  while (x) {
    const unsigned bit = std::countr_zero(x);
    x &= x - 1;
    answer &= relation[bit];
  }
  return answer;
}

u32 closure(u32 x, const std::array<u32, 30>& relation) {
  return derivation(derivation(x, relation), relation);
}

std::string bits(bool x, bool y, bool z) {
  return std::string{x ? '1' : '0', y ? '1' : '0', z ? '1' : '0'};
}

}  // namespace

#ifndef LRC_TWO_MASK_NO_MAIN
int main(int argc, char** argv) {
  try {
    require(argc == 2, "usage: two_mask_concepts FAILURE_CSV");
    const auto failures = read_tagged_failures(argv[1]);
    const std::array<Geometry, 3> geometry = {
        build_full_geometry(96, 595), build_full_geometry(100, 595),
        build_full_geometry(210, 595)};

    std::array<u32, 30> cooccurrence{};
    u32 all_union = 0;
    for (const auto& family : failures)
      for (u32 body : family) {
        all_union |= body;
        u32 rest = body;
        while (rest) {
          const unsigned bit = std::countr_zero(rest);
          rest &= rest - 1;
          cooccurrence[bit] |= body;
        }
      }
    std::array<u32, 30> relation{};
    for (unsigned bit = 0; bit < 30; ++bit)
      relation[bit] = kFull ^ cooccurrence[bit];

    std::cout << "LRC14_TWO_MASK_FORMAL_CONCEPT_AUDIT_V1\n";
    std::cout << "FULL_CLASSES " << geometry[0].classes.size() << ','
              << geometry[1].classes.size() << ','
              << geometry[2].classes.size() << '\n';
    std::cout << "ALL_UNION " << hex8(all_union) << " ABSENT "
              << hex8(kFull ^ all_union) << '\n';
    for (int p = 0; p < 3; ++p) {
      const u32 allowed = kFull ^ family_union(failures[p]);
      const auto sig = active_signature(allowed, geometry);
      std::cout << "WHOLE_MAX PAIR " << kQ[p] << ",595 MASK "
                << hex8(allowed) << " RANK " << std::popcount(allowed)
                << " ACTIVE_SIGNATURE " << bits(sig[0], sig[1], sig[2])
                << " OWN_ACTIVE " << sig[p] << " OWN_TICKS "
                << decimal(margin(geometry[p], allowed).ticks) << '\n';
    }

    u64 concepts = 0, active_96_210 = 0, active_all = 0, exact = 0;
    Fnv concept_ledger, active_ledger;
    u32 best_a = 0, best_b = 0;
    unsigned best_total = 0;
    std::array<unsigned, 3> best_union{};

    u32 a = closure(0, relation);
    while (true) {
      const u32 b = derivation(a, relation);
      require(closure(a, relation) == a && derivation(b, relation) == a,
              "formal concept invariant failed");
      ++concepts;
      concept_ledger.add(a);
      concept_ledger.add(b);
      const auto sa = active_signature(a, geometry);
      const auto sb = active_signature(b, geometry);
      if (sa[0] && sa[2] && sb[0] && sb[2]) {
        ++active_96_210;
        active_ledger.add(a);
        active_ledger.add(b);
        if (sa[1] && sb[1]) ++active_all;
        const auto ca = coverage(a, failures, geometry);
        const auto cb = coverage(b, failures, geometry);
        std::array<unsigned, 3> joined{};
        unsigned total = 0;
        for (int p = 0; p < 3; ++p) {
          for (u32 body : failures[p])
            joined[p] += ((sa[p] && !(a & body)) ||
                          (sb[p] && !(b & body)));
          total += joined[p];
        }
        if (total > best_total ||
            (total == best_total && std::pair{a, b} < std::pair{best_a, best_b})) {
          best_total = total;
          best_a = a;
          best_b = b;
          best_union = joined;
        }
        if (covers_all(a, b, failures, geometry)) {
          ++exact;
          std::cout << "EXACT_CONCEPT A " << hex8(a) << " RANK "
                    << std::popcount(a) << " SIG "
                    << bits(sa[0], sa[1], sa[2]) << " B " << hex8(b)
                    << " RANK " << std::popcount(b) << " SIG "
                    << bits(sb[0], sb[1], sb[2]) << '\n';
        }
      }
      if (a == kFull) break;
      bool advanced = false;
      for (int i = 29; i >= 0; --i) {
        const u32 bit = u32{1} << i;
        if (a & bit) continue;
        const u32 prefix = a & (bit - 1);
        const u32 candidate = closure(prefix | bit, relation);
        if ((candidate & (bit - 1)) == prefix) {
          a = candidate;
          advanced = true;
          break;
        }
      }
      require(advanced, "NextClosure failed before full set");
    }

    // If exactly one mask is active at q=100, the active mask must be
    // disjoint from the whole q=100 family.  Enumerate its complete 2^15
    // support universe and maximize the other side by derivation.
    const u32 allowed100 = kFull ^ family_union(failures[1]);
    u64 one_sided_tested = 0, one_sided_candidates = 0,
        one_sided_exact = 0;
    Fnv one_sided_ledger;
    for (u32 x = allowed100;; x = (x - 1) & allowed100) {
      ++one_sided_tested;
      const auto sx = active_signature(x, geometry);
      if (sx[0] && sx[1] && sx[2]) {
        const u32 y = derivation(x, relation);
        const auto sy = active_signature(y, geometry);
        if (sy[0] && sy[2]) {
          ++one_sided_candidates;
          one_sided_ledger.add(x);
          one_sided_ledger.add(y);
          if (covers_all(x, y, failures, geometry)) {
            ++one_sided_exact;
            std::cout << "EXACT_ONE_SIDED A " << hex8(x) << " RANK "
                      << std::popcount(x) << " B " << hex8(y) << " RANK "
                      << std::popcount(y) << " B_SIG "
                      << bits(sy[0], sy[1], sy[2]) << '\n';
          }
        }
      }
      if (x == 0) break;
    }

    std::cout << "CONCEPTS " << concepts << " FNV " << std::hex
              << concept_ledger.state << std::dec << " ACTIVE_96_210 "
              << active_96_210 << " ACTIVE_FNV " << std::hex
              << active_ledger.state << std::dec << " ACTIVE_ALL "
              << active_all << " EXACT " << exact << '\n';
    std::cout << "BEST_CONCEPT A " << hex8(best_a) << " B " << hex8(best_b)
              << " COVER " << best_union[0] << ',' << best_union[1] << ','
              << best_union[2] << " TOTAL " << best_total << '\n';
    std::cout << "ONE_SIDED_TESTED " << one_sided_tested
              << " CANDIDATES " << one_sided_candidates << " FNV "
              << std::hex << one_sided_ledger.state << std::dec << " EXACT "
              << one_sided_exact << '\n';
    require(exact == 0 && one_sided_exact == 0,
            "two-mask cover unexpectedly found");
    std::cout << "SCOPE RANK_FREE_FULL_LITERAL_CLASSES_EXACT_TWO_MASK_OBSTRUCTION_ON_FIXED_145_FAILURES_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14\n"
              << "VERDICT PASS_NO_TWO_MASK_COVER\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "TWO_MASK_CONCEPT_ERROR " << error.what() << '\n';
    return 1;
  }
}
#endif
