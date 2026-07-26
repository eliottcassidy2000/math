#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

namespace {

constexpr int L = 91;

using Mask = bitset<L>;
using Signature = vector<int>;
using SignatureSet = set<Signature>;

array<uint64_t, 2> key91(const Mask& a) {
  array<uint64_t, 2> key{0, 0};
  for (int x = 0; x < L; ++x) {
    if (a[x]) key[x / 64] |= uint64_t{1} << (x % 64);
  }
  return key;
}

void require(bool condition, const string& message) {
  if (!condition) throw runtime_error(message);
}

struct Universe {
  Mask guard;
  vector<Mask> ap;
  vector<int> step;
  vector<int> center;
  vector<vector<int>> through;
  array<Mask, 13> row;
  map<array<uint64_t, 2>, int> id;

  Universe() : through(L) {
    for (int x = 0; x < 26; ++x) guard[x] = 1;
    for (int h = 0; h < 13; ++h) {
      for (int x = h; x < L; x += 13) row[h][x] = 1;
    }

    for (int d = 1; d < L; ++d) {
      if (gcd(d, L) != 1) continue;
      for (int s = 0; s < L; ++s) {
        Mask a;
        for (int j = 0; j < 13; ++j) a[(s + j * d) % L] = 1;
        const auto key = key91(a);
        if (id.count(key)) continue;
        id[key] = static_cast<int>(ap.size());
        ap.push_back(a);
        step.push_back(min(d, L - d));
        center.push_back((s + 6 * d) % L);
      }
    }

    require(ap.size() == 3276, "unit AP13 universe is not 3276");
    for (int i = 0; i < static_cast<int>(ap.size()); ++i) {
      for (int x = 0; x < L; ++x) {
        if (ap[i][x]) through[x].push_back(i);
      }
    }
  }

  int reflected_id(int i) const {
    Mask reflected;
    for (int x = 0; x < L; ++x) {
      if (ap[i][x]) reflected[(25 - x + L) % L] = 1;
    }
    return id.at(key91(reflected));
  }
};

struct Atlas {
  SignatureSet solutions;
  size_t flow_words = 0;
  size_t step_spectra = 0;
  uint64_t flat = 0;
  uint64_t bad = 0;
  uint64_t repeated_support = 0;
  uint64_t no_repeated_step = 0;
  uint64_t h_independent_nonzero = 0;
  uint64_t repeated_duplicate_target = 0;
  uint64_t any_double_row = 0;
  uint64_t max_double_rows = 0;
  uint64_t max_base_multiplicity = 0;
  uint64_t affine_arrow_words = 0;
  uint64_t chi4_failures = 0;
  map<int, uint64_t> total_holes_hist;
  map<int, uint64_t> arrow_hist;
  map<int, set<int>> center_bank;
  map<unsigned, uint64_t> exact_hole_residue_hist;
};

unsigned source_mask(int s0, int s1) {
  return (1u << s0) | (1u << s1);
}

unsigned hole_residue_mask(const Universe& u, const Signature& v) {
  unsigned answer = 0;
  for (int h = 0; h < 13; ++h) {
    array<int, 7> m{};
    for (int x = h; x < L; x += 13) {
      if (u.guard[x]) ++m[x % 7];
    }
    for (int i : v) {
      for (int x = h; x < L; x += 13) {
        if (u.ap[i][x]) ++m[x % 7];
      }
    }
    for (int r = 0; r < 7; ++r) {
      if (m[r] == 0) answer |= 1u << r;
    }
  }
  return answer;
}

Atlas enumerate_atlas(const Universe& u, int s0, int s1) {
  Mask blockers;
  for (int x = 0; x < L; ++x) {
    if (x % 7 == s0 || x % 7 == s1) blockers[x] = 1;
  }
  const Mask required = ~(u.guard | blockers);

  auto feasible = [&](const Mask& covered, int remaining) {
    const Mask missing = required & ~covered;
    if (missing.count() > static_cast<size_t>(13 * remaining)) return false;
    for (int h = 0; h < 13; ++h) {
      if ((missing & u.row[h]).count() > static_cast<size_t>(remaining)) {
        return false;
      }
    }
    return true;
  };

  Signature chosen;
  SignatureSet solutions;
  function<void(const Mask&)> search = [&](const Mask& covered) {
    const int depth = static_cast<int>(chosen.size());
    const int remaining = 5 - depth;
    if (!feasible(covered, remaining)) return;
    const Mask missing = required & ~covered;

    if (depth == 5) {
      if (missing.none()) {
        Signature signature = chosen;
        sort(signature.begin(), signature.end());
        solutions.insert(signature);
      }
      return;
    }

    if (depth == 4) {
      if (missing.none()) {
        for (int i = 0; i < static_cast<int>(u.ap.size()); ++i) {
          Signature signature = chosen;
          signature.push_back(i);
          sort(signature.begin(), signature.end());
          solutions.insert(signature);
        }
      } else {
        int pivot = -1;
        for (int x = 0; x < L && pivot < 0; ++x) {
          if (missing[x]) pivot = x;
        }
        for (int i : u.through[pivot]) {
          if ((missing & ~u.ap[i]).none()) {
            Signature signature = chosen;
            signature.push_back(i);
            sort(signature.begin(), signature.end());
            solutions.insert(signature);
          }
        }
      }
      return;
    }

    int pivot = -1;
    for (int x = 0; x < L && pivot < 0; ++x) {
      if (missing[x]) pivot = x;
    }
    require(pivot >= 0, "premature complete cover before depth four");
    for (int i : u.through[pivot]) {
      const Mask next = covered | u.ap[i];
      if (!feasible(next, remaining - 1)) continue;
      chosen.push_back(i);
      search(next);
      chosen.pop_back();
    }
  };
  search(Mask{});

  Atlas atlas;
  atlas.solutions = std::move(solutions);
  set<string> flow_words;
  set<vector<int>> step_spectra;

  for (const Signature& v : atlas.solutions) {
    if (adjacent_find(v.begin(), v.end()) != v.end()) {
      ++atlas.repeated_support;
    }

    vector<array<int, 7>> defect_rows;
    vector<int> duplicate_targets;
    string flow_word;
    int total_holes = 0;
    int double_rows = 0;
    bool all_flat = true;

    for (int h = 0; h < 13; ++h) {
      array<int, 7> m{};
      for (int x = h; x < L; x += 13) {
        if (u.guard[x]) ++m[x % 7];
      }
      for (int i : v) {
        for (int x = h; x < L; x += 13) {
          if (u.ap[i][x]) ++m[x % 7];
        }
      }

      vector<int> holes;
      vector<int> duplicates;
      array<int, 7> defect{};
      for (int r = 0; r < 7; ++r) {
        atlas.max_base_multiplicity =
            max<uint64_t>(atlas.max_base_multiplicity, m[r]);
        defect[r] = 1 - m[r];
        if (m[r] == 0) holes.push_back(r);
        for (int copy = 1; copy < m[r]; ++copy) duplicates.push_back(r);
      }

      if (holes.size() != duplicates.size()) ++atlas.bad;
      for (int r : holes) {
        if (r != s0 && r != s1) ++atlas.bad;
      }
      if (!holes.empty()) all_flat = false;
      if (holes.size() == 2) ++double_rows;
      require(holes.size() <= 2, "more than two holes in one row");
      total_holes += static_cast<int>(holes.size());
      duplicate_targets.insert(
          duplicate_targets.end(), duplicates.begin(), duplicates.end());
      defect_rows.push_back(defect);

      flow_word.push_back('[');
      for (int r : holes) flow_word.push_back(char('0' + r));
      flow_word.push_back('>');
      for (int r : duplicates) flow_word.push_back(char('0' + r));
      flow_word.push_back(']');
    }

    atlas.total_holes_hist[total_holes]++;
    atlas.any_double_row += (double_rows != 0);
    atlas.max_double_rows =
        max<uint64_t>(atlas.max_double_rows, double_rows);
    if (all_flat) ++atlas.flat;
    flow_words.insert(flow_word);

    {
      set<int> targets(duplicate_targets.begin(), duplicate_targets.end());
      if (targets.size() < duplicate_targets.size()) {
        ++atlas.repeated_duplicate_target;
      }
    }

    bool h_independent = true;
    for (int h = 1; h < 13; ++h) {
      if (defect_rows[h] != defect_rows[0]) h_independent = false;
    }
    if (h_independent && !all_flat) ++atlas.h_independent_nonzero;

    vector<int> steps;
    for (int i : v) steps.push_back(u.step[i]);
    sort(steps.begin(), steps.end());
    step_spectra.insert(steps);
    bool repeated_step = false;
    for (int j = 1; j < 5; ++j) {
      if (steps[j] == steps[j - 1]) repeated_step = true;
    }
    if (!repeated_step) ++atlas.no_repeated_step;

    if (s0 == s1) {
      vector<pair<int, int>> arrows;
      for (int h = 0; h < 13; ++h) {
        int hole = -1;
        int duplicate = -1;
        for (int r = 0; r < 7; ++r) {
          if (defect_rows[h][r] == 1) hole = r;
          if (defect_rows[h][r] == -1) duplicate = r;
        }
        if (hole >= 0) {
          require(hole == s0 && duplicate >= 0,
                  "one-source defect is not an arrow");
          arrows.push_back({h, duplicate});
        }
      }
      atlas.arrow_hist[static_cast<int>(arrows.size())]++;

      bool affine = false;
      for (int a = 0; a < 7; ++a) {
        for (int b = 0; b < 7; ++b) {
          bool works = true;
          for (auto [h, target] : arrows) {
            if ((a * h + b) % 7 != target) works = false;
          }
          if (works) affine = true;
        }
      }
      if (affine) ++atlas.affine_arrow_words;
      if (arrows.size() >= 4) {
        bool positive = false;
        bool negative = false;
        for (auto [h, target] : arrows) {
          (void)h;
          const int delta = (target - s0 + 7) % 7;
          require(delta != 0, "arrow has zero target difference");
          if (delta == 1 || delta == 2 || delta == 4) {
            positive = true;
          } else {
            negative = true;
          }
        }
        if (!positive || !negative) ++atlas.chi4_failures;
      }

    }

    for (int i = 0; i < 5; ++i) {
      for (int j = i + 1; j < 5; ++j) {
        const int p = v[i];
        const int q = v[j];
        if (u.step[p] != u.step[q]) continue;
        const int d = u.step[p];
        const int delta = (u.center[p] - u.center[q] + L) % L;
        atlas.center_bank[d].insert(delta);
        atlas.center_bank[d].insert((L - delta) % L);
      }
    }

    const unsigned holes = hole_residue_mask(u, v);
    atlas.exact_hole_residue_hist[holes]++;
  }

  atlas.flow_words = flow_words.size();
  atlas.step_spectra = step_spectra.size();
  require(atlas.bad == 0, "invalid hole/duplicate flow");
  return atlas;
}

SignatureSet independent_flat_exact_cover(const Universe& u) {
  const Mask required = ~u.guard;
  vector<int> eligible;
  vector<vector<int>> through(L);
  for (int i = 0; i < static_cast<int>(u.ap.size()); ++i) {
    if ((u.ap[i] & u.guard).any()) continue;
    eligible.push_back(i);
    for (int x = 0; x < L; ++x) {
      if (u.ap[i][x]) through[x].push_back(i);
    }
  }
  require(eligible.size() == 182, "guard-avoiding AP count is not 182");

  Signature chosen;
  SignatureSet solutions;
  function<void(const Mask&)> search = [&](const Mask& covered) {
    if (chosen.size() == 5) {
      if (covered == required) {
        Signature v = chosen;
        sort(v.begin(), v.end());
        solutions.insert(v);
      }
      return;
    }
    const Mask missing = required & ~covered;
    int pivot = -1;
    int best = 1'000'000'000;
    for (int x = 0; x < L; ++x) {
      if (!missing[x]) continue;
      int count = 0;
      for (int i : through[x]) {
        if ((u.ap[i] & covered).none()) ++count;
      }
      if (count < best) {
        best = count;
        pivot = x;
      }
    }
    if (pivot < 0) return;
    for (int i : through[pivot]) {
      if ((u.ap[i] & covered).any()) continue;
      chosen.push_back(i);
      search(covered | u.ap[i]);
      chosen.pop_back();
    }
  };
  search(Mask{});
  return solutions;
}

Signature reflected_signature(const Universe& u, const Signature& v) {
  Signature answer;
  for (int i : v) answer.push_back(u.reflected_id(i));
  sort(answer.begin(), answer.end());
  return answer;
}

int reflected_source(int r) {
  return (4 - r + 7) % 7;
}

int inverse_mod91(int a) {
  for (int x = 1; x < L; ++x) {
    if ((a * x) % L == 1) return x;
  }
  throw runtime_error("nonunit inverse request");
}

string double_row_witness(const Universe& u, const Atlas& atlas) {
  for (const Signature& v : atlas.solutions) {
    vector<string> row_words;
    for (int h = 0; h < 13; ++h) {
      array<int, 7> m{};
      for (int x = h; x < L; x += 13) {
        if (u.guard[x]) ++m[x % 7];
      }
      for (int i : v) {
        for (int x = h; x < L; x += 13) {
          if (u.ap[i][x]) ++m[x % 7];
        }
      }
      vector<int> holes;
      vector<int> duplicates;
      for (int r = 0; r < 7; ++r) {
        if (m[r] == 0) holes.push_back(r);
        for (int copy = 1; copy < m[r]; ++copy) duplicates.push_back(r);
      }
      if (holes.size() == 2) {
        string word = to_string(h) + ":[";
        for (int r : holes) word += to_string(r);
        word += ">";
        for (int r : duplicates) word += to_string(r);
        word += "]";
        row_words.push_back(word);
      }
    }
    if (row_words.size() != 2) continue;
    string answer = "supports=";
    for (int i : v) {
      answer += to_string(u.center[i]) + ":" + to_string(u.step[i]) + ",";
    }
    answer.pop_back();
    answer += " rows=" + row_words[0] + "," + row_words[1];
    return answer;
  }
  return "";
}

vector<int> unsigned_step_spectrum(const Universe& u, const Signature& v) {
  vector<int> steps;
  for (int i : v) steps.push_back(u.step[i]);
  sort(steps.begin(), steps.end());
  return steps;
}

bool only_repeated_step_is_one(const vector<int>& steps) {
  bool one_repeats = false;
  for (int i = 0; i < 5;) {
    int j = i + 1;
    while (j < 5 && steps[j] == steps[i]) ++j;
    if (j - i >= 2) {
      if (steps[i] != 1) return false;
      one_repeats = true;
    }
    i = j;
  }
  return one_repeats;
}

tuple<size_t, size_t, map<vector<int>, tuple<uint64_t, size_t, size_t>>>
step_one_exceptional_bank(
    const Universe& u,
    const map<pair<int, int>, Atlas>& atlases,
    const vector<pair<int, int>>& source_configurations) {
  uint64_t assignments = 0;
  SignatureSet unique;
  map<vector<int>, uint64_t> spectrum_assignments;
  map<vector<int>, set<int>> spectrum_centers;
  for (const auto& sources : source_configurations) {
    for (const Signature& v : atlases.at(sources).solutions) {
      const vector<int> steps = unsigned_step_spectrum(u, v);
      if (!only_repeated_step_is_one(steps)) continue;
      ++assignments;
      unique.insert(v);
      spectrum_assignments[steps]++;
      for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
          if (u.step[v[i]] != 1 || u.step[v[j]] != 1) continue;
          const int delta = (u.center[v[i]] - u.center[v[j]] + L) % L;
          spectrum_centers[steps].insert(delta);
          spectrum_centers[steps].insert((L - delta) % L);
        }
      }
    }
  }
  map<vector<int>, tuple<uint64_t, size_t, size_t>> result;
  for (const auto& [steps, centers] : spectrum_centers) {
    set<int> sharp = centers;
    for (int x : centers) sharp.insert((x + L - 1) % L);
    result[steps] = {
        spectrum_assignments.at(steps), centers.size(), sharp.size()};
  }
  return {assignments, unique.size(), result};
}

}  // namespace

int main() {
  try {
    const Universe u;
    cout << "THM-2436 exact companion\n";
    cout << "unoriented_unit_AP13=" << u.ap.size() << "\n";

    for (int h = 0; h < 13; ++h) {
      set<int> columns;
      for (int x = h; x < L; x += 13) {
        if (u.guard[x]) columns.insert(x % 7);
      }
      require(columns == set<int>{h % 7, (h + 6) % 7},
              "guard row is not {h,h-1}");
    }
    cout << "guard_row_law=G_h_{h,h-1}:PASS\n";

    const SignatureSet flat = independent_flat_exact_cover(u);
    require(flat.size() == 62, "independent flat count is not 62");
    cout << "independent_flat_exact_tilings=" << flat.size() << "\n";

    map<pair<int, int>, Atlas> atlases;
    map<unsigned, SignatureSet> exact_hole_bank;
    SignatureSet unique_support_configurations;

    for (int s0 = 0; s0 < 7; ++s0) {
      for (int s1 = s0; s1 < 7; ++s1) {
        Atlas atlas = enumerate_atlas(u, s0, s1);
        require(atlas.flat == 62, "flat subatlas is not 62");
        require(atlas.repeated_support == 0, "repeated AP support survived");
        require(atlas.no_repeated_step == 0,
                "a punctured cover has no repeated unsigned step");
        require(atlas.h_independent_nonzero == 0,
                "nonzero h-independent defect survived");
        require(atlas.max_base_multiplicity <= 2,
                "base multiplicity three survived");
        require(atlas.max_double_rows <= 2,
                "more than two double-hole rows survived");
        require(atlas.chi4_failures == 0,
                "four-arrow word misses one chi_7 sign");

        SignatureSet flat_here;
        for (const Signature& v : atlas.solutions) {
          const unsigned holes = hole_residue_mask(u, v);
          exact_hole_bank[holes].insert(v);
          unique_support_configurations.insert(v);
          if (holes == 0) flat_here.insert(v);
        }
        require(flat_here == flat, "flat subatlas disagrees with Algorithm X");
        atlases[{s0, s1}] = std::move(atlas);
      }
    }

    for (const auto& [sources, atlas] : atlases) {
      const unsigned allowed = source_mask(sources.first, sources.second);
      SignatureSet reconstructed;
      for (const auto& [holes, bank] : exact_hole_bank) {
        if ((holes & ~allowed) == 0) {
          reconstructed.insert(bank.begin(), bank.end());
        }
      }
      require(reconstructed == atlas.solutions,
              "source-superset reconstruction disagrees");
    }

    for (const auto& [sources, atlas] : atlases) {
      const int t0 = reflected_source(sources.first);
      const int t1 = reflected_source(sources.second);
      const pair<int, int> target = minmax(t0, t1);
      SignatureSet reflected;
      for (const Signature& v : atlas.solutions) {
        reflected.insert(reflected_signature(u, v));
      }
      require(reflected == atlases.at(target).solutions,
              "reflection does not transport an atlas exactly");
    }
    cout << "reflection_and_source_superset_controls=PASS\n";

    const array<uint64_t, 7> expected_one{
        264, 257, 229, 257, 264, 679, 679};
    const array<uint64_t, 7> expected_words{
        81, 74, 64, 74, 81, 232, 232};
    const array<uint64_t, 7> expected_affine{
        232, 232, 218, 232, 232, 481, 481};
    cout << "ONE_SOURCE_ATLAS\n";
    cout << "r covers nonflat flow_words arrow_hist affine_words "
            "repeated_target step_spectra\n";
    uint64_t one_total = 0;
    map<int, set<int>> global_center_bank;
    for (int r = 0; r < 7; ++r) {
      const Atlas& a = atlases.at({r, r});
      require(a.solutions.size() == expected_one[r], "one-source count drift");
      require(a.flow_words == expected_words[r], "one-source word count drift");
      require(a.affine_arrow_words == expected_affine[r],
              "one-source affine count drift");
      require(a.repeated_duplicate_target == 0,
              "one-source duplicate target repeats across rows");
      one_total += a.solutions.size();
      cout << r << " " << a.solutions.size() << " "
           << a.solutions.size() - 62 << " " << a.flow_words << " ";
      bool first = true;
      for (auto [arrows, count] : a.arrow_hist) {
        if (!first) cout << ",";
        cout << arrows << ":" << count;
        first = false;
      }
      cout << " " << a.affine_arrow_words << " "
           << a.repeated_duplicate_target << " " << a.step_spectra << "\n";
      for (const auto& [d, bank] : a.center_bank) {
        global_center_bank[d].insert(bank.begin(), bank.end());
      }
    }
    require(one_total == 2629, "one-source total drift");

    uint64_t one_fixed = 0;
    for (const Signature& v : atlases.at({2, 2}).solutions) {
      if (reflected_signature(u, v) == v) ++one_fixed;
    }
    require(one_fixed == 15, "one-source reflection fixed count drift");
    require((one_total + one_fixed) / 2 == 1322,
            "one-source orbit count drift");
    cout << "one_source_total=2629 reflection_fixed=15 "
            "reflection_orbits=1322\n";

    cout << "ONE_SOURCE_REPEATED_PAIR_BANK\n";
    cout << "step center_differences sharp_cells coarse_cells\n";
    const map<int, tuple<int, int, int>> expected_banks{
        {1, {72, 75, 76}}, {2, {52, 57, 62}}, {3, {14, 20, 26}},
        {4, {12, 15, 18}}, {5, {18, 22, 26}}, {44, {12, 14, 16}},
        {45, {16, 20, 24}}};
    for (const auto& [d, bank] : global_center_bank) {
      const int inverse = inverse_mod91(d);
      set<int> q;
      for (int x : bank) q.insert((inverse * x) % L);
      set<int> sharp = q;
      set<int> coarse = q;
      for (int x : q) {
        sharp.insert((x + L - 1) % L);
        coarse.insert((x + L - 1) % L);
        coarse.insert((x + 1) % L);
      }
      require(expected_banks.at(d) ==
                  tuple<int, int, int>{
                      static_cast<int>(bank.size()),
                      static_cast<int>(sharp.size()),
                      static_cast<int>(coarse.size())},
              "repeated-pair bank drift");
      cout << d << " " << bank.size() << " " << sharp.size() << " "
           << coarse.size() << "\n";
    }

    const map<pair<int, int>, uint64_t> expected_two{
        {{0, 1}, 1019}, {{0, 2}, 602},  {{0, 3}, 717},
        {{0, 4}, 912},  {{0, 5}, 2600}, {{0, 6}, 2168},
        {{1, 2}, 944},  {{1, 3}, 583},  {{1, 4}, 717},
        {{1, 5}, 1790}, {{1, 6}, 2096}, {{2, 3}, 944},
        {{2, 4}, 602},  {{2, 5}, 1677}, {{2, 6}, 1677},
        {{3, 4}, 1019}, {{3, 5}, 2096}, {{3, 6}, 1790},
        {{4, 5}, 2168}, {{4, 6}, 2600}, {{5, 6}, 10029}};
    const map<pair<int, int>, uint64_t> expected_two_words{
        {{0, 1}, 448},  {{0, 2}, 258},  {{0, 3}, 315},
        {{0, 4}, 460},  {{0, 5}, 1500}, {{0, 6}, 1021},
        {{1, 2}, 404},  {{1, 3}, 234},  {{1, 4}, 315},
        {{1, 5}, 882},  {{1, 6}, 1042}, {{2, 3}, 404},
        {{2, 4}, 258},  {{2, 5}, 798},  {{2, 6}, 798},
        {{3, 4}, 448},  {{3, 5}, 1042}, {{3, 6}, 882},
        {{4, 5}, 1021}, {{4, 6}, 1500}, {{5, 6}, 5098}};

    cout << "TWO_DISTINCT_SOURCE_ATLAS\n";
    cout << "sources covers exact_both flow_words double_row_covers "
            "max_double_rows max_total_holes step_spectra\n";
    uint64_t two_total = 0;
    uint64_t global_max_double_rows = 0;
    map<int, set<int>> global_two_center_bank;
    for (int s0 = 0; s0 < 7; ++s0) {
      for (int s1 = s0 + 1; s1 < 7; ++s1) {
        const Atlas& a = atlases.at({s0, s1});
        require(a.solutions.size() == expected_two.at({s0, s1}),
                "two-source count drift");
        require(a.flow_words == expected_two_words.at({s0, s1}),
                "two-source word count drift");
        const unsigned both = (1u << s0) | (1u << s1);
        const uint64_t exact_both = exact_hole_bank[both].size();
        const uint64_t inclusion =
            a.solutions.size() - atlases.at({s0, s0}).solutions.size() -
            atlases.at({s1, s1}).solutions.size() + flat.size();
        require(exact_both == inclusion, "two-source inclusion drift");
        const int max_holes = a.total_holes_hist.rbegin()->first;
        cout << s0 << "," << s1 << " " << a.solutions.size() << " "
             << exact_both << " " << a.flow_words << " "
             << a.any_double_row << " " << a.max_double_rows << " "
             << max_holes << " "
             << a.step_spectra << "\n";
        global_max_double_rows =
            max(global_max_double_rows, a.max_double_rows);
        if (a.max_double_rows == 2) {
          cout << "double2_witness " << s0 << "," << s1 << " "
               << double_row_witness(u, a) << "\n";
        }
        two_total += a.solutions.size();
        for (const auto& [d, bank] : a.center_bank) {
          global_two_center_bank[d].insert(bank.begin(), bank.end());
        }
      }
    }
    require(two_total == 38750, "two-distinct-source total drift");
    for (int r = 0; r < 7; ++r) {
      for (const auto& [d, bank] : atlases.at({r, r}).center_bank) {
        global_two_center_bank[d].insert(bank.begin(), bank.end());
      }
    }

    uint64_t two_fixed = 0;
    for (pair<int, int> sources : {
             pair<int, int>{0, 4}, {1, 3}, {5, 6}}) {
      for (const Signature& v : atlases.at(sources).solutions) {
        if (reflected_signature(u, v) == v) ++two_fixed;
      }
    }
    require(two_fixed == 100, "two-source reflection fixed drift");
    require((two_total + two_fixed) / 2 == 19425,
            "two-source orbit count drift");
    cout << "two_distinct_total=38750 reflection_fixed=100 "
            "reflection_orbits=19425\n";
    cout << "TWO_SOURCE_REPEATED_PAIR_BANK\n";
    cout << "step center_differences sharp_cells coarse_cells\n";
    const map<int, tuple<int, int, int>> expected_two_banks{
        {1, {80, 81, 82}},  {2, {66, 71, 74}},  {3, {56, 62, 66}},
        {4, {22, 25, 28}},  {5, {26, 30, 34}},  {6, {6, 7, 8}},
        {30, {8, 10, 12}},  {31, {4, 5, 6}},    {36, {4, 6, 8}},
        {44, {20, 22, 24}}, {45, {38, 44, 46}}};
    for (const auto& [d, bank] : global_two_center_bank) {
      const int inverse = inverse_mod91(d);
      set<int> q;
      for (int x : bank) q.insert((inverse * x) % L);
      set<int> sharp = q;
      set<int> coarse = q;
      for (int x : q) {
        sharp.insert((x + L - 1) % L);
        coarse.insert((x + L - 1) % L);
        coarse.insert((x + 1) % L);
      }
      require(expected_two_banks.at(d) ==
                  tuple<int, int, int>{
                      static_cast<int>(bank.size()),
                      static_cast<int>(sharp.size()),
                      static_cast<int>(coarse.size())},
              "two-source repeated-pair bank drift");
      cout << d << " " << bank.size() << " " << sharp.size() << " "
           << coarse.size() << "\n";
    }

    vector<pair<int, int>> one_source_configurations;
    vector<pair<int, int>> all_two_label_configurations;
    for (int s0 = 0; s0 < 7; ++s0) {
      one_source_configurations.push_back({s0, s0});
      for (int s1 = s0; s1 < 7; ++s1) {
        all_two_label_configurations.push_back({s0, s1});
      }
    }
    const auto [one_exceptional_assignments, one_exceptional_unique,
                one_exceptional] =
        step_one_exceptional_bank(u, atlases, one_source_configurations);
    const auto [two_exceptional_assignments, two_exceptional_unique,
                two_exceptional] =
        step_one_exceptional_bank(u, atlases, all_two_label_configurations);
    require(one_exceptional_assignments == 512 &&
                one_exceptional_unique == 464,
            "one-source exceptional total drift");
    require(two_exceptional_assignments == 8864 &&
                two_exceptional_unique == 5912,
            "two-label exceptional total drift");
    const map<vector<int>, tuple<uint64_t, size_t, size_t>>
        expected_one_exceptional{
            {{1, 1, 1, 1, 1}, {144, 40, 46}},
            {{1, 1, 1, 1, 30}, {106, 36, 41}},
            {{1, 1, 1, 1, 45}, {262, 48, 57}}};
    const map<vector<int>, tuple<uint64_t, size_t, size_t>>
        expected_two_exceptional{
            {{1, 1, 1, 1, 1}, {2681, 60, 66}},
            {{1, 1, 1, 1, 15}, {85, 30, 35}},
            {{1, 1, 1, 1, 18}, {13, 26, 38}},
            {{1, 1, 1, 1, 30}, {1611, 54, 59}},
            {{1, 1, 1, 1, 45}, {4314, 62, 67}},
            {{1, 1, 1, 18, 45}, {16, 10, 16}},
            {{1, 1, 23, 30, 45}, {144, 6, 8}}};
    require(one_exceptional == expected_one_exceptional,
            "one-source exceptional spectrum drift");
    require(two_exceptional == expected_two_exceptional,
            "two-label exceptional spectrum drift");

    cout << "STEP_ONE_ONLY_EXCEPTIONAL_SPECTRA_ONE_SOURCE\n";
    cout << "spectrum assignments center_differences sharp_cells\n";
    size_t one_exceptional_max = 0;
    for (const auto& [steps, data] : one_exceptional) {
      for (int x : steps) cout << x << ",";
      const auto [assignments, centers, sharp] = data;
      cout << " " << assignments << " " << centers << " " << sharp << "\n";
      one_exceptional_max = max(one_exceptional_max, sharp);
    }
    cout << "one_source_exceptional_assignments="
         << one_exceptional_assignments << " unique="
         << one_exceptional_unique << " max_sharp_cells="
         << one_exceptional_max << "\n";

    cout << "STEP_ONE_ONLY_EXCEPTIONAL_SPECTRA_TWO_LABEL\n";
    cout << "spectrum assignments center_differences sharp_cells\n";
    size_t two_exceptional_max = 0;
    for (const auto& [steps, data] : two_exceptional) {
      for (int x : steps) cout << x << ",";
      const auto [assignments, centers, sharp] = data;
      cout << " " << assignments << " " << centers << " " << sharp << "\n";
      two_exceptional_max = max(two_exceptional_max, sharp);
    }
    cout << "two_label_exceptional_assignments="
         << two_exceptional_assignments << " unique="
         << two_exceptional_unique << " max_sharp_cells="
         << two_exceptional_max << "\n";

    const uint64_t assignment_total = one_total + two_total;
    const uint64_t fixed_total = one_fixed + two_fixed;
    require(assignment_total == 41379, "all-source assignment total drift");
    require(unique_support_configurations.size() == 26535,
            "unique support configuration count drift");
    require((assignment_total + fixed_total) / 2 == 20747,
            "all-source orbit count drift");
    cout << "all_source_multiset_assignments=41379\n";
    cout << "unique_ordinary_support_configurations=26535\n";
    cout << "all_reflection_fixed=115 all_reflection_orbits=20747\n";
    cout << "h_independent_nonzero=0\n";
    cout << "one_source_duplicate_target_repetitions=0\n";
    cout << "four_arrow_words_have_both_chi7_signs=PASS\n";
    cout << "max_base_multiplicity=2 max_double_hole_rows="
         << global_max_double_rows << "\n";
    cout << "every_cover_repeats_unsigned_step=PASS\n";
    cout << "ALL CHECKS PASSED\n";
  } catch (const exception& error) {
    cerr << "FAIL: " << error.what() << "\n";
    return 1;
  }
}
