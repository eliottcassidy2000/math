// Exact truncated-Radon census over the audited THM-2436 punctured atlas.
//
// This translation unit deliberately reuses the pre-existing exhaustive atlas
// engine, but not any THM-2507 implementation.  Renaming the old main exposes
// its Universe/enumerate_atlas machinery in this translation unit; the code
// below independently derives the 13 x 7 defect matrices and their truncated
// Radon zero slopes.

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type"
#endif
#define main thm2436_full_check_main
#include "lrc14_punctured_91_stalk_mixed_mode_thm2436.cpp"
#undef main
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace {

using Defect = array<array<int, 7>, 13>;

struct Census {
  uint64_t assignments = 0;
  uint64_t flat_assignments = 0;
  uint64_t nonflat_assignments = 0;
  map<unsigned, uint64_t> assignment_zero_masks;
  map<Defect, unsigned> patterns;
};

Defect defect_of(const Universe& u, const Signature& signature) {
  Defect defect{};
  for (int h = 0; h < 13; ++h) {
    array<int, 7> multiplicity{};
    for (int x = h; x < L; x += 13) {
      if (u.guard[x]) ++multiplicity[x % 7];
    }
    for (int i : signature) {
      for (int x = h; x < L; x += 13) {
        if (u.ap[i][x]) ++multiplicity[x % 7];
      }
    }
    for (int r = 0; r < 7; ++r) defect[h][r] = 1 - multiplicity[r];
  }
  return defect;
}

bool is_zero(const Defect& defect) {
  for (const auto& row : defect) {
    for (int x : row) {
      if (x != 0) return false;
    }
  }
  return true;
}

unsigned zero_slope_mask(const Defect& defect) {
  unsigned mask = 0;
  for (int kappa = 0; kappa < 13; ++kappa) {
    array<int, 13> radon{};
    for (int v = 0; v < 13; ++v) {
      for (int r = 0; r < 7; ++r) {
        const int h = (v - kappa * r % 13 + 13) % 13;
        radon[v] += defect[h][r];
      }
    }
    if (all_of(radon.begin(), radon.end(), [](int x) { return x == 0; })) {
      mask |= 1u << kappa;
    }
  }
  return mask;
}

void insert(Census& census, const Defect& defect) {
  ++census.assignments;
  for (const auto& row : defect) {
    require(accumulate(row.begin(), row.end(), 0) == 0,
            "defect row sum is not zero");
  }
  if (is_zero(defect)) {
    ++census.flat_assignments;
    return;
  }

  ++census.nonflat_assignments;
  const unsigned mask = zero_slope_mask(defect);
  require(mask & 1u, "row sums did not force the kappa=0 zero slope");
  const int bad = __builtin_popcount(mask);
  require(bad <= 6, "degree-six Radon root floor failed");
  ++census.assignment_zero_masks[mask];
  const auto [it, fresh] = census.patterns.emplace(defect, mask);
  require(fresh || it->second == mask,
          "one defect matrix acquired two Radon zero masks");
}

void merge_patterns(Census& target, const Census& source) {
  target.assignments += source.assignments;
  target.flat_assignments += source.flat_assignments;
  target.nonflat_assignments += source.nonflat_assignments;
  for (const auto& [mask, count] : source.assignment_zero_masks) {
    target.assignment_zero_masks[mask] += count;
  }
  for (const auto& [defect, mask] : source.patterns) {
    const auto [it, fresh] = target.patterns.emplace(defect, mask);
    require(fresh || it->second == mask,
            "merged defect matrix acquired two Radon zero masks");
  }
}

string mask_string(unsigned mask) {
  string answer = "{";
  bool first = true;
  for (int kappa = 0; kappa < 13; ++kappa) {
    if (!((mask >> kappa) & 1u)) continue;
    if (!first) answer += ",";
    answer += to_string(kappa);
    first = false;
  }
  answer += "}";
  return answer;
}

map<unsigned, uint64_t> pattern_mask_histogram(const Census& census) {
  map<unsigned, uint64_t> answer;
  for (const auto& [defect, mask] : census.patterns) {
    (void)defect;
    ++answer[mask];
  }
  return answer;
}

void print_histogram(const string& label,
                     const map<unsigned, uint64_t>& histogram) {
  cout << label;
  for (const auto& [mask, count] : histogram) {
    cout << " " << mask_string(mask) << ":" << count;
  }
  cout << "\n";
}

void print_summary(const string& label, const Census& census) {
  int min_survivors = 13;
  int max_survivors = 0;
  for (const auto& [mask, count] : census.assignment_zero_masks) {
    (void)count;
    const int survivors = 13 - __builtin_popcount(mask);
    min_survivors = min(min_survivors, survivors);
    max_survivors = max(max_survivors, survivors);
  }
  cout << label << " assignments=" << census.assignments
       << " flat=" << census.flat_assignments
       << " nonflat=" << census.nonflat_assignments
       << " unique_nonflat_defects=" << census.patterns.size()
       << " survivor_range=" << min_survivors << ".." << max_survivors
       << "\n";
  print_histogram(label + "_assignment_zero_masks",
                  census.assignment_zero_masks);
  print_histogram(label + "_unique_zero_masks",
                  pattern_mask_histogram(census));
}

}  // namespace

int main() {
  try {
    const Universe universe;
    Census coincident;
    Census distinct;
    Census all;

    for (int s0 = 0; s0 < 7; ++s0) {
      for (int s1 = s0; s1 < 7; ++s1) {
        const Atlas atlas = enumerate_atlas(universe, s0, s1);
        Census local;
        for (const Signature& signature : atlas.solutions) {
          insert(local, defect_of(universe, signature));
        }
        if (s0 == s1) {
          merge_patterns(coincident, local);
        } else {
          merge_patterns(distinct, local);
        }
        merge_patterns(all, local);
      }
    }

    require(coincident.assignments == 2629,
            "coincident-source assignment count drift");
    require(distinct.assignments == 38750,
            "distinct-source assignment count drift");
    require(all.assignments == 41379, "all assignment count drift");
    require(all.flat_assignments == 1736, "flat assignment count drift");
    require(all.nonflat_assignments == 39643,
            "nonflat assignment count drift");
    require(all.patterns.size() == 14952,
            "unique nonflat defect count drift");

    const map<unsigned, uint64_t> expected_assignments{
        {1u << 0, 38810},
        {(1u << 0) | (1u << 1), 112},
        {(1u << 0) | (1u << 2), 28},
        {(1u << 0) | (1u << 3), 16},
        {(1u << 0) | (1u << 4), 77},
        {(1u << 0) | (1u << 5), 51},
        {(1u << 0) | (1u << 6), 93},
        {(1u << 0) | (1u << 7), 106},
        {(1u << 0) | (1u << 8), 49},
        {(1u << 0) | (1u << 9), 42},
        {(1u << 0) | (1u << 10), 85},
        {(1u << 0) | (1u << 11), 71},
        {(1u << 0) | (1u << 12), 103}};
    require(all.assignment_zero_masks == expected_assignments,
            "assignment Radon zero-mask histogram drift");

    const map<unsigned, uint64_t> expected_unique{
        {1u << 0, 14711},
        {(1u << 0) | (1u << 1), 18},
        {(1u << 0) | (1u << 2), 14},
        {(1u << 0) | (1u << 3), 9},
        {(1u << 0) | (1u << 4), 21},
        {(1u << 0) | (1u << 5), 15},
        {(1u << 0) | (1u << 6), 36},
        {(1u << 0) | (1u << 7), 17},
        {(1u << 0) | (1u << 8), 12},
        {(1u << 0) | (1u << 9), 17},
        {(1u << 0) | (1u << 10), 27},
        {(1u << 0) | (1u << 11), 19},
        {(1u << 0) | (1u << 12), 36}};
    require(pattern_mask_histogram(all) == expected_unique,
            "unique Radon zero-mask histogram drift");
    for (const auto& [defect, mask] : all.patterns) {
      (void)defect;
      require(__builtin_popcount(mask) == 1 ||
                  __builtin_popcount(mask) == 2,
              "atlas defect has more than one additional zero slope");
    }

    Defect sharp{};
    sharp[0][0] = -1;
    sharp[0][1] = 1;
    sharp[12][1] = 1;
    sharp[12][2] = -1;
    require(all.patterns.count(sharp) == 1,
            "sharp extra-zero witness left the atlas");
    require(all.patterns.at(sharp) == ((1u << 0) | (1u << 1)),
            "sharp witness zero-mask drift");

    cout << "THM-2507 TRUNCATED-RADON COMPLETE THM-2436-ATLAS CENSUS\n";
    cout << "universe=guard_[0,25]_plus_five_unit_AP13_supports_covering"
            "_with_blocker_columns;recursion_allows_repeats_but_none_survive"
            ";source_multisets=all_28_(s0<=s1)_in_F7\n";
    print_summary("coincident_source", coincident);
    print_summary("distinct_sources", distinct);
    print_summary("all_sources", all);
    cout << "sharp_extra_zero_witness="
            "d(0)=-e0+e1,d(12)=e1-e2;zero_slopes={0,1}\n";
    cout << "ALL CHECKS PASSED\n";
  } catch (const exception& error) {
    cerr << "FAIL: " << error.what() << "\n";
    return 1;
  }
}
