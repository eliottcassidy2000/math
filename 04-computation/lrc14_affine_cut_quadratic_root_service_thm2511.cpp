// Exact THM-2511 census: the positive cut-bundle energy retains a
// nonconstant C_13 current on every nonflat THM-2436 atlas defect for at
// least 64 of the 72 nonzero (tau,a) pairs.

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
    for (int value : row) {
      if (value != 0) return false;
    }
  }
  return true;
}

using Energy = array<long long, 13>;

Energy cut_energy(const Defect& defect, int tau, int a) {
  array<long long, 13> energy{};
  for (int v = 0; v < 13; ++v) {
    for (int c = 0; c < 7; ++c) {
      long long value = 0;
      for (int r = 0; r < 7; ++r) {
        const int representative = (a * r + c) % 7;
        const int h = (v - tau * representative % 13 + 13) % 13;
        value += defect[h][r];
      }
      energy[v] += value * value;
    }
  }
  return energy;
}

bool energy_is_constant(const Energy& energy) {
  return all_of(energy.begin() + 1, energy.end(),
                [&](long long value) { return value == energy[0]; });
}

}  // namespace

int main() {
  try {
    const Universe universe;
    map<Defect, unsigned> patterns;
    uint64_t assignments = 0;
    uint64_t flat_assignments = 0;
    for (int s0 = 0; s0 < 7; ++s0) {
      for (int s1 = s0; s1 < 7; ++s1) {
        const Atlas atlas = enumerate_atlas(universe, s0, s1);
        for (const Signature& signature : atlas.solutions) {
          ++assignments;
          const Defect defect = defect_of(universe, signature);
          if (is_zero(defect)) {
            ++flat_assignments;
          } else {
            patterns.emplace(defect, 0);
          }
        }
      }
    }

    uint64_t constant_pairs = 0;
    uint64_t defects_with_constant_pair = 0;
    uint64_t defects_with_all_pairs_constant = 0;
    map<int, uint64_t> nonconstant_pair_histogram;
    int minimum_nonconstant_pairs = 72;
    long long maximum_total_energy = 0;
    Defect minimum_witness{};
    for (const auto& [defect, ignored] : patterns) {
      (void)ignored;
      int defect_l1 = 0;
      for (const auto& row : defect) {
        require(accumulate(row.begin(), row.end(), 0) == 0,
                "defect row sum is not zero");
        for (int value : row) defect_l1 += abs(value);
      }
      require(defect_l1 <= 18, "inherited L1 bound failed");
      int nonconstant_pairs = 0;
      for (int tau = 1; tau < 13; ++tau) {
        for (int a = 1; a < 7; ++a) {
          const Energy energy = cut_energy(defect, tau, a);
          const long long total_energy =
              accumulate(energy.begin(), energy.end(), 0LL);
          require(total_energy > 0, "nonflat defect acquired zero cut energy");
          require(total_energy <= 7LL * 18 * 18,
                  "quadratic energy exceeded 7*18^2");
          maximum_total_energy = max(maximum_total_energy, total_energy);
          if (energy_is_constant(energy)) {
            ++constant_pairs;
          } else {
            ++nonconstant_pairs;
          }
        }
      }
      if (nonconstant_pairs < 72) ++defects_with_constant_pair;
      if (nonconstant_pairs == 0) ++defects_with_all_pairs_constant;
      ++nonconstant_pair_histogram[nonconstant_pairs];
      if (nonconstant_pairs < minimum_nonconstant_pairs) {
        minimum_nonconstant_pairs = nonconstant_pairs;
        minimum_witness = defect;
      }
    }

    const map<int, uint64_t> expected_histogram{
        {64, 4}, {66, 107}, {68, 1498}, {70, 1977}, {72, 11366}};
    require(assignments == 41379, "assignment count drifted");
    require(flat_assignments == 1736, "flat assignment count drifted");
    require(patterns.size() == 14952, "unique nonflat defect count drifted");
    require(constant_pairs == 10620, "constant-energy incidence count drifted");
    require(defects_with_constant_pair == 3586,
            "constant-pair hostile count drifted");
    require(defects_with_all_pairs_constant == 0,
            "all-pairs quadratic escape failed");
    require(minimum_nonconstant_pairs == 64,
            "64-of-72 quadratic invoice drifted");
    require(nonconstant_pair_histogram == expected_histogram,
            "nonconstant-pair histogram drifted");

    cout << "THM-2510 AFFINE-CUT QUADRATIC ROOT-SERVICE ATLAS CENSUS\n";
    cout << "universe=THM-2436_guard_[0,25]_plus_five_unit_AP13_supports"
            "_and_blocker_columns;all_28_source_multisets\n";
    cout << "assignments=" << assignments << " flat=" << flat_assignments
         << " nonflat=" << assignments - flat_assignments << "\n";
    cout << "unique_nonflat_defects=" << patterns.size() << "\n";
    cout << "constant_energy_pairs=" << constant_pairs << "\n";
    cout << "defects_with_constant_pair=" << defects_with_constant_pair << "\n";
    cout << "defects_with_all_pairs_constant="
         << defects_with_all_pairs_constant << "\n";
    cout << "minimum_nonconstant_pairs=" << minimum_nonconstant_pairs << "/72\n";
    cout << "nonconstant_pair_histogram=";
    for (const auto& [count, defects] : nonconstant_pair_histogram) {
      cout << " " << count << ":" << defects;
    }
    cout << "\n";
    cout << "minimum_witness_constant_pairs=";
    for (int tau = 1; tau < 13; ++tau) {
      for (int a = 1; a < 7; ++a) {
        if (energy_is_constant(cut_energy(minimum_witness, tau, a))) {
          cout << " (" << tau << "," << a << ")";
        }
      }
    }
    cout << "\n";
    cout << "minimum_witness_rows=";
    for (int h = 0; h < 13; ++h) {
      bool nonzero = any_of(minimum_witness[h].begin(),
                            minimum_witness[h].end(),
                            [](int x) { return x != 0; });
      if (!nonzero) continue;
      cout << " h" << h << "=[";
      for (int r = 0; r < 7; ++r) {
        if (r) cout << ",";
        cout << minimum_witness[h][r];
      }
      cout << "]";
    }
    cout << "\n";
    cout << "maximum_total_energy=" << maximum_total_energy
         << ";universal_bound=2268\n";
    cout << "fixed-pair_Fubini_floor=64/72=8/9;"
            "parent_floors=16/63,8/21\n";
    cout << "nonconstant rational C13 vector iff all 12 nonzero colours;"
            "algebraic-integer floor=2268^-11\n";
    cout << "ALL EXACT CHECKS PASSED\n";
  } catch (const exception& error) {
    cerr << "FAIL: " << error.what() << "\n";
    return 1;
  }
}
