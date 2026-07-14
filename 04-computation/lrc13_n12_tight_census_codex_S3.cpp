// Exact bounded census for primitive tight 12-speed LRC(13) instances.
// codex-2026-07-14-S3, companion audit for THM-765.
//
// Build/run:
//   c++ -O3 -std=c++17 04-computation/lrc13_n12_tight_census_codex_S3.cpp \
//     -o /tmp/lrc13_n12_tight_census
//   /tmp/lrc13_n12_tight_census 30
//
// Exactness.  If S subset {1,...,N} has M(S)>1/13, a maximizer of the lower
// envelope min_v ||vt|| occurs at a pair crossing (denominator v+w or |v-w|)
// or a self-cusp (denominator 2v), hence at p/q with q<=2N.  Conversely, any
// p/q at which all clearances are >1/13 is already a strict witness.  Thus it
// is exact to enumerate all 1<=p<q<=2N.  Keeping maximal safe-speed masks is a
// lossless acceleration.
//
// The explicit 2v/self-cusp clause corrects the denominator omission in the
// local S107/S108 scripts.  (The consolidated lrc14_certificates.M_exact
// already has the correct clause.)
//
// Tournament-analysis note.  Runner-level mask dominance is transitive and
// contributes no structural information, so this census does not promote it
// as a tournament.  THM-765 uses the meaningful vertices: gcd-deck lifts;
// signed killer-phase displacement is the pairwise observable, circle
// reversal is the gauge, and deck order is the antipodal tie path.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

using std::uint64_t;
using std::vector;

int main(int argc, char** argv) {
  const int N = argc > 1 ? std::atoi(argv[1]) : 25;
  if (N < 12 || N > 62) {
    std::cerr << "usage: lrc13_n12_tight_census N, with 12 <= N <= 62\n";
    return 2;
  }

  vector<uint64_t> masks;
  for (int q = 2; q <= 2 * N; ++q) {
    for (int p = 1; p < q; ++p) {
      uint64_t safe = 0;
      for (int v = 1; v <= N; ++v) {
        const int r = static_cast<int>((1LL * v * p) % q);
        const int distance_numerator = std::min(r, q - r);
        if (13 * distance_numerator > q) safe |= uint64_t{1} << (v - 1);
      }
      if (__builtin_popcountll(safe) >= 12) masks.push_back(safe);
    }
  }

  std::sort(masks.begin(), masks.end());
  masks.erase(std::unique(masks.begin(), masks.end()), masks.end());
  const std::size_t raw_mask_count = masks.size();

  std::sort(masks.begin(), masks.end(), [](uint64_t a, uint64_t b) {
    const int pa = __builtin_popcountll(a), pb = __builtin_popcountll(b);
    return pa != pb ? pa > pb : a < b;
  });
  vector<uint64_t> maximal_masks;
  for (uint64_t safe : masks) {
    bool contained = false;
    for (uint64_t larger : maximal_masks) {
      if ((safe & ~larger) == 0) {
        contained = true;
        break;
      }
    }
    if (!contained) maximal_masks.push_back(safe);
  }

  uint64_t total = 0, primitive = 0, hereditary_primitive = 0, tight = 0;
  vector<uint64_t> tight_sets;
  uint64_t subset = (uint64_t{1} << 12) - 1;
  const uint64_t limit = uint64_t{1} << N;

  while (subset < limit) {
    ++total;
    vector<int> speeds;
    speeds.reserve(12);
    int all_gcd = 0;
    for (int bit = 0; bit < N; ++bit) {
      if ((subset >> bit) & 1) {
        speeds.push_back(bit + 1);
        all_gcd = std::gcd(all_gcd, bit + 1);
      }
    }

    if (all_gcd == 1) {
      ++primitive;
      bool hereditary = true;
      for (int skip = 0; skip < 12 && hereditary; ++skip) {
        int core_gcd = 0;
        for (int j = 0; j < 12; ++j) {
          if (j != skip) core_gcd = std::gcd(core_gcd, speeds[j]);
        }
        hereditary = core_gcd == 1;
      }
      if (hereditary) ++hereditary_primitive;

      bool has_strict_witness = false;
      for (uint64_t safe : maximal_masks) {
        if ((subset & ~safe) == 0) {
          has_strict_witness = true;
          break;
        }
      }
      if (!has_strict_witness) {
        ++tight;
        tight_sets.push_back(subset);
      }
    }

    // Gosper's hack: next 12-subset of an N-set.
    const uint64_t low = subset & -subset;
    const uint64_t ripple = subset + low;
    if (ripple == 0 || ripple >= limit) break;
    subset = (((ripple ^ subset) >> 2) / low) | ripple;
  }

  std::cout << "N=" << N << "\n";
  std::cout << "raw_safe_masks=" << raw_mask_count
            << " maximal_safe_masks=" << maximal_masks.size() << "\n";
  std::cout << "total=" << total << " primitive=" << primitive
            << " hereditary_primitive=" << hereditary_primitive
            << " tight=" << tight << "\n";
  for (uint64_t set : tight_sets) {
    std::cout << "{";
    bool first = true;
    for (int bit = 0; bit < N; ++bit) {
      if ((set >> bit) & 1) {
        if (!first) std::cout << ",";
        std::cout << bit + 1;
        first = false;
      }
    }
    std::cout << "}\n";
  }
}
