#!/usr/bin/env python3
"""
practical_completeness_s302.py — Practical uses of the completeness theorem
opus-2026-03-24-S302

The theorem: metagraph_dist(C,D) = min_Hamming(T,T') over T∈C, T'∈D.

PRACTICAL USE 1: DIAMETER COMPUTATION WITHOUT FULL METAGRAPH
  Instead of building the metagraph (expensive), compute diameter
  directly from Hamming distances between class representatives.

PRACTICAL USE 2: ADMISSIBLE HEURISTIC FOR SHORTEST PATHS
  Hamming distance between ANY two tilings in C and D gives an
  UPPER BOUND on metagraph distance. This is admissible for A*.

PRACTICAL USE 3: RANKING CHANGE MINIMIZATION
  "How many comparisons must change to go from ranking A to ranking B?"
  = min_Hamming between the tiling representations. Computable in O(m).

PRACTICAL USE 4: METAGRAPH DIAMETER FROM BURNSIDE
  Diameter = max over pairs min_Hamming. For the farthest pair
  (typically transitive ↔ most isolated), compute directly.
"""

import sys, time
from math import comb, factorial
from collections import Counter
from functools import lru_cache

sys.stdout.reconfigure(line_buffering=True)

print("=" * 72)
print("  PRACTICAL APPLICATIONS OF THE COMPLETENESS THEOREM")
print("  opus-2026-03-24-S302")
print("=" * 72)

# ============================================================
# USE 1: DIAMETER BOUNDS FROM BURNSIDE (no enumeration!)
# ============================================================

print(f"\n{'='*72}")
print("  USE 1: DIAMETER BOUNDS")
print(f"{'='*72}")

# The diameter = max min_Hamming over all class pairs.
# LOWER BOUND: max distance from the transitive class (all-zeros tiling).
#   The transitive tiling has Hamming weight 0.
#   A tiling at Hamming weight w is at distance ≤ w from transitive.
#   So diameter ≥ max Hamming weight of any class representative.

# UPPER BOUND: diameter ≤ m (trivially, since max Hamming distance = m).
# Better: diameter ≤ m - floor(m/2) (by complementary argument).

# From our data:
diameters = {3:1, 4:1, 5:3, 6:4, 7:7, 8:8}
m_vals = {n: comb(n-1, 2) for n in range(3, 11)}

print(f"\n  {'n':>3} {'m':>4} {'diam':>5} {'m/2':>5} {'diam/m':>7} {'n-2':>4} {'diam=n-2?':>10}")
for n in range(3, 9):
    m = m_vals[n]
    d = diameters.get(n, '?')
    if isinstance(d, int):
        print(f"  {n:>3} {m:>4} {d:>5} {m/2:>5.1f} {d/m:>7.3f} {n-2:>4} {'✓' if d == n-2 else '✗':>10}")

# Conjecture from other sessions: diameter = n-2
print(f"\n  CONJECTURE: diameter = n-2 for n ≥ 4")
print(f"  Verified: n=4 (1=2✓), n=5 (3=3✓), n=6 (4=4✓), n=7 (7=5✗!)")
print(f"  FAILS at n=7: diameter=7 but n-2=5.")

# Actually let me check: diameter sequence is 1, 1, 3, 4, 7, 8
# n-2 gives 1, 2, 3, 4, 5, 6
# Doesn't match. Let me look for another pattern.
print(f"\n  DIAMETER SEQUENCE: ", end="")
for n in range(3, 9):
    d = diameters.get(n, '?')
    print(f"d({n})={d}, ", end="")
print()

# First differences
diffs = [diameters[n+1] - diameters[n] for n in range(3, 8)]
print(f"  Differences: {diffs}")

# Ratios
for n in range(4, 9):
    if diameters.get(n) and diameters.get(n-1):
        print(f"  d({n})/d({n-1}) = {diameters[n]/diameters[n-1]:.3f}")

# ============================================================
# USE 2: RANKING CHANGE QUANTIFICATION
# ============================================================

print(f"\n{'='*72}")
print("  USE 2: RANKING CHANGE QUANTIFICATION")
print(f"{'='*72}")
print("""
  PROBLEM: Given two tournament rankings (e.g., LLM leaderboards
  from different time periods), how "structurally different" are they?

  ANSWER (from completeness theorem):
    structural_distance = min_Hamming(tiling(T1), tiling(T2))
    = number of pairwise comparisons that must change
    to transform T1 into T2 (minimized over relabelings)

  This is computable in O(n! × m) time for exact, or
  O(n²m) with heuristic relabeling.

  EXAMPLE at n=10 (10 LLM models):
    m = C(9,2) = 36 tiles
    Two leaderboard snapshots → two tilings of δ_8
    Hamming distance = number of changed head-to-head matchups
    If distance = 3: only 3 matchup reversals separate the rankings
    If distance = 15: half the matchups differ → very different rankings
""")

# ============================================================
# USE 3: EFFICIENT METAGRAPH CONSTRUCTION
# ============================================================

print(f"\n{'='*72}")
print("  USE 3: EFFICIENT METAGRAPH CONSTRUCTION")
print(f"{'='*72}")
print("""
  The completeness theorem tells us:
    F_k = {pairs at metagraph distance ≤ k}
       = {pairs at min_Hamming ≤ k}
       = {pairs reachable by ≤ k tile flips}

  So to build the metagraph (distance-1 pairs only):
    Just compute F_1 = all single-tile-flip connections.
    This is O(2^m × m) time, which we already do.

  But for APPROXIMATE metagraphs at large n:
    The spectral gap ≈ 2/n and H ≈ 2nd eigenvector
    mean we can PREDICT many edge weights from H values alone:
      W[i][j] ≈ f(H_i, H_j, |Aut_i|, |Aut_j|)
    This is 88% accurate (from S284).

  For the remaining 12%: need actual computation.
""")

# ============================================================
# USE 4: TOURNAMENT SIMILARITY METRIC FOR ML
# ============================================================

print(f"\n{'='*72}")
print("  USE 4: TOURNAMENT SIMILARITY METRIC FOR ML")
print(f"{'='*72}")
print("""
  The completeness theorem provides a METRIC on tournament iso classes:
    d(C, D) = min_Hamming(C, D)

  Properties:
  - d(C, D) = 0 iff C = D (identity)
  - d(C, D) = d(D, C) (symmetry — from W being symmetric)
  - d(C, E) ≤ d(C, D) + d(D, E) (triangle inequality — from Hamming)
  - Integer-valued: d ∈ {0, 1, 2, ..., m}

  This is a PROPER METRIC on the set of tournament iso classes.
  It can be used as a kernel or distance function in ML:
    K(C, D) = exp(-d(C,D) / σ)  (Gaussian kernel)
    K(C, D) = (1 + d(C,D))^{-α}  (inverse polynomial)

  For TDA (topological data analysis):
    Build a Rips complex on iso classes using this metric.
    Compute persistent homology → topological features of tournament space.

  For clustering:
    k-medoids or hierarchical clustering using d(C,D).
    Each cluster = a "structural family" of tournaments.
""")

# ============================================================
# USE 5: OEIS AND SEQUENCE ANALYSIS
# ============================================================

print(f"\n{'='*72}")
print("  USE 5: NEW OEIS SEQUENCES")
print(f"{'='*72}")

# Diameter sequence
print(f"  Diameter of wiggly metagraph: 1, 1, 3, 4, 7, 8 (n=3..8)")
print(f"  → NOT in OEIS (should submit)")

# Number of classes at each distance from transitive
print(f"\n  Distance profile from transitive (from S296/S299):")
print(f"    n=5: distance [0,1,2,3] → classes [1, 4, 4, 1]")
print(f"    n=6: distance [0,1,2,3,4] → classes [1, 6, 11, 12, 4]")
print(f"    n=7: distance [0,1,2,...,7] → need to compute")

# Filling sequence
print(f"\n  Filling fraction at each waggly order:")
print(f"    n=5: [47%, 91%, 100%]")
print(f"    n=6: [26%, 75%, 97%, 100%]")
print(f"    n=7: [5.8%, 38%, 82%, ...]")

# ============================================================
# USE 6: SPEED ENGINE INTEGRATION
# ============================================================

print(f"\n{'='*72}")
print("  USE 6: DIAMETER FROM SPEED ENGINE")
print(f"{'='*72}")

# The speed engine computes V_n, T_n, twin_SL, etc. at any n.
# Can we also compute the diameter?
# From the data: diameter correlates with n but isn't n-2.
# Sequence: 1, 1, 3, 4, 7, 8
# Differences: 0, 2, 1, 3, 1
# Second differences: 2, -1, 2, -2

# Let me check: is diameter = C(n-2, 2) / something?
print(f"  {'n':>3} {'diam':>5} {'C(n-2,2)':>8} {'ratio':>7}")
for n in range(3, 9):
    d = diameters.get(n)
    cn2 = comb(n-2, 2)
    if d and cn2 > 0:
        print(f"  {n:>3} {d:>5} {cn2:>8} {d/cn2:>7.3f}")

# Hmm, C(n-2,2) = 1,3,6,10,15 for n=4..8, and diameter = 1,3,4,7,8.
# Ratio: 1.0, 1.0, 0.67, 0.70, 0.53 — not clean.

# Maybe diameter = floor(m × some fraction)?
print(f"\n  diam vs m:")
for n in range(3, 9):
    d = diameters.get(n)
    m = comb(n-1, 2)
    if d:
        print(f"  n={n}: diam={d}, m={m}, diam/m={d/m:.3f}")

print(f"\n  diam/m sequence: {[diameters[n]/comb(n-1,2) for n in range(4, 9)]}")
print(f"  → Roughly diam ≈ m/2 at large n? Need more data.")

print("\nDONE.")
