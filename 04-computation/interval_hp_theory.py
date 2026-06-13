#!/usr/bin/env python3
"""
Why does the cyclic interval tournament maximize H?

The interval tournament C_p has S = {1,...,m} where m = (p-1)/2.
This means: i -> j iff 0 < (j-i) mod p <= m.

Equivalently: vertex i beats the "next" m vertices in the cyclic order.

STRUCTURAL INSIGHT: In C_p, the natural cyclic ordering 0,1,2,...,p-1
is "almost" a Hamiltonian path. Every arc i -> i+1 exists (since 1 is in S).
In fact, EVERY k-step forward arc (i -> i+k for k <= m) exists.

A Hamiltonian path in C_p is a permutation (v_0, v_1, ..., v_{p-1}) where
v_i -> v_{i+1} for all i, meaning 0 < (v_{i+1} - v_i) mod p <= m.

This is equivalent to: the SUCCESSIVE DIFFERENCES d_i = (v_{i+1} - v_i) mod p
must all be in {1,...,m}.

So H(C_p) = #{sequences (d_0,...,d_{p-2}) with d_i in {1,...,m},
             sum(d_i) = 0 mod p, all partial sums distinct mod p}

The constraint "all partial sums distinct" makes this a RESTRICTED composition.

For PALEY: the condition is d_i in QR_p instead of {1,...,m}.
The QR set is "spread out" (elements from 1 to p-1 with no pattern),
while {1,...,m} is a CONTIGUOUS block.

HYPOTHESIS: The contiguous block {1,...,m} allows more valid compositions
because the "small step" options (d_i = 1,2,...) are available, which
create tighter packing of partial sums and avoid collisions better.

This script tests this composition counting interpretation.
"""

import numpy as np
from itertools import product

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p-1)//2, p) == 1

def count_hp_composition(S, p):
    """Count Hamiltonian paths as restricted compositions.

    HP(C_p, S) = #{(d_0,...,d_{p-2}) : d_i in S,
                   partial sums s_k = sum_{i<k} d_i distinct mod p,
                   s_{p-1} = 0 mod p}

    Use DP: state = (last partial sum, set of visited partial sums)
    """
    n = p
    m = n - 1  # number of steps

    # Start at partial sum 0 (vertex v_0 = 0 WLOG)
    # dp[mask][s] = # ways to reach partial sum s with visited set mask
    # mask has bit k set iff partial sum k has been visited
    # Starting: mask = {0}, s = 0

    dp = {}
    dp[(1, 0)] = 1  # mask=1 (only vertex 0 visited), current sum=0

    for step in range(m):
        new_dp = {}
        for (mask, s), count in dp.items():
            for d in S:
                s_new = (s + d) % p
                bit = 1 << s_new
                if mask & bit:
                    continue  # already visited
                new_mask = mask | bit
                key = (new_mask, s_new)
                new_dp[key] = new_dp.get(key, 0) + count
        dp = new_dp

    # Final: all vertices visited (mask = 2^p - 1), end at sum = 0 mod p
    # Wait: the last partial sum s_{p-1} doesn't need to be 0.
    # The partial sums s_0=0, s_1, ..., s_{p-1} must all be distinct mod p.
    # Since there are p partial sums and p residues, they must be a PERMUTATION
    # of Z_p. So s_{p-1} is determined (the missing element).
    # But we also need v_0 arbitrary, not fixed at 0.

    # Actually, for CIRCULANT tournaments, fixing v_0 = 0 gives H/p of the paths
    # (by rotational symmetry). So total H = p * (# paths starting at 0).

    # For paths starting at 0: the partial sums s_k = v_k are the vertices visited.
    # We need them to be a permutation of Z_p.

    full_mask = (1 << p) - 1
    total = sum(count for (mask, s), count in dp.items() if mask == full_mask)

    return total * p  # multiply by p for all starting vertices

# But this DP is essentially the same as the adjacency matrix DP!
# Let me just verify the connection and then do the analysis.

# For small p, compare:
for p in [5, 7]:
    m = (p-1)//2
    S_int = set(range(1, m+1))
    S_pal = set(j for j in range(1,p) if is_qr(j,p))

    H_int = count_hp_composition(S_int, p)
    H_pal = count_hp_composition(S_pal, p)

    print(f"p={p}: H(interval)={H_int}, H(Paley)={H_pal}")

# Now analyze WHY the interval produces more paths.
# The key difference: {1,...,m} vs QR_p.

print(f"\n{'='*60}")
print(f"STEP SIZE DISTRIBUTION ANALYSIS")
print(f"{'='*60}")

for p in [7, 11, 19]:
    m = (p-1)//2
    S_int = set(range(1, m+1))
    S_pal = set(j for j in range(1,p) if is_qr(j,p))

    # Average step size
    avg_int = sum(S_int) / len(S_int)
    avg_pal = sum(S_pal) / len(S_pal)

    print(f"\np={p}:")
    print(f"  Interval: S = {sorted(S_int)}")
    print(f"    avg step = {avg_int:.2f}, max step = {max(S_int)}")
    print(f"  Paley: S = {sorted(S_pal)}")
    print(f"    avg step = {avg_pal:.2f}, max step = {max(S_pal)}")

    # "Small step" count: steps of size 1 or 2
    small_int = sum(1 for s in S_int if s <= 2)
    small_pal = sum(1 for s in S_pal if s <= 2)
    print(f"  Small steps (<=2): interval={small_int}, Paley={small_pal}")

    # "Consecutive" property: how many consecutive integers?
    consecutive_int = max(S_int) - min(S_int) + 1 == len(S_int)
    consecutive_pal = max(S_pal) - min(S_pal) + 1 == len(S_pal)
    print(f"  Consecutive? interval={consecutive_int}, Paley={consecutive_pal}")

    # Gap structure
    sorted_pal = sorted(S_pal)
    gaps = [sorted_pal[i+1] - sorted_pal[i] for i in range(len(sorted_pal)-1)]
    print(f"  Paley gaps: {gaps}")

# Theoretical analysis
print(f"\n{'='*60}")
print(f"THEORETICAL INSIGHT: CONSECUTIVE vs SCATTERED STEPS")
print(f"{'='*60}")
print("""
The cyclic interval S = {1,...,m} is a CONSECUTIVE block of step sizes.
The Paley set QR_p is SCATTERED across {1,...,p-1}.

Key advantage of consecutive steps:
1. FINE-GRAINED CONTROL: steps of size 1 and 2 allow precise navigation
   around the cycle, avoiding collisions with previously visited vertices.

2. GREEDY PACKING: with small steps available, the HP can fill in gaps
   between large jumps. A path that makes a large jump (step m) can then
   use small steps (1,2,...) to fill the nearby vertices.

3. NUMBER OF VALID COMPOSITIONS: for a random walk on Z_p with steps in S,
   the probability that all partial sums are distinct depends on how
   "fine-grained" the step set is. Consecutive steps give maximum
   fine-grainedness.

At small p (7, 11): the QR set still has good coverage (e.g., QR_7 = {1,2,4}
   includes steps 1 and 2). The additional structure of QR (balanced phases,
   Gauss sum properties) compensates for the non-consecutive gaps.

At large p (19+): the QR set has larger gaps and misses critical small steps.
   The QR_19 = {1,4,5,6,7,9,11,16,17} includes 1 but not {2,3}. The gaps
   {2,3,8,10,12,13,14,15,18} are significant, and many fine-grained
   navigation options are lost.

PREDICTION: The crossover happens when the "gap penalty" of QR exceeds
the "phase alignment bonus" of the Gauss sum structure.
""")

# Count the NUMBER of missing consecutive integers from QR
print(f"Missing consecutive pairs (k, k+1 both missing):")
for p in [7, 11, 19, 23, 29, 31]:
    if p % 4 != 3:
        continue
    S_pal = set(j for j in range(1,p) if is_qr(j,p))
    missing = set(range(1,p)) - S_pal
    missing_pairs = sum(1 for k in range(1,p-1)
                       if k in missing and k+1 in missing)
    print(f"  p={p}: QR={sorted(S_pal)}")
    print(f"    Missing from {1,...,{(p-1)//2}}: "
          f"{sorted(k for k in range(1,(p-1)//2+1) if k not in S_pal)}")
    print(f"    Total missing consecutive pairs: {missing_pairs}")
