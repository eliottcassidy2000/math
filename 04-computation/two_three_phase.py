"""
two_three_phase.py -- kind-pasteur-2026-03-13-S62

The b1 PHASE TRANSITION between n=7 and n=9.

At n=7: b1 = a1 - 2*a2 > 0 ALWAYS (400/400 samples)
At n=9: b1 < 0 ALWAYS (just from 3-cycles: 30/30 samples)

Questions:
1. What happens at n=8 (even n)?
2. What's the EXACT transition point?
3. How does this relate to Omega graph density?
4. What does b1 < 0 mean for the 3-adic tower?

The key structural reason:
  At n=7: two 3-cycles use 6 vertices, leaving 1. There's limited room
  for disjoint pairs.
  At n=9: two 3-cycles use 6 vertices, leaving 3. The unused 3 vertices
  can form ANOTHER 3-cycle, giving alpha_3 > 0 and many more disjoint pairs.

The transition is at n=2k where two k-cycles can first be disjoint.
For k=3: n=6 is the first time two 3-cycles can be disjoint.
But at n=6, C(6,3)=20 triples, and the ~10 cyclic ones have ~7 disjoint pairs.
With a1 ~ 10, a2 ~ 7, b1 = 10 - 14 = -4. So b1 < 0 might ALREADY happen at n=6!

Wait: but at n=5, a2=0 always (two triples need 6 > 5).
At n=6: two triples CAN be disjoint (6=3+3). So a2 > 0 possible.
At n=7: same, but with 7 > 6, so "more room" for overlap.

The question is the RATIO a2/a1. At small n, a2/a1 is small (few disjoint pairs).
As n grows, a2/a1 grows because there's more room for disjoint placements.

Let me check this systematically.
"""

import numpy as np
from itertools import combinations
from collections import Counter

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def get_c3_alpha(A, n):
    """Get alpha_1 and alpha_2 counting only 3-cycles (fast)."""
    cycles = []
    for combo in combinations(range(n), 3):
        sub = np.zeros((3, 3), dtype=int)
        verts = list(combo)
        for a in range(3):
            for b in range(3):
                sub[a][b] = A[verts[a]][verts[b]]
        c = int(np.trace(sub @ sub @ sub)) // 3
        if c > 0:
            cycles.append(frozenset(combo))

    a1 = len(cycles)
    a2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                a2 += 1

    return a1, a2

def rand_bits(total_bits):
    """Generate random bits for large total_bits."""
    if total_bits <= 30:
        return np.random.randint(0, 1 << total_bits)
    # Build bits in chunks
    bits = 0
    remaining = total_bits
    shift = 0
    while remaining > 0:
        chunk = min(remaining, 30)
        bits |= int(np.random.randint(0, 1 << chunk)) << shift
        shift += chunk
        remaining -= chunk
    return bits

# ============================================================
# PART 1: b1 sign across n = 3 to 11
# ============================================================
print("=" * 70)
print("PART 1: b1 sign transition across n = 3 to 11")
print("=" * 70)

np.random.seed(42)

for n in range(3, 12):
    total_bits = n*(n-1)//2
    n_samples = min(2**total_bits, 500) if n <= 6 else 200

    b1_neg = 0
    b1_zero = 0
    b1_pos = 0
    a1_sum = 0
    a2_sum = 0

    for trial in range(n_samples):
        if n <= 6 and trial < 2**total_bits:
            bits = trial
        else:
            bits = rand_bits(total_bits)

        A = bits_to_adj(bits, n)
        a1, a2 = get_c3_alpha(A, n)

        b1 = a1 - 2*a2
        a1_sum += a1
        a2_sum += a2

        if b1 < 0:
            b1_neg += 1
        elif b1 == 0:
            b1_zero += 1
        else:
            b1_pos += 1

    total = b1_neg + b1_zero + b1_pos
    avg_a1 = a1_sum / total
    avg_a2 = a2_sum / total
    ratio = avg_a2 / avg_a1 if avg_a1 > 0 else 0

    print(f"  n={n:>2}: b1<0: {b1_neg:>4} ({100*b1_neg/total:>5.1f}%), "
          f"b1=0: {b1_zero:>4} ({100*b1_zero/total:>5.1f}%), "
          f"b1>0: {b1_pos:>4} ({100*b1_pos/total:>5.1f}%), "
          f"avg(a1)={avg_a1:>6.1f}, avg(a2)={avg_a2:>7.1f}, "
          f"a2/a1={ratio:>6.3f}")

# ============================================================
# PART 2: b0 sign across n
# ============================================================
print("\n" + "=" * 70)
print("PART 2: b0 sign transition across n = 3 to 11")
print("=" * 70)

np.random.seed(42)

for n in range(3, 12):
    total_bits = n*(n-1)//2
    n_samples = min(2**total_bits, 500) if n <= 6 else 200

    b0_neg = 0
    b0_zero = 0
    b0_pos = 0

    for trial in range(n_samples):
        if n <= 6 and trial < 2**total_bits:
            bits = trial
        else:
            bits = rand_bits(total_bits)

        A = bits_to_adj(bits, n)
        a1, a2 = get_c3_alpha(A, n)
        b0 = 1 - a1 + a2

        if b0 < 0:
            b0_neg += 1
        elif b0 == 0:
            b0_zero += 1
        else:
            b0_pos += 1

    total = b0_neg + b0_zero + b0_pos
    print(f"  n={n:>2}: b0<0: {b0_neg:>4} ({100*b0_neg/total:>5.1f}%), "
          f"b0=0: {b0_zero:>4} ({100*b0_zero/total:>5.1f}%), "
          f"b0>0: {b0_pos:>4} ({100*b0_pos/total:>5.1f}%)")

# ============================================================
# PART 3: Omega graph density across n
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Omega(T) density (3-cycle conflict graph)")
print("=" * 70)

np.random.seed(42)

for n in range(3, 12):
    total_bits = n*(n-1)//2
    n_samples = min(2**total_bits, 300) if n <= 6 else 100

    densities = []

    for trial in range(n_samples):
        if n <= 6 and trial < 2**total_bits:
            bits = trial
        else:
            bits = rand_bits(total_bits)

        A = bits_to_adj(bits, n)
        a1, a2 = get_c3_alpha(A, n)

        max_edges = a1*(a1-1)//2 if a1 >= 2 else 0
        if max_edges > 0:
            e_omega = max_edges - a2
            density = e_omega / max_edges
            densities.append(density)

    if densities:
        print(f"  n={n:>2}: avg density = {np.mean(densities):.4f}, "
              f"min = {min(densities):.4f}, max = {max(densities):.4f}")
    else:
        print(f"  n={n:>2}: no data (all transitive?)")

# ============================================================
# PART 4: Why the transition? Counting argument
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Why the b1 transition happens")
print("=" * 70)

print("""
The key is the ratio a2/a1 (disjoint pairs / total cycles).

b1 = a1 - 2*a2 > 0  iff  a2/a1 < 1/2

For 3-cycles only:
  a1 = c3 = # directed 3-cycles (vertex sets with a Ham cycle)
  a2 = #{(C1, C2) : C1, C2 are 3-cycles on disjoint vertex sets}

  Expected a1 = C(n,3) * P(cyclic) where P(cyclic) = 3/4 for random T
  (3 of 4 tournament types on 3 vertices are acyclic? No, 2 of 4...
   Actually for 3-vertex tournament, P(has 3-cycle) = 1/4.
   Wait: C(3,2)=3 edges, each with 2 orientations = 8 tournaments.
   Cyclic ones: 2 (clockwise and counterclockwise). So P = 2/8 = 1/4.)

  Expected c3 = C(n,3) / 4

  Expected a2: for random T, two disjoint triples each have P=1/4 of being cyclic.
  # disjoint triple pairs = C(n,3) * C(n-3,3) / 2 (choosing 2 disjoint triples)

  Expected a2 = C(n,3) * C(n-3,3) / 2 * (1/4)^2 = C(n,3)*C(n-3,3)/32

  Ratio: a2/a1 ~ (C(n-3,3)/32) / (1/4) = C(n-3,3)/8

  For n=5: C(2,3)/8 = 0  (can't pick 3 from 2)
  For n=6: C(3,3)/8 = 1/8 = 0.125
  For n=7: C(4,3)/8 = 4/8 = 0.5
  For n=8: C(5,3)/8 = 10/8 = 1.25
  For n=9: C(6,3)/8 = 20/8 = 2.5

  TRANSITION at n=7: expected a2/a1 ~ 0.5, so b1 ~ 0!
  At n=8: a2/a1 ~ 1.25 > 0.5, so b1 < 0!

  This predicts:
  n <= 6: b1 > 0 (almost always)
  n = 7: b1 ~ 0 (transition point)
  n >= 8: b1 < 0 (almost always)
""")

# Verify the prediction
print("  Verification:")
for n in range(3, 12):
    if n >= 6:
        from math import comb
        predicted_ratio = comb(n-3, 3) / 8
        print(f"    n={n}: predicted a2/a1 ~ {predicted_ratio:.3f}, "
              f"b1 {'> 0' if predicted_ratio < 0.5 else '~ 0' if abs(predicted_ratio-0.5) < 0.1 else '< 0'}")

# ============================================================
# PART 5: The full alpha independence polynomial at n=8
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Full independence structure at n=8 (3-cycles only)")
print("=" * 70)

n = 8
total_bits = n*(n-1)//2
np.random.seed(42)

# At n=8: max independent set of 3-cycles
# Two disjoint 3-cycles use 6 vertices, leaving 2. Cannot have 3rd cycle.
# So alpha_3 = 0 (for 3-cycles) at n=8. Need n=9 for alpha_3.
# alpha_2 = # disjoint 3-cycle pairs (uses 6 of 8 vertices)

for trial in range(10):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_c3_alpha(A, n)

    b0 = 1 - a1 + a2
    b1 = a1 - 2*a2
    b2 = a2

    H_approx = b0 + 3*b1 + 9*b2  # NOT correct H -- only 3-cycle contribution

    print(f"  trial {trial}: a1={a1:>3}, a2={a2:>4}, "
          f"b0={b0:>5}, b1={b1:>5}, b2={a2:>4}, "
          f"b1 {'> 0' if b1 > 0 else '<= 0'}")

# ============================================================
# PART 6: The meaning of b1 < 0 for the 3-adic tower
# ============================================================
print("\n" + "=" * 70)
print("PART 6: What b1 < 0 means for the 3-adic decomposition")
print("=" * 70)

print("""
H = b0 + 3*b1 + 9*b2

When b1 < 0:
  The "derivative layer" SUBTRACTS from H.
  H = (topology) + 3*(negative correction) + 9*(disjoint pairs)

  This means: at large n, the 3-adic decomposition has:
  - b0 = I(Omega, -1) << 0 (very negative topology)
  - 3*b1 < 0 (correction also negative)
  - 9*b2 >> 0 (disjoint pairs dominate)

  So H is DRIVEN by the curvature term 9*b2 = 9*alpha_2 at large n!

  The hierarchy INVERTS:
  Small n (n <= 7): H ~ b0 + 3*b1 (topology + derivative, curvature small)
  Large n (n >= 9): H ~ 9*b2 (curvature dominates, topology is noise)

  At the transition (n ~ 7-8): all three terms comparable.

This is a PHASE TRANSITION in the 3-adic structure of H!

The 2-adic structure has no such transition:
  H mod 2 = 1 (always, Redei)
  H mod 4 = 1 + 2*(a1 mod 2) (bit 1 = a1 parity, always works)

But the 3-adic structure changes character:
  Small n: H mod 3 tells you about TOPOLOGY (Euler char)
  Large n: H mod 3 still = b0 mod 3, but b0 is now a HUGE
           number with essentially random mod-3 residue.
           The topological content is "washed out" by the growth.
""")

# At n=7: what fraction of H comes from each layer?
print("  Layer contributions as fraction of H:")
n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

for trial in range(10):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_c3_alpha(A, n)
    H = 1 + 2*a1 + 4*a2  # OCF

    b0 = 1 - a1 + a2
    b1 = a1 - 2*a2
    b2 = a2

    # Check: H = b0 + 3*b1 + 9*b2
    assert H == b0 + 3*b1 + 9*b2

    pct_b0 = abs(b0) / H * 100
    pct_b1 = abs(3*b1) / H * 100
    pct_b2 = abs(9*b2) / H * 100

    # NOTE: at n=7, this uses only 3-cycles, so H is approximate
    # But the structure is clear
    print(f"  n=7: H={H:>5}, |b0|/H={pct_b0:>5.1f}%, "
          f"|3*b1|/H={pct_b1:>5.1f}%, |9*b2|/H={pct_b2:>5.1f}%")

print("\n  --- n=9 layer fractions ---")
n = 9
total_bits = n*(n-1)//2
np.random.seed(42)

for trial in range(10):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_c3_alpha(A, n)

    b0 = 1 - a1 + a2
    b1 = a1 - 2*a2
    b2 = a2

    H_c3 = b0 + 3*b1 + 9*b2  # Only 3-cycle contribution
    # H_c3 = 1 + 2*a1 + 4*a2

    if H_c3 > 0:
        pct_b0 = abs(b0) / H_c3 * 100
        pct_b1 = abs(3*b1) / H_c3 * 100
        pct_b2 = abs(9*b2) / H_c3 * 100
        print(f"  n=9: H_c3={H_c3:>5}, |b0|/H={pct_b0:>5.1f}%, "
              f"|3*b1|/H={pct_b1:>5.1f}%, |9*b2|/H={pct_b2:>5.1f}%, "
              f"b0={b0:>5}, b1={b1:>5}")

print("\n\nDone.")
