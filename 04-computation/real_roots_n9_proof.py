#!/usr/bin/env python3
"""
Proof attempt: I(Omega_3(T), x) has all real roots for n=9.

At n=9, I(Omega_3, x) has degree 3 (alpha(Omega_3) = 3 = floor(9/3)).
Write I = 1 + a1*x + a2*x^2 + a3*x^3 where:
  a1 = c3 = number of 3-cycle vertex sets
  a2 = number of pairs of vertex-disjoint 3-cycles
  a3 = number of triples of pairwise vertex-disjoint 3-cycles

For all real roots of a degree-3 polynomial with positive coefficients,
the discriminant must be non-negative:
  disc = 18*a1*a2*a3 - 4*a1^3*a3 + a1^2*a2^2 - 4*a2^3 - 27*a3^2 >= 0

Strategy:
1. Use Turan's theorem: the "disjoint pair graph" (3-cycles as vertices,
   edge iff disjoint) is triangle-free at n=9. So a2 <= a1^2/4.
2. Bound a3 in terms of a1 and a2.
3. Show disc >= 0 using these bounds.

Key combinatorial facts at n=9:
- a1 = c3 in range [0, C(9,3)/... ]. Max c3 is C(9,3) for completely
  non-transitive tournament. Actually c3 = C(n,3) - sum_v C(d+(v),2)
  by Moon's formula.
- Three pairwise disjoint 3-cycles use 9 vertices = ALL of n=9.
  So a3 = number of ways to partition [9] into three 3-cycles, each forming
  a directed 3-cycle. This is very constrained.
- a2 <= a1^2/4 by Turan (triangle-free graph).

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, random_tournament
from collections import defaultdict
from math import comb
import numpy as np

def find_3cycles(T):
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append((i,j,k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append((i,j,k))
    return cycles

def count_a2_a3(cycles):
    m = len(cycles)
    cycle_sets = [frozenset(c) for c in cycles]
    disjoint_after = [[] for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if not (cycle_sets[i] & cycle_sets[j]):
                disjoint_after[i].append(j)

    a2 = sum(len(d) for d in disjoint_after)

    a3 = 0
    for i in range(m):
        ui = cycle_sets[i]
        for j in disjoint_after[i]:
            uij = ui | cycle_sets[j]
            for k in disjoint_after[j]:
                if k > j and not (cycle_sets[k] & uij):
                    a3 += 1
    return a2, a3

print("=" * 70)
print("REAL-ROOTEDNESS PROOF ATTEMPT FOR n=9")
print("=" * 70)

# Part 1: Verify Turan bound and collect statistics
print("\n--- Part 1: Turan bound verification ---")
n = 9
samples = 2000
data = []

for _ in range(samples):
    T = random_tournament(n)
    c3_list = find_3cycles(T)
    a1 = len(c3_list)
    if a1 < 3:
        continue
    a2, a3 = count_a2_a3(c3_list)
    if a3 == 0:
        continue
    data.append((a1, a2, a3))

print(f"  Collected {len(data)} tournaments with a3 > 0")

# Verify Turan: a2 <= a1^2/4
turan_holds = sum(1 for a1,a2,a3 in data if a2 <= a1**2/4)
print(f"  Turan (a2 <= a1^2/4): {turan_holds}/{len(data)}")

# Part 2: What's the relationship between a3 and a1, a2?
print("\n--- Part 2: a3 bounds ---")
# a3 counts triples of pairwise disjoint 3-cycles.
# At n=9, three disjoint 3-cycles partition all 9 vertices.
# The number of such partitions depends on the tournament structure.

a3_over_a2 = [a3/a2 for a1,a2,a3 in data if a2 > 0]
a3_over_a1 = [a3/a1 for a1,a2,a3 in data]
a3_over_a1sq = [a3/a1**2 for a1,a2,a3 in data]

print(f"  a3/a2: min={min(a3_over_a2):.4f}, max={max(a3_over_a2):.4f}, mean={np.mean(a3_over_a2):.4f}")
print(f"  a3/a1: min={min(a3_over_a1):.4f}, max={max(a3_over_a1):.4f}, mean={np.mean(a3_over_a1):.4f}")
print(f"  a3/a1^2: min={min(a3_over_a1sq):.6f}, max={max(a3_over_a1sq):.6f}, mean={np.mean(a3_over_a1sq):.6f}")

# Part 3: Check discriminant with Turan bound
print("\n--- Part 3: Discriminant analysis ---")
# disc = 18*a1*a2*a3 - 4*a1^3*a3 + a1^2*a2^2 - 4*a2^3 - 27*a3^2

# With Turan: a2 <= a1^2/4, let's write a2 = t*a1^2/4 where 0 <= t <= 1
# and a3 = s*a1^3/... need to find the right parameterization

discs = []
for a1, a2, a3 in data:
    disc = 18*a1*a2*a3 - 4*a1**3*a3 + a1**2*a2**2 - 4*a2**3 - 27*a3**2
    discs.append(disc)

print(f"  Discriminant: min={min(discs)}, max={max(discs)}")
print(f"  Negative discriminants: {sum(1 for d in discs if d < 0)}")

# Part 4: Try to find the tightest case
print("\n--- Part 4: Tightest discriminant cases ---")
sorted_cases = sorted(zip(data, discs), key=lambda x: x[1])
for (a1, a2, a3), disc in sorted_cases[:10]:
    t = 4*a2 / a1**2 if a1 > 0 else 0  # Turan parameter
    print(f"  a1={a1}, a2={a2}, a3={a3}: disc={disc}, t=4a2/a1^2={t:.4f}")

# Part 5: Algebraic bound on a3
print("\n--- Part 5: Bounding a3 ---")
# a3 = number of partitions of [9] into three 3-cycles.
# For a 3-cycle on {a,b,c}: the tournament must have a directed 3-cycle
# on these vertices. Probability = 1/4 (uniform random tournament).
# Number of ways to partition [9] into three triples: C(9,3)*C(6,3)*C(3,3)/3! = 280
# Expected a3 for random tournament: 280 * (1/4)^3 = 280/64 = 4.375
# But conditional on each triple forming a 3-cycle: the triples are not independent.

# Let's compute the exact maximum a3 for n=9
max_a3 = 0
max_a3_tournament = None
for _ in range(5000):
    T = random_tournament(n)
    c3_list = find_3cycles(T)
    a1 = len(c3_list)
    if a1 < 3:
        continue
    _, a3 = count_a2_a3(c3_list)
    if a3 > max_a3:
        max_a3 = a3
        max_a3_tournament = [row[:] for row in T]

print(f"  Max a3 found in 5000 random n=9 tournaments: {max_a3}")
print(f"  280 total triple-partitions, so a3 <= 280")
print(f"  But constrained: 3-cycle probability ~1/4 per triple")
print(f"  Theoretical max: each of 84 triples forms a 3-cycle (regular tournament)")

# For the regular tournament on 9 vertices (if it exists):
# Each vertex has out-degree 4. c3 = C(9,3) - sum_v C(4,2) = 84 - 9*6 = 84 - 54 = 30.
# Wait, Moon's formula: c3 = C(n,3) - sum_v C(d+(v), 2)
# For regular: c3 = 84 - 9*C(4,2) = 84 - 54 = 30
print(f"  Regular tournament on 9: c3 = {comb(9,3)} - 9*{comb(4,2)} = {comb(9,3) - 9*comb(4,2)}")

# Part 6: Can we prove disc >= 0 analytically?
print("\n--- Part 6: Analytical bound attempt ---")
print("  disc = 18*a1*a2*a3 - 4*a1^3*a3 + a1^2*a2^2 - 4*a2^3 - 27*a3^2")
print("  = a1^2*a2^2 - 4*a2^3 + a3*(18*a1*a2 - 4*a1^3 - 27*a3)")
print()
print("  Factor 1: a1^2*a2^2 - 4*a2^3 = a2^2*(a1^2 - 4*a2)")
print("  By Turan: a2 <= a1^2/4, so a1^2 - 4*a2 >= 0. Factor 1 >= 0.")
print()
print("  Factor 2: 18*a1*a2 - 4*a1^3 - 27*a3")
print("  = a1*(18*a2 - 4*a1^2) - 27*a3")
print("  By Turan: a2 <= a1^2/4, so 18*a2 <= 18*a1^2/4 = 4.5*a1^2")
print("  So 18*a2 - 4*a1^2 <= 0.5*a1^2")
print("  Factor 2 <= 0.5*a1^3 - 27*a3")
print()
print("  So disc = Factor1 + a3*Factor2")
print("  >= 0 + a3*(0.5*a1^3 - 27*a3) ... wait, Factor2 can be negative")
print("  Need: a2^2*(a1^2 - 4*a2) >= a3*(4*a1^3 - 18*a1*a2 + 27*a3)")

# Let's check: is 4*a1^3 - 18*a1*a2 + 27*a3 positive or negative?
for a1, a2, a3 in data[:20]:
    f2 = 4*a1**3 - 18*a1*a2 + 27*a3
    print(f"  a1={a1}, a2={a2}, a3={a3}: 4a1^3-18a1a2+27a3 = {f2}")

# Part 7: Alternative approach — use a2/a1 and a3/a1 ratios
print("\n--- Part 7: Normalized analysis ---")
print("  Let r = a2/a1, s = a3/a1. Then:")
print("  disc/a1^2 = a1^2*r^2 - 4*a1*r^3 + 18*a1*r*s - 4*a1^2*s + a1^2*(r/a1)^2...")
print("  Actually, let's use t = a2/a1^2 and u = a3/a1^3.")
print("  disc/a1^6 = t^2 - 4*t^3/a1^2 + 18*t*u/a1 - 4*u + t^2*... ")
print("  This parameterization is messy. Let me try numerically.")

# For each data point, compute disc / a1^6
for (a1, a2, a3), disc in sorted_cases[:5]:
    t = a2 / a1**2
    u = a3 / a1**3
    print(f"  a1={a1}, t={t:.4f}, u={u:.6f}, disc={disc}, disc/a1^6={disc/a1**6:.6f}")

# Part 8: Direct bound
print("\n--- Part 8: Direct sufficient condition ---")
print("  Sufficient for disc >= 0:")
print("  a1^2*a2^2 >= 4*a2^3 + 4*a1^3*a3 + 27*a3^2 - 18*a1*a2*a3")
print("  LHS = a2^2*(a1^2)")
print("  RHS = 4*a2^3 + a3*(4*a1^3 + 27*a3 - 18*a1*a2)")
print()
print("  At n=9 with Turan (a2 <= a1^2/4):")
print("  If a3 is small enough relative to a1, the quadratic in a3 wins.")
print()

# The discriminant as a quadratic in a3:
# disc = -27*a3^2 + (18*a1*a2 - 4*a1^3)*a3 + (a1^2*a2^2 - 4*a2^3)
# = -27*a3^2 + B*a3 + C
# where B = 18*a1*a2 - 4*a1^3, C = a2^2*(a1^2 - 4*a2)
# This is non-negative when a3 <= (B + sqrt(B^2 + 108*C)) / 54

print("  disc as quadratic in a3: -27*a3^2 + B*a3 + C >= 0")
print("  where B = 18*a1*a2 - 4*a1^3, C = a2^2*(a1^2 - 4*a2)")
print("  Note: C >= 0 by Turan.")
print("  disc >= 0 iff a3 <= (B + sqrt(B^2 + 108*C)) / 54")

for (a1, a2, a3), disc in sorted_cases[:5]:
    B = 18*a1*a2 - 4*a1**3
    C = a2**2 * (a1**2 - 4*a2)
    bound = (B + np.sqrt(B**2 + 108*C)) / 54
    print(f"  a1={a1}, a2={a2}: B={B}, C={C}, a3_max={bound:.1f}, actual a3={a3}")

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
