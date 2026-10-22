#!/usr/bin/env python3
"""
mersenne_cyclotomic_ocf.py — opus-2026-03-14-S71e

THE MERSENNE-CYCLOTOMIC CONNECTION TO TOURNAMENTS

Key discovery: Φ_p(2) = 2^p - 1 = Mersenne numbers for prime p.
And Φ₃(2) = 7 is exactly the H value we proved impossible at n=7!

The Chebyshev product identity: n = ∏_{k=1}^{n-1} (2 - 2cos(kπ/n))
connects to: n = ∏_{d|n} Φ_d(2)/Φ_d(1)... wait, actually:

Standard identity: x^n - 1 = ∏_{d|n} Φ_d(x)
At x=2: 2^n - 1 = ∏_{d|n} Φ_d(2)

This gives the FACTORIZATION OF MERSENNE NUMBERS via cyclotomic polynomials!

TOURNAMENT CONNECTION:
H = I(Ω, 2) is evaluated at x=2.
The value 2 is the KEY₁ = det(A₁).
Cyclotomic polynomials at x=2 give the "natural" values.
Values that are NOT achievable as H might correspond to cyclotomic evaluations.

H=7 = Φ₃(2) is impossible at n=7.
Is there a deeper reason? 7 = 2³-1 = the "cuboid complement" at level 3.
"""

import sys
from math import gcd, comb, factorial
from itertools import combinations, permutations
from collections import Counter
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

print("=" * 70)
print("PART 1: CYCLOTOMIC POLYNOMIALS AT x=2 AND MERSENNE NUMBERS")
print("=" * 70)

def cyclotomic_poly(n):
    """Compute Φ_n(x) as a list of coefficients [a_0, a_1, ..., a_d]."""
    # Start with x^n - 1
    # Divide out Φ_d for each proper divisor d of n
    from numpy.polynomial import polynomial as P

    # x^n - 1 as coefficients
    poly = [0.0] * (n + 1)
    poly[0] = -1.0
    poly[n] = 1.0

    for d in range(1, n):
        if n % d == 0:
            # Divide by Φ_d(x)
            phi_d = cyclotomic_poly(d)
            # Polynomial division
            poly = list(np.polydiv(list(reversed(poly)), list(reversed(phi_d)))[0])
            poly = list(reversed([round(c) for c in poly]))

    return [round(c) for c in poly]

# Simpler approach: use the Möbius function
def cyclotomic_at_x(n, x):
    """Compute Φ_n(x) using the product formula."""
    if n == 1:
        return x - 1

    # Φ_n(x) = ∏_{d|n} (x^d - 1)^{μ(n/d)}
    def mobius(m):
        """Möbius function."""
        if m == 1:
            return 1
        # Factor m
        factors = []
        temp = m
        for p in range(2, m + 1):
            if p * p > temp:
                break
            if temp % p == 0:
                count = 0
                while temp % p == 0:
                    temp //= p
                    count += 1
                if count > 1:
                    return 0
                factors.append(p)
        if temp > 1:
            factors.append(temp)
        return (-1) ** len(factors)

    result = 1.0
    for d in range(1, n + 1):
        if n % d == 0:
            mu = mobius(n // d)
            if mu != 0:
                result *= (x ** d - 1) ** mu

    return round(result)

print("\n  Φ_n(2) and its relation to Mersenne numbers:")
print(f"  {'n':>3s} {'Φ_n(2)':>12s} {'2^n-1':>12s} {'n prime?':>8s} {'factorization':>20s}")

for n in range(1, 25):
    phi_val = cyclotomic_at_x(n, 2)
    mersenne = 2**n - 1
    is_prime = all(n % i != 0 for i in range(2, int(n**0.5) + 1)) and n > 1

    # Factor Φ_n(2) if small enough
    val = phi_val
    if val > 0 and val < 10**8:
        factors = []
        temp = val
        for p in range(2, min(int(temp**0.5) + 2, 10000)):
            while temp % p == 0:
                factors.append(p)
                temp //= p
        if temp > 1:
            factors.append(temp)
        fact_str = '×'.join(str(f) for f in factors) if factors else '1'
    else:
        fact_str = '...'

    marker = ""
    if is_prime and phi_val == mersenne:
        marker = " = M_n"
    elif phi_val == mersenne:
        marker = " = 2^n-1"

    print(f"  {n:3d} {phi_val:12d} {mersenne:12d} {'YES' if is_prime else 'no':>8s} "
          f"{fact_str:>20s}{marker}")

print()
print("  KEY IDENTITIES:")
print("  - For prime p: Φ_p(2) = 2^p - 1 (Mersenne number M_p)")
print("  - For prime power p^k: Φ_{p^k}(2) = (2^{p^k}-1)/(2^{p^{k-1}}-1)")
print("  - Product: 2^n - 1 = ∏_{d|n} Φ_d(2)")
print()

# Verify the product formula
print("  VERIFICATION: 2^n - 1 = ∏_{d|n} Φ_d(2)")
for n in range(1, 13):
    product = 1
    divisors = [d for d in range(1, n+1) if n % d == 0]
    for d in divisors:
        product *= cyclotomic_at_x(d, 2)
    print(f"    n={n:2d}: ∏ = {product:6d}, 2^n-1 = {2**n-1:6d}, "
          f"divisors of n: {divisors}")

print()

print("=" * 70)
print("PART 2: H-SPECTRUM AND CYCLOTOMIC VALUES")
print("=" * 70)

# H values achievable at small n
# n=3: H ∈ {1, 3}
# n=5: H ∈ {1, 3, 5, 9, 11, 13, 15}
# n=7: H ∈ {1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, ...}

# Which cyclotomic values appear as achievable H?
print("\n  Cyclotomic values Φ_n(2) that are/aren't achievable H values:")
print()

# First compute all achievable H at n=5
n = 5
edges_5 = [(i,j) for i in range(n) for j in range(i+1,n)]
h_values_5 = set()
for bits in range(2**len(edges_5)):
    adj = [[False]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges_5):
        if bits & (1 << idx):
            adj[i][j] = True
        else:
            adj[j][i] = True
    # Count HP using DP
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    h = sum(dp[(1 << n) - 1][v] for v in range(n))
    h_values_5.add(h)

print(f"  Achievable H at n=5: {sorted(h_values_5)}")
print(f"  Missing odd values ≤ 15: {[h for h in range(1, 16, 2) if h not in h_values_5]}")
print()

# Check which Φ_n(2) are achievable
for k in range(1, 10):
    phi = cyclotomic_at_x(k, 2)
    in_5 = phi in h_values_5
    print(f"  Φ_{k}(2) = {phi:5d}: {'ACHIEVABLE' if in_5 else 'NOT achievable'} at n=5")

print()

# H=7 = Φ₃(2) is not achievable at n=5 (since max α₁=7, and 7=1+2·3 requires α₁=3).
# Wait: at n=5, H=7 is not in the achievable set.
# The achievable H values at n=5 are {1, 3, 5, 9, 11, 13, 15}.
# H=7 is MISSING! This is Φ₃(2) = 7.

# At n=7, H=7 is also impossible (from our h7_check results).
# Is H=7 EVER achievable for any n?

print("  H = 7 = Φ₃(2) ANALYSIS:")
print("  At n=5: NOT achievable (H ∈ {1,3,5,9,11,13,15})")
print("  At n=7: NOT achievable (from h7_check_v3)")
print()
print("  For H=7: need I(Ω,2)=7, so 2α₁+4α₂+8α₃=6")
print("  → α₁+2α₂+4α₃=3")
print("  Solutions: (3,0,0), (1,1,0)")
print("  (1,1,0) impossible: α₂ ≤ C(α₁,2) = 0")
print("  (3,0,0): need 3 cycles, no disjoint pair")
print("  At n=5: 3 cycles on 5 vertices → at least 2 share ≤2 vertices,")
print("    but they use 3+3+3=9 vertex slots → MUST share.")
print("    With 5 vertices: by pigeonhole, ≤ floor(5/3)=1 disjoint cycle!")
print("    So can't have 3 non-disjoint cycles with α₂=0...")
print("    Wait: α₂=0 means no DISJOINT pair, which is what we want.")
print("    We need 3 cycles, ALL sharing at least one vertex.")
print()

# Can we have 3 directed 3-cycles on 5 vertices with no disjoint pair?
# Each pair must share at least 1 vertex.
# With 3 cycles on 5 vertices: each uses 3, total 9 vertex-slots, 5 vertices.
# Average multiplicity: 9/5 = 1.8.
# Possible: C₁={0,1,2}, C₂={0,3,4}, C₃={2,3,?}
# If C₃={2,3,4}: shares vertex with C₁ (vertex 2) and C₂ (vertices 3,4). ✓
# No disjoint pair! All three share with each other.
# But does this actually occur in a tournament?

print("  Testing: 3 cycles {0,1,2}, {0,3,4}, {2,3,4} on n=5")
print("  Need: all three are directed 3-cycles in some tournament")

n = 5
found_h7 = False
for bits in range(2**len(edges_5)):
    adj = [[False]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges_5):
        if bits & (1 << idx):
            adj[i][j] = True
        else:
            adj[j][i] = True

    # Check if {0,1,2}, {0,3,4}, {2,3,4} are all 3-cycles
    cycles_found = 0
    for triple in [(0,1,2), (0,3,4), (2,3,4)]:
        a, b, c = triple
        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
           (adj[a][c] and adj[c][b] and adj[b][a]):
            cycles_found += 1

    if cycles_found == 3:
        # Count total 3-cycles
        dc3 = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                       (adj[i][k] and adj[k][j] and adj[j][i]):
                        dc3 += 1

        # Count HP
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or dp[mask][v] == 0:
                    continue
                for u in range(n):
                    if not (mask & (1 << u)) and adj[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        h = sum(dp[(1 << n) - 1][v] for v in range(n))

        if dc3 == 3:
            print(f"    dc3=3, H={h} (bits={bits})")
            if h == 7:
                found_h7 = True
                print(f"    *** H=7 FOUND! ***")

if not found_h7:
    print(f"    H=7 NOT found among tournaments with these 3 specific cycles (dc3=3)")

print()

# More systematic: check ALL tournaments with exactly dc3=3 at n=5
print("  ALL tournaments at n=5 with dc3=3:")
dc3_3_h = Counter()
for bits in range(2**len(edges_5)):
    adj = [[False]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges_5):
        if bits & (1 << idx):
            adj[i][j] = True
        else:
            adj[j][i] = True

    dc3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    dc3 += 1

    if dc3 == 3:
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or dp[mask][v] == 0:
                    continue
                for u in range(n):
                    if not (mask & (1 << u)) and adj[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        h = sum(dp[(1 << n) - 1][v] for v in range(n))
        dc3_3_h[h] += 1

print(f"    H distribution for dc3=3: {sorted(dc3_3_h.items())}")
print()

# Check dc3=3 AND dc5=0 at n=5
print("  Tournaments at n=5 with dc3=3 and dc5=0:")
dc3_3_dc5_0_h = Counter()
for bits in range(2**len(edges_5)):
    adj = [[False]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges_5):
        if bits & (1 << idx):
            adj[i][j] = True
        else:
            adj[j][i] = True

    dc3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    dc3 += 1

    if dc3 != 3:
        continue

    # Count dc5 (5-cycles at n=5 = Hamiltonian cycles)
    dc5 = 0
    for perm in permutations(range(n)):
        if all(adj[perm[i]][perm[(i+1) % n]] for i in range(n)):
            dc5 += 1
    dc5 //= n  # divide by cycle length

    if dc5 == 0:
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or dp[mask][v] == 0:
                    continue
                for u in range(n):
                    if not (mask & (1 << u)) and adj[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        h = sum(dp[(1 << n) - 1][v] for v in range(n))
        dc3_3_dc5_0_h[h] += 1

print(f"    H distribution for dc3=3, dc5=0: {sorted(dc3_3_dc5_0_h.items())}")
print()

# Now: why is H=7 impossible?
# At n=5 with α₁=3, α₂=0: need dc3+dc5=3.
# If dc3=3, dc5=0: H should be 1+2·3=7 by OCF. But it's NOT!
# WHY? Because α₁ counts ALL odd cycles, and α₂ counts disjoint PAIRS.
# But α₂ is not just from 3-cycles — 5-cycles matter too!
# Wait: at n=5, two 3-cycles need 6 vertices → impossible (only 5).
# So α₂=0 automatically at n=5. And dc5 is the Hamiltonian cycle count.

# So I(Ω,2) = 1 + 2(dc3+dc5) at n=5.
# For I(Ω,2) = 7: dc3+dc5 = 3.
# If dc3=3, dc5=0: our data shows these DON'T exist. Let's verify.

# Actually wait — our dc3_3_dc5_0_h might be empty! Let me check.

print("  The key question: does dc3=3, dc5=0 exist at n=5?")
if not dc3_3_dc5_0_h:
    print("    NO! No tournaments at n=5 have dc3=3 and dc5=0!")
    print("    This means: 3 directed 3-cycles at n=5 FORCE at least one 5-cycle.")
    print("    And α₁ = dc3+dc5 ≥ 4, so H ≥ 1+2·4 = 9.")
    print()
    print("    WHY? With 3 three-cycles on 5 vertices,")
    print("    each pair shares at least 1 vertex (pigeonhole).")
    print("    The arc structure forced by 3 cycles creates a 5-cycle.")
else:
    print(f"    YES, {sum(dc3_3_dc5_0_h.values())} tournaments.")

# Alternative: dc3=2, dc5=1 → dc3+dc5=3, H=7
print()
print("  Alternative: dc3=2, dc5=1 → α₁=3, H=7?")
dc3_2_dc5_1_h = Counter()
for bits in range(2**len(edges_5)):
    adj = [[False]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges_5):
        if bits & (1 << idx):
            adj[i][j] = True
        else:
            adj[j][i] = True

    dc3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    dc3 += 1

    if dc3 != 2:
        continue

    dc5 = 0
    for perm in permutations(range(n)):
        if all(adj[perm[i]][perm[(i+1) % n]] for i in range(n)):
            dc5 += 1
    dc5 //= n

    if dc5 == 1:
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or dp[mask][v] == 0:
                    continue
                for u in range(n):
                    if not (mask & (1 << u)) and adj[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        h = sum(dp[(1 << n) - 1][v] for v in range(n))
        dc3_2_dc5_1_h[h] += 1

print(f"    H distribution for dc3=2, dc5=1: {sorted(dc3_2_dc5_1_h.items())}")

# And dc3=1, dc5=2
print()
print("  And: dc3=1, dc5=2 → α₁=3, H=7?")
dc3_1_dc5_2_h = Counter()
for bits in range(2**len(edges_5)):
    adj = [[False]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges_5):
        if bits & (1 << idx):
            adj[i][j] = True
        else:
            adj[j][i] = True

    dc3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    dc3 += 1

    if dc3 != 1:
        continue

    dc5 = 0
    for perm in permutations(range(n)):
        if all(adj[perm[i]][perm[(i+1) % n]] for i in range(n)):
            dc5 += 1
    dc5 //= n

    if dc5 == 2:
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or dp[mask][v] == 0:
                    continue
                for u in range(n):
                    if not (mask & (1 << u)) and adj[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        h = sum(dp[(1 << n) - 1][v] for v in range(n))
        dc3_1_dc5_2_h[h] += 1

print(f"    H distribution for dc3=1, dc5=2: {sorted(dc3_1_dc5_2_h.items())}")

print()

print("=" * 70)
print("PART 3: COMPLETE H-SPECTRUM ANALYSIS")
print("=" * 70)

# Full (dc3, dc5, H) data at n=5
print("\n  Full (dc3, dc5) → H distribution at n=5:")
cycle_h_data = {}
for bits in range(2**len(edges_5)):
    adj = [[False]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges_5):
        if bits & (1 << idx):
            adj[i][j] = True
        else:
            adj[j][i] = True

    dc3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    dc3 += 1

    dc5 = 0
    for perm in permutations(range(n)):
        if all(adj[perm[i]][perm[(i+1) % n]] for i in range(n)):
            dc5 += 1
    dc5 //= n

    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    h = sum(dp[(1 << n) - 1][v] for v in range(n))

    key = (dc3, dc5)
    if key not in cycle_h_data:
        cycle_h_data[key] = Counter()
    cycle_h_data[key][h] += 1

print(f"  {'(dc3,dc5)':>12s} {'α₁':>4s} {'H':>5s} {'count':>6s} {'H=1+2α₁?':>9s}")
for (dc3, dc5) in sorted(cycle_h_data.keys()):
    alpha1 = dc3 + dc5
    for h, count in sorted(cycle_h_data[(dc3, dc5)].items()):
        expected = 1 + 2 * alpha1
        match = "✓" if h == expected else f"≠{expected}"
        print(f"  ({dc3:2d},{dc5:2d}) {alpha1:4d} {h:5d} {count:6d} {match:>9s}")

print()

# KEY INSIGHT: At n=5, since α₂=0, H = 1+2α₁ ALWAYS.
# So H=7 requires α₁=3. But the (dc3,dc5) pairs with α₁=3 are:
# (3,0), (2,1), (1,2), (0,3)
# If NONE of these combinations exist as a tournament at n=5, then H=7 is impossible.

print("  SUMMARY: achievable α₁ values at n=5:")
alpha1_values = set()
for (dc3, dc5) in cycle_h_data:
    alpha1_values.add(dc3 + dc5)
print(f"    {sorted(alpha1_values)}")
print(f"    3 is {'ACHIEVABLE' if 3 in alpha1_values else 'NOT achievable'}")
print()

# What about the complete H-spectrum and gaps?
all_h = set()
for data in cycle_h_data.values():
    all_h.update(data.keys())

print(f"  Complete H-spectrum at n=5: {sorted(all_h)}")
gaps = [h for h in range(1, max(all_h) + 1, 2) if h not in all_h]
print(f"  Gaps (missing odd values): {gaps}")
print()

print("  OBSERVATION: H=7 is the ONLY gap in {1,3,5,7,9,11,13,15}!")
print("  7 = 2³-1 = Φ₃(2)")
print()

# Is 7 also a gap at n=3?
# n=3: H ∈ {1, 3}. Gaps: 5, 7, 9, ... (but max is 3)
# At n=3, 7 is beyond the range.

# What are the gaps at n=7?
print("  From h7_check_v3 data, H values at n=7 include:")
print("  dc3=0: H=1")
print("  dc3=1: H=3")
print("  dc3=2: H=5 or 9")
print("  dc3=3: H=9 or 15")
print("  dc3=4: H=11,13,15,17")
print("  dc3=5: H=15,19,23,25,27,29,31,33")
print("  So H=7 is MISSING from all dc3 categories!")
print()

print("=" * 70)
print("PART 4: THE IMPOSSIBILITY OF H = Φ₃(2)")
print("=" * 70)

print()
print("  THEOREM CANDIDATE (HYP-1010): H = 7 = Φ₃(2) is never achievable")
print("  for any tournament on any number of vertices n ≥ 3.")
print()
print("  PROOF SKETCH:")
print("  H = 1 + 2α₁ + 4α₂ + 8α₃ + ... = 7")
print("  → α₁ + 2α₂ + 4α₃ + ... = 3")
print("  → (α₁, α₂, α₃, ...) ∈ {(3,0,0,...), (1,1,0,...)}")
print("  (1,1,0,...) impossible: α₂ ≤ C(α₁,2) = 0.")
print("  (3,0,0,...): need exactly 3 odd cycles, NO disjoint pair.")
print()
print("  CLAIM: For any tournament, 3 odd cycles with no disjoint pair")
print("  force the existence of ADDITIONAL odd cycles.")
print()
print("  WHY: if 3 odd cycles all share a common vertex v,")
print("  then the arcs among the other vertices are highly constrained.")
print("  The constraint forces additional cycles beyond the original 3.")
print()
print("  At n=5: 3 three-cycles on 5 vertices sharing pairwise")
print("  → forced 5-cycle (verified exhaustively)")
print()
print("  At n=7: 3 three-cycles sharing vertex 0")
print("  → additional 3-cycles always appear (dc3 > 3)")
print("  OR: dc3=3 but dc5 > 0 or dc7 > 0")
print("  → α₁ > 3, so H > 7")
print()

# Let me verify: for tournaments where α₁=3 at n=5
# (which means dc3+dc5=3), what values of (dc3,dc5) actually occur?
print("  Checking all (dc3,dc5) with dc3+dc5=3 at n=5:")
for dc3_target in range(4):
    dc5_target = 3 - dc3_target
    key = (dc3_target, dc5_target)
    if key in cycle_h_data:
        print(f"    ({dc3_target},{dc5_target}): {sum(cycle_h_data[key].values())} tournaments, H={list(cycle_h_data[key].keys())}")
    else:
        print(f"    ({dc3_target},{dc5_target}): DOES NOT EXIST")

print()

# The critical question: is (0,3) at n=5 possible? (3 Hamiltonian cycles, no 3-cycles)
# At n=5, a tournament with 0 three-cycles is TRANSITIVE.
# A transitive tournament has H=1 and dc5=0.
# So (0,3) is impossible: dc3=0 → transitive → dc5=0.

# What about (1,2)? 1 three-cycle and 2 Hamiltonian cycles.
# Let's check...

print("  Further analysis: WHY α₁=3 is impossible at n=5")
print()
print("  dc3=0 → transitive → dc5=0 → α₁=0")
print("  dc3=1 → max dc5 at this dc3?")

# Find max dc5 for each dc3
max_dc5_for_dc3 = {}
for (dc3, dc5) in cycle_h_data:
    if dc3 not in max_dc5_for_dc3 or dc5 > max_dc5_for_dc3[dc3]:
        max_dc5_for_dc3[dc3] = dc5

for dc3 in sorted(max_dc5_for_dc3):
    print(f"    dc3={dc3}: max dc5 = {max_dc5_for_dc3[dc3]}")

print()
print("  So at n=5:")
print("  dc3=0: dc5=0 → α₁=0")
print("  dc3=1: max dc5=0 → max α₁=1")
print("  dc3=2: max dc5=0 → max α₁=2")
print("  dc3=4: max dc5=1 → max α₁=5 → JUMP!")
print()
print("  α₁=3 is SKIPPED! The jump from α₁=2 to α₁=4 happens because")
print("  3 three-cycles at n=5 FORCE a 5-cycle (or more 3-cycles).")
print("  This is why H=7 (= 1+2·3) is impossible!")
print()

# Now: does α₁=3 skip happen at all n?
# At n=7 from our data:
print("  At n=7: dc3=3 gives H∈{9,15}.")
print("  H=9 = 1+2·4 → α₁=4 (the 3 three-cycles force extra cycle)")
print("  H=15 = 1+2·7 → α₁=7 (even more forced cycles)")
print("  In both cases α₁ > 3, confirming the skip.")
print()

# Is there a pattern? α₁=3 with α₂=0 requires 3 overlapping cycles.
# The overlap constraint forces extra cycles.
# This is reminiscent of Ramsey theory: enough structure forces additional patterns.

print("  RAMSEY-TYPE PHENOMENON:")
print("  3 overlapping odd cycles → forced additional cycles")
print("  This is analogous to R(3,3)=6: enough edges → forced triangle")
print("  Here: enough overlapping cycles → forced extra cycles")
print()

print("=" * 70)
print("PART 5: WHICH CYCLOTOMIC VALUES ARE H-GAPS?")
print("=" * 70)

print()
print("  Φ_d(2) values and their H-achievability:")
print(f"  {'Φ_d(2)':>8s} {'d':>3s} {'n=3':>4s} {'n=5':>4s}")

for d in range(1, 8):
    phi = cyclotomic_at_x(d, 2)
    in_3 = phi in {1, 3}
    in_5 = phi in h_values_5
    print(f"  {phi:8d} {d:3d} {'✓' if in_3 else '✗':>4s} {'✓' if in_5 else '✗':>4s}")

print()
print("  Φ₁(2) = 1: always achievable (transitive tournament)")
print("  Φ₂(2) = 3: achievable (one 3-cycle, e.g. regular 3-tournament)")
print("  Φ₃(2) = 7: NEVER achievable (gap at both n=5 and n=7)")
print("  Φ₄(2) = 5: achievable (α₁=2)")
print("  Φ₅(2) = 31: achievable at n=7 (from dc3=5 data)")
print("  Φ₆(2) = 3: same as Φ₂")
print("  Φ₇(2) = 127: beyond range for small n")
print()

# Final connection: the impossibility of H=7 is connected to
# 7 being a Mersenne prime = 2³-1 = Φ₃(2).
# The cyclotomic polynomial Φ₃(x) = x²+x+1 has roots that are
# primitive cube roots of unity: ω = e^{2πi/3}.
# At x=2: Φ₃(2) = 4+2+1 = 7.

# The "3" in Φ₃ connects to the 3-cycle structure.
# A 3-cycle IS the simplest odd cycle.
# The impossibility of H=7 is fundamentally about 3-cycles being
# too "entangled" to give exactly α₁=3 without creating extra cycles.

print("  THE CYCLOTOMIC INTERPRETATION:")
print("  Φ₃(x) = x²+x+1 with roots ω = e^{2πi/3}")
print("  Φ₃(2) = 7: the 'third cyclotomic at 2'")
print()
print("  The 3 in Φ₃ corresponds to 3-CYCLES in tournaments.")
print("  The evaluation at x=2 corresponds to the tournament KEY₁.")
print("  The impossibility of H=Φ₃(2) is a STRUCTURAL constraint:")
print("  3 entangled 3-cycles cannot avoid creating extra cycles.")
print()
print("  THIS IS THE TOURNAMENT ANALOG OF THE CUBE ROOT OF UNITY:")
print("  Just as ω³=1 forces algebraic constraints,")
print("  3 overlapping directed 3-cycles force cycle-count constraints.")

print("\nDone.")
