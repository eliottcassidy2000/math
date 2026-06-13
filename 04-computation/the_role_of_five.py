#!/usr/bin/env python3
"""
the_role_of_five.py — opus-2026-03-14-S73

"Pay particular attention to 5."

5 is the BRIDGE number in tournament theory. It connects:
  - 2 and 3 (5 = 2+3)
  - Fibonacci and Jacobsthal (discriminant √5 vs √9)
  - 3-cycles and 7-cycles (5-cycles are the intermediate level)
  - The (2,3) pair and the (7,8) pair (5 sits between)
  - Binary: 101 (palindrome, = x²+1 at x=2)

PART 1:  5 = 2+3: the sum of the keys
PART 2:  √5 and the Fibonacci-Jacobsthal bridge
PART 3:  5 as the Galois distance at k=3
PART 4:  5-cycles: the second odd cycle
PART 5:  5 in binary: 101 as palindrome and x²+1
PART 6:  5 in the polynomial ring Z[x] at x=2
PART 7:  The pentic recurrence: 5-nacci and weighted 5-nacci
PART 8:  The role of 5 in cycle packing
PART 9:  5 and the discriminant ladder
PART 10: 5 as Fermat prime and its consequences
PART 11: The complete role of 5 in the hierarchy
"""

import numpy as np
import math
from itertools import combinations, permutations
from fractions import Fraction
import time

def adj_matrix(bits, n):
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

def hamiltonian_paths(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                count = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        count += dp.get((prev_mask, u), 0)
                if count:
                    dp[(mask, v)] = count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_directed_k_cycles(A, k):
    """Count directed k-cycles in tournament."""
    n = len(A)
    count = 0
    for verts in combinations(range(n), k):
        for perm in permutations(verts):
            if perm[0] != min(perm):
                continue
            if perm[1] > perm[-1]:
                continue
            valid = True
            for i in range(k):
                if not A[perm[i]][perm[(i+1) % k]]:
                    valid = False
                    break
            if valid:
                count += 1
    return count


print("=" * 70)
print("PART 1: 5 = 2 + 3 — THE SUM OF THE KEYS")
print("=" * 70)
print()

print("The number 5 is the SUM of the two fundamental primes.")
print("  5 = 2 + 3")
print()
print("In the Z[x] polynomial ring at x=2:")
print("  5 = x + (1+x) = 2x + 1 = 2·2 + 1")
print("  So 5 = 2x+1, which is the Galois distance at k=x=2:")
print("  2k-1 = 2·2-1 = 3... wait, that's 3, not 5.")
print()
print("  Actually: 5 = x² + 1 = 4 + 1 = 2² + 1")
print("  Also: 5 = x + (1+x) = the sum of the pair")
print()
print("  In the hierarchy {1, 2, 3, 4, 5, 6, 7, 8, ...}:")
print("  5 sits between (2,3) and (7,8)")
print("  It's the FIRST number not directly in the (x, 1+x) pair")
print("  but not yet at the x³ level either.")
print()
print("  5 = x²+1 = (1+x)² - (1+x) + 1 = Φ_3(1+x) where Φ_3 is cyclotomic!")
print("  Actually Φ_3(y) = y²-y+1. At y=3: 9-3+1 = 7, not 5.")
print("  Let me just note: 5 = x²+1 is a Fermat prime (2^{2^1}+1).")
print()

print("THE FIVE DECOMPOSITIONS OF 5:")
print("  5 = 2 + 3     (sum of keys)")
print("  5 = 2² + 1    (Fermat prime)")
print("  5 = 2·3 - 1   (product minus 1)")
print("  5 = 3² - 2²   (difference of squares of keys)")
print("  5 = 2³ - 3     (cube of first key minus second key)")
print()
print("  Each decomposition has meaning:")
print("  2+3: the combined weight")
print("  2²+1: x² plus identity")
print("  2·3-1: product of Jacobsthal levels minus identity")
print("  3²-2²: 9-4 = (1+x)² - x²")
print("  2³-3: x³ - (1+x) = 8-3 = the GAP from level 2 back to level 1")
print()

print("=" * 70)
print("PART 2: √5 AND THE FIBONACCI-JACOBSTHAL BRIDGE")
print("=" * 70)
print()

print("The Fibonacci recurrence f(n) = f(n-1) + f(n-2) [x=1]:")
print("  Roots: (1±√5)/2")
print("  φ = (1+√5)/2 ≈ 1.618")
print("  The discriminant is √(1+4·1) = √5")
print()
print("The Jacobsthal recurrence f(n) = f(n-1) + 2f(n-2) [x=2]:")
print("  Roots: (1±√9)/2 = (1±3)/2 = 2, -1")
print("  The discriminant is √(1+4·2) = √9 = 3 (RATIONAL!)")
print()
print("THE BRIDGE:")
print("  x=1: discriminant = √5   (irrational, involves 5)")
print("  x=2: discriminant = √9=3 (rational, involves 3)")
print()
print("  As x goes from 1 to 2:")
print("  discriminant = √(1+4x) goes from √5 to 3")
print("  At what x does the discriminant equal an ODD INTEGER?")
print()

# 1+4x = m² for odd m: x = (m²-1)/4
print("  1+4x = m² requires m odd: m=1,3,5,7,9,...")
print("  x = (m²-1)/4:")
print()
for m in range(1, 16, 2):
    x_val = (m*m - 1) // 4
    k = (m + 1) // 2
    print(f"    m={m:2d}: x = {x_val:3d},  roots = ({k}, {1-k}),  x = {k}·{k-1}")
print()
print("  So the integer-discriminant points are x = k(k-1): 0, 2, 6, 12, 20, ...")
print("  The discriminants are: 1, 3, 5, 7, 9, ...")
print("  And 5 IS THE FIBONACCI DISCRIMINANT!")
print()
print("  √5 = the discriminant at x=1 (Fibonacci)")
print("  3 = the discriminant at x=2 (Jacobsthal)")
print("  5 = the discriminant at x=6 (third level)")
print()
print("  Wait: at x=6, disc = √(1+24) = √25 = 5 ✓")
print("  So √5 appears TWICE: as disc(Fibonacci) and as disc(x=6 level)!")
print()
print("  The BRIDGE: 5 connects x=1 (Fibonacci) to x=6 (third level)")
print("  through the discriminant sequence.")
print()

# The discriminant at each x = k(k-1)
print("DISCRIMINANT TABLE:")
print(f"  {'k':>3s}  {'x=k(k-1)':>8s}  {'disc=2k-1':>10s}  {'name':>20s}")
print("  " + "-" * 45)
names = {1: "trivial", 2: "Jacobsthal/OCF", 3: "6-spectrum",
         4: "12-spectrum", 5: "20-level", 6: "30-level"}
for k in range(1, 8):
    x_val = k*(k-1)
    disc = 2*k-1
    name = names.get(k, "")
    print(f"  {k:3d}  {x_val:8d}  {disc:10d}  {name:>20s}")
print()
print("  The discriminants ARE the odd numbers: 1, 3, 5, 7, 9, ...")
print("  And 5 sits at position k=3, x=6 — the THIRD level.")
print("  The first truly nontrivial discriminant after 3.")
print()

print("=" * 70)
print("PART 3: 5 AS THE GALOIS DISTANCE AT k=3")
print("=" * 70)
print()

print("At x=6 (k=3): roots are 3 and -2")
print("  Galois distance = 3-(-2) = 5")
print("  J₆(n) = (3^n - (-2)^n) / 5")
print()
print("  So 5 is the DENOMINATOR of the second Jacobsthal sequence.")
print("  This means 5 | (3^n - (-2)^n) for all n.")
print()

# Verify
print("  Verification: 3^n - (-2)^n mod 5")
for n in range(1, 12):
    val = 3**n - (-2)**n
    print(f"    n={n:2d}: 3^n-(-2)^n = {val:8d}, mod 5 = {val % 5}")
print()

print("  The sequence J₆(n) = (3^n-(-2)^n)/5:")
for n in range(1, 12):
    j6 = (3**n - (-2)**n) // 5
    print(f"    J₆({n:2d}) = {j6:8d}")
print()
print("  J₆ starts: 1, 1, 7, 13, 55, 133, 463, 1261, ...")
print("  Note J₆(3) = 7! The Mersenne prime at the (7,8) threshold!")
print()
print("  So 5 (via J₆) GENERATES 7, which is the next level's threshold.")
print("  The chain: 5 → J₆(3) = 7 → (7,8) transition")
print()

print("=" * 70)
print("PART 4: 5-CYCLES — THE SECOND ODD CYCLE")
print("=" * 70)
print()

print("In tournaments:")
print("  3-cycles: appear at n ≥ 3 (smallest non-transitive structure)")
print("  5-cycles: appear at n ≥ 5 (first 'higher' cycle)")
print("  7-cycles: appear at n ≥ 7 (first 'cross-level' interaction)")
print()

# Compute 5-cycle statistics for n=5 (the critical case)
n = 5
total = 2**(n*(n-1)//2)
print(f"n=5: All {total} tournaments")
print()

dc3_list = []
dc5_list = []
H_list = []
for bits in range(total):
    A = adj_matrix(bits, n)
    dc3 = count_directed_k_cycles(A, 3)
    dc5 = count_directed_k_cycles(A, 5)
    H = hamiltonian_paths(A)
    dc3_list.append(dc3)
    dc5_list.append(dc5)
    H_list.append(H)

print("  5-cycle distribution:")
from collections import Counter
dc5_dist = Counter(dc5_list)
for k in sorted(dc5_dist.keys()):
    print(f"    dc5={k}: {dc5_dist[k]:4d} tournaments ({dc5_dist[k]/total*100:.1f}%)")

print()
print("  KEY: at n=5, a 5-cycle uses ALL 5 vertices.")
print("  So a 5-cycle at n=5 is a HAMILTONIAN CYCLE.")
print("  And dc5 = #{Hamiltonian cycles}/2 (each undirected cycle = 2 directed)")
print()
print("  dc5=0: 480 tournaments (no Hamiltonian cycle)")
print("  dc5≥1: 544 tournaments (have Hamiltonian cycle)")
print()

# Relationship between H and dc5
print("  Correlation between H and dc5:")
arr_H = np.array(H_list, dtype=float)
arr_dc5 = np.array(dc5_list, dtype=float)
arr_dc3 = np.array(dc3_list, dtype=float)
corr_H_dc5 = np.corrcoef(arr_H, arr_dc5)[0, 1]
corr_H_dc3 = np.corrcoef(arr_H, arr_dc3)[0, 1]
print(f"    corr(H, dc5) = {corr_H_dc5:.6f}")
print(f"    corr(H, dc3) = {corr_H_dc3:.6f}")
print()

# The formula H = 1 + 2α₁ + 4α₂
# α₁ = dc3 + dc5 (total directed odd cycles)
# α₂ = #{disjoint pairs} = 0 at n=5 (can't fit two disjoint 3-cycles in 5 vertices)
print("  At n=5: α₁ = dc3 + dc5 (total directed odd cycles)")
print("  α₂ = 0 always (can't fit two disjoint cycles in 5 vertices)")
print("  H = 1 + 2(dc3 + dc5)")
print()

# Verify
mismatches = 0
for i in range(total):
    a1 = dc3_list[i] + dc5_list[i]
    H_check = 1 + 2 * a1
    if H_check != H_list[i]:
        mismatches += 1
print(f"  H = 1 + 2(dc3+dc5) verified: {mismatches} mismatches out of {total}")
print()

print("  THE 5-CYCLE CONTRIBUTION:")
print("  For each directed 5-cycle, H increases by 2.")
print("  For each directed 3-cycle, H increases by 2.")
print("  At n=5, 3-cycles and 5-cycles are ON EQUAL FOOTING.")
print("  Each contributes +2 to H, and they don't interfere")
print("  (can't be disjoint at n=5, so α₂=0).")
print()

print("  But at n=8: 3-cycles and 5-cycles CAN be disjoint (3+5=8).")
print("  This is the cross-level coupling, and 5 is essential to it.")
print("  Without 5-cycles, the first cross-level would be at n=6 (3+3).")
print("  But 5-cycles delay the INTERESTING cross-level to n=8.")
print()

print("=" * 70)
print("PART 5: 5 IN BINARY — 101 AS PALINDROME")
print("=" * 70)
print()

print("Binary representations of key numbers:")
for num in [1, 2, 3, 4, 5, 6, 7, 8, 10, 11]:
    b = bin(num)[2:]
    is_palindrome = b == b[::-1]
    print(f"  {num:3d} = {b:>5s}  {'(palindrome)' if is_palindrome else ''}")

print()
print("  5 = 101₂ is a PALINDROME in binary!")
print("  The only palindromes in {1,...,11}: 1, 3, 5, 7")
print("  These are exactly the ODD NUMBERS ≤ 7!")
print("  And they are exactly the CYCLE LENGTHS in tournaments ≤ 7.")
print()
print("  Palindromes in binary = numbers of the form:")
print("  1, 11, 101, 111, 1001, 1111, ...")
print("  = 1, 3, 5, 7, 9, 15, 21, 27, ...")
print()
print("  5 = 101₂: the FIRST non-trivial binary palindrome.")
print("  It has the structure '1_1' — a frame around a gap.")
print("  Compare: 10 = 1010₂ (interleaved 2s)")
print("           11 = 1011₂ (interleaved 2 and 3)")
print("           5  = 101₂  (a single 2 framed by 1s)")
print()

print("=" * 70)
print("PART 6: 5 IN THE POLYNOMIAL RING Z[x] AT x=2")
print("=" * 70)
print()

print("Multiple representations of 5 in Z[x]:")
print(f"  x²+1 = {2**2+1} = 5  ✓")
print(f"  x+(1+x) = {2+(1+2)} = 5  ✓  (sum of pair)")
print(f"  x(1+x)-1 = {2*3-1} = 5  ✓  (product minus 1)")
print(f"  (1+x)²-x² = {3**2-2**2} = 5  ✓  (difference of squares)")
print(f"  x³-(1+x) = {8-3} = 5  ✓  (cube minus next key)")
print()

print("  The 'simplest' representation: 5 = x² + 1")
print("  This makes 5 the NORM of the Gaussian integer 2+i:")
print("  |2+i|² = 2²+1² = 5")
print()
print("  In the Eisenstein integers Z[ω] (ω = e^{2πi/3}):")
print("  5 = (2+ω)(2+ω²) = |2+ω|²")
print("  where ω = (-1+√3i)/2, so 2+ω = (3+√3i)/2")
print("  |2+ω|² = (3²+3)/4 = 3  (NOT 5 in Eisenstein norm)")
print()
print("  Actually in Eisenstein: 5 is inert (doesn't split)")
print("  because 5 ≡ 2 (mod 3), not ≡ 1.")
print("  So 5 is a PRIME in Z[ω] but splits in Z[i].")
print()
print("  This means: 5 'knows about' i=√(-1) (Gaussian) but")
print("  is BLIND to ω=cube root of unity (Eisenstein).")
print("  In tournament theory: 5 relates to the x²+1 structure")
print("  (mod 4 residues) but not to the cubic (mod 3) structure.")
print()

print("=" * 70)
print("PART 7: THE 5-NACCI AND WEIGHTED 5-NACCI")
print("=" * 70)
print()

def knacci_root(k, iters=200):
    x = 2.0 - 0.5**(k-1)
    for _ in range(iters):
        f = x**(k+1) - 2*x**k + 1
        fp = (k+1)*x**k - 2*k*x**(k-1)
        if abs(fp) < 1e-30: break
        x -= f/fp
    return x

def weighted_knacci_root(k, iters=200):
    x = 3.0 - (2/3)**(k-1)
    for _ in range(iters):
        val = x**(k+1) - 3*x**k + 2**k
        dval = (k+1)*x**k - 3*k*x**(k-1)
        if abs(dval) < 1e-30: break
        x -= val/dval
    return x

phi5 = knacci_root(5)
psi5 = weighted_knacci_root(5)

print(f"5-nacci dominant root: φ₅ = {phi5:.10f}")
print(f"  Gap to 2: {2-phi5:.10e}")
print(f"  ≈ 1/2⁵ = {1/32:.10e}")
print()
print(f"Weighted 5-nacci dominant root: ψ₅ = {psi5:.10f}")
print(f"  Gap to 3: {3-psi5:.10e}")
print(f"  ≈ (2/3)⁵ = {(2/3)**5:.10e}")
print()
print(f"  Ratio ψ₅/φ₅ = {psi5/phi5:.10f} (→ 3/2 = 1.5)")
print()

# The 5-nacci sequence
print("5-nacci sequence: f(n) = f(n-1) + f(n-2) + f(n-3) + f(n-4) + f(n-5)")
f = [0]*5 + [1]  # f(0)=...=f(3)=0, f(4)=0, f(5)=1...
# Standard: f(1)=...=f(4)=1, f(5)=f(1)+...+f(5)
f = [1, 1, 2, 4, 8]
for i in range(15):
    f.append(sum(f[-5:]))
print(f"  {f[:15]}")
print()

# Ratios
print("  Consecutive ratios → φ₅:")
for i in range(5, 14):
    r = f[i] / f[i-1] if f[i-1] > 0 else 0
    print(f"    f({i})/f({i-1}) = {f[i]}/{f[i-1]} = {r:.8f}")
print()

print("  The 5-nacci ratio stabilizes to φ₅ ≈ 1.9659")
print("  which is already 98.3% of the way to 2.")
print("  By k=5, the k-nacci is 'almost Jacobsthal'.")
print()

print("=" * 70)
print("PART 8: 5 AND CYCLE PACKING — THE CRITICAL ROLE")
print("=" * 70)
print()

print("CYCLE PACKING THRESHOLDS:")
print("  n=3: first 3-cycle (3)")
print("  n=5: first 5-cycle (5)")
print("  n=6: first disjoint pair of 3-cycles (3+3)")
print("  n=7: first 7-cycle (7)")
print("  n=8: first mixed pair (3+5)")
print("  n=9: first triple of 3-cycles (3+3+3)")
print("  n=10: first pair of 5-cycles (5+5)")
print("  n=11: first triple including 5 (3+3+5)")
print()
print("  The 5-cycle creates NEW packing types at:")
print("  n=5: α₁ type (5)")
print("  n=8: α₂ type (3,5)  ← THE (7,8) TRANSITION!")
print("  n=10: α₂ type (5,5)")
print("  n=11: α₃ type (3,3,5)")
print("  n=13: α₂ type (5,8) ... no, 8 is even. α₂ type (5,7) at n=12")
print()
print("  Without 5-cycles, the packing structure would be:")
print("  Only 3-cycles: n=3 (single), n=6 (pair), n=9 (triple), ...")
print("  Only 7-cycles: n=7, ...")
print("  The 5-cycle FILLS THE GAP between n=6 and n=7.")
print()

# Count tournaments with specific packing types at n=7
n = 7
print(f"n={n}: 5-cycle role in cycle structure (200 random samples)")
import random
random.seed(73)
t0 = time.time()

has_5cycle = 0
has_disjoint_35 = 0
total_samples = 200

for trial in range(total_samples):
    bits = random.randint(0, 2**(n*(n-1)//2) - 1)
    A = adj_matrix(bits, n)
    dc5 = count_directed_k_cycles(A, 5)
    if dc5 > 0:
        has_5cycle += 1
    # Check for disjoint 3+5 pair: impossible at n=7 (3+5=8>7)... wait, vertex-disjoint needs 8 vertices
    # Actually at n=7: disjoint (3,3) is possible (uses 6 of 7), but (3,5) needs 8 > 7
    # So α₂^{35} = 0 at n=7!

print(f"  {has_5cycle}/{total_samples} have at least one directed 5-cycle ({has_5cycle/total_samples*100:.1f}%)")
print(f"  Time: {time.time()-t0:.1f}s")
print()
print("  At n=7: disjoint (3,5) pair IMPOSSIBLE (3+5=8>7)")
print("  So 5-cycles contribute to α₁ but NOT to α₂ at n≤7.")
print("  The 5-cycle's contribution to independence structure")
print("  only becomes 'interesting' (pairwise) at n≥8.")
print()
print("  This is why n=8 (= 2³) is the cross-level threshold:")
print("  it's the SMALLEST n where a 3-cycle and 5-cycle can coexist disjointly.")
print()

print("=" * 70)
print("PART 9: 5 AND THE DISCRIMINANT LADDER")
print("=" * 70)
print()

print("The discriminant of f(n) = f(n-1) + x·f(n-2) is Δ = 1+4x.")
print()
print("  At x=1: Δ = 5   → Fibonacci (the √5 recurrence)")
print("  At x=2: Δ = 9   → Jacobsthal (the √9 = 3 recurrence)")
print("  At x=6: Δ = 25  → 6-spectrum (the √25 = 5 recurrence!)")
print("  At x=12: Δ = 49 → 12-spectrum (the √49 = 7 recurrence)")
print()
print("THE LADDER:")
print("  Δ = 5, 9, 25, 49, 81, 121, ...")
print("  = 5, 3², 5², 7², 9², 11², ...")
print()
print("  Wait: at x=1, Δ=5 is NOT a perfect square!")
print("  That's EXACTLY why Fibonacci has IRRATIONAL roots.")
print("  And for x=2,6,12,..., Δ IS a perfect square → rational roots.")
print()
print("  So 5 is the discriminant that separates IRRATIONAL (Fibonacci)")
print("  from RATIONAL (Jacobsthal and beyond).")
print()
print("  Fibonacci: Δ=5, roots involve √5 (irrational)")
print("  Jacobsthal: Δ=9=3², roots are 2,-1 (rational)")
print()
print("  The NEXT non-square discriminant: at x=3, Δ=13.")
print("  And 13 is prime, so √13 is irrational too.")
print("  At x=5: Δ=21 = 3·7, irrational.")
print("  At x=10: Δ=41, prime, irrational.")
print("  At x=11: Δ=45 = 9·5, so √45 = 3√5. Still involves √5!")
print()
print("  REMARKABLE: x=11 gives discriminant 45 = 9·5 = 3²·5")
print("  So √45 = 3√5, meaning the roots at x=11 are:")
print("  (1 ± 3√5)/2")
print("  The x=11 recurrence has roots in the SAME FIELD Q(√5)")
print("  as Fibonacci!")
print()

# The x values where discriminant involves √5
print("  Which x values have discriminant involving √5?")
print("  Need 1+4x = 5·m² for some integer m.")
print("  4x = 5m²-1, x = (5m²-1)/4")
for m in range(1, 8):
    val = 5*m*m - 1
    if val % 4 == 0:
        x_val = val // 4
        disc = 1 + 4*x_val
        sqrt_disc = math.sqrt(disc)
        print(f"    m={m}: x = {x_val}, Δ = {disc} = 5·{m}², √Δ = {m}√5 = {sqrt_disc:.6f}")
    else:
        print(f"    m={m}: x = {val}/4 (not integer)")
print()
print("  So: x=1 (Δ=5), x=6 (Δ=25=5·1²... wait, 25=5²)")
print("  Hmm, x=6: Δ=25=5², so √25=5, which is RATIONAL.")
print("  x=11: Δ=45=5·9=5·3², √45=3√5.")
print("  x=31: Δ=125=5³=5·25, √125=5√5.")
print()
print("  x=1: √5 (m=1)")
print("  x=11: 3√5 (m=3)")
print("  x=31: 5√5 (m=5)")
print()
print("  The x-values with √5 discriminant: x = (5m²-1)/4 for ODD m.")
print("  m=1: x=1 (Fibonacci)")
print("  m=3: x=11 (the 'decimal pair' number!)")
print("  m=5: x=31 (a Mersenne prime!)")
print()
print("  DISCOVERY: x=1 (Fibonacci) and x=11 share the √5 field!")
print("  The roots at x=11 are (1±3√5)/2, which is φ-related:")
print(f"  φ = (1+√5)/2 = {(1+math.sqrt(5))/2:.6f}")
print(f"  (1+3√5)/2 = {(1+3*math.sqrt(5))/2:.6f} = 3φ + (-1)")
print(f"  Since 3φ = 3·{(1+math.sqrt(5))/2:.6f} = {3*(1+math.sqrt(5))/2:.6f}")
print(f"  And 3φ-1 = {3*(1+math.sqrt(5))/2-1:.6f}")
print(f"  (1+3√5)/2 = {(1+3*math.sqrt(5))/2:.6f} ← close!")
print()

# Exact: (1+3√5)/2 = (1+3√5)/2. And 3φ = (3+3√5)/2.
# 3φ - 1 = (3+3√5)/2 - 1 = (1+3√5)/2 ✓
print("  EXACT: (1+3√5)/2 = 3φ - 1 = 3·(1+√5)/2 - 1")
print("  So root(x=11) = 3φ - 1  ← beautiful!")
print()
print("  AND the conjugate: (1-3√5)/2 = 3(-1/φ) - 1 + 1 = ...")
print("  (1-3√5)/2 = -(3√5-1)/2 = -(3φ-2) = 2-3φ")
print("  Conjugate(x=11) = 2 - 3φ = 2 - 3·1.618... = -2.854...")
print()
print("  So at x=11:")
print("    r = 3φ-1 ≈ 3.854")
print("    s = 2-3φ ≈ -2.854")
print("    r·s = (3φ-1)(2-3φ) = 6φ-9φ²-2+3φ = 9φ-9φ²-2")
print("    = 9φ(1-φ)-2 = 9φ(-1/φ²)-2... let me just compute:")
print(f"    r·s = {(3*(1+math.sqrt(5))/2-1) * (2-3*(1+math.sqrt(5))/2):.6f}")
print(f"    Should be -x = -11: {-11}")
print()

print("=" * 70)
print("PART 10: 5 AS A FERMAT PRIME")
print("=" * 70)
print()

print("5 = 2^{2^1} + 1 is a Fermat prime (F₁).")
print("The Fermat primes are: 3, 5, 17, 257, 65537, ...")
print("  F₀ = 3 = 2^1 + 1 (the other key!)")
print("  F₁ = 5 = 2^2 + 1")
print("  F₂ = 17 = 2^4 + 1")
print("  F₃ = 257 = 2^8 + 1")
print()
print("So both 3 and 5 are FERMAT PRIMES.")
print("  3 = F₀ = 2^{2^0} + 1 = 1+x at x=2")
print("  5 = F₁ = 2^{2^1} + 1 = x²+1 at x=2")
print()
print("  The pattern: F_k = x^{2^k} + 1")
print("  F₀ = x+1 = 3")
print("  F₁ = x²+1 = 5")
print("  F₂ = x⁴+1 = 17")
print()
print("  In the Z[x] ring: Fermat primes are x^{2^k}+1.")
print("  These are special because they relate to CONSTRUCTIBILITY")
print("  (regular polygon constructibility requires Fermat primes).")
print()
print("  In tournament theory: n-gons relate to circulant tournaments")
print("  on Z_n. A regular n-gon is constructible iff n is a product")
print("  of a power of 2 and distinct Fermat primes.")
print("  Constructible n: 3, 4, 5, 6, 8, 10, 12, 15, 16, 17, ...")
print("  Note: 5, 10, 15 are ALL constructible (involving F₁=5).")
print()
print("  The 5-cycle in a tournament on Z_5 (the only cyclic tournament")
print("  on 5 vertices) corresponds to the regular pentagon — the")
print("  simplest Fermat-prime polygon beyond the triangle (F₀=3).")
print()

print("=" * 70)
print("PART 11: THE COMPLETE ROLE OF 5 IN THE HIERARCHY")
print("=" * 70)
print()

print("5 plays SIX distinct roles in tournament theory:")
print()
print("1. SUM OF KEYS: 5 = 2+3 = x+(1+x)")
print("   The combined weight of counting + topology.")
print()
print("2. FIBONACCI DISCRIMINANT: Δ(x=1) = 1+4 = 5")
print("   The 'irrationality' that makes φ transcendental-looking.")
print("   At x=2: Δ=9=3², eliminating the irrationality.")
print("   5 is the OBSTRUCTION that Jacobsthal overcomes.")
print()
print("3. SECOND CYCLE LENGTH: 5-cycles in tournaments")
print("   First appears at n=5, enables cross-level at n=8.")
print("   The 5-cycle is to 3-cycle as 3 is to 1:")
print("   the first non-trivial extension.")
print()
print("4. JACOBSTHAL DENOMINATOR: J₆(n) = (3^n-(-2)^n)/5")
print("   The Galois distance at k=3, connecting the")
print("   x=6 level to cycle packing.")
print()
print("5. FERMAT PRIME: 5 = 2²+1 = F₁")
print("   Relates to constructibility and circulant structure.")
print("   The regular pentagon (5-gon) is the simplest")
print("   non-trivial constructible polygon.")
print()
print("6. √5 FIELD BRIDGE: Q(√5) contains BOTH Fibonacci roots")
print("   and x=11 roots (via 3φ-1).")
print("   So 5 connects the (2,3) level to the (10,11) level")
print("   through the number field Q(√5).")
print()
print("THE HIERARCHY WITH 5:")
print()
print("  Level 0: 1 (identity)")
print("  Level 1: (2, 3) — the two keys")
print("  Level 1.5: 5 = 2+3 — the BRIDGE")
print("  Level 2: (7, 8) — the transition (= 2³-1, 2³)")
print("  Level 3: (10, 11) — the interleaving")
print()
print("  5 sits at 'Level 1.5' — halfway between the fundamental")
print("  pair and the transition pair.")
print()
print("  In the recurrence tower:")
print("  k=2: root → φ ≈ 1.618    (Fibonacci, discriminant 5)")
print("  k=3: root → 1.839...      (tribonacci)")
print("  k=4: root → 1.928...      (tetranacci)")
print("  k=5: root → 1.966...      (pentanacci — the '5-nacci')")
print(f"       gap to 2: {2-knacci_root(5):.6e}")
print()
print("  The 5-nacci is 98.3% of the way from 1 to 2.")
print("  Only 1.7% of the asymptotic gap remains.")
print()
print("  And in the weighted tower:")
print(f"  ψ₅ = {weighted_knacci_root(5):.10f}")
print(f"  gap to 3: {3-weighted_knacci_root(5):.6e}")
print(f"  This is (2/3)^5 = {(2/3)**5:.6e} of the total gap.")
print()

# Final synthesis
print("=" * 70)
print("FINAL: WHY 5 MATTERS")
print("=" * 70)
print()
print("  Without 5:")
print("  - Fibonacci would have rational roots (no golden ratio)")
print("  - The cross-level transition would happen at n=6 (3+3), not n=8")
print("  - The Jacobsthal denominator sequence would skip from 3 to 7")
print("  - The x=11 level would not connect back to Fibonacci via Q(√5)")
print("  - There would be no regular pentagon")
print()
print("  With 5:")
print("  - 5 = 2+3 bridges the counting and topology keys")
print("  - √5 makes Fibonacci irrational, creating the golden ratio")
print("  - 5-cycles delay cross-level coupling to n=8=2³")
print("  - J₆(3) = 7, generating the transition threshold")
print("  - The number field Q(√5) links Fibonacci to x=11")
print("  - F₁=5 enables pentagon constructibility")
print()
print("  5 is the MORTAR between the bricks (2,3).")
print("  You need 2 and 3 to build. You need 5 to hold them together.")
print()
print("Done.")
