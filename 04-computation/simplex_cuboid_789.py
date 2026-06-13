#!/usr/bin/env python3
"""
simplex_cuboid_789.py — opus-2026-03-14-S71e

HOW 7, 8, 9, 10 FIT IN THE SIMPLEX-CUBOID FRAMEWORK

The user's framing:
  Simplices = (x+1)^n → evaluate at x=1 gives 2^n
  Cuboids = (x+2)^n → evaluate at x=1 gives 3^n
  k-nacci → 2 (simplex limit)
  weighted/doubled k-nacci → 3 (cuboid limit)

KEY QUESTION: What happens at 7 and 8? At 9=3² and 10=5*2?

Connection to tournament theory:
  n=7: QR_7 tournament, first β_3>0 homology
  n=8: first even tournament size with rich structure
  n=9: 3² — first squared cycle length, α_3 first appears
  n=10: 5*2 — first time 5-cycles can pair up

The (x+1)/(x+2) framework at these special values of n.
"""

import sys
import numpy as np
from math import comb, factorial
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

print("=" * 70)
print("PART 1: THE (x+1)^n vs (x+2)^n RATIO AT SPECIAL n")
print("=" * 70)

print("\nThe ratio (x+1)^n / (x+2)^n as a function of x:")
print("  At x=1: (2/3)^n → 0 as n → ∞")
print("  At x=0: (1/2)^n → 0")
print("  At x=-1: (0/1)^n = 0 for n>0")
print("  At x=∞: → 1")
print()

for n in [3, 5, 7, 8, 9, 10, 11, 12, 15]:
    r = (2/3)**n
    print(f"  n={n:2d}: (2/3)^n = {r:.8f}, diff 3^n - 2^n = {3**n - 2**n:>10d}")

print("\n" + "=" * 70)
print("PART 2: DIFFERENCE (x+2)^n - (x+1)^n AT EACH n")
print("=" * 70)

print("""
  (x+2)^n - (x+1)^n evaluated at specific x values:

  x=0: 2^n - 1^n = 2^n - 1 (Mersenne numbers!)
  x=1: 3^n - 2^n
  x=2: 4^n - 3^n
  x=3: 5^n - 4^n

  At x=0, these give Mersenne numbers:
    n=2: 3, n=3: 7, n=5: 31, n=7: 127 — ALL MERSENNE PRIMES!
    n=4: 15, n=6: 63 — not prime.

  At x=1:
    3^n - 2^n: 1, 5, 19, 65, 211, 665, 2059, 6305, 19171, ...
    (OEIS A001047: 3^n - 2^n)

  TOURNAMENT CONNECTION:
  I(Omega, 3) - I(Omega, 2) = alpha_1 + 5*alpha_2 + 19*alpha_3 + ...
  The coefficients ARE the Mersenne-like numbers 3^k - 2^k!
""")

print(f"  {'n':>3s} {'2^n-1':>8s} {'3^n-2^n':>10s} {'4^n-3^n':>12s} {'5^n-4^n':>14s}")
for n in range(1, 16):
    print(f"  {n:3d} {2**n-1:8d} {3**n-2**n:10d} {4**n-3**n:12d} {5**n-4**n:14d}")

print("\n" + "=" * 70)
print("PART 3: THE ROLE OF 7 AND 8")
print("=" * 70)

print("""
  7 = (x+2)^n - (x+1)^n at n=3, x=0 = 2^3 - 1 = Mersenne prime!
  Also: 7 = Φ_3(2) = cyclotomic polynomial of order 3 at 2.
  Also: 7 = C(7,1) = first prime that enters at n=18 in the packing dictionary.
  Also: H=7 is IMPOSSIBLE for n=5 tournaments (the Forbidden 3 / S74).

  8 = 2^3 = (x+1)^n at n=3, x=1
  Also: 8 = first n where even-tournament structure appears
  Also: H_8 Hadamard matrix: creates the n=7 regular simplex in {0,1}^7

  HOW 7 AND 8 CHANGE THINGS in the (x+1)/(x+2) framework:

  At x=1 (tournament evaluation point):
    (x+1)^3 = 8 = 2^3    ← the simplex "fills" 8 slots
    (x+2)^3 = 27 = 3^3   ← the cuboid has 27 slots
    Complement: 27 - 8 = 19 = 3^3 - 2^3

  The number 19 appears as:
    - 3^3 - 2^3 = 19
    - Coefficient of alpha_3 in I(3)-I(2)
    - "The 19 corner pieces at level 3"

  7 + 8 = 15 = C(6,2) = number of edges in T_6
  7 * 8 = 56 = C(8,3) = number of 3-element subsets of 8

  THE (5,6) → (7,8) SHIFT:
  Tournament polynomial: z^2 - 5z + 6 = (z-2)(z-3)
  Roots: 2, 3. Sum=5, Product=6.
  The "next level": z^2 - 7z + 8... doesn't factor nicely.
  But z^2 - 7z + 12 = (z-3)(z-4) — shift up by 1!
  Sum=7, Product=12 = 3*4.
""")

# What polynomial has roots 2 and 3 with multiplicity reflecting the packing?
# z^2 - 5z + 6 = 0 at z=2,3
# The sequence p_k = 2^k + 3^k satisfies p_k = 5p_{k-1} - 6p_{k-2}

print("  The (5,6)-recurrence: p_k = 5p_{k-1} - 6p_{k-2}")
print(f"  {'k':>3s} {'2^k+3^k':>10s} {'5*prev-6*prev2':>16s} {'match':>6s}")
p = [2, 5]  # p_0 = 2^0+3^0 = 2, p_1 = 2+3 = 5
for k in range(2, 12):
    next_p = 5*p[-1] - 6*p[-2]
    actual = 2**k + 3**k
    print(f"  {k:3d} {actual:10d} {next_p:16d} {'✓' if actual == next_p else '✗':>6s}")
    p.append(next_p)

print("\n  NOW: what recurrence does 7^k + 8^k satisfy?")
print("  (z-7)(z-8) = z^2 - 15z + 56. Sum=15, Product=56.")
print(f"  {'k':>3s} {'7^k+8^k':>12s} {'15*prev-56*prev2':>18s} {'match':>6s}")
q = [2, 15]  # 7^0+8^0=2, 7^1+8^1=15
for k in range(2, 10):
    next_q = 15*q[-1] - 56*q[-2]
    actual = 7**k + 8**k
    print(f"  {k:3d} {actual:12d} {next_q:18d} {'✓' if actual == next_q else '✗':>6s}")
    q.append(next_q)

print("\n" + "=" * 70)
print("PART 4: 9=3² AND 10=5*2 IN THE POLYNOMIAL FRAMEWORK")
print("=" * 70)

print("""
  9 = 3² = (x+1)^2 at x=2 = simplex squared
  10 = 5*2 = (x+2)*(x+1) at x=1? No, (1+2)*(1+1) = 6.
       10 = C(5,2) = triangular number
       10 = (x+2)^n - (x+1)^n at... let's check.

  At x=1: 3^n - 2^n = 10 has no integer solution (3^2-2^2=5, 3^3-2^3=19).
  At x=0: 2^n - 1 = 10 → 2^n = 11, no solution.

  So 10 doesn't appear directly as a simplex-cuboid difference.

  But 9 = 3^2 DOES:
  (x+2)^n at x=1, n=2: 3^2 = 9. The 2-cuboid has 9 faces.
  (x+1)^n at x=2, n=2: 3^2 = 9. The 2-simplex at x=2 has 9 elements.

  The DOUBLE APPEARANCE of 9:
    9 = (1+2)^2 = 3^2 = cuboid^2
    9 = (2+1)^2 = 3^2 = simplex-at-x=2 squared

  This is because (x+2)^n at x=1 EQUALS (x+1)^n at x=2! Both give 3^n.
  It's the SAME evaluation: 3^n = I(Omega, 3) in tournament theory.

  THE BEAUTIFUL SYMMETRY:
    (x+1)^n|_{x=1} = 2^n = I(2) = "simplex"
    (x+1)^n|_{x=2} = 3^n = I(3) = "cuboid"
    (x+2)^n|_{x=1} = 3^n = I(3) = "cuboid"

  So I(3) can be computed EITHER as (x+1)^n at x=2 OR (x+2)^n at x=1.
  The simplex polynomial at a HIGHER evaluation point equals the
  cuboid polynomial at the STANDARD point!
""")

print("  GENERALIZATION: (x+r)^n at x=1 gives (r+1)^n")
print(f"  {'r':>3s}  {'(r+1)^1':>8s} {'(r+1)^2':>8s} {'(r+1)^3':>8s} {'(r+1)^4':>8s} {'name':>20s}")
for r in range(0, 8):
    vals = [(r+1)**n for n in range(1, 5)]
    name = {0: "simplex at x=0", 1: "simplex", 2: "cuboid",
            3: "4-hedroid", 4: "5-hedroid"}.get(r, f"(r={r})")
    print(f"  {r:3d}  {vals[0]:8d} {vals[1]:8d} {vals[2]:8d} {vals[3]:8d} {name:>20s}")

print("\n" + "=" * 70)
print("PART 5: THE 9=3² PACKING — SIMPLEX INSIDE SIMPLEX")
print("=" * 70)

print("""
  9 = 3² means: the 3-CYCLE can be SQUARED.
  In tournament theory: 3 disjoint 3-cycles cover all 9 vertices.
  This is alpha_3 first appearing.

  In the (x+1)/(x+2) framework:
  (x+1)^9 = ((x+1)^3)^3 = simplex^3
  (x+2)^9 = ((x+2)^3)^3 = cuboid^3

  The RATIO: ((x+2)/(x+1))^9 = ((x+2)/(x+1))^{3*3}
  At x=1: (3/2)^9 = 19683/512 ≈ 38.44

  But also: (3/2)^9 = ((3/2)^3)^3 = (27/8)^3 = 19683/512

  THE TRIPLE NESTING:
  (x+2)^9 = ((x+1)+1)^9 = sum C(9,k)(x+1)^k
  But ALSO:
  (x+2)^9 = (((x+1)+1)^3)^3

  Can we decompose the packing hierarchically?
  Level 1: cuboid^3 = simplex^3 + complement^3
  Level 2: each complement^3 = simplex^2 * complement + ...
  Level 3: ...

  This is the ITERATED PACKING of simplices inside cuboids.
  At 9=3², we get a DOUBLY NESTED structure.
  At 27=3³, a TRIPLY NESTED structure.
  At 3^k, a k-FOLD NESTED structure.
""")

# Verify iterated decomposition
print("  ITERATED DECOMPOSITION:")
print(f"  (x+2)^9 = ((x+2)^3)^3")
print(f"  At x=1: 3^9 = {3**9} = (3^3)^3 = {27**3}")
print()

# (y+1)^3 where y = (x+1)^3... no that's not right.
# (x+2)^9 = ((x+2)^3)^3 is just 3^9 at x=1.
# But the algebraic structure:
# (x+2)^3 = (x+1)^3 + 3(x+1)^2 + 3(x+1) + 1
# So (x+2)^9 = ((x+1)^3 + 3(x+1)^2 + 3(x+1) + 1)^3

print("  Let y = (x+1). Then (x+2) = y+1.")
print("  (x+2)^3 = (y+1)^3 = y^3 + 3y^2 + 3y + 1")
print("  (x+2)^9 = ((y+1)^3)^3 = (y^3 + 3y^2 + 3y + 1)^3")
print()
print("  But also: (x+2)^9 = (y+1)^9 = sum C(9,k) y^k")
print("  Both expressions are equal (of course!).")
print()
print("  The HIERARCHICAL view: ((y+1)^3)^3")
print("  First level: y^3 + (3y^2 + 3y + 1) = simplex + complement")
print("  Second level: (simplex + complement)^3")
print("    = simplex^3 + 3*simplex^2*complement + 3*simplex*complement^2 + complement^3")
print()
print("  At y=2 (x=1):")
print(f"    simplex = y^3 = 8")
print(f"    complement = 3y^2 + 3y + 1 = {3*4+3*2+1} = 19")
print(f"    simplex^3 = {8**3}")
print(f"    3*simplex^2*complement = {3*64*19}")
print(f"    3*simplex*complement^2 = {3*8*361}")
print(f"    complement^3 = {19**3}")
print(f"    Total = {8**3 + 3*64*19 + 3*8*361 + 19**3} = 3^9 = {3**9}")

print("\n" + "=" * 70)
print("PART 6: 10=5*2 — THE CROSS-PRODUCT")
print("=" * 70)

print("""
  10 = 5 * 2 = C(5,2)

  In tournament terms: n=10 is where 5-cycles can first pair up:
    Two disjoint 5-cycles covering all 10 vertices → alpha_2^{55}
    Also: (3,7) pair first appears → alpha_2^{37}

  In the polynomial framework:
  10 is NOT a power of a single base. It's a PRODUCT of distinct primes.

  (x+1)^{10} = ((x+1)^5)^2 = ((x+1)^2)^5
  (x+2)^{10} = ((x+2)^5)^2 = ((x+2)^2)^5

  Both factorizations are available: 10 = 5*2 = 2*5.

  The CROSS-PRODUCT structure:
  Let y = (x+1). Then (x+2)^10 = (y+1)^10.

  (y+1)^10 = sum C(10,k) y^k
  = y^10 + 10y^9 + 45y^8 + 120y^7 + 210y^6 + 252y^5 + 210y^4 + 120y^3 + 45y^2 + 10y + 1

  The LARGEST binomial coefficient is C(10,5) = 252 at the CENTER.
  This is EXACTLY where the "cuboid excess" peaks!

  At y=2: each term C(10,k)*2^k:
""")

n = 10
y = 2
print(f"  (y+1)^{n} at y={y}:")
total = 0
for k in range(n, -1, -1):
    term = comb(n, k) * y**k
    total += term
    pct = term / 3**n * 100
    bar = '#' * int(pct / 2)
    print(f"    C({n},{k})*{y}^{k} = {comb(n,k):>4d}*{y**k:>6d} = {term:>8d}  ({pct:5.1f}%) {bar}")
print(f"    Total = {total} = 3^{n} = {3**n}")

print(f"\n  Simplex part y^{n} = {y**n}")
print(f"  Complement = {3**n} - {y**n} = {3**n - y**n}")
print(f"  Complement/simplex = {(3**n - y**n)/(y**n):.4f}")

print("\n" + "=" * 70)
print("PART 7: THE CONNECTION — k-NACCI, SIMPLICES, AND TOURNAMENTS")
print("=" * 70)

print("""
  UNIFIED PICTURE:

  1. k-NACCI LIMIT = 2:
     The ratio of consecutive k-nacci numbers → 2.
     2 = (x+1)|_{x=1} = simplex evaluation point.
     Error: ~1/2^k = inverse simplex at dimension k.

  2. DOUBLED k-NACCI LIMIT = 3:
     The ratio of doubled k-nacci → 3.
     3 = (x+2)|_{x=1} = cuboid evaluation point.
     Error: ~2/3^k = 2 * inverse cuboid at dimension k.

  3. SIMPLEX-IN-CUBOID:
     (x+2)^n = sum C(n,k)(x+1)^k
     The cuboid is a binomial sum of simplex pieces.
     The complement has corner pieces of size 3^k - 2^k.

  4. TOURNAMENT INDEPENDENCE POLYNOMIAL:
     I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2 + ...
     I(2) = H (Hamiltonian paths = simplex evaluation)
     I(3) = "cuboid evaluation"
     I(3) - I(2) = sum alpha_k(3^k - 2^k) = corner piece count

  5. THE KEYS 2 AND 3:
     2 = smallest prime = simplex = k-nacci limit
     3 = smallest odd prime = cuboid = doubled k-nacci limit
     z^2 - 5z + 6 = (z-2)(z-3) = tournament polynomial
     5 = sum, 6 = product = 3! = 2*3

  6. HOW 7 AND 8 MODIFY:
     7 = Φ_3(2) = 2^3 - 1 = Mersenne prime = cuboid complement at n=3
     8 = 2^3 = simplex at n=3
     Together: (7,8) represents the "first nontrivial corner" —
     the first dimension where corner pieces have structure.

     7 + 8 = 15 = (3/2)^n at n=... no.
     15 = 2^4 - 1 = Mersenne number (not prime)
     15 = (x+2)^4 - (x+1)^4 at... no, that's 65.
     15 = C(6,2) = edges in 6-vertex tournament.
     15 = 3*5 = both odd primes multiplied.

  7. HOW 9=3² AND 10=5*2 FIT:
     9 = simplex evaluation squared = 3^2
     → Doubly-nested packing: ((y+1)^3)^3
     → Tournament: alpha_3 first appears (3 disjoint 3-cycles)
     → The SELF-SIMILAR level of the hierarchy

     10 = 5*2 = odd prime * even prime
     → Cross-product packing: ((y+1)^5)^2 or ((y+1)^2)^5
     → Tournament: alpha_2^{55} first appears (paired 5-cycles)
     → The MIXED level: odd and even interact

  THE FINAL SYNTHESIS:
  The simplex-cuboid framework IS the independence polynomial framework.
  Every tournament lives on a spectrum from:
    - TRANSITIVE (I=1, simplex=cuboid, no corners, pure order)
    - REGULAR (I maximal, simplex << cuboid, many corners, pure chaos)

  The "packing ratio" I(3)/I(2) is a single number measuring this:
    I(3)/I(2) = 1: transitive (simplex fills cuboid)
    I(3)/I(2) → (3/2)^{n/3}: maximally cyclic (empty simplex)

  The k-nacci limit (2) and doubled k-nacci limit (3) are the
  BOUNDARIES of this spectrum: the simplex wall and cuboid wall.
""")

# Compute the maximal I(3)/I(2) ratio for small n
print("  MAXIMAL I(3)/I(2) RATIOS vs (3/2)^{n/3}:")
print(f"  {'n':>3s} {'max I3/I2':>10s} {'(3/2)^(n/3)':>12s} {'ratio':>8s}")

# For n=5, we computed earlier: max ratio ≈ 1.4667
# Let's compute for n=5 and n=7
for n in [3, 5, 7]:
    num_edges = n*(n-1)//2
    max_ratio = 0
    for bits in range(min(2**num_edges, 2**15)):  # limit for large n
        # Build tournament
        A = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        # Count dc3
        dc3 = 0
        for i in range(n):
            for j in range(n):
                if i == j: continue
                for k in range(n):
                    if k == i or k == j: continue
                    if A[i][j] and A[j][k] and A[k][i]:
                        dc3 += 1
        dc3 //= 3

        # Count dc5 (for n >= 5)
        dc5 = 0
        if n >= 5:
            from itertools import combinations, permutations
            for combo in combinations(range(n), 5):
                for perm in permutations(combo):
                    if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                        dc5 += 1
            dc5 //= 5

        alpha1 = dc3 + dc5

        # Count HP
        dp = {}
        for v in range(n):
            dp[(1 << v, v)] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                if (mask, v) not in dp:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    if A[v][u]:
                        new_mask = mask | (1 << u)
                        dp[(new_mask, u)] = dp.get((new_mask, u), 0) + dp[(mask, v)]
        full_mask = (1 << n) - 1
        H = sum(dp.get((full_mask, v), 0) for v in range(n))

        alpha2 = (H - 1 - 2*alpha1) // 4 if H >= 1 + 2*alpha1 else 0
        I3 = 1 + 3*alpha1 + 9*alpha2
        ratio = I3 / H if H > 0 else 0
        if ratio > max_ratio:
            max_ratio = ratio

    theoretical = (3/2)**(n/3)
    print(f"  {n:3d} {max_ratio:10.4f} {theoretical:12.4f} {max_ratio/theoretical:8.4f}")

print("\nDone.")
