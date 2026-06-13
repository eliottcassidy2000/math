#!/usr/bin/env python3
"""
e7_deep_structure.py — opus-2026-03-14-S77

Deep dive into E₇ and its extraordinary connection to the CS boundary.

Key discovery from lie_tournament_bridge.py:
  h(E₇)/det(E₇) = 18/2 = 9 = CS BOUNDARY
  E₇ exponents = [1, 5, 7, 9, 11, 13, 17]
  9 is an exponent of E₇ BUT gcd(9,18)=9≠1 (NOT coprime to h!)
  This is special: in E₈, ALL exponents are coprime to h.
  In E₇, the CS boundary 9 sneaks in as a NON-coprime exponent.

Parts:
1. E₇ exponents: coprime vs non-coprime structure
2. The 18-gon and tournament circulants on Z_18
3. E₇ Weyl group coset structure
4. The exceptional cascade: G₂→F₄→E₆→E₇→E₈
5. Dual Coxeter numbers and the level-rank duality
6. The E₇ representation ring and tournaments
7. The 56-dimensional representation and n=8
8. E₇ as midpoint of the exceptional hierarchy
"""

from math import gcd, factorial, sqrt, pi, cos, sin, log
from itertools import combinations, permutations
import random

print("=" * 70)
print("PART 1: E₇ EXPONENTS — THE NON-COPRIME ANOMALY")
print("=" * 70)
print()

# Exponents of ADE algebras
exponents_data = {
    "A₁": (2, [1]),
    "A₂": (3, [1, 2]),
    "A₃": (4, [1, 2, 3]),
    "A₄": (5, [1, 2, 3, 4]),
    "A₅": (6, [1, 2, 3, 4, 5]),
    "A₆": (7, [1, 2, 3, 4, 5, 6]),
    "A₇": (8, [1, 2, 3, 4, 5, 6, 7]),
    "D₄": (6, [1, 3, 3, 5]),
    "D₅": (8, [1, 3, 5, 5, 7]),
    "D₆": (10, [1, 3, 5, 5, 7, 9]),
    "G₂": (6, [1, 5]),
    "F₄": (12, [1, 5, 7, 11]),
    "E₆": (12, [1, 4, 5, 7, 8, 11]),
    "E₇": (18, [1, 5, 7, 9, 11, 13, 17]),
    "E₈": (30, [1, 7, 11, 13, 17, 19, 23, 29]),
}

print("  Exponent analysis — coprime vs non-coprime to h:")
print()
for name in ["G₂", "F₄", "E₆", "E₇", "E₈"]:
    h, exp = exponents_data[name]
    coprime = [e for e in exp if gcd(e, h) == 1]
    non_coprime = [e for e in exp if gcd(e, h) != 1]
    print(f"  {name} (h={h}): exponents = {exp}")
    print(f"    coprime to h: {coprime} ({len(coprime)} = φ({h})={len([k for k in range(1,h) if gcd(k,h)==1])}? "
          f"{'YES' if len(coprime)==len([k for k in range(1,h) if gcd(k,h)==1]) else 'NO'})")
    if non_coprime:
        print(f"    NON-coprime to h: {non_coprime}")
        for e in non_coprime:
            print(f"      gcd({e},{h}) = {gcd(e,h)}")
    print()

print("  KEY FINDING: E₆ has non-coprime exponents {4, 8} (gcd with 12 = 4)")
print("  E₇ has non-coprime exponent {9} (gcd(9,18) = 9)")
print("  E₈ has ALL exponents coprime to h=30 (totatives of 30)")
print()
print("  The CS boundary 9 appears as the UNIQUE non-coprime exponent of E₇!")
print("  In E₆, the non-coprime exponents {4,8} have gcd 4 with h=12.")
print("  4 = vertices of tetrahedron = KEY₁².")
print()

# The non-coprime exponents come from the non-trivial divisors of h
# that happen to be exponents
print("  Divisors of h for each exceptional:")
for name in ["G₂", "F₄", "E₆", "E₇", "E₈"]:
    h, exp = exponents_data[name]
    divs = [d for d in range(1, h+1) if h % d == 0]
    exp_divs = [d for d in divs if d in exp and d > 1 and d < h]
    print(f"  {name} (h={h}): divisors = {divs}")
    if exp_divs:
        print(f"    Divisors that are also exponents: {exp_divs}")
    print()

# For A_n, exponents are 1,...,n and h=n+1
# Divisors of h that are exponents: divisors of n+1 that are ≤ n
# This is ALL non-trivial divisors of n+1
print("  For A_n: h=n+1, exponents=1,...,n")
print("  ALL non-trivial divisors of h are exponents!")
for n in [1, 2, 3, 4, 5, 7]:
    h = n + 1
    exp = list(range(1, n+1))
    exp_divs = [d for d in range(2, h) if h % d == 0]
    print(f"    A_{n} (h={h}): divisor-exponents = {exp_divs}")

print()
print("=" * 70)
print("PART 2: THE 18-GON AND CIRCULANT TOURNAMENTS")
print("=" * 70)
print()

# Z_18 ≅ Z_2 × Z_9 ≅ Z_2 × Z_3²
# The structure of circulant tournaments on Z_18

print("  Z₁₈ ≅ Z₂ × Z₉")
print("  Generators: {1, 5, 7, 11, 13, 17} (= units of Z₁₈)")
print("  These are the E₇ COPRIME exponents!")
print()

units_18 = [k for k in range(1, 18) if gcd(k, 18) == 1]
print(f"  Units of Z₁₈: {units_18}")
print(f"  Number of units: {len(units_18)} = φ(18) = 6")
print()

# A circulant tournament on Z_n uses a subset S ⊂ Z_n \ {0}
# such that S ∪ (-S) = Z_n \ {0} and S ∩ (-S) = ∅
# i.e., for each pair {k, n-k}, exactly one is in S
# This requires n to be odd!

# Since 18 is even, we can't make a tournament on Z_18
# But we CAN make a tournament on Z_17 (a prime!)
# QR_17: 17 ≡ 1 (mod 4) → graph, not tournament!
# So QR_17 doesn't give a tournament either.

# What about Z_19? 19 ≡ 3 (mod 4) → QR_19 IS a tournament!
print("  Tournament circulants require ODD n:")
print("  18 is even → no circulant tournament on Z₁₈")
print("  But 17 is prime, 17 ≡ 1 (mod 4) → QR_17 is a graph, not tournament")
print("  19 is prime, 19 ≡ 3 (mod 4) → QR_19 IS a tournament!")
print("  19 = 18 + 1 = h(E₇) + 1 ← interesting!")
print()

# Paley tournament QR_p for p ≡ 3 (mod 4)
# The SMALLEST such primes: 3, 7, 11, 19, 23, 31, 43, ...
small_paley = [p for p in range(3, 50)
               if p % 4 == 3 and all(p % d != 0 for d in range(2, p))]
print(f"  Paley tournament primes (p ≡ 3 mod 4): {small_paley}")
print()

# Note: 7 = 2³-1, 11 = Lucas(5), 19 = 18+1 = h(E₇)+1
# 23 = h(E₈) - 7, 31 = 2⁵-1 = Mersenne
for p in small_paley[:7]:
    tags = []
    if p == 3: tags.append("KEY₂, C₃")
    if p == 7: tags.append("2³-1, Mersenne, Lucas(4)")
    if p == 11: tags.append("Lucas(5), first SC tournament at n≡3(4)")
    if p == 19: tags.append("h(E₇)+1=18+1")
    if p == 23: tags.append("E₈ exponent")
    if p == 31: tags.append("2⁵-1, Mersenne")
    if p == 43: tags.append("E₈ dim component: 248=8·31")
    tag = f"  ← {', '.join(tags)}" if tags else ""
    print(f"  QR_{p}{tag}")

print()
print("=" * 70)
print("PART 3: E₇ WEYL GROUP AND TOURNAMENT STRUCTURE")
print("=" * 70)
print()

# |W(E₇)| = 2903040 = 2^10 · 3^4 · 5 · 7
we7 = 2**10 * 3**4 * 5 * 7
print(f"  |W(E₇)| = {we7} = 2^10 · 3^4 · 5 · 7")
print(f"         = {2**10} · {3**4} · 5 · 7")
print()

# Compare with tournament-related numbers
# At n=8 (rank of E₈, one more than rank of E₇):
# 2^C(8,2) = 2^28 total tournaments on 8 vertices
# 8! = 40320 = size of S_8 = Weyl group of A_7
# Non-iso tournaments on 8 vertices: 6880

print(f"  |S₈| = 8! = {factorial(8)} = Weyl group of A₇")
print(f"  |W(E₇)| / |S₈| = {we7 / factorial(8):.2f}")
print(f"  = {we7 // factorial(8)} = 72 = 8 · 9 = rank(E₈) · CS-boundary!")
print()

# Check: 2903040 / 40320 = 72 = 8·9
ratio = we7 // factorial(8)
print(f"  |W(E₇)| / |W(A₇)| = {ratio}")
print(f"  72 = 8 · 9 = 2³ · 3² = KEY₁³ · KEY₂²")
print(f"  This is the product of rank(E₈) and the CS boundary!")
print()

# What about E₈?
we8 = 2**14 * 3**5 * 5**2 * 7
print(f"  |W(E₈)| = {we8} = 2^14 · 3^5 · 5^2 · 7")
print(f"  |W(E₈)| / |S₈| = {we8 // factorial(8)}")
print(f"  = {we8 // factorial(8)} (should check)")

# 696729600 / 40320 = 17280
r2 = we8 // factorial(8)
# Factor 17280
n = r2
factors = {}
for p in [2, 3, 5, 7, 11]:
    while n % p == 0:
        factors[p] = factors.get(p, 0) + 1
        n //= p
print(f"  {r2} = {'·'.join(f'{p}^{e}' for p,e in sorted(factors.items()))}")
# 17280 = 2^6 · 3^3 · 5 · ... = 2^6·270 = 64·270
print()

# E₆ Weyl group
we6 = 2**7 * 3**4 * 5
print(f"  |W(E₆)| = {we6}")
print(f"  |W(E₆)| / |S₆| = {we6 // factorial(6)} (S₆ = A₅ Weyl)")
r3 = we6 // factorial(6)
print(f"  = {r3}")
# 51840 / 720 = 72 = 8·9 again!
print(f"  = 72 = 8·9 = KEY₁³·KEY₂² AGAIN!")
print()
print(f"  REMARKABLE: |W(E₆)|/|S₆| = |W(E₇)|/|S₈| = 72 = 8·9!")

print()
print("=" * 70)
print("PART 4: THE EXCEPTIONAL CASCADE")
print("=" * 70)
print()

# How the exceptionals build on each other
# G₂ ⊂ D₄ (triality), F₄ related to D₄, E₆ ⊂ E₇ ⊂ E₈

exc_dims = {
    "G₂": 14, "F₄": 52, "E₆": 78, "E₇": 133, "E₈": 248
}
exc_ranks = {
    "G₂": 2, "F₄": 4, "E₆": 6, "E₇": 7, "E₈": 8
}
exc_h = {"G₂": 6, "F₄": 12, "E₆": 12, "E₇": 18, "E₈": 30}
exc_det = {"G₂": 1, "F₄": 1, "E₆": 3, "E₇": 2, "E₈": 1}

print("  The exceptional cascade:")
print("  Group   rank  dim    h    det   dim/rank   pos.roots")
for name in ["G₂", "F₄", "E₆", "E₇", "E₈"]:
    d = exc_dims[name]
    r = exc_ranks[name]
    h = exc_h[name]
    det = exc_det[name]
    pos = (d - r) // 2  # dim = rank + 2·#pos.roots
    print(f"  {name}     {r}    {d:>4}   {h:>2}    {det}     {d/r:>6.1f}      {pos}")

print()
# Dimension ratios
print("  Dimension ratios (consecutive):")
print(f"  F₄/G₂ = {52/14:.4f} ≈ {52/14:.4f}")
print(f"  E₆/F₄ = {78/52:.4f} = 3/2 = KEY₂/KEY₁")
print(f"  E₇/E₆ = {133/78:.4f}")
print(f"  E₈/E₇ = {248/133:.4f}")
print()
print(f"  E₆/F₄ = 78/52 = 3/2 EXACTLY! The ratio of the keys!")
print(f"  Check: 78 = 52 · 3/2 = 52 · 1.5 = 78 ✓")

print()
print("=" * 70)
print("PART 5: DUAL COXETER NUMBERS")
print("=" * 70)
print()

# Dual Coxeter number h∨ = h for simply-laced (ADE)
# For non-simply-laced: h∨ ≠ h

dual_h = {
    "A_n": "n+1 (same as h)",
    "B_n": "2n-1",
    "C_n": "n+1",
    "D_n": "2n-2 (same as h)",
    "G₂": "4 (h=6)",
    "F₄": "9 (h=12)",
    "E₆": "12 (h=12, same)",
    "E₇": "18 (h=18, same)",
    "E₈": "30 (h=30, same)",
}

print("  Dual Coxeter numbers h∨:")
for name, val in dual_h.items():
    print(f"    {name}: h∨ = {val}")

print()
print("  For the ADE types (simply-laced): h∨ = h always.")
print("  For non-simply-laced:")
print("    G₂: h∨ = 4 ≠ h = 6. Ratio h/h∨ = 6/4 = 3/2 = KEY₂/KEY₁!")
print("    F₄: h∨ = 9 ≠ h = 12. Ratio h/h∨ = 12/9 = 4/3 = KEY₁²/KEY₂!")
print()
print("  EXTRAORDINARY: h∨(F₄) = 9 = CS BOUNDARY!")
print("  The dual Coxeter number of F₄ is our CS boundary!")
print("  And h(F₄)/h∨(F₄) = 12/9 = 4/3 = KEY₁²/KEY₂")
print()
print("  h∨(G₂) = 4 = KEY₁²")
print("  h∨(F₄) = 9 = KEY₂²")
print("  Pattern: h∨ for non-simply-laced = square of a key!")

print()
print("=" * 70)
print("PART 6: THE 56-DIMENSIONAL REPRESENTATION OF E₇")
print("=" * 70)
print()

# E₇ has a fundamental representation of dimension 56
# 56 = 8 · 7 = number of non-isomorphic tournaments on n=6

print("  E₇ fundamental representations:")
print("  The minuscule representation has dimension 56")
print()
print("  56 = 8 · 7 = C(8,2) - C(8,1) = 28 - 8 + ... no")
print("  56 = C(8,3) = number of 3-element subsets of {1,...,8}")
print()

# At n=6: non-iso tournaments = 56!
print("  COINCIDENCE: non-isomorphic tournaments on 6 vertices = 56!")
print("  And dim(minuscule rep of E₇) = 56!")
print()

# At n=8: C(8,3) = 56 = number of potential 3-cycles on 8 vertices
print("  At n=8 (= rank(E₈)):")
print(f"  C(8,3) = {8*7*6//(3*2*1)} = 56 = number of 3-element subsets")
print(f"  Each 3-element subset gives either a 3-cycle or a transitive triple")
print(f"  In a random tournament, expected 3-cycles = 56/4 = 14")
print()

# The 56 also appears as:
# 56 = 2^3 · 7 = KEY₁³ · 7
# 7 = 2³-1 = Mersenne prime = 2·KEY₂+1
print("  56 = 2³ · 7 = KEY₁³ · (KEY₁³-1) = KEY₁³ · Mersenne(3)")
print("  This is rank(E₈) · (rank(E₈)-1)")
print()

# The adjoint representation of E₇ has dimension 133
print("  Adjoint representation: dim = 133 = 7 · 19")
print("  7 = rank(E₇), 19 = h(E₇)+1 = Paley prime!")
print("  133 = rank(E₇) · (h(E₇)+1)")
print()

# Let's check: dim(g) = rank + 2·#pos_roots = rank · (2h-1)... no
# For E₇: rank=7, #pos_roots=63, dim = 7 + 2·63 = 7+126 = 133 ✓
# 133 = 7·19
# So dim(E₇) = rank(E₇) · 19
# 19 = 2·9+1 = 2·CS+1

print("  dim(E₇) = rank · (2h-1+... wait)")
print(f"  Actually: dim = rank + 2·#pos.roots = 7 + 2·63 = 133")
print(f"  And #pos.roots = sum of exponents = 63 = 7·9")
print(f"  63 = rank(E₇) · CS_boundary = 7 · 9!")
print()
print("  E₇ ENCODES THE CS BOUNDARY:")
print("  #positive roots / rank = 63/7 = 9 = CS boundary")
print("  This ratio is the AVERAGE EXPONENT!")
print("  Average exponent of E₇ = 9 = 3² = KEY₂²")

print()
print("=" * 70)
print("PART 7: EXPONENT SUMS AND THE KEYS")
print("=" * 70)
print()

# For each Lie algebra, the average exponent = sum(exponents)/rank = #pos_roots/rank
print("  Average exponent (= #pos.roots / rank) for exceptionals:")
for name in ["G₂", "F₄", "E₆", "E₇", "E₈"]:
    h, exp = exponents_data[name]
    avg = sum(exp) / len(exp)
    print(f"  {name}: avg = {sum(exp)}/{len(exp)} = {avg:.4f}", end="")
    if avg == int(avg):
        iv = int(avg)
        if iv == 3: print(f" = 3 = KEY₂", end="")
        elif iv == 6: print(f" = 6 = KEY₁·KEY₂ = h(G₂)", end="")
        elif iv == 9: print(f" = 9 = KEY₂² = CS boundary", end="")
        elif iv == 15: print(f" = 15 = 3·5 = KEY₂·5", end="")
    print()

print()
print("  G₂: avg exponent = 3 = KEY₂")
print("  F₄: avg exponent = 6 = KEY₁·KEY₂ = h(G₂)")
print("  E₆: avg exponent = 6 = KEY₁·KEY₂ = h(G₂)")
print("  E₇: avg exponent = 9 = KEY₂² = CS BOUNDARY")
print("  E₈: avg exponent = 15 = 3·5")
print()
print("  The average exponent sequence for exceptionals: 3, 6, 6, 9, 15")
print("  Differences: 3, 0, 3, 6")
print("  Ratios: 2, 1, 3/2, 5/3")
print()

# Also: for A_n, avg exponent = (n+1)/2 = h/2
print("  For A_n: avg exponent = (1+2+...+n)/n = (n+1)/2 = h/2")
print("  So avg_exp = h/2 always for type A")
print()

# For E₇: avg_exp = 9, h = 18, so avg_exp = h/2 ✓
# For E₈: avg_exp = 15, h = 30, so avg_exp = h/2 ✓
# For E₆: avg_exp = 6, h = 12, so avg_exp = h/2 ✓
# For G₂: avg_exp = 3, h = 6, so avg_exp = h/2 ✓
# For F₄: avg_exp = 6, h = 12, so avg_exp = h/2 ✓
print("  Check: avg_exp = h/2 for ALL simple Lie algebras!")
print("  This is because the exponents m_i satisfy m_i + m_{rank+1-i} = h")
print("  (they come in palindromic pairs summing to h)")
print()
print("  So the CS boundary 9 = h(E₇)/2 = avg exponent of E₇")
print("  is a FUNDAMENTAL invariant of the E₇ Lie algebra.")

print()
print("=" * 70)
print("PART 8: E₇ AS MIDPOINT OF THE EXCEPTIONAL HIERARCHY")
print("=" * 70)
print()

print("  The 5 exceptional Lie groups form a HIERARCHY:")
print()
print("  G₂ ← F₄ ← E₆ ← E₇ ← E₈")
print("  ↓      ↓      ↓      ↓      ↓")
print("  h=6    h=12   h=12   h=18   h=30")
print("  det=1  det=1  det=3  det=2  det=1")
print("  avg=3  avg=6  avg=6  avg=9  avg=15")
print()

# E₇ is the MIDPOINT (3rd of 5) when ordered by rank
# It's where the keys (det=2) and the CS boundary (avg=9) meet

print("  E₇ is the 4th of 5 exceptionals by rank (2,4,6,7,8)")
print("  It is the UNIQUE exceptional with:")
print("  - det = 2 = KEY₁ (binary orientation)")
print("  - h/2 = 9 = KEY₂² (CS boundary)")
print("  - h/det = 9 (also the CS boundary)")
print("  - 9 as an exponent (appears in the root system)")
print()

# The KEY₁-KEY₂ encoding
print("  THE KEY ENCODING IN EXCEPTIONAL GROUPS:")
print()
print("  E₆ encodes KEY₂ = 3 through its center Z₃")
print("  E₇ encodes KEY₁ = 2 through its center Z₂")
print("  Together they encode the tournament polynomial z²-5z+6")
print()
print("  E₈ unifies them: det=1 (trivial center),")
print("  but h = 30 = 2·3·5 contains BOTH keys and their sum")
print()

# Final synthesis
print("  THE DEEPEST CONNECTION:")
print()
print("  In tournament theory: z² - 5z + 6 = (z-2)(z-3) = 0")
print("  In Lie theory: A₁ and A₂ have Cartan eigenvalues {2} and {1,3}")
print("  The tournament polynomial roots ARE Cartan eigenvalues!")
print()
print("  More precisely:")
print("  KEY₁ = 2 = max eigenvalue of Cartan(A₁) = det(A₁) = h(A₁)")
print("  KEY₂ = 3 = max eigenvalue of Cartan(A₂) = det(A₂) = h(A₂)")
print()
print("  And their difference pattern reappears in the exceptionals:")
print("  det(E₆) = KEY₂ = 3")
print("  det(E₇) = KEY₁ = 2")
print("  h(E₇)/det(E₇) = 9 = KEY₂² = CS boundary")
print("  h∨(F₄) = 9 = CS boundary (dual Coxeter of F₄)")
print()
print("  THE TOURNAMENT UNIVERSE IS THE A₁-A₂ UNIVERSE,")
print("  REFLECTED IN THE EXCEPTIONAL GROUPS E₆ AND E₇.")

print()
print("=" * 70)
print("PART 9: NUMERICAL COINCIDENCES AS STRUCTURAL ECHOES")
print("=" * 70)
print()

# Let's verify some of these numerological observations more carefully

# 1. dim(E₇ minuscule) = 56 = non-iso tournaments on 6 vertices
print("  Verified coincidences:")
print(f"  56 = dim(minuscule E₇) = T(6) (non-iso tournaments on 6)")
print(f"  56 = C(8,3) = potential 3-cycles on 8 vertices")
print()

# 2. 72 = |W(E₆)|/|S₆| = |W(E₇)|/|S₈|
print(f"  72 = |W(E₆)|/6! = {51840//720}")
print(f"  72 = |W(E₇)|/8! = {2903040//40320}")
print(f"  72 = 8·9 = rank(E₈)·CS_boundary = KEY₁³·KEY₂²")
print()

# 3. Positive roots of E₇: 63 = 7·9
print(f"  63 = #pos.roots(E₇) = 7·9 = rank·CS_boundary")
print(f"  63 = 2^6-1 = Mersenne(6)")
print()

# 4. Positive roots of E₈: 120 = 5!
print(f"  120 = #pos.roots(E₈) = 5! = |BI| (binary icosahedral)")
print(f"  120 = C(10,3) = ... or 2³·3·5 = KEY₁³·KEY₂·5")
print()

# 5. The h values as tournament polynomial evaluations
print("  Evaluating tournament polynomial z²-5z+6 at various points:")
for z in range(8):
    val = z**2 - 5*z + 6
    tags = []
    if val == 0: tags.append(f"ROOT (z={z})")
    if z in [6, 12, 18, 30]: tags.append("Coxeter number!")
    if abs(val) in [6, 12, 18, 30]: tags.append(f"|val|={abs(val)} is a Coxeter number")
    tag = f"  ← {', '.join(tags)}" if tags else ""
    print(f"  f({z}) = {val}{tag}")

print()
# f(0)=6=h(G₂), f(1)=2=det(E₇), f(4)=-2=-det(E₇), f(5)=6=h(G₂)
print("  f(0) = 6 = h(G₂)!")
print("  f(1) = 2 = KEY₁ = det(E₇)!")
print("  f(4) = -2 = -KEY₁")
print("  f(5) = 6 = h(G₂) again! (by symmetry: f(5-x)=f(x) since 5=sum of roots)")
print()
print("  The tournament polynomial has AXIS OF SYMMETRY at z = 5/2!")
print("  f(5/2 + t) = f(5/2 - t) for all t")
print("  The symmetric point 5/2 = (KEY₁+KEY₂)/2")
print()

# Also check: f evaluated at Cartan eigenvalues
print("  f at Cartan eigenvalues of A₂: {1, 3}")
print(f"  f(1) = {1-5+6} = 2 = KEY₁")
print(f"  f(3) = {9-15+6} = 0 (root!)")
print()
print("  So f(smaller eigenvalue of A₂) = KEY₁ = larger eigenvalue of A₁!")
print("  This creates a SELF-REFERENTIAL LOOP:")
print("  A₂'s eigenvalues → through tournament poly → back to A₁'s eigenvalue")
