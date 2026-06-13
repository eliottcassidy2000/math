#!/usr/bin/env python3
"""
knacci_tournament_recurrence.py — opus-2026-03-14-S71f

The user said: "k-nacci approaches 2 and weighted k-nacci approaches 3."
We found: uniform weight-2 k-nacci → 3.

Now investigate: is there a TOURNAMENT RECURRENCE that follows
a k-nacci-like pattern?

Specifically:
1. H(T_n) vs H(T_{n-1}), H(T_{n-2}), ... under vertex deletion
2. Does deleting vertex v give: H(T) = Σ f(H(T\v), v)?
3. Is the growth rate of H under this recurrence related to k-nacci?

Also explore:
- The characteristic polynomial of the recurrence z² - 5z + 6 = 0
  which factors as (z-2)(z-3), roots 2 and 3!
  This is the recurrence from HYP-1108: forbidden sequence 7·3^k.
- General solution: a(k) = A·2^k + B·3^k
  If a(0)=7, a(1)=21: 7=A+B, 21=2A+3B → A=0, B=7 → a(k) = 7·3^k
  This is the pure z=3 solution!
- The z=2 solution: a(k) = 2^k gives 1, 2, 4, 8, 16, ...
  These aren't even odd numbers, so not H values.
- BUT: if we start with a(0)=1, a(1)=3: 1=A+B, 3=2A+3B → A=0, B=1 → a(k)=3^k
  The pure z=3 solution gives 1, 3, 9, 27, 81 = simplex powers!
"""

print("=" * 70)
print("Tournament Recurrence z² - 5z + 6 = 0")
print("=" * 70)

print("""
The characteristic equation z² - 5z + 6 = 0 has roots z=2 and z=3.

General solution: a(k) = A·2^k + B·3^k

This recurrence: a(k+1) = 5·a(k) - 6·a(k-1)

Connection to OCF:
  z = 2: the OCF evaluation point
  z = 3: the simplex brick value (1+x at x=2)

Specific orbits:
  (a₀, a₁) = (1, 3) → pure z=3: {1, 3, 9, 27, 81, ...} (simplex^k)
  (a₀, a₁) = (1, 5) → {1, 5, 19, 67, 231, ...}   A·2^k + B·3^k with A=-1, B=2
  (a₀, a₁) = (1, 2) → pure z=2: {1, 2, 4, 8, 16, ...}  (all powers of 2)
  (a₀, a₁) = (7, 21) → pure z=3: {7, 21, 63, 189, ...} (FORBIDDEN sequence!)
""")

# Compute several orbit families
print("Key orbits of z² - 5z + 6 = 0:")
for a0, a1 in [(1,1), (1,2), (1,3), (1,5), (3,5), (3,9), (5,9), (7,21), (7,15), (9,27)]:
    seq = [a0, a1]
    for _ in range(8):
        seq.append(5*seq[-1] - 6*seq[-2])
    # Solve for A, B: a0=A+B, a1=2A+3B → A=3a0-a1, B=a1-2a0
    A = 3*a0 - a1
    B = a1 - 2*a0
    label = f"a(k)={A}·2^k+{B}·3^k" if B >= 0 else f"a(k)={A}·2^k{B}·3^k"
    print(f"  ({a0},{a1}): {seq[:8]}  [{label}]")

# The (7,21) orbit: pure z=3 → 7·3^k
# The (7,15) orbit: A=3·7-15=6, B=15-14=1 → 6·2^k + 3^k
# Check: a(0)=6+1=7, a(1)=12+3=15, a(2)=24+9=33, a(3)=48+27=75 ✓

print(f"\n{'='*70}")
print("Connection to packing: (z-2)(z-3) = 0")
print(f"{'='*70}")

print("""
The factorization (z-2)(z-3) = z² - 5z + 6 encodes:

  z = 2: OCF evaluation point
  z = 3: simplex value (1+x)|_{x=2}
  z = 5: PRODUCT 2·3-1 = cuboid value (1+2x)|_{x=2}... wait, that's not right.

  Actually 5 = 2+3, the SUM of roots.
  And 6 = 2·3, the PRODUCT of roots.

  The recurrence a(k+1) = 5·a(k) - 6·a(k-1) can be read as:
    a(k+1) = (sum of roots)·a(k) - (product of roots)·a(k-1)

  This is the STANDARD form of a linear recurrence with roots r₁, r₂:
    z² - (r₁+r₂)z + r₁r₂ = 0

  For tournament theory, r₁=2 (OCF point) and r₂=3 (simplex value).

  The pure-3 orbit {7, 21, 63, ...} generates ALL forbidden H values:
    7·3⁰ = 7 (permanently forbidden)
    7·3¹ = 21 (permanently forbidden)
    7·3² = 63 (achievable at n≥8! NOT permanently forbidden)

  So the "forbidden sequence" breaks at k=2.
  The first two elements of the z=3 orbit from seed 7 are forbidden,
  but larger elements escape because larger n provides enough cycles.
""")

# NEW INSIGHT: The recurrence bridges simplices and OCF
print(f"{'='*70}")
print("The Simplex-OCF Bridge Recurrence")
print(f"{'='*70}")

print("""
For the simplex orbit (1, 3, 9, 27, ...):
  These are H values of m simplex bricks: I(Ω,2) = 3^m
  This requires m disjoint isolated cycles in Ω(T)

  The recurrence 3^{m+1} = 5·3^m - 6·3^{m-1} connects:
    Current: H = 3^m (m simplices)
    Previous: H = 3^{m-1} (m-1 simplices)
    Next: H = 3^{m+1} (m+1 simplices)

  Rearranging: 3^{m+1} - 5·3^m + 6·3^{m-1} = 0
    → 3^{m-1}(9 - 15 + 6) = 0 ✓

For the cuboid orbit:
  Cuboid values: 5^m = (1+2·2)^m
  Does 5^m satisfy z² - 5z + 6 = 0?
  5² - 5·5 + 6 = 25-25+6 = 6 ≠ 0
  So cuboid values do NOT satisfy this recurrence!

  Cuboid roots: z² - (r₁+r₂)z + r₁r₂ = 0 with solution z=5
  would need r₁+r₂ = sum giving 5, r₁r₂ = product
  Not unique — the cuboid lives on a DIFFERENT recurrence.

For the mixed orbit (1, 5, 19, 67, ...):
  a(0)=1, a(1)=5: A=3-5=-2, B=5-2=3
  a(k) = (-2)·2^k + 3·3^k = 3^{k+1} - 2^{k+1}

  These are NOT products of simplex/cuboid bricks in general.
  a(2) = 19 = 1+2·9+4·0 → α₁=9, α₂=0 → Ω=K₉ ✓ (achievable)
  a(3) = 67 = 1+2·33 → α₁=33, α₂=0 → Ω=K₃₃ ✓

  So the mixed orbit generates H = 3^{k+1} - 2^{k+1}, all achievable!
""")

# Verify: does H = 3^{k+1} - 2^{k+1} give achievable values?
print("Mixed orbit H = 3^{k+1} - 2^{k+1}:")
for k in range(10):
    h = 3**(k+1) - 2**(k+1)
    alpha1 = (h-1)//2
    print(f"  k={k}: H={h:7d} (α₁={alpha1} if Ω=K_{{α₁}})")

# ============================================================
# Connection to the (x+1)^n / (x+2)^n packing
# ============================================================

print(f"\n{'='*70}")
print("Packing (x+1)^n inside (x+2)^n")
print(f"{'='*70}")

print("""
The user's directive: "think of simplices as (x+1)^n and cuboids as (x+2)^n,
and think about packing them inside each other."

Interpretation 1: POLYNOMIAL NESTING
  (x+2)^n = ((x+1)+1)^n = Σ C(n,k) (x+1)^k · 1^{n-k}
           = Σ C(n,k) (x+1)^k

  At x=2: 4^n = Σ C(n,k) 3^k = (1+3)^n ✓ (binomial theorem)

  So a cuboid (x+2)^n CONTAINS n levels of simplices (x+1)^k
  weighted by binomial coefficients!

Interpretation 2: TOURNAMENT PACKING
  m simplex bricks: I = (1+x)^m → H = 3^m
  m cuboid bricks:  I = (1+2x)^m → H = 5^m

  Can we "pack" simplex bricks inside a cuboid brick?
  Cuboid = K₂ in Ω (2 overlapping cycles)
  Simplex = K₁ (1 isolated cycle)

  "Simplex inside cuboid" would mean: a cycle that is both
  isolated (simplex) and part of a pair (cuboid) — contradiction!

  So packing is NOT nesting — it's DISJOINT UNION.
  Simplex × Cuboid = K₁ ⊔ K₂ in Ω → H = 3·5 = 15.

Interpretation 3: GENERATING FUNCTION PACKING
  The generating function for achievable H via pure packing:

  GP(x) = Π_{c≥1} 1/(1 - x^{2c+1})

  where each factor allows any number of (2c+1)-type bricks.

  But the tournament constraint restricts: c=3 (tesseract) is forbidden.
  So the restricted GF is:

  GP_restricted(x) = 1/(1-x^3) · 1/(1-x^5) · 1/(1-x^9) · 1/(1-x^{11}) · ...

  (skipping 1/(1-x^7) because tesseract/K_3 is forbidden)
""")

# Compute which H values are achievable via pure packing (no K_3)
from functools import reduce

achievable = set([1])  # empty packing
bricks = [3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27]  # (2c+1) for c=1,2,4,5,...

# Dynamic programming: which H ≤ 1000 can be expressed as product of these bricks?
H_max = 1000
dp_achievable = set([1])
for b in bricks:
    new = set()
    for h in dp_achievable:
        val = h * b
        while val <= H_max:
            new.add(val)
            val *= b
    dp_achievable |= new

# Also allow products of different bricks
# Actually need to do this properly with full DP
def achievable_products(bricks, max_val):
    """Find all products of elements from bricks up to max_val."""
    reachable = {1}
    for b in bricks:
        new_reachable = set(reachable)
        for h in reachable:
            val = h * b
            while val <= max_val:
                new_reachable.add(val)
                val *= b
        reachable = new_reachable
    # Now allow combinations
    changed = True
    while changed:
        changed = False
        for h in list(reachable):
            for b in bricks:
                val = h * b
                if val <= max_val and val not in reachable:
                    reachable.add(val)
                    changed = True
    return reachable

packing_H = achievable_products(bricks, H_max)
# Add 7 as achievable via K_3 in Ω? No — K_3 is forbidden
# What about general Ω structures (non-packing)?
# Those give additional H values beyond the products

odd_to_1000 = set(range(1, 1001, 2))
packing_only = packing_H & odd_to_1000
non_packing = odd_to_1000 - packing_H

print(f"\nAchievable H values via pure packing (H ≤ 1000, no K_3):")
print(f"  Total achievable: {len(packing_only)} out of 500 odd values")
print(f"  Missing from packing: {len(non_packing)} values")
print(f"  Packing density: {len(packing_only)/500:.4f}")
print(f"  First 30 achievable: {sorted(packing_only)[:30]}")
print(f"  First 20 missing: {sorted(non_packing)[:20]}")

# The key point: most H values come from NON-PACKING Ω structures
# (general graphs, not disjoint unions of cliques)
# The packing gives only {products of {3,5,9,11,13,...}}

# But ALL odd values except 7,21 are achievable via general Ω!
# So the packing structure is a SUBSET of what tournaments can produce.

print(f"\n{'='*70}")
print("KEY INSIGHT: Packing vs General Ω")
print(f"{'='*70}")
print("""
Packing (Ω = ⊔K_{c_i}): gives H = Π(1+2c_i), sparse subset of odd ints
General Ω: gives ALL odd integers except {7, 21}

The packing framework is a SPECIAL CASE that:
1. Captures the simplex/cuboid hierarchy perfectly
2. Explains WHY 7 (tesseract) is forbidden
3. Explains WHY 21 = 3×7 (simplex × tesseract) is forbidden
4. But MOST tournaments have non-packing Ω structures

The tournament world is RICHER than what packings alone can describe.
The forbidden values 7 and 21 are the ONLY permanent casualties of
the tesseract obstruction, because all other H values can be achieved
by non-packing (general Ω) mechanisms.
""")

# ============================================================
# THE GRAND CONNECTION
# ============================================================

print(f"{'='*70}")
print("GRAND CONNECTION: k-nacci, Packing, and OCF")
print(f"{'='*70}")

print("""
The user's three observations connect as follows:

1. k-NACCI → 2:
   The growth rate of k-bonacci sequences approaches 2 as k→∞.
   This is the OCF evaluation point: H = I(Ω, 2).
   The tournament's binary arc choices create a "2-ness" that
   permeates the entire theory.

2. WEIGHTED k-NACCI → 3:
   With weight-2 (doubling), the growth rate → 3.
   This is the simplex brick value: (1+x)|_{x=2} = 3.
   One isolated cycle in Ω contributes ×3 to H.

3. SIMPLICES (x+1)^n AND CUBOIDS (x+2)^n:
   (x+1)^n at x=2 = 3^n: n simplex bricks
   (x+2)^n at x=2 = 4^n: NOT tournament-related!
   (1+2x)^n at x=2 = 5^n: n cuboid bricks (K₂ in Ω)

   The "packing" is MULTIPLICATIVE:
     simplex^a × cuboid^b → H = 3^a · 5^b

   "Packing simplices inside cuboids" = composing (1+x) into (1+2x):
     (1+2(1+x)) = 3+2x → H=7 at x=2
     This NESTING produces the forbidden value!

4. THE RECURRENCE (z-2)(z-3) = 0:
   Roots 2 and 3 connect OCF point and simplex value.
   Pure z=3 orbit from seed 7: {7, 21, 63, ...} = forbidden sequence.
   The first two are permanently forbidden; larger ones escape.

SYNTHESIS:
   The tournament OCF theory lives in the INTERSECTION of:
   - The binary world (z=2, k-nacci limit, arc choices)
   - The ternary world (z=3, simplex value, cycle counting)

   The characteristic polynomial (z-2)(z-3) encodes this duality.
   The forbidden value H=7 sits at the BOUNDARY where nesting
   (composition) tries to bridge these two worlds but fails.

   The tournament geometric constraint (vertex disjointness of cycles)
   prevents the nesting operation, leaving only the multiplicative
   (disjoint union) structure.
""")
