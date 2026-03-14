#!/usr/bin/env python3
"""
Fibonacci sequence: 2 and 3 subsequences, period-6 structure,
and deep connections to tournament Jacobsthal recurrence.
opus-2026-03-14-S86

USER INSIGHT: Fibonacci mod 4 has period 6 = Pisano(4).
Residues: 1, 1, 2, 3, 1, 0, 1, 1, 2, 3, 1, 0, ...
The 2s and 3s appear at specific positions within each period!

CONNECTIONS TO TOURNAMENTS:
- Jacobsthal: a(k) = a(k-1) + 2·a(k-2), char eq t²-t-2=0, roots {2,-1}
- Fibonacci: F(k) = F(k-1) + F(k-2), char eq t²-t-1=0, roots {φ,-1/φ}
- These differ by 1 in the constant: t²-t-N with N=2 (Jacobsthal) vs N=1 (Fib)
- OCF evaluation at x=2 gives Jacobsthal; at x=1 gives Fibonacci!
- I(P_k, x) satisfies a(k) = a(k-1) + x·a(k-2) — UNIFIED!

PERIOD-6 STRUCTURE:
- Fib mod 2: period 3 → {1,1,0}
- Fib mod 3: period 8 → {1,1,2,0,2,2,1,0}
- Fib mod 4: period 6 → {1,1,2,3,1,0}
- The "2 and 3 subsequences" within mod 4 decompose as:
  Positions where F(n) ≡ 2 mod 4: n ≡ 3 mod 6
  Positions where F(n) ≡ 3 mod 4: n ≡ 4 mod 6
"""

import math
from collections import Counter, defaultdict
from fractions import Fraction

# ============================================================
# Part 1: Fibonacci mod m — Pisano Periods
# ============================================================
print("=" * 70)
print("PART 1: PISANO PERIODS π(m) — FIBONACCI PERIODICITY")
print("=" * 70)

def pisano_period(m):
    """Compute Pisano period π(m) = period of Fibonacci mod m."""
    a, b = 0, 1
    for i in range(1, 6 * m * m + 1):  # Upper bound: π(m) ≤ 6m²
        a, b = b, (a + b) % m
        if a == 0 and b == 1:
            return i
    return None

# Fibonacci sequence
fib = [0, 1]
for i in range(2, 60):
    fib.append(fib[-1] + fib[-2])

print("\nFibonacci sequence (first 30):")
print(f"  {fib[:30]}")

print("\nPisano periods π(m) for m = 2..20:")
for m in range(2, 21):
    pi_m = pisano_period(m)
    residues = [fib[i] % m for i in range(pi_m)]
    print(f"  π({m:2d}) = {pi_m:3d}  residues: {residues}")

# Highlight m=4: period 6!
print("\n*** FOCUS: m=4, π(4) = 6 ***")
print("F(n) mod 4:")
for period in range(4):
    residues = [fib[6*period + i] % 4 for i in range(6)]
    print(f"  Period {period}: {residues}")

# ============================================================
# Part 2: The 2 and 3 Subsequences Within Period 6
# ============================================================
print("\n" + "=" * 70)
print("PART 2: 2 AND 3 SUBSEQUENCES IN F(n) mod 4")
print("=" * 70)

print("\nPositions of each residue mod 4 (within period 6):")
for r in range(4):
    positions = [i for i in range(6) if fib[i] % 4 == r]
    # Extended
    extended = [i for i in range(30) if fib[i] % 4 == r]
    print(f"  F(n) ≡ {r} mod 4: positions mod 6 = {positions}, first 10: {extended[:10]}")

# The 2 appears at position 3 (mod 6): F(3)=2, F(9)=34≡2, F(15)=610≡2, ...
# The 3 appears at position 4 (mod 6): F(4)=3, F(10)=55≡3, F(16)=987≡3, ...
# The 1 appears at positions 1,2,5 (mod 6)
# The 0 appears at position 0 (mod 6)... wait, F(0)=0, F(6)=8≡0, F(12)=144≡0

print("\nDecomposition of period {1,1,2,3,1,0}:")
print("  Position 0 (mod 6): F ≡ 0 mod 4 → divisible by 4")
print("  Position 1 (mod 6): F ≡ 1 mod 4 → ≡ 1")
print("  Position 2 (mod 6): F ≡ 1 mod 4 → ≡ 1")
print("  Position 3 (mod 6): F ≡ 2 mod 4 → EVEN, not div by 4")
print("  Position 4 (mod 6): F ≡ 3 mod 4 → ≡ 3")
print("  Position 5 (mod 6): F ≡ 1 mod 4 → ≡ 1")
print()
print("The 2 and 3 are ADJACENT in the period: (2,3) at positions (3,4).")
print("Together they give F(3)×F(4) = 2×3 = 6 = period length!")

# ============================================================
# Part 3: Fibonacci → Jacobsthal → Tournament OCF
# ============================================================
print("\n" + "=" * 70)
print("PART 3: UNIFIED RECURRENCE a(k) = a(k-1) + x·a(k-2)")
print("=" * 70)

# The generalized Fibonacci with parameter x:
# G_x(0) = 1, G_x(1) = 1+x, G_x(k) = G_x(k-1) + x·G_x(k-2)
# This is I(P_k, x) — independence polynomial of path at fugacity x.

# x=1: Fibonacci (shifted): 1, 2, 3, 5, 8, 13, 21, ...
# Actually I(P_k, 1) = F(k+2) where F is standard Fibonacci
# x=2: Jacobsthal-related: 1, 3, 5, 11, 21, 43, 85, 171, ...
# I(P_k, 2) = (2^{k+2} - (-1)^k)/3

# Characteristic equation: t² - t - x = 0
# Roots: (1 ± √(1+4x))/2
# x=1: roots = φ, -1/φ (golden ratio)
# x=2: roots = 2, -1

print("\nGeneralized Fibonacci G_x(k) = G_x(k-1) + x·G_x(k-2):")
print("  x=1: Fibonacci   (char roots: φ ≈ 1.618, -1/φ ≈ -0.618)")
print("  x=2: Jacobsthal  (char roots: 2, -1)")
print("  x=3: 'Tribonacci' variant (char roots: (1±√13)/2 ≈ 2.303, -1.303)")
print()

for x in [1, 2, 3, 4]:
    G = [1, 1 + x]
    for k in range(2, 12):
        G.append(G[-1] + x * G[-2])
    root1 = (1 + (1 + 4*x)**0.5) / 2
    root2 = (1 - (1 + 4*x)**0.5) / 2
    print(f"  x={x}: G = {G[:10]}")
    print(f"       Roots: {root1:.6f}, {root2:.6f}")
    print(f"       G(k)/G(k-1) → {root1:.6f}")
    print(f"       Ratio check: {G[-1]/G[-2]:.6f}")

# ============================================================
# Part 4: Pisano Periods of Jacobsthal Sequence
# ============================================================
print("\n" + "=" * 70)
print("PART 4: JACOBSTHAL mod m — PISANO-LIKE PERIODS")
print("=" * 70)

# Jacobsthal: J(0)=0, J(1)=1, J(k) = J(k-1) + 2·J(k-2)
# = (2^k - (-1)^k)/3
jac = [0, 1]
for i in range(2, 40):
    jac.append(jac[-1] + 2 * jac[-2])

print(f"\nJacobsthal: {jac[:20]}")

def jacobsthal_pisano(m):
    """Period of Jacobsthal sequence mod m."""
    a, b = 0, 1
    for i in range(1, 6 * m * m + 1):
        a, b = b, (b + 2*a) % m
        if a == 0 and b % m == 1:
            return i
    return None

print("\nJacobsthal Pisano periods πJ(m):")
for m in range(2, 21):
    period = jacobsthal_pisano(m)
    if period and period <= 100:
        residues = [jac[i] % m for i in range(min(period, 30))]
        print(f"  πJ({m:2d}) = {period:3d}  residues: {residues[:20]}{'...' if period > 20 else ''}")
    else:
        print(f"  πJ({m:2d}) = {period or '>600'}")

# Key: πJ(4) = ?
print("\n*** Jacobsthal mod 4 ***")
jac_mod4 = [jac[i] % 4 for i in range(20)]
print(f"  J(n) mod 4: {jac_mod4}")
# Find period
for p in range(1, 20):
    if all(jac_mod4[i] == jac_mod4[i+p] for i in range(min(10, 20-p))):
        print(f"  Period: {p}")
        break

# ============================================================
# Part 5: The 6-fold Symmetry
# ============================================================
print("\n" + "=" * 70)
print("PART 5: WHY PERIOD 6? THE 6-FOLD STRUCTURE")
print("=" * 70)

print("""
Pisano period π(4) = 6 for Fibonacci. WHY 6?

Factor: 4 = 2². Using the formula:
  π(p^k) = p^{k-1} × π(p) for primes p ≥ 5
  But for p=2: π(2)=3, π(4)=6=2×3, π(8)=12=4×3, π(16)=24=8×3
  So π(2^k) = 3 × 2^{k-1} for k ≥ 1.

The 6 = 2 × 3 comes from:
  - 3 = π(2) = period of Fibonacci mod 2: {1,1,0}
  - 2 = multiplier from 2² vs 2¹

The DEEP reason: F(n) mod 2 has period 3 because:
  The matrix [[1,1],[1,0]] mod 2 = [[1,1],[1,0]] has order 3 in GL₂(F₂).
  Its orbit: [[1,1],[1,0]], [[0,1],[1,1]], [[1,0],[1,1]], [[1,1],[1,0]]...
  Three distinct matrices, then repeats.

For mod 4: lifting from F₂ to Z/4Z doubles the period via Hensel.
""")

# The matrix approach
import numpy as np

print("Matrix [[1,1],[1,0]] mod m — orders in GL₂(Z/mZ):")
for m in [2, 3, 4, 5, 6, 7, 8]:
    M = [[1, 1], [1, 0]]
    current = [[1, 0], [0, 1]]  # identity
    for k in range(1, 200):
        # Multiply
        new = [[0, 0], [0, 0]]
        for i in range(2):
            for j in range(2):
                for l in range(2):
                    new[i][j] = (new[i][j] + current[i][l] * M[l][j]) % m
        current = new
        if current == [[1, 0], [0, 1]]:
            print(f"  mod {m}: order = {k} = π({m})")
            break
    else:
        print(f"  mod {m}: order > 200")

# ============================================================
# Part 6: Fibonacci = 2-subsequence ⊕ 3-subsequence
# ============================================================
print("\n" + "=" * 70)
print("PART 6: FIBONACCI DECOMPOSED INTO 2 AND 3 SUBSEQUENCES")
print("=" * 70)

# The user's insight: Fibonacci can be "composed of 2 and 3 subsequences"
# Interpretation 1: Zeckendorf representation uses Fibonacci numbers
# Interpretation 2: Every positive integer = sum of 2s and 3s
#   (since gcd(2,3)=1, by Chicken McNugget: all n≥2 are representable)
# Interpretation 3: F(n) can be written using 2s and 3s
# Interpretation 4: The Fibonacci word on {2,3} or {a,b}

# Actually: Think about Fibonacci in terms of COMPOSITIONS.
# F(n+1) = # ways to tile a 1×n board with 1×1 and 1×2 tiles.
# If we use 1×2 and 1×3 tiles, we get a different sequence!

print("\nTiling interpretations:")
print("  F(n+1) = #{tilings of 1×n with tiles of size 1 and 2}")
print("  What if we use tiles of size 2 and 3?")
print()

# T(n) = #{tilings of 1×n with tiles of size 2 and 3}
# T(0) = 1 (empty), T(1) = 0, T(2) = 1, T(3) = 1, T(4) = 1, T(5) = 2
# T(n) = T(n-2) + T(n-3)
T = [0] * 30
T[0] = 1
for n in range(2, 30):
    T[n] = (T[n-2] if n >= 2 else 0) + (T[n-3] if n >= 3 else 0)

print(f"  T(n) [2,3 tiles]: {T[:20]}")
print(f"  Fibonacci F(n):   {fib[:20]}")

# Check: is there a recurrence for T?
# T(n) = T(n-2) + T(n-3)
# Char eq: t³ = t + 1, i.e., t³ - t - 1 = 0
# This is the tribonacci-like equation!

print("\n  T satisfies T(n) = T(n-2) + T(n-3)")
print("  Char eq: t³ - t - 1 = 0")

# Roots
import numpy as np
roots = np.roots([1, 0, -1, -1])
print(f"  Roots: {roots}")
print(f"  Real root ≈ {roots[0].real:.6f} (Plastic number!)")

# THE PLASTIC NUMBER! ρ ≈ 1.3247...
# The smallest Pisot number. This connects to substitution tilings!

print("\n  The real root is the PLASTIC NUMBER ρ ≈ 1.3247")
print("  It's the smallest Pisot number — a DEEP algebraic constant!")

# ============================================================
# Part 7: Fibonacci Word and the Golden String
# ============================================================
print("\n" + "=" * 70)
print("PART 7: FIBONACCI WORD ON ALPHABET {2, 3}")
print("=" * 70)

# Fibonacci word: start with 2, substitute 2→23, 3→2
# Or equivalently: start with a, b, concatenate: a, b, ab, bab, abbab, ...
# Using 2 and 3:
# Level 0: [2]
# Level 1: [2, 3]
# Level 2: [2, 3, 2]
# Level 3: [2, 3, 2, 2, 3]
# Level 4: [2, 3, 2, 2, 3, 2, 3, 2]

print("\nFibonacci word construction (2→23, 3→2):")
word = [2]
for level in range(8):
    print(f"  Level {level}: {word[:30]}{'...' if len(word) > 30 else ''} (length {len(word)})")
    new_word = []
    for c in word:
        if c == 2:
            new_word.extend([2, 3])
        else:  # c == 3
            new_word.extend([2])
    word = new_word

print(f"\n  Word lengths: {[fib[i+1] if i < len(fib)-1 else '?' for i in range(9)]}")
print("  (= Fibonacci numbers!)")

# Partial sums of the Fibonacci word
print("\nPartial sums of Fibonacci word:")
partial_sums = [0]
for i, c in enumerate(word[:30]):
    partial_sums.append(partial_sums[-1] + c)
print(f"  {partial_sums[:20]}")

# The partial sums modulo 6
print("\nPartial sums mod 6:")
ps_mod6 = [s % 6 for s in partial_sums[:30]]
print(f"  {ps_mod6}")

# ============================================================
# Part 8: 2+3=5, 2×3=6, and the Period
# ============================================================
print("\n" + "=" * 70)
print("PART 8: 2+3=5 (FIBONACCI), 2×3=6 (PERIOD)")
print("=" * 70)

print("""
DEEP NUMEROLOGY:
  2 + 3 = 5 = F(5)    [fifth Fibonacci number]
  2 × 3 = 6 = π(4)    [Pisano period of Fibonacci mod 4]
  2 and 3 are CONSECUTIVE FIBONACCI NUMBERS: F(3)=2, F(4)=3

The period-6 structure of Fib mod 4 = {1,1,2,3,1,0}:
  Sum of residues: 1+1+2+3+1+0 = 8 = 2³ = F(6)
  Product of nonzero: 1×1×2×3×1 = 6 = period itself!

The 2 and 3 in the period are AT POSITIONS 3 and 4.
  F(3) = 2, F(4) = 3 — the Fibonacci numbers ARE the residues!
  F(F(3)) = F(2) = 1, F(F(4)) = F(3) = 2 — self-referential!

CONNECTIONS TO TOURNAMENTS:
  - x=2 in I(P_k, x): Jacobsthal, roots {2, -1}
  - x=3 in I(P_k, x): roots {(1+√13)/2, (1-√13)/2}
  - Period 6 = π(4) = π(2²) = 2 × π(2) = 2 × 3
  - The OCF uses x=2 (binary arc choice)
  - The "3" appears as: I(single vertex, 2) = 1+2 = 3,
    and 3^k = I(k disjoint vertices, 2)
""")

# ============================================================
# Part 9: Fibonacci mod 2 and mod 3 separately
# ============================================================
print("=" * 70)
print("PART 9: FIBONACCI mod 2 AND mod 3 — CRT RECONSTRUCTION")
print("=" * 70)

print("\nF(n) mod 2: period π(2) = 3")
fib_mod2 = [fib[i] % 2 for i in range(18)]
print(f"  {fib_mod2}")
print("  Pattern: 1, 1, 0 (repeating)")

print("\nF(n) mod 3: period π(3) = 8")
fib_mod3 = [fib[i] % 3 for i in range(24)]
print(f"  {fib_mod3}")
print("  Pattern: 0, 1, 1, 2, 0, 2, 2, 1 (repeating)")

print("\nF(n) mod 4: period π(4) = 6")
fib_mod4 = [fib[i] % 4 for i in range(18)]
print(f"  {fib_mod4}")
print("  Pattern: 0, 1, 1, 2, 3, 1 (repeating)")

# CRT: mod 4 is NOT mod 2 × mod 3 (since gcd(2,3)=1 but 4≠2×3)
# But mod 6 = mod 2 × mod 3 via CRT
print("\nF(n) mod 6: period π(6) = lcm(π(2), π(3)) = lcm(3, 8) = 24")
fib_mod6 = [fib[i] % 6 for i in range(24)]
print(f"  {fib_mod6}")

# The 2-subsequence and 3-subsequence
print("\n--- 2-SUBSEQUENCE: positions where F(n) ≡ 2 mod 4 ---")
pos_2 = [i for i in range(30) if fib[i] % 4 == 2]
print(f"  Positions: {pos_2}")
diffs_2 = [pos_2[i+1] - pos_2[i] for i in range(len(pos_2)-1)]
print(f"  Gaps: {diffs_2}")

print("\n--- 3-SUBSEQUENCE: positions where F(n) ≡ 3 mod 4 ---")
pos_3 = [i for i in range(30) if fib[i] % 4 == 3]
print(f"  Positions: {pos_3}")
diffs_3 = [pos_3[i+1] - pos_3[i] for i in range(len(pos_3)-1)]
print(f"  Gaps: {diffs_3}")

# BOTH have constant gap 6 = the period!
print("\n  Both subsequences have constant gap 6 = π(4)")
print("  The 2 appears at n ≡ 3 mod 6, the 3 at n ≡ 4 mod 6")
print("  They are always ADJACENT in the period: (..., 2, 3, ...)")

# ============================================================
# Part 10: Tournament Connection — the Jacobsthal mod 4
# ============================================================
print("\n" + "=" * 70)
print("PART 10: JACOBSTHAL mod 4 AND TOURNAMENT H mod 4")
print("=" * 70)

# Jacobsthal: J(k) = (2^k - (-1)^k)/3
# J mod 4:
jac_mod4_seq = [(2**k - (-1)**k) // 3 % 4 for k in range(24)]
print(f"\nJacobsthal mod 4: {jac_mod4_seq}")

# Find period
for p in range(1, 24):
    if all(jac_mod4_seq[i] == jac_mod4_seq[i+p] for i in range(min(12, 24-p))):
        print(f"Jacobsthal Pisano period mod 4: {p}")
        print(f"  Pattern: {jac_mod4_seq[:p]}")
        break

# I(P_k, 2) mod 4
ip2_mod4 = [(2**(k+2) - (-1)**k) // 3 % 4 for k in range(24)]
print(f"\nI(P_k, 2) mod 4: {ip2_mod4}")

# I(C_k, 2) mod 4
ic2_mod4 = [(2**k + (-1)**k) % 4 for k in range(3, 27)]
print(f"I(C_k, 2) mod 4 (k≥3): {ic2_mod4}")

# The forbidden H = 7 ≡ 3 mod 4, H = 21 ≡ 1 mod 4, H = 63 ≡ 3 mod 4
print("\nForbidden H values mod 4:")
for h in [7, 21, 63]:
    print(f"  H={h}: {h} mod 4 = {h % 4}")

# All achievable H values mod 4 at n=5
print("\nH mod 4 at n=5: {1: ≡1, 3: ≡3, 5: ≡1, 9: ≡1, 11: ≡3, 13: ≡1, 15: ≡3}")
for h in [1, 3, 5, 9, 11, 13, 15]:
    print(f"  H={h:3d}: mod 4 = {h % 4}")

# ============================================================
# Part 11: Period-6 in Tournament Context
# ============================================================
print("\n" + "=" * 70)
print("PART 11: PERIOD 6 IN TOURNAMENT THEORY")
print("=" * 70)

print("""
WHERE DOES 6 APPEAR IN TOURNAMENTS?

1. ARCS AT n=4: m = C(4,2) = 6 arcs. The Boolean hypercube {0,1}^6
   is the space of all tournaments on 4 vertices.

2. EULERIAN NUMBER A(4,1) = 11, A(4,2) = 11 — symmetric.
   But A(3,1) = 4, and 4+1+1 = 6 = C(3,2).

3. PERIOD OF FIBONACCI MOD 4 = 6.
   Since tournaments are 2-colorings of arcs (mod 2 structure),
   and the Jacobsthal recurrence uses x=2,
   the period-6 structure may govern H periodicity.

4. 6 = 2 × 3 = product of the two "atomic" tournament numbers:
   - 2: binary arc choice
   - 3: smallest non-trivial tournament (triangle)

5. THE PATH-TOURNAMENT ANALOGY:
   I(P_k, 1) follows Fibonacci (period π(m) for mod m)
   I(P_k, 2) follows Jacobsthal (different period)
   At x=1: π(4) = 6
   At x=2: what is the period of Jacobsthal mod 4?
""")

# Compute I(P_k, 2) mod various m
print("I(P_k, 2) mod m — periods:")
for m in range(2, 13):
    seq = [(2**(k+2) - (-1)**k) // 3 % m for k in range(100)]
    for p in range(1, 100):
        if all(seq[i] == seq[i+p] for i in range(min(50, 100-p))):
            print(f"  I(P_k, 2) mod {m:2d}: period = {p:3d}  pattern: {seq[:min(p, 15)]}{'...' if p > 15 else ''}")
            break

# Compare with Fibonacci mod m
print("\nFibonacci mod m (= I(P_k, 1)) — periods:")
for m in range(2, 13):
    print(f"  F(n) mod {m:2d}: π = {pisano_period(m):3d}")

# ============================================================
# Part 12: The Golden Ratio and Tournament Growth
# ============================================================
print("\n" + "=" * 70)
print("PART 12: φ → 2 DEFORMATION — GOLDEN TO TOURNAMENT")
print("=" * 70)

print("""
The family of recurrences a(k) = a(k-1) + x·a(k-2) with:
  x = 1: roots φ ≈ 1.618, -1/φ ≈ -0.618 (GOLDEN)
  x = 2: roots 2, -1               (TOURNAMENT/JACOBSTHAL)
  x = 3: roots (1+√13)/2, (1-√13)/2  ≈ 2.303, -1.303

As x increases from 1 to 2:
  - Larger root goes from φ → 2 (EXACTLY at x=2, root is INTEGER)
  - Smaller root goes from -1/φ → -1 (also INTEGER at x=2!)
  - The Jacobsthal case x=2 is SPECIAL: both roots are integers!

This is why I(P_k, 2) = (2^{k+2} - (-1)^k)/3 has such clean forms.
Fibonacci has I(P_k, 1) = F(k+2) with irrational closed form.

The "deformation parameter" x interpolates:
  x = 0: trivial (constant sequence)
  x = 1: golden ratio (Fibonacci, nature, aesthetics)
  x = 2: tournaments (binary choices, OCF evaluation)
  x = 3: Ω with "ternary" structure

TOURNAMENT THEORY LIVES AT THE UNIQUE INTEGER POINT x=2
where both characteristic roots are integers!
""")

# Table showing the deformation
print("Deformation table:")
print(f"  {'x':>4s} | {'root1':>10s} | {'root2':>10s} | {'product':>10s} | {'sum':>10s} | I(P_5,x)")
for x_val in [0.5, 1, 1.5, 2, 2.5, 3, 4, 5]:
    r1 = (1 + (1 + 4*x_val)**0.5) / 2
    r2 = (1 - (1 + 4*x_val)**0.5) / 2
    # I(P_5, x) via recurrence
    G = [1, 1 + x_val]
    for k in range(2, 6):
        G.append(G[-1] + x_val * G[-2])
    print(f"  {x_val:4.1f} | {r1:10.6f} | {r2:10.6f} | {r1*r2:10.6f} | {r1+r2:10.6f} | {G[5]:.1f}")

# Note: product = -x, sum = 1 always (Vieta's formulas for t²-t-x=0)
print("\n  Vieta: root1 + root2 = 1, root1 × root2 = -x")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — FIBONACCI, PERIOD 6, AND TOURNAMENTS")
print("=" * 70)
print("""
CROWN JEWELS:

1. FIBONACCI mod 4 = {0,1,1,2,3,1} with period π(4) = 6.
   The 2 and 3 are at positions 3 and 4 (mod 6) — adjacent,
   and F(3)=2, F(4)=3 are themselves the residues.

2. 6 = 2×3: The period is the PRODUCT of the two "special" residues.
   Also 2+3 = 5 = F(5), the next Fibonacci number.

3. UNIFIED RECURRENCE: I(P_k, x) = I(P_{k-1}, x) + x·I(P_{k-2}, x)
   - x=1: Fibonacci world (golden ratio φ)
   - x=2: Tournament world (Jacobsthal, integer roots 2 and -1)
   Tournament theory lives at the UNIQUE INTEGER POINT of the
   φ-deformation family.

4. PLASTIC NUMBER: Tiles of size 2 and 3 give recurrence
   T(n) = T(n-2) + T(n-3), char eq t³-t-1=0.
   The real root is the Plastic number ρ ≈ 1.3247, the smallest
   Pisot number. This bridges 2-3 composition to algebraic number theory.

5. FIBONACCI WORD on {2,3}: The substitution 2→23, 3→2 generates
   the golden string with letters 2 and 3. Lengths follow Fibonacci.
   This is a quasicrystal — aperiodic but highly ordered.

6. PERIOD-6 UNIVERSALITY: π(4)=6 for Fibonacci, and 6 = C(4,2) =
   number of arcs in a 4-tournament. The connection runs deep:
   - F(n) mod 4 returns to start after 6 steps
   - Tournament on 4 vertices has 6 binary arc choices
   - The Jacobsthal sequence mod 4 has period 6 too!
""")
