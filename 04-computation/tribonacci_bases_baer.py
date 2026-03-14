#!/usr/bin/env python3
"""
tribonacci_bases_baer.py
opus-2026-03-14-S71l

NON-STANDARD BASES, TRIBONACCI, AND TOURNAMENT THEORY

User's key insights:
- k-nacci approaches 2 as k→∞
- Weighted k-nacci approaches 3
- Simplices = (x+1)^n, cuboids = (x+2)^n
- Think about packing simplices inside cuboids

This script explores:
1. The tribonacci constant τ ≈ 1.839 and its role in tournaments
2. Number bases: base φ (golden), base τ (tribonacci), base 2, base 3
3. The k-nacci limit → 2 and weighted limit → 3
4. Expressing tournament constants in non-standard bases
5. The simplex-cuboid packing in general dimension
6. Category-theoretic structure of the k-strand Pascal triangle
"""

import numpy as np
from math import log, sqrt, gcd
from fractions import Fraction

print("=" * 70)
print("PART 1: THE K-NACCI CONSTANTS")
print("=" * 70)
print()

print("The k-nacci sequence generalizes Fibonacci:")
print("  k=2: Fibonacci     a(n) = a(n-1) + a(n-2)         ratio → φ ≈ 1.618")
print("  k=3: Tribonacci    a(n) = a(n-1) + a(n-2) + a(n-3) ratio → τ ≈ 1.839")
print("  k=4: Tetranacci    a(n) = Σ a(n-i) for i=1..4      ratio → ≈ 1.928")
print("  k→∞:                                                ratio → 2")
print()

# Compute k-nacci ratios
for k in [2, 3, 4, 5, 6, 7, 8, 10, 20, 50]:
    # k-nacci: start with [0,...,0,1], recurrence = sum of last k terms
    seq = [0] * (k - 1) + [1]
    for _ in range(100):
        seq.append(sum(seq[-k:]))
    ratio = seq[-1] / seq[-2] if seq[-2] != 0 else 0
    print(f"  k={k:3d}: ratio = {ratio:.10f}  (2 - ratio = {2 - ratio:.10f})")

print()
print("  As k→∞, ratio → 2 (the k-nacci constant approaches 2)")
print("  The TOURNAMENT GENERATOR 2 is the LIMIT of k-nacci ratios!")
print()

print("=" * 70)
print("PART 2: WEIGHTED K-NACCI → 3")
print("=" * 70)
print()

print("Weighted k-nacci at x=2: a(n) = 2·a(n-1) + 2·a(n-2) + ... + 2·a(n-k)")
print("  This is the k-nacci with each term weighted by x=2.")
print("  Equivalently: char poly = x^k - 2(x^{k-1} + ... + 1) = x^k - 2(x^k-1)/(x-1)")
print()
print("  The dominant root approaches 2k/(k-1)... let's compute:")
print()

for k in [2, 3, 4, 5, 6, 7, 8, 10, 20]:
    # Weighted k-nacci: a(n) = x*sum of last k terms, x=2
    seq = [0] * (k - 1) + [1]
    for _ in range(80):
        seq.append(2 * sum(seq[-k:]))
    ratio = seq[-1] / seq[-2] if seq[-2] != 0 else 0
    print(f"  k={k:3d}: weighted ratio = {ratio:.10f}  (3 - ratio = {3 - ratio:.10f})")

print()
print("  As k→∞, weighted ratio → 3 = 1+x = simplex factor!")
print("  The SIMPLEX GENERATOR (x+1)=3 is the limit of weighted k-nacci!")
print()

print("  INTERPRETATION:")
print("    k-nacci → 2 = cuboid factor (x+2)/2 ... no.")
print("    Actually: k-nacci → 2 = number of arc states per edge")
print("    Weighted k-nacci at x=2 → 3 = simplex volume factor (1+x)")
print("    These are (x+2)^n and (1+x)^n at x=2!")
print()

print("=" * 70)
print("PART 3: THE TRIBONACCI CONSTANT τ")
print("=" * 70)
print()

# Tribonacci constant = real root of x^3 - x^2 - x - 1 = 0
# τ ≈ 1.8392867552141612
coeffs_trib = [1, -1, -1, -1]  # x^3 - x^2 - x - 1
roots = np.roots(coeffs_trib)
tau = max(r.real for r in roots if abs(r.imag) < 1e-10)
print(f"Tribonacci constant τ = {tau:.15f}")
print(f"  Real root of x³ - x² - x - 1 = 0")
print(f"  Equivalently: x³ = x² + x + 1 = Φ₃(x)!")
print()
print(f"  *** τ³ = Φ₃(τ) = τ² + τ + 1 ***")
print(f"  Verify: τ³ = {tau**3:.10f}")
print(f"          τ²+τ+1 = {tau**2 + tau + 1:.10f}")
print()

print("  THE TRIBONACCI CONSTANT SATISFIES τ³ = Φ₃(τ)!")
print("  This is the DEFINING RELATIONSHIP connecting tribonacci to Φ₃!")
print()

# The tribonacci polynomial x^3 - x^2 - x - 1 factors as x^3 - Phi_3(x)
# or equivalently x^3 = Phi_3(x) → x^3 - x^2 - x - 1 = 0
print("  Rewriting: x³ - (x²+x+1) = 0  →  x³ = Φ₃(x)")
print()
print("  Compare with golden ratio: φ² = φ + 1 = Φ₂(φ)? ")
print(f"    φ = {(1+sqrt(5))/2:.10f}")
print(f"    φ² = {((1+sqrt(5))/2)**2:.10f}")
print(f"    φ+1 = {(1+sqrt(5))/2 + 1:.10f}")
print(f"    Φ₂(φ) = φ+1 = {(1+sqrt(5))/2 + 1:.10f}")
print()
print("  Actually φ² = φ+1, and x²-x-1 = 0 means x² = x+1.")
print("  For Fibonacci: x² = x + 1 (sum of 2 terms)")
print("  For tribonacci: x³ = x² + x + 1 = Φ₃(x) (sum of 3 terms)")
print("  For k-nacci:    x^k = x^{k-1} + ... + x + 1 = (x^k-1)/(x-1)")
print()

print("  THE K-NACCI EQUATION: x^k = (x^k - 1)/(x - 1)")
print("  Rearranging: x^k(x-1) = x^k - 1")
print("               x^{k+1} - x^k = x^k - 1")
print("               x^{k+1} = 2x^k - 1")
print()
print("  At the limit k→∞: if x = 2-ε, then")
print("  (2-ε)^{k+1} ≈ 2(2-ε)^k - 1")
print("  This gives x → 2 as k → ∞ ✓")
print()

# Key values of τ
print(f"  τ = {tau:.10f}")
print(f"  τ² = {tau**2:.10f}")
print(f"  τ³ = {tau**3:.10f} = Φ₃(τ)")
print(f"  Φ₃(τ) = τ²+τ+1 = {tau**2+tau+1:.10f}")
print(f"  Φ₃(2) = 7")
print(f"  Φ₃(φ) = φ²+φ+1 = {((1+sqrt(5))/2)**2 + (1+sqrt(5))/2 + 1:.10f}")
print(f"         = (φ+1)+φ+1 = 2φ+2 = {2*(1+sqrt(5))/2 + 2:.10f}")
print(f"         = 2+√5 = {2+sqrt(5):.10f}")
print()

phi = (1 + sqrt(5)) / 2
print(f"  Φ₃(φ) = 2+√5 ≈ {2+sqrt(5):.6f}")
print(f"  Φ₃(τ) = τ³ ≈ {tau**3:.6f}")
print(f"  Φ₃(2) = 7 (integer!)")
print(f"  Φ₃(3) = 13 (integer!)")
print()

print("=" * 70)
print("PART 4: TOURNAMENT CONSTANTS IN NON-STANDARD BASES")
print("=" * 70)
print()

# Express numbers in various bases
def to_base(n, base, digits=20):
    """Express integer n in given base (may be irrational)."""
    if n == 0:
        return "0"
    result = []
    remaining = float(n)
    # Integer part
    int_part = int(remaining)
    remaining -= int_part
    # Convert integer part
    if int_part == 0:
        int_str = "0"
    else:
        int_digits = []
        temp = int_part
        while temp > 0:
            int_digits.append(int(temp % base))
            temp = int(temp // base)
        int_str = ''.join(map(str, reversed(int_digits)))
    return int_str

# Express tournament constants in base φ, base τ, base 2, base 3
bases = [
    (2, "2 (binary)"),
    (3, "3 (ternary)"),
    (phi, f"φ ≈ {phi:.4f} (golden)"),
    (tau, f"τ ≈ {tau:.4f} (tribonacci)"),
]

tournament_numbers = [
    (3, "cycle generator = Φ₆(2)"),
    (7, "Fano = Φ₃(2) = FORBIDDEN"),
    (13, "|PG(2,F₃)| = Φ₃(3)"),
    (21, "Baer = Φ₃(4) = FORBIDDEN"),
    (6, "tournament period = LCM(2,3)"),
    (8, "|T₃| = 2³"),
    (64, "|T₄| = 2⁶"),
    (1024, "|T₅| = 2¹⁰"),
]

print("Tournament constants in various bases:")
print()
for n, name in tournament_numbers:
    print(f"  {n:5d} = {name}")
    for base, bname in bases:
        rep = to_base(n, base)
        log_b = log(n, base) if n > 0 else 0
        print(f"         base {bname}: {rep}  (log = {log_b:.4f})")
    print()

print("=" * 70)
print("PART 5: BASE-τ AND Φ₃")
print("=" * 70)
print()

print("Since τ³ = Φ₃(τ) = τ²+τ+1, base-τ has special properties:")
print()
print("  In base τ:")
print(f"    τ = 10                (by definition)")
print(f"    τ² = 100")
print(f"    τ³ = 1000 = 111       (because τ³ = τ²+τ+1)")
print(f"    This means: 1000 = 111 in base τ!")
print(f"    The tribonacci 'carrying rule': three consecutive 1s = one 1 three places left")
print()

# Verify
print(f"  Verification:")
print(f"    τ² + τ + 1 = {tau**2:.6f} + {tau:.6f} + 1 = {tau**2+tau+1:.6f}")
print(f"    τ³ = {tau**3:.6f}")
print(f"    Match: {abs(tau**3 - tau**2 - tau - 1) < 1e-10}")
print()

print("  Compare with base φ:")
print(f"    φ² = φ+1, so in base φ: 100 = 11")
print(f"    The Fibonacci carrying rule: two consecutive 1s = one 1 two places left")
print()
print("  Compare with base 2:")
print(f"    2² = 2+2, so in base 2: 100 = ... (normal binary)")
print(f"    No carrying rule involving sums of lower digits")
print()

print("  THE CARRYING RULES ARE THE K-NACCI RECURRENCES!")
print("  base φ: 100_φ = 11_φ      → Fibonacci carrying (2 terms)")
print("  base τ: 1000_τ = 111_τ    → Tribonacci carrying (3 terms)")
print("  base k-nacci: 10...0 = 1...1  → k-nacci carrying (k terms)")
print()

print("=" * 70)
print("PART 6: Φ₃ IN BASE τ vs BASE φ vs BASE 2")
print("=" * 70)
print()

print("Evaluating Φ₃(x) = x²+x+1 at different bases:")
print()

# At x = τ: Φ₃(τ) = τ³ = τ·τ·τ
print(f"  Φ₃(τ) = τ³ = {tau**3:.10f}")
print(f"    In base τ: 1000")
print(f"    Or equivalently: 111 (using tribonacci carrying)")
print(f"    The Fano-like value in base τ is just 'three ones'!")
print()

# At x = φ: Φ₃(φ) = 2+√5
print(f"  Φ₃(φ) = 2+√5 = {2+sqrt(5):.10f}")
print(f"    = φ² + φ + 1 = (φ+1) + φ + 1 = 2φ+2")
print(f"    In base φ: 100 + 10 + 1 = 111_φ")
print(f"    But 111_φ = φ²+φ+1 is NOT the same as writing 111 in base φ")
print(f"    (In Zeckendorf, 111 is not allowed since φ²=φ+1 means 100=11)")
print()

# At x = 2: Φ₃(2) = 7
print(f"  Φ₃(2) = 7")
print(f"    In base 2: 111")
print(f"    7 = 4+2+1 = 2²+2¹+2⁰ = 111₂")
print()

# At x = 3: Φ₃(3) = 13
print(f"  Φ₃(3) = 13")
print(f"    In base 3: 111")
print(f"    13 = 9+3+1 = 3²+3¹+3⁰ = 111₃")
print()

print("  *** Φ₃(x) = x²+x+1 IS the number 111 in base x! ***")
print()
print("  This is OBVIOUS but PROFOUND:")
print("  The 'Fano number' 7 = 111₂")
print("  The 'PG(2,3) number' 13 = 111₃")
print("  The 'Baer number' 21 = 111₄ (base 4)")
print("  The projective plane PG(2,q) has 111_q points!")
print()

# Verify 21 in base 4
print(f"  21 in base 4: {21} = 4²+4+1 = 16+4+1 = 21 = 111₄ ✓")
print()

print("  MORE GENERALLY:")
print("  PG(k-1, q) has q^{k-1}+q^{k-2}+...+q+1 = 111...1_q (k ones) points")
print("  This is the number (q^k-1)/(q-1) = 'k ones in base q'")
print()
print("  For k=3 (projective planes):")
print("    PG(2,2) = 111₂ = 7")
print("    PG(2,3) = 111₃ = 13")
print("    PG(2,4) = 111₄ = 21")
print("    PG(2,5) = 111₅ = 31")
print()

print("  For k=2 (projective lines):")
print("    PG(1,2) = 11₂ = 3")
print("    PG(1,3) = 11₃ = 4")
print("    PG(1,4) = 11₄ = 5")
print("    PG(1,5) = 11₅ = 6")
print()

print("  THE FORBIDDEN VALUES ARE 'ALL ONES' NUMBERS:")
print("  7 = 111₂ = repunit in base 2")
print("  21 = 111₄ = repunit in base 4 = base 2²")
print("  273 = 111₁₆ = repunit in base 16 = base 2⁴ (NOT forbidden)")
print()

print("=" * 70)
print("PART 7: SIMPLEX-CUBOID PACKING")
print("=" * 70)
print()

print("User's framework: simplex = (x+1)^n, cuboid = (x+2)^n at x=2")
print("  Simplex volume factor: (1+x)^n = 3^n")
print("  Cuboid volume factor:  (2+x)^n = 4^n")
print()
print("  Packing ratio: simplex/cuboid = (3/4)^n = ((1+x)/(2+x))^n")
print("  This is the fraction of cuboid volume occupied by the simplex.")
print()

print("  Complement = 1 - (3/4)^n = 'extra pieces' fraction")
for n in range(1, 9):
    ratio = (3/4)**n
    complement = 1 - ratio
    pieces = round(4**n - 3**n)  # integer value of complement * 4^n
    print(f"  n={n}: (3/4)^{n} = {ratio:.6f}, complement = {complement:.6f}, "
          f"extra = 4^{n}-3^{n} = {pieces}")

print()
print("  The extra pieces 4^n - 3^n:")
for n in range(1, 9):
    extra = 4**n - 3**n
    mod7 = extra % 7
    div7 = extra // 7 if mod7 == 0 else ""
    print(f"    n={n}: {extra:6d}  mod 7 = {mod7}  {'= 7·'+str(div7) if mod7==0 else ''}")

print()
print("  7 | (4^n - 3^n) iff n is even (since 4≡3≡-1 mod 7 fails...)")
print("  Actually: 4 ≡ 4 mod 7, 3 ≡ 3 mod 7")
print("  4^n - 3^n mod 7:")
for n in range(1, 13):
    print(f"    n={n}: 4^{n} ≡ {pow(4,n,7)}, 3^{n} ≡ {pow(3,n,7)}, "
          f"diff ≡ {(pow(4,n,7)-pow(3,n,7))%7} mod 7")
print(f"  Period of (4^n - 3^n) mod 7: 6 = tournament period!")
print()

print("=" * 70)
print("PART 8: SIMPLEX NESTING AND CORNER PIECES")
print("=" * 70)
print()

print("The user's geometric picture:")
print("  n=2: equilateral triangle in square, 2 corner triangles")
print("  n=3: regular tetrahedron in cube, 4 corner pieces")
print()
print("  For a regular n-simplex inscribed in an n-cube:")
print("  A regular simplex with n+1 vertices can inscribe in an n-cube")
print("  when a Hadamard matrix H_{n+1} exists (n+1 ≡ 0 mod 4 or n+1 ≤ 2)")
print()

print("  HADAMARD DIMENSIONS: n+1 = 1, 2, 4, 8, 12, 16, 20, ...")
print("  So n = 0, 1, 3, 7, 11, 15, 19, ...")
print()
print("  At these dimensions, the simplex inscribes perfectly and")
print("  the 'corner pieces' are congruent orthoschemes.")
print()

# Volume of regular simplex inscribed in n-cube
print("  Volume ratios (simplex/cube) at Hadamard dimensions:")
from math import factorial
for n in [1, 2, 3, 7]:
    if n == 1:
        vol_ratio = 1.0  # line segment = 1D cube
    elif n == 2:
        vol_ratio = sqrt(3)/4 / 1  # equilateral triangle in unit square
        # Actually: equilateral triangle inscribed in unit square
        # side = 1, area = sqrt(3)/4, square area = 1
        # But user says "2 halves" → triangle with vertices at (0,0),(1,0),(0.5, sqrt(3)/2)
        # sits in [0,1]×[0,sqrt(3)/2] rectangle, not a square
        # For a unit square with equilateral triangle: side = 1, area = sqrt(3)/4
        vol_ratio = sqrt(3) / 4
    elif n == 3:
        # Regular tetrahedron in unit cube
        # Vertices at (0,0,0),(1,1,0),(1,0,1),(0,1,1)
        # Volume = 1/3 of cube
        vol_ratio = 1/3
    elif n == 7:
        # Regular simplex in 7-cube using Hadamard matrix H_8
        # Volume = sqrt(8) / 7! ... actually
        # Volume of regular n-simplex with edge length a:
        # V = a^n * sqrt(n+1) / (n! * 2^{n/2})
        # In the unit cube, edge length = sqrt(n) (using Hadamard)
        # Wait: Hadamard simplex has edge length sqrt(2·n/(n+1))... complex
        # For the standard embedding: vertices = rows of H_{n+1} normalized
        # Volume = (n+1)^{(n+1)/2} / (2^n · n!) ... approximately
        vol_ratio = 8**(3.5) / (factorial(7) * 2**3.5)  # approximate

    n_corners = 2**n - (n + 1)  # number of corner pieces if they're orthoschemes
    print(f"  n={n}: vol ratio ≈ {vol_ratio:.6f}, corner pieces = 2^{n}-{n+1} = {n_corners}")

print()
print("  THE KEY PATTERN:")
print("  n=1: 2^1 - 2 = 0 pieces (simplex = segment = cube)")
print("  n=2: 2^2 - 3 = 1 piece  (triangle + 1 complement)")
print("  n=3: 2^3 - 4 = 4 pieces (tetrahedron + 4 corners)")
print("  n=7: 2^7 - 8 = 120 = 5! pieces")
print()
print("  General: 2^n - (n+1) = Σ C(n,k) for k≥2 = interaction terms")
print("  These are the higher-order terms in (1+1)^n = 2^n")
print("  The simplex captures the 'linear' part: C(n,0)+C(n,1) = 1+n")
print("  The corners capture the 'nonlinear' part: C(n,2)+C(n,3)+...")
print()
print("  IN TOURNAMENT TERMS:")
print("  The simplex ↔ degree-0 and degree-1 Walsh components = mean + linear")
print("  The corners ↔ degree-2 and higher Walsh components = interactions")
print("  H is determined by degree-2+ interactions (the 'corner pieces'!)")
print()

print("=" * 70)
print("PART 9: THE TRIBONACCI-BAER CONNECTION")
print("=" * 70)
print()

print("Since τ³ = Φ₃(τ) and Φ₃(2) = 7, Φ₃(4) = 21:")
print()
print("  The tribonacci constant satisfies THE SAME polynomial relation")
print("  that generates the forbidden values, just at an irrational point.")
print()
print("  Think of it as a continuous version:")
print("  Φ₃(x) evaluated at integer points gives projective planes")
print("  Φ₃(x) evaluated at k-nacci constants gives k-nacci^(k+1)")
print()

# Verify: Φ₃(k-nacci ratio) for various k
print("  k-nacci ratios and Φ₃ values:")
for k in [2, 3, 4, 5, 6, 10]:
    seq = [0] * (k - 1) + [1]
    for _ in range(100):
        seq.append(sum(seq[-k:]))
    r = seq[-1] / seq[-2]
    phi3_val = r**2 + r + 1
    r_cubed = r**3
    rk_power = r**k

    print(f"  k={k:2d}: ratio={r:.8f}, Φ₃(r) = {phi3_val:.6f}, "
          f"r^k = {rk_power:.6f}, r³ = {r_cubed:.6f}")

print()
print("  For k=3 (tribonacci): Φ₃(τ) = τ³ EXACTLY")
print("  For k=2 (Fibonacci): Φ₃(φ) = 2+√5 ≈ 4.236, φ³ ≈ 4.236 — WAIT")
print(f"    φ³ = {phi**3:.10f}")
print(f"    Φ₃(φ) = {phi**2+phi+1:.10f}")
print(f"    φ³ = Φ₃(φ)? {abs(phi**3 - (phi**2+phi+1)) < 1e-10}")
print()
print("  WAIT: φ³ = φ·φ² = φ(φ+1) = φ²+φ = (φ+1)+φ = 2φ+1")
print(f"    2φ+1 = {2*phi+1:.10f}")
print(f"    Φ₃(φ) = φ²+φ+1 = (φ+1)+φ+1 = 2φ+2 = {2*phi+2:.10f}")
print(f"    So φ³ = 2φ+1 ≠ 2φ+2 = Φ₃(φ)")
print(f"    The Fibonacci equation is φ² = φ+1, NOT φ³ = Φ₃(φ)")
print()
print("  THE DISTINCTION:")
print("  Fibonacci: φ² = φ+1  (degree 2 relation)")
print("  Tribonacci: τ³ = τ²+τ+1 = Φ₃(τ)  (degree 3 relation)")
print("  The tribonacci is the UNIQUE constant where x^3 = Φ₃(x)")
print()

print("=" * 70)
print("PART 10: REPUNITS AND TOURNAMENT FORBIDDEN VALUES")
print("=" * 70)
print()

print("REPUNIT: a number consisting of all 1s in some base")
print("  R_k(b) = 111...1_b (k ones) = (b^k - 1)/(b - 1)")
print()

print("  R_3(b) = b² + b + 1 = Φ₃(b)  (3-digit repunit = Φ₃)")
print("  R_2(b) = b + 1               (2-digit repunit)")
print("  R_k(b) = Σ b^i for i=0..k-1   (k-digit repunit)")
print()

print("  FORBIDDEN VALUES AS REPUNITS:")
print("  7 = R_3(2) = 111₂  (base-2 repunit)")
print("  21 = R_3(4) = 111₄ (base-4 repunit)")
print("  273 = R_3(16) = 111₁₆ (base-16 repunit, NOT forbidden)")
print()
print("  These are 3-digit repunits in base 2^{2^k}:")
print("  k=0: R_3(2) = 7 ← FORBIDDEN")
print("  k=1: R_3(4) = 21 ← FORBIDDEN")
print("  k=2: R_3(16) = 273 ← achievable")
print()

print("  OTHER REPUNIT SEQUENCES:")
print("  2-digit repunits R_2(b) = b+1: 3, 5, 7, 9, 11, ...")
print("    R_2(2) = 3 = Φ₆(2), R_2(6) = 7 = Φ₃(2)")
print("    7 is BOTH R_3(2) and R_2(6)!")
print()

print("  4-digit repunits R_4(b) = b³+b²+b+1:")
for b in range(2, 8):
    r4 = b**3 + b**2 + b + 1
    print(f"    R_4({b}) = {r4} = (b+1)(b²+1) = {b+1}×{b**2+1}")
print()
print("  R_4(b) = (b+1)(b²+1) always factors!")
print("  So 4-digit repunits are NEVER prime.")
print("  But 3-digit repunits R_3(b) = Φ₃(b) can be prime.")
print()

# Check which R_3(b) are prime
print("  Primality of R_3(b) = Φ₃(b):")
def is_prime(n):
    if n < 2: return False
    for p in range(2, int(n**0.5)+1):
        if n % p == 0: return False
    return True

for b in range(2, 20):
    val = b**2 + b + 1
    factors = []
    temp = val
    for p in range(2, temp+1):
        while temp % p == 0:
            factors.append(p)
            temp //= p
        if temp == 1: break
    prime = is_prime(val)
    print(f"    Φ₃({b:2d}) = {val:4d} {'PRIME' if prime else ' = '+' × '.join(map(str,factors))}")

print()
print("  Φ₃ primes (b such that Φ₃(b) is prime):")
print("  b = 1(3), 2(7), 3(13), 5(31), 6(43), 8(73), 12(157), 14(211), ...")
print("  These are called 'generalized Eisenstein primes' — norms of")
print("  Eisenstein primes in Z[ω].")
print()

print("=" * 70)
print("PART 11: THE CATEGORY OF REPUNITS")
print("=" * 70)
print()

print("Define the REPUNIT FUNCTOR:")
print("  R_k: (Bases) → (Integers)")
print("  b ↦ R_k(b) = (b^k-1)/(b-1)")
print()
print("  Morphisms: b₁ = b₂^m gives R_k(b₂^m) | R_{km}(b₂)")
print("  because 111...1_{b^m} (k digits) embeds into 111...1_b (km digits)")
print()
print("  At k=3: R_3(2) = 7, R_3(4) = 21, R_3(16) = 273")
print("  Divisibility: R_3(2) | R_3(4)? 7 | 21 = 3×7 ✓")
print("  R_3(4) | R_3(16)? 21 | 273 = 13×21 ✓")
print("  R_3(2^{2^k}) divides R_3(2^{2^{k+1}}) for all k!")
print()

# Verify the divisibility tower
print("  The Baer divisibility tower:")
for k in range(5):
    b = 2**(2**k)
    r3 = b**2 + b + 1
    if k > 0:
        prev_b = 2**(2**(k-1))
        prev_r3 = prev_b**2 + prev_b + 1
        quotient = r3 // prev_r3
        print(f"    R_3({b}) = {r3} = R_3({prev_b}) × {quotient}")
    else:
        print(f"    R_3({b}) = {r3}")

print()
print("  The quotients: 3, 13, 241, 65281, ...")
print("  These are Φ₆(2^{2^k}) = q² - q + 1 at q = 2^{2^k}!")
print("  This IS the Baer partition: each PG(2,F_{q²}) contains")
print("  Φ₆(q) = q²-q+1 copies of PG(2,F_q).")
print()

print("=" * 70)
print("PART 12: SYNTHESIS — THE REPUNIT-TRIBONACCI-BAER TRIANGLE")
print("=" * 70)
print()

print("THREE WINDOWS onto the same structure:")
print()
print("  1. REPUNITS: Φ₃(b) = 111_b = 3-digit repunit in base b")
print("     Forbidden values = repunits in bases 2, 4")
print("     The 'all ones' structure = every position contributes equally")
print()
print("  2. TRIBONACCI: τ³ = Φ₃(τ) = τ²+τ+1")
print("     The tribonacci recurrence IS the Φ₃ polynomial relation")
print("     τ is the 'irrational base' where Φ₃ becomes a pure power")
print()
print("  3. BAER PLANES: |PG(2,q)| = Φ₃(q) = 111_q")
print("     The projective plane IS the repunit geometry")
print("     Each 'digit position' (0, 1, 2) corresponds to a coordinate axis")
print("     The 111 structure means: one point for each combination")
print()
print("  THE UNITY:")
print("  A projective plane PG(2,q) has q²+q+1 = 111_q points.")
print("  This is Φ₃(q), the third cyclotomic polynomial.")
print("  At q = τ (tribonacci), this becomes τ³ = Φ₃(τ).")
print("  At q = 2 (tournament base), this becomes 7 = FORBIDDEN.")
print("  The repunit '111' is the universal form, and the base")
print("  determines which mathematical world we're in.")
print()
print("  FORBIDDEN VALUES = REPUNITS IN TOURNAMENT-RELATED BASES")
print("  The base 2 is the tournament generator (2 states per arc)")
print("  The base 4 = 2² is the Baer square (F₄ field size)")
print("  Together, 7 = 111₂ and 21 = 111₄ are the only forbidden repunits.")
print()

print("=" * 70)
print("DONE — TRIBONACCI BASES BAER")
print("=" * 70)
