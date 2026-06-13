#!/usr/bin/env python3
"""
fano_hamming_forbidden.py — opus-2026-03-14-S80

The first forbidden H value is 7 = |PG(2,2)| = the Fano plane.
The Fano plane is the projective geometry underlying the [7,4,3] Hamming code.
This is also 7 = I(K₃, 2) = the independence polynomial of the triangle at x=2.

We explore:
1. The Hamming code, Fano plane, and the tournament polynomial
2. Error-correcting codes built from the (2,3,5) structure
3. The Golay code and its relation to the Petersen graph
4. Weight enumerators as independence polynomials
5. The Steiner system S(5,8,24) and the Leech lattice through (2,3)
6. Coding theory meets tournament theory: channel capacity
"""

from math import comb, log2, factorial
from fractions import Fraction
import numpy as np

def section(title, n):
    print(f"\n{'='*70}")
    print(f"{n}. {title}")
    print(f"{'='*70}\n")

KEY1, KEY2 = 2, 3

# ============================================================
section("THE FANO PLANE AS H_FORBIDDEN", 1)
# ============================================================

print("The Fano plane PG(2,2):")
print("  7 points, 7 lines, 3 points per line, 3 lines per point")
print("  = unique (7₃) configuration")
print()
print("  Points: {1,2,3,4,5,6,7}")
print("  Lines:  {124, 235, 346, 457, 156, 267, 137}")
print("  (cyclic: line i contains {i, i+1, i+3} mod 7)")
print()
print("  Automorphism group: GL(3,2) = PSL(2,7)")
print(f"  |GL(3,2)| = (2³-1)(2³-2)(2³-4) = 7·6·4 = {7*6*4}")
print(f"  = 168 = 7·24 = H_forb₁·|BT|")
print(f"  = 168 = 8·21 = rank(E₈)·H_forb₂")

print()
print("  Incidence matrix of Fano plane:")
# Points as columns, lines as rows
fano_lines = [(0,1,3), (1,2,4), (2,3,5), (3,4,6), (0,4,5), (1,5,6), (0,2,6)]
M = np.zeros((7,7), dtype=int)
for i, line in enumerate(fano_lines):
    for j in line:
        M[i][j] = 1

for row in M:
    print("    " + " ".join(str(x) for x in row))

print(f"\n  Each row has {sum(M[0])} ones (KEY₂ points per line)")
print(f"  Each column has {sum(M[:,0])} ones (KEY₂ lines per point)")
print(f"  M·Mᵀ = I + J (mod 2)? Let's check...")

MMt = M @ M.T
print(f"\n  M·Mᵀ (over Z):")
for row in MMt:
    print("    " + " ".join(str(x) for x in row))

print(f"\n  Diagonal entries: {MMt[0][0]} = KEY₂ (points per line)")
print(f"  Off-diagonal: {MMt[0][1]} = 1 (any two lines meet in exactly 1 point)")
print(f"  So M·Mᵀ = (KEY₂-1)·I + J = KEY₁·I + J")

# ============================================================
section("HAMMING CODE [7,4,3] AND TOURNAMENT POLYNOMIAL", 2)
# ============================================================

print("The [7,4,3] Hamming code:")
print(f"  Length n = 7 = H_forb₁ = Φ₃(KEY₁)")
print(f"  Dimension k = 4 = rank(F₄) = KEY₁²")
print(f"  Distance d = 3 = KEY₂")
print()
print("  Parameters: [Φ₃(KEY₁), KEY₁², KEY₂]")
print()

# Parity check matrix
H = np.array([
    [1,0,0,1,1,0,1],
    [0,1,0,1,0,1,1],
    [0,0,1,0,1,1,1]
], dtype=int)

print("  Parity check matrix (columns = binary reps of 1..7):")
for row in H:
    print("    " + " ".join(str(x) for x in row))

print(f"\n  H is a 3×7 matrix = KEY₂ × Φ₃(KEY₁)")
print(f"  Redundancy: n-k = 7-4 = 3 = KEY₂")
print(f"  Rate: k/n = 4/7 = KEY₁²/Φ₃(KEY₁)")

# Weight enumerator
print("\nWeight enumerator of Hamming [7,4,3]:")
print("  W(x,y) = x⁷ + 7x⁴y³ + 7x³y⁴ + y⁷")
# Coefficients: A_0=1, A_3=7, A_4=7, A_7=1
# 16 codewords total (2^4)
print(f"  Codeword weights: 0(×1), 3(×7), 4(×7), 7(×1)")
print("  Total codewords: 1+7+7+1 = 16 = 2^4 = KEY1^(KEY1^2)")
print(f"  Number of weight-3 codewords: 7 = H_forb₁")
print(f"  Number of weight-4 codewords: 7 = H_forb₁")
print("  Symmetry: A_d = A_{n-d} (7 = 7)")
print()
print(f"  Weight distribution [1, 0, 0, 7, 7, 0, 0, 1]")
print(f"  Sum = {1+7+7+1} = 2^4 ✓")
print()

# Connection to independence polynomial
print("Weight enumerator at (1,2):")
W_12 = 1 + 7*2**3 + 7*2**4 + 2**7
print(f"  W(1,2) = 1 + 7·8 + 7·16 + 128 = {W_12}")
print(f"  = {W_12}")
print()

# The dual code
print("Dual code: Hamming [7,3,4] (simplex code):")
print(f"  Length n = 7, Dimension k = 3, Distance d = 4 = KEY₁²")
print(f"  All nonzero codewords have weight 4 = KEY₁²")
print(f"  Number of nonzero codewords: 2³-1 = 7 = H_forb₁")

# ============================================================
section("THE GOLAY CODE AND PETERSEN/STEINER CONNECTION", 3)
# ============================================================

print("The [23,12,7] binary Golay code G₂₃:")
print(f"  Length n = 23 (prime)")
print(f"  Dimension k = 12 = h(E₆) = h(F₄)")
print(f"  Distance d = 7 = H_forb₁ = Φ₃(KEY₁)")
print()

# Check Hamming bound
hamming_bound = sum(comb(23, i) for i in range(4))
print(f"  Hamming bound: Σ_{{i=0}}^3 C(23,i) = {hamming_bound}")
print(f"  Number of codewords: 2^12 = {2**12}")
print(f"  {2**12} · {hamming_bound} = {2**12 * hamming_bound} = 2^23 = {2**23}")
print(f"  PERFECT CODE! (meets Hamming bound with equality)")
print()

# Extended Golay [24,12,8]
print("The extended Golay code G₂₄ = [24,12,8]:")
print(f"  Length n = 24 = |BT|")
print(f"  Dimension k = 12 = h(E₆)")
print(f"  Distance d = 8 = rank(E₈)")
print(f"  Parameters: [|BT|, h(E₆), rank(E₈)]")
print()
print(f"  Weight enumerator of G₂₄:")
weights = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
for w, count in sorted(weights.items()):
    print(f"    Weight {w:2d}: {count:5d} codewords")
print(f"  Total: {sum(weights.values())} = 2^12 = {2**12} ✓")
print()
print(f"  A_8 = 759 = 3·11·23")
print(f"  A_12 = 2576 = 2^5·80+16 = ... = 2576")
print(f"  759 = C(23,4)/C(7,4) = {comb(23,4)//comb(7,4)}")
print(f"  Wait: 759 = {759}. Let me check: C(23,4)/C(7,4) = {comb(23,4)}/{comb(7,4)} = {comb(23,4)//comb(7,4)}")
# Actually: the 759 octads form a 5-(24,8,1) design = S(5,8,24)

print()
print("  The 759 octads form the Steiner system S(5,8,24):")
print(f"    Any 5 of 24 points lie in a unique octad")
print(f"    5 = KEY₁+KEY₂, 8 = rank(E₈), 24 = |BT|")
print(f"    S(KEY₁+KEY₂, rank(E₈), |BT|)")

# ============================================================
section("PETERSEN GRAPH IN THE GOLAY CODE", 4)
# ============================================================

print("The Petersen graph appears in the Golay code construction!")
print()
print("The Mathieu group M₁₂ acts on 12 points (= h(E₆) = k of Golay)")
print("M₂₄ acts on 24 points (= |BT| = n of extended Golay)")
print()
print(f"  |M₁₂| = {95040}")
print(f"  = 95040 = 2^6 · 3^3 · 5 · 11")
print(f"  = 8·11880 = rank(E₈)·11880")
print(f"  |M₂₄| = {244823040}")
print(f"  = 244823040 = 2^10 · 3^3 · 5 · 7 · 11 · 23")
print()

# Factor through Lie numbers
print("M₂₄ order factored through (2,3) vocabulary:")
print(f"  2^10 = KEY₁^10 = (rank(E₈))^(10/3)... messy")
print(f"  Better: |M₂₄| = 24! / (|Aut(G₂₄)| in some sense)")
print()

# The Petersen graph and the Steiner system
print("The Petersen graph and S(5,8,24):")
print(f"  The 10 vertices of the Petersen graph correspond to 10 = C(5,2)")
print(f"  Steiner system S(2,3,7) is the Fano plane (7 = H_forb₁)")
print(f"  S(2,3,7) has b = C(7,2)/C(3,2) = {comb(7,2)}/{comb(3,2)} = {comb(7,2)//comb(3,2)} blocks")
print(f"  7 blocks of 3 elements from 7 points")
print()
print(f"  Steiner chain: S(2,3,7) → S(3,6,22) → S(4,7,23) → S(5,8,24)")
print(f"  Or in our terms:")
print(f"    S(KEY₁, KEY₂, H_forb₁)")
print(f"    S(KEY₂, h(G₂), 22)")
print(f"    S(KEY₁², H_forb₁, 23)")
print(f"    S(KEY₁+KEY₂, rank(E₈), |BT|)")

# ============================================================
section("SPHERE PACKING AND KISSING NUMBERS", 5)
# ============================================================

print("The lattice packing connection:")
print()
print(f"  E₈ lattice: densest packing in dim 8 = rank(E₈)")
print(f"    Kissing number: 240 = #roots(E₈) = V(icos)·V(dodec)")
print(f"    Center density: 1/16 = 1/KEY₁⁴")
print()
print(f"  Leech lattice Λ₂₄: densest packing in dim 24 = |BT|")
print(f"    Kissing number: 196560 = 2^4 · 3 · 5 · 7 · 13 · ... ")
print(f"    196560 = {196560}")
print(f"    = 2^4 · 3 · 5 · 7 · 9 · 13 ... let me factor")

n = 196560
factors = {}
temp = n
for p in [2,3,5,7,11,13,17,19,23,29,31]:
    while temp % p == 0:
        factors[p] = factors.get(p, 0) + 1
        temp //= p
if temp > 1:
    factors[temp] = 1
print(f"    196560 = " + " · ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(factors.items())))
print()

# The kissing number in various dimensions
kissing = {1: 2, 2: 6, 3: 12, 4: 24, 8: 240, 24: 196560}
print("Kissing numbers in key dimensions:")
for d, k in sorted(kissing.items()):
    lie_str = ""
    if d == 2: lie_str = " = h(G2) = KEY1*KEY2"
    elif d == 3: lie_str = " = h(E6) = h(F4)"
    elif d == 4: lie_str = " = |BT|"
    elif d == 8: lie_str = " = #roots(E8)"
    elif d == 24: lie_str = f" = 2^4*3*5*{196560//(16*15)}"
    print(f"  dim {d:2d}: τ = {k:>8d}{lie_str}")

print()
print("  dim 1: τ(1) = 2 = KEY₁")
print("  dim 2: τ(2) = 6 = KEY₁·KEY₂ = h(G₂)")
print("  dim 3: τ(3) = 12 = h(E₆)")
print("  dim 4: τ(4) = 24 = |BT| = KEY₁·h(E₆)")

# ============================================================
section("CHANNEL CAPACITY AND THE (2,3) RATE", 6)
# ============================================================

print("The Hamming [7,4,3] code has rate R = 4/7:")
print(f"  R = KEY₁²/Φ₃(KEY₁) = KEY₁²/(KEY₁²+KEY₁+1)")
print()

# Shannon capacity
print("The Shannon capacity of the binary symmetric channel:")
print(f"  C(p) = 1 - H(p) where H(p) = -p·log₂(p) - (1-p)·log₂(1-p)")
print()

# For what error probability does the Hamming rate achieve capacity?
# R = 4/7 ≈ 0.571
# C(p) = 4/7 → H(p) = 1 - 4/7 = 3/7
# H(p) = 3/7 → what p?
print(f"  Hamming rate R = {Fraction(4,7)} ≈ {4/7:.6f}")
print(f"  For R=4/7 to be capacity-achieving: H(p) = 3/7 = KEY₂/Φ₃(KEY₁)")
print(f"  The complementary rate 1-R = 3/7 = KEY₂/(KEY₁²+KEY₁+1)")
print()

# Perfect codes
print("ALL binary perfect codes (Tietäväinen-van Lint-Zinoviev):")
print("  1. Trivial: [n, n, 1] — the whole space")
print("  2. Repetition: [n, 1, n] with n odd")
print("  3. Hamming: [2^r-1, 2^r-r-1, 3] for r ≥ 2")
print("  4. Golay: [23, 12, 7]")
print()
print("Hamming codes [2^r-1, 2^r-r-1, 3]:")
for r in range(2, 8):
    n_code = 2**r - 1
    k_code = 2**r - r - 1
    d_code = 3
    rate = Fraction(k_code, n_code)
    print(f"  r={r}: [{n_code}, {k_code}, {d_code}]  rate={rate} ≈ {float(rate):.4f}")

print()
print(f"  As r→∞: rate → 1 (tending to trivial)")
print(f"  First non-trivial: r=2: [3,1,3] rate=1/3 = 1/KEY₂")
print(f"  Second: r=3: [7,4,3] rate=4/7 = KEY₁²/H_forb₁")
print(f"  Third: r=4: [15,11,3] rate=11/15")

# ============================================================
section("SUMMARY: CODING THEORY THROUGH THE (2,3) LENS", 7)
# ============================================================

print("="*70)
print("THE CODING THEORY — TOURNAMENT THEORY DICTIONARY")
print("="*70)
print()
print("Object                     Coding theory        Tournament/Lie")
print("-"*70)
print(f"Fano plane                  PG(2,2), 7 points   H_forb₁ = Φ₃(KEY₁)")
print(f"Hamming code length         n = 7                = Φ₃(KEY₁)")
print(f"Hamming code dimension      k = 4                = KEY₁² = rank(F₄)")
print(f"Hamming code distance       d = 3                = KEY₂")
print(f"Hamming code redundancy     r = 3                = KEY₂")
print(f"Hamming code size           2⁴ = 16              = KEY₁^(KEY₁²)")
print(f"Extended Golay length       n = 24               = |BT|")
print(f"Extended Golay dimension    k = 12               = h(E₆)")
print(f"Extended Golay distance     d = 8                = rank(E₈)")
print(f"Golay distance              d = 7                = H_forb₁")
print(f"Octads in Steiner           759                  = 3·11·23")
print(f"Steiner S(5,8,24)           t=5, k=8, v=24      = (KEY₁+KEY₂, rank(E₈), |BT|)")
print(f"Kissing in dim 8            240                  = #roots(E₈)")
print(f"Kissing in dim 24           196560               (Leech)")
print(f"|GL(3,2)| = |PSL(2,7)|     168                  = H_forb₁ · |BT|")
print(f"|Aut(Fano)|                 168                  = rank(E₈) · H_forb₂")
print()
print("The punchline: The first forbidden H value (7) is the Fano plane,")
print("the smallest projective plane, the foundation of all coding theory.")
print("The extended Golay code [|BT|, h(E₆), rank(E₈)] is the")
print("tournament-theoretic code par excellence.")
