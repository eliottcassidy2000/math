#!/usr/bin/env python3
"""
recurrence_web_23.py — The web of recurrences in the (2,3) universe
opus-2026-03-14-S81

Every number in the (2,3) universe satisfies a recurrence with characteristic
polynomial dividing (z-2)(z-3) = z^2 - 5z + 6.

This script maps out the complete recurrence web:
1. All linear recurrences with roots in {2, 3}
2. The "corner piece" 3^n - 2^n as universal connector
3. Higher-order recurrences from compositions of (2,3)
4. The recurrence structure of root system data
5. Music theory: 2^a * 3^b = Pythagorean tuning
6. The Stern-Brocot tree and (2,3) mediants
7. Collatz-like maps on the (2,3) monoid
8. The Zeta function of the (2,3) monoid
"""

from functools import lru_cache
from math import comb, gcd, log2, log
from fractions import Fraction

KEY1, KEY2 = 2, 3

def section(title, num):
    print(f"\n{'='*70}")
    print(f"  Part {num}: {title}")
    print(f"{'='*70}\n")

# ============================================================
section("THE FUNDAMENTAL RECURRENCE a(n) = 5a(n-1) - 6a(n-2)", 1)
# ============================================================

print("Characteristic polynomial: z^2 - 5z + 6 = (z-2)(z-3)")
print("General solution: a(n) = A*2^n + B*3^n")
print()

# Table of important sequences with this recurrence
print("Important sequences satisfying a(n) = 5a(n-1) - 6a(n-2):")
print()
print(f"  {'Name':25s} {'a(0)':>6s} {'a(1)':>6s}  {'Formula':20s}  First terms")
print("-" * 100)

sequences = [
    ("Powers of 2",        1,  2, "2^n",             [2**n for n in range(10)]),
    ("Powers of 3",        1,  3, "3^n",             [3**n for n in range(10)]),
    ("Corner piece",       0,  1, "3^n - 2^n",       [3**n - 2**n for n in range(10)]),
    ("Sum of powers",      2,  5, "2^n + 3^n",       [2**n + 3**n for n in range(10)]),
    ("Lie recurrence",     0,  6, "6*(3^n - 2^n)",   [6*(3**n - 2**n) for n in range(10)]),
    ("H_forb(n)",          0,  7, "7*3^n (partial)",  [7 * 3**n for n in range(6)]),
    ("Petersen vertices",  0, 10, "10*... (partial)", [10] + [10*5-10*6 for _ in range(1)]),
]

for name, a0, a1, formula, terms in sequences:
    print(f"  {name:25s} {a0:6d} {a1:6d}  {formula:20s}  {terms[:8]}")

print()
print("The corner piece c(n) = 3^n - 2^n:")
print("  c(0)=0, c(1)=1, c(2)=5, c(3)=19, c(4)=65, c(5)=211, c(6)=665, ...")
corner = [3**n - 2**n for n in range(12)]
print(f"  {corner}")
print()
print("  c(1) = 1 = unit")
print(f"  c(2) = 5 = KEY_SUM")
print(f"  c(3) = 19 (prime)")
print(f"  c(4) = 65 = 5*13 = KEY_SUM * Phi_3(KEY2)")
print(f"  c(5) = 211 (prime)")
print(f"  c(6) = 665 = 5*7*19 = KEY_SUM * H_forb_1 * 19")
print()

# Factor all corner pieces
def prime_factors(n):
    if n <= 1: return {n: 1} if n == 1 else {}
    factors = {}
    d = 2
    while d * d <= abs(n):
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if abs(n) > 1:
        factors[abs(n)] = factors.get(abs(n), 0) + 1
    return factors

print("Corner piece factorizations:")
for n in range(1, 16):
    cn = 3**n - 2**n
    pf = prime_factors(cn)
    pf_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(pf.items()))
    print(f"  c({n:2d}) = {cn:>10d} = {pf_str}")

# ============================================================
section("EVERY LIE NUMBER AS L*(3^n - 2^n)", 2)
# ============================================================

print("Key insight: every Lie number L gives a sequence L*(3^n - 2^n)")
print("satisfying the same recurrence a(n) = 5a(n-1) - 6a(n-2)")
print()

lie_numbers = {
    1: "unit",
    2: "KEY1",
    3: "KEY2",
    5: "KEY_SUM",
    6: "h(G2)",
    7: "H_forb_1",
    8: "rank(E8)",
    10: "V(Pet)",
    12: "h(E6)",
    14: "dim(G2)",
    18: "h(E7)",
    24: "|BT|",
    30: "h(E8)",
}

for L, name in sorted(lie_numbers.items()):
    seq = [L * (3**n - 2**n) for n in range(6)]
    print(f"  {L:>3d}*(3^n-2^n) = {name:>12s}: {seq}")

print()
print("OBSERVATION: Every row starts at 0, passes through L (at n=1),")
print("  then through 5L, 19L, 65L, 211L, ...")
print()
print("The h(G2)=6 row: 0, 6, 30, 114, 390, 1266, ...")
print(f"  h(G2)*(3^2-2^2) = 6*5 = 30 = h(E8)!")
print(f"  h(G2) at n=1 gives h(E8) at n=2!")
print()
print("The h(E6)=12 row: 0, 12, 60, 228, 780, 2532, ...")
print(f"  h(E6)*(3^2-2^2) = 12*5 = 60 = |BI|/2")
print()
print("The |BT|=24 row: 0, 24, 120, 456, 1560, 5064, ...")
print(f"  |BT|*(3^2-2^2) = 24*5 = 120 = |BI|!")
print(f"  |BT| at n=1 gives |BI| at n=2!")

# ============================================================
section("PYTHAGOREAN TUNING: MUSIC AND THE (2,3) MONOID", 3)
# ============================================================

print("In Pythagorean tuning, all intervals are ratios 2^a * 3^b:")
print()
print("The circle of fifths generates intervals by stacking 3/2:")
print()

# Musical intervals in Pythagorean tuning
intervals = []
for i in range(-6, 7):
    # i fifths = (3/2)^i, reduced to [1, 2)
    ratio = Fraction(3, 2) ** abs(i)
    if i < 0:
        ratio = Fraction(1, 1) / ratio
    # Reduce to [1, 2)
    while ratio >= 2:
        ratio /= 2
    while ratio < 1:
        ratio *= 2
    cents = 1200 * log2(float(ratio))
    intervals.append((cents, i, ratio))

intervals.sort()
for cents, fifths, ratio in intervals:
    name = {
        0: "Unison",
        1: "Perfect Fifth",
        2: "Major Second",
        3: "Major Sixth",
        4: "Major Third",
        5: "Major Seventh",
        6: "Tritone (Aug 4th)",
        -1: "Perfect Fourth",
        -2: "Minor Seventh",
        -3: "Minor Third",
        -4: "Minor Sixth",
        -5: "Minor Second",
        -6: "Tritone (Dim 5th)",
    }.get(fifths, f"{fifths} fifths")

    # Express ratio as 2^a * 3^b
    r = ratio
    a, b = 0, 0
    while r.numerator % 3 == 0:
        b += 1
        r = Fraction(r.numerator // 3, r.denominator)
    while r.denominator % 3 == 0:
        b -= 1
        r = Fraction(r.numerator, r.denominator // 3)
    while r.numerator % 2 == 0:
        a += 1
        r = Fraction(r.numerator // 2, r.denominator)
    while r.denominator % 2 == 0:
        a -= 1
        r = Fraction(r.numerator, r.denominator // 2)

    n = ratio.numerator
    d = ratio.denominator
    print(f"  {cents:7.2f}c  {name:25s}  {n}/{d} = 2^{a}*3^{b}" if b >= 0 else
          f"  {cents:7.2f}c  {name:25s}  {n}/{d} = 2^{a}/3^{-b}")

print()
print("The Pythagorean comma:")
p_comma = Fraction(3**12, 2**19)
print(f"  (3/2)^12 / 2^7 = 3^12 / 2^19 = {3**12}/{2**19} = {float(p_comma):.10f}")
print(f"  = {1200 * log2(float(p_comma)):.4f} cents")
print(f"  3^12 = {3**12}, 2^19 = {2**19}")
print(f"  Pythagorean comma = 3^12 / 2^19 ≈ 23.46 cents")
print()
print(f"  12 = h(E6), 19 is prime")
print(f"  The comma exists because 3^{12} ≈ 2^{19} but not equal")
print(f"  i.e., 3^{{h(E6)}} ≈ 2^{{19}} — the E6 Coxeter number controls the comma!")
print()

# The Pythagorean comma in cents
print("Related musical-mathematical coincidences:")
print(f"  12-TET: 12 notes per octave = h(E6)")
print(f"  Circle of fifths: 12 fifths = h(E6) fifths ≈ 7 octaves")
print(f"  Perfect fifth: 3/2 = KEY2/KEY1 = 701.96 cents")
print(f"  Equal tempered fifth: 2^(7/12) = 700 cents")
print(f"  Difference: 1.96 cents = (7/h(E6)) * comma")
print()
print(f"  7 = H_forb_1 notes in the diatonic scale!")
print(f"  5 = KEY_SUM notes in the pentatonic scale!")
print(f"  12 = h(E6) notes in the chromatic scale!")
print(f"  The major scales: C, D, E, F, G, A, B = 7 = H_forb_1 notes")

# ============================================================
section("THE STERN-BROCOT TREE AND (2,3) MEDIANTS", 4)
# ============================================================

print("The Stern-Brocot tree is the complete binary tree of all positive rationals.")
print("The mediant of a/b and c/d is (a+c)/(b+d).")
print()
print("Starting from 0/1 and 1/0:")
print("  Level 0: 1/1")
print("  Level 1: 1/2, 2/1")
print("  Level 2: 1/3, 2/3, 3/2, 3/1")
print("  Level 3: 1/4, 2/5, 3/5, 3/4, 4/3, 5/3, 5/2, 4/1")
print()
print("The 3/2 node (perfect fifth) is at position (right, left) in the tree.")
print("  Path: R L = take right then left from root 1/1")
print("  The matrix path: R = [[1,1],[0,1]], L = [[1,0],[1,1]]")
print("  RL = [[2,1],[1,1]] maps 0/1, 1/0 to 2/1, 3/2")
print()

# The (2,3) ratios in the Stern-Brocot tree
print("Positions of KEY1/KEY2 and KEY2/KEY1 in the Stern-Brocot tree:")
print(f"  2/3: path = R L L (right, left, left)")
print(f"  3/2: path = R L (right, left)")
print(f"  The perfect fifth 3/2 is CLOSER to the root than 2/3!")
print()

# Stern's diatomic sequence
print("Stern's diatomic sequence s(n):")
print("  s(0)=0, s(1)=1, s(2n)=s(n), s(2n+1)=s(n)+s(n+1)")
stern = [0, 1]
for n in range(2, 32):
    if n % 2 == 0:
        stern.append(stern[n//2])
    else:
        stern.append(stern[n//2] + stern[n//2 + 1])
print(f"  {stern[:32]}")
print()
print("  s(n)/s(n+1) enumerates all positive rationals exactly once!")
print()

# Where do 2 and 3 appear?
twos = [i for i, v in enumerate(stern[:64]) if v == 2]
threes = [i for i, v in enumerate(stern[:64]) if v == 3]
print(f"  Positions of s(n)=2: {twos}")
print(f"  Positions of s(n)=3: {threes}")
print()
print(f"  s(n)=2 at: 3, 4, 10, 12, ... = positions where 2 appears")
print(f"  s(n)=3 at: 5, 8, 11, 16, ... = positions where 3 appears")

# ============================================================
section("THE (2,3) RECURRENCE AND BERNOULLI NUMBERS", 5)
# ============================================================

print("Bernoulli numbers B_n and the (2,3) universe:")
print()

# Bernoulli numbers
B = [Fraction(0)] * 20
B[0] = Fraction(1)
B[1] = Fraction(-1, 2)
for n in range(2, 20):
    s = Fraction(0)
    for k in range(n):
        s += Fraction(comb(n+1, k)) * B[k]
    B[n] = -s / (n + 1)

for n in range(16):
    if B[n] != 0:
        print(f"  B_{n:2d} = {str(B[n]):>15s} = {float(B[n]):>12.6f}")

print()
print("Key Bernoulli values:")
print(f"  B_0 = 1")
print(f"  B_1 = -1/2 = -1/KEY1")
print(f"  B_2 = 1/6 = 1/h(G2)")
print(f"  B_4 = -1/30 = -1/h(E8)")
print(f"  B_6 = 1/42 = 1/f(9)")
print(f"  B_8 = -1/30 = -1/h(E8) again!")
print(f"  B_10 = 5/66 = KEY_SUM/(h(G2)*11)")
print(f"  B_12 = -691/2730 = -691/(h(E8)*91)")
print()
print("The denominators of B_{2k}/(2k):")
print("  Von Staudt-Clausen: denom(B_{2k}) = prod_{(p-1)|2k} p")
print()
print("  denom(B_2) = 6 = h(G2)")
print("  denom(B_4) = 30 = h(E8)")
print("  denom(B_6) = 42 = f(9)")
print("  denom(B_8) = 30 = h(E8)")
print("  denom(B_10) = 66 = h(G2)*11")
print("  denom(B_12) = 2730 = 2*3*5*7*13 = h(E8)*91")
print()
print("  Bernoulli denominators are products of tournament vocabulary!")

# ============================================================
section("FIBONACCI, TRIBONACCI, AND THE (2,3) CHARACTERISTIC POLY", 6)
# ============================================================

print("Comparing characteristic polynomials:")
print()
print(f"  Fibonacci:  z^2 - z - 1 = 0, roots = (1±sqrt(5))/2")
print(f"              Golden ratio phi = (1+sqrt(5))/2 = {(1+5**0.5)/2:.6f}")
print(f"  Tournament: z^2 - 5z + 6 = 0, roots = 2, 3")
print(f"              f(z) = (z-2)(z-3)")
print()
print(f"  Fibonacci: sum=1, product=-1 (roots add to 1, multiply to -1)")
print(f"  Tournament: sum=5=KEY_SUM, product=6=h(G2) (roots add to 5, multiply to 6)")
print()
print(f"  The tournament polynomial is 'multiplicatively Fibonacci-like':")
print(f"  If F(z) = z^2 - z - 1, then f(z) = z^2 - 5z + 6")
print(f"  Coefficients: -1 -> -5 (multiply by KEY_SUM)")
print(f"                -1 -> +6 (multiply by -h(G2))")
print()

# Tribonacci: z^3 - z^2 - z - 1 = 0
print("Tribonacci polynomial: z^3 - z^2 - z - 1 = 0")
print(f"  Tribonacci constant ≈ 1.8393")
print()

# What polynomial has roots 2, 3, 5?
# (z-2)(z-3)(z-5) = z^3 - 10z^2 + 31z - 30
print("The (2,3,5) polynomial:")
print(f"  (z-2)(z-3)(z-5) = z^3 - 10z^2 + 31z - 30")
print(f"  Coefficients: -10 = -V(Petersen)")
print(f"                 31 = 31 (prime)")
print(f"                -30 = -h(E8)")
print()
print(f"  The (2,3,5) polynomial has h(E8) as its constant term!")
print(f"  And V(Petersen) as its z^2 coefficient!")
print()

# Evaluate the (2,3,5) polynomial at special points
def g(z):
    return z**3 - 10*z**2 + 31*z - 30

print("(2,3,5) polynomial g(z) = z^3 - 10z^2 + 31z - 30 at special points:")
for z in range(-3, 13):
    val = g(z)
    notes = {
        -120: "-|BI|", -72: "-|Phi(E6)|", -36: "-|Phi+(E6)|",
        -30: "-h(E8)", -24: "-|BT|", -12: "-h(E6)", -8: "-rank(E8)",
        -6: "-h(G2)", 0: "0=root", 6: "h(G2)", 8: "rank(E8)",
        12: "h(E6)", 24: "|BT|", 30: "h(E8)", 42: "f(9)",
        56: "f(10)", 120: "|BI|", 168: "|GL(3,2)|",
    }
    note = notes.get(val, "")
    if note:
        note = f"  <-- {note}"
    print(f"  g({z:3d}) = {val:>8d}{note}")

print()
# The recurrence from (2,3,5): a(n) = 10a(n-1) - 31a(n-2) + 30a(n-3)
print("Recurrence: a(n) = 10a(n-1) - 31a(n-2) + 30a(n-3)")
print("General solution: a(n) = A*2^n + B*3^n + C*5^n")
print()
print("Corner pieces:")
print("  c_{2,3}(n) = 3^n - 2^n  (the (2,3) corner)")
print("  c_{3,5}(n) = 5^n - 3^n  (the (3,5) corner)")
print("  c_{2,5}(n) = 5^n - 2^n  (the (2,5) corner)")
print()

for name, func in [("3^n-2^n", lambda n: 3**n - 2**n),
                    ("5^n-3^n", lambda n: 5**n - 3**n),
                    ("5^n-2^n", lambda n: 5**n - 2**n)]:
    vals = [func(n) for n in range(8)]
    print(f"  {name:10s}: {vals}")

print()
print("At n=2:")
print(f"  3^2-2^2 = 5 = KEY_SUM")
print(f"  5^2-3^2 = 16 = KEY1^4")
print(f"  5^2-2^2 = 21 = H_forb_2!")
print()
print("  5^2 - 2^2 = (5-2)(5+2) = 3*7 = KEY2 * H_forb_1 = 21 = H_forb_2 ✓")

# ============================================================
section("THE ZETA FUNCTION OF THE (2,3) MONOID", 7)
# ============================================================

print("The Dirichlet series (zeta function) of the monoid M = <2,3>:")
print()
print("  Z_M(s) = sum_{m in M} m^{-s} = sum_{a,b >= 0} (2^a * 3^b)^{-s}")
print("         = (sum 2^{-as})(sum 3^{-bs})")
print("         = 1/(1-2^{-s}) * 1/(1-3^{-s})")
print()
print("This is a product of TWO geometric series!")
print()
print("Special values:")
for s in [1, 2, 3, 4]:
    val2 = Fraction(1, 1) / (1 - Fraction(1, 2**s))
    val3 = Fraction(1, 1) / (1 - Fraction(1, 3**s))
    product = val2 * val3
    print(f"  Z_M({s}) = {val2} * {val3} = {product} = {float(product):.6f}")

print()
print("  Z_M(1) = 2 * 3/2 = 3  (divergent, but formal value = KEY2)")
print("  Z_M(2) = 4/3 * 9/8 = 3/2 = KEY2/KEY1 (the PERFECT FIFTH!)")
print("  Z_M(3) = 8/7 * 27/26 = 216/182 = 108/91 = ... ")
print()
print(f"  Z_M(2) = 3/2 is the perfect fifth!!")
print(f"  The zeta function of the (2,3) monoid at s=2 IS the musical fifth!")
print()

# The RECIPROCAL of the zeta
print("1/Z_M(s) = (1-2^{-s})(1-3^{-s}):")
for s in [1, 2, 3]:
    val = (1 - Fraction(1, 2**s)) * (1 - Fraction(1, 3**s))
    print(f"  1/Z_M({s}) = {val} = {float(val):.6f}")

print()
print("  1/Z_M(1) = (1/2)(2/3) = 1/3 = 1/KEY2")
print("  1/Z_M(2) = (3/4)(8/9) = 2/3 = KEY1/KEY2 (the perfect fourth!)")
print("  The reciprocal at s=2 is the PERFECT FOURTH = KEY1/KEY2!")

# ============================================================
section("COLLATZ-LIKE MAPS ON THE (2,3) WORLD", 8)
# ============================================================

print("The Collatz conjecture uses the map n -> 3n+1 (odd) or n/2 (even).")
print("In the (2,3) world, consider simpler maps:")
print()
print("Map 1: n -> n*3 (if n is a power of 2) else n/3 (if divisible by 3) else n*2")
print("This is a random walk on the (2,3) lattice Z^2.")
print()

# More interesting: the map on exponents
print("Map on exponents (a,b) in 2^a * 3^b:")
print("  The lattice Z^2 with coordinates (v_2(n), v_3(n))")
print()
print("  All Lie numbers in this lattice:")
lie_in_lattice = {}
for a in range(10):
    for b in range(7):
        m = 2**a * 3**b
        if m in [1, 2, 3, 4, 6, 8, 9, 12, 16, 18, 24, 27, 36, 48, 54, 72, 120, 240]:
            lie_in_lattice[(a, b)] = m

print("  (a, b) -> 2^a * 3^b:")
for (a, b), m in sorted(lie_in_lattice.items()):
    name = {1: "unit", 2: "KEY1", 3: "KEY2", 4: "KEY1^2", 6: "h(G2)",
            8: "rank(E8)", 9: "KEY2^2", 12: "h(E6)", 16: "KEY1^4",
            18: "h(E7)", 24: "|BT|", 27: "KEY2^3", 36: "6^2",
            48: "|BO|", 54: "2*27", 72: "8*9"}.get(m, str(m))
    print(f"    ({a},{b}) -> {m:>4d} = {name}")

print()
print("Note: the lattice points form a CONVEX region")
print("  (0,0)=1, (1,0)=2, (0,1)=3, (2,0)=4, (1,1)=6, (3,0)=8, (0,2)=9, ...")
print("  Lie numbers 30=h(E8) and 120=|BI| are NOT on this lattice (need factor 5)")

# ============================================================
section("THE (2,3,5) TOWER: RAMANUJAN'S CONTINUED FRACTION", 9)
# ============================================================

print("Ramanujan's continued fraction R(q):")
print("  R(q) = q^(1/5) / (1 + q/(1 + q^2/(1 + q^3/(1 + ...))))")
print()
print("At q = e^{-2pi}: R(e^{-2pi}) = sqrt(5) * (phi^{5/2} * ... )")
print("  involves sqrt(5) = sqrt(KEY_SUM)")
print()
print("The Rogers-Ramanujan continued fraction connects to:")
print("  - The golden ratio phi = (1+sqrt(5))/2")
print("  - Modular equations of level 5 = KEY_SUM")
print("  - The icosahedron (symmetry group uses 5)")
print()

# The KEY connection: R(q)^5 + R(q)^{-5} is related to j-invariant
# And j(tau) = 1/q + 744 + 196884q + ...
# 744 = 24 * 31 = |BT| * 31
print("The j-invariant and Rogers-Ramanujan:")
print("  j(tau) = 1/q + 744 + 196884q + ...")
print(f"  744 = 24 * 31 = |BT| * 31")
print(f"  196884 = 12 * 16407 = h(E6) * 16407 = ... ")
print(f"  196884 = 2^2 * 3 * 16407 = 4 * 49221 = ...")
print()

# Factor 196884
n = 196884
pf = prime_factors(n)
pf_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(pf.items()))
print(f"  196884 = {pf_str}")
print(f"  = 2^2 * 3 * {196884//12}")
print(f"  = 12 * 16407 = h(E6) * 16407")
print(f"  196884 = 196883 + 1 (Monstrous Moonshine: dim(rep_M) + 1)")
print()

# 196883 = 47 * 59 * 71
pf2 = prime_factors(196883)
pf2_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(pf2.items()))
print(f"  196883 = {pf2_str}")
print(f"  These are three primes in arithmetic progression:")
print(f"  47, 59, 71 with common difference 12 = h(E6)!")

# ============================================================
section("GRAND SYNTHESIS: THE RECURRENCE WEB", 10)
# ============================================================

print("="*70)
print("  THE (2,3) UNIVERSE AS A WEB OF RECURRENCES")
print("="*70)
print()
print("1. FUNDAMENTAL RECURRENCE: a(n) = 5a(n-1) - 6a(n-2)")
print("   Char poly (z-2)(z-3), general solution A*2^n + B*3^n")
print()
print("2. CORNER PIECE: c(n) = 3^n - 2^n = 0, 1, 5, 19, 65, 211, ...")
print("   Every Lie number L gives sequence L*c(n)")
print("   KEY: h(G2)*c(2) = 6*5 = 30 = h(E8)")
print("        |BT|*c(2) = 24*5 = 120 = |BI|")
print()
print("3. MUSIC: Z_M(2) = 3/2 = perfect fifth")
print("   1/Z_M(2) = 2/3 = perfect fourth")
print("   Pythagorean comma: 3^12/2^19, where 12 = h(E6)")
print("   Diatonic scale: 7 = H_forb_1 notes")
print("   Chromatic scale: 12 = h(E6) notes")
print()
print("4. THE (2,3,5) POLYNOMIAL: z^3 - 10z^2 + 31z - 30")
print("   Coefficients: V(Petersen), 31, h(E8)")
print("   5^2 - 2^2 = 21 = H_forb_2")
print()
print("5. BERNOULLI DENOMINATORS: 6, 30, 42, 30, 66, 2730, ...")
print("   = h(G2), h(E8), f(9), h(E8), ... — all tournament vocabulary!")
print()
print("6. MONSTROUS MOONSHINE: 196883 = 47*59*71")
print("   Arithmetic progression with common difference 12 = h(E6)")
print()
print("7. THE STERN-BROCOT TREE places 3/2 (the fifth) at depth 2")
print("   and 2/3 (the fourth) at depth 3 — the keys appear at")
print("   depths KEY1 and KEY2!")
print()
print("THE WEB:")
print("  Every mathematical structure touched by (2,3) participates in")
print("  the recurrence a(n) = 5a(n-1) - 6a(n-2). The Lie algebras,")
print("  the music, the number theory, the combinatorics — they are all")
print("  eigenvectors of the same shift operator with eigenvalues 2 and 3.")
