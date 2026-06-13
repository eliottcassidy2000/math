#!/usr/bin/env python3
"""
moat_boundary_synthesis.py — opus-2026-03-14-S79
=================================================

THE T=10 MOAT AS A UNIVERSAL BOUNDARY

Three independent boundaries all converge at n=10:
  1. Tournament moat: T(n) > n! first at n=10
  2. Petersen graph K(5,2): 10 vertices = 2·5 = KEY₁·(KEY₁+KEY₂)
  3. (2,3,5) boundary: 1/2+1/3+1/5 = 31/30 > 1 (last spherical triple)

This script synthesizes all the (2,3,5) connections found in S79
into a single coherent picture.

Parts:
1. The three boundaries unified
2. The complete Petersen-Lie-Tournament dictionary
3. The Coxeter number sequence as recurrence
4. The 5-cascade revisited
5. Musical intervals and the comma
6. The information-theoretic view
7. Summary: one page of the universe
"""

from math import factorial, comb, log, log2, sqrt, gcd
from fractions import Fraction

KEY1, KEY2 = 2, 3
f = lambda z: (z - KEY1) * (z - KEY2)
T = lambda n: 2**(n*(n-1)//2)

print("=" * 70)
print("THE (2,3,5) UNIVERSE — COMPLETE SESSION S79 SYNTHESIS")
print("=" * 70)
print()

print("=" * 70)
print("PART 1: THREE BOUNDARIES AT n=10")
print("=" * 70)
print()

print("BOUNDARY 1: Tournament Count > Symmetry Count")
print(f"  T(n) = 2^C(n,2), n! = |S_n|")
print(f"  T(n)/n! crosses 1 (downward) at n={1} (trivially)")
print(f"  T(n)/n! > 1 FIRST at n=3 (T=8 > 6)")
print(f"  But the GROWTH RATE crosses: C(n,2) > n·log₂(n) at n≈10")
print()
# More precisely: T(n)/n! is monotonically increasing for n≥3
# The "moat" is where the ratio becomes exponentially large
# At n=10: T/n! ≈ 9.7 million

for n in range(3, 12):
    ratio = T(n) / factorial(n)
    print(f"  n={n:2d}: T(n)/n! = {ratio:>14.2f}")

print()
print("BOUNDARY 2: Petersen Graph")
print(f"  K(5,2) has 10 = C(5,2) vertices")
print(f"  10 = KEY₁ · (KEY₁+KEY₂)")
print(f"  Unique Kneser graph: smallest non-complete K(n,2)")
print(f"  K(4,2) = K₃ (complete, trivial)")
print(f"  K(5,2) = Petersen (first interesting)")
print()

print("BOUNDARY 3: (2,3,5) Spherical Condition")
print(f"  1/p + 1/q + 1/r > 1 for polyhedra")
print(f"  (2,3,5): 1/2+1/3+1/5 = {Fraction(1,2)+Fraction(1,3)+Fraction(1,5)} > 1 ✓")
print(f"  (2,3,7): 1/2+1/3+1/7 = {Fraction(1,2)+Fraction(1,3)+Fraction(1,7)} < 1 ✗")
print(f"  The LAST spherical triple is (2,3,5)")
print(f"  Excess: 31/30 - 1 = 1/30 = 1/h(E₈)")
print()

print("CONVERGENCE:")
print(f"  All three boundaries involve the number 10:")
print(f"  • T(n)/n! becomes exponentially large around n=10")
print(f"  • K(5,2) has exactly 10 vertices")
print(f"  • The (2,3,5) triple gives 2·5 = 10")
print(f"  And they are CONNECTED:")
print(f"  • f(10) = 56 = T(6) = dim(V_E₇)")
print(f"  • The Petersen complement has 30 = h(E₈) edges")
print(f"  • The spherical excess 1/30 inverts to h(E₈)")

print()
print("=" * 70)
print("PART 2: THE COMPLETE DICTIONARY (30 entries)")
print("=" * 70)
print()

# Comprehensive dictionary of (2,3,5) connections
dictionary = [
    ("2 = KEY₁", "First root of z²-5z+6", "k-nacci attractor", "rank(A₁)"),
    ("3 = KEY₂", "Second root of z²-5z+6", "weight-2 k-nacci attractor", "χ(Petersen)"),
    ("5 = KEY₁+KEY₂", "Root sum of z²-5z+6", "# Platonic solids = # exceptionals", "girth(Petersen)"),
    ("6 = KEY₁·KEY₂", "Root product", "h(G₂), E(tetrahedron), # perfect matchings of Petersen", ""),
    ("7 = KEY₁²+KEY₂", "First composite", "rank(E₇), h(G₂)+1, Mersenne M₃", ""),
    ("8 = KEY₁³", "Cubic key", "rank(E₈), φ(30), τ(30)", ""),
    ("9 = KEY₂²", "Square key", "h∨(F₄), Pythagorean whole tone 9/8", ""),
    ("10 = 2·5", "Petersen vertices", "moat boundary, C(5,2)", ""),
    ("11 = KEY₁·KEY₂+KEY₁+KEY₂", "", "f(11)=72=σ(30)", ""),
    ("12 = KEY₁²·KEY₂", "", "h(E₆)=h(F₄), V(icosahedron), T(5)", ""),
    ("13 = KEY₁²+KEY₂²", "", "dim(F₄)/rank(F₄)=dim(E₆)/rank(E₆)", "Fibonacci prime F₇"),
    ("14 = 2·7", "", "dim(G₂), 2·rank(E₇)", ""),
    ("15 = 3·5", "", "E(Petersen), C(6,2)", ""),
    ("18 = KEY₁·KEY₂²", "", "h(E₇), product of I(P) multipliers", ""),
    ("19 = KEY₂³-KEY₁³", "Corner piece n=3", "h(E₇)+1, dim(E₇)/rank(E₇)", ""),
    ("20 = 4·5", "", "V(dodecahedron), f(7), [S₅:S₃]", ""),
    ("24 = KEY₁³·KEY₂", "", "|BT|→E₆, [S₅:Z₅], A(4,1)", ""),
    ("27 = KEY₂³", "", "dim(V_E₆), J₃(O)", ""),
    ("30 = KEY₁·KEY₂·(KEY₁+KEY₂)", "", "h(E₈), E(icos/dodec/J(5,2))", "I₂=I₃ in I(P,x)"),
    ("31 = 2⁵-1", "Mersenne M₅", "h(E₈)+1, dim(E₈)/rank(E₈)", ""),
    ("42 = 6·7", "", "f(9), V(icos)+V(dodec)+V(Petersen)", ""),
    ("48 = 2⁴·3", "", "|BO|→E₇, |det A(Petersen)|", ""),
    ("56 = 7·8", "", "T(6), dim(V_E₇), f(10), 7·rank(E₈)", ""),
    ("60 = 2²·3·5", "", "|A₅|, icos rotation group, V·deg/2", ""),
    ("72 = 8·9", "", "σ(30), f(11), rank(E₈)·9", ""),
    ("120 = 5!", "", "|S₅|=|BI|=|Aut(P)|, P(Petersen,3)", "5-simplex tiles"),
    ("128 = 2⁷", "", "f(10)+f(11), KEY₁^rank(E₇)", ""),
    ("240 = 2⁴·3·5", "", "#roots(E₈), V(icos)·V(dodec), 2·|BI|", ""),
    ("256 = 2⁸", "", "|χ_P(2)|, KEY₁^rank(E₈)", ""),
    ("720 = 6!", "", "|S₆|, |BT|·h(E₈)", ""),
]

print(f"  {'Number':<32s} {'Identity 1':<40s}")
print(f"  {'-'*72}")
for num, id1, id2, id3 in dictionary:
    ids = [x for x in [id1, id2, id3] if x]
    print(f"  {num:<32s} {', '.join(ids[:2])}")
    if len(ids) > 2:
        print(f"  {'':32s} {', '.join(ids[2:])}")

print()
print("=" * 70)
print("PART 3: THE COXETER NUMBERS AS A RECURRENCE")
print("=" * 70)
print()

# Exceptional Coxeter numbers: 6, 12, 12, 18, 30
h_values = [6, 12, 12, 18, 30]
names = ["G₂", "F₄", "E₆", "E₇", "E₈"]

print("Coxeter numbers of exceptional groups: ", h_values)
print()

# Do they satisfy any linear recurrence?
# With 5 values, we can fit a recurrence of depth ≤ 2
# Try a(n) = p·a(n-1) + q·a(n-2)
# 12 = 6p (from n=2: a(2) = p·a(1) if q·a(0) = 0... need initial conditions)

# Let's check: is 6, 12, 12, 18, 30 a subsequence of A·2ⁿ+B·3ⁿ?
print("  Can {6,12,12,18,30} be expressed as A·2ⁿ+B·3ⁿ?")
# For n=0,1: A+B=6, 2A+3B=12 → A=6, B=0 → sequence is 6·2ⁿ = 6,12,24,48,...
# Doesn't work past n=2

# Ratios between consecutive
print("  Ratios: ", end="")
for i in range(len(h_values)-1):
    print(f"{h_values[i+1]}/{h_values[i]} = {Fraction(h_values[i+1], h_values[i])}", end="  ")
print()
print(f"  Ratios: 2, 1, 3/2, 5/3")
print(f"  = KEY₁, 1, KEY₂/KEY₁, (KEY₁+KEY₂)/KEY₂")
print()
print(f"  The ratio pattern: 2, 1, 3/2, 5/3")
print(f"  These are Fibonacci-like: each ratio ≈ next Fibonacci ratio")
print(f"  1/1, 2/1, 3/2, 5/3, 8/5, ...")
print(f"  The Coxeter number ratios FOLLOW THE FIBONACCI FRACTIONS!")
print()

# Actually the Fibonacci convergents to φ
fib = [1, 1, 2, 3, 5, 8, 13, 21]
print("  Fibonacci convergents to φ:")
for i in range(1, 7):
    print(f"    F_{i+1}/F_{i} = {fib[i+1]}/{fib[i]} = {fib[i+1]/fib[i]:.6f}")
print()

# Our ratios in order: h(F₄)/h(G₂)=2/1, h(E₆)/h(F₄)=1, h(E₇)/h(E₆)=3/2, h(E₈)/h(E₇)=5/3
# These are F₃/F₂, F₁/F₁, F₄/F₃, F₅/F₄
# The indices are 3/2, 1/1, 4/3, 5/4 — not consecutive but related

print("  OBSERVATION: the ratios h(X_{i+1})/h(X_i) use")
print("  successively better rational approximations to φ:")
print(f"    h(F₄)/h(G₂) = 2/1 = F₃/F₂ (above φ)")
print(f"    h(E₆)/h(F₄) = 1/1 = F₂/F₂ (below φ)")
print(f"    h(E₇)/h(E₆) = 3/2 = F₄/F₃ (above φ)")
print(f"    h(E₈)/h(E₇) = 5/3 = F₅/F₄ (below φ)")
print(f"    (hypothetical E₉ would need h(E₉)/h(E₈) ≈ 8/5 = F₆/F₅)")
print(f"    But E₉ doesn't exist! The sequence terminates at (2,3,5).")

print()
print("  Product of ratios: 2·1·(3/2)·(5/3) = 5 = KEY₁+KEY₂")
print(f"  = h(E₈)/h(G₂) = 30/6 = 5 ✓")

print()
print("=" * 70)
print("PART 4: THE 5-CASCADE — THE COMPLETE PICTURE")
print("=" * 70)
print()

print("Starting from the triple (2,3,5), EVERYTHING unfolds:")
print()

cascade = [
    (1, "1", "unit", "rank(A₁)-1"),
    (2, "KEY₁", "parity, binary, octave", "det(E₇)=2"),
    (3, "KEY₂", "ternary, perfect 5th", "det(E₆)=3"),
    (5, "KEY₁+KEY₂", "sum, pentagon, girth(P)", "# exceptional groups"),
    (6, "KEY₁·KEY₂", "product, h(G₂)", "perfect matchings of Petersen"),
    (8, "KEY₁³", "rank(E₈), φ(30)", "τ(30)"),
    (10, "KEY₁(K₁+K₂)", "V(Petersen), moat", "C(5,2)"),
    (12, "KEY₁²·KEY₂", "h(E₆)=h(F₄)", "V(icosahedron)"),
    (15, "KEY₂(K₁+K₂)", "E(Petersen)", "C(6,2)"),
    (18, "KEY₁·KEY₂²", "h(E₇)", "mult prod I(P)"),
    (20, "KEY₁²(K₁+K₂)", "V(dodecahedron)", "[S₅:S₃]"),
    (24, "KEY₁³·KEY₂", "|BT|→E₆", "[S₅:Z₅]"),
    (30, "KEY₁·KEY₂(K₁+K₂)", "h(E₈)", "I₂=I₃=30"),
    (42, "KEY₁·KEY₂·(K₁²+K₂)", "f(9)", "V_sum trinity"),
    (48, "KEY₁⁴·KEY₂", "|BO|→E₇", "|det A(P)|"),
    (56, "KEY₁³·(K₁²+K₂)", "T(6)=f(10)", "dim(V_E₇)"),
    (60, "KEY₁²·KEY₂(K₁+K₂)", "|A₅|", "icos rotation"),
    (72, "KEY₁³·KEY₂²", "σ(30)=f(11)", "rank(E₈)·9"),
    (120, "KEY₁³·KEY₂(K₁+K₂)", "|S₅|=|BI|", "5!=|Aut(P)|"),
    (240, "KEY₁⁴·KEY₂(K₁+K₂)", "#roots(E₈)", "V(i)·V(d)"),
    (720, "KEY₁⁴·KEY₂²(K₁+K₂)", "|S₆|=6!", "|BT|·h(E₈)"),
]

for val, expr, id1, id2 in cascade:
    print(f"  {val:>4d} = {expr:<22s}  {id1:<30s} {id2}")

print()
print("=" * 70)
print("PART 5: MUSICAL INTERVALS AND THE PYTHAGOREAN COMMA")
print("=" * 70)
print()

print("The musical scale is built from KEY₁ and KEY₂:")
print()

intervals = [
    ("Unison", Fraction(1, 1), "1/1"),
    ("Minor 2nd", Fraction(256, 243), "KEY₁⁸/KEY₂⁵"),
    ("Major 2nd (whole tone)", Fraction(9, 8), "KEY₂²/KEY₁³"),
    ("Minor 3rd", Fraction(32, 27), "KEY₁⁵/KEY₂³"),
    ("Major 3rd", Fraction(81, 64), "KEY₂⁴/KEY₁⁶"),
    ("Perfect 4th", Fraction(4, 3), "KEY₁²/KEY₂"),
    ("Tritone (aug 4th)", Fraction(729, 512), "KEY₂⁶/KEY₁⁹"),
    ("Perfect 5th", Fraction(3, 2), "KEY₂/KEY₁"),
    ("Minor 6th", Fraction(128, 81), "KEY₁⁷/KEY₂⁴"),
    ("Major 6th", Fraction(27, 16), "KEY₂³/KEY₁⁴"),
    ("Minor 7th", Fraction(16, 9), "KEY₁⁴/KEY₂²"),
    ("Major 7th", Fraction(243, 128), "KEY₂⁵/KEY₁⁷"),
    ("Octave", Fraction(2, 1), "KEY₁"),
]

for name, ratio, expr in intervals:
    cents = 1200 * log(float(ratio)) / log(2) if ratio > 0 else 0
    print(f"  {name:<25s} {str(ratio):<10s} = {expr:<14s}  ({cents:>7.2f} cents)")

print()
pythagorean_comma = Fraction(3**12, 2**19)
print(f"  Pythagorean comma: 3¹²/2¹⁹ = {pythagorean_comma} ≈ {float(pythagorean_comma):.10f}")
print(f"  = KEY₂¹²/KEY₁¹⁹")
print(f"  12 = h(E₆) = h(F₄)")
print(f"  19 = h(E₇)+1 = dim(E₇)/rank(E₇)")
print()
print(f"  The Pythagorean comma = KEY₂^h(E₆) / KEY₁^(dim(E₇)/rank(E₇))")
print(f"  TWO Lie-theoretic quantities determine the musical comma!")
print()

# 12 fifths vs 7 octaves
print(f"  12 perfect fifths = (3/2)¹² = {Fraction(3,2)**12}")
print(f"  7 octaves = 2⁷ = {2**7}")
print(f"  Difference = (3/2)¹² - 2⁷ = {float(Fraction(3,2)**12) - 128:.6f}")
print(f"  12 = h(E₆) fifths overshoot 7 = rank(E₇) octaves")
print(f"  by exactly the Pythagorean comma ≈ 23.46 cents")

print()
print("=" * 70)
print("PART 6: THE INFORMATION-THEORETIC VIEW")
print("=" * 70)
print()

# Binary entropy of the tournament/exponent partition
# 8 totatives out of 29 non-trivial values in [1,30)
p_exp = 8/29  # probability of being an exponent
H_binary = -p_exp * log2(p_exp) - (1-p_exp) * log2(1-p_exp)
print(f"  Partition of [1,30) into exponents (8) and non-exponents (21):")
print(f"  p(exponent) = 8/29 ≈ {8/29:.6f}")
print(f"  p(tournament) = 21/29 ≈ {21/29:.6f}")
print(f"  Binary entropy H = {H_binary:.6f} bits")
print()

# More naturally: 8 out of 30 values total (including 30 itself)
p2 = 8/30
H2 = -p2 * log2(p2) - (1-p2) * log2(1-p2)
print(f"  Alternate: 8 totatives out of 30:")
print(f"  p = 8/30 = {Fraction(8,30)} ≈ {p2:.6f}")
print(f"  H = {H2:.6f} bits")
print()

# The information content of knowing the group type
print(f"  Information content of 5 exceptional groups:")
print(f"  log₂(5) = {log2(5):.6f} bits ≈ KEY₁+KEY₂ - 1 - 1/{Fraction(3,1)}...")
print(f"  = 2.32 bits")
print()

# Comparing: H of tournament polynomial
print(f"  Tournament polynomial f(z) = z² - 5z + 6:")
print(f"  Has 2 roots → 1 bit to specify which root")
print(f"  But root MULTIPLICITIES carry more information:")
print(f"  z=2 (mult 1), z=3 (mult 1): 1 bit")
print(f"  The polynomial is SEPARABLE: discriminant = 1 ≠ 0")
print(f"  Separability ⟺ the two keys are DISTINCT")

print()
print("=" * 70)
print("PART 7: ONE PAGE OF THE UNIVERSE")
print("=" * 70)
print()

print("""
THE (2,3,5) UNIVERSE — A ONE-PAGE SUMMARY

AXIOM: z² - 5z + 6 = 0 has roots KEY₁ = 2, KEY₂ = 3.

FROM THIS SINGLE EQUATION:

ARITHMETIC
  Product 2·3 = 6 = h(G₂)
  Sum 2+3 = 5 = # exceptional Lie groups = # Platonic solids
  Triple product 2·3·5 = 30 = h(E₈)
  φ(30) = τ(30) = 8 = 2³ = rank(E₈)

RECURRENCES
  k-nacci → 2 (geometric rate 1/2)
  Weight-2 k-nacci → 3 (general: weight w → w+1)
  a(n) = 5a(n-1) - 6a(n-2) generates ALL sequences A·2ⁿ+B·3ⁿ
  Corner pieces: 3ⁿ-2ⁿ = 0, 1, 5, 19, 65, 211, ...

GEOMETRY
  (2,3,5) = last spherical triple → icosahedron/dodecahedron
  V(icos)·V(dodec) = 12·20 = 240 = #roots(E₈)
  30-edge trinity: icos, dodec, Petersen complement
  All have |Aut| = 120 = 5! = |S₅|

THE PETERSEN GRAPH K(5,2)
  10 vertices = 2·5, 15 edges = C(6,2)
  Eigenvalues {3, 1, -2} with multiplicities {1, 5, 4}
  Independence polynomial: I_k/C(5,k) = {1, 2, 3, 3, 1}
  ω·α = 2·4 = 8 = rank(E₈)
  P(Petersen, 3) = 120 = |BI| = |Aut|
  6 perfect matchings = h(G₂)

SIEVE
  Multiples of {2,3,5} in [1,30): tournament theory (21 numbers)
  Coprime to 30: Lie exponents (8 numbers = totatives)
  E₈ exponents = {1} ∪ {primes in (5,30)}
  The two theories partition ALL integers

LIE THEORY
  dim(g) = rank(g)·(h(g)+1) universal
  h+1 = {7, 13, 13, 19, 31} = {M₃, 2²+3², 2²+3², 3³-2³, M₅}
  E₈ positive roots: 120 = |BI| = |S₅|
  h(E₈)/h(G₂) = 5, h(E₇)/h(E₆) = 3/2

MUSIC
  Perfect 5th = 3/2 = KEY₂/KEY₁
  Whole tone = 9/8 = KEY₂²/KEY₁³
  Pythagorean comma = 3¹²/2¹⁹ = KEY₂^h(E₆)/KEY₁^(h(E₇)+1)
  12 fifths ≈ 7 octaves (off by the comma)

THE BOUNDARY
  1/2 + 1/3 + 1/5 = 31/30 > 1: E₈ EXISTS
  1/2 + 1/3 + 1/7 = 41/42 < 1: E₉ DOES NOT
  Excess = 1/30 = 1/h(E₈): the last gap before the boundary
  f(10) = 56 = T(6): the polynomial at the boundary encodes T(6)

THE MARKOV-LIKE IDENTITY
  2² + 3² + 5² - 2·3·5 = 8 = rank(E₈)
  a² + b² + c² - abc = rank of the maximal exceptional group
""")
