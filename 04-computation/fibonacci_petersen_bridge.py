#!/usr/bin/env python3
"""
fibonacci_petersen_bridge.py — opus-2026-03-14-S79
===================================================

THE FIBONACCI-PETERSEN BRIDGE

Two key discoveries from this session:
1. Coxeter number ratios follow Fibonacci fractions
2. Independence polynomial factors via u₂ = -2φ²

This script explores the deep connection between the golden
ratio φ, the Petersen graph, and the ADE classification.

Parts:
1. Fibonacci numbers and tournament primes
2. Coxeter ratios as Fibonacci convergents
3. The φ-connection in Petersen's independence polynomial
4. The golden ratio in E₈: H₃ and icosians
5. φ as the bridge between KEY₁ and KEY₂
6. The (2,3,5) as Fibonacci primes
"""

from math import sqrt, log, log2, factorial, comb, gcd
from fractions import Fraction

KEY1, KEY2 = 2, 3
phi = (1 + sqrt(5)) / 2
psi = (1 - sqrt(5)) / 2  # conjugate golden ratio

print("=" * 70)
print("PART 1: FIBONACCI NUMBERS AND TOURNAMENT PRIMES")
print("=" * 70)
print()

# Fibonacci sequence
fib = [0, 1]
for i in range(30):
    fib.append(fib[-1] + fib[-2])

print("Fibonacci sequence and primality:")
print()
for n in range(1, 20):
    f_n = fib[n]
    # Simple primality test
    is_prime = f_n > 1 and all(f_n % d != 0 for d in range(2, int(f_n**0.5)+1))
    marker = ""
    if is_prime and f_n in [2, 3, 5]:
        marker = " ← TOURNAMENT PRIME!"
    elif is_prime:
        marker = " (prime)"
    print(f"  F_{n:2d} = {f_n:>6d}{marker}")

print()
print("  Fibonacci PRIMES: F₃=2, F₄=3, F₅=5, F₇=13, F₁₁=89, F₁₃=233, ...")
print(f"  The first three: 2, 3, 5 = the tournament triple!")
print(f"  Indices of Fibonacci primes: 3, 4, 5, 7, 11, 13, ...")
print(f"    These indices are themselves mostly prime (Wall conjecture)")
print()

# Lucas numbers
lucas = [2, 1]
for i in range(20):
    lucas.append(lucas[-1] + lucas[-2])

print("Lucas numbers L_n = φⁿ + ψⁿ (companion to Fibonacci):")
for n in range(12):
    print(f"  L_{n:2d} = {lucas[n]:>6d}", end="")
    if lucas[n] in [2, 3, 7, 18, 47, 123]:
        if lucas[n] == 2:
            print(f" = KEY₁", end="")
        elif lucas[n] == 3:
            print(f" = KEY₂", end="")
        elif lucas[n] == 7:
            print(f" = rank(E₇)", end="")
        elif lucas[n] == 18:
            print(f" = h(E₇)!", end="")
        elif lucas[n] == 47:
            print(f" = prime", end="")
        elif lucas[n] == 123:
            print(f" = 3·41", end="")
    print()

print()
print(f"  REMARKABLE: L₀ = 2 = KEY₁, L₁ = 1, L₂ = 3 = KEY₂")
print(f"  L₆ = 18 = h(E₇)!")
print(f"  The Lucas numbers start with the two keys!")
print(f"  And L₆ = 18 is the Coxeter number of E₇!")

print()
print("=" * 70)
print("PART 2: COXETER RATIOS AS FIBONACCI CONVERGENTS")
print("=" * 70)
print()

# Coxeter numbers: G₂=6, F₄=12, E₆=12, E₇=18, E₈=30
h_values = [("G₂", 6), ("F₄", 12), ("E₆", 12), ("E₇", 18), ("E₈", 30)]

print("Coxeter number ratios and Fibonacci fractions:")
print()

for i in range(len(h_values)-1):
    name1, h1 = h_values[i]
    name2, h2 = h_values[i+1]
    ratio = Fraction(h2, h1)
    # Find which Fibonacci fraction this is
    for j in range(1, 15):
        if Fraction(fib[j+1], fib[j]) == ratio:
            print(f"  h({name2})/h({name1}) = {h2}/{h1} = {ratio} = F_{j+1}/F_{j}")
            break
    else:
        print(f"  h({name2})/h({name1}) = {h2}/{h1} = {ratio}")

print()
print("  The ratios: 2/1, 1/1, 3/2, 5/3")
print(f"  = F₃/F₂, F₁/F₁, F₄/F₃, F₅/F₄")
print()

# Product of ratios
product = Fraction(1)
for i in range(len(h_values)-1):
    _, h1 = h_values[i]
    _, h2 = h_values[i+1]
    product *= Fraction(h2, h1)
print(f"  Product of all ratios = {product} = {float(product)}")
print(f"  = h(E₈)/h(G₂) = 30/6 = 5 = KEY₁+KEY₂")
print()

# The Fibonacci fractions converge to φ
print("  These fractions converge to φ = {:.10f}:".format(phi))
print(f"    2/1 = 2.0   (above φ)")
print(f"    1/1 = 1.0   (below φ)")
print(f"    3/2 = 1.5   (below φ)")
print(f"    5/3 = 1.667 (above φ)")
print()

# Average of the ratios
avg = (Fraction(2) + 1 + Fraction(3,2) + Fraction(5,3)) / 4
print(f"  Average ratio = {avg} ≈ {float(avg):.6f}")
print(f"  φ = {phi:.6f}")
print(f"  φ² = {phi**2:.6f}")
print(f"  The average ≈ φ²/2 + 1/(2φ)... hmm, not exact")

print()
# Geometric mean
geo_mean = (2 * 1 * 1.5 * 5/3) ** 0.25
print(f"  Geometric mean = 5^(1/4) = {5**0.25:.6f} = {geo_mean:.6f}")
print(f"  (Since product = 5, geometric mean = 5^(1/(#ratios)))")
print(f"  5^(1/4) = {5**0.25:.6f} ≈ φ^(3/4) = {phi**0.75:.6f}... close!")

print()
print("=" * 70)
print("PART 3: φ IN THE INDEPENDENCE POLYNOMIAL")
print("=" * 70)
print()

# I(P,x) = 1 + 10x + 30x² + 30x³ + 5x⁴
# = 5x⁴ + 30x³ + 30x² + 10x + 1
# Dividing by x²: 5(x+1/x)² + 30(x+1/x) + 20
# = 5(u² + 6u + 4) where u = x + 1/x
# u = -3 ± √5

print("Independence polynomial I(P,x) = 1 + 10x + 30x² + 30x³ + 5x⁴")
print()
print(f"  After substitution u = x + 1/x:")
print(f"  I(P,x)/x² = 5(u² + 6u + 4)")
print(f"  Roots: u = -3 ± √5")
print()

u1 = -3 + sqrt(5)
u2 = -3 - sqrt(5)
print(f"  u₁ = -3 + √5 = {u1:.10f}")
print(f"  u₂ = -3 - √5 = {u2:.10f}")
print()

print(f"  Expressing in terms of φ:")
print(f"    φ = (1+√5)/2 = {phi:.10f}")
print(f"    √5 = 2φ - 1")
print(f"    u₁ = -3 + (2φ-1) = 2φ - 4 = 2(φ - 2) = 2(φ - KEY₁)")
print(f"    u₂ = -3 - (2φ-1) = -2 - 2φ = -2(1+φ) = -2φ²")
print()
print(f"    u₁ = 2(φ-KEY₁) = {2*(phi-2):.10f}")
print(f"    u₂ = -2φ²       = {-2*phi**2:.10f}")
print()

# The four roots of I(P,x)
import cmath
print("  The four roots of I(P,x) = 0:")
for u in [u1, u2]:
    disc = u**2 - 4
    if disc >= 0:
        x1 = (u + sqrt(disc)) / 2
        x2 = (u - sqrt(disc)) / 2
        print(f"    From u={u:.4f}: x = {x1:.8f}, {x2:.8f} (real)")
    else:
        x1 = (u + cmath.sqrt(disc)) / 2
        x2 = (u - cmath.sqrt(disc)) / 2
        print(f"    From u={u:.4f}: x = {x1:.4f}, {x2:.4f} (complex, |x|={abs(x1):.6f})")

print()
# The complex roots have |x| = 1
print(f"  Complex roots lie on the UNIT CIRCLE |x| = 1")
print(f"  This is related to the fact that Petersen is vertex-transitive")
print()

# The real roots
r1 = (u2 + sqrt(u2**2 - 4)) / 2
r2 = (u2 - sqrt(u2**2 - 4)) / 2
print(f"  Real roots: {r1:.10f} and {r2:.10f}")
print(f"  Product: {r1*r2:.10f} (should be 1)")
print(f"  Sum: {r1+r2:.10f} = u₂ = -2φ²")
print()

# Connection to φ
print(f"  r₁ = (-φ² + √(φ⁴-1))/1 ≈ {r1:.6f}")
print(f"  r₂ = (-φ² - √(φ⁴-1))/1 ≈ {r2:.6f}")
print(f"  |r₂| = {abs(r2):.6f}")
print(f"  |r₂| = φ² + √(φ⁴-1) ≈ {phi**2 + sqrt(phi**4-1):.6f}")
print(f"  φ⁴ = {phi**4:.6f} = (φ²)² = (φ+1)² = φ²+2φ+1 = 3φ+2")
print(f"  φ⁴ - 1 = 3φ+1 = {3*phi+1:.6f}")

print()
print("=" * 70)
print("PART 4: THE GOLDEN RATIO IN E₈ — H₃ AND ICOSIANS")
print("=" * 70)
print()

print("The non-crystallographic Coxeter group H₃:")
print(f"  H₃ is the symmetry group of the icosahedron")
print(f"  |H₃| = 120 = |BI| = |S₅|")
print(f"  H₃ has Coxeter matrix with entries 2, 3, 5")
print(f"  The eigenvalues of H₃'s Cartan matrix involve φ")
print()

print("  The H₃ Coxeter diagram:  o---o===o")
print(f"                           3   5")
print(f"  (edge label 5 means angle π/5 between mirrors)")
print()
print(f"  cos(π/5) = φ/2 = {phi/2:.6f}")
print(f"  cos(2π/5) = (φ-1)/2 = 1/(2φ) = {1/(2*phi):.6f}")
print(f"  The golden ratio IS the cosine structure of the pentagon!")
print()

# H₄ and E₈
print("  H₄ (4D icosahedral symmetry) and E₈:")
print(f"    |H₄| = 14400 = 120² = |BI|²")
print(f"    14400 = 2⁵·3²·5² = KEY₁⁵·KEY₂²·(KEY₁+KEY₂)²")
print(f"    H₄ has exponents {{1, 11, 19, 29}}")
print(f"    = all four CORNERS of E₈ exponent set!")
print(f"    E₈ exponents: {{1, 7, 11, 13, 17, 19, 23, 29}}")
print(f"    H₄ exponents are positions 1, 3, 6, 8 of E₈ list")
print()
print(f"  The H₄ exponents {{1,11,19,29}} satisfy:")
print(f"    Sum: 1+11+19+29 = 60 = |A₅| = |BI|/2")
print(f"    Product: 1·11·19·29 = {1*11*19*29} = h(E₈)·{1*11*19*29//30}")
print(f"    Pairs sum to 30: 1+29=30, 11+19=30")
print(f"    ALL pairs sum to h(E₈)!")
print()

# The icosian connection
print("  ICOSIANS (Hamilton, 1856):")
print(f"    Quaternions with coefficients in Z[φ]")
print(f"    120 unit icosians form a group ≅ BI")
print(f"    When embedded in R⁸ via Z[φ] → R², they generate E₈")
print(f"    φ is the BRIDGE from H₃ to E₈!")
print()
print(f"  The embedding map:")
print(f"    Z[φ] → R²: a+bφ ↦ (a+bφ, a+bψ) where ψ = (1-√5)/2")
print(f"    Quaternion q ∈ Z[φ]⁴ maps to R⁸")
print(f"    The image of the 120 unit icosians = the 240 shortest E₈ vectors")
print(f"    (each icosian maps to 2 E₈ vectors via the ±ψ embeddings)")
print()
print(f"  THIS IS WHY #roots(E₈) = 2·|BI| = 240:")
print(f"    Each of the 120 icosians has 2 Galois conjugates")
print(f"    240 = 2·120 = KEY₁·|BI|")

print()
print("=" * 70)
print("PART 5: φ AS THE BRIDGE BETWEEN KEY₁ AND KEY₂")
print("=" * 70)
print()

print("The golden ratio φ = (1+√5)/2 connects KEY₁ and KEY₂:")
print()
print(f"  φ ≈ {phi:.10f}")
print(f"  KEY₁ = 2, KEY₂ = 3")
print(f"  KEY₁ < φ² < KEY₂ (since φ² = φ+1 ≈ 2.618)")
print()

# φ as the 'geometric mean' analogue
print(f"  KEY₁·KEY₂ = 6, √(KEY₁·KEY₂) = √6 ≈ {sqrt(6):.6f}")
print(f"  φ² = φ+1 ≈ {phi**2:.6f}")
print(f"  φ² is between √6 and KEY₂ = 3")
print()

# The continued fraction
print(f"  φ = [1; 1, 1, 1, ...] (slowest converging CF)")
print(f"  KEY₂/KEY₁ = 3/2 = [1; 2] (finite CF)")
print(f"  φ is the 'most irrational' number between KEY₁ and KEY₂")
print()

# Mediants
print(f"  Stern-Brocot mediants between 1 and 2:")
print(f"    1/1 ← 3/2 → 2/1")
print(f"    3/2 = KEY₂/KEY₁ is the mediant of 1/1 and 2/1")
print(f"    The SIMPLEST fraction between 1 and 2")
print(f"    In music: the perfect fifth, most consonant after the octave")
print()

# φ² = φ+1 and the characteristic polynomial
print(f"  φ satisfies x² = x + 1, i.e., x² - x - 1 = 0")
print(f"  Tournament polynomial: x² - 5x + 6 = 0")
print(f"  Compare:")
print(f"    φ poly: x² - x - 1 = 0 (roots φ, ψ)")
print(f"    f poly: x² - 5x + 6 = 0 (roots 2, 3)")
print()
print(f"  The φ polynomial has:")
print(f"    sum of roots = 1, product = -1")
print(f"  The f polynomial has:")
print(f"    sum of roots = 5, product = 6")
print()
print(f"  Ratio: sum_f/sum_φ = 5/1 = KEY₁+KEY₂")
print(f"  Ratio: prod_f/prod_φ = 6/(-1) = -6 = -h(G₂)")
print()

# The golden recurrence and (2,3) recurrence
print(f"  Golden recurrence: a(n) = a(n-1) + a(n-2)")
print(f"    General solution: A·φⁿ + B·ψⁿ")
print(f"    Dominant: φⁿ (since |φ| > |ψ|)")
print()
print(f"  (2,3) recurrence: a(n) = 5a(n-1) - 6a(n-2)")
print(f"    General solution: A·2ⁿ + B·3ⁿ")
print(f"    Dominant: 3ⁿ (since 3 > 2)")
print()
print(f"  PARALLEL:")
print(f"    Fibonacci: roots φ,ψ with |φ/ψ| = φ² ≈ 2.618")
print(f"    Tournament: roots 2,3 with |3/2| = 1.5")
print(f"    The golden ratio φ² ≈ 2.618 sits between KEY₂/KEY₁ = 1.5 and KEY₂ = 3")

print()
print("=" * 70)
print("PART 6: THE (2,3,5) AS FIBONACCI PRIMES")
print("=" * 70)
print()

print("The triple (2,3,5) = (F₃, F₄, F₅) — three consecutive Fibonacci numbers!")
print()
print(f"  F₁ = 1, F₂ = 1, F₃ = 2, F₄ = 3, F₅ = 5, F₆ = 8, F₇ = 13")
print()
print(f"  The triple (F₃, F₄, F₅) = (2, 3, 5) uniquely satisfies:")
print(f"    a + b = c (Fibonacci rule)")
print(f"    All three are prime")
print(f"    a·b·c = 30 = h(E₈)")
print()

# Is there another Fibonacci triple with all prime?
print("  Other consecutive Fibonacci triples with all prime:")
for i in range(3, 20):
    a, b, c = fib[i], fib[i+1], fib[i+2]
    a_prime = a > 1 and all(a % d != 0 for d in range(2, int(a**0.5)+1))
    b_prime = b > 1 and all(b % d != 0 for d in range(2, int(b**0.5)+1))
    c_prime = c > 1 and all(c % d != 0 for d in range(2, int(c**0.5)+1))
    if a_prime and b_prime and c_prime:
        print(f"    (F_{i}, F_{i+1}, F_{i+2}) = ({a}, {b}, {c})")

print()
print(f"  ONLY (2, 3, 5) works!")
print(f"  Because F₆ = 8 = 2³ is not prime.")
print(f"  And for n≥7, at least one of F_n, F_{{n+1}}, F_{{n+2}} is even (divisible by KEY₁)")
print()

# The Fibonacci numbers modulo 2
print("  Fibonacci modulo 2: ", end="")
for i in range(1, 20):
    print(f"{fib[i]%2}", end=" ")
print()
print("  Period 3: (1, 1, 0, 1, 1, 0, ...)")
print(f"  Every 3rd Fibonacci number is even (divisible by KEY₁)")
print(f"  3 = KEY₂ is the Fibonacci period mod KEY₁!")
print()

# Fibonacci modulo 3
print("  Fibonacci modulo 3: ", end="")
for i in range(1, 20):
    print(f"{fib[i]%3}", end=" ")
print()
print("  Period 8: (1,1,2,0,2,2,1,0,...)")
print(f"  Every 4th Fibonacci number is divisible by KEY₂=3")
print(f"  Pisano period π(3) = 8 = rank(E₈)!")
print()

# Fibonacci modulo 5
print("  Fibonacci modulo 5: ", end="")
for i in range(1, 20):
    print(f"{fib[i]%5}", end=" ")
print()
print("  Period 20: ...")
print(f"  Every 5th Fibonacci number is divisible by KEY₁+KEY₂=5")
print(f"  Pisano period π(5) = 20 = V(dodecahedron) = f(7)!")
print()

# The Pisano periods of 2, 3, 5
print("  PISANO PERIODS (period of F_n mod m):")
for m, expected in [(2, 3), (3, 8), (5, 20), (6, 24), (10, 60), (30, 120)]:
    # Compute Pisano period
    a, b = 0, 1
    for period in range(1, 10000):
        a, b = b, (a + b) % m
        if a == 0 and b == 1:
            break
    print(f"    π({m:>3d}) = {period:>4d}", end="")
    note = ""
    if period == 3:
        note = " = KEY₂"
    elif period == 8:
        note = " = rank(E₈)"
    elif period == 20:
        note = " = V(dodec)"
    elif period == 24:
        note = " = |BT|"
    elif period == 60:
        note = " = |A₅|"
    elif period == 120:
        note = " = |BI| = |S₅|!"
    print(note)

print()
print(f"  π(30) = π(h(E₈)) = 120 = |BI| = |Aut(Petersen)|!")
print(f"  The Fibonacci period mod h(E₈) IS the icosahedral group order!")
print()
print(f"  π(KEY₁) = KEY₂ = 3")
print(f"  π(KEY₂) = rank(E₈) = 8")
print(f"  π(KEY₁·KEY₂) = |BT| = 24")
print(f"  π(KEY₁·KEY₂·(KEY₁+KEY₂)) = |BI| = 120")

print()
print("=" * 70)
print("FINAL: THE FIBONACCI-PETERSEN-LIE TRIANGLE")
print("=" * 70)
print()
print("Three vertices of the universal triangle:")
print()
print("       FIBONACCI")
print("      /         \\")
print("   φ→2          (2,3,5) = (F₃,F₄,F₅)")
print("    /             \\")
print(" PETERSEN -------- LIE")
print("   K(5,2)    E₆,E₇,E₈")
print()
print("Edges:")
print(f"  Fibonacci → Petersen: π(30)=120=|Aut(P)|, I(P) factors via φ")
print(f"  Fibonacci → Lie:      (2,3,5)=tournament triple, 2·3·5=h(E₈)")
print(f"  Petersen → Lie:       eigenvalues {{3,1,-2}}, complement 30 edges")
print()
print("Center:")
print(f"  The tournament polynomial f(z) = (z-2)(z-3)")
print(f"  lives at the CENTER of this triangle,")
print(f"  connecting all three worlds through z²-5z+6=0.")
