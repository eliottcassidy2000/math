#!/usr/bin/env python3
"""
five_six_to_seven_eight.py — opus-2026-03-14-S74
The (5,6) → (7,8) transition: sum/product → quadratic/cubic.

Key framing from user:
  5 = 2 + 3     (sum)
  6 = 2 * 3     (product)
  7 = 2² + 3    (quadratic + additive)
  8 = 2³        (cubic)

How do these relate as recurrences?
What is the structural meaning of (5,6) as a pair?
How does (5,6) generate (7,8)?
"""

import numpy as np
from math import sqrt, gcd, comb
from fractions import Fraction
from itertools import permutations, combinations
from collections import Counter
import random, time

def banner(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")

phi = (1 + sqrt(5)) / 2

# ─────────────────────────────────────────────────────────────────────
# PART 1: The (5,6) pair as sum and product
# ─────────────────────────────────────────────────────────────────────
banner("PART 1: (5,6) = (2+3, 2·3) — SUM AND PRODUCT")

print("For ANY two numbers a, b:")
print("  s = a + b  (sum)")
print("  p = a · b  (product)")
print("  a and b are roots of z² - sz + p = 0")
print()
print("At (a,b) = (2,3):")
print("  s = 5, p = 6")
print("  z² - 5z + 6 = 0")
print("  (z-2)(z-3) = 0  ✓")
print()
print("So 2 and 3 ARE THE ROOTS of z² - 5z + 6 = 0!")
print("The polynomial z²-5z+6 with coefficients 5 and 6 ENCODES the keys.")
print()
print("Compare to the Fibonacci polynomial z²-z-1:")
print("  z²-z-1: sum=1, product=-1, roots=φ,ψ")
print("  z²-5z+6: sum=5, product=6, roots=2,3")
print()
print("The 'tournament polynomial' z²-5z+6 relates to Fibonacci z²-z-1:")
print("  Replace z by z+2: (z+2)² - 5(z+2) + 6 = z² - z = z(z-1)")
print("  So shifting by 2 turns z²-5z+6 into z²-z = z(z-1)")
print("  Not directly Fibonacci, but close.")
print()

# The discriminant
disc = 5**2 - 4*6
print(f"Discriminant of z²-5z+6: Δ = 25-24 = {disc}")
print(f"√Δ = 1")
print(f"So the 'gap' between the keys is √Δ = 1.")
print(f"Compare: Fibonacci Δ = 5, Jacobsthal (z²-z-2) Δ = 9.")
print()

# The characteristic polynomial of the 2×2 companion matrix
print("Companion matrix of z²-5z+6:")
C = np.array([[5, -6], [1, 0]])
print(C)
eigs = np.linalg.eigvals(C)
print(f"Eigenvalues: {sorted(eigs, reverse=True)}")
print(f"= 3, 2  (the keys!) ✓")
print()

# The recurrence x(n) = 5x(n-1) - 6x(n-2) has roots 2,3
print("Recurrence: x(n) = 5x(n-1) - 6x(n-2)")
print("With x(0)=0, x(1)=1: closed form x(n) = (3^n - 2^n)/(3-2) = 3^n - 2^n")
print()
seq = [0, 1]
for i in range(2, 15):
    seq.append(5*seq[-1] - 6*seq[-2])
print("  Sequence: ", seq[:15])
print("  = 3^n - 2^n: ", [3**n - 2**n for n in range(15)])
print()
print("With x(0)=2, x(1)=5: closed form x(n) = 3^n + 2^n")
seq2 = [2, 5]
for i in range(2, 15):
    seq2.append(5*seq2[-1] - 6*seq2[-2])
print("  Sequence: ", seq2[:15])
print("  = 3^n + 2^n: ", [3**n + 2**n for n in range(15)])
print()
print("So the 'Lucas-type' for (2,3) gives: 2, 5, 13, 35, 97, ...")
print("And the 'Fibonacci-type' for (2,3) gives: 0, 1, 5, 19, 65, 211, ...")

# ─────────────────────────────────────────────────────────────────────
# PART 2: How (5,6) generates (7,8)
# ─────────────────────────────────────────────────────────────────────
banner("PART 2: (5,6) GENERATES (7,8)")

print("Starting from 5 and 6:")
print(f"  5 + 2 = 7   (adding 2)")
print(f"  6 + 2 = 8   (adding 2)")
print(f"  5 + 1 + 1 = 7")
print(f"  6 + 1 + 1 = 8")
print()
print("But more structurally:")
print(f"  7 = 5 + 2 = (2+3) + 2 = 2·2 + 3 = 2² + 3")
print(f"  8 = 6 + 2 = (2·3) + 2 = 2³")
print(f"  Also: 8 = 5 + 3 = (2+3) + 3")
print(f"  Also: 7 = 6 + 1 = 2·3 + 1")
print()
print("In the polynomial ring Z[x] at x=2:")
print(f"  5 = x + (1+x) = 2x + 1 = x² + 1")
print(f"  6 = x(1+x) = x² + x")
print(f"  7 = x³ - 1")
print(f"  8 = x³")
print()
print("KEY: 6 = x(1+x) = x² + x  and  7 = x² + x + 1 = 6 + 1")
print("So 7 = 2·3 + 1 (product of keys plus identity)")
print()
print("RECURRENCE CONNECTION:")
print("  The recurrence x(n) = 5x(n-1) - 6x(n-2) with x(0)=1, x(1)=2:")
seq3 = [1, 2]
for i in range(2, 10):
    seq3.append(5*seq3[-1] - 6*seq3[-2])
print(f"  Sequence: {seq3[:10]}")
print(f"  = 2^n: {[2**n for n in range(10)]}")
print(f"  So x(n) = 5x(n-1) - 6x(n-2) with x(0)=1, x(1)=2 gives POWERS OF 2!")
print()

seq4 = [1, 3]
for i in range(2, 10):
    seq4.append(5*seq4[-1] - 6*seq4[-2])
print(f"  With x(0)=1, x(1)=3: {seq4[:10]}")
print(f"  = 3^n: {[3**n for n in range(10)]}")
print(f"  So x(n) = 5x(n-1) - 6x(n-2) with x(0)=1, x(1)=3 gives POWERS OF 3!")
print()
print("THEOREM: The recurrence with coefficients (5,-6) generates both")
print("power sequences 2^n and 3^n, because 2 and 3 are its roots.")
print()
print("At n=3:")
print(f"  2^3 = 8 = 5·4 - 6·2 = 20-12")
print(f"  3^3 = 27 = 5·9 - 6·3 = 45-18")
print()
print("So (5,6) literally GENERATES 8=2³ via the recurrence.")
print("And 7 = 8-1 = 2³-1 (the Mersenne prime just below).")

# ─────────────────────────────────────────────────────────────────────
# PART 3: The (5,6) pair in tournament theory
# ─────────────────────────────────────────────────────────────────────
banner("PART 3: (5,6) IN TOURNAMENT THEORY")

print("In the OCF: H = I(CG, 2)")
print()
print("I(CG, x) at x=5 and x=6:")
print("  I(5) evaluates at the SUM of the keys")
print("  I(6) evaluates at the PRODUCT of the keys")
print()
print("Since z²-5z+6 = (z-2)(z-3):")
print("  I(CG, 5) = 1 + 5α₁ + 25α₂ + 125α₃ + ...")
print("  I(CG, 6) = 1 + 6α₁ + 36α₂ + 216α₃ + ...")
print()

# Compute I(5) and I(6) at n=5
def adj_matrix(n, idx):
    A = [[0]*n for _ in range(n)]
    bits = idx
    for i in range(n):
        for j in range(i+1, n):
            if bits & 1: A[i][j] = 1
            else: A[j][i] = 1
            bits >>= 1
    return A

def count_ham_paths_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size: continue
            for v in range(n):
                if not (S & (1 << v)): continue
                S_prev = S ^ (1 << v)
                for u in range(n):
                    if (S_prev & (1 << u)) and A[u][v]:
                        dp[(S, v)] = dp.get((S, v), 0) + dp.get((S_prev, u), 0)
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# At n=5: H = 1 + 2α₁, so I(x) = 1 + α₁·x
# I(5) = 1 + 5α₁, I(6) = 1 + 6α₁
# I(5) - I(6) = -α₁
# I(6) - H = 4α₁

# Verify and find relationship
print("At n=5 (α₂=0): I(x) = 1 + α₁·x")
print("  I(5) = 1 + 5α₁ = I(6) - α₁")
print("  I(6) = 1 + 6α₁ = I(5) + α₁")
print("  H = 1 + 2α₁")
print("  I(3) = 1 + 3α₁")
print()
print("Relations:")
print("  I(5) = 3H - 2      (= 3(1+2α₁)-2 = 1+6α₁... wait)")
# I(5) = 1 + 5α₁ = 1 + 5(H-1)/2 = 1 + (5H-5)/2 = (2+5H-5)/2 = (5H-3)/2
print(f"  I(5) = (5H - 3) / 2")
print(f"  I(6) = (6H - 4) / 2 = 3H - 2")
print(f"  I(6) = 3H - 2  ← simple!")
print()
print("  So I(product) = 3·I(2) - 2 = 3H - 2")
print("  And I(sum) = (5H-3)/2")
print()

# Verify
for idx in [0, 40, 341, 500, 800]:
    A = adj_matrix(5, idx)
    h = count_ham_paths_dp(A, 5)
    alpha1 = (h - 1) // 2
    i5 = 1 + 5*alpha1
    i6 = 1 + 6*alpha1
    i5_formula = (5*h - 3) // 2
    i6_formula = 3*h - 2
    print(f"  idx={idx:3d}: H={h:3d}, α₁={alpha1}, I(5)={i5}={i5_formula}, I(6)={i6}={i6_formula}")

# ─────────────────────────────────────────────────────────────────────
# PART 4: The Vandermonde at (5,6) = (sum, product)
# ─────────────────────────────────────────────────────────────────────
banner("PART 4: VANDERMONDE AT (5,6)")

print("From the S71d work: V(a,b) = ab(b-a) for degree-2 extraction.")
print(f"  V(2,3) = 2·3·1 = 6")
print(f"  V(5,6) = 5·6·1 = 30")
print(f"  V(2,5) = 2·5·3 = 30")
print(f"  V(2,6) = 2·6·4 = 48")
print(f"  V(3,5) = 3·5·2 = 30")
print(f"  V(3,6) = 3·6·3 = 54")
print()
print("V(2,3) = 6 = THE PRODUCT OF THE KEYS = the minimum.")
print("V(5,6) = 30 = V(2,5) = V(3,5) = 5·6 = 5·V(2,3)")
print()
print("The Vandermonde at the sum/product pair (5,6) equals")
print("5 times the Vandermonde at the original pair (2,3).")
print("This factor of 5 = the gap between sum and product levels.")

# ─────────────────────────────────────────────────────────────────────
# PART 5: 6 as the multiplicative bridge
# ─────────────────────────────────────────────────────────────────────
banner("PART 5: 6 = 2·3 — THE MULTIPLICATIVE BRIDGE")

print("6 plays the multiplicative role that 5 plays additively:")
print()
print("  5 = 2 + 3  (additive combination)")
print("  6 = 2 · 3  (multiplicative combination)")
print()
print("Key properties of 6:")
print("  6 = 3! (factorial)")
print("  6 = C(4,2) (binomial coefficient)")
print("  6 = T₃ (3rd triangular number)")
print("  6 is the smallest perfect number (1+2+3=6)")
print("  6 = V(2,3) (Vandermonde determinant)")
print()
print("In tournament theory:")
print("  6 = number of edges in K₄ (the 4-vertex complete graph)")
print("  C(4,2) = 6 → tournament on 4 vertices has 6 arcs")
print("  2⁶ = 64 = number of tournaments on 4 vertices")
print()
print("The (5,6) pair has a UNIQUE property:")
print("  5 · 6 = 30")
print("  5 + 6 = 11  ← THE DECIMAL PAIR NUMBER!")
print()
print("So 5 + 6 = 11 = L(5) = J(5), connecting back to the bridge!")

# ─────────────────────────────────────────────────────────────────────
# PART 6: From (5,6) to (7,8) via operations
# ─────────────────────────────────────────────────────────────────────
banner("PART 6: FROM (5,6) TO (7,8)")

print("The four operations from (2,3) to (7,8):")
print()
print("  2 + 3 = 5   →  5 + 2 = 7")
print("  2 · 3 = 6   →  6 + 2 = 8")
print()
print("Both jumps are +2 (adding the first key)!")
print()
print("  (5, 6)  +2→  (7, 8)")
print("  (sum, product) → (sum+2, product+2)")
print()
print("Or in terms of the keys:")
print("  7 = (2+3) + 2 = 2·2 + 3 = 2² + 3")
print("  8 = (2·3) + 2 = 2³")
print()
print("7 combines BOTH keys (2² + 3)")
print("8 is PURELY the first key cubed (2³)")
print()
print("ALTERNATIVE: (5,6) → (7,8) via (5+2, 6+2)")
print("  Both shifted by x=2 from the (sum,product) level")
print()
print("AND: 5 + 6 = 11,  7 + 8 = 15")
print("  11 = L(5) = J(5) (the bridge)")
print("  15 = 2⁴ - 1 = max H at n=5")
print()
print("  5 · 6 = 30")
print("  7 · 8 = 56 = C(8,3) = number of 3-element subsets in 8 vertices")
print("  Also: 56 = edges in K₈ = C(8,2)... no, C(8,2)=28.")
print("  56 = 7·8 = C(8,3) = number of potential 3-cycles at n=8!")
print("  At n=8, there are C(8,3) = 56 triples, each potentially a 3-cycle.")

# ─────────────────────────────────────────────────────────────────────
# PART 7: The recurrence hierarchy {5,6,7,8} 
# ─────────────────────────────────────────────────────────────────────
banner("PART 7: RECURRENCE HIERARCHY")

print("Second-order recurrence f(n) = f(n-1) + x·f(n-2):")
print()
print(f"{'x':>4} {'Δ=1+4x':>8} {'roots':>20} {'roots rational?':>16} {'field':>12}")
for x in range(1, 13):
    delta = 1 + 4*x
    sq = int(sqrt(delta))
    rational = (sq*sq == delta)
    r1 = (1 + sqrt(delta)) / 2
    r2 = (1 - sqrt(delta)) / 2
    if rational:
        field = "Q"
        roots = f"{int(r1)}, {int(r2)}"
    else:
        roots = f"({1}±√{delta})/2"
        # squarefree part
        sf = delta
        for p in [4, 9, 25, 49]:
            while sf % p == 0:
                sf //= p
        field = f"Q(√{sf})"
    print(f"{x:4d} {delta:8d} {roots:>20} {'YES' if rational else 'no':>16} {field:>12}")

print()
print("KEY PATTERN:")
print("  x=5: Δ=21 = 3·7,  field Q(√21)")
print("  x=6: Δ=25 = 5²,   field Q (RATIONAL! roots 3,-2)")
print("  x=7: Δ=29,         field Q(√29)")
print("  x=8: Δ=33 = 3·11,  field Q(√33)")
print()
print("x=6 is the THIRD rational point (after x=0 and x=2)!")
print("At x=6: roots are 3 and -2. So f(n) = A·3ⁿ + B·(-2)ⁿ.")
print("  J₆(n) = (3ⁿ - (-2)ⁿ) / 5")
print("  L₆(n) = 3ⁿ + (-2)ⁿ")
print()
print("The rational points x = k(k-1): 0, 2, 6, 12, 20, 30, ...")
print("  x=0: roots 1, 0   (trivial)")
print("  x=2: roots 2, -1  (Jacobsthal — THE tournament level)")
print("  x=6: roots 3, -2  (the next level — denominator 5!)")
print("  x=12: roots 4, -3")
print()
print("The gap between rational points: 2, 4, 6, 8, 10, ...")
print("  = 2k for k=1,2,3,4,5,...")
print("  So x=2 → x=6 is a gap of 4 = 2·2")
print("  And x=6 → x=12 is a gap of 6 = 2·3")
print()
print("x=6 sits EXACTLY between x=2 and x=12:")
print("  Gap 2→6 = 4 = 2²")
print("  Gap 6→12 = 6 = 2·3")
print("  So the 'multiplicative bridge' 6 is also a bridge in the rational point sequence!")

# ─────────────────────────────────────────────────────────────────────
# PART 8: The 6-spectrum: J₆ and the number 7
# ─────────────────────────────────────────────────────────────────────
banner("PART 8: THE 6-SPECTRUM AND 7")

print("J₆(n) = (3ⁿ - (-2)ⁿ) / 5:")
j6 = [(3**n - (-2)**n) // 5 for n in range(15)]
print(f"  {j6}")
print()
print("J₆(3) = (27 + 8) / 5 = 35 / 5 = 7  ← THE TRANSITION NUMBER!")
print("J₆(4) = (81 - 16) / 5 = 65 / 5 = 13")
print("J₆(5) = (243 + 32) / 5 = 275 / 5 = 55 = F(10)")
print()
print("L₆(n) = 3ⁿ + (-2)ⁿ:")
l6 = [3**n + (-2)**n for n in range(15)]
print(f"  {l6}")
print()
print("L₆(3) = 27 - 8 = 19")
print("L₆(4) = 81 + 16 = 97")
print()

# The (5,6) recurrence x(n) = 5x(n-1) - 6x(n-2)
print("The (5,6)-recurrence x(n) = 5x(n-1) - 6x(n-2):")
print("  With x(0)=0, x(1)=1 → x(n) = 3ⁿ-2ⁿ")
diff = [3**n - 2**n for n in range(12)]
print(f"  {diff}")
print()
print("x(3) = 3³-2³ = 27-8 = 19")
print("x(4) = 3⁴-2⁴ = 81-16 = 65")
print()
print("Now: 3ⁿ-2ⁿ = 5·J₆(n) when n odd, and = 5·J₆(n) always?")
for n in range(1, 10):
    dn = 3**n - 2**n
    j6n = (3**n - (-2)**n) // 5
    print(f"  n={n}: 3ⁿ-2ⁿ={dn}, J₆(n)={j6n}, ratio={dn/j6n if j6n else 'inf':.4f}, diff={dn - 5*j6n}")

print()
print("3ⁿ-2ⁿ ≠ 5·J₆(n) in general (because (-2)ⁿ ≠ -2ⁿ for even n).")
print("But: 3ⁿ-2ⁿ for POSITIVE 2 vs J₆ using NEGATIVE (-2).")
print()
print("The (5,-6) recurrence uses positive roots (2,3).")
print("The x=6 recurrence uses roots (3,-2).")
print("Same ABSOLUTE values, different SIGNS for the smaller root.")

# ─────────────────────────────────────────────────────────────────────
# PART 9: 7 = 2² + 3 and 8 = 2³ — the quadratic-cubic pair
# ─────────────────────────────────────────────────────────────────────
banner("PART 9: 7 = 2² + 3  AND  8 = 2³")

print("The user's key insight: 7 = 2·2 + 3 = 2² + 3")
print()
print("In the polynomial ring Z[x] at x=2:")
print(f"  7 = x² + (1+x) = x² + x + 1  (cyclotomic Φ₃!)")
print(f"  8 = x³")
print()
print("Φ₃(x) = x²+x+1 = (x³-1)/(x-1)")
print("So 7 = Φ₃(2) = (2³-1)/(2-1) = 7/1 = 7 ✓")
print()
print("THIS IS THE KEY: 7 is the VALUE of the 3rd cyclotomic polynomial at x=2!")
print("  Φ₁(2) = 2-1 = 1")
print("  Φ₂(2) = 2+1 = 3")
print("  Φ₃(2) = 2²+2+1 = 7")
print("  Φ₄(2) = 2²+1 = 5")
print("  Φ₅(2) = 2⁴+2³+2²+2+1 = 31")
print("  Φ₆(2) = 2²-2+1 = 3")
print()
print("So the hierarchy IS the cyclotomic evaluation at x=2!")
print()

# Compute cyclotomic polynomials evaluated at 2
def cyclotomic_at_2(n):
    """Φ_n(2) using the formula: Φ_n(x) = prod_{d|n} (x^d - 1)^{μ(n/d)}"""
    from functools import reduce
    # Mobius function
    def mobius(k):
        if k == 1: return 1
        n = k
        factors = set()
        d = 2
        while d*d <= n:
            while n % d == 0:
                factors.add(d)
                n //= d
            d += 1
        if n > 1:
            factors.add(n)
        # Check squarefree
        n2 = k
        for p in factors:
            if n2 % (p*p) == 0:
                return 0
        return (-1)**len(factors)
    
    # Φ_n(2) = prod_{d|n} (2^d - 1)^{μ(n/d)}
    result = Fraction(1)
    for d in range(1, n+1):
        if n % d == 0:
            mu = mobius(n // d)
            if mu == 1:
                result *= (2**d - 1)
            elif mu == -1:
                result /= (2**d - 1)
    return int(result)

print(f"{'n':>4} {'Φ_n(2)':>10} {'factored':>20} {'in hierarchy?':>15}")
for n in range(1, 20):
    val = cyclotomic_at_2(n)
    # Factor
    factors = []
    v = val
    for p in range(2, 100):
        while v % p == 0:
            factors.append(p)
            v //= p
    if v > 1: factors.append(v)
    fstr = '·'.join(str(f) for f in factors) if factors else '1'
    
    hier = ""
    if val in [1, 2, 3, 5, 7, 8, 10, 11]:
        hier = f"  ← {val}"
    print(f"{n:4d} {val:10d} {fstr:>20} {hier:>15}")

print()
print("THE CYCLOTOMIC DICTIONARY:")
print("  Φ₁(2) = 1  (identity)")
print("  Φ₂(2) = 3  (second key)")
print("  Φ₃(2) = 7  (transition)")
print("  Φ₄(2) = 5  (bridge = sum of keys)")
print("  Φ₅(2) = 31 (Mersenne prime)")
print("  Φ₆(2) = 3  (= Φ₂, composite cyclotomic)")
print()
print("  And 8 = 2³ = x³ is NOT cyclotomic — it's a pure power.")
print("  7 = Φ₃(2) = cyclotomic")
print("  8 = 2³ = exponential")
print("  This is why (7,8) is a transition: cyclotomic meets exponential.")

# ─────────────────────────────────────────────────────────────────────
# PART 10: The cyclotomic factorization of 2^n - 1
# ─────────────────────────────────────────────────────────────────────
banner("PART 10: 2^n - 1 = PRODUCT OF CYCLOTOMIC VALUES")

print("2^n - 1 = ∏_{d|n} Φ_d(2)")
print()
for n in range(1, 13):
    prod = 1
    divs = [d for d in range(1, n+1) if n % d == 0]
    factors = []
    for d in divs:
        val = cyclotomic_at_2(d)
        prod *= val
        factors.append(f"Φ_{d}={val}")
    mn = 2**n - 1
    match = "✓" if prod == mn else "✗"
    print(f"  2^{n:2d}-1 = {mn:6d} = {'·'.join(factors)} = {prod} {match}")

print()
print("At n=3: 2³-1 = 7 = Φ₁·Φ₃ = 1·7")
print("At n=6: 2⁶-1 = 63 = Φ₁·Φ₂·Φ₃·Φ₆ = 1·3·7·3")
print()
print("The 6th Mersenne number 2⁶-1 = 63 = 3·3·7 = 9·7")
print("Contains BOTH 3 and 7 as factors!")
print("This is because 6 = 2·3, and both Φ₂(2)=3 and Φ₃(2)=7 appear.")

# ─────────────────────────────────────────────────────────────────────
# PART 11: Synthesis — the complete (5,6)→(7,8) story
# ─────────────────────────────────────────────────────────────────────
banner("PART 11: SYNTHESIS")

print("""
THE COMPLETE PICTURE:

(2, 3) are roots of z² - 5z + 6 = 0.
  5 = sum,  6 = product.

The recurrence x(n) = 5x(n-1) - 6x(n-2) generates:
  2ⁿ (with initial (1,2))
  3ⁿ (with initial (1,3))
  3ⁿ-2ⁿ (with initial (0,1))  
  3ⁿ+2ⁿ (with initial (2,5))

At n=3: x(3) = 2³ = 8  and  3³ = 27
  The (5,6)-recurrence generates 8 in 3 steps.

Cyclotomic structure:
  7 = Φ₃(2) = (2³-1)/(2-1) = the 3rd cyclotomic at 2
  8 = 2³ = pure power
  So (7,8) = (cyclotomic, exponential) pair.

From (5,6) to (7,8):
  5 + 2 = 7    (sum + first key = cyclotomic)
  6 + 2 = 8    (product + first key = exponential)
  Both shift by +2.

The number 6 is special:
  6 = 2·3 (product of keys)
  6 = V(2,3) (Vandermonde determinant)
  6 = 3! (factorial of second key)
  6 = C(4,2) (edges in K₄)
  
  And at x=6 in the recurrence f(n)=f(n-1)+6f(n-2):
  Roots are 3, -2 (the keys with sign flip on 2)
  J₆(3) = 7 (generates the transition number!)
  J₆(5) = 55 = F(10) (generates the doubled Fibonacci)

5 + 6 = 11 = L(5) = J(5) (the bridge to the decimal pair!)
5 · 6 = 30 = V(5,6) (the Vandermonde at sum/product level)

Everything connects through the polynomial z² - 5z + 6.
""")

print("Done.")
