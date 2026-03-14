#!/usr/bin/env python3
"""
Oblong Fugacity Landscape — integer eigenvalue points of the independence polynomial.
opus-2026-03-14-S86

DISCOVERY: The transfer matrix T_x = [[1,x],[1,0]] has integer eigenvalues
exactly when x = k(k+1) for some non-negative integer k.

x = 0:  eigenvalues (1, 0)   — trivial
x = 2:  eigenvalues (2, -1)  — TOURNAMENT WORLD
x = 6:  eigenvalues (3, -2)  — ???
x = 12: eigenvalues (4, -3)  — ???
x = 20: eigenvalues (5, -4)  — ???

General: x = k(k+1), eigenvalues (k+1, -k).

At each integer fugacity, I(P_n, x) satisfies a "generalized Fibonacci"
recurrence a(n) = a(n-1) + k(k+1)·a(n-2) with integer roots.

This script explores what these higher fugacity worlds look like
and whether they have tournament-like interpretations.
"""

import math
from collections import Counter, defaultdict
from fractions import Fraction

def I_path(n, x):
    """I(P_n, x) where P_n has n vertices."""
    if n == 0: return 1
    if n == 1: return 1 + x
    a, b = 1, 1 + x
    for _ in range(n - 1):
        a, b = b, b + x * a
    return b

def I_cycle(n, x):
    """I(C_n, x) where C_n has n vertices, n >= 3."""
    if n < 3: return None
    return I_path(n-1, x) + x * I_path(n-3, x)

def gen_fib(n, x):
    """Generalized Fibonacci at fugacity x: a(0)=0, a(1)=1, a(n)=a(n-1)+x*a(n-2)"""
    if n == 0: return 0
    if n == 1: return 1
    a, b = 0, 1
    for _ in range(n - 1):
        a, b = b, b + x * a
    return b

# ============================================================
# PART 1: THE OBLONG FUGACITY LADDER
# ============================================================
print("=" * 70)
print("PART 1: OBLONG FUGACITY LADDER — x = k(k+1)")
print("=" * 70)

print("\n  k | x=k(k+1) | eigenvalues | gen_Fib sequence (first 10)")
for k in range(0, 8):
    x = k * (k + 1)
    seq = [gen_fib(n, x) for n in range(10)]
    print(f"  {k} | {x:5d}    | ({k+1:2d}, {-k:3d})   | {seq}")

# ============================================================
# PART 2: CLOSED FORMS AT EACH FUGACITY
# ============================================================
print("\n" + "=" * 70)
print("PART 2: CLOSED FORMS — gen_Fib_k(n) = ((k+1)^n - (-k)^n) / (2k+1)")
print("=" * 70)

# For x = k(k+1), eigenvalues are (k+1) and -k.
# gen_Fib_k(n) = ((k+1)^n - (-k)^n) / ((k+1) - (-k)) = ((k+1)^n - (-k)^n) / (2k+1)

for k in range(0, 6):
    x = k * (k + 1)
    print(f"\nk={k}, x={x}: gen_Fib(n) = (({k+1})^n - ({-k})^n) / {2*k+1}")
    for n in range(1, 10):
        actual = gen_fib(n, x)
        if 2*k + 1 != 0:
            formula = ((k+1)**n - (-k)**n) // (2*k+1)
        else:
            formula = n  # k=0 case: x=0, gen_fib = n? No, 0,1,1,1,...
        print(f"  n={n}: actual={actual:10d}, formula={formula:10d}, match={actual==formula}")

# ============================================================
# PART 3: I(C_n, x) AT INTEGER EIGENVALUE POINTS
# ============================================================
print("\n" + "=" * 70)
print("PART 3: CYCLE INDEPENDENCE — I(C_n, k(k+1)) = (k+1)^n + (-k)^n")
print("=" * 70)

# I(C_n, x) with eigenvalues α, β satisfies I(C_n, x) = α^n + β^n
# (generalized Lucas)

for k in range(0, 6):
    x = k * (k + 1)
    alpha, beta = k + 1, -k
    print(f"\nk={k}, x={x}: I(C_n,x) = ({alpha})^n + ({beta})^n")
    for n in range(3, 10):
        actual = I_cycle(n, x)
        formula = alpha**n + beta**n
        print(f"  n={n}: I(C_n,{x}) = {actual:10d}, formula = {formula:10d}, match = {actual==formula}")

# ============================================================
# PART 4: FORBIDDEN VALUES AT EACH FUGACITY
# ============================================================
print("\n" + "=" * 70)
print("PART 4: 'FORBIDDEN' CYCLE VALUES AT EACH FUGACITY")
print("=" * 70)

# At x=2 (k=1): I(C_3, 2) = 7 = FORBIDDEN because K_3 can't be Ω
# At x=6 (k=2): I(C_3, 6) = 3^3 + (-2)^3 = 27 - 8 = 19
# At x=12 (k=3): I(C_3, 12) = 4^3 + (-3)^3 = 64 - 27 = 37
# At x=20 (k=4): I(C_3, 20) = 5^3 + (-4)^3 = 125 - 64 = 61

print("\n'Forbidden C_3 value' at each fugacity:")
print("  k | x=k(k+1) | I(C_3,x) = (k+1)^3 + (-k)^3 = (2k+1)(k^2+k+1)")
for k in range(0, 10):
    x = k * (k + 1)
    ic3 = (k+1)**3 + (-k)**3
    factor = (2*k+1) * (k*k + k + 1)
    print(f"  {k} | {x:5d}    | {ic3:10d} = {2*k+1} x {k*k+k+1}")

# Note: (k+1)^3 - k^3 = 3k^2 + 3k + 1 = 3k(k+1) + 1 = 3x + 1
# And (k+1)^3 + (-k)^3 = (k+1-k)((k+1)^2 + k(k+1) + k^2)
# = 1 · (k^2+2k+1 + k^2+k + k^2) = 3k^2 + 3k + 1 = 3x + 1

print(f"\n  PATTERN: I(C_3, k(k+1)) = 3k(k+1) + 1 = 3x + 1")
print(f"  At x=2: I(C_3,2) = 7 = 3(2)+1")
print(f"  At x=6: I(C_3,6) = 19 = 3(6)+1")
print(f"  At x=12: I(C_3,12) = 37 = 3(12)+1")
print(f"  BEAUTIFUL: the 'forbidden' C_3 value at fugacity x is always 3x+1!")

# Verify
for x_test in range(0, 30):
    ic3 = I_cycle(3, x_test)
    if ic3 != 3*x_test + 1:
        print(f"  COUNTEREXAMPLE: I(C_3,{x_test}) = {ic3} != {3*x_test+1}")
        break
else:
    print(f"\n  VERIFIED: I(C_3, x) = 3x + 1 for ALL x = 0..29")
    print(f"  (This is just 1 + 3x = I({{triangle}}, x) since C_3 has 3 vertices,")
    print(f"   3 edges, and independence sets: {{}}, {{0}}, {{1}}, {{2}} giving 1 + 3x.)")
    print(f"  Wait, C_3 has NO independent sets of size >= 2 (it's a triangle)!")
    print(f"  So I(C_3, x) = 1 + 3x exactly. Not 3k^2+3k+1!")

# Let me recheck...
# C_3 = triangle (K_3). Independent sets: {}, {0}, {1}, {2}. That's all.
# I(K_3, x) = 1 + 3x.
# At x=2: I(K_3, 2) = 1 + 6 = 7. ✓
# But earlier we said I(C_3, x) = (k+1)^3 + (-k)^3 when x=k(k+1)?
# That would give I(C_3, 2) = 3^3 + (-1)^3 = 27-1 = 26 ≠ 7!
# ERROR! The cycle formula I(C_n, x) = α^n + β^n uses the transfer matrix,
# but C_3 = K_3 is special.

# Wait — C_3 as a GRAPH is the triangle, which is K_3.
# As a cycle graph, C_3 has 3 vertices and 3 edges forming a cycle.
# This IS K_3. So I(C_3, x) = 1 + 3x because α(K_3) = 1.

# But the cycle recurrence gives I(C_n, x) = I(P_{n-1}, x) + x*I(P_{n-3}, x)
# For n=3: I(C_3, x) = I(P_2, x) + x*I(P_0, x) = (1+2x) + x(1) = 1+3x ✓

# So where does α^n + β^n come in?
# I(C_n, x) = I(P_{n-1}, x) + x*I(P_{n-3}, x)
# For LARGE n, I(C_n, x) ~ α^n + β^n via the characteristic equation
# But for small n, this isn't exactly right.

# Let's check: I(C_n, x) vs α^n + β^n
print(f"\n  CORRECTION: I(C_n, x) vs alpha^n + beta^n:")
for k in [1]:  # x=2
    x = k * (k + 1)
    alpha, beta = k + 1, -k
    print(f"\n  x={x}: alpha={alpha}, beta={beta}")
    for n in range(3, 12):
        ic = I_cycle(n, x)
        ab = alpha**n + beta**n
        print(f"    I(C_{n},{x}) = {ic:8d}, alpha^n+beta^n = {ab:8d}, match = {ic==ab}")

# AH HA! They DON'T match for C_3! I(C_3,2) = 7 but 2^3+(-1)^3 = 7. They DO match!
# Wait... 2^3 = 8, (-1)^3 = -1, sum = 7. And I(C_3,2) = 1+3(2) = 7. ✓ They agree!
# Let me recheck my earlier error...

# Actually my closed form in the earlier script said I(C_k, 2) = 2^k + (-1)^k
# For k=3: 2^3 + (-1)^3 = 8-1 = 7 ✓
# For k=4: 2^4 + (-1)^4 = 16+1 = 17 ✓
# So the formula IS correct! My confusion was thinking (k+1)^3+(-k)^3 for k=1
# gives 2^3 + (-1)^3 = 7, which equals 3x+1 = 7. Both are right!

print(f"\n  RESOLUTION: Both formulas agree!")
print(f"  I(C_3, x) = 1 + 3x (from direct counting)")
print(f"  I(C_3, k(k+1)) = (k+1)^3 + (-k)^3 = 3k^2+3k+1 = 3k(k+1)+1 = 3x+1")
print(f"  These are THE SAME formula: 1 + 3x = 3x + 1 ✓")

# ============================================================
# PART 5: THE JACOBSTHAL-LUCAS FAMILY
# ============================================================
print("\n" + "=" * 70)
print("PART 5: GENERALIZED FIBONACCI/LUCAS AT EACH OBLONG FUGACITY")
print("=" * 70)

# At x=k(k+1):
# Gen Fibonacci: G_k(n) = ((k+1)^n - (-k)^n) / (2k+1)
# Gen Lucas: L_k(n) = (k+1)^n + (-k)^n

print("\nFibonacci family (paths):")
print("  n\\k |" + "".join(f" k={k:>4d}  " for k in range(6)))
for n in range(0, 10):
    print(f"  {n:3d} |", end="")
    for k in range(6):
        x = k * (k + 1)
        val = gen_fib(n, x)
        print(f" {val:6d} ", end="")
    print()

print("\nLucas family (cycles): L_k(n) = (k+1)^n + (-k)^n")
print("  n\\k |" + "".join(f" k={k:>4d}  " for k in range(6)))
for n in range(3, 10):
    print(f"  {n:3d} |", end="")
    for k in range(6):
        lk = (k+1)**n + (-k)**n
        print(f" {lk:6d} ", end="")
    print()

# ============================================================
# PART 6: WHAT IS x=6 THEORY? (k=2)
# ============================================================
print("\n" + "=" * 70)
print("PART 6: THE x=6 WORLD — EIGENVALUES (3, -2)")
print("=" * 70)

print("""
At x=6 (k=2):
  Eigenvalues: 3 and -2
  Gen_Fib: G_2(n) = (3^n - (-2)^n) / 5
  Gen_Lucas: L_2(n) = 3^n + (-2)^n

This is the "base-3 minus base-(-2)" sequence divided by 5.
G_2: 0, 1, 1, 7, 13, 55, 133, 463, 1261, ...
L_2: 2, 1, 13, 19, 97, 211, 793, 1981, ...

INTERPRETATION:
  If tournaments use binary orientations (2 choices per arc),
  what uses TERNARY structure (3 choices per edge)?
  Answer: SIGNED GRAPHS! Edge can be +, -, or 0.
  Or: COLORED TOURNAMENTS with 3 arc types.
  Or: HYPERGRAPHS with 3-uniform structure.

At x=6: each independent set element contributes 6 = 2×3 weight.
This is the product of binary and ternary worlds!
""")

# Sequence G_2
G2 = [gen_fib(n, 6) for n in range(15)]
print(f"G_2(n) = (3^n - (-2)^n) / 5: {G2}")

# Check if in OEIS? (3^n-(-2)^n)/5 is related to Jacobsthal of order 3?
# Actually, gen_fib at x=6 is sometimes called "Jacobsthal of order 2"
# or related to Narayana's cows? Let me check values.

# Is G_2 related to any well-known sequence?
# 0, 1, 1, 7, 13, 55, 133, 463, 1261, 3895
# Check OEIS: this is A015518 (or similar)
print(f"\nG_2(n) values: {G2[:12]}")
print(f"Differences: {[G2[i+1] - G2[i] for i in range(11)]}")
print(f"Ratios: {[G2[i+1]/G2[i] if G2[i] else 'inf' for i in range(11)]}")

# Check: G_2(n) mod small numbers
print(f"\nG_2(n) mod 6: {[g % 6 for g in G2]}")
print(f"G_2(n) mod 5: {[g % 5 for g in G2]}")
print(f"G_2(n) mod 7: {[g % 7 for g in G2]}")

# ============================================================
# PART 7: THE TRIBONACCI CONNECTION (x=3, NOT k(k+1))
# ============================================================
print("\n" + "=" * 70)
print("PART 7: TRIBONACCI AND THE PLASTIC NUMBER (x=3)")
print("=" * 70)

# x=3 is NOT an oblong number (3 ≠ k(k+1) for any integer k).
# At x=3, eigenvalues are (1±√13)/2 ≈ 2.303 and -1.303.
# These are NOT integers.

# But the plastic number ρ ≈ 1.3247 satisfies ρ³ = ρ + 1.
# This is related to the tribonacci recurrence, not the Fibonacci one.

# The 2/3 TILE RECURRENCE from earlier sessions:
# T(n) = T(n-2) + T(n-3) counts tilings with tiles of length 2 and 3.
# This gives the Padovan/plastic number sequence.

# Transfer matrix for T(n) = T(n-2) + T(n-3):
# State: (a(n), a(n-1), a(n-2))
# a(n+1) = a(n-1) + a(n-2) = 0*a(n) + 1*a(n-1) + 1*a(n-2)
# Matrix: [[0,1,1],[1,0,0],[0,1,0]]

import numpy as np
M23 = np.array([[0, 1, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
eigs = np.linalg.eigvals(M23)
print(f"\n2/3 tile transfer matrix eigenvalues: {eigs}")
print(f"Largest = {max(abs(eigs)):.6f} (plastic number = {1.324717957:.6f})")

# The plastic number is the REAL root of t³ - t - 1 = 0.
# At tournament fugacity x=2, the path recurrence has roots 2 and -1.
# The 2/3 tile recurrence has root ρ ≈ 1.3247.
# These are DIFFERENT recurrences!

# But there IS a connection: the Fibonacci word on {2,3} uses tiles of length 1
# (for "2") and length 1 (for "3"), but the INTERPRETATION is:
# - "2" means a binary choice (arc orientation)
# - "3" means a ternary cycle (directed 3-cycle)

# The 2/3 COMPOSITION NUMBER from fibonacci_23_composition.py:
# Number of compositions of n using parts 2 and 3:
# C(0)=1, C(1)=0, C(2)=1, C(3)=1, C(4)=1, C(5)=2, C(6)=2, C(7)=3, C(8)=4, ...
comp23 = [0] * 20
comp23[0] = 1
for n in range(2, 20):
    if n >= 2: comp23[n] += comp23[n-2]
    if n >= 3: comp23[n] += comp23[n-3]

print(f"\n2/3 compositions: {comp23[:15]}")
print(f"Ratios: {[comp23[n]/comp23[n-1] if comp23[n-1] > 0 else 'inf' for n in range(4, 15)]}")

# ============================================================
# PART 8: THE k-TOURNAMENT GENERALIZATION
# ============================================================
print("\n" + "=" * 70)
print("PART 8: k-TOURNAMENTS — WHAT ARE THEY?")
print("=" * 70)

print("""
At each oblong fugacity x = k(k+1), we get a "k-tournament" theory:

k=0 (x=0): Trivial — all I values are 1.
k=1 (x=2): TOURNAMENTS — binary arc orientations.
  Each arc has 2 directions. H = I(Omega, 2).
k=2 (x=6): 3-TOURNAMENTS — ternary edge types?
  Each edge has 6 = 2*3 "states" (e.g., oriented + colored).
  H_6 = I(Omega, 6). Each independent set member contributes 6.
k=3 (x=12): 4-TOURNAMENTS — quaternary edge types?
  Each edge has 12 = 3*4 "states".
  H_12 = I(Omega, 12).

PATTERN: k-tournament has edges with k(k+1) states.
  The number of states factorizes as k × (k+1).
  At k=1: 1×2 = 2 states (forward/backward)
  At k=2: 2×3 = 6 states (like 2 colors × 3 orientations)
  At k=3: 3×4 = 12 states

The eigenvalue pair (k+1, -k) encodes:
  - k+1: "forward" growth rate
  - k: "backward" decay rate
  - Sum = 1 (always! trace of T_x)
  - Product = -k(k+1) = -x (determinant)

CONJECTURE: k-tournaments correspond to oriented (k+1)-uniform
hypergraphs with k colors per edge, giving k(k+1) total states.
""")

# ============================================================
# PART 9: CROSS-FUGACITY IDENTITIES
# ============================================================
print("\n" + "=" * 70)
print("PART 9: CROSS-FUGACITY IDENTITIES")
print("=" * 70)

# Is there a relationship between gen_Fib at different fugacities?
# At x=k(k+1): G_k(n) = ((k+1)^n - (-k)^n) / (2k+1)
# Cross: G_k(n) * (2k+1) = (k+1)^n - (-k)^n

# For k=1: 3*J(n) = 2^n - (-1)^n (Jacobsthal identity)
# For k=2: 5*G_2(n) = 3^n - (-2)^n

# Relation between G_1 and G_2?
# G_1(n) = (2^n - (-1)^n)/3 = Jacobsthal
# G_2(n) = (3^n - (-2)^n)/5

# Is G_2(n) = f(G_1(n))?
# Not obviously. But there might be a convolution identity.

print("\nCross-fugacity comparison:")
print(f"  n | G_0=n | G_1=J(n) | G_2 | G_3 | G_4")
for n in range(10):
    g = [gen_fib(n, k*(k+1)) for k in range(5)]
    print(f"  {n} | {g[0]:5d} | {g[1]:8d} | {g[2]:8d} | {g[3]:10d} | {g[4]:12d}")

# Check: G_k(n) satisfies G_k(n) = G_k(n-1) + k(k+1)*G_k(n-2)
# What if we compose? G_2(G_1(n))?
print("\nNested evaluation G_2(G_1(n)):")
for n in range(1, 8):
    j = gen_fib(n, 2)
    if j < 20:
        g2j = gen_fib(j, 6)
        print(f"  G_1({n}) = J({n}) = {j}, G_2({j}) = {g2j}")

# ============================================================
# PART 10: UNIVERSALITY — ALL INTEGER EIGENVALUE SEQUENCES
# ============================================================
print("\n" + "=" * 70)
print("PART 10: THE GRAND FAMILY — ALL INTEGER-EIGENVALUE FIBONACCI SEQUENCES")
print("=" * 70)

print("""
THE COMPLETE CLASSIFICATION:

All sequences satisfying a(n) = a(n-1) + x*a(n-2) with integer roots
have x = k(k+1) for some integer k >= 0.

The corresponding sequences form a HIERARCHY:

  k=0: a(n) = 0, 1, 1, 1, 1, ...  (eventually constant)
  k=1: a(n) = Jacobsthal           (tournament counts at x=2)
  k=2: a(n) = (3^n-(-2)^n)/5       (x=6 generalization)
  k=3: a(n) = (4^n-(-3)^n)/7       (x=12)
  ...
  k→∞: a(n) ~ (k+1)^n / (2k+1)    (exponential growth)

LUCAS PARTNERS:
  k=0: L(n) = 2, 1, 1, 1, ...
  k=1: L(n) = 2^n + (-1)^n = Jacobsthal-Lucas
  k=2: L(n) = 3^n + (-2)^n
  k=3: L(n) = 4^n + (-3)^n

THE TOURNAMENT IS k=1 in this infinite family.
The "forbidden value" I(C_3, k(k+1)) = 3k(k+1)+1 = 3x+1 grows linearly.

DEEPER STRUCTURE: The Fibonacci number itself (k=0.5..., x=3/4...no)
is NOT an integer-eigenvalue point! φ is irrational.
The Fibonacci sequence is the "continuous" version;
Jacobsthal (k=1) is its "integer shadow" at the nearest oblong point.
""")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS — OBLONG FUGACITY LANDSCAPE")
print("=" * 70)
print("""
CROWN JEWELS:

1. INTEGER EIGENVALUE FUGACITIES: x = k(k+1) (oblong numbers)
   0, 2, 6, 12, 20, 30, 42, 56, ...
   Eigenvalues: (k+1, -k). Tournament = k=1.

2. FORBIDDEN C_3 VALUE: I(C_3, x) = 3x + 1 for ALL x.
   At x=2: I(C_3,2) = 7 = forbidden H.
   At x=6: I(C_3,6) = 19.
   The triangle's independence polynomial is ALWAYS 1+3x, trivially.

3. CLOSED FORMS: G_k(n) = ((k+1)^n - (-k)^n) / (2k+1)
   This generalizes Jacobsthal to all oblong fugacities.
   L_k(n) = (k+1)^n + (-k)^n generalizes Jacobsthal-Lucas.

4. THE k-TOURNAMENT INTERPRETATION: At fugacity k(k+1),
   each edge has k(k+1) states = k × (k+1).
   For k=1: 2 states (orientations). For k=2: 6 states.

5. FIBONACCI IS NOT AN INTEGER POINT: φ ≈ 1.618 lives between
   k=0 (x=0) and k=1 (x=2). The golden ratio is the irrational
   "seed" that integer eigenvalue points approximate.

6. THE 2/3 COMPOSITION RECURRENCE: T(n) = T(n-2) + T(n-3)
   with plastic number root ρ ≈ 1.325 is a DIFFERENT family
   (3-step recurrence, not 2-step). The (2,3) tile counting
   is NOT the same as the independence polynomial recurrence.
   They share the atoms {2,3} but combine them differently.
""")
