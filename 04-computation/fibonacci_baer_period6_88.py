#!/usr/bin/env python3
"""
fibonacci_baer_period6_88.py — opus-2026-03-14-S88

THE FIBONACCI-BAER-PERIOD-6 TRIANGLE

Three threads to unify:
1. Fibonacci sequence decomposed into 2s and 3s
2. Period 6 structure (transfer matrix M(-1) has order 6 in SL(2,Z))
3. Baer subplanes: PG(2,F_2) ⊂ PG(2,F_4)

Key insight: F_4 = F_2[α]/(α²+α+1) and x²+x+1 = Φ₃(x),
the 3rd cyclotomic polynomial. The primitive 3rd roots of unity
are the elements of F_4 \ F_2.

And Φ₃(x) divides x⁶-1 (since 3|6), connecting to period 6!

The chain: PERIOD 6 → CYCLOTOMIC Φ₃ → FIELD EXTENSION F_4/F_2 → BAER
"""

from collections import Counter, defaultdict
import sys

# ══════════════════════════════════════════════════════════════════
# PART 1: Fibonacci Zeckendorf decomposition — the 2-3 composition
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("PART 1: FIBONACCI AS A SEQUENCE OF 2s AND 3s")
print("=" * 70)

# Every positive integer has a unique representation as a sum of
# non-consecutive Fibonacci numbers (Zeckendorf). But that's not what
# we want. We want to see the Fibonacci RECURRENCE as generating
# intervals of size 2 and 3.

# The Fibonacci word: f₁ = "b", f₂ = "a", f_n = f_{n-1} f_{n-2}
# This gives: a, b, ab, bab, abbab, bababbab, ...
# If we replace a→2, b→3:
# The "musical" Fibonacci word with "long" (3) and "short" (2) intervals

# Actually: the CONTINUED FRACTION of φ = [1;1,1,1,...] means
# the Fibonacci sequence can be seen as:
# F(n) = F(n-1) + F(n-2) = "one of size F(n-1) and one of size F(n-2)"

# The Stern-Brocot / Fibonacci substitution rule:
# L → LS, S → L
# Starting: L
# Step 1: LS
# Step 2: LSL
# Step 3: LSLLS
# Step 4: LSLLSLSL
# Step 5: LSLLSLSLLSLLS

word = "L"
for i in range(10):
    new = ""
    for c in word:
        if c == "L":
            new += "LS"
        else:  # S
            new += "L"
    word = new

# If L = 3, S = 2:
values = [3 if c == "L" else 2 for c in word]
print(f"Fibonacci word (L=3, S=2), first 60 terms:")
print(f"  {''.join(str(v) for v in values[:60])}")
print(f"  Length: {len(values)}")

# Running sums = Fibonacci-like sequence
partial_sums = [0]
for v in values:
    partial_sums.append(partial_sums[-1] + v)

print(f"\nPartial sums (first 20): {partial_sums[:21]}")

# These partial sums should be related to Fibonacci numbers
fib = [1, 2]
while fib[-1] < 1000:
    fib.append(fib[-1] + fib[-2])
print(f"Fibonacci numbers: {fib[:15]}")

# Count 2s and 3s in first k terms
for k in [5, 8, 13, 21, 34, 55]:
    if k > len(values):
        break
    sub = values[:k]
    n2 = sub.count(2)
    n3 = sub.count(3)
    total = sum(sub)
    print(f"  First {k:2d}: n₂={n2}, n₃={n3}, ratio 3s/2s={n3/n2 if n2 > 0 else '∞':.4f}, "
          f"sum={total}")

# The ratio of 3s to 2s should approach φ = (1+√5)/2
phi = (1 + 5**0.5) / 2
print(f"\n  φ = {phi:.6f}")
print(f"  Ratio of Ls to Ss → φ (golden ratio)")

# ══════════════════════════════════════════════════════════════════
# PART 2: Period 6 modular structure
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 2: PERIOD 6 — THE UNIVERSAL TOURNAMENT PERIOD")
print("=" * 70)

# M(-1) = [[1,-1],[1,0]] has order 6 in SL(2,Z)
# M(2) = [[1,2],[1,0]] has order 6 mod 3 (because 2≡-1 mod 3)
# Fibonacci mod m has Pisano period dividing 6 for m=2,3,4,8

print("Fibonacci mod m periods (Pisano period π(m)):")
for m in range(2, 25):
    fib_mod = [0, 1]
    for i in range(200):
        fib_mod.append((fib_mod[-1] + fib_mod[-2]) % m)
    # Find period
    period = None
    for p in range(1, 201):
        if fib_mod[p] == 0 and fib_mod[p+1] == 1:
            period = p
            break
    divides_6 = period is not None and 6 % period == 0
    print(f"  m={m:2d}: π(m)={period:4d}" +
          (f"  (divides 6!)" if divides_6 else "") +
          (f"  = 6" if period == 6 else ""))

# Jacobsthal mod m periods
print("\nJacobsthal mod m periods:")
for m in [3, 5, 7, 9, 11, 13]:
    if m % 2 == 0:
        continue
    j = [0, 1]
    for i in range(200):
        j.append((j[-1] + 2*j[-2]) % m)
    period = None
    for p in range(1, 201):
        if j[p] == 0 and j[p+1] == 1:
            period = p
            break
    print(f"  m={m:2d}: π_J(m)={period}" +
          (f"  = 6" if period == 6 else ""))

# ══════════════════════════════════════════════════════════════════
# PART 3: Cyclotomic polynomials and the Baer connection
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 3: CYCLOTOMIC Φ₃ AND F_4/F_2")
print("=" * 70)

print("""
THE CYCLOTOMIC CONNECTION:

  Φ₃(x) = x² + x + 1 (3rd cyclotomic polynomial)

  Key properties:
  1. Φ₃(x) = 0 defines the primitive cube roots of unity
  2. F_4 = F_2[x] / Φ₃(x) — the field extension defining PG(2,F_4)!
  3. Φ₃(2) = 4 + 2 + 1 = 7 = H_forb_1 = |PG(2,F_2)|
  4. Φ₃(4) = 16 + 4 + 1 = 21 = H_forb_2 = |PG(2,F_4)|
  5. Φ₃(8) = 64 + 8 + 1 = 73 = |PG(2,F_8)| (NOT forbidden)

  So: H_forb_k = Φ₃(2^k) = |PG(2, F_{2^k})|
  for k=1,2 (but not k=3)

  The period connection: the 6th roots of unity include
  the 3rd roots (since 3|6). And:
  x⁶ - 1 = (x³-1)(x³+1) = Φ₁Φ₃(x²+1)(x²-x+1)...
  In characteristic 2: x⁶+1 = (x³+1)² = (x+1)²(x²+x+1)²
  So Φ₃ divides x⁶-1, and the period-6 structure
  CONTAINS the Baer subplane structure!
""")

# Verify: Φ₃(2^k) for k=1,2,3,...
print("Φ₃(2^k) = (2^k)² + 2^k + 1:")
for k in range(1, 8):
    val = (1 << (2*k)) + (1 << k) + 1
    print(f"  k={k}: Φ₃(2^{k}) = {1<<(2*k)} + {1<<k} + 1 = {val}" +
          (f" = H_forb_1" if val == 7 else "") +
          (f" = H_forb_2" if val == 21 else ""))

# ══════════════════════════════════════════════════════════════════
# PART 4: The Fibonacci word on the Baer plane
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 4: FIBONACCI WORD ON THE PROJECTIVE LINE")
print("=" * 70)

# The projective line PG(1, F_q) has q+1 points.
# PG(1, F_2) = 3 points (the "triangle")
# PG(1, F_4) = 5 points
# PG(1, F_8) = 9 points

# The Fibonacci word substitution L→LS, S→L has a connection
# to the automorphism of PG(1, F_4):
# |PGL(2, F_4)| = (4²-1)(4²-4)/(...) = 60 = |A_5|
# And A_5 is the icosahedral group, connected to φ!

print("Projective line sizes:")
for k in range(1, 6):
    q = 1 << k
    print(f"  |PG(1, F_{q})| = {q+1}")

print(f"\n  |PGL(2, F_4)| = (16-1)(16-4)/(4-1) = 15×12/3 = 60 = |A₅|")
print(f"  A₅ = icosahedral group, connected to golden ratio φ!")
print(f"  The icosahedron has golden ratio proportions.")
print(f"  |PGL(2, F_2)| = (4-1)(4-2)/(2-1) = 3×2 = 6")
print(f"  6 = the period of M(-1) in SL(2,Z)!")

# The chain: 6 → 60 → ...
# 60/6 = 10 = C(5,2) = Petersen vertices = BIBD points!

print(f"\n  60/6 = 10 = C(5,2) = Petersen vertices = BIBD points")
print(f"  The index [PGL(2,F_4) : PGL(2,F_2)] = 10")
print(f"  There are 10 copies of PGL(2,F_2) inside PGL(2,F_4)")
print(f"  And there are 10 partition pairs of {{0,...,5}}!")

# ══════════════════════════════════════════════════════════════════
# PART 5: The 360 = |A₆| connection
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 5: 360 BAER SUBPLANES = |A₆|")
print("=" * 70)

# We found 360 Baer subplanes in PG(2, F_4)
# 360 = 6!/2 = |A_6|
# S_6 has the unique outer automorphism!

print(f"  360 Baer subplanes in PG(2, F_4)")
print(f"  360 = 6!/2 = |A₆| (alternating group on 6)")
print(f"  S₆ is the ONLY symmetric group with outer automorphism!")
print(f"  |Aut(PG(2,F_4))| = |PΓL(3,F_4)| = |PGL(3,F_4)| × |Gal(F_4/F_2)|")
print(f"  |PGL(3,F_4)| = |GL(3,F_4)|/(4-1) = 181440/3 = 60480")
print(f"  |Gal(F_4/F_2)| = 2 (Frobenius x→x²)")
print(f"  |PΓL(3,F_4)| = 60480 × 2 = 120960")
print(f"  Orbits of Baer subplanes: |PΓL(3,F_4)| acts transitively")
print(f"  Stabilizer: 120960 / 360 = {120960 // 360}")
print(f"  336 = |PΓL(3,F_2)| = 2 × |PGL(3,F_2)| = 2 × 168 = 336")

print(f"\n  BEAUTIFUL: the stabilizer of a Baer subplane is |PΓL(3,F_2)| = 336!")
print(f"  This is 2 × |GL(3,F_2)| = 2 × 168 = 2 × 8 × 21")
print(f"  = 2 × 8 × H_forb_2")

# ══════════════════════════════════════════════════════════════════
# PART 6: The Fibonacci word in the Baer decomposition
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 6: FIBONACCI ↔ BAER DECOMPOSITION")
print("=" * 70)

# The Fibonacci word decomposes into Long (3) and Short (2).
# The projective plane decomposes into Baer (7) and Exterior (14).
# 7 = Φ₃(2) and 14 = 2 × 7.
# The ratio 14/7 = 2 (generator of tournament theory).
# But in the Fibonacci word, the ratio L/S → φ.

# Where do these meet?
# At the level of the INDEPENDENCE POLYNOMIAL:
# I(CP(n), x) = (1+2x)^n (cocktail party)
# I(Petersen, x) = ? (the Petersen graph, related to our BIBD)

# Actually: let me compute I(Petersen, x)
print("Independence polynomial of the Petersen graph:")

# Petersen: 10 vertices, 15 edges, 3-regular, girth 5
# Independence number = 4
# I(P, x) = 1 + 10x + 30x² + 30x³ + 5x⁴

# Let me compute from scratch
from itertools import combinations

petersen_verts = list(combinations(range(5), 2))
petersen_edges = []
for i, v1 in enumerate(petersen_verts):
    for j, v2 in enumerate(petersen_verts):
        if j > i and set(v1).isdisjoint(set(v2)):
            petersen_edges.append((i, j))

adj = [[False]*10 for _ in range(10)]
for a, b in petersen_edges:
    adj[a][b] = True
    adj[b][a] = True

# Count independent sets by size
indep_by_size = Counter()
for mask in range(1 << 10):
    # Check independence
    verts = [i for i in range(10) if mask & (1 << i)]
    is_indep = True
    for i, j in combinations(verts, 2):
        if adj[i][j]:
            is_indep = False
            break
    if is_indep:
        indep_by_size[len(verts)] += 1

print(f"  I(Petersen, x) = ", end="")
terms = []
for k in sorted(indep_by_size.keys()):
    coeff = indep_by_size[k]
    if k == 0:
        terms.append(f"{coeff}")
    elif k == 1:
        terms.append(f"{coeff}x")
    else:
        terms.append(f"{coeff}x^{k}")
print(" + ".join(terms))

# Evaluate at key points
for x in [-1, 2, -1/4]:
    val = sum(indep_by_size[k] * x**k for k in indep_by_size)
    print(f"  I(Petersen, {x:5.2f}) = {val:.4f}")

# I(Petersen, 2):
val_2 = sum(indep_by_size[k] * 2**k for k in indep_by_size)
print(f"\n  I(Petersen, 2) = {int(val_2)}")
print(f"  Compare: H_forb_2 = 21, H_forb_1 = 7")
print(f"  I(Petersen, 2) = {int(val_2)} = ", end="")
v = int(val_2)
# Factor
factors = []
n = v
for p in [2,3,5,7,11,13]:
    while n % p == 0:
        factors.append(p)
        n //= p
if n > 1: factors.append(n)
print(" × ".join(map(str, factors)))

# ══════════════════════════════════════════════════════════════════
# PART 7: The Grand Unification
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PART 7: THE FIBONACCI-BAER-PERIOD-6 TRIANGLE")
print("=" * 70)

print("""
THREE THREADS, ONE STRUCTURE:

THREAD 1: FIBONACCI PERIOD-6
  M(-1) = [[1,-1],[1,0]] has order 6 in SL(2,Z)
  M(2) ≡ M(-1) mod 3 (the mod-3 unification)
  Fibonacci mod 8 has period 12 = 2 × 6
  Jacobsthal mod 3 has period 6
  The fundamental period of tournament parity is 6 = 2 × 3

THREAD 2: BAER SUBPLANES
  PG(2,F_2) ⊂ PG(2,F_4): Fano inside 21-point plane
  F_4 = F_2[x]/Φ₃(x) where Φ₃(x) = x²+x+1
  7 = Φ₃(2), 21 = Φ₃(4), 73 = Φ₃(8)
  360 Baer subplanes = |A₆| (alternating on 6)
  Stabilizer = |PΓL(3,F_2)| = 336 = 2 × 168 = 2 × 8 × 21

THREAD 3: FIBONACCI 2-3 DECOMPOSITION
  Fibonacci word: L→LS, S→L (L=3, S=2)
  Ratio L/S → φ (golden ratio)
  Partial sums approximate Fibonacci numbers
  Connected to continued fraction [1;1,1,1,...]

THE UNIFICATION:

  ┌─────────────────┐
  │   PERIOD 6      │
  │ = LCM(2, 3)     │
  │ = ord(M(-1))    │
  └────────┬────────┘
           │
     ┌─────┴──────┐
     │             │
  ┌──┴───┐   ┌────┴────┐
  │  2   │   │    3    │
  │ gen  │   │  cycle  │
  └──┬───┘   └────┬────┘
     │             │
     └──────┬──────┘
            │
  ┌─────────┴──────────┐
  │  Φ₃(x) = x²+x+1   │
  │  = cyclotomic of 3  │
  │  divides x⁶-1      │
  └─────────┬──────────┘
            │
     ┌──────┴───────┐
     │              │
  ┌──┴───┐   ┌─────┴─────┐
  │ F_4   │   │ PG(2,F_4) │
  │ =F_2² │   │ =21 pts   │
  │via Φ₃ │   │ =H_forb_2 │
  └──┬───┘   └─────┬─────┘
     │              │
     └──────┬───────┘
            │
  ┌─────────┴──────────┐
  │   BAER SUBPLANE    │
  │  7 ⊂ 21 (Fano ⊂ PG)│
  │  H_forb_1 ⊂ H_forb_2│
  └─────────┬──────────┘
            │
  ┌─────────┴──────────┐
  │ FIBONACCI WORD      │
  │ L=3, S=2            │
  │ L/S → φ             │
  │ φ = root of x²-x-1 │
  │   = Φ₃ twin:        │
  │   x²+x+1 vs x²-x-1 │
  └────────────────────┘

THE DEEPEST CONNECTION:
  Φ₃(x) = x² + x + 1 (defines F_4, Baer, forbidden values)
  Fibonacci: x² - x - 1 = 0 (defines φ, golden ratio)

  These differ by: signs of x and constant!
  Φ₃(x) at x=2: 7 (forbidden)
  x²-x-1 at x=2: 4-2-1 = 1 (trivial)

  But: Φ₃(-x) = x² - x + 1 (which IS close to x²-x-1)
  The difference: Φ₃(-x) - (x²-x-1) = 2
  THE GENERATOR 2 IS THE GAP BETWEEN BAER AND FIBONACCI!
""")
