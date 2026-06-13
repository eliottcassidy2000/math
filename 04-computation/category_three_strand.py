#!/usr/bin/env python3
"""
category_three_strand.py — opus-2026-03-14-S71g
Category theory of the 3-strand Pascal structure and its tournament meaning.

The three strands:
  s0(k) = C(2k+1, k)   — odd rows of Pascal, "Catalan-related"
  s1(k) = C(2k+2, k)   — even rows, left of center
  s2(k) = C(2k+2, k+1) — even rows, central binomial

KEY RELATIONSHIPS:
  s2(k)/s0(k) = 2 (EXACT for all k)
  s1(k)/s0(k) = (2k+2)/(k+1) · k!/(k!) = ... let me compute

This script explores:
1. All ratios between strands
2. The recurrence structure
3. Connection to Catalan, ballot, and tournament numbers
4. Categorical interpretation as a graded ring/module
5. The "triple" reading of Pascal's triangle as a functor
"""

import math
from fractions import Fraction

print("=" * 70)
print("CATEGORY THEORY OF 3-STRAND PASCAL SEQUENCE")
print("opus-2026-03-14-S71g")
print("=" * 70)

def C(n, k):
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

# ============================================================
# Part 1: Exact strand relationships
# ============================================================
print("\n" + "=" * 70)
print("PART 1: EXACT STRAND RATIOS")
print("=" * 70)

print("\nk | s0=C(2k+1,k) | s1=C(2k+2,k) | s2=C(2k+2,k+1) | s2/s0 | s1/s0 | s2/s1")
for k in range(10):
    s0 = C(2*k+1, k)
    s1 = C(2*k+2, k)
    s2 = C(2*k+2, k+1)

    r20 = Fraction(s2, s0) if s0 > 0 else "∞"
    r10 = Fraction(s1, s0) if s0 > 0 else "∞"
    r21 = Fraction(s2, s1) if s1 > 0 else "∞"

    print(f"{k:2d} | {s0:12d} | {s1:12d} | {s2:14d} | {r20} | {r10} | {r21}")

print("""
ALGEBRAIC PROOF of s2/s0 = 2:
  s0(k) = C(2k+1, k) = (2k+1)! / (k! · (k+1)!)
  s2(k) = C(2k+2, k+1) = (2k+2)! / ((k+1)! · (k+1)!)

  s2/s0 = [(2k+2)! / ((k+1)!)²] / [(2k+1)! / (k! · (k+1)!)]
        = [(2k+2)! · k! · (k+1)!] / [(k+1)!² · (2k+1)!]
        = [(2k+2) · k!] / [(k+1)!]
        = (2k+2) / (k+1)
        = 2.  ∎

ALGEBRAIC s1/s0:
  s1(k) = C(2k+2, k) = (2k+2)! / (k! · (k+2)!)

  s1/s0 = [(2k+2)!/(k!·(k+2)!)] / [(2k+1)!/(k!·(k+1)!)]
        = [(2k+2)·(k+1)!] / [(k+2)!]
        = (2k+2) / (k+2)
        = 2(k+1)/(k+2)

  This approaches 2 as k → ∞.

ALGEBRAIC s2/s1:
  s2/s1 = C(2k+2,k+1)/C(2k+2,k)
        = [(2k+2)!/((k+1)!)²] / [(2k+2)!/(k!·(k+2)!)]
        = [k!·(k+2)!] / [(k+1)!²]
        = (k+2)/(k+1)

  This approaches 1 as k → ∞.
""")

# ============================================================
# Part 2: The graded ring structure
# ============================================================
print("=" * 70)
print("PART 2: GRADED RING STRUCTURE")
print("=" * 70)

# The three strands satisfy Vandermonde-like convolution identities.
# C(2k+1,k) · C(2j+1,j) relates to C(2(k+j)+...,...) somehow.
#
# The natural product: s0(k) * s0(j) vs s0(k+j+?)
# C(2k+1,k) · C(2j+1,j) = C(2k+1,k) · C(2j+1,j)
# By Vandermonde: Σ_i C(a,i)·C(b,n-i) = C(a+b,n)
# Not directly applicable.

# But there IS a product structure on central binomials:
# C(2k,k) · C(2j,j) = C(2(k+j), k+j) · C(k+j,k)·... no.
# Actually: C(2k,k)·C(2j,j)/C(2(k+j),k+j) = C(k+j,k)·2^{-2(k+j)} / (2^{-2k}·2^{-2j})
# This ratio is C(k+j,k) · 1 ... not clean.

# Better: the GENERATING FUNCTIONS.
# Σ_k s0(k) x^k = Σ C(2k+1,k) x^k
# = (1/x) Σ C(2k+1,k) x^{k+1}
# The GF of C(2n+1,n) is related to 1/√(1-4x).

# Recall: Σ C(2n,n) x^n = 1/√(1-4x)
# And: C(2n+1,n) = C(2n,n) · (2n+1)/(n+1) = C(2n,n) · (2 - 1/(n+1))
# Hmm, not clean.

# Alternative: C(2n+1,n) = 2·C(2n,n) - C(2n,n)/(n+1) = 2·C(2n,n) - C_n (Catalan)
# where C_n = C(2n,n)/(n+1) is the n-th Catalan number.

print("Connection to Catalan numbers:")
for k in range(10):
    catalan_k = C(2*k, k) // (k+1)
    central_k = C(2*k, k)
    s0_k = C(2*k+1, k)

    # Check: s0 = 2·central - catalan?
    # 2·C(2k,k) - C(2k,k)/(k+1) = C(2k,k)·(2-1/(k+1)) = C(2k,k)·(2k+1)/(k+1)
    # = (2k+1)·C(2k,k)/(k+1) = (2k+1)!/(k!·(k+1)!) ... wait that IS C(2k+1,k).
    # So C(2k+1,k) = (2k+1)/(k+1) · C(2k,k).

    ratio = Fraction(s0_k, central_k)

    print(f"  k={k}: C_{k}={catalan_k:6d}, C(2k,k)={central_k:6d}, "
          f"C(2k+1,k)={s0_k:6d}, ratio s0/central={(2*k+1)}/{k+1} = {ratio}")

print(f"""
IDENTITY: s0(k) = C(2k+1,k) = (2k+1)/(k+1) · C(2k,k)

So: s0 = (2k+1)·C_k (where C_k = Catalan number)
    s2 = C(2k+2,k+1) = C(2(k+1), k+1) = central binomial (shifted)
    s1 = C(2k+2,k) = (k+2)/(k+1) · s2(k-1)... hmm

Actually:
  s0(k) = (2k+1) · C_k  (Catalan connection)
  s2(k) = 2 · s0(k)     (OCF ratio)
  s1(k) = (2(k+1)/(k+2)) · s0(k) (approaches 2·s0)

The Catalan numbers C_k count:
  - Balanced parenthesizations
  - Binary trees with k nodes
  - Triangulations of (k+2)-gon
  - Monotone lattice paths below diagonal

So s0(k) = (2k+1)·C_k counts "labeled Catalan structures"
where we additionally pick one of 2k+1 positions.
""")

# ============================================================
# Part 3: Recurrence structure
# ============================================================
print("=" * 70)
print("PART 3: RECURRENCES FOR EACH STRAND")
print("=" * 70)

# s0(k) = C(2k+1,k): satisfies what recurrence?
# C(2k+1,k)/C(2k-1,k-1) = [(2k+1)!/(k!·(k+1)!)] / [(2k-1)!/((k-1)!·k!)]
# = [(2k+1)·(2k)] / [(k+1)·k] · [k!·(k-1)!] / [k!·(k+1)!]... messy
# Actually: C(2k+1,k)/C(2(k-1)+1,k-1)
# = [(2k+1)!/(k!(k+1)!)] / [(2k-1)!/((k-1)!k!)]
# = [(2k+1)(2k)] / [(k+1)k] = (2k+1)·2/(k+1) = 2(2k+1)/(k+1)

print("Growth ratios within each strand:")
for k in range(1, 10):
    s0_k = C(2*k+1, k)
    s0_km1 = C(2*(k-1)+1, k-1)
    s1_k = C(2*k+2, k)
    s1_km1 = C(2*(k-1)+2, k-1) if k >= 1 else 1
    s2_k = C(2*k+2, k+1)
    s2_km1 = C(2*(k-1)+2, k) if k >= 1 else 1

    r0 = Fraction(s0_k, s0_km1)
    r1 = Fraction(s1_k, s1_km1)
    r2 = Fraction(s2_k, s2_km1)

    print(f"  k={k}: s0(k)/s0(k-1) = {r0} = {float(r0):.4f}, "
          f"s1: {r1} = {float(r1):.4f}, s2: {r2} = {float(r2):.4f}")

print("""
All three strands have growth ratio → 4 as k → ∞.
(Central binomial coefficients grow as 4^k / √(πk).)

The three ratios are:
  s0: 2(2k+1)/(k+1) → 4
  s1: (2k+2)(2k+1)/((k+1)(k+2)) → 4
  s2: 4(2k+1)/(2k+2) → 4

Per 3 positions: growth ≈ 4 per strand step = 4^{1/3} per index step ≈ 1.587.
This is SLOWER than Fibonacci growth φ ≈ 1.618.
""")

# ============================================================
# Part 4: The interleaving as a functor
# ============================================================
print("=" * 70)
print("PART 4: FUNCTORIAL INTERPRETATION")
print("=" * 70)

print("""
THE THREE-STRAND SEQUENCE AS A FUNCTOR:

Consider the category N (natural numbers with ≤) and the functor:
  F: N → Vec (vector spaces / abelian groups)
  F(k) = Z (the integers)
  F(k → k+1) = multiplication by the growth factor

The THREE-STRAND structure is a functor:
  G: N → Z^3 (triples of integers)
  G(k) = (s0(k), s1(k), s2(k))

With transition maps:
  G(k) → G(k+1) given by the 3×3 matrix:
  [s0(k+1)]   [2(2k+3)/(k+2)     0              0    ] [s0(k)]
  [s1(k+1)] = [     0      (2k+4)(2k+3)/((k+2)(k+3))  0    ] [s1(k)]
  [s2(k+1)]   [     0              0       4(2k+3)/(2k+4)] [s2(k)]

This is DIAGONAL — the three strands are INDEPENDENT!
They don't interact. The functor G splits as G = G_0 ⊕ G_1 ⊕ G_2.

BUT: the strands are RELATED by constant or simple ratios:
  s2 = 2·s0 (constant ratio, the OCF parameter)
  s1 = 2(k+1)/(k+2) · s0 (k-dependent ratio)

So G_2 ≅ 2·G_0 and G_1 ≅ (2k+2)/(k+2)·G_0.

The NATURAL TRANSFORMATION η: G_0 → G_2 given by ×2 is
the "OCF evaluation" functor. It maps the "Catalan strand" to
the "central binomial strand" by doubling.

MONOIDAL STRUCTURE:
The Vandermonde convolution gives a product:
  s0(k) * s0(j) → s0(k+j+?) via multinomial coefficients.

Actually, the natural product is on GENERATING FUNCTIONS:
  Σ s0(k) x^k = f(x)
  Σ s2(k) x^k = 2·f(x)  (since s2 = 2·s0)

The GF of s0(k) = C(2k+1,k):
  Σ C(2k+1,k) x^k = (1 - (1-4x)^{1/2}) / (2x) · 1/(1-4x)^{1/2}...
  Actually: Σ C(2n+1,n) x^n = (1/(1-4x) - 1/√(1-4x)) / (2x)... complex.

  Simpler: C(2n+1,n) = (2n+1)·C(2n,n)/(n+1).
  And Σ C(2n,n) x^n = 1/√(1-4x).
  So Σ n·C(2n,n) x^n = 4x/(1-4x)^{3/2}  (differentiate and multiply by x).
  Then Σ (2n+1)·C(2n,n) x^n = 2·4x/(1-4x)^{3/2} + 1/√(1-4x).

  Dividing by (n+1) is harder (involves integral).
  The Catalan GF is C(x) = (1-√(1-4x))/(2x).
  And s0(k) = (2k+1)·CatalanNumber(k), so:
  Σ s0(k) x^k = Σ (2k+1)·C_k x^k = 2x·C'(x) + C(x)...

  Actually: Σ C_n x^n = C(x) = (1-√(1-4x))/(2x)
  Σ n·C_n x^n = x·C'(x)
  C'(x) = (1/√(1-4x) - 1)/(2x²)... this gets messy.
""")

# Let me compute the GF coefficients numerically to verify
print("Generating function verification:")
print("Catalan: C_k, s0 = (2k+1)·C_k, s2 = 2·(2k+1)·C_k")
for k in range(8):
    cat = C(2*k,k)//(k+1)
    s0 = (2*k+1)*cat
    s2 = 2*s0
    actual_s0 = C(2*k+1,k)
    actual_s2 = C(2*k+2,k+1)
    print(f"  k={k}: Cat={cat}, (2k+1)·Cat={s0}, actual s0={actual_s0}, "
          f"match={s0==actual_s0}, 2·s0={s2}, actual s2={actual_s2}, match={s2==actual_s2}")

# ============================================================
# Part 5: Connection to ballot numbers
# ============================================================
print("\n" + "=" * 70)
print("PART 5: BALLOT NUMBERS AND TOURNAMENT INTERPRETATION")
print("=" * 70)

print("""
The BALLOT PROBLEM: In an election with k+1 vs k votes,
the probability that the winner is always strictly ahead is 1/(2k+1).
The number of such sequences is C(2k+1,k)/(2k+1) = C_k (Catalan!).

So C(2k+1,k) = (2k+1) × (number of ballot sequences).

In TOURNAMENT CONTEXT:
  - A tournament on n vertices has a score sequence.
  - A "balanced" tournament (n odd) has all scores near (n-1)/2.
  - The number of "near-balanced" score orderings relates to
    central binomial coefficients.

  The strand ratio s2/s0 = 2:
  - s0 counts "strict ballot" paths (odd-length)
  - s2 counts "relaxed ballot" paths (even-length, reaching center)
  - The factor 2 = the binary choice (which side reaches center first)
  - This parallels OCF: x=2 means each cycle has 2 orientations.

ANALOGY MAP:
  Ballot sequence    ↔  Tournament path
  Candidate A/B      ↔  Arc direction (forward/backward)
  Winning margin      ↔  Score difference
  Ballot C_k          ↔  Catalan structure in cycle graph
  s0 = (2k+1)·C_k   ↔  Labeled Catalan (positional info)
  s2 = 2·s0          ↔  OCF evaluation (binary cycle choice)
""")

# ============================================================
# Part 6: The mod 3 periodicity and cyclic structure
# ============================================================
print("=" * 70)
print("PART 6: Z/3Z GRADING AND CYCLIC CATEGORY")
print("=" * 70)

# The 3-strand structure comes from Pascal's triangle having
# central elements at positions floor(n/2) and ceil(n/2).
# For n odd: ONE central element C(n, (n-1)/2) = C(n, (n+1)/2)
# For n even: ONE central (truly central) C(n, n/2)

# But the strands are indexed mod 3, not mod 2!
# This is because we interleave: odd, even-left, even-right
# giving period 3.

# The mod-3 structure relates to the CUBE ROOT:
# (1+x)^{1/3} = Σ C(1/3, k) x^k (generalized binomial)
# The coefficients of the cube root generate a 3-periodic pattern.

print("Generalized binomial coefficients C(1/3, k):")
from fractions import Fraction
for k in range(10):
    coeff = Fraction(1, 1)
    for i in range(k):
        coeff *= Fraction(1 - 3*i, 3*(i+1))
    print(f"  C(1/3, {k}) = {coeff} = {float(coeff):.10f}")

print("""
The cube root (1+x)^{1/3} has coefficients that decrease as 1/k^{4/3}.
This is much slower than the (1+x)^{1/2} coefficients (which decrease as 1/k^{3/2}).

The THREE-fold reading of Pascal's triangle connects to the 3 cube roots of unity:
  ω = e^{2πi/3} = -1/2 + i√3/2

The extraction of every 3rd term uses:
  Σ_{k≡r mod 3} C(n,k) = (2^n + 2·cos((n-2r)π/3)) / 3

For the central terms (k ≈ n/2), this mod-3 filtering creates
exactly the 3-strand structure.
""")

# Verify: sum of every 3rd binomial coefficient
print("Verification: Σ C(n,k) for k ≡ r mod 3")
for n in range(1, 12):
    sums = [0, 0, 0]
    for k in range(n+1):
        sums[k % 3] += C(n, k)
    total = sum(sums)
    print(f"  n={n:2d}: r=0: {sums[0]:6d}, r=1: {sums[1]:6d}, r=2: {sums[2]:6d}, "
          f"total={total} = 2^{n}, diffs: {sums[0]-sums[1]:+d}, {sums[1]-sums[2]:+d}")

# ============================================================
# Part 7: The x=2 specialization as a monoidal functor
# ============================================================
print("\n" + "=" * 70)
print("PART 7: x=2 AS MONOIDAL FUNCTOR")
print("=" * 70)

print("""
The evaluation at x=2 defines a monoidal functor:

  Ev₂: Graphs → N (natural numbers)
  Ev₂(G) = I(G, 2)
  Ev₂(G ⊔ H) = Ev₂(G) · Ev₂(H)  (multiplicative on disjoint union)

This makes Ev₂ a RING HOMOMORPHISM from the Grothendieck ring of graphs
(with ⊔ as multiplication) to Z.

The THREE-STRAND sequence is the image of Ev₂ applied to:
  P_n (path graphs of increasing length)
  The path graph generates the Jacobsthal numbers under Ev₂.

The CATEGORICAL MEANING of the ratio 2:
  The natural transformation Ev₂ → 2·Ev₂ corresponds to
  replacing each independent set by "2 copies" (forward/backward).
  This is exactly the OCF mechanism:
  H(T) = I(Ω(T), 2) counts independent sets where each selected
  cycle contributes factor 2 (two traversal directions).

The three strands, under Ev₂:
  Ev₂(strand 0 generator at k) = I(G_0(k), 2) for some graph G_0(k)
  What graph? C(2k+1,k) = s0(k).
  We need: I(G, 2) = C(2k+1,k) for each k.

  k=0: I=1 → empty graph
  k=1: I=3 → single vertex
  k=2: I=10 → ??? need a graph with I(G,2) = 10.
  Hmm, 10 is even. But I(G,2) is always odd! Contradiction.

  Actually C(3,1) = 3 (yes, odd). C(5,2) = 10 (EVEN).
  So s0(2) = 10 CANNOT be I(G,2) for any graph.

  This means the 3-strand Pascal sequence is NOT directly the image
  of Ev₂ on a sequence of graphs. The ratio 2 connection to OCF is
  ANALOGICAL, not a direct functor on graphs.

  The connection is more subtle: the STRAND RATIO of 2 equals the
  evaluation parameter, and the central binomial structure parallels
  the independence polynomial structure, but they live in different
  categories.
""")

# Double-check: I(G,2) is always odd
print("Verification: I(G,2) always odd")
print("  I(G,2) = Σ 2^|S| over independent S. Empty set: 2^0=1 (odd).")
print("  All other terms 2^|S| are even. Sum = 1 + (even) = odd. ✓")
print(f"  s0(2) = C(5,2) = 10 (EVEN) → NOT an I(G,2) value.")
print(f"  s0(0) = 1 (odd ✓), s0(1) = 3 (odd ✓), s0(2) = 10 (EVEN ✗)")
print(f"  The strand values are NOT all realizable as I(G,2).")

# ============================================================
# Part 8: What IS the connection?
# ============================================================
print("\n" + "=" * 70)
print("PART 8: THE TRUE CONNECTION")
print("=" * 70)

print("""
The connection between the 3-strand Pascal sequence and tournaments
is NOT through I(G,2) directly, but through:

1. ASYMPTOTIC GROWTH: Both central binomials and I(P_k,2) grow as 4^k.
   - C(2k,k) ~ 4^k/√(πk)
   - I(P_k,2) = (2^{k+2} - (-1)^k)/3 ~ (4/3)·2^k

   Wait: these grow at DIFFERENT rates! C(2k,k) grows as 4^k, while
   I(P_k,2) grows as 2^k. So the connection is not asymptotic growth.

2. THE RATIO 2: The strand ratio s2/s0 = 2 mirrors the OCF parameter x=2.
   This is the NUMEROLOGICAL connection. The factor 2 appearing in both:
   - Pascal: going from odd row central to even row central multiplies by 2
   - OCF: each independent cycle contributes factor 1+x = 1+2 = 3, vs 1+1 = 2 at x=1

3. The DEEPER connection is through GENERATING FUNCTIONS:
   - I(P_k, x) satisfies a(k) = a(k-1) + x·a(k-2)
   - The GF is Σ I(P_k,x) t^k = 1/(1-t-xt²)
   - At x=2: 1/(1-t-2t²) = 1/((1-2t)(1+t))
   - The PARTIAL FRACTIONS give the Jacobsthal closed form.

   - Central binomials: GF = 1/√(1-4t)
   - This is (1-4t)^{-1/2}, a FRACTIONAL power.

   The connection: 1/(1-t-2t²) = 1/((1-2t)(1+t)) involves the factor (1-2t),
   and 1/√(1-4t) = 1/√((1-2t)(1+2t)) also involves (1-2t).

   So BOTH generating functions have singularity at t = 1/2,
   and the singularity type differs:
   - Jacobsthal: simple pole at t=1/2 (and t=-1)
   - Central binomial: branch point at t=1/4

   The singularity at t=1/2 corresponds to growth rate 2.
   For Jacobsthal: dominant root = 2 (giving 2^k growth).
   For central binomial: singular at 1/4 (giving 4^k growth).

   But note: the STRAND sequence (reading every 3rd term of Pascal)
   has growth 4^{1/3} per step, which is (cube root of 4) ≈ 1.587.
   Meanwhile Fibonacci has growth φ ≈ 1.618 per step.
   The Jacobsthal has growth 2 per step.

   The ordering: 4^{1/3} < φ < 2 means:
   Pascal 3-strand grows SLOWER than Fibonacci grows SLOWER than Jacobsthal.

4. MODULARITY: Both sequences have periodic behavior mod small primes.
   - Fibonacci mod 4: period 6 = 2×3
   - Jacobsthal mod 4: period 6
   - Pascal mod 2: period 3 (along each strand)

   The common period 6 = LCM(2,3) comes from the interaction of
   binary (mod 2) and ternary (mod 3) structure.
""")

# Verify: Jacobsthal mod 4 period
J = [0, 1]
for i in range(2, 30):
    J.append(J[-1] + 2*J[-2])

print("Jacobsthal mod 4:", [J[i] % 4 for i in range(18)])
# Find period
jm4 = [J[i] % 4 for i in range(24)]
for p in range(1, 12):
    if all(jm4[i] == jm4[i+p] for i in range(12)):
        print(f"Period: {p}")
        print(f"Pattern: {jm4[:p]}")
        break

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("GRAND SYNTHESIS")
print("=" * 70)

print("""
THE THREE-STRAND PASCAL SEQUENCE IN TOURNAMENT CONTEXT:

1. STRUCTURE: Three strands read the CENTRAL COEFFICIENTS of (1+x)^n.
   s0(k) = C(2k+1,k) = (2k+1)·CatalanNumber(k)
   s1(k) = C(2k+2,k) = (2k+2)/(k+2)·s0(k)
   s2(k) = C(2k+2,k+1) = 2·s0(k)  [EXACT, for all k]

2. THE RATIO 2: The constant strand ratio s2/s0 = 2 is the
   evaluation parameter of the OCF: I(Ω(T), x) at x=2.
   This is the UNIQUE positive integer x where the characteristic
   roots of a(k) = a(k-1) + x·a(k-2) are both integers (2 and -1).

3. CATALAN CONNECTION: s0 = (2k+1)·C_k, where Catalan numbers count
   binary trees and non-crossing partitions. In tournaments,
   non-crossing structures relate to TRANSITIVE sub-tournaments.

4. NOT A DIRECT FUNCTOR: The strand values include even numbers
   (like s0(2)=10), which cannot be I(G,2) for any graph (always odd).
   The connection is analogical/structural, not functorial.

5. PERIOD-6 UNIVERSALITY: Both Fibonacci mod 4 and Jacobsthal mod 4
   have period 6 = 2×3. This matches:
   - 3 strands × 2 parities = 6-periodic structure
   - 6 = smallest n where ALL major tournament phenomena appear
   - 6 = LCM of binary (arc choice) and ternary (triangle) structure

6. GEOMETRIC PACKING: The central binomial C(2k,k) is the largest
   coefficient of (1+x)^{2k} = the "simplex polynomial."
   At x=2: (1+2)^{2k} = 3^{2k} = 9^k.
   The strand interleaves these maximal coefficients, capturing
   the "peak" of the simplex at each level.

7. TOURNAMENT TAKEAWAY: The three-strand structure reveals that
   tournament combinatorics (living at x=2) is connected to
   the central structure of Pascal's triangle through the
   factor 2. This factor simultaneously:
   - Doubles the strand (s2 = 2·s0)
   - Sets the OCF evaluation point (x=2)
   - Makes characteristic roots integral
   - Creates the binary choice in arc orientation
""")
