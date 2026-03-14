#!/usr/bin/env python3
"""
packing_both_directions.py — opus-2026-03-14-S71e

PACKING IN BOTH DIRECTIONS: Simplex in cuboid AND cuboid in simplex

Key identities:
  (x+2)^n = sum C(n,k) (x+1)^k         [cuboid = sum of simplex pieces]
  (x+1)^n = sum C(n,k) (-1)^{n-k} (x+2)^k   [simplex = alternating sum of cuboid pieces!]

The second identity follows from the binomial inversion formula.

Also: (x+r)^n for general r gives the "r-hedroid" family:
  r=1: simplex (x+1)^n
  r=2: cuboid (x+2)^n
  r=3: "3-hedroid" (x+3)^n
  ...

The ENTIRE family is related by:
  (x+r+1)^n = sum C(n,k) (x+r)^k

This is the PASCAL TRIANGLE operating on the hedroid family!

CONNECTION TO TOURNAMENTS:
  I(Omega, x) is evaluated at x = 2 (simplex) and x = 3 (cuboid).
  The independence polynomial I itself plays the role of a
  "generalized binomial" connecting the simplex and cuboid evaluations.
"""

import sys
from math import comb, factorial
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

print("=" * 70)
print("PART 1: BINOMIAL INVERSION — SIMPLEX FROM CUBOID PIECES")
print("=" * 70)

print("""
  FORWARD: (x+2)^n = sum_{k=0}^n C(n,k) (x+1)^k
  INVERSE: (x+1)^n = sum_{k=0}^n C(n,k) (-1)^{n-k} (x+2)^k

  Proof: Let y = x+1. Then x+2 = y+1.
  (y+1)^n = sum C(n,k) y^k  ← forward
  y^n = sum C(n,k) (-1)^{n-k} (y+1)^k  ← binomial inversion

  This means: SIMPLEX = ALTERNATING SUM of cuboid pieces!
""")

print("  VERIFICATION:")
for n in range(1, 8):
    # Check: (x+1)^n = sum C(n,k)(-1)^{n-k}(x+2)^k at x=1
    # LHS = 2^n
    lhs = 2**n
    # RHS = sum C(n,k)(-1)^{n-k} 3^k
    rhs = sum(comb(n, k) * (-1)**(n-k) * 3**k for k in range(n+1))
    print(f"  n={n}: (x+1)^{n}|_{{x=1}} = {lhs}, sum C({n},k)(-1)^{{{n}-k}} 3^k = {rhs}, match = {lhs == rhs}")

print(f"\n  KEY: The alternating signs mean INCLUSION-EXCLUSION!")
print(f"  The simplex is obtained from cuboid pieces via Möbius inversion.")
print(f"  In tournament terms: H = I(2) can be recovered from I(3) and lower I values.")

print("\n" + "=" * 70)
print("PART 2: THE GENERAL r-HEDROID FAMILY")
print("=" * 70)

print("""
  Define F_r(x) = (x+r)^n for fixed n.

  Then: F_{r+1}(x) = sum C(n,k) F_r^{(k)}(x) ... no, simpler:
  F_{r+1}(x) = (x+r+1)^n = ((x+r)+1)^n = sum C(n,k) (x+r)^k
             = sum C(n,k) F_r(x)^{k/n} ... no, that's wrong.

  Actually: (x+r+1)^n = sum C(n,k) (x+r)^k * 1^{n-k}
  This gives F_{r+1} in terms of POWERS of the ARGUMENT (x+r), not F_r itself.

  The correct relation is simpler:
  (x+r+1)^n = sum C(n,k) (x+r)^k

  At x=1: (r+2)^n = sum C(n,k) (r+1)^k

  This IS the binomial theorem: (a+1)^n = sum C(n,k) a^k where a = r+1.

  For the HEDROID FAMILY at x=1:
    r=0: F_0 = 1^n = 1
    r=1: F_1 = 2^n (simplex)
    r=2: F_2 = 3^n (cuboid)
    r=3: F_3 = 4^n (tesseract)
    r=4: F_4 = 5^n (penteract)

  And: F_{r+1} = sum C(n,k) (r+1)^k = (r+2)^n ✓ (trivially true)

  The INTERESTING structure is the DIFFERENCE OPERATOR:
  Delta F_r = F_{r+1} - F_r = (r+2)^n - (r+1)^n

  This gives: Delta F_0 = 2^n - 1 (Mersenne numbers)
              Delta F_1 = 3^n - 2^n (simplex-cuboid gap)
              Delta F_2 = 4^n - 3^n
              Delta F_3 = 5^n - 4^n
""")

print("  FINITE DIFFERENCE TABLE (n fixed, varying r):")
for n in [3, 5, 7, 9]:
    print(f"\n  n={n}:")
    # Build the difference table
    f = [(r+1)**n for r in range(8)]
    print(f"    F:  {f}")

    diffs = [f]
    for level in range(1, 6):
        d = [diffs[-1][i+1] - diffs[-1][i] for i in range(len(diffs[-1])-1)]
        diffs.append(d)
        if len(d) > 0:
            print(f"    Δ^{level}: {d}")
        if all(x == d[0] for x in d):
            print(f"    *** CONSTANT at level {level}: Δ^{level} = {d[0]} = {n}! = {factorial(n)}? {'YES' if d[0] == factorial(n) else 'NO'} ***")
            break

print("""
  OBSERVATION: The n-th finite difference Δ^n F_r is ALWAYS n!
  This is because F_r = (r+1)^n is a polynomial of degree n in r,
  and its leading coefficient is 1, so Δ^n = n! * (leading coeff) = n!.

  THIS IS THE FINITE DIFFERENCE CALCULUS CONNECTION:
  n! = Δ^n[(r+1)^n]|_{r=0} = sum_{k=0}^n (-1)^{n-k} C(n,k) (k+1)^n

  This is the STIRLING NUMBER formula!
  n! = sum_{k=0}^n (-1)^{n-k} C(n,k) (k+1)^n

  TOURNAMENT INTERPRETATION:
  The number of Hamiltonian paths in a transitive tournament = 1.
  n! = number of permutations of n.
  The finite difference Δ^n at the hedroid family gives n!.
  So: the PERMUTATION GROUP (n!) is the "top-level finite difference"
  of the hedroid counting sequence.
""")

print("  VERIFICATION: n! = sum (-1)^{n-k} C(n,k) (k+1)^n")
for n in range(1, 10):
    val = sum((-1)**(n-k) * comb(n,k) * (k+1)**n for k in range(n+1))
    print(f"  n={n}: sum = {val}, n! = {factorial(n)}, match = {val == factorial(n)}")

print("\n" + "=" * 70)
print("PART 3: THE INVERSE — EXPRESSING I(x) IN THE HEDROID BASIS")
print("=" * 70)

print("""
  The independence polynomial I(Omega, x) = sum alpha_k x^k.

  In the HEDROID BASIS, we can write I(Omega, x) as:
  I(Omega, x) = sum_{j=0}^{d} c_j (x+1)^j   [simplex basis]
  or:
  I(Omega, x) = sum_{j=0}^{d} d_j (x+2)^j   [cuboid basis]

  The simplex basis coefficients c_j:
  I(Omega, x) = sum alpha_k x^k = sum alpha_k ((x+1)-1)^k
              = sum alpha_k sum C(k,j) (-1)^{k-j} (x+1)^j
              = sum_j (sum_k alpha_k C(k,j) (-1)^{k-j}) (x+1)^j

  So: c_j = sum_{k=j}^{d} alpha_k C(k,j) (-1)^{k-j}

  This is the ALTERNATING BINOMIAL TRANSFORM of the alpha sequence!
""")

# Compute for some example tournaments (n=5)
print("  EXAMPLE: n=5 tournaments in simplex and cuboid bases")
print("  (Using alpha_1, alpha_2 from known cases)")
print()

cases = [
    ("Transitive", 0, 0),
    ("Near-trans", 1, 0),
    ("Moderate", 3, 0),
    ("Moderate+", 4, 1),
    ("Regular C5", 5, 2),  # Circulant
]

for name, a1, a2 in cases:
    H = 1 + 2*a1 + 4*a2
    I3 = 1 + 3*a1 + 9*a2

    # Standard basis: I(x) = 1 + a1*x + a2*x^2
    # Simplex basis: I(x) = c0 + c1*(x+1) + c2*(x+1)^2
    # c2 = alpha_2
    # c1 = alpha_1 - 2*alpha_2
    # c0 = 1 - alpha_1 + alpha_2
    c0 = 1 - a1 + a2
    c1 = a1 - 2*a2
    c2 = a2

    # Cuboid basis: I(x) = d0 + d1*(x+2) + d2*(x+2)^2
    # x = (x+2) - 2
    # I(x) = 1 + a1*(x+2-2) + a2*(x+2-2)^2
    #       = 1 + a1*(x+2) - 2*a1 + a2*(x+2)^2 - 4*a2*(x+2) + 4*a2
    #       = (1 - 2*a1 + 4*a2) + (a1 - 4*a2)*(x+2) + a2*(x+2)^2
    d0 = 1 - 2*a1 + 4*a2
    d1 = a1 - 4*a2
    d2 = a2

    # Verify
    # At x=2: c0 + c1*3 + c2*9 should = H
    h_check = c0 + c1*3 + c2*9
    # At x=2: d0 + d1*4 + d2*16 should = H
    h_check2 = d0 + d1*4 + d2*16

    print(f"  {name:12s}: α=({a1},{a2}), H={H}, I(3)={I3}")
    print(f"    Standard:  1 + {a1}x + {a2}x²")
    print(f"    Simplex:   {c0} + {c1}(x+1) + {c2}(x+1)²")
    print(f"    Cuboid:    {d0} + {d1}(x+2) + {d2}(x+2)²")
    print(f"    H check (simplex): {h_check} {'✓' if h_check == H else '✗'}")
    print(f"    H check (cuboid):  {h_check2} {'✓' if h_check2 == H else '✗'}")
    print()

print("\n" + "=" * 70)
print("PART 4: THE FALLING FACTORIAL / STIRLING CONNECTION")
print("=" * 70)

print("""
  The change of basis from x^k to (x+1)^j involves Stirling numbers!

  x^k = sum_j S(k,j) * x^{(j)}  where x^{(j)} = x(x-1)...(x-j+1) = falling factorial

  And: (x+1)^n = sum_k C(n,k) x^k = sum_k C(n,k) sum_j S(k,j) x^{(j)}

  The SIMPLEX BASIS expansion of I(Omega, x) involves the STIRLING TRANSFORM
  of the alpha_k sequence.

  For tournament theory:
    I(x) in standard basis: alpha_0=1, alpha_1, alpha_2, ...
    I(x) in falling factorial basis: involves Stirling numbers of the first kind
    I(x) in simplex basis (x+1)^j: involves the shifted Stirling transform

  THE DEEP CONNECTION:
  The Stirling numbers S(n,k) count set partitions.
  The alpha_k count independent cycle sets.
  The SIMPLEX BASIS coefficients c_j mix these via Stirling numbers.

  This means: THE HEDROID PACKING IS A STIRLING TRANSFORM OF THE CYCLE STRUCTURE.
""")

# Compute Stirling numbers of the second kind (for small values)
def stirling2(n, k):
    """Stirling numbers of the second kind."""
    if n == 0 and k == 0:
        return 1
    if n == 0 or k == 0:
        return 0
    return k * stirling2(n-1, k) + stirling2(n-1, k-1)

print("  Stirling numbers S(n,k) of the second kind:")
print("  " + "".join(f"  k={k}" for k in range(7)))
for n in range(7):
    print(f"  n={n}: ", end="")
    for k in range(7):
        s = stirling2(n, k)
        print(f"  {s:3d}", end="")
    print()

print("""
  KEY OBSERVATION:
  x^n = sum_k S(n,k) x(x-1)...(x-k+1) = sum_k S(n,k) x^{(k)}

  At x=2 (simplex point):
  2^n = sum_k S(n,k) * 2^{(k)} = sum_k S(n,k) * 2*1*0*(-1)... = S(n,1)*2 + S(n,2)*2

  Wait: 2^{(0)} = 1, 2^{(1)} = 2, 2^{(2)} = 2*1 = 2, 2^{(3)} = 2*1*0 = 0.

  So: 2^n = S(n,0)*1 + S(n,1)*2 + S(n,2)*2 + 0 + 0 + ...
         = 1 + 2*S(n,1) + 2*S(n,2)... no, S(n,0)=0 for n>0.
         = 2*1 + 2*S(n,2)    [since S(n,1)=1 always]
         = 2 + 2*S(n,2)

  Check: n=3: 2^3 = 8 = 2 + 2*S(3,2) = 2 + 2*3 = 8. ✓
         n=4: 2^4 = 16 = 2 + 2*S(4,2) = 2 + 2*7 = 16. ✓

  So: S(n,2) = (2^n - 2)/2 = 2^{n-1} - 1.

  This is: S(n,2) = 2^{n-1} - 1 = MERSENNE NUMBER!
  And: 2^n = 2 + 2*(2^{n-1}-1) = 2 + 2^n - 2 = 2^n. Trivially true but illuminating.

  For I(Omega, 2) = H:
  H = sum alpha_k * 2^k = sum alpha_k * (2 + 2*S(k,2) + higher terms with 0)
    = sum alpha_k * 2 + sum alpha_k * 2 * S(k,2)
  But S(k,2) = 2^{k-1}-1, so:
  H = 2*sum(alpha_k) + 2*sum(alpha_k*(2^{k-1}-1))
    = 2*sum(alpha_k) + sum(alpha_k*(2^k - 2))
    = 2*sum(alpha_k) + H - 2*sum(alpha_k) = H. Circular.
""")

print("\n" + "=" * 70)
print("PART 5: THE INCLUSION-EXCLUSION INTERPRETATION")
print("=" * 70)

print("""
  THE MOST IMPORTANT IDENTITY:

  (x+1)^n = sum C(n,k) (-1)^{n-k} (x+2)^k

  At x=1: 2^n = sum C(n,k) (-1)^{n-k} 3^k

  This is INCLUSION-EXCLUSION!

  Think of 3^n as counting functions from [n] to {a,b,c} (3 colors).
  Then 2^n counts functions using only {a,b} (2 colors).
  By inclusion-exclusion:
    |functions using ≤2 colors| = sum (-1)^k C(3,k) (3-k)^n
    But that's for choosing which colors to EXCLUDE.

  Actually, the simpler interpretation:
  2^n = (3-1)^n = sum C(n,k) 3^k (-1)^{n-k}

  This is just the binomial theorem with a = 3 and b = -1.

  But in the TOURNAMENT context:
  H = I(2) = I(3-1) = sum alpha_k (3-1)^k
            = sum alpha_k sum C(k,j) 3^j (-1)^{k-j}
            = sum_j (sum_k alpha_k C(k,j) (-1)^{k-j}) 3^j

  This expresses H as a LINEAR COMBINATION of I(3)^{powers}!
  (Not exactly, because I(3) ≠ 3, but the STRUCTURE is analogous.)

  THE PACKING IN REVERSE:
  To go from cuboid (I(3)) to simplex (I(2)=H), we use
  alternating signs — INCLUSION-EXCLUSION of the cycle structure.

  The corners that we ADD when going simplex→cuboid,
  we SUBTRACT when going cuboid→simplex.

  GEOMETRICALLY: the simplex is the cuboid MINUS corners.
  But the corners overlap, so we need inclusion-exclusion.
  The alternating signs in C(n,k)(-1)^{n-k} ARE the
  inclusion-exclusion correction terms.
""")

print("  NUMERICAL VERIFICATION: 2^n = sum C(n,k)(-1)^{n-k} 3^k")
for n in range(1, 10):
    lhs = 2**n
    terms = [(comb(n,k) * (-1)**(n-k) * 3**k, k) for k in range(n+1)]
    rhs = sum(t[0] for t in terms)
    print(f"  n={n}: {lhs} = ", end="")
    print(" + ".join(f"({t[0]:+d})" for t in terms), end="")
    print(f" = {rhs} {'✓' if lhs == rhs else '✗'}")

print("\n" + "=" * 70)
print("PART 6: THE EULER TRANSFORM AND TOURNAMENT GENERATING FUNCTIONS")
print("=" * 70)

print("""
  The Euler transform of a sequence {a_n} is:
  E(x) = prod_{k=1}^∞ (1-x^k)^{-a_k}

  For tournaments, the "partition function" is:
  Z(x) = sum_T x^{|T|} (sum over all tournaments, weighted by size)

  The independence polynomial I(Omega, x) can be seen as a
  "local partition function" for the odd-cycle conflict graph.

  The GENERATING FUNCTION for the hedroid family:
  F(r, n) = (r+1)^n

  Its EGF in n: sum_n F(r,n) x^n/n! = sum_n (r+1)^n x^n/n! = e^{(r+1)x}

  So: F(r, ·) has EGF e^{(r+1)x}.
  And: F(0,·) = e^x (simplex-at-0)
       F(1,·) = e^{2x} (simplex)
       F(2,·) = e^{3x} (cuboid)

  The RATIO: F(2,·)/F(1,·) = e^{3x}/e^{2x} = e^x
  This IS the EGF of the constant sequence 1,1,1,...!

  In other words: the cuboid-to-simplex ratio at each n is
  (coefficient of x^n/n! in e^x) / (coefficient of x^n/n! in e^{2x})
  = 1/2^n (after normalization)... no, that's not quite right.

  The point is: e^{3x} = e^{2x} * e^x
  So: 3^n = sum_{k=0}^n C(n,k) 2^k * 1^{n-k}
  = sum C(n,k) 2^k  ← binomial theorem ✓

  The EGF factorization e^{3x} = e^{2x} * e^x corresponds to:
  cuboid = simplex ⊗ complement
  where ⊗ is the BINOMIAL CONVOLUTION (Hadamard product of EGFs).

  THE TOURNAMENT ANALOG:
  If I(Omega, x) has EGF-like structure, then
  I(3) = I(2) * "correction factor"
  where the correction factor involves the cycle structure.
""")

# The EGF ratio
print("  EGF RATIO e^{3x}/e^{2x} = e^x:")
print("  The coefficient of x^n/n! in e^x is 1.")
print("  So the 'correction factor' at each n is: 3^n/2^n = (3/2)^n.")
print("  This is the simplex-cuboid ratio at dimension n.")
print()

for n in range(1, 12):
    ratio = (3/2)**n
    # Express as sum: 3^n = sum C(n,k) 2^k
    # The sum has n+1 terms, dominated by the k=n term (2^n)
    # The correction: sum_{k=0}^{n-1} C(n,k) 2^k = 3^n - 2^n
    correction = 3**n - 2**n
    print(f"  n={n:2d}: (3/2)^n = {ratio:10.4f}, correction = {correction:10d}, correction/2^n = {correction/2**n:.4f}")

print("\n  The correction/2^n = (3/2)^n - 1 → ∞ as n → ∞.")
print("  The simplex becomes NEGLIGIBLE inside the cuboid.")
print("  But the RATIO (3/2)^n grows EXPONENTIALLY.")

print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print("""
  THE SIMPLEX-CUBOID PACKING FRAMEWORK SUMMARIZED:

  1. POLYNOMIALS: (x+1)^n ↔ simplex, (x+2)^n ↔ cuboid
  2. FORWARD: cuboid = sum C(n,k) simplex^k (binomial expansion)
  3. INVERSE: simplex = alternating sum of cuboid^k (Möbius inversion)
  4. GEOMETRY: Regular simplex inscribable in cube at Hadamard dimensions (n=3,7,11,...)
     At n=3: 4 corner pieces, each = half the simplex volume
  5. TOURNAMENTS: I(2)=H (simplex), I(3) (cuboid), gap = cycle count
  6. k-NACCI: limit 2 (simplex), doubled limit 3 (cuboid), ratio → 3/2
  7. CORNER PIECES: 3^k - 2^k = 1, 5, 19, 65, 211, ... (level-k excess)
  8. 9=3²: nested packing (simplex of simplices)
  9. 10=5*2: cross-product packing (mixed prime structure)
  10. n! = Δ^n[(r+1)^n]: permutations = top finite difference of hedroid family
  11. EGF: e^{3x} = e^{2x} · e^x ↔ cuboid = simplex ⊗ correction
  12. STIRLING: S(n,2) = 2^{n-1}-1 = Mersenne numbers = simplex Stirling value
""")

print("Done.")
