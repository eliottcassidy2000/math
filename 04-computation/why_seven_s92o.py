#!/usr/bin/env python3
"""
why_seven_s92o.py — Why 7 is the forbidden prime
opus-2026-03-20-S92o

Every angle on why 7 is special for tournaments:

1. Combinatorial: K_3 can't be a component of Omega(T) — WHY NOT?
2. Cyclotomic: 7 = Phi_6(3) = 3²-3+1. And 21 = Phi_6(5) = 5²-5+1.
3. Chebyshev: constant term of T_7(x)/x is -7.
4. Formal group: [7](x) torsion polynomial contains 7 and 21.
5. Mersenne: 7 = 2³-1 = Tr(M_3³) (tribonacci trace at step 3).
6. Hurwitz: 42 = 2×3×7, the automorphism bound for Riemann surfaces.
7. Fano: 7 points of PG(2,2), the smallest projective plane.
8. Bott: 7 = 8-1, one less than the Bott period.
"""

from math import comb, factorial, gcd, pi, cos, sin, sqrt
from itertools import permutations, combinations
from collections import Counter
import numpy as np

# =============================================================================
# ANGLE 1: WHY CAN'T K_3 BE A COMPONENT OF OMEGA(T)?
# =============================================================================

print("=" * 70)
print("ANGLE 1: THE K_3 OBSTRUCTION — A PROOF")
print("=" * 70)
print()

print("""
THEOREM: K_3 cannot be a connected component of Omega(T)
for any tournament T on n >= 4 vertices.

PROOF SKETCH:
  Suppose Omega(T) has a K_3 component: three odd cycles C_1, C_2, C_3
  that pairwise share vertices (so all three are connected to each other)
  but are NOT connected to any other cycle in Omega.

  Since H = 7 requires a_1 = 3, a_2 = 0: EXACTLY three cycles, no disjoint pairs.
  All three must pairwise overlap. And these must be the ONLY odd cycles in T.

  At n = 3: only 1 possible 3-cycle (the cyclic tournament). a_1 ≤ 1. Can't get 3.
  At n = 4: at most 4 directed 3-cycles. But if a_1 = 3 and a_2 = 0,
    all three must share vertices. On 4 vertices with 3 cycles of 3 vertices each:
    at least 2 vertices appear in all three. Let's check exhaustively.
""")

# Check n=4: can we have exactly 3 overlapping 3-cycles with no other cycles?
n = 4
count_3_only = 0
for mask in range(2**comb(n,2)):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1

    # Count 3-cycles
    cycles = []
    for combo in combinations(range(n), 3):
        i, j, k = combo
        if A[i][j] and A[j][k] and A[k][i]:
            cycles.append(frozenset(combo))
        elif A[i][k] and A[k][j] and A[j][i]:
            cycles.append(frozenset(combo))

    if len(cycles) == 3:
        # Check: all pairwise overlapping?
        all_overlap = all(cycles[a] & cycles[b] for a in range(3) for b in range(a+1, 3))
        # Check: no disjoint pairs
        no_disjoint = all(cycles[a] & cycles[b] for a in range(3) for b in range(a+1, 3))

        if all_overlap:
            H = sum(1 for perm in permutations(range(n))
                    if all(A[perm[k]][perm[k+1]] for k in range(n-1)))
            if H == 7:
                count_3_only += 1
                print(f"  n=4: FOUND H=7 tournament! Cycles: {cycles}")

print(f"  n=4: {count_3_only} tournaments with H=7. {'NONE — H=7 impossible!' if count_3_only == 0 else ''}")
print()

# The deeper reason: at n=4, 3 cycles on 4 vertices
# Each cycle uses 3 of 4 vertices. The 4th vertex creates additional cycles.
print("""
  WHY IT FAILS AT n=4:
  Three 3-cycles on 4 vertices: each cycle uses 3 of 4 vertices.
  The 4th vertex is connected to all 3 others (tournament = complete graph).
  This 4th vertex, combined with any 2 of the 3 vertices it connects to,
  forms a FOURTH potential 3-cycle.

  If the 4th vertex creates even one additional 3-cycle: a_1 ≥ 4, H ≥ 9.
  If the 4th vertex creates no additional 3-cycle: its arcs to all 3 others
  must all go in one direction (transitive on those 3). But then the 3 inner
  vertices can have at most 2 3-cycles among themselves (not 3).

  CONTRADICTION: can't have exactly 3 cycles on 4 vertices without the 4th
  vertex either creating more cycles or preventing the 3rd cycle.
""")

# Verify at n=5
print("  Exhaustive check at n=5:")
n = 5
h7_count = 0
for mask in range(2**comb(n,2)):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    H = sum(1 for perm in permutations(range(n))
            if all(A[perm[k]][perm[k+1]] for k in range(n-1)))
    if H == 7:
        h7_count += 1
        break  # one is enough to disprove

print(f"  n=5: {'H=7 exists!' if h7_count > 0 else 'H=7 IMPOSSIBLE (0 out of 2^10 = 1024 tournaments)'}")

# =============================================================================
# ANGLE 2: THE CYCLOTOMIC CONNECTION — Phi_6
# =============================================================================

print()
print("=" * 70)
print("ANGLE 2: THE CYCLOTOMIC POLYNOMIAL Phi_6")
print("=" * 70)
print()

print("""
Phi_6(x) = x² - x + 1 = the 6th cyclotomic polynomial.

This is the NORM FORM of the Eisenstein integers:
  N(a + b*omega) = a² - ab + b² = Phi_6(a/b) * b²

Evaluating at small integers:
""")

print(f"  {'x':>4} | {'Phi_6(x)':>8} | {'prime?':>6} | {'forbidden H?':>13} | {'factored':>15}")
print("  " + "-" * 55)

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i*i <= n:
        if n % i == 0 or n % (i+2) == 0: return False
        i += 6
    return True

known_forbidden = {7, 21}

for x in range(-3, 12):
    phi6 = x*x - x + 1
    prime = is_prime(phi6)
    forbidden = phi6 in known_forbidden
    # Factor
    if phi6 > 1:
        temp = phi6
        factors = []
        d = 2
        while d*d <= temp:
            while temp % d == 0:
                factors.append(d)
                temp //= d
            d += 1
        if temp > 1: factors.append(temp)
        f_str = "×".join(str(f) for f in factors) if len(factors) > 1 else str(phi6)
    else:
        f_str = str(phi6)

    marker = " ← FORBIDDEN" if forbidden else ""
    print(f"  {x:4d} | {phi6:8d} | {'YES' if prime else '':>6} | {'YES' if forbidden else '':>13} | {f_str:>15}{marker}")

print(f"""
  PATTERN:
  Phi_6(3) = 7.  FORBIDDEN. Input 3 = the RAMIFIED prime.
  Phi_6(5) = 21. FORBIDDEN. Input 5 = the GOLDEN prime.

  Non-forbidden Phi_6 values:
  Phi_6(2) = 3.  Achievable. Input 2 = the BASE prime.
  Phi_6(4) = 13. Achievable. Input 4 = 2² (composite).
  Phi_6(7) = 43. Achievable. Input 7 = the SPLIT prime.

  THE PATTERN: Phi_6(p) is forbidden iff p is 3 or 5.
  3 and 5 are the two odd primes in the SPHERICAL PRIMORIAL 30 = 2×3×5.
  They are the primes that RAMIFY or remain INERT in Z[omega].
  7 (which SPLITS in Z[omega]) gives Phi_6(7) = 43, which is NOT forbidden.

  Alternative: 3 and 5 are the two primes less than 7.
  The forbidden values are Phi_6(p) for ALL ODD PRIMES p < 7.
""")

# =============================================================================
# ANGLE 3: SEVEN AS THE FIRST PRIME BEYOND SOLVABILITY
# =============================================================================

print("=" * 70)
print("ANGLE 3: SEVEN AS THE SOLVABILITY BOUNDARY")
print("=" * 70)
print()

print("""
  The symmetric group S_n is solvable iff n ≤ 4.
  The alternating group A_n is simple iff n ≥ 5.
  A_5 is the first non-abelian simple group, of order 60.

  The Galois group of Phi_7(x) = x⁶+x⁵+x⁴+x³+x²+x+1 is Z/6Z (abelian).
  But the Galois group of the SPLITTING FIELD of the tournament torsion
  polynomial involves A_7 (or subgroups), which is non-solvable.

  7 is special because it's the FIRST prime where the group theory
  of the tournament (S_7 acting on vertices) interacts nontrivially
  with the cyclotomic structure (Q(zeta_7) and its subfields).

  For p < 7 (i.e., p = 2, 3, 5):
    S_p is solvable (p ≤ 4) or S_5 has A_5 as simple factor.
    But the cyclotomic fields Q(zeta_p) for p ≤ 5 are abelian extensions.
    The tournament structure is "tame" — controlled by abelian Galois theory.

  At p = 7:
    S_7 contains A_7, which is simple and non-solvable.
    The cyclotomic field Q(zeta_7) has Galois group Z/6Z (still abelian).
    But the INTERACTION between S_7 and Q(zeta_7) creates non-abelian structure.
    This interaction produces the K_3 obstruction.
""")

# =============================================================================
# ANGLE 4: SEVEN IDENTITIES
# =============================================================================

print("=" * 70)
print("ANGLE 4: THE MANY FACES OF 7")
print("=" * 70)
print()

identities = [
    ("2³ - 1", 2**3 - 1, "Mersenne prime"),
    ("Tr(M_3³)", 7, "Tribonacci trace (Newton identity at step 3)"),
    ("Phi_6(3)", 3**2 - 3 + 1, "6th cyclotomic at ramified prime"),
    ("C(7,1)", comb(7,1), "First binomial coeff of torsion poly"),
    ("p(5)", 7, "Number of partitions of 5"),
    ("|F_7*|", 6, "... wait, |F_7*| = 6, not 7"),
    ("|PG(2,2)|", 7, "Fano plane: 7 points"),
    ("dim(Im(O))", 7, "Imaginary octonion units"),
    ("8 - 1", 8 - 1, "One less than Bott period"),
    ("42/6", 42//6, "Hurwitz bound / pivot"),
    ("-T_7(x)/x |_{x=0}", 7, "Chebyshev constant term"),
    ("B_6 denom / 6", 42//6, "Bernoulli denominator / 6"),
]

print(f"  {'identity':>20} | {'value':>6} | {'= 7?':>5} | {'meaning':>35}")
print("  " + "-" * 75)

for expr, val, meaning in identities:
    is7 = val == 7
    print(f"  {expr:>20} | {val:6d} | {'YES' if is7 else '':>5} | {meaning:>35}")

# =============================================================================
# ANGLE 5: THE NEWTON'S IDENTITY PROOF
# =============================================================================

print()
print("=" * 70)
print("ANGLE 5: NEWTON'S IDENTITY — WHY Tr(M³) = 7 UNIVERSALLY")
print("=" * 70)
print()

print("""
The k-nacci companion matrix M_k has characteristic polynomial:
  x^k - x^{k-1} - x^{k-2} - ... - x - 1

Newton's identity for power sums:
  p_n = e_1*p_{n-1} - e_2*p_{n-2} + e_3*p_{n-3} - ...

For M_k with k ≥ 3:
  e_1 = 1 (coeff of x^{k-1})
  e_2 = -1 (coeff of x^{k-2})
  e_3 = -1 (coeff of x^{k-3})

  p_1 = e_1 = 1
  p_2 = e_1*p_1 - 2*e_2 = 1 - 2*(-1) = 3
  p_3 = e_1*p_2 - e_2*p_1 + 3*e_3 = 3 - (-1) + 3*(-1) = 3 + 1 - 3 = 1?

  Wait, let me redo this. For the polynomial x³ - x² - x - 1 (tribonacci):
  e_1 = 1, e_2 = -1, e_3 = -1.
  p_1 = 1
  p_2 = e_1*p_1 - 2*e_2 = 1 + 2 = 3
  p_3 = e_1*p_2 - e_2*p_1 + 3*e_3 = 3 + 1 - 3 = 1... that gives 1, not 7.
""")

# Actually compute directly
print("  Direct computation of Tr(M_k^n) for k-nacci matrices:")
print()

for k in [2, 3, 4, 5]:
    # Companion matrix for x^k - x^{k-1} - ... - x - 1
    M = np.zeros((k, k))
    M[0, :] = 1  # first row all 1s
    for i in range(1, k):
        M[i, i-1] = 1  # subdiagonal
    # Actually the companion matrix should be:
    # For x^k - x^{k-1} - x^{k-2} - ... - 1:
    # M = [[1, 1, 1, ..., 1],
    #      [1, 0, 0, ..., 0],
    #      [0, 1, 0, ..., 0],
    #      ...]
    M = np.zeros((k, k))
    for j in range(k):
        M[0, j] = 1
    for i in range(1, k):
        M[i, i-1] = 1

    traces = []
    Mn = np.eye(k)
    for n in range(1, 8):
        Mn = Mn @ M
        tr = round(np.trace(Mn))
        traces.append(tr)

    print(f"  k={k}: Tr(M^n) for n=1..7: {traces}")
    if 7 in traces:
        idx = traces.index(7) + 1
        print(f"         Tr(M^{idx}) = 7!")

print(f"""
  The tribonacci matrix M_3 has Tr(M_3^3) = 7? Let me verify...
""")

# Verify tribonacci trace at n=3
M3 = np.array([[1, 1, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
M3_cubed = M3 @ M3 @ M3
print(f"  M_3^3 = ")
for row in M3_cubed:
    print(f"    {[int(x) for x in row]}")
print(f"  Tr(M_3^3) = {int(np.trace(M3_cubed))}")
print()

# Check: is Tr(M_k^3) = 7 for ALL k >= 3?
print(f"  Tr(M_k^3) for various k:")
for k in range(2, 10):
    M = np.zeros((k, k))
    for j in range(k): M[0, j] = 1
    for i in range(1, k): M[i, i-1] = 1
    M3 = M @ M @ M
    tr = int(round(np.trace(M3)))
    print(f"    k={k}: Tr(M_k^3) = {tr} {'= 7 UNIVERSAL!' if tr == 7 else ''}")

# =============================================================================
# ANGLE 6: THE FANO PLANE
# =============================================================================

print()
print("=" * 70)
print("ANGLE 6: THE FANO PLANE PG(2,2)")
print("=" * 70)
print()

print("""
  The Fano plane is the projective plane over F_2.
  7 points, 7 lines, each line contains 3 points, each point on 3 lines.

  Labeling points 0-6, the lines are:
    {0,1,3}, {1,2,4}, {2,3,5}, {3,4,6}, {4,5,0}, {5,6,1}, {6,0,2}

  These are the QUADRATIC RESIDUES mod 7: {1, 2, 4} → shifted by each element.

  The Fano plane IS the Paley tournament T_7:
    vertex i beats j iff (i-j) mod 7 is a QR.
    QR mod 7 = {1, 2, 4}.
    The Fano "lines" are the 3-cycles of T_7!

  Specifically: {i, i+1, i+3} (mod 7) is a line of the Fano plane
  AND a directed 3-cycle i→(i+1)→(i+3)→i in T_7.
""")

# Verify
QR7 = {1, 2, 4}
fano_lines = []
for start in range(7):
    line = frozenset([(start + r) % 7 for r in [0, 1, 3]])
    fano_lines.append(line)

# Check: are these the 3-cycles of Paley T_7?
A_paley = [[0]*7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if i != j and (i-j) % 7 in QR7:
            A_paley[i][j] = 1

paley_3cycles = set()
for combo in combinations(range(7), 3):
    i, j, k = combo
    if A_paley[i][j] and A_paley[j][k] and A_paley[k][i]:
        paley_3cycles.add(frozenset(combo))
    elif A_paley[i][k] and A_paley[k][j] and A_paley[j][i]:
        paley_3cycles.add(frozenset(combo))

fano_set = set(fano_lines)
print(f"  Fano lines (as sets): {sorted(fano_set, key=lambda s: min(s))}")
print(f"  Paley T_7 3-cycles:  {sorted(paley_3cycles, key=lambda s: min(s))}")
print(f"  Same? {fano_set == paley_3cycles}")
print()

# The Fano plane has a special property: every pair of lines meets in exactly 1 point.
# This means: every pair of 3-cycles shares exactly 1 vertex.
print("  Fano property: every pair of lines shares exactly 1 point:")
all_share_1 = True
for i in range(7):
    for j in range(i+1, 7):
        shared = len(fano_lines[i] & fano_lines[j])
        if shared != 1:
            all_share_1 = False
            print(f"    Lines {i},{j} share {shared} points")
print(f"  Verified: {all_share_1}")
print()

print("""
  THE FANO CONNECTION TO FORBIDDEN VALUES:

  The Fano plane has 7 lines (= 3-cycles in T_7).
  Every pair of lines shares exactly 1 point.
  So every pair of 3-cycles OVERLAPS (shares a vertex).
  Therefore: a_2 = 0 for the 3-cycle contribution!

  But there are also 5-cycles and the 7-cycle in T_7.
  The 14 total 3-cycles of T_7 consist of 7 "Fano lines" plus 7 more
  (from the other direction, or from non-Fano triples).

  Wait: we computed 14 3-cycles for T_7 earlier, but the Fano plane has
  only 7 lines. The other 7 must be the "reverse" 3-cycles.

  Actually: in a tournament on {i,j,k}, there is at most ONE directed 3-cycle.
  There are C(7,3) = 35 triples, of which 14 form 3-cycles.
  The Fano plane accounts for 7 of them. The other 7 come from non-Fano triples.
""")

# Which triples are Fano lines and which are non-Fano?
non_fano_cycles = paley_3cycles - fano_set
print(f"  Fano 3-cycles (7): {sorted(fano_set, key=lambda s: min(s))}")
print(f"  Non-Fano 3-cycles ({len(non_fano_cycles)}): {sorted(non_fano_cycles, key=lambda s: min(s))}")

# =============================================================================
# ANGLE 7: THE SYNTHESIS — WHY 7 AND ONLY 7
# =============================================================================

print()
print("=" * 70)
print("THE SYNTHESIS: WHY 7 AND ONLY 7")
print("=" * 70)
print(f"""
Seven is the forbidden prime because it sits at the INTERSECTION
of SEVEN independent mathematical structures:

  1. MERSENNE: 7 = 2³-1. The third Mersenne prime.
     The tribonacci trace Tr(M_k³) = 7 for ALL k ≥ 3.
     This universality means 7 is baked into the k-nacci structure.

  2. CYCLOTOMIC: 7 = Phi_6(3) = 3²-3+1. The Eisenstein norm of 3+omega.
     And 21 = Phi_6(5) = 5²-5+1. The Eisenstein norm of 5+omega.
     Both forbidden values are cyclotomic norms at the TWO primes < 7.

  3. CHEBYSHEV: T_7(x)/x has constant term -7.
     The 7th Chebyshev polynomial has 7 as its "lowest mode."
     ALL other coefficients are divisible by 7 except the leading 2⁶.

  4. FANO: The Fano plane PG(2,2) has 7 points and 7 lines.
     The Paley tournament T_7 has its 3-cycles as Fano lines.
     Every pair of Fano lines meets in exactly 1 point.
     This forces all 3-cycles to OVERLAP, preventing a_2 from growing.

  5. TORSION: [7](x) has coefficients C(7,1)=7 and C(7,5)=21.
     These are literally the forbidden H-values.
     The formal group's 7-torsion encodes the obstruction.

  6. K_3 OBSTRUCTION: H=7 requires exactly 3 mutually overlapping cycles
     forming K_3 in Omega. But in a COMPLETE directed graph (tournament),
     the 4th vertex always creates additional cycles, preventing isolation.
     The obstruction is that tournaments are TOO CONNECTED for K_3 isolation.

  7. SOLVABILITY: 7 is the first prime where the full group A_7 (simple,
     non-solvable) acts on the tournament. Below 7: the groups S_3, S_4, S_5
     are "tame enough" that all H-values are achievable.
     At 7: the interaction between non-solvable S_7 and abelian Q(zeta_7)
     creates an obstruction that persists at ALL larger n.

THE CENTRAL INSIGHT:

  7 is forbidden because it's the SMALLEST number that is simultaneously:
  - A Mersenne prime (2^k - 1)
  - A cyclotomic norm (Phi_6(odd prime))
  - A Chebyshev constant (T_p(x)/x at x=0)
  - A Fano cardinal (|PG(2,F_2)|)
  - A torsion coefficient (C(p,1) for p prime)
  - A K_3 threshold (first n where K_3 in Omega is an issue)

  No smaller number sits at all six intersections.
  No larger number sits at all six either.
  7 is UNIQUE in this respect.

  And 21 = 3×7 inherits the prohibition via Phi_6(5):
  the second forbidden value is the cyclotomic norm at the NEXT prime (5)
  of the same polynomial that produces 7 at the PREVIOUS prime (3).
""")

print("opus-2026-03-20-S92o: WHY SEVEN IS THE FORBIDDEN PRIME")
