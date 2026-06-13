#!/usr/bin/env python3
"""
PROJECTIVE & ALGEBRAIC GEOMETRY — DEEP EXTENSION
opus-2026-03-14-S71n

Extending the projective/algebraic geometry investigation from S71m.
The user requests HEAVY consideration of:
- Projective geometry
- Algebraic geometry
- The unity in their duality

NEW INVESTIGATIONS:
1. The Fano matroid and tournament matroids — when is the tournament matroid
   isomorphic to PG(2,2)?
2. Veronese embedding: tournaments as points on the Veronese surface
3. Segre variety: the product structure of tournament space
4. Schubert calculus on the Grassmannian of M[a,b]
5. The tournament scheme: Spec of the H-polynomial ring
6. Blow-up at forbidden values: what happens if we "resolve" H=7?
7. Étale cohomology of the tournament variety over F_2
8. Tropical geometry of the H landscape
9. The projective dual of the H hypersurface
10. Baer subplanes and sub-tournaments
11. Algebraic curves defined by the F-polynomial
12. The moduli stack of tournaments
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import math
from fractions import Fraction
import cmath

def adj_matrix(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_hp(n, A):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def compute_F_polynomial(n, A):
    F = [0] * n
    for perm in permutations(range(n)):
        valid = True
        for k in range(n-1):
            if not A[perm[k]][perm[k+1]]:
                valid = False
                break
        if valid:
            asc = sum(1 for k in range(n-1) if perm[k] < perm[k+1])
            F[asc] += 1
    return F

def score_seq(n, A):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

def get_all_tournaments(n):
    m = n * (n-1) // 2
    results = []
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        results.append((bits, A, H))
    return results

print("=" * 70)
print("PROJECTIVE & ALGEBRAIC GEOMETRY — DEEP EXTENSION")
print("opus-2026-03-14-S71n")
print("=" * 70)

# ======================================================================
# PART 1: THE VERONESE EMBEDDING OF TOURNAMENT SPACE
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: THE VERONESE EMBEDDING")
print("=" * 70)

print("""
  The VERONESE MAP v_d: P^n → P^{C(n+d,d)-1} sends a point
  [x_0 : ... : x_n] to all degree-d monomials.

  For tournaments in F_2^m:
  The degree-2 Veronese map sends x = (x_1,...,x_m) to all
  products x_i*x_j (including x_i^2 = x_i over F_2).

  Since H is degree 2 at n=3,4, the H-level sets live on
  the Veronese surface v_2(P^{m-1}) ⊂ P^{C(m+2,2)-1}.

  This means the H hypersurface {H(x)=h} becomes a HYPERPLANE SECTION
  of the Veronese surface after lifting to the Veronese space!

  QUESTION: What are the H-level sets as intersections of
  the Veronese variety with hyperplanes?

  For n=3 (m=3): H = 1 + 2x₁ - 2x₀x₁ + 2x₀x₂ - 2x₁x₂
  This is degree 2, so H is a QUADRATIC FORM + linear + constant.

  The degree-2 Veronese of F_2^3 maps to F_2^{C(3+2,2)} = F_2^{10}?
  No — over F_2, x_i^2 = x_i, so degree-2 monomials are:
  {x_0, x_1, x_2, x_0x_1, x_0x_2, x_1x_2} plus constant 1.
  That's 7 = C(3,1) + C(3,2) + 1 = 3 + 3 + 1 monomials.
  Wait, that's 7 = |PG(2,2)|!!! The Fano plane appears again!

  Actually: the number of multilinear monomials in m variables of
  degree ≤ 2 (including constant) is 1 + m + C(m,2).
  At m=3: 1 + 3 + 3 = 7 = Φ_3(2).
  At m=6: 1 + 6 + 15 = 22.
  At m=10: 1 + 10 + 45 = 56.

  GENERAL: 1 + C(n,2) + C(C(n,2), 2) = 1 + m + m(m-1)/2.
""")

# Compute the number of degree-≤d monomials for tournament space
for n in [3, 4, 5, 6, 7]:
    m = n * (n-1) // 2
    for d in [1, 2, 3, 4]:
        count = sum(math.comb(m, k) for k in range(d+1))
        print(f"  n={n} (m={m}): degree-≤{d} monomials = {count}")
        if d == 2 and count == 7:
            print(f"    *** = 7 = |PG(2,2)| = Φ_3(2) ***")
    print()

# ======================================================================
# PART 2: THE QUADRATIC FORM STRUCTURE OF H
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: THE QUADRATIC FORM STRUCTURE OF H")
print("=" * 70)

print("""
  H at n=3 is a QUADRATIC FORM over Z (not F_2!):
  H = 1 + 2x₁ - 2x₀x₁ + 2x₀x₂ - 2x₁x₂
    = 1 + 2(x₁ - x₀x₁ + x₀x₂ - x₁x₂)
    = 1 + 2(x₁(1-x₀) + x₂(x₀-x₁))

  The QUADRATIC PART is: Q = -2x₀x₁ + 2x₀x₂ - 2x₁x₂
  The LINEAR PART is: L = 2x₁
  The CONSTANT is: 1

  So H = 1 + L + Q.

  The MATRIX of Q (the Gram matrix) is:
  Q(x) = x^T M_Q x where M_Q = [[ 0, -1,  1],
                                   [-1,  0, -1],
                                   [ 1, -1,  0]]
  (upper triangle × 2, since Q = sum_{i<j} 2a_{ij} x_i x_j)

  Wait, that's the SKEW-SYMMETRIC matrix! M_Q = skew-adjacency matrix
  of the "reference tournament" (the one with all arcs oriented 0→1→2→0).

  No — the multilinear coefficients are a_{01}=-2, a_{02}=+2, a_{12}=-2.
  These are NOT the skew-adjacency ±1 entries.

  Let me be more careful.
""")

# Compute the quadratic form matrix for n=3
n = 3
m = 3
H_at = {}
for bits in range(1 << m):
    A = adj_matrix(n, bits)
    H_at[bits] = count_hp(n, A)

# Multilinear expansion
coeffs = {}
for S in range(1 << m):
    a_S = 0
    for T in range(1 << m):
        if T & S == T:
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            a_S += sign * H_at[T]
    coeffs[S] = a_S

# Extract quadratic part
print(f"  n=3: H multilinear coefficients:")
for S in range(1 << m):
    if coeffs[S] != 0:
        S_bits = [i for i in range(m) if S & (1 << i)]
        deg = len(S_bits)
        monomial = '*'.join(f'x_{i}' for i in S_bits) if S_bits else '1'
        print(f"    deg {deg}: a = {coeffs[S]:+d}  ({monomial})")

# Build the symmetric bilinear form B(x,y) associated to Q
# Q(x) = sum_{i<j} a_{ij} x_i x_j
# B(x,y) = sum_{i<j} a_{ij} (x_i y_j + x_j y_i) / 2
# But over Z, we just record the coefficients

print("\n  Quadratic form matrix (B[i,j] = coefficient of x_i*x_j):")
B = [[0]*m for _ in range(m)]
for S, a in coeffs.items():
    bits = [i for i in range(m) if S & (1 << i)]
    if len(bits) == 2:
        i, j = bits
        B[i][j] = a
        B[j][i] = a

for row in B:
    print(f"    {row}")

# Determinant of B
det_B = (B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1])
        -B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0])
        +B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]))
print(f"  det(B) = {det_B}")

# The discriminant of the quadratic form
# For a binary quadratic form over Z: discriminant = -4*det(B)?
# Actually for 3 variables: discriminant is more complex
print(f"  Rank of B over Q: ", end="")
if det_B != 0:
    print(f"3 (full rank)")
else:
    # Check 2x2 minors
    rank = 0
    for i in range(m):
        for j in range(m):
            if B[i][j] != 0:
                rank = max(rank, 1)
    for i1, i2 in combinations(range(m), 2):
        for j1, j2 in combinations(range(m), 2):
            if B[i1][j1]*B[i2][j2] - B[i1][j2]*B[i2][j1] != 0:
                rank = 2
    print(f"{rank}")

# Same for n=4
n = 4
m = 6
H_at_4 = {}
for bits in range(1 << m):
    A = adj_matrix(n, bits)
    H_at_4[bits] = count_hp(n, A)

coeffs_4 = {}
for S in range(1 << m):
    a_S = 0
    for T in range(1 << m):
        if T & S == T:
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            a_S += sign * H_at_4[T]
    coeffs_4[S] = a_S

print(f"\n  n=4: Quadratic form matrix (6×6):")
B4 = [[0]*m for _ in range(m)]
for S, a in coeffs_4.items():
    bits = [i for i in range(m) if S & (1 << i)]
    if len(bits) == 2:
        i, j = bits
        B4[i][j] = a
        B4[j][i] = a

for row in B4:
    print(f"    {row}")

# Show the arc labels
print(f"\n  Arc labels for n=4:")
idx = 0
for i in range(4):
    for j in range(i+1, 4):
        print(f"    x_{idx} = arc ({i},{j})")
        idx += 1

# ======================================================================
# PART 3: THE F-POLYNOMIAL CURVE IN PROJECTIVE SPACE
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: THE F-POLYNOMIAL AS AN ALGEBRAIC CURVE")
print("=" * 70)

print("""
  The F-polynomial F(T, q) = sum_k F_k(T) q^k is a polynomial in q
  with nonneg integer coefficients summing to H(T).

  For fixed T, F defines a POINT in P^{n-1} (projective space):
  [F_0 : F_1 : ... : F_{n-1}]

  The SET of all possible F-polynomial vectors {F(T) : T tournament on n}
  is a FINITE subset of P^{n-1}.

  QUESTION: What is the CONVEX HULL of these points?
  What is the ALGEBRAIC VARIETY they span?

  For n=5: F lives in P^4. The achievable F-vectors are a finite set.
  The SECANT VARIETY connects pairs of F-vectors.
  The TANGENT SPACE at a point [F] tells us how F changes under arc flips.
""")

n = 5
tournaments = get_all_tournaments(n)

# Compute all F-polynomials
F_vecs = defaultdict(list)
for bits, A, H in tournaments:
    F = compute_F_polynomial(n, A)
    F_vecs[tuple(F)].append(bits)

print(f"\n  n={n}: {len(F_vecs)} distinct F-polynomial vectors")

# Show the vectors
for F_tuple in sorted(F_vecs.keys()):
    count = len(F_vecs[F_tuple])
    H_val = sum(F_tuple)
    print(f"    F = {list(F_tuple)}, H={H_val}, count={count}")

# Check: do the F-vectors lie on a linear subspace?
# Reduce over Q (just check rank)
F_list = list(F_vecs.keys())
# Center by subtracting mean
# Actually just check if they span all of R^n or a proper subspace

# F_0 + F_1 + ... + F_{n-1} = H (always)
# So all F-vectors lie on the hyperplane sum = H.
# But H varies, so they DON'T all lie on one hyperplane.

# What RELATIONS do they satisfy?
print(f"\n  Linear relations among F-vector components:")
print(f"    F_0 + F_1 + ... + F_{n-1} = H (always)")

# Check: F_0 = F_{n-1}? (complement duality)
# F_k(T) = F_{n-1-k}(T^op), so if T = T^op (self-complementary),
# then F_k = F_{n-1-k} (palindromic)
palindromic_count = 0
for F_tuple in F_vecs.keys():
    F = list(F_tuple)
    if F == F[::-1]:
        palindromic_count += 1
print(f"    Palindromic F (F_k = F_{{n-1-k}}): {palindromic_count}/{len(F_vecs)}")

# F evaluated at various points
print(f"\n  F evaluated at special points:")
for F_tuple in sorted(F_vecs.keys()):
    F = list(F_tuple)
    H = sum(F)
    F_neg1 = sum(F[k] * (-1)**k for k in range(n))

    omega = cmath.exp(2j * cmath.pi / 3)
    F_omega = sum(F[k] * omega**k for k in range(n))
    norm_omega = abs(F_omega)**2

    # F at golden ratio
    phi = (1 + 5**0.5) / 2
    F_phi = sum(F[k] * phi**k for k in range(n))

    print(f"    F={F}: F(1)={H}, F(-1)={F_neg1}, "
          f"|F(ω)|²={norm_omega:.0f}, F(φ)={F_phi:.2f}")

# ======================================================================
# PART 4: SEGRE PRODUCT AND TOURNAMENT DECOMPOSITION
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: SEGRE PRODUCT — COMPOSITION OF TOURNAMENTS")
print("=" * 70)

print("""
  The SEGRE EMBEDDING σ: P^a × P^b → P^{(a+1)(b+1)-1} realizes the
  product of two projective spaces as a projective variety.

  For tournaments: if T₁ is a tournament on n₁ and T₂ on n₂,
  can we define a PRODUCT T₁ × T₂?

  COMPOSITION (Lexicographic product):
  T₁[T₂] has n₁*n₂ vertices. Vertex (i,j) beats (k,l) iff:
  - T₁: i→k, OR
  - i=k and T₂: j→l.

  H(T₁[T₂]) = H(T₁) * H(T₂)^{n₁}?
  No — it's more complex. The HP must traverse all n₁*n₂ vertices.

  SIMPLER: Direct sum T₁ ⊕ T₂ (disjoint union, all arcs from T₁ to T₂).
  This is a tournament on n₁+n₂ vertices.
  H(T₁ ⊕ T₂) depends on how many HPs cross between T₁ and T₂.

  For the ORDERED direct sum (all arcs from V₁ to V₂):
  Any HP must visit all of V₁ before entering V₂ (or vice versa)?
  No — it can interleave.

  Let's compute: for small T₁, T₂, what is H(T₁ ⊕ T₂)?
""")

# Compute direct sum of two tournaments
def direct_sum(n1, A1, n2, A2, direction='forward'):
    """Direct sum: tournament on n1+n2 with all arcs from V1 to V2."""
    n = n1 + n2
    A = [[0]*n for _ in range(n)]
    # Copy A1
    for i in range(n1):
        for j in range(n1):
            A[i][j] = A1[i][j]
    # Copy A2
    for i in range(n2):
        for j in range(n2):
            A[n1+i][n1+j] = A2[i][j]
    # Cross arcs: all from V1 to V2 (forward) or V2 to V1 (backward)
    if direction == 'forward':
        for i in range(n1):
            for j in range(n2):
                A[i][n1+j] = 1
    else:
        for i in range(n2):
            for j in range(n1):
                A[n1+i][j] = 1
    return A

# All 3-tournaments
t3 = get_all_tournaments(3)
H_3 = {bits: H for bits, _, H in t3}

# Direct sum of two 3-tournaments
print(f"\n  Direct sum T₁ ⊕ T₂ (n₁=n₂=3, all arcs V₁→V₂):")
print(f"  {'H(T₁)':>6s} {'H(T₂)':>6s} {'H(T₁⊕T₂)':>10s} {'Product':>8s} {'Ratio':>8s}")

for bits1, A1, H1 in t3:
    if bits1 > 0 and bits1 != 5:  # Pick one from each H class
        continue
    for bits2, A2, H2 in t3:
        if bits2 > 0 and bits2 != 5:
            continue
        A_sum = direct_sum(3, A1, 3, A2)
        H_sum = count_hp(6, A_sum)
        product = H1 * H2
        ratio = H_sum / product if product > 0 else float('inf')
        print(f"  {H1:6d} {H2:6d} {H_sum:10d} {product:8d} {ratio:8.2f}")

# ======================================================================
# PART 5: THE TOURNAMENT MATROID
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: THE TOURNAMENT MATROID AND FANO CONNECTION")
print("=" * 70)

print("""
  A MATROID M(E, I) on ground set E with independent sets I
  captures the notion of "linear independence" abstractly.

  For tournaments: consider the CYCLE MATROID of the underlying
  undirected graph (always K_n), which is the graphic matroid M(K_n).

  M(K_n) has rank n-1 and is the same for ALL tournaments on n vertices.
  So this doesn't distinguish tournaments.

  MORE INTERESTING: The ORIENTED MATROID of T.
  The signed incidence matrix of T (±1 entries) defines an oriented matroid.
  The chirotope (sign pattern of maximal minors) encodes the orientation.

  For n=3: The signed incidence matrix is 3×3 (skew-symmetric).
  Its chirotope has 1 element (the Pfaffian sign... but n=3 is odd).

  Actually: the oriented matroid of a tournament T on n vertices
  is determined by the SIGNS of all (n-1)×(n-1) minors of the
  signed incidence matrix.

  FANO MATROID: The Fano matroid F_7 is the matroid of PG(2,2).
  It has 7 elements, rank 3, and is the smallest non-representable
  matroid over R (but IS representable over F_2).

  QUESTION: Does the Fano matroid appear in tournament theory?

  At n=3, m=3 arcs. The cycle space of K_3 has dimension C(3,2)-3+1=1.
  The Fano matroid has 7 elements, so it can't appear at n=3.

  At n=4, m=6 arcs. The graphic matroid M(K_4) has rank 3.
  The Fano matroid also has rank 3, on 7 elements.
  But m=6 ≠ 7, so no direct match.

  HOWEVER: if we ADD the H function as a 7th "element" to the
  matroid on m=6 arcs, we get a matroid on 7 elements of rank 3.
  Could THIS be the Fano matroid?

  The Fano matroid's circuits (minimal dependent sets) are:
  {1,2,4}, {2,3,5}, {3,4,6}, {4,5,7}, {5,6,1}, {6,7,2}, {7,1,3}
  (the 7 lines of PG(2,2)).

  For K_4: the graphic matroid circuits are the cycles of K_4:
  {01,12,02}, {01,13,03}, {02,23,03}, {12,23,13},
  plus the whole graph {01,02,03,12,13,23}.

  These are 3-element circuits (triangles) and the 6-element circuit.
  The Fano matroid has only 3-element circuits. So M(K_4) is NOT F_7.

  But the DUAL of M(K_4) has rank m-n+1 = 6-3 = 3, same as F_7.
  The dual matroid M*(K_4) = the bond matroid of K_4.
  Its circuits are the minimal edge cuts of K_4.
  Min cuts of K_4: removing all edges incident to one vertex (3 edges each).
  There are 4 such cuts, each of size 3.

  Circuits of M*(K_4): {01,02,03}, {01,12,13}, {02,12,23}, {03,13,23}.
  Four circuits of size 3 — same as the graphic matroid of K_4
  (which is self-dual for K_4? No, K_4 is planar so M*(K_4) = M(K_4*)
  where K_4* is the dual graph = K_4 again... actually that gives the
  same matroid since K_4 is self-dual as a planar graph!).

  Nope: M(K_4) is uniform U_{3,6} only if K_4 were "generic," but
  K_4 has structure. Actually M(K_4) is isomorphic to the graphic
  matroid of K_4, which has rank 3 on 6 elements with the 4 triangle circuits.
""")

# Let's just verify the structure
print("  Cycle matroid of K_4 (m=6 edges, rank 3):")
print("  Edges: 01, 02, 03, 12, 13, 23")
print("  Circuits (3-cycles): {01,02,12}, {01,03,13}, {02,03,23}, {12,13,23}")
print("  4-cycles: {01,12,23,03}, {01,13,23,02}, {02,12,13,03}")
print()

# But the KEY connection to PG(2,2) is different:
# PG(2,2) has 7 POINTS and 7 LINES.
# The points of PG(2,2) can be labeled by nonzero vectors in F_2^3:
# {001, 010, 011, 100, 101, 110, 111}
# The lines are the 2-dimensional subspaces:
# Each line has 3 points, and each pair of points determines a unique line.

print("  PG(2,2) = Fano plane:")
print("  Points = nonzero vectors in F_2^3:")
points = [(i,j,k) for i in range(2) for j in range(2) for k in range(2) if (i,j,k) != (0,0,0)]
for i, p in enumerate(points):
    print(f"    P_{i+1} = {p}")

print("  Lines (3-point subsets where sum = 0 mod 2):")
lines = []
for triple in combinations(range(7), 3):
    p1, p2, p3 = [points[i] for i in triple]
    if all((p1[k] + p2[k] + p3[k]) % 2 == 0 for k in range(3)):
        lines.append(triple)
        pt_names = [f"P_{i+1}" for i in triple]
        print(f"    {{{', '.join(pt_names)}}}")

# ======================================================================
# PART 6: TOURNAMENT SPACE AS TORIC VARIETY
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: TORIC GEOMETRY OF TOURNAMENT SPACE")
print("=" * 70)

print("""
  The tournament hypercube Q_m = {0,1}^m is the set of LATTICE POINTS
  of the m-dimensional unit cube [0,1]^m.

  In TORIC GEOMETRY, the unit cube corresponds to the projective space P^m.
  More precisely:
  - The vertices {0,1}^m are the torus-fixed points of P^m.
  - The faces of the cube correspond to the torus orbits.
  - The H function is a "toric polynomial" — a Laurent polynomial
    on the torus (C*)^m.

  Actually, since H is a multilinear polynomial, it IS a section of
  a line bundle on the toric variety X_Cube associated to the cube.

  The toric variety X = TV(Δ) where Δ is the normal fan of [0,1]^m.
  For the cube, Δ is the fan over the faces of the CROSS-POLYTOPE.
  And X_Cube ≅ (P^1)^m = the product of m copies of P^1.

  Each tournament is a point on this toric variety!
  And H is a section of O(1,1,...,1) (the "multilinear" line bundle).

  For n=3: X = (P^1)^3. H is a section of O(1,1,1).
  sections of O(1,1,1) on (P^1)^3 form a space of dimension 2^3 = 8.
  H is one specific section.

  The ZERO LOCUS of H - h (as a section of O(1,1,1)) is a divisor in X.
  Over C: this divisor has degree (1,1,1) = total degree 3 in the
  Segre embedding of (P^1)^3 → P^7.

  THIS IS THE PROJECTIVE-ALGEBRAIC DUALITY:
  - PROJECTIVELY: H-level sets are hyperplane sections of the Segre variety
  - ALGEBRAICALLY: H is a polynomial function on affine space
  - These UNIFY because the Segre embedding realizes (P^1)^m as a
    projective variety, and multilinear functions become linear on the
    ambient projective space

  THE SEGRE-VERONESE DUALITY:
  - Segre: product P^1 × ... × P^1 → P^{2^m - 1} (multilinear)
  - Veronese: P^m → P^{C(m+d,d)-1} (degree-d homogeneous)
  - For d=1: Veronese = identity, Segre gives the full structure
  - H being degree 2 means: H factors through the degree-2 Veronese
    of the linear span, OR equivalently, H is LINEAR on the Segre variety
""")

# Verify: count of sections of O(1,...,1) on (P^1)^m
for n in [3, 4, 5]:
    m = n * (n-1) // 2
    sections = 2**m
    print(f"  n={n} (m={m}): dim H^0((P^1)^{m}, O(1,...,1)) = {sections}")
    print(f"    = total functions {0,1}^{m} → Z  (multilinear)")
    print(f"    H is ONE such section")

# ======================================================================
# PART 7: BAER SUBPLANES AND SUB-TOURNAMENTS
# ======================================================================
print("\n" + "=" * 70)
print("PART 7: BAER SUBPLANES AND SUB-TOURNAMENTS")
print("=" * 70)

print("""
  A BAER SUBPLANE of PG(2,q²) is a copy of PG(2,q) embedded within it.
  |PG(2,q²)| = q⁴ + q² + 1
  |PG(2,q)| = q² + q + 1
  Ratio: (q⁴+q²+1)/(q²+q+1) = q²-q+1 = Φ_6(q)

  At q=2: PG(2,4) has 21 points, PG(2,2) has 7 points.
  Ratio = 21/7 = 3 = Φ_6(2) = 2² - 2 + 1.

  TOURNAMENT CONNECTION:
  The FORBIDDEN VALUES are 7 = |PG(2,2)| and 21 = |PG(2,4)|.
  A Baer subplane of PG(2,4) is a copy of PG(2,2).
  So PG(2,4) contains Φ_6(2) = 3 copies of PG(2,2)? No — each point
  is in exactly one Baer subplane of the Baer partition.

  PG(2,4) has exactly q²-q+1 = 3 disjoint Baer subplanes that
  partition all the points when q=2: 21 = 3 × 7.

  THIS IS 21 = 3 × 7 = Φ_6(2) × Φ_3(2).
  The factorization of the second forbidden value involves the first!

  For SUB-TOURNAMENTS: a tournament on k vertices embedded in
  a tournament on n vertices (induced sub-tournament).
  The number of sub-tournaments of size k is C(n, k).

  QUESTION: Is there a "Baer sub-tournament" concept?
  A sub-tournament T' ⊂ T such that H(T') relates to H(T)
  in a way analogous to Baer subplanes?

  For the transitive tournament T_n (H=1):
  Every sub-tournament is transitive: H(T'_k) = 1 for all k ≤ n.
  Not interesting.

  For the regular tournament C_5 (H=15):
  Sub-tournaments of size 3: C(5,3) = 10 induced sub-tournaments.
  How many have H=1 and how many have H=3?
""")

# Compute sub-tournament structure for n=5
n = 5
t5 = get_all_tournaments(n)

# Pick specific tournaments
trans_5 = None
reg_5 = None
for bits, A, H in t5:
    if H == 1 and trans_5 is None:
        trans_5 = (bits, A, H)
    if H == 15:
        ss = score_seq(n, A)
        if ss == (2,2,2,2,2) and reg_5 is None:
            reg_5 = (bits, A, H)

for name, (bits, A, H) in [("Transitive (H=1)", trans_5), ("Regular (H=15)", reg_5)]:
    print(f"\n  {name}:")
    print(f"    Score: {score_seq(n, A)}")

    # All induced 3-sub-tournaments
    sub_H = Counter()
    for triple in combinations(range(n), 3):
        # Extract induced tournament
        sub_A = [[0]*3 for _ in range(3)]
        for i, vi in enumerate(triple):
            for j, vj in enumerate(triple):
                sub_A[i][j] = A[vi][vj]
        sub_h = count_hp(3, sub_A)
        sub_H[sub_h] += 1

    print(f"    Induced 3-sub-tournaments: {dict(sorted(sub_H.items()))}")

    # 4-sub-tournaments
    sub_H4 = Counter()
    for quad in combinations(range(n), 4):
        sub_A = [[0]*4 for _ in range(4)]
        for i, vi in enumerate(quad):
            for j, vj in enumerate(quad):
                sub_A[i][j] = A[vi][vj]
        sub_h = count_hp(4, sub_A)
        sub_H4[sub_h] += 1

    print(f"    Induced 4-sub-tournaments: {dict(sorted(sub_H4.items()))}")

# ======================================================================
# PART 8: THE HILBERT POLYNOMIAL OF THE H-SCHEME
# ======================================================================
print("\n" + "=" * 70)
print("PART 8: THE SCHEME Spec(Z[H]) AND ITS FIBERS")
print("=" * 70)

print("""
  Consider the map H: {0,1}^m → Z.
  This defines a SCHEME morphism π: Spec(Z[x_1,...,x_m]/(x_i²-x_i)) → Spec(Z[H]).

  The fiber π^{-1}(h) is the tournament variety {T : H(T) = h}.
  This is a 0-dimensional scheme over Z.

  The HILBERT FUNCTION of this fiber:
  h_F(d) = #{monomials of degree ≤ d constant on the fiber}

  But for us, the more interesting question is:
  What is the GENUS of the curve {H(x) = h} in (P^1)^m?

  For n=3, H has degree (1,1,1) in (P^1)^3.
  By adjunction formula: the generic fiber has genus
  g = (d₁-1)(d₂-1) + ... = 0 (for a trilinear hypersurface in (P^1)^3).
  Actually for a smooth divisor D of class (1,1,1) in (P^1)^3:
  D is isomorphic to a del Pezzo surface of degree 6.

  Wait — (P^1)^3 is 3-dimensional, and a divisor D of class (1,1,1)
  is a surface (2-dimensional), not a curve.

  The CURVE arises when we intersect D with another divisor.
  E.g., {H = h₁} ∩ {H = h₂} is empty (since H takes one value at each point).

  Better: the H-level set is a FINITE set of points in (P^1)^3.
  These points form a 0-dimensional subscheme.
  Its DEGREE = |{T : H(T) = h}| = multiplicity of h in the H-spectrum.
""")

# Compute degrees (= multiplicities) of H-level sets
for n in [3, 4, 5]:
    m = n * (n-1) // 2
    t = get_all_tournaments(n)
    H_counter = Counter(H for _, _, H in t)

    print(f"\n  n={n} (m={m}):")
    for h in sorted(H_counter.keys()):
        count = H_counter[h]
        # As a scheme: degree = count / 2^m of the full space?
        # No: degree = count (number of points in the fiber)
        print(f"    H={h:3d}: degree = {count}, "
              f"frac of (P^1)^{m} = {count / 2**m:.6f}")

# ======================================================================
# PART 9: THE ALGEBRAIC-PROJECTIVE DUALITY MADE EXPLICIT
# ======================================================================
print("\n" + "=" * 70)
print("PART 9: THE DUALITY MADE EXPLICIT")
print("=" * 70)

print("""
  THE CENTRAL THEOREM (of this investigation):

  Let T_n be the space of tournaments on n vertices.
  T_n ≅ F_2^m as an affine space (m = C(n,2)).
  T_n ≅ (P^1)^m as a toric variety.
  T_n ↪ P^{2^m-1} via the Segre embedding.

  H: T_n → Z is a multilinear function of degree ≤ n-1.

  PROJECTIVE VIEWPOINT:
  The Segre embedding maps (P^1)^m → P^{2^m-1}.
  H becomes a LINEAR function on P^{2^m-1}.
  The H-level sets are HYPERPLANE SECTIONS of the Segre variety.
  Projective duality exchanges H-level sets with their dual hyperplanes.

  ALGEBRAIC VIEWPOINT:
  H is a polynomial in F_2[x_1,...,x_m]/(x_i²-x_i).
  The H-level sets are varieties V(H-h) in affine space.
  Their ideal I(V) lives in the Boolean polynomial ring.
  The Hilbert function of I(V) encodes the algebraic complexity.

  THE UNIFICATION:
  The Segre embedding realizes T_n as a projective variety.
  Multilinear polynomials become linear forms on the ambient space.
  H-level sets become hyperplane sections.
  Projective duality exchanges points ↔ hyperplanes.
  Algebraic duality exchanges varieties ↔ ideals.

  THESE ARE THE SAME DUALITY: the hyperplane {H = h} in P^{2^m-1}
  defines both a projective variety (the section) and an algebraic ideal
  (the kernel of the linear form H-h).

  The FANO PLANE enters because:
  - At n=3 (m=3), the Segre variety Σ ⊂ P^7 has 2^3 = 8 points
  - The number of degree-≤2 multilinear monomials is 7 = |PG(2,2)|
  - H is one of these 7 "directions" in the multilinear space
  - The forbidden value 7 = Φ_3(2) is the COUNT of these directions
  - The 7 Fano lines correspond to the 7 degree-2 equations that
    constrain the tournament structure

  THE BRIDGE BETWEEN PROJECTIVE AND ALGEBRAIC:
  In classical geometry, "projective" = "linear algebra over a field"
  and "algebraic" = "polynomial equations."
  The Segre/Veronese embeddings LINEARIZE polynomial equations,
  making algebraic geometry into projective geometry in a higher-dimensional space.

  For tournaments: H is algebraic (polynomial) in T_n but PROJECTIVE (linear)
  in the Segre-embedded P^{2^m-1}. This is the unity of the duality.
""")

# Final computation: dimension of the linear span of each H-level set
# in the Segre-embedded space
for n in [3, 4, 5]:
    m = n * (n-1) // 2
    t = get_all_tournaments(n)
    H_set = sorted(set(H for _, _, H in t))

    print(f"\n  n={n}: H-level sets in Segre-embedded space P^{2**m - 1}:")

    for h in H_set:
        level = [bits for bits, _, H in t if H == h]

        # In the Segre embedding, each tournament bits maps to a point
        # whose coordinates are indexed by subsets S ⊆ {0,...,m-1}:
        # coord_S(bits) = product_{i ∈ S} x_i = 1 if S ⊆ supp(bits), else 0
        # Actually: coord_S(bits) = prod_{i∈S} bits_i = 1 iff all bits in S are 1

        # Compute the Segre vectors
        segre_vecs = []
        for bits in level:
            vec = []
            for S in range(1 << m):
                if (bits & S) == S:
                    vec.append(1)
                else:
                    vec.append(0)
            segre_vecs.append(vec)

        # Compute rank over F_2
        # (Gaussian elimination on the binary matrix)
        mat = [list(v) for v in segre_vecs]
        rows = len(mat)
        cols = len(mat[0])
        rank = 0
        pivot_col = 0
        for r in range(rows):
            if pivot_col >= cols:
                break
            # Find pivot
            found = False
            for rr in range(r, rows):
                if mat[rr][pivot_col] == 1:
                    mat[r], mat[rr] = mat[rr], mat[r]
                    found = True
                    break
            if not found:
                pivot_col += 1
                continue
            # Eliminate
            for rr in range(rows):
                if rr != r and mat[rr][pivot_col] == 1:
                    for cc in range(cols):
                        mat[rr][cc] ^= mat[r][cc]
            rank += 1
            pivot_col += 1

        print(f"    H={h:3d}: |level|={len(level):5d}, "
              f"Segre rank (F_2) = {rank}/{2**m}")

print("\n" + "=" * 70)
print("DONE — PROJECTIVE/ALGEBRAIC GEOMETRY DEEP EXTENSION")
print("=" * 70)
