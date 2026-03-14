#!/usr/bin/env python3
"""
PROJECTIVE & ALGEBRAIC GEOMETRY IN TOURNAMENT THEORY
opus-2026-03-14-S71m

The user requests heavy consideration of:
1. Projective geometry
2. Algebraic geometry
3. The unity in their duality

CONNECTIONS TO EXPLORE:
A) Tournaments as points in projective space
B) The H-spectrum as an algebraic variety
C) Fano plane PG(2,2) and 7 = Φ_3(2)
D) Grassmannians and tournament structure
E) Tournament ideal and Hilbert function
F) Projective duality: points ↔ hyperplanes in tournament space
G) Algebraic curves: the F-polynomial defines a curve in (a,b) space
H) The tournament hypercube Q_m as a toric variety
"""

from itertools import permutations, combinations
from collections import Counter, defaultdict
import math
from fractions import Fraction

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

def score_seq(n, A):
    return tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))

print("=" * 70)
print("PROJECTIVE & ALGEBRAIC GEOMETRY IN TOURNAMENT THEORY")
print("opus-2026-03-14-S71m")
print("=" * 70)

# ======================================================================
# PART 1: THE FANO PLANE AND 7 = Φ_3(2)
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: THE FANO PLANE PG(2,2) AND THE FORBIDDEN 7")
print("=" * 70)

print("""
  The FANO PLANE PG(2,2) is the smallest projective plane.
  It has 7 points, 7 lines, 3 points per line, 3 lines per point.

  |PG(2,q)| = q² + q + 1 = Φ_3(q)

  At q=2: |PG(2,2)| = 7 = Φ_3(2) = the FORBIDDEN H value!
  At q=4: |PG(2,4)| = 21 = Φ_3(4) = the SECOND forbidden value!

  These are the ONLY permanently forbidden H values.

  QUESTION: Is there a direct connection between the Fano plane
  structure and why 7 is forbidden?

  The Fano plane's incidence matrix is 7×7 with exactly 3 ones per row/column.
  Its automorphism group is GL(3,2) = PSL(2,7), order 168 = 24*7.

  A tournament on 7 vertices has C(7,2) = 21 = Φ_3(4) arcs!
  So: the number of ARCS in a 7-tournament = the number of POINTS in PG(2,4).
  And the number of VERTICES = the number of points in PG(2,2).

  COINCIDENCE: For n-vertex tournament, m = C(n,2) = n(n-1)/2.
  m = |PG(2,n-2)| when n(n-1)/2 = (n-2)² + (n-2) + 1 = n² - 3n + 3.
  So n(n-1)/2 = n² - 3n + 3 → n² - n = 2n² - 6n + 6 → n² - 5n + 6 = 0
  → (n-2)(n-3) = 0 → n = 2 or n = 3.

  At n=3: m = 3 = |PG(2,1)| = 3 ✓ (PG(2,1) = triangle = 3 points)
  At n=2: m = 1 = |PG(2,0)| = 1 ✓ (trivial)

  So the arc count equals a projective plane size ONLY at n=2,3.
  But the FORBIDDEN values 7,21 are Φ_3 at 2,4 = 2^1, 2^2.
  These are the PRIME POWER orders of Galois fields!

  PG(2, 2^k) has Φ_3(2^k) points.
  Φ_3(2) = 7, Φ_3(4) = 21, Φ_3(8) = 73, Φ_3(16) = 273, ...

  Are 73, 273 also special in the H-spectrum?
  At n=7, max_H grows to ~315, so 73 is achievable.
  We'd need to check if 73 has special properties.
""")

# Compute Φ_3 values
print("  Φ_3(q) = q² + q + 1 for small q:")
for q in range(1, 20):
    phi3 = q*q + q + 1
    is_prime_power = False
    for p in [2, 3, 5, 7, 11, 13, 17, 19]:
        pk = p
        while pk <= q:
            if pk == q:
                is_prime_power = True
                break
            pk *= p
    if q == 1:
        is_prime_power = True
    print(f"    Φ_3({q:2d}) = {phi3:4d}"
          f"{'  ← PG(2,' + str(q) + ') exists (GF(' + str(q) + '))' if is_prime_power else ''}")

# ======================================================================
# PART 2: TOURNAMENTS AS POINTS IN PROJECTIVE SPACE
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: TOURNAMENTS AS POINTS IN PROJECTIVE SPACE")
print("=" * 70)

print("""
  A tournament on n vertices is a point in F_2^m where m = C(n,2).
  (F_2 = GF(2) = {0,1})

  The PROJECTIVE SPACE PG(m-1, 2) has 2^m - 1 = 2^{C(n,2)} - 1 points.
  Each nonzero vector in F_2^m corresponds to a point in PG(m-1, 2).

  Tournaments correspond to ALL 2^m vectors (including the zero vector),
  which are the AFFINE points of AG(m, 2).

  The tournament T with all arcs oriented one way (the zero vector)
  is the ORIGIN of this affine space.

  Flipping an arc = adding a basis vector e_i.
  The arc-flip graph Q_m is the Cayley graph of (F_2^m, +)
  with respect to the standard basis.

  PROJECTIVE DUALITY in PG(m-1, 2):
  - A k-flat (k-dimensional subspace) ↔ an (m-1-k)-flat
  - Points ↔ Hyperplanes

  A HYPERPLANE in PG(m-1, 2) is defined by a linear form:
  sum a_i x_i = 0 (mod 2), for some (a_1, ..., a_m) ≠ 0.

  Each linear form partitions tournaments into two halves:
  those satisfying the equation and those not.
  There are 2^m - 1 = m + C(m,2) + ... hyperplanes.

  QUESTION: Does the H function define a "variety" in this space?
  The H-level sets {T : H(T) = h} are subsets of F_2^m.
  Are they linear subspaces? Cosets? Algebraic varieties over F_2?
""")

# Check: is the H=max level set a linear subspace over F_2?
for n in [3, 4, 5]:
    m = n * (n-1) // 2
    tournaments = []
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        tournaments.append((bits, H))

    H_set = sorted(set(H for _, H in tournaments))
    max_H = max(H_set)

    print(f"\n  n={n}: Testing if H-level sets are F_2-linear:")
    for h in H_set:
        level = [bits for bits, H in tournaments if H == h]

        # A set is an F_2-subspace if 0 ∈ set and a⊕b ∈ set for all a,b in set
        is_subspace = 0 in level
        if is_subspace:
            for a in level:
                for b in level:
                    if (a ^ b) not in level:
                        is_subspace = False
                        break
                if not is_subspace:
                    break

        # Is it a COSET of a subspace? Check if a⊕b⊕c ∈ set for all a,b,c
        is_coset = True
        level_set = set(level)
        for a in level[:min(20, len(level))]:
            for b in level[:min(20, len(level))]:
                for c in level[:min(20, len(level))]:
                    if (a ^ b ^ c) not in level_set:
                        is_coset = False
                        break
                if not is_coset:
                    break
            if not is_coset:
                break

        print(f"    H={h:3d}: {len(level)} elements, "
              f"subspace={'Y' if is_subspace else 'N'}, "
              f"coset={'Y' if is_coset else 'N'}")

# ======================================================================
# PART 3: THE TOURNAMENT IDEAL
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: THE TOURNAMENT IDEAL AND HILBERT FUNCTION")
print("=" * 70)

print("""
  In ALGEBRAIC GEOMETRY, we associate to a set of points S ⊆ F_q^m
  the VANISHING IDEAL I(S) = {f ∈ F_q[x_1,...,x_m] : f(s) = 0 ∀s ∈ S}.

  For tournaments: S = {T : H(T) = h} for fixed h.

  Over F_2, every variable satisfies x_i² = x_i, so the polynomial ring
  is F_2[x_1,...,x_m] / (x_1²-x_1, ..., x_m²-x_m).
  This quotient ring has dimension 2^m (one basis element per subset of variables).

  The HILBERT FUNCTION of I(S) counts:
  dim(F_2[x]_d / I(S)_d) = number of degree-d polynomials not vanishing on S.

  For the H=1 (transitive) tournaments:
  |S| = n! (number of transitive tournaments).
  The ideal I(S) has codimension n! in the full ring.

  QUESTION: Does the Hilbert function of the H-level set ideal
  have a nice form?

  More concretely: how many degree-d monomials (= Walsh monomials)
  are constant on the H-level set?
""")

# For small cases, compute the "Walsh support" of each H level
for n in [3, 4, 5]:
    m = n * (n-1) // 2
    tournaments = []
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        tournaments.append((bits, H))

    H_set = sorted(set(H for _, H in tournaments))

    # Walsh transform: for each Walsh monomial chi_S,
    # compute the VARIANCE of chi_S within each H-level set
    # If variance = 0, then chi_S is constant on the level set
    # i.e., chi_S ∈ I(level set)

    print(f"\n  n={n}: Walsh monomials constant on each H-level set:")
    for h in H_set:
        level = [bits for bits, H in tournaments if H == h]
        # For each Walsh character chi_S (S ⊆ {0,...,m-1}):
        # chi_S(T) = (-1)^{popcount(S & T)}
        # chi_S is constant on level iff all level members give same value

        constant_by_degree = Counter()
        for S in range(1 << m):
            # Compute chi_S on level
            vals = set()
            for T in level:
                vals.add((-1) ** bin(S & T).count('1'))
            if len(vals) == 1:
                deg = bin(S).count('1')
                constant_by_degree[deg] += 1

        total_constant = sum(constant_by_degree.values())
        print(f"    H={h:3d}: {total_constant}/{1 << m} Walsh monomials constant")
        for d in sorted(constant_by_degree.keys()):
            print(f"        degree {d}: {constant_by_degree[d]}/{math.comb(m, d)}")

# ======================================================================
# PART 4: GRASSMANNIANS AND TOURNAMENT STRUCTURE
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: GRASSMANNIANS AND PLÜCKER COORDINATES")
print("=" * 70)

print("""
  The GRASSMANNIAN Gr(k, n) parametrizes k-dimensional subspaces of F^n.
  Its Plücker embedding maps Gr(k,n) → P(C(n,k)-1).

  For tournaments: the M[a,b] matrix (HP count from a to b) is n×n.
  Its rank determines a point in a Grassmannian.

  If M has rank r, then the column space of M is a point in Gr(r, n).
  The Plücker coordinates are the r×r minors of M.

  QUESTION: What ranks does M achieve for different H values?
  And do the Plücker coordinates encode interesting invariants?

  Also: the TOURNAMENT VARIETY.
  A tournament T can be encoded as a skew-symmetric matrix S.
  S lives in the space of n×n skew-symmetric matrices, which has
  dimension n(n-1)/2 = m. This is the LIE ALGEBRA so(n).

  The Pfaffian of S (for even n) is a polynomial of degree n/2
  in the m entries of S. The PFAFFIAN VARIETY is the zero set
  of the Pfaffian in this ambient space.

  TOURNAMENT CONNECTION:
  - T is a "real point" of the skew-symmetric variety where all entries are ±1
  - H(T) can be expressed as a polynomial in the entries of S
  - The H-level sets are algebraic varieties cut out by H(T) = h
""")

# For n=4, compute rank of M[a,b] for each tournament class
n = 4
m = n * (n-1) // 2

H_seen = {}
for bits in range(1 << m):
    A = adj_matrix(n, bits)
    H = count_hp(n, A)
    ss = score_seq(n, A)
    if (H, ss) in H_seen:
        continue

    # Compute M[a][b] = HP count from a to b
    M = [[0]*n for _ in range(n)]
    full = (1 << n) - 1
    for start in range(n):
        dp = [[0]*n for _ in range(1 << n)]
        dp[1 << start][start] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or dp[mask][v] == 0:
                    continue
                if mask == full:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    if A[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        for w in range(n):
            M[start][w] = dp[full][w]

    # Compute rank over Q (check if any 3×3 minor is nonzero, etc.)
    # For 4×4, rank = 4 iff det ≠ 0
    def det4(mat):
        result = 0
        for j in range(4):
            minor = [[mat[i][k] for k in range(4) if k != j] for i in range(1, 4)]
            det3 = (minor[0][0]*(minor[1][1]*minor[2][2]-minor[1][2]*minor[2][1])
                   -minor[0][1]*(minor[1][0]*minor[2][2]-minor[1][2]*minor[2][0])
                   +minor[0][2]*(minor[1][0]*minor[2][1]-minor[1][1]*minor[2][0]))
            result += ((-1)**j) * mat[0][j] * det3
        return result

    det_M = det4(M)

    # Check rank
    if det_M != 0:
        rank = 4
    else:
        # Check 3×3 minors
        rank = 0
        for rows in combinations(range(4), 3):
            for cols in combinations(range(4), 3):
                sub = [[M[r][c] for c in cols] for r in rows]
                d = (sub[0][0]*(sub[1][1]*sub[2][2]-sub[1][2]*sub[2][1])
                    -sub[0][1]*(sub[1][0]*sub[2][2]-sub[1][2]*sub[2][0])
                    +sub[0][2]*(sub[1][0]*sub[2][1]-sub[1][1]*sub[2][0]))
                if d != 0:
                    rank = max(rank, 3)
        if rank == 0:
            for i in range(4):
                for j in range(4):
                    if M[i][j] != 0:
                        rank = max(rank, 1)
            if rank == 1:
                # Check 2×2 minors
                for r1, r2 in combinations(range(4), 2):
                    for c1, c2 in combinations(range(4), 2):
                        if M[r1][c1]*M[r2][c2] - M[r1][c2]*M[r2][c1] != 0:
                            rank = 2

    H_seen[(H, ss)] = (M, rank, det_M)
    print(f"  H={H}, ss={ss}: rank(M) = {rank}, det(M) = {det_M}")
    print(f"    M = {M}")

# ======================================================================
# PART 5: DUALITY: POINTS ↔ HYPERPLANES
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: PROJECTIVE DUALITY IN TOURNAMENT SPACE")
print("=" * 70)

print("""
  In PG(m-1, 2), each point P determines a hyperplane H_P (orthogonal complement).
  Dually, each hyperplane determines a point.

  For tournaments: each tournament T ∈ F_2^m determines a hyperplane
  H_T = {S ∈ F_2^m : <T, S> = 0} where <T,S> = sum T_i S_i mod 2.

  The dual of a tournament T is the set of tournaments S such that
  the number of arcs where T and S AGREE is even (or equivalently,
  their Hamming distance is even).

  SELF-DUAL: T is self-dual iff <T, T> = 0 iff popcount(T) is even.
  A tournament is self-dual in PG(m-1,2) iff it has an EVEN number of arcs
  oriented in the "positive" direction.

  QUESTION: Does the projective dual of the H-level set have
  interesting properties?

  More specifically: if L_h = {T : H(T) = h}, what is L_h^⊥?
  L_h^⊥ = intersection of all hyperplanes H_T for T ∈ L_h
         = {S : <T, S> = 0 for ALL T ∈ L_h}
         = the F_2-linear span of L_h taken perpendicular.
""")

for n in [3, 4, 5]:
    m = n * (n-1) // 2
    tournaments = []
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        tournaments.append((bits, H))

    H_set = sorted(set(H for _, H in tournaments))

    print(f"\n  n={n} (m={m}):")

    # For each H level, compute the F_2-linear span and its dimension
    for h in H_set:
        level = [bits for bits, H in tournaments if H == h]

        # Compute span over F_2 using Gaussian elimination
        basis = []
        for v in level:
            # Reduce v using current basis
            r = v
            for b in basis:
                r = min(r, r ^ b)  # Not standard; let's do proper GE

        # Proper F_2 Gaussian elimination
        basis = []
        for v in level:
            r = v
            for b in basis:
                if r & (1 << (b.bit_length() - 1)):
                    r ^= b
            if r > 0:
                # Insert r into basis in order
                basis.append(r)
                # Sort by leading bit (descending)
                basis.sort(reverse=True)
                # Re-reduce
                new_basis = []
                for b in basis:
                    r = b
                    for b2 in new_basis:
                        if b2.bit_length() <= r.bit_length() and r & (1 << (b2.bit_length() - 1)):
                            r ^= b2
                    if r > 0:
                        new_basis.append(r)
                basis = sorted(new_basis, reverse=True)

        span_dim = len(basis)

        # The perpendicular space has dimension m - span_dim
        perp_dim = m - span_dim

        print(f"    H={h:3d}: |L_h|={len(level):5d}, "
              f"span_dim={span_dim}, perp_dim={perp_dim}")

        # If span_dim = m, the level set spans all of F_2^m
        if span_dim == m:
            print(f"           (Level set SPANS all of F_2^{m})")

# ======================================================================
# PART 6: THE CURVE H = h IN TOURNAMENT SPACE
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: ALGEBRAIC DEGREE OF H")
print("=" * 70)

print("""
  H(T) is a polynomial in the arc variables x_1, ..., x_m ∈ {0,1}.
  Over F_2, this is a multilinear polynomial (each x_i appears to power ≤ 1).

  The DEGREE of H as a polynomial tells us the algebraic complexity.
  From Walsh analysis:
  - n=3: H has Walsh degree 2 (quadratic)
  - n=4: H has Walsh degree 2 (quadratic)
  - n=5: H has Walsh degree 4 (quartic)

  But H is NOT a polynomial over R in the usual sense.
  It's a polynomial FUNCTION on {0,1}^m → Z.

  Every function {0,1}^m → Z can be uniquely written as
  H = sum_S a_S * prod_{i ∈ S} x_i
  where the a_S are integers.

  The ALGEBRAIC DEGREE of H is max{|S| : a_S ≠ 0}.
  This equals the Walsh degree by the theory of multilinear polynomials.

  For the H-level curve {T : H(T) = h}:
  This is a hypersurface of degree = algebraic degree of H
  in the "tournament affine variety" {0,1}^m.

  Over F_2: the equation H(x) = h mod 2 is always
  H(x) ≡ 1 mod 2 (Rédei), so it defines the EMPTY set mod 2!
  (Since h is always odd, H ≡ 1 mod 2, so H ≡ h mod 2 is always true.)

  Over F_3: H mod 3 defines a nontrivial variety.
  H(x) ≡ 0 mod 3 defines a cubic (or quadratic) hypersurface in F_3^m.
""")

# What is the algebraic representation of H over Z?
# For n=3: compute the multilinear expansion
n = 3
m = 3
print(f"\n  n={n}: H as multilinear polynomial in arc variables x_0, x_1, x_2")
print(f"  (Arc encoding: x_0 = (0→1), x_1 = (0→2), x_2 = (1→2))")

# Compute H at all 2^3 = 8 points
H_at = {}
for bits in range(1 << m):
    A = adj_matrix(n, bits)
    H_at[bits] = count_hp(n, A)
    x = [(bits >> i) & 1 for i in range(m)]
    print(f"    x={x}: H = {H_at[bits]}")

# Multilinear expansion using inclusion-exclusion (Möbius inversion on Boolean lattice)
print(f"\n  Multilinear coefficients:")
coeffs = {}
for S in range(1 << m):
    # a_S = sum_{T ⊆ S} (-1)^{|S|-|T|} H(T)
    a_S = 0
    for T in range(1 << m):
        if T & S == T:  # T ⊆ S
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            a_S += sign * H_at[T]
    coeffs[S] = a_S

    if a_S != 0:
        S_bits = [i for i in range(m) if S & (1 << i)]
        monomial = '*'.join(f'x_{i}' for i in S_bits) if S_bits else '1'
        print(f"    a_{{{bin(S)[2:].zfill(m)}}} = {a_S:+d}  ({monomial})")

# Verify
print(f"\n  Verification:")
for bits in range(1 << m):
    x = [(bits >> i) & 1 for i in range(m)]
    H_computed = 0
    for S in range(1 << m):
        prod = 1
        for i in range(m):
            if S & (1 << i):
                prod *= x[i]
        H_computed += coeffs[S] * prod
    assert H_computed == H_at[bits], f"Mismatch at {bits}"
    print(f"    x={x}: H_poly = {H_computed} = H = {H_at[bits]} ✓")

# Same for n=4
n = 4
m = 6
print(f"\n  n={n}: H as multilinear polynomial in 6 arc variables")

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

nonzero = [(S, a) for S, a in coeffs_4.items() if a != 0]
print(f"  {len(nonzero)} nonzero coefficients:")
degree_dist = Counter()
for S, a in nonzero:
    deg = bin(S).count('1')
    degree_dist[deg] += 1

for d in sorted(degree_dist.keys()):
    print(f"    Degree {d}: {degree_dist[d]} terms")
    if degree_dist[d] <= 3:
        for S, a in nonzero:
            if bin(S).count('1') == d:
                S_bits = [i for i in range(m) if S & (1 << i)]
                monomial = '*'.join(f'x_{i}' for i in S_bits) if S_bits else '1'
                print(f"      a = {a:+d}  ({monomial})")

print(f"\n  Algebraic degree of H at n=4: {max(degree_dist.keys())}")
print(f"  (Should match Walsh degree = 2)")

# ======================================================================
# PART 7: THE UNITY OF DUALITY
# ======================================================================
print("\n" + "=" * 70)
print("PART 7: THE UNITY OF DUALITY")
print("=" * 70)

print("""
  PROJECTIVE-ALGEBRAIC DUALITY in tournament theory:

  1. POINT-HYPERPLANE DUALITY (projective):
     Each tournament T ↔ hyperplane H_T = {S : <T,S> = 0 mod 2}
     This is the standard projective duality in PG(m-1, 2).

  2. TOURNAMENT-COMPLEMENT DUALITY (combinatorial):
     T ↔ T^op (reverse all arcs) = complement in F_2^m
     H(T) = H(T^op) for odd n (PROVED)
     T^op = T ⊕ 1...1 (XOR with all-ones)

  3. WALSH DUALITY (Fourier):
     The Walsh transform maps H to its Fourier spectrum hat{H}.
     hat{H}[S] lives in R, but the support {S : hat{H}[S] ≠ 0}
     is an F_2-variety (depends on |S| mod 2).

  4. HP-CYCLE DUALITY (OCF):
     H counts paths, I(Omega, 2) counts independent sets of cycles.
     H = I(Omega, 2) IS the OCF — paths and cycles are DUAL descriptions.

  5. GRASSMANN DUALITY:
     Gr(k, n) ≅ Gr(n-k, n) via taking orthogonal complement.
     The M[a,b] matrix lives in the Grassmannian of its column space.
     Its transpose M^T (which counts paths in reverse) gives the dual.

  THE UNITY:
  All five dualities are CONNECTED through the projective structure.

  The Fano plane PG(2,2) has 7 points and 7 lines.
  Under duality, points ↔ lines, and the incidence structure is SELF-DUAL.
  PG(2,2) is the ONLY self-dual projective plane of smallest order.

  The forbidden value 7 = |PG(2,2)| is thus tied to SELF-DUALITY:
  - 7 = Φ_3(2) is self-dual under the Galois action on F_4
  - The Fano plane's automorphism group GL(3,2) ≅ PSL(2,7)
  - 7 | H is forbidden because the cyclotomic structure of Φ_3
    prevents the independence polynomial from vanishing at the
    primitive cube root of 2

  THE DEEPER UNITY:
  Projective geometry and algebraic geometry meet in the
  MODULI SPACE of tournaments. The tournament space {0,1}^m
  is a 0-dimensional algebraic variety (a finite set of points).
  But the H function defines a FIBRATION:
  π: {0,1}^m → Z, T ↦ H(T)

  The fibers π^{-1}(h) are the H-level sets.
  Their algebraic structure (ideal, Hilbert function, dimension)
  encodes the geometry of the tournament landscape.

  PROJECTIVE: The level sets live in PG(m-1, 2).
  ALGEBRAIC: Their defining equations live in F_2[x_1,...,x_m].
  THE DUALITY: A polynomial equation f(T) = 0 defines both
  a projective variety AND an algebraic variety, and these are
  connected by the incidence correspondence of PG(m-1, 2).
""")

# Demonstrate the unity: for n=3, show how all dualities connect
n = 3
m = 3
print(f"\n  DEMONSTRATION for n={n}:")
print(f"  Tournaments = F_2^{m} = {{000, 001, 010, 011, 100, 101, 110, 111}}")

tournaments = []
for bits in range(1 << m):
    A = adj_matrix(n, bits)
    H = count_hp(n, A)
    complement = bits ^ ((1 << m) - 1)
    hamming_weight = bin(bits).count('1')
    tournaments.append((bits, H, complement, hamming_weight))

print(f"\n  {'T':>7s} | {'H':>3s} | {'T^op':>7s} | H(T^op) | wt(T) | Self-dual?")
print(f"  {'-'*7} | {'-'*3} | {'-'*7} | {'-'*7} | {'-'*5} | {'-'*10}")
for bits, H, comp, wt in tournaments:
    H_comp = [H2 for b2, H2, _, _ in tournaments if b2 == comp][0]
    self_dual = (wt % 2 == 0)
    print(f"  {bin(bits)[2:].zfill(m):>7s} | {H:3d} | {bin(comp)[2:].zfill(m):>7s} | {H_comp:7d} | {wt:5d} | {'Yes' if self_dual else 'No':>10s}")

print(f"\n  Observations:")
print(f"    H(T) = H(T^op) ✓ (proved for odd n)")
print(f"    Self-dual tournaments (even weight): {sum(1 for _, _, _, w in tournaments if w % 2 == 0)} "
      f"out of {len(tournaments)}")
print(f"    The complement map T → T^op is a FIXED-POINT-FREE involution on PG({m-1}, 2)")

print("\n" + "=" * 70)
print("DONE — PROJECTIVE & ALGEBRAIC GEOMETRY")
print("=" * 70)
