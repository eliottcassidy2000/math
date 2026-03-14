#!/usr/bin/env python3
"""
PROJECTIVE & ALGEBRAIC GEOMETRY — DEEPEST EXTENSION (PART 3)
opus-2026-03-14-S71n

The deepest connections: hyperdeterminants, flag varieties, derived categories,
and the categorical unification of all dualities.

Building on:
- Walsh basis = Segre basis (PROVED)
- H is multiplicative under direct sum (PROVED)
- Gram matrix eigenvalues {4,-2,-2} (n=3)
- Uniform gradient only for cyclic tournaments
- Six dualities identified

NOW exploring:
1. Cayley's hyperdeterminant and tournament critical points
2. Flag variety and tournament filtrations
3. The spectral curve of the tournament matrix
4. Resultant and discriminant of H-level polynomials
5. The moment map and symplectic geometry
6. Derived category perspective: tournaments as sheaves
7. Hodge structure on tournament cohomology
8. The period map: tournaments → period domain
9. Motivic integration over tournament space
10. The GRAND SYNTHESIS
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict
import math

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

def get_all_tournaments(n):
    m = n*(n-1)//2
    results = {}
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        H = count_hp(n, A)
        results[bits] = H
    return results

print("=" * 70)
print("PROJECTIVE & ALGEBRAIC GEOMETRY — DEEPEST EXTENSION (PART 3)")
print("opus-2026-03-14-S71n")
print("=" * 70)

# ======================================================================
# PART 1: CAYLEY'S HYPERDETERMINANT
# ======================================================================
print(f"\n{'='*70}")
print("PART 1: CAYLEY'S HYPERDETERMINANT AND TOURNAMENT SINGULARITIES")
print(f"{'='*70}")

print("""
  For a 2x2x2 tensor T = (t_ijk), Cayley's hyperdeterminant is:
  Det(T) = t000^2 t111^2 + t001^2 t110^2 + t010^2 t101^2 + t011^2 t100^2
           - 2(t000 t001 t110 t111 + t000 t010 t101 t111
              + t000 t011 t100 t111 + t001 t010 t101 t110
              + t001 t011 t100 t110 + t010 t011 t100 t101)
           + 4(t000 t011 t101 t110 + t001 t010 t100 t111)

  For tournament H viewed as a tensor on {0,1}^m:
  At n=3 (m=3): H is a 2x2x2 tensor!
  H[x0,x1,x2] = number of Hamiltonian paths in tournament T(x0,x1,x2).

  The SINGULARITIES of the Segre variety under H are detected
  by the hyperdeterminant.
""")

# Compute H as a 2x2x2 tensor at n=3
all_T3 = get_all_tournaments(3)
tensor = {}
for bits in range(8):
    x = [(bits >> i) & 1 for i in range(3)]
    tensor[tuple(x)] = all_T3[bits]

print("  H tensor at n=3:")
for x0 in range(2):
    for x1 in range(2):
        for x2 in range(2):
            print(f"    H[{x0},{x1},{x2}] = {tensor[(x0,x1,x2)]}")

# Compute hyperdeterminant
t = tensor
Det = (t[(0,0,0)]**2 * t[(1,1,1)]**2 + t[(0,0,1)]**2 * t[(1,1,0)]**2
     + t[(0,1,0)]**2 * t[(1,0,1)]**2 + t[(0,1,1)]**2 * t[(1,0,0)]**2
     - 2*(t[(0,0,0)]*t[(0,0,1)]*t[(1,1,0)]*t[(1,1,1)]
        + t[(0,0,0)]*t[(0,1,0)]*t[(1,0,1)]*t[(1,1,1)]
        + t[(0,0,0)]*t[(0,1,1)]*t[(1,0,0)]*t[(1,1,1)]
        + t[(0,0,1)]*t[(0,1,0)]*t[(1,0,1)]*t[(1,1,0)]
        + t[(0,0,1)]*t[(0,1,1)]*t[(1,0,0)]*t[(1,1,0)]
        + t[(0,1,0)]*t[(0,1,1)]*t[(1,0,0)]*t[(1,0,1)])
     + 4*(t[(0,0,0)]*t[(0,1,1)]*t[(1,0,1)]*t[(1,1,0)]
        + t[(0,0,1)]*t[(0,1,0)]*t[(1,0,0)]*t[(1,1,1)]))

print(f"\n  Cayley hyperdeterminant Det(H) = {Det}")
print(f"  Det(H) = 0 means H has a SINGULAR point on the Segre variety.")
print(f"  Since Det(H) = {Det}, H is {'SINGULAR' if Det == 0 else 'NON-SINGULAR'} on (P^1)^3.")

# ======================================================================
# PART 2: THE SPECTRAL CURVE OF THE TOURNAMENT MATRIX
# ======================================================================
print(f"\n{'='*70}")
print("PART 2: SPECTRAL CURVE OF THE TOURNAMENT MATRIX")
print(f"{'='*70}")

print("""
  For a tournament T with adjacency matrix A, the SPECTRAL CURVE is
  C_T: det(xI - A) = 0  in  P^1 x P^{n-1}.

  The eigenvalues of A lie on C_T. For tournaments:
  - tr(A) = 0 (diagonal is 0)
  - A + A^T = J - I (every off-diagonal entry sums to 1)
  - So A = (J-I)/2 + B/2 where B is skew-symmetric
  - Eigenvalues of A = (n-1)/2 + eigenvalues of B/2

  The SKEW-SYMMETRIC MATRIX B = A - A^T = 2A - J + I has:
  - B[i,j] = +1 if i->j, -1 if j->i, 0 if i=j
  - B is the SIGNED ADJACENCY MATRIX

  Eigenvalues of B: purely imaginary (B is real skew-symmetric).
  B = i*H_B where H_B is Hermitian.
  Eigenvalues of B: {0} union {+/- i*lambda_k} for some lambda_k > 0.
  (0 always appears when n is odd.)
""")

# Compute eigenvalues of B for small tournaments
import cmath

def char_poly_2x2(M):
    """Characteristic polynomial of 2x2 matrix: x^2 - tr*x + det"""
    tr = M[0][0] + M[1][1]
    det = M[0][0]*M[1][1] - M[0][1]*M[1][0]
    return tr, det

def eigenvalues_3x3_skew(B):
    """Eigenvalues of 3x3 skew-symmetric matrix: 0, ±i*sqrt(B01^2+B02^2+B12^2)"""
    val = B[0][1]**2 + B[0][2]**2 + B[1][2]**2
    return 0, val**0.5

for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    spectral_classes = Counter()

    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        B = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    B[i][j] = 2*A[i][j] - 1

        if n == 3:
            _, val = eigenvalues_3x3_skew(B)
            spectral_classes[round(val, 6)] += 1
        else:
            # Compute B^T B eigenvalues (they're lambda_k^2)
            BtB = [[sum(B[k][i]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
            # Trace and other invariants
            tr_BtB = sum(BtB[i][i] for i in range(n))
            spectral_classes[tr_BtB] += 1

    print(f"\n  n={n}: Spectral invariant distribution:")
    if n == 3:
        print(f"  (eigenvalue magnitude |lambda|):")
    else:
        print(f"  (tr(B^T B) = sum of squared eigenvalue magnitudes):")
    for val, count in sorted(spectral_classes.items()):
        H_vals = set()
        for bits in range(1 << m):
            A = adj_matrix(n, bits)
            B = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    if i != j:
                        B[i][j] = 2*A[i][j] - 1
            if n == 3:
                _, v = eigenvalues_3x3_skew(B)
                if round(v, 6) == val:
                    H_vals.add(all_T[bits])
            else:
                BtB = [[sum(B[k][i]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
                if sum(BtB[i][i] for i in range(n)) == val:
                    H_vals.add(all_T[bits])
        print(f"    {val}: {count} tournaments, H values = {sorted(H_vals)}")

# ======================================================================
# PART 3: FLAG VARIETY AND TOURNAMENT FILTRATIONS
# ======================================================================
print(f"\n{'='*70}")
print("PART 3: FLAG VARIETY AND TOURNAMENT FILTRATIONS")
print(f"{'='*70}")

print("""
  A tournament defines a TOTAL ORDER on vertices (its transitive closure),
  but only if the tournament is transitive.

  For general tournaments: the tournament defines a PARTIAL ORDER
  via reachability: v >= w if there's a directed path v -> ... -> w.

  But EVERY tournament is reachable from every other vertex!
  (Strong connectivity for n >= 2.)

  So the partial order is trivial. Instead, consider:

  The SCORE SEQUENCE (s_1 <= s_2 <= ... <= s_n) defines a FLAG:
  F_1 ⊂ F_2 ⊂ ... ⊂ F_n where F_k = {vertices with score <= s_k}.

  But vertices with the same score are not distinguishable in this flag.
  The flag type is determined by the PARTITION of n into groups of
  equal scores.

  The FLAG VARIETY Fl(d_1,...,d_k; n) parametrizes chains of subspaces.

  QUESTION: Does H factor through the flag variety?
  Answer: NO — H depends on more than the score sequence.
  But H restricted to the flag stratum is simpler.
""")

# Group tournaments by score sequence and check H variance
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    score_to_H = defaultdict(list)
    for bits, H in all_T.items():
        A = adj_matrix(n, bits)
        ss = score_seq(n, A)
        score_to_H[ss].append(H)

    print(f"\n  n={n}: H variance within score sequence strata:")
    for ss in sorted(score_to_H.keys()):
        H_vals = score_to_H[ss]
        H_set = sorted(set(H_vals))
        if len(H_set) > 1:
            mean_H = sum(H_vals) / len(H_vals)
            print(f"    {ss}: |stratum|={len(H_vals)}, H={H_set}, mean={mean_H:.1f}")
        else:
            print(f"    {ss}: |stratum|={len(H_vals)}, H={H_set[0]} (constant!)")

# ======================================================================
# PART 4: RESULTANT OF H(x) - h₁ AND H(x) - h₂
# ======================================================================
print(f"\n{'='*70}")
print("PART 4: RESULTANT AND DISCRIMINANT STRUCTURE")
print(f"{'='*70}")

print("""
  The DISCRIMINANT of H is the set of values h where the
  fiber H^{-1}(h) is singular (has a multiple root).

  But over F_2, EVERY fiber is "singular" in the algebraic sense.

  More interesting: the TROPICAL DISCRIMINANT.
  This is the set of h where the level set {H=h} has
  "extra" structure (e.g., unusual F_2-rank, symmetry, etc.).

  Even more interesting: the NEWTON POLYTOPE of H.
  As a multilinear polynomial, H's Newton polytope is a
  subpolytope of the unit cube [0,1]^m.

  The support of H = set of monomials with nonzero coefficient.
""")

# Compute Newton polytope for small n
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    # Multilinear coefficients via Mobius inversion
    ml_coeffs = {}
    for S in range(1 << m):
        val = 0
        for T in range(1 << m):
            if T & S != T:
                continue
            sign = (-1) ** (bin(S).count('1') - bin(T).count('1'))
            val += sign * all_T[T]
        if val != 0:
            ml_coeffs[S] = val

    support_sizes = Counter()
    for S in ml_coeffs:
        support_sizes[bin(S).count('1')] += 1

    print(f"\n  n={n} (m={m}): Newton polytope support")
    print(f"    Total nonzero monomials: {len(ml_coeffs)}")
    print(f"    By degree: {dict(sorted(support_sizes.items()))}")

    # The degree of H as a multilinear polynomial
    max_deg = max(bin(S).count('1') for S in ml_coeffs)
    print(f"    Max degree: {max_deg}")

    # Distribution of coefficient magnitudes
    coeff_dist = Counter(abs(v) for v in ml_coeffs.values())
    print(f"    Coefficient magnitudes: {dict(sorted(coeff_dist.items()))}")

    # Are coefficients always ±2^k?
    all_powers_of_2 = all(abs(v) & (abs(v)-1) == 0 for v in ml_coeffs.values())
    print(f"    All coefficients are +-2^k: {all_powers_of_2}")

# ======================================================================
# PART 5: THE MOMENT MAP — SYMPLECTIC GEOMETRY
# ======================================================================
print(f"\n{'='*70}")
print("PART 5: THE MOMENT MAP — SYMPLECTIC GEOMETRY")
print(f"{'='*70}")

print("""
  The torus T = (S^1)^m acts on (P^1)^m by rotation.
  The MOMENT MAP mu: (P^1)^m -> R^m sends a point to its
  "angular momentum" in each factor.

  For (P^1)^m: mu(x) = (|x_1|^2/(|x_0|^2+|x_1|^2), ..., |x_m|^2/...).
  On the Boolean lattice {0,1}^m: mu(x) = x itself!

  So the tournament hypercube {0,1}^m IS the image of the moment map.
  The SYMPLECTIC REDUCTION at level h gives the symplectic quotient.

  For tournaments: mu^{-1}(x) = {x} (the moment map is injective on F_2).
  But the H-fibers mu^{-1}(H=h) are the level sets of H restricted
  to the moment polytope.

  The DUISTERMAAT-HECKMAN theorem says:
  The push-forward of the symplectic volume by H is piecewise polynomial.
  For the tournament hypercube: this is the H-distribution!

  Let's check: is the H-distribution piecewise polynomial in h?
""")

# Check H-distribution for polynomial structure
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)
    H_counter = Counter(all_T.values())

    H_vals = sorted(H_counter.keys())
    print(f"\n  n={n}: H-distribution as function of h:")
    for h in H_vals:
        count = H_counter[h]
        # Normalize to fraction of total
        frac = count / (1 << m)
        print(f"    h={h:3d}: count={count:5d}, frac={frac:.6f}")

    # Check if counts satisfy any polynomial relation
    # At n=3: counts are [6, 2] at h=[1, 3]
    # At n=4: counts are [24, 16, 24] at h=[1, 3, 5]
    # At n=5: counts are [120, 120, 240, 240, 120, 120, 64] at h=[1,3,5,9,11,13,15]

    # Check symmetry: count(h) = count(n! - h)?
    # At n=5: max H = 15 = 5!/8, not 5!
    # What about complement: count(h) = count(H_max + H_min - h)?
    h_max = max(H_vals)
    h_min = min(H_vals)
    h_sum = h_max + h_min
    sym_pairs = [(h, h_sum - h) for h in H_vals if h <= h_sum/2]
    is_symmetric = all(H_counter.get(h, 0) == H_counter.get(h_sum - h, 0) for h, _ in sym_pairs)
    print(f"    Symmetric about (h_min+h_max)/2 = {h_sum/2}: {is_symmetric}")

    # Check: does h_min + h_max work with complement?
    print(f"    h_max + h_min = {h_sum}")

    # Actually H(T) + H(T^op) is not constant...
    # But the level sets |{T: H(T)=h}| might be related to |{T: H(T^op)=h}|
    # Since H(T^op) = H(T) for odd n (complement preserves H for odd n)...
    # Wait, that's not true. H(T) = H(T^op) only for self-complementary T.

    # Check H(T) vs H(T^op)
    complement_mask = (1 << m) - 1
    same_count = sum(1 for bits in range(1 << m) if all_T[bits] == all_T[bits ^ complement_mask])
    print(f"    H(T) = H(T^op): {same_count}/{1<<m}")

# ======================================================================
# PART 6: THE HESSIAN AND MORSE THEORY
# ======================================================================
print(f"\n{'='*70}")
print("PART 6: HESSIAN AND MORSE THEORY ON THE TOURNAMENT HYPERCUBE")
print(f"{'='*70}")

print("""
  Morse theory on a discrete space: the FORMAN DISCRETE MORSE THEORY.

  On the tournament hypercube Q_m with height function H:
  A vertex is CRITICAL if it's a local extremum.
  The MORSE INEQUALITY says:
  #{critical points of index k} >= b_k(Q_m)

  The Betti numbers of Q_m = {0,1}^m (as a simplicial complex of the cube):
  b_k = C(m, k) (the hypercube has the topology of the m-torus!)

  Wait — Q_m as a GRAPH has b_0 = 1 (connected), b_1 = C(m,2) - m + 1.
  But as a CUBE COMPLEX: Q_m ~ (S^1)^m, so b_k = C(m,k).

  For Morse theory on the graph Q_m with height H:
  The HESSIAN of H at a critical point is the matrix of
  second differences: D_{ij} = H(x + e_i + e_j) - H(x + e_i) - H(x + e_j) + H(x).

  The INDEX of a critical point = number of negative eigenvalues of D.
""")

# Compute Hessian and Morse index at critical points
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    # Find critical points (local maxima and minima)
    critical = []
    for bits in range(1 << m):
        H = all_T[bits]
        is_max = True
        is_min = True
        for arc in range(m):
            nbr = bits ^ (1 << arc)
            if all_T[nbr] > H:
                is_max = False
            if all_T[nbr] < H:
                is_min = False
        if is_max or is_min:
            # Compute Hessian
            D = [[0]*m for _ in range(m)]
            for i in range(m):
                for j in range(m):
                    # D[i,j] = H(x+ei+ej) - H(x+ei) - H(x+ej) + H(x)
                    xi = bits ^ (1 << i)
                    xj = bits ^ (1 << j)
                    xij = bits ^ (1 << i) ^ (1 << j)
                    D[i][j] = all_T[xij] - all_T[xi] - all_T[xj] + H

                # Diagonal: D[i,i] = H(x) - H(x+ei) (wait, x+ei+ei = x for binary)
                # So D[i,i] = H(x) - 2*H(x+ei) + H(x) = 2H - 2H(x+ei)
                # Actually D[i,i] is not well-defined for binary variables
                # Since e_i + e_i = 0 in F_2
                D[i][i] = 0  # Set diagonal to 0

            # Count positive/negative off-diagonal elements
            neg_count = sum(1 for i in range(m) for j in range(i+1,m) if D[i][j] < 0)
            pos_count = sum(1 for i in range(m) for j in range(i+1,m) if D[i][j] > 0)
            zero_count = sum(1 for i in range(m) for j in range(i+1,m) if D[i][j] == 0)

            critical.append((H, "MAX" if is_max else "MIN", neg_count, pos_count, zero_count))

    # Summarize
    max_critical = [c for c in critical if c[1] == "MAX"]
    min_critical = [c for c in critical if c[1] == "MIN"]

    print(f"\n  n={n}: Critical points of H on Q_{m}:")
    print(f"    Local maxima: {len(max_critical)}")
    for H, _, neg, pos, zero in sorted(set(max_critical)):
        count = max_critical.count((H, "MAX", neg, pos, zero))
        print(f"      H={H}: {count} points, Hessian off-diag: {neg} neg, {pos} pos, {zero} zero")
    print(f"    Local minima: {len(min_critical)}")
    for H, _, neg, pos, zero in sorted(set(min_critical)):
        count = min_critical.count((H, "MIN", neg, pos, zero))
        print(f"      H={H}: {count} points, Hessian off-diag: {neg} neg, {pos} pos, {zero} zero")

# ======================================================================
# PART 7: HODGE STRUCTURE ON TOURNAMENT COHOMOLOGY
# ======================================================================
print(f"\n{'='*70}")
print("PART 7: HODGE DECOMPOSITION OF THE H FUNCTION SPACE")
print(f"{'='*70}")

print("""
  The space of functions {0,1}^m -> R has dimension 2^m.
  The Walsh decomposition gives a grading by degree:
  V = V_0 + V_1 + ... + V_m where V_k has dimension C(m,k).

  This is EXACTLY the Hodge decomposition H^k(T^m) = C^{C(m,k)}
  of the cohomology of the m-torus T^m = (S^1)^m.

  The HODGE NUMBERS h^{p,q} of (P^1)^m are:
  h^{p,q} = C(m,p) * C(m,q) / C(m,p+q) * something...

  Actually, for (P^1)^m: H^{p,q} = C(m,p) * delta_{p,q}.
  So the Hodge diamond is concentrated on the diagonal.

  For tournaments: the Walsh decomposition places H in specific
  Hodge components. Since H has Walsh components only at EVEN degrees
  (for odd n), H lives in H^{0,0} + H^{2,2} + H^{4,4} + ...

  The F-polynomial further decomposes: F_k(T) lives in
  the degree-k part of some auxiliary Hodge structure.
""")

# Walsh degree decomposition of H
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    # Walsh transform
    walsh_hat = {}
    for S in range(1 << m):
        val = 0
        for bits in range(1 << m):
            sign = 1
            for i in range(m):
                if (S >> i) & 1 and (bits >> i) & 1:
                    sign *= -1
            val += sign * all_T[bits]
        walsh_hat[S] = val

    # Group by Walsh degree
    degree_energy = defaultdict(int)
    degree_count = defaultdict(int)
    for S, val in walsh_hat.items():
        deg = bin(S).count('1')
        degree_energy[deg] += val**2
        if val != 0:
            degree_count[deg] += 1

    total_energy = sum(degree_energy.values())
    print(f"\n  n={n}: Walsh-Hodge decomposition of H:")
    print(f"    Degree  C(m,deg)  NonzeroCoeffs  Energy      FractionOfTotal")
    for deg in sorted(degree_energy.keys()):
        if degree_energy[deg] > 0:
            frac = degree_energy[deg] / total_energy
            print(f"    {deg:5d}  {math.comb(m,deg):7d}  {degree_count[deg]:13d}  {degree_energy[deg]:10d}  {frac:.6f}")

# ======================================================================
# PART 8: THE TOURNAMENT ZETA FUNCTION
# ======================================================================
print(f"\n{'='*70}")
print("PART 8: THE TOURNAMENT ZETA FUNCTION")
print(f"{'='*70}")

print("""
  Define the TOURNAMENT ZETA FUNCTION:
  Z_n(s) = sum_T H(T)^{-s} = sum_h mult(h) * h^{-s}

  This is an analogue of the RIEMANN ZETA FUNCTION where
  "primes" are replaced by tournament H-values.

  Properties:
  - Z_n(0) = 2^m (total number of tournaments)
  - Z_n(1) = sum_T 1/H(T) (harmonic sum)
  - The ABSCISSA OF CONVERGENCE is sigma_c = 0
    (since there are finitely many terms)

  More interesting: the EULER PRODUCT.
  If H is multiplicative under direct sum, then:
  Z_n(s) = prod_{indecomposable T} (1 - H(T)^{-s})^{-1} ???

  No — that's for a free abelian monoid of indecomposables.
  The tournament direct sum is not free (different decompositions possible).

  But the DIRICHLET SERIES is still well-defined:
  sum_{h >= 1} a_h * h^{-s} where a_h = |{T : H(T) = h}|.
""")

for n in [3, 4, 5, 6]:
    m = n*(n-1)//2
    if n <= 5:
        all_T = get_all_tournaments(n)
    else:
        all_T = None

    if all_T:
        H_counter = Counter(all_T.values())

        # Compute Z_n(s) for s = 1, 2
        Z1 = sum(count / h for h, count in H_counter.items())
        Z2 = sum(count / h**2 for h, count in H_counter.items())

        print(f"\n  n={n}: Tournament zeta function")
        print(f"    Z_n(0) = {1 << m} (total tournaments)")
        print(f"    Z_n(1) = {Z1:.6f} (harmonic sum)")
        print(f"    Z_n(2) = {Z2:.6f}")
        print(f"    Z_n(1) / 2^m = {Z1 / (1<<m):.6f} (normalized)")

        # The "mean of 1/H"
        mean_inv_H = Z1 / (1 << m)
        mean_H = sum(all_T.values()) / (1 << m)
        print(f"    E[1/H] = {mean_inv_H:.6f}")
        print(f"    E[H] = {mean_H:.6f}")
        print(f"    E[H] * E[1/H] = {mean_H * mean_inv_H:.6f}")
        print(f"    (= 1 if H and 1/H were 'independent', >1 by Jensen's)")

# ======================================================================
# PART 9: THE PERIOD MAP
# ======================================================================
print(f"\n{'='*70}")
print("PART 9: THE PERIOD MAP — TOURNAMENTS AS HODGE STRUCTURES")
print(f"{'='*70}")

print("""
  In algebraic geometry, the PERIOD MAP associates to a variety
  its HODGE STRUCTURE (the decomposition of cohomology into
  H^{p,q} pieces).

  For tournaments: the "period" of a tournament T is the vector
  of Walsh coefficients hat{H}[S] for all S.

  This embeds the tournament in the PERIOD DOMAIN:
  Omega = {sequences (a_S)_S : sum_S a_S W_S(x) >= 1 for all x}

  The PERIOD MAP phi: T_n -> Omega sends T to its Walsh spectrum.

  TORELLI THEOREM ANALOGUE: Is phi injective?
  (Does the Walsh spectrum uniquely determine the tournament?)

  Answer: YES! The Walsh spectrum determines H completely,
  and H determines the tournament (? No — different T can have same H).

  But the FULL Walsh spectrum of ALL functions (H, t3, t5, ...)
  DOES determine T. The tournament IS its Walsh spectrum.

  Actually: the Walsh spectrum of H alone determines
  T only if H is injective on tournaments (rarely true).

  The Walsh spectrum of the INDICATOR FUNCTION 1_T determines T
  (trivially — it IS T in Walsh basis).
""")

# Check: how many tournaments have the same Walsh-H spectrum?
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    # Walsh transform of H (same as before)
    # Actually H is one function; the Walsh coefficients are determined
    # by the single number H at each point.
    # The "Walsh spectrum of H" is the set of hat{H}[S] values.
    # Different n give different H functions.

    # Better question: for a FIXED tournament T, what is the Walsh
    # expansion of delta_T (the indicator function)?
    # delta_T(x) = product_{i: T has bit i set} x_i * product_{i: T has bit i unset} (1-x_i)
    # In Walsh: delta_T = (1/2^m) * sum_S (-1)^{<S,T>} W_S
    # So the Walsh spectrum of delta_T is just (-1)^{<S,T>} / 2^m.

    # This means: the Walsh spectrum of delta_T uniquely identifies T.
    # But this is trivial.

    # More interesting: two tournaments T1, T2 have the SAME H-spectrum
    # iff H(T1) = H(T2). How does the Walsh CONTENT of H relate?

    # The Walsh content at each INDIVIDUAL tournament:
    # hat{H}[S] is determined by ALL tournaments, not by one.
    # So we can't define a "Walsh spectrum of T" from H alone.

    print(f"  n={n}: H-spectrum determines {len(set(all_T.values()))} out of {1<<m} tournaments")
    print(f"    (only {len(set(all_T.values()))} distinct H values)")

# ======================================================================
# PART 10: THE GRAND SYNTHESIS
# ======================================================================
print(f"\n{'='*70}")
print("PART 10: THE GRAND SYNTHESIS — ALL STRUCTURES UNIFIED")
print(f"{'='*70}")

print("""
  THE GRAND PICTURE:

  A tournament on n vertices is simultaneously:

  1. A POINT on the toric variety (P^1)^m
     (projective geometry — Segre embedding)

  2. A VECTOR in F_2^m
     (algebraic geometry — affine coordinates)

  3. A TENSOR in (C^2)^{otimes m}
     (multilinear algebra — hyperdeterminant)

  4. A FUNCTION on the Walsh spectrum
     (harmonic analysis — Fourier dual)

  5. A PARTITION FUNCTION for a statistical system
     (physics — multiplicativity under direct sum)

  6. An ELEMENT of the tournament monoid
     (algebra — direct sum operation)

  7. A CRITICAL POINT of the Segre variety under H
     (Morse theory — Hessian and index)

  8. A POINT on the period domain
     (Hodge theory — Walsh-Hodge decomposition)

  AND the H function is simultaneously:

  1. A SECTION of O(1,...,1) on (P^1)^m
     (line bundle — projective)

  2. A MULTILINEAR POLYNOMIAL on F_2^m
     (algebraic — polynomial ring)

  3. A LINEAR FUNCTIONAL on P^{2^m-1}
     (Segre-linearized — projective)

  4. A HERMITIAN OPERATOR on L^2(F_2^m)
     (spectral — eigenvalues are Walsh coefficients)

  5. A HEIGHT FUNCTION on the moment polytope
     (symplectic — Duistermaat-Heckman)

  6. A PFAFFIAN/HYPERDETERMINANT of a tensor
     (tensor algebra — Cayley)

  ALL THESE VIEWPOINTS ARE CONNECTED BY THE SEGRE EMBEDDING.
  The Segre embedding is the ROSETTA STONE of tournament theory.

  The SIX DUALITIES are unified:

  1. Projective (T <-> H-hyperplane): via Segre linearization
  2. Complement (T <-> T^op): via F_2 involution x -> 1-x
  3. Walsh-Fourier (tournament <-> spectrum): Segre coords = Walsh chars
  4. OCF (paths <-> cycles): H = I(Omega, 2) [Grinberg-Stanley]
  5. Grassmann (Gr(k,n) <-> Gr(n-k,n)): Plucker on tournament arcs
  6. Segre-Veronese (product <-> polynomial): the bridge itself

  And now a SEVENTH DUALITY:
  7. Statistical mechanics (tournament <-> partition function):
     H(T1 + T2) = H(T1) * H(T2) makes H a multiplicative character
     on the tournament monoid. Its LOG is an additive character —
     i.e., a HOMOMORPHISM to (R, +).
     The dual space of characters IS the "spectrum" of the monoid.

     This connects tournaments to NUMBER THEORY:
     - The tournament monoid under direct sum ~ the positive integers under multiplication
     - H is like a multiplicative arithmetic function
     - The tournament zeta function Z_n(s) is the Dirichlet series
     - "Prime tournaments" = indecomposable tournaments (not direct sums)
     - The PRIME FACTORIZATION of a tournament is its decomposition
       into indecomposable sub-tournaments
""")

# Count indecomposable tournaments at each n
print("  Indecomposable tournament count:")
for n in [3, 4, 5]:
    m = n*(n-1)//2
    all_T = get_all_tournaments(n)

    # A tournament is decomposable if there exists a subset S with
    # |S| in {1,...,n-1} such that all arcs go from S to V\S.
    # (This means T = T[S] + T[V\S])

    decomposable = 0
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        is_decomposable = False
        for size in range(1, n):
            for subset in combinations(range(n), size):
                S = set(subset)
                complement = set(range(n)) - S
                all_from_S_to_C = all(A[i][j] == 1 for i in S for j in complement)
                all_from_C_to_S = all(A[j][i] == 1 for i in S for j in complement)
                if all_from_S_to_C or all_from_C_to_S:
                    is_decomposable = True
                    break
            if is_decomposable:
                break
        if is_decomposable:
            decomposable += 1

    indecomposable = (1 << m) - decomposable
    print(f"  n={n}: {indecomposable}/{1<<m} indecomposable ({100*indecomposable/(1<<m):.1f}%)")

    # H values of decomposable vs indecomposable
    decomp_H = Counter()
    indecomp_H = Counter()
    for bits in range(1 << m):
        A = adj_matrix(n, bits)
        is_decomposable = False
        for size in range(1, n):
            for subset in combinations(range(n), size):
                S = set(subset)
                complement = set(range(n)) - S
                if all(A[i][j] == 1 for i in S for j in complement) or \
                   all(A[j][i] == 1 for i in S for j in complement):
                    is_decomposable = True
                    break
            if is_decomposable:
                break
        if is_decomposable:
            decomp_H[all_T[bits]] += 1
        else:
            indecomp_H[all_T[bits]] += 1

    print(f"    Decomposable H: {dict(sorted(decomp_H.items()))}")
    print(f"    Indecomposable H: {dict(sorted(indecomp_H.items()))}")

# Final: the connection to number theory
print(f"\n  ARITHMETIC TOURNAMENT ANALOGY:")
print(f"  ============================================")
print(f"  Integers           | Tournaments")
print(f"  ============================================")
print(f"  n ∈ Z+             | T tournament")
print(f"  n*m (mult)         | T1 + T2 (direct sum)")
print("  sigma_0(n)=#div    | H(T) = #HP")
print(f"  p prime            | T indecomposable")
print(f"  p^a (prime power)  | iterared sums of same type")
print("  zeta(s)=sum 1/n^s  | Z_n(s)=sum H(T)^{-s}")
print(f"  Euler product      | Product over indecomposables")
print(f"  Mobius function     | Inclusion-exclusion on sums")
print(f"  ============================================")

print("\n" + "=" * 70)
print("DONE — PROJECTIVE/ALGEBRAIC GEOMETRY DEEPEST EXTENSION (PART 3)")
print("=" * 70)
