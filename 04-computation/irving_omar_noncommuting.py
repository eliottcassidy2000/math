#!/usr/bin/env python3
"""
Irving-Omar noncommuting walk generating function and its connection
to the transfer matrix M[a,b].

Irving-Omar (arXiv:2412.10572):
  W_D(z) = det(I + z*X*A_bar) / det(I - z*X*A)
where X = diag(x_1, ..., x_n) with NONCOMMUTING x_i.

The Hamiltonian path count is:
  ham(D) = L_n(W_D(1))
where L_n extracts the coefficient of x_{i_1}*...*x_{i_n} (each x_i once).

KEY QUESTION: Can we express M[a,b] via L_n extraction?
  M[a,b] = L_{n-2}( something involving W restricted to paths a->...->b )?

For tournaments: A_bar = A^T, so:
  W_D(z) = det(I + z*X*A^T) / det(I - z*X*A)

The even r-powers of M[a,b] (THM-030) should correspond to some
symmetry property of W_D under z -> -z (or X -> -X).

Let's explore this connection computationally using SYMBOLIC noncommuting
variables (represented as strings/tuples tracking variable order).
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

# ===================================================================
# APPROACH 1: Direct path decomposition via inclusion-exclusion
# ===================================================================
# Instead of implementing full noncommuting algebra, observe:
#
# The Irving-Omar det/per formula (Prop 2) states:
#   ham(D) = sum_S det(A_bar[S,S]) * per(A[S^c, S^c])
#
# For tournaments: A_bar = A^T, so det(A_bar[S,S]) = det(A^T[S,S]) = det(A[S,S])
# Thus: ham(D) = sum_S det(A[S,S]) * per(A[S^c, S^c])
#
# Our transfer matrix:
#   M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
#
# E_a(S+a) = paths through S ending at a
# B_b(R+b) = paths through R starting at b
#
# The Irving-Omar decomposition is over subsets S of [n],
# while our M decomposition is over subsets S of [n]\{a,b}.
# And IO uses det*per (cycle covers), while we use E*B (paths).
#
# Can we BRIDGE these? Yes, via the Matrix Tree / Permanent connection.

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start:
            continue
        if end is not None and p[-1] != end:
            continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def ham_path_count(A):
    n = len(A)
    return count_paths_subset(A, list(range(n)))

def transfer_matrix(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_subset(A, S_verts, end=a)
                    bb = count_paths_subset(A, R_verts, start=b)
                    total += ((-1)**k) * ea * bb

            M[a][b] = total
    return M

def permanent(M):
    n = len(M)
    if n == 0:
        return 1
    total = 0
    for perm in permutations(range(n)):
        prod = 1
        for i in range(n):
            prod *= M[i][perm[i]]
        total += prod
    return total

def det_int(M):
    n = len(M)
    if n == 0:
        return 1
    return round(np.linalg.det(np.array(M, dtype=float)))

# ===================================================================
# APPROACH 2: The "endpoint-refined" Irving-Omar decomposition
# ===================================================================
#
# Irving-Omar's walk generating function with noncommuting variables
# can be expanded as:
#   W_D(z) = sum_k z^k * sum_{walks w of length k} x_{w_0}*...*x_{w_k}
#
# The L_n extraction picks out walks that visit EVERY vertex exactly once,
# i.e., Hamiltonian PATHS.
#
# If we condition on the endpoints a (start) and b (end), we get:
#   ham(D) = sum_{a,b} L_n(walks from a to b)
#
# This "endpoint-refined" count is exactly:
#   H(a,b) = number of Hamiltonian paths from a to b
#
# Our transfer matrix M[a,b] is NOT the same as H(a,b)!
# M[a,b] involves inclusion-exclusion over subset splits.
# But they're related: H = sum_{a,b} H(a,b) = tr(M) (odd n).
#
# Let's compute H(a,b) and compare with M[a,b].

def endpoint_ham_paths(A):
    """Compute H(a,b) = #Ham paths from a to b."""
    n = len(A)
    H_ab = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            H_ab[a][b] = count_paths_subset(A, list(range(n)), start=a, end=b)
    return H_ab

print("=" * 70)
print("M[a,b] vs H(a,b): TRANSFER MATRIX vs ENDPOINT PATHS")
print("=" * 70)

# Test tournaments
A_cyc = [[0,1,0],[0,0,1],[1,0,0]]  # n=3 cyclic
A_trans = [[0,1,1],[0,0,1],[0,0,0]]  # n=3 transitive
A_4a = [[0,1,1,0],[0,0,1,1],[0,0,0,1],[1,0,0,0]]  # n=4

for name, A in [("n=3 cyclic", A_cyc), ("n=3 trans", A_trans), ("n=4", A_4a)]:
    n = len(A)
    M = transfer_matrix(A)
    H_ab = endpoint_ham_paths(A)
    H = ham_path_count(A)

    print(f"\n--- {name} (H={H}) ---")
    print(f"  M[a,b] = {M.tolist()}")
    print(f"  H(a,b) = {H_ab.tolist()}")
    print(f"  M = H_ab? {np.allclose(M, H_ab)}")
    print(f"  tr(M) = {np.trace(M)}, sum(H_ab) = {H_ab.sum()}")
    print(f"  H_ab symmetric? {np.allclose(H_ab, H_ab.T)}")
    print(f"  M symmetric? {np.allclose(M, M.T)}")

# ===================================================================
# APPROACH 3: Irving-Omar det/per REFINED by endpoints
# ===================================================================
#
# Can we write ham(a->b) using an IO-like formula?
#
# For a Hamiltonian path from a to b:
#   ham(a->b) = sum over perm sigma with sigma(0)=a, sigma(n-1)=b
#               of prod A[sigma(i), sigma(i+1)]
#
# This is per(A with row a and col b fixed) in some sense.
# More precisely:
#   ham(a->b) = per(A[{a}∪U, U∪{b}])  where U = [n]\{a,b}
#            ... but this isn't quite right either.
#
# Actually, ham(a->b) counts directed Hamiltonian paths from a to b,
# which is the (a,b) entry of the "path permanent" of A.

print()
print("=" * 70)
print("IRVING-OMAR REFINED: per/det by endpoint")
print("=" * 70)

# For the IO formula, we need to understand what the
# NONCOMMUTING extraction L_n does.
#
# L_n picks out monomials x_{i_1}*x_{i_2}*...*x_{i_n} where
# {i_1,...,i_n} = [n]. So L_n counts WALKS that are HAMILTONIAN.
#
# In the NONCOMMUTING setting:
#   L_n(W_D(1)) = sum of products over all Hamiltonian walks
#               = ham(D) when we set all x_i = 1 AFTER extraction

# For the ENDPOINT-REFINED version:
#   L_n(x_a * ... * x_b) = Hamiltonian paths from a to b
#
# The IO formula at the matrix level:
#   W_D(z) = det(I + zXA^T) / det(I - zXA)
#
# At the noncommuting level, this is a formal power series in z
# whose coefficients are polynomials in the noncommuting x_i.
#
# The coefficient of z^{n-1} (for n vertices) after L_n extraction
# gives exactly ham(D).

# KEY INSIGHT: The transfer matrix M[a,b] has an inclusion-exclusion
# structure that looks like the EXPANSION of the det/per formula!
#
# det(A[S,S]) = sum over permutations of S, with sign
# per(A[S^c,S^c]) = sum over permutations of S^c, without sign
#
# Our M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
#
# The (-1)^|S| in our formula plays a role similar to the det in IO.
# The E_a * B_b product plays a role similar to the per in IO.
# But E_a and B_b count PATHS, while det and per count CYCLE COVERS.

print()
print("COMPARISON: IO det*per vs Transfer M inclusion-exclusion")
print()

for name, A in [("n=3 cyclic", A_cyc), ("n=3 trans", A_trans)]:
    n = len(A)
    A_np = np.array(A, dtype=float)
    A_T = A_np.T

    print(f"\n--- {name} ---")

    # IO decomposition by subset
    io_by_S = {}
    for mask in range(2**n):
        S = [i for i in range(n) if (mask >> i) & 1]
        Sc = [i for i in range(n) if not ((mask >> i) & 1)]

        if len(S) > 0:
            det_S = det_int(A_T[np.ix_(S, S)])  # det(A^T[S,S]) = det(A[S,S])
        else:
            det_S = 1

        if len(Sc) > 0:
            per_Sc = permanent(A_np[np.ix_(Sc, Sc)])
        else:
            per_Sc = 1

        io_by_S[tuple(S)] = (det_S, per_Sc, det_S * per_Sc)
        print(f"  S={S}: det(A[S])={det_S}, per(A[Sc])={per_Sc}, contribution={det_S*per_Sc}")

    total_io = sum(v[2] for v in io_by_S.values())
    print(f"  Total ham(IO) = {total_io}")

    # Our M decomposition
    M = transfer_matrix(A)
    print(f"  Transfer M = {M.tolist()}")

    # Show M[0,0] decomposition vs IO restricted to paths through 0
    print(f"\n  M[0,0] breakdown:")
    U = list(range(1, n))
    for k in range(len(U)+1):
        for S in combinations(U, k):
            S_set = set(S)
            R = [v for v in U if v not in S_set]
            S_verts = sorted(list(S) + [0])
            R_verts = sorted(R + [0])
            ea = count_paths_subset(A, S_verts, end=0)
            bb = count_paths_subset(A, R_verts, start=0)
            if ea * bb != 0:
                print(f"    S={list(S)}: (-1)^{k} * E_0({S_verts})={ea} * B_0({R_verts})={bb} = {(-1)**k * ea * bb}")

# ===================================================================
# APPROACH 4: Even r-powers in the Irving-Omar framework
# ===================================================================
print()
print("=" * 70)
print("EVEN R-POWERS AND THE IO FRAMEWORK")
print("=" * 70)
print()
print("THM-030: M[a,b] has only even powers of r in c-tournament.")
print("This means: M[a,b](r,s) = M[a,b](-r,s) for all r.")
print()
print("In the IO framework, the generating function is:")
print("  W_D(z) = det(I + zXA^T) / det(I - zXA)")
print()
print("For a c-tournament with t_{ij} = r + s_{ij}:")
print("  A has entries r + s_{ij} (for i->j)")
print("  A^T has entries r - s_{ij} (for j->i, since t_{ji} = r - s_{ij})")
print()
print("Under r -> -r:")
print("  A becomes -r + s_{ij} = -(r - s_{ij})")
print("  A^T becomes -r - s_{ij} = -(r + s_{ij})")
print("So A(-r) = -A^T(r) and A^T(-r) = -A(r)")
print()
print("Therefore:")
print("  W_D(-z, -r) = det(I + (-z)X(-A)) / det(I - (-z)X(-A^T))")
print("             = det(I - zXA) / det(I + zXA^T)")
print("             = 1/W_D(z, r)")
print()
print("This gives W(-z,-r) = 1/W(z,r).")
print("At the Hamiltonian coefficient level (z^{n-1}):")
print("  The coefficient of z^{n-1} in W(-z,-r) = (-1)^{n-1} * [z^{n-1}]W(z,-r)")
print("  And the coefficient of z^{n-1} in 1/W(z,r) involves convolution.")
print()
print("This is the IO-level manifestation of the THM-030 symmetry!")

# Let's verify: for c-tournament with parameter r, compute W(z) and check
# the even-r-power property

def T_weighted_from_to(s_vals, r_val, n, start, end):
    """Weighted Hamiltonian path count from start to end."""
    total = 0.0
    for p in permutations(range(n)):
        if p[0] != start or p[-1] != end:
            continue
        w = 1.0
        for i in range(n-1):
            a, b = p[i], p[i+1]
            t_ab = r_val + s_vals.get((a,b), 0)
            w *= t_ab
        total += w
    return total

def tournament_to_s(A):
    n = len(A)
    s = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                s[(i,j)] = A[i][j] - 0.5
    return s

# Verify even r-powers at the IO W(z) level
print()
print("--- Verification: W(z,r) * W(-z,-r) = 1 ---")

for name, A in [("n=3 cyclic", A_cyc), ("n=3 trans", A_trans)]:
    n = len(A)
    s = tournament_to_s(A)

    print(f"\n  {name}:")
    for r in [0.3, 0.5, 0.7, 1.0]:
        # Build A matrix for c-tournament
        A_r = np.zeros((n,n))
        A_r_T = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    A_r[i][j] = r + s.get((i,j), 0)
                    A_r_T[i][j] = r + s.get((j,i), 0)  # = r - s_{ij}

        I_n = np.eye(n)
        for z in [0.3, 0.7]:
            det_num = np.linalg.det(I_n + z * A_r_T)
            det_den = np.linalg.det(I_n - z * A_r)
            if abs(det_den) > 1e-10:
                W = det_num / det_den
            else:
                W = float('inf')

            # W(-z, -r): replace z->-z and r->-r
            A_mr = np.zeros((n,n))
            A_mr_T = np.zeros((n,n))
            for i in range(n):
                for j in range(n):
                    if i != j:
                        A_mr[i][j] = -r + s.get((i,j), 0)
                        A_mr_T[i][j] = -r + s.get((j,i), 0)

            det_num2 = np.linalg.det(I_n + (-z) * A_mr_T)
            det_den2 = np.linalg.det(I_n - (-z) * A_mr)
            if abs(det_den2) > 1e-10:
                W2 = det_num2 / det_den2
            else:
                W2 = float('inf')

            product = W * W2
            print(f"    r={r}, z={z}: W(z,r)={W:.4f}, W(-z,-r)={W2:.4f}, product={product:.6f}")

# ===================================================================
# APPROACH 5: Skeleton and the IO determinant structure
# ===================================================================
print()
print("=" * 70)
print("IO DET STRUCTURE ALONG THE TILING SKELETON")
print("=" * 70)
print()
print("When we flip a single tile (i,j) in the tiling:")
print("  A' = A + delta_A where delta_A has only 2 nonzero entries:")
print("  delta_A[i,j] changes by -1, delta_A[j,i] changes by +1 (or vice versa)")
print()
print("The IO generating function changes:")
print("  W' = det(I + zX(A')^T) / det(I - zXA')")
print()
print("Since delta_A is rank-2 (at most), we can use the matrix determinant lemma:")
print("  det(I - zXA') = det(I - zX(A + delta_A))")
print("               = det(I - zXA) * det(I - zX(I-zXA)^{-1} delta_A)")
print()
print("This gives a RATIO formula for how W changes under a single tile flip!")

# Compute for n=3
n = 3
A = A_cyc
A_np = np.array(A, dtype=float)
I_n = np.eye(n)

# Flip arc (2,0): currently 1, flip to (0,2)
A_flipped = [row[:] for row in A]
A_flipped[2][0] = 0
A_flipped[0][2] = 1
A_f = np.array(A_flipped, dtype=float)

print(f"\nn=3 cyclic: flip arc (2,0)")
print(f"  Original A = {A}")
print(f"  Flipped A' = {A_flipped}")

z = 0.5
# Original W
W_orig = np.linalg.det(I_n + z * A_np.T) / np.linalg.det(I_n - z * A_np)
W_flip = np.linalg.det(I_n + z * A_f.T) / np.linalg.det(I_n - z * A_f)

delta_A = A_f - A_np
print(f"  delta_A = {delta_A.tolist()}")
print(f"  z={z}: W_orig={W_orig:.6f}, W_flip={W_flip:.6f}, ratio={W_flip/W_orig:.6f}")

# The delta_A for a tournament arc flip is always rank 2:
# delta_A[i,j] = -1, delta_A[j,i] = +1
# (or the reverse)
rank = np.linalg.matrix_rank(delta_A)
print(f"  rank(delta_A) = {rank}")

# Sherman-Morrison for det change
# det(I - zX(A+dA)) / det(I - zXA) = det(I + (I-zXA)^{-1}*(-zX*dA))
R_inv = np.linalg.inv(I_n - z * A_np)
factor_den = np.linalg.det(I_n + R_inv @ (-z * delta_A))
print(f"  det ratio for denominator: {factor_den:.6f}")

R_inv2 = np.linalg.inv(I_n + z * A_np.T)
factor_num = np.linalg.det(I_n + R_inv2 @ (z * delta_A.T))
print(f"  det ratio for numerator: {factor_num:.6f}")
print(f"  W_flip/W_orig via det ratios: {factor_num/factor_den:.6f}")
print(f"  Match? {abs(W_flip/W_orig - factor_num/factor_den) < 1e-10}")

# ===================================================================
# APPROACH 6: What does M[a,b] count in terms of IO?
# ===================================================================
print()
print("=" * 70)
print("WHAT DOES M[a,b] COUNT IN THE IO FRAMEWORK?")
print("=" * 70)
print()

# M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
#
# E_a(S+a) = #paths through S+{a} ending at a
# B_b(R+b) = #paths through R+{b} starting at b
#
# The PRODUCT E_a * B_b counts ordered pairs of paths:
#   (path through S+{a} ending at a, path through R+{b} starting at b)
# These two paths together cover ALL vertices (a,b included).
# They form a "path decomposition" of the tournament.
#
# The inclusion-exclusion with (-1)^|S| converts this from an
# over-count to the exact count of... what?
#
# CLAIM: M[a,b] = ham(a->b) when a=b (loop, but only for odd n)
# and M[a,b] involves cancelation for a != b.
#
# Actually, for M[a,a] (diagonal):
# U = [n]\{a}, we split U into S and R=U\S.
# M[a,a] = sum_S (-1)^|S| E_a(S+a) B_a(R+a)
# Both paths go through vertex a.
# At |S|=0: E_a({a}) * B_a(U+a) = 1 * (all paths starting at a) = out-ham(a)
# At |S|=n-1: E_a([n]) * B_a({a}) = (all paths ending at a) * 1 = in-ham(a)
# The inclusion-exclusion corrects for overlaps.

for name, A in [("n=3 cyclic", A_cyc), ("n=3 trans", A_trans)]:
    n = len(A)
    M = transfer_matrix(A)
    H_ab = endpoint_ham_paths(A)

    print(f"\n--- {name} ---")
    print(f"  M[a,b]  = {M.tolist()}")
    print(f"  H(a->b) = {H_ab.tolist()}")

    # Compare row sums and column sums
    print(f"  M row sums = {[sum(M[i]) for i in range(n)]}")
    print(f"  M col sums = {[sum(M[:,j]) for j in range(n)]}")
    print(f"  H row sums = {[sum(H_ab[i]) for i in range(n)]}")
    print(f"  H col sums = {[sum(H_ab[:,j]) for j in range(n)]}")

    # Is there a simple matrix relationship?
    # M = H_ab @ P for some permutation matrix P?
    # Or M = f(H_ab) for some function f?

    # Check: does M = H_ab - something?
    diff = M - H_ab
    print(f"  M - H(a->b) = {diff.tolist()}")

# For n=4, let's check the same thing
print(f"\n--- n=4 ---")
A = A_4a
n = len(A)
M = transfer_matrix(A)
H_ab = endpoint_ham_paths(A)
H = ham_path_count(A)
print(f"  H = {H}")
print(f"  M[a,b]  = {M.tolist()}")
print(f"  H(a->b) = {H_ab.tolist()}")
print(f"  M - H(a->b) = {(M - H_ab).tolist()}")
print(f"  tr(M) = {np.trace(M)}, tr(H_ab) = {np.trace(H_ab)}")
print(f"  sum(M) = {M.sum()}, sum(H_ab) = {H_ab.sum()}")

# At n=4 (even): tr(M) should = 0, sum(M) = 2H
print(f"  Expected at even n: tr(M)=0, sum(M)=2H={2*H}")

print()
print("=" * 70)
print("KEY STRUCTURAL OBSERVATION")
print("=" * 70)
print("""
M[a,b] != H(a->b) in general!

M[a,b] is the "signed path decomposition" count, while H(a->b) is
the simple endpoint-refined Hamiltonian path count.

For the IO connection:
- ham(D) = L_n(W_D(1)) extracts Hamiltonian paths from the walk GF
- M[a,b] involves an inclusion-exclusion that seems related to the
  EXPANSION of det/per in IO, but with path-endpoint constraints

The even r-powers (THM-030) manifest in the IO framework as:
  W_D(z, -r) = 1/W_D(-z, r)  (for c-tournaments)

This means the IO generating function has a RECIPROCITY under
(z, r) -> (-z, -r), which is exactly the T <-> T^op symmetry
that underlies THM-030.
""")
