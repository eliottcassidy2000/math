#!/usr/bin/env python3
"""
Irving-Omar det/per formula and its connection to the transfer matrix.

Irving-Omar (arXiv:2412.10572, Proposition 2):
  ham(D) = sum_{S subset [n]} det(A_bar[S]) * per(A[S^c])

where A is the adjacency matrix of D, A_bar = J - I - A is the complement.

For tournaments: A_bar = A^T (the converse/transpose).

Our transfer matrix:
  M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)

Question: What is the PRECISE relationship between these two decompositions?

Also explore: Irving-Omar's walk generating function
  W_D(z) = det(I + z*X*A_bar) / det(I - z*X*A)
where X = diag(x_1, ..., x_n) are noncommuting variables.

At the commutative specialization x_i = 1, what does W_D(z) give?
"""

import numpy as np
from itertools import permutations, combinations
from math import factorial
from collections import defaultdict

def adj_matrix(T):
    """Tournament adjacency matrix."""
    return np.array(T)

def complement(A):
    """Complement adjacency: A_bar = J - I - A. For tournaments: A_bar = A^T."""
    n = len(A)
    J = np.ones((n,n)) - np.eye(n)
    return J - A

def submatrix(A, rows, cols):
    """Extract submatrix A[rows, cols]."""
    return A[np.ix_(rows, cols)]

def permanent(M):
    """Compute permanent of a matrix (brute force)."""
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

def det_np(M):
    """Determinant."""
    if len(M) == 0:
        return 1
    return round(np.linalg.det(M))

def count_paths_subset(A, verts, start=None, end=None):
    """Count Hamiltonian paths through verts."""
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

def transfer_matrix(A):
    """Our transfer matrix M[a,b]."""
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

def irving_omar_decomposition(A):
    """
    Compute ham(D) using Irving-Omar det/per formula:
    ham(D) = sum_{S subset [n]} det(A_bar[S]) * per(A[S^c])

    Returns: total ham count AND the per-subset contributions.
    """
    n = len(A)
    A = np.array(A, dtype=float)
    A_bar = complement(A)

    total = 0
    contributions = {}

    for mask in range(2**n):
        S = [i for i in range(n) if (mask >> i) & 1]
        Sc = [i for i in range(n) if not ((mask >> i) & 1)]

        if len(S) == 0:
            det_S = 1
        else:
            det_S = det_np(submatrix(A_bar, S, S))

        if len(Sc) == 0:
            per_Sc = 1
        else:
            per_Sc = permanent(submatrix(A, Sc, Sc))

        contrib = det_S * per_Sc
        contributions[tuple(S)] = (det_S, per_Sc, contrib)
        total += contrib

    return total, contributions

def irving_omar_walk_gf(A, z_val):
    """
    Evaluate W_D(z) = det(I + z*A_bar) / det(I - z*A) at scalar z.
    (Commutative specialization: X = I)
    """
    n = len(A)
    A = np.array(A, dtype=float)
    A_bar = complement(A)

    num = np.linalg.det(np.eye(n) + z_val * A_bar)
    den = np.linalg.det(np.eye(n) - z_val * A)

    if abs(den) < 1e-10:
        return float('inf')
    return num / den

print("=" * 70)
print("IRVING-OMAR DET/PER FORMULA vs TRANSFER MATRIX")
print("=" * 70)

# Test tournaments
# n=3 cyclic: 0->1->2->0
A_cyc = [[0,1,0],[0,0,1],[1,0,0]]
# n=3 transitive: 0->1, 0->2, 1->2
A_trans = [[0,1,1],[0,0,1],[0,0,0]]
# n=4 tournament
A_4 = [[0,1,1,0],[0,0,1,1],[0,0,0,1],[1,0,0,0]]

for name, A in [("n=3 cyclic", A_cyc), ("n=3 trans", A_trans), ("n=4", A_4)]:
    n = len(A)
    H_direct = count_paths_subset(A, list(range(n)))
    H_io, contribs = irving_omar_decomposition(A)
    M = transfer_matrix(A)

    print(f"\n--- {name} (H={H_direct}) ---")
    print(f"  Irving-Omar: ham = {H_io}")
    print(f"  Transfer matrix M = {M.tolist()}")
    print(f"  tr(M) = {np.trace(M)}, Sigma = {M.sum()-np.trace(M)}")

    # Show per-subset contributions
    print(f"  IO contributions by |S|:")
    by_size = defaultdict(float)
    for S, (det_s, per_sc, c) in contribs.items():
        by_size[len(S)] += c
    for k in sorted(by_size.keys()):
        print(f"    |S|={k}: total contribution = {by_size[k]}")

print()
print("=" * 70)
print("CONNECTING IO DECOMPOSITION TO M[a,b]")
print("=" * 70)
print()
print("The IO formula sums over ALL subsets S of [n].")
print("Our M[a,b] sums over subsets of [n]\\{a,b}.")
print("Can we express M[a,b] as an IO-like decomposition?")
print()

# For our M[a,b], let's rewrite using permanents/determinants
# E_a(S+a) = per of path matrix from S+a ending at a
# B_b(R+b) = per of path matrix from R+b starting at b
# These are NOT permanents of submatrices of A, but path permanents

# Actually, E_a(S+a) = sum over permutations of S+a that end at a
# = number of Hamiltonian paths through S+{a} ending at a
# In matrix terms: E_a(S) = sum_{sigma in S_{S+a}, sigma ends at a} prod A[sigma(i), sigma(i+1)]

# The permanent per(A[S,S]) counts CYCLE COVERS, not PATHS.
# So there's a fundamental difference: IO uses cycle covers (det/per),
# while our M uses path counts (E_a, B_b).

# However, the "walk generating function" approach might bridge this gap.
# W_D(z) at the noncommuting level tracks WALKS, and the coefficient
# extraction L_n picks out Hamiltonian paths.

print("KEY OBSERVATION:")
print("IO det/per uses CYCLE COVERS of subgraphs.")
print("Our M[a,b] uses PATH COUNTS through subsets.")
print("These are different combinatorial objects!")
print()
print("The connection is through the walk generating function W_D(z),")
print("which tracks all walks and whose appropriate coefficient gives ham(D).")

print()
print("=" * 70)
print("WALK GENERATING FUNCTION W(z) AT COMMUTATIVE LEVEL")
print("=" * 70)

# W_D(z) = det(I + z*A^T) / det(I - z*A) for tournaments
# At z=0: W(0) = 1
# At z=1: W(1) = det(I + A^T) / det(I - A)

# For tournaments, I - A has diagonal 1 and off-diagonal -(A_ij) or 0
# I + A^T = I + complement of A (since A_bar = A^T for tournaments)

for name, A in [("n=3 cyclic", A_cyc), ("n=3 trans", A_trans), ("n=4", A_4)]:
    n = len(A)
    A_np = np.array(A, dtype=float)
    A_bar = complement(A_np)

    H = count_paths_subset(A, list(range(n)))

    print(f"\n--- {name} ---")

    # Characteristic polynomials
    I_n = np.eye(n)

    # det(I - z*A) and det(I + z*A^T) as polynomials in z
    # For small n, compute at several z values
    z_vals = [0, 0.1, 0.2, 0.5, 1.0]
    for z in z_vals:
        W = irving_omar_walk_gf(A, z)
        det_num = np.linalg.det(I_n + z * A_bar)
        det_den = np.linalg.det(I_n - z * A_np)
        print(f"  z={z}: W(z)={W:.4f}, det(I+zA^T)={det_num:.4f}, det(I-zA)={det_den:.4f}")

    # At z=1 for tournament:
    # det(I - A) = det of matrix with 1s on diagonal, -A[i][j] off-diagonal
    # I - A for tournament has row sums = 1 - (outdegree of i)
    W_1 = irving_omar_walk_gf(A, 1.0)
    print(f"  W(1) = {W_1:.4f} (compare H = {H})")

print()
print("=" * 70)
print("WHAT DOES det(I - zA) ENCODE?")
print("=" * 70)
print()

# For a tournament A, det(I - zA) = 1 - z*tr(A) + z^2*...
# The coefficients are elementary symmetric functions of eigenvalues of A
# For tournaments on n vertices, tr(A) = 0 always

for name, A in [("n=3 cyclic", A_cyc), ("n=3 trans", A_trans), ("n=4", A_4)]:
    n = len(A)
    A_np = np.array(A, dtype=float)

    # Eigenvalues of A
    evals = np.linalg.eigvals(A_np)
    print(f"  {name}: eigenvalues of A = {[f'{e:.3f}' for e in evals]}")

    # Characteristic polynomial of A: det(zI - A)
    # = z^n - tr(A)*z^{n-1} + ...
    # For tournament: tr(A) = 0

    # det(I - zA) = det(I - zA) = prod(1 - z*lambda_i)
    # This is the RECIPROCAL characteristic polynomial evaluated at 1/z

    # For n=3 cyclic: eigenvalues are cube roots of unity scaled
    # A has eigenvalues summing to 0 (tr=0)

print()
print("=" * 70)
print("RATIO det(I+zA^T) / det(I-zA) FOR TOURNAMENTS")
print("=" * 70)
print()
print("For tournaments, A^T = J - I - A (complement).")
print("So I + z*A^T = I + z*(J - I - A) = (1-z)*I + z*J - z*A")
print("             = (1-z)*I + z*(J - A)")
print("And I - z*A is straightforward.")
print()
print("The ratio W(z) = det((1-z)I + z(J-A)) / det(I - zA)")

# For a regular tournament (doubly regular), the eigenvalues of A are
# {(n-1)/2, (-1 + i*sqrt(n))/2, (-1 - i*sqrt(n))/2, ...}

# Let's compute W(z) as a rational function in z for n=3
print()
print("--- Exact W(z) for n=3 cyclic ---")
# A = [[0,1,0],[0,0,1],[1,0,0]]
# eigenvalues of A: 1, omega, omega^2 where omega = e^{2pi*i/3}
# det(I - zA) = (1-z)(1-omega*z)(1-omega^2*z) = (1-z)(1-z*omega)(1-z*omega^2)
# = (1-z)(1 - z*(omega+omega^2) + z^2*omega^3) = (1-z)(1 + z + z^2)
# = 1 - z^3

# A^T = [[0,0,1],[1,0,0],[0,1,0]]
# eigenvalues same as A: 1, omega, omega^2
# det(I + zA^T) = (1+z)(1+omega*z)(1+omega^2*z) = (1+z)(1+z+z^2) = 1+z^3

# So W(z) = (1+z^3)/(1-z^3) for cyclic n=3

# Verify:
A_np = np.array(A_cyc, dtype=float)
for z in [0.1, 0.5, 0.9]:
    W = irving_omar_walk_gf(A_cyc, z)
    W_exact = (1 + z**3) / (1 - z**3)
    print(f"  z={z}: W(z)={W:.6f}, (1+z^3)/(1-z^3)={W_exact:.6f}, match={abs(W-W_exact)<1e-6}")

print(f"\n  W(z) = (1+z^3)/(1-z^3) for cyclic 3-tournament")
print("  Taylor: W(z) = 1 + 2z^3 + 2z^6 + ... = 1 + 2*sum_k z^(3k)")
print(f"  Coefficient of z^{2} (= n-1 = 2): 0")
print(f"  But ham(D) should be extracted as L_n W(1) in noncommuting setting")
print(f"  At commutative level, z^{2} coefficient is 0, but H = 3")
print(f"  => The commutative W(z) does NOT directly give H as a coefficient!")

print()
print("--- Exact W(z) for n=3 transitive ---")
# A = [[0,1,1],[0,0,1],[0,0,0]]
# A is nilpotent: A^3 = 0. Eigenvalues all 0.
# det(I - zA) = 1 (since A nilpotent)
# det(I + zA^T) = 1 (A^T also nilpotent)
# W(z) = 1/1 = 1 for all z

for z in [0.1, 0.5, 0.9]:
    W = irving_omar_walk_gf(A_trans, z)
    print(f"  z={z}: W(z)={W:.6f}")

print(f"  W(z) = 1 for transitive tournament (both A and A^T nilpotent)")
print(f"  But H = 1. So again, W(z) at commutative level doesn't give H directly.")

print()
print("=" * 70)
print("KEY INSIGHT: W(z) needs NONCOMMUTING variables X")
print("=" * 70)
print()
print("The commutative W(z) = det(I+zA^T)/det(I-zA) encodes CYCLE COVER")
print("information, not path information. To extract Hamiltonian paths,")
print("Irving-Omar use NONCOMMUTING variables X = diag(x_1,...,x_n)")
print("and extract the coefficient L_n (product of all x_i in some order).")
print()
print("Our transfer matrix M[a,b] directly encodes path endpoints.")
print("The connection may be: M[a,b] = coefficient of x_a * ... * x_b")
print("in some expansion of W_D, with appropriate noncommuting extraction.")
print()
print("This would give a MATRIX-ALGEBRAIC interpretation of M[a,b]!")

# Let me check: for n=3, does tr(M) relate to any expansion of W?
print()
for name, A in [("n=3 cyclic", A_cyc), ("n=3 trans", A_trans)]:
    n = len(A)
    M = transfer_matrix(A)
    H = count_paths_subset(A, list(range(n)))
    print(f"{name}: H={H}, tr(M)={np.trace(M)}, M={M.tolist()}")

print()
print("=" * 70)
print("DONE — Irving-Omar analysis")
print("=" * 70)
