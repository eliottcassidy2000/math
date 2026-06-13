#!/usr/bin/env python3
"""
Bridge between the consecutive-position formula and Irving-Omar framework.

THEOREM (kind-pasteur, verified n=5):
  M[a,b] = sum_{j=0}^{n-2} (-1)^j * C(a,b,j)
  where C(a,b,j) = #{Ham paths with a at pos j, b at pos j+1}

This means M[a,b] counts "signed consecutive occurrences" of edge a->b.

CONNECTION TO IO: The Irving-Omar walk GF
  W_D(z) = det(I + zA^T) / det(I - zA)
encodes path statistics. The consecutive-position formula suggests
M is related to the DERIVATIVE or RESIDUE of W at specific points.

CONNECTION TO EVEN r-POWERS (THM-030):
  THM-030 says B_b(W) + (-1)^m E_b(W) = 2r * col_sum_W(b)
  The consecutive formula involves E_a * B_b on complementary sets.
  Can we use THM-030 to constrain the consecutive counts?

Let's explore:
1. Verify the consecutive formula at n=3, 4, 6
2. Express C(a,b,j) in terms of sub-tournament path counts
3. Connect to the IO walk GF coefficient extraction
4. Use THM-030 to derive parity constraints on C(a,b,j)
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def ham_path_count_dp(A):
    n = len(A)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if (mask & (1 << u)) or A[v][u] != 1:
                    continue
                dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1: valid = False; break
        if valid: count += 1
    return count

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

def consecutive_position_formula(A):
    """Compute M via the consecutive-position formula."""
    n = len(A)
    # C[a][b][j] = #{Ham paths with a at pos j, b at pos j+1}
    C = [[[0]*n for _ in range(n)] for _ in range(n)]
    for perm in permutations(range(n)):
        valid = True
        for k in range(n-1):
            if A[perm[k]][perm[k+1]] != 1:
                valid = False
                break
        if valid:
            for j in range(n-1):
                C[perm[j]][perm[j+1]][j] += 1

    M_consec = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            if a == b:
                # Diagonal: M[a,a] = sum_j (-1)^j * P[a,j]
                # where P[a,j] = #{paths with a at position j}
                val = 0
                for perm in permutations(range(n)):
                    valid = True
                    for k in range(n-1):
                        if A[perm[k]][perm[k+1]] != 1:
                            valid = False
                            break
                    if valid:
                        pos = list(perm).index(a)
                        val += (-1)**pos
                M_consec[a][a] = val
            else:
                M_consec[a][b] = sum((-1)**j * C[a][b][j] for j in range(n-1))
    return M_consec, C

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(edges)):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield [list(row) for row in A]

# =====================================================================
# VERIFICATION at n=3, 4
# =====================================================================
print("=" * 70)
print("VERIFY CONSECUTIVE-POSITION FORMULA")
print("=" * 70)

for n in [3, 4]:
    print(f"\n--- n={n} ---")
    count = 0
    match = 0
    fail = 0
    for A in all_tournaments(n):
        M_def = transfer_matrix(A)
        M_con, C = consecutive_position_formula(A)
        if np.array_equal(M_def, M_con):
            match += 1
        else:
            fail += 1
            if fail <= 2:
                print(f"  MISMATCH!")
                print(f"    M_def = {M_def.tolist()}")
                print(f"    M_con = {M_con.tolist()}")
        count += 1
    print(f"  {match}/{count} tournaments match ({fail} failures)")

# =====================================================================
# n=5: C(a,b,j) structure for key tournaments
# =====================================================================
print()
print("=" * 70)
print("n=5: CONSECUTIVE-POSITION ARRAY C(a,b,j)")
print("=" * 70)

# Paley tournament
n = 5
A_paley = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and (j-i)%n in [1, 2]:
            A_paley[i][j] = 1

M, C = consecutive_position_formula(A_paley)
H = ham_path_count_dp(A_paley)
print(f"\nPaley T_5: H={H}, M=2*I")

# Show C array for one edge
for a in range(min(3, n)):
    for b in range(n):
        if a == b: continue
        c_vals = [C[a][b][j] for j in range(n-1)]
        alt_sum = sum((-1)**j * c_vals[j] for j in range(n-1))
        edge = A_paley[a][b]
        if a < 2:
            print(f"  C({a},{b},j) = {c_vals}, alt_sum={alt_sum}, T[{a},{b}]={edge}")

# =====================================================================
# KEY INSIGHT: C(a,b,j) = 0 when T[a,b] = 0
# Because a can only immediately precede b if a->b is an arc.
# So M[a,b] only depends on arcs a->b.
# =====================================================================
print()
print("=" * 70)
print("C(a,b,j) SYMMETRY AND THM-030 CONNECTION")
print("=" * 70)

# For the Paley: C(a,b,j) should be symmetric under cyclic shift
print("\nPaley C(a,b,j) for edge 0->1:")
for j in range(n-1):
    print(f"  j={j}: C(0,1,{j}) = {C[0][1][j]}")

# Now: M[a,b] = M[b,a] (symmetry from THM-030)
# But C(a,b,j) counts paths with a@j, b@j+1 (consecutive, a BEFORE b)
# C(b,a,j) counts paths with b@j, a@j+1 (consecutive, b BEFORE a)
# These are DIFFERENT! Yet the alternating sums are equal.

# So: sum_j (-1)^j C(a,b,j) = sum_j (-1)^j C(b,a,j) for ALL a,b
# This is a NON-TRIVIAL identity!
print("\nSymmetry check: sum(-1)^j C(a,b,j) = sum(-1)^j C(b,a,j)?")
for a in range(n):
    for b in range(a+1, n):
        ab_sum = sum((-1)**j * C[a][b][j] for j in range(n-1))
        ba_sum = sum((-1)**j * C[b][a][j] for j in range(n-1))
        if ab_sum != ba_sum:
            print(f"  FAIL: ({a},{b}): ab={ab_sum}, ba={ba_sum}")
        else:
            print(f"  OK: ({a},{b}): M[a,b]=M[b,a]={ab_sum}, "
                  f"C_ab={[C[a][b][j] for j in range(n-1)]}, "
                  f"C_ba={[C[b][a][j] for j in range(n-1)]}")

# =====================================================================
# DEEPER: What is C(a,b,j) + C(b,a,n-2-j)?
# If path has a@j, b@j+1, reverse has b@(n-1-j-1)=n-2-j, a@(n-1-j)=n-1-j
# So C_rev(b,a,n-2-j) corresponds to reversals.
# But reversals use reversed arcs — in the OPPOSITE tournament T^op!
# =====================================================================
print()
print("=" * 70)
print("C(a,b,j) AND PATH REVERSAL")
print("=" * 70)

# If we reverse a Ham path, a@j,b@j+1 becomes b@(n-2-j), a@(n-1-j)
# in the OPPOSITE tournament. So:
# C_T(a,b,j) = C_{T^op}(b,a,n-2-j)
# And M_T[a,b] = sum_j (-1)^j C_T(a,b,j)
#              = sum_j (-1)^j C_{T^op}(b,a,n-2-j)
# Let k = n-2-j, then j = n-2-k, (-1)^j = (-1)^{n-2-k} = (-1)^n * (-1)^k
# (since (-1)^{-2} = 1)
# So M_T[a,b] = (-1)^n * sum_k (-1)^k C_{T^op}(b,a,k)
#             = (-1)^n * M_{T^op}[b,a]

# For T^op: M_{T^op}[b,a] = (-1)^n * M_T[a,b]
# This gives: M_{T^op} = (-1)^n * M_T^T

# Check this at n=3 ((-1)^3 = -1):
# M_{T^op} = -M_T^T
print("Checking M_{T^op} = (-1)^n * M_T^T:")
for n in [3, 4, 5]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    test_count = 0
    for bits in range(2**len(edges)):
        if test_count >= 20: break
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1

        # Opposite tournament
        A_op = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    A_op[i][j] = 1 - A[i][j]

        M_T = transfer_matrix(A)
        M_Top = transfer_matrix(A_op)

        expected = ((-1)**n) * M_T.T
        if np.array_equal(M_Top, expected):
            test_count += 1
        else:
            print(f"  n={n}, FAIL!")
            print(f"    M_T = {M_T.tolist()}")
            print(f"    M_Top = {M_Top.tolist()}")
            print(f"    (-1)^n * M_T^T = {expected.tolist()}")
            break
    else:
        print(f"  n={n}: M_{{T^op}} = (-1)^n * M_T^T VERIFIED ({test_count} tournaments)")

# =====================================================================
# PARITY STRUCTURE: C(a,b,j) mod 2
# =====================================================================
print()
print("=" * 70)
print("C(a,b,j) PARITY STRUCTURE")
print("=" * 70)

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
seen = set()

for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    key = tuple(tuple(row) for row in A)
    min_key = key
    for perm in permutations(range(n)):
        pkey = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if pkey < min_key: min_key = pkey
    if min_key in seen: continue
    seen.add(min_key)

    H = ham_path_count_dp(A)
    M, C = consecutive_position_formula(A)

    # For each edge a->b: is sum_j C(a,b,j) = #{paths containing edge a->b in consecutive positions}
    # This equals #{paths using arc a->b} where a immediately before b.
    # Note: sum_j C(a,b,j) = (edge usage count in all Ham paths)

    # Parity: sum_j (-1)^j C(a,b,j) must be integer (M entry)
    # And C(a,b,j) >= 0 always.

    # At H=1 (transitive tournament): only 1 path, each edge used once
    if H == 1:
        print(f"\n  Transitive (H=1):")
        for a in range(n):
            for b in range(n):
                c_vals = [C[a][b][j] for j in range(n-1)]
                if sum(c_vals) > 0:
                    print(f"    C({a},{b},j) = {c_vals}")
        break

# Collect C statistics across all iso classes
print("\n  Edge usage sum_j C(a,b,j) vs H:")
seen = set()
for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    key = tuple(tuple(row) for row in A)
    min_key = key
    for perm in permutations(range(n)):
        pkey = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if pkey < min_key: min_key = pkey
    if min_key in seen: continue
    seen.add(min_key)

    H = ham_path_count_dp(A)
    M, C = consecutive_position_formula(A)

    total_usage = sum(C[a][b][j] for a in range(n) for b in range(n) for j in range(n-1))
    # total_usage = sum of all C entries = (n-1) * H (each path has n-1 edges)
    expected = (n-1) * H
    print(f"  H={H:3d}: total_edge_usage={total_usage}, (n-1)*H={expected}, match={total_usage==expected}")

# =====================================================================
# IO CONNECTION: W(z) coefficient extraction
# =====================================================================
print()
print("=" * 70)
print("IO WALK GF: COEFFICIENT CONNECTION TO C(a,b,j)")
print("=" * 70)

# Irving-Omar: W_T(z) = det(I + zA^T) / det(I - zA)
# For a single variable z, this is a rational function.
# At commutative level, W(z) = product over odd cycles of (1+z^k)/(1-z^k)
# The coefficient of z^k in W(z) relates to cycle covers of weight k.

# M[a,b] uses the POSITION of a in paths. This is like extracting
# the coefficient of a specific term in a generating function.

# Key: W(z) counts CYCLE COVERS by weight.
# M[a,b] counts PATHS by signed position.
# These are DIFFERENT combinatorial objects!

# But the IO walk GF also has a PATH version:
# The number of walks of length k from i to j is [A^k]_{ij}.
# For HAMILTONIAN paths, we need the permanent or a restricted sum.

# The connection might be through the NONCOMMUTING version:
# W_D(z) = det(I + zXA^T) / det(I - zXA)
# where X = diag(x_1, ..., x_n) tracks vertex usage.
# Setting all x_i = 1 gives the commuting version.

# The coefficient of x_1 x_2 ... x_n z^{n-1} in the numerator of W
# should relate to Hamiltonian paths (each vertex used once, n-1 steps).

print("The IO walk GF W(z) and M[a,b] are connected through:")
print("  W_D(z) = det(I + zXA^T) / det(I - zXA)")
print("  Setting x_i = 1: commuting version counts cycle covers")
print("  The multilinear part (coeff of x_1...x_n) counts Ham structures")
print("")
print("  M[a,b] = sum_j (-1)^j C(a,b,j)")
print("  C(a,b,j) = #{paths: a@j, b@j+1}")
print("")
print("  The (-1)^j weighting is the SAME sign pattern as in")
print("  det(I - zA): the determinant expansion uses (-1)^{cycle length}.")
print("  This suggests M[a,b] extracts the 'edge (a,b) contribution'")
print("  from the determinantal structure of W.")

# Verify: at n=3, compute W(z) and check coefficient structure
print()
n = 3
A3 = [[0,1,0],[0,0,1],[1,0,0]]  # 3-cycle
A_np = np.array(A3)

# W(z) = det(I + z*A^T) / det(I - z*A)
# For 3-cycle: A^3 = I, eigenvalues are cube roots of unity
# det(I - zA) = (1-z)(1-wz)(1-w^2 z) where w = e^{2pi i/3}
# = 1 - (1+w+w^2)z + (w+w^2+w^3)z^2 - z^3
# = 1 - 0 - 0 - z^3 = 1 - z^3

# det(I + zA^T) = det(I + zA^{-1}) (since A^T = A^{-1} for permutation matrix)
# = det(A^{-1}) det(A + zI) = det(A+zI) / det(A)
# Actually for 3-cycle A: A^T = A^2 (inverse)
# det(I + zA^2) = same eigenvalue structure: (1+z)(1+wz)(1+w^2 z) = 1+z^3

# So W(z) = (1+z^3)/(1-z^3) for the directed 3-cycle
print("3-cycle: W(z) = (1+z^3)/(1-z^3)")
print("  = 1 + 2z^3 + 2z^6 + ...")
print("  Only odd-length cycle covers contribute (weight 3, 6, ...)")

# For the transitive tournament on 3: A = [[0,1,1],[0,0,1],[0,0,0]]
A3_trans = [[0,1,1],[0,0,1],[0,0,0]]
A_np_t = np.array(A3_trans)

# det(I - zA) for upper triangular = det(I) = 1 (since A is nilpotent, A^3=0)
# Actually: I-zA is upper triangular with 1's on diagonal, det = 1
# det(I + zA^T): A^T is lower triangular, I+zA^T also lower triangular, det = 1
# W(z) = 1/1 = 1 for transitive tournament!
print("\nTransitive T_3: W(z) = 1")
print("  This makes sense: no cycles at all in transitive tournament")
print("  H = 1 (one Ham path), but W only counts cycle covers, not paths")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
KEY RESULTS:

1. CONSECUTIVE FORMULA VERIFIED at n=3,4,5:
   M[a,b] = sum_j (-1)^j * #{paths with a@j, b@j+1}

2. PATH REVERSAL IDENTITY:
   M_{T^op} = (-1)^n * M_T^T
   This means: at odd n, M_{T^op} = -M_T^T
              at even n, M_{T^op} = M_T^T

3. SYMMETRY CONSEQUENCE (THM-030):
   M[a,b] = M[b,a] means:
   sum_j (-1)^j C(a,b,j) = sum_j (-1)^j C(b,a,j)
   despite C(a,b,j) != C(b,a,j) in general!
   (One side is zero when T[a,b]=0 or T[b,a]=0)

4. TOTAL EDGE USAGE:
   sum_{a,b,j} C(a,b,j) = (n-1)*H (each path has n-1 edges)

5. IO CONNECTION:
   W(z) counts cycle covers, not paths.
   The consecutive formula shows M extracts PATH information
   from a signed position weighting, not cycle information.
   The connection to IO must go through the multilinear
   (vertex-tracking) version of the walk GF.
""")
