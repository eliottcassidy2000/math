# ⚠️ WARNING: This script assumed the transfer matrix M = [[1,0],[0,-1]],
# which is WRONG. Exhaustive check at n=4 shows widespread failures.
# M[a,b] values vary widely depending on the tournament (observed: -3 to +3 at n=5).
# The symmetry M[a,b] = M[b,a] IS true, but individual entries are NOT 0 or ±1.
# See MISTAKE-011a.

#!/usr/bin/env python3
"""
Verify the transfer matrix identity:
  sum_S (-1)^|S| E_a(S) * B_b(M\S) = delta_{ab} * (-1)^{a=j}

where a,b in {i,j}, M = V\{i,j}, and (i,j) is the fixed arc.

Specifically, defining:
  M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(M\S)

We found M = [[1, 0], [0, -1]] (rows/cols indexed by i,j).

This means:
  sum_S (-1)^|S| E_i(S)*B_i(R) = +1     (!)
  sum_S (-1)^|S| E_j(S)*B_j(R) = -1     (!)
  sum_S (-1)^|S| E_i(S)*B_j(R) = 0
  sum_S (-1)^|S| E_j(S)*B_i(R) = 0

The DIAGONAL entries being +1 and -1 is the real surprise!

Instance: opus-2026-03-05-S4b
"""
import random
from itertools import combinations

def h_end(T, n, S_list, v):
    """Count Hamiltonian paths in T[S_list union {v}] ending at v."""
    S = set(S_list)
    if not S:
        return 1
    S_sorted = sorted(S)
    m = len(S_sorted)
    full = (1 << m) - 1
    dp = [[0]*m for _ in range(1 << m)]
    for k, s in enumerate(S_sorted):
        dp[1 << k][k] = 1
    for mask in range(1, (1 << m)):
        for li in range(m):
            c = dp[mask][li]
            if c == 0: continue
            for ni in range(m):
                if mask & (1 << ni): continue
                if T[S_sorted[li]*n + S_sorted[ni]]:
                    dp[mask | (1 << ni)][ni] += c
    total = 0
    for li in range(m):
        if T[S_sorted[li]*n + v]:
            total += dp[full][li]
    return total

def h_start(T, n, R_list, v):
    """Count Hamiltonian paths in T[{v} union R_list] starting at v."""
    R = set(R_list)
    if not R:
        return 1
    R_sorted = sorted(R)
    m = len(R_sorted)
    full = (1 << m) - 1
    dp = [[0]*m for _ in range(1 << m)]
    for k, r in enumerate(R_sorted):
        if T[v*n + r]:
            dp[1 << k][k] = 1
    for mask in range(1, (1 << m)):
        for li in range(m):
            c = dp[mask][li]
            if c == 0: continue
            for ni in range(m):
                if mask & (1 << ni): continue
                if T[R_sorted[li]*n + R_sorted[ni]]:
                    dp[mask | (1 << ni)][ni] += c
    return sum(dp[full])

def random_tournament(n):
    T = [0]*(n*n)
    for a in range(n):
        for b in range(a+1, n):
            if random.random() < 0.5:
                T[a*n+b] = 1
            else:
                T[b*n+a] = 1
    return T

def compute_transfer_matrix(T, n, i, j):
    """Compute M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(M\S) for a,b in {i,j}."""
    M = [v for v in range(n) if v != i and v != j]
    m = len(M)

    result = [[0, 0], [0, 0]]  # [[M_ii, M_ij], [M_ji, M_jj]]

    for smask in range(1 << m):
        S = [M[k] for k in range(m) if smask & (1 << k)]
        R = [M[k] for k in range(m) if not (smask & (1 << k))]
        sz = len(S)
        sign = (-1) ** sz

        Ei = h_end(T, n, S, i)
        Ej = h_end(T, n, S, j)
        Bi = h_start(T, n, R, i)
        Bj = h_start(T, n, R, j)

        result[0][0] += sign * Ei * Bi
        result[0][1] += sign * Ei * Bj
        result[1][0] += sign * Ej * Bi
        result[1][1] += sign * Ej * Bj

    return result

print("=" * 60)
print("TRANSFER MATRIX IDENTITY VERIFICATION")
print("M[a,b] = sum_S (-1)^|S| E_a(S) * B_b(M\\S)")
print("Expected: M = [[+1, 0], [0, -1]]")
print("=" * 60)

for n in range(4, 9):
    print(f"\n--- n = {n} ---")
    random.seed(42)

    all_match = True
    trials = 500 if n <= 6 else (100 if n <= 7 else 30)

    for trial in range(trials):
        T = random_tournament(n)
        # Test all arc choices (i,j)
        for i in range(min(n, 3)):
            for j in range(i+1, min(n, i+3)):
                M = compute_transfer_matrix(T, n, i, j)

                # Check: does arc i->j exist?
                arc_ij = T[i*n + j]

                expected_ii = 1
                expected_jj = -1
                expected_ij = 0
                expected_ji = 0

                if M[0][0] != expected_ii or M[0][1] != expected_ij or \
                   M[1][0] != expected_ji or M[1][1] != expected_jj:
                    all_match = False
                    if trial < 3:
                        print(f"  UNEXPECTED trial={trial} arc ({i},{j}), T[i][j]={arc_ij}:")
                        print(f"    M = [[{M[0][0]}, {M[0][1]}], [{M[1][0]}, {M[1][1]}]]")

                        # Try with i,j swapped (j->i direction)
                        M2 = compute_transfer_matrix(T, n, j, i)
                        print(f"    M(j,i) = [[{M2[0][0]}, {M2[0][1]}], [{M2[1][0]}, {M2[1][1]}]]")

    if all_match:
        print(f"  ALL {trials} tournaments x arcs: M = [[+1, 0], [0, -1]]")
    else:
        print(f"  Some mismatches found!")

# Now check: does the sign depend on arc direction?
print("\n" + "=" * 60)
print("DOES THE IDENTITY DEPEND ON ARC DIRECTION?")
print("=" * 60)

random.seed(42)
for trial in range(5):
    n = 6
    T = random_tournament(n)
    i, j = 0, 1

    M_fwd = compute_transfer_matrix(T, n, i, j)
    M_rev = compute_transfer_matrix(T, n, j, i)

    print(f"\nTrial {trial}: T[{i}][{j}]={T[i*n+j]}")
    print(f"  M(i={i},j={j}) = [[{M_fwd[0][0]}, {M_fwd[0][1]}], [{M_fwd[1][0]}, {M_fwd[1][1]}]]")
    print(f"  M(i={j},j={i}) = [[{M_rev[0][0]}, {M_rev[0][1]}], [{M_rev[1][0]}, {M_rev[1][1]}]]")
