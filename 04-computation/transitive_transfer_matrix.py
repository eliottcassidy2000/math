#!/usr/bin/env python3
"""
Transfer matrix structure for TRANSITIVE tournaments.

The transitive tournament T_n has unique Ham path 0->1->2->...->n-1.
H(T_n) = 1.

M is tridiagonal with alternating diagonal: M[k,k] = (-1)^k.
Off-diagonal: M[k,k+1] = M[k+1,k] = (-1)^k (by symmetry).

What are the eigenvalues? What is the pattern?

kind-pasteur-2026-03-06-S25 (continuation)
"""

from itertools import permutations
import numpy as np

def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod != 0:
                pos = list(perm).index(a)
                val += (-1)**pos * prod
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S)
            S_set = set(S) | {a}
            R_set = set(R) | {b}
            ea = 0
            for p in permutations(sorted(S_set)):
                if p[-1] != a: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                ea += prod
            if len(S_set) == 1: ea = 1
            bb2 = 0
            for p in permutations(sorted(R_set)):
                if p[0] != b: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                bb2 += prod
            if len(R_set) == 1: bb2 = 1
            val += sign * ea * bb2
        return val


# ============================================================
# Compute M for transitive tournaments
# ============================================================
print("=" * 70)
print("Transfer matrix for transitive tournaments")
print("=" * 70)

for n in range(3, 8):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if i < j else 0

    M = np.zeros((n,n), dtype=int)
    for a in range(n):
        for b in range(n):
            M[a,b] = compute_M_entry(T, n, a, b)

    print(f"\n  n={n}:")
    print(f"    M =")
    for row in M:
        print(f"      {list(row)}")

    # Check tridiagonal
    is_tridiag = True
    for i in range(n):
        for j in range(n):
            if abs(i-j) > 1 and M[i,j] != 0:
                is_tridiag = False
    print(f"    Tridiagonal: {is_tridiag}")

    if is_tridiag:
        diag = [int(M[i,i]) for i in range(n)]
        super_diag = [int(M[i,i+1]) for i in range(n-1)]
        sub_diag = [int(M[i+1,i]) for i in range(n-1)]
        print(f"    Diagonal: {diag}")
        print(f"    Super-diagonal: {super_diag}")
        print(f"    Sub-diagonal: {sub_diag}")
        print(f"    Symmetric: {super_diag == sub_diag}")

        # Pattern: diagonal = (-1)^i, off-diagonal = (-1)^min(i,j)
        # Let's verify
        diag_pattern = [(-1)**i for i in range(n)]
        off_pattern = [(-1)**i for i in range(n-1)]
        print(f"    diag = (-1)^i: {diag == diag_pattern}")
        print(f"    off-diag = (-1)^i: {super_diag == off_pattern}")

    evals = sorted(np.linalg.eigvalsh(M))[::-1]
    print(f"    Eigenvalues: {[round(e, 6) for e in evals]}")
    print(f"    tr(M) = {int(np.trace(M))}")
    print(f"    det(M) = {int(round(np.linalg.det(M)))}")


# ============================================================
# OBSERVE: M for transitive tournament is related to Chebyshev
# ============================================================
print("\n" + "=" * 70)
print("Eigenvalue pattern analysis")
print("=" * 70)

print("""
For transitive T_n:
  M[i,i] = (-1)^i
  M[i,i+1] = M[i+1,i] = (-1)^i

Let D = diag((-1)^0, (-1)^1, ..., (-1)^{n-1}).
Then D * M * D should have a nice form since D^{-1} = D.

Actually, let M' = D * M:
  M'[i,j] = (-1)^i * M[i,j]

For tridiag: M'[i,i] = (-1)^i * (-1)^i = 1
             M'[i,i+1] = (-1)^i * (-1)^i = 1
             M'[i+1,i] = (-1)^{i+1} * (-1)^i = -1

So M' has diagonal 1, super-diagonal 1, sub-diagonal -1.
This is I + U - L where U = strict upper bidiagonal, L = strict lower bidiagonal.
""")

for n in [3, 5, 7]:
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if i < j else 0

    M = np.zeros((n,n), dtype=int)
    for a in range(n):
        for b in range(n):
            M[a,b] = compute_M_entry(T, n, a, b)

    D = np.diag([(-1)**i for i in range(n)])
    Mp = D @ M  # D * M

    print(f"\n  n={n}: D*M =")
    for row in Mp:
        print(f"    {[int(x) for x in row]}")

    # Check: is D*M = I + U - L?
    expected = np.eye(n, dtype=int)
    for i in range(n-1):
        expected[i, i+1] = 1
        expected[i+1, i] = -1
    print(f"    D*M = I + U - L: {np.array_equal(Mp, expected)}")

    # Eigenvalues of D*M
    evals = sorted(np.linalg.eigvals(Mp), key=lambda x: -abs(x))
    print(f"    Eigenvalues of D*M: {[round(complex(e).real, 4) + round(complex(e).imag, 4)*1j for e in evals]}")


# ============================================================
# Connection to skew-adjacency matrix
# ============================================================
print("\n" + "=" * 70)
print("Connection to skew-adjacency matrix")
print("=" * 70)

print("""
The skew-adjacency matrix S of a tournament T has:
  S[i,j] = 1 if i->j, S[i,j] = -1 if j->i, S[i,i] = 0

For the transitive tournament: S[i,j] = sgn(j-i) for i != j.
This is an antisymmetric matrix.

The matrix I + U - L (where U,L are bidiagonals) looks like:
  I + (upper bidiag 1's) - (lower bidiag 1's)
  = I + tridiag part of S

Is D*M equal to I + tridiag(S)?
""")

for n in [3, 5, 7]:
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if i < j else 0

    S = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i < j: S[i,j] = 1
            elif i > j: S[i,j] = -1

    M = np.zeros((n,n), dtype=int)
    for a in range(n):
        for b in range(n):
            M[a,b] = compute_M_entry(T, n, a, b)

    D = np.diag([(-1)**i for i in range(n)])

    # Tridiagonal part of S
    S_tridiag = np.zeros((n,n), dtype=int)
    for i in range(n):
        S_tridiag[i,i] = S[i,i]
        if i < n-1: S_tridiag[i,i+1] = S[i,i+1]
        if i > 0: S_tridiag[i,i-1] = S[i,i-1]

    print(f"\n  n={n}:")
    print(f"    D*M = I + tridiag(S): {np.array_equal(D @ M, np.eye(n, dtype=int) + S_tridiag)}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
