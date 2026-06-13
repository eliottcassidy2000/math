#!/usr/bin/env python3
"""
COROLLARY: Transfer matrix of vertex-transitive tournaments.

For vertex-transitive T on n vertices (odd n):
  M = (H(T)/n) * I

For vertex-transitive T on n vertices (even n):
  M = (2H(T))/(n(n-1)) * (J - I)

PROOF: By THM-030, M is symmetric. By THM-027, tr(M) = H for odd n, 0 for even n.
For odd n, off-diag sum = 0 (from even r-powers + Sigma identity).
For even n, off-diag sum = 2H.
Vertex-transitivity forces all diagonal entries equal and all off-diag entries equal.

kind-pasteur-2026-03-06-S25
"""

from itertools import permutations
import numpy as np

def make_circulant(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(len(perm)-1):
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

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(len(perm)-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count


# ============================================================
# ODD n: vertex-transitive => M = (H/n)*I
# ============================================================
print("=" * 70)
print("ODD n: Vertex-transitive tournaments => M = (H/n)*I")
print("=" * 70)

# n=3: only Paley T_3
T = make_circulant(3, {1})
n = 3
H = count_H(T, n)
M00 = compute_M_entry(T, n, 0, 0)
M01 = compute_M_entry(T, n, 0, 1)
print(f"  n=3 (Paley): H={H}, M[0,0]={M00}, M[0,1]={M01}, H/n={H/n}")
print(f"    M = {M00}*I: {M01 == 0 and M00 == H/n}")

# n=5: two circulant classes
for S_name, S in [("{1,2}", {1,2}), ("{1,3}", {1,3})]:
    T = make_circulant(5, S)
    n = 5
    H = count_H(T, n)
    M00 = compute_M_entry(T, n, 0, 0)
    M01 = compute_M_entry(T, n, 0, 1)
    print(f"  n=5 S={S_name}: H={H}, M[0,0]={M00}, M[0,1]={M01}, H/n={H/n}")
    print(f"    M = {M00}*I: {M01 == 0 and M00 == H/n}")

# n=7: three circulant classes
for S_name, S in [("{1,2,4} (Paley)", {1,2,4}), ("{1,3,5}", {1,3,5}), ("{2,3,4}", {2,3,4})]:
    T = make_circulant(7, S)
    n = 7
    H = count_H(T, n)
    M00 = compute_M_entry(T, n, 0, 0)
    M01 = compute_M_entry(T, n, 0, 1)
    print(f"  n=7 S={S_name}: H={H}, M[0,0]={M00}, M[0,1]={M01}, H/n={H/n:.1f}")
    print(f"    M = {M00}*I: {M01 == 0 and abs(M00 - H/n) < 0.001}")


# ============================================================
# EVEN n: vertex-transitive => M = c*(J-I)
# ============================================================
print("\n" + "=" * 70)
print("EVEN n: Vertex-transitive tournaments")
print("Expected: M = (2H)/(n(n-1)) * (J - I)")
print("=" * 70)

# n=4: circulant tournament (only one VT class: S={1,2})
for S_name, S in [("{1,2}", {1,2}), ("{1,3}", {1,3})]:
    T = make_circulant(4, S)
    n = 4
    H = count_H(T, n)
    M_full = np.zeros((n,n))
    for a in range(n):
        for b in range(n):
            M_full[a,b] = compute_M_entry(T, n, a, b)
    c_expected = 2*H / (n*(n-1))
    print(f"\n  n=4 S={S_name}: H={H}")
    print(f"    M =")
    for row in M_full:
        print(f"      {[int(x) for x in row]}")
    print(f"    Expected c = 2H/(n(n-1)) = {c_expected:.4f}")
    print(f"    Diagonal: {set(int(M_full[i,i]) for i in range(n))}")
    print(f"    Off-diag: {set(int(M_full[i,j]) for i in range(n) for j in range(n) if i != j)}")
    # Check M = c*(J-I)
    J_minus_I = np.ones((n,n)) - np.eye(n)
    is_match = np.allclose(M_full, c_expected * J_minus_I)
    print(f"    M = c*(J-I)? {is_match}")

# n=6: circulant tournament
for S_name, S in [("{1,2,3}", {1,2,3}), ("{1,2,4}", {1,2,4}), ("{1,3,5}", {1,3,5})]:
    T = make_circulant(6, S)
    n = 6
    H = count_H(T, n)
    # Just check a few entries
    M00 = compute_M_entry(T, n, 0, 0)
    M01 = compute_M_entry(T, n, 0, 1)
    M02 = compute_M_entry(T, n, 0, 2)
    c_expected = 2*H / (n*(n-1))
    print(f"\n  n=6 S={S_name}: H={H}")
    print(f"    M[0,0]={M00}, M[0,1]={M01}, M[0,2]={M02}")
    print(f"    Expected: diag=0, off-diag={c_expected:.4f}")
    print(f"    Match: diag=0? {M00==0}, all off-diag equal? {M01==M02}")


# ============================================================
# The eigenvalue formula for VT tournaments
# ============================================================
print("\n" + "=" * 70)
print("Eigenvalue formulas for vertex-transitive transfer matrix")
print("=" * 70)

print("""
For vertex-transitive T on n vertices:

ODD n:
  M = (H/n) * I
  Eigenvalue: H/n with multiplicity n
  det(M) = (H/n)^n

EVEN n:
  M = c * (J - I)  where c = 2H/(n(n-1))
  Eigenvalues: c*(n-1) = 2H/n (multiplicity 1)
               -c = -2H/(n(n-1)) (multiplicity n-1)
  det(M) = (-1)^{n-1} * c^n * (n-1)
         = (-1)^{n-1} * [2H/(n(n-1))]^n * (n-1)
""")

# Verify eigenvalue formulas
for n, S_name, S in [(3, "{1}", {1}), (5, "{1,2}", {1,2}), (7, "{1,2,4}", {1,2,4}),
                      (4, "{1,2}", {1,2}), (6, "{1,2,3}", {1,2,3})]:
    T = make_circulant(n, S)
    H = count_H(T, n)
    M_full = np.zeros((n,n))
    for a in range(n):
        for b in range(n):
            M_full[a,b] = compute_M_entry(T, n, a, b)

    evals = sorted(np.linalg.eigvalsh(M_full))[::-1]

    if n % 2 == 1:
        expected_eval = H / n
        print(f"  n={n} (odd): H={H}, expected eval={expected_eval}, actual={evals}")
    else:
        c = 2 * H / (n * (n-1))
        e1 = c * (n-1)
        e2 = -c
        print(f"  n={n} (even): H={H}, c={c:.4f}")
        print(f"    Expected: {e1:.4f} (x1), {e2:.4f} (x{n-1})")
        print(f"    Actual: {[round(e,4) for e in evals]}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
