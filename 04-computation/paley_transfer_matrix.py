#!/usr/bin/env python3
"""
Transfer matrix M for Paley tournaments.

DISCOVERY: For Paley T_7, M = 27*I (scalar matrix!).
Question: Does M = (H/n)*I for ALL vertex-transitive tournaments?
Or is this specific to Paley?

kind-pasteur-2026-03-06-S25
"""

from itertools import permutations
import numpy as np

def make_paley(p):
    """Paley tournament on Z_p for prime p = 3 mod 4."""
    qr = set()
    for k in range(1, p):
        qr.add((k*k) % p)
    T = {}
    for i in range(p):
        for j in range(p):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % p) in qr else 0
    return T

def compute_M_numeric(T, n):
    """Compute full transfer matrix M[a,b] for a numeric tournament."""
    M = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            if a == b:
                val = 0
                for perm in permutations(range(n)):
                    prod = 1
                    for k in range(len(perm)-1):
                        prod *= T.get((perm[k], perm[k+1]), 0)
                    if prod != 0:
                        pos = list(perm).index(a)
                        val += (-1)**pos * prod
                M[a,a] = val
            else:
                val = 0
                U = [v for v in range(n) if v != a and v != b]
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
                M[a,b] = val
    return M

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(len(perm)-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count


# ============================================================
# Test Paley tournaments
# ============================================================
print("=" * 70)
print("Transfer matrix for Paley tournaments")
print("=" * 70)

for p in [3, 5, 7]:
    if p == 5:
        print(f"\n  p=5: p = 1 mod 4, not Paley tournament (skip)")
        continue

    T = make_paley(p)
    M = compute_M_numeric(T, p)
    H = count_H(T, p)

    print(f"\n  Paley T_{p}: H = {H}")
    print(f"    M = ")
    for row in M:
        print(f"      {[int(x) for x in row]}")

    # Check if M = c*I
    diag = M[0, 0]
    is_scalar = np.allclose(M, diag * np.eye(p))
    print(f"    M = {int(diag)}*I ? {is_scalar}")
    if is_scalar:
        print(f"    H/n = {H}/{p} = {H/p}, M[a,a] = {int(diag)}")
        print(f"    Match: {abs(H/p - diag) < 0.001}")

    # Eigenvalues
    evals = np.linalg.eigvalsh(M)
    print(f"    eigenvalues = {sorted(set(int(round(e)) for e in evals))}")


# ============================================================
# Test ALL regular tournaments at n=5 and n=7
# ============================================================
print("\n" + "=" * 70)
print("Transfer matrix for regular tournaments")
print("=" * 70)

import random
random.seed(123)

# n=5: all regular tournaments have score (2,2,2,2,2)
# There are 24 regular tournaments on 5 vertices (up to labeling more)
# Let's generate some and check

def is_regular(T, n):
    for i in range(n):
        out_deg = sum(T.get((i,j), 0) for j in range(n) if j != i)
        if out_deg != (n-1)//2:
            return False
    return True

# Generate random regular n=5 tournaments
print("\nn=5 regular tournaments:")
found = set()
for _ in range(100000):
    T = {}
    pairs = [(i,j) for i in range(5) for j in range(i+1, 5)]
    for (i,j) in pairs:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    if not is_regular(T, 5):
        continue

    # Canonical form: tuple of edge directions
    key = tuple(T.get((i,j), 0) for (i,j) in pairs)
    if key in found:
        continue
    found.add(key)

    M = compute_M_numeric(T, 5)
    H = count_H(T, 5)

    # Check if M = (H/n)*I
    diag = M[0, 0]
    is_scalar = np.allclose(M, diag * np.eye(5))
    off_diag = set(int(round(M[i,j])) for i in range(5) for j in range(5) if i != j)

    if is_scalar:
        print(f"  H={H}: M = {int(diag)}*I (SCALAR!)")
    else:
        print(f"  H={H}: M is NOT scalar. Diagonal: {set(int(round(M[i,i])) for i in range(5))}, Off-diag: {off_diag}")

print(f"  Total regular tournaments found: {len(found)}")


# ============================================================
# Vertex-transitive but non-Paley at n=7
# ============================================================
print("\n" + "=" * 70)
print("Non-Paley vertex-transitive tournaments at n=7")
print("=" * 70)

# The doubly-regular tournament on Z_7 with QR = {1,2,4} is the Paley one.
# The other DRT has QR = {1,3,5} (non-QR set for mod 7... actually {1,2,4} is QR mod 7)
# Let me construct the other circulant tournament with S = {1,3,5}

for S_set_name, S_set in [("QR={1,2,4} (Paley)", {1,2,4}), ("S={1,3,5}", {1,3,5}), ("S={2,3,4}", {2,3,4})]:
    T = {}
    n = 7
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S_set else 0

    if not is_regular(T, n):
        continue

    M = compute_M_numeric(T, n)
    H = count_H(T, n)

    diag = M[0, 0]
    is_scalar = np.allclose(M, diag * np.eye(n))

    print(f"\n  Circulant {S_set_name}: H = {H}")
    if is_scalar:
        print(f"    M = {int(diag)}*I (SCALAR!)")
    else:
        diag_vals = sorted(set(int(round(M[i,i])) for i in range(n)))
        off_vals = sorted(set(int(round(M[i,j])) for i in range(n) for j in range(n) if i != j))
        print(f"    M diagonal: {diag_vals}, off-diag: {off_vals}")
        evals = sorted(np.linalg.eigvalsh(M))[::-1]
        print(f"    eigenvalues: {[round(e, 2) for e in evals]}")


# ============================================================
# What about n=3 Paley?
# ============================================================
print("\n" + "=" * 70)
print("Paley T_3")
print("=" * 70)

T = make_paley(3)
M = compute_M_numeric(T, 3)
H = count_H(T, 3)
print(f"  H = {H}")
print(f"  M = ")
for row in M:
    print(f"    {[int(x) for x in row]}")
is_scalar = np.allclose(M, M[0,0] * np.eye(3))
print(f"  M = {int(M[0,0])}*I ? {is_scalar}")


# ============================================================
# Test: ALL vertex-transitive n=5 tournaments
# ============================================================
print("\n" + "=" * 70)
print("Vertex-transitive (circulant) n=5 tournaments")
print("=" * 70)

# On Z_5, tournament defined by S subset {1,2,3,4} with |S|=2 and j in S iff 5-j not in S
# Possible: S={1,2}, S={1,3}, S={2,4}, S={3,4}
# But S={1,2} and S={3,4} are reverses (T^op). Same for S={1,3} and S={2,4}.

for S_name, S in [("S={1,2}", {1,2}), ("S={1,3}", {1,3})]:
    T = {}
    n = 5
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j-i) % n) in S else 0

    M = compute_M_numeric(T, n)
    H = count_H(T, n)

    diag = M[0, 0]
    is_scalar = np.allclose(M, diag * np.eye(n))

    print(f"\n  Circulant {S_name}: H = {H}")
    if is_scalar:
        print(f"    M = {int(diag)}*I (SCALAR!)")
    else:
        print(f"    M =")
        for row in M:
            print(f"      {[int(x) for x in row]}")
        evals = sorted(np.linalg.eigvalsh(M))[::-1]
        print(f"    eigenvalues: {[round(e, 2) for e in evals]}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
