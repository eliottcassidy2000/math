#!/usr/bin/env python3
"""
pfaffian_H_connection.py — Does Pfaffian of B = 2A-J+I relate to H(T)?

At n=6 (even n), B is 6x6 skew-symmetric, so Pf(B)^2 = det(B).
The Pfaffian counts perfect matchings of a graph weighted by B.

B[i][j] = 2*A[i][j] - 1 for i≠j (so B[i][j] = +1 if i→j, -1 if j→i)
B[i][i] = 0 (skew-symmetric diagonal)

Wait: B = 2A - J + I means B[i][j] = 2*A[i][j] - 1 for i≠j and B[i][i] = 0.
So B is the SIGNED adjacency matrix: B[i][j] = +1 if i→j, -1 if j→i.

Pf(B) = sum over perfect matchings M of prod sign(M) * prod B[i_k][j_k]
       = sum over perfect matchings of sgn * prod (±1)
       = ±(some integer)

For odd n, det(B) = 0 (skew-symmetric odd matrix).

QUESTION: Is Pf(B) ≡ H mod something?

Also: the HAFNIAN of B: Haf(B) = sum over perfect matchings (without sign).
For skew-symmetric ±1 matrices, Haf counts matchings weighted by ±1.

Author: opus-2026-03-07-S44
"""
import numpy as np
from itertools import permutations, combinations
import math
import random

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def count_t3(A, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            t3 += 1
    return t3

def pfaffian(M):
    """Compute Pfaffian of skew-symmetric matrix M of even size."""
    n = len(M)
    if n == 0:
        return 1
    if n == 2:
        return M[0][1]
    if n % 2 == 1:
        return 0

    # Expansion along first row
    result = 0
    for j in range(1, n):
        if M[0][j] == 0:
            continue
        # Remove rows/cols 0 and j
        sub = []
        for r in range(1, n):
            if r == j:
                continue
            row = []
            for c in range(1, n):
                if c == j:
                    continue
                row.append(M[r][c])
            sub.append(row)
        sign = (-1) ** (j - 1)
        result += sign * M[0][j] * pfaffian(sub)
    return result

# ============================================================
# PFAFFIAN vs H at n=4, 6
# ============================================================
print("=" * 60)
print("PFAFFIAN OF SIGNED ADJACENCY vs H(T)")
print("=" * 60)

for n in [4, 6]:
    m = n*(n-1)//2
    print(f"\n--- n={n} ---")

    seen = set()
    data = []

    for bits in range(1 << m):
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        H = F[n-1]
        t3 = count_t3(A, n)

        # Build signed adjacency B
        B = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    B[i][j] = 1 if A[i][j] else -1

        pf = pfaffian(B)
        data.append((H, pf, t3))

    data.sort()
    print(f"  {'H':>3} {'Pf':>4} {'t3':>3} {'H mod Pf':>8} {'H-Pf':>6}")
    for H, pf, t3 in data:
        h_mod_pf = H % abs(pf) if pf != 0 else 'N/A'
        print(f"  {H:3d} {pf:4d} {t3:3d} {str(h_mod_pf):>8} {H-abs(pf):6d}")

# ============================================================
# PFAFFIAN vs H at n=8 (sampling)
# ============================================================
print("\n" + "=" * 60)
print("PFAFFIAN vs H at n=8 (sampling)")
print("=" * 60)

n = 8
m = n*(n-1)//2
random.seed(42)
seen = set()
data = []

for trial in range(100):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]

    B = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                B[i][j] = 1 if A[i][j] else -1

    pf = pfaffian(B)
    data.append((H, pf))

data.sort()
print(f"  H range: [{data[0][0]}, {data[-1][0]}]")
print(f"  Pf range: [{min(d[1] for d in data)}, {max(d[1] for d in data)}]")

# Check: is Pf always odd?
all_odd_pf = all(abs(d[1]) % 2 == 1 for d in data)
print(f"  Pf always odd: {all_odd_pf}")

# Check: Pf mod small numbers
for p in [2, 3, 4, 5, 7]:
    vals = sorted(set(d[1] % p for d in data))
    print(f"  Pf mod {p}: {vals}")

# Show examples
print(f"\n  Examples:")
for H, pf in data[:15]:
    print(f"    H={H:5d}, Pf={pf:6d}, H mod |Pf|={H % abs(pf) if pf != 0 else 'N/A'}")

# ============================================================
# NOVEL: H mod Pf — is it always 0?
# ============================================================
print("\n" + "=" * 60)
print("DOES Pf(B) DIVIDE H(T)?")
print("=" * 60)

for n in [4, 6]:
    m = n*(n-1)//2
    print(f"\n  n={n}:")
    seen = set()
    all_divide = True

    for bits in range(1 << m):
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        H = F[n-1]
        B = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    B[i][j] = 1 if A[i][j] else -1
        pf = pfaffian(B)

        if pf != 0 and H % abs(pf) != 0:
            print(f"    FAILS: H={H}, Pf={pf}, H mod |Pf| = {H % abs(pf)}")
            all_divide = False

    if all_divide:
        print(f"    Pf(B) ALWAYS divides H(T)!")
    else:
        print(f"    Divisibility FAILS for some tournaments")

# ============================================================
# NOVEL: Pf(B)^2 = det(B). What is det(B) in terms of t3?
# ============================================================
print("\n" + "=" * 60)
print("det(B) vs t3")
print("=" * 60)

for n in [4, 6]:
    m = n*(n-1)//2
    print(f"\n  n={n}:")
    seen = set()
    t3_det = {}

    for bits in range(1 << m):
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        H = F[n-1]
        t3 = count_t3(A, n)
        B = np.array([[1 if A[i][j] else (-1 if A[j][i] else 0) for j in range(n)] for i in range(n)], dtype=float)
        det_B = round(np.linalg.det(B))

        if t3 not in t3_det:
            t3_det[t3] = set()
        t3_det[t3].add((H, det_B))

    for t3 in sorted(t3_det.keys()):
        vals = sorted(t3_det[t3])
        print(f"    t3={t3}: {vals}")
