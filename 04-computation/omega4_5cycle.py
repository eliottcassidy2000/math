#!/usr/bin/env python3
"""
omega4_5cycle.py — opus-2026-03-13-S70

The two (1,2,2,2,3) types with t3=4, c4=12 have different omega4_local (0 vs 1).
What invariant distinguishes them? Check 5-cycle count, Omega_2, Omega_3 locally.

Also: express omega4_local as a polynomial in (t3, t5) or similar.
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict

def adj_matrix_5(bits):
    A = np.zeros((5,5), dtype=int)
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def omega4_local(A):
    count = 0
    for perm in permutations(range(5)):
        v0, v1, v2, v3, v4 = perm
        if (A[v0][v1] and A[v1][v2] and A[v2][v3] and A[v3][v4] and
            A[v0][v2] and A[v1][v3] and A[v2][v4]):
            count += 1
    return count

def count_3cycles(A):
    n = A.shape[0]
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                    count += 1
    return count

def count_5cycles(A):
    """Count oriented 5-cycles a→b→c→d→e→a on 5 vertices."""
    count = 0
    for perm in permutations(range(5)):
        a, b, c, d, e = perm
        if A[a][b] and A[b][c] and A[c][d] and A[d][e] and A[e][a]:
            count += 1
    return count

def count_oriented_4cycles(A):
    n = A.shape[0]
    count = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                for d in range(n):
                    if d == a or d == b or d == c: continue
                    if A[c][d] and A[d][a]:
                        count += 1
    return count

def omega2_local(A):
    """Count regular 2-paths = transitive triples."""
    n = A.shape[0]
    count = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    count += 1
    return count

def omega3_local(A):
    """Count regular 3-paths."""
    n = A.shape[0]
    count = 0
    for v0 in range(n):
        for v1 in range(n):
            if v1 == v0 or not A[v0][v1]: continue
            for v2 in range(n):
                if v2 == v0 or v2 == v1: continue
                if not A[v1][v2] or not A[v0][v2]: continue
                for v3 in range(n):
                    if v3 == v0 or v3 == v1 or v3 == v2: continue
                    if A[v2][v3] and A[v1][v3]:
                        count += 1
    return count

def canon_form(A):
    n = A.shape[0]
    best = None
    for perm in permutations(range(n)):
        enc = []
        for i in range(n):
            for j in range(i+1, n):
                enc.append(A[perm[i]][perm[j]])
        enc = tuple(enc)
        if best is None or enc < best:
            best = enc
    return best

# Enumerate all 12 types at n=5
print("="*70)
print("ALL 12 TOURNAMENT TYPES ON 5 VERTICES")
print("="*70)

types = {}
for bits in range(2**10):
    A = adj_matrix_5(bits)
    cf = canon_form(A)
    if cf not in types:
        types[cf] = bits  # representative

print(f"  Found {len(types)} types\n")

print(f"  {'Type':>4s} {'scores':20s} {'t3':>3s} {'c4':>4s} {'c5':>4s} {'O2':>4s} {'O3':>4s} {'O4':>4s} {'C4-C5':>5s}")
print(f"  {'----':>4s} {'---':20s} {'---':>3s} {'----':>4s} {'----':>4s} {'----':>4s} {'----':>4s} {'----':>4s} {'-----':>5s}")

results = []
for i, (cf, bits) in enumerate(sorted(types.items())):
    A = adj_matrix_5(bits)
    scores = tuple(sorted(sum(A[v][w] for w in range(5) if w != v) for v in range(5)))
    t3 = count_3cycles(A)
    c4 = count_oriented_4cycles(A)
    c5 = count_5cycles(A)
    o2 = omega2_local(A)
    o3 = omega3_local(A)
    o4 = omega4_local(A)

    results.append((scores, t3, c4, c5, o2, o3, o4))
    print(f"  {i+1:4d} {str(scores):20s} {t3:3d} {c4:4d} {c5:4d} {o2:4d} {o3:4d} {o4:4d} {c4-c5:5d}")

# Check: is omega4_local = f(t3, c5)?
print(f"\n  Is omega4_local = f(t3, c5)?")
by_t3_c5 = defaultdict(set)
for scores, t3, c4, c5, o2, o3, o4 in results:
    by_t3_c5[(t3, c5)].add(o4)
all_det = True
for (t3, c5), o4s in sorted(by_t3_c5.items()):
    det = "YES" if len(o4s) == 1 else "NO"
    print(f"    t3={t3}, c5={c5}: omega4_local ∈ {sorted(o4s)} [{det}]")
    if len(o4s) > 1:
        all_det = False
if all_det:
    print("    => YES, omega4_local is determined by (t3, c5)")

# Check: is omega4_local = f(c4, c5)?
print(f"\n  Is omega4_local = f(c4, c5)?")
by_c4_c5 = defaultdict(set)
for scores, t3, c4, c5, o2, o3, o4 in results:
    by_c4_c5[(c4, c5)].add(o4)
for (c4, c5), o4s in sorted(by_c4_c5.items()):
    det = "YES" if len(o4s) == 1 else "NO"
    print(f"    c4={c4}, c5={c5}: omega4_local ∈ {sorted(o4s)} [{det}]")

# Check: is omega4_local = f(o3)?
print(f"\n  Is omega4_local = f(o3)?")
by_o3 = defaultdict(set)
for scores, t3, c4, c5, o2, o3, o4 in results:
    by_o3[o3].add(o4)
for o3_val, o4s in sorted(by_o3.items()):
    det = "YES" if len(o4s) == 1 else "NO"
    print(f"    Omega_3={o3_val}: omega4_local ∈ {sorted(o4s)} [{det}]")

# Try to fit omega4_local as polynomial in (t3, c5)
print(f"\n  Fitting omega4_local = a*t3 + b*c5 + c*t3^2 + d*c5^2 + e*t3*c5 + f:")
from numpy.linalg import lstsq
X = []
y = []
for scores, t3, c4, c5, o2, o3, o4 in results:
    X.append([t3, c5, t3**2, c5**2, t3*c5, 1])
    y.append(o4)
X = np.array(X, dtype=float)
y = np.array(y, dtype=float)
coeffs, residuals, _, _ = lstsq(X, y, rcond=None)
print(f"    coeffs = {coeffs}")
print(f"    residuals = {residuals}")
predicted = X @ coeffs
for i, (scores, t3, c4, c5, o2, o3, o4) in enumerate(results):
    print(f"    scores={str(scores):20s} actual={o4}, predicted={predicted[i]:.4f}")

# Check: omega4_local vs Omega_3 + t3 + c5 combo
print(f"\n  Trying omega4_local = a*Omega_3 + b*t3 + c*c5 + d:")
X2 = []
for scores, t3, c4, c5, o2, o3, o4 in results:
    X2.append([o3, t3, c5, 1])
X2 = np.array(X2, dtype=float)
coeffs2, res2, _, _ = lstsq(X2, y, rcond=None)
print(f"    coeffs = {coeffs2}")
predicted2 = X2 @ coeffs2
for i, (scores, t3, c4, c5, o2, o3, o4) in enumerate(results):
    print(f"    actual={o4}, predicted={predicted2[i]:.4f}")

# The key insight: the FULL homology profile (O2, O3, O4) is all determined
# by isomorphism type. What's the relationship?
print(f"\n  Full Omega profile by type:")
print(f"  {'scores':20s} {'t3':>3s} {'c5':>4s} | {'O0':>3s} {'O1':>3s} {'O2':>3s} {'O3':>3s} {'O4':>3s} | {'chi':>4s}")
for scores, t3, c4, c5, o2, o3, o4 in results:
    # O0 = n = 5, O1 = n*(n-1)/2 = 20 always? No, O1 = number of 1-paths = edges = C(n,2)
    o0 = 5
    o1 = 10  # C(5,2) edges, all directed, all are 1-paths
    chi = o0 - o1 + o2 - o3 + o4
    print(f"  {str(scores):20s} {t3:3d} {c5:4d} | {o0:3d} {o1:3d} {o2:3d} {o3:3d} {o4:3d} | {chi:4d}")

print("\nDONE.")
