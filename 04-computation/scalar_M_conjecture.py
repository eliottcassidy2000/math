#!/usr/bin/env python3
"""
CONJECTURE: M = (H/n)*I for vertex-transitive tournaments at odd n.

Test: is this true for ALL vertex-transitive tournaments, or just circulant ones?
What about non-vertex-transitive regular tournaments?

kind-pasteur-2026-03-06-S25
"""

from itertools import permutations
import numpy as np
import random

def compute_M_entry(T, n, a, b):
    """Compute M[a,b] for numeric tournament T."""
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

def compute_full_M(T, n):
    M = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            M[a,b] = compute_M_entry(T, n, a, b)
    return M

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(len(perm)-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def is_regular(T, n):
    for i in range(n):
        out_deg = sum(T.get((i,j), 0) for j in range(n) if j != i)
        if out_deg != (n-1)//2:
            return False
    return True

# ============================================================
# Test ALL tournaments at n=5 (enumerate)
# ============================================================
print("=" * 70)
print("n=5: ALL tournaments - check M structure")
print("=" * 70)

# Classify by score sequence
from collections import defaultdict
score_classes = defaultdict(list)

random.seed(42)
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
seen = set()

# Sample random tournaments
for trial in range(50000):
    T = {}
    for (i,j) in pairs:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    key = tuple(T.get((i,j), 0) for (i,j) in pairs)
    if key in seen:
        continue
    seen.add(key)

    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
    H = count_H(T, n)

    # Only compute M for a few per class
    if len(score_classes[scores]) < 5:
        M = compute_full_M(T, n)
        diag_vals = sorted(set(int(round(M[i,i])) for i in range(n)))
        off_diag_vals = sorted(set(int(round(M[i,j])) for i in range(n) for j in range(n) if i != j))
        is_scalar = np.allclose(M, M[0,0] * np.eye(n))
        score_classes[scores].append((H, is_scalar, diag_vals, off_diag_vals))

print(f"  Found {len(seen)} distinct labeled tournaments")
for scores in sorted(score_classes.keys()):
    entries = score_classes[scores]
    print(f"\n  Score {scores}:")
    for H, is_sc, dv, ov in entries:
        if is_sc:
            print(f"    H={H}: M = {dv[0]}*I (SCALAR)")
        else:
            print(f"    H={H}: diag={dv}, off-diag={ov}")


# ============================================================
# n=7: Check non-circulant regular tournaments
# ============================================================
print("\n" + "=" * 70)
print("n=7: Random regular tournaments")
print("=" * 70)

n = 7
pairs7 = [(i,j) for i in range(n) for j in range(i+1, n)]
seen7 = set()
results7 = []

for trial in range(200000):
    T = {}
    for (i,j) in pairs7:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    if not is_regular(T, n):
        continue

    key = tuple(T.get((i,j), 0) for (i,j) in pairs7)
    if key in seen7:
        continue
    seen7.add(key)

    H = count_H(T, n)

    # Compute just M[0,0], M[0,1], M[1,0] to check
    M00 = compute_M_entry(T, n, 0, 0)
    M01 = compute_M_entry(T, n, 0, 1)
    M10 = compute_M_entry(T, n, 1, 0)
    M11 = compute_M_entry(T, n, 1, 1)

    is_scalar_hint = (M01 == 0 and M10 == 0 and M00 == M11)
    results7.append((H, M00, M01, M10, M11, is_scalar_hint))

    if len(results7) >= 30:
        break

print(f"  Found {len(results7)} distinct regular n=7 tournaments")
for H, M00, M01, M10, M11, scalar in results7:
    if scalar:
        print(f"  H={H}: M[0,0]={M00}, M[0,1]={M01} -> scalar M={M00}*I")
    else:
        print(f"  H={H}: M[0,0]={M00}, M[0,1]={M01}, M[1,0]={M10}, M[1,1]={M11} -> NOT scalar")


# ============================================================
# Deeper check: does the non-scalar structure relate to Aut(T)?
# ============================================================
print("\n" + "=" * 70)
print("n=7: Non-scalar M structure analysis")
print("=" * 70)

# For a non-vertex-transitive tournament, what does M look like?
# Pick a non-regular tournament
T = {}
# Transitive tournament: 0->1->2->3->4->5->6
for i in range(7):
    for j in range(7):
        if i != j:
            T[(i,j)] = 1 if i < j else 0

n = 7
H = count_H(T, n)
M = compute_full_M(T, n)
print(f"\n  Transitive T_7: H = {H}")
print(f"    M =")
for row in M:
    print(f"      {[int(x) for x in row]}")
print(f"    trace = {int(np.trace(M))}")
evals = sorted(np.linalg.eigvalsh(M))[::-1]
print(f"    eigenvalues: {[round(e, 2) for e in evals]}")

# Almost-transitive: flip one arc
T2 = dict(T)
T2[(0,6)] = 0; T2[(6,0)] = 1  # 6 now beats 0
H2 = count_H(T2, n)
M2 = compute_full_M(T2, n)
print(f"\n  Almost-transitive (flip 0->6): H = {H2}")
print(f"    M =")
for row in M2:
    print(f"      {[int(x) for x in row]}")
print(f"    trace = {int(np.trace(M2))}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
