#!/usr/bin/env python3
"""
Determinant of the transfer matrix M for various tournaments.

For VT tournaments at odd n: M = (H/n)*I, so det(M) = (H/n)^n.
For transitive tournaments: |det(M)| = F(n+1) (Fibonacci).

What about general tournaments?

kind-pasteur-2026-03-06-S25 (continuation)
"""

from itertools import permutations
import numpy as np
import random

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

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count


# ============================================================
# det(M) for various tournaments
# ============================================================
print("=" * 70)
print("det(M) for various tournaments")
print("=" * 70)

# Collect (H, det(M)) pairs
results = {}

# n=3
n = 3
print(f"\n  n={n}:")

# Transitive
T = {}
for i in range(n):
    for j in range(n):
        if i != j:
            T[(i,j)] = 1 if i < j else 0
M = np.zeros((n,n));
for a in range(n):
    for b in range(n):
        M[a,b] = compute_M_entry(T, n, a, b)
H = count_H(T, n)
d = int(round(np.linalg.det(M)))
print(f"    Transitive: H={H}, det(M)={d}")
results.setdefault(n, []).append((H, d, "trans"))

# Cyclic
T = make_circulant(n, {1})
M = np.zeros((n,n))
for a in range(n):
    for b in range(n):
        M[a,b] = compute_M_entry(T, n, a, b)
H = count_H(T, n)
d = int(round(np.linalg.det(M)))
print(f"    Cyclic: H={H}, det(M)={d}")
results.setdefault(n, []).append((H, d, "cyclic"))


# n=5
n = 5
print(f"\n  n={n}:")

# Transitive
T = {}
for i in range(n):
    for j in range(n):
        if i != j:
            T[(i,j)] = 1 if i < j else 0
M = np.zeros((n,n))
for a in range(n):
    for b in range(n):
        M[a,b] = compute_M_entry(T, n, a, b)
H = count_H(T, n)
d = int(round(np.linalg.det(M)))
print(f"    Transitive: H={H}, det(M)={d}")
results.setdefault(n, []).append((H, d, "trans"))

# Paley
T = make_circulant(n, {1,2})
M = np.zeros((n,n))
for a in range(n):
    for b in range(n):
        M[a,b] = compute_M_entry(T, n, a, b)
H = count_H(T, n)
d = int(round(np.linalg.det(M)))
print(f"    Paley: H={H}, det(M)={d}, (H/n)^n={(H//n)**n}")
results.setdefault(n, []).append((H, d, "paley"))

# circ {1,3}
T = make_circulant(n, {1,3})
M = np.zeros((n,n))
for a in range(n):
    for b in range(n):
        M[a,b] = compute_M_entry(T, n, a, b)
H = count_H(T, n)
d = int(round(np.linalg.det(M)))
print(f"    circ {{1,3}}: H={H}, det(M)={d}")
results.setdefault(n, []).append((H, d, "circ13"))

# Some random non-regular
random.seed(42)
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
seen = set()
count = 0
for trial in range(10000):
    T = {}
    for (i,j) in pairs:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    key = tuple(T.get((i,j),0) for (i,j) in pairs)
    if key in seen: continue
    seen.add(key)

    H = count_H(T, n)
    if H in [1, 15]: continue  # Skip trans and regular

    M = np.zeros((n,n))
    for a in range(n):
        for b in range(n):
            M[a,b] = compute_M_entry(T, n, a, b)
    d = int(round(np.linalg.det(M)))

    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
    print(f"    scores={scores}: H={H}, det(M)={d}")
    results.setdefault(n, []).append((H, d, f"scores={scores}"))
    count += 1
    if count >= 8: break


# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY: det(M) patterns")
print("=" * 70)

for n in sorted(results.keys()):
    print(f"\n  n={n}:")
    for H, d, label in results[n]:
        # Check if det(M) = (H/n)^n for scalar M
        if n % 2 == 1:
            scalar_det = round((H/n)**n)
            is_scalar = abs(d - scalar_det) < 0.5
        else:
            is_scalar = False
        print(f"    {label}: H={H}, det(M)={d}" + (f" = (H/n)^n = {scalar_det}" if is_scalar else ""))


# ============================================================
# det(M) vs H: is det(M) always odd?
# ============================================================
print("\n" + "=" * 70)
print("Is det(M) always odd?")
print("=" * 70)

random.seed(123)
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
seen = set()
odd_count = 0
total = 0

for trial in range(50000):
    T = {}
    for (i,j) in pairs:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    key = tuple(T.get((i,j),0) for (i,j) in pairs)
    if key in seen: continue
    seen.add(key)

    H = count_H(T, n)
    M = np.zeros((n,n))
    for a in range(n):
        for b in range(n):
            M[a,b] = compute_M_entry(T, n, a, b)
    d = int(round(np.linalg.det(M)))

    total += 1
    if d % 2 == 1:
        odd_count += 1
    else:
        print(f"  EVEN det: H={H}, det(M)={d}")
        break

    if total >= 200:
        break

print(f"  Tested {total} tournaments at n={n}: {odd_count} odd, {total-odd_count} even")
if odd_count == total:
    print("  ALL determinants are odd! (consistent with Redei)")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
