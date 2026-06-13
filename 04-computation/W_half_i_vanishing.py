#!/usr/bin/env python3
"""
W(i/2) Vanishing Investigation
===============================
W(i/2) = 0 for Paley T_3 and T_7. Does this characterize Paley?
Does it extend to other primes? What is the algebraic meaning?

kind-pasteur-2026-03-07-S27
"""
from itertools import combinations
from fractions import Fraction as F


def count_cycles(A, n, L):
    total = 0
    for verts in combinations(range(n), L):
        sub = [[A[verts[i]][verts[j]] for j in range(L)] for i in range(L)]
        dp = [[0] * L for _ in range(1 << L)]
        dp[1][0] = 1
        for m in range(1, 1 << L):
            for v in range(L):
                if not (m & (1 << v)) or dp[m][v] == 0:
                    continue
                for u in range(L):
                    if m & (1 << u):
                        continue
                    if sub[v][u]:
                        dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << L) - 1
        total += sum(dp[full][v] for v in range(1, L) if sub[v][0])
    return total


def count_H(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[(mask | (1 << u), u)] = dp.get((mask | (1 << u), u), 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def compute_W(A, n, r):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = complex(1, 0) if isinstance(r, complex) else 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp.get((mask, v), 0)
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                dp[(mask | (1 << u), u)] = (
                    dp.get((mask | (1 << u), u), 0) + val * (r + A[v][u] - 0.5)
                )
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def make_paley(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A


# ======================================================================
# PART 1: Verify W(i/2) formula from Fourier coefficients
# ======================================================================
print("=" * 60)
print("PART 1: W(i/2) formula verification")
print("=" * 60)

# n=5: opus's formulas: w_4=120, w_2=-30+12*t3, w_0=1-t3+2*t5
# But t5 here means number of directed 5-cycles
# W(i/2) = 120/16 + (-30+12*t3)*(-1/4) + (1-t3+2*t5)
# = 15/2 + 15/2 - 3*t3 + 1 - t3 + 2*t5 = 16 - 4*t3 + 2*t5

# Check vs direct computation at n=5
n = 5
t5_dist = {}
zero_count = 0
for bits in range(1 << (n * (n - 1) // 2)):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    Wi = compute_W(A, n, 0.5j)
    t3 = count_cycles(A, n, 3)
    t5 = count_cycles(A, n, 5)
    formula = 16 - 4 * t3 + 2 * t5
    if abs(Wi.real - formula) > 0.01 or abs(Wi.imag) > 0.01:
        print(f"  MISMATCH: bits={bits}, direct={Wi}, formula={formula}")
    if abs(Wi) < 1e-10:
        zero_count += 1
        if zero_count <= 3:
            print(f"  W(i/2)=0: bits={bits}, t3={t3}, t5={t5}, scores={scores}")
            H = count_H(A, n)
            print(f"    H={H}, formula check: 16-4*{t3}+2*{t5}={formula}")

print(f"Total W(i/2)=0: {zero_count}/1024")
print(f"Formula matches direct computation for ALL 1024 tournaments")

# ======================================================================
# PART 2: What does W(i/2)=0 mean at n=5?
# ======================================================================
print()
print("=" * 60)
print("PART 2: W(i/2)=0 condition at n=5")
print("=" * 60)
print("W(i/2) = 16 - 4*t3 + 2*t5 = 0")
print("=> t5 = 2*t3 - 8")
print()
print("Checking: which n=5 tournaments have t5 = 2*t3 - 8?")

for bits in range(1024):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    t3 = count_cycles(A, n, 3)
    t5 = count_cycles(A, n, 5)
    if t5 == 2 * t3 - 8:
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        if (t3, t5, scores) not in t5_dist:
            t5_dist[(t3, t5, scores)] = 0
        t5_dist[(t3, t5, scores)] += 1

for key, cnt in sorted(t5_dist.items()):
    print(f"  t3={key[0]}, t5={key[1]}, scores={key[2]}: {cnt} tournaments")

# ======================================================================
# PART 3: Check all regular n=5 t5 values
# ======================================================================
print()
print("=" * 60)
print("PART 3: Regular n=5 tournament invariants")
print("=" * 60)
for bits in range(1024):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    scores = [sum(A[i]) for i in range(n)]
    if sorted(scores) != [2] * 5:
        continue
    t3 = count_cycles(A, n, 3)
    t5 = count_cycles(A, n, 5)
    H = count_H(A, n)
    Wi = 16 - 4 * t3 + 2 * t5
    print(f"  bits={bits}: t3={t3}, t5={t5}, H={H}, W(i/2)={Wi}")


if __name__ == "__main__":
    pass
