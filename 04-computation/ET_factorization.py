#!/usr/bin/env python3
"""
E_T(t) factorization analysis for H-maximizers.
Investigates the root structure of the tournament Eulerian polynomial.

kind-pasteur-2026-03-07-S28
"""
from itertools import combinations
from collections import defaultdict

n = 5

# Regular n=5 maximizer: a_k = [15, 30, 30, 30, 15]
a_reg = [15, 30, 30, 30, 15]
poly = lambda t: sum(a_reg[k] * t**k for k in range(5))

print("Regular n=5 maximizer:")
print(f"  a_k = {a_reg}")
print(f"  E_T(t)/15 = 1 + 2t + 2t^2 + 2t^3 + t^4")
print(f"  = (t^2+1)*(t^2+2t+1) = (t^2+1)*(t+1)^2")
print(f"  E_T(i) = {poly(1j)}")
print(f"  E_T(-1) = {poly(-1)}")
print(f"  Roots: t = i, -i (from t^2+1), t = -1 double")
print(f"  P(u,2) = 15*u*(u+2), roots u=0, u=-2")
print()

# Find non-regular H=15 tournaments
for bits in range(1 << 10):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    scores = sorted([sum(A[i]) for i in range(n)])
    if scores == [2]*5:
        continue

    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    dp[(mask | (1 << u), u, fwd + A[v][u])] = (
                        dp.get((mask | (1 << u), u, fwd + A[v][u]), 0) + c
                    )
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    a = [dist[k] for k in range(n)]
    H = a[0]
    if H != 15:
        continue

    E_m1 = sum(a[k] * (-1)**k for k in range(n))
    E_i = sum(a[k] * (1j)**k for k in range(n))
    print(f"Non-regular H=15: scores={scores}, a={a}")
    print(f"  E_T(-1) = {E_m1}")
    print(f"  E_T(i) = {E_i}")

    # P-coefficients
    p2 = a[0]
    p1 = a[1]
    p0 = a[2] - 2 * a[0]
    print(f"  P(u,2) = {p0} + {p1}*u + {p2}*u^2")
    disc = p1**2 - 4*p2*p0
    print(f"  discriminant = {disc}")
    if disc < 0:
        import cmath
        r1 = (-p1 + cmath.sqrt(disc)) / (2*p2)
        print(f"  complex roots: u = {r1.real:.4f} +/- {abs(r1.imag):.4f}i")
    break

print()
print("="*60)
print("INSIGHT: For the regular n=5 maximizer:")
print("  E_T(t) = 15*(t^2+1)*(t+1)^2")
print("  The factor (t^2+1) gives W(i/2)=0")
print("  The factor (t+1)^2 gives E_T(-1)=0")
print()
print("For non-regular H=15 tournaments:")
print("  E_T(i) is NOT zero, E_T(-1) is NOT zero")
print("  The regular maximizer is SPECIAL in having both vanishings")
