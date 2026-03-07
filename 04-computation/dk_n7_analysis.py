#!/usr/bin/env python3
"""
dk_n7_analysis.py — Analyze D_k dependencies at n=7.

At n=5: D_2 = 150 + 12*t3, D_3 = 30 + 12*t3 (both exactly linear in t3).
What about n=7?

Author: opus-2026-03-07-S43b
"""
from itertools import permutations, combinations
from collections import defaultdict
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

def compute_all(A, n):
    """Compute H, S, D_k, t3, t5, t7."""
    D = [0] * n
    S = 0

    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        for k in range(n):
            D[k] += math.comb(fwd, k)
        prod_b = 1
        for i in range(n-1):
            prod_b *= (2*A[P[i]][P[i+1]] - 1)
        S += prod_b

    H = D[n-1]

    # Count directed 3-cycles
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            t3 += 1

    # Count directed 5-cycles
    t5 = 0
    for quint in combinations(range(n), 5):
        for perm in permutations(quint):
            if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                t5 += 1
    t5 //= 10  # 5 rotations * 2 directions

    # Count directed 7-cycles (only at n=7)
    t7 = 0
    if n == 7:
        for perm in permutations(range(n)):
            if all(A[perm[i]][perm[(i+1)%7]] for i in range(7)):
                t7 += 1
        t7 //= 14  # 7 rotations * 2 directions

    return H, S, D, t3, t5, t7

n = 7
m = n*(n-1)//2
random.seed(123)

print(f"=== n={n}: D_k analysis ===")
print(f"{'t3':>3} {'t5':>3} {'t7':>4} {'D2':>6} {'D3':>6} {'D4':>6} {'D5':>6} {'H':>5} {'S':>7}")

results = []
for trial in range(50):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    H, S, D, t3, t5, t7 = compute_all(A, n)
    results.append((t3, t5, t7, D, H, S))
    if trial < 15:
        print(f"{t3:3d} {t5:3d} {t7:4d} {D[2]:6d} {D[3]:6d} {D[4]:6d} {D[5]:6d} {H:5d} {S:7d}")

# Check linearity of D_2, D_3 in t3
print("\n=== D_2 vs t3 ===")
# D_2 = a + b*t3 + c*t5 + d*t7?
import numpy as np

t3s = np.array([r[0] for r in results], dtype=float)
t5s = np.array([r[1] for r in results], dtype=float)
t7s = np.array([r[2] for r in results], dtype=float)
D2s = np.array([r[3][2] for r in results], dtype=float)
D3s = np.array([r[3][3] for r in results], dtype=float)
D4s = np.array([r[3][4] for r in results], dtype=float)
D5s = np.array([r[3][5] for r in results], dtype=float)
Hs = np.array([r[4] for r in results], dtype=float)
Ss = np.array([r[5] for r in results], dtype=float)

# D_2 = a + b*t3
X1 = np.column_stack([np.ones_like(t3s), t3s])
for name, Y in [("D_2", D2s), ("D_3", D3s), ("D_4", D4s), ("D_5", D5s)]:
    c, res, _, _ = np.linalg.lstsq(X1, Y, rcond=None)
    pred = X1 @ c
    maxerr = max(abs(Y - pred))
    print(f"  {name} = {c[0]:.1f} + {c[1]:.1f}*t3, max_err={maxerr:.1f}")

# D_2 = a + b*t3 + c*t5
X2 = np.column_stack([np.ones_like(t3s), t3s, t5s])
for name, Y in [("D_2", D2s), ("D_3", D3s), ("D_4", D4s), ("D_5", D5s)]:
    c, res, _, _ = np.linalg.lstsq(X2, Y, rcond=None)
    pred = X2 @ c
    maxerr = max(abs(Y - pred))
    print(f"  {name} = {c[0]:.1f} + {c[1]:.1f}*t3 + {c[2]:.1f}*t5, max_err={maxerr:.1f}")

# D_2 = a + b*t3 + c*t5 + d*t7
X3 = np.column_stack([np.ones_like(t3s), t3s, t5s, t7s])
for name, Y in [("D_2", D2s), ("D_3", D3s), ("D_4", D4s), ("D_5", D5s)]:
    c, res, _, _ = np.linalg.lstsq(X3, Y, rcond=None)
    pred = X3 @ c
    maxerr = max(abs(Y - pred))
    print(f"  {name} = {c[0]:.1f} + {c[1]:.1f}*t3 + {c[2]:.1f}*t5 + {c[3]:.1f}*t7, max_err={maxerr:.1f}")

# Check S/64 formula
print("\n=== S/64 regression ===")
c064 = Ss / 64.0
for name, X in [("t3", X1), ("t3,t5", X2), ("t3,t5,t7", X3)]:
    c, res, _, _ = np.linalg.lstsq(X, c064, rcond=None)
    pred = X @ c
    maxerr = max(abs(c064 - pred))
    print(f"  S/64 = {' + '.join(f'{v:.4f}*{n}' for n, v in zip(['1','t3','t5','t7'][:len(c)], c))}, max_err={maxerr:.4f}")

# Check: S/64 = H - 3*t3 - 5*t5 - 7*t7? (conjecture from pattern)
print("\n=== Testing S/64 = H - sum (2k+1)*t_{2k+1} patterns ===")
for a, b, c_coef in [(3, 5, 7), (3, 0, 0), (0, 5, 0), (3, 5, 0)]:
    pred = Hs - a*t3s - b*t5s - c_coef*t7s
    maxerr = max(abs(c064 - pred))
    print(f"  S/64 = H - {a}*t3 - {b}*t5 - {c_coef}*t7? max_err={maxerr:.4f}")

# What about S/64 = I(Omega, 2) - 3*(I(Omega,1)-1)?
# I(Omega,1) = 1 + alpha1, so I(2) - 3*alpha1 = H - 3*alpha1
# alpha1 = t3 + t5 + t7
alpha1s = t3s + t5s + t7s
pred_ia = Hs - 3*alpha1s
maxerr = max(abs(c064 - pred_ia))
print(f"\n  S/64 = H - 3*alpha1? max_err={maxerr:.4f}")

# Try: S/64 = I(Omega, -2)/4^? or some evaluation
# I(Omega, -2) = 1 - 2*alpha1 + 4*alpha2 - 8*alpha3
# This involves alpha2, alpha3 which we haven't computed. Let's try.

# Actually let's check: does D_2 = a + b*t3 exactly?
# From the skeleton doc: D_2 = 16800 + 240*t3 at n=7
print("\n=== Checking D_2 = 16800 + 240*t3 at n=7 ===")
exact = True
for r in results:
    t3, t5, t7, D, H, S = r
    pred = 16800 + 240*t3
    if D[2] != pred:
        print(f"  MISMATCH: t3={t3}, D_2={D[2]}, pred={pred}, diff={D[2]-pred}")
        exact = False
        break
if exact:
    print("  CONFIRMED for all 50 samples!")

# Similarly: D_3 = 8400 + 480*t3?
print("\n=== Checking D_3 = 8400 + 480*t3 at n=7 ===")
exact = True
for r in results:
    t3, t5, t7, D, H, S = r
    pred = 8400 + 480*t3
    if D[3] != pred:
        print(f"  MISMATCH: t3={t3}, D_3={D[3]}, pred={pred}, diff={D[3]-pred}")
        exact = False
        break
if exact:
    print("  CONFIRMED for all 50 samples!")
