#!/usr/bin/env python3
"""
Verify and explore: w_0 = H - 3*t_3 at n=5.
Then check if w_0 = H - C(n)*t_3 at n=7.

The identity w_0 = H - 3*t_3 means:
  W(0) = H(T) - 3*t_3(T)

Since W(0) = sum_P prod(s_e) = sum_P prod(A[e] - 1/2)
= (1/2)^{n-1} sum_P (-1)^{b(P)}

where b(P) = number of backward arcs in path P.

So: sum_P (-1)^{b(P)} = 2^{n-1} * (H - 3*t_3) at n=5.

This signed count has deep structure!

kind-pasteur-2026-03-06-S25g
"""

from itertools import permutations, combinations
from math import factorial, comb

# Verify at n=5 exhaustively
n = 5
print("=" * 70)
print(f"w_0 = H - C*t_3 at n={n}")
print("=" * 70)

for bits in range(1 << 10):
    A = [[0]*5 for _ in range(5)]
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    H = sum(1 for p in permutations(range(5))
            if all(A[p[i]][p[i+1]] == 1 for i in range(4)))

    t3 = sum(1 for v in combinations(range(5), 3)
             for p in permutations(v)
             if all(A[p[i]][p[(i+1)%3]] == 1 for i in range(3))) // 3

    # W(0)
    w0 = 0.0
    for p in permutations(range(5)):
        prod = 1.0
        for i in range(4):
            prod *= (A[p[i]][p[i+1]] - 0.5)
        w0 += prod

    predicted = H - 3*t3
    if abs(w0 - predicted) > 0.01:
        scores = tuple(sorted(sum(A[i]) for i in range(5)))
        print(f"FAIL: bits={bits}, H={H}, t3={t3}, w0={w0:.4f}, predicted={predicted}")
        break
else:
    print(f"ALL 1024 tournaments: w_0 = H - 3*t_3  [VERIFIED]")
    print(f"Coefficient: C(5) = 3")

# Now check n=7 (sampling)
n = 7
print(f"\n{'='*70}")
print(f"w_0 = H - C*t_3 at n={n}?")
print(f"{'='*70}")

import random
random.seed(42)

data = []
for trial in range(30):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = sum(1 for v in combinations(range(n), 3)
             for p in permutations(v)
             if all(A[p[i]][p[(i+1)%3]] == 1 for i in range(3))) // 3

    # H via DP
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size:
                continue
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                S_prev = S ^ (1 << v)
                total = 0
                for u in range(n):
                    if (S_prev & (1 << u)) and A[u][v]:
                        total += dp.get((S_prev, u), 0)
                if total > 0:
                    dp[(S, v)] = total
    full = (1 << n) - 1
    H = sum(dp.get((full, v), 0) for v in range(n))

    # W(0) via full enumeration
    w0 = 0.0
    for p in permutations(range(n)):
        prod = 1.0
        for i in range(n-1):
            prod *= (A[p[i]][p[i+1]] - 0.5)
        w0 += prod

    data.append({'H': H, 't3': t3, 'w0': w0})
    diff = w0 - H + 3*t3  # Check if C=3 works at n=7
    if trial < 5:
        print(f"  trial={trial}: H={H}, t3={t3}, w0={w0:.4f}, H-3*t3={H-3*t3}, diff={diff:.4f}")

# Fit w0 = H - C*t3
import numpy as np
y = np.array([d['w0'] for d in data])
X = np.array([[d['H'], d['t3'], 1] for d in data])
c = np.linalg.lstsq(X, y, rcond=None)[0]
err = max(abs(y - X @ c))
print(f"\n  w0 = {c[0]:.6f}*H + {c[1]:.6f}*t3 + {c[2]:.6f}  (max_err={err:.4f})")

# Try w0 = H - C*t3 (no constant)
X2 = np.array([[d['H'], d['t3']] for d in data])
c2 = np.linalg.lstsq(X2, y, rcond=None)[0]
err2 = max(abs(y - X2 @ c2))
print(f"  w0 = {c2[0]:.6f}*H + {c2[1]:.6f}*t3  (max_err={err2:.4f})")

# Try w0 = a*H + b*t3 + c*t5
t5_data = []
for trial, d in enumerate(data):
    A = [[0]*n for _ in range(n)]
    random.seed(42)  # Reset seed
    for _ in range(trial):
        # Skip past previous trials
        for i in range(n):
            for j in range(i+1, n):
                random.random()
    random.seed(42 + trial * 1000)  # Use different seed
    # Recompute... actually let me just compute t5 fresh
    pass

# Instead, recompute with t5
random.seed(42)
data2 = []
for trial in range(30):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = sum(1 for v in combinations(range(n), 3)
             for p in permutations(v)
             if all(A[p[i]][p[(i+1)%3]] == 1 for i in range(3))) // 3

    t5 = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts):
            if all(A[p[i]][p[(i+1)%5]] == 1 for i in range(5)):
                t5 += 1
    t5 //= 5

    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size:
                continue
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                S_prev = S ^ (1 << v)
                total = 0
                for u in range(n):
                    if (S_prev & (1 << u)) and A[u][v]:
                        total += dp.get((S_prev, u), 0)
                if total > 0:
                    dp[(S, v)] = total
    full = (1 << n) - 1
    H = sum(dp.get((full, v), 0) for v in range(n))

    w0 = 0.0
    for p in permutations(range(n)):
        prod = 1.0
        for i in range(n-1):
            prod *= (A[p[i]][p[i+1]] - 0.5)
        w0 += prod

    data2.append({'H': H, 't3': t3, 't5': t5, 'w0': w0})

y2 = np.array([d['w0'] for d in data2])
X3 = np.array([[d['H'], d['t3'], d['t5'], 1] for d in data2])
c3 = np.linalg.lstsq(X3, y2, rcond=None)[0]
err3 = max(abs(y2 - X3 @ c3))
print(f"\n  w0 = {c3[0]:.6f}*H + {c3[1]:.6f}*t3 + {c3[2]:.6f}*t5 + {c3[3]:.6f}  (max_err={err3:.4f})")

# Check: at n=5, w0 = H - 3*t3 = 1*(H) + (-3)*t3 + 0*t5 + 0
# At n=7, what is the coefficient of H?

# Also check: is w0 = H - (n-2)*t3 + ... a general pattern?
print(f"\n  At n=5: C = 3 = n-2")
print(f"  At n=7: coefficient of t3 = {c3[1]:.4f}")
print(f"  n-2 = 5. Does C = n-2?")

# The signed backward count interpretation
print(f"\n{'='*70}")
print("SIGNED BACKWARD-PARITY COUNT")
print(f"{'='*70}")
print(f"""
  W(0) = (1/2)^{{n-1}} * sum_P (-1)^{{b(P)}}
  where b(P) = number of backward arcs.

  At n=5: W(0) = (1/16) * [E - O] where E = # even-b paths, O = # odd-b paths.
  Also: E + O = n! = 120.

  So: E - O = 16 * w_0 = 16*(H - 3*t_3).
  And: E = (120 + 16*(H - 3*t_3))/2 = 60 + 8*(H - 3*t_3)
       O = (120 - 16*(H - 3*t_3))/2 = 60 - 8*(H - 3*t_3)
""")

# Verify at n=5
for bits in [0, 1, 4, 15, 31]:
    A = [[0]*5 for _ in range(5)]
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    H = sum(1 for p in permutations(range(5))
            if all(A[p[i]][p[i+1]] == 1 for i in range(4)))
    t3 = sum(1 for v in combinations(range(5), 3)
             for p in permutations(v)
             if all(A[p[i]][p[(i+1)%3]] == 1 for i in range(3))) // 3

    even_b = 0
    odd_b = 0
    for p in permutations(range(5)):
        b = sum(1 for i in range(4) if A[p[i]][p[i+1]] == 0)
        if b % 2 == 0:
            even_b += 1
        else:
            odd_b += 1

    w0 = H - 3*t3
    predicted_diff = 16 * w0

    print(f"  bits={bits:2d}: H={H}, t3={t3}, w0={w0}, E={even_b}, O={odd_b}, E-O={even_b-odd_b}, 16*w0={predicted_diff}")

print("\nDONE")
