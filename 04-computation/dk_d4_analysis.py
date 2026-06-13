#!/usr/bin/env python3
"""
dk_d4_analysis.py — What does D_4 depend on at n=7?

D_2, D_3 are exactly linear in t3. D_4 is not.
D_4 depends on t3, t5, t7, and probably some intersection statistic.

Candidate invariants:
- t3: number of 3-cycles
- t5: number of 5-cycles
- t7: number of 7-cycles
- bc33: number of vertex-disjoint 3-cycle pairs (alpha_2 of Omega_3)
- sum_v in(v)*out(v): score-based invariant
- sum_v in(v)^2: score variance

Author: opus-2026-03-07-S43b
"""
from itertools import permutations, combinations
import math
import random
import numpy as np

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

def compute_Dk(A, n):
    D = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        for k in range(n):
            D[k] += math.comb(fwd, k)
    return D

def count_cycles(A, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            t3 += 1

    t5 = 0
    cycles5 = []
    for quint in combinations(range(n), 5):
        cnt = 0
        for perm in permutations(quint):
            if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                cnt += 1
        t5 += cnt // 10

    t7 = 0
    if n == 7:
        cnt = 0
        for perm in permutations(range(n)):
            if all(A[perm[i]][perm[(i+1)%7]] for i in range(7)):
                cnt += 1
        t7 = cnt // 14

    return t3, t5, t7

def count_bc33(A, n):
    """Count vertex-disjoint 3-cycle pairs."""
    # First enumerate all 3-cycles as vertex sets
    three_cycles = []
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            three_cycles.append(set(triple))

    bc33 = 0
    for i in range(len(three_cycles)):
        for j in range(i+1, len(three_cycles)):
            if not (three_cycles[i] & three_cycles[j]):
                bc33 += 1
    return bc33

def score_stats(A, n):
    scores = [sum(A[i][j] for j in range(n)) for i in range(n)]
    sum_io = sum(s * (n-1-s) for s in scores)
    sum_s2 = sum(s*s for s in scores)
    return sum_io, sum_s2, tuple(sorted(scores))

n = 7
m = n*(n-1)//2
random.seed(456)

print(f"=== n={n}: D_4 dependency analysis ===")

data = []
for trial in range(100):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    D = compute_Dk(A, n)
    t3, t5, t7 = count_cycles(A, n)
    bc33 = count_bc33(A, n)
    sum_io, sum_s2, scores = score_stats(A, n)
    data.append((D[4], D[5], D[6], t3, t5, t7, bc33, sum_io, sum_s2))

D4 = np.array([d[0] for d in data], dtype=float)
D5 = np.array([d[1] for d in data], dtype=float)
H  = np.array([d[2] for d in data], dtype=float)
t3 = np.array([d[3] for d in data], dtype=float)
t5 = np.array([d[4] for d in data], dtype=float)
t7 = np.array([d[5] for d in data], dtype=float)
bc = np.array([d[6] for d in data], dtype=float)
sio= np.array([d[7] for d in data], dtype=float)
ss2= np.array([d[8] for d in data], dtype=float)

# Try various models for D_4
models = {
    "t3": np.column_stack([np.ones_like(t3), t3]),
    "t3,t5": np.column_stack([np.ones_like(t3), t3, t5]),
    "t3,t5,t7": np.column_stack([np.ones_like(t3), t3, t5, t7]),
    "t3,bc33": np.column_stack([np.ones_like(t3), t3, bc]),
    "t3,t5,bc33": np.column_stack([np.ones_like(t3), t3, t5, bc]),
    "t3,t5,t7,bc33": np.column_stack([np.ones_like(t3), t3, t5, t7, bc]),
    "t3,t5,t7,bc33,sio": np.column_stack([np.ones_like(t3), t3, t5, t7, bc, sio]),
    "t3,t5,t7,bc33,ss2": np.column_stack([np.ones_like(t3), t3, t5, t7, bc, ss2]),
    "t3,bc33,sio": np.column_stack([np.ones_like(t3), t3, bc, sio]),
}

print("\n--- D_4 models ---")
for name, X in sorted(models.items(), key=lambda x: x[0]):
    c, _, _, _ = np.linalg.lstsq(X, D4, rcond=None)
    pred = X @ c
    maxerr = max(abs(D4 - pred))
    rmse = np.sqrt(np.mean((D4-pred)**2))
    names = ['1'] + name.split(',')
    formula = ' + '.join(f'{v:.2f}*{n}' for n, v in zip(names, c))
    print(f"  {name:25s}: maxerr={maxerr:7.2f}, rmse={rmse:6.2f}  | {formula}")

print("\n--- D_5 models ---")
for name, X in sorted(models.items(), key=lambda x: x[0]):
    c, _, _, _ = np.linalg.lstsq(X, D5, rcond=None)
    pred = X @ c
    maxerr = max(abs(D5 - pred))
    rmse = np.sqrt(np.mean((D5-pred)**2))
    names = ['1'] + name.split(',')
    formula = ' + '.join(f'{v:.2f}*{n}' for n, v in zip(names, c))
    print(f"  {name:25s}: maxerr={maxerr:7.2f}, rmse={rmse:6.2f}  | {formula}")

# Check: D_4 - D_5 = D_4 - D_5 pattern
print("\n--- D_4 - D_5 ---")
diff45 = D4 - D5
for name, X in [("t3", np.column_stack([np.ones_like(t3), t3])),
                ("t3,t5", np.column_stack([np.ones_like(t3), t3, t5])),
                ("t3,t5,t7", np.column_stack([np.ones_like(t3), t3, t5, t7])),
                ("t3,t5,t7,bc33", np.column_stack([np.ones_like(t3), t3, t5, t7, bc]))]:
    c, _, _, _ = np.linalg.lstsq(X, diff45, rcond=None)
    pred = X @ c
    maxerr = max(abs(diff45 - pred))
    names = ['1'] + name.split(',')
    formula = ' + '.join(f'{v:.2f}*{n}' for n, v in zip(names, c))
    print(f"  {name:25s}: maxerr={maxerr:7.2f} | {formula}")

# Key observation from skeleton: D_4 - D_5 depends on more than t3
# because D_4 has the t5, bc33 contributions while D_5 also does

# Check: does H = D_6 follow a pattern?
print("\n--- H models ---")
for name, X in [("t3", np.column_stack([np.ones_like(t3), t3])),
                ("t3,t5", np.column_stack([np.ones_like(t3), t3, t5])),
                ("t3,t5,t7", np.column_stack([np.ones_like(t3), t3, t5, t7])),
                ("t3,t5,t7,bc33", np.column_stack([np.ones_like(t3), t3, t5, t7, bc]))]:
    c, _, _, _ = np.linalg.lstsq(X, H, rcond=None)
    pred = X @ c
    maxerr = max(abs(H - pred))
    names = ['1'] + name.split(',')
    formula = ' + '.join(f'{v:.2f}*{n}' for n, v in zip(names, c))
    print(f"  {name:25s}: maxerr={maxerr:7.2f} | {formula}")

# OCF says H = 1 + 2*alpha1 + 4*alpha2 + 8*alpha3
# where alpha_k = #independent k-sets in Omega(T)
# alpha_1 = t3 + t5 + t7
# So H - 2*(t3+t5+t7) - 1 = 4*alpha2 + 8*alpha3
# Let's check
alpha1 = t3 + t5 + t7
residual = H - 1 - 2*alpha1
print(f"\n--- H - 1 - 2*alpha1 = 4*alpha2 + 8*alpha3 ---")
print(f"  min={min(residual):.0f}, max={max(residual):.0f}")
print(f"  All divisible by 4? {all(r % 4 == 0 for r in residual)}")

# alpha2 includes 3-5 pairs, 3-7 pairs, 5-7 pairs, and 3-3 pairs (=bc33)
# Let's check: is residual/4 linear in bc33?
res4 = residual / 4
X_bc = np.column_stack([np.ones_like(bc), bc])
c, _, _, _ = np.linalg.lstsq(X_bc, res4, rcond=None)
pred = X_bc @ c
maxerr = max(abs(res4 - pred))
print(f"  alpha2 + 2*alpha3 = {c[0]:.2f} + {c[1]:.2f}*bc33, maxerr={maxerr:.2f}")
