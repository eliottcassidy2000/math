#!/usr/bin/env python3
"""Verify Hertzsprung recurrence and find bivariate generalization for N(n,j)."""
from fractions import Fraction as F

# N(n,j) data
Nnj = {
    1: [1], 2: [0, 1], 3: [0, 2, 1], 4: [2, 5, 3, 1],
    5: [14, 20, 14, 4, 1], 6: [90, 115, 72, 26, 5, 1],
    7: [646, 790, 467, 168, 41, 6, 1],
    8: [5242, 6217, 3557, 1285, 319, 59, 7, 1],
    9: [47622, 55160, 30968, 11120, 2834, 536, 80, 8, 1],
    10: [479306, 545135, 301970, 107918, 27752, 5432, 830, 104, 9, 1],
    11: [5296790, 5938490, 3255789, 1158624, 299454, 59580, 9450, 1212, 131, 10, 1],
    12: [63779034, 70686805, 38434607, 13626345, 3530514, 710034, 114894, 15312, 1693, 161, 11, 1],
    13: [831283558, 912660508, 492864354, 174148604, 45172191, 9147096, 1501020, 204456, 23495, 2284, 194, 12, 1],
    14: [11661506218, 12702694075, 6820525996, 2402630714, 623426927, 126812523, 21010320, 2907630, 341809, 34529, 2996, 230, 13, 1],
    15: [175203184374, 189579135710, 101292697671, 35583328840, 9231722935, 1883522306, 314193675, 43976556, 5258009, 543586, 48997, 3840, 269, 14, 1],
}

def N(n, j):
    if n not in Nnj: return 0
    d = Nnj[n]
    return d[j] if 0 <= j < len(d) else 0

# Verify Hertzsprung recurrence: a(n) = (n+1)a(n-1) - (n-2)a(n-2) - (n-5)a(n-3) + (n-3)a(n-4)
print("=== Verify Hertzsprung recurrence for j=0 column ===")
for n in range(5, 16):
    pred = (n+1)*N(n-1,0) - (n-2)*N(n-2,0) - (n-5)*N(n-3,0) + (n-3)*N(n-4,0)
    actual = N(n, 0)
    print(f"  n={n}: pred={pred}, actual={actual}, {'OK' if pred==actual else 'FAIL'}")

# Now check if this recurrence works for other j values
print("\n=== Check Hertzsprung recurrence for each j ===")
for j in range(12):
    ok = True
    for n in range(5, 16):
        pred = (n+1)*N(n-1,j) - (n-2)*N(n-2,j) - (n-5)*N(n-3,j) + (n-3)*N(n-4,j)
        actual = N(n, j)
        if pred != actual:
            ok = False
            first_fail_n = n
            diff = pred - actual
            break
    if ok:
        print(f"  j={j}: OK for all n=5..15")
    else:
        print(f"  j={j}: FAIL at n={first_fail_n}, diff={diff}")

# The recurrence probably fails for j>0. Let me compute the correction.
print("\n=== Corrections: N(n,j) - [(n+1)N(n-1,j) - (n-2)N(n-2,j) - (n-5)N(n-3,j) + (n-3)N(n-4,j)] ===")
for n in range(5, 14):
    corr = []
    for j in range(min(n, 10)):
        pred = (n+1)*N(n-1,j) - (n-2)*N(n-2,j) - (n-5)*N(n-3,j) + (n-3)*N(n-4,j)
        corr.append(N(n,j) - pred)
    while corr and corr[-1] == 0:
        corr.pop()
    print(f"  n={n}: {corr}")

# Maybe the correction involves N(n-k, j-1) terms?
# Try: N(n,j) = Hertzsprung_recurrence(N(·,j)) + α*(something involving j-1)
# Let me check if correction = c₁·N(n-1,j-1) + c₂·N(n-2,j-1) + ...
print("\n=== Correction analysis: is it proportional to N(·,j-1)? ===")
for n in range(5, 14):
    for j in range(1, min(n-1, 6)):
        pred = (n+1)*N(n-1,j) - (n-2)*N(n-2,j) - (n-5)*N(n-3,j) + (n-3)*N(n-4,j)
        corr = N(n,j) - pred
        # Check: corr / N(n-1, j-1)
        if N(n-1, j-1) != 0:
            ratio = F(corr, N(n-1, j-1))
            print(f"  n={n}, j={j}: corr={corr}, N(n-1,j-1)={N(n-1,j-1)}, ratio={float(ratio):.6f}")

# Try: corr = α·N(n-1,j-1) + β·N(n-2,j-1) + γ·N(n-3,j-1) + δ·N(n-4,j-1)
# with α,β,γ,δ constants (not depending on n or j)
print("\n=== Fit: corr = α·N(n-1,j-1) + β·N(n-2,j-1) + γ·N(n-3,j-1) + δ·N(n-4,j-1) ===")
# Collect equations
eqs = []
for n in range(5, 16):
    for j in range(1, min(n-1, 8)):
        pred = (n+1)*N(n-1,j) - (n-2)*N(n-2,j) - (n-5)*N(n-3,j) + (n-3)*N(n-4,j)
        corr = N(n,j) - pred
        row = [F(N(n-1,j-1)), F(N(n-2,j-1)), F(N(n-3,j-1)), F(N(n-4,j-1))]
        eqs.append((row, F(corr)))

# Solve from first 4 equations
mat = [list(eqs[i][0]) + [eqs[i][1]] for i in range(4)]
nv = 4
for col in range(nv):
    pivot = None
    for r in range(col, nv):
        if mat[r][col] != 0:
            pivot = r; break
    if pivot is None:
        print(f"  Degenerate at col {col}")
        continue
    mat[col], mat[pivot] = mat[pivot], mat[col]
    for r in range(nv):
        if r != col and mat[r][col] != 0:
            fac = F(mat[r][col], mat[col][col])
            for c in range(nv+1):
                mat[r][c] -= fac * mat[col][c]

sol_const = [F(mat[i][nv], mat[i][i]) for i in range(nv)]
print(f"  α={sol_const[0]}, β={sol_const[1]}, γ={sol_const[2]}, δ={sol_const[3]}")

# Verify
fail_count = 0
for i, (row, rhs) in enumerate(eqs):
    pred = sum(sol_const[k] * row[k] for k in range(4))
    if pred != rhs:
        if fail_count < 3:
            # Find which (n,j) this corresponds to
            print(f"    FAIL at eq {i}: pred={pred}, rhs={rhs}")
        fail_count += 1
if fail_count == 0:
    print(f"  ALL {len(eqs)} equations verified! BIVARIATE RECURRENCE FOUND!")
else:
    print(f"  {fail_count} failures out of {len(eqs)}")

# Try with n-dependent coefficients: corr = (α₁n+β₁)·N(n-1,j-1) + (α₂n+β₂)·N(n-2,j-1) + ...
print("\n=== Fit with n-linear coeffs: corr = (α₁n+β₁)N(n-1,j-1) + (α₂n+β₂)N(n-2,j-1) + (α₃n+β₃)N(n-3,j-1) + (α₄n+β₄)N(n-4,j-1) ===")
eqs2 = []
for n in range(5, 16):
    for j in range(1, min(n-1, 8)):
        pred = (n+1)*N(n-1,j) - (n-2)*N(n-2,j) - (n-5)*N(n-3,j) + (n-3)*N(n-4,j)
        corr = N(n,j) - pred
        row = [F(n*N(n-1,j-1)), F(N(n-1,j-1)),
               F(n*N(n-2,j-1)), F(N(n-2,j-1)),
               F(n*N(n-3,j-1)), F(N(n-3,j-1)),
               F(n*N(n-4,j-1)), F(N(n-4,j-1))]
        eqs2.append((row, F(corr)))

mat2 = [list(eqs2[i][0]) + [eqs2[i][1]] for i in range(8)]
nv = 8
for col in range(nv):
    pivot = None
    for r in range(col, nv):
        if mat2[r][col] != 0:
            pivot = r; break
    if pivot is None:
        print(f"  Degenerate at col {col}")
        continue
    mat2[col], mat2[pivot] = mat2[pivot], mat2[col]
    for r in range(nv):
        if r != col and mat2[r][col] != 0:
            fac = F(mat2[r][col], mat2[col][col])
            for c in range(nv+1):
                mat2[r][c] -= fac * mat2[col][c]

sol_lin = [F(mat2[i][nv], mat2[i][i]) for i in range(nv)]
labels = ['α₁', 'β₁', 'α₂', 'β₂', 'α₃', 'β₃', 'α₄', 'β₄']
for label, s in zip(labels, sol_lin):
    print(f"  {label} = {s} = {float(s):.6f}")

fail_count = 0
for i, (row, rhs) in enumerate(eqs2):
    pred = sum(sol_lin[k] * row[k] for k in range(8))
    if pred != rhs:
        if fail_count < 3:
            print(f"    FAIL at eq {i}: diff={pred - rhs}")
        fail_count += 1
if fail_count == 0:
    print(f"  ALL {len(eqs2)} equations verified! BIVARIATE RECURRENCE WITH n-LINEAR COEFFS FOUND!")
else:
    print(f"  {fail_count} failures out of {len(eqs2)}")

# Try mixed: some n-dependent, some not.
# Actually let me try the simplest possible: correction = α·N(n-1,j-1) + β·N(n-2,j-1)
print("\n=== Simplest: corr = α·N(n-1,j-1) + β·N(n-2,j-1) ===")
eqs3 = []
for n in range(5, 16):
    for j in range(1, min(n-1, 8)):
        pred = (n+1)*N(n-1,j) - (n-2)*N(n-2,j) - (n-5)*N(n-3,j) + (n-3)*N(n-4,j)
        corr = N(n,j) - pred
        row = [F(N(n-1,j-1)), F(N(n-2,j-1))]
        eqs3.append((row, F(corr)))

mat3 = [list(eqs3[i][0]) + [eqs3[i][1]] for i in range(2)]
nv = 2
for col in range(nv):
    pivot = None
    for r in range(col, nv):
        if mat3[r][col] != 0:
            pivot = r; break
    if pivot is None: continue
    mat3[col], mat3[pivot] = mat3[pivot], mat3[col]
    for r in range(nv):
        if r != col and mat3[r][col] != 0:
            fac = F(mat3[r][col], mat3[col][col])
            for c in range(nv+1):
                mat3[r][c] -= fac * mat3[col][c]

sol3 = [F(mat3[i][nv], mat3[i][i]) for i in range(nv)]
print(f"  α={sol3[0]}, β={sol3[1]}")

fail_count = 0
for i, (row, rhs) in enumerate(eqs3):
    pred = sum(sol3[k] * row[k] for k in range(2))
    if pred != rhs:
        if fail_count < 3:
            print(f"    FAIL at eq {i}")
        fail_count += 1
if fail_count == 0:
    print(f"  ALL {len(eqs3)} equations verified!")
else:
    print(f"  {fail_count} failures out of {len(eqs3)}")
