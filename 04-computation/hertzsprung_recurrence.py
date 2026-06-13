#!/usr/bin/env python3
"""Find the recurrence for Hertzsprung numbers A002464 from data,
then check if N(n,j) satisfies a generalization."""
from fractions import Fraction as F

# Hertzsprung: H(n) = N(n,0) = perms of [n] with no unit ascent or descent
H = {1: 1, 2: 0, 3: 0, 4: 2, 5: 14, 6: 90, 7: 646, 8: 5242,
     9: 47622, 10: 479306, 11: 5296790, 12: 63779034,
     13: 831283558, 14: 11661506218, 15: 175203184374}

# Known recurrence for A002464 (from OEIS):
# a(n) = (n+1)*a(n-1) - (n-2)*a(n-2) - 2*(n-2)*a(n-3) + (n-3)*a(n-4) + (n-3)*a(n-5)
# Let me verify this
print("=== Verify known Hertzsprung recurrence ===")
print("a(n) = (n+1)*a(n-1) - (n-2)*a(n-2) - 2*(n-2)*a(n-3) + (n-3)*a(n-4) + (n-3)*a(n-5)")
for n in range(6, 16):
    pred = (n+1)*H[n-1] - (n-2)*H[n-2] - 2*(n-2)*H[n-3] + (n-3)*H[n-4] + (n-3)*H[n-5]
    match = pred == H[n]
    if not match:
        print(f"  n={n}: FAIL (pred={pred}, actual={H[n]}, diff={pred-H[n]})")
    else:
        print(f"  n={n}: OK")

# That doesn't work. Let me find the actual recurrence.
# Try: a(n) = Σ_{k=1}^{d} p_k(n) * a(n-k) where p_k is linear in n
# Try depth 3
print("\n=== Search for Hertzsprung recurrence, depth 3, quadratic coeffs ===")
# a(n) = (α₁n²+β₁n+γ₁)a(n-1) + (α₂n²+β₂n+γ₂)a(n-2) + (α₃n²+β₃n+γ₃)a(n-3)
# 9 unknowns, use n=4..12

rows = []
for n in range(4, 13):
    row = [F(n**2*H[n-1]), F(n*H[n-1]), F(H[n-1]),
           F(n**2*H[n-2]), F(n*H[n-2]), F(H[n-2]),
           F(n**2*H[n-3]), F(n*H[n-3]), F(H[n-3])]
    rows.append((row, F(H[n])))

mat = [list(row) + [rhs] for row, rhs in rows]
nv = 9
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

sol = [F(mat[i][nv], mat[i][i]) for i in range(nv)]
for i, s in enumerate(sol):
    print(f"  sol[{i}] = {s} = {float(s):.10f}")

# Verify on n=13,14,15
for n in [13, 14, 15]:
    pred = F(0)
    for k in range(3):
        pk = sol[3*k]*n**2 + sol[3*k+1]*n + sol[3*k+2]
        pred += pk * H[n-1-k]
    match = pred == H[n]
    print(f"  n={n}: {'OK' if match else 'FAIL'}")

# Now try 5-term, linear coefficients
print("\n=== 5-term recurrence with linear coefficients ===")
# a(n) = (a₁n+b₁)a(n-1) + (a₂n+b₂)a(n-2) + ... + (a₅n+b₅)a(n-5)
# 10 unknowns, use n=6..15
rows = []
for n in range(6, 16):
    row = []
    for k in range(1, 6):
        row.extend([F(n*H[n-k]), F(H[n-k])])
    rows.append((row, F(H[n])))

mat = [list(row) + [rhs] for row, rhs in rows]
nv = 10
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

sol5 = [F(mat[i][nv], mat[i][i]) for i in range(nv)]
print("  Coefficients:")
for k in range(5):
    a_coeff, b_coeff = sol5[2*k], sol5[2*k+1]
    print(f"    P_{k+1}(n) = {a_coeff}*n + {b_coeff} = {float(a_coeff):.6f}*n + {float(b_coeff):.6f}")

# Check if these are simple rationals
print("\n  Exact coefficients:")
for k in range(5):
    a_coeff, b_coeff = sol5[2*k], sol5[2*k+1]
    print(f"    ({a_coeff})n + ({b_coeff})")

# Now try the same 5-term recurrence for the full N(n,j) triangle
print("\n=== Check: does the same recurrence work for general j? ===")
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

def get_nnj(n, j):
    if n not in Nnj:
        return 0
    d = Nnj[n]
    return d[j] if j < len(d) else 0

# For each j, check the 5-term Hertzsprung recurrence
for j in range(8):
    ok = True
    for n in range(6, 16):
        pred = F(0)
        for k in range(5):
            pk = sol5[2*k]*n + sol5[2*k+1]
            pred += pk * get_nnj(n-1-k, j)
        actual = F(get_nnj(n, j))
        if pred != actual:
            ok = False
            break
    print(f"  j={j}: {'OK' if ok else f'FAIL at n={n}, diff={pred - actual}'}")

# The recurrence probably fails for j>0. Let me check if there's a j-shifted version.
# Try: N(n,j) = Σ_{k=1}^{5} P_k(n) N(n-k, j) + Σ_{k=1}^{5} Q_k(n) N(n-k, j-1)
# This doubles the unknowns but might capture the bivariate structure.

print("\n=== Bivariate 3-term recurrence with j and j-1 ===")
# N(n,j) = (a₁n+b₁)N(n-1,j) + (a₂n+b₂)N(n-2,j) + (c₁n+d₁)N(n-1,j-1) + (c₂n+d₂)N(n-2,j-1)
# 8 unknowns. Use multiple (n,j) pairs.

# Collect equations
eqs = []
for n in range(4, 16):
    max_j = min(n-1, 8)
    for j in range(0, max_j):
        # N(n,j) = (a₁n+b₁)N(n-1,j) + (a₂n+b₂)N(n-2,j) + (c₁n+d₁)N(n-1,j-1) + (c₂n+d₂)N(n-2,j-1)
        lhs = F(get_nnj(n, j))
        row = [
            F(n*get_nnj(n-1, j)), F(get_nnj(n-1, j)),
            F(n*get_nnj(n-2, j)), F(get_nnj(n-2, j)),
            F(n*get_nnj(n-1, j-1)), F(get_nnj(n-1, j-1)),
            F(n*get_nnj(n-2, j-1)), F(get_nnj(n-2, j-1)),
        ]
        eqs.append((row, lhs))

# Use first 8 equations, verify rest
nv = 8
mat = [list(eqs[i][0]) + [eqs[i][1]] for i in range(nv)]
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

sol_biv = [F(mat[i][nv], mat[i][i]) for i in range(nv)]
print("  Coefficients:")
labels = ['a1', 'b1', 'a2', 'b2', 'c1', 'd1', 'c2', 'd2']
for i, (label, s) in enumerate(zip(labels, sol_biv)):
    print(f"    {label} = {s} = {float(s):.6f}")

# Verify on all equations
fail_count = 0
for i, (row, lhs) in enumerate(eqs):
    pred = sum(sol_biv[k] * row[k] for k in range(nv))
    if pred != lhs:
        if fail_count < 5:
            print(f"    FAIL at eq {i}: pred={pred}, actual={lhs}, diff={pred-lhs}")
        fail_count += 1
if fail_count == 0:
    print(f"    ALL {len(eqs)} equations verified!")
else:
    print(f"    {fail_count} failures out of {len(eqs)}")
