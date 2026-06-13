#!/usr/bin/env python3
"""Find bivariate recurrence for N(n,j) or W_u(n) = Σ u^j N(n,j).
Try: does W_u satisfy a recurrence with coefficients depending on n and u?

Key observation: NUD recurrence is NUD(n) = (n-1)*NUD(n-1) + (n-2)*NUD(n-2).
W_u might satisfy a modified version with u-dependent coefficients.
"""
from fractions import Fraction as F
import sys

# Pre-computed N(n,j) data from nud_adj1_dist.py
Nnj = {
    1: [1],
    2: [0, 1],
    3: [0, 2, 1],
    4: [2, 5, 3, 1],
    5: [14, 20, 14, 4, 1],
    6: [90, 115, 72, 26, 5, 1],
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

def Wu(n, u):
    """W_u(n) = Σ_j u^j N(n,j)"""
    return sum(F(u)**j * F(d) for j, d in enumerate(Nnj[n]))

# Check NUD recurrence at u=1:
print("=== Verify NUD recurrence at u=1 ===")
for n in range(3, 16):
    lhs = Wu(n, 1)
    rhs = (n-1)*Wu(n-1, 1) + (n-2)*Wu(n-2, 1)
    print(f"  n={n}: LHS={lhs}, RHS={rhs}, match={lhs==rhs}")

# Try: W_u(n) = (n-1+a*u)*W_u(n-1) + (n-2+b*u)*W_u(n-2)
# At u=0: W_0(n) = N(n,0) = Hertzsprung(n). Does Hertzsprung satisfy (n-1)H(n-1) + (n-2)H(n-2)?
print("\n=== Does Hertzsprung (j=0 column) satisfy NUD recurrence? ===")
for n in range(3, 16):
    h = Wu(n, 0)
    rhs = (n-1)*Wu(n-1, 0) + (n-2)*Wu(n-2, 0)
    print(f"  n={n}: H(n)={h}, (n-1)H(n-1)+(n-2)H(n-2)={rhs}, diff={h-rhs}")

# The difference is the correction we already computed. Try adding more terms.
# Try: W_u(n) = P(n,u)*W_u(n-1) + Q(n,u)*W_u(n-2) + R(n,u)*W_u(n-3)
# with P,Q,R linear in n and polynomial in u.

# Systematic: P = α + βn + γu + δnu, Q = ε + ζn + ηu + θnu, etc.
# Use u = 0, 1, 2, 3, 4 to get 5 sets of recurrences

print("\n=== Check 2-term recurrence W_u(n) = P(n,u)W_u(n-1) + Q(n,u)W_u(n-2) ===")
print("    for each u, find P(n), Q(n) by solving from 2 equations")
for u in range(5):
    print(f"\n  u={u}:")
    # For each n, W_u(n) = P*W_u(n-1) + Q*W_u(n-2)
    # From n=3,4: solve for P,Q
    # W(3) = P*W(2) + Q*W(1)
    # W(4) = P*W(3) + Q*W(2)
    w = [Wu(n, u) for n in range(1, 16)]
    # w[i] = W_u(i+1)

    for n in range(5, 16):
        # Use n-1 and n to solve for P, Q
        # W(n) = P*W(n-1) + Q*W(n-2)
        # W(n-1) = P*W(n-2) + Q*W(n-3)
        wn = w[n-1]; wn1 = w[n-2]; wn2 = w[n-3]; wn3 = w[n-4] if n >= 4 else 0

        det = wn1*wn3 - wn2*wn2
        if det != 0:
            P = F(wn*wn3 - wn1*wn2, det)  # Cramer's rule... actually let me redo this

        # W(n) = P*W(n-1) + Q*W(n-2) ... (1)
        # W(n+1) = P*W(n) + Q*W(n-1) ... (2)
        # Better: just solve from two consecutive n values
        pass

    # Just compute the ratio: if Q=0, then P = W(n)/W(n-1)
    for n in range(3, 16):
        if w[n-2] != 0:
            r = F(w[n-1], w[n-2])
            print(f"    W_{u}({n})/W_{u}({n-1}) = {float(r):.6f}")

# Let's try the 3-term recurrence more carefully for u=2
print("\n=== 3-term recurrence search for W_2(n) ===")
# W(n) = (an+b)W(n-1) + (cn+d)W(n-2) + (en+f)W(n-3)
# 6 unknowns, use n=4..9
w2 = [Wu(n, 2) for n in range(1, 16)]

import numpy as np

# Build matrix for n=4..9
rows = []
for n in range(4, 10):
    i = n - 1  # index into w2
    row = [n*w2[i-1], w2[i-1], n*w2[i-2], w2[i-2], n*w2[i-3], w2[i-3]]
    rows.append(row)
rhs = [w2[n-1] for n in range(4, 10)]

# Convert to Fraction for exact arithmetic
mat_f = []
for i, row in enumerate(rows):
    mat_f.append([F(x) for x in row] + [F(rhs[i])])

# Gaussian elimination
nv = 6
for col in range(nv):
    pivot = None
    for r in range(col, nv):
        if mat_f[r][col] != 0:
            pivot = r; break
    if pivot is None: continue
    mat_f[col], mat_f[pivot] = mat_f[pivot], mat_f[col]
    for r in range(nv):
        if r != col and mat_f[r][col] != 0:
            fac = F(mat_f[r][col], mat_f[col][col])
            for c in range(nv+1):
                mat_f[r][c] -= fac * mat_f[col][c]

sol = [F(mat_f[i][nv], mat_f[i][i]) for i in range(nv)]
a, b, c, d, e, f_val = sol
print(f"  W(n) = ({float(a):.6f}n + {float(b):.6f})W(n-1) + ({float(c):.6f}n + {float(d):.6f})W(n-2) + ({float(e):.6f}n + {float(f_val):.6f})W(n-3)")
print(f"  Exact: a={a}, b={b}, c={c}, d={d}, e={e}, f={f_val}")

# Verify
for n in range(4, 16):
    i = n - 1
    pred = (a*n+b)*w2[i-1] + (c*n+d)*w2[i-2] + (e*n+f_val)*w2[i-3]
    match = pred == w2[i]
    if not match:
        print(f"  n={n}: FAIL (diff={pred - w2[i]})")
    else:
        print(f"  n={n}: OK")

# Now try for general u: same form W_u(n) = (an+b)W_u(n-1) + (cn+d)W_u(n-2) + (en+f)W_u(n-3)
# where a,b,c,d,e,f may depend on u
print("\n=== 3-term recurrence for each u ===")
for u in range(6):
    wu = [Wu(n, u) for n in range(1, 16)]
    rows = []
    for n in range(4, 10):
        i = n - 1
        row = [F(n*wu[i-1]), F(wu[i-1]), F(n*wu[i-2]), F(wu[i-2]), F(n*wu[i-3]), F(wu[i-3])]
        rows.append(row)
    rhs = [F(wu[n-1]) for n in range(4, 10)]

    mat_f = [list(row) + [rhs[j]] for j, row in enumerate(rows)]
    for col in range(6):
        pivot = None
        for r in range(col, 6):
            if mat_f[r][col] != 0:
                pivot = r; break
        if pivot is None:
            print(f"  u={u}: degenerate at col {col}")
            continue
        mat_f[col], mat_f[pivot] = mat_f[pivot], mat_f[col]
        for r in range(6):
            if r != col and mat_f[r][col] != 0:
                fac = F(mat_f[r][col], mat_f[col][col])
                for c in range(7):
                    mat_f[r][c] -= fac * mat_f[col][c]

    sol = [F(mat_f[i][6], mat_f[i][i]) for i in range(6)]
    a, b, c, d, e, f_val = sol

    # Verify
    ok = True
    for n in range(4, 16):
        i = n - 1
        pred = (a*n+b)*wu[i-1] + (c*n+d)*wu[i-2] + (e*n+f_val)*wu[i-3]
        if pred != wu[i]:
            ok = False
            break

    status = "OK (all n=4..15)" if ok else f"FAIL at n={n}"
    print(f"  u={u}: ({float(a):.4f}n+{float(b):.4f})W(n-1) + ({float(c):.4f}n+{float(d):.4f})W(n-2) + ({float(e):.4f}n+{float(f_val):.4f})W(n-3) [{status}]")

# Try 4-term recurrence with constant coefficients (not n-dependent)
print("\n=== 4-term constant-coeff recurrence ===")
for u in [0, 1, 2]:
    wu = [Wu(n, u) for n in range(1, 16)]
    # W(n) = a*W(n-1) + b*W(n-2) + c*W(n-3) + d*W(n-4)
    # Use n=5..8
    rows = []
    for n in range(5, 9):
        i = n - 1
        row = [F(wu[i-1]), F(wu[i-2]), F(wu[i-3]), F(wu[i-4])]
        rows.append(row)
    rhs = [F(wu[n-1]) for n in range(5, 9)]

    mat_f = [list(row) + [rhs[j]] for j, row in enumerate(rows)]
    for col in range(4):
        pivot = None
        for r in range(col, 4):
            if mat_f[r][col] != 0: pivot = r; break
        if pivot is None: continue
        mat_f[col], mat_f[pivot] = mat_f[pivot], mat_f[col]
        for r in range(4):
            if r != col and mat_f[r][col] != 0:
                fac = F(mat_f[r][col], mat_f[col][col])
                for c in range(5):
                    mat_f[r][c] -= fac * mat_f[col][c]
    sol = [F(mat_f[i][4], mat_f[i][i]) for i in range(4)]

    ok = True
    for n in range(5, 16):
        i = n - 1
        pred = sum(sol[k]*wu[i-1-k] for k in range(4))
        if pred != wu[i]: ok = False; break

    print(f"  u={u}: {[float(s) for s in sol]} {'OK' if ok else f'FAIL at n={n}'}")
