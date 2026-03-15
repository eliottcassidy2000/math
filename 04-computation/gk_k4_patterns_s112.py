#!/usr/bin/env python3
"""
gk_k4_patterns_s112.py — Extract patterns from k>=4 correction structure
kind-pasteur-2026-03-15-S112

From gk_k4_investigation_s112.py:
- Delta_4(m) = 8*C(m-1,4) exactly (clean, degree 4)
- Delta_k for k>=5 has components at ALL degrees 3 through k
- The alternating tail alpha_{k,j} for j >> k approaches +-|g_k(0)|
- alpha_{k,3} sequence: -160, -24960, -6314240, -2144267840, ...

Questions to investigate:
1. Exact formula for alpha_{k,j} in terms of k and j
2. Relationship between alpha_{k,j} and eigenvalues of M(x)
3. Generating function / recurrence for the correction terms
"""

from fractions import Fraction
from math import factorial, comb

def transfer_gk_values(k, m_max):
    results = {}
    for m in range(0, m_max + 1):
        n = m + 2*k
        num_edges = n - 2
        k_max = (n - 1) // 2
        if k > k_max:
            continue
        state = [[Fraction(1)] + [Fraction(0)] * k_max,
                 [Fraction(0)] * (k_max + 1),
                 [Fraction(0)] * (k_max + 1)]
        for step in range(num_edges):
            A, B, C = state
            nA = [A[i] + C[i] for i in range(k_max + 1)]
            nB = [Fraction(0)] + [2*A[i] + C[i] for i in range(k_max)]
            nC = list(B)
            state = [nA, nB, nC]
        total = [state[0][i] + state[1][i] + state[2][i] for i in range(k_max + 1)]
        if k < len(total) and total[k] != 0:
            results[m] = total[k] / 2
        else:
            results[m] = Fraction(0)
    return results

def opus_gk(k, m):
    C = {3:(2,0,1,0), 4:(10,-33,50,-24), 5:(388,-2040,3431,-1776),
         6:(69660,-380445,653748,-342960), 7:(19826270,-109486152,189674605,-100014720),
         8:(7309726742,-40641958545,70757788486,-37425556680),
         9:(3262687720240,-18232387983408,31858349908595,-16888649645424)}
    if k == 1: return Fraction(m)
    if k == 2: return Fraction(m * m)
    if k not in C: return None
    a,b,c,d = C[k]
    return Fraction(a*m**3 + b*m**2 + c*m + d, 3)

# ============================================================
# PATTERN 1: The alternating tail and g_k(0)
# ============================================================

print("="*60)
print("PATTERN 1: Alternating tail and g_k(0)")
print("="*60)

for k in range(4, 10):
    gk0 = opus_gk(k, 0)
    print(f"\nk={k}: g_k(0) = {gk0}")

    # Compute alpha_{k,j} for high j
    tm = transfer_gk_values(k, k + 6)
    deltas = {}
    for m in range(0, k + 5):
        if m not in tm:
            continue
        op = opus_gk(k, m)
        if op is not None:
            deltas[m] = tm[m] - op

    vals = [deltas.get(m, Fraction(0)) for m in range(k + 5)]
    cur = list(vals)
    alphas = []
    for j in range(len(cur)):
        alphas.append(cur[0])
        cur = [cur[i+1] - cur[i] for i in range(len(cur) - 1)]
        if not cur:
            break

    # Show high-j alphas and compare to -g_k(0)
    neg_gk0 = -gk0
    for j in range(max(0, k-2), len(alphas)):
        if alphas[j] != 0:
            ratio = alphas[j] / neg_gk0 if neg_gk0 != 0 else "?"
            sign = "+" if alphas[j] > 0 else "-"
            print(f"  alpha_{j} = {alphas[j]}, -g_k(0) = {neg_gk0}, "
                  f"ratio = {ratio}, (-1)^j = {(-1)**j}")

# ============================================================
# PATTERN 2: The k=4 case in detail
# ============================================================

print("\n" + "="*60)
print("PATTERN 2: k=4 correction — why 8*C(m-1,4)?")
print("="*60)

# The transfer matrix at k=4 counts 4-matchings of path P_{m+7}.
# The dominant eigenvalue gives g_4 ~ 2^3/4! * m^4 + lower = m^4/3 + lower
# The opus g_4 = (10m^3-33m^2+50m-24)/3
# So the degree-4 part of TM g_4 is the correction.

# Let's check: TM g_4(m) = opus_g_4(m) + 8*C(m-1,4)
# opus_g_4(m) = (10m^3-33m^2+50m-24)/3
# 8*C(m-1,4) = 8*(m-1)(m-2)(m-3)(m-4)/24 = (m-1)(m-2)(m-3)(m-4)/3

# TM g_4 = [(10m^3-33m^2+50m-24) + (m-1)(m-2)(m-3)(m-4)] / 3
# Expand (m-1)(m-2)(m-3)(m-4) = m^4 - 10m^3 + 35m^2 - 50m + 24
# TM g_4 = [m^4 - 10m^3 + 35m^2 - 50m + 24 + 10m^3 - 33m^2 + 50m - 24] / 3
#         = [m^4 + 2m^2] / 3 = m^2(m^2+2)/3

tm_g4_formula = lambda m: m*m*(m*m+2)//3
print("TM g_4(m) = m^2(m^2+2)/3:")
for m in range(1, 10):
    tm_val = transfer_gk_values(4, m)[m]
    formula_val = Fraction(m*m*(m*m+2), 3)
    match = "OK" if tm_val == formula_val else f"FAIL ({tm_val} vs {formula_val})"
    print(f"  m={m}: {match}")

# ============================================================
# PATTERN 3: Try to find TM g_k formula for k=5
# ============================================================

print("\n" + "="*60)
print("PATTERN 3: TM g_5(m) formula search")
print("="*60)

# TM g_4 = m^2(m^2+2)/3. Clean!
# Can we find TM g_5 similarly?
# TM g_5(m): 1, 10, 51, 180, 501, 1182, 2471, 4712, 8361, 14002, ...
# These are: [x^5 in transfer matrix] / 2

# Check: is TM g_5(m) = polynomial of degree 5 in m?
# 5th differences should be constant = 16
tm5 = transfer_gk_values(5, 15)

vals = [tm5[m] for m in range(1, 13)]
print("TM g_5 values:", vals)

# Try to find the degree-5 polynomial
# g_5(m) = am^5 + bm^4 + cm^3 + dm^2 + em + f
# Use 6 values to determine 6 coefficients

from fractions import Fraction

# Set up system for m=1..6
# Actually, let's use the fact that the 5th difference is 16
# and work backwards
m_vals = list(range(1, 12))
v = [tm5[m] for m in m_vals]

# Forward differences
diffs = [list(v)]
for d in range(10):
    prev = diffs[-1]
    diffs.append([prev[i+1] - prev[i] for i in range(len(prev)-1)])

print("Forward difference table:")
for d in range(6):
    print(f"  d{d}: {diffs[d][:8]}")

# 5th diff = 16 for all. So leading coeff = 16/(5!) = 16/120 = 2/15
# g_5(m) = (2/15)*m^5 + ...

# Let's compute the polynomial by solving:
# 2*g_5(m)/3 = (2/15)*m^5*3/2 ... no, let me use a matrix approach

# Actually, use Lagrange or Newton interpolation
# g_5(m) = sum_{j=0}^{5} delta^j * C(m-1, j) where delta^j = j-th forward diff at m=1

# Newton forward differences at m=1:
newton_coeffs = [diffs[d][0] for d in range(6)]
print(f"\nNewton basis (at m=1): {newton_coeffs}")

# So g_5(m) = sum_{j=0}^5 newton_coeffs[j] * C(m-1, j)
# = 1 + 9*C(m-1,1) + 22*C(m-1,2) + 16*C(m-1,3) + 0*C(m-1,4) + 16*C(m-1,5)
# Wait: newton_coeffs[4] should be the 4th forward diff at m=1

print("\nVerify g_5(m) = 1 + 9*C(m-1,1) + 22*C(m-1,2) + ...")
print("  " + " + ".join(f"{newton_coeffs[j]}*C(m-1,{j})" for j in range(6) if newton_coeffs[j] != 0))

for m in range(1, 12):
    pred = sum(newton_coeffs[j] * comb(m-1, j) for j in range(6))
    actual = tm5[m]
    match = "OK" if pred == actual else f"FAIL"
    if m <= 8:
        print(f"  m={m}: pred={pred}, actual={actual} {match}")

# Convert to standard polynomial basis
# C(m-1,0) = 1
# C(m-1,1) = m-1
# C(m-1,2) = (m-1)(m-2)/2
# C(m-1,3) = (m-1)(m-2)(m-3)/6
# C(m-1,4) = (m-1)(m-2)(m-3)(m-4)/24
# C(m-1,5) = (m-1)(m-2)(m-3)(m-4)(m-5)/120

# 3*g_5(m) = 3 + 27*(m-1) + 33*(m-1)(m-2) + 8*(m-1)(m-2)(m-3) + 0 + 8/5*(m-1)(m-2)(m-3)(m-4)(m-5)
# Hmm, not integer. Let me check if 3*g_5 works.

# Actually let me find 15*g_5:
print("\n15*g_5(m):")
for m in range(1, 10):
    print(f"  m={m}: 15*g_5 = {15*tm5[m]}")

# 15*g_5: 15, 150, 765, 2700, 7515, 17730, 37065, 70680, 125415
# Check: 2*m^5 + lower
check = [2*m**5 for m in range(1, 10)]
print(f"  2*m^5:  {check}")
diff_from_m5 = [15*int(tm5[m]) - 2*m**5 for m in range(1, 10)]
print(f"  15*g5 - 2*m^5: {diff_from_m5}")

# 15*g5 - 2*m^5: 13, 86, 279, 620, 1265, 2502, 4649, 8056, 13107
# Is this a degree-4 polynomial?
d1 = [diff_from_m5[i+1] - diff_from_m5[i] for i in range(len(diff_from_m5)-1)]
d2 = [d1[i+1]-d1[i] for i in range(len(d1)-1)]
d3 = [d2[i+1]-d2[i] for i in range(len(d2)-1)]
d4 = [d3[i+1]-d3[i] for i in range(len(d3)-1)]
print(f"  d1: {d1}")
print(f"  d2: {d2}")
print(f"  d3: {d3}")
print(f"  d4: {d4}")

# So 15*g_5 = 2*m^5 + degree-4 polynomial
# Let's find: 15*g_5 = 2*m^5 + am^4 + bm^3 + cm^2 + dm + e
# From the differences: d4 constant -> degree 4 poly confirmed
# d4 all = 48? Let me check
# d3: [100, 148, 196, 244, 292] -> d4: [48, 48, 48, 48]
# Leading coeff of degree-4 part: 48/24 = 2. So 15*g5 = 2m^5 + 2m^4 + ...

r2 = [15*int(tm5[m]) - 2*m**5 - 2*m**4 for m in range(1, 10)]
print(f"\n  15*g5 - 2m^5 - 2m^4: {r2}")
# Check degree 3
d1r = [r2[i+1]-r2[i] for i in range(len(r2)-1)]
d2r = [d1r[i+1]-d1r[i] for i in range(len(d1r)-1)]
d3r = [d2r[i+1]-d2r[i] for i in range(len(d2r)-1)]
print(f"  d1: {d1r}")
print(f"  d2: {d2r}")
print(f"  d3: {d3r}")
# d3 constant -> degree 3. d3 = ?

# So 15*g_5 = 2*m^5 + 2*m^4 + am^3 + bm^2 + cm + d
r3 = [r2[i] for i in range(len(r2))]
print(f"  residual after m^5,m^4: {r3}")

# d3 is constant = -24? Let me see
# r3: 11, 54, 135, 260, 445, 706, 1059, 1520, 2105
# d1: 43, 81, 125, 185, 261, 353, 461, 585
# d2: 38, 44, 60, 76, 92, 108, 124  -- NOT constant at d2!
# Hmm wait let me recompute
d1_r3 = [r3[i+1]-r3[i] for i in range(len(r3)-1)]
d2_r3 = [d1_r3[i+1]-d1_r3[i] for i in range(len(d1_r3)-1)]
d3_r3 = [d2_r3[i+1]-d2_r3[i] for i in range(len(d2_r3)-1)]
print(f"  d1: {d1_r3}")
print(f"  d2: {d2_r3}")
print(f"  d3: {d3_r3}")

# Let me just solve for the polynomial directly
# 15*g_5(m) = 2m^5 + am^4 + bm^3 + cm^2 + dm + e
# m=1: 15 = 2 + a + b + c + d + e
# m=2: 150 = 64 + 16a + 8b + 4c + 2d + e
# m=3: 765 = 486 + 81a + 27b + 9c + 3d + e
# m=4: 2700 = 2048 + 256a + 64b + 16c + 4d + e
# m=5: 7515 = 6250 + 625a + 125b + 25c + 5d + e

# From m=1: a+b+c+d+e = 13
# From m=2: 16a+8b+4c+2d+e = 86
# From m=3: 81a+27b+9c+3d+e = 279
# From m=4: 256a+64b+16c+4d+e = 652
# From m=5: 625a+125b+25c+5d+e = 1265

# Subtract consecutive:
# 15a+7b+3c+d = 73
# 65a+19b+5c+d = 193
# 175a+37b+7c+d = 373
# 369a+61b+9c+d = 613

# Subtract again:
# 50a+12b+2c = 120
# 110a+18b+2c = 180
# 194a+24b+2c = 240

# Subtract:
# 60a+6b = 60 => 10a+b = 10
# 84a+6b = 60 => 14a+b = 10

# 14a+b - (10a+b) = 0 => 4a = 0 => a = 0! Wait...
# 10a+b=10, a=0, b=10. Hmm but we had 2m^4 earlier from d4=48.
# Let me recheck: 15*g_5(m) at m=1..5: 15, 150, 765, 2700, 7515
# minus 2m^5: 15-2=13, 150-64=86, 765-486=279, 2700-2048=652, 7515-6250=1265

# I said 2m^4 coefficient from d4=48, but actually the d4 of (15g5-2m^5) is:
# r2 = [13, 86, 279, 652, 1265, 2178, 3449, 5136, 7297]... wait I need to recompute
r2_exact = [15*int(tm5[m]) - 2*m**5 for m in range(1, 10)]
print(f"\nExact 15*g5 - 2m^5: {r2_exact}")

# OK so: [13, 86, 279, 620, 1265, 2502, 4649, 8056, 13107]
# d1: [73, 193, 341, 645, 1237, 2147, 3407, 5051]... let me compute
d1e = [r2_exact[i+1]-r2_exact[i] for i in range(len(r2_exact)-1)]
d2e = [d1e[i+1]-d1e[i] for i in range(len(d1e)-1)]
d3e = [d2e[i+1]-d2e[i] for i in range(len(d2e)-1)]
d4e = [d3e[i+1]-d3e[i] for i in range(len(d3e)-1)]
print(f"d1: {d1e}")
print(f"d2: {d2e}")
print(f"d3: {d3e}")
print(f"d4: {d4e}")

# From the equations: a+b+c+d+e=13, with a=0, b=10:
# c+d+e = 3. And 15*0+7*10+3c+d = 73 => 3c+d = 3. And 50*0+12*10+2c = 120 => 2c = 0 => c=0.
# So d=3, e=0. Check: 15*g_5(m) = 2m^5 + 10m^3 + 3m^2? No wait, we have degree-4 issues.

# Let me just solve with all 5 equations using exact arithmetic
from fractions import Fraction as F
# 15*g_5(m) = sum_{j=0}^5 c_j * m^j where c_5 = 2

# Vandermonde system for m=1..5:
# [1^j for j=0..4] * [c_0,...,c_4] = [15*g5(m) - 2*m^5 for m=1..5]

rhs = [F(15*int(tm5[m]) - 2*m**5) for m in range(1, 6)]
# Matrix: m^j for j=0..4, m=1..5
import numpy as np

# Solve with fractions
mat = [[F(m**j) for j in range(5)] for m in range(1, 6)]
# Gaussian elimination
for col in range(5):
    # Find pivot
    for row in range(col, 5):
        if mat[row][col] != 0:
            mat[col], mat[row] = mat[row], mat[col]
            rhs[col], rhs[row] = rhs[row], rhs[col]
            break
    pivot = mat[col][col]
    for j in range(5):
        mat[col][j] /= pivot
    rhs[col] /= pivot
    for row in range(5):
        if row != col and mat[row][col] != 0:
            factor = mat[row][col]
            for j in range(5):
                mat[row][j] -= factor * mat[col][j]
            rhs[row] -= factor * rhs[col]

print(f"\n15*g_5(m) = 2*m^5 + {rhs[4]}*m^4 + {rhs[3]}*m^3 + {rhs[2]}*m^2 + {rhs[1]}*m + {rhs[0]}")

# Verify
for m in range(1, 10):
    pred = 2*m**5 + sum(rhs[j]*m**j for j in range(5))
    actual = F(15*int(tm5[m]))
    match = "OK" if pred == actual else f"FAIL ({pred} vs {actual})"
    if m <= 6:
        print(f"  m={m}: {match}")

# ============================================================
# PATTERN 4: All TM g_k formulas
# ============================================================
print("\n" + "="*60)
print("PATTERN 4: Exact TM g_k polynomial formulas")
print("="*60)

for k in range(1, 8):
    tm = transfer_gk_values(k, k + 6)

    # Solve for polynomial of degree k
    # D * g_k(m) = sum_{j=0}^k c_j * m^j
    # We need to find the right common denominator D

    # Use Lagrange interpolation with k+1 points
    points = [(m, tm[m]) for m in range(1, k + 2)]
    if len(points) < k + 1:
        continue

    # Solve Vandermonde
    n_pts = k + 1
    ms = [p[0] for p in points]
    vs = [p[1] for p in points]

    mat = [[Fraction(m**j) for j in range(n_pts)] for m in ms]
    rhs = list(vs)

    for col in range(n_pts):
        for row in range(col, n_pts):
            if mat[row][col] != 0:
                mat[col], mat[row] = mat[row], mat[col]
                rhs[col], rhs[row] = rhs[row], rhs[col]
                break
        pivot = mat[col][col]
        for j in range(n_pts):
            mat[col][j] /= pivot
        rhs[col] /= pivot
        for row in range(n_pts):
            if row != col and mat[row][col] != 0:
                factor = mat[row][col]
                for j in range(n_pts):
                    mat[row][j] -= factor * mat[col][j]
                rhs[row] -= factor * rhs[col]

    # Find common denominator
    from math import gcd
    denoms = [c.denominator for c in rhs if c != 0]
    lcm = denoms[0]
    for d in denoms[1:]:
        lcm = lcm * d // gcd(lcm, d)

    scaled = [int(c * lcm) for c in rhs]
    print(f"\nk={k}: {lcm}*g_{k}(m) = ", end="")
    terms = []
    for j in range(n_pts-1, -1, -1):
        if scaled[j] != 0:
            terms.append(f"{scaled[j]}*m^{j}" if j > 0 else f"{scaled[j]}")
    print(" + ".join(terms))

    # Verify at more points
    ok = True
    for m in range(1, min(k + 6, 12)):
        if m in tm:
            pred = sum(rhs[j] * Fraction(m**j) for j in range(n_pts))
            if pred != tm[m]:
                ok = False
                print(f"  FAIL at m={m}: pred={pred}, actual={tm[m]}")
    if ok:
        print(f"  Verified at m=1..{min(k+5, 11)}")

print("\nDone!")
