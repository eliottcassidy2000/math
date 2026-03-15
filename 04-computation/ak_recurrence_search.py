#!/usr/bin/env python3
"""Search for recurrence relations for the a_k sequence.
a_k = leading coeff of 3*g_k(m), for k >= 1.
Known: 0, 0, 2, 10, 388, 69660, 19826270, 7309726742, 3262687720240, 1707771898925208, 1030345293948358170
"""
from fractions import Fraction

# a_k for k=1..11
a = [None, 0, 0, 2, 10, 388, 69660, 19826270, 7309726742, 3262687720240, 1707771898925208, 1030345293948358170]

# g_k(0) for k=1..11
g0 = [None, 0, 0, 0, -8, -592, -114320, -33338240, -12475185560, -5629549881808, -2973062116837472, -1807627677927630720]

# 3*g_k coefficients (a, b, c, d) for k=1..11
coeff = [None,
    (0, 0, 3, 0),
    (0, 3, 0, 0),
    (2, 0, 1, 0),
    (10, -33, 50, -24),
    (388, -2040, 3431, -1776),
    (69660, -380445, 653748, -342960),
    (19826270, -109486152, 189674605, -100014720),
    (7309726742, -40641958545, 70757788486, -37425556680),
    (3262687720240, -18232387983408, 31858349908595, -16888649645424),
    (1707771898925208, -9582908872031805, 16794323323619016, -8919186350512416),
    (1030345293948358170, -5802477398736520560, 10195015138571054553, -5422883033782892160),
]

print("=== Ratio analysis ===")
for k in range(4, 12):
    if a[k-1] != 0:
        r = Fraction(a[k], a[k-1])
        print(f"a_{k}/a_{k-1} = {float(r):.6f}, /k^2 = {float(r)/k**2:.6f}")

print("\n=== Looking for 2-term recurrence: a_k = P(k) * a_{k-1} + Q(k) * a_{k-2} ===")
# Try P(k) = α*k^2 + β*k + γ, Q(k) = δ*k^2 + ε*k + ζ
# Use k=5,6,7,8,9,10 to set up system (6 eqns, 6 unknowns)

from fractions import Fraction as F

# a_k = (α k² + β k + γ) a_{k-1} + (δ k² + ε k + ζ) a_{k-2}
# For k=5..10, this gives 6 equations in 6 unknowns
rows = []
for k in range(5, 11):
    # a_k = α k² a_{k-1} + β k a_{k-1} + γ a_{k-1} + δ k² a_{k-2} + ε k a_{k-2} + ζ a_{k-2}
    row = [F(k**2 * a[k-1]), F(k * a[k-1]), F(a[k-1]),
           F(k**2 * a[k-2]), F(k * a[k-2]), F(a[k-2])]
    rhs = F(a[k])
    rows.append((row, rhs))

# Gaussian elimination
import copy
mat = [row + [rhs] for row, rhs in rows]
n = 6
for col in range(n):
    # Find pivot
    pivot = None
    for r in range(col, n):
        if mat[r][col] != 0:
            pivot = r
            break
    if pivot is None:
        print(f"  No pivot at column {col}")
        continue
    mat[col], mat[pivot] = mat[pivot], mat[col]
    for r in range(n):
        if r != col and mat[r][col] != 0:
            factor = F(mat[r][col], mat[col][col])
            for c in range(n + 1):
                mat[r][c] -= factor * mat[col][c]

sol = [F(mat[i][n], mat[i][i]) for i in range(n)]
alpha, beta, gamma, delta, epsilon, zeta = sol
print(f"  α={float(alpha):.6f}, β={float(beta):.6f}, γ={float(gamma):.6f}")
print(f"  δ={float(delta):.6f}, ε={float(epsilon):.6f}, ζ={float(zeta):.6f}")

# Verify on k=11
pred_11 = (alpha * 11**2 + beta * 11 + gamma) * a[10] + (delta * 11**2 + epsilon * 11 + zeta) * a[9]
print(f"  Predicted a_11 = {pred_11}")
print(f"  Actual a_11    = {a[11]}")
print(f"  Match: {pred_11 == a[11]}")

# Also try with a_k starting from k=3 (first nonzero)
print("\n=== Same but using k=4..9 ===")
rows2 = []
for k in range(4, 10):
    row = [F(k**2 * a[k-1]), F(k * a[k-1]), F(a[k-1]),
           F(k**2 * a[k-2]), F(k * a[k-2]), F(a[k-2])]
    rhs = F(a[k])
    rows2.append((row, rhs))

mat2 = [row + [rhs] for row, rhs in rows2]
for col in range(n):
    pivot = None
    for r in range(col, n):
        if mat2[r][col] != 0:
            pivot = r
            break
    if pivot is None:
        print(f"  No pivot at column {col}")
        continue
    mat2[col], mat2[pivot] = mat2[pivot], mat2[col]
    for r in range(n):
        if r != col and mat2[r][col] != 0:
            factor = F(mat2[r][col], mat2[col][col])
            for c in range(n + 1):
                mat2[r][c] -= factor * mat2[col][c]

sol2 = [F(mat2[i][n], mat2[i][i]) for i in range(n)]
alpha2, beta2, gamma2, delta2, epsilon2, zeta2 = sol2
print(f"  α={float(alpha2):.6f}, β={float(beta2):.6f}, γ={float(gamma2):.6f}")
print(f"  δ={float(delta2):.6f}, ε={float(epsilon2):.6f}, ζ={float(zeta2):.6f}")

# Verify on k=10, 11
for check_k in [10, 11]:
    pred = (alpha2 * check_k**2 + beta2 * check_k + gamma2) * a[check_k-1] + (delta2 * check_k**2 + epsilon2 * check_k + zeta2) * a[check_k-2]
    print(f"  Predicted a_{check_k} = {pred}")
    print(f"  Actual a_{check_k}    = {a[check_k]}")
    print(f"  Match: {pred == a[check_k]}")

# Try 3-term recurrence: a_k = P(k) a_{k-1} + Q(k) a_{k-2} + R(k) a_{k-3}
# with P, Q, R all linear in k
print("\n=== 3-term recurrence, linear coeffs: a_k = (αk+β)a_{k-1} + (γk+δ)a_{k-2} + (εk+ζ)a_{k-3} ===")
rows3 = []
for k in range(6, 12):
    row = [F(k * a[k-1]), F(a[k-1]),
           F(k * a[k-2]), F(a[k-2]),
           F(k * a[k-3]), F(a[k-3])]
    rhs = F(a[k])
    rows3.append((row, rhs))

mat3 = [row + [rhs] for row, rhs in rows3]
for col in range(6):
    pivot = None
    for r in range(col, 6):
        if mat3[r][col] != 0:
            pivot = r
            break
    if pivot is None:
        print(f"  No pivot at column {col}")
        continue
    mat3[col], mat3[pivot] = mat3[pivot], mat3[col]
    for r in range(6):
        if r != col and mat3[r][col] != 0:
            factor = F(mat3[r][col], mat3[col][col])
            for c in range(7):
                mat3[r][c] -= factor * mat3[col][c]

sol3 = [F(mat3[i][7-1], mat3[i][i]) for i in range(6)]
print(f"  Coeffs: {[float(s) for s in sol3]}")

# Check if any are rational with small denominators
for i, s in enumerate(sol3):
    print(f"  sol[{i}] = {s} ≈ {float(s):.10f}")

# Now try the b_k sequence (second coeff of 3*g_k)
print("\n=== b_k sequence analysis ===")
b = [None] + [c[1] for c in coeff[1:]]
print(f"b_k = {b[1:]}")
print("\nb_k / a_k ratios:")
for k in range(4, 12):
    if a[k] != 0:
        print(f"  b_{k}/a_{k} = {Fraction(b[k], a[k])} ≈ {b[k]/a[k]:.6f}")

# Check if b_k has a simple relation to a_k
# From THM-217: D1 = 1 - D0, D2 = D0 + 2(k-1) where D0=g_k(0), D3=2*a_k
# 3*g_k(m) = a_k*m^3 + b_k*m^2 + c_k*m + d_k
# g_k(0) = d_k/3
# g_k(1) = 1 => a_k + b_k + c_k + d_k = 3
# g_k(2) = 2k => 8a_k + 4b_k + 2c_k + d_k = 6k
# From these:
# (1) a + b + c + d = 3
# (2) 8a + 4b + 2c + d = 6k
# (2)-(1): 7a + 3b + c = 6k - 3
# Also d = 3*g_k(0)
# So c = 3 - a - b - d = 3 - a - b - 3*g_k(0)
# And 7a + 3b + (3 - a - b - 3*g_k(0)) = 6k - 3
# 6a + 2b - 3*g_k(0) = 6k - 6
# b = (6k - 6 - 6a + 3*g_k(0)) / 2 = 3k - 3 - 3a + 3*g_k(0)/2

print("\n=== Verify b_k formula from constraints ===")
for k in range(1, 12):
    ak, bk, ck, dk = coeff[k]
    predicted_b = 3*k - 3 - 3*ak + F(3*g0[k], 2)  # 3*g_k(0)/2 = dk/2
    # Actually dk = 3*g_k(0), so dk/2 = 3*g_k(0)/2
    actual_b = bk
    print(f"  k={k}: predicted b = {predicted_b}, actual = {actual_b}, match = {predicted_b == actual_b}")

print("\n=== So b_k is determined by a_k and g_k(0). Checking g_k(0)/a_k ratio ===")
for k in range(4, 12):
    if a[k] != 0:
        print(f"  g_{k}(0)/a_{k} = {Fraction(g0[k], a[k])} ≈ {g0[k]/a[k]:.10f}")

# Check c_k and d_k patterns
print("\n=== c_k sequence ===")
c_seq = [coeff[k][2] for k in range(1, 12)]
print(f"c_k = {c_seq}")

print("\n=== d_k = 3*g_k(0) sequence ===")
d_seq = [coeff[k][3] for k in range(1, 12)]
print(f"d_k = {d_seq}")

# Let's look at a_k from a different angle: EGF or OGF coefficients
# If W(n) has an EGF, the g_k structure might come from coefficient extraction
print("\n=== W(n)/n! analysis ===")
import math
wn = [None, 1, 2, 8, 32, 158, 928, 6350, 49752, 439670, 4327904, 46963358, 556953448,
      7166360054, 99428495088, 1479600188798, 23506712352248, 397095175477430,
      7107209383674112, 134345623603516190, 2674426516381764744,
      55925062706620208438, 1225582324497129993488, 28088043261491650347134,
      671901551280362054066968, 16746599265666151628174198,
      434187457363955400414175008, 11692423738081050318010736030]

print("W(n)/n! - 1:")
for n in range(3, 28):
    cv2 = Fraction(wn[n], math.factorial(n)) - 1
    print(f"  n={n}: {float(cv2):.15f}")

# Look at (W(n)/n! - 1) * n to see convergence to 2
print("\nCV^2 * n:")
for n in range(3, 28):
    cv2 = Fraction(wn[n], math.factorial(n)) - 1
    print(f"  n={n}: {float(cv2 * n):.10f}")

# Look at (CV^2 * n - 2) * n to find next correction term
print("\n(CV^2 * n - 2) * n (should converge to something):")
for n in range(5, 28):
    cv2 = Fraction(wn[n], math.factorial(n)) - 1
    val = (cv2 * n - 2) * n
    print(f"  n={n}: {float(val):.10f}")
