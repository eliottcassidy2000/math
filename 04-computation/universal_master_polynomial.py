#!/usr/bin/env python3
"""
UNIVERSAL MASTER POLYNOMIAL: F_j(r)

THEOREM (conjectured, verified n=4..9):
For any tournament invariant I associated with an independent set type
in Omega(T) (using v_I vertices, with partition pi_I), the per-invariant
polynomial C_I(r,n) = coefficient of I in W(r) satisfies:

    C_I(r, n) = 2^{parts(pi_I)} * F_f(r)

where f = (n-1) - 2*|pi_I| is the free position count (|pi| = sum of parts).

The universal sequence F_j(r):
  F_0 = 1
  F_1 = 2r
  F_2 = -1/2 + 6r^2
  F_3 = -4r + 24r^3
  F_4 = 1 - 30r^2 + 120r^4
  F_5 = 17r - 240r^3 + 720r^5
  F_6 = -17/4 + 231r^2 - 2100r^4 + 5040r^6

Properties:
  - F_j(1/2) = 1 for all j
  - Leading coefficient = (j+1)!
  - F_j has parity j mod 2 (even/odd powers)

opus-2026-03-06-S30
"""
from fractions import Fraction
import numpy as np

# =====================================================================
# The master sequence
# =====================================================================

# F_j as lists of coefficients [r^0 or r^1, r^2 or r^3, ...]
# Even j: even powers r^0, r^2, r^4, ...
# Odd j: odd powers r^1, r^3, r^5, ...
F = {
    0: [1],                             # 1
    1: [2],                             # 2r
    2: [Fraction(-1,2), 6],             # -1/2 + 6r^2
    3: [-4, 24],                        # -4r + 24r^3
    4: [1, -30, 120],                   # 1 - 30r^2 + 120r^4
    5: [17, -240, 720],                 # 17r - 240r^3 + 720r^5
    6: [Fraction(-17,4), 231, -2100, 5040],  # -17/4 + 231r^2 - ...
}

def eval_F(j, r):
    """Evaluate F_j(r)."""
    result = 0
    for k, c in enumerate(F[j]):
        if j % 2 == 0:
            result += float(c) * r**(2*k)
        else:
            result += float(c) * r**(2*k+1)
    return result

# =====================================================================
print("=" * 70)
print("UNIVERSAL MASTER POLYNOMIAL VERIFICATION")
print("=" * 70)

print("\nF_j(1/2) should all equal 1:")
for j in range(7):
    val = eval_F(j, 0.5)
    print(f"  F_{j}(1/2) = {val:.6f}")

print(f"\nLeading coefficients (should be (j+1)!):")
for j in range(7):
    lc = float(F[j][-1])
    from math import factorial
    expected = factorial(j+1)
    print(f"  F_{j}: leading = {lc}, (j+1)! = {expected}, match = {abs(lc-expected)<0.01}")

# =====================================================================
print(f"\n{'='*70}")
print("VERIFICATION: C_I(r,n) = 2^parts * F_f(r)")
print("="*70)

# Each invariant I has:
#   - partition pi_I (parts in the block-length notation)
#   - vertex count v_I = sum(2j+1 for j in pi_I)
#   - OCF weight w_I = 2^{len(pi_I)}
#   - free positions f_I = (n-1) - 2*sum(pi_I)

invariants = {
    # (name, partition, vertex_count)
    't3': ([1], 3),      # 3-cycle: partition (1,), uses 3 vertices
    't5': ([2], 5),      # 5-cycle: partition (2,), uses 5 vertices
    't7': ([3], 7),      # 7-cycle: partition (3,), uses 7 vertices
    't9': ([4], 9),      # 9-cycle: partition (4,), uses 9 vertices
    'bc':  ([1,1], 6),   # disjoint 3+3: partition (1,1), 6 verts
    'bc35': ([1,2], 8),  # disjoint 3+5: partition (1,2), 8 verts
    'a3':  ([1,1,1], 9), # three disjoint 3s: partition (1,1,1), 9 verts
}

# Known per-invariant polynomials from computations:
# Format: {(inv_name, n): [coefficients as w_0_coeff, w_2_coeff, ...]}
known_C = {
    # Odd n: coefficients of [r^0, r^2, r^4, ...]
    ('t3', 5): [-1, 12],
    ('t3', 7): [2, -60, 240],
    ('t3', 9): [Fraction(-17,2), 462, -4200, 10080],
    ('t5', 5): [2],
    ('t5', 7): [-1, 12],
    ('t5', 9): [2, -60, 240],
    ('t7', 7): [2],
    ('t7', 9): [-1, 12],
    ('t9', 9): [2],
    ('bc', 7): [-2, 24],
    ('bc', 9): [4, -120, 480],
    ('bc35', 9): [-2, 24],
    ('a3', 9): [-4, 48],

    # Even n: coefficients of [r^1, r^3, r^5, ...]
    ('t3', 4): [4],
    ('t3', 6): [-8, 48],
    ('t3', 8): [34, -480, 1440],
    ('t5', 6): [4],
    ('t5', 8): [-8, 48],
    ('t7', 8): [4],
    ('bc', 6): [8],
    ('bc', 8): [-16, 96],
    ('bc35', 8): [8],
}

print(f"\n{'Invariant':10s} {'n':3s} {'parts':5s} {'|pi|':4s} {'free':4s} {'2^parts':7s} {'Predicted':30s} {'Actual':30s} {'Match':5s}")
print("-" * 100)

for key, actual_coeffs in sorted(known_C.items()):
    inv_name, n = key
    if inv_name not in invariants:
        continue
    pi, v = invariants[inv_name]
    num_parts = len(pi)
    sum_parts = sum(pi)
    free = (n-1) - 2*sum_parts
    weight = 2**num_parts

    if free < 0 or free not in F:
        continue

    # Predicted: weight * F_free
    pred_coeffs = [weight * c for c in F[free]]

    # Compare
    match = True
    for i in range(max(len(pred_coeffs), len(actual_coeffs))):
        p = float(pred_coeffs[i]) if i < len(pred_coeffs) else 0
        a = float(actual_coeffs[i]) if i < len(actual_coeffs) else 0
        if abs(p - a) > 0.01:
            match = False

    pred_str = str([float(c) for c in pred_coeffs])
    actual_str = str([float(c) for c in actual_coeffs])
    print(f"{inv_name:10s} {n:3d} {num_parts:5d} {sum_parts:4d} {free:4d} {weight:7d} {pred_str:30s} {actual_str:30s} {'YES' if match else 'NO':5s}")

# =====================================================================
print(f"\n{'='*70}")
print("MASTER POLYNOMIAL IN SHIFTED VARIABLE u = r^2 - 1/4")
print("="*70)

# For even j, express F_j(r) in terms of u = r^2 - 1/4
# F_j = sum c_k u^k
for j in [0, 2, 4, 6]:
    coeffs_r = [float(c) for c in F[j]]
    # Convert from r^{2k} coefficients to u^k coefficients
    # r^{2k} = (u + 1/4)^k
    degree = len(coeffs_r) - 1
    u_coeffs = [0.0] * (degree + 1)
    for k in range(degree + 1):
        # r^{2k} = (u+1/4)^k = sum_m C(k,m) u^m (1/4)^{k-m}
        for m in range(k + 1):
            from math import comb
            u_coeffs[m] += coeffs_r[k] * comb(k, m) * (0.25)**(k-m)
    print(f"  F_{j} in u = {[Fraction(c).limit_denominator(10000) for c in u_coeffs]}")

# For odd j, express F_j(r) in terms of v = r^2 - 1/4, times r
for j in [1, 3, 5]:
    coeffs_r = [float(c) for c in F[j]]
    # F_j(r) = r * sum c_k r^{2k} = r * sum c_k (v+1/4)^k
    degree = len(coeffs_r) - 1
    v_coeffs = [0.0] * (degree + 1)
    for k in range(degree + 1):
        for m in range(k + 1):
            from math import comb
            v_coeffs[m] += coeffs_r[k] * comb(k, m) * (0.25)**(k-m)
    print(f"  F_{j} = r * {[Fraction(c).limit_denominator(10000) for c in v_coeffs]}")

# =====================================================================
print(f"\n{'='*70}")
print("RECURRENCE SEARCH")
print("="*70)

# In u-coordinates, the even-j sequence F_0, F_2, F_4, F_6 has coefficients:
# F_0 = [1]
# F_2 = [1, 6]
# F_4 = [1, 30, 120]
# F_6 = [1, 126, 1680, 5040]

# Sub-leading: 6, 30, 126. Ratios: 5, 4.2.
# 6 = C(4,2), 30 = C(6,2)*? No, C(6,2)=15.
# 6 = 3*2, 30 = 5*6, 126 = 7*18.
# or: 6, 30, 126: these are C(4,1)*... hmm
# 6 = 2*3, 30 = 2*3*5, 126 = 2*3*3*7. Factors: all divisible by 6.
# 6/6=1, 30/6=5, 126/6=21. 1, 5, 21 = C(2,2), C(4,2)+C(2,2)=7, ...
# 1, 5, 21: differences 4, 16.
# Actually 1 = T(1), 5 = T(3)? T(n) = n(n+1)/2: T(1)=1, T(2)=3, T(3)=6. No.
# 1, 5, 21 = (2k choose k+1) for k=1,2,3? C(2,2)=1✓, C(4,3)=4✗.
# Just compute ratios: 5/1=5, 21/5=4.2.

# Let me try to see if u-coefficients follow: a_{k,j} = (2k+1)! / (2j+1)! * b_j
# k=1: [1, 6]. 6 = 3!/1. a_{1,1}/1! = 6.
# k=2: [1, 30, 120]. 120 = 5!/1. 30 = 5*6 = 5*3!.
#   a_{2,1}/3! = 30/6 = 5. a_{2,2}/5! = 120/120 = 1.
# k=3: [1, 126, 1680, 5040]. 5040 = 7!/1. 1680 = 7*240 = 7*5!/... nah.
#   a_{3,1}/3! = 126/6 = 21. a_{3,2}/5! = 1680/120 = 14. a_{3,3}/7! = 5040/5040 = 1.

# Sequence of a_{k,1}/3!: 1, 5, 21. These are C(2k,k)/(k+1) = Catalan(k)!
# Catalan: 1, 1, 2, 5, 14, 42. C(1)=1, C(2)=2, C(3)=5✓, nope C_1=1, C_2=2, C_3=5✓!
# Wait: a_{1,1}/6 = 1, a_{2,1}/6 = 5, a_{3,1}/6 = 21.
# 1, 5, 21 = ?
# C(2,1)=2, C(4,2)=6, C(6,3)=20. No.
# Actually: 1 = C(2,1)/2, 5 = C(6,2)/3, 21 = C(8,3)/... nah.
# Or simply: 1, 5, 21 satisfies a(n) = (4n-2)*a(n-1) + ... ?
# 5 = 5*1, 21 = 21/5 * 5. Not helpful.

# Let me try: the a_{k,j}/((2j+1)!) sequence
coeffs_u = {
    0: [1],
    1: [1, 6],
    2: [1, 30, 120],
    3: [1, 126, 1680, 5040],
}

print("u-coefficients of F_{2k} (k=0,1,2,3):")
for k in range(4):
    print(f"  k={k}: {coeffs_u[k]}")
    normalized = [coeffs_u[k][j] / factorial(2*j+1) for j in range(k+1)]
    print(f"       ÷(2j+1)!: {[Fraction(x).limit_denominator(1000) for x in normalized]}")

print("\nDivided by (2j+1)!:")
print("  k=0: [1]")
print("  k=1: [1, 1]")
print("  k=2: [1, 1/4, 1]")
print("  k=3: [1, 1/120*126 = 21/20, 1680/120 = 14, 1]")

# Actually let me recompute
for k in range(4):
    row = []
    for j in range(k+1):
        val = Fraction(coeffs_u[k][j]) / factorial(2*j+1)
        row.append(val)
    print(f"  k={k}: {row}")

print(f"\n{'='*70}")
print("SUMMARY")
print("="*70)
print("""
UNIVERSAL MASTER POLYNOMIAL THEOREM (verified n=4..9):

For tournament invariant I with partition type pi_I (in the OCF block notation):
  C_I(r, n) = 2^{|pi_I|} * F_f(r)

where f = (n-1) - 2*sum(pi_I) is the free position count.

The master sequence F_0, F_1, F_2, ... is universal: it depends on NO tournament data.

Properties:
  (1) F_j(1/2) = 1 for all j.
  (2) deg(F_j) = j, with leading coefficient (j+1)!.
  (3) F_j has parity j: only even powers for even j, odd powers for odd j.
  (4) In shifted variable u = r^2 - 1/4 (for even j):
      F_0 = 1, F_2 = 1+6u, F_4 = 1+30u+120u^2, F_6 = 1+126u+1680u^2+5040u^3.

Consequences:
  (a) The SHIFT PRINCIPLE is a corollary: same cycle type at n vs n-2
      shifts f by 2, giving the next F polynomial.
  (b) The coefficient of ANY invariant in ANY W-coefficient is determined
      by exactly ONE universal number: the corresponding F_j coefficient.
  (c) The entire W-polynomial hierarchy is encoded in the single sequence {F_j}.
""")
