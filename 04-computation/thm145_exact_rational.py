#!/usr/bin/env python3
"""
thm145_exact_rational.py — opus-2026-03-12-S69

Use exact rational arithmetic (fractions) to find the affine formula
Ω_m = c_0 + c_1*e_2 + ... + c_5*e_6 at p=13.

Since all values are integers, the coefficients MUST be rational.
"""

from fractions import Fraction

# p=13 data (all integers)
esfs = [
    (70, 84, 45, 11, 1),
    (122, 305, 357, 180, 27),
    (148, 448, 604, 362, 79),
    (161, 552, 812, 375, 53),
    (161, 552, 851, 570, 131),
    (174, 721, 1566, 1701, 729),
]

omegas = [
    [1,6,30,140,556,1823,4983,11241,19730,24780,20344,9642,1989],
    [1,6,30,135,517,1622,4129,8436,13292,15078,11316,4988,977],
    [1,6,30,135,520,1627,3915,6965,9012,8325,5157,1906,321],
    [1,6,30,135,524,1675,4176,7707,10295,9681,6081,2298,399],
    [1,6,30,135,524,1663,4001,6876,8409,7266,4237,1493,238],
    [1,6,30,135,528,1707,4245,7503,9177,7770,4374,1488,237],
]

def solve_exact(A, b):
    """Solve Ax=b exactly using Fraction arithmetic (Gaussian elimination)."""
    n = len(A)
    # Augmented matrix
    M = [[Fraction(A[i][j]) for j in range(n)] + [Fraction(b[i])] for i in range(n)]

    for col in range(n):
        # Find pivot
        pivot = None
        for row in range(col, n):
            if M[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            return None
        M[col], M[pivot] = M[pivot], M[col]

        # Eliminate
        for row in range(n):
            if row != col and M[row][col] != 0:
                factor = M[row][col] / M[col][col]
                for j in range(n+1):
                    M[row][j] -= factor * M[col][j]

    return [M[i][n] / M[i][i] for i in range(n)]

# Build design matrix [1, e_2, e_3, e_4, e_5, e_6]
A = [[1] + list(e) for e in esfs]

print("="*80)
print("EXACT RATIONAL AFFINE FORMULA: Ω_m = c_0 + Σ c_j * e_{j+1}")
print("="*80)

for mm in range(13):
    omega_vals = [omegas[i][mm] for i in range(6)]
    if len(set(omega_vals)) == 1:
        print(f"\nΩ_{mm} = {omega_vals[0]} (universal)")
        continue

    x = solve_exact(A, omega_vals)
    if x is None:
        print(f"\nΩ_{mm}: singular system!")
        continue

    # Verify
    for i in range(6):
        val = sum(x[j] * Fraction(A[i][j]) for j in range(6))
        assert val == omega_vals[i], f"Verification failed at row {i}: {val} != {omega_vals[i]}"

    print(f"\nΩ_{mm}:")
    labels = ['const', 'e_2', 'e_3', 'e_4', 'e_5', 'e_6']
    for j, label in enumerate(labels):
        print(f"  {label:6s} = {x[j]} = {float(x[j]):.6f}")

    # Find common denominator
    from math import lcm
    common_d = 1
    for f in x:
        common_d = lcm(common_d, f.denominator)
    nums = [int(f * common_d) for f in x]
    print(f"  Common denom = {common_d}")
    print(f"  Ω_{mm} = ({' + '.join(f'{nums[j]}*{labels[j]}' for j in range(6))}) / {common_d}")

# chi_per formula
print(f"\n{'='*80}")
print(f"chi_per EXACT FORMULA")
print(f"{'='*80}")

chi_vals = [sum((-1)**mm * omegas[i][mm] for mm in range(13)) for i in range(6)]
print(f"chi_per values: {chi_vals}")

x_chi = solve_exact(A, chi_vals)
from math import lcm
labels = ['const', 'e_2', 'e_3', 'e_4', 'e_5', 'e_6']
common_d = 1
for f in x_chi:
    common_d = lcm(common_d, f.denominator)
nums = [int(f * common_d) for f in x_chi]

print(f"\nchi_per = ({' + '.join(f'{nums[j]}*{labels[j]}' for j in range(6))}) / {common_d}")
for j, label in enumerate(labels):
    print(f"  {label:6s} = {x_chi[j]}")

# === Now verify this formula predicts the data ===
print(f"\n  Verification:")
for i in range(6):
    pred = sum(x_chi[j] * Fraction(A[i][j]) for j in range(6))
    print(f"  e_j={esfs[i]}: chi_per = {pred} (actual {chi_vals[i]}) {'✓' if pred == chi_vals[i] else '✗'}")

# === Look for patterns in the denominators ===
print(f"\n{'='*80}")
print(f"DENOMINATOR ANALYSIS")
print(f"{'='*80}")

for mm in range(3, 13):
    omega_vals = [omegas[i][mm] for i in range(6)]
    x = solve_exact(A, omega_vals)
    common_d = 1
    for f in x:
        common_d = lcm(common_d, f.denominator)
    print(f"  Ω_{mm}: common denominator = {common_d}, factored = ", end="")
    d = common_d
    factors = []
    for p in [2,3,5,7,11,13,17,19,23]:
        while d % p == 0:
            factors.append(p)
            d //= p
    if d > 1:
        factors.append(d)
    print(factors)

print(f"\n  chi_per: common denominator = {common_d}")

print("\nDONE.")
