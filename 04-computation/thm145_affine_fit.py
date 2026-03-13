#!/usr/bin/env python3
"""
thm145_affine_fit.py — opus-2026-03-12-S69

At p=13: 6 data points, 5 ESF variables.
Affine fit: Ω_m = c_0 + c_1*e_2 + c_2*e_3 + c_3*e_4 + c_4*e_5 + c_5*e_6
has 6 unknowns for 6 equations — always solvable.

Question: Are the coefficients rational with small denominators?
If so, we have an explicit formula.
"""

import numpy as np
from fractions import Fraction

# p=13 data
data_13 = [
    ((70, 84, 45, 11, 1),     [1,6,30,140,556,1823,4983,11241,19730,24780,20344,9642,1989]),
    ((122, 305, 357, 180, 27), [1,6,30,135,517,1622,4129,8436,13292,15078,11316,4988,977]),
    ((148, 448, 604, 362, 79), [1,6,30,135,520,1627,3915,6965,9012,8325,5157,1906,321]),
    ((161, 552, 812, 375, 53), [1,6,30,135,524,1675,4176,7707,10295,9681,6081,2298,399]),
    ((161, 552, 851, 570, 131),[1,6,30,135,524,1663,4001,6876,8409,7266,4237,1493,238]),
    ((174, 721, 1566, 1701, 729),[1,6,30,135,528,1707,4245,7503,9177,7770,4374,1488,237]),
]

# Build affine design matrix [1, e_2, e_3, e_4, e_5, e_6]
A = np.array([[1] + list(d[0]) for d in data_13], dtype=float)

print("="*80)
print("AFFINE FIT: Ω_m = c_0 + c_1*e_2 + c_2*e_3 + c_3*e_4 + c_4*e_5 + c_5*e_6")
print("="*80)

for mm in range(13):
    omega_vals = np.array([d[1][mm] for d in data_13], dtype=float)
    if len(set(omega_vals)) == 1:
        print(f"\n  Ω_{mm} = {int(omega_vals[0])} (universal)")
        continue

    x = np.linalg.solve(A, omega_vals)
    pred = A @ x
    err = max(abs(pred - omega_vals))

    print(f"\n  Ω_{mm}: max_err = {err:.2e}")
    print(f"    c_0 = {x[0]:.6f}")
    for j in range(5):
        print(f"    c_{j+1} (e_{j+2}) = {x[j+1]:.6f}")

    # Try to rationalize
    # Use continued fraction approach
    for denom_limit in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                        26, 39, 52, 78, 132, 143, 156, 286, 429, 858, 1716]:
        xi = [round(c * denom_limit) for c in x]
        pred2 = np.array([sum(xi[j]*A[i,j] for j in range(6)) for i in range(6)]) / denom_limit
        err2 = max(abs(pred2 - omega_vals))
        if err2 < 0.01:
            # Express as fraction
            fracs = [Fraction(xi[j], denom_limit) for j in range(6)]
            # Find common denominator
            from math import lcm
            common_d = 1
            for f in fracs:
                common_d = lcm(common_d, f.denominator)
            nums = [int(f * common_d) for f in fracs]
            print(f"    RATIONAL! denominator = {common_d}")
            print(f"    Ω_{mm} = ({nums[0]} + {nums[1]}*e_2 + {nums[2]}*e_3 + {nums[3]}*e_4 + {nums[4]}*e_5 + {nums[5]}*e_6) / {common_d}")
            break

# === Now try: does Ω_m depend on POWER SUMS rather than ESFs? ===
print(f"\n{'='*80}")
print(f"POWER SUM BASIS: p_k = Σ Q_i^k")
print(f"{'='*80}")

# Compute power sums from e_j using Newton's identities
# For p=13, m=6, e_1=21 (universal)
m = 6
e1 = m*(13-m)//2  # = 21

def esf_to_power_sums(esf_nonuniv, m=6, e1=21):
    """Convert (e_2,...,e_m) to power sums (p_1,...,p_m)."""
    e = [1, e1] + list(esf_nonuniv)
    pk = [m]  # p_0 = m
    # Newton: p_k = Σ_{j=1}^{k-1} (-1)^{j-1} e_j p_{k-j} + (-1)^{k-1} k e_k
    for k in range(1, m+1):
        val = 0
        for j in range(1, k):
            val += (-1)**(j-1) * e[j] * pk[k-j]
        val += (-1)**(k-1) * k * e[k]
        pk.append(val)
    return pk[1:]  # p_1,...,p_m

for d in data_13:
    pk = esf_to_power_sums(d[0])
    print(f"  e_j={d[0]}: p_k = {pk}")

# Affine fit using power sums
P = np.zeros((6, 7))  # 6 data points, 7 unknowns (const + p_1...p_6)
for i, d in enumerate(data_13):
    pk = esf_to_power_sums(d[0])
    P[i, 0] = 1
    P[i, 1:] = pk

# This is 7 unknowns for 6 equations — underdetermined
# Try just p_2,...,p_6 (p_1 is universal)
P5 = np.zeros((6, 6))
for i, d in enumerate(data_13):
    pk = esf_to_power_sums(d[0])
    P5[i, 0] = 1
    P5[i, 1:] = pk[1:]  # p_2,...,p_6

print(f"\n  Power sums: p_1 = {e1} (universal)")
print(f"\n  Affine fit in power sums (const + p_2...p_6):")

for mm in range(3, 13):
    omega_vals = np.array([d[1][mm] for d in data_13], dtype=float)
    x = np.linalg.solve(P5, omega_vals)
    pred = P5 @ x
    err = max(abs(pred - omega_vals))

    # Try to rationalize
    best_denom = None
    for denom_limit in range(1, 200):
        xi = [round(c * denom_limit) for c in x]
        pred2 = np.array([sum(xi[j]*P5[i,j] for j in range(6)) for i in range(6)]) / denom_limit
        err2 = max(abs(pred2 - omega_vals))
        if err2 < 0.01:
            best_denom = denom_limit
            break

    if best_denom and best_denom <= 100:
        xi = [round(c * best_denom) for c in x]
        print(f"  Ω_{mm}: denom = {best_denom}, coeffs = {xi}")
    else:
        print(f"  Ω_{mm}: no small rational form found")

# === LOOK AT chi_per formula ===
print(f"\n{'='*80}")
print(f"chi_per = Σ(-1)^m Ω_m as function of e_j")
print(f"{'='*80}")

chi_vals = np.array([sum((-1)**mm * d[1][mm] for mm in range(13)) for d in data_13], dtype=float)
print(f"  chi_per values: {[int(c) for c in chi_vals]}")

# Affine fit for chi_per
x_chi = np.linalg.solve(A, chi_vals)
pred_chi = A @ x_chi
err_chi = max(abs(pred_chi - chi_vals))
print(f"  Affine fit: max_err = {err_chi:.2e}")
print(f"    c_0 = {x_chi[0]:.6f}")
for j in range(5):
    print(f"    c_{j+1} (e_{j+2}) = {x_chi[j+1]:.6f}")

# Try rational
for denom_limit in range(1, 500):
    xi = [round(c * denom_limit) for c in x_chi]
    pred2 = np.array([sum(xi[j]*A[i,j] for j in range(6)) for i in range(6)]) / denom_limit
    err2 = max(abs(pred2 - chi_vals))
    if err2 < 0.01:
        fracs = [Fraction(xi[j], denom_limit) for j in range(6)]
        from math import lcm
        common_d = 1
        for f in fracs:
            common_d = lcm(common_d, f.denominator)
        nums = [int(f * common_d) for f in fracs]
        print(f"  RATIONAL chi_per formula (denom = {common_d}):")
        print(f"    chi_per = ({nums[0]} + {nums[1]}*e_2 + {nums[2]}*e_3 + {nums[3]}*e_4 + {nums[4]}*e_5 + {nums[5]}*e_6) / {common_d}")
        break

print("\nDONE.")
