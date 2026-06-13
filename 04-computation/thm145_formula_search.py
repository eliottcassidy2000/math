#!/usr/bin/env python3
"""
thm145_formula_search.py — opus-2026-03-12-S69

Search for the explicit formula Ω_m = f(e_2,...,e_m).

Strategy: At p=13, we have 6 data points and 5 ESF variables.
Try polynomial fits of various degrees.

Also look for RATIOS and DIFFERENCES that simplify.
"""

import numpy as np
from numpy.linalg import lstsq

# Data from p=13 verification (sorted by e_2)
# Format: (e_2, e_3, e_4, e_5, e_6), Ω profile
data_13 = [
    ((70, 84, 45, 11, 1),     [1,6,30,140,556,1823,4983,11241,19730,24780,20344,9642,1989]),
    ((122, 305, 357, 180, 27), [1,6,30,135,517,1622,4129,8436,13292,15078,11316,4988,977]),
    ((148, 448, 604, 362, 79), [1,6,30,135,520,1627,3915,6965,9012,8325,5157,1906,321]),
    ((161, 552, 812, 375, 53), [1,6,30,135,524,1675,4176,7707,10295,9681,6081,2298,399]),
    ((161, 552, 851, 570, 131),[1,6,30,135,524,1663,4001,6876,8409,7266,4237,1493,238]),
    ((174, 721, 1566, 1701, 729),[1,6,30,135,528,1707,4245,7503,9177,7770,4374,1488,237]),
]

# Data from p=11
data_11 = [
    ((35,28,9,1),    [1,5,20,74,224,522,926,1222,1115,611,148]),
    ((68,127,97,23),  [1,5,20,70,200,439,711,827,648,301,64]),
    ((79,171,130,23), [1,5,20,70,201,430,620,596,384,151,26]),
    ((90,270,405,243),[1,5,20,70,205,460,700,690,450,180,30]),
]

# Data from p=7
data_7 = [
    ((5, 1),  [1,3,6,11,14,9,2]),
    ((12, 8), [1,3,6,9,9,6,3]),
]

print("="*80)
print("SEARCHING FOR Ω_m = f(e_j) FORMULA")
print("="*80)

# === p=7: 2 points, 2 variables ===
print("\n--- p=7: Linear fit (2 points, 2 unknowns) ---")
E = np.array([[d[0][0], d[0][1]] for d in data_7], dtype=float)
for mm in range(7):
    omega_vals = [d[1][mm] for d in data_7]
    if omega_vals[0] == omega_vals[1]:
        print(f"  Ω_{mm} = {omega_vals[0]} (universal)")
        continue
    # Ω_m = a + b*e_2 + c*e_3
    A = np.array([[1, data_7[0][0][0], data_7[0][0][1]],
                   [1, data_7[1][0][0], data_7[1][0][1]]], dtype=float)
    # Underdetermined: 2 equations, 3 unknowns. Try without constant:
    # Ω_m = b*e_2 + c*e_3
    A2 = np.array([[data_7[0][0][0], data_7[0][0][1]],
                    [data_7[1][0][0], data_7[1][0][1]]], dtype=float)
    b = np.array(omega_vals, dtype=float)
    x = np.linalg.solve(A2, b)
    print(f"  Ω_{mm} = {x[0]:.6f}*e_2 + {x[1]:.6f}*e_3 → {x[0]}*e_2 + {x[1]}*e_3")

# === p=11: 4 points, 4 variables ===
print("\n--- p=11: Linear fit Ω_m = a*e_2 + b*e_3 + c*e_4 + d*e_5 ---")
E11 = np.array([[d[0][0], d[0][1], d[0][2], d[0][3]] for d in data_11], dtype=float)
for mm in range(11):
    omega_vals = np.array([d[1][mm] for d in data_11], dtype=float)
    if len(set(omega_vals)) == 1:
        print(f"  Ω_{mm} = {int(omega_vals[0])} (universal)")
        continue
    x, res, rank, sv = lstsq(E11, omega_vals, rcond=None)
    pred = E11 @ x
    err = max(abs(pred - omega_vals))
    # Check if coefficients are rational with small denominators
    print(f"  Ω_{mm}: coeffs = [{', '.join(f'{c:.6f}' for c in x)}], max_err = {err:.6f}")
    # Check if integer or simple rational
    for denom in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 55]:
        xi = [round(c*denom) for c in x]
        pred2 = np.array([sum(xi[j]*E11[i,j] for j in range(4)) for i in range(4)]) / denom
        err2 = max(abs(pred2 - omega_vals))
        if err2 < 0.01:
            print(f"    → Ω_{mm} = ({' + '.join(f'{xi[j]}*e_{j+2}' for j in range(4))}) / {denom}")
            break

# === p=13: 6 points, 5 variables ===
print("\n--- p=13: Linear fit Ω_m = a*e_2 + b*e_3 + c*e_4 + d*e_5 + f*e_6 ---")
E13 = np.array([[d[0][0], d[0][1], d[0][2], d[0][3], d[0][4]] for d in data_13], dtype=float)
for mm in range(13):
    omega_vals = np.array([d[1][mm] for d in data_13], dtype=float)
    if len(set(omega_vals)) == 1:
        print(f"  Ω_{mm} = {int(omega_vals[0])} (universal)")
        continue
    x, res, rank, sv = lstsq(E13, omega_vals, rcond=None)
    pred = E13 @ x
    err = max(abs(pred - omega_vals))
    is_linear = err < 0.1
    print(f"  Ω_{mm}: max_err = {err:.4f} {'✓ LINEAR' if is_linear else '✗ NOT LINEAR'}")
    if is_linear:
        for denom in range(1, 100):
            xi = [round(c*denom) for c in x]
            pred2 = np.array([sum(xi[j]*E13[i,j] for j in range(5)) for i in range(6)]) / denom
            err2 = max(abs(pred2 - omega_vals))
            if err2 < 0.01:
                print(f"    → ({' + '.join(f'{xi[j]}*e_{j+2}' for j in range(5))}) / {denom}")
                break

# === Check if Ω_m - (universal value) depends on product e_m only ===
print("\n--- p=13: Does Ω_m depend mainly on e_6 (= prod Q)? ---")
for mm in range(3, 13):
    omega_vals = [d[1][mm] for d in data_13]
    e6_vals = [d[0][4] for d in data_13]
    if len(set(omega_vals)) <= 1:
        continue
    corr = np.corrcoef(omega_vals, e6_vals)[0,1]
    print(f"  Ω_{mm}: corr(Ω, e_6) = {corr:+.4f}")

# === Look at TOTAL Ω ===
print("\n--- Total Ω = Σ Ω_m ---")
for dataset, p_val in [(data_7, 7), (data_11, 11), (data_13, 13)]:
    print(f"\n  p={p_val}:")
    for d in dataset:
        total = sum(d[1])
        print(f"    e_j = {d[0]}: total Ω = {total}, chi_per = {sum((-1)**mm * w for mm, w in enumerate(d[1]))}")

# === Check Ω differences between orbits ===
print("\n--- p=13: Ω RATIOS between consecutive degrees ---")
for d in data_13:
    ratios = [d[1][mm+1]/d[1][mm] if d[1][mm] > 0 else 0 for mm in range(12)]
    print(f"  e_j={d[0]}: ratios = [{', '.join(f'{r:.3f}' for r in ratios)}]")

# === Newton's identities: power sums p_k = Σ Q_i^k ===
print("\n--- p=13: Power sums of Q vs Ω ---")
# At p=13, we can compute power sums from ESF via Newton's identities
# p_1 = e_1, p_2 = e_1*p_1 - 2*e_2, etc.
# But we already know e_1 = m(p-m)/2 = 6*7/2 = 21
e1 = 21  # universal at p=13
for d in data_13:
    e = (e1,) + d[0]  # (e_1, e_2, ..., e_6)
    # Newton: p_k = e_1*p_{k-1} - e_2*p_{k-2} + ... + (-1)^{k-1}*k*e_k
    pk = [0, e1]  # p_0=m=6, p_1=e_1=21
    pk[0] = 6  # p_0 = m
    for k in range(2, 7):
        val = sum((-1)**(j-1) * e[j] * pk[k-j] for j in range(1, k))
        val += (-1)**(k-1) * k * e[k]
        pk.append(val)
    print(f"  e_j={d[0]}: p_k = {pk[1:]}, chi_per = {sum((-1)**mm * w for mm, w in enumerate(d[1]))}")

print("\nDONE.")
