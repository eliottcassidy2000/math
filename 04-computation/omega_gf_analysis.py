"""
Omega generating function and product formula analysis.
opus-2026-03-13-S71b
"""
import numpy as np
from math import comb, factorial

omega_7 = [1, 3, 6, 9, 9, 6, 3]
omega_11 = [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]

print("O(-1) = chi/p:")
for name, omega in [("P_7", omega_7), ("P_11", omega_11)]:
    alt_sum = sum((-1)**d * omega[d] for d in range(len(omega)))
    print(f"  {name}: {alt_sum}")

print("\nRatios Omega_{d+1}/Omega_d:")
for name, omega, m, p in [("P_7", omega_7, 3, 7), ("P_11", omega_11, 5, 11)]:
    ratios = [omega[d+1]/omega[d] for d in range(len(omega)-1)]
    print(f"  {name}: {[f'{r:.4f}' for r in ratios]}")

print("\nA_d/p analysis:")
A_7 = [1, 3, 9, 21, 39, 45, 27]
A_11 = [1, 5, 25, 110, 430, 1430, 3970, 8735, 14395, 15745, 8645]

for name, A, m, p in [("P_7", A_7, 3, 7), ("P_11", A_11, 5, 11)]:
    ratios = [A[d+1]/A[d] for d in range(len(A)-1)]
    print(f"  {name}: A_{d+1}/A_d ratios: {[f'{r:.4f}' for r in ratios]}")

# CRITICAL: Check if Omega = convolution of two sequences
print("\nConvolution analysis:")
for name, omega, m, p in [("P_7", omega_7, 3, 7), ("P_11", omega_11, 5, 11)]:
    f = [comb(m, k) for k in range(m+1)]
    g_conv = np.convolve(f, f)
    print(f"  {name}: C(m,k)*C(m,k) convolution = {list(g_conv.astype(int))}")
    print(f"  {name}: Omega = {omega}")
    if len(g_conv) == len(omega):
        ratios = [omega[d] / g_conv[d] if g_conv[d] > 0 else 0 for d in range(len(omega))]
        print(f"    ratio: {[f'{r:.4f}' for r in ratios]}")

# Check Omega vs C(2m, d)
print("\nOmega vs C(2m, d):")
for name, omega, m in [("P_7", omega_7, 3), ("P_11", omega_11, 5)]:
    binom = [comb(2*m, d) for d in range(2*m+1)]
    print(f"  {name}: C(2m,d) = {binom}")
    ratios = [omega[d] / binom[d] if binom[d] > 0 else 0 for d in range(len(omega))]
    print(f"    Omega/C(2m,d) = {[f'{r:.4f}' for r in ratios]}")

# Check: Omega_d = product_{j=0}^{d-1} (m - f(j)) for some f
print("\nProduct form: Omega_d = Omega_{d-1} * r_d")
for name, omega, m, p in [("P_7", omega_7, 3, 7), ("P_11", omega_11, 5, 11)]:
    ratios = [omega[d]/omega[d-1] if omega[d-1] > 0 else 0 for d in range(1, len(omega))]
    print(f"  {name}: r_d = {[f'{r:.4f}' for r in ratios]}")
    # Try r_d = (a*m - b*d + c) / d
    for d, r in enumerate(ratios, 1):
        # r_d = (stuff) / d
        print(f"    d={d}: r*d = {r*d:.4f}, r*(d+1) = {r*(d+1):.4f}", end="")
        # Check: is r_d = (2m - d + 1) * something?
        check1 = (2*m - d + 1)
        check2 = (p - d)
        print(f"  (2m-d+1)={check1}, (p-d)={check2}, r/(2m-d+1)={r/check1:.4f} r/((p-d)/d)={r*d/(p-d):.4f}" if check2 > 0 else "")
