# Sato-Tate moments of the four new seven-vertex motives m1 = N02, m2 = N09, m3 = N18, m4 = N05 + N18 (normalised new
# parts: residuals of the q^-5-weighted fit on the base, divided by q^(5/2)).  References for E[x^2], E[x^4], E[x^6]:
# SU(2) 1, 2, 5; USp(4) 1, 3, 14; SU(2) x SU(2) 1, 4, 25.  mac-mini-2026-10-25-pyth.
import sys
import numpy as np
exec(open('analyze_new.py').read().split('R = (Yw - Bw @ coef)')[0].replace('print(f"  {nm}:', 'pass  # print(f"  {nm}:').replace('print(f"{len(qs)} primes', 'pass  # print('))
R = (Yw - Bw @ coef)
m = np.stack([R[:, 0], R[:, 2], R[:, 3], R[:, 1] + R[:, 3]], axis=1)
print(f"{len(qs)} primes in [{int(qs.min())}, {int(qs.max())}]")
for j, nm in enumerate(['m1 (Y02)', 'm2 (Y09)', 'm3 (Y18)', 'm4 (Y05+Y18)']):
    x = m[:, j]
    print(f"  {nm}: E[x^2] = {np.mean(x**2):.3f}, E[x^4] = {np.mean(x**4):.3f}, E[x^6] = {np.mean(x**6):.3f}, max|x| = {np.max(np.abs(x)):.3f}, mean = {np.mean(x):+.3f}")
G = m.T @ m / len(qs)
np.set_printoptions(precision=3, suppress=True)
print("Gram of (m1..m4):\n", G)
