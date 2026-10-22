# Audit of THM-4755 (4) / HYP-9364 (a): tower normalised counts, Richardson in 1/n, deviation from rho(S_k).
# kappa values for k <= 5 recomputed by the auditor (b_tower_H.out, b_tower_k5_modp.out); k = 6 is the session's value.
import re
from fractions import Fraction
from math import factorial
from mpmath import mp, mpf, matrix, lu_solve, exp, log, atanh
mp.dps = 40
s = open('/Users/e/Documents/GitHub/math-wt-fanopsp/04-computation/experiments/sierpinski_ocf_20261022/ternary_fast_k6.out').read()
H = {int(m.group(1)): int(m.group(3)) for m in re.finditer(r'k=(\d) n=(\d+): H\(S_k\) = (\d+)', s)}
K = {}
for k, h in H.items():
    n = 3 ** k
    f = Fraction(h * 2 ** (n - 1), factorial(n))
    K[k] = mpf(f.numerator) / f.denominator
rho = lambda k: exp(sum(2 * mpf(3) ** i * atanh(mpf(3) ** (-(2 * i + 1))) for i in range(k)))
kS = exp(sum(2 * mpf(3) ** i * atanh(mpf(3) ** (-(2 * i + 1))) for i in range(80)))
def fit0(ks):
    A = matrix(len(ks), len(ks)); b = matrix(len(ks), 1)
    for i, k in enumerate(ks):
        for j in range(len(ks)): A[i, j] = mpf(3) ** (-k * j)
        b[i] = K[k]
    return lu_solve(A, b)[0]
print("kappa(S_k):", {k: mp.nstr(K[k], 8) for k in sorted(K)})
print("kappa_S =", mp.nstr(kS, 22))
for ks in [(3, 4, 5, 6), (4, 5, 6), (2, 3, 4, 5, 6), (2, 3, 4, 5), (3, 4, 5), (5, 6)]:
    v = fit0(ks)
    print(f"Richardson in 1/n through k={ks}: {mp.nstr(v, 13)}   minus kappa_S = {mp.nstr(v - kS, 3)}")
print("differences ratios:", [mp.nstr((K[k + 1] - K[k]) / (K[k] - K[k - 1]), 4) for k in range(2, 6)])
for k in range(1, 7):
    d = 1 - K[k] / rho(k)
    print(f"k={k} n={3**k}: rho(S_k)={mp.nstr(rho(k), 10)}  1-kappa/rho = {mp.nstr(d, 4)}  n*(1-kappa/rho) = {mp.nstr(d * 3 ** k, 4)}")
