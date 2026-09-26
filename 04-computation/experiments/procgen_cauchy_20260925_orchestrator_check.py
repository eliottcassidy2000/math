#!/usr/bin/env python3
"""Orchestrator's independent audit of THM-4476 (Cauchy-Schwarz price bound), written without the lane's code.
  C1  exact M2(8) = E_Haar[W_8^2] = 501975/4096 from the backward tree over v mod 3^7 (exact rationals).
  C2  Perron roots of the class-mass majorant Phi at levels r = 1..10 by converged power iteration (max = min ratio
      on the non-zero classes): theta_1 = 1.34307, theta_10 = 1.04911.
  C3  rigorous exact-rational certificate at r = 10 on the non-zero classes (upper-rounded square roots; zero classes
      evolve exactly, ratio 1/4): theta_10 <= 1.0491582, hence the exponent 2(1-h) + log2(theta_10) <= 0.1694.
"""
import math
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
@lru_cache(maxsize=None)
def Z(k, v):
    if k == 0: return F(1)
    s = Z(k-1, 2*v) / 2
    if v % 3 == 2: s += F(3, 2) * Z(k-1, (2*v - 1)//3)
    return s
L = 8; M = 3**(L-1)
M2 = sum((sum(Z(k, v) for k in range(L)))**2 for v in range(M, 2*M)) / M
assert M2 == F(501975, 4096); print("C1  M2(8) = 501975/4096 exactly: ok")
Z.cache_clear()
SQ = math.sqrt(3)/2
def phi(e, r):
    N = 3**r; Nb = 3**(r-1); c = np.arange(N); twoc = (2*c) % N
    t = np.arange(Nb); Ec = (2*t + 1) % Nb
    U = e.reshape(3, Nb).sum(axis=0); new = e[twoc]*0.25
    b = U[Ec]; new[2::3] += np.sqrt(e[twoc[2::3]]*b)*SQ + 0.75*b
    return new
e = None; th = {}
for r in range(1, 11):
    N = 3**r; e = np.ones(N) if e is None else np.tile(e, 3)/3.0
    for _ in range(3000):
        f = phi(e, r); e = f/f.sum()
    nz = (np.arange(N) % 3 != 0); ratio = phi(e, r)[nz]/e[nz]
    assert ratio.max() - ratio.min() < 1e-9
    th[r] = ratio.max()
assert abs(th[1] - 1.3430703) < 1e-6 and abs(th[10] - 1.0491107) < 1e-6
print(f"C2  Perron roots: theta_1 = {th[1]:.7f}, theta_5 = {th[5]:.7f}, theta_10 = {th[10]:.7f}: ok")
R = 10; N = 3**R; Nb = 3**(R-1)
nz = (np.arange(N) % 3 != 0); eps = 1e-9*e[nz].min(); e = np.where(nz, e, eps)
S = 10**15
def sqrt_up(x):
    n = x.numerator*S*S//x.denominator + 1
    return F(math.isqrt(n) + 1, S)
s3 = sqrt_up(F(3))/2; eF = [F(float(x)) for x in e]; best = F(0)
for c in range(N):
    if c % 3 == 0: continue
    tc = (2*c) % N; val = eF[tc]/4
    if c % 3 == 2:
        d = (2*((c-2)//3) + 1) % Nb; UU = eF[d] + eF[d+Nb] + eF[d+2*Nb]
        val += F(3, 4)*UU + s3*sqrt_up(eF[tc]*UU)
    best = max(best, val/eF[c])
eta = 1 - (-(math.log(2)/math.log(3))*math.log2(math.log(2)/math.log(3)) - (1-math.log(2)/math.log(3))*math.log2(1-math.log(2)/math.log(3)))
assert float(best) < 1.04916
print(f"C3  rigorous theta_10 <= {float(best):.7f}; exponent 2*eta + log2(theta_10) = {2*eta + math.log2(float(best)):.4f}: ok")
print("ALL CHECKS PASSED")
