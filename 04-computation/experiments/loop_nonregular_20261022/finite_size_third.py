"""Third-order finite-size term under the one-loop forest measure.  Session opus-2026-10-22-S54.

log g_n(L) = sum_(i < L) -log(1 - i/n) = a/n + b/n^2 + c/n^3 + ..., a = s_1, b = s_2/2, c = s_3/3, s_j = sum_(i < L) i^j.
log E[g] = E a/n + (Var a/2 + E b)/n^2 + (k3(a)/6 + Cov(a, b) + E c)/n^3 + O(n^-4)  under any measure on L.
With the one-loop measure w_L = [z^L] rho(z)/rho(1) (Taylor coefficients of exp X0), the n^-3 coefficient B3_0 is the
leading part of the finite-size remainder B - B_f after second order (the first-order cluster shift of the n^-2 terms
is of the same order and is reported separately as the derivative term).  Compared here with the exact n^3 (B - B_f)
from budget files.  Usage: python finite_size_third.py n BUDGETFILE"""
import sys, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))

n = int(sys.argv[1]); path = sys.argv[2]
iu = np.triu_indices(n, 1)
K = 128
zs = np.exp(2j * np.pi * np.arange(K) / K)
Ls = np.arange(K)
s1 = Ls * (Ls - 1) / 2
s2 = np.array([sum(i * i for i in range(L)) for L in Ls], dtype=float)
s3 = s1 ** 2
a, b, c = s1, s2 / 2, s3 / 3
B3, exact = [], []
rows = [l.split() for l in open(path) if l.strip()]
head = rows[0]
for r in rows[1:]:
    d = dict(zip(head, r))
    code = d['code']
    A = np.zeros((n, n)); bb = np.array([1 if ch == '1' else 0 for ch in code])
    A[iu[0], iu[1]] = bb; A[iu[1], iu[0]] = 1 - bb
    x2 = np.linalg.eigvalsh(1j * (A - A.T) / n) ** 2
    X0 = np.array([0.5 * np.sum(np.log1p(z * z * x2) - np.log1p(-z * z * x2)) for z in zs])
    w = (np.fft.fft(np.exp(X0)) / K).real
    w = w / w.sum()
    Ea = w @ a; Va = w @ (a - Ea) ** 2; k3 = w @ (a - Ea) ** 3
    Eb = w @ b; Cab = w @ ((a - Ea) * (b - Eb)); Ec = w @ c
    B3.append(k3 / 6 + Cab + Ec)
    exact.append(n ** 3 * (float(d['B']) - float(d['Bf'])))
B3 = np.array(B3); exact = np.array(exact)
print(f'n={n}: {len(B3)} tournaments: one-loop third-order finite-size coefficient B3_0 mean {B3.mean():+.2f} (sd {B3.std():.2f}); '
      f'exact n^3 (B - B_f) mean {exact.mean():+.2f} (sd {exact.std():.2f}); corr {np.corrcoef(B3, exact)[0, 1]:+.3f}')
