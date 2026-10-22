# Audit of note section 10 (uncommitted): kappa/rho over ALL regular n-tournaments, n = 9, 11, and Spearman(kappa, rho).
import sys
import numpy as np
from math import factorial
sys.path.insert(0, '/Users/e/Documents/GitHub/math-wt-fanopsp/04-computation/experiments/sierpinski_ocf_20261022/audit')
from d_circ_band import spearman
def run(n, fname):
    rows = [l.split() for l in open(fname) if l.strip()]
    kap, rho = [], []
    for code, h in rows:
        A = np.zeros((n, n)); b = 0
        for i in range(n):
            for j in range(i + 1, n):
                if code[b] == '1': A[i, j] = 1
                else: A[j, i] = 1
                b += 1
        ev = np.linalg.eigvals(A - A.T)
        r = np.prod([(n * n + s * s) / (n * n - s * s) for s in np.abs(ev.imag) if s > 1e-9]) ** 0.5
        kap.append(int(h) * 2 ** (n - 1) / factorial(n)); rho.append(r)
    q = [k / r for k, r in zip(kap, rho)]
    print(f"n={n}: regular classes={len(rows)}  kappa/rho in [{min(q):.4f}, {max(q):.4f}]  spearman(kappa, rho) = {spearman(kap, rho):+.3f}")
run(9, sys.argv[1]); run(11, sys.argv[2])
