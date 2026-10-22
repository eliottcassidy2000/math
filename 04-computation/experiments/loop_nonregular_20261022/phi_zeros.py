"""Zeros of the exact forest polynomial and the reach of the cluster expansion.  Session opus-2026-10-22-S54.

Phi(z) = sum_L N_L (z/n)^L is a polynomial (L <= n - 1, odd L vanish), with N_L = c_L/(n - L)! from the exact forest
numbers (ar_tool forest).  log Phi(z) = X0(z) + (all cluster classes); X0 = log rho(z) is analytic for
|z| < 1/max|x_j| (>= sqrt(2n/(n - 1)) > 1.41), so the radius of convergence of the z-series of the cluster sum
log Phi - X0 is min(|zero of Phi|, 1/max|x_j|).  Prints, per file, the distribution of the smallest zero modulus and
the exact cluster part A = log Phi(1) - log rho.
Usage: python phi_zeros.py n FILE [MAX]"""
import sys, os, math, subprocess
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from predict_nonregular import load

n = int(sys.argv[1]); path = sys.argv[2]; MAX = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
codes = [line.split()[0] for i, line in enumerate(open(path)) if line.strip() and i < MAX]
tool = os.path.join(HERE, 'ar_tool2.exe')
cl = subprocess.run([tool, 'forest', str(n)], input='\n'.join(codes) + '\n', capture_output=True, text=True).stdout.split('\n')
rows = []
for (code, H, S), line in zip(load(n, path), cl[:len(codes)]):
    c = [int(t) for t in line.split()]
    N = [c[L] / math.factorial(n - L) for L in range(n)]
    assert all(abs(N[L]) < 1e-9 for L in range(1, n, 2))
    # polynomial in w = z^2: P(w) = sum_k N_2k (w/n^2)^k
    coef = [N[2 * k] / n ** (2 * k) for k in range((n - 1) // 2 + 1)]
    roots = np.roots(coef[::-1]) if len(coef) > 1 else np.array([])
    zmin = float(np.min(np.sqrt(np.abs(roots)))) if len(roots) else np.inf
    w0 = roots[np.argmin(np.abs(roots))] if len(roots) else 0
    argw = float(np.angle(w0)) if len(roots) else 0.0
    x = np.linalg.eigvals(S / n).imag
    logrho = 0.5 * np.sum(np.log((1 + x ** 2) / (1 - x ** 2)))
    Phi1 = sum(coef)
    rows.append(dict(zmin=zmin, xmax=float(np.max(np.abs(x))), A=math.log(Phi1) - logrho if Phi1 > 0 else np.nan,
                     neg=int(np.sum((roots.real > 0) & (abs(roots.imag) < 1e-12))), argw=argw))
Z = np.array([r['zmin'] for r in rows]); A = np.array([r['A'] for r in rows])
print(f'n={n}: {len(rows)} tournaments ({os.path.basename(path)})')
print(f'  smallest zero modulus of Phi: min {Z.min():.4f}, 5% {np.quantile(Z, .05):.4f}, median {np.median(Z):.4f}, max {Z.max():.4f}; '
      f'fraction with a zero in |z| <= 1: {np.mean(Z <= 1):.4f}; in |z| <= 1.2: {np.mean(Z <= 1.2):.4f}')
print(f'  spectral singularity 1/max|x|: min {1 / max(r["xmax"] for r in rows):.4f}; exact cluster part A: mean {np.nanmean(A):+.5f}, sd {np.nanstd(A):.5f}')
print(f'  zeros on the positive w = z^2 axis (real z): {sum(r["neg"] for r in rows)} in total')
AW = np.array([r['argw'] for r in rows])
print(f'  argument of the smallest zero in the w = z^2 plane: w > 0 (real z) {np.mean(np.abs(AW) < 1e-9):.3f}, w < 0 (imaginary z) '
      f'{np.mean(np.abs(np.abs(AW) - np.pi) < 1e-9):.3f}, complex {np.mean((np.abs(AW) > 1e-9) & (np.abs(np.abs(AW) - np.pi) > 1e-9)):.3f}; '
      f'median |arg w| {np.median(np.abs(AW)):.4f}')
