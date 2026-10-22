"""Are the zeros of the forest polynomial Szego zeros of the truncated one-loop series?  Session opus-2026-10-22-S54.

Compare the smallest zero of the exact Phi(z) (exact forest numbers) with the smallest zero of the Taylor truncation of
rho(z) = exp(X0(z)) at the same degree n - 1, and of the truncated Gaussian exp((n - 1) z^2 / n) (the flat limit of rho).
All three are polynomials in w = z^2 of degree floor((n - 1)/2).  Usage: python phi_zeros_model.py n FILE [MAX]"""
import sys, os, math, subprocess
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from predict_nonregular import load


def exp_series(a, m):
    """coefficients of exp(sum_k a_k w^k) to degree m (a_0 = 0)."""
    b = [1.0] + [0.0] * m
    for k in range(1, m + 1):
        b[k] = sum(j * a[j] * b[k - j] for j in range(1, k + 1)) / k
    return b


def smallest_root(coef):
    r = np.roots(coef[::-1])
    i = np.argmin(np.abs(r))
    return r[i]


n = int(sys.argv[1]); path = sys.argv[2]; MAX = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
codes = [line.split()[0] for i, line in enumerate(open(path)) if line.strip() and i < MAX]
cl = subprocess.run([os.path.join(HERE, 'ar_tool2.exe'), 'forest', str(n)], input='\n'.join(codes) + '\n',
                    capture_output=True, text=True).stdout.split('\n')
m = (n - 1) // 2
out = []
for (code, H, S), line in zip(load(n, path), cl[:len(codes)]):
    c = [int(t) for t in line.split()]
    exact = [c[2 * k] / math.factorial(n - 2 * k) / n ** (2 * k) for k in range(m + 1)]
    x2 = np.linalg.eigvals(S / n).imag ** 2
    a = [0.0] + [(np.sum(x2 ** k) / k if k % 2 else 0.0) for k in range(1, m + 1)]   # log rho = sum over all eigenvalues of atanh(w x_j^2)
    trunc = exp_series(a, m)
    gauss = [((n - 1) / n) ** k / math.factorial(k) for k in range(m + 1)]
    out.append((smallest_root(exact), smallest_root(trunc), smallest_root(gauss)))
E = np.array([o[0] for o in out]); T = np.array([o[1] for o in out]); G = np.array([o[2] for o in out])
f = lambda v: f'|z| median {np.median(np.sqrt(np.abs(v))):.4f} (range {np.min(np.sqrt(np.abs(v))):.4f}..{np.max(np.sqrt(np.abs(v))):.4f}), median arg(w) {np.median(np.abs(np.angle(v))):.3f}'
print(f'n={n} ({len(out)} tournaments, {os.path.basename(path)}), degree {m} in w = z^2:')
print(f'  exact Phi:             {f(E)}')
print(f'  truncated rho(z):      {f(T)};  median |w_exact/w_trunc| {np.median(np.abs(E / T)):.4f}')
print(f'  truncated exp(c z^2):  {f(G)}')
