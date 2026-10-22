"""Coefficient-wise loop prediction: exact finite-size factors on the truncated cluster series.  Session
opus-2026-10-22-S54.

kappa = sum_(L <= n-1) phi_L g_n(L) exactly, g_n(L) = n^L/(n)_L, Phi = sum phi_L z^L = exp(sum of all classes).
Replace log Phi by a truncation X_trunc(z) = X_(0,0) + chosen leaf classes, take psi_L = [z^L] exp(X_trunc) (Cauchy sums
on |z| = r) and predict kappa_cw = sum_(L <= n-1) psi_L g_n(L): no 1/n expansion of g_n and the support L <= n - 1 of
the exact forest measure is respected. Truncations:
  T0 = cycles only;  T1 = T0 + X(-1,2) + X(1,0);  T2 = T1 + X(0,1) + X(0,2) + X(2,0);  T3 = T2 + X(-1,3) + X(-1,4) + X(0,3) + X(1,1).
Usage: python predict_coeffwise.py n FILE [MAX] [K]"""
import sys, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leaf_engine as LE
from predict_nonregular import X1_closed, tree_closed, load

TR = {'T0': [], 'T1': [(-1, 2), (1, 0)], 'T2': [(-1, 2), (1, 0), (0, 1), (0, 2), (2, 0)],
      'T3': [(-1, 2), (1, 0), (0, 1), (0, 2), (2, 0), (-1, 3), (-1, 4), (0, 3), (1, 1)]}


def class_at(Sh, c, z, P):
    if c == (-1, 2):
        return tree_closed(Sh.astype(complex), z)
    if c == (1, 0):
        return X1_closed(Sh.astype(complex), z)
    return LE.X(Sh, c[0], c[1], z, P)


def predict_cw(S, K=64, r=1.0):
    n = S.shape[0]; Sh = S / n
    x = np.linalg.eigvals(Sh).imag
    rho = float(np.sqrt(np.prod((1 + x ** 2) / (1 - x ** 2))))
    zs = r * np.exp(2j * np.pi * np.arange(K) / K)
    vals = {}
    need = TR['T3']
    for kk, z in enumerate(zs):
        P = LE.LEF.stack(Sh.astype(complex), z) if any(c not in ((-1, 2), (1, 0)) for c in need) else None
        X0 = 0.5 * np.sum(np.log((1 + (z * x) ** 2) / (1 - (z * x) ** 2)))
        vals.setdefault('X0', []).append(X0)
        for c in need:
            vals.setdefault(c, []).append(class_at(Sh, c, z, P))
    g = np.array([n ** L / math.prod(range(n - L + 1, n + 1)) if L else 1.0 for L in range(n)])
    out = {}
    for name, cls in TR.items():
        Xz = np.array(vals['X0']) + sum(np.array(vals[c]) for c in cls)
        coef = np.fft.fft(np.exp(Xz)) / K / r ** np.arange(K)
        psi = coef.real[:n]
        out[name] = float(psi @ g)
    return rho, out


if __name__ == '__main__':
    n = int(sys.argv[1]); path = sys.argv[2]
    MAX = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
    K = int(sys.argv[4]) if len(sys.argv) > 4 else 64
    Y = []; PR = {k: [] for k in TR}
    for i, (code, H, S) in enumerate(load(n, path)):
        if i >= MAX:
            break
        rho, out = predict_cw(S, K)
        kap = H * 2 ** (n - 1) / math.factorial(n)
        Y.append(n * (1 - kap / rho))
        for k, v in out.items():
            PR[k].append(n * (1 - v / rho))
    Y = np.array(Y)
    print(f'n={n}: {len(Y)} tournaments ({os.path.basename(path)}), coefficient-wise prediction with exact g_n(L), L <= n-1')
    for k, v in PR.items():
        v = np.array(v); R = Y - v
        print(f'  {k}: corr {np.corrcoef(Y, v)[0, 1]:+.4f}, mean(y - pred) {R.mean():+.4f}, sd {R.std():.4f}, max {np.abs(R).max():.4f}', flush=True)
