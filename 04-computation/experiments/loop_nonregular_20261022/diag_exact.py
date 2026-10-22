"""Exact error budget of the loop expansion on non-regular tournaments.  Session opus-2026-10-22-S54.

log(kappa/rho) = [log Phi(1) - log rho] + log(kappa/Phi(1))   (cluster part A, finite-size part B), both exact from
the forest numbers c_L = (n - L)! N_L (ar_tool forest): Phi(1) = sum c_L/((n - L)! n^L), kappa = sum c_L/n!.
  A = sum of all non-cycle classes X_(e,k)(1);  truncations A1 = X(-1,2) + X(1,0); A2 = A1 + X(0,1) + X(0,2) + X(2,0);
      A3 = A2 + X(-1,3) + X(-1,4) + X(0,3) + X(1,1).
  B = log E[g_n(L)] under the exact measure phi_L/Phi(1); expansions B1 = E0[a]/n, B2 = (E0[a] + dE[a])/n + (Var0(a)/2 +
      E0[b])/n^2 (one-loop measure E0 plus the first-order cumulant shift), and Bx = log E[g_n(L)] with the exact E.
Prints mean and sd of each error over the file, and the exact-vs-predicted y in log form,
y2L = n (1 - exp(A2 + B2)).  Usage: python diag_exact.py n FILE"""
import sys, os, math, subprocess
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leaf_engine as LE
from predict_nonregular import cumulants, X1_closed, tree_closed, load

n = int(sys.argv[1]); path = sys.argv[2]
codes = [line.split()[0] for line in open(path) if line.strip()]
tool = os.path.join(HERE, 'ar_tool2.exe')
cl = subprocess.run([tool, 'forest', str(n)], input='\n'.join(codes) + '\n', capture_output=True, text=True).stdout.split('\n')
rows = []
for (code, H, S), line in zip(load(n, path), cl):
    c = [int(t) for t in line.split()]
    assert sum(c) == H * 2 ** (n - 1), (code, sum(c), H)
    Sh = S / n
    x = np.linalg.eigvals(Sh).imag
    logrho = 0.5 * np.sum(np.log((1 + x ** 2) / (1 - x ** 2)))
    phi = np.array([c[L] / (math.factorial(n - L) * n ** L) for L in range(n)])
    Phi1 = phi.sum(); kap = sum(c) / math.factorial(n)
    A = math.log(Phi1) - logrho; B = math.log(kap / Phi1)
    Ls = np.arange(n); w = phi / Phi1
    g = np.array([n ** L / math.prod(range(n - L + 1, n + 1)) if L else 1.0 for L in range(n)])
    aL = Ls * (Ls - 1) / 2
    Ea_exact = float(w @ aL)
    k = cumulants(lambda z: 0.5 * np.sum(np.log((1 + (z * x) ** 2) / (1 - (z * x) ** 2))))
    k1, k2, k3, k4 = k[1], k[2], k[3], k[4]
    EL = k1; EL2 = k2 + k1 ** 2; EL3 = k3 + 3 * k2 * k1 + k1 ** 3
    EL4 = k4 + 4 * k3 * k1 + 3 * k2 ** 2 + 6 * k2 * k1 ** 2 + k1 ** 4
    Ea = (EL2 - EL) / 2; Ea2 = (EL4 - 2 * EL3 + EL2) / 4; Eb = (2 * EL3 - 3 * EL2 + EL) / 12
    P = LE.LEF.stack(Sh.astype(complex), 1.0)
    X = {cc: LE.X(Sh, cc[0], cc[1], 1.0, P).real for cc in [(-1, 2), (1, 0), (0, 1), (0, 2), (2, 0), (-1, 3), (-1, 4), (0, 3), (1, 1)]}
    Dl = cumulants(lambda z: tree_closed(Sh.astype(complex), z) + X1_closed(Sh.astype(complex), z), kmax=2)
    dEa = (Dl[2] + (2 * k1 - 1) * Dl[1]) / 2
    A1 = X[(-1, 2)] + X[(1, 0)]
    A2 = A1 + X[(0, 1)] + X[(0, 2)] + X[(2, 0)]
    A3 = A2 + X[(-1, 3)] + X[(-1, 4)] + X[(0, 3)] + X[(1, 1)]
    B1 = Ea / n
    B2 = (Ea + dEa) / n + ((Ea2 - Ea ** 2) / 2 + Eb) / n ** 2
    y = n * (1 - kap / math.exp(logrho))
    rows.append(dict(A=A, A1=A1, A2=A2, A3=A3, B=B, B1=B1, B2=B2, Ea_exact=Ea_exact, Ea0=Ea, Ea1=Ea + dEa, y=y,
                     y1L=n * (1 - math.exp(A1 + B1)), y2L=n * (1 - math.exp(A2 + B2)), y3L=n * (1 - math.exp(A3 + B2))))
R = {k: np.array([r[k] for r in rows]) for k in rows[0]}
print(f'n={n}: {len(rows)} tournaments ({os.path.basename(path)})')
print(f'  cluster part A = log Phi(1) - log rho: mean {R["A"].mean():+.5f}; truncation errors (n x): '
      f'A - A1 {n * (R["A"] - R["A1"]).mean():+.4f} (sd {n * (R["A"] - R["A1"]).std():.4f}), '
      f'A - A2 {n * (R["A"] - R["A2"]).mean():+.4f} (sd {n * (R["A"] - R["A2"]).std():.4f}), '
      f'A - A3 {n * (R["A"] - R["A3"]).mean():+.4f} (sd {n * (R["A"] - R["A3"]).std():.4f})')
print(f'  finite-size part B = log(kappa/Phi(1)): mean {R["B"].mean():+.5f}; errors (n x): B - B1 {n * (R["B"] - R["B1"]).mean():+.4f} '
      f'(sd {n * (R["B"] - R["B1"]).std():.4f}), B - B2 {n * (R["B"] - R["B2"]).mean():+.4f} (sd {n * (R["B"] - R["B2"]).std():.4f})')
print(f'  E[a]: exact {R["Ea_exact"].mean():.4f}, one-loop E0[a] {R["Ea0"].mean():.4f}, with first-order shift {R["Ea1"].mean():.4f}')
for key in ('y1L', 'y2L', 'y3L'):
    res = R['y'] - R[key]
    print(f'  {key}: corr(y, pred) {np.corrcoef(R["y"], R[key])[0, 1]:+.4f}, mean(y - pred) {res.mean():+.4f}, sd {res.std():.4f}')
