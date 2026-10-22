# Loop expansion beyond regular tournaments: for almost-regular even-n tournaments (scores n/2 and n/2 - 1) the tree
# level (path trees) X_(-1)(z) = -z^2 delta^T (I + z^2 Shat^2)^-1 delta, delta = Shat 1, is O(1/n), the same order as two loops.
# Test: residual y = n(1 - kappa/rho) against y_pred (two-loop formula, as if regular) and t = -n X_(-1)(1).
# Usage: almost_regular_tree.py n HK.  mac-mini-2026-10-22-sierp.
import subprocess, sys, math
import numpy as np
sys.path.insert(0, '.')
from two_loop_regular import two_loop
n = int(sys.argv[1]); HK = sys.argv[2]; h = n // 2
gen = subprocess.run(['gentourng', '-q', f'-d{h - 1}', f'-D{h}', str(n)], capture_output=True, text=True).stdout
codes = [c for c in gen.split() if c]
Hs = [int(x) for x in subprocess.run([HK], input=gen, capture_output=True, text=True).stdout.split()]
iu = np.triu_indices(n, 1)
Y, P, T = [], [], []
for c, H in zip(codes, Hs):
    A = np.zeros((n, n)); b = np.array([1 if ch == '1' else 0 for ch in c]); A[iu[0], iu[1]] = b; A[iu[1], iu[0]] = 1 - b
    S = A - A.T; Sh = S / n
    rho, yp, *_ = two_loop(S)
    delta = Sh @ np.ones(n)
    t = n * delta @ np.linalg.solve(np.eye(n) + Sh @ Sh, delta)
    kap = H * 2 ** (n - 1) / math.factorial(n)
    Y.append(n * (1 - kap / rho)); P.append(yp); T.append(t)
Y, P, T = map(np.array, (Y, P, T))
print(f"n={n}: {len(Y)} almost-regular tournaments; path-tree term t = n delta^T (I + Shat^2)^-1 delta in [{T.min():.4f}, {T.max():.4f}]")
print(f"  corr(y, y_pred) = {np.corrcoef(Y, P)[0, 1]:+.4f};  corr(y, y_pred + t) = {np.corrcoef(Y, P + T)[0, 1]:+.4f};  "
      f"mean(y - y_pred) = {np.mean(Y - P):+.4f}, mean(y - y_pred - t) = {np.mean(Y - P - T):+.4f}; sd {np.std(Y - P):.4f} -> {np.std(Y - P - T):.4f}")
