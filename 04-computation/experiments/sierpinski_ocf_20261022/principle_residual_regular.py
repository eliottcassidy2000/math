# Second-order law on ALL regular tournaments at order n: y = n(1 - kappa/rho) against M4 = sum_k sigma_k^4/n^4.
# usage: principle_residual_regular.py n HK_BINARY.  mac-mini-2026-10-22-sierp.
import subprocess, sys, math
import numpy as np
n = int(sys.argv[1]); HK = sys.argv[2]; m = (n - 1) // 2
gen = subprocess.run(['gentourng', f'-d{m}', f'-D{m}', str(n)], capture_output=True, text=True).stdout
codes = [c for c in gen.split() if not c.startswith('>')]
Hs = np.array([int(x) for x in subprocess.run([HK], input=gen, capture_output=True, text=True).stdout.split()], dtype=float)
iu = np.triu_indices(n, 1)
B = np.zeros((len(codes), n, n)); bits = np.array([[1 if ch == '1' else 0 for ch in c] for c in codes])
B[:, iu[0], iu[1]] = bits; B[:, iu[1], iu[0]] = 1 - bits
s2 = np.sort(np.abs(np.linalg.eigvals(B - np.transpose(B, (0, 2, 1))).imag) ** 2, axis=1)[:, 1::2]
t = s2 / n ** 2
rho = np.prod((1 + t) / (1 - t), axis=1); M4 = np.sum(t ** 2, axis=1); M6 = np.sum(t ** 3, axis=1)
kap = Hs * 2 ** (n - 1) / math.factorial(n)
y = n * (1 - kap / rho)
b, a = np.polyfit(M4, y, 1)
k2 = rho * (1 - (a + b * M4) / n)
rk = lambda v: np.argsort(np.argsort(v))
print(f"n={n} ({len(codes)} regular classes): corr(y, M4) = {np.corrcoef(y, M4)[0,1]:+.3f}; fit y = {a:.4f} + {b:.4f} M4; "
      f"spearman(kappa, rho) = {np.corrcoef(rk(kap), rk(rho))[0,1]:+.3f} -> spearman(kappa, kappa2) = {np.corrcoef(rk(kap), rk(k2))[0,1]:+.3f}; "
      f"max rel err {np.max(np.abs(k2/kap - 1)):.2e}; argmax agree: {int(np.argmax(kap)) == int(np.argmax(k2))}")
