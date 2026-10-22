# Residual of the self-blow-up principle on circulant classes: y = n(1 - kappa/rho).  Correlate y with candidate
# second-order statistics: log rho, M4 = sum_j sigma_j^4 / n^4, M6 = sum_j sigma_j^6/n^6, and the largest sigma^2/n^2.
# mac-mini-2026-10-22-sierp.
import re, glob, math
import numpy as np
from math import factorial
from collections import defaultdict
rows = defaultdict(dict)
for f in sorted(glob.glob('../tournament_pseudoprimes_20261021/circ_*.out')):
    for line in open(f):
        m = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
        if m: rows[int(m.group(1))][tuple(sorted(int(x) for x in m.group(2).split(',')))] = int(m.group(4))
for line in open('../lattice_bridges_20261021/circ25_H.txt'):
    m = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
    if m: rows[25][tuple(int(x) for x in m.group(1).split(','))] = int(m.group(2))
for n in (17, 19, 21, 23, 25):
    Y, LR, M4, M6, TOP = [], [], [], [], []
    for S, H in rows[n].items():
        sig = np.array([2 * sum(math.sin(2 * math.pi * j * x / n) for x in S) for j in range(1, (n - 1) // 2 + 1)])
        t = sig ** 2 / n ** 2
        rho = float(np.prod((1 + t) / (1 - t)))
        kap = H * 2 ** (n - 1) / factorial(n)
        Y.append(n * (1 - kap / rho)); LR.append(math.log(rho)); M4.append(float(np.sum(t ** 2))); M6.append(float(np.sum(t ** 3))); TOP.append(float(t.max()))
    Y = np.array(Y)
    def r(x): x = np.array(x); return np.corrcoef(Y, x)[0, 1]
    # best linear fit on M4 and M6
    X = np.column_stack([np.ones_like(Y), M4, M6])
    coef, res, *_ = np.linalg.lstsq(X, Y, rcond=None)
    pred = X @ coef
    print(f"n={n} ({len(Y)} classes): y in [{Y.min():.3f}, {Y.max():.3f}]; corr(y, log rho) = {r(LR):+.3f}, corr(y, M4) = {r(M4):+.3f}, corr(y, M6) = {r(M6):+.3f}, corr(y, max t) = {r(TOP):+.3f}; "
          f"fit y ~ a + b M4 + c M6: R^2 = {1 - np.sum((Y - pred) ** 2) / np.sum((Y - Y.mean()) ** 2):.3f}", flush=True)

# Second stage: fit y = a + b*M4 per n; then test whether the corrected predictor kappa2 = rho*(1 - (a+b*M4)/n)
# ranks the classes (Spearman with kappa) and finds the argmax.
print()
for n in (13, 15, 17, 19, 21, 23, 25):
    data = []
    for S, H in rows[n].items():
        sig = np.array([2 * sum(math.sin(2 * math.pi * j * x / n) for x in S) for j in range(1, (n - 1) // 2 + 1)])
        t = sig ** 2 / n ** 2
        rho = float(np.prod((1 + t) / (1 - t)))
        kap = H * 2 ** (n - 1) / factorial(n)
        data.append((kap, rho, float(np.sum(t ** 2))))
    kap = np.array([d[0] for d in data]); rho = np.array([d[1] for d in data]); M4 = np.array([d[2] for d in data])
    Y = n * (1 - kap / rho)
    b, a = np.polyfit(M4, Y, 1)
    k2 = rho * (1 - (a + b * M4) / n)
    rk = lambda v: np.argsort(np.argsort(v))
    sp1 = np.corrcoef(rk(kap), rk(rho))[0, 1]; sp2 = np.corrcoef(rk(kap), rk(k2))[0, 1]
    relerr = np.max(np.abs(k2 / kap - 1))
    print(f"n={n}: y = {a:.4f} + {b:.4f} M4; spearman(kappa, rho) = {sp1:+.3f} -> spearman(kappa, rho*(1-(a+bM4)/n)) = {sp2:+.3f}; max rel. error {relerr:.2e}; argmax agree: {int(np.argmax(kap)) == int(np.argmax(k2))}", flush=True)
