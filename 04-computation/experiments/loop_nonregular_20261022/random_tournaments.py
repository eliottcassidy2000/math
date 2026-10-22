"""Loop expansion on uniformly random tournaments: the tree level pays back the one-loop gain.
Session opus-2026-10-22-S54.

For every tournament sum_k x_k^2 = -tr Shat^2 = (n - 1)/n, so the one-loop gain is rho ~ e^(1 - 1/n) (times higher
cycle terms), while the path-tree term X(-1,2) = -delta^T (I + Shat^2)^-1 delta ~ -|delta|^2 with
|delta|^2 = sum_u ((2 s_u - n + 1)/n)^2, of mean (n - 1)/n for a uniformly random tournament (and 0 for regular ones).
So the leading loop expansion predicts log kappa ~ log rho + X(-1,2): regular tournaments keep the factor e, random ones
pay it back through score irregularity, consistent with E[kappa] = 1.
Prints means, correlations and the first/second-order log-form predictions against exact H.
Usage: python random_tournaments.py n FILE"""
import sys, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leaf_engine as LE
import leaf_closed as LC
from predict_nonregular import cumulants, X1_closed, tree_closed, load

n = int(sys.argv[1]); path = sys.argv[2]
rows = []
for code, H, S in load(n, path):
    Sh = S / n
    x = np.linalg.eigvals(Sh).imag
    logrho = 0.5 * np.sum(np.log((1 + x ** 2) / (1 - x ** 2)))
    kap = H * 2 ** (n - 1) / math.factorial(n)
    d = Sh @ np.ones(n)
    k = cumulants(lambda z: 0.5 * np.sum(np.log((1 + (z * x) ** 2) / (1 - (z * x) ** 2))))
    k1, k2 = k[1], k[2]
    Ea = (k2 + k1 ** 2 - k1) / 2
    t2 = LC.X_tree2(Sh); x10 = X1_closed(Sh, 1.0).real
    rows.append(dict(logk=math.log(kap) if kap > 0 else -np.inf, logrho=logrho, d2=float(d @ d), t2=t2, x10=x10,
                     lol=LC.X_lollipop(Sh), star=LC.X_star(Sh), x02=LC.X_02(Sh), Ea=Ea, maxx=float(np.max(np.abs(x)))))
R = {k: np.array([r[k] for r in rows]) for k in rows[0]}
ok = np.isfinite(R['logk']) & (R['maxx'] < 0.9)
print(f'n={n}: {len(rows)} uniformly random tournaments ({ok.sum()} with H > 0 and spectral radius < 0.9)')
print(f'  mean kappa = {np.mean(np.exp(R["logk"][np.isfinite(R["logk"])])):.4f} (exact expectation 1); mean log rho = {R["logrho"].mean():.4f}; '
      f'mean |delta|^2 = {R["d2"].mean():.4f} (expected (n-1)/n = {(n - 1) / n:.4f}); mean X(-1,2) = {R["t2"].mean():+.4f}')
y = R['logk'][ok] - R['logrho'][ok]
for label, pred in (('X(-1,2)', R['t2']), ('X(-1,2) + X(1,0) + E0[a]/n', R['t2'] + R['x10'] + R['Ea'] / n),
                    ('+ second-order leaves (lollipop, X(0,2)) + star', R['t2'] + R['x10'] + R['Ea'] / n + R['lol'] + R['x02'] + R['star'])):
    p = pred[ok]
    print(f'  log(kappa/rho) vs {label:45s}: corr {np.corrcoef(y, p)[0, 1]:+.4f}, mean diff {np.mean(y - p):+.4f}, sd diff {np.std(y - p):.4f} (sd of log(kappa/rho) {np.std(y):.4f})')
print(f'  corr(log kappa, log rho) = {np.corrcoef(R["logk"][ok], R["logrho"][ok])[0, 1]:+.4f}; corr(log kappa, -|delta|^2) = '
      f'{np.corrcoef(R["logk"][ok], -R["d2"][ok])[0, 1]:+.4f}; corr(log kappa, log rho + X(-1,2)) = {np.corrcoef(R["logk"][ok], R["logrho"][ok] + R["t2"][ok])[0, 1]:+.4f}')
