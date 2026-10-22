# Audit of HYP-9364 addendum (f): y = n(1 - kappa/rho) vs M4 = sum_k sigma_k^4/n^4 on circulant classes (n = 17..25).
import re, glob, io, contextlib
from fractions import Fraction
from math import factorial
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import d_circ_band as cb
for n in (17, 19, 21, 23, 25):
    Ds = list(cb.data[n]); H = [cb.data[n][D] for D in Ds]
    kap = np.array([h * 2 ** (n - 1) / factorial(n) for h in H])
    sig = [np.array([2 * sum(np.sin(2 * np.pi * j * x / n) for x in D) for j in range(1, (n - 1) // 2 + 1)]) for D in Ds]
    rho = np.array([np.prod((n * n + s * s) / (n * n - s * s)) for s in sig])
    M4 = np.array([np.sum(s ** 4) / n ** 4 for s in sig])
    y = n * (1 - kap / rho)
    b, a = np.polyfit(M4, y, 1)
    pred = rho * (1 - (a + b * M4) / n)
    print(f"n={n}: corr(y,M4)={np.corrcoef(y, M4)[0,1]:.3f}  fit a={a:.3f} b={b:.3f}  spearman(kappa, rho)={cb.spearman(list(kap), list(rho)):+.3f}  spearman(kappa, pred)={cb.spearman(list(kap), list(pred)):+.3f}  max rel err={np.max(np.abs(pred / kap - 1)):.2e}  argmax pred = argmax kappa: {np.argmax(pred) == np.argmax(kap)}")
d = 27
print("rho(R_27) closed form =", float(Fraction((d + 1) ** d + (d - 1) ** d) / Fraction(int((((d + 1j) ** d) + ((d - 1j) ** d)).real))) if False else None)
from mpmath import mp, mpf, mpc
mp.dps = 30
print("rho(R_27) =", ((mpf(d + 1) ** d + mpf(d - 1) ** d) / (mpc(d, 1) ** d + mpc(d, -1) ** d)).real, " (14/13)^13 =", mpf(14) ** 13 / mpf(13) ** 13)
