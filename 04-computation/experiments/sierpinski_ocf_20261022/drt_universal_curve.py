# Do all doubly regular tournaments lie on one universal curve kappa = e (1 - 3/(2n) + c2/n^2 + c3/n^3 + ...)?
# Data: Paley P_q (q = 3, 7, 11, 19, 23, 27, 31) and the non-Paley classes D_15, D_19, D_23 (drt_results.out).
# Fit c2, c3, c4 by least squares on z(n) = n^2 (kappa/e - 1 + 3/(2n)) against 1/n.  mac-mini-2026-10-22-sierp.
from mpmath import mp, mpf, e, factorial
import numpy as np
mp.dps = 30
data = [(3, 3, 'P'), (7, 189, 'P'), (11, 95095, 'P'), (15, 198335025, 'D'), (19, 1172695746915, 'P'), (19, 1172732961771, 'D'),
        (23, 15760206976379349, 'P'), (23, 15759953146384333, 'D'), (23, 15759901654866253, 'D'), (23, 15759902090984757, 'D'),
        (27, 418425945688463432031, 'P'), (31, 19876729650849384218015865, 'P')]
rows = []
for n, H, lab in data:
    k = mpf(H) * mpf(2) ** (n - 1) / factorial(n)
    z = n * n * (k / e - 1 + mpf(3) / (2 * n))
    rows.append((n, lab, float(k), float(z)))
    print(f"n={n:2d} {lab}: kappa = {mp.nstr(k, 12)}   n(1 - kappa/e) = {mp.nstr(n * (1 - k / e), 8)}   z = n^2(kappa/e - 1 + 3/(2n)) = {mp.nstr(z, 8)}")
sel = [r for r in rows if r[0] >= 11]
x = np.array([1 / r[0] for r in sel]); z = np.array([r[3] for r in sel])
for deg in (1, 2):
    c = np.polyfit(x, z, deg); res = z - np.polyval(c, x)
    print(f"fit z = c2 + c3/n{' + c4/n^2' if deg == 2 else ''} on n >= 11: c2 = {c[-1]:.4f}, c3 = {c[-2]:.4f}{', c4 = %.3f' % c[0] if deg == 2 else ''}; max |residual| = {np.max(np.abs(res)):.4f}")
