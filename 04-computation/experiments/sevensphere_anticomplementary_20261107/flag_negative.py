#!/usr/bin/env python3
# THM-4910 statement 6 (mac-mini-2026-11-07-sevensphere): NEGATIVE test of a flag extension of the corner-power rule.
# For a random n-simplex and the flag F_n > F_(n-1) > ... > F_2 (F_k = conv(V_(n-k), ..., V_n)), ask for one sphere with
# Pow(O_k) = (k+1)(R^2 - |OO_k|^2) and Pow(G_k) = k(R^2 - |OG_k|^2) for all k. Pow(x) - |x|^2 must then be affine; we fit
# it by least squares on the 2(n-1) corners and report the worst relative residual (0 = consistent). Output: flag_negative.out
import numpy as np
rng = np.random.default_rng(11)
def circumcenter(P):
    p0 = P[0]; Dm = P[1:] - p0
    lam = np.linalg.solve(Dm @ Dm.T, 0.5 * (Dm**2).sum(1)); return p0 + lam @ Dm
for n in (3, 4, 5, 6):
    worst = 0.0
    for trial in range(200):
        V = rng.normal(size=(n + 1, n)); O = circumcenter(V); R2 = ((V[0] - O)**2).sum()
        pts, vals = [], []
        for k in range(n, 1, -1):
            F = V[n - k:]; Ok = circumcenter(F); Gk = F.mean(0)
            pts += [Ok, Gk]; vals += [(k + 1) * (R2 - ((Ok - O)**2).sum()), k * (R2 - ((Gk - O)**2).sum())]
        P = np.array(pts); y = np.array(vals) - (P**2).sum(1)
        X = np.hstack([-2 * P, np.ones((len(P), 1))])
        sol, *_ = np.linalg.lstsq(X, y, rcond=None)
        worst = max(worst, np.abs(X @ sol - y).max() / np.abs(y).max())
    print(f"n={n}: {2*(n-1)} flag corners; worst relative residual of an affine fit over 200 simplices: {worst:.2e}"
          + ("  (n = 3: four corners always fit; this is THM-4910(1))" if n == 3 else "  (no common sphere)"))
