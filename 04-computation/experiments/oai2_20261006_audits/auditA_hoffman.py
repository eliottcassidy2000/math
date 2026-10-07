"""Audit A: Hoffman ratio bound for kappa(p) = chi(Cay(F_{p^2}, mu_{p+1})), p odd prime.
F_{p^2} = F_p(t), t^2 = n (n a non-residue). Element a + b t; norm a^2 - n b^2. mu = norm-1 circle (p+1 points).
Eigenvalues: lambda_c = sum_{u in mu} exp(2 pi i Tr(c u)/p), Tr(x + y t) = 2x. c ranges over F_{p^2} \ 0.
Computed with numpy (double precision) AND re-verified exactly-ish with mpmath at 50 digits for the minimising c."""
import numpy as np, math, mpmath
from sympy import primerange, legendre_symbol
mpmath.mp.dps = 50
rows = []
for p in primerange(3, 140):
    n = next(a for a in range(2, p) if legendre_symbol(a, p) == -1)
    circ = np.array([(a, b) for a in range(p) for b in range(p) if (a * a - n * b * b) % p == 1])
    assert len(circ) == p + 1
    A, B = circ[:, 0], circ[:, 1]
    best = (1e9, None)
    for c in range(p):
        d = np.arange(p)
        # c*u = (c0 + c1 t)(a + b t) = (c0 a + n c1 b) + (...) t ; Tr = 2 (c0 a + n c1 b)
        X = (2 * (c * A[None, :] + n * d[:, None] * B[None, :])) % p
        lam = np.cos(2 * np.pi * X / p).sum(axis=1)
        if c == 0: lam[0] = 1e9
        i = int(np.argmin(lam))
        if lam[i] < best[0]: best = (float(lam[i]), (c, int(d[i])))
    c0, c1 = best[1]
    exact = mpmath.fsum(mpmath.cos(2 * mpmath.pi * ((2 * (c0 * int(a) + n * c1 * int(b))) % p) / p) for a, b in circ)
    deg = p + 1
    hof = 1 + deg / (-exact)
    rows.append((p, float(exact), float(hof)))
    print(f"p={p:3d} lambda_min={float(exact):10.6f}  -2sqrt(p)={-2*math.sqrt(p):8.3f}  Hoffman 1+d/|l|={float(hof):.6f}  => kappa >= {math.ceil(float(hof) - 1e-12)}", flush=True)
print()
print("kappa >= 5 (Hoffman) for p:", [p for p, l, h in rows if h > 4])
print("kappa >= 6 (Hoffman) for p:", [p for p, l, h in rows if h > 5])
print("margin at p=31:", [h - 4 for p, l, h in rows if p == 31])
