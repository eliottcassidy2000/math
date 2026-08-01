"""Is 'nu is a positive measure on [0,inf)' ENOUGH to force Phi_n squarefree?

Phi_n(z) = int (z-u)^n dnu(u).  Multiple root  <=>  Phi_n(z0)=Phi_{n-1}(z0)=0.

TWO ATOMS: impossible.  p1 r^n = -p2 and p1 r^{n-1} = -p2 with r=(z-u1)/(z-u2)
force r=1, i.e. u1=u2.  PROVED, no computation needed.

THREE ATOMS: search.  nu = p1 d_{u1} + p2 d_{u2} + p3 d_{u3}, all p_i > 0.
"""
import numpy as np
from scipy.optimize import fsolve

def eqs(v, n):
    a, b, c, x, y = v
    u = np.array([0.0, 1.0, a]); p = np.array([1.0, b, c])
    z = complex(x, y)
    F1 = np.sum(p * (z - u) ** n)
    F2 = np.sum(p * (z - u) ** (n - 1))
    return [F1.real, F1.imag, F2.real, F2.imag, 0.0]   # 4 eqs, 5 unknowns -> family

found = []
rng = np.random.default_rng(20260731)
for n in (3, 4, 5, 6):
    hits = 0
    for _ in range(4000):
        g = [rng.uniform(1.5, 6), rng.uniform(.1, 5), rng.uniform(.1, 5),
             rng.uniform(-3, 5), rng.uniform(.1, 4)]
        try:
            sol, info, ier, _ = fsolve(eqs, g, args=(n,), full_output=True, xtol=1e-13)
        except Exception:
            continue
        if ier != 1:
            continue
        a, b, c, x, y = sol
        if not (a > 1.05 and b > 1e-3 and c > 1e-3 and abs(y) > 1e-4):
            continue
        r = np.array(eqs(sol, n))
        if np.max(np.abs(r)) < 1e-9:
            hits += 1
            if len(found) < 4:
                found.append((n, a, b, c, x, y, float(np.max(np.abs(r)))))
    print(f"  n={n}: {hits} valid positive-weight solutions found")

print("\nWITNESSES (u = (0,1,a), p = (1,b,c) all > 0, z0 = x+iy a MULTIPLE root of Phi_n):")
for (n, a, b, c, x, y, res) in found:
    u = np.array([0.0, 1.0, a]); p = np.array([1.0, b, c]); z = complex(x, y)
    print(f"  n={n}  u=(0, 1, {a:.6f})  p=(1, {b:.6f}, {c:.6f})  z0={z:.6f}")
    for k in (n, n - 1):
        print(f"      Phi_{k}(z0) = {np.sum(p*(z-u)**k):+.3e}")
