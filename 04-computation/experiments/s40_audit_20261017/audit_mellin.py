"""S40 audit: Gamma(s) L(s, chi_-p) = int_0^oo u^(s-1) f_p(e^-u) / (1 - e^(-p u)) du (THM-4666 (5)), double precision.
L(s, chi) from the Hurwitz-zeta form p^(-s) sum_a chi(a) zeta(s, a/p) (mpmath, 20 digits); integral by scipy quad."""
import math
import numpy as np
import mpmath as mp
from scipy.integrate import quad
mp.mp.dps = 20
for p in (23, 43, 47, 163):
    chi = [0] + [1 if pow(a, (p - 1) // 2, p) == 1 else -1 for a in range(1, p)]
    coeffs = np.array(chi, dtype=float)
    def g(u, s):
        x = math.exp(-u)
        if u < 1e-6:
            # f_p(e^-u)/(1-e^-pu) -> -f_p'(1)/p as u -> 0 (f_p(1) = 0)
            val = -sum(t * chi[t] for t in range(p)) / p
        else:
            val = np.polyval(coeffs[::-1], x) / (1 - math.exp(-p * u))
        return u ** (s - 1) * val
    for s in (0.5, 0.8, 1.0, 2.0):
        L = p ** (-s) * mp.fsum(chi[a] * mp.zeta(s, mp.mpf(a) / p) for a in range(1, p))
        lhs = float(mp.gamma(s) * L)
        rhs = quad(g, 0, 1, args=(s,), limit=400)[0] + quad(g, 1, 60, args=(s,), limit=400)[0]
        print('p=%d s=%.1f Gamma(s)L(s,chi) = %.10f  integral = %.10f  diff %.1e' % (p, s, lhs, rhs, abs(lhs - rhs)))
    xs = np.linspace(1e-4, 1 - 1e-4, 20001)
    vals = np.polyval(coeffs[::-1], xs)
    print('   p=%d: min of f_p on (0,1) grid = %.5f (sign change: %s)' % (p, vals.min(), bool((vals < 0).any())))
