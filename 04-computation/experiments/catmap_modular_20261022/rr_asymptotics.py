"""Growth constants and q -> 1 ratios for the Rogers-Ramanujan-type product sides.  Session opus-2026-10-22-S55.
c_eff from part densities (exact); numerical check from log f(exp(-t)) ~ pi^2 c_eff/(6t); q -> 1 ratios of companion
products compared with sin(pi a/m) ratios (quantum dimensions at q = exp(i pi/m)).
Usage: python rr_asymptotics.py"""
import math
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 30


def logprod(t, factors):
    """log of prod over factors (power, sign, a, m) = ((sign) q^a; q^m)_inf^power at q = exp(-t), by direct summation."""
    import numpy as np
    t = float(t)
    tot = 0.0
    for power, sign, a, m in factors:
        K = int(60 / (t * m)) + 10
        n = a + m * np.arange(K, dtype=float)
        tot += power * np.sum(np.log1p(sign * np.exp(-t * n)))
    return tot


def ceff(factors):
    c = Fraction(0)
    for power, sign, a, m in factors:
        d = Fraction(1, m)
        c += (-power) * d if sign == -1 else power * d / 2            # (q^a;q^m)^power: power*(-d)... see below
    return c


# a factor (power, sign, a, m) means ((sign) q^a; q^m)_inf ^ power with sign = -1 for (q^a;q^m)_inf, +1 for (-q^a;q^m)_inf
# log (q^a;q^m)_inf ~ -(1/m) pi^2/(6t); log (-q^a;q^m)_inf ~ +(1/(2m)) pi^2/(6t)
def ceff2(factors):
    c = Fraction(0)
    for power, sign, a, m in factors:
        c += power * (Fraction(-1, m) if sign == -1 else Fraction(1, 2 * m))
    return c


def triple(a, m):
    return [(1, -1, a, m), (1, -1, m - a, m), (1, -1, m, m)]


fam = {
    'RR (mod 5)': {'G': triple(2, 5) + [(-1, -1, 1, 1)], 'H': triple(1, 5) + [(-1, -1, 1, 1)]},
    'Rogers-Selberg (mod 7)': {'A': triple(3, 7) + [(-1, -1, 2, 2)], 'B': triple(2, 7) + [(-1, -1, 2, 2)], 'C': triple(1, 7) + [(-1, -1, 2, 2)]},
    'Rogers (mod 14)': {'A': triple(6, 14) + [(-1, -1, 1, 1)], 'B': triple(4, 14) + [(-1, -1, 1, 1)], 'C': triple(2, 14) + [(-1, -1, 1, 1)]},
    'Goellnitz-Gordon (mod 8)': {'1': [(-1, -1, 1, 8), (-1, -1, 4, 8), (-1, -1, 7, 8)], '2': [(-1, -1, 3, 8), (-1, -1, 4, 8), (-1, -1, 5, 8)]},
    'Bailey (mod 9)': {'A': triple(4, 9) + [(-1, -1, 3, 3)], 'B': triple(2, 9) + [(-1, -1, 3, 3)], 'C': triple(1, 9) + [(-1, -1, 3, 3)]},
    'Dyson (mod 27)': {'A': triple(12, 27) + [(-1, -1, 1, 1)], 'B': [(1, -1, 9, 9), (-1, -1, 1, 1)], 'C': triple(6, 27) + [(-1, -1, 1, 1)], 'D': triple(3, 27) + [(-1, -1, 1, 1)]},
    'Jackson-Slater (Ising vacuum)': {'JS': [(1, -1, 1, 8), (1, -1, 7, 8), (1, -1, 8, 8), (1, -1, 6, 16), (1, -1, 10, 16), (-1, -1, 1, 1)]},
}
for name, fs in fam.items():
    print(name)
    base = None
    for lab, f in fs.items():
        c = ceff2(f)
        t1, t2 = 0.004, 0.002
        l1, l2 = logprod(t1, f), logprod(t2, f)
        cnum = (l2 - l1) / (math.pi ** 2 / 6 * (1 / t2 - 1 / t1))
        val = math.exp(logprod(0.0005, f) - logprod(0.0005, list(fs.values())[0]))
        print(f'   {lab:3s}: c_eff = {str(c):6s} (numeric {float(cnum):+.4f}); 6/c_eff = {(str(6 / c) if c else "inf"):6s}; ratio to first at t = 5e-4: {val:.6f}')
phi = (1 + 5 ** 0.5) / 2
print('reference ratios: 1/phi =', 1 / phi, '; sin(k pi/7)/sin(3 pi/7):', [round(math.sin(k * math.pi / 7) / math.sin(3 * math.pi / 7), 6) for k in (2, 1)],
      '; sin(k pi/9)/sin(4 pi/9):', [round(math.sin(k * math.pi / 9) / math.sin(4 * math.pi / 9), 6) for k in (2, 1)],
      '; sin(k pi/27)/sin(12 pi/27):', [round(math.sin(k * math.pi / 27) / math.sin(12 * math.pi / 27), 6) for k in (9, 6, 3)],
      '; sin(k pi/14)/sin(6 pi/14):', [round(math.sin(k * math.pi / 14) / math.sin(6 * math.pi / 14), 6) for k in (4, 2)])
