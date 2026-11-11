#!/usr/bin/env python3
"""The mod-4 race in the free monoid of tournaments.
Every tournament is uniquely an ordered sum (transitive concatenation) of strong tournaments, H is multiplicative over
the components, and chi(H) = (-1)^((H-1)/2) is a character mod 4.  With labelled exponential generating functions
    T(x) = 1/(1 - S(x)),   T_chi(x) = 1/(1 - S_chi(x)),
so the signed counts over STRONG labelled tournaments follow from the signed counts over all labelled tournaments:
    s_chi(n) = t_chi(n) - sum_{k=1}^{n-1} C(n,k) s_chi(k) t_chi(n-k).
Input: t_chi(n) = 2^C(n,2) D_n for n <= 11 (THM-5127, FINITE-EXACT)."""
from fractions import Fraction
from math import comb
tchi = {0: 1, 1: 1, 2: 2, 3: 4, 4: 32, 5: 416, 6: 8192, 7: 218368, 8: 11149312, 9: 516939776, 10: 85258534912, 11: 17075203538944}
N = 11
t = {n: 2 ** (n * (n - 1) // 2) for n in range(N + 1)}
def invert(tt):
    s = {}
    for n in range(1, N + 1):
        s[n] = tt[n] - sum(comb(n, k) * s[k] * tt[n - k] for k in range(1, n))
    return s
s = invert(t); schi = invert(tchi)
print(" n   strong labelled s(n)        signed strong s_chi(n)      D_n (all)        D_n^str (strong)   D^str/D    share of D_n from decomposable")
for n in range(1, N + 1):
    D = Fraction(tchi[n], t[n]); Ds = Fraction(schi[n], s[n]) if s[n] else None
    dec = Fraction(tchi[n] - schi[n], tchi[n])
    print("%2d   %-26d  %-26d  %.9f      %s   %s    %.6f" % (n, s[n], schi[n], float(D), ("%.9f" % float(Ds)) if Ds is not None else "    -      ",
          ("%.5f" % float(Ds / D)) if Ds is not None else "   -   ", float(dec)))
print()
print("scaled: D_n * 4^n and D_n^str * 4^n")
for n in range(3, N + 1):
    print("%2d  %10.2f  %10.2f" % (n, float(Fraction(tchi[n], t[n])) * 4 ** n, float(Fraction(schi[n], s[n])) * 4 ** n))
# decomposition of t_chi(n) by the number of strong components
print()
print("t_chi(n) by number of strong components r (labelled):")
for n in range(3, N + 1):
    # coefficient extraction of S_chi(x)^r
    from functools import lru_cache
    @lru_cache(None)
    def part(m, r):
        if r == 0: return 1 if m == 0 else 0
        return sum(comb(m, k) * schi[k] * part(m - k, r - 1) for k in range(1, m - r + 2))
    row = [part(n, r) for r in range(1, n + 1)]
    assert sum(row) == tchi[n]
    print("%2d " % n, [ "%.3g" % (v / tchi[n]) for v in row[:5]])
