#!/usr/bin/env python3
"""Loneliness of the Fibonacci leapers (chessboard-weave, 2026-10-06).

For speeds a, b (the (a,b)-leaper read as two runners) put
    delta(a,b) = max_t min(||t a||, ||t b||),   ||.|| = distance to the nearest integer.
Claim: for coprime 1 <= a <= b, delta(a,b) = floor((a+b)/2)/(a+b), attained at t = k/(a+b).
Exact check for (a,b) = (F_n, F_{n+1}), n = 0..30:
  * over the candidate times k/(a+b) (the times the line t(a,b) mod 1 meets u+v = 0 mod 1);
  * over ALL breakpoints of the piecewise-linear function min(||ta||,||tb||) on [0,1]:
    m/(2a), m/(2b) (kinks of each tent) and k/(a+b), k/(b-a) (crossings ||ta|| = ||tb||),
    so the maximum over this finite set is the true maximum.
Integer arithmetic only (t = k/D: ||t a|| = min(ka mod D, D - ka mod D)/D).
Run:  python3 chessboard_weave_20261006_fib_lonely.py
"""
from fractions import Fraction
import numpy as np


def fib(k):
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def family_max(a, b, D):
    """max over t = k/D, k = 0..D, of min(||ta||, ||tb||); exact, returned as Fraction."""
    k = np.arange(0, D + 1, dtype=np.int64)
    ra = (k * a) % D
    rb = (k * b) % D
    m = np.minimum(np.minimum(ra, D - ra), np.minimum(rb, D - rb))
    i = int(m.argmax())
    return Fraction(int(m[i]), D), i


names = {(0, 1): 'rook', (1, 1): 'bishop', (1, 2): 'knight', (2, 3): 'zebra'}
print(" n   (F_n, F_n+1)        leaper    delta (exact)           = floor((a+b)/2)/(a+b)  attained at k/(a+b)  colour-pres.  delta = 1/2")
allok = True
for n in range(0, 31):
    a, b = fib(n), fib(n + 1)
    dens = [D for D in {2 * a, 2 * b, a + b, b - a} if D > 0]
    best = max(family_max(a, b, D)[0] for D in dens)
    cand, kstar = family_max(a, b, a + b)
    formula = Fraction((a + b) // 2, a + b)
    colour_pres = (a + b) % 2 == 0
    ok = best == formula == cand and ((best == Fraction(1, 2)) == colour_pres)
    allok &= ok
    print("%2d (%8d,%8d) %-8s %-22s %-6s %-20s %-13s %s"
          % (n, a, b, names.get((a, b), ''), str(best), best == formula, "%s/%d" % (kstar, a + b) if n else '-',
             colour_pres, best == Fraction(1, 2)))
print("all n = 0..30: global max over all breakpoints == max over k/(a+b) == floor((a+b)/2)/(a+b), "
      "and delta = 1/2 exactly iff a+b even (colour-preserving):", allok)
print("for a+b odd: delta = 1/2 - 1/(2(a+b)); e.g. knight 1/3, zebra 2/5, (5,8): 6/13, (8,13): 10/21")
