#!/usr/bin/env python3
"""Audit F2 (typing check for 'the start 5 pairs with a window of 7, the start 7 with a window of 5'):
for which a is the (a,1) Gibonacci sequence a window sum of r >= 2 consecutive Fibonacci (theta = +-phi^k (phi^r - 1))
or Lucas numbers (theta = +-sqrt5 phi^k (phi^r - 1))?  theta(a,1) = 1 + a(phi - 1) = a phi - (a - 1)."""
def zm(p, q): return (p[0]*q[0] + p[1]*q[1], p[0]*q[1] + p[1]*q[0] + p[1]*q[1])
def zp(p, k):
    r = (1, 0)
    for _ in range(k): r = zm(r, p)
    return r
PHI = (0, 1); PHIINV = (-1, 1); SQ5 = (-1, 2)
def phik(k): return zp(PHI, k) if k >= 0 else zp(PHIINV, -k)
reps = {}
for a in range(1, 41):
    th = (-(a - 1), a)
    for r in range(2, 40):
        base = phik(r); base = (base[0] - 1, base[1])
        for k in range(-45, 45):
            for kind, x in (('Fib', zm(phik(k), base)), ('Luc', zm(zm(phik(k), base), SQ5))):
                for s in (1, -1):
                    if (s * x[0], s * x[1]) == th:
                        reps.setdefault(a, []).append((kind, r, k, s))
for a in range(1, 41):
    print('a = %2d, N = -(a^2+a-1) = %5d :' % (a, -(a*a + a - 1)), reps.get(a, 'no window representation'))
