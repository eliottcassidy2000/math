#!/usr/bin/env python3
"""Pythagorean ternary tree and unit walks in the rational plane Q(i) (mac-mini-2026-10-18-unitpath).
(1) Berggren tree from (3,4,5) with A, B, C; every hypotenuse has all prime factors = 1 mod 4 (Fermat), so these are
    exactly the denominators of the unit vectors (a + b i)/c, and Gamma_(Q(i)) = Z[i][1/p : p = 1 mod 4] (B2 for Q(i):
    the conjugation-fixed primes are (1 + i) and the inert p = 3 mod 4).
(2) Two-step reciprocals never exist in Q(i): L(1/n) = 2 would need 4n^2 - 1 to be a square.
(3) The symmetric three-step walk 1/c = delta + conj(delta) -+ 1, delta = (a + b i)/c, needs 2a = c +- 1, i.e. a leg
    (c +- 1)/2: the near-60-degree triples. These are ((3c -+ 1)/2)^2 - 3 b^2 = 1 (Pell in Z[sqrt3]) and form the
    alternating branch C, A, C, A, ... of the tree: (3,4,5) -> (15,8,17) -> (33,56,65) -> (209,120,241) -> ...
(4) The 45-degree branch B, B, ...: legs differing by one, (2a +- 1)^2 - 2 c^2 = -1 type (Pell in Z[sqrt2]).
All checks exact (integers / Gaussian rationals)."""
from fractions import Fraction as Fr
from sympy import factorint
A = ((1, -2, 2), (2, -1, 2), (2, -2, 3)); B = ((1, 2, 2), (2, 1, 2), (2, 2, 3)); C = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
app = lambda M, t: tuple(sum(M[i][k] * t[k] for k in range(3)) for i in range(3))
root = (3, 4, 5)
# (1) BFS to depth 7
level = [root]; allt = [root]
for dep in range(7):
    level = [app(M, t) for t in level for M in (A, B, C)]; allt += level
ok1 = all(all(p % 4 == 1 for p in factorint(t[2])) and t[0] ** 2 + t[1] ** 2 == t[2] ** 2 for t in allt)
print('(1) %d triples to depth 7: all a^2 + b^2 = c^2 with every prime factor of c = 1 mod 4: %s' % (len(allt), ok1))
# (2)
from math import isqrt
print('(2) no n in 2..10^6 has 4n^2 - 1 a perfect square:', all(isqrt(4 * n * n - 1) ** 2 != 4 * n * n - 1 for n in range(2, 10 ** 6)))
# (3) the CACA branch
t = root; branch = [t]
for k in range(8):
    t = app(C if k % 2 == 0 else A, t); branch.append(t)
def walk3(tr):
    a, b, c = tr
    for leg in (a, b):
        for sgn in (1, -1):
            if 2 * leg - sgn == c or 2 * leg + sgn == c:
                # delta + conj delta = 2 leg / c ; 1/c = (delta + conj delta) - 1 if 2 leg - 1 = c ; = 1 - (...) if 2 leg + 1 = c
                return True
    return False
pell = []
for tr in branch:
    a, b, c = sorted(tr[:2]) + [tr[2]]
    leg = a if abs(2 * a - c) == 1 else b; other = b if leg == a else a
    X2 = 3 * c - (1 if 2 * leg - c == 1 else -1)        # 3c -+ 1 = 2X
    pell.append((c, X2 // 2, other, (X2 // 2) ** 2 - 3 * other ** 2))
print('(3) branch C,A,C,A,...:', [tr for tr in branch[:6]])
print('    hypotenuses', [p[0] for p in pell], '; each has a leg (c +- 1)/2 (symmetric 3-step walk):', all(walk3(tr) for tr in branch))
print('    Pell data (c, X, Y = other leg, X^2 - 3Y^2):', pell[:6])
# exact Gaussian-rational check of the walks
def gmul(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
okw = True
for tr in branch[:7]:
    a, b, c = tr
    for leg, oth in ((a, b), (b, a)):
        if abs(2 * leg - c) == 1:
            d = (Fr(leg, c), Fr(oth, c)); dc = (d[0], -d[1])
            assert d[0] ** 2 + d[1] ** 2 == 1
            s = (d[0] + dc[0], d[1] + dc[1])
            val = (s[0] - 1, s[1]) if 2 * leg - c == 1 else (1 - s[0], -s[1])
            okw &= val == (Fr(1, c), 0)
print('    exact: 1/c = delta + conj(delta) - 1 (or 1 - delta - conj(delta)) for the first seven branch triples:', okw)
# is every near-60 triple on this branch? search all triples to depth 9 for a leg (c +- 1)/2
level = [root]; near60 = [root]
for dep in range(9):
    level = [app(M, t) for t in level for M in (A, B, C)]
    near60 += [t for t in level if abs(2 * t[0] - t[2]) == 1 or abs(2 * t[1] - t[2]) == 1]
print('    near-60 triples (leg = (c +- 1)/2) among all triples to depth 9:', sorted(t[2] for t in near60), '; branch hypotenuses', sorted(set(p[0] for p in pell)))
# (4) the BBB branch
t = root; bb = [t]
for k in range(6): t = app(B, t); bb.append(t)
print('(4) branch B,B,B,...:', bb[:5], '; legs differ by one:', all(abs(x[0] - x[1]) == 1 for x in bb))
