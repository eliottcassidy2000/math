#!/usr/bin/env python3
"""Exact checks for the 19-adic obstruction to connecting the coefficient copies of the cubic recipe (mac-mini-2026-10-18).
The tower order Z[zeta15][j][w] has discriminant supported on {2, 3, 5, 11, 37} (cubic disc = 592 = 2^4 * 37), so it is
19-maximal and an element with a 19 in a coordinate denominator is non-integral at some prime above 19. Checks:
  1. coordinate denominators of b and c (19 divides them);
  2. every point of C and every unit direction used by copy_connectors.py has 19-free denominators;
  3. N_{F/Q}(m - m') for all spindle pairs m != m' is prime to 19 (so b(m - m') is non-integral at the same prime).
Hence no path of unit steps along 19-integral directions joins C + b m to C + b m' (sums of 19-integral elements are
19-integral), and likewise for c."""
import os, sys, itertools
from fractions import Fraction as Q
from math import lcm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from golden_forcing_search_20261017 import Geo, recipe
geo = Geo(); f = geo.f; A, C = recipe(geo, 'base')
d, db = geo.d, geo.bar(geo.d)
m2 = f.scale(f.add(f.scalar(20), f.scale(geo.j, 11)), Q(-1, 19)); m1 = f.scale(f.sub(f.scalar(20), f.scale(geo.j, 11)), Q(1, 19))
db2 = f.mul(db, db)
cc = f.mul(f.add(f.sub(f.add(f.mul(db, m2), f.one), f.mul(db2, m1)), f.mul(db2, db)), f.inv(f.sub(f.mul(d, db), f.one)))
cb = f.add(f.mul(cc, d), f.sub(f.mul(db, m1), db2))
den = lambda x: lcm(*[c.denominator for c in x])
print('denominator of b: %d (19 | it: %s); of c: %d (19 | it: %s)' % (den(cb), den(cb) % 19 == 0, den(cc), den(cc) % 19 == 0))
EC = geo.unit_edges(C)
U0 = sorted({f.sub(C[b2], C[a2]) for a2, b2 in EC} | {f.sub(C[a2], C[b2]) for a2, b2 in EC})
z15 = f.gen(0); roots = [f.one]
for _ in range(14): roots.append(f.mul(roots[-1], z15))
roots = roots + [f.scale(r, -1) for r in roots]
U1 = sorted({f.mul(r, q) for r in roots for q in U0 + [f.one]})
print('points of C 19-integral:', all(den(x) % 19 for x in C), '; all %d directions 19-integral:' % len(U1), all(den(q) % 19 for q in U1))
# exact norm to Q via the regular representation (Gaussian elimination over Q)
def norm(x):
    n = f.n; cols = [f.mul(x, tuple(Q(1) if i == k else Q(0) for i in range(n))) for k in range(n)]
    Mx = [[cols[k][i] for k in range(n)] for i in range(n)]; det = Q(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if Mx[r][i] != 0), None)
        if piv is None: return Q(0)
        if piv != i: Mx[i], Mx[piv] = Mx[piv], Mx[i]; det = -det
        det *= Mx[i][i]; inv = 1 / Mx[i][i]
        for r in range(i + 1, n):
            if Mx[r][i] != 0:
                fac = Mx[r][i] * inv; Mx[r] = [a - fac * b for a, b in zip(Mx[r], Mx[i])]
    return det
bad = []
for (i, m), (k, mm) in itertools.combinations(enumerate(geo.M), 2):
    nm = norm(f.sub(mm, m))
    if nm.numerator % 19 == 0 or nm.denominator % 19 == 0: bad.append((i, k))
print('spindle differences m - m\' with 19 in N_{F/Q}:', bad if bad else 'none (all 21 pairs are 19-units)')
nb = norm(cb); nc = norm(cc)
print('v_19(N(b)) and v_19(N(c)) (negative = non-integral above 19):', end=' ')
def v19(q):
    v = 0; a, b = q.numerator, q.denominator
    while a % 19 == 0 and a: a //= 19; v += 1
    while b % 19 == 0: b //= 19; v -= 1
    return v
print(v19(nb), v19(nc))
