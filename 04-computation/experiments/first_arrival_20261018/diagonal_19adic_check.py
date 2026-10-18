#!/usr/bin/env python3
"""Can DIAGONAL moves across the two 19-adic coefficient factors of the cubic recipe be 19-integral? (mac-mini-2026-10-18)
A point of the recipe sum is alpha + a m1 + gamma2 + b m2 + gamma3 + c m3 (alpha in A, gamma_i in C, m_i in the spindle M,
a = conj(d)). If A, C and a are 19-integral, a path of 19-integral unit steps can change the sheet index (m2, m3) by
(dm2, dm3) only if b dm2 + c dm3 is 19-integral (integral at every prime above 19). THM-4690(5) is the case dm3 = 0 (or
dm2 = 0). Here every pair (dm2, dm3) in (M - M)^2 is tested; the tower order is 19-maximal, so 19-integrality is the
absence of 19 from the coordinate denominators."""
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
i19 = lambda x: den(x) % 19 != 0
print('19-integral: a = conj(d): %s ; all points of A: %s ; all points of C: %s ; spindle points: %s' % (i19(db), all(i19(x) for x in A), all(i19(x) for x in C), all(i19(x) for x in geo.M)))
D = sorted({f.sub(y, x) for x in geo.M for y in geo.M})
print('|M - M| =', len(D))
zero = f.sub(f.one, f.one)
hits = []
for x2 in D:
    bx = f.mul(cb, x2)
    for x3 in D:
        if x2 == zero and x3 == zero: continue
        if i19(f.add(bx, f.mul(cc, x3))): hits.append((D.index(x2), D.index(x3)))
print('nonzero (dm2, dm3) with b dm2 + c dm3 19-integral:', len(hits), hits[:20])
# ratio b / c: if it were a 19-unit times a spindle-difference ratio, diagonal resonances would exist
r = f.mul(cb, f.inv(cc))
print('b / c 19-integral: %s ; c / b 19-integral: %s ; denominators %d, %d' % (i19(r), i19(f.inv(r)), den(r), den(f.inv(r))))
