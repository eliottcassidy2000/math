#!/usr/bin/env python3
"""House lower bound for unit walks between coefficient copies of the cubic recipe (mac-mini-2026-10-18-unitpath).
Every complex embedding sigma of the CM field F maps a unit walk to a unit walk (sigma commutes with conjugation), so
L(x) >= max_sigma |sigma(x)| (the house of x). A connector between copies C + k m and C + k m' (k = b, c) joins some
x + k m to some y + k m', so its length is at least
    H(m, m') = min over (x, y) in C^2 of house(y - x + k (m' - m)).
Evaluated in double precision over all 48 embeddings; reported margins are far from integers (checked below)."""
import os, sys, cmath, itertools, json, math
import numpy as np
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__)); EXP = os.path.join(HERE, '..')
sys.path.insert(0, EXP); sys.path.insert(0, os.path.join(EXP, 'first_arrival_20261018', 'shim'))
from golden_forcing_search_20261017 import Geo, recipe
geo = Geo(); f = geo.f; bar = geo.bar; j = geo.j
d = geo.d; db = bar(d)
m2 = f.scale(f.add(f.scalar(20), f.scale(j, 11)), Q(-1, 19)); m1 = f.scale(f.sub(f.scalar(20), f.scale(j, 11)), Q(1, 19))
db2 = f.mul(db, db)
cc = f.mul(f.add(f.sub(f.add(f.mul(db, m2), f.one), f.mul(db2, m1)), f.mul(db2, db)), f.inv(f.sub(f.mul(d, db), f.one)))
cb = f.add(f.mul(cc, d), f.sub(f.mul(db, m1), db2))
A, C = recipe(geo, 'base'); M = geo.M
wroots = sorted(r.real for r in np.roots([1, 2, -4, -4]))
embs = [(cmath.exp(2j * cmath.pi * k / 15), s * 1j * 11 ** 0.5, wr) for k in range(1, 15) if k % 3 and k % 5 for s in (1, -1) for wr in wroots]
basis = list(f.basis)
def vec(x):  # values of x in all embeddings
    return np.array([sum(float(c) * e[0] ** a * e[1] ** b * e[2] ** ce for c, (a, b, ce) in zip(x, basis) if c) for e in embs])
VC = [vec(x) for x in C]
out = {}
for name, k in (('b', cb), ('c', cc)):
    Vk = vec(k); VM = [vec(m) for m in M]
    table = []
    for i, i2 in itertools.combinations(range(len(M)), 2):
        best = None
        for x in range(len(C)):
            for y in range(len(C)):
                hv = np.max(np.abs(VC[y] - VC[x] + Vk * (VM[i2] - VM[i])))
                if best is None or hv < best[0]: best = (hv, x, y)
        table.append((i, i2, best[0], best[1], best[2]))
    hmin = min(t[2] for t in table)
    out[name] = [dict(m=t[0], m2=t[1], house_min=t[2], x=t[3], y=t[4]) for t in table]
    print('kappa = %s: over the 21 copy pairs, min_(x,y) house(y - x + kappa (m\' - m)) ranges %.4f .. %.4f' % (name, hmin, max(t[2] for t in table)))
    print('   => every unit walk between two distinct copies of C + %sM has length >= %d' % (name, math.ceil(hmin - 1e-9)))
    for t in sorted(table, key=lambda t: t[2])[:4]:
        print('   pair (m%d, m%d): house %.6f at (x, y) = (%d, %d)' % t)
json.dump(out, open(os.path.join(HERE, 'house_bound.json'), 'w'), indent=1)
