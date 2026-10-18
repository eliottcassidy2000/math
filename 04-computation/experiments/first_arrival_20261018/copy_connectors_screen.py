#!/usr/bin/env python3
"""Re-screen the three-step connector candidates of copy_connectors.py at many more degree-1 places, and try to certify
an exact middle point by solving s^2 = D over F with gp (nffactor) when the screen passes."""
import os, sys, json, itertools
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from golden_forcing_search_20261017 import Geo, recipe
from golden_vsquare_unit_census_20261017 import places, red
geo = Geo(); f = geo.f; A, C = recipe(geo, 'base')
d, db = geo.d, geo.bar(geo.d)
m2 = f.scale(f.add(f.scalar(20), f.scale(geo.j, 11)), Q(-1, 19)); m1 = f.scale(f.sub(f.scalar(20), f.scale(geo.j, 11)), Q(1, 19))
db2 = f.mul(db, db)
cc = f.mul(f.add(f.sub(f.add(f.mul(db, m2), f.one), f.mul(db2, m1)), f.mul(db2, db)), f.inv(f.sub(f.mul(d, db), f.one)))
cb = f.add(f.mul(cc, d), f.sub(f.mul(db, m1), db2))
EC = geo.unit_edges(C)
U0 = sorted({f.sub(C[b2], C[a2]) for a2, b2 in EC} | {f.sub(C[a2], C[b2]) for a2, b2 in EC})
SUFFIX = ''
if os.environ.get('ROOTS_OF_UNITY'):
    z15 = f.gen(0); roots = [f.one]
    for _ in range(14): roots.append(f.mul(roots[-1], z15))
    roots = roots + [f.scale(r, -1) for r in roots]
    U0 = sorted({f.mul(r, q) for r in roots for q in U0 + [f.one]}); SUFFIX = '_roots'
cands = json.load(open(os.path.join(HERE, 'copy_connectors%s.json' % SUFFIX)))['b']['three_step_survivors']
pl = places(f, int(sys.argv[1]) if len(sys.argv) > 1 else 80, start=1000)
print('places', len(pl), 'up to', pl[-1][0])
for (i1, i2, a1, a2, qi, nt) in cands:
    x = f.add(C[a1], f.mul(cb, geo.M[i1])); y = f.add(C[a2], f.mul(cb, geo.M[i2]))
    t = f.sub(f.sub(y, x), U0[qi]); Nt = f.mul(t, geo.bar(t)); X = f.mul(Nt, f.sub(Nt, f.scalar(4)))
    res = [pow(red(X, p, emb), (p - 1) // 2, p) for p, emb, embc in pl]
    print('candidate copies', (i1, i2), 'points', (a1, a2), 'q index', qi, '|t|^2 = %.6f' % nt,
          'residue at all places' if all(r in (0, 1) for r in res) else 'NON-square at place %d' % next(p for (p, _, _), r in zip(pl, res) if r not in (0, 1)))
