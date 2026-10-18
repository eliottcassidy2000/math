#!/usr/bin/env python3
"""Direct translation matchings between coefficient copies of clock_duality's cubic recipe (mac-mini-2026-10-18-unitpath).

Setting (golden_contact_circuit_20261017i.md section 4): F = Q(zeta15, sqrt-11, w) (degree 48, CM, conjugation c),
u^2 = a + b u + c Y with a = conj(d), d = 1 + zeta5, and the factors A + aM, C + bM, C + cM on axes 0, 1, 4.

Lemma (elementary, PROVED in THM-4700(3)-(4)): for kappa in F^* and units q != q' the copies C + kappa q and
C + kappa q' are joined by the translation matching x + kappa q ~ x + kappa q' (10 unit edges) iff kappa (q' - q) is a
unit, iff the unit r = q'/q is a root of T^2 - (2 - 1/N) T + 1 with N = kappa conj(kappa) (discriminant (1 - 4N)/N^2),
iff 4N - 1 lies in 3 (F+)^{*2} (r must be non-real: its square root of the discriminant lies in sqrt-3 F+).
Necessary condition tested here: 1 - 4N is a square in F; one quadratic non-residue excludes the matching.
This script computes N for a, b, c (and for the spindle-scaled kappa (m' - m)), screens 1 - 4N for squareness by
quadratic characters at split degree-1 places, and records the result. One non-residue at a place where the element is
a unit proves non-square.
"""
import os, sys, json, time
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__)); EXP = os.path.join(HERE, '..')
sys.path.insert(0, EXP); sys.path.insert(0, os.path.join(EXP, 'first_arrival_20261018', 'shim'))
from golden_forcing_search_20261017 import Geo
from golden_vsquare_unit_census_20261017 import places, red

T0 = time.time()
def log(*a): print('[%5.0f s]' % (time.time() - T0), *a, flush=True)

geo = Geo(); f = geo.f; bar = geo.bar
d = geo.d; db = bar(d)
w = f.gen(2); j = geo.j
m2 = f.scale(f.add(f.scalar(20), f.scale(j, 11)), Q(-1, 19))
m1 = f.scale(f.sub(f.scalar(20), f.scale(j, 11)), Q(1, 19))
db2 = f.mul(db, db)
cc = f.mul(f.add(f.sub(f.add(f.mul(db, m2), f.one), f.mul(db2, m1)), f.mul(db2, db)), f.inv(f.sub(f.mul(d, db), f.one)))
cb = f.add(f.mul(cc, d), f.sub(f.mul(db, m1), db2))
ca = db
# sanity: the cubic unit u and relation u^2 = a + b u + c Y
h = f.add(w, f.scalar(3)); u = f.mul(f.sub(h, j), f.inv(f.add(h, j)))
assert f.mul(u, bar(u)) == f.one
Z = f.add(d, u); Y = f.scale(f.mul(Z, f.inv(bar(Z))), -1)
assert f.mul(u, u) == f.add(ca, f.add(f.mul(cb, u), f.mul(cc, Y))), 'u^2 = a + b u + c Y'
log('field degree', f.n, '; u^2 = a + b u + c Y verified exactly')

def den(x):
    from math import lcm
    D = 1
    for c in x:
        if c: D = lcm(D, c.denominator)
    return D

def factor_small(n):
    out = []; p = 2
    while p * p <= n:
        while n % p == 0: out.append(p); n //= p
        p += 1
    if n > 1: out.append(n)
    return out

pl = places(f, 60)
log('screening places:', [p for p, _, _ in pl][:12], '... (%d places)' % len(pl))

def screen(x):
    """Return ('NON-square', p) at the first place where x is a unit non-residue, else ('residue at all', count)."""
    cnt = 0
    for p, emb, embc in pl:
        if any(c and c.denominator % p == 0 for c in x): continue
        r = red(x, p, emb)
        if r == 0: continue
        cnt += 1
        if pow(r, (p - 1) // 2, p) == p - 1: return ('NON-square', p)
    return ('residue at all', cnt)

res = {}
four = f.scalar(4)
for name, k in (('a = conj(d)', ca), ('b', cb), ('c', cc)):
    N = f.mul(k, bar(k)); X = f.sub(f.one, f.scale(N, 4))
    sX = screen(X); sm3X = screen(f.scale(X, -3))
    log('%-12s den(kappa) = %s = %s ; N = kappa conj(kappa) real: %s ; 1 - 4N: %s' % (name, den(k), factor_small(den(k)), bar(N) == N, sX))
    res[name] = {'den': den(k), 'den_factors': factor_small(den(k)), 'one_minus_4N': list(map(str, sX))}
# scaled by spindle differences: kappa (m' - m) for the 21 pairs of M -- the direct copy edge needs kappa (m' - m) unit
M = geo.M
log('spindle coefficient set M has', len(M), 'points')
pairs = []
for i in range(len(M)):
    for k in range(i + 1, len(M)):
        dm = f.sub(M[k], M[i])
        for name, kap in (('b', cb), ('c', cc)):
            x = f.mul(kap, dm)
            pairs.append((i, k, name, geo.is_unit(x), den(x)))
log('kappa (m\' - m) unit for any pair:', any(p[3] for p in pairs), '; denominators:', sorted({p[4] for p in pairs}))
res['pair_denominators'] = sorted({p[4] for p in pairs})
res['pair_unit'] = any(p[3] for p in pairs)
json.dump(res, open(os.path.join(HERE, 'matching_criterion.json'), 'w'), indent=1)
log('done')
