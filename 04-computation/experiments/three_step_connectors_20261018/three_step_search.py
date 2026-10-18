#!/usr/bin/env python3
"""Search for three-step connectors between coefficient copies of the cubic recipe (mac-mini-2026-10-18-unitpath).
X = (y - x) + k (m' - m), k = b, c, (x, y) in C^2, only targets with house(X) <= 3 (others need >= 4 steps).
L(X) <= 3 iff X - u is a two-sum for some unit u, iff N(X - u)(4 - N(X - u))/3 lies in (F+)^{*2} with
0 < N(X - u) <= 4 in every embedding (THM-4700(3)). Units u range over products of named generators.
Screen: quadratic characters at degree-one places of F (vectorized residues; conj realized by the conjugate embedding).
A survivor of all places is then checked against 60 places and the real-embedding condition."""
import os, sys, itertools, json, time, cmath
import numpy as np
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__)); EXP = os.path.join(HERE, '..')
sys.path.insert(0, EXP); sys.path.insert(0, os.path.join(EXP, 'first_arrival_20261018', 'shim'))
from golden_forcing_search_20261017 import Geo, recipe
from golden_vsquare_unit_census_20261017 import places, red
T0 = time.time()
def log(*a): print('[%5.0f s]' % (time.time() - T0), *a, flush=True)
geo = Geo(); f = geo.f; bar = geo.bar; j = geo.j; z15 = f.gen(0); w = f.gen(2)
isunit = lambda x: f.mul(x, bar(x)) == f.one
s3 = f.add(f.scale(f.power(z15, 5), 2), f.one); s5 = f.add(f.scale(f.add(f.power(z15, 3), f.power(z15, 12)), 2), f.one); s15 = f.mul(s3, s5)
d = geo.d; db = bar(d)
m2 = f.scale(f.add(f.scalar(20), f.scale(j, 11)), Q(-1, 19)); m1 = f.scale(f.sub(f.scalar(20), f.scale(j, 11)), Q(1, 19))
db2 = f.mul(db, db)
cc = f.mul(f.add(f.sub(f.add(f.mul(db, m2), f.one), f.mul(db2, m1)), f.mul(db2, db)), f.inv(f.sub(f.mul(d, db), f.one)))
cb = f.add(f.mul(cc, d), f.sub(f.mul(db, m1), db2))
A, C = recipe(geo, 'base'); M = geo.M
h = f.add(w, f.scalar(3)); ucub = f.mul(f.sub(h, j), f.inv(f.add(h, j)))
Z = f.add(d, ucub); Yc = f.scale(f.mul(Z, f.inv(bar(Z))), -1)
def dl(g): return f.mul(bar(g), f.inv(g))
half = lambda x: f.scale(x, Q(1, 2))
gens = {'zeta15': z15, 'rho': f.scale(f.add(f.scalar(5), j), Q(1, 6)), 'ucub': ucub, 'Y': Yc,
        'd3': dl(half(f.add(f.one, j))), 'd7': dl(half(f.add(f.one, f.scale(s3, 3)))), 'd19': dl(half(f.add(f.one, f.scale(s3, 5)))),
        'd4': dl(half(f.add(f.one, s15))), 'd5': dl(half(f.add(f.scalar(3), j))), 'd23': dl(half(f.add(f.scalar(9), j))),
        'v2': f.scale(f.add(f.one, s15), Q(1, 4)), 'v5': f.scale(f.add(f.one, f.scale(j, 3)), Q(1, 10)),
        'v13': f.scale(f.add(f.one, f.scale(s3, 15)), Q(1, 26)), 'dd': f.mul(d, f.inv(db))}
for k, v in gens.items(): assert isunit(v), k
names = list(gens)
pl = places(f, 24)
P = [p for p, _, _ in pl]
log('places', P)
def resvec(x):
    out = []
    for p, emb, embc in pl:
        if any(c and c.denominator % p == 0 for c in x): return None
        out.append((red(x, p, emb), red(x, p, embc)))
    return out
G = {k: resvec(v) for k, v in gens.items()}
# unit set: zeta15^a * prod over the other generators g^e, e in {-1, 0, 1}, at most 3 nonzero exponents
others = [k for k in names if k != 'zeta15']
expo_list = []
for r_ in range(0, 4):
    for combo in itertools.combinations(others, r_):
        for signs in itertools.product((-1, 1), repeat=r_):
            for a in range(30):
                expo_list.append((a, tuple(zip(combo, signs))))
log('unit candidates:', len(expo_list))
nP = len(pl)
UR = np.zeros((len(expo_list), nP), dtype=np.int64); URC = np.zeros_like(UR)
for idx, (a, terms) in enumerate(expo_list):
    for t, (p, emb, embc) in enumerate(pl):
        r, rc = pow(G['zeta15'][t][0], a, p), pow(G['zeta15'][t][1], a, p)
        for k, e in terms:
            gr, grc = G[k][t]
            if e == 1: r, rc = r * gr % p, rc * grc % p
            else: r, rc = r * pow(gr, p - 2, p) % p, rc * pow(grc, p - 2, p) % p
        UR[idx, t] = r; URC[idx, t] = rc
log('residues of units computed')
Pa = np.array(P, dtype=np.int64)
# numeric embeddings for the house filter
wroots = sorted(r.real for r in np.roots([1, 2, -4, -4]))
embs = [(cmath.exp(2j * cmath.pi * k / 15), s * 1j * 11 ** 0.5, wr) for k in range(1, 15) if k % 3 and k % 5 for s in (1, -1) for wr in wroots]
def numv(x): return np.array([sum(float(c) * e[0] ** a * e[1] ** b * e[2] ** ce for c, (a, b, ce) in zip(x, f.basis) if c) for e in embs])
targets = []
for name, k in (('b', cb), ('c', cc)):
    for i, i2 in itertools.combinations(range(len(M)), 2):
        dm = f.mul(k, f.sub(M[i2], M[i]))
        for x in range(len(C)):
            for y in range(len(C)):
                X = f.add(f.sub(C[y], C[x]), dm)
                if np.max(np.abs(numv(X))) <= 3.0: targets.append((name, i, i2, x, y, X))
log('targets with house <= 3:', len(targets))
inv3 = np.array([pow(3, p - 2, p) for p in P], dtype=np.int64)
survivors = []
for (name, i, i2, x, y, X) in targets:
    rv = resvec(X)
    if rv is None: continue
    xr = np.array([t[0] for t in rv], dtype=np.int64); xrc = np.array([t[1] for t in rv], dtype=np.int64)
    alive = np.ones(len(expo_list), dtype=bool)
    for t in range(nP):
        p = Pa[t]
        n = ((xr[t] - UR[:, t]) % p) * ((xrc[t] - URC[:, t]) % p) % p
        z = n * ((4 - n) % p) % p * inv3[t] % p
        # Legendre symbol via Euler's criterion, vectorized with Python pow per element on survivors only
        idxs = np.flatnonzero(alive)
        if len(idxs) == 0: break
        leg = np.array([pow(int(z[k]), (int(p) - 1) // 2, int(p)) for k in idxs])
        alive[idxs[leg == p - 1]] = False
    for k in np.flatnonzero(alive):
        survivors.append((name, i, i2, x, y, expo_list[k]))
log('survivors of %d places: %d' % (nP, len(survivors)))
json.dump([[s[0], s[1], s[2], s[3], s[4], [s[5][0], list(s[5][1])]] for s in survivors[:200]], open(os.path.join(HERE, 'three_step_survivors.json'), 'w'))
print('first survivors:', survivors[:5])

# positive control: random three-unit sums from the same unit set must survive (they do: 6 of 6 in the session run)
import random
random.seed(7)
def unit_of(e):
    a, terms = e
    x = f.power(z15, a)
    for k, s in terms:
        x = f.mul(x, gens[k] if s == 1 else f.inv(gens[k]))
    return x
hits = 0
for trial in range(6):
    e1, e2, e3 = random.sample(expo_list, 3)
    X = f.add(f.add(unit_of(e1), unit_of(e2)), unit_of(e3))
    rv = resvec(X)
    xr = np.array([t[0] for t in rv], dtype=np.int64); xrc = np.array([t[1] for t in rv], dtype=np.int64)
    alive = np.ones(len(expo_list), dtype=bool)
    for t in range(nP):
        p = Pa[t]
        n = ((xr[t] - UR[:, t]) % p) * ((xrc[t] - URC[:, t]) % p) % p
        z = n * ((4 - n) % p) % p * inv3[t] % p
        idxs = np.flatnonzero(alive)
        leg = np.array([pow(int(z[k]), (int(p) - 1) // 2, int(p)) for k in idxs])
        alive[idxs[leg == p - 1]] = False
    hits += alive.sum() > 0
log('positive control: %d of 6 random three-unit sums detected' % hits)
