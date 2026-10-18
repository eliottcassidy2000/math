#!/usr/bin/env python3
"""Heron census of coefficient triangles in the six-colour field (mac-mini-2026-10-18-unitpath).

F = Q(zeta15, sqrt-11, w) is a CM field, F = F+(sqrt-3), conjugation c.
Lemma (Heron criterion, elementary): for totally positive R1, R2, R3 in F+, there are z1, z2 in F with
  z1 conj(z1) = R1, z2 conj(z2) = R2, (z1 - z2) conj(z1 - z2) = R3, z1/z2 not real
(a nondegenerate triangle with squared sides R1, R2, R3 and a vertex at 0) iff some element of norm R1 exists and
  H = 2(R1 R2 + R2 R3 + R3 R1) - (R1^2 + R2^2 + R3^2)   (= 16 Area^2)
lies in 3 (F+)^{*2}. [y = z2/z1 has y + conj y = (R1 + R2 - R3)/R1 = T and y conj y = R2/R1; y in F iff
T^2 - 4 R2/R1 = -H/R1^2 is a square in F; z^2 in F+ forces z in F+ or z in sqrt-3 F+, and y non-real needs the second.]
Coefficient pairs (m, m') of the cubic recipe interact at |m - m'| = 1 (u-square family), d_a = 1/|a|, d_b = 1/|b|, d_c = 1/|c| (axis-0,
axis-1, axis-4 translation matchings), where u^2 = a + b u + c Y, a = conj(1 + zeta5).
We test H/3 for squareness in F by quadratic characters (one non-residue proves non-square) and, when it survives,
check total positivity of H numerically in all 24 real embeddings of F+ (then H/3 in F^2 cap F+ positive => in (F+)^2).
"""
import os, sys, json, time, itertools, cmath
from fractions import Fraction as Q
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); EXP = os.path.join(HERE, '..')
sys.path.insert(0, EXP); sys.path.insert(0, os.path.join(EXP, 'first_arrival_20261018', 'shim'))
from golden_forcing_search_20261017 import Geo
from golden_vsquare_unit_census_20261017 import places, red

geo = Geo(); f = geo.f; bar = geo.bar
d = geo.d; db = bar(d); w = f.gen(2); j = geo.j
m2 = f.scale(f.add(f.scalar(20), f.scale(j, 11)), Q(-1, 19)); m1 = f.scale(f.sub(f.scalar(20), f.scale(j, 11)), Q(1, 19))
db2 = f.mul(db, db)
cc = f.mul(f.add(f.sub(f.add(f.mul(db, m2), f.one), f.mul(db2, m1)), f.mul(db2, db)), f.inv(f.sub(f.mul(d, db), f.one)))
cb = f.add(f.mul(cc, d), f.sub(f.mul(db, m1), db2))
ca = db
pl = places(f, 60)

# all complex embeddings: zeta15 -> e(k/15), k in (Z/15)^*, j -> +- i sqrt 11, w -> each real root
wroots = sorted(r.real for r in np.roots([1, 2, -4, -4]))
embs = [(cmath.exp(2j * cmath.pi * k / 15), s * 1j * 11 ** 0.5, wr) for k in range(1, 15) if k % 3 and k % 5 for s in (1, -1) for wr in wroots]
def ev(x, e):
    zc, jc, wc = e
    return sum(float(c) * zc ** a * jc ** b * wc ** cexp for c, (a, b, cexp) in zip(x, f.basis) if c)

def screen(x):
    cnt = 0
    for p, emb, embc in pl:
        if any(c and c.denominator % p == 0 for c in x): continue
        r = red(x, p, emb)
        if r == 0: continue
        cnt += 1
        if pow(r, (p - 1) // 2, p) == p - 1: return 'NON-square (p = %d)' % p
    return 'square at all %d places' % cnt

R = {'1': f.one}
for name, k in (('a', ca), ('b', cb), ('c', cc)):
    N = f.mul(k, bar(k)); R['d_%s' % name] = f.inv(N)      # side d_k = 1/|k|, squared length 1/N_k
geom = embs[0]
print('field degree', f.n, '; %d complex embeddings' % len(embs))
for name, x in R.items():
    vals = [ev(x, e).real for e in embs]
    print('  side %-4s squared length: min/max over the 48 embeddings %.5f / %.5f ; value at the first listed embedding %.6f' % (name, min(vals), max(vals), vals[0]))
rows = []
for combo in itertools.combinations_with_replacement(list(R), 3):
    R1, R2, R3 = (R[n] for n in combo)
    S = f.add(f.add(f.mul(R1, R2), f.mul(R2, R3)), f.mul(R3, R1))
    H = f.sub(f.scale(S, 2), f.add(f.add(f.mul(R1, R1), f.mul(R2, R2)), f.mul(R3, R3)))
    Hv = [ev(H, e).real for e in embs]
    tp = min(Hv) > 0
    s = screen(f.scale(H, Q(1, 3)))
    verdict = 'TRIANGLE EXISTS in F' if (s.startswith('square') and tp) else 'no'
    rows.append(dict(sides=combo, H_min=min(Hv), H_max=max(Hv), totally_positive=tp, H_over_3=s, verdict=verdict))
    print('%-28s H range [%9.4f, %9.4f] tot.pos %-5s H/3: %-28s -> %s' % (' , '.join(combo), min(Hv), max(Hv), tp, s, verdict))
json.dump(rows, open(os.path.join(HERE, 'heron_census.json'), 'w'), indent=1)
