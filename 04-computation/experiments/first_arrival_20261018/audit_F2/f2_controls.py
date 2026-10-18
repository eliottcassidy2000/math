#!/usr/bin/env python3
"""Audit F2: controls for the residue certificates of f2_connectors.py.
  1. homomorphism control: red(x y) = red(x) red(y) and red(conj x) = red_conj(x) at every embedding of 3 primes;
  2. negative control: X = s^2 for s in F is never certified as a non-square (all residues are squares or 0);
  3. the survivors' X = N(t)(N(t)-4) is nonzero and not a square at a prime of a different residue characteristic class
     (an extra prime p = 1 mod 30 above 10^4), so the certificates do not depend on the first 40 primes.
"""
import sys, os, json, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f2_field import E, ONE, J, W, N, build, roots_of_unity, conj, norm_cc
from f2_connectors import places_all, red, red_conj

random.seed(20261018)
O = build()
b, c, u, M, C = O['b'], O['c'], O['u'], O['M'], O['C']
PL = places_all(3, start=1000)
def rnd():
    return E([random.randint(-5, 5) for _ in range(N)], random.choice([1, 2, 3, 4, 6, 19]))
ok_hom = True; ok_conj = True
for _ in range(20):
    x, y = rnd(), rnd()
    if x.d % 19 == 0 or y.d % 19 == 0:
        pass
    for p, embs in PL:
        for emb in embs:
            if red(x * y, p, emb) != red(x, p, emb) * red(y, p, emb) % p: ok_hom = False
            if red(conj(x), p, emb) != red_conj(x, p, emb): ok_conj = False
print('homomorphism control (20 random pairs, all embeddings of %s): %s ; conjugate-place control: %s'
      % ([p for p, _ in PL], ok_hom, ok_conj))
bad = 0; tests = 0
for _ in range(30):
    s = rnd() + b * rnd()
    X = s * s
    for p, embs in PL:
        for emb in embs:
            r = red(X, p, emb); tests += 1
            if r and pow(r, (p - 1) // 2, p) == p - 1: bad += 1
print('negative control: squares X = s^2 certified as non-squares in %d of %d residue tests (must be 0)' % (bad, tests))
# extra prime above 10^4
PX = places_all(1, start=10000)
p, embs = PX[0]
HERE = os.path.dirname(os.path.abspath(__file__)); SESSION = os.path.dirname(HERE)
EC = [(i, k) for i in range(len(C)) for k in range(i + 1, len(C)) if norm_cc(C[k] - C[i]) == ONE]
U0 = sorted({C[k] - C[i] for i, k in EC} | {C[i] - C[k] for i, k in EC}, key=lambda e: e.fracs())
roots = roots_of_unity()
D90 = sorted({r * q for r in roots for q in U0 + [ONE]}, key=lambda e: e.fracs())
ui = u.inv()
D120 = sorted({r * q for r in roots for q in (u, u * u, ui, ui * ui)}, key=lambda e: e.fracs())
res = []
for fname, D in (('copy_connectors.json', U0), ('copy_connectors_roots.json', D90), ('copy_connectors_udirs.json', D120)):
    Jn = json.load(open(os.path.join(SESSION, fname)))
    for name, co in (('b', b), ('c', c)):
        for (i1, i2, a1, a2, qi, nt) in Jn[name]['three_step_survivors']:
            t = (C[a2] + co * M[i2]) - (C[a1] + co * M[i1]) - D[qi]
            Nt = norm_cc(t); X = Nt * (Nt - 4)
            nres = sum(1 for emb in embs if red(X, p, emb) and pow(red(X, p, emb), (p - 1) // 2, p) == p - 1)
            res.append(nres > 0)
print('extra prime p = %d (%d embeddings): every recorded survivor certified there too: %s (%d survivors)' % (p, len(embs), all(res), len(res)))
