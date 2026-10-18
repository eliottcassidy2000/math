#!/usr/bin/env python3
"""Valuations at the 24 primes above 19 of F = Q(zeta15)(sqrt(-11))(w) (mac-mini-2026-10-18-sixcolour).
The tower order Z[zeta15][j][w] is 19-maximal and 19 is unramified, so the primes above 19 are the Frobenius orbits of the
48 ring homomorphisms Z[zeta15, j, w] -> F_361 (zeta -> a primitive 15th root of unity, j -> a square root of -11, w -> a
root of w^3 + 2w^2 - 4w - 4), all orbits of size 2, residue fields F_361. For x in F with 19^k exactly dividing the
coordinate denominators, y = 19^k x is 19-integral and v_p(x) = -k if y does not vanish at p, else v_p(x) >= 1 - k.
Elements: the circuit coefficients b, c; the contact unit u = (h - j)/(h + j), h = w + 3, and the mixed-golden unit
u2 = (h2 - j)/(h2 + j), h2 = w + 3 + phi (escaped_ports make()); a = conj(d); the spindle differences."""
import os, sys, itertools
from fractions import Fraction as Q
from math import lcm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from golden_forcing_search_20261017 import Geo, recipe
from escaped_ports_20261017d import make
P = 19; NR = (-11) % P            # F_361 = F_19[t]/(t^2 - 8), t = sqrt(-11)
def fm(x, y): return ((x[0] * y[0] + NR * x[1] * y[1]) % P, (x[0] * y[1] + x[1] * y[0]) % P)
def fpow(x, k):
    r = (1, 0)
    for _ in range(k): r = fm(r, x)
    return r
def fadd(x, y): return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)
els = [(a, b) for a in range(P) for b in range(P)]
z15 = [x for x in els if fpow(x, 15) == (1, 0) and fpow(x, 5) != (1, 0) and fpow(x, 3) != (1, 0)]
jr = [(0, 1), (0, P - 1)]
def cub(x): return fadd(fadd(fpow(x, 3), fm((2, 0), fpow(x, 2))), fadd(fm((P - 4, 0), x), (P - 4, 0)))
wr = [x for x in els if cub(x) == (0, 0)]
geo = Geo(); f = geo.f; A, C = recipe(geo, 'base')
for poly, rs in zip(f.polys, (z15, jr, wr)):
    for r in rs:
        v = (0, 0)
        for c in reversed(poly): v = fadd(fm(v, r), (c % P, 0))
        assert v == (0, 0)
homs = [(z, j, w) for z in z15 for j in jr for w in wr]
frob = lambda x: fpow(x, P)
primes, seen = [], set()
for hm in homs:
    if hm in seen: continue
    orb = (hm, tuple(frob(x) for x in hm)); seen.update(orb); primes.append(orb[0])
print('primitive 15th roots in F_361: %d ; sqrt(-11): 2 ; cubic roots: %s ; homomorphisms %d ; primes above 19: %d' % (len(z15), wr, len(homs), len(primes)))
def hval(x, hm):
    z, j, w = hm; s = (0, 0)
    for c, e in zip(x, f.basis):
        if c:
            assert c.denominator % P
            cc = c.numerator * pow(c.denominator, -1, P) % P
            s = fadd(s, fm((cc, 0), fm(fpow(z, e[0]), fm(fpow(j, e[1]), fpow(w, e[2])))))
    return s
def v19den(x):
    D = lcm(*[c.denominator for c in x]); k = 0
    while D % P == 0: D //= P; k += 1
    return k
def val(x, hm):
    k = v19den(x); y = f.scale(x, Q(P ** k))
    return -k if hval(y, hm) != (0, 0) else '>=%d' % (1 - k)
d, db = geo.d, geo.bar(geo.d)
m2c = f.scale(f.add(f.scalar(20), f.scale(geo.j, 11)), Q(-1, 19)); m1c = f.scale(f.sub(f.scalar(20), f.scale(geo.j, 11)), Q(1, 19))
db2 = f.mul(db, db)
cc = f.mul(f.add(f.sub(f.add(f.mul(db, m2c), f.one), f.mul(db2, m1c)), f.mul(db2, db)), f.inv(f.sub(f.mul(d, db), f.one)))
cb = f.add(f.mul(cc, d), f.sub(f.mul(db, m1c), db2))
_, _, _, _, _, units, _, _, _, _ = make('sqrt3')
u, u2 = units
w0 = f.gen(2); h = f.add(w0, f.scalar(3)); assert f.mul(f.sub(h, geo.j), f.inv(f.add(h, geo.j))) == u
print('prime | h residue type | v(b) v(c) | v(u) v(u2) | v(a)')
tab = []
for hm in primes:
    hres = fadd(hm[2], (3, 0)); typ = 'h=+j' if hres == hm[1] else ('h=-j' if hres == ((-hm[1][0]) % P, (-hm[1][1]) % P) else 'h=%s' % (hres,))
    row = (val(cb, hm), val(cc, hm), val(u, hm), val(u2, hm), val(db, hm))
    tab.append((typ, row)); print(' %s | %s | %s %s | %s %s | %s' % (hm, typ, *row))
blocked = [i for i, (typ, (vb, vc, vu, vu2, va)) in enumerate(tab) if vb == -1 and vu == 0 and vu2 == 0]
print('primes with v(b) = -1 where u and u2 are both units:', len(blocked))
# per-prime diagonal resonances: (dm2, dm3) with b dm2 + c dm3 integral at that single prime
Dl = sorted({f.sub(y, x) for x in geo.M for y in geo.M}); zero = f.sub(f.one, f.one)
sm = all(hval(f.scale(x, Q(P ** v19den(x))), hm) != (0, 0) for x in Dl if x != zero for hm in primes)
print('every nonzero spindle difference is a unit at every prime above 19:', sm, '; all spindle differences 19-integral:', all(v19den(x) == 0 for x in Dl))
# --- refinement: exact valuations of u at the three prime types, and diagonal sheet moves at the h = 7 primes ---
def typ_of(hm):
    hres = fadd(hm[2], (3, 0))
    return '+j' if hres == hm[1] else ('-j' if hres == ((-hm[1][0]) % P, (-hm[1][1]) % P) else '7')
num, den_ = f.sub(h, geo.j), f.add(h, geo.j)          # u = num/den_, both 19-integral
for t in ('+j', '-j', '7'):
    ps = [hm for hm in primes if typ_of(hm) == t]
    nz = lambda x: [hval(x, hm) != (0, 0) for hm in ps]
    print('type h = %2s: %d primes; h - j nonzero at all: %s ; h + j nonzero at all: %s' % (t, len(ps), all(nz(num)), all(nz(den_))))
P7 = [hm for hm in primes if typ_of(hm) == '7']
res = []
for x2 in Dl:
    for x3 in Dl:
        if x2 == zero and x3 == zero: continue
        y = f.add(f.mul(cb, x2), f.mul(cc, x3)); k = v19den(y); yy = f.scale(y, Q(P ** k))
        integral_at = sum(1 for hm in P7 if k == 0 or hval(yy, hm) == (0, 0))
        res.append(integral_at)
from collections import Counter
print('nonzero (dm2, dm3): number of h = 7 primes at which b dm2 + c dm3 is integral -> count:', sorted(Counter(res).items()))
print('pairs integral at ALL eight h = 7 primes:', sum(1 for r in res if r == len(P7)))
