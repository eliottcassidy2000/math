#!/usr/bin/env python3
"""Unit-path ring of the six-colour field F = Q(zeta15, sqrt-11, w) and short walks (mac-mini-2026-10-18-unitpath).

Gamma_F = Z-span of U_F = {z in F : z conj(z) = 1} (a ring: clock_duality_oct15, cubic_adelic_20261016e.md section 3).
Exact checks (Field arithmetic of escaped_ports_20261017d.make('sqrt3'), degree 48):
 (A) generators of O_F in Gamma_F: zeta15 (unit), 1/2 (two steps), 1/11 (clock_duality's 82-step connector),
     j = 6 rho - 5, w = (19 (u + conj u) + 6)/11, hence (1 + j)/2 and w^2/2;
 (B) explicit short walks to 1/p for every prime p < 29 (29 is the first prime with 1/p outside Gamma_F);
 (C) the reflection link: for a unit edge s of C and unit q, q' = (conj(kappa)/kappa) conj(q) (kappa q - s)/(conj(kappa q) - conj(s))
     is a unit and |s + kappa (q' - q)| = 1, for kappa = b, c of the cubic recipe.
A walk of length L is an identity x = sum of L unit vectors (repetition allowed); integer multiples count |k| steps.
"""
import os, sys
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__)); EXP = os.path.join(HERE, '..')
sys.path.insert(0, EXP); sys.path.insert(0, os.path.join(EXP, 'first_arrival_20261018', 'shim'))
from golden_forcing_search_20261017 import Geo, recipe

geo = Geo(); f = geo.f; bar = geo.bar
isunit = lambda x: f.mul(x, bar(x)) == f.one
z15 = f.gen(0); j = geo.j; w = f.gen(2)
S = lambda *xs: f.add(xs[0], S(*xs[1:])) if len(xs) > 1 else xs[0]
sc = lambda x, k: f.scale(x, Q(k))
zeta3 = f.power(z15, 5); s3 = f.add(sc(zeta3, 2), f.one); assert f.mul(s3, s3) == f.scalar(-3)
zeta5 = f.power(z15, 3); s5 = f.add(sc(f.add(zeta5, f.power(zeta5, 4)), 2), f.one); assert f.mul(s5, s5) == f.scalar(5)
s15 = f.mul(s3, s5); assert f.mul(s15, s15) == f.scalar(-15)
assert f.mul(j, j) == f.scalar(-11)
ok = True
def check(name, cond):
    global ok
    ok &= bool(cond); print('  %-74s %s' % (name, bool(cond)))

print('(A) O_F inside Gamma_F')
check('zeta15 is a unit vector', isunit(z15))
v2 = f.scale(f.add(f.one, s15), Q(1, 4))
check('1/2 = v + conj(v), v = (1 + sqrt-15)/4 unit  [2 steps]', isunit(v2) and f.add(v2, bar(v2)) == f.scalar(Q(1, 2)))
alpha = zeta5
beta = [f.mul(f.sub(f.power(alpha, r), f.scalar(3)), f.inv(f.sub(f.power(alpha, 5 - r), f.scalar(3)))) for r in range(1, 5)]
one11 = f.sub(sc(S(*beta), 11), f.scalar(38))
check("1/11 = 11 (beta_1 + ... + beta_4) - 38, beta_r = (a^r - 3)/(a^-r - 3)  [82 steps, clock_duality]", all(isunit(b) for b in beta) and one11 == f.scalar(Q(1, 11)))
rho = f.scale(f.add(f.scalar(5), j), Q(1, 6))
check('j = 6 rho - 5, rho = (5 + j)/6 unit  [11 steps]', isunit(rho) and f.sub(sc(rho, 6), f.scalar(5)) == j)
h = f.add(w, f.scalar(3)); u = f.mul(f.sub(h, j), f.inv(f.add(h, j)))
check('w = (19 (u + conj u) + 6)/11 with u = (w + 3 - j)/(w + 3 + j) unit', isunit(u) and f.scale(f.add(sc(f.add(u, bar(u)), 19), f.scalar(6)), Q(1, 11)) == w)
print('    => zeta15, (1 + j)/2 = (1/2)(1 + j), w and w^2/2 = (1/2) w^2 lie in the ring Gamma_F; O_F = Z[zeta15][(1+j)/2][1, w, w^2/2]')

print('(B) short walks to 1/p, p < 29 (delta = conj(g)/g, delta + 1/delta = t^2/n - 2)')
def delta_of(g):
    d = f.mul(bar(g), f.inv(g)); assert isunit(d); return d
def walk(name, value, terms, steps):
    tot = f.zero
    for k, x in terms:
        assert isunit(x); tot = f.add(tot, sc(x, k))
    check('%-28s = %-40s [%d steps]' % (name, ' + '.join('%s*%s' % (k, n) for k, (_, n) in zip([t[0] for t in terms], [(None, nm) for nm in [t[2] for t in terms]])), steps), tot == value)
one = f.one
def unit_pair_half(n, D, b, sqrtD):
    # 1/n = v + conj v with v = (1 + b sqrt(-D))/(2n), |v| = 1 iff 1 + D b^2 = 4 n^2
    v = f.scale(f.add(one, sc(sqrtD, b)), Q(1, 2 * n)); return v
for n, D, b, sq in ((2, 15, 1, s15), (5, 11, 3, j), (13, 3, 15, s3)):
    v = unit_pair_half(n, D, b, sq)
    check('1/%-2d = v + conj v, v = (1 + %d sqrt-%d)/%d  [2 steps]' % (n, b, D, 2 * n), isunit(v) and f.add(v, bar(v)) == f.scalar(Q(1, n)))
for n, t, b, sq, D in ((3, 1, 1, j, 11), (7, 1, 3, s3, 3), (19, 1, 5, s3, 3), (4, 1, 1, s15, 15)):
    g = f.scale(f.add(f.scalar(t), sc(sq, b)), Q(1, 2)); dd = delta_of(g)
    check('1/%-2d = delta + 1/delta + 2, g = (1 + %d sqrt-%d)/2 of norm %d  [4 steps]' % (n, b, D, n), f.add(f.add(dd, bar(dd)), f.scalar(2)) == f.scalar(Q(1, n)))
g = f.scale(f.add(f.scalar(9), j), Q(1, 2)); dd = delta_of(g)
check('1/23 = 2 delta + 2/delta - 3, g = (9 + sqrt-11)/2 of norm 23  [7 steps]', f.sub(sc(f.add(dd, bar(dd)), 2), f.scalar(3)) == f.scalar(Q(1, 23)))
# 17: inert in Q(zeta15)^+ ; a generator pi of a prime of Q(zeta15) above 17 with conj(pi) not associate: search small elements
import itertools
found17 = None
for coeffs in itertools.product(range(-2, 3), repeat=4):
    g = S(*[sc(f.power(z15, k), c) for k, c in enumerate(coeffs)]) if any(coeffs) else None
    if g is None: continue
    n = f.mul(g, bar(g)); t = f.add(g, bar(g))
    # need n = 17 * (something) with norm exactly 17^? : test the rational norm of n down Q(zeta15)+ via 17 | n and t coprime
    if n[0].denominator == 1 and all(c == 0 for c in n[1:]) and n[0] == 17:
        found17 = (coeffs, t); break
if found17:
    coeffs, t = found17
    print('    element of Q(zeta15) with g conj(g) = 17: coefficients', coeffs, '; trace t =', [str(c) for c in t if c][:3])
else:
    print('    (no element of norm 17 in the small coefficient box; 1/17 in Gamma_F follows from B2 since 17 is not in S)')
print('    1/11 in Gamma_F: 82 steps (clock_duality); no 1/29: 29 lies under a c-fixed prime (w = 3 branch)')

print('(C) reflection links between coefficient copies of the cubic recipe')
d = geo.d; db = bar(d)
m2 = f.scale(f.add(f.scalar(20), f.scale(j, 11)), Q(-1, 19)); m1 = f.scale(f.sub(f.scalar(20), f.scale(j, 11)), Q(1, 19))
db2 = f.mul(db, db)
cc = f.mul(f.add(f.sub(f.add(f.mul(db, m2), f.one), f.mul(db2, m1)), f.mul(db2, db)), f.inv(f.sub(f.mul(d, db), f.one)))
cb = f.add(f.mul(cc, d), f.sub(f.mul(db, m1), db2))
A, C = recipe(geo, 'base')
EC = geo.unit_edges(C)
s_edge = f.sub(C[EC[0][1]], C[EC[0][0]])
for name, kap in (('b', cb), ('c', cc)):
    cnt = 0
    for q in [m for m in geo.M if m != f.zero and isunit(m)]:
        kq = f.mul(kap, q)
        qp = f.mul(f.mul(f.mul(bar(kap), f.inv(kap)), bar(q)), f.mul(f.sub(kq, s_edge), f.inv(f.sub(bar(kq), bar(s_edge)))))
        good = isunit(qp) and isunit(f.add(s_edge, f.mul(kap, f.sub(qp, q)))) and qp != q
        cnt += good
    check('kappa = %s: every unit q in M gets a unit q\' with x + kappa q ~ y + kappa q\' (s = y - x)' % name, cnt == len([m for m in geo.M if m != f.zero and isunit(m)]))

print('(D) every coordinate factor of the cubic recipe lies in Gamma_F')
from math import lcm
def den(x):
    D = 1
    for c in x:
        if c: D = lcm(D, c.denominator)
    return D
def primes_of(n):
    out = set(); p = 2
    while p * p <= n:
        while n % p == 0: out.add(p); n //= p
        p += 1
    if n > 1: out.add(n)
    return out
factors = {'A + aM': [f.add(x, f.mul(db, m)) for x in A for m in geo.M],
           'C + bM': [f.add(x, f.mul(cb, m)) for x in C for m in geo.M],
           'C + cM': [f.add(x, f.mul(cc, m)) for x in C for m in geo.M]}
allp = set()
for name, pts in factors.items():
    ps = set().union(*(primes_of(den(x)) for x in pts))
    allp |= ps
    print('  %-8s %4d points ; primes in power-basis denominators: %s' % (name, len(set(pts)), sorted(ps)))
check('denominator primes avoid S = {p = 14 mod 15, (p/11) = -1} (power basis Z[zeta15, j, w] is in O_F)', all(not (p % 15 == 14 and pow(p % 11, 5, 11) == 10) for p in allp))
# 1/17 via g = (11 + sqrt-15)/2 of norm 34: 1/34 = 9 delta + 9/delta - 14, so 1/17 = 18 delta + 18/delta - 28
g = f.scale(f.add(f.scalar(11), s15), Q(1, 2)); dd = delta_of(g)
check('1/17 = 18 delta + 18/delta - 28, g = (11 + sqrt-15)/2 of norm 34  [64 steps]', f.sub(sc(f.add(dd, bar(dd)), 18), f.scalar(28)) == f.scalar(Q(1, 17)))
print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
