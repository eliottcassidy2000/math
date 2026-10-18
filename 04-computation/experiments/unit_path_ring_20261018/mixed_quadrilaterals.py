#!/usr/bin/env python3
"""Explicit mixed quadrilaterals (two unit sides, two sides d_k = 1/|k|) in the six-colour field (mac-mini-2026-10-18-unitpath).

Mixed triangles {1, d_b} or {1, d_c} do not exist (THM-4700(4)). A closed walk u1 + u2 + (e1 + e2)/k = 0 with units
u_i, e_i and z = u1 + u2 != 0 exists iff L(z) <= 2 and L(kz) <= 2. With mu = N(kz) = 4/(1 + 3t^2) (a two-sum norm)
and nu = N(z) = mu/N, N = k conj(k), the remaining condition is the conic
    sigma^2 - N t^2 = (N - 1)/3      over F+.
Its quaternion algebra is (N, (N-1)/3) = (N, -3), split because N = N_(F/F+)(k). Explicitly, writing k = A + B sqrt-3
(A, B in F+, B != 0): t = (1 + A)/(3B), sigma = (A + N)/(3B), since sigma + t sqrt N = (1 + sqrt N)(A + sqrt N)/(3B) and
N(1 + sqrt N) = 1 - N, N((A + sqrt N)/B) = -3.
Then z = 2/(conj(k)(1 + t sqrt-3)), and the two-sum splittings x = (x/2)(1 + sqrt-3 s/N(x)) + (x/2)(1 - sqrt-3 s/N(x)),
s^2 = N(x)(4 - N(x))/3, give the four steps. Everything is checked in exact arithmetic."""
import os, sys
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__)); EXP = os.path.join(HERE, '..')
sys.path.insert(0, EXP); sys.path.insert(0, os.path.join(EXP, 'first_arrival_20261018', 'shim'))
from golden_forcing_search_20261017 import Geo
geo = Geo(); f = geo.f; bar = geo.bar; j = geo.j
z15 = f.gen(0); s3 = f.add(f.scale(f.power(z15, 5), 2), f.one); assert f.mul(s3, s3) == f.scalar(-3)
d = geo.d; db = bar(d)
m2 = f.scale(f.add(f.scalar(20), f.scale(j, 11)), Q(-1, 19)); m1 = f.scale(f.sub(f.scalar(20), f.scale(j, 11)), Q(1, 19))
db2 = f.mul(db, db)
cc = f.mul(f.add(f.sub(f.add(f.mul(db, m2), f.one), f.mul(db2, m1)), f.mul(db2, db)), f.inv(f.sub(f.mul(d, db), f.one)))
cb = f.add(f.mul(cc, d), f.sub(f.mul(db, m1), db2))
isunit = lambda x: f.mul(x, bar(x)) == f.one
isreal = lambda x: bar(x) == x
inv_s3 = f.inv(s3)
ok = True
def two_sum(x, s):
    """x = v1 + v2 with v1, v2 units, given s in F+ with s^2 = N(x)(4 - N(x))/3."""
    Nx = f.mul(x, bar(x)); r = f.mul(f.mul(s3, s), f.inv(Nx))
    half = f.scale(x, Q(1, 2))
    return f.mul(half, f.add(f.one, r)), f.mul(half, f.sub(f.one, r))
for name, k in (('a', db), ('b', cb), ('c', cc)):
    N = f.mul(k, bar(k))
    A = f.scale(f.add(k, bar(k)), Q(1, 2))                        # k = A + B sqrt-3
    B = f.mul(f.scale(f.sub(k, bar(k)), Q(1, 2)), inv_s3)
    assert isreal(A) and isreal(B) and f.add(A, f.mul(B, s3)) == k
    iB = f.inv(B)
    t = f.scale(f.mul(f.add(f.one, A), iB), Q(1, 3)); sig = f.scale(f.mul(f.add(A, N), iB), Q(1, 3))
    conic = f.sub(f.mul(sig, sig), f.mul(N, f.mul(t, t))) == f.scale(f.sub(N, f.one), Q(1, 3))
    one3t2 = f.add(f.one, f.scale(f.mul(t, t), 3))
    z = f.scale(f.inv(f.mul(bar(k), f.add(f.one, f.mul(t, s3)))), 2)
    kz = f.mul(k, z)
    mu = f.mul(kz, bar(kz)); nu = f.mul(z, bar(z))
    s_mu = f.scale(f.mul(t, f.inv(one3t2)), 4)                     # s^2 = mu (4 - mu)/3
    s_nu = f.mul(f.scale(sig, 4), f.inv(f.mul(N, one3t2)))          # s^2 = nu (4 - nu)/3
    chk_mu = f.mul(s_mu, s_mu) == f.scale(f.mul(mu, f.sub(f.scalar(4), mu)), Q(1, 3))
    chk_nu = f.mul(s_nu, s_nu) == f.scale(f.mul(nu, f.sub(f.scalar(4), nu)), Q(1, 3))
    u1, u2 = two_sum(z, s_nu)
    e1, e2 = two_sum(f.scale(kz, -1), s_mu)
    ik = f.inv(k)
    closed = f.add(f.add(u1, u2), f.mul(f.add(e1, e2), ik)) == f.zero
    units = all(isunit(x) for x in (u1, u2, e1, e2))
    nonpar = f.add(u1, u2) != f.zero
    # the four vertices 0, u1, z, z + e1/k are distinct
    verts = [f.zero, u1, z, f.add(z, f.mul(e1, ik))]
    distinct = len(set(verts)) == 4
    good = conic and chk_mu and chk_nu and closed and units and nonpar and distinct
    ok &= good
    print('k = %s: conic point %s ; two-sum norms %s/%s ; closed walk u1 + u2 + (e1 + e2)/k = 0: %s ; steps are units: %s ; non-parallelogram: %s ; 4 distinct vertices: %s'
          % (name, conic, chk_nu, chk_mu, closed, units, nonpar, distinct))
print('ALL MIXED QUADRILATERALS VERIFIED' if ok else 'FAILURE')
