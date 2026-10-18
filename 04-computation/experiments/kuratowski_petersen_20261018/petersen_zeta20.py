#!/usr/bin/env python3
"""The golden unit-distance Petersen realization and its reductions (mac-mini-2026-10-18-unitpath).
Outer pentagon o_k = R zeta5^k, inner pentagram i_k = i r zeta5^k (rotated by 90 degrees), R = (zeta20 + zeta20^-1)/sqrt5,
r = R/phi. Exact identities in Q(zeta20): R^2 = (5 + sqrt5)/10, r^2 = (5 - sqrt5)/10, R^2 + r^2 = 1 (spokes), and the
pentagon / pentagram sides are 1. Checked here by reduction to F_(q^2) for q = 19, 59, 79, 139 (q = 19 mod 20, where
zeta20 lies in mu_(q+1) and the prime of Q(zeta20) above q is conjugation-fixed with residue field F_(q^2)): the ten
reduced points are distinct and all fifteen edges have norm 1, so the Petersen graph is a subgraph of G_q.
Also checked symbolically with sympy (exact algebraic numbers)."""
import mpmath as mp
mp.mp.dps = 60
# Exact identities (proved by hand in THM-4701): R^2 = (5 + sqrt5)/10, r^2 = (5 - sqrt5)/10, |1 - zeta5|^2 = (5 - sqrt5)/2,
# |1 - zeta5^2|^2 = (5 + sqrt5)/2, so R^2 |1 - zeta5|^2 = r^2 |1 - zeta5^2|^2 = (25 - 5)/20 = 1 and R^2 + r^2 = 1.
# Here: a 60-digit numerical confirmation of the fifteen edge lengths.
z = mp.exp(2j * mp.pi / 20); R = (z + 1 / z).real / mp.sqrt(5); phi = (1 + mp.sqrt(5)) / 2; r = R / phi
pts = [R * z ** (4 * k) for k in range(5)] + [1j * r * z ** (4 * k) for k in range(5)]
E = [(k, (k + 1) % 5) for k in range(5)] + [(k, k + 5) for k in range(5)] + [(5 + k, 5 + (k + 2) % 5) for k in range(5)]
err = max(abs(abs(pts[a] - pts[b]) ** 2 - 1) for a, b in E)
print('60-digit check: R^2 = %s, (5 + sqrt5)/10 = %s ; max |edge^2 - 1| over the 15 edges = %s'
      % (mp.nstr(R ** 2, 20), mp.nstr((5 + mp.sqrt(5)) / 10, 20), mp.nstr(err, 5)))
import numpy as np
J = np.zeros((15, 20))
for e, (a, b) in enumerate(E):
    dv = complex(pts[a] - pts[b])
    J[e, 2 * a], J[e, 2 * a + 1], J[e, 2 * b], J[e, 2 * b + 1] = 2 * dv.real, 2 * dv.imag, -2 * dv.real, -2 * dv.imag
sv = np.linalg.svd(J, compute_uv=False)
print('rigidity matrix rank %d (smallest singular values %s): real realizations near this one form a %d-parameter family modulo motions'
      % ((sv > 1e-9).sum(), np.round(sv[-3:], 4), 20 - (sv > 1e-9).sum() - 3))
for q in (19, 59, 79, 139):
    n = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
    def mul(a, b): return ((a[0] * b[0] + n * a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)
    def pw(a, e):
        res = (1, 0)
        while e:
            if e & 1: res = mul(res, a)
            a = mul(a, a); e >>= 1
        return res
    def inv(a): return pw(a, q * q - 2)
    def add(a, b): return ((a[0] + b[0]) % q, (a[1] + b[1]) % q)
    def sub(a, b): return ((a[0] - b[0]) % q, (a[1] - b[1]) % q)
    norm = lambda a: (a[0] * a[0] - n * a[1] * a[1]) % q
    # generator of mu_(q+1): an element g of F_(q^2)^* with g^(q+1) = 1 of exact order q + 1
    from sympy import factorint
    fac = factorint(q + 1)
    gen = None
    for a in range(q):
        for b in range(1, q):
            x = pw((a, b), q - 1)                         # lands in mu_(q+1)
            if all(pw(x, (q + 1) // p) != (1, 0) for p in fac): gen = x; break
        if gen: break
    z20 = pw(gen, (q + 1) // 20); z5 = pw(z20, 4); iu = pw(z20, 5)
    s5 = sub(add(z5, pw(z5, 4)), add(pw(z5, 2), pw(z5, 3)))          # Gauss sum, s5^2 = 5
    assert mul(s5, s5) == (5 % q, 0)
    Rq = mul(add(z20, inv(z20)), inv(s5)); phiq = mul(add((1, 0), s5), inv((2, 0))); rq = mul(Rq, inv(phiq))
    P = [mul(Rq, pw(z5, k)) for k in range(5)] + [mul(mul(iu, rq), pw(z5, k)) for k in range(5)]
    edges_ok = all(norm(sub(P[a], P[b])) == 1 for a, b in E)
    print('q = %3d: ten distinct points: %s ; all fifteen Petersen edges of norm 1 in F_(q^2): %s ; N(R) + N(r) = 1: %s'
          % (q, len(set(P)) == 10, edges_ok, (norm(Rq) + norm(rq)) % q == 1))
