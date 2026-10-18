#!/usr/bin/env python3
"""Referee check B2 (second S41 audit): [F(s, h) : F] = 4, i.e. D, R and DR are non-squares in F = Q(zeta15, j, w).

The proof of THM-4668 compares the four coefficients of 1, s, h, sh over F, so it needs {1, s, h, sh} to be an F-basis.
At a degree-one place of F (p = 1 mod 15, (-11/p) = 1, the cubic has a root mod p; p not in {2, 3, 5, 11, 37}, so the
order Z[zeta15, j, w] is p-maximal and the kernel is a prime of O_F with residue field F_p) a square of F reduces to a
square of F_p. A non-residue therefore certifies a non-square. Own evaluation from the defining formulas:
u = (w + 3 - j)/(w + 3 + j), conj u = 1/u, d = 1 + zeta5, conj d = 1 + zeta5^-1, N = Z conj(Z), D = 1 - 4/N, R = 4 - N.
Also replays clock_duality's two places at 31 (golden_square_contexts_20261017h.md, G2 table).
Usage: python b2_degree192.py
"""
import sympy as sp


def place_values(p, z0, j0, w0):
    inv = lambda x: pow(x % p, -1, p)
    u = (w0 + 3 - j0) * inv(w0 + 3 + j0) % p
    ub = (w0 + 3 + j0) * inv(w0 + 3 - j0) % p
    assert u * ub % p == 1
    a0 = pow(z0, 3, p)
    d = (1 + a0) % p; db = (1 + inv(a0)) % p
    Z = (d + u) % p; Zb = (db + ub) % p
    N = Z * Zb % p
    D = (1 - 4 * inv(N)) % p
    R = (4 - N) % p
    return N, D, R


def leg(x, p):
    return sp.legendre_symbol(x % p, p) if x % p else 0


print('clock_duality G2 rows at p = 31 (zeta15, j, w) -> (N, D, R, chi(D), chi(R), chi(DR)):')
for z0 in (7, 10):
    N, D, R = place_values(31, z0, 19, 25)
    print('  (%d, 19, 25): N = %d, D = %d, R = %d, chi = (%d, %d, %d)' % (z0, N, D, R, leg(D, 31), leg(R, 31), leg(D * R, 31)))

found = {'D': [], 'R': [], 'DR': []}
rows = 0
for p in sp.primerange(29, 3000):
    if p % 15 != 1 or p in (2, 3, 5, 11, 37):
        continue
    js = sp.sqrt_mod(-11 % p, p, all_roots=True)
    ws = [x for x in range(p) if (x ** 3 + 2 * x * x - 4 * x - 4) % p == 0]
    zs = [x for x in range(2, p) if pow(x, 15, p) == 1 and pow(x, 3, p) != 1 and pow(x, 5, p) != 1]
    if not (js and ws and zs):
        continue
    for z0 in zs:
        for j0 in js:
            for w0 in ws:
                try:
                    N, D, R = place_values(p, z0, j0, w0)
                except ValueError:
                    continue
                if N == 0 or D == 0 or R == 0:
                    continue
                rows += 1
                for name, val in (('D', D), ('R', R), ('DR', D * R)):
                    if leg(val, p) == -1 and len(found[name]) < 3:
                        found[name].append((p, z0, j0, w0))
print('degree-one places examined (31 <= p < 3000):', rows)
for name in ('D', 'R', 'DR'):
    print('  places where %s is a non-residue (first three): %s' % (name, found[name]))
ok = all(found[nm] for nm in found)
print('=> D, R, DR are all non-squares in F; [F(s,h):F] = 4 and [E\':Q] = 192:', ok)
