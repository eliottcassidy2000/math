#!/usr/bin/env python3
"""THM-4604 checks: Collatz U-words as upper-triangular matrices G_w = [[3^|w|, B_w],[0, 2^sum w]] (chronological product):
cocycle identity, coboundary on cyclic submonoids (cycle points), trace formula, Cayley cubic x^2+y^2+z^2-xyz = 4 (exact)."""
import random
from fractions import Fraction as Fr
def G(word):
    a, b, d = 1, 0, 1          # [[a, b],[0, d]]
    for x in word:
        # new = G_x * old, G_x = [[3, 1],[0, 2^x]]
        a, b, d = 3*a, 3*b + d, d * 2**x
    return a, b, d
rnd = random.Random(4)
CH = 0
for _ in range(2000):
    u = [rnd.randint(1, 6) for _ in range(rnd.randint(1, 8))]
    v = [rnd.randint(1, 6) for _ in range(rnd.randint(1, 8))]
    au, bu, du = G(u); av, bv, dv = G(v); a, b, d = G(u + v)
    assert (a, d) == (au*av, du*dv)
    assert b == 3**len(v) * bu + 2**sum(u) * bv          # cocycle identity
    assert b > 0
    # Cayley cubic with SL2-normalised traces (exact): x^2 = tr^2/det etc., xyz = tr_u tr_v tr_uv / det_uv
    tu, tv, t = au + du, av + dv, a + d
    x2, y2, z2 = Fr(tu*tu, au*du), Fr(tv*tv, av*dv), Fr(t*t, a*d)
    xyz = Fr(tu*tv*t, a*d)
    assert x2 + y2 + z2 - xyz == 4
    # coboundary on <u> when 2^sum != 3^len
    if 2**sum(u) != 3**len(u):
        c = Fr(bu, 2**sum(u) - 3**len(u))
        for j in (2, 3):
            aj, bj, dj = G(u*j)
            assert bj == c * (2**(j*sum(u)) - 3**(j*len(u)))
    CH += 1
print(f"cocycle identity, positivity, Cayley cubic x^2+y^2+z^2-xyz=4 and cyclic coboundaries: {CH} random word pairs")
# the odd/even generating pair
O = (3, 1, 2); E = (1, 0, 2)
x2, y2 = Fr(5*5, 6), Fr(3*3, 2)
z2 = Fr(7*7, 12); xyz = Fr(5*3*7, 12)
assert x2 + y2 + z2 - xyz == 4
print("odd/even pair (5/sqrt6, 3/sqrt2, 7/sqrt12) lies on the Cayley cubic")
print("ALL CHECKS PASSED")
