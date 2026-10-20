#!/usr/bin/env python3
"""p = 211 (p-1 = 210 = 2*3*5*7, p = 3 mod 4): the slice-mixed construction gives a tournament set that is traceless
(all power sums sigma_k, gcd(k,210) = 1, vanish) with trivial multiplier stabiliser, so the power-sum detector fails
even for p = 3 mod 4 once p-1 has three odd prime factors. Also p = 31 (= 2*3*5, one slice class) for contrast."""
from math import gcd
from sympy import primitive_root
def test(p, f):
    g = primitive_root(p); log = {pow(g, t, p): t for t in range(p - 1)}
    S = frozenset(x for x in range(1, p) if f(log[x]) == 1)
    assert all((x in S) != ((p - x) in S) for x in range(1, p)), "not a tournament set"
    units = [k for k in range(1, p - 1) if gcd(k, p - 1) == 1]
    traceless = all(sum(pow(s, k, p) for s in S) % p == 0 for k in units)
    st = [a for a in range(1, p) if frozenset(a * x % p for x in S) == S]
    return traceless, st
p = 211
phi = {z: {y: (1 if (y * (z + 2)) % 3 == 0 else -1) for y in range(5)} for z in range(7)}   # phi_z(y), nonconstant in y
psi = {z: {x: (1 if (x + z) % 3 == 0 else -1) for x in range(3)} for z in range(7)}         # psi_z(x), nonconstant in x
Z1 = {0, 1, 2}; Z2 = {3, 4, 5, 6}
def f211(t):
    t2, x, y, z = t % 2, t % 3, t % 5, t % 7
    A = phi[z][y] if z in Z1 else psi[z][x]
    return A if t2 == 0 else -A
print("p = 211 slice-mixed set: traceless =", test(p, f211)[0], ", multiplier stabiliser =", test(p, f211)[1])
# p = 31 analogue (Z_30 = Z_2 x Z_3 x Z_5): the only option is a single slice type; mixing is impossible, so test a
# generic U + V function and confirm it is never +-1 unless one part is constant: brute force all U, V in {0,+-1}
p = 31; cnt = 0; found = 0
import itertools
for U in itertools.product([-1, 0, 1], repeat=5):          # U(t5)
    for V in itertools.product([-1, 0, 1], repeat=3):      # V(t3)
        if all(U[y] + V[x] in (1, -1) for x in range(3) for y in range(5)):
            cnt += 1
            if len(set(U)) > 1 and len(set(V)) > 1: found += 1
print(f"p = 31: +-1-valued arrays U(t5)+V(t3): {cnt}; with both U and V nonconstant: {found} (the argument predicts 0)")
