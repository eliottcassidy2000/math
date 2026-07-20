#!/usr/bin/env python3
"""
death-star-2026-07-19-S59m -- verify the claimed Jacobian Conjecture counterexample
F: C^3 -> C^3,
  F1 = (1+xy)^3 z + y^2 (1+xy)(4+3xy)
  F2 = y + 3x(1+xy)^2 z + 3xy^2 (4+3xy)
  F3 = 2x - 3x^2 y - x^3 z
Claims: det JF = -2 (constant); F maps (0,0,-1/4), (1,-3/2,13/2), (-1,3/2,13/2)
all to (-1/4, 0, 0).  Exact arithmetic, no dependencies.
"""
from fractions import Fraction
from itertools import permutations

def pmul(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = (ka[0]+kb[0], ka[1]+kb[1], ka[2]+kb[2])
            r[k] = r.get(k, 0) + ca*cb
    return {k: c for k, c in r.items() if c != 0}

def padd(*ps):
    r = {}
    for p in ps:
        for k, c in p.items():
            r[k] = r.get(k, 0) + c
    return {k: c for k, c in r.items() if c != 0}

def pscale(p, s):
    return {k: c*s for k, c in p.items() if c*s != 0}

def pdiff(p, i):
    r = {}
    for k, c in p.items():
        if k[i] > 0:
            k2 = list(k); k2[i] -= 1
            r[tuple(k2)] = r.get(tuple(k2), 0) + c*k[i]
    return {k: c for k, c in r.items() if c != 0}

def peval(p, pt):
    s = Fraction(0)
    for (i, j, k), c in p.items():
        s += Fraction(c) * pt[0]**i * pt[1]**j * pt[2]**k
    return s

X = {(1,0,0): 1}; Y = {(0,1,0): 1}; Z = {(0,0,1): 1}; ONE = {(0,0,0): 1}
U = padd(ONE, pmul(X, Y))                      # 1 + xy
U2 = pmul(U, U); U3 = pmul(U2, U)
W = padd(pscale(ONE, 4), pscale(pmul(X, Y), 3))  # 4 + 3xy

F1 = padd(pmul(U3, Z), pmul(pmul(pmul(Y, Y), U), W))
F2 = padd(Y, pscale(pmul(pmul(X, U2), Z), 3), pscale(pmul(pmul(X, pmul(Y, Y)), W), 3))
F3 = padd(pscale(X, 2), pscale(pmul(pmul(X, X), Y), -3), pscale(pmul(pmul(pmul(X, X), X), Z), -1))
F = [F1, F2, F3]

J = [[pdiff(Fi, v) for v in range(3)] for Fi in F]
det = {}
for perm in permutations(range(3)):
    sgn = 1 if perm in [(0,1,2),(1,2,0),(2,0,1)] else -1
    t = pscale(pmul(pmul(J[0][perm[0]], J[1][perm[1]]), J[2][perm[2]]), sgn)
    det = padd(det, t)
print("det JF =", det, " (constant -2:", det == {(0,0,0): -2}, ")")

pts = [(Fraction(0), Fraction(0), Fraction(-1,4)),
       (Fraction(1), Fraction(-3,2), Fraction(13,2)),
       (Fraction(-1), Fraction(3,2), Fraction(13,2))]
target = (Fraction(-1,4), Fraction(0), Fraction(0))
allok = True
for pt in pts:
    img = tuple(peval(Fi, pt) for Fi in F)
    ok = img == target
    allok &= ok
    print(f"F({pt[0]},{pt[1]},{pt[2]}) = {img}  -> target: {ok}")
print("degrees:", [max(sum(k) for k in Fi) for Fi in F])
print("VERDICT: constant det -2 AND triple collision =", (det == {(0,0,0): -2}) and allok)
