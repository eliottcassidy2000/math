#!/usr/bin/env python3
"""Exact rational verification of the claimed counterexample to Erdos Problem 106 (f(17) > 4).

Source (UNREFEREED): "The Fibonacci geometry of eleven squares" (PINGYOU LTD, 2026), Appendix D, Theorem D.1, building
on Sprite143's seventeen-square layout: sixteen axis-parallel squares (Table 12) and one square R rotated by the
Pythagorean angle (c, s) = (143/145, 24/145), with (a, h, b, r) = (1016863, 449731, 567132, 390485)/2306608.
Claim: the seventeen squares lie in [0,1]^2 with pairwise disjoint interiors and total side 4 + 101/1153304.
This script checks, in exact rational arithmetic:
  (1) the parameter identities c^2 + s^2 = 1, b = a - h, r(c + s) = h;
  (2) containment of every square in [0,1]^2 (vertices of R included);
  (3) pairwise interior-disjointness by the separating-axis theorem for convex polygons (axes: the edge normals of
      both squares; interiors are disjoint iff some axis has projection intervals meeting in at most a point);
  (4) the side sum.
Session opus-2026-10-14-S36.
"""
from fractions import Fraction as Fr
from itertools import combinations

c, s = Fr(143, 145), Fr(24, 145)
D = 2306608
a, h, b, r = Fr(1016863, D), Fr(449731, D), Fr(567132, D), Fr(390485, D)
assert c * c + s * s == 1 and b == a - h and r * (c + s) == h, 'parameter identities'
d = (1 - a) / 2
x0, y0 = 1 - 3 * h, 1 - a


def U(x, y, l):
    return ('axis', [(x, y), (x + l, y), (x + l, y + l), (x, y + l)], l)


squares = {}
for j in range(4):
    squares[f'L{j}'] = U(0, Fr(j, 4), Fr(1, 4))
for j in range(3):
    squares[f'B{j}'] = U(Fr(1, 4), j * (a - Fr(1, 4)), a - Fr(1, 4))
for i in range(2):
    for j in range(2):
        squares[f'C{i}{j}'] = U(a + i * d, j * d, d)
squares['D'] = U(Fr(1, 4), Fr(1, 4) + 2 * b, Fr(3, 4) - 2 * b)
for i in range(2):
    squares[f'E{i}'] = U(1 - 2 * b + i * b, 1 - b, b)
    squares[f'F{i}'] = U(1 - 2 * h + i * h, 1 - a, h)
o = (x0 + r * s, y0)
Rv = [(o[0] + v * c - w * s, o[1] + v * s + w * c) for (v, w) in [(0, 0), (r, 0), (r, r), (0, r)]]
squares['R'] = ('rot', Rv, r)


def axes(sq):
    kind, V, l = sq
    if kind == 'axis':
        return [(1, 0), (0, 1)]
    return [(c, s), (-s, c)]


def proj(V, ax):
    vals = [p[0] * ax[0] + p[1] * ax[1] for p in V]
    return min(vals), max(vals)


def interiors_disjoint(A, B):
    for ax in axes(A) + axes(B):
        a0, a1 = proj(A[1], ax)
        b0, b1 = proj(B[1], ax)
        if a1 <= b0 or b1 <= a0:
            return True
    return False


ok_side = all(sq[2] > 0 for sq in squares.values())
ok_contain = all(0 <= x <= 1 and 0 <= y <= 1 for sq in squares.values() for (x, y) in sq[1])
bad = [(p, q) for p, q in combinations(squares, 2) if not interiors_disjoint(squares[p], squares[q])]
total = sum(sq[2] for sq in squares.values())
print('number of squares:', len(squares))
print('all side lengths positive:', ok_side)
print('all squares inside [0,1]^2 (exact):', ok_contain)
print('pairs with overlapping interiors (separating-axis test, exact):', bad)
print('total side length =', total, '= 4 +', total - 4, ' ; exceeds 4:', total > 4)
print('R vertices:', [(str(x), str(y)) for x, y in Rv])
