#!/usr/bin/env python3
"""Audit B, item 4: exact recount of unit (sqrt 7) distances in the Eisenstein-integer sets of MISTAKE-575.
Eisenstein integers x + y*w with w = e^{i pi/3} (w^2 = w - 1), stored as integer pairs (x, y).
Norm N(x + y w) = x^2 + x y + y^2 (since |x + y w|^2 = x^2 + 2xy Re(w) + y^2 = x^2 + xy + y^2)."""
import itertools
from fractions import Fraction


def mul(a, b):  # (a0 + a1 w)(b0 + b1 w) = a0b0 + (a0b1 + a1b0) w + a1b1 w^2,  w^2 = w - 1
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0] + a[1] * b[1])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def norm(a):
    return a[0] ** 2 + a[0] * a[1] + a[1] ** 2


def conj(a):  # conj(w) = 1 - w
    return (a[0] + a[1], -a[1])


w = (0, 1)
one = (1, 0)
zero = (0, 0)
units = [one, w, mul(w, w), (-1, 0), mul((-1, 0), w), mul((-1, 0), mul(w, w))]
assert all(norm(u) == 1 for u in units) and len(set(units)) == 6
alpha = (2, 1)            # 2 + w
alphab = (3, -1)          # 3 - w
assert conj(alpha) == alphab and norm(alpha) == norm(alphab) == 7
# cos of the angle between alpha and alphab: Re(alpha * conj(alphab)) / 7
prod = mul(alpha, conj(alphab))
re2 = 2 * prod[0] + prod[1]          # 2*Re(x + y w) = 2x + y
print("alpha*conj(alphab) =", prod, "-> cos angle =", Fraction(re2, 2 * 7))

tri = [zero, one, w]
W6 = [zero] + units


def count(P, d2=7):
    P = list(P)
    assert len(set(P)) == len(P), "points not distinct"
    return sum(1 for p, q in itertools.combinations(P, 2) if norm(add(p, (-q[0], -q[1]))) == d2)


S21 = [add(mul(alpha, p), mul(alphab, q)) for p in tri for q in W6]
print("21-point set: distinct =", len(set(S21)), " pairs at distance sqrt7 =", count(S21))
S49 = [add(mul(alpha, p), mul(alphab, q)) for p in W6 for q in W6]
print("49-point set W6 x W6: distinct =", len(set(S49)), " pairs at distance sqrt7 =", count(S49), " 3N =", 3 * 49)
# extra (non-product) edges?
def product_edges(A, B):
    eA = sum(1 for p, q in itertools.combinations(A, 2) if norm(add(p, (-q[0], -q[1]))) == 1)
    eB = sum(1 for p, q in itertools.combinations(B, 2) if norm(add(p, (-q[0], -q[1]))) == 1)
    return eA * len(B) + eB * len(A)
print("generic product counts: tri x W6 =", product_edges(tri, W6), " W6 x W6 =", product_edges(W6, W6))
# all distance multiplicities in the 21-point set
from collections import Counter
print("distance^2 multiset (21 pts), top:", Counter(norm(add(p, (-q[0], -q[1]))) for p, q in itertools.combinations(S21, 2)).most_common(6))
