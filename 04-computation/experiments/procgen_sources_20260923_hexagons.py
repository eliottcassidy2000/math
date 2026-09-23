#!/usr/bin/env python3
"""procgen_sources_20260923_hexagons.py

Exact check (arithmetic in Q(sqrt 3), no floating point in any decision) that
seven non-overlapping unit regular hexagons fit in a regular hexagon of side
5/sqrt(3), and of the rotated-flower family N_k = 3k(k-1)+1 hexagons in side
(3k-1)/sqrt(3).

Configuration (Morandi's record as displayed on Friedman's Packing Center,
"Hexagons in Hexagons", n = 6-7): the honeycomb flower of 7 "pointy-top" unit
hexagons (centres 0 and sqrt(3) e^(i pi k/3)) inside a "flat-top" container
(vertices at angles 0, 60, ..., 300 degrees), i.e. the container is rotated by
30 degrees relative to the tiles.

Checks:
  H1 containment: every tile vertex v satisfies n_k . v <= apothem for the six
     outward container normals n_k (convexity => tile inside container);
  H2 disjoint interiors: separating-axis test on every pair of tiles, exact;
  H3 contacts: which tile vertices lie exactly on the container boundary;
  H4 rigidity under shrinking: for side 5/sqrt3 - t (t > 0) this configuration
     (scaled container, same tiles) is infeasible, exactly;
  H5 the area bound s >= sqrt(7) and the gap;
  H6 the family k = 2, 3, 4 (7, 19, 37 tiles in side 5/sqrt3, 8/sqrt3, 11/sqrt3),
     compared with the values listed on Friedman's page (numerically).
"""
import sys
from fractions import Fraction as Fr
from itertools import combinations


def say(*a):
    print(*a)
    sys.stdout.flush()


class Q3:
    """a + b sqrt(3), a, b rational."""
    __slots__ = ('a', 'b')

    def __init__(self, a, b=0):
        self.a = Fr(a)
        self.b = Fr(b)

    def __add__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return Q3(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return Q3(-self.a, -self.b)

    def __sub__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return Q3(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return Q3(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, Q3) else Q3(o)
        return Q3(self.a * o.a + 3 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def sign(self):
        a, b = self.a, self.b
        if a >= 0 and b >= 0:
            return 0 if (a == 0 and b == 0) else 1
        if a <= 0 and b <= 0:
            return -1
        # opposite signs: compare a^2 with 3 b^2
        d = a * a - 3 * b * b
        if d == 0:
            return 0  # impossible unless a = b = 0 (sqrt 3 irrational)
        return (1 if a > 0 else -1) if d > 0 else (1 if b > 0 else -1)

    def __float__(self):
        return float(self.a) + float(self.b) * 3 ** 0.5

    def __repr__(self):
        return f'({self.a}+{self.b}r3)'

    def eq(self, o):
        return (self - o).sign() == 0

    def le(self, o):
        return (self - o).sign() <= 0

    def lt(self, o):
        return (self - o).sign() < 0


R3 = Q3(0, 1)
HALF = Fr(1, 2)

# pointy-top unit hexagon vertices: angles 30 + 60 k
TILE = [(Q3(0, HALF), Q3(HALF)), (Q3(0), Q3(1)), (Q3(0, -HALF), Q3(HALF)),
        (Q3(0, -HALF), Q3(-HALF)), (Q3(0), Q3(-1)), (Q3(0, HALF), Q3(-HALF))]
# outward normals of a flat-top container (angles 30 + 60 k), unit length
NORMALS = [(Q3(0, HALF), Q3(HALF)), (Q3(0), Q3(1)), (Q3(0, -HALF), Q3(HALF)),
           (Q3(0, -HALF), Q3(-HALF)), (Q3(0), Q3(-1)), (Q3(0, HALF), Q3(-HALF))]
# edge normals of a pointy-top tile: angles 0, 60, 120 (and opposite)
TILE_AXES = [(Q3(1), Q3(0)), (Q3(HALF), Q3(0, HALF)), (Q3(-HALF), Q3(0, HALF))]


def flower_centres(k):
    """centres of the k-ring honeycomb flower of pointy-top unit hexagons.
    Lattice basis: sqrt3*(1,0) and sqrt3*(1/2, sqrt3/2) = (sqrt3/2, 3/2)."""
    cs = []
    for i in range(-(k - 1), k):
        for j in range(-(k - 1), k):
            if abs(i) <= k - 1 and abs(j) <= k - 1 and abs(i + j) <= k - 1:
                x = Q3(0, i) + Q3(0, Fr(j, 2))
                y = Q3(Fr(3 * j, 2))
                cs.append((x, y))
    return cs


def tile_vertices(c):
    return [(c[0] + vx, c[1] + vy) for vx, vy in TILE]


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1]


def check(k, side):
    """side is a Q3. Returns (inside, disjoint, contacts, n_tiles)."""
    apothem = side * R3 * Q3(HALF)  # s sqrt3 / 2
    cs = flower_centres(k)
    inside = True
    contacts = 0
    for c in cs:
        for v in tile_vertices(c):
            for n in NORMALS:
                d = dot(n, v)
                if not d.le(apothem):
                    inside = False
                if d.eq(apothem):
                    contacts += 1
    disjoint = True
    for c1, c2 in combinations(cs, 2):
        V1, V2 = tile_vertices(c1), tile_vertices(c2)
        sep = False
        for ax in TILE_AXES:
            p1 = [dot(ax, v) for v in V1]
            p2 = [dot(ax, v) for v in V2]
            mx1 = max(p1, key=float)
            mn1 = min(p1, key=float)
            mx2 = max(p2, key=float)
            mn2 = min(p2, key=float)
            # exact confirmation of the float-chosen extremes
            assert all(x.le(mx1) for x in p1) and all(mn1.le(x) for x in p1)
            assert all(x.le(mx2) for x in p2) and all(mn2.le(x) for x in p2)
            if mx1.le(mn2) or mx2.le(mn1):
                sep = True
                break
        if not sep:
            disjoint = False
    return inside, disjoint, contacts, len(cs)


def main():
    say('procgen_sources_20260923_hexagons.py')
    say('=' * 78)
    s7 = Q3(0, Fr(5, 3))  # 5/sqrt3 = 5 sqrt3 / 3
    inside, disjoint, contacts, n = check(2, s7)
    say(f'H1-H3  k = 2: {n} unit hexagons, container side 5/sqrt3 = {float(s7):.12f}')
    say(f'   every tile inside the container (exact): {inside}')
    say(f'   pairwise disjoint interiors (exact separating axes): {disjoint}')
    say(f'   (tile vertex, container edge) incidences with equality: {contacts}')
    # H4: shrink by t = 1/10^6 and by an exact tiny amount
    for t in (Fr(1, 10 ** 6), Fr(1, 10 ** 30)):
        ins, _, _, _ = check(2, s7 - Q3(t))
        say(f'H4  side 5/sqrt3 - {t}: this configuration fits: {ins}')
    # H5 area bound
    say(f'H5  area bound: 7 * (3 sqrt3/2) <= (3 sqrt3/2) s^2  =>  s >= sqrt7 = {7 ** 0.5:.12f}; '
        f'record/bound = {float(s7) / 7 ** 0.5:.6f}')
    # H6 family
    say('H6  rotated-flower family: N_k = 3k(k-1)+1 tiles in side (3k-1)/sqrt3 = (3k-1) sqrt3/3')
    friedman = {7: '6-7: 5/sqrt3 = 2.886+', 19: '18-19: 8/sqrt3 = 4.618+',
                37: '(page ends at n = 34: 6.35085+)'}
    for k in (2, 3, 4):
        side = Q3(0, Fr(3 * k - 1, 3))
        ins, dis, con, n = check(k, side)
        tighter, _, _, _ = check(k, side - Q3(Fr(1, 10 ** 9)))
        say(f'   k = {k}: N = {n:2d}, side = {float(side):.9f}, inside {ins}, disjoint {dis}, contacts {con}, '
            f'fits at side - 1e-9: {tighter}; Friedman: {friedman[n]}')
    say(f'   11/sqrt3 = {11 / 3 ** 0.5:.9f} against the listed n = 34 value 6.35085+')
    say('=' * 78)
    say('hexagons: done')


if __name__ == '__main__':
    main()
