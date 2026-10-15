#!/usr/bin/env python3
"""S36 referee audit, part G: exact check of the seventeen-square packing of the eleven-squares paper, Appendix D.

Built directly from the paper's (D.1), (D.2) and Table 12.  Disjointness of interiors is tested by a method different
from the session's separating-axis test: the exact area of the intersection of every pair of squares is computed by
Sutherland-Hodgman clipping in rational arithmetic (interiors are disjoint iff the intersection has zero area, the
squares being convex with nonempty interior).  Also checks containment, the parameter identities, (D.4), the two
supporting-contact identities (D.5) and the side sum.
"""
from fractions import Fraction as Fr
from itertools import combinations

c, s = Fr(143, 145), Fr(24, 145)
D0 = 2306608
a, h, b, r = (Fr(x, D0) for x in (1016863, 449731, 567132, 390485))
print('c^2 + s^2 = 1:', c * c + s * s == 1, '; b = a - h:', b == a - h, '; r(c + s) = h:', r * (c + s) == h)
d = (1 - a) / 2
x0, y0 = 1 - 3 * h, 1 - a
D4 = [Fr(1, 4) < a < Fr(1, 2), 0 < h < b < Fr(3, 8), 1 - 3 * h > Fr(1, 4), Fr(3, 4) < a + 2 * b < 1, a + 2 * h < 1]
print('(D.4) inequalities:', D4)


def U(x, y, l):
    return [(x, y), (x + l, y), (x + l, y + l), (x, y + l)], l


sq = {}
for j in range(4):
    sq[f'L{j}'] = U(0, Fr(j, 4), Fr(1, 4))
for j in range(3):
    sq[f'B{j}'] = U(Fr(1, 4), j * (a - Fr(1, 4)), a - Fr(1, 4))
for i in range(2):
    for j in range(2):
        sq[f'C{i}{j}'] = U(a + i * d, j * d, d)
sq['D'] = U(Fr(1, 4), Fr(1, 4) + 2 * b, Fr(3, 4) - 2 * b)
for i in range(2):
    sq[f'E{i}'] = U(1 - 2 * b + i * b, 1 - b, b)
    sq[f'F{i}'] = U(1 - 2 * h + i * h, 1 - a, h)
base = (x0 + r * s, y0)
sq['R'] = ([(base[0] + v * c - w * s, base[1] + v * s + w * c) for (v, w) in ((0, 0), (r, 0), (r, r), (0, r))], r)


def area(P):
    return sum(P[i][0] * P[(i + 1) % len(P)][1] - P[(i + 1) % len(P)][0] * P[i][1] for i in range(len(P))) / 2


def clip(subject, clipper):
    """Sutherland-Hodgman: intersection of convex polygons (both counter-clockwise), exact."""
    out = subject
    for i in range(len(clipper)):
        A, B = clipper[i], clipper[(i + 1) % len(clipper)]
        inp = out
        out = []
        if not inp:
            break
        def side(P):
            return (B[0] - A[0]) * (P[1] - A[1]) - (B[1] - A[1]) * (P[0] - A[0])
        for j in range(len(inp)):
            P, Q = inp[j], inp[(j + 1) % len(inp)]
            sp, sq_ = side(P), side(Q)
            if sp >= 0:
                out.append(P)
            if (sp > 0 and sq_ < 0) or (sp < 0 and sq_ > 0):
                t = sp / (sp - sq_)
                out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    return out


for k, (V, l) in sq.items():
    assert area(V) == l * l and l > 0, k
inside = all(0 <= x <= 1 and 0 <= y <= 1 for V, _ in sq.values() for (x, y) in V)
overl = []
touching = 0
for p_, q_ in combinations(sq, 2):
    I = clip(sq[p_][0], sq[q_][0])
    A = area(I) if len(I) >= 3 else Fr(0)
    if A != 0:
        overl.append((p_, q_, A))
    elif len(I) >= 1:
        touching += 1
total = sum(l for _, l in sq.values())
print('squares:', len(sq), '; all inside [0,1]^2:', inside)
print('pairs with positive-area intersection:', overl, '; pairs touching (boundary contact):', touching)
print('total side =', total, '= 4 +', total - 4, '; > 4:', total > 4)
P = (a, 3 * a - Fr(3, 4))
Qp = (1 - 2 * b, Fr(1, 4) + 2 * b)
e1 = c * a + s * (3 * a - Fr(3, 4)) == c * x0 + s * y0 + c * r * s
e2 = -s * (1 - 2 * b) + c * (Fr(1, 4) + 2 * b) == -s * x0 + c * y0 + r * c * c
print('(D.5) supporting-contact identities:', e1, e2)
# Corollary D.2 constant: (4/k)(4 + 101/1153304) + (k^2 - 16)/k - k = 404/(1153304 k) = 101/(288326 k)
print('Corollary D.2 surplus coefficient: 4 * 101/1153304 =', Fr(4 * 101, 1153304))
# minimal clearance of R to the B column and D (should be exactly 0: contact) and to the others
def proj(V, ax):
    vals = [x * ax[0] + y * ax[1] for x, y in V]
    return min(vals), max(vals)
RB = min(proj(sq['R'][0], (c, s))[0] - proj(sq[k][0], (c, s))[1] for k in ('B0', 'B1', 'B2'))
RD = proj(sq['D'][0], (-s, c))[0] - proj(sq['R'][0], (-s, c))[1]
print('R vs B column clearance along (c, s):', RB, '; D vs R clearance along (-s, c):', RD)
