"""procgen_fence2_20260926_constr.py -- constructions and families for Part C of lane "fence2"
(session collatz-procgen-20260922, 2026-09-26).

- pythagorean(a): the Pythagorean pinwheel tiling as a torus configuration (2 fences per period).
- flipped_square(s): the pentagon with sides (1, s, 1, 1-s, 1), exact rational coordinates, area 1.
- pent_roundabout_family(c): the O-pinwheel pentagon whose three corners are partnered by a tiny
  I-triangle (circumdiameter c) and two by mirror pentagons; maximise a_P + a_T with a_P <= 1.
"""
import itertools
import math
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import minimize


def pythagorean(a):
    """representative fences and lattice of the Pythagorean tiling with squares a and b = 1 - a."""
    a = Fr(a)
    b = 1 - a
    fences = [((-b, a), (a, a)), ((a, Fr(0)), (a, Fr(1)))]
    lattice = ((a, b), (-b, a))
    return fences, lattice


def reflect(p, a, b):
    """reflection of p in the line through a, b (exact)."""
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)
    fx, fy = a[0] + t * dx, a[1] + t * dy
    return (2 * fx - p[0], 2 * fy - p[1])


def flipped_square(s):
    """unit square with the corner triangle (1,0),(1,1),(1-s,1) replaced by the congruent right
    triangle on the other side of its hypotenuse with the legs exchanged (legs s at (1,0), 1 at
    (1-s,1)): the pentagon (0,0),(1,0),c2,(1-s,1),(0,1) with sides (1, s, 1, 1-s, 1), area 1."""
    s = Fr(s)
    c0, c1, c3, c4 = (Fr(0), Fr(0)), (Fr(1), Fr(0)), (1 - s, Fr(1)), (Fr(0), Fr(1))
    c2 = reflect((1 - s, Fr(0)), c1, c3)
    P = [c0, c1, c2, c3, c4]
    A2 = sum(P[i][0] * P[(i + 1) % 5][1] - P[(i + 1) % 5][0] * P[i][1] for i in range(5))
    l2 = [(P[(i + 1) % 5][0] - P[i][0]) ** 2 + (P[(i + 1) % 5][1] - P[i][1]) ** 2 for i in range(5)]
    # convexity: all cross products of consecutive edges positive
    conv = all(((P[(i + 1) % 5][0] - P[i][0]) * (P[(i + 2) % 5][1] - P[(i + 1) % 5][1]) -
                (P[(i + 1) % 5][1] - P[i][1]) * (P[(i + 2) % 5][0] - P[(i + 1) % 5][0])) > 0 for i in range(5))
    return dict(P=P, area=A2 / 2, len2=l2, convex=conv)


def _pent(ang, sides):
    d = 0.0
    P = [np.zeros(2)]
    for i in range(5):
        P.append(P[-1] + sides[i] * np.array([math.cos(d), math.sin(d)]))
        d += math.pi - ang[(i + 1) % 5]
    P = np.array(P)
    return 0.5 * np.sum(P[:-1, 0] * P[1:, 1] - P[1:, 0] * P[:-1, 1])


def pent_roundabout_family(c, starts=40, seed=0, Pset=(2, 4)):
    """max of a_P + a_T over the family with fixed triangle circumdiameter c (sides t_k = c sin(opp.)):
    O-pinwheel pentagon, corners T=[0,1,3] partnered by the tiny I-triangle (angle 180 - alpha_k,
    preceding side 1 - t_k), corners Pset by mirror pentagons (angles th, 180-th; sides s, 1-s)."""
    rng = np.random.default_rng(seed)
    T = [q for q in range(5) if q not in Pset]
    best = None
    for perm in itertools.permutations(range(3)):
        for which in (0, 1):
            def build(v):
                a1, a2, th = v
                al = [a1, a2, 180 - a1 - a2]
                t = [c * math.sin(math.radians(al[(k + 2) % 3])) for k in range(3)]
                ang = [0.0] * 5
                sides = [None] * 5
                for idx, cc in enumerate(T):
                    k = perm[idx]
                    ang[cc] = math.radians(180 - al[k])
                    sides[(cc - 1) % 5] = 1 - t[k]
                ang[Pset[which]] = math.radians(th)
                ang[Pset[1 - which]] = math.radians(180 - th)
                d = 0.0
                dirs = []
                for i in range(5):
                    dirs.append(d)
                    d += math.pi - ang[(i + 1) % 5]
                u = [np.array([math.cos(x), math.sin(x)]) for x in dirs]
                fr = [(p - 1) % 5 for p in Pset]
                rhs = -sum(sides[i] * u[i] for i in range(5) if sides[i] is not None)
                M = np.array([u[fr[0]], u[fr[1]]]).T
                if abs(np.linalg.det(M)) < 1e-12:
                    return None
                x = np.linalg.solve(M, rhs)
                sides[fr[0]], sides[fr[1]] = x
                aT = 0.5 * c * c * math.sin(math.radians(al[0])) * math.sin(math.radians(al[1])) * \
                    math.sin(math.radians(al[2]))
                aP = _pent(ang, sides)
                return aP, aT, sides, fr

            def obj(v):
                p = build(v)
                return 10.0 if p is None else -(p[0] + p[1])

            def eq(v):
                p = build(v)
                return 10.0 if p is None else p[2][p[3][0]] + p[2][p[3][1]] - 1

            def ineq(v):
                p = build(v)
                if p is None:
                    return -1.0
                return min(1 - p[0], min(p[2]) - 1e-6, 1 - max(p[2]), v[0] - 0.01, v[1] - 0.01,
                           179.99 - v[0] - v[1], v[2] - 0.01, 179.99 - v[2])
            for st in range(starts):
                v0 = [rng.uniform(5, 120), rng.uniform(5, 120), rng.uniform(20, 160)]
                if v0[0] + v0[1] > 170:
                    continue
                try:
                    r = minimize(obj, v0, method='SLSQP', constraints=[{'type': 'eq', 'fun': eq},
                                 {'type': 'ineq', 'fun': ineq}], options={'maxiter': 500, 'ftol': 1e-14})
                except Exception:
                    continue
                p = build(r.x)
                if p is None or abs(eq(r.x)) > 1e-9 or ineq(r.x) < -1e-9:
                    continue
                if best is None or p[0] + p[1] > best[0]:
                    best = (p[0] + p[1], p[0], p[1], r.x.copy(), perm, which)
    return best
