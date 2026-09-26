"""procgen_smallgraph_20260926_fence6.py -- T2.H (optional): numerical fence optimizer.

Session collatz-procgen-20260922, lane "smallgraph" (2026-09-26).  Called by the runner.

Configurations of the form "outer k-gon with unit sides + unit chords whose ends lie on sides or on
earlier chords" (T-junctions).  Fields are computed combinatorially (half-edge faces), each field is
constrained to area <= 1, fences to length 1, and the total area is maximized by multi-start SLSQP.
This is a numerical search: nothing here is claimed as a new record unless it beats the table,
and then only after an exact re-verification (none occurred).
"""
import math

import numpy as np
from scipy.optimize import minimize

from procgen_smallgraph_20260926_lib import check


def faces(points, fences):
    """points: list of (x,y); fences: list of (i, j, interior_ids).  Returns list of face areas
    (positive = bounded, counterclockwise) or None if the half-edge structure is degenerate."""
    pieces = set()
    for i, j, inner in fences:
        a, b = np.array(points[i]), np.array(points[j])
        d = b - a
        seq = sorted([i, j] + list(inner), key=lambda k: float(np.dot(np.array(points[k]) - a, d)))
        for p, q in zip(seq, seq[1:]):
            if p == q:
                return None
            pieces.add((min(p, q), max(p, q)))
    out = {}
    for u, v in pieces:
        out.setdefault(u, []).append(v)
        out.setdefault(v, []).append(u)
    for u in out:
        pu = points[u]
        out[u].sort(key=lambda w: math.atan2(points[w][1] - pu[1], points[w][0] - pu[0]))
    visited, areas = set(), []
    for u in out:
        for v in out[u]:
            if (u, v) in visited:
                continue
            face, a, b, guard = [], u, v, 0
            while (a, b) not in visited:
                visited.add((a, b))
                face.append(a)
                lst = out[b]
                c = lst[(lst.index(a) - 1) % len(lst)]
                a, b = b, c
                guard += 1
                if guard > 1000:
                    return None
            s = 0.0
            for k in range(len(face)):
                p, q = points[face[k]], points[face[(k + 1) % len(face)]]
                s += p[0] * q[1] - q[0] * p[1]
            areas.append(s / 2)
    return areas


def proper_cross(p1, p2, p3, p4, eps=1e-9):
    def o(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    d1, d2, d3, d4 = o(p3, p4, p1), o(p3, p4, p2), o(p1, p2, p3), o(p1, p2, p4)
    return ((d1 > eps and d2 < -eps) or (d1 < -eps and d2 > eps)) and ((d3 > eps and d4 < -eps) or (d3 < -eps and d4 > eps))


class Topology:
    """outer k-gon (unit sides) + chords.  chord spec: (host1, host2) where a host is ('s', side) or ('c', chord)."""

    def __init__(self, k, chords):
        self.k, self.chords = k, chords
        self.nv = (k - 1) + 2 * len(chords)   # k-1 side angles + one parameter per chord end

    def build(self, x):
        k = self.k
        ph = x[:k - 1]
        V = [np.zeros(2)]
        for i in range(k - 1):
            V.append(V[-1] + np.array([math.cos(ph[i]), math.sin(ph[i])]))
        pts = [tuple(v) for v in V]
        segs = [(i, (i + 1) % k) for i in range(k)]          # host segments: sides then chords
        inner = [[] for _ in range(k)]
        chord_ends = []
        for ci, (h1, h2) in enumerate(self.chords):
            ends = []
            for h, t in ((h1, x[k - 1 + 2 * ci]), (h2, x[k - 1 + 2 * ci + 1])):
                hid = h[1] if h[0] == 's' else k + h[1]
                a, b = np.array(pts[segs[hid][0]]), np.array(pts[segs[hid][1]])
                pts.append(tuple(a + t * (b - a)))
                inner[hid].append(len(pts) - 1)
                ends.append(len(pts) - 1)
            segs.append((ends[0], ends[1]))
            inner.append([])
            chord_ends.append(ends)
        fences = [(segs[i][0], segs[i][1], inner[i]) for i in range(len(segs))]
        return pts, fences

    def lengths(self, x):
        """closing side and chords (the first k-1 sides are unit by construction)"""
        pts, fences = self.build(x)
        return np.array([math.dist(pts[i], pts[j]) - 1.0 for i, j, _ in fences[self.k - 1:]])

    def valid(self, x):
        pts, fences = self.build(x)
        segs = [(pts[i], pts[j]) for i, j, _ in fences]
        for a in range(len(segs)):
            for b in range(a + 1, len(segs)):
                if proper_cross(*segs[a], *segs[b]):
                    return False
        return True


def optimise(top, seeds, nstart, rng):
    best = (-1.0, None)
    k = top.k
    for _ in range(nstart):
        base = rng.uniform(0, 2 * math.pi)
        x0 = np.concatenate([base + np.cumsum(np.full(k - 1, 2 * math.pi / k)) + rng.normal(0, 0.25, k - 1),
                             rng.uniform(0.15, 0.85, 2 * len(top.chords))])

        def area_list(x):
            pts, fences = top.build(x)
            A = faces(pts, fences)
            return A

        def obj(x):
            A = area_list(x)
            if A is None:
                return 10.0
            return -sum(a for a in A if a > 0)

        def cons_face(x):
            A = area_list(x)
            if A is None:
                return -np.ones(len(top.chords) + 1)
            pos = sorted([a for a in A if a > 0], reverse=True)
            pos = (pos + [0.0] * (len(top.chords) + 1))[:len(top.chords) + 1]
            return np.array([1.0 - a for a in pos])
        cons = [{'type': 'eq', 'fun': top.lengths},
                {'type': 'ineq', 'fun': cons_face},
                {'type': 'ineq', 'fun': lambda x: np.concatenate([x[k - 1:] - 0.02, 0.98 - x[k - 1:]])}]
        try:
            r = minimize(obj, x0, constraints=cons, method='SLSQP', options={'maxiter': 400, 'ftol': 1e-13})
        except Exception:
            continue
        x = r.x
        if np.max(np.abs(top.lengths(x))) > 1e-9 or not top.valid(x):
            continue
        A = area_list(x)
        if A is None:
            continue
        pos = [a for a in A if a > 1e-9]
        neg = [a for a in A if a < -1e-9]
        if len(neg) != 1 or len(pos) != len(top.chords) + 1 or max(pos) > 1 + 1e-9:
            continue
        tot = sum(pos)
        if abs(tot + neg[0]) > 1e-9:       # union of fields = outer polygon
            continue
        if tot > best[0]:
            best = (tot, (x.copy(), sorted(pos)))
    return best


def t2h():
    print('== T2.H (optional) numerical fence optimizer ==')
    rng = np.random.default_rng(20260926)
    # n = 6: pentagon + one chord between two sides (the record's topology), all side pairs
    best6 = (-1, None, None)
    for a in range(5):
        for b in range(a + 1, 5):
            top = Topology(5, [(('s', a), ('s', b))])
            v, info = optimise(top, None, 25, rng)
            if v > best6[0]:
                best6 = (v, (a, b), info)
    v, sides, info = best6
    a, b = sides
    adjacent = (b - a) % 5 in (1, 4)
    check('T2.H1 n=6: pentagon of unit sides + unit chord, optimum reproduces the record 1.47585',
          1.47585 - 1e-6 <= v <= 1.47586 and not adjacent,
          f'best {v:.7f} (fields {info[1][0]:.6f} + {info[1][1]:.6f}), chord between sides {a} and {b} (non-adjacent); '
          f'record 1.47585+ (Lagache 2026) matched to the printed digits, not beaten')
    # n = 8: hexagon + two chords (second chord may end on the first)
    tops = []
    for a in range(6):
        for b in range(a + 2, 6):
            if (a, b) == (0, 5):
                continue
            for c in list(range(6)) + ['c0']:
                for d in range(6):
                    h1 = ('c', 0) if c == 'c0' else ('s', c)
                    if h1 == ('s', d):
                        continue
                    tops.append(Topology(6, [(('s', a), ('s', b)), (h1, ('s', d))]))
    rng8 = np.random.default_rng(8)
    sample = [tops[i] for i in sorted(rng8.choice(len(tops), size=min(40, len(tops)), replace=False))]
    best8 = -1
    for top in sample:
        v8, info8 = optimise(top, None, 6, rng8)
        best8 = max(best8, v8)
    check('T2.H2 n=8: hexagon of unit sides + two unit chords (40 random combinatorial types, 6 starts each)',
          best8 < 2.10306 + 1e-4, f'best found {best8:.5f} < record 2.10306 (Gomes 2026); no improvement '
          '(the record uses a fence that runs from the boundary into the interior, outside this family)')
