#!/usr/bin/env python3
"""Finite exact regression controls, not an all-level proof."""
from fractions import Fraction as F
from itertools import combinations
import json

DEPTH = 7


def require(predicate, label):
    if not predicate:
        raise RuntimeError(label)


def point(v):
    k, j = v
    h = F(1, 2**k)
    return (j + F(1, 2)) * h, h


def cosh_distance(p, q):
    x, y = p
    xx, yy = q
    return 1 + ((x-xx)**2 + (y-yy)**2) / (2*y*yy)


def adjacent(u, v):
    if u[0] > v[0]:
        u, v = v, u
    return v[0] == u[0] + 1 and v[1] // 2 == u[1]


def run():
    vertices = [(k, j) for k in range(DEPTH+1) for j in range(2**k)]
    points = {v: point(v) for v in vertices}
    edges = [(v, (v[0]+1, 2*v[1]+b))
             for v in vertices if v[0] < DEPTH for b in (0, 1)]
    require(len(set(points.values())) == len(vertices), "vertex injectivity")
    pair_count = equality_count = 0
    for u, v in combinations(vertices, 2):
        value = cosh_distance(points[u], points[v])
        require(value >= F(21, 16), "sharp separation")
        require((value == F(21, 16)) == adjacent(u, v), "equality iff edge")
        delta = abs(u[0]-v[0])
        lower = F(3, 2) if delta == 0 else (F(2**delta)+F(1, 2**delta))/2
        require(value >= lower, "same-level or height-ratio bound")
        equality_count += value == F(21, 16)
        pair_count += 1

    rectangles = {}
    for p, q in edges:
        a, h = points[p]
        x, y = points[q]
        s = 2*(q[1] % 2)-1
        require((x, y) == (a+s*h/4, h/2), "child offset")
        require(cosh_distance(points[p], points[q]) == F(21, 16), "edge length")
        center = a-s*11*h/8
        r2 = F(185, 64)*h*h
        require((a-center)**2+h*h == r2, "parent circle equation")
        require((x-center)**2+y*y == r2, "child circle equation")
        require(s*(a-center) > 0 and s*(x-center) > 0, "same monotone side")
        lo, hi = sorted((a, x))
        # Rational squared-height samples, with no square-root rounding.
        for t in (F(i, 8) for i in range(9)):
            xx = (1-t)*a+t*x
            yy2 = r2-(xx-center)**2
            require(h*h/4 <= yy2 <= h*h, "arc height-band samples")
            require(lo <= xx <= hi, "arc x-interval samples")
        rectangles[(p, q)] = (lo, hi, h/2, h)

    edge_pairs = incident_pairs = 0
    for e, f in combinations(edges, 2):
        if e[0][0] > f[0][0]:
            e, f = f, e
        k, l = e[0][0], f[0][0]
        e_lo, e_hi, e_bottom, e_top = rectangles[e]
        f_lo, f_hi, f_bottom, f_top = rectangles[f]
        incident = bool(set(e) & set(f))
        if k == l:
            if e[0] == f[0]:
                require(max(e_lo, f_lo) == min(e_hi, f_hi)
                        == points[e[0]][0], "sibling x-overlap only at parent")
                require(incident, "sibling incidence")
            else:
                require(e_hi < f_lo or f_hi < e_lo, "other parent x-disjoint")
                require(not incident, "different parents same-level incidence")
        elif l == k+1:
            require(e_bottom == f_top, "consecutive band boundary")
            require(incident == (e[1] == f[0]), "boundary endpoints incidence")
        else:
            require(f_top < e_bottom and not incident, "nonconsecutive bands")
        edge_pairs += 1
        incident_pairs += incident

    for k in range(1, DEPTH+1):
        left, right = (k, 2**(k-1)-1), (k, 2**(k-1))
        require(cosh_distance(points[left], points[right]) == F(3, 2),
                "root-split hyperbolic distance constant")
        require(points[right][0]-points[left][0] == F(1, 2**k),
                "Euclidean crowding")
        a, b = left, right
        steps = 0
        while a != b:
            a, b = (a[0]-1, a[1]//2), (b[0]-1, b[1]//2)
            steps += 2
        require(a == (0, 0) and steps == 2*k, "root-split graph distance")

    require(equality_count == len(edges), "equality count")
    print(json.dumps({
        "status": "FINITE-EXACT PASS",
        "universe": "all vertices at levels 0 through 7; every pair and edge pair",
        "vertices": len(vertices), "edges": len(edges),
        "vertex_pairs": pair_count, "edge_pairs": edge_pairs,
        "incident_edge_pairs": incident_pairs,
        "minimum_cosh_distance": "21/16",
        "minimum_attained_exactly_on_tree_edges": True,
        "circle_radius_squared_over_parent_height_squared": "185/64",
        "negative_controls": ["Euclidean separation tends to zero",
                              "not quasi-isometric: root-split graph distance 2k, cosh H-distance 3/2"],
        "all_level_status": "proved in PROOF.md, not by finite enumeration"
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    run()
