#!/usr/bin/env python3
"""Exact controls for the five-interior extension theorem; see FIVE_INTERIORS.md."""
from fractions import Fraction as F
from itertools import combinations, product
import json

CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def chromatic(adj):
    n = len(adj)
    order = sorted(range(n), key=lambda v: -len(adj[v]))
    colours = [-1]*n
    def visit(i, k):
        if i == n:
            return True
        v = order[i]; banned = {colours[w] for w in adj[v]}
        for c in range(k):
            if c not in banned:
                colours[v] = c
                if visit(i+1, k):
                    return True
        colours[v] = -1
        return False
    return next(k for k in range(1, n+1) if visit(0, k)) if n else 0


def complement(parts):
    n = sum(s for _, s in parts); H = [set() for _ in range(n)]
    off = 0
    for typ, s in parts:
        for i in range(s-1):
            H[off+i].add(off+i+1); H[off+i+1].add(off+i)
        if typ == 'C':
            H[off].add(off+s-1); H[off+s-1].add(off)
        off += s
    return [set(range(n))-{v}-H[v] for v in range(n)]


def delete(adj, v):
    keep = [w for w in range(len(adj)) if w != v]
    return [{j for j, x in enumerate(keep) if x in adj[w]} for w in keep]


def critical_controls():
    types = [('P', n) for n in range(1, 10)]+[('C', n) for n in range(3, 10)]
    def parts(n, start=0):
        if n == 0:
            yield ()
        for i in range(start, len(types)):
            t, s = types[i]
            if s <= n:
                for tail in parts(n-s, i):
                    yield ((t, s),)+tail
    tested = 0; critical = []
    for n in (7, 8, 9):
        for p in parts(n):
            G = complement(p)
            # Minimum degree six is necessary for seven-criticality.
            if min(map(len, G)) < 6:
                continue
            tested += 1
            if chromatic(G) == 7 and all(chromatic(delete(G, v)) == 6 for v in range(n)):
                critical.append(p)
    expected = [(('P', 1),)*7, (('P', 1),)*4+(('C', 5),)]
    need(critical == expected, ('critical quotient classification', critical))
    return dict(complement_component_types_tested=tested, seven_critical_complements=critical)


def sub(a, b):
    return a[0]-b[0], a[1]-b[1]


def norm(p):
    x, y = p
    return x*x+x*y+y*y


def centre(p, q, r):
    # Physical basis (1,0),(1/2,sqrt(3)/2), solved entirely over Q.
    a, b = sub(q, p), sub(r, p)
    a0, a1 = 2*a[0]+a[1], a[0]+2*a[1]
    b0, b1 = 2*b[0]+b[1], b[0]+2*b[1]
    det = a0*b1-a1*b0
    if not det:
        return None
    c, d = norm(q)-norm(p), norm(r)-norm(p)
    z = F(c*b1-a1*d, det), F(a0*d-c*b0, det)
    return z if norm(sub(z, p)) == 1 else None


def geometry_controls():
    fan = [(0, 0), (1, 0), (0, 1), (-1, 1), (-1, 0)]
    external = set()
    for triple in combinations(fan, 3):
        z = centre(*triple)
        if z is not None and z not in fan:
            external.add(z)
    need(not external, 'unit K1 join P4 has no external triple-contact centre')
    directions = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    cycles = 0; induced = 0
    for steps in product(directions, repeat=3):
        vs = [(0, 0), (1, 0)]
        for d in steps:
            vs.append((vs[-1][0]+d[0], vs[-1][1]+d[1]))
        if len(set(vs)) != 5 or norm(sub(vs[-1], vs[0])) != 1:
            continue
        cycles += 1
        diagonal = any(norm(sub(vs[j], vs[i])) == 1
                       for i, j in combinations(range(5), 2) if j-i not in (1, 4))
        induced += not diagonal
    need(induced == 0 and cycles > 0, 'triangular lattice has no induced five-cycle')
    # Positive control: a unit hexagon DOES have an external centre meeting
    # all six vertices. The pentagon argument must not extend by counting.
    need(centre(directions[0], directions[1], directions[2]) == (0, 0), 'hexagon hostile to a six-cycle extension')
    return dict(fan_external_triple_centres=0, rooted_directed_lattice_five_cycles=cycles,
                induced_lattice_five_cycles=induced, hexagon_common_centre=[0, 0])


def main():
    critical = critical_controls(); geometry = geometry_controls()
    return dict(status='PROVED in FIVE_INTERIORS.md; finite exact independent controls',
                theorem='Every proper six-colouring constant on four disjoint port pairs extends across at most five private interior vertices in the Euclidean unit-distance graph',
                critical=critical, geometry=geometry, checks=CHECKS,
                exclusions='No conclusion for six interiors, coincident ports, arbitrary abstract graphs, or global plane colouring')


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
