#!/usr/bin/env python3
"""Exact triangular-lattice port controls, with all equality partitions retained."""
from itertools import combinations
from collections import Counter
import json

CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def matchings(xs):
    if not xs:
        yield ()
        return
    a = xs[0]
    for i in range(1, len(xs)):
        b = xs[i]
        for rest in matchings(xs[1:i]+xs[i+1:]):
            yield ((a, b),)+rest


def partitions(n):
    def rec(xs):
        if len(xs) == n:
            yield tuple(xs)
        else:
            for c in range(max(xs, default=-1)+2):
                yield from rec(xs+[c])
    yield from rec([])


def unit(p, q):
    a, b = p[0]-q[0], p[1]-q[1]
    return a*a+a*b+b*b == 1


def extends(rows, boundary, inner):
    colours = dict(boundary)
    order = sorted(inner, key=lambda v: -len(rows[v]))
    def rec(i):
        if i == len(order):
            return True
        v = order[i]
        banned = {colours[w] for w in rows[v] if w in colours}
        for c in range(6):
            if c not in banned:
                colours[v] = c
                if rec(i+1):
                    return True
                del colours[v]
        return False
    return rec(0)


def main():
    Ms = list(matchings(tuple(range(8))))
    Ps = list(partitions(4))
    need(len(Ms) == 105 and len(Ps) == 15, 'universe counts')
    special = {M for M in Ms if len({a ^ b for a, b in M}) == 4}
    need(len(special) == 56, 'four-direction matching bridge')
    points = [(a, b) for a in range(4) for b in range(3)]
    rows = {i: {j for j in range(12) if unit(points[i], points[j])} for i in range(12)}
    count = Counter()
    for inner in combinations(range(12), 4):
        ports = [i for i in range(12) if i not in inner]
        for M in Ms:
            for P in Ps:
                boundary = {ports[v]: P[k] for k, pair in enumerate(M) for v in pair}
                count['candidate_patterns'] += 1
                if any(boundary[v] == boundary[w] for v in ports for w in rows[v] if w in boundary):
                    count['improper_boundary_rejected'] += 1
                    continue
                need(extends(rows, boundary, inner), ('unit selector violated lemma', inner, M, P))
                count['proper_patterns_extend'] += 1
                if M in special:
                    count['proper_patterns_in_56_family'] += 1
    # Exact geometric reason: three exterior double-contact points, no fourth.
    triangle = [(0, 0), (1, 0), (0, 1)]
    doubles = [p for p in ((a, b) for a in range(-3, 4) for b in range(-3, 4))
               if p not in triangle and sum(unit(p, q) for q in triangle) >= 2]
    need(sorted(doubles) == [(-1, 1), (1, -1), (1, 1)], 'triangle contact certificate')
    # Hostile abstract graph: K3 joined to eight independent ports. Four
    # distinct prescribed pair colours leave just two colours for its triangle.
    abstract = {v: set() for v in range(11)}
    for v in range(8, 11):
        abstract[v] = set(range(11))-{v}
    for v in range(8):
        abstract[v] = set(range(8, 11))
    need(not extends(abstract, {i: i//2 for i in range(8)}, range(8, 11)),
         'abstract three-interior hostile should fail')
    # Exhaust all complements on eight vertices whose matching number <=1:
    # they are stars or triangles. These are exactly the only candidates after
    # the two-disjoint-nonedges six-colouring test.
    shapes = {frozenset()}
    edges = list(combinations(range(8), 2))
    for centre in range(8):
        leaves = [v for v in range(8) if v != centre]
        for mask in range(1 << 7):
            shapes.add(frozenset(tuple(sorted((centre, leaves[k]))) for k in range(7) if mask >> k & 1))
    for tri in combinations(range(8), 3):
        shapes.add(frozenset(combinations(tri, 2)))
    for nonedges in shapes:
        triangle_missing = any(set(combinations(t, 2)) <= nonedges for t in combinations(range(8), 3))
        clique7 = any(not (set(combinations(t, 2)) & nonedges) for t in combinations(range(8), 7))
        need(triangle_missing or clique7, 'seven-chromatic small graph lemma')
    return dict(status='FINITE-EXACT', points=points, interior_choices=495,
                matchings=105, oriented_bridge_matchings=56, equality_partitions=15,
                **dict(count), exterior_triangle_double_contacts=doubles,
                abstract_three_interior_hostile_rejected=True,
                small_complement_controls=len(shapes), checks=CHECKS)


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
