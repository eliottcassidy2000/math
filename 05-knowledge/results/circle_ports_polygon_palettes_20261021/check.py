#!/usr/bin/env python3
"""Exact, standalone controls. No assertions disappear under python -O."""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def mod1(x):
    return x % 1


def phase_graph(phases, alpha=F(1, 6)):
    need(len(set(phases)) == len(phases), "Physical points must be distinct")
    adj = [set() for _ in phases]
    for i, j in combinations(range(len(phases)), 2):
        if mod1(phases[i] - phases[j]) in {alpha, 1 - alpha}:
            adj[i].add(j)
            adj[j].add(i)
    return adj


def components(adj):
    unseen = set(range(len(adj)))
    out = []
    while unseen:
        start = min(unseen)
        seen, stack = {start}, [start]
        while stack:
            u = stack.pop()
            for v in adj[u] - seen:
                seen.add(v)
                stack.append(v)
        unseen -= seen
        out.append(sorted(seen))
    return out


def is_colouring(adj, lists, colour):
    return (colour is not None and len(colour) == len(adj)
            and all(colour[v] in lists[v] for v in range(len(adj)))
            and all(colour[u] != colour[v]
                    for u in range(len(adj)) for v in adj[u]))


def exact_list_colouring(adj, lists):
    """Independent MRV backtracker, including empty and singleton lists."""
    colours = [None] * len(adj)

    def rec():
        free = [v for v in range(len(adj)) if colours[v] is None]
        if not free:
            return colours.copy()
        avail = {v: set(lists[v]) - {colours[u] for u in adj[v]} for v in free}
        v = min(free, key=lambda v: (len(avail[v]), -len(adj[v]), v))
        for c in sorted(avail[v]):
            colours[v] = c
            found = rec()
            if found is not None:
                return found
        colours[v] = None
        return None

    return rec()


def two_list_colouring(adj, lists):
    """Constructive path/even-cycle proof, not a general graph solver."""
    need(all(len(L) >= 2 for L in lists), "Need two colours in every list")
    need(all(len(a) <= 2 for a in adj), "Need maximum degree two")
    out = [None] * len(adj)
    for comp in components(adj):
        endpoints = [v for v in comp if len(adj[v]) < 2]
        if endpoints:
            order, prev, curr = [], None, min(endpoints)
            while curr is not None:
                order.append(curr)
                nxt = adj[curr] - ({prev} if prev is not None else set())
                prev, curr = curr, min(nxt) if nxt else None
            need(set(order) == set(comp), "Path traversal failed")
            for v in order:
                available = set(lists[v]) - {out[u] for u in adj[v]}
                need(bool(available), "Greedy path failed")
                out[v] = min(available)
            continue

        need(len(comp) % 2 == 0, "Odd cycle is not 2-choosable")
        start = min(comp)
        order = [start]
        prev, curr = start, min(adj[start])
        while curr != start:
            order.append(curr)
            nxt = (adj[curr] - {prev}).pop()
            prev, curr = curr, nxt
        need(set(order) == set(comp), "Cycle traversal failed")
        if all(set(lists[v]) == set(lists[start]) for v in comp):
            a, b = sorted(lists[start])[:2]
            for i, v in enumerate(order):
                out[v] = a if i % 2 == 0 else b
            continue

        oriented = None
        for u in comp:
            for v in sorted(adj[u]):
                difference = set(lists[u]) - set(lists[v])
                if difference:
                    oriented = u, v, min(difference)
                    break
            if oriented is not None:
                break
        need(oriented is not None, "Unequal cycle lists need a directed difference")
        u, final, c = oriented
        out[u] = c
        prev, curr = u, (adj[u] - {final}).pop()
        while curr != u:
            available = set(lists[curr]) - {out[w] for w in adj[curr]}
            need(bool(available), "Constructive cycle proof failed")
            out[curr] = min(available)
            prev, curr = curr, (adj[curr] - {prev}).pop()
    need(is_colouring(adj, lists, out), "Invalid constructive witness")
    return out


def polygon_controls():
    tested = 0
    for n in range(2, 11):
        A = [F(i, n) for i in range(n)]
        for m in range(2, 11):
            for shift in range(12):
                delta = F(shift, 12)
                B = [mod1(delta + F(j, m)) for j in range(m)]
                counts = {F(0): 0, F(1, 6): 0, F(5, 6): 0}
                for a, b in product(A, B):
                    d = mod1(b-a)
                    if d in counts:
                        counts[d] += 1
                L, g = lcm(n, m), gcd(n, m)
                for angle, count in counts.items():
                    expected = g if (L * (delta-angle)).denominator == 1 else 0
                    need(count == expected, "Polygon fibre-count disagreement")
                # Merge literal coincident points before constructing the graph.
                adj = phase_graph(sorted(set(A+B)))
                for comp in components(adj):
                    need(max(len(adj[v]) for v in comp) <= 2, "Bad degree")
                    if all(len(adj[v]) == 2 for v in comp):
                        need(len(comp) == 6, "Unit-radius cycle must have size six")
                tested += 1
    return tested


def rational_circle_controls():
    tested = 0
    for q in range(2, 13):
        for p in range(1, q//2 + 1):
            if gcd(p, q) != 1:
                continue
            alpha = F(p, q)
            for n in range(2, 37):
                adj = phase_graph([F(i, n) for i in range(n)], alpha)
                sizes = sorted(map(len, components(adj)))
                expect = [q]*(n//q) if n % q == 0 else [1]*n
                need(sizes == expect, "Single polygon classification failed")
                tested += 1
    return tested


def list_controls():
    # Every induced subset of C6, with every exact two-list from three colours.
    cases = 0
    pairs3 = list(combinations(range(3), 2))
    for mask in range(64):
        phases = [F(i, 6) for i in range(6) if mask >> i & 1]
        adj = phase_graph(phases)
        for lists in product(pairs3, repeat=len(phases)):
            got = two_list_colouring(adj, lists)
            other = exact_list_colouring(adj, lists)
            need(is_colouring(adj, lists, got), "Invalid constructive colouring")
            need(is_colouring(adj, lists, other), "Independent list solver disagrees")
            cases += 1
    need(cases == 4**6, "Wrong induced-subset list universe")

    # Tight five-colour, fixed-centre lists have four available colour names.
    cycle = phase_graph([F(i, 6) for i in range(6)])
    pairs4 = list(combinations(range(4), 2))
    cycle_cases = 0
    for lists in product(pairs4, repeat=6):
        got = two_list_colouring(cycle, lists)
        need(is_colouring(cycle, lists, got), "Four-colour list witness invalid")
        cycle_cases += 1
    need(cycle_cases == 6**6, "Wrong full-cycle list universe")

    # The odd-circle hostile: alpha=1/5 represents a different radius.
    odd = phase_graph([F(i, 5) for i in range(5)], F(1, 5))
    need(exact_list_colouring(odd, [(0, 1)]*5) is None, "Odd-cycle hostile failed")
    need(exact_list_colouring(odd, [(0, 1, 2)]*5) is not None, "C5 needs only three")
    # A proper precolouring can fail when its rim ports are not free.
    need(sum(map(len, phase_graph([F(i, 5) for i in range(5)]))) == 0,
         "Radius-one regular pentagon must have no rim edges")
    star = [{1, 2, 3, 4, 5}, {0}, {0}, {0}, {0}, {0}]
    need(exact_list_colouring(star, [tuple(range(5))] + [(i,) for i in range(5)])
         is None, "Fixed-rim-port hostile failed")
    return cases, cycle_cases


def dist2(a, b):
    # (x,y) represents the Euclidean point (x, sqrt(3)*y), exactly.
    return (a[0]-b[0])**2 + 3*(a[1]-b[1])**2


def sharp_geometry():
    pts = [(F(0), F(0)), (F(1), F(0)), (F(1, 2), F(1, 2))]
    edges = {(0, 1), (0, 2), (1, 2)}
    leaves, parameters = [], []
    for parent in (1, 2):
        for colour in range(3):
            t = 1
            while True:
                denominator = 1 + 3*t*t
                r = (F(1-3*t*t, denominator), F(2*t, denominator))
                candidate = (pts[parent][0]+r[0], pts[parent][1]+r[1])
                if (candidate not in pts and dist2(candidate, pts[parent]) == 1
                        and all(dist2(candidate, old) != 1
                                for j, old in enumerate(pts) if j != parent)):
                    break
                t += 1
                need(t < 1000, "Unexpected finite-avoidance failure")
            idx = len(pts)
            pts.append(candidate)
            edges.add(tuple(sorted((parent, idx))))
            leaves.append((idx, colour))
            parameters.append(t)
    actual = {ij for ij in combinations(range(len(pts)), 2)
              if dist2(pts[ij[0]], pts[ij[1]]) == 1}
    need(actual == edges, "Extra or missing unit edge in sharpness example")
    adj = [set() for _ in pts]
    for i, j in actual:
        adj[i].add(j)
        adj[j].add(i)
    fixed = {i: c for i, c in leaves}
    for centre_colour, expected in ((3, False), (0, True)):
        lists = [tuple(range(5)) for _ in pts]
        lists[0] = (centre_colour,)
        for i, c in fixed.items():
            lists[i] = (c,)
        found = exact_list_colouring(adj, lists)
        need((found is not None) == expected, "Sharp chosen-centre claim failed")
        if expected:
            repair = found
    three = exact_list_colouring(adj, [tuple(range(3)) for _ in pts])
    need(is_colouring(adj, [tuple(range(3)) for _ in pts], three), "Not 3-colourable")
    return {"points_x_y_sqrt3": [[str(a), str(b)] for a, b in pts],
            "unit_edges": sorted(map(list, actual)), "rotation_parameters": parameters,
            "off_circle_fixed_colours": leaves, "repaired_five_colouring": repair,
            "ordinary_three_colouring": three,
            "all_unordered_distances_checked": len(pts)*(len(pts)-1)//2}


def clock_controls():
    cases = 0
    speeds = range(13)
    for q in range(1, 37):
        for a in range(1, q+1):
            if gcd(a, q) != 1:
                continue
            for u, v in combinations(speeds, 2):
                angular = mod1(F(a*(u-v), q)) in {F(1, 6), F(5, 6)}
                congruence = q % 6 == 0 and (a*(u-v)) % q in {q//6, 5*q//6}
                need(angular == congruence, "Rational runner-clock criterion failed")
                cases += 1
    ap14 = phase_graph([F(i, 14) for i in range(14)])
    ap12 = phase_graph([F(i, 12) for i in range(12)])
    scaled14 = phase_graph([F(i, 14) for i in range(14)], F(1, 14))
    need(sum(map(len, ap14)) == 0, "AP14 radius-one hostile failed")
    need(sorted(map(len, components(ap12))) == [6, 6], "AP12 positive control failed")
    need(sorted(map(len, components(scaled14))) == [14], "Rescaled AP14 must be C14")
    return {"pair_clock_cases": cases, "AP14_radius_one_edges": 0,
            "AP12_component_sizes": [6, 6], "AP14_side_unit_component_sizes": [14]}


def main():
    subset_cases, full_cases = list_controls()
    receipt = {
        "status": "PASS",
        "scope": "exact finite controls of stated all-size proofs; not a graph census",
        "phased_polygon_cases": polygon_controls(),
        "rational_single_polygon_cases": rational_circle_controls(),
        "induced_C6_two_lists_independently_checked": subset_cases,
        "full_C6_four_colour_two_lists_with_witnesses": full_cases,
        "sharp_fixed_centre_example": sharp_geometry(),
        "runner_clock_controls": clock_controls(),
        "hostiles": ["odd cycle at a different radius is not bipartite",
                     "fixed rim ports cannot be treated as free two-lists",
                     "conditional failure does not imply high chromatic number"]}
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
