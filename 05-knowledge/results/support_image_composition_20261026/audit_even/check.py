#!/usr/bin/env python3
"""Independent root-first audit. No producer or repository imports."""
import itertools as it
import json
import math
from fractions import Fraction
from functools import lru_cache

checks = 0


def check(test, explanation):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(explanation)


@lru_cache(None)
def pairs(n):
    return tuple(it.combinations(range(n), 2))


@lru_cache(None)
def permutations(n):
    return tuple(it.permutations(range(n)))


def bits(n, mask):
    return frozenset(e for k, e in enumerate(pairs(n)) if mask >> k & 1)


def encode(n, edges):
    return sum((e in edges) << k for k, e in enumerate(pairs(n)))


def image_graph(n, edges, p):
    return encode(n, {tuple(sorted((p[u], p[v]))) for u, v in edges})


def sign(edges, p):
    return (-1) ** sum(p[u] > p[v] for u, v in edges)


def sign_product(edges, p):
    # Primary-source product of integer differences, computed literally over Q.
    result = Fraction(1)
    for u, v in edges:
        result *= Fraction(p[v] - p[u], v - u)
    return result


@lru_cache(None)
def automorphisms(n, mask):
    edges = bits(n, mask)
    return tuple(p for p in permutations(n) if image_graph(n, edges, p) == mask)


@lru_cache(None)
def is_even(n, mask):
    edges = bits(n, mask)
    return all(sign(edges, p) == 1 for p in automorphisms(n, mask))


def graph_representatives(n):
    # Remove whole relabelling orbits from the complete labelled-mask universe.
    remaining = set(range(1 << len(pairs(n))))
    reps = []
    while remaining:
        mask = min(remaining)
        edges = bits(n, mask)
        orbit = {image_graph(n, edges, p) for p in permutations(n)}
        remaining.difference_update(orbit)
        reps.append(mask)
    return reps


def vertex_orbits(n, group):
    remaining = set(range(n))
    result = []
    while remaining:
        v = min(remaining)
        orbit = {p[v] for p in group}
        remaining.difference_update(orbit)
        result.append(orbit)
    return result


def induced(n, edges, vertices):
    inverse = {v: i for i, v in enumerate(vertices)}
    return encode(len(vertices), {tuple(sorted((inverse[u], inverse[v])))
                                  for u, v in edges if u in inverse and v in inverse})


def components(n, edges):
    remaining = set(range(n))
    result = []
    while remaining:
        comp = {min(remaining)}
        while True:
            expanded = comp | {v for u, v in edges if u in comp} | {u for u, v in edges if v in comp}
            if expanded == comp:
                break
            comp = expanded
        result.append(sorted(comp))
        remaining.difference_update(comp)
    return result


def marked_parent_key(n, edges, root):
    vertices = [v for v in range(n) if v != root]
    g = bits(n - 1, induced(n, edges, vertices))
    neighbours = {j for j, v in enumerate(vertices)
                  if tuple(sorted((root, v))) in edges}
    return min((image_graph(n - 1, g, p), sum(1 << p[v] for v in neighbours))
               for p in permutations(n - 1))


def card_transports(n, edges, root, target):
    # Construct each order-(n-1) card map on original vertex names, then check
    # graph incidences and the deleted vertex's neighbourhood separately.
    domain = [v for v in range(n) if v != root]
    codomain = [v for v in range(n) if v != target]
    ans = []
    for destinations in it.permutations(codomain):
        f = dict(zip(domain, destinations))
        if any(((u, v) in edges) != (tuple(sorted((f[u], f[v]))) in edges)
               for u, v in it.combinations(domain, 2)):
            continue
        if any((tuple(sorted((root, u))) in edges) !=
               (tuple(sorted((target, f[u]))) in edges) for u in domain):
            continue
        f[root] = target
        ans.append(tuple(f[v] for v in range(n)))
    return tuple(ans)


def tournament_census(n):
    remaining = set(range(1 << len(pairs(n))))
    count = rooted = 0
    while remaining:
        mask = min(remaining)
        arcs = [(u, v) if mask >> k & 1 else (v, u)
                for k, (u, v) in enumerate(pairs(n))]
        images, group = set(), []
        for p in permutations(n):
            image = encode(n, {(p[u], p[v]) for u, v in arcs if p[u] < p[v]})
            images.add(image)
            if image == mask:
                group.append(p)
        remaining.difference_update(images)
        count += 1
        rooted += len(vertex_orbits(n, group))
    return count, rooted


def main():
    graphs = [graph_representatives(n) for n in range(7)]
    check([len(x) for x in graphs] == [1, 1, 2, 4, 11, 34, 156], 'graph universe')
    even = [0] * 7
    connected_defect = [0] * 7
    defect = [0] * 7
    rooted_even = [0] * 7
    recurrence = [Fraction(1)] + [Fraction(0)] * 6
    even_parent_only = [Fraction(1)] + [Fraction(0)] * 6
    rooted_rows = [0] * 7
    transport_maps = 0
    for n, reps in enumerate(graphs):
        parent_keys = set()
        for mask in reps:
            edges = bits(n, mask)
            group = automorphisms(n, mask)
            char = {p: sign(edges, p) for p in group}
            for p in group:
                check(sign_product(edges, p) == char[p], 'edge-character product')
            glob_even = all(c == 1 for c in char.values())
            even[n] += glob_even
            check(sum(char.values()) == (len(group) if glob_even else 0), 'global cancellation')
            comps = components(n, edges)
            comp_info = [(c, induced(n, edges, c)) for c in comps]
            comp_even = [is_even(len(c), cmask) for c, cmask in comp_info]
            check(glob_even == all(comp_even), 'component-even equivalence')
            for orbit in vertex_orbits(n, group):
                root = min(orbit)
                rooted_rows[n] += 1
                rooted_even[n] += glob_even
                key = marked_parent_key(n, edges, root)
                check(key not in parent_keys, 'two rooted children share one marked parent')
                parent_keys.add(key)
                stabilizer = tuple(p for p in group if p[root] == root)
                rel_even = all(char[p] == 1 for p in stabilizer)
                charge = Fraction(0)
                fibre_signs = {}
                for target in range(n):
                    cards = card_transports(n, edges, root, target)
                    transport_maps += len(cards)
                    expected = {p for p in group if p[root] == target}
                    check(set(cards) == expected, 'marked card / full-child map mismatch')
                    check(len(cards) in (0, len(stabilizer)), 'coset fibre size')
                    fibre_signs[target] = {char[p] for p in cards}
                    charge += Fraction(sum(char[p] for p in cards), len(stabilizer))
                    if rel_even:
                        check(len(fibre_signs[target]) == (target in orbit), 'transporter choice dependence')
                    else:
                        check(sum(char[p] for p in cards) == 0, 'nonrelative coset failed to cancel')
                check(charge == (len(orbit) if glob_even else 0), 'individual root-card charge')
                recurrence[n] += charge / n
                if is_even(n - 1, key[0]):
                    even_parent_only[n] += charge / n

                # Independent component criterion using the component's own
                # automorphisms, not just restrictions of the full group.
                ci = next(j for j, c in enumerate(comps) if root in c)
                c, cmask = comp_info[ci]
                local_root = c.index(root)
                local_edges = bits(len(c), cmask)
                local_group = automorphisms(len(c), cmask)
                local_rel = all(sign(local_edges, p) == 1 for p in local_group if p[local_root] == local_root)
                criterion = local_rel and all(ok for j, ok in enumerate(comp_even) if j != ci)
                check(rel_even == criterion, 'E3 rooted component criterion')
                if rel_even:
                    eta = {}
                    for u in orbit:
                        for w in orbit:
                            vals = {char[p] for p in group if p[u] == w}
                            check(len(vals) == 1, 'orbit sign not well-defined')
                            eta[u, w] = next(iter(vals))
                    for u, w, t in it.product(orbit, repeat=3):
                        check(eta[u, w] * eta[w, t] == eta[u, t], 'transport composition')
                    for u, w in it.product(orbit, repeat=2):
                        check(eta[u, w] == eta[root, u] * eta[root, w], 'rank-one identity')
                    check(sum(eta[root, w] for w in orbit) == (len(orbit) if glob_even else 0),
                          'balanced root signs')
                    if not glob_even:
                        defect[n] += 1
                        connected_defect[n] += len(comps) == 1
                        check(sum(not ok for ok in comp_even) == 1 and not comp_even[ci],
                              'unique odd root component')
                        local_orbit = {c[p[local_root]] for p in local_group}
                        check(local_orbit == orbit, 'root orbit changed under even decorations')
                        inverse = {v: i for i, v in enumerate(c)}
                        for p in group:
                            restriction = tuple(inverse[p[v]] for v in c)
                            check(char[p] == sign(local_edges, restriction), 'odd-component sign localization')

        if n:
            # Check completeness against ALL marked parent orbits, separately
            # from the root-first card enumeration used above.
            attachment_keys = set()
            for parent in graphs[n - 1]:
                pedges = bits(n - 1, parent)
                for s in range(1 << (n - 1)):
                    attachment_keys.add(min((image_graph(n - 1, pedges, p),
                                             sum(1 << p[v] for v in range(n - 1) if s >> v & 1))
                                            for p in permutations(n - 1)))
            check(parent_keys == attachment_keys, 'rooted child / attachment orbit bijection')
            check(recurrence[n] == even[n], 'E1 recurrence')

    check(even == [1, 1, 1, 2, 4, 12, 56], 'even-graph sequence')
    check(connected_defect == [0, 0, 1, 0, 4, 9, 58], 'connected defect sequence')
    check(defect == [0, 0, 1, 1, 5, 15, 75], 'defect sequence')
    for n in range(7):
        check(defect[n] == sum(even[j] * connected_defect[n - j] for j in range(n + 1)),
              'E4 convolution')
    tournaments = [tournament_census(n) for n in range(7)]
    check([a for a, _ in tournaments] == even, 'direct tournament / even count equality')
    check([r for _, r in tournaments] == [a + b for a, b in zip(rooted_even, defect)],
          'relative-even pointed identity')
    check([r for _, r in tournaments] == [0, 1, 2, 4, 12, 48, 296], 'pointed tournament sequence')
    check(even_parent_only[3] == Fraction(4, 3), 'odd parent hostile')
    check(rooted_rows == [0, 1, 2, 6, 20, 90, 544], 'rooted attachment census')
    print(json.dumps({'status': 'PASS', 'checks': checks,
                      'method': 'full labelled graph orbit removal, root-first marked-card maps, independent tournament orbit removal',
                      'graph_types_by_order': [len(g) for g in graphs],
                      'rooted_attachment_types_by_order': rooted_rows,
                      'transport_maps_checked': transport_maps,
                      'even_graphs': even,
                      'connected_rooted_odd_defect': connected_defect,
                      'all_rooted_odd_defect': defect,
                      'pointed_tournaments': [r for _, r in tournaments],
                      'even_parent_only_order_three': str(even_parent_only[3])},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
