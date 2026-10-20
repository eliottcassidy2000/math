"""Independent C2 kernel audit using adjacency matrices and direct permutations.

Universe: every ordered deficiency-three primitive symbolic record at residual
budget j=2.  The six integer menus satisfy 2q+delta+u=3; defect graphs are
simple, have r=j-q edges, and have no isolated endpoint.  Endpoint weights
are one, with the allowed single degree-one endpoint of weight two when q=1.
Exceptional gap symbols are 0,1,2,L=3; H=99 is a prohibitive heavy symbol.

This file imports no repository module.  Pair-freeness is checked directly
from equality of adjacency rows outside the pair.  Ownership is checked by
every permutation inside the intrinsic heavy gaps, with direct weighted
backward-cost calculation and numerical record reconstruction.  This is a
finite-exact audit of C2, not a proof of the all-budget kernel formula.

Run from the repository root with Python, or Python -O, followed by this path.
"""

from collections import Counter
from itertools import combinations, permutations, product


H = 99


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def primitive(defects, weights, gaps):
    """Close marked endpoints under actual defect adjacency."""
    endpoint_count = len(weights)
    neighbors = [set() for _ in weights]
    for a, b in defects:
        neighbors[a].add(b)
        neighbors[b].add(a)
    active = {v for v in range(endpoint_count)
              if len(neighbors[v]) != 1 or weights[v] != 1}
    for i, weight in enumerate(gaps):
        if weight != H:
            if i:
                active.add(i - 1)
            if i < endpoint_count:
                active.add(i)
    while True:
        expanded = active.copy()
        for vertex in active:
            expanded.update(neighbors[vertex])
        if expanded == active:
            return len(active) == endpoint_count
        active = expanded


def has_homogeneous_pair(adjacency):
    size = len(adjacency)
    return any(all(adjacency[a][v] == adjacency[b][v]
                   for v in range(size) if v not in (a, b))
               for a, b in combinations(range(size), 2))


def accepted(defects, endpoint_weights, gaps, budget):
    weights, endpoints = [], []
    for i, weight in enumerate(gaps):
        if weight:
            weights.append(weight)
        if i < len(endpoint_weights):
            endpoints.append(len(weights))
            weights.append(endpoint_weights[i])
    size = len(weights)
    adjacency = [[int(a < b) for b in range(size)] for a in range(size)]
    for a, b in defects:
        x, y = endpoints[a], endpoints[b]
        adjacency[x][y], adjacency[y][x] = 0, 1
    if has_homogeneous_pair(adjacency):
        return False

    heavy = [v for v in range(size) if weights[v] == H]
    need(len(heavy) == 2 * budget - 2, "unexpected heavy count")
    light_gaps = [[] for _ in range(len(heavy) + 1)]
    for v in range(size):
        if weights[v] == H:
            continue
        signature = tuple(adjacency[z][v] for z in heavy)
        cut = sum(signature)
        need(signature == (1,) * cut + (0,) * (len(heavy) - cut),
             "light vertex lacks an intrinsic heavy cut")
        light_gaps[cut].append(v)

    source = (budget - len(defects),
              2 * len(defects) - len(endpoint_weights),
              sum(weight != H for weight in gaps),
              endpoint_weights, gaps, defects)
    source_seen = False
    minimum = None
    for parts in product(*(permutations(gap) for gap in light_gaps)):
        order = []
        for i, part in enumerate(parts):
            order.extend(part)
            if i < len(heavy):
                order.append(heavy[i])
        backward = [(a, b) for a, b in combinations(range(size), 2)
                    if adjacency[order[b]][order[a]]]
        cost = sum(weights[order[a]] * weights[order[b]]
                   for a, b in backward)
        if cost > budget:
            continue
        ordered_endpoints = sorted({v for edge in backward for v in edge})
        endpoint_index = {v: i for i, v in enumerate(ordered_endpoints)}
        gap_weights = []
        for a, b in zip([-1] + ordered_endpoints,
                        ordered_endpoints + [size]):
            need(b - a <= 2, "pair-free core has two unmarked gap vertices")
            gap_weights.append(weights[order[a + 1]] if b - a == 2 else 0)
        record = (
            budget - len(backward),
            2 * len(backward) - len(ordered_endpoints),
            sum(weight != H for weight in gap_weights),
            tuple(weights[order[v]] for v in ordered_endpoints),
            tuple(gap_weights),
            tuple((endpoint_index[a], endpoint_index[b]) for a, b in backward),
        )
        need(2 * record[0] + record[1] + record[2] == 3,
             "reconstructed record violates deficiency identity")
        source_seen = source_seen or record == source
        if minimum is None or record < minimum:
            minimum = record
    need(source_seen, "displayed source order absent from good-order search")
    return minimum == source


def enumerate_c2():
    budget = 2
    total, polynomial = 0, Counter()
    for q in (0, 1):
        edge_count = budget - q
        for delta in range(4 - 2 * q):
            exceptional_count = 3 - 2 * q - delta
            endpoint_count = 2 * edge_count - delta
            if endpoint_count < 2:
                continue
            possible_edges = tuple(combinations(range(endpoint_count), 2))
            for defects in combinations(possible_edges, edge_count):
                degrees = [sum(v in edge for edge in defects)
                           for v in range(endpoint_count)]
                if 0 in degrees:
                    continue
                assignments = [(1,) * endpoint_count]
                if q:
                    assignments += [tuple(1 + int(v == a)
                                          for v in range(endpoint_count))
                                    for a in range(endpoint_count)
                                    if degrees[a] == 1]
                for endpoint_weights in assignments:
                    need(sum(endpoint_weights[a] * endpoint_weights[b]
                             for a, b in defects) <= budget,
                         "source weighted cost exceeds budget")
                    for positions in combinations(range(endpoint_count + 1),
                                                  exceptional_count):
                        for values in product(range(4), repeat=exceptional_count):
                            gaps = [H] * (endpoint_count + 1)
                            for i, weight in zip(positions, values):
                                gaps[i] = weight
                            gaps = tuple(gaps)
                            if not primitive(defects, endpoint_weights, gaps):
                                continue
                            total += 1
                            if accepted(defects, endpoint_weights, gaps, budget):
                                polynomial[gaps.count(3)] += 1
    return total, polynomial


def main():
    triangle = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    need(not has_homogeneous_pair(triangle), "directed triangle positive control")
    need(has_homogeneous_pair([[0, 1], [0, 0]]), "transitive pair hostile control")
    total, polynomial = enumerate_c2()
    expected = {0: 756, 1: 792, 2: 270, 3: 30}
    need(total == 2244, "unexpected primitive C2 universe")
    need(dict(polynomial) == expected, "independent C2 coefficients disagree")
    print("Universe: all primitive symbolic D3 records at residual budget j=2.")
    print("Independent adjacency-pair tests and exhaustive intrinsic-gap permutations.")
    print("Primitive records:", total)
    print("C2 coefficient list:", [polynomial[i] for i in range(4)])
    print("PASS: C2(y)=756+792y+270y^2+30y^3; no repository imports or assertions.")


if __name__ == "__main__":
    main()
