#!/usr/bin/env python3
"""Independent finite audit of whole-side ED carriers; no geometry solver claim."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def check_involution(p):
    need(type(p) is tuple, "involution must be a tuple")
    n = len(p)
    need(all(type(x) is int and 0 <= x < n for x in p), "bad image")
    need(all(p[p[i]] == i for i in range(n)), "not an involution")


def involutions(n):
    def visit(todo, p):
        if not todo:
            yield tuple(p)
            return
        a, *rest = todo
        p[a] = a
        yield from visit(rest, p)
        for b in rest:
            p[a], p[b] = b, a
            yield from visit([x for x in rest if x != b], p)
    yield from visit(list(range(n)), list(range(n)))


def incidence(a, b):
    check_involution(a)
    check_involution(b)
    need(len(a) == len(b), "different named side universes")
    nodes = {}
    edges = []
    for s in range(len(a)):
        endpoints = tuple((color, tuple(sorted({s, p[s]})))
                          for color, p in enumerate((a, b)))
        edges.append(endpoints)
        for node in endpoints:
            nodes.setdefault(node, []).append(s)
    return nodes, edges


def decode(nodes, edges):
    n = len(edges)
    decoded = [list(range(n)), list(range(n))]
    seen = [set(), set()]
    for node, links in nodes.items():
        color, block = node
        need(color in (0, 1), "bad shore")
        need(len(links) in (1, 2), "isolated or overfull ED node")
        need(len(set(links)) == len(links), "duplicate incidence")
        need(tuple(sorted(links)) == block, "block/link mismatch")
        need(all(type(s) is int and 0 <= s < n for s in links), "bad side")
        for s in links:
            need(s not in seen[color], "side has two endpoints on one shore")
            seen[color].add(s)
            need(edges[s][color] == node, "wrong incidence")
        if len(links) == 2:
            x, y = links
            decoded[color][x], decoded[color][y] = y, x
    need(all(z == set(range(n)) for z in seen), "missing endpoint")
    return tuple(map(tuple, decoded))


def side_components(a, b):
    remaining = set(range(len(a)))
    result = []
    while remaining:
        start = min(remaining)
        todo, visited = [start], set()
        while todo:
            x = todo.pop()
            if x in visited:
                continue
            visited.add(x)
            todo += [a[x], b[x]]
        remaining -= visited
        result.append(tuple(sorted(visited)))
    return sorted(result)


def graph_components(nodes, edges):
    remaining = set(nodes)
    result, shapes = [], []
    while remaining:
        todo, visited, links = [min(remaining)], set(), set()
        while todo:
            node = todo.pop()
            if node in visited:
                continue
            visited.add(node)
            for side in nodes[node]:
                links.add(side)
                todo += list(edges[side])
        remaining -= visited
        degrees = [len(nodes[x]) for x in visited]
        leaves = degrees.count(1)
        need(leaves in (0, 2), "neither path nor cycle")
        need(len(links) == len(visited) - (leaves == 2), "wrong Euler count")
        if not leaves:
            need(len(links) % 2 == 0, "odd bipartite cycle")
        result.append(tuple(sorted(links)))
        shapes.append((len(links), "path" if leaves else "cycle"))
    return sorted(result), sorted(shapes)


def overlap_area(p, q):
    return max(0, min(p[0] + 1, q[0] + 1) - max(p[0], q[0])) * \
        max(0, min(p[1] + 1, q[1] + 1) - max(p[1], q[1]))


def geometry_checks():
    origins = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
    deltas = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    outgoing = ["E", "N", "W", "S"]
    opposite = {"E": "W", "W": "E", "N": "S", "S": "N"}
    used = [[] for _ in origins]
    for i, (dx, dy) in enumerate(deltas):
        need(origins[i + 1] == (origins[i][0] + dx, origins[i][1] + dy),
             "incorrect translation transport")
        need(overlap_area(origins[i], origins[i + 1]) == 0, "adjacent overlap")
        used[i].append(outgoing[i])
        used[i + 1].append(opposite[outgoing[i]])
    need(all(len(set(x)) == len(x) for x in used), "port reused on a piece")
    overlaps = [(i, j, overlap_area(origins[i], origins[j]))
                for i in range(5) for j in range(i + 1, 5)
                if overlap_area(origins[i], origins[j]) > 0]
    need(overlaps == [(0, 4, 1)], "overlap witness lost")
    need(sum(x for x, y in deltas) == sum(y for x, y in deltas) == 0,
         "return displacement changed")
    # Reflection of the right square: its intrinsic E side maps to the shared
    # physical E side of the left square with endpoint indices unchanged.
    east = [(1, 0), (1, 1)]
    reflected = [(2 - x, y) for x, y in east]
    need(east == reflected, "reflection endpoint hostile changed")
    return {"five_square_origins": origins, "ports": used,
            "positive_area_overlaps": overlaps,
            "holonomy_scope": "tree consistency; not nonoverlap"}


def endpoint_parity_checks():
    tested = 0
    for steps in range(9):
        for signs in itertools.product((0, 1), repeat=steps + 1):
            for start in (0, 1):
                value = start
                for i in range(steps):
                    value ^= 1 ^ signs[i] ^ signs[i + 1]
                need(value == (start ^ (steps % 2) ^ signs[0] ^ signs[-1]),
                     "endpoint-sign telescoping failed")
                tested += 1
    # With flips between assemblies, an interior side can have different
    # orientation signs in its two occurrences. Do not telescope those away.
    for steps in range(7):
        for contact_signs in itertools.product((0, 1), repeat=2 * steps):
            value, parity = 0, steps % 2
            for i in range(steps):
                left, right = contact_signs[2 * i:2 * i + 2]
                value ^= 1 ^ left ^ right
                parity ^= left ^ right
            need(value == parity, "assembly-dependent endpoint parity failed")
            tested += 1
    need((1 ^ 0 ^ 0) ^ (1 ^ 1 ^ 0) == 1,
         "even-step mixed-orientation hostile lost")
    return tested


def expect_reject(fn):
    try:
        fn()
    except ValueError:
        return
    raise RuntimeError("hostile input accepted")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutate", choices=["involution", "overlap", "product"])
    args = parser.parse_args()
    census = []
    for n in range(8):
        bank = list(involutions(n))
        for a, b in itertools.product(bank, repeat=2):
            nodes, edges = incidence(a, b)
            need(decode(nodes, edges) == (a, b), "roundtrip failed")
            components, shapes = graph_components(nodes, edges)
            need(components == side_components(a, b), "component map failed")
            # Reverse-label equivariance of the named incidence carrier.
            rev_a = tuple(n - 1 - a[n - 1 - i] for i in range(n))
            rev_b = tuple(n - 1 - b[n - 1 - i] for i in range(n))
            rn, re = incidence(rev_a, rev_b)
            rc, rs = graph_components(rn, re)
            need(rs == shapes, "relabeling changed graph shapes")
        census.append({"sides": n, "involutions": len(bank), "pairs": len(bank) ** 2})
    id2, swap = (0, 1), (1, 0)
    product1 = tuple(id2[id2[x]] for x in range(2))
    product2 = tuple(swap[swap[x]] for x in range(2))
    need(product1 == product2 == id2, "product hostile changed")
    path_shapes = graph_components(*incidence(id2, id2))[1]
    cycle_shapes = graph_components(*incidence(swap, swap))[1]
    need(path_shapes == [(1, "path"), (1, "path")], "wrong singleton model")
    need(cycle_shapes == [(2, "cycle")], "parallel-edge cycle lost")
    # Flat T-junction: the long side and two short sides are all singletons
    # in paper ED; atomization splits long into two matching occurrences.
    paper = (0, 1, 2)
    atomized = (2, 3, 0, 1)
    need(graph_components(*incidence(paper, paper))[1] == [(1, "path")] * 3,
         "flat singleton distinction lost")
    need(graph_components(*incidence(atomized, atomized))[1] == [(2, "cycle")] * 2,
         "atomized seam carrier distinction lost")
    for bad in ((1, 2, 0), (True,), (1,), [0], (0.0,)):
        expect_reject(lambda bad=bad: check_involution(bad))
    nodes, edges = incidence(id2, id2)
    isolated = dict(nodes)
    isolated[(0, ())] = []
    expect_reject(lambda: decode(isolated, edges))
    duplicate = dict(nodes)
    duplicate[(0, (0,))] = [0, 0]
    expect_reject(lambda: decode(duplicate, edges))
    geo = geometry_checks()
    parity_count = endpoint_parity_checks()
    if args.mutate == "involution":
        check_involution((1, 2, 0))
    if args.mutate == "overlap":
        need(not geo["positive_area_overlaps"], "CORRUPTION: overlaps erased")
    if args.mutate == "product":
        need(path_shapes == cycle_shapes, "CORRUPTION: product claimed faithful")
    result = {"status": "PASS", "universe": "all ordered pairs of totalized partial involutions on n named sides, 0<=n<=7",
              "census": census, "pairs_total": sum(x["pairs"] for x in census),
              "active_checks": CHECKS, "endpoint_parity_cases": parity_count,
              "product_only_collision": {"product": product1,
                 "first_shapes": path_shapes, "second_shapes": cycle_shapes},
              "flat_control": "whole-side ED != atomized seam incidence",
              "geometry": geo,
              "scope": "finite carrier audit and exact hostile witnesses, not a dissection enumeration or Penrose extension solver"}
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
