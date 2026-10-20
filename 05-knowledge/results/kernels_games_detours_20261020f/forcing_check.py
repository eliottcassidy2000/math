"""Exact geometric and list-colouring checks; standard library only.

Q(sqrt(3),sqrt(11)) arithmetic uses the ordered basis 1,sqrt3,sqrt11,sqrt33.
All guards are explicit exceptions, including in python -O.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def number(a=0, b=0, c=0, d=0):
    return tuple(F(x) for x in (a, b, c, d))


ZERO, ONE = number(), number(1)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(a, q):
    return tuple(x * q for x in a)


def mul(a, b):
    out = [F(0)] * 4
    for i in range(4):
        for j in range(4):
            factor = (3 if i & j & 1 else 1) * (11 if i & j & 2 else 1)
            out[i ^ j] += a[i] * b[j] * factor
    return tuple(out)


def inverse(a):
    need(a != ZERO, "division by zero")
    conjugates = [tuple(x * (-1 if (i & mask).bit_count() % 2 else 1)
                        for i, x in enumerate(a)) for mask in (1, 2, 3)]
    numerator = mul(mul(conjugates[0], conjugates[1]), conjugates[2])
    norm = mul(a, numerator)
    need(norm[1:] == (0, 0, 0) and norm[0] != 0, "field norm")
    out = scale(numerator, 1 / norm[0])
    need(mul(a, out) == ONE, "inverse round trip")
    return out


def point_add(a, b):
    return tuple(add(x, y) for x, y in zip(a, b))


def point_sub(a, b):
    return tuple(sub(x, y) for x, y in zip(a, b))


def norm2(p):
    return add(mul(p[0], p[0]), mul(p[1], p[1]))


def distance2(p, q):
    return norm2(point_sub(p, q))


def adjacency(points):
    need(len(set(points)) == len(points), "coincident points")
    graph = [set() for _ in points]
    for i, j in combinations(range(len(points)), 2):
        if distance2(points[i], points[j]) == ONE:
            graph[i].add(j)
            graph[j].add(i)
    return graph


def edges(graph):
    return [(v, w) for v in range(len(graph)) for w in sorted(graph[v]) if v < w]


def list_colour(graph, lists):
    colours = {}

    def go():
        if len(colours) == len(graph):
            return dict(colours)
        choices = [(set(lists[v]) - {colours[w] for w in graph[v] if w in colours}, v)
                   for v in range(len(graph)) if v not in colours]
        allowed, v = min(choices, key=lambda item: (len(item[0]), -len(graph[item[1]]), item[1]))
        for c in sorted(allowed):
            colours[v] = c
            result = go()
            if result is not None:
                return result
        colours.pop(v, None)
        return None

    return go()


def check_colouring(graph, lists, colouring):
    need(colouring is not None, "missing colouring")
    need(all(colouring[v] in lists[v] for v in range(len(graph))), "list violation")
    need(all(colouring[v] != colouring[w] for v, w in edges(graph)), "edge violation")


def build_geometry():
    # Two unit diamonds about 0; the second is rotated by (5+i sqrt11)/6.
    p = [
        (ZERO, ZERO),
        (ONE, ZERO),
        (number(F(1, 2)), number(0, F(1, 2))),
        (number(F(3, 2)), number(0, F(1, 2))),
        (number(F(5, 6)), number(0, 0, F(1, 6))),
        (number(F(5, 12), 0, 0, F(-1, 12)), number(0, F(5, 12), F(1, 12))),
        (number(F(15, 12), 0, 0, F(-1, 12)), number(0, F(5, 12), F(3, 12))),
    ]
    core = adjacency(p)
    expected = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3),
                (2, 3), (3, 6), (4, 5), (4, 6), (5, 6)]
    need(edges(core) == expected, "exact spindle edge set")
    # Each triple (i,j,k) supplies the other common unit neighbour p_i+p_j-p_k.
    groups = [
        [(1, 3, 2), (2, 4, 0), (5, 6, 4)],
        [(1, 4, 0), (2, 6, 3), (3, 5, 6)],
        [(1, 5, 0), (2, 3, 1), (4, 6, 5)],
    ]
    port_colours = {}
    for c, group in enumerate(groups, 3):
        need(sorted(i for triple in group for i in triple[:2]) == list(range(1, 7)),
             "perfect matching cover")
        for i, j, k in group:
            need(k in core[i] & core[j], "reflection source")
            v = len(p)
            p.append(point_sub(point_add(p[i], p[j]), p[k]))
            port_colours[v] = c
    graph16 = adjacency(p)
    need([graph16[v] & set(range(7)) for v in range(7, 16)] ==
         [set(t[:2]) for group in groups for t in group], "port core contacts")
    need(all(port_colours[v] != port_colours[w]
             for v, w in edges(graph16) if v >= 7), "proper nine-port boundary")
    need([len(graph16[v]) for v in range(7)] == [4, 6, 6, 6, 6, 6, 6],
         "gate interior degrees")
    # Private central ports: deterministic rational unit directions, all contacts tested.
    tangents = []
    for c in range(3, 6):
        for t in range(1, 513):
            candidate = (number(F(1 - t * t, 1 + t * t)),
                         number(F(2 * t, 1 + t * t)))
            if candidate in p:
                continue
            if any(distance2(candidate, old) == ONE for old in p[1:]):
                continue
            p.append(candidate)
            port_colours[len(p) - 1] = c
            tangents.append(t)
            break
        else:
            raise RuntimeError("finite private-port construction exhausted")
    graph19 = adjacency(p)
    need(edges(graph19) == sorted(edges(graph16) + [(0, 16), (0, 17), (0, 18)]),
         "no accidental new unit contacts")
    need(all(port_colours[v] != port_colours[w]
             for v, w in edges(graph19) if v >= 7), "proper twelve-port boundary")
    need([len(graph19[v]) for v in range(7)] == [7, 6, 6, 6, 6, 6, 6],
         "one excess degree")
    return p, core, graph16, graph19, port_colours, groups, tangents


def check_circumcentres(core_points):
    histogram = Counter()
    unit_centres = []
    for triple in combinations(range(7), 3):
        p0, p1, p2 = (core_points[i] for i in triple)
        a, b = point_sub(p1, p0), point_sub(p2, p0)
        determinant = sub(mul(a[0], b[1]), mul(a[1], b[0]))
        if determinant == ZERO:
            histogram["collinear"] += 1
            continue
        r1 = scale(sub(norm2(p1), norm2(p0)), F(1, 2))
        r2 = scale(sub(norm2(p2), norm2(p0)), F(1, 2))
        inv = inverse(determinant)
        centre = (mul(sub(mul(r1, b[1]), mul(a[1], r2)), inv),
                  mul(sub(mul(a[0], r2), mul(r1, b[0])), inv))
        need(distance2(centre, p0) == distance2(centre, p1) == distance2(centre, p2),
             "circumcentre equal distances")
        if distance2(centre, p0) == ONE:
            need(centre in core_points, ("external triple neighbour", triple))
            histogram["unit centre is core vertex"] += 1
            unit_centres.append((triple, core_points.index(centre)))
        else:
            histogram["circumradius not one"] += 1
    need(sum(histogram.values()) == 35, "all triples covered")
    return dict(histogram), unit_centres


def diamond_colour(lists):
    # 0,1 are nonadjacent; 2,3 are adjacent and see both.
    for a in sorted(lists[0]):
        for b in sorted(lists[1]):
            for c in sorted(set(lists[2]) - {a, b}):
                for d in sorted(set(lists[3]) - {a, b, c}):
                    return {0: a, 1: b, 2: c, 3: d}
    return None


def list_controls():
    diamond = [{2, 3}, {2, 3}, {0, 1, 3}, {0, 1, 2}]
    total = 0
    for q in range(3, 7):
        pairs = list(combinations(range(q), 2))
        triples = list(combinations(range(q), 3))
        for a, b, c, d in product(pairs, pairs, triples, triples):
            lists = list(map(set, (a, b, c, d)))
            check_colouring(diamond, lists, diamond_colour(lists))
            total += 1
    need(total == 100585, "diamond universe")
    # All available-colour set pairs, including empty lists, check the two-hub formula.
    palette = set(range(6))
    subsets = [{c for c in palette if mask & (1 << c)} for mask in range(64)]
    for a, b in product(subsets, repeat=2):
        actual = sum(x != y for x in a for y in b)
        need(actual == len(a) * len(b) - len(a & b), "two-hub extension count")
    return total


def main():
    points, core, g16, g19, fixed, groups, tangents = build_geometry()
    hist, centres = check_circumcentres(points[:7])
    need(list_colour(core, [set(range(3)) for _ in core]) is None, "spindle not three-colourable")
    witness4 = list_colour(g19, [set(range(4)) for _ in g19])
    check_colouring(g19, [set(range(4)) for _ in g19], witness4)
    # Exact central-output relation of the 16-point gate under the displayed boundary.
    output_counts = []
    for c in range(6):
        count = 0
        for assignment in product(range(3), repeat=6):
            word = (c,) + assignment
            count += all(word[v] != word[w] for v, w in edges(core))
        output_counts.append(count)
    need(output_counts == [0, 0, 0, 24, 24, 24], "palette-complement gate")
    lists19 = [set(range(6)) - {fixed[w] for w in g19[v] if w in fixed} for v in range(7)]
    need(lists19 == [set(range(3))] * 7, "explicit blocked lists")
    need(list_colour(core, lists19) is None, "nineteen-point blocked boundary")
    # The 18-point central-port deletion has an explicit positive control.
    lists18 = [set(x) for x in lists19]
    lists18[0].add(5)
    witness18 = list_colour(core, lists18)
    check_colouring(core, lists18, witness18)
    # Every individual displayed port-colour constraint is necessary for this witness:
    # removing a private port frees one core colour, removing a paired one frees two.
    deletions = {}
    for port in range(7, 19):
        lists = [set(range(6)) - {fixed[w] for w in g19[v] if w in fixed and w != port}
                 for v in range(7)]
        colouring = list_colour(core, lists)
        check_colouring(core, lists, colouring)
        deletions[port] = colouring
    diamond_cases = list_controls()
    result = {
        "status": "PASS",
        "field_basis": ["1", "sqrt3", "sqrt11", "sqrt33"],
        "spindle_edges": edges(core),
        "reflection_groups_colours_3_4_5": groups,
        "private_central_tangents": tangents,
        "points": [[[str(c) for c in x] for x in p] for p in points],
        "graph19_edges": edges(g19),
        "boundary_colouring": fixed,
        "gate_central_extension_counts": output_counts,
        "global_four_colouring": witness4,
        "eighteen_point_core_colouring": witness18,
        "all_twelve_port_deletion_extensions": deletions,
        "circumcentre_histogram": hist,
        "unit_circumcentres": centres,
        "diamond_minimum_list_assignments": diamond_cases,
        "two_hub_available_list_pairs": 4096,
    }
    target = Path(__file__).with_suffix(".json")
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("Exact field: Q(sqrt3,sqrt11); distinct geometric points:", len(points))
    print("Spindle edges:", len(edges(core)), "; gate16 edges:", len(edges(g16)),
          "; full19 edges:", len(edges(g19)))
    print("Core triple circumcentres:", hist)
    print("Private central rational direction parameters:", tangents)
    print("Gate output multiplicities at colours 0..5:", output_counts)
    print("19-point boundary blocked; all 12 single-port deletions extend.")
    print("Whole 19-point graph is 4-colourable, and its spindle is not 3-colourable.")
    print("100585 degree-list diamond cases and 4096 two-hub list pairs pass.")
    print("PASS exact contacts, minimum-port certificate, lists and hostile controls.")


if __name__ == "__main__":
    main()
