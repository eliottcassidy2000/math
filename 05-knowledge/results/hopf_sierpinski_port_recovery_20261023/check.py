"""Exact checks for connected three-boundary affine-F3 section carriers.

Default replays a frozen JSON receipt. --write creates it. No external modules.
The observer is parallel-owner compatibility, not tournament reconstruction.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb, factorial
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
COUNT = 0


def require(condition, message):
    global COUNT
    COUNT += 1
    if not condition:
        raise RuntimeError(message)


F = range(3)
AFF = tuple(product((1, -1), F))


def compose(g, h):
    a, b = g
    c, d = h
    return (a * c, (a * d + b) % 3)


def inv(g):
    a, b = g
    return (a, (-a * b) % 3)


def act(g, x):
    return (g[0] * x + g[1]) % 3


def state_key(R):
    return sum(1 << (9*a + 3*b + c) for a, b, c in R)


EMPTY = frozenset()
states = [EMPTY]
models = [((1, 0), (1, 0), (1, 1))]
for a, b, c in product(F, repeat=3):
    states.append(frozenset(((a, b, c),)))
    models.append(((1, (b-a) % 3), (1, (c-a) % 3), (-1, 2*a % 3)))
for g, h in product(AFF, repeat=2):
    states.append(frozenset((x, act(g, x), act(h, x)) for x in F))
    models.append((g, h, (1, 0)))
index = {R: i for i, R in enumerate(states)}
require(len(index) == len(states) == 64, "64 distinct states")


def edges(model, labels):
    g, h, loop = model
    a, b, c, u, v = labels
    return ((a, b, g), (a, c, h), (a, u, (1, 0)),
            (u, v, (1, 0)), (v, a, loop))


def graph_relation(es, order, boundary):
    """Independent propagation through literal simple support edges."""
    adj = [[] for _ in range(order)]
    for a, b, g in es:
        adj[a].append((b, g))
        adj[b].append((a, inv(g)))
    result = set()
    for source in F:
        values = [None] * order
        values[0] = source
        stack = [0]
        valid = True
        while stack and valid:
            a = stack.pop()
            for b, g in adj[a]:
                y = act(g, values[a])
                if values[b] is None:
                    values[b] = y
                    stack.append(b)
                elif values[b] != y:
                    valid = False
                    break
        if valid:
            if any(x is None for x in values):
                raise RuntimeError("support must be connected")
            result.add(tuple(values[a] for a in boundary))
    return frozenset(result)


def join(R, S, T):
    result = set()
    for a, x, z in R:
        for xx, b, y in S:
            if x != xx:
                continue
            for zz, yy, c in T:
                if z == zz and y == yy:
                    result.add((a, b, c))
    return frozenset(result)


def run():
    for R, model in zip(states, models):
        es = edges(model, (0, 1, 2, 3, 4))
        require(len({frozenset((a, b)) for a, b, _ in es}) == 5,
                "realization support is simple")
        brute = frozenset(xs[:3] for xs in product(F, repeat=5)
                          if all(xs[b] == act(g, xs[a]) for a, b, g in es))
        require(brute == R, "literal 3^5 realization")
        require(graph_relation(es, 5, (0, 1, 2)) == R, "root propagation")

    # Every two states are separated by one connected singleton pin context.
    for R, S in combinations(states, 2):
        witness = next(iter(R ^ S))
        pin = frozenset((witness,))
        require(pin in index, "pin context is realizable")
        require(bool(R & pin) != bool(S & pin), "context distinguishes states")
        g, h, loop = models[index[pin]]
        # Boundary-to-boundary context edges are subdivided through fresh
        # vertices, so simultaneous gluing cannot introduce parallel edges.
        pin_edges = ((0, 7, (1, 0)), (7, 1, g), (0, 8, (1, 0)), (8, 2, h),
                     (0, 5, (1, 0)), (5, 6, (1, 0)), (6, 0, loop))
        for Q in (R, S):
            union = edges(models[index[Q]], (0, 1, 2, 3, 4)) + pin_edges
            require(len({frozenset((a, b)) for a, b, _ in union}) == len(union),
                    "pin union is a simple support")
            require(graph_relation(union, 9, (0, 1, 2)) == Q & pin,
                    "literal subdivided pin context")

    # Compare all 64^3 compositions against actual 12-cell support propagation.
    hist, table = Counter(), bytearray()
    placements = ((0, 3, 5, 6, 7), (3, 1, 4, 8, 9), (5, 4, 2, 10, 11))
    placed = [[edges(model, placement) for model in models]
              for placement in placements]
    for i, R in enumerate(states):
        for j, S in enumerate(states):
            left = placed[0][i] + placed[1][j]
            for k, T in enumerate(states):
                result = join(R, S, T)
                require(result in index, "closure under triangle gluing")
                actual = graph_relation(left + placed[2][k], 12, (0, 1, 2))
                require(result == actual, "relation join equals literal support")
                hist[len(result)] += 1
                table.append(index[result])

    diagonal = frozenset((x, x, x) for x in F)
    children = {
        "three": diagonal,
        "one": frozenset((x, -x % 3, x) for x in F),
        "zero": frozenset((x, (x+1) % 3, x) for x in F),
    }
    outcomes = {}
    for key, T in children.items():
        require(len(T) == 3, "all child relations have same cardinality")
        outcomes[key] = sorted(join(diagonal, diagonal, T))
    require({key: len(value) for key, value in outcomes.items()} ==
            {"three": 3, "one": 1, "zero": 0}, "3,1,0 scalar hostile")

    # D3 -> O(2): angles theta -> epsilon*theta + delta/3.
    for g, h, theta in product(AFF, AFF, (Fraction(i, 6) for i in range(6))):
        lhs = (g[0] * ((h[0]*theta + Fraction(h[1], 3)) % 1)
               + Fraction(g[1], 3)) % 1
        gh = compose(g, h)
        rhs = (gh[0]*theta + Fraction(gh[1], 3)) % 1
        require(lhs == rhs, "faithful triangle-circle composition")
    subgroups = []
    for bits in range(1 << len(AFF)):
        H = {AFF[i] for i in range(len(AFF)) if bits >> i & 1}
        if (1, 0) not in H or any(compose(g, h) not in H for g in H for h in H):
            continue
        fixed_triangle = [x for x in F if all(act(g, x) == x for g in H)]
        if H == {(1, 0)}:
            fixed_circle = "whole circle"
            circle_nonempty = True
        else:
            # Any nonidentity rotation has no fixed point. Otherwise any
            # reflection confines fixed points to its two antipodal solutions.
            if any(a == 1 and b != 0 for a, b in H):
                circle_points = []
            else:
                g = next(g for g in H if g[0] == -1)
                candidates = [Fraction(g[1], 6), (Fraction(g[1], 6)+Fraction(1, 2)) % 1]
                circle_points = [x for x in candidates
                                 if all((a*x+Fraction(b, 3)) % 1 == x for a, b in H)]
            fixed_circle = [str(x) for x in circle_points]
            circle_nonempty = bool(circle_points)
        require(bool(fixed_triangle) == circle_nonempty, "section existence transfers")
        subgroups.append({"elements": sorted(H), "triangle_fixed": fixed_triangle,
                          "circle_fixed": fixed_circle})
    require(len(subgroups) == 6, "all six subgroups of S3")

    # General boundary size formula; do not extrapolate powers of four.
    totals = {b: 1 + 3**b + 6**(b-1) for b in range(1, 6)}
    require(totals == {1: 5, 2: 16, 3: 64, 4: 298, 5: 1540}, "boundary count")
    gauge_orbits = set()
    for R in states:
        orbit = frozenset(state_key(frozenset(tuple((x+a) % 3 for x, a in zip(xs, shift))
                                                   for xs in R))
                          for shift in product(F, repeat=3))
        gauge_orbits.add(orbit)
    require(sorted(map(len, gauge_orbits)) == [1, 9, 9, 9, 9, 27],
            "six coarse states after independent boundary gauge")

    def boundary_count(r, b):
        return (sum(comb(r, k)*(factorial(r)//factorial(r-k))**(b-1)
                    for k in range(r-1)) + factorial(r)**(b-1))

    # Independent small literal enumeration of all permitted partial bijections.
    permutation_states = set()
    for k in (0, 1, 2, 4):
        for domain in combinations(range(4), k):
            for values in permutations(range(4), k):
                permutation_states.add(frozenset(zip(domain, values)))
    require(len(permutation_states) == boundary_count(4, 2) == 113,
            "four-label two-boundary partial-bijection count")
    seven_counts = {b: boundary_count(7, b) for b in (1, 2, 3)}
    require(seven_counts == {1: 121, 2: 95642, 3: 185036888}, "seven-label formula")
    require(2**27 < seven_counts[3] <= 2**28, "seven-label 28-bit minimum")
    return {"scope": "connected affine-F3 matching supports; actual boundary coordinates",
            "state_count": len(states), "state_cardinalities": dict(sorted(Counter(map(len, states)).items())),
            "state_masks": [state_key(R) for R in states],
            "all_compositions": len(table), "composition_cardinalities": dict(sorted(hist.items())),
            "composition_table_sha256": sha256(table).hexdigest(),
            "literal_support_cells_per_composition": 12,
            "scalar_hostile": outcomes, "triangle_circle_subgroups": subgroups,
            "general_boundary_counts": totals,
            "independent_boundary_gauge_orbit_sizes": sorted(map(len, gauge_orbits)),
            "seven_label_boundary_counts_by_proved_formula": seven_counts,
            "four_label_two_boundary_enumerated_states": len(permutation_states), "checks": COUNT}


if __name__ == "__main__":
    receipt = run()
    data = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    output = HERE / "receipt.json"
    if "--write" in sys.argv:
        output.write_text(data, encoding="utf-8")
    else:
        require(json.loads(output.read_text(encoding="utf-8")) == json.loads(data), "frozen receipt replay")
    print(json.dumps({"checks": COUNT, "states": 64, "compositions": receipt["all_compositions"],
                      "composition_cardinalities": receipt["composition_cardinalities"],
                      "table_sha256": receipt["composition_table_sha256"]}, sort_keys=True))
