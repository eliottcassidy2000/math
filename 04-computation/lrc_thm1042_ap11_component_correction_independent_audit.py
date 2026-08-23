#!/usr/bin/env python3
"""Independent exact wall-graph audit of THM-1042's B={1,...,11} row.

This implementation does not import the canonical component-length script.  It
distinguishes strict-safe open cells from the connected components of the
closed safe set, including zero-length components supported at wall points.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256


DELTA = Q(1, 14)
BASE = tuple(range(1, 12))
CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def circle_distance(x):
    r = x % 1
    return min(r, 1 - r)


def safe_closed(x):
    return all(circle_distance(v * x) >= DELTA for v in BASE)


def safe_strict(x):
    return all(circle_distance(v * x) > DELTA for v in BASE)


def frac(x):
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def cyclic_midpoint(a, b):
    """Midpoint of the positively oriented circle cell a -> b."""
    lift_b = b if b > a else b + 1
    return ((a + lift_b) / 2) % 1


def main():
    # All equality walls ||v t||=1/14, reduced modulo one.  Zero is inserted
    # only as a display/cyclic cut, not counted as an equality wall.
    equality_walls = {
        ((Q(k, 1) + sign * DELTA) / v) % 1
        for v in BASE
        for k in range(v)
        for sign in (-1, 1)
    }
    walls = sorted(equality_walls)
    n = len(walls)

    # Cyclic open cells between consecutive equality walls.
    cell_safe = []
    for i, a in enumerate(walls):
        b = walls[(i + 1) % n]
        m = cyclic_midpoint(a, b)
        # No wall lies inside a cell, so closed/strict tests coincide there.
        require(safe_closed(m) == safe_strict(m), f"cell convention {i}")
        cell_safe.append(safe_strict(m))

    wall_safe = [safe_closed(x) for x in walls]

    # Bipartite incidence graph: safe wall vertices and safe open-cell
    # vertices; a safe cell is incident to its two closed endpoints.  Its
    # connected components are exactly those of the closed safe subset.
    vertices = [("w", i) for i in range(n) if wall_safe[i]]
    vertices += [("c", i) for i in range(n) if cell_safe[i]]
    adjacency = {v: set() for v in vertices}
    for i in range(n):
        if not cell_safe[i]:
            continue
        for wi in (i, (i + 1) % n):
            require(wall_safe[wi], f"safe cell endpoint {i}:{wi}")
            adjacency[("c", i)].add(("w", wi))
            adjacency[("w", wi)].add(("c", i))

    components = []
    unseen = set(vertices)
    while unseen:
        seed = min(unseen)
        stack = [seed]
        comp = set()
        unseen.remove(seed)
        while stack:
            x = stack.pop()
            comp.add(x)
            for y in adjacency[x]:
                if y in unseen:
                    unseen.remove(y)
                    stack.append(y)
        components.append(comp)

    isolated = sorted(
        walls[i]
        for i in range(n)
        if wall_safe[i] and not adjacency[("w", i)]
    )

    # Extract the nondegenerate components by their extreme incident cell
    # endpoints.  In this census each has exactly one safe cell, but the code
    # computes its cyclic length from all cells and checks that structure.
    arcs = []
    for comp in components:
        cell_ids = sorted(i for typ, i in comp if typ == "c")
        if not cell_ids:
            continue
        # Sum cyclic cell lengths; shared safe walls have zero measure.
        length = Q(0)
        for i in cell_ids:
            a, b = walls[i], walls[(i + 1) % n]
            length += (b - a) % 1
        arcs.append((walls[cell_ids[0]], walls[(cell_ids[-1] + 1) % n], length))

    lengths = sorted(length for _, _, length in arcs)
    measure = sum(lengths, Q(0))

    # Machine-readable semantic packet independent of presentation order.
    semantic = "\n".join(
        [
            "base=" + ",".join(map(str, BASE)),
            "endpoint=closed-ge",
            "walls=" + ",".join(map(frac, walls)),
            "safe_cells=" + ",".join(str(i) for i, x in enumerate(cell_safe) if x),
            "safe_walls=" + ",".join(frac(walls[i]) for i, x in enumerate(wall_safe) if x),
            "isolated=" + ",".join(map(frac, isolated)),
            "lengths=" + ",".join(map(frac, lengths)),
        ]
    )

    require(n == 128, "equality wall count")
    require(sum(cell_safe) == 14, "strict safe cell count")
    require(sum(wall_safe) == 32, "closed safe wall count")
    require(len(components) == 18, "closed component count")
    require(len(arcs) == 14, "nondegenerate arc count")
    require(isolated == [Q(3, 14), Q(5, 14), Q(9, 14), Q(11, 14)], "isolated walls")
    require(measure == Q(10931, 194040), "safe measure")
    require(max(lengths) == Q(1, 77), "maximum component")
    require(Counter(lengths) == Counter(
        {
            Q(1, 588): 2,
            Q(1, 560): 2,
            Q(1, 504): 2,
            Q(1, 420): 2,
            Q(1, 308): 2,
            Q(1, 245): 2,
            Q(1, 77): 2,
        }
    ), "positive component spectrum")

    print("THM1042_AP11_INDEPENDENT_WALL_GRAPH_20260823")
    print("universe=B={1,...,11};circle=R/Z;delta=1/14")
    print("danger=open:{t:||vt||<1/14};safe=closed:{t:all ||vt||>=1/14}")
    print(f"equality_walls={n};closed_safe_walls={sum(wall_safe)}")
    print(f"strict_safe_open_cells={sum(cell_safe)}")
    print(f"closed_topological_components={len(components)}")
    print(f"closed_nondegenerate_arcs={len(arcs)};closed_isolated_points={len(isolated)}")
    print("isolated_points=" + ",".join(map(frac, isolated)))
    print("positive_component_lengths=" + ",".join(map(frac, lengths)))
    print(f"measure={frac(measure)};max_positive_component={frac(max(lengths))};threshold={frac(1/max(lengths))}")
    print("positive_length_convention_components=14;closed_topological_convention_components=18")
    print("semantic_sha256=" + sha256(semantic.encode("ascii")).hexdigest())
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
