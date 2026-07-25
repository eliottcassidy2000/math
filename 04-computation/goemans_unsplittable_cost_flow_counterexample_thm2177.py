#!/usr/bin/env python3
"""Exact referee for THM-2177.

The checker starts from the directed graph, not from a prescribed path list.
It independently enumerates every simple source--terminal path, every
unsplittable routing, its arc loads, its cost, and every violated additive
capacity constraint.  All arithmetic is integral or rational.

Reproduce:
    python3 04-computation/goemans_unsplittable_cost_flow_counterexample_thm2177.py
    python3 -O 04-computation/goemans_unsplittable_cost_flow_counterexample_thm2177.py
"""

from fractions import Fraction
from itertools import product


Vertex = str
Arc = tuple[Vertex, Vertex]
Path = tuple[Arc, ...]

SOURCE = "s"
TERMINALS = ("t1", "t2", "t3")
DEMAND = {"t1": 15, "t2": 10, "t3": 15}
MAX_DEMAND = max(DEMAND.values())

# arc: (fractional load, nonnegative per-unit cost)
ARC_DATA: dict[Arc, tuple[int, int]] = {
    ("s", "t1"): (10, 2),
    ("s", "t2"): (6, 3),
    ("s", "u"): (24, 0),
    ("u", "t3"): (10, 2),
    ("u", "v"): (14, 0),
    ("v", "t1"): (5, 0),
    ("v", "w"): (9, 0),
    ("w", "t2"): (4, 0),
    ("w", "t3"): (5, 0),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def all_simple_paths(source: Vertex, target: Vertex) -> tuple[Path, ...]:
    adjacency: dict[Vertex, list[Vertex]] = {}
    for tail, head in ARC_DATA:
        adjacency.setdefault(tail, []).append(head)

    paths: list[Path] = []

    def visit(vertex: Vertex, seen: frozenset[Vertex], path: Path) -> None:
        if vertex == target:
            paths.append(path)
            return
        for head in sorted(adjacency.get(vertex, ())):
            if head not in seen:
                visit(head, seen | {head}, path + ((vertex, head),))

    visit(source, frozenset({source}), ())
    return tuple(paths)


def fractional_flow_conservation() -> None:
    vertices = {endpoint for arc in ARC_DATA for endpoint in arc}
    for vertex in sorted(vertices):
        inflow = sum(
            load for (tail, head), (load, _) in ARC_DATA.items() if head == vertex
        )
        outflow = sum(
            load for (tail, head), (load, _) in ARC_DATA.items() if tail == vertex
        )
        expected_divergence = (
            sum(DEMAND.values())
            if vertex == SOURCE
            else -DEMAND.get(vertex, 0)
        )
        require(
            outflow - inflow == expected_divergence,
            f"flow conservation failed at {vertex}",
        )


def path_cost(path: Path, demand: int) -> int:
    return demand * sum(ARC_DATA[arc][1] for arc in path)


def routing_loads(paths: tuple[Path, ...]) -> dict[Arc, int]:
    loads = {arc: 0 for arc in ARC_DATA}
    for terminal, path in zip(TERMINALS, paths):
        for arc in path:
            loads[arc] += DEMAND[terminal]
    return loads


def routing_cost(loads: dict[Arc, int]) -> int:
    return sum(loads[arc] * ARC_DATA[arc][1] for arc in ARC_DATA)


def violations(loads: dict[Arc, int]) -> tuple[tuple[Arc, int], ...]:
    return tuple(
        (arc, loads[arc] - (fractional_load + MAX_DEMAND))
        for arc, (fractional_load, _) in ARC_DATA.items()
        if loads[arc] > fractional_load + MAX_DEMAND
    )


def path_vertices(path: Path) -> tuple[Vertex, ...]:
    if not path:
        return (SOURCE,)
    return (path[0][0],) + tuple(head for _, head in path)


def family_specialization_control() -> None:
    """Check the rational point in the three-parameter family exactly."""
    b = Fraction(2, 3)
    r = Fraction(1, 3)
    q = Fraction(2, 5)
    require(2 * r + q > 1, "family cost inequality failed")
    require(b * (1 - q) > r, "family outer-conflict inequality failed")
    require(2 * r + b * q < 1, "family middle-conflict inequality failed")

    normalized_loads: dict[Arc, Fraction] = {
        ("s", "t1"): 1 - r,
        ("s", "t2"): b * (1 - q),
        ("s", "u"): 1 + r + b * q,
        ("u", "t3"): 1 - r,
        ("u", "v"): 2 * r + b * q,
        ("v", "t1"): r,
        ("v", "w"): r + b * q,
        ("w", "t2"): b * q,
        ("w", "t3"): r,
    }
    require(
        {
            arc: 15 * load for arc, load in normalized_loads.items()
        }
        == {
            arc: Fraction(load)
            for arc, (load, _) in ARC_DATA.items()
        },
        "integer instance is not the displayed family specialization",
    )
    normalized_cost = 3 - (2 * r + q)
    require(normalized_cost == Fraction(29, 15), "family cost changed")
    require(15 * 2 * normalized_cost == 58, "scaled family cost changed")


def main() -> None:
    require(MAX_DEMAND == 15, "maximum demand changed")
    require(all(cost >= 0 for _, cost in ARC_DATA.values()), "negative cost")
    fractional_flow_conservation()
    family_specialization_control()

    paths_by_terminal = {
        terminal: all_simple_paths(SOURCE, terminal) for terminal in TERMINALS
    }
    require(
        tuple(len(paths_by_terminal[t]) for t in TERMINALS) == (2, 2, 2),
        "the graph does not have exactly two paths per terminal",
    )

    # E_i is the unique positive-cost path and Z_i the zero-cost path.
    labelled: dict[str, tuple[Vertex, ...]] = {}
    ordered_paths: list[tuple[Path, Path]] = []
    for index, terminal in enumerate(TERMINALS, start=1):
        paths = paths_by_terminal[terminal]
        expensive = tuple(
            path for path in paths if path_cost(path, DEMAND[terminal]) > 0
        )
        zero = tuple(
            path for path in paths if path_cost(path, DEMAND[terminal]) == 0
        )
        require(
            len(expensive) == len(zero) == 1,
            f"path cost split failed for {terminal}",
        )
        require(
            path_cost(expensive[0], DEMAND[terminal]) == 30,
            f"expensive path cost changed for {terminal}",
        )
        labelled[f"E{index}"] = path_vertices(expensive[0])
        labelled[f"Z{index}"] = path_vertices(zero[0])
        ordered_paths.append((expensive[0], zero[0]))

    expected_paths = {
        "E1": ("s", "t1"),
        "Z1": ("s", "u", "v", "t1"),
        "E2": ("s", "t2"),
        "Z2": ("s", "u", "v", "w", "t2"),
        "E3": ("s", "u", "t3"),
        "Z3": ("s", "u", "v", "w", "t3"),
    }
    require(labelled == expected_paths, "simple-path enumeration changed")

    fractional_cost = sum(
        load * cost for load, cost in ARC_DATA.values()
    )
    require(fractional_cost == 58, "fractional cost changed")

    records = []
    for choices in product((0, 1), repeat=3):
        selected_paths = tuple(
            ordered_paths[index][choice]
            for index, choice in enumerate(choices)
        )
        names = tuple(
            ("E" if choice == 0 else "Z") + str(index + 1)
            for index, choice in enumerate(choices)
        )
        loads = routing_loads(selected_paths)
        bad = violations(loads)
        records.append((names, routing_cost(loads), bad))

    good = tuple(record for record in records if not record[2])
    bad = tuple(record for record in records if record[2])
    require(len(records) == 8, "routing universe is not 2^3")
    require(len(good) == 4 and len(bad) == 4, "good/bad routing split changed")
    require(
        all(sum(name.startswith("Z") for name in names) <= 1 for names, _, _ in good),
        "a good routing uses two zero-cost paths",
    )
    require(
        all(sum(name.startswith("Z") for name in names) >= 2 for names, _, _ in bad),
        "a routing with at most one zero-cost path was rejected",
    )

    minimum_good_cost = min(cost for _, cost, _ in good)
    require(minimum_good_cost == 60, "minimum good cost changed")
    require(
        minimum_good_cost > fractional_cost,
        "the cost conjecture is not separated",
    )

    exact_bad_witnesses = {
        ("E1", "Z2", "Z3"): {(("v", "w"), 1)},
        ("Z1", "E2", "Z3"): {(("u", "v"), 1)},
        ("Z1", "Z2", "E3"): {(("s", "u"), 1)},
        ("Z1", "Z2", "Z3"): {
            (("s", "u"), 1),
            (("u", "v"), 11),
            (("v", "w"), 1),
        },
    }
    require(
        {
            names: set(witnesses)
            for names, _, witnesses in bad
        }
        == exact_bad_witnesses,
        "violating-arc certificate changed",
    )

    # Independent convex-hull obstruction: the fractional path weights of
    # Z_1,Z_2,Z_3 violate the triangle stable-set inequality.
    cheap_probabilities = (
        Fraction(5, 15),
        Fraction(4, 10),
        Fraction(5, 15),
    )
    require(
        sum(cheap_probabilities) == Fraction(16, 15) > 1,
        "triangle stable-set separation changed",
    )

    # Derive the underlying undirected graph from ARC_DATA, suppress its three
    # degree-two terminal vertices, and recover K_4 on s,u,v,w.
    underlying_edges = {frozenset(arc) for arc in ARC_DATA}
    for terminal in TERMINALS:
        incident = {
            edge for edge in underlying_edges if terminal in edge
        }
        require(
            len(incident) == 2,
            f"{terminal} is not a degree-two subdivision vertex",
        )
        neighbours = {
            next(iter(edge - {terminal})) for edge in incident
        }
        require(
            len(neighbours) == 2,
            f"{terminal} does not have two distinct neighbours",
        )
        underlying_edges.difference_update(incident)
        underlying_edges.add(frozenset(neighbours))

    k4_edges = {
        frozenset(("s", "u")),
        frozenset(("u", "v")),
        frozenset(("v", "w")),
        frozenset(("s", "v")),
        frozenset(("s", "w")),
        frozenset(("u", "w")),
    }
    require(
        underlying_edges == k4_edges,
        "suppressing the terminals does not recover K4",
    )

    print("THM-2177 exact planar Goemans counterexample")
    print(f"vertices=7, arcs={len(ARC_DATA)}, terminals=3")
    print(f"demands={(15, 10, 15)}, d_max={MAX_DEMAND}")
    print("simple_paths_per_terminal=(2, 2, 2)")
    print("all_unsplittable_routings=8")
    print("capacity_good_routings=4, capacity_bad_routings=4")
    print(f"fractional_cost={fractional_cost}")
    print(f"minimum_capacity_good_cost={minimum_good_cost}")
    print("cost_gap=2")
    print("cheap_path_probability_sum=16/15")
    print("family_point=(b,r,q)=(2/3,1/3,2/5)")
    print("underlying_undirected_graph=K4_subdivision, hence planar")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
