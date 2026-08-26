#!/usr/bin/env python3
"""Exact audit for THM-4221's cycle-left source-padding response sector.

The inherited THM-4184 engine constructs capacities from literal exposed
two-path covers.  This script additionally rebuilds both ordinal remainders
from the literal child tournaments, enumerates marked Hamilton-path arcs, and
checks the endpoint-chirality improvement used in the all-order proof.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations

import tournament_ordinal_cocycle_parity_thm4184 as base


def need(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def tournament_from_word(order: int, word: int) -> tuple[int, ...]:
    out = [0] * order
    cursor = 0
    for left in range(order):
        for right in range(left + 1, order):
            if (word >> cursor) & 1:
                out[left] |= 1 << right
            else:
                out[right] |= 1 << left
            cursor += 1
    return tuple(out)


def cycle3() -> tuple[int, ...]:
    return (1 << 1, 1 << 2, 1 << 0)


def tower_one(order: int) -> tuple[int, ...]:
    """The THM-4219 strong hostile T(order,1)."""
    need(order >= 5, "tower order")
    last = order - 1
    special = {0, order - 3}
    out = [0] * order
    for left in range(last):
        for right in range(left + 1, last):
            out[left] |= 1 << right
    for vertex in range(last):
        if vertex in special:
            out[last] |= 1 << vertex
        else:
            out[vertex] |= 1 << last
    return tuple(out)


def hamilton_paths(out: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        path
        for path in permutations(range(len(out)))
        if all(out[path[k]] & (1 << path[k + 1]) for k in range(len(out) - 1))
    )


def literal_remainder(
    left: base.TournamentData, right: base.TournamentData
) -> int:
    child = base.tournament_data(base.ordinal_out(left.out, right.out))
    return (
        child.gate
        - right.hamilton**2 * left.gate
        - left.hamilton**2 * right.gate
    )


def coordinates(
    data: base.TournamentData,
) -> tuple[
    int,
    int,
    tuple[int, ...],
    tuple[int, ...],
    tuple[int, ...],
    tuple[int, ...],
    tuple[int, ...],
]:
    order = len(data.out)
    degree = tuple(sum(data.capacities[i]) for i in range(order))
    incoming = tuple(
        sum(
            data.capacities[i][j]
            for j in range(order)
            if j != i and data.out[j] & (1 << i)
        )
        for i in range(order)
    )
    v = tuple(sum(data.ends[i]) for i in range(order))
    starts = tuple(data.ends[i][1] - data.ends[i][0] for i in range(order))
    ends = tuple(data.starts[i][1] - data.starts[i][0] for i in range(order))
    need(all(v[i] == incoming[i] + starts[i] for i in range(order)), "fan collapse")
    return data.hamilton, data.mass, v, starts, ends, degree, incoming


def increment_formula(
    data: base.TournamentData,
) -> tuple[int, int, int, int, int, int, tuple[int, ...], tuple[int, ...]]:
    order = len(data.out)
    h, z, v, starts, ends, degree, incoming = coordinates(data)
    mass_square = sum(value * value for value in v)
    tau = sum(value * value for value in starts)
    delta = (z + h) * (z + 3 * h) - 3 * mass_square
    auxiliary = sum(
        (z - degree[i]) * (27 * v[i] - 9 * starts[i]) for i in range(order)
    )
    first = (
        324 * z * z
        + 1188 * z * h
        + 744 * h * h
        - 945 * mass_square
        - 15 * tau
        - 27 * sum(degree[i] * v[i] for i in range(order))
        + 9 * sum(degree[i] * starts[i] for i in range(order))
    )
    second = (
        -18 * z * z
        - 90 * z * h
        - 201 * h * h
        + 315 * delta
        - 15 * tau
        + auxiliary
    )
    need(first == second, "two increment forms")
    return first, delta, auxiliary, mass_square, tau, h, v, starts


def marked_arc_audit(
    data: base.TournamentData,
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    out = data.out
    order = len(out)
    paths = hamilton_paths(out)
    starts = tuple(sum(path[0] == i for path in paths) for i in range(order))
    ends = tuple(sum(path[-1] == i for path in paths) for i in range(order))
    marked = [[0] * order for _ in range(order)]
    for path in paths:
        for k in range(order - 1):
            left, right = sorted((path[k], path[k + 1]))
            marked[left][right] += 1
    for left in range(order):
        for right in range(left + 1, order):
            need(
                marked[left][right] <= data.capacities[left][right],
                "marked arc injection",
            )
    avoiding = tuple(
        sum(
            marked[left][right]
            for left in range(order)
            for right in range(left + 1, order)
            if left != vertex and right != vertex
        )
        for vertex in range(order)
    )
    expected = tuple(
        (order - 3) * len(paths) + starts[i] + ends[i] for i in range(order)
    )
    need(avoiding == expected, "marked count identity")
    return starts, ends, avoiding


def padding_audit(
    data: base.TournamentData, one: base.TournamentData
) -> base.TournamentData:
    padded_out = base.ordinal_out(one.out, data.out)
    literal = base.tournament_data(padded_out)
    transferred = base.ordinal_data(one, data)
    need(literal == transferred, "literal source-padding transfer")
    h, z, v, starts, _ends, degree, incoming = coordinates(data)
    ph, pz, pv, _pstarts, _pends, pdegree, pincoming = coordinates(literal)
    need(ph == h and pz == 2 * (z + h), "padded scalar transfer")
    need(literal.ends[0] == (0, h), "padded source V state")
    need(
        all(literal.ends[i + 1] == (v[i], v[i]) for i in range(len(v))),
        "padded old V states",
    )
    need(pincoming[0] == 0 and pdegree[0] == z + 2 * h, "padded source fan")
    need(
        all(pincoming[i + 1] == 2 * v[i] for i in range(len(v))),
        "padded old incoming fans",
    )
    need(
        all(
            pdegree[i + 1] == degree[i] + incoming[i] + 2 * starts[i]
            for i in range(len(v))
        ),
        "padded old degrees",
    )
    need(pv[0] == h and all(pv[i + 1] == 2 * v[i] for i in range(len(v))), "padded v")
    return literal


def audit_labelled() -> None:
    one = base.tournament_data(base.transitive(1))
    cycle = base.tournament_data(cycle3())
    total = 0
    no_sink = 0
    by_order: list[int] = []
    sector_by_order: list[int] = []
    conditional_by_order: list[int] = []
    minimum_ratio: dict[int, Fraction] = {}
    equality = 0
    for order in range(1, 6):
        rows = 1 << (order * (order - 1) // 2)
        order_no_sink = 0
        order_sector = 0
        order_conditional = 0
        for word in range(rows):
            data = base.tournament_data(tournament_from_word(order, word))
            padded = padding_audit(data, one)
            increment, delta, auxiliary, _mass_square, tau, h, v, starts = increment_formula(data)
            literal_before = literal_remainder(cycle, data)
            literal_after = literal_remainder(cycle, padded)
            need(literal_after - literal_before == increment, "literal response increment")
            need(literal_before == base.remainder(cycle, data), "literal base remainder")
            need(literal_after == base.remainder(cycle, padded), "literal padded remainder")
            total += 1
            if any(row == 0 for row in data.out):
                continue

            no_sink += 1
            order_no_sink += 1
            path_starts, path_ends, avoiding = marked_arc_audit(data)
            _h, z, _v, coordinate_starts, coordinate_ends, degree, _incoming = coordinates(data)
            need(path_starts == coordinate_starts, "Hamilton starts")
            need(path_ends == coordinate_ends, "Hamilton ends")
            need(all(v[i] >= h for i in range(order)), "split floor")
            need(delta >= order * (order - 1) * h * h, "endpoint-energy floor")
            need(
                all(z - degree[i] >= avoiding[i] for i in range(order)),
                "avoidance capacity bound",
            )

            full = (1 << order) - 1
            sources = [i for i in range(order) if data.out[i] == (full ^ (1 << i))]
            need(all(2 * path_ends[i] <= h for i in range(order)), "no-sink terminal chirality")
            if sources:
                need(len(sources) == 1, "unique universal source")
                source = sources[0]
                need(path_starts[source] == h and path_ends[source] == 0, "source endpoint boundary")
            else:
                need(all(2 * path_starts[i] <= h for i in range(order)), "initial chirality")
            endpoint_cross = sum(path_starts[i] * path_ends[i] for i in range(order))
            need(tau + endpoint_cross <= h * h, "joint endpoint-energy sidecar")

            auxiliary_floor = 9 * (order - 3) * h * (3 * z + 2 * h) + 45 * h * h
            need(auxiliary >= auxiliary_floor, "improved auxiliary floor")
            lower = (
                -18 * z * z
                + (27 * order - 171) * z * h
                + (315 * order * order - 297 * order - 225) * h * h
            )
            need(increment >= lower, "unconditional response lower bound")
            need(z >= (order - 1) * h, "one-defect mass floor")
            ratio = Fraction(increment, h * h)
            minimum_ratio[order] = min(minimum_ratio.get(order, ratio), ratio)

            if order >= 4 and lower >= 1480 * h * h:
                order_sector += 1
                need(increment >= 1480 * h * h, "exact radical sector")
            if order >= 4 and z <= (5 * order - 11) * h:
                need(increment > 1480 * h * h, "clean strict sector")
            if 8 * delta >= 3 * (z + 2 * h) ** 2:
                order_conditional += 1
                need(increment >= 1480 * h * h, "conditional gap consequence")
            if increment == 1480 * h * h:
                equality += 1
                need(order == 3, "response equality only C3 in finite audit")
        by_order.append(order_no_sink)
        sector_by_order.append(order_sector)
        conditional_by_order.append(order_conditional)

    need(total == 1099, "labelled universe size")
    need(no_sink == 738, "labelled no-sink size")
    need(equality == 2, "two labelled C3 equalities")
    shown = {order: str(value) for order, value in sorted(minimum_ratio.items())}
    print(f"labelled_orders_1_to_5={total} source_padding_and_literal_response=PASS")
    print(f"no_sink={no_sink} no_sink_by_order={by_order}")
    print(f"radical_sector_by_order={sector_by_order} conditional_rows_by_order={conditional_by_order}")
    print(f"min_increment_over_H2={shown} equality_labelled={equality}")
    print("marked_arc_injection=PASS endpoint_chirality_45H2=PASS lower_sector=PASS")


def controls() -> None:
    cycle = base.tournament_data(cycle3())
    c_increment = increment_formula(cycle)[0]
    need(c_increment == 13_320 == 1480 * cycle.hamilton**2, "C3 equality control")

    transitive = base.tournament_data(base.transitive(3))
    transitive_increment = increment_formula(transitive)[0]
    need(any(row == 0 for row in transitive.out), "transitive filter hostile")

    tower = base.tournament_data(tower_one(9))
    tower_increment = increment_formula(tower)[0]
    need(not any(row == 0 for row in tower.out), "tower no-sink control")
    need(tower.mass > (5 * len(tower.out) - 11) * tower.hamilton, "tower outside clean sector")
    print(
        "controls="
        f"C3_increment={c_increment};"
        f"transitive3_filtered_increment={transitive_increment};"
        f"T1_n9_outside_clean_increment={tower_increment}"
    )


def algebra_controls() -> None:
    for order in range(4, 200):
        clean_slack = 531 * order - 2002
        need(clean_slack > 0, "clean linear sector algebra")
        conditional_numerator = 1017 * order * order + 738 * order + 369
        need(conditional_numerator > 8 * 1480, "conditional all-order algebra")
        radicand = 2601 * order * order - 3402 * order - 10391
        need(radicand > 0, "radical sector is real")
    print("symbolic_integer_reductions=clean_slack_531n-2002;conditional_numerator_1017n2+738n+369")


def main() -> None:
    audit_labelled()
    controls()
    algebra_controls()
    print("THM4221_EXACT_AUDIT=PASS")


if __name__ == "__main__":
    main()
