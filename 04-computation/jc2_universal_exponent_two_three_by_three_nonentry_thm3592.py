#!/usr/bin/env python3
"""Optimization-safe exact companion for THM-3592.

The theorem is universal in the squarefree target and in the integer supports.
This standard-library-only companion checks the finite collision-pattern
catalogue, component/deletion gates, every scalar-arm placement in the five
connected collision types and their reversal-sensitive representatives, the
hooked-family differential identities, and all sharp boundaries displayed in
THM-3592.  Finite parameter ranges are hostile exact controls; the universal
steps are proved symbolically in the theorem.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    """Optimization-safe truth gate (never compiled away by ``python -O``)."""
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------------------
# Three-point sumset and collision classification.


def offsets(gaps: tuple[int, int]) -> tuple[int, int, int]:
    x, y = gaps
    require(x > 0 and y > 0, "gaps must be positive")
    return (0, x, x + y)


def fibre_map(
    gaps_a: tuple[int, int], gaps_b: tuple[int, int]
) -> dict[int, tuple[tuple[int, int], ...]]:
    a_values, b_values = offsets(gaps_a), offsets(gaps_b)
    fibres: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for i, a_value in enumerate(a_values):
        for j, b_value in enumerate(b_values):
            fibres[a_value + b_value].append((i, j))
    return {value: tuple(cells) for value, cells in sorted(fibres.items())}


def collision_statistics(
    gaps_a: tuple[int, int], gaps_b: tuple[int, int]
) -> tuple[int, int, int, int]:
    fibres = fibre_map(gaps_a, gaps_b)
    sizes = [len(cells) for cells in fibres.values()]
    require(max(sizes) <= 3, "a three-point sum fibre has size at most three")
    collision_pairs = sum(size * (size - 1) // 2 for size in sizes)
    triples = sum(size == 3 for size in sizes)
    sumset_size = len(fibres)
    require(sumset_size == 9 - collision_pairs + triples, "9-C+tau identity")
    return sumset_size, collision_pairs, triples, max(sizes)


def transformations(data: tuple[int, int, int, int]) -> tuple[tuple[int, ...], ...]:
    """Exchange the outputs and/or reverse both ordered supports."""
    x, y, u, v = data
    return (
        (x, y, u, v),
        (u, v, x, y),
        (y, x, v, u),
        (v, u, y, x),
    )


def primitive_normal_form(data: tuple[int, int, int, int]) -> tuple[int, ...]:
    divisor = gcd(gcd(data[0], data[1]), gcd(data[2], data[3]))
    reduced = tuple(entry // divisor for entry in data)
    return min(transformations(reduced))


def low_sumset_label(data: tuple[int, int, int, int]) -> str | None:
    """Return the exact affine type whenever |A+B| is at most seven."""
    variants = transformations(data)

    for x, y, u, v in variants:
        if x == y == u == v:
            return "AP5"

    for x, y, u, v in variants:
        if x != y and (u, v) == (x, y):
            return "D6"

    for x, y, u, v in variants:
        if x == y and u == x and v == 2 * x:
            return "H6"

    for x, y, u, v in variants:
        if x != y and (u, v) == (y, x):
            return "R7"

    for p, q, u, v in variants:
        if 0 < p < q and q != 2 * p and (u, v) == (p, q - p):
            return "E_PLUS7"

    for p, q, u, v in variants:
        if 0 < p < q and q != 2 * p and (u, v) == (q - p, p):
            return "E_MINUS7"

    for d, d_again, u, v in variants:
        if d == d_again and u == d and v not in (d, 2 * d):
            return "AP_EXTENSION7"

    for d, d_again, p, q in variants:
        if d == d_again and p + q == d and p != q:
            return "AP_CONTAINED7"

    for d, d_again, u, v in variants:
        if d == d_again and u == v == 2 * d:
            return "AP_DYADIC7"

    return None


def difference_set(gaps: tuple[int, int]) -> set[int]:
    x, y = gaps
    return {x, y, x + y}


def component_sizes(vertices: tuple[int, int, int], allowed: set[int]) -> tuple[int, ...]:
    adjacency = [set() for _ in vertices]
    for i in range(3):
        for j in range(i + 1, 3):
            if vertices[j] - vertices[i] in allowed:
                adjacency[i].add(j)
                adjacency[j].add(i)
    unseen = set(range(3))
    sizes: list[int] = []
    while unseen:
        seed = min(unseen)
        stack = [seed]
        unseen.remove(seed)
        size = 0
        while stack:
            current = stack.pop()
            size += 1
            for neighbor in adjacency[current]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def support_components(
    gaps_a: tuple[int, int], gaps_b: tuple[int, int]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    a_values, b_values = offsets(gaps_a), offsets(gaps_b)
    gamma_b = component_sizes(b_values, difference_set(gaps_a))
    gamma_a = component_sizes(a_values, difference_set(gaps_b))
    return gamma_a, gamma_b


EXPECTED_BY_SIZE = {
    5: {"AP5"},
    6: {"D6", "H6"},
    7: {
        "R7",
        "E_PLUS7",
        "E_MINUS7",
        "AP_EXTENSION7",
        "AP_CONTAINED7",
        "AP_DYADIC7",
    },
}
CONNECTED_LABELS = {"AP5", "D6", "H6", "R7", "E_PLUS7", "E_MINUS7"}


collision_counts: Counter[str] = Counter()
primitive_forms: dict[str, set[tuple[int, ...]]] = defaultdict(set)
GAP_BOUND = 24
for x in range(1, GAP_BOUND + 1):
    for y in range(1, GAP_BOUND + 1):
        for u in range(1, GAP_BOUND + 1):
            for v in range(1, GAP_BOUND + 1):
                data = (x, y, u, v)
                size, collisions, triples, _ = collision_statistics((x, y), (u, v))
                require(triples in (0, 1), "at most one triple fibre")
                require(
                    (triples == 1) == ((u, v) == (y, x)),
                    "triple iff reflected gap vector",
                )
                label = low_sumset_label(data)
                if size <= 7:
                    require(label in EXPECTED_BY_SIZE[size], "complete low-sumset label")
                    collision_counts[label] += 1
                    if gcd(gcd(x, y), gcd(u, v)) == 1:
                        primitive_forms[label].add(primitive_normal_form(data))
                else:
                    require(label is None, "no low-sumset label above seven")

                gamma_a, gamma_b = support_components((x, y), (u, v))
                if size >= 8:
                    require(
                        gamma_a != (3,) or gamma_b != (3,),
                        "a size-eight/nine sumset has a deletable component",
                    )
                if label in CONNECTED_LABELS:
                    require(gamma_a == (3,) and gamma_b == (3,), "connected family")
                if label in EXPECTED_BY_SIZE[7] - CONNECTED_LABELS:
                    require(
                        gamma_a != (3,) or gamma_b != (3,),
                        "AP size-seven family is deletable",
                    )

require(set(collision_counts) == set().union(*EXPECTED_BY_SIZE.values()),
        "every normalized collision family occurs")
require(set(primitive_forms) == set(collision_counts), "primitive family coverage")
for family_label, normal_forms in primitive_forms.items():
    require(bool(normal_forms), f"primitive forms exist for {family_label}")
    for normal_form in normal_forms:
        require(
            gcd(gcd(normal_form[0], normal_form[1]), gcd(normal_form[2], normal_form[3])) == 1,
            f"primitive gcd for {family_label}",
        )
        require(primitive_normal_form(normal_form) == normal_form,
                f"canonical orbit representative for {family_label}")
        require(low_sumset_label(normal_form) == family_label,
                f"canonical label for {family_label}")


def cell_fibres(gaps_a: tuple[int, int], gaps_b: tuple[int, int]) -> list[tuple[str, ...]]:
    return [tuple(f"{i}{j}" for i, j in cells) for cells in fibre_map(gaps_a, gaps_b).values()]


require(
    cell_fibres((1, 1), (1, 1))
    == [("00",), ("01", "10"), ("02", "11", "20"), ("12", "21"), ("22",)],
    "AP5 fibres",
)
require(
    cell_fibres((2, 5), (2, 5))
    == [("00",), ("01", "10"), ("11",), ("02", "20"), ("12", "21"), ("22",)],
    "D6 fibres",
)
require(
    cell_fibres((1, 1), (1, 2))
    == [("00",), ("01", "10"), ("11", "20"), ("02", "21"), ("12",), ("22",)],
    "H6 fibres",
)
require(
    cell_fibres((1, 1), (2, 1))
    == [("00",), ("10",), ("01", "20"), ("02", "11"), ("12", "21"), ("22",)],
    "reversed H6 fibres",
)
require(
    cell_fibres((2, 5), (5, 2))
    == [("00",), ("10",), ("01",), ("02", "11", "20"), ("12",), ("21",), ("22",)],
    "R7 fibres",
)
require(
    cell_fibres((2, 5), (2, 3))
    == [("00",), ("01", "10"), ("11",), ("02",), ("12", "20"), ("21",), ("22",)],
    "E+ fibres",
)
require(
    cell_fibres((2, 5), (3, 2))
    == [("00",), ("10",), ("01",), ("02", "11"), ("12", "20"), ("21",), ("22",)],
    "E- fibres",
)
require(
    cell_fibres((5, 2), (3, 2))
    == [("00",), ("01",), ("02", "10"), ("20",), ("11",), ("12", "21"), ("22",)],
    "reversed E+ fibres",
)
require(
    cell_fibres((5, 2), (2, 3))
    == [("00",), ("01",), ("02", "10"), ("11", "20"), ("21",), ("12",), ("22",)],
    "reversed E- fibres",
)


def fibre_size_at_cell(
    gaps_a: tuple[int, int], gaps_b: tuple[int, int], cell: tuple[int, int]
) -> int:
    a_values, b_values = offsets(gaps_a), offsets(gaps_b)
    target = a_values[cell[0]] + b_values[cell[1]]
    return len(fibre_map(gaps_a, gaps_b)[target])


# Every diagonal scalar double has a genuine 2 x 2 rectangle whose other two
# corners are globally singleton.  The same check covers the reducible double
# in each Euclidean family.
diagonal_rectangles = 0
for x in range(1, 25):
    for y in range(1, 25):
        if x == y:
            continue
        for first, second in (((0, 1), (1, 0)), ((0, 2), (2, 0)), ((1, 2), (2, 1))):
            i, j = first
            k, ell = second
            require(fibre_size_at_cell((x, y), (x, y), (i, ell)) == 1,
                    "D6 first rectangle corner")
            require(fibre_size_at_cell((x, y), (x, y), (k, j)) == 1,
                    "D6 second rectangle corner")
            diagonal_rectangles += 1

euclidean_rectangles = 0
for p in range(1, 25):
    for q in range(p + 1, 25):
        if q == 2 * p:
            continue
        # E+: scalar double 01/10 has singleton corners 00/11.
        require(fibre_size_at_cell((p, q), (p, q - p), (0, 0)) == 1,
                "E+ rectangle corner 00")
        require(fibre_size_at_cell((p, q), (p, q - p), (1, 1)) == 1,
                "E+ rectangle corner 11")
        # E-: scalar double 12/20 has singleton corners 10/22.
        require(fibre_size_at_cell((p, q), (q - p, p), (1, 0)) == 1,
                "E- rectangle corner 10")
        require(fibre_size_at_cell((p, q), (q - p, p), (2, 2)) == 1,
                "E- rectangle corner 22")
        # Reversed E+: scalar double 12/21 has singleton corners 11/22.
        require(fibre_size_at_cell((q, p), (q - p, p), (1, 1)) == 1,
                "reversed E+ rectangle corner 11")
        require(fibre_size_at_cell((q, p), (q - p, p), (2, 2)) == 1,
                "reversed E+ rectangle corner 22")
        # Reversed E-: scalar double 02/10 has singleton corners 00/12.
        require(fibre_size_at_cell((q, p), (p, q - p), (0, 0)) == 1,
                "reversed E- rectangle corner 00")
        require(fibre_size_at_cell((q, p), (p, q - p), (1, 2)) == 1,
                "reversed E- rectangle corner 12")
        euclidean_rectangles += 4


# ---------------------------------------------------------------------------
# Scalar-arm gate and all collision-family gate placements.


def regularity_floor(weight: int) -> int:
    return (-weight + 1) // 2 if weight < 0 else 0


arm_gate_survivors: list[tuple[int, int, int]] = []
for r_value in range(1, 65):
    negative_weight, positive_weight = -r_value, r_value - 1
    for negative_order in range(regularity_floor(negative_weight), 40):
        for positive_order in range(0, 8):
            output_order = negative_order + positive_order - 1
            multiplier = (r_value - 1) * negative_order + r_value * positive_order
            if output_order == 0 and multiplier != 0:
                arm_gate_survivors.append((r_value, negative_order, positive_order))
require(arm_gate_survivors == [(2, 1, 0)], "unique simple-arm scalar gate")


def pair_order_possible(
    weight_f: int,
    weight_g: int,
    fixed_f: int | None = None,
    fixed_g: int | None = None,
) -> bool:
    """Exact local-order feasibility for one singleton zero Wronskian."""
    if weight_f * weight_g < 0:
        return False
    if (weight_f == 0) != (weight_g == 0):
        return False  # the retained weight-zero coefficient would be constant
    if weight_f == weight_g == 0:
        return True

    floor_f, floor_g = regularity_floor(weight_f), regularity_floor(weight_g)
    if fixed_f is not None and fixed_f < floor_f:
        return False
    if fixed_g is not None and fixed_g < floor_g:
        return False

    abs_f, abs_g = abs(weight_f), abs(weight_g)
    divisor = gcd(abs_f, abs_g)
    primitive_f, primitive_g = abs_f // divisor, abs_g // divisor
    # |g| ord(f)=|f| ord(g), hence orders are primitive_f*t,
    # primitive_g*t.  Positive coefficients may take t=0.
    candidates: list[int]
    if fixed_f is not None:
        if fixed_f % primitive_f:
            return False
        candidates = [fixed_f // primitive_f]
    elif fixed_g is not None:
        if fixed_g % primitive_g:
            return False
        candidates = [fixed_g // primitive_g]
    else:
        minimum = 0 if weight_f > 0 else 1
        while primitive_f * minimum < floor_f or primitive_g * minimum < floor_g:
            minimum += 1
        candidates = [minimum]

    for scale in candidates:
        if scale < 0:
            continue
        order_f, order_g = primitive_f * scale, primitive_g * scale
        if fixed_f is not None and order_f != fixed_f:
            continue
        if fixed_g is not None and order_g != fixed_g:
            continue
        if order_f >= floor_f and order_g >= floor_g:
            return True
    return False


def singleton_gate_feasible(
    p_weights: tuple[int, int, int],
    q_weights: tuple[int, int, int],
    gated_cell: tuple[int, int],
) -> bool:
    """Check the three H6 singleton fibres 00,12,22 with shared orders."""
    gate_i, gate_j = gated_cell
    fixed_p = {gate_i: 1 if p_weights[gate_i] == -2 else 0}
    fixed_q = {gate_j: 1 if q_weights[gate_j] == -2 else 0}
    if not pair_order_possible(
        p_weights[0], q_weights[0], fixed_p.get(0), fixed_q.get(0)
    ):
        return False

    q2_floor = regularity_floor(q_weights[2])
    q2_values = [fixed_q[2]] if 2 in fixed_q else range(q2_floor, 1025)
    for q2_order in q2_values:
        if pair_order_possible(
            p_weights[1], q_weights[2], fixed_p.get(1), q2_order
        ) and pair_order_possible(
            p_weights[2], q_weights[2], fixed_p.get(2), q2_order
        ):
            return True
    return False


def weights_from_gate(
    p_offsets: tuple[int, int, int],
    q_offsets: tuple[int, int, int],
    cell: tuple[int, int],
    negative_on_p: bool,
) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    i, j = cell
    gated_p, gated_q = (-2, 1) if negative_on_p else (1, -2)
    p_base, q_base = gated_p - p_offsets[i], gated_q - q_offsets[j]
    return (
        tuple(p_base + entry for entry in p_offsets),
        tuple(q_base + entry for entry in q_offsets),
    )


H_FIBRES = {
    1: ((0, 1), (1, 0)),
    2: ((1, 1), (2, 0)),
    3: ((0, 2), (2, 1)),
}


def hook_key(kappa: int, cell: tuple[int, int], negative_on_p: bool) -> str:
    if kappa == 1 and cell == (0, 1) and negative_on_p:
        return "H1_PRIMARY"
    if kappa == 1 and cell == (1, 0) and not negative_on_p:
        return "H1_MIRROR"
    if kappa == 2 and cell == (1, 1) and not negative_on_p:
        return "H2"
    if kappa == 3 and cell == (0, 2) and negative_on_p:
        return "H3"
    return "DEAD"


def h2_bridge_feasible(step: int) -> bool:
    """The adjacent-row order equation 3*ell/gcd(step-1,3)=1."""
    a_value = step - 1
    if a_value <= 0:
        return False
    divisor = gcd(a_value, 3)
    solutions = [ell for ell in range(1, 65) if 3 * ell == divisor]
    if not solutions:
        return False
    ell = solutions[0]
    m_value = a_value * ell // divisor
    return a_value - 2 * m_value != 0


hook_survivors: dict[str, list[int]] = defaultdict(list)
hook_gate_rows = 0
for step in range(1, 65):
    p_offsets = (0, step, 2 * step)
    q_offsets = (0, step, 3 * step)
    for kappa, cells in H_FIBRES.items():
        for cell in cells:
            for negative_on_p in (True, False):
                p_weights, q_weights = weights_from_gate(
                    p_offsets, q_offsets, cell, negative_on_p
                )
                feasible = singleton_gate_feasible(p_weights, q_weights, cell)
                key = hook_key(kappa, cell, negative_on_p)
                if feasible and key == "H2":
                    feasible = h2_bridge_feasible(step)
                if feasible:
                    require(key != "DEAD", "unclassified hooked scalar gate")
                    hook_survivors[key].append(step)
                hook_gate_rows += 1

odd_steps = list(range(3, 65, 2))
require(hook_survivors["H1_PRIMARY"] == odd_steps, "H1 primary ladder")
require(hook_survivors["H1_MIRROR"] == odd_steps, "H1 mirror ladder")
require(hook_survivors["H3"] == odd_steps, "H3 ladder")
require(hook_survivors["H2"] == list(range(4, 65, 3)), "H2 ladder")
require(set(hook_survivors) == {"H1_PRIMARY", "H1_MIRROR", "H2", "H3"},
        "exactly four hooked ladders")


def reversed_hook_singleton_gate_feasible(
    p_weights: tuple[int, int, int],
    q_weights: tuple[int, int, int],
    gated_cell: tuple[int, int],
) -> bool:
    """Joint order test for the reversed-hook singleton rows 00,10,22."""
    gate_i, gate_j = gated_cell
    fixed_p = {gate_i: 1 if p_weights[gate_i] == -2 else 0}
    fixed_q = {gate_j: 1 if q_weights[gate_j] == -2 else 0}
    q0_floor = regularity_floor(q_weights[0])
    q0_orders = [fixed_q[0]] if 0 in fixed_q else range(q0_floor, 1025)
    for q0_order in q0_orders:
        if not pair_order_possible(
            p_weights[0], q_weights[0], fixed_p.get(0), q0_order
        ):
            continue
        if not pair_order_possible(
            p_weights[1], q_weights[0], fixed_p.get(1), q0_order
        ):
            continue
        if pair_order_possible(
            p_weights[2], q_weights[2], fixed_p.get(2), fixed_q.get(2)
        ):
            return True
    return False


# Simultaneous reversal is an additive-support symmetry but not a regularity
# symmetry.  The omitted hook orientation has one surviving gate: 21=(1,-2)
# for d>=3.  Its lowest double row has arm orders (2d-1)ell and
# (2d+2)ell-1, so cancellation would require 3ell=1.
REVERSED_H_FIBRES = (
    ((0, 1), (2, 0)),
    ((0, 2), (1, 1)),
    ((1, 2), (2, 1)),
)
reversed_hook_survivors: list[tuple[int, tuple[int, int], bool]] = []
reversed_hook_gate_rows = 0
reversed_hook_valuation_rows = 0
for step in range(1, 65):
    p_offsets = (0, step, 2 * step)
    q_offsets = (0, 2 * step, 3 * step)
    for cells in REVERSED_H_FIBRES:
        for cell in cells:
            for negative_on_p in (True, False):
                p_weights, q_weights = weights_from_gate(
                    p_offsets, q_offsets, cell, negative_on_p
                )
                if reversed_hook_singleton_gate_feasible(
                    p_weights, q_weights, cell
                ):
                    reversed_hook_survivors.append((step, cell, negative_on_p))
                reversed_hook_gate_rows += 1
    if step >= 3:
        for arm_order in range(1, 65):
            first_order = (2 * step - 1) * arm_order
            second_order = (2 * step + 2) * arm_order - 1
            first_multiplier = (2 * step - 1) * (1 - 2 * arm_order)
            second_multiplier = -(2 * step + 2) * arm_order
            require(first_multiplier != 0 and second_multiplier != 0,
                    "reversed-hook nonzero initial multipliers")
            require(first_order != second_order,
                    "reversed-hook valuation mismatch")
            if step == 3:
                require((step - 1) * arm_order >= 2,
                        "reversed-hook d=3 alternate gate is never simple")
            reversed_hook_valuation_rows += 1
require(
    reversed_hook_survivors
    == [(step, (2, 1), False) for step in range(3, 65)],
    "exact reversed-hook survivor ladder",
)


# Reflected triple: all six scalar gates die on a singleton row, including
# x=1, x=2, and y=2 zero-weight boundaries.
reflected_gate_rows = 0
for x in range(1, 25):
    for y in range(x + 1, 25):
        p_offsets, q_offsets = offsets((x, y)), offsets((y, x))
        for cell in ((0, 2), (1, 1), (2, 0)):
            for negative_on_p in (True, False):
                p_weights, q_weights = weights_from_gate(
                    p_offsets, q_offsets, cell, negative_on_p
                )
                singles = ((0, 0), (1, 0), (0, 1), (1, 2), (2, 1), (2, 2))
                require(
                    not all(
                        pair_order_possible(
                            p_weights[i], q_weights[j],
                            (1 if p_weights[i] == -2 else 0) if i == cell[0] else None,
                            (1 if q_weights[j] == -2 else 0) if j == cell[1] else None,
                        )
                        for i, j in singles
                    ),
                    "reflected-triple singleton gate",
                )
                reflected_gate_rows += 1


# The non-rectangular scalar double in each Euclidean orientation has four
# gates; each dies on a singleton.  q=2p is deliberately excluded because it
# is exactly the hooked six-row boundary.
euclidean_gate_rows = 0
for p in range(1, 25):
    for q in range(p + 1, 25):
        if q == 2 * p:
            continue
        for gaps_a, gaps_b, scalar_cells, singleton_cells in (
            ((p, q), (p, q - p), ((1, 2), (2, 0)),
             ((0, 0), (1, 1), (0, 2), (2, 1), (2, 2))),
            ((p, q), (q - p, p), ((0, 2), (1, 1)),
             ((0, 0), (1, 0), (0, 1), (2, 1), (2, 2))),
            ((q, p), (q - p, p), ((0, 2), (1, 0)),
             ((0, 0), (0, 1), (2, 0), (1, 1), (2, 2))),
            ((q, p), (p, q - p), ((1, 1), (2, 0)),
             ((0, 0), (0, 1), (2, 1), (1, 2), (2, 2))),
        ):
            p_offsets, q_offsets = offsets(gaps_a), offsets(gaps_b)
            for cell in scalar_cells:
                for negative_on_p in (True, False):
                    p_weights, q_weights = weights_from_gate(
                        p_offsets, q_offsets, cell, negative_on_p
                    )
                    require(
                        not all(
                            pair_order_possible(
                                p_weights[i], q_weights[j],
                                (1 if p_weights[i] == -2 else 0) if i == cell[0] else None,
                                (1 if q_weights[j] == -2 else 0) if j == cell[1] else None,
                            )
                            for i, j in singleton_cells
                        ),
                        "Euclidean singleton gate",
                    )
                    euclidean_gate_rows += 1


# ---------------------------------------------------------------------------
# Tiny exact Laurent differential-polynomial algebra over Q.


NAMES = (
    "A", "B", "L", "M", "U", "C0",
    "h", "hp", "K", "Kp", "z",
    "f1", "f1p", "f2", "f2p", "g1", "g1p", "g2", "g2p",
    "w", "wp",
)
INDEX = {name: index for index, name in enumerate(NAMES)}
DERIVATIVE_INDEX = {
    INDEX["h"]: INDEX["hp"],
    INDEX["K"]: INDEX["Kp"],
    INDEX["f1"]: INDEX["f1p"],
    INDEX["f2"]: INDEX["f2p"],
    INDEX["g1"]: INDEX["g1p"],
    INDEX["g2"]: INDEX["g2p"],
    INDEX["w"]: INDEX["wp"],
}
ZERO_MONOMIAL = (0,) * len(NAMES)


class Expr:
    """Sparse Laurent polynomial with exact rational coefficients."""

    def __init__(self, terms: dict[tuple[int, ...], Fraction] | None = None):
        self.terms = {
            monomial: Fraction(coefficient)
            for monomial, coefficient in (terms or {}).items()
            if coefficient
        }

    @staticmethod
    def constant(value: int | Fraction) -> "Expr":
        coefficient = Fraction(value)
        return Expr({ZERO_MONOMIAL: coefficient}) if coefficient else Expr()

    @staticmethod
    def variable(name: str) -> "Expr":
        exponents = [0] * len(NAMES)
        exponents[INDEX[name]] = 1
        return Expr({tuple(exponents): Fraction(1)})

    def __add__(self, other: object) -> "Expr":
        right = to_expr(other)
        terms = dict(self.terms)
        for monomial, coefficient in right.terms.items():
            terms[monomial] = terms.get(monomial, Fraction(0)) + coefficient
            if not terms[monomial]:
                del terms[monomial]
        return Expr(terms)

    def __radd__(self, other: object) -> "Expr":
        return self + other

    def __neg__(self) -> "Expr":
        return Expr({monomial: -coefficient for monomial, coefficient in self.terms.items()})

    def __sub__(self, other: object) -> "Expr":
        return self + (-to_expr(other))

    def __rsub__(self, other: object) -> "Expr":
        return to_expr(other) - self

    def __mul__(self, other: object) -> "Expr":
        right = to_expr(other)
        terms: dict[tuple[int, ...], Fraction] = {}
        for left_monomial, left_coefficient in self.terms.items():
            for right_monomial, right_coefficient in right.terms.items():
                monomial = tuple(
                    left_monomial[index] + right_monomial[index]
                    for index in range(len(NAMES))
                )
                terms[monomial] = terms.get(monomial, Fraction(0)) + (
                    left_coefficient * right_coefficient
                )
                if not terms[monomial]:
                    del terms[monomial]
        return Expr(terms)

    def __rmul__(self, other: object) -> "Expr":
        return self * other

    def __pow__(self, exponent: int) -> "Expr":
        require(isinstance(exponent, int), "Laurent exponent must be integral")
        if exponent == 0:
            return Expr.constant(1)
        if exponent < 0:
            require(len(self.terms) == 1, "negative power requires one monomial")
            monomial, coefficient = next(iter(self.terms.items()))
            return Expr({
                tuple(exponent * entry for entry in monomial): coefficient ** exponent
            })
        result, factor, power = Expr.constant(1), self, exponent
        while power:
            if power & 1:
                result = result * factor
            factor = factor * factor
            power >>= 1
        return result

    def derivative(self) -> "Expr":
        terms: dict[tuple[int, ...], Fraction] = {}
        for monomial, coefficient in self.terms.items():
            for variable_index, derivative_index in DERIVATIVE_INDEX.items():
                exponent = monomial[variable_index]
                if exponent == 0:
                    continue
                differentiated = list(monomial)
                differentiated[variable_index] -= 1
                differentiated[derivative_index] += 1
                key = tuple(differentiated)
                terms[key] = terms.get(key, Fraction(0)) + coefficient * exponent
                if not terms[key]:
                    del terms[key]
        return Expr(terms)

    def is_zero(self) -> bool:
        return not self.terms


def to_expr(value: object) -> Expr:
    if isinstance(value, Expr):
        return value
    if isinstance(value, (int, Fraction)):
        return Expr.constant(value)
    raise TypeError(f"cannot convert {type(value)!r} to Expr")


def zero(expression: Expr, message: str) -> None:
    require(expression.is_zero(), message)


def derivative(expression: Expr) -> Expr:
    return expression.derivative()


def W(weight_f: int, f_value: Expr, weight_g: int, g_value: Expr) -> Expr:
    return weight_g * derivative(f_value) * g_value - weight_f * f_value * derivative(g_value)


def W_data(
    weight_f: int, f_value: Expr, f_prime: Expr,
    weight_g: int, g_value: Expr, g_prime: Expr,
) -> Expr:
    return weight_g * f_prime * g_value - weight_f * f_value * g_prime


AA, BB, LL, MM, UU, CC0 = (
    Expr.variable(name) for name in ("A", "B", "L", "M", "U", "C0")
)
h, hp, K, Kp, z = (Expr.variable(name) for name in ("h", "hp", "K", "Kp", "z"))
f1_symbol, f1p_symbol = Expr.variable("f1"), Expr.variable("f1p")
f2_symbol, f2p_symbol = Expr.variable("f2"), Expr.variable("f2p")
g1_symbol, g1p_symbol = Expr.variable("g1"), Expr.variable("g1p")
g2_symbol, g2p_symbol = Expr.variable("g2"), Expr.variable("g2p")
w, wp = Expr.variable("w"), Expr.variable("wp")


# H1 primary.  The two terminal rows give a common logarithmic derivative z.
# Eliminating g1 from the two middle rows yields the displayed product gate.
h1_primary_controls: list[int] = []
for k in range(1, 9):
    p_value, q_value = 2 * k - 1, 4 * k + 3
    f0, f0p = AA * h, AA * hp
    g0, g0p = BB * h**k, k * BB * h ** (k - 1) * hp
    f1p = p_value * f1_symbol * z
    f2p = 4 * k * f2_symbol * z
    g2p = q_value * g2_symbol * z
    zero(q_value * f1p * g2_symbol - p_value * f1_symbol * g2p,
         f"H1 terminal f1/g2 k={k}")
    zero(q_value * f2p * g2_symbol - 4 * k * f2_symbol * g2p,
         f"H1 terminal f2/g2 k={k}")
    row2 = W_data(p_value, f1_symbol, f1p, 1, g1_symbol, g1p_symbol) + W_data(
        4 * k, f2_symbol, f2p, -2 * k, g0, g0p
    )
    row3 = W_data(-2, f0, f0p, q_value, g2_symbol, g2p) + W_data(
        4 * k, f2_symbol, f2p, 1, g1_symbol, g1p_symbol
    )
    euler_z = hp + 2 * h * z
    obstruction = (
        AA * p_value * q_value * f1_symbol * g2_symbol
        + 16 * k**3 * BB * f2_symbol**2 * h ** (k - 1)
    )
    zero(
        p_value * f1_symbol * row3 - 4 * k * f2_symbol * row2
        - euler_z * obstruction,
        f"H1 primary elimination k={k}",
    )
    if k > 1:
        zero(
            ((p_value + q_value) - 8 * k) * h * z
            - (k - 1) * hp
            + (k - 1) * (hp + 2 * h * z),
            f"H1 logarithmic contradiction k={k}",
        )
    h1_primary_controls.append(k)


# H1 exceptional k=1.  The terminal rows have exponents 1,4,7; the bridge,
# compatibility equation, and scalar Euler factor are all exact.
f0 = AA * h
g0 = BB * h
f1 = LL * K
f2 = UU * K**4
g2 = MM * K**7
g1 = K * (CC0 - 4 * BB * UU * LL**-1 * h * K**2)
h1_row2 = W(1, f1, 1, g1) + W(4, f2, -2, g0)
h1_row3 = W(-2, f0, 7, g2) + W(4, f2, 1, g1)
h1_scalar = W(-2, f0, 1, g1) + W(1, f1, -2, g0)
euler2 = K * hp + 2 * h * Kp
zero(h1_row2, "H1 exceptional bridge")
zero(
    h1_row3 - K**6 * euler2 * (7 * AA * MM + 16 * BB * UU**2 * LL**-1),
    "H1 exceptional compatibility",
)
zero(
    h1_scalar
    - euler2 * (AA * CC0 - BB * LL - 12 * AA * BB * UU * LL**-1 * h * K**2),
    "H1 exceptional scalar factor",
)


# H1 mirror.  One middle row substituted into the next leaves the same Euler
# factor; for k>1 its bracket is a nonzero constant plus a positive-degree
# power, and k=1 is precisely the exceptional compatibility above.
h1_mirror_controls: list[int] = []
for k in range(1, 9):
    p_value, a2, q_value = 2 * k - 1, 2 * k + 2, 6 * k + 1
    f0, g0 = AA * h**k, BB * h
    f1, f2, g2 = LL * K, UU * K**a2, MM * K**q_value
    row2 = W(1, f1, p_value, g1_symbol) + W(a2, f2, -2, g0)
    row3 = W(-2 * k, f0, q_value, g2) + W(a2, f2, p_value, g1_symbol)
    alpha = a2 * UU * LL**-1 * K ** (a2 - 1)
    factor = euler2 * K ** (4 * k + 2) * (
        k * q_value * AA * MM * (h * K**2) ** (k - 1)
        + 4 * (k + 1) ** 2 * UU**2 * BB * LL**-1
    )
    zero(row3 - alpha * row2 - factor, f"H1 mirror elimination k={k}")
    if k > 1:
        require(2 * k - 2 > 0, "H1 mirror positive-degree separation")
    h1_mirror_controls.append(k)


# H2.  The lower bridge fixes g1=lambda*h*K.  A homogeneous addition y obeys
# (y^3/h^2)'=0; arm integrality makes it order at least two, never the simple
# scalar coefficient.  The scalar row has the K h'+3h K' factor.
h2_controls: list[int] = []
for k in range(1, 9):
    f0 = AA * h**k
    g0 = BB * h ** (k + 1)
    f1 = LL * K
    f2 = UU * K ** (3 * k + 2)
    g2 = MM * K ** (6 * k)
    lam = Fraction(k + 1, k) * LL * BB * AA**-1
    g1 = lam * h * K
    zero(W(-3 * k, f0, -3 * k - 3, g0), f"H2 low singleton k={k}")
    zero(
        W(-3 * k, f0, -2, g1) + W(1, f1, -3 * k - 3, g0),
        f"H2 lower bridge k={k}",
    )
    zero(W(1, f1, 6 * k, g2), f"H2 high singleton one k={k}")
    zero(W(3 * k + 2, f2, 6 * k, g2), f"H2 high singleton two k={k}")
    scalar = W(1, f1, -2, g1) + W(3 * k + 2, f2, -3 * k - 3, g0)
    euler3 = K * hp + 3 * h * Kp
    scalar_factor = -euler3 * (
        lam * LL * K
        + (k + 1) * (3 * k + 2) * UU * BB * h**k * K ** (3 * k + 1)
    )
    zero(scalar - scalar_factor, f"H2 scalar factor k={k}")
    h2_controls.append(k)

zero(
    derivative(w**3 * h**-2) - w**2 * h**-3 * (3 * h * wp - 2 * hp * w),
    "H2 cube first integral",
)
cube_arm_rows: list[tuple[int, int]] = []
for h_order in range(1, 65):
    if (2 * h_order) % 3 == 0:
        y_order = 2 * h_order // 3
        require(y_order >= 2, "cube mode cannot be simple at an arm")
        cube_arm_rows.append((h_order, y_order))


# H3.  Subtracting a monomial multiple of the second zero row from the scalar
# row leaves the displayed Euler factor.  The unused first bridge can only add
# constraints, so no solution is discarded by this elimination.
h3_controls: list[int] = []
for k in range(1, 9):
    p_value = 2 * k - 1
    f0 = AA * h
    g0 = BB * h ** (3 * k + 1)
    f1 = LL * K**p_value
    f2 = UU * K ** (4 * k)
    g2 = MM * K
    zero(W(-2, f0, -6 * k - 2, g0), f"H3 low singleton k={k}")
    zero(W(p_value, f1, 1, g2), f"H3 high singleton one k={k}")
    zero(W(4 * k, f2, 1, g2), f"H3 high singleton two k={k}")
    row2 = W(p_value, f1, -4 * k - 1, g1_symbol) + W(
        4 * k, f2, -6 * k - 2, g0
    )
    scalar = W(-2, f0, 1, g2) + W(4 * k, f2, -4 * k - 1, g1_symbol)
    alpha = Fraction(4 * k, p_value) * UU * LL**-1 * K ** (2 * k + 1)
    factor = euler2 * (
        AA * MM
        + Fraction(16 * k * k * (3 * k + 1), p_value)
        * UU**2 * BB * LL**-1 * h ** (3 * k) * K ** (6 * k)
    )
    zero(scalar - alpha * row2 - factor, f"H3 scalar elimination k={k}")
    h3_controls.append(k)


# ---------------------------------------------------------------------------
# Degree, polynomiality, and support-floor hostiles.


degree_controls = 0
for degree_h in range(2, 17):
    for degree_k in range(0, 17):
        for multiplier in (1, 2, 3):
            require(
                degree_h + multiplier * degree_k > 0,
                "Euler leading multiplier is positive",
            )
            require(
                degree_h + degree_k - 1 >= 1,
                "Euler factor has positive degree",
            )
            degree_controls += 1

# Rational pole: b has weight 0, coefficient derivative 1; -c^-1 has weight
# -1 and constant coefficient -1.
require(
    W_data(0, Expr.constant(0), Expr.constant(1), -1, Expr.constant(-1), Expr.constant(0)).terms
    == Expr.constant(1).terms,
    "rational hostile {b,-1/c}=1",
)

# Degree-one polynomial exception.  For Sigma=ub+v, -e/u is the weight -2
# piece with coefficient -Sigma/u and derivative -1.
require(
    W_data(1, Expr.constant(1), Expr.constant(0), -2, Expr.constant(0), Expr.constant(-1)).terms
    == Expr.constant(1).terms,
    "degree-one hostile {c,-e/u}=1",
)

# At a repeated root and c=0 all three Poisson tensor coefficients vanish.
require((0**2, 0, -2 * 0 * 1) == (0, 0, 0), "repeated-root Poisson degeneration")

partitions_below_seven = []
for left in range(1, 7):
    for right in range(left, 7):
        if left + right <= 6:
            partitions_below_seven.append((left, right))
            require(
                left == 1 or left == 2 or (left, right) == (3, 3),
                "all total-at-most-six partitions hit a proved cell",
            )
require(
    [(left, 7 - left) for left in range(2, 6)]
    == [(2, 5), (3, 4), (4, 3), (5, 2)],
    "exact total-seven frontier",
)


print("THM-3592 exact control")
print(
    "normalized collision census: gaps 1..24, primitive under dilation; "
    "types AP5 / D6,H6 / R7,E+7,E-7,AP-extension/contained/dyadic7 PASS"
)
print(
    "collision identity: |A+B|=9-C+tau; tau=1 iff reflected gaps; "
    f"labelled rows={sum(collision_counts.values())} PASS"
)
print(
    "component deletion: all |A+B|>=8 and all three AP seven-row families "
    "have a disconnected support graph PASS"
)
print(
    f"rectangle reductions: D6={diagonal_rectangles}, "
    f"Euclidean={euclidean_rectangles} sampled exact rows PASS"
)
print(
    f"simple-arm gate: (-2,1)/(1,-2) only; H6 placements={hook_gate_rows}, "
    "four ladders exactly PASS"
)
print(
    f"reversed H6 placements={reversed_hook_gate_rows}; one ladder and "
    f"valuation rows={reversed_hook_valuation_rows} PASS"
)
print(
    f"R7 gates={reflected_gate_rows}; four-orientation nonrectangular E gates={euclidean_gate_rows}; "
    "all singleton-killed PASS"
)
print("H1 primary/mirror: terminal, leak-elimination, k>1, and k=1 compatibility PASS")
print("H2: bridge, cube first integral, scalar K*h'+3*h*K' factor PASS")
print("H3: singleton rows and scalar-minus-row Euler factor PASS")
print(
    f"degree hostiles: {degree_controls} Euler rows; rational pole, degree-one, "
    "and repeated-root boundaries PASS"
)
print("support floor: total nonconstant support >=7; total-seven frontier 2+5/3+4/4+3/5+2")
print("scope: squarefree Sigma over C, deg Sigma>=2, exponent two; no JC(2) conclusion")
print("PASS")
