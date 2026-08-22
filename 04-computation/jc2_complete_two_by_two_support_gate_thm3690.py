#!/usr/bin/env python3
"""Exact all-activity closure of the normalized 2x2 sparse Keller chart.

This companion enumerates every activity mask and every possible equality
partition of the surviving Jacobian contributions.  It first performs an
affine-linear exit pass, then uses exact Groebner saturation for the residual
parallel/axis strata.  The final sixteen algebraic packets are printed and
checked directly against the nonnegative-integral exponent requirement.
"""

from __future__ import annotations

import ast
import hashlib
from collections import Counter
from itertools import combinations
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def no_singleton_partitions(items: tuple[int, ...]):
    """Enumerate each set partition once, including the empty partition."""

    blocks: list[list[int]] = []

    def visit(index: int):
        if index == len(items):
            if all(len(block) >= 2 for block in blocks):
                yield tuple(tuple(block) for block in blocks)
            return
        value = items[index]
        for block_index in range(len(blocks)):
            blocks[block_index].append(value)
            yield from visit(index + 1)
            blocks[block_index].pop()
        blocks.append([value])
        yield from visit(index + 1)
        blocks.pop()

    yield from visit(0)


def bounded_no_singleton_census(max_total_degree: int) -> tuple[int, int, int]:
    """Independent direct bucket census, with no symbolic partition solver."""

    monomials = [
        (i, j)
        for i in range(max_total_degree + 1)
        for j in range(max_total_degree + 1 - i)
        if i + j >= 2
    ]
    support_pairs = list(combinations(monomials, 2))
    survivors = 0

    for left in support_pairs:
        for right in support_pairs:
            multiplicities: dict[tuple[int, int], int] = {}

            def add(bucket: tuple[int, int]) -> None:
                multiplicities[bucket] = multiplicities.get(bucket, 0) + 1

            for alpha, beta in left:
                if alpha != 0:
                    add((alpha - 1, beta))
            for gamma, delta in right:
                if delta != 0:
                    add((gamma, delta - 1))
            for alpha, beta in left:
                for gamma, delta in right:
                    if alpha * delta - beta * gamma != 0:
                        add((alpha + gamma - 1, beta + delta - 1))

            # An empty debt dictionary would mean Jac(P,Q)=1 and therefore is
            # a survivor too; `all` gives that vacuous case the correct value.
            if all(value >= 2 for value in multiplicities.values()):
                survivors += 1

    return len(monomials), len(support_pairs) ** 2, survivors


def main() -> None:
    ax, ay, cx, cy, bx, by, dx, dy = variables = sp.symbols(
        "ax ay cx cy bx by dx dy"
    )
    sat = sp.symbols("sat")

    # Label order: P_A,P_C,Q_B,Q_D,J_AB,J_AD,J_CB,J_CD.
    buckets = (
        (ax - 1, ay),
        (cx - 1, cy),
        (bx, by - 1),
        (dx, dy - 1),
        (ax + bx - 1, ay + by - 1),
        (ax + dx - 1, ay + dy - 1),
        (cx + bx - 1, cy + by - 1),
        (cx + dx - 1, cy + dy - 1),
    )
    divergence_coordinates = (ax, cx, by, dy)
    cross_determinants = (
        ax * by - ay * bx,
        ax * dy - ay * dx,
        cx * by - cy * bx,
        cx * dy - cy * dx,
    )

    linear_exits: Counter[str] = Counter()
    residual = []
    activity_systems = 0

    # A bit is one precisely when the corresponding labelled contribution is
    # declared active.  Inactive divergence labels add a coordinate-zero
    # equation.  Inactive cross labels are imposed at the saturation stage.
    for divergence_mask in range(16):
        for cross_mask in range(16):
            active_labels = tuple(
                [i for i in range(4) if divergence_mask & (1 << i)]
                + [4 + j for j in range(4) if cross_mask & (1 << j)]
            )

            for partition in no_singleton_partitions(active_labels):
                activity_systems += 1
                equations = [
                    divergence_coordinates[i]
                    for i in range(4)
                    if not divergence_mask & (1 << i)
                ]
                for block in partition:
                    base = buckets[block[0]]
                    for label in block[1:]:
                        equations.extend(
                            (
                                buckets[label][0] - base[0],
                                buckets[label][1] - base[1],
                            )
                        )

                matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
                solution_set = sp.linsolve((matrix, rhs), variables)
                if solution_set is sp.EmptySet:
                    linear_exits["inconsistent"] += 1
                    continue

                solution = tuple(next(iter(solution_set)))
                substitution = dict(zip(variables, solution))
                vectors = (
                    solution[0:2],
                    solution[2:4],
                    solution[4:6],
                    solution[6:8],
                )

                reason = None
                if vectors[0] == vectors[1]:
                    reason = "equal_P_supports"
                elif vectors[2] == vectors[3]:
                    reason = "equal_Q_supports"
                else:
                    for name, vector in zip(
                        ("A_low", "C_low", "B_low", "D_low"), vectors
                    ):
                        total = sp.expand(vector[0] + vector[1])
                        if total.is_number and total < 2:
                            reason = name
                            break

                if reason is None:
                    active_divergences = (
                        vectors[0][0], vectors[1][0], vectors[2][1], vectors[3][1]
                    )
                    if any(
                        divergence_mask & (1 << i) and sp.expand(value) == 0
                        for i, value in enumerate(active_divergences)
                    ):
                        reason = "active_divergence_vanishes"

                specialized_determinants = tuple(
                    sp.factor(det.subs(substitution)) for det in cross_determinants
                )
                if reason is None:
                    for index, determinant in enumerate(specialized_determinants):
                        active = bool(cross_mask & (1 << index))
                        if active and determinant == 0:
                            reason = "active_cross_vanishes"
                            break
                        if not active and determinant.is_number and determinant != 0:
                            reason = "inactive_cross_nonzero"
                            break

                if reason is not None:
                    linear_exits[reason] += 1
                    continue

                residual.append(
                    (
                        divergence_mask,
                        cross_mask,
                        partition,
                        solution,
                        specialized_determinants,
                    )
                )

    expected_linear_exits = Counter(
        {
            "inconsistent": 1698,
            "equal_P_supports": 1805,
            "equal_Q_supports": 284,
            "A_low": 128,
            "C_low": 96,
            "B_low": 16,
            "D_low": 16,
        }
    )
    gate(activity_systems == 4140, "wrong activity/partition system count")
    gate(linear_exits == expected_linear_exits, "wrong affine exit ledger")
    gate(len(residual) == 97, "wrong residual-system count")

    saturation_attempts = 0
    saturation_skips = 0
    raw_packets = []
    nonunit_basis_shapes: Counter[tuple[str, ...]] = Counter()

    for divergence_mask, cross_mask, partition, solution, determinants in residual:
        parameters = sorted(
            set().union(*(entry.free_symbols for entry in solution)), key=str
        )
        vectors = (
            solution[0:2],
            solution[2:4],
            solution[4:6],
            solution[6:8],
        )
        active_divergences = (
            vectors[0][0], vectors[1][0], vectors[2][1], vectors[3][1]
        )
        inactive_determinants = [
            determinants[index]
            for index in range(4)
            if not cross_mask & (1 << index) and determinants[index] != 0
        ]
        active_determinants = [
            determinants[index]
            for index in range(4)
            if cross_mask & (1 << index)
        ]

        # Saturating by total*(total-1) is lossless for nonnegative integral
        # exponent vectors of total degree >=2.  Distinct supports are the
        # union of two principal opens on each side, so all four witness pairs
        # are tested separately below.
        nonzero_product = sp.prod(
            value
            for index, value in enumerate(active_divergences)
            if divergence_mask & (1 << index)
        ) * sp.prod(active_determinants)
        for vector in vectors:
            total = sp.expand(vector[0] + vector[1])
            nonzero_product *= total * (total - 1)

        left_witnesses = (
            sp.expand(vectors[0][0] - vectors[1][0]),
            sp.expand(vectors[0][1] - vectors[1][1]),
        )
        right_witnesses = (
            sp.expand(vectors[2][0] - vectors[3][0]),
            sp.expand(vectors[2][1] - vectors[3][1]),
        )

        for left_index, left_witness in enumerate(left_witnesses):
            for right_index, right_witness in enumerate(right_witnesses):
                saturation_factor = sp.factor(
                    nonzero_product * left_witness * right_witness
                )
                if saturation_factor == 0:
                    saturation_skips += 1
                    continue

                saturation_attempts += 1
                basis = sp.groebner(
                    inactive_determinants + [sat * saturation_factor - 1],
                    *(parameters + [sat]),
                    order="grevlex",
                    method="f5b",
                )
                if basis.contains(sp.Integer(1)):
                    continue

                basis_expressions = [poly.as_expr() for poly in basis.polys]
                generators = parameters + [sat]
                gate(basis.is_zero_dimensional, "nonunit saturation is not zero-dimensional")
                gate(
                    len(basis_expressions) == len(generators)
                    and all(
                        sp.Poly(expression, *generators).total_degree() == 1
                        for expression in basis_expressions
                    ),
                    "nonunit saturation did not reduce to a square affine-linear basis",
                )
                nonunit_basis_shapes[
                    tuple(str(expression) for expression in basis_expressions)
                ] += 1

                # Every surviving reduced Groebner basis is a square affine-
                # linear system.  Solve it by linsolve, not by a heuristic
                # nonlinear root extractor, and demand its unique exact point.
                point_set = sp.linsolve(basis_expressions, generators)
                gate(point_set != sp.EmptySet, "nonunit saturation has no point")
                points = list(point_set)
                gate(len(points) == 1, "nonunit saturation is not a singleton")
                point = points[0]
                gate(
                    all(not entry.free_symbols for entry in point),
                    "nonunit saturation point retains a free parameter",
                )
                parameter_solution = dict(zip(generators, point))
                defining_equations = inactive_determinants + [
                    sat * saturation_factor - 1
                ]
                gate(
                    all(
                        sp.factor(expression.subs(parameter_solution)) == 0
                        for expression in defining_equations
                    ),
                    "extracted point does not satisfy its saturation ideal",
                )
                packet = tuple(
                    sp.factor(entry.subs(parameter_solution)) for entry in solution
                )
                gate(
                    all(not entry.free_symbols for entry in packet),
                    "saturation point did not specialize all exponents",
                )
                packet_substitution = dict(zip(variables, packet))
                packet_divergences = (
                    packet[0], packet[2], packet[5], packet[7]
                )
                packet_determinants = tuple(
                    sp.factor(det.subs(packet_substitution))
                    for det in cross_determinants
                )
                gate(
                    all(
                        (value != 0) == bool(divergence_mask & (1 << index))
                        for index, value in enumerate(packet_divergences)
                    ),
                    "packet divergence activity does not match its mask",
                )
                gate(
                    all(
                        (value != 0) == bool(cross_mask & (1 << index))
                        for index, value in enumerate(packet_determinants)
                    ),
                    "packet cross activity does not match its mask",
                )
                packet_vectors = (
                    packet[0:2], packet[2:4], packet[4:6], packet[6:8]
                )
                gate(
                    all(sum(vector) not in (0, 1) for vector in packet_vectors),
                    "packet violates the saturated nonlinear-total condition",
                )
                gate(
                    packet_vectors[0][left_index]
                    != packet_vectors[1][left_index]
                    and packet_vectors[2][right_index]
                    != packet_vectors[3][right_index],
                    "packet violates its selected support-distinctness witnesses",
                )
                packet_buckets = tuple(
                    tuple(sp.factor(entry.subs(packet_substitution)) for entry in bucket)
                    for bucket in buckets
                )
                gate(
                    all(
                        all(
                            packet_buckets[label] == packet_buckets[block[0]]
                            for label in block[1:]
                        )
                        for block in partition
                    ),
                    "packet violates its bucket-equality partition",
                )
                raw_packets.append(
                    (
                        divergence_mask,
                        cross_mask,
                        partition,
                        left_index,
                        right_index,
                        packet,
                    )
                )

    packet_multiplicities = Counter(row[-1] for row in raw_packets)
    gate(saturation_attempts == 369, "wrong saturation-attempt count")
    gate(saturation_skips == 19, "wrong identically-empty open count")
    gate(len(raw_packets) == 64, "wrong raw algebraic packet count")
    gate(sum(nonunit_basis_shapes.values()) == 64, "wrong nonunit-basis count")
    gate(len(nonunit_basis_shapes) == 48, "wrong affine-basis-shape count")
    gate(len(packet_multiplicities) == 16, "wrong unique algebraic packet count")
    gate(
        set(packet_multiplicities.values()) == {4},
        "distinctness-witness multiplicity is not uniformly four",
    )

    invalidity_rows = []
    for packet in sorted(packet_multiplicities, key=str):
        invalid_coordinates = []
        for index, entry in enumerate(packet):
            gate(entry.is_Rational, "final exponent coordinate is not rational")
            if entry.q != 1:
                invalid_coordinates.append(f"{index}:noninteger={entry}")
            elif entry < 0:
                invalid_coordinates.append(f"{index}:negative={entry}")
        gate(
            bool(invalid_coordinates),
            f"algebraic packet survived in nonnegative integer exponents: {packet}",
        )
        invalidity_rows.append(
            f"{packet}::" + ";".join(invalid_coordinates)
        )

    # Independent bounded representation: direct bucket multiplicities, no
    # activity masks, partitions, linsolve, Groebner bases, or saturation.
    monomial_count, bounded_pairs, bounded_survivors = bounded_no_singleton_census(8)
    gate(monomial_count == 42, "wrong bounded monomial count")
    gate(bounded_pairs == 741321, "wrong bounded support-pair universe")
    gate(bounded_survivors == 0, "bounded direct census found a survivor")

    # Positive tame support-drop control.
    x, y, alpha, delta = sp.symbols("x y alpha delta")
    P_control = x + alpha * y**2
    Q_control = y + delta * P_control**2
    J_control = sp.expand(
        sp.diff(P_control, x) * sp.diff(Q_control, y)
        - sp.diff(P_control, y) * sp.diff(Q_control, x)
    )
    gate(J_control == 1, "tame support-drop control failed")

    # Noninjective boundary hostile: collision is exact, but the Jacobian debt
    # remains visible.  This keeps nonproperness separate from Keller support.
    a, b, d = sp.symbols("a b d")
    P_hostile = x + a * x**2 - x**3
    Q_hostile = y + b * x**2 + d * x**4
    gate(
        sp.expand(P_hostile.subs({x: 1, y: 0}) - P_hostile.subs({x: -1, y: 0}))
        == 0,
        "hostile P collision failed",
    )
    gate(
        sp.expand(Q_hostile.subs({x: 1, y: 0}) - Q_hostile.subs({x: -1, y: 0}))
        == 0,
        "hostile Q collision failed",
    )
    J_hostile = sp.expand(
        sp.diff(P_hostile, x) * sp.diff(Q_hostile, y)
        - sp.diff(P_hostile, y) * sp.diff(Q_hostile, x)
    )
    gate(J_hostile == 1 + 2 * a * x - 3 * x**2, "hostile Jacobian debt changed")

    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    gate(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "inactive Python assert found",
    )

    semantic_rows = [
        f"linear_exits={sorted(linear_exits.items())}",
        f"residual={len(residual)}",
        f"saturation={saturation_attempts},{saturation_skips}",
        *(
            f"basis={count}:{shape}"
            for shape, count in sorted(nonunit_basis_shapes.items())
        ),
        *invalidity_rows,
    ]
    semantic_digest = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("THM-3690 exact companion -- complete normalized 2x2 sparse-support closure")
    print(f"activity_partition_systems={activity_systems}")
    print("linear_exits=" + ",".join(f"{key}:{linear_exits[key]}" for key in sorted(linear_exits)))
    print(f"residual_systems={len(residual)}")
    print(f"saturation_attempts={saturation_attempts}")
    print(f"saturation_skips={saturation_skips}")
    print(f"nonunit_affine_bases={sum(nonunit_basis_shapes.values())}")
    print(f"unique_affine_basis_shapes={len(nonunit_basis_shapes)}")
    print(f"raw_algebraic_packets={len(raw_packets)}")
    print(f"unique_algebraic_packets={len(packet_multiplicities)}")
    print("algebraic_packets_begin")
    for row in invalidity_rows:
        print(row)
    print("algebraic_packets_end")
    print(
        f"bounded_degree8=monomials:{monomial_count},support_pair_pairs:{bounded_pairs},survivors:{bounded_survivors}"
    )
    print(f"semantic_sha256={semantic_digest}")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
