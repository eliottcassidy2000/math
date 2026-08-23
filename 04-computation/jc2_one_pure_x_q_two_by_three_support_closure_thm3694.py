#!/usr/bin/env python3
"""Exact companion for the one-pure-x-Q normalized 2x3 closure.

The finite equality-partition core is followed by exact Groebner saturation
of every determinant, nonlinear-support, and distinctness condition.  Every
nonunit stratum forces the pure-x exponent to be 1/2 or 3/2.
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


def determinant(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def bounded_axis_census(max_total_degree: int) -> tuple[int, int, int]:
    """Direct census with two active P seeds and exactly one pure-x Q term."""

    monomials = [
        (x_degree, total - x_degree)
        for total in range(2, max_total_degree + 1)
        for x_degree in range(total + 1)
    ]
    left_monomials = [monomial for monomial in monomials if monomial[0] > 0]
    active_right = [monomial for monomial in monomials if monomial[1] > 0]
    pure_x_right = [monomial for monomial in monomials if monomial[1] == 0]
    universe = 0
    survivors = 0

    for left in combinations(left_monomials, 2):
        for active_pair in combinations(active_right, 2):
            for pure_x in pure_x_right:
                right = active_pair + (pure_x,)
                universe += 1
                buckets: Counter[tuple[int, int]] = Counter()
                for px, py in left:
                    buckets[(px - 1, py)] += 1
                for qx, qy in active_pair:
                    buckets[(qx, qy - 1)] += 1
                for p in left:
                    for q in right:
                        if determinant(p, q) != 0:
                            buckets[(p[0] + q[0] - 1, p[1] + q[1] - 1)] += 1
                if all(multiplicity >= 2 for multiplicity in buckets.values()):
                    survivors += 1
    return len(monomials), universe, survivors


def no_singleton_partitions(items: tuple[int, ...]):
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


def main() -> None:
    p1x, p1y, p2x, p2y, q1x, q1y, q2x, q2y, q3x, q3y = variables = sp.symbols(
        "p1x p1y p2x p2y q1x q1y q2x q2y q3x q3y"
    )
    sat = sp.symbols("sat")
    left = ((p1x, p1y), (p2x, p2y))
    right = ((q1x, q1y), (q2x, q2y), (q3x, q3y))
    buckets = (
        (p1x - 1, p1y),
        (p2x - 1, p2y),
        (q1x, q1y - 1),
        (q2x, q2y - 1),
        (q3x, q3y - 1),
        *((p[0] + q[0] - 1, p[1] + q[1] - 1) for p in left for q in right),
    )
    determinants = tuple(
        p[0] * q[1] - p[1] * q[0] for p in left for q in right
    )

    # Up to relabelling Q supports: q3 is the unique pure-x nonlinear term.
    active_divergence_labels = (0, 1, 2, 3)
    linear_exits: Counter[str] = Counter()
    residual = []
    system_count = 0

    for cross_mask in range(64):
        active_labels = active_divergence_labels + tuple(
            5 + index for index in range(6) if cross_mask & (1 << index)
        )
        for partition in no_singleton_partitions(active_labels):
            system_count += 1
            equations = [q3y]
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
            if solution_set == sp.EmptySet:
                linear_exits["inconsistent"] += 1
                continue
            solution = tuple(next(iter(solution_set)))
            substitution = dict(zip(variables, solution))
            vectors = (
                solution[0:2],
                solution[2:4],
                solution[4:6],
                solution[6:8],
                solution[8:10],
            )
            reason = None
            if vectors[0] == vectors[1]:
                reason = "equal_P"
            elif len({vectors[2], vectors[3], vectors[4]}) < 3:
                reason = "equal_Q"
            else:
                for name, vector in zip(("P1", "P2", "Q1", "Q2", "Q3"), vectors):
                    total = sp.expand(sum(vector))
                    if total.is_number and total < 2:
                        reason = f"{name}_low"
                        break
            active_divergences = (
                vectors[0][0], vectors[1][0], vectors[2][1], vectors[3][1]
            )
            if reason is None and any(value == 0 for value in active_divergences):
                reason = "active_divergence_zero"

            specialized_determinants = tuple(
                sp.factor(value.subs(substitution)) for value in determinants
            )
            if reason is None:
                for index, value in enumerate(specialized_determinants):
                    if cross_mask & (1 << index) and value == 0:
                        reason = "active_cross_zero"
                        break
                    if (
                        not cross_mask & (1 << index)
                        and value.is_number
                        and value != 0
                    ):
                        reason = "inactive_cross_nonzero"
                        break
            if reason is not None:
                linear_exits[reason] += 1
                continue
            residual.append((cross_mask, partition, solution, specialized_determinants))

    expected_linear_exits = Counter(
        {
            "inconsistent": 42_447,
            "equal_P": 7_782,
            "equal_Q": 1_813,
            "P1_low": 226,
            "P2_low": 226,
            "Q1_low": 108,
            "Q2_low": 108,
        }
    )
    gate(system_count == 52_922, "wrong activity/partition system count")
    gate(linear_exits == expected_linear_exits, "wrong affine exit ledger")
    gate(sum(linear_exits.values()) == 52_710, "wrong affine exit total")
    gate(len(residual) == 212, "wrong residual-system count")

    saturation_attempts = 0
    saturation_skips = 0
    unit_ideals = 0
    groebner_cache: dict[tuple[tuple[str, ...], tuple[str, ...], str], sp.GroebnerBasis] = {}
    groebner_cache_hits = 0
    nonunit_profiles: Counter[tuple[int, int, tuple[int, ...]]] = Counter()
    nonunit_dimensions: Counter[str] = Counter()
    nonunit_basis_shapes: Counter[tuple[str, ...]] = Counter()
    obstruction = sp.expand((2 * q3x - 1) * (2 * q3x - 3))

    for cross_mask, partition, solution, specialized_determinants in residual:
        parameters = sorted(
            set().union(*(entry.free_symbols for entry in solution)), key=str
        )
        vectors = (
            solution[0:2],
            solution[2:4],
            solution[4:6],
            solution[6:8],
            solution[8:10],
        )
        inactive_equations = [
            specialized_determinants[index]
            for index in range(6)
            if not cross_mask & (1 << index)
            and specialized_determinants[index] != 0
        ]
        active_factors = [
            specialized_determinants[index]
            for index in range(6)
            if cross_mask & (1 << index)
        ]
        active_factors.extend(
            (vectors[0][0], vectors[1][0], vectors[2][1], vectors[3][1])
        )
        for vector in vectors:
            total = sp.expand(sum(vector))
            active_factors.append(total * (total - 1))

        left_witnesses = (
            sp.expand(vectors[0][0] - vectors[1][0]),
            sp.expand(vectors[0][1] - vectors[1][1]),
        )
        # Q1 and Q2 are automatically distinct from pure-x Q3 because their
        # active y coordinates are saturated.  Only Q1!=Q2 needs a witness.
        right_witnesses = (
            sp.expand(vectors[2][0] - vectors[3][0]),
            sp.expand(vectors[2][1] - vectors[3][1]),
        )

        base_product = sp.prod(active_factors)
        for left_witness in left_witnesses:
            for right_witness in right_witnesses:
                saturation_factor = sp.factor(
                    base_product * left_witness * right_witness
                )
                if saturation_factor == 0:
                    saturation_skips += 1
                    continue
                saturation_attempts += 1
                normalized_inactive = tuple(
                    sorted({str(sp.factor(value)) for value in inactive_equations})
                )
                cache_key = (
                    tuple(str(parameter) for parameter in parameters),
                    normalized_inactive,
                    str(saturation_factor),
                )
                if cache_key in groebner_cache:
                    groebner_cache_hits += 1
                    basis = groebner_cache[cache_key]
                else:
                    basis = sp.groebner(
                        [sp.sympify(value) for value in normalized_inactive]
                        + [sat * saturation_factor - 1],
                        *(parameters + [sat]),
                        order="grevlex",
                        method="f5b",
                    )
                    groebner_cache[cache_key] = basis
                if basis.contains(sp.Integer(1)):
                    unit_ideals += 1
                    continue
                degrees = tuple(
                    sp.Poly(poly.as_expr(), *(parameters + [sat])).total_degree()
                    for poly in basis.polys
                )
                profile = (len(parameters), len(basis.polys), degrees)
                nonunit_profiles[profile] += 1
                nonunit_dimensions[
                    "zero_dimensional" if basis.is_zero_dimensional else "positive_dimensional"
                ] += 1
                expressions = tuple(poly.as_expr() for poly in basis.polys)
                nonunit_basis_shapes[tuple(str(value) for value in expressions)] += 1
                specialized_obstruction = sp.expand(
                    obstruction.subs(dict(zip(variables, solution)))
                )
                _, remainder = basis.reduce(specialized_obstruction)
                gate(
                    sp.factor(remainder) == 0,
                    "nonunit stratum does not force the half-integral obstruction",
                )

    expected_profiles = Counter(
        {
            (1, 2, (1, 1)): 32,
            (3, 3, (18, 1, 1)): 16,
        }
    )
    gate(saturation_attempts == 848, "wrong saturation-attempt count")
    gate(saturation_skips == 0, "unexpected identically empty witness chart")
    gate(unit_ideals == 800, "wrong unit-ideal count")
    gate(sum(nonunit_profiles.values()) == 48, "wrong nonunit-stratum count")
    gate(nonunit_profiles == expected_profiles, "wrong nonunit basis profiles")
    gate(
        nonunit_dimensions
        == Counter({"zero_dimensional": 32, "positive_dimensional": 16}),
        "wrong nonunit dimension ledger",
    )
    gate(
        sp.solve(obstruction, q3x) == [sp.Rational(1, 2), sp.Rational(3, 2)],
        "wrong roots of the terminal exponent obstruction",
    )

    bounded_monomials, bounded_universe, bounded_survivors = bounded_axis_census(6)
    gate(bounded_monomials == 25, "wrong bounded monomial count")
    gate(bounded_universe == 180_500, "wrong bounded axis universe")
    gate(bounded_survivors == 0, "bounded axis census found a survivor")

    # A canonical near-survivor: every bucket except y^2 can be paired, but
    # its coefficient 2ab is forced nonzero in characteristic zero.
    x, y, a, b, c, d, e = sp.symbols("x y a b c d e")
    P_hostile = x + a * x * y + c * x**3
    Q_hostile = y + b * y**2 + d * x**2 * y + e * x**4
    hostile_debt = sp.Poly(
        sp.expand(
            sp.diff(P_hostile, x) * sp.diff(Q_hostile, y)
            - sp.diff(P_hostile, y) * sp.diff(Q_hostile, x)
            - 1
        ),
        x,
        y,
    )
    gate(hostile_debt.coeff_monomial(y**2) == 2 * a * b, "hostile terminal debt changed")
    gate(
        Counter(hostile_debt.monoms())
        == Counter(((4, 0), (2, 1), (2, 0), (0, 2), (0, 1))),
        "hostile debt support changed",
    )
    gate(
        sp.expand(
            hostile_debt.as_expr().subs(
                {b: -a / 2, d: -3 * c, e: -9 * c**2 / (4 * a)}
            )
        )
        == -a**2 * y**2,
        "near-survivor coefficient cancellation changed",
    )

    # Positive support-drop control.
    alpha, delta = sp.symbols("alpha delta")
    P_control = x + alpha * y**2
    Q_control = y + delta * P_control**2
    J_control = sp.expand(
        sp.diff(P_control, x) * sp.diff(Q_control, y)
        - sp.diff(P_control, y) * sp.diff(Q_control, x)
    )
    gate(J_control == 1, "one-by-three tame control failed")

    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    gate(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "inactive Python assert found",
    )

    semantic_rows = [
        f"linear={sorted(linear_exits.items())}",
        f"residual={len(residual)}",
        f"saturation={saturation_attempts},{unit_ideals}",
        f"cache={len(groebner_cache)},{groebner_cache_hits}",
        f"profiles={sorted(nonunit_profiles.items(), key=str)}",
        f"dimensions={sorted(nonunit_dimensions.items())}",
        *(f"basis={count}:{shape}" for shape, count in sorted(nonunit_basis_shapes.items())),
        f"obstruction={obstruction}",
        f"bounded={bounded_monomials},{bounded_universe},{bounded_survivors}",
    ]
    semantic_digest = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("THM-3694 exact companion -- one-pure-x-Q normalized 2x3 closure")
    print(f"systems={system_count}")
    print("linear_exits=" + repr(sorted(linear_exits.items())))
    print(f"linear_exit_total={sum(linear_exits.values())}")
    print(f"residual={len(residual)}")
    print(f"saturation_attempts={saturation_attempts}")
    print(f"saturation_skips={saturation_skips}")
    print(f"unit_ideals={unit_ideals}")
    print(f"unique_groebner_systems={len(groebner_cache)}")
    print(f"groebner_cache_hits={groebner_cache_hits}")
    print(f"nonunit={saturation_attempts-unit_ideals}")
    print("nonunit_profiles=" + repr(sorted(nonunit_profiles.items(), key=str)))
    print("nonunit_dimensions=" + repr(sorted(nonunit_dimensions.items())))
    print(f"unique_nonunit_basis_shapes={len(nonunit_basis_shapes)}")
    print(f"forced_exponent_obstruction={obstruction}")
    print("forced_exponent_roots=1/2,3/2")
    print(
        f"bounded_degree6=monomials:{bounded_monomials},universe:{bounded_universe},survivors:{bounded_survivors}"
    )
    print("controls=near_survivor_terminal:PASS,one_by_three_tame:PASS")
    print(f"semantic_sha256={semantic_digest}")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
