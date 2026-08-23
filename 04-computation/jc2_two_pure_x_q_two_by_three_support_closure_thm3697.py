#!/usr/bin/env python3
"""Exact companion for the two-pure-x-Q normalized 2x3 closure.

Every equality partition of the active Jacobian labels is solved over Q.
The remaining determinant, nonlinearity, and distinctness conditions are
encoded by principal-open saturation.  All surviving affine systems have
empty saturated support strata over C.
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
    """Direct typed census: one y-active and two pure-x Q supports."""

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
        for active in active_right:
            for pure_pair in combinations(pure_x_right, 2):
                right = (active,) + pure_pair
                universe += 1
                multiplicities: Counter[tuple[int, int]] = Counter()
                for px, py in left:
                    multiplicities[(px - 1, py)] += 1
                multiplicities[(active[0], active[1] - 1)] += 1
                for p in left:
                    for q in right:
                        if determinant(p, q) != 0:
                            multiplicities[
                                (p[0] + q[0] - 1, p[1] + q[1] - 1)
                            ] += 1
                if all(value >= 2 for value in multiplicities.values()):
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

    # Up to relabelling Q supports, q2 and q3 are precisely the pure-x terms.
    active_divergence_labels = (0, 1, 2)
    linear_exits: Counter[str] = Counter()
    residual = []
    system_count = 0

    for cross_mask in range(64):
        active_labels = active_divergence_labels + tuple(
            5 + index for index in range(6) if cross_mask & (1 << index)
        )
        for partition in no_singleton_partitions(active_labels):
            system_count += 1
            equations = [q2y, q3y]
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
            if reason is None and any(
                value == 0 for value in (vectors[0][0], vectors[1][0], vectors[2][1])
            ):
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
            "inconsistent": 9_504,
            "equal_P": 1_233,
            "equal_Q": 294,
            "Q1_low": 34,
            "P1_low": 24,
            "P2_low": 24,
            "active_cross_zero": 8,
        }
    )
    gate(system_count == 11_155, "wrong activity/partition system count")
    gate(linear_exits == expected_linear_exits, "wrong affine exit ledger")
    gate(sum(linear_exits.values()) == 11_121, "wrong affine exit total")
    gate(len(residual) == 34, "wrong residual-system count")

    saturation_attempts = 0
    saturation_skips = 0
    unit_ideals = 0
    groebner_cache: dict[tuple[tuple[str, ...], tuple[str, ...], str], sp.GroebnerBasis] = {}
    groebner_cache_hits = 0
    nonunit_profiles: Counter[tuple[int, int, tuple[int, ...]]] = Counter()

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
        active_factors.extend((vectors[0][0], vectors[1][0], vectors[2][1]))
        for vector in vectors:
            total = sp.expand(sum(vector))
            active_factors.append(total * (total - 1))

        # B differs from the two pure-x terms because B_y is active.  The two
        # pure-x terms differ exactly when their x coordinates differ.
        active_factors.append(sp.expand(vectors[3][0] - vectors[4][0]))
        left_witnesses = (
            sp.expand(vectors[0][0] - vectors[1][0]),
            sp.expand(vectors[0][1] - vectors[1][1]),
        )

        base_product = sp.prod(active_factors)
        for left_witness in left_witnesses:
            saturation_factor = sp.factor(base_product * left_witness)
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
            nonunit_profiles[(len(parameters), len(basis.polys), degrees)] += 1

    gate(saturation_attempts == 68, "wrong saturation-attempt count")
    gate(saturation_skips == 0, "unexpected identically empty witness chart")
    gate(unit_ideals == 68, "a saturated support stratum survived")
    gate(not nonunit_profiles, "unexpected nonunit Groebner profile")

    bounded_monomials, bounded_universe, bounded_survivors = bounded_axis_census(6)
    gate(bounded_monomials == 25, "wrong bounded monomial count")
    gate(bounded_universe == 38_000, "wrong bounded typed universe")
    gate(bounded_survivors == 0, "bounded axis census found a survivor")

    x, y, a, b, c, d, e = sp.symbols("x y a b c d e")

    # A canonical near-survivor.  The two pure-x Q terms are invisible to the
    # bracket because P is pure x; cancelling the first two product layers
    # forces the nonzero terminal x^3 debt.
    P_hostile = x + a * x**2 + c * x**3
    Q_hostile = y + b * x * y + d * x**2 + e * x**3
    hostile_debt = sp.expand(
        sp.diff(P_hostile, x) * sp.diff(Q_hostile, y)
        - sp.diff(P_hostile, y) * sp.diff(Q_hostile, x)
        - 1
    )
    gate(
        hostile_debt
        == sp.expand(
            (2 * a + b) * x + (3 * c + 2 * a * b) * x**2 + 3 * b * c * x**3
        ),
        "hostile debt expansion changed",
    )
    gate(
        sp.expand(hostile_debt.subs({b: -2 * a, c: sp.Rational(4, 3) * a**2}))
        == -8 * a**3 * x**3,
        "near-survivor terminal debt changed",
    )

    # Positive controls on the support boundary.
    alpha, delta = sp.symbols("alpha delta")
    P_control = x + alpha * y**2
    Q_control = y + delta * P_control**2
    J_control = sp.expand(
        sp.diff(P_control, x) * sp.diff(Q_control, y)
        - sp.diff(P_control, y) * sp.diff(Q_control, x)
    )
    gate(J_control == 1, "one-by-three tame control failed")
    P_axis_control = x
    Q_axis_control = y + d * x**2 + e * x**3
    J_axis_control = sp.expand(
        sp.diff(P_axis_control, x) * sp.diff(Q_axis_control, y)
        - sp.diff(P_axis_control, y) * sp.diff(Q_axis_control, x)
    )
    gate(J_axis_control == 1, "two-pure-x support-drop control failed")

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
        f"bounded={bounded_monomials},{bounded_universe},{bounded_survivors}",
        f"hostile={sp.factor(hostile_debt)}",
    ]
    semantic_digest = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("THM-3697 exact companion -- two-pure-x-Q normalized 2x3 closure")
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
    print(
        f"bounded_degree6=monomials:{bounded_monomials},universe:{bounded_universe},survivors:{bounded_survivors}"
    )
    print("controls=near_survivor_terminal:PASS,one_by_three_tame:PASS,two_pure_x_drop:PASS")
    print(f"semantic_sha256={semantic_digest}")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
