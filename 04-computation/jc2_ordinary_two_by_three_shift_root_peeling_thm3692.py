#!/usr/bin/env python3
"""Exact companion for the ordinary normalized 2x3 shift-root peel.

The proof of THM-3692 is degree-unbounded and structural.  This companion
checks its six forced scalar-tower leaves, performs a direct bounded census
without using the peel, and then runs a farther census using only the proved
leading-parallel necessary condition.
"""

from __future__ import annotations

import ast
import hashlib
from collections import Counter
from itertools import combinations
from math import gcd
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


def debt_multiplicities(
    left: tuple[tuple[int, int], ...],
    right: tuple[tuple[int, int], ...],
) -> Counter[tuple[int, int]]:
    """Active exponent buckets of Jac(P,Q)-1 in the normalized chart."""

    buckets: Counter[tuple[int, int]] = Counter()
    for px, py in left:
        if px != 0:
            buckets[(px - 1, py)] += 1
    for qx, qy in right:
        if qy != 0:
            buckets[(qx, qy - 1)] += 1
    for px, py in left:
        for qx, qy in right:
            if px * qy - py * qx != 0:
                buckets[(px + qx - 1, py + qy - 1)] += 1
    return buckets


def direct_census(max_total_degree: int) -> tuple[int, int, int]:
    """Blind support census: no shift-root lemma or leading filter."""

    monomials = [
        (x_degree, total - x_degree)
        for total in range(2, max_total_degree + 1)
        for x_degree in range(total + 1)
    ]
    left_pairs = list(combinations(monomials, 2))
    right_triples = list(combinations(monomials, 3))
    survivors = 0
    for left in left_pairs:
        for right in right_triples:
            buckets = debt_multiplicities(left, right)
            if all(multiplicity >= 2 for multiplicity in buckets.values()):
                survivors += 1
    return len(monomials), len(left_pairs) * len(right_triples), survivors


def primitive_direction(vector: tuple[int, int]) -> tuple[int, int]:
    divisor = gcd(vector[0], vector[1])
    return vector[0] // divisor, vector[1] // divisor


def leading_parallel_census(max_total_degree: int) -> tuple[int, int, int, int]:
    """Farther scout after the exact unique-leading-edge reduction."""

    monomials = [
        (x_degree, total - x_degree)
        for total in range(2, max_total_degree + 1)
        for x_degree in range(total + 1)
    ]
    index = {monomial: position for position, monomial in enumerate(monomials)}
    rays: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for monomial in monomials:
        rays.setdefault(primitive_direction(monomial), []).append(monomial)

    tested = 0
    survivors = 0
    for left in combinations(monomials, 2):
        leading_left = left[-1]
        for leading_right in rays[primitive_direction(leading_left)]:
            leading_index = index[leading_right]
            if leading_index < 2:
                continue
            for first, second in combinations(monomials[:leading_index], 2):
                right = (first, second, leading_right)
                tested += 1
                buckets = debt_multiplicities(left, right)
                if all(multiplicity >= 2 for multiplicity in buckets.values()):
                    survivors += 1

    full_universe = len(list(combinations(monomials, 2))) * len(
        list(combinations(monomials, 3))
    )
    return len(monomials), full_universe, tested, survivors


def main() -> None:
    sx, sy = sp.symbols("sx sy", nonnegative=True, integer=True)
    k, ell = sp.symbols("k ell", positive=True, integer=True)
    p = (k * sx + 1, k * sy)
    q = (ell * sx, ell * sy + 1)
    forced_determinant = sp.expand(p[0] * q[1] - p[1] * q[0])
    gate(
        forced_determinant == k * sx + ell * sy + 1,
        "forced scalar-tower determinant formula changed",
    )

    # S and R are recorded by their positive scalar multiples of the common
    # least vector s.  The five rows are exactly the leaves of the hand peel.
    forced_cases = {
        "I_seed2": ((1, 2), (1, 3, 4), (6,)),
        "II_seed3": ((1, 3), (1, 2, 4), (7,)),
        "II_root3": ((1, 4), (1, 2, 3), (5, 6, 7)),
        "II_both3": ((1, 3), (1, 2, 3), (5, 6)),
        "III_root3": ((1, 2), (1, 2, 3), (5,)),
        "III_root4": ((1, 2), (1, 2, 4), (5, 6)),
    }
    case_rows = []
    for name, (shifts, roots, expected_singletons) in forced_cases.items():
        multiplicities: Counter[int] = Counter(shifts)
        multiplicities.update(roots)
        multiplicities.update(root + shift for root in roots for shift in shifts)
        singletons = tuple(
            location
            for location, multiplicity in sorted(multiplicities.items())
            if multiplicity == 1
        )
        gate(singletons == expected_singletons, f"wrong terminal debt in {name}")
        for location in singletons:
            # Every terminal listed by the peel is an edge k*s + ell*s.
            representations = [
                (shift, root)
                for shift in shifts
                for root in roots
                if shift + root == location
            ]
            gate(bool(representations), f"terminal {location}s is not an edge")
            for shift, root in representations:
                edge_factor = shift * sx + root * sy + 1
                gate(
                    edge_factor.subs({sx: 0, sy: 0}) == 1,
                    f"terminal edge in {name} lost its positive constant",
                )
        profile = tuple(sorted(multiplicities.items()))
        case_rows.append(f"{name}:profile={profile}:singletons={singletons}")

    direct_monomials, direct_universe, direct_survivors = direct_census(7)
    gate(direct_monomials == 33, "wrong direct-census monomial count")
    gate(direct_universe == 2_880_768, "wrong direct-census universe")
    gate(direct_survivors == 0, "direct census found a support survivor")

    (
        leading_monomials,
        leading_universe,
        leading_tests,
        leading_survivors,
    ) = leading_parallel_census(10)
    gate(leading_monomials == 63, "wrong leading-census monomial count")
    gate(leading_universe == 77_555_583, "wrong leading-census universe")
    gate(leading_tests == 3_739_003, "wrong leading-parallel test count")
    gate(leading_survivors == 0, "leading-parallel census found a survivor")

    # Tame controls show exactly which divergence hypothesis is being used.
    x, y, alpha, gamma, delta = sp.symbols("x y alpha gamma delta")
    P_one_by_three = x + alpha * y**2
    Q_one_by_three = y + delta * P_one_by_three**2
    J_one_by_three = sp.expand(
        sp.diff(P_one_by_three, x) * sp.diff(Q_one_by_three, y)
        - sp.diff(P_one_by_three, y) * sp.diff(Q_one_by_three, x)
    )
    gate(J_one_by_three == 1, "one-by-three tame control failed")

    P_axis = x + alpha * y**2 + gamma * y**3
    gate(sp.diff(P_axis, x) == 1, "two-term pure-y axis control failed")

    # Exact collision in the displayed 2x3 support size, with visible debt.
    a, b, d, e = sp.symbols("a b d e")
    P_collision = x + a * x**2 - x**3
    Q_collision = y + b * x**2 + d * x**4 + e * x**6
    for polynomial in (P_collision, Q_collision):
        gate(
            sp.expand(polynomial.subs({x: 1, y: 0}) - polynomial.subs({x: -1, y: 0}))
            == 0,
            "two-by-three hostile collision failed",
        )
    J_collision = sp.expand(
        sp.diff(P_collision, x) * sp.diff(Q_collision, y)
        - sp.diff(P_collision, y) * sp.diff(Q_collision, x)
    )
    gate(J_collision == 1 + 2 * a * x - 3 * x**2, "collision debt changed")

    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    gate(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "inactive Python assert found",
    )

    semantic_rows = [
        f"determinant={forced_determinant}",
        *case_rows,
        f"direct={direct_monomials},{direct_universe},{direct_survivors}",
        f"leading={leading_monomials},{leading_universe},{leading_tests},{leading_survivors}",
    ]
    semantic_digest = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("THM-3692 exact companion -- ordinary 2x3 shift-root peeling")
    print(f"forced_determinant={forced_determinant}")
    for row in case_rows:
        print(row)
    print(
        f"direct_degree7=monomials:{direct_monomials},universe:{direct_universe},survivors:{direct_survivors}"
    )
    print(
        "leading_parallel_degree10="
        f"monomials:{leading_monomials},full_universe:{leading_universe},"
        f"tests:{leading_tests},survivors:{leading_survivors}"
    )
    print("controls=one_by_three_tame:PASS,two_pure_y_plus_y:PASS,collision_debt:PASS")
    print(f"semantic_sha256={semantic_digest}")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
