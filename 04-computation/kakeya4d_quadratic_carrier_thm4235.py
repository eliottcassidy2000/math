#!/usr/bin/env python3
"""Exact F_61^4 quadratic Kakeya carrier and boundary-chart audit.

The THM-4035 direction torus is

    D* = {[1:u:v:w] : u,v,w in F_61^*}.

For k=0,1,2,3 we attach one line to each direction by using
``s^2+t*s`` in the first k spatial coordinates and ``t*s`` in the
remaining coordinates.  The fully quadratic carrier is audited twice:
by its one-coordinate fibre law and by direct enumeration of all incidences.

We then complete the projective boundary.  If the first nonzero homogeneous
coordinate of a direction is coordinate j, normalize it to one and use

    (0,...,0,t,s_(j+1)^2+t*s_(j+1),...,s_3^2+t*s_3).

The four disjoint charts contain every direction of P^3(F_61) exactly once.
Their union and multiplicity histogram are again checked by two independent
paths: a pointwise inverse-fibre formula and direct line enumeration.

This is finite-field geometry only.  It contains no Euclidean angles, tube
thickness, shadings, two-ends, or induction on scale.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import comb
import sys


Q = 61
PHI = 44
NONZERO = tuple(range(1, Q))
FIELD = tuple(range(Q))

sys.stdout.reconfigure(newline="\n")


def require(condition: bool, label) -> None:
    if not condition:
        raise RuntimeError(label)


def multiplicative_order(value: int, prime: int) -> int:
    current = 1
    for exponent in range(1, prime):
        current = current * value % prime
        if current == 1:
            return exponent
    raise RuntimeError("order search exhausted")


def det_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    work = [list(row) for row in matrix]
    determinant = 1
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column] % prime),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column] % prime
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        for row in range(column + 1, len(work)):
            multiplier = work[row][column] * inverse % prime
            for entry in range(column, len(work)):
                work[row][entry] = (
                    work[row][entry] - multiplier * work[column][entry]
                ) % prime
    return determinant % prime


def zero_sum_count(nonzero_terms: int, target_is_zero: bool) -> int:
    """Number of nonzero r-tuples with prescribed zero/nonzero sum."""
    if nonzero_terms == 0:
        return int(target_is_zero)
    if target_is_zero:
        return (
            (Q - 1) ** nonzero_terms
            + (Q - 1) * (-1) ** nonzero_terms
        ) // Q
    return ((Q - 1) ** nonzero_terms - (-1) ** nonzero_terms) // Q


def hyperplane_intersection(constant_nonzero: int, active_variables: int) -> int:
    if not constant_nonzero and active_variables == 0:
        raise RuntimeError("zero linear form is not a hyperplane")
    if active_variables == 0:
        return 0
    active = zero_sum_count(active_variables, not constant_nonzero)
    return active * (Q - 1) ** (3 - active_variables)


def torus_hyperplane_histogram() -> dict[int, int]:
    histogram: Counter[int] = Counter()
    for active in range(1, 4):
        number = comb(3, active) * (Q - 1) ** (active - 1)
        histogram[hyperplane_intersection(0, active)] += number
    for active in range(4):
        number = comb(3, active) * (Q - 1) ** active
        histogram[hyperplane_intersection(1, active)] += number
    return dict(sorted(histogram.items()))


def placement_histogram(quadratic_coordinates: int) -> dict[int, int]:
    """Fibre-product formula for the torus placement with k folds."""
    k = quadratic_coordinates
    require(0 <= k <= 3, (k, "quadratic coordinate count"))
    histogram: Counter[int] = Counter()
    if k == 0:
        histogram[(Q - 1) ** 3] += 1
    else:
        histogram[(2**k) * (Q - 1) ** (3 - k)] += ((Q - 1) // 2) ** k
    for doubled in range(k + 1):
        multiplicity = 2**doubled
        points_per_nonzero_slice = (
            comb(k, doubled)
            * ((Q - 3) // 2) ** doubled
            * 2 ** (k - doubled)
            * (Q - 1) ** (3 - k)
        )
        histogram[multiplicity] += (Q - 1) * points_per_nonzero_slice
    return dict(sorted(histogram.items()))


def direct_torus_quadratic() -> tuple[dict[int, int], str]:
    """Enumerate the 216000*61 quadratic torus line incidences."""
    counts = bytearray(Q**4)
    values = tuple(
        tuple((s * s + t * s) % Q for s in NONZERO)
        for t in FIELD
    )
    q2 = Q * Q
    q3 = q2 * Q
    for t in FIELD:
        image = values[t]
        t_offset = t * q3
        for x1 in image:
            x1_offset = t_offset + x1 * q2
            for x2 in image:
                x2_offset = x1_offset + x2 * Q
                for x3 in image:
                    counts[x2_offset + x3] += 1
    histogram = Counter(counts)
    positive = {
        multiplicity: number
        for multiplicity, number in sorted(histogram.items())
        if multiplicity
    }
    return positive, sha256(counts).hexdigest()


def cube_determinant_histogram() -> dict[int, int]:
    vertices = tuple((1,) + signs for signs in product((-1, 1), repeat=3))
    histogram: Counter[int] = Counter()
    for indices in combinations(range(8), 4):
        determinant = det_mod(tuple(vertices[index] for index in indices), Q)
        signed = determinant if determinant <= Q // 2 else determinant - Q
        histogram[signed] += 1
    return dict(sorted(histogram.items()))


def all_field_fibre_table() -> tuple[tuple[tuple[int, ...], ...], ...]:
    table = []
    for t in FIELD:
        counts = [0] * Q
        for s in FIELD:
            counts[(s * s + t * s) % Q] += 1
        table.append(tuple(counts))
    return tuple(table)


def full_chart_formula_histogram(
    fibres: tuple[tuple[tuple[int, ...], ...], ...]
) -> dict[int, int]:
    """Pointwise inverse-fibre formula, independent of line enumeration."""
    histogram: Counter[int] = Counter()
    for x0 in FIELD:
        f0 = fibres[x0]
        for x1 in FIELD:
            f1 = fibres[x1]
            for x2 in FIELD:
                f2 = fibres[x2]
                base = f0[x1] * f0[x2]
                for x3 in FIELD:
                    multiplicity = base * f0[x3]
                    if x0 == 0:
                        multiplicity += f1[x2] * f1[x3]
                        if x1 == 0:
                            multiplicity += f2[x3]
                            if x2 == 0:
                                multiplicity += 1
                    histogram[multiplicity] += 1
    return {
        multiplicity: number
        for multiplicity, number in sorted(histogram.items())
        if multiplicity
    }


def direct_full_chart_carrier() -> tuple[dict[int, int], str]:
    """Enumerate every line-point incidence in the four boundary charts."""
    counts = bytearray(Q**4)
    values = tuple(
        tuple((s * s + t * s) % Q for s in FIELD)
        for t in FIELD
    )
    q2 = Q * Q
    q3 = q2 * Q

    # [1:s1:s2:s3]
    for t in FIELD:
        image = values[t]
        t_offset = t * q3
        for x1 in image:
            x1_offset = t_offset + x1 * q2
            for x2 in image:
                x2_offset = x1_offset + x2 * Q
                for x3 in image:
                    counts[x2_offset + x3] += 1

    # [0:1:s2:s3]
    for t in FIELD:
        image = values[t]
        t_offset = t * q2
        for x2 in image:
            x2_offset = t_offset + x2 * Q
            for x3 in image:
                counts[x2_offset + x3] += 1

    # [0:0:1:s3]
    for t in FIELD:
        for x3 in values[t]:
            counts[t * Q + x3] += 1

    # [0:0:0:1]
    for t in FIELD:
        counts[t] += 1

    histogram = Counter(counts)
    positive = {
        multiplicity: number
        for multiplicity, number in sorted(histogram.items())
        if multiplicity
    }
    return positive, sha256(counts).hexdigest()


def main() -> None:
    require(PHI * PHI % Q == (PHI + 1) % Q, "golden root")
    require(multiplicative_order(PHI, Q) == 60, "clock order")
    require({pow(PHI, phase, Q) for phase in range(60)} == set(NONZERO), "clock chart")

    torus_directions = (Q - 1) ** 3
    all_directions = Q**3 + Q**2 + Q + 1
    boundary_debt = all_directions - torus_directions
    require((torus_directions, all_directions, boundary_debt) == (216_000, 230_764, 14_764), "direction ledger")

    plane_histogram = torus_hyperplane_histogram()
    require(
        plane_histogram == {0: 4, 3540: 14_400, 3541: 216_000, 3600: 360},
        "torus plane histogram",
    )
    require(sum(plane_histogram.values()) == all_directions, "dual projective count")

    incidence_mass = torus_directions * Q
    placements = {}
    for k in range(4):
        histogram = placement_histogram(k)
        require(
            sum(multiplicity * number for multiplicity, number in histogram.items())
            == incidence_mass,
            (k, "torus incidence mass"),
        )
        placements[k] = (sum(histogram.values()), histogram)
    require(
        tuple(placements[k][0] for k in range(4))
        == (12_960_001, 6_696_030, 3_460_500, 1_814_460),
        "placement ladder",
    )

    torus_formula = placements[3][1]
    require(
        torus_formula == {1: 480, 2: 20_880, 4: 302_760, 8: 1_490_340},
        "torus multiplicity formula",
    )
    torus_direct, torus_digest = direct_torus_quadratic()
    require(torus_direct == torus_formula, "torus formula/direct mismatch")

    determinant_histogram = cube_determinant_histogram()
    require(
        determinant_histogram == {-16: 1, -8: 28, 0: 12, 8: 28, 16: 1},
        "cube determinant histogram",
    )
    transverse_per_cube = sum(
        number for determinant, number in determinant_histogram.items() if determinant
    )
    transverse_incidences = transverse_per_cube * torus_formula[8]
    require((transverse_per_cube, transverse_incidences) == (58, 86_439_720), "transverse quartets")

    fibres = all_field_fibre_table()
    for t in FIELD:
        require(Counter(fibres[t]) == Counter({0: 30, 1: 1, 2: 30}), (t, "full-field fibre law"))
    full_formula = full_chart_formula_histogram(fibres)
    expected_full = {
        1: 150,
        2: 9_450,
        4: 210_151,
        6: 120,
        7: 30,
        8: 1_641_540,
        9: 60,
        10: 1_260,
        12: 5_880,
    }
    require(full_formula == expected_full, "full-chart formula histogram")
    full_direct, full_digest = direct_full_chart_carrier()
    require(full_direct == full_formula, "full-chart formula/direct mismatch")
    full_union = sum(full_formula.values())
    full_incidence_mass = sum(
        multiplicity * number for multiplicity, number in full_formula.items()
    )
    require((full_union, full_incidence_mass) == (1_868_641, all_directions * Q), "full carrier ledger")

    ambient = Q**4
    asymptotic_chart_upper = Q * sum(((Q + 1) // 2) ** power for power in range(4))
    print("THM-4235 FINITE FOUR-DIMENSIONAL QUADRATIC KAKEYA CARRIER")
    print(
        "universe=(field:61,ambient_points:13845841,projective_directions:230764,"
        "torus_directions:216000,boundary_directions:14764)"
    )
    print(
        f"torus_plane_concentration=(histogram:{plane_histogram},maximum:3600,"
        f"maximum_fraction:{Fraction(3600, torus_directions)})"
    )
    print(
        "torus_placement_ladder=(concurrent:12960001,one_fold:6696030,"
        "two_folds:3460500,three_folds:1814460)"
    )
    print(
        f"torus_quadratic=(multiplicity_histogram:{torus_formula},"
        f"direct_sha256:{torus_digest},cube_determinants:{determinant_histogram},"
        f"transverse_quartets:{transverse_incidences})"
    )
    print(
        f"full_projective_carrier=(union:{full_union},ambient_fraction:"
        f"{Fraction(full_union, ambient)},multiplicity_histogram:{full_formula},"
        f"maximum_multiplicity:{max(full_formula)},incidences:{full_incidence_mass},"
        f"direct_sha256:{full_digest},chart_sum_upper:{asymptotic_chart_upper})"
    )
    print(
        "euclidean_polynomial_shadow=(map:(t,s1^2+t*s1,s2^2+t*s2,s3^2+t*s3),"
        "focal_determinant:(2*s1+t)*(2*s2+t)*(2*s3+t),"
        "det_DF_sign:-1,focal_degree:3)"
    )
    print(
        "scope_boundary=(finite_field_kakeya_set:True,euclidean_transfer:False,"
        "angles:False,tube_thickness:False,shadings:False,two_ends:False,scales:False)"
    )
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
