#!/usr/bin/env python3
"""Exact referee for THM-2457.

The script is deliberately dependency-free.  Every truth-bearing check uses
``require`` rather than ``assert``, so ``python`` and ``python -O`` execute the
same verification.
"""

from fractions import Fraction
from functools import reduce
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def valuation(n, prime):
    require(n != 0, "valuation requires a nonzero integer")
    n = abs(n)
    value = 0
    while n % prime == 0:
        n //= prime
        value += 1
    return value


def centered_norm(value):
    """Distance of a rational number to the nearest integer."""
    floor_value = value.numerator // value.denominator
    residue = value - floor_value
    return min(residue, 1 - residue)


def guard_safe(speed, point):
    return centered_norm(speed * point) > Fraction(1, 7)


def unit_safe(speed, point):
    return centered_norm(speed * point) > Fraction(1, 14)


def danger(speed, point):
    return centered_norm(speed * point) < Fraction(1, 14)


def cyclotomic_reduction_13(values, character):
    """Reduce the normalized C_13 transform modulo Phi_13.

    The return value is the coefficient vector in the basis
    1,zeta,...,zeta^11.
    """
    require(len(values) == 13, "C_13 vector must have thirteen entries")
    require(character % 13 != 0, "charged character must be nonzero")
    coefficients = [Fraction(0) for _ in range(13)]
    for shift, value in enumerate(values):
        exponent = (-character * shift) % 13
        coefficients[exponent] += value
    top = coefficients[12]
    return tuple((coefficients[index] - top) / 13 for index in range(12))


def graph_referee():
    atoms = 128
    directed_cells = atoms * atoms

    # Sharp realization: one rational base stratum for every ordered pair.
    matrix = [
        [Fraction(1, directed_cells) for _ in range(atoms)]
        for _ in range(atoms)
    ]
    total_mass = sum(
        (sum(row, Fraction(0)) for row in matrix),
        Fraction(0),
    )
    maximum_edge = max(max(row) for row in matrix)
    require(total_mass == 1, "sharp matrix must have initial mass one")
    require(
        maximum_edge == Fraction(1, directed_cells),
        "sharp directed-edge floor mismatch",
    )

    # A loop is exact-atom service.
    loop_matrix = [[Fraction(0) for _ in range(atoms)] for _ in range(atoms)]
    loop_matrix[5][5] = 1
    require(loop_matrix[5][5] > 0, "loop positive control failed")

    # A directed edge is two-atom service; an unrelated selected atom makes
    # the service union have at most three atoms.
    edge_matrix = [[Fraction(0) for _ in range(atoms)] for _ in range(atoms)]
    edge_matrix[5][9] = 1
    two_atom = {5, 9}
    two_atom_mass = sum(
        edge_matrix[left][right]
        for left in two_atom
        for right in two_atom
    )
    require(two_atom_mass == 1, "two-atom service identity failed")
    selected_atom = 17
    three_atom = {selected_atom, 5, 9}
    three_atom_mass = sum(
        edge_matrix[left][right]
        for left in three_atom
        for right in three_atom
    )
    require(len(three_atom) == 3, "three-atom cardinality control failed")
    require(three_atom_mass == 1, "three-atom service identity failed")
    selected_table = {selected_atom: Fraction(1, 3)}
    coarsened_table_mass = sum(
        selected_table.get(atom, Fraction(0)) for atom in three_atom
    )
    require(
        coarsened_table_mass == Fraction(1, 3),
        "positive selected table was not retained",
    )

    # Separate marginals are not simultaneous mass.
    marginal_a = Fraction(1, 2)
    marginal_f = Fraction(1, 2)
    simultaneous = Fraction(0)
    require(marginal_a > 0 and marginal_f > 0, "marginal hostile malformed")
    require(simultaneous == 0, "separate marginals spuriously coupled")

    # No drift-to-co-support floor without an explicit coupling.  The
    # two-state table (1,0) has constant-projection drift 1/4, while an
    # unrelated semantic co-support stratum can have arbitrary rational mass.
    drift_table = (Fraction(1), Fraction(0))
    drift_mean = sum(drift_table, Fraction(0)) / 2
    drift = sum(
        ((value - drift_mean) ** 2 for value in drift_table),
        Fraction(0),
    ) / 2
    epsilon_cosupport = Fraction(1, 1000)
    require(drift == Fraction(1, 4), "two-state drift control failed")
    require(
        epsilon_cosupport < drift,
        "no-floor hostile needs smaller co-support than drift",
    )

    return {
        "atoms": atoms,
        "directed_cells": directed_cells,
        "total_mass": total_mass,
        "maximum_edge": maximum_edge,
        "two_atom_mass": two_atom_mass,
        "three_atom_mass": three_atom_mass,
        "drift": drift,
        "epsilon_cosupport": epsilon_cosupport,
    }


def charged_energy_referee():
    # Sharp off-zero-uniform correlation with total mass one.
    correlation = [Fraction(0)] + [Fraction(1, 12)] * 12
    require(sum(correlation, Fraction(0)) == 1, "correlation mass mismatch")

    correlation_energy = (
        sum((value * value for value in correlation), Fraction(0)) / 13
        - Fraction(1, 169)
    )
    joint_energy = correlation_energy / 169
    maximum_joint_square = joint_energy / 12
    require(
        correlation_energy == Fraction(1, 2028),
        "charged correlation energy constant failed",
    )
    require(
        joint_energy == Fraction(1, 342732),
        "joint root energy constant failed",
    )
    require(
        maximum_joint_square == Fraction(1, 2028**2),
        "maximum joint coefficient constant failed",
    )

    expected_reduction = (Fraction(-1, 156),) + (Fraction(0),) * 11
    for character in range(1, 13):
        reduced = cyclotomic_reduction_13(correlation, character)
        require(
            reduced == expected_reduction,
            f"charged cyclotomic reduction failed at {character}",
        )
        joint_reduced = tuple(value / 13 for value in reduced)
        require(
            joint_reduced[0] == Fraction(-1, 2028)
            and all(value == 0 for value in joint_reduced[1:]),
            f"sharp joint coefficient failed at {character}",
        )

    global_root_denominator = 16384 * 2028
    global_energy_denominator = 342732 * 16384**2
    word_root_denominator = 65536 * 2028
    adaptive_denominators = tuple(
        2028 * count for count in (1, 16, 8, 4, 2, 1)
    )
    require(
        global_root_denominator == 33226752,
        "global maximum-root denominator failed",
    )
    require(
        global_energy_denominator == 92001420705792,
        "global root-energy denominator failed",
    )
    require(
        word_root_denominator == 132907008,
        "word maximum-root denominator failed",
    )
    require(
        adaptive_denominators == (2028, 32448, 16224, 8112, 4056, 2028),
        "adaptive root denominators failed",
    )

    return {
        "correlation_energy": correlation_energy,
        "joint_energy": joint_energy,
        "maximum_joint": Fraction(1, 2028),
        "global_root_denominator": global_root_denominator,
        "global_energy_denominator": global_energy_denominator,
        "word_root_denominator": word_root_denominator,
        "adaptive_denominators": adaptive_denominators,
    }


def hostile_margins(point, terminal_index):
    guard = 9
    units = (7, 8, 10, 11, 12)
    source = 13
    other = 338
    deepest = 1113879
    clock = 169
    terminal = clock * point

    margins = [centered_norm(guard * point) - Fraction(1, 7)]
    margins.extend(
        centered_norm(speed * point) - Fraction(1, 14) for speed in units
    )
    margins.extend(
        (
            Fraction(1, 14) - centered_norm(source * point),
            centered_norm(other * point) - Fraction(1, 14),
            centered_norm(deepest * point) - Fraction(1, 14),
            Fraction(1, 14)
            - centered_norm(deepest * point - Fraction(2, 13)),
        )
    )
    margins.append(centered_norm(guard * terminal) - Fraction(1, 7))
    margins.extend(
        centered_norm(speed * terminal) - Fraction(1, 14)
        for speed in units
    )
    margins.extend(
        (
            centered_norm(source * terminal) - Fraction(1, 14),
            Fraction(1, 14) - centered_norm(other * terminal),
        )
    )
    if terminal_index == 0:
        margins.append(
            Fraction(1, 14) - centered_norm(deepest * terminal)
        )
    else:
        margins.append(
            centered_norm(deepest * terminal) - Fraction(1, 14)
        )
    require(all(margin > 0 for margin in margins), "hostile hit a seam")
    return min(margins)


def semantic_word_hostile_referee():
    guard = 9
    units = (7, 8, 10, 11, 12)
    q_star = 7
    source = 13
    other = 338
    deepest = 1113879
    clock = 169
    speeds = (guard,) + units + (source, other, deepest)
    points = (
        Fraction(1041874, 14480427),
        Fraction(135443621, 1882455510),
    )

    require(deepest == 3 * 13**5, "deepest blocker factorization failed")
    require(clock == 13**2, "source clock failed")
    require(len(speeds) == 9, "speed count failed")
    require(len(set(speeds)) == 9, "speeds must be distinct")
    require(reduce(gcd, speeds) == 1, "speed row must be primitive")
    require(
        (valuation(source, 13), valuation(other, 13), valuation(deepest, 13))
        == (1, 2, 5),
        "strict blocker profile failed",
    )
    require(valuation(q_star, 7) == 1, "q_star septimal depth failed")
    require(valuation(deepest, 7) == 0, "deepest septimal depth failed")

    masks = []
    words = []
    minimum_margins = []
    deep_norms = []
    deep_probe_norms = []
    terminal_source_norms = []
    terminal_other_norms = []
    terminal_deep_norms = []

    for index, point in enumerate(points):
        require(0 < point < 1, "source point must lie on the base circle")
        require(guard_safe(guard, point), "source guard must be safe")
        require(
            all(unit_safe(speed, point) for speed in units),
            "source ordinary role must be safe",
        )
        require(danger(source, point), "source owner must be dangerous")
        require(not danger(other, point), "other blocker must be safe")
        require(not danger(deepest, point), "deepest blocker must be safe")
        require(unit_safe(q_star, point), "q_star graft must be safe")

        split_roles = (guard, 8, 10, 11, 12)
        split_bits = tuple(
            0
            if (
                guard_safe(speed, point)
                if speed == guard
                else unit_safe(speed, point)
            )
            else 1
            for speed in split_roles
        )
        blocker_bits = (
            1 if danger(source, point) else 0,
            1 if danger(other, point) else 0,
        )
        mask = split_bits + blocker_bits
        require(mask == (0, 0, 0, 0, 0, 1, 0), "local mask failed")
        masks.append(mask)

        deep_norm = centered_norm(deepest * point)
        deep_probe_norm = centered_norm(
            deepest * point - Fraction(2, 13)
        )
        require(deep_norm > Fraction(1, 14), "deep safe label failed")
        require(
            deep_probe_norm < Fraction(1, 14),
            "deep shifted probe failed",
        )
        deep_norms.append(deep_norm)
        deep_probe_norms.append(deep_probe_norm)

        terminal = clock * point
        require(guard_safe(guard, terminal), "terminal guard must be safe")
        require(
            all(unit_safe(speed, terminal) for speed in units),
            "terminal ordinary role must be safe",
        )
        require(not danger(source, terminal), "source must expire")
        require(danger(other, terminal), "other blocker must be dangerous")
        terminal_source_norms.append(centered_norm(source * terminal))
        terminal_other_norms.append(centered_norm(other * terminal))
        terminal_deep_norms.append(centered_norm(deepest * terminal))

        word = ["a"]
        if danger(deepest, terminal):
            word.append("c3")
        words.append(tuple(word))
        minimum_margins.append(hostile_margins(point, index))

    require(masks[0] == masks[1], "the hostile masks must coincide")
    require(
        deep_norms == [Fraction(2, 13), Fraction(261, 1690)],
        "deep safe norms failed",
    )
    require(
        deep_probe_norms == [Fraction(0), Fraction(1, 1690)],
        "deep probe norms failed",
    )
    require(
        terminal_source_norms
        == [Fraction(496, 6591), Fraction(64481, 856830)],
        "terminal source norms failed",
    )
    require(
        terminal_other_norms
        == [Fraction(22, 507), Fraction(1429, 32955)],
        "terminal other-blocker norms failed",
    )
    require(
        terminal_deep_norms == [Fraction(0), Fraction(1, 10)],
        "terminal deepest norms failed",
    )
    require(
        words == [("a", "c3"), ("a",)],
        "pure/fork semantic split failed",
    )
    require(
        minimum_margins
        == [Fraction(353, 92274), Fraction(11476, 2998905)],
        "exact hostile margins failed",
    )

    radii = tuple(
        margin / (2 * clock * deepest) for margin in minimum_margins
    )
    require(all(radius > 0 for radius in radii), "open hostile radius failed")

    return {
        "speeds": speeds,
        "profile": (1, 2, 5),
        "qstar_depth": valuation(q_star, 7),
        "deep_depth": valuation(deepest, 7),
        "mask": masks[0],
        "deep_norms": tuple(deep_norms),
        "deep_probe_norms": tuple(deep_probe_norms),
        "terminal_source_norms": tuple(terminal_source_norms),
        "terminal_other_norms": tuple(terminal_other_norms),
        "terminal_deep_norms": tuple(terminal_deep_norms),
        "words": tuple(words),
        "minimum_margins": tuple(minimum_margins),
        "radii": radii,
    }


def main():
    graph = graph_referee()
    charged = charged_energy_referee()
    hostile = semantic_word_hostile_referee()

    print("THM2457 COMPLETE-ATOM ROOT CO-SUPPORT REFEREE")
    print(f"atoms={graph['atoms']}")
    print(f"directed_matrix_cells={graph['directed_cells']}")
    print(f"sharp_matrix_total_mass={graph['total_mass']}")
    print(f"sharp_max_edge={graph['maximum_edge']}")
    print(f"two_atom_mass={graph['two_atom_mass']}")
    print(f"three_atom_mass={graph['three_atom_mass']}")
    print(
        "split_base_no_floor="
        f"drift:{graph['drift']},cosupport:{graph['epsilon_cosupport']}"
    )
    print(f"sharp_correlation_energy={charged['correlation_energy']}")
    print(f"sharp_joint_root_energy={charged['joint_energy']}")
    print(f"sharp_max_joint_coefficient={charged['maximum_joint']}")
    print(f"global_root_denominator={charged['global_root_denominator']}")
    print(
        "global_energy_denominator="
        f"{charged['global_energy_denominator']}"
    )
    print(f"word_root_denominator={charged['word_root_denominator']}")
    print(f"adaptive_denominators={charged['adaptive_denominators']}")
    print(f"hostile_speeds={hostile['speeds']}")
    print(f"hostile_profile={hostile['profile']}")
    print(
        "hostile_septimal_depths="
        f"(qstar:{hostile['qstar_depth']},c3:{hostile['deep_depth']})"
    )
    print(f"hostile_complete_mask={hostile['mask']}")
    print(f"hostile_deep_norms={hostile['deep_norms']}")
    print(f"hostile_deep_probe_norms={hostile['deep_probe_norms']}")
    print(
        "hostile_terminal_source_norms="
        f"{hostile['terminal_source_norms']}"
    )
    print(
        "hostile_terminal_other_norms="
        f"{hostile['terminal_other_norms']}"
    )
    print(f"hostile_terminal_deep_norms={hostile['terminal_deep_norms']}")
    print(f"hostile_words={hostile['words']}")
    print(f"hostile_minimum_margins={hostile['minimum_margins']}")
    print(f"hostile_open_radii={hostile['radii']}")
    print("fraction_only_truth_checks=PASS")
    print("ALL_CHECKS_PASSED")


if __name__ == "__main__":
    main()
