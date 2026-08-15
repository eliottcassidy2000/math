#!/usr/bin/env python3
"""Exact referee for THM-3396.

The analytic theorem classifies laws on four Rademacher bits whose Walsh
moments of degrees one and two vanish.  This companion checks the Fourier
inversion and positivity facets on an exact rational grid, exhausts every
0/1-supported law on the sixteen atoms, checks the coordinatewise-convolution
law, and freezes the three higher-chaos sidecars extracted from the Hadamard
certificate bank.  It uses no floats, external packages, or ``assert`` gates.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
import math


Atom = tuple[int, int, int, int]
Counts = tuple[int, ...]

ATOMS: tuple[Atom, ...] = tuple(product((-1, 1), repeat=4))
ATOM_INDEX = {atom: index for index, atom in enumerate(ATOMS)}
LOW_SUBSETS = tuple((i,) for i in range(4)) + tuple(combinations(range(4), 2))


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def character(atom: Atom, subset: tuple[int, ...]) -> int:
    return math.prod(atom[index] for index in subset)


def parity(atom: Atom) -> int:
    return math.prod(atom)


def low_moments(counts: Counts) -> tuple[int, ...]:
    return tuple(
        sum(count * character(atom, subset) for count, atom in zip(counts, ATOMS))
        for subset in LOW_SUBSETS
    )


def high_moments(counts: Counts) -> tuple[tuple[int, ...], int]:
    cubic = tuple(
        sum(
            count * math.prod(atom[j] for j in range(4) if j != i)
            for count, atom in zip(counts, ATOMS)
        )
        for i in range(4)
    )
    quartic = sum(count * parity(atom) for count, atom in zip(counts, ATOMS))
    return cubic, quartic


def reconstruct_counts(
    mass: int | Fraction,
    cubic: tuple[int | Fraction, ...],
    quartic: int | Fraction,
) -> tuple[Fraction, ...]:
    require(len(cubic) == 4, "wrong cubic packet length")
    return tuple(
        Fraction(
            mass
            + parity(atom)
            * (quartic + sum(coefficient * bit for coefficient, bit in zip(cubic, atom)))
        , 16)
        for atom in ATOMS
    )


def cone_inequalities(
    cubic: tuple[Fraction, ...], quartic: Fraction
) -> bool:
    even_max = max(
        abs(sum(coefficient * bit for coefficient, bit in zip(cubic, atom)))
        for atom in ATOMS
        if parity(atom) == 1
    )
    odd_max = max(
        abs(sum(coefficient * bit for coefficient, bit in zip(cubic, atom)))
        for atom in ATOMS
        if parity(atom) == -1
    )
    return even_max <= 1 + quartic and odd_max <= 1 - quartic


def coordinatewise_convolution(left: Counts, right: Counts) -> Counts:
    result = [0] * len(ATOMS)
    for left_count, left_atom in zip(left, ATOMS):
        for right_count, right_atom in zip(right, ATOMS):
            target = tuple(a * b for a, b in zip(left_atom, right_atom))
            result[ATOM_INDEX[target]] += left_count * right_count
    return tuple(result)


def integer_tuple(values: tuple[Fraction, ...], label: str) -> Counts:
    require(all(value.denominator == 1 for value in values), (label, "nonintegral"))
    return tuple(value.numerator for value in values)


def main() -> None:
    # Exact character orthogonality on the four-cube.
    all_subsets = tuple(
        subset for size in range(5) for subset in combinations(range(4), size)
    )
    for left in all_subsets:
        for right in all_subsets:
            inner = sum(
                character(atom, left) * character(atom, right) for atom in ATOMS
            )
            require(inner == (16 if left == right else 0),
                    ("Walsh orthogonality", left, right, inner))

    # On a 1/2-grid, pointwise nonnegativity is exactly the two parity-facet
    # inequalities.  The analytic proof has no grid restriction.
    grid = tuple(Fraction(value, 2) for value in range(-2, 3))
    grid_cells = 0
    feasible_cells = 0
    for packet in product(grid, repeat=5):
        cubic = packet[:4]
        quartic = packet[4]
        reconstructed = reconstruct_counts(Fraction(1), cubic, quartic)
        pointwise = all(value >= 0 for value in reconstructed)
        facets = cone_inequalities(cubic, quartic)
        require(pointwise == facets, ("facet iff failed", packet))
        grid_cells += 1
        feasible_cells += int(pointwise)

    # Exhaust every subset of the sixteen atoms.  Exactly eleven nonempty
    # subsets give unbiased pairwise-independent laws: ten half-cubes and the
    # full cube.  Each must round-trip through only four cubics and one quartic.
    pairwise_laws: list[Counts] = []
    support_histogram: Counter[int] = Counter()
    for mask in range(1, 1 << len(ATOMS)):
        counts = tuple((mask >> index) & 1 for index in range(len(ATOMS)))
        if any(low_moments(counts)):
            continue
        mass = sum(counts)
        cubic, quartic = high_moments(counts)
        reconstructed = integer_tuple(
            reconstruct_counts(mass, cubic, quartic), "subset round trip"
        )
        require(reconstructed == counts, ("subset round trip changed", mask))
        normalized_cubic = tuple(Fraction(value, mass) for value in cubic)
        normalized_quartic = Fraction(quartic, mass)
        require(cone_inequalities(normalized_cubic, normalized_quartic),
                ("law outside cone", mask))
        pairwise_laws.append(counts)
        support_histogram[mass] += 1
    require(support_histogram == Counter({8: 10, 16: 1}),
            ("unexpected pairwise-law census", support_histogram))

    # Coordinatewise multiplication of independent rows multiplies every
    # Fourier coefficient.  Check all 11^2 exact subset-law pairs.
    convolution_pairs = 0
    for left in pairwise_laws:
        left_mass = sum(left)
        left_cubic, left_quartic = high_moments(left)
        for right in pairwise_laws:
            right_mass = sum(right)
            right_cubic, right_quartic = high_moments(right)
            convolved = coordinatewise_convolution(left, right)
            require(not any(low_moments(convolved)), "convolution lost pairwise independence")
            cubic, quartic = high_moments(convolved)
            require(sum(convolved) == left_mass * right_mass,
                    "convolution mass changed")
            require(cubic == tuple(a * b for a, b in zip(left_cubic, right_cubic)),
                    "cubic convolution law failed")
            require(quartic == left_quartic * right_quartic,
                    "quartic convolution law failed")
            convolution_pairs += 1

    # Three exact OA(N,4,2,2) sidecars occurring inside the Hadamard bank.
    # Coordinates are A_i=sum product_(j!=i) x_j and D=sum product_j x_j.
    puzzle_packets = (
        (48, (-16, -8, 0, 0), -8, Counter({1: 2, 2: 4, 3: 4, 4: 4, 5: 2})),
        (120, (0, 0, -8, 8), 24, Counter({5: 2, 6: 4, 7: 2, 8: 2, 9: 4, 10: 2})),
        (896, (-224, 0, 0, -224), 448,
         Counter({0: 2, 28: 4, 56: 4, 84: 4, 112: 2})),
    )
    packet_rows = []
    for mass, cubic, quartic, expected_histogram in puzzle_packets:
        counts = integer_tuple(
            reconstruct_counts(mass, cubic, quartic), f"packet {mass}"
        )
        require(all(value >= 0 for value in counts), ("negative packet", mass))
        require(sum(counts) == mass, ("packet mass", mass))
        require(not any(low_moments(counts)), ("packet low moment", mass))
        require(high_moments(counts) == (cubic, quartic), ("packet high moment", mass))
        require(Counter(counts) == expected_histogram, ("packet histogram", mass))
        normalized_cubic = tuple(Fraction(value, mass) for value in cubic)
        normalized_quartic = Fraction(quartic, mass)
        even_max = max(
            abs(sum(a * x for a, x in zip(normalized_cubic, atom)))
            for atom in ATOMS if parity(atom) == 1
        )
        odd_max = max(
            abs(sum(a * x for a, x in zip(normalized_cubic, atom)))
            for atom in ATOMS if parity(atom) == -1
        )
        require(cone_inequalities(normalized_cubic, normalized_quartic),
                ("packet facets", mass))
        packet_rows.append({
            "mass": mass,
            "cubic": list(cubic),
            "quartic": quartic,
            "count_histogram": sorted(expected_histogram.items()),
            "even_margin": str(1 + normalized_quartic - even_max),
            "odd_margin": str(1 - normalized_quartic - odd_max),
        })

    require(packet_rows[-1]["odd_margin"] == "0", "896 facet is no longer sharp")

    # Hostiles: one parity family of inequalities cannot replace both, and
    # vanishing degree-one moments alone cannot replace pairwise independence.
    one_sided_plus = (Fraction(1), Fraction(0), Fraction(0), Fraction(0))
    require(
        max(abs(sum(a * x for a, x in zip(one_sided_plus, atom)))
            for atom in ATOMS if parity(atom) == 1) <= Fraction(3, 2),
        "plus-side hostile lost its passing side",
    )
    require(
        any(value < 0 for value in reconstruct_counts(1, one_sided_plus, Fraction(1, 2))),
        "plus-side hostile no longer violates positivity",
    )
    correlated = tuple(int(atom[0] == atom[1]) for atom in ATOMS)
    require(all(low_moments(correlated)[i] == 0 for i in range(4)),
            "correlated hostile is not marginally unbiased")
    require(low_moments(correlated)[4] != 0,
            "correlated hostile lost its quadratic obstruction")
    correlated_cubic, correlated_quartic = high_moments(correlated)
    require(
        reconstruct_counts(sum(correlated), correlated_cubic, correlated_quartic)
        != tuple(Fraction(value) for value in correlated),
        "quadratic hostile was accidentally reconstructed",
    )

    semantic = {
        "theorem": "THM-3396",
        "grid_cells": grid_cells,
        "grid_feasible": feasible_cells,
        "subset_support_histogram": sorted(support_histogram.items()),
        "convolution_pairs": convolution_pairs,
        "puzzle_packets": packet_rows,
        "one_sided_hostile": "plus_passes_odd_fails",
        "quadratic_hostile": "x0_equals_x1",
    }
    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("THM-3396 FOUR-BIT PAIRWISE-INDEPENDENT FOURIER-CONE AUDIT")
    print(f"rational_half_grid={grid_cells} feasible={feasible_cells} facet_iff=PASS")
    print("binary_atom_subsets=65535 pairwise_laws=11 support_histogram=8:10,16:1")
    print(f"coordinatewise_convolution_pairs={convolution_pairs} moment_product=PASS")
    for row in packet_rows:
        print(
            "puzzle_OA=" + str(row["mass"])
            + " cubic=" + ",".join(map(str, row["cubic"]))
            + " quartic=" + str(row["quartic"])
            + " even_margin=" + row["even_margin"]
            + " odd_margin=" + row["odd_margin"]
        )
    print("puzzle_896_odd_parity_facet=SHARP two_atoms_absent")
    print("one_parity_only_hostile=PASS")
    print("quadratic_hypothesis_hostile=PASS")
    print(f"semantic_sha256={digest}")


if __name__ == "__main__":
    main()
