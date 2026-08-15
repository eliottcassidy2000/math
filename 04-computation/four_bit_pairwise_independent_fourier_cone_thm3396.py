#!/usr/bin/env python3
"""Exact referee for THM-3396.

The analytic theorem classifies laws on four Rademacher bits whose Walsh
moments of degrees one and two vanish.  This companion checks the Fourier
inversion and positivity facets on an exact rational grid, exhausts every
0/1-supported law on the sixteen atoms, checks the coordinatewise-convolution
law, reconstructs the complete face lattice of the polar odd 5-demicube, and
freezes the three higher-chaos sidecars extracted from the Hadamard certificate
bank.  It uses no floats, external packages, or ``assert`` gates.
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
Sign5 = tuple[int, int, int, int, int]
Vector5 = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]

ATOMS: tuple[Atom, ...] = tuple(product((-1, 1), repeat=4))
ATOM_INDEX = {atom: index for index, atom in enumerate(ATOMS)}
LOW_SUBSETS = tuple((i,) for i in range(4)) + tuple(combinations(range(4), 2))
ODD_NORMALS: tuple[Sign5, ...] = tuple(
    sign for sign in product((-1, 1), repeat=5) if math.prod(sign) == -1
)


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


def dot5(left: tuple[int | Fraction, ...], right: Vector5) -> Fraction:
    return sum((Fraction(a) * b for a, b in zip(left, right)), Fraction(0))


def solve_square(
    matrix: tuple[Sign5, ...], right_hand_side: tuple[int, ...]
) -> Vector5 | None:
    """Solve one 5 by 5 system over Q, returning None when singular."""
    augmented = [
        [Fraction(entry) for entry in row] + [Fraction(value)]
        for row, value in zip(matrix, right_hand_side)
    ]
    for column in range(5):
        pivot = next(
            (row for row in range(column, 5) if augmented[row][column]), None
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [entry / scale for entry in augmented[column]]
        for row in range(5):
            if row == column or not augmented[row][column]:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                entry - scale * pivot_entry
                for entry, pivot_entry in zip(augmented[row], augmented[column])
            ]
    return tuple(row[-1] for row in augmented)  # type: ignore[return-value]


def matrix_rank(rows: tuple[tuple[Fraction, ...], ...], width: int = 5) -> int:
    matrix = [list(row) for row in rows]
    rank = 0
    for column in range(width):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]), None
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = matrix[rank][column]
        matrix[rank] = [entry / scale for entry in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                entry - scale * pivot_entry
                for entry, pivot_entry in zip(matrix[row], matrix[rank])
            ]
        rank += 1
    return rank


def affine_dimension(points: tuple[Vector5, ...]) -> int:
    require(points, "empty point set has no affine dimension here")
    base = points[0]
    differences = tuple(
        tuple(entry - base_entry for entry, base_entry in zip(point, base))
        for point in points[1:]
    )
    return matrix_rank(differences)


def points_from_mask(vertices: tuple[Vector5, ...], mask: int) -> tuple[Vector5, ...]:
    return tuple(vertex for index, vertex in enumerate(vertices) if mask & (1 << index))


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
    feasible_grid_packets: list[tuple[Fraction, ...]] = []
    for packet in product(grid, repeat=5):
        cubic = packet[:4]
        quartic = packet[4]
        reconstructed = reconstruct_counts(Fraction(1), cubic, quartic)
        pointwise = all(value >= 0 for value in reconstructed)
        facets = cone_inequalities(cubic, quartic)
        require(pointwise == facets, ("facet iff failed", packet))
        grid_cells += 1
        feasible_cells += int(pointwise)
        if pointwise:
            feasible_grid_packets.append(packet)

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

    # The moment polytope is the polar of the odd 5-demicube.  Enumerate all
    # rank-five intersections of its sixteen atom facets over Q and compare
    # them with the analytic list: ten coordinate vertices and sixteen
    # one-third even-sign vertices.
    enumerated_vertices: set[Vector5] = set()
    for facet_indices in combinations(range(len(ODD_NORMALS)), 5):
        matrix = tuple(ODD_NORMALS[index] for index in facet_indices)
        candidate = solve_square(matrix, (1, 1, 1, 1, 1))
        if candidate is None:
            continue
        if all(dot5(normal, candidate) <= 1 for normal in ODD_NORMALS):
            enumerated_vertices.add(candidate)

    coordinate_vertices: set[Vector5] = set()
    for index in range(5):
        for sign in (-1, 1):
            coordinate_vertices.add(tuple(
                Fraction(sign if coordinate == index else 0)
                for coordinate in range(5)
            ))  # type: ignore[arg-type]
    simplex_vertices: set[Vector5] = {
        tuple(Fraction(sign, 3) for sign in signs)  # type: ignore[misc]
        for signs in product((-1, 1), repeat=5)
        if math.prod(signs) == 1
    }
    expected_vertices = coordinate_vertices | simplex_vertices
    require(enumerated_vertices == expected_vertices,
            ("odd-demicube polar vertex mismatch", enumerated_vertices ^ expected_vertices))
    vertices = tuple(sorted(enumerated_vertices))

    active_facets = tuple(
        frozenset(
            index for index, normal in enumerate(ODD_NORMALS)
            if dot5(normal, vertex) == 1
        )
        for vertex in vertices
    )
    for vertex, active in zip(vertices, active_facets):
        active_rows = tuple(
            tuple(Fraction(entry) for entry in ODD_NORMALS[index])
            for index in active
        )
        require(matrix_rank(active_rows) == 5, ("nonvertex candidate", vertex, active))

    # Intersections of subsets of the sixteen facets recover the entire face
    # lattice.  Store faces by their vertex-incidence bitmask, so duplicates
    # from different defining subsets disappear exactly.
    all_vertices_mask = (1 << len(vertices)) - 1
    facet_vertex_masks = tuple(
        sum(
            1 << vertex_index
            for vertex_index, vertex in enumerate(vertices)
            if dot5(normal, vertex) == 1
        )
        for normal in ODD_NORMALS
    )
    face_masks = {all_vertices_mask}
    for facet_subset in range(1, 1 << len(ODD_NORMALS)):
        common = all_vertices_mask
        remaining = facet_subset
        while remaining and common:
            low_bit = remaining & -remaining
            facet_index = low_bit.bit_length() - 1
            common &= facet_vertex_masks[facet_index]
            remaining ^= low_bit
        if common:
            face_masks.add(common)
    face_dimensions = {
        mask: affine_dimension(points_from_mask(vertices, mask)) for mask in face_masks
    }
    face_histogram = Counter(face_dimensions.values())
    require(face_histogram == Counter({0: 26, 1: 120, 2: 160, 3: 80, 4: 16, 5: 1}),
            ("unexpected face lattice", face_histogram))

    # Absolute support on the moment polytope is the maximum of the coordinate
    # and one-third l1 norms.  The ten coordinate vertices realize the first
    # term; in odd dimension one of sign(c), -sign(c) is an even sign vector
    # and realizes the second term in absolute value.
    support_grid = tuple(product(range(-2, 3), repeat=5))
    for functional in support_grid:
        vertex_maximum = max(abs(dot5(functional, vertex)) for vertex in vertices)
        predicted = max(
            Fraction(max(abs(value) for value in functional)),
            Fraction(sum(abs(value) for value in functional), 3),
        )
        require(vertex_maximum == predicted,
                ("absolute support norm", functional, vertex_maximum, predicted))

    # The two vertex orbits have literal OA compilers.  Coordinate vertices
    # are uniform eight-run half-cubes.  One-third even-sign vertices are
    # twelve-run arrays with atom multiplicities 0^5,1^10,2^1.
    vertex_packet_histogram: Counter[tuple[int, tuple[tuple[int, int], ...]]] = Counter()
    for vertex in vertices:
        mass = 8 if vertex in coordinate_vertices else 12
        cubic = tuple(mass * value for value in vertex[:4])
        quartic = mass * vertex[4]
        counts = integer_tuple(
            reconstruct_counts(mass, cubic, quartic), "vertex OA compiler"
        )
        require(not any(low_moments(counts)), ("vertex OA low moments", vertex))
        vertex_packet_histogram[
            (mass, tuple(sorted(Counter(counts).items())))
        ] += 1
    require(
        vertex_packet_histogram
        == Counter({
            (8, ((0, 8), (1, 8))): 10,
            (12, ((0, 5), (1, 10), (2, 1))): 16,
        }),
        ("vertex OA orbit census", vertex_packet_histogram),
    )

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

    # Repeated coordinatewise multiplication is group convolution.  Its exact
    # chi-square distance from uniform is the Parseval sum of the five
    # surviving Fourier coefficients.  Check six powers on every feasible
    # half-grid packet, together with the sharp total-variation bound.
    mixing_checks = 0
    sharp_tv_checks = 0
    for packet in feasible_grid_packets:
        for power_index in range(1, 7):
            powered = tuple(value**power_index for value in packet)
            probabilities = reconstruct_counts(
                Fraction(1), powered[:4], powered[4]
            )
            chi_square = sum(
                (16 * probability - 1) ** 2 for probability in probabilities
            ) / 16
            expected_chi_square = sum(
                value ** (2 * power_index) for value in packet
            )
            require(chi_square == expected_chi_square,
                    ("convolution Parseval", packet, power_index))
            total_variation = sum(
                abs(probability - Fraction(1, 16))
                for probability in probabilities
            ) / 2
            require(4 * total_variation**2 <= expected_chi_square,
                    ("convolution TV bound", packet, power_index))
            if sum(value != 0 for value in packet) == 1:
                require(4 * total_variation**2 == expected_chi_square,
                        ("one-coordinate TV sharpness", packet, power_index))
                sharp_tv_checks += 1
            mixing_checks += 1

    nonmixing_vertices = {
        vertex for vertex in vertices if max(abs(value) for value in vertex) == 1
    }
    require(nonmixing_vertices == coordinate_vertices,
            ("nonmixing vertex classification", nonmixing_vertices))

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

    # The 896 packet has exactly two zero atoms.  Under the atom-to-facet map
    # x -> -q(x)(x,1), those are adjacent odd-demicube vertices.  Their common
    # face in the polar has five vertices and f-vector (5,9,6): a triangular
    # bipyramid.  Exactly two active facets place the packet in its relative
    # interior, rather than on one of the six triangular boundary faces.
    packet_896 = (
        Fraction(-1, 4), Fraction(0), Fraction(0), Fraction(-1, 4), Fraction(1, 2)
    )
    packet_896_active = tuple(
        index for index, normal in enumerate(ODD_NORMALS)
        if dot5(normal, packet_896) == 1
    )
    require(len(packet_896_active) == 2,
            ("896 active atom facets", packet_896_active))
    active_normals = {ODD_NORMALS[index] for index in packet_896_active}
    expected_active_normals = {
        (-1, -1, 1, -1, 1),
        (-1, 1, -1, -1, 1),
    }
    require(active_normals == expected_active_normals,
            ("896 active normals", active_normals))
    ridge_mask = all_vertices_mask
    for index in packet_896_active:
        ridge_mask &= facet_vertex_masks[index]
    ridge_points = points_from_mask(vertices, ridge_mask)
    expected_ridge_points: set[Vector5] = {
        (Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1)),
        (Fraction(-1), Fraction(0), Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0), Fraction(-1), Fraction(0)),
        tuple(Fraction(value, 3) for value in (-1, -1, -1, -1, 1)),
        tuple(Fraction(value, 3) for value in (-1, 1, 1, -1, 1)),
    }  # type: ignore[assignment]
    require(set(ridge_points) == expected_ridge_points,
            ("896 ridge vertices", ridge_points))
    require(affine_dimension(ridge_points) == 3, "896 face is not a ridge")
    ridge_face_histogram = Counter(
        dimension for mask, dimension in face_dimensions.items()
        if not (mask & ~ridge_mask)
    )
    require(ridge_face_histogram == Counter({0: 5, 1: 9, 2: 6, 3: 1}),
            ("896 ridge face lattice", ridge_face_histogram))
    require(
        all(
            len(points_from_mask(vertices, mask)) == 3
            for mask, dimension in face_dimensions.items()
            if dimension == 2 and not (mask & ~ridge_mask)
        ),
        "896 ridge has a nontriangular facet",
    )

    # The same 896-row law has two exact vertex decompositions, one entirely
    # in the H8 orbit and one mixing H8 with the two H12 apices.  This freezes
    # the nonuniqueness that a Gram or moment packet cannot remember.
    e5 = (Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(1))
    minus_e1 = (Fraction(-1), Fraction(0), Fraction(0), Fraction(0), Fraction(0))
    minus_e4 = (Fraction(0), Fraction(0), Fraction(0), Fraction(-1), Fraction(0))
    simplex_minus = tuple(Fraction(value, 3) for value in (-1, -1, -1, -1, 1))
    simplex_plus = tuple(Fraction(value, 3) for value in (-1, 1, 1, -1, 1))

    def weighted_sum(
        terms: tuple[tuple[Fraction, Vector5], ...]
    ) -> Vector5:
        return tuple(
            sum((weight * vertex[index] for weight, vertex in terms), Fraction(0))
            for index in range(5)
        )  # type: ignore[return-value]

    coordinate_decomposition = (
        (Fraction(1, 2), e5),
        (Fraction(1, 4), minus_e1),
        (Fraction(1, 4), minus_e4),
    )
    mixed_decomposition = (
        (Fraction(1, 4), e5),
        (Fraction(3, 8), simplex_minus),
        (Fraction(3, 8), simplex_plus),
    )
    require(weighted_sum(coordinate_decomposition) == packet_896,
            "896 coordinate-vertex decomposition")
    require(weighted_sum(mixed_decomposition) == packet_896,
            "896 mixed-vertex decomposition")

    packet_896_counts = integer_tuple(
        reconstruct_counts(896, (-224, 0, 0, -224), 448), "896 packet"
    )
    for label, blocks in (
        ("coordinate", ((448, e5), (224, minus_e1), (224, minus_e4))),
        ("mixed", ((224, e5), (336, simplex_minus), (336, simplex_plus))),
    ):
        combined = [0] * len(ATOMS)
        for mass, vertex in blocks:
            block = integer_tuple(
                reconstruct_counts(
                    mass,
                    tuple(mass * value for value in vertex[:4]),
                    mass * vertex[4],
                ),
                f"896 {label} block",
            )
            combined = [left + right for left, right in zip(combined, block)]
        require(tuple(combined) == packet_896_counts,
                ("896 row-level decomposition", label))

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
        "polar_vertices": len(vertices),
        "polar_f_vector": [face_histogram[dimension] for dimension in range(5)],
        "polar_absolute_support": "max(linf,l1/3)",
        "polar_support_grid": len(support_grid),
        "vertex_OA_orbits": sorted(vertex_packet_histogram.items()),
        "convolution_pairs": convolution_pairs,
        "convolution_mixing_checks": mixing_checks,
        "convolution_TV_sharp_checks": sharp_tv_checks,
        "convolution_nonmixing_vertices": len(nonmixing_vertices),
        "puzzle_packets": packet_rows,
        "puzzle_896_active_normals": sorted(active_normals),
        "puzzle_896_face_f_vector": [
            ridge_face_histogram[dimension] for dimension in range(3)
        ],
        "puzzle_896_vertex_decompositions": ["H8:448,224,224", "H8/H12:224,336,336"],
        "one_sided_hostile": "plus_passes_odd_fails",
        "quadratic_hostile": "x0_equals_x1",
    }
    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("THM-3396 FOUR-BIT PAIRWISE-INDEPENDENT FOURIER-CONE AUDIT")
    print(f"rational_half_grid={grid_cells} feasible={feasible_cells} facet_iff=PASS")
    print("binary_atom_subsets=65535 pairwise_laws=11 support_histogram=8:10,16:1")
    print("odd_5_demicube_polar_vertices=26 f_vector=26,120,160,80,16")
    print("polar_absolute_support=max(linf,l1/3) grid=3125 PASS")
    print("polar_vertex_OA_orbits=H8:10,H12:16")
    print(f"coordinatewise_convolution_pairs={convolution_pairs} moment_product=PASS")
    print(
        f"convolution_mixing_checks={mixing_checks} "
        + f"tv_sharp_checks={sharp_tv_checks} nonmixing_signed_halfcubes=10"
    )
    for row in packet_rows:
        print(
            "puzzle_OA=" + str(row["mass"])
            + " cubic=" + ",".join(map(str, row["cubic"]))
            + " quartic=" + str(row["quartic"])
            + " even_margin=" + row["even_margin"]
            + " odd_margin=" + row["odd_margin"]
        )
    print("puzzle_896_face=triangular_bipyramid vertices=5 edges=9 facets=6 active_atoms=2")
    print("puzzle_896_vertex_decompositions=H8(448,224,224)=H8/H12(224,336,336)")
    print("one_parity_only_hostile=PASS")
    print("quadratic_hypothesis_hostile=PASS")
    print(f"semantic_sha256={digest}")


if __name__ == "__main__":
    main()
