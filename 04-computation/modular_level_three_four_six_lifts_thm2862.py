#!/usr/bin/env python3
"""Exact modular level-three/four and inequivalent six-lift audit.

Dependency-free companion for THM-2862.  All truth-bearing checks use
explicit exceptions so that ``python -O`` performs the same verification.
"""

from collections import Counter
from itertools import combinations, permutations, product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def canonical_matrix(matrix, modulus):
    matrix = tuple(entry % modulus for entry in matrix)
    negative = tuple((-entry) % modulus for entry in matrix)
    return min(matrix, negative)


def determinant(matrix, modulus):
    return (
        matrix[0] * matrix[3] - matrix[1] * matrix[2]
    ) % modulus


def matrix_product(left, right, modulus):
    return canonical_matrix(
        (
            left[0] * right[0] + left[1] * right[2],
            left[0] * right[1] + left[1] * right[3],
            left[2] * right[0] + left[3] * right[2],
            left[2] * right[1] + left[3] * right[3],
        ),
        modulus,
    )


def matrix_inverse(matrix, modulus):
    return canonical_matrix(
        (matrix[3], -matrix[1], -matrix[2], matrix[0]),
        modulus,
    )


def sl_size(modulus):
    return sum(
        determinant(matrix, modulus) == 1 % modulus
        for matrix in product(range(modulus), repeat=4)
    )


def projective_sl(modulus):
    return tuple(
        sorted(
            {
                canonical_matrix(matrix, modulus)
                for matrix in product(range(modulus), repeat=4)
                if determinant(matrix, modulus) == 1 % modulus
            }
        )
    )


def generated_subgroup(generators, modulus):
    identity = canonical_matrix((1, 0, 0, 1), modulus)
    subgroup = {identity}
    pending = [identity]
    while pending:
        element = pending.pop()
        for generator in generators:
            candidate = matrix_product(element, generator, modulus)
            if candidate not in subgroup:
                subgroup.add(candidate)
                pending.append(candidate)
    return frozenset(subgroup)


def element_order(element, modulus):
    identity = canonical_matrix((1, 0, 0, 1), modulus)
    power = identity
    for exponent in range(1, 101):
        power = matrix_product(power, element, modulus)
        if power == identity:
            return exponent
    raise RuntimeError("element order exceeded exact cap")


def conjugate_element(conjugator, element, modulus):
    return matrix_product(
        matrix_product(conjugator, element, modulus),
        matrix_inverse(conjugator, modulus),
        modulus,
    )


def canonical_point(vector, modulus):
    units = tuple(
        scalar
        for scalar in range(modulus)
        if gcd(scalar, modulus) == 1
    )
    return min(
        tuple((scalar * coordinate) % modulus for coordinate in vector)
        for scalar in units
    )


def projective_line(modulus):
    return tuple(
        sorted(
            {
                canonical_point((x_coordinate, y_coordinate), modulus)
                for x_coordinate, y_coordinate
                in product(range(modulus), repeat=2)
                if gcd(gcd(x_coordinate, y_coordinate), modulus) == 1
            }
        )
    )


def matrix_on_point(matrix, point, modulus):
    return canonical_point(
        (
            (
                matrix[0] * point[0]
                + matrix[1] * point[1]
            ) % modulus,
            (
                matrix[2] * point[0]
                + matrix[3] * point[1]
            ) % modulus,
        ),
        modulus,
    )


def point_action(group, points, modulus):
    index = {point: position for position, point in enumerate(points)}
    return {
        element: tuple(
            index[matrix_on_point(element, point, modulus)]
            for point in points
        )
        for element in group
    }


def compose_permutations(left, right):
    return tuple(left[right[position]] for position in range(len(left)))


def cycle_type(permutation):
    visited = set()
    lengths = []
    for start in range(len(permutation)):
        if start in visited:
            continue
        current = start
        length = 0
        while current not in visited:
            visited.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def permutation_sign(permutation):
    return -1 if (
        sum(length - 1 for length in cycle_type(permutation)) % 2
    ) else 1


def action_census(action):
    return Counter(cycle_type(permutation) for permutation in action.values())


def check_action(group, action, modulus, label):
    identity = canonical_matrix((1, 0, 0, 1), modulus)
    identity_permutation = tuple(range(len(next(iter(action.values())))))
    require(
        action[identity] == identity_permutation,
        f"{label}: identity action failed",
    )
    for left in group:
        for right in group:
            require(
                action[matrix_product(left, right, modulus)]
                == compose_permutations(action[left], action[right]),
                f"{label}: action law failed",
            )


def action_orbit(action, point=0):
    return {permutation[point] for permutation in action.values()}


def all_two_generated_subgroups(group, modulus):
    return {
        generated_subgroup((left, right), modulus)
        for left in group
        for right in group
    }


def conjugate_subgroup(conjugator, subgroup, modulus):
    return frozenset(
        conjugate_element(conjugator, element, modulus)
        for element in subgroup
    )


def subgroup_conjugation_action(group, subgroups, modulus):
    index = {
        subgroup: position
        for position, subgroup in enumerate(subgroups)
    }
    return {
        element: tuple(
            index[conjugate_subgroup(element, subgroup, modulus)]
            for subgroup in subgroups
        )
        for element in group
    }


def element_conjugation_action(group, elements, modulus):
    index = {
        element: position
        for position, element in enumerate(elements)
    }
    return {
        conjugator: tuple(
            index[conjugate_element(conjugator, element, modulus)]
            for element in elements
        )
        for conjugator in group
    }


def subset_action(group, base_action, subset_size):
    base_size = len(next(iter(base_action.values())))
    subsets = tuple(combinations(range(base_size), subset_size))
    index = {
        subset: position
        for position, subset in enumerate(subsets)
    }
    action = {
        element: tuple(
            index[
                tuple(
                    sorted(
                        base_action[element][point]
                        for point in subset
                    )
                )
            ]
            for subset in subsets
        )
        for element in group
    }
    return subsets, action


def block_action(group, action, blocks):
    blocks = tuple(frozenset(block) for block in blocks)
    index = {
        block: position
        for position, block in enumerate(blocks)
    }
    induced = {}
    for element in group:
        permutation = action[element]
        induced[element] = tuple(
            index[
                frozenset(
                    permutation[position]
                    for position in block
                )
            ]
            for block in blocks
        )
    return induced


def equivariant_bijections(group, source_action, target_action):
    source_size = len(next(iter(source_action.values())))
    target_size = len(next(iter(target_action.values())))
    if source_size != target_size:
        return ()
    witnesses = []
    for candidate in permutations(range(target_size)):
        if all(
            candidate[source_action[element][position]]
            == target_action[element][candidate[position]]
            for element in group
            for position in range(source_size)
        ):
            witnesses.append(candidate)
    return tuple(witnesses)


def orbit_map_from_base(
    group,
    source_action,
    source_base,
    target_base,
    target_transport,
):
    """Build and verify g*x -> g*target_base is well defined."""
    source_size = len(next(iter(source_action.values())))
    images = [set() for _ in range(source_size)]
    for element in group:
        source_point = source_action[element][source_base]
        images[source_point].add(target_transport(element, target_base))
    require(
        all(len(image) == 1 for image in images),
        "orbit map is not well defined",
    )
    return tuple(next(iter(image)) for image in images)


def subgroup_element_orders(subgroup, modulus):
    return Counter(element_order(element, modulus) for element in subgroup)


def character_with_kernel(group, kernel):
    return {
        element: 1 if element in kernel else -1
        for element in group
    }


def commutator(left, right, modulus):
    return matrix_product(
        matrix_product(
            matrix_product(left, right, modulus),
            matrix_inverse(left, modulus),
            modulus,
        ),
        matrix_inverse(right, modulus),
        modulus,
    )


def check_character(group, character, modulus, label):
    for left in group:
        for right in group:
            require(
                character[matrix_product(left, right, modulus)]
                == character[left] * character[right],
                f"{label}: character law failed",
            )


def format_census(census):
    return ", ".join(
        f"{cycle}:{census[cycle]}"
        for cycle in sorted(census, reverse=True)
    )


def main():
    s_integral = (0, -1, 1, 0)
    c_integral = (0, -1, 1, 1)

    # Level three: the projective four-set is the four-C3-complement torsor.
    group_three = projective_sl(3)
    identity_three = canonical_matrix((1, 0, 0, 1), 3)
    s_three = canonical_matrix(s_integral, 3)
    c_three = canonical_matrix(c_integral, 3)
    t_three = matrix_product(
        matrix_inverse(s_three, 3),
        c_three,
        3,
    )
    require(
        (
            sl_size(3),
            len(group_three),
            len(generated_subgroup((s_three, c_three), 3)),
        ) == (24, 12, 12),
        "level-three group sizes failed",
    )
    require(
        (
            element_order(s_three, 3),
            element_order(c_three, 3),
            element_order(t_three, 3),
        ) == (2, 3, 3),
        "level-three triangle orders failed",
    )
    line_three = projective_line(3)
    projective_three = point_action(group_three, line_three, 3)
    check_action(
        group_three,
        projective_three,
        3,
        "level-three projective",
    )
    require(
        len(line_three) == 4
        and len(set(projective_three.values())) == 12
        and len(action_orbit(projective_three)) == 4,
        "level-three projective action failed",
    )
    expected_a4_census = Counter(
        {
            (3, 1): 8,
            (2, 2): 3,
            (1, 1, 1, 1): 1,
        }
    )
    require(
        action_census(projective_three) == expected_a4_census,
        "level-three A4 census failed",
    )
    klein_three = frozenset(
        element
        for element in group_three
        if element_order(element, 3) <= 2
    )
    require(
        len(klein_three) == 4
        and subgroup_element_orders(klein_three, 3)
        == Counter({2: 3, 1: 1}),
        "level-three Klein four failed",
    )
    subgroups_three = all_two_generated_subgroups(group_three, 3)
    complements_three = tuple(
        sorted(
            (
                subgroup
                for subgroup in subgroups_three
                if len(subgroup) == 3
                and subgroup.intersection(klein_three)
                == {identity_three}
            ),
            key=lambda subgroup: tuple(sorted(subgroup)),
        )
    )
    require(
        len(complements_three) == 4,
        "level-three C3 complement count failed",
    )
    complement_three = subgroup_conjugation_action(
        group_three,
        complements_three,
        3,
    )
    check_action(
        group_three,
        complement_three,
        3,
        "level-three complement",
    )
    require(
        len(
            equivariant_bijections(
                group_three,
                projective_three,
                complement_three,
            )
        ) == 1,
        "level-three projective/complement identification failed",
    )

    # Level four: reduction has normal V4 kernel and four S3 complements.
    group_four = projective_sl(4)
    identity_four = canonical_matrix((1, 0, 0, 1), 4)
    identity_two = canonical_matrix((1, 0, 0, 1), 2)
    s_four = canonical_matrix(s_integral, 4)
    c_four = canonical_matrix(c_integral, 4)
    t_four = matrix_product(
        matrix_inverse(s_four, 4),
        c_four,
        4,
    )
    require(
        (
            sl_size(4),
            len(group_four),
            len(generated_subgroup((s_four, c_four), 4)),
        ) == (48, 24, 24),
        "level-four group sizes failed",
    )
    require(
        (
            element_order(s_four, 4),
            element_order(c_four, 4),
            element_order(t_four, 4),
        ) == (2, 3, 4),
        "level-four triangle orders failed",
    )
    reduction = {
        element: canonical_matrix(
            tuple(entry % 2 for entry in element),
            2,
        )
        for element in group_four
    }
    quotient_two = set(reduction.values())
    klein_four = frozenset(
        element
        for element in group_four
        if reduction[element] == identity_two
    )
    require(
        len(quotient_two) == 6
        and len(klein_four) == 4
        and subgroup_element_orders(klein_four, 4)
        == Counter({2: 3, 1: 1}),
        "level-four reduction exact sequence failed",
    )
    for left in group_four:
        for right in group_four:
            require(
                reduction[matrix_product(left, right, 4)]
                == matrix_product(reduction[left], reduction[right], 2),
                "level-four reduction homomorphism failed",
            )

    subgroups_four = all_two_generated_subgroups(group_four, 4)
    complements_four = tuple(
        sorted(
            (
                subgroup
                for subgroup in subgroups_four
                if len(subgroup) == 6
                and subgroup.intersection(klein_four)
                == {identity_four}
                and len({reduction[element] for element in subgroup})
                == 6
            ),
            key=lambda subgroup: tuple(sorted(subgroup)),
        )
    )
    require(
        len(complements_four) == 4,
        "level-four S3 complement count failed",
    )
    complement_four = subgroup_conjugation_action(
        group_four,
        complements_four,
        4,
    )
    check_action(
        group_four,
        complement_four,
        4,
        "level-four complement",
    )
    expected_s4_census = Counter(
        {
            (3, 1): 8,
            (2, 1, 1): 6,
            (4,): 6,
            (2, 2): 3,
            (1, 1, 1, 1): 1,
        }
    )
    require(
        len(set(complement_four.values())) == 24
        and action_census(complement_four) == expected_s4_census,
        "level-four complement S4 action failed",
    )

    # The six edges of the four-complement torsor.
    edges, edge_action = subset_action(
        group_four,
        complement_four,
        2,
    )
    check_action(group_four, edge_action, 4, "six-edge")
    expected_edge_census = Counter(
        {
            (2, 2, 1, 1): 9,
            (3, 3): 8,
            (4, 2): 6,
            (1, 1, 1, 1, 1, 1): 1,
        }
    )
    require(
        action_census(edge_action) == expected_edge_census
        and len(action_orbit(edge_action)) == 6,
        "six-edge census failed",
    )

    # The natural projective six-set and the conjugacy class of T.
    line_four = projective_line(4)
    projective_four = point_action(group_four, line_four, 4)
    check_action(group_four, projective_four, 4, "level-four projective")
    expected_projective_census = Counter(
        {
            (3, 3): 8,
            (2, 2, 2): 6,
            (4, 1, 1): 6,
            (2, 2, 1, 1): 3,
            (1, 1, 1, 1, 1, 1): 1,
        }
    )
    require(
        len(line_four) == 6
        and len(set(projective_four.values())) == 24
        and len(action_orbit(projective_four)) == 6
        and action_census(projective_four)
        == expected_projective_census,
        "level-four projective census failed",
    )
    infinity = line_four.index((1, 0))
    infinity_stabilizer = frozenset(
        element
        for element in group_four
        if projective_four[element][infinity] == infinity
    )
    cyclic_four = generated_subgroup((t_four,), 4)
    require(
        infinity_stabilizer == cyclic_four
        and subgroup_element_orders(cyclic_four, 4)
        == Counter({4: 2, 2: 1, 1: 1}),
        "projective stabilizer is not <T> = C4",
    )
    four_cycles = tuple(
        sorted(
            element
            for element in group_four
            if element_order(element, 4) == 4
        )
    )
    require(
        len(four_cycles) == 6,
        "order-four conjugacy class size failed",
    )
    four_cycle_action = element_conjugation_action(
        group_four,
        four_cycles,
        4,
    )
    check_action(
        group_four,
        four_cycle_action,
        4,
        "four-cycle conjugation",
    )
    t_index = four_cycles.index(t_four)
    point_to_four_cycle = orbit_map_from_base(
        group_four,
        projective_four,
        infinity,
        t_index,
        lambda element, base: four_cycles.index(
            conjugate_element(
                element,
                four_cycles[base],
                4,
            )
        ),
    )
    require(
        len(set(point_to_four_cycle)) == 6,
        "projective-to-four-cycle map is not bijective",
    )
    for element in group_four:
        for point in range(6):
            require(
                point_to_four_cycle[
                    projective_four[element][point]
                ]
                == four_cycle_action[element][
                    point_to_four_cycle[point]
                ],
                "projective-to-four-cycle equivariance failed",
            )
    projective_four_cycle_bijections = equivariant_bijections(
        group_four,
        projective_four,
        four_cycle_action,
    )
    require(
        len(projective_four_cycle_bijections) == 2
        and point_to_four_cycle in projective_four_cycle_bijections,
        "projective/four-cycle oriented equivariant bijection failed",
    )
    inverse_point_to_four_cycle = tuple(
        four_cycles.index(
            matrix_inverse(
                four_cycles[cycle_position],
                4,
            )
        )
        for cycle_position in point_to_four_cycle
    )
    require(
        set(projective_four_cycle_bijections)
        == {
            point_to_four_cycle,
            inverse_point_to_four_cycle,
        },
        "the two projective/four-cycle gauges are not T and T inverse",
    )

    # Squaring the six four-cycles produces the common matching quotient.
    square_elements = tuple(
        sorted(
            {
                matrix_product(cycle, cycle, 4)
                for cycle in four_cycles
            }
        )
    )
    require(
        len(square_elements) == 3
        and set(square_elements) == set(klein_four) - {identity_four},
        "four-cycle squares are not the nonzero normal V4 elements",
    )
    square_action = element_conjugation_action(
        group_four,
        square_elements,
        4,
    )
    check_action(group_four, square_action, 4, "square quotient")
    four_cycle_to_square = tuple(
        square_elements.index(
            matrix_product(cycle, cycle, 4)
        )
        for cycle in four_cycles
    )
    require(
        Counter(four_cycle_to_square) == Counter({0: 2, 1: 2, 2: 2}),
        "four-cycle square fibres failed",
    )

    line_two = projective_line(2)
    line_four_reduction = {
        point: canonical_point(
            tuple(coordinate % 2 for coordinate in point),
            2,
        )
        for point in line_four
    }
    projective_blocks = tuple(
        frozenset(
            position
            for position, point in enumerate(line_four)
            if line_four_reduction[point] == reduced_point
        )
        for reduced_point in line_two
    )
    require(
        all(len(block) == 2 for block in projective_blocks),
        "projective reduction fibres failed",
    )
    inverse_point_permutation = tuple(
        point_to_four_cycle.index(cycle_position)
        for cycle_position in inverse_point_to_four_cycle
    )
    require(
        all(
            inverse_point_permutation[position] != position
            for position in range(6)
        )
        and {
            frozenset(
                (
                    position,
                    inverse_point_permutation[position],
                )
            )
            for position in range(6)
        } == set(projective_blocks),
        "T/T-inverse gauge does not swap the mod-two fibres",
    )
    for block in projective_blocks:
        require(
            len(
                {
                    four_cycle_to_square[
                        point_to_four_cycle[position]
                    ]
                    for position in block
                }
            ) == 1,
            "gamma-square is not constant on a projective fibre",
        )
    require(
        {
            frozenset(
                position
                for position in range(6)
                if four_cycle_to_square[
                    point_to_four_cycle[position]
                ] == square
            )
            for square in range(3)
        } == set(projective_blocks),
        "gamma-square fibres do not equal mod-two fibres",
    )

    matchings = (
        frozenset(
            (
                edges.index((0, 1)),
                edges.index((2, 3)),
            )
        ),
        frozenset(
            (
                edges.index((0, 2)),
                edges.index((1, 3)),
            )
        ),
        frozenset(
            (
                edges.index((0, 3)),
                edges.index((1, 2)),
            )
        ),
    )
    matching_action = block_action(
        group_four,
        edge_action,
        matchings,
    )
    projective_block_action = block_action(
        group_four,
        projective_four,
        projective_blocks,
    )
    require(
        len(
            equivariant_bijections(
                group_four,
                square_action,
                matching_action,
            )
        ) == 1
        and len(
            equivariant_bijections(
                group_four,
                projective_block_action,
                matching_action,
            )
        ) == 1,
        "common three-point quotient failed",
    )
    kernel_projective_blocks = frozenset(
        element
        for element in group_four
        if projective_block_action[element] == (0, 1, 2)
    )
    kernel_matching_blocks = frozenset(
        element
        for element in group_four
        if matching_action[element] == (0, 1, 2)
    )
    require(
        kernel_projective_blocks
        == kernel_matching_blocks
        == klein_four,
        "common three-block kernel is not the normal V4",
    )

    # Exact parity comparison.
    require(
        all(
            permutation_sign(projective_four[element])
            == permutation_sign(complement_four[element])
            for element in group_four
        ),
        "projective ambient sign does not equal quartic sign",
    )
    require(
        all(
            permutation_sign(edge_action[element]) == 1
            for element in group_four
        ),
        "six-edge ambient sign is not identically positive",
    )

    # One square/matching block has a common D8 stabilizer.
    square_base = four_cycle_to_square[t_index]
    square_element = square_elements[square_base]
    matching_by_square = {}
    for square_position, square in enumerate(square_elements):
        permutation = complement_four[square]
        moved_pairs = tuple(
            sorted(
                tuple(sorted((position, permutation[position])))
                for position in range(4)
                if position < permutation[position]
            )
        )
        require(
            len(moved_pairs) == 2,
            "normal V4 element does not define a matching",
        )
        matching = frozenset(edges.index(pair) for pair in moved_pairs)
        require(
            matching in matchings,
            "normal V4 matching is not in the matching bank",
        )
        matching_by_square[square_position] = matchings.index(matching)
    for element in group_four:
        for square_position in range(3):
            require(
                matching_by_square[
                    square_action[element][square_position]
                ]
                == matching_action[element][
                    matching_by_square[square_position]
                ],
                "square-to-matching equivariance failed",
            )
    matching_base = matching_by_square[square_base]
    common_dihedral = frozenset(
        element
        for element in group_four
        if conjugate_element(
            element,
            square_element,
            4,
        ) == square_element
    )
    matching_stabilizer = frozenset(
        element
        for element in group_four
        if matching_action[element][matching_base] == matching_base
    )
    require(
        common_dihedral == matching_stabilizer
        and len(common_dihedral) == 8
        and subgroup_element_orders(common_dihedral, 4)
        == Counter({2: 5, 4: 2, 1: 1}),
        "common block stabilizer is not D8",
    )

    edge_base = min(matchings[matching_base])
    edge_stabilizer = frozenset(
        element
        for element in group_four
        if edge_action[element][edge_base] == edge_base
    )
    require(
        edge_stabilizer.issubset(common_dihedral)
        and subgroup_element_orders(edge_stabilizer, 4)
        == Counter({2: 3, 1: 1})
        and edge_stabilizer != klein_four,
        "edge stabilizer is not the nonnormal V4",
    )
    index_two_subgroups = {
        subgroup
        for subgroup in subgroups_four
        if len(subgroup) == 4
        and subgroup.issubset(common_dihedral)
    }
    require(
        index_two_subgroups
        == {klein_four, cyclic_four, edge_stabilizer},
        "D8 does not have exactly the three typed index-two subgroups",
    )
    common_center = frozenset(
        element
        for element in common_dihedral
        if all(
            matrix_product(element, other, 4)
            == matrix_product(other, element, 4)
            for other in common_dihedral
        )
    )
    commutator_subgroup = generated_subgroup(
        tuple(
            commutator(left, right, 4)
            for left in common_dihedral
            for right in common_dihedral
        ),
        4,
    )
    triple_kernel = (
        klein_four
        .intersection(cyclic_four)
        .intersection(edge_stabilizer)
    )
    pairwise_intersections = (
        klein_four.intersection(cyclic_four),
        cyclic_four.intersection(edge_stabilizer),
        edge_stabilizer.intersection(klein_four),
    )
    require(
        common_center
        == commutator_subgroup
        == triple_kernel
        == generated_subgroup(
            (
                matrix_product(t_four, t_four, 4),
            ),
            4,
        )
        and all(
            intersection == triple_kernel
            for intersection in pairwise_intersections
        ),
        "D8 center/commutator/character-kernel triangle failed",
    )
    pairwise_joins = (
        generated_subgroup(
            tuple(klein_four) + tuple(cyclic_four),
            4,
        ),
        generated_subgroup(
            tuple(cyclic_four) + tuple(edge_stabilizer),
            4,
        ),
        generated_subgroup(
            tuple(edge_stabilizer) + tuple(klein_four),
            4,
        ),
    )
    require(
        all(join == common_dihedral for join in pairwise_joins),
        "two D8 index-two kernels do not generate D8",
    )

    character_kernel = character_with_kernel(
        common_dihedral,
        klein_four,
    )
    character_projective = character_with_kernel(
        common_dihedral,
        cyclic_four,
    )
    character_edge = character_with_kernel(
        common_dihedral,
        edge_stabilizer,
    )
    check_character(
        common_dihedral,
        character_kernel,
        4,
        "normal-V4 character",
    )
    check_character(
        common_dihedral,
        character_projective,
        4,
        "C4 character",
    )
    check_character(
        common_dihedral,
        character_edge,
        4,
        "edge-V4 character",
    )
    for element in common_dihedral:
        require(
            character_kernel[element]
            * character_projective[element]
            == character_edge[element],
            "D8 character product chi_K chi_P = chi_E failed",
        )
        require(
            character_projective[element]
            * character_edge[element]
            == character_kernel[element],
            "D8 character product chi_P chi_E = chi_K failed",
        )
        require(
            character_edge[element]
            * character_kernel[element]
            == character_projective[element],
            "D8 character product chi_E chi_K = chi_P failed",
        )
    character_triples = Counter(
        (
            character_kernel[element],
            character_projective[element],
            character_edge[element],
        )
        for element in common_dihedral
    )
    require(
        character_triples
        == Counter(
            {
                (1, 1, 1): 2,
                (1, -1, -1): 2,
                (-1, 1, -1): 2,
                (-1, -1, 1): 2,
            }
        ),
        "D8 character-triple census failed",
    )

    require(
        len(
            equivariant_bijections(
                group_four,
                projective_four,
                edge_action,
            )
        ) == 0,
        "the two six-actions were incorrectly identified",
    )
    infinity_orders = sorted(
        element_order(element, 4)
        for element in infinity_stabilizer
    )
    edge_orders = sorted(
        element_order(element, 4)
        for element in edge_stabilizer
    )
    require(
        infinity_orders == [1, 2, 4, 4]
        and edge_orders == [1, 2, 2, 2],
        "six-action stabilizer hostile failed",
    )

    # Class-by-class comparison of the two six-actions.
    class_rows = {}
    for element in group_four:
        abstract_type = cycle_type(complement_four[element])
        row = (
            cycle_type(projective_four[element]),
            cycle_type(edge_action[element]),
        )
        class_rows.setdefault(abstract_type, set()).add(row)
    require(
        all(len(rows) == 1 for rows in class_rows.values()),
        "six-action types are not constant on quartic classes",
    )

    print("MODULAR LEVEL THREE/FOUR SIX-LIFT AUDIT - exact")
    print(
        "group sizes: "
        f"SL3={sl_size(3)} PSL3={len(group_three)}; "
        f"SL4={sl_size(4)} PSL4={len(group_four)}"
    )
    print(
        "triangle orders: "
        f"level3={element_order(s_three,3)},{element_order(c_three,3)},{element_order(t_three,3)}; "
        f"level4={element_order(s_four,4)},{element_order(c_four,4)},{element_order(t_four,4)}"
    )
    print(
        "level3: "
        f"P1={len(line_three)} V4={len(klein_three)} "
        f"C3_complements={len(complements_three)}"
    )
    print(
        "level3 census: "
        f"{format_census(action_census(projective_three))}"
    )
    print(
        "level4 reduction: "
        f"quotient={len(quotient_two)} kernel={len(klein_four)} "
        f"kernel_orders={sorted(element_order(x,4) for x in klein_four)}"
    )
    print(
        "level4 complements: "
        f"count={len(complements_four)} "
        f"census={format_census(action_census(complement_four))}"
    )
    print(
        "P1(Z/4): "
        f"points={len(line_four)} four_cycles={len(four_cycles)} "
        "equivariant_bijections=2 oriented_base=T stabilizer=C4"
    )
    print(
        "P1 census: "
        f"{format_census(action_census(projective_four))}"
    )
    print(
        "edge census: "
        f"{format_census(action_census(edge_action))}"
    )
    print(
        "gamma-square quotient: "
        f"image={len(square_elements)} fibre_sizes={sorted(Counter(four_cycle_to_square).values())} "
        "matching_equivariant=1"
    )
    print(
        "sign census: "
        "projective_equals_quartic=24/24 edge_positive=24/24"
    )
    print(
        "six stabilizers: "
        f"projective_orders={infinity_orders} edge_orders={edge_orders}"
    )
    print(
        "common block stabilizer: "
        f"size={len(common_dihedral)} "
        f"orders={sorted(element_order(x,4) for x in common_dihedral)}"
    )
    print(
        "D8 index-two bank: "
        f"count={len(index_two_subgroups)} "
        "kernels=normal_V4,C4,edge_V4"
    )
    print(
        "D8 characters: "
        "chi_K*chi_P=chi_E chi_P*chi_E=chi_K chi_E*chi_K=chi_P"
    )
    print(
        "D8 fixed-field core: "
        "center=commutator=triple_kernel=<T^2> size=2; "
        "pair_intersections=2,2,2 pair_joins=8,8,8"
    )
    print(
        "D8 character triples: "
        "(+++):2 (+--):2 (-+-):2 (--+):2"
    )
    print(
        "six actions: equivariant=0; "
        "three-block quotients equivariant=1; block_kernel=V4"
    )
    print("PASS")


if __name__ == "__main__":
    main()
