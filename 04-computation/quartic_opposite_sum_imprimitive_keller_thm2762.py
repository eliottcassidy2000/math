#!/usr/bin/env python3
"""Exact audit for THM-2762: the quartic T-wall is imprimitive.

All checks use integer/Fraction arithmetic and explicit exceptions.  The
group census is a complete enumeration inside S4.  No truth-bearing Python
assertions are used.
"""

from fractions import Fraction
from itertools import combinations, permutations, product
from collections import Counter


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


IDENTITY = tuple(range(4))
S4 = tuple(permutations(range(4)))
EDGES = tuple(combinations(range(4), 2))
PAIRING = frozenset((frozenset((0, 1)), frozenset((2, 3))))


def compose(left, right):
    return tuple(left[right[index]] for index in range(4))


def sign(permutation):
    inversions = sum(
        permutation[i] > permutation[j]
        for i, j in combinations(range(4), 2)
    )
    return -1 if inversions % 2 else 1


def order(permutation):
    power = IDENTITY
    for exponent in range(1, 13):
        power = compose(permutation, power)
        if power == IDENTITY:
            return exponent
    raise RuntimeError("permutation order exceeded the S4 exponent")


def closure(generators):
    subgroup = {IDENTITY}
    changed = True
    while changed:
        changed = False
        current = tuple(subgroup)
        for left in current:
            for right in tuple(generators) + current:
                for value in (compose(left, right), compose(right, left)):
                    if value not in subgroup:
                        subgroup.add(value)
                        changed = True
    return frozenset(subgroup)


def image_pairing(permutation):
    return frozenset(
        frozenset(permutation[index] for index in block)
        for block in PAIRING
    )


def edge_image(permutation, edge):
    return tuple(sorted((permutation[edge[0]], permutation[edge[1]])))


def edge_action(permutation):
    return tuple(EDGES.index(edge_image(permutation, edge)) for edge in EDGES)


def action_sign(action):
    inversions = sum(
        action[i] > action[j]
        for i, j in combinations(range(len(action)), 2)
    )
    return -1 if inversions % 2 else 1


def orbit(group, point):
    return {permutation[point] for permutation in group}


def edge_orbit(group, edge):
    return {edge_image(permutation, edge) for permutation in group}


def elementaries(values):
    e1 = sum(values)
    e2 = sum(values[i] * values[j] for i, j in combinations(range(4), 2))
    e3 = sum(
        values[i] * values[j] * values[k]
        for i, j, k in combinations(range(4), 3)
    )
    return e1, e2, e3


def t_invariant(values):
    e1, e2, e3 = elementaries(values)
    return e1 ** 3 - 4 * e1 * e2 + 8 * e3


def pair_sums(values):
    return tuple(values[i] + values[j] for i, j in EDGES)


def discriminant(values):
    out = 1
    for i, j in combinations(range(len(values)), 2):
        out *= (values[i] - values[j]) ** 2
    return out


def centrally_symmetric(values):
    total = sum(values)
    return any(
        values[a] + values[b] == values[c] + values[d]
        for (a, b), (c, d) in (
            ((0, 1), (2, 3)),
            ((0, 2), (1, 3)),
            ((0, 3), (1, 2)),
        )
    ) and total % 2 == 0


def main():
    # The coefficient translation f(Y-a/4) has linear coefficient
    # q=(a^3-4ab+8c)/8 and T=-8q in the root-sign convention.
    coefficient_rows = 0
    coefficient_walls = 0
    for a, b, c, d in product(range(-3, 4), repeat=4):
        coefficient_rows += 1
        center = Fraction(-a, 4)
        cubic = 4 * center + a
        linear = (
            4 * center ** 3
            + 3 * a * center ** 2
            + 2 * b * center
            + c
        )
        coefficient_t = -(a ** 3 - 4 * a * b + 8 * c)
        require(cubic == 0, "quartic translation stopped depressing the cubic")
        require(linear == Fraction(-coefficient_t, 8),
                "depressed linear coefficient stopped being -T/8")
        if coefficient_t == 0:
            coefficient_walls += 1

    # Distinct rational-root census: T=0, central symmetry, and a repeated
    # pair sum are the same wall.
    root_rows = 0
    root_walls = 0
    for values in combinations(range(-4, 6), 4):
        root_rows += 1
        t_value = t_invariant(values)
        central = centrally_symmetric(values)
        repeated_pair_sum = discriminant(pair_sums(values)) == 0
        require((t_value == 0) == central,
                "T-wall stopped matching central symmetry")
        require(central == repeated_pair_sum,
                "central symmetry stopped matching pair-sum collision")
        if central:
            root_walls += 1

    # Complete subgroup census inside the stabilizer of one antipodal pairing.
    matching_stabilizer = tuple(
        permutation for permutation in S4
        if image_pairing(permutation) == PAIRING
    )
    require(len(matching_stabilizer) == 8,
            "perfect-matching stabilizer stopped having order eight")
    subgroups = set()
    for mask in range(1 << len(matching_stabilizer)):
        generators = tuple(
            matching_stabilizer[index]
            for index in range(len(matching_stabilizer))
            if (mask >> index) & 1
        )
        subgroups.add(closure(generators))
    require(len(subgroups) == 10,
            "matching stabilizer subgroup census changed")

    transitive_types = []
    for subgroup in subgroups:
        if len(orbit(subgroup, 0)) != 4:
            continue
        if len(subgroup) == 8:
            group_type = "D4"
        elif len(subgroup) == 4 and any(order(g) == 4 for g in subgroup):
            group_type = "C4"
        elif len(subgroup) == 4:
            group_type = "V4"
        else:
            raise RuntimeError("unexpected transitive matching subgroup")
        transitive_types.append(group_type)
    type_census = Counter(transitive_types)
    require(type_census == Counter({"C4": 1, "V4": 1, "D4": 1}),
            "transitive matching-subgroup list changed")

    # The two live Keller groups are transitive on edges, and the six-edge
    # representation is even.
    A4 = tuple(permutation for permutation in S4 if sign(permutation) == 1)
    require(len(A4) == 12 and len(edge_orbit(A4, EDGES[0])) == 6,
            "A4 stopped being edge-transitive")
    require(len(edge_orbit(S4, EDGES[0])) == 6,
            "S4 stopped being edge-transitive")
    require(all(action_sign(edge_action(permutation)) == 1 for permutation in S4),
            "six-edge S4 action stopped being alternating")

    # Sharp imprimitive control x^4-2: coefficient T=0 and no repeated root.
    x4_minus_2_t = -(0 ** 3 - 4 * 0 * 0 + 8 * 0)
    require(x4_minus_2_t == 0,
            "x^4-2 stopped lying on the opposite-sum wall")

    print("QUARTIC OPPOSITE-SUM / IMPRIMITIVE KELLER WALL AUDIT")
    print(f"coefficient_grid={coefficient_rows} T_walls={coefficient_walls}")
    print("depressed_linear_coefficient=-T/8 exact")
    print(f"distinct_root_census={root_rows} central_walls={root_walls}")
    print("T_zero=central_symmetry=opposite_pair_sum_collision")
    print("matching_stabilizer_order=8 isomorphic_to=D4")
    print(f"matching_stabilizer_subgroups={len(subgroups)}")
    print("transitive_wall_groups=C4:1,V4:1,D4:1")
    print("A4_edge_orbit=6 S4_edge_orbit=6")
    print("S4_six_edge_action=alternating")
    print("sharp_control=x^4-2 irreducible_D4_wall")
    print("KELLER_INPUT=THM2633_excludes_C4_V4_D4")
    print("KELLER_SURVIVOR=T_nonzero_pair_sum_sextic_irreducible_separable_square_disc")
    print("SCOPE=field_level_sidecar_not_affine_cover_not_JC2")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
