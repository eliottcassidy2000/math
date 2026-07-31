#!/usr/bin/env python3
"""Exact signed-wreath/Kummer audit for THM-2766.

No floating point and no truth-bearing Python assertions are used.
"""

from itertools import combinations, permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


S3 = tuple(permutations(range(3)))
MASKS = tuple(range(8))
OMEGA = tuple(mask for mask in MASKS if mask.bit_count() % 2 == 0)


def permute_mask(permutation, mask):
    out = 0
    for index in range(3):
        if (mask >> index) & 1:
            out |= 1 << permutation[index]
    return out


def signed_compose(left, right):
    """Compose (mask, permutation) actions, left after right."""
    left_mask, left_permutation = left
    right_mask, right_permutation = right
    permutation = tuple(
        left_permutation[right_permutation[index]] for index in range(3)
    )
    mask = left_mask ^ permute_mask(left_permutation, right_mask)
    return mask, permutation


def sign_of_permutation(permutation):
    inversions = sum(
        permutation[i] > permutation[j]
        for i, j in combinations(range(len(permutation)), 2)
    )
    return -1 if inversions % 2 else 1


def six_root_action(element):
    mask, permutation = element
    action = []
    for index in range(3):
        for sign_bit in range(2):
            moved_index = permutation[index]
            moved_sign = sign_bit ^ ((mask >> moved_index) & 1)
            action.append(2 * moved_index + moved_sign)
    return tuple(action)


def omega_action(element):
    mask, permutation = element
    return tuple(
        OMEGA.index(mask ^ permute_mask(permutation, state))
        for state in OMEGA
    )


def span(vectors):
    values = {0}
    for vector in vectors:
        values |= {value ^ vector for value in tuple(values)}
    return frozenset(values)


def dot(left, right):
    return (left & right).bit_count() % 2


def annihilator(relations):
    return frozenset(
        vector for vector in MASKS
        if all(dot(vector, relation) == 0 for relation in relations)
    )


def discriminant(values):
    out = 1
    for i, j in combinations(range(len(values)), 2):
        out *= (values[i] - values[j]) ** 2
    return out


def quartic_t(values):
    e1 = sum(values)
    e2 = sum(values[i] * values[j] for i, j in combinations(range(4), 2))
    e3 = sum(
        values[i] * values[j] * values[k]
        for i, j, k in combinations(range(4), 3)
    )
    return e1 ** 3 - 4 * e1 * e2 + 8 * e3


def opposite_differences(values):
    return (
        values[0] + values[1] - values[2] - values[3],
        values[0] + values[2] - values[1] - values[3],
        values[0] + values[3] - values[1] - values[2],
    )


def main():
    signed_group = tuple((mask, permutation) for mask in MASKS for permutation in S3)
    require(len(set(signed_group)) == 48, "B3 stopped having order 48")

    # Six-root parity is exactly mask parity; block permutations occur twice.
    for element in signed_group:
        mask, _ = element
        require(
            sign_of_permutation(six_root_action(element))
            == (-1 if mask.bit_count() % 2 else 1),
            "six-root sign stopped equalling flip parity",
        )

    even_group = tuple(element for element in signed_group if element[0] in OMEGA)
    require(len(even_group) == 24, "W(D3) stopped having order 24")
    even_actions = {omega_action(element) for element in even_group}
    require(len(even_actions) == 24,
            "four-state action of W(D3) stopped being faithful/S4")

    three_cycle = (1, 2, 0)
    C3 = {
        (0, 1, 2),
        three_cycle,
        tuple(three_cycle[three_cycle[index]] for index in range(3)),
    }
    cyclic_preimage = tuple(
        element for element in even_group if element[1] in C3
    )
    require(len(cyclic_preimage) == 12,
            "cyclic-base inverse image stopped having order 12")
    cyclic_actions = {omega_action(element) for element in cyclic_preimage}
    require(len(cyclic_actions) == 12
            and all(sign_of_permutation(action) == 1 for action in cyclic_actions),
            "cyclic-base inverse image stopped being A4")

    # All subspaces of F2^3 and the exact relation/annihilator dictionary.
    subspaces = {
        span(vector for vector in MASKS if (choice >> vector) & 1)
        for choice in range(1 << 8)
    }
    require(len(subspaces) == 16, "F2^3 subspace census changed")
    product_relation_space = span((0b111,))
    product_annihilator = annihilator(product_relation_space)
    require(product_annihilator == frozenset(OMEGA),
            "product relation stopped cutting out the even sign plane")
    require(len(product_annihilator) == 4,
            "rank-two product Kummer kernel stopped being V4")

    # Exact pullback discriminant formula on all 120 nonzero distinct triples.
    triple_rows = 0
    for values in combinations(range(1, 11), 3):
        triple_rows += 1
        pullback_roots = []
        # Use square integer t-values to test the literal six roots.
        square_values = tuple(value * value for value in values)
        for value in values:
            pullback_roots.extend((value, -value))
        left = discriminant(tuple(pullback_roots))
        right = (
            4 ** 3
            * square_values[0] * square_values[1] * square_values[2]
            * discriminant(square_values) ** 2
        )
        require(left == right, "quadratic-cubic pullback discriminant changed")
    require(triple_rows == 120, "triple census changed")

    # Quartic specialization t_i=d_i^2/4 and product=(T/8)^2.
    quartic_rows = 0
    nonzero_rows = 0
    for values in combinations(range(-4, 6), 4):
        quartic_rows += 1
        differences = opposite_differences(values)
        t_value = quartic_t(values)
        require(differences[0] * differences[1] * differences[2] == t_value,
                "quartic opposite product stopped being T")
        require(
            differences[0] ** 2 * differences[1] ** 2 * differences[2] ** 2
            == t_value ** 2,
            "quartic pullback product-square relation changed",
        )
        if t_value != 0:
            nonzero_rows += 1
    require(quartic_rows == 210 and nonzero_rows == 160,
            "quartic specialization census changed")

    print("QUADRATIC-CUBIC PULLBACK / EVEN-SIGN WEYL D3 AUDIT")
    print("ambient_signed_wreath_order=48")
    print("six_root_sign=sign_flip_parity")
    print("even_sign_subgroup_order=24")
    print("even_sign_four_state_action=S4=W(D3)=W(A3)")
    print("cyclic_base_preimage_order=12 action=A4")
    print("F2^3_subspaces=16")
    print("product_relation_111_annihilator=000,011,101,110=V4")
    print(f"pullback_triple_census={triple_rows} discriminant_formula=exact")
    print("disc_RV2=4^3*product_roots*disc_R^2")
    print(f"quartic_census={quartic_rows} T_nonzero={nonzero_rows}")
    print("quartic_product_ti=(T/8)^2")
    print("rank2+C3=A4 rank2+S3=S4")
    print("SCOPE=semidirect_binary_over_ternary_not_free_product_not_JC2")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
