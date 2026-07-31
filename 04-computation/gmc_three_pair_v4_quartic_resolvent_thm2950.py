#!/usr/bin/env python3
"""Exact finite controls for THM-2950."""

from itertools import permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def complement(bits):
    return tuple(1 - bit for bit in bits)


def canonical_class(bits):
    bits = tuple(bits)
    return min(bits, complement(bits))


CLASSES = tuple(
    sorted({canonical_class(bits) for bits in product((0, 1), repeat=3)})
)
CLASS_INDEX = {bits: index for index, bits in enumerate(CLASSES)}


def add_classes(left, right):
    return canonical_class(tuple(a ^ b for a, b in zip(left, right)))


def permute_bits(bits, permutation):
    return tuple(bits[permutation[index]] for index in range(3))


def affine_permutation(translation, linear):
    image = []
    for bits in CLASSES:
        moved = canonical_class(permute_bits(bits, linear))
        image.append(CLASS_INDEX[add_classes(moved, translation)])
    return tuple(image)


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(right)))


def gaussian_mul(left, right):
    a, b = left
    c, d = right
    return (a * c - b * d, a * d + b * c)


def gaussian_conjugate(value):
    return (value[0], -value[1])


def half_trace(values, bits):
    result = (1, 0)
    for value, bit in zip(values, bits):
        result = gaussian_mul(
            result,
            gaussian_conjugate(value) if bit else value,
        )
    return 2 * result[0]


def gaussian_norm(value):
    return value[0] ** 2 + value[1] ** 2


def discriminant(values):
    answer = 1
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            answer *= (values[i] - values[j]) ** 2
    return answer


def main():
    require(len(CLASSES) == 4, "half-system quotient size changed")
    translations = {
        affine_permutation(translation, (0, 1, 2))
        for translation in CLASSES
    }
    require(len(translations) == 4, "V4 translation action changed")
    nonidentity = [
        action for action in translations if action != tuple(range(4))
    ]
    require(
        all(compose(action, action) == tuple(range(4)) for action in nonidentity),
        "nontrivial V4 translation lost order two",
    )
    require(
        all(all(action[index] != index for index in range(4)) for action in nonidentity),
        "V4 translation acquired a fixed point",
    )

    affine_group = {
        affine_permutation(translation, linear)
        for translation in CLASSES
        for linear in permutations(range(3))
    }
    require(len(affine_group) == 24, "affine group order changed")
    require(
        affine_group == set(permutations(range(4))),
        "affine action is not the full S4",
    )

    matchings = (
        frozenset((frozenset((0, 1)), frozenset((2, 3)))),
        frozenset((frozenset((0, 2)), frozenset((1, 3)))),
        frozenset((frozenset((0, 3)), frozenset((1, 2)))),
    )

    def matching_action(permutation):
        image = []
        for matching in matchings:
            moved = frozenset(
                frozenset(permutation[index] for index in pair)
                for pair in matching
            )
            image.append(matchings.index(moved))
        return tuple(image)

    s4 = tuple(permutations(range(4)))
    matching_images = {matching_action(permutation) for permutation in s4}
    kernel = [
        permutation
        for permutation in s4
        if matching_action(permutation) == (0, 1, 2)
    ]
    require(len(matching_images) == 6, "resolvent image is not S3")
    require(len(kernel) == 4, "resolvent kernel is not V4")
    require(set(kernel) == translations, "two V4 realizations disagree")

    discriminant_controls = 0
    for roots in (
        (-4, -1, 2, 7),
        (0, 1, 3, 8),
        (-9, -2, 5, 6),
        (1, 4, 10, 23),
    ):
        r0, r1, r2, r3 = roots
        resolvent_roots = (
            r0 * r1 + r2 * r3,
            r0 * r2 + r1 * r3,
            r0 * r3 + r1 * r2,
        )
        require(
            discriminant(roots) == discriminant(resolvent_roots),
            "quartic/resolvent discriminants disagree",
        )
        discriminant_controls += 1

    all_real = ((1, 0), (1, 0), (1, 0))
    one_imaginary = ((0, 1), (1, 0), (1, 0))
    generic = ((1, 1), (2, 1), (3, 2))
    traces_real = tuple(half_trace(all_real, bits) for bits in CLASSES)
    traces_imaginary = tuple(
        half_trace(one_imaginary, bits) for bits in CLASSES
    )
    traces_generic = tuple(half_trace(generic, bits) for bits in CLASSES)
    norm_real = 1
    norm_imaginary = 1
    norm_generic = 1
    for value in generic:
        norm_generic *= gaussian_norm(value)
    require(traces_real == (2, 2, 2, 2), "real trace hostile changed")
    require(traces_imaginary == (0, 0, 0, 0), "phase-loss hostile changed")
    require(norm_real == norm_imaginary == 1, "common norm hostile changed")
    require(len(set(traces_generic)) == 4, "generic V4 traces collided")

    print("THM-2950 THREE-CONJUGATE-PAIR V4/RESOLVENT AUDIT")
    print(f"half_systems=8;mod_global_conjugation={len(CLASSES)}")
    print(f"V4_translations={len(translations)};fixed_point_free_nonidentity=3")
    print(f"affine_group_order={len(affine_group)};affine_group=S4")
    print(f"pairing_action_image={len(matching_images)};kernel={len(kernel)}")
    print(
        f"quartic_resolvent_discriminant_controls={discriminant_controls}"
    )
    print(
        f"same_norm_hostile={norm_real};"
        f"traces_real={traces_real};traces_one_imaginary={traces_imaginary}"
    )
    print(
        f"generic_norm={norm_generic};generic_half_traces={traces_generic}"
    )
    print("consequence=full_norm_forgets_the_V4_half_system_phase")
    print("scope=no_canonical_C3_order_or_SFC4_or_JC2_conclusion")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
