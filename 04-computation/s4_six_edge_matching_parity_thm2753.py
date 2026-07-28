#!/usr/bin/env python3
"""Exact S4 six-edge parity erasure and matching restoration census.

The action on two-subsets is faithful but has trivial ambient S6 sign.  Its
unlabelled cycle type identifies transpositions with double transpositions.
The action on the three perfect matchings has kernel V4 and restores the
original S4 sign.  All checks are exact and use explicit exceptions rather
than Python assertions.
"""

from itertools import combinations, permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


VERTICES = tuple(range(4))
IDENTITY = VERTICES
EDGES = tuple(combinations(VERTICES, 2))
MATCHINGS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


def canonical_edge(a, b):
    return tuple(sorted((a, b)))


def compose(left, right):
    """Return left after right."""
    return tuple(left[right[index]] for index in VERTICES)


def inverse(permutation):
    out = [0] * len(permutation)
    for index, image in enumerate(permutation):
        out[image] = index
    return tuple(out)


def induced_on_edges(permutation):
    return tuple(
        EDGES.index(canonical_edge(permutation[a], permutation[b]))
        for a, b in EDGES
    )


def induced_on_matchings(permutation):
    result = []
    canonical_matchings = tuple(
        frozenset(canonical_edge(*edge) for edge in matching)
        for matching in MATCHINGS
    )
    for matching in MATCHINGS:
        image = frozenset(
            canonical_edge(permutation[a], permutation[b])
            for a, b in matching
        )
        result.append(canonical_matchings.index(image))
    return tuple(result)


def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def permutation_sign(permutation):
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def generated_group(generators):
    group = {IDENTITY}
    frontier = [IDENTITY]
    generators = tuple(generators) + tuple(inverse(g) for g in generators)
    while frontier:
        element = frontier.pop()
        for generator in generators:
            product = compose(generator, element)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return frozenset(group)


def matching_image(group):
    return frozenset(induced_on_matchings(element) for element in group)


def format_type(parts):
    counts = {}
    for part in parts:
        counts[part] = counts.get(part, 0) + 1
    return " ".join(
        str(part) if multiplicity == 1 else f"{part}^{multiplicity}"
        for part, multiplicity in sorted(counts.items())
    )


def main():
    s4 = tuple(permutations(VERTICES))
    edge_actions = {g: induced_on_edges(g) for g in s4}
    matching_actions = {g: induced_on_matchings(g) for g in s4}

    require(len(set(edge_actions.values())) == 24,
            "the six-edge action stopped being faithful")
    require(all(permutation_sign(edge_actions[g]) == 1 for g in s4),
            "the exterior two-subset action acquired ambient odd parity")
    require(all(
        permutation_sign(matching_actions[g]) == permutation_sign(g)
        for g in s4
    ), "the matching quotient stopped restoring quartic parity")

    matching_kernel = frozenset(
        g for g in s4 if matching_actions[g] == (0, 1, 2)
    )
    expected_v4 = frozenset((
        IDENTITY,
        (1, 0, 3, 2),
        (2, 3, 0, 1),
        (3, 2, 1, 0),
    ))
    require(matching_kernel == expected_v4,
            "the matching kernel stopped being the normal V4")
    require(len(set(matching_actions.values())) == 6,
            "the matching quotient stopped being S3")

    rows = {}
    for g in s4:
        key = cycle_type(g)
        value = (
            permutation_sign(g),
            cycle_type(edge_actions[g]),
            cycle_type(matching_actions[g]),
        )
        if key in rows:
            require(rows[key][0] == value,
                    "one S4 conjugacy class split in the action table")
            rows[key] = (value, rows[key][1] + 1)
        else:
            rows[key] = (value, 1)

    expected_rows = {
        (1, 1, 1, 1): ((1, (1, 1, 1, 1, 1, 1), (1, 1, 1)), 1),
        (2, 1, 1): ((-1, (2, 2, 1, 1), (2, 1)), 6),
        (2, 2): ((1, (2, 2, 1, 1), (1, 1, 1)), 3),
        (3, 1): ((1, (3, 3), (3,)), 8),
        (4,): ((-1, (4, 2), (2, 1)), 6),
    }
    require(rows == expected_rows, "the complete S4 action table changed")

    transposition = (3, 1, 2, 0)       # (0 3)
    double_translation = (1, 0, 3, 2)  # (0 1)(2 3)
    ternary = (1, 2, 0, 3)             # (0 1 2)
    require(cycle_type(edge_actions[transposition]) == (2, 2, 1, 1)
            == cycle_type(edge_actions[double_translation]),
            "the sharp six-edge cycle collision disappeared")
    require(cycle_type(matching_actions[transposition]) == (2, 1)
            and matching_actions[double_translation] == (0, 1, 2),
            "the matching quotient stopped separating the binary classes")

    a4_group = generated_group((double_translation, ternary))
    s4_group = generated_group((transposition, ternary))
    require(len(a4_group) == 12
            and all(permutation_sign(g) == 1 for g in a4_group)
            and len(matching_image(a4_group)) == 3,
            "the marked A4/C3 branch changed")
    require(len(s4_group) == 24 and len(matching_image(s4_group)) == 6,
            "the marked S4/S3 branch changed")

    witness_edge = EDGES.index((0, 2))
    transposition_image = EDGES[edge_actions[transposition][witness_edge]]
    double_image = EDGES[edge_actions[double_translation][witness_edge]]
    require(transposition_image == (2, 3)
            and double_image == (1, 3)
            and transposition_image != double_image,
            "the labelled-edge hostile stopped separating the generators")
    transposition_joint_type = cycle_type(induced_on_edges(
        compose(transposition, ternary)
    ))
    double_joint_type = cycle_type(induced_on_edges(
        compose(double_translation, ternary)
    ))
    require(transposition_joint_type == (4, 2)
            and double_joint_type == (3, 3),
            "the labelled joint-word hostile changed")

    edge_type_to_matching_types = {}
    for g in s4:
        edge_type_to_matching_types.setdefault(
            cycle_type(edge_actions[g]), set()
        ).add(cycle_type(matching_actions[g]))
    require(edge_type_to_matching_types[(2, 2, 1, 1)]
            == {(2, 1), (1, 1, 1)},
            "the no-factor-through witness changed")

    print("S4 SIX-EDGE PARITY ERASURE / MATCHING RESTORATION AUDIT")
    print("class size sign edge_cycle matching_cycle")
    for source_type in (
        (1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,)
    ):
        (sign, edge_type, match_type), size = rows[source_type]
        print(
            f"{format_type(source_type):7s} {size:2d} {sign:+d} "
            f"{format_type(edge_type):7s} {format_type(match_type)}"
        )
    print("edge_action=faithful image_subset=A6 ambient_sign=trivial")
    print("matching_action=quotient_S3 kernel=V4 sign_pullback=quartic_sign")
    print("binary_collision=edge_cycle(03)=edge_cycle((01)(23))=1^2 2^2")
    print("matching_separation=reflection_vs_identity")
    print("marked_images=double+C3:A4->C3 transposition+C3:S4->S3")
    print(
        "labelled_hostile=edge02 maps under (03) to "
        f"{transposition_image}, under (01)(23) to {double_image}"
    )
    print(
        "joint_hostile=edge_cycle((03)*C3)=2 4 "
        "edge_cycle((01)(23)*C3)=3^2"
    )
    print("factor_boundary=edge_cycle_class_does_not_determine_matching_class")
    print("SCOPE: finite permutation information-loss theorem; no Keller or LRC exclusion")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
