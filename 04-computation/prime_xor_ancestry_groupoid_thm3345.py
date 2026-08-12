#!/usr/bin/env python3
"""Exact ancestry-path audit for THM-3345.

Universe
--------
* the complete Berggren tree through hypotenuse 1,105;
* the four vertices and six edges of the first spine-selected affine V4 fibre;
* all 12 directed prime-XOR arrows, all 24 ordered affine squares, and all
  60 directed simple closed walks of lengths 2, 3, and 4 on that K4.

All arithmetic and free-group reductions are exact.  Validation uses explicit
exceptions so normal and optimized Python execute the same checks.
"""

from itertools import permutations


U = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
A = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
D = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
MATRICES = {"U": U, "A": A, "D": D}
ROOT = (3, 4, 5)
LIMIT = 1105

VERTEX_ORDER = ("000", "100", "010", "001")
EXPECTED = {
    "000": ((943, 576, 1105), "DAUD"),
    "100": ((817, 744, 1105), "UDUA"),
    "010": ((1073, 264, 1105), "DADDD"),
    "001": ((47, 1104, 1105), "U" * 22),
}
DIRECTIONS = (
    (5, "100", (5, 221), ((47, 1073), (817, 943))),
    (13, "010", (13, 85), ((47, 817), (943, 1073))),
    (17, "001", (17, 65), ((47, 943), (817, 1073))),
)
EXPECTED_COSTS = {
    5: ((8, 0), (27, 17)),
    13: ((5, 1), (24, 18)),
    17: ((26, 18), (9, 1)),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def matrix_vector(matrix, vector):
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(3))
        for i in range(3)
    )


def apply_word(word):
    value = ROOT
    for letter in word:
        value = matrix_vector(MATRICES[letter], value)
    return value


def build_tree(limit):
    addresses = {ROOT: ""}
    frontier = [ROOT]
    cursor = 0
    while cursor < len(frontier):
        parent = frontier[cursor]
        cursor += 1
        for letter, matrix in MATRICES.items():
            child = matrix_vector(matrix, parent)
            require(child[2] > parent[2],
                    f"hypotenuse failed to increase at {parent}->{child}")
            if child[2] <= limit:
                require(child not in addresses,
                        f"duplicate Berggren address at {child}")
                addresses[child] = addresses[parent] + letter
                frontier.append(child)
    return addresses


def bit_tuple(label):
    return tuple(int(char) for char in label)


def bit_label(bits):
    return "".join(str(bit) for bit in bits)


def quotient_representative(bits):
    """Return the weight-zero-or-one representative modulo 111."""
    require(all(bit in (0, 1) for bit in bits), "nonbinary V4 coordinate")
    if sum(bits) <= 1:
        return bits
    return tuple(1 - bit for bit in bits)


def quotient_add(left, right):
    bits = tuple(x ^ y for x, y in zip(bit_tuple(left), bit_tuple(right)))
    return bit_label(quotient_representative(bits))


def inverse(tokens):
    return tuple((letter, -sign) for letter, sign in reversed(tokens))


def reduce_free(tokens):
    stack = []
    for token in tokens:
        if stack and stack[-1] == (token[0], -token[1]):
            stack.pop()
        else:
            stack.append(token)
    return tuple(stack)


def word_tokens(word):
    return tuple((letter, 1) for letter in word)


def ancestry_path(source, target):
    source_word = EXPECTED[source][1]
    target_word = EXPECTED[target][1]
    return reduce_free(inverse(word_tokens(source_word)) + word_tokens(target_word))


def common_prefix_length(left, right):
    length = 0
    for x, y in zip(left, right):
        if x != y:
            break
        length += 1
    return length


def compress_tokens(tokens):
    if not tokens:
        return "1"
    pieces = []
    cursor = 0
    while cursor < len(tokens):
        token = tokens[cursor]
        end = cursor + 1
        while end < len(tokens) and tokens[end] == token:
            end += 1
        exponent = token[1] * (end - cursor)
        pieces.append(token[0] if exponent == 1 else f"{token[0]}^{exponent}")
        cursor = end
    return " ".join(pieces)


def matching_pairs(direction):
    seen = set()
    pairs = []
    for source in VERTEX_ORDER:
        if source in seen:
            continue
        target = quotient_add(source, direction)
        require(target in VERTEX_ORDER and target != source,
                "direction failed to act freely on V4")
        seen.update((source, target))
        pairs.append((source, target))
    require(len(pairs) == 2 and len(seen) == 4,
            "direction did not produce a perfect matching")
    return tuple(pairs)


def path_signature(path):
    reverse = inverse(path)
    return min(path, reverse)


def compose_paths(vertices):
    tokens = ()
    for source, target in zip(vertices, vertices[1:]):
        tokens = reduce_free(tokens + ancestry_path(source, target))
    return tokens


def main():
    addresses = build_tree(LIMIT)
    fibre = {triple: word for triple, word in addresses.items()
             if triple[2] == LIMIT}
    require(len(fibre) == 4, "c=1105 fibre did not have four vertices")

    for bit, (triple, word) in EXPECTED.items():
        require(apply_word(word) == triple,
                f"word evaluation failed at bit {bit}")
        require(addresses.get(triple) == word,
                f"bounded-tree address failed at bit {bit}")
        require(triple[0] ** 2 + triple[1] ** 2 == triple[2] ** 2,
                f"Pythagorean identity failed at bit {bit}")

    require(set(fibre) == {row[0] for row in EXPECTED.values()},
            "unexpected c=1105 vertex")
    require(addresses[(13, 84, 85)] == "U" * 5
            and addresses[(77, 36, 85)] == "AD",
            "c=85 positive ancestry control failed")
    c85_cost = len(ancestry_path_for_words("U" * 5, "AD"))
    require(c85_cost == 7, "c=85 path cost failed")

    rows = []
    matching_summary = {}
    for prime, direction, folded_weight, expected_legs in DIRECTIONS:
        pairs = matching_pairs(direction)
        odd_leg_pairs = tuple(sorted(tuple(sorted((EXPECTED[left][0][0],
                                                   EXPECTED[right][0][0])))
                                     for left, right in pairs))
        require(odd_leg_pairs == tuple(sorted(expected_legs)),
                f"prime-{prime} odd-leg matching failed")

        costs = []
        signatures = []
        cross_depth = False
        for left, right in pairs:
            left_word = EXPECTED[left][1]
            right_word = EXPECTED[right][1]
            prefix = common_prefix_length(left_word, right_word)
            path = ancestry_path(left, right)
            distance = len(path)
            depth_jump = abs(len(left_word) - len(right_word))
            require(distance == len(left_word) + len(right_word) - 2 * prefix,
                    "tree-distance formula failed")
            require(ancestry_path(right, left) == inverse(path),
                    "path reversal failed")
            require(not left_word.startswith(right_word)
                    and not right_word.startswith(left_word),
                    "equal-hypotenuse endpoints became ancestor-related")
            cross_depth = cross_depth or depth_jump != 0
            costs.append((distance, depth_jump))
            signatures.append(path_signature(path))
            rows.append((prime, folded_weight, left, right, prefix,
                         distance, depth_jump, compress_tokens(path)))

        require(tuple(costs) == EXPECTED_COSTS[prime],
                f"prime-{prime} ancestry costs changed")
        require(costs[0][0] != costs[1][0],
                f"prime-{prime} colour unexpectedly fixed path cost")
        require(signatures[0] != signatures[1],
                f"prime-{prime} colour unexpectedly fixed a path label")
        require(cross_depth,
                f"prime-{prime} matching escaped the depth obstruction")
        matching_summary[prime] = tuple(costs)

    # The source-dependent arrows are involutive and every affine square is flat.
    involution_checks = 0
    square_checks = 0
    direction_labels = tuple(direction for _, direction, _, _ in DIRECTIONS)
    for source in VERTEX_ORDER:
        for direction in direction_labels:
            target = quotient_add(source, direction)
            require(reduce_free(ancestry_path(source, target)
                                + ancestry_path(target, source)) == (),
                    "direction involution acquired holonomy")
            involution_checks += 1

        for first in direction_labels:
            for second in direction_labels:
                if first == second:
                    continue
                via_first = quotient_add(source, first)
                target = quotient_add(via_first, second)
                via_second = quotient_add(source, second)
                left_path = reduce_free(ancestry_path(source, via_first)
                                        + ancestry_path(via_first, target))
                right_path = reduce_free(ancestry_path(source, via_second)
                                         + ancestry_path(via_second, target))
                require(left_path == right_path == ancestry_path(source, target),
                        "affine square failed flatness")
                square_checks += 1

    # Every closed walk is null after retaining exact endpoints and tree paths.
    loop_checks = 0
    for length in (2, 3, 4):
        for vertices in permutations(VERTEX_ORDER, length):
            closed = vertices + (vertices[0],)
            require(compose_paths(closed) == (),
                    f"nontrivial ambient-tree holonomy at {closed}")
            loop_checks += 1
    require(loop_checks == 60, "closed-walk census changed")

    # Freezing all arrows at the basepoint does not respect the V4 relation.
    basepoint = "000"
    base_arrows = {
        prime: ancestry_path(basepoint, quotient_add(basepoint, direction))
        for prime, direction, _, _ in DIRECTIONS
    }
    prime_for_direction = {direction: prime for prime, direction, _, _ in DIRECTIONS}
    naive_defects = []
    for first_index in range(len(DIRECTIONS)):
        for second_index in range(first_index + 1, len(DIRECTIONS)):
            first_prime, first_direction, _, _ = DIRECTIONS[first_index]
            second_prime, second_direction, _, _ = DIRECTIONS[second_index]
            target_direction = quotient_add(first_direction, second_direction)
            target_prime = prime_for_direction[target_direction]
            defect = reduce_free(base_arrows[first_prime]
                                 + base_arrows[second_prime]
                                 + inverse(base_arrows[target_prime]))
            require(defect, "frozen-source colour labels accidentally formed V4")
            naive_defects.append((first_prime, second_prime, target_prime,
                                  len(defect), compress_tokens(defect)))

    print("THM-3345 VERIFIED-EXACT")
    print(f"berggren_tree_nodes_through_1105={len(addresses)}")
    print("vertices(bit,depth,word,odd_leg):")
    for bit in VERTEX_ORDER:
        triple, word = EXPECTED[bit]
        shown_word = "U^22" if word == "U" * 22 else word
        print(f"  {bit}: depth={len(word)}, word={shown_word}, odd_leg={triple[0]}")
    print(f"c85_positive_control_distance={c85_cost}")
    print("prime_matchings(weight,source,target,lcp,distance,abs_depth_jump,path):")
    for row in rows:
        print(f"  p={row[0]}, K={row[1]}, {row[2]}->{row[3]}, lcp={row[4]}, "
              f"distance={row[5]}, jump={row[6]}, path={row[7]}")
    print(f"matching_costs={matching_summary}")
    print(f"directed_involutions={involution_checks}; affine_square_checks={square_checks}")
    print(f"closed_walk_holonomy_checks={loop_checks}; all_identity=True")
    print("frozen_basepoint_V4_defects(first,second,target,length,path):")
    for row in naive_defects:
        print(f"  {row}")
    print("rooted_tree_automorphism_controls=3/3 blocked_by_depth")
    print("colour_only_path_controls=3/3 blocked_by_distinct_labels")
    print("ALL CHECKS PASSED")


def ancestry_path_for_words(source_word, target_word):
    return reduce_free(inverse(word_tokens(source_word)) + word_tokens(target_word))


if __name__ == "__main__":
    main()
