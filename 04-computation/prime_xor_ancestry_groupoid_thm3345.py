#!/usr/bin/env python3
"""Exact ancestry-path audit for THM-3345.

Universe
--------
* the complete Berggren tree through hypotenuse 99,905;
* the four vertices and six edges of the first spine-selected affine V4 fibre;
* the eight vertices and 28 edges of the first spine-selected affine F_2^3
  fibre;
* all 12 directed prime-XOR arrows, all 24 ordered affine squares, and all
  60 directed closed walks with 2, 3, or 4 distinct nonterminal vertices on
  that K4;
* all 56 directed arrows and 336 ordered affine squares on the eight-parent
  scale control;
* the first 2,000 rows of the explicit infinite prime-5 toggle family and five
  simultaneous Boolean-rank CRT controls.

All arithmetic and free-group reductions are exact.  Validation uses explicit
exceptions so normal and optimized Python execute the same checks.
"""

from itertools import permutations
from math import gcd


U = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
A = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
D = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
MATRICES = {"U": U, "A": A, "D": D}
IDENTITY3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
ROOT = (3, 4, 5)
FIRST_LIMIT = 1105
SECOND_LIMIT = 99905
LIMIT = SECOND_LIMIT

VERTEX_ORDER = ("000", "100", "010", "001")
EXPECTED = {
    "000": ((943, 576, 1105), "DAUD"),
    "100": ((817, 744, 1105), "UDUA"),
    "010": ((1073, 264, 1105), "DADDD"),
    "001": ((47, 1104, 1105), "U" * 22),
}
EXPECTED_99905 = {
    "0000": ((96033, 27544, 99905), "DA" + "U" * 8 + "D" * 3),
    "0001": ((67137, 73984, 99905), "DADDDUDA"),
    "0010": ((32193, 94576, 99905), "A" + "U" * 11 + "AUU"),
    "0011": ((70623, 70664, 99905), "AUAAAA"),
    "0100": ((48063, 87584, 99905), "U" * 4 + "A" + "U" * 4 + "DU"),
    "0101": ((99807, 4424, 99905), "U" * 6 + "D" * 22),
    "0110": ((89823, 43736, 99905), "U" * 22 + "AAD"),
    "0111": ((447, 99904, 99905), "U" * 222),
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
SECOND_DIRECTIONS = ("0001", "0010", "0011", "0100", "0101", "0110", "0111")
SECOND_WEIGHTS = {
    "0001": (53, 1885),
    "0010": (29, 3445),
    "0011": (65, 1537),
    "0100": (13, 7685),
    "0101": (145, 689),
    "0110": (265, 377),
    "0111": (5, 19981),
}
SECOND_COSTS = {
    "0001": (17, 17, 31, 203),
    "0010": (28, 14, 28, 238),
    "0011": (19, 23, 225, 41),
    "0100": (24, 36, 40, 228),
    "0101": (41, 19, 237, 31),
    "0110": (38, 230, 26, 34),
    "0111": (235, 33, 43, 17),
}
SECOND_JUMPS = {
    "0001": (5, 9, 17, 197),
    "0010": (2, 2, 14, 194),
    "0011": (7, 7, 211, 3),
    "0100": (2, 20, 10, 216),
    "0101": (15, 3, 207, 19),
    "0110": (12, 214, 4, 22),
    "0111": (209, 17, 13, 5),
}
SECOND_GAUSSIAN_FACTORS = ((2, 1), (3, 2), (5, 2), (7, 2))
INFINITE_FAMILY_BOUND = 2000
CRT_PRIMES = (13, 17, 29, 37, 41)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def matrix_vector(matrix, vector):
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(3))
        for i in range(3)
    )


def matrix_multiply(left, right):
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(3))
              for j in range(3))
        for i in range(3)
    )


def matrix_power(matrix, exponent):
    require(exponent >= 0, "negative matrix exponent")
    result = IDENTITY3
    base = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = matrix_multiply(result, base)
        base = matrix_multiply(base, base)
        remaining //= 2
    return result


def gaussian_multiply(left, right):
    return (left[0] * right[0] - left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


def gaussian_parameters(bit_label_value):
    value = (1, 0)
    for bit, factor in zip(bit_tuple(bit_label_value), SECOND_GAUSSIAN_FACTORS):
        chosen = factor if bit == 0 else (factor[0], -factor[1])
        value = gaussian_multiply(value, chosen)
    return tuple(sorted((abs(value[0]), abs(value[1])), reverse=True))


def euclid_triple(parameters):
    m, n = parameters
    return (m * m - n * n, 2 * m * n, m * m + n * n)


def spine_c(t):
    return 2 * t * (t + 1) + 1


def roots_of_spine_c(prime):
    return tuple(t for t in range(prime) if spine_c(t) % prime == 0)


def crt_pair(residue, modulus, new_residue, new_modulus):
    step = ((new_residue - residue) * pow(modulus, -1, new_modulus)) % new_modulus
    return ((residue + modulus * step) % (modulus * new_modulus),
            modulus * new_modulus)


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
    """Return the chosen chart representative modulo the all-ones vector."""
    require(all(bit in (0, 1) for bit in bits), "nonbinary affine coordinate")
    opposite = tuple(1 - bit for bit in bits)
    if len(bits) == 3:
        return tuple(bits) if sum(bits) <= 1 else opposite
    return min(tuple(bits), opposite)


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
    return ancestry_path_for_words(source_word, target_word)


def ancestry_path_in(table, source, target):
    return ancestry_path_for_words(table[source][1], table[target][1])


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


def matching_pairs(direction, vertices=VERTEX_ORDER):
    seen = set()
    pairs = []
    for source in vertices:
        if source in seen:
            continue
        target = quotient_add(source, direction)
        require(target in vertices and target != source,
                "direction failed to act freely on the affine quotient")
        seen.update((source, target))
        pairs.append((source, target))
    require(2 * len(pairs) == len(vertices) and len(seen) == len(vertices),
            "direction did not produce a perfect matching")
    return tuple(pairs)


def audit_flat_action(vertices, directions, path_function):
    involution_checks = 0
    square_checks = 0
    for source in vertices:
        for direction in directions:
            target = quotient_add(source, direction)
            require(reduce_free(path_function(source, target)
                                + path_function(target, source)) == (),
                    "direction involution acquired holonomy")
            involution_checks += 1

        for first in directions:
            for second in directions:
                if first == second:
                    continue
                via_first = quotient_add(source, first)
                target = quotient_add(via_first, second)
                via_second = quotient_add(source, second)
                left_path = reduce_free(path_function(source, via_first)
                                        + path_function(via_first, target))
                right_path = reduce_free(path_function(source, via_second)
                                         + path_function(via_second, target))
                require(left_path == right_path == path_function(source, target),
                        "affine square failed flatness")
                square_checks += 1
    return involution_checks, square_checks


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
             if triple[2] == FIRST_LIMIT}
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

    second_fibre = {triple: word for triple, word in addresses.items()
                    if triple[2] == SECOND_LIMIT}
    require(len(second_fibre) == 8,
            "c=99905 fibre did not have eight vertices")
    for bit, (triple, word) in EXPECTED_99905.items():
        require(euclid_triple(gaussian_parameters(bit)) == triple,
                f"c=99905 Gaussian allocation failed at bit {bit}")
        require(apply_word(word) == triple,
                f"c=99905 word evaluation failed at bit {bit}")
        require(addresses.get(triple) == word,
                f"c=99905 bounded-tree address failed at bit {bit}")
        require(triple[0] ** 2 + triple[1] ** 2 == triple[2] ** 2,
                f"c=99905 Pythagorean identity failed at bit {bit}")
    require(set(second_fibre) == {row[0] for row in EXPECTED_99905.values()},
            "unexpected c=99905 vertex")
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
    direction_labels = tuple(direction for _, direction, _, _ in DIRECTIONS)
    involution_checks, square_checks = audit_flat_action(
        VERTEX_ORDER, direction_labels, ancestry_path
    )

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

    # Scale control: the first eight-parent record fibre at c=99905.
    second_vertices = tuple(sorted(EXPECTED_99905))
    second_costs = {}
    second_jumps = {}
    second_primes = (5, 13, 29, 53)
    for direction in SECOND_DIRECTIONS:
        pairs = matching_pairs(direction, second_vertices)
        costs = []
        jumps = []
        signatures = []
        for left, right in pairs:
            left_word = EXPECTED_99905[left][1]
            right_word = EXPECTED_99905[right][1]
            prefix = common_prefix_length(left_word, right_word)
            path = ancestry_path_in(EXPECTED_99905, left, right)
            distance = len(path)
            jump = abs(len(left_word) - len(right_word))
            require(distance == len(left_word) + len(right_word) - 2 * prefix,
                    "c=99905 tree-distance formula failed")
            require(ancestry_path_in(EXPECTED_99905, right, left) == inverse(path),
                    "c=99905 path reversal failed")
            require(not left_word.startswith(right_word)
                    and not right_word.startswith(left_word),
                    "c=99905 endpoints became ancestor-related")
            costs.append(distance)
            jumps.append(jump)
            signatures.append(path_signature(path))

        product = 1
        for bit, prime in zip(bit_tuple(direction), second_primes):
            if bit:
                product *= prime
        weight = tuple(sorted((product, SECOND_LIMIT // product)))
        require(weight == tuple(sorted(SECOND_WEIGHTS[direction])),
                f"c=99905 folded weight failed at {direction}")
        require(tuple(costs) == SECOND_COSTS[direction]
                and tuple(jumps) == SECOND_JUMPS[direction],
                f"c=99905 path census changed at {direction}")
        require(len(set(costs)) > 1 and len(set(signatures)) > 1,
                f"c=99905 direction {direction} became source-independent")
        require(any(jump != 0 for jump in jumps),
                f"c=99905 direction {direction} escaped depth obstruction")
        second_costs[direction] = tuple(costs)
        second_jumps[direction] = tuple(jumps)

    second_involutions, second_squares = audit_flat_action(
        second_vertices,
        SECOND_DIRECTIONS,
        lambda source, target: ancestry_path_in(EXPECTED_99905, source, target),
    )
    require((second_involutions, second_squares) == (56, 336),
            "c=99905 flat-action census changed")

    # Infinite prime-5 toggle: t=25s+1 has a short explicit partner word.
    infinite_checks = 0
    for s in range(1, INFINITE_FAMILY_BOUND + 1):
        t = 25 * s + 1
        hypotenuse = spine_c(t)
        spine = (2 * t + 1, 2 * t * (t + 1), hypotenuse)
        partner_parameters = (35 * s + 2, 5 * s + 1)
        partner = (
            1200 * s * s + 130 * s + 3,
            350 * s * s + 90 * s + 4,
            1250 * s * s + 150 * s + 5,
        )
        require(matrix_vector(matrix_power(U, 25 * s), ROOT) == spine,
                f"infinite spine word failed at s={s}")
        require(euclid_triple(partner_parameters) == partner
                and partner[2] == hypotenuse,
                f"prime-5 partner formula failed at s={s}")
        require(gcd(*partner_parameters) == 1
                and (partner_parameters[0] - partner_parameters[1]) % 2 == 1,
                f"prime-5 partner lost primitivity at s={s}")

        quotient = (15 * s + 1, 5 * s)
        require(gaussian_multiply((2, 1), quotient) == (t + 1, t),
                f"2+i factorization failed at s={s}")
        require(gaussian_multiply((2, -1), quotient)
                == (partner_parameters[0], -partner_parameters[1]),
                f"prime-5 allocation toggle failed at s={s}")
        require(hypotenuse == 5 * (250 * s * s + 30 * s + 1)
                and (hypotenuse // 5) % 5 == 1,
                f"exact 5-adic valuation failed at s={s}")

        word_value = ROOT
        for matrix in (D, D):
            word_value = matrix_vector(matrix, word_value)
        word_value = matrix_vector(matrix_power(U, s - 1), word_value)
        for matrix in (A, D, D):
            word_value = matrix_vector(matrix, word_value)
        require(word_value == partner,
                f"partner Berggren word failed at s={s}")
        require(abs(25 * s - (s + 4)) == 24 * s - 4,
                f"depth dispersion failed at s={s}")
        require(25 * s + (s + 4) == 26 * s + 4,
                f"path cost failed at s={s}")
        infinite_checks += 1

    # CRT controls show compatibility with arbitrarily high Boolean rank.
    crt_residue = 1
    crt_modulus = 25
    crt_controls = []
    selected_primes = []
    for prime in CRT_PRIMES:
        roots = roots_of_spine_c(prime)
        require(len(roots) == 2,
                f"split prime {prime} did not have two spine roots")
        crt_residue, crt_modulus = crt_pair(
            crt_residue, crt_modulus, roots[0], prime
        )
        selected_primes.append(prime)
        t = crt_residue if crt_residue > 1 else crt_residue + crt_modulus
        require(t % 25 == 1 and all(spine_c(t) % p == 0 for p in selected_primes),
                f"simultaneous rank CRT failed through prime {prime}")
        require((spine_c(t) // 5) % 5 == 1,
                "CRT control lost the exact prime-5 toggle")
        crt_controls.append((len(selected_primes), t, crt_modulus))

    print("THM-3345 VERIFIED-EXACT")
    print(f"berggren_tree_nodes_through_99905={len(addresses)}")
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
    print("c99905_vertices(bit,depth,word,odd_leg):")
    for bit in second_vertices:
        triple, word = EXPECTED_99905[bit]
        shown_word = f"U^{len(word)}" if set(word) == {"U"} else word
        print(f"  {bit}: depth={len(word)}, word={shown_word}, odd_leg={triple[0]}")
    print(f"c99905_direction_weights={SECOND_WEIGHTS}")
    print(f"c99905_direction_costs={second_costs}")
    print(f"c99905_direction_jumps={second_jumps}")
    print(f"c99905_involutions={second_involutions}; affine_square_checks={second_squares}")
    print("c99905_rooted_automorphism_controls=7/7 blocked_by_depth")
    print("c99905_colour_only_path_controls=7/7 blocked_by_source_dispersion")
    print(f"prime5_infinite_family_checks=1..{infinite_checks}; "
          f"max_depth_jump={24 * INFINITE_FAMILY_BOUND - 4}; "
          f"max_path_cost={26 * INFINITE_FAMILY_BOUND + 4}")
    print(f"prime5_rank_CRT_controls={crt_controls}")
    print("ALL CHECKS PASSED")


def ancestry_path_for_words(source_word, target_word):
    return reduce_free(inverse(word_tokens(source_word)) + word_tokens(target_word))


if __name__ == "__main__":
    main()
