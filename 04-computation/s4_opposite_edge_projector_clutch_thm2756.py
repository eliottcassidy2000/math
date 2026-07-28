#!/usr/bin/env python3
"""Exact opposite-edge projector and integral-clutch audit for THM-2756.

The six-edge permutation module of K4 splits over Q into the +1 and -1
eigenspaces of edge complementation.  The positive block is the permutation
module on the three perfect matchings; the negative block is the standard
three-dimensional S4 module.  Both block determinants are the quartic sign,
so their product is the trivial ambient six-edge sign.  Integrally the two
blocks have index eight and are glued by an F2^3 matching-permutation clutch.
All checks use explicit exceptions rather than Python assertions.
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


def induced_on_edges(permutation):
    return tuple(
        EDGES.index(canonical_edge(permutation[a], permutation[b]))
        for a, b in EDGES
    )


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


def permutation_matrix(action):
    size = len(action)
    matrix = [[0] * size for _ in range(size)]
    for source, target in enumerate(action):
        matrix[target][source] = 1
    return tuple(tuple(row) for row in matrix)


def matrix_multiply(left, right):
    return tuple(tuple(
        sum(left[i][k] * right[k][j] for k in range(len(right)))
        for j in range(len(right[0]))
    ) for i in range(len(left)))


def matrix_trace(matrix):
    return sum(matrix[i][i] for i in range(len(matrix)))


def determinant(matrix):
    """Fraction-free Bareiss determinant."""
    work = [list(row) for row in matrix]
    size = len(work)
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        pivot_row = next(
            (row for row in range(pivot_index, size)
             if work[row][pivot_index] != 0),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = (
                work[pivot_row], work[pivot_index]
            )
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for i in range(pivot_index + 1, size):
            for j in range(pivot_index + 1, size):
                numerator = (
                    work[i][j] * pivot
                    - work[i][pivot_index] * work[pivot_index][j]
                )
                require(numerator % previous == 0,
                        "Bareiss exact division failed")
                work[i][j] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def rank_mod_two(matrix):
    work = [[entry % 2 for entry in row] for row in matrix]
    rows = len(work)
    columns = len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for row in range(rows):
            if row != rank and work[row][column]:
                work[row] = [
                    a ^ b for a, b in zip(work[row], work[rank])
                ]
        rank += 1
    return rank


def block_diagonal(first, second):
    n_first = len(first)
    n_second = len(second)
    out = [[0] * (n_first + n_second)
           for _ in range(n_first + n_second)]
    for i in range(n_first):
        for j in range(n_first):
            out[i][j] = first[i][j]
    for i in range(n_second):
        for j in range(n_second):
            out[n_first + i][n_first + j] = second[i][j]
    return tuple(tuple(row) for row in out)


def block_actions(permutation):
    """Return the + matching block and the oriented - block."""
    plus = [[0] * 3 for _ in range(3)]
    minus = [[0] * 3 for _ in range(3)]
    canonical_matchings = tuple(
        frozenset(canonical_edge(*edge) for edge in matching)
        for matching in MATCHINGS
    )
    for source, (first, second) in enumerate(MATCHINGS):
        image_first = canonical_edge(
            permutation[first[0]], permutation[first[1]]
        )
        image_second = canonical_edge(
            permutation[second[0]], permutation[second[1]]
        )
        image_set = frozenset((image_first, image_second))
        target = canonical_matchings.index(image_set)
        target_first, target_second = MATCHINGS[target]
        sign = 1 if (image_first, image_second) == (
            target_first, target_second
        ) else -1
        plus[target][source] = 1
        minus[target][source] = sign
    return (
        tuple(tuple(row) for row in plus),
        tuple(tuple(row) for row in minus),
    )


def eigenbasis_matrix():
    """Columns are p_m=e+e^opp, then n_m=e-e^opp."""
    columns = []
    for sign in (1, -1):
        for first, second in MATCHINGS:
            column = [0] * 6
            column[EDGES.index(first)] = 1
            column[EDGES.index(second)] = sign
            columns.append(column)
    return tuple(tuple(columns[j][i] for j in range(6)) for i in range(6))


def clutch(vector):
    return tuple(
        (vector[EDGES.index(first)] + vector[EDGES.index(second)]) % 2
        for first, second in MATCHINGS
    )


def act_vector(permutation, vector):
    action = induced_on_edges(permutation)
    out = [0] * 6
    for source, target in enumerate(action):
        out[target] = vector[source]
    return tuple(out)


def act_matching_mod_two(plus, vector):
    return tuple(
        sum(plus[i][j] * vector[j] for j in range(3)) % 2
        for i in range(3)
    )


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
    basis = eigenbasis_matrix()
    require(abs(determinant(basis)) == 8,
            "the integral plus/minus lattice index stopped being eight")
    require(rank_mod_two(basis) == 3,
            "the integral clutch stopped having three binary coordinates")

    opposite = {}
    for first, second in MATCHINGS:
        opposite[EDGES.index(first)] = EDGES.index(second)
        opposite[EDGES.index(second)] = EDGES.index(first)
    require(tuple(opposite[opposite[i]] for i in range(6))
            == tuple(range(6)), "edge complementation stopped being involutive")

    rows = {}
    matching_kernel = []
    plus_actions = []
    minus_actions = []
    for g in s4:
        edge_action = induced_on_edges(g)
        edge_matrix = permutation_matrix(edge_action)
        plus, minus = block_actions(g)
        plus_actions.append(plus)
        minus_actions.append(minus)

        require(all(
            edge_action[opposite[index]] == opposite[edge_action[index]]
            for index in range(6)
        ), "opposite-edge involution stopped commuting with S4")
        require(matrix_multiply(edge_matrix, basis)
                == matrix_multiply(basis, block_diagonal(plus, minus)),
                "the rational opposite-edge block diagonalization failed")
        require(determinant(plus) == permutation_sign(g)
                and determinant(minus) == permutation_sign(g),
                "one block stopped carrying the quartic sign")
        require(determinant(edge_matrix)
                == determinant(plus) * determinant(minus) == 1,
                "the two determinant characters stopped cancelling")
        fixed_vertices = sum(g[i] == i for i in range(4))
        require(matrix_trace(minus) == fixed_vertices - 1,
                "the negative block stopped being the standard module")

        # Equivariance of the integral parity clutch on the six basis vectors.
        for edge_index in range(6):
            vector = tuple(int(i == edge_index) for i in range(6))
            require(
                clutch(act_vector(g, vector))
                == act_matching_mod_two(plus, clutch(vector)),
                "the F2^3 clutch stopped transforming as the matching module",
            )

        if plus == ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            matching_kernel.append(g)

        key = cycle_type(g)
        value = (
            matrix_trace(plus), determinant(plus),
            matrix_trace(minus), determinant(minus),
            matrix_trace(edge_matrix), determinant(edge_matrix),
        )
        if key in rows:
            require(rows[key][0] == value,
                    "a conjugacy class split in the block-character table")
            rows[key] = (value, rows[key][1] + 1)
        else:
            rows[key] = (value, 1)

    expected_rows = {
        (1, 1, 1, 1): ((3, 1, 3, 1, 6, 1), 1),
        (2, 1, 1): ((1, -1, 1, -1, 2, 1), 6),
        (2, 2): ((3, 1, -1, 1, 2, 1), 3),
        (3, 1): ((0, 1, 0, 1, 0, 1), 8),
        (4,): ((1, -1, -1, -1, 0, 1), 6),
    }
    require(rows == expected_rows, "the block-character table changed")
    require(len(set(plus_actions)) == 6 and len(set(minus_actions)) == 24,
            "the quotient/faithful block image sizes changed")

    # The matching block is 1+[22]: subtracting its invariant line gives the
    # irreducible degree-two character (2,0,2,-1,0), of norm one.
    class_order = ((1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,))
    matching_zero_sum_character = tuple(rows[key][0][0] - 1
                                        for key in class_order)
    require(matching_zero_sum_character == (2, 0, 2, -1, 0),
            "the matching zero-sum block stopped being [22]")
    class_sizes = tuple(rows[key][1] for key in class_order)
    character_norm_numerator = sum(
        size * character * character
        for size, character in zip(class_sizes, matching_zero_sum_character)
    )
    require(character_norm_numerator == 24,
            "the [22] character stopped having norm one")

    require(len(matching_kernel) == 4
            and {cycle_type(g) for g in matching_kernel}
            == {(1, 1, 1, 1), (2, 2)},
            "the positive block kernel stopped being V4")

    # Direct congruence model for the integral sum L_+ + L_-.
    for residues in range(4 ** 6):
        value = residues
        vector = []
        for _ in range(6):
            vector.append(value % 4)
            value //= 4
        same_parity_on_pairs = all(
            (vector[EDGES.index(first)] - vector[EDGES.index(second)]) % 2 == 0
            for first, second in MATCHINGS
        )
        require((clutch(vector) == (0, 0, 0)) == same_parity_on_pairs,
                "the integral clutch kernel congruence changed")

    print("S4 OPPOSITE-EDGE PROJECTOR / INTEGRAL CLUTCH AUDIT")
    print("class size tr_plus det_plus tr_minus det_minus tr_edge det_edge")
    for source_type in class_order:
        value, size = rows[source_type]
        print(
            f"{format_type(source_type):7s} {size:2d} "
            f"{value[0]:+d} {value[1]:+d} "
            f"{value[2]:+d} {value[3]:+d} "
            f"{value[4]:+d} {value[5]:+d}"
        )
    print("opposition=central_involution rational_split=E_plus+E_minus")
    print("E_plus=matching_permutation=1+[22] kernel=V4 quotient=S3")
    print("E_minus=standard_[31] faithful_character=fixed_vertices-1")
    print("determinants=det_plus=det_minus=quartic_sign product=+1")
    print("trace_collision=transposition:(1+1) double:(3-1) total=2")
    print("integral_basis_index=8 mod2_rank=3 smith=1^3,2^3")
    print("integral_clutch=L/(L_plus+L_minus)=F2^3 matching_permutation_module")
    print("SCOPE: exact finite representation/lattice theorem; no Keller or LRC exclusion")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
