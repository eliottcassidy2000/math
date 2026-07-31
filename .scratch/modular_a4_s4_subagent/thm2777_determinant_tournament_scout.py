#!/usr/bin/env python3
"""Exact scout for the proposed marked D3 determinant tournament.

The scout separates the genuinely preserved weighted chirotope from the
chamber-dependent tournament representative and checks its blindness to the
110 affine Kummer boundary word.
"""

from itertools import combinations, permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def determinant(columns):
    a, b, c = columns
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def matvec(matrix, vector):
    return tuple(sum(matrix[i][j] * vector[j] for j in range(3))
                 for i in range(3))


def determinant_three(matrix):
    return determinant(tuple(tuple(matrix[row][column] for row in range(3))
                             for column in range(3)))


def all_weyl_d3():
    out = set()
    for coordinate_permutation in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            if signs[0] * signs[1] * signs[2] != 1:
                continue
            matrix = []
            for row_index in range(3):
                row = [0, 0, 0]
                row[coordinate_permutation[row_index]] = signs[row_index]
                matrix.append(tuple(row))
            out.add(tuple(matrix))
    return frozenset(out)


def tournament(roots, h):
    adjacency = [[False] * len(roots) for _ in roots]
    weights = {}
    for first, second in combinations(range(len(roots)), 2):
        value = determinant((roots[first], roots[second], h))
        require(value != 0, "determinant orientation acquired a tie")
        adjacency[first][second] = value > 0
        adjacency[second][first] = value < 0
        weights[(first, second)] = abs(value)
    return adjacency, weights


def score_sequence(adjacency):
    return tuple(sorted(sum(row) for row in adjacency))


def cyclic_triangle_count(adjacency):
    count = 0
    for first, second, third in combinations(range(len(adjacency)), 3):
        if ((adjacency[first][second]
             and adjacency[second][third]
             and adjacency[third][first])
                or (adjacency[second][first]
                    and adjacency[first][third]
                    and adjacency[third][second])):
            count += 1
    return count


def strongly_connected(adjacency):
    for source in range(len(adjacency)):
        reached = {source}
        frontier = [source]
        while frontier:
            vertex = frontier.pop()
            for target in range(len(adjacency)):
                if adjacency[vertex][target] and target not in reached:
                    reached.add(target)
                    frontier.append(target)
        if len(reached) != len(adjacency):
            return False
    return True


def hamiltonian_path_count(adjacency):
    return sum(
        all(adjacency[path[index]][path[index + 1]]
            for index in range(len(path) - 1))
        for path in permutations(range(len(adjacency)))
    )


def switched(adjacency, signs):
    return tuple(tuple(
        False if first == second
        else adjacency[first][second] ^ (signs[first] != signs[second])
        for second in range(len(adjacency))
    ) for first in range(len(adjacency)))


def relabelled(adjacency, permutation):
    return tuple(tuple(
        adjacency[permutation[first]][permutation[second]]
        for second in range(len(adjacency))
    ) for first in range(len(adjacency)))


def reversed_tournament(adjacency):
    return tuple(tuple(
        False if first == second else not adjacency[first][second]
        for second in range(len(adjacency))
    ) for first in range(len(adjacency)))


def weight_at(weights, first, second):
    return weights[tuple(sorted((first, second)))]


def line_index(vector, roots):
    if vector in roots:
        return roots.index(vector)
    negative = tuple(-entry for entry in vector)
    require(negative in roots, "W(D3) image left the six root lines")
    return roots.index(negative)


def line_action(matrix, roots):
    return tuple(line_index(matvec(matrix, root), roots) for root in roots)


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


def main():
    h = (1, 1, 1)
    names = ("a12", "a23", "a13", "s12", "s13", "s23")
    roots = (
        (1, -1, 0),
        (0, 1, -1),
        (1, 0, -1),
        (1, 1, 0),
        (1, 0, 1),
        (0, 1, 1),
    )
    adjacency, weights = tournament(roots, h)
    weight_counts = {
        weight: sum(value == weight for value in weights.values())
        for weight in sorted(set(weights.values()))
    }
    require(weight_counts == {1: 9, 2: 3, 3: 3},
            "weighted determinant spectrum changed")

    a2_vertices = {0, 1, 2}
    opposition = {(0, 3), (1, 5), (2, 4)}
    weight_three = {
        pair for pair, weight in weights.items() if weight == 3
    }
    weight_two = {
        pair for pair, weight in weights.items() if weight == 2
    }
    require(weight_three == set(combinations(sorted(a2_vertices), 2)),
            "weight-three stratum stopped being the A2 triangle")
    require(weight_two == opposition,
            "weight-two stratum stopped being the opposition matching")
    require(score_sequence(adjacency) == (1, 2, 2, 3, 3, 4),
            "marked tournament score sequence changed")
    require(cyclic_triangle_count(adjacency) == 6,
            "marked tournament cyclic-triangle count changed")
    require(strongly_connected(adjacency),
            "marked determinant tournament stopped being strongly connected")
    require(hamiltonian_path_count(adjacency) == 23,
            "marked tournament Hamiltonian-path count changed")
    adjacency_tuple = tuple(tuple(row) for row in adjacency)
    reverse = reversed_tournament(adjacency)
    automorphisms = [
        permutation for permutation in permutations(range(6))
        if relabelled(adjacency, permutation) == adjacency_tuple
    ]
    anti_isomorphisms = [
        permutation for permutation in permutations(range(6))
        if relabelled(adjacency, permutation) == reverse
    ]
    require(automorphisms == [tuple(range(6))],
            "marked tournament stopped being rigid")
    require(not anti_isomorphisms,
            "marked tournament unexpectedly became self-converse")

    # Six A2 chambers are the six coordinate orders.  They only choose signs
    # on the three A2 root lines.  The resulting tournaments are vertex
    # switchings of the baseline, so the chamber-free object is a switching
    # class rather than one canonical tournament.
    difference_lines = {
        frozenset((0, 1)): 0,
        frozenset((1, 2)): 1,
        frozenset((0, 2)): 2,
    }
    chamber_tournaments = set()
    for order in permutations(range(3)):
        chosen = [None, None, None]
        for early, late in combinations(range(3), 2):
            first = order[early]
            second = order[late]
            vector = tuple(int(index == first) - int(index == second)
                           for index in range(3))
            slot = difference_lines[frozenset((first, second))]
            chosen[slot] = vector
        chamber_roots = tuple(chosen) + roots[3:]
        chamber_adjacency, _ = tournament(chamber_roots, h)
        sign_switch = tuple(
            chamber_roots[index] != roots[index] if index < 3 else False
            for index in range(6)
        )
        require(tuple(tuple(row) for row in chamber_adjacency)
                == switched(adjacency, sign_switch),
                "A2 chamber change stopped being vertex switching")
        chamber_tournaments.add(tuple(tuple(row) for row in chamber_adjacency))
    require(len(chamber_tournaments) == 6,
            "six A2 chambers stopped giving six marked representatives")

    # The marked tournament is chiral, but its reverse returns after one
    # allowed chamber switch and a weight-preserving coordinate relabelling.
    # Quartic sign is therefore lost in the chamber-free switching quotient.
    reverse_switch_mask = (True, False, False, False, False, False)
    reverse_relabelling = (0, 2, 1, 3, 5, 4)
    require(relabelled(switched(adjacency, reverse_switch_mask),
                       reverse_relabelling) == reverse,
            "reverse stopped lying in the chamber switching class")
    require(all(
        weight_at(weights, first, second)
        == weight_at(weights,
                     reverse_relabelling[first],
                     reverse_relabelling[second])
        for first, second in combinations(range(6), 2)
    ), "reverse switching witness stopped preserving determinant weights")

    # Covariance is the actual invariant content:
    # det(g alpha,g beta,g h)=det(g)det(alpha,beta,h).
    weyl = all_weyl_d3()
    require(len(weyl) == 24, "W(D3) census changed")
    for matrix in weyl:
        matrix_sign = determinant_three(matrix)
        require(matrix_sign in (-1, 1), "W(D3) determinant left +/-1")
        moved_h = matvec(matrix, h)
        for first, second in combinations(range(6), 2):
            before = determinant((roots[first], roots[second], h))
            after = determinant((
                matvec(matrix, roots[first]),
                matvec(matrix, roots[second]),
                moved_h,
            ))
            require(after == matrix_sign * before,
                    "determinant-tournament covariance failed")

    x_s = (
        (1, 0, 0),
        (0, 0, -1),
        (0, -1, 0),
    )
    y_s = (
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 0),
    )
    face_square = (
        (-1, 0, 0),
        (0, -1, 0),
        (0, 0, 1),
    )
    require(determinant_three(x_s) == -1
            and determinant_three(y_s) == 1
            and determinant_three(face_square) == 1,
            "binary/ternary/face-square determinant characters changed")
    require(cycle_type(line_action(x_s, roots)) == (2, 2, 1, 1)
            and cycle_type(line_action(face_square, roots))
            == (2, 2, 1, 1),
            "binary reflection and 110 face square stopped colliding on lines")
    require(cycle_type(line_action(y_s, roots)) == (3, 3),
            "ternary line spectrum changed")

    # The tournament inputs contain no divisor valuation row.  In particular
    # both 000 and 110 are even, while every V4 element has determinant +1.
    boundary_words = {(0, 0, 0), (1, 1, 0)}
    require(all(sum(word) % 2 == 0 for word in boundary_words)
            and determinant_three(identity_face := (
                (1, 0, 0), (0, 1, 0), (0, 0, 1)
            )) == determinant_three(face_square) == 1,
            "boundary-blindness control changed")
    require(identity_face != face_square,
            "110 face square unexpectedly became the identity")

    print("THM-2777 MARKED D3 DETERMINANT-TOURNAMENT SCOUT")
    print("tie_count=0 weight_spectrum={1:9,2:3,3:3}")
    print("weight3=A2_triangle weight2=three_opposition_pairs")
    print("tournament_scores=(1,2,2,3,3,4)")
    print("cyclic_triangles=6 strongly_connected=True Hamiltonian_paths=23")
    print("marked_tournament=rigid_chiral")
    print("A2_chambers=6 distinct_tournaments=6 one_switching_class")
    print("W(D3)_covariance=even_preserves_odd_globally_reverses")
    print("marked_preserved_character=quartic_sign_via_global_reversal")
    print("chamber_free_weighted_switching_class=self_converse")
    print("X_line_cycle=(2,2,1,1) face110_line_cycle=(2,2,1,1)")
    print("Y_line_cycle=(3,3)")
    print("orientation_character_blindness=000_and_110_both_determinant_plus")
    print("labelled_line_action_distinguishes_identity_from_face110")
    print("VERDICT=NONCOSMETIC_SIGN_ONLY_WITH_MARKED_CHAMBER")
    print("CHAMBER_FREE_SURVIVOR=WEIGHTED_A2_OPPOSITION_INCIDENCE")
    print("NO_NEW_AFFINE_KUMMER_OR_LRC_COORDINATE")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
