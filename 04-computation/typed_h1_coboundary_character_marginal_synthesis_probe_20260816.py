#!/usr/bin/env python3
"""Independent finite probe for the typed H^1/coboundary synthesis.

This file deliberately imports none of the four source companions.  It checks
only finite consequences used by the accompanying reflection:

* incidence ranks and Betti numbers for K4, C7, and the 13-edge role graph;
* functorial exactness of the primitive-view and role-potential gradients;
* the tetrahedral opposite-face transform, including its characteristic-two
  collapse and the audited degree-three THM-3494 four-view index packet;
* the cyclic-cover seam/deck-defect construction at p=13 and p=2;
* graph H^1 versus group characters for the affine Berggren D4;
* the variable- and fixed-drift automaton character gates; and
* the 13 by 13 marginal tensor kernel and its checkerboard hostile.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from json import dumps


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rank_mod(matrix, prime):
    a = [[entry % prime for entry in row] for row in matrix]
    if not a:
        return 0
    rows = len(a)
    cols = len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [(inv * x) % prime for x in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col]:
                scalar = a[r][col]
                a[r] = [
                    (x - scalar * y) % prime
                    for x, y in zip(a[r], a[rank])
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def incidence(vertex_count, edges):
    matrix = []
    for tail, head in edges:
        row = [0] * vertex_count
        row[tail] = -1
        row[head] = 1
        matrix.append(row)
    return matrix


def matrix_vector_mod(matrix, vector, prime):
    return tuple(
        sum(entry * value for entry, value in zip(row, vector)) % prime
        for row in matrix
    )


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def determinant_integer(matrix):
    """Fraction-free Bareiss determinant for a square integer matrix."""
    a = [row[:] for row in matrix]
    size = len(a)
    require(all(len(row) == size for row in a), "determinant matrix is not square")
    if size == 0:
        return 1
    sign = 1
    previous_pivot = 1
    for column in range(size - 1):
        pivot = next((row for row in range(column, size) if a[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            sign = -sign
        pivot_value = a[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = (
                    a[row][target] * pivot_value
                    - a[row][column] * a[column][target]
                )
                require(numerator % previous_pivot == 0, "Bareiss nonintegral division")
                a[row][target] = numerator // previous_pivot
            a[row][column] = 0
        previous_pivot = pivot_value
    return sign * a[-1][-1]


def triangle_boundary(vertices, edge_index, edge_count):
    a, b, c = vertices
    vector = [0] * edge_count
    for tail, head, coefficient in ((b, c, 1), (a, c, -1), (a, b, 1)):
        edge = (min(tail, head), max(tail, head))
        sign = 1 if tail < head else -1
        vector[edge_index[edge]] += coefficient * sign
    return vector


def opposite_face_matrix(ordered_vertices, all_edges, total_vertices):
    """Vertex values -> boundaries of opposite faces of an oriented K4."""
    edge_index = {edge: i for i, edge in enumerate(all_edges)}
    matrix = [[0] * total_vertices for _ in all_edges]
    for omitted, vertex in enumerate(ordered_vertices):
        face = tuple(v for i, v in enumerate(ordered_vertices) if i != omitted)
        boundary = triangle_boundary(face, edge_index, len(all_edges))
        face_sign = -1 if omitted % 2 else 1
        for edge, coefficient in enumerate(boundary):
            matrix[edge][vertex] += face_sign * coefficient
    return matrix


def edge_action_matrix(vertex_permutation, edges):
    edge_index = {edge: i for i, edge in enumerate(edges)}
    matrix = [[0] * len(edges) for _ in edges]
    for source, (tail, head) in enumerate(edges):
        moved_tail = vertex_permutation[tail]
        moved_head = vertex_permutation[head]
        target_edge = (min(moved_tail, moved_head), max(moved_tail, moved_head))
        sign = 1 if moved_tail < moved_head else -1
        matrix[edge_index[target_edge]][source] = sign
    return matrix


def graph_profile(vertex_count, edges):
    d0 = incidence(vertex_count, edges)
    rank = rank_mod(d0, 1000003)
    require(rank == vertex_count - 1, "graph must be connected")
    return {
        "vertices": vertex_count,
        "edges": len(edges),
        "B1": rank,
        "beta1": len(edges) - rank,
    }


def perm_compose(left, right):
    """Return left after right."""
    return tuple(left[right[i]] for i in range(len(right)))


def perm_inverse(p):
    out = [0] * len(p)
    for i, value in enumerate(p):
        out[value] = i
    return tuple(out)


def perm_parity(p):
    inversions = sum(
        p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p))
    )
    return inversions % 2


def cycle_type(p):
    seen = set()
    lengths = []
    for start in range(len(p)):
        if start in seen:
            continue
        cur = start
        length = 0
        while cur not in seen:
            seen.add(cur)
            cur = p[cur]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        old = frontier.pop()
        for generator in generators:
            new = perm_compose(generator, old)
            if new not in group:
                group.add(new)
                frontier.append(new)
    return tuple(sorted(group))


def all_f2_characters(group):
    identity = tuple(range(len(group[0])))
    index = {g: i for i, g in enumerate(group)}
    multiplication = [
        [index[perm_compose(a, b)] for b in group]
        for a in group
    ]
    chars = []
    for values in product((0, 1), repeat=len(group)):
        if values[index[identity]]:
            continue
        if all(
            values[multiplication[i][j]] == (values[i] ^ values[j])
            for i in range(len(group))
            for j in range(len(group))
        ):
            chars.append({g: values[index[g]] for g in group})
    return chars


def matching_action(p):
    matchings = (
        frozenset((frozenset((0, 1)), frozenset((2, 3)))),
        frozenset((frozenset((0, 2)), frozenset((1, 3)))),
        frozenset((frozenset((0, 3)), frozenset((1, 2)))),
    )
    lookup = {matching: i for i, matching in enumerate(matchings)}
    image = []
    for matching in matchings:
        moved = frozenset(
            frozenset((p[min(edge)], p[max(edge)])) for edge in matching
        )
        image.append(lookup[moved])
    return tuple(image)


def affine_perm(linear_bit, translation):
    # V4={0,p,q,r} is encoded as two bits, with p=1, q=2, r=3.
    def swap_bits(x):
        return ((x & 1) << 1) | ((x & 2) >> 1)

    return tuple(
        ((swap_bits(x) if linear_bit else x) ^ translation)
        for x in range(4)
    )


def cyclic_transgression(base_length, prime, coefficient_dim):
    """Check every edge/coefficient basis class on the prime-degree cover."""
    checks = 0
    for edge in range(base_length):
        for coordinate in range(coefficient_dim):
            base = [[0] * coefficient_dim for _ in range(base_length)]
            base[edge][coordinate] = 1
            seam = [sum(row[j] for row in base) % prime for j in range(coefficient_dim)]
            pulled = [base[i % base_length] for i in range(base_length * prime)]
            primitive = [[0] * coefficient_dim]
            for value in pulled:
                primitive.append(
                    [
                        (primitive[-1][j] + value[j]) % prime
                        for j in range(coefficient_dim)
                    ]
                )
            require(primitive[-1] == [0] * coefficient_dim, "cover did not exactify")
            cover_vertices = primitive[:-1]
            for i in range(base_length * prime):
                shifted = (i + base_length) % (base_length * prime)
                defect = [
                    (cover_vertices[shifted][j] - cover_vertices[i][j]) % prime
                    for j in range(coefficient_dim)
                ]
                require(defect == seam, "deck defect is not the seam")
            checks += 1
    return checks


def marginal_map_rank(size, prime):
    # Rows followed by columns.  Their totals satisfy one relation.
    matrix = [[0] * (size * size) for _ in range(2 * size)]
    for i in range(size):
        for j in range(size):
            matrix[i][i * size + j] = 1
            matrix[size + j][i * size + j] = 1
    return rank_mod(matrix, prime)


def main():
    k4_edges = tuple(combinations(range(4), 2))
    c6_edges = tuple((i, (i + 1) % 6) for i in range(6))
    c7_edges = tuple((i, (i + 1) % 7) for i in range(7))
    role_edges = tuple(
        sorted(
            set(combinations((0, 3, 4, 5), 2))
            | set(combinations((1, 2, 4, 7), 2))
            | {(4, 6)}
        )
    )
    k4 = graph_profile(4, k4_edges)
    c6 = graph_profile(6, c6_edges)
    c7 = graph_profile(7, c7_edges)
    role = graph_profile(8, role_edges)
    require(k4 == {"vertices": 4, "edges": 6, "B1": 3, "beta1": 3}, "K4 profile")
    require(c6 == {"vertices": 6, "edges": 6, "B1": 5, "beta1": 1}, "C6 profile")
    require(c7 == {"vertices": 7, "edges": 7, "B1": 6, "beta1": 1}, "C7 profile")
    require(role == {"vertices": 8, "edges": 13, "B1": 7, "beta1": 6}, "role profile")

    # Oriented opposite-face map on K4.  It is not the graph gradient.
    omega_k4 = opposite_face_matrix((0, 1, 2, 3), k4_edges, 4)
    k4_boundary = transpose(incidence(4, k4_edges))
    for column in transpose(omega_k4):
        require(matrix_vector_mod(k4_boundary, column, 13) == (0, 0, 0, 0), "face boundary")
    require(rank_mod(omega_k4, 13) == 3, "opposite-face rank")
    require(
        rank_mod([cut + curl for cut, curl in zip(incidence(4, k4_edges), omega_k4)], 13) == 6,
        "K4 cut/cycle direct sum",
    )
    require(
        matrix_vector_mod(omega_k4, (1, 1, 1, 1), 13) == (0,) * 6,
        "constant potential must vanish",
    )
    k4_split_minor = [
        incidence_row[:3] + omega_row[:3]
        for incidence_row, omega_row in zip(incidence(4, k4_edges), omega_k4)
    ]
    require(determinant_integer(k4_split_minor) == -16, "K4 split lattice index")

    # In characteristic two all four oriented face boundaries have one common
    # nonzero cohomology class.  Thus Omega remembers exactly the sum of the
    # four vertex square classes, i.e. the square class of their product.
    k4_incidence = incidence(4, k4_edges)
    require(rank_mod(k4_incidence, 2) == 3, "K4 cut rank in characteristic two")
    require(rank_mod(omega_k4, 2) == 3, "K4 face rank in characteristic two")
    require(
        rank_mod([cut + curl for cut, curl in zip(k4_incidence, omega_k4)], 2) == 4,
        "K4 characteristic-two quotient rank",
    )
    first_face = tuple(row[0] for row in omega_k4)
    require(
        rank_mod([row + [first_face[index]] for index, row in enumerate(k4_incidence)], 2)
        == 4,
        "K4 face class should be nonzero in characteristic two",
    )
    for vertex in range(1, 4):
        face_difference = tuple(
            row[vertex] - row[0]
            for row in omega_k4
        )
        require(
            rank_mod(
                [row + [face_difference[index]] for index, row in enumerate(k4_incidence)],
                2,
            )
            == 3,
            "K4 face classes should agree in characteristic two",
        )

    # THM-3494's audited n=3 primitive views (w,x,y,z) have index roots
    # 1, C^3*f, C^3*(27Q), C^6*I_z.  At the prime divisor
    # f=9P-27Q-2, I_z does not vanish: setting Q=0 and P=2/9 gives
    # 32/2187.  Hence f has odd valuation in their product and the resulting
    # characteristic-two opposite-face H^1 class is genuinely nonzero.
    p_on_index_divisor = Fraction(2, 9)
    iz_on_index_divisor = (
        972 * p_on_index_divisor**7
        + 756 * p_on_index_divisor**6
        - 189 * p_on_index_divisor**5
    )
    require(iz_on_index_divisor == Fraction(32, 2187), "THM-3494 z-index divisor test")
    require(iz_on_index_divisor != 0, "THM-3494 four-view product became square")
    for vertex_permutation in permutations(range(4)):
        edge_action = edge_action_matrix(vertex_permutation, k4_edges)
        orientation_sign = -1 if perm_parity(vertex_permutation) else 1
        for vertex in range(4):
            old_column = tuple(row[vertex] for row in omega_k4)
            moved_column = matrix_vector_mod(edge_action, old_column, 13)
            target_column = tuple(
                orientation_sign * row[vertex_permutation[vertex]] % 13
                for row in omega_k4
            )
            require(moved_column == target_column, "tetrahedral sign covariance")

    # The two-K4-plus-bridge role graph has a rank-six oriented cycle map.
    wing_left = (0, 3, 4, 5)
    wing_right = (1, 2, 4, 7)
    omega_left = opposite_face_matrix(wing_left, role_edges, 8)
    omega_right = opposite_face_matrix(wing_right, role_edges, 8)
    omega_role = [
        [left + right for left, right in zip(left_row, right_row)]
        for left_row, right_row in zip(omega_left, omega_right)
    ]
    role_incidence = incidence(8, role_edges)
    role_boundary = transpose(role_incidence)
    for column in transpose(omega_role):
        require(matrix_vector_mod(role_boundary, column, 13) == (0,) * 8, "role cycle")
    require(rank_mod(omega_role, 13) == 6, "role opposite-face rank")
    require(
        rank_mod([cut + curl for cut, curl in zip(role_incidence, omega_role)], 13) == 13,
        "role cut/cycle direct sum",
    )
    bridge_row = [0] * 8
    bridge_row[4] = -1
    bridge_row[6] = 1
    require(rank_mod(omega_role + [bridge_row], 13) == 7, "H1 plus bridge decomposition")
    require(
        matrix_vector_mod(omega_role, (0, 0, 0, 0, 0, 0, 1, 0), 13) == (0,) * 13,
        "leaf direction must be the sharp H1 kernel",
    )
    role_split_minor = [
        incidence_row[:7] + omega_row[:6]
        for incidence_row, omega_row in zip(role_incidence, omega_role)
    ]
    require(determinant_integer(role_split_minor) == -256, "role split lattice index")

    role_edge_set = set(role_edges)
    role_automorphisms = tuple(
        vertex_permutation
        for vertex_permutation in permutations(range(8))
        if {
            (min(vertex_permutation[a], vertex_permutation[b]),
             max(vertex_permutation[a], vertex_permutation[b]))
            for a, b in role_edges
        } == role_edge_set
    )
    require(len(role_automorphisms) == 72, "role automorphism count")
    invariant_constraints = [row[:] for row in role_boundary]
    identity_edges = tuple(range(len(role_edges)))
    for vertex_permutation in role_automorphisms:
        action = edge_action_matrix(vertex_permutation, role_edges)
        for row in range(len(role_edges)):
            invariant_constraints.append(
                [
                    action[row][column] - (1 if identity_edges[row] == column else 0)
                    for column in range(len(role_edges))
                ]
            )
    invariant_cycle_dimension = len(role_edges) - rank_mod(invariant_constraints, 13)
    require(invariant_cycle_dimension == 0, "role H1 should have no invariant scalar line")

    certified_prime = 572252886246508880869
    role_values = {
        "c1": 405336876493642499425,
        "c2": 503604956476841920373,
        "c3": 518539850465495448196,
        "H": 320618948602619577408,
        "q2": 15703541686881447885,
        "q3": 503604956476841920373,
        "q4": 503604956476841920373,
        "q5": 503604956476841920373,
    }
    require(
        rank_mod([cut + curl for cut, curl in zip(role_incidence, omega_role)], certified_prime)
        == 13,
        "role cut/cycle split at the certified prime",
    )
    left_outer = tuple(v for v in wing_left if v != 4)
    right_outer = tuple(v for v in wing_right if v != 4)
    role_class_rows = []
    both_wings_nonzero = 0
    bridge_nonzero = 0
    for swap_wings in (False, True):
        c_vertices, q_vertices = (
            (right_outer, left_outer) if swap_wings else (left_outer, right_outer)
        )
        for c_order in permutations(c_vertices):
            for q_order in permutations(q_vertices):
                potential = [0] * 8
                potential[4] = role_values["H"]
                potential[6] = role_values["q5"]
                for role_name, vertex in zip(("c1", "c2", "c3"), c_order):
                    potential[vertex] = role_values[role_name]
                for role_name, vertex in zip(("q2", "q3", "q4"), q_order):
                    potential[vertex] = role_values[role_name]
                left_class = matrix_vector_mod(omega_left, potential, certified_prime)
                right_class = matrix_vector_mod(omega_right, potential, certified_prime)
                full_class = matrix_vector_mod(omega_role, potential, certified_prime)
                if any(left_class) and any(right_class):
                    both_wings_nonzero += 1
                bridge = (potential[6] - potential[4]) % certified_prime
                if bridge:
                    bridge_nonzero += 1
                role_class_rows.append((full_class, bridge))
    require(len(role_class_rows) == 72, "role chart count")
    require(both_wings_nonzero == 72, "both H1 wing classes must fire")
    require(bridge_nonzero == 72, "bridge sidecar must fire")
    role_class_digest = sha256(
        dumps(role_class_rows, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    # The universal primitive-view cochain has edge vector e_head-e_tail.
    # Every pullback around C7 telescopes, for every one of the 4^7 vertex maps.
    pullback_maps = 0
    for vertex_map in product(range(4), repeat=7):
        total = [0] * 4
        for i in range(7):
            total[vertex_map[(i + 1) % 7]] += 1
            total[vertex_map[i]] -= 1
        require(total == [0] * 4, "exact K4 cochain acquired a C7 seam")
        pullback_maps += 1
    require(pullback_maps == 4**7, "pullback census")

    transgression_13 = cyclic_transgression(7, 13, 1)
    transgression_2 = cyclic_transgression(6, 2, 2)
    require(transgression_13 == 7, "p=13 transgression census")
    require(transgression_2 == 12, "p=2 transgression census")

    # Affine Berggren D4.
    G = affine_perm(1, 1)  # Ju+p, a four-cycle
    T = affine_perm(0, 1)  # u+p, a double transposition
    d4 = generated_group((G, T))
    require(len(d4) == 8, "affine group is not D4 of order eight")
    chars = all_f2_characters(d4)
    require(len(chars) == 4, "H^1_grp(D4;F2) should have four elements")
    translations = {t: affine_perm(0, t) for t in range(4)}
    restriction_rows = sorted(
        tuple(character[translations[t]] for t in (1, 2, 3))
        for character in chars
    )
    require(
        restriction_rows == [(0, 0, 0), (0, 0, 0), (1, 1, 0), (1, 1, 0)],
        "D4 character restrictions to V4",
    )
    common_translation_kernel = tuple(
        t
        for t in range(4)
        if all(character[translations[t]] == 0 for character in chars)
    )
    require(common_translation_kernel == (0, 3), "commutator translation should be r")

    eta = next(
        character for character in chars if character[G] == 1 and character[T] == 1
    )
    sign_d4 = {g: perm_parity(g) for g in d4}
    require(sign_d4[G] == 1 and sign_d4[T] == 0, "D4 sign character")
    require(eta[translations[3]] == 0, "period character must miss r=G^2")

    # Variable-translation language: independently test all 48 (S4, linear-bit) pairs.
    s4 = tuple(permutations(range(4)))
    variable_formula_checks = 0
    variable_states = set()
    variable_accept_states = set()
    variable_kernel_states = set()
    for source in s4:
        mu = matching_action(source)
        for epsilon in (0, 1):
            passes = any(
                cycle_type(source) == cycle_type(affine_perm(epsilon, t))
                for t in range(4)
            )
            formula = (
                (epsilon == 0 and mu == (0, 1, 2))
                or (epsilon == 1 and cycle_type(mu) == (1, 2))
            )
            require(passes == formula, "variable-language quotient formula")
            state = (mu, epsilon)
            variable_states.add(state)
            if passes:
                variable_accept_states.add(state)
            if perm_parity(mu) == epsilon:
                variable_kernel_states.add(state)
            variable_formula_checks += 1
    require(len(variable_states) == 12, "variable state count")
    require(len(variable_accept_states) == 4, "variable accept count")
    require(len(variable_kernel_states) == 6, "variable character-kernel count")
    require(variable_accept_states <= variable_kernel_states, "variable sign gate")

    r3 = (1, 2, 0)
    s3 = (1, 0, 2)
    variable_generators = {"A": (r3, 1), "B": (s3, 1), "C": (r3, 0)}
    a_parity_values = {
        letter: perm_parity(mu) ^ epsilon
        for letter, (mu, epsilon) in variable_generators.items()
    }
    require(a_parity_values == {"A": 1, "B": 0, "C": 0}, "variable A parity")
    require((r3, 0) in variable_kernel_states, "C hostile should pass character")
    require((r3, 0) not in variable_accept_states, "C hostile should fail language")

    # Fixed-drift language on the full S4 x D4 state group.
    fixed_accept = []
    fixed_kernel = []
    fixed_accept_by_eta = {0: 0, 1: 0}
    for source in s4:
        for target in d4:
            passes = cycle_type(source) == cycle_type(target)
            sign_gate = perm_parity(source) == perm_parity(target)
            if passes:
                require(sign_gate, "cycle type did not preserve sign")
                fixed_accept.append((source, target))
                fixed_accept_by_eta[eta[target]] += 1
            if sign_gate:
                fixed_kernel.append((source, target))
    require(len(fixed_accept) == 34, "fixed-drift accept count")
    require(len(fixed_kernel) == 96, "fixed-drift sign-kernel count")
    require(fixed_accept_by_eta == {0: 16, 1: 18}, "fixed eta split")

    A3 = (1, 3, 2, 0)
    B3 = (1, 3, 0, 2)
    C3 = (3, 1, 0, 2)
    fixed_generators = {"A": (A3, G), "B": (B3, G), "C": (C3, T)}
    fixed_sign_values = {
        letter: perm_parity(source) ^ perm_parity(target)
        for letter, (source, target) in fixed_generators.items()
    }
    fixed_eta_values = {
        letter: eta[target] for letter, (_, target) in fixed_generators.items()
    }
    require(fixed_sign_values == {"A": 1, "B": 0, "C": 0}, "fixed A parity")
    require(fixed_eta_values == {"A": 1, "B": 1, "C": 1}, "fixed length parity")
    require(cycle_type(B3) == cycle_type(G), "B positive control")
    require(
        perm_parity(C3) == perm_parity(T) and cycle_type(C3) != cycle_type(T),
        "C sign-gate hostile",
    )

    marginal_rank = marginal_map_rank(13, 13)
    marginal_kernel = 13 * 13 - marginal_rank
    require(marginal_rank == 25 and marginal_kernel == 144, "13x13 marginal kernel")
    checkerboard = [[0] * 13 for _ in range(13)]
    checkerboard[0][0] = 1
    checkerboard[0][1] = -1
    checkerboard[1][0] = -1
    checkerboard[1][1] = 1
    row_sums = [sum(row) % 13 for row in checkerboard]
    col_sums = [sum(checkerboard[i][j] for i in range(13)) % 13 for j in range(13)]
    require(not all(entry == 0 for row in checkerboard for entry in row), "checkerboard zero")
    require(row_sums == [0] * 13 and col_sums == [0] * 13, "checkerboard marginals")

    # A hom C_m -> C_n is determined by x with m*x=0 mod n.
    hom_2_to_13 = [x for x in range(13) if (2 * x) % 13 == 0]
    hom_13_to_2 = [x for x in range(2) if (13 * x) % 2 == 0]
    require(hom_2_to_13 == [0] and hom_13_to_2 == [0], "coprime coefficient no-go")

    print("TYPED H1 / COBOUNDARY / CHARACTER / MARGINAL SYNTHESIS")
    print("implementation=independent_standard_library")
    print()
    print("GRAPH PROFILES (connected graphs)")
    print(f"K4: V={k4['vertices']} E={k4['edges']} dim_B1={k4['B1']} beta1={k4['beta1']}")
    print(
        f"C6: V={c6['vertices']} E={c6['edges']} dim_B1={c6['B1']} "
        f"beta1={c6['beta1']} V4_dim_B1=10 V4_dim_H1=2"
    )
    print(f"C7: V={c7['vertices']} E={c7['edges']} dim_B1={c7['B1']} beta1={c7['beta1']}")
    print(
        f"role_graph: V={role['vertices']} E={role['edges']} "
        f"dim_B1={role['B1']} beta1={role['beta1']}"
    )
    print(f"all_K4_to_C7_vertex_pullbacks_exact={pullback_maps}")
    print()
    print("ORIENTED TETRAHEDRAL H1 TRANSPORT")
    print("K4: opposite-face map C0/constants -> H1 has rank 3; split_index=16")
    print("K4_COVARIANCE: pi*Omega=sign(pi)*Omega*pi on all 24 vertex permutations")
    print("K4_char2: quotient_rank=1 class=product_of_four_vertex_square_classes")
    print(
        "THM3494_n3_four_view_class_nonzero=True "
        f"Iz(P=2/9,Q=0)={iz_on_index_divisor}"
    )
    print("role_graph: opposite-face map rank=6; split_index=256; H1 plus bridge rank=7")
    print("role_graph_kernel: constants plus the leaf-only direction; bridge restores the leaf")
    print(
        f"role_charts: both_H1_wings_nonzero={both_wings_nonzero}/72 "
        f"bridge_nonzero={bridge_nonzero}/72"
    )
    print(f"role_H1_chart_digest={role_class_digest}")
    print(
        f"role_automorphisms={len(role_automorphisms)} "
        f"invariant_H1_dimension_mod13={invariant_cycle_dimension}"
    )
    print()
    print("CYCLIC SEAM -> DECK DEFECT")
    print(f"p=13,m=7,basis_checks={transgression_13}, kill=13, restore=14")
    print(f"p=2,m=6,V4_basis_checks={transgression_2}, kill=2, restore=3")
    print()
    print("BERGGREN GROUP CHARACTERS")
    print(f"|D4|={len(d4)} |Hom(D4,F2)|={len(chars)} dim_H1_group=2")
    print(f"restrictions_on_(p,q,r)={restriction_rows}")
    print(f"common_character_kernel_on_translations={common_translation_kernel}")
    print("nonzero_seam_r_is_invisible_to_every_D4_character=True")
    print()
    print("WORD AUTOMATA")
    print(
        f"variable: pair_checks={variable_formula_checks} states={len(variable_states)} "
        f"character_kernel={len(variable_kernel_states)} accepted={len(variable_accept_states)}"
    )
    print(f"variable_character_on_(A,B,C)={a_parity_values}")
    print("variable_C_hostile=character_zero_but_rejected")
    print(
        f"fixed: states={len(s4) * len(d4)} sign_kernel={len(fixed_kernel)} "
        f"accepted={len(fixed_accept)} eta_split={fixed_accept_by_eta}"
    )
    print(f"fixed_sign_obstruction_on_(A,B,C)={fixed_sign_values}")
    print(f"fixed_period_character_on_(A,B,C)={fixed_eta_values}")
    print("fixed_C_hostile=sign_gate_passes_but_cycle_type_fails")
    print()
    print("MARGINAL TENSOR KERNEL")
    print(f"13x13_joint_dim=169 marginal_image_dim={marginal_rank} kernel_dim={marginal_kernel}")
    print("two_by_two_checkerboard_nonzero_with_zero_marginals=True")
    print()
    print("COEFFICIENT HOSTILE")
    print(
        f"nonzero_Hom(C2,C13)={len(hom_2_to_13) - 1} "
        f"nonzero_Hom(C13,C2)={len(hom_13_to_2) - 1}"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
