#!/usr/bin/env python3
"""Exact finite checks for THM-2315.

All load-bearing tests use explicit exceptions, so optimized mode performs
the same audit.  The span and marginal arguments are symbolic; their minimal
finite hostile controls are checked here.
"""

from itertools import combinations, product


P = 13
INF = "inf"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def inv(value):
    return pow(value % P, -1, P)


def projective_points():
    return list(range(P)) + [INF]


def vector(point):
    if point == INF:
        return (0, 1)
    return (1, point)


def support_word(point):
    if point == 0:
        return "a"
    if point == INF:
        return "b"
    return "ab"


def torus_action(c, point):
    if point == INF:
        return INF
    return (c * point) % P


def swap_action(point):
    if point == 0:
        return INF
    if point == INF:
        return 0
    return inv(point)


def antidiagonal_action(c, point):
    swapped = swap_action(point)
    return torus_action(c, swapped)


def permutation(action):
    points = projective_points()
    return tuple(action(point) for point in points)


def chi(value):
    value %= P
    if value == 0:
        return 0
    return 1 if pow(value, (P - 1) // 2, P) == 1 else -1


def determinant(left, right):
    x1, y1 = vector(left)
    x2, y2 = vector(right)
    return (x1 * y2 - y1 * x2) % P


def marginal(word_mass):
    p_a, p_b, fork = word_mass
    return (p_a + fork, p_b + fork)


def reconstruct_mass(total, marginals):
    m_a, m_b = marginals
    return (total - m_b, total - m_a, m_a + m_b - total)


def hafnian4(matrix):
    return (
        matrix[0][1] * matrix[2][3]
        + matrix[0][2] * matrix[1][3]
        + matrix[0][3] * matrix[1][2]
    )


def star_with_context(m_a, m_b, context):
    """Vertex order j,a,b,z; context is 'read_a' or 'read_b'."""
    matrix = [[0] * 4 for _ in range(4)]

    def edge(i, j, weight):
        matrix[i][j] = weight
        matrix[j][i] = weight

    edge(0, 1, m_a)
    edge(0, 2, m_b)
    if context == "read_a":
        edge(2, 3, 1)
    elif context == "read_b":
        edge(1, 3, 1)
    else:
        raise RuntimeError("unknown context")
    return matrix


def check_projective_carrier():
    points = projective_points()
    require(len(points) == 14, "wrong projective-line size")
    counts = {"a": 0, "b": 0, "ab": 0}
    for point in points:
        counts[support_word(point)] += 1
    require(counts == {"a": 1, "b": 1, "ab": 12}, "wrong support fibres")

    torus_permutations = {
        permutation(lambda point, c=c: torus_action(c, point))
        for c in range(1, P)
    }
    normalizer_permutations = set(torus_permutations)
    normalizer_permutations.update(
        permutation(lambda point, c=c: antidiagonal_action(c, point))
        for c in range(1, P)
    )
    require(len(torus_permutations) == 12, "wrong labelled-boundary torus")
    require(len(normalizer_permutations) == 24, "wrong boundary normalizer")

    mixed = set(range(1, P))
    orbit = {torus_action(c, 1) for c in range(1, P)}
    require(orbit == mixed, "mixed gains are not one torus orbit")
    fixed_by_all = {
        gain
        for gain in mixed
        if all(torus_action(c, gain) == gain for c in range(1, P))
    }
    require(not fixed_by_all, "a canonical torus-fixed fork gain exists")

    for gain in mixed:
        require(swap_action(swap_action(gain)) == gain, "swap is not involutive")
    fixed_inversion = {gain for gain in mixed if swap_action(gain) == gain}
    require(fixed_inversion == {1, 12}, "wrong inversion-fixed gains")

    for left, right in combinations(points, 2):
        require(determinant(left, right) != 0, "projective pair is dependent")
    return counts, len(torus_permutations), len(normalizer_permutations), fixed_inversion


def check_head_selectors_and_characters():
    mixed = set(range(1, P))
    fixed = {gain for gain in mixed if inv(gain) == gain}
    remaining = mixed - fixed
    inverse_pairs = []
    while remaining:
        gain = min(remaining)
        partner = inv(gain)
        inverse_pairs.append((gain, partner))
        remaining.remove(gain)
        remaining.remove(partner)
    require(len(inverse_pairs) == 5, "wrong inverse-pair count")

    selectors = 0
    for choices in product((0, 1), repeat=len(inverse_pairs)):
        head = {1: "tie", 12: "tie"}
        for choice, (gain, partner) in zip(choices, inverse_pairs):
            head[gain] = "a" if choice == 0 else "b"
            head[partner] = "b" if choice == 0 else "a"
        for gain in range(1, P):
            swapped = {"a": "b", "b": "a", "tie": "tie"}[head[gain]]
            require(head[inv(gain)] == swapped, "head selector breaks swap")
        require(sum(value == "tie" for value in head.values()) == 2, "extra tie")
        selectors += 1
    require(selectors == 32, "wrong minimal-tie selector count")

    squares = {gain for gain in range(1, P) if chi(gain) == 1}
    nonsquares = set(range(1, P)) - squares
    require(len(squares) == len(nonsquares) == 6, "wrong character classes")
    require(all(chi(inv(gain)) == chi(gain) for gain in range(1, P)), "chi inversion")
    require(chi(1) == chi(4) == 1, "signed-gain collision failed")
    require((1 - 1) % P == 0 and (4 - 1) % P != 0, "linear selector failed")

    for left, right in combinations(projective_points(), 2):
        det_lr = determinant(left, right)
        det_rl = determinant(right, left)
        require(det_rl == (-det_lr) % P, "determinant is not alternating")
        require(chi(det_lr) == chi(det_rl), "chi_13 unexpectedly orients")
    return len(inverse_pairs), selectors, len(squares), len(nonsquares)


def check_marginal_and_hafnian_boundary():
    two_pure = (1, 1, 0)
    one_fork = (0, 0, 1)
    require(marginal(two_pure) == marginal(one_fork) == (1, 1), "marginal collision")

    kernel_direction = (-1, -1, 1)
    require(marginal(kernel_direction) == (0, 0), "wrong marginal kernel")

    checked = 0
    for p_a, p_b, fork in product(range(6), repeat=3):
        masses = (p_a, p_b, fork)
        total = sum(masses)
        require(reconstruct_mass(total, marginal(masses)) == masses, "mass repair")
        checked += 1
    require(checked == 216, "wrong mass-repair census")

    m_a, m_b = marginal(one_fork)
    read_a = hafnian4(star_with_context(m_a, m_b, "read_a"))
    read_b = hafnian4(star_with_context(m_a, m_b, "read_b"))
    require((read_a, read_b) == (m_a, m_b), "forced context misses marginals")

    require((1 * 2) % P == (3 * 5) % P == 2, "pair-product collision failed")
    gain_12 = (2 * inv(1)) % P
    gain_35 = (5 * inv(3)) % P
    require((gain_12, gain_35) == (2, 6), "wrong collided gains")
    return marginal(two_pure), checked, (read_a, read_b), (gain_12, gain_35)


def check_composition_boundary():
    universe = {1, 2, 3, 4}
    incoming_arrival = {1, 2}
    outgoing_overlap = {1, 2}
    outgoing_disjoint = {3, 4}
    require(len(incoming_arrival) == len(outgoing_overlap) == len(outgoing_disjoint) == 2,
            "edge masses differ")
    positive_pullback = incoming_arrival & outgoing_overlap
    empty_pullback = incoming_arrival & outgoing_disjoint
    require(len(positive_pullback) == 2, "positive pullback has wrong size")
    require(not empty_pullback, "disjoint pullback is nonempty")
    return len(universe), len(positive_pullback), len(empty_pullback)


def main():
    support_counts, torus_size, normalizer_size, fixed = check_projective_carrier()
    inverse_pairs, selectors, squares, nonsquares = check_head_selectors_and_characters()
    marginal_collision, mass_census, forced_reads, collided_gains = (
        check_marginal_and_hafnian_boundary()
    )
    universe_size, positive_size, empty_size = check_composition_boundary()

    print("THM-2315 exact companion")
    print(f"projective_support_counts={support_counts}")
    print(f"labelled_boundary_torus_size={torus_size}")
    print(f"unordered_boundary_normalizer_size={normalizer_size}")
    print(f"inversion_fixed_gains={sorted(fixed)}")
    print(f"inverse_gain_pairs={inverse_pairs}")
    print(f"minimal_tie_head_selectors={selectors}")
    print(f"quadratic_character_classes={squares}+{nonsquares}")
    print(f"pair_marginal_collision={marginal_collision}")
    print(f"mass_repair_census={mass_census}")
    print(f"forced_hafnian_reads={forced_reads}")
    print(f"pair_product_collided_gains={collided_gains}")
    print(f"composition_universe_size={universe_size}")
    print(f"positive_vs_empty_pullback={positive_size},{empty_size}")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
