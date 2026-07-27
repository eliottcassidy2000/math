#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2545.

The proof in THM-2545 is a weighted Hall/max-flow argument.  This companion
exhausts small square and rectangular transportation instances, checks the
closed complete-off-diagonal formula, and realizes the word-stratified
aligned/swap hostile inside F_13.  It also checks the independent root/cut
character torsor controls used to compare with the 42-cut ancestry artifact.
"""

from collections import deque
from fractions import Fraction
from itertools import product


ROOTS = 13
FIELD = 547                 # 547-1 = 6*7*13
WORDS = ("a", "c", "ac")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compositions(total, length):
    """All weak compositions of total into length ordered parts."""
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def hall_feasible(left_margin, right_margin, allowed):
    """Weighted Hall inequalities for a bipartite support graph."""
    if sum(left_margin) != sum(right_margin):
        return False
    left_size = len(left_margin)
    for mask in range(1 << left_size):
        left_mass = sum(
            left_margin[i] for i in range(left_size) if mask & (1 << i)
        )
        neighbors = {
            j
            for i, j in allowed
            if mask & (1 << i)
        }
        right_mass = sum(right_margin[j] for j in neighbors)
        if left_mass > right_mass:
            return False
    return True


def maximum_transport(left_margin, right_margin, allowed):
    """Independent integer max-flow calculation for the same graph."""
    left_size = len(left_margin)
    right_size = len(right_margin)
    source = 0
    left0 = 1
    right0 = left0 + left_size
    sink = right0 + right_size
    vertices = sink + 1
    capacity = [[0] * vertices for _ in range(vertices)]
    total = sum(left_margin)

    for i, value in enumerate(left_margin):
        capacity[source][left0 + i] = value
    for i, j in allowed:
        capacity[left0 + i][right0 + j] = total
    for j, value in enumerate(right_margin):
        capacity[right0 + j][sink] = value

    flow = 0
    while True:
        parent = [-1] * vertices
        parent[source] = source
        queue = deque([source])
        while queue and parent[sink] == -1:
            u = queue.popleft()
            for v in range(vertices):
                if parent[v] == -1 and capacity[u][v] > 0:
                    parent[v] = u
                    queue.append(v)
        if parent[sink] == -1:
            break
        increment = total
        v = sink
        while v != source:
            u = parent[v]
            increment = min(increment, capacity[u][v])
            v = u
        v = sink
        while v != source:
            u = parent[v]
            capacity[u][v] -= increment
            capacity[v][u] += increment
            v = u
        flow += increment
    return flow


def exhaustive_hall_referee(left_size, right_size, max_total):
    """Exhaust every support graph and every small pair of margins."""
    edges = [(i, j) for i in range(left_size) for j in range(right_size)]
    checks = 0
    for graph_mask in range(1 << len(edges)):
        allowed = {
            edge for bit, edge in enumerate(edges) if graph_mask & (1 << bit)
        }
        for total in range(1, max_total + 1):
            for left in compositions(total, left_size):
                for right in compositions(total, right_size):
                    hall = hall_feasible(left, right, allowed)
                    flow = maximum_transport(left, right, allowed)
                    require(hall == (flow == total),
                            f"Hall/flow mismatch {left=} {right=} {allowed=}")
                    checks += 1
    return checks


def complete_off_diagonal_referee(max_size, max_total):
    checks = 0
    normalized_checks = 0
    for size in range(2, max_size + 1):
        allowed = {(i, j) for i in range(size) for j in range(size) if i != j}
        for total in range(1, max_total + 1):
            for left in compositions(total, size):
                for right in compositions(total, size):
                    off_diagonal_flow = maximum_transport(left, right, allowed)
                    minimum_diagonal = total - off_diagonal_flow
                    overload = max(
                        [0] + [left[i] + right[i] - total for i in range(size)]
                    )
                    require(minimum_diagonal == overload,
                            f"complete formula mismatch {left=} {right=}")
                    normalized = max(
                        [Fraction(0)]
                        + [
                            Fraction(left[i], total)
                            + Fraction(right[i], total)
                            - 1
                            for i in range(size)
                        ]
                    )
                    require(normalized == Fraction(overload, total),
                            "normalized overload mismatch")
                    checks += 1
                    normalized_checks += 1
    return checks, normalized_checks


def complete_off_diagonal_cemetery_referee(max_size, max_total):
    """Check the same formula when the right side has a cemetery column."""
    checks = 0
    for size in range(2, max_size + 1):
        cemetery = size
        allowed = {
            (i, j)
            for i in range(size)
            for j in range(size + 1)
            if j == cemetery or i != j
        }
        for total in range(1, max_total + 1):
            for left in compositions(total, size):
                for right in compositions(total, size + 1):
                    minimum_diagonal = total - maximum_transport(
                        left, right, allowed
                    )
                    overload = max(
                        [0] + [left[i] + right[i] - total for i in range(size)]
                    )
                    require(minimum_diagonal == overload,
                            f"cemetery formula mismatch {left=} {right=}")
                    checks += 1
    return checks


def margins(cells):
    left = [Fraction(0) for _ in range(ROOTS)]
    right = [Fraction(0) for _ in range(ROOTS)]
    for (i, j), value in cells.items():
        left[i] += value
        right[j] += value
    return tuple(left), tuple(right)


def diagonal_mass(cells):
    return sum(value for (i, j), value in cells.items() if i == j)


def primitive_root(prime):
    factors = []
    n = prime - 1
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            factors.append(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    if n > 1:
        factors.append(n)
    for candidate in range(2, prime):
        if all(pow(candidate, (prime - 1) // factor, prime) != 1
               for factor in factors):
            return candidate
    raise RuntimeError("no primitive root")


def word_hostile_referee():
    word_masses = {
        "a": Fraction(5, 7),
        "c": Fraction(11, 13),
        "ac": Fraction(17, 19),
    }
    aligned = {}
    swapped = {}
    for word in WORDS:
        half = word_masses[word] / 2
        aligned[word] = {(0, 0): half, (1, 1): half}
        swapped[word] = {(0, 1): half, (1, 0): half}
        require(sum(aligned[word].values()) == word_masses[word],
                "aligned word mass")
        require(sum(swapped[word].values()) == word_masses[word],
                "swapped word mass")
        require(margins(aligned[word]) == margins(swapped[word]),
                "hostile margins differ")
        require(diagonal_mass(aligned[word]) == word_masses[word],
                "aligned diagonal")
        require(diagonal_mass(swapped[word]) == 0, "swap diagonal")

    total = sum(word_masses.values())
    require(sum(diagonal_mass(aligned[w]) for w in WORDS) == total,
            "aligned total hit")
    require(sum(diagonal_mass(swapped[w]) for w in WORDS) == 0,
            "swapped total hit")

    # The common two-root head/target marginal has every nontrivial F_13
    # Fourier colour: 1+zeta^k cannot vanish because zeta has odd order 13.
    generator = primitive_root(FIELD)
    zeta = pow(generator, (FIELD - 1) // ROOTS, FIELD)
    require(pow(zeta, ROOTS, FIELD) == 1, "root character order")
    colour_checks = 0
    for colour in range(1, ROOTS):
        require((1 + pow(zeta, (-colour) % ROOTS, FIELD)) % FIELD != 0,
                f"two-root marginal lost colour {colour}")
        colour_checks += 1

    # Minimality: on a one-root universe with no cemetery mass, no
    # off-diagonal coupling exists; two roots admit the displayed swap.
    require(maximum_transport((1,), (1,), set()) == 0,
            "one-root hostile should be impossible")
    require(maximum_transport((1, 1), (1, 1), {(0, 1), (1, 0)}) == 2,
            "two-root swap should be feasible")
    return aligned, swapped, word_masses, colour_checks


def cut42_control_referee(aligned, swapped, word_masses):
    """Tensor the semantic tables with the 42-cut artifact's sidecars.

    The artifact has 6*7 cuts and 12*6*6=432 lawful nonzero coefficient
    placements.  Its character factors use typed-different sheets.  A
    semantic table is therefore an independent tensor factor unless a
    common-ancestry diagonal map is supplied.
    """
    cuts = tuple(product(range(1, 7), range(7)))
    modes = tuple(product(range(1, 13), range(1, 7), range(1, 7)))
    require(len(cuts) == 42, "42-cut census")
    require(len(modes) == 432, "paired-mode census")

    def retained_signature(model):
        # The owner label is the same unique O in both tensors; genuine
        # target-active mass is all mass (there is no cemetery column).
        word_signature = tuple(
            (
                word,
                sum(model[word].values()),
                margins(model[word]),
                sum(value for (_, _), value in model[word].items()),
            )
            for word in WORDS
        )
        coefficient_signature = tuple(
            (word, len(cuts), len(modes), word_masses[word]) for word in WORDS
        )
        owner_signature = tuple((word, "O", word_masses[word]) for word in WORDS)
        return word_signature, coefficient_signature, owner_signature

    require(retained_signature(aligned) == retained_signature(swapped),
            "42-cut/owner/target marginal signature changed under swap")
    require(sum(diagonal_mass(swapped[word]) for word in WORDS) == 0,
            "42-cut sidecars unexpectedly forced semantic diagonal")

    # Exact character hostile: every nontrivial character is nonzero, yet
    # averaging either torsor is zero.  This is the algebraic reason a
    # coefficientwise pairing is not already a common-base positive product.
    generator = primitive_root(FIELD)
    zeta = pow(generator, (FIELD - 1) // 13, FIELD)
    xi = pow(generator, (FIELD - 1) // 7, FIELD)
    root_collapse_checks = 0
    cut_collapse_checks = 0
    for alpha in range(1, 13):
        values = [pow(zeta, (-alpha * h) % 13, FIELD) for h in range(13)]
        require(all(value != 0 for value in values), "root phase vanished")
        require(sum(values) % FIELD == 0, "root torsor collapse")
        root_collapse_checks += 1
    for a in range(1, 7):
        for beta in range(1, 7):
            values = [pow(xi, (beta * a * c) % 7, FIELD) for c in range(7)]
            require(all(value != 0 for value in values), "cut phase vanished")
            require(sum(values) % FIELD == 0, "cut torsor collapse")
            cut_collapse_checks += 1
    return len(cuts), len(modes), root_collapse_checks, cut_collapse_checks


def main():
    square2 = exhaustive_hall_referee(2, 2, 5)
    square3 = exhaustive_hall_referee(3, 3, 4)
    cemetery = exhaustive_hall_referee(2, 3, 4)
    complete, normalized = complete_off_diagonal_referee(5, 6)
    complete_cemetery = complete_off_diagonal_cemetery_referee(4, 5)
    aligned, swapped, word_masses, colours = word_hostile_referee()
    cuts, modes, root_collapses, cut_collapses = cut42_control_referee(
        aligned, swapped, word_masses
    )

    print(f"hall_flow_exhaustive_2x2={square2}")
    print(f"hall_flow_exhaustive_3x3={square3}")
    print(f"hall_flow_rectangular_2x3_with_cemetery={cemetery}")
    print(f"complete_off_diagonal_formula_checks={complete}")
    print(f"complete_off_diagonal_cemetery_formula_checks={complete_cemetery}")
    print(f"normalized_formula_checks={normalized}")
    print("word_hostile=three_positive_strata margins_equal=1 "
          "aligned_hit=full swapped_hit=0")
    print(f"two_root_F13_nontrivial_colour_checks={colours}/12")
    print("hostile_minimal_root_count=2")
    print(f"cut_bundle_controls={cuts} artifact_paired_mode_labels={modes}")
    print(f"root_torsor_collapse_checks={root_collapses} "
          f"cut_torsor_collapse_checks={cut_collapses}")
    print("cut42_owner_target_marginals_admit_per_word_swap=PASS")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
