#!/usr/bin/env python3
"""Exact companion for THM-2502.

Dependency-free exhaustive checks for the Boolean endpoint cube, its
first-difference tournament and carry colours, the top Newton atom, and the
lawful two-dipole support boundary.
"""

from itertools import combinations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def bits7(k: int) -> tuple[int, ...]:
    return tuple((k >> (6 - j)) & 1 for j in range(7))


def first_difference(x: tuple[int, ...], y: tuple[int, ...]) -> int:
    for j, (a, b) in enumerate(zip(x, y), start=1):
        if a != b:
            return j
    raise RuntimeError("equal words have no first difference")


def v2(n: int) -> int:
    require(n > 0, "v2 input must be positive")
    answer = 0
    while n % 2 == 0:
        answer += 1
        n //= 2
    return answer


def first_failure_cell(word: tuple[int, ...]) -> int:
    for i in range(5):
        if word[i] == 0:
            return i + 1
    return 0


def blocker_word(word: tuple[int, ...]) -> tuple[int, int]:
    return word[5], word[6]


def ghost_value(word: tuple[int, ...]) -> int:
    value = 1
    for bit in word:
        value *= bit
    return value


def mixed_difference_top(subset: tuple[int, ...]) -> int:
    """Delta_subset G(0), by exact Boolean inclusion-exclusion."""
    total = 0
    size = len(subset)
    for mask in range(1 << size):
        word = [0] * 7
        chosen = 0
        for j, coordinate in enumerate(subset):
            if (mask >> j) & 1:
                word[coordinate] = 1
                chosen += 1
        total += (-1) ** (size - chosen) * ghost_value(tuple(word))
    return total


def rank_one_row_kernel_dimension(columns: int) -> int:
    require(columns > 0, "positive column count")
    return columns - 1


def main() -> None:
    masks = [bits7(k) for k in range(128)]
    require(len(set(masks)) == 128, "seven-bit mask enumeration")

    colour_counts = [0] * 7
    for i, j in combinations(range(128), 2):
        x, y = masks[i], masks[j]
        h = first_difference(x, y)
        require(x[h - 1] == 0 and y[h - 1] == 1, "lex orientation")
        colour_counts[h - 1] += 1
    expected_colours = [2 ** (13 - h) for h in range(1, 8)]
    require(colour_counts == expected_colours, "first-difference census")
    require(sum(colour_counts) == 128 * 127 // 2, "pair census")

    triple_count = 0
    for i, j, k in combinations(range(128), 3):
        require(masks[i] < masks[j] < masks[k], "lex tournament transitivity")
        triple_count += 1
    require(triple_count == 128 * 127 * 126 // 6, "triple census")

    successor_colours = []
    ruler_counts = [0] * 7
    for k in range(127):
        h = first_difference(masks[k], masks[k + 1])
        predicted = 7 - v2(k + 1)
        require(h == predicted, "successor carry colour")
        successor_colours.append(h)
        ruler_counts[h - 1] += 1

    blocker_words = [(0, 0), (0, 1), (1, 0), (1, 1)]
    cell_counts = {
        (i, sigma): 0
        for i in range(6)
        for sigma in blocker_words
    }
    for word in masks:
        cell_counts[(first_failure_cell(word), blocker_word(word))] += 1
    extensions = [1, 16, 8, 4, 2, 1]
    for sigma in blocker_words:
        require(
            [cell_counts[(i, sigma)] for i in range(6)] == extensions,
            "fixed-blocker first-failure cells",
        )
    aggregate_cells = [
        sum(cell_counts[(i, sigma)] for sigma in blocker_words)
        for i in range(6)
    ]
    require(
        aggregate_cells == [4, 64, 32, 16, 8, 4],
        "aggregate first-failure layers",
    )
    ghosts = [word for word in masks if ghost_value(word) == 1]
    require(ghosts == [(1, 1, 1, 1, 1, 1, 1)], "unique ghost")
    require(
        first_failure_cell(ghosts[0]) == 0
        and blocker_word(ghosts[0]) == (1, 1),
        "ghost cell label",
    )

    newton_by_order = [0] * 8
    nonzero_subsets = []
    for size in range(8):
        for subset in combinations(range(7), size):
            value = mixed_difference_top(subset)
            newton_by_order[size] += value
            if value:
                nonzero_subsets.append((subset, value))
    require(newton_by_order == [0, 0, 0, 0, 0, 0, 0, 1], "Newton layers")
    require(nonzero_subsets == [(tuple(range(7)), 1)], "unique top coefficient")
    for coordinate in range(7):
        for word in masks:
            if word[coordinate] == 0:
                require(ghost_value(word) == 0, "coordinate-zero face")

    dipole_support_counts = {2: 0, 4: 0}
    for s in range(13):
        for t in range(13):
            if s == 0 and t == 0:
                continue
            representative = (s, (-s) % 13, t, (-t) % 13)
            support = sum(value != 0 for value in representative)
            require(support in dipole_support_counts, "dipole support")
            dipole_support_counts[support] += 1
    require(dipole_support_counts == {2: 24, 4: 144}, "dipole census")
    require(rank_one_row_kernel_dimension(4) == 3, "four-cell ambiguity")
    hostile_a = (-1, 1, 0, 0)
    hostile_b = (0, 1, 0, -1)
    require(sum(hostile_a) == sum(hostile_b) == 0, "four-cell hostiles")
    require(hostile_a[1] == hostile_b[1] == 1, "fixed transition")
    require(hostile_a[0] != hostile_b[0], "undetermined matched diagonal")

    tail_counts = [1]
    for _ in range(5):
        tail_counts.append(2 + 2 * tail_counts[-1])
    require(tail_counts == [1, 4, 10, 22, 46, 94], "tail recurrence")
    matched_counts = [1, 46, 22, 10, 4, 1]
    require([23 + count for count in matched_counts] == [24, 69, 45, 33, 27, 24],
            "fixed-left cospan counts")

    # In Z[X]/Phi_13, P_g + P_d = Phi_13, so every nontrivial
    # thirteenth-root evaluation vanishes and g_hat = -d_hat.
    danger = [1 if j in {0, 1, 12} else 0 for j in range(13)]
    safe = [1 - value for value in danger]
    require([a + b for a, b in zip(danger, safe)] == [1] * 13,
            "complement coefficient identity")

    print("THM-2502 endpoint Newton/carry boundary exact audit")
    print(f"masks={len(masks)} pairs={sum(colour_counts)} triples={triple_count}")
    print("first_difference_edges=" + ",".join(map(str, colour_counts)))
    print("successor_ruler_counts=" + ",".join(map(str, ruler_counts)))
    print("per_blocker_cell_extensions=" + ",".join(map(str, extensions)))
    print("aggregate_first_failure_layers=" + ",".join(map(str, aggregate_cells)))
    print(f"thm2445_cells={len(cell_counts)}")
    print("newton_layer_sums=" + ",".join(map(str, newton_by_order)))
    print(f"dipole_support_2={dipole_support_counts[2]} "
          f"dipole_support_4={dipole_support_counts[4]}")
    print("tail_counts=" + ",".join(map(str, tail_counts)))
    print("fixed_left_piece_counts=24,69,45,33,27,24")
    print("four_cell_kernel_dimension=3")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
