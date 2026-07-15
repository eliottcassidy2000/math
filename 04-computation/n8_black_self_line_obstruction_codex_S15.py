#!/usr/bin/env python3
"""Exact obstruction to the proposed all-n black self-line law (THM-849).

For a fixed-path tiling t, kappa reverses every non-path arc and sigma is
staircase reflection.  This script counts

    Q_K(n) = #{t : sigma(t) != t and T(kappa(t)) isomorphic to T(t)}.

Thus Q_K(n) = 2*selfK(n), where selfK is the number of ordinary black
complement self-lines.  Isomorphism is decided by common directed equitable
refinement followed by exact individualization/backtracking.  An independent
signed cycle-index calculation supplies the ordinary and self-converse class
counts, including T(8)=6880 and SC(8)=176.

Tournament Analysis deliberately uses the four Klein operations as vertices,
not runners or arcs.  On every black quasi-fixed orbit, H is tied on all four
vertices; lexicographic endpoint order is the tie gauge and gives the
transitive four-vertex tournament.  This records that the quotient preserves
the Klein sheet but destroys the class/path incidence needed by Q_K.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import permutations
from math import factorial
from time import perf_counter


def staircase(n: int) -> tuple[list[tuple[int, int]], list[int]]:
    """Return zero-based staircase cells and their reflection permutation."""
    cells = [
        (x, y)
        for y in range(n - 2)
        for x in range(n - 1, y + 1, -1)
        if x - y >= 2
    ]
    index = {cell: i for i, cell in enumerate(cells)}
    reflection = [index[(n - 1 - y, n - 1 - x)] for x, y in cells]
    return cells, reflection


def tournament_rows(n: int, cells: list[tuple[int, int]], word: int) -> tuple[int, ...]:
    """Build row bitsets; the observer path is n-1 -> ... -> 0."""
    rows = [0] * n
    for x in range(1, n):
        rows[x] |= 1 << (x - 1)
    for i, (x, y) in enumerate(cells):
        if word >> i & 1:
            rows[x] |= 1 << y
        else:
            rows[y] |= 1 << x
    return tuple(rows)


def paired_equitable_colors(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[list[int], list[int]] | None:
    """Refine both tournaments in one color namespace, or reject them."""
    n = len(left)
    left_colors = [row.bit_count() for row in left]
    right_colors = [row.bit_count() for row in right]
    if Counter(left_colors) != Counter(right_colors):
        return None

    while True:
        signatures: list[tuple[int, tuple[int, ...], tuple[int, ...]]] = []
        for rows, colors in ((left, left_colors), (right, right_colors)):
            for vertex in range(n):
                out_colors = tuple(
                    sorted(colors[w] for w in range(n) if rows[vertex] >> w & 1)
                )
                in_colors = tuple(
                    sorted(colors[w] for w in range(n) if rows[w] >> vertex & 1)
                )
                signatures.append((colors[vertex], out_colors, in_colors))

        color_of = {signature: i for i, signature in enumerate(sorted(set(signatures)))}
        new_left = [color_of[signatures[v]] for v in range(n)]
        new_right = [color_of[signatures[n + v]] for v in range(n)]
        if Counter(new_left) != Counter(new_right):
            return None

        old_class_count = len(set(left_colors + right_colors))
        new_class_count = len(set(new_left + new_right))
        left_colors, right_colors = new_left, new_right
        if new_class_count == old_class_count:
            return left_colors, right_colors


def isomorphism_count(
    left: tuple[int, ...], right: tuple[int, ...], stop_after_one: bool = False
) -> int:
    """Count all directed isomorphisms, with an optional existence shortcut."""
    colors = paired_equitable_colors(left, right)
    if colors is None:
        return 0
    left_colors, right_colors = colors
    n = len(left)
    candidates = {
        v: [w for w in range(n) if right_colors[w] == left_colors[v]]
        for v in range(n)
    }
    image = [-1] * n

    def search(used: int, assigned: int) -> int:
        if assigned == n:
            return 1

        best_vertex = -1
        best_options: list[int] | None = None
        for v in range(n):
            if image[v] >= 0:
                continue
            options = []
            for w in candidates[v]:
                if used >> w & 1:
                    continue
                if all(
                    image[u] < 0
                    or ((left[v] >> u) & 1) == ((right[w] >> image[u]) & 1)
                    for u in range(n)
                ):
                    options.append(w)
            if not options:
                return 0
            if best_options is None or len(options) < len(best_options):
                best_vertex, best_options = v, options

        assert best_options is not None
        total = 0
        for w in best_options:
            image[best_vertex] = w
            total += search(used | (1 << w), assigned + 1)
            image[best_vertex] = -1
            if stop_after_one and total:
                return 1
        return total

    return search(0, 0)


def brute_isomorphic(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    """Permutation oracle used only for the independent n<=6 audit."""
    n = len(left)
    left_degrees = [row.bit_count() for row in left]
    right_degrees = [row.bit_count() for row in right]
    for image in permutations(range(n)):
        if any(left_degrees[v] != right_degrees[image[v]] for v in range(n)):
            continue
        if all(
            ((left[v] >> u) & 1) == ((right[image[v]] >> image[u]) & 1)
            for v in range(n)
            for u in range(n)
            if u != v
        ):
            return True
    return False


def sigma_word(word: int, reflection: list[int]) -> int:
    return sum(((word >> reflection[i]) & 1) << i for i in range(len(reflection)))


def endpoint_face_word(
    n: int,
    cells: list[tuple[int, int]],
    word: int,
    remove_low: bool,
) -> int:
    """Delete path endpoint 0 or n-1 and return the size-(n-1) face word."""
    lower_cells, _ = staircase(n - 1)
    lower_index = {cell: i for i, cell in enumerate(lower_cells)}
    face = 0
    for i, (x, y) in enumerate(cells):
        if remove_low:
            if y == 0:
                continue
            lower_cell = (x - 1, y - 1)
        else:
            if x == n - 1:
                continue
            lower_cell = (x, y)
        face |= ((word >> i) & 1) << lower_index[lower_cell]
    return face


def lower_line_type(n: int, word: int) -> str:
    """Classify a lower line as ordinary self, merged-only self, or neither."""
    cells, _ = staircase(n)
    full = (1 << len(cells)) - 1
    tournament = tournament_rows(n, cells, word)
    flipped = tournament_rows(n, cells, word ^ full)
    if isomorphism_count(tournament, flipped, stop_after_one=True):
        return "ordinary"
    off_diagonal = (1 << n) - 1
    converse = tuple(
        off_diagonal ^ (1 << vertex) ^ tournament[vertex]
        for vertex in range(n)
    )
    if isomorphism_count(flipped, converse, stop_after_one=True):
        return "merged_only"
    return "neither"


def cyclic_triangles(rows: tuple[int, ...]) -> int:
    n = len(rows)
    total = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                total += int(
                    ((rows[i] >> j) & 1)
                    and ((rows[j] >> k) & 1)
                    and ((rows[k] >> i) & 1)
                )
                total += int(
                    ((rows[i] >> k) & 1)
                    and ((rows[k] >> j) & 1)
                    and ((rows[j] >> i) & 1)
                )
    return total


def hamiltonian_paths(rows: tuple[int, ...]) -> int:
    n = len(rows)
    dp = {(1 << v, v): 1 for v in range(n)}
    for size in range(1, n):
        next_dp: dict[tuple[int, int], int] = {}
        for (mask, v), count in dp.items():
            if mask.bit_count() != size:
                continue
            for w in range(n):
                if not (mask >> w & 1) and rows[v] >> w & 1:
                    key = (mask | (1 << w), w)
                    next_dp[key] = next_dp.get(key, 0) + count
        dp.update(next_dp)
    full = (1 << n) - 1
    return sum(count for (mask, _), count in dp.items() if mask == full)


def signed_pair_burnside(n: int, anti: bool) -> int:
    """Ordinary (anti=False) or converse-twisted tournament Burnside count."""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    index = {pair: i for i, pair in enumerate(pairs)}
    fixed_sum = 0
    for permutation in permutations(range(n)):
        edge_image = []
        edge_sign = []
        for i, j in pairs:
            a, b = permutation[i], permutation[j]
            edge_image.append(index[tuple(sorted((a, b)))])
            edge_sign.append((a > b) ^ anti)

        seen = 0
        cycle_count = 0
        consistent = True
        for edge in range(len(pairs)):
            if seen >> edge & 1:
                continue
            cycle_count += 1
            parity = 0
            current = edge
            while not (seen >> current & 1):
                seen |= 1 << current
                parity ^= edge_sign[current]
                current = edge_image[current]
            if parity:
                consistent = False
                break
        if consistent:
            fixed_sum += 1 << cycle_count
    assert fixed_sum % factorial(n) == 0
    return fixed_sum // factorial(n)


def census(n: int, brute_oracle: bool) -> dict[str, object]:
    cells, reflection = staircase(n)
    dimension = len(cells)
    full = (1 << dimension) - 1
    qualifying = blue = black = 0
    score_survivors = 0
    realizer_pairs = blue_realizer_pairs = black_realizer_pairs = 0
    automorphisms: Counter[tuple[str, int]] = Counter()
    blue_orbits: set[int] = set()
    black_orbits: set[int] = set()
    black_orbit_words: dict[int, tuple[int, int, int, int]] = {}

    for word in range(1 << dimension):
        left = tournament_rows(n, cells, word)
        right = tournament_rows(n, cells, word ^ full)
        if sorted(row.bit_count() for row in left) != sorted(
            row.bit_count() for row in right
        ):
            continue
        score_survivors += 1
        exact = bool(isomorphism_count(left, right, stop_after_one=True))
        if brute_oracle and n <= 6:
            assert exact == brute_isomorphic(left, right)
        if not exact:
            continue

        aut = isomorphism_count(left, left)
        assert aut > 0
        reflected = sigma_word(word, reflection)
        is_blue = reflected == word
        orbit = (word, word ^ full, reflected, reflected ^ full)
        representative = min(orbit)

        qualifying += 1
        realizer_pairs += aut
        if is_blue:
            blue += 1
            blue_realizer_pairs += aut
            blue_orbits.add(representative)
            automorphisms[("blue", aut)] += 1
            assert len(set(orbit)) == 2
        else:
            black += 1
            black_realizer_pairs += aut
            black_orbits.add(representative)
            black_orbit_words.setdefault(representative, orbit)
            automorphisms[("black", aut)] += 1
            assert len(set(orbit)) == 4

    assert qualifying == blue + black
    assert blue % 2 == 0 and black % 4 == 0
    assert len(blue_orbits) == blue // 2
    assert len(black_orbits) == black // 4

    ta_ties = 0
    for orbit in black_orbit_words.values():
        rows = [tournament_rows(n, cells, word) for word in orbit]
        payloads = {
            (
                tuple(sorted(row.bit_count() for row in tournament)),
                cyclic_triangles(tournament),
                hamiltonian_paths(tournament),
            )
            for tournament in rows
        }
        if len(payloads) == 1:
            ta_ties += 1
    assert ta_ties == len(black_orbits)

    face_recursion: Counter[str] = Counter()
    if n == 8:
        for orbit in black_orbit_words.values():
            sheet_types = []
            for word in orbit:
                low = lower_line_type(
                    n - 1, endpoint_face_word(n, cells, word, remove_low=True)
                )
                high = lower_line_type(
                    n - 1, endpoint_face_word(n, cells, word, remove_low=False)
                )
                sheet_types.append((low, high))
            ordinary_per_sheet = [pair.count("ordinary") for pair in sheet_types]
            merged_only = sum(pair.count("merged_only") for pair in sheet_types)
            if ordinary_per_sheet == [1, 1, 1, 1] and merged_only == 0:
                face_recursion["one_ordinary_face_per_sheet"] += 1
            elif ordinary_per_sheet == [0, 0, 0, 0] and merged_only == 0:
                face_recursion["no_self_face"] += 1
            else:
                face_recursion["other"] += 1
        assert face_recursion == {
            "one_ordinary_face_per_sheet": 55,
            "no_self_face": 46,
        }

    return {
        "n": n,
        "dimension": dimension,
        "tilings": 1 << dimension,
        "classes": signed_pair_burnside(n, anti=False),
        "SC": signed_pair_burnside(n, anti=True),
        "score_survivors": score_survivors,
        "Q": qualifying,
        "Q_blue": blue,
        "Q_black": black,
        "selfB": blue // 2,
        "selfK": black // 2,
        "blue_V4_orbits": len(blue_orbits),
        "black_V4_orbits": len(black_orbits),
        "realizer_pairs": realizer_pairs,
        "blue_realizer_pairs": blue_realizer_pairs,
        "black_realizer_pairs": black_realizer_pairs,
        "automorphisms": automorphisms,
        "TA_tied_black_orbits": ta_ties,
        "face_recursion": face_recursion,
    }


def format_histogram(histogram: Counter[tuple[str, int]]) -> str:
    return ", ".join(
        f"{color}/Aut{aut}:{count}"
        for (color, aut), count in sorted(histogram.items())
    ) or "empty"


def main() -> None:
    if not __debug__:
        raise RuntimeError("this exact verifier requires Python assertions")
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=8, choices=range(3, 9))
    parser.add_argument(
        "--no-brute-oracle",
        action="store_true",
        help="skip the independent permutation oracle at n<=6",
    )
    args = parser.parse_args()

    print("THM-849: exact black self-line obstruction")
    print("Q_black = # non-grid-symmetric t with T(kappa t) ~= T(t) = 2*selfK")
    print()
    print(
        " n  dim  tilings  classes   SC  score-pass   Q  Qblue  Qblack  "
        "selfB  selfK  V4black  pair-sum"
    )
    rows = []
    started = perf_counter()
    for n in range(3, args.max_n + 1):
        row = census(n, brute_oracle=not args.no_brute_oracle)
        rows.append(row)
        print(
            f"{n:2d} {row['dimension']:4d} {row['tilings']:8d} "
            f"{row['classes']:8d} {row['SC']:4d} {row['score_survivors']:11d} "
            f"{row['Q']:3d} {row['Q_blue']:6d} {row['Q_black']:7d} "
            f"{row['selfB']:6d} {row['selfK']:6d} "
            f"{row['black_V4_orbits']:8d} {row['realizer_pairs']:9d}"
        )
        print(f"   automorphism histogram: {format_histogram(row['automorphisms'])}")

    expected = {
        3: (2, 2, 0, 0, 0),
        4: (4, 2, 2, 2, 0),
        5: (12, 8, 8, 0, 8),
        6: (56, 12, 16, 4, 12),
        7: (456, 88, 88, 0, 88),
        8: (6880, 176, 412, 8, 404),
    }
    for row in rows:
        assert (
            row["classes"],
            row["SC"],
            row["Q"],
            row["Q_blue"],
            row["Q_black"],
        ) == expected[row["n"]]

    if args.max_n == 8:
        row8 = rows[-1]
        print()
        print("n=8 obstruction:")
        print(f"  2*selfK(8) = Q_black(8) = {row8['Q_black']}")
        print(f"  SC(8) = {row8['SC']}")
        print(f"  defect = {row8['Q_black'] - row8['SC']} = classes(7)/2")
        assert row8["Q_black"] - row8["SC"] == rows[-2]["classes"] // 2
        print(
            "  black quasi-fixed set = 101 free Klein-four orbits; "
            "there are 202 ordinary black self-lines"
        )
        print(
            "  endpoint-face audit: 55/101 orbits have exactly one ordinary "
            "size-7 self-line face on every Klein sheet; 46/101 have neither "
            "face; merged-only faces = 0"
        )
        print(
            "  Tournament Analysis on {1,kappa,sigma,kappa*sigma}: H/score/C3 "
            "tie on all 101 black orbits; lex gauge has score histogram "
            "{0:1,1:1,2:1,3:1}, 0 directed cycles, SCCs 1+1+1+1, "
            "and 1 Hamiltonian path"
        )

    print()
    print(
        "Independent checks: signed pair-cycle Burnside gives classes(8)=6880 "
        "and SC(8)=176; exact backtracking agrees with brute permutations through n=6."
    )
    print(f"all assertions passed ({perf_counter() - started:.3f}s)")


if __name__ == "__main__":
    main()
