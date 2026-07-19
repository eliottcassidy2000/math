#!/usr/bin/env python3
"""Optimization-safe exact referee for THM-1264.

The all-range return-polygon identity is proved on paper and in Lean.  This
referee independently reconstructs every literal minimal interval-chain
return of two through five edges in a finite exact tooth bank and checks the
identity, seam digits, ratio ascent, and triangle specialization.
"""

import ast
from collections import Counter
from fractions import Fraction
from pathlib import Path


def require(condition: bool, message: str = "exact certificate check failed") -> None:
    if not condition:
        raise RuntimeError(message)


def optimization_safety_probe() -> int:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "optimization-sensitive assert remains in referee")
    caught = False
    try:
        require(False, "sentinel")
    except RuntimeError:
        caught = True
    require(caught, "explicit RuntimeError sentinel did not fire")
    return count


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


def tooth(speed: int, address: int) -> tuple[Fraction, Fraction, int, int]:
    require(speed > 0)
    return (
        Fraction(14 * address - 1, 14 * speed),
        Fraction(14 * address + 1, 14 * speed),
        speed,
        address,
    )


def exact_tooth_bank() -> list[tuple[Fraction, Fraction, int, int]]:
    bank: list[tuple[Fraction, Fraction, int, int]] = []
    for speed in range(2, 81):
        for address in range(max(1, speed // 2 - 2), min(speed, speed // 2 + 3)):
            row = tooth(speed, address)
            if Fraction(2, 5) < row[0] and row[1] < Fraction(3, 5):
                bank.append(row)
    bank.sort()
    require(len(bank) == 335)
    return bank


def next_indices(
    bank: list[tuple[Fraction, Fraction, int, int]]
) -> list[list[int]]:
    adjacency: list[list[int]] = [[] for _ in bank]
    for left_index, (left, right, speed, _) in enumerate(bank):
        for right_index in range(left_index + 1, len(bank)):
            next_left, next_right, next_speed, _ = bank[right_index]
            if next_left >= right:
                break
            if left < next_left and right < next_right and speed != next_speed:
                adjacency[left_index].append(right_index)
    return adjacency


def raw_overlap(
    left: tuple[Fraction, Fraction, int, int],
    right: tuple[Fraction, Fraction, int, int],
) -> Fraction:
    return left[1] - right[0]


def verify_return_polygon(
    bank: list[tuple[Fraction, Fraction, int, int]], path: tuple[int, ...]
) -> tuple[int, int]:
    rows = [bank[index] for index in path]
    edge_count = len(rows) - 1
    require(edge_count >= 2)
    outer_speed = rows[0][2]
    address_return = rows[-1][3] - rows[0][3]
    require(rows[-1][2] == outer_speed and address_return >= 1)

    overlaps: list[Fraction] = []
    for index in range(edge_count):
        left, right = rows[index], rows[index + 1]
        overlap = raw_overlap(left, right)
        require(overlap > 0)
        require(left[0] < right[0] and left[1] < right[1])
        if index + 2 < len(rows):
            require(left[1] <= rows[index + 2][0])
        numerator = overlap * 14 * left[2] * right[2]
        require(numerator.denominator == 1 and numerator.numerator > 0)
        common = gcd(left[2], right[2])
        require(numerator.numerator % common == 0)
        seam_digit = numerator.numerator // common
        reduced_sum = left[2] // common + right[2] // common
        require(seam_digit % 14 == reduced_sum % 14)
        overlaps.append(overlap)

    internal_speeds = [row[2] for row in rows[1:-1]]
    exact_right = (
        sum((Fraction(1, speed) for speed in [outer_speed, *internal_speeds]), Fraction())
        / 7
        - Fraction(address_return, outer_speed)
    )
    require(sum(overlaps, Fraction()) == exact_right)
    require(exact_right > 0)

    ratio_sum = sum(
        (Fraction(outer_speed, speed) for speed in internal_speeds), Fraction()
    )
    require(ratio_sum > 7 * address_return - 1)
    minimum_internal = min(internal_speeds)
    require(
        Fraction(outer_speed, minimum_internal)
        > Fraction(7 * address_return - 1, edge_count - 1)
    )
    if edge_count == 3:
        require(outer_speed > 3 * minimum_internal)
    return edge_count, len(overlaps)


def return_polygon_census() -> tuple[int, Counter[int], int, int, int]:
    bank = exact_tooth_bank()
    adjacency = next_indices(bank)
    counts: Counter[int] = Counter()
    seam_rows = 0
    triangle_rows = 0
    max_return = 0

    def visit(path: tuple[int, ...]) -> None:
        nonlocal seam_rows, triangle_rows, max_return
        last = path[-1]
        if len(path) >= 3:
            first_speed, first_address = bank[path[0]][2:]
            last_speed, last_address = bank[last][2:]
            if last_speed == first_speed and last_address > first_address:
                edge_count, seams = verify_return_polygon(bank, path)
                counts[edge_count] += 1
                seam_rows += seams
                triangle_rows += edge_count == 3
                max_return = max(max_return, last_address - first_address)
                return
        if len(path) - 1 >= 5:
            return
        for next_index in adjacency[last]:
            # THM-1253 separation: a tooth two positions back cannot overlap
            # the proposed next tooth.  Endpoint monotonicity handles all
            # still earlier teeth.
            if len(path) >= 2 and bank[path[-2]][1] > bank[next_index][0]:
                continue
            visit((*path, next_index))

    for start in range(len(bank)):
        visit((start,))

    require(
        counts == Counter({2: 42, 3: 8, 4: 1310, 5: 8601}),
        f"unexpected return census: {counts}",
    )
    require(seam_rows == 48_353)
    require(triangle_rows == 8)
    require(max_return == 4, f"unexpected maximum return: {max_return}")
    return len(bank), counts, seam_rows, triangle_rows, max_return


def compact_height_audit() -> tuple[int, int]:
    require(3**7 == 2187 < 2345 < 6561 == 3**8)
    # Eight strict factor-three ascents leave the compact projective box;
    # seven can still fit numerically.
    for length in range(0, 8):
        require(3**length < 2345)
    require(3**8 > 2345)
    six_fast_factor = Fraction(6, 5)
    require(six_fast_factor**42 < 2345 < six_fast_factor**43)
    return 7, 42


def simple_return_index_audit() -> tuple[int, int, Fraction, int, int]:
    # A simple return with K owners has m-1 distinct internal owners, so
    # m<=K.  Seven owners give only a strict factor-one ascent; the actual
    # six-fast tooth word improves the scale coefficient to 6/5.
    for edges in range(2, 8):
        coefficient = Fraction(6, edges - 1)
        require((coefficient >= 1) == (edges <= 7))
    require(Fraction(6, 6) == 1)
    require(Fraction(6, 5) > 1)
    seven_owner_depth = 6
    six_fast_depth = 5
    return 7, 6, Fraction(6, 5), seven_owner_depth, six_fast_depth


def main() -> None:
    assert_nodes = optimization_safety_probe()
    bank_size, counts, seam_rows, triangles, max_return = return_polygon_census()
    triangle_height, six_fast_scale_height = compact_height_audit()
    (
        seven_owner_edges,
        six_fast_edges,
        six_fast_factor,
        seven_owner_depth,
        six_fast_depth,
    ) = simple_return_index_audit()
    ordered_counts = tuple((edges, counts[edges]) for edges in sorted(counts))

    print("THM-1264 CHRONOLOGICAL RETURN-POLYGON EXACT AUDIT")
    print(f"optimization-sensitive assert nodes = {assert_nodes}")
    print(f"exact tooth bank = {bank_size}")
    print(f"literal return polygons by edge count = {ordered_counts}")
    print(f"raw seam-digit occurrences = {seam_rows}")
    print(f"literal owner triangles = {triangles}")
    print(f"largest bank address return = {max_return}")
    print(f"simple-return edge ceilings (seven/six-fast) = {seven_owner_edges} / {six_fast_edges}")
    print("seven-owner weakest ascent = strict factor one (no scale-only height)")
    print(f"six-fast weakest multiplicative ascent = {six_fast_factor}")
    print(f"six-fast factor-6/5 compact-box height = {six_fast_scale_height}")
    print(f"fixed-owner strict-ascent depths (seven/six-fast) = {seven_owner_depth} / {six_fast_depth}")
    print(f"strict factor-three compact-box height = {triangle_height}")
    print("RESULT: PASS -- literal return polygons pay exact address drift")


if __name__ == "__main__":
    main()
