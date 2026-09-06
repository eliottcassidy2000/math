#!/usr/bin/env python3
"""Independent exact audit for THM-4448 general shore attachment.

This program is self-contained: it imports no primary THM-4448 code.  Every
interval endpoint and every decision is integral or Fraction-exact.
"""

from __future__ import annotations

from fractions import Fraction as Q
from math import gcd


THRESHOLD = Q(1, 14)
EXPECTED_FAILURE = {
    (1, 5, 11): (
        (Q(25, 154), Q(31, 154), (0, 1, 2)),
        (Q(123, 154), Q(129, 154), (2, 1, 0)),
    ),
    (2, 11, 20): (
        (Q(5, 56), Q(3, 28), (0, 1, 2)),
        (Q(123, 280), Q(129, 280), (1, 2, 0)),
        (Q(151, 280), Q(157, 280), (1, 0, 2)),
        (Q(25, 28), Q(51, 56), (2, 1, 0)),
    ),
}


def dist(x: Q) -> Q:
    x %= 1
    return min(x, 1 - x)


def killed_sheets(t: int, y: Q) -> tuple[int, ...]:
    """Literal fibre mask after representing y in [0,1)."""
    return tuple(j for j in range(3) if dist(Q(t) * (y + j) / 3) < THRESHOLD)


def tail_failure_cells(tails: tuple[int, int, int]):
    walls = {Q(0), Q(1)}
    for t in tails:
        for n in range(t + 1):
            for sign in (-1, 1):
                y = Q(n, t) + sign * Q(3, 14 * t)
                if 0 <= y <= 1:
                    walls.add(y)
    walls = sorted(walls)
    answer = []
    for lo, hi in zip(walls, walls[1:]):
        y = (lo + hi) / 2
        masks = [killed_sheets(t, y) for t in tails]
        if set().union(*map(set, masks)) == {0, 1, 2}:
            assert all(len(mask) == 1 for mask in masks)
            answer.append((lo, hi, tuple(mask[0] for mask in masks)))
    return tuple(answer)


def integer_teeth(speeds: tuple[int, ...]):
    """Strict danger teeth on a common exact integer grid."""
    product = 1
    for v in speeds:
        product *= v // gcd(product, v)
    denominator = 14 * product
    intervals = []
    for v in speeds:
        scale = product // v
        for k in range(v + 1):
            lo = max(0, scale * (14 * k - 1))
            hi = min(denominator, scale * (14 * k + 1))
            if lo < hi:
                intervals.append((lo, hi))
    return denominator, sorted(intervals)


def merged_blocks(speeds: tuple[int, ...], merge_touch: bool = False):
    denominator, intervals = integer_teeth(speeds)
    blocks: list[list[int]] = []
    for lo, hi in intervals:
        if not blocks:
            blocks.append([lo, hi])
            continue
        separated = lo > blocks[-1][1] if merge_touch else lo >= blocks[-1][1]
        if separated:
            blocks.append([lo, hi])
        elif hi > blocks[-1][1]:
            blocks[-1][1] = hi
    return denominator, blocks


def longest_circle_danger(speeds: tuple[int, ...], merge_touch: bool = False) -> Q:
    denominator, blocks = merged_blocks(speeds, merge_touch)
    lengths = [hi - lo for lo, hi in blocks]
    if len(blocks) > 1 and blocks[0][0] == 0 and blocks[-1][1] == denominator:
        lengths = [blocks[0][1] + denominator - blocks[-1][0]] + [
            hi - lo for lo, hi in blocks[1:-1]
        ]
    return Q(max(lengths), denominator)


def safe_components(speeds: tuple[int, ...]):
    denominator, blocks = merged_blocks(speeds)
    answer = []
    cursor = 0
    for lo, hi in blocks:
        if cursor <= lo and cursor != 0:
            answer.append((Q(cursor, denominator), Q(lo, denominator)))
        cursor = max(cursor, hi)
    if cursor < denominator:
        answer.append((Q(cursor, denominator), Q(1)))
    return tuple(answer)


def intersection_measure(left, right) -> Q:
    total = Q(0)
    for a, b in left:
        for c, d, *_ in right:
            total += max(Q(0), min(b, d) - max(a, c))
    return total


def admissible_sum(n: int) -> bool:
    p = 2
    while p * p <= n:
        if n % p:
            p += 1
            continue
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if p % 3 != 2 or exponent > 2:
            return False
        p += 1
    return n == 1 or n % 3 == 2


def pair_atlases():
    all_rows = []
    filtered_rows = []
    for total in range(3, 357):
        filtered = admissible_sum(total)
        for p in range(1, (total + 1) // 2):
            q = total - p
            if gcd(p, q) != 1:
                continue
            row = (longest_circle_danger((p, q)), p, q)
            all_rows.append(row)
            if filtered:
                filtered_rows.append(row)
    return all_rows, filtered_rows


def cross_height(body: tuple[int, ...], pair: tuple[int, int]) -> int:
    shore = set(pair)
    return min(
        max(x // gcd(x, z), z // gcd(x, z))
        for x in pair
        for z in body
        if z not in shore
    )


def hostile_checks():
    small = tuple(range(1, 9))
    controls = (
        ((1, 5, 11), 53, 2310, Q(2, 11), Q(3, 154)),
        ((2, 11, 20), 121, 210, Q(1, 10), Q(1, 140)),
    )
    rows = []
    for tails, base, step, y, boundary_distance in controls:
        # These congruences prove gcd(N,1*...*8)=1 and freeze N*y mod one
        # for every k, not only for the diagnostic samples below.
        assert step % 210 == 0 and gcd(base, 210) == 1
        assert (step * y.numerator) % y.denominator == 0
        assert all(gcd(base, a) == 1 for a in small)
        for k in (0, 1, 2, 17, 199):
            n = base + step * k
            body = small + (n, 4 * n)
            assert all(gcd(n, a) == 1 for a in small)
            assert cross_height(body, (n, 4 * n)) == n
            assert all(dist(Q(c) * y) > THRESHOLD for c in body)
            assert set().union(*(set(killed_sheets(t, y)) for t in tails)) == {0, 1, 2}
            assert Q(6, 7 * n) < boundary_distance
            assert 14 * n >= 87 * max(small)
            assert 14 * n >= 29 * max(tails)
        rows.append((tails, base, step, y, boundary_distance))
    return tuple(rows)


def main() -> None:
    print("LRC14_GENERAL_SHORE_ATTACHMENT_THM4448_INDEPENDENT")
    print("STATUS=FINITE_EXACT_REPLAY_PLUS_EXACT_STRUCTURAL_CONTROLS;LRC14_OPEN")

    for tails, expected in EXPECTED_FAILURE.items():
        actual = tail_failure_cells(tails)
        assert actual == expected
        mass = sum((hi - lo for lo, hi, _ in actual), Q(0))
        print(f"tail={tails} cells={actual} mass={mass}")

    # Equality walls are safe.  Closed merging erases the two singleton good
    # components of (1,13) and creates the stated false 15/91 danger block.
    singletons = tuple(j for j in safe_components((1, 13)) if j[0] == j[1])
    assert singletons == ((Q(1, 14), Q(1, 14)), (Q(13, 14), Q(13, 14)))
    strict_gap = longest_circle_danger((1, 13))
    closed_merge_gap = longest_circle_danger((1, 13), merge_touch=True)
    assert strict_gap == Q(1, 7) and closed_merge_gap == Q(15, 91)
    print(
        f"strict_touch_control=(1,13) singleton_good={singletons} "
        f"strict_gap={strict_gap} closed_merge_spurious={closed_merge_gap}"
    )

    equality_controls = {
        (1, 5, 11): (1, 2, 3, 4, 7, 8, 9, 13, 14, 19),
        (2, 11, 20): (1, 2, 3, 4, 5, 6, 7, 8, 12, 14),
    }
    for tails, body in equality_controls.items():
        failure = EXPECTED_FAILURE[tails]
        target = sum((hi - lo for lo, hi, _ in failure), Q(0))
        good = safe_components(body)
        overlap = intersection_measure(good, failure)
        assert overlap == target
        assert sum((hi - lo for lo, hi in good), Q(0)) > target
        print(f"full_failure_overlap_control=T{tails} body={body} overlap={overlap}")

    all_rows, filtered_rows = pair_atlases()
    all_max = max(row[0] for row in all_rows)
    filtered_max = max(row[0] for row in filtered_rows)
    all_leaders = [(p, q) for gap, p, q in all_rows if gap == all_max]
    filtered_leaders = [(p, q) for gap, p, q in filtered_rows if gap == filtered_max]
    assert len(all_rows) == 19314
    assert len(filtered_rows) == 5855
    assert all_max == Q(15, 98) and all_leaders == [(1, 14)]
    assert filtered_max == Q(29, 196) and filtered_leaders == [(1, 28)]
    print(f"unfiltered_pairs={len(all_rows)} max={all_max} leaders={all_leaders}")
    print(
        f"decoder_filtered_pairs={len(filtered_rows)} max={filtered_max} "
        f"leaders={filtered_leaders}"
    )
    print("decoder_seam_scope=46837_positive_scale_triples_plus_5855_scale_one_equals_52692")

    # Exact finite controls for the pullback scaling and the LRC margin used
    # by the general r-shore proof.
    scale_tests = 0
    for r in range(1, 10):
        u = tuple(range(1, r + 1))
        delta = longest_circle_danger(u)
        margin = Q(r, 14 * (14 - r))
        assert margin == Q(1, 14 - r) - Q(1, 14) and margin > 0
        for h in range(1, 6):
            assert longest_circle_danger(tuple(h * v for v in u)) == delta / h
            scale_tests += 1
    print(f"general_r_margin_and_pullback_scale_tests={scale_tests} r=1..9 h=1..5")

    # A divisibility/overlap control.  For t=3s the mask is all three sheets
    # or empty, while the selected local inverse branch remains meaningful.
    assert killed_sheets(3, Q(0)) == (0, 1, 2)
    assert killed_sheets(3, Q(1, 2)) == ()
    body = tuple(range(1, 9))
    tails = (3, 6, 9)
    x0 = Q(1, 27)
    y0 = 3 * x0
    rho = min(Q(1, 84 * max(body)), Q(1, 28 * max(tails)))
    base = tuple(sorted(set(3 * v for v in body) | set(tails)))
    assert len(base) == 8 <= 11
    assert min(dist(Q(v) * x0) for v in base) >= Q(1, 12)
    for y in (y0 - rho, y0, y0 + rho):
        assert all(dist(Q(v) * y) >= THRESHOLD for v in body)
        assert all(dist(Q(t) * y / 3) >= THRESHOLD for t in tails)
    print(
        f"nonunit_tail_control=T{tails} distinct_base_speeds={len(base)} "
        f"rho={rho} local_branch_endpoints_safe=True"
    )

    hostile = hostile_checks()
    print(f"hostile_progressions={hostile} sampled_k=(0,1,2,17,199) exact=True")
    print("CONCLUSION=INDEPENDENT_CHECKS_PASS;GENERAL_ATTACHMENT_TYPED;LRC14_OPEN")


if __name__ == "__main__":
    main()
