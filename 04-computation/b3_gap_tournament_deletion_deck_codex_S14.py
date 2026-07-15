#!/usr/bin/env python3
"""Exact gap-tournament/deletion-deck audit for THM-830.

The vertices used here are deliberately not runners or original tournament
vertices.  A free chord (a,b) of a fixed Hamiltonian path is sent to the arc
(a-1,b) between two path gaps.  This identifies the full fixed-path tiling
cube X_n with the full labelled tournament cube on n-1 ordered gaps.

The audit verifies:

* positive-composition and simplex forms of the A/B/C seven-term recursion;
* the gap-tournament bijection, its seam row, and actual deletion cards;
* the B face as the Cech consensus of internal repaired cards;
* weighted deletion roles, incidence/Gram/rank laws, and mirror current;
* the half-cube fibre product and exact blue/black defect normal form;
* the exact defect relation algebra and its Hamming/Krawtchouk fusion;
* the complete radius-one node/line/tiling classification;
* the invariant-ring construction of the even and odd half recursions;
* the sign correction in the anchored c3 deletion sieve;
* Tournament Analysis on radius-one merged nodes.

Tournament Analysis uses radius-one merged node addresses (gap, mirror depth)
as vertices.  The pair observable is lexicographic comparison, with either
cyclic depth or mirror depth chosen as the first gauge.  Ties use the other
coordinate, producing the displayed Hamiltonian paths.  These quotients
preserve one-flip C3, H, absolute mirror depth, and line colour; they destroy
the signed converse sheet and all multi-flip packet interactions, and
therefore do not preserve LRC loneliness.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache


@dataclass(frozen=True)
class Root:
    a: int
    b: int
    u: int
    g: int
    w: int
    d: int
    x: int
    y: int
    z: int


@lru_cache(maxsize=None)
def roots(n: int) -> tuple[Root, ...]:
    result = []
    for a in range(3, n + 1):
        for b in range(1, a - 1):
            u = b
            g = a - b - 1
            w = n - a + 1
            result.append(
                Root(
                    a=a,
                    b=b,
                    u=u,
                    g=g,
                    w=w,
                    d=w - u,
                    x=w - 1,
                    y=u - 1,
                    z=g - 1,
                )
            )
    return tuple(result)


def m(n: int) -> int:
    return math.comb(n - 1, 2) if n >= 3 else 0


def fixed_cells(n: int) -> int:
    return (n - 1) // 2


def half_cells(n: int) -> int:
    return ((n - 1) ** 2) // 4


def hamming_intersection(r: int, i: int, j: int, k: int) -> int:
    twice_on_support = i - j + k
    twice_off_support = i + j - k
    if twice_on_support % 2 or twice_off_support % 2:
        return 0
    on_support = twice_on_support // 2
    off_support = twice_off_support // 2
    if not (0 <= on_support <= k and 0 <= off_support <= r - k):
        return 0
    return math.comb(k, on_support) * math.comb(r - k, off_support)


def krawtchouk(r: int, distance: int, weight: int) -> int:
    return sum(
        (-1) ** chosen_on_support
        * math.comb(weight, chosen_on_support)
        * math.comb(r - weight, distance - chosen_on_support)
        for chosen_on_support in range(distance + 1)
        if chosen_on_support <= weight and distance - chosen_on_support <= r - weight
    )


def reflection(root: Root, n: int) -> tuple[int, int]:
    return n + 1 - root.b, n + 1 - root.a


def face_word(root: Root) -> str:
    return "".join(
        name
        for name, coordinate in (("A", root.x), ("B", root.z), ("C", root.y))
        if coordinate > 0
    )


@lru_cache(maxsize=None)
def root_index(n: int) -> dict[tuple[int, int], int]:
    return {(root.a, root.b): index for index, root in enumerate(roots(n))}


def reflect_mask(mask: int, n: int) -> int:
    index = root_index(n)
    value = 0
    for bit, root in enumerate(roots(n)):
        value |= ((mask >> bit) & 1) << index[reflection(root, n)]
    return value


def complement_mask(mask: int, n: int) -> int:
    return mask ^ ((1 << m(n)) - 1)


@lru_cache(maxsize=None)
def gap_pairs(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((i, j) for i in range(2, n) for j in range(1, i))


def gap_mask(mask: int, n: int) -> int:
    """The full labelled tournament on the n-1 ordered path gaps."""
    pairs = gap_pairs(n)
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    value = 0
    for bit, root in enumerate(roots(n)):
        value |= ((mask >> bit) & 1) << pair_index[(root.a - 1, root.b)]
    return value


def b_face_and_seams(mask: int, n: int) -> tuple[int, tuple[int, ...]]:
    lower_index = root_index(n - 1)
    lower = 0
    seams = []
    upper_index = root_index(n)
    for root in roots(n - 1):
        upper = (root.a + 1, root.b)
        lower |= ((mask >> upper_index[upper]) & 1) << lower_index[(root.a, root.b)]
    for j in range(2, n):
        seams.append((mask >> upper_index[(j + 1, j - 1)]) & 1)
    return lower, tuple(seams)


def card(mask: int, n: int, deleted: int) -> tuple[int, int | None]:
    """Delete one ordered path vertex and retain the shortcut seam separately."""
    lower_index = root_index(n - 1)
    upper_index = root_index(n)
    value = 0
    for lower in roots(n - 1):
        a = lower.a + (lower.a >= deleted)
        b = lower.b + (lower.b >= deleted)
        value |= ((mask >> upper_index[(a, b)]) & 1) << lower_index[(lower.a, lower.b)]
    seam = None
    if 1 < deleted < n:
        seam = (mask >> upper_index[(deleted + 1, deleted - 1)]) & 1
    return value, seam


def d_b(mask: int, n: int) -> int:
    return b_face_and_seams(mask, n)[0]


def u_window(n: int, deleted: int) -> set[tuple[int, int]]:
    return {
        (root.a, root.b)
        for root in roots(n - 1)
        if root.b < deleted <= root.a
    }


def composite_face_support(n: int, alpha: int, beta: int, gamma: int) -> set[tuple[int, int]]:
    """Source gap-tournament arcs retained by A^alpha B^beta C^gamma."""
    size = n - 1 - alpha - beta - gamma
    return {
        (p + beta + gamma, q + gamma)
        for p in range(2, size + 1)
        for q in range(1, p)
    }


@lru_cache(maxsize=None)
def half_keys(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((root.g, root.d) for root in roots(n) if root.d >= 0))


@lru_cache(maxsize=None)
def positive_keys(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((root.g, root.d) for root in roots(n) if root.d > 0))


def normal_form(mask: int, n: int) -> tuple[int, int]:
    by_key = {(root.g, root.d): index for index, root in enumerate(roots(n))}
    x_word = 0
    defect = 0
    for bit, (g, d) in enumerate(half_keys(n)):
        x_word |= ((mask >> by_key[(g, d)]) & 1) << bit
    for bit, (g, d) in enumerate(positive_keys(n)):
        pair_xor = ((mask >> by_key[(g, d)]) ^ (mask >> by_key[(g, -d)])) & 1
        defect |= pair_xor << bit
    return x_word, defect


def from_normal_form(x_word: int, defect: int, n: int) -> int:
    by_key = {(root.g, root.d): index for index, root in enumerate(roots(n))}
    value = 0
    for bit, (g, d) in enumerate(half_keys(n)):
        x_bit = (x_word >> bit) & 1
        value |= x_bit << by_key[(g, d)]
        if d > 0:
            defect_bit = (defect >> positive_keys(n).index((g, d))) & 1
            value |= (x_bit ^ defect_bit) << by_key[(g, -d)]
    return value


def gf2_rank(rows: list[int]) -> int:
    rows = rows[:]
    rank = 0
    while rows:
        pivot = max(rows)
        if pivot == 0:
            break
        bit = pivot.bit_length() - 1
        rows.remove(pivot)
        rows = [row ^ pivot if (row >> bit) & 1 else row for row in rows]
        rank += 1
    return rank


def deletion_incidence_rows(n: int) -> list[int]:
    result = []
    for deleted in range(1, n + 1):
        row = 0
        for bit, root in enumerate(roots(n)):
            if deleted not in (root.a, root.b):
                row |= 1 << bit
        result.append(row)
    return result


def arc(mask: int, n: int, i: int, j: int) -> bool:
    """Return whether i -> j in the fixed-path tournament."""
    if i == j:
        return False
    if i < j:
        return not arc(mask, n, j, i)
    if i == j + 1:
        return True
    bit = root_index(n)[(i, j)]
    return ((mask >> bit) & 1) == 0


def c3(mask: int, n: int) -> int:
    total = 0
    for i, j, k in itertools.combinations(range(1, n + 1), 3):
        total += (
            arc(mask, n, i, j)
            and arc(mask, n, j, k)
            and arc(mask, n, k, i)
        ) or (
            arc(mask, n, i, k)
            and arc(mask, n, k, j)
            and arc(mask, n, j, i)
        )
    return total


def hp_count(mask: int, n: int) -> int:
    dp = [[0] * n for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
    for subset in range(1 << n):
        for last in range(n):
            count = dp[subset][last]
            if not count:
                continue
            for nxt in range(n):
                if (subset >> nxt) & 1:
                    continue
                if arc(mask, n, last + 1, nxt + 1):
                    dp[subset | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def score_signature(mask: int, n: int) -> tuple[int, ...]:
    scores = [sum(arc(mask, n, i, j) for j in range(1, n + 1) if i != j) for i in range(1, n + 1)]
    return tuple(sorted(scores))


def merged_score_signature(mask: int, n: int) -> tuple[int, ...]:
    scores = score_signature(mask, n)
    converse = tuple(sorted(n - 1 - score for score in scores))
    return min(scores, converse)


def deck_vector(mask: int, n: int) -> tuple[int, int, int, int]:
    values = [0, 0, 0, 0]
    for bit, root in enumerate(roots(n)):
        if not ((mask >> bit) & 1):
            continue
        values[0] += root.x
        values[2] += root.y
        if root.g >= 2:
            values[1] += root.g
        else:
            values[3] += 1
    return tuple(values)  # type: ignore[return-value]


def mask_from_tiles(n: int, tiles: tuple[tuple[int, int], ...]) -> int:
    index = root_index(n)
    return sum(1 << index[tile] for tile in tiles)


def packet_eta(left: Root, right: Root) -> int:
    if left.a == right.a or left.b == right.b:
        return -1
    if left.a == right.b or right.a == left.b:
        return 1
    return 0


def c3_quadratic(mask: int, n: int) -> int:
    cells = roots(n)
    value = sum(root.g * ((mask >> bit) & 1) for bit, root in enumerate(cells))
    for left in range(len(cells)):
        if not ((mask >> left) & 1):
            continue
        for right in range(left):
            if (mask >> right) & 1:
                value += packet_eta(cells[left], cells[right])
    return value


def arbitrary_tournament_c3(mask: int, n: int, vertices: tuple[int, ...] | None = None) -> int:
    if vertices is None:
        vertices = tuple(range(n))
    pairs = list(itertools.combinations(range(n), 2))
    index = {pair: bit for bit, pair in enumerate(pairs)}

    def beats(i: int, j: int) -> bool:
        if i < j:
            return ((mask >> index[(i, j)]) & 1) == 0
        return not beats(j, i)

    total = 0
    for i, j, k in itertools.combinations(vertices, 3):
        total += (beats(i, j) and beats(j, k) and beats(k, i)) or (
            beats(i, k) and beats(k, j) and beats(j, i)
        )
    return total


def anchored_sieve_c3(mask: int, n: int, anchor: int = 0) -> int:
    deletable = tuple(vertex for vertex in range(n) if vertex != anchor)
    total = 0
    for size in range(1, len(deletable) + 1):
        sign = 1 if size % 2 else -1
        for deleted in itertools.combinations(deletable, size):
            survivors = tuple(vertex for vertex in range(n) if vertex not in deleted)
            total += sign * arbitrary_tournament_c3(mask, n, survivors)
    return total


def shift(counter: Counter[tuple[int, int]], dg: int, dd: int, factor: int = 1) -> Counter[tuple[int, int]]:
    return Counter({(g + dg, d + dd): factor * value for (g, d), value in counter.items()})


def add_counters(*parts: Counter[tuple[int, int]]) -> Counter[tuple[int, int]]:
    result: Counter[tuple[int, int]] = Counter()
    for part in parts:
        result.update(part)
    return Counter({key: value for key, value in result.items() if value != 0})


def root_spectrum(n: int, half: bool = False) -> Counter[tuple[int, int]]:
    return Counter((root.g, root.d) for root in roots(n) if not half or root.d >= 0)


def invariant_monomials(degree: int) -> set[tuple[int, int, int]]:
    """Monomials s^i Z^j p^k with degrees 1,1,2."""
    return {
        (i, j, k)
        for k in range(degree // 2 + 1)
        for i in range(degree - 2 * k + 1)
        for j in (degree - 2 * k - i,)
    }


def ta_fingerprint(n: int) -> dict[str, object]:
    vertices = sorted((root.g, root.d) for root in roots(n) if root.d >= 0)
    by_gap = sorted(vertices, key=lambda item: (item[0], item[1]))
    by_depth = sorted(vertices, key=lambda item: (item[1], item[0]))
    position_gap = {vertex: index for index, vertex in enumerate(by_gap)}
    position_depth = {vertex: index for index, vertex in enumerate(by_depth)}

    def order_fingerprint(order: list[tuple[int, int]]) -> dict[str, object]:
        position = {vertex: index for index, vertex in enumerate(order)}

        def beats(left: tuple[int, int], right: tuple[int, int]) -> bool:
            return position[left] < position[right]

        scores = Counter(
            sum(beats(vertex, other) for other in order if other != vertex)
            for vertex in order
        )
        triangles = sum(
            (beats(a, b) and beats(b, c) and beats(c, a))
            or (beats(a, c) and beats(c, b) and beats(b, a))
            for a, b, c in itertools.combinations(order, 3)
        )
        reach = [[beats(left, right) for right in order] for left in order]
        for index in range(len(order)):
            reach[index][index] = True
        for pivot in range(len(order)):
            for left in range(len(order)):
                if reach[left][pivot]:
                    reach[left] = [
                        old or reach[pivot][right]
                        for right, old in enumerate(reach[left])
                    ]
        unseen = set(range(len(order)))
        scc_sizes = []
        while unseen:
            seed = next(iter(unseen))
            component = {other for other in unseen if reach[seed][other] and reach[other][seed]}
            scc_sizes.append(len(component))
            unseen -= component
        assert all(beats(order[index], order[index + 1]) for index in range(len(order) - 1))
        assert scores == Counter({score: 1 for score in range(len(order))})
        assert triangles == 0
        assert scc_sizes == [1] * len(order)
        return {
            "score_histogram": {str(score): count for score, count in sorted(scores.items())},
            "directed_triangles": triangles,
            "scc_sizes": scc_sizes,
            "hamiltonian_path_count": 1,
        }

    gap_fingerprint = order_fingerprint(by_gap)
    depth_fingerprint = order_fingerprint(by_depth)
    assert gap_fingerprint == depth_fingerprint
    flips = sum(
        (position_gap[u] < position_gap[v]) != (position_depth[u] < position_depth[v])
        for u, v in itertools.combinations(vertices, 2)
    )
    size = len(vertices)
    return {
        "vertices": size,
        **gap_fingerprint,
        "gap_first_path": by_gap,
        "depth_first_path": by_depth,
        "edge_flips": flips,
    }


def verify_size(n: int, exhaustive: bool) -> dict[str, object]:
    cells = roots(n)
    assert len(cells) == m(n)
    assert all(root.u + root.g + root.w == n for root in cells)
    assert all(root.x + root.y + root.z == n - 3 for root in cells)
    assert all(reflection(root, n) in root_index(n) for root in cells)
    assert len(gap_pairs(n)) == m(n)
    assert {(root.a - 1, root.b) for root in cells} == set(gap_pairs(n))

    face_a = composite_face_support(n, 1, 0, 0)
    face_b = composite_face_support(n, 0, 1, 0)
    face_c = composite_face_support(n, 0, 0, 1)
    assert face_a | face_b | face_c == set(gap_pairs(n))
    assert face_a & face_b == composite_face_support(n, 1, 1, 0)
    assert face_a & face_c == composite_face_support(n, 1, 0, 1)
    assert face_b & face_c == composite_face_support(n, 0, 1, 1)
    assert face_a & face_b & face_c == composite_face_support(n, 1, 1, 1)
    for depth in range(0, n - 2):
        support = composite_face_support(n, 0, depth, 0)
        assert len(support) == math.comb(n - 1 - depth, 2)
        assert len(support) + sum(n - 1 - layer for layer in range(1, depth + 1)) == m(n)

    strata = Counter(face_word(root) for root in cells)
    expected_strata = Counter(
        {
            "A": 1,
            "B": 1,
            "C": 1,
            "AB": n - 4,
            "AC": n - 4,
            "BC": n - 4,
            "ABC": math.comb(n - 4, 2) if n >= 6 else 0,
        }
    )
    assert strata == +expected_strata

    total_a = sum(root.x for root in cells)
    total_c = sum(root.y for root in cells)
    total_internal = sum(root.g for root in cells)
    total_b = sum(root.g for root in cells if root.g >= 2)
    total_s = sum(root.g == 1 for root in cells)
    assert total_a == total_c == math.comb(n - 1, 3)
    assert total_internal == math.comb(n, 3)
    assert total_b == math.comb(n, 3) - (n - 2)
    assert total_s == n - 2
    assert sum(root.g * root.g for root in cells) == (n - 2) * (n - 1) ** 2 * n // 12
    assert sum(root.d * root.d for root in cells) == 2 * math.comb(n, 4)
    assert sum(root.g * root.d for root in cells) == 0

    for root in cells:
        ja = root.x
        jc = root.y
        jb = root.g if root.g >= 2 else 0
        js = 1 if root.g == 1 else 0
        assert ja + jb + jc + js == n - 2

    windows = {j: u_window(n, j) for j in range(2, n)}
    for j, window in windows.items():
        assert len(window) == (j - 1) * (n - j) - 1
    assert set().union(*windows.values()) == {(root.a, root.b) for root in roots(n - 1)}
    for j, k in itertools.combinations(range(2, n), 2):
        assert len(windows[j] & windows[k]) == (j - 1) * (n - k)
    assert sum(map(len, windows.values())) == math.comb(n, 3) - (n - 2)

    rows = deletion_incidence_rows(n)
    row_sums = [row.bit_count() for row in rows]
    assert row_sums[0] == row_sums[-1] == m(n - 1)
    assert all(value == m(n - 1) + 1 for value in row_sums[1:-1])
    assert sum(row_sums) == (n - 2) * m(n)
    expected_rank = 1 if n == 3 else (n if n % 2 else n - 1)
    assert gf2_rank(rows) == expected_rank
    for j, k in itertools.combinations(range(n), 2):
        endpoint_j = int(j in (0, n - 1))
        endpoint_k = int(k in (0, n - 1))
        expected = m(n - 2) + 2 - endpoint_j - endpoint_k - int(k == j + 1)
        assert (rows[j] & rows[k]).bit_count() == expected

    assert len(half_keys(n)) == half_cells(n)
    assert len(positive_keys(n)) == half_cells(n - 1)
    assert m(n) == half_cells(n) + half_cells(n - 1)
    assert fixed_cells(n) == half_cells(n) - half_cells(n - 1)

    line_count = 2 ** (m(n) - 1)
    blue_count = 2 ** (half_cells(n) - 1)
    black_count = line_count - blue_count
    lower_line_count = 2 ** (m(n - 1) - 1) if n > 4 else 1
    lower_blue_count = 2 ** (half_cells(n - 1) - 1) if n > 4 else 1
    lower_black_count = lower_line_count - lower_blue_count
    bb = 2 ** fixed_cells(n) * lower_blue_count
    bk = 0
    kb = (2 ** (n - 2) - 2 ** fixed_cells(n)) * lower_blue_count
    kk = 2 ** (n - 2) * lower_black_count
    assert (bb, bk) == (blue_count, 0)
    assert bb + kb + kk == line_count
    assert kb + kk == black_count
    assert sum(
        2 ** ((n - 2 + 1) // 2) * math.comb((n - 2) // 2, defect)
        for defect in range((n - 2) // 2 + 1)
    ) == 2 ** (n - 2)

    spectrum = root_spectrum(n)
    if n >= 6:
        recurrence = add_counters(
            shift(root_spectrum(n - 1), 1, 0),
            shift(root_spectrum(n - 1), 0, 1),
            shift(root_spectrum(n - 1), 0, -1),
            shift(root_spectrum(n - 2), 1, 1, -1),
            shift(root_spectrum(n - 2), 0, 0, -1),
            shift(root_spectrum(n - 2), 1, -1, -1),
            shift(root_spectrum(n - 3), 1, 0),
        )
        assert spectrum == recurrence

    half_spectrum = root_spectrum(n, half=True)
    if n % 2 == 0:
        even_recurrence = add_counters(
            shift(root_spectrum(n - 1, half=True), 1, 0),
            shift(root_spectrum(n - 1, half=True), 0, 1),
            shift(root_spectrum(n - 2, half=True), 1, 1, -1),
        )
        assert half_spectrum == even_recurrence
    elif n >= 7:
        odd_recurrence = add_counters(
            shift(root_spectrum(n - 1, half=True), 1, 0),
            shift(root_spectrum(n - 1, half=True), 0, 1),
            root_spectrum(n - 2, half=True),
            shift(root_spectrum(n - 2, half=True), 1, 1, -1),
            shift(root_spectrum(n - 3, half=True), 1, 0, -1),
            shift(root_spectrum(n - 3, half=True), 0, 1, -1),
            shift(root_spectrum(n - 4, half=True), 1, 1),
        )
        assert half_spectrum == odd_recurrence

    score_cells: dict[tuple[int, ...], list[Root]] = defaultdict(list)
    merged_cells: dict[tuple[int, ...], list[Root]] = defaultdict(list)
    index = root_index(n)
    for root in cells:
        mask = 1 << index[(root.a, root.b)]
        score_cells[score_signature(mask, n)].append(root)
        merged_cells[merged_score_signature(mask, n)].append(root)
        if n <= 9:
            assert c3(mask, n) == root.g
            assert hp_count(mask, n) == 1 + 2**root.g
    assert len(score_cells) == m(n)
    assert len(merged_cells) == half_cells(n)
    assert all(
        {(root.g, abs(root.d)) for root in fibre} == {(fibre[0].g, abs(fibre[0].d))}
        and len(fibre) == (1 if fibre[0].d == 0 else 2)
        for fibre in merged_cells.values()
    )

    packet_counts = Counter()
    for left, root in enumerate(cells):
        row_sum = 0
        for right, other in enumerate(cells):
            if left == right:
                continue
            row_sum += packet_eta(root, other)
            if left > right:
                packet_counts[packet_eta(root, other)] += 1
        boundary_weight = int(root.b == 1) + int(root.a == n)
        assert row_sum == boundary_weight - 2 * root.g
    assert packet_counts[-1] == 2 * math.comb(n - 1, 3)
    assert packet_counts[1] == math.comb(n - 2, 3)

    checked_masks = 0
    if exhaustive:
        all_mask = (1 << m(n)) - 1
        c3_values = []
        line_steps = []
        line_midpoints = []
        for mask in range(1 << m(n)):
            checked_masks += 1
            assert gap_mask(mask, n) == mask
            bmask, seams = b_face_and_seams(mask, n)
            assert d_b(complement_mask(mask, n), n) == complement_mask(bmask, n - 1)
            assert b_face_and_seams(complement_mask(mask, n), n)[1] == tuple(1 - bit for bit in seams)
            free_sum = 0
            seam_sum = 0
            for deleted in range(1, n + 1):
                card_mask, seam = card(mask, n, deleted)
                free_sum += card_mask.bit_count()
                seam_sum += 0 if seam is None else seam
                if seam is not None:
                    assert seam == seams[deleted - 2]
            assert free_sum + seam_sum == (n - 2) * mask.bit_count()

            for deleted, window in windows.items():
                card_mask, _ = card(mask, n, deleted)
                for lower_root in roots(n - 1):
                    if (lower_root.a, lower_root.b) in window:
                        bit = root_index(n - 1)[(lower_root.a, lower_root.b)]
                        assert ((card_mask >> bit) & 1) == ((bmask >> bit) & 1)

            x_word, defect = normal_form(mask, n)
            assert from_normal_form(x_word, defect, n) == mask
            rx, rd = normal_form(reflect_mask(mask, n), n)
            assert rd == defect
            cx, cd = normal_form(mask ^ all_mask, n)
            assert cx == x_word ^ ((1 << half_cells(n)) - 1)
            assert cd == defect

            direct_c3 = c3(mask, n)
            complement_c3 = c3(mask ^ all_mask, n)
            assert c3_quadratic(mask, n) == direct_c3
            boundary_sum = sum(
                ((mask >> bit) & 1) * (int(root.b == 1) + int(root.a == n))
                for bit, root in enumerate(cells)
            )
            assert complement_c3 - direct_c3 == n - 2 - boundary_sum
            c3_values.append(direct_c3)
            if mask < (mask ^ all_mask):
                line_steps.append(complement_c3 - direct_c3)
                line_midpoints.append(Fraction(complement_c3 + direct_c3, 2))

        count = len(c3_values)
        mean = Fraction(sum(c3_values), count)
        variance = sum((Fraction(value) - mean) ** 2 for value in c3_values) / count
        assert variance == Fraction(n**3 - 7 * n**2 + 20 * n - 16, 32)
        step_second_moment = Fraction(sum(step * step for step in line_steps), len(line_steps))
        assert step_second_moment == Fraction(n - 1, 2)
        midpoint_mean = sum(line_midpoints) / len(line_midpoints)
        midpoint_variance = sum((value - midpoint_mean) ** 2 for value in line_midpoints) / len(line_midpoints)
        assert midpoint_variance == Fraction((n - 3) * (n - 2) ** 2, 32)

    return {
        "n": n,
        "tiles": m(n),
        "half_cells": half_cells(n),
        "fixed_cells": fixed_cells(n),
        "strata": dict(strata),
        "seams": total_s,
        "role_totals": {"A": total_a, "B": total_b, "C": total_c, "S": total_s},
        "deletion_rank_f2": expected_rank,
        "radius_one_merged_nodes": len(merged_cells),
        "radius_one_blue_lines": sum(root.d == 0 for root in cells),
        "radius_one_black_lines": sum(root.d != 0 for root in cells),
        "checked_tilings": checked_masks,
    }


def verify_defect_relation_algebra(n: int) -> None:
    r = half_cells(n - 1)
    components = 2 ** (fixed_cells(n) - 1)
    blue = 2 ** (half_cells(n) - 1)
    assert blue == components * 2**r

    intersection = [
        [
            [hamming_intersection(r, i, j, k) for k in range(r + 1)]
            for j in range(r + 1)
        ]
        for i in range(r + 1)
    ]
    eigenvalues = [
        [krawtchouk(r, i, weight) for weight in range(r + 1)]
        for i in range(r + 1)
    ]

    for i in range(r + 1):
        valency_i = math.comb(r, i)
        for j in range(r + 1):
            valency_j = math.comb(r, j)
            assert sum(
                intersection[i][j][k] * math.comb(r, k) for k in range(r + 1)
            ) == valency_i * valency_j

    # The distance-one relation generates the Hamming Bose-Mesner algebra.
    for j in range(r + 1):
        for weight in range(r + 1):
            assert sum(
                intersection[1][j][k] * eigenvalues[k][weight]
                for k in range(r + 1)
            ) == eigenvalues[1][weight] * eigenvalues[j][weight]

    for i in range(r + 1):
        for j in range(r + 1):
            inner_product = sum(
                math.comb(r, weight) * eigenvalues[i][weight] * eigenvalues[j][weight]
                for weight in range(r + 1)
            )
            assert inner_product == (2**r * math.comb(r, i) if i == j else 0)

    black_square_blue_coefficient = 2**r - 1
    black_square_black_coefficient = 2**r - 2
    if r <= 8:
        assert black_square_blue_coefficient == sum(
            1 for delta in range(1, 2**r) if delta ^ delta == 0
        )
        for result in range(1, 2**r):
            assert black_square_black_coefficient == sum(
                1
                for left in range(1, 2**r)
                if (left ^ result) != 0
            )

    q_packets = blue + sum(
        blue // 2 * math.comb(r, weight) for weight in range(1, r + 1)
    )
    assert q_packets == 2 ** (m(n) - 2) + 2 ** (half_cells(n) - 2)


def verify_half_ring(max_n: int) -> None:
    for n in range(4, max_n + 1):
        degree = n - 3
        universe = invariant_monomials(degree)
        assert len(universe) == half_cells(n)
        s_corner = {monomial for monomial in universe if monomial[0] > 0}
        z_corner = {monomial for monomial in universe if monomial[1] > 0}
        if n % 2 == 0:
            assert universe == s_corner | z_corner
            assert len(s_corner) == len(z_corner) == half_cells(n - 1)
            assert len(s_corner & z_corner) == half_cells(n - 2)
        else:
            p_corner = {monomial for monomial in universe if monomial[2] > 0}
            assert universe == s_corner | z_corner | p_corner
            assert len(p_corner) == half_cells(n - 2)
            assert len(s_corner & z_corner) == half_cells(n - 2)
            assert len(s_corner & p_corner) == len(z_corner & p_corner) == half_cells(n - 3)
            assert len(s_corner & z_corner & p_corner) == half_cells(n - 4)
            assert (0, 0, degree // 2) in p_corner - s_corner - z_corner


def verify_anchored_sieve() -> dict[str, object]:
    result = {}
    for n in (4, 5):
        histogram = Counter()
        total_masks = 1 << math.comb(n, 2)
        deletable = tuple(range(1, n))
        for mask in range(total_masks):
            direct = arbitrary_tournament_c3(mask, n)
            sieve = anchored_sieve_c3(mask, n)
            residual = arbitrary_tournament_c3(mask, n, deletable) if len(deletable) == 3 else 0
            assert sieve == direct - residual
            histogram[direct - sieve] += 1
        result[str(n)] = {
            "tournaments": total_masks,
            "direct_minus_sieve": dict(histogram),
        }
    return result


def verify_preservation_boundaries() -> dict[str, object]:
    examples = [
        (
            5,
            ((4, 2), (5, 2)),
            ((4, 1), (5, 3)),
            (1, 2, 2, 1),
            (2, 3),
        ),
        (
            6,
            ((6, 2), (5, 3)),
            ((6, 4), (5, 1)),
            (1, 3, 3, 1),
            (4, 4),
        ),
    ]
    result = []
    for n, left_tiles, right_tiles, expected_deck, expected_c3 in examples:
        left = mask_from_tiles(n, left_tiles)
        right = mask_from_tiles(n, right_tiles)
        assert deck_vector(left, n) == deck_vector(right, n) == expected_deck
        assert (c3(left, n), c3(right, n)) == expected_c3
        result.append(
            {
                "n": n,
                "left_tiles": left_tiles,
                "right_tiles": right_tiles,
                "deck": expected_deck,
                "c3": expected_c3,
                "score_signatures": (score_signature(left, n), score_signature(right, n)),
            }
        )
    assert result[1]["score_signatures"][0] != result[1]["score_signatures"][1]
    return {"examples": result}


def n14_shell_table() -> list[dict[str, int]]:
    n = 14
    rows = []
    for g in range(1, n - 1):
        literal = n - g - 1
        blue = int((n - g) % 2 == 0)
        merged = (n - g) // 2
        rows.append(
            {
                "gap": g,
                "C3": g,
                "H": 1 + 2**g,
                "literal_lines": literal,
                "blue_lines": blue,
                "black_lines": literal - blue,
                "merged_nodes_or_Q_packets": merged,
                "black_packets": merged - blue,
            }
        )
    return rows


def run(max_n: int, exhaustive_n: int) -> dict[str, object]:
    assert max_n >= 6
    assert 4 <= exhaustive_n <= 7
    reports = [verify_size(n, n <= exhaustive_n) for n in range(4, max_n + 1)]
    for n in range(4, max_n + 1):
        verify_defect_relation_algebra(n)
    verify_half_ring(max_n)

    recurrence = []
    for n in range(4, max_n + 1):
        hn = half_cells(n)
        curvature = hn - 2 * half_cells(n - 1) + half_cells(n - 2)
        assert curvature == (n % 2)
        if n % 2 == 0:
            assert hn == 2 * half_cells(n - 1) - half_cells(n - 2)
        if n >= 5:
            assert hn == 2 * half_cells(n - 1) - 2 * half_cells(n - 3) + half_cells(n - 4)
        recurrence.append(
            {
                "n": n,
                "m": m(n),
                "h": hn,
                "f": fixed_cells(n),
                "delta2_h": curvature,
            }
        )

    ta = {str(n): ta_fingerprint(n) for n in range(4, max_n + 1)}
    n14 = {
        "m": m(14),
        "h": half_cells(14),
        "h_previous": half_cells(13),
        "f": fixed_cells(14),
        "literal_lines": 2 ** (m(14) - 1),
        "blue_lines": 2 ** (half_cells(14) - 1),
        "black_lines": (2 ** half_cells(13) - 1) * 2 ** (half_cells(14) - 1),
        "reflection_orbit_lines": 2 ** (m(14) - 2) + 2 ** (half_cells(14) - 2),
        "blue_fraction": "2^-36",
        "defect_slice_formula": "2^41 * binom(36,d)",
        "boundary_components": 2 ** (fixed_cells(14) - 1),
        "objects_per_component": 2 ** half_cells(13),
        "exact_black_Q_packets": 2 ** (half_cells(14) - 2),
        "weight_Q_formula": "2^40 * binom(36,d), d>0",
        "gap_shell": n14_shell_table(),
    }
    assert n14["m"] == n14["h"] + n14["h_previous"]

    return {
        "method": {
            "alternate_vertices": "radius-one merged node addresses (gap, mirror depth)",
            "pairwise_observable": "lexicographic comparison",
            "gauges": ["gap first", "mirror depth first"],
            "tie_hamiltonian_path": "the corresponding lexicographic order",
            "preserves": ["one-flip C3", "one-flip H", "absolute mirror depth", "blue/black packet"],
            "destroys": ["signed converse sheet", "multi-flip packet interactions", "general node fibre", "LRC metric loneliness"],
            "challenged_assumption": "tournament vertices need not be original runners or arcs; path gaps and information carriers work",
        },
        "size_reports": reports,
        "recurrence": recurrence,
        "anchored_c3_sieve": verify_anchored_sieve(),
        "preservation_boundaries": verify_preservation_boundaries(),
        "tournament_analysis": ta,
        "n14": n14,
    }


def print_report(report: dict[str, object]) -> None:
    print("THM-830 GAP-TOURNAMENT / DELETION-DECK EXACT AUDIT")
    print("=" * 68)
    print("method:", json.dumps(report["method"], sort_keys=True))
    print("\nsize ledger")
    for row in report["size_reports"]:  # type: ignore[assignment]
        print(
            "n={n:2d} m={tiles:3d} h={half_cells:3d} f={fixed_cells:2d} "
            "rank(deck F2)={deletion_rank_f2:2d} root-nodes={radius_one_merged_nodes:3d} "
            "root B/K={radius_one_blue_lines:2d}/{radius_one_black_lines:3d} "
            "checked={checked_tilings}".format(**row)
        )
    print("\nrecurrence curvature")
    for row in report["recurrence"]:  # type: ignore[assignment]
        print("n={n:2d} m={m:3d} h={h:3d} f={f:2d} Delta2(h)={delta2_h}".format(**row))
    print("\nanchored c3 sieve (direct - sieve; n=4 must UNDERSHOOT)")
    print(json.dumps(report["anchored_c3_sieve"], sort_keys=True))
    print("\nTournament Analysis edge flips: gap-first vs depth-first")
    ta = report["tournament_analysis"]
    print(" ".join(f"n={n}:{ta[n]['edge_flips']}" for n in sorted(ta, key=int)))  # type: ignore[index]
    print("\nn=14 exact ledger")
    n14 = report["n14"]  # type: ignore[assignment]
    for key in (
        "m",
        "h",
        "h_previous",
        "f",
        "literal_lines",
        "blue_lines",
        "black_lines",
        "reflection_orbit_lines",
        "blue_fraction",
        "defect_slice_formula",
        "boundary_components",
        "objects_per_component",
        "exact_black_Q_packets",
        "weight_Q_formula",
    ):
        print(f"{key}: {n14[key]}")
    print("gap C3 H literal blue black merged blackQ")
    for row in n14["gap_shell"]:
        print(
            "{gap:2d} {C3:3d} {H:5d} {literal_lines:7d} {blue_lines:4d} "
            "{black_lines:5d} {merged_nodes_or_Q_packets:6d} {black_packets:6d}".format(**row)
        )
    print("\nALL ASSERTIONS PASSED")


def main() -> None:
    if not __debug__:
        raise RuntimeError("this exact audit requires Python assertions; do not run with -O")
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=14)
    parser.add_argument("--exhaustive-n", type=int, default=7)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    report = run(args.max_n, args.exhaustive_n)
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print_report(report)


if __name__ == "__main__":
    main()
