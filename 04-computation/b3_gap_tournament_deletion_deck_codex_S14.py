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
* the complete radius-one node/line/tiling classification;
* the invariant-ring construction of the even and odd half recursions;
* the sign correction in the anchored c3 deletion sieve;
* Tournament Analysis on radius-one merged nodes.

Tournament Analysis uses radius-one merged node addresses (gap, mirror depth)
as vertices.  The pair observable is lexicographic comparison, with either
cyclic depth or mirror depth chosen as the first gauge.  Ties use the other
coordinate, producing the displayed Hamiltonian paths.  These quotients
preserve one-flip C3, H, converse sheet, and line colour; they destroy all
multi-flip packet interactions and therefore do not preserve LRC loneliness.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
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
    return +result


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
    flips = sum(
        (position_gap[u] < position_gap[v]) != (position_depth[u] < position_depth[v])
        for u, v in itertools.combinations(vertices, 2)
    )
    size = len(vertices)
    return {
        "vertices": size,
        "score_histogram": {str(score): 1 for score in range(size)},
        "directed_triangles": 0,
        "scc_sizes": [1] * size,
        "hamiltonian_path_count": 1,
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

    checked_masks = 0
    if exhaustive:
        all_mask = (1 << m(n)) - 1
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
        "gap_shell": n14_shell_table(),
    }
    assert n14["m"] == n14["h"] + n14["h_previous"]

    return {
        "method": {
            "alternate_vertices": "radius-one merged node addresses (gap, mirror depth)",
            "pairwise_observable": "lexicographic comparison",
            "gauges": ["gap first", "mirror depth first"],
            "tie_hamiltonian_path": "the corresponding lexicographic order",
            "preserves": ["one-flip C3", "one-flip H", "converse sheet", "blue/black packet"],
            "destroys": ["multi-flip packet interactions", "general node fibre", "LRC metric loneliness"],
            "challenged_assumption": "tournament vertices need not be original runners or arcs; path gaps and information carriers work",
        },
        "size_reports": reports,
        "recurrence": recurrence,
        "anchored_c3_sieve": verify_anchored_sieve(),
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
