#!/usr/bin/env python3
"""Exact referee for THM-2562: canonical duty-commutator holotopy."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations


P = 13
HALF = 6
BASE = 2_316_060
A = 210_552
B = 12
C = 19_128
N7 = 1_439_676
FIBRE = 27_825_593_350_009
CHECKS = 0

Point = tuple[int, int]


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def add(q: Point, r: Point) -> Point:
    return ((q[0] + r[0]) % P, (q[1] + r[1]) % P)


def scale(k: int, q: Point) -> Point:
    return ((k * q[0]) % P, (k * q[1]) % P)


def duty(q: Point) -> int:
    x, y = q
    return (BASE + A * (int(x == 0) + int(y == 0))
            + B * int((x + y) % P == 0)
            + C * int(x == 0 and y == 0))


def rank_q(matrix: list[list[int | Fraction]]) -> int:
    """Exact Gaussian rank over Q for the tiny matrices used here."""
    if not matrix:
        return 0
    work = [[Fraction(value) for value in row] for row in matrix]
    rows = len(work)
    cols = len(work[0])
    pivot_row = 0
    for col in range(cols):
        pivot = next((i for i in range(pivot_row, rows)
                      if work[i][col]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][col]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for i in range(rows):
            if i == pivot_row or not work[i][col]:
                continue
            factor = work[i][col]
            work[i] = [work[i][j] - factor * work[pivot_row][j]
                       for j in range(cols)]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


H = [[int((i - j) % P <= HALF) for j in range(P)]
     for i in range(P)]


def commutator_block(values: list[int]) -> list[list[int]]:
    return [[H[i][j] * (values[i] - values[j]) for j in range(P)]
            for i in range(P)]


def lines_for_gain(g: Point) -> list[list[Point]]:
    remaining = {(x, y) for x in range(P) for y in range(P)}
    lines: list[list[Point]] = []
    while remaining:
        q = min(remaining)
        line = [add(q, scale(k, g)) for k in range(P)]
        require(len(set(line)) == P, "gain line is not a 13-cycle")
        for point in line:
            remaining.remove(point)
        lines.append(line)
    require(len(lines) == P, "gain did not partition target plane")
    return lines


def punctures(values: list[int]) -> tuple[list[int], list[int]]:
    baseline = min(values)
    support = [i for i, value in enumerate(values) if value != baseline]
    weights = [values[i] - baseline for i in support]
    return support, weights


def factor_block(support: list[int], weights: list[int]) -> list[list[int]]:
    weight_at = {s: weights[i] for i, s in enumerate(support)}
    return [[H[i][j] * (weight_at.get(i, 0) - weight_at.get(j, 0))
             for j in range(P)] for i in range(P)]


def constraint_matrix(support: list[int]) -> list[list[int]]:
    rows: list[list[int]] = []
    for s in support:
        row = [0] * P
        row[s] = 1
        rows.append(row)
    rows.extend(H[s][:] for s in support)
    return rows


def cyclic_gaps(support: tuple[int, int, int]) -> tuple[int, int, int]:
    a, b, c = support
    return tuple(sorted((b - a, c - b, P - c + a)))


def gain_class(g: Point) -> str:
    x, y = g
    if x == 0 or y == 0 or (x + y) % P == 0:
        return "axis_anti"
    slope = y * pow(x, -1, P) % P
    if slope == 1:
        return "diagonal"
    if slope in (6, 11):
        return "quadratic"
    return "generic"


def flat_arc(g: Point, q: Point) -> bool:
    return len({duty(add(q, scale(j, g))) for j in range(HALF + 1)}) == 1


def main() -> None:
    points = [(x, y) for x in range(P) for y in range(P)]
    require(duty((0, 0)) == 2_756_304, "origin duty failed")

    # Independent finite incidence census for every 3-puncture shape.
    gap_types = {
        (1, 1, 11), (1, 2, 10), (1, 3, 9), (1, 4, 8),
        (1, 5, 7), (1, 6, 6), (2, 2, 9), (2, 3, 8),
        (2, 4, 7), (2, 5, 6), (3, 3, 7), (3, 4, 6),
        (3, 5, 5), (4, 4, 5),
    }
    seen_gap_types: set[tuple[int, int, int]] = set()
    exceptional_sets = 0
    for support in combinations(range(P), 3):
        target = [i for i in range(P) if i not in support]
        h_ts = [[H[i][j] for j in support] for i in target]
        h_st = [[H[i][j] for j in target] for i in support]
        gap = cyclic_gaps(support)
        seen_gap_types.add(gap)
        expected = 2 if gap == (1, 6, 6) else 3
        require(rank_q(h_ts) == expected and rank_q(h_st) == expected,
                "cyclic incidence rank classification failed")
        exceptional_sets += int(gap == (1, 6, 6))
    require(seen_gap_types == gap_types, "wrong cyclic gap universe")
    require(exceptional_sets == 13, "wrong exceptional gap-set count")

    equal_values = [0] * P
    unequal_values = [0] * P
    for i, value in {0: A, 6: B, 12: A}.items():
        equal_values[i] = value
    for i, value in {0: A, 6: A, 12: B}.items():
        unequal_values[i] = value
    require(rank_q(commutator_block(equal_values)) == 4,
            "equal Schur boundary is not rank four")
    require(rank_q(commutator_block(unequal_values)) == 5,
            "unequal Schur boundary is not rank five")
    signed_circuit = [0] * P
    signed_circuit[0] = -1
    signed_circuit[1] = -(A // B - 1)
    signed_circuit[12] = 1
    require(all(sum(commutator_block(equal_values)[i][j]
                        * signed_circuit[j] for j in range(P)) == 0
                for i in range(P)), "rank-four signed circuit failed")

    rank_histogram: Counter[int] = Counter()
    nullity_histogram: Counter[int] = Counter()
    flat_histogram: Counter[int] = Counter()
    class_counts: Counter[str] = Counter()
    line_blocks = 0
    anchor_rows = 0
    rank4_blocks = 0
    rank5_blocks = 0

    expected_patterns = {
        "axis_anti": Counter({4: 12, 2: 1}),
        "diagonal": Counter({6: 10, 4: 2, 2: 1}),
        "quadratic": Counter({6: 10, 5: 2, 2: 1}),
        "generic": Counter({6: 12, 2: 1}),
    }
    expected_ranks = {
        "axis_anti": 50,
        "diagonal": 70,
        "quadratic": 72,
        "generic": 74,
    }

    for g in points[1:]:
        label = gain_class(g)
        class_counts[label] += 1
        block_ranks: Counter[int] = Counter()
        total_rank = 0
        flat_points: set[Point] = set()

        for line in lines_for_gain(g):
            values = [duty(q) for q in line]
            block = commutator_block(values)
            support, weights = punctures(values)
            require(block == factor_block(support, weights),
                    "puncture factorization failed")
            block_rank = rank_q(block)
            require(block_rank <= 2 * len(support),
                    "puncture rank bound failed")
            block_ranks[block_rank] += 1
            total_rank += block_rank
            line_blocks += 1

            constraints = constraint_matrix(support)
            constraint_rank = rank_q(constraints)
            stacked_rank = rank_q(constraints + block)
            require(stacked_rank == constraint_rank,
                    "commutator rows left the puncture constraint span")
            if block_rank == constraint_rank:
                if len(support) == 3 and block_rank == 5:
                    rank5_blocks += 1
            elif block_rank == 4:
                require(len(support) == 3 and constraint_rank == 5,
                        "unexpected rank-four kernel extension")
                adjacent = [(s, t) for s in support for t in support
                            if (t - s) % P == 1 and values[s] == values[t]]
                require(len(adjacent) == 1,
                        "rank-four block lacks its equal adjacent pair")
                s, t = adjacent[0]
                other = next(u for u in support if u not in (s, t))
                row_other = [0] * P
                row_other[other] = 1
                row_sum = [0] * P
                row_sum[s] = row_sum[t] = 1
                require(rank_q(block + [row_other]) == block_rank,
                        "rank-four kernel does not kill third puncture")
                require(rank_q(block + [row_sum]) == block_rank,
                        "rank-four kernel lacks opposite-sign restriction")
                rank4_blocks += 1
            else:
                require(block_rank == 5 and constraint_rank == 5,
                        "unexpected exceptional block")
                rank5_blocks += 1

            # Coordinate zero columns are exactly the forward flat arcs.
            for column, q in enumerate(line):
                zero_column = all(block[row][column] == 0 for row in range(P))
                is_flat = flat_arc(g, q)
                require(zero_column == is_flat,
                        "flat arc does not match coordinate kernel")
                if is_flat:
                    flat_points.add(q)

        require(block_ranks == expected_patterns[label],
                "wrong affine-line block pattern")
        require(total_rank == expected_ranks[label],
                "wrong full commutator rank")
        rank_histogram[total_rank] += 1
        nullity_histogram[P * P - total_rank] += 1

        flat_count = len(flat_points)
        flat_histogram[flat_count] += 1
        require((0, 0) not in flat_points, "origin entered flat-arc cone")
        if label == "axis_anti":
            require(flat_count == 36, "wrong tangent flat-arc count")
        else:
            slope = g[1] * pow(g[0], -1, P) % P
            expected_flat = 18 if slope in (1, 3, 6, 9, 11) else 16
            require(flat_count == expected_flat,
                    "wrong transverse flat-arc count")

        # The six uncontaminated zero-anchor rows.
        n0 = duty((0, 0))
        ng = duty(g)
        defect = n0 - ng
        expected_defect = (229_692 if g[0] * g[1] == 0 else
                           440_232 if (g[0] + g[1]) % P == 0 else
                           440_244)
        require(defect == expected_defect, "wrong anchor defect class")
        origin_line = [scale(k, g) for k in range(P)]
        origin_block = commutator_block([duty(q) for q in origin_line])
        for j in range(1, HALF + 1):
            nonzero = [(column, value)
                       for column, value in enumerate(origin_block[j])
                       if value]
            require(nonzero == [(0, -defect)],
                    "anchor row is contaminated")
            anchor_rows += 1

    require(line_blocks == 2_184, "wrong affine-line block count")
    require(class_counts == Counter({
        "axis_anti": 36, "diagonal": 12,
        "quadratic": 24, "generic": 96,
    }), "wrong gain-class census")
    require(rank_histogram == Counter({50: 36, 70: 12, 72: 24, 74: 96}),
            "wrong rank histogram")
    require(nullity_histogram
            == Counter({119: 36, 99: 12, 97: 24, 95: 96}),
            "wrong nullity histogram")
    require(flat_histogram == Counter({36: 36, 18: 60, 16: 72}),
            "wrong flat-arc histogram")
    require(rank4_blocks == 24 and rank5_blocks == 48,
            "wrong exceptional block counts")
    require(anchor_rows == 1_008, "wrong anchor-row count")
    require(Fraction(N7, FIBRE) == Fraction(205_668, 3_975_084_764_287),
            "duty normalization failed")

    print("THM-2562 duty-commutator holotopy referee")
    print("line_blocks=2184 triple_incidence_sets=286 exceptional_gap_sets=13")
    print("rank_histogram=50:36,70:12,72:24,74:96")
    print("nullity_histogram=119:36,99:12,97:24,95:96")
    print("block_patterns=axis_anti:2+12x4,diagonal:2+2x4+10x6,quadratic:2+2x5+10x6,generic:2+12x6")
    print("exceptional_slopes=1,6,11 schur_equal_rank=4 schur_unequal_rank=5")
    print("flat_arc_histogram=36:36,18:60,16:72 flat_arc_minimum=16 origin_flat=0")
    print("nonnegative_kernel=flat_arc_coordinate_cone")
    print("anchor_rows=1008 defects=229692,440232,440244")
    print("kappa=1439676/27825593350009 reduced=205668/3975084764287")
    print(f"explicit_require_checks={CHECKS}")
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
