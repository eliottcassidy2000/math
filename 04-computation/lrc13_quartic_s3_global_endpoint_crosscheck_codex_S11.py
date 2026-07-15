#!/usr/bin/env python3
"""Independent endpoint-cell cross-check for THM-816.

The analytic recursive comb bound in THM-816 proves that every hypothetical
covering row has all four off-sheet speeds at most 398.  This second program
does not use the recursive interval traversal.  For each of the twelve
coset/parity contexts it cuts the circle at every threshold endpoint arising
from the core and every permitted exception speed through 398.  Threshold
membership is constant on each resulting open cell.  Python integer bitsets
then test every one of the 132,510 four-way congruence rows directly.

Fractions are exact; no floating point or sampled grid is used.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import product


P = 13
DELTA = F(1, P)
CAP = 398
H = frozenset((1, 5, 8, 12))
BASE = tuple(range(1, P))


def circle_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def cosets() -> tuple[tuple[int, ...], ...]:
    rows: list[tuple[int, ...]] = []
    for a in BASE:
        row = tuple(sorted({a * h % P for h in H}))
        if row not in rows:
            rows.append(row)
    assert len(rows) == 3
    return tuple(rows)


def negative_pairs(labels: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    remaining = set(labels)
    rows: list[tuple[int, int]] = []
    while remaining:
        r = min(remaining)
        pair = tuple(sorted((r, -r % P)))
        assert len(pair) == 2 and set(pair) <= remaining
        rows.append(pair)
        remaining -= set(pair)
    assert len(rows) == 2
    return tuple(rows)


def crt_base(label: int, parity: int) -> int:
    values = [
        u
        for u in range(1, 3 * P)
        if u % P == 3 * label % P and u % 3 == parity
    ]
    assert len(values) == 1
    return values[0]


def threshold_endpoints(speed: int) -> tuple[F, ...]:
    return tuple(
        endpoint
        for k in range(speed)
        for endpoint in (
            F(P * k + 1, P * speed),
            F(P * (k + 1) - 1, P * speed),
        )
    )


def main() -> None:
    contexts = 0
    total_rows = 0
    total_covers = 0
    cell_counts: list[int] = []
    row_counts: list[int] = []
    lines = [
        "quartic_s3 independent global endpoint crosscheck",
        f"analytic_global_cap={CAP}",
    ]
    digest = sha256()

    for labels in cosets():
        pairs = negative_pairs(labels)
        core = tuple(3 * q for q in BASE if q not in labels)
        for pair_bits in product((1, 2), repeat=2):
            parity = {
                r: pair_bits[index]
                for index, pair in enumerate(pairs)
                for r in pair
            }
            choices: dict[int, tuple[int, ...]] = {}
            for r in labels:
                first = crt_base(r, parity[r])
                choices[r] = tuple(range(first, CAP + 1, 39))
            exception_bank = tuple(u for r in labels for u in choices[r])

            endpoints = {F(0), F(1)}
            for speed in (*core, *exception_bank):
                endpoints.update(threshold_endpoints(speed))
            points = sorted(endpoints)
            midpoints = tuple(
                (left + right) / 2
                for left, right in zip(points, points[1:])
                if left < right
            )

            core_safe = 0
            danger = {speed: 0 for speed in exception_bank}
            for index, time in enumerate(midpoints):
                bit = 1 << index
                if all(circle_distance(speed * time) > DELTA for speed in core):
                    core_safe |= bit
                for speed in exception_bank:
                    if circle_distance(speed * time) <= DELTA:
                        danger[speed] |= bit
            assert core_safe

            local_rows = 0
            local_covers = 0
            for speeds in product(*(choices[r] for r in labels)):
                local_rows += 1
                cover = 0
                for speed in speeds:
                    cover |= danger[speed]
                residual = core_safe & ~cover
                if not residual:
                    local_covers += 1
                first_uncovered = (residual & -residual).bit_length() - 1
                digest.update(
                    repr((labels, pair_bits, speeds, first_uncovered)).encode()
                )
                digest.update(b"\n")

            contexts += 1
            total_rows += local_rows
            total_covers += local_covers
            cell_counts.append(len(midpoints))
            row_counts.append(local_rows)
            lines.append(
                f"context={contexts} labels={labels} parity={pair_bits} "
                f"cells={len(midpoints)} core_safe_cells={core_safe.bit_count()} "
                f"rows={local_rows} covers={local_covers}"
            )

    assert contexts == 12
    assert row_counts == [
        10_000,
        11_000,
        10_000,
        11_000,
        12_100,
        11_000,
        11_000,
        10_000,
        11_000,
        13_310,
        10_000,
        12_100,
    ]
    assert cell_counts == [
        15_617,
        15_273,
        16_025,
        15_681,
        16_405,
        16_033,
        16_165,
        15_793,
        16_205,
        16_817,
        15_989,
        16_601,
    ]
    assert total_rows == 132_510 and total_covers == 0
    certificate_digest = digest.hexdigest()
    assert certificate_digest == "a02c30e784969c5865606455a38843f44ff49fb6a1b96faee2548a64204f8b83"

    lines.extend(
        [
            f"contexts={contexts}",
            f"rows={total_rows}",
            f"covers={total_covers}",
            f"cell_count_range=({min(cell_counts)},{max(cell_counts)})",
            f"certificate_digest={certificate_digest}",
        ]
    )
    output_digest = sha256(("\n".join(lines) + "\n").encode()).hexdigest()
    lines.append(f"sha256={output_digest}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
