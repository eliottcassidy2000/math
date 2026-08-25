#!/usr/bin/env python3
"""Fraction-exact primary audit for THM-4098.

The proof is analytic.  This referee checks its finite boundary data and a
declared hostile census: every row-compatible bank drawn from {1,...,10}, at
the first sufficient dilation and the following (even) dilation.  Literal
open danger intervals are merged exactly; no floating point or ``assert`` is
used.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations


DELTA = F(1, 14)
ROWS = {
    # core: (left, right, number of outliers, odd outliers, first q)
    7: (F(4, 35), F(13, 98), 4, 3, 55),
    6: (F(4, 35), F(1, 7), 5, 2, 35),
    5: (F(4, 35), F(1, 7), 6, 1, 35),
}
SAMPLE_BANKS = {
    7: (8, 9, 11, 13),
    6: (7, 8, 9, 10, 12),
    5: (6, 8, 10, 12, 14, 15),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def omega(speed: int) -> int:
    return 1 if speed % 2 == 0 else 2


def circle_norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def two_phase_safe(speed: int, theta: F) -> bool:
    return all(
        circle_norm(speed * (theta + F(phase, 2))) >= DELTA
        for phase in (0, 1)
    )


def two_phase_teeth(speed: int, left: F, right: F) -> list[tuple[F, F, int]]:
    """Lifted literal open teeth of U_speed meeting [left,right]."""
    count = omega(speed) * speed
    radius = DELTA / speed
    lo = count * (left - radius)
    hi = count * (right + radius)
    first = lo.numerator // lo.denominator - 2
    last = hi.numerator // hi.denominator + 3
    return [
        (F(index, count) - radius, F(index, count) + radius, speed)
        for index in range(first, last + 1)
        if F(index, count) + radius > left
        and F(index, count) - radius < right
    ]


def single_phase_teeth(speed: int, left: F, right: F) -> list[tuple[F, F, int]]:
    radius = DELTA / speed
    lo = speed * (left - radius)
    hi = speed * (right + radius)
    first = lo.numerator // lo.denominator - 2
    last = hi.numerator // hi.denominator + 3
    return [
        (F(index, speed) - radius, F(index, speed) + radius, speed)
        for index in range(first, last + 1)
        if F(index, speed) + radius > left
        and F(index, speed) - radius < right
    ]


def merged_union(
    intervals: list[tuple[F, F, int]], left: F, right: F
) -> list[tuple[F, F]]:
    clipped = sorted(
        (max(left, start), min(right, end))
        for start, end, _ in intervals
        if end > left and start < right
    )
    merged: list[tuple[F, F]] = []
    for start, end in clipped:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def safe_components(
    bank: tuple[int, ...], *, single_phase: bool = False
) -> list[tuple[F, F]]:
    teeth = single_phase_teeth if single_phase else two_phase_teeth
    danger: list[tuple[F, F, int]] = []
    for speed in bank:
        danger.extend(teeth(speed, F(0), F(1)))
    merged = merged_union(danger, F(0), F(1))
    safe: list[tuple[F, F]] = []
    cursor = F(0)
    for start, end in merged:
        if cursor < start:
            safe.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < 1:
        safe.append((cursor, F(1)))
    return safe


def uncovered_measure(bank: tuple[int, ...], left: F, right: F) -> F:
    danger: list[tuple[F, F, int]] = []
    for speed in bank:
        danger.extend(two_phase_teeth(speed, left, right))
    covered = sum((end - start for start, end in merged_union(danger, left, right)), F(0))
    return right - left - covered


def literally_covers(bank: tuple[int, ...], left: F, right: F) -> bool:
    """Whether literal open teeth contain every point of [left,right]."""
    intervals: list[tuple[F, F, int]] = []
    for speed in bank:
        intervals.extend(two_phase_teeth(speed, left, right))
    active: F | None = None
    for start, end, _ in sorted(intervals):
        if active is None:
            if start < left:
                active = end
            continue
        if start < active:
            active = max(active, end)
        elif active <= right:
            return False
        if active > right:
            return True
    return active is not None and active > right


def even_presentation(value: F) -> int:
    return value.denominator if value.denominator % 2 == 0 else 2 * value.denominator


def arrangement_witness(
    core: int, bank: tuple[int, ...], left: F, right: F
) -> tuple[F, int, str]:
    body = tuple(range(1, core + 1)) + bank
    candidates: list[tuple[F, int, str]] = [
        (left, even_presentation(left), "core-left"),
        (right, even_presentation(right), "core-right"),
    ]
    for speed in bank:
        for start, end, _ in two_phase_teeth(speed, left, right):
            if left <= start <= right:
                candidates.append((start, 14 * speed, f"{speed}-left"))
            if left <= end <= right:
                candidates.append((end, 14 * speed, f"{speed}-right"))
    for theta, clock, source in candidates:
        if all(two_phase_safe(speed, theta) for speed in body):
            require((theta * clock).denominator == 1, "bad endpoint presentation")
            return theta, clock, source
    raise RuntimeError(f"no arrangement witness for AP{core}, bank={bank}")


def mod_norm(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def main() -> None:
    # Exact inherited core intervals and arithmetic threshold gates.
    for core, (left, right, count, odd_count, first_q) in ROWS.items():
        length = right - left
        require(first_q * length >= 1, f"AP{core} mesh threshold fails")
        require((first_q - 1) * length < 1, f"AP{core} q is not first integer threshold")
        require(count + odd_count == 7, f"AP{core} is not the W=7 row")
        for speed in range(1, core + 1):
            require(two_phase_safe(speed, left), f"AP{core} left endpoint unsafe")
            require(two_phase_safe(speed, right), f"AP{core} right endpoint unsafe")

    # The exact THM-4092 hostile really covers AP7 before scale escape.
    hostile = SAMPLE_BANKS[7]
    left7, right7, *_ = ROWS[7]
    require(uncovered_measure(hostile, left7, right7) == 0, "AP7 hostile stopped covering")
    require(literally_covers(hostile, left7, right7), "AP7 hostile left a wall survivor")

    census_banks = 0
    census_scaled_rows = 0
    smallest_margin: F | None = None
    base_interval_floor: F | None = None
    for core, (left, right, count, odd_count, first_q) in ROWS.items():
        length = right - left
        for bank in combinations(range(1, 11), count):
            if sum(speed % 2 for speed in bank) != odd_count:
                continue
            require(sum(omega(speed) for speed in bank) == 7, "census weight changed")

            odd_safe = safe_components(bank)
            even_safe = safe_components(bank, single_phase=True)
            require(odd_safe, "weight-seven base has no two-phase interval")
            require(even_safe, "sub-seven single-phase base has no interval")
            parity_interval_floor = min(
                max(end - start for start, end in odd_safe),
                max(end - start for start, end in even_safe),
            )
            base_interval_floor = (
                parity_interval_floor
                if base_interval_floor is None
                else min(base_interval_floor, parity_interval_floor)
            )

            for q in (first_q, first_q + 1):
                scaled = tuple(q * speed for speed in bank)
                margin = uncovered_measure(scaled, left, right)
                require(margin > 0, f"scaled row covered AP{core}: {bank}, q={q}")
                smallest_margin = margin if smallest_margin is None else min(smallest_margin, margin)
                # Direct mesh-side check using the parity-appropriate base interval.
                components = odd_safe if q % 2 else even_safe
                ell = max(end - start for start, end in components)
                require((1 - ell) / q < length, "mesh gap did not clear core interval")
                census_scaled_rows += 1
            census_banks += 1

    # Compile explicit endpoint clocks and replay the complementary odd-tail gate.
    witness_rows: list[tuple[int, int, tuple[int, ...], F, int, str]] = []
    owner_gates = 0
    for core, base in SAMPLE_BANKS.items():
        left, right, count, odd_count, first_q = ROWS[core]
        require(len(base) == count and sum(v % 2 for v in base) == odd_count, "sample typing")
        for q in (first_q, first_q + 1):
            bank = tuple(q * speed for speed in base)
            theta, clock, source = arrangement_witness(core, bank, left, right)
            require(clock % 2 == 0, "owner clock is odd")
            require(clock <= 14 * max(bank), "owner clock exceeds endpoint bound")
            label = int(theta * clock) % clock
            opposite = (label + clock // 2) % clock
            for tail in range(1, 2 * clock, 2):
                require(
                    not (
                        7 * mod_norm(tail * label, clock) < clock
                        and 7 * mod_norm(tail * opposite, clock) < clock
                    ),
                    "odd tail survived both complementary labels",
                )
                owner_gates += 1
            witness_rows.append((core, q, bank, theta, clock, source))

    require(smallest_margin is not None and smallest_margin > 0, "missing census margin")
    require(base_interval_floor is not None and base_interval_floor > 0, "missing base interval")

    print("THM-4098 primary Fraction-exact audit: PASS")
    print(f"declared base-bank census={census_banks}; scaled rows={census_scaled_rows}")
    print(f"smallest direct uncovered margin={smallest_margin}")
    print(f"smallest parity-branch longest base interval in census={base_interval_floor}")
    print(f"owner complement gates={owner_gates}")
    print("row thresholds: " + ", ".join(
        f"AP{core}:q>={data[4]} (L={data[1]-data[0]})" for core, data in sorted(ROWS.items(), reverse=True)
    ))
    print("sample endpoint witnesses:")
    for core, q, bank, theta, clock, source in witness_rows:
        print(f"  AP{core} q={q} bank={bank} theta={theta} N={clock} source={source}")


if __name__ == "__main__":
    main()
