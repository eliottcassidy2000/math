#!/usr/bin/env python3
"""Independent exact checks for THM-2186's octagon-needle profile."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256


SWEEP_LIMIT = 512
EXCEPTIONAL = (2, 3, 4, 5, 6, 10, 11, 12, 18)
HOSTILE = (1024, 4095, 4096, 4097, 10000, 10001, 10002)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def speeds(c: int) -> tuple[int, ...]:
    return (1,) + tuple(
        value for k in range(1, 7) for value in (k * c, k * c + 1)
    )


def circle_distance(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def profile_at(c: int, time: Fraction) -> Fraction:
    return min(circle_distance(a * time) for a in speeds(c))


def deficit(c: int) -> Fraction:
    residue = c % 8
    if residue == 1:
        return Fraction(0)
    if residue == 2:
        return Fraction(3, 8 * (3 * c + 1))
    numerator = {0: 1, 3: 6, 4: 5, 5: 4, 6: 3, 7: 2}[residue]
    return Fraction(numerator, 8 * (7 * c + 1))


def predicted(c: int) -> Fraction:
    return Fraction(1, 8) - deficit(c)


def explicit_witness(c: int) -> Fraction:
    residue = c % 8
    margin_deficit = deficit(c)
    if residue == 0:
        return Fraction(7, 8) - 7 * margin_deficit
    if residue == 1:
        return Fraction(1, 8)
    if residue == 2:
        return Fraction(7, 8) + margin_deficit
    if residue == 3:
        return Fraction(1, 8) + 7 * margin_deficit
    if residue == 4:
        return Fraction(3, 8) - Fraction(7, 5) * margin_deficit
    if residue == 5:
        return Fraction(1, 8) + 7 * margin_deficit
    if residue == 6:
        return Fraction(5, 8) - Fraction(7, 3) * margin_deficit
    return Fraction(7, 8) - 7 * margin_deficit


def affine_cell_maximum(c: int) -> tuple[Fraction, Fraction]:
    """Independent evaluator: exact affine cells and in-cell line crossings."""

    row = speeds(c)
    boundaries = sorted(
        {Fraction(k, 2 * a) for a in row for k in range(2 * a + 1)}
    )
    best = (Fraction(-1), Fraction(0))

    for left, right in zip(boundaries, boundaries[1:]):
        midpoint = (left + right) / 2
        lines: list[tuple[int, int]] = []
        for a in row:
            scaled = a * midpoint
            integer_part = scaled.numerator // scaled.denominator
            fractional_part = scaled - integer_part
            if fractional_part < Fraction(1, 2):
                lines.append((a, -integer_part))
            else:
                lines.append((-a, integer_part + 1))

        candidates = {left, right}
        for index, (slope, intercept) in enumerate(lines):
            for other_slope, other_intercept in lines[index + 1 :]:
                if slope == other_slope:
                    continue
                crossing = Fraction(
                    other_intercept - intercept, slope - other_slope
                )
                if left <= crossing <= right:
                    candidates.add(crossing)

        for time in candidates:
            candidate = (profile_at(c, time), time)
            if candidate > best:
                best = candidate

    return best


def integer_breakpoint_maximum(c: int) -> tuple[Fraction, Fraction]:
    """Independent evaluator using only integer residues at all breakpoints."""

    row = speeds(c)
    denominators = {2 * a for a in row}
    for index, a in enumerate(row):
        for b in row[index + 1 :]:
            denominators.add(a + b)
            if a != b:
                denominators.add(abs(a - b))

    best_numerator = 0
    best_denominator = 1
    best_time = Fraction(0)

    for denominator in sorted(denominators):
        for numerator in range(denominator + 1):
            minimum_residue = denominator
            for a in row:
                residue = (a * numerator) % denominator
                residue = min(residue, denominator - residue)
                if residue < minimum_residue:
                    minimum_residue = residue
                if (
                    minimum_residue * best_denominator
                    < best_numerator * denominator
                ):
                    break

            comparison = (
                minimum_residue * best_denominator
                - best_numerator * denominator
            )
            if comparison > 0:
                best_numerator = minimum_residue
                best_denominator = denominator
                best_time = Fraction(numerator, denominator)
            elif comparison == 0:
                time = Fraction(numerator, denominator)
                if time > best_time:
                    best_time = time

    return Fraction(best_numerator, best_denominator), best_time


def projective_height(c: int) -> Fraction | None:
    gap = deficit(c)
    if gap == 0:
        return None
    return predicted(c) / gap


def main() -> None:
    exception_lines: list[str] = []
    for c in EXCEPTIONAL:
        arrangement_value, arrangement_time = affine_cell_maximum(c)
        breakpoint_value, breakpoint_time = integer_breakpoint_maximum(c)
        require(
            arrangement_value == breakpoint_value,
            f"independent evaluators disagree at c={c}",
        )
        require(
            arrangement_value == predicted(c),
            f"exceptional formula failure at c={c}",
        )
        require(
            profile_at(c, explicit_witness(c)) == predicted(c),
            f"exceptional witness failure at c={c}",
        )
        exception_lines.append(
            f"{c}:{arrangement_value}:{arrangement_time}:"
            f"{breakpoint_time}"
        )

    transcript = sha256()
    for c in range(1, SWEEP_LIMIT + 1):
        value, time = integer_breakpoint_maximum(c)
        expected = predicted(c)
        witness = explicit_witness(c)
        require(value == expected, f"sweep formula failure at c={c}")
        require(
            profile_at(c, witness) == expected,
            f"sweep witness failure at c={c}",
        )
        transcript.update(f"{c}:{value}:{time}:{witness}\n".encode())

    hostile_lines: list[str] = []
    for c in HOSTILE:
        value, time = integer_breakpoint_maximum(c)
        require(value == predicted(c), f"hostile formula failure at c={c}")
        require(
            profile_at(c, explicit_witness(c)) == value,
            f"hostile witness failure at c={c}",
        )
        hostile_lines.append(f"{c}:{value}:{time}")

    for c in range(8, SWEEP_LIMIT + 1, 8):
        margin = predicted(c)
        height = projective_height(c)
        require(height == 7 * c, f"projective height failure at c={c}")
        for scale in (2, 3):
            scaled_margin = predicted(scale * c)
            mobius = Fraction(scale) * margin / (
                1 + 8 * (scale - 1) * margin
            )
            require(
                scaled_margin == mobius,
                f"Mobius scale failure at c={c}, q={scale}",
            )

    expected_heights = {
        0: lambda c: Fraction(7 * c),
        2: lambda c: Fraction(3 * c - 2, 3),
        3: lambda c: Fraction(7 * c - 5, 6),
        4: lambda c: Fraction(7 * c - 4, 5),
        5: lambda c: Fraction(7 * c - 3, 4),
        6: lambda c: Fraction(7 * c - 2, 3),
        7: lambda c: Fraction(7 * c - 1, 2),
    }
    for c in range(1, SWEEP_LIMIT + 1):
        if c % 8 == 1:
            require(projective_height(c) is None, "critical H must be infinite")
        else:
            require(
                projective_height(c) == expected_heights[c % 8](c),
                f"residue H failure at c={c}",
            )

    proposed_matches = tuple(
        c
        for c in range(1, SWEEP_LIMIT + 1)
        if predicted(c) == Fraction(7 * c, 8 * (7 * c + 1))
    )
    require(
        proposed_matches == tuple(range(8, SWEEP_LIMIT + 1, 8)),
        "claimed-formula range is not exactly 8|c",
    )

    print("THM-2186 exact octagon-needle checks")
    print(f"exceptional={EXCEPTIONAL}")
    print(f"exceptional_exact={';'.join(exception_lines)}")
    print(f"sweep=1..{SWEEP_LIMIT}")
    print(f"sweep_sha256={transcript.hexdigest()}")
    print(f"hostile_exact={';'.join(hostile_lines)}")
    print("proposed_formula_iff=8|c")
    print("dyadic_values=c2:1/14,c4:3/29,c8+:7c/[8(7c+1)]")
    print("projective_height_on_8Z=7c")
    print("mobius_scale=M(qc)=qM(c)/[1+8(q-1)M(c)]")
    print("independent_evaluators=affine_cells,integer_breakpoints")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
