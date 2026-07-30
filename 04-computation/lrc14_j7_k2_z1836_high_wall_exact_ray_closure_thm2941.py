#!/usr/bin/env python3
"""Exact all-label closure of the two projected k=2 z1=1836 high walls.

The horizon scan at this frontier left two rows whose forced largest label
was represented by an analytic ``HIGH-TAIL`` placeholder.  This verifier
removes the horizon altogether.  If ``C_E`` is the six-body safe carrier,
``L=14*lcm(E)``, and

    delta_E(z) = mu(C_E intersect D_z) - mu(C_E)/7,

then the endpoints of ``C_E`` lie on the ``1/L`` grid and periodicity of the
primitive of ``1_D1`` gives the exact residue-ray law

    delta_E(r+mL) = A_E(r)/(r+mL),        1 <= r < L.       (1)

The program checks (1) directly for every residue.  On a positive ray, only
its first four labels strictly after z1 can occur in the global top four:
every later point has four distinct, strictly larger predecessors on its own
ray.  Nonpositive rays cannot enter because the computed fourth value is
positive.  Likewise the best label satisfying the forced high wall is the
first high point on one of the positive rays.

For k=2 a countercover must have one of the four later drifts at least

    floor(13 L / 150) + 1,

and must satisfy the scalar necessary condition

    delta_E(z1) + ... + delta_E(z5) >= mu(C_E)/91.         (2)

For each assigned row the unconstrained top four lie below the high wall.
Thus the exact constrained maximum in (2) is the sum of the unconstrained
top three and the best high label.  Both exact maxima are strictly below the
right side of (2), so both rows are scalar-empty for *all* later labels.

As an independent geometric control, the script also reconstructs the
lossless projected residual of THM-2941(25i) for each maximizing packet.  In
both cases the common danger intersection is already empty after a short
prefix of the body-safe cells, so the projected residual has exact mass one
(and in particular at least 25/91).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import inspect
import math
from fractions import Fraction as F
from pathlib import Path


def repo_root(start: Path) -> Path:
    for parent in (start, *start.parents):
        if (parent / "04-computation").is_dir() and (parent / "01-canon").is_dir():
            return parent
    raise RuntimeError("cannot locate repository root")


HERE = Path(__file__).resolve().parent
ROOT = repo_root(HERE)
SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_critical_scalar_wall_independent_thm2941.py"
)
BAND_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1836_1836_thm2941.out"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1836_high_wall_exact_ray_closure_thm2941.out"
)
EXPECTED_SOURCE_SHA256 = (
    "5d25a955fe184d6c1a3d8b632b4bbf901dc996ee46ad67c5748836fcc7134404"
)
EXPECTED_SUPPORT_SHA256 = (
    "5482e10635ecf72840bc0c083360fd7ddad65c2885d743820061bcba58cd5609"
)
EXPECTED_BAND_OUTPUT_SHA256 = (
    "fdbeee8d14057edc7b6d1a5a24c660899846aef5eb2a7c6a3004367caf3aa73b"
)
EXPECTED_SEMANTIC_SHA256 = (
    "8b303e85838d5a8850e3550887ffade1499de07a0a0799d7121ddcb311f9f8ef"
)

FIRST = 1_836
SUFFIX_SLOTS = 4
HIGH_WALL_RATIO = F(13, 150)
SCALAR_ETA = F(1, 91)
PROJECTED_THRESHOLD = F(25, 91)
BODIES = (
    (1, 4, 10, 12, 13, 14),
    (1, 8, 10, 12, 13, 14),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def normalized_sha256(path: Path) -> str:
    payload = path.read_bytes()
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: lone carriage return",
    )
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


spec = importlib.util.spec_from_file_location("z1836_high_wall_source", SOURCE)
require(spec is not None and spec.loader is not None, "cannot load source")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)

SUPPORT_NAMES = (
    "merge_intervals",
    "danger_intervals",
    "carrier_for",
    "danger_primitive",
    "singleton_coverage",
)
support_payload = "\n".join(
    inspect.getsource(getattr(A, name)) for name in SUPPORT_NAMES
).encode()
require(
    normalized_sha256(SOURCE) == EXPECTED_SOURCE_SHA256,
    "pinned THM-2941 source changed",
)
require(
    normalized_sha256(BAND_OUTPUT) == EXPECTED_BAND_OUTPUT_SHA256,
    "pinned z1836 scalar slice changed",
)
require(
    hashlib.sha256(support_payload).hexdigest() == EXPECTED_SUPPORT_SHA256,
    "pinned carrier/singleton support changed",
)
require(A.RULER == 5_045_040 and A.BASE_LABEL == 15, "source constants changed")
require(HIGH_WALL_RATIO == F(13, 150), "projected wall changed")
require(SCALAR_ETA == F(1, 91), "k=2 scalar surplus changed")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def excess_amplitude(
    carrier: tuple[tuple[int, int], ...],
    mass_numerator: int,
    label: int,
) -> F:
    """Return ``label * delta_E(label)`` exactly."""

    covered_numerator = sum(
        A.danger_primitive(label * right) - A.danger_primitive(label * left)
        for left, right in carrier
    )
    return F(7 * covered_numerator - mass_numerator * label, 7 * A.RULER)


def first_strictly_after(residue: int, bound: int, period: int) -> int:
    if residue > bound:
        return residue
    return residue + ((bound - residue) // period + 1) * period


def first_at_least(residue: int, bound: int, period: int) -> int:
    if residue >= bound:
        return residue
    return residue + ((bound - residue + period - 1) // period) * period


def merge_fraction(
    intervals: list[tuple[F, F]],
) -> tuple[tuple[F, F], ...]:
    rows = sorted((left, right) for left, right in intervals if left < right)
    if not rows:
        return ()
    merged: list[list[F]] = [[rows[0][0], rows[0][1]]]
    for left, right in rows[1:]:
        if left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    return tuple((left, right) for left, right in merged)


def intersect_fraction(
    first: tuple[tuple[F, F], ...],
    second: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    rows: list[tuple[F, F]] = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            rows.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(rows)


def interval_mass(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def phase_danger(cell: int, speed: int, period: int) -> tuple[tuple[F, F], ...]:
    """Danger subset in the normalized coordinate on one body cell."""

    start = F(speed * cell, period)
    finish = F(speed * (cell + 1), period)
    rows: list[tuple[F, F]] = []
    for tooth in range(floor_fraction(start) - 1, ceil_fraction(finish) + 2):
        left = max(
            F(0),
            F(period, speed) * (F(tooth) - F(1, 14)) - cell,
        )
        right = min(
            F(1),
            F(period, speed) * (F(tooth) + F(1, 14)) - cell,
        )
        if left < right:
            rows.append((left, right))
    return merge_fraction(rows)


def body_cells(
    carrier: tuple[tuple[int, int], ...], period: int
) -> tuple[int, ...]:
    require(A.RULER % period == 0, "body period does not divide ruler")
    scale = A.RULER // period
    cells: list[int] = []
    for left, right in carrier:
        require(
            left % scale == 0 and right % scale == 0,
            "carrier endpoint left the body grid",
        )
        cells.extend(range(left // scale, right // scale))
    require(len(cells) == len(set(cells)) > 0, "body-safe cell ledger changed")
    return tuple(cells)


def projected_residual(
    carrier: tuple[tuple[int, int], ...],
    period: int,
    packet: tuple[int, ...],
) -> tuple[F, int, F, int, int]:
    """Exact residual mass and first prefix attaining the k=2 threshold."""

    cells = body_cells(carrier, period)
    common: tuple[tuple[F, F], ...] = ((F(0), F(1)),)
    threshold_at = 0
    threshold_mass = F(0)
    exhausted_at = 0
    for used, cell in enumerate(cells, 1):
        local_union = merge_fraction(
            [
                interval
                for speed in packet
                for interval in phase_danger(cell, speed, period)
            ]
        )
        common = intersect_fraction(common, local_union)
        residual_mass = 1 - interval_mass(common)
        if threshold_at == 0 and residual_mass >= PROJECTED_THRESHOLD:
            threshold_at = used
            threshold_mass = residual_mass
        if not common:
            exhausted_at = used
            break
    require(threshold_at > 0, "projected threshold was not attained")
    if exhausted_at == 0:
        exhausted_at = len(cells)
    return 1 - interval_mass(common), threshold_at, threshold_mass, exhausted_at, len(cells)


def profile(body: tuple[int, ...]) -> tuple[object, ...]:
    carrier = A.carrier_for(body)
    mass_numerator = sum(right - left for left, right in carrier)
    h = F(mass_numerator, A.RULER)
    components = len(carrier)
    period = 14 * math.lcm(*body)
    high_floor = max(A.BASE_LABEL, floor_fraction(HIGH_WALL_RATIO * period) + 1)

    scale = A.RULER // period
    require(
        A.RULER % period == 0
        and all(left % scale == right % scale == 0 for left, right in carrier),
        (body, "carrier is off the body ruler"),
    )

    amplitudes: list[F] = [F(0)] * period
    recurrence_checks = 0
    for residue in range(1, period):
        amplitude = excess_amplitude(carrier, mass_numerator, residue)
        next_amplitude = excess_amplitude(carrier, mass_numerator, residue + period)
        require(next_amplitude == amplitude, (body, "ray recurrence", residue))
        amplitudes[residue] = amplitude
        recurrence_checks += 1
    require(
        excess_amplitude(carrier, mass_numerator, period) == 0,
        (body, "aligned amplitude is nonzero"),
    )
    for residue in range(1, period):
        require(
            amplitudes[period - residue] == -amplitudes[residue],
            (body, "residue antipode", residue),
        )

    arbitrary: list[tuple[F, int, int, int]] = []
    high: list[tuple[F, int, int]] = []
    fifth_by_ray: list[tuple[F, int, int]] = []
    positive_residues = 0
    for residue in range(1, period):
        amplitude = amplitudes[residue]
        if amplitude <= 0:
            continue
        positive_residues += 1
        first_label = first_strictly_after(residue, FIRST, period)
        for step in range(SUFFIX_SLOTS):
            label = first_label + step * period
            arbitrary.append((amplitude / label, label, residue, step))
        fifth_label = first_label + SUFFIX_SLOTS * period
        fifth_by_ray.append((amplitude / fifth_label, fifth_label, residue))
        high_label = first_at_least(residue, high_floor, period)
        high.append((amplitude / high_label, high_label, residue))

    rank4 = tuple(sorted(arbitrary, key=lambda row: (-row[0], row[1:]))[:4])
    require(len(rank4) == 4 and rank4[-1][0] > 0, (body, "positive top four missing"))
    omitted_max = min(fifth_by_ray, key=lambda row: (-row[0], row[1:]))
    require(
        omitted_max[0] <= rank4[-1][0],
        (body, "first-four-per-ray truncation failed", omitted_max, rank4[-1]),
    )
    best_high = min(high, key=lambda row: (-row[0], row[1:]))

    # Here the unconstrained maximizer violates the wall.  Therefore every
    # feasible four-set has one high item bounded by best_high and three
    # remaining items bounded by the unconstrained top-three sum.
    require(
        all(label < high_floor for _value, label, _residue, _step in rank4),
        (body, "unconstrained top four already meet the wall"),
    )
    require(best_high[1] >= high_floor, (body, "high label below wall"))
    constrained = (*rank4[:3], (best_high[0], best_high[1], best_high[2], 0))
    require(len({row[1] for row in constrained}) == 4, (body, "label collision"))

    first_delta = amplitudes[FIRST % period] / FIRST
    scalar_upper = first_delta + sum((row[0] for row in constrained), F(0))
    scalar_lower = h * SCALAR_ETA
    scalar_gap = scalar_upper - scalar_lower
    require(scalar_gap < 0, (body, "high-wall scalar row survived", scalar_gap))

    packet = tuple(sorted((FIRST, *(row[1] for row in constrained))))
    residual = projected_residual(carrier, period, packet)
    require(residual[0] == 1, (body, "maximizing projected residual changed"))

    return (
        body,
        period,
        h,
        components,
        high_floor,
        recurrence_checks,
        positive_residues,
        rank4,
        omitted_max,
        best_high,
        first_delta,
        tuple(constrained),
        packet,
        scalar_upper,
        scalar_lower,
        scalar_gap,
        residual,
    )


def render(rows: tuple[tuple[object, ...], ...]) -> str:
    require(tuple(row[0] for row in rows) == BODIES, "body universe changed")
    semantic_payload = (
        EXPECTED_SOURCE_SHA256,
        EXPECTED_BAND_OUTPUT_SHA256,
        FIRST,
        SUFFIX_SLOTS,
        HIGH_WALL_RATIO,
        SCALAR_ETA,
        PROJECTED_THRESHOLD,
        rows,
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_hash == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest changed", semantic_hash),
        )

    lines = [
        "LRC14 projected k=2 z1=1836 high-wall exact ray closure",
        f"source_sha256={normalized_sha256(SOURCE)}",
        f"scalar_band_output_sha256={normalized_sha256(BAND_OUTPUT)}",
        (
            "scope=two inherited HIGH-TAIL body rows;all distinct later labels;"
            "no finite label horizon"
        ),
        (
            "ray_law=delta(r+mL)=A(r)/(r+mL);"
            "search=first four points of every positive ray plus first high point"
        ),
        (
            "forced_high_wall=floor(13L/150)+1;"
            "scalar_necessary=sum(delta)>=h/91;projected_threshold=25/91"
        ),
    ]
    for row in rows:
        (
            body,
            period,
            h,
            components,
            high_floor,
            recurrence_checks,
            positive_residues,
            rank4,
            omitted_max,
            best_high,
            first_delta,
            constrained,
            packet,
            scalar_upper,
            scalar_lower,
            scalar_gap,
            residual,
        ) = row
        lines.extend(
            (
                (
                    f"ROW;E={','.join(map(str, body))};L={period};h={ftext(h)};"
                    f"r={components};high_floor={high_floor};"
                    f"positive_residues={positive_residues};"
                    f"recurrence_checks={recurrence_checks}"
                ),
                "  unconstrained_top4="
                + ",".join(f"{item[1]}:{ftext(item[0])}" for item in rank4),
                (
                    f"  first_omitted_ray_point={omitted_max[1]}:"
                    f"{ftext(omitted_max[0])};"
                    f"fourth_top={rank4[-1][1]}:{ftext(rank4[-1][0])}"
                ),
                f"  best_high={best_high[1]}:{ftext(best_high[0])}",
                (
                    "  exact_constrained_suffix="
                    + ",".join(f"{item[1]}:{ftext(item[0])}" for item in constrained)
                ),
                (
                    f"  delta1={ftext(first_delta)};upper={ftext(scalar_upper)};"
                    f"lower={ftext(scalar_lower)};gap={ftext(scalar_gap)}"
                ),
                (
                    f"  projected_packet={packet};exact_residual={ftext(residual[0])};"
                    f"threshold_prefix_cells={residual[1]};"
                    f"threshold_prefix_mass={ftext(residual[2])};"
                    f"common_danger_empty_at={residual[3]};body_safe_cells={residual[4]}"
                ),
                "  conclusion=SCALAR-EMPTY",
            )
        )
    lines.extend(
        (
            "global_high_wall_rows=2;scalar_empty=2;survivors=0",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    rows = tuple(profile(body) for body in BODIES)
    output = render(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
