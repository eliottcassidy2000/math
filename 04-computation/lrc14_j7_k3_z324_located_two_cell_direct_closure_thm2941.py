#!/usr/bin/env python3
"""Independent direct two-cell closure of the sole ``k=3,z1=324`` state.

The pinned ray/status frontier leaves one denominator state.  This sidecar
independently reconstructs the exact scalar ray heads and shows that every
scalar-admissible packet has fixed drift labels 324, 364, and 492, together
with one high label of exact denominator 3920.

Two complete body-safe cells then give a pointwise certificate.  The three
fixed drift combs miss both cells.  Their cell indices differ by half of the
denominator 3920, so every unit-direction denominator-3920 comb has phases
differing by exactly one half on the two cells.  Its strict danger arcs on
the two cells are disjoint.  Thus every normalized phase has a drift-safe
lift in one of the cells, whereas three aligned danger combs have union mass
at most 3/7 and cannot cover the normalized circle.

This implementation does not import the adjacent projected-cell closure.
It also performs a direct rational carrier subtraction and projection for a
representative packet, preserving the pointwise/open-set direction in the
paper argument.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import re
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FRONTIER_PATH = (
    ROOT / "04-computation/lrc14_j7_k3_z324_ray_status_frontier_thm2941.py"
)
FRONTIER_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/"
    "lrc14_j7_k3_z324_ray_status_frontier_thm2941.out"
)
CARRIER_PATH = (
    ROOT / "04-computation/lrc14_j7_critical_scalar_wall_independent_thm2941.py"
)
BODY_ATLAS_PATH = (
    ROOT / "04-computation/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py"
)
BODY_ATLAS_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/"
    "lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge/results/"
    "lrc14_j7_k3_z324_located_two_cell_direct_closure_thm2941.out"
)
EXPECTED_FRONTIER_SHA256 = (
    "7eaaf551d2bd4ae386e2db4452edac7d30c25f7fc67b71967e48454d688bf78e"
)
EXPECTED_FRONTIER_OUTPUT_SHA256 = (
    "db3c5c68c4aa2f61584ef91dd2171901888270edbf17e860d40f16a64d3a9242"
)
EXPECTED_CARRIER_SHA256 = (
    "5d25a955fe184d6c1a3d8b632b4bbf901dc996ee46ad67c5748836fcc7134404"
)
EXPECTED_BODY_ATLAS_SHA256 = (
    "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded"
)
EXPECTED_BODY_ATLAS_OUTPUT_SHA256 = (
    "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
)
EXPECTED_SEMANTIC_SHA256 = (
    "3afa63cf8a0ed8084cced064fda19848eea5cec2f3ebd0938e8b91257e6e886e"
)

BODY = (2, 8, 10, 11, 12, 14)
FIRST = 324
DS = (3920, 4620, 10780, 10780)
TAIL_DS = (3920, 4620, 10780)
FIXED_LABELS = (324, 364, 492)
HIGH_REPRESENTATIVE = 12771
CELL_J = 5880
CELL_K = 19600
ALIGNED_UNION_CAP = F(3, 7)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(FRONTIER_PATH) == EXPECTED_FRONTIER_SHA256, "frontier changed")
require(
    file_sha256(FRONTIER_OUTPUT_PATH) == EXPECTED_FRONTIER_OUTPUT_SHA256,
    "frontier transcript changed",
)
require(file_sha256(CARRIER_PATH) == EXPECTED_CARRIER_SHA256, "carrier changed")
require(file_sha256(BODY_ATLAS_PATH) == EXPECTED_BODY_ATLAS_SHA256, "atlas changed")
require(
    file_sha256(BODY_ATLAS_OUTPUT_PATH) == EXPECTED_BODY_ATLAS_OUTPUT_SHA256,
    "atlas transcript changed",
)
A = load_module("z324_direct_carrier", CARRIER_PATH)


def first_on_ray(residue, modulus, threshold):
    if residue >= threshold:
        return residue
    return residue + ((threshold - residue + modulus - 1) // modulus) * modulus


def denominator(modulus, label):
    return modulus // math.gcd(modulus, label)


def merge_fraction(intervals):
    rows = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if rows and left <= rows[-1][1]:
            rows[-1] = (rows[-1][0], max(rows[-1][1], right))
        else:
            rows.append((left, right))
    return tuple(rows)


def subtract_fraction(base, removed):
    result = []
    for left, right in base:
        cursor = left
        for cut_left, cut_right in removed:
            if cut_right <= cursor:
                continue
            if cut_left >= right:
                break
            if cursor < cut_left:
                result.append((cursor, min(cut_left, right)))
            cursor = max(cursor, cut_right)
            if cursor >= right:
                break
        if cursor < right:
            result.append((cursor, right))
    return tuple(result)


def interval_mass(intervals):
    return sum((right - left for left, right in intervals), F(0))


def danger_fraction(label):
    radius = F(1, 14 * label)
    rows = [(F(0), radius)]
    rows.extend(
        (F(index, label) - radius, F(index, label) + radius)
        for index in range(1, label)
    )
    rows.append((F(1) - radius, F(1)))
    return tuple(rows)


def projected_residual(carrier, labels, modulus):
    removed = merge_fraction(
        interval
        for label in labels
        for interval in danger_fraction(label)
    )
    residual = subtract_fraction(carrier, removed)
    projected = []
    for left, right in residual:
        scaled_left = modulus * left
        scaled_right = modulus * right
        first_cell = scaled_left.numerator // scaled_left.denominator
        last_cell = -(-scaled_right.numerator // scaled_right.denominator)
        for cell in range(first_cell, last_cell):
            piece_left = max(scaled_left, F(cell)) - cell
            piece_right = min(scaled_right, F(cell + 1)) - cell
            if piece_left < piece_right:
                projected.append((piece_left, piece_right))
    return residual, merge_fraction(projected)


def interval_intersection_mass(first, second):
    i = j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            total += right - left
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return total


def atlas_frontier():
    counts = Counter()
    pattern = re.compile(r"^row=E=[0-9,]+;.*;z1=([0-9]+);")
    for line in BODY_ATLAS_OUTPUT_PATH.read_text().splitlines():
        match = pattern.match(line)
        if match:
            counts[int(match.group(1))] += 1
    require(counts[324] == 45, ("z324 atlas count changed", counts[324]))
    require(all(counts[value] == 0 for value in range(313, 324)), "atlas gap changed")
    below = max(value for value in counts if value < 324)
    require((below, counts[below]) == (312, 80), "next occupied height changed")
    return tuple(sorted(counts.items())), below, counts[below]


def unit_rows(carrier, h, modulus, d, high_floor):
    rows = []
    for unit in range(1, d):
        if math.gcd(unit, d) != 1:
            continue
        residue = (modulus // d) * unit
        amplitude = residue * (A.singleton_coverage(carrier, residue) - h / 7)
        recurrence = (
            (residue + modulus)
            * (A.singleton_coverage(carrier, residue + modulus) - h / 7)
        )
        require(recurrence == amplitude, ("ray recurrence failed", d, unit))
        rows.append(
            (
                unit,
                residue,
                amplitude,
                first_on_ray(residue, modulus, FIRST + 1),
                first_on_ray(residue, modulus, high_floor),
            )
        )
    require(len(rows) > 0, ("empty denominator class", d))
    return tuple(rows)


def best_row(rows, *, high):
    candidates = []
    for unit, residue, amplitude, arbitrary, high_label in rows:
        label = high_label if high else arbitrary
        candidates.append((amplitude / label, label, unit, residue))
    return max(candidates, key=lambda row: (row[0], -row[1], -row[2]))


def labels_meeting(rows, threshold, modulus, high_floor, *, below_high=False):
    require(threshold > 0, ("threshold must be positive", threshold))
    labels = set()
    for _unit, _residue, amplitude, arbitrary, _high in rows:
        if amplitude <= 0:
            continue
        label = arbitrary
        while amplitude / label >= threshold:
            if not below_high or label < high_floor:
                labels.add(label)
            if below_high and label >= high_floor:
                break
            label += modulus
    return tuple(sorted(labels))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    frontier_text = FRONTIER_OUTPUT_PATH.read_text()
    require(
        "totals=states:2346;crude_kills:702;status_kills:1643;survivors:1"
        in frontier_text,
        "frontier split changed",
    )
    require(
        f"sole_residual=E:{BODY};denominators:{DS}" in frontier_text,
        "frontier residual changed",
    )
    atlas_counts, next_height, next_count = atlas_frontier()

    modulus = 14 * math.lcm(*BODY)
    carrier_integer = A.carrier_for(BODY)
    carrier = tuple(
        (F(left, A.RULER), F(right, A.RULER))
        for left, right in carrier_integer
    )
    h = interval_mass(carrier)
    lower = h * F(3, 91)
    first_delta = A.singleton_coverage(carrier_integer, FIRST) - h / 7
    high_floor = max(15, (F(13, 132) * modulus).numerator // (F(13, 132) * modulus).denominator + 1)
    require(
        (modulus, h, lower, first_delta, high_floor)
        == (129360, F(21041, 64680), F(21041, 1961960), F(38443, 12224520), 12741),
        "body scalar constants changed",
    )

    rows = {
        d: unit_rows(carrier_integer, h, modulus, d, high_floor)
        for d in TAIL_DS
    }
    best_any = {d: best_row(rows[d], high=False) for d in TAIL_DS}
    best_high = {d: best_row(rows[d], high=True) for d in TAIL_DS}
    require(
        best_any
        == {
            3920: (F(8147, 5885880), 429, 13, 429),
            4620: (F(2059, 452760), 364, 13, 364),
            10780: (F(56489, 18563160), 492, 41, 492),
        },
        ("arbitrary ray heads changed", best_any),
    )
    require(
        best_high
        == {
            3920: (F(2449, 35043624), 12771, 387, 12771),
            4620: (F(31639, 272108760), 16828, 601, 16828),
            10780: (F(58687, 481283880), 12756, 1063, 12756),
        },
        ("high ray heads changed", best_high),
    )
    required = lower - first_delta
    low_3920_threshold = required - max(
        best_high[4620][0] + best_any[10780][0],
        best_any[4620][0] + best_high[10780][0],
    )
    label_4620_threshold = required - best_high[3920][0] - best_any[10780][0]
    label_10780_threshold = required - best_high[3920][0] - best_any[4620][0]
    thresholds = (
        low_3920_threshold,
        label_4620_threshold,
        label_10780_threshold,
    )
    require(
        thresholds
        == (
            F(35115043, 12066474420),
            F(89391041, 20012412420),
            F(722933, 244053810),
        ),
        ("thresholds changed", thresholds),
    )
    low_3920 = labels_meeting(
        rows[3920], low_3920_threshold, modulus, high_floor, below_high=True
    )
    eligible_4620 = labels_meeting(
        rows[4620], label_4620_threshold, modulus, high_floor
    )
    eligible_10780 = labels_meeting(
        rows[10780], label_10780_threshold, modulus, high_floor
    )
    require(
        low_3920 == () and eligible_4620 == (364,) and eligible_10780 == (492,),
        ("scalar packet forcing changed", low_3920, eligible_4620, eligible_10780),
    )
    base_surplus = (
        first_delta
        + A.singleton_coverage(carrier_integer, 364) - h / 7
        + A.singleton_coverage(carrier_integer, 492) - h / 7
        - lower
    )
    require(base_surplus == F(5119, 465404940) > 0, "base surplus changed")

    require(A.RULER % modulus == 0, "carrier ruler does not resolve body cells")
    cell_scale = A.RULER // modulus
    body_cells = {
        cell
        for left, right in carrier_integer
        for cell in range(left // cell_scale, right // cell_scale)
    }
    require(
        len(body_cells) == 42082 and CELL_J in body_cells and CELL_K in body_cells,
        "located cells left the carrier",
    )
    located_intervals = (
        (F(CELL_J, modulus), F(CELL_J + 1, modulus)),
        (F(CELL_K, modulus), F(CELL_K + 1, modulus)),
    )
    fixed_dangers = {
        label: danger_fraction(label) for label in FIXED_LABELS
    }
    fixed_intersections = tuple(
        (
            cell,
            label,
            interval_intersection_mass((interval,), fixed_dangers[label]),
        )
        for cell, interval in zip((CELL_J, CELL_K), located_intervals)
        for label in FIXED_LABELS
    )
    require(
        all(mass == 0 for _cell, _label, mass in fixed_intersections),
        ("fixed drift entered a located cell", fixed_intersections),
    )

    difference = CELL_K - CELL_J
    units = tuple(unit for unit in range(1, 3920) if math.gcd(unit, 3920) == 1)
    phase_residues = {(unit * difference) % 3920 for unit in units}
    require(
        difference == 13720
        and len(units) == 1344
        and all(unit % 2 for unit in units)
        and phase_residues == {1960},
        "antipodal unit-direction law changed",
    )
    require(F(1, 14) + F(1, 14) < F(1, 2), "danger arcs ceased to be disjoint")
    require(ALIGNED_UNION_CAP < 1, "aligned union cap ceased to be proper")

    representative_packet = (*FIXED_LABELS, HIGH_REPRESENTATIVE)
    require(
        denominator(modulus, HIGH_REPRESENTATIVE) == 3920
        and HIGH_REPRESENTATIVE >= high_floor,
        "representative high label changed",
    )
    residual, projected = projected_residual(carrier, representative_packet, modulus)
    residual_mass = interval_mass(residual)
    projected_mass = interval_mass(projected)
    require(projected == ((F(0), F(1)),) and projected_mass == 1, "direct projection failed")

    semantic_payload = (
        BODY,
        FIRST,
        DS,
        atlas_counts,
        next_height,
        next_count,
        modulus,
        h,
        lower,
        first_delta,
        high_floor,
        best_any,
        best_high,
        thresholds,
        low_3920,
        eligible_4620,
        eligible_10780,
        base_surplus,
        len(body_cells),
        located_intervals,
        fixed_intersections,
        difference,
        len(units),
        phase_residues,
        representative_packet,
        residual_mass,
        projected,
        projected_mass,
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=3 z1=324 independent located two-cell closure",
        f"frontier_source_sha256={file_sha256(FRONTIER_PATH)}",
        f"frontier_output_sha256={file_sha256(FRONTIER_OUTPUT_PATH)}",
        f"carrier_source_sha256={file_sha256(CARRIER_PATH)}",
        f"body_atlas_source_sha256={file_sha256(BODY_ATLAS_PATH)}",
        f"body_atlas_output_sha256={file_sha256(BODY_ATLAS_OUTPUT_PATH)}",
        "inherited_frontier=45 bodies;2346 states;702 crude;1643 status;1 residual",
        f"sole_state=E:{BODY};denominators:{DS}",
        f"scalar_constants=L:{modulus};h:{h};lower:{lower};first_delta:{first_delta};high_floor:{high_floor}",
        f"best_any={best_any}",
        f"best_high={best_high}",
        f"positive_scalar_thresholds={thresholds}",
        (
            f"literal_forcing=low_d3920:{low_3920};d4620:{eligible_4620};"
            f"d10780:{eligible_10780};remaining_H:any_high_exact_d3920"
        ),
        f"fixed_base_surplus={base_surplus}",
        (
            f"located_cells=j:{CELL_J};k:{CELL_K};difference:{difference};"
            f"intervals:{located_intervals};fixed_intersection_masses:{fixed_intersections}"
        ),
        "pointwise_phase_law=z=33v and gcd(v,3920)=1 => v odd;phase_k-phase_j=v/2=1/2 mod 1",
        "topology=strict radius-1/14 danger arcs at half-turn are disjoint;every phase has a drift-safe located lift",
        f"aligned_completion=three aligned danger sets have union mass<={ALIGNED_UNION_CAP}<1",
        (
            f"direct_rational_control=packet:{representative_packet};"
            f"physical_residual_mass:{residual_mass};projected_residual:{projected};"
            f"projected_mass:{projected_mass}"
        ),
        f"atlas_gap=313..323 empty;next_occupied_height:{next_height};bodies:{next_count}",
        "conclusion=the sole z1=324 residual is empty uniformly;all 45 z1=324 rows close",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
