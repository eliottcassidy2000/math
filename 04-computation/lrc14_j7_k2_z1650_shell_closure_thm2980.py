#!/usr/bin/env python3
"""Exact projected k=2 descent through the twelve-row z1=1650 shell.

This finite referee first runs the complete scalar atlas on 1600..1655.  The
top occupied height is 1650, with exactly twelve rows and an empty gap
1651..1655.  Three ordinary rows close by exact common-status/projected
certificates.  Eight forced-high rows are scalar-empty.  The ninth forced-high
row has 167 literal denominator classes, all sharing one low triple; located
torsion closes 166 classes and an exact projected-cell prefix closes d=143.

The consequence is the projected k=2 bound z1<=1649.  This is a finite-exact
boundary input for reserved THM-2980, not the full 1600..1679 theorem and not
a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, gcd, lcm
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
REF_NAME = "lrc14_j7_k2_z1656_boundary_closure_thm2980"
REF_SOURCE = ROOT / f"04-computation/{REF_NAME}.py"
REF_OUTPUT = ROOT / f"05-knowledge/results/{REF_NAME}.out"
BASE_NAME = "lrc14_j7_k2_z1660_status_descent_closure_thm2941"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1650_shell_closure_thm2980.out"

EXPECTED_REF_SOURCE_SHA256 = "593988b1963fc9f8406d7f7d0ae0f5cc0e9a8bd2604a5554522787fd2d0fa6b3"
EXPECTED_REF_OUTPUT_SHA256 = "c93f0fd172bb3a2a2148686566d4645fdd2db2a73e77c5048ba302a7c1bb1207"
EXPECTED_BASE_SOURCE_SHA256 = "b7d3344888034c1abc4cc5da72b73c8a81992e287ba7fc843b35cfbc6d5b305a"
EXPECTED_BASE_OUTPUT_SHA256 = "1a1a85083ebb382d478139ce5206a20131ed412c5e47c61b26b58975cb2fed2e"
EXPECTED_ATLAS_PROFILE_SHA256 = "218d46b572fce7696f4d878fae03e7248e24c2c6e06e75d2fec06682ec2a6949"
EXPECTED_ATLAS_SURVIVOR_SHA256 = "2b567b151dbc6b616b01ac71cb7b94877631ef19e16facab4bf065ba84a4d701"
EXPECTED_CLOSURE_PROFILE_SHA256 = "c60de012a865086415c757741e7d473ab977fb38092d1906ab72140f09187a03"
EXPECTED_HARD_CLASS_SHA256 = "5bef285d533d0fd37e141ccb478359ba6550ac486afcf6f2c993440c3e8bc1f6"
EXPECTED_SEMANTIC_SHA256 = "a38d0d4ae78bbf0b8a31411c32861390ce3d60cf52237152d01c2f038e66998f"

START = 1600
END = 1655
EXPECTED_CANDIDATE_ROWS = 168_168
EXPECTED_HEIGHTS = (
    (1604, 1), (1608, 2), (1612, 19), (1616, 1), (1618, 1),
    (1620, 8), (1625, 1), (1630, 1), (1631, 1), (1642, 1),
    (1645, 1), (1650, 12),
)
CASES = tuple((1650, body) for body in (
    (1, 2, 10, 11, 12, 14),
    (1, 2, 10, 12, 13, 14),
    (1, 4, 8, 10, 12, 14),
    (1, 4, 10, 11, 12, 14),
    (1, 4, 10, 12, 13, 14),
    (1, 8, 10, 12, 13, 14),
    (1, 10, 11, 12, 13, 14),
    (2, 3, 10, 11, 12, 14),
    (2, 4, 6, 8, 10, 14),
    (2, 4, 8, 10, 12, 14),
    (2, 4, 10, 11, 12, 14),
    (2, 8, 10, 11, 12, 14),
))
HARD_CASE = (1650, (1, 10, 11, 12, 13, 14))
ORDINARY_CASES = tuple(case for case in CASES if 14 * lcm(*case[1]) == 11760)
EMPTY_DIGEST = sha256(b"()").hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(REF_SOURCE) == EXPECTED_REF_SOURCE_SHA256, "z1656 referee source changed")
require(file_sha256(REF_OUTPUT) == EXPECTED_REF_OUTPUT_SHA256, "z1656 referee output changed")
require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1660 atlas source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1660 atlas output changed")
SPEC = spec_from_file_location("k2_z1650_referee_base", REF_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1656 referee")
REF = module_from_spec(SPEC)
SPEC.loader.exec_module(REF)
ENGINE = REF.ENGINE
TORSION = REF.BASE.BASE.TORSION
EXACT_STATUS_AUDIT = REF.EXACT_STATUS_AUDIT


def configure_atlas():
    bank = REF.BASE.BASE.BASE.ATLAS.B
    bank.START = START
    bank.END = END
    bank.LABELS = np.arange(START, bank.HORIZON + 1, dtype=np.int64)


def atlas_profile(body):
    return REF.BASE.BASE.BASE.ATLAS.B.profile(body)


def ordinary_profile(case):
    first, body = case
    carrier = ENGINE.U.suffix.A.carrier_for(body)
    require(ENGINE.P.A.carrier_for(body) == carrier, (case, "carrier engines disagree"))
    h = F(sum(right - left for left, right in carrier), ENGINE.U.suffix.A.RULER)
    lower = h * ENGINE.U.suffix.ETAS[2]
    modulus = 14 * lcm(*body)
    result = ENGINE.ray_and_status(first, body, carrier, h, lower, modulus)
    counts = tuple(len(result[index]) for index in (6, 7, 8, 9))
    certificate_count = EXACT_STATUS_AUDIT(case, modulus, result[8], result[10][2])
    projected = (
        ENGINE.projected_packets(
            first, body, carrier, h, lower, modulus,
            result[0], result[4], result[9],
        )
        if result[9]
        else (0, 0, 0, None, 0, None, None, EMPTY_DIGEST, ())
    )
    require(certificate_count == counts[2], (case, "certificate ledger changed"))
    require(projected[1] == projected[2], (case, "projected packet survived"))
    return (
        "ALL-LABEL", first, body, h, len(carrier), modulus, lower,
        result[4], result[5], result[1], result[2], result[3], counts,
        result[10], certificate_count, *projected[:8],
    )


def high_profile(case):
    first, body = case
    modulus = 14 * lcm(*body)
    carrier, mass_numerator, amplitudes = TORSION.amplitude_data(body, modulus)
    h = F(mass_numerator, TORSION.RULER)
    first_delta = TORSION.delta(amplitudes, first)
    lower = h / 91
    high_floor = (13 * modulus) // 150 + 1
    low_rows = sorted(
        ((TORSION.delta(amplitudes, label), label) for label in range(first + 1, high_floor)),
        key=lambda item: (-item[0], item[1]),
    )
    required = lower - first_delta
    maxima = TORSION.high_ray_maxima(amplitudes, modulus, high_floor)
    best_d, best_high = max(maxima.items(), key=lambda item: item[1][0])
    top3 = tuple(low_rows[:3])
    two_gap = required - (
        2 * best_high[0]
        + max(best_high[0], top3[0][0])
        + max(best_high[0], top3[1][0])
    )
    exact_gap = required - sum((value for value, _label in top3), F()) - best_high[0]
    require(two_gap > 0, (case, "two-high scalar branch survived"))
    return (
        "FORCED-HIGH", first, body, h, len(carrier), modulus, high_floor,
        first_delta, lower, required, two_gap, exact_gap,
        tuple(label for _value, label in top3), best_d, best_high[0], best_high[1],
        sha256(repr(tuple(int(value) for value in amplitudes)).encode()).hexdigest(),
    )


def closure_profile(case):
    return ordinary_profile(case) if case in ORDINARY_CASES else high_profile(case)


def hard_class_audit():
    first, body = HARD_CASE
    modulus = 14 * lcm(*body)
    carrier, mass_numerator, amplitudes = TORSION.amplitude_data(body, modulus)
    h = F(mass_numerator, TORSION.RULER)
    lower = h / 91
    first_delta = TORSION.delta(amplitudes, first)
    required = lower - first_delta
    high_floor = (13 * modulus) // 150 + 1
    low_rows = sorted(
        ((TORSION.delta(amplitudes, label), label) for label in range(first + 1, high_floor)),
        key=lambda item: (-item[0], item[1]),
    )
    maxima = TORSION.high_ray_maxima(amplitudes, modulus, high_floor)
    fixed_labels = (first, 1736, 1800, 2340)
    cells = TORSION.complete_body_cells(carrier, modulus)
    safe = TORSION.fixed_safe_cells(cells, modulus, fixed_labels)
    classes = []
    torsion_failures = []
    order_histogram = Counter()
    for denominator, high in sorted(maxima.items()):
        high_value, high_label, high_residue, high_amplitude = high
        if high_label is None:
            continue
        for low_labels, low_total in TORSION.finite_triples(low_rows, required - high_value):
            require(low_labels == fixed_labels[1:], (denominator, "second low triple appeared", low_labels))
            unit = high_residue // (modulus // denominator)
            require(high_amplitude > 0 and gcd(unit, denominator) == 1, (denominator, "high ray typing"))
            torsions = TORSION.torsion_pairs(safe, denominator)
            if not torsions:
                torsion_failures.append((denominator, high_label, unit))
            for order, shift, first_cell, second_cell in torsions:
                require(
                    (second_cell - first_cell) % denominator == shift
                    and order * shift == denominator
                    and gcd(unit, order) == 1,
                    (denominator, unit, "located torsion typing"),
                )
                order_histogram[order] += 1
            excess = low_total + high_value - required
            classes.append((denominator, high_label, high_residue, unit, excess, torsions))
    classes = tuple(classes)
    require(len(classes) == 167, "hard literal class count changed")
    require(torsion_failures == [(143, 105840, 18)], ("torsion failure ledger changed", torsion_failures))
    require(min(row[4] for row in classes) == F(83, 10377542175), "minimum hard excess changed")
    require(len(cells) == 256232 and len(safe) == 135766, "hard fixed-safe cell ledger changed")
    require(
        order_histogram == Counter({2: 128, 3: 86, 4: 88, 5: 87, 6: 65, 7: 128}),
        ("hard torsion histogram changed", order_histogram),
    )
    class_hash = sha256(repr(classes).encode()).hexdigest()
    require(class_hash == EXPECTED_HARD_CLASS_SHA256, "hard class digest changed")

    low_total = sum((TORSION.delta(amplitudes, label) for label in fixed_labels[1:]), F())
    high_threshold = required - low_total
    eligible_high = []
    denominator = 143
    step = modulus // denominator
    for unit in range(1, denominator):
        if gcd(unit, denominator) != 1:
            continue
        residue = step * unit
        amplitude = int(amplitudes[residue - 1])
        if amplitude <= 0:
            continue
        label = residue if residue >= high_floor else residue + modulus
        while F(amplitude, 7 * TORSION.RULER * label) >= high_threshold:
            value = F(amplitude, 7 * TORSION.RULER * label)
            eligible_high.append((value, label, residue, unit, value - high_threshold))
            label += modulus
    eligible_high.sort(key=lambda row: (-row[0], row[1:]))
    require(
        tuple(eligible_high) == (
            (F(1, 245245), 105840, 105840, 18, F(3754, 1482506025)),
            (F(97, 58858800), 117600, 117600, 20, F(487, 4744019280)),
            (F(97, 58858800), 235200, 235200, 40, F(487, 4744019280)),
        ),
        ("d=143 eligible ray ledger changed", eligible_high),
    )
    projected_cells = ENGINE.P.body_cells(ENGINE.P.A.carrier_for(body), modulus)
    projected_packets = []
    for high_value, high_label, _residue, _unit, excess in eligible_high:
        hostile_labels = (*fixed_labels, high_label)
        projected_lower, cells_used, common = ENGINE.projected_safe_lower_bound(
            projected_cells, modulus, hostile_labels
        )
        projected_margin = projected_lower - REF.ALIGNED_TWO_UNION_CAP
        require(
            len(projected_cells) == len(cells)
            and cells_used == 146
            and projected_lower == F(2, 5)
            and projected_margin == F(57, 455)
            and len(common) == 1,
            ("d=143 projected hostile changed", high_label),
        )
        projected_packets.append(
            (hostile_labels, high_value, excess, cells_used, projected_lower, projected_margin)
        )
    return (
        fixed_labels, len(cells), len(safe), len(classes),
        min(row[4] for row in classes), tuple(sorted(order_histogram.items())),
        tuple(torsion_failures), tuple(projected_packets), class_hash,
    )


def render(atlas_profiles, closure_profiles, hard_audit):
    require(len(atlas_profiles) == comb(14, 6) == 3003, "atlas body universe changed")
    candidate_rows = sum(row[7] for row in atlas_profiles)
    survivors = tuple(sorted(
        (survivor for row in atlas_profiles for survivor in row[10]),
        key=lambda row: (row[1], row[0]),
    ))
    heights = tuple(sorted(Counter(row[1] for row in survivors).items()))
    keys = tuple((row[1], row[0]) for row in survivors)
    top_keys = tuple(key for key in keys if key[0] == max(first for first, _body in keys))
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "candidate ledger changed")
    require(heights == EXPECTED_HEIGHTS, ("height ledger changed", heights))
    require(len(survivors) == 49 and top_keys == CASES, ("z1650 universe changed", top_keys))
    require(not any(1651 <= first <= 1655 for first, _body in keys), "descent gap occupied")
    atlas_hash = sha256(repr(tuple(atlas_profiles)).encode()).hexdigest()
    survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    require(atlas_hash == EXPECTED_ATLAS_PROFILE_SHA256, "atlas profile digest changed")
    require(survivor_hash == EXPECTED_ATLAS_SURVIVOR_SHA256, "atlas survivor digest changed")

    closure_hash = sha256(repr(tuple(closure_profiles)).encode()).hexdigest()
    require(closure_hash == EXPECTED_CLOSURE_PROFILE_SHA256, "closure profile digest changed")
    ordinary = tuple(row for row in closure_profiles if row[0] == "ALL-LABEL")
    high = tuple(row for row in closure_profiles if row[0] == "FORCED-HIGH")
    require((len(ordinary), len(high)) == (3, 9), "closure partition changed")
    ordinary_totals = tuple(sum(row[12][index] for row in ordinary) for index in range(4))
    require(ordinary_totals == (593, 71, 365, 157), ("ordinary status ledger changed", ordinary_totals))
    require(sum(row[14] for row in ordinary) == 365, "certificate count changed")
    require(sum(row[16] for row in ordinary) == sum(row[17] for row in ordinary) == 329, "packet ledger changed")
    scalar_empty_high = tuple(row for row in high if row[11] > 0)
    scalar_live_high = tuple(row for row in high if row[11] <= 0)
    require(len(scalar_empty_high) == 8, "forced-high scalar-empty split changed")
    require(tuple((row[1], row[2]) for row in scalar_live_high) == (HARD_CASE,), "hard row changed")
    require(scalar_live_high[0][11] == F(-29147179, 1081241060900), "hard scalar gap changed")

    semantic_payload = (
        EXPECTED_REF_SOURCE_SHA256, EXPECTED_REF_OUTPUT_SHA256,
        EXPECTED_BASE_SOURCE_SHA256, EXPECTED_BASE_OUTPUT_SHA256,
        START, END, candidate_rows, heights, top_keys, ordinary_totals,
        365, 329, hard_audit, 1649, atlas_hash, survivor_hash, closure_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact z1=1650 shell closure",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"z1656_source_sha256={file_sha256(REF_SOURCE)}",
        f"z1656_output_sha256={file_sha256(REF_OUTPUT)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        f"atlas_band={START}..{END};candidate_rows={candidate_rows};survivors={len(survivors)};height_counts={heights}",
        "descent_gap=1651..1655 empty;next_occupied_height=1650;z1650_rows=12",
        "closure_partition=ordinary_all_label:3;forced_high_scalar_empty:8;forced_high_literal:1;empty:12;survivors:0",
        f"ordinary_ledger=scalar:{ordinary_totals[0]};crude:{ordinary_totals[1]};status:{ordinary_totals[2]};status_survivors:{ordinary_totals[3]};exact_status_certificates:365;literal_packets:329;projected_kills:329",
    ]
    for survivor in survivors:
        body, first, delta1, suffix, upper, lower, gap, h, carrier_count, modulus, high_floor = survivor
        route = "HIGH-TAIL" if any(kind == "HIGH-TAIL" for _value, _label, kind in suffix) else "EXACT"
        lines.append(
            f"ATLAS;z1={first};E={','.join(map(str, body))};route={route};h={ftext(h)};"
            f"r={carrier_count};L={modulus};high_floor={high_floor};delta1={ftext(delta1)};"
            f"lower={ftext(lower)};upper={ftext(upper)};gap={ftext(gap)};suffix={suffix}"
        )
    for row in closure_profiles:
        if row[0] == "FORCED-HIGH":
            conclusion = "SCALAR-EMPTY" if row[11] > 0 else "LITERAL-CLOSURE-BELOW"
            lines.append(
                f"CASE;route=FORCED-HIGH;z1={row[1]};E={','.join(map(str,row[2]))};"
                f"h={ftext(row[3])};r={row[4]};L={row[5]};high_floor={row[6]};"
                f"delta1={ftext(row[7])};lower={ftext(row[8])};two_high_gap={ftext(row[10])};"
                f"exact_one_high_gap={ftext(row[11])};top_lows={row[12]};best_high_d={row[13]};"
                f"best_high_value={ftext(row[14])};best_high_label={row[15]};"
                f"amplitude_sha256={row[16]};conclusion={conclusion}"
            )
        else:
            lines.append(
                f"CASE;route=ALL-LABEL;z1={row[1]};E={','.join(map(str,row[2]))};"
                f"h={ftext(row[3])};r={row[4]};L={row[5]};lower={ftext(row[6])};"
                f"delta1={ftext(row[7])};first_d={row[8]};ray_sha256={row[9]};"
                f"divisors={row[10]};trials={row[11]};counts={row[12]};"
                f"stage_sha256={row[13]};exact_status_certificates={row[14]};"
                f"cells={row[15]};packets={row[16]};kills={row[17]};"
                f"minimum_margin={ftext(row[18])};maximum_prefix_cells={row[19]};"
                f"direct_mass={ftext(row[21])};state_sha256={row[22]};conclusion=EMPTY"
            )
    lines.extend((
        f"HARD;z1=1650;E=1,10,11,12,13,14;fixed_labels={hard_audit[0]};complete_cells={hard_audit[1]};fixed_safe_cells={hard_audit[2]};literal_denominator_classes={hard_audit[3]};minimum_excess={ftext(hard_audit[4])};torsion_order_histogram={hard_audit[5]};torsion_failures={hard_audit[6]};class_sha256={hard_audit[8]}",
        f"HARD-PROJECTED;d=143;eligible_literal_packets={hard_audit[7]};conclusion=all-three-EMPTY",
        "hostile_controls=the unrestricted one-high scalar gap is positive for the hard row;d143 is the unique located-torsion failure and requires projected cells",
        "consequence=projected k=2 first drift label z1<=1649",
        "scope_boundary=finite-exact THM-2980 input only;full reserved theorem and LRC(14) remain open",
        f"atlas_profile_sha256={atlas_hash}",
        f"atlas_survivor_sha256={survivor_hash}",
        f"closure_profile_sha256={closure_hash}",
        f"hard_class_sha256={hard_audit[8]}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(4, mp.cpu_count() or 1))
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    os.environ["PYTHONHASHSEED"] = "0"
    configure_atlas()
    bodies = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        atlas_profiles = tuple(atlas_profile(body) for body in bodies)
        closure_profiles = tuple(closure_profile(case) for case in CASES)
    else:
        with mp.get_context("spawn").Pool(args.workers, initializer=configure_atlas) as pool:
            atlas_profiles = tuple(pool.imap(atlas_profile, bodies, chunksize=4))
        with mp.get_context("spawn").Pool(args.workers) as pool:
            closure_profiles = tuple(pool.imap(closure_profile, CASES))
    payload = render(atlas_profiles, closure_profiles, hard_class_audit())
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
