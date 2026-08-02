#!/usr/bin/env python3
"""Exact closure of all eight projected k=2 rows at z1=1620.

Seven exact-suffix rows close by common status and projected packets.  The
single genuine high-tail row has a surviving one-high scalar invoice, but the
THM-2984 complete-cell section gate closes all 147 exact denominator/triple
classes.  Located torsion closes the same 147 classes independently.
This finite-exact boundary input gives z1<=1619; THM-2980 remains reserved.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
BASE_NAME = "lrc14_j7_k2_z1625_forced_high_closure_thm2980"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1620_mixed_shell_closure_thm2980.out"

EXPECTED_BASE_SOURCE_SHA256 = "ec98feabf60fe138df56647fc436ddec54b3f36b0704043c46af21e648e43935"
EXPECTED_BASE_OUTPUT_SHA256 = "85947083c3bfdaeeaafd02aee9dac4a651b16c97fd80c1b84759d9380217fe7a"
EXPECTED_PROFILE_SHA256 = "a0ba343ada8e718e5a277954578c189f263a571bd8a5dbf905ec5e087f5b3df1"
EXPECTED_CLASS_SHA256 = "aab1b5f7f99717e6afbb163350d80a24b0927fe630370223f244ddc12eb9ec78"
EXPECTED_SECTION_SHA256 = "a96f2b6df34e40f9bceafa7fb1ffe2c0a57ca67bcb7f05214588341f4dae056f"
EXPECTED_SEMANTIC_SHA256 = "46cd5ec80a0d214a218d9cded1f42147f32f5ee92d1e206f36d45332352089cb"

HIGH_CASE = (1620, (2, 8, 9, 10, 12, 14))
CASES = tuple((1620, body) for body in (
    (1, 2, 8, 10, 12, 14),
    (1, 4, 8, 10, 12, 14),
    (1, 4, 9, 10, 12, 14),
    (2, 4, 8, 10, 12, 14),
    (2, 4, 9, 10, 12, 14),
    (2, 8, 9, 10, 12, 14),
    (3, 4, 8, 10, 12, 14),
    (4, 8, 9, 10, 12, 14),
))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1625 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1625 output changed")
require(
    "consequence=projected k=2 first drift label z1<=1624" in BASE_OUTPUT.read_text(encoding="utf-8"),
    "z1625 handoff changed",
)
SPEC = spec_from_file_location("k2_z1620_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1625 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
ENGINE = BASE.BASE.BASE.BASE
TORSION = ENGINE.TORSION


def profile(case):
    return ENGINE.high_profile(case) if case == HIGH_CASE else ENGINE.ordinary_profile(case)


def high_class_audit():
    first, body = HIGH_CASE
    modulus = 14 * lcm(*body)
    carrier, mass_numerator, amplitudes = TORSION.amplitude_data(body, modulus)
    lower = F(mass_numerator, TORSION.RULER) / 91
    required = lower - TORSION.delta(amplitudes, first)
    high_floor = (13 * modulus) // 150 + 1
    low_rows = sorted(
        ((TORSION.delta(amplitudes, label), label) for label in range(first + 1, high_floor)),
        key=lambda item: (-item[0], item[1]),
    )
    maxima = TORSION.high_ray_maxima(amplitudes, modulus, high_floor)
    cells = TORSION.complete_body_cells(carrier, modulus)
    safe_cache = {}
    classes = []
    failures = []
    order_histogram = Counter()
    active_denominators = Counter()
    section_rows = []
    for denominator, high in sorted(maxima.items()):
        high_value, high_label, high_residue, high_amplitude = high
        if high_label is None:
            continue
        unit = high_residue // (modulus // denominator)
        require(high_amplitude > 0 and gcd(unit, denominator) == 1, "high ray typing changed")
        for low_labels, low_total in TORSION.finite_triples(low_rows, required - high_value):
            active_denominators[denominator] += 1
            if low_labels not in safe_cache:
                safe_cache[low_labels] = TORSION.fixed_safe_cells(
                    cells, modulus, (first, *low_labels)
                )
            safe = safe_cache[low_labels]
            residue_count = len({cell % denominator for cell in safe})
            short_order = max(order for order in range(1, 8) if denominator % order == 0)
            section_threshold = denominator // short_order
            section_surplus = residue_count - section_threshold
            require(section_surplus > 0, (denominator, low_labels, "complete-cell section gate failed"))
            section_rows.append(
                (denominator, low_labels, residue_count, short_order, section_threshold, section_surplus)
            )
            torsions = TORSION.torsion_pairs(safe, denominator)
            if not torsions:
                failures.append((denominator, low_labels, high_label, unit))
            for order, shift, first_cell, second_cell in torsions:
                require(
                    (second_cell - first_cell) % denominator == shift
                    and order * shift == denominator
                    and gcd(unit, order) == 1,
                    (denominator, unit, "located torsion typing"),
                )
                order_histogram[order] += 1
            excess = low_total + high_value - required
            classes.append(
                (denominator, low_labels, high_label, high_residue, unit, excess, len(safe), torsions)
            )
    classes = tuple(classes)
    section_rows = tuple(section_rows)
    class_hash = sha256(repr(classes).encode()).hexdigest()
    section_hash = sha256(repr(section_rows).encode()).hexdigest()
    require(len(cells) == 11958 and len(classes) == 147 and len(safe_cache) == 59, "class universe changed")
    require(not failures, ("located torsion failure", failures))
    require(len(section_rows) == 147 and min(row[5] for row in section_rows) == 90, "section ledger changed")
    require(min(row[5] for row in classes) == F(464015191, 1855060847738280), "minimum excess changed")
    require((min(map(len, safe_cache.values())), max(map(len, safe_cache.values()))) == (4780, 5236), "safe-cell range changed")
    require(
        order_histogram == Counter({2: 143, 3: 139, 4: 138, 5: 49, 6: 135, 7: 147}),
        ("torsion histogram changed", order_histogram),
    )
    require(tuple(active_denominators.items()) == (
        (105, 1), (196, 3), (420, 9), (504, 23), (588, 59),
        (630, 5), (840, 5), (980, 5), (1176, 2), (1260, 5),
        (1764, 9), (2205, 3), (2520, 14), (2940, 2), (3528, 2),
    ), "active denominator ledger changed")
    require(class_hash == EXPECTED_CLASS_SHA256, "class digest changed")
    if EXPECTED_SECTION_SHA256 is not None:
        require(section_hash == EXPECTED_SECTION_SHA256, "section digest changed")
    return (
        len(cells), len(classes), len(safe_cache),
        (min(map(len, safe_cache.values())), max(map(len, safe_cache.values()))),
        tuple(active_denominators.items()), min(row[5] for row in classes),
        tuple(sorted(order_histogram.items())), min(row[5] for row in section_rows),
        section_hash, class_hash,
    )


def render(rows, hard):
    require(tuple((row[1], row[2]) for row in rows) == CASES, "z1620 universe changed")
    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    ordinary = tuple(row for row in rows if row[0] == "ALL-LABEL")
    high = tuple(row for row in rows if row[0] == "FORCED-HIGH")
    require((len(ordinary), len(high)) == (7, 1), "route partition changed")
    totals = tuple(sum(row[12][index] for row in ordinary) for index in range(4))
    require(totals == (1647, 281, 1202, 164), ("ordinary ledger changed", totals))
    require(sum(row[14] for row in ordinary) == 1202, "certificate ledger changed")
    require(sum(row[16] for row in ordinary) == sum(row[17] for row in ordinary) == 418, "packet ledger changed")
    require(high[0][10] > 0 and high[0][11] == F(-73815531457, 339041265184440) < 0, "high hostile changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256, EXPECTED_BASE_OUTPUT_SHA256, CASES,
        totals, 1202, 418, hard, 1619, profile_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 exact z1620 mixed-shell closure",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1624;atlas_gap=1621..1624 empty;boundary_rows=z1620:8",
        "closure_partition=all_label:7;forced_high_literal:1;empty:8;survivors:0",
        f"ordinary_ledger=scalar:{totals[0]};crude:{totals[1]};status:{totals[2]};status_survivors:{totals[3]};exact_status_certificates:1202;literal_packets:418;projected_kills:418",
    ]
    for row in rows:
        if row[0] == "FORCED-HIGH":
            lines.append(
                f"CASE;route=FORCED-HIGH;z1={row[1]};E={','.join(map(str,row[2]))};"
                f"h={ftext(row[3])};r={row[4]};L={row[5]};high_floor={row[6]};"
                f"delta1={ftext(row[7])};lower={ftext(row[8])};two_high_gap={ftext(row[10])};"
                f"exact_one_high_gap={ftext(row[11])};top_lows={row[12]};best_high_d={row[13]};"
                f"best_high_value={ftext(row[14])};best_high_label={row[15]};"
                f"amplitude_sha256={row[16]};conclusion=LITERAL-TORSION-EMPTY"
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
        f"HIGH-LITERAL;complete_cells={hard[0]};classes={hard[1]};low_triples={hard[2]};safe_cell_range={hard[3]};active_denominators={hard[4]};minimum_excess={ftext(hard[5])};section_gate=147/147;minimum_section_surplus={hard[7]};section_sha256={hard[8]};torsion_order_histogram={hard[6]};torsion_failures=0;class_sha256={hard[9]};conclusion=EMPTY",
        "hostile_controls=one-high scalar gap is negative;THM-2984 complete-cell sections close every class without pair selection, while located torsion independently closes the same 147 classes",
        "consequence=projected k=2 first drift label z1<=1619",
        "scope_boundary=finite-exact THM-2980 input only;full reserved theorem and LRC(14) remain open",
        f"profile_sha256={profile_hash}",
        f"section_sha256={hard[8]}",
        f"class_sha256={hard[9]}",
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
    if args.workers == 1:
        rows = tuple(profile(case) for case in CASES)
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = tuple(pool.imap(profile, CASES))
    payload = render(rows, high_class_audit())
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
