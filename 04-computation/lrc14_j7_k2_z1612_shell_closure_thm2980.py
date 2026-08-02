#!/usr/bin/env python3
"""Exact projected k=2 closure of all nineteen rows at z1=1612.

Eight exact-suffix rows close by common status and projected packets.  Seven
high-tail rows are scalar-empty.  For each of the four genuinely surviving
one-high rows, this referee retains complete fixed-safe cells and applies the
THM-2984 translated-band threshold |S mod d|>ceil(d/7) on every exact
denominator/low-triple class.  A cell-count gate is tried first and exact
residue counting is used only at its boundary.

This is a finite-exact boundary input for reserved THM-2980, not its proof.
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

import numpy as np


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
BASE_NAME = "lrc14_j7_k2_z1618_z1616_status_descent_thm2980"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1612_shell_closure_thm2980.out"

EXPECTED_BASE_SOURCE_SHA256 = "3d8a0a3daa4775826beeee7d98cd26c8f3483e410239274c4f6d41c4d82faf8a"
EXPECTED_BASE_OUTPUT_SHA256 = "20fdd1070210a4a61686dd13dc3bb1bee0c1ebc896b47400853593f44fa38fd2"
EXPECTED_PROFILE_SHA256 = "d83d42c00ba02afd7d7329be7318689d1d80a618a03a4f49aa90f3fb4cd59a56"
EXPECTED_SECTION_SHA256 = "8eb6544449075408dbc49ca46d7e614042786f8f053fab72adc8188157b7c51f"
EXPECTED_SEMANTIC_SHA256 = "1e3bf172d0d4844890456d72b1025c1aceb0ff00cba11f9f3010ae46ac82e971"

ORDINARY_BODIES = (
    (1, 2, 8, 10, 12, 14),
    (1, 4, 6, 8, 10, 14),
    (1, 4, 8, 10, 12, 13),
    (1, 4, 8, 10, 12, 14),
    (1, 6, 8, 10, 12, 14),
    (2, 4, 8, 10, 12, 14),
    (2, 6, 8, 10, 12, 14),
    (4, 6, 8, 10, 12, 14),
)
HIGH_BODIES = (
    (1, 2, 10, 12, 13, 14),
    (1, 3, 10, 12, 13, 14),
    (1, 4, 10, 12, 13, 14),
    (1, 6, 8, 10, 13, 14),
    (1, 6, 10, 12, 13, 14),
    (1, 8, 10, 12, 13, 14),
    (1, 10, 11, 12, 13, 14),
    (2, 3, 10, 12, 13, 14),
    (2, 4, 10, 12, 13, 14),
    (2, 8, 10, 12, 13, 14),
    (4, 8, 10, 12, 13, 14),
)
HARD_BODIES = (
    (1, 2, 10, 12, 13, 14),
    (1, 4, 10, 12, 13, 14),
    (1, 8, 10, 12, 13, 14),
    (1, 10, 11, 12, 13, 14),
)
CASES = tuple((1612, body) for body in sorted((*ORDINARY_BODIES, *HIGH_BODIES)))
HARD_CASES = tuple((1612, body) for body in HARD_BODIES)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1618/z1616 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1618/z1616 output changed")
require(
    "consequence=projected k=2 first drift label z1<=1615"
    in BASE_OUTPUT.read_text(encoding="utf-8"),
    "z1618/z1616 handoff changed",
)
SPEC = spec_from_file_location("k2_z1612_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1618/z1616 referee")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
ENGINE = BASE.ENGINE
TORSION = ENGINE.TORSION


def closure_profile(case):
    return ENGINE.ordinary_profile(case) if case[1] in ORDINARY_BODIES else ENGINE.high_profile(case)


def high_section_audit(case):
    first, body = case
    modulus = 14 * lcm(*body)
    carrier, mass_numerator, amplitudes = TORSION.amplitude_data(body, modulus)
    lower = F(mass_numerator, TORSION.RULER) / 91
    required = lower - TORSION.delta(amplitudes, first)
    high_floor = (13 * modulus) // 150 + 1
    lows = sorted(
        ((TORSION.delta(amplitudes, label), label) for label in range(first + 1, high_floor)),
        key=lambda item: (-item[0], item[1]),
    )
    maxima = TORSION.high_ray_maxima(amplitudes, modulus, high_floor)
    scalar_classes = []
    active_denominators = Counter()
    zero_supremum_classes = 0
    for denominator, high in sorted(maxima.items()):
        high_value, high_label, high_residue, high_amplitude = high
        if high_label is None:
            # Exact signed-ray attainment may make a zero/negative ray relevant
            # when the fixed low labels already reach the floor.  Retaining the
            # zero supremum is a safe finite over-approximation even at the
            # unattained negative equality boundary.
            unit = None
        else:
            unit = high_residue // (modulus // denominator)
            require(high_amplitude > 0 and gcd(unit, denominator) == 1, "high ray typing changed")
        for low_labels, low_total in TORSION.finite_triples(lows, required - high_value):
            excess = low_total + high_value - required
            require(excess >= 0, "finite triple below threshold")
            scalar_classes.append(
                (denominator, low_labels, high_label, high_residue, unit, excess)
            )
            active_denominators[denominator] += 1
            zero_supremum_classes += high_label is None
    scalar_classes = tuple(scalar_classes)
    require(scalar_classes, (case, "hard class universe vanished"))

    cell_tuple = TORSION.complete_body_cells(carrier, modulus)
    cells = np.asarray(cell_tuple, dtype=np.int64)
    require(cells.size > 0, "complete-cell universe vanished")
    mask_cache = {}

    def safe_mask(label):
        cached = mask_cache.get(label)
        if cached is None:
            residue = (cells * label) % modulus
            cached = (14 * residue >= modulus) & (14 * (residue + label) <= 13 * modulus)
            mask_cache[label] = cached
        return cached

    first_mask = safe_mask(first)
    safe_cache = {}
    section_rows = []
    failures = []
    count_passes = 0
    residue_passes = 0
    for denominator, low_labels, high_label, high_residue, unit, excess in scalar_classes:
        safe = safe_cache.get(low_labels)
        if safe is None:
            mask = first_mask.copy()
            for label in low_labels:
                mask &= safe_mask(label)
            safe = cells[mask]
            safe_cache[low_labels] = safe
        kappa = (denominator + 6) // 7
        count_threshold = kappa * (modulus // denominator)
        count_surplus = int(safe.size) - count_threshold
        residue_count = None
        residue_surplus = None
        if count_surplus > 0:
            route = "COUNT"
            count_passes += 1
        else:
            residue_count = int(np.unique(safe % denominator).size)
            residue_surplus = residue_count - kappa
            if residue_surplus > 0:
                route = "RESIDUE"
                residue_passes += 1
            else:
                route = "OPEN"
                failures.append((denominator, low_labels, high_label, unit, residue_count, kappa))
        section_rows.append((
            denominator, low_labels, high_label, high_residue, unit, excess,
            int(safe.size), kappa, count_threshold, count_surplus,
            residue_count, residue_surplus, route,
        ))
    section_rows = tuple(section_rows)
    safe_sizes = tuple(int(values.size) for values in safe_cache.values())
    control_rows = []
    control_labels = sorted(safe_cache)
    for low_labels in dict.fromkeys((control_labels[0], control_labels[-1])):
        scalar_safe = TORSION.fixed_safe_cells(cell_tuple, modulus, (first, *low_labels))
        vector_safe = tuple(int(cell) for cell in safe_cache[low_labels])
        require(vector_safe == scalar_safe, (case, low_labels, "vector/scalar cell mismatch"))
        control_rows.append((low_labels, len(scalar_safe), sha256(repr(scalar_safe).encode()).hexdigest()))
    control_hash = sha256(repr(tuple(control_rows)).encode()).hexdigest()
    section_hash = sha256(repr(section_rows).encode()).hexdigest()
    return (
        case, modulus, int(cells.size), len(scalar_classes), len(safe_cache),
        (min(safe_sizes), max(safe_sizes)), min(row[5] for row in scalar_classes),
        count_passes, residue_passes, tuple(failures),
        len(active_denominators), tuple(active_denominators.items()),
        min(row[9] for row in section_rows),
        min((row[11] for row in section_rows if row[11] is not None), default=None),
        section_hash,
        sha256(repr(tuple(int(value) for value in amplitudes)).encode()).hexdigest(),
        zero_supremum_classes, control_hash,
    )


def render(rows, audits):
    require(tuple((row[1], row[2]) for row in rows) == CASES, "z1612 universe changed")
    ordinary = tuple(row for row in rows if row[0] == "ALL-LABEL")
    high = tuple(row for row in rows if row[0] == "FORCED-HIGH")
    require((len(ordinary), len(high)) == (8, 11), "closure partition changed")
    require(all(row[16] == row[17] for row in ordinary), "ordinary projected packet survived")
    scalar_high = tuple(row for row in high if row[11] > 0)
    hard_high = tuple(row for row in high if row[11] <= 0)
    require(len(scalar_high) == 7, "scalar-high partition changed")
    require(tuple((row[1], row[2]) for row in hard_high) == HARD_CASES, "hard-high universe changed")
    require(tuple(audit[0] for audit in audits) == HARD_CASES, "section audit universe changed")
    require(all(not audit[9] for audit in audits), ("complete-section failures", audits))
    require(all(audit[7] + audit[8] == audit[3] for audit in audits), "section ledger mismatch")

    ordinary_totals = tuple(sum(row[12][index] for row in ordinary) for index in range(4))
    certificate_total = sum(row[14] for row in ordinary)
    packet_total = sum(row[16] for row in ordinary)
    profile_hash = sha256(repr(tuple(rows)).encode()).hexdigest()
    section_hash = sha256(repr(tuple(audits)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    if EXPECTED_SECTION_SHA256 is not None:
        require(section_hash == EXPECTED_SECTION_SHA256, "section digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256, EXPECTED_BASE_OUTPUT_SHA256, CASES,
        ordinary_totals, certificate_total, packet_total,
        tuple((row[1], row[2], row[11]) for row in scalar_high),
        tuple((audit[0], audit[3], audit[7], audit[8], audit[14], audit[17]) for audit in audits),
        1611, profile_hash, section_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact z1612 shell closure",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        "inherited_cap=1615;atlas_gap=1613..1615 empty;boundary_rows=z1612:19",
        "closure_partition=all_label:8;forced_high_scalar:7;forced_high_section:4;empty:19;survivors:0",
        f"ordinary_ledger=scalar:{ordinary_totals[0]};crude:{ordinary_totals[1]};status:{ordinary_totals[2]};status_survivors:{ordinary_totals[3]};exact_status_certificates:{certificate_total};literal_packets:{packet_total};projected_kills:{packet_total}",
    ]
    for row in rows:
        if row[0] == "ALL-LABEL":
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
        else:
            conclusion = "SCALAR-EMPTY" if row[11] > 0 else "SECTION-CLOSURE-BELOW"
            lines.append(
                f"CASE;route=FORCED-HIGH;z1={row[1]};E={','.join(map(str,row[2]))};"
                f"h={ftext(row[3])};r={row[4]};L={row[5]};high_floor={row[6]};"
                f"delta1={ftext(row[7])};lower={ftext(row[8])};two_high_gap={ftext(row[10])};"
                f"exact_one_high_gap={ftext(row[11])};top_lows={row[12]};best_high_d={row[13]};"
                f"best_high_value={ftext(row[14])};best_high_label={row[15]};"
                f"amplitude_sha256={row[16]};conclusion={conclusion}"
            )
    for audit in audits:
        lines.append(
            f"HIGH-SECTION;z1={audit[0][0]};E={','.join(map(str,audit[0][1]))};"
            f"L={audit[1]};complete_cells={audit[2]};classes={audit[3]};"
            f"low_triples={audit[4]};safe_cell_range={audit[5]};"
            f"minimum_excess={ftext(audit[6])};count_gate={audit[7]};"
            f"residue_gate={audit[8]};failures={len(audit[9])};"
            f"active_denominators={audit[10]};minimum_count_surplus={audit[12]};"
            f"minimum_residue_surplus={audit[13]};section_sha256={audit[14]};"
            f"amplitude_sha256={audit[15]};zero_supremum_classes={audit[16]};"
            f"scalar_vector_control_sha256={audit[17]};conclusion=EMPTY"
        )
    lines.extend((
        "hostile_controls=four one-high scalar gaps are nonpositive;every exact class closes pointwise by the THM-2984 complete-section gate",
        "consequence=projected k=2 first drift label z1<=1611",
        "scope_boundary=finite-exact THM-2980 input only;full reserved theorem and LRC(14) remain open",
        f"profile_sha256={profile_hash}",
        f"section_sha256={section_hash}",
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
        rows = tuple(closure_profile(case) for case in CASES)
        audits = tuple(high_section_audit(case) for case in HARD_CASES)
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = tuple(pool.imap(closure_profile, CASES))
            audits = tuple(pool.imap(high_section_audit, HARD_CASES))
    payload = render(rows, audits)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
