#!/usr/bin/env python3
"""Exact common-torsion and finite-height reduction for THM-3078.

The two THM-3078 bodies with a nonpositive two-high scalar threshold contain
1,699 residual denominator masks.  This companion first applies a pointwise
common-modulus capacity certificate.  For the unresolved masks it then uses a
complete physical 1/L-cell, the sharp three-aligned union cap, and the exact
one-comb interval discrepancy to bound both high labels.

The mask bank is a text projection of the pinned THM-3078 screen checkpoints.
Use ``--extract-bank`` after a full THM-3078 run to rebuild it.  The default
verification reads the tracked bank and recomputes every new certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import pickle
import tempfile
from bisect import bisect_left
from collections import defaultdict
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3078 = (
    ROOT
    / "04-computation/lrc14_j7_k3_z234_direct_farkas_four_two_high_boundary_thm3078.py"
)
OUTPUT_3078 = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k3_z234_direct_farkas_four_two_high_boundary_thm3078.out"
)
BANK = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k3_z234_two_high_mask_bank_thm3087.tsv"
)
OUTPUT = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k3_z234_two_high_common_torsion_finite_height_thm3087.out"
)

SOURCE_3078_SHA256 = "2a051babe109f56056fe61476870f8e2e13cfc99b2f9bb7ac122b8780c8fa168"
OUTPUT_3078_SHA256 = "d5fc52e922b7e083f5cbe5aa5a6066304c4e5c8c7963904dd3b649477caf5e42"
SEMANTIC_3078_SHA256 = "749a2292da58925c9653ce919d5bab6b374f64f8b713ac52efa1f18fba74c918"
BANK_SHA256 = "b3d424991e816bd5ddc0540554483cd3bc40efbaf44a3f3064197c337ae6ca8d"
EXPECTED_SEMANTIC_SHA256 = "e7f27276343813afa2b34471e8d3b98a28f32491b2be49204e6ca5acc3bc144d"

LEVEL = 234
ALIGNED_UNION_CAP = F(36, 91)
INFINITE = {
    (1, 5, 9, 11, 12, 14): {
        "low": 243,
        "excess": F(34729, 119189070),
        "masks": 886,
        "clean": 44062,
        "torsion_pass": 663,
        "torsion_fail": 223,
        "positive_slack": (1575, 17434),
        "worst_slack": -80678,
        "finite_candidates": 36200669,
    },
    (1, 9, 10, 11, 12, 14): {
        "low": 260,
        "excess": F(3307, 4414410),
        "masks": 813,
        "clean": 44616,
        "torsion_pass": 382,
        "torsion_fail": 431,
        "positive_slack": (700, 13482),
        "worst_slack": -80124,
        "finite_candidates": 72765829,
    },
}
EXPECTED_EXCEPTIONS = {
    97020: (F(1, 7), 14850, F(2162160, 29), F(47, 637)),
    129360: (F(3, 28), 14850, F(2882880, 43), F(435, 2548)),
    145530: (F(2, 21), 14851, F(324324, 5), F(388, 1911)),
}
DISPLAY = {
    (1, 5, 6, 9, 12, 14): (
        (LEVEL, 8820, 324, 1939),
        F(191990179, 1000465830),
        502,
    ),
    (1, 5, 9, 11, 12, 14): (
        (LEVEL, 20174, 19308, 243),
        F(7803086659, 46609560270),
        9504,
    ),
    (1, 9, 10, 11, 12, 14): (
        (LEVEL, 20370, 260, 19521),
        F(7485939823, 44226712530),
        13482,
    ),
    (2, 5, 9, 11, 12, 14): (
        (LEVEL, 20174, 19308, 243),
        F(302562628919, 1794468070395),
        8566,
    ),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes()
    require(b"\r" not in payload.replace(b"\r\n", b""), f"bare CR in {path}")
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(sha(SOURCE_3078) == SOURCE_3078_SHA256, "THM-3078 source changed")
require(sha(OUTPUT_3078) == OUTPUT_3078_SHA256, "THM-3078 output changed")
require(
    f"semantic_sha256={SEMANTIC_3078_SHA256}" in OUTPUT_3078.read_text(),
    "THM-3078 semantic changed",
)
thm = load("thm3078_two_high_finite_reduction_base", SOURCE_3078)
eng = thm.eng


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def body_text(body):
    return ",".join(map(str, body))


def parse_body(text):
    body = tuple(map(int, text.split(",")))
    require(len(body) == 6, text)
    return body


def cell_clean(cell, label, L):
    residue = (cell * label) % L
    return 14 * residue >= L and 14 * (residue + label) <= 13 * L


def mask_rows_from_residual(body, residual):
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(body)
    needed = {
        d for ds in residual for d in eng.suffix_slots(ds, stream.first_d)
    }
    low, _high, _signs, _checks = eng.build_literal_tables(stream, needed)
    required = stream.lower - stream.first_delta
    rows = []
    for ds in residual:
        slots = eng.suffix_slots(ds, stream.first_d)
        for low_index, low_d in enumerate(slots):
            high_ds = tuple(
                d for index, d in enumerate(slots) if index != low_index
            )
            for value, label in low[low_d]:
                if value >= required:
                    rows.append((body, high_ds[0], high_ds[1], label,
                                 value - required))
    rows.sort()
    return tuple(rows)


def extract_bank(checkpoint_dir, target):
    rows_by_body = {}
    pattern = f"screen-{thm.CACHE_FINGERPRINT[:16]}-*.pickle"
    for path in checkpoint_dir.glob(pattern):
        fingerprint, task, row = pickle.loads(path.read_bytes())
        require(fingerprint == thm.CACHE_FINGERPRINT, (path, "fingerprint"))
        require(task[:2] == row[:2], (path, "task/row"))
        if task[1] in INFINITE:
            rows_by_body[task[1]] = mask_rows_from_residual(task[1], row[13])
    require(set(rows_by_body) == set(INFINITE), ("missing bodies", rows_by_body))
    rows = tuple(
        row for body in sorted(rows_by_body) for row in rows_by_body[body]
    )
    require(len(rows) == 1699, len(rows))
    lines = ["body\td1\td2\tlow\texcess"]
    lines.extend(
        f"{body_text(body)}\t{d1}\t{d2}\t{low}\t{ftext(excess)}"
        for body, d1, d2, low, excess in rows
    )
    target.write_text("\n".join(lines) + "\n", newline="\n")
    print(f"bank_rows={len(rows)};bank_sha256={sha(target)}")


def read_bank():
    if BANK_SHA256 is not None:
        require(sha(BANK) == BANK_SHA256, "mask bank changed")
    lines = BANK.read_text().splitlines()
    require(lines and lines[0] == "body\td1\td2\tlow\texcess", "bank header")
    rows = []
    for line in lines[1:]:
        body_s, d1_s, d2_s, low_s, excess_s = line.split("\t")
        numerator, denominator = map(int, excess_s.split("/"))
        rows.append(
            (
                parse_body(body_s),
                int(d1_s),
                int(d2_s),
                int(low_s),
                F(numerator, denominator),
            )
        )
    require(rows == sorted(rows), "bank order")
    require(len(rows) == len(set(rows)) == 1699, "bank cardinality")
    return tuple(rows)


def fixed_safe_cells(stream, labels):
    cells = []
    for left, right in stream.ranges:
        for cell in range(left, right):
            if all(cell_clean(cell, label, stream.L) for label in labels):
                cells.append(cell)
    return tuple(cells)


def common_torsion_record(cells, d1, d2):
    Q = lcm(d1, d2)
    support = len({cell % Q for cell in cells})
    cap1 = ((d1 + 6) // 7) * (Q // d1)
    cap2 = ((d2 + 6) // 7) * (Q // d2)
    return (d1, d2, Q, support, cap1 + cap2, support - cap1 - cap2)


def normalized_danger_mass(cell, label, L):
    """Exact normalized load of D_label on [cell/L,(cell+1)/L]."""
    require(0 < label < L and L % 14 == 0, (label, L))
    left = (cell * label) % L
    right = left + label
    radius = L // 14
    hit = 0
    first_k = (left - radius) // L - 1
    last_k = (right + radius) // L + 1
    for index in range(first_k, last_k + 1):
        danger_left = index * L - radius
        danger_right = index * L + radius
        hit += max(0, min(right, danger_right) - max(left, danger_left))
    return F(hit, label)


def count_pairs(first, second, excess, same):
    second_values = [value for value, _label in second]
    total = 0
    for value, _label in first:
        index = bisect_left(second_values, -excess - value)
        total += len(second) - index
        if same and 2 * value + excess >= 0:
            total -= 1
    return total


def verify_display_controls():
    rows = []
    for body, (labels, expected_mass, expected_slack) in DISPLAY.items():
        eng.FIRST = LEVEL
        eng.ray.FIRST = LEVEL
        stream = eng.ray.Stream(body)
        low = tuple(label for label in labels[1:] if label < stream.high_floor)
        high = tuple(label for label in labels[1:] if label >= stream.high_floor)
        require(len(low) == 1 and len(high) == 2, (body, labels))
        cells = fixed_safe_cells(stream, (LEVEL, *low))
        denominators = tuple(stream.L // gcd(stream.L, label) for label in high)
        torsion = common_torsion_record(cells, *denominators)
        require(torsion[-1] == expected_slack > 0, (body, torsion))
        mass, ruler, body_components, danger_components = (
            thm.exact_maximizer_safe_mass(body, labels)
        )
        require(mass == expected_mass > 0, (body, mass))
        rows.append(
            (body, stream.L, labels, denominators, torsion, mass, ruler,
             body_components, danger_components)
        )
    return tuple(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--extract-bank", action="store_true")
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "lrc14-k3-z234-thm3078-v1",
    )
    args = parser.parse_args()
    if args.extract_bank:
        extract_bank(args.checkpoint_dir, BANK)
        return

    bank = read_bank()
    grouped = defaultdict(list)
    for row in bank:
        grouped[row[0]].append(row)
    require(set(grouped) == set(INFINITE), grouped.keys())

    body_records = []
    unresolved = []
    clean_cache = {}
    for body in sorted(INFINITE):
        expected = INFINITE[body]
        eng.FIRST = LEVEL
        eng.ray.FIRST = LEVEL
        stream = eng.ray.Stream(body)
        require((stream.L, stream.high_floor) == (194040, 19111), body)
        low = expected["low"]
        require(
            eng.delta(stream, low)
            - (stream.lower - stream.first_delta)
            == expected["excess"],
            (body, "low excess"),
        )
        rows = tuple(grouped[body])
        require(len(rows) == expected["masks"], (body, len(rows)))
        require(all(row[3:] == (low, expected["excess"]) for row in rows), body)
        cells = fixed_safe_cells(stream, (LEVEL, low))
        require(len(cells) == expected["clean"], (body, len(cells)))
        clean_cache[body] = cells
        torsion = tuple(
            common_torsion_record(cells, row[1], row[2]) for row in rows
        )
        passed = tuple(record for record in torsion if record[-1] > 0)
        failed = tuple(record for record in torsion if record[-1] <= 0)
        positive_slacks = tuple(record[-1] for record in passed)
        require(len(passed) == expected["torsion_pass"], (body, "pass"))
        require(len(failed) == expected["torsion_fail"], (body, "fail"))
        require(
            (min(positive_slacks), max(positive_slacks))
            == expected["positive_slack"],
            (body, "positive slack"),
        )
        require(min(record[-1] for record in failed) == expected["worst_slack"],
                (body, "worst slack"))
        require(sum(record[-1] == 0 for record in failed) == 1,
                (body, "equality census"))
        require(
            tuple(record[:2] for record in failed if record[-1] == 0)
            == ((2, 2),),
            (body, "equality pair"),
        )
        for row, record in zip(rows, torsion):
            if record[-1] <= 0:
                unresolved.append((row, record))
        body_records.append((body, stream.L, stream.high_floor, low,
                             expected["excess"], len(cells), torsion))

    require(len(unresolved) == 654, len(unresolved))

    # Uniform two-stage physical height bound.  The first stage needs only one
    # complete base cell and therefore applies before the torsion split.
    L = 194040
    high_floor = 19111
    initial_boundary = F(1092 * L, 1421)
    initial_ceiling = initial_boundary.numerator // initial_boundary.denominator
    require(initial_boundary == F(4324320, 29), initial_boundary)
    require(initial_ceiling == 149114, initial_ceiling)
    final_ceiling = F(13 * L, 49)
    require(final_ceiling.denominator == 1 and final_ceiling == 51480,
            final_ceiling)

    exception_records = []
    for body in sorted(INFINITE):
        cells = clean_cache[body]
        exceptions = []
        for label in range(high_floor, initial_ceiling + 1):
            if not any(cell_clean(cell, label, L) for cell in cells):
                exceptions.append(label)
        require(tuple(exceptions) == tuple(EXPECTED_EXCEPTIONS),
                (body, exceptions))
        for label in exceptions:
            best_mass, best_cell = min(
                (normalized_danger_mass(cell, label, L), cell)
                for cell in cells
            )
            expected_mass, expected_cell, cutoff, margin = (
                EXPECTED_EXCEPTIONS[label]
            )
            require((best_mass, best_cell) == (expected_mass, expected_cell),
                    (body, label, best_mass, best_cell))
            derived_cutoff = F(6 * L, 49) / (F(6, 13) - best_mass)
            require(derived_cutoff == cutoff < label,
                    (body, label, derived_cutoff))
            endpoint_rhs = (
                ALIGNED_UNION_CAP + best_mass + F(1, 7)
                + F(6 * L, 49 * label)
            )
            require(1 - endpoint_rhs == margin > 0,
                    (body, label, endpoint_rhs))
            exception_records.append(
                (body, label, best_mass, best_cell, cutoff, margin)
            )

    # Exact finite candidate census after the common-torsion closures.
    finite_records = []
    for body in sorted(INFINITE):
        expected = INFINITE[body]
        stream = eng.ray.Stream(body)
        body_unresolved = [row for row, record in unresolved if row[0] == body]
        needed = {d for row in body_unresolved for d in row[1:3]}
        values = defaultdict(list)
        for label in range(high_floor, int(final_ceiling) + 1):
            d = L // gcd(L, label)
            if d in needed:
                value = eng.delta(stream, label)
                require(
                    (label + L) * eng.delta(stream, label + L) == label * value,
                    (body, label, "ray law"),
                )
                values[d].append((value, label))
        for d in values:
            values[d].sort()
        total = 0
        raw_total = 0
        zero_masks = 0
        per_mask = []
        for row in body_unresolved:
            _body, d1, d2, _low, excess = row
            first = values[d1]
            second = values[d2]
            count = count_pairs(first, second, excess, d1 == d2)
            raw = len(first) * len(second) - (len(first) if d1 == d2 else 0)
            require(count == raw, (body, d1, d2, count, raw))
            total += count
            raw_total += raw
            zero_masks += count == 0
            per_mask.append((d1, d2, count))
        require(total == raw_total == expected["finite_candidates"],
                (body, total, raw_total))
        finite_records.append(
            (body, len(body_unresolved), total, zero_masks, tuple(per_mask))
        )

    display = verify_display_controls()
    require(sum(row[2] for row in finite_records) == 108966498,
            "finite total")
    require(sum(INFINITE[body]["torsion_pass"] for body in INFINITE) == 1045,
            "torsion total")

    bank_hash = sha(BANK)
    semantic_packet = (
        bank_hash,
        tuple(body_records),
        tuple(exception_records),
        tuple(finite_records),
        display,
        initial_boundary,
        final_ceiling,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 z234 two-high common-torsion and finite-height reduction",
        f"dependency=THM3078_source:{SOURCE_3078_SHA256};output:{OUTPUT_3078_SHA256};semantic:{SEMANTIC_3078_SHA256}",
        f"mask_bank=rows:{len(bank)};sha256:{bank_hash};bodies:2;nonpositive_threshold:1699",
        "common_torsion=closed:1045;unresolved:654;criterion:|S_mod_Q|>ceil(d1/7)Q/d1+ceil(d2/7)Q/d2;pointwise_projection:T",
    ]
    for body, _L, _H, low, excess, clean, torsion in body_records:
        positive = [row[-1] for row in torsion if row[-1] > 0]
        failed = [row[-1] for row in torsion if row[-1] <= 0]
        lines.append(
            f"BODY;E={body};low={low};excess={ftext(excess)};clean_cells={clean};"
            f"masks={len(torsion)};torsion_closed={len(positive)};unresolved={len(failed)};"
            f"positive_slack={min(positive)}..{max(positive)};worst={min(failed)};equalities={sum(x == 0 for x in failed)}@high_pair(2,2)"
        )
    lines.extend(
        [
            f"height_stage1=if_both_high>{ftext(initial_boundary)}_then_36/91+2/7+12L/(49N)<1;therefore_min_high<={initial_ceiling}",
            f"height_stage2=complete_cell_safe_for_every_{high_floor}<=z<={initial_ceiling}_except:{','.join(map(str, EXPECTED_EXCEPTIONS))};ordinary_cutoff:other_high<={int(final_ceiling)}",
            "exception_controls="
            + ";".join(
                f"z{label}:load={ftext(mass)},cell={cell},far_cutoff={ftext(cutoff)},margin_at_w=z={ftext(margin)}"
                for _body, label, mass, cell, cutoff, margin
                in exception_records[:3]
            ),
            f"finite_contract=both_high_in_[{high_floor},{int(final_ceiling)}];unresolved_masks:654;ordered_scalar_candidates:108966498",
        ]
    )
    for body, masks, total, zero_masks, per_mask in finite_records:
        lines.append(
            f"FINITE;E={body};masks={masks};ordered_candidates={total};zero_masks={zero_masks};record_sha256={hashlib.sha256(repr(per_mask).encode()).hexdigest()}"
        )
    lines.append(
        "displayed_maximizer_control="
        + ";".join(
            f"E={body}:P=T,slack={torsion[-1]},raw_safe={ftext(mass)}"
            for body, _L, _labels, _ds, torsion, mass, _ruler, _bc, _dc
            in display
        )
    )
    lines.extend(
        [
            "first_failed_implication=finite_height_or_positive_raw_safe_mass=>mu(P_(E,Z))>36/91_for_every_remaining_packet",
            "missing_coordinate=the_three-aligned_tail_section_U_A_inside_the_exact_four-drift_projection_P_(E,Z);the_108966498_candidates_are_not_closed",
            "scope=no_ledger_or_projected_cap_decrement;the_common-torsion successes_close_mask_families_but_654_finite_masks_remain",
            f"semantic_sha256={semantic}",
        ]
    )
    args.output.write_text("\n".join(lines) + "\n", newline="\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
