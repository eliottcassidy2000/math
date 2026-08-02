#!/usr/bin/env python3
"""Exact closure of the final two THM-3078 z1=234 boundary bodies.

The tracked bank is the sorted projection of the two relevant THM-3078
screen checkpoint residuals.  The proof has two different finite exits:

* the 16-mask body has no denominator-2 high label in its exact height box;
* the 1,138-mask body has 1,018 literal scalar candidates, and every one has
  a common complete body/234/low/high/high cell.

Use ``--extract-bank`` after a complete THM-3078 run to rebuild the bank.
The default verification is checkpoint-independent and recomputes every
truth-bearing count with explicit checks that remain active under ``-O``.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import pickle
import tempfile
from collections import defaultdict
from fractions import Fraction as F
from math import gcd
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3078 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z234_direct_farkas_four_two_high_boundary_thm3078.py"
)
OUTPUT_3078 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z234_direct_farkas_four_two_high_boundary_thm3078.out"
)
SOURCE_3094 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z234_two_high_complete_cell_intersection_closure_thm3094.py"
)
OUTPUT_3094 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z234_two_high_complete_cell_intersection_closure_thm3094.out"
)
SOURCE_3071 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z236_z235_compositional_descent_thm3071.py"
)
OUTPUT_3071 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z236_z235_compositional_descent_thm3071.out"
)
BANK = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z234_final_two_body_mask_bank_thm3098.tsv"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z234_final_two_body_height_complete_cell_closure_thm3098.out"
)

SOURCE_3078_SHA256 = "2a051babe109f56056fe61476870f8e2e13cfc99b2f9bb7ac122b8780c8fa168"
OUTPUT_3078_SHA256 = "d5fc52e922b7e083f5cbe5aa5a6066304c4e5c8c7963904dd3b649477caf5e42"
SEMANTIC_3078_SHA256 = "749a2292da58925c9653ce919d5bab6b374f64f8b713ac52efa1f18fba74c918"
SOURCE_3094_SHA256 = "c456e28c88cac8d0b687c1636463aa79642b2fb39ea187eafc836f46e3ab3e0e"
OUTPUT_3094_SHA256 = "6752dbac509aba6c739eb9ce3cd59da85d39ae47b6f4f60b1657ecaf4cf36a95"
SEMANTIC_3094_SHA256 = "7c3fa428224d59b717770b936bb2d35aea2642dec83c77df8395c538a8015682"
SOURCE_3071_SHA256 = "c8735a0e1328b08e98e9f27b86f901d2b0c491f832381d98f7a8a14f11e4f345"
OUTPUT_3071_SHA256 = "c04e6f57ac7025100645f1a0f546e3e5f79f6444fa2269469c09b38e746e772f"
SEMANTIC_3071_SHA256 = "d4cfb4c7e497007f044041e9f7e56d3b109107d033421d592f7f3f74e0995f02"
BANK_SHA256 = "a1b3959152cc87b05d352492870803cc23e04fb6b38dd7cbf09d2e9605850144"
EXPECTED_SEMANTIC_SHA256 = "f5427dc8392ed287162f076df9b5e418afcec2a2f3a63a3b94a840c313a5b6bd"

LEVEL = 234
ALIGNED_UNION_CAP = F(36, 91)
LEDGER_BEFORE = 374768
LAYER_ROWS = 381
LEDGER_AFTER = 374387

BODY_A = (1, 5, 6, 9, 12, 14)
BODY_B = (2, 5, 9, 11, 12, 14)
OTHER_CLOSED_BODIES = {
    (1, 5, 9, 11, 12, 14),
    (1, 9, 10, 11, 12, 14),
}

EXPECTED = {
    BODY_A: {
        "L": 17640,
        "H": 1738,
        "masks": 16,
        "bank_digest": "fb53e6a6c3e2098560a8971a079aa595a04a95c1fcdff77e3be3a8691bfff6d5",
        "three_gap": F(13734121, 4337473140),
        "three_ds": (2, 980, 1764, 1960),
        "raw_expansions": 2748,
        "survivors": 3,
        "survivor_digest": "1fd3d18d3a4fd592a2391643fe7fdef1ed196e2a237916093e611c3bba76a87f",
        "low": 324,
        "deficit": F(37, 3611790),
        "survivor_margin": (F(31396, 122542875), F(26554, 71461845)),
        "clean": 4446,
        "initial_boundary": F(393120, 29),
        "initial_ceiling": 13555,
        "final_ceiling": 4680,
        "exceptions": (
            (8820, F(1, 7), 1350, F(196560, 29), F(47, 637)),
            (11760, F(3, 28), 1350, F(262080, 43), F(435, 2548)),
            (13230, F(2, 21), 1591, F(29484, 5), F(388, 1911)),
        ),
        "literal_candidates": 0,
        "raw_pairs": 0,
        "endpoint": (14, 1350),
    },
    BODY_B: {
        "L": 194040,
        "H": 19111,
        "masks": 1138,
        "bank_digest": "c5c1949cc5a2037ab6af71f8d41d36b6f2161fce0fd710a84315a3a3f4f1218f",
        "three_gap": F(3429147350803, 807271369268364),
        "three_ds": (1260, 10780, 16170, 21560),
        "raw_expansions": 1512196,
        "survivors": 69,
        "survivor_digest": "ae7d6b1b307992b48ade07543354798f85ea4d4136eb1b1de01b58b959151a77",
        "low": 243,
        "deficit": F(4117, 39729690),
        "survivor_margin": (F(1063, 30830239440), F(138056623, 3349673731404)),
        "clean": 44364,
        "initial_boundary": F(4324320, 29),
        "initial_ceiling": 149114,
        "final_ceiling": 51480,
        "exceptions": (
            (97020, F(1, 7), 6930, F(2162160, 29), F(47, 637)),
            (129360, F(3, 28), 6930, F(2882880, 43), F(435, 2548)),
            (145530, F(2, 21), 6931, F(324324, 5), F(388, 1911)),
        ),
        "literal_candidates": 1018,
        "raw_pairs": 7366342,
        "packet_digest": "3277b423430c9fe955d78e056942b0e73e0f8ec740f2f240e108451d46ceefe0",
        "pair_digest": "b723c2705b0f25c04e5b07307908ccdfdd666e9ef4501d4baff1cccf41b22aff",
        "candidate_labels": 342,
        "label_range": (19166, 48540),
        "minimum_bonferroni": (15470, (20174, 48540)),
        "minimum_intersection": (19908, (20174, 48540), 6930),
        "endpoint": (2, 6930),
    },
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


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def body_text(body):
    return ",".join(map(str, body))


def parse_body(text):
    body = tuple(map(int, text.split(",")))
    require(len(body) == 6, text)
    return body


require(sha(SOURCE_3078) == SOURCE_3078_SHA256, "THM-3078 source changed")
require(sha(OUTPUT_3078) == OUTPUT_3078_SHA256, "THM-3078 output changed")
require(
    f"semantic_sha256={SEMANTIC_3078_SHA256}" in OUTPUT_3078.read_text(),
    "THM-3078 semantic changed",
)
require(sha(SOURCE_3094) == SOURCE_3094_SHA256, "THM-3094 source changed")
require(sha(OUTPUT_3094) == OUTPUT_3094_SHA256, "THM-3094 output changed")
require(
    f"semantic_sha256={SEMANTIC_3094_SHA256}" in OUTPUT_3094.read_text(),
    "THM-3094 semantic changed",
)
require(sha(SOURCE_3071) == SOURCE_3071_SHA256, "THM-3071 source changed")
require(sha(OUTPUT_3071) == OUTPUT_3071_SHA256, "THM-3071 output changed")
require(
    f"semantic_sha256={SEMANTIC_3071_SHA256}" in OUTPUT_3071.read_text(),
    "THM-3071 semantic changed",
)

thm = load("thm3078_final_two_body_base", SOURCE_3078)
eng = thm.eng


def extract_bank(checkpoint_dir, target):
    residual_by_body = {}
    pattern = f"screen-{thm.CACHE_FINGERPRINT[:16]}-*.pickle"
    for path in checkpoint_dir.glob(pattern):
        fingerprint, task, row = pickle.loads(path.read_bytes())
        require(fingerprint == thm.CACHE_FINGERPRINT, (path, "fingerprint"))
        require(task[:2] == row[:2], (path, "task/row"))
        if task[1] in EXPECTED:
            residual_by_body[task[1]] = tuple(sorted(row[13]))
    require(set(residual_by_body) == set(EXPECTED), residual_by_body.keys())
    lines = ["body\td1\td2\td3\td4"]
    for body in sorted(residual_by_body):
        rows = residual_by_body[body]
        expected = EXPECTED[body]
        require(len(rows) == expected["masks"], (body, len(rows)))
        require(
            hashlib.sha256(repr(rows).encode()).hexdigest()
            == expected["bank_digest"],
            (body, "residual digest"),
        )
        lines.extend(
            f"{body_text(body)}\t" + "\t".join(map(str, ds)) for ds in rows
        )
    target.write_text("\n".join(lines) + "\n", newline="\n")
    print(f"bank_rows={sum(x['masks'] for x in EXPECTED.values())};sha256={sha(target)}")


def read_bank():
    if BANK_SHA256 is not None:
        require(sha(BANK) == BANK_SHA256, "THM-3098 mask bank changed")
    lines = BANK.read_text().splitlines()
    require(lines and lines[0] == "body\td1\td2\td3\td4", "bank header")
    grouped = defaultdict(list)
    flat = []
    for line in lines[1:]:
        fields = line.split("\t")
        require(len(fields) == 5, line)
        body = parse_body(fields[0])
        ds = tuple(map(int, fields[1:]))
        require(ds == tuple(sorted(ds)), (body, ds, "unsorted mask"))
        grouped[body].append(ds)
        flat.append((body, ds))
    require(flat == sorted(flat), "bank order")
    require(len(flat) == len(set(flat)) == 1154, "bank cardinality")
    require(set(grouped) == set(EXPECTED), grouped.keys())
    for body in EXPECTED:
        rows = tuple(grouped[body])
        require(len(rows) == EXPECTED[body]["masks"], (body, len(rows)))
        require(
            hashlib.sha256(repr(rows).encode()).hexdigest()
            == EXPECTED[body]["bank_digest"],
            (body, "bank body digest"),
        )
    return {body: tuple(grouped[body]) for body in grouped}, tuple(flat)


def cell_clean(cell, label, L):
    residue = (cell * label) % L
    return 14 * residue >= L and 14 * (residue + label) <= 13 * L


def fixed_safe_cells(stream, low):
    cells = thm.vector_fixed_safe_cells(stream, (low,))
    require(
        all(
            cell_clean(int(cell), label, stream.L)
            for cell in cells[: min(128, len(cells))]
            for label in (LEVEL, low)
        ),
        (stream.body, "vector fixed-safe hostile sample"),
    )
    return cells


def normalized_danger_mass(cell, label, L):
    """Exact value of L times the danger mass on one complete 1/L cell."""
    require(0 < label < L and L % 14 == 0, (label, L))
    left = (cell * label) % L
    right = left + label
    radius = L // 14
    hit = 0
    first_index = (left - radius) // L - 1
    last_index = (right + radius) // L + 1
    for index in range(first_index, last_index + 1):
        danger_left = index * L - radius
        danger_right = index * L + radius
        hit += max(0, min(right, danger_right) - max(left, danger_left))
    return F(hit, label)


def safe_mask(cells, label, L):
    residues = (cells * int(label)) % int(L)
    return (14 * residues >= L) & (14 * (residues + int(label)) <= 13 * L)


def scalar_reduction(body, residual):
    expected = EXPECTED[body]
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(body)
    require(
        (stream.L, stream.high_floor) == (expected["L"], expected["H"]),
        (body, stream.L, stream.high_floor),
    )
    needed = {d for ds in residual for d in eng.suffix_slots(ds, stream.first_d)}
    low, high, low_signs, ray_checks = eng.build_literal_tables(stream, needed)
    required = stream.lower - stream.first_delta

    three_records = []
    for ds in residual:
        slots = eng.suffix_slots(ds, stream.first_d)
        rows = tuple((d, high[d][0], high[d][1]) for d in slots)
        best = sum((row[1] for row in rows), F())
        gap = required - best
        require(gap > 0, (body, ds, "three-high upper survived", gap))
        three_records.append((ds, gap, rows))
    three_records = tuple(sorted(three_records))
    minimum_three = min(three_records, key=lambda row: (row[1], row[0]))
    require(
        (minimum_three[1], minimum_three[0])
        == (expected["three_gap"], expected["three_ds"]),
        (body, "three-high minimum", minimum_three),
    )

    raw_expansions = 0
    survivors = []
    for ds in residual:
        slots = eng.suffix_slots(ds, stream.first_d)
        for low_index, low_d in enumerate(slots):
            high_ds = tuple(
                d for index, d in enumerate(slots) if index != low_index
            )
            best = high[high_ds[0]][0] + high[high_ds[1]][0]
            for value, label in low[low_d]:
                raw_expansions += 1
                deficit = required - value
                margin = best - deficit
                if margin >= 0:
                    survivors.append(
                        (
                            ds,
                            low_index,
                            low_d,
                            high_ds,
                            label,
                            deficit,
                            best,
                            margin,
                        )
                    )
    survivors = tuple(sorted(survivors))
    require(raw_expansions == expected["raw_expansions"], (body, raw_expansions))
    require(len(survivors) == expected["survivors"], (body, len(survivors)))
    require(
        hashlib.sha256(repr(survivors).encode()).hexdigest()
        == expected["survivor_digest"],
        (body, "survivor digest"),
    )
    require(
        all(
            row[4] == expected["low"] and row[5] == expected["deficit"]
            for row in survivors
        ),
        (body, "low literal/deficit"),
    )
    require(
        (min(row[7] for row in survivors), max(row[7] for row in survivors))
        == expected["survivor_margin"],
        (body, "survivor margins"),
    )
    require(not any(row[3][0] == row[3][1] for row in survivors),
            (body, "same high denominator survivor"))
    return stream, low_signs, ray_checks, three_records, survivors


def height_reduction(body, stream, cells):
    expected = EXPECTED[body]
    L = stream.L
    H = stream.high_floor
    initial_boundary = F(1092 * L, 1421)
    initial_ceiling = initial_boundary.numerator // initial_boundary.denominator
    final_ceiling = F(13 * L, 49)
    require(
        (initial_boundary, initial_ceiling, final_ceiling)
        == (
            expected["initial_boundary"],
            expected["initial_ceiling"],
            F(expected["final_ceiling"]),
        ),
        (body, "height constants"),
    )

    exceptions = tuple(
        label
        for label in range(H, initial_ceiling + 1)
        if not any(cell_clean(int(cell), label, L) for cell in cells)
    )
    require(exceptions == tuple(row[0] for row in expected["exceptions"]),
            (body, "height exceptions", exceptions))
    exception_records = []
    for label, expected_mass, expected_cell, expected_cutoff, expected_margin in (
        expected["exceptions"]
    ):
        mass, cell = min(
            (normalized_danger_mass(int(candidate), label, L), int(candidate))
            for candidate in cells
        )
        cutoff = F(6 * L, 49) / (F(6, 13) - mass)
        endpoint_rhs = ALIGNED_UNION_CAP + mass + F(1, 7) + F(6 * L, 49 * label)
        margin = 1 - endpoint_rhs
        require(
            (mass, cell, cutoff, margin)
            == (expected_mass, expected_cell, expected_cutoff, expected_margin),
            (body, label, mass, cell, cutoff, margin),
        )
        require(cutoff < label and margin > 0, (body, label, "exception exit"))
        exception_records.append((label, mass, cell, cutoff, margin))
    return (
        initial_boundary,
        initial_ceiling,
        int(final_ceiling),
        tuple(exception_records),
    )


def literal_candidates(body, stream, survivors, final_ceiling):
    expected = EXPECTED[body]
    needed = {d for row in survivors for d in row[3]}
    values = defaultdict(list)
    for label in range(stream.high_floor, final_ceiling + 1):
        d = stream.L // gcd(stream.L, label)
        if d in needed:
            value = eng.delta(stream, label)
            require(
                (label + stream.L) * eng.delta(stream, label + stream.L)
                == label * value,
                (body, label, "high-ray law"),
            )
            values[d].append((value, label))
    for d in values:
        values[d].sort(key=lambda row: (row[1], row[0]))

    packets = []
    per_survivor = []
    raw_total = 0
    for row in survivors:
        ds, low_index, low_d, high_ds, low_label, deficit, _best, _margin = row
        first = values[high_ds[0]]
        second = values[high_ds[1]]
        raw = len(first) * len(second)
        if high_ds[0] == high_ds[1]:
            raw -= len(first)
        count = 0
        for first_value, first_label in first:
            for second_value, second_label in second:
                if first_label == second_label:
                    continue
                if first_value + second_value >= deficit:
                    packets.append(
                        (
                            ds,
                            low_index,
                            low_d,
                            high_ds,
                            low_label,
                            first_label,
                            second_label,
                        )
                    )
                    count += 1
        raw_total += raw
        per_survivor.append(
            (ds, low_index, low_d, high_ds, low_label, deficit, count, raw)
        )
    packets = tuple(sorted(packets))
    per_survivor = tuple(sorted(per_survivor))
    require(len(packets) == expected["literal_candidates"],
            (body, len(packets)))
    require(raw_total == expected["raw_pairs"], (body, raw_total))
    require(all(row[-2] != row[-1] for row in packets),
            (body, "repeated high label"))
    if body == BODY_A:
        require(not values[2], (body, "unexpected denominator-2 high label"))
        require(all(2 in row[3] for row in survivors),
                (body, "missing denominator-2 obstruction"))
    else:
        require(
            hashlib.sha256(repr(packets).encode()).hexdigest()
            == expected["packet_digest"],
            (body, "packet digest"),
        )
        pairs = tuple(sorted({(row[-2], row[-1]) for row in packets}))
        require(len(pairs) == len(packets), (body, "duplicate packet pair"))
        require(
            hashlib.sha256(repr(pairs).encode()).hexdigest()
            == expected["pair_digest"],
            (body, "pair digest"),
        )
        labels = tuple(sorted({label for pair in pairs for label in pair}))
        require(
            (len(labels), (labels[0], labels[-1]))
            == (expected["candidate_labels"], expected["label_range"]),
            (body, "candidate label census"),
        )
    return packets, per_survivor, raw_total


def complete_cell_closure(body, stream, cells, packets):
    expected = EXPECTED[body]
    if not packets:
        return (), None, None
    pairs = tuple(sorted({(row[-2], row[-1]) for row in packets}))
    labels = tuple(sorted({label for pair in pairs for label in pair}))
    masks = {label: safe_mask(cells, label, stream.L) for label in labels}
    counts = {label: int(mask.sum()) for label, mask in masks.items()}
    records = []
    for first, second in pairs:
        intersection = masks[first] & masks[second]
        intersection_count = int(intersection.sum())
        bonferroni = counts[first] + counts[second] - len(cells)
        require(bonferroni > 0 and intersection_count >= bonferroni,
                (body, first, second, bonferroni, intersection_count))
        witness = int(cells[int(np.flatnonzero(intersection)[0])])
        records.append(
            (
                first,
                second,
                counts[first],
                counts[second],
                bonferroni,
                intersection_count,
                witness,
            )
        )
    records = tuple(records)
    minimum_bonferroni = min(records, key=lambda row: (row[4], row[:2]))
    minimum_intersection = min(records, key=lambda row: (row[5], row[:2]))
    require(
        (minimum_bonferroni[4], minimum_bonferroni[:2])
        == expected["minimum_bonferroni"],
        (body, "minimum Bonferroni", minimum_bonferroni),
    )
    require(
        (
            minimum_intersection[5],
            minimum_intersection[:2],
            minimum_intersection[6],
        )
        == expected["minimum_intersection"],
        (body, "minimum intersection", minimum_intersection),
    )
    return records, minimum_bonferroni, minimum_intersection


def composition_record():
    tasks, neighbor_census = thm.atlas_tasks()
    require(neighbor_census[234] == (381, 330, 51), neighbor_census)
    require(neighbor_census[233] == (62, 45, 17), neighbor_census)
    task_bodies = tuple(task[1] for task in tasks)
    boundary = set(thm.EXPECTED_SURVIVORS)
    require(boundary == set(EXPECTED) | OTHER_CLOSED_BODIES, boundary)
    require(all(body in task_bodies for body in boundary), "boundary row absent")
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")
    require(
        f"promotion_consequence=ledger 374780-12=374768;projected cap z1<=234"
        in OUTPUT_3071.read_text(),
        "THM-3071 ledger formula changed",
    )
    return (
        len(tasks),
        neighbor_census[234],
        neighbor_census[233],
        tuple(sorted(boundary)),
        LEDGER_BEFORE,
        LAYER_ROWS,
        LEDGER_AFTER,
        233,
    )


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

    grouped, flat_bank = read_bank()
    body_records = []
    for body in sorted(EXPECTED):
        stream, low_signs, ray_checks, three_records, survivors = scalar_reduction(
            body, grouped[body]
        )
        expected = EXPECTED[body]
        cells = fixed_safe_cells(stream, expected["low"])
        require(len(cells) == expected["clean"], (body, len(cells)))
        height = height_reduction(body, stream, cells)
        packets, per_survivor, raw_total = literal_candidates(
            body, stream, survivors, height[2]
        )
        intersections, minimum_bonferroni, minimum_intersection = (
            complete_cell_closure(body, stream, cells, packets)
        )
        endpoint_label, endpoint_cell = expected["endpoint"]
        endpoint_residue = endpoint_cell * endpoint_label % stream.L
        require(14 * endpoint_residue - stream.L == 0,
                (body, "strict-open endpoint control"))
        body_records.append(
            (
                body,
                stream.L,
                stream.high_floor,
                expected["low"],
                expected["deficit"],
                len(grouped[body]),
                low_signs,
                ray_checks,
                three_records,
                survivors,
                len(cells),
                height,
                packets,
                per_survivor,
                raw_total,
                intersections,
                minimum_bonferroni,
                minimum_intersection,
                (endpoint_label, endpoint_cell, endpoint_residue),
            )
        )

    composition = composition_record()
    bank_hash = sha(BANK)
    semantic_packet = (
        SOURCE_3078_SHA256,
        OUTPUT_3078_SHA256,
        SEMANTIC_3078_SHA256,
        SOURCE_3094_SHA256,
        OUTPUT_3094_SHA256,
        SEMANTIC_3094_SHA256,
        SOURCE_3071_SHA256,
        OUTPUT_3071_SHA256,
        SEMANTIC_3071_SHA256,
        bank_hash,
        flat_bank,
        tuple(body_records),
        composition,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 z234 final two-body height and complete-cell closure",
        f"dependency=THM3078_source:{SOURCE_3078_SHA256};output:{OUTPUT_3078_SHA256};semantic:{SEMANTIC_3078_SHA256}",
        f"dependency=THM3094_source:{SOURCE_3094_SHA256};output:{OUTPUT_3094_SHA256};semantic:{SEMANTIC_3094_SHA256}",
        f"ledger_dependency=THM3071_source:{SOURCE_3071_SHA256};output:{OUTPUT_3071_SHA256};semantic:{SEMANTIC_3071_SHA256};atlas:{thm.ATLAS_SHA256}",
        f"mask_bank=rows:{len(flat_bank)};bodies:2;sha256:{bank_hash}",
    ]
    for row in body_records:
        (
            body,
            L,
            H,
            low,
            deficit,
            masks,
            _low_signs,
            ray_checks,
            three_records,
            survivors,
            clean,
            height,
            packets,
            per_survivor,
            raw_total,
            intersections,
            minimum_bonferroni,
            minimum_intersection,
            endpoint,
        ) = row
        minimum_three = min(three_records, key=lambda item: (item[1], item[0]))
        lines.extend(
            [
                f"BODY;E={body};L={L};high={H};masks={masks};low={low};deficit={ftext(deficit)};clean_cells={clean};ray_checks={ray_checks}",
                f"THREE_HIGH;E={body};closed={masks};minimum_gap={ftext(minimum_three[1])};witness_mask={minimum_three[0]};record_sha256={hashlib.sha256(repr(three_records).encode()).hexdigest()}",
                f"TWO_HIGH_UPPER;E={body};raw_expansions={EXPECTED[body]['raw_expansions']};survivors={len(survivors)};same_high_denominator=0;margin={ftext(min(x[7] for x in survivors))}..{ftext(max(x[7] for x in survivors))};record_sha256={hashlib.sha256(repr(survivors).encode()).hexdigest()}",
                f"HEIGHT;E={body};stage1_boundary={ftext(height[0])};smaller_high_ceiling={height[1]};exceptions={','.join(str(x[0]) for x in height[3])};final_box=[{H},{height[2]}]",
                "EXCEPTIONS;E=" + str(body) + ";" + ";".join(
                    f"z={label},load={ftext(mass)},cell={cell},far_cutoff={ftext(cutoff)},margin_at_w=z={ftext(margin)}"
                    for label, mass, cell, cutoff, margin in height[3]
                ),
                f"FINITE;E={body};survivors={len(per_survivor)};raw_box_pairs={raw_total};ordered_literal_candidates={len(packets)};empty_masks={sum(x[-2] == 0 for x in per_survivor)};per_mask_sha256={hashlib.sha256(repr(per_survivor).encode()).hexdigest()};packet_sha256={hashlib.sha256(repr(packets).encode()).hexdigest()}",
            ]
        )
        if body == BODY_A:
            lines.append(
                "EXIT;E=" + str(body)
                + ";mechanism=all_three_scalar_survivors_require_high_denominator_2_but_no_label_in_height_box_has_denominator_2"
            )
        else:
            lines.append(
                f"EXIT;E={body};mechanism=literal_complete_cell_intersection;candidate_pairs={len(intersections)};minimum_bonferroni={minimum_bonferroni[4]}@{minimum_bonferroni[:2]};minimum_actual_intersection={minimum_intersection[5]}@{minimum_intersection[:2]};witness={minimum_intersection[6]};intersection_sha256={hashlib.sha256(repr(intersections).encode()).hexdigest()}"
            )
        lines.append(
            f"STRICT_OPEN;E={body};body_label={endpoint[0]};cell={endpoint[1]};left_endpoint_slack=0;weak_endpoint_inequalities_are_load_bearing"
        )
    lines.extend(
        [
            "projection_logic=a_common_complete_1/L_cell_is_contained_in_S_(E,Z);phi_L(t)=Lt_mod_1_maps_it_surjectively_to_T;therefore_P_(E,Z)=T",
            "THM2941_direction=a_literal_completion_forces_P_(E,Z)_subset_U_A;THM1166_gives_mu(U_A)<=36/91<1;contradiction",
            f"composition=THM3078_closed_377_plus_boundary_4;THM3094_closed_2;THM3098_closes_final_2;all_{LAYER_ROWS}_z1=234_rows_empty",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<=233;next_layer:z1=233_rows:62",
            "scope=projected_k3_necessary_atlas_only;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    args.output.write_text("\n".join(lines) + "\n", newline="\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
