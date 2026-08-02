#!/usr/bin/env python3
"""Exact complete-cell closure for THM-3087's finite two-high pair bank."""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation/lrc14_j7_k3_z234_two_high_common_torsion_finite_height_thm3087.py"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z234_two_high_common_torsion_finite_height_thm3087.out"
BANK = ROOT / "05-knowledge/results/lrc14_j7_k3_z234_two_high_mask_bank_thm3087.tsv"

SOURCE_SHA256 = "440077a9fda3f9f5cb8898659d0e5006b6e2529b6bbb8c12ced37d7a67e5c69a"
OUTPUT_SHA256 = "cd3059f6a83b2b5cd63176ea359ec76bd3f46fe75e144cef95752f6b462f5bff"
BANK_SHA256 = "b3d424991e816bd5ddc0540554483cd3bc40efbaf44a3f3064197c337ae6ca8d"
EXPECTED = {
    (1, 5, 9, 11, 12, 14): (44062, (22020, 48510), (25066, 48512), 3024),
    (1, 9, 10, 11, 12, 14): (44616, (22308, 48510), (24816, 48512), 2508),
}
FIRST = 19111
LAST = 51480


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha(path: Path) -> str:
    payload = path.read_bytes()
    require(b"\r" not in payload.replace(b"\r\n", b""), f"bare CR: {path}")
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(lf_sha(SOURCE) == SOURCE_SHA256, "THM-3087 source changed")
require(lf_sha(OUTPUT) == OUTPUT_SHA256, "THM-3087 output changed")
require(lf_sha(BANK) == BANK_SHA256, "THM-3087 bank changed")
require("ordered_scalar_candidates:108966498" in OUTPUT.read_text(), "candidate census changed")
thm = load("thm3087_pair_cell_closure", SOURCE)


def exact_safe_counts(cells: np.ndarray, L: int, labels: np.ndarray) -> np.ndarray:
    counts = np.empty(len(labels), dtype=np.int64)
    for start in range(0, len(labels), 128):
        block = labels[start : start + 128]
        residues = (cells[:, None] * block[None, :]) % L
        safe = (14 * residues >= L) & (
            14 * (residues + block[None, :]) <= 13 * L
        )
        counts[start : start + len(block)] = np.count_nonzero(safe, axis=0)
    return counts


def safe_mask(cells: np.ndarray, label: int, L: int) -> np.ndarray:
    residues = (cells * label) % L
    return (14 * residues >= L) & (14 * (residues + label) <= 13 * L)


def main() -> None:
    labels = np.arange(FIRST, LAST + 1, dtype=np.int64)
    records = []
    print("LRC14 z234 two-high complete-cell intersection closure THM3094")
    print(
        f"dependency=THM3087_source:{SOURCE_SHA256};output:{OUTPUT_SHA256};bank:{BANK_SHA256}"
    )
    for body in sorted(EXPECTED):
        clean_expected, minimum_expected, second_expected, slack_expected = EXPECTED[body]
        stream = thm.eng.ray.Stream(body)
        low = thm.INFINITE[body]["low"]
        cells_tuple = thm.fixed_safe_cells(stream, (thm.LEVEL, low))
        require(len(cells_tuple) == clean_expected, (body, "clean count"))
        cells = np.asarray(cells_tuple, dtype=np.int64)
        counts = exact_safe_counts(cells, stream.L, labels)

        order = np.argsort(counts, kind="stable")
        first_index, second_index = map(int, order[:2])
        minimum = (int(counts[first_index]), int(labels[first_index]))
        second = (int(counts[second_index]), int(labels[second_index]))
        require(minimum == minimum_expected, (body, "minimum", minimum))
        require(second == second_expected, (body, "second", second))
        require(np.count_nonzero(counts == minimum[0]) == 1, (body, "minimum uniqueness"))

        slack = minimum[0] + second[0] - len(cells)
        require(slack == slack_expected > 0, (body, "pair slack", slack))

        # A scalar implementation independently checks both extremal columns.
        for count, label in (minimum, second):
            direct = sum(thm.cell_clean(cell, label, stream.L) for cell in cells_tuple)
            require(direct == count, (body, label, direct, count))

        first_mask = safe_mask(cells, minimum[1], stream.L)
        second_mask = safe_mask(cells, second[1], stream.L)
        intersection = np.flatnonzero(first_mask & second_mask)
        require(len(intersection) >= slack > 0, (body, "extremal intersection"))
        witness = int(cells[int(intersection[0])])
        require(
            thm.cell_clean(witness, minimum[1], stream.L)
            and thm.cell_clean(witness, second[1], stream.L),
            (body, "witness"),
        )

        count_digest = hashlib.sha256(
            repr(tuple(zip(map(int, labels), map(int, counts)))).encode("ascii")
        ).hexdigest()
        record = (
            body,
            low,
            stream.L,
            len(cells),
            minimum,
            second,
            slack,
            len(intersection),
            witness,
            count_digest,
        )
        records.append(record)
        print(
            f"BODY;E={body};low={low};L={stream.L};clean={len(cells)};"
            f"minimum={minimum[1]}:{minimum[0]};second={second[1]}:{second[0]};"
            f"distinct_pair_slack={slack};extremal_intersection={len(intersection)};"
            f"witness_cell={witness};count_sha256={count_digest}"
        )

    semantic = hashlib.sha256(repr(tuple(records)).encode("ascii")).hexdigest()
    print(
        "consequence=every_distinct_pair_in_[19111,51480]^2_has_a_common_complete_safe_cell;"
        "P_(E,Z)=T;both_scalar-unbounded_THM3078_rows_close"
    )
    print(
        "scope=108966498_THM3087_candidates_removed;the_other_two_THM3078_rows_remain;"
        "no_projected_cap_or_ledger_or_LRC_claim"
    )
    print(f"semantic_sha256={semantic}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
