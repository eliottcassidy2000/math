#!/usr/bin/env python3
"""Independent HiGHS MIP audit of the 101-bit maximal-response atlas."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("atlas", type=Path)
    args = parser.parse_args()
    raw = args.atlas.read_bytes()
    rows = []
    with args.atlas.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(
                (
                    int(row["mask_hex"], 16),
                    int(row["word_lo"], 16) | (int(row["word_hi"], 16) << 64),
                )
            )
    assert rows
    obligation_count = max(pattern.bit_length() for _, pattern in rows)
    assert 0 < obligation_count <= 128
    matrix = np.zeros((obligation_count, len(rows)), dtype=np.float64)
    for j, (_, pattern) in enumerate(rows):
        for i in range(obligation_count):
            matrix[i, j] = (pattern >> i) & 1
    result = milp(
        c=np.ones(len(rows)),
        integrality=np.ones(len(rows)),
        bounds=Bounds(np.zeros(len(rows)), np.ones(len(rows))),
        constraints=LinearConstraint(
            matrix, np.ones(obligation_count), np.full(obligation_count, np.inf)
        ),
        options={"time_limit": 300.0, "mip_rel_gap": 0.0, "presolve": True},
    )
    assert result.success and result.fun is not None, result.message
    chosen = [rows[i][0] for i, value in enumerate(result.x) if value > 0.5]
    covered = 0
    for i, value in enumerate(result.x):
        if value > 0.5:
            covered |= rows[i][1]
    assert covered == (1 << obligation_count) - 1 and len(chosen) == round(result.fun)
    print(
        f"ATLAS_SHA256 {hashlib.sha256(raw).hexdigest()} ROWS {len(rows)} "
        f"OBLIGATIONS {obligation_count}"
    )
    print(f"STATUS {result.status} SUCCESS {result.success} OBJECTIVE {result.fun:.0f}")
    print(f"MIP_GAP {getattr(result, 'mip_gap', None)} NODE_COUNT {getattr(result, 'mip_node_count', None)}")
    print("MASKS " + " ".join(f"{mask:08x}" for mask in chosen))
    print("VERDICT PASS INDEPENDENT_HIGHS_INTEGER_SET_COVER")


if __name__ == "__main__":
    main()
