#!/usr/bin/env python3
"""Solve and exactly verify the 84-obligation k=350 retention quotient."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csc_matrix


def main() -> None:
    path = Path(sys.argv[1])
    responses: list[int] = []
    masks: list[int] = []
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["w0", "w1", "count", "least_mask"]:
            raise RuntimeError("atlas header changed")
        for row in reader:
            responses.append(int(row["w0"], 16) | (int(row["w1"], 16) << 64))
            masks.append(int(row["least_mask"], 16))
    if len(responses) != 53 or len(set(responses)) != 53:
        raise RuntimeError("response universe changed")
    rows: list[int] = []
    cols: list[int] = []
    for col, response in enumerate(responses):
        bits = response
        while bits:
            least = bits & -bits
            rows.append(least.bit_length() - 1)
            cols.append(col)
            bits ^= least
    incidence = csc_matrix(
        (np.ones(len(rows)), (np.asarray(rows), np.asarray(cols))),
        shape=(84, len(responses)),
    )
    cover = milp(
        c=np.ones(len(responses)),
        integrality=np.ones(len(responses)),
        bounds=Bounds(np.zeros(len(responses)), np.ones(len(responses))),
        constraints=LinearConstraint(incidence, np.ones(84), np.full(84, np.inf)),
        options={"mip_rel_gap": 0.0},
    )
    if not cover.success:
        raise RuntimeError(f"cover MILP failed: {cover.message}")
    chosen = [i for i, value in enumerate(cover.x) if value > 0.5]
    union = 0
    for i in chosen:
        union |= responses[i]
    if union != (1 << 84) - 1 or len(chosen) != round(cover.fun):
        raise RuntimeError("cover failed exact integer replay")

    packing = milp(
        c=-np.ones(84),
        integrality=np.ones(84),
        bounds=Bounds(np.zeros(84), np.ones(84)),
        constraints=LinearConstraint(
            incidence.T, np.full(len(responses), -np.inf), np.ones(len(responses))
        ),
        options={"mip_rel_gap": 0.0},
    )
    if not packing.success:
        raise RuntimeError(f"packing MILP failed: {packing.message}")
    support = [i for i, value in enumerate(packing.x) if value > 0.5]
    for response in responses:
        if sum((response >> i) & 1 for i in support) > 1:
            raise RuntimeError("packing failed exact integer replay")

    print("LRC14_K350_RESIDUAL_RETENTION_SOLVER_V1")
    print(
        f"TYPES {len(responses)} OBLIGATIONS 84 COVER_MIN {len(chosen)} "
        f"PACKING_MAX {len(support)} COVER_NODES {cover.mip_node_count} "
        f"PACKING_NODES {packing.mip_node_count}"
    )
    for i in chosen:
        w0 = responses[i] & ((1 << 64) - 1)
        w1 = responses[i] >> 64
        print(f"RETAIN {masks[i]:08x} RESPONSE {w0:016x}:{w1:016x}")
    print("PACKING " + " ".join(map(str, support)))
    print(f"COVER_UNION {union:x}")
    print(
        "VERDICT " + ("PASS" if len(chosen) == len(support) else "NO_MATCHING_PACKING")
    )
    if len(chosen) != len(support):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
