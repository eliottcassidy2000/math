#!/usr/bin/env python3
"""Find and bitwise-verify a small cover of the 100 endpoint-590 failures."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csc_matrix


FULL = (1 << 100) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_signatures(path: Path) -> list[tuple[int, int, int]]:
    records: list[tuple[int, int, int]] = []
    seen: set[int] = set()
    with path.open(newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle):
            reply = int(row["w0"], 16) | (int(row["w1"], 16) << 64)
            require(reply != 0 and reply < 1 << 100,
                    "response escaped the 100-obligation universe")
            require(reply not in seen, "duplicate response signature")
            seen.add(reply)
            candidates: list[tuple[int, int]] = []
            if row["least8"]:
                candidates.append((int(row["least8"], 16), 8))
            if row["least9"]:
                candidates.append((int(row["least9"], 16), 9))
            require(bool(candidates), "signature has no representative")
            mask, rank = min(candidates)
            require(mask.bit_count() == rank, "representative rank changed")
            records.append((reply, mask, rank))
    require(len(records) == 14368, "signature count changed")
    return records


def build_matrix(records: list[tuple[int, int, int]]) -> csc_matrix:
    rows: list[int] = []
    columns: list[int] = []
    for column, (reply, _, _) in enumerate(records):
        work = reply
        while work:
            least = work & -work
            rows.append(least.bit_length() - 1)
            columns.append(column)
            work ^= least
    return csc_matrix(
        (np.ones(len(rows)), (np.asarray(rows), np.asarray(columns))),
        shape=(100, len(records)),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signatures", type=Path, required=True)
    parser.add_argument("--cover", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    records = read_signatures(args.signatures)
    matrix = build_matrix(records)
    result = milp(
        c=np.ones(len(records)),
        integrality=np.ones(len(records)),
        bounds=Bounds(np.zeros(len(records)), np.ones(len(records))),
        constraints=LinearConstraint(
            matrix, np.ones(100), np.full(100, np.inf)
        ),
        options={"mip_rel_gap": 0.0, "time_limit": 300, "presolve": True},
    )
    require(result.x is not None, f"MILP produced no incumbent: {result.message}")
    chosen = [index for index, value in enumerate(result.x) if value > 0.5]
    union = 0
    for index in chosen:
        union |= records[index][0]
    require(union == FULL, "chosen response family does not cover exactly")
    require(abs(float(result.fun) - len(chosen)) < 1e-7,
            "MILP objective differs from exact chosen count")

    args.cover.parent.mkdir(parents=True, exist_ok=True)
    with args.cover.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["mask_hex", "rank", "weight", "w0", "w1"])
        for index in chosen:
            reply, mask, rank = records[index]
            writer.writerow([
                f"{mask:08x}", rank, reply.bit_count(),
                f"{reply & ((1 << 64) - 1):016x}",
                f"{reply >> 64:016x}",
            ])

    summary = {
        "status": int(result.status),
        "success": bool(result.success),
        "message": str(result.message),
        "signature_count": len(records),
        "cover_size": len(chosen),
        "exact_union_hex": f"{union:025x}",
        "mip_dual_bound": (
            None if getattr(result, "mip_dual_bound", None) is None
            else float(result.mip_dual_bound)
        ),
        "mip_gap": (
            None if getattr(result, "mip_gap", None) is None
            else float(result.mip_gap)
        ),
        "mip_node_count": getattr(result, "mip_node_count", None),
        "scope": (
            "chosen_cover_bitwise_exact; minimum_only_solver_reported; "
            "no independently certified integer lower bound"
        ),
    }
    args.summary.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )
    print("ENDPOINT590_RESPONSE_COVER_SOLVER_V1")
    print(f"SIGNATURES {len(records)} COVER_SIZE {len(chosen)}")
    print(f"SOLVER_SUCCESS {int(result.success)} STATUS {result.status}")
    print(f"DUAL_BOUND {summary['mip_dual_bound']} GAP {summary['mip_gap']}")
    print("SCOPE BITWISE_EXACT_CHOSEN_COVER_MINIMUM_ONLY_SOLVER_REPORTED")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
