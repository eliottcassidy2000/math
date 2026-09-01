#!/usr/bin/env python3
"""Solve and exactly certify the endpoint-595 pair-tagged response covers."""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp
from scipy.sparse import csc_matrix


ROW_SLICES = {
    "all": tuple(range(145)),
    "96,595": tuple(range(116)),
    "100,595": tuple(range(116, 129)),
    "210,595": tuple(range(129, 145)),
}
BIT_POSITIONS = tuple(range(116)) + tuple(range(128, 141)) + tuple(range(192, 208))
FROZEN_MIXED_ALL_DUAL = {
    4: 361, 11: 49, 12: 153, 13: 52, 16: 664, 17: 714, 26: 202,
    27: 1, 28: 184, 32: 345, 33: 168, 36: 157, 38: 91, 42: 160,
    43: 449, 49: 119, 50: 295, 52: 386, 53: 440, 55: 167,
    64: 254, 66: 128, 69: 649, 73: 160, 81: 312, 91: 553,
    94: 31, 97: 626, 98: 304, 101: 715, 107: 49, 117: 649,
    118: 354, 120: 265, 123: 153, 124: 136, 125: 135, 127: 729,
    128: 97, 129: 128, 130: 8, 132: 248, 133: 176, 135: 537,
    136: 112, 137: 312, 138: 394, 139: 649, 141: 706, 142: 304,
}


@dataclass(frozen=True)
class Source:
    global_reply: int
    mask: int
    rank: int


def read_atlas(path: Path):
    rows = []
    with path.open(newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle):
            words = [int(row[f"w{i}"], 16) for i in range(4)]
            reply = words[0] | (words[1] << 64) | (words[2] << 128) | (words[3] << 192)
            rows.append(
                (
                    reply,
                    int(row["count8"]),
                    int(row["least8"], 16) if row["least8"] else None,
                    int(row["count9"]),
                    int(row["least9"], 16) if row["least9"] else None,
                )
            )
    return rows


def project(reply: int, obligations: tuple[int, ...]) -> int:
    answer = 0
    for local, global_index in enumerate(obligations):
        if reply >> BIT_POSITIONS[global_index] & 1:
            answer |= 1 << local
    return answer


def allowed_sources(rows, mode: str, obligations: tuple[int, ...]):
    unique: dict[int, Source] = {}
    for reply, count8, least8, count9, least9 in rows:
        candidates = []
        if mode in ("mixed", "rank8") and count8:
            candidates.append((least8, 8))
        if mode == "mixed" and count9:
            candidates.append((least9, 9))
        if not candidates:
            continue
        projected = project(reply, obligations)
        if projected == 0:
            continue
        mask, rank = min(candidates)
        source = Source(reply, mask, rank)
        if projected not in unique or mask < unique[projected].mask:
            unique[projected] = source
    return unique


def maximal_responses(unique: dict[int, Source]):
    maximal: list[int] = []
    for reply in sorted(unique, key=lambda value: (-value.bit_count(), value)):
        if any(reply & keep == reply for keep in maximal):
            continue
        maximal.append(reply)
    return maximal


def incidence(responses: list[int], obligation_count: int):
    row_index = []
    column_index = []
    for column, reply in enumerate(responses):
        bits = reply
        while bits:
            least = bits & -bits
            row_index.append(least.bit_length() - 1)
            column_index.append(column)
            bits ^= least
    data = np.ones(len(row_index), dtype=float)
    return csc_matrix(
        (data, (np.asarray(row_index), np.asarray(column_index))),
        shape=(obligation_count, len(responses)),
    )


def exact_cover(responses: list[int], obligation_count: int):
    matrix = incidence(responses, obligation_count)
    result = milp(
        c=np.ones(len(responses)),
        integrality=np.ones(len(responses)),
        bounds=Bounds(np.zeros(len(responses)), np.ones(len(responses))),
        constraints=LinearConstraint(
            matrix, np.ones(obligation_count), np.full(obligation_count, np.inf)
        ),
        options={"mip_rel_gap": 0.0, "time_limit": 600},
    )
    if not result.success:
        raise RuntimeError(f"cover MILP failed: {result.message}")
    chosen = [index for index, value in enumerate(result.x) if value > 0.5]
    union = 0
    for index in chosen:
        union |= responses[index]
    if union != (1 << obligation_count) - 1:
        raise AssertionError("floating cover did not verify exactly")
    optimum = len(chosen)
    if abs(result.fun - optimum) > 1e-6:
        raise AssertionError("cover objective is not integral")
    return optimum, chosen, result


def integer_dual_certificate(
    all_responses: list[int], obligation_count: int, lower: int
):
    matrix = incidence(all_responses, obligation_count)
    # Seek small-denominator nonnegative integral weights z_i/D with
    # sum_{i in response} z_i <= D.  Any total > D*(lower-1) certifies lower.
    for denominator in range(1, 13):
        result = milp(
            c=-np.ones(obligation_count),
            integrality=np.ones(obligation_count),
            bounds=Bounds(
                np.zeros(obligation_count),
                np.full(obligation_count, denominator),
            ),
            constraints=LinearConstraint(
                matrix.T,
                np.full(len(all_responses), -np.inf),
                np.full(len(all_responses), denominator),
            ),
            options={"mip_rel_gap": 0.0, "time_limit": 600},
        )
        if not result.success:
            continue
        weights = np.rint(result.x).astype(np.int64)
        loads = np.asarray(matrix.T @ weights).reshape(-1)
        if np.any(loads > denominator):
            raise AssertionError("dual load failed exact integer check")
        total = int(weights.sum())
        if total > denominator * (lower - 1):
            support = [(index, int(value)) for index, value in enumerate(weights) if value]
            return denominator, total, int(loads.max(initial=0)), support, result
    # Some exact LP vertices have a larger but still transparent common
    # denominator.  Recover it from HiGHS, then verify the scaled certificate
    # with integer arithmetic against every complete realized response type.
    lp = linprog(
        -np.ones(obligation_count),
        A_ub=matrix.T,
        b_ub=np.ones(len(all_responses)),
        bounds=(0, None),
        method="highs",
    )
    if lp.success and -lp.fun > lower - 1 + 1e-9:
        for denominator in range(13, 10001):
            scaled = lp.x * denominator
            weights = np.rint(scaled).astype(np.int64)
            if np.max(np.abs(scaled - weights), initial=0.0) > 1e-7:
                continue
            loads = np.asarray(matrix.T @ weights).reshape(-1)
            total = int(weights.sum())
            if np.any(loads > denominator) or total <= denominator * (lower - 1):
                continue
            support = [(index, int(value)) for index, value in enumerate(weights) if value]
            return denominator, total, int(loads.max(initial=0)), support, lp
    return None


def frozen_mixed_all_certificate(all_responses: list[int]):
    denominator = 1667
    weights = np.zeros(145, dtype=np.int64)
    for index, value in FROZEN_MIXED_ALL_DUAL.items():
        weights[index] = value
    matrix = incidence(all_responses, 145)
    loads = np.asarray(matrix.T @ weights).reshape(-1)
    total = int(weights.sum())
    if total != 15030 or np.any(loads > denominator) or int(loads.max()) != denominator:
        raise AssertionError("frozen D=1667 dual certificate failed exact check")
    support = sorted(FROZEN_MIXED_ALL_DUAL.items())
    return (
        denominator,
        total,
        int(loads.max()),
        support,
        SimpleNamespace(mip_node_count="FROZEN_EXACT_INTEGER_CHECK"),
    )


def response_hex(reply: int) -> str:
    return ":".join(f"{(reply >> (64 * word)) & ((1 << 64) - 1):016x}" for word in range(4))


def run(path: Path):
    rows = read_atlas(path)
    print("LRC14_ENDPOINT595_PAIR_TAGGED_SET_COVER_SOLVER_V1")
    print(f"ATLAS_ROWS {len(rows)}")
    for mode in ("mixed", "rank8"):
        for scope, obligations in ROW_SLICES.items():
            unique = allowed_sources(rows, mode, obligations)
            maximal = maximal_responses(unique)
            optimum, chosen, cover_result = exact_cover(maximal, len(obligations))
            if mode == "mixed" and scope == "all":
                dual = frozen_mixed_all_certificate(list(unique))
            elif mode == "rank8" and scope == "all":
                # The rank-eight LP optimum is 167/8 < 21, so an obligation-
                # weight dual cannot certify the integer scout value 22.
                dual = None
            else:
                dual = integer_dual_certificate(list(unique), len(obligations), optimum)
            if dual is None:
                dual_text = "DUAL_CERTIFICATE none_denominator_le_12"
                support = []
            else:
                denominator, total, max_load, support, dual_result = dual
                dual_text = (
                    f"DUAL_DENOMINATOR {denominator} DUAL_TOTAL {total} "
                    f"DUAL_VALUE {total}/{denominator} DUAL_MAX_LOAD {max_load} "
                    f"DUAL_MIP_NODES {getattr(dual_result, 'mip_node_count', 'LP')}"
                )
            optimum_label = (
                "MIP_OPTIMUM_SCOUT" if mode == "rank8" and scope == "all"
                else "EXACT_MINIMUM"
            )
            print(
                f"CASE MODE {mode} SCOPE {scope} OBLIGATIONS {len(obligations)} "
                f"UNIQUE_TYPES {len(unique)} MAXIMAL_TYPES {len(maximal)} "
                f"{optimum_label} {optimum} COVER_MIP_NODES {cover_result.mip_node_count} "
                + dual_text
            )
            for ordinal, index in enumerate(chosen):
                reply = maximal[index]
                source = unique[reply]
                print(
                    f"COVER {mode} {scope} {ordinal} MASK {source.mask:08x} "
                    f"RANK {source.rank} PROJECTED {reply:x} "
                    f"GLOBAL {response_hex(source.global_reply)}"
                )
            print(
                "DUAL_SUPPORT "
                + mode
                + " "
                + scope
                + " "
                + " ".join(f"{index}:{value}" for index, value in support)
            )
    print("SCOPE EXACT_FINITE_SET_COVER_ON_FROZEN_COMPLETE_RESPONSE_ATLAS")
    print("VERDICT PASS")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: solve_cover.py RESPONSE_ATLAS.csv")
    run(Path(sys.argv[1]))
