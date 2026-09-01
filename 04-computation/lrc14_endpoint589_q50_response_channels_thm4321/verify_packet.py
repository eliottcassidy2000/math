#!/usr/bin/env python3
"""Hardened semantic verifier for the independent q50 19/95 audit."""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
EXPECTED_FILES = {
    "q50_response_common.hpp",
    "audit_capacity_residual_19_95.cpp",
    "input_endpoint589_failures.csv",
    "input_packing19.csv",
    "input_cover95.csv",
    "audit_O2.out",
    "audit_O3.out",
    "freeze_manifest.py",
    "q50_cover_two_for_one_direct.py",
    "q50_cover_two_for_one_direct.out",
    "q50_cover_two_for_one_direct_opt.out",
    "verify_packet.py",
    "RESULT.md",
    "REPRODUCTION.md",
}
EXPECTED_OUTPUT = [
    "ENDPOINT589_Q50_CAPACITY_19_95_INDEPENDENT_AUDIT_V1",
    "FAILURES 20025 PEEL_COVERED 891 RESIDUAL 19134 RESIDUAL_FNV e67c635fbc71d3b",
    "UPPER COVER_SIZE 95 COVER_FNV fa44f9bfad76cfe7 COVERED 19134 MIN_ACTIVITY_TICKS 1550209054968 MIN_MASK 0a624049",
    "LOWER PACKING_SIZE 19 RESPONDERS_RANK8 172 RESPONDERS_RANK9 17049 MAX_LOAD 1 HIT_FNV 430fdb51e2ee1fa1",
    "STREAM ACTIVE_RANK8 480409 ACTIVE_RANK9 4112383 ACTIVE_FNV b79e52255c5b3522",
    "CONSEQUENCE 19_LE_RESIDUAL_RESPONSE_COVER_NUMBER_LE_95",
    "SCOPE FIXED_Q96_PEEL_AND_Q50_RESPONSE_REPRESENTATION_ONLY",
    "VERDICT PASS",
]
EXPECTED_LOCAL_OUTPUT = [
    "Q50_COVER95_TWO_FOR_ONE_DIRECT_AUDIT_V1",
    "RESIDUAL 19134 COVER 95 PAIRS 4465",
    "AVAILABLE_WIDTH_HIST 0:2;1:25;2:106;3:326;4:642;5:835;6:885;7:704;8:507;9:244;10:131;11:39;12:15;13:3;14:1",
    "ENUMERATED 40468 ACTIVE_FOUND 0 EXCHANGES 0",
    "VERDICT TWO_FOR_ONE_LOCAL_MINIMUM",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fnv_add(state: int, word: int) -> int:
    for byte in range(8):
        state ^= (word >> (8 * byte)) & 0xFF
        state = (state * PRIME) & MASK64
    return state


def hash_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_manifest(packet: Path) -> None:
    entries: dict[str, str] = {}
    for line in (packet / "SHA256SUMS").read_text(encoding="ascii").splitlines():
        digest, name = line.split("  ", 1)
        require(len(digest) == 64 and name not in entries,
                "malformed/duplicate manifest row")
        entries[name] = digest
    require(set(entries) == EXPECTED_FILES, "manifest file set changed")
    for name, digest in entries.items():
        require(hash_file(packet / name) == digest, f"SHA mismatch: {name}")


def build_low_classes() -> tuple[int, list[tuple[int, int]]]:
    speeds = (*POOL, 50, 589)
    grid = 1
    for speed in speeds:
        grid = math.lcm(grid, 14 * speed)
    require(grid == 2827379709554400, "q50 grid changed")
    walls = {0, grid}
    for speed in speeds:
        require(grid % (14 * speed) == 0, "nonintegral wall unit")
        unit = grid // (14 * speed)
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    ordered = sorted(walls)

    def safe(speed: int, left: int, right: int) -> bool:
        residue = speed * (left + right) % (2 * grid)
        return grid <= 7 * residue <= 13 * grid

    classes: dict[int, int] = {}
    for left, right in zip(ordered, ordered[1:]):
        if not (safe(50, left, right) and safe(589, left, right)):
            continue
        failure = 0
        for bit, speed in enumerate(POOL):
            if not safe(speed, left, right):
                failure |= 1 << bit
        if failure.bit_count() <= 9:
            classes[failure] = classes.get(failure, 0) + right - left
    result = sorted(classes.items())
    ledger = OFFSET
    for failure, width in result:
        ledger = fnv_add(fnv_add(ledger, failure), width)
    require(len(ordered) - 1 == 8389 and len(result) == 2383 and
            ledger == 0x88D3EB2D7A477232,
            "q50 wall/class geometry changed")
    return grid, result


def main() -> None:
    packet = Path(__file__).resolve().parent
    parse_manifest(packet)
    require((packet / "audit_O2.out").read_bytes() ==
            (packet / "audit_O3.out").read_bytes(),
            "O2/O3 transcript mismatch")
    require((packet / "audit_O2.out").read_text(encoding="ascii").splitlines()
            == EXPECTED_OUTPUT, "audit transcript semantics changed")
    require((packet / "q50_cover_two_for_one_direct.out").read_bytes() ==
            (packet / "q50_cover_two_for_one_direct_opt.out").read_bytes(),
            "normal/optimized local-exchange transcript mismatch")
    require((packet / "q50_cover_two_for_one_direct.out").read_text(
                encoding="ascii").splitlines() == EXPECTED_LOCAL_OUTPUT,
            "local-exchange transcript semantics changed")
    for name in ("verify_packet.py", "audit_capacity_residual_19_95.cpp",
                 "q50_response_common.hpp", "q50_cover_two_for_one_direct.py"):
        text = (packet / name).read_text(encoding="utf-8")
        forbidden = "ass" + "ert("
        require(forbidden not in text, f"runtime assertion found in {name}")

    q50: list[int] = []
    q96: list[int] = []
    q50_fnv = OFFSET
    with (packet / "input_endpoint589_failures.csv").open(
            newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["q", "r", "body_hex"],
                "failure header changed")
        for row in reader:
            q, r, body = int(row["q"]), int(row["r"]), int(row["body_hex"], 16)
            require(r == 589 and q in (50, 96) and body.bit_count() == 9,
                    "failure row escaped declared universe")
            if q == 50:
                q50.append(body)
                q50_fnv = fnv_add(q50_fnv, body)
            else:
                q96.append(body)
    require(len(q50) == 20025 and len(q96) == 11 and
            len(set(q50)) == len(q50) and q50_fnv == 0xFF421454F02D9099,
            "failure universe changed")

    active_peel = 0x0220932C
    peeled_ledger = OFFSET
    residual_ledger = OFFSET
    residual: list[int] = []
    for index, body in enumerate(q50):
        if active_peel & body == 0:
            peeled_ledger = fnv_add(fnv_add(peeled_ledger, index), body)
        else:
            residual.append(body)
            residual_ledger = fnv_add(fnv_add(residual_ledger, index), body)
    require(len(residual) == 19134 and
            peeled_ledger == 0xF766352D228791ED and
            residual_ledger == 0x0E67C635FBC71D3B,
            "fixed-peel residual changed")
    residual_set = set(residual)

    with (packet / "input_packing19.csv").open(
            newline="", encoding="ascii") as handle:
        packing_rows = list(csv.DictReader(handle))
    packing = [int(row["body_hex"], 16) for row in packing_rows]
    require(len(packing) == 19 and len(set(packing)) == 19 and
            all(body in residual_set for body in packing),
            "packing certificate changed")
    packing_fnv = OFFSET
    for index, row in enumerate(packing_rows):
        require(int(row["solution_index"]) == index,
                "packing order changed")
        packing_fnv = fnv_add(packing_fnv, index)
        packing_fnv = fnv_add(packing_fnv, int(row["core_index"]))
        packing_fnv = fnv_add(packing_fnv, packing[index])
    require(packing_fnv == 0x195708723B22FE7D,
            "packing semantic FNV changed")

    with (packet / "input_cover95.csv").open(
            newline="", encoding="ascii") as handle:
        cover_rows = list(csv.DictReader(handle))
    cover = [int(row["mask_hex"], 16) for row in cover_rows]
    require(len(cover) == 95 and len(set(cover)) == 95 and
            all(mask.bit_count() == 9 for mask in cover),
            "cover mask universe changed")
    cover_fnv = OFFSET
    for mask in cover:
        cover_fnv = fnv_add(cover_fnv, mask)
    require(cover_fnv == 0xFA44F9BFAD76CFE7, "cover FNV changed")

    grid, classes = build_low_classes()

    def ticks(mask: int) -> int:
        mass = sum(width for failure, width in classes
                   if failure & (((1 << 30) - 1) ^ mask) == 0)
        return 63 * mass - 4 * grid

    require(ticks(0x0000B3A5) < 0 and ticks(active_peel) >= 0,
            "fixed-peel activity changed")
    uncovered = set(residual)
    total = 0
    activity = []
    for step, (mask, row) in enumerate(zip(cover, cover_rows), start=1):
        value = ticks(mask)
        require(value >= 0, "cover contains inactive response")
        activity.append((value, mask))
        hit = {body for body in uncovered if mask & body == 0}
        uncovered.difference_update(hit)
        total += len(hit)
        require(int(row["step"]) == step and
                int(row["newly_covered"]) == len(hit) and
                int(row["total_covered"]) == total and
                int(row["residual"]) == len(residual) - total,
                "declared greedy gain changed")
    require(not uncovered and total == 19134 and
            min(activity) == (1550209054968, 0x0A624049),
            "95-mask direct cover replay failed")

    print("ENDPOINT589_Q50_CAPACITY_19_95_INDEPENDENT_PACKET_VERIFY PASS")
    print("q50_failures=20025 peeled=891 residual=19134")
    print("packing=19 full_stream_max_load=1 cover=95 covered=19134")
    print("cover95_two_for_one_exchanges=0")
    print("scope=response_representation_diagnostic_not_endpoint_closure")


if __name__ == "__main__":
    main()
