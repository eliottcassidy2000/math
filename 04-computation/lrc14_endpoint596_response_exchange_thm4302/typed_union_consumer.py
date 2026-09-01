#!/usr/bin/env python3
"""Independent typed row-set consequence consumer for THM-4302."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

Pair = tuple[int, int]
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def read_pairs(path: Path) -> set[Pair]:
    raw = path.read_bytes()
    rows = [tuple(map(int, line.split(b","))) for line in raw.splitlines()]
    require(rows == sorted(set(rows)), f"unordered or duplicate: {path}")
    require(all(0 < q < r for q, r in rows), f"invalid pair: {path}")
    require(raw == payload(set(rows)), f"noncanonical bytes: {path}")
    return set(rows)


def payload(rows: set[Pair]) -> bytes:
    return "".join(f"{q},{r}\n" for q, r in sorted(rows)).encode("ascii")


def fnv(rows: set[Pair]) -> str:
    state = FNV_OFFSET
    for q, r in sorted(rows):
        for value in (q, r):
            for byte in value.to_bytes(8, "little"):
                state = ((state ^ byte) * FNV_PRIME) & MASK64
    return f"{state:016x}"


def identity(rows: set[Pair]) -> tuple[int, str, str]:
    raw = payload(rows)
    return len(rows), fnv(rows), hashlib.sha256(raw).hexdigest()


def describe(label: str, rows: set[Pair]) -> None:
    count, fingerprint, digest = identity(rows)
    maximum = max((r for _, r in rows), default=0)
    print(
        f"{label} rows={count} fnv={fingerprint} sha256={digest} "
        f"max_endpoint={maximum}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("universe", type=Path)
    parser.add_argument("old_union", type=Path)
    parser.add_argument("old_residual", type=Path)
    args = parser.parse_args()

    universe = read_pairs(args.universe)
    old_union = read_pairs(args.old_union)
    old_residual = read_pairs(args.old_residual)
    require(identity(universe) == (
        22647,
        "df5374d4aca67677",
        "14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317",
    ), "universe changed")
    require(identity(old_union) == (
        1624,
        "11414a33ab91fef6",
        "ef9102553cd030f67ab1bdb7d6965c3efaf4b0d8aa85daa1092354c9703caf26",
    ), "THM-4300 union changed")
    require(identity(old_residual) == (
        21023,
        "e93e8089e9dc58c0",
        "ce215cb53a742e5e0d0d4f16e344687da9583ce5dc1417c7ea70399fb0bf70ba",
    ), "THM-4300 residual changed")
    require(universe - old_union == old_residual,
            "THM-4300 row partition changed")

    k596 = {pair for pair in universe if pair[1] >= 596}
    k597 = {pair for pair in universe if pair[1] >= 597}
    old_overlap = old_union & k596
    new_rows = k596 - old_union
    typed_union = old_union | k596
    residual = universe - typed_union
    require(identity(k596) == (
        363,
        "fbf4e6bd7a593649",
        "6fa377dd16aaa12e62fb5e4c6ec36ba30701c9055a1d503b227cb76546f13960",
    ), "K596 changed")
    require(identity(old_overlap) == (
        354,
        "33b0069ca994b786",
        "e653857c9bfe1ef50e7724cfad05232b3695f88534ca289844a46d914ca52df5",
    ), "old/K596 overlap changed")
    require(old_overlap == k597,
            "old/K596 overlap is not the complete K597 prefix")
    expected_new_rows = {
        (96, 596), (100, 596), (192, 596), (210, 596), (256, 596),
        (260, 596), (294, 596), (306, 596), (384, 596),
    }
    require(new_rows == expected_new_rows,
            "strictly new rows are not the complete layer-596 boundary")
    require(identity(new_rows) == (
        9,
        "1f9be79f9df5bf16",
        "1632cdf580e32ceaf627b0b9bb9d2580fabe18f5000f6a5d4c6b2d1f86b830c3",
    ), "new layer-596 rows changed")
    require(identity(typed_union) == (
        1633,
        "b1c8ecf1dd4a71c5",
        "28084fded429f407188e471183d5645d28eded621967e10e386261bbe52844c0",
    ), "new typed union changed")
    require(identity(residual) == (
        21014,
        "7da11cd038486887",
        "2a3ee951deb5b7cfbb4b86aabd4058c8073aae713b42afebabc15e3159deb3b6",
    ), "new residual changed")
    require(typed_union.isdisjoint(residual) and typed_union | residual == universe,
            "new row partition failed")

    maximum = max(r for _, r in residual)
    top = {pair for pair in residual if pair[1] == maximum}
    expected_top = {
        (96, 595), (100, 595), (147, 595), (186, 595), (192, 595),
        (206, 595), (210, 595), (220, 595), (244, 595), (256, 595),
        (260, 595), (265, 595), (294, 595), (296, 595), (306, 595),
        (320, 595), (332, 595), (338, 595), (346, 595), (366, 595),
        (370, 595), (372, 595), (384, 595), (416, 595), (420, 595),
        (440, 595), (512, 595), (520, 595),
    }
    require(maximum == 595 and top == expected_top,
            "new residual boundary changed")
    require(identity(top) == (
        28,
        "47981ce64825ef2a",
        "c607dab04e4f6849a2226f518771e43b1301d4fc582b47bfa5752c4643c93702",
    ), "top layer identity changed")

    print("THM4302_TYPED_UNION_CONSUMER_V1")
    describe("UNIVERSE", universe)
    describe("OLD_TYPED_UNION", old_union)
    describe("K596_COMPLETE_PREFIX", k596)
    describe("OLD_K596_OVERLAP", old_overlap)
    describe("STRICTLY_NEW_ROWS", new_rows)
    describe("TYPED_UNION", typed_union)
    describe("FINAL_RESIDUAL", residual)
    describe("FINAL_TOP_LAYER", top)
    print("FINAL_TOP_ROWS " + " ".join(f"({q},{r})" for q, r in sorted(top)))
    print("TYPE_GUARD row_consequences_only_carrier_not_physical_entry")
    print("SCOPE FINITE_EXACT_TYPED_UNION_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
