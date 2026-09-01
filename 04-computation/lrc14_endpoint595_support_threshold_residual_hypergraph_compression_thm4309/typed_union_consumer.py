#!/usr/bin/env python3
"""Typed set consumer for the repaired complete endpoint-595 layer."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

Pair = tuple[int, int]
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(ok: bool, message: str):
    if not ok:
        raise RuntimeError(message)


def payload(rows: set[Pair]) -> bytes:
    return "".join(f"{q},{r}\n" for q, r in sorted(rows)).encode("ascii")


def read_pairs(path: Path) -> set[Pair]:
    raw = path.read_bytes()
    rows = [tuple(map(int, line.split(b","))) for line in raw.splitlines()]
    require(rows == sorted(set(rows)), f"unordered/duplicate pair file: {path}")
    require(all(0 < q < r for q, r in rows), f"invalid pair file: {path}")
    answer = set(rows)
    require(raw == payload(answer), f"noncanonical pair bytes: {path}")
    return answer


def fnv(rows: set[Pair]) -> str:
    state = FNV_OFFSET
    for q, r in sorted(rows):
        for value in (q, r):
            for byte in value.to_bytes(8, "little"):
                state = ((state ^ byte) * FNV_PRIME) & MASK64
    return f"{state:016x}"


def identity(rows: set[Pair]):
    raw = payload(rows)
    return len(rows), fnv(rows), hashlib.sha256(raw).hexdigest()


def describe(label: str, rows: set[Pair]):
    count, fingerprint, digest = identity(rows)
    print(
        f"{label} rows={count} fnv={fingerprint} sha256={digest} "
        f"max_endpoint={max((r for _, r in rows), default=0)}"
    )


def read_closed_audit(path: Path) -> set[Pair]:
    expected = [
        "q", "r", "active", "active_fnv", "active_joint", "active_nonjoint",
        "exposed", "exposed_fnv", "minimum_hits", "maximum_hits", "failures",
        "failure_fnv",
    ]
    rows = set()
    previous = None
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == expected, "pair-audit header changed")
        for record in reader:
            pair = int(record["q"]), int(record["r"])
            require(previous is None or previous < pair, "pair audit unordered")
            require(pair not in rows and int(record["failures"]) == 0,
                    "duplicate/failing repaired pair")
            require(int(record["active"]) ==
                    int(record["active_joint"]) + int(record["active_nonjoint"]),
                    "active partition failed")
            require(int(record["minimum_hits"]) >= 1, "closed row has no witness")
            rows.add(pair)
            previous = pair
    require(path.read_bytes().endswith(b"\n"), "pair audit lacks LF")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("universe", type=Path)
    parser.add_argument("old_union", type=Path)
    parser.add_argument("old_residual", type=Path)
    parser.add_argument("pair_audit", type=Path)
    parser.add_argument("union_out", type=Path)
    parser.add_argument("residual_out", type=Path)
    parser.add_argument("top_out", type=Path)
    args = parser.parse_args()

    universe = read_pairs(args.universe)
    old_union = read_pairs(args.old_union)
    old_residual = read_pairs(args.old_residual)
    require(identity(universe) == (
        22647, "df5374d4aca67677",
        "14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317",
    ), "universe changed")
    require(identity(old_union) == (
        1624, "11414a33ab91fef6",
        "ef9102553cd030f67ab1bdb7d6965c3efaf4b0d8aa85daa1092354c9703caf26",
    ), "THM-4300 union changed")
    require(identity(old_residual) == (
        21023, "e93e8089e9dc58c0",
        "ce215cb53a742e5e0d0d4f16e344687da9583ce5dc1417c7ea70399fb0bf70ba",
    ) and universe - old_union == old_residual, "THM-4300 partition changed")

    k596 = {pair for pair in universe if pair[1] >= 596}
    thm4302_union = old_union | k596
    thm4302_residual = universe - thm4302_union
    require(identity(thm4302_union) == (
        1633, "b1c8ecf1dd4a71c5",
        "28084fded429f407188e471183d5645d28eded621967e10e386261bbe52844c0",
    ), "THM-4302 union changed")
    require(identity(thm4302_residual) == (
        21014, "7da11cd038486887",
        "2a3ee951deb5b7cfbb4b86aabd4058c8073aae713b42afebabc15e3159deb3b6",
    ), "THM-4302 residual changed")
    maximum = max(r for _, r in thm4302_residual)
    top28 = {pair for pair in thm4302_residual if pair[1] == maximum}
    require(maximum == 595 and identity(top28) == (
        28, "47981ce64825ef2a",
        "c607dab04e4f6849a2226f518771e43b1301d4fc582b47bfa5752c4643c93702",
    ), "THM-4302 top layer changed")

    audited_target = read_closed_audit(args.pair_audit)
    target391 = k596 | top28
    require(audited_target == target391,
            "raw replay is not exactly the complete 391-row target")
    audited = audited_target & top28
    typed_union = thm4302_union | audited
    residual = universe - typed_union
    require(identity(typed_union) == (
        1661, "5bdd2ebf09e9404a",
        "de00493a80ca88eb4ed802be00fce19967f0978508439bd07afc7393708a4b62",
    ), "new typed union changed")
    require(identity(residual) == (
        20986, "606bf18913a49a14",
        "67561f7f0c5c3a32155811e9978b42b2393c10ea8964387eb712cea6a6683f50",
    ), "new residual changed")
    require(typed_union.isdisjoint(residual) and typed_union | residual == universe,
            "new partition failed")
    residual_maximum = max(r for _, r in residual)
    residual_top = {pair for pair in residual if pair[1] == residual_maximum}
    require(residual_maximum == 594 and identity(residual_top) == (
        25, "cce015c81f7121d9",
        "920638d6fb23a8f6492d34cf50e7dc247c2eddfe7ba3f2088c59155e1a56167e",
    ), "new residual maximum layer changed")

    args.union_out.write_bytes(payload(typed_union))
    args.residual_out.write_bytes(payload(residual))
    args.top_out.write_bytes(payload(residual_top))
    print("LRC14_ENDPOINT595_REPAIRED_TYPED_UNION_CONSUMER_V1")
    describe("UNIVERSE", universe)
    describe("THM4302_TYPED_UNION", thm4302_union)
    describe("AUDITED_ENDPOINT595_LAYER", audited)
    describe("TYPED_UNION", typed_union)
    describe("FINAL_RESIDUAL", residual)
    describe("NEW_TOP_LAYER", residual_top)
    print("NEW_TOP_ROWS " + " ".join(f"({q},{r})" for q, r in sorted(residual_top)))
    print("TYPE_GUARD row_consequences_only_fixed_carrier_not_physical_entry")
    print("SCOPE FINITE_EXACT_TYPED_UNION_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
