#!/usr/bin/env python3
"""Typed row-set consumer for the THM-4305 five-mask endpoint-595 repair."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

Pair = tuple[int, int]
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def payload(rows: set[Pair]) -> bytes:
    return "".join(f"{q},{r}\n" for q, r in sorted(rows)).encode("ascii")


def read_pairs(path: Path) -> set[Pair]:
    raw = path.read_bytes()
    rows = [tuple(map(int, line.split(b","))) for line in raw.splitlines()]
    require(rows == sorted(set(rows)), f"unordered or duplicate: {path}")
    require(all(0 < q < r for q, r in rows), f"invalid pair: {path}")
    answer = set(rows)
    require(raw == payload(answer), f"noncanonical bytes: {path}")
    return answer


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


def read_pair_audit(path: Path) -> dict[Pair, tuple[int, str]]:
    expected_header = [
        "q", "r", "active", "active_fnv", "active_joint",
        "active_nonjoint", "exposed", "exposed_fnv", "minimum_hits",
        "maximum_hits", "failures", "failure_fnv",
    ]
    answer: dict[Pair, tuple[int, str]] = {}
    previous: Pair | None = None
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == expected_header, "pair-audit header changed")
        for row in reader:
            pair = (int(row["q"]), int(row["r"]))
            require(0 < pair[0] < pair[1], "invalid audited pair")
            require(previous is None or previous < pair,
                    "pair-audit order changed")
            require(pair not in answer, "duplicate audited pair")
            active = int(row["active"])
            active_joint = int(row["active_joint"])
            active_nonjoint = int(row["active_nonjoint"])
            exposed = int(row["exposed"])
            minimum_hits = int(row["minimum_hits"])
            maximum_hits = int(row["maximum_hits"])
            failures = int(row["failures"])
            failure_fnv = f"{int(row['failure_fnv'], 16):016x}"
            require(active == active_joint + active_nonjoint,
                    "active partition changed")
            require(exposed >= failures and 0 <= minimum_hits <= maximum_hits,
                    "invalid raw-audit statistics")
            require(len(row["active_fnv"]) <= 16 and
                    len(row["exposed_fnv"]) <= 16 and
                    len(row["failure_fnv"]) <= 16,
                    "non-u64 audit fingerprint")
            answer[pair] = failures, failure_fnv
            previous = pair
    require(path.read_bytes().endswith(b"\n"), "pair audit lacks final LF")
    return answer


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("universe", type=Path)
    parser.add_argument("old_union", type=Path)
    parser.add_argument("old_residual", type=Path)
    parser.add_argument("pair_audit", type=Path)
    parser.add_argument("five_cover", type=Path)
    parser.add_argument("output_union", type=Path)
    parser.add_argument("output_residual", type=Path)
    parser.add_argument("output_top", type=Path)
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
            "THM-4300 partition changed")

    k596 = {pair for pair in universe if pair[1] >= 596}
    thm4302_union = old_union | k596
    thm4302_residual = universe - thm4302_union
    require(identity(thm4302_union) == (
        1633,
        "b1c8ecf1dd4a71c5",
        "28084fded429f407188e471183d5645d28eded621967e10e386261bbe52844c0",
    ), "THM-4302 union changed")
    require(identity(thm4302_residual) == (
        21014,
        "7da11cd038486887",
        "2a3ee951deb5b7cfbb4b86aabd4058c8073aae713b42afebabc15e3159deb3b6",
    ), "THM-4302 residual changed")
    maximum = max(r for _, r in thm4302_residual)
    top28 = {pair for pair in thm4302_residual if pair[1] == maximum}
    require(identity(top28) == (
        28,
        "47981ce64825ef2a",
        "c607dab04e4f6849a2226f518771e43b1301d4fc582b47bfa5752c4643c93702",
    ) and maximum == 595, "THM-4302 top layer changed")

    audit = read_pair_audit(args.pair_audit)
    require(set(audit) == top28,
            "raw audit is not the complete derived top layer")
    expected_failures = {
        (96, 595): (116, "fedacdbff3f31981"),
        (100, 595): (13, "3ac9ac8b4b9ad93f"),
        (210, 595): (16, "a6a226f12c168d3a"),
    }
    survivors = {pair for pair, value in audit.items() if value[0] != 0}
    closed = set(audit) - survivors
    require({pair: audit[pair] for pair in survivors} == expected_failures,
            "endpoint-595 failure identities changed")
    require(identity(closed) == (
        25,
        "1ad4f3c2ab6ea09d",
        "6839c30eafed0d6d73dcd6a0eff6ed7b4751798ff88dce33bcf282397ae247d5",
    ), "closed-row identity changed")
    require(identity(survivors) == (
        3,
        "9853e590efc73022",
        "5f521627485bc829bc6ec2c9d25ce3fb899cfe8d9b04a8f96a225da7e9d84dde",
    ), "surviving-row identity changed")

    typed_union = thm4302_union | closed
    residual = universe - typed_union
    require(identity(typed_union) == (
        1658,
        "43317f1aee06e8bd",
        "bfeb739dcad61dd42bdd9a8b295b6058f3ecee5cc30acd64c469e3b8132393c7",
    ), "THM-4303 union changed")
    require(identity(residual) == (
        20989,
        "b0fbaa28440a118f",
        "bbf2dbba58a5492f6e5d136f14940c3bac8b3ddea604b19e3e2b926abb8bad00",
    ), "THM-4303 residual changed")
    require(typed_union.isdisjoint(residual) and
            typed_union | residual == universe,
            "THM-4303 partition failed")
    residual_maximum = max(r for _, r in residual)
    residual_top = {pair for pair in residual if pair[1] == residual_maximum}
    require(residual_maximum == 595 and residual_top == survivors,
            "THM-4303 residual boundary changed")

    cover_text = args.five_cover.read_text(encoding="ascii")
    require(
        "UNION JOINED 116,13,16 TOTAL 145 UNCOVERED 0 "
        "UNCOVERED_FNV cbf29ce484222325 INCIDENCES 244" in cover_text and
        "WITNESS_FNV ebea2511eb7fa46f RESPONSE_FNV 8e13cc00c5bca917"
        in cover_text and
        "VERDICT PASS FIVE_MASK_UPPER_WITNESS" in cover_text,
        "five-mask response certificate changed",
    )

    repaired_union = typed_union | survivors
    repaired_residual = universe - repaired_union
    require(identity(repaired_union) == (
        1661,
        "5bdd2ebf09e9404a",
        "de00493a80ca88eb4ed802be00fce19967f0978508439bd07afc7393708a4b62",
    ), "THM-4305 union changed")
    require(identity(repaired_residual) == (
        20986,
        "606bf18913a49a14",
        "67561f7f0c5c3a32155811e9978b42b2393c10ea8964387eb712cea6a6683f50",
    ), "THM-4305 residual changed")
    require(repaired_union.isdisjoint(repaired_residual) and
            repaired_union | repaired_residual == universe,
            "THM-4305 partition failed")
    repaired_maximum = max(r for _, r in repaired_residual)
    repaired_top = {
        pair for pair in repaired_residual if pair[1] == repaired_maximum
    }
    require(repaired_maximum == 594 and identity(repaired_top) == (
        25,
        "cce015c81f7121d9",
        "920638d6fb23a8f6492d34cf50e7dc247c2eddfe7ba3f2088c59155e1a56167e",
    ), "THM-4305 residual boundary changed")

    for path, rows in (
        (args.output_union, repaired_union),
        (args.output_residual, repaired_residual),
        (args.output_top, repaired_top),
    ):
        path.write_bytes(payload(rows))
        require(path.read_bytes() == payload(rows), f"write failed: {path}")

    print("THM4305_TYPED_UNION_CONSUMER_V1")
    describe("UNIVERSE", universe)
    describe("THM4302_TYPED_UNION", thm4302_union)
    describe("THM4302_RESIDUAL", thm4302_residual)
    describe("AUDITED_TOP_LAYER", top28)
    describe("CLOSED_ROWS", closed)
    describe("SURVIVING_ROWS", survivors)
    print("SURVIVING_TOP_ROWS " +
          " ".join(f"({q},{r})" for q, r in sorted(survivors)))
    describe("THM4303_TYPED_UNION", typed_union)
    describe("THM4303_RESIDUAL", residual)
    describe("REPAIRED_ROWS", survivors)
    describe("TYPED_UNION", repaired_union)
    describe("FINAL_RESIDUAL", repaired_residual)
    describe("FINAL_RESIDUAL_TOP", repaired_top)
    print("FINAL_RESIDUAL_TOP_ROWS " +
          " ".join(f"({q},{r})" for q, r in sorted(repaired_top)))
    print("TYPE_GUARD row_consequences_only_fixed_carrier_not_physical_entry")
    print("SCOPE FINITE_EXACT_TYPED_UNION_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
