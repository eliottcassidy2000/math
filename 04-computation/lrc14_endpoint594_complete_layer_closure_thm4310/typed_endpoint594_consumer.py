#!/usr/bin/env python3
"""Typed consequence consumer for THM-4310 endpoint-594 closure.

This consumer treats only proved row closure as a typed consequence;
it never identifies the THM-4309 carrier with the THM-4306 rebuilt deck.
"""

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
    parsed = [tuple(map(int, line.split(b","))) for line in raw.splitlines()]
    require(parsed == sorted(set(parsed)), f"unordered/duplicate rows: {path}")
    rows = set(parsed)
    require(raw == payload(rows), f"noncanonical pair bytes: {path}")
    require(all(0 < q < r for q, r in rows), f"invalid pair: {path}")
    return rows


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
    print(f"{label} rows={count} fnv={fingerprint} sha256={digest} "
          f"max_endpoint={maximum}")


def write_rows(path: Path, rows: set[Pair]) -> None:
    path.write_bytes(payload(rows))


def audit_pair_csv(path: Path, target: set[Pair]) -> None:
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == [
            "q", "r", "active", "active_fnv", "active_joint",
            "active_nonjoint", "exposed", "exposed_fnv", "minimum_hits",
            "maximum_hits", "failures", "failure_fnv",
        ], "pair audit header changed")
        rows = list(reader)
    observed = {(int(row["q"]), int(row["r"])) for row in rows}
    require(len(rows) == 25 and observed == target, "pair target changed")
    require(all(int(row["failures"]) == 0 for row in rows),
            "endpoint-594 carrier failure present")
    require(sum(int(row["exposed"]) for row in rows) == 46_178,
            "exposed total changed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--prior-union", type=Path, required=True)
    parser.add_argument("--prior-residual", type=Path, required=True)
    parser.add_argument("--top594", type=Path, required=True)
    parser.add_argument("--carrier-target", type=Path, required=True)
    parser.add_argument("--pair-audit", type=Path, required=True)
    parser.add_argument("--failures", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    universe = read_pairs(args.universe)
    prior_union = read_pairs(args.prior_union)
    prior_residual = read_pairs(args.prior_residual)
    top594 = read_pairs(args.top594)
    carrier_target = read_pairs(args.carrier_target)
    require(identity(universe) == (
        22647, "df5374d4aca67677",
        "14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317",
    ), "fixed universe changed")
    require(identity(prior_union) == (
        2014, "ee275cd5d460d153",
        "836b7ccd0f93268ae039d749c6d778778dd0fcf706dd58cb6dfbd3c340fcbcd1",
    ), "THM-4306 union changed")
    require(identity(prior_residual) == (
        20633, "3acd694e62cb7841",
        "77e54cc8614f2750a7e3a46fd9ec5a3a9b4f86db38121d62802d7387934e7e7f",
    ), "THM-4306 residual changed")
    require(identity(top594) == (
        22, "8413f0d2282e4cd6",
        "2a46ac360974ee95b5c468f1f76fb9ddd6b5165fa6e410dd3b6bad02ca93dd54",
    ), "THM-4306 top layer changed")
    require(identity(carrier_target) == (
        25, "cce015c81f7121d9",
        "920638d6fb23a8f6492d34cf50e7dc247c2eddfe7ba3f2088c59155e1a56167e",
    ), "THM-4309 top layer changed")
    require(prior_union.isdisjoint(prior_residual) and
            prior_union | prior_residual == universe,
            "inherited typed partition changed")
    require(top594 <= prior_residual and
            top594 == {row for row in prior_residual if row[1] == 594},
            "top-594 layer not complete")
    require(top594 <= carrier_target and carrier_target - top594 == {
        (173, 594), (381, 594), (383, 594),
    }, "THM-4306/4307 endpoint-594 seam changed")

    audit_pair_csv(args.pair_audit, carrier_target)
    require(args.failures.read_bytes().splitlines() == [b"q,r,body_hex"],
            "failure ledger is nonempty/noncanonical")

    enlarged = prior_union | top594
    residual = prior_residual - top594
    require(enlarged.isdisjoint(residual) and enlarged | residual == universe,
            "enlarged typed partition failed")
    maximum = max(r for _, r in residual)
    next_top = {row for row in residual if row[1] == maximum}
    require(maximum == 593, "next endpoint is not 593")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.output_dir / "typed_union2036.csv", enlarged)
    write_rows(args.output_dir / "final_residual20611.csv", residual)
    write_rows(args.output_dir / "residual_top593.csv", next_top)

    print("LRC14_ENDPOINT594_C3925_TYPED_CONSUMER_V1")
    describe("CONTROL_THM4306_UNION", prior_union)
    describe("CONTROL_THM4306_RESIDUAL", prior_residual)
    describe("CONTROL_TOP594", top594)
    describe("CONTROL_CARRIER_TARGET594", carrier_target)
    describe("ENLARGED_UNION", enlarged)
    describe("ENLARGED_RESIDUAL", residual)
    describe("NEXT_TOP593", next_top)
    print("NEXT_TOP_ROWS " +
          " ".join(f"({q},{r})" for q, r in sorted(next_top)))
    print("TYPE_GUARD row_consequences_only_separate_decks_no_physical_entry")
    print("SCOPE FINITE_EXACT_TYPED_CONSEQUENCE_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
