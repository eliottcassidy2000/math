#!/usr/bin/env python3
"""Canonical typed consumer for T_4303/T_4305 joined with H_265.

All evidence and output locations are CLI arguments.  The consumer derives
the typed sets from frozen predicates and audits the supplied H_265 response
quotient; it never resolves a scratch directory implicitly.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import deque
from pathlib import Path

Pair = tuple[int, int]
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
F3 = {(96, 595), (100, 595), (210, 595)}
BODIES = (
    0x09392104, 0x0D30E080, 0x0D382104, 0x15386080,
    0x186C9080, 0x19786000, 0x1D489080, 0x1F087000,
)
WITNESSES = (0x22020E09, 0x00868489)


def require(ok: bool, message: str) -> None:
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


def identity(rows: set[Pair]) -> tuple[int, str, str]:
    raw = payload(rows)
    return len(rows), fnv(rows), hashlib.sha256(raw).hexdigest()


def describe(label: str, rows: set[Pair]) -> None:
    count, fingerprint, digest = identity(rows)
    maximum = max((r for _, r in rows), default=0)
    print(f"{label} rows={count} fnv={fingerprint} sha256={digest} "
          f"max_endpoint={maximum}")


def response(mask: int) -> int:
    answer = 0
    for index, body in enumerate(BODIES):
        if mask & body == 0:
            answer |= 1 << index
    return answer


def read_and_audit_quotient(path: Path) -> dict[int, tuple[int, int, bool]]:
    require(hashlib.sha256(path.read_bytes()).hexdigest() ==
            "586e3b735298be94c7c58cae44de8a4fb1bb6aef2ec99b427c3d0278068001ac",
            "H_265 quotient byte identity changed")
    classes: dict[int, tuple[int, int, bool]] = {}
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames ==
                ["response_hex", "count", "least_mask", "maximal"],
                "quotient header changed")
        for row in reader:
            pattern = int(row["response_hex"], 16)
            require(pattern not in classes, "duplicate response class")
            classes[pattern] = (
                int(row["count"]), int(row["least_mask"], 16),
                bool(int(row["maximal"])),
            )
    require(len(classes) == 38, "response-class count changed")
    require(sum(count for count, _, _ in classes.values()) == 1_494_889,
            "common-active count changed")
    require(sum(count for pattern, (count, _, _) in classes.items() if pattern)
            == 57_752, "nonempty response count changed")
    require(0xFF not in classes, "unexpected full responder")
    for pattern, (_, least, _) in classes.items():
        require(response(least) == pattern,
                f"least representative response changed at {pattern:02x}")
    maximal = {
        pattern for pattern in classes if pattern and not any(
            other != pattern and pattern & ~other == 0 for other in classes
        )
    }
    require(maximal == {0x7F, 0xA5, 0xD0}, "maximal antichain changed")
    require({pattern for pattern, (_, _, flag) in classes.items() if flag}
            == maximal, "maximal flags changed")

    distance = [999] * 256
    distance[0] = 0
    queue: deque[int] = deque([0])
    while queue:
        state = queue.popleft()
        for pattern in maximal:
            target = state | pattern
            if distance[target] <= distance[state] + 1:
                continue
            distance[target] = distance[state] + 1
            queue.append(target)
    require(distance[0xFF] == 2, "response DP minimum changed")
    require(all((pattern & 0x82).bit_count() <= 1 for pattern in classes),
            "lower packing support 0x82 changed")
    require(tuple(response(mask) for mask in WITNESSES) == (0x7F, 0xA5),
            "explicit response cover changed")
    return classes


def write_set(output_dir: Path, name: str, rows: set[Pair]) -> None:
    (output_dir / name).write_bytes(payload(rows))
    describe(name.removesuffix(".csv").upper(), rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--universe", required=True, type=Path)
    parser.add_argument("--thm4300-union", required=True, type=Path)
    parser.add_argument("--thm4300-residual", required=True, type=Path)
    parser.add_argument("--ideal", required=True, type=Path)
    parser.add_argument("--quotient", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    universe = read_pairs(args.universe)
    old_union = read_pairs(args.thm4300_union)
    old_residual = read_pairs(args.thm4300_residual)
    ideal = read_pairs(args.ideal)
    read_and_audit_quotient(args.quotient)

    require(identity(universe) == (
        22647, "df5374d4aca67677",
        "14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317",
    ), "fixed-pool universe changed")
    require(identity(old_union) == (
        1624, "11414a33ab91fef6",
        "ef9102553cd030f67ab1bdb7d6965c3efaf4b0d8aa85daa1092354c9703caf26",
    ), "THM-4300 union changed")
    require(identity(old_residual) == (
        21023, "e93e8089e9dc58c0",
        "ce215cb53a742e5e0d0d4f16e344687da9583ce5dc1417c7ea70399fb0bf70ba",
    ) and universe - old_union == old_residual,
            "THM-4300 partition changed")
    require(identity(ideal) == (
        367, "d422b161d94ebae4",
        "750a31e2f0ebe6573835cf9d2cd43e83403d94fe802616e539ee9d330a3dab65",
    ), "H_265 ideal changed")
    require(ideal <= universe, "H_265 escapes fixed-pool universe")

    k596 = {pair for pair in universe if pair[1] >= 596}
    t4302 = old_union | k596
    r4302 = universe - t4302
    layer595 = {pair for pair in r4302 if pair[1] == 595}
    require(identity(t4302) == (
        1633, "b1c8ecf1dd4a71c5",
        "28084fded429f407188e471183d5645d28eded621967e10e386261bbe52844c0",
    ), "THM-4302 union changed")
    require(identity(layer595) == (
        28, "47981ce64825ef2a",
        "c607dab04e4f6849a2226f518771e43b1301d4fc582b47bfa5752c4643c93702",
    ) and F3 <= layer595, "endpoint-595 layer changed")

    # Preserve the prior theorem state as an explicit control.
    t4303 = t4302 | (layer595 - F3)
    r4303 = universe - t4303
    require(identity(t4303) == (
        1658, "43317f1aee06e8bd",
        "bfeb739dcad61dd42bdd9a8b295b6058f3ecee5cc30acd64c469e3b8132393c7",
    ), "T_4303 control union changed")
    require(identity(r4303) == (
        20989, "b0fbaa28440a118f",
        "bbf2dbba58a5492f6e5d136f14940c3bac8b3ddea604b19e3e2b926abb8bad00",
    ), "T_4303 control residual changed")
    require({pair for pair in r4303 if pair[1] == 595} == F3,
            "T_4303 control top layer changed")

    # The caller-specified continuation is T_4305=T_4303 union F3.
    t4305 = t4303 | F3
    r4305 = universe - t4305
    require(identity(t4305) == (
        1661, "5bdd2ebf09e9404a",
        "de00493a80ca88eb4ed802be00fce19967f0978508439bd07afc7393708a4b62",
    ), "T_4305 union changed")
    require(identity(r4305) == (
        20986, "606bf18913a49a14",
        "67561f7f0c5c3a32155811e9978b42b2393c10ea8964387eb712cea6a6683f50",
    ), "T_4305 residual changed")
    top4305 = {pair for pair in r4305
               if pair[1] == max(r for _, r in r4305)}
    require(identity(top4305) == (
        25, "cce015c81f7121d9",
        "920638d6fb23a8f6492d34cf50e7dc247c2eddfe7ba3f2088c59155e1a56167e",
    ), "T_4305 top layer changed")

    overlap4303 = ideal & t4303
    overlap4305 = ideal & t4305
    require(F3.isdisjoint(ideal) and overlap4303 == overlap4305,
            "new endpoint-595 consequences change H_265 overlap")
    require(identity(overlap4305) == (
        14, "ee5919716c115125",
        "355e379f9c43a0a4c51a0c71c58f20b090fa40e2a4368b2c18ce88403a533dd2",
    ), "H_265 typed overlap changed")
    new265 = ideal - t4305
    require(identity(new265) == (
        353, "5aa9e290e62d23d0",
        "63eb9a881c254aa294d3836f9fb938d44aaa727e0f1081bb787cfaf3b08fa981",
    ), "H_265 honest new rows changed")

    combined = t4305 | ideal
    residual = universe - combined
    require(combined.isdisjoint(residual) and combined | residual == universe,
            "combined partition failed")
    maximum = max(r for _, r in residual)
    top = {pair for pair in residual if pair[1] == maximum}
    removed_top = ideal & top4305
    require(removed_top == {(173, 594), (381, 594), (383, 594)},
            "H_265 intersection with T_4305 top layer changed")
    require(identity(combined) == (
        2014, "ee275cd5d460d153",
        "836b7ccd0f93268ae039d749c6d778778dd0fcf706dd58cb6dfbd3c340fcbcd1",
    ), "combined typed union changed")
    require(identity(residual) == (
        20633, "3acd694e62cb7841",
        "77e54cc8614f2750a7e3a46fd9ec5a3a9b4f86db38121d62802d7387934e7e7f",
    ), "combined residual changed")
    require(maximum == 594 and identity(top) == (
        22, "8413f0d2282e4cd6",
        "2a46ac360974ee95b5c468f1f76fb9ddd6b5165fa6e410dd3b6bad02ca93dd54",
    ), "combined maximum layer changed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print("LRC14_T4305_INDEX265_TYPED_CONSUMER_CANDIDATE_V1")
    write_set(args.output_dir, "control_t4303_union1658.csv", t4303)
    write_set(args.output_dir, "control_t4303_residual20989.csv", r4303)
    write_set(args.output_dir, "control_t4303_top595_3.csv", F3)
    write_set(args.output_dir, "t4305_union1661.csv", t4305)
    write_set(args.output_dir, "t4305_residual20986.csv", r4305)
    write_set(args.output_dir, "t4305_top594_25.csv", top4305)
    write_set(args.output_dir, "h265_overlap14.csv", overlap4305)
    write_set(args.output_dir, "h265_new353.csv", new265)
    write_set(args.output_dir, "h265_removed_from_top594_3.csv", removed_top)
    write_set(args.output_dir, "combined_union2014.csv", combined)
    write_set(args.output_dir, "combined_residual20633.csv", residual)
    write_set(args.output_dir, "combined_top594_22.csv", top)
    print("T4305_ADDITIONS " + " ".join(f"({q},{r})" for q, r in sorted(F3)))
    print("T4305_ADDITIONS_INTERSECT_H265 0")
    print("H265_OVERLAP_ROWS " +
          " ".join(f"({q},{r})" for q, r in sorted(overlap4305)))
    print(f"COMBINED_MAX_ENDPOINT {maximum}")
    print("COMBINED_TOP_ROWS " +
          " ".join(f"({q},{r})" for q, r in sorted(top)))
    print("QUOTIENT exact_minimum=2 lower_support=82 "
          "cover=22020e09:7f,00868489:a5")
    print("TYPE_GUARD row_consequences_only_separate_decks_no_physical_entry")
    print("SCOPE FINITE_EXACT_TYPED_UNION_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
