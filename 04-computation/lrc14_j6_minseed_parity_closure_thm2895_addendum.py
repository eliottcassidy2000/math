#!/usr/bin/env python3
"""Locked THM-2895 addendum for the 15 minimum scalar-bootstrap seeds.

The committed minimum-seed activation battery proves that scalar bootstrap
closes every nonseed apex once the listed ordered seeds are certified on
their *actual* prefixes.  This probe applies the locked suffix-H4 verifier
to exactly those 15 branch obligations.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARITY_PATH = (
    ROOT
    / "04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895.py"
)
PARITY_SHA = "970d77503f8d56d737e223dabb3c3562d7b19cd018ca75398e3deb054715e5f6"
PARITY_OUTPUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_suffix_parity_flag_closure_thm2895.out"
)
PARITY_OUTPUT_SHA = "c11260f6544a319e1cc1862c9221b188a4314860422470e465b82e7ce492b1b4"
MINSEED_PATH = (
    ROOT
    / "04-computation/lrc14_j6_suffix_bootstrap_minseed_battery_codex_20260729.py"
)
MINSEED_SHA = "a7a77dc433b21d94a54524064ccd62e553ed67ae8a3d3364bf79c41c36849d04"
MINSEED_OUTPUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_suffix_bootstrap_minseed_battery_codex_20260729.out"
)
MINSEED_OUTPUT_SHA = "bdcf896d152d206b3ae77a3568609887190e1ba991909481769d7ad560a68835"

SEEDS = (
    ((2, 8, 9, 10, 11, 13, 14), (17, 23, 46, 24)),
    ((1, 3, 9, 10, 11, 12, 14), (39, 23, 16)),
    ((2, 5, 9, 11, 12, 13, 14), (16, 23, 19, 40, 46)),
    ((2, 3, 4, 5, 6, 7, 8), (22, 26, 33)),
)
EXPECTED_AGGREGATES = (15, 15, 0, 1464, 1461, 1444, 1462, 2, 17_971)
EXPECTED_H_SIZES = (30, 17, 24, 15, 12, 7, 7, 16, 12, 12, 10, 7, 11, 5, 7)
EXPECTED_RECURSIVE = (
    (
        (2, 5, 9, 11, 12, 13, 14),
        1,
        16,
        (23, 25),
        (20, 37),
        ((20, 37),),
        ("33/7252",),
        (31,),
        13,
    ),
    (
        (2, 5, 9, 11, 12, 13, 14),
        1,
        16,
        (25, 34),
        (20, 23),
        ((20, 23),),
        ("19/4508",),
        (33,),
        15,
    ),
)
EXPECTED_LEDGER_DIGEST: str | None = (
    "08712e5aff4875c8fd2ba336f21b29d4ac2694ddf69b1aaa258f6df74073a05c"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path, expected: str, name: str):
    require(file_sha256(path) == expected, f"{path.name} changed")
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(
    file_sha256(PARITY_OUTPUT) == PARITY_OUTPUT_SHA,
    "locked parity-descent output changed",
)
require(
    file_sha256(MINSEED_OUTPUT) == MINSEED_OUTPUT_SHA,
    "minimum-seed output changed",
)
Q = load(PARITY_PATH, PARITY_SHA, "j6_minseed_parity")
M = load(MINSEED_PATH, MINSEED_SHA, "j6_minseed_source")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def seed_branches() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    locked_seeds = tuple(
        row[5] for row in M.EXPECTED_ROOT_RESULTS
    )
    require(
        locked_seeds == tuple(seed for _, seed in SEEDS),
        "deterministic minimum seeds changed",
    )
    for body, seed in SEEDS:
        root = Q.S.global_root(body)
        gate = tuple(speed for _, speed in root["top"][: root["K"]])
        require(
            all(apex in gate for apex in seed),
            f"seed left root gate: {body}",
        )
        for seed_index, apex in enumerate(seed, start=1):
            prefix = seed[:seed_index]
            excluded = set(prefix)
            carrier = Q.H.R.subtract_local(root["good"], apex)
            carrier_mass = mass(carrier)
            direct, direct_r, direct_m = Q.H.T.CORE.good_norm(
                tuple(sorted((*body, apex)))
            )
            require(
                carrier == direct
                and len(carrier) == direct_r
                and carrier_mass == direct_m > 0,
                f"seed carrier reconstruction failed: {body}, {apex}",
            )
            ranked, tail_single = Q.H.globally_ranked(carrier, excluded)
            q1 = ranked[0][0]
            singleton_margin = F(3, 7) * carrier_mass - q1
            if singleton_margin <= 0:
                level = None
                tail_first = None
                core: tuple[int, ...] = ()
            else:
                level = (carrier_mass - q1) / 4
                threshold = (
                    Q.H.S2
                    * direct_r
                    / (7 * (level - carrier_mass / 7))
                )
                tail_first = max(Q.H.FIRST_EXTERNAL, ceiling(threshold))
                by_speed = {speed: value for value, speed in ranked}
                if tail_first > Q.H.HORIZON + 1:
                    by_speed.update(
                        {
                            speed: value
                            for value, speed in Q.H.T.coverages_many(
                                carrier,
                                [
                                    speed
                                    for speed in range(
                                        Q.H.HORIZON + 1,
                                        tail_first,
                                    )
                                    if speed not in excluded
                                ],
                            )
                        }
                    )
                require(
                    carrier_mass / 7
                    + Q.H.S2 * direct_r / (7 * tail_first)
                    <= level,
                    "seed H4 tail did not seal",
                )
                core = tuple(
                    sorted(
                        speed
                        for speed, value in by_speed.items()
                        if speed not in excluded and value >= level
                    )
                )
            root_rank = gate.index(apex) + 1
            rows.append(
                {
                    "root": root,
                    "branch": {
                        "body": body,
                        "rank": seed_index,
                        "root_rank": root_rank,
                        "apex": apex,
                        "excluded_prefix": prefix,
                        "m": carrier_mass,
                        "r": direct_r,
                    },
                    "carrier": carrier,
                    "excluded": excluded,
                    "q1": q1,
                    "tail_single": tail_single,
                    "singleton_margin": singleton_margin,
                    "level": level,
                    "Htail": tail_first,
                    "H": core,
                }
            )
    require(len(rows) == 15, "minimum-seed branch universe changed")
    return rows


def main() -> None:
    branches = seed_branches()
    unavailable = [
        row for row in branches if row["singleton_margin"] <= 0
    ]
    print("J6 MINIMUM-SEED PARITY CLOSURE VERIFIER")
    print(
        f"seed_branches={len(branches)};"
        f"B1_below_3h7={len(branches)-len(unavailable)};"
        f"B1_failures={len(unavailable)}"
    )
    for row in branches:
        branch = row["branch"]
        print(
            f"BRANCH E={branch['body']};seed_index={branch['rank']};"
            f"root_rank={branch['root_rank']};apex={branch['apex']};"
            f"prefix={branch['excluded_prefix']};"
            f"3h7-q1={ftext(row['singleton_margin'])};"
            f"H={len(row['H'])};Htail={row['Htail']}"
        )
    if unavailable:
        print("scope=singleton-complement obstruction found")
        return

    all_pairs: list[tuple[dict[str, object], dict[str, object]]] = []
    for row in branches:
        local = [
            Q.pair_residual(row, hpair)
            for hpair in combinations(row["H"], 2)
        ]
        all_pairs.extend((row, pair) for pair in local)
        branch = row["branch"]
        print(
            f"PAIR_SUMMARY E={branch['body']};seed_index={branch['rank']};"
            f"apex={branch['apex']};pairs={len(local)};"
            f"top3={sum(pair['direct_margin'] > 0 for pair in local)};"
            f"B2q1={sum(pair['pair_margin'] > 0 for pair in local)};"
            f"union={sum(pair['adaptive_closed'] for pair in local)}"
        )

    hard = [
        (row, pair)
        for row, pair in all_pairs
        if not pair["adaptive_closed"]
    ]
    aggregates = (
        len(branches),
        len(branches) - len(unavailable),
        len(unavailable),
        len(all_pairs),
        sum(pair["direct_margin"] > 0 for _, pair in all_pairs),
        sum(pair["pair_margin"] > 0 for _, pair in all_pairs),
        len(all_pairs) - len(hard),
        len(hard),
        sum(pair["cap"]["paid"] for _, pair in all_pairs),
    )
    require(
        aggregates == EXPECTED_AGGREGATES,
        f"minimum-seed aggregate census changed: {aggregates}",
    )
    require(
        tuple(len(row["H"]) for row in branches) == EXPECTED_H_SIZES,
        "minimum-seed H4 sizes changed",
    )
    recursive_rows = [
        (row, pair, Q.recursive_k3_close(row, pair))
        for row, pair in hard
    ]
    recursive_data = tuple(
        (
            row["branch"]["body"],
            row["branch"]["rank"],
            row["branch"]["apex"],
            pair["hpair"],
            recursive["H"],
            recursive["heavy"],
            tuple(ftext(value) for value in recursive["longest"]),
            recursive["horizons"],
            recursive["checks"],
        )
        for row, pair, recursive in recursive_rows
    )
    require(
        recursive_data == EXPECTED_RECURSIVE,
        f"minimum-seed recursive rows changed: {recursive_data}",
    )
    for row, pair, recursive in recursive_rows:
        branch = row["branch"]
        print(
            f"RECURSIVE E={branch['body']};seed_index={branch['rank']};"
            f"apex={branch['apex']};H4pair={pair['hpair']};"
            f"H2={recursive['H']};heavy={recursive['heavy']};"
            f"ell={tuple(ftext(value) for value in recursive['longest'])};"
            f"W={recursive['horizons']};checks={recursive['checks']};"
            f"covers={recursive['covers']}"
        )
    print(
        f"totals=pairs:{len(all_pairs)},"
        f"top3:{sum(pair['direct_margin'] > 0 for _, pair in all_pairs)},"
        f"B2q1:{sum(pair['pair_margin'] > 0 for _, pair in all_pairs)},"
        f"union:{len(all_pairs)-len(hard)},hard:{len(hard)},"
        f"paid_pairs:{sum(pair['cap']['paid'] for _, pair in all_pairs)},"
        f"recursive_closed:{sum(recursive['closed'] for _, _, recursive in recursive_rows)},"
        f"recursive_open:{sum(not recursive['closed'] for _, _, recursive in recursive_rows)},"
        f"heavy_edges:{sum(len(recursive['heavy']) for _, _, recursive in recursive_rows)},"
        f"singleton_checks:{sum(recursive['checks'] for _, _, recursive in recursive_rows)}"
    )
    require(
        all(recursive["closed"] for _, _, recursive in recursive_rows),
        "minimum-seed recursive residual remains open",
    )
    ledger = hashlib.sha256()
    ledger.update(b"LRC14/j6/minimum-seed-H4-parity-descent/v1\n")
    for row in branches:
        branch = row["branch"]
        ledger.update(
            (
                f"seed_index={branch['rank']};"
                f"root_rank={branch['root_rank']};"
            ).encode()
        )
        ledger.update(Q.branch_ledger_line(row).encode())
    for row, pair in all_pairs:
        ledger.update(Q.pair_ledger_line(row, pair).encode())
    for row, pair, recursive in recursive_rows:
        ledger.update(Q.recursive_ledger_line(row, pair, recursive).encode())
    ledger_digest = ledger.hexdigest()
    if EXPECTED_LEDGER_DIGEST is not None:
        require(
            ledger_digest == EXPECTED_LEDGER_DIGEST,
            f"minimum-seed ledger digest changed: {ledger_digest}",
        )
    print(f"canonical_ledger_sha256={ledger_digest}")
    print(
        "mode="
        + ("DISCOVERY" if EXPECTED_LEDGER_DIGEST is None else "LOCKED")
    )
    print(
        "consequence=all 15 ordered seed branches close;"
        "locked scalar bootstrap closes all four root gates"
    )
    print("scope=four-root minimum-seed atlas;not uniform;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
