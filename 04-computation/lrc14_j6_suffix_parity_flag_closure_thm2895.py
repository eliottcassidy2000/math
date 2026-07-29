#!/usr/bin/env python3
"""Locked suffix-aware parity-flag verifier for the four-root j=6 battery.

The committed ranked-suffix scalar battery has 25 open branches on four
seven-body roots.  This script preserves each branch's full excluded root
prefix, computes its singleton complement cap and H4 core, then tests every
H4-pair residual against two incomparable triple-cover certificates:

  (1) the exact global top-three singleton sum;
  (2) the exact globally sealed B2 + q1 cap.

Only residuals failing both are sent to exact triple maximization.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUFFIX_PATH = (
    ROOT
    / "04-computation/lrc14_j6_rank_first_suffix_scalar_battery_codex_20260729.py"
)
SUFFIX_SHA = "6434f020c5aa4000ac81fa081881d93ac0b4190516f854fbd9d8493475baf539"
SUFFIX_OUTPUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_rank_first_suffix_scalar_battery_codex_20260729.out"
)
SUFFIX_OUTPUT_SHA = "d03f4be7ead138135447b4d720e91d215b1668c2182a099da92129289e605ee9"
H4_PATH = (
    ROOT
    / "04-computation/lrc14_j6_h4_pair_residual_exact_kernel_codex_20260729.py"
)
H4_SHA = "b82f318bf89ffd3ab4c918c87736461d068e03f25941aa25a0961d0f74b4d70a"

EXPECTED_AGGREGATES = (25, 25, 0, 784, 771, 773, 779, 5, 7_551)
EXPECTED_RECURSIVE_AGGREGATES = (5, 5, 0, 5, 86, 5, 86)
EXPECTED_RECURSIVE_ROWS = (
    (
        (2, 8, 9, 10, 11, 13, 14),
        1,
        19,
        (37, 108),
        (17, 23, 46),
        ((17, 23),),
        ("121/34776",),
        (41,),
        23,
    ),
    (
        (2, 8, 9, 10, 11, 13, 14),
        1,
        19,
        (37, 125),
        (17, 23, 46),
        ((17, 23),),
        ("13/3500",),
        (38,),
        20,
    ),
    (
        (2, 8, 9, 10, 11, 13, 14),
        1,
        19,
        (108, 125),
        (17, 23, 46),
        ((17, 23),),
        ("125/28728",),
        (32,),
        15,
    ),
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
    "c541e41f688dd9873e57b5316a8c9b28c496e9bd808a65c62f7d20b7a3b87d4e"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path, expected: str | None, name: str):
    if expected is not None:
        require(file_sha256(path) == expected, f"{path.name} changed")
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(
    file_sha256(SUFFIX_OUTPUT) == SUFFIX_OUTPUT_SHA,
    "ranked-suffix transcript changed",
)
S = load(SUFFIX_PATH, SUFFIX_SHA, "j6_suffix_h4_source")
H = load(H4_PATH, H4_SHA, "j6_suffix_h4_helpers")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def sealed_top(
    carrier: list[tuple[F, F]],
    excluded: set[int],
    count: int,
    base_ranked: list[tuple[F, int]],
) -> tuple[tuple[tuple[F, int], ...], int]:
    """Globally seal the first ``count`` singleton ranks."""

    carrier_mass = mass(carrier)
    components = len(carrier)
    q_count = base_ranked[count - 1][0]
    require(q_count > carrier_mass / 7, "requested rank misses limit")
    threshold = H.S2 * components / (7 * (q_count - carrier_mass / 7))
    tail_first = max(H.HORIZON + 1, ceiling(threshold))
    rows = list(base_ranked)
    if tail_first > H.HORIZON + 1:
        rows.extend(
            H.T.coverages_many(
                carrier,
                [
                    speed
                    for speed in range(H.HORIZON + 1, tail_first)
                    if speed not in excluded
                ],
            )
        )
    require(
        carrier_mass / 7
        + H.S2 * components / (7 * tail_first)
        <= q_count,
        "top-rank tail did not seal",
    )
    return (
        tuple(sorted(rows, key=lambda item: (-item[0], item[1]))[:count]),
        tail_first,
    )


def open_branch_profiles() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for body in S.BATTERY_ROOTS:
        root = S.global_root(body)
        for rank in range(1, root["K"] + 1):
            branch = S.suffix_branch(root, rank)
            if branch["closed"]:
                continue
            carrier = H.R.subtract_local(root["good"], branch["apex"])
            require(
                mass(carrier) == branch["m"]
                and len(carrier) == branch["r"],
                f"suffix carrier mismatch: {body}, rank {rank}",
            )
            excluded = set(branch["excluded_prefix"])
            ranked, tail_single = H.globally_ranked(carrier, excluded)
            q1 = ranked[0][0]
            singleton_margin = F(3, 7) * branch["m"] - q1
            if singleton_margin <= 0:
                core: tuple[int, ...] = ()
                level = None
                tail_first = None
            else:
                level = (branch["m"] - q1) / 4
                require(level > branch["m"] / 7, "H4 level is not finite")
                threshold = (
                    H.S2
                    * branch["r"]
                    / (7 * (level - branch["m"] / 7))
                )
                tail_first = max(H.FIRST_EXTERNAL, ceiling(threshold))
                by_speed = {speed: value for value, speed in ranked}
                if tail_first > H.HORIZON + 1:
                    by_speed.update(
                        {
                            speed: value
                            for value, speed in H.T.coverages_many(
                                carrier,
                                [
                                    speed
                                    for speed in range(
                                        H.HORIZON + 1,
                                        tail_first,
                                    )
                                    if speed not in excluded
                                ],
                            )
                        }
                    )
                require(
                    branch["m"] / 7
                    + H.S2 * branch["r"] / (7 * tail_first)
                    <= level,
                    "H4 tail did not seal",
                )
                core = tuple(
                    sorted(
                        speed
                        for speed, value in by_speed.items()
                        if speed not in excluded and value >= level
                    )
                )
            rows.append(
                {
                    "root": root,
                    "branch": branch,
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
    require(len(rows) == 25, "open suffix-branch universe changed")
    return rows


def pair_residual(
    branch: dict[str, object],
    hpair: tuple[int, int],
) -> dict[str, object]:
    residual = H.R.subtract_local_multi(branch["carrier"], hpair)
    residual_mass = mass(residual)
    require(residual_mass > 0, "H4 pair already covers suffix carrier")
    family = tuple(
        sorted(
            (
                *branch["branch"]["body"],
                branch["branch"]["apex"],
                *hpair,
            )
        )
    )
    direct, direct_r, direct_m = H.T.CORE.good_norm(family)
    require(
        residual == direct
        and len(residual) == direct_r
        and residual_mass == direct_m,
        f"literal/direct residual mismatch: {family}",
    )
    excluded = set(branch["excluded"]) | set(hpair)
    cap = H.pair_cap(residual, excluded)
    top3, top3_tail = sealed_top(
        residual,
        excluded,
        3,
        cap["ranked"],
    )
    direct_margin = residual_mass - sum(
        (value for value, _ in top3),
        F(0),
    )
    pair_margin = residual_mass - cap["cap"] - cap["q1"]
    return {
        "hpair": hpair,
        "residual": residual,
        "m": residual_mass,
        "r": len(residual),
        "cap": cap,
        "top3": top3,
        "top3_tail": top3_tail,
        "direct_margin": direct_margin,
        "pair_margin": pair_margin,
        "adaptive_closed": direct_margin > 0 or pair_margin > 0,
    }


def recursive_k3_close(
    branch: dict[str, object],
    pair: dict[str, object],
) -> dict[str, object]:
    """Apply THM-2893 with (k,s,ell)=(3,2,2) to one hard residual."""

    residual = pair["residual"]
    residual_mass = pair["m"]
    q1 = pair["cap"]["q1"]
    singleton_margin = F(5, 7) * residual_mass - q1
    require(singleton_margin > 0, "recursive singleton cap reaches 5h/7")
    theta = residual_mass - q1
    level = theta / 2
    require(level > residual_mass / 7, "recursive H2 is not finite")
    threshold = (
        H.S2 * pair["r"] / (7 * (level - residual_mass / 7))
    )
    tail_first = max(H.FIRST_EXTERNAL, ceiling(threshold))
    excluded = set(branch["excluded"]) | set(pair["hpair"])
    by_speed = {
        speed: value
        for value, speed in pair["cap"]["ranked"]
    }
    if tail_first > H.HORIZON + 1:
        by_speed.update(
            {
                speed: value
                for value, speed in H.T.coverages_many(
                    residual,
                    [
                        speed
                        for speed in range(H.HORIZON + 1, tail_first)
                        if speed not in excluded
                    ],
                )
            }
        )
    require(
        residual_mass / 7
        + H.S2 * pair["r"] / (7 * tail_first)
        <= level,
        "recursive H2 tail did not seal",
    )
    core = tuple(
        sorted(
            speed
            for speed, value in by_speed.items()
            if speed not in excluded and value >= level
        )
    )
    heavy: list[tuple[int, int]] = []
    heavy_residuals: list[list[tuple[F, F]]] = []
    heavy_hash = hashlib.sha256(b"LRC14/j6/suffix-H2-heavy/v1\n")
    for edge in combinations(core, 2):
        after = H.R.subtract_local_multi(residual, edge)
        union = residual_mass - mass(after)
        heavy_hash.update(
            (
                f"P={edge[0]},{edge[1]};U={ftext(union)};"
                f"L={ftext(mass(after))};r={len(after)}\n"
            ).encode()
        )
        if union >= theta:
            heavy.append(edge)
            heavy_residuals.append(after)

    checks = 0
    covers: list[tuple[tuple[int, int], int]] = []
    horizons: list[int] = []
    longest_components: list[F] = []
    direct_controls = 0
    noncontainment_controls = 0
    singleton_hash = hashlib.sha256(
        b"LRC14/j6/suffix-H2-singleton-noncontainment/v1\n"
    )
    for edge, after in zip(heavy, heavy_residuals):
        after_mass = mass(after)
        require(after_mass > 0, "recursive heavy edge covers residual")
        family = tuple(
            sorted(
                (
                    *branch["branch"]["body"],
                    branch["branch"]["apex"],
                    *pair["hpair"],
                    *edge,
                )
            )
        )
        direct, direct_r, direct_m = H.T.CORE.good_norm(family)
        require(
            after == direct
            and len(after) == direct_r
            and after_mass == direct_m,
            f"recursive literal/direct reconstruction failed: {family}",
        )
        direct_controls += 1

        # A component of D_w has length at most 1/(7w).  Hence the longest
        # connected interval of the residual gives the geometric singleton
        # horizon w <= floor((1/7)/ell), independent of discrepancy.
        longest = max(right - left for left, right in after)
        longest_components.append(longest)
        horizon_fraction = F(1, 7) / longest
        horizon = horizon_fraction.numerator // horizon_fraction.denominator
        horizons.append(horizon)
        edge_excluded = excluded | set(edge)
        for speed in range(H.FIRST_EXTERNAL, horizon + 1):
            if speed in edge_excluded:
                continue
            value = H.T.coverage(after, speed)
            survivor = H.R.subtract_local(after, speed)
            survivor_mass = mass(survivor)
            require(
                survivor_mass == after_mass - value,
                "singleton coverage/subtraction control failed",
            )
            singleton_hash.update(
                (
                    f"E={edge[0]},{edge[1]};w={speed};"
                    f"c={ftext(value)};L={ftext(survivor_mass)}\n"
                ).encode()
            )
            checks += 1
            noncontainment_controls += 1
            if value == after_mass:
                covers.append((edge, speed))
        require(
            F(1, 7 * (horizon + 1)) < longest,
            "geometric singleton horizon did not close",
        )
    return {
        "singleton_margin": singleton_margin,
        "theta": theta,
        "level": level,
        "Htail": tail_first,
        "H": core,
        "heavy": tuple(heavy),
        "heavy_digest": heavy_hash.hexdigest(),
        "longest": tuple(longest_components),
        "horizons": tuple(horizons),
        "checks": checks,
        "direct_controls": direct_controls,
        "noncontainment_controls": noncontainment_controls,
        "singleton_digest": singleton_hash.hexdigest(),
        "covers": tuple(covers),
        "closed": not covers,
    }


def branch_ledger_line(row: dict[str, object]) -> str:
    branch = row["branch"]
    return (
        f"B:E={','.join(map(str, branch['body']))};rank={branch['rank']};"
        f"a={branch['apex']};"
        f"prefix={','.join(map(str, branch['excluded_prefix']))};"
        f"h={ftext(branch['m'])};r={branch['r']};"
        f"q1={ftext(row['q1'])};"
        f"d1={ftext(row['singleton_margin'])};"
        f"level={ftext(row['level'])};Htail={row['Htail']};"
        f"H={','.join(map(str, row['H']))}\n"
    )


def pair_ledger_line(
    row: dict[str, object],
    pair: dict[str, object],
) -> str:
    branch = row["branch"]
    cap = pair["cap"]
    return (
        f"P:E={','.join(map(str, branch['body']))};rank={branch['rank']};"
        f"a={branch['apex']};L={pair['hpair'][0]},{pair['hpair'][1]};"
        f"h={ftext(pair['m'])};r={pair['r']};"
        f"direct={ftext(pair['direct_margin'])};"
        f"B2q1={ftext(pair['pair_margin'])};"
        f"q1={ftext(cap['q1'])};"
        f"B2head={ftext(cap['head'])};B2tail={ftext(cap['tail'])};"
        f"B2={ftext(cap['cap'])};"
        f"maxP={cap['maximizer'][0]},{cap['maximizer'][1]};"
        f"paid={cap['paid']};paidsha={cap['paid_digest']};"
        f"top3tail={pair['top3_tail']};top3="
        + ",".join(
            f"{speed}:{ftext(value)}"
            for value, speed in pair["top3"]
        )
        + "\n"
    )


def recursive_ledger_line(
    row: dict[str, object],
    pair: dict[str, object],
    recursive: dict[str, object],
) -> str:
    branch = row["branch"]
    return (
        f"R:E={','.join(map(str, branch['body']))};rank={branch['rank']};"
        f"a={branch['apex']};L={pair['hpair'][0]},{pair['hpair'][1]};"
        f"d1={ftext(recursive['singleton_margin'])};"
        f"theta={ftext(recursive['theta'])};"
        f"level={ftext(recursive['level'])};"
        f"Htail={recursive['Htail']};"
        f"H={','.join(map(str, recursive['H']))};"
        + "heavy="
        + ",".join(f"{first}:{second}" for first, second in recursive["heavy"])
        + ";ell="
        + ",".join(ftext(value) for value in recursive["longest"])
        + ";W="
        + ",".join(map(str, recursive["horizons"]))
        + f";checks={recursive['checks']};"
        f"heavysha={recursive['heavy_digest']};"
        f"singlesha={recursive['singleton_digest']};"
        f"covers={recursive['covers']}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h-only", action="store_true")
    parser.add_argument("--skip-exact", action="store_true")
    parser.add_argument("--exact-fallback", action="store_true")
    args = parser.parse_args()

    branches = open_branch_profiles()
    unavailable = [
        row for row in branches if row["singleton_margin"] <= 0
    ]
    print("J6 SUFFIX-AWARE PARITY-FLAG CLOSURE VERIFIER")
    print(
        f"open_suffix_branches={len(branches)};"
        f"B1_below_3h7={len(branches)-len(unavailable)};"
        f"B1_failures={len(unavailable)}"
    )
    for row in branches:
        branch = row["branch"]
        print(
            f"BRANCH E={branch['body']};rank={branch['rank']};"
            f"apex={branch['apex']};prefix={branch['excluded_prefix']};"
            f"3h7-q1={ftext(row['singleton_margin'])};"
            f"H={len(row['H'])};Htail={row['Htail']}"
        )
    if args.h_only or unavailable:
        if unavailable:
            for row in unavailable:
                branch = row["branch"]
                print(
                    f"B1_FAILURE E={branch['body']};rank={branch['rank']};"
                    f"apex={branch['apex']};"
                    f"margin={ftext(row['singleton_margin'])}"
                )
        print("scope=H4 structural layer only")
        print("all_exact_controls=PASS")
        return

    all_pairs: list[tuple[dict[str, object], dict[str, object]]] = []
    for row in branches:
        local = [
            pair_residual(row, hpair)
            for hpair in combinations(row["H"], 2)
        ]
        all_pairs.extend((row, pair) for pair in local)
        branch = row["branch"]
        print(
            f"PAIR_SUMMARY E={branch['body']};rank={branch['rank']};"
            f"apex={branch['apex']};pairs={len(local)};"
            f"top3={sum(pair['direct_margin'] > 0 for pair in local)};"
            f"B2q1={sum(pair['pair_margin'] > 0 for pair in local)};"
            f"union={sum(pair['adaptive_closed'] for pair in local)}"
        )

    adaptive_failures = [
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
        len(all_pairs) - len(adaptive_failures),
        len(adaptive_failures),
        sum(pair["cap"]["paid"] for _, pair in all_pairs),
    )
    require(
        aggregates == EXPECTED_AGGREGATES,
        f"suffix-H4 aggregate census changed: {aggregates}",
    )
    print(
        f"pair_totals={len(all_pairs)};"
        f"top3={sum(pair['direct_margin'] > 0 for _, pair in all_pairs)};"
        f"B2q1={sum(pair['pair_margin'] > 0 for _, pair in all_pairs)};"
        f"adaptive_union={len(all_pairs)-len(adaptive_failures)};"
        f"adaptive_failures={len(adaptive_failures)};"
        f"paid_pairs={sum(pair['cap']['paid'] for _, pair in all_pairs)}"
    )

    recursive_rows = []
    for row, pair in adaptive_failures:
        recursive = recursive_k3_close(row, pair)
        recursive_rows.append((row, pair, recursive))
        branch = row["branch"]
        print(
            f"RECURSIVE E={branch['body']};rank={branch['rank']};"
            f"apex={branch['apex']};H4pair={pair['hpair']};"
            f"H2={recursive['H']};heavy={recursive['heavy']};"
            f"ell={tuple(ftext(value) for value in recursive['longest'])};"
            f"horizons={recursive['horizons']};"
            f"checks={recursive['checks']};"
            f"direct_controls={recursive['direct_controls']};"
            f"noncontainment_controls={recursive['noncontainment_controls']};"
            f"covers={recursive['covers']}"
        )
    print(
        f"recursive_rows={len(recursive_rows)};"
        f"recursive_closed={sum(recursive['closed'] for _, _, recursive in recursive_rows)};"
        f"recursive_open={sum(not recursive['closed'] for _, _, recursive in recursive_rows)};"
        f"recursive_heavy_edges={sum(len(recursive['heavy']) for _, _, recursive in recursive_rows)};"
        f"recursive_singleton_checks={sum(recursive['checks'] for _, _, recursive in recursive_rows)};"
        f"recursive_direct_controls={sum(recursive['direct_controls'] for _, _, recursive in recursive_rows)};"
        f"recursive_noncontainment_controls={sum(recursive['noncontainment_controls'] for _, _, recursive in recursive_rows)}"
    )
    recursive_aggregates = (
        len(recursive_rows),
        sum(recursive["closed"] for _, _, recursive in recursive_rows),
        sum(not recursive["closed"] for _, _, recursive in recursive_rows),
        sum(len(recursive["heavy"]) for _, _, recursive in recursive_rows),
        sum(recursive["checks"] for _, _, recursive in recursive_rows),
        sum(
            recursive["direct_controls"]
            for _, _, recursive in recursive_rows
        ),
        sum(
            recursive["noncontainment_controls"]
            for _, _, recursive in recursive_rows
        ),
    )
    require(
        recursive_aggregates == EXPECTED_RECURSIVE_AGGREGATES,
        f"recursive aggregate census changed: {recursive_aggregates}",
    )
    recursive_row_data = tuple(
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
        recursive_row_data == EXPECTED_RECURSIVE_ROWS,
        f"recursive hostile rows changed: {recursive_row_data}",
    )

    ledger = hashlib.sha256()
    ledger.update(b"LRC14/j6/suffix-H4-parity-descent/v1\n")
    for row in branches:
        ledger.update(branch_ledger_line(row).encode())
    for row, pair in all_pairs:
        ledger.update(pair_ledger_line(row, pair).encode())
    for row, pair, recursive in recursive_rows:
        ledger.update(recursive_ledger_line(row, pair, recursive).encode())
    ledger_digest = ledger.hexdigest()
    if EXPECTED_LEDGER_DIGEST is not None:
        require(
            ledger_digest == EXPECTED_LEDGER_DIGEST,
            f"canonical ledger digest changed: {ledger_digest}",
        )

    exact_rows = []
    recursive_open = [
        (row, pair)
        for row, pair, recursive in recursive_rows
        if not recursive["closed"]
    ]
    if args.exact_fallback and not args.skip_exact:
        for row, pair in recursive_open:
            exact = H.exact_triple_cap(pair["residual"], pair["cap"])
            exact_rows.append((row, pair, exact))
            branch = row["branch"]
            print(
                f"EXACT E={branch['body']};rank={branch['rank']};"
                f"apex={branch['apex']};Hpair={pair['hpair']};"
                f"margin={ftext(exact['margin'])};"
                f"maxtriple={exact['maximizer']};paid={exact['paid']};"
                f"tail_active={int(exact['tail'] >= exact['head'])}"
            )
        print(
            f"exact_triples={len(exact_rows)};"
            f"exact_closed={sum(exact['margin'] > 0 for _, _, exact in exact_rows)};"
            f"exact_open={sum(exact['margin'] <= 0 for _, _, exact in exact_rows)};"
            f"exact_paid={sum(exact['paid'] for _, _, exact in exact_rows)}"
        )
    print(f"canonical_ledger_sha256={ledger_digest}")
    print(
        "mode="
        + (
            "DISCOVERY"
            if EXPECTED_LEDGER_DIGEST is None
            else "LOCKED"
        )
    )
    print("scope=four-root/25-open-branch battery;not uniform")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
