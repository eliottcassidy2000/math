#!/usr/bin/env python3
"""Close the two K>=20 pair-flag survivors through H4 pair children."""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CENSUS_PATH = (
    ROOT / ".scratch/lrc_sparse_j6_gate_audit_20260729/"
    "all_root_q5_pair_suffix_census.py"
)
CENSUS_SHA = "02e8802eb60f6d598903b84e47f454592a78001f5143c02450533c0899e06c07"
SURVIVORS = (
    ((1, 3, 7, 8, 10, 11, 13), 1, 18),
    ((1, 4, 9, 10, 12, 13, 14), 1, 22),
)
EXPECTED_TOTALS: tuple[int, ...] | None = (
    2,
    651,
    412,
    602,
    602,
    49,
    4883,
    165,
    124,
    2325,
    0,
)
EXPECTED_ROOT_ROWS: tuple[tuple[object, ...], ...] | None = (
    (
        (1, 3, 7, 8, 10, 11, 13),
        18,
        F(21237, 140140),
        F(66259, 1331330),
        F(282883, 18638620),
        F(379, 14896),
        1171,
        (
            17, 19, 21, 23, 25, 29, 34, 37, 38, 42, 46, 47, 53,
            68, 75, 80, 84, 100, 104, 125, 130, 143, 156, 173, 216,
        ),
        300,
        225,
        287,
        287,
    ),
    (
        (1, 4, 9, 10, 12, 13, 14),
        22,
        F(4643, 25740),
        F(29257, 504504),
        F(48721, 2522520),
        F(308729, 10090080),
        1422,
        (
            16, 17, 19, 23, 24, 29, 32, 34, 38, 40, 43, 46, 48,
            62, 69, 74, 75, 79, 86, 108, 115, 125, 130, 140, 156,
            175, 182,
        ),
        351,
        187,
        315,
        315,
    ),
)
EXPECTED_MINIMUM: tuple[object, ...] | None = (
    (1, 4, 9, 10, 12, 13, 14),
    22,
    (19, 34),
    "-555503/37934820",
    "-59209/2408560",
    (16, 23),
)
EXPECTED_LEDGER_DIGEST: str | None = (
    "e84b6f0472307458b85150c3f545a57e8b96cc2f15108d6138bb436484c5c988"
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


C = load(CENSUS_PATH, CENSUS_SHA, "j6_two_survivor_h4")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def parent_branch(
    body: tuple[int, ...],
    rank: int,
    expected_apex: int,
) -> dict[str, object]:
    root = C.A.profile_body(body)
    apex = root["top"][rank - 1][1]
    require(apex == expected_apex, f"survivor apex changed: {body}")
    prefix = tuple(speed for _, speed in root["top"][:rank])
    excluded = set(prefix)
    root_good, _, _ = C.A.T.CORE.good_norm(body)
    carrier = C.H.R.subtract_local(root_good, apex)
    direct, direct_r, direct_h = C.H.T.CORE.good_norm(
        tuple(sorted((*body, apex)))
    )
    h = C.mass(carrier)
    require(
        carrier == direct and len(carrier) == direct_r and h == direct_h > 0,
        f"parent literal/direct mismatch: {body}",
    )
    ranked, _, _ = C.sealed_ranking(carrier, excluded)
    q1 = ranked[0][0]
    singleton_margin = F(3, 7) * h - q1
    require(singleton_margin > 0, f"H4 singleton gate fails: {body}")
    level = (h - q1) / 4
    eta = level - h / 7
    require(eta > 0, f"H4 level is not finite: {body}")
    cutoff = C.ceiling(C.H.S2 * direct_r / (7 * eta)) - 1
    rows = C.H.T.coverages_many(
        carrier,
        [
            speed
            for speed in range(C.H.FIRST_EXTERNAL, cutoff + 1)
            if speed not in excluded
        ],
    )
    core = tuple(
        sorted(speed for value, speed in rows if value >= level)
    )
    require(len(core) >= 2, f"H4 core too small: {body}")
    return {
        "body": body,
        "rank": rank,
        "apex": apex,
        "prefix": prefix,
        "excluded": excluded,
        "carrier": carrier,
        "h": h,
        "r": direct_r,
        "q1": q1,
        "singleton_margin": singleton_margin,
        "level": level,
        "cutoff": cutoff,
        "H4": core,
    }


def child_profile(
    parent: dict[str, object],
    hpair: tuple[int, int],
) -> dict[str, object]:
    residual = C.H.R.subtract_local_multi(parent["carrier"], hpair)
    h = C.mass(residual)
    require(h > 0, f"H4 pair covers parent: {parent['body']}, {hpair}")
    family = tuple(
        sorted((*parent["body"], parent["apex"], *hpair))
    )
    direct, direct_r, direct_h = C.H.T.CORE.good_norm(family)
    require(
        residual == direct and len(residual) == direct_r and h == direct_h,
        f"child literal/direct mismatch: {family}",
    )
    excluded = set(parent["excluded"]) | set(hpair)
    pair = C.H.pair_cap(residual, excluded)
    top3 = tuple(pair["ranked"][:3])
    require(
        pair["tail_single"] <= top3[2][0],
        f"child rank three did not seal: {family}",
    )
    q3 = top3[2][0]
    top3_margin = h - sum((value for value, _ in top3), F(0))
    rank_pair_margin = h - q3 - pair["cap"]
    return {
        "body": parent["body"],
        "rank": parent["rank"],
        "apex": parent["apex"],
        "hpair": hpair,
        "family": family,
        "residual": residual,
        "excluded": excluded,
        "pair_data": pair,
        "h": h,
        "r": direct_r,
        "top3": top3,
        "q3": q3,
        "b2": pair["cap"],
        "b2_head": pair["head"],
        "b2_tail": pair["tail"],
        "b2_max": pair["maximizer"],
        "top3_margin": top3_margin,
        "rank_pair_margin": rank_pair_margin,
        "closed": top3_margin > 0 or rank_pair_margin > 0,
        "paid": pair["paid"],
        "paid_digest": pair["paid_digest"],
    }


def recursive_k3_close(row: dict[str, object]) -> dict[str, object]:
    """Apply THM-2895's fresh 3->1 parity descent to one hard child."""

    residual = row["residual"]
    residual_mass = row["h"]
    pair_data = row["pair_data"]
    q1 = pair_data["q1"]
    singleton_margin = F(5, 7) * residual_mass - q1
    require(
        singleton_margin > 0,
        f"recursive singleton cap reaches 5h/7: {row['family']}",
    )
    theta = residual_mass - q1
    level = theta / 2
    require(
        level > residual_mass / 7,
        f"recursive H2 is not finite: {row['family']}",
    )
    tail_first = max(
        C.H.FIRST_EXTERNAL,
        C.ceiling(C.H.S2 * row["r"] / (7 * (level - residual_mass / 7))),
    )
    by_speed = {
        speed: value for value, speed in pair_data["ranked"]
    }
    if tail_first > C.H.HORIZON + 1:
        by_speed.update(
            {
                speed: value
                for value, speed in C.H.T.coverages_many(
                    residual,
                    [
                        speed
                        for speed in range(C.H.HORIZON + 1, tail_first)
                        if speed not in row["excluded"]
                    ],
                )
            }
        )
    require(
        residual_mass / 7
        + C.H.S2 * row["r"] / (7 * tail_first)
        <= level,
        f"recursive H2 tail did not seal: {row['family']}",
    )
    core = tuple(
        sorted(
            speed
            for speed, value in by_speed.items()
            if speed not in row["excluded"] and value >= level
        )
    )

    heavy: list[tuple[int, int]] = []
    heavy_residuals: list[list[tuple[F, F]]] = []
    heavy_digest = hashlib.sha256(
        b"LRC14/j6/two-survivor-H2-heavy/v1\n"
    )
    for edge in combinations(core, 2):
        after = C.H.R.subtract_local_multi(residual, edge)
        union = residual_mass - C.mass(after)
        heavy_digest.update(
            (
                f"E={edge[0]},{edge[1]};U={ftext(union)};"
                f"L={ftext(C.mass(after))};r={len(after)}\n"
            ).encode()
        )
        if union >= theta:
            heavy.append(edge)
            heavy_residuals.append(after)

    checks = 0
    covers: list[tuple[tuple[int, int], int]] = []
    horizons: list[int] = []
    longest: list[F] = []
    singleton_digest = hashlib.sha256(
        b"LRC14/j6/two-survivor-H2-singleton/v1\n"
    )
    for edge, after in zip(heavy, heavy_residuals):
        after_mass = C.mass(after)
        require(
            after_mass > 0,
            f"recursive H2 edge covers residual: {row['family']}, {edge}",
        )
        family = tuple(sorted((*row["family"], *edge)))
        direct, direct_r, direct_h = C.H.T.CORE.good_norm(family)
        require(
            after == direct and len(after) == direct_r and after_mass == direct_h,
            f"recursive literal/direct mismatch: {family}",
        )
        ell = max(right - left for left, right in after)
        longest.append(ell)
        horizon_fraction = F(1, 7) / ell
        horizon = horizon_fraction.numerator // horizon_fraction.denominator
        horizons.append(horizon)
        edge_excluded = set(row["excluded"]) | set(edge)
        for speed in range(C.H.FIRST_EXTERNAL, horizon + 1):
            if speed in edge_excluded:
                continue
            value = C.H.T.coverage(after, speed)
            survivor = C.H.R.subtract_local(after, speed)
            survivor_mass = C.mass(survivor)
            require(
                survivor_mass == after_mass - value,
                f"recursive singleton subtraction mismatch: {family}, {speed}",
            )
            singleton_digest.update(
                (
                    f"E={edge[0]},{edge[1]};w={speed};"
                    f"c={ftext(value)};L={ftext(survivor_mass)}\n"
                ).encode()
            )
            checks += 1
            if value == after_mass:
                covers.append((edge, speed))
        require(
            F(1, 7 * (horizon + 1)) < ell,
            f"recursive geometric horizon failed: {family}",
        )
    return {
        "singleton_margin": singleton_margin,
        "theta": theta,
        "level": level,
        "tail_first": tail_first,
        "H2": core,
        "heavy": tuple(heavy),
        "heavy_digest": heavy_digest.hexdigest(),
        "longest": tuple(longest),
        "horizons": tuple(horizons),
        "checks": checks,
        "covers": tuple(covers),
        "singleton_digest": singleton_digest.hexdigest(),
        "closed": not covers,
    }


def ledger_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};rank={row['rank']};a={row['apex']};"
        f"L={row['hpair']};F={row['family']};h={ftext(row['h'])};r={row['r']};"
        f"q3={ftext(row['q3'])};B2={ftext(row['b2'])};"
        f"Bhead={ftext(row['b2_head'])};Btail={ftext(row['b2_tail'])};"
        f"Bmax={row['b2_max']};top3margin={ftext(row['top3_margin'])};"
        f"q3B2margin={ftext(row['rank_pair_margin'])};"
        f"closed={int(row['closed'])};paid={row['paid']};"
        f"paid_digest={row['paid_digest']};top3="
        + ",".join(f"{speed}:{ftext(value)}" for value, speed in row["top3"])
        + "\n"
    )


def recursive_ledger_line(
    row: dict[str, object],
    recursive: dict[str, object],
) -> str:
    return (
        f"R:E={row['body']};a={row['apex']};L={row['hpair']};"
        f"d1={ftext(recursive['singleton_margin'])};"
        f"theta={ftext(recursive['theta'])};"
        f"level={ftext(recursive['level'])};"
        f"Htail={recursive['tail_first']};H2={recursive['H2']};"
        f"heavy={recursive['heavy']};"
        f"ell={tuple(ftext(value) for value in recursive['longest'])};"
        f"W={recursive['horizons']};checks={recursive['checks']};"
        f"covers={recursive['covers']};"
        f"heavy_digest={recursive['heavy_digest']};"
        f"singleton_digest={recursive['singleton_digest']}\n"
    )


def main() -> None:
    parents = [
        parent_branch(body, rank, apex) for body, rank, apex in SURVIVORS
    ]
    local_children = [
        tuple(
            child_profile(parent, hpair)
            for hpair in combinations(parent["H4"], 2)
        )
        for parent in parents
    ]
    children = [row for local in local_children for row in local]
    failures = [row for row in children if not row["closed"]]
    recursive_rows = [
        (row, recursive_k3_close(row)) for row in failures
    ]
    recursive_failures = [
        (row, recursive)
        for row, recursive in recursive_rows
        if not recursive["closed"]
    ]
    root_rows = tuple(
        (
            parent["body"],
            parent["apex"],
            parent["h"],
            parent["q1"],
            parent["singleton_margin"],
            parent["level"],
            parent["cutoff"],
            parent["H4"],
            len(local),
            sum(row["top3_margin"] > 0 for row in local),
            sum(row["rank_pair_margin"] > 0 for row in local),
            sum(row["closed"] for row in local),
        )
        for parent, local in zip(parents, local_children)
    )
    totals = (
        len(parents),
        len(children),
        sum(row["top3_margin"] > 0 for row in children),
        sum(row["rank_pair_margin"] > 0 for row in children),
        sum(row["closed"] for row in children),
        len(failures),
        sum(row["paid"] for row in children),
        sum(len(recursive["H2"]) for _, recursive in recursive_rows),
        sum(len(recursive["heavy"]) for _, recursive in recursive_rows),
        sum(recursive["checks"] for _, recursive in recursive_rows),
        len(recursive_failures),
    )
    minimum = min(
        children,
        key=lambda row: (
            row["rank_pair_margin"],
            row["body"],
            row["hpair"],
        ),
    )
    minimum_tuple = (
        minimum["body"],
        minimum["apex"],
        minimum["hpair"],
        ftext(minimum["rank_pair_margin"]),
        ftext(minimum["top3_margin"]),
        minimum["b2_max"],
    )
    digest = hashlib.sha256(b"LRC14/j6/two-survivor-H4-q3B2/v1\n")
    for row in children:
        digest.update(ledger_line(row).encode())
    for row, recursive in recursive_rows:
        digest.update(recursive_ledger_line(row, recursive).encode())
    digest_text = digest.hexdigest()
    if EXPECTED_TOTALS is not None:
        require(totals == EXPECTED_TOTALS, "aggregate totals changed")
    if EXPECTED_ROOT_ROWS is not None:
        require(root_rows == EXPECTED_ROOT_ROWS, "parent rows changed")
    if EXPECTED_MINIMUM is not None:
        require(minimum_tuple == EXPECTED_MINIMUM, "minimum row changed")
    if EXPECTED_LEDGER_DIGEST is not None:
        require(digest_text == EXPECTED_LEDGER_DIGEST, "ledger changed")

    print("LRC14 j6 two-survivor H4 child q3+B2 probe")
    print(f"totals={totals}")
    print(f"parent_rows={root_rows}")
    print(f"minimum={minimum_tuple}")
    print(f"failure_count={len(failures)}")
    for row, recursive in recursive_rows:
        print(
            f"FAIL E={row['body']};a={row['apex']};L={row['hpair']};"
            f"top3={ftext(row['top3_margin'])};"
            f"q3B2={ftext(row['rank_pair_margin'])};"
            f"Bmax={row['b2_max']};"
            f"H2={recursive['H2']};heavy={recursive['heavy']};"
            f"ell={tuple(ftext(value) for value in recursive['longest'])};"
            f"W={recursive['horizons']};checks={recursive['checks']};"
            f"covers={recursive['covers']}"
        )
    print(f"recursive_failure_count={len(recursive_failures)}")
    print(f"canonical_ledger_sha256={digest_text}")
    print("scope=two marked first-apex branches;651 H4 children")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
