#!/usr/bin/env python3
"""Exact THM-2897 Hunter-star scalar envelope on all K>=20 branches."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CENSUS_PATH = (
    ROOT / ".scratch/lrc_sparse_j6_gate_audit_20260729/"
    "all_root_q5_pair_suffix_census.py"
)
CENSUS_SHA = "6472d6e52333a4e08688d220e53a383c84f3ec924e369e4c76977d5b2e2e7e34"
EXPECTED_COUNTS: tuple[int, ...] | None = (
    62,
    1289,
    840,
    888,
    895,
    394,
    921,
    26,
    921,
    368,
    0,
)
EXPECTED_PER_K: tuple[tuple[int, tuple[int, ...]], ...] | None = (
    (20, (600, 420, 432, 12, 168)),
    (21, (420, 286, 293, 7, 127)),
    (22, (198, 140, 145, 5, 53)),
    (23, (46, 30, 31, 1, 15)),
    (25, (25, 19, 20, 1, 5)),
)
EXPECTED_EXTREMA: tuple[tuple[object, ...], ...] | None = (
    (
        (1, 2, 7, 8, 10, 11, 12),
        20,
        7,
        38,
        "1279/709320920",
        "3238437679/14895739320",
        "76382921/1752439920",
    ),
    (
        (2, 8, 9, 10, 11, 12, 14),
        21,
        6,
        52,
        "-107/2802800",
        "5101363/25225200",
        "1050893/25225200",
    ),
    (
        (1, 4, 9, 10, 12, 13, 14),
        20,
        1,
        22,
        "-82927571/1315073760",
        "960424609/3945221280",
        "2423647/46414368",
    ),
    (
        (1, 4, 9, 10, 12, 13, 14),
        20,
        1,
        22,
        "-82927571/1315073760",
        "960424609/3945221280",
        "2423647/46414368",
    ),
)
EXPECTED_FLAG_COUNTS: tuple[int, ...] | None = (366, 282, 132, 366, 2)
EXPECTED_FLAG_HISTOGRAM: tuple[
    tuple[tuple[bool, bool, bool], int], ...
] | None = (
    ((False, False, False), 2),
    ((True, False, False), 84),
    ((True, True, False), 150),
    ((True, True, True), 132),
)
EXPECTED_NOFLAG_SURVIVORS: tuple[tuple[object, ...], ...] | None = (
    (
        (1, 3, 7, 8, 10, 11, 13),
        21,
        1,
        18,
        "-26374006621/543454231320",
        "108729775927/543454231320",
        "144619549/3123300180",
        (
            "-1848989/5703417720",
            "-1265297389/65589303780",
            "-125320907/2851708860",
        ),
    ),
    (
        (1, 4, 9, 10, 12, 13, 14),
        20,
        1,
        22,
        "-82927571/1315073760",
        "960424609/3945221280",
        "2423647/46414368",
        (
            "-157867/116035920",
            "-58997/2109744",
            "-3147959/58017960",
        ),
    ),
)
EXPECTED_MATCHING_COUNTS: tuple[int, ...] | None = (
    394,
    254,
    140,
    4,
    136,
    60,
    25,
    1,
    35,
    61,
    333,
    0,
    0,
    9871,
    1188,
    2122,
)
EXPECTED_MATCHING_STATUS: tuple[tuple[str, int], ...] | None = (
    ("ENUMERATED", 136),
    ("NONPOSITIVE", 254),
    ("TOO_LARGE", 4),
)
EXPECTED_FULLY_CLOSED_ROOTS: tuple[tuple[int, ...], ...] | None = ()
EXPECTED_SURVIVOR_DIGEST: str | None = (
    "26812ec972338364aee7ad2c80785c7e6f12a8eb5e66231952a74f880dd30588"
)
EXPECTED_COMBINED_SURVIVOR_DIGEST: str | None = (
    "6ebdf3ada4d2c4275baf12cca60324590aff9d30fcd22a180ecc173a46805af7"
)
EXPECTED_LEDGER_DIGEST: str | None = (
    "71eac8b09cd7aeeab87b213acbc423568cf4de275b5dc4e5584834b3097d0453"
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


C = load(CENSUS_PATH, CENSUS_SHA, "j6_highK_hunter_star")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def hunter_star(qs: tuple[F, ...], b2: F) -> tuple[F, F]:
    """Return G5 and a deterministic maximizing scalar breakpoint."""

    require(len(qs) == 5, "Hunter-star envelope needs five ranks")
    require(
        all(qs[index] >= qs[index + 1] for index in range(4)),
        "singleton ranks are not decreasing",
    )
    upper = min(qs[0], b2)
    candidates = {F(0), upper, b2 / 2}
    for q in qs[1:]:
        candidates.add(q)
        candidates.add(b2 - q)
    candidates = {a for a in candidates if 0 <= a <= upper}
    require(candidates, "empty Hunter-star breakpoint set")

    def invoice(a: F) -> F:
        return a + sum((min(a, q, b2 - a) for q in qs[1:]), F(0))

    value, maximizer = max((invoice(a), a) for a in candidates)
    return value, maximizer


def matching_repair(row: dict[str, object], max_cutoff: int = 10_000):
    """Replay the prior q5+M_(2,2) heavy-edge certificate from this profile."""

    if row["union_closed"]:
        return {
            "matching_status": "PRIOR_CLOSED",
            "matching_delta": None,
            "matching_cutoff": None,
            "matching_candidates": 0,
            "matching_edges": 0,
            "matching_hostile_pairs": 0,
            "hostile_matching": None,
            "matching_closed": False,
            "matching_digest": "-",
        }
    h = row["h"]
    qs = tuple(value for value, _ in row["top5"])
    b2 = row["b2"]
    target = h - qs[4]
    edge_floor = target - b2
    gamma = C.H.S2 * row["r"] / 7
    delta = edge_floor - qs[0] - h / 7
    base = {
        "matching_status": "NONPOSITIVE",
        "matching_delta": delta,
        "matching_cutoff": None,
        "matching_candidates": 0,
        "matching_edges": 0,
        "matching_hostile_pairs": 0,
        "hostile_matching": None,
        "matching_closed": False,
        "matching_digest": "-",
    }
    if delta <= 0:
        return base
    cutoff = C.ceiling(gamma / delta) - 1
    require(cutoff >= C.H.FIRST_EXTERNAL, "matching cutoff below suffix")
    base["matching_cutoff"] = cutoff
    if cutoff > max_cutoff:
        base["matching_status"] = "TOO_LARGE"
        return base

    excluded = set(row["prefix"])
    by_speed = {speed: value for value, speed in row["ranked"]}
    if cutoff >= row["tail_first"]:
        by_speed.update(
            {
                speed: value
                for value, speed in C.H.T.coverages_many(
                    row["carrier"],
                    [
                        speed
                        for speed in range(row["tail_first"], cutoff + 1)
                        if speed not in excluded
                    ],
                )
            }
        )
    vertices = sorted(
        (
            (by_speed[speed], speed)
            for speed in range(C.H.FIRST_EXTERNAL, cutoff + 1)
            if speed not in excluded
        ),
        key=lambda item: (-item[0], item[1]),
    )
    edges: list[tuple[F, int, int]] = []
    candidate_pairs = 0
    edge_hash = hashlib.sha256(b"LRC14/j6/highK-matching-heavy-edges/v1\n")
    for first_index, (first_value, first) in enumerate(vertices[:-1]):
        if first_value + vertices[first_index + 1][0] < edge_floor:
            break
        after_first = C.H.R.subtract_local(row["carrier"], first)
        for second_value, second in vertices[first_index + 1 :]:
            if first_value + second_value < edge_floor:
                break
            candidate_pairs += 1
            survivor = C.H.R.subtract_local(after_first, second)
            weight = h - C.mass(survivor)
            if weight >= edge_floor:
                x, y = sorted((first, second))
                edges.append((weight, x, y))
                edge_hash.update(f"E={x},{y};W={ftext(weight)}\n".encode())
    edges.sort(key=lambda edge: (-edge[0], edge[1], edge[2]))
    hostile_pairs = 0
    hostile_matching = None
    for first_index, first_edge in enumerate(edges[:-1]):
        if first_edge[0] + edges[first_index + 1][0] < target:
            break
        for second_edge in edges[first_index + 1 :]:
            total = first_edge[0] + second_edge[0]
            if total < target:
                break
            hostile_pairs += 1
            if {first_edge[1], first_edge[2]}.isdisjoint(
                {second_edge[1], second_edge[2]}
            ):
                witness = (
                    total,
                    (first_edge[1], first_edge[2]),
                    (second_edge[1], second_edge[2]),
                )
                if hostile_matching is None or witness > hostile_matching:
                    hostile_matching = witness
    base.update(
        {
            "matching_status": "ENUMERATED",
            "matching_candidates": candidate_pairs,
            "matching_edges": len(edges),
            "matching_hostile_pairs": hostile_pairs,
            "hostile_matching": hostile_matching,
            "matching_closed": hostile_matching is None,
            "matching_digest": edge_hash.hexdigest(),
        }
    )
    return base


def profile_body(body: tuple[int, ...]) -> tuple[dict[str, object], ...] | None:
    item = C.profile_body_star((body, 20, "all"))
    if item is None:
        return None
    _, rows = item
    profiled: list[dict[str, object]] = []
    for row in rows:
        matching = matching_repair(row)
        qs = tuple(value for value, _ in row["top5"])
        b2 = row["b2"]
        require(b2 is not None, "pair cap missing in all-pair mode")
        g5, maximizing_a = hunter_star(qs, b2)
        require(
            g5 <= sum(qs, F(0)),
            f"Hunter-star lost scalar dominance: {row['body']}, {row['rank']}",
        )
        require(
            g5 <= qs[4] + 2 * b2,
            f"Hunter-star lost pair-ladder dominance: {row['body']}, {row['rank']}",
        )
        margin = row["h"] - g5
        flag_margins = (
            F(4, 7) * row["h"] - b2,
            F(5, 7) * row["h"] - qs[2] - b2,
            F(6, 7) * row["h"] - 2 * b2,
        )
        profiled.append(
            {
                **row,
                **matching,
                "qs": qs,
                "g5": g5,
                "maximizing_a": maximizing_a,
                "g5_margin": margin,
                "g5_closed": margin > 0,
                "flag_margins": flag_margins,
                "flag_bits": tuple(value > 0 for value in flag_margins),
            }
        )
    return tuple(profiled)


def ledger_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};a={row['apex']};"
        f"P={row['prefix']};h={ftext(row['h'])};"
        f"q={tuple(ftext(q) for q in row['qs'])};B2={ftext(row['b2'])};"
        f"G5={ftext(row['g5'])};astar={ftext(row['maximizing_a'])};"
        f"margin={ftext(row['g5_margin'])};"
        f"flags={row['flag_bits']};"
        f"Mstatus={row['matching_status']};Mdelta="
        f"{ftext(row['matching_delta']) if row['matching_delta'] is not None else '-'};"
        f"MN={row['matching_cutoff']};Mclosed={int(row['matching_closed'])};"
        f"Mmatch={row['hostile_matching']};Mdigest={row['matching_digest']};"
        f"S={ftext(row['scalar_margin'])};Q={ftext(row['pair_margin'])};"
        f"prior={int(row['union_closed'])};G={int(row['g5_closed'])}\n"
    )


def row_key(row: dict[str, object]) -> tuple[object, ...]:
    return (row["body"], row["rank"], row["apex"])


def extremum_tuple(row: dict[str, object]) -> tuple[object, ...]:
    return (
        row["body"],
        row["K"],
        row["rank"],
        row["apex"],
        ftext(row["g5_margin"]),
        ftext(row["g5"]),
        ftext(row["maximizing_a"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers", type=int, default=min(os.cpu_count() or 1, 4)
    )
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    bodies = tuple(C.A.BODIES)
    if args.workers == 1:
        raw = [profile_body(body) for body in bodies]
    else:
        context = mp.get_context("spawn")
        with context.Pool(args.workers) as pool:
            raw = list(pool.imap(profile_body, bodies, chunksize=4))
    selected = [local for local in raw if local is not None]
    rows = [row for local in selected for row in local]
    require(len(selected) == 62 and len(rows) == 1289, "high-K universe changed")
    residual = [row for row in rows if not row["union_closed"]]
    survivors = [
        row for row in residual if not row["g5_closed"]
    ]
    gains = [row for row in residual if row["g5_closed"]]
    g5_closed = [row for row in rows if row["g5_closed"]]
    g5_open = [row for row in rows if not row["g5_closed"]]
    counts = (
        len(selected),
        len(rows),
        sum(row["scalar_margin"] > 0 for row in rows),
        sum(row["pair_margin"] > 0 for row in rows),
        sum(row["union_closed"] for row in rows),
        len(residual),
        len(g5_closed),
        len(gains),
        sum(row["union_closed"] or row["g5_closed"] for row in rows),
        len(survivors),
        sum(row["union_closed"] and not row["g5_closed"] for row in rows),
    )
    per_k = tuple(
        (
            k,
            (
                sum(row["K"] == k for row in rows),
                sum(row["K"] == k and row["union_closed"] for row in rows),
                sum(row["K"] == k and row["g5_closed"] for row in rows),
                sum(
                    row["K"] == k
                    and not row["union_closed"]
                    and row["g5_closed"]
                    for row in rows
                ),
                sum(row["K"] == k for row in survivors),
            ),
        )
        for k in sorted(Counter(row["K"] for row in rows))
    )
    closest_closed = min(
        g5_closed, key=lambda row: (row["g5_margin"], row_key(row))
    )
    closest_open = max(
        g5_open, key=lambda row: (row["g5_margin"], row_key(row))
    )
    worst = min(rows, key=lambda row: (row["g5_margin"], row_key(row)))
    worst_residual = min(
        residual, key=lambda row: (row["g5_margin"], row_key(row))
    )
    extrema = tuple(
        extremum_tuple(row)
        for row in (closest_closed, closest_open, worst, worst_residual)
    )
    flag_counts = (
        *(sum(row["flag_bits"][index] for row in survivors) for index in range(3)),
        sum(any(row["flag_bits"]) for row in survivors),
        sum(not any(row["flag_bits"]) for row in survivors),
    )
    flag_histogram = tuple(
        sorted(Counter(row["flag_bits"] for row in survivors).items())
    )
    noflag_survivors = tuple(
        (
            *extremum_tuple(row),
            tuple(ftext(value) for value in row["flag_margins"]),
        )
        for row in survivors
        if not any(row["flag_bits"])
    )
    matching_status = tuple(
        sorted(Counter(row["matching_status"] for row in residual).items())
    )
    matching_closed = [row for row in residual if row["matching_closed"]]
    matching_enumerated = [
        row for row in residual if row["matching_status"] == "ENUMERATED"
    ]
    combined_survivors = [
        row
        for row in residual
        if not row["g5_closed"] and not row["matching_closed"]
    ]
    g5_only = [
        row for row in residual if row["g5_closed"] and not row["matching_closed"]
    ]
    matching_only = [
        row for row in residual if not row["g5_closed"] and row["matching_closed"]
    ]
    survivor_bodies = {row["body"] for row in combined_survivors}
    fully_closed_roots = tuple(
        sorted(
            local[0]["body"]
            for local in selected
            if local[0]["body"] not in survivor_bodies
        )
    )
    delta_nonpositive = [
        row for row in residual if row["matching_status"] == "NONPOSITIVE"
    ]
    too_large = [
        row for row in residual if row["matching_status"] == "TOO_LARGE"
    ]
    matching_counts = (
        len(residual),
        len(delta_nonpositive),
        sum(row["matching_delta"] is not None and row["matching_delta"] > 0 for row in residual),
        len(too_large),
        len(matching_enumerated),
        len(matching_closed),
        sum(row["g5_closed"] and row["matching_closed"] for row in residual),
        sum(row["g5_closed"] and not row["matching_closed"] for row in residual),
        sum(not row["g5_closed"] and row["matching_closed"] for row in residual),
        len(residual) - len(combined_survivors),
        len(combined_survivors),
        sum(row["g5_closed"] for row in delta_nonpositive),
        sum(row["g5_closed"] for row in too_large),
        sum(row["matching_candidates"] for row in matching_enumerated),
        sum(row["matching_edges"] for row in matching_enumerated),
        sum(row["matching_hostile_pairs"] for row in matching_enumerated),
    )
    survivor_hash = hashlib.sha256(
        b"LRC14/j6/highK-Hunter-star-survivors/v1\n"
    )
    for row in survivors:
        survivor_hash.update(ledger_line(row).encode())
    survivor_digest = survivor_hash.hexdigest()
    combined_hash = hashlib.sha256(
        b"LRC14/j6/highK-Hunter-star-matching-survivors/v1\n"
    )
    for row in combined_survivors:
        combined_hash.update(ledger_line(row).encode())
    combined_digest = combined_hash.hexdigest()
    digest = hashlib.sha256(b"LRC14/j6/highK-Hunter-star/v1\n")
    for row in rows:
        digest.update(ledger_line(row).encode())
    digest_text = digest.hexdigest()
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "Hunter-star counts changed")
    if EXPECTED_PER_K is not None:
        require(per_k == EXPECTED_PER_K, "Hunter-star K strata changed")
    if EXPECTED_EXTREMA is not None:
        require(extrema == EXPECTED_EXTREMA, "Hunter-star extrema changed")
    if EXPECTED_FLAG_COUNTS is not None:
        require(flag_counts == EXPECTED_FLAG_COUNTS, "flag counts changed")
    if EXPECTED_FLAG_HISTOGRAM is not None:
        require(
            flag_histogram == EXPECTED_FLAG_HISTOGRAM,
            "flag histogram changed",
        )
    if EXPECTED_NOFLAG_SURVIVORS is not None:
        require(
            noflag_survivors == EXPECTED_NOFLAG_SURVIVORS,
            "nonflag survivors changed",
        )
    if EXPECTED_MATCHING_COUNTS is not None:
        require(
            matching_counts == EXPECTED_MATCHING_COUNTS,
            "matching overlap counts changed",
        )
    if EXPECTED_MATCHING_STATUS is not None:
        require(
            matching_status == EXPECTED_MATCHING_STATUS,
            "matching status distribution changed",
        )
    if EXPECTED_FULLY_CLOSED_ROOTS is not None:
        require(
            fully_closed_roots == EXPECTED_FULLY_CLOSED_ROOTS,
            "fully closed root set changed",
        )
    if EXPECTED_SURVIVOR_DIGEST is not None:
        require(
            survivor_digest == EXPECTED_SURVIVOR_DIGEST,
            "Hunter-star survivor ledger changed",
        )
    if EXPECTED_COMBINED_SURVIVOR_DIGEST is not None:
        require(
            combined_digest == EXPECTED_COMBINED_SURVIVOR_DIGEST,
            "combined survivor ledger changed",
        )
    if EXPECTED_LEDGER_DIGEST is not None:
        require(digest_text == EXPECTED_LEDGER_DIGEST, "Hunter-star ledger changed")

    print("LRC14 j6 K>=20 Hunter-star G5 envelope census")
    print(f"parameters=workers:{args.workers}")
    print(f"counts={counts}")
    print(f"per_K={per_k}")
    print(f"extrema={extrema}")
    print(f"survivor_flag_counts_H3_H2_H1_union_none={flag_counts}")
    print(f"survivor_flag_histogram={flag_histogram}")
    print(f"nonflag_survivors={noflag_survivors}")
    print(f"matching_overlap_counts={matching_counts}")
    print(f"matching_status={matching_status}")
    print(f"fully_closed_roots={fully_closed_roots}")
    print(f"G5_only_rows={tuple(extremum_tuple(row) for row in g5_only)}")
    print(
        "matching_only_rows="
        + repr(tuple(extremum_tuple(row) for row in matching_only))
    )
    print(f"gain_rows={tuple(extremum_tuple(row) for row in gains)}")
    print(
        "survivor_rows="
        + repr(tuple(extremum_tuple(row) for row in survivors))
    )
    print(f"survivor_ledger_sha256={survivor_digest}")
    print(f"combined_survivor_ledger_sha256={combined_digest}")
    print(f"canonical_ledger_sha256={digest_text}")
    print("scope=K>=20 exact scalar Hunter-star envelope;certificate only")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
