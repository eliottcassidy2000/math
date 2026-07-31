#!/usr/bin/env python3
"""Probe THM-2897's q5 + M_(2,2) repair on the K>=20 failure core.

The weighted graph is undirected.  Every edge {x,y} retains its exact
union-coverage weight; ties are preserved and no orientation is introduced.
"""

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
CENSUS_SHA = "02e8802eb60f6d598903b84e47f454592a78001f5143c02450533c0899e06c07"


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


C = load(CENSUS_PATH, CENSUS_SHA, "j6_highK_matching_census")


def ftext(value: F | None) -> str:
    if value is None:
        return "-"
    return f"{value.numerator}/{value.denominator}"


def matching_repair(
    root: dict[str, object],
    root_good: list[tuple[F, F]],
    rank: int,
    max_cutoff: int,
) -> dict[str, object]:
    body = root["body"]
    apex = root["top"][rank - 1][1]
    prefix = tuple(speed for _, speed in root["top"][:rank])
    excluded = set(prefix)
    carrier = C.H.R.subtract_local(root_good, apex)
    direct, direct_r, direct_m = C.H.T.CORE.good_norm(
        tuple(sorted((*body, apex)))
    )
    h = C.mass(carrier)
    require(
        carrier == direct and len(carrier) == direct_r and h == direct_m > 0,
        f"literal/direct carrier mismatch: {body}, rank {rank}",
    )
    ranked, tail_single, tail_first = C.sealed_ranking(carrier, excluded)
    top5 = tuple(ranked[:5])
    q1 = ranked[0][0]
    q5 = top5[4][0]
    scalar_margin = h - sum((value for value, _ in top5), F(0))
    pair = C.pair_cap_from_ranking(carrier, ranked, tail_single)
    b2 = pair["cap"]
    pair_margin = h - q5 - 2 * b2
    union_closed = scalar_margin > 0 or pair_margin > 0
    base = {
        "body": body,
        "K": root["adaptive_k"],
        "rank": rank,
        "apex": apex,
        "prefix": prefix,
        "h": h,
        "r": direct_r,
        "q1": q1,
        "q5": q5,
        "b2": b2,
        "scalar_margin": scalar_margin,
        "pair_margin": pair_margin,
        "union_closed": union_closed,
        "tail_first": tail_first,
    }
    if union_closed:
        return base

    target = h - q5
    edge_floor = target - b2
    gamma = C.H.S2 * direct_r / 7
    delta = edge_floor - q1 - h / 7
    base.update(
        {
            "target": target,
            "edge_floor": edge_floor,
            "delta": delta,
            "cutoff": None,
            "cutoff_status": "NONPOSITIVE",
            "candidate_pairs": 0,
            "heavy_edges": 0,
            "hostile_edge_pairs": 0,
            "hostile_matching": None,
            "matching_closed": False,
            "edge_digest": "-",
        }
    )
    if delta <= 0:
        return base

    cutoff = C.ceiling(gamma / delta) - 1
    require(cutoff >= C.H.FIRST_EXTERNAL, "positive cutoff below external range")
    base["cutoff"] = cutoff
    if cutoff > max_cutoff:
        base["cutoff_status"] = "TOO_LARGE"
        return base

    by_speed = {speed: value for value, speed in ranked}
    if cutoff >= tail_first:
        by_speed.update(
            {
                speed: value
                for value, speed in C.H.T.coverages_many(
                    carrier,
                    [
                        speed
                        for speed in range(tail_first, cutoff + 1)
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
        after_first = C.H.R.subtract_local(carrier, first)
        for second_value, second in vertices[first_index + 1 :]:
            if first_value + second_value < edge_floor:
                break
            candidate_pairs += 1
            survivor = C.H.R.subtract_local(after_first, second)
            weight = h - C.mass(survivor)
            if weight >= edge_floor:
                x, y = sorted((first, second))
                edges.append((weight, x, y))
                edge_hash.update(
                    f"E={x},{y};W={ftext(weight)}\n".encode()
                )
    edges.sort(key=lambda edge: (-edge[0], edge[1], edge[2]))
    hostile_edge_pairs = 0
    hostile_matching: tuple[F, tuple[int, int], tuple[int, int]] | None = None
    for first_index, first_edge in enumerate(edges[:-1]):
        if first_edge[0] + edges[first_index + 1][0] < target:
            break
        for second_edge in edges[first_index + 1 :]:
            total = first_edge[0] + second_edge[0]
            if total < target:
                break
            hostile_edge_pairs += 1
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
            "cutoff_status": "ENUMERATED",
            "candidate_pairs": candidate_pairs,
            "heavy_edges": len(edges),
            "hostile_edge_pairs": hostile_edge_pairs,
            "hostile_matching": hostile_matching,
            "matching_closed": hostile_matching is None,
            "edge_digest": edge_hash.hexdigest(),
        }
    )
    return base


def profile_body_star(task: tuple[tuple[int, ...], int, int]):
    body, min_k, max_cutoff = task
    root = C.A.profile_body(body)
    if root["adaptive_k"] < min_k:
        return None
    root_good, root_r, root_h = C.A.T.CORE.good_norm(body)
    require(
        root_r == root["r"] and root_h == root["m"],
        f"root carrier mismatch: {body}",
    )
    return tuple(
        matching_repair(root, root_good, rank, max_cutoff)
        for rank in range(1, root["adaptive_k"] + 1)
    )


def ledger_line(row: dict[str, object]) -> str:
    if row["union_closed"]:
        return (
            f"E={row['body']};K={row['K']};rank={row['rank']};"
            f"a={row['apex']};U=1;S={ftext(row['scalar_margin'])};"
            f"Q={ftext(row['pair_margin'])}\n"
        )
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};"
        f"a={row['apex']};P={row['prefix']};h={ftext(row['h'])};r={row['r']};"
        f"q1={ftext(row['q1'])};q5={ftext(row['q5'])};B2={ftext(row['b2'])};"
        f"T={ftext(row['target'])};L={ftext(row['edge_floor'])};"
        f"delta={ftext(row['delta'])};N={row['cutoff']};"
        f"status={row['cutoff_status']};candidates={row['candidate_pairs']};"
        f"edges={row['heavy_edges']};hostile_pairs={row['hostile_edge_pairs']};"
        f"matching={row['hostile_matching']};closed={int(row['matching_closed'])};"
        f"edge_digest={row['edge_digest']}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-k", type=int, default=20)
    parser.add_argument("--max-cutoff", type=int, default=10_000)
    parser.add_argument(
        "--workers", type=int, default=min(os.cpu_count() or 1, 4)
    )
    args = parser.parse_args()
    require(
        args.min_k >= 0 and args.max_cutoff >= 15 and args.workers >= 1,
        "bad probe arguments",
    )
    tasks = tuple(
        (body, args.min_k, args.max_cutoff) for body in C.A.BODIES
    )
    if args.workers == 1:
        raw = [profile_body_star(task) for task in tasks]
    else:
        context = mp.get_context("spawn")
        with context.Pool(args.workers) as pool:
            raw = list(pool.imap(profile_body_star, tasks, chunksize=4))
    selected = [local for local in raw if local is not None]
    rows = [row for local in selected for row in local]
    open_rows = [row for row in rows if not row["union_closed"]]
    cutoff_fail = [row for row in open_rows if row["delta"] <= 0]
    too_large = [
        row for row in open_rows if row["cutoff_status"] == "TOO_LARGE"
    ]
    enumerated = [
        row for row in open_rows if row["cutoff_status"] == "ENUMERATED"
    ]
    matching_closed = [row for row in enumerated if row["matching_closed"]]
    matching_open = [row for row in enumerated if not row["matching_closed"]]
    counts = (
        len(selected),
        len(rows),
        len(open_rows),
        len(cutoff_fail),
        len(too_large),
        len(enumerated),
        len(matching_closed),
        len(matching_open),
        sum(row["candidate_pairs"] for row in enumerated),
        sum(row["heavy_edges"] for row in enumerated),
        sum(row["hostile_edge_pairs"] for row in enumerated),
    )
    cutoff_distribution = tuple(
        sorted(
            Counter(
                row["cutoff_status"]
                for row in open_rows
            ).items()
        )
    )
    positive = [row for row in open_rows if row["delta"] > 0]
    max_cutoff_row = max(
        positive, key=lambda row: (row["cutoff"], row["body"], row["rank"])
    ) if positive else None
    closest_delta = max(
        cutoff_fail, key=lambda row: (row["delta"], row["body"], row["rank"])
    ) if cutoff_fail else None
    digest = hashlib.sha256(b"LRC14/j6/highK-matching-repair/v1\n")
    for row in rows:
        digest.update(ledger_line(row).encode())

    print("LRC14 j6 K>=20 matching-repair probe")
    print(
        f"parameters=min_k:{args.min_k},max_cutoff:{args.max_cutoff},"
        f"workers:{args.workers}"
    )
    print(f"counts={counts}")
    print(f"cutoff_distribution={cutoff_distribution}")
    if max_cutoff_row is not None:
        print(
            "maximum_positive_cutoff="
            f"{max_cutoff_row['cutoff']};E={max_cutoff_row['body']};"
            f"rank={max_cutoff_row['rank']};a={max_cutoff_row['apex']};"
            f"delta={ftext(max_cutoff_row['delta'])}"
        )
    if closest_delta is not None:
        print(
            "closest_nonpositive_delta="
            f"{ftext(closest_delta['delta'])};E={closest_delta['body']};"
            f"rank={closest_delta['rank']};a={closest_delta['apex']};"
            f"T={ftext(closest_delta['target'])};"
            f"L={ftext(closest_delta['edge_floor'])};"
            f"q1={ftext(closest_delta['q1'])}"
        )
    unique = [
        row for row in open_rows if row["body"] == C.UNIQUE_MAX_ROOT
    ]
    print(
        "unique_K25="
        + repr(
            tuple(
                (
                    row["rank"],
                    row["apex"],
                    row["cutoff_status"],
                    row["cutoff"],
                    row["matching_closed"],
                    row["hostile_matching"],
                )
                for row in unique
            )
        )
    )
    for label, group in (
        ("NONPOSITIVE", cutoff_fail),
        ("TOO_LARGE", too_large),
        ("MATCHING_OPEN", matching_open),
    ):
        print(f"{label}_count={len(group)}")
        for row in group[:40]:
            print(
                f"{label} E={row['body']};K={row['K']};"
                f"rank={row['rank']};a={row['apex']};"
                f"delta={ftext(row['delta'])};N={row['cutoff']};"
                f"T={ftext(row['target'])};L={ftext(row['edge_floor'])};"
                f"q1={ftext(row['q1'])};B2={ftext(row['b2'])};"
                f"edges={row['heavy_edges']};matching={row['hostile_matching']}"
            )
    print(f"canonical_ledger_sha256={digest.hexdigest()}")
    print("graph=undirected exact weighted graph;ties retained;no orientation")
    print("scope=K>=20 q5/top5 failure core;matching certificate only")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
