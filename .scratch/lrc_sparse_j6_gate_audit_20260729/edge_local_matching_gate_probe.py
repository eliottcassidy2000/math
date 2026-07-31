#!/usr/bin/env python3
"""Edge-local residual cutoff for THM-2897 matching repair.

For an L-heavy edge {x,y}, at least one endpoint has c(x)>=L/2.  Enumerate
that finite high core.  After fixing x, use the literal residual C\\D_x:

    U_C({x,y}) = c(x) + c_(C\\D_x)(y).

The residual discrepancy cutoff is often finite when the coarser parent
singleton cutoff is not.
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
MATCHING_PATH = (
    ROOT / ".scratch/lrc_sparse_j6_gate_audit_20260729/"
    "highK_matching_repair_probe.py"
)
MATCHING_SHA = "52e409f24610b6241af72885a600a6242bb89809fa77cd348292e706545b3149"


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


M = load(MATCHING_PATH, MATCHING_SHA, "j6_edge_local_matching")
C = M.C
EXPECTED_LOCKED_COUNTS = (
    62,
    1_289,
    394,
    254,
    76,
    0,
    76,
    14_603,
    11_548,
    36_843,
)
EXPECTED_LOCKED_STATUS = (
    ("ENUMERATED", 76),
    ("HIGH_CUTOFF_TOO_LARGE", 6),
    ("HIGH_LEVEL_NONPOSITIVE", 59),
    ("PARTNER_GATE_FAILED", 113),
)
EXPECTED_LOCKED_DIGEST = (
    "d113e4e9a5037dc846eb5259c59b1b8188a44ec9e9a7ee3c8e9b5cd6d20ee5e7"
)
EXPECTED_HOSTILE_CONTROL = (
    4,
    46,
    1_115,
    15,
    1_746,
    84,
    313,
    (F(210257, 1381380), (34, 44), (17, 18)),
    False,
    "dab2481bc344cf9e16a04633a43e55156fef41894b2a98648a68b64aae227d85",
)
EXPECTED_POSITIVE_CONTROL = (
    6,
    17,
    682,
    2,
    696,
    3,
    2,
    None,
    True,
    "6c0105d477d23e27f3cff4ed4f4521dbe896ad4340ac9208e18713b5da207bca",
)


def ftext(value: F | None) -> str:
    return M.ftext(value)


def edge_local_graph(
    root: dict[str, object],
    root_good: list[tuple[F, F]],
    rank: int,
    max_high_cutoff: int,
    max_partner_cutoff: int,
    parent_delta: str,
) -> dict[str, object]:
    base = M.matching_repair(root, root_good, rank, 15)
    if base["union_closed"]:
        return base
    if parent_delta == "nonpositive" and base["delta"] > 0:
        base.update(
            {
                "local_status": "PARENT_POSITIVE_SKIPPED",
                "high_level": None,
                "high_eta": None,
                "high_cutoff": None,
                "high_core": (),
                "bad_high": (),
                "max_partner_cutoff": None,
                "local_candidate_pairs": 0,
                "local_edges": (),
                "local_hostile_pairs": 0,
                "local_hostile_matching": None,
                "local_closed": False,
                "local_edge_digest": "-",
            }
        )
        return base
    body = root["body"]
    apex = base["apex"]
    excluded = set(base["prefix"])
    carrier = C.H.R.subtract_local(root_good, apex)
    h = C.mass(carrier)
    level = base["edge_floor"] / 2
    gamma = C.H.S2 * base["r"] / 7
    eta = level - h / 7
    base.update(
        {
            "local_status": "HIGH_LEVEL_NONPOSITIVE",
            "high_level": level,
            "high_eta": eta,
            "high_cutoff": None,
            "high_core": (),
            "bad_high": (),
            "max_partner_cutoff": None,
            "local_candidate_pairs": 0,
            "local_edges": (),
            "local_hostile_pairs": 0,
            "local_hostile_matching": None,
            "local_closed": False,
            "local_edge_digest": "-",
        }
    )
    if eta <= 0:
        return base
    high_cutoff = C.ceiling(gamma / eta) - 1
    base["high_cutoff"] = high_cutoff
    if high_cutoff > max_high_cutoff:
        base["local_status"] = "HIGH_CUTOFF_TOO_LARGE"
        return base

    high_rows = C.H.T.coverages_many(
        carrier,
        [
            speed
            for speed in range(C.H.FIRST_EXTERNAL, high_cutoff + 1)
            if speed not in excluded
        ],
    )
    high_core = tuple(
        sorted(speed for value, speed in high_rows if value >= level)
    )
    by_high = {speed: value for value, speed in high_rows}
    base["high_core"] = high_core
    edge_bank: dict[tuple[int, int], F] = {}
    bad_high: list[tuple[int, F | None, int | None, str]] = []
    partner_cutoffs: list[int] = []
    candidate_pairs = 0
    for x in high_core:
        cx = by_high[x]
        residual = C.H.R.subtract_local(carrier, x)
        residual_h = C.mass(residual)
        require(residual_h == h - cx, f"edge residual mass mismatch: {body}, {x}")
        beta = base["edge_floor"] - cx - residual_h / 7
        if beta <= 0:
            bad_high.append((x, beta, None, "NONPOSITIVE"))
            continue
        partner_gamma = C.H.S2 * len(residual) / 7
        partner_cutoff = C.ceiling(partner_gamma / beta) - 1
        partner_cutoffs.append(partner_cutoff)
        if partner_cutoff > max_partner_cutoff:
            bad_high.append((x, beta, partner_cutoff, "TOO_LARGE"))
            continue
        partners = [
            speed
            for speed in range(C.H.FIRST_EXTERNAL, partner_cutoff + 1)
            if speed not in excluded and speed != x
        ]
        for residual_value, y in C.H.T.coverages_many(residual, partners):
            if cx + residual_value < base["edge_floor"]:
                continue
            candidate_pairs += 1
            edge = tuple(sorted((x, y)))
            weight = cx + residual_value
            if edge in edge_bank:
                require(
                    edge_bank[edge] == weight,
                    f"edge weight depends on chosen endpoint: {body}, {edge}",
                )
            else:
                edge_bank[edge] = weight

    base["bad_high"] = tuple(bad_high)
    base["max_partner_cutoff"] = max(partner_cutoffs, default=0)
    if bad_high:
        base["local_status"] = "PARTNER_GATE_FAILED"
        return base

    edges = sorted(
        ((weight, edge[0], edge[1]) for edge, weight in edge_bank.items()),
        key=lambda edge: (-edge[0], edge[1], edge[2]),
    )
    edge_hash = hashlib.sha256(b"LRC14/j6/edge-local-heavy-graph/v1\n")
    for weight, x, y in edges:
        edge_hash.update(f"E={x},{y};W={ftext(weight)}\n".encode())
    hostile_pairs = 0
    hostile_matching: tuple[F, tuple[int, int], tuple[int, int]] | None = None
    target = base["target"]
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
            "local_status": "ENUMERATED",
            "local_candidate_pairs": candidate_pairs,
            "local_edges": tuple(edges),
            "local_hostile_pairs": hostile_pairs,
            "local_hostile_matching": hostile_matching,
            "local_closed": hostile_matching is None,
            "local_edge_digest": edge_hash.hexdigest(),
        }
    )
    return base


def profile_body_star(
    task: tuple[tuple[int, ...], int, int, int, str],
):
    body, min_k, max_high_cutoff, max_partner_cutoff, parent_delta = task
    root = C.A.profile_body(body)
    if root["adaptive_k"] < min_k:
        return None
    root_good, root_r, root_h = C.A.T.CORE.good_norm(body)
    require(
        root_r == root["r"] and root_h == root["m"],
        f"root carrier mismatch: {body}",
    )
    return tuple(
        edge_local_graph(
            root,
            root_good,
            rank,
            max_high_cutoff,
            max_partner_cutoff,
            parent_delta,
        )
        for rank in range(1, root["adaptive_k"] + 1)
    )


def ledger_line(row: dict[str, object]) -> str:
    if row["union_closed"]:
        return (
            f"E={row['body']};K={row['K']};rank={row['rank']};U=1;"
            f"S={ftext(row['scalar_margin'])};Q={ftext(row['pair_margin'])}\n"
        )
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};a={row['apex']};"
        f"P={row['prefix']};T={ftext(row['target'])};L={ftext(row['edge_floor'])};"
        f"level={ftext(row['high_level'])};eta={ftext(row['high_eta'])};"
        f"Nh={row['high_cutoff']};H={row['high_core']};bad={row['bad_high']};"
        f"Nmax={row['max_partner_cutoff']};status={row['local_status']};"
        f"candidates={row['local_candidate_pairs']};"
        f"edges={row['local_edges']};hostile={row['local_hostile_pairs']};"
        f"matching={row['local_hostile_matching']};"
        f"closed={int(row['local_closed'])};"
        f"digest={row['local_edge_digest']}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", choices=("unique", "high"), default="unique")
    parser.add_argument("--max-high-cutoff", type=int, default=20_000)
    parser.add_argument("--max-partner-cutoff", type=int, default=20_000)
    parser.add_argument(
        "--parent-delta",
        choices=("nonpositive", "all"),
        default="nonpositive",
    )
    parser.add_argument(
        "--workers", type=int, default=min(os.cpu_count() or 1, 4)
    )
    args = parser.parse_args()
    require(
        args.max_high_cutoff >= 15
        and args.max_partner_cutoff >= 15
        and args.workers >= 1,
        "bad probe arguments",
    )
    bodies = (
        (C.UNIQUE_MAX_ROOT,)
        if args.scope == "unique"
        else C.A.BODIES
    )
    min_k = 20
    tasks = tuple(
        (
            body,
            min_k,
            args.max_high_cutoff,
            args.max_partner_cutoff,
            args.parent_delta,
        )
        for body in bodies
    )
    if args.workers == 1 or len(tasks) == 1:
        raw = [profile_body_star(task) for task in tasks]
    else:
        context = mp.get_context("spawn")
        with context.Pool(args.workers) as pool:
            raw = list(pool.imap(profile_body_star, tasks, chunksize=4))
    selected = [local for local in raw if local is not None]
    rows = [row for local in selected for row in local]
    open_rows = [row for row in rows if not row["union_closed"]]
    target_rows = [
        row
        for row in open_rows
        if row["local_status"] != "PARENT_POSITIVE_SKIPPED"
    ]
    status = tuple(
        sorted(Counter(row["local_status"] for row in target_rows).items())
    )
    enumerated = [
        row for row in target_rows if row["local_status"] == "ENUMERATED"
    ]
    closed = [row for row in enumerated if row["local_closed"]]
    matching_open = [row for row in enumerated if not row["local_closed"]]
    counts = (
        len(selected),
        len(rows),
        len(open_rows),
        len(target_rows),
        len(enumerated),
        len(closed),
        len(matching_open),
        sum(row["local_candidate_pairs"] for row in enumerated),
        sum(len(row["local_edges"]) for row in enumerated),
        sum(row["local_hostile_pairs"] for row in enumerated),
    )
    digest = hashlib.sha256(b"LRC14/j6/edge-local-matching-gate/v1\n")
    for row in rows:
        digest.update(ledger_line(row).encode())
    digest_text = digest.hexdigest()
    locked_mode = (
        args.scope == "high"
        and args.parent_delta == "nonpositive"
        and args.max_high_cutoff == 20_000
        and args.max_partner_cutoff == 20_000
    )
    controls = None
    if locked_mode:
        require(counts == EXPECTED_LOCKED_COUNTS, "locked counts changed")
        require(status == EXPECTED_LOCKED_STATUS, "locked statuses changed")
        require(digest_text == EXPECTED_LOCKED_DIGEST, "locked ledger changed")
        hostile_row = next(
            row
            for row in target_rows
            if row["body"] == C.UNIQUE_MAX_ROOT and row["rank"] == 4
        )
        unique_root = C.A.profile_body(C.UNIQUE_MAX_ROOT)
        unique_good, _, _ = C.A.T.CORE.good_norm(C.UNIQUE_MAX_ROOT)
        positive_row = edge_local_graph(
            unique_root,
            unique_good,
            6,
            20_000,
            20_000,
            "all",
        )

        def control_tuple(row: dict[str, object]) -> tuple[object, ...]:
            return (
                row["rank"],
                row["apex"],
                row["high_cutoff"],
                len(row["high_core"]),
                row["max_partner_cutoff"],
                len(row["local_edges"]),
                row["local_hostile_pairs"],
                row["local_hostile_matching"],
                row["local_closed"],
                row["local_edge_digest"],
            )

        controls = (control_tuple(hostile_row), control_tuple(positive_row))
        require(
            controls
            == (EXPECTED_HOSTILE_CONTROL, EXPECTED_POSITIVE_CONTROL),
            "hostile/positive controls changed",
        )

    print("LRC14 j6 edge-local matching cutoff probe")
    print(
        f"parameters=scope:{args.scope},max_high:{args.max_high_cutoff},"
        f"max_partner:{args.max_partner_cutoff},"
        f"parent_delta:{args.parent_delta},workers:{args.workers}"
    )
    print(f"counts={counts}")
    print(f"status_distribution={status}")
    if controls is not None:
        print(f"hostile_positive_controls={controls}")
    for row in target_rows:
        print(
            f"ROW E={row['body']};rank={row['rank']};a={row['apex']};"
            f"status={row['local_status']};eta={ftext(row['high_eta'])};"
            f"Nh={row['high_cutoff']};H={len(row['high_core'])};"
            f"bad={row['bad_high']};Nmax={row['max_partner_cutoff']};"
            f"edges={len(row['local_edges'])};"
            f"hostile={row['local_hostile_pairs']};"
            f"matching={row['local_hostile_matching']};"
            f"closed={int(row['local_closed'])}"
        )
    print(f"canonical_ledger_sha256={digest_text}")
    print("graph=undirected exact weighted graph;ties retained;no orientation")
    print("scope=edge-local finite-tail matching certificate only")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
