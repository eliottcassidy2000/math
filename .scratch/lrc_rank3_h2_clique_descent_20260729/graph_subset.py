#!/usr/bin/env python3
"""Run the exact heavy graph/K4 census on a sealed-core size slice."""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PILOT_PATH = ROOT / ".scratch/lrc_rank3_h2_clique_descent_20260729/pilot.py"
PILOT_SHA256 = (
    "9ec1b8d8c697f54fdbf2836638fa6c7dc5284f4ea2644849d98d71398c1f3520"
)
CORE_PATH = (
    ROOT
    / ".scratch/lrc_rank3_h2_clique_descent_20260729/core.discovery.out"
)
CORE_SHA256 = (
    "5569d8d34b59e5eed2bfa82148648dcb5e63515146ff55b95234a002d0c87ba2"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_pilot():
    require(hashlib.sha256(PILOT_PATH.read_bytes()).hexdigest() == PILOT_SHA256, "pilot changed")
    require(hashlib.sha256(CORE_PATH.read_bytes()).hexdigest() == CORE_SHA256, "core ledger changed")
    spec = importlib.util.spec_from_file_location("rank3_h2_subset_pilot", PILOT_PATH)
    require(spec is not None and spec.loader is not None, "cannot load pilot")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


P = load_pilot()


def profile_task(task):
    return P.profile_target(task)


def core_sizes() -> dict[tuple[tuple[int, ...], int], int]:
    out = {}
    for line in CORE_PATH.read_text().splitlines():
        if not line.startswith("H2;"):
            continue
        fields = dict(part.split("=", 1) for part in line.split(";")[1:])
        key = (ast.literal_eval(fields["E"]), int(fields["rank"]))
        out[key] = int(fields["H"])
    require(len(out) == 143, "core key universe changed")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-core", type=int, default=40)
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 1, 8),
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="also reconstruct every K4 residual and run singleton seal",
    )
    parser.add_argument("--emit-ledger", action="store_true")
    args = parser.parse_args()
    require(args.max_core >= 4 and args.workers >= 1, "bad parameters")
    sizes = core_sizes()
    selected = tuple(
        row
        for row in P.TARGET_ROWS
        if sizes[(row["body"], row["rank"])] <= args.max_core
    )
    tasks = tuple((row, not args.full, False) for row in selected)
    if args.workers == 1:
        rows = [profile_task(task) for task in tasks]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = list(pool.imap(profile_task, tasks, chunksize=1))
    require(
        tuple((row["body"], row["rank"]) for row in rows)
        == tuple((row["body"], row["rank"]) for row in selected),
        "worker order changed",
    )
    require(
        all(len(row["core"]) == sizes[(row["body"], row["rank"])] for row in rows),
        "replayed core size changed",
    )
    digest = hashlib.sha256(
        f"LRC14/j6/highK-r3-H2-graph/maxH={args.max_core}/v1\n".encode()
    )
    for row in rows:
        line = P.core_line(row)
        digest.update(line.encode())
        if args.emit_ledger:
            print("CORE;" + line.rstrip())
        if args.full:
            for item in row["clique_profiles"]:
                leaf = P.clique_line(row, item)
                digest.update(leaf.encode())
                if args.emit_ledger:
                    print("LEAF;" + leaf.rstrip())
    print("LRC14 j6 K>=20 strongest-r3 exact heavy-graph slice")
    print(
        f"max_core={args.max_core};branches={len(rows)};"
        f"vertices={sum(len(row['core']) for row in rows)};"
        f"pairs={sum(row['pairs'] for row in rows)};"
        f"paid={sum(row['paid_pairs'] for row in rows)};"
        f"edges={sum(len(row['edges']) for row in rows)};"
        f"ties={sum(len(row['equality_edges']) for row in rows)};"
        f"K4={sum(len(row['cliques']) for row in rows)}"
    )
    print(
        f"edge_histogram={tuple(sorted(Counter(len(row['edges']) for row in rows).items()))}"
    )
    print(
        f"degeneracy_histogram={tuple(sorted(Counter(row['degeneracy'] for row in rows).items()))}"
    )
    print(
        f"K4_histogram={tuple(sorted(Counter(len(row['cliques']) for row in rows).items()))}"
    )
    print(
        f"K4_free={sum(not row['cliques'] for row in rows)};"
        f"maximum_edges={max((len(row['edges']),row['body'],row['rank']) for row in rows)};"
        f"maximum_degeneracy={max((row['degeneracy'],row['body'],row['rank']) for row in rows)};"
        f"maximum_K4={max((len(row['cliques']),row['body'],row['rank']) for row in rows)}"
    )
    print(
        f"ordered_work=degeneracy_forward_triples:"
        f"{sum(row['forward_triples'] for row in rows)},"
        f"numeric_forward_triples:"
        f"{sum(row['natural_forward_triples'] for row in rows)},"
        f"triangle_checks:{sum(row['triangle_candidates'] for row in rows)}"
    )
    if args.full:
        survivors = tuple(
            (
                row["body"],
                row["rank"],
                row["apex"],
                row["empty_residuals"],
                row["witnesses"],
            )
            for row in rows
            if not row["closed"]
        )
        print(
            f"closure=closed:{sum(row['closed'] for row in rows)},"
            f"survivors:{len(survivors)},"
            f"K4_free:{sum(not row['cliques'] for row in rows)}"
        )
        print(
            f"singleton=scans:{sum(row['singleton_scans'] for row in rows)},"
            f"controls:{sum(row['singleton_controls'] for row in rows)},"
            f"empty:{sum(len(row['empty_residuals']) for row in rows)},"
            f"witnesses:{sum(len(row['witnesses']) for row in rows)}"
        )
        print(f"survivors={survivors}")
    print(f"canonical_ledger_sha256={digest.hexdigest()}")
    print("orientation=degeneracy enumeration gauge on symmetric heavy graph;not tournament")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
