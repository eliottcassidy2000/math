#!/usr/bin/env python3
r"""Refine the coarse heavy-pair graph on the five j=5 pair hostiles.

For a first-apex carrier ``C`` of mass ``m``, let ``B<5m/7`` be the exact
nonexceptional global pair cap from the pair-hostile audit and put
``theta=m-B>2m/7``.  Every four-comb cover avoiding the unique exceptional
pair would make all six of its pairs ``theta``-heavy.  The finite set

    H = {x : |C intersect D_x| >= theta/2}

is a vertex cover of the heavy graph, so every heavy ``K4`` contains a
triangle in ``H``.

The coarse graph is known to contain many false ``K4`` witnesses because its
edge weights double-count.  This scout keeps an edge ``xy`` only when an exact
global pair cap on the literal residual

    C \ (D_x union D_y)

fails to prove that two further combs leave positive measure.  Thus any actual
four-comb cover still induces a clique in the refined graph.  If the refined
graph on ``H`` is triangle-free, the entire nonexceptional branch closes.

This is a targeted five-carrier scout, not a theorem companion for all 30,030
first-apex carriers.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = (
    ROOT
    / "04-computation/lrc14_j5_pair_hostiles_recursive_audit_codex_20260729.py"
)
AUDIT_SHA256 = (
    "aa1880b84a3d54693dffd29d2a687c622d360aee3a975eef7757ada2828389bf"
)
HOSTILES = (
    ((1, 2, 3, 7, 8, 10, 11, 13), 5, 19, (17, 18)),
    ((1, 2, 3, 7, 8, 10, 11, 13), 7, 17, (18, 19)),
    ((1, 4, 6, 7, 9, 10, 11, 13), 6, 23, (16, 17)),
    ((2, 3, 4, 8, 10, 12, 13, 14), 6, 16, (18, 22)),
    ((2, 4, 6, 8, 10, 12, 13, 14), 4, 16, (18, 22)),
)
EXPECTED_ROWS = (
    (
        (1, 2, 3, 7, 8, 10, 11, 13),
        19,
        1031,
        24,
        184,
        F(107083894, 9534924945),
        (17, 84),
        (18, 37),
    ),
    (
        (1, 2, 3, 7, 8, 10, 11, 13),
        17,
        782,
        16,
        87,
        F(346626179567, 43981430462970),
        (19, 37),
        (18, 125),
    ),
    (
        (1, 4, 6, 7, 9, 10, 11, 13),
        23,
        2045,
        47,
        709,
        F(13681079513, 2840499545520),
        (16, 19),
        (17, 108),
    ),
    (
        (2, 3, 4, 8, 10, 12, 13, 14),
        16,
        779,
        15,
        61,
        F(26683, 1261260),
        (17, 18),
        (22, 26),
    ),
    (
        (2, 4, 6, 8, 10, 12, 13, 14),
        16,
        716,
        13,
        48,
        F(356, 35035),
        (18, 20),
        (22, 26),
    ),
)
EXPECTED_DIGEST = (
    "1d682c78bb1b2b9a9b907b3d25d62b15f86f85044334d0a19ea3b7839b6426f5"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_audit():
    require(
        hashlib.sha256(AUDIT_PATH.read_bytes()).hexdigest() == AUDIT_SHA256,
        "pair-hostile audit changed",
    )
    spec = importlib.util.spec_from_file_location("j5_refined_k4_audit", AUDIT_PATH)
    require(spec is not None and spec.loader is not None, "cannot load audit")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = load_audit()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def hostile_profile(
    hostile: tuple[tuple[int, ...], int, int, tuple[int, int]],
) -> dict[str, object]:
    body, apex_rank, apex, exceptional = hostile
    pair_row = A.pair_threshold_profile(
        {"body": body, "prefix": (apex,), "apex_rank": apex_rank}
    )
    require(
        len(pair_row["crossing"]) == 1
        and pair_row["crossing"][0][0] == exceptional,
        f"exceptional pair changed: {body}, {apex}",
    )
    carrier = A.reconstruct_carrier(body, (apex,))
    mass = pair_row["m"]
    theta = mass - pair_row["nonexception_cap"]
    require(theta > 2 * mass / 7, f"heavy threshold not finite: {body}, {apex}")
    level = theta / 2
    threshold = A.S2 * pair_row["r"] / (7 * (level - mass / 7))
    tail_first = A.ceiling(threshold)
    speeds = [
        speed
        for speed in range(A.FIRST_EXTERNAL, tail_first)
        if speed != apex
    ]
    rows = A.P1.THM2885.coverages_many(carrier, speeds)
    head = tuple(sorted(speed for value, speed in rows if value >= level))
    coarse_edges = []
    coarse_rows = []
    for x, y in combinations(head, 2):
        if (x, y) == exceptional:
            continue
        residual = A.P1.THM2883.subtract_local_multi(carrier, (x, y))
        union = mass - sum((right - left for left, right in residual), F(0))
        if union >= theta:
            coarse_edges.append((x, y))
            coarse_rows.append(
                (
                    body,
                    apex_rank,
                    apex,
                    exceptional,
                    theta,
                    x,
                    y,
                    union,
                )
            )
    return {
        "body": body,
        "apex_rank": apex_rank,
        "apex": apex,
        "exceptional": exceptional,
        "m": mass,
        "theta": theta,
        "head_threshold": threshold,
        "head_tail_first": tail_first,
        "head": head,
        "coarse_edges": tuple(coarse_edges),
        "coarse_rows": tuple(coarse_rows),
    }


def refine_edge(
    task: tuple[
        tuple[int, ...],
        int,
        int,
        tuple[int, int],
        F,
        int,
        int,
        F,
    ],
) -> dict[str, object]:
    body, apex_rank, apex, exceptional, theta, x, y, union = task
    first_carrier = A.reconstruct_carrier(body, (apex,))
    residual = A.P1.THM2883.subtract_local_multi(first_carrier, (x, y))
    residual_mass = sum((right - left for left, right in residual), F(0))
    require(residual_mass > 0, f"coarse edge already covers carrier: {body}, {apex}, {x}, {y}")
    cap = A.global_pair_cap(body, (apex, x, y), residual)
    margin = residual_mass - cap["cap"]
    return {
        "body": body,
        "apex_rank": apex_rank,
        "apex": apex,
        "exceptional": exceptional,
        "theta": theta,
        "edge": (x, y),
        "union": union,
        "residual_mass": residual_mass,
        "residual_cap": cap["cap"],
        "residual_margin": margin,
        "residual_maximizer": cap["finite_pair"],
        "residual_paid": cap["paid"],
        "survives": margin <= 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    profiles = [hostile_profile(hostile) for hostile in HOSTILES]
    tasks = [task for profile in profiles for task in profile["coarse_rows"]]
    context = mp.get_context("spawn")
    if args.workers == 1:
        refined = list(map(refine_edge, tasks))
    else:
        with context.Pool(args.workers) as pool:
            refined = pool.map(refine_edge, tasks, chunksize=1)

    by_root: dict[tuple[tuple[int, ...], int], list[dict[str, object]]] = {}
    for row in refined:
        by_root.setdefault((row["body"], row["apex"]), []).append(row)
    ledger_rows = []
    print("LRC14 J=5 FIVE-HOSTILE REFINED K4 SCOUT")
    print("status=FINITE-EXACT-SCOUT;refined_graph_is_necessary_not_equivalent")
    total_coarse = 0
    total_refined = 0
    total_triangles = 0
    total_k4 = 0
    observed_rows = []
    for profile in profiles:
        key = (profile["body"], profile["apex"])
        rows = by_root[key]
        require(
            tuple(row["edge"] for row in rows) == profile["coarse_edges"],
            f"parallel edge order changed: {key}",
        )
        refined_edges = {
            row["edge"]
            for row in rows
            if row["survives"]
        }
        head = profile["head"]

        def is_edge(x: int, y: int) -> bool:
            return tuple(sorted((x, y))) in refined_edges

        triangles = tuple(
            triple
            for triple in combinations(head, 3)
            if all(is_edge(x, y) for x, y in combinations(triple, 2))
        )
        k4s = tuple(
            quad
            for quad in combinations(head, 4)
            if all(is_edge(x, y) for x, y in combinations(quad, 2))
        )
        total_coarse += len(rows)
        total_refined += len(refined_edges)
        total_triangles += len(triangles)
        total_k4 += len(k4s)
        minimum_removed = min(
            (
                row["residual_margin"],
                row["edge"],
                row["residual_maximizer"],
            )
            for row in rows
            if not row["survives"]
        )
        observed_rows.append(
            (
                profile["body"],
                profile["apex"],
                profile["head_tail_first"],
                len(head),
                len(rows),
                minimum_removed[0],
                minimum_removed[1],
                minimum_removed[2],
            )
        )
        ledger_rows.extend(
            "E="
            + ",".join(map(str, row["body"]))
            + f";rank={row['apex_rank']};a={row['apex']};"
            + f"edge={row['edge'][0]},{row['edge'][1]};"
            + f"U={ftext(row['union'])};"
            + f"L={ftext(row['residual_mass'])};"
            + f"B2={ftext(row['residual_cap'])};"
            + f"margin={ftext(row['residual_margin'])};"
            + f"maxpair={row['residual_maximizer']};"
            + f"survives={int(row['survives'])}\n"
            for row in rows
        )
        print(
            f"hostile=E{profile['body']};rank={profile['apex_rank']};"
            f"apex={profile['apex']};exceptional={profile['exceptional']};"
            f"theta={ftext(profile['theta'])};"
            f"head_tail={profile['head_tail_first']};"
            f"head={len(head)};coarse_edges={len(rows)};"
            f"refined_edges={len(refined_edges)};"
            f"triangles={len(triangles)};head_K4={len(k4s)};"
            f"first_triangle={triangles[0] if triangles else None};"
            f"first_K4={k4s[0] if k4s else None};"
            f"minimum_removed_margin={ftext(minimum_removed[0])};"
            f"minimum_removed_edge={minimum_removed[1]};"
            f"minimum_removed_maxpair={minimum_removed[2]}"
        )
    digest = hashlib.sha256(
        (
            "LRC14/j5/five-hostile-refined-residual-pair-graph/v1\n"
            + "".join(ledger_rows)
        ).encode()
    ).hexdigest()
    require(
        tuple(observed_rows) == EXPECTED_ROWS,
        "five-hostile refined row census changed",
    )
    require(
        (total_coarse, total_refined, total_triangles, total_k4)
        == (1089, 0, 0, 0)
        and digest == EXPECTED_DIGEST,
        "refined residual-edge closure changed",
    )
    print(
        f"coarse_edges={total_coarse};refined_edges={total_refined};"
        f"triangles={total_triangles};head_K4={total_k4};"
        f"refined_digest={digest}"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
