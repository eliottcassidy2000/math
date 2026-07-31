#!/usr/bin/env python3
"""Scalable marked-suffix q5 + 2 B2 census behind THM-2896 gates.

The default high-gate slice profiles every root with K >= 20.  Use
``--min-k 0 --pair-mode scalar-open`` for the staged full 41,415-branch
screen: exact pair caps are then paid only where the cheaper top-five
singleton sum has not already closed the marked branch.
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
GATE_PATH = (
    ROOT / "04-computation/lrc14_j6_all_root_adaptive_gate_atlas_codex_20260729.py"
)
GATE_SHA = "fc36f26d4c8da5b005465696b954eec700c080376eef9ee5ba74a7111def99d7"
GATE_OUTPUT = (
    ROOT / "05-knowledge/results/lrc14_j6_all_root_adaptive_gate_atlas_codex_20260729.out"
)
GATE_OUTPUT_SHA = "3081b93a870faacb31d205e43f7ca87872d7a9f196f4774d8740ced6a314d80b"
PAIR_PATH = (
    ROOT / "04-computation/lrc14_j6_h4_pair_residual_exact_kernel_codex_20260729.py"
)
PAIR_SHA = "b82f318bf89ffd3ab4c918c87736461d068e03f25941aa25a0961d0f74b4d70a"
EXPECTED_HIGH_COUNTS = (
    62,
    1_289,
    840,
    1_289,
    888,
    895,
    394,
    0,
    62,
    30_075,
    2_839,
)
EXPECTED_HIGH_K_DISTRIBUTION = ((20, 30), (21, 20), (22, 9), (23, 2), (25, 1))
EXPECTED_HIGH_DIGEST = (
    "20d8c1fb216815ee9ba4138ab12bd7b2869d4d0bef79baa3d64fdde312e309c7"
)
EXPECTED_HIGH_PER_K_DIGESTS: tuple[tuple[int, str], ...] | None = (
    (20, "0819247438299b5046954a822afe79fc554a9986177dd74ec0aa19b52e74849f"),
    (21, "cfc3b71b5b6782cc400349a2564b110de0184d724c907e26559965357a5cf502"),
    (22, "83cfffbf2abe46344af6608369be70d5e03e12920ccabc69d5ae7b6bf2b2fa68"),
    (23, "b989fcea52136e96085985ef4f0e7a17ac09a14aa4857195df8370e4262de5cb"),
    (25, "558f19c39bd19da434006e4e60e8416797dcb331c8a85604a2867340a1f395be"),
)
UNIQUE_MAX_ROOT = (1, 8, 10, 11, 12, 13, 14)
EXPECTED_UNIQUE_MAX_OPEN_PATH = (23, 27, 19, 46, 18, 17)


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


require(file_sha256(GATE_OUTPUT) == GATE_OUTPUT_SHA, "THM-2896 output changed")
A = load(GATE_PATH, GATE_SHA, "j6_allroot_q5_gate")
H = load(PAIR_PATH, PAIR_SHA, "j6_allroot_q5_pair")


def ftext(value: F | None) -> str:
    if value is None:
        return "-"
    return f"{value.numerator}/{value.denominator}"


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def sealed_ranking(
    carrier: list[tuple[F, F]],
    excluded: set[int],
) -> tuple[list[tuple[F, int]], F, int]:
    """Seal the exact fifth order statistic and the singleton tail."""

    h = mass(carrier)
    r = len(carrier)
    speeds = [
        speed
        for speed in range(H.FIRST_EXTERNAL, H.HORIZON + 1)
        if speed not in excluded
    ]
    rows = H.T.coverages_many(carrier, speeds)
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    require(len(ranked) >= 5, "fewer than five allowed labels")
    q5_base = ranked[4][0]
    require(q5_base > h / 7, "finite rank five does not clear limiting density")
    threshold = H.S2 * r / (7 * (q5_base - h / 7))
    tail_first = max(H.HORIZON + 1, ceiling(threshold))
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
    tail_single = h / 7 + H.S2 * r / (7 * tail_first)
    require(tail_single <= q5_base, "rank-five tail did not seal")
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    require(ranked[4][0] == q5_base, "rank five changed after tail extension")
    return ranked, tail_single, tail_first


def pair_cap_from_ranking(
    carrier: list[tuple[F, F]],
    ranked: list[tuple[F, int]],
    tail_single: F,
) -> dict[str, object]:
    """Exact finite pair maximum with a one-finite/one-tail seal."""

    h = mass(carrier)
    q1 = ranked[0][0]
    head = F(0)
    maximizer: tuple[int, int] | None = None
    paid = 0
    paid_hash = hashlib.sha256(b"LRC14/j6/all-root-q5-pair/paid/v1\n")
    for first_index, (first_value, first) in enumerate(ranked[:-1]):
        if first_value + ranked[first_index + 1][0] <= head:
            break
        after_first = H.R.subtract_local(carrier, first)
        for second_value, second in ranked[first_index + 1 :]:
            if first_value + second_value <= head:
                break
            survivor = H.R.subtract_local(after_first, second)
            survivor_mass = mass(survivor)
            union = h - survivor_mass
            paid += 1
            paid_hash.update(
                (
                    f"P={min(first, second)},{max(first, second)};"
                    f"U={ftext(union)};L={ftext(survivor_mass)}\n"
                ).encode()
            )
            if union > head:
                head = union
                maximizer = tuple(sorted((first, second)))
    require(paid > 0 and maximizer is not None, "empty pair maximization")
    tail = q1 + tail_single
    return {
        "head": head,
        "tail": tail,
        "cap": max(head, tail),
        "maximizer": maximizer,
        "paid": paid,
        "paid_digest": paid_hash.hexdigest(),
    }


def branch_profile(
    root: dict[str, object],
    root_good: list[tuple[F, F]],
    rank: int,
    pair_mode: str,
) -> dict[str, object]:
    body = root["body"]
    apex = root["top"][rank - 1][1]
    prefix = tuple(speed for _, speed in root["top"][:rank])
    excluded = set(prefix)
    carrier = H.R.subtract_local(root_good, apex)
    direct, direct_r, direct_m = H.T.CORE.good_norm(tuple(sorted((*body, apex))))
    h = mass(carrier)
    require(
        carrier == direct and len(carrier) == direct_r and h == direct_m > 0,
        f"literal/direct carrier mismatch: {body}, rank {rank}",
    )
    ranked, tail_single, tail_first = sealed_ranking(carrier, excluded)
    top5 = tuple(ranked[:5])
    q5 = top5[4][0]
    scalar_margin = h - sum((value for value, _ in top5), F(0))
    pay_pair = pair_mode == "all" or scalar_margin <= 0
    pair = (
        pair_cap_from_ranking(carrier, ranked, tail_single)
        if pay_pair
        else None
    )
    pair_margin = h - q5 - 2 * pair["cap"] if pair is not None else None
    return {
        "body": body,
        "K": root["adaptive_k"],
        "rank": rank,
        "apex": apex,
        "prefix": prefix,
        "carrier": carrier,
        "ranked": ranked,
        "tail_single": tail_single,
        "h": h,
        "r": direct_r,
        "top5": top5,
        "q5": q5,
        "b2": pair["cap"] if pair is not None else None,
        "tail_first": tail_first,
        "scalar_margin": scalar_margin,
        "pair_margin": pair_margin,
        "pair_head": pair["head"] if pair is not None else None,
        "pair_tail": pair["tail"] if pair is not None else None,
        "pair_maximizer": pair["maximizer"] if pair is not None else None,
        "paid": pair["paid"] if pair is not None else 0,
        "paid_digest": pair["paid_digest"] if pair is not None else "-",
        "union_closed": scalar_margin > 0
        or (pair_margin is not None and pair_margin > 0),
    }


def profile_body_star(task: tuple[tuple[int, ...], int, str]):
    body, min_k, pair_mode = task
    root = A.profile_body(body)
    if root["adaptive_k"] < min_k:
        return None
    root_good, root_r, root_h = A.T.CORE.good_norm(body)
    require(
        root_r == root["r"] and root_h == root["m"],
        f"root carrier changed: {body}",
    )
    return (
        root,
        tuple(
            branch_profile(root, root_good, rank, pair_mode)
            for rank in range(1, root["adaptive_k"] + 1)
        ),
    )


def branch_ledger_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};"
        f"a={row['apex']};P={row['prefix']};h={ftext(row['h'])};r={row['r']};"
        f"q5={ftext(row['q5'])};tail={row['tail_first']};"
        f"S={ftext(row['scalar_margin'])};Q={ftext(row['pair_margin'])};"
        f"Bhead={ftext(row['pair_head'])};Btail={ftext(row['pair_tail'])};"
        f"Bmax={row['pair_maximizer']};paid={row['paid']};"
        f"paid_digest={row['paid_digest']};U={int(row['union_closed'])}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-k", type=int, default=20)
    parser.add_argument(
        "--workers", type=int, default=min(os.cpu_count() or 1, 4)
    )
    parser.add_argument(
        "--pair-mode", choices=("all", "scalar-open"), default="all"
    )
    args = parser.parse_args()
    require(args.min_k >= 0 and args.workers >= 1, "bad census arguments")
    tasks = tuple((body, args.min_k, args.pair_mode) for body in A.BODIES)
    if args.workers == 1:
        raw = [profile_body_star(task) for task in tasks]
    else:
        context = mp.get_context("spawn")
        with context.Pool(args.workers) as pool:
            raw = list(pool.imap(profile_body_star, tasks, chunksize=4))
    selected = [item for item in raw if item is not None]
    roots = [item[0] for item in selected]
    branches = [row for _, local in selected for row in local]
    require(
        len(branches) == sum(root["adaptive_k"] for root in roots),
        "branch count does not equal total gate size",
    )

    root_rows = []
    for root, local_tuple in selected:
        local = list(local_tuple)
        root_rows.append(
            (
                root["body"],
                root["adaptive_k"],
                sum(row["scalar_margin"] > 0 for row in local),
                sum(
                    row["pair_margin"] is not None
                    and row["pair_margin"] > 0
                    for row in local
                ),
                sum(row["union_closed"] for row in local),
            )
        )
    fully_closed = tuple(
        row for row in root_rows if row[4] == row[1]
    )
    open_roots = tuple(row for row in root_rows if row[4] < row[1])
    open_branches = sorted(
        (row for row in branches if not row["union_closed"]),
        key=lambda row: (
            row["pair_margin"]
            if row["pair_margin"] is not None
            else row["scalar_margin"],
            row["body"],
            row["rank"],
        ),
    )
    counts = (
        len(roots),
        len(branches),
        sum(row["scalar_margin"] > 0 for row in branches),
        sum(row["pair_margin"] is not None for row in branches),
        sum(
            row["pair_margin"] is not None and row["pair_margin"] > 0
            for row in branches
        ),
        sum(row["union_closed"] for row in branches),
        len(open_branches),
        len(fully_closed),
        len(open_roots),
        sum(row["paid"] for row in branches),
        max((row["tail_first"] for row in branches), default=0),
    )
    k_distribution = tuple(
        sorted(Counter(root["adaptive_k"] for root in roots).items())
    )
    digest = hashlib.sha256(b"LRC14/j6/all-root-q5-pair-suffix/v1\n")
    per_k_hashes = {
        k: hashlib.sha256(
            f"LRC14/j6/all-root-q5-pair-suffix/K={k}/v1\n".encode()
        )
        for k, _ in k_distribution
    }
    for row in branches:
        encoded = branch_ledger_line(row).encode()
        digest.update(encoded)
        per_k_hashes[row["K"]].update(encoded)
    digest_text = digest.hexdigest()
    per_k_digests = tuple(
        (k, per_k_hashes[k].hexdigest()) for k, _ in k_distribution
    )
    if args.min_k == 20 and args.pair_mode == "all":
        require(counts == EXPECTED_HIGH_COUNTS, "high-K counts changed")
        require(
            k_distribution == EXPECTED_HIGH_K_DISTRIBUTION,
            "high-K distribution changed",
        )
        require(digest_text == EXPECTED_HIGH_DIGEST, "high-K ledger changed")
        if EXPECTED_HIGH_PER_K_DIGESTS is not None:
            require(
                per_k_digests == EXPECTED_HIGH_PER_K_DIGESTS,
                "high-K stratified ledgers changed",
            )
        unique_open = tuple(
            row["apex"]
            for row in branches
            if row["body"] == UNIQUE_MAX_ROOT and not row["union_closed"]
        )
        require(
            unique_open == EXPECTED_UNIQUE_MAX_OPEN_PATH,
            f"unique maximum-gate paid path changed: {unique_open}",
        )

    print("LRC14 j6 all-root marked-suffix q5+2B2 census")
    print(
        f"parameters=min_k:{args.min_k},workers:{args.workers},"
        f"pair_mode:{args.pair_mode}"
    )
    print(f"counts={counts}")
    print(f"K_distribution={k_distribution}")
    print(f"fully_closed_roots={fully_closed}")
    print(f"open_root_count={len(open_roots)}")
    print(f"open_roots={open_roots}")
    print(f"open_branch_count={len(open_branches)}")
    for row in open_branches[:40]:
        print(
            f"OPEN E={row['body']};K={row['K']};rank={row['rank']};"
            f"apex={row['apex']};prefix={row['prefix']};"
            f"scalar={ftext(row['scalar_margin'])};"
            f"pair={ftext(row['pair_margin'])};q5={ftext(row['q5'])};"
            f"Bhead={ftext(row['pair_head'])};"
            f"Btail={ftext(row['pair_tail'])};"
            f"Bmax={row['pair_maximizer']};paid={row['paid']}"
        )
    print(f"canonical_ledger_sha256={digest_text}")
    print(f"per_K_ledger_sha256={per_k_digests}")
    print("scope=THM-2896 marked branches;certificate census only")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
