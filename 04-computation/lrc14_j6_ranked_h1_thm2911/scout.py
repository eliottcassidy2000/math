#!/usr/bin/env python3
"""Exact earliest-label scout for the THM-2897 r=4 ranked H1 core.

The fixed battery is the exact K>=20 slice of the THM-2896 atlas.  For every
scalar-open marked suffix, this scout tests

    R4=q1+q2+q3+q4 < 6h/7.

When it holds, THM-2893 forces all five residual labels into

    H1={w not in P : c_C(w) >= h-R4}.

The strict discrepancy tail makes H1 finite.  Order H1 by decreasing literal
coverage.  Every hypothetical five-set has a unique earliest label x.  The
scout checks, for every x, that the best four labels in its suffix have total
literal coverage on C\D_x strictly below |C\D_x|.  This is a lossless
depth-one certificate.  Any exception falls back to an exact residual
branch-and-bound; any witness is directly reconstructed.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ENGINE = (
    ROOT
    / "04-computation"
    / "lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.py"
)
ENGINE_SHA = "e0ff69252870f194549bba61289c1c5b15bef451e37a72836d71f9e71b1016e9"

# Exact K>=20 root slice from the locked THM-2896 atlas.
HIGH_K_BODIES = (
    (1, 2, 5, 7, 8, 9, 11),
    (1, 2, 5, 7, 9, 11, 12),
    (1, 2, 5, 9, 11, 12, 14),
    (1, 2, 7, 8, 10, 11, 12),
    (1, 2, 7, 9, 10, 11, 12),
    (1, 2, 8, 9, 10, 11, 14),
    (1, 2, 8, 9, 10, 12, 14),
    (1, 2, 8, 10, 11, 12, 14),
    (1, 2, 9, 10, 11, 12, 14),
    (1, 2, 9, 10, 12, 13, 14),
    (1, 2, 10, 11, 12, 13, 14),
    (1, 3, 5, 7, 8, 9, 11),
    (1, 3, 7, 8, 10, 11, 13),
    (1, 3, 7, 9, 10, 11, 12),
    (1, 3, 8, 10, 11, 12, 14),
    (1, 3, 8, 10, 12, 13, 14),
    (1, 3, 9, 10, 11, 12, 14),
    (1, 3, 10, 11, 12, 13, 14),
    (1, 4, 5, 9, 11, 13, 14),
    (1, 4, 6, 8, 10, 11, 14),
    (1, 4, 6, 9, 10, 12, 14),
    (1, 4, 8, 9, 10, 12, 14),
    (1, 4, 9, 10, 11, 13, 14),
    (1, 4, 9, 10, 12, 13, 14),
    (1, 5, 6, 7, 8, 9, 11),
    (1, 5, 6, 7, 8, 11, 13),
    (1, 5, 6, 8, 9, 10, 14),
    (1, 5, 7, 8, 9, 11, 12),
    (1, 5, 7, 8, 9, 11, 13),
    (1, 5, 7, 8, 11, 12, 13),
    (1, 5, 8, 9, 11, 12, 14),
    (1, 5, 8, 9, 11, 13, 14),
    (1, 5, 8, 11, 12, 13, 14),
    (1, 6, 8, 9, 10, 11, 14),
    (1, 6, 8, 9, 10, 12, 14),
    (1, 6, 8, 10, 11, 13, 14),
    (1, 6, 8, 10, 12, 13, 14),
    (1, 6, 9, 10, 11, 12, 14),
    (1, 6, 9, 10, 12, 13, 14),
    (1, 7, 8, 10, 11, 12, 13),
    (1, 7, 9, 10, 11, 12, 13),
    (1, 8, 9, 10, 11, 12, 14),
    (1, 8, 9, 10, 11, 13, 14),
    (1, 8, 9, 10, 12, 13, 14),
    (1, 8, 10, 11, 12, 13, 14),
    (1, 9, 10, 11, 12, 13, 14),
    (2, 4, 6, 10, 11, 13, 14),
    (2, 4, 8, 9, 10, 12, 14),
    (2, 5, 6, 8, 9, 10, 14),
    (2, 5, 6, 8, 9, 12, 14),
    (2, 5, 7, 8, 9, 11, 12),
    (2, 5, 8, 9, 11, 12, 14),
    (2, 5, 8, 9, 11, 13, 14),
    (2, 5, 9, 11, 12, 13, 14),
    (2, 6, 7, 8, 10, 11, 14),
    (2, 6, 8, 9, 10, 12, 14),
    (2, 6, 8, 10, 11, 12, 14),
    (2, 6, 8, 10, 11, 13, 14),
    (2, 8, 9, 10, 11, 12, 14),
    (2, 8, 10, 11, 12, 13, 14),
    (2, 9, 10, 11, 12, 13, 14),
    (3, 6, 8, 9, 10, 11, 14),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_engine():
    require(
        hashlib.sha256(ENGINE.read_bytes()).hexdigest() == ENGINE_SHA,
        "ranked-suffix engine changed",
    )
    spec = importlib.util.spec_from_file_location("ranked_h1_engine", ENGINE)
    require(spec is not None and spec.loader is not None, "cannot load engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_engine()


def ftext(value: F | None) -> str:
    if value is None:
        return "-"
    return f"{value.numerator}/{value.denominator}"


def residual_mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def nearest_quantiles(
    values: list[int],
    percentages: tuple[int, ...] = (0, 25, 50, 75, 90, 95, 99, 100),
) -> tuple[tuple[int, int], ...]:
    if not values:
        return ()
    values = sorted(values)
    size = len(values)
    return tuple(
        (
            pct,
            values[
                0 if pct == 0 else math.ceil(pct * size / 100) - 1
            ],
        )
        for pct in percentages
    )


def exact_five_cover(
    carrier: list[tuple[F, F]],
    core_rows: tuple[tuple[F, int], ...],
) -> dict[str, object]:
    """Search exactly for a cover by at most five H1 labels."""

    rows = tuple(sorted(core_rows, key=lambda item: (-item[0], item[1])))
    speeds = tuple(speed for _, speed in rows)
    coverages = tuple(value for value, _ in rows)
    size = len(rows)
    stats = Counter()
    memo: set[tuple[int, int, tuple[tuple[F, F], ...]]] = set()
    found: tuple[int, ...] | None = None

    def pad(chosen: tuple[int, ...]) -> tuple[int, ...]:
        padded = list(chosen)
        for speed in speeds:
            if speed not in padded:
                padded.append(speed)
                if len(padded) == 5:
                    break
        require(len(padded) == 5, "cover witness cannot be padded")
        return tuple(padded)

    def dfs(
        start: int,
        slots: int,
        residual: list[tuple[F, F]],
        chosen: tuple[int, ...],
    ) -> tuple[int, ...] | None:
        stats["nodes"] += 1
        if not residual:
            stats["early_covers"] += 1
            return pad(chosen)
        if slots == 0 or start == size:
            stats["leaves"] += 1
            return None

        key = (start, slots, tuple(residual))
        if key in memo:
            stats["memo_prunes"] += 1
            return None
        memo.add(key)

        mass = residual_mass(residual)
        upper = sum(
            coverages[start : start + min(slots, size - start)],
            F(0),
        )
        if upper < mass:
            stats["parent_bound_prunes"] += 1
            return None

        children = []
        for index in range(start, size):
            child = S.R.subtract_local(residual, speeds[index])
            gain = mass - residual_mass(child)
            children.append((index, gain, child))
        local_upper = sum(
            sorted(
                (gain for _, gain, _ in children),
                reverse=True,
            )[:slots],
            F(0),
        )
        if local_upper < mass:
            stats["local_bound_prunes"] += 1
            return None

        for index, gain, child in children:
            stats["candidate_edges"] += 1
            if gain == 0:
                stats["zero_gain_skips"] += 1
                continue
            witness = dfs(
                index + 1,
                slots - 1,
                child,
                (*chosen, speeds[index]),
            )
            if witness is not None:
                return witness
        return None

    if size < 5:
        witness = None
        stats["small_core_closures"] = 1
    else:
        witness = dfs(0, 5, carrier, ())

    if witness is not None:
        residual = carrier
        for speed in witness:
            residual = S.R.subtract_local(residual, speed)
        require(not residual, "returned witness does not cover literal carrier")
    return {
        "witness": witness,
        "nodes": stats["nodes"],
        "leaves": stats["leaves"],
        "candidate_edges": stats["candidate_edges"],
        "parent_prunes": stats["parent_bound_prunes"],
        "local_prunes": stats["local_bound_prunes"],
        "memo_prunes": stats["memo_prunes"],
        "zero_gain_skips": stats["zero_gain_skips"],
        "memo_states": len(memo),
    }


def depth_one_certificate(
    carrier: list[tuple[F, F]],
    core_rows: tuple[tuple[F, int], ...],
) -> dict[str, object]:
    """Certify noncoverage by assigning each five-set to its earliest label.

    ``core_rows`` is sorted by decreasing coverage on ``carrier``.  If a
    five-set exists, its earliest row i leaves four distinct choices in the
    strict suffix i+1,... .  First use the inherited carrier coverages as an
    upper bound.  Only when that fails, recompute the exact four largest
    coverages on the literal residual C\D_x.
    """

    rows = tuple(sorted(core_rows, key=lambda item: (-item[0], item[1])))
    size = len(rows)
    if size < 5:
        return {
            "depth1_closed": True,
            "depth1_checks": size,
            "depth1_parent": 0,
            "depth1_local": 0,
            "depth1_short": size,
            "depth1_unresolved": (),
            "depth1_min_margin": None,
            "depth1_min_row": None,
            "depth1_digest": hashlib.sha256(
                f"SMALL:{size}\n".encode()
            ).hexdigest(),
        }

    certifications: list[tuple[int, int, str, F, F, F]] = []
    unresolved: list[tuple[int, int, F, F, F]] = []
    parent_count = 0
    local_count = 0
    short_count = 0
    for index, (_, speed) in enumerate(rows):
        residual = S.R.subtract_local(carrier, speed)
        mass = residual_mass(residual)
        suffix = rows[index + 1 :]
        if len(suffix) < 4:
            # No five-set can have this row as its earliest member.
            short_count += 1
            certifications.append(
                (index, speed, "SHORT", F(0), mass, F(0))
            )
            continue
        if not residual:
            # A single H1 label covers C.  This is necessarily a genuine
            # padded five-cover because size>=5.
            unresolved.append((index, speed, mass, F(0), F(0)))
            continue

        inherited_upper = sum((value for value, _ in suffix[:4]), F(0))
        inherited_margin = mass - inherited_upper
        if inherited_margin > 0:
            parent_count += 1
            certifications.append(
                (
                    index,
                    speed,
                    "PARENT",
                    inherited_margin,
                    mass,
                    inherited_upper,
                )
            )
            continue

        local_rows = []
        for _, other in suffix:
            child = S.R.subtract_local(residual, other)
            gain = mass - residual_mass(child)
            local_rows.append((gain, other))
        local_rows.sort(key=lambda item: (-item[0], item[1]))
        local_upper = sum((value for value, _ in local_rows[:4]), F(0))
        local_margin = mass - local_upper
        if local_margin > 0:
            local_count += 1
            certifications.append(
                (
                    index,
                    speed,
                    "LOCAL",
                    local_margin,
                    mass,
                    local_upper,
                )
            )
        else:
            unresolved.append(
                (index, speed, mass, inherited_upper, local_upper)
            )

    positive = [
        row for row in certifications if row[2] in ("PARENT", "LOCAL")
    ]
    minimum = min(
        positive,
        key=lambda row: (row[3], row[0], row[1]),
        default=None,
    )
    digest = hashlib.sha256(b"LRC14/j6/H1-earliest-label/v1\n")
    for row in certifications:
        digest.update(
            (
                f"{row[0]};{row[1]};{row[2]};"
                f"{ftext(row[3])};{ftext(row[4])};{ftext(row[5])}\n"
            ).encode()
        )
    for row in unresolved:
        digest.update(
            (
                f"OPEN;{row[0]};{row[1]};{ftext(row[2])};"
                f"{ftext(row[3])};{ftext(row[4])}\n"
            ).encode()
        )
    return {
        "depth1_closed": not unresolved,
        "depth1_checks": len(certifications) + len(unresolved),
        "depth1_parent": parent_count,
        "depth1_local": local_count,
        "depth1_short": short_count,
        "depth1_unresolved": tuple(unresolved),
        "depth1_min_margin": None if minimum is None else minimum[3],
        "depth1_min_row": minimum,
        "depth1_digest": digest.hexdigest(),
    }


def profile_root_task(
    task: tuple[tuple[int, ...], int, int, bool],
) -> dict[str, object]:
    body, max_cutoff, max_combinations, require_high_k = task
    root = S.G.profile_body(body)
    if require_high_k:
        require(root["adaptive_k"] >= 20, f"battery root lost K>=20: {body}")
    good, components, mass = S.G.T.CORE.good_norm(body)
    require(
        components == root["r"] and mass == root["m"],
        f"root reconstruction changed: {body}",
    )
    root["good"] = good
    branches = [
        S.profile_branch(root, rank)
        for rank in range(1, root["adaptive_k"] + 1)
    ]
    hard = [row for row in branches if not row["closed"]]
    rows = []
    for row in hard:
        margins = S.ranked_core_margins(row)
        if margins[3] <= 0:
            continue
        r4 = sum((value for value, _ in row["top5"][:4]), F(0))
        level = row["m"] - r4
        epsilon = level - row["m"] / 7
        require(epsilon == margins[3] > 0, "r4 threshold identity failed")
        cutoff = S.ceiling(S.S2 * row["r"] / (7 * epsilon)) - 1
        require(cutoff >= S.FIRST_EXTERNAL, "H1 cutoff below label range")
        record = {
            "body": body,
            "K": root["adaptive_k"],
            "rank": row["rank"],
            "apex": row["apex"],
            "prefix": row["prefix"],
            "h": row["m"],
            "components": row["r"],
            "r4": r4,
            "margin": margins[3],
            "level": level,
            "cutoff": cutoff,
            "core": (),
            "core_rows": (),
            "combinations": None,
            "status": "CUTOFF_SKIP",
            "witness": None,
            "nodes": 0,
            "leaves": 0,
            "candidate_edges": 0,
            "parent_prunes": 0,
            "local_prunes": 0,
            "memo_prunes": 0,
            "zero_gain_skips": 0,
            "memo_states": 0,
            "depth1_closed": False,
            "depth1_checks": 0,
            "depth1_parent": 0,
            "depth1_local": 0,
            "depth1_short": 0,
            "depth1_unresolved": (),
            "depth1_min_margin": None,
            "depth1_min_row": None,
            "depth1_digest": "-",
        }
        if cutoff > max_cutoff:
            rows.append(record)
            continue

        carrier = S.R.subtract_local(good, row["apex"])
        direct, direct_r, direct_h = S.G.T.CORE.good_norm(
            tuple(sorted((*body, row["apex"])))
        )
        require(
            carrier == direct
            and len(carrier) == direct_r == row["r"]
            and residual_mass(carrier) == direct_h == row["m"],
            f"literal H1 carrier changed: {body}, rank {row['rank']}",
        )
        excluded = set(row["prefix"])
        speeds = [
            speed
            for speed in range(S.FIRST_EXTERNAL, cutoff + 1)
            if speed not in excluded
        ]
        coverages = S.G.T.coverages_many(carrier, speeds)
        core_rows = tuple(
            sorted(
                (
                    (value, speed)
                    for value, speed in coverages
                    if value >= level
                ),
                key=lambda item: (-item[0], item[1]),
            )
        )
        require(
            all(speed not in excluded for _, speed in core_rows),
            "excluded suffix label entered H1",
        )
        core = tuple(speed for _, speed in core_rows)
        workload = math.comb(len(core), 5) if len(core) >= 5 else 0
        record.update(
            {
                "core": core,
                "core_rows": core_rows,
                "combinations": workload,
                "status": "DEPTH1_PENDING",
            }
        )
        depth1 = depth_one_certificate(carrier, core_rows)
        record.update(depth1)
        if depth1["depth1_closed"]:
            record["status"] = "DEPTH1_CLOSED"
            rows.append(record)
            continue

        record["status"] = "COMBINATION_SKIP"
        if workload > max_combinations:
            rows.append(record)
            continue

        search = exact_five_cover(carrier, core_rows)
        record.update(search)
        record["status"] = (
            "WITNESS" if search["witness"] is not None else "CLOSED"
        )
        if search["witness"] is not None:
            family = tuple(
                sorted((*body, row["apex"], *search["witness"]))
            )
            direct_residual, _, direct_mass = S.G.T.CORE.good_norm(family)
            require(
                not direct_residual and direct_mass == 0,
                "witness direct reconstruction did not cover the circle",
            )
        rows.append(record)

    closed_r4 = {
        row["rank"]
        for row in rows
        if row["status"] in ("DEPTH1_CLOSED", "CLOSED")
    }
    unresolved = tuple(
        row["rank"]
        for row in hard
        if row["rank"] not in closed_r4
    )
    return {
        "body": body,
        "K": root["adaptive_k"],
        "scalar_open": len(hard),
        "eligible": len(rows),
        "rows": tuple(rows),
        "unresolved": unresolved,
        "root_closed": not unresolved,
    }


def branch_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};"
        f"a={row['apex']};P={row['prefix']};h={ftext(row['h'])};"
        f"r={row['components']};R4={ftext(row['r4'])};"
        f"margin={ftext(row['margin'])};ell={ftext(row['level'])};"
        f"N={row['cutoff']};H={len(row['core'])};"
        f"C5={row['combinations']};status={row['status']};"
        f"d1checks={row['depth1_checks']};"
        f"d1parent={row['depth1_parent']};"
        f"d1local={row['depth1_local']};"
        f"d1short={row['depth1_short']};"
        f"d1min={ftext(row['depth1_min_margin'])};"
        f"d1open={row['depth1_unresolved']};"
        f"d1digest={row['depth1_digest']};"
        f"witness={row['witness']};nodes={row['nodes']};"
        f"leaves={row['leaves']};edges={row['candidate_edges']};"
        f"pprune={row['parent_prunes']};lprune={row['local_prunes']};"
        f"mprune={row['memo_prunes']};zero={row['zero_gain_skips']};"
        f"memo={row['memo_states']};Hlabels={row['core']}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 1, 4),
    )
    parser.add_argument("--max-cutoff", type=int, default=15_000)
    parser.add_argument("--max-combinations", type=int, default=2_000_000)
    parser.add_argument(
        "--scope",
        choices=("high-k", "all"),
        default="high-k",
    )
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    args = parser.parse_args()
    require(
        args.workers >= 1
        and args.max_cutoff >= S.FIRST_EXTERNAL
        and args.max_combinations >= 0,
        "bad scout parameters",
    )
    require(
        args.shard_count >= 1
        and 0 <= args.shard_index < args.shard_count,
        "bad shard parameters",
    )
    require(
        len(HIGH_K_BODIES) == 62
        and len(set(HIGH_K_BODIES)) == 62,
        "K>=20 battery changed",
    )
    all_bodies = (
        HIGH_K_BODIES if args.scope == "high-k" else S.G.BODIES
    )
    bodies = all_bodies[args.shard_index :: args.shard_count]
    require(bodies, "empty shard")
    tasks = tuple(
        (
            body,
            args.max_cutoff,
            args.max_combinations,
            args.scope == "high-k",
        )
        for body in bodies
    )
    if args.workers == 1:
        roots = [profile_root_task(task) for task in tasks]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            roots = list(pool.imap(profile_root_task, tasks, chunksize=1))
    rows = [row for root in roots for row in root["rows"]]
    statuses = Counter(row["status"] for row in rows)
    exact = [
        row
        for row in rows
        if row["status"] in ("DEPTH1_CLOSED", "CLOSED", "WITNESS")
    ]
    closed = [
        row
        for row in rows
        if row["status"] in ("DEPTH1_CLOSED", "CLOSED")
    ]
    witnesses = [row for row in rows if row["status"] == "WITNESS"]
    finite = [row for row in rows if row["combinations"] is not None]
    root_closures = tuple(
        root["body"] for root in roots if root["root_closed"]
    )
    counts = (
        len(roots),
        sum(root["K"] for root in roots),
        sum(root["scalar_open"] for root in roots),
        len(rows),
        len(finite),
        len(exact),
        len(closed),
        len(witnesses),
        len(root_closures),
    )
    digest = hashlib.sha256(b"LRC14/j6/ranked-H1-core-scout/v1\n")
    for row in rows:
        digest.update(branch_line(row).encode())

    print("LRC14 j6 ranked r4/H1 exact scout")
    print(
        f"parameters=workers:{args.workers},max_cutoff:{args.max_cutoff},"
        f"max_combinations:{args.max_combinations},scope:{args.scope},"
        f"shard:{args.shard_index}/{args.shard_count}"
    )
    print(f"counts={counts}")
    print(f"statuses={tuple(sorted(statuses.items()))}")
    print(
        "K_distribution="
        + repr(tuple(sorted(Counter(root["K"] for root in roots).items())))
    )
    print(
        "cutoff_quantiles="
        + repr(nearest_quantiles([row["cutoff"] for row in rows]))
    )
    print(
        "H1_size_histogram="
        + repr(
            tuple(
                sorted(
                    Counter(len(row["core"]) for row in finite).items()
                )
            )
        )
    )
    print(
        "H1_size_quantiles="
        + repr(nearest_quantiles([len(row["core"]) for row in finite]))
    )
    print(
        "C5_quantiles="
        + repr(
            nearest_quantiles(
                [row["combinations"] for row in finite]
            )
        )
    )
    print(
        "node_quantiles="
        + repr(nearest_quantiles([row["nodes"] for row in exact]))
    )
    print(
        "depth1_check_quantiles="
        + repr(
            nearest_quantiles(
                [row["depth1_checks"] for row in exact]
            )
        )
    )
    print(
        "depth1_local_quantiles="
        + repr(
            nearest_quantiles(
                [row["depth1_local"] for row in exact]
            )
        )
    )
    print(
        "depth1_min_margin="
        + ftext(
            min(
                (
                    row["depth1_min_margin"]
                    for row in exact
                    if row["depth1_min_margin"] is not None
                ),
                default=None,
            )
        )
    )
    print(
        "root_closures="
        + repr(root_closures)
    )
    print(
        "unresolved_per_root="
        + repr(
            tuple(
                sorted(
                    Counter(len(root["unresolved"]) for root in roots).items()
                )
            )
        )
    )
    hardest = sorted(
        finite,
        key=lambda row: (
            row["combinations"],
            row["cutoff"],
            row["body"],
            row["rank"],
        ),
        reverse=True,
    )[:20]
    print("HIGHEST_WORKLOAD")
    for row in hardest:
        print(branch_line(row).rstrip())
    print("WITNESSES")
    for row in witnesses:
        print(branch_line(row).rstrip())
    print("DEPTH1_OPEN")
    for row in rows:
        if (
            row["combinations"] is not None
            and not row["depth1_closed"]
        ):
            print(branch_line(row).rstrip())
    print("ROOT_CLOSURES")
    for root in roots:
        if root["root_closed"]:
            print(
                f"E={root['body']};K={root['K']};"
                f"scalar_open={root['scalar_open']};"
                f"eligible={root['eligible']}"
            )
    print(f"ledger_sha256={digest.hexdigest()}")
    print(
        f"scope={args.scope};shard={args.shard_index}/{args.shard_count};"
        "bounded H1 earliest-label scout;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
