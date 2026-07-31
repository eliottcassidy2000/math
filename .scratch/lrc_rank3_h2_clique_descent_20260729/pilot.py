#!/usr/bin/env python3
"""Exact K>=20 pilot for the ranked r=3 H2-heavy K4 descent.

The input universe is the locked scalar-hard ledger behind THM-2896.  On a
branch whose strongest eligible THM-2897 ranked-complement flag is r=3, put

    R3 = q1 + q2 + q3,     theta = h - R3.

THM-2893 then forces at least four labels of a hypothetical five-cover into

    H2 = {v : c_C(v) >= theta/2},

and every pair among those four labels has exact union coverage at least
theta.  Hence it suffices to enumerate the K4s of the exact theta-heavy graph
on the globally sealed H2 core and rule out one-label covers of each literal
four-label residual.

This file is deliberately scratch-only until its universe, exact controls,
ordinary/optimized replay, and workload have been audited.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ENGINE_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.py"
)
ENGINE_SHA256 = (
    "e0ff69252870f194549bba61289c1c5b15bef451e37a72836d71f9e71b1016e9"
)
HARD_LEDGER_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j6_all_root_ranked_suffix_scalar_hard_ledger_codex_20260729.out"
)
HARD_LEDGER_SHA256 = (
    "6be9a6c9218f3b42b2eea733c9050f5d35160664af0f19390337b3c5be57cb37"
)
CHUNK = 10_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_engine():
    require(file_sha256(ENGINE_PATH) == ENGINE_SHA256, "scalar engine changed")
    spec = importlib.util.spec_from_file_location("rank3_h2_scalar_engine", ENGINE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load scalar engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_engine()
require(file_sha256(HARD_LEDGER_PATH) == HARD_LEDGER_SHA256, "hard ledger changed")


def ftext(value: F | None) -> str:
    if value is None:
        return "-"
    return f"{value.numerator}/{value.denominator}"


def parse_fraction_speed(text: str) -> tuple[F, int]:
    speed, value = text.split(":", 1)
    return F(value), int(speed)


def parse_hard_row(line: str) -> dict[str, object]:
    require(line.startswith("HARD;"), "malformed hard-ledger row")
    fields = dict(part.split("=", 1) for part in line.rstrip().split(";")[1:])
    row = {
        "stratum": fields["S"],
        "body": tuple(map(int, fields["E"].split(","))),
        "K": int(fields["K"]),
        "rank": int(fields["rank"]),
        "apex": int(fields["apex"]),
        "prefix": tuple(map(int, fields["prefix"].split(","))),
        "h": F(fields["m"]),
        "components": int(fields["r"]),
        "scalar_margin": F(fields["margin"]),
        "tail_first": int(fields["tail_first"]),
        "top5": tuple(
            parse_fraction_speed(entry) for entry in fields["top5"].split(",")
        ),
    }
    require(row["scalar_margin"] <= 0, "hard ledger contains scalar closure")
    require(len(row["body"]) == 7 and len(row["top5"]) == 5, "row arity changed")
    require(row["prefix"][-1] == row["apex"], "apex/prefix mismatch")
    return row


def ranked_margins(row: dict[str, object]) -> tuple[F, ...]:
    h = row["h"]
    partial = F(0)
    margins = []
    for rank, (value, _) in enumerate(row["top5"][:4], start=1):
        partial += value
        margins.append(F(rank + 2, 7) * h - partial)
    require(
        all(margins[index + 1] < margins[index] for index in range(3)),
        "ranked margins lost strict nesting",
    )
    return tuple(margins)


def strongest_rank(row: dict[str, object]) -> int:
    return max(
        (
            rank
            for rank, margin in enumerate(ranked_margins(row), start=1)
            if margin > 0
        ),
        default=0,
    )


ALL_HIGH_K_HARD_ROWS = tuple(
    row
    for row in map(parse_hard_row, HARD_LEDGER_PATH.read_text().splitlines())
    if row["K"] >= 20
)
TARGET_ROWS = tuple(row for row in ALL_HIGH_K_HARD_ROWS if strongest_rank(row) == 3)
require(len(ALL_HIGH_K_HARD_ROWS) == 449, "K>=20 scalar-hard universe changed")
require(len(TARGET_ROWS) == 143, "K>=20 strongest-r3 universe changed")


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def interval_is_tooth_covered(left: F, right: F, speed: int) -> bool:
    """Whether one lifted interval lies in one radius-1/(14w) tooth."""

    lower_center = ceil_fraction(speed * right - F(1, 14))
    upper_center = floor_fraction(speed * left + F(1, 14))
    return lower_center <= upper_center


def singleton_covers(carrier: list[tuple[F, F]], speed: int) -> bool:
    """Exact early-exit predicate for carrier subset D_speed a.e."""

    return all(
        interval_is_tooth_covered(left, right, speed)
        for left, right in sorted(
            carrier,
            key=lambda interval: interval[1] - interval[0],
            reverse=True,
        )
    )


def sealed_coverages(
    carrier: list[tuple[F, F]],
    excluded: set[int],
    cutoff: int,
) -> tuple[tuple[F, int], ...]:
    """Evaluate every allowed H2 candidate in bounded vector chunks."""

    rows: list[tuple[F, int]] = []
    for first in range(S.FIRST_EXTERNAL, cutoff + 1, CHUNK):
        last = min(cutoff + 1, first + CHUNK)
        speeds = [speed for speed in range(first, last) if speed not in excluded]
        rows.extend(S.G.T.coverages_many(carrier, speeds))
    return tuple(rows)


def degeneracy_order(
    vertices: tuple[int, ...],
    edges: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, ...], int, dict[int, set[int]]]:
    """Return deterministic minimum-current-degree order and symmetric graph."""

    adjacency = {vertex: set() for vertex in vertices}
    for x, y in edges:
        require(x < y and x in adjacency and y in adjacency, "bad heavy edge")
        adjacency[x].add(y)
        adjacency[y].add(x)
    remaining = set(vertices)
    order = []
    degeneracy = 0
    while remaining:
        vertex = min(
            remaining,
            key=lambda item: (len(adjacency[item].intersection(remaining)), item),
        )
        current_degree = len(adjacency[vertex].intersection(remaining))
        degeneracy = max(degeneracy, current_degree)
        order.append(vertex)
        remaining.remove(vertex)
    return tuple(order), degeneracy, adjacency


def four_cliques(
    vertices: tuple[int, ...],
    edges: tuple[tuple[int, int], ...],
    brute_control: bool,
) -> dict[str, object]:
    """Enumerate each K4 as a forward-neighborhood triangle.

    The forward direction is induced solely by a deterministic degeneracy
    order.  It is an enumeration gauge on the symmetric heavy graph, not a
    theorem-bearing orientation and not a tournament.
    """

    order, degeneracy, adjacency = degeneracy_order(vertices, edges)
    position = {vertex: index for index, vertex in enumerate(order)}
    forward = {
        vertex: {
            neighbor
            for neighbor in adjacency[vertex]
            if position[neighbor] > position[vertex]
        }
        for vertex in vertices
    }
    forward_triples = sum(math.comb(len(forward[vertex]), 3) for vertex in vertices)
    natural_forward_triples = sum(
        math.comb(sum(neighbor > vertex for neighbor in adjacency[vertex]), 3)
        for vertex in vertices
    )
    edge_set = set(edges)
    pivoted = []
    triangle_candidates = 0
    for pivot in order:
        neighbors = sorted(forward[pivot], key=position.__getitem__)
        for first_index, first in enumerate(neighbors[:-2]):
            for second_index in range(first_index + 1, len(neighbors) - 1):
                second = neighbors[second_index]
                if tuple(sorted((first, second))) not in edge_set:
                    continue
                for third in neighbors[second_index + 1 :]:
                    triangle_candidates += 1
                    if (
                        tuple(sorted((first, third))) in edge_set
                        and tuple(sorted((second, third))) in edge_set
                    ):
                        pivoted.append(
                            (pivot, tuple(sorted((pivot, first, second, third))))
                        )
    cliques = tuple(clique for _, clique in pivoted)
    require(len(cliques) == len(set(cliques)), "degeneracy K4 enumeration duplicated")
    if brute_control:
        brute = tuple(
            quad
            for quad in combinations(vertices, 4)
            if all(tuple(sorted(pair)) in edge_set for pair in combinations(quad, 2))
        )
        require(
            set(cliques) == set(brute) and len(cliques) == len(brute),
            "degeneracy/brute K4 census differs",
        )
    return {
        "cliques": cliques,
        "pivoted": tuple(pivoted),
        "order": order,
        "degeneracy": degeneracy,
        "forward_degrees": tuple(len(forward[vertex]) for vertex in order),
        "forward_triples": forward_triples,
        "natural_forward_triples": natural_forward_triples,
        "triangle_candidates": triangle_candidates,
    }


def singleton_residual_profile(
    row: dict[str, object],
    clique: tuple[int, int, int, int],
    residual: list[tuple[F, F]],
) -> dict[str, object]:
    """Use the exact longest-component horizon and containment test."""

    require(residual, "singleton profile received empty residual")
    residual_h = mass(residual)
    longest = max(right - left for left, right in residual)
    geometric_last = floor_fraction(F(1, 7) / longest)
    tail_first = geometric_last + 1
    excluded = set((*row["prefix"], *clique))
    scanned = tuple(
        speed
        for speed in range(S.FIRST_EXTERNAL, tail_first)
        if speed not in excluded
    )
    covering = tuple(speed for speed in scanned if singleton_covers(residual, speed))

    controls = []
    if scanned:
        controls.extend((scanned[0], scanned[-1]))
    controls.extend(covering)
    for speed in sorted(set(controls)):
        literal_empty = not S.R.subtract_local(residual, speed)
        scalar_full = S.G.T.coverage(residual, speed) == residual_h
        require(
            singleton_covers(residual, speed) == literal_empty == scalar_full,
            f"singleton control mismatch: {row['body']}, {row['rank']}, "
            f"{clique}, {speed}",
        )

    witnesses = []
    for speed in covering:
        family = tuple(sorted((*row["body"], row["apex"], *clique, speed)))
        direct, direct_r, direct_h = S.G.T.CORE.good_norm(family)
        require(
            not direct and direct_r == 0 and direct_h == 0,
            f"direct five-cover witness failed: {family}",
        )
        witnesses.append(speed)
    return {
        "mass": residual_h,
        "components": len(residual),
        "longest": longest,
        "tail_first": tail_first,
        "scanned": len(scanned),
        "controls": len(set(controls)),
        "covering": tuple(witnesses),
    }


def profile_target(
    task: tuple[dict[str, object], bool, bool],
) -> dict[str, object]:
    row, profile_only, core_only = task
    body = row["body"]
    apex = row["apex"]
    excluded = set(row["prefix"])
    root_good, _, _ = S.G.T.CORE.good_norm(body)
    carrier = S.R.subtract_local(root_good, apex)
    direct, direct_r, direct_h = S.G.T.CORE.good_norm(
        tuple(sorted((*body, apex)))
    )
    h = mass(carrier)
    require(
        carrier == direct
        and len(carrier) == direct_r == row["components"]
        and h == direct_h == row["h"] > 0,
        f"literal/direct carrier mismatch: {body}, rank {row['rank']}",
    )

    top5_speeds = [speed for _, speed in row["top5"]]
    top5_control = tuple(S.G.T.coverages_many(carrier, top5_speeds))
    require(
        top5_control == row["top5"],
        f"top-five control mismatch: {body}, rank {row['rank']}",
    )
    margins = ranked_margins(row)
    require(margins[2] > 0 >= margins[3], "row is not strongest r=3")
    r3 = sum((value for value, _ in row["top5"][:3]), F(0))
    theta = h - r3
    level = theta / 2
    epsilon = level - h / 7
    require(epsilon == margins[2] / 2 > 0, "H2 threshold identity failed")
    gamma = S.S2 * row["components"] / 7
    cutoff = S.ceiling(gamma / epsilon) - 1
    tail_first = cutoff + 1
    require(
        h / 7 + gamma / tail_first <= level,
        f"H2 tail did not seal: {body}, rank {row['rank']}",
    )

    coverage_rows = sealed_coverages(carrier, excluded, cutoff)
    core_rows = tuple(
        sorted(
            (
                (value, speed)
                for value, speed in coverage_rows
                if value >= level
            ),
            key=lambda item: item[1],
        )
    )
    core = tuple(speed for _, speed in core_rows)
    c_by_speed = {speed: value for value, speed in core_rows}
    require(len(core) == len(set(core)), "duplicate H2 label")
    require(all(speed not in excluded for speed in core), "prefix entered H2")
    common = {
        **row,
        "margins": margins,
        "r3": r3,
        "theta": theta,
        "level": level,
        "epsilon": epsilon,
        "cutoff": cutoff,
        "core_rows": core_rows,
        "core": core,
    }
    if core_only:
        return common

    paid_pairs = 0
    singleton_pruned_pairs = 0
    equality_edges = []
    weighted_edges = []
    closest_pruned_gap: tuple[F, tuple[int, int]] | None = None
    after_first: dict[int, list[tuple[F, F]]] = {}
    for x, y in combinations(core, 2):
        singleton_sum = c_by_speed[x] + c_by_speed[y]
        if singleton_sum < theta:
            singleton_pruned_pairs += 1
            gap = theta - singleton_sum
            candidate = (gap, (x, y))
            if closest_pruned_gap is None or candidate < closest_pruned_gap:
                closest_pruned_gap = candidate
            continue
        if x not in after_first:
            after_first[x] = S.R.subtract_local(carrier, x)
        sequential = S.R.subtract_local(after_first[x], y)
        simultaneous = S.R.subtract_local_multi(carrier, (x, y))
        require(
            sequential == simultaneous,
            f"pair residual path mismatch: {body}, {row['rank']}, {(x, y)}",
        )
        union = h - mass(sequential)
        paid_pairs += 1
        if union >= theta:
            weighted_edges.append((x, y, union))
            if union == theta:
                equality_edges.append((x, y))

    edges = tuple((x, y) for x, y, _ in weighted_edges)
    clique_graph = four_cliques(core, edges, brute_control=len(core) <= 80)
    cliques = clique_graph["cliques"]
    base = {
        **common,
        "pairs": math.comb(len(core), 2),
        "paid_pairs": paid_pairs,
        "singleton_pruned_pairs": singleton_pruned_pairs,
        "edges": tuple(weighted_edges),
        "equality_edges": tuple(equality_edges),
        "closest_pruned_gap": closest_pruned_gap,
        "cliques": cliques,
        "pivoted_cliques": clique_graph["pivoted"],
        "degeneracy_order": clique_graph["order"],
        "degeneracy": clique_graph["degeneracy"],
        "forward_degrees": clique_graph["forward_degrees"],
        "forward_triples": clique_graph["forward_triples"],
        "natural_forward_triples": clique_graph["natural_forward_triples"],
        "triangle_candidates": clique_graph["triangle_candidates"],
        "clique_profiles": (),
        "empty_residuals": (),
        "witnesses": (),
        "singleton_scans": 0,
        "singleton_controls": 0,
        "closed": None if profile_only else True,
    }
    if profile_only:
        return base

    clique_profiles = []
    empty_residuals = []
    witnesses = []
    singleton_scans = 0
    singleton_controls = 0
    for clique in cliques:
        residual = S.R.subtract_local_multi(carrier, clique)
        sequential = carrier
        for speed in clique:
            sequential = S.R.subtract_local(sequential, speed)
        require(
            residual == sequential,
            f"K4 residual path mismatch: {body}, {row['rank']}, {clique}",
        )
        if not residual:
            empty_residuals.append(clique)
            clique_profiles.append((clique, None))
            continue
        singleton = singleton_residual_profile(row, clique, residual)
        singleton_scans += singleton["scanned"]
        singleton_controls += singleton["controls"]
        if singleton["covering"]:
            witnesses.extend((clique, speed) for speed in singleton["covering"])
        clique_profiles.append((clique, singleton))
    base.update(
        {
            "clique_profiles": tuple(clique_profiles),
            "empty_residuals": tuple(empty_residuals),
            "witnesses": tuple(witnesses),
            "singleton_scans": singleton_scans,
            "singleton_controls": singleton_controls,
            "closed": not empty_residuals and not witnesses,
        }
    )
    return base


def core_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};a={row['apex']};"
        f"P={row['prefix']};h={ftext(row['h'])};r={row['components']};"
        f"R3={ftext(row['r3'])};theta={ftext(row['theta'])};"
        f"level={ftext(row['level'])};eps={ftext(row['epsilon'])};"
        f"N={row['cutoff']};H={len(row['core'])};"
        f"pairs={row['pairs']};paid={row['paid_pairs']};"
        f"pruned={row['singleton_pruned_pairs']};"
        f"edges={len(row['edges'])};ties={len(row['equality_edges'])};"
        f"deg={row['degeneracy']};"
        f"F3={row['forward_triples']};N3={row['natural_forward_triples']};"
        f"tri_checks={row['triangle_candidates']};"
        f"K4={len(row['cliques'])};closed={row['closed']};"
        f"order={row['degeneracy_order']};"
        f"forward_degrees={row['forward_degrees']};"
        f"Hlabels={tuple(speed for _, speed in row['core_rows'])};"
        f"weighted_edges={tuple((x,y,ftext(w)) for x,y,w in row['edges'])};"
        f"tie_edges={row['equality_edges']}\n"
    )


def sealed_core_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};a={row['apex']};"
        f"P={row['prefix']};h={ftext(row['h'])};r={row['components']};"
        f"R3={ftext(row['r3'])};theta={ftext(row['theta'])};"
        f"level={ftext(row['level'])};eps={ftext(row['epsilon'])};"
        f"N={row['cutoff']};H={len(row['core'])};"
        f"Hrows={tuple((speed,ftext(value)) for value,speed in row['core_rows'])}\n"
    )


def clique_line(row: dict[str, object], item: tuple[object, object]) -> str:
    clique, profile = item
    prefix = (
        f"E={row['body']};rank={row['rank']};a={row['apex']};"
        f"P={row['prefix']};K4={clique};"
    )
    if profile is None:
        return prefix + "residual=EMPTY\n"
    return (
        prefix
        + f"m={ftext(profile['mass'])};r={profile['components']};"
        + f"lambda={ftext(profile['longest'])};"
        + f"tail={profile['tail_first']};scan={profile['scanned']};"
        + f"controls={profile['controls']};cover={profile['covering']}\n"
    )


def nearest_quantiles(
    values: list[int],
    percentages: tuple[int, ...] = (0, 25, 50, 75, 90, 95, 99, 100),
) -> tuple[tuple[int, int], ...]:
    if not values:
        return ()
    ordered = sorted(values)
    size = len(ordered)
    return tuple(
        (
            percentage,
            ordered[
                0
                if percentage == 0
                else ceil_fraction(F(percentage * size, 100)) - 1
            ],
        )
        for percentage in percentages
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 1, 4),
    )
    parser.add_argument(
        "--core-only",
        action="store_true",
        help="stop after globally sealing H2, before any pair evaluation",
    )
    parser.add_argument(
        "--profile-only",
        action="store_true",
        help="stop after globally sealed H2 graph and K4 workload census",
    )
    parser.add_argument(
        "--emit-ledger",
        action="store_true",
        help="emit every literal core, weighted heavy edge, and K4 residual",
    )
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(
        not (args.core_only and args.profile_only),
        "choose at most one of --core-only and --profile-only",
    )
    tasks = tuple((row, args.profile_only, args.core_only) for row in TARGET_ROWS)
    if args.workers == 1:
        rows = [profile_target(task) for task in tasks]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = list(pool.imap(profile_target, tasks, chunksize=1))

    require(
        tuple((row["body"], row["rank"]) for row in rows)
        == tuple((row["body"], row["rank"]) for row in TARGET_ROWS),
        "worker order changed",
    )
    if args.core_only:
        if args.emit_ledger:
            for row in rows:
                print("H2;" + sealed_core_line(row).rstrip())
        ledger = hashlib.sha256(b"LRC14/j6/highK-r3-sealed-H2/v1\n")
        for row in rows:
            ledger.update(sealed_core_line(row).encode())
        print("LRC14 j6 K>=20 strongest-r3 globally sealed H2 core census")
        print(
            f"hard_branches={len(ALL_HIGH_K_HARD_ROWS)};"
            f"strongest_r3_branches={len(rows)}"
        )
        print(
            f"cutoff_quantiles={nearest_quantiles([row['cutoff'] for row in rows])};"
            f"maximum_cutoff="
            f"{max((row['cutoff'],row['body'],row['rank']) for row in rows)}"
        )
        print(
            "H2_size_histogram="
            + repr(tuple(sorted(Counter(len(row["core"]) for row in rows).items())))
        )
        print(
            f"H2_size_quantiles="
            f"{nearest_quantiles([len(row['core']) for row in rows])};"
            f"maximum_H2="
            f"{max((len(row['core']),row['body'],row['rank']) for row in rows)};"
            f"total_pairs={sum(math.comb(len(row['core']),2) for row in rows)}"
        )
        print(f"canonical_ledger_sha256={ledger.hexdigest()}")
        print("scope=exact K>=20 strongest-r3 H2 cores;graph not yet evaluated")
        print("all_exact_controls=PASS")
        return

    if args.emit_ledger:
        for row in rows:
            print("CORE;" + core_line(row).rstrip())
            if not args.profile_only:
                for item in row["clique_profiles"]:
                    print("LEAF;" + clique_line(row, item).rstrip())

    per_root_hard: dict[tuple[int, ...], set[int]] = defaultdict(set)
    for row in ALL_HIGH_K_HARD_ROWS:
        per_root_hard[row["body"]].add(row["rank"])
    closed_by_root: dict[tuple[int, ...], set[int]] = defaultdict(set)
    if not args.profile_only:
        for row in rows:
            if row["closed"]:
                closed_by_root[row["body"]].add(row["rank"])
    root_closures = tuple(
        body
        for body in sorted(per_root_hard)
        if per_root_hard[body] <= closed_by_root[body]
    )

    ledger = hashlib.sha256(b"LRC14/j6/highK-r3-H2-K4-singleton/v1\n")
    for row in rows:
        ledger.update(core_line(row).encode())
        if not args.profile_only:
            for item in row["clique_profiles"]:
                ledger.update(clique_line(row, item).encode())

    print("LRC14 j6 K>=20 strongest-r3 H2-heavy K4 pilot")
    print(
        f"mode={'PROFILE' if args.profile_only else 'FULL'};"
        f"hard_branches={len(ALL_HIGH_K_HARD_ROWS)};"
        f"strongest_r3_branches={len(rows)}"
    )
    print(
        "strongest_rank_histogram="
        + repr(
            tuple(
                sorted(Counter(strongest_rank(row) for row in ALL_HIGH_K_HARD_ROWS).items())
            )
        )
    )
    print(
        f"cutoff_quantiles={nearest_quantiles([row['cutoff'] for row in rows])};"
        f"maximum_cutoff={max((row['cutoff'],row['body'],row['rank']) for row in rows)}"
    )
    print(
        "H2_size_histogram="
        + repr(tuple(sorted(Counter(len(row["core"]) for row in rows).items())))
    )
    print(
        f"H2_size_quantiles={nearest_quantiles([len(row['core']) for row in rows])};"
        f"maximum_H2={max((len(row['core']),row['body'],row['rank']) for row in rows)}"
    )
    print(
        f"graph_totals=pairs:{sum(row['pairs'] for row in rows)},"
        f"paid:{sum(row['paid_pairs'] for row in rows)},"
        f"singleton_pruned:{sum(row['singleton_pruned_pairs'] for row in rows)},"
        f"edges:{sum(len(row['edges']) for row in rows)},"
        f"ties:{sum(len(row['equality_edges']) for row in rows)},"
        f"K4:{sum(len(row['cliques']) for row in rows)}"
    )
    print(
        f"edge_quantiles={nearest_quantiles([len(row['edges']) for row in rows])};"
        f"degeneracy_quantiles={nearest_quantiles([row['degeneracy'] for row in rows])};"
        f"K4_quantiles={nearest_quantiles([len(row['cliques']) for row in rows])};"
        f"maximum_K4={max((len(row['cliques']),row['body'],row['rank']) for row in rows)}"
    )
    print(
        f"ordered_pivot_work=degeneracy_forward_triples:"
        f"{sum(row['forward_triples'] for row in rows)},"
        f"numeric_forward_triples:"
        f"{sum(row['natural_forward_triples'] for row in rows)},"
        f"triangle_checks:{sum(row['triangle_candidates'] for row in rows)},"
        f"K4_free:{sum(not row['cliques'] for row in rows)}"
    )
    if not args.profile_only:
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
            f"survivors:{len(survivors)},root_closures:{len(root_closures)}"
        )
        print(
            f"pre_singleton_closure=K4_free:"
            f"{sum(not row['cliques'] for row in rows)},"
            f"singleton_tested:{sum(bool(row['cliques']) for row in rows)}"
        )
        print(
            f"singleton_totals=scans:{sum(row['singleton_scans'] for row in rows)},"
            f"controls:{sum(row['singleton_controls'] for row in rows)},"
            f"empty:{sum(len(row['empty_residuals']) for row in rows)},"
            f"witnesses:{sum(len(row['witnesses']) for row in rows)}"
        )
        print(f"survivors={survivors}")
        print(f"root_closures={root_closures}")
    print(f"canonical_ledger_sha256={ledger.hexdigest()}")
    print(
        "controls=source-hashes;literal-vs-direct;top5-replay;"
        "strict-tail;pair-sequential-vs-simultaneous;"
        "degeneracy-forward-triangle-vs-brute-K4;"
        "K4-sequential-vs-simultaneous;"
        "singleton-containment-vs-literal-vs-scalar;direct-witness"
    )
    print("scope=exact K>=20 strongest-r3 marked suffixes;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
