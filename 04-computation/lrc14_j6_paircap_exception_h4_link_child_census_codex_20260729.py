#!/usr/bin/env python3
r"""Exact heavy-link child census on THM-2901's 52 pair-cap exceptions.

For one of the 52 exception carriers C, let h=|C| and let q1 be its attained
global allowed singleton cap.  Since q1<3h/7, every hypothetical five-cover
contains at least two vertices of

    H4(C)={w:c_C(w)>=(h-q1)/4}.

For every unordered pair L from the actual H4 core, form the literal residual
R=C\D_L.  If T is the remaining three labels of a parent cover, every
two-subset A of T obeys the exact heavy-link condition

    U_R(A)>=theta_L,                 theta_L=|R|-q1.

When theta_L/2>|R|/7, every heavy edge meets the finite core

    H2(R)={w:c_R(w)>=theta_L/2}.

Consequently a heavy triangle T contains at least two H2 vertices.  The script
enumerates those core pairs, retains only actual heavy edges, and uses exact
literal singleton completion (with the longest-component cutoff) for the last
vertex, while checking all three link edges.

If theta_L/2<=|R|/7, or if the analytic H2 head/core is deliberately routed
away from an impractically large finite search, the script uses a generic but
exact cover recursion.  Every three-cover of R has a vertex in

    G3(R)={w:c_R(w)>=|R|/3};

behind such a vertex every residual two-cover has a vertex of coverage at
least half the new residual.  The final label is again decided by literal
singleton completion.  Heavy-link sidecars are checked throughout, so this
fallback loses efficiency, not correctness.

All labels respect the inherited forbidden prefix and every selected flag.
This is a scoped exact child computation, not LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MEMBERSHIP_PATH = (
    ROOT
    / "04-computation/lrc14_j6_paircap_exception_h4_membership_census_codex_20260729.py"
)
MEMBERSHIP_SHA256 = (
    "63a80908a6380a877345f0cc4aba7a5e0ef2bb3d59b1b10d58367444ed406b75"
)

FIRST_EXTERNAL = 15
S2 = F(99, 70)

# Discovery values are locked only after ordinary and optimized replay.
EXPECTED_COUNTS: tuple[object, ...] | None = None
EXPECTED_DIGEST: str | None = None


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(MEMBERSHIP_PATH) == MEMBERSHIP_SHA256, "H4 census changed")
M = load_module("lrc14_h4_membership", MEMBERSHIP_PATH)
R = M.E.R
V = M.E.T


def ftext(value: F | None) -> str:
    if value is None:
        return "None"
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def interval_mass(carrier: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in carrier), F(0))


def exact_coverages(
    carrier: list[tuple[F, F]],
    labels: list[int],
    controls: bool = True,
) -> list[tuple[F, int]]:
    rows = V.coverages_many(carrier, labels)
    require(len(rows) == len(labels), "vector coverage length changed")
    if rows and controls:
        controls = tuple(
            dict.fromkeys((labels[0], labels[-1], labels[len(labels) // 2]))
        )
        by_label = {label: value for value, label in rows}
        for label in controls:
            require(
                by_label[label] == R.coverage(carrier, label),
                f"scalar/vector mismatch at {label}",
            )
    return rows


def finite_core(
    carrier: list[tuple[F, F]],
    threshold: F,
    excluded: frozenset[int],
) -> tuple[list[tuple[F, int]], int]:
    """Return exact core rows and the globally sealed inclusive cutoff."""

    mass = interval_mass(carrier)
    require(mass > 0 and threshold > mass / 7, "non-finite core request")
    gamma = S2 * len(carrier) / 7
    cutoff = ceiling(gamma / (threshold - mass / 7)) - 1
    require(
        mass / 7 + gamma / (cutoff + 1) <= threshold,
        "discrepancy tail did not seal",
    )
    labels = [
        label
        for label in range(FIRST_EXTERNAL, cutoff + 1)
        if label not in excluded
    ]
    rows = [
        (value, label)
        for value, label in exact_coverages(carrier, labels)
        if value >= threshold
    ]
    rows.sort(key=lambda item: (-item[0], item[1]))
    require(
        all(label not in excluded for _, label in rows),
        "excluded label entered finite core",
    )
    return rows, cutoff


def singleton_cover_candidates(
    carrier: list[tuple[F, F]],
    excluded: frozenset[int],
) -> tuple[list[int], int, int]:
    """Return every allowed label whose danger comb covers carrier exactly."""

    if not carrier:
        return [], 0, 0
    mass = interval_mass(carrier)
    longest = max(right - left for left, right in carrier)
    bound_ratio = F(1, 7) / longest
    cutoff = bound_ratio.numerator // bound_ratio.denominator
    labels = [
        label
        for label in range(FIRST_EXTERNAL, cutoff + 1)
        if label not in excluded
    ]
    # The pinned vector primitive is exact.  Candidate equalities receive the
    # stronger independent literal-residual check below; repeating three slow
    # scalar controls at every terminal leaf would dominate this census.
    rows = exact_coverages(carrier, labels, controls=False)
    candidates = [label for value, label in rows if value == mass]
    for label in candidates:
        require(
            not R.subtract_local(carrier, label),
            "coverage equality did not empty literal residual",
        )
    return candidates, cutoff, len(labels)


def union_mass(
    carrier: list[tuple[F, F]],
    labels: tuple[int, ...],
) -> F:
    return interval_mass(carrier) - interval_mass(
        R.subtract_local_multi(carrier, labels)
    )


def link_triangle_ok(
    carrier: list[tuple[F, F]],
    triangle: tuple[int, int, int],
    threshold: F,
) -> bool:
    return all(
        union_mass(carrier, tuple(sorted(pair))) >= threshold
        for pair in combinations(triangle, 2)
    )


def finite_h2_route(
    residual: list[tuple[F, F]],
    threshold: F,
    excluded: frozenset[int],
    h2_rows: list[tuple[F, int]],
) -> dict[str, object]:
    """Decide heavy triangles via two high-core vertices and one completion."""

    pair_unions = 0
    heavy_edges = 0
    singleton_heads = 0
    max_singleton_cutoff = 0
    attempted_core_pairs = 0
    coverage_by_label = {label: value for value, label in h2_rows}
    for first, second in combinations(sorted(coverage_by_label), 2):
        attempted_core_pairs += 1
        if coverage_by_label[first] + coverage_by_label[second] < threshold:
            continue
        edge_union = union_mass(residual, (first, second))
        pair_unions += 1
        if edge_union < threshold:
            continue
        heavy_edges += 1
        leaf = R.subtract_local_multi(residual, (first, second))
        if not leaf:
            return {
                "closed": False,
                "witness": (first, second, None),
                "partial_witness": True,
                "core_pairs": attempted_core_pairs,
                "pair_unions": pair_unions,
                "heavy_edges": heavy_edges,
                "singleton_heads": singleton_heads,
                "max_singleton_cutoff": max_singleton_cutoff,
            }
        leaf_excluded = excluded | frozenset((first, second))
        candidates, cutoff, paid = singleton_cover_candidates(
            leaf,
            leaf_excluded,
        )
        singleton_heads += paid
        max_singleton_cutoff = max(max_singleton_cutoff, cutoff)
        for third in candidates:
            triangle = tuple(sorted((first, second, third)))
            if link_triangle_ok(residual, triangle, threshold):
                return {
                    "closed": False,
                    "witness": triangle,
                    "partial_witness": False,
                    "core_pairs": attempted_core_pairs,
                    "pair_unions": pair_unions,
                    "heavy_edges": heavy_edges,
                    "singleton_heads": singleton_heads,
                    "max_singleton_cutoff": max_singleton_cutoff,
                }
    return {
        "closed": True,
        "witness": None,
        "partial_witness": False,
        "core_pairs": attempted_core_pairs,
        "pair_unions": pair_unions,
        "heavy_edges": heavy_edges,
        "singleton_heads": singleton_heads,
        "max_singleton_cutoff": max_singleton_cutoff,
    }


def generic_three_route(
    residual: list[tuple[F, F]],
    link_threshold: F,
    excluded: frozenset[int],
) -> dict[str, object]:
    """Exact generic three-cover recursion, filtered by all heavy-link edges."""

    mass = interval_mass(residual)
    g3_rows, g3_cutoff = finite_core(residual, mass / 3, excluded)
    pair_unions = 0
    singleton_heads = 0
    max_g2_cutoff = 0
    max_singleton_cutoff = 0
    g2_vertices = 0
    for _, first in g3_rows:
        first_residual = R.subtract_local(residual, first)
        if not first_residual:
            return {
                "closed": False,
                "witness": (first, None, None),
                "partial_witness": True,
                "g3_size": len(g3_rows),
                "g3_cutoff": g3_cutoff,
                "g2_vertices": g2_vertices,
                "max_g2_cutoff": max_g2_cutoff,
                "pair_unions": pair_unions,
                "singleton_heads": singleton_heads,
                "max_singleton_cutoff": max_singleton_cutoff,
            }
        first_excluded = excluded | frozenset((first,))
        first_mass = interval_mass(first_residual)
        g2_rows, g2_cutoff = finite_core(
            first_residual,
            first_mass / 2,
            first_excluded,
        )
        g2_vertices += len(g2_rows)
        max_g2_cutoff = max(max_g2_cutoff, g2_cutoff)
        for _, second in g2_rows:
            if union_mass(residual, tuple(sorted((first, second)))) < link_threshold:
                pair_unions += 1
                continue
            pair_unions += 1
            leaf = R.subtract_local(first_residual, second)
            if not leaf:
                return {
                    "closed": False,
                    "witness": (first, second, None),
                    "partial_witness": True,
                    "g3_size": len(g3_rows),
                    "g3_cutoff": g3_cutoff,
                    "g2_vertices": g2_vertices,
                    "max_g2_cutoff": max_g2_cutoff,
                    "pair_unions": pair_unions,
                    "singleton_heads": singleton_heads,
                    "max_singleton_cutoff": max_singleton_cutoff,
                }
            second_excluded = first_excluded | frozenset((second,))
            candidates, cutoff, paid = singleton_cover_candidates(
                leaf,
                second_excluded,
            )
            singleton_heads += paid
            max_singleton_cutoff = max(max_singleton_cutoff, cutoff)
            for third in candidates:
                triangle = tuple(sorted((first, second, third)))
                if link_triangle_ok(residual, triangle, link_threshold):
                    return {
                        "closed": False,
                        "witness": triangle,
                        "partial_witness": False,
                        "g3_size": len(g3_rows),
                        "g3_cutoff": g3_cutoff,
                        "g2_vertices": g2_vertices,
                        "max_g2_cutoff": max_g2_cutoff,
                        "pair_unions": pair_unions,
                        "singleton_heads": singleton_heads,
                        "max_singleton_cutoff": max_singleton_cutoff,
                    }
    return {
        "closed": True,
        "witness": None,
        "partial_witness": False,
        "g3_size": len(g3_rows),
        "g3_cutoff": g3_cutoff,
        "g2_vertices": g2_vertices,
        "max_g2_cutoff": max_g2_cutoff,
        "pair_unions": pair_unions,
        "singleton_heads": singleton_heads,
        "max_singleton_cutoff": max_singleton_cutoff,
    }


def branch_rows() -> list[dict[str, object]]:
    rows = [M.exact_row(fields) for fields in M.parse_exceptions()]
    require(len(rows) == 52, "exception universe changed")
    for row in rows:
        carrier, components, mass = R.CORE.good_norm(
            tuple(sorted((*row["body"], row["apex"])))
        )
        require(
            components == row["components"] and mass == row["mass"],
            "parent carrier reconstruction changed",
        )
        row["carrier"] = carrier
        row["forbidden"] = frozenset(row["prefix"])
    return rows


def child_task(
    task: tuple[dict[str, object], tuple[int, int], int, int],
) -> dict[str, object]:
    branch, pair, h2_head_limit, h2_core_limit = task
    carrier = branch["carrier"]
    parent_mass = branch["mass"]
    q1 = branch["q1"]
    residual = R.subtract_local_multi(carrier, pair)
    residual_mass = interval_mass(residual)
    require(residual_mass > 0, "H4 pair covered the parent carrier")
    pair_union = parent_mass - residual_mass
    direct, direct_components, direct_mass = R.CORE.good_norm(
        tuple(sorted((*branch["body"], branch["apex"], *pair)))
    )
    require(
        residual == direct
        and residual_mass == direct_mass
        and len(residual) == direct_components,
        "literal/direct H4-pair child mismatch",
    )
    link_threshold = residual_mass - q1
    h2_threshold = link_threshold / 2
    delta = h2_threshold - residual_mass / 7
    excluded = branch["forbidden"] | frozenset(pair)
    gamma = S2 * len(residual) / 7
    analytic_h2_cutoff = (
        ceiling(gamma / delta) - 1 if delta > 0 else None
    )

    if (
        delta > 0
        and analytic_h2_cutoff is not None
        and analytic_h2_cutoff <= h2_head_limit
    ):
        h2_rows, h2_cutoff = finite_core(
            residual,
            h2_threshold,
            excluded,
        )
        require(h2_cutoff == analytic_h2_cutoff, "H2 cutoff changed")
        h2_size = len(h2_rows)
        if h2_size <= h2_core_limit:
            result = finite_h2_route(
                residual,
                link_threshold,
                excluded,
                h2_rows,
            )
            route = "finite-H2"
            generic_reason = None
        else:
            result = generic_three_route(
                residual,
                link_threshold,
                excluded,
            )
            route = "generic-three"
            generic_reason = "large-H2-core"
    else:
        result = generic_three_route(
            residual,
            link_threshold,
            excluded,
        )
        route = "generic-three"
        generic_reason = "delta-nonpositive" if delta <= 0 else "large-H2-head"
        h2_size = None

    return {
        "body": branch["body"],
        "stratum": branch["stratum"],
        "rank": branch["rank"],
        "apex": branch["apex"],
        "prefix": branch["prefix"],
        "pair": pair,
        "parent_mass": parent_mass,
        "q1": q1,
        "pair_union": pair_union,
        "residual_mass": residual_mass,
        "residual_components": len(residual),
        "link_threshold": link_threshold,
        "h2_delta": delta,
        "analytic_h2_cutoff": analytic_h2_cutoff,
        "h2_size": h2_size,
        "route": route,
        "generic_reason": generic_reason,
        **result,
    }


def row_key(row: dict[str, object]) -> tuple[object, ...]:
    return (row["body"], row["rank"], row["apex"], row["pair"])


def row_line(row: dict[str, object]) -> str:
    return (
        f"E={','.join(map(str, row['body']))};S={row['stratum']};"
        f"rank={row['rank']};a={row['apex']};"
        f"P={','.join(map(str, row['prefix']))};"
        f"L={','.join(map(str, row['pair']))};"
        f"h={ftext(row['parent_mass'])};q1={ftext(row['q1'])};"
        f"UL={ftext(row['pair_union'])};hR={ftext(row['residual_mass'])};"
        f"rR={row['residual_components']};theta={ftext(row['link_threshold'])};"
        f"delta2={ftext(row['h2_delta'])};"
        f"N2={row['analytic_h2_cutoff']};H2={row['h2_size']};"
        f"route={row['route']};reason={row['generic_reason']};"
        f"closed={int(row['closed'])};witness={row['witness']};"
        f"partial={int(row['partial_witness'])};"
        f"corepairs={row.get('core_pairs')};"
        f"heavyedges={row.get('heavy_edges')};"
        f"G3={row.get('g3_size')};N3={row.get('g3_cutoff')};"
        f"G2total={row.get('g2_vertices')};N2child={row.get('max_g2_cutoff')};"
        f"pairunions={row['pair_unions']};"
        f"singletonheads={row['singleton_heads']};"
        f"W1={row['max_singleton_cutoff']}\n"
    )


def nearest_rank(values: list[int]) -> tuple[tuple[int, int], ...]:
    if not values:
        return ()
    ordered = sorted(values)
    return tuple(
        (
            p,
            ordered[0 if p == 0 else (p * len(ordered) + 99) // 100 - 1],
        )
        for p in (0, 25, 50, 75, 90, 95, 99, 100)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--h2-head-limit", type=int, default=2000)
    parser.add_argument("--h2-core-limit", type=int, default=12)
    parser.add_argument("--ledger", type=Path)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.h2_head_limit >= FIRST_EXTERNAL, "H2 head limit too small")
    require(args.h2_core_limit >= 2, "H2 core limit too small")

    branches = branch_rows()
    tasks = [
        (branch, pair, args.h2_head_limit, args.h2_core_limit)
        for branch in branches
        for pair in combinations(branch["core"], 2)
    ]
    require(len(tasks) == 18_290, "H4 pair-flag count changed")
    context = mp.get_context("spawn")
    if args.workers == 1:
        rows = list(map(child_task, tasks))
    else:
        with context.Pool(args.workers) as pool:
            rows = pool.map(child_task, tasks, chunksize=8)
    rows.sort(key=row_key)

    routes = Counter(row["route"] for row in rows)
    reasons = Counter(row["generic_reason"] for row in rows if row["generic_reason"])
    open_rows = [row for row in rows if not row["closed"]]
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_body[row["body"]].append(row)
    closed_bodies = tuple(
        sorted(
            body
            for body, body_rows_ in by_body.items()
            if all(row["closed"] for row in body_rows_)
        )
    )
    open_bodies = tuple(sorted(set(by_body) - set(closed_bodies)))
    finite_rows = [row for row in rows if row["route"] == "finite-H2"]
    generic_rows = [row for row in rows if row["route"] == "generic-three"]
    counts = (
        len(rows),
        routes["finite-H2"],
        routes["generic-three"],
        reasons["delta-nonpositive"],
        reasons["large-H2-head"],
        reasons["large-H2-core"],
        sum(row["closed"] for row in rows),
        len(open_rows),
        len(closed_bodies),
        len(open_bodies),
        sum(
            row["h2_size"]
            for row in rows
            if row["h2_size"] is not None
        ),
        max(
            (
                row["h2_size"]
                for row in rows
                if row["h2_size"] is not None
            ),
            default=None,
        ),
        sum(row["pair_unions"] for row in rows),
        sum(row["singleton_heads"] for row in rows),
        max(row["max_singleton_cutoff"] for row in rows),
        sum(row.get("heavy_edges", 0) for row in finite_rows),
        sum(row.get("g3_size", 0) for row in generic_rows),
        sum(row.get("g2_vertices", 0) for row in generic_rows),
        max((row.get("g3_cutoff", 0) for row in generic_rows), default=None),
        max((row.get("max_g2_cutoff", 0) for row in generic_rows), default=None),
    )
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "aggregate counts changed")

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/paircap-exception/H4-link-child/v1\n")
    for row in rows:
        digest.update(row_line(row).encode())
    semantic_digest = digest.hexdigest()
    if EXPECTED_DIGEST is not None:
        require(semantic_digest == EXPECTED_DIGEST, "semantic digest changed")

    if args.ledger is not None:
        args.ledger.write_text(
            "LRC14 j6 pair-cap-exception H4 link-child ledger\n"
            + "".join("CHILD;" + row_line(row) for row in rows)
            + f"counts={counts}\n"
            + f"semantic_digest={semantic_digest}\n"
            + "scope=52 exact-B2 exceptions;18,290 actual H4 pair flags;"
            "exact heavy-link child recursion;not LRC14\n"
        )

    print("LRC14 j6 pair-cap-exception H4 link-child census")
    print(f"h2_head_limit={args.h2_head_limit}")
    print(f"h2_core_limit={args.h2_core_limit}")
    print(f"counts={counts}")
    print(f"routes={tuple(sorted(routes.items()))}")
    print(f"generic_reasons={tuple(sorted(reasons.items()))}")
    print(f"closed_bodies={closed_bodies}")
    print(f"open_bodies={open_bodies}")
    print(
        "h2_cutoff_quantiles="
        f"{nearest_rank([row['analytic_h2_cutoff'] for row in finite_rows])}"
    )
    print(
        "h2_size_quantiles="
        f"{nearest_rank([row['h2_size'] for row in rows if row['h2_size'] is not None])}"
    )
    print(
        "generic_g3_size_quantiles="
        f"{nearest_rank([row['g3_size'] for row in generic_rows])}"
    )
    print(
        "open_rows="
        f"{tuple((row['body'], row['rank'], row['apex'], row['pair'], row['route'], row['witness'], row['partial_witness']) for row in open_rows)}"
    )
    print(f"semantic_digest={semantic_digest}")
    print(
        "mode=DISCOVERY"
        if EXPECTED_COUNTS is None or EXPECTED_DIGEST is None
        else "mode=LOCKED"
    )
    print(
        "scope=52 exact-B2 exceptions;18,290 actual H4 pair flags;"
        "finite-H2 and generic exact heavy-link recursion;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
