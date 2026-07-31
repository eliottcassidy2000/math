#!/usr/bin/env python3
r"""Exact H3 pair-link census on the 61 unresolved one-hard LRC14 roots.

THM-2901 leaves 61 roots having exactly one scalar-hard marked suffix and
failing its direct q5+2 beta2 certificate.  On each such literal carrier C,

    H3 = {w : c_C(w) >= (h-beta2)/3}

is finite and every hypothetical five-cover contains at least three H3
labels.  Enumerate a pair L from the *actual* H3.  If T is the remaining
three labels, every parent triple L union {z}, z in T, has coverage at least
h-beta2 because its complementary pair has coverage at most beta2.  For the
literal child R=C minus D_L this gives the exact link restriction

    c_R(z) >= h_R-beta2.

Thus all of T lies in

    J_L={z allowed outside L : c_R(z)>=h_R-beta2}.

When delta_L=6h_R/7-beta2>0, discrepancy seals J_L at

    N_L=ceil(gamma_R/delta_L)-1,  gamma_R=(99/70)r_R/7.

The program first tests the sum of the three largest J_L singleton
coverages.  Survivors receive the rank-selective restricted certificate

    q_3(J_L)+B_2(J_L)<h_R,

where B_2(J_L) is the exact maximum pair union over distinct vertices of
J_L.  Literal carriers, ordered prefixes, H3 pairs, and link cores are
retained.  This is a scoped finite census, not LRC(14).
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
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAIR_LEDGER = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.ledger.out"
)
PAIR_LEDGER_SHA256 = (
    "5dea0eaa45dd52fbf1bef7cfcc328899a4789bc277b6e1e8ac2f4bdf192b85e4"
)
RESIDUAL_PATH = (
    ROOT
    / "04-computation/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.py"
)
RESIDUAL_SHA256 = (
    "a5f3dcc1a23defea4b3dc067675d83141f1866022d6d01946617a97de69e5b0e"
)
VECTOR_PATH = (
    ROOT
    / "04-computation/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.py"
)
VECTOR_SHA256 = (
    "dff97f67b1104c25589802a6a2f216b6e7bfedd58eebfa1bcce615d59c1e872f"
)

FIRST_EXTERNAL = 15
S2 = F(99, 70)

# Lock after discovery and replay under ordinary and optimized Python.
EXPECTED_H3_COUNTS: tuple[object, ...] | None = (
    61,
    24,
    461,
    1961,
    2,
    23,
    920,
    (
        (5, 10),
        (7, 9),
        (9, 8),
        (6, 6),
        (4, 6),
        (8, 5),
        (10, 4),
        (3, 3),
        (15, 3),
        (12, 2),
        (2, 2),
        (14, 1),
        (23, 1),
        (17, 1),
    ),
)
EXPECTED_LINK_COUNTS: tuple[object, ...] | None = (
    1961,
    1908,
    53,
    1409,
    484,
    14,
    1,
    86,
    60333,
    1260,
    10104,
    ((False, 1273, 42, 1), (True, 688, 11, 0)),
    51,
    10,
    ((False, 37, 31, 6), (True, 24, 20, 4)),
)
EXPECTED_LEDGER_SHA256: str | None = (
    "b27ca638baf6194be3ee436916233389fe5b7ace7f49020bacc1bd1221c92ebb"
)
EXPECTED_TRIANGLE_COUNTS: tuple[object, ...] | None = (
    10,
    54,
    2,
    (
        ((1, 3, 6, 9, 11, 12, 13), (15, 24, 33)),
        ((1, 3, 6, 9, 11, 12, 13), (24, 30, 33)),
    ),
    2,
    2,
    12,
    61,
    0,
    (),
    ((False, 37), (True, 24)),
)
EXPECTED_PROOF_SHA256: str | None = (
    "08dc4a539544a417ff884ef4631b10d6eb14fd63d7b62b06e76d92e3d4d9b162"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(PAIR_LEDGER) == PAIR_LEDGER_SHA256, "pair ledger changed")
require(file_sha256(RESIDUAL_PATH) == RESIDUAL_SHA256, "residual engine changed")
require(file_sha256(VECTOR_PATH) == VECTOR_SHA256, "vector engine changed")
R = load_module("lrc14_link_residual", RESIDUAL_PATH)
V = load_module("lrc14_link_vector", VECTOR_PATH)


def interval_mass(carrier: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in carrier), F(0))


def parse_fraction(text: str) -> F:
    numerator, denominator = text.split("/")
    return F(int(numerator), int(denominator))


def parse_pair_ledger() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for line in PAIR_LEDGER.read_text().splitlines():
        if not line.startswith("PAIR;"):
            continue
        fields = {
            item.split("=", 1)[0]: item.split("=", 1)[1]
            for item in line.split(";")[1:]
            if "=" in item
        }
        h3_cutoff, h3_count, h3_raw = fields["H3"].split(":")
        rows.append(
            {
                "stratum": fields["S"],
                "body": tuple(map(int, fields["E"].split(","))),
                "K": int(fields["K"]),
                "rank": int(fields["rank"]),
                "apex": int(fields["a"]),
                "prefix": tuple(map(int, fields["P"].split(","))),
                "mass": parse_fraction(fields["h"]),
                "components": int(fields["r"]),
                "pair_cap": parse_fraction(fields["B2"]),
                "pair_margin": parse_fraction(fields["mB2"]),
                "direct_margin": parse_fraction(fields["mdirect"]),
                "h3_cutoff": None if h3_cutoff == "None" else int(h3_cutoff),
                "h3_cutoff_count": None if h3_count == "None" else int(h3_count),
                "h3_cutoff_raw": None if h3_raw == "None" else int(h3_raw),
            }
        )
    require(len(rows) == 14_806, "pair-ledger row count changed")
    return rows


def one_hard_inputs() -> list[dict[str, object]]:
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in parse_pair_ledger():
        by_body[row["body"]].append(row)
    one_hard = [rows[0] for rows in by_body.values() if len(rows) == 1]
    inputs = [row for row in one_hard if row["direct_margin"] <= 0]
    require(
        (len(by_body), len(one_hard), len(inputs)) == (3427, 65, 61),
        "one-hard root universe changed",
    )
    require(
        all(row["pair_margin"] > 0 and row["h3_cutoff"] is not None for row in inputs),
        "one-hard input lost H3 eligibility",
    )
    return sorted(inputs, key=lambda row: row["body"])


def exact_coverages(
    carrier: list[tuple[F, F]],
    labels: list[int],
) -> list[tuple[F, int]]:
    rows = V.coverages_many(carrier, labels)
    require(len(rows) == len(labels), "vector coverage length changed")
    if rows:
        controls = tuple(dict.fromkeys((labels[0], labels[-1], labels[len(labels) // 2])))
        by_label = {label: value for value, label in rows}
        for label in controls:
            require(
                by_label[label] == R.coverage(carrier, label),
                f"scalar/vector mismatch at {label}",
            )
    return rows


def branch_h3(row: dict[str, object]) -> dict[str, object]:
    body = row["body"]
    apex = row["apex"]
    carrier, components, mass = R.CORE.good_norm(tuple(sorted((*body, apex))))
    require(
        components == row["components"] and mass == row["mass"] and mass > 0,
        "literal branch reconstruction changed",
    )
    forbidden = frozenset(row["prefix"])
    cutoff = row["h3_cutoff"]
    require(isinstance(cutoff, int), "missing H3 cutoff")
    labels = [
        label
        for label in range(FIRST_EXTERNAL, cutoff + 1)
        if label not in forbidden
    ]
    require(len(labels) == row["h3_cutoff_count"], "H3 cutoff count changed")
    threshold = (mass - row["pair_cap"]) / 3
    beta = threshold - mass / 7
    gamma = S2 * components / 7
    require(beta > 0, "H3 threshold not above limiting density")
    require(
        mass / 7 + gamma / (cutoff + 1) <= threshold,
        "H3 discrepancy cutoff failed",
    )
    coverages = exact_coverages(carrier, labels)
    core = tuple(label for value, label in coverages if value >= threshold)
    require(
        all(label not in forbidden for label in core),
        "excluded prefix entered actual H3",
    )
    return {
        **row,
        "carrier": carrier,
        "forbidden": forbidden,
        "h3_threshold": threshold,
        "H3": core,
        "actual_h3_size": len(core),
        "actual_h3_pairs": comb(len(core), 2),
        "omits6": 6 not in body,
    }


def finite_pair_cap(
    carrier: list[tuple[F, F]],
    ranked: list[tuple[F, int]],
) -> tuple[F, tuple[int, int] | None, int]:
    if len(ranked) < 2:
        return F(0), None, 0
    mass = interval_mass(carrier)
    best = F(-1)
    witness: tuple[int, int] | None = None
    paid = 0
    for first_index, (first_value, first) in enumerate(ranked[:-1]):
        if first_value + ranked[first_index + 1][0] <= best:
            break
        for second_value, second in ranked[first_index + 1 :]:
            if first_value + second_value <= best:
                break
            pair = tuple(sorted((first, second)))
            survivor = R.subtract_local_multi(carrier, pair)
            union = mass - interval_mass(survivor)
            paid += 1
            if (union, tuple(-x for x in pair)) > (
                best,
                tuple(-x for x in witness) if witness is not None else (),
            ):
                best = union
                witness = pair
    require(witness is not None and paid > 0, "finite pair search stayed empty")
    direct = mass - interval_mass(R.subtract_local_multi(carrier, witness))
    require(best == direct, "restricted pair winner failed direct check")
    return best, witness, paid


def globally_sealed_singleton_top(
    carrier: list[tuple[F, F]],
    excluded: set[int],
    initial_tail: int = 1601,
) -> tuple[F, int, int]:
    """Return the attained global singleton cap and its strict tail start."""

    mass = interval_mass(carrier)
    gamma = S2 * len(carrier) / 7
    labels = [
        label
        for label in range(FIRST_EXTERNAL, initial_tail)
        if label not in excluded
    ]
    rows = exact_coverages(carrier, labels)
    require(rows, "empty singleton head")
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    q1, witness = ranked[0]
    require(q1 > mass / 7, "singleton head did not beat limiting density")
    tail_first = max(initial_tail, ceiling(gamma / (q1 - mass / 7)))
    if tail_first > initial_tail:
        extra = [
            label
            for label in range(initial_tail, tail_first)
            if label not in excluded
        ]
        rows.extend(exact_coverages(carrier, extra))
        q1, witness = min(rows, key=lambda item: (-item[0], item[1]))
    require(
        mass / 7 + gamma / tail_first <= q1,
        "singleton discrepancy tail did not seal",
    )
    return q1, witness, tail_first


def bad_triangle_child(
    branch: dict[str, object],
    triangle: tuple[int, int, int],
) -> dict[str, object]:
    """Exclude the literal two-cover behind one legal bad H3 triangle."""

    carrier = branch["carrier"]
    mass = branch["mass"]
    pair_cap = branch["pair_cap"]
    residual = R.subtract_local_multi(carrier, triangle)
    residual_mass = interval_mass(residual)
    require(residual_mass > 0, "bad triangle covered its parent carrier")
    triple_union = mass - residual_mass
    heavy_margin = triple_union - (mass - pair_cap)
    excluded = set(branch["forbidden"]) | set(triangle)
    base = {
        "body": branch["body"],
        "rank": branch["rank"],
        "apex": branch["apex"],
        "triangle": triangle,
        "residual_mass": residual_mass,
        "residual_components": len(residual),
        "triple_union": triple_union,
        "heavy_margin": heavy_margin,
    }
    if heavy_margin < 0:
        return {
            **base,
            "legal": False,
            "q1": None,
            "q1_witness": None,
            "singleton_tail": None,
            "pair_cutoff": None,
            "pair_head": None,
            "pair_witness": None,
            "paid_pairs": 0,
            "tail_cap": None,
            "margin": None,
            "closed": True,
        }

    direct, direct_components, direct_mass = R.CORE.good_norm(
        tuple(sorted((*branch["body"], branch["apex"], *triangle)))
    )
    require(
        residual == direct
        and residual_mass == direct_mass
        and len(residual) == direct_components,
        "bad-triangle literal/direct child mismatch",
    )
    q1, q1_witness, singleton_tail = globally_sealed_singleton_top(
        residual,
        excluded,
    )
    delta = 6 * residual_mass / 7 - q1
    require(delta > 0, "two-cover target has no one-slot discrepancy gap")
    gamma = S2 * len(residual) / 7
    ratio = gamma / delta
    pair_cutoff = ratio.numerator // ratio.denominator + 1
    labels = [
        label
        for label in range(FIRST_EXTERNAL, pair_cutoff)
        if label not in excluded
    ]
    rows = sorted(
        exact_coverages(residual, labels),
        key=lambda item: (-item[0], item[1]),
    )
    pair_head, pair_witness, paid = finite_pair_cap(residual, rows)
    tail_cap = q1 + residual_mass / 7 + gamma / pair_cutoff
    require(tail_cap < residual_mass, "two-cover discrepancy tail not strict")
    pair_upper = max(pair_head, tail_cap)
    margin = residual_mass - pair_upper
    require(margin > 0, "bad-triangle child pair cap did not close")
    return {
        **base,
        "legal": True,
        "q1": q1,
        "q1_witness": q1_witness,
        "singleton_tail": singleton_tail,
        "pair_cutoff": pair_cutoff,
        "pair_head": pair_head,
        "pair_witness": pair_witness,
        "paid_pairs": paid,
        "tail_cap": tail_cap,
        "margin": margin,
        "closed": True,
    }


def link_child(task: tuple[dict[str, object], tuple[int, int]]) -> dict[str, object]:
    branch, pair = task
    carrier = branch["carrier"]
    mass = branch["mass"]
    pair_cap = branch["pair_cap"]
    residual = R.subtract_local_multi(carrier, pair)
    residual_mass = interval_mass(residual)
    require(residual_mass > 0, "H3 pair covered the parent carrier")
    union_pair = mass - residual_mass
    require(union_pair <= pair_cap, "actual H3 pair exceeded global pair cap")
    direct, direct_components, direct_mass = R.CORE.good_norm(
        tuple(sorted((*branch["body"], branch["apex"], *pair)))
    )
    require(
        residual == direct
        and residual_mass == direct_mass
        and len(residual) == direct_components,
        "literal/direct H3-pair child mismatch",
    )
    delta = 6 * residual_mass / 7 - pair_cap
    base = {
        "body": branch["body"],
        "stratum": branch["stratum"],
        "rank": branch["rank"],
        "apex": branch["apex"],
        "prefix": branch["prefix"],
        "pair": pair,
        "omits6": branch["omits6"],
        "parent_mass": mass,
        "parent_pair_cap": pair_cap,
        "pair_union": union_pair,
        "residual_mass": residual_mass,
        "residual_components": len(residual),
        "delta": delta,
    }
    if delta <= 0:
        return {
            **base,
            "finite": False,
            "link_cutoff": None,
            "link_size": None,
            "top3_margin": None,
            "rank_pair_margin": None,
            "pair_witness": None,
            "paid_pairs": 0,
            "route": "delta-failure",
            "J": (),
        }

    gamma = S2 * len(residual) / 7
    cutoff = ceiling(gamma / delta) - 1
    threshold = residual_mass - pair_cap
    require(
        threshold == residual_mass / 7 + delta and threshold > residual_mass / 7,
        "link threshold identity changed",
    )
    forbidden = set(branch["forbidden"]) | set(pair)
    labels = [
        label
        for label in range(FIRST_EXTERNAL, cutoff + 1)
        if label not in forbidden
    ]
    require(
        residual_mass / 7 + gamma / (cutoff + 1) <= threshold,
        "link discrepancy cutoff failed",
    )
    rows = exact_coverages(residual, labels)
    link_rows = sorted(
        (
            (value, label)
            for value, label in rows
            if value >= threshold
        ),
        key=lambda item: (-item[0], item[1]),
    )
    J = tuple(sorted(label for _, label in link_rows))
    require(
        all(label not in forbidden for label in J),
        "excluded label entered link core",
    )

    # Check the theorem-bearing link identity on every retained vertex.
    for value, label in link_rows:
        full_union = mass - interval_mass(
            R.subtract_local_multi(carrier, (*pair, label))
        )
        require(
            full_union == union_pair + value,
            "parent/child link identity failed",
        )

    if len(link_rows) < 3:
        top3_margin = residual_mass
        rank_pair_margin = residual_mass
        route = "cardinality"
        witness = None
        paid = 0
    else:
        top3_sum = sum((value for value, _ in link_rows[:3]), F(0))
        top3_margin = residual_mass - top3_sum
        if top3_margin > 0:
            rank_pair_margin = None
            route = "top3"
            witness = None
            paid = 0
        else:
            pair_bound, witness, paid = finite_pair_cap(residual, link_rows)
            q3 = link_rows[2][0]
            rank_pair_margin = residual_mass - q3 - pair_bound
            route = "q3+B2" if rank_pair_margin > 0 else "open"

    return {
        **base,
        "finite": True,
        "link_cutoff": cutoff,
        "link_size": len(J),
        "top3_margin": top3_margin,
        "rank_pair_margin": rank_pair_margin,
        "pair_witness": witness,
        "paid_pairs": paid,
        "route": route,
        "J": J,
    }


def row_key(row: dict[str, object]) -> tuple[object, ...]:
    return (row["body"], row["rank"], row["apex"], row["pair"])


def link_line(row: dict[str, object]) -> str:
    return (
        f"E={','.join(map(str, row['body']))};S={row['stratum']};"
        f"rank={row['rank']};a={row['apex']};"
        f"P={','.join(map(str, row['prefix']))};"
        f"L={','.join(map(str, row['pair']))};omit6={int(row['omits6'])};"
        f"h={ftext(row['parent_mass'])};B2={ftext(row['parent_pair_cap'])};"
        f"UL={ftext(row['pair_union'])};hR={ftext(row['residual_mass'])};"
        f"rR={row['residual_components']};delta={ftext(row['delta'])};"
        f"N={row['link_cutoff']};J={','.join(map(str, row['J']))};"
        f"top3margin={None if row['top3_margin'] is None else ftext(row['top3_margin'])};"
        f"rankpairmargin={None if row['rank_pair_margin'] is None else ftext(row['rank_pair_margin'])};"
        f"B2w={row['pair_witness']};paid={row['paid_pairs']};"
        f"route={row['route']}\n"
    )


def triangle_line(row: dict[str, object]) -> str:
    return (
        f"E={','.join(map(str, row['body']))};rank={row['rank']};"
        f"a={row['apex']};T={','.join(map(str, row['triangle']))};"
        f"hR={ftext(row['residual_mass'])};rR={row['residual_components']};"
        f"UT={ftext(row['triple_union'])};"
        f"heavy_margin={ftext(row['heavy_margin'])};legal={int(row['legal'])};"
        f"q1={None if row['q1'] is None else ftext(row['q1'])};"
        f"q1w={row['q1_witness']};Stail={row['singleton_tail']};"
        f"W2={row['pair_cutoff']};"
        f"H2={None if row['pair_head'] is None else ftext(row['pair_head'])};"
        f"H2w={row['pair_witness']};paid={row['paid_pairs']};"
        f"tail={None if row['tail_cap'] is None else ftext(row['tail_cap'])};"
        f"margin={None if row['margin'] is None else ftext(row['margin'])};"
        f"closed={int(row['closed'])}\n"
    )


def nearest_rank(values: list[int]) -> tuple[tuple[int, int], ...]:
    require(values, "empty quantile population")
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
    parser.add_argument("--workload-only", action="store_true")
    parser.add_argument("--ledger", type=Path)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")

    inputs = one_hard_inputs()
    context = mp.get_context("spawn")
    if args.workers == 1:
        branches = list(map(branch_h3, inputs))
    else:
        with context.Pool(args.workers) as pool:
            branches = pool.map(branch_h3, inputs)
    branches.sort(key=lambda row: row["body"])
    h3_counts = (
        len(branches),
        sum(row["omits6"] for row in branches),
        sum(row["actual_h3_size"] for row in branches),
        sum(row["actual_h3_pairs"] for row in branches),
        min(row["actual_h3_size"] for row in branches),
        max(row["actual_h3_size"] for row in branches),
        max(row["h3_cutoff"] for row in branches),
        tuple(Counter(row["actual_h3_size"] for row in branches).most_common()),
    )
    if EXPECTED_H3_COUNTS is not None:
        require(h3_counts == EXPECTED_H3_COUNTS, "actual H3 census changed")
    print("LRC14 one-hard H3 actual-core and link-child census")
    print(f"h3_counts={h3_counts}")
    print(
        "h3_size_quantiles="
        f"{nearest_rank([row['actual_h3_size'] for row in branches])}"
    )
    if args.workload_only:
        print("mode=WORKLOAD_ONLY")
        print("scope=61 non-direct one-hard roots;actual H3 membership;not LRC14")
        return

    tasks = [
        (branch, pair)
        for branch in branches
        for pair in combinations(branch["H3"], 2)
    ]
    require(len(tasks) == h3_counts[3], "actual H3 pair task count changed")
    if args.workers == 1:
        links = list(map(link_child, tasks))
    else:
        with context.Pool(args.workers) as pool:
            links = pool.map(link_child, tasks)
    links.sort(key=row_key)
    finite = [row for row in links if row["finite"]]
    routes = Counter(row["route"] for row in links)
    links_by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in links:
        links_by_body[row["body"]].append(row)
    closed_bodies = tuple(
        sorted(
            body
            for body, body_rows in links_by_body.items()
            if all(
                row["route"] not in {"delta-failure", "open"}
                for row in body_rows
            )
        )
    )
    failed_bodies = tuple(sorted(set(links_by_body) - set(closed_bodies)))
    link_counts = (
        len(links),
        len(finite),
        len(links) - len(finite),
        routes["cardinality"],
        routes["top3"],
        routes["q3+B2"],
        routes["open"],
        sum(row["paid_pairs"] for row in links),
        max((row["link_cutoff"] for row in finite), default=None),
        max((row["link_size"] for row in finite), default=None),
        sum(row["link_size"] for row in finite),
        tuple(
            (
                omit,
                sum(row["omits6"] == omit for row in links),
                sum(row["omits6"] == omit and row["route"] == "delta-failure" for row in links),
                sum(row["omits6"] == omit and row["route"] == "open" for row in links),
            )
            for omit in (False, True)
        ),
        len(closed_bodies),
        len(failed_bodies),
        tuple(
            (
                omit,
                sum((6 not in body) == omit for body in links_by_body),
                sum((6 not in body) == omit for body in closed_bodies),
                sum((6 not in body) == omit for body in failed_bodies),
            )
            for omit in (False, True)
        ),
    )
    if EXPECTED_LINK_COUNTS is not None:
        require(link_counts == EXPECTED_LINK_COUNTS, "link-child census changed")

    branch_by_body = {branch["body"]: branch for branch in branches}
    bad_edges_by_body: dict[tuple[int, ...], set[tuple[int, int]]] = defaultdict(set)
    for row in links:
        if row["route"] in {"delta-failure", "open"}:
            bad_edges_by_body[row["body"]].add(tuple(sorted(row["pair"])))
    bad_triangles: list[tuple[tuple[int, ...], tuple[int, int, int]]] = []
    for body, edges in sorted(bad_edges_by_body.items()):
        vertices = sorted(set().union(*(set(edge) for edge in edges)))
        for triangle in combinations(vertices, 3):
            if all(tuple(sorted(edge)) in edges for edge in combinations(triangle, 2)):
                bad_triangles.append((body, triangle))
    triangle_rows = [
        bad_triangle_child(branch_by_body[body], triangle)
        for body, triangle in bad_triangles
    ]
    graph_closed_bodies = tuple(
        sorted(
            body
            for body in links_by_body
            if all(
                row["closed"]
                for row in triangle_rows
                if row["body"] == body
            )
        )
    )
    graph_failed_bodies = tuple(
        sorted(set(links_by_body) - set(graph_closed_bodies))
    )
    triangle_counts = (
        len(bad_edges_by_body),
        sum(len(edges) for edges in bad_edges_by_body.values()),
        len(bad_triangles),
        tuple(bad_triangles),
        sum(row["legal"] for row in triangle_rows),
        sum(row["closed"] for row in triangle_rows),
        sum(row["paid_pairs"] for row in triangle_rows),
        len(graph_closed_bodies),
        len(graph_failed_bodies),
        graph_failed_bodies,
        tuple(
            (
                omit,
                sum((6 not in body) == omit for body in graph_closed_bodies),
            )
            for omit in (False, True)
        ),
    )
    if EXPECTED_TRIANGLE_COUNTS is not None:
        require(
            triangle_counts == EXPECTED_TRIANGLE_COUNTS,
            "bad-triangle census changed",
        )

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/one-hard-H3-link/v1\n")
    for row in links:
        digest.update(link_line(row).encode())
    ledger_sha256 = digest.hexdigest()
    if EXPECTED_LEDGER_SHA256 is not None:
        require(ledger_sha256 == EXPECTED_LEDGER_SHA256, "link ledger changed")
    if args.ledger is not None:
        args.ledger.write_text(
            "LRC14 j6 one-hard actual H3 link-child ledger\n"
            + "".join("LINK;" + link_line(row) for row in links)
            + f"ledger_sha256={ledger_sha256}\n"
            + "scope=61 non-direct one-hard roots;actual H3 pairs and exact link cores;not LRC14\n"
        )

    proof_digest = hashlib.sha256()
    proof_digest.update(b"LRC14/j6/one-hard-H3-link-proof/v1\n")
    for row in links:
        proof_digest.update(link_line(row).encode())
    for row in triangle_rows:
        proof_digest.update(("TRI;" + triangle_line(row)).encode())
    proof_sha256 = proof_digest.hexdigest()
    if EXPECTED_PROOF_SHA256 is not None:
        require(proof_sha256 == EXPECTED_PROOF_SHA256, "proof digest changed")

    failures = [row for row in links if row["route"] in {"delta-failure", "open"}]
    print(f"link_counts={link_counts}")
    print(f"failed_bodies={failed_bodies}")
    print(f"triangle_counts={triangle_counts}")
    print(f"triangle_rows={tuple(triangle_line(row).strip() for row in triangle_rows)}")
    print(
        "link_cutoff_quantiles="
        f"{nearest_rank([row['link_cutoff'] for row in finite])}"
    )
    print(
        "link_size_quantiles="
        f"{nearest_rank([row['link_size'] for row in finite])}"
    )
    print(
        "failure_rows="
        f"{tuple((row['body'], row['rank'], row['apex'], row['pair'], row['route'], ftext(row['delta'])) for row in failures)}"
    )
    print(f"ledger_sha256={ledger_sha256}")
    print(f"proof_sha256={proof_sha256}")
    print(
        "mode=DISCOVERY"
        if any(
            value is None
            for value in (
                EXPECTED_H3_COUNTS,
                EXPECTED_LINK_COUNTS,
                EXPECTED_LEDGER_SHA256,
                EXPECTED_TRIANGLE_COUNTS,
                EXPECTED_PROOF_SHA256,
            )
        )
        else "mode=LOCKED"
    )
    print(
        "scope=61 non-direct one-hard roots;actual H3 pair flags;"
        "delta-sealed residual links and restricted child caps;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
