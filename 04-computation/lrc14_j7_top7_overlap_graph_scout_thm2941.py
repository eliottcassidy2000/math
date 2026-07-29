#!/usr/bin/env python3
"""Restricted-overlap graph scout for the critical j=7 seam.

This is a sidecar experiment, not a closure theorem.  It imports only the
independent reconstruction in ``lrc_j7_independent_audit_20260729.py`` and
retains the actual restricted pair intersections among each carrier's seven
largest singleton coverages.  It compares three compatible pair objects:

* a three-edge matching forest;
* the best six-edge star;
* the maximum-weight spanning tree.

The scalar G7 envelope reuses one global pair-union cap independently for six
leaves.  These objects instead require the pair credits to coexist on one
labelled graph.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
AUDIT_PATH = HERE / "lrc14_j7_critical_scalar_wall_independent_thm2941.py"
EXPECTED_AUDIT_SHA256 = (
    "0fde0cd17e93a0cbc16de69e42f6dde25e79eaab1c256372825553087b098b92"
)
if (
    hashlib.sha256(AUDIT_PATH.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    != EXPECTED_AUDIT_SHA256
):
    raise RuntimeError("independent THM2941 verifier changed")
SPEC = importlib.util.spec_from_file_location("lrc_j7_independent", AUDIT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load independent j7 audit")
A = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(A)


def perfect_matching_sums(
    vertices: tuple[int, ...],
    edge: dict[tuple[int, int], F],
) -> tuple[tuple[F, tuple[tuple[int, int], ...]], ...]:
    if not vertices:
        return ((F(0), ()),)
    first = vertices[0]
    rows: list[tuple[F, tuple[tuple[int, int], ...]]] = []
    for index in range(1, len(vertices)):
        second = vertices[index]
        pair = tuple(sorted((first, second)))
        rest = vertices[1:index] + vertices[index + 1 :]
        for value, matching in perfect_matching_sums(rest, edge):
            rows.append((edge[pair] + value, (pair,) + matching))
    return tuple(rows)


def maximum_three_matching(
    labels: tuple[int, ...],
    edge: dict[tuple[int, int], F],
) -> tuple[F, tuple[tuple[int, int], ...]]:
    rows: list[tuple[F, tuple[tuple[int, int], ...]]] = []
    for omitted in labels:
        vertices = tuple(label for label in labels if label != omitted)
        rows.extend(perfect_matching_sums(vertices, edge))
    return min(
        rows,
        key=lambda row: (-row[0], tuple(sorted(row[1]))),
    )


def maximum_spanning_tree(
    labels: tuple[int, ...],
    edge: dict[tuple[int, int], F],
) -> tuple[F, tuple[tuple[int, int], ...]]:
    parent = {label: label for label in labels}

    def find(label: int) -> int:
        while parent[label] != label:
            parent[label] = parent[parent[label]]
            label = parent[label]
        return label

    tree: list[tuple[int, int]] = []
    total = F(0)
    ordered = sorted(
        edge,
        key=lambda pair: (-edge[pair], pair),
    )
    for pair in ordered:
        first, second = pair
        root_first = find(first)
        root_second = find(second)
        if root_first == root_second:
            continue
        parent[root_second] = root_first
        tree.append(pair)
        total += edge[pair]
        if len(tree) == len(labels) - 1:
            break
    A.require(len(tree) == len(labels) - 1, "tree did not span")
    return total, tuple(tree)


def owner_bound(
    carrier: tuple[tuple[int, int], ...],
) -> tuple[F, int]:
    longest = F(max(right - left for left, right in carrier), A.RULER)
    # A closed positive-length carrier component contained in the open tooth
    # of D_w has length strictly below 1/(7w).
    bound = A.ceil_fraction(1 / (7 * longest)) - 1
    return longest, bound


def literal_uncovered(
    carrier: tuple[tuple[int, int], ...],
    labels: tuple[int, ...],
) -> F:
    """Exact common-ruler complement, independent of Hunter inequalities."""
    ruler = math.lcm(A.RULER, *(14 * label for label in labels))
    scale = ruler // A.RULER
    carrier_scaled = tuple(
        (left * scale, right * scale) for left, right in carrier
    )
    danger = A.merge_intervals(
        [
            interval
            for label in labels
            for interval in A.danger_intervals(label, ruler)
        ]
    )
    covered = A.intersect_intervals(carrier_scaled, danger)
    uncovered_numerator = (
        sum(right - left for left, right in carrier_scaled)
        - sum(right - left for left, right in covered)
    )
    return F(uncovered_numerator, ruler)


def profile(body: tuple[int, ...]) -> dict[str, object]:
    scalar = A.profile_root(body)
    carrier = A.carrier_for(body)
    labels = tuple(label for label, _ in scalar["top7"])
    coverage = dict(scalar["top7"])
    edges = {
        tuple(sorted(pair)): A.restricted_pair_intersection(carrier, *pair)
        for pair in combinations(labels, 2)
    }
    base = scalar["mass"] - sum(coverage.values(), F(0))

    matching_weight, matching = maximum_three_matching(labels, edges)
    star_rows = tuple(
        (
            sum(
                (
                    edges[tuple(sorted((center, leaf)))]
                    for leaf in labels
                    if leaf != center
                ),
                F(0),
            ),
            center,
        )
        for center in labels
    )
    star_weight, star_center = min(
        star_rows,
        key=lambda row: (-row[0], row[1]),
    )
    tree_weight, tree = maximum_spanning_tree(labels, edges)
    longest, critical_owner_bound = owner_bound(carrier)
    true_gap = literal_uncovered(carrier, labels)
    union_mass = scalar["mass"] - true_gap
    overlap_redundancy = sum(coverage.values(), F(0)) - union_mass
    simplex_delta_sum = sum(
        (value - scalar["mass"] / 7 for value in coverage.values()),
        F(0),
    )
    # For a cover the exact multiplicity identity is
    # integral_C(m-1)=sum_i c_i-h=sum_i delta_i.  For this positive control,
    # overlap_redundancy also contains the uncovered correction.
    A.require(
        overlap_redundancy == simplex_delta_sum + true_gap,
        "multiplicity/uncovered identity failed",
    )

    return {
        **scalar,
        "labels": labels,
        "base": base,
        "edges": edges,
        "matching_weight": matching_weight,
        "matching": matching,
        "matching_margin": base + matching_weight,
        "star_weight": star_weight,
        "star_center_label": star_center,
        "star_margin": base + star_weight,
        "tree_weight": tree_weight,
        "tree": tree,
        "tree_margin": base + tree_weight,
        "longest_component": longest,
        "critical_owner_bound": critical_owner_bound,
        "true_gap": true_gap,
        "overlap_redundancy": overlap_redundancy,
        "simplex_delta_sum": simplex_delta_sum,
    }


def sign_histogram(rows: list[dict[str, object]], field: str) -> tuple[int, int, int]:
    return (
        sum(row[field] < 0 for row in rows),
        sum(row[field] == 0 for row in rows),
        sum(row[field] > 0 for row in rows),
    )


def extreme(
    rows: list[dict[str, object]],
    field: str,
    largest: bool = False,
) -> dict[str, object]:
    if largest:
        return min(rows, key=lambda row: (-row[field], row["body"]))
    return min(rows, key=lambda row: (row[field], row["body"]))


def row_text(row: dict[str, object], field: str) -> str:
    value = row[field]
    value_text = A.ftext(value) if isinstance(value, F) else str(value)
    edge_text = tuple(
        (
            pair,
            A.ftext(row["edges"][pair]),
            math.gcd(*pair),
            tuple(label // math.gcd(*pair) for label in pair),
        )
        for pair in row["tree"]
    )
    return (
        f"value={value_text};body={row['body']};labels={row['labels']};"
        f"base={A.ftext(row['base'])};matching={row['matching']};"
        f"star_center={row['star_center_label']};tree={edge_text};"
        f"longest={A.ftext(row['longest_component'])};"
        f"critical_owner_bound={row['critical_owner_bound']}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(8, mp.cpu_count() or 1),
    )
    args = parser.parse_args()
    A.require(args.workers >= 1, "worker count must be positive")
    roots = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        rows = [profile(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = list(pool.imap(profile, roots, chunksize=1))
    rows.sort(key=lambda row: row["body"])

    A.require(len(rows) == 3_003, "root universe changed")
    print("LRC14 j7 restricted-overlap graph scout")
    print(
        "object=complete labelled graph on the actual top-seven singleton "
        "labels; edge=mu(G_E cap D_u cap D_v)"
    )
    for field in ("matching_margin", "star_margin", "tree_margin"):
        print(f"{field}_signs={sign_histogram(rows, field)}")
        print(f"minimum_{field}=({row_text(extreme(rows, field), field)})")
        print(
            f"maximum_{field}=({row_text(extreme(rows, field, True), field)})"
        )

    print(f"true_gap_signs={sign_histogram(rows, 'true_gap')}")
    print(f"minimum_true_gap=({row_text(extreme(rows, 'true_gap'), 'true_gap')})")
    print(
        f"maximum_true_gap=({row_text(extreme(rows, 'true_gap', True), 'true_gap')})"
    )
    owner_hist: dict[int, int] = {}
    for row in rows:
        bound = row["critical_owner_bound"]
        owner_hist[bound] = owner_hist.get(bound, 0) + 1
    min_longest = extreme(rows, "longest_component")
    print(f"critical_owner_bound_histogram={tuple(sorted(owner_hist.items()))}")
    print(
        "minimum_longest_component="
        f"({row_text(min_longest, 'longest_component')})"
    )
    print(
        "critical_partition_consequence=for a pointwise cover at zero excess, "
        "m=1 a.e.; on each actual closed positive-length carrier component K, "
        "every pair intersection is relatively open and a nonempty one has "
        "positive length, so the relatively open owners are disjoint and "
        "clopen in connected K; one owner contains all of K; closed K inside "
        "an open D_w tooth gives ell<1/(7w); all owner bounds are <=6<15"
    )
    print(
        "critical_gram=for f_i=1_(D_wi)|C-1/7, an exact critical partition "
        "has Gram (h/49)(7I-J), rank 6, and sum_i f_i=0; near-critical "
        "deformations are delta_i=c_i-h/7 and "
        "p_ij=mu(C cap D_wi cap D_wj)"
    )
    print(
        "scope=top-seven compatibility scout and exact critical-boundary "
        "obstruction only; arbitrary seven-label packets remain open"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
