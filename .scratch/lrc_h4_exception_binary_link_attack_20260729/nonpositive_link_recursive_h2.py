#!/usr/bin/env python3
"""Close the nonfinite binary-link flags by an exact recursive H2 descent.

Scope: the 179 H4 pairs in the 52 THM-2901 pair-cap exceptions for which

    (h_R-q1(C))/2 <= h_R/7.

At this boundary the inherited binary-link edge high core is not made
finite by discrepancy.  This does *not* make the child intrinsically hard.
For the literal three-label child R, compute its attained singleton maximum
q1(R).  If q1(R)<5h_R/7, THM-2893 applied afresh with
``(k,s,ell)=(3,2,2)`` gives the finite core

    H2(R)={w:c_R(w)>=(h_R-q1(R))/2}.

Every three-cover contains a heavy H2 pair.  Behind each such pair, an
exact longest-component bound reduces the final singleton to a finite
integer interval.  We also retain the parent binary-link sidecar: a leaf
cover counts as unresolved only when all three of its pairs lie in the
inherited link graph.

All computations are exact.  Failed sufficient tests are never reported as
cover witnesses.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
H4_PATH = (
    ROOT
    / "04-computation/lrc14_j6_paircap_exception_h4_membership_census_codex_20260729.py"
)
H4_SHA256 = "63a80908a6380a877345f0cc4aba7a5e0ef2bb3d59b1b10d58367444ed406b75"
BASE_HORIZON = 600


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_h4():
    require(file_sha256(H4_PATH) == H4_SHA256, "H4 membership engine changed")
    spec = importlib.util.spec_from_file_location("h4_nonpositive_h2_base", H4_PATH)
    require(spec is not None and spec.loader is not None, "cannot load H4 engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H4 = load_h4()
E = H4.E


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def sealed_singletons(
    residual: list[tuple[F, F]],
    excluded: frozenset[int],
) -> dict[str, object]:
    """Seal the top three singleton ranks and retain all scanned values."""

    mass = E.interval_mass(residual)
    components = len(residual)
    labels = [
        label
        for label in range(E.FIRST_EXTERNAL, BASE_HORIZON + 1)
        if label not in excluded
    ]
    rows = E.T.coverages_many(residual, labels)
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    require(len(ranked) >= 3, "rank-three head has fewer than three labels")
    q3_base = ranked[2][0]
    gamma = E.S2 * components / 7
    gap = q3_base - mass / 7
    require(gap > 0, "rank-three head does not clear the limiting density")
    tail_first = max(BASE_HORIZON + 1, E.ceiling(gamma / gap))
    if tail_first > BASE_HORIZON + 1:
        extra = [
            label
            for label in range(BASE_HORIZON + 1, tail_first)
            if label not in excluded
        ]
        rows.extend(E.T.coverages_many(residual, extra))
    tail_cap = mass / 7 + gamma / tail_first
    require(tail_cap <= q3_base, "rank-three tail did not seal")
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    require(ranked[2][0] == q3_base, "rank three changed after tail extension")
    by_label = {label: value for value, label in rows}
    for label in dict.fromkeys(
        (labels[0], labels[-1], ranked[0][1], ranked[1][1], ranked[2][1])
    ):
        require(
            by_label[label] == E.T.coverage(residual, label),
            f"singleton scalar/vector mismatch at {label}",
        )
    return {
        "ranked": ranked,
        "by_label": by_label,
        "q1": ranked[0][0],
        "q3": ranked[2][0],
        "top3": tuple(ranked[:3]),
        "tail_first": tail_first,
        "tail_cap": tail_cap,
        "gamma": gamma,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", choices=("nonpositive", "all"), default="nonpositive")
    parser.add_argument("--parent-start", type=int, default=0)
    parser.add_argument("--parent-stop", type=int)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    require(args.parent_start >= 0, "bad parent start")
    require(args.limit is None or args.limit >= 1, "bad limit")
    parent_rows = [H4.exact_row(fields) for fields in H4.parse_exceptions()]
    parent_stop = len(parent_rows) if args.parent_stop is None else args.parent_stop
    require(
        args.parent_start <= parent_stop <= len(parent_rows),
        "bad parent slice",
    )
    totals = Counter()
    parent_rows_out: list[tuple[object, ...]] = []
    minimum_entry = None
    maximum_h2_size = None
    maximum_leaf_horizon = None
    unresolved_rows: list[tuple[object, ...]] = []
    semantic = hashlib.sha256(
        (
            "LRC14/j6/paircap-exception/H4-recursive-H2/v2;"
            f"scope={args.scope};parents={args.parent_start}:{parent_stop}\n"
        ).encode()
    )

    stop = False
    for row in parent_rows[args.parent_start : parent_stop]:
        carrier, components, mass = E.R.CORE.good_norm(
            tuple(sorted((*row["body"], row["apex"])))
        )
        require(
            components == row["components"] and mass == row["mass"],
            "parent carrier reconstruction changed",
        )
        forbidden = frozenset(row["prefix"])
        coverages = {
            label: value
            for value, label in E.T.coverages_many(carrier, list(row["core"]))
        }
        parent_ledger = E.PairLedger(carrier, coverages)
        local = Counter()

        for hpair in combinations(row["core"], 2):
            if args.limit is not None and totals["flags"] >= args.limit:
                stop = True
                break
            pair_union = parent_ledger.edge(*hpair)
            residual = E.R.subtract_local_multi(carrier, hpair)
            residual_mass = E.interval_mass(residual)
            require(residual_mass == mass - pair_union > 0, "H4 residual mismatch")
            link_threshold = residual_mass - row["q1"]
            link_delta = 5 * residual_mass / 7 - row["q1"]
            if args.scope == "nonpositive" and link_delta > 0:
                continue
            totals["flags"] += 1
            local["flags"] += 1
            excluded = frozenset((*forbidden, *hpair))
            singles = sealed_singletons(residual, excluded)
            direct_margin = residual_mass - sum(
                (value for value, _ in singles["top3"]), F(0)
            )
            direct_closed = direct_margin > 0
            totals["direct_closed"] += direct_closed
            local["direct_closed"] += direct_closed
            if direct_closed:
                semantic.update(
                    (
                        f"{row['body']};{row['rank']};{hpair};"
                        f"h={ftext(residual_mass)};d={ftext(link_delta)};"
                        f"S={ftext(direct_margin)};route=direct\n"
                    ).encode()
                )
                continue

            q1_child = singles["q1"]
            entry = 5 * residual_mass / 7 - q1_child
            entry_record = (
                entry,
                row["body"],
                row["rank"],
                hpair,
                q1_child,
                residual_mass,
            )
            if minimum_entry is None or entry_record < minimum_entry:
                minimum_entry = entry_record
            if entry <= 0:
                totals["h2_nonfinite"] += 1
                local["h2_nonfinite"] += 1
                semantic.update(
                    (
                        f"{row['body']};{row['rank']};{hpair};"
                        f"h={ftext(residual_mass)};d={ftext(link_delta)};"
                        f"S={ftext(direct_margin)};entry={ftext(entry)};"
                        "route=H2-nonfinite\n"
                    ).encode()
                )
                continue

            level = (residual_mass - q1_child) / 2
            gap = level - residual_mass / 7
            require(gap == entry / 2 > 0, "H2 gap mismatch")
            cutoff = E.ceiling(singles["gamma"] / gap) - 1
            require(cutoff >= E.FIRST_EXTERNAL, "H2 cutoff below label floor")
            by_label = dict(singles["by_label"])
            scanned_last = singles["tail_first"] - 1
            if cutoff > scanned_last:
                extra = [
                    label
                    for label in range(scanned_last + 1, cutoff + 1)
                    if label not in excluded
                ]
                by_label.update(
                    {
                        label: value
                        for value, label in E.T.coverages_many(residual, extra)
                    }
                )
            require(
                residual_mass / 7 + singles["gamma"] / (cutoff + 1) <= level,
                "H2 discrepancy tail did not seal",
            )
            h2_core = tuple(
                sorted(
                    label
                    for label, value in by_label.items()
                    if label not in excluded and value >= level
                )
            )
            require(len(h2_core) >= 2, "open child has fewer than two H2 labels")
            h2_record = (
                len(h2_core),
                row["body"],
                row["rank"],
                hpair,
                cutoff,
            )
            if maximum_h2_size is None or h2_record > maximum_h2_size:
                maximum_h2_size = h2_record

            heavy_pairs = 0
            leaf_checks = 0
            singleton_covers = 0
            compatible_covers: list[tuple[int, int, int]] = []
            heavy_digest = hashlib.sha256(
                b"LRC14/j6/H4-exception/nonpositive/H2-heavy/v1\n"
            )
            leaf_digest = hashlib.sha256(
                b"LRC14/j6/H4-exception/nonpositive/H2-leaf/v1\n"
            )
            for edge in combinations(h2_core, 2):
                after = E.R.subtract_local_multi(residual, edge)
                after_mass = E.interval_mass(after)
                union = residual_mass - after_mass
                heavy_digest.update(
                    (
                        f"{edge};U={ftext(union)};"
                        f"L={ftext(after_mass)};r={len(after)}\n"
                    ).encode()
                )
                if union < residual_mass - q1_child:
                    continue
                heavy_pairs += 1
                require(after_mass > 0, "H2 pair already covers child")
                longest = max(right - left for left, right in after)
                horizon_fraction = F(1, 7) / longest
                horizon = horizon_fraction.numerator // horizon_fraction.denominator
                horizon_record = (
                    horizon,
                    row["body"],
                    row["rank"],
                    hpair,
                    edge,
                    longest,
                )
                if (
                    maximum_leaf_horizon is None
                    or horizon_record > maximum_leaf_horizon
                ):
                    maximum_leaf_horizon = horizon_record
                leaf_excluded = frozenset((*excluded, *edge))
                labels = [
                    label
                    for label in range(E.FIRST_EXTERNAL, horizon + 1)
                    if label not in leaf_excluded
                ]
                leaf_rows = E.T.coverages_many(after, labels)
                leaf_checks += len(leaf_rows)
                for value, label in leaf_rows:
                    leaf_digest.update(
                        (
                            f"{edge};w={label};c={ftext(value)};"
                            f"h={ftext(after_mass)}\n"
                        ).encode()
                    )
                    if value != after_mass:
                        continue
                    singleton_covers += 1
                    survivor = E.R.subtract_local(after, label)
                    require(
                        not survivor,
                        "singleton coverage equality did not empty the leaf",
                    )
                    # Retain the inherited binary-link triangle sidecar.
                    union_first = residual_mass - E.interval_mass(
                        E.R.subtract_local_multi(residual, (edge[0], label))
                    )
                    union_second = residual_mass - E.interval_mass(
                        E.R.subtract_local_multi(residual, (edge[1], label))
                    )
                    if (
                        union_first >= link_threshold
                        and union_second >= link_threshold
                    ):
                        compatible_covers.append((*edge, label))
                require(
                    F(1, 7 * (horizon + 1)) < longest,
                    "longest-component horizon did not seal",
                )

            closed = not compatible_covers
            if not closed:
                unresolved_rows.append(
                    (
                        row["body"],
                        row["rank"],
                        row["apex"],
                        hpair,
                        residual_mass,
                        q1_child,
                        direct_margin,
                        entry,
                        cutoff,
                        h2_core,
                        tuple(compatible_covers),
                    )
                )
            totals["h2_routes"] += 1
            totals["h2_closed"] += closed
            totals["heavy_pairs"] += heavy_pairs
            totals["leaf_checks"] += leaf_checks
            totals["singleton_covers"] += singleton_covers
            totals["compatible_covers"] += len(compatible_covers)
            local["h2_routes"] += 1
            local["h2_closed"] += closed
            local["heavy_pairs"] += heavy_pairs
            local["leaf_checks"] += leaf_checks
            local["singleton_covers"] += singleton_covers
            local["compatible_covers"] += len(compatible_covers)
            semantic.update(
                (
                    f"{row['body']};{row['rank']};{hpair};"
                    f"h={ftext(residual_mass)};d={ftext(link_delta)};"
                    f"S={ftext(direct_margin)};q1={ftext(q1_child)};"
                    f"entry={ftext(entry)};cutoff={cutoff};H2={h2_core};"
                    f"heavy={heavy_pairs};checks={leaf_checks};"
                    f"covers={tuple(compatible_covers)};"
                    f"heavy_digest={heavy_digest.hexdigest()};"
                    f"leaf_digest={leaf_digest.hexdigest()}\n"
                ).encode()
            )

        parent_rows_out.append(
            (
                row["body"],
                row["rank"],
                local["flags"],
                local["direct_closed"],
                local["h2_routes"],
                local["h2_closed"],
                local["h2_nonfinite"],
                local["heavy_pairs"],
                local["leaf_checks"],
                local["singleton_covers"],
                local["compatible_covers"],
            )
        )
        if stop:
            break

    if (
        args.scope == "nonpositive"
        and args.parent_start == 0
        and parent_stop == len(parent_rows)
        and args.limit is None
    ):
        require(totals["flags"] == 179, "nonpositive-link flag count changed")
    print(f"LRC14 j6 H4 {args.scope} recursive H2 census")
    print(f"totals={tuple(sorted(totals.items()))}")
    print(
        "parent_partition="
        + repr(tuple(sorted(Counter(row[2:] for row in parent_rows_out).items())))
    )
    print(f"minimum_H2_entry={minimum_entry}")
    print(f"maximum_H2_size={maximum_h2_size}")
    print(f"maximum_leaf_horizon={maximum_leaf_horizon}")
    print(f"unresolved_rows={tuple(unresolved_rows)}")
    print(f"semantic_digest={semantic.hexdigest()}")
    print(
        f"scope={args.scope};parents={args.parent_start}:{parent_stop};"
        "direct top-three or recursive H2+singleton leaf;"
        "parent binary-link sidecar retained"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
