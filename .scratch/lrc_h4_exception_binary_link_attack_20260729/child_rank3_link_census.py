#!/usr/bin/env python3
"""Exact rank-three and binary-link scout for the 52 pair-cap exceptions.

The parent carrier ``C`` has a globally attained singleton maximum ``q1``
and pair cap ``B2``.  THM-2893 with ``(k,s,ell)=(5,4,2)`` says that every
hypothetical five-cover contains an H4 pair ``L``.  On the literal child

    R = C \ (D_x union D_y),      h_R = |R|,

the remaining three labels form a triangle in the binary link

    yz is an edge iff U_R({y,z}) >= h_R-q1.

This first exact screen retains two independent consequences:

* ``q3(R)+B2<h_R`` excludes every three-cover of R.  The inherited B2 is
  valid because R is a subset of C, and q3(R) is sealed exactly against the
  discrepancy tail.
* if ``(h_R-q1)/2>h_R/7``, every link edge meets the finite high core
  ``A_R={w:c_R(w)>=(h_R-q1)/2}``; hence every link triangle has at least two
  vertices in A_R.  This script records the exact cutoff, but deliberately
  does not enumerate A_R unless ``--enumerate-link-cores`` is supplied.

Failure of either sufficient test is not a covering witness.
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
    spec = importlib.util.spec_from_file_location("h4_exception_link_base", H4_PATH)
    require(spec is not None and spec.loader is not None, "cannot load H4 engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H4 = load_h4()
E = H4.E


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def exact_rank_three(
    residual: list[tuple[F, F]],
    excluded: frozenset[int],
) -> dict[str, object]:
    """Seal the three largest allowed singleton coverages exactly."""

    mass = E.interval_mass(residual)
    components = len(residual)
    require(mass > 0 and components > 0, "empty child residual")
    labels = [
        label
        for label in range(E.FIRST_EXTERNAL, BASE_HORIZON + 1)
        if label not in excluded
    ]
    rows = E.T.coverages_many(residual, labels)
    require(len(rows) == len(labels) >= 3, "rank-three base head changed")
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    q3_base = ranked[2][0]
    gap = q3_base - mass / 7
    require(gap > 0, "rank-three base value does not clear limiting density")
    gamma = E.S2 * components / 7
    tail_first = max(BASE_HORIZON + 1, E.ceiling(gamma / gap))
    if tail_first > BASE_HORIZON + 1:
        extra_labels = [
            label
            for label in range(BASE_HORIZON + 1, tail_first)
            if label not in excluded
        ]
        rows.extend(E.T.coverages_many(residual, extra_labels))
    tail_cap = mass / 7 + gamma / tail_first
    require(tail_cap <= q3_base, "rank-three discrepancy tail did not seal")
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    require(ranked[2][0] == q3_base, "rank three changed after tail extension")
    by_label = {label: value for value, label in rows}
    for label in dict.fromkeys((labels[0], labels[-1], ranked[0][1], ranked[2][1])):
        require(
            by_label[label] == E.T.coverage(residual, label),
            f"rank-three scalar/vector mismatch at {label}",
        )
    return {
        "top3": tuple(ranked[:3]),
        "q3": ranked[2][0],
        "tail_first": tail_first,
        "tail_cap": tail_cap,
        "scanned": len(rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int)
    parser.add_argument("--parent-start", type=int, default=0)
    parser.add_argument("--parent-stop", type=int)
    parser.add_argument("--enumerate-link-cores", action="store_true")
    args = parser.parse_args()
    require(args.limit is None or args.limit >= 1, "bad limit")
    require(args.parent_start >= 0, "bad parent start")

    exception_fields = H4.parse_exceptions()
    parent_rows = [H4.exact_row(fields) for fields in exception_fields]
    pair_caps = {
        (
            tuple(map(int, fields["E"].split(","))),
            int(fields["rank"]),
        ): F(fields["B2"])
        for fields in exception_fields
    }
    counters = Counter()
    per_parent: list[tuple[object, ...]] = []
    extrema: dict[str, tuple[object, ...] | None] = {
        "minimum_rank_pair_margin": None,
        "minimum_scalar_margin": None,
        "largest_rank3_tail": None,
        "largest_link_cutoff": None,
    }
    semantic = hashlib.sha256(
        b"LRC14/j6/paircap-exception/H4-child-rank3-link/v1\n"
    )

    stop = args.limit
    parent_stop = len(parent_rows) if args.parent_stop is None else args.parent_stop
    require(
        args.parent_start <= parent_stop <= len(parent_rows),
        "bad parent slice",
    )
    for parent_index, row in enumerate(
        parent_rows[args.parent_start : parent_stop],
        start=args.parent_start,
    ):
        carrier, components, mass = E.R.CORE.good_norm(
            tuple(sorted((*row["body"], row["apex"])))
        )
        require(
            components == row["components"] and mass == row["mass"],
            "parent carrier reconstruction changed",
        )
        forbidden = frozenset(row["prefix"])
        parent_coverages = {
            label: value
            for value, label in E.T.coverages_many(carrier, list(row["core"]))
        }
        require(set(parent_coverages) == set(row["core"]), "H4 core changed")
        parent_pair_ledger = E.PairLedger(carrier, parent_coverages)
        B2 = pair_caps[(row["body"], row["rank"])]
        local = Counter()

        for pair_index, pair in enumerate(combinations(row["core"], 2)):
            if stop is not None and counters["pairs"] >= stop:
                break
            pair_union = parent_pair_ledger.edge(*pair)
            residual = E.R.subtract_local_multi(carrier, pair)
            residual_mass = E.interval_mass(residual)
            require(residual_mass == mass - pair_union > 0, "pair residual mismatch")
            direct, direct_components, direct_mass = E.R.CORE.good_norm(
                tuple(sorted((*row["body"], row["apex"], *pair)))
            )
            require(
                residual == direct
                and len(residual) == direct_components
                and residual_mass == direct_mass,
                "literal/direct H4 child mismatch",
            )
            excluded = frozenset((*forbidden, *pair))
            rank3 = exact_rank_three(residual, excluded)
            scalar_margin = residual_mass - sum(
                (value for value, _ in rank3["top3"]), F(0)
            )
            rank_pair_margin = residual_mass - B2 - rank3["q3"]
            scalar_closed = scalar_margin > 0
            rank_pair_closed = rank_pair_margin > 0
            closed = scalar_closed or rank_pair_closed
            link_delta = 5 * residual_mass / 7 - row["q1"]
            finite_link = link_delta > 0
            link_cutoff = None
            link_core_size = None
            high_link_edges = None
            link_closed = False
            if finite_link:
                threshold = (residual_mass - row["q1"]) / 2
                gap = threshold - residual_mass / 7
                require(gap == link_delta / 2 > 0, "link-core gap mismatch")
                gamma = E.S2 * len(residual) / 7
                link_cutoff = E.ceiling(gamma / gap) - 1
                require(link_cutoff >= E.FIRST_EXTERNAL, "link cutoff below floor")
                if args.enumerate_link_cores and not closed:
                    link_labels = [
                        label
                        for label in range(E.FIRST_EXTERNAL, link_cutoff + 1)
                        if label not in excluded
                    ]
                    link_rows = E.T.coverages_many(residual, link_labels)
                    link_core = tuple(
                        label for value, label in link_rows if value >= threshold
                    )
                    link_core_size = len(link_core)
                    require(
                        residual_mass / 7 + gamma / (link_cutoff + 1)
                        <= threshold,
                        "link-core discrepancy tail did not seal",
                    )
                    link_core_set = set(link_core)
                    link_singletons = {
                        label: value
                        for value, label in link_rows
                        if label in link_core_set
                    }
                    link_ledger = E.PairLedger(residual, link_singletons)
                    high_link_edges = 0
                    first_edge = None
                    for high_pair in combinations(link_core, 2):
                        counters["high_pairs_tested"] += 1
                        local["high_pairs_tested"] += 1
                        union = link_ledger.edge(*high_pair)
                        if union >= residual_mass - row["q1"]:
                            high_link_edges += 1
                            if first_edge is None:
                                first_edge = (high_pair, union)
                    if first_edge is not None:
                        high_pair, union = first_edge
                        require(
                            union
                            == residual_mass
                            - E.interval_mass(
                                E.R.subtract_local_multi(residual, high_pair)
                            ),
                            "high-link edge/direct subtraction mismatch",
                        )
                    link_closed = high_link_edges == 0
                    counters["link_core_flags"] += 1
                    counters["link_core_vertices"] += link_core_size
                    counters["high_link_edges"] += high_link_edges
                    counters["link_closed_zero_high_edge"] += link_closed
                    local["link_core_flags"] += 1
                    local["link_core_vertices"] += link_core_size
                    local["high_link_edges"] += high_link_edges
                    local["link_closed_zero_high_edge"] += link_closed

            counters["pairs"] += 1
            counters["scalar_closed"] += scalar_closed
            counters["rank_pair_closed"] += rank_pair_closed
            counters["closed_union"] += closed
            counters["closed_with_link"] += (closed or link_closed)
            counters["open"] += not closed
            counters["open_after_link"] += not (closed or link_closed)
            counters["finite_link"] += finite_link
            counters["open_finite_link"] += (not closed and finite_link)
            counters["open_nonfinite_link"] += (not closed and not finite_link)
            counters["rank3_scanned"] += rank3["scanned"]
            local["pairs"] += 1
            local["closed"] += closed
            local["closed_with_link"] += (closed or link_closed)
            local["open"] += not closed
            local["open_after_link"] += not (closed or link_closed)
            local["finite_link"] += finite_link
            local["open_finite_link"] += (not closed and finite_link)
            local["open_nonfinite_link"] += (not closed and not finite_link)

            margin_record = (
                rank_pair_margin,
                row["body"],
                row["rank"],
                pair,
                rank3["top3"],
                B2,
            )
            if (
                extrema["minimum_rank_pair_margin"] is None
                or margin_record < extrema["minimum_rank_pair_margin"]
            ):
                extrema["minimum_rank_pair_margin"] = margin_record
            scalar_record = (
                scalar_margin,
                row["body"],
                row["rank"],
                pair,
                rank3["top3"],
            )
            if (
                extrema["minimum_scalar_margin"] is None
                or scalar_record < extrema["minimum_scalar_margin"]
            ):
                extrema["minimum_scalar_margin"] = scalar_record
            tail_record = (
                rank3["tail_first"],
                row["body"],
                row["rank"],
                pair,
            )
            if (
                extrema["largest_rank3_tail"] is None
                or tail_record > extrema["largest_rank3_tail"]
            ):
                extrema["largest_rank3_tail"] = tail_record
            if link_cutoff is not None:
                link_record = (
                    link_cutoff,
                    row["body"],
                    row["rank"],
                    pair,
                    link_core_size,
                )
                if (
                    extrema["largest_link_cutoff"] is None
                    or link_record > extrema["largest_link_cutoff"]
                ):
                    extrema["largest_link_cutoff"] = link_record

            semantic.update(
                (
                    f"{row['body']};{row['rank']};{pair};"
                    f"h={ftext(residual_mass)};"
                    f"top3={','.join(f'{label}:{ftext(value)}' for value, label in rank3['top3'])};"
                    f"B2={ftext(B2)};S={ftext(scalar_margin)};"
                    f"Q={ftext(rank_pair_margin)};d={ftext(link_delta)};"
                    f"rt={rank3['tail_first']};lt={link_cutoff};"
                    f"A={link_core_size};eA={high_link_edges};"
                    f"LC={int(link_closed)}\n"
                ).encode()
            )
        per_parent.append(
            (
                row["body"],
                row["rank"],
                local["pairs"],
                local["closed"],
                local["closed_with_link"],
                local["open"],
                local["open_after_link"],
                local["finite_link"],
                local["open_finite_link"],
                local["open_nonfinite_link"],
                local["link_core_flags"],
                local["link_core_vertices"],
                local["high_link_edges"],
                local["link_closed_zero_high_edge"],
            )
        )
        if stop is not None and counters["pairs"] >= stop:
            break

    print("LRC14 j6 pair-cap-exception H4 child rank3/link census")
    print(f"counts={tuple(sorted(counters.items()))}")
    print(
        "parent_partition="
        + repr(tuple(sorted(Counter(row[2:] for row in per_parent).items())))
    )
    for name in sorted(extrema):
        print(f"{name}={extrema[name]}")
    print(f"semantic_digest={semantic.hexdigest()}")
    print(
        f"mode={'LINK-CORE' if args.enumerate_link_cores else 'RANK3'};"
        f"base_horizon={BASE_HORIZON};"
        f"parent_slice={args.parent_start}:{parent_stop};"
        "failed certificates are not cover witnesses"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
