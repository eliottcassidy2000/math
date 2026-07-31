#!/usr/bin/env python3
"""Exact first diagnostics for the 52-exception H4 pair-link route.

For each actual H4 pair L, compute the literal residual mass h_L and test:

* the parent-singleton scalar closure 3*q1 < h_L;
* the inherited pair-partition closure q1+B2 < h_L; and
* finiteness of the binary link high core,
      (h_L-q1)/2 > h_L/7  iff  5*h_L/7-q1 > 0.

This is a discovery census only.  Failure of every printed test is not a
cover witness.
"""

from __future__ import annotations

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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_h4():
    require(file_sha256(H4_PATH) == H4_SHA256, "H4 membership engine changed")
    spec = importlib.util.spec_from_file_location("h4_link_membership", H4_PATH)
    require(spec is not None and spec.loader is not None, "cannot load H4 engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H4 = load_h4()
E = H4.E


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    exception_fields = H4.parse_exceptions()
    parent_rows = [H4.exact_row(fields) for fields in exception_fields]
    pair_caps = {
        (
            tuple(map(int, fields["E"].split(","))),
            int(fields["rank"]),
        ): F(fields["B2"])
        for fields in exception_fields
    }
    route_counts = Counter()
    sign_by_branch = []
    closest_positive = None
    closest_nonpositive = None
    finite_cutoffs = []
    finite_components = []
    semantic = hashlib.sha256(
        b"LRC14/j6/paircap-exception/H4-pair-link-delta/v1\n"
    )

    for row in parent_rows:
        carrier, components, mass = E.R.CORE.good_norm(
            tuple(sorted((*row["body"], row["apex"])))
        )
        require(
            components == row["components"] and mass == row["mass"],
            "parent carrier reconstruction changed",
        )
        core = tuple(row["core"])
        coverages = {
            label: value
            for value, label in E.T.coverages_many(carrier, list(core))
        }
        require(set(coverages) == set(core), "H4 singleton universe changed")
        pair_ledger = E.PairLedger(carrier, coverages)
        branch_signs = Counter()

        for pair in combinations(core, 2):
            union = pair_ledger.edge(*pair)
            residual_mass = mass - union
            require(residual_mass >= 0, "negative literal residual mass")
            scalar_margin = residual_mass - 3 * row["q1"]
            # The exact parent B2 is not part of the H4 membership row, so
            # recover it from the pinned THM-2901 exception ledger.
            inherited_pair_margin = (
                residual_mass
                - row["q1"]
                - pair_caps[(row["body"], row["rank"])]
            )
            link_delta = 5 * residual_mass / 7 - row["q1"]
            scalar_closed = scalar_margin > 0
            pair_closed = inherited_pair_margin > 0
            finite_link = link_delta > 0
            equality = link_delta == 0
            route_counts["pairs"] += 1
            route_counts["scalar_closed"] += scalar_closed
            route_counts["pair_closed"] += pair_closed
            route_counts["finite_link"] += finite_link
            route_counts["link_equality"] += equality
            branch_signs["positive" if finite_link else "nonpositive"] += 1

            record = (
                abs(link_delta),
                row["body"],
                row["rank"],
                pair,
                residual_mass,
                link_delta,
            )
            if finite_link:
                residual = E.R.subtract_local_multi(carrier, pair)
                require(
                    E.interval_mass(residual) == residual_mass,
                    "literal pair residual mass mismatch",
                )
                residual_components = len(residual)
                beta = link_delta / 2
                gamma = E.S2 * residual_components / 7
                cutoff = E.ceiling(gamma / beta) - 1
                require(cutoff >= E.FIRST_EXTERNAL, "link cutoff below label floor")
                finite_cutoffs.append(cutoff)
                finite_components.append(residual_components)
                if closest_positive is None or record < closest_positive:
                    closest_positive = record
            elif closest_nonpositive is None or record < closest_nonpositive:
                closest_nonpositive = record

            semantic.update(
                (
                    f"{row['body']};{row['rank']};{pair};"
                    f"{ftext(residual_mass)};{ftext(scalar_margin)};"
                    f"{ftext(inherited_pair_margin)};{ftext(link_delta)}\n"
                ).encode()
            )

        sign_by_branch.append(
            (
                row["body"],
                row["rank"],
                len(core),
                branch_signs["positive"],
                branch_signs["nonpositive"],
            )
        )

    require(route_counts["pairs"] == 18_290, "H4 pair count changed")
    finite_cutoffs.sort()
    finite_components.sort()
    quantile_indices = (0, 25, 50, 75, 90, 95, 99, 100)
    cutoff_quantiles = tuple(
        (
            percentile,
            finite_cutoffs[
                min(
                    len(finite_cutoffs) - 1,
                    percentile * (len(finite_cutoffs) - 1) // 100,
                )
            ],
        )
        for percentile in quantile_indices
    )
    component_quantiles = tuple(
        (
            percentile,
            finite_components[
                min(
                    len(finite_components) - 1,
                    percentile * (len(finite_components) - 1) // 100,
                )
            ],
        )
        for percentile in quantile_indices
    )
    print("LRC14 j6 pair-cap-exception H4 pair-link delta scout")
    print(f"counts={tuple(sorted(route_counts.items()))}")
    print(
        "branch_partition="
        + repr(
            tuple(
                sorted(
                    Counter(
                        (positive, nonpositive)
                        for _, _, _, positive, nonpositive in sign_by_branch
                    ).items()
                )
            )
        )
    )
    print(f"closest_positive={closest_positive}")
    print(f"closest_nonpositive={closest_nonpositive}")
    print(f"finite_cutoff_quantiles={cutoff_quantiles}")
    print(f"finite_component_quantiles={component_quantiles}")
    print(f"semantic_digest={semantic.hexdigest()}")
    print(
        "scope=52 pair-cap exceptions;all 18,290 actual H4 pair flags;"
        "delta diagnostics only;failed tests are not cover witnesses"
    )


if __name__ == "__main__":
    main()
