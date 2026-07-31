#!/usr/bin/env python3
"""Exact projected-safe prefix filter on the subcritical k=5 pair bank."""

from __future__ import annotations

import argparse
import hashlib
import multiprocessing as mp
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
CHILD = HERE.parent / "lrc_j7_two_drift_transition_child_20260729"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CHILD))

import projected_safe_set_pair_filter as P  # noqa: E402
import subcritical_exact_pair_bank as B  # noqa: E402


def filtered_profile(body: tuple[int, ...]) -> dict[str, object]:
    profile = B.profile(body)
    pairs = profile["exact_pairs"]
    if not pairs:
        return {
            **profile,
            "killed": 0,
            "survivors": [],
            "maximum_prefix": 0,
            "minimum_certificate": None,
        }

    carrier_i = B.A.integer_carrier(body)
    canonical_l = pairs[0][3]
    cells = P.body_cells(carrier_i, canonical_l)
    killed = 0
    survivors: list[tuple[object, ...]] = []
    maximum_prefix = 0
    minimum_certificate: tuple[object, ...] | None = None

    for pair in pairs:
        (
            pair_body,
            first,
            second,
            pair_l,
            h,
            components,
            delta_first,
            gap,
            second_floor,
            second_cap,
            delta_second,
        ) = pair
        B.require(pair_body == body and pair_l == canonical_l, "pair mistyped")
        projected_lower, cells_used = P.projected_safe_lower_bound(
            cells,
            canonical_l,
            first,
            second,
        )
        maximum_prefix = max(maximum_prefix, cells_used)
        certificate = (
            body,
            first,
            second,
            canonical_l,
            h,
            components,
            delta_first,
            gap,
            second_floor,
            second_cap,
            delta_second,
            projected_lower,
            cells_used,
        )
        if projected_lower >= P.ALIGNED_UNION_CAP:
            killed += 1
            margin = projected_lower - P.ALIGNED_UNION_CAP
            candidate = (margin, certificate)
            if minimum_certificate is None or candidate < minimum_certificate:
                minimum_certificate = candidate
        else:
            survivors.append(certificate)

    return {
        **profile,
        "killed": killed,
        "survivors": survivors,
        "maximum_prefix": maximum_prefix,
        "minimum_certificate": minimum_certificate,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    args = parser.parse_args()
    B.require(args.workers >= 1, "worker count must be positive")
    roots = tuple(B.combinations(range(1, 15), 6))
    if args.workers == 1:
        profiles = [filtered_profile(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(filtered_profile, roots, chunksize=1))
    profiles.sort(key=lambda row: row["body"])

    analytic_rows = sum(row["analytic_rows"] for row in profiles)
    pair_candidates = sum(row["pair_candidates"] for row in profiles)
    suffix_rows = sorted(
        row
        for profile in profiles
        for row in profile["suffix_rows"]
    )
    exact_pairs = sorted(
        row
        for profile in profiles
        for row in profile["exact_pairs"]
    )
    killed = sum(row["killed"] for row in profiles)
    survivors = sorted(
        row
        for profile in profiles
        for row in profile["survivors"]
    )
    maximum_prefix = max(row["maximum_prefix"] for row in profiles)
    certificates = [
        row["minimum_certificate"]
        for row in profiles
        if row["minimum_certificate"] is not None
    ]
    minimum_certificate = min(certificates, default=None)

    bank_digest = hashlib.sha256(
        b"LRC14/k5/subcritical-exact-pair-bank/v1\n"
        + repr(tuple(suffix_rows)).encode()
        + b"\n"
        + repr(tuple(exact_pairs)).encode()
    ).hexdigest()
    survivor_digest = hashlib.sha256(
        b"LRC14/k5/projected-safe-subcritical-pairs/v1\n"
        + repr(tuple(survivors)).encode()
    ).hexdigest()

    print("LRC14 k5 subcritical exact-pair projected-safe prefix filter")
    print(
        f"analytic_first_rows={analytic_rows};"
        f"projected_suffix_rows={len(suffix_rows)};"
        f"suffix_roots={len({row[0] for row in suffix_rows})}"
    )
    print(
        f"finite_z2_candidates={pair_candidates};"
        f"exact_excess_admissible_pairs={len(exact_pairs)};"
        f"pair_roots={len({row[0] for row in exact_pairs})};"
        f"pair_(E,z1)_rows={len({(row[0],row[1]) for row in exact_pairs})}"
    )
    print(
        f"killed_by_projected_prefix={killed};"
        f"surviving_pairs={len(survivors)};"
        f"surviving_roots={len({row[0] for row in survivors})}"
    )
    print(f"maximum_cell_prefix={maximum_prefix}")
    print(f"minimum_certificate={minimum_certificate}")
    print(f"bank_digest={bank_digest}")
    print(f"survivor_digest={survivor_digest}")
    if survivors:
        print(f"survivors={tuple(survivors)}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
