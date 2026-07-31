#!/usr/bin/env python3
"""Compose the 52 discharged H4-exception branches into whole-root ledgers.

Two baselines are audited:

1. the current proved 82-root union through THM-2905, with G5 as the other
   branch certificate;
2. the candidate 179-root branchwise finite-H1/G5 union, with the same 52
   exception discharges added as a third branch certificate.

This script performs set joins only.  It relies on the separate exact H4 child
census and boundary-tiler audit for the soundness of the 52 exception keys.
"""

from __future__ import annotations

import hashlib
import importlib.util
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
G5_PATH = (
    ROOT
    / "04-computation/lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.py"
)
G5_SHA256 = "794111b992e912ec8471c8334a867d7b2db1d248f4b08f744f52faf7f50b86c3"
H4_PATH = (
    ROOT
    / "04-computation/lrc14_j6_paircap_exception_h4_membership_census_codex_20260729.py"
)
H4_SHA256 = "63a80908a6380a877345f0cc4aba7a5e0ef2bb3d59b1b10d58367444ed406b75"
S2 = F(99, 70)
MAX_H1_CUTOFF = 15_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name: str, path: Path, expected: str):
    require(sha256(path) == expected, f"{path.name} changed")
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


G5 = load("h4_root_join_g5", G5_PATH, G5_SHA256)
H4 = load("h4_root_join_membership", H4_PATH, H4_SHA256)


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def h1_cutoff(row: dict[str, object]) -> int | None:
    r4 = sum(row["qs"][:4], F(0))
    epsilon = row["mass"] - r4 - row["mass"] / 7
    if epsilon <= 0:
        return None
    return ceiling(S2 * row["components"] / (7 * epsilon)) - 1


def body_digest(bodies: set[tuple[int, ...]], tag: str) -> str:
    digest = hashlib.sha256((tag + "\n").encode())
    for body in sorted(bodies):
        digest.update((",".join(map(str, body)) + "\n").encode())
    return digest.hexdigest()


def main() -> None:
    thm2895, thm2898, thm2899, thm2901, thm2902 = G5.canonical_root_sets()
    G5.thm2903_controls()
    rows = G5.parse_pair_rows(G5.parse_hard_rows())
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    by_key: dict[tuple[tuple[int, ...], int, int], dict[str, object]] = {}
    for row in rows:
        key = (row["body"], row["rank"], row["apex"])
        require(key not in by_key, "duplicate hard branch key")
        by_key[key] = row
        by_body[row["body"]].append(row)

    exception_keys = {
        (
            tuple(map(int, fields["E"].split(","))),
            int(fields["rank"]),
            int(fields["a"]),
        )
        for fields in H4.parse_exceptions()
    }
    require(
        len(exception_keys) == 52
        and len({body for body, _, _ in exception_keys}) == 52,
        "exception key/body count changed",
    )
    require(
        exception_keys <= set(by_key)
        and all(by_key[key]["pair_margin"] <= 0 for key in exception_keys),
        "exception key does not match the all-hard pair-cap ledger",
    )

    g5_keys = {
        key for key, row in by_key.items() if row["margin"] > 0
    }
    h1_keys = {
        key
        for key, row in by_key.items()
        if (
            (cutoff := h1_cutoff(row)) is not None
            and cutoff <= MAX_H1_CUTOFF
        )
    }
    require(
        (len(g5_keys), len(h1_keys), len(g5_keys | h1_keys))
        == (2964, 5999, 6054),
        "G5/H1 branch counts changed",
    )
    require(not (g5_keys & exception_keys), "G5 unexpectedly closed an exception")

    def hard_roots(keys: set[tuple[tuple[int, ...], int, int]]):
        return {
            body
            for body, body_rows in by_body.items()
            if all((row["body"], row["rank"], row["apex"]) in keys for row in body_rows)
        }

    # G5 is a certificate only on the scalar-hard ledger.  The five
    # scalar-terminal roots in THM-2899 enter through the proved root
    # baselines, not through G5 itself.  The finite-H1/G5 join's published
    # all-root count does adjoin those five vacuous scalar roots.
    scalar_roots = set(thm2899)
    g5_roots = hard_roots(g5_keys)
    g5_exception_roots = hard_roots(g5_keys | exception_keys)
    h1_g5_roots = hard_roots(h1_keys | g5_keys) | scalar_roots
    h1_g5_exception_roots = (
        hard_roots(h1_keys | g5_keys | exception_keys) | scalar_roots
    )
    require(
        (len(g5_roots), len(h1_g5_roots)) == (16, 134),
        "baseline recomposed root counts changed",
    )

    prior_fifteen = thm2895 | thm2898 | thm2899 | thm2901
    one_hard = {
        body
        for body, body_rows in by_body.items()
        if len(body_rows) == 1 and body_rows[0]["direct_margin"] <= 0
    }
    through_2903 = prior_fifteen | thm2902 | one_hard
    through_2905 = through_2903 | g5_roots
    candidate_2904 = through_2905 | h1_g5_roots
    require(
        (len(through_2903), len(through_2905), len(candidate_2904))
        == (76, 82, 179),
        "proved/candidate baseline changed",
    )

    current_additions = g5_exception_roots - through_2905
    current_union = through_2905 | g5_exception_roots
    candidate_additions = h1_g5_exception_roots - candidate_2904
    candidate_union = candidate_2904 | h1_g5_exception_roots

    print("LRC14 H4-exception whole-root consequence join")
    print(
        "branch_counts="
        f"(hard={len(rows)},G5={len(g5_keys)},H1={len(h1_keys)},"
        f"H1G5={len(h1_keys | g5_keys)},exceptions={len(exception_keys)},"
        f"exception_in_H1={len(exception_keys & h1_keys)},"
        f"G5_exception_union={len(g5_keys | exception_keys)},"
        f"H1G5_exception_union={len(h1_keys | g5_keys | exception_keys)})"
    )
    print(
        "root_certificate_counts="
        f"(G5={len(g5_roots)},G5_exception={len(g5_exception_roots)},"
        f"H1G5={len(h1_g5_roots)},"
        f"H1G5_exception={len(h1_g5_exception_roots)})"
    )
    print(
        "against_T2905="
        f"(baseline={len(through_2905)},"
        f"certificate_intersection={len(g5_exception_roots & through_2905)},"
        f"additive={len(current_additions)},union={len(current_union)},"
        f"residual={3432-len(current_union)})"
    )
    print(f"T2905_additive_roots={tuple(sorted(current_additions))}")
    print(
        "after_candidate_T2904="
        f"(baseline={len(candidate_2904)},"
        f"certificate_intersection={len(h1_g5_exception_roots & candidate_2904)},"
        f"additive={len(candidate_additions)},union={len(candidate_union)},"
        f"residual={3432-len(candidate_union)})"
    )
    print(f"T2904_additive_roots={tuple(sorted(candidate_additions))}")
    print(
        "exception_body_root_intersections="
        f"(T2905={len({b for b,_,_ in exception_keys} & through_2905)},"
        f"T2904={len({b for b,_,_ in exception_keys} & candidate_2904)},"
        f"G5_exception={len({b for b,_,_ in exception_keys} & g5_exception_roots)},"
        f"H1G5_exception={len({b for b,_,_ in exception_keys} & h1_g5_exception_roots)})"
    )
    print(
        "digests="
        f"(G5E={body_digest(g5_exception_roots, 'G5E')},"
        f"H1G5E={body_digest(h1_g5_exception_roots, 'H1G5E')},"
        f"current={body_digest(current_union, 'current')},"
        f"candidate={body_digest(candidate_union, 'candidate')})"
    )
    print(
        "scope=52 discharged pair-cap-exception branches;"
        "whole-root joins against current and candidate branch baselines"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
