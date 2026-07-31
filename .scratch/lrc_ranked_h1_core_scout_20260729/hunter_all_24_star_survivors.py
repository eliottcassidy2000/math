#!/usr/bin/env python3
"""Exact relative-Hunter audit of all 24 finite-core star survivors."""

from __future__ import annotations

import hashlib
import importlib.util
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
AUDIT = HERE / "hunter_two_star_exceptions.py"
AUDIT_SHA = "c452781c7cf6a2ada6be8984b9c9cbe7aab8369c5a9333721ef8ecc8e0207393"
TARGETS = (
    ((2, 3, 4, 5, 10, 12, 13), 2, 16),
    ((1, 2, 3, 5, 6, 8, 13), 1, 21),
    ((1, 2, 3, 4, 6, 7, 9), 3, 24),
    ((2, 3, 5, 8, 9, 11, 12), 4, 28),
    ((1, 4, 5, 8, 9, 11, 13), 3, 28),
    ((1, 2, 5, 7, 9, 11, 12), 3, 32),
    ((1, 3, 4, 8, 9, 10, 11), 3, 28),
    ((1, 3, 6, 9, 11, 13, 14), 3, 17),
    ((1, 2, 5, 7, 9, 11, 13), 6, 48),
    ((1, 2, 4, 6, 10, 11, 13), 4, 21),
    ((1, 4, 6, 7, 8, 11, 12), 2, 19),
    ((1, 3, 4, 8, 9, 12, 13), 3, 28),
    ((1, 2, 3, 5, 8, 10, 13), 2, 36),
    ((1, 2, 5, 9, 11, 12, 13), 2, 35),
    ((1, 2, 3, 4, 6, 9, 11), 2, 24),
    ((1, 2, 3, 4, 11, 12, 13), 1, 28),
    ((1, 2, 3, 6, 9, 11, 12), 2, 42),
    ((1, 2, 5, 7, 8, 9, 11), 5, 25),
    ((1, 2, 3, 6, 7, 9, 14), 1, 24),
    ((1, 3, 5, 6, 7, 11, 13), 5, 24),
    ((1, 2, 4, 11, 12, 13, 14), 2, 15),
    ((1, 3, 4, 9, 10, 11, 12), 4, 28),
    ((1, 3, 4, 6, 7, 11, 13), 3, 18),
    ((2, 4, 5, 6, 8, 10, 11), 1, 18),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_audit():
    require(
        hashlib.sha256(AUDIT.read_bytes()).hexdigest() == AUDIT_SHA,
        "Hunter audit source changed",
    )
    spec = importlib.util.spec_from_file_location("hunter_base", AUDIT)
    require(spec is not None and spec.loader is not None, "cannot load audit")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H = load_audit()
ORIGINAL_PROFILE = H.S.profile_root_task


def profile_without_high_k_requirement(task):
    body, max_cutoff, max_combinations, _ = task
    return ORIGINAL_PROFILE(
        (body, max_cutoff, max_combinations, False)
    )


H.S.profile_root_task = profile_without_high_k_requirement


def main() -> None:
    require(
        len(TARGETS) == 24 and len(set(TARGETS)) == 24,
        "star-survivor target battery changed",
    )
    rows = [H.audit_target(*target) for target in TARGETS]
    statuses = Counter(
        "HUNTER_CLOSED" if row["hunter_hard_sets"] == 0
        else "HUNTER_OPEN"
        for row in rows
    )
    aggregate = hashlib.sha256(
        b"LRC14/j6/H1-relative-Hunter/all-star-survivors/v1\n"
    )
    print("All-root finite H1 relative-Hunter star-survivor audit")
    print(f"targets={len(rows)};statuses={tuple(sorted(statuses.items()))}")
    for row in rows:
        aggregate.update((row["ledger"] + "\n").encode())
        print(
            f"E={row['body']};rank={row['rank']};a={row['apex']};"
            f"P={row['prefix']};h={H.ftext(row['h'])};H={row['H']};"
            f"pivot={row['pivot']};nodes={row['nodes']};"
            f"star_prunes={row['star_bound_prunes']};"
            f"star_safe_leaves={row['star_safe_leaves']};"
            f"star_hostile={row['star_hostile_sets']};"
            f"hunter_repairs={row['hunter_repairs']};"
            f"hunter_hard={row['hunter_hard_sets']};"
            f"pairs={row['pair_evaluations']};"
            f"maxPsi={H.ftext(row['max_hunter'])};"
            f"margin={H.ftext(row['margin'])};"
            f"min_extra_credit={H.ftext(row['min_extra_credit'])};"
            f"maxrow={row['max_row']};ledger={row['ledger']}"
        )
    print(f"aggregate_sha256={aggregate.hexdigest()}")
    print("scope=24 finite-core ordered-star survivors;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
