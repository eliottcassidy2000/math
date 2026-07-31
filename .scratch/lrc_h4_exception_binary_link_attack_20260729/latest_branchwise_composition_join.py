#!/usr/bin/env python3
"""Join H4 exceptions with G5, THM-2904 pivots, and finite-H1 candidates."""

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
PIVOT_SOURCE = (
    ROOT
    / "04-computation/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.py"
)
PIVOT_SOURCE_SHA256 = (
    "99f1938f264d90c2b34ec3c64566605cc8fd12520424ad2f5cd0957342202ba0"
)
PIVOT_LEDGER = Path(__file__).resolve().parent / "latest_2904_pivot.ledger.out"
PIVOT_LEDGER_SHA256 = (
    "bec35518329b5d9e6ba2c9a8c87bfb20234a0c07dc1a5c5f2babec21888d452a"
)
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


G5 = load("latest_join_g5", G5_PATH, G5_SHA256)
H4 = load("latest_join_h4", H4_PATH, H4_SHA256)
require(sha256(PIVOT_SOURCE) == PIVOT_SOURCE_SHA256, "pivot source changed")
require(sha256(PIVOT_LEDGER) == PIVOT_LEDGER_SHA256, "pivot ledger changed")


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def h1_cutoff(row: dict[str, object]) -> int | None:
    r4 = sum(row["qs"][:4], F(0))
    epsilon = row["mass"] - r4 - row["mass"] / 7
    if epsilon <= 0:
        return None
    return ceiling(S2 * row["components"] / (7 * epsilon)) - 1


def fields(line: str) -> dict[str, str]:
    return {
        item.split("=", 1)[0]: item.split("=", 1)[1]
        for item in line.split(";")[1:]
        if "=" in item
    }


def full_key(row: dict[str, object]):
    return (
        row["body"],
        row["stratum"],
        row["rank"],
        row["apex"],
        row["prefix"],
    )


def pivot_closed_keys():
    keys = set()
    total = 0
    for line in PIVOT_LEDGER.read_text().splitlines():
        if not line.startswith("H1;"):
            continue
        total += 1
        data = fields(line)
        pivots = data["pivot"].split(",")
        require(pivots and all(item for item in pivots), "empty pivot list")
        if all(item.rsplit(":", 2)[1] == "1" for item in pivots):
            keys.add(
                (
                    tuple(map(int, data["E"].split(","))),
                    data["S"],
                    int(data["rank"]),
                    int(data["a"]),
                    tuple(map(int, data["P"].split(","))) if data["P"] else (),
                )
            )
    require(total == 11_842 and len(keys) == 279, "pivot key counts changed")
    return keys


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
    by_key = {}
    for row in rows:
        key = full_key(row)
        require(key not in by_key, "duplicate full hard key")
        by_key[key] = row
        by_body[row["body"]].append(row)

    g5_keys = {full_key(row) for row in rows if row["margin"] > 0}
    h1_keys = {
        full_key(row)
        for row in rows
        if (cutoff := h1_cutoff(row)) is not None and cutoff <= MAX_H1_CUTOFF
    }
    exception_keys = {
        (
            tuple(map(int, data["E"].split(","))),
            data["S"],
            int(data["rank"]),
            int(data["a"]),
            tuple(map(int, data["P"].split(","))) if data["P"] else (),
        )
        for data in H4.parse_exceptions()
    }
    pivot_keys = pivot_closed_keys()
    require(
        len(g5_keys) == 2_964
        and len(h1_keys) == 5_999
        and len(exception_keys) == 52
        and len(pivot_keys) == 279,
        "branch certificate counts changed",
    )
    require(
        pivot_keys <= set(by_key)
        and exception_keys <= set(by_key)
        and not (g5_keys & pivot_keys)
        and not (g5_keys & exception_keys)
        and not (pivot_keys & exception_keys),
        "branch certificate alignment changed",
    )

    def hard_roots(keys):
        return {
            body
            for body, body_rows in by_body.items()
            if all(full_key(row) in keys for row in body_rows)
        }

    routes = {
        "G5": g5_keys,
        "G5P": g5_keys | pivot_keys,
        "G5PE": g5_keys | pivot_keys | exception_keys,
        "H1G5": h1_keys | g5_keys,
        "H1G5P": h1_keys | g5_keys | pivot_keys,
        "H1G5PE": h1_keys | g5_keys | pivot_keys | exception_keys,
    }
    roots = {name: hard_roots(keys) for name, keys in routes.items()}
    scalar_roots = set(thm2899)
    prior_fifteen = thm2895 | thm2898 | thm2899 | thm2901
    one_hard = {
        body
        for body, body_rows in by_body.items()
        if len(body_rows) == 1 and body_rows[0]["direct_margin"] <= 0
    }
    through_2903 = prior_fifteen | thm2902 | one_hard
    through_2905 = through_2903 | roots["G5"]
    through_2904 = through_2905 | roots["G5P"]
    require(
        (len(through_2903), len(through_2905), len(through_2904))
        == (76, 82, 88),
        "proved root baselines changed",
    )
    old_candidate = through_2905 | (roots["H1G5"] | scalar_roots)
    require(len(old_candidate) == 179, "old finite-H1/G5 candidate changed")

    proved_h4 = through_2904 | roots["G5PE"]
    old_candidate_pivot = (
        through_2904 | roots["H1G5P"] | scalar_roots
    )
    old_candidate_pivot_h4 = (
        through_2904 | roots["H1G5PE"] | scalar_roots
    )
    print("LRC14 latest branchwise H4 composition join")
    print(
        "branch_counts="
        f"(G5={len(g5_keys)},pivot={len(pivot_keys)},"
        f"exception={len(exception_keys)},H1={len(h1_keys)},"
        f"G5P={len(routes['G5P'])},G5PE={len(routes['G5PE'])},"
        f"H1G5P={len(routes['H1G5P'])},H1G5PE={len(routes['H1G5PE'])})"
    )
    print(
        "hard_root_counts="
        + repr(tuple((name, len(roots[name])) for name in routes))
    )
    print(
        "proved_T2904_plus_H4="
        f"(baseline={len(through_2904)},"
        f"additive={len(proved_h4-through_2904)},"
        f"union={len(proved_h4)},residual={3432-len(proved_h4)})"
    )
    print(f"proved_H4_additive_roots={tuple(sorted(proved_h4-through_2904))}")
    print(
        "old_candidate_finite_H1_G5="
        f"(baseline={len(old_candidate)},"
        f"plus_pivot={len(old_candidate_pivot)},"
        f"plus_pivot_H4={len(old_candidate_pivot_h4)},"
        f"H4_additive={len(old_candidate_pivot_h4-old_candidate_pivot)},"
        f"residual={3432-len(old_candidate_pivot_h4)})"
    )
    print(
        "old_candidate_H4_additive_roots="
        f"{tuple(sorted(old_candidate_pivot_h4-old_candidate_pivot))}"
    )
    print(
        "digests="
        f"(provedH4={body_digest(proved_h4,'provedH4')},"
        f"candidateP={body_digest(old_candidate_pivot,'candidateP')},"
        f"candidatePE={body_digest(old_candidate_pivot_h4,'candidatePE')})"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
