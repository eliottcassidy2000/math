#!/usr/bin/env python3
"""Exact branchwise join of finite r4/H1 closure with THM-2905 G5.

This postprocessor does not repeat the interval scan.  It uses:

* the hash-pinned all-hard ledgers parsed by the canonical THM-2905 verifier;
* the proved r4/H1 fact that every eligible branch whose strict discrepancy
  cutoff is at most 15,000 closes; and
* the literal THM-2905 test G5 < h on the same branch keys.

It recomposes whole roots against both the 76-root union through THM-2903 and
the current 82-root union through THM-2905.
"""

from __future__ import annotations

import hashlib
import importlib.util
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
G5_SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.py"
)
G5_SOURCE_SHA256 = (
    "794111b992e912ec8471c8334a867d7b2db1d248f4b08f744f52faf7f50b86c3"
)
H1_SCOUT = Path(__file__).resolve().parent / "scout.py"
H1_SCOUT_SHA256 = (
    "6abc06972cda64bb4f53db26cca62bbea5362ca178bdbf5bf7398e8b0f28317a"
)
H1_MERGE_OUTPUT = (
    Path(__file__).resolve().parent / "merge_all_32_shards.ordinary.out"
)
H1_MERGE_OUTPUT_SHA256 = (
    "dd2d08d523cd481173f7c0d75fb23e458eea935f0f0bd55a5c5d9cac2461f70d"
)
HUNTER_OUTPUT = (
    Path(__file__).resolve().parent / "hunter_all_24_star_survivors.out"
)
HUNTER_OUTPUT_SHA256 = (
    "c142f6389a38549f3b2096b1d4b785f2c19b4062c8879eeeebd2ae04d7be5c7f"
)
S2 = F(99, 70)
MAX_CUTOFF = 15_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_g5():
    require(sha256(G5_SOURCE) == G5_SOURCE_SHA256, "THM-2905 verifier changed")
    require(sha256(H1_SCOUT) == H1_SCOUT_SHA256, "r4/H1 scout changed")
    require(
        sha256(H1_MERGE_OUTPUT) == H1_MERGE_OUTPUT_SHA256,
        "r4/H1 locked merge output changed",
    )
    require(
        sha256(HUNTER_OUTPUT) == HUNTER_OUTPUT_SHA256,
        "r4/H1 Hunter output changed",
    )
    spec = importlib.util.spec_from_file_location("thm2905_g5_join", G5_SOURCE)
    require(spec is not None and spec.loader is not None, "cannot load THM-2905")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def h1_cutoff(row: dict[str, object]) -> int | None:
    r4 = sum(row["qs"][:4], F(0))
    epsilon = row["mass"] - r4 - row["mass"] / 7
    if epsilon <= 0:
        return None
    return ceiling(S2 * row["components"] / (7 * epsilon)) - 1


def body_digest(bodies: set[tuple[int, ...]]) -> str:
    payload = "".join(
        ",".join(map(str, body)) + "\n" for body in sorted(bodies)
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def key_digest(
    keys: set[tuple[tuple[int, ...], int, int]], header: str
) -> str:
    digest = hashlib.sha256((header + "\n").encode())
    for body, rank, apex in sorted(keys):
        digest.update(
            f"{','.join(map(str, body))};{rank};{apex}\n".encode()
        )
    return digest.hexdigest()


def main() -> None:
    g5 = load_g5()
    thm2895, thm2898, thm2899, thm2901, thm2902 = (
        g5.canonical_root_sets()
    )
    g5.thm2903_controls()
    rows = g5.parse_pair_rows(g5.parse_hard_rows())
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_body[row["body"]].append(row)

    h1_keys = {
        (row["body"], row["rank"], row["apex"])
        for row in rows
        if (
            (cutoff := h1_cutoff(row)) is not None
            and cutoff <= MAX_CUTOFF
        )
    }
    g5_keys = {
        (row["body"], row["rank"], row["apex"])
        for row in rows
        if row["margin"] > 0
    }

    def key(row: dict[str, object]) -> tuple[tuple[int, ...], int, int]:
        return row["body"], row["rank"], row["apex"]

    h1_hard_roots = {
        body
        for body, body_rows in by_body.items()
        if all(key(row) in h1_keys for row in body_rows)
    }
    g5_roots = {
        body
        for body, body_rows in by_body.items()
        if all(key(row) in g5_keys for row in body_rows)
    }
    joined_hard_roots = {
        body
        for body, body_rows in by_body.items()
        if all(key(row) in h1_keys or key(row) in g5_keys for row in body_rows)
    }

    # These five roots have no scalar-hard branch and are vacuously closed in
    # the all-root r4/H1 recomposition.
    scalar_roots = set(thm2899)
    h1_roots = h1_hard_roots | scalar_roots
    joined_roots = joined_hard_roots | scalar_roots
    merge_roots = {
        tuple(map(int, line.split("=", 1)[1].split(",")))
        for line in H1_MERGE_OUTPUT.read_text().splitlines()
        if line.startswith("CLOSED_ROOT=")
    }
    require(merge_roots == h1_roots, "recomputed H1 root set differs from merge")

    prior_fifteen = thm2895 | thm2898 | thm2899 | thm2901
    one_hard = {
        body
        for body, body_rows in by_body.items()
        if len(body_rows) == 1 and body_rows[0]["direct_margin"] <= 0
    }
    through_2903 = prior_fifteen | thm2902 | one_hard
    through_2905 = through_2903 | g5_roots

    require(len(rows) == 14_806, "hard branch count changed")
    require(
        (
            len(h1_keys),
            len(g5_keys),
            len(h1_keys & g5_keys),
            len(h1_keys | g5_keys),
        )
        == (5_999, 2_964, 2_909, 6_054),
        "branchwise join counts changed",
    )
    require(
        (
            len(h1_hard_roots),
            len(h1_roots),
            len(g5_roots),
            len(joined_hard_roots),
            len(joined_roots),
        )
        == (127, 132, 16, 129, 134),
        "whole-root join counts changed",
    )
    require(
        joined_roots - h1_roots
        == {
            (1, 2, 3, 4, 5, 8, 9),
            (5, 7, 8, 10, 11, 13, 14),
        },
        "two cross-route roots changed",
    )
    require(
        len(through_2903) == 76 and len(through_2905) == 82,
        "proved baseline changed",
    )
    require(
        (
            len(h1_roots & prior_fifteen),
            len(h1_roots & one_hard),
            len(joined_roots & prior_fifteen),
            len(joined_roots & one_hard),
            len(h1_roots & (through_2905 - through_2903)),
            len(joined_roots & (through_2905 - through_2903)),
        )
        == (9, 21, 10, 21, 6, 6),
        "proved-baseline overlap decomposition changed",
    )
    require(
        (
            len(h1_roots & through_2903),
            len(h1_roots - through_2903),
            len(h1_roots | through_2903),
            len(joined_roots & through_2903),
            len(joined_roots - through_2903),
            len(joined_roots | through_2903),
        )
        == (30, 102, 178, 31, 103, 179),
        "THM-2903 comparison changed",
    )
    require(
        (
            len(h1_roots & through_2905),
            len(h1_roots - through_2905),
            len(h1_roots | through_2905),
            len(joined_roots & through_2905),
            len(joined_roots - through_2905),
            len(joined_roots | through_2905),
        )
        == (36, 96, 178, 37, 97, 179),
        "current comparison changed",
    )

    mixed = []
    mixed_details = []
    for body in sorted(joined_hard_roots):
        classes = []
        for row in by_body[body]:
            row_key = key(row)
            h1 = row_key in h1_keys
            star = row_key in g5_keys
            r4 = sum(row["qs"][:4], F(0))
            epsilon = row["mass"] - r4 - row["mass"] / 7
            classes.append(
                (
                    "both" if h1 and star else "H1" if h1 else "G5",
                    row["rank"],
                    row["apex"],
                )
            )
            if body == (1, 2, 3, 4, 5, 8, 9):
                mixed_details.append(
                    (
                        row["rank"],
                        row["apex"],
                        epsilon,
                        h1_cutoff(row),
                        row["margin"],
                    )
                )
        if any(item[0] == "H1" for item in classes) and any(
            item[0] == "G5" for item in classes
        ):
            mixed.append((body, tuple(classes)))
    require(
        mixed
        == [
            (
                (1, 2, 3, 4, 5, 8, 9),
                (("G5", 1, 22), ("H1", 2, 33)),
            )
        ],
        "genuinely mixed root changed",
    )
    require(
        tuple(mixed_details)
        == (
            (1, 22, F(-193, 2_522_520), None, F(43, 84_084)),
            (2, 33, F(26_111, 2_522_520), 585, F(-1, 10_010)),
        ),
        "genuinely mixed branch margins changed",
    )

    print("LRC14 finite-r4/H1 plus Hunter-star G5 branchwise join")
    print(
        "branch_counts="
        f"(hard={len(rows)},H1={len(h1_keys)},G5={len(g5_keys)},"
        f"both={len(h1_keys & g5_keys)},union={len(h1_keys | g5_keys)},"
        f"H1_only={len(h1_keys - g5_keys)},G5_only={len(g5_keys - h1_keys)})"
    )
    print(
        "root_counts="
        f"(H1_hard={len(h1_hard_roots)},H1_all={len(h1_roots)},"
        f"G5={len(g5_roots)},join_hard={len(joined_hard_roots)},"
        f"join_all={len(joined_roots)})"
    )
    print(f"join_minus_H1={tuple(sorted(joined_roots - h1_roots))}")
    print(f"genuinely_mixed={tuple(mixed)}")
    print(f"mixed_branch_details={tuple(mixed_details)}")
    print(
        "against_T2903="
        f"(baseline={len(through_2903)},H1_intersection={len(h1_roots & through_2903)},"
        f"H1_additive={len(h1_roots - through_2903)},H1_union={len(h1_roots | through_2903)},"
        f"join_intersection={len(joined_roots & through_2903)},"
        f"join_additive={len(joined_roots - through_2903)},"
        f"join_union={len(joined_roots | through_2903)},"
        f"residual={3432 - len(joined_roots | through_2903)})"
    )
    print(
        "baseline_overlap_decomposition="
        f"(H1_prior15={len(h1_roots & prior_fifteen)},"
        f"H1_one_hard={len(h1_roots & one_hard)},"
        f"join_prior15={len(joined_roots & prior_fifteen)},"
        f"join_one_hard={len(joined_roots & one_hard)},"
        f"H1_T2905_additions={len(h1_roots & (through_2905 - through_2903))},"
        f"join_T2905_additions={len(joined_roots & (through_2905 - through_2903))})"
    )
    print(
        "against_T2905="
        f"(baseline={len(through_2905)},H1_intersection={len(h1_roots & through_2905)},"
        f"H1_additive={len(h1_roots - through_2905)},H1_union={len(h1_roots | through_2905)},"
        f"join_intersection={len(joined_roots & through_2905)},"
        f"join_additive={len(joined_roots - through_2905)},"
        f"join_union={len(joined_roots | through_2905)},"
        f"residual={3432 - len(joined_roots | through_2905)})"
    )
    print(f"H1_root_digest={body_digest(h1_roots)}")
    print(f"join_root_digest={body_digest(joined_roots)}")
    print(f"T2903_digest={body_digest(through_2903)}")
    print(f"T2905_digest={body_digest(through_2905)}")
    print(f"H1_additive_T2903_digest={body_digest(h1_roots - through_2903)}")
    print(f"H1_union_T2903_digest={body_digest(h1_roots | through_2903)}")
    print(f"H1_additive_T2905_digest={body_digest(h1_roots - through_2905)}")
    print(f"H1_union_T2905_digest={body_digest(h1_roots | through_2905)}")
    print(f"join_additive_T2903_digest={body_digest(joined_roots - through_2903)}")
    print(f"join_additive_T2905_digest={body_digest(joined_roots - through_2905)}")
    print(f"proved_union_digest={body_digest(joined_roots | through_2905)}")
    print(
        "H1_branch_digest="
        + key_digest(h1_keys, "LRC14/j6/r4-H1-finite-branch-keys/v1")
    )
    print(
        "G5_branch_digest="
        + key_digest(g5_keys, "LRC14/j6/G5-positive-branch-keys/v1")
    )
    print(
        "branch_union_digest="
        + key_digest(
            h1_keys | g5_keys,
            "LRC14/j6/G5-or-finite-H1-branch-keys/v1",
        )
    )
    print("scope=14806 hard branches plus five scalar-terminal roots;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
