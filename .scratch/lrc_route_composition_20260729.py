#!/usr/bin/env python3
"""Discovery join of the four live exact j=6 branch certificates."""

from __future__ import annotations

import hashlib
import importlib.util
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
G5_PATH = ROOT / "04-computation/lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.py"
H4_PATH = ROOT / "04-computation/lrc14_j6_paircap_exception_h4_membership_census_codex_20260729.py"
PIVOT_PATH = ROOT / ".scratch/lrc_h4_exception_binary_link_attack_20260729/latest_2904_pivot.ledger.out"
G5_SHA = "794111b992e912ec8471c8334a867d7b2db1d248f4b08f744f52faf7f50b86c3"
H4_SHA = "63a80908a6380a877345f0cc4aba7a5e0ef2bb3d59b1b10d58367444ed406b75"
PIVOT_SHA = "bec35518329b5d9e6ba2c9a8c87bfb20234a0c07dc1a5c5f2babec21888d452a"
S2 = F(99, 70)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load(name: str, path: Path, digest: str):
    require(hashlib.sha256(path.read_bytes()).hexdigest() == digest, f"{name} hash")
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"load {name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def h1_cutoff(row: dict[str, object]) -> int | None:
    epsilon = row["mass"] - sum(row["qs"][:4], F(0)) - row["mass"] / 7
    if epsilon <= 0:
        return None
    return ceiling(S2 * row["components"] / (7 * epsilon)) - 1


def short_key(row: dict[str, object]):
    return row["body"], row["rank"], row["apex"], row["prefix"]


def full_key(row: dict[str, object]):
    return row["body"], row["gate_size"], row["rank"], row["apex"], row["prefix"]


def digest_keys(keys, tag: str) -> str:
    digest = hashlib.sha256((tag + "\n").encode())
    for key in sorted(keys):
        digest.update((repr(key) + "\n").encode())
    return digest.hexdigest()


def main() -> None:
    g5 = load("route_join_g5", G5_PATH, G5_SHA)
    h4 = load("route_join_h4", H4_PATH, H4_SHA)
    require(hashlib.sha256(PIVOT_PATH.read_bytes()).hexdigest() == PIVOT_SHA, "pivot ledger hash")
    rows = g5.parse_pair_rows(g5.parse_hard_rows())
    by_short = {short_key(row): row for row in rows}
    require(len(by_short) == len(rows) == 14_806, "hard key collision")
    by_body = defaultdict(list)
    for row in rows:
        by_body[row["body"]].append(row)

    pivot = set()
    ledger_rows = 0
    for line in PIVOT_PATH.read_text().splitlines():
        if not line.startswith("H1;"):
            continue
        ledger_rows += 1
        fields = dict(part.split("=", 1) for part in line.split(";")[1:])
        entries = fields["pivot"].split(",")
        closed = all(entry.rsplit(":", 2)[-2] == "1" for entry in entries)
        if not closed:
            continue
        key = (
            tuple(map(int, fields["E"].split(","))),
            int(fields["rank"]),
            int(fields["a"]),
            tuple(map(int, fields["P"].split(","))),
        )
        require(key in by_short, f"unmatched pivot {key}")
        pivot.add(full_key(by_short[key]))
    require((ledger_rows, len(pivot)) == (11_842, 279), "pivot counts")

    g5_keys = {full_key(row) for row in rows if row["margin"] > 0}
    h1_keys = {
        full_key(row)
        for row in rows
        if (cutoff := h1_cutoff(row)) is not None and cutoff <= 15_000
    }
    exception_short = {
        (
            tuple(map(int, fields["E"].split(","))),
            int(fields["rank"]),
            int(fields["a"]),
        )
        for fields in h4.parse_exceptions()
    }
    exception_keys = {
        full_key(row)
        for row in rows
        if (row["body"], row["rank"], row["apex"]) in exception_short
    }
    require(len(exception_keys) == len(exception_short) == 52, "exception keys")

    thm2895, thm2898, thm2899, thm2901, thm2902 = g5.canonical_root_sets()
    g5.thm2903_controls()
    prior15 = thm2895 | thm2898 | thm2899 | thm2901
    one_hard = {
        body
        for body, body_rows in by_body.items()
        if len(body_rows) == 1 and body_rows[0]["direct_margin"] <= 0
    }
    through2905 = prior15 | thm2902 | one_hard
    g5_terminal = {
        body
        for body, body_rows in by_body.items()
        if all(full_key(row) in g5_keys for row in body_rows)
    } | set(thm2899)
    through2905 |= g5_terminal
    pivot_additions = {
        (1, 2, 3, 4, 5, 6, 11),
        (1, 2, 3, 4, 5, 8, 9),
        (3, 4, 6, 7, 9, 11, 12),
        (5, 6, 7, 8, 10, 11, 12),
        (5, 7, 9, 10, 11, 12, 14),
        (6, 7, 8, 9, 11, 12, 13),
    }
    proved88 = through2905 | pivot_additions
    require((len(through2905), len(proved88)) == (82, 88), "proved baseline")

    scalar_roots = set(thm2899)

    def roots(keys):
        return {
            body
            for body, body_rows in by_body.items()
            if all(full_key(row) in keys for row in body_rows)
        } | scalar_roots

    routes = {
        "G": g5_keys,
        "H": h1_keys,
        "P": pivot,
        "E": exception_keys,
    }
    membership = Counter(
        "".join(name for name, keys in routes.items() if key in keys)
        for key in set().union(*routes.values())
    )
    for names in ("G", "H", "P", "E", "GH", "GHP", "GHPE"):
        key_union = set().union(*(routes[name] for name in names))
        terminal = roots(key_union)
        union = proved88 | terminal
        print(
            f"{names}:branches={len(key_union)};roots={len(terminal)};"
            f"intersection={len(terminal & proved88)};additions={len(terminal - proved88)};"
            f"proved_union={len(union)};residual={3432-len(union)};"
            f"digest={digest_keys(union, names)}"
        )
        if names in ("GHP", "GHPE"):
            print(f"{names}_ADDITIONS={tuple(sorted(terminal - proved88))}")
    gh_roots = roots(g5_keys | h1_keys)
    ghp_roots = roots(g5_keys | h1_keys | pivot)
    ghpe_roots = roots(g5_keys | h1_keys | pivot | exception_keys)
    print(f"GHP_NEW_BEYOND_GH_AND_BASELINE={tuple(sorted(ghp_roots - gh_roots - proved88))}")
    print(f"GHPE_NEW_BEYOND_GHP_AND_BASELINE={tuple(sorted(ghpe_roots - ghp_roots - proved88))}")
    for body in sorted(ghp_roots - gh_roots - proved88):
        anatomy = []
        for row in by_body[body]:
            key = full_key(row)
            anatomy.append(
                (
                    row["rank"],
                    row["apex"],
                    row["prefix"],
                    "".join(name for name, keys in routes.items() if key in keys),
                    h1_cutoff(row),
                    row["margin"],
                )
            )
        print(f"GHP_MIXED={body};anatomy={tuple(anatomy)}")
    print(f"atomic_membership={tuple(sorted(membership.items()))}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
