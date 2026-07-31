#!/usr/bin/env python3
"""Exact H4 membership census on THM-2901's 52 pair-cap exceptions.

For a scalar-hard carrier C of mass h, THM-2899 gives the attained global
singleton maximum q1<3h/7.  Hence every hypothetical five-cover contains
at least two labels in

    H4(C)={w allowed: c_C(w)>=(h-q1)/4}.

The discrepancy bound seals this core globally.  This script reconstructs
the 52 exact-B2 exceptions in THM-2901, computes actual H4 membership rather
than cutoff-universe binomial counts, and records the unordered pair-flag
workload for the heavy-link child census.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENGINE = (
    ROOT
    / "04-computation/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.py"
)
ENGINE_SHA256 = (
    "7ba8244d8fc78ebc0d9381e05d69ca53c849d6008ff9cfb43f0efcbb4b394f81"
)
LEDGER = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.ledger.out"
)
LEDGER_SHA256 = (
    "5dea0eaa45dd52fbf1bef7cfcc328899a4789bc277b6e1e8ac2f4bdf192b85e4"
)

EXPECTED_COUNTS: tuple[int, ...] | None = (
    52,
    51,
    1,
    52,
    12,
    44,
    18_290,
    2_026,
    1_348,
)
EXPECTED_SIZE_DISTRIBUTION: tuple[tuple[int, int], ...] | None = (
    (12, 1),
    (14, 1),
    (16, 1),
    (17, 2),
    (18, 4),
    (20, 2),
    (21, 3),
    (22, 5),
    (23, 3),
    (24, 4),
    (25, 3),
    (26, 3),
    (27, 2),
    (28, 4),
    (29, 2),
    (30, 3),
    (32, 1),
    (33, 1),
    (39, 1),
    (40, 1),
    (41, 2),
    (42, 1),
    (43, 1),
    (44, 1),
)
EXPECTED_DIGEST: str | None = (
    "f75c1440a58c8724cb6bcebd41343f81daa6d3bca302a6122ba3721373479cb6"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_engine():
    require(file_sha256(ENGINE) == ENGINE_SHA256, "THM-2901 engine changed")
    spec = importlib.util.spec_from_file_location("h4_exception_engine", ENGINE)
    require(spec is not None and spec.loader is not None, "cannot load engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E = load_engine()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_ints(text: str) -> tuple[int, ...]:
    return tuple(map(int, text.split(","))) if text else ()


def parse_exceptions() -> list[dict[str, str]]:
    require(file_sha256(LEDGER) == LEDGER_SHA256, "THM-2901 ledger changed")
    rows: list[dict[str, str]] = []
    for line in LEDGER.read_text().splitlines():
        if not line.startswith("PAIR;"):
            continue
        fields = dict(
            field.split("=", 1)
            for field in line.removeprefix("PAIR;").split(";")
        )
        if F(fields["mB2"]) <= 0:
            rows.append(fields)
    require(len(rows) == 52, "pair-cap exception count changed")
    return rows


def exact_row(fields: dict[str, str]) -> dict[str, object]:
    body = parse_ints(fields["E"])
    apex = int(fields["a"])
    prefix = parse_ints(fields["P"])
    mass = F(fields["h"])
    components = int(fields["r"])
    q1 = F(fields["q1"])
    carrier, reconstructed_components, reconstructed_mass = E.R.CORE.good_norm(
        tuple(sorted((*body, apex)))
    )
    require(
        reconstructed_components == components and reconstructed_mass == mass,
        "literal carrier reconstruction changed",
    )
    level = (mass - q1) / 4
    delta = level - mass / 7
    require(delta > 0 and q1 < 3 * mass / 7, "H4 core is not finite")
    gamma = E.S2 * components / 7
    ratio = gamma / delta
    cutoff = (ratio.numerator - 1) // ratio.denominator
    require(
        mass / 7 + gamma / (cutoff + 1) <= level,
        "H4 discrepancy tail did not seal",
    )
    forbidden = frozenset(prefix)
    labels = [
        label
        for label in range(E.FIRST_EXTERNAL, cutoff + 1)
        if label not in forbidden
    ]
    coverages = E.T.coverages_many(carrier, labels)
    ranked = sorted(coverages, key=lambda item: (-item[0], item[1]))
    require(ranked[0][0] == q1, "ledger q1/reconstructed maximum changed")
    core = tuple(label for value, label in ranked if value >= level)
    require(len(core) >= 2, "H4 core has fewer than two labels")
    for label in dict.fromkeys((core[0], core[-1], labels[-1])):
        vector_value = next(value for value, speed in coverages if speed == label)
        require(
            vector_value == E.T.coverage(carrier, label),
            f"scalar/vector mismatch at {label}",
        )
    return {
        "stratum": fields["S"],
        "body": body,
        "rank": int(fields["rank"]),
        "apex": apex,
        "prefix": prefix,
        "mass": mass,
        "components": components,
        "q1": q1,
        "level": level,
        "cutoff": cutoff,
        "core": core,
        "flags": comb(len(core), 2),
    }


def row_line(row: dict[str, object]) -> str:
    return (
        f"H4;S={row['stratum']};E={','.join(map(str, row['body']))};"
        f"rank={row['rank']};a={row['apex']};"
        f"P={','.join(map(str, row['prefix']))};"
        f"h={ftext(row['mass'])};r={row['components']};"
        f"q1={ftext(row['q1'])};level={ftext(row['level'])};"
        f"cutoff={row['cutoff']};size={len(row['core'])};"
        f"flags={row['flags']};core={','.join(map(str, row['core']))}\n"
    )


def main() -> None:
    rows = [exact_row(fields) for fields in parse_exceptions()]
    sizes = [len(row["core"]) for row in rows]
    distribution = tuple(
        (size, sizes.count(size))
        for size in sorted(set(sizes))
    )
    counts = (
        len(rows),
        sum(row["rank"] == 1 for row in rows),
        sum(row["rank"] == 2 for row in rows),
        len({row["body"] for row in rows}),
        min(sizes),
        max(sizes),
        sum(row["flags"] for row in rows),
        max(row["cutoff"] for row in rows),
        sum(sizes),
    )
    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/paircap-exception/H4-membership/v1\n")
    for row in rows:
        digest.update(row_line(row).encode())
    semantic_digest = digest.hexdigest()
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "aggregate counts changed")
    if EXPECTED_SIZE_DISTRIBUTION is not None:
        require(
            distribution == EXPECTED_SIZE_DISTRIBUTION,
            "H4 size distribution changed",
        )
    if EXPECTED_DIGEST is not None:
        require(semantic_digest == EXPECTED_DIGEST, "row digest changed")

    print("LRC14 j6 pair-cap-exception exact H4 membership census")
    print(f"counts={counts}")
    print(f"size_distribution={distribution}")
    for row in rows:
        print(row_line(row).rstrip())
    print(f"semantic_digest={semantic_digest}")
    print(
        "mode="
        + ("DISCOVERY" if EXPECTED_DIGEST is None else "LOCKED")
    )
    print(
        "scope=52 exact-B2 exceptions;actual H4 membership and pair flags;"
        "no child closures;not LRC14"
    )


if __name__ == "__main__":
    main()
