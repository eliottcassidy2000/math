#!/usr/bin/env python3
r"""Actual hostile-centre census behind the all-hard Hunter-star envelope.

For one THM-2905 row, write

    g(a)=a+sum_{r=2}^5 min(a,q_r,B_2-a).

If a five-cover exists and a is its largest singleton coverage, the Hunter
star proof gives g(a)>=h.  Let lambda be the least a in the legal centre
domain with g(a)>=h.  Every possible maximum singleton then belongs to

    H_1^star={w allowed:c_C(w)>=lambda}.

On every row surviving G_5<h, this script proves lambda>h/7, seals H_1^star
by discrepancy, and computes the actual core exactly.  The result is an
ordered-pivot workload and proof sidecar, not a closure and not LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HARD_LEDGER = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_root_ranked_suffix_scalar_hard_ledger_codex_20260729.out"
)
PAIR_LEDGER = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.ledger.out"
)
RESIDUAL_PATH = (
    ROOT
    / "04-computation/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.py"
)
VECTOR_PATH = (
    ROOT
    / "04-computation/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.py"
)
HARD_LEDGER_SHA256 = (
    "6be9a6c9218f3b42b2eea733c9050f5d35160664af0f19390337b3c5be57cb37"
)
PAIR_LEDGER_SHA256 = (
    "5dea0eaa45dd52fbf1bef7cfcc328899a4789bc277b6e1e8ac2f4bdf192b85e4"
)
RESIDUAL_SHA256 = (
    "a5f3dcc1a23defea4b3dc067675d83141f1866022d6d01946617a97de69e5b0e"
)
VECTOR_SHA256 = (
    "dff97f67b1104c25589802a6a2f216b6e7bfedd58eebfa1bcce615d59c1e872f"
)

FIRST_EXTERNAL = 15
S2 = F(99, 70)

# Locked after discovery; replay under ordinary and optimized Python.
EXPECTED_COUNTS: tuple[object, ...] | None = (
    11_842,
    52,
    4_797_677,
    149,
    1_013,
    55_293,
    1,
    13,
    0,
    (
        (0, 1),
        (25, 3),
        (50, 5),
        (75, 6),
        (90, 7),
        (95, 8),
        (99, 9),
        (100, 13),
    ),
    (
        (5, 2_566),
        (4, 2_488),
        (3, 1_906),
        (6, 1_782),
        (2, 1_028),
        (7, 1_016),
        (8, 454),
        (1, 283),
        (9, 209),
        (10, 82),
        (11, 20),
        (12, 7),
        (13, 1),
    ),
)
EXPECTED_LEDGER_SHA256: str | None = (
    "05b2ce123e995b13f4ef93c948e5dff2ff190ef50dd0be4e893ffe452d25a791"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_fraction(text: str) -> F:
    numerator, denominator = text.split("/")
    return F(int(numerator), int(denominator))


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def fields(line: str) -> dict[str, str]:
    return {
        item.split("=", 1)[0]: item.split("=", 1)[1]
        for item in line.split(";")[1:]
        if "=" in item
    }


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(HARD_LEDGER) == HARD_LEDGER_SHA256, "hard ledger changed")
require(file_sha256(PAIR_LEDGER) == PAIR_LEDGER_SHA256, "pair ledger changed")
require(file_sha256(RESIDUAL_PATH) == RESIDUAL_SHA256, "residual engine changed")
require(file_sha256(VECTOR_PATH) == VECTOR_SHA256, "vector engine changed")
R = load_module("lrc14_ranked_h1_residual", RESIDUAL_PATH)
V = load_module("lrc14_ranked_h1_vector", VECTOR_PATH)


def parse_hard_rows() -> dict[tuple[object, ...], dict[str, object]]:
    rows: dict[tuple[object, ...], dict[str, object]] = {}
    for line in HARD_LEDGER.read_text().splitlines():
        if not line.startswith("HARD;"):
            continue
        row = fields(line)
        body = tuple(map(int, row["E"].split(",")))
        prefix = tuple(map(int, row["prefix"].split(",")))
        key = (body, int(row["rank"]), int(row["apex"]), prefix)
        top5 = tuple(
            (int(item.split(":", 1)[0]), parse_fraction(item.split(":", 1)[1]))
            for item in row["top5"].split(",")
        )
        require(len(top5) == 5, "hard row lost top-five data")
        require(
            all(top5[index][1] >= top5[index + 1][1] for index in range(4)),
            "hard singleton ranks are not ordered",
        )
        require(key not in rows, "duplicate hard row")
        rows[key] = {
            "body": body,
            "rank": key[1],
            "apex": key[2],
            "prefix": prefix,
            "stratum": row["S"],
            "gate_size": int(row["K"]),
            "mass": parse_fraction(row["m"]),
            "components": int(row["r"]),
            "top5": top5,
        }
    require(len(rows) == 14_806, "hard row count changed")
    return rows


def clipped(value: F, upper: F) -> F:
    return max(F(0), min(upper, value))


def star_data(qs: tuple[F, ...], pair_cap: F, mass: F) -> tuple[F, F]:
    upper = min(qs[0], pair_cap)
    breakpoints = {F(0), upper, clipped(pair_cap / 2, upper)}
    for singleton in qs[1:]:
        breakpoints.add(clipped(singleton, upper))
        breakpoints.add(clipped(pair_cap - singleton, upper))
    breakpoints = sorted(breakpoints)

    def objective(center: F) -> F:
        return center + sum(
            (min(center, singleton, pair_cap - center) for singleton in qs[1:]),
            F(0),
        )

    envelope = max(objective(center) for center in breakpoints)
    if envelope < mass:
        return envelope, F(-1)
    for left, right in zip(breakpoints, breakpoints[1:]):
        left_value = objective(left)
        right_value = objective(right)
        midpoint = (left + right) / 2
        require(
            2 * objective(midpoint) == left_value + right_value,
            "star breakpoint list is incomplete",
        )
        if left_value >= mass:
            return envelope, left
        if left_value < mass <= right_value:
            threshold = (
                left
                + (mass - left_value) * (right - left) / (right_value - left_value)
            )
            require(objective(threshold) == mass, "hostile threshold solve failed")
            return envelope, threshold
    require(objective(breakpoints[-1]) >= mass, "hostile star set disappeared")
    return envelope, breakpoints[-1]


def survivor_inputs() -> list[dict[str, object]]:
    hard_rows = parse_hard_rows()
    inputs: list[dict[str, object]] = []
    seen: set[tuple[object, ...]] = set()
    for line in PAIR_LEDGER.read_text().splitlines():
        if not line.startswith("PAIR;"):
            continue
        row = fields(line)
        body = tuple(map(int, row["E"].split(",")))
        prefix = tuple(map(int, row["P"].split(",")))
        key = (body, int(row["rank"]), int(row["a"]), prefix)
        require(key in hard_rows and key not in seen, "hard/pair join changed")
        seen.add(key)
        hard = hard_rows[key]
        require(row["S"] == hard["stratum"], "joined stratum changed")
        require(int(row["K"]) == hard["gate_size"], "joined gate size changed")
        require(int(row["r"]) == hard["components"], "joined component count changed")
        mass = parse_fraction(row["h"])
        require(mass == hard["mass"], "joined mass changed")
        qs = tuple(value for _, value in hard["top5"])
        require(
            (
                parse_fraction(row["q1"]),
                parse_fraction(row["q2"]),
                parse_fraction(row["q3"]),
                parse_fraction(row["q5"]),
            )
            == (qs[0], qs[1], qs[2], qs[4]),
            "joined singleton ranks changed",
        )
        pair_cap = parse_fraction(row["B2"])
        envelope, threshold = star_data(qs, pair_cap, mass)
        if envelope < mass:
            continue
        require(threshold >= 0, "star survivor lost hostile threshold")
        delta = threshold - mass / 7
        require(delta > 0, "hostile centre threshold is not discrepancy-finite")
        gamma = S2 * hard["components"] / 7
        cutoff = ceiling(gamma / delta) - 1
        require(
            mass / 7 + gamma / (cutoff + 1) <= threshold,
            "hostile-centre cutoff failed",
        )
        inputs.append(
            {
                **hard,
                "qs": qs,
                "pair_cap": pair_cap,
                "pair_exception": parse_fraction(row["mB2"]) <= 0,
                "star_envelope": envelope,
                "threshold": threshold,
                "delta": delta,
                "cutoff": cutoff,
            }
        )
    require(len(seen) == len(hard_rows), "pair ledger lost hard rows")
    require(len(inputs) == 11_842, "G5 survivor universe changed")
    return sorted(inputs, key=lambda row: (row["body"], row["rank"], row["apex"]))


def exact_coverages(
    carrier: list[tuple[F, F]], labels: list[int]
) -> list[tuple[F, int]]:
    rows = V.coverages_many(carrier, labels)
    require(len(rows) == len(labels), "vector coverage length changed")
    if rows:
        controls = tuple(dict.fromkeys((labels[0], labels[-1], labels[len(labels) // 2])))
        by_label = {label: value for value, label in rows}
        for label in controls:
            require(
                by_label[label] == R.coverage(carrier, label),
                f"scalar/vector mismatch at {label}",
            )
    return rows


def actual_core(row: dict[str, object]) -> dict[str, object]:
    carrier, components, mass = R.CORE.good_norm(
        tuple(sorted((*row["body"], row["apex"])))
    )
    require(
        components == row["components"] and mass == row["mass"] and mass > 0,
        "literal carrier reconstruction changed",
    )
    forbidden = frozenset(row["prefix"])
    labels = [
        label
        for label in range(FIRST_EXTERNAL, row["cutoff"] + 1)
        if label not in forbidden
    ]
    coverages = exact_coverages(carrier, labels)
    core_rows = tuple(
        sorted(
            (
                (value, label)
                for value, label in coverages
                if value >= row["threshold"]
            ),
            key=lambda item: (-item[0], item[1]),
        )
    )
    require(
        all(label not in forbidden for _, label in core_rows),
        "forbidden label entered hostile centre core",
    )
    return {
        **row,
        "actual_size": len(core_rows),
        "core_rows": core_rows,
        "scanned": len(labels),
    }


def ledger_line(row: dict[str, object]) -> str:
    return (
        f"E={','.join(map(str, row['body']))};S={row['stratum']};"
        f"rank={row['rank']};a={row['apex']};"
        f"P={','.join(map(str, row['prefix']))};"
        f"h={ftext(row['mass'])};B2={ftext(row['pair_cap'])};"
        f"G5={ftext(row['star_envelope'])};lambda={ftext(row['threshold'])};"
        f"delta={ftext(row['delta'])};N={row['cutoff']};scan={row['scanned']};"
        f"exception={int(row['pair_exception'])};"
        + "H1="
        + ",".join(f"{label}:{ftext(value)}" for value, label in row["core_rows"])
        + "\n"
    )


def nearest_rank(values: list[int]) -> tuple[tuple[int, int], ...]:
    require(values, "empty quantile population")
    ordered = sorted(values)
    return tuple(
        (
            p,
            ordered[0 if p == 0 else (p * len(ordered) + 99) // 100 - 1],
        )
        for p in (0, 25, 50, 75, 90, 95, 99, 100)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument("--ledger", type=Path)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    inputs = survivor_inputs()
    context = mp.get_context("spawn")
    if args.workers == 1:
        rows = list(map(actual_core, inputs))
    else:
        with context.Pool(args.workers) as pool:
            rows = pool.map(actual_core, inputs)
    rows.sort(key=lambda row: (row["body"], row["rank"], row["apex"]))

    sizes = [row["actual_size"] for row in rows]
    counts = (
        len(rows),
        sum(row["pair_exception"] for row in rows),
        sum(row["scanned"] for row in rows),
        min(row["cutoff"] for row in rows),
        max(row["cutoff"] for row in rows),
        sum(sizes),
        min(sizes),
        max(sizes),
        sum(size == 0 for size in sizes),
        nearest_rank(sizes),
        tuple(Counter(size for size in sizes).most_common(20)),
    )
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "ranked-H1 counts changed")

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/all-hard/ranked-H1-Hunter-pivots/v1\n")
    for row in rows:
        digest.update(ledger_line(row).encode())
    ledger_sha256 = digest.hexdigest()
    if EXPECTED_LEDGER_SHA256 is not None:
        require(ledger_sha256 == EXPECTED_LEDGER_SHA256, "ranked-H1 ledger changed")
    if args.ledger is not None:
        args.ledger.write_text(
            "LRC14 j6 all-hard ranked H1 Hunter-pivot ledger\n"
            + "".join("H1;" + ledger_line(row) for row in rows)
            + f"ledger_sha256={ledger_sha256}\n"
            + "scope=11842 G5-surviving hard rows;actual hostile centre cores;not LRC14\n"
        )

    print("LRC14 all-hard ranked H1 Hunter-pivot census")
    print(f"counts={counts}")
    print(f"cutoff_quantiles={nearest_rank([row['cutoff'] for row in rows])}")
    print(f"size_quantiles={nearest_rank(sizes)}")
    print(f"ledger_sha256={ledger_sha256}")
    print(
        "mode=DISCOVERY"
        if EXPECTED_COUNTS is None or EXPECTED_LEDGER_SHA256 is None
        else "mode=LOCKED"
    )
    print(
        "scope=11842 G5-surviving scalar-hard rows;actual hostile-centre "
        "cores;ordered-pivot workload only;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
