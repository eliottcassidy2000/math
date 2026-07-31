#!/usr/bin/env python3
r"""Actual hostile-centre census behind the all-hard Hunter-star envelope.

For one THM-2905 row, write

    g(a)=a+sum_{r=2}^5 min(a,q_r,B_2-a).

If a five-cover exists and a is its largest singleton coverage, the Hunter
star proof gives g(a)>=h.  Let lambda be the least a in the legal centre
domain with g(a)>=h.  Every possible maximum singleton then belongs to

    H_1^star={w allowed:c_C(w)>=lambda}.

On every row not closed by G_5<h, this script proves lambda>h/7, seals
H_1^star by discrepancy, and computes the actual core exactly.  Ordering
the possible maxima then gives a pivot certificate on 4,071 centre
flags and six new whole roots.  This is not a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import ast
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
THM2895_OUT = (
    ROOT / "05-knowledge/results/lrc14_j6_suffix_parity_flag_closure_thm2895.out"
)
THM2898_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_k25_five_parity_matching_closure_thm2898.out"
)
THM2899_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.out"
)
THM2901_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_hard_global_pair_cap_census_codex_20260729.out"
)
THM2902_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_omission6_ranked_h1_depth1_closure_thm2902.out"
)
THM2903_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_one_hard_h3_link_core_census_codex_20260729.out"
)
THM2905_OUT = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_hard_hunter_star_envelope_census_codex_20260729.out"
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
THEOREM_OUTPUT_SHA256 = {
    THM2895_OUT: "c11260f6544a319e1cc1862c9221b188a4314860422470e465b82e7ce492b1b4",
    THM2898_OUT: "41f5e443f6d1ee2553c332da7709bd0c89f400b9ca154ddb6047f8ca724c6a40",
    THM2899_OUT: "dbd3dc5a8c44a55957a6e1ce660ca0e89fcd70e6c0d06d5ba47dc3a22f40c680",
    THM2901_OUT: "98b8ba171be1d38980e7271ef82e2bc1bde536afcf9864fa39138dbfbc93a3eb",
    THM2902_OUT: "cff46490d4a904947ec0fbe0cedfa59484c6b7974923656ba2459a55781192d7",
    THM2903_OUT: "5719083a83b275206907f47141fee8da2ba489194e31ba7c119f5313e3dfe73d",
    THM2905_OUT: "c346cbce451b4d0104707b071c9874798e2cadc853102038b229be9ad8a6afe4",
}

FIRST_EXTERNAL = 15
S2 = F(99, 70)

# Locked after ordered-pivot discovery; replay under ordinary/optimized Python.
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
    55_293,
    4_071,
    51_222,
    4,
    4,
    0,
    279,
    0,
    11_563,
    3_411,
    10,
    4,
    6,
    88,
    3_344,
    (
        (0, 279),
        (1, 737),
        (2, 1_428),
        (3, 1_822),
        (4, 2_033),
        (5, 2_183),
        (6, 1_629),
        (7, 969),
        (8, 446),
        (9, 206),
        (10, 82),
        (11, 20),
        (12, 7),
        (13, 1),
    ),
)
EXPECTED_ROOTS: tuple[tuple[int, ...], ...] | None = (
    (1, 2, 3, 4, 5, 6, 11),
    (1, 2, 3, 4, 5, 8, 9),
    (1, 2, 3, 5, 7, 10, 14),
    (1, 2, 4, 8, 11, 12, 13),
    (1, 2, 5, 6, 10, 12, 13),
    (2, 4, 6, 8, 11, 12, 13),
    (3, 4, 6, 7, 9, 11, 12),
    (5, 6, 7, 8, 10, 11, 12),
    (5, 7, 9, 10, 11, 12, 14),
    (6, 7, 8, 9, 11, 12, 13),
)
EXPECTED_ADDITIVE_ROOTS: tuple[tuple[int, ...], ...] | None = (
    (1, 2, 3, 4, 5, 6, 11),
    (1, 2, 3, 4, 5, 8, 9),
    (3, 4, 6, 7, 9, 11, 12),
    (5, 6, 7, 8, 10, 11, 12),
    (5, 7, 9, 10, 11, 12, 14),
    (6, 7, 8, 9, 11, 12, 13),
)
EXPECTED_EQUALITY_PIVOTS: tuple[tuple[object, ...], ...] | None = (
    ((1, 2, 3, 4, 7, 9, 14), 1, 24, 20),
    ((1, 2, 3, 4, 9, 12, 14), 4, 44, 15),
    ((1, 2, 6, 8, 9, 11, 12), 1, 20, 52),
    ((1, 3, 4, 6, 8, 9, 13), 2, 21, 20),
)
EXPECTED_OPEN_EQUALITY_PIVOTS: tuple[tuple[object, ...], ...] | None = ()
EXPECTED_LEDGER_SHA256: str | None = (
    "ec878244b922ba5f48633614a86a1f9706c1fbdd0ebd6c61f020291cfd737bab"
)
RAW_PATH_READ_BYTES = Path.read_bytes


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    """Hash repository text independently of LF/CRLF checkout policy."""

    payload = RAW_PATH_READ_BYTES(path)
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: unexpected lone carriage return",
    )
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


def lf_read_bytes(path: Path) -> bytes:
    """Read repository text on the canonical LF byte basis."""

    payload = RAW_PATH_READ_BYTES(path)
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: unexpected lone carriage return",
    )
    return payload.replace(b"\r\n", b"\n")


def unique_output_text(path: Path, prefix: str) -> str:
    matches = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, f"{path.name}: expected one {prefix!r} line")
    return matches[0]


def literal_output_value(path: Path, prefix: str) -> object:
    return ast.literal_eval(unique_output_text(path, prefix))


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
    # The imported residual/vector engines predate checkout-independent
    # text hashes.  Canonicalize their single-threaded import-time reads,
    # then retain the same normalized helper for any runtime hash checks.
    original_read_bytes = Path.read_bytes
    Path.read_bytes = lf_read_bytes
    try:
        spec.loader.exec_module(module)
    finally:
        Path.read_bytes = original_read_bytes
    if hasattr(module, "file_sha256"):
        module.file_sha256 = file_sha256
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
    pivot_rows: list[tuple[int, F, F, bool, bool]] = []
    for index, (coverage, label) in enumerate(core_rows):
        leaf_cap = row["pair_cap"] - coverage
        require(leaf_cap >= 0, "pair cap fell below pivot singleton")
        later_bounds = [
            min(later_coverage, leaf_cap)
            for later_coverage, _ in core_rows[index + 1 :]
        ]
        # Four formal noncore leaves suffice because a child has four slots.
        noncore_bound = min(row["threshold"], leaf_cap)
        leaf_bounds = later_bounds + [noncore_bound] * 4
        top_four_bound = sum(sorted(leaf_bounds, reverse=True)[:4], F(0))
        pivot_margin = mass - coverage - top_four_bound
        # When rho>=lambda, a noncore leaf has strict coverage <lambda.
        # At zero margin this closes exactly when no four later-core bounds
        # can attain the same top-four sum without a noncore contribution.
        later_top_four = (
            sum(sorted(later_bounds, reverse=True)[:4], F(0))
            if len(later_bounds) >= 4
            else None
        )
        equality_strict = (
            pivot_margin == 0
            and noncore_bound == row["threshold"]
            and (later_top_four is None or later_top_four < top_four_bound)
        )
        pivot_closed = pivot_margin > 0 or equality_strict
        pivot_rows.append(
            (label, top_four_bound, pivot_margin, pivot_closed, equality_strict)
        )
    return {
        **row,
        "actual_size": len(core_rows),
        "core_rows": core_rows,
        "pivot_rows": tuple(pivot_rows),
        "branch_closed": all(closed for _, _, _, closed, _ in pivot_rows),
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
        + ";pivot="
        + ",".join(
            f"{label}:{ftext(bound)}:{ftext(margin)}:{int(closed)}:{int(repaired)}"
            for label, bound, margin, closed, repaired in row["pivot_rows"]
        )
        + "\n"
    )


def current_proved_union() -> set[tuple[int, ...]]:
    for path, expected_hash in THEOREM_OUTPUT_SHA256.items():
        require(file_sha256(path) == expected_hash, f"{path.name} changed")

    thm2895 = {
        ast.literal_eval(line.removeprefix("BRANCH E=").split(";", 1)[0])
        for line in THM2895_OUT.read_text().splitlines()
        if line.startswith("BRANCH E=")
    }
    thm2898 = {
        ast.literal_eval(unique_output_text(THM2898_OUT, "root=").split(";", 1)[0])
    }
    thm2899 = set(literal_output_value(THM2899_OUT, "terminal_bodies="))
    thm2901 = set(literal_output_value(THM2901_OUT, "direct_terminal_bodies="))
    thm2902 = set(literal_output_value(THM2902_OUT, "closed_roots="))
    root_groups = (thm2895, thm2898, thm2899, thm2901, thm2902)
    require(
        tuple(map(len, root_groups)) == (4, 1, 5, 5, 2),
        "prior theorem root counts changed",
    )
    require(
        len(set().union(*root_groups)) == sum(map(len, root_groups)) == 17,
        "prior theorem root groups are no longer disjoint",
    )

    thm2903_controls = {
        "closed_roots": literal_output_value(THM2903_OUT, "closed_roots="),
        "overlap": literal_output_value(THM2903_OUT, "overlap_THM2902="),
        "proved_union": literal_output_value(THM2903_OUT, "proved_union="),
        "official_residual": literal_output_value(
            THM2903_OUT, "official_residual="
        ),
        "root_digest": unique_output_text(THM2903_OUT, "root61_digest="),
        "proof_digest": unique_output_text(THM2903_OUT, "proof_digest="),
        "mode": unique_output_text(THM2903_OUT, "mode="),
        "sentinel": unique_output_text(THM2903_OUT, "all_exact_controls="),
    }
    require(
        thm2903_controls
        == {
            "closed_roots": 61,
            "overlap": (
                (1, 2, 3, 4, 5, 10, 12),
                (1, 2, 3, 4, 5, 12, 13),
            ),
            "proved_union": 76,
            "official_residual": 3_356,
            "root_digest": (
                "f58c4143f329d215ff9bb7ec594d172e831fbccbab713a2398dc6cb53c60b8b7"
            ),
            "proof_digest": (
                "08dc4a539544a417ff884ef4631b10d6eb14fd63d7b62b06e76d92e3d4d9b162"
            ),
            "mode": "LOCKED ordinary_and_optimized_replays_agree",
            "sentinel": "PASS",
        },
        "THM2903 proof controls changed",
    )

    thm2905_counts = literal_output_value(THM2905_OUT, "counts=")
    require(
        thm2905_counts[:17]
        == (
            14_806,
            2_964,
            1_835,
            1_129,
            0,
            52,
            0,
            3_427,
            16,
            5,
            11,
            61,
            5,
            76,
            6,
            82,
            3_350,
        ),
        "THM2905 root-composition counts changed",
    )
    star_roots = set(literal_output_value(THM2905_OUT, "star_roots="))
    thm2905_additive = set(
        literal_output_value(THM2905_OUT, "additive_roots=")
    )
    require(
        (
            len(star_roots),
            len(thm2905_additive),
            unique_output_text(THM2905_OUT, "row_digest="),
            unique_output_text(THM2905_OUT, "mode="),
            unique_output_text(THM2905_OUT, "all_exact_controls="),
        )
        == (
            16,
            6,
            "c1f60bdbc3669202dde60d8f44c9db887807b179cb42e02137d475c0f282d066",
            "LOCKED",
            "PASS",
        ),
        "THM2905 proof controls changed",
    )

    by_body: dict[tuple[int, ...], list[F]] = {}
    for line in PAIR_LEDGER.read_text().splitlines():
        if not line.startswith("PAIR;"):
            continue
        row = fields(line)
        body = tuple(map(int, row["E"].split(",")))
        by_body.setdefault(body, []).append(parse_fraction(row["mdirect"]))
    one_hard = {
        body
        for body, margins in by_body.items()
        if len(margins) == 1 and margins[0] <= 0
    }
    prior_fifteen = thm2895 | thm2898 | thm2899 | thm2901
    require(len(prior_fifteen) == 15, "prior fifteen-root union changed")
    require(
        len(one_hard) == thm2903_controls["closed_roots"],
        "THM2903 one-hard bank changed",
    )
    require(
        one_hard & thm2902 == set(thm2903_controls["overlap"]),
        "THM2903/THM2902 overlap changed",
    )
    through_2903 = prior_fifteen | thm2902 | one_hard
    require(
        len(through_2903) == thm2903_controls["proved_union"],
        "proved union through THM2903 changed",
    )
    require(
        3_432 - len(through_2903) == thm2903_controls["official_residual"],
        "THM2903 residual identity changed",
    )
    require(
        star_roots - through_2903 == thm2905_additive,
        "THM2905 additive-root difference changed",
    )
    through_2905 = through_2903 | thm2905_additive
    require(
        len(through_2905) == thm2905_counts[15]
        and 3_432 - len(through_2905) == thm2905_counts[16],
        "proved union through THM2905 changed",
    )
    return through_2905


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
    pivots = [
        (row, label, bound, margin, closed, repaired)
        for row in rows
        for label, bound, margin, closed, repaired in row["pivot_rows"]
    ]
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = {}
    for row in rows:
        by_body.setdefault(row["body"], []).append(row)
    closed_roots = tuple(
        sorted(
            body
            for body, body_rows in by_body.items()
            if all(row["branch_closed"] for row in body_rows)
        )
    )
    proved_union = current_proved_union()
    additive_roots = tuple(sorted(set(closed_roots) - proved_union))
    new_union = proved_union | set(closed_roots)
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
        len(pivots),
        sum(closed for _, _, _, _, closed, _ in pivots),
        sum(not closed for _, _, _, _, closed, _ in pivots),
        sum(margin == 0 for _, _, _, margin, _, _ in pivots),
        sum(repaired for _, _, _, _, _, repaired in pivots),
        sum(
            margin == 0 and not closed
            for _, _, _, margin, closed, _ in pivots
        ),
        sum(row["branch_closed"] for row in rows),
        sum(row["pair_exception"] and row["branch_closed"] for row in rows),
        sum(not row["branch_closed"] for row in rows),
        len(by_body),
        len(closed_roots),
        len(set(closed_roots) & proved_union),
        len(additive_roots),
        len(new_union),
        3_432 - len(new_union),
        tuple(
            sorted(
                Counter(
                    sum(not closed for _, _, _, closed, _ in row["pivot_rows"])
                    for row in rows
                ).items()
            )
        ),
    )
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "ranked-H1 counts changed")
    if EXPECTED_ROOTS is not None:
        require(closed_roots == EXPECTED_ROOTS, "ordered-H1 root list changed")
    if EXPECTED_ADDITIVE_ROOTS is not None:
        require(
            additive_roots == EXPECTED_ADDITIVE_ROOTS,
            "additive ordered-H1 roots changed",
        )
    equality_pivots = tuple(
        (
            row["body"],
            row["rank"],
            row["apex"],
            label,
        )
        for row, label, _, margin, _, _ in pivots
        if margin == 0
    )
    if EXPECTED_EQUALITY_PIVOTS is not None:
        require(
            equality_pivots == EXPECTED_EQUALITY_PIVOTS,
            "ordered-H1 equality boundary changed",
        )
    open_equality_pivots = tuple(
        (
            row["body"],
            row["rank"],
            row["apex"],
            label,
        )
        for row, label, _, margin, closed, _ in pivots
        if margin == 0 and not closed
    )
    if EXPECTED_OPEN_EQUALITY_PIVOTS is not None:
        require(
            open_equality_pivots == EXPECTED_OPEN_EQUALITY_PIVOTS,
            "ordered-H1 open equality boundary changed",
        )

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/all-hard/ranked-H1-Hunter-pivots/v2\n")
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
    print(f"closed_roots={closed_roots}")
    print(f"additive_roots={additive_roots}")
    print(f"strict_equality_pivots={equality_pivots}")
    print(f"open_equality_pivots={open_equality_pivots}")
    positive_pivots = [pivot for pivot in pivots if pivot[3] > 0]
    failed_pivots = [pivot for pivot in pivots if not pivot[4]]
    closest_positive = min(
        positive_pivots,
        key=lambda pivot: (
            pivot[3],
            pivot[0]["body"],
            pivot[0]["rank"],
            pivot[1],
        ),
    )
    closest_failure = max(
        failed_pivots,
        key=lambda pivot: (
            pivot[3],
            tuple(-value for value in pivot[0]["body"]),
            -pivot[0]["rank"],
            -pivot[1],
        ),
    )
    print(
        "closest_positive_pivot="
        f"{closest_positive[0]['body']};rank={closest_positive[0]['rank']};"
        f"a={closest_positive[0]['apex']};x={closest_positive[1]};"
        f"margin={ftext(closest_positive[3])}"
    )
    print(
        "closest_failed_pivot="
        f"{closest_failure[0]['body']};rank={closest_failure[0]['rank']};"
        f"a={closest_failure[0]['apex']};x={closest_failure[1]};"
        f"margin={ftext(closest_failure[3])}"
    )
    print(f"ledger_sha256={ledger_sha256}")
    print(
        "mode=DISCOVERY"
        if any(
            value is None
            for value in (
                EXPECTED_COUNTS,
                EXPECTED_ROOTS,
                EXPECTED_ADDITIVE_ROOTS,
                EXPECTED_EQUALITY_PIVOTS,
                EXPECTED_OPEN_EQUALITY_PIVOTS,
                EXPECTED_LEDGER_SHA256,
            )
        )
        else "mode=LOCKED"
    )
    print(
        "scope=11842 G5-surviving scalar-hard rows;actual hostile-centre "
        "cores and ordered-pivot upper bounds;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
