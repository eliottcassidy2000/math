#!/usr/bin/env python3
r"""Exact sharp-discrepancy extension of the THM-2911 ranked-H1 census.

THM-2911 used the strict THM-735 interval tail

    c_C(w) < |C|/7 + (99/70) r(C)/(7w).

Section 2 of THM-1094 proves the sharper one-interval estimate

    |I intersect D_w| <= |I|/7 + 6/(49w).

Summing over the ``r(C)`` interval components (endpoints have measure zero)
therefore gives the globally valid non-strict tail

    c_C(w) <= |C|/7 + 6 r(C)/(49w).

This program changes only that tail constant, from ``S2=99/70`` to
``S2=6/7`` in the existing THM-2911 exact scout.  It:

1. reconstructs all 14,806 hard branch rows;
2. identifies exactly the branches newly brought under the fixed 15,000
   cutoff;
3. reruns the literal ranked-H1 ordered-child certificate on those rows;
4. composes the certified keys with the already-proved G5, THM-2904 pivot,
   and THM-2907 H4 exception keys; and
5. takes exact set difference against both the historical 351-root union
   through THM-2913 and the live 1,041-root union through THM-2920.

No sampled residue or floating-point comparison is used.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import multiprocessing as mp
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM1094 = (
    ROOT
    / "01-canon/theorems/THM-1094-exact-two-comb-component-theorem.md"
)
SCOUT_PATH = ROOT / "04-computation/lrc14_j6_ranked_h1_thm2911/scout.py"
VERIFY_PATH = ROOT / "04-computation/lrc14_j6_ranked_h1_thm2911/verify.py"
HUNTER_PATH = (
    ROOT
    / "04-computation/lrc14_j6_ranked_h1_thm2911/"
    "hunter_two_star_exceptions.py"
)
THM2911_OUT = (
    ROOT / "05-knowledge/results/lrc14_j6_ranked_h1_thm2911/locked.out"
)
THM2912_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_one_h3_row_child_top4_census_codex_20260729.out"
)
THM2913_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_one_h3_row_pair_hunter_toothpick_closure_codex_20260729.out"
)
THM2916_SOURCE = (
    ROOT
    / "04-computation/"
    "lrc14_j6_two_h3_row_child_top4_census_codex_20260729.py"
)
THM2916_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_two_h3_row_child_top4_census_codex_20260729.out"
)
THM2920_SOURCE = (
    ROOT
    / "04-computation/"
    "lrc14_j6_two_h3_row_pair_hunter_recursive_toothpick_closure_codex_20260729.py"
)
THM2920_OUT = (
    ROOT
    / "05-knowledge/results/"
    "lrc14_j6_two_h3_row_pair_hunter_recursive_toothpick_closure_codex_20260729.out"
)

THM1094_SHA256 = (
    "b59393bdd726c88bd3692fefb9bde1c52d0760e7bc438878d7fd9157d3ed29c2"
)
SCOUT_SHA256 = (
    "6abc06972cda64bb4f53db26cca62bbea5362ca178bdbf5bf7398e8b0f28317a"
)
VERIFY_SHA256 = (
    "e0ac67539f7ff09376645a62beef0a9ac7d0766a2e749666f94d1fd4d6487b15"
)
HUNTER_SHA256 = (
    "c452781c7cf6a2ada6be8984b9c9cbe7aab8369c5a9333721ef8ecc8e0207393"
)
THM2911_OUT_SHA256 = (
    "e5c58cc2eb325928c00839c2593450ea7cce8945b3835898ec83c6c5f42fac9b"
)
THM2912_OUT_SHA256 = (
    "454d87c8beeb81405b031cce4b40bdda0d385cfcd9c48e6fcf4eb810cfc00c5a"
)
THM2913_OUT_SHA256 = (
    "3604644a9691b13e7fa245249b68c9586ec2775996834f7761f32eb0b89f3e64"
)
THM2916_SOURCE_SHA256 = (
    "7d3c4e82abb0ce3af13c43c1e03f09d4be1ee3dbfd496ca05273672cd2a462a6"
)
THM2916_OUT_SHA256 = (
    "6b31abaadd4f089a9f98a9eea49845c0ed0123a810bb4bf4f3c0309c519e7df9"
)
THM2920_SOURCE_SHA256 = (
    "049be1e331fb0c4fc46e703ffb37d61be2b3ec3b781835f480e55ba89bdf894e"
)
THM2920_OUT_SHA256 = (
    "1a38fd441dfd77a4f5d30d45d3160febc33d2d4eeb6247b223f10a1e31a8aefb"
)

OLD_S2 = F(99, 70)
SHARP_S2 = F(6, 7)
MAX_CUTOFF = 15_000
MAX_COMBINATIONS = 50_000_000_000
EXTENSION_BODY = (1, 3, 5, 7, 10, 11, 14)
EXTENSION_CUTOFF = 18_869
HUNTER_BODY = (3, 6, 7, 9, 11, 12, 13)
HUNTER_RANK = 3
HUNTER_PIVOT = 25

# Locked after an independent ordinary replay.
EXPECTED_COUNTS: tuple[object, ...] | None = (
    14_806,
    5_999,
    6_180,
    6_079,
    80,
    79,
    1,
    80,
    81,
    3,
    15_028,
    24_701,
    9_107,
    14_970,
    315,
    755,
    525,
    1,
    F(167, 1_891_890),
    2_964,
    279,
    52,
    6_170,
    6_240,
    6_241,
    6_331,
    138,
    139,
    140,
    146,
    2,
    351,
    2,
    1_041,
    1,
    1_042,
    2_390,
    101,
    7,
    6,
    2,
    1,
    99,
)
EXPECTED_TARGET_KEY_DIGEST: str | None = (
    "a6fe73463c89a889137e0cf25c6d1badbb792654862676c3186c5c9adc12018f"
)
EXPECTED_CLOSED_KEY_DIGEST: str | None = (
    "a6fe73463c89a889137e0cf25c6d1badbb792654862676c3186c5c9adc12018f"
)
EXPECTED_CERTIFIED_KEY_DIGEST: str | None = (
    "3f6280598288370f5acd7ce51e7439faa730edc0af46980a65da6095d8c25b57"
)
EXPECTED_LEDGER_SHA256: str | None = (
    "6de74387474e14bbfd622c5371d2d246588532e29764938f6210277478950831"
)
EXPECTED_THROUGH_2913_ROOT_DIGEST = (
    "34c21b0a217aaf6bc73fa80f1076c5ca58d96cd287d49f461c0fc6d19814b451"
)
EXPECTED_CURRENT_ROOT_DIGEST = (
    "ad3f17c4f8732831ff4810c4fedeb40ef299468db35e83d029301ed7d9eb69f2"
)
EXPECTED_FIXED_ROUTE_ROOT_DIGEST = (
    "b474b90ece71b49545794e971feae11e9e923bac679aa518fe8bd1d149cec92e"
)
EXPECTED_SHARP_ROUTE_ROOT_DIGEST: str | None = (
    "e27fa211742afd735facab5e6e3abfbc1bb9bd8b130171cdee1148e22c0e72ab"
)
EXPECTED_SATURATION_ROOT_DIGEST: str | None = (
    "fb3c52072666c4860032ebef2c50a561dfe3dcfbdc66bfb0325b6f27370b1a41"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    """Hash repository text independently of LF/CRLF checkout policy."""

    payload = path.read_bytes()
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: unexpected lone carriage return",
    )
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def literal_output_value(path: Path, prefix: str) -> object:
    matches = [
        line.removeprefix(prefix)
        for line in path.read_text().splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, f"{path.name}: expected one {prefix!r} line")
    return ast.literal_eval(matches[0])


def ftext(value: F | None) -> str:
    if value is None:
        return "-"
    return f"{value.numerator}/{value.denominator}"


def tuple_digest(items: set[object] | tuple[object, ...]) -> str:
    return hashlib.sha256(repr(tuple(sorted(items))).encode()).hexdigest()


def tail_ratio(row: dict[str, object], s2: F) -> F | None:
    r4 = sum(row["qs"][:4], F(0))
    epsilon = row["mass"] - r4 - row["mass"] / 7
    if epsilon <= 0:
        return None
    return s2 * row["components"] / (7 * epsilon)


def strict_cutoff(verifier, row: dict[str, object], s2: F) -> int | None:
    """Largest label forced out by a strict ``c(w)<...`` tail."""

    ratio = tail_ratio(row, s2)
    if ratio is None:
        return None
    return verifier.ceiling(ratio) - 1


def nonstrict_cutoff(row: dict[str, object], s2: F) -> int | None:
    """Largest label retained by a non-strict ``c(w)<=...`` tail."""

    ratio = tail_ratio(row, s2)
    if ratio is None:
        return None
    return ratio.numerator // ratio.denominator


def dependency_controls() -> None:
    for path, expected in (
        (THM1094, THM1094_SHA256),
        (SCOUT_PATH, SCOUT_SHA256),
        (VERIFY_PATH, VERIFY_SHA256),
        (HUNTER_PATH, HUNTER_SHA256),
        (THM2911_OUT, THM2911_OUT_SHA256),
        (THM2912_OUT, THM2912_OUT_SHA256),
        (THM2913_OUT, THM2913_OUT_SHA256),
        (THM2916_SOURCE, THM2916_SOURCE_SHA256),
        (THM2916_OUT, THM2916_OUT_SHA256),
        (THM2920_SOURCE, THM2920_SOURCE_SHA256),
        (THM2920_OUT, THM2920_OUT_SHA256),
    ):
        require(file_sha256(path) == expected, f"{path.name} changed")
    theorem = THM1094.read_text()
    require(
        "|J intersect D_k| <= ell/7 + 6/(49k)." in theorem
        and "Endpoint conventions do not affect any interval length" in theorem,
        "THM1094 sharp interval discrepancy statement changed",
    )
    require(
        "mode=LOCKED_FINITE_EXACT" in THM2911_OUT.read_text()
        and "all_exact_controls=PASS" in THM2911_OUT.read_text()
        and "mode=LOCKED" in THM2912_OUT.read_text()
        and "all_exact_controls=PASS" in THM2912_OUT.read_text()
        and "mode=LOCKED" in THM2913_OUT.read_text()
        and "all_exact_controls=PASS" in THM2913_OUT.read_text(),
        "dependency output lock changed",
    )
    require(
        "mode=LOCKED" in THM2916_OUT.read_text()
        and "all_exact_controls=PASS" in THM2916_OUT.read_text(),
        "THM2916 dependency output lock changed",
    )
    require(
        "mode=LOCKED" in THM2920_OUT.read_text()
        and "all_exact_controls=PASS" in THM2920_OUT.read_text(),
        "THM2920 dependency output lock changed",
    )


def current_proved_union(
    verifier,
    g5,
    rows: list[dict[str, object]],
    pivot_keys: set[tuple[object, ...]],
    pivot_additive: set[tuple[int, ...]],
) -> tuple[
    set[tuple[int, ...]],
    set[tuple[int, ...]],
    set[tuple[int, ...]],
    set[tuple[int, ...]],
]:
    thm2895, thm2898, thm2899, thm2901, thm2902 = (
        g5.canonical_root_sets()
    )
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_body[row["body"]].append(row)

    def hard_roots(keys: set[tuple[object, ...]]) -> set[tuple[int, ...]]:
        return {
            body
            for body, body_rows in by_body.items()
            if all(verifier.branch_key(row) in keys for row in body_rows)
        }

    g5_keys = {
        verifier.branch_key(row) for row in rows if row["margin"] > 0
    }
    prior_fifteen = set(thm2895) | set(thm2898) | set(thm2899) | set(thm2901)
    one_hard = {
        body
        for body, body_rows in by_body.items()
        if len(body_rows) == 1 and body_rows[0]["direct_margin"] <= 0
    }
    through_2903 = prior_fifteen | set(thm2902) | one_hard
    through_2905 = through_2903 | hard_roots(g5_keys)
    through_2904 = through_2905 | pivot_additive
    require(
        (
            len(prior_fifteen),
            len(one_hard),
            len(through_2903),
            len(through_2905),
            len(through_2904),
            len(pivot_keys),
        )
        == (15, 61, 76, 82, 88, 279),
        "proved union through THM2904 changed",
    )
    thm2911_additive = {
        tuple(map(int, line.removeprefix("ADDITIVE_ROOT=").split(",")))
        for line in THM2911_OUT.read_text().splitlines()
        if line.startswith("ADDITIVE_ROOT=")
    }
    thm2912_closed = set(literal_output_value(THM2912_OUT, "closed_roots="))
    thm2913_closed = set(literal_output_value(THM2913_OUT, "closed_roots="))
    through_2913 = (
        through_2904
        | thm2911_additive
        | thm2912_closed
        | thm2913_closed
    )
    require(
        (
            len(thm2911_additive),
            len(thm2912_closed),
            len(thm2913_closed),
            len((through_2904 | thm2911_additive) & thm2912_closed),
            len(
                (through_2904 | thm2911_additive | thm2912_closed)
                & thm2913_closed
            ),
            len(through_2913),
        )
        == (93, 172, 38, 39, 1, 351),
        "proved union through THM2913 changed",
    )
    require(
        tuple_digest(through_2913) == EXPECTED_THROUGH_2913_ROOT_DIGEST,
        "historical 351-root digest changed",
    )
    thm2916_closed = set(
        literal_output_value(THM2916_OUT, "closed_roots=")
    )
    through_2916 = through_2913 | thm2916_closed
    require(
        len(thm2916_closed) == 394
        and not (through_2913 & thm2916_closed)
        and len(through_2916) == 745
        and tuple_digest(thm2916_closed)
        == "c68d09676683f6204df3b04353a3b3107ebbb4285d13a3b6001446372e351e1b",
        "proved union through THM2916 changed",
    )
    thm2920_closed = set(
        literal_output_value(THM2920_OUT, "closed_roots=")
    )
    current = through_2916 | thm2920_closed
    require(
        len(thm2920_closed) == 296
        and not (through_2916 & thm2920_closed)
        and len(current) == 1_041
        and tuple_digest(thm2920_closed)
        == "e3045198e08804c78025bd532111377309882911e08bc50604aa7119ac266c71"
        and tuple_digest(current) == EXPECTED_CURRENT_ROOT_DIGEST,
        "live proved union through THM2920 changed",
    )
    return through_2913, current, thm2916_closed, thm2920_closed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    dependency_controls()

    verifier = load_module("sharp_h1_recompose_verify", VERIFY_PATH)
    scout = load_module("sharp_h1_recompose_scout", SCOUT_PATH)
    g5 = verifier.load_module("sharp_h1_recompose_g5", verifier.G5_SOURCE)
    rows = g5.parse_pair_rows(g5.parse_hard_rows())
    pivot_keys, _pivot_roots, pivot_additive = verifier.verify_pivot_ledger(rows)

    old_h1_keys = {
        verifier.branch_key(row)
        for row in rows
        if (
            (old_cutoff := strict_cutoff(verifier, row, OLD_S2)) is not None
            and old_cutoff <= MAX_CUTOFF
        )
    }
    sharp_all_keys = {
        verifier.branch_key(row)
        for row in rows
        if nonstrict_cutoff(row, SHARP_S2) is not None
    }
    sharp_finite_keys = {
        verifier.branch_key(row)
        for row in rows
        if (
            (sharp_cutoff := nonstrict_cutoff(row, SHARP_S2)) is not None
            and sharp_cutoff <= MAX_CUTOFF
        )
    }
    target_keys = sharp_finite_keys - old_h1_keys
    integral_rows = [
        row
        for row in rows
        if (
            (ratio := tail_ratio(row, SHARP_S2)) is not None
            and ratio.denominator == 1
        )
    ]
    integral_controls = tuple(
        sorted(
            (
                verifier.branch_key(row),
                tail_ratio(row, SHARP_S2).numerator,
                nonstrict_cutoff(row, SHARP_S2),
            )
            for row in integral_rows
        )
    )
    require(
        len(rows) == 14_806
        and len(old_h1_keys) == 5_999
        and len(sharp_all_keys) == 6_180
        and len(sharp_finite_keys) == 6_079
        and len(target_keys) == 80
        and len({key[0] for key in target_keys}) == 80,
        "sharp cutoff census changed",
    )
    require(
        integral_controls
        == (
            (((1, 3, 6, 7, 9, 11, 14), 9, 2, 39, (15, 39)), 3718, 3718),
            (((2, 4, 6, 8, 11, 12, 13), 4, 1, 20, (20,)), 154440, 154440),
            (((3, 5, 6, 8, 9, 11, 12), 10, 2, 30, (21, 30)), 1536, 1536),
        )
        and not any(key in target_keys for key, _ratio, _cutoff in integral_controls),
        "sharp non-strict equality layer changed",
    )

    # The scout reads its discrepancy constant from the imported exact engine.
    # Fork workers inherit this deliberate one-coordinate substitution.
    require(scout.S.S2 == OLD_S2, "scout input discrepancy constant changed")
    scout.S.S2 = SHARP_S2
    tasks = [
        (body, MAX_CUTOFF, MAX_COMBINATIONS, False)
        for body in sorted({key[0] for key in target_keys})
    ]
    tasks.append(
        (EXTENSION_BODY, EXTENSION_CUTOFF, MAX_COMBINATIONS, False)
    )
    if args.workers == 1:
        profiles = list(map(scout.profile_root_task, tasks))
    else:
        require("fork" in mp.get_all_start_methods(), "fork replay unavailable")
        with mp.get_context("fork").Pool(args.workers) as pool:
            profiles = pool.map(scout.profile_root_task, tasks)

    profile_by_body = {profile["body"]: profile for profile in profiles}
    require(len(profile_by_body) == 81, "target profile body count changed")
    exact_rows = []
    for key in sorted(target_keys):
        body, gate_size, rank, apex, prefix = key
        candidates = [
            row
            for row in profile_by_body[body]["rows"]
            if (
                row["K"],
                row["rank"],
                row["apex"],
                row["prefix"],
            )
            == (gate_size, rank, apex, prefix)
        ]
        require(len(candidates) == 1, f"target row reconstruction changed: {key}")
        row = candidates[0]
        parent = next(
            item for item in rows if verifier.branch_key(item) == key
        )
        require(
            tail_ratio(parent, SHARP_S2).denominator != 1
            and row["cutoff"] == nonstrict_cutoff(parent, SHARP_S2),
            f"sharp cutoff mismatch: {key}",
        )
        exact_rows.append((key, row))

    depth1_closed_keys = {
        key
        for key, row in exact_rows
        if row["status"] in ("DEPTH1_CLOSED", "CLOSED")
    }
    depth1_open_rows = [
        (key, row) for key, row in exact_rows if key not in depth1_closed_keys
    ]
    require(
        len(depth1_closed_keys) == 79
        and len(depth1_open_rows) == 1
        and depth1_open_rows[0][0]
        == ((3, 6, 7, 9, 11, 12, 13), 13, 3, 19, (15, 30, 19))
        and depth1_open_rows[0][1]["status"] == "COMBINATION_SKIP",
        "79/1 exact target partition changed",
    )

    # Reuse the independently hash-pinned generalized relative-Hunter engine.
    # Its original wrapper restricted attention to K>=20, a battery choice
    # rather than a mathematical hypothesis, so disable only that wrapper.
    hunter = load_module("sharp_h1_recompose_hunter", HUNTER_PATH)
    require(hunter.S.S.S2 == OLD_S2, "Hunter scout constant changed")
    hunter.S.S.S2 = SHARP_S2
    original_hunter_profile = hunter.S.profile_root_task

    def hunter_profile_without_high_k(task):
        body, max_cutoff, max_combinations, _require_high_k = task
        return original_hunter_profile(
            (body, max_cutoff, max_combinations, False)
        )

    hunter.S.profile_root_task = hunter_profile_without_high_k
    hunter_repair = hunter.audit_target(
        HUNTER_BODY, HUNTER_RANK, HUNTER_PIVOT
    )
    hunter_key = depth1_open_rows[0][0]
    require(
        hunter_repair["H"] == 525
        and hunter_repair["pivot"] == 25
        and hunter_repair["star_hostile_sets"] == 1
        and hunter_repair["hunter_repairs"] == 1
        and hunter_repair["hunter_hard_sets"] == 0
        and hunter_repair["max_hunter"] == F(7_543_517, 37_519_300)
        and hunter_repair["margin"] == F(1_932_131, 28_139_475)
        and hunter_repair["max_row"][1] == (16, 25, 31, 32, 64)
        and hunter_repair["ledger"]
        == "df8fceb7a21cd9773709a8b1f31c911046312290e756f72ae1ff7642bfc42925",
        "sharp survivor Hunter repair changed",
    )
    fixed_closed_keys = depth1_closed_keys | {hunter_key}
    require(
        fixed_closed_keys == target_keys and len(fixed_closed_keys) == 80,
        "fixed-window sharp branch closure changed",
    )

    extension_parent_rows = [
        row
        for row in rows
        if (
            row["body"] == EXTENSION_BODY
            and row["rank"] == 1
            and row["apex"] == 25
            and row["prefix"] == (25,)
        )
    ]
    require(len(extension_parent_rows) == 1, "extension parent key changed")
    extension_parent = extension_parent_rows[0]
    extension_key = verifier.branch_key(extension_parent)
    extension_candidates = [
        row
        for row in profile_by_body[EXTENSION_BODY]["rows"]
        if (
            row["K"],
            row["rank"],
            row["apex"],
            row["prefix"],
        )
        == (
            extension_key[1],
            extension_key[2],
            extension_key[3],
            extension_key[4],
        )
    ]
    require(len(extension_candidates) == 1, "extension replay row changed")
    extension = extension_candidates[0]
    require(
        extension_key
        == ((1, 3, 5, 7, 10, 11, 14), 9, 1, 25, (25,))
        and tail_ratio(extension_parent, SHARP_S2).denominator != 1
        and nonstrict_cutoff(extension_parent, SHARP_S2) == EXTENSION_CUTOFF
        and extension["cutoff"] == EXTENSION_CUTOFF
        and len(extension["core"]) == 831
        and extension["status"] == "DEPTH1_CLOSED"
        and (
            extension["depth1_parent"],
            extension["depth1_local"],
            extension["depth1_short"],
            len(extension["depth1_unresolved"]),
        )
        == (824, 3, 4, 0)
        and extension["depth1_min_margin"] == F(89_371, 111_411_300)
        and extension["depth1_digest"]
        == "6659430a59abd0c529177b424a784763ea312c8f99dc9e0ac1d1fb249964aee7",
        "isolated sharp extension certificate changed",
    )
    certified_sharp_keys = fixed_closed_keys | {extension_key}

    ledger_lines = []
    for key, row in exact_rows:
        parent = next(item for item in rows if verifier.branch_key(item) == key)
        old_cutoff = strict_cutoff(verifier, parent, OLD_S2)
        final_status = (
            "HUNTER_CLOSED" if key == hunter_key else row["status"]
        )
        ledger_lines.append(
            "E="
            + ",".join(map(str, key[0]))
            + f";K={key[1]};rank={key[2]};a={key[3]};"
            + "P="
            + ",".join(map(str, key[4]))
            + f";oldN={old_cutoff};sharpN={row['cutoff']};"
            + f"H={len(row['core'])};status={final_status};"
            + f"d1min={ftext(row['depth1_min_margin'])};"
            + f"d1open={row['depth1_unresolved']};"
            + f"d1digest={row['depth1_digest']}"
        )
    ledger_lines.append(
        "HUNTER="
        + ",".join(map(str, HUNTER_BODY))
        + f";rank={HUNTER_RANK};pivot={HUNTER_PIVOT};"
        + f"star_hostile={hunter_repair['star_hostile_sets']};"
        + f"hunter_repairs={hunter_repair['hunter_repairs']};"
        + f"hunter_hard={hunter_repair['hunter_hard_sets']};"
        + f"maxPsi={ftext(hunter_repair['max_hunter'])};"
        + f"margin={ftext(hunter_repair['margin'])};"
        + f"vertices={hunter_repair['max_row'][1]};"
        + f"digest={hunter_repair['ledger']}"
    )
    ledger_lines.append(
        "EXTENSION="
        + ",".join(map(str, EXTENSION_BODY))
        + f";K={extension_key[1]};rank=1;a=25;P=25;"
        + f"oldN={strict_cutoff(verifier, extension_parent, OLD_S2)};"
        + f"sharpN={extension['cutoff']};H={len(extension['core'])};"
        + f"status={extension['status']};"
        + f"d1min={ftext(extension['depth1_min_margin'])};"
        + f"d1digest={extension['depth1_digest']}"
    )
    ledger = "\n".join(ledger_lines) + "\n"
    ledger_sha256 = hashlib.sha256(ledger.encode()).hexdigest()

    g5_keys = {
        verifier.branch_key(row) for row in rows if row["margin"] > 0
    }
    h4_keys = {
        verifier.branch_key(row) for row in rows if row["pair_margin"] <= 0
    }
    by_body: dict[tuple[int, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_body[row["body"]].append(row)

    def hard_roots(keys: set[tuple[object, ...]]) -> set[tuple[int, ...]]:
        return {
            body
            for body, body_rows in by_body.items()
            if all(verifier.branch_key(row) in keys for row in body_rows)
        }

    scalar_roots = set(g5.canonical_root_sets()[2])
    old_route_keys = old_h1_keys | g5_keys | pivot_keys | h4_keys
    fixed_route_keys = old_route_keys | fixed_closed_keys
    sharp_route_keys = old_route_keys | certified_sharp_keys
    saturation_route_keys = sharp_all_keys | g5_keys | pivot_keys | h4_keys
    old_route_roots = hard_roots(old_route_keys) | scalar_roots
    fixed_route_roots = hard_roots(fixed_route_keys) | scalar_roots
    sharp_route_roots = hard_roots(sharp_route_keys) | scalar_roots
    saturation_route_roots = (
        hard_roots(saturation_route_keys) | scalar_roots
    )
    (
        through_2913_roots,
        current_roots,
        thm2916_roots,
        _thm2920_roots,
    ) = current_proved_union(
        verifier, g5, rows, pivot_keys, set(pivot_additive)
    )
    additive_route_roots = sharp_route_roots - old_route_roots
    additive_through_2913_roots = sharp_route_roots - through_2913_roots
    additive_current_roots = sharp_route_roots - current_roots
    proved_union = current_roots | sharp_route_roots
    saturation_gain_over_fixed = saturation_route_roots - fixed_route_roots
    saturation_gain_over_sharp = saturation_route_roots - sharp_route_roots
    high_cutoff_keys = sharp_all_keys - sharp_finite_keys
    saturation_missing = {
        body: tuple(
            sorted(
                verifier.branch_key(row)
                for row in by_body[body]
                if verifier.branch_key(row) not in fixed_route_keys
            )
        )
        for body in saturation_gain_over_fixed
    }
    saturation_missing_cutoffs = {
        body: tuple(
            nonstrict_cutoff(
                next(
                    row
                    for row in rows
                    if verifier.branch_key(row) == key
                ),
                SHARP_S2,
            )
            for key in keys
        )
        for body, keys in saturation_missing.items()
    }
    saturation_new_beyond_351 = (
        saturation_gain_over_fixed - through_2913_roots
    )
    saturation_new_beyond_current = (
        saturation_gain_over_fixed - current_roots
    )
    historically_additive_high_keys = {
        key
        for body in saturation_new_beyond_351
        for key in saturation_missing[body]
    }
    nonadditive_high_keys = high_cutoff_keys - historically_additive_high_keys
    require(
        len(high_cutoff_keys) == 101
        and saturation_gain_over_fixed
        == {
            (1, 2, 3, 4, 5, 6, 14),
            (1, 2, 3, 4, 6, 7, 12),
            (1, 2, 3, 4, 6, 9, 12),
            (1, 2, 4, 8, 9, 10, 11),
            (1, 2, 4, 8, 11, 13, 14),
            (1, 3, 5, 6, 7, 9, 14),
            EXTENSION_BODY,
        }
        and saturation_missing[EXTENSION_BODY] == (extension_key,)
        and saturation_missing_cutoffs[EXTENSION_BODY] == (18_869,)
        and len(saturation_missing[(1, 2, 3, 4, 5, 6, 14)]) == 1
        and saturation_missing_cutoffs[(1, 2, 3, 4, 5, 6, 14)]
        == (34_222,)
        and saturation_new_beyond_351
        == {(1, 2, 3, 4, 5, 6, 14), EXTENSION_BODY}
        and saturation_new_beyond_current == {EXTENSION_BODY}
        and saturation_gain_over_sharp
        == {
            (1, 2, 3, 4, 5, 6, 14),
            (1, 2, 3, 4, 6, 7, 12),
            (1, 2, 3, 4, 6, 9, 12),
            (1, 2, 4, 8, 9, 10, 11),
            (1, 2, 4, 8, 11, 13, 14),
            (1, 3, 5, 6, 7, 9, 14),
        }
        and len(historically_additive_high_keys) == 2
        and len(nonadditive_high_keys) == 99,
        "sharp route-saturation ceiling changed",
    )

    atoms = Counter()
    for key in sharp_route_keys:
        atom = "".join(
            label
            for label, keys in (
                ("E", h4_keys),
                ("G", g5_keys),
                ("H", old_h1_keys),
                ("P", pivot_keys),
                ("S", certified_sharp_keys),
            )
            if key in keys
        )
        atoms[atom] += 1

    closed_core_sizes = [
        len(row["core"]) for key, row in exact_rows if key in depth1_closed_keys
    ]
    positive_depth1_margins = [
        row["depth1_min_margin"]
        for key, row in exact_rows
        if key in depth1_closed_keys and row["depth1_min_margin"] is not None
    ]
    old_cutoffs = [
        strict_cutoff(
            verifier,
            next(item for item in rows if verifier.branch_key(item) == key),
            OLD_S2,
        )
        for key, _row in exact_rows
    ]
    sharp_cutoffs = [row["cutoff"] for _key, row in exact_rows]
    survivor = depth1_open_rows[0][1]
    counts = (
        len(rows),
        len(old_h1_keys),
        len(sharp_all_keys),
        len(sharp_finite_keys),
        len(target_keys),
        len(depth1_closed_keys),
        len(depth1_open_rows),
        len(fixed_closed_keys),
        len(certified_sharp_keys),
        len(integral_controls),
        min(old_cutoffs),
        max(old_cutoffs),
        min(sharp_cutoffs),
        max(sharp_cutoffs),
        min(closed_core_sizes),
        max(closed_core_sizes),
        len(survivor["core"]),
        len(survivor["depth1_unresolved"]),
        min(positive_depth1_margins),
        len(g5_keys),
        len(pivot_keys),
        len(h4_keys),
        len(old_route_keys),
        len(fixed_route_keys),
        len(sharp_route_keys),
        len(saturation_route_keys),
        len(old_route_roots),
        len(fixed_route_roots),
        len(sharp_route_roots),
        len(saturation_route_roots),
        len(additive_route_roots),
        len(through_2913_roots),
        len(additive_through_2913_roots),
        len(current_roots),
        len(additive_current_roots),
        len(proved_union),
        3_432 - len(proved_union),
        len(high_cutoff_keys),
        len(saturation_gain_over_fixed),
        len(saturation_gain_over_sharp),
        len(saturation_new_beyond_351),
        len(saturation_new_beyond_current),
        len(nonadditive_high_keys),
    )

    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "locked count tuple changed")
    target_digest = tuple_digest(target_keys)
    closed_digest = tuple_digest(fixed_closed_keys)
    certified_digest = tuple_digest(certified_sharp_keys)
    if EXPECTED_TARGET_KEY_DIGEST is not None:
        require(
            target_digest == EXPECTED_TARGET_KEY_DIGEST,
            "target-key digest changed",
        )
    if EXPECTED_CLOSED_KEY_DIGEST is not None:
        require(
            closed_digest == EXPECTED_CLOSED_KEY_DIGEST,
            "closed-key digest changed",
        )
    if EXPECTED_CERTIFIED_KEY_DIGEST is not None:
        require(
            certified_digest == EXPECTED_CERTIFIED_KEY_DIGEST,
            "certified-key digest changed",
        )
    if EXPECTED_LEDGER_SHA256 is not None:
        require(
            ledger_sha256 == EXPECTED_LEDGER_SHA256,
            "exact target ledger changed",
        )
    require(
        tuple_digest(fixed_route_roots) == EXPECTED_FIXED_ROUTE_ROOT_DIGEST,
        "fixed-window route-root digest changed",
    )
    if EXPECTED_SHARP_ROUTE_ROOT_DIGEST is not None:
        require(
            tuple_digest(sharp_route_roots)
            == EXPECTED_SHARP_ROUTE_ROOT_DIGEST,
            "sharp route-root digest changed",
        )
    if EXPECTED_SATURATION_ROOT_DIGEST is not None:
        require(
            tuple_digest(saturation_route_roots)
            == EXPECTED_SATURATION_ROOT_DIGEST,
            "saturation route-root digest changed",
        )
    require(
        additive_route_roots
        == additive_through_2913_roots
        == {
            (2, 4, 8, 10, 11, 13, 14),
            EXTENSION_BODY,
        }
        and additive_current_roots == {EXTENSION_BODY}
        and (sharp_route_roots & thm2916_roots)
        - (old_route_roots & thm2916_roots)
        == {(2, 4, 8, 10, 11, 13, 14)}
        and len(proved_union) == 1_042,
        "sharp whole-root addition changed",
    )

    new_root_profiles = []
    for new_root in sorted(additive_through_2913_roots):
        profile = []
        for row in sorted(by_body[new_root], key=lambda item: item["rank"]):
            key = verifier.branch_key(row)
            routes = "".join(
                label
                for label, keys in (
                    ("E", h4_keys),
                    ("G", g5_keys),
                    ("H", old_h1_keys),
                    ("P", pivot_keys),
                    ("S", certified_sharp_keys),
                )
                if key in keys
            )
            profile.append(
                (
                    row["rank"],
                    row["apex"],
                    row["prefix"],
                    routes,
                    strict_cutoff(verifier, row, OLD_S2),
                    nonstrict_cutoff(row, SHARP_S2),
                )
            )
        new_root_profiles.append((new_root, tuple(profile)))

    print("LRC14 sharp ranked-H1 tail recomposition")
    print("sharp_tail=component_union_of_THM1094_section2;S2=6/7")
    print(f"parameters=workers:{args.workers},max_cutoff:{MAX_CUTOFF}")
    print(f"counts={counts}")
    print(f"integral_equality_controls={integral_controls}")
    print(f"branch_atoms={tuple(sorted(atoms.items()))}")
    print(f"target_key_digest={target_digest}")
    print(f"closed_key_digest={closed_digest}")
    print(f"certified_key_digest={certified_digest}")
    print(f"ledger_sha256={ledger_sha256}")
    print(f"through_2913_root_digest={tuple_digest(through_2913_roots)}")
    print(f"current_root_digest={tuple_digest(current_roots)}")
    print(f"fixed_route_root_digest={tuple_digest(fixed_route_roots)}")
    print(f"sharp_route_root_digest={tuple_digest(sharp_route_roots)}")
    print(f"saturation_root_digest={tuple_digest(saturation_route_roots)}")
    print(f"depth1_open_key={depth1_open_rows[0][0]}")
    print(f"depth1_open_unresolved={survivor['depth1_unresolved']}")
    print(
        "hunter_repair="
        f"{hunter_key};maxPsi={ftext(hunter_repair['max_hunter'])};"
        f"margin={ftext(hunter_repair['margin'])};"
        f"vertices={hunter_repair['max_row'][1]};"
        f"ledger={hunter_repair['ledger']}"
    )
    print(f"extension_key={extension_key}")
    print(f"additive_route_roots={tuple(sorted(additive_route_roots))}")
    print(
        "additive_through_2913_roots="
        f"{tuple(sorted(additive_through_2913_roots))}"
    )
    print(f"additive_current_roots={tuple(sorted(additive_current_roots))}")
    print(f"new_root_profiles={tuple(new_root_profiles)}")
    print(f"high_cutoff_keys={len(high_cutoff_keys)}")
    print(
        "saturation_gain_over_fixed="
        f"{tuple(sorted(saturation_gain_over_fixed))}"
    )
    print(
        "saturation_new_beyond_351="
        f"{tuple(sorted(saturation_new_beyond_351))}"
    )
    print(
        "saturation_new_beyond_current="
        f"{tuple(sorted(saturation_new_beyond_current))}"
    )
    print(f"nonadditive_high_keys={len(nonadditive_high_keys)}")
    print(
        "saturation_missing="
        f"{tuple(sorted(saturation_missing.items()))}"
    )
    print(
        "saturation_missing_cutoffs="
        f"{tuple(sorted(saturation_missing_cutoffs.items()))}"
    )
    print(f"ledger_rows={len(ledger_lines)}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
