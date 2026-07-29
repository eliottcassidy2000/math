#!/usr/bin/env python3
r"""Exact pair-cap atlas for the eight-body / five-slot LRC(14) rung.

THM-2885 assigns ten globally ranked first apices to every
``E in C({1,...,14},8)``.  For each of the resulting 30,030 literal carriers

    C(E;a) = G_E \ D_a,          h=|C(E;a)|,
    c(w) = |C(E;a) intersect D_w|,

this verifier constructs a global cap for the exact two-comb union

    U(x,y) = |C(E;a) intersect (D_x union D_y)|.

The strict discrepancy estimate

    c(w) < h/7 + (99/70) r/(7w)

and an exact scan through 2500 seal both the leading scalar coverages and the
pair tail.  Branch-and-bound computes the maximum exact pair union in the
finite head.  If the resulting global cap ``Ubar`` is below ``5h/7``, then
with

    delta = 5h/7-Ubar,
    W = floor(2(99/70)r/(7 delta)),

every possible four-speed obstruction has at most one speed above ``W``.
If ``2Ubar<h``, an arbitrary 2+2 partition closes all four slots directly.

Exactly five carriers fail ``Ubar<5h/7``.  Each has exactly one head pair
whose exact union reaches ``5h/7``.  The literal residual after this forced
pair has a global two-comb union cap below its mass, so every branch
containing that pair is closed.  Branches avoiding it, and the one-tail/all-
head obligations of the positive pair caps, remain separate obligations.
This script does not prove the five-slot rung or LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM2885_PATH = (
    ROOT
    / "04-computation/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.py"
)
THM2885_SHA256 = (
    "dff97f67b1104c25589802a6a2f216b6e7bfedd58eebfa1bcce615d59c1e872f"
)
THM2883_PATH = (
    ROOT
    / "04-computation/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.py"
)
THM2883_SHA256 = (
    "a5f3dcc1a23defea4b3dc067675d83141f1866022d6d01946617a97de69e5b0e"
)

FIRST_EXTERNAL = 15
PAIR_HORIZON = 2500
S2 = F(99, 70)

EXPECTED_SCALAR_COUNTS = (13_969, 10_939, 5_122)
EXPECTED_PAIR_COUNTS = (30_025, 5, 13_802, 53)
EXPECTED_PAIR_WORK = (191_306, 353)
EXPECTED_MINIMUM_POSITIVE = (
    F(709, 8_828_820),
    (2, 4, 6, 7, 8, 10, 12, 14),
    3,
    26,
    130_827,
    (18, 22),
)
EXPECTED_HOSTILES = (
    (
        (1, 2, 3, 7, 8, 10, 11, 13),
        5,
        19,
        17,
        18,
        F(-50_887, 285_170_886),
    ),
    (
        (1, 2, 3, 7, 8, 10, 11, 13),
        7,
        17,
        18,
        19,
        F(-3_650_923, 2_851_708_860),
    ),
    (
        (1, 4, 6, 7, 9, 10, 11, 13),
        6,
        23,
        16,
        17,
        F(-480_313, 920_551_632),
    ),
    (
        (2, 3, 4, 8, 10, 12, 13, 14),
        6,
        16,
        18,
        22,
        F(-23, 1_070_160),
    ),
    (
        (2, 4, 6, 8, 10, 12, 13, 14),
        4,
        16,
        18,
        22,
        F(-3_553, 1_070_160),
    ),
)
EXPECTED_FORCED_RESIDUALS = (
    (
        (1, 2, 3, 7, 8, 10, 11, 13),
        19,
        17,
        18,
        F(2_597_104, 212_952_285),
        21,
        23,
        38,
    ),
    (
        (1, 2, 3, 7, 8, 10, 11, 13),
        17,
        18,
        19,
        F(2_597_104, 212_952_285),
        21,
        23,
        38,
    ),
    (
        (1, 4, 6, 7, 9, 10, 11, 13),
        23,
        16,
        17,
        F(389_293, 63_488_880),
        19,
        108,
        41,
    ),
    (
        (2, 3, 4, 8, 10, 12, 13, 14),
        16,
        18,
        22,
        F(15_287, 840_840),
        26,
        40,
        1,
    ),
    (
        (2, 4, 6, 8, 10, 12, 13, 14),
        16,
        18,
        22,
        F(356, 35_035),
        20,
        26,
        4,
    ),
)

# Filled and made mandatory after the first exact discovery replay.
EXPECTED_PROFILE_DIGEST: str | None = None
EXPECTED_PAID_PAIR_DIGEST: str | None = None
EXPECTED_THRESHOLD_DIGEST: str | None = None
EXPECTED_FORCED_DIGEST: str | None = None
EXPECTED_NORMALIZATION_DIGEST: str | None = None
EXPECTED_ROOT_GATE_DIGEST: str | None = None


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(path: Path, expected_sha256: str, name: str):
    require(file_sha256(path) == expected_sha256, f"{path.name} hash changed")
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


THM2885 = load_module(THM2885_PATH, THM2885_SHA256, "thm2888_thm2885")
THM2883 = load_module(THM2883_PATH, THM2883_SHA256, "thm2888_thm2883")
CORE = THM2885.CORE
BODIES = THM2885.BODIES


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def carrier_mass(good: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in good), F(0))


def pair_maximum(
    good: list[tuple[F, F]],
    mass: F,
    ranked: list[tuple[F, int]],
    direct_control_family: tuple[int, ...] | None,
    excluded_pair: tuple[int, int] | None = None,
) -> dict[str, object]:
    """Branch-and-bound the exact pair maximum in the finite head."""

    head_cap = F(0)
    maximizing_pair: tuple[int, int] | None = None
    paid = 0
    rows: list[str] = []
    first_control_done = False
    for first_index, (first_value, first) in enumerate(ranked[:-1]):
        if first_value + ranked[first_index + 1][0] <= head_cap:
            break
        after_first = THM2883.subtract_local(good, first)
        for second_value, second in ranked[first_index + 1 :]:
            if first_value + second_value <= head_cap:
                break
            pair = tuple(sorted((first, second)))
            if pair == excluded_pair:
                continue
            after_pair = THM2883.subtract_local(after_first, second)
            survivor = carrier_mass(after_pair)
            union = mass - survivor
            rows.append(
                f"P={pair[0]},{pair[1]};U={ftext(union)};"
                f"h={ftext(survivor)};r={len(after_pair)}\n"
            )
            paid += 1
            if union > head_cap:
                head_cap = union
                maximizing_pair = pair
            if not first_control_done:
                require(
                    after_pair
                    == THM2883.subtract_local(
                        THM2883.subtract_local(good, second),
                        first,
                    )
                    == THM2883.subtract_local_multi(good, pair),
                    f"pair subtraction paths disagree: {pair}",
                )
                if direct_control_family is not None:
                    family = tuple(sorted((*direct_control_family, *pair)))
                    direct_good, direct_r, direct_m = CORE.good_norm(family)
                    require(
                        direct_good == after_pair
                        and direct_r == len(after_pair)
                        and direct_m == survivor,
                        f"direct pair reconstruction disagrees: {family}",
                    )
                first_control_done = True
    require(
        paid > 0 and maximizing_pair is not None and first_control_done,
        "empty pair branch-and-bound",
    )
    return {
        "head_cap": head_cap,
        "maximizing_pair": maximizing_pair,
        "paid": paid,
        "paid_digest": hashlib.sha256(
            (
                "THM2888/paid-pairs/v1;"
                + (
                    "excluded=none\n"
                    if excluded_pair is None
                    else f"excluded={excluded_pair[0]},{excluded_pair[1]}\n"
                )
                + "".join(rows)
            ).encode()
        ).hexdigest(),
    }


def threshold_pairs(
    good: list[tuple[F, F]],
    mass: F,
    ranked: list[tuple[F, int]],
    threshold: F,
) -> tuple[tuple[tuple[int, int, F], ...], int, str]:
    """Enumerate every pair whose scalar sum can reach ``threshold``."""

    candidates: list[str] = []
    heavy: list[tuple[int, int, F]] = []
    checks = 0
    for first_index, (first_value, first) in enumerate(ranked[:-1]):
        if first_value + ranked[first_index + 1][0] < threshold:
            break
        after_first = THM2883.subtract_local(good, first)
        for second_value, second in ranked[first_index + 1 :]:
            if first_value + second_value < threshold:
                break
            after_pair = THM2883.subtract_local(after_first, second)
            union = mass - carrier_mass(after_pair)
            pair = tuple(sorted((first, second)))
            candidates.append(
                f"P={pair[0]},{pair[1]};U={ftext(union)}\n"
            )
            checks += 1
            if union >= threshold:
                heavy.append((pair[0], pair[1], union))
    digest = hashlib.sha256(
        ("THM2888/threshold-pairs/v1\n" + "".join(candidates)).encode()
    ).hexdigest()
    return tuple(heavy), checks, digest


def profile_apex(
    body: tuple[int, ...],
    root_good: list[tuple[F, F]],
    apex_rank: int,
    apex: int,
) -> dict[str, object]:
    carrier = THM2883.subtract_local(root_good, apex)
    components = len(carrier)
    mass = carrier_mass(carrier)
    direct_r, direct_m, direct_good = CORE.subtract(root_good, apex)
    require(
        direct_good == carrier
        and direct_r == components
        and direct_m == mass
        and mass > 0,
        f"first-apex carrier mismatch: {body}, {apex}",
    )

    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, PAIR_HORIZON + 1)
        if speed != apex
    ]
    rows = THM2885.coverages_many(carrier, speeds)
    by_speed = {speed: value for value, speed in rows}
    for control_speed in (speeds[0], speeds[-1]):
        require(
            by_speed[control_speed]
            == THM2885.coverage(carrier, control_speed),
            f"vector/scalar coverage mismatch: {body}, {apex}, {control_speed}",
        )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    q1 = ranked[0][0]
    q4 = ranked[3][0]
    require(
        q4 > mass / 7
        and mass / 7 + S2 * components / (7 * (PAIR_HORIZON + 1)) <= q4,
        f"global top-four tail seal failed: {body}, {apex}",
    )
    rank1_threshold = S2 * components / (7 * (q1 - mass / 7))
    require(
        rank1_threshold < PAIR_HORIZON + 1,
        f"global rank-one tail seal failed: {body}, {apex}",
    )
    q1_margin = F(4, 7) * mass - q1
    require(q1_margin > 0, f"q1 reaches 4h/7: {body}, {apex}")
    pair_entry = S2 * components / (7 * q1_margin)
    tail_single_cap = mass / 7 + S2 * components / (
        7 * (PAIR_HORIZON + 1)
    )
    tail_pair_cap = q1 + tail_single_cap
    require(
        pair_entry < PAIR_HORIZON + 1
        and tail_pair_cap < F(5, 7) * mass,
        f"pair tail misses 5h/7: {body}, {apex}",
    )

    direct_control_family = (*body, apex) if apex_rank == 1 else None
    pair = pair_maximum(
        carrier,
        mass,
        ranked,
        direct_control_family,
    )
    global_cap = max(pair["head_cap"], tail_pair_cap)
    margin = F(5, 7) * mass - global_cap

    heavy, threshold_checks, threshold_digest = threshold_pairs(
        carrier,
        mass,
        ranked,
        F(5, 7) * mass,
    )
    if margin > 0:
        require(not heavy, f"positive cap has a heavy pair: {body}, {apex}")
        cutoff = max(
            14,
            floor_fraction(
                2 * S2 * components / (7 * margin)
            ),
        )
        require(
            global_cap
            + 2
            * (
                mass / 7
                + S2 * components / (7 * (cutoff + 1))
            )
            < mass,
            f"two-tail cutoff failed: {body}, {apex}",
        )
    else:
        cutoff = None
        require(heavy, f"nonpositive cap has no heavy pair: {body}, {apex}")

    direct_margin = mass - sum(
        (value for value, _ in ranked[:4]),
        F(0),
    )
    residual3 = mass - sum(
        (value for value, _ in ranked[:3]),
        F(0),
    )
    scalar_class = (
        "direct"
        if direct_margin > 0
        else ("rank3" if residual3 > mass / 7 else "failure")
    )
    return {
        "body": body,
        "apex_rank": apex_rank,
        "apex": apex,
        "carrier": carrier
        if (
            margin <= 0
            or (
                scalar_class == "failure"
                and cutoff is not None
                and cutoff > PAIR_HORIZON
            )
        )
        else None,
        "m": mass,
        "r": components,
        "top4": tuple(ranked[:4]),
        "direct_margin": direct_margin,
        "residual3": residual3,
        "scalar_class": scalar_class,
        "rank1_threshold": rank1_threshold,
        "pair_entry": pair_entry,
        "head_cap": pair["head_cap"],
        "tail_pair_cap": tail_pair_cap,
        "global_cap": global_cap,
        "margin": margin,
        "cutoff": cutoff,
        "pairpair_direct": 2 * global_cap < mass,
        "maximizing_pair": pair["maximizing_pair"],
        "paid": pair["paid"],
        "paid_digest": pair["paid_digest"],
        "heavy": heavy,
        "threshold_checks": threshold_checks,
        "threshold_digest": threshold_digest,
    }


def profile_body(body: tuple[int, ...]) -> list[dict[str, object]]:
    root = THM2885.profile_body(body)
    root_good, root_r, root_m = CORE.good_norm(body)
    require(
        root["m"] == root_m and root["r"] == root_r,
        f"root reconstruction changed: {body}",
    )
    return [
        profile_apex(body, root_good, rank, apex)
        for rank, (_, apex) in enumerate(root["top15"][:10], start=1)
    ]


def forced_residual(
    row: dict[str, object],
    forced_pair: tuple[int, int] | None = None,
) -> dict[str, object]:
    carrier = row["carrier"]
    heavy = row["heavy"]
    require(
        isinstance(carrier, list),
        f"bad hostile carrier: {row['body']}, {row['apex']}",
    )
    if forced_pair is None:
        require(len(heavy) == 1, "hostile carrier does not have one heavy pair")
        first, second, _ = heavy[0]
    else:
        first, second = forced_pair
    residual = THM2883.subtract_local_multi(carrier, (first, second))
    reverse = THM2883.subtract_local(
        THM2883.subtract_local(carrier, second),
        first,
    )
    require(residual == reverse, "forced residual order changed")
    mass = carrier_mass(residual)
    components = len(residual)
    family = tuple(sorted((*row["body"], row["apex"], first, second)))
    direct_good, direct_r, direct_m = CORE.good_norm(family)
    require(
        direct_good == residual
        and direct_r == components
        and direct_m == mass,
        f"forced residual direct reconstruction failed: {family}",
    )

    excluded = {row["apex"], first, second}
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, PAIR_HORIZON + 1)
        if speed not in excluded
    ]
    coverage_rows = THM2885.coverages_many(residual, speeds)
    ranked = sorted(coverage_rows, key=lambda item: (-item[0], item[1]))
    q1 = ranked[0][0]
    require(q1 > mass / 7, f"forced residual q1 misses limit: {family}")
    tail_single_cap = mass / 7 + S2 * components / (
        7 * (PAIR_HORIZON + 1)
    )
    require(
        tail_single_cap <= q1,
        f"forced residual rank-one tail failed: {family}",
    )
    pair = pair_maximum(residual, mass, ranked, None)
    tail_pair_cap = q1 + tail_single_cap
    global_cap = max(pair["head_cap"], tail_pair_cap)
    margin = mass - global_cap
    require(margin > 0, f"forced pair residual remains open: {family}")
    return {
        "body": row["body"],
        "apex": row["apex"],
        "pair": (first, second),
        "m": mass,
        "r": components,
        "head_cap": pair["head_cap"],
        "tail_pair_cap": tail_pair_cap,
        "global_cap": global_cap,
        "margin": margin,
        "maximizing_pair": pair["maximizing_pair"],
        "paid": pair["paid"],
        "paid_digest": pair["paid_digest"],
    }


def normalize_heavy_edge(row: dict[str, object]) -> dict[str, object]:
    """Close the maximizing edge and recompute the cap with it deleted."""

    carrier = row["carrier"]
    edge = row["maximizing_pair"]
    require(
        isinstance(carrier, list)
        and isinstance(edge, tuple)
        and len(edge) == 2,
        "heavy-edge normalization input changed",
    )
    forced = forced_residual(row, edge)
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, PAIR_HORIZON + 1)
        if speed != row["apex"]
    ]
    coverage_rows = THM2885.coverages_many(carrier, speeds)
    ranked = sorted(coverage_rows, key=lambda item: (-item[0], item[1]))
    deleted = pair_maximum(
        carrier,
        row["m"],
        ranked,
        None,
        excluded_pair=edge,
    )
    deleted_cap = max(deleted["head_cap"], row["tail_pair_cap"])
    deleted_margin = F(5, 7) * row["m"] - deleted_cap
    require(deleted_margin > 0, "heavy-edge deletion leaves nonpositive cap")
    deleted_cutoff = max(
        14,
        floor_fraction(
            2 * S2 * row["r"] / (7 * deleted_margin)
        ),
    )
    require(
        deleted_cap
        + 2
        * (
            row["m"] / 7
            + S2 * row["r"] / (7 * (deleted_cutoff + 1))
        )
        < row["m"],
        "deleted-edge two-tail cutoff failed",
    )
    return {
        "body": row["body"],
        "apex_rank": row["apex_rank"],
        "apex": row["apex"],
        "edge": edge,
        "old_margin": row["margin"],
        "old_cutoff": row["cutoff"],
        "forced": forced,
        "deleted_head_cap": deleted["head_cap"],
        "deleted_cap": deleted_cap,
        "deleted_margin": deleted_margin,
        "deleted_cutoff": deleted_cutoff,
        "deleted_maximizing_pair": deleted["maximizing_pair"],
        "deleted_paid": deleted["paid"],
        "deleted_digest": deleted["paid_digest"],
    }


def digest_text(header: str, rows: list[str]) -> str:
    return hashlib.sha256((header + "".join(rows)).encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(6, os.cpu_count() or 1),
    )
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(
        len(BODIES) == 3003
        and len(set(BODIES)) == 3003
        and S2**2 > 2,
        "pair atlas universe or sqrt(2) majorant changed",
    )

    if args.workers == 1:
        nested = list(map(profile_body, BODIES))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            nested = pool.map(profile_body, BODIES, chunksize=1)
    rows = [row for body_rows in nested for row in body_rows]
    require(
        len(rows) == 30_030
        and tuple(row["body"] for row in rows[::10]) == BODIES,
        "first-apex atlas order changed",
    )

    scalar_counts = tuple(
        sum(row["scalar_class"] == name for row in rows)
        for name in ("direct", "rank3", "failure")
    )
    positive = [row for row in rows if row["margin"] > 0]
    hostile = [row for row in rows if row["margin"] <= 0]
    pair_counts = (
        len(positive),
        len(hostile),
        sum(row["pairpair_direct"] for row in rows),
        sum(
            row["pairpair_direct"] and row["scalar_class"] == "failure"
            for row in rows
        ),
    )
    require(
        scalar_counts == EXPECTED_SCALAR_COUNTS,
        "scalar first-apex census changed",
    )
    require(pair_counts == EXPECTED_PAIR_COUNTS, "pair-cap census changed")
    require(
        (
            sum(row["paid"] for row in rows),
            sum(row["threshold_checks"] for row in rows),
        )
        == EXPECTED_PAIR_WORK,
        "pair work census changed",
    )

    hostile_ledger = tuple(
        (
            row["body"],
            row["apex_rank"],
            row["apex"],
            row["heavy"][0][0],
            row["heavy"][0][1],
            row["margin"],
        )
        for row in hostile
    )
    require(hostile_ledger == EXPECTED_HOSTILES, "hostile atlas changed")

    minimum_positive_row = min(
        positive,
        key=lambda row: (
            row["margin"],
            row["body"],
            row["apex_rank"],
            row["apex"],
        ),
    )
    minimum_positive = (
        minimum_positive_row["margin"],
        minimum_positive_row["body"],
        minimum_positive_row["apex_rank"],
        minimum_positive_row["apex"],
        minimum_positive_row["cutoff"],
        minimum_positive_row["maximizing_pair"],
    )
    require(
        minimum_positive == EXPECTED_MINIMUM_POSITIVE,
        "minimum positive pair cap changed",
    )

    high_cutoff = [
        row
        for row in rows
        if row["scalar_class"] == "failure"
        and row["margin"] > 0
        and row["cutoff"] > PAIR_HORIZON
    ]
    require(len(high_cutoff) == 29, "high-cutoff scalar-failure census changed")
    normalization_sources = [
        row
        for row in rows
        if row["scalar_class"] == "failure"
        and (
            row["margin"] <= 0
            or (
                row["cutoff"] is not None
                and row["cutoff"] > PAIR_HORIZON
            )
        )
    ]
    require(
        len(normalization_sources) == 34,
        "heavy-edge normalization universe changed",
    )
    normalized = [normalize_heavy_edge(row) for row in normalization_sources]
    normalized_by_key = {
        (row["body"], row["apex"]): row
        for row in normalized
    }
    forced = [
        normalized_by_key[(row["body"], row["apex"])]["forced"]
        for row in hostile
    ]
    forced_ledger = tuple(
        (
            row["body"],
            row["apex"],
            row["pair"][0],
            row["pair"][1],
            row["margin"],
            row["maximizing_pair"][0],
            row["maximizing_pair"][1],
            row["paid"],
        )
        for row in forced
    )
    require(
        forced_ledger == EXPECTED_FORCED_RESIDUALS,
        "forced residual atlas changed",
    )

    profile_lines = [
        "E="
        + ",".join(map(str, row["body"]))
        + f";rank={row['apex_rank']};a={row['apex']};"
        + f"h={ftext(row['m'])};r={row['r']};"
        + f"class={row['scalar_class']};"
        + "top4="
        + ",".join(
            f"{speed}:{ftext(value)}"
            for value, speed in row["top4"]
        )
        + f";T1={ftext(row['rank1_threshold'])};"
        + f"T2={ftext(row['pair_entry'])};"
        + f"H={ftext(row['head_cap'])};"
        + f"tail={ftext(row['tail_pair_cap'])};"
        + f"Ubar={ftext(row['global_cap'])};"
        + f"delta={ftext(row['margin'])};"
        + f"W={row['cutoff'] if row['cutoff'] is not None else 'NA'};"
        + f"pairpair={int(row['pairpair_direct'])};"
        + "maxpair="
        + ",".join(map(str, row["maximizing_pair"]))
        + f";paid={row['paid']};heavy={len(row['heavy'])};"
        + f"paidsha={row['paid_digest']};"
        + f"thresholdsha={row['threshold_digest']}\n"
        for row in rows
    ]
    profile_digest = digest_text(
        "THM2888/first-apex-pair-profile/v1\n",
        profile_lines,
    )
    paid_digest = digest_text(
        "THM2888/hierarchical-paid-pair-ledger/v1\n",
        [row["paid_digest"] + "\n" for row in rows],
    )
    threshold_digest = digest_text(
        "THM2888/hierarchical-threshold-pair-ledger/v1\n",
        [row["threshold_digest"] + "\n" for row in rows],
    )
    forced_lines = [
        "E="
        + ",".join(map(str, row["body"]))
        + f";a={row['apex']};P={row['pair'][0]},{row['pair'][1]};"
        + f"h={ftext(row['m'])};r={row['r']};"
        + f"H={ftext(row['head_cap'])};"
        + f"tail={ftext(row['tail_pair_cap'])};"
        + f"cap={ftext(row['global_cap'])};"
        + f"margin={ftext(row['margin'])};"
        + "maxpair="
        + ",".join(map(str, row["maximizing_pair"]))
        + f";paid={row['paid']};sha={row['paid_digest']}\n"
        for row in forced
    ]
    forced_digest = digest_text(
        "THM2888/forced-heavy-pair-residual/v1\n",
        forced_lines,
    )
    normalization_lines = [
        "E="
        + ",".join(map(str, row["body"]))
        + f";rank={row['apex_rank']};a={row['apex']};"
        + f"edge={row['edge'][0]},{row['edge'][1]};"
        + f"old_delta={ftext(row['old_margin'])};"
        + f"old_W={row['old_cutoff'] if row['old_cutoff'] is not None else 'NA'};"
        + f"forced_margin={ftext(row['forced']['margin'])};"
        + f"forced_paid={row['forced']['paid']};"
        + f"deleted_H={ftext(row['deleted_head_cap'])};"
        + f"deleted_cap={ftext(row['deleted_cap'])};"
        + f"deleted_delta={ftext(row['deleted_margin'])};"
        + f"deleted_W={row['deleted_cutoff']};"
        + "deleted_maxpair="
        + ",".join(map(str, row["deleted_maximizing_pair"]))
        + f";deleted_paid={row['deleted_paid']};"
        + f"deleted_sha={row['deleted_digest']}\n"
        for row in normalized
    ]
    normalization_digest = digest_text(
        "THM2888/heavy-edge-normalization/v1\n",
        normalization_lines,
    )
    if EXPECTED_PROFILE_DIGEST is not None:
        require(profile_digest == EXPECTED_PROFILE_DIGEST, "profile digest changed")
    if EXPECTED_PAID_PAIR_DIGEST is not None:
        require(paid_digest == EXPECTED_PAID_PAIR_DIGEST, "paid-pair digest changed")
    if EXPECTED_THRESHOLD_DIGEST is not None:
        require(
            threshold_digest == EXPECTED_THRESHOLD_DIGEST,
            "threshold-pair digest changed",
        )
    if EXPECTED_FORCED_DIGEST is not None:
        require(forced_digest == EXPECTED_FORCED_DIGEST, "forced digest changed")
    if EXPECTED_NORMALIZATION_DIGEST is not None:
        require(
            normalization_digest == EXPECTED_NORMALIZATION_DIGEST,
            "normalization digest changed",
        )

    scalar_failures = [
        row for row in rows if row["scalar_class"] == "failure"
    ]
    positive_scalar_failures = [
        row for row in scalar_failures if row["margin"] > 0
    ]
    w_distribution = tuple(
        (
            cutoff,
            sum(row["cutoff"] <= cutoff for row in positive_scalar_failures),
        )
        for cutoff in (
            500,
            1000,
            2500,
            5000,
            10_000,
            25_000,
            50_000,
            100_000,
        )
    )
    maximum_paid = max(
        (row["paid"], row["body"], row["apex_rank"], row["apex"])
        for row in rows
    )
    maximum_rank1 = max(
        (
            row["rank1_threshold"],
            row["body"],
            row["apex_rank"],
            row["apex"],
        )
        for row in rows
    )
    maximum_entry = max(
        (
            row["pair_entry"],
            row["body"],
            row["apex_rank"],
            row["apex"],
        )
        for row in rows
    )

    hostile_keys = {
        (body, apex)
        for body, _, apex, _, _, _ in EXPECTED_HOSTILES
    }
    root_gate_lines: list[str] = []
    root_gate_base_closed = 0
    root_gate_refined_closed = 0
    active_base_closed = 0
    active_refined_closed = 0
    residual_refined: list[tuple[int, ...]] = []
    for body_index, body in enumerate(BODIES):
        body_rows = rows[10 * body_index : 10 * body_index + 10]
        require(
            len(body_rows) == 10
            and all(row["body"] == body for row in body_rows),
            f"root gate row block changed: {body}",
        )
        root = THM2885.profile_body(body)
        terminal_base = {
            row["apex"]
            for row in body_rows
            if row["scalar_class"] == "direct"
            or row["pairpair_direct"]
        }
        terminal_refined = terminal_base | {
            row["apex"]
            for row in body_rows
            if (body, row["apex"]) in hostile_keys
        }

        def complement_margin(terminals: set[int]):
            allowed = [
                (value, speed)
                for rank, (value, speed) in enumerate(
                    root["top15"],
                    start=1,
                )
                if rank >= 11 or speed not in terminals
            ]
            require(len(allowed) >= 5, f"short root complement: {body}")
            top5 = tuple(allowed[:5])
            return (
                root["m"]
                - sum((value for value, _ in top5), F(0)),
                top5,
            )

        base_margin, base_allowed = complement_margin(terminal_base)
        refined_margin, refined_allowed = complement_margin(terminal_refined)
        base_ok = base_margin > 0
        refined_ok = refined_margin > 0
        require(
            not base_ok or refined_ok,
            f"refined terminal gate regressed: {body}",
        )
        root_gate_base_closed += base_ok
        root_gate_refined_closed += refined_ok
        if max(body) >= 13:
            active_base_closed += base_ok
            active_refined_closed += refined_ok
            if not refined_ok:
                residual_refined.append(body)
        root_gate_lines.append(
            "E="
            + ",".join(map(str, body))
            + ";base="
            + ",".join(map(str, sorted(terminal_base)))
            + ";refined="
            + ",".join(map(str, sorted(terminal_refined)))
            + f";base_margin={ftext(base_margin)};"
            + "base_allowed="
            + ",".join(
                f"{speed}:{ftext(value)}"
                for value, speed in base_allowed
            )
            + f";refined_margin={ftext(refined_margin)};"
            + "refined_allowed="
            + ",".join(
                f"{speed}:{ftext(value)}"
                for value, speed in refined_allowed
            )
            + "\n"
        )
    root_gate_digest = digest_text(
        "THM2888/genuine-terminal-root-weighted-gate/v1\n",
        root_gate_lines,
    )
    if EXPECTED_ROOT_GATE_DIGEST is not None:
        require(
            root_gate_digest == EXPECTED_ROOT_GATE_DIGEST,
            "genuine-terminal root gate digest changed",
        )

    print("THM-2888 EIGHT-BODY FIRST-APEX GLOBAL PAIR-CAP ATLAS")
    print("status=FINITE-EXACT+GLOBAL-PAIR-TAIL-SEALED")
    print(
        "universe=3003 bodies x 10 THM-2885 apices=30030 literal carriers"
    )
    print(
        f"scalar_classes=direct:{scalar_counts[0]},"
        f"rank3_finite:{scalar_counts[1]},failure:{scalar_counts[2]}"
    )
    print(
        f"pair_cap_below_5h7={pair_counts[0]};"
        f"nonpositive={pair_counts[1]};"
        f"pairpair_direct_all={pair_counts[2]};"
        f"pairpair_direct_scalar_failures={pair_counts[3]}"
    )
    print(
        f"exact_pair_max_checks={sum(row['paid'] for row in rows)};"
        f"maximum_per_carrier={maximum_paid[0]};"
        f"maximum_body={maximum_paid[1]};"
        f"maximum_rank={maximum_paid[2]};maximum_apex={maximum_paid[3]};"
        f"threshold_pair_checks={sum(row['threshold_checks'] for row in rows)}"
    )
    print(
        f"maximum_rank1_threshold={ftext(maximum_rank1[0])};"
        f"body={maximum_rank1[1]};rank={maximum_rank1[2]};"
        f"apex={maximum_rank1[3]};horizon={PAIR_HORIZON}"
    )
    print(
        f"maximum_pair_entry_threshold={ftext(maximum_entry[0])};"
        f"body={maximum_entry[1]};rank={maximum_entry[2]};"
        f"apex={maximum_entry[3]};horizon={PAIR_HORIZON}"
    )
    print(
        f"minimum_positive_pair_margin={ftext(minimum_positive[0])};"
        f"body={minimum_positive[1]};rank={minimum_positive[2]};"
        f"apex={minimum_positive[3]};W={minimum_positive[4]};"
        f"maxpair={minimum_positive[5][0]},{minimum_positive[5][1]}"
    )
    print(
        "scalar_failure_W_distribution="
        + ",".join(f"le{cutoff}:{count}" for cutoff, count in w_distribution)
    )
    print(
        "five_hostiles="
        + "|".join(
            ",".join(map(str, row["body"]))
            + f";rank={row['apex_rank']};a={row['apex']};"
            + f"P={row['heavy'][0][0]},{row['heavy'][0][1]};"
            + f"delta={ftext(row['margin'])}"
            for row in hostile
        )
    )
    print(
        "forced_pair_residuals=5/5 globally pair-closed;"
        f"minimum_margin={ftext(min(row['margin'] for row in forced))};"
        f"exact_pair_checks={sum(row['paid'] for row in forced)}"
    )
    maximum_deleted_cutoff = max(
        (
            row["deleted_cutoff"],
            row["body"],
            row["apex_rank"],
            row["apex"],
        )
        for row in normalized
    )
    minimum_normalized_forced = min(
        (
            row["forced"]["margin"],
            row["body"],
            row["apex_rank"],
            row["apex"],
            row["edge"],
        )
        for row in normalized
    )
    print(
        f"heavy_edge_normalizations={len(normalized)};"
        f"positive_high_W={len(high_cutoff)};hostile={len(hostile)};"
        f"max_deleted_W={maximum_deleted_cutoff[0]};"
        f"max_deleted_body={maximum_deleted_cutoff[1]};"
        f"max_deleted_rank={maximum_deleted_cutoff[2]};"
        f"max_deleted_apex={maximum_deleted_cutoff[3]};"
        f"minimum_forced_margin={ftext(minimum_normalized_forced[0])};"
        f"minimum_forced_body={minimum_normalized_forced[1]};"
        f"minimum_forced_apex={minimum_normalized_forced[3]}"
    )
    print(
        f"genuine_terminal_root_gate_all="
        f"base:{root_gate_base_closed}/3003,"
        f"plus_five_refined:{root_gate_refined_closed}/3003;"
        f"active_base:{active_base_closed}/2508;"
        f"active_plus_five_refined:{active_refined_closed}/2508;"
        f"active_residual={len(residual_refined)}"
    )
    print(f"profile_digest_sha256={profile_digest}")
    print(f"paid_pair_digest_sha256={paid_digest}")
    print(f"threshold_pair_digest_sha256={threshold_digest}")
    print(f"forced_residual_digest_sha256={forced_digest}")
    print(f"normalization_digest_sha256={normalization_digest}")
    print(f"root_gate_digest_sha256={root_gate_digest}")
    print("vector_scalar_controls=60060;pair_path_controls=30030")
    print(
        "scope=pair atlas and forced-heavy-edge subbranches only;"
        "5069 non-pair-direct scalar failures retain triangle/literal obligations"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
