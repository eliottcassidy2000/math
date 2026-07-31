#!/usr/bin/env python3
r"""Exact structural scout for the eight-body / five-tail LRC(14) rung.

This is a research artifact, not a canonical theorem companion.  It starts
from THM-2885's global top-ten hitting gate for every
``E in C({1,...,14},8)``.  For each of the 30,030 selected first apices it
forms the literal carrier

    C(E;a) = G_E \ D_a

and globally seals its four largest remaining individual tooth coverages.
The carrier is then classified as

* ``direct``: the four largest coverages sum to less than ``|C|``;
* ``rank3``: direct closure fails, but the residual after the three largest
  coverages is strictly greater than ``|C|/7`` and hence has a finite exact
  fourth-speed head by the strict THM-735 discrepancy cap;
* ``failure``: neither scalar certificate applies.

The first weighted-complement gate treats ``direct`` and ``rank3`` as
structural terminals (the latter is a finite reduction, not a positivity
proof).  It asks whether a dangerous five-set can avoid every terminal apex.
The script reports the resulting exact residual and is intended to recurse
only there.

THM-885 separately closes the 495 bodies contained in ``{1,...,12}``, because
adding five speeds at least 15 gives exactly five outliers relative to
``{1,...,12}``.  We nevertheless classify all 30,030 apices so that the
arithmetic census and the inherited shortcut are visible separately.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
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
THM885_PATH = (
    ROOT / "01-canon/theorems/THM-885-covering-case-decomposition-j56-sweep.md"
)
THM885_SHA256 = (
    "3d86a3610cfd2c366833eb991226efa378e06d8a1458c1fd9730b0b37cd4f56b"
)

FIRST_EXTERNAL = 15
BASE_HORIZON = 600
S2 = F(99, 70)
BODIES = tuple(combinations(range(1, 15), 8))


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


THM2885 = load_module(THM2885_PATH, THM2885_SHA256, "j5_scout_thm2885")
THM2883 = load_module(THM2883_PATH, THM2883_SHA256, "j5_scout_thm2883")
require(file_sha256(THM885_PATH) == THM885_SHA256, "THM-885 source changed")
CORE = THM2885.CORE


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def digest_rows(header: str, rows: list[str]) -> str:
    return hashlib.sha256((header + "".join(rows)).encode()).hexdigest()


def root_profile(body: tuple[int, ...]) -> dict[str, object]:
    row = THM2885.profile_body(body)
    good, components, mass = CORE.good_norm(body)
    require(
        row["m"] == mass and row["r"] == components,
        f"THM-2885/root reconstruction mismatch: {body}",
    )
    return {
        "body": body,
        "good": good,
        "m": mass,
        "r": components,
        "top15": row["top15"],
        "threshold15": row["threshold"],
        "tail_first15": row["tail_first"],
    }


def profile_carrier(
    body: tuple[int, ...],
    prefix: tuple[int, ...],
    carrier: list[tuple[F, F]],
    remaining_slots: int,
) -> dict[str, object]:
    """Globally seal the leading ``remaining_slots`` scalar coverages."""

    require(
        2 <= remaining_slots <= 5
        and len(prefix) + remaining_slots == 5
        and len(set(prefix)) == len(prefix),
        f"bad recursion state: {body}, {prefix}, {remaining_slots}",
    )
    mass = sum((right - left for left, right in carrier), F(0))
    components = len(carrier)
    require(mass > 0 and components > 0, f"empty carrier: {body}, {prefix}")
    base_speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, BASE_HORIZON + 1)
        if speed not in prefix
    ]
    rows = THM2885.coverages_many(carrier, base_speeds)
    by_speed = {speed: value for value, speed in rows}
    require(
        by_speed[base_speeds[0]] == THM2885.coverage(carrier, base_speeds[0]),
        f"base vector/scalar mismatch: {body}, {prefix}",
    )
    ranked600 = sorted(rows, key=lambda item: (-item[0], item[1]))
    qk600 = ranked600[remaining_slots - 1][0]
    require(
        qk600 > mass / 7,
        f"rank-{remaining_slots} base value misses limit: {body}, {prefix}",
    )
    threshold_k = S2 * components / (7 * (qk600 - mass / 7))
    tail_first_k = ceiling(threshold_k)
    tail_speeds = [
        speed
        for speed in range(BASE_HORIZON + 1, tail_first_k)
        if speed not in prefix
    ]
    rows.extend(THM2885.coverages_many(carrier, tail_speeds))
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    qk = ranked[remaining_slots - 1][0]
    require(
        mass / 7 + S2 * components / (7 * tail_first_k) <= qk600,
        f"rank-{remaining_slots} strict tail seal failed: {body}, {prefix}",
    )
    control_speed = (
        tail_speeds[-1] if tail_speeds else base_speeds[min(8, len(base_speeds) - 1)]
    )
    by_speed = {speed: value for value, speed in rows}
    require(
        by_speed[control_speed] == THM2885.coverage(carrier, control_speed),
        f"tail vector/scalar mismatch: {body}, {prefix}",
    )

    direct_margin = mass - sum(
        (value for value, _ in ranked[:remaining_slots]),
        F(0),
    )
    residual = mass - sum(
        (value for value, _ in ranked[: remaining_slots - 1]),
        F(0),
    )
    classification = (
        "direct"
        if direct_margin > 0
        else ("rank_feasible" if residual > mass / 7 else "failure")
    )
    head_threshold = (
        S2 * components / (7 * (residual - mass / 7))
        if classification == "rank_feasible"
        else None
    )
    return {
        "body": body,
        "prefix": prefix,
        "m": mass,
        "r": components,
        "remaining_slots": remaining_slots,
        "topk": tuple(ranked[:remaining_slots]),
        "direct_margin": direct_margin,
        "residual": residual,
        "classification": classification,
        "threshold_k": threshold_k,
        "tail_first_k": tail_first_k,
        "head_threshold": head_threshold,
        "head_tail_first": ceiling(head_threshold)
        if head_threshold is not None
        else None,
    }


def profile_first_apices(root: dict[str, object]) -> list[dict[str, object]]:
    body = root["body"]
    good = root["good"]
    top15 = root["top15"]
    require(
        isinstance(body, tuple)
        and isinstance(good, list)
        and isinstance(top15, tuple),
        "root profile types changed",
    )
    out = []
    for apex_rank, (_, apex) in enumerate(top15[:10], start=1):
        carrier = THM2883.subtract_local(good, apex)
        direct_r, direct_m, direct_good = CORE.subtract(good, apex)
        require(
            direct_good == carrier
            and direct_r == len(carrier)
            and direct_m
            == sum((right - left for left, right in carrier), F(0)),
            f"literal/full first-apex mismatch: {body}, {apex}",
        )
        row = profile_carrier(body, (apex,), carrier, 4)
        row["apex_rank"] = apex_rank
        out.append(row)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(
        len(BODIES) == 3003
        and len(set(BODIES)) == 3003
        and sum(max(body) <= 12 for body in BODIES) == 495,
        "eight-body universe changed",
    )

    context = mp.get_context("spawn")
    if args.workers == 1:
        roots = list(map(root_profile, BODIES))
    else:
        with context.Pool(args.workers) as pool:
            roots = pool.map(root_profile, BODIES, chunksize=2)
    require(
        tuple(row["body"] for row in roots) == BODIES,
        "parallel root order changed",
    )
    print("LRC14 J=5 WEIGHTED-TRANSVERSAL STRUCTURAL SCOUT", flush=True)
    print(
        "status=FINITE-EXACT-SCOUT;not_a_positivity_proof;"
        "rank_feasible_means_finite_head_only",
        flush=True,
    )
    print(
        "root_universe=C(14,8)=3003;"
        "thm885_terminal_bodies=C(12,8)=495;"
        "active_bodies_meeting_13_or_14=2508",
        flush=True,
    )

    if args.workers == 1:
        nested = list(map(profile_first_apices, roots))
    else:
        with context.Pool(args.workers) as pool:
            nested = pool.map(profile_first_apices, roots, chunksize=1)
    carriers = [row for group in nested for row in group]
    require(len(carriers) == 30030, "first-apex task universe changed")
    counts = {
        name: sum(row["classification"] == name for row in carriers)
        for name in ("direct", "rank_feasible", "failure")
    }
    active_carriers = [row for row in carriers if max(row["body"]) >= 13]
    active_counts = {
        name: sum(row["classification"] == name for row in active_carriers)
        for name in ("direct", "rank_feasible", "failure")
    }
    classification_rows = [
        "E="
        + ",".join(map(str, row["body"]))
        + f";rank={row['apex_rank']};a={row['prefix'][0]};"
        + f"class={row['classification']};m={ftext(row['m'])};r={row['r']};"
        + "top4="
        + ",".join(
            f"{speed}:{ftext(value)}" for value, speed in row["topk"]
        )
        + f";direct={ftext(row['direct_margin'])};"
        + f"R3={ftext(row['residual'])};"
        + f"T4={ftext(row['threshold_k'])};tail4={row['tail_first_k']};"
        + (
            f"H3={ftext(row['head_threshold'])};"
            f"head_tail={row['head_tail_first']}"
            if row["head_threshold"] is not None
            else "H3=NA;head_tail=NA"
        )
        + "\n"
        for row in carriers
    ]
    classification_digest = digest_rows(
        "LRC14/j5/first-apex-scalar-classification/v1\n",
        classification_rows,
    )
    maximum_top4_tail = max(
        (
            row["tail_first_k"],
            row["body"],
            row["apex_rank"],
            row["prefix"][0],
            row["threshold_k"],
        )
        for row in carriers
    )
    feasible = [row for row in carriers if row["classification"] == "rank_feasible"]
    maximum_head_tail = max(
        (
            row["head_tail_first"],
            row["body"],
            row["apex_rank"],
            row["prefix"][0],
            row["head_threshold"],
        )
        for row in feasible
    )
    direct_equalities = sum(row["direct_margin"] == 0 for row in carriers)
    residual_equalities = sum(row["residual"] == row["m"] / 7 for row in carriers)
    print(
        f"first_apices={len(carriers)};"
        f"direct={counts['direct']};rank3_feasible={counts['rank_feasible']};"
        f"failed={counts['failure']};"
        f"direct_equalities={direct_equalities};"
        f"R3_equalities={residual_equalities};"
        f"classification_digest={classification_digest}",
        flush=True,
    )
    print(
        f"active_first_apices={len(active_carriers)};"
        f"direct={active_counts['direct']};"
        f"rank3_feasible={active_counts['rank_feasible']};"
        f"failed={active_counts['failure']}",
        flush=True,
    )
    print(
        f"max_top4_tail={maximum_top4_tail[0]};"
        f"body={maximum_top4_tail[1]};rank={maximum_top4_tail[2]};"
        f"apex={maximum_top4_tail[3]};"
        f"threshold={ftext(maximum_top4_tail[4])}",
        flush=True,
    )
    print(
        f"max_rank3_head_tail={maximum_head_tail[0]};"
        f"body={maximum_head_tail[1]};rank={maximum_head_tail[2]};"
        f"apex={maximum_head_tail[3]};"
        f"threshold={ftext(maximum_head_tail[4])}",
        flush=True,
    )

    by_key = {
        (row["body"], row["prefix"][0]): row
        for row in carriers
    }
    gate_rows = []
    for root in roots:
        body = root["body"]
        if max(body) <= 12:
            continue
        failed = {
            row["prefix"][0]
            for row in carriers
            if row["body"] == body and row["classification"] == "failure"
        }
        allowed = [
            (value, speed)
            for rank, (value, speed) in enumerate(root["top15"], start=1)
            if rank >= 11 or speed in failed
        ]
        require(len(allowed) >= 5, f"short first complement: {body}")
        allowed_top5 = tuple(allowed[:5])
        margin = root["m"] - sum(
            (value for value, _ in allowed_top5),
            F(0),
        )
        gate_rows.append(
            {
                "body": body,
                "failed": tuple(sorted(failed)),
                "allowed_top5": allowed_top5,
                "margin": margin,
                "closed": margin > 0,
            }
        )
    require(len(gate_rows) == 2508, "active root gate universe changed")
    residual_roots = [row for row in gate_rows if not row["closed"]]
    minimum_gate = min(
        (row["margin"], row["body"], row["allowed_top5"])
        for row in gate_rows
    )
    gate_text = [
        "E="
        + ",".join(map(str, row["body"]))
        + ";failed="
        + ",".join(map(str, row["failed"]))
        + ";allowed="
        + ",".join(
            f"{speed}:{ftext(value)}" for value, speed in row["allowed_top5"]
        )
        + f";margin={ftext(row['margin'])};closed={int(row['closed'])}\n"
        for row in gate_rows
    ]
    gate_digest = digest_rows(
        "LRC14/j5/first-weighted-complement-gate/v1\n",
        gate_text,
    )
    residual_digest = digest_rows(
        "LRC14/j5/first-weighted-complement-residual/v1\n",
        [
            ",".join(map(str, row["body"])) + "\n"
            for row in residual_roots
        ],
    )
    print(
        f"first_weighted_gate_roots={len(gate_rows)};"
        f"structurally_terminal={len(gate_rows)-len(residual_roots)};"
        f"residual={len(residual_roots)};"
        f"equalities={sum(row['margin']==0 for row in gate_rows)};"
        f"minimum_margin={ftext(minimum_gate[0])};"
        f"minimum_body={minimum_gate[1]};"
        "minimum_allowed="
        + ",".join(str(speed) for _, speed in minimum_gate[2])
        + f";gate_digest={gate_digest};residual_digest={residual_digest}",
        flush=True,
    )
    print("phase1_exact_controls=PASS", flush=True)


if __name__ == "__main__":
    main()
