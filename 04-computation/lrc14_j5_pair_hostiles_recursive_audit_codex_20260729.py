#!/usr/bin/env python3
r"""Independent pair-threshold audit and recursive probe for the j=5 scout.

This companion deliberately keeps two statements separate.

1.  It replays the full 30,030 first-apex scalar classification from the
    phase-one scout and, only on its 5,122 scalar failures, applies an exact
    pair-threshold branch-and-bound.  If ``C`` has mass ``h`` and global
    largest individual coverage ``q1 < 4h/7``, a pair with union at least
    ``5h/7`` must have both individual coverages at least

        R2 = 5h/7 - q1 > h/7.

    The strict THM-735 cap therefore makes the candidate-speed set finite.
    Every candidate pair is subtracted literally from the interval carrier.

2.  For the resulting exceptional first-apex carriers it independently seals
    the carrier top fourteen, classifies the ten second-apex branches with
    three slots left, tests their weighted complement, and also profiles the
    residual left by each exceptional pair.

The pair threshold is a structural reduction, not by itself a proof that all
literal four-comb survivors are positive.  A carrier with no exceptional pair
has an at-most-one-asymptotic-tail reduction; a carrier with an exceptional
pair additionally needs the explicitly profiled pair branch.
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
PHASE1_PATH = (
    ROOT
    / "04-computation/lrc14_j5_weighted_transversal_recursion_scout_codex_20260729.py"
)
PHASE1_SHA256 = (
    "a942068031d1d389af34af18148996fd8fc346dd4a0a093cf4c3b2acb9119ca1"
)
FIRST_EXTERNAL = 15
BASE_HORIZON = 600
PAIR_HORIZON = 2500
S2 = F(99, 70)

EXPECTED_PAIR_HOSTILES = (
    ((1, 2, 3, 7, 8, 10, 11, 13), 5, 19, (17, 18), F(-50887, 285170886)),
    ((1, 2, 3, 7, 8, 10, 11, 13), 7, 17, (18, 19), F(-3650923, 2851708860)),
    ((1, 4, 6, 7, 9, 10, 11, 13), 6, 23, (16, 17), F(-480313, 920551632)),
    ((2, 3, 4, 8, 10, 12, 13, 14), 6, 16, (18, 22), F(-23, 1070160)),
    ((2, 4, 6, 8, 10, 12, 13, 14), 4, 16, (18, 22), F(-3553, 1070160)),
)
EXPECTED_THRESHOLD_SUMMARY = (
    353,
    (5, (1, 4, 6, 8, 9, 10, 13, 14), 3, 23),
    (457, (1, 2, 6, 8, 9, 10, 13, 14), 6, 46, F(1559025468, 3415457)),
    (6, (1, 4, 6, 8, 9, 10, 13, 14), 3, 23),
)
EXPECTED_CAP_SUMMARY = (
    35699,
    (765, (2, 3, 5, 7, 8, 10, 11, 13), 5, 46),
    53,
    (
        130828,
        (2, 4, 6, 7, 8, 10, 12, 14),
        3,
        26,
        F(709, 8828820),
        (18, 22),
    ),
)
EXPECTED_PAIR_DIGEST = (
    "4ed0d61f76fa20ebb6fe6876747146c36405e6a81a53581c57e1caec4c77ab1a"
)
EXPECTED_FORCED_PAIR_MARGINS = (
    F(2597104, 212952285),
    F(2597104, 212952285),
    F(389293, 63488880),
    F(15287, 840840),
    F(356, 35035),
)
EXPECTED_AUDIT_DIGEST = (
    "47bef0aa54683e857abc1bff25f66281689ac142220ea7878b592037ea7c023e"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_phase1():
    require(file_sha256(PHASE1_PATH) == PHASE1_SHA256, "phase-one script changed")
    spec = importlib.util.spec_from_file_location("j5_pair_audit_phase1", PHASE1_PATH)
    require(spec is not None and spec.loader is not None, "cannot load phase one")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


P1 = load_phase1()
CORE = P1.CORE


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def digest_rows(header: str, rows: list[str]) -> str:
    return hashlib.sha256((header + "".join(rows)).encode()).hexdigest()


def reconstruct_carrier(
    body: tuple[int, ...],
    prefix: tuple[int, ...],
) -> list[tuple[F, F]]:
    good, _, _ = CORE.good_norm(body)
    carrier = P1.THM2883.subtract_local_multi(good, prefix)
    family = tuple(sorted((*body, *prefix)))
    direct_good, direct_r, direct_m = CORE.good_norm(family)
    mass = sum((right - left for left, right in carrier), F(0))
    require(
        direct_good == carrier
        and direct_r == len(carrier)
        and direct_m == mass,
        f"recursive/full carrier mismatch: {body}, {prefix}",
    )
    return carrier


def seal_top(
    body: tuple[int, ...],
    prefix: tuple[int, ...],
    carrier: list[tuple[F, F]],
    rank_count: int,
) -> dict[str, object]:
    """Globally seal the requested leading scalar coverage ranks."""

    mass = sum((right - left for left, right in carrier), F(0))
    components = len(carrier)
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, BASE_HORIZON + 1)
        if speed not in prefix
    ]
    rows = P1.THM2885.coverages_many(carrier, speeds)
    ranked600 = sorted(rows, key=lambda item: (-item[0], item[1]))
    q_base = ranked600[rank_count - 1][0]
    require(
        q_base > mass / 7,
        f"rank-{rank_count} base value misses limit: {body}, {prefix}",
    )
    threshold = S2 * components / (7 * (q_base - mass / 7))
    tail_first = ceiling(threshold)
    tail_speeds = [
        speed
        for speed in range(BASE_HORIZON + 1, tail_first)
        if speed not in prefix
    ]
    rows.extend(P1.THM2885.coverages_many(carrier, tail_speeds))
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    require(
        mass / 7 + S2 * components / (7 * tail_first) <= q_base,
        f"rank-{rank_count} strict tail seal failed: {body}, {prefix}",
    )
    controls = [speeds[0], tail_speeds[-1] if tail_speeds else speeds[8]]
    by_speed = {speed: value for value, speed in rows}
    for speed in controls:
        require(
            by_speed[speed] == P1.THM2885.coverage(carrier, speed),
            f"top-rank vector/scalar mismatch: {body}, {prefix}, {speed}",
        )
    return {
        "m": mass,
        "r": components,
        "top": tuple(ranked[:rank_count]),
        "threshold": threshold,
        "tail_first": tail_first,
    }


def scalar_failure_body(body: tuple[int, ...]) -> list[dict[str, object]]:
    root = P1.root_profile(body)
    rows = P1.profile_first_apices(root)
    return [row for row in rows if row["classification"] == "failure"]


def global_pair_cap(
    body: tuple[int, ...],
    prefix: tuple[int, ...],
    carrier: list[tuple[F, F]],
    excluded_pairs: frozenset[tuple[int, int]] = frozenset(),
) -> dict[str, object]:
    """Bound every pair union by exact head search plus a strict tail cap."""

    mass = sum((right - left for left, right in carrier), F(0))
    components = len(carrier)
    sealed1 = seal_top(body, prefix, carrier, 1)
    q1, q1_speed = sealed1["top"][0]
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, PAIR_HORIZON + 1)
        if speed not in prefix
    ]
    rows = P1.THM2885.coverages_many(carrier, speeds)
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    finite_cap = F(0)
    finite_pair = None
    paid = 0
    for first in range(len(ranked) - 1):
        if ranked[first][0] + ranked[first + 1][0] < finite_cap:
            break
        x = ranked[first][1]
        after_x = P1.THM2883.subtract_local(carrier, x)
        for second in range(first + 1, len(ranked)):
            if ranked[first][0] + ranked[second][0] < finite_cap:
                break
            y = ranked[second][1]
            pair = tuple(sorted((x, y)))
            if pair in excluded_pairs:
                continue
            survivor_carrier = P1.THM2883.subtract_local(after_x, y)
            survivor = sum(
                (right - left for left, right in survivor_carrier),
                F(0),
            )
            union = mass - survivor
            paid += 1
            if (
                union > finite_cap
                or (
                    union == finite_cap
                    and finite_pair is not None
                    and pair < finite_pair
                )
            ):
                finite_cap = union
                finite_pair = pair
    tail_cap = (
        q1
        + mass / 7
        + S2 * components / (7 * (PAIR_HORIZON + 1))
    )
    cap = max(finite_cap, tail_cap)
    return {
        "m": mass,
        "r": components,
        "q1": q1,
        "q1_speed": q1_speed,
        "q1_tail_first": sealed1["tail_first"],
        "finite_cap": finite_cap,
        "finite_pair": finite_pair,
        "tail_cap": tail_cap,
        "cap": cap,
        "paid": paid,
    }


def pair_threshold_profile(row: dict[str, object]) -> dict[str, object]:
    body = row["body"]
    apex = row["prefix"][0]
    apex_rank = row["apex_rank"]
    require(
        isinstance(body, tuple)
        and isinstance(apex, int)
        and isinstance(apex_rank, int),
        "scalar failure types changed",
    )
    carrier = reconstruct_carrier(body, (apex,))
    sealed1 = seal_top(body, (apex,), carrier, 1)
    mass = sealed1["m"]
    components = sealed1["r"]
    q1, q1_speed = sealed1["top"][0]
    require(
        isinstance(mass, F)
        and isinstance(components, int)
        and isinstance(q1, F),
        "pair carrier types changed",
    )
    residual = F(5, 7) * mass - q1
    require(
        residual > mass / 7,
        f"pair-threshold residual is not finite: {body}, {apex}",
    )
    threshold = S2 * components / (7 * (residual - mass / 7))
    tail_first = ceiling(threshold)
    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, tail_first)
        if speed != apex
    ]
    coverage_rows = P1.THM2885.coverages_many(carrier, speeds)
    head = tuple(
        sorted(
            (
                (value, speed)
                for value, speed in coverage_rows
                if value >= residual
            ),
            key=lambda item: (-item[0], item[1]),
        )
    )
    if speeds:
        controls = [speeds[0], speeds[-1]]
        by_speed = {speed: value for value, speed in coverage_rows}
        for speed in controls:
            require(
                by_speed[speed] == P1.THM2885.coverage(carrier, speed),
                f"pair-head vector/scalar mismatch: {body}, {apex}, {speed}",
            )
    crossing = []
    paid = 0
    threshold_union = F(5, 7) * mass
    for first in range(len(head) - 1):
        if head[first][0] + head[first + 1][0] < threshold_union:
            break
        x = head[first][1]
        after_x = P1.THM2883.subtract_local(carrier, x)
        for second in range(first + 1, len(head)):
            if head[first][0] + head[second][0] < threshold_union:
                break
            y = head[second][1]
            survivor_carrier = P1.THM2883.subtract_local(after_x, y)
            survivor = sum(
                (right - left for left, right in survivor_carrier),
                F(0),
            )
            union = mass - survivor
            paid += 1
            if union >= threshold_union:
                simultaneous = P1.THM2883.subtract_local_multi(carrier, (x, y))
                require(
                    simultaneous == survivor_carrier,
                    f"pair cached/simultaneous mismatch: {body}, {apex}, {x}, {y}",
                )
                crossing.append(
                    (
                        tuple(sorted((x, y))),
                        union,
                        threshold_union - union,
                        survivor,
                    )
                )
    require(
        len(crossing) == len({entry[0] for entry in crossing}),
        f"duplicate exceptional pair: {body}, {apex}",
    )
    cap_row = global_pair_cap(body, (apex,), carrier)
    crossing_pairs = frozenset(entry[0] for entry in crossing)
    nonexception_row = (
        global_pair_cap(body, (apex,), carrier, crossing_pairs)
        if crossing_pairs
        else cap_row
    )
    threshold_union = F(5, 7) * mass
    require(
        nonexception_row["cap"] < threshold_union,
        f"nonexceptional pair cap failed: {body}, {apex}",
    )
    pair_delta = threshold_union - nonexception_row["cap"]
    one_tail_threshold = 2 * S2 * components / (7 * pair_delta)
    return {
        "body": body,
        "apex_rank": apex_rank,
        "apex": apex,
        "m": mass,
        "r": components,
        "q1": q1,
        "q1_speed": q1_speed,
        "q1_tail_first": sealed1["tail_first"],
        "pair_residual": residual,
        "pair_threshold": threshold,
        "pair_tail_first": tail_first,
        "head_size": len(head),
        "paid": paid,
        "crossing": tuple(crossing),
        "pair_cap": cap_row["cap"],
        "pair_cap_finite": cap_row["finite_cap"],
        "pair_cap_finite_pair": cap_row["finite_pair"],
        "pair_cap_tail": cap_row["tail_cap"],
        "pair_cap_paid": cap_row["paid"],
        "nonexception_cap": nonexception_row["cap"],
        "nonexception_pair": nonexception_row["finite_pair"],
        "nonexception_paid": nonexception_row["paid"],
        "pair_delta": pair_delta,
        "one_tail_threshold": one_tail_threshold,
        "one_tail_first": ceiling(one_tail_threshold),
        "two_pair_direct": 2 * cap_row["cap"] < mass,
    }


def audit_hostile(row: dict[str, object]) -> dict[str, object]:
    body = row["body"]
    apex = row["apex"]
    apex_rank = row["apex_rank"]
    crossing = row["crossing"]
    require(len(crossing) == 1, f"hostile pair is not unique: {body}, {apex}")
    pair = crossing[0][0]
    first_carrier = reconstruct_carrier(body, (apex,))
    top14 = seal_top(body, (apex,), first_carrier, 14)
    raw_hitting_margin = top14["m"] - sum(
        (value for value, _ in top14["top"][10:14]),
        F(0),
    )

    second_rows = []
    for second_rank, (_, second) in enumerate(top14["top"][:10], start=1):
        second_carrier = P1.THM2883.subtract_local(first_carrier, second)
        profile = P1.profile_carrier(
            body,
            (apex, second),
            second_carrier,
            3,
        )
        profile["second_rank"] = second_rank
        second_rows.append(profile)
    second_counts = {
        name: sum(profile["classification"] == name for profile in second_rows)
        for name in ("direct", "rank_feasible", "failure")
    }
    failed_second = {
        profile["prefix"][1]
        for profile in second_rows
        if profile["classification"] == "failure"
    }
    allowed = [
        (value, speed)
        for rank, (value, speed) in enumerate(top14["top"], start=1)
        if rank >= 11 or speed in failed_second
    ]
    require(len(allowed) >= 4, f"short second-level complement: {body}, {apex}")
    allowed_top4 = tuple(allowed[:4])
    weighted_margin = top14["m"] - sum(
        (value for value, _ in allowed_top4),
        F(0),
    )

    pair_prefix = (apex, *pair)
    pair_carrier = reconstruct_carrier(body, pair_prefix)
    pair_profile = P1.profile_carrier(body, pair_prefix, pair_carrier, 2)
    pair_union_profile = global_pair_cap(body, pair_prefix, pair_carrier)
    require(
        pair_union_profile["cap"] < pair_union_profile["m"],
        f"forced-pair residual union cap failed: {body}, {pair_prefix}",
    )
    return {
        "body": body,
        "apex_rank": apex_rank,
        "apex": apex,
        "pair": pair,
        "pair_margin": crossing[0][2],
        "pair_survivor": crossing[0][3],
        "top14": top14,
        "raw_hitting_margin": raw_hitting_margin,
        "second_rows": tuple(second_rows),
        "second_counts": second_counts,
        "failed_second": tuple(sorted(failed_second)),
        "allowed_top4": allowed_top4,
        "weighted_margin": weighted_margin,
        "pair_profile": pair_profile,
        "pair_union_profile": pair_union_profile,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    context = mp.get_context("spawn")
    if args.workers == 1:
        nested = list(map(scalar_failure_body, P1.BODIES))
    else:
        with context.Pool(args.workers) as pool:
            nested = pool.map(scalar_failure_body, P1.BODIES, chunksize=1)
    failures = [row for body_rows in nested for row in body_rows]
    require(len(failures) == 5122, "phase-one scalar failure census changed")
    if args.workers == 1:
        pair_rows = list(map(pair_threshold_profile, failures))
    else:
        with context.Pool(args.workers) as pool:
            pair_rows = pool.map(pair_threshold_profile, failures, chunksize=1)

    hostiles = [row for row in pair_rows if row["crossing"]]
    pair_rows_text = [
        "E="
        + ",".join(map(str, row["body"]))
        + f";rank={row['apex_rank']};a={row['apex']};"
        + f"m={ftext(row['m'])};r={row['r']};"
        + f"q1={row['q1_speed']}:{ftext(row['q1'])};"
        + f"T1={row['q1_tail_first']};"
        + f"R2={ftext(row['pair_residual'])};"
        + f"T2={ftext(row['pair_threshold'])};tail2={row['pair_tail_first']};"
        + f"K2={row['head_size']};paid={row['paid']};cross="
        + ",".join(
            f"{pair[0]}-{pair[1]}:{ftext(margin)}"
            for pair, _, margin, _ in row["crossing"]
        )
        + f";U2={ftext(row['pair_cap'])};"
        + f"U2pair={row['pair_cap_finite_pair']};"
        + f"U2tail={ftext(row['pair_cap_tail'])};"
        + f"U2paid={row['pair_cap_paid']};"
        + f"U2nonexception={ftext(row['nonexception_cap'])};"
        + f"delta={ftext(row['pair_delta'])};"
        + f"W={row['one_tail_first']};"
        + f"twopair={int(row['two_pair_direct'])}"
        + "\n"
        for row in pair_rows
    ]
    pair_digest = digest_rows(
        "LRC14/j5/scalar-failure-pair-threshold/v1\n",
        pair_rows_text,
    )
    maximum_threshold = max(
        (
            row["pair_tail_first"],
            row["body"],
            row["apex_rank"],
            row["apex"],
            row["pair_threshold"],
        )
        for row in pair_rows
    )
    maximum_head = max(
        (
            row["head_size"],
            row["body"],
            row["apex_rank"],
            row["apex"],
        )
        for row in pair_rows
    )
    maximum_paid = max(
        (
            row["paid"],
            row["body"],
            row["apex_rank"],
            row["apex"],
        )
        for row in pair_rows
    )
    maximum_cap_paid = max(
        (
            row["pair_cap_paid"],
            row["body"],
            row["apex_rank"],
            row["apex"],
        )
        for row in pair_rows
    )
    maximum_one_tail = max(
        (
            row["one_tail_first"],
            row["body"],
            row["apex_rank"],
            row["apex"],
            row["pair_delta"],
            row["nonexception_pair"],
        )
        for row in pair_rows
    )
    minimum_delta = min(
        (
            row["pair_delta"],
            row["body"],
            row["apex_rank"],
            row["apex"],
            row["nonexception_pair"],
        )
        for row in pair_rows
    )
    hostile_signature = tuple(
        (
            row["body"],
            row["apex_rank"],
            row["apex"],
            row["crossing"][0][0],
            row["crossing"][0][2],
        )
        for row in hostiles
    )
    require(
        hostile_signature == EXPECTED_PAIR_HOSTILES,
        "five pair-hostile identities changed",
    )
    require(
        (
            sum(row["paid"] for row in pair_rows),
            maximum_paid,
            maximum_threshold,
            maximum_head,
        )
        == EXPECTED_THRESHOLD_SUMMARY,
        "pair-threshold summary changed",
    )
    require(
        (
            sum(row["pair_cap_paid"] for row in pair_rows),
            maximum_cap_paid,
            sum(row["two_pair_direct"] for row in pair_rows),
            maximum_one_tail,
        )
        == EXPECTED_CAP_SUMMARY,
        "global pair-cap summary changed",
    )
    require(
        minimum_delta[:4]
        == (
            F(709, 8828820),
            (2, 4, 6, 7, 8, 10, 12, 14),
            3,
            26,
        )
        and pair_digest == EXPECTED_PAIR_DIGEST,
        "pair-cap extremum or ledger changed",
    )

    print("LRC14 J=5 PAIR-HOSTILE RECURSIVE AUDIT")
    print(
        "status=FINITE-EXACT-SCOUT;pair_reduction_is_structural;"
        "not_literal_four_comb_positivity"
    )
    print(
        f"scalar_failures={len(pair_rows)};"
        f"pair_finite_reduced={len(pair_rows)-len(hostiles)};"
        f"pair_hostiles={len(hostiles)};"
        f"threshold_paid_pairs={sum(row['paid'] for row in pair_rows)};"
        f"max_paid={maximum_paid[0]};max_paid_body={maximum_paid[1]};"
        f"max_paid_rank={maximum_paid[2]};max_paid_apex={maximum_paid[3]};"
        f"pair_digest={pair_digest}"
    )
    print(
        f"max_pair_tail={maximum_threshold[0]};"
        f"body={maximum_threshold[1]};rank={maximum_threshold[2]};"
        f"apex={maximum_threshold[3]};"
        f"threshold={ftext(maximum_threshold[4])};"
        f"max_pair_head={maximum_head[0]};"
        f"max_head_body={maximum_head[1]};"
        f"max_head_rank={maximum_head[2]};max_head_apex={maximum_head[3]}"
    )
    print(
        f"global_pair_cap_paid={sum(row['pair_cap_paid'] for row in pair_rows)};"
        f"max_cap_paid={maximum_cap_paid[0]};"
        f"max_cap_paid_body={maximum_cap_paid[1]};"
        f"max_cap_paid_rank={maximum_cap_paid[2]};"
        f"max_cap_paid_apex={maximum_cap_paid[3]};"
        f"two_pair_direct={sum(row['two_pair_direct'] for row in pair_rows)};"
        f"one_tail_reductions={len(pair_rows)};"
        f"max_one_tail_first={maximum_one_tail[0]};"
        f"max_one_tail_body={maximum_one_tail[1]};"
        f"max_one_tail_rank={maximum_one_tail[2]};"
        f"max_one_tail_apex={maximum_one_tail[3]};"
        f"max_one_tail_delta={ftext(maximum_one_tail[4])};"
        f"max_one_tail_pair={maximum_one_tail[5]};"
        f"minimum_delta={ftext(minimum_delta[0])};"
        f"minimum_delta_body={minimum_delta[1]};"
        f"minimum_delta_rank={minimum_delta[2]};"
        f"minimum_delta_apex={minimum_delta[3]};"
        f"minimum_delta_pair={minimum_delta[4]}"
    )
    for row in hostiles:
        crossing_text = ",".join(
            f"{pair}:union={ftext(union)}:margin={ftext(margin)}"
            for pair, union, margin, _ in row["crossing"]
        )
        print(
            f"pair_hostile=E{row['body']};rank={row['apex_rank']};"
            f"apex={row['apex']};{crossing_text}"
        )

    audits = [audit_hostile(row) for row in hostiles]
    audit_rows = []
    for row in audits:
        pair_profile = row["pair_profile"]
        pair_union_profile = row["pair_union_profile"]
        audit_rows.append(
            "E="
            + ",".join(map(str, row["body"]))
            + f";rank={row['apex_rank']};a={row['apex']};"
            + f"pair={row['pair'][0]},{row['pair'][1]};"
            + f"pair_margin={ftext(row['pair_margin'])};"
            + f"top14_tail={row['top14']['tail_first']};"
            + f"raw_hit={ftext(row['raw_hitting_margin'])};"
            + f"second={row['second_counts']};failed="
            + ",".join(map(str, row["failed_second"]))
            + ";allowed="
            + ",".join(str(speed) for _, speed in row["allowed_top4"])
            + f";weighted={ftext(row['weighted_margin'])};"
            + f"pair_class={pair_profile['classification']};"
            + f"pair_m={ftext(pair_profile['m'])};"
            + "pair_top2="
            + ",".join(
                f"{speed}:{ftext(value)}"
                for value, speed in pair_profile["topk"]
            )
            + f";pair_direct={ftext(pair_profile['direct_margin'])};"
            + f"pair_R1={ftext(pair_profile['residual'])};"
            + (
                f"pair_head_tail={pair_profile['head_tail_first']}"
                if pair_profile["head_tail_first"] is not None
                else "pair_head_tail=NA"
            )
            + f";pair_U2={ftext(pair_union_profile['cap'])};"
            + f"pair_U2_pair={pair_union_profile['finite_pair']};"
            + f"pair_U2_tail={ftext(pair_union_profile['tail_cap'])};"
            + f"pair_U2_margin="
            + f"{ftext(pair_union_profile['m']-pair_union_profile['cap'])}"
            + "\n"
        )
        print(
            f"recursive_hostile=E{row['body']};rank={row['apex_rank']};"
            f"apex={row['apex']};pair={row['pair']};"
            f"top14_tail={row['top14']['tail_first']};"
            f"raw_top10_hit_margin={ftext(row['raw_hitting_margin'])};"
            f"second_direct={row['second_counts']['direct']};"
            f"second_rank2={row['second_counts']['rank_feasible']};"
            f"second_failed={row['second_counts']['failure']};"
            f"weighted_margin={ftext(row['weighted_margin'])};"
            f"pair_residual_class={pair_profile['classification']};"
            f"pair_residual_margin={ftext(pair_profile['direct_margin'])};"
            f"pair_residual_R1_minus_m7="
            f"{ftext(pair_profile['residual']-pair_profile['m']/7)};"
            f"pair_union_cap_margin="
            f"{ftext(pair_union_profile['m']-pair_union_profile['cap'])};"
            f"pair_union_maximizer={pair_union_profile['finite_pair']}"
        )
    audit_digest = digest_rows(
        "LRC14/j5/five-pair-hostile-recursion/v1\n",
        audit_rows,
    )
    forced_pair_margins = tuple(
        row["pair_union_profile"]["m"] - row["pair_union_profile"]["cap"]
        for row in audits
    )
    require(
        forced_pair_margins == EXPECTED_FORCED_PAIR_MARGINS
        and audit_digest == EXPECTED_AUDIT_DIGEST
        and all(row["raw_hitting_margin"] > 0 for row in audits)
        and all(row["pair_union_profile"]["cap"] < row["pair_union_profile"]["m"] for row in audits),
        "five-hostile recursive audit changed",
    )
    print(
        f"recursive_audits={len(audits)};"
        f"raw_top10_hitting_positive="
        f"{sum(row['raw_hitting_margin']>0 for row in audits)};"
        f"second_weighted_positive="
        f"{sum(row['weighted_margin']>0 for row in audits)};"
        f"pair_residual_direct="
        f"{sum(row['pair_profile']['classification']=='direct' for row in audits)};"
        f"pair_residual_rank1="
        f"{sum(row['pair_profile']['classification']=='rank_feasible' for row in audits)};"
        f"pair_residual_failed="
        f"{sum(row['pair_profile']['classification']=='failure' for row in audits)};"
        f"forced_pair_union_closed="
        f"{sum(row['pair_union_profile']['cap']<row['pair_union_profile']['m'] for row in audits)};"
        f"audit_digest={audit_digest}"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
