#!/usr/bin/env python3
"""Exact closure of the THM-2941 five-aligned/two-drift branch.

This verifier reconstructs the literal six-body boundary from the defining
danger combs.  It starts with the inclusive first-drift universe forced by
the THM-1094 component discrepancy bound, then splits the 4,702 projected-
suffix rows into two exact banks.

High-excess bank
    ``delta_E(z1) >= (88/1365) h_E``.  Delete ``D_z1`` and apply the
    six-tail first-apex gate to the residual.  This closes 3,827 of 4,084
    rows and leaves 42,912 integral second-drift candidates.

Subcritical bank
    ``delta_E(z1) < (88/1365) h_E``.  The remaining scalar gap and the
    component discrepancy bound give a nonempty finite ``z2`` interval on
    2,290 first rows, containing 7,218,110 row-labelled candidates.  Exact
    singleton integration together with the projected-suffix predicate leaves
    194,073 admissible pairs on 590 of the 618 suffix rows.

For each remaining pair put

    P = phi_L(C_E minus (D_z1 union D_z2)),  phi_L(t)=Lt mod 1.

Five aligned combs could complete the cover only if the compact set ``P``
were contained in an open normalized union of measure at most
``887/1365``.  On a body-safe ``1/L`` cell ``j`` let

    E_z(j) = {u in [0,1] : ||z(j+u)/L|| < 1/14}.

The exact De Morgan identity

    T minus P = intersection_j (E_z1(j) union E_z2(j))

turns this into rational interval arithmetic.  A prefix of the intersection
gives a rigorous lower bound for ``mu(P)``.  Every pair reaches
``mu(P) >= 887/1365`` and is therefore impossible; equality would already
contradict compact containment in a proper open set.

No finite search horizon is used as an exhaustive cutoff.  ``HORIZON`` only
reconstructs the already-rigorous suffix predicate: omitted labels obey
``delta_E(w) <= 6 r_E/(49 w)``.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import inspect
import math
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


def repo_root(start: Path) -> Path:
    for parent in (start, *start.parents):
        if (parent / "04-computation").is_dir() and (parent / "01-canon").is_dir():
            return parent
    raise RuntimeError("cannot locate repository root")


HERE = Path(__file__).resolve().parent
ROOT = repo_root(HERE)
SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_critical_scalar_wall_independent_thm2941.py"
)
EXPECTED_SOURCE_SHA256 = (
    "5d25a955fe184d6c1a3d8b632b4bbf901dc996ee46ad67c5748836fcc7134404"
)
EXPECTED_SUPPORT_SHA256 = (
    "5482e10635ecf72840bc0c083360fd7ddad65c2885d743820061bcba58cd5609"
)
DEFAULT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.out"
)

HORIZON = 7_000
ETA = F(88, 1_365)
U5 = F(478, 1_365)
ALIGNED_UNION_CAP = 1 - U5
PROJECTED_WALL_RATIO = F(2_275, 18_627)
BASE_LABEL = 15


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def normalized_sha256(path: Path) -> str:
    payload = path.read_bytes()
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: lone carriage return",
    )
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


spec = importlib.util.spec_from_file_location("thm2941_projected_source", SOURCE)
require(spec is not None and spec.loader is not None, "cannot load source engine")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)

SUPPORT_NAMES = (
    "merge_intervals",
    "danger_intervals",
    "carrier_for",
    "danger_primitive",
    "singleton_coverage",
)
support_payload = "\n".join(
    inspect.getsource(getattr(A, name)) for name in SUPPORT_NAMES
).encode()
require(
    normalized_sha256(SOURCE) == EXPECTED_SOURCE_SHA256,
    "pinned THM-2941 source changed",
)
require(
    hashlib.sha256(support_payload).hexdigest() == EXPECTED_SUPPORT_SHA256,
    "pinned carrier/singleton support changed",
)
require(
    A.BASE_LABEL == BASE_LABEL and A.RULER == 5_045_040,
    "source constants changed",
)
require(ETA == F(88, 1_365), "safe-surplus constant changed")
require(ALIGNED_UNION_CAP == F(887, 1_365), "aligned union cap changed")
require(
    PROJECTED_WALL_RATIO == F(2_275, 18_627),
    "projected wall changed",
)


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def merge_fraction(
    rows: list[tuple[F, F]],
) -> tuple[tuple[F, F], ...]:
    rows = sorted((left, right) for left, right in rows if left < right)
    if not rows:
        return ()
    out: list[list[F]] = [[rows[0][0], rows[0][1]]]
    for left, right in rows[1:]:
        if left <= out[-1][1]:
            out[-1][1] = max(out[-1][1], right)
        else:
            out.append([left, right])
    return tuple((left, right) for left, right in out)


def intersect_fraction(
    first: tuple[tuple[F, F], ...],
    second: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    out: list[tuple[F, F]] = []
    i = 0
    j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            out.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def interval_mass(rows: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in rows), F(0))


def danger_fraction(speed: int) -> tuple[tuple[F, F], ...]:
    """The open danger comb modulo one; endpoints are measure-null."""

    radius = F(1, 14 * speed)
    rows: list[tuple[F, F]] = []
    for index in range(speed):
        center = F(index, speed)
        left = center - radius
        right = center + radius
        if left < 0:
            rows.extend(((F(0), right), (left + 1, F(1))))
        elif right > 1:
            rows.extend(((left, F(1)), (F(0), right - 1)))
        else:
            rows.append((left, right))
    return merge_fraction(rows)


def subtract_fraction(
    carrier: tuple[tuple[F, F], ...],
    removed: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    out: list[tuple[F, F]] = []
    first_removed = 0
    for left, right in carrier:
        cursor = left
        while first_removed < len(removed) and removed[first_removed][1] <= left:
            first_removed += 1
        index = first_removed
        while index < len(removed) and removed[index][0] < right:
            bad_left, bad_right = removed[index]
            if cursor < bad_left:
                out.append((cursor, min(right, bad_left)))
            cursor = max(cursor, bad_right)
            if cursor >= right:
                break
            index += 1
        if cursor < right:
            out.append((cursor, right))
    return tuple(out)


def excess_numerator(
    carrier: tuple[tuple[int, int], ...],
    carrier_numerator: int,
    speed: int,
) -> int:
    """Return N with delta_E(speed)=N/(7*RULER*speed)."""

    covered = sum(
        A.danger_primitive(speed * right) - A.danger_primitive(speed * left)
        for left, right in carrier
    )
    return 7 * covered - carrier_numerator * speed


def excess_meets(numerator: int, speed: int, target: F) -> bool:
    return (
        numerator * target.denominator
        >= target.numerator * 7 * A.RULER * speed
    )


def body_cells(
    carrier: tuple[tuple[int, int], ...],
    canonical_l: int,
) -> tuple[int, ...]:
    require(A.RULER % canonical_l == 0, "body grid does not divide ruler")
    scale = A.RULER // canonical_l
    cells: list[int] = []
    for left, right in carrier:
        require(
            left % scale == 0 and right % scale == 0,
            "carrier endpoint left the body grid",
        )
        cells.extend(range(left // scale, right // scale))
    require(cells, "body carrier has no safe cells")
    return tuple(cells)


def phase_danger(
    cell: int,
    speed: int,
    canonical_l: int,
) -> tuple[tuple[F, F], ...]:
    """Danger subset in normalized coordinate u on body cell ``cell``."""

    start = F(speed * cell, canonical_l)
    finish = F(speed * (cell + 1), canonical_l)
    rows: list[tuple[F, F]] = []
    for tooth in range(floor_fraction(start) - 1, ceil_fraction(finish) + 2):
        left = max(
            F(0),
            F(canonical_l, speed) * (F(tooth) - F(1, 14)) - cell,
        )
        right = min(
            F(1),
            F(canonical_l, speed) * (F(tooth) + F(1, 14)) - cell,
        )
        if left < right:
            rows.append((left, right))
    return merge_fraction(rows)


def projected_safe_lower_bound(
    cells: tuple[int, ...],
    canonical_l: int,
    first: int,
    second: int,
) -> tuple[F, int]:
    """Return a prefix lower bound for projected-safe mass and prefix size."""

    common_danger: tuple[tuple[F, F], ...] = ((F(0), F(1)),)
    cells_used = 0
    for cell in cells:
        local_union = merge_fraction(
            [
                *phase_danger(cell, first, canonical_l),
                *phase_danger(cell, second, canonical_l),
            ]
        )
        common_danger = intersect_fraction(common_danger, local_union)
        cells_used += 1
        if interval_mass(common_danger) <= U5:
            break
    return 1 - interval_mass(common_danger), cells_used


def update_minimum(
    current: tuple[F, tuple[object, ...]] | None,
    candidate: tuple[F, tuple[object, ...]],
) -> tuple[F, tuple[object, ...]]:
    if current is None or candidate < current:
        return candidate
    return current


def profile(body: tuple[int, ...]) -> dict[str, object]:
    carrier = A.carrier_for(body)
    carrier_numerator = sum(right - left for left, right in carrier)
    h = F(carrier_numerator, A.RULER)
    components = len(carrier)
    canonical_l = 14 * math.lcm(*body)
    cells = body_cells(carrier, canonical_l)
    carrier_f = tuple(
        (F(left, A.RULER), F(right, A.RULER)) for left, right in carrier
    )
    first_bound = F(12 * components, 49) / (ETA * h)
    first_cap = floor_fraction(first_bound)
    projected_wall = PROJECTED_WALL_RATIO * canonical_l
    second_wall = max(BASE_LABEL, floor_fraction(projected_wall) + 1)

    candidate_first_rows = 0
    equality_cap_rows: list[tuple[tuple[int, ...], int]] = []

    high_rows = 0
    high_roots = False
    high_max_first = 0
    high_first_apex: list[tuple[object, ...]] = []
    high_typed_counts = {
        (False, False): 0,
        (False, True): 0,
        (True, False): 0,
        (True, True): 0,
    }
    high_pairs: list[tuple[object, ...]] = []
    high_killed = 0
    high_survivors: list[tuple[object, ...]] = []
    high_maximum_prefix = 0
    high_minimum: tuple[F, tuple[object, ...]] | None = None

    sub_analytic_rows = 0
    sub_pair_candidates = 0
    sub_suffix_rows: list[tuple[object, ...]] = []
    sub_pairs: list[tuple[object, ...]] = []
    sub_killed = 0
    sub_survivors: list[tuple[object, ...]] = []
    sub_maximum_prefix = 0
    sub_minimum: tuple[F, tuple[object, ...]] | None = None

    if first_bound.denominator == 1 and first_cap % canonical_l:
        equality_cap_rows.append((body, first_cap))

    lower = ETA * h
    for first in range(BASE_LABEL, first_cap + 1):
        if first % canonical_l == 0:
            continue
        candidate_first_rows += 1
        first_numerator = excess_numerator(carrier, carrier_numerator, first)
        first_delta = F(first_numerator, 7 * A.RULER * first)

        if first_delta >= lower:
            high_rows += 1
            high_roots = True
            high_max_first = max(high_max_first, first)
            first_coverage = h / 7 + first_delta
            residual = subtract_fraction(carrier_f, danger_fraction(first))
            residual_h = interval_mass(residual)
            residual_r = len(residual)
            require(
                residual_h == h - first_coverage,
                f"{body},z1={first}: residual mass mismatch",
            )
            require(
                residual_h > 0 and residual_r > 0,
                f"{body},z1={first}: residual disappeared",
            )
            apex_cap = floor_fraction(F(36 * residual_r, 7) / residual_h)
            second_floor = max(first + 1, second_wall)
            aligned_floor = canonical_l
            drift_possible = second_floor <= apex_cap
            aligned_possible = aligned_floor <= apex_cap
            if drift_possible or aligned_possible:
                first_row = (
                    body,
                    first,
                    h,
                    components,
                    canonical_l,
                    first_cap,
                    first_delta,
                    residual_h,
                    residual_r,
                    apex_cap,
                    second_floor,
                    aligned_floor,
                    drift_possible,
                    aligned_possible,
                )
                high_first_apex.append(first_row)
                high_typed_counts[(drift_possible, aligned_possible)] += 1
                for second in range(second_floor, apex_cap + 1):
                    if second % canonical_l == 0:
                        continue
                    pair = (
                        body,
                        first,
                        second,
                        canonical_l,
                        h,
                        components,
                        residual_h,
                        residual_r,
                        apex_cap,
                        second_floor,
                    )
                    high_pairs.append(pair)
                    projected_lower, cells_used = projected_safe_lower_bound(
                        cells,
                        canonical_l,
                        first,
                        second,
                    )
                    high_maximum_prefix = max(high_maximum_prefix, cells_used)
                    certificate = (*pair, projected_lower, cells_used)
                    if projected_lower >= ALIGNED_UNION_CAP:
                        high_killed += 1
                        high_minimum = update_minimum(
                            high_minimum,
                            (
                                projected_lower - ALIGNED_UNION_CAP,
                                certificate,
                            ),
                        )
                    else:
                        high_survivors.append(certificate)
            continue

        gap = lower - first_delta
        second_floor = max(first + 1, second_wall)
        second_cap = floor_fraction(F(6 * components, 49) / gap)
        if second_floor > second_cap:
            continue
        sub_analytic_rows += 1
        exact_suffix_hit = False
        row_pairs: list[tuple[object, ...]] = []
        for second in range(second_floor, second_cap + 1):
            if second % canonical_l == 0:
                continue
            sub_pair_candidates += 1
            second_numerator = excess_numerator(
                carrier,
                carrier_numerator,
                second,
            )
            if not excess_meets(second_numerator, second, gap):
                continue
            if second <= HORIZON:
                exact_suffix_hit = True
            second_delta = F(second_numerator, 7 * A.RULER * second)
            pair = (
                body,
                first,
                second,
                canonical_l,
                h,
                components,
                first_delta,
                gap,
                second_floor,
                second_cap,
                second_delta,
            )
            row_pairs.append(pair)
            projected_lower, cells_used = projected_safe_lower_bound(
                cells,
                canonical_l,
                first,
                second,
            )
            sub_maximum_prefix = max(sub_maximum_prefix, cells_used)
            certificate = (*pair, projected_lower, cells_used)
            if projected_lower >= ALIGNED_UNION_CAP:
                sub_killed += 1
                sub_minimum = update_minimum(
                    sub_minimum,
                    (
                        projected_lower - ALIGNED_UNION_CAP,
                        certificate,
                    ),
                )
            else:
                sub_survivors.append(certificate)

        high_tail_start = max(HORIZON + 1, second_floor)
        high_tail = F(6 * components, 49 * high_tail_start)
        tail_suffix_hit = high_tail >= gap
        if exact_suffix_hit or tail_suffix_hit:
            sub_suffix_rows.append(
                (
                    body,
                    first,
                    canonical_l,
                    h,
                    components,
                    first_delta,
                    gap,
                    second_floor,
                    second_cap,
                    exact_suffix_hit,
                    tail_suffix_hit,
                    high_tail_start,
                    high_tail,
                    len(row_pairs),
                )
            )
            sub_pairs.extend(row_pairs)
        else:
            require(
                not row_pairs,
                f"{body},z1={first}: exact pair escaped suffix predicate",
            )

    return {
        "body": body,
        "candidate_first_rows": candidate_first_rows,
        "equality_cap_rows": equality_cap_rows,
        "high_rows": high_rows,
        "high_roots": high_roots,
        "high_max_first": high_max_first,
        "high_first_apex": high_first_apex,
        "high_typed_counts": high_typed_counts,
        "high_pairs": high_pairs,
        "high_killed": high_killed,
        "high_survivors": high_survivors,
        "high_maximum_prefix": high_maximum_prefix,
        "high_minimum": high_minimum,
        "sub_analytic_rows": sub_analytic_rows,
        "sub_pair_candidates": sub_pair_candidates,
        "sub_suffix_rows": sub_suffix_rows,
        "sub_pairs": sub_pairs,
        "sub_killed": sub_killed,
        "sub_survivors": sub_survivors,
        "sub_maximum_prefix": sub_maximum_prefix,
        "sub_minimum": sub_minimum,
    }


def flatten(
    profiles: list[dict[str, object]],
    key: str,
) -> list[tuple[object, ...]]:
    return sorted(row for profile_row in profiles for row in profile_row[key])


def render(profiles: list[dict[str, object]]) -> str:
    require(len(profiles) == math.comb(14, 6) == 3_003, "root universe changed")
    profiles.sort(key=lambda row: row["body"])

    candidate_first_rows = sum(
        row["candidate_first_rows"] for row in profiles
    )
    equality_cap_rows = flatten(profiles, "equality_cap_rows")
    high_rows = sum(row["high_rows"] for row in profiles)
    high_first_apex = flatten(profiles, "high_first_apex")
    high_pairs = flatten(profiles, "high_pairs")
    high_killed = sum(row["high_killed"] for row in profiles)
    high_survivors = flatten(profiles, "high_survivors")
    high_maximum_prefix = max(row["high_maximum_prefix"] for row in profiles)
    high_minimum = min(
        (
            row["high_minimum"]
            for row in profiles
            if row["high_minimum"] is not None
        ),
        default=None,
    )
    high_typed_counts = {
        key: sum(row["high_typed_counts"][key] for row in profiles)
        for key in ((False, False), (False, True), (True, False), (True, True))
    }

    sub_analytic_rows = sum(row["sub_analytic_rows"] for row in profiles)
    sub_pair_candidates = sum(row["sub_pair_candidates"] for row in profiles)
    sub_suffix_rows = flatten(profiles, "sub_suffix_rows")
    sub_pairs = flatten(profiles, "sub_pairs")
    sub_killed = sum(row["sub_killed"] for row in profiles)
    sub_survivors = flatten(profiles, "sub_survivors")
    sub_maximum_prefix = max(row["sub_maximum_prefix"] for row in profiles)
    sub_minimum = min(
        (
            row["sub_minimum"]
            for row in profiles
            if row["sub_minimum"] is not None
        ),
        default=None,
    )

    first_apex_digest = hashlib.sha256(
        b"LRC14/five-aligned/two-drift/residual-first-apex/v1\n"
        + repr(tuple(high_first_apex)).encode()
    ).hexdigest()
    high_bank_digest = hashlib.sha256(
        b"LRC14/k5/high-excess-second-drift-bank/v1\n"
        + repr(tuple(high_pairs)).encode()
    ).hexdigest()
    high_survivor_digest = hashlib.sha256(
        b"LRC14/k5/projected-safe-set-pair-filter/v1\n"
        + repr(tuple(high_survivors)).encode()
    ).hexdigest()
    sub_bank_digest = hashlib.sha256(
        b"LRC14/k5/subcritical-exact-pair-bank/v1\n"
        + repr(tuple(sub_suffix_rows)).encode()
        + b"\n"
        + repr(tuple(sub_pairs)).encode()
    ).hexdigest()
    sub_survivor_digest = hashlib.sha256(
        b"LRC14/k5/projected-safe-subcritical-pairs/v1\n"
        + repr(tuple(sub_survivors)).encode()
    ).hexdigest()
    combined_survivor_digest = hashlib.sha256(
        b"LRC14/k5/five-aligned-two-drift-projected-closure/v1\n"
        + repr((tuple(high_survivors), tuple(sub_survivors))).encode()
    ).hexdigest()

    expected_equality_rows = [
        ((1, 2, 4, 6, 12, 14), 182),
        ((1, 3, 4, 9, 11, 12), 200),
        ((1, 4, 5, 6, 12, 14), 210),
        ((2, 3, 4, 6, 8, 14), 182),
        ((4, 5, 6, 7, 12, 14), 189),
    ]
    require(candidate_first_rows == 626_787, "inclusive first universe changed")
    require(equality_cap_rows == expected_equality_rows, "equality caps changed")
    require(high_rows == 4_084, "high-excess row count changed")
    require(
        sum(row["high_roots"] for row in profiles) == 2_309,
        "high-excess root count changed",
    )
    require(
        max(row["high_max_first"] for row in profiles) == 66,
        "high-excess maximum z1 changed",
    )
    require(len(high_first_apex) == 257, "first-apex row count changed")
    require(
        len({row[0] for row in high_first_apex}) == 183,
        "first-apex root count changed",
    )
    require(
        max(row[1] for row in high_first_apex) == 44,
        "first-apex maximum z1 changed",
    )
    require(
        high_typed_counts
        == {
            (False, False): 0,
            (False, True): 0,
            (True, False): 257,
            (True, True): 0,
        },
        "typed first-apex split changed",
    )
    require(
        first_apex_digest
        == "ac39be269cbd8695fb1cdefdd5d2bb8ce5512f636fd93626795c1aba949c3362",
        "first-apex digest changed",
    )
    require(len(high_pairs) == 42_912, "high pair bank changed")
    require(
        len({row[0] for row in high_pairs}) == 183
        and len({(row[0], row[1]) for row in high_pairs}) == 257,
        "high pair bank support changed",
    )
    require(
        high_bank_digest
        == "b488499d27625ccab7bf32d9fedfb50b091bf08db6dfdfc4b92ec98d66bef7ae",
        "high pair bank digest changed",
    )
    require(high_killed == len(high_pairs), "high bank no longer closes")
    require(not high_survivors, "high projected survivor bank is nonempty")
    require(high_maximum_prefix == 22, "high maximum prefix changed")
    require(
        high_minimum
        == (
            F(1, 378_105),
            (
                (1, 2, 3, 5, 9, 10),
                24,
                277,
                1_260,
                F(136, 315),
                14,
                F(43, 126),
                26,
                391,
                154,
                F(180, 277),
                1,
            ),
        ),
        "high minimum certificate changed",
    )
    require(
        high_survivor_digest
        == "3f001616fc494fca31891589993823d9aab957b9a5b73b96159cc899361d27bd",
        "high survivor digest changed",
    )

    require(sub_analytic_rows == 2_290, "subcritical analytic rows changed")
    require(len(sub_suffix_rows) == 618, "subcritical suffix rows changed")
    require(
        len({row[0] for row in sub_suffix_rows}) == 435,
        "subcritical suffix roots changed",
    )
    require(
        high_rows + len(sub_suffix_rows) == 4_702,
        "projected suffix split changed",
    )
    require(
        max(row[1] for row in sub_suffix_rows) == 55,
        "subcritical maximum z1 changed",
    )
    require(
        max(row[8] for row in sub_suffix_rows) == 688_956,
        "subcritical maximum z2 cap changed",
    )
    require(
        sub_pair_candidates == 7_218_110,
        "subcritical finite z2 universe changed",
    )
    require(len(sub_pairs) == 194_073, "subcritical exact pair bank changed")
    require(
        len({row[0] for row in sub_pairs}) == 407,
        "subcritical exact pair roots changed",
    )
    require(
        len({(row[0], row[1]) for row in sub_pairs}) == 590,
        "subcritical exact (E,z1) rows changed",
    )
    require(
        max(row[2] for row in sub_pairs) == 327_502,
        "subcritical maximum admissible z2 changed",
    )
    require(
        sub_bank_digest
        == "4827dabbc24f07bdf9e628e6b706fe66e7ed451c893cece993499a4cf13394c0",
        "subcritical bank digest changed",
    )
    require(sub_killed == len(sub_pairs), "subcritical bank no longer closes")
    require(not sub_survivors, "subcritical projected survivor bank is nonempty")
    require(sub_maximum_prefix == 871, "subcritical maximum prefix changed")
    require(
        sub_minimum
        == (
            F(1, 378_105),
            (
                (2, 3, 4, 5, 6, 9),
                22,
                554,
                2_520,
                F(13, 30),
                16,
                F(191, 6_930),
                F(13, 34_650),
                308,
                5_221,
                F(4, 9_695),
                F(180, 277),
                1,
            ),
        ),
        "subcritical minimum certificate changed",
    )
    require(
        sub_survivor_digest
        == "340938a1ca495881500b8a704aeb88c2e4fa5b976eb031fb675f23c3b028c33a",
        "subcritical survivor digest changed",
    )

    suffix_mode_counts = {
        (exact_hit, tail_hit): sum(
            row[9] == exact_hit and row[10] == tail_hit
            for row in sub_suffix_rows
        )
        for exact_hit in (False, True)
        for tail_hit in (False, True)
    }
    pair_root_counts = (
        len({row[0] for row in high_pairs}),
        len({(row[0], row[1]) for row in high_pairs}),
    )
    require(
        suffix_mode_counts
        == {
            (False, False): 0,
            (False, True): 56,
            (True, False): 452,
            (True, True): 110,
        },
        "subcritical suffix modes changed",
    )
    require(
        len(high_pairs) + len(sub_pairs) == 236_985
        and high_killed + sub_killed == 236_985,
        "combined pair closure count changed",
    )
    require(
        combined_survivor_digest
        == "e728c27d7cec406eb3c61923d47e74d11232e2e3f4b88a268ec20d512f9a53ae",
        "combined survivor digest changed",
    )

    lines = [
        "LRC14 THM-2941 five-aligned/two-drift projected closure",
        f"source_sha256={normalized_sha256(SOURCE)}",
        f"support_sha256={hashlib.sha256(support_payload).hexdigest()}",
        "universe=(six_body_roots=3003,body_labels=1..14,"
        "z1,z2>=15,strictly_increasing,nonmultiples_of_L)",
        f"eta5={ftext(ETA)};u5={ftext(U5)};"
        f"aligned_union_cap={ftext(ALIGNED_UNION_CAP)};"
        "projected_wall_ratio=2275/18627(=325/2661)",
        f"inclusive_first_rows={candidate_first_rows};"
        f"inclusive_boundary_rows={len(equality_cap_rows)};"
        f"projected_suffix_rows={high_rows + len(sub_suffix_rows)}",
        f"inclusive_boundary_bank={tuple(equality_cap_rows)}",
        "branch=high_excess;"
        f"first_rows={high_rows};first_apex_closed={high_rows-len(high_first_apex)};"
        f"first_apex_rows={len(high_first_apex)};"
        f"first_apex_roots={len({row[0] for row in high_first_apex})};"
        f"typed_counts={high_typed_counts}",
        f"branch=high_excess;candidate_pairs={len(high_pairs)};"
        f"pair_roots={pair_root_counts[0]};"
        f"pair_(E,z1)_rows={pair_root_counts[1]};"
        f"killed={high_killed};survivors={len(high_survivors)};"
        f"maximum_cell_prefix={high_maximum_prefix}",
        f"branch=high_excess;minimum_certificate={high_minimum}",
        f"branch=high_excess;first_apex_digest={first_apex_digest}",
        f"branch=high_excess;pair_bank_digest={high_bank_digest}",
        f"branch=high_excess;survivor_digest={high_survivor_digest}",
        "branch=subcritical;"
        f"analytic_first_rows={sub_analytic_rows};"
        f"projected_suffix_rows={len(sub_suffix_rows)};"
        f"suffix_roots={len({row[0] for row in sub_suffix_rows})};"
        f"suffix_mode_counts={suffix_mode_counts}",
        f"branch=subcritical;finite_z2_candidates={sub_pair_candidates};"
        f"exact_admissible_pairs={len(sub_pairs)};"
        f"pair_roots={len({row[0] for row in sub_pairs})};"
        f"pair_(E,z1)_rows={len({(row[0],row[1]) for row in sub_pairs})};"
        f"killed={sub_killed};survivors={len(sub_survivors)};"
        f"maximum_cell_prefix={sub_maximum_prefix}",
        f"branch=subcritical;minimum_certificate={sub_minimum}",
        f"branch=subcritical;pair_bank_digest={sub_bank_digest}",
        f"branch=subcritical;survivor_digest={sub_survivor_digest}",
        f"combined_candidate_pairs={len(high_pairs)+len(sub_pairs)};"
        f"combined_killed={high_killed+sub_killed};"
        f"combined_survivors={len(high_survivors)+len(sub_survivors)}",
        f"combined_survivor_digest={combined_survivor_digest}",
        "all_exact_controls=PASS",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    roots = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        profiles = [profile(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(profile, roots, chunksize=1))
    output = render(profiles)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
