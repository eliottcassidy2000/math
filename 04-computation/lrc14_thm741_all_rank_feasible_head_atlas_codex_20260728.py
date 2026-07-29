#!/usr/bin/env python3
r"""Exact ranked-head closure for every rank-feasible THM-741 root.

For a nine-speed root ``E`` in ``{1,...,14}``, write

    c_E(w)=|G(E) intersect D_w|,        m_E=|G(E)|.

Let ``q_1>=q_2>=q_3`` be the three largest coverages over all ``w>=15`` and
put ``R_E=m_E-q_1-q_2-q_3``.  THM-735(ii) gives

    c_E(w) < m_E/7 + (99/70) r_E/(7w).

When ``R_E>m_E/7``, this makes the ranked head

    H_E={w>=15:c_E(w)>=R_E}

finite.  If ``K=|H_E|``, every four-speed set not contained in ``H_E`` is
strictly closed by the union bound: its fourth coverage is below ``R_E``.
It remains only to inspect the finite head quadruples whose individual
coverage sum is at least ``m_E``.

This companion starts from all ``C(14,9)=2002`` roots, reconstructs every
finite global head, selects exactly the 816 nonpositive roots with
``R_E>m_E/7``, and finds 63,265 union-bound-dangerous head quadruples.  It
checks all of them by cached exact twelve-speed carriers and sparse
fourth-comb subtraction.  As a positive control it retains the older
exhaustive check of every head quadruple when ``K<=10``, for 79,930 literal
carrier checks in total.  Every survivor is positive.  The computation
itself treats added speeds at least 15; THM-741 composes it with THM-738 for
the small-speed chamber.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
from math import lcm
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
FIRST_EXTERNAL = 15
BASE_HORIZON = 600
BODIES = tuple(combinations(range(1, 15), 9))
FLOOD_CORE = tuple(range(8, 15))

EXPECTED_TOP4_POSITIVE = 584
EXPECTED_TOP4_NONPOSITIVE = 1418
EXPECTED_RANK_FEASIBLE = 816
EXPECTED_RANK_IMPOSSIBLE = 602
EXPECTED_SELECTED = 816
EXPECTED_K_HISTOGRAM = {
    4: 87,
    5: 74,
    6: 62,
    7: 67,
    8: 43,
    9: 37,
    10: 36,
    11: 32,
    12: 20,
    13: 16,
    14: 18,
    15: 19,
    16: 12,
    17: 8,
    18: 10,
    19: 8,
    20: 10,
    21: 4,
    22: 10,
    23: 12,
    24: 6,
    25: 6,
    26: 7,
    27: 6,
    28: 8,
    29: 7,
    30: 6,
    31: 7,
    32: 5,
    33: 2,
    34: 6,
    35: 3,
    36: 5,
    37: 2,
    38: 1,
    39: 8,
    40: 4,
    41: 3,
    42: 1,
    43: 2,
    44: 5,
    45: 3,
    46: 2,
    47: 1,
    48: 1,
    49: 2,
    50: 7,
    51: 3,
    52: 1,
    53: 1,
    54: 3,
    55: 3,
    56: 2,
    57: 1,
    58: 1,
    59: 1,
    61: 2,
    62: 1,
    63: 2,
    64: 4,
    65: 1,
    66: 1,
    67: 3,
    68: 3,
    69: 2,
    70: 4,
    71: 1,
    72: 1,
    74: 3,
    76: 1,
    77: 1,
    78: 1,
    79: 1,
    80: 2,
    81: 1,
    83: 1,
    84: 1,
    85: 1,
    86: 2,
    91: 2,
    94: 2,
    96: 1,
    98: 1,
    100: 1,
    105: 1,
    106: 1,
    111: 1,
    116: 1,
    120: 1,
    127: 1,
    130: 1,
    131: 1,
    133: 1,
    135: 1,
    139: 1,
    140: 1,
    148: 1,
    152: 1,
    158: 1,
    175: 1,
    176: 1,
    178: 1,
    188: 1,
    197: 1,
    206: 1,
    208: 1,
    216: 1,
    228: 1,
    241: 1,
    250: 1,
    260: 1,
    275: 1,
    280: 1,
    288: 1,
    293: 1,
    299: 1,
    315: 1,
    318: 1,
    328: 1,
    365: 1,
    369: 1,
    370: 1,
    382: 1,
    392: 1,
    407: 1,
    502: 1,
    515: 1,
    554: 1,
    559: 1,
    574: 1,
    742: 1,
    1320: 1,
    1321: 1,
    1763: 1,
    5975: 1,
    10519: 1,
}
EXPECTED_DANGEROUS = 63265
EXPECTED_LITERAL_CHECKS = 79930
EXPECTED_MINIMUM = (
    F(362, 30107),
    (1, 4, 5, 6, 7, 10, 11, 12, 13, 16, 17, 18, 23),
)
EXPECTED_MINIMUM_COMPONENTS = (
    6,
    F(577, 22610),
    (2, 3, 5, 7, 8, 9, 11, 12, 13, 15, 17, 19, 20),
)
EXPECTED_FLOOD_OVERLAP = ((2, 3), (2, 5), (2, 6), (2, 7), (3, 7))
EXPECTED_MAXIMUM_THRESHOLD = (
    (2, 3, 4, 6, 7, 8, 9, 10, 14),
    F(42810768, 145),
    295247,
)
EXPECTED_MAXIMUM_HEAD_SPEED = 90850
EXPECTED_FEASIBLE_BODY_DIGEST = (
    "4db4ae0257114fc91e1a4cd5ad05bdd563f576a7d5dacb15f2d3ba98c9abd47e"
)
EXPECTED_IMPOSSIBLE_BODY_DIGEST = (
    "a4ad9f9b5a8ce16103450ac05684cd13ec33637e8e3737218a9186086cd639d4"
)
EXPECTED_PROVED_UNION_DIGEST = (
    "e79dc2d7e2a77adfe0901c9fdcc198c6d3f2338bda9d36ba0d5f02e2b8c46133"
)
EXPECTED_RESIDUAL_BODY_DIGEST = (
    "cd7d255516269e2ccc57a43dfc4d3e5eed2f7d364e9fce82d09cd7d651b7a0b9"
)
EXPECTED_CLASSIFICATION_DIGEST = (
    "d89f6ec420e8b2418b97658a3b4499fa4eb0c489be674c8dd6744061805aad71"
)
EXPECTED_K10_CLOSED_BODY_DIGEST = (
    "321d427048021d633af245ddc990a8a4d594fad9746fbe8f0ec50c72ab0f6101"
)
EXPECTED_K10_HEAD_DIGEST = (
    "b4bfd1858e32edc188b0d161933c55a9a8b3bb27a5c5c110ec7e641a5a6761e6"
)
EXPECTED_K10_CANDIDATE_DIGEST = (
    "29343d70edfab1eaedd64a53fcf555715c0338f661d81a2f9729c68d4d00179a"
)
EXPECTED_K10_LITERAL_DIGEST = (
    "24baf0abcbedd117b7c729b645126f69a2ccd7a93adc7692b7619b63bbea4974"
)
EXPECTED_K11_20_CLOSED_BODY_DIGEST = (
    "9c1fa71e4f50797987d0afd2765ccccf9cf75428215f158e8cba8bbda376576d"
)
EXPECTED_K11_20_HEAD_DIGEST = (
    "80da08f837460543ea5182a31be61f3323e73ca93800d50e4d0fbe283700a919"
)
EXPECTED_K11_20_CANDIDATE_DIGEST = (
    "8d87bb818a31fa857511f4ab2c31d9987ea9a4914f439e31ebce41b93805077d"
)
EXPECTED_K11_20_LITERAL_DIGEST = (
    "50e6c6960ac07d785ad8562a0fa376e2c6aa3e70304a2cf328ef426e0ce574a5"
)
EXPECTED_K21_50_CLOSED_BODY_DIGEST = (
    "4c45f25dbbfffff8888ed0e020ad1cc7b8d79ec08e4895f5a46bb768342c0e06"
)
EXPECTED_K21_50_HEAD_DIGEST = (
    "b157da9a19a2297e1f21a3eacd6ffba55743d4a85f94a74e5a997184d3cd158d"
)
EXPECTED_K21_50_CANDIDATE_DIGEST = (
    "165af80d7828c12b63c3d990c82dc7ffaa52bf69799c8ced45146624a541c454"
)
EXPECTED_K21_50_LITERAL_DIGEST = (
    "369f8c8eb980aa8af9001d875932b4fe5608d2165f91214a059cce68af1d9a29"
)
EXPECTED_K51_100_CLOSED_BODY_DIGEST = (
    "1410cfe933af10390a1f5b3ad2a7eb5639f274c34c596f4732752d0b3d44e2a9"
)
EXPECTED_K51_100_HEAD_DIGEST = (
    "1a0430bfcbe0d28f890dad5bc07ebc1d06077657e7fa0b17f1c75a2f341bc471"
)
EXPECTED_K51_100_CANDIDATE_DIGEST = (
    "cf6e77325c44077fdcec70de12bb97076849b8f15f0da1134a63c48590e57eca"
)
EXPECTED_K51_100_LITERAL_DIGEST = (
    "db283faac44bc9aeb7227165d48ea4ef986c89870e5f064ff68a2d1cb9d68558"
)
EXPECTED_K101_200_CLOSED_BODY_DIGEST = (
    "a489d381868f0268f228260ea28488eae1029f1094d27daef0a772ab6603e9a7"
)
EXPECTED_K101_200_HEAD_DIGEST = (
    "bbc882664d6d7695d54010af6682eb2881db16e44a2eb9eec1b536f8ad9b3b20"
)
EXPECTED_K101_200_CANDIDATE_DIGEST = (
    "72a836d68a3985593cf8537d5be178141f50560a260fc8de0b30919b7aeb3f4e"
)
EXPECTED_K101_200_LITERAL_DIGEST = (
    "faf042c1dc790e57ef5f62dd7d759b68a656d218cf4e337af62e656ddb137c3f"
)
EXPECTED_K201_500_CLOSED_BODY_DIGEST = (
    "bbef092014d583a921872bfa3623c1d595eea3847e5da90a03b21777fa69a6f3"
)
EXPECTED_K201_500_HEAD_DIGEST = (
    "f3b23a394cfe7406a10c596caac3b9deb6ac269c55e6964dd80b8796c2b18bf5"
)
EXPECTED_K201_500_CANDIDATE_DIGEST = (
    "8401dcc6b2cb17155b515c280957229745e52d798a6f135a44907b1f855d5eb4"
)
EXPECTED_K201_500_LITERAL_DIGEST = (
    "899bc5cbf46e6625f50a9a1333bc3e8a32a51b510099978ad4a08a2a4f58cd35"
)
EXPECTED_K501_1000_CLOSED_BODY_DIGEST = (
    "06fd145190a33fea71309700b5479f959dbc3e68365e96465139f37466fc88ea"
)
EXPECTED_K501_1000_HEAD_DIGEST = (
    "8f90d8a4b2d4a5b3396bea18c605942e90426c55b36c50bc0783f9405085525d"
)
EXPECTED_K501_1000_CANDIDATE_DIGEST = (
    "8acb1f85fe03fd3b44e556fec3e6f0e70dd1d88ce46ac8a4099166e15135163a"
)
EXPECTED_K501_1000_LITERAL_DIGEST = (
    "25186b7742eb0584dba901157840f3e863225487ce20c8df38336cc9ca3199b6"
)
EXPECTED_K1001_10000_CLOSED_BODY_DIGEST = (
    "1e32a269f947ab5f72198a1d0f9d5b087cd4434c6776287dcd52f6a5968405ac"
)
EXPECTED_K1001_10000_HEAD_DIGEST = (
    "7a7d754f96fa8f47d66dfb70e684e773acd31b389aeaad22c814e16932945cff"
)
EXPECTED_K1001_10000_CANDIDATE_DIGEST = (
    "c82d6213251e982303735a34b37c09f6b11639426e2de80e3bb94bbdf0ce91fd"
)
EXPECTED_K1001_10000_LITERAL_DIGEST = (
    "fd6707b7cf179028bce7e5769d9aded2a61ebbb3759d3ad25779340642348504"
)
EXPECTED_K10001_PLUS_CLOSED_BODY_DIGEST = (
    "8a66f7362df0d9c63a26e257b75fe704b00a6626330c372171abbaa0c9a39029"
)
EXPECTED_K10001_PLUS_HEAD_DIGEST = (
    "5f96bb4708166bd1e2c95b197b9abf6ae856e42aba7dc7c568dcf3344b115c66"
)
EXPECTED_K10001_PLUS_CANDIDATE_DIGEST = (
    "b2eabf92144a059e8b7c61bf44bd4df3bd5db0db61d497be9bf29856e296b973"
)
EXPECTED_K10001_PLUS_LITERAL_DIGEST = (
    "6f52ad7e988b307e6d9fbf5faddea50ad66512b749b442313fb33953f1bf4b8f"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(file_sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location(
        "thm741_ranked_head_k10_dependency",
        CORE_PATH,
    )
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CORE = load_core()


def fraction_text(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def tooth_primitive(value: F) -> F:
    """Primitive of the unit-periodic radius-1/14 tooth indicator."""

    integer = value.numerator // value.denominator
    remainder = value - integer
    return (
        F(integer, 7)
        + min(remainder, F(1, 14))
        + max(F(0), remainder - F(13, 14))
    )


def coverage(good: list[tuple[F, F]], root_m: F, speed: int) -> F:
    del root_m
    return sum(
        (
            tooth_primitive(speed * right)
            - tooth_primitive(speed * left)
        )
        / speed
        for left, right in good
    )


def coverages_many(
    good: list[tuple[F, F]],
    speeds: range | list[int],
    batch_size: int = 20000,
) -> list[tuple[F, int]]:
    """Evaluate the periodic primitive exactly in bounded int64 batches."""

    speed_list = list(speeds)
    if not speed_list:
        return []
    denominator = lcm(
        *(endpoint.denominator for interval in good for endpoint in interval)
    )
    endpoint_ints = np.array(
        [
            endpoint.numerator * (denominator // endpoint.denominator)
            for interval in good
            for endpoint in interval
        ],
        dtype=np.int64,
    )
    require(
        4 * len(good) * denominator * max(speed_list) < 2**62,
        "vector primitive would overflow int64",
    )
    out: list[tuple[F, int]] = []
    for start in range(0, len(speed_list), batch_size):
        speed_batch = speed_list[start : start + batch_size]
        speed_array = np.array(speed_batch, dtype=np.int64)
        scaled = endpoint_ints[:, None] * speed_array[None, :]
        quotient, remainder = np.divmod(scaled, denominator)
        primitive = (
            2 * denominator * quotient
            + np.minimum(14 * remainder, denominator)
            + np.maximum(0, 14 * remainder - 13 * denominator)
        )
        numerators = (primitive[1::2] - primitive[0::2]).sum(axis=0)
        out.extend(
            (F(int(numerator), 14 * denominator * speed), speed)
            for numerator, speed in zip(numerators, speed_batch)
        )
    return out


def local_bad_pieces(left: F, right: F, speed: int) -> list[tuple[F, F]]:
    """Build only teeth of ``D_speed`` meeting one carrier component."""

    radius = F(1, 14 * speed)
    first = (speed * (left - radius)).numerator // (
        speed * (left - radius)
    ).denominator - 1
    upper = speed * (right + radius)
    last = -((-upper.numerator) // upper.denominator) + 1
    pieces: list[tuple[F, F]] = []
    for tooth in range(first, last + 1):
        center = F(tooth, speed)
        lo = max(left, center - radius)
        hi = min(right, center + radius)
        if lo < hi:
            pieces.append((lo, hi))
    return pieces


def subtract_local(
    carrier: list[tuple[F, F]],
    speeds: tuple[int, ...],
) -> list[tuple[F, F]]:
    """Remove several combs simultaneously using only local endpoint events."""

    out: list[tuple[F, F]] = []
    for left, right in carrier:
        pieces = sorted(
            piece
            for speed in speeds
            for piece in local_bad_pieces(left, right, speed)
        )
        cursor = left
        for lo, hi in pieces:
            if hi <= cursor:
                continue
            if lo > cursor:
                out.append((cursor, lo))
            cursor = max(cursor, hi)
            if cursor >= right:
                break
        if cursor < right:
            out.append((cursor, right))
    return out


def profile_body(body: tuple[int, ...]) -> dict[str, object]:
    """Build one exact global ranked-head profile."""

    good, root_r, root_m = CORE.good_norm(body)
    require(root_r == len(good) and root_m > 0, f"bad root carrier {body}")

    rows = coverages_many(good, range(FIRST_EXTERNAL, BASE_HORIZON + 1))
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    top3 = tuple(ranked[:3])
    fourth = ranked[3][0]
    require(fourth > root_m / 7, f"fourth rank misses limiting cap {body}")
    fourth_threshold = CORE.S2 * root_r / (7 * (fourth - root_m / 7))
    require(
        fourth_threshold < BASE_HORIZON + 1,
        f"W=600 does not seal the global top three {body}",
    )
    top4_margin = root_m - sum((value for value, _ in ranked[:4]), F(0))
    residual = root_m - sum((value for value, _ in top3), F(0))

    result: dict[str, object] = {
        "body": body,
        "r": root_r,
        "m": root_m,
        "top3": top3,
        "top4_margin": top4_margin,
        "residual": residual,
        "rank_feasible": residual > root_m / 7,
        "selected": False,
    }
    if top4_margin > 0 or residual <= root_m / 7:
        return result

    threshold = CORE.S2 * root_r / (7 * (residual - root_m / 7))
    tail_first = ceiling(threshold)
    head = [(value, speed) for value, speed in rows if value >= residual]

    if tail_first > BASE_HORIZON + 1:
        head.extend(
            (value, speed)
            for value, speed in coverages_many(
                good,
                range(BASE_HORIZON + 1, tail_first),
            )
            if value >= residual
        )

    # Every omitted speed is below residual: scanned omissions by exact
    # comparison, and the unscanned tail by the strict THM-735(ii) cap.
    head.sort(key=lambda item: (-item[0], item[1]))
    require(len(head) >= 4, f"nonpositive row has fewer than four head speeds {body}")
    require(
        all(value >= residual for value, _ in head),
        f"bad head membership {body}",
    )
    result.update(
        {
            "threshold": threshold,
            "tail_first": tail_first,
            "K": len(head),
            "head": tuple(head),
            "selected": True,
        }
    )
    return result


def close_selected(row: dict[str, object]) -> dict[str, object]:
    """Check every dangerous quadruple, plus the legacy small-head controls."""

    body = row["body"]
    root_m = row["m"]
    head = row["head"]
    require(isinstance(body, tuple), "body type changed")
    require(isinstance(root_m, F), "root measure type changed")
    require(isinstance(head, tuple), "head type changed")
    good, root_r, replay_m = CORE.good_norm(body)
    require(root_r == row["r"] and replay_m == root_m, f"root replay changed {body}")

    ranked_head = sorted(head, key=lambda item: (-item[0], item[1]))
    exhaustive_small_head = len(ranked_head) <= 10
    checked_records: list[
        tuple[tuple[int, int, int, int], tuple[int, int, int], int]
    ] = []
    for first in range(len(ranked_head) - 3):
        if not exhaustive_small_head and (
            ranked_head[first][0]
            + ranked_head[first + 1][0]
            + ranked_head[first + 2][0]
            + ranked_head[first + 3][0]
            < root_m
        ):
            break
        for second in range(first + 1, len(ranked_head) - 2):
            if not exhaustive_small_head and (
                ranked_head[first][0]
                + ranked_head[second][0]
                + ranked_head[second + 1][0]
                + ranked_head[second + 2][0]
                < root_m
            ):
                break
            for third in range(second + 1, len(ranked_head) - 1):
                prefix = (
                    ranked_head[first][0]
                    + ranked_head[second][0]
                    + ranked_head[third][0]
                )
                if (
                    not exhaustive_small_head
                    and prefix + ranked_head[third + 1][0] < root_m
                ):
                    break
                for fourth in range(third + 1, len(ranked_head)):
                    if (
                        not exhaustive_small_head
                        and prefix + ranked_head[fourth][0] < root_m
                    ):
                        break
                    rank_triple = (
                        ranked_head[first][1],
                        ranked_head[second][1],
                        ranked_head[third][1],
                    )
                    fourth_speed = ranked_head[fourth][1]
                    quadruple = tuple(sorted((*rank_triple, fourth_speed)))
                    checked_records.append(
                        (quadruple, tuple(sorted(rank_triple)), fourth_speed)
                    )
    checked_records.sort()
    require(
        len(checked_records)
        == len({record[0] for record in checked_records}),
        f"duplicate checked quadruple {body}",
    )
    coverage_by_speed = {speed: value for value, speed in ranked_head}
    dangerous = sum(
        sum((coverage_by_speed[speed] for speed in record[0]), F(0)) >= root_m
        for record in checked_records
    )
    checked = len(checked_records)
    minimum: tuple[F, tuple[int, ...], int, tuple[int, ...]] | None = None
    minimum_components: tuple[int, F, tuple[int, ...], tuple[int, ...]] | None = None
    candidate_rows: list[str] = []
    literal_rows: list[str] = []

    triple_carriers: dict[tuple[int, int, int], list[tuple[F, F]]] = {}
    checked_triples: set[tuple[int, int, int]] = set()
    for speeds, rank_triple, fourth_speed in checked_records:
        triple_carrier = triple_carriers.get(rank_triple)
        if triple_carrier is None:
            triple_carrier = subtract_local(good, rank_triple)
            triple_carriers[rank_triple] = triple_carrier
        nested = subtract_local(triple_carrier, (fourth_speed,))
        nested_m = sum((right - left for left, right in nested), F(0))

        family = tuple(sorted((*body, *speeds)))
        direct_r = len(nested)
        direct_m = nested_m
        require(len(family) == 13 and len(set(family)) == 13, f"bad family {family}")
        if rank_triple not in checked_triples:
            require(
                subtract_local(good, speeds) == nested,
                f"cached-triple/direct paths disagree {body},{speeds}",
            )
            checked_triples.add(rank_triple)
        require(direct_m > 0, f"nonpositive dangerous survivor {body},{speeds}")

        record = (direct_m, family, direct_r, speeds)
        if minimum is None or record < minimum:
            minimum = record
        component_record = (direct_r, direct_m, family, speeds)
        if minimum_components is None or component_record < minimum_components:
            minimum_components = component_record
        body_text = ",".join(map(str, body))
        speed_text = ",".join(map(str, speeds))
        candidate_rows.append(f"E={body_text};Q={speed_text}\n")
        literal_rows.append(
            f"E={body_text};Q={speed_text};"
            f"L={fraction_text(direct_m)};r={direct_r};tight=-\n"
        )

    return {
        "body": body,
        "K": row["K"],
        "dangerous": dangerous,
        "checked": checked,
        "minimum": minimum,
        "minimum_components": minimum_components,
        "candidate_rows": tuple(candidate_rows),
        "literal_rows": tuple(literal_rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(10, os.cpu_count() or 1))
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(len(BODIES) == 2002 and len(set(BODIES)) == 2002, "root census changed")
    require(CORE.S2 == F(99, 70) and CORE.S2**2 > 2, "sqrt(2) majorant changed")

    control_body = BODIES[0]
    control_good, _, control_m = CORE.good_norm(control_body)
    for speed in (15, 23, 601, 5003):
        _, _, core_carrier = CORE.subtract(control_good, speed)
        require(
            subtract_local(control_good, (speed,)) == core_carrier,
            f"local/core subtraction control changed at {speed}",
        )
        require(
            coverage(control_good, control_m, speed)
            == control_m - CORE.subtract_sparse(control_good, speed),
            f"periodic-primitive coverage control changed at {speed}",
        )
    vector_controls = coverages_many(control_good, [15, 23, 601, 5003])
    require(
        vector_controls
        == [
            (coverage(control_good, control_m, speed), speed)
            for speed in (15, 23, 601, 5003)
        ],
        "vector/scalar periodic-primitive controls changed",
    )
    control_speeds = (15, 16, 17, 18)
    control_direct, _, _ = CORE.good_norm((*control_body, *control_speeds))
    require(
        subtract_local(control_good, control_speeds) == control_direct,
        "simultaneous local/direct-union control changed",
    )

    context = mp.get_context("spawn")
    if args.workers == 1:
        profiles = list(map(profile_body, BODIES))
    else:
        with context.Pool(args.workers) as pool:
            profiles = pool.map(profile_body, BODIES, chunksize=2)
    require(
        tuple(row["body"] for row in profiles) == BODIES,
        "parallel profile order changed",
    )

    positive = [row for row in profiles if row["top4_margin"] > 0]
    nonpositive = [row for row in profiles if row["top4_margin"] <= 0]
    feasible = [row for row in nonpositive if row["rank_feasible"]]
    impossible = [row for row in nonpositive if not row["rank_feasible"]]
    selected = [row for row in feasible if row["selected"]]
    require(
        len(positive) == EXPECTED_TOP4_POSITIVE
        and len(nonpositive) == EXPECTED_TOP4_NONPOSITIVE,
        "top-four census changed",
    )
    require(
        len(feasible) == EXPECTED_RANK_FEASIBLE
        and len(impossible) == EXPECTED_RANK_IMPOSSIBLE,
        "rank-feasibility census changed",
    )
    require(len(selected) == EXPECTED_SELECTED, "rank-feasible root count changed")

    histogram = {}
    for size in sorted({int(row["K"]) for row in selected}):
        count = sum(row["K"] == size for row in selected)
        if count:
            histogram[size] = count
    require(histogram == EXPECTED_K_HISTOGRAM, "K histogram changed")

    if args.workers == 1:
        closures = list(map(close_selected, selected))
    else:
        with context.Pool(args.workers) as pool:
            closures = pool.map(close_selected, selected, chunksize=1)
    require(
        tuple(row["body"] for row in closures)
        == tuple(row["body"] for row in selected),
        "parallel closure order changed",
    )

    dangerous = sum(int(row["dangerous"]) for row in closures)
    literal_checks = sum(int(row["checked"]) for row in closures)
    dangerous_band_counts = {
        tag: sum(
            int(closure["dangerous"])
            for closure in closures
            if lower <= int(closure["K"]) <= upper
        )
        for lower, upper, tag in (
            (4, 10, "KLE10"),
            (11, 20, "K11-20"),
            (21, 50, "K21-50"),
            (51, 100, "K51-100"),
            (101, 200, "K101-200"),
            (201, 500, "K201-500"),
            (501, 1000, "K501-1000"),
            (1001, 10000, "K1001-10000"),
            (10001, 10**9, "K10001-plus"),
        )
    }
    require(
        dangerous == EXPECTED_DANGEROUS,
        f"dangerous-quadruple count changed: {dangerous}; {dangerous_band_counts}",
    )
    require(
        literal_checks == EXPECTED_LITERAL_CHECKS,
        f"literal-check count changed: {literal_checks}",
    )
    minima = [row["minimum"] for row in closures if row["minimum"] is not None]
    require(minima, "no dangerous quadruples were checked")
    global_minimum = min(minima)
    require(
        (global_minimum[0], global_minimum[1]) == EXPECTED_MINIMUM,
        "minimum dangerous survivor changed",
    )
    component_minima = [
        row["minimum_components"]
        for row in closures
        if row["minimum_components"] is not None
    ]
    global_minimum_components = min(component_minima)
    require(
        global_minimum_components[:3] == EXPECTED_MINIMUM_COMPONENTS,
        "minimum component-count survivor changed",
    )

    def band_digests(
        lower: int,
        upper: int,
        tag: str,
    ) -> tuple[str, str, str, str]:
        band_profiles = [
            row for row in selected if lower <= int(row["K"]) <= upper
        ]
        band_bodies = {row["body"] for row in band_profiles}
        band_closures = [
            row for row in closures if row["body"] in band_bodies
        ]
        closed_text = "\n".join(
            ",".join(map(str, row["body"])) for row in band_profiles
        )
        head_text = f"THM741/{tag}-global-head-ledger/v1\n" + "".join(
            "E="
            + ",".join(map(str, row["body"]))
            + f";K={row['K']};R={fraction_text(row['residual'])};head="
            + ",".join(map(str, sorted(speed for _, speed in row["head"])))
            + "\n"
            for row in band_profiles
        )
        # Compatibility note: the KLE10 v1 ledger predates the dangerous-sum
        # prune and its historical header names all exhaustive head controls
        # "dangerous".  The separate exact count above records the true
        # predicate; retaining this byte string preserves its audited digest.
        candidate_text = f"THM741/{tag}-dangerous-quadruples/v1\n" + "".join(
            line
            for closure in band_closures
            for line in closure["candidate_rows"]
        )
        literal_text = f"THM741/{tag}-literal-endpoint-ledger/v1\n" + "".join(
            line
            for closure in band_closures
            for line in closure["literal_rows"]
        )
        return tuple(
            hashlib.sha256(text.encode()).hexdigest()
            for text in (closed_text, head_text, candidate_text, literal_text)
        )

    k10_digests = band_digests(4, 10, "KLE10")
    require(
        k10_digests
        == (
            EXPECTED_K10_CLOSED_BODY_DIGEST,
            EXPECTED_K10_HEAD_DIGEST,
            EXPECTED_K10_CANDIDATE_DIGEST,
            EXPECTED_K10_LITERAL_DIGEST,
        ),
        "K<=10 canonical ledgers changed",
    )
    k11_20_digests = band_digests(11, 20, "K11-20")
    require(
        k11_20_digests
        == (
            EXPECTED_K11_20_CLOSED_BODY_DIGEST,
            EXPECTED_K11_20_HEAD_DIGEST,
            EXPECTED_K11_20_CANDIDATE_DIGEST,
            EXPECTED_K11_20_LITERAL_DIGEST,
        ),
        "K=11..20 canonical ledgers changed",
    )
    k21_50_digests = band_digests(21, 50, "K21-50")
    require(
        k21_50_digests
        == (
            EXPECTED_K21_50_CLOSED_BODY_DIGEST,
            EXPECTED_K21_50_HEAD_DIGEST,
            EXPECTED_K21_50_CANDIDATE_DIGEST,
            EXPECTED_K21_50_LITERAL_DIGEST,
        ),
        "K=21..50 canonical ledgers changed",
    )
    k51_100_digests = band_digests(51, 100, "K51-100")
    require(
        k51_100_digests
        == (
            EXPECTED_K51_100_CLOSED_BODY_DIGEST,
            EXPECTED_K51_100_HEAD_DIGEST,
            EXPECTED_K51_100_CANDIDATE_DIGEST,
            EXPECTED_K51_100_LITERAL_DIGEST,
        ),
        "K=51..100 canonical ledgers changed",
    )
    k101_200_digests = band_digests(101, 200, "K101-200")
    require(
        k101_200_digests
        == (
            EXPECTED_K101_200_CLOSED_BODY_DIGEST,
            EXPECTED_K101_200_HEAD_DIGEST,
            EXPECTED_K101_200_CANDIDATE_DIGEST,
            EXPECTED_K101_200_LITERAL_DIGEST,
        ),
        "K=101..200 canonical ledgers changed",
    )
    k201_500_digests = band_digests(201, 500, "K201-500")
    require(
        k201_500_digests
        == (
            EXPECTED_K201_500_CLOSED_BODY_DIGEST,
            EXPECTED_K201_500_HEAD_DIGEST,
            EXPECTED_K201_500_CANDIDATE_DIGEST,
            EXPECTED_K201_500_LITERAL_DIGEST,
        ),
        "K=201..500 canonical ledgers changed",
    )
    k501_1000_digests = band_digests(501, 1000, "K501-1000")
    require(
        k501_1000_digests
        == (
            EXPECTED_K501_1000_CLOSED_BODY_DIGEST,
            EXPECTED_K501_1000_HEAD_DIGEST,
            EXPECTED_K501_1000_CANDIDATE_DIGEST,
            EXPECTED_K501_1000_LITERAL_DIGEST,
        ),
        "K=501..1000 canonical ledgers changed",
    )
    k1001_10000_digests = band_digests(1001, 10000, "K1001-10000")
    require(
        k1001_10000_digests
        == (
            EXPECTED_K1001_10000_CLOSED_BODY_DIGEST,
            EXPECTED_K1001_10000_HEAD_DIGEST,
            EXPECTED_K1001_10000_CANDIDATE_DIGEST,
            EXPECTED_K1001_10000_LITERAL_DIGEST,
        ),
        "K=1001..10000 canonical ledgers changed",
    )
    k10001_plus_digests = band_digests(
        10001,
        max(int(row["K"]) for row in selected),
        "K10001-plus",
    )
    require(
        k10001_plus_digests
        == (
            EXPECTED_K10001_PLUS_CLOSED_BODY_DIGEST,
            EXPECTED_K10001_PLUS_HEAD_DIGEST,
            EXPECTED_K10001_PLUS_CANDIDATE_DIGEST,
            EXPECTED_K10001_PLUS_LITERAL_DIGEST,
        ),
        "K>=10001 canonical ledgers changed",
    )

    selected_set = {row["body"] for row in selected}
    flood_overlap = []
    for edge in combinations(range(1, 8), 2):
        body = tuple(sorted((*FLOOD_CORE, *edge)))
        if body in selected_set:
            flood_overlap.append(edge)
    require(tuple(flood_overlap) == EXPECTED_FLOOD_OVERLAP, "flood overlap changed")

    positive_set = {row["body"] for row in positive}
    feasible_set = {row["body"] for row in feasible}
    impossible_set = {row["body"] for row in impossible}
    flood_bodies = {
        tuple(sorted((*FLOOD_CORE, *edge)))
        for edge in combinations(range(1, 8), 2)
    }
    extra_flood = impossible_set & flood_bodies
    require(len(extra_flood) == 6, "rank-impossible flood count changed")
    proved_union = positive_set | feasible_set | extra_flood
    residual_set = impossible_set - extra_flood
    require(
        len(proved_union) == 1406 and len(residual_set) == 596,
        "proved/residual root partition changed",
    )

    def body_digest(bodies: set[tuple[int, ...]]) -> str:
        text = "\n".join(
            ",".join(map(str, body)) for body in BODIES if body in bodies
        )
        return hashlib.sha256(text.encode()).hexdigest()

    feasible_digest = body_digest(feasible_set)
    impossible_digest = body_digest(impossible_set)
    proved_digest = body_digest(proved_union)
    residual_digest = body_digest(residual_set)
    require(
        (
            feasible_digest,
            impossible_digest,
            proved_digest,
            residual_digest,
        )
        == (
            EXPECTED_FEASIBLE_BODY_DIGEST,
            EXPECTED_IMPOSSIBLE_BODY_DIGEST,
            EXPECTED_PROVED_UNION_DIGEST,
            EXPECTED_RESIDUAL_BODY_DIGEST,
        ),
        "global root-partition ledgers changed",
    )
    classification_text = (
        "THM741/all-root-ranked-head-classification/v1\n"
        + "".join(
            "E="
            + ",".join(map(str, body))
            + ";class="
            + (
                "positive"
                if body in positive_set
                else "feasible"
                if body in feasible_set
                else "impossible"
            )
            + f";flood={int(body in flood_bodies)}\n"
            for body in BODIES
        )
    )
    classification_digest = hashlib.sha256(
        classification_text.encode()
    ).hexdigest()
    require(
        classification_digest == EXPECTED_CLASSIFICATION_DIGEST,
        "global classification ledger changed",
    )

    maximum_threshold = max(selected, key=lambda row: row["threshold"])
    maximum_tail_first = max(int(row["tail_first"]) for row in selected)
    require(
        (
            maximum_threshold["body"],
            maximum_threshold["threshold"],
            maximum_tail_first,
        )
        == EXPECTED_MAXIMUM_THRESHOLD,
        "maximum selected tail threshold changed",
    )
    maximum_head_speed = max(speed for row in selected for _, speed in row["head"])
    require(
        maximum_head_speed == EXPECTED_MAXIMUM_HEAD_SPEED,
        "maximum selected head speed changed",
    )

    print("THM-741 ALL RANK-FEASIBLE RANKED-HEAD ATLAS")
    print("status=FINITE-EXACT+GLOBAL-HEAD-SEALED")
    print("root_universe=C(14,9)=2002; added_speeds=15<=a<b<c<d")
    print(
        "top4="
        f"positive:{len(positive)},nonpositive:{len(nonpositive)};"
        f"rank_feasible:{len(feasible)},rank_impossible:{len(impossible)}"
    )
    print(
        "selected_K_histogram="
        + ",".join(f"{size}:{count}" for size, count in histogram.items())
    )
    print(f"selected_roots={len(selected)}")
    print(
        f"union_bound_dangerous={dangerous};"
        f"literal_carrier_checks={literal_checks};"
        f"positive={literal_checks};tight=0"
    )
    print(
        "minimum_survivor="
        + fraction_text(global_minimum[0])
        + ":family="
        + ",".join(map(str, global_minimum[1]))
        + f":components={global_minimum[2]}"
    )
    print(
        "minimum_component_count="
        f"{global_minimum_components[0]}:"
        + fraction_text(global_minimum_components[1])
        + ":family="
        + ",".join(map(str, global_minimum_components[2]))
    )
    print(
        "maximum_selected_tail_threshold="
        + ",".join(map(str, maximum_threshold["body"]))
        + ":"
        + fraction_text(maximum_threshold["threshold"])
        + f";tail_first_max={maximum_tail_first};maximum_head_speed={maximum_head_speed}"
    )
    print(
        "flood_overlap="
        + ";".join(",".join(map(str, edge)) for edge in flood_overlap)
        + ";new_beyond_595=811;proved_union=1406/2002;residual=596"
    )
    print("KLE10_digests=" + ",".join(k10_digests))
    print("K11_20_digests=" + ",".join(k11_20_digests))
    print("K21_50_digests=" + ",".join(k21_50_digests))
    print("K51_100_digests=" + ",".join(k51_100_digests))
    print("K101_200_digests=" + ",".join(k101_200_digests))
    print("K201_500_digests=" + ",".join(k201_500_digests))
    print("K501_1000_digests=" + ",".join(k501_1000_digests))
    print("K1001_10000_digests=" + ",".join(k1001_10000_digests))
    print("K10001_plus_digests=" + ",".join(k10001_plus_digests))
    print(
        "global_partition_digests="
        f"feasible:{feasible_digest},impossible:{impossible_digest},"
        f"proved:{proved_digest},residual:{residual_digest},"
        f"classification:{classification_digest}"
    )
    print("scope=816 ranked-head roots; global THM-741 and LRC14 remain open")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
