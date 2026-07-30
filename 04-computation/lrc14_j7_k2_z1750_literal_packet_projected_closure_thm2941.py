#!/usr/bin/env python3
"""Exact projected k=2 closure of the first-drift slice z1=1750.

This verifier starts from all C(14,6)=3003 literal six-body carriers.  It
reconstructs the rigorous scalar upper envelope over every later label, rather
than importing an exploratory survivor list.  Four rows whose projected wall
forces a high later drift are then closed by the exact residue-ray law

    delta_E(r+mL) = A_E(r)/(r+mL).

The eight remaining exact-suffix rows are reconstructed over every denominator
state and passed through the all-divisor/common-K5 status screen.  Four rows
are status-empty.  Exactly 100 states remain on four bodies.

For each remaining state, its positive scalar slack gives a finite literal
label gate on every residue ray.  All resulting distinct five-drift packets
are enumerated.  For every packet Z the program computes the exact lossless
projected residual

    P_(E,Z) = phi_L(C_E minus union_(z in Z) D_z).

If its mass is at least 25/91, it cannot be contained in the union of the two
remaining aligned danger combs.  A direct global interval subtraction and
projection independently checks the cellwise De Morgan computation on the
minimum-margin packet of every surviving body.

There is no finite label horizon in the exact-ray, status, packet, or
projected-residual stages.  The horizon 7000 occurs only in the rigorous global
scalar upper envelope, whose four omitted slots are bounded separately by
6r/[49(H+1)].
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCALAR_BASE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_scalar_band_1836_1836_thm2941.py"
)
PROJECTED = (
    ROOT
    / "04-computation"
    / "lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
)
RAY_STATUS = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_z1932_ray_status_closure_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1750_literal_packet_projected_closure_thm2941.out"
)
EXPECTED_SCALAR_BASE_SHA256 = (
    "99a8170413b34c2887eb7deac1ab038e6492c74b6d2f36ef8b9d0127afccd19e"
)
EXPECTED_PROJECTED_SHA256 = (
    "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662"
)
EXPECTED_RAY_STATUS_SHA256 = (
    "77c848a2b18b6fd0f60b31bd1a2b9c9972793e3fa1fe3838633b46211a64d037"
)

FIRST = 1750
ALIGNED_TWO_UNION_CAP = F(25, 91)
HIGH_WALL_RATIO = F(13, 150)
SCALAR_ETA = F(1, 91)
SUFFIX_SLOTS = 4

FORCED_HIGH_BODIES = (
    (1, 2, 10, 11, 12, 14),
    (1, 2, 10, 12, 13, 14),
    (1, 4, 10, 12, 13, 14),
    (1, 6, 10, 12, 13, 14),
)
EXACT_BODIES = (
    (1, 4, 6, 8, 10, 14),
    (1, 4, 8, 10, 12, 14),
    (1, 6, 8, 10, 12, 14),
    (2, 4, 6, 8, 10, 14),
    (2, 4, 8, 10, 12, 14),
    (2, 6, 8, 10, 12, 14),
    (3, 4, 8, 10, 12, 14),
    (4, 6, 8, 10, 12, 14),
)
RESIDUAL_BODIES = (
    (1, 4, 8, 10, 12, 14),
    (2, 4, 8, 10, 12, 14),
    (2, 6, 8, 10, 12, 14),
    (4, 6, 8, 10, 12, 14),
)
EXPECTED_SCALAR_KEYS = tuple(
    (body, FIRST) for body in sorted((*FORCED_HIGH_BODIES, *EXACT_BODIES))
)
EXPECTED_GLOBAL_PROFILE_SHA256 = (
    "b0c4283b61cf20c19a327537b82ee701bc43e68ead1930d564f33070d26854a2"
)
EXPECTED_GLOBAL_SURVIVOR_SHA256 = (
    "b2bf3a4df0857587074a58194a7bbfae52b47c04298d2d4a842437675b6f8ba1"
)
EXPECTED_HIGH_GAPS = {
    FORCED_HIGH_BODIES[0]: F(-122114680703, 813986523850500),
    FORCED_HIGH_BODIES[1]: F(-5167406441, 47765952848250),
    FORCED_HIGH_BODIES[2]: F(-27564560809, 95531905696500),
    FORCED_HIGH_BODIES[3]: F(-23299726123, 95531905696500),
}

# scalar, scalar digest, crude kills, crude digest, status kills, status
# digest, survivors, survivor digest.  The empty crude stage is intentional at
# z1=1750; the stronger common-status table supplies all pre-projection kills.
EXPECTED_EXACT_STAGES = {
    EXACT_BODIES[0]: (
        1,
        "3902ea5927438aa75fb23971f9e006d6cdff076f4694b0e8047239c5e7ba605c",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        1,
        "37de3656c0dbd7dddcf11a3fd45acfb701e322b92489515b16b4da8782daada8",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ),
    EXACT_BODIES[1]: (
        474,
        "a4547c928f9fc7602e3cf9dbee2122e0d76fc4deed95ee187b99a1381415f321",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        471,
        "ebc7448c7d8b3196159c5abdfded0a4e799754c8c13c29fbad6b3f3c2c18ff51",
        3,
        "6f36b3ba5b851a76a629e75aaa479ead9caf0491bd0424112712c99376e80c6b",
    ),
    EXACT_BODIES[2]: (
        7,
        "c77a2141e0b3d3224b51a1e2fd438c7604df2df8ec35327ae80994ef2b204694",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        7,
        "2d5a5128006454e1f3696ba0d29c79b71c95b818ef8b181d92d5e5530094e7a1",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ),
    EXACT_BODIES[3]: (
        5,
        "5e9df93ae08046a6d699dc3e12676e92e499e85c179cd790856c745a9075a844",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        5,
        "0c734a03ed364ec7f0101157b9ea6b5720a152818f6cd739d5da457ac01e5d5d",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ),
    EXACT_BODIES[4]: (
        186,
        "0731c408fd674d6af106980c022b14a208eacaef27b219a6a110e5aa786528ef",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        94,
        "810f82af8d5398d037dbcc97d73872f1fbbaa38e7aa4f92e75a88ee130088997",
        92,
        "424ca96366625136417d1b705b514cf92b35b141e234b8fd5a9213515c90c334",
    ),
    EXACT_BODIES[5]: (
        5,
        "da5b0e592ed6a23467172cd250b3cc3264cff69d14483722ad0131b637397f55",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        1,
        "790bb3c706ca2fa6b2670349cbb60df25bcd9d12bfe5a743916ff2b593235f95",
        4,
        "91a2f8f9d4a97a346ad7c49559d3fa5d57dd2fa0981d4fa174c992d9994a9370",
    ),
    EXACT_BODIES[6]: (
        2,
        "e8aec5003ad194517745dda8f01b5543f5cb3e006e46b1f8c2ec3101c94eea9b",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        2,
        "fa65b693b9011270ca20740e0681c1b07eab3f3b7ed262f34f5bb43b586df633",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ),
    EXACT_BODIES[7]: (
        2,
        "c026e698cd5e161bef160db161cfba6e74294591469e75e03f6a7796559a0619",
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        1,
        "742fabb93b9eeb620506f69b224bb2b7ecf57be8ef5c356e4f30a7bce0a8bf51",
        1,
        "475eef6bca4067a12390d49befb62b600b73d8088733be0fb026ff317fc83a74",
    ),
}

EXPECTED_PACKET_SUMMARIES = (
    (
        RESIDUAL_BODIES[0],
        (
            4,
            F(79223, 234325),
            2,
            F(2591645947, 9539068470750),
            2,
            2,
            (
                F(79223, 234325),
                (1750, 1805, 1836, 2060, 2172),
                F(1578, 2575),
                2,
            ),
            F(1),
            "7ab12308da5f054968b4ea854d7b9eb046d8aad65ce70c81ef2a952dfa16c37d",
            "754956cff396ef141b8c6d5c120cae283a265293a5a1e0d98409307fb69c68f7",
        ),
    ),
    (
        RESIDUAL_BODIES[1],
        (
            219,
            F(4085, 54691),
            2,
            F(2968487, 6174154350),
            7,
            7,
            (
                F(4085, 54691),
                (1750, 1810, 1836, 2060, 2404),
                F(210, 601),
                1,
            ),
            F(1),
            "2d79d6a74888c29f5ba824c864b0af77f6b0089791406ac8ec6019126332baff",
            "e913b30536def41c03729ef57044d0ee78fd6ccab5164c615e42976f6c506598",
        ),
    ),
    (
        RESIDUAL_BODIES[2],
        (
            5,
            F(723886, 3632447),
            2,
            F(2027349047, 2860691515500),
            2,
            2,
            (
                F(723886, 3632447),
                (1750, 1784, 1790, 1810, 1836),
                F(18921, 39917),
                2,
            ),
            F(1),
            "904a17a45c71df67542d914c8543dd13d761ee8f20273afe4d00039d59018439",
            "3660c5a2d777e4d3855936ae03b3819b0b547f2d5c497a9dd6dde058285a6b7e",
        ),
    ),
    (
        RESIDUAL_BODIES[3],
        (
            1,
            F(8215, 16471),
            1,
            F(7094632, 9261231525),
            2,
            1,
            (
                F(8215, 16471),
                (1750, 1810, 1836, 2060, 2172),
                F(140, 181),
                1,
            ),
            F(1),
            "e704b6f78edb390cd4a8740e1cacf76251440fef944fcc821d4d6e7815977420",
            "d5d8de04f260aa92e06635d4637058210aeed91bd1dae4e1266e1404f2f8ac33",
        ),
    ),
)
EXPECTED_PROFILE_SHA256 = (
    "0c911c7222dc6686735e878a2905a0dcf42d3bd00ecc56c7835b5447741969da"
)
EXPECTED_SEMANTIC_SHA256 = (
    "71dfba440f4573d82de685adb7f4684909e99b7fafcfe729215ac2fc0573e972"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def load_module(path: Path, name: str, expected_sha256: str):
    require(file_sha256(path) == expected_sha256, (path.name, "dependency changed"))
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (path.name, "cannot load"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load_module(SCALAR_BASE, "z1750_scalar_base", EXPECTED_SCALAR_BASE_SHA256)
P = load_module(PROJECTED, "z1750_projected_support", EXPECTED_PROJECTED_SHA256)
K = load_module(RAY_STATUS, "z1750_ray_status_support", EXPECTED_RAY_STATUS_SHA256)
U = K.U
require(
    B.RULER == U.suffix.A.RULER == P.A.RULER,
    "dependency rulers disagree",
)
B.START = B.END = FIRST
B.LABELS = np.arange(FIRST, B.HORIZON + 1, dtype=np.int64)


def global_scalar_profile(body: tuple[int, ...]):
    return B.profile(body)


def projected_safe_lower_bound(
    cells: tuple[int, ...],
    period: int,
    labels: tuple[int, ...],
    stop_at_cap: bool = True,
):
    common_danger = ((F(0), F(1)),)
    cells_used = 0
    for cell in cells:
        local_union = P.merge_fraction(
            [
                interval
                for label in labels
                for interval in P.phase_danger(cell, label, period)
            ]
        )
        common_danger = P.intersect_fraction(common_danger, local_union)
        cells_used += 1
        if stop_at_cap and P.interval_mass(common_danger) <= 1 - ALIGNED_TWO_UNION_CAP:
            break
    return 1 - P.interval_mass(common_danger), cells_used, common_danger


def direct_projected_mass(
    body: tuple[int, ...], period: int, labels: tuple[int, ...]
) -> F:
    carrier = tuple(
        (F(left, P.A.RULER), F(right, P.A.RULER))
        for left, right in P.A.carrier_for(body)
    )
    removed = P.merge_fraction(
        [
            interval
            for label in labels
            for interval in P.danger_fraction(label)
        ]
    )
    residual = P.subtract_fraction(carrier, removed)
    projected = []
    for left, right in residual:
        scaled_left = period * left
        scaled_right = period * right
        for integer in range(P.floor_fraction(scaled_left), P.ceil_fraction(scaled_right)):
            piece_left = max(scaled_left, F(integer)) - integer
            piece_right = min(scaled_right, F(integer + 1)) - integer
            if piece_left < piece_right:
                projected.append((piece_left, piece_right))
    return P.interval_mass(P.merge_fraction(projected))


def ray_data(body: tuple[int, ...], carrier, h: F, period: int):
    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, period):
        amplitude = residue * U.delta(carrier, h, residue)
        require(
            (residue + period) * U.delta(carrier, h, residue + period) == amplitude,
            (body, "ray recurrence failed", residue),
        )
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(
            amplitudes[period - residue] == -amplitudes[residue]
            for residue in range(1, period)
        ),
        (body, "ray antipode failed"),
    )
    return tuple(amplitudes), tuple(sorted(signs.items())), sha256(
        repr(tuple(amplitudes)).encode()
    ).hexdigest()


def first_strictly_after(residue: int, bound: int, period: int) -> int:
    if residue > bound:
        return residue
    return residue + ((bound - residue) // period + 1) * period


def first_at_least(residue: int, bound: int, period: int) -> int:
    if residue >= bound:
        return residue
    return residue + ((bound - residue + period - 1) // period) * period


def forced_high_profile(body: tuple[int, ...]):
    carrier = U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    period = 14 * lcm(*body)
    high_floor = max(15, (HIGH_WALL_RATIO * period).numerator // (HIGH_WALL_RATIO * period).denominator + 1)
    require(FIRST < high_floor, (body, "row does not carry a forced-high wall"))
    amplitudes, signs, ray_digest = ray_data(body, carrier, h, period)

    arbitrary = []
    high = []
    fifth = []
    for residue in range(1, period):
        amplitude = amplitudes[residue]
        if amplitude <= 0:
            continue
        first_label = first_strictly_after(residue, FIRST, period)
        arbitrary.extend(
            (amplitude / (first_label + offset * period), first_label + offset * period, residue)
            for offset in range(SUFFIX_SLOTS)
        )
        fifth.append((amplitude / (first_label + SUFFIX_SLOTS * period), first_label + SUFFIX_SLOTS * period, residue))
        high_label = first_at_least(residue, high_floor, period)
        high.append((amplitude / high_label, high_label, residue))
    arbitrary.sort(key=lambda row: (-row[0], row[1], row[2]))
    high.sort(key=lambda row: (-row[0], row[1], row[2]))
    fifth.sort(key=lambda row: (-row[0], row[1], row[2]))
    rank4 = tuple(arbitrary[:4])
    require(
        len(rank4) == 4
        and rank4[-1][0] > 0
        and fifth[0][0] <= rank4[-1][0]
        and all(row[1] < high_floor for row in rank4),
        (body, "finite top-four ray gate failed"),
    )
    best_high = high[0]
    constrained = (*rank4[:3], best_high)
    require(len({row[1] for row in constrained}) == 4, (body, "label collision"))
    first_delta = amplitudes[FIRST % period] / FIRST
    upper = first_delta + sum((row[0] for row in constrained), F())
    lower = h * SCALAR_ETA
    gap = upper - lower
    require(gap == EXPECTED_HIGH_GAPS[body] < 0, (body, "forced-high gap changed", gap))
    packet = (FIRST, *sorted(row[1] for row in constrained))
    cells = P.body_cells(P.A.carrier_for(body), period)
    prefix_mass, prefix_cells, _ = projected_safe_lower_bound(cells, period, packet)
    require(prefix_mass >= ALIGNED_TWO_UNION_CAP, (body, "projected control failed"))
    return (
        body,
        h,
        len(carrier),
        period,
        high_floor,
        signs,
        ray_digest,
        rank4,
        fifth[0],
        best_high,
        first_delta,
        tuple(constrained),
        upper,
        lower,
        gap,
        packet,
        prefix_mass,
        prefix_cells,
    )


def exact_status_frontier(body, carrier, h, lower, period, first_delta, first_d):
    amplitudes, signs, ray_digest = ray_data(body, carrier, h, period)
    divisors = tuple(d for d in U.support.divisors(period) if d > 1)

    @lru_cache(maxsize=None)
    def top_values(d: int, multiplicity: int):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (period // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = first_strictly_after(residue, FIRST, period)
            candidates.extend(
                (amplitude / (label + offset * period), label + offset * period, residue)
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        require(len(candidates) >= multiplicity, (body, "missing ray", d))
        return tuple(candidates[:multiplicity])

    scalar = []
    for tail in combinations_with_replacement(divisors, SUFFIX_SLOTS):
        upper = first_delta
        labels = []
        for d, multiplicity in Counter(tail).items():
            chosen = top_values(d, multiplicity)
            upper += sum((value for value, _label, _residue in chosen), F())
            labels.extend(label for _value, label, _residue in chosen)
        if upper >= lower:
            scalar.append((tuple(sorted((first_d, *tail))), upper, (FIRST, *sorted(labels))))
    require(len({row[0] for row in scalar}) == len(scalar), (body, "duplicate state"))

    actual_period, ranges = U.support.safe_cell_ranges(body)
    require(actual_period == period, (body, "safe-cell ruler changed"))
    arcs_cache = {}
    crude_kills = []
    crude_survivors = []
    for ds, upper, labels in scalar:
        D = lcm(*ds)
        arcs = arcs_cache.setdefault(D, U.fibre.projected_support_arcs(D, ranges))
        witness = None
        for q in U.support.divisors(D):
            histogram = U.fibre.residue_load_histogram(arcs, q)
            target = max(load for load, count in histogram if count)
            capacity = sum(U.fibre_cap(D, d, q) for d in ds)
            if target > capacity:
                witness = (q, D // q, target, capacity)
                break
        row = (ds, upper, labels, witness)
        (crude_survivors if witness is None else crude_kills).append(row)

    status_kills = []
    states = []
    for ds, upper, labels, _ in crude_survivors:
        D = lcm(*ds)
        arcs = arcs_cache[D]
        witness = None
        for M in U.support.divisors(D):
            q = D // M
            marginals, capacities = K.hunter_status_data5(D, ds, q)
            histogram = U.fibre.residue_load_histogram(arcs, q)
            feasible, certificate = K.common_status_feasible5(q, marginals, capacities, histogram)
            if not feasible:
                require(certificate is not None, (body, ds, "missing Farkas witness"))
                witness = (
                    q,
                    M,
                    marginals,
                    tuple(sorted(set(capacities))),
                    histogram,
                    certificate,
                )
                break
        row = (ds, upper, labels, witness)
        (states if witness is None else status_kills).append(row)

    stage_rows = (tuple(scalar), tuple(crude_kills), tuple(status_kills), tuple(states))
    stage_digests = tuple(sha256(repr(rows).encode()).hexdigest() for rows in stage_rows)
    actual = (
        len(scalar),
        stage_digests[0],
        len(crude_kills),
        stage_digests[1],
        len(status_kills),
        stage_digests[2],
        len(states),
        stage_digests[3],
    )
    require(actual == EXPECTED_EXACT_STAGES[body], (body, "status frontier changed", actual))
    return (
        amplitudes,
        signs,
        ray_digest,
        len(divisors),
        comb(len(divisors) + SUFFIX_SLOTS - 1, SUFFIX_SLOTS),
        stage_rows,
        stage_digests,
    )


def close_literal_packets(body, period, amplitudes, cells, lower, first_delta, states):
    @lru_cache(maxsize=None)
    def ray_heads(d: int, multiplicity: int):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (period // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = first_strictly_after(residue, FIRST, period)
            candidates.extend(
                (amplitude / (label + offset * period), label + offset * period, residue)
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        return tuple(candidates[:multiplicity])

    @lru_cache(maxsize=None)
    def eligible_group(d: int, multiplicity: int, slack: F):
        heads = ray_heads(d, multiplicity)
        top_sum = sum((row[0] for row in heads), F())
        threshold = heads[-1][0] - slack
        require(threshold > 0, (body, "nonfinite slack gate", d, multiplicity, slack))
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (period // d) * direction
            amplitude = amplitudes[residue]
            if amplitude <= 0:
                continue
            label = first_strictly_after(residue, FIRST, period)
            while amplitude / label >= threshold:
                candidates.append((amplitude / label, label, residue))
                label += period
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        choices = []
        for chosen in combinations(candidates, multiplicity):
            labels = tuple(sorted(row[1] for row in chosen))
            if len(set(labels)) != multiplicity:
                continue
            value = sum((row[0] for row in chosen), F())
            deficit = top_sum - value
            if deficit <= slack:
                choices.append((deficit, labels, value))
        choices.sort()
        require(choices and choices[0][0] == 0, (body, "maximizer lost", d, multiplicity))
        return top_sum, threshold, tuple(candidates), tuple(choices)

    packet_audit = []
    state_audit = []
    minimum_row = None
    maximum_prefix_cells = 0
    for state_index, (ds, upper, maximizing_labels, witness) in enumerate(states, 1):
        require(witness is None, (body, ds, "status survivor carries witness"))
        slack = upper - lower
        tail = list(ds)
        tail.remove(period // gcd(period, FIRST))
        groups = []
        for d, multiplicity in sorted(Counter(tail).items()):
            top_sum, threshold, candidates, choices = eligible_group(d, multiplicity, slack)
            groups.append((d, multiplicity, top_sum, threshold, len(candidates), choices))
        require(
            first_delta + sum((group[2] for group in groups), F()) == upper,
            (body, ds, "state envelope changed"),
        )
        state_packets = 0
        state_minimum = None
        for selection in product(*(group[5] for group in groups)):
            deficit = sum((choice[0] for choice in selection), F())
            if deficit > slack:
                continue
            labels = (FIRST, *sorted(label for choice in selection for label in choice[1]))
            require(len(labels) == 5 and len(set(labels)) == 5, (body, "nondistinct packet", labels))
            actual_upper = first_delta + sum((choice[2] for choice in selection), F())
            require(actual_upper == upper - deficit >= lower, (body, "packet scalar gate", labels))
            projected_mass, used, _ = projected_safe_lower_bound(cells, period, labels)
            margin = projected_mass - ALIGNED_TWO_UNION_CAP
            require(margin >= 0, (body, "projected survivor", state_index, ds, labels, margin))
            packet_audit.append((state_index, ds, labels, actual_upper, projected_mass, used))
            state_packets += 1
            maximum_prefix_cells = max(maximum_prefix_cells, used)
            state_minimum = margin if state_minimum is None else min(state_minimum, margin)
            candidate = (margin, labels, projected_mass, used)
            if minimum_row is None or candidate < minimum_row:
                minimum_row = candidate
        require(state_packets > 0, (body, ds, "state has no literal packet"))
        state_audit.append(
            (
                state_index,
                ds,
                upper,
                maximizing_labels,
                slack,
                tuple(
                    (d, multiplicity, top_sum, threshold, candidate_count, len(choices))
                    for d, multiplicity, top_sum, threshold, candidate_count, choices in groups
                ),
                state_packets,
                state_minimum,
            )
        )
    require(minimum_row is not None, (body, "minimum packet missing"))
    require(
        len({row[2] for row in packet_audit}) == len(packet_audit),
        (body, "literal packet appeared in two denominator states"),
    )
    minimum_threshold = min(
        group[3] for state in state_audit for group in state[5]
    )
    maximum_ray_candidates = max(
        group[4] for state in state_audit for group in state[5]
    )
    maximum_group_choices = max(
        group[5] for state in state_audit for group in state[5]
    )
    require(minimum_threshold > 0, (body, "literal ray universe is not finite"))
    _margin, control_labels, prefix_mass, _used = minimum_row
    full_mass, full_cells, _ = projected_safe_lower_bound(cells, period, control_labels, False)
    direct_mass = direct_projected_mass(body, period, control_labels)
    require(
        full_cells == len(cells) and direct_mass == full_mass and direct_mass >= prefix_mass,
        (body, "direct projected subtraction control failed"),
    )
    return (
        len(packet_audit),
        minimum_row[0],
        maximum_prefix_cells,
        minimum_threshold,
        maximum_ray_candidates,
        maximum_group_choices,
        minimum_row,
        direct_mass,
        sha256(repr(tuple(packet_audit)).encode()).hexdigest(),
        sha256(repr(tuple(state_audit)).encode()).hexdigest(),
        tuple(state_audit),
    )


def exact_body_profile(body: tuple[int, ...]):
    carrier = U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * U.suffix.ETAS[2]
    period = 14 * lcm(*body)
    require(period == 11760 and FIRST >= 1020, (body, "ordinary exact geometry changed"))
    first_delta = U.delta(carrier, h, FIRST)
    first_d = period // gcd(period, FIRST)
    amplitudes, signs, ray_digest, divisor_count, trials, stage_rows, stage_digests = exact_status_frontier(
        body, carrier, h, lower, period, first_delta, first_d
    )
    states = stage_rows[3]
    packet_summary = None
    cell_count = 0
    if states:
        require(body in RESIDUAL_BODIES, (body, "unexpected residual body"))
        cells = P.body_cells(P.A.carrier_for(body), period)
        cell_count = len(cells)
        packet_summary = close_literal_packets(
            body, period, amplitudes, cells, lower, first_delta, states
        )
    else:
        require(body not in RESIDUAL_BODIES, (body, "expected residual body vanished"))
    return (
        body,
        h,
        len(carrier),
        period,
        lower,
        first_delta,
        first_d,
        signs,
        ray_digest,
        divisor_count,
        trials,
        tuple(len(rows) for rows in stage_rows),
        stage_digests,
        cell_count,
        packet_summary,
    )


def pool_map(function, rows, workers):
    if workers == 1:
        return [function(row) for row in rows]
    with mp.get_context("spawn").Pool(workers) as pool:
        return list(pool.imap(function, rows, chunksize=4 if function is global_scalar_profile else 1))


def render(global_profiles, high_profiles, exact_profiles) -> str:
    require(len(global_profiles) == comb(14, 6) == 3003, "global body universe changed")
    global_profiles = tuple(sorted(global_profiles, key=lambda row: row[0]))
    scalar_survivors = tuple(
        sorted(
            (survivor for row in global_profiles for survivor in row[10]),
            key=lambda row: (row[0], row[1]),
        )
    )
    scalar_keys = tuple((row[0], row[1]) for row in scalar_survivors)
    require(sum(row[7] for row in global_profiles) == 3003, "candidate-row count changed")
    require(scalar_keys == EXPECTED_SCALAR_KEYS, ("global survivor keys changed", scalar_keys))
    global_profile_hash = sha256(repr(global_profiles).encode()).hexdigest()
    global_survivor_hash = sha256(repr(scalar_survivors).encode()).hexdigest()
    require(global_profile_hash == EXPECTED_GLOBAL_PROFILE_SHA256, "global profile digest changed")
    require(global_survivor_hash == EXPECTED_GLOBAL_SURVIVOR_SHA256, "global survivor digest changed")

    high_profiles = tuple(sorted(high_profiles, key=lambda row: row[0]))
    exact_profiles = tuple(sorted(exact_profiles, key=lambda row: row[0]))
    require(tuple(row[0] for row in high_profiles) == FORCED_HIGH_BODIES, "high universe changed")
    require(tuple(row[0] for row in exact_profiles) == EXACT_BODIES, "exact universe changed")
    residual_profiles = tuple(row for row in exact_profiles if row[14] is not None)
    require(tuple(row[0] for row in residual_profiles) == RESIDUAL_BODIES, "residual universe changed")
    packet_summaries = tuple((row[0], row[14][:-1]) for row in residual_profiles)
    if EXPECTED_PACKET_SUMMARIES is not None:
        require(packet_summaries == EXPECTED_PACKET_SUMMARIES, ("packet summaries changed", packet_summaries))
    status_survivors = sum(row[11][3] for row in exact_profiles)
    literal_packets = sum(row[14][0] for row in residual_profiles)
    minimum_margin = min(row[14][1] for row in residual_profiles)
    maximum_prefix = max(row[14][2] for row in residual_profiles)

    profile_payload = (
        EXPECTED_SCALAR_BASE_SHA256,
        EXPECTED_PROJECTED_SHA256,
        EXPECTED_RAY_STATUS_SHA256,
        FIRST,
        ALIGNED_TWO_UNION_CAP,
        global_profile_hash,
        global_survivor_hash,
        high_profiles,
        exact_profiles,
    )
    profile_hash = sha256(repr(profile_payload).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "complete profile digest changed")
    semantic_payload = (
        FIRST,
        tuple(range(1, 15)),
        6,
        5,
        2,
        B.HORIZON,
        HIGH_WALL_RATIO,
        SCALAR_ETA,
        ALIGNED_TWO_UNION_CAP,
        scalar_keys,
        tuple((row[0], row[14]) for row in exact_profiles),
        status_survivors,
        literal_packets,
        minimum_margin,
        maximum_prefix,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 z1=1750 literal-packet projected closure",
        f"scalar_base_sha256={file_sha256(SCALAR_BASE)}",
        f"projected_source_sha256={file_sha256(PROJECTED)}",
        f"ray_status_source_sha256={file_sha256(RAY_STATUS)}",
        (
            "universe=all C(14,6)=3003 bodies;first nonaligned label exactly 1750;"
            "four distinct later drifts;two aligned labels"
        ),
        (
            f"global_scalar_horizon={B.HORIZON};omitted_slot_bound=6r/[49(H+1)];"
            f"global_scalar_survivors={len(scalar_survivors)}"
        ),
        (
            "ray_law=delta(r+mL)=A(r)/(r+mL);forced_high_wall=floor(13L/150)+1;"
            "scalar_necessary=sum(delta)>=h/91"
        ),
        (
            "projected_identity=P_(E,Z)=phi_L(C_E minus union_z D_z);"
            "two_aligned_open_union_cap=25/91"
        ),
        (
            f"forced_high_rows={len(high_profiles)};forced_high_scalar_empty={len(high_profiles)};"
            f"exact_suffix_rows={len(exact_profiles)};status_survivors={status_survivors}"
        ),
        (
            f"all_scalar_eligible_literal_packets={literal_packets};projected_kills={literal_packets};"
            f"survivors=0;global_minimum_margin={ftext(minimum_margin)};"
            f"max_prefix_cells={maximum_prefix}"
        ),
    ]
    for row in high_profiles:
        lines.append(
            f"HIGH;E={','.join(map(str, row[0]))};h={ftext(row[1])};r={row[2]};L={row[3]};"
            f"high_floor={row[4]};ray_signs={row[5]};ray_sha256={row[6]};"
            f"upper={ftext(row[12])};lower={ftext(row[13])};gap={ftext(row[14])};"
            f"maximizing_packet={row[15]};projected_control_mass={ftext(row[16])};"
            f"projected_prefix_cells={row[17]};conclusion=SCALAR-EMPTY"
        )
    for row in exact_profiles:
        packet = row[14]
        line = (
            f"EXACT;E={','.join(map(str, row[0]))};h={ftext(row[1])};r={row[2]};L={row[3]};"
            f"lower={ftext(row[4])};delta1={ftext(row[5])};first_d={row[6]};"
            f"ray_signs={row[7]};ray_sha256={row[8]};denominators={row[9]};trials={row[10]};"
            f"scalar={row[11][0]};crude={row[11][1]};status={row[11][2]};"
            f"residual_states={row[11][3]};stage_sha256={row[12]}"
        )
        if packet is None:
            line += ";conclusion=STATUS-EMPTY"
        else:
            line += (
                f";body_cells={row[13]};literal_packets={packet[0]};"
                f"projected_kills={packet[0]};min_margin={ftext(packet[1])};"
                f"max_prefix_cells={packet[2]};min_ray_cutoff={ftext(packet[3])};"
                f"max_ray_candidates={packet[4]};max_group_choices={packet[5]};"
                f"minimum_row={packet[6]};direct_control_mass={ftext(packet[7])};"
                f"packet_sha256={packet[8]};state_audit_sha256={packet[9]};"
                "conclusion=PROJECTED-EMPTY"
            )
        lines.append(line)
    lines.extend(
        (
            "independent_control=direct interval subtraction/projection agrees with full cell De Morgan identity on each residual body's minimum row",
            "conclusion=the complete projected k=2 first-drift slice z1=1750 is empty",
            f"global_profile_sha256={global_profile_hash}",
            f"global_survivor_sha256={global_survivor_hash}",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(4, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    bodies = tuple(combinations(range(1, 15), 6))
    global_profiles = pool_map(global_scalar_profile, bodies, args.workers)
    high_profiles = pool_map(forced_high_profile, FORCED_HIGH_BODIES, args.workers)
    exact_profiles = pool_map(exact_body_profile, EXACT_BODIES, args.workers)
    output = render(global_profiles, high_profiles, exact_profiles)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
