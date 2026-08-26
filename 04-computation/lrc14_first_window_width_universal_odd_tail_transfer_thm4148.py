#!/usr/bin/env python3
"""Exact arithmetic controls for THM-4148's first-window transfer."""

from __future__ import annotations

import hashlib
from fractions import Fraction
from math import comb


EXPECTED_SEMANTIC_SHA256 = "6f4c7b6fd62bca62f0a8baccb63822e6f19df831d2ece22001c2412cfce44e2e"
EXPECTED_SEMANTIC_FNV64 = "1f67e51174709013"


def gap(speed: int, phase: Fraction) -> Fraction:
    value = speed * phase
    residue = value - value.numerator // value.denominator
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Fraction) -> Fraction:
    return min(gap(speed, phase) for speed in speeds)


def width(first: int, last: int) -> Fraction:
    return Fraction(13, 14 * last) - Fraction(1, 14 * first)


def width_gate(first: int, last: int) -> bool:
    return 27 * (13 * first - last) >= 4 * first * last


def fnv64(data: bytes) -> str:
    value = 14695981039346656037
    for byte in data:
        value ^= byte
        value = value * 1099511628211 & ((1 << 64) - 1)
    return f"{value:016x}"


def main() -> None:
    delta = Fraction(1, 14)
    gate = Fraction(2, 189)

    assert 3 * gate == Fraction(2, 63)
    assert Fraction(2, 7 * 27) == gate

    residual_min = min(
        Fraction(1, 2) - Fraction(r, 28 * m)
        for m in range(3, 501)
        for r in range(1, 26, 2)
    )
    assert residual_min == Fraction(17, 84) > delta

    valid_46: list[int] = []
    for first in range(3, 1001):
        if width_gate(first, first + 45):
            valid_46.append(first)
        assert not width_gate(first, first + 46)
    assert valid_46 == list(range(14, 23))
    assert 140 * 140 - 16 * 1242 == -272
    assert 144 * 144 - 16 * 1215 == 1296

    block = tuple(range(15, 61))
    body_phase_left = Fraction(1, 210)
    body_phase_right = Fraction(13, 840)
    block_width = body_phase_right - body_phase_left
    assert block_width == Fraction(3, 280)
    assert block_width - gate == Fraction(1, 7560)
    for phase in (body_phase_left, body_phase_right):
        assert clearance(block, phase) >= delta

    residual_phase = Fraction(211, 420)
    residual_tail_min = min(gap(r, residual_phase) for r in range(1, 26, 2))
    assert residual_tail_min == Fraction(37, 84)
    assert clearance(tuple(2 * h for h in block), residual_phase) >= delta

    low_lift = Fraction(1, 420)
    high_lift = Fraction(211, 420)
    assert gap(1, low_lift) == Fraction(1, 420)
    assert gap(211, high_lift) == Fraction(1, 420)
    assert (211 * 211) % 420 == 1

    rescue_body_phase = Fraction(1, 105)
    rescue_phase = Fraction(53, 105)
    assert 2 * rescue_phase - 1 == rescue_body_phase
    superbody = tuple(2 * h for h in block) + (1, 211)
    rescue_clearance = clearance(superbody, rescue_phase)
    assert rescue_clearance == Fraction(1, 7)
    assert gap(1, rescue_phase) == Fraction(52, 105)
    assert gap(211, rescue_phase) == Fraction(52, 105)

    family_count = comb(46, 11)
    assert family_count == 13_340_783_196

    ledger = (
        "gate=2/189;scale3=2/63;q27=2/189;"
        "residual=17/84;valid46=14..22;max47disc=-272;"
        "block15_60=3/280;surplus=1/7560;"
        "clock=211/420;clock_tail=37/84;"
        "hostile=1/420,1/420;rescue=1/7;"
        "families=13340783196"
    )
    semantic = hashlib.sha256(ledger.encode()).hexdigest()
    semantic_fnv = fnv64(ledger.encode())
    assert semantic == EXPECTED_SEMANTIC_SHA256
    assert semantic_fnv == EXPECTED_SEMANTIC_FNV64

    print("THM4148_FIRST_WINDOW_WIDTH_ODD_TAIL_TRANSFER_20260825")
    print("status=PASS;scope=exact rational arithmetic controls")
    print("width_gate=2/189;scale_gates=(2/63,2/189)")
    print(f"residual_clock_min={residual_min}")
    print("consecutive_46_starts=(14,15,16,17,18,19,20,21,22)")
    print("consecutive_47_discriminant=-272;consecutive_46_discriminant=1296")
    print(f"block_15_60_interval=({body_phase_left},{body_phase_right},{block_width})")
    print(f"block_15_60_surplus={block_width-gate}")
    print(f"residual_clock=(211/420,{residual_tail_min})")
    print("moving_base_hostile=((1/420,1),(211/420,211),1/420)")
    print(f"moving_base_rescue=(1/105,53/105,{rescue_clearance})")
    print(f"eleven_subsets={family_count}")
    print(f"semantic_sha256={semantic}")
    print(f"semantic_fnv64={semantic_fnv}")


if __name__ == "__main__":
    main()
