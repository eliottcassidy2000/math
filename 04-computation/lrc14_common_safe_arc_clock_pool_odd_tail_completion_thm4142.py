#!/usr/bin/env python3
"""Exact primary audit for THM-4142.

The THM-4136 two-sheet theorem uses only one closed body-safe arc and three
low-ratio clocks.  This audit computes the complete positive-integer alphabet
sharing those four certificates, then checks that every eleven-subset inherits
the infinite cross-comb proof and all sixty-eight residual ratios.
"""

from __future__ import annotations

import ast
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from math import comb, gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/lrc14_fixed_body_universal_odd_tail_completion_thm4136.py"
EXPECTED_BASE_SHA256 = "50c40de4756362d05ffe7f6064a0990e062a08b8b317708febd554e780de8c23"
EXPECTED_SEMANTIC = "fe13698f821b5328bf7fc56267c0745280a7607261770f1b99424627f6686a3e"

DELTA = Q(1, 14)
J = (Q(33, 70), Q(27, 56))
CLOCKS = (Q(89, 1176), Q(181, 4704), Q(431, 4480))
J_POOL = (
    1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22,
    24, 26, 28, 30, 32, 33, 35, 37, 39, 41, 43, 45, 47, 49, 60, 62, 64, 66,
)
POOL = (
    1, 3, 4, 6, 8, 9, 10, 11, 12, 14, 15, 16, 18, 22, 24, 28, 30, 32,
    35, 37, 41, 43, 45, 49, 60, 64,
)
EXCLUDED_J_SAFE = (5, 7, 13, 20, 26, 33, 39, 47, 62, 66)
DIRECT_CASES = (
    ((1, 4, 6, 8, 10, 12, 14, 15, 16, 18, 22), (1, 9)),
    ((1, 3, 4, 6, 8, 9, 10, 11, 12, 14, 15), (13, 23)),
    ((28, 30, 32, 35, 37, 41, 43, 45, 49, 60, 64), (1, 13)),
    ((1, 4, 8, 10, 12, 15, 18, 24, 30, 35, 41), (1, 27)),
    ((3, 6, 9, 11, 14, 22, 28, 37, 45, 60, 64), (3, 27)),
    ((28, 30, 32, 35, 37, 41, 43, 45, 49, 60, 64), (2001, 18009)),
)

CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def file_digest(path: Path) -> str:
    hasher = sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


require(file_digest(BASE_PATH) == EXPECTED_BASE_SHA256, "pinned THM-4136 source")
SPEC = importlib.util.spec_from_file_location("thm4136_base_for_4142", BASE_PATH)
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)


def pair(value: Q) -> tuple[int, int]:
    return value.numerator, value.denominator


def fmt(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def circle_distance(value: Q) -> Q:
    residue = value % 1
    return min(residue, 1 - residue)


def interval_integer_distance(left: Q, right: Q) -> Q:
    """Minimum distance from the real interval [left,right] to Z."""
    require(left <= right, ("ordered interval", left, right))
    first_integer = -((-left.numerator) // left.denominator)
    if Q(first_integer) <= right:
        return Q(0)
    return min(circle_distance(left), circle_distance(right))


def body_interval_minimum(speed: int) -> Q:
    require(speed > 0, ("positive speed", speed))
    return interval_integer_distance(speed * J[0], speed * J[1])


def common_pool_audit() -> dict[str, object]:
    require(J[1] - J[0] == Q(3, 280), "body arc length")
    # For s>=94 the image interval has length >1 and therefore contains an
    # integer.  Thus 1..93 is the complete finite universe.
    require(94 * (J[1] - J[0]) > 1, "all larger speeds hit an integer")
    j_pool = tuple(speed for speed in range(1, 94) if body_interval_minimum(speed) >= DELTA)
    require(j_pool == J_POOL, "complete J-safe alphabet")
    common = tuple(
        speed for speed in j_pool
        if all(circle_distance(2 * speed * phase) > DELTA for phase in CLOCKS)
    )
    require(common == POOL, "complete common safe-arc/clock pool")
    require(tuple(speed for speed in j_pool if speed not in common) == EXCLUDED_J_SAFE,
            "clock-excluded J-safe speeds")

    interval_minima = tuple((speed, pair(body_interval_minimum(speed))) for speed in common)
    clock_minima = []
    clock_owners = []
    for phase in CLOCKS:
        rows = tuple((circle_distance(2 * speed * phase), speed) for speed in common)
        minimum = min(value for value, _ in rows)
        require(minimum > DELTA, ("strict common clock margin", phase, minimum))
        clock_minima.append(pair(minimum))
        clock_owners.append(tuple(speed for value, speed in rows if value == minimum))
    require(tuple(clock_minima) == ((4, 49), (11, 147), (17, 224)),
            "common body clock minima")
    require(tuple(clock_owners) == ((60,), (64,), (10,)), "common body clock owners")

    excluded_clock_rows = tuple(
        (
            speed,
            tuple(pair(circle_distance(2 * speed * phase)) for phase in CLOCKS),
        )
        for speed in EXCLUDED_J_SAFE
    )
    require(comb(len(common), 11) == 7_726_160, "body count")
    return {
        "finite_universe": (1, 93),
        "large_speed_gate": (94, pair(94 * (J[1] - J[0]))),
        "J_pool": j_pool,
        "common_pool": common,
        "excluded_J_safe": excluded_clock_rows,
        "interval_minima": interval_minima,
        "clock_minima": tuple(clock_minima),
        "clock_owners": tuple(clock_owners),
        "body_count": comb(len(common), 11),
    }


def clock_index(p: int, q: int) -> int:
    require(0 < p < q <= 25 and p % 2 == q % 2 == 1 and gcd(p, q) == 1,
            ("residual ratio", p, q))
    if 13 not in (p, q):
        return 0
    if (p, q) in ((1, 13), (13, 25)):
        return 2
    return 1


def quotient_and_residual_audit(pool: dict[str, object]) -> dict[str, object]:
    scan_count = 0
    scan_max = Q(0)
    scan_maximizers: list[tuple[int, int]] = []
    for q in range(3, 102, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            scan_count += 1
            components = BASE.two_sheet_components(p, q)
            beta = max((right - left for left, right in components), default=Q(0))
            if q >= 9:
                require(beta <= Q(2, 7 * q), ("cross-comb q-tooth bound", p, q, beta))
            require(beta <= Q(2, 63), ("universal cross-comb bound", p, q, beta))
            if beta > scan_max:
                scan_max = beta
                scan_maximizers = [(p, q)]
            elif beta == scan_max:
                scan_maximizers.append((p, q))
    require(scan_count == 1053, "quotient hostile count")
    require(scan_max == Q(2, 63) and tuple(scan_maximizers) == ((1, 9),),
            "sharp quotient hostile")
    require(3 * (J[1] - J[0]) - Q(2, 63) == Q(1, 2520), "odd gcd gate")
    require((J[1] - J[0]) - Q(2, 189) == Q(1, 7560), "large q gate")

    pairs = tuple(
        (p, q)
        for q in range(3, 26, 2)
        for p in range(1, q, 2)
        if gcd(p, q) == 1
    )
    require(len(pairs) == 68, "residual ratio count")
    groups: list[list[tuple[int, int]]] = [[], [], []]
    tail_minima = [Q(1), Q(1), Q(1)]
    for p, q in pairs:
        index = clock_index(p, q)
        phase = CLOCKS[index]
        tail_gap = min(circle_distance(p * phase), circle_distance(q * phase))
        require(tail_gap > DELTA, ("tail clock", p, q, phase, tail_gap))
        groups[index].append((p, q))
        tail_minima[index] = min(tail_minima[index], tail_gap)
    require(tuple(len(group) for group in groups) == (56, 10, 2), "residual categories")
    require(tuple(tail_minima) == (Q(89, 1176), Q(541, 4704), Q(431, 4480)),
            "tail-only minima")
    body_minima = tuple(Q(*entry) for entry in pool["clock_minima"])
    full_minima = tuple(min(tail, body) for tail, body in zip(tail_minima, body_minima))
    require(full_minima == (Q(89, 1176), Q(11, 147), Q(17, 224)),
            "all-body residual minima")
    require(all(value > DELTA for value in full_minima), "strict residual closure")
    return {
        "cross_comb_scan": (scan_count, 101, pair(scan_max), tuple(scan_maximizers)),
        "scale_gates": (pair(Q(1, 2520)), pair(Q(1, 7560))),
        "residual_count": len(pairs),
        "category_counts": tuple(len(group) for group in groups),
        "tail_minima": tuple(pair(value) for value in tail_minima),
        "full_minima": tuple(pair(value) for value in full_minima),
        "pair_digest": digest(pairs),
    }


def direct_full_row(body: tuple[int, ...], tails: tuple[int, int]) -> tuple[Q, Q, int]:
    require(len(body) == len(set(body)) == 11, ("eleven-body", body))
    require(set(body) <= set(POOL), ("body in common pool", body))
    a, b = tails
    require(0 < a < b and a % 2 == b % 2 == 1, ("odd tails", tails))
    speeds = tuple(2 * speed for speed in body) + tails
    require(len(speeds) == len(set(speeds)) == 13, ("distinct full row", body, tails))
    walls = BASE.strict_danger_walls(speeds)
    for left, right in zip(walls, walls[1:]):
        phase = (left + right) / 2
        gap = BASE.clearance(speeds, phase)
        if gap > DELTA:
            return phase, gap, len(walls)
    raise RuntimeError(f"no direct full-row witness for {body}, {tails}")


def main() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "no assertions")
    require(tuple(CLOCKS) == tuple(BASE.CLOCKS), "inherited clock bank")
    require(J == tuple(BASE.J), "inherited body arc")
    pool = common_pool_audit()
    quotient = quotient_and_residual_audit(pool)
    direct = tuple(
        (body, tails, pair(phase), pair(gap), wall_count)
        for body, tails in DIRECT_CASES
        for phase, gap, wall_count in (direct_full_row(body, tails),)
    )
    ledger = {
        "theorem": "THM-4142",
        "statement": (
            "for every eleven-subset B of the 26-speed common pool and every distinct "
            "positive odd a,b, 2B union {a,b} is 1/14-safe"
        ),
        "base_sha256": EXPECTED_BASE_SHA256,
        "delta": pair(DELTA),
        "J": (pair(J[0]), pair(J[1]), pair(J[1] - J[0])),
        "clocks": tuple(pair(value) for value in CLOCKS),
        "pool": pool,
        "quotient": quotient,
        "direct_controls": direct,
        "scope": (
            "7,726,160 bodies in one common certificate pool; bodies outside the pool, "
            "physical entry, mixed tail parity, and LRC14 remain open"
        ),
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")
    print("THM4142_COMMON_SAFE_ARC_CLOCK_POOL_ODD_TAIL_COMPLETION_20260825")
    print("status=PASS")
    print(f"base_sha256={EXPECTED_BASE_SHA256}")
    print(f"J=[{fmt(J[0])},{fmt(J[1])}];length={fmt(J[1]-J[0])};clocks={tuple(fmt(x) for x in CLOCKS)}")
    print(f"J_pool_count={len(J_POOL)};common_pool_count={len(POOL)};common_pool={POOL}")
    print(f"excluded_J_safe={EXCLUDED_J_SAFE};body_count={comb(len(POOL),11)}")
    print(f"clock_body_minima={pool['clock_minima']};clock_body_owners={pool['clock_owners']}")
    print(f"cross_comb_scan={quotient['cross_comb_scan']};scale_gates={quotient['scale_gates']}")
    print(f"residual_categories={quotient['category_counts']};tail_minima={quotient['tail_minima']};full_minima={quotient['full_minima']}")
    print(f"direct_controls={direct}")
    print("scope=common 26-speed certificate pool only; arbitrary bodies, entry, mixed parity, and LRC14 open")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
