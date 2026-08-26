#!/usr/bin/env python3
"""Clean-room exact referee for THM-4142.

This file imports neither the primary audit nor THM-4136.  It reconstructs
the safe alphabet from real-interval/integer separation, rebuilds the
two-sheet quotient from strict wall cells, verifies every residual ratio, and
solves a disjoint bank of complete thirteen-speed rows directly.
"""

from __future__ import annotations

import ast
from fractions import Fraction as F
from hashlib import sha256
import json
from math import ceil, comb, gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

DELTA = F(1, 14)
LEFT = F(33, 70)
RIGHT = F(27, 56)
CLOCKS = (F(89, 1176), F(181, 4704), F(431, 4480))
EXPECTED_POOL = (
    1, 3, 4, 6, 8, 9, 10, 11, 12, 14, 15, 16, 18, 22, 24, 28, 30, 32,
    35, 37, 41, 43, 45, 49, 60, 64,
)
EXPECTED_SEMANTIC = "dfddeafede97f47d6b8c38a4676cd6bc24c16559833690f51400addae8a5801b"
DIRECT_ROWS = (
    ((3, 6, 8, 9, 11, 16, 22, 28, 37, 49, 64), (5, 7)),
    ((1, 10, 12, 14, 18, 24, 30, 32, 41, 45, 60), (13, 25)),
    ((4, 6, 8, 15, 16, 22, 28, 35, 43, 49, 64), (101, 103)),
    ((1, 3, 9, 11, 14, 18, 24, 32, 37, 45, 60), (999, 1001)),
)

CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def pair(value: F) -> tuple[int, int]:
    return value.numerator, value.denominator


def digest(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(raw).hexdigest()


def distance(value: F) -> F:
    residue = value % 1
    return residue if residue <= F(1, 2) else 1 - residue


def interval_gap(speed: int) -> F:
    """Exact min_{x in [LEFT,RIGHT]} ||speed*x||."""
    lo = speed * LEFT
    hi = speed * RIGHT
    first = ceil(lo)
    if F(first) <= hi:
        return F(0)
    return min(distance(lo), distance(hi))


def reconstruct_pool() -> dict[str, object]:
    require(94 * (RIGHT - LEFT) > 1, "finite speed cutoff")
    arc_safe = tuple(s for s in range(1, 94) if interval_gap(s) >= DELTA)
    require(len(arc_safe) == 36, "arc-safe count")
    pool = tuple(
        s for s in arc_safe
        if min(distance(2 * s * x) for x in CLOCKS) > DELTA
    )
    require(pool == EXPECTED_POOL, "independent common-pool reconstruction")
    failure_masks = tuple(
        (s, tuple(distance(2 * s * x) <= DELTA for x in CLOCKS))
        for s in arc_safe if s not in pool
    )
    require(all(any(mask) for _, mask in failure_masks), "every exclusion witnessed")
    margins = tuple(
        min(distance(2 * s * x) for s in pool)
        for x in CLOCKS
    )
    require(margins == (F(4, 49), F(11, 147), F(17, 224)), "clock margins")
    return {
        "arc_safe": arc_safe,
        "pool": pool,
        "failure_masks": failure_masks,
        "margins": tuple(pair(x) for x in margins),
        "body_count": comb(len(pool), 11),
    }


def bad(speed: int, phase: F) -> bool:
    return distance(speed * phase) < DELTA


def both_roots_bad(p: int, q: int, w: F) -> bool:
    z = (w % 1) / 2
    first = bad(p, z) or bad(q, z)
    second = bad(p, z + F(1, 2)) or bad(q, z + F(1, 2))
    return first and second


def quotient_walls(p: int, q: int) -> tuple[F, ...]:
    walls = {F(0), F(1)}
    for shift in (F(0), F(1, 2)):
        for speed in (p, q):
            for k in range(speed):
                walls.add((2 * ((F(k) - DELTA) / speed - shift)) % 1)
                walls.add((2 * ((F(k) + DELTA) / speed - shift)) % 1)
    return tuple(sorted(walls))


def quotient_components(p: int, q: int) -> tuple[tuple[F, F], ...]:
    require(0 < p < q and p % 2 == q % 2 == 1 and gcd(p, q) == 1,
            ("primitive odd ratio", p, q))
    walls = quotient_walls(p, q)
    cells = tuple(
        (left, right)
        for left, right in zip(walls, walls[1:])
        if both_roots_bad(p, q, (left + right) / 2)
    )
    components: list[tuple[F, F]] = []
    for left, right in cells:
        if components and components[-1][1] == left and both_roots_bad(p, q, left):
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    require(not both_roots_bad(p, q, F(0)), ("zero outside bad quotient", p, q))
    return tuple(components)


def quotient_referee() -> dict[str, object]:
    scanned = 0
    maximum = F(0)
    maximizers: list[tuple[int, int]] = []
    for q in range(3, 152, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            scanned += 1
            components = quotient_components(p, q)
            beta = max((right - left for left, right in components), default=F(0))
            if q >= 9:
                require(beta <= F(2, 7 * q), ("q-tooth bound", p, q, beta))
            require(beta <= F(2, 63), ("uniform quotient bound", p, q, beta))
            if beta > maximum:
                maximum = beta
                maximizers = [(p, q)]
            elif beta == maximum:
                maximizers.append((p, q))
    require(maximum == F(2, 63) and tuple(maximizers) == ((1, 9),),
            "independent sharp quotient")
    hostiles = ((1, 999), (997, 999), (1, 1001), (999, 1001))
    hostile_rows = []
    for p, q in hostiles:
        if gcd(p, q) != 1:
            continue
        beta = max((r - l for l, r in quotient_components(p, q)), default=F(0))
        require(beta <= F(2, 7 * q), ("far quotient hostile", p, q, beta))
        hostile_rows.append(((p, q), pair(beta), pair(F(2, 7 * q))))
    return {
        "scan": (scanned, 151, pair(maximum), tuple(maximizers)),
        "far_hostiles": tuple(hostile_rows),
    }


def residual_referee(pool: tuple[int, ...]) -> dict[str, object]:
    pairs = tuple(
        (p, q)
        for q in range(3, 26, 2)
        for p in range(1, q, 2)
        if gcd(p, q) == 1
    )
    require(len(pairs) == 68, "residual ratio count")
    counts = [0, 0, 0]
    minima = [F(1), F(1), F(1)]
    for p, q in pairs:
        index = 0 if 13 not in (p, q) else (2 if (p, q) in ((1, 13), (13, 25)) else 1)
        phase = CLOCKS[index]
        body_gap = min(distance(2 * s * phase) for s in pool)
        tail_gap = min(distance(p * phase), distance(q * phase))
        gap = min(body_gap, tail_gap)
        require(gap > DELTA, ("all-pool residual clock", p, q, gap))
        counts[index] += 1
        minima[index] = min(minima[index], gap)
    require(tuple(counts) == (56, 10, 2), "residual split")
    require(tuple(minima) == (F(89, 1176), F(11, 147), F(17, 224)),
            "residual full minima")
    return {"counts": tuple(counts), "minima": tuple(pair(x) for x in minima)}


def danger_walls(speeds: tuple[int, ...]) -> tuple[F, ...]:
    walls = {F(0), F(1)}
    for speed in speeds:
        for k in range(speed):
            walls.add(((F(k) - DELTA) / speed) % 1)
            walls.add(((F(k) + DELTA) / speed) % 1)
    return tuple(sorted(walls))


def row_gap(speeds: tuple[int, ...], phase: F) -> F:
    return min(distance(speed * phase) for speed in speeds)


def solve_row(body: tuple[int, ...], tails: tuple[int, int]) -> tuple[F, F, int]:
    require(len(body) == len(set(body)) == 11 and set(body) <= set(EXPECTED_POOL),
            ("direct body", body))
    require(tails[0] < tails[1] and all(x > 0 and x % 2 == 1 for x in tails),
            ("direct tails", tails))
    speeds = tuple(2 * s for s in body) + tails
    walls = danger_walls(speeds)
    for left, right in zip(walls, walls[1:]):
        phase = (left + right) / 2
        gap = row_gap(speeds, phase)
        if gap > DELTA:
            return phase, gap, len(walls)
    raise RuntimeError(f"direct row failed: {body}, {tails}")


def main() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "no assertions")
    pool = reconstruct_pool()
    require(pool["body_count"] == 7_726_160, "eleven-subset count")
    quotient = quotient_referee()
    residual = residual_referee(tuple(pool["pool"]))
    require(3 * (RIGHT - LEFT) > F(2, 63), "gcd scale surplus")
    require(RIGHT - LEFT > F(2, 189), "large primitive-ratio surplus")
    direct = tuple(
        (body, tails, pair(phase), pair(gap), wall_count)
        for body, tails in DIRECT_ROWS
        for phase, gap, wall_count in (solve_row(body, tails),)
    )
    ledger = {
        "theorem": "THM-4142-independent",
        "pool": pool,
        "quotient": quotient,
        "residual": residual,
        "surpluses": (pair(3 * (RIGHT - LEFT) - F(2, 63)),
                      pair((RIGHT - LEFT) - F(2, 189))),
        "direct": direct,
        "scope": "common-certificate pool only; no arbitrary-body or LRC14 closure",
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")
    print("THM4142_COMMON_POOL_INDEPENDENT_AUDIT_20260825")
    print("status=PASS;implementation=clean-room standard-library exact arithmetic")
    print(f"arc_safe_count={len(pool['arc_safe'])};pool={pool['pool']};body_count={pool['body_count']}")
    print(f"failure_masks={pool['failure_masks']};clock_margins={pool['margins']}")
    print(f"quotient_scan={quotient['scan']};far_hostiles={quotient['far_hostiles']}")
    print(f"residual={residual};surpluses={ledger['surpluses']}")
    print(f"direct_rows={direct}")
    print("scope=26-speed common certificate alphabet only; outside bodies, entry, mixed parity, LRC14 open")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
