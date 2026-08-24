#!/usr/bin/env python3
"""Exact THM-3910 auxiliary-center erosion response for THM-3878.

For a scale-one survivor (p,q), an auxiliary multiplier a, and t>=U,
cited LRC through twelve speeds applied to u union {a*t} gives a centre y0
at which all eleven body speeds and a*t have clearance at least 1/13.  Put
w0=t*y0.  The body alone remains 1/14-safe on the closed w-arc of radius
1/182 about w0; importantly, the auxiliary runner need not remain safe on
that entire arc.  A counterexample therefore requires the a-deep centre w0
to lie in the strict radius-1/182 erosion of D_p union D_q.

This program checks that necessary condition by exact rational arithmetic,
freezes the 41 killed scale-one types, and exhausts all other multipliers by
the largest-eroded-component bound.  It also checks the scale-two (1,9)
quotient obstruction.
"""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
import math
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "cyclic_probe",
    ROOT / "04-computation" / "lrc14_eleven_plus_two_cyclic_slack_thm3878.py",
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load THM-3878 cyclic probe")
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)


SCALE1 = (
    (1, 3), (1, 4), (1, 9), (1, 10),
    (2, 3), (2, 9), (2, 15), (2, 21), (2, 23),
    (3, 7), (3, 8), (3, 14), (3, 17), (3, 19), (3, 20),
    (3, 22), (3, 26), (3, 31), (3, 38),
    (4, 7), (4, 13), (4, 19), (4, 21), (4, 25), (4, 37),
    (4, 43), (4, 49), (4, 51),
    (5, 6), (5, 12), (5, 17), (5, 18), (5, 24), (5, 29),
    (5, 36), (5, 39), (5, 41), (5, 42), (5, 48), (5, 53),
    (5, 54), (5, 63),
    (6, 11), (6, 17), (6, 19), (6, 23), (6, 41), (6, 47),
    (6, 53), (6, 65),
    (7, 10), (7, 13), (7, 15), (7, 22),
    (8, 9), (8, 21), (9, 11),
)

RHO = Q(1, 182)
CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def raw_danger(frequency: int) -> list[tuple[Q, Q]]:
    """Open D_frequency pieces clipped to [0,1], cutting the circle at 0."""
    radius = Q(1, 14 * frequency)
    pieces = [(Q(0), radius), (Q(1) - radius, Q(1))]
    for k in range(1, frequency):
        centre = Q(k, frequency)
        pieces.append((centre - radius, centre + radius))
    return pieces


def circle_open_components(frequencies: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    """Lifted open components; the component through zero may end above 1."""
    merged: list[list[Q]] = []
    pieces = sum((raw_danger(frequency) for frequency in frequencies), [])
    for left, right in sorted(pieces):
        if not merged or left >= merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    if len(merged) >= 2 and merged[0][0] == 0 and merged[-1][1] == 1:
        wrap = (merged[-1][0], merged[0][1] + 1)
        middle = tuple((left, right) for left, right in merged[1:-1])
        return (wrap,) + middle
    return tuple((left, right) for left, right in merged)


def has_deep_eroded_centre(
    components: tuple[tuple[Q, Q], ...], auxiliary: int
) -> bool:
    """Test (erosion interior) intersect {||auxiliary*w||>=1/13}."""
    for left, right in components:
        core_left = left + RHO
        core_right = right - RHO
        if core_left >= core_right:
            continue
        k_lo = math.floor(auxiliary * core_left) - 2
        k_hi = math.ceil(auxiliary * core_right) + 2
        for k in range(k_lo, k_hi + 1):
            deep_left = (Q(k) + Q(1, 13)) / auxiliary
            deep_right = (Q(k) + Q(12, 13)) / auxiliary
            # The erosion core is open.  Equality at either core boundary is
            # not enough for a closed radius-RHO arc to lie in the open set.
            if max(core_left, deep_left) < min(core_right, deep_right):
                return True
    return False


def finite_auxiliary_limit(components: tuple[tuple[Q, Q], ...]) -> int:
    """First a forced positive by the largest eroded component length."""
    longest_core = max((right - left - 2 * RHO for left, right in components), default=Q(0))
    if longest_core <= 0:
        return 1
    # If a*longest_core>2/13, the core cannot fit in a connected component
    # of {||a*w||<1/13}; hence it contains an a-deep point.
    return math.floor(Q(2, 13) / longest_core) + 1


def fmt(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    require(len(SCALE1) == 57 and len(set(SCALE1)) == 57, "57 survivor pairs")

    closed = []
    survivors = []
    auxiliary_tests = 0
    for p, q in SCALE1:
        components = circle_open_components((p, q))
        beta = max(right - left for left, right in components)
        require(beta == P.widths(p, q)[0], f"component width {(p, q)}")
        limit = finite_auxiliary_limit(components)
        require(limit >= 1, f"positive search limit {(p, q)}")
        longest_core = beta - 2 * RHO
        require(limit * longest_core > Q(2, 13), f"large-a exhaustion {(p, q)}")
        killers = []
        for auxiliary in range(1, limit + 1):
            auxiliary_tests += 1
            if not has_deep_eroded_centre(components, auxiliary):
                killers.append(auxiliary)
        record = (p, q, beta, limit, tuple(killers))
        if killers:
            closed.append(record)
        else:
            survivors.append(record)

    require(len(closed) == 41 and len(survivors) == 16, "41/16 split")
    require(all(killers == (p,) for p, _, _, _, killers in closed), "canonical a=p killers")

    # Audit the older THM-3878 one-auxiliary cutoff statistic.  Its frozen
    # `at most 78` loop bound is safe, but the exact attained maximum is 77.
    legacy_cutoffs = tuple(
        (p, q, math.floor(Q(13, 7) / beta) + 1)
        for p, q, beta, _, _ in closed + survivors
    )
    legacy_max = max(cutoff for _, _, cutoff in legacy_cutoffs)
    legacy_max_rows = tuple((p, q) for p, q, cutoff in legacy_cutoffs if cutoff == legacy_max)
    require(legacy_max == 77, "THM-3878 exact legacy cutoff")
    require(legacy_max_rows == ((6, 19), (6, 47), (8, 21)), "legacy cutoff equality rows")

    scale2_components = ((Q(2, 21), Q(8, 63)), (Q(55, 63), Q(19, 21)))
    scale2_beta = Q(2, 63)
    require(
        all(right - left == scale2_beta for left, right in scale2_components),
        "scale-two component widths",
    )
    scale2_limit = finite_auxiliary_limit(scale2_components)
    require(
        scale2_limit * (scale2_beta - 2 * RHO) > Q(2, 13),
        "scale-two large-a exhaustion",
    )
    scale2_killers = tuple(
        auxiliary
        for auxiliary in range(1, scale2_limit + 1)
        if not has_deep_eroded_centre(scale2_components, auxiliary)
    )
    require(not scale2_killers, "scale-two has no auxiliary-centre closure")

    semantic = {
        "scope": "THM3878 t>=U auxiliary-centre radius-1/182 erosion",
        "closed": [
            [p, q, fmt(beta), limit, list(killers)]
            for p, q, beta, limit, killers in closed
        ],
        "survivors": [
            [p, q, fmt(beta), limit]
            for p, q, beta, limit, _ in survivors
        ],
        "scale2_limit": scale2_limit,
        "scale2_killers": list(scale2_killers),
        "legacy_cutoff": [legacy_max, [list(pair) for pair in legacy_max_rows]],
    }
    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("LRC14_AUXILIARY_CENTER_EROSION_RESPONSE_THM3910")
    print("scope=THM3878_11+2;t>=U;necessary_certificate;LRC14=OPEN")
    print(f"scale1_input={len(SCALE1)};closed={len(closed)};remaining={len(survivors)}")
    print("closed_by_a=p=" + repr(tuple((p, q) for p, q, *_ in closed)))
    print("remaining=" + repr(tuple((p, q) for p, q, *_ in survivors)))
    print(f"auxiliary_tests={auxiliary_tests};max_limit={max(row[3] for row in closed + survivors)}")
    print(f"scale2_limit={scale2_limit};scale2_killers={scale2_killers}")
    print(f"THM3878_legacy_auxiliary_cutoff_max={legacy_max};rows={legacy_max_rows};at_most_78_remains_safe")
    print("large_a_exhaustion=a*(beta-1/91)>2/13")
    print("semantic_sha256=" + digest)
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
