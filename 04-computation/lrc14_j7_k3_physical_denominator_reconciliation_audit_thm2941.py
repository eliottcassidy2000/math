#!/usr/bin/env python3
"""Independent controls for the LRC14 k=3 physical denominator bridge.

The primary bridge extracts exact-lcm coefficients by divisor-Mobius
inversion.  This audit instead propagates the current lcm while adjoining
each possible later denominator symbol.  It compares the full retained
feature distribution, not merely the total coefficient.  It also replays a
deterministic collection of physical bodies through both the vectorized
bridge and THM-2941's original scalar implementation.
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import lcm
from pathlib import Path


HERE = Path(__file__).resolve().parent
BRIDGE_PATH = HERE / "lrc14_j7_k3_physical_denominator_reconciliation_thm2941.py"
EXPECTED_BRIDGE_SHA256 = (
    "9cf4cab69d9fa6f65ece80160eb244a5773ff7c7943ca81697e0236150da6033"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def load_bridge():
    require(
        hashlib.sha256(BRIDGE_PATH.read_bytes()).hexdigest()
        == EXPECTED_BRIDGE_SHA256,
        "physical denominator reconciliation changed",
    )
    spec = importlib.util.spec_from_file_location("physical_denominator_bridge", BRIDGE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load bridge")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load_bridge()


def lcm_state_distribution(D: int, d1: int) -> Counter:
    """Direct symbol-by-symbol multiset DP with current lcm retained."""
    fixed = B.denominator_feature(D, d1)
    # state=(later_used,current_lcm,m2,m3,m4,large,uniform)
    states = Counter({(0, d1, *fixed): 1})
    for d in (value for value in B.divisors_of(D) if value > 1):
        unit = B.denominator_feature(D, d)
        updated = Counter(states)
        for state, multiplicity in states.items():
            used, current, m2, m3, m4, large, uniform = state
            for copies in range(1, 4 - used):
                updated[
                    (
                        used + copies,
                        lcm(current, d),
                        m2 + copies * unit[0],
                        m3 + copies * unit[1],
                        m4 + copies * unit[2],
                        large + copies * unit[3],
                        uniform + copies * unit[4],
                    )
                ] += multiplicity
        states = updated
    result = Counter()
    for state, multiplicity in states.items():
        used, current, m2, m3, m4, large, uniform = state
        if used == 3 and current == D:
            result[((m2, m3, m4), large, uniform)] += multiplicity
    return result


def recurrence_controls() -> tuple[int, str]:
    cases = 0
    semantic = hashlib.sha256()
    for D in (14, 28, 42, 56, 84, 168, 420, 840, 2520, 11760):
        alphabet = tuple(d for d in B.divisors_of(D) if d > 1)
        # All fixed denominators at small/moderate D; a deterministic spanning
        # sample at the two largest controls keeps the audit quick.
        if D <= 840:
            fixed_denominators = alphabet
        else:
            fixed_denominators = tuple(
                sorted({alphabet[0], alphabet[len(alphabet) // 3], alphabet[-2], alphabet[-1]})
            )
        for d1 in fixed_denominators:
            direct = lcm_state_distribution(D, d1)
            mobius = B.conditional_feature_distribution(D, d1)
            require(direct == mobius, ("lcm-state/Mobius mismatch", D, d1))
            semantic.update(f"{D}|{d1}|{tuple(sorted(direct.items()))}\n".encode())
            cases += 1
    return cases, semantic.hexdigest()


def physical_controls() -> tuple[int, str]:
    bodies = tuple(combinations(range(1, 15), 6))
    indices = (0, 1, 17, 149, 377, 751, 1201, 1801, 2401, 3002)
    semantic = hashlib.sha256()
    for index in indices:
        body = bodies[index]
        vectorized = B.physical_rows_for_body(body)
        scalar = scalar_physical_rows(body)
        require(vectorized[2] == scalar, ("exact physical row mismatch", body))
        semantic.update(
            f"{body}|{vectorized[1]}|{vectorized[2]}|{vectorized[3:]}\n".encode()
        )
    return len(indices), semantic.hexdigest()


def scalar_physical_rows(body: tuple[int, ...]) -> tuple[int, ...]:
    """Literal Fraction replay of the canonical projected-suffix loop."""
    carrier = B.physical.A.carrier_for(body)
    h = Q(sum(right - left for left, right in carrier), B.physical.A.RULER)
    components = len(carrier)
    L = 14 * math.lcm(*body)
    delta = {
        label: B.physical.A.singleton_coverage(carrier, label) - h / 7
        for label in range(B.physical.A.BASE_LABEL, B.physical.HORIZON + 1)
        if label % L
    }
    bound = Q(24 * components, 49) / (h * B.ETA3)
    cap = bound.numerator // bound.denominator
    wall = B.PROJECTED_RATIO3 * L
    high_floor = max(15, wall.numerator // wall.denominator + 1)
    ordinary_tail = Q(6 * components, 49 * (B.physical.HORIZON + 1))
    high_tail = Q(
        6 * components,
        49 * max(B.physical.HORIZON + 1, high_floor),
    )
    arbitrary_top = []
    high_top = []
    survivors = []
    for first in range(B.physical.HORIZON, B.physical.A.BASE_LABEL - 1, -1):
        if first % L == 0:
            continue
        if first <= cap:
            suffix = B.physical.suffix_upper(
                arbitrary_exact=arbitrary_top,
                high_exact=high_top,
                need=3,
                tail=ordinary_tail,
                high_tail=high_tail,
                constrained=first < high_floor,
            )
            upper = delta[first] + sum((value for value, _label, _kind in suffix), Q(0))
            if upper >= h * B.ETA3:
                survivors.append(first)
        item = (delta[first], first, "EXACT")
        B.physical.top_insert(arbitrary_top, item, limit=3)
        if first >= high_floor:
            B.physical.top_insert(high_top, item, limit=3)
    return tuple(sorted(survivors))


def main() -> None:
    recurrence_cases, recurrence_semantic = recurrence_controls()
    physical_cases, physical_semantic = physical_controls()
    print("LRC14 k3 physical denominator bridge independent audit")
    print(f"bridge_sha256={hashlib.sha256(BRIDGE_PATH.read_bytes()).hexdigest()}")
    print(f"lcm_state_recurrence_cases={recurrence_cases}")
    print(f"lcm_state_semantic_sha256={recurrence_semantic}")
    print(f"canonical_physical_replay_cases={physical_cases}")
    print(f"physical_semantic_sha256={physical_semantic}")
    print("all_independent_controls=PASS")


if __name__ == "__main__":
    main()
