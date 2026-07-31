#!/usr/bin/env python3
"""Exact geometric audit of the sharpest THM-2892 discrepancy row.

For a connected residual component I=[a,b], a danger comb D_w can cover I
up to measure zero only if I is contained in the closure of one of its teeth.
Equivalently, some integer k must satisfy

    w*b - 1/14 <= k <= w*a + 1/14.

Thus noncontainment is the exact integer condition

    ceil(w*b - 1/14) > floor(w*a + 1/14).

In particular a singleton cover of a residual having a component of length
ell requires w*ell <= 1/7.  This gives a geometric finite horizon independent
of the discrepancy tail used by the main locked replay.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAIN_PATH = (
    ROOT
    / "04-computation/lrc14_j5_postgate_heavy_triangle_census_codex_20260729.py"
)
MAIN_SHA256 = (
    "0a386864dee44144130d25060301ed0f6d8c3cd02136b6aecf58efa3ae3a790d"
)

BODY = (1, 5, 6, 7, 8, 9, 11, 13)
APEX = 19
TRIANGLE = (37, 121, 130)
EXPECTED_MASS = F(1_099_999, 24_444_420)
EXPECTED_COMPONENTS = 42
EXPECTED_LONGEST = (F(911, 220_220), F(939, 1_820), F(881, 1_694))
EXPECTED_BOUND = F(31_460, 911)
EXPECTED_HORIZON = 34
EXPECTED_CANDIDATES = (27,)
EXPECTED_FINITE_MAX = (F(7_111, 459_800), 25)
EXPECTED_FINITE_MARGIN = F(137_171_599, 4_644_439_800)
EXPECTED_CANDIDATE_27 = (F(77_101, 8_368_360), F(33_241_751, 928_887_960))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_main():
    require(
        hashlib.sha256(MAIN_PATH.read_bytes()).hexdigest() == MAIN_SHA256,
        "locked THM-2892 verifier changed",
    )
    spec = importlib.util.spec_from_file_location("j5_component_audit_main", MAIN_PATH)
    require(spec is not None and spec.loader is not None, "cannot load main verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ceil_fraction(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def closed_tooth_can_contain(a: F, b: F, speed: int) -> bool:
    return (
        ceil_fraction(speed * b - F(1, 14))
        <= floor_fraction(speed * a + F(1, 14))
    )


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    module = load_main()
    root_good, _, _ = module.T.CORE.good_norm(BODY)
    residual = module.T.THM2883.subtract_local_multi(
        root_good,
        (APEX, *TRIANGLE),
    )
    direct, components, mass = module.T.CORE.good_norm(
        tuple(sorted((*BODY, APEX, *TRIANGLE)))
    )
    require(
        residual == direct
        and mass == module.interval_mass(residual) == EXPECTED_MASS
        and components == len(residual) == EXPECTED_COMPONENTS,
        "literal residual reconstruction changed",
    )

    longest = max((right - left, left, right) for left, right in residual)
    require(longest == EXPECTED_LONGEST, "longest component changed")
    ell, left, right = longest
    bound = F(1, 1) / (7 * ell)
    horizon = floor_fraction(bound)
    require(
        bound == EXPECTED_BOUND and horizon == EXPECTED_HORIZON,
        "geometric singleton horizon changed",
    )
    require(
        F(1, 7 * (horizon + 1)) < ell
        and F(1, 7 * horizon) >= ell,
        "horizon equality boundary changed",
    )

    candidates = tuple(
        speed
        for speed in range(15, horizon + 1)
        if speed not in {APEX, *TRIANGLE}
        and closed_tooth_can_contain(left, right, speed)
    )
    require(candidates == EXPECTED_CANDIDATES, "containment candidates changed")

    speeds = [
        speed
        for speed in range(15, horizon + 1)
        if speed not in {APEX, *TRIANGLE}
    ]
    coverages = module.T.THM2885.coverages_many(residual, speeds)
    require(
        all(
            value == module.T.THM2885.coverage(residual, speed)
            for value, speed in coverages
        ),
        "vector/scalar coverage control failed",
    )
    finite_max = min(
        ((-value, speed) for value, speed in coverages),
        key=lambda item: (item[0], item[1]),
    )
    finite_max = (-finite_max[0], finite_max[1])
    require(finite_max == EXPECTED_FINITE_MAX, "finite maximum changed")
    margin = mass - finite_max[0]
    require(margin == EXPECTED_FINITE_MARGIN, "finite margin changed")

    coverage_27 = module.T.THM2885.coverage(residual, 27)
    require(
        (coverage_27, mass - coverage_27) == EXPECTED_CANDIDATE_27,
        "speed-27 candidate control changed",
    )

    print("LRC14 J=5 SINGLETON COMPONENT SEAL AUDIT")
    print(f"body={BODY};apex={APEX};triangle={TRIANGLE}")
    print(f"residual_mass={ftext(mass)};components={components}")
    print(
        f"longest_component={ftext(ell)};"
        f"left={ftext(left)};right={ftext(right)}"
    )
    print(
        f"singleton_speed_bound={ftext(bound)};"
        f"geometric_horizon={horizon};"
        f"closed_tooth_candidates={','.join(map(str, candidates))}"
    )
    print(
        f"finite_head_max={ftext(finite_max[0])};"
        f"maximizer={finite_max[1]};"
        f"exact_finite_margin={ftext(margin)}"
    )
    print(
        f"candidate_27_coverage={ftext(coverage_27)};"
        f"candidate_27_margin={ftext(mass - coverage_27)}"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
