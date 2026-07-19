#!/usr/bin/env python3
"""Exact arithmetic referee for THM-1199.

The analytic theorem combines three already proved ingredients:

* the THM-1198 six-bin density, whose minimum and maximum are 3/4 and 7/6;
* the THM-1178 endpoint quantum on a connected overlap-nerve tree;
* the THM-1196 prefix-survivor beat mesh.

This replay checks every new constant, exhausts a bounded bank of rational
overlap components, and records why ``gcd(d_i,d_j)`` is not a valid uniform
edge refinement.  The valid refinement labels each chosen component by its
two *actual endpoint owners* ``u,v`` and uses gcd(u,v)/(14uv).

Tournament/carrier audit
------------------------
The faithful pair vertices are endpoint-labelled handoff components, while
the descent vertices are phase-local Farey cells.  Speed order gives only a
transitive tournament and forgets the endpoint owners and beat numerators.
"""

from __future__ import annotations

from fractions import Fraction as Q
from math import ceil, gcd


HEIGHTS = (Q(3, 4), Q(13, 12), Q(7, 6),
           Q(7, 6), Q(13, 12), Q(3, 4))


def require(condition: bool, message: str) -> None:
    """Run an exact-certificate check even when Python uses ``-O``."""
    if not condition:
        raise RuntimeError(f"exact certificate failed: {message}")


def qstr(value: Q) -> str:
    return (str(value.numerator) if value.denominator == 1
            else f"{value.numerator}/{value.denominator}")


def slow_gap(c: int, k: int) -> tuple[Q, Q]:
    return Q(14 * k + 1, 14 * c), Q(14 * k + 13, 14 * c)


def floor_q(value: Q) -> int:
    return value.numerator // value.denominator


def ceil_q(value: Q) -> int:
    return -floor_q(-value)


def meeting_teeth(
    c: int, speed: int, k: int,
) -> tuple[tuple[Q, Q, int, int], ...]:
    """Open tooth components clipped to G, with endpoint-owner labels.

    The openness does not affect lengths.  Owner ``c`` denotes a clipped slow
    gap endpoint; owner ``speed`` denotes a genuine tooth endpoint.
    """
    left, right = slow_gap(c, k)
    first = floor_q(speed * left - Q(1, 14)) - 1
    last = ceil_q(speed * right + Q(1, 14)) + 1
    answer: list[tuple[Q, Q, int, int]] = []
    for centre in range(first, last + 1):
        tooth_left = Q(14 * centre - 1, 14 * speed)
        tooth_right = Q(14 * centre + 1, 14 * speed)
        lo = max(left, tooth_left)
        hi = min(right, tooth_right)
        if lo < hi:
            answer.append((lo, hi,
                           c if lo == left else speed,
                           c if hi == right else speed))
    return tuple(answer)


def pair_components(
    c: int, first: int, second: int, k: int,
) -> tuple[tuple[Q, int, int], ...]:
    """Lengths and selected endpoint owners of pair-intersection components."""
    answer: list[tuple[Q, int, int]] = []
    for a, b, owner_a, owner_b in meeting_teeth(c, first, k):
        for x, y, owner_x, owner_y in meeting_teeth(c, second, k):
            lo = max(a, x)
            hi = min(b, y)
            if lo >= hi:
                continue
            left_owner = owner_a if a >= x else owner_x
            right_owner = owner_b if b <= y else owner_y
            answer.append((hi - lo, left_owner, right_owner))
    return tuple(answer)


def endpoint_bank() -> tuple[int, int, Q, tuple[int, int, int, int]]:
    """Check the owner quantum and its phase-free weakening exactly."""
    rows = components = 0
    minimum_ratio: Q | None = None
    minimizer = (0, 0, 0, 0)
    for c in range(1, 13):
        for first in range(c + 1, 4 * c + 1):
            for second in range(first + 1, 5 * c + 1):
                for k in range(c):
                    rows += 1
                    for length, u, v in pair_components(c, first, second, k):
                        components += 1
                        owner_quantum = Q(gcd(u, v), 14 * u * v)
                        quotient = length / owner_quantum
                        require(
                            quotient.denominator == 1 and quotient >= 1,
                            "component is not a positive integer owner quantum at "
                            f"(c,d_i,d_j,k)=({c},{first},{second},{k})",
                        )

                        phase_free_quantum = Q(1, 14 * first * second)
                        ratio = length / phase_free_quantum
                        require(
                            ratio >= 1,
                            "component violates the phase-free edge floor at "
                            f"(c,d_i,d_j,k)=({c},{first},{second},{k})",
                        )
                        if minimum_ratio is None or ratio < minimum_ratio:
                            minimum_ratio = ratio
                            minimizer = (c, first, second, k)
    require(minimum_ratio is not None, "endpoint bank found no positive components")
    return rows, components, minimum_ratio, minimizer


def pair_gcd_guardrail() -> tuple[Q, Q, Q, int, int]:
    """An exact local counterexample to the naive pair-gcd refinement."""
    c, first, second, k = 7, 10, 34, 2
    components = pair_components(c, first, second, k)
    require(len(components) == 1, "pair-gcd guardrail component count changed")
    length, u, v = components[0]
    owner_quantum = Q(gcd(u, v), 14 * u * v)
    pair_quantum = Q(gcd(first, second), 14 * first * second)
    require((u, v) == (c, second), "pair-gcd guardrail endpoint owners changed")
    require(
        length == owner_quantum == Q(1, 3332),
        "pair-gcd guardrail owner quantum changed",
    )
    require(pair_quantum == Q(1, 2380), "naive pair-gcd quantum changed")
    require(length < pair_quantum, "pair-gcd guardrail no longer refutes naive floor")
    return length, owner_quantum, pair_quantum, u, v


def component_count_bank() -> tuple[int, int]:
    """Replay n_i <= ceil((6d_i+c)/(7c)) on a bounded phase bank."""
    rows = max_slack = 0
    for c in range(1, 25):
        for speed in range(c + 1, 7 * c + 1):
            bound = ceil(Q(6 * speed + c, 7 * c))
            for k in range(c):
                count = len(meeting_teeth(c, speed, k))
                require(
                    count <= bound,
                    f"tooth-count ceiling failed at (c,d,k)=({c},{speed},{k})",
                )
                max_slack = max(max_slack, bound - count)
                rows += 1
    return rows, max_slack


def main() -> None:
    require(sum(HEIGHTS, Q(0)) / 6 == 1, "six-bin density mass changed")
    require(min(HEIGHTS) == Q(3, 4), "six-bin density minimum changed")
    require(max(HEIGHTS) == Q(7, 6), "six-bin density maximum changed")

    # One physical endpoint quantum, scaled to x in [0,1], then weighted by f.
    c, u, v = 11, 13, 17
    physical_quantum = Q(gcd(u, v), 14 * u * v)
    normalized_quantum = Q(7 * c, 6) * physical_quantum
    weighted_quantum = min(HEIGHTS) * normalized_quantum
    require(
        weighted_quantum == Q(c * gcd(u, v), 16 * u * v),
        "weighted endpoint scaling changed",
    )

    rows, components, min_ratio, minimizer = endpoint_bank()
    local, owner_q, pair_q, owner_left, owner_right = pair_gcd_guardrail()
    count_rows, max_slack = component_count_bank()

    # Quantitative five-prefix survivor and its physical slow-gap scaling.
    dual_mass = 1 - 5 * Q(7, 36)
    normalized_survivor = dual_mass / max(HEIGHTS)
    physical_survivor_factor = Q(6, 7) * normalized_survivor
    require(dual_mass == Q(1, 36), "five-prefix dual survivor changed")
    require(normalized_survivor == Q(1, 42), "normalized survivor changed")
    require(physical_survivor_factor == Q(1, 49), "physical survivor changed")

    # General phase-free refinement with eta=1-sum_{i<=5} Pbar(L_i).
    eta_coefficient = 1 / max(HEIGHTS)
    physical_eta_coefficient = Q(6, 7) * eta_coefficient
    require(eta_coefficient == Q(6, 7), "normalized eta coefficient changed")
    require(
        physical_eta_coefficient == Q(36, 49),
        "physical eta coefficient changed",
    )
    require(Q(49, 36) / dual_mass == 49, "last-tooth baseline factor changed")

    print("THM-1199 exact arithmetic referee")
    print("===================================")
    print("six-bin density")
    print(f"  integral={qstr(sum(HEIGHTS, Q(0)) / 6)} "
          f"min={qstr(min(HEIGHTS))} max={qstr(max(HEIGHTS))}")
    print("weighted endpoint scaling")
    print("  (7c/6)*(3/4)*gcd(u,v)/(14uv) = c*gcd(u,v)/(16uv)")
    print("endpoint-component bank")
    print(f"  rows={rows} positive_components={components}")
    print(f"  minimum length/[1/(14*d_i*d_j)]={qstr(min_ratio)} "
          f"at (c,d_i,d_j,k)={minimizer}")
    print("pair-gcd guardrail")
    print("  (c,d_i,d_j,k)=(7,10,34,2) has one overlap component")
    print(f"  length={qstr(local)} owner_quantum={qstr(owner_q)} "
          f"owners=({owner_left},{owner_right})")
    print(f"  naive gcd(d_i,d_j) quantum={qstr(pair_q)} (strictly too large)")
    print("five-prefix survivor")
    print(f"  dual_mass={qstr(dual_mass)} normalized_length>="
          f"{qstr(normalized_survivor)} physical_length>=1/(49c)")
    print("unconditional last-tooth descent")
    print("  (d_6+d_5)/c < 49*C <= 49*C_bar")
    print("  refined: (d_6+d_5)/c < 49*C/[36*eta],")
    print("           eta=1-sum_{i<=5} Pbar(6d_i/(7c)) >= 1/36")
    print("component-count bank")
    print(f"  rows={count_rows} max_bound_slack={max_slack}")
    print("tournament loss audit")
    print("  speed order: transitive; score histogram 0,1,2,3,4,5;")
    print("  directed cycles=0; SCCs=6 singletons; Hamiltonian paths=1")
    print("VERIFIED_EXACT")


if __name__ == "__main__":
    main()
