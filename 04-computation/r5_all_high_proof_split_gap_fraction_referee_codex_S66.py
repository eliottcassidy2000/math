#!/usr/bin/env python3
"""Exact Fraction referee for an all-high gap in THM-1101's proof split.

This is not a Lonely Runner counterexample.  It proves only that the two
branches stated in THM-1101 do not cover every five-killer tuple: the finite
horn was restricted to killers below 235, while the stated measure horn needs
the fifth killer to exceed the threshold computed after four removals.

All set construction and every decision use fractions.Fraction.  Decimal
values are presentation-only.  The interval carrier retains exact endpoint
coordinates and tooth-owner labels.

Tournament Analysis audit.  Vertices are challenged rather than assumed to
be runners: runners, danger combs, core components, wall-crossing endpoints,
residues, and proof obligations are all considered.  The useful pairwise
observable on endpoints is exact left/right order; the coordinate is the
gauge that makes it a binary relation, exact equality is the tie, and the tie
Hamiltonian path is the sorted endpoint word after coincident endpoints are
coalesced.  This tournament is transitive.  Order alone destroys the metric
lengths, endpoint owners, removal stages, and residue-witness clearance, so
the faithful carrier is the coordinate-plus-owner interval word, not a plain
runner or endpoint tournament.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, getcontext
from fractions import Fraction as F


getcontext().prec = 24

CORE = (1, 2, 4, 5, 7, 9, 11, 12)
REMOVED = (294, 298, 299, 303)
FIFTH = 320
FINITE_BOUND = 235


@dataclass(frozen=True)
class Boundary:
    value: F
    owner: str


Interval = tuple[Boundary, Boundary]
Region = list[Interval]


def decimal(value: F) -> str:
    return str(Decimal(value.numerator) / Decimal(value.denominator))


def remove_danger(region: Region, speed: int) -> Region:
    """Subtract D_speed using exact tooth endpoints and retain their owners."""

    result: Region = []
    denominator = 14 * speed
    for left, right in region:
        cursor = left
        j_lo = int(left.value * speed) - 1
        j_hi = int(right.value * speed) + 1
        for j in range(j_lo, j_hi + 1):
            tooth_left = F(14 * j - 1, denominator)
            tooth_right = F(14 * j + 1, denominator)
            if tooth_right <= left.value or right.value <= tooth_left:
                continue

            clipped_left = max(left.value, tooth_left)
            clipped_right = min(right.value, tooth_right)
            if cursor.value < clipped_left:
                result.append(
                    (
                        cursor,
                        Boundary(clipped_left, f"D_{speed}:tooth_{j}:left"),
                    )
                )
            if cursor.value < clipped_right:
                cursor = Boundary(
                    clipped_right, f"D_{speed}:tooth_{j}:right"
                )

        if cursor.value < right.value:
            result.append((cursor, right))
    return result


def safe_region(speeds: tuple[int, ...]) -> Region:
    region = [
        (Boundary(F(0), "unit:0"), Boundary(F(1), "unit:1"))
    ]
    for speed in speeds:
        region = remove_danger(region, speed)
    return region


def length(interval: Interval) -> F:
    return interval[1].value - interval[0].value


def least_absolute_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def main() -> None:
    core_safe = safe_region(CORE)
    remainder = core_safe
    for killer in REMOVED:
        remainder = remove_danger(remainder, killer)

    component_count = len(remainder)
    measure = sum((length(interval) for interval in remainder), F(0))
    longest = max(length(interval) for interval in remainder)
    longest_intervals = [
        interval for interval in remainder if length(interval) == longest
    ]

    threshold_mass = F(component_count, 1) / (6 * measure)
    threshold_component = F(1, 1) / (3 * longest)
    threshold = min(threshold_mass, threshold_component)
    ratio = threshold / REMOVED[-1]

    assert component_count == 198
    assert measure == F(14_258_767_904, 152_794_649_007)
    assert longest == F(431, 415_716)
    assert threshold_mass == F(5_042_223_417_231, 14_258_767_904)
    assert threshold_component == F(138_572, 431)
    assert threshold == threshold_component
    assert FIFTH <= threshold
    assert all(killer >= FINITE_BOUND for killer in REMOVED + (FIFTH,))

    family = CORE + REMOVED + (FIFTH,)
    covering_owner: dict[int, int] = {}
    for modulus in range(2, 15):
        owners = [speed for speed in family if speed % modulus == 0]
        assert owners
        covering_owner[modulus] = owners[0]

    witness_modulus = 22
    witness_numerator = 7
    witness_threshold = (witness_modulus + 13) // 14
    witness_clearances = tuple(
        least_absolute_residue(speed * witness_numerator, witness_modulus)
        for speed in family
    )
    assert min(witness_clearances) >= witness_threshold

    after_fifth = remove_danger(remainder, FIFTH)
    witness_time = F(witness_numerator, witness_modulus)
    witness_components = [
        interval
        for interval in after_fifth
        if interval[0].value <= witness_time <= interval[1].value
    ]
    assert len(witness_components) == 1

    endpoint_values = [
        boundary.value for interval in remainder for boundary in interval
    ]
    endpoint_owners = [
        boundary.owner for interval in remainder for boundary in interval
    ]
    endpoint_count = 2 * component_count
    assert len(endpoint_values) == endpoint_count
    assert len(set(endpoint_values)) == endpoint_count
    assert endpoint_values == sorted(endpoint_values)
    owner_speed_histogram = Counter(
        owner.split(":", 1)[0] for owner in endpoint_owners
    )

    print("r=5 all-high proof-split gap: exact Fraction referee")
    print("arithmetic=fractions.Fraction; decimal fields are display-only")
    print(f"core={list(CORE)}")
    print(f"removed_four={REMOVED}")
    print(f"fifth={FIFTH}")
    print(f"finite_horn_bound={FINITE_BOUND} (script scope: every killer < bound)")
    print(f"finite_horn_in_scope={all(k < FINITE_BOUND for k in REMOVED + (FIFTH,))}")
    print(f"components_N={component_count}")
    print(f"measure_mu={measure} decimal={decimal(measure)}")
    print(f"longest_L={longest} decimal={decimal(longest)}")
    print(
        f"Tmass=N/(6mu)={threshold_mass} "
        f"decimal={decimal(threshold_mass)}"
    )
    print(
        f"Tcomp=1/(3L)={threshold_component} "
        f"decimal={decimal(threshold_component)}"
    )
    print(f"T=min(Tmass,Tcomp)={threshold} decimal={decimal(threshold)}")
    print(f"R=T/k4={ratio} decimal={decimal(ratio)}")
    print(f"measure_horn_condition_k5_gt_T={FIFTH > threshold}")
    print(f"exact_gap_inequality=k5={FIFTH} <= T={threshold}")
    for index, interval in enumerate(longest_intervals, start=1):
        print(
            f"longest_interval_{index}="
            f"[{interval[0].value},{interval[1].value}] "
            f"owners=({interval[0].owner},{interval[1].owner})"
        )

    covering_text = ",".join(
        f"{modulus}:{covering_owner[modulus]}" for modulus in range(2, 15)
    )
    print("covering_q2_through_q14=True")
    print(f"covering_owner_q_to_speed={covering_text}")
    print(
        f"explicit_safe_witness=(q,a)=({witness_modulus},{witness_numerator}) "
        f"time={witness_time}"
    )
    print(f"witness_required_least_residue={witness_threshold}")
    print(f"witness_least_residues={witness_clearances}")
    print(f"witness_min_clearance={min(witness_clearances)}")
    witness_component = witness_components[0]
    print(
        "witness_component_after_fifth="
        f"[{witness_component[0].value},{witness_component[1].value}]"
    )
    print("classification=proof-gap, not an LRC counterexample")
    print(
        "reason=finite horn is out of scope and its stated measure sufficient "
        "condition fails, while the explicit q=22 witness proves safety"
    )

    print("tournament_analysis:")
    print(
        "  alternate_vertices=runners|danger_combs|core_components|"
        "wall_crossing_endpoints|residues|proof_obligations"
    )
    print("  chosen_vertices=surviving exact endpoints")
    print("  pairwise_observable=exact rational left/right comparison")
    print("  switch_gauge=orient x->y iff coordinate(x)<coordinate(y)")
    print(
        "  tie_hamiltonian_path=coalesce exact equalities, then use the sorted "
        "endpoint word"
    )
    print(f"  endpoint_vertices={endpoint_count}")
    print(f"  score_histogram={{0..{endpoint_count - 1}: multiplicity 1}}")
    print("  directed_cycles=0")
    print(f"  strongly_connected_components={endpoint_count} singleton")
    print("  hamiltonian_path_count=1")
    print("  edge_flip_locus=exact endpoint coincidence")
    print(f"  endpoint_owner_speed_histogram={dict(sorted(owner_speed_histogram.items()))}")
    print(
        "  preserved_by_full_carrier=component count, mass, longest length, "
        "removal owner, witness membership"
    )
    print(
        "  destroyed_by_order_only=metric lengths, endpoint owners, removal "
        "stage, residue-witness clearance"
    )
    print(
        "  challenged_assumption=tournament vertices need not be runners; "
        "endpoint order is faithful only with coordinate and owner sidecars"
    )
    print("certificate=PASS")


if __name__ == "__main__":
    main()
