#!/usr/bin/env python3
"""Independent exact probe for the next THM-3818 two-component implication.

This file tests only consequences of the proved rank-eleven quotient

    n = s*u direct-sum t*v,  gcd(s,t)=1.

It checks the cyclic branch-gluing conductor, the exact two-detuned tariff on
all 5,855 decoder ratios, and the singleton packet-fibre collapse.  The
general proofs are in THM-3818.  Nothing here claims LRC(14).
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import product
from json import dumps
from math import gcd, lcm
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


Q = 91**6
B = Q**2
REQUIREMENTS = 0


def require(condition: bool, label: object) -> None:
    global REQUIREMENTS
    REQUIREMENTS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def factor(n: int) -> dict[int, int]:
    require(n >= 1, ("positive factor input", n))
    answer: dict[int, int] = {}
    candidate = 2
    while candidate * candidate <= n:
        while n % candidate == 0:
            answer[candidate] = answer.get(candidate, 0) + 1
            n //= candidate
        candidate = 3 if candidate == 2 else candidate + 2
    if n > 1:
        answer[n] = answer.get(n, 0) + 1
    return answer


def admissible_sum(total: int) -> bool:
    factors = factor(total)
    return bool(factors) and all(prime % 3 == 2 and exponent <= 2 for prime, exponent in factors.items())


def inert_scale(scale: int) -> bool:
    return all(prime % 3 == 2 for prime in factor(scale))


def ratio_atlas() -> tuple[tuple[int, int], ...]:
    ratios: list[tuple[int, int]] = []
    for total in range(3, 357):
        if not admissible_sum(total):
            continue
        for p in range(1, (total + 1) // 2):
            q = total - p
            if p < q and gcd(p, q) == 1:
                ratios.append((p, q))
    return tuple(ratios)


def reduced_pair(x: int, y: int) -> tuple[int, int]:
    common = gcd(x, y)
    return tuple(sorted((x // common, y // common)))


def decoder_edges(shape: tuple[int, ...], ratio_set: set[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    return tuple(
        (i, j)
        for i in range(len(shape))
        for j in range(i + 1, len(shape))
        if reduced_pair(shape[i], shape[j]) in ratio_set and inert_scale(gcd(shape[i], shape[j]))
    )


def branch_danger_bound(scale: int, speed: int) -> int:
    """Sharp maximum bad branches for one speed at radius 1/14."""

    require(scale >= 1 and speed >= 1, ("positive branch data", scale, speed))
    common = gcd(scale, speed)
    order = scale // common
    distinct_bad = order // 7 if order % 7 == 0 else order // 7 + 1
    return common * distinct_bad


def danger_bound_controls() -> dict[str, object]:
    singleton_min_safe = None
    for scale in range(2, 513):
        for detuned in range(1, 513):
            if gcd(scale, detuned) != 1:
                continue
            bad = branch_danger_bound(scale, detuned)
            require(bad < scale, ("primitive singleton leaves a branch", scale, detuned, bad))
            safe = scale - bad
            singleton_min_safe = safe if singleton_min_safe is None else min(singleton_min_safe, safe)

    y = Fraction(1, 100)
    core = tuple(range(8, 20))
    require(all(min((u * y) % 1, 1 - (u * y) % 1) >= Fraction(1, 14) for u in core), "nonprimitive core control")
    nonprimitive_phases = tuple((2 * (y + branch) / 2) % 1 for branch in range(2))
    require(nonprimitive_phases == (y, y), ("nonprimitive phase collapse", nonprimitive_phases))
    require(all(min(phase, 1 - phase) < Fraction(1, 14) for phase in nonprimitive_phases), "nonprimitive hostile")

    ap_time = Fraction(1, 13)
    require(all(min((u * ap_time) % 1, 1 - (u * ap_time) % 1) >= Fraction(1, 13) for u in range(1, 13)), "AP12 cited-time control")
    require((13 * ap_time) % 1 == 0, "scale-one branch hostile")
    full_ap_time = Fraction(1, 14)
    require(all(min((u * full_ap_time) % 1, 1 - (u * full_ap_time) % 1) >= Fraction(1, 14) for u in range(1, 14)), "AP13 positive control")

    return {
        "primitive_singleton_scales_checked": [2, 512],
        "primitive_detuned_values_checked": [1, 512],
        "minimum_certified_safe_branches": singleton_min_safe,
        "nonprimitive_hostile_phases": tuple(str(value) for value in nonprimitive_phases),
        "scale_one_method_hostile": "AP12@1/13 leaves speed13 at zero; AP13@1/14 is safe",
    }


def danger_intervals(speed: int) -> tuple[tuple[Fraction, Fraction], ...]:
    radius = Fraction(1, 14 * speed)
    intervals: list[tuple[Fraction, Fraction]] = [(Fraction(0), radius), (Fraction(1) - radius, Fraction(1))]
    for k in range(1, speed):
        centre = Fraction(k, speed)
        intervals.append((centre - radius, centre + radius))
    return tuple(intervals)


def longest_safe_arc(speeds: tuple[int, ...]) -> Fraction:
    """Longest component length of {x:min ||speed*x|| >= 1/14}."""

    intervals = sorted(interval for speed in speeds for interval in danger_intervals(speed))
    merged: list[tuple[Fraction, Fraction]] = []
    for left, right in intervals:
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    gaps = tuple(merged[index + 1][0] - merged[index][1] for index in range(len(merged) - 1))
    require(bool(gaps), ("nonempty safe-gap ledger", speeds))
    answer = max(gaps)
    require(answer > 0, ("positive safe arc", speeds))
    return answer


def ceiling_reciprocal(value: Fraction) -> int:
    return (value.denominator + value.numerator - 1) // value.numerator


def two_detuned_atlas(ratios: tuple[tuple[int, int], ...]) -> dict[str, object]:
    """Exact 11+2 gluing/tariff census over every THM-3818 ratio."""

    conductor_records: list[tuple[int, Fraction, int, int]] = []
    tariff_survivors: list[tuple[int, int, int]] = []
    for p, q in ratios:
        safe_arc = longest_safe_arc((p, q))
        conductor = ceiling_reciprocal(safe_arc)
        conductor_records.append((conductor, safe_arc, p, q))
        for scale in range(2, conductor):
            bad_bound = branch_danger_bound(scale, p) + branch_danger_bound(scale, q)
            if bad_bound >= scale:
                tariff_survivors.append((p, q, scale))
                require(
                    scale == 2 or p % scale == 0 or q % scale == 0,
                    ("two-detuned survivor must be the scale-two/divisor seam", p, q, scale, bad_bound),
                )
            else:
                require(
                    scale != 2 and p % scale != 0 and q % scale != 0,
                    ("classification converse", p, q, scale, bad_bound),
                )

    max_conductor = max(record[0] for record in conductor_records)
    max_conductor_records = tuple(record for record in conductor_records if record[0] == max_conductor)
    require(max_conductor == 415, ("two-speed conductor maximum", max_conductor))
    require(max_conductor_records == ((415, Fraction(6, 2485), 1, 355),), ("unique conductor extremizer", max_conductor_records))
    require(len(tariff_survivors) == 46837, ("two-detuned tariff survivor census", len(tariff_survivors)))
    require(len(tariff_survivors) + len(ratios) == 52692, "two-detuned census including scale one")
    require(max(scale for _, _, scale in tariff_survivors) == 355, "two-detuned survivor scale ceiling")

    seam_types = {
        "scale_two": sum(1 for p, q, scale in tariff_survivors if scale == 2),
        "divides_p": sum(1 for p, q, scale in tariff_survivors if p % scale == 0),
        "divides_q": sum(1 for p, q, scale in tariff_survivors if q % scale == 0),
    }
    return {
        "ratio_count": len(ratios),
        "maximum_exact_grid_conductor": max_conductor,
        "unique_grid_conductor_extremizer": {"ratio": [1, 355], "longest_safe_arc": "6/2485"},
        "tariff_survivor_triples_s_ge_2": len(tariff_survivors),
        "tariff_residual_triples_including_s_1": len(tariff_survivors) + len(ratios),
        "tariff_survivor_scale_max_s_ge_2": 355,
        "tariff_survivor_exact_form_s_ge_2": "scale=2 or scale|p or scale|q",
        "seam_incidence_counts_nonexclusive": seam_types,
    }


def positive_congruence_count(bound: int, modulus: int, residue: int) -> int:
    require(bound >= 0 and modulus >= 1, ("congruence-count input", bound, modulus))
    reduced = residue % modulus
    first = reduced if reduced > 0 else modulus
    return 0 if first > bound else 1 + (bound - first) // modulus


def singleton_fibre(ratio_set: set[tuple[int, int]]) -> dict[str, object]:
    """Collapse the canonical 12+1 hostile after the gluing obstruction."""

    shape = tuple(4**power for power in range(12))
    edges = decoder_edges(shape, ratio_set)
    require(bool(edges), "powers-of-four decoder component")
    incident = {vertex for edge in edges for vertex in edge}
    require(len(incident) == 12, ("powers-of-four decoder connected support", incident))
    packet_lcm = 1
    for i, j in edges:
        packet_lcm = lcm(packet_lcm, shape[i] + shape[j])
    require(packet_lcm == 22906142720, ("canonical packet lcm", packet_lcm))

    maximum = max(shape)
    require(packet_lcm > maximum, "connected packet lcm exceeds component maximum")
    counterexample_cutoff = 13 * maximum
    singleton_0 = 2 * Q * maximum + 1
    residue = singleton_0 % packet_lcm
    candidate_count = positive_congruence_count(counterexample_cutoff, packet_lcm, residue)
    require(candidate_count == 0, ("canonical hostile counterexample fibre must be empty", candidate_count))

    dominance_tail_count = (B - sum(shape) - singleton_0) // packet_lcm + 1
    require(dominance_tail_count == 14077914720208, ("inherited raw packet fibre", dominance_tail_count))
    require(13 * maximum <= Q, "canonical component is below the crossing-row cutoff")

    sharp_control = positive_congruence_count(13 * 100, 101, 1)
    require(sharp_control == 13, ("abstract thirteen-candidate sharp control", sharp_control))

    return {
        "component_shape": "(1,4,...,4^11)",
        "component_maximum": maximum,
        "packet_lcm": packet_lcm,
        "singleton_residue": residue,
        "raw_dominance_packet_fibre": dominance_tail_count,
        "counterexample_cutoff": counterexample_cutoff,
        "counterexample_candidates_in_same_packet": candidate_count,
        "general_candidate_bound": 13,
        "crossing_row_shape_floor": Q // 13 + 1,
        "abstract_bound_sharp_control": {"U": 100, "L": 101, "residue": 1, "candidates": sharp_control},
    }


def main() -> None:
    ratios = ratio_atlas()
    require(len(ratios) == 5855, ("THM-3818 ratio count", len(ratios)))
    require(len(set(ratios)) == len(ratios), "ratio uniqueness")

    branch = danger_bound_controls()
    pair_atlas = two_detuned_atlas(ratios)
    singleton = singleton_fibre(set(ratios))

    semantic = {"Q": Q, "B": B, "branch": branch, "pair_atlas": pair_atlas, "singleton": singleton}
    semantic_hash = sha256(dumps(semantic, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    print("THM-3818 NEXT CYCLIC-GLUING / SINGLETON-FIBRE PROBE")
    print("status=PROVED_ANALYTIC_IMPLICATIONS+FINITE_EXACT_CONTROLS;LRC14_OPEN")
    print(f"Q={Q}")
    print(f"B={B}")
    print("cyclic_gluing=if_component_good_arc_length_lambda_and_pack_scale_s_satisfy_s*lambda>=1_then_row_is_lonely")
    print("universal_scale_gate=for_small_component_size_b_and_maximum_V,_s<7*(b+1)*V/(13-b)_in_any_counterexample")
    print("singleton_consequence=12+1_counterexample_forces_nonsingleton_scale_s=1")
    print("singleton_killer_gate=singleton_speed_t<=13*U")
    print("singleton_crossing_gate=U>Q/13=43682250157")
    print("singleton_packet_gate=L>U_and_t_mod_L_known_imply_at_most_13_candidates")
    print("two_detuned_consequence=s=1_or_(s>=2_seam);s>=3_forces_s_divides_p_or_s_divides_q;s>=2_seams=46837;including_s1=52692;positive_scale_max=355")
    print("branch_controls=" + dumps(branch, sort_keys=True, separators=(",", ":")))
    print("two_detuned_atlas=" + dumps(pair_atlas, sort_keys=True, separators=(",", ":")))
    print("singleton=" + dumps(singleton, sort_keys=True, separators=(",", ":")))
    print(f"requirements={REQUIREMENTS}")
    print(f"semantic_sha256={semantic_hash}")
    print("scope=rank11_two_component_quotient_only;no_owner_or_global_LRC14_closure")
    print("PASS")


if __name__ == "__main__":
    main()
