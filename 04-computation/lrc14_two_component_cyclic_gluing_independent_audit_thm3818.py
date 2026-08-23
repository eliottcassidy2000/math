#!/usr/bin/env python3
"""Independent hostile audit of the THM-3818 cyclic-gluing corollary.

This checker does not import the primary companion.  It uses a sieve atlas, a
wall-cell construction of good arcs, direct shifted-grid sampling on a bounded
hostile universe, and a divisor-union count for the 11+2 residual.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from json import dumps
from math import gcd, isqrt, lcm
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


Q = 91**6
B = Q**2
GATES = 0


def gate(condition: bool, label: object) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"audit gate failed: {label}")


def least_residue_distance(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def smallest_prime_factors(limit: int) -> list[int]:
    spf = list(range(limit + 1))
    for prime in range(2, isqrt(limit) + 1):
        if spf[prime] != prime:
            continue
        for multiple in range(prime * prime, limit + 1, prime):
            if spf[multiple] == multiple:
                spf[multiple] = prime
    return spf


SPF = smallest_prime_factors(356)


def sieve_factorization(n: int) -> tuple[tuple[int, int], ...]:
    factors: list[tuple[int, int]] = []
    while n > 1:
        prime = SPF[n]
        exponent = 0
        while n % prime == 0:
            exponent += 1
            n //= prime
        factors.append((prime, exponent))
    return tuple(factors)


def trial_prime_divisors(n: int) -> tuple[int, ...]:
    """Prime support for values beyond the small ratio-atlas sieve."""

    answer: list[int] = []
    candidate = 2
    while candidate * candidate <= n:
        if n % candidate == 0:
            answer.append(candidate)
            while n % candidate == 0:
                n //= candidate
        candidate = 3 if candidate == 2 else candidate + 2
    if n > 1:
        answer.append(n)
    return tuple(answer)


def ratio_atlas_sieve() -> tuple[tuple[int, int], ...]:
    answer: list[tuple[int, int]] = []
    for total in range(3, 357):
        factors = sieve_factorization(total)
        if not factors or any(prime % 3 != 2 or exponent > 2 for prime, exponent in factors):
            continue
        answer.extend(
            (p, total - p)
            for p in range(1, (total + 1) // 2)
            if p < total - p and gcd(p, total - p) == 1
        )
    return tuple(answer)


def safe(value: Fraction, speeds: tuple[int, ...]) -> bool:
    return all(least_residue_distance(speed * value) >= Fraction(1, 14) for speed in speeds)


def longest_good_arc_from_walls(speeds: tuple[int, ...]) -> Fraction:
    """Independent wall-cell computation, not danger-interval merging."""

    walls = {Fraction(0), Fraction(1)}
    for speed in speeds:
        for k in range(speed):
            walls.add(Fraction(14 * k + 1, 14 * speed) % 1)
            walls.add(Fraction(14 * k - 1, 14 * speed) % 1)
    ordered = sorted(walls)
    gate(ordered[0] == 0 and ordered[-1] == 1, ("wall endpoints", speeds))

    best = Fraction(0)
    running = Fraction(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if safe(midpoint, speeds):
            running += right - left
            best = max(best, running)
        else:
            running = Fraction(0)
    gate(best > 0, ("positive wall-cell good arc", speeds))
    return best


def conductor(arc: Fraction) -> int:
    return (arc.denominator + arc.numerator - 1) // arc.numerator


def closed_form_bad_bound(scale: int, speed: int) -> int:
    common = gcd(scale, speed)
    order = scale // common
    distinct = order // 7 if order % 7 == 0 else (order + 6) // 7
    return common * distinct


def brute_bad_bound(scale: int, speed: int) -> int:
    """Max bad branches by sampling every shift cell exactly."""

    phase_steps = tuple(Fraction((speed * branch) % scale, scale) for branch in range(scale))
    walls = set()
    for phase in phase_steps:
        walls.add((-phase + Fraction(1, 14)) % 1)
        walls.add((-phase - Fraction(1, 14)) % 1)
    ordered = sorted(walls)
    probes = list(ordered)
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += 1
        probes.append(((left + right) / 2) % 1)
    return max(
        sum(least_residue_distance(shift + phase) < Fraction(1, 14) for phase in phase_steps)
        for shift in probes
    )


def divisors_at_least_two(n: int) -> set[int]:
    answer: set[int] = set()
    for divisor in range(1, isqrt(n) + 1):
        if n % divisor:
            continue
        if divisor >= 2:
            answer.add(divisor)
        partner = n // divisor
        if partner >= 2:
            answer.add(partner)
    return answer


def positive_congruence_count(bound: int, modulus: int, residue: int) -> int:
    representative = residue % modulus
    first = representative if representative else modulus
    return 0 if first > bound else 1 + (bound - first) // modulus


def audit_shifted_grids() -> dict[str, object]:
    brute_pairs = 0
    for scale in range(1, 29):
        for speed in range(1, 29):
            brute = brute_bad_bound(scale, speed)
            formula = closed_form_bad_bound(scale, speed)
            gate(brute == formula, ("brute/formula branch count", scale, speed, brute, formula))
            brute_pairs += 1

    unit_grids = 0
    worst_best_clearance = Fraction(1, 2)
    for scale in range(2, 65):
        for detuned in range(1, scale + 1):
            if gcd(scale, detuned) != 1:
                continue
            phase_steps = tuple(Fraction((detuned * k) % scale, scale) for k in range(scale))
            shifts = {(-phase) % 1 for phase in phase_steps}
            shifts.update(((-phase + Fraction(1, 2 * scale)) % 1) for phase in phase_steps)
            best_for_shift = tuple(max(least_residue_distance(shift + phase) for phase in phase_steps) for shift in shifts)
            minimum_best = min(best_for_shift)
            gate(minimum_best >= Fraction(scale - 1, 2 * scale), ("uniform-grid antipodal radius", scale, detuned, minimum_best))
            gate(minimum_best >= Fraction(1, 4), ("singleton grid is 1/14-safe", scale, detuned, minimum_best))
            worst_best_clearance = min(worst_best_clearance, minimum_best)
            unit_grids += 1

    y = Fraction(1, 100)
    core = tuple(range(8, 20))
    gate(all(safe(y, (speed,)) for speed in core), "nonprimitive hostile core")
    collapsed = tuple((2 * (y + branch) / 2) % 1 for branch in range(2))
    gate(collapsed == (y, y), ("nonprimitive collapsed orbit", collapsed))
    gate(all(not safe(phase, (1,)) for phase in collapsed), "nonprimitive singleton remains bad")

    ap_time = Fraction(1, 13)
    gate(safe(ap_time, tuple(range(1, 13))), "scale-one AP12 core")
    gate(not safe(ap_time, (13,)), "scale-one singleton hostile")
    gate(safe(Fraction(1, 14), tuple(range(1, 14))), "AP13 positive control")

    return {
        "brute_shift_count_pairs": brute_pairs,
        "unit_shifted_grids": unit_grids,
        "worst_unit_grid_best_clearance": str(worst_best_clearance),
        "hostiles": ["gcd(2,2)>1 collapses both lifts", "scale-one AP12 witness has speed13 at zero"],
    }


def audit_two_detuned(ratios: tuple[tuple[int, int], ...]) -> dict[str, object]:
    records: list[tuple[int, Fraction, int, int]] = []
    divisor_survivors: list[tuple[int, int, int]] = []
    for p, q in ratios:
        arc = longest_good_arc_from_walls((p, q))
        height = conductor(arc)
        records.append((height, arc, p, q))

        seams = {2} | divisors_at_least_two(p) | divisors_at_least_two(q)
        gate(all(scale < height for scale in seams), ("all divisor seams precede grid conductor", p, q, height, seams))
        for scale in seams:
            divisor_survivors.append((p, q, scale))
            gate(
                closed_form_bad_bound(scale, p) + closed_form_bad_bound(scale, q) >= scale,
                ("seam defeats separate danger tariff", p, q, scale),
            )

        for scale in range(2, height):
            if scale in seams:
                continue
            gate(
                closed_form_bad_bound(scale, p) + closed_form_bad_bound(scale, q) < scale,
                ("nonseam leaves a common branch", p, q, scale),
            )

    maximum = max(record[0] for record in records)
    extremizers = tuple(record for record in records if record[0] == maximum)
    gate(maximum == 415, ("independent conductor maximum", maximum))
    gate(extremizers == ((415, Fraction(6, 2485), 1, 355),), ("independent conductor extremizer", extremizers))
    gate(len(divisor_survivors) == 46837, ("independent divisor-union count", len(divisor_survivors)))
    gate(len(divisor_survivors) + len(ratios) == 52692, "independent divisor count including scale one")
    gate(max(scale for _, _, scale in divisor_survivors) == 355, "independent scale ceiling")

    incidence = {
        "scale_two": sum(scale == 2 for _, _, scale in divisor_survivors),
        "divides_p": sum(p % scale == 0 for p, _, scale in divisor_survivors),
        "divides_q": sum(q % scale == 0 for _, q, scale in divisor_survivors),
    }
    gate(incidence == {"scale_two": 5855, "divides_p": 19320, "divides_q": 25866}, ("seam incidence", incidence))
    return {
        "ratio_count": len(ratios),
        "grid_conductor_max": maximum,
        "grid_conductor_extremizer": [1, 355],
        "extremal_longest_arc": "6/2485",
        "divisor_union_survivors_s_ge_2": len(divisor_survivors),
        "residual_including_s_1": len(divisor_survivors) + len(ratios),
        "survivor_scale_max_s_ge_2": 355,
        "incidence_nonexclusive": incidence,
    }


def audit_singleton(ratio_set: set[tuple[int, int]]) -> dict[str, object]:
    gate(Q == 567869252041 and B == 322475487413604782665681, "Q/B constants")
    gate(Q % 13 == 0 and Q // 13 == 43682250157, "exact crossing cutoff")

    abstract_cases = 0
    abstract_max = 0
    for maximum in range(1, 81):
        for modulus in range(maximum + 1, 2 * maximum + 4):
            for residue in range(modulus):
                count = positive_congruence_count(13 * maximum, modulus, residue)
                gate(count <= 13, ("abstract packet candidate cap", maximum, modulus, residue, count))
                abstract_max = max(abstract_max, count)
                abstract_cases += 1
    gate(abstract_max == 13, ("abstract cap attained", abstract_max))
    gate(positive_congruence_count(1300, 101, 1) == 13, "named cap hostile")

    shape = tuple(4**power for power in range(12))
    edges = []
    for i in range(12):
        for j in range(i + 1, 12):
            common = gcd(shape[i], shape[j])
            reduced = tuple(sorted((shape[i] // common, shape[j] // common)))
            inert = all(prime % 3 == 2 for prime in trial_prime_divisors(common))
            if reduced in ratio_set and inert:
                edges.append((i, j))
    gate({vertex for edge in edges for vertex in edge} == set(range(12)), "canonical decoder support")
    packet_lcm = 1
    for i, j in edges:
        packet_lcm = lcm(packet_lcm, shape[i] + shape[j])
    gate(packet_lcm == 22906142720, ("independent packet lcm", packet_lcm))

    maximum = shape[-1]
    raw_singleton = 2 * Q * maximum + 1
    residue = raw_singleton % packet_lcm
    raw_fibre = (B - sum(shape) - raw_singleton) // packet_lcm + 1
    cutoff = 13 * maximum
    candidates = positive_congruence_count(cutoff, packet_lcm, residue)
    gate(maximum == 4194304 and cutoff == 54525952, "canonical maximum/cutoff")
    gate(residue == 9906946049, ("independent singleton residue", residue))
    gate(raw_fibre == 14077914720208, ("independent raw fibre", raw_fibre))
    gate(candidates == 0, ("canonical counterexample fibre empty", candidates))

    boundary_u = Q // 13
    for component_speed in (1, boundary_u // 2, boundary_u):
        for singleton_speed in (1, Q // 2, Q):
            common = gcd(component_speed, singleton_speed)
            height = max(component_speed // common, singleton_speed // common)
            gate(height <= Q, ("crossing boundary relation", component_speed, singleton_speed, height))
    above_u = boundary_u + 1
    above_t = 13 * above_u
    gate(above_t == Q + 13 and max(above_t, 1) > Q, "first-above crude-height hostile")

    return {
        "abstract_candidate_cases": abstract_cases,
        "abstract_candidate_max": abstract_max,
        "Q_over_13": Q // 13,
        "first_allowed_shape_maximum": Q // 13 + 1,
        "first_above_size_only_hostile_height": above_t,
        "canonical": {
            "U": maximum,
            "L": packet_lcm,
            "residue": residue,
            "raw_fibre": raw_fibre,
            "cutoff": cutoff,
            "remaining_candidates": candidates,
        },
    }


def main() -> None:
    ratios = ratio_atlas_sieve()
    gate(len(ratios) == 5855 and len(set(ratios)) == 5855, "independent ratio atlas")
    shifted = audit_shifted_grids()
    two_detuned = audit_two_detuned(ratios)
    singleton = audit_singleton(set(ratios))

    semantic = {"Q": Q, "B": B, "shifted": shifted, "singleton": singleton, "two_detuned": two_detuned}
    semantic_hash = sha256(dumps(semantic, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    print("THM-3818 CYCLIC-GLUING INDEPENDENT HOSTILE AUDIT")
    print("verdict=PASS;analytic_implications_sound;finite_constants_independently_rederived;LRC14_OPEN")
    print("method=sieve_ratio_atlas+wall_cell_good_arcs+direct_shift_cells+divisor_union+independent_congruence_count")
    print("shifted=" + dumps(shifted, sort_keys=True, separators=(",", ":")))
    print("two_detuned=" + dumps(two_detuned, sort_keys=True, separators=(",", ":")))
    print("singleton=" + dumps(singleton, sort_keys=True, separators=(",", ":")))
    print(f"gates={GATES}")
    print(f"semantic_sha256={semantic_hash}")
    print("scope=corollary_inside_THM3818_dimW11_two_component_branch;no_owner_arrival_or_global_closure")
    print("PASS")


if __name__ == "__main__":
    main()
