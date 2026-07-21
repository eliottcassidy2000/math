#!/usr/bin/env python3
"""Exact audit for THM-2043: period-14 parity-Hasse packets.

The characteristic-seven decomposition

    F_7[C_14] = F_7[X]/(X-1)^7 x F_7[X]/(X+1)^7

identifies a reduced phase function with fourteen Hasse jets.  The audit checks
that this coordinate map is invertible, tests the jets on three exact phase
functions (danger count, weak-safe mask, and boundary-contact count), and
exhibits an infinite family on which every period-14 phase statistic is blind
to the distinction between the tight AP and a strict denominator witness.

The constructive repair is a bounded owner-labelled CRT chart: residues modulo
13 and 14 recover every positive speed at most 181.  This is a finite-bank
tomography statement, not an LRC proof.

Tournament Analysis uses proof carriers, not runners.  The observable records
retention of strict-witness status, endpoint labels, magnitude, the full
phase vector, and the bad-prime filtration.  Lexicographic comparison orients
the tournament and declaration order breaks ties.  This preserves the proof
obligations and destroys raw circle geometry.  The challenged assumption is
that deeper jets can recover information absent from the underlying mod-14
phase function: completeness shows they cannot.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from math import comb, gcd
PERIOD = 14
JET_PRIME = 7


def phase_functions(
    speeds: tuple[int, ...], period: int = PERIOD
) -> tuple[tuple[int, ...], ...]:
    danger: list[int] = []
    weak_safe: list[int] = []
    boundary: list[int] = []
    for a in range(period):
        distances = []
        for speed in speeds:
            residue = (speed * a) % period
            distances.append(min(residue, period - residue))
        danger.append(sum(distance == 0 for distance in distances))
        weak_safe.append(int(all(distance >= 1 for distance in distances)))
        boundary.append(sum(distance == 1 for distance in distances))
    return tuple(danger), tuple(weak_safe), tuple(boundary)


def hasse_jets(values: tuple[int, ...], depth: int = 6) -> tuple[int, ...]:
    assert len(values) == PERIOD and 0 <= depth <= 6
    out = []
    for center in (1, -1):
        for order in range(depth + 1):
            value = sum(
                comb(exponent, order)
                * coefficient
                * center ** (exponent - order)
                for exponent, coefficient in enumerate(values)
                if exponent >= order
            )
            out.append(value % JET_PRIME)
    return tuple(out)


def jet_packet(speeds: tuple[int, ...], depth: int = 6) -> tuple[int, ...]:
    return sum((hasse_jets(values, depth) for values in phase_functions(speeds)), ())


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    work = [[entry % prime for entry in row] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if work[r][col]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inv = pow(work[rank][col], -1, prime)
        work[rank] = [(inv * entry) % prime for entry in work[rank]]
        for r in range(rows):
            if r == rank or not work[r][col]:
                continue
            factor = work[r][col]
            work[r] = [
                (left - factor * right) % prime
                for left, right in zip(work[r], work[rank])
            ]
        rank += 1
        if rank == rows:
            break
    return rank


def polynomial_product(left: list[int], right: list[int], prime: int) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] = (out[i + j] + a * b) % prime
    return out


def hasse_algebra_checks() -> tuple[int, int]:
    matrix: list[list[int]] = []
    for center in (1, -1):
        for order in range(7):
            matrix.append(
                [
                    (
                        comb(exponent, order)
                        * center ** (exponent - order)
                        if exponent >= order
                        else 0
                    )
                    % JET_PRIME
                    for exponent in range(PERIOD)
                ]
            )
    rank = rank_mod(matrix, JET_PRIME)
    assert rank == PERIOD

    minus = [(-1) ** (7 - j) * comb(7, j) % 7 for j in range(8)]
    plus = [comb(7, j) % 7 for j in range(8)]
    product = polynomial_product(minus, plus, 7)
    target = [0] * 15
    target[0] = -1 % 7
    target[14] = 1
    assert product == target
    return rank, len(matrix)


def mobius(n: int) -> int:
    if n == 1:
        return 1
    factors = 0
    d = 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            factors += 1
            if n % d == 0:
                return 0
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors += 1
    return -1 if factors % 2 else 1


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def ramanujan_14(n: int) -> int:
    return sum(d * mobius(14 // d) for d in divisors(gcd(14, n)))


def cyclic_convolution(left: tuple[int, ...], right: tuple[int, ...], prime: int) -> tuple[int, ...]:
    out = [0] * PERIOD
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[(i + j) % PERIOD] = (out[(i + j) % PERIOD] + a * b) % prime
    return tuple(out)


def exact_period_component(values: tuple[int, ...], prime: int) -> tuple[int, ...]:
    inverse_period = pow(PERIOD, -1, prime)
    projector = tuple(inverse_period * ramanujan_14(x) % prime for x in range(PERIOD))
    return cyclic_convolution(projector, tuple(value % prime for value in values), prime)


def energy_mod(values: tuple[int, ...], prime: int) -> int:
    return sum(
        values[a] * values[b] * ramanujan_14(a - b)
        for a in range(PERIOD)
        for b in range(PERIOD)
    ) % prime


def char3_component_checks() -> dict[str, object]:
    ap = tuple(range(1, 14))
    covering = tuple(list(range(1, 12)) + [13, 84])
    ap_danger, ap_weak, _ = phase_functions(ap)
    covering_danger, _, _ = phase_functions(covering)
    assert not any(exact_period_component(ap_danger, 3))
    assert not any(exact_period_component(covering_danger, 3))
    ap_weak_component = exact_period_component(ap_weak, 3)
    assert any(ap_weak_component)
    assert energy_mod(ap_weak, 3) == 0
    return {
        "danger_collision_components_mod3": "both_zero",
        "ap_weak_energy_mod3": 0,
        "ap_weak_component_mod3": ap_weak_component,
        "ap_weak_component_nonzero": True,
    }


def q_threshold(speeds: tuple[int, ...], ceiling: int = PERIOD) -> int | None:
    """First denominator q for which a=1 has no runner at zero."""
    return next(
        (q for q in range(2, ceiling + 1) if all(speed % q for speed in speeds)),
        None,
    )


def first_strict_certificate(
    speeds: tuple[int, ...], ceiling: int = 42
) -> tuple[int, int, int] | None:
    """Return (q,a,minimum integer slack) for the first strict rational phase."""
    for q in range(2, ceiling + 1):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            slacks = []
            for speed in speeds:
                residue = (a * speed) % q
                slacks.append(14 * min(residue, q - residue) - q)
            if min(slacks) > 0:
                return q, a, min(slacks)
    return None


def strict_margin(speeds: tuple[int, ...], q: int, a: int) -> int:
    return min(
        14 * min((a * speed) % q, q - ((a * speed) % q)) - q
        for speed in speeds
    )


def blocked_mask(speeds: tuple[int, ...], ceiling: int = 13) -> tuple[int, ...]:
    return tuple(int(any(speed % q == 0 for speed in speeds)) for q in range(2, ceiling + 1))


def far_alias(index: int) -> tuple[int, ...]:
    """The q-threshold-blind infinite family with a fixed q=41 exit."""
    return replaced(12, 96 + 3444 * index)


def infinite_far_alias_checks(test_depth: int = 40) -> dict[str, object]:
    ap = tuple(range(1, 14))
    ap_mask = blocked_mask(ap)
    for index in range(test_depth + 1):
        speeds = far_alias(index)
        assert phase_functions(speeds) == phase_functions(ap)
        assert jet_packet(speeds) == jet_packet(ap)
        assert blocked_mask(speeds) == ap_mask
        assert q_threshold(speeds) == 14
        assert strict_margin(speeds, 41, 17) == 1

    replacement = far_alias(1)[-1]
    height = (replacement - 12) // 14
    assert replacement == 3540 and height == 252 and height % 7 == 0
    return {
        "family": "T_n=AP with 12 replaced by 96+3444n",
        "symbolic_periods": "3444=lcm(12,14,41)",
        "q2_to_q13_blocked_mask_equal_AP": True,
        "q_threshold": 14,
        "uniform_strict_certificate": "(q,a,integer_margin)=(41,17,1)",
        "checked_indices": f"0..{test_depth}",
        "mod7_height_collision": "n=1 has lift height 252=0 mod7",
    }


def finite_seven_adic_height_no_go(max_precision: int = 6) -> dict[str, object]:
    """Construct an alias invisible to every requested finite height precision."""
    ap = tuple(range(1, 14))
    rows = []
    for precision in range(max_precision + 1):
        power = 7**precision
        inverse = pow(power, -1, 41)
        replacement = 12 + 84 * power * inverse
        speeds = replaced(12, replacement)
        height = (replacement - 12) // 14
        assert phase_functions(speeds) == phase_functions(ap)
        assert blocked_mask(speeds) == blocked_mask(ap)
        assert q_threshold(speeds) == 14
        assert height % power == 0
        assert replacement % 41 == 14
        assert strict_margin(speeds, 41, 17) == 1
        rows.append((precision, inverse, replacement))
    return {
        "precision_rows_k_inverse_replacement": tuple(rows),
        "all_height_mod_7_to_k_equal_AP": True,
        "all_uniform_q41_margin": 1,
    }


def strict_primitive_numerators(speeds: tuple[int, ...], q: int) -> tuple[int, ...]:
    return tuple(
        a
        for a in range(1, q)
        if gcd(a, q) == 1 and strict_margin(speeds, q, a) > 0
    )


def alias_audit() -> dict[str, object]:
    """Audit one full 14-step lift window above each AP speed."""
    ap = tuple(range(1, 14))
    histogram: Counter[int] = Counter()
    strict_histogram: Counter[int] = Counter()
    aliases = []
    for drop in range(1, 14):
        for height in range(1, 13):
            speeds = replaced(drop, drop + 14 * height)
            assert phase_functions(speeds) == phase_functions(ap)
            assert jet_packet(speeds) == jet_packet(ap)
            threshold = q_threshold(speeds)
            assert threshold is not None
            histogram[threshold] += 1
            strict = first_strict_certificate(speeds)
            assert strict is not None
            strict_histogram[strict[0]] += 1
            aliases.append((drop, height, speeds, threshold, strict))

    direct = [row for row in aliases if row[3] <= 13]
    assert len(aliases) == 156
    assert len(direct) == 63
    return {
        "aliases": len(aliases),
        "q_threshold_histogram": dict(sorted(histogram.items())),
        "direct_q_le_13": len(direct),
        "q14_only_in_window": len(aliases) - len(direct),
        "strict_certificate_q_le_42": len(aliases),
        "first_strict_q_histogram": dict(sorted(strict_histogram.items())),
    }


def owner_crt_checks(bound: int = 181) -> dict[str, object]:
    """Owner-labelled residues modulo 13 and 14 recover bounded heights."""
    chart = {(speed % 13, speed % 14): speed for speed in range(1, bound + 1)}
    assert len(chart) == bound
    assert all(chart[(speed % 13, speed % 14)] == speed for speed in range(1, bound + 1))
    return {
        "height_bound": bound,
        "owner_residue_pairs": len(chart),
        "collisions": bound - len(chart),
        "modulus_product": 13 * 14,
    }


CARRIERS = {
    # strict-witness status, endpoint labels, magnitude, full phase, Hasse depth
    "labelled_packet_sheaf": (1, 1, 1, 1, 7),
    "owner_labelled_13x14_chart": (0, 1, 1, 1, 7),
    "full_jets_plus_height_owner": (0, 1, 1, 1, 7),
    "full_unlabelled_hasse_packet": (0, 0, 0, 1, 7),
    "depth3_unlabelled_hasse_packet": (0, 0, 0, 0, 4),
    "depth1_unlabelled_hasse_packet": (0, 0, 0, 0, 2),
    "parity_values_only": (0, 0, 0, 0, 1),
    "scalar_ramanujan_energy": (0, 0, 0, 0, 0),
}


def tournament_fingerprint() -> dict[str, object]:
    names = list(CARRIERS)
    wins = {name: set() for name in names}
    for left, right in combinations(names, 2):
        winner, loser = (
            (left, right) if CARRIERS[left] >= CARRIERS[right] else (right, left)
        )
        wins[winner].add(loser)
    score_hist = dict(sorted(Counter(map(len, wins.values())).items()))
    cycles = 0
    for a, b, c in combinations(names, 3):
        cycles += int(
            (b in wins[a] and c in wins[b] and a in wins[c])
            or (c in wins[a] and b in wins[c] and a in wins[b])
        )
    order = sorted(
        names,
        key=lambda name: (CARRIERS[name], -names.index(name)),
        reverse=True,
    )
    assert all(order[j + 1] in wins[order[j]] for j in range(len(order) - 1))
    return {
        "score_hist": score_hist,
        "directed_3cycles": cycles,
        "scc_sizes": [1] * len(names),
        "hamiltonian_path_count": 1,
        "tie_hamiltonian_path": " > ".join(order),
    }


def replaced(drop: int, add: int) -> tuple[int, ...]:
    return tuple(sorted((set(range(1, 14)) - {drop}) | {add}))


def replace_many(drops: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(range(1, 14)) - set(drops)) | set(adds)))


def named_route_fiber_audit() -> dict[str, object]:
    """Cheap exact audit of the eleven named HYP-2979 controls."""
    ap = tuple(range(1, 14))
    rows = (
        ("AP", "BOUNDARY-AP-GW", ap),
        ("GW12to24", "BOUNDARY-AP-GW", replaced(12, 24)),
        ("residue_liar12to26", "Q-WITNESS", replaced(12, 26)),
        ("near12to36", "K33-STATE-LIFT", replaced(12, 36)),
        ("petal10to20", "BOUNDARY-PETAL", replaced(10, 20)),
        ("petal13to26", "BOUNDARY-PETAL", replaced(13, 26)),
        ("P10plusGW", "BOUNDARY-PETAL", replace_many((10, 12), (20, 24))),
        ("P10plusK33", "K33-STATE-LIFT", replace_many((10, 12), (20, 36))),
        ("covering12to84", "COVERING-MOMENT", replaced(12, 84)),
        ("covering12to168", "COVERING-MOMENT", replaced(12, 168)),
        ("covering6to98", "COVERING-MOMENT", replaced(6, 98)),
    )

    def fibers(signature):
        grouped: dict[tuple, list[tuple[str, str]]] = {}
        for name, route, speeds in rows:
            grouped.setdefault(signature(speeds), []).append((name, route))
        mixed = [fiber for fiber in grouped.values() if len({route for _, route in fiber}) > 1]
        return grouped, mixed

    exact_fibers, exact_mixed = fibers(phase_functions)
    assert len(exact_fibers) == 3 and len(exact_mixed) == 1
    ap_fiber = next(fiber for fiber in exact_fibers.values() if ("AP", "BOUNDARY-AP-GW") in fiber)
    assert len(ap_fiber) == 7
    depth_census = []
    for depth in range(7):
        depth_fibers, depth_mixed = fibers(lambda speeds, d=depth: jet_packet(speeds, d))
        depth_census.append((depth, len(depth_fibers), len(depth_mixed), sum(map(len, depth_mixed))))
    assert all(row[1:] == (3, 1, 7) for row in depth_census)
    return {
        "named_rows": len(rows),
        "exact_phase_fibers": len(exact_fibers),
        "mixed_exact_phase_fibers": len(exact_mixed),
        "ap_fiber_size": len(ap_fiber),
        "ap_fiber_route_counts": dict(sorted(Counter(route for _, route in ap_fiber).items())),
        "ap_fiber_names": tuple(name for name, _ in ap_fiber),
        "jet_depth_census": tuple(depth_census),
    }


def main() -> None:
    rank, coordinate_count = hasse_algebra_checks()
    print("THM-2043 PERIOD-14 PARITY-HASSE JET AUDIT")
    print(f"hasse_coordinate_count={coordinate_count}")
    print(f"hasse_matrix_rank_mod7={rank}")
    print("factorization_mod7=(X-1)^7*(X+1)^7=X^14-1 PASS")

    ap = tuple(range(1, 14))
    liar = replaced(12, 26)
    assert phase_functions(ap) == phase_functions(liar)
    assert jet_packet(ap, 6) == jet_packet(liar, 6)
    assert q_threshold(ap) == 14
    assert q_threshold(liar) == 12
    assert first_strict_certificate(ap) is None
    assert first_strict_certificate(liar) == (12, 1, 2)
    assert phase_functions(ap, 12) != phase_functions(liar, 12)
    print("SHARP MAGNITUDE-BLIND COLLISION")
    print("AP_route=BOUNDARY-TIGHT (pigeonhole upper bound and t=1/14 equality)")
    print("replace12by26_route=STRICT-Q12-WITNESS")
    print("all_three_raw_phase_functions_equal=PASS")
    print("all_fourteen_hasse_coordinates_per_function_equal=PASS")
    print("period12_phase_functions_differ=PASS")

    print("INFINITE COLLISION FAMILY")
    print("S_h=AP with 12 replaced by 12+14h")
    print("period14_packet_equal_for_every_h=by congruence")
    print("strict_q12_witness_when_6_does_not_divide_h=since (12+14h) mod 12=2h")

    print("Q-THRESHOLD-BLIND INFINITE COLLISION FAMILY")
    for key, value in infinite_far_alias_checks().items():
        print(f"{key}={value}")
    print(f"AP_strict_primitive_numerators_q41={strict_primitive_numerators(ap, 41)}")
    print(f"T0_strict_primitive_numerators_q41={strict_primitive_numerators(far_alias(0), 41)}")

    print("FINITE 7-ADIC HEIGHT TRUNCATION NO-GO")
    for key, value in finite_seven_adic_height_no_go().items():
        print(f"{key}={value}")

    print("BOUNDED SAME-RESIDUE ALIAS AUDIT")
    for key, value in alias_audit().items():
        print(f"{key}={value}")

    print("NAMED HYP-2979 ROUTE-FIBER AUDIT")
    for key, value in named_route_fiber_audit().items():
        print(f"{key}={value}")

    print("OWNER-LABELLED CRT REPAIR")
    for key, value in owner_crt_checks().items():
        print(f"{key}={value}")
    print("owner_labelled_mod13xmod14_recovers_height=PASS")

    named = {
        "AP": ap,
        "GW12to24": replaced(12, 24),
        "residue_liar12to26": liar,
        "near12to36": replaced(12, 36),
        "petal10to20": replaced(10, 20),
        "covering12to84": replaced(12, 84),
    }
    print("NAMED FIRST PERIOD-14 JET DEPTH DIFFERING FROM AP")
    for name, speeds in named.items():
        first = next(
            (depth for depth in range(7) if jet_packet(speeds, depth) != jet_packet(ap, depth)),
            None,
        )
        print(f"{name}: first_depth={first}")

    print("CHARACTERISTIC-3 COMPONENT VERSUS ENERGY")
    for key, value in char3_component_checks().items():
        print(f"{key}={value}")

    print("TOURNAMENT ANALYSIS")
    for key, value in tournament_fingerprint().items():
        print(f"{key}={value}")
    print(
        "assumption_challenge=complete raw Hasse jets cannot recover magnitude or "
        "endpoint ownership absent from the phase function"
    )
    print(
        "verdict=PASS; Hasse jets are complete coordinates after mod-7 reduction, "
        "but a single period cannot supply LRC magnitude or exit data"
    )


if __name__ == "__main__":
    main()
