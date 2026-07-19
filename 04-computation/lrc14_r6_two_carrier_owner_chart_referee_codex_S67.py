#!/usr/bin/env python3
"""Exact referee for THM-1154's residue-owner chart and two-carrier closure.

The general chart based at a normalized carrier ``s`` is

    t_a = a/13 + c/(182s),          a=1,...,12,

under the broad cap ``c*M <= 13*s``.  Its core obstruction is exactly the
set of distinct residue owners ``{-p^(-1): c*p>s}``.  For an exceptional
normalized carrier pair ``13u < 13v`` the closing specialization is

    t_a = a/13 + 15/(182v),       a=1,...,12.

All decisions below use integers or ``Fraction``.  The finite lift box is a
replay rather than the proof: the script separately enumerates every cell of
an arbitrarily translated thirteen-grid and thereby verifies the universal
two-point noncarrier capacity.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import comb


THRESHOLD = Fraction(1, 14)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circle_distance(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def safe(speed: int, time: Fraction) -> bool:
    return circle_distance(speed * time) >= THRESHOLD


def exceptional_triples() -> list[tuple[int, int, int]]:
    return [
        (M, M + 1, v)
        for M in range(8, 13)
        for v in range(13 * M + 14, 15 * M)
    ]


def general_chart_time(a: int, s: int, c: int) -> Fraction:
    require(1 <= a <= 12, "multiplier outside F_13^*")
    require(s > 0 and c > 0, "source and switch must be positive")
    return Fraction(a, 13) + Fraction(c, 182 * s)


def chart_time(a: int, v: int) -> Fraction:
    return general_chart_time(a, v, 15)


def general_core_bad_mask(P: tuple[int, ...], s: int, c: int) -> int:
    mask = 0
    for a in range(1, 13):
        if any(not safe(p, general_chart_time(a, s, c)) for p in P):
            mask |= 1 << (a - 1)
    return mask


def predicted_owner_mask(P: tuple[int, ...], s: int, c: int) -> int:
    mask = 0
    for p in P:
        if c * p > s:
            owner = next(a for a in range(1, 13) if (a * p) % 13 == 12)
            mask |= 1 << (owner - 1)
    return mask


def core_bad_mask(P: tuple[int, ...], v: int) -> int:
    return general_core_bad_mask(P, v, 15)


def speed_bad_mask(k: int, v: int) -> int:
    mask = 0
    for a in range(1, 13):
        if not safe(k, chart_time(a, v)):
            mask |= 1 << (a - 1)
    return mask


def normalized_carrier_safe(c: int, o: int, s: int) -> bool:
    """Return dist(c*o/s,14Z)>=1 after clearing the positive denominator."""
    modulus = 14 * s
    remainder = (c * o) % modulus
    return min(remainder, modulus - remainder) >= s


def translated_grid_bad_mask(theta: Fraction) -> int:
    """Bad b in {1,...,12} for b/13+theta, independent of a labels."""
    mask = 0
    for b in range(1, 13):
        if circle_distance(Fraction(b, 13) + theta) < THRESHOLD:
            mask |= 1 << (b - 1)
    return mask


def translated_grid_cell_audit() -> tuple[int, int, int]:
    """Enumerate every endpoint and open cell of the translation circle."""
    breaks = {
        frac_part(sign * THRESHOLD - Fraction(b, 13))
        for b in range(1, 13)
        for sign in (-1, 1)
    }
    ordered = sorted(breaks)
    probes = set(ordered)
    for i, left in enumerate(ordered):
        right = ordered[(i + 1) % len(ordered)]
        if i + 1 == len(ordered):
            right += 1
        probes.add(frac_part((left + right) / 2))

    masks = {translated_grid_bad_mask(theta) for theta in probes}
    max_load = max(mask.bit_count() for mask in masks)
    require(len(breaks) == 24, "unexpected translated-grid breakpoint collision")
    require(max_load == 2, "a translated thirteen-grid used more than two charts")
    return len(breaks), len(probes), len(masks)


def abstract_hyperedge_audit() -> tuple[tuple[int, int, int, tuple[int, ...]], ...]:
    """Exhaust owner/noncarrier capacities for rho=2,...,6 carriers."""
    masks = [0]
    masks.extend(1 << i for i in range(12))
    masks.extend((1 << i) | (1 << j) for i, j in combinations(range(12), 2))
    require(len(masks) == 79, "size-at-most-two mask universe is wrong")

    ledger: list[tuple[int, int, int, tuple[int, ...]]] = []
    for rho in range(2, 7):
        edge_count = 6 - rho
        owner_masks = tuple((1 << d) - 1 for d in range(8))
        maxima = [0] * 8
        assignments = 0
        for indices in combinations_with_replacement(range(len(masks)), edge_count):
            edge_union = 0
            for index in indices:
                edge_union |= masks[index]
            for d, owner in enumerate(owner_masks):
                maxima[d] = max(maxima[d], (owner | edge_union).bit_count())
            assignments += 1

        expected_assignments = comb(79 + edge_count - 1, edge_count)
        require(assignments == expected_assignments, "multiset assignment count mismatch")
        for d, maximum in enumerate(maxima):
            expected_maximum = min(12, d + 2 * edge_count)
            require(maximum == expected_maximum, "abstract owner/edge capacity mismatch")
            if d < 2 * rho:
                require(maximum < 12, "the d<2rho survival condition failed")
        ledger.append((rho, edge_count, assignments, tuple(maxima)))

    require(ledger[0][2] == 1749060, "four-edge assignment count changed")
    require(ledger[0][3][1] == 9, "two-carrier one-owner capacity changed")
    return tuple(ledger)


def general_owner_cap_audit(s_max: int = 60) -> tuple[int, int, int, int]:
    """Replay every broad-cap core chart in a normalized finite box."""
    cap_rows = 0
    core_rows = 0
    phase_checks = 0
    maximum_owners = 0

    for M in range(7, 13):
        cores = [tuple(low) + (M,) for low in combinations(range(1, M), 6)]
        for s in range(M + 1, s_max + 1):
            for c in range(1, (13 * s) // M + 1):
                require(c * M <= 13 * s, "broad core cap rounding failed")
                cap_rows += 1
                individual_masks: dict[int, int] = {}

                for p in range(1, M + 1):
                    actual = 0
                    for a in range(1, 13):
                        residue = (a * p) % 13
                        numerator = 14 * s * residue + c * p
                        require(14 * s <= numerator <= 181 * s,
                                "broad-cap phase left the no-wrap strip")
                        integer_safe = 13 * s <= numerator <= 169 * s
                        metric_safe = safe(p, general_chart_time(a, s, c))
                        require(integer_safe == metric_safe, "cleared core phase disagrees")
                        if not metric_safe:
                            actual |= 1 << (a - 1)
                        phase_checks += 1

                    if c * p > s:
                        owner = next(a for a in range(1, 13) if (a * p) % 13 == 12)
                        predicted = 1 << (owner - 1)
                    else:
                        predicted = 0
                    require(actual == predicted, "individual residue-owner law failed")
                    individual_masks[p] = actual

                crossed = [p for p in range(1, M + 1) if c * p > s]
                crossed_masks = [individual_masks[p] for p in crossed]
                require(len(set(crossed_masks)) == len(crossed_masks),
                        "distinct core speeds shared a boundary owner")
                crossed_union = 0
                for mask in crossed_masks:
                    require((crossed_union & mask) == 0, "boundary-owner masks overlap")
                    crossed_union |= mask
                require(crossed_union.bit_count() == len(crossed),
                        "owner union cardinality is not the number crossed")

                for P in cores:
                    actual = 0
                    for p in P:
                        actual |= individual_masks[p]
                    predicted = predicted_owner_mask(P, s, c)
                    require(actual == predicted, "seven-core owner union failed")
                    d = sum(c * p > s for p in P)
                    require(actual.bit_count() == d, "core owner cardinality is not d")
                    maximum_owners = max(maximum_owners, d)
                    core_rows += 1

    require(maximum_owners == 7, "finite broad-cap box did not realize seven owners")
    return cap_rows, core_rows, phase_checks, maximum_owners


def general_carrier_chart_audit(o_max: int = 20) -> tuple[int, int, tuple[int, ...]]:
    """Check normalized carrier condition (G3) against literal metric phases."""
    carrier_sets = 0
    candidate_charts = 0
    safe_histogram = [0] * 7

    for M in range(7, 13):
        universe = range(M + 1, o_max + 1)
        for rho in range(1, min(6, o_max - M) + 1):
            for carriers in combinations(universe, rho):
                carrier_sets += 1
                for s in carriers:
                    for c in range(1, (13 * s) // M + 1):
                        candidate_charts += 1
                        predicted = all(normalized_carrier_safe(c, o, s) for o in carriers)
                        t = general_chart_time(1, s, c)
                        actual = all(safe(13 * o, t) for o in carriers)
                        require(actual == predicted, "normalized carrier condition disagrees")
                        if actual:
                            safe_histogram[rho] += 1

    require(all(safe_histogram[rho] > 0 for rho in range(1, 7)),
            "finite carrier box missed a carrier cardinality")
    return carrier_sets, candidate_charts, tuple(safe_histogram[1:])


def core_and_carrier_audit() -> tuple[int, int, int]:
    triples = exceptional_triples()
    require(len(triples) == 30, "exceptional triple count is not thirty")
    core_rows = 0
    chart_checks = 0

    for M, u, v in triples:
        z = 13 * u
        h = 13 * v
        require(u == M + 1, "exception has the wrong small carrier")
        require(13 * u < v < 15 * u, "carrier ratio orientation failed")
        ratio = Fraction(15 * u, v)
        require(1 < ratio < 2, "small-carrier phase is outside (1,2)")
        require(15 * (M - 1) <= v < 15 * M, "owner threshold is misoriented")
        require(15 * M < 14 * v, "residue -1 phase could wrap through one")

        for a in range(1, 13):
            t = chart_time(a, v)
            dz = circle_distance(z * t)
            dh = circle_distance(h * t)
            require(THRESHOLD < dz < Fraction(1, 7), "small carrier sign/wrap failed")
            require(dh == THRESHOLD, "large carrier is not at the safe endpoint")
            chart_checks += 1

        for low in combinations(range(1, M), 6):
            P = tuple(low) + (M,)
            core_rows += 1
            actual = core_bad_mask(P, v)
            owner_a = next(a for a in range(1, 13) if (a * M) % 13 == 12)
            predicted = 1 << (owner_a - 1)
            require(actual == predicted, f"core owner mismatch at {(M, u, v, P)}")

            for a in range(1, 13):
                t = chart_time(a, v)
                for p in P:
                    r = (a * p) % 13
                    y = Fraction(15 * p, 182 * v)
                    require(r != 0, "core acquired residue zero")
                    require(Fraction(r, 13) + y < 1, "core phase wrapped through one")
                    is_bad = not safe(p, t)
                    predicted_bad = p == M and r == 12
                    require(is_bad == predicted_bad, "core sign/orientation ledger failed")

    require(core_rows == 6930, "exception/core row count mismatch")
    require(chart_checks == 360, "carrier chart count mismatch")
    return len(triples), core_rows, chart_checks


def noncarrier_lift_audit(k_max: int = 5000) -> tuple[int, int, int]:
    rows = 0
    maximum = 0
    realized_masks: set[int] = set()
    for M, _u, v in exceptional_triples():
        for k in range(13 * M + 1, k_max + 1):
            if k % 13 == 0:
                continue
            mask = speed_bad_mask(k, v)
            load = mask.bit_count()
            require(load <= 2, f"noncarrier used {load} charts at {(M, v, k)}")
            rows += 1
            maximum = max(maximum, load)
            realized_masks.add(mask)
    require(maximum == 2, "finite lift box did not realize the sharp load two")
    return rows, maximum, len(realized_masks)


def displayed_example_audit() -> tuple[Fraction, tuple[int, ...], tuple[int, ...]]:
    P = (1, 2, 7, 9, 10, 11, 12)
    K = (157, 158, 159, 160, 169, 2210)
    v = 170
    survivors: list[int] = []
    for a in range(1, 13):
        t = chart_time(a, v)
        if all(safe(w, t) for w in P + K):
            survivors.append(a)
    require(tuple(survivors) == (2, 5, 7, 8, 9, 10, 11), "example survivors changed")

    t = chart_time(2, v)
    distances = tuple(circle_distance(w * t) for w in P + K)
    minimum = min(distances)
    owners = tuple(w for w, distance in zip(P + K, distances) if distance == minimum)
    require(t == Fraction(955, 6188), "displayed time changed")
    require(minimum == THRESHOLD, "displayed minimum changed")
    require(owners == (2210,), "displayed equality owner changed")
    return t, tuple(survivors), owners


def tournament_audit() -> tuple[tuple[int, ...], int]:
    """Load tournament for the displayed chart; tie reversal is telemetry only."""
    P = (1, 2, 7, 9, 10, 11, 12)
    noncarriers = (157, 158, 159, 160)
    v = 170
    core = core_bad_mask(P, v)
    edges = [speed_bad_mask(k, v) for k in noncarriers]
    loads = tuple(
        ((core >> i) & 1) + sum((edge >> i) & 1 for edge in edges)
        for i in range(12)
    )

    order_up = sorted(range(12), key=lambda i: (loads[i], i))
    order_down = sorted(range(12), key=lambda i: (loads[i], -i))
    position_up = {vertex: index for index, vertex in enumerate(order_up)}
    position_down = {vertex: index for index, vertex in enumerate(order_down)}
    flips = sum(
        (position_up[i] < position_up[j]) != (position_down[i] < position_down[j])
        for i in range(12)
        for j in range(i + 1, 12)
    )
    require(sorted(position_up.values()) == list(range(12)), "tournament order failed")
    require(flips == sum(comb(loads.count(level), 2) for level in set(loads)),
            "tie-gauge flip count failed")
    return loads, flips


def main() -> None:
    triples, core_rows, carrier_checks = core_and_carrier_audit()
    breakpoints, probes, geometric_masks = translated_grid_cell_audit()
    hyperedge_ledger = abstract_hyperedge_audit()
    assignments = hyperedge_ledger[0][2]
    abstract_maximum = hyperedge_ledger[0][3][1]
    cap_rows, general_core_rows, phase_checks, maximum_owners = general_owner_cap_audit()
    carrier_sets, candidate_charts, safe_histogram = general_carrier_chart_audit()
    lift_rows, lift_maximum, lift_masks = noncarrier_lift_audit()
    time, survivors, owners = displayed_example_audit()
    loads, tie_flips = tournament_audit()

    print("THM-1154 exact two-carrier residue-owner closure")
    print(f"exceptional triples checked: {triples}")
    print(f"seven-core/exception rows checked: {core_rows}")
    print(f"carrier chart checks: {carrier_checks}")
    print("carrier orientation: 1 < 15u/v < 2; h is safe at equality 1/14")
    print("core orientation: only aM = -1 (mod 13) is bad; exactly 1 of 12 charts")
    print("general owner chart:")
    print(f"  broad-cap rows={cap_rows}, seven-core rows={general_core_rows}, phase checks={phase_checks}")
    print(f"  maximum distinct core owners realized={maximum_owners}")
    print(f"  carrier sets={carrier_sets}, candidate charts={candidate_charts}, safe-by-rho={safe_histogram}")
    for rho, edge_count, assignment_count, maxima in hyperedge_ledger:
        print(f"  rho={rho}: noncarriers={edge_count}, assignments={assignment_count}, max-unions-by-d={maxima}")
    print("  exact survival criterion: d + 2(6-rho) < 12 iff d < 2rho")
    print(f"translated-grid cells: breakpoints={breakpoints}, probes={probes}, masks={geometric_masks}")
    print("universal noncarrier capacity: at most 2 of 12 charts")
    print(f"finite noncarrier lift replay: rows={lift_rows}, max_load={lift_maximum}, masks={lift_masks}")
    print(f"abstract four-edge assignments: {assignments}; max with core owner={abstract_maximum}")
    print("survival lower bound: 12 - 1 - 4*2 = 3")
    print(f"displayed witness: t={time}, survivors={survivors}, equality_owners={owners}")

    print("Tournament Analysis:")
    print("  two-carrier order tournament: scores=(0,1), cycles=0, SCCs=(1,1), HP=1")
    print(f"  displayed chart obstruction loads={loads}")
    print("  load tournament is transitive: scores=0..11, cycles=0, singleton SCCs, HP=1")
    print(f"  reversing chart-label tie gauge flips {tie_flips} edges")
    print("  challenged vertices: runners, carriers, ratios, arcs, owners, charts, obligations")
    print("  faithful object: 12 chart vertices plus one core singleton and four <=2 hyperedges")
    print("  pair orientation destroys hyperedge union capacity; killer order is unnecessary")
    print("VERDICT: all 30 switch obstructions close; the exactly-two-carrier branch is empty")


if __name__ == "__main__":
    main()
