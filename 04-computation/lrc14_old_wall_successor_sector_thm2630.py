#!/usr/bin/env python3
"""Exact companion for THM-2630: old-wall clutching versus successor carry.

There are two different translated deep gates on the canonical common-x
packet.  The old THM-2587 gate has speed two in the route-five coordinate;
the later THM-2600 probe has speed ``2*13**5``.  This script first checks the
affine old-wall formula on the proved THM-2586/2600 rails, and then splits
every later probe tooth into its exact left/right wall halves on THM-2616's
raw q=h diagonal.  The latter support relation is emphatically not a map.
"""

from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_deep_root_coshift_incidence_wall_thm2587 as wallmod
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old


P = 13
Q7 = 7
SHARDS = ((0, 41), (41, 82), (82, 123), (123, 162))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def bit_count(mask):
    return mask.bit_count()


def selected_theta(s, ell):
    """THM-2600's uniform q=0 rail selector."""
    return int(s in (6, 11) or (s, ell) == (8, 2))


def old_wall_census():
    """Check the old speed-two wall before and after the full pullback."""
    E = wallmod.base.build_set(wallmod.base.PAT_E, wallmod.base.ZELL)
    Q = wallmod.base.build_set(wallmod.host.PAT_QB, wallmod.base.ZELL)
    _, wall, _, _ = wallmod.root_chart_route(E, Q)

    # THM-2586's priority theta-zero selector is root zero except for the
    # three exact fallback cells, where it uses root twelve.
    fallback = {(7, 4), (7, 5), (7, 6)}
    priority_patterns = Counter()
    priority_middle_zero = []
    for s in range(1, P):
        for ell in range(Q7):
            root = 12 if (s, ell) in fallback else 0
            pattern = tuple(wall[0][zone][s][root][ell] > 0
                            for zone in range(3))
            priority_patterns[pattern] += 1
            if not pattern[1]:
                priority_middle_zero.append((s, ell))
    require(
        priority_patterns
        == Counter({(True, True, False): 44,
                    (False, True, False): 36,
                    (True, False, False): 4}),
        "THM-2586 priority wall census changed",
    )
    require(
        priority_middle_zero == [(7, 0), (7, 1), (7, 2), (7, 3)],
        "THM-2586 pure-low boundary changed",
    )

    # THM-2600 keeps arrival v=6 and chooses t=12 on theta zero or t=0 on
    # theta one.  Its old wall has a positive middle stratum in all 84 cells.
    selector_patterns = Counter()
    naive_failures = []
    corrected_failures = []
    for s in range(1, P):
        for ell in range(Q7):
            theta = selected_theta(s, ell)
            root = 0 if theta else 12
            pattern = tuple(wall[theta][zone][s][root][ell] > 0
                            for zone in range(3))
            selector_patterns[(theta, pattern)] += 1
            # h=v=6 and inverse(2)=7 in F_13.
            if 7 * root % P != 6:
                naive_failures.append((s, ell))
            if 7 * (root - theta) % P != 6:
                corrected_failures.append((s, ell))
    require(
        selector_patterns
        == Counter({(0, (False, True, True)): 39,
                    (0, (False, True, False)): 30,
                    (1, (True, True, False)): 8,
                    (1, (False, True, False)): 7}),
        "THM-2600 selected old-wall census changed",
    )
    require(len(naive_failures) == 15 and not corrected_failures,
            "theta correction boundary changed")

    # Algebraic equivariance of h=7(tau-epsilon-theta).  Normalize tau by
    # tau_bar=7*tau; under v -> v+a, both h and tau_bar shift by a.
    for theta in (0, 1):
        for v in range(P):
            root = (2 * v + theta) % P
            for epsilon in (0, 1):
                tau = (root + epsilon) % P
                require(7 * (tau - epsilon - theta) % P == v,
                        "old-wall affine formula failed")
                for a in range(P):
                    tau_shift = (tau + 2 * a) % P
                    require(
                        7 * (tau_shift - epsilon - theta) % P
                        == (v + a) % P,
                        "normalized old-wall equivariance failed",
                    )

    # Strict danger seams.  At z=1/14 only epsilon zero survives; at z=13/14
    # only epsilon one survives.  In the open middle both survive.
    def incident_offsets(z):
        return tuple(
            epsilon for epsilon in (0, 1)
            if (epsilon == 0 and z < Fraction(13, 14))
            or (epsilon == 1 and z > Fraction(1, 14))
        )

    require(incident_offsets(Fraction(1, 14)) == (0,),
            "lower strict seam changed")
    require(incident_offsets(Fraction(1, 2)) == (0, 1),
            "open middle lost an edge")
    require(incident_offsets(Fraction(13, 14)) == (1,),
            "upper strict seam changed")

    # Refine THM-2600's selected q=0, h=6 physical carrier by the old middle
    # wall before integrating.  Positivity here is not inferred from separate
    # marginals.  We do not claim that the restricted Bockstein remains a unit.
    (module, prefixes, _, _, rails, present, starts) = cross.build_carrier_data()
    rail_by_key = {(s, ell, t): pieces for s, ell, t, pieces in rails}
    old_danger = [wallmod.danger_intervals(2, tau) for tau in range(P)]
    full_middle_support = []
    for s in range(1, P):
        for ell4 in range(Q7):
            theta = selected_theta(s, ell4)
            root = 0 if theta else 12
            pieces = rail_by_key[s, ell4, root]
            middle = wallmod.root_tensor.intersect_lists(
                old_danger[root], old_danger[(root + 1) % P]
            )
            middle_pieces = old.intersect_weighted_union(
                pieces, middle, [left for left, _ in middle]
            )
            support = 0
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    middle_pieces, present[ell5, 0], starts[ell5, 0]
                )
                for probe in range(1, P):
                    probed = old.intersect_weighted_comb(
                        overlap, module.C3, 182,
                        14 * probe - 13, 14 * probe + 13,
                    )
                    value = old.delayed_weighted_numerator(
                        probed, prefixes[ell5][6]
                    )
                    support += int(value > 0)
            require(support > 0, "full common-x old middle wall vanished")
            full_middle_support.append(support)
    full_middle_hist = Counter(full_middle_support)
    require(
        full_middle_hist == Counter({36: 37, 24: 23, 12: 12, 48: 12}),
        "full common-x old-middle support census changed",
    )
    return (
        priority_patterns,
        tuple(priority_middle_zero),
        selector_patterns,
        tuple(naive_failures),
        full_middle_hist,
    )


def later_probe_shard(bounds):
    """Split later Delta_r teeth and retain support on the raw q=h diagonal."""
    start, stop = bounds
    (module, prefixes, _, _, rails, present, starts) = cross.build_carrier_data()
    require(0 <= start < stop <= len(rails), "invalid later-probe shard")

    records = []
    clock_records = []
    ambiguous_rails = 0
    partition_checks = 0
    for s, ell4, root, pieces in rails[start:stop]:
        theta = (root - 12) % P
        require(theta in (0, 1), "middle rail left its theta chart")
        support = defaultdict(int)
        clock_support = defaultdict(int)
        for h in range(P):
            shift = (-h) % P
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    pieces, present[ell5, shift], starts[ell5, shift]
                )
                if not overlap:
                    continue
                for probe in range(1, P):
                    whole = old.intersect_weighted_comb(
                        overlap, module.C3, 182,
                        14 * probe - 13, 14 * probe + 13,
                    )
                    whole_value = old.delayed_weighted_numerator(
                        whole, prefixes[ell5][h]
                    )
                    split_values = []
                    # epsilon=0 is the right half with absolute root=probe;
                    # epsilon=1 is the left half with absolute root=probe-1.
                    for epsilon, left, right in (
                        (0, 14 * probe, 14 * probe + 13),
                        (1, 14 * probe - 13, 14 * probe),
                    ):
                        half = old.intersect_weighted_comb(
                            overlap, module.C3, 182, left, right
                        )
                        value = old.delayed_weighted_numerator(
                            half, prefixes[ell5][h]
                        )
                        split_values.append(value)
                        if value:
                            support[probe, epsilon] |= 1 << h
                            clock_support[ell5, probe, epsilon] |= 1 << h
                    require(sum(split_values) == whole_value,
                            "left/right later tooth did not partition")
                    partition_checks += 1
        rail_records = tuple(
            ((s, ell4, theta, probe, epsilon), mask)
            for (probe, epsilon), mask in sorted(support.items())
        )
        ambiguous_rails += int(any(bit_count(mask) > 1
                                   for _, mask in rail_records))
        records.extend(rail_records)
        clock_records.extend(
            ((s, ell4, theta, ell5, probe, epsilon), mask)
            for (ell5, probe, epsilon), mask
            in sorted(clock_support.items())
        )
    return (bounds, tuple(records), tuple(clock_records), ambiguous_rails,
            partition_checks)


def strict_local_hostiles():
    """Check the two displayed strict hostiles and the seven-sector minimum."""
    examples = []
    for y, expected_h in ((Fraction(11, 260), 7),
                          (Fraction(3, 52), 9)):
        phase = (2 * y) % 1
        absolute = int(P * phase)
        z = P * phase - absolute
        probe = 1
        epsilon = (probe - absolute) % P
        u = (P * y) % 1
        predecessor = int(P * y)
        kappa = int(2 * u)
        h = int(P * u)
        require(
            (probe, epsilon, kappa, h)
            == (1, 0, 1, expected_h),
            "strict local hostile labels changed",
        )
        require(Fraction(1, 14) < z < Fraction(13, 14),
                "strict local hostile left the open middle")
        require(absolute == (2 * predecessor + kappa) % P,
                "later predecessor-carry formula failed")
        require(abs(phase - Fraction(probe, P)) < Fraction(1, 14),
                "strict local hostile left the translated danger tooth")
        require(h == int(P * (kappa + z) / 2),
                "successor-sector formula failed")
        examples.append((y, z, h))

    # Even after kappa and the three-state middle wall are retained, the fixed
    # later signature (probe,epsilon)=(1,0) realizes seven h sectors on each
    # half.  Therefore no additional bit can recover h; the seven-sector
    # refinement (equivalently h itself after gluing the halves) is minimal.
    half_sector_support = []
    for kappa in (0, 1):
        values = set()
        for h in (range(0, 7) if kappa == 0 else range(6, 13)):
            lower = max(
                Fraction(h, P),
                Fraction(kappa, 2),
                Fraction(kappa, 2) + Fraction(1, 28),
            )
            upper = min(
                Fraction(h + 1, P),
                Fraction(kappa + 1, 2),
                Fraction(kappa, 2) + Fraction(13, 28),
            )
            require(lower < upper, "expected middle successor sector empty")
            u = (lower + upper) / 2
            v = 7 if kappa == 0 else 0
            y = (v + u) / P
            phase = (2 * y) % 1
            absolute = int(P * phase)
            z = P * phase - absolute
            require(absolute == 1 and int(2 * u) == kappa,
                    "seven-sector witness changed absolute root")
            require(absolute == (2 * v + kappa) % P,
                    "seven-sector predecessor-carry formula failed")
            require(Fraction(1, 14) < z < Fraction(13, 14),
                    "seven-sector witness left the middle wall")
            require(int(P * u) == h,
                    "seven-sector witness changed successor digit")
            values.add(h)
        require(len(values) == 7, "half-wall lost a successor sector")
        half_sector_support.append(tuple(sorted(values)))
    return tuple(examples), tuple(half_sector_support)


def main():
    (priority_patterns, priority_middle_zero, selector_patterns,
     naive_failures, full_middle_hist) = old_wall_census()
    hostiles, half_sector_support = strict_local_hostiles()

    with ProcessPoolExecutor(max_workers=len(SHARDS)) as pool:
        shard_results = list(pool.map(later_probe_shard, SHARDS))
    require(tuple(result[0] for result in shard_results) == SHARDS,
            "later-probe shards returned out of order")

    signatures = {}
    clock_signatures = {}
    ambiguous_rails = 0
    partition_checks = 0
    for _, records, clock_records, ambiguous, checks in shard_results:
        ambiguous_rails += ambiguous
        partition_checks += checks
        for key, mask in records:
            require(key not in signatures, "duplicate later-probe signature")
            signatures[key] = mask
        for key, mask in clock_records:
            require(key not in clock_signatures,
                    "duplicate clock-retained later-probe signature")
            clock_signatures[key] = mask

    signature_hist = Counter(
        (key[2], bit_count(mask)) for key, mask in signatures.items()
    )
    ambiguous_hist = Counter(
        (key[2], bit_count(mask))
        for key, mask in signatures.items() if bit_count(mask) > 1
    )
    ambiguous_signatures = sum(bit_count(mask) > 1
                               for mask in signatures.values())
    functional_signatures = len(signatures) - ambiguous_signatures
    ambiguous_cells = {
        (s, ell) for (s, ell, theta, probe, epsilon), mask
        in signatures.items() if bit_count(mask) > 1
    }
    future_mask = sum(1 << h for h in range(1, 12))

    # Retaining the complete future-owner clock still leaves almost every
    # later tooth signature nonfunctional.  This is a strictly finer test
    # than the clock-united census above.
    clock_hist = Counter(bit_count(mask)
                         for mask in clock_signatures.values())
    clock_ambiguous = sum(bit_count(mask) > 1
                          for mask in clock_signatures.values())
    clock_ambiguous_cells = {
        (s, ell4) for (s, ell4, theta, ell5, probe, epsilon), mask
        in clock_signatures.items() if bit_count(mask) > 1
    }

    # Test every affine coefficient graph r=alpha*h+beta against this fully
    # clock-labelled physical support.  The identically-zero graph is the
    # forbidden absent probe and is excluded.  More importantly, no affine
    # bijection (alpha nonzero) is functional at even one visible signature.
    visible = sorted({
        (s, ell4, theta, ell5, epsilon)
        for s, ell4, theta, ell5, probe, epsilon in clock_signatures
    })
    graph_rows = []
    for alpha in range(P):
        for beta in range(P):
            if (alpha, beta) == (0, 0):
                continue
            sizes = []
            for s, ell4, theta, ell5, epsilon in visible:
                mask = 0
                for h in range(P):
                    probe = (alpha * h + beta) % P
                    if probe and clock_signatures.get(
                        (s, ell4, theta, ell5, probe, epsilon), 0
                    ) >> h & 1:
                        mask |= 1 << h
                sizes.append(bit_count(mask))
            graph_rows.append((
                alpha, beta, max(sizes), min(sizes),
                tuple(sorted(Counter(sizes).items())),
            ))
    bijective_graph_rows = [row for row in graph_rows if row[0] != 0]
    graph_max_hist = Counter(row[2] for row in bijective_graph_rows)
    opposite_graph = next(row for row in graph_rows
                          if row[:2] == (P - 1, P - 1))

    require(partition_checks == 61_248,
            "later-tooth partition-check universe changed")
    require(len(signatures) == 3_792,
            "later-probe signature census changed")
    require(ambiguous_signatures == 3_476
            and functional_signatures == 316,
            "later-probe ambiguity census changed")
    require(ambiguous_rails == 158,
            "ambiguous later-rail census changed")
    require(len(ambiguous_cells) == 84,
            "a base cell unexpectedly acquired a later-probe map")
    require(all(mask & ~future_mask == 0 for mask in signatures.values()),
            "later signature left THM-2616 future support")
    require(
        ambiguous_hist
        == Counter({(1, 11): 1452, (0, 11): 1254,
                    (0, 10): 484, (1, 10): 286}),
        "later-probe ambiguous support histogram changed",
    )
    require(
        len(clock_signatures) == 10_480
        and clock_hist
        == Counter({10: 2948, 9: 2904, 11: 1980,
                    7: 1452, 1: 712, 8: 484}),
        "clock-retained signature census changed",
    )
    require(clock_ambiguous == 9_768
            and len(clock_ambiguous_cells) == 84,
            "future clock unexpectedly repaired the successor map")
    require(len(visible) == 888 and len(graph_rows) == 168
            and len(bijective_graph_rows) == 156,
            "affine graph universe changed")
    require(all(row[3] >= 5 for row in bijective_graph_rows)
            and graph_max_hist == Counter({9: 74, 10: 74, 11: 8}),
            "an affine bijection became locally functional")
    require(
        opposite_graph
        == (12, 12, 11, 7,
            ((7, 132), (8, 143), (9, 266), (10, 257), (11, 90))),
        "opposite affine graph spectrum changed",
    )

    print("THM-2630 exact old-wall/successor-sector controls")
    print(f"priority_wall_patterns={sorted(priority_patterns.items())}")
    print(f"priority_middle_positive=80 pure_low_cells={priority_middle_zero}")
    print(f"thm2600_selector_patterns={sorted(selector_patterns.items())}")
    print(f"naive_theta0_formula_failures={len(naive_failures)} corrected_failures=0")
    print(f"full_common_x_middle_support_hist={sorted(full_middle_hist.items())}")
    print("old_wall_formula=h=qprime=7*(tau-epsilon-theta)_mod13")
    print(f"strict_hostiles={hostiles}")
    print(f"middle_successor_sectors_by_kappa={half_sector_support}")
    print(f"later_partition_checks={partition_checks} ambiguous_rails={ambiguous_rails}/162")
    print(f"later_signature_hist={sorted(signature_hist.items())}")
    print(f"later_ambiguous_hist={sorted(ambiguous_hist.items())}")
    print(f"later_signatures={len(signatures)} ambiguous={ambiguous_signatures} functional={functional_signatures}")
    print(f"later_ambiguous_base_cells={len(ambiguous_cells)}/84")
    print(f"clock_retained_hist={sorted(clock_hist.items())}")
    print(f"clock_retained_signatures={len(clock_signatures)} ambiguous={clock_ambiguous} base_cells={len(clock_ambiguous_cells)}/84")
    print(f"affine_bijections={len(bijective_graph_rows)} visible_signatures={len(visible)} support_floor={min(row[3] for row in bijective_graph_rows)} max_support_hist={sorted(graph_max_hist.items())}")
    print(f"opposite_graph_physical_spectrum={opposite_graph[4]}")
    print("verdict=PASS: old wall clutches affinely; later probe needs successor-sector/carry data")


if __name__ == "__main__":
    main()
