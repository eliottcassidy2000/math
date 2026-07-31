#!/usr/bin/env python3
"""Exact D-handoffs on the full supported-arrival LRC rail bank.

Retain every positive THM-2584 route-two rail as a weighted half-open set

    R_(m,s,b,t)(x),

where ``m`` is the arrival digit, ``s`` the source, ``b`` the owner clock,
and ``t`` the deep label.  The complete inherited support is

    m=0,t=0: 81 rails;  m=6,t in {0,12}: 162 rails;
    m=12,t=12: 81 rails.

For every ordered pair of supported arrivals this script first computes

    R_(m0,s0,b,t0)(x) AND R_(m1,s1,c,t1)(D x),  D(x)={13x},

on the exact grid ``13*T``.  It then cuts by the current shallow clock ``a``.
The covariance ``shallow(Dx)=owner(x)`` types the two stored clock edges as

    a -> b -> c.

Raw rail-key pairs exist only when ``m0=m1``.  After requiring both edges to
be nonconstant, the exact legal-object counts are 2,376, 4,224, and 2,376
for arrivals 0, 6, and 12.  The endpoint objects are identified piecewise by
reflection.  All six varying-arrival cells of the projected phase graph are
therefore absent before any independent phase marginal is multiplied.

The global reason is unexpectedly simple.  The unions of all positive rail
pieces at the three arrivals are the exact envelopes

    B0=[0,1/28), B6=[13/28,15/28), B12=[27/28,1).

Their D-images are ``[0,13/28)``, ``[1/28,27/28)``, and ``[15/28,1)``;
they meet the supported envelopes only on the same arrival.  The three raw
same-arrival rail unions are respectively

    [0,1/4732), [2365/4732,2367/4732), [4731/4732,1),

where ``4732=28*13^2``.  On each one the first stored edge satisfies
``c7(Dx)=c7(D^2x)`` and is diagonal.  Thus none of the 27 arrival words in
``{0,6,12}^3`` extends a legal two-edge object to a third rail, even before
the third clock edge is tested.  Endpoint alternation dies at the first
varying handoff, while all three same-arrival loops die by the same diagonal
clock mechanism.  Arrival-only and arrival/deep cuts are retained as a
diagnostic ledger locating the phase slivers discarded by the envelopes.

This is not the false statement that every raw three-rail intersection is
empty.  Explicit reflected endpoint triples have positive raw support but
force the first clock edge ``0->0``.  Nor is this a full sharp-graph event:
present packets, tooth/carry labels, delayed words, Bockstein units, endpoint
sections, and scalar exclusion are not inserted.  Product-profile weighted
lengths are exact support certificates, not canonical transition masses.
Only positive-length half-open components are counted; isolated endpoints
have measure zero and are deliberately ignored.
"""

from bisect import bisect_left
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product

import lrc14_dilation_reversed_clock_fibre_product_probe as two
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old
import lrc_central_arrival_clock_return_classification as phase


P = 13
Q7 = 7
ARRIVALS = (0, 6, 12)
DEEP = (0, 12)

EXPECTED_TWO = {
    (0, 0): (3744, 972, 540, 144, 2376, 132, ((18, 132),),
             4392, 6, 7_741_745_231_220,
             777_238_009_789_549_502_903_819_799_148_387_200),
    (0, 6): (0, 0, 0, 0, 0, 0, (), 0, 0, 0, 0),
    (0, 12): (0, 0, 0, 0, 0, 0, (), 0, 0, 0, 0),
    (6, 0): (0, 0, 0, 0, 0, 0, (), 0, 0, 0, 0),
    (6, 6): (6696, 1944, 816, 288, 4224, 144,
             ((16, 24), (32, 120)), 7726, 10, 12_714_085_664_280,
             1_274_494_381_119_049_317_956_010_741_084_837_120),
    (6, 12): (0, 0, 0, 0, 0, 0, (), 0, 0, 0, 0),
    (12, 0): (0, 0, 0, 0, 0, 0, (), 0, 0, 0, 0),
    (12, 6): (0, 0, 0, 0, 0, 0, (), 0, 0, 0, 0),
    (12, 12): (3744, 972, 540, 144, 2376, 132, ((18, 132),),
               4392, 6, 7_741_745_231_220,
               777_238_009_789_549_502_903_819_799_148_387_200),
}

EXPECTED_ARRIVAL_EXTENSIONS = {
    (0, 0, 0): (132, ((0, 1, 0),)),
    (0, 0, 6): (132, ((0, 3, 0),)),
    (6, 6, 0): (660, (
        (3, 0, 3), (3, 0, 4), (3, 0, 5), (3, 0, 6), (3, 1, 0),
    )),
    (6, 6, 12): (660, (
        (4, 0, 1), (4, 0, 2), (4, 0, 3), (4, 0, 4), (4, 6, 0),
    )),
    (12, 12, 6): (132, ((0, 4, 0),)),
    (12, 12, 12): (132, ((0, 6, 0),)),
}

EXPECTED_DEEP_EXTENSIONS = {
    (0, 0, 6): (132, ((0, 3, 0, 12),)),
    (6, 6, 0): (132, ((3, 0, 3, 0),)),
    (6, 6, 12): (132, ((4, 0, 4, 12),)),
    (12, 12, 6): (132, ((0, 4, 0, 0),)),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def weighted_mass(pieces):
    return sum(weight * (right - left) for left, right, weight in pieces)


def intersect_weighted_union(pieces, intervals):
    """Grid-independent weighted/unweighted half-open intersection."""
    out = []
    i = j = 0
    while i < len(pieces) and j < len(intervals):
        left = max(pieces[i][0], intervals[j][0])
        right = min(pieces[i][1], intervals[j][1])
        if left < right:
            out.append((left, right, pieces[i][2]))
        if pieces[i][1] <= intervals[j][1]:
            i += 1
        else:
            j += 1
    return out


def merge_support(intervals):
    """Merge a family of unweighted half-open intervals."""
    out = []
    for left, right in sorted(intervals):
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return out


def sparse_intersect_weighted(left, right, right_ends):
    """Intersect a tiny weighted list with a cached large sorted list."""
    out = []
    for a, b, wa in left:
        index = bisect_left(right_ends, a + 1)
        while index < len(right) and right[index][0] < b:
            c, d, wb = right[index]
            lo = max(a, c)
            hi = min(b, d)
            if lo < hi:
                out.append((lo, hi, wa * wb))
            index += 1
    return out


def build_rail_bank():
    e4 = old.base.build_set(old.base.PAT_E, old.base.ZELL)
    qb = old.base.build_set(old.host.PAT_QB, old.base.ZELL)
    ust, uv, vst, vv = old.rail.packet_profiles(e4, qb)
    _, _, tensor = old.rail.exact_tensor(e4, qb)
    owner = old.base.clock_cells(P, Q7, old.T, P * P)
    deep = old.rail.deep_cells()
    bank = {}
    for m in ARRIVALS:
        arrival = [(m * old.T // P, (m + 1) * old.T // P)]
        for s in range(1, P):
            rvst, rvv = old.base.rotate_profile(
                vst, vv, s * old.T // P, old.T
            )
            profile_starts, profile_values, _, _ = old.base.product_cum(
                ust, uv, rvst, rvv, old.T
            )
            for ell in range(Q7):
                for t in DEEP:
                    cell = old.intersect_sorted(
                        old.intersect_sorted(owner[ell], deep[t]), arrival
                    )
                    pieces = old.profile_on_intervals(
                        cell, profile_starts, profile_values
                    )
                    numerator = P * weighted_mass(pieces)
                    require(numerator == tensor[s][m][t][ell],
                            "rail reconstruction disagrees with tensor")
                    if numerator:
                        bank[m, s, ell, t] = pieces
    require(Counter(key[0] for key in bank) == Counter({6: 162, 0: 81, 12: 81}),
            "alternate-arrival rail census changed")
    return bank


def build_shallow_cells(refined_grid_factor=1):
    module = old.load_present_module()
    cells = []
    for ell in range(Q7):
        intervals = module.make_comb(
            module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        cells.append([
            (refined_grid_factor * left, refined_grid_factor * right)
            for left, right in intervals
        ])
    return cells


def pullback_d2_weighted(pieces):
    """Represent ``x -> {P^2*x}`` pullback on the grid ``P^2*T``."""
    return sorted(
        (left + branch * old.T, right + branch * old.T, weight)
        for left, right, weight in pieces
        for branch in range(P * P)
    )


def pullback_d2_intervals(intervals):
    return sorted(
        (left + branch * old.T, right + branch * old.T)
        for left, right in intervals
        for branch in range(P * P)
    )


def audit_clock_covariance():
    """Check ``D^-1(shallow_ell)=owner_ell`` as exact interval unions."""
    owner = old.base.clock_cells(P, Q7, old.T, P * P)
    shallow = build_shallow_cells()
    checks = 0
    for ell in range(Q7):
        scaled_owner = two.merge_adjacent(
            (P * left, P * right) for left, right in owner[ell]
        )
        pulled_shallow = two.merge_adjacent(
            (left + branch * old.T, right + branch * old.T)
            for left, right in shallow[ell]
            for branch in range(P)
        )
        require(scaled_owner == pulled_shallow,
                "D clock covariance changed")
        checks += 1
    return checks


def audit_envelopes(bank):
    """Extract the three aggregate rail envelopes and their D geometry."""
    require(old.T % 28 == 0, "rail grid stopped resolving denominator 28")
    unit = old.T // 28
    expected = {
        0: [(0, unit)],
        6: [(13 * unit, 15 * unit)],
        12: [(27 * unit, old.T)],
    }
    envelopes = {}
    for m in ARRIVALS:
        envelopes[m] = merge_support(
            (left, right)
            for key, pieces in bank.items()
            if key[0] == m
            for left, right, _ in pieces
        )
    require(envelopes == expected,
            "aggregate supported-arrival rail envelopes changed")

    # Direct images have no wrap on these three short envelopes.
    image_envelopes = {
        0: [(0, 13 * unit)],
        6: [(unit, 27 * unit)],
        12: [(15 * unit, old.T)],
    }
    diagonal = []
    for m0 in ARRIVALS:
        scaled = [(P * left, P * right, 1)
                  for left, right in envelopes[m0]]
        for m1 in ARRIVALS:
            pulled = two.pullback_dilation_weighted(
                [(left, right, 1) for left, right in envelopes[m1]], old.T
            )
            positive = bool(two.intersect_weighted(scaled, pulled))
            require(positive == (m0 == m1),
                    "envelope D-handoff stopped being arrival-diagonal")
            if positive:
                diagonal.append((m0, m1))

    # The union over all raw rail-key triples distributes over intersection,
    # so it is computed exactly from the three aggregate envelopes.
    grid2 = P * P * old.T
    expected_raw_triples = {
        0: [(0, unit)],
        6: [(2365 * unit, 2367 * unit)],
        12: [(4731 * unit, grid2)],
    }
    raw_triples = {}
    raw_clock_pairs = {}
    shallow2 = build_shallow_cells(P * P)
    owner2 = [
        [(P * P * left, P * P * right) for left, right in cell]
        for cell in old.base.clock_cells(P, Q7, old.T, P * P)
    ]
    for m in ARRIVALS:
        base = [(left, right, 1) for left, right in envelopes[m]]
        first = two.scale_weighted(base, P * P)
        second = two.scale_weighted(
            two.pullback_dilation_weighted(base, old.T), P
        )
        third = pullback_d2_weighted(base)
        joint = two.intersect_weighted(
            two.intersect_weighted(first, second), third
        )
        raw_triples[m] = merge_support(
            (left, right) for left, right, _ in joint
        )
        pairs = []
        for shallow_clock in range(Q7):
            for owner_clock in range(Q7):
                cell = old.intersect_sorted(
                    old.intersect_sorted(
                        raw_triples[m], shallow2[shallow_clock]
                    ),
                    owner2[owner_clock],
                )
                if cell:
                    pairs.append((shallow_clock, owner_clock))
        raw_clock_pairs[m] = tuple(pairs)
    require(raw_triples == expected_raw_triples,
            "raw same-arrival three-rail envelopes changed")
    require(raw_clock_pairs == {
        0: ((0, 0),),
        6: ((3, 3), (4, 4)),
        12: ((0, 0),),
    }, "raw three-rail envelope acquired a nonconstant first edge")
    return envelopes, image_envelopes, tuple(diagonal), raw_triples, raw_clock_pairs


def raw_constituent_triple_controls(bank):
    """Positive raw triples showing why the legal qualifier is essential."""
    grid2 = P * P * old.T
    expected = {
        (0, 1, 0, 0): [(0, 125_882_036_280,
                         20_996_345_682_699_618_893_895_762_446_016)],
        (12, 1, 0, 12): [(
            50_334_309_852_666_840,
            50_334_435_734_703_120,
            20_996_345_682_699_618_893_895_762_446_016,
        )],
    }
    shallow0 = build_shallow_cells(P * P)[0]
    checks = 0
    for key, target in expected.items():
        pieces = bank[key]
        joint = two.intersect_weighted(
            two.intersect_weighted(
                two.scale_weighted(pieces, P * P),
                two.scale_weighted(
                    two.pullback_dilation_weighted(pieces, old.T), P
                ),
            ),
            pullback_d2_weighted(pieces),
        )
        require(joint == target,
                "raw endpoint constituent triple control changed")
        require(intersect_weighted_union(joint, shallow0) == joint,
                "raw endpoint triple left shallow clock zero")
        require(key[2] == 0,
                "raw endpoint positive control stopped being clock-diagonal")
        require(0 <= joint[0][0] < joint[0][1] <= grid2,
                "raw endpoint control left refined circle")
        checks += 1
    return checks


def two_rail_handoffs(bank):
    shallow = build_shallow_cells(P)
    by_arrival = {
        m: tuple((key, pieces) for key, pieces in sorted(bank.items())
                 if key[0] == m)
        for m in ARRIVALS
    }
    summary = {}
    all_objects = defaultdict(list)
    for m0 in ARRIVALS:
        for m1 in ARRIVALS:
            source_counts = Counter()
            clock_counts = Counter()
            raw_pairs = 0
            clock_refined_objects = 0
            zero_first = 0
            zero_second = 0
            zero_both = 0
            legal_objects = 0
            total_weighted_length = 0
            all_support = []
            for key0, rail0 in by_arrival[m0]:
                _, s0, owner0, t0 = key0
                scaled0 = two.scale_weighted(rail0, P)
                for key1, rail1 in by_arrival[m1]:
                    _, s1, owner1, t1 = key1
                    joint = two.intersect_weighted(
                        scaled0,
                        two.pullback_dilation_weighted(rail1, old.T),
                    )
                    if not joint:
                        continue
                    raw_pairs += 1
                    for shallow0 in range(Q7):
                        typed = intersect_weighted_union(
                            joint, shallow[shallow0]
                        )
                        if not typed:
                            continue
                        clock_refined_objects += 1
                        # E0 has shallow0 -> owner0; E1 has
                        # owner0 -> owner1 under the D covariance.
                        first_zero = shallow0 == owner0
                        second_zero = owner0 == owner1
                        zero_first += first_zero
                        zero_second += second_zero
                        zero_both += first_zero and second_zero
                        if first_zero or second_zero:
                            continue
                        legal_objects += 1
                        mass = weighted_mass(typed)
                        require(mass > 0, "positive rail handoff lost mass")
                        total_weighted_length += mass
                        all_support.extend(
                            (left, right) for left, right, _ in typed
                        )
                        source_counts[s0, s1] += 1
                        clock_counts[shallow0, owner0, owner1] += 1
                        all_objects[m0, m1].append((
                            key0, key1, shallow0, typed
                        ))
            components = merge_support(all_support)
            summary[m0, m1] = {
                "raw_pairs": raw_pairs,
                "clock_refined_objects": clock_refined_objects,
                "zero_first": zero_first,
                "zero_second": zero_second,
                "zero_both": zero_both,
                "legal_objects": legal_objects,
                "source_hist": Counter(source_counts.values()),
                "positive_source_pairs": len(source_counts),
                "clock_paths": tuple(sorted(clock_counts)),
                "clock_occurrence_hist": Counter(clock_counts.values()),
                "weighted_length": total_weighted_length,
                "piece_occurrences": len(all_support),
                "geometric_components": len(components),
                "geometric_length": sum(
                    right - left for left, right in components
                ),
            }
    return summary, all_objects


def third_arrival_deep_filters(bank, two_objects):
    """Locate the layer at which a possible third supported rail dies."""
    supported_deep = {
        m: tuple(sorted({key[3] for key in bank if key[0] == m}))
        for m in ARRIVALS
    }
    require(supported_deep == {0: (0,), 6: (0, 12), 12: (12,)},
            "supported arrival/deep ledger changed")
    deep_cells = old.rail.deep_cells()
    summary = {}
    for triple in product(ARRIVALS, repeat=3):
        m0, m1, m2 = triple
        arrival_interval = [(
            m2 * old.T // P,
            (m2 + 1) * old.T // P,
        )]
        pulled_arrival = pullback_d2_intervals(arrival_interval)
        arrival_objects = 0
        arrival_paths = Counter()
        deep_objects = 0
        deep_paths = Counter()
        for key0, key1, shallow0, pair_pieces in two_objects[m0, m1]:
            scaled_pair = two.scale_weighted(pair_pieces, P)
            arrival_cut = intersect_weighted_union(
                scaled_pair, pulled_arrival
            )
            if not arrival_cut:
                continue
            arrival_objects += 1
            path = (shallow0, key0[2], key1[2])
            arrival_paths[path] += 1
            # Since m1 is the arrival of z=Dx and m2 that of Dz, the
            # current edge (shallow0,owner0) must lie in this exact phase
            # support.  This is an independent typing control.
            require(
                (shallow0, key0[2])
                in phase.forward_cross_arrival_support(Q7, m1, m2),
                "third-arrival cut disagrees with phase-return typing",
            )
            for deep_label in supported_deep[m2]:
                base_cell = old.intersect_sorted(
                    arrival_interval, deep_cells[deep_label]
                )
                deep_cut = intersect_weighted_union(
                    scaled_pair, pullback_d2_intervals(base_cell)
                )
                if deep_cut:
                    deep_objects += 1
                    deep_paths[path + (deep_label,)] += 1
        summary[triple] = {
            "arrival_objects": arrival_objects,
            "arrival_paths": tuple(sorted(arrival_paths)),
            "arrival_path_hist": Counter(arrival_paths.values()),
            "deep_objects": deep_objects,
            "deep_paths": tuple(sorted(deep_paths)),
            "deep_path_hist": Counter(deep_paths.values()),
        }
    actual_arrival = {
        triple: (value["arrival_objects"], value["arrival_paths"])
        for triple, value in summary.items()
        if value["arrival_objects"]
    }
    actual_deep = {
        triple: (value["deep_objects"], value["deep_paths"])
        for triple, value in summary.items()
        if value["deep_objects"]
    }
    require(actual_arrival == EXPECTED_ARRIVAL_EXTENSIONS,
            "third-arrival extension atlas changed")
    require(actual_deep == EXPECTED_DEEP_EXTENSIONS,
            "third-arrival/deep extension atlas changed")
    require(all(
        value["arrival_path_hist"] == Counter({132: len(value["arrival_paths"])})
        for value in summary.values() if value["arrival_objects"]
    ), "third-arrival path multiplicities changed")
    require(all(
        value["deep_path_hist"] == Counter({132: len(value["deep_paths"])})
        for value in summary.values() if value["deep_objects"]
    ), "third-deep path multiplicities changed")
    return summary


def three_rail_handoffs(bank, two_objects, arrival_triples):
    """Extend physical two-rail objects by one further D handoff.

    A two-rail piece is represented on the grid ``P*T``.  Scaling it by P
    and intersecting with the full D^2-pullback of a third base rail gives
    the exact common-x support on ``P^2*T``.  No independent phase-return
    interval is multiplied here.
    """
    by_arrival = {
        m: tuple((key, pieces) for key, pieces in sorted(bank.items())
                 if key[0] == m)
        for m in ARRIVALS
    }
    pull2_cache = {}
    pull2_identity_checks = 0
    for key, pieces in bank.items():
        pulled = pullback_d2_weighted(pieces)
        sequential = two.pullback_dilation_weighted(
            two.pullback_dilation_weighted(pieces, old.T), P * old.T
        )
        require(pulled == sequential,
                "direct and sequential D^2 pullbacks disagree")
        pull2_identity_checks += 1
        pull2_cache[key] = (
            pulled,
            tuple(right for _, right, _ in pulled),
        )

    summary = {}
    all_objects = defaultdict(list)
    sparse_controls = 0
    for m0, m1, m2 in arrival_triples:
        source_counts = Counter()
        clock_counts = Counter()
        extensions_of_legal_two = 0
        legal_objects = 0
        total_weighted_length = 0
        all_support = []
        for key0, key1, shallow0, pair_pieces in two_objects[m0, m1]:
            _, s0, owner0, t0 = key0
            _, s1, owner1, t1 = key1
            scaled_pair = two.scale_weighted(pair_pieces, P)
            for key2, _ in by_arrival[m2]:
                _, s2, owner2, t2 = key2
                pulled2, pulled2_ends = pull2_cache[key2]
                joint = sparse_intersect_weighted(
                    scaled_pair, pulled2, pulled2_ends
                )
                if sparse_controls < 1200:
                    require(
                        joint == two.intersect_weighted(
                            scaled_pair, pulled2
                        ),
                        "sparse third-rail intersection disagrees with sweep",
                    )
                    sparse_controls += 1
                if not joint:
                    continue
                extensions_of_legal_two += 1
                # The third event has shallow=owner1 by D covariance.
                if owner1 == owner2:
                    continue
                legal_objects += 1
                mass = weighted_mass(joint)
                require(mass > 0, "positive three-rail handoff lost mass")
                total_weighted_length += mass
                source_counts[s0, s1, s2] += 1
                clock_counts[shallow0, owner0, owner1, owner2] += 1
                all_support.extend((left, right) for left, right, _ in joint)
                all_objects[m0, m1, m2].append((
                    key0, key1, key2, shallow0, joint
                ))
        components = merge_support(all_support)
        summary[m0, m1, m2] = {
            "extensions_of_legal_two": extensions_of_legal_two,
            "legal_objects": legal_objects,
            "positive_source_triples": len(source_counts),
            "source_hist": Counter(source_counts.values()),
            "clock_paths": tuple(sorted(clock_counts)),
            "clock_occurrence_hist": Counter(clock_counts.values()),
            "weighted_length": total_weighted_length,
            "piece_occurrences": len(all_support),
            "geometric_components": len(components),
            "geometric_length": sum(right - left for left, right in components),
        }
    require(pull2_identity_checks == len(bank) == 324,
            "D^2 pullback identity census changed")
    require(sparse_controls == 1200,
            "sparse/ordinary intersection control census changed")
    return summary, all_objects, pull2_identity_checks, sparse_controls


def audit_two_rail_atlas(summary, objects):
    expected_paths = {
        (0, 0): tuple(
            (0, owner, following)
            for owner in (1, 2, 3)
            for following in range(Q7)
            if following != owner
        ),
        (6, 6): (
            (3, 0, 3), (3, 0, 4), (3, 0, 5), (3, 0, 6),
            (3, 1, 0), (3, 1, 2), (3, 1, 3), (3, 1, 4),
            (3, 1, 5), (3, 1, 6),
            (3, 2, 0), (3, 2, 1), (3, 2, 3), (3, 2, 4),
            (3, 2, 5), (3, 2, 6),
            (4, 0, 1), (4, 0, 2), (4, 0, 3), (4, 0, 4),
            (4, 5, 0), (4, 5, 1), (4, 5, 2), (4, 5, 3),
            (4, 5, 4), (4, 5, 6),
            (4, 6, 0), (4, 6, 1), (4, 6, 2), (4, 6, 3),
            (4, 6, 4), (4, 6, 5),
        ),
        (12, 12): tuple(
            (0, owner, following)
            for owner in (4, 5, 6)
            for following in range(Q7)
            if following != owner
        ),
    }
    source_missing = {}
    for key, value in summary.items():
        actual = (
            value["raw_pairs"],
            value["zero_first"],
            value["zero_second"],
            value["zero_both"],
            value["legal_objects"],
            value["positive_source_pairs"],
            tuple(sorted(value["source_hist"].items())),
            value["piece_occurrences"],
            value["geometric_components"],
            value["geometric_length"],
            value["weighted_length"],
        )
        require(actual == EXPECTED_TWO[key],
                f"two-rail atlas changed at arrivals {key}")
        require(value["clock_refined_objects"] == value["raw_pairs"],
                "a raw rail-key pair crossed two shallow clock cells")
        require(
            value["legal_objects"]
            == value["raw_pairs"]
            - value["zero_first"] - value["zero_second"]
            + value["zero_both"],
            "two-edge rejection inclusion-exclusion failed",
        )
        expected = expected_paths.get(key, ())
        require(value["clock_paths"] == expected,
                f"clock-path support changed at arrivals {key}")
        require(value["clock_occurrence_hist"]
                == (Counter({132: len(expected)}) if expected else Counter()),
                "clock-path multiplicity changed")

        positive_sources = {
            (obj[0][1], obj[1][1]) for obj in objects[key]
        }
        source_missing[key] = tuple(
            (s0, s1)
            for s0 in range(1, P)
            for s1 in range(1, P)
            if (s0, s1) not in positive_sources
        )
    require(source_missing[0, 0]
            == tuple((s, 7) for s in range(1, P)),
            "arrival-zero following-source hole changed")
    require(source_missing[6, 6] == (),
            "central source-pair coverage stopped being complete")
    require(source_missing[12, 12]
            == tuple((s, 6) for s in range(1, P)),
            "arrival-twelve following-source hole changed")

    # Exact a.e. reflection on the refined grid, including every weight.
    refined_grid = P * old.T

    def reflected_key(key):
        m, source, owner, deep = key
        return 12 - m, P - source, (-owner) % Q7, 12 - deep

    def reflected_object(obj):
        key0, key1, shallow0, pieces = obj
        return (
            reflected_key(key0),
            reflected_key(key1),
            (-shallow0) % Q7,
            tuple(sorted(
                (refined_grid - right, refined_grid - left, weight)
                for left, right, weight in pieces
            )),
        )

    reflection_checks = 0
    for key in ((0, 0), (6, 6), (12, 12)):
        target = (12 - key[0], 12 - key[1])
        reflected = {reflected_object(obj) for obj in objects[key]}
        target_objects = {
            (key0, key1, shallow0, tuple(pieces))
            for key0, key1, shallow0, pieces in objects[target]
        }
        require(reflected == target_objects,
                f"two-rail reflection failed at arrivals {key}")
        reflection_checks += len(objects[key])
    require(reflection_checks == 8976,
            "two-rail reflection object census changed")
    return expected_paths, source_missing, reflection_checks


def main():
    bank = build_rail_bank()
    clock_covariance_checks = audit_clock_covariance()
    (envelopes, image_envelopes, diagonal_arrivals,
     raw_triples, raw_clock_pairs) = audit_envelopes(bank)
    raw_constituent_checks = raw_constituent_triple_controls(bank)
    summary, objects = two_rail_handoffs(bank)
    expected_paths, source_missing, reflection_checks = audit_two_rail_atlas(
        summary, objects
    )
    third_filters = third_arrival_deep_filters(bank, objects)
    arrival_triples = tuple(product(ARRIVALS, repeat=3))
    (three_summary, three_objects,
     pull2_identity_checks, sparse_controls) = three_rail_handoffs(
        bank, objects, arrival_triples
    )
    require(len(three_summary) == 27,
            "supported three-arrival universe changed")
    require(all(
        value["extensions_of_legal_two"] == 0
        and value["legal_objects"] == 0
        and value["weighted_length"] == 0
        for value in three_summary.values()
    ), "a legal two-edge rail object acquired a third rail")
    require(not three_objects,
            "zero three-rail atlas unexpectedly stored an object")

    pair_rows = tuple(
        (
            key,
            summary[key]["raw_pairs"],
            summary[key]["zero_first"],
            summary[key]["zero_second"],
            summary[key]["zero_both"],
            summary[key]["legal_objects"],
            summary[key]["positive_source_pairs"],
            tuple(sorted(summary[key]["source_hist"].items())),
            summary[key]["piece_occurrences"],
            summary[key]["geometric_components"],
        )
        for key in sorted(summary)
    )
    weighted_rows = tuple(
        (key, summary[key]["weighted_length"])
        for key in sorted(summary)
    )
    geometric_rows = tuple(
        (
            key,
            Fraction(summary[key]["geometric_length"], P * old.T),
        )
        for key in sorted(summary)
    )
    clock_rows = tuple(
        (key, expected_paths.get(key, ())) for key in sorted(summary)
    )
    arrival_filter_rows = tuple(sorted(
        (key, value) for key, value in EXPECTED_ARRIVAL_EXTENSIONS.items()
    ))
    deep_filter_rows = tuple(sorted(
        (key, value) for key, value in EXPECTED_DEEP_EXTENSIONS.items()
    ))
    require(all(
        third_filters[key]["arrival_objects"] == value[0]
        for key, value in EXPECTED_ARRIVAL_EXTENSIONS.items()
    ) and all(
        third_filters[key]["deep_objects"] == value[0]
        for key, value in EXPECTED_DEEP_EXTENSIONS.items()
    ), "printed third-filter atlas disagrees with its audit")

    print("LRC14 exact supported-arrival physical rail D-handoff atlas")
    print("scope=all 324 positive THM-2584 route-two weighted rail profiles on arrivals {0,6,12}")
    print("handoff=D(x)={13x}; rail_key=(arrival,source,owner,deep); typed_edges=(shallow0,owner0,owner1)")
    print(
        "aggregate_rail_envelopes="
        "B0=[0,1/28),B6=[13/28,15/28),B12=[27/28,1)"
    )
    print(
        "D_images="
        "D(B0)=[0,13/28),D(B6)=[1/28,27/28),"
        "D(B12)=[15/28,1)"
    )
    print(f"positive_envelope_arrival_handoffs={diagonal_arrivals}")
    print(
        "two_rail_rows_"
        "(arrivals,raw_keys,zero_first,zero_second,zero_both,"
        "legal_objects,source_pairs,source_hist,piece_occurrences,"
        f"geometric_components)={pair_rows}"
    )
    print(f"two_rail_clock_path_atlas={clock_rows}")
    print(
        "same_arrival_missing_source_pairs="
        f"{tuple((key, source_missing[key]) for key in ((0,0),(6,6),(12,12)))}"
    )
    print(
        f"legal_two_rail_weighted_length_numerators_over_{P * old.T}="
        f"{weighted_rows}"
    )
    print(f"legal_two_rail_geometric_union_measures={geometric_rows}")
    print(f"piecewise_reflection_object_checks={reflection_checks}")
    print(
        "raw_same_arrival_three_rail_unions="
        "000:[0,1/4732),666:[2365/4732,2367/4732),"
        "12^3:[4731/4732,1)"
    )
    print(f"raw_three_rail_first_clock_pairs={tuple(sorted(raw_clock_pairs.items()))}")
    print(f"positive_clock_diagonal_raw_constituent_triple_controls={raw_constituent_checks}")
    print(f"third_arrival_only_extensions={arrival_filter_rows}")
    print(f"third_arrival_deep_extensions={deep_filter_rows}")
    print(
        "full_supported_three_arrival_words=27; "
        "extensions_of_legal_two_edge_objects_before_third_edge=0"
    )
    print(
        f"controls=clock_D_covariance:{clock_covariance_checks},"
        f"D2_direct_equals_sequential:{pull2_identity_checks},"
        f"sparse_equals_standard:{sparse_controls}"
    )
    print(
        "verdict=PASS: the supported rail bank is arrival-diagonal under D, "
        "and every raw same-arrival three-rail loop forces a diagonal first edge"
    )
    print(
        "semantics=positive-length a.e. rail-clock nilpotence only; raw clock-"
        "illegal triples exist, and no present/tooth/word/unit/endpoint/scalar "
        "conclusion is asserted"
    )


if __name__ == "__main__":
    main()
