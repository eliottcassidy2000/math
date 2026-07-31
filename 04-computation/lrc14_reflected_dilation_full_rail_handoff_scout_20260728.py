#!/usr/bin/env python3
"""Exact reflected-dilation scout on the complete THM-2584 rail bank.

Let ``D(x)={13x}``, let ``R(x)={-x}`` be circle reflection, and put
``T=R D=D R``.  On the THM-2684 rail bank, reflection acts by

    (m,s,b,t) -> (12-m,13-s,-b,12-t).

If ``S`` and ``O=S o D`` are the stored shallow and owner clocks, then

    S(Tx)=-O(x).

Thus ``T`` is not a same-label replacement for the canonical ``D`` handoff:
the existing boundary equation ``O_i=S_(i+1)`` first fails by a sign.  It is
an exact handoff only after declaring the reflection involution
``J(b)=-b`` on every interface.  Equivalently, reflect every odd event.
That alternating gauge sends a ``T`` chain to a ``D`` chain exactly.

This script checks those statements using exact half-open weighted rail
pieces, classifies all nine arrival pairs and all 27 arrival triples, and
separately records the accidental fixed-clock-zero part on which the old
untwisted interface equation holds.  It does not add present packets,
delayed words, units, semantic endpoints, or scalar rows.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
import hashlib

import lrc14_alternate_arrival_physical_rail_handoff as full


P = full.P
Q7 = full.Q7
ARRIVALS = full.ARRIVALS
old = full.old
two = full.two
EXPECTED_OBJECT_DIGEST = (
    "25711579f77d93db7e8fa9cb99484194dfbca6b570b5c88a146a334570c8e019"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def reflected_clock(label):
    return (-label) % Q7


def reflected_arrival(label):
    return 12 - label


def reflected_key(key):
    arrival, source, owner, deep = key
    return (
        reflected_arrival(arrival),
        P - source,
        reflected_clock(owner),
        12 - deep,
    )


def reflect_weighted(pieces, grid):
    """A.e. circle reflection of sorted half-open weighted pieces."""
    return sorted(
        (grid - right, grid - left, weight)
        for left, right, weight in pieces
    )


def reflect_intervals(intervals, grid):
    return sorted((grid - right, grid - left) for left, right in intervals)


def pullback_t_weighted(pieces, grid):
    """Represent ``T^-1`` on the refined grid ``13*grid``."""
    return two.pullback_dilation_weighted(
        reflect_weighted(pieces, grid), grid
    )


def pullback_t_intervals(intervals, grid):
    weighted = pullback_t_weighted(
        [(left, right, 1) for left, right in intervals], grid
    )
    return [(left, right) for left, right, _ in weighted]


def audit_reflection_and_covariance(bank):
    rail_checks = 0
    for key, pieces in bank.items():
        target = reflected_key(key)
        require(target in bank, "reflection left the complete rail bank")
        require(reflect_weighted(pieces, old.T) == bank[target],
                "weighted rail reflection changed")
        require(reflected_key(target) == key,
                "rail reflection stopped being involutive")
        rail_checks += 1

    shallow = full.build_shallow_cells()
    owner = old.base.clock_cells(P, Q7, old.T, P * P)
    shallow_reflection_checks = 0
    owner_reflection_checks = 0
    t_covariance_checks = 0
    for label in range(Q7):
        target = reflected_clock(label)
        require(reflect_intervals(shallow[label], old.T) == shallow[target],
                "shallow clock reflection changed")
        require(reflect_intervals(owner[label], old.T) == owner[target],
                "owner clock reflection changed")
        shallow_reflection_checks += 1
        owner_reflection_checks += 1

        # S(Tx)=-O(x), equivalently T^-1(S_-b)=O_b.
        scaled_owner = two.merge_adjacent(
            (P * left, P * right) for left, right in owner[label]
        )
        pulled_reflected_shallow = two.merge_adjacent(
            pullback_t_intervals(shallow[target], old.T)
        )
        require(scaled_owner == pulled_reflected_shallow,
                "reflected-dilation clock covariance changed")
        t_covariance_checks += 1

    t2_checks = 0
    for pieces in bank.values():
        sequential = pullback_t_weighted(
            pullback_t_weighted(pieces, old.T), P * old.T
        )
        direct = full.pullback_d2_weighted(pieces)
        require(sequential == direct,
                "T^2 stopped agreeing with D^2")
        t2_checks += 1
    return (
        rail_checks,
        shallow_reflection_checks,
        owner_reflection_checks,
        t_covariance_checks,
        t2_checks,
    )


def audit_envelopes(bank):
    (envelopes, d_images, _, d_raw_triples,
     d_raw_clock_pairs) = full.audit_envelopes(bank)

    expected_t_images = {
        0: [(15 * old.T // 28, old.T)],
        6: [(old.T // 28, 27 * old.T // 28)],
        12: [(0, 13 * old.T // 28)],
    }
    t_images = {
        arrival: reflect_intervals(d_images[arrival], old.T)
        for arrival in ARRIVALS
    }
    require(t_images == expected_t_images,
            "reflected-dilation envelope images changed")

    positive_pairs = []
    for first in ARRIVALS:
        scaled = [(P * left, P * right, 1)
                  for left, right in envelopes[first]]
        for second in ARRIVALS:
            pulled = pullback_t_weighted(
                [(left, right, 1) for left, right in envelopes[second]],
                old.T,
            )
            positive = bool(two.intersect_weighted(scaled, pulled))
            expected = second == reflected_arrival(first)
            require(positive == expected,
                    "T arrival adjacency stopped being anti-diagonal")
            if positive:
                positive_pairs.append((first, second))

    raw_triples = {}
    raw_clock_pairs = {}
    grid2 = P * P * old.T
    shallow2 = full.build_shallow_cells(P * P)
    owner2 = [
        [(P * P * left, P * P * right) for left, right in cell]
        for cell in old.base.clock_cells(P, Q7, old.T, P * P)
    ]
    for word in product(ARRIVALS, repeat=3):
        first, second, third = word
        first_pieces = two.scale_weighted(
            [(left, right, 1) for left, right in envelopes[first]],
            P * P,
        )
        second_pieces = two.scale_weighted(
            pullback_t_weighted(
                [(left, right, 1) for left, right in envelopes[second]],
                old.T,
            ),
            P,
        )
        third_pieces = pullback_t_weighted(
            pullback_t_weighted(
                [(left, right, 1) for left, right in envelopes[third]],
                old.T,
            ),
            P * old.T,
        )
        joint = two.intersect_weighted(
            two.intersect_weighted(first_pieces, second_pieces),
            third_pieces,
        )
        support = full.merge_support(
            (left, right) for left, right, _ in joint
        )
        if not support:
            continue
        raw_triples[word] = support
        pairs = []
        for shallow_label in range(Q7):
            for owner_label in range(Q7):
                typed = old.intersect_sorted(
                    old.intersect_sorted(support, shallow2[shallow_label]),
                    owner2[owner_label],
                )
                if typed:
                    pairs.append((shallow_label, owner_label))
        raw_clock_pairs[word] = tuple(pairs)

    expected_words = ((0, 12, 0), (6, 6, 6), (12, 0, 12))
    require(tuple(raw_triples) == expected_words,
            "positive T triple-word support changed")
    for first, second, third in expected_words:
        require(raw_triples[first, second, third] == d_raw_triples[first],
                "alternating T triple stopped matching the D return cylinder")
        require(raw_clock_pairs[first, second, third]
                == d_raw_clock_pairs[first],
                "T raw triple acquired a non-D first clock edge")
    require(raw_clock_pairs == {
        (0, 12, 0): ((0, 0),),
        (6, 6, 6): ((3, 3), (4, 4)),
        (12, 0, 12): ((0, 0),),
    }, "a positive T triple acquired a nonconstant first stored edge")
    return (
        envelopes,
        t_images,
        tuple(positive_pairs),
        raw_triples,
        raw_clock_pairs,
    )


def two_rail_t_handoffs(bank, d_summary, d_objects):
    shallow = full.build_shallow_cells(P)
    by_arrival = {
        arrival: tuple(
            (key, pieces) for key, pieces in sorted(bank.items())
            if key[0] == arrival
        )
        for arrival in ARRIVALS
    }
    summary = {}
    objects = defaultdict(list)
    bijection_checks = 0
    for first in ARRIVALS:
        for second in ARRIVALS:
            raw_pairs = 0
            clock_refined = 0
            zero_first = 0
            zero_second = 0
            zero_both = 0
            twisted_legal = 0
            twisted_paths = Counter()
            source_counts = Counter()
            weighted_length = 0
            all_support = []
            untwisted_interface_raw = 0
            untwisted_interface_legal = 0
            untwisted_paths = Counter()
            untwisted_source_counts = Counter()

            for key0, rail0 in by_arrival[first]:
                _, source0, owner0, _ = key0
                scaled0 = two.scale_weighted(rail0, P)
                for key1, rail1 in by_arrival[second]:
                    _, source1, owner1, _ = key1
                    joint = two.intersect_weighted(
                        scaled0, pullback_t_weighted(rail1, old.T)
                    )
                    if not joint:
                        continue
                    raw_pairs += 1

                    reflected1 = reflected_key(key1)
                    d_joint = two.intersect_weighted(
                        scaled0,
                        two.pullback_dilation_weighted(
                            bank[reflected1], old.T
                        ),
                    )
                    require(joint == d_joint,
                            "T pair failed its odd-event D gauge identity")
                    bijection_checks += 1

                    for shallow0 in range(Q7):
                        typed = full.intersect_weighted_union(
                            joint, shallow[shallow0]
                        )
                        if not typed:
                            continue
                        clock_refined += 1
                        shallow1 = reflected_clock(owner0)
                        reflected_owner1 = reflected_clock(owner1)
                        first_zero = shallow0 == owner0
                        second_zero = shallow1 == owner1
                        require(second_zero
                                == (owner0 == reflected_owner1),
                                "twisted second-edge typing changed")
                        zero_first += first_zero
                        zero_second += second_zero
                        zero_both += first_zero and second_zero

                        untwisted_interface = owner0 == shallow1
                        require(untwisted_interface == (owner0 == 0),
                                "odd-clock fixed locus changed")
                        untwisted_interface_raw += int(untwisted_interface)

                        if first_zero or second_zero:
                            continue
                        twisted_legal += 1
                        mass = full.weighted_mass(typed)
                        require(mass > 0,
                                "positive T rail handoff lost mass")
                        weighted_length += mass
                        all_support.extend(
                            (left, right) for left, right, _ in typed
                        )
                        source_counts[source0, source1] += 1
                        twisted_paths[
                            shallow0, owner0, reflected_owner1
                        ] += 1
                        objects[first, second].append((
                            key0, key1, shallow0, typed
                        ))
                        if untwisted_interface:
                            untwisted_interface_legal += 1
                            untwisted_paths[
                                shallow0, owner0, shallow1, owner1
                            ] += 1
                            untwisted_source_counts[source0, source1] += 1

            components = full.merge_support(all_support)
            summary[first, second] = {
                "raw_pairs": raw_pairs,
                "clock_refined": clock_refined,
                "zero_first": zero_first,
                "zero_second": zero_second,
                "zero_both": zero_both,
                "twisted_legal": twisted_legal,
                "twisted_paths": tuple(sorted(twisted_paths)),
                "twisted_path_hist": Counter(twisted_paths.values()),
                "positive_source_pairs": len(source_counts),
                "source_hist": Counter(source_counts.values()),
                "weighted_length": weighted_length,
                "piece_occurrences": len(all_support),
                "geometric_components": len(components),
                "geometric_length": sum(
                    right - left for left, right in components
                ),
                "untwisted_interface_raw": untwisted_interface_raw,
                "untwisted_interface_legal": untwisted_interface_legal,
                "untwisted_paths": tuple(sorted(untwisted_paths)),
                "untwisted_path_hist": Counter(untwisted_paths.values()),
                "untwisted_positive_source_pairs": len(
                    untwisted_source_counts
                ),
                "untwisted_source_hist": Counter(
                    untwisted_source_counts.values()
                ),
            }

            d_key = (first, reflected_arrival(second))
            reference = d_summary[d_key]
            require(
                (
                    raw_pairs,
                    clock_refined,
                    zero_first,
                    zero_second,
                    zero_both,
                    twisted_legal,
                    len(source_counts),
                    Counter(source_counts.values()),
                    weighted_length,
                    len(all_support),
                    len(components),
                    sum(right - left for left, right in components),
                    tuple(sorted(twisted_paths)),
                    Counter(twisted_paths.values()),
                )
                == (
                    reference["raw_pairs"],
                    reference["clock_refined_objects"],
                    reference["zero_first"],
                    reference["zero_second"],
                    reference["zero_both"],
                    reference["legal_objects"],
                    reference["positive_source_pairs"],
                    reference["source_hist"],
                    reference["weighted_length"],
                    reference["piece_occurrences"],
                    reference["geometric_components"],
                    reference["geometric_length"],
                    reference["clock_paths"],
                    reference["clock_occurrence_hist"],
                ),
                "T two-rail atlas stopped being alternating-gauge D",
            )

            mapped_objects = {
                (key0, reflected_key(key1), shallow0, tuple(pieces))
                for key0, key1, shallow0, pieces
                in objects[first, second]
            }
            reference_objects = {
                (key0, key1, shallow0, tuple(pieces))
                for key0, key1, shallow0, pieces in d_objects[d_key]
            }
            require(mapped_objects == reference_objects,
                    "T object bank failed exact alternating-gauge descent")

    return summary, objects, bijection_checks


def audit_t_semantics(summary):
    positive_twisted = tuple(
        key for key in sorted(summary) if summary[key]["twisted_legal"]
    )
    require(positive_twisted == ((0, 12), (6, 6), (12, 0)),
            "twisted T pair support stopped being anti-diagonal")

    positive_untwisted = tuple(
        key for key in sorted(summary)
        if summary[key]["untwisted_interface_legal"]
    )
    require(positive_untwisted == ((6, 6),),
            "untwisted fixed-clock-zero subcarrier changed")

    expected_untwisted_paths = (
        (3, 0, 0, 1), (3, 0, 0, 2),
        (3, 0, 0, 3), (3, 0, 0, 4),
        (4, 0, 0, 3), (4, 0, 0, 4),
        (4, 0, 0, 5), (4, 0, 0, 6),
    )
    require(all(
        summary[key]["untwisted_interface_raw"] == 0
        and summary[key]["untwisted_interface_legal"] == 0
        for key in summary
        if key not in ((0, 12), (6, 6), (12, 0))
    ), "an arrival-empty T cell acquired an untwisted interface")
    require(
        (
            summary[0, 12]["untwisted_interface_raw"],
            summary[0, 12]["untwisted_interface_legal"],
            summary[12, 0]["untwisted_interface_raw"],
            summary[12, 0]["untwisted_interface_legal"],
        ) == (972, 0, 972, 0),
        "endpoint fixed-clock-zero hostile changed",
    )
    require(
        (
            summary[6, 6]["untwisted_interface_raw"],
            summary[6, 6]["untwisted_interface_legal"],
            summary[6, 6]["untwisted_positive_source_pairs"],
            summary[6, 6]["untwisted_source_hist"],
            summary[6, 6]["untwisted_paths"],
            summary[6, 6]["untwisted_path_hist"],
        ) == (
            1056,
            1056,
            144,
            Counter({4: 24, 8: 120}),
            expected_untwisted_paths,
            Counter({132: 8}),
        ),
        "central untwisted fixed-clock-zero atlas changed",
    )

    # Two consecutive untwisted interfaces force b=-b and c=-c, hence
    # b=c=0.  The middle reflected-gauge edge b->c is then diagonal.
    nominal_twisted_paths = 0
    nominal_untwisted_three_paths = 0
    for a, b, c, d in product(range(Q7), repeat=4):
        if a != b and b != c and c != d:
            nominal_twisted_paths += 1
            if b == reflected_clock(b) and c == reflected_clock(c):
                nominal_untwisted_three_paths += 1
    require(nominal_twisted_paths == Q7 * 6**3,
            "nominal nonconstant clock-path universe changed")
    require(nominal_untwisted_three_paths == 0,
            "two untwisted T interfaces acquired a legal clock path")
    return (
        positive_twisted,
        positive_untwisted,
        nominal_twisted_paths,
        nominal_untwisted_three_paths,
    )


def main():
    bank = full.build_rail_bank()
    reflection_controls = audit_reflection_and_covariance(bank)
    (envelopes, t_images, arrival_pairs,
     raw_triples, raw_clock_pairs) = audit_envelopes(bank)

    d_summary, d_objects = full.two_rail_handoffs(bank)
    full.audit_two_rail_atlas(d_summary, d_objects)
    t_summary, t_objects, bijection_checks = two_rail_t_handoffs(
        bank, d_summary, d_objects
    )
    (positive_twisted, positive_untwisted,
     nominal_twisted_paths,
     nominal_untwisted_three_paths) = audit_t_semantics(t_summary)
    require(reflection_controls == (324, 7, 7, 7, 324),
            "reflection/covariance control census changed")
    require(bijection_checks == 14_184,
            "T-to-D positive raw-pair bijection census changed")

    object_serialization = ";".join(
        repr((pair, key0, key1, shallow0, tuple(pieces)))
        for pair in sorted(t_objects)
        for key0, key1, shallow0, pieces in t_objects[pair]
    )
    object_digest = hashlib.sha256(
        object_serialization.encode("ascii")
    ).hexdigest()
    require(object_digest == EXPECTED_OBJECT_DIGEST,
            "T rail-object serialization digest changed")

    pair_rows = tuple(
        (
            key,
            value["raw_pairs"],
            value["zero_first"],
            value["zero_second"],
            value["zero_both"],
            value["twisted_legal"],
            value["positive_source_pairs"],
            tuple(sorted(value["source_hist"].items())),
            value["piece_occurrences"],
            value["geometric_components"],
        )
        for key, value in sorted(t_summary.items())
    )
    twisted_clock_rows = tuple(
        (key, value["twisted_paths"])
        for key, value in sorted(t_summary.items())
        if value["twisted_paths"]
    )
    untwisted_rows = tuple(
        (
            key,
            value["untwisted_interface_raw"],
            value["untwisted_interface_legal"],
            value["untwisted_positive_source_pairs"],
            tuple(sorted(value["untwisted_source_hist"].items())),
            value["untwisted_paths"],
            tuple(sorted(value["untwisted_path_hist"].items())),
        )
        for key, value in sorted(t_summary.items())
        if value["untwisted_interface_raw"]
    )
    weighted_rows = tuple(
        (key, value["weighted_length"])
        for key, value in sorted(t_summary.items())
    )
    geometric_rows = tuple(
        (
            key,
            Fraction(value["geometric_length"], P * old.T),
        )
        for key, value in sorted(t_summary.items())
    )
    raw_triple_rows = tuple(
        (word, tuple(
            (Fraction(left, P * P * old.T),
             Fraction(right, P * P * old.T))
            for left, right in intervals
        ))
        for word, intervals in raw_triples.items()
    )

    print("LRC14 reflected-dilation full-rail handoff scout")
    print("status=VERIFIED-EXACT SCOUT; not a theorem dependency")
    print("scope=all 324 positive THM-2584 rails; T(x)={-13x}=R o D; positive-length a.e. half-open supports")
    print("reflection_key=(m,s,b,t)->(12-m,13-s,-b,12-t)")
    print("clock_covariance=S(Tx)=-O(x); a.e._boundary=T^-1(S_-b)=O_b")
    print("semantic_first_failure=existing untwisted boundary O_i=S_(i+1) becomes O_i=-S_(i+1); equality survives only at clock 0")
    print("lawful_reframe=J-twisted interface with J(b)=-b, equivalently reflect every odd event; this maps T chains exactly to D chains")
    print("rail_envelopes=B0:[0,1/28),B6:[13/28,15/28),B12:[27/28,1)")
    print("T_images=B0:[15/28,1),B6:[1/28,27/28),B12:[0,13/28)")
    print(f"positive_arrival_pairs={arrival_pairs}")
    print("two_rail_columns=(arrivals,raw_keys,zero_first,zero_second,zero_both,J_twisted_legal,source_pairs,source_hist,piece_occurrences,geometric_components)")
    print(f"two_rail_rows={pair_rows}")
    print(f"twisted_clock_path_atlas={twisted_clock_rows}")
    print(f"untwisted_fixed_clock_rows={untwisted_rows}")
    print(f"twisted_weighted_length_numerators_over_{P * old.T}={weighted_rows}")
    print(f"twisted_geometric_union_measures={geometric_rows}")
    print(f"positive_raw_arrival_triples={raw_triple_rows}")
    print(f"raw_triple_first_clock_pairs={tuple(raw_clock_pairs.items())}")
    print("twisted_three_event_legal=0: every positive raw triple has a diagonal first stored edge")
    print(f"untwisted_three_event_clock_paths={nominal_untwisted_three_paths}/{nominal_twisted_paths}: two fixed interfaces force the middle edge 0->0")
    print(f"positive_twisted_pair_types={positive_twisted}; positive_untwisted_pair_types={positive_untwisted}")
    print(f"T_rail_object_digest={object_digest}")
    print("controls=rail_reflection:%d,shallow_reflection:%d,owner_reflection:%d,T_clock_covariance:%d,T2_equals_D2:%d,T_pair_to_D_pair:%d" % (*reflection_controls, bijection_checks))
    print("verdict=NO NEW HANDOFF: untwisted T fails the state identification; the lawful J-twist is an alternating gauge copy of D and inherits THM-2684 length-three nilpotence")
    print("boundary=no present/tooth/word/unit/semantic-endpoint/scalar consequence; the clock-zero untwisted two-event subcarrier is positive but cannot iterate")


if __name__ == "__main__":
    main()
