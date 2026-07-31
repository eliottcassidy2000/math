#!/usr/bin/env python3
"""Exact first scout for THM-2670's physical D-pullback handoff.

This keeps the source, clock, predecessor, future, rail, half-tooth, and
binary-carry labels until after the two physical edge events have been
constructed.  For

    D(x) = {13 x},       y = {13^6 x},

an edge has the form ``B(x) 1_Q(y)`` and the next pulled-back edge has the
form ``B'(13x) 1_Q'(13y)``.  Both products are resolved exactly on the
denominator ``13*T`` grid.  The existing prefix sweep then decides the
intersection without expanding the 13^6 inverse images.

The script checks the richest two-edge clock/source cell in each guard
sector, exhausts the entire danger/danger bank, and resolves all four cospan
legs in the canonical positive cell.  It records a raw labelled atom witness
rather than inferring positivity from a Boolean matrix product.
"""

from collections import Counter, defaultdict

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old
import lrc14_guard_cospan_successor_private_clock_collapse as guard
import lrc14_successor_private_sharp_graph_clock_collapse as sharp


P = 13
Q7 = 7
T = old.T
T2 = P * T
R = P**6
INV2 = 7

EXPECTED_EVENT_SUMMARIES = {
    "safe": (6_855, 1_797_403, 407),
    "danger": (1_536, 385_929, 407),
    "guard_free": (7_679, 2_001_882, 407),
}
EXPECTED_RICHEST = {
    "safe": ((3, 1, 0, 1, 1), 36, 20, 2, 14),
    "danger": ((0, 1, 0, 1, 1), 2, 0, 2, 0),
    "guard_free": ((3, 1, 0, 1, 1), 43, 25, 7, 11),
}
EXPECTED_FIRST_WITNESSES = {
    "safe": (
        (5, 2, 6),
        166_706_319_780,
        (2, 12, 0, 0),
        (0, 12, 1, 1),
        126_816_337_986_097_204_341_478_787_325_120,
        23,
    ),
    "guard_free": (
        (12, 1, 9),
        260_931_630_960,
        (2, 12, 0, 0),
        (0, 12, 0, 1),
        198_496_045_051_817_934_627_393_332_745_600,
        12,
    ),
}
EXPECTED_COSPAN_HISTS = {
    ("safe", "safe"): Counter({(36, 20): 132, (0, 0): 12}),
    ("safe", "danger"): Counter({(5, 2): 132, (0, 0): 12}),
    ("danger", "safe"): Counter({(10, 3): 121, (6, 3): 11,
                                    (0, 0): 12}),
    ("danger", "danger"): Counter({(2, 0): 121, (0, 0): 23}),
}
EXPECTED_COSPAN_SUPPORT_SIZES = {
    ("safe", "safe"): 20,
    ("safe", "danger"): 2,
    ("danger", "safe"): 3,
    ("danger", "danger"): 0,
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def merge_intervals(intervals):
    """Return the Boolean union of half-open intervals."""
    if not intervals:
        return ()
    ordered = sorted(intervals)
    out = [list(ordered[0])]
    for left, right in ordered[1:]:
        if left <= out[-1][1]:
            if right > out[-1][1]:
                out[-1][1] = right
        else:
            out.append([left, right])
    return tuple((left, right) for left, right in out)


def intersect_sorted(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def lift_first(intervals):
    """View a T-grid set on the common T2=P*T grid."""
    return tuple((P * a, P * b) for a, b in intervals)


def pullback_D(intervals):
    """Return {x:{13x} lies in intervals} on the T2 grid."""
    return tuple(
        (sheet * T + left, sheet * T + right)
        for sheet in range(P)
        for left, right in intervals
    )


def make_prefix(intervals):
    starts = [left for left, _ in intervals]
    lens = [right - left for left, right in intervals]
    pref = [0]
    for length in lens:
        pref.append(pref[-1] + length)
    return starts, lens, pref


def phi_at(x, starts, lens, pref):
    """Distribution function of an interval union on [0,T2)."""
    from bisect import bisect_right

    i = bisect_right(starts, x) - 1
    if i < 0:
        return 0
    return pref[i] + min(x - starts[i], lens[i])


def delayed_overlap_numerator(base_intervals, delayed_prefix):
    """Numerator over R*T2 of int w_B(x)1_Q({Rx}) dx."""
    starts, lens, pref = delayed_prefix
    q_length = pref[-1]
    weighted_length = 0
    acc_r = 0
    acc_phi = 0
    rred = R % T2
    for piece in base_intervals:
        if len(piece) == 2:
            left, right = piece
            weight = 1
        else:
            left, right, weight = piece
        rleft = left * rred % T2
        rright = right * rred % T2
        weighted_length += weight * (right - left)
        acc_r += weight * (rright - rleft)
        acc_phi += weight * (
            phi_at(rright, starts, lens, pref)
            - phi_at(rleft, starts, lens, pref)
        )
    floor_numerator = R * weighted_length - acc_r
    require(floor_numerator % T2 == 0,
            "T2 prefix sweep has a nonintegral floor count")
    value = q_length * (floor_numerator // T2) + acc_phi
    require(value >= 0, "T2 prefix sweep produced negative mass")
    return value


def prefix_intervals(prefix):
    starts, lens, _ = prefix
    return tuple((left, left + length)
                 for left, length in zip(starts, lens))


def build_physical_events():
    """Rebuild every positive labelled edge occurrence before unioning."""
    (module, inherited_safe, _, _, rails,
     present, present_starts) = cross.build_carrier_data()
    prefixes, _, _, _ = guard.build_guard_cospan(module, inherited_safe)

    occurrences = {sector: defaultdict(list) for sector in guard.SECTORS}
    positive_occurrences = Counter()
    for alpha, (source, rail_clock, rail_digit, pieces) in enumerate(rails):
        for step in sharp.CLOCK_STEPS:
            future_clock = (rail_clock + step) % Q7
            for future_digit in guard.ADMISSIBLE_GRAPH_DIGITS:
                deep_digit = (-future_digit - 1) % P
                overlap = old.intersect_weighted_union(
                    pieces,
                    present[future_clock, (-future_digit) % P],
                    present_starts[future_clock, (-future_digit) % P],
                )
                for half, left, right in (
                    (0, 14 * deep_digit, 14 * deep_digit + 13),
                    (1, 14 * deep_digit - 13, 14 * deep_digit),
                ):
                    half_tooth = old.intersect_weighted_comb(
                        overlap, module.C3, 182, left, right
                    )
                    for carry in (0, 1):
                        predecessor = INV2 * (
                            deep_digit - half - carry
                        ) % P
                        predecessor_piece = old.intersect_weighted_comb(
                            half_tooth,
                            guard.PREDECESSOR_SPEED,
                            P,
                            predecessor,
                            predecessor + 1,
                        )
                        carry_piece = old.intersect_weighted_comb(
                            predecessor_piece,
                            guard.SUCCESSOR_SPEED,
                            2,
                            carry,
                            carry + 1,
                        )
                        if not carry_piece:
                            continue
                        support = merge_intervals(
                            (left0, right0)
                            for left0, right0, weight in carry_piece
                            if weight > 0
                        )
                        if not support:
                            continue
                        key = (
                            future_clock,
                            rail_clock,
                            source,
                            predecessor,
                            future_digit,
                        )
                        label = (alpha, rail_digit, half, carry)
                        for sector in guard.SECTORS:
                            value = old.delayed_weighted_numerator(
                                carry_piece,
                                prefixes[sector][future_clock][future_digit],
                            )
                            if value:
                                occurrences[sector][key].append(
                                    (label, support, tuple(carry_piece), value)
                                )
                                positive_occurrences[sector] += 1

    require(
        all(positive_occurrences[sector]
            == guard.EXPECTED_POSITIVE_ATOMS[sector]
            for sector in guard.SECTORS),
        "positive occurrence census disagrees with THM-2670",
    )

    events = {sector: {} for sector in guard.SECTORS}
    labels = {sector: {} for sector in guard.SECTORS}
    for sector in guard.SECTORS:
        for key, atoms in occurrences[sector].items():
            events[sector][key] = merge_intervals(
                interval
                for _, support, _, _ in atoms
                for interval in support
            )
            labels[sector][key] = tuple(
                label for label, _, _, _ in atoms
            )
    return prefixes, events, labels, occurrences, positive_occurrences


def build_edge_index(events, sector):
    index = defaultdict(lambda: defaultdict(list))
    for key in events[sector]:
        future, rail, source, predecessor, output = key
        index[future, rail, source][output].append((predecessor, key))
    return index


def formal_chains(index0, a, b, c, source0, source1, index1=None):
    """Return formal (j,h,k) chains for E_(b,a), E_(c,b)."""
    if index1 is None:
        index1 = index0
    first = index0.get((a, b, source0), {})
    raw_second = index1.get((b, c, source1), {})
    second = defaultdict(list)
    for output, entries in raw_second.items():
        for predecessor, key in entries:
            second[predecessor].append((output, key))
    return tuple(
        (j, h, k, key0, key1)
        for h in range(P)
        for j, key0 in first.get(h, ())
        for k, key1 in second.get(h, ())
    )


def delayed_pair_prefix(prefixes, sector0, a, h, b, k, sector1=None):
    if sector1 is None:
        sector1 = sector0
    q0 = prefix_intervals(prefixes[sector0][a][h])
    q1 = prefix_intervals(prefixes[sector1][b][k])
    joint = intersect_sorted(lift_first(q0), pullback_D(q1))
    return make_prefix(joint)


def physical_chain_value(prefixes, events, sector0, chain, sector1=None,
                         delayed_cache=None, lift_cache=None,
                         pullback_cache=None):
    if sector1 is None:
        sector1 = sector0
    if delayed_cache is None:
        delayed_cache = {}
    if lift_cache is None:
        lift_cache = {}
    if pullback_cache is None:
        pullback_cache = {}
    j, h, k, key0, key1 = chain
    a, b, _, _, _ = key0
    b1, c, _, _, _ = key1
    require(b1 == b, "formal middle clock mismatch")
    lift_key = (sector0, key0)
    pullback_key = (sector1, key1)
    if lift_key not in lift_cache:
        lift_cache[lift_key] = lift_first(events[sector0][key0])
    if pullback_key not in pullback_cache:
        pullback_cache[pullback_key] = pullback_D(events[sector1][key1])
    base = intersect_sorted(
        lift_cache[lift_key],
        pullback_cache[pullback_key],
    )
    if not base:
        return 0, 0, 0
    delayed_key = (sector0, a, h, sector1, b, k)
    if delayed_key not in delayed_cache:
        delayed_cache[delayed_key] = delayed_pair_prefix(
            prefixes, sector0, a, h, b, k, sector1
        )
    delayed = delayed_cache[delayed_key]
    if delayed[2][-1] == 0:
        return 0, len(base), 0
    value = delayed_overlap_numerator(base, delayed)
    return value, len(base), len(delayed[0])


def intersect_weighted_D(first, second):
    """Multiply two weighted base packets after pulling back the second."""
    left = tuple((P * a, P * b, w) for a, b, w in first)
    right = tuple(
        (sheet * T + a, sheet * T + b, w)
        for sheet in range(P)
        for a, b, w in second
    )
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b, left[i][2] * right[j][2]))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def atom_pair_witnesses(prefixes, occurrences, sector0, chain, sector1=None):
    """Return every positive labelled atom pair beneath one union chain."""
    if sector1 is None:
        sector1 = sector0
    _, h, k, key0, key1 = chain
    a, b, _, _, _ = key0
    witnesses = []
    delayed = delayed_pair_prefix(
        prefixes, sector0, a, h, b, k, sector1
    )
    for label0, _, weighted0, _ in occurrences[sector0][key0]:
        for label1, _, weighted1, _ in occurrences[sector1][key1]:
            base = intersect_weighted_D(weighted0, weighted1)
            if not base:
                continue
            value = delayed_overlap_numerator(base, delayed)
            if value:
                witnesses.append((label0, label1, value, len(base)))
    return tuple(witnesses)


def richest_formal_cell(index):
    best = None
    histogram = Counter()
    for a in range(Q7):
        for b in range(Q7):
            if a == b:
                continue
            for c in range(Q7):
                if b == c:
                    continue
                for source0 in range(1, P):
                    for source1 in range(1, P):
                        chains = formal_chains(
                            index, a, b, c, source0, source1
                        )
                        histogram[len(chains)] += 1
                        candidate = (
                            len(chains), -a, -b, -c,
                            -source0, -source1,
                        )
                        if best is None or candidate > best[0]:
                            best = (
                                candidate,
                                (a, b, c, source0, source1, chains),
                            )
    return best[1], histogram


def exhaustive_sector_census(prefixes, events, occurrences, sector):
    """Exhaust every formal two-edge chain in one sector."""
    index = build_edge_index(events, sector)
    totals = Counter()
    state_support = set()
    first_witness = None
    for a in range(Q7):
        for b in range(Q7):
            if a == b:
                continue
            for c in range(Q7):
                if b == c:
                    continue
                for source0 in range(1, P):
                    for source1 in range(1, P):
                        chains = formal_chains(
                            index, a, b, c, source0, source1
                        )
                        totals["formal"] += len(chains)
                        for chain in chains:
                            value, base_parts, _ = physical_chain_value(
                                prefixes, events, sector, chain
                            )
                            if not value:
                                totals[
                                    "base_zero" if base_parts == 0
                                    else "skew_zero"
                                ] += 1
                                continue
                            totals["physical"] += 1
                            state_support.add(chain[:3])
                            if first_witness is None:
                                atom_witnesses = atom_pair_witnesses(
                                    prefixes, occurrences, sector, chain
                                )
                                require(atom_witnesses,
                                        "positive union chain has no atom witness")
                                first_witness = (
                                    (a, b, c, source0, source1),
                                    chain[:3],
                                    atom_witnesses[0],
                                )
    require(
        totals["formal"]
        == totals["physical"] + totals["base_zero"] + totals["skew_zero"],
        "exhaustive physical census lost a formal chain",
    )
    return totals, tuple(sorted(state_support)), first_witness


def fixed_cell_cross_census(prefixes, events, occurrences,
                            sector0, sector1, cell=(3, 1, 0)):
    """Exhaust both source labels in one ordered clock triple."""
    a, b, c = cell
    index0 = build_edge_index(events, sector0)
    index1 = build_edge_index(events, sector1)
    delayed_cache = {}
    lift_cache = {}
    pullback_cache = {}
    by_source = {}
    support = set()
    witness = None
    for source0 in range(1, P):
        for source1 in range(1, P):
            chains = formal_chains(
                index0, a, b, c, source0, source1, index1
            )
            physical = []
            for chain in chains:
                value, _, _ = physical_chain_value(
                    prefixes, events, sector0, chain, sector1,
                    delayed_cache, lift_cache, pullback_cache,
                )
                if value:
                    physical.append(chain)
                    support.add(chain[:3])
                    if witness is None:
                        atom_witnesses = atom_pair_witnesses(
                            prefixes, occurrences,
                            sector0, chain, sector1,
                        )
                        require(atom_witnesses,
                                "cross-sector union has no atom witness")
                        witness = (
                            (source0, source1), chain[:3],
                            atom_witnesses[0],
                        )
            by_source[source0, source1] = (
                len(chains), len(physical)
            )
    return Counter(by_source.values()), tuple(
        pair for pair, counts in sorted(by_source.items()) if counts[1] == 0
    ), tuple(sorted(support)), witness


def main():
    prefixes, events, labels, occurrences, positive_occurrences = (
        build_physical_events()
    )
    print("THM-2670 physical D-pullback first scout")
    print(f"T={T} T2={T2} R={R}")
    print(f"positive_occurrences={tuple(sorted(positive_occurrences.items()))}")
    richest_rows = {}
    for sector in guard.SECTORS:
        index = build_edge_index(events, sector)
        interval_counts = tuple(map(len, events[sector].values()))
        summary = (
            len(events[sector]), sum(interval_counts), max(interval_counts)
        )
        require(summary == EXPECTED_EVENT_SUMMARIES[sector],
                f"{sector} event support summary changed")
        print(
            f"sector={sector} event_summary="
            f"(keys,total_intervals,max_intervals){summary}"
        )
        (a, b, c, source0, source1, chains), histogram = (
            richest_formal_cell(index)
        )
        positive = []
        base_zeros = delayed_zeros = 0
        for chain in chains:
            value, base_parts, delayed_parts = physical_chain_value(
                prefixes, events, sector, chain
            )
            if value:
                positive.append((chain[:3], value, base_parts, delayed_parts))
            elif base_parts == 0:
                base_zeros += 1
            else:
                delayed_zeros += 1
        richest = (
            (a, b, c, source0, source1), len(chains), len(positive),
            base_zeros, delayed_zeros,
        )
        require(richest == EXPECTED_RICHEST[sector],
                f"{sector} richest physical cell changed")
        require(sum(count * multiplicity
                    for count, multiplicity in histogram.items()) > 0,
                f"{sector} formal atlas unexpectedly empty")
        richest_rows[sector] = richest
        print(f"sector={sector} richest={richest}")
        if positive:
            state = positive[0][0]
            selected = next(chain for chain in chains if chain[:3] == state)
            witnesses = atom_pair_witnesses(
                prefixes, occurrences, sector, selected
            )
            require(witnesses,
                    "displayed physical chain has no labelled atom witness")
            displayed = (
                state, positive[0][1], witnesses[0][0], witnesses[0][1],
                witnesses[0][2], witnesses[0][3],
            )
            require(displayed == EXPECTED_FIRST_WITNESSES[sector],
                    f"{sector} first labelled witness changed")
            print(f"sector={sector} first_atom_witness={displayed}")
    danger_totals, danger_states, danger_witness = exhaustive_sector_census(
        prefixes, events, occurrences, "danger"
    )
    print(f"danger_exhaustive_totals={tuple(sorted(danger_totals.items()))}")
    print(f"danger_exhaustive_state_support={danger_states}")
    print(f"danger_exhaustive_first_witness={danger_witness}")
    require(danger_totals == Counter({"formal": 27_640,
                                      "base_zero": 27_640}),
            "danger/danger exhaustive zero atlas changed")
    require(not danger_states and danger_witness is None,
            "danger/danger acquired a physical state")
    for sector0, sector1 in (
        ("safe", "safe"),
        ("safe", "danger"),
        ("danger", "safe"),
        ("danger", "danger"),
    ):
        histogram, zero_sources, support, witness = fixed_cell_cross_census(
            prefixes, events, occurrences, sector0, sector1
        )
        leg = (sector0, sector1)
        require(histogram == EXPECTED_COSPAN_HISTS[leg],
                f"{leg} source histogram changed")
        require(len(support) == EXPECTED_COSPAN_SUPPORT_SIZES[leg],
                f"{leg} state support changed")
        if leg != ("danger", "danger"):
            require(zero_sources == tuple((source0, 6)
                                          for source0 in range(1, P)),
                    f"{leg} inactive source column changed")
            require(witness is not None, f"{leg} lost its atom witness")
        else:
            require(len(zero_sources) == 144 and witness is None,
                    "danger/danger zero-source atlas changed")
        print(
            f"fixed_cell_310_legs={sector0}->{sector1} "
            f"source_hist={tuple(sorted(histogram.items()))} "
            f"zero_source_count={len(zero_sources)} "
            f"state_support_size={len(support)} witness={witness}"
        )
    print("verdict=PASS: positive atomwise D-pullback exists; danger-square is physically zero")
    print("scope=two-edge cospan; global safe atlas and three-edge chronology audited separately")


if __name__ == "__main__":
    main()
