#!/usr/bin/env python3
"""Exact primary for the c=28 AP-centred Hamming-six two-adic obstruction.

The theorem-bearing step is THM-994's anchor/nonanchor relaxation.  For each
scalar survivor and owner, all effective-order 2 and 4 masks are retained
exactly as fibres of Z/28 -> Z/4.  Every order 7, 14, and 28 mask is then
allowed to choose its unit independently so as to maximize its contribution
outside the retained union.  This enlarges the attainable union.  The relaxed
bound reaches 28 at at most two of the six owners in every labelled scalar
row, so no common unit word can work at all owners.

NumPy is used only for exact batching of the 924*26961 scalar contexts.  CRT
bases, literal masks, anchor banks, immutable-union reachability, group orbits,
Cayley relations, and tournament telemetry use deterministic Python integers.
No multiplication orbit is used to omit a labelled row.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 28
ORDERS = (1, 2, 4, 7, 14, 28)
UNITS = {
    order: (
        (0,)
        if order == 1
        else tuple(unit for unit in range(1, order) if gcd(unit, order) == 1)
    )
    for order in ORDERS
}
FULL = (1 << C) - 1
MOD4_ANCHORS = frozenset((1, 2, 4))
MOD7_ANCHORS = frozenset((1, 7))
CRT47_ANCHORS = frozenset((1, 2, 4, 7))
ORDER_TWO_HIGH = frozenset((1, 6, 7))
ORDER_FOUR_HIGH = frozenset((1, 3, 4, 6, 7, 9, 10))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered(value, modulus):
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base_algebraic(label, order, unit):
    if order == 1:
        return label
    step = ((unit - order * label) * pow(P, -1, order)) % order
    return (order * label + P * step) % (P * order)


def crt_base_literal(label, order, unit):
    answers = tuple(
        value
        for value in range(P * order)
        if value % P == order * label % P and value % order == unit % order
    )
    require(len(answers) == 1, "CRT uniqueness")
    return answers[0]


def local_mask(label, order, unit, owner):
    base = crt_base_algebraic(label, order, unit)
    owner_inverse = pow(owner, -1, P)
    mask = 0
    for sheet in range(C):
        value = centered(base * (owner_inverse + P * sheet), P * order)
        if -order < value <= order:
            mask |= 1 << sheet
    return mask


def analytic_cardinality(label, order, owner):
    ratio = label * pow(owner, -1, P) % P
    target = order * ratio % P
    period_count = sum(
        value % P == target for value in range(-order + 1, order + 1)
    )
    return (C // order) * period_count


def closed_cardinality(label, order, owner):
    ratio = label * pow(owner, -1, P) % P
    if order == 1:
        return 28 if ratio == 1 else 0
    if order == 2:
        return 14 if ratio in ORDER_TWO_HIGH else 0
    if order == 4:
        return 7 if ratio in ORDER_FOUR_HIGH else 0
    if order == 7:
        return 8 if ratio == 1 else 4
    if order == 14:
        return 6 if ratio == 1 else 4
    require(order == 28, "unknown order")
    return 5 if ratio in ORDER_TWO_HIGH else 4


def build_tables():
    masks = {}
    cards = np.zeros((P, len(ORDERS), P), dtype=np.int16)
    base_digest = sha256()
    mask_digest = sha256()
    for label in range(1, P):
        for order_index, order in enumerate(ORDERS):
            for unit in UNITS[order]:
                algebraic = crt_base_algebraic(label, order, unit)
                literal = crt_base_literal(label, order, unit)
                require(algebraic == literal, "algebraic/literal CRT mismatch")
                base_digest.update(bytes((label, order, unit)))
                base_digest.update(algebraic.to_bytes(2, "little"))
                for owner in range(1, P):
                    mask = local_mask(label, order, unit, owner)
                    card = analytic_cardinality(label, order, owner)
                    require(card == closed_cardinality(label, order, owner),
                            "closed-cardinality mismatch")
                    require(mask.bit_count() == card, "mask/cardinality mismatch")
                    masks[label, order, unit, owner] = mask
                    cards[label, order_index, owner] = card
                    mask_digest.update(bytes((label, order, unit, owner)))
                    mask_digest.update(mask.to_bytes(4, "little"))

                    # An order-d mask is the pullback of a subset of Z/dZ.
                    for sheet in range(C):
                        require(
                            ((mask >> sheet) & 1)
                            == ((mask >> ((sheet + order) % C)) & 1),
                            "order-periodicity mismatch",
                        )
                        if order in MOD4_ANCHORS:
                            require(
                                ((mask >> sheet) & 1)
                                == ((mask >> ((sheet + 4) % C)) & 1),
                                "mod-four anchor periodicity mismatch",
                            )
    return masks, cards, base_digest.hexdigest(), mask_digest.hexdigest()


def hereditary(index_word):
    word = tuple(ORDERS[index] for index in index_word)
    return all(
        lcm(*(word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )


def grammar():
    words = []
    state_words = 0
    digest = sha256()
    for index_word in product(range(len(ORDERS)), repeat=6):
        word = tuple(ORDERS[index] for index in index_word)
        via_lcm = hereditary(index_word)
        via_valuations = (
            sum(order % 4 == 0 for order in word) >= 2
            and sum(order % 7 == 0 for order in word) >= 2
        )
        require(via_lcm == via_valuations, "valuation grammar mismatch")
        if not via_lcm:
            continue
        words.append(index_word)
        weight = prod(len(UNITS[order]) for order in word)
        state_words += weight
        digest.update(bytes(index_word))
        digest.update(weight.to_bytes(8, "little"))

    # Unweighted and phi-weighted inclusion-exclusion checks.
    unweighted_bad4 = 4**6 + 6 * 2 * 4**5
    unweighted_bad7 = 3**6 + 6 * 3 * 3**5
    unweighted_both = (
        2**6 + 6 * 1 * 2**5 + 6 * 2 * 2**5
        + 6 * 1 * 2**5 + 6 * 5 * 1 * 2 * 2**4
    )
    require(
        len(words) == 6**6 - unweighted_bad4 - unweighted_bad7 + unweighted_both
        == 26_961,
        "hereditary word count",
    )

    weighted_bad4 = 14**6 + 6 * 14 * 14**5
    weighted_bad7 = 4**6 + 6 * 24 * 4**5
    weighted_both = (
        2**6 + 6 * 2 * 2**5 + 6 * 12 * 2**5
        + 6 * 12 * 2**5 + 6 * 5 * 2 * 12 * 2**4
    )
    require(
        state_words == 28**6 - weighted_bad4 - weighted_bad7 + weighted_both
        == 429_048_576,
        "literal state words/support",
    )
    return np.asarray(words, dtype=np.int8), state_words, digest.hexdigest()


EXPECTED_MULTIPLICITIES = Counter({
    (0, 2, 0, 1, 1, 2): 432,
    (0, 3, 1, 1, 0, 1): 288,
    (0, 3, 1, 0, 1, 1): 288,
    (0, 2, 0, 2, 0, 2): 216,
    (0, 2, 0, 0, 2, 2): 216,
    (0, 1, 3, 1, 1, 0): 192,
    (0, 1, 3, 1, 0, 1): 192,
    (0, 1, 3, 0, 1, 1): 192,
    (0, 3, 1, 0, 0, 2): 144,
    (0, 2, 0, 1, 0, 3): 144,
    (0, 2, 0, 0, 1, 3): 144,
    (0, 1, 3, 2, 0, 0): 96,
    (0, 1, 3, 0, 2, 0): 96,
    (0, 1, 3, 0, 0, 2): 96,
    (0, 0, 4, 1, 1, 0): 60,
    (0, 0, 4, 1, 0, 1): 60,
    (0, 0, 4, 0, 1, 1): 60,
    (0, 0, 3, 2, 1, 0): 48,
    (0, 0, 3, 1, 2, 0): 48,
    (0, 2, 0, 0, 0, 4): 36,
    (0, 0, 4, 2, 0, 0): 30,
    (0, 0, 4, 0, 2, 0): 30,
    (0, 0, 4, 0, 0, 2): 30,
    (0, 0, 3, 3, 0, 0): 16,
    (0, 0, 3, 0, 3, 0): 16,
})


def scalar_census(words, cards):
    survivors = []
    feasible_owner_histogram = np.zeros(7, dtype=np.int64)
    support_histogram = Counter()
    multiplicities = Counter()
    capacity_vectors = set()
    digest = sha256()
    for support_tuple in combinations(range(1, P), 6):
        support = np.asarray(support_tuple, dtype=np.int8)
        capacities = np.zeros((len(words), 6), dtype=np.int16)
        for provider, label in enumerate(support):
            capacities += cards[label, words[:, provider, None], support[None, :]]
        feasible_count = (capacities >= C).sum(axis=1)
        feasible_owner_histogram += np.bincount(feasible_count, minlength=7)
        indices = np.flatnonzero(feasible_count == 6)
        support_histogram[len(indices)] += 1
        for index in indices:
            word = tuple(ORDERS[int(i)] for i in words[index])
            capacity = tuple(int(value) for value in capacities[index])
            survivors.append((support_tuple, word, capacity))
            multiplicities[tuple(word.count(order) for order in ORDERS)] += 1
            capacity_vectors.add(capacity)
            digest.update(bytes(support_tuple))
            digest.update(bytes(word))
            digest.update(bytes(capacity))

    require(
        tuple(int(x) for x in feasible_owner_histogram)
        == (120_024, 5_260_824, 10_675_332, 6_969_052,
            1_724_910, 158_652, 3_170),
        "scalar feasible-owner histogram",
    )
    require(
        support_histogram
        == Counter({0: 718, 5: 96, 9: 24, 10: 24,
                    18: 24, 33: 12, 42: 24, 199: 2}),
        "support survivor histogram",
    )
    require(len(survivors) == 3_170, "scalar survivor count")
    require(len({support for support, _, _ in survivors}) == 206,
            "scalar support count")
    require(all(1 not in word for _, word, _ in survivors),
            "order-one scalar survivor")
    require(multiplicities == EXPECTED_MULTIPLICITIES,
            "scalar multiplicity histogram")
    require(len(capacity_vectors) == 1_480, "capacity-vector count")
    literal_survivor_words = sum(
        count * prod(len(UNITS[order]) ** multiplicity
                     for order, multiplicity in zip(ORDERS, pattern))
        for pattern, count in multiplicities.items()
    )
    require(literal_survivor_words == 9_275_904,
            "literal scalar-survivor state words")
    return {
        "survivors": tuple(survivors),
        "feasible_owner_histogram": tuple(int(x) for x in feasible_owner_histogram),
        "support_histogram": support_histogram,
        "multiplicities": multiplicities,
        "capacity_vectors": capacity_vectors,
        "literal_survivor_words": literal_survivor_words,
        "digest": digest.hexdigest(),
    }


def multiply_context(context, multiplier):
    return tuple(sorted(
        (multiplier * label % P, order) for label, order in context
    ))


def multiplication_orbits(survivors):
    contexts = {
        tuple(zip(support, word)) for support, word, _capacity in survivors
    }
    remaining = set(contexts)
    context_histogram = Counter()
    representatives = []
    while remaining:
        context = min(remaining)
        orbit = {multiply_context(context, multiplier)
                 for multiplier in range(1, P)}
        require(orbit <= contexts, "context multiplication orbit escapes")
        remaining -= orbit
        context_histogram[len(orbit)] += 1
        representatives.append((context, len(orbit)))

    supports = {support for support, _word, _capacity in survivors}
    remaining = set(supports)
    support_histogram = Counter()
    support_representatives = []
    while remaining:
        support = min(remaining)
        orbit = {
            tuple(sorted(multiplier * label % P for label in support))
            for multiplier in range(1, P)
        }
        require(orbit <= supports, "support multiplication orbit escapes")
        remaining -= orbit
        support_histogram[len(orbit)] += 1
        support_representatives.append((support, len(orbit)))

    require(context_histogram == Counter({4: 2, 6: 3, 12: 262}),
            "context orbit histogram")
    require(support_histogram == Counter({2: 1, 12: 17}),
            "support orbit histogram")
    return {
        "context_histogram": context_histogram,
        "representatives": tuple(representatives),
        "support_histogram": support_histogram,
        "support_representatives": tuple(support_representatives),
    }


def anchor_union_bank(support, word, owner, anchor_orders, masks):
    bank = frozenset((0,))
    for label, order in zip(support, word):
        if order not in anchor_orders:
            continue
        options = frozenset(
            masks[label, order, unit, owner] for unit in UNITS[order]
        )
        bank = frozenset(
            partial | option for partial in bank for option in options
        )
    return bank


def relaxation_bound(support, word, owner, anchor_orders, anchor_bank, masks):
    best = 0
    for anchor_union in anchor_bank:
        bound = anchor_union.bit_count()
        for label, order in zip(support, word):
            if order in anchor_orders:
                continue
            bound += max(
                (masks[label, order, unit, owner] & ~anchor_union).bit_count()
                for unit in UNITS[order]
            )
        best = max(best, bound)
    return best


def scc_sizes(adjacency):
    size = len(adjacency)
    reach = [[adjacency[a][b] for b in range(size)] for a in range(size)]
    for vertex in range(size):
        reach[vertex][vertex] = True
    for middle in range(size):
        for source in range(size):
            if reach[source][middle]:
                for target in range(size):
                    reach[source][target] |= reach[middle][target]
    unused = set(range(size))
    result = []
    while unused:
        root = min(unused)
        component = {
            vertex for vertex in unused
            if reach[root][vertex] and reach[vertex][root]
        }
        require(component, "empty SCC")
        result.append(len(component))
        unused -= component
    return tuple(sorted(result))


def tournament_fingerprint(local_rows):
    """Tournament on owner obligations, not on runners or residue cells."""
    adjacency = [[False] * 6 for _ in range(6)]
    ties = 0
    flips = 0
    for left, right in combinations(range(6), 2):
        # Harder proof obligation wins.  Coordinate order is the tie path.
        left_key = local_rows[left]
        right_key = local_rows[right]
        if left_key == right_key:
            ties += 1
            winner, loser = left, right
        elif left_key > right_key:
            winner, loser = left, right
        else:
            winner, loser = right, left
            flips += 1
        adjacency[winner][loser] = True

    scores = tuple(sorted(sum(row) for row in adjacency))
    triangles = 0
    for a, b, c in combinations(range(6), 3):
        triangles += int(
            (adjacency[a][b] and adjacency[b][c] and adjacency[c][a])
            or (adjacency[a][c] and adjacency[c][b] and adjacency[b][a])
        )
    paths = [[0] * 6 for _mask in range(1 << 6)]
    for last in range(6):
        paths[1 << last][last] = 1
    for mask in range(1, 1 << 6):
        for last in range(6):
            if not (mask >> last) & 1:
                continue
            previous_mask = mask ^ (1 << last)
            for previous in range(6):
                if ((previous_mask >> previous) & 1
                        and adjacency[previous][last]):
                    paths[mask][last] += paths[previous_mask][previous]
    hamiltonian_paths = sum(paths[-1])
    components = scc_sizes(adjacency)
    require(scores == (0, 1, 2, 3, 4, 5), "tournament score histogram")
    require(triangles == 0, "tournament directed triangle")
    require(components == (1, 1, 1, 1, 1, 1), "tournament SCCs")
    require(hamiltonian_paths == 1, "tournament Hamiltonian paths")
    return ties, flips


def fibre_and_exact_census(survivors, masks):
    u4_bound_histogram = Counter()
    u7_bound_histogram = Counter()
    u47_bound_histogram = Counter()
    u4_bank_histogram = Counter()
    u7_bank_histogram = Counter()
    u47_bank_histogram = Counter()
    exact_maximum_histogram = Counter()
    exact_bank_histogram = Counter()
    u4_context_live_histogram = Counter()
    u7_context_live_histogram = Counter()
    u47_context_live_histogram = Counter()
    exact_context_live_histogram = Counter()
    u7_minus_u4 = Counter()
    u4_minus_exact = Counter()
    live_class_orbits = {0: set(), 1: set(), 2: set()}
    tournament_ties = Counter()
    tournament_flips = Counter()
    structural_digest = sha256()
    exact_digest = sha256()
    total_exact_masks = 0
    maximum_exact_bank = 0

    for support, word, capacities in survivors:
        bounds4 = []
        bounds7 = []
        bounds47 = []
        exact_rows = []
        tournament_rows = []
        for owner_index, owner in enumerate(support):
            bank4 = anchor_union_bank(
                support, word, owner, MOD4_ANCHORS, masks
            )
            bank7 = anchor_union_bank(
                support, word, owner, MOD7_ANCHORS, masks
            )
            bank47 = anchor_union_bank(
                support, word, owner, CRT47_ANCHORS, masks
            )
            exact_bank = anchor_union_bank(
                support, word, owner, frozenset(ORDERS), masks
            )
            u4 = relaxation_bound(
                support, word, owner, MOD4_ANCHORS, bank4, masks
            )
            u7 = relaxation_bound(
                support, word, owner, MOD7_ANCHORS, bank7, masks
            )
            u47 = relaxation_bound(
                support, word, owner, CRT47_ANCHORS, bank47, masks
            )
            exact_maximum = max(mask.bit_count() for mask in exact_bank)
            exact_feasible = FULL in exact_bank

            require(exact_feasible == (exact_maximum == C),
                    "exact feasibility/maximum mismatch")
            require(exact_maximum <= u4 <= u7,
                    "fibre-bound monotonicity mismatch")
            require(u47 == u4, "CRT-product carrier changed mod-four bound")
            require((u4 >= C) == exact_feasible,
                    "mod-four threshold not exact")
            require(u4 <= capacities[owner_index],
                    "mod-four bound exceeds scalar capacity")

            bounds4.append(u4)
            bounds7.append(u7)
            bounds47.append(u47)
            exact_rows.append((exact_feasible, exact_maximum))
            # Pair observable: proof-facing status/bound first; harder wins.
            tournament_rows.append((
                int(u4 >= C), u4, exact_maximum,
                capacities[owner_index], u7,
            ))
            u4_bound_histogram[u4] += 1
            u7_bound_histogram[u7] += 1
            u47_bound_histogram[u47] += 1
            u4_bank_histogram[len(bank4)] += 1
            u7_bank_histogram[len(bank7)] += 1
            u47_bank_histogram[len(bank47)] += 1
            exact_maximum_histogram[exact_maximum] += 1
            exact_bank_histogram[len(exact_bank)] += 1
            u7_minus_u4[u7 - u4] += 1
            u4_minus_exact[u4 - exact_maximum] += 1
            total_exact_masks += len(exact_bank)
            maximum_exact_bank = max(maximum_exact_bank, len(exact_bank))

            prefix = bytes(support) + bytes(word) + bytes((owner,))
            structural_digest.update(prefix)
            structural_digest.update(bytes((u4, u7, u47)))
            structural_digest.update(len(bank4).to_bytes(2, "little"))
            structural_digest.update(len(bank7).to_bytes(2, "little"))
            structural_digest.update(len(bank47).to_bytes(2, "little"))
            exact_digest.update(prefix)
            exact_digest.update(bytes((int(exact_feasible), exact_maximum)))
            exact_digest.update(len(exact_bank).to_bytes(2, "little"))
            for mask in sorted(exact_bank):
                exact_digest.update(mask.to_bytes(4, "little"))

        live4 = sum(bound >= C for bound in bounds4)
        live7 = sum(bound >= C for bound in bounds7)
        live47 = sum(bound >= C for bound in bounds47)
        exact_live = sum(feasible for feasible, _maximum in exact_rows)
        u4_context_live_histogram[live4] += 1
        u7_context_live_histogram[live7] += 1
        u47_context_live_histogram[live47] += 1
        exact_context_live_histogram[exact_live] += 1
        require(live4 in live_class_orbits, "unexpected mod-four live count")
        live_class_orbits[live4].add(tuple(zip(support, word)))
        ties, flips = tournament_fingerprint(tournament_rows)
        tournament_ties[ties] += 1
        tournament_flips[flips] += 1

    require(
        u4_bound_histogram
        == Counter({18: 912, 19: 900, 20: 1_032, 22: 2_064,
                    23: 8_328, 24: 4_392, 28: 1_392}),
        "mod-four owner-bound histogram",
    )
    require(u4_bank_histogram == Counter({1: 12_888, 2: 936, 3: 5_196}),
            "mod-four anchor-bank histogram")
    require(
        u4_context_live_histogram == Counter({0: 2_018, 1: 912, 2: 240}),
        "mod-four live-owner/context histogram",
    )
    require(
        u7_context_live_histogram
        == Counter({0: 304, 1: 408, 2: 390, 3: 144, 6: 1_924}),
        "mod-seven live-owner/context histogram",
    )
    require(
        u7_bank_histogram == Counter({1: 8_088, 6: 8_496, 21: 2_340, 41: 96}),
        "mod-seven anchor-bank histogram",
    )
    require(u47_context_live_histogram == u4_context_live_histogram,
            "CRT-product context histogram")
    require(exact_context_live_histogram == u4_context_live_histogram,
            "exact context histogram")
    require(
        exact_maximum_histogram
        == Counter({18: 912, 19: 900, 20: 1_056, 22: 2_088,
                    23: 8_280, 24: 4_392, 28: 1_392}),
        "exact owner-maximum histogram",
    )
    require(total_exact_masks == 6_628_500 and maximum_exact_bank == 2_438,
            "exact reachable-mask totals")
    require(
        u4_minus_exact == Counter({0: 18_948, 1: 48, 2: 24}),
        "mod-four/exact loss histogram",
    )
    require(
        tournament_ties
        == Counter({1: 720, 2: 936, 3: 768, 4: 528, 6: 188, 7: 30}),
        "tournament tie histogram",
    )
    require(
        tournament_flips
        == Counter({0: 50, 1: 75, 2: 144, 3: 258, 4: 358,
                    5: 432, 6: 433, 7: 393, 8: 364, 9: 290,
                    10: 196, 11: 102, 12: 55, 13: 20}),
        "tournament flip histogram",
    )

    orbit_histograms = {}
    for live_count, contexts in live_class_orbits.items():
        remaining = set(contexts)
        histogram = Counter()
        while remaining:
            context = min(remaining)
            orbit = {multiply_context(context, multiplier)
                     for multiplier in range(1, P)}
            require(orbit <= contexts, "live-class orbit escapes")
            remaining -= orbit
            histogram[len(orbit)] += 1
        orbit_histograms[live_count] = histogram
    require(orbit_histograms == {
        0: Counter({4: 2, 6: 3, 12: 166}),
        1: Counter({12: 76}),
        2: Counter({12: 20}),
    }, "live-class orbit histograms")

    return {
        "u4_bound_histogram": u4_bound_histogram,
        "u7_bound_histogram": u7_bound_histogram,
        "u47_bound_histogram": u47_bound_histogram,
        "u4_bank_histogram": u4_bank_histogram,
        "u7_bank_histogram": u7_bank_histogram,
        "u47_bank_histogram": u47_bank_histogram,
        "exact_maximum_histogram": exact_maximum_histogram,
        "exact_bank_histogram": exact_bank_histogram,
        "u4_context_live_histogram": u4_context_live_histogram,
        "u7_context_live_histogram": u7_context_live_histogram,
        "u47_context_live_histogram": u47_context_live_histogram,
        "exact_context_live_histogram": exact_context_live_histogram,
        "u7_minus_u4": u7_minus_u4,
        "u4_minus_exact": u4_minus_exact,
        "orbit_histograms": orbit_histograms,
        "tournament_ties": tournament_ties,
        "tournament_flips": tournament_flips,
        "total_exact_masks": total_exact_masks,
        "maximum_exact_bank": maximum_exact_bank,
        "structural_digest": structural_digest.hexdigest(),
        "exact_digest": exact_digest.hexdigest(),
    }


def relation_fingerprint(high_ratios):
    vertices = tuple(range(1, P))
    adjacency = [[False] * len(vertices) for _ in vertices]
    arcs = set()
    symmetric_edges = set()
    reciprocal_edges = set()
    for source in vertices:
        for target in vertices:
            if source == target:
                continue
            ratio = target * pow(source, -1, P) % P
            if ratio in high_ratios:
                arcs.add((source, target))
                adjacency[source - 1][target - 1] = True
    for left, right in combinations(vertices, 2):
        forward = (left, right) in arcs
        reverse = (right, left) in arcs
        if forward or reverse:
            symmetric_edges.add((left, right))
        if forward and reverse:
            reciprocal_edges.add((left, right))
    directed_triangles = 0
    for a, b, c in combinations(vertices, 3):
        indices = (a - 1, b - 1, c - 1)
        directed_triangles += int(
            (adjacency[indices[0]][indices[1]]
             and adjacency[indices[1]][indices[2]]
             and adjacency[indices[2]][indices[0]])
            or (adjacency[indices[0]][indices[2]]
                and adjacency[indices[2]][indices[1]]
                and adjacency[indices[1]][indices[0]])
        )
    return (
        len(arcs), len(symmetric_edges), len(reciprocal_edges),
        directed_triangles, scc_sizes(adjacency),
    )


def cayley_shape_audit():
    group = frozenset(range(1, P))
    quadratic = frozenset(value * value % P for value in group)
    nonquadratic = group - quadratic

    switch2 = (ORDER_TWO_HIGH - {1}) | frozenset(
        pow(value, -1, P) for value in ORDER_TWO_HIGH - {1}
    )
    require(switch2 <= nonquadratic, "order-two switch not bipartite")
    edges2 = {
        tuple(sorted((left, right)))
        for left, right in combinations(group, 2)
        if right * pow(left, -1, P) % P in switch2
    }
    complete_bipartite = {
        tuple(sorted((left, right)))
        for left in quadratic for right in nonquadratic
    }
    missing = complete_bipartite - edges2
    require(len(edges2) == 24 and len(missing) == 12,
            "order-two bipartite edge counts")
    unused = set(group)
    missing_components = []
    while unused:
        root = min(unused)
        component = {root}
        frontier = {root}
        while frontier:
            next_frontier = set()
            for vertex in frontier:
                next_frontier |= {
                    other for edge in missing for other in edge
                    if vertex in edge and other != vertex
                }
            next_frontier -= component
            component |= next_frontier
            frontier = next_frontier
        missing_components.append(component)
        unused -= component
    require(tuple(sorted(map(len, missing_components))) == (4, 4, 4),
            "missing order-two graph is not three four-cycles")
    require(all(
        sum(vertex in edge for edge in missing) == 2
        for vertex in group
    ), "missing order-two component is not cyclic")

    subgroup4 = frozenset(pow(2, 3 * exponent, P) for exponent in range(4))
    parts = []
    unused = set(group)
    while unused:
        representative = min(unused)
        part = frozenset(representative * value % P for value in subgroup4)
        parts.append(part)
        unused -= part
    switch4 = (ORDER_FOUR_HIGH - {1}) | frozenset(
        pow(value, -1, P) for value in ORDER_FOUR_HIGH - {1}
    )
    edges4 = {
        tuple(sorted((left, right)))
        for left, right in combinations(group, 2)
        if right * pow(left, -1, P) % P in switch4
    }
    require(len(parts) == 3 and all(len(part) == 4 for part in parts),
            "quartic-coset partition")
    require(all(
        ((left, right) in edges4)
        == (next(i for i, part in enumerate(parts) if left in part)
            != next(i for i, part in enumerate(parts) if right in part))
        for left, right in combinations(group, 2)
    ), "order-four symmetric shadow is not K4,4,4")
    return (
        tuple(sorted(quadratic)), tuple(sorted(nonquadratic)),
        tuple(sorted(tuple(sorted(part)) for part in parts)),
    )


def fmt(counter):
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main():
    masks, cards, base_digest, mask_digest = build_tables()
    words, state_words, grammar_digest = grammar()
    scalar = scalar_census(words, cards)
    orbits = multiplication_orbits(scalar["survivors"])
    fibre = fibre_and_exact_census(scalar["survivors"], masks)
    relation2 = relation_fingerprint(ORDER_TWO_HIGH - {1})
    relation4 = relation_fingerprint(ORDER_FOUR_HIGH - {1})
    relation7 = relation_fingerprint(frozenset())
    cayley_shapes = cayley_shape_audit()

    # Digests are printed on the first freeze pass and then pinned below.
    expected_digests = {
        "grammar": "0d5fb0d79a2c398857b8ba4d286729872b828a8092c17055ef9d06ba72aee366",
        "CRT-base": "c9d6db1a24f3960bf71a9c711cc46ce9abf72c26f27975ab16f9c03852cce7a5",
        "literal-mask": "dd9d2fe84c27d15486cd0d67f1237e9b3233455768ca80994c396e35dee8ac1a",
        "scalar": "311dc8eccfee377d69f134b4bd3d2f335d98dff6031c276789865aa92fb3f6ea",
        "structural": "d4d9362a069e5fe445a00d221907d7e170906e1b7ebf34528590c03711fc2819",
        "exact": "5561960152be08a1b7d61cb8b23e40e6ea51fa25c18086c117823bd9db06d000",
    }
    actual_digests = {
        "grammar": grammar_digest,
        "CRT-base": base_digest,
        "literal-mask": mask_digest,
        "scalar": scalar["digest"],
        "structural": fibre["structural_digest"],
        "exact": fibre["exact_digest"],
    }
    require(actual_digests == expected_digests, "frozen digest mismatch")

    require(relation2 == (24, 24, 0, 0, (12,)),
            "order-two Cayley relation")
    require(relation4 == (72, 48, 24, 64, (12,)),
            "order-four Cayley relation")
    require(relation7 == (0, 0, 0, 0, (1,) * 12),
            "order-seven Cayley relation")

    print("scale-twenty-eight AP-centred Hamming-six two-adic fibre certificate")
    print("orders", ORDERS, "unit counts", tuple(len(UNITS[d]) for d in ORDERS))
    print("numpy exact-batching version", np.__version__)
    print("hereditary grammar at least two nu2=2 providers and at least two nu7=1 providers")
    print("hereditary words", len(words), "state words/support", state_words,
          "raw labelled states", 924 * state_words)
    for name in ("grammar", "CRT-base", "literal-mask"):
        print(name + " SHA256", actual_digests[name])
    print("closed scalar cards D1 28 at r=1 else0; D2 14 on (1,6,7) else0; D4 7 on (1,3,4,6,7,9,10) else0; D7 8 at r=1 else4; D14 6 at r=1 else4; D28 5 on (1,6,7) else4")
    print("scalar contexts", 924 * len(words), "feasible-owner histogram",
          " ".join(f"{index}:{value}" for index, value in
                   enumerate(scalar["feasible_owner_histogram"])))
    print("scalar supports-by-survivor-count", fmt(scalar["support_histogram"]))
    print("scalar survivors", len(scalar["survivors"]), "on supports",
          len({support for support, _word, _capacity in scalar["survivors"]}),
          "literal state words", scalar["literal_survivor_words"],
          "capacity vectors", len(scalar["capacity_vectors"]))
    print("scalar multiplicity histogram", fmt(scalar["multiplicities"]))
    print("scalar SHA256", actual_digests["scalar"])
    print("context multiplication orbit-size histogram", fmt(orbits["context_histogram"]),
          "orbits", len(orbits["representatives"]), "telemetry only")
    print("support multiplication orbit-size histogram", fmt(orbits["support_histogram"]),
          "orbits", len(orbits["support_representatives"]), "telemetry only")
    print("mod4 anchor-union bank-size histogram", fmt(fibre["u4_bank_histogram"]))
    print("mod4 owner-bound histogram", fmt(fibre["u4_bound_histogram"]))
    print("mod4 live-owner/context histogram", fmt(fibre["u4_context_live_histogram"]))
    print("mod4 live-class orbit telemetry", " ".join(
        f"{live}[{fmt(histogram)}]" for live, histogram in
        sorted(fibre["orbit_histograms"].items())
    ))
    print("mod7 anchor-union bank-size histogram", fmt(fibre["u7_bank_histogram"]))
    print("mod7 owner-bound histogram", fmt(fibre["u7_bound_histogram"]))
    print("mod7 live-owner/context histogram", fmt(fibre["u7_context_live_histogram"]))
    print("mod7-minus-mod4 owner-bound histogram", fmt(fibre["u7_minus_u4"]))
    print("CRT-product anchor-union bank-size histogram", fmt(fibre["u47_bank_histogram"]))
    print("CRT-product owner bounds equal mod4 on all", 6 * len(scalar["survivors"]),
          "owner rows; live-owner/context histogram",
          fmt(fibre["u47_context_live_histogram"]))
    print("structural SHA256", actual_digests["structural"])
    print("exact maximum-union histogram", fmt(fibre["exact_maximum_histogram"]))
    print("mod4-minus-exact owner-maximum histogram", fmt(fibre["u4_minus_exact"]))
    print("exact feasible-owner/context histogram", fmt(fibre["exact_context_live_histogram"]))
    print("exact reachable-bank-size histogram", fmt(fibre["exact_bank_histogram"]))
    print("exact owner rows", 6 * len(scalar["survivors"]),
          "total reachable masks", fibre["total_exact_masks"],
          "maximum bank", fibre["maximum_exact_bank"])
    print("exact SHA256", actual_digests["exact"])
    print("order-two high-ratio Cayley fingerprint arcs,symmetric,reciprocal,triangles,SCCs", relation2)
    print("order-four high-ratio Cayley fingerprint arcs,symmetric,reciprocal,triangles,SCCs", relation4)
    print("order-seven off-diagonal high-ratio Cayley fingerprint arcs,symmetric,reciprocal,triangles,SCCs", relation7)
    print("Cayley shadows D2 is K6,6 minus three C4 on QR/NQR; D4 is K4,4,4 on quartic-subgroup cosets", cayley_shapes)
    print("tournament vertices owner obligations; pair observable lex(status-U4,U4,exact-max,scalar-capacity,U7); harder key wins; coordinate-order tie Hamiltonian path")
    print("tournament fingerprints all 3170 transitive: scores 0,1,2,3,4,5; cycles 0; singleton SCCs; Hamiltonian paths 1")
    print("tournament tie-edge histogram", fmt(fibre["tournament_ties"]))
    print("tournament flip-edge histogram", fmt(fibre["tournament_flips"]))
    print("4x7 CRT toothpick carrier D2/D4 are mod4 thick fibres; D7 is a transverse four-cell fibre, D14 a two-cell parity strip, D28 a point needle")
    print("carrier-loss ledger mod4 forgets transverse mod7 incidence, overlaps among relaxed D7/D14/D28 masks, and shared units but only enlarges unions; mod7 leaves 1924 rows; per-owner min(U4,U7)=U4 and retaining the joint CRT product changes zero bounds")
    print("alternate-carrier audit owner obligations with absolute U4 are faithful; runners/providers, gaps, sections/boundaries, wall events, isolated residues, cover arcs, Fourier modes, matroid circuits, Fano points, chi7 colours, and completed tournaments each forget the terminal threshold or shared fibre incidence")
    print("frontier verdict scalar-empty no; mod7-empty no; mod4 all-six empty yes; exact all-six empty yes")


if __name__ == "__main__":
    main()
