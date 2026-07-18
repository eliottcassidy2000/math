#!/usr/bin/env python3
"""Exact primary for the c=30 AP-centred Hamming-six two-flag obstruction.

The proof-facing certificate uses THM-994's sound anchor/nonanchor relaxation.
It retains every effective-order divisor of 6 exactly on the first pass.  The
complete labelled scalar bank then has only 120 rows whose six owner bounds
still reach 30.  On those rows it instead retains every effective-order
divisor of 10; every resulting owner bound is strictly below 30.

NumPy is used only to batch the 924*185193 scalar order contexts.  CRT bases,
literal masks, quotient-anchor banks, exact residual union banks, covariance,
orbits, Cayley relations, and tournaments use deterministic Python integers.
No multiplication orbit is used to omit a labelled row.  The exact immutable-
union DP is a sharpness sidecar on the complete 120-row structural residual;
the theorem-bearing implication is the two sound upper relaxations.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 30
ORDERS = (1, 2, 3, 5, 6, 10, 15, 30)
UNITS = {
    order: (
        (0,)
        if order == 1
        else tuple(unit for unit in range(1, order) if gcd(unit, order) == 1)
    )
    for order in ORDERS
}
FULL = (1 << C) - 1
FLAGS = (6, 10, 15)

HIGH_RATIOS = {
    1: frozenset((1,)),
    2: frozenset((1, 6, 7)),
    3: frozenset((1, 4, 5, 8, 9)),
    5: frozenset((1, 2, 3, 5, 6, 7, 8, 10, 11)),
    6: frozenset(range(1, 12)),
    10: frozenset((1, 2, 3, 6, 7, 10, 11)),
    15: frozenset((1, 6, 7)),
    30: frozenset((1, 3, 4, 6, 7, 9, 10)),
}

EXPECTED_SUPPORT_HISTOGRAM = Counter({
    0: 152,
    2: 112,
    4: 24,
    8: 24,
    12: 120,
    14: 6,
    17: 96,
    22: 96,
    30: 24,
    54: 96,
    68: 24,
    122: 24,
    137: 24,
    208: 24,
    253: 24,
    294: 12,
    296: 12,
    545: 24,
    549: 6,
})

EXPECTED_RESIDUAL_MULTIPLICITIES = Counter({
    (0, 0, 0, 0, 0, 2, 3, 1): 24,
    (0, 0, 0, 0, 0, 4, 0, 2): 48,
    (0, 0, 0, 1, 0, 3, 0, 2): 48,
})


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fmt(counter: Counter) -> str:
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base_algebraic(label: int, order: int, unit: int) -> int:
    if order == 1:
        return label
    step = ((unit - order * label) * pow(P, -1, order)) % order
    return (order * label + P * step) % (P * order)


def crt_base_literal(label: int, order: int, unit: int) -> int:
    answers = tuple(
        value
        for value in range(P * order)
        if value % P == order * label % P and value % order == unit % order
    )
    require(len(answers) == 1, "CRT uniqueness")
    return answers[0]


def local_mask(label: int, order: int, unit: int, owner: int) -> int:
    base = crt_base_algebraic(label, order, unit)
    owner_inverse = pow(owner, -1, P)
    mask = 0
    for sheet in range(C):
        value = centered(base * (owner_inverse + P * sheet), P * order)
        if -order < value <= order:
            mask |= 1 << sheet
    return mask


def analytic_cardinality(label: int, order: int, owner: int) -> int:
    ratio = label * pow(owner, -1, P) % P
    target = order * ratio % P
    period_count = sum(
        value % P == target for value in range(-order + 1, order + 1)
    )
    return (C // order) * period_count


def closed_cardinality(label: int, order: int, owner: int) -> int:
    ratio = label * pow(owner, -1, P) % P
    if order == 1:
        return 30 if ratio == 1 else 0
    if order == 2:
        return 15 if ratio in HIGH_RATIOS[2] else 0
    if order == 3:
        return 10 if ratio in HIGH_RATIOS[3] else 0
    if order == 5:
        return 6 if ratio in HIGH_RATIOS[5] else 0
    if order == 6:
        return 5 if ratio in HIGH_RATIOS[6] else 0
    if order == 10:
        return 6 if ratio in HIGH_RATIOS[10] else 3
    if order == 15:
        return 6 if ratio in HIGH_RATIOS[15] else 4
    require(order == 30, "unknown order")
    return 5 if ratio in HIGH_RATIOS[30] else 4


def build_tables():
    masks = {}
    choices = {}
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
                card = analytic_cardinality(label, order, owner)
                require(card == closed_cardinality(label, order, owner),
                        "closed-cardinality mismatch")
                options = []
                for unit in UNITS[order]:
                    mask = local_mask(label, order, unit, owner)
                    require(mask.bit_count() == card,
                            "mask/cardinality mismatch")
                    masks[label, order, unit, owner] = mask
                    options.append(mask)
                    mask_digest.update(bytes((label, order, unit, owner)))
                    mask_digest.update(mask.to_bytes(4, "little"))

                    # An order-D mask is a pullback from Z/DZ.  Consequently
                    # it is also a valid q-anchor whenever D divides q.
                    for sheet in range(C):
                        require(
                            ((mask >> sheet) & 1)
                            == ((mask >> ((sheet + order) % C)) & 1),
                            "order-periodicity mismatch",
                        )
                        for flag in FLAGS:
                            if flag % order == 0:
                                require(
                                    ((mask >> sheet) & 1)
                                    == ((mask >> ((sheet + flag) % C)) & 1),
                                    "flag-periodicity mismatch",
                                )
                cards[label, order_index, owner] = card
                choices[label, order, owner] = tuple(sorted(set(options)))

    return {
        "masks": masks,
        "choices": choices,
        "cards": cards,
        "base_digest": base_digest.hexdigest(),
        "mask_digest": mask_digest.hexdigest(),
    }


def hereditary_lcm(index_word) -> bool:
    word = tuple(ORDERS[index] for index in index_word)
    return all(
        lcm(*(word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )


def hereditary_squarefree(index_word) -> bool:
    word = tuple(ORDERS[index] for index in index_word)
    return all(sum(order % prime == 0 for order in word) >= 2
               for prime in (2, 3, 5))


def build_grammar():
    words = []
    state_words = 0
    digest = sha256()
    for index_word in product(range(len(ORDERS)), repeat=6):
        via_lcm = hereditary_lcm(index_word)
        via_squarefree = hereditary_squarefree(index_word)
        require(via_lcm == via_squarefree, "squarefree grammar mismatch")
        if not via_lcm:
            continue
        word = tuple(ORDERS[index] for index in index_word)
        weight = prod(len(UNITS[order]) for order in word)
        words.append(index_word)
        state_words += weight
        digest.update(bytes(index_word))
        digest.update(weight.to_bytes(8, "little"))

    # For each prime, choose independently a coordinate subset of size >=2.
    require(len(words) == (2**6 - 1 - 6) ** 3 == 185_193,
            "hereditary word count")
    weighted_factors = tuple(
        prime**6 - 1 - 6 * (prime - 1) for prime in (2, 3, 5)
    )
    require(weighted_factors == (57, 716, 15_600),
            "weighted squarefree factors")
    require(state_words == prod(weighted_factors) == 636_667_200,
            "literal state words/support")
    return {
        "words": np.asarray(words, dtype=np.int8),
        "state_words": state_words,
        "weighted_factors": weighted_factors,
        "digest": digest.hexdigest(),
    }


def scalar_census(words: np.ndarray, cards: np.ndarray):
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
            capacities += cards[
                label, words[:, provider, None], support[None, :]
            ]
        feasible_count = (capacities >= C).sum(axis=1)
        feasible_owner_histogram += np.bincount(feasible_count, minlength=7)
        indices = np.flatnonzero(feasible_count == 6)
        support_histogram[len(indices)] += 1
        for row_index in indices:
            word = tuple(ORDERS[int(index)] for index in words[row_index])
            capacity = tuple(int(value) for value in capacities[row_index])
            survivors.append((support_tuple, word, capacity))
            multiplicities[tuple(word.count(order) for order in ORDERS)] += 1
            capacity_vectors.add(capacity)
            digest.update(bytes(support_tuple))
            digest.update(bytes(word))
            for value in capacity:
                digest.update(value.to_bytes(2, "little"))

    require(
        tuple(int(value) for value in feasible_owner_histogram)
        == (1_401_966, 36_143_640, 66_874_158, 49_478_260,
            15_326_622, 1_839_636, 54_050),
        "scalar feasible-owner histogram",
    )
    require(support_histogram == EXPECTED_SUPPORT_HISTOGRAM,
            "scalar support histogram")
    require(len(survivors) == 54_050, "scalar survivor count")
    require(len({support for support, _word, _cap in survivors}) == 772,
            "scalar survivor support count")
    require(len(multiplicities) == 244, "scalar multiplicity count")
    require(len(capacity_vectors) == 28_965, "capacity-vector count")

    literal_survivor_words = sum(
        count * prod(
            len(UNITS[order]) ** multiplicity
            for order, multiplicity in zip(ORDERS, pattern)
        )
        for pattern, count in multiplicities.items()
    )
    require(literal_survivor_words == 64_678_912,
            "literal scalar-survivor words")
    return {
        "survivors": tuple(survivors),
        "feasible_owner_histogram": tuple(
            int(value) for value in feasible_owner_histogram
        ),
        "support_histogram": support_histogram,
        "multiplicities": multiplicities,
        "capacity_vectors": capacity_vectors,
        "literal_survivor_words": literal_survivor_words,
        "digest": digest.hexdigest(),
    }


def context_of(support, word):
    return tuple(zip(support, word))


def multiply_context(context, multiplier):
    return tuple(sorted(
        (multiplier * label % P, order) for label, order in context
    ))


def anchor_union_bank(support, word, owner, flag, choices):
    bank = frozenset((0,))
    for label, order in zip(support, word):
        if flag % order != 0:
            continue
        bank = frozenset(
            partial | option
            for partial in bank
            for option in choices[label, order, owner]
        )
    return bank


def relaxation_bound(support, word, owner, flag, anchor_bank, choices):
    best = 0
    for anchor_union in anchor_bank:
        outside = FULL ^ anchor_union
        bound = anchor_union.bit_count()
        for label, order in zip(support, word):
            if flag % order == 0:
                continue
            bound += max(
                (mask & outside).bit_count()
                for mask in choices[label, order, owner]
            )
        best = max(best, bound)
    return best


def exact_union_bank(support, word, owner, choices):
    bank = frozenset((0,))
    for label, order in zip(support, word):
        bank = frozenset(
            partial | option
            for partial in bank
            for option in choices[label, order, owner]
        )
    return bank


def scc_sizes(adjacency):
    size = len(adjacency)
    reach = [[bool(adjacency[a][b]) for b in range(size)] for a in range(size)]
    for vertex in range(size):
        reach[vertex][vertex] = True
    for middle in range(size):
        for source in range(size):
            if reach[source][middle]:
                for target in range(size):
                    reach[source][target] |= reach[middle][target]
    unused = set(range(size))
    components = []
    while unused:
        root = min(unused)
        component = {
            vertex for vertex in unused
            if reach[root][vertex] and reach[vertex][root]
        }
        require(component, "empty SCC")
        components.append(len(component))
        unused -= component
    return tuple(sorted(components))


def tournament_fingerprint(keys):
    size = len(keys)
    adjacency = [[False] * size for _ in range(size)]
    ties = 0
    flips = 0
    for left, right in combinations(range(size), 2):
        if keys[left] == keys[right]:
            ties += 1
            winner, loser = left, right
        elif keys[left] > keys[right]:
            winner, loser = left, right
        else:
            winner, loser = right, left
            flips += 1
        adjacency[winner][loser] = True

    scores = tuple(sorted(sum(row) for row in adjacency))
    triangles = 0
    for a, b, c in combinations(range(size), 3):
        triangles += int(
            (adjacency[a][b] and adjacency[b][c] and adjacency[c][a])
            or (adjacency[a][c] and adjacency[c][b] and adjacency[b][a])
        )
    paths = [[0] * size for _mask in range(1 << size)]
    for last in range(size):
        paths[1 << last][last] = 1
    for mask in range(1, 1 << size):
        for last in range(size):
            if not ((mask >> last) & 1):
                continue
            previous_mask = mask ^ (1 << last)
            for previous in range(size):
                if ((previous_mask >> previous) & 1
                        and adjacency[previous][last]):
                    paths[mask][last] += paths[previous_mask][previous]
    hamiltonian_paths = sum(paths[-1])
    return (
        scores,
        triangles,
        scc_sizes(adjacency),
        hamiltonian_paths,
        ties,
        flips,
    )


def flag_tournament(values):
    """Smaller quotient upper bound wins; flag order is the tie path."""
    adjacency = [[False] * 3 for _ in range(3)]
    ties = 0
    flips = 0
    for left, right in combinations(range(3), 2):
        if values[left] == values[right]:
            ties += 1
        if values[left] <= values[right]:
            winner, loser = left, right
        else:
            winner, loser = right, left
            flips += 1
        adjacency[winner][loser] = True
    scores = tuple(sorted(sum(row) for row in adjacency))
    cycle = int(
        (adjacency[0][1] and adjacency[1][2] and adjacency[2][0])
        or (adjacency[0][2] and adjacency[2][1] and adjacency[1][0])
    )
    return scores, cycle, ties, flips


def two_flag_census(survivors, choices):
    u6_bound_histogram = Counter()
    u6_bank_histogram = Counter()
    u6_live_histogram = Counter()
    residual = []
    q6_digest = sha256()

    # First pass: exact D|6 anchors on every labelled scalar survivor.
    for support, word, capacities in survivors:
        bounds = []
        for owner_index, owner in enumerate(support):
            bank = anchor_union_bank(support, word, owner, 6, choices)
            raw_bound = relaxation_bound(
                support, word, owner, 6, bank, choices
            )
            require(raw_bound <= capacities[owner_index],
                    "U6 exceeds scalar capacity")
            bound = min(C, raw_bound)
            bounds.append(bound)
            u6_bound_histogram[bound] += 1
            u6_bank_histogram[len(bank)] += 1
            q6_digest.update(bytes(support))
            q6_digest.update(bytes(word))
            q6_digest.update(bytes((owner, bound)))
            q6_digest.update(len(bank).to_bytes(2, "little"))
        live = sum(bound >= C for bound in bounds)
        u6_live_histogram[live] += 1
        if live == 6:
            residual.append((support, word, capacities, tuple(bounds)))

    require(
        u6_bound_histogram
        == Counter({16: 144, 17: 144, 18: 384, 19: 840,
                    20: 7_632, 21: 9_264, 22: 26_496, 23: 60_972,
                    24: 104_388, 25: 31_368, 26: 42_480,
                    27: 22_464, 28: 6_600, 29: 300, 30: 10_824}),
        "U6 owner-bound histogram",
    )
    require(
        u6_live_histogram
        == Counter({0: 45_110, 1: 7_536, 2: 1_284, 6: 120}),
        "U6 live-owner/context histogram",
    )
    require(
        u6_bank_histogram
        == Counter({1: 114_588, 2: 123_024, 3: 42_996,
                    4: 35_760, 5: 5_244, 7: 2_688}),
        "U6 anchor-bank histogram",
    )
    require(len(residual) == 120, "U6 residual size")

    residual_multiplicities = Counter(
        tuple(word.count(order) for order in ORDERS)
        for _support, word, _capacities, _bounds in residual
    )
    require(residual_multiplicities == EXPECTED_RESIDUAL_MULTIPLICITIES,
            "U6 residual multiplicities")
    require(all(all(order in (5, 10, 15, 30) for order in word)
                for _support, word, _capacities, _bounds in residual),
            "U6 residual contains a divisor-of-six anchor")
    residual_literal_words = sum(
        count * prod(
            len(UNITS[order]) ** multiplicity
            for order, multiplicity in zip(ORDERS, pattern)
        )
        for pattern, count in residual_multiplicities.items()
    )
    require(residual_literal_words == 3_145_728,
            "U6 residual literal words")

    u10_bound_histogram = Counter()
    u15_bound_histogram = Counter()
    u10_bank_histogram = Counter()
    u15_bank_histogram = Counter()
    u10_live_histogram = Counter()
    u15_live_histogram = Counter()
    exact_maximum_histogram = Counter()
    exact_bank_histogram = Counter()
    exact_live_histogram = Counter()
    u10_minus_exact = Counter()
    joint_contingency = Counter()
    profile_u10 = {}
    owner_ties = Counter()
    owner_flips = Counter()
    flag_scores = Counter()
    flag_cycles = 0
    flag_ties = Counter()
    flag_flips = Counter()
    q10_le_u15 = Counter()
    residual_digest = sha256()
    flag_digest = sha256()
    exact_digest = sha256()
    total_exact_masks = 0
    maximum_exact_bank = 0
    owner_summaries = {}

    for support, word, capacities, bounds6 in residual:
        pattern = tuple(word.count(order) for order in ORDERS)
        profile_u10.setdefault(pattern, Counter())
        bounds10 = []
        bounds15 = []
        exact_rows = []
        tournament_keys = []
        context = context_of(support, word)
        residual_digest.update(bytes(support))
        residual_digest.update(bytes(word))
        for value in capacities:
            residual_digest.update(value.to_bytes(2, "little"))

        for owner_index, owner in enumerate(support):
            bank10 = anchor_union_bank(support, word, owner, 10, choices)
            bank15 = anchor_union_bank(support, word, owner, 15, choices)
            raw10 = relaxation_bound(
                support, word, owner, 10, bank10, choices
            )
            raw15 = relaxation_bound(
                support, word, owner, 15, bank15, choices
            )
            exact_bank = exact_union_bank(
                support, word, owner, choices
            )
            exact_maximum = max(mask.bit_count() for mask in exact_bank)
            exact_feasible = FULL in exact_bank
            require(exact_feasible == (exact_maximum == C),
                    "exact residual feasibility/maximum mismatch")
            require(exact_maximum <= raw10 <= capacities[owner_index],
                    "U10 soundness mismatch")
            require(exact_maximum <= raw15 <= capacities[owner_index],
                    "U15 soundness mismatch")

            bound10 = min(C, raw10)
            bound15 = min(C, raw15)
            require(bound10 < C, "U10 residual owner reaches threshold")
            require(bound10 <= bound15,
                    "U10 does not dominate U15 on residual")
            bounds10.append(bound10)
            bounds15.append(bound15)
            exact_rows.append((exact_feasible, exact_maximum))
            profile_u10[pattern][bound10] += 1
            u10_bound_histogram[bound10] += 1
            u15_bound_histogram[bound15] += 1
            u10_bank_histogram[len(bank10)] += 1
            u15_bank_histogram[len(bank15)] += 1
            exact_maximum_histogram[exact_maximum] += 1
            exact_bank_histogram[len(exact_bank)] += 1
            u10_minus_exact[bound10 - exact_maximum] += 1
            joint_contingency[
                (bounds6[owner_index], bound10, bound15, exact_maximum)
            ] += 1
            q10_le_u15[(bound10 < bound15)] += 1
            total_exact_masks += len(exact_bank)
            maximum_exact_bank = max(maximum_exact_bank, len(exact_bank))

            # Harder owner proof obligation wins.  Support coordinate order
            # supplies the tie Hamiltonian path.
            tournament_keys.append((
                int(bound10 >= C), bound10, exact_maximum,
                capacities[owner_index], bound15, len(exact_bank),
            ))
            score, cycle, ties, flips = flag_tournament(
                (bounds6[owner_index], bound10, bound15)
            )
            flag_scores[score] += 1
            flag_cycles += cycle
            flag_ties[ties] += 1
            flag_flips[flips] += 1

            owner_summaries[context, owner] = (
                capacities[owner_index],
                bounds6[owner_index], bound10, bound15,
                exact_maximum, int(exact_feasible),
                len(bank10), len(bank15), len(exact_bank),
            )
            prefix = bytes(support) + bytes(word) + bytes((owner,))
            flag_digest.update(prefix)
            flag_digest.update(bytes((bounds6[owner_index], bound10, bound15)))
            flag_digest.update(len(bank10).to_bytes(2, "little"))
            flag_digest.update(len(bank15).to_bytes(2, "little"))
            exact_digest.update(prefix)
            exact_digest.update(bytes((int(exact_feasible), exact_maximum)))
            exact_digest.update(len(exact_bank).to_bytes(2, "little"))
            for mask in sorted(exact_bank):
                exact_digest.update(mask.to_bytes(4, "little"))

        u10_live_histogram[sum(bound >= C for bound in bounds10)] += 1
        u15_live_histogram[sum(bound >= C for bound in bounds15)] += 1
        exact_live_histogram[
            sum(feasible for feasible, _maximum in exact_rows)
        ] += 1
        fingerprint = tournament_fingerprint(tournament_keys)
        require(
            fingerprint[:4]
            == ((0, 1, 2, 3, 4, 5), 0, (1, 1, 1, 1, 1, 1), 1),
            "owner tournament fingerprint",
        )
        owner_ties[fingerprint[4]] += 1
        owner_flips[fingerprint[5]] += 1

    require(
        u10_bound_histogram
        == Counter({23: 48, 24: 120, 25: 192,
                    26: 240, 27: 72, 28: 48}),
        "U10 residual bound histogram",
    )
    require(u10_live_histogram == Counter({0: 120}),
            "U10 residual live-owner histogram")
    require(
        u10_bank_histogram
        == Counter({10: 144, 34: 48, 52: 96,
                    70: 96, 86: 48, 94: 288}),
        "U10 anchor-bank histogram",
    )
    require(
        u15_bound_histogram == Counter({27: 48, 28: 144, 29: 144, 30: 384}),
        "U15 residual bound histogram",
    )
    require(u15_live_histogram == Counter({0: 24, 2: 48, 6: 48}),
            "U15 residual live-owner histogram")
    require(
        u15_bank_histogram
        == Counter({1: 336, 4: 240, 132: 24,
                    172: 48, 204: 48, 260: 24}),
        "U15 anchor-bank histogram",
    )
    require(
        exact_maximum_histogram
        == Counter({23: 72, 24: 144, 25: 192, 26: 264, 27: 48}),
        "exact residual maximum histogram",
    )
    require(
        exact_bank_histogram
        == Counter({1608: 48, 2234: 48, 2528: 48, 3680: 48,
                    3772: 96, 3820: 96, 3856: 96, 4264: 96,
                    7116: 24, 8572: 48, 10784: 48, 12936: 24}),
        "exact residual bank histogram",
    )
    require(exact_live_histogram == Counter({0: 120}),
            "exact residual live-owner histogram")
    require(total_exact_masks == 3_401_088 and maximum_exact_bank == 12_936,
            "exact residual reachable-mask totals")
    require(u10_minus_exact == Counter({0: 576, 1: 48, 2: 96}),
            "U10/exact loss histogram")
    require(q10_le_u15 == Counter({True: 672, False: 48}),
            "U10/U15 strictness histogram")
    require(flag_scores == Counter({(0, 1, 2): 720}) and flag_cycles == 0,
            "flag tournament score/cycle fingerprint")
    require(flag_ties == Counter({0: 288, 1: 432}),
            "flag tournament tie histogram")
    require(flag_flips == Counter({1: 384, 2: 336}),
            "flag tournament flip histogram")
    require(
        owner_ties == Counter({1: 96, 2: 24}),
        "owner tournament tie histogram",
    )
    require(
        owner_flips
        == Counter({1: 1, 2: 5, 3: 6, 4: 11, 5: 17, 6: 18,
                    7: 9, 8: 15, 9: 20, 10: 6, 11: 6, 12: 5, 13: 1}),
        "owner tournament flip histogram",
    )

    require(
        joint_contingency
        == Counter({
            (30, 23, 30, 23): 48,
            (30, 24, 28, 23): 24,
            (30, 24, 29, 24): 48,
            (30, 24, 30, 24): 48,
            (30, 25, 28, 25): 48,
            (30, 25, 29, 25): 96,
            (30, 25, 30, 25): 48,
            (30, 26, 27, 24): 48,
            (30, 26, 30, 26): 192,
            (30, 27, 28, 26): 24,
            (30, 27, 30, 27): 48,
            (30, 28, 28, 26): 48,
        }),
        "joint flag/exact contingency",
    )

    expected_profile_u10 = {
        (0, 0, 0, 0, 0, 2, 3, 1):
            Counter({24: 24, 26: 48, 27: 24, 28: 48}),
        (0, 0, 0, 0, 0, 4, 0, 2):
            Counter({23: 48, 24: 48, 25: 48, 26: 144}),
        (0, 0, 0, 1, 0, 3, 0, 2):
            Counter({24: 48, 25: 144, 26: 48, 27: 48}),
    }
    require(profile_u10 == expected_profile_u10,
            "profile-wise U10 histogram")

    return {
        "residual": tuple(residual),
        "residual_multiplicities": residual_multiplicities,
        "residual_literal_words": residual_literal_words,
        "u6_bound_histogram": u6_bound_histogram,
        "u6_bank_histogram": u6_bank_histogram,
        "u6_live_histogram": u6_live_histogram,
        "u10_bound_histogram": u10_bound_histogram,
        "u15_bound_histogram": u15_bound_histogram,
        "u10_bank_histogram": u10_bank_histogram,
        "u15_bank_histogram": u15_bank_histogram,
        "u10_live_histogram": u10_live_histogram,
        "u15_live_histogram": u15_live_histogram,
        "exact_maximum_histogram": exact_maximum_histogram,
        "exact_bank_histogram": exact_bank_histogram,
        "exact_live_histogram": exact_live_histogram,
        "u10_minus_exact": u10_minus_exact,
        "joint_contingency": joint_contingency,
        "profile_u10": profile_u10,
        "owner_ties": owner_ties,
        "owner_flips": owner_flips,
        "flag_scores": flag_scores,
        "flag_ties": flag_ties,
        "flag_flips": flag_flips,
        "q10_le_u15": q10_le_u15,
        "total_exact_masks": total_exact_masks,
        "maximum_exact_bank": maximum_exact_bank,
        "owner_summaries": owner_summaries,
        "q6_digest": q6_digest.hexdigest(),
        "residual_digest": residual_digest.hexdigest(),
        "flag_digest": flag_digest.hexdigest(),
        "exact_digest": exact_digest.hexdigest(),
    }


def orbit_and_covariance_census(survivors, two_flag):
    contexts = {
        context_of(support, word): tuple(zip(support, capacities))
        for support, word, capacities in survivors
    }
    remaining = set(contexts)
    context_orbits = Counter()
    while remaining:
        context = min(remaining)
        orbit = {
            multiply_context(context, multiplier)
            for multiplier in range(1, P)
        }
        require(orbit <= contexts.keys(), "scalar context orbit escapes")
        remaining -= orbit
        context_orbits[len(orbit)] += 1

    # Check labelled scalar capacities, not merely membership, under every
    # multiplier.  This orbit action is telemetry and deletes no primary row.
    for context, capacity_pairs in contexts.items():
        capacity_map = dict(capacity_pairs)
        for multiplier in range(1, P):
            target = multiply_context(context, multiplier)
            target_capacity = dict(contexts[target])
            for owner, capacity in capacity_map.items():
                require(
                    target_capacity[multiplier * owner % P] == capacity,
                    "scalar capacity covariance mismatch",
                )

    supports = {tuple(label for label, _order in context) for context in contexts}
    remaining_supports = set(supports)
    support_orbits = Counter()
    while remaining_supports:
        support = min(remaining_supports)
        orbit = {
            tuple(sorted(multiplier * label % P for label in support))
            for multiplier in range(1, P)
        }
        require(orbit <= supports, "scalar support orbit escapes")
        remaining_supports -= orbit
        support_orbits[len(orbit)] += 1

    residual_contexts = {
        context_of(support, word)
        for support, word, _capacities, _bounds in two_flag["residual"]
    }
    remaining_residual = set(residual_contexts)
    residual_orbits = Counter()
    residual_representatives = []
    while remaining_residual:
        context = min(remaining_residual)
        orbit = {
            multiply_context(context, multiplier)
            for multiplier in range(1, P)
        }
        require(orbit <= residual_contexts, "residual context orbit escapes")
        remaining_residual -= orbit
        residual_orbits[len(orbit)] += 1
        residual_representatives.append((min(orbit), len(orbit)))
    require(residual_orbits == Counter({12: 10}),
            "residual context orbit histogram")

    normalized_owner_keys = set()
    summaries = two_flag["owner_summaries"]
    for context in residual_contexts:
        owners = tuple(label for label, _order in context)
        for owner in owners:
            normalized_owner_keys.add(tuple(sorted(
                (label * pow(owner, -1, P) % P, order)
                for label, order in context
            )))
        for multiplier in range(1, P):
            target = multiply_context(context, multiplier)
            for owner in owners:
                require(
                    summaries[target, multiplier * owner % P]
                    == summaries[context, owner],
                    "residual owner-summary covariance mismatch",
                )
    require(len(normalized_owner_keys) == 60,
            "owner-normalized residual key count")
    require(context_orbits == Counter({4: 2, 6: 9, 12: 4_499}),
            "scalar context orbit histogram")
    require(support_orbits == Counter({4: 1, 6: 2, 12: 63}),
            "scalar support orbit histogram")

    return {
        "context_orbits": context_orbits,
        "support_orbits": support_orbits,
        "residual_orbits": residual_orbits,
        "residual_representatives": tuple(residual_representatives),
        "normalized_owner_keys": normalized_owner_keys,
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
    triangles = 0
    for a, b, c in combinations(vertices, 3):
        x, y, z = a - 1, b - 1, c - 1
        triangles += int(
            (adjacency[x][y] and adjacency[y][z] and adjacency[z][x])
            or (adjacency[x][z] and adjacency[z][y] and adjacency[y][x])
        )
    return (
        len(arcs), len(symmetric_edges), len(reciprocal_edges),
        triangles, scc_sizes(adjacency),
    )


def main() -> None:
    tables = build_tables()
    grammar = build_grammar()
    scalar = scalar_census(grammar["words"], tables["cards"])
    two_flag = two_flag_census(scalar["survivors"], tables["choices"])
    orbit = orbit_and_covariance_census(scalar["survivors"], two_flag)

    expected_digests = {
        "grammar": "d5ec0f6ac4c8f74ed399fbce0e8b50f232789f82619f556f7606b2e41fb6d55f",
        "bases": "c94c62de38c0238f686e34cfc77293359e8638418485e59e49bed509f992dc3e",
        "masks": "dcff92077f5e320927a98c7f5a86df40e9c6baed2075e7cad46542d40d3e8da9",
        "scalar": "6be8cfd5291e955df7cffd5112471bc0d2e232e816cc68d25de3031034e4dd93",
        "u6": "bf5be40ba8f998c194319593abbddfbe277a15df3805c2ebff0f24d3b81e3342",
        "residual": "e8c698b8d62f804c67925977c4dbc397b26af00a5908b553fb5fb4f0eb10d59d",
        "flags": "c1607cc6959e5bd1aea3bbcf63353dbee6eadce5b30709ee98864c108ca3ea32",
        "exact": "955669e80f11fc3239ea076416966fc4f27ad6fb5e9dc9453f1e692979c959d8",
    }
    require(grammar["digest"] == expected_digests["grammar"],
            "grammar stream digest")
    require(tables["base_digest"] == expected_digests["bases"],
            "CRT-base stream digest")
    require(tables["mask_digest"] == expected_digests["masks"],
            "literal-mask stream digest")
    require(scalar["digest"] == expected_digests["scalar"],
            "scalar-bank stream digest")
    require(two_flag["q6_digest"] == expected_digests["u6"],
            "U6 stream digest")
    require(two_flag["residual_digest"] == expected_digests["residual"],
            "residual stream digest")
    require(two_flag["flag_digest"] == expected_digests["flags"],
            "flag stream digest")
    require(two_flag["exact_digest"] == expected_digests["exact"],
            "exact-bank stream digest")

    print("scale-thirty AP-centred Hamming-six two-flag primary")
    print("scope primitive proper common-scale H6 owner-local gate only")
    print("divisor grammar", ",".join(map(str, ORDERS)),
          "; unit counts", ",".join(str(len(UNITS[d])) for d in ORDERS))
    print("supports 924; hereditary order words", len(grammar["words"]),
          "; labelled order contexts", 924 * len(grammar["words"]))
    print("squarefree weighted factors", grammar["weighted_factors"])
    print("state words/support", grammar["state_words"],
          "; raw labelled states", 924 * grammar["state_words"])
    print("numpy exact-batching version", np.__version__)
    print("scalar feasible-owner histogram", " ".join(
        f"{index}:{value}"
        for index, value in enumerate(scalar["feasible_owner_histogram"])
    ))
    print("scalar contexts", len(scalar["survivors"]), "on",
          len({s for s, _w, _c in scalar["survivors"]}),
          "supports; multiplicity profiles", len(scalar["multiplicities"]),
          "; capacity vectors", len(scalar["capacity_vectors"]))
    print("scalar support survivor histogram", fmt(scalar["support_histogram"]))
    print("literal scalar-survivor unit words", scalar["literal_survivor_words"])
    print("U6 saturated owner-bound histogram", fmt(two_flag["u6_bound_histogram"]))
    print("U6 anchor-bank histogram", fmt(two_flag["u6_bank_histogram"]))
    print("U6 live owners/context", fmt(two_flag["u6_live_histogram"]))
    print("U6 residual contexts", len(two_flag["residual"]),
          "; literal unit words", two_flag["residual_literal_words"])
    print("U6 residual multiplicities", fmt(two_flag["residual_multiplicities"]))
    print("U10 residual owner-bound histogram", fmt(two_flag["u10_bound_histogram"]))
    print("U10 anchor-bank histogram", fmt(two_flag["u10_bank_histogram"]))
    print("U10 live owners/context", fmt(two_flag["u10_live_histogram"]))
    print("U15 competing owner-bound histogram", fmt(two_flag["u15_bound_histogram"]))
    print("U15 anchor-bank histogram", fmt(two_flag["u15_bank_histogram"]))
    print("U15 live owners/context", fmt(two_flag["u15_live_histogram"]))
    print("U10 strictly beats U15", fmt(two_flag["q10_le_u15"]))
    print("profile-wise U10 histograms")
    for pattern in sorted(two_flag["profile_u10"]):
        print(" ", pattern, fmt(two_flag["profile_u10"][pattern]),
              "maximum", max(two_flag["profile_u10"][pattern]))
    print("exact residual maximum histogram", fmt(two_flag["exact_maximum_histogram"]))
    print("exact residual live owners/context", fmt(two_flag["exact_live_histogram"]))
    print("exact residual bank histogram", fmt(two_flag["exact_bank_histogram"]))
    print("exact residual reachable masks", two_flag["total_exact_masks"],
          "; largest bank", two_flag["maximum_exact_bank"])
    print("U10 minus exact maximum", fmt(two_flag["u10_minus_exact"]))
    print("joint U6,U10,U15,exact contingency")
    for key in sorted(two_flag["joint_contingency"]):
        print(" ", key, two_flag["joint_contingency"][key])
    print("scalar context orbit histogram", fmt(orbit["context_orbits"]),
          "(telemetry; no quotient)")
    print("scalar support orbit histogram", fmt(orbit["support_orbits"]),
          "(telemetry; no quotient)")
    print("residual context orbit histogram", fmt(orbit["residual_orbits"]),
          "; owner-normalized keys", len(orbit["normalized_owner_keys"]))
    print("residual orbit representatives")
    for representative, size in orbit["residual_representatives"]:
        print(" ", " ".join(f"{label}:{order}" for label, order in representative),
              "orbit", size)
    print("owner-obligation tournament observable",
          "(U10-live,U10,exact,capacity,U15,bank); harder wins;",
          "support-coordinate tie Hamiltonian path")
    print("owner tournament all 120 transitive: scores 0,1,2,3,4,5;",
          "cycles 0; SCCs 6; Hamiltonian paths 1")
    print("owner tournament tie histogram", fmt(two_flag["owner_ties"]))
    print("owner tournament edge-flip histogram", fmt(two_flag["owner_flips"]))
    print("flag tournament vertices 6,10,15; smaller bound wins;",
          "flag-order tie Hamiltonian path")
    print("flag tournament scores", fmt(two_flag["flag_scores"]),
          "; cycles 0; tie histogram", fmt(two_flag["flag_ties"]),
          "; edge-flip histogram", fmt(two_flag["flag_flips"]))
    for order in (2, 3, 10, 30):
        print("Cayley high-cardinality relation D", order,
              "arcs,edges,reciprocal,triangles,SCCs",
              relation_fingerprint(HIGH_RATIOS[order] - {1}))
    print("internal SHA256 grammar", grammar["digest"])
    print("internal SHA256 CRT bases", tables["base_digest"])
    print("internal SHA256 literal masks", tables["mask_digest"])
    print("internal SHA256 scalar bank", scalar["digest"])
    print("internal SHA256 U6 rows", two_flag["q6_digest"])
    print("internal SHA256 residual rows", two_flag["residual_digest"])
    print("internal SHA256 U10/U15 flags", two_flag["flag_digest"])
    print("internal SHA256 exact residual banks", two_flag["exact_digest"])
    print("carrier loss U6/U10 retain full masks for D|q; nonanchors choose",
          "independent best units outside the anchor union, forgetting their",
          "shared units and pair overlaps; this only enlarges coverage")
    print("challenged vertices owner obligations and q-fibres preserve the",
          "absolute terminal bound; runners/providers lose fibre overlap;",
          "the completed tournaments lose the threshold and are telemetry")
    print("primary verdict scale 30 structurally empty conditional on",
          "independent replay; theorem status remains CLAIMED + FROZEN PRIMARY")


if __name__ == "__main__":
    main()
