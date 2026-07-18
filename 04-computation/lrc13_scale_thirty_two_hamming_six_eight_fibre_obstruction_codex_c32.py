#!/usr/bin/env python3
"""Exact primary for the c=32 AP-centred Hamming-six eight-fibre obstruction.

The theorem-bearing step is THM-994's anchor/nonanchor relaxation.  After the
complete hereditary and scalar scans, every mask of effective order dividing
8 is retained exactly as a fibre of Z/32 -> Z/8.  Each order-16 or order-32
provider is then allowed to choose its unit independently so as to maximize
its contribution outside the retained union.  This forgets shared units and
nonanchor overlaps and therefore enlarges attainable coverage.  The relaxed
bound reaches 32 at no more than two of the six owners in every scalar row.

NumPy is used only to batch the 924*12281 exact scalar contexts.  CRT bases,
literal masks, fibre banks, covariance, and tournament telemetry use Python
integers.  No orbit quotient omits a labelled row.

Tournament Analysis deliberately uses owner obligations, not runners, as
vertices.  The pairwise observable is the lexicographic proof cost
(relaxed bound, scalar capacity); coordinate order is the tie Hamiltonian
path.  This quotient preserves which owner-local obligations remain hardest,
but destroys the shared global unit word and all nonanchor overlaps.  It is
telemetry only; soundness comes from the relaxation inequality.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 32
ORDERS = (1, 2, 4, 8, 16, 32)
UNITS = {
    order: (
        (0,)
        if order == 1
        else tuple(unit for unit in range(1, order) if gcd(unit, order) == 1)
    )
    for order in ORDERS
}
ANCHORS = frozenset((1, 2, 4, 8))
FULL = (1 << C) - 1

EXPECTED_SCALAR_FEASIBLE = (
    76_548, 2_800_212, 4_692_582, 2_946_408, 743_040, 85_404, 3_450
)
EXPECTED_SUPPORT_HISTOGRAM = Counter({
    0: 640, 1: 96, 2: 48, 3: 30, 17: 24,
    22: 24, 33: 12, 34: 24, 38: 24, 54: 2,
})
EXPECTED_BANK_HISTOGRAM = Counter({
    1: 8_064, 2: 4_008, 3: 636, 4: 3_576, 6: 2_088, 7: 240,
    8: 192, 9: 36, 10: 684, 12: 984, 14: 168, 26: 24,
})
EXPECTED_BOUND_HISTOGRAM = Counter({
    20: 24, 21: 396, 22: 420, 23: 564, 24: 1_416, 25: 2_352,
    26: 3_252, 27: 3_780, 28: 3_708, 29: 2_208, 30: 1_260,
    31: 480, 32: 624, 33: 192, 34: 24,
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


def build_tables():
    choices = {}
    cards = np.zeros((P, len(ORDERS), P), dtype=np.int16)
    base_digest = sha256()
    mask_digest = sha256()
    cardinality_rows = {}

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
                options = []
                for unit in UNITS[order]:
                    mask = local_mask(label, order, unit, owner)
                    require(mask.bit_count() == card, "mask/cardinality mismatch")
                    options.append(mask)
                    mask_digest.update(bytes((label, order, unit, owner)))
                    mask_digest.update(mask.to_bytes(4, "little"))

                    # Every order-D mask is pulled back from Z/DZ.  If D|8,
                    # it is consequently a sound Z/8 anchor.
                    for sheet in range(C):
                        require(
                            ((mask >> sheet) & 1)
                            == ((mask >> ((sheet + order) % C)) & 1),
                            "order-periodicity mismatch",
                        )
                        if order in ANCHORS:
                            require(
                                ((mask >> sheet) & 1)
                                == ((mask >> ((sheet + 8) % C)) & 1),
                                "eight-fibre periodicity mismatch",
                            )
                choices[label, order, owner] = tuple(sorted(set(options)))
                cards[label, order_index, owner] = card

        # Owner normalization makes the cardinality depend only on a/b.
        cardinality_rows[label] = tuple(
            analytic_cardinality(label, order, 1) for order in ORDERS
        )

    for label in range(1, P):
        for owner in range(1, P):
            ratio = label * pow(owner, -1, P) % P
            require(
                tuple(int(cards[label, i, owner]) for i in range(len(ORDERS)))
                == cardinality_rows[ratio],
                "ratio covariance of scalar table",
            )

    return {
        "choices": choices,
        "cards": cards,
        "base_digest": base_digest.hexdigest(),
        "mask_digest": mask_digest.hexdigest(),
        "cardinality_rows": cardinality_rows,
    }


def hereditary_lcm(index_word) -> bool:
    word = tuple(ORDERS[index] for index in index_word)
    return all(
        lcm(*(word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )


def build_grammar():
    words = []
    state_words = 0
    digest = sha256()
    for index_word in product(range(len(ORDERS)), repeat=6):
        via_lcm = hereditary_lcm(index_word)
        via_valuation = sum(ORDERS[index] == 32 for index in index_word) >= 2
        require(via_lcm == via_valuation, "2-adic grammar mismatch")
        if not via_lcm:
            continue
        word = tuple(ORDERS[index] for index in index_word)
        weight = prod(len(UNITS[order]) for order in word)
        words.append(index_word)
        state_words += weight
        digest.update(bytes(index_word))
        digest.update(weight.to_bytes(8, "little"))

    require(
        len(words) == 6**6 - 5**6 - 6 * 5**5 == 12_281,
        "hereditary word count",
    )
    require(sum(len(UNITS[order]) for order in ORDERS[:-1]) == 16,
            "proper-divisor unit mass")
    require(len(UNITS[32]) == 16, "top-order unit mass")
    require(
        state_words == 32**6 - 7 * 16**6 == 956_301_312,
        "literal state words/support",
    )
    return {
        "words": np.asarray(words, dtype=np.int8),
        "state_words": state_words,
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
            pattern = tuple(word.count(order) for order in ORDERS)
            multiplicities[pattern] += 1
            capacity_vectors.add(capacity)
            digest.update(bytes(support_tuple))
            digest.update(bytes(word))
            digest.update(bytes(capacity))

    require(len(survivors) == 3_450, "scalar survivor count")
    require(len({support for support, _word, _cap in survivors}) == 284,
            "scalar survivor support count")
    literal_survivor_words = sum(
        count * prod(
            len(UNITS[order]) ** multiplicity
            for order, multiplicity in zip(ORDERS, pattern)
        )
        for pattern, count in multiplicities.items()
    )
    require(tuple(int(value) for value in feasible_owner_histogram)
            == EXPECTED_SCALAR_FEASIBLE,
            "scalar feasible-owner histogram")
    require(support_histogram == EXPECTED_SUPPORT_HISTOGRAM,
            "scalar support histogram")
    require(len(multiplicities) == 23, "scalar multiplicity-profile count")
    require(len(capacity_vectors) == 1_649, "scalar capacity-vector count")
    require(literal_survivor_words == 621_084_672,
            "literal scalar-survivor state words")
    require(digest.hexdigest()
            == "aba8c8c9d589103d83c135f512bec5375ded7997adc166c77b9051cc2b2c106a",
            "scalar digest")
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


def anchor_union_bank(support, word, owner, choices):
    bank = frozenset((0,))
    for label, order in zip(support, word):
        if order not in ANCHORS:
            continue
        bank = frozenset(
            partial | option
            for partial in bank
            for option in choices[label, order, owner]
        )
    return bank


def relaxation_bound(support, word, owner, anchor_bank, choices):
    best = 0
    for anchor_union in anchor_bank:
        outside = FULL ^ anchor_union
        bound = anchor_union.bit_count()
        for label, order in zip(support, word):
            if order in ANCHORS:
                continue
            bound += max(
                (mask & outside).bit_count()
                for mask in choices[label, order, owner]
            )
        best = max(best, bound)
    return best


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
        components.append(len(component))
        unused -= component
    return tuple(sorted(components))


def tournament_fingerprint(bounds, capacities):
    """Cost tournament on owner obligations; coordinate order breaks ties."""
    adjacency = [[False] * 6 for _ in range(6)]
    ties = 0
    flips = 0
    for left, right in combinations(range(6), 2):
        left_key = (bounds[left], capacities[left])
        right_key = (bounds[right], capacities[right])
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
            or (adjacency[b][a] and adjacency[c][b] and adjacency[a][c])
        )
    components = scc_sizes(adjacency)
    require(scores == (0, 1, 2, 3, 4, 5), "cost tournament not transitive")
    require(triangles == 0 and components == (1, 1, 1, 1, 1, 1),
            "cost tournament cycle/SCC mismatch")
    return ties, flips, scores, triangles, components


def flag_census(survivors, choices):
    bound_histogram = Counter()
    bank_size_histogram = Counter()
    live_owner_histogram = Counter()
    tournament_histogram = Counter()
    digest = sha256()

    for support, word, capacities in survivors:
        bounds = []
        for owner in support:
            bank = anchor_union_bank(support, word, owner, choices)
            bound = relaxation_bound(support, word, owner, bank, choices)
            bank_size_histogram[len(bank)] += 1
            bound_histogram[bound] += 1
            bounds.append(bound)
            digest.update(bytes(support))
            digest.update(bytes(word))
            digest.update(bytes((owner, bound)))
            digest.update(len(bank).to_bytes(4, "little"))
        live_owner_histogram[sum(bound >= C for bound in bounds)] += 1
        tournament_histogram[tournament_fingerprint(bounds, capacities)] += 1

    require(live_owner_histogram == Counter({0: 2_802, 1: 456, 2: 192}),
            "eight-fibre live-owner histogram")
    require(max(live_owner_histogram) == 2, "three-owner flag survivor")
    require(bank_size_histogram == EXPECTED_BANK_HISTOGRAM,
            "anchor-bank histogram")
    require(bound_histogram == EXPECTED_BOUND_HISTOGRAM,
            "relaxed owner-bound histogram")
    require(sum(tournament_histogram.values()) == 3_450,
            "tournament context count")
    require(digest.hexdigest()
            == "d6c6f54aa1efbfc0968a1bc11583bab9de59ef589fc509164d3940ae67c145e1",
            "flag digest")
    return {
        "bound_histogram": bound_histogram,
        "bank_size_histogram": bank_size_histogram,
        "live_owner_histogram": live_owner_histogram,
        "tournament_histogram": tournament_histogram,
        "digest": digest.hexdigest(),
    }


def multiply_context(context, multiplier):
    return tuple(sorted(
        (multiplier * label % P, order) for label, order in context
    ))


def orbit_census(survivors):
    contexts = {
        tuple(zip(support, word)) for support, word, _capacity in survivors
    }
    remaining = set(contexts)
    context_histogram = Counter()
    while remaining:
        context = min(remaining)
        orbit = {
            multiply_context(context, multiplier) for multiplier in range(1, P)
        }
        require(orbit <= contexts, "multiplication orbit escapes scalar bank")
        context_histogram[len(orbit)] += 1
        remaining -= orbit

    supports = {support for support, _word, _capacity in survivors}
    remaining_supports = set(supports)
    support_orbit_histogram = Counter()
    while remaining_supports:
        support = min(remaining_supports)
        orbit = {
            tuple(sorted(multiplier * label % P for label in support))
            for multiplier in range(1, P)
        }
        require(orbit <= supports, "support orbit escapes scalar bank")
        support_orbit_histogram[len(orbit)] += 1
        remaining_supports -= orbit
    return context_histogram, support_orbit_histogram


def main() -> None:
    tables = build_tables()
    grammar = build_grammar()
    scalar = scalar_census(grammar["words"], tables["cards"])
    flag = flag_census(scalar["survivors"], tables["choices"])
    context_orbits, support_orbits = orbit_census(scalar["survivors"])
    require(context_orbits == Counter({6: 3, 12: 286}),
            "context orbit histogram")
    require(support_orbits == Counter({2: 1, 6: 1, 12: 23}),
            "support orbit histogram")
    require(tables["base_digest"]
            == "68b88f331468047800d59975eaabd0bdf02af4567bc70ffe8a4afb0984a12b74",
            "CRT-base digest")
    require(tables["mask_digest"]
            == "f715f217aae1be03edc31d51a4eb4380024e815c67cc73983bd321beff8c5b81",
            "literal-mask digest")

    print("### c=32 AP-centred H6 exact eight-fibre primary ###")
    print("scope: primitive proper common-scale H6 owner-local gate only")
    print("orders:", ORDERS)
    print("unit counts:", tuple(len(UNITS[order]) for order in ORDERS))
    print("grammar equivalence: leave-one-out lcm=32 iff at least two D=32")
    print("hereditary order words:", len(grammar["words"]))
    print("literal state words/support:", grammar["state_words"])
    print("raw labelled state words:", 924 * grammar["state_words"])
    print("grammar digest:", grammar["digest"])
    print("CRT-base digest:", tables["base_digest"])
    print("literal-mask digest:", tables["mask_digest"])
    print("cardinality vectors by ratio r=1..12; columns D=1,2,4,8,16,32:")
    for ratio in range(1, P):
        print(f"  r={ratio:2d}: {tables['cardinality_rows'][ratio]}")
    print("scalar feasible-owner histogram:",
          " ".join(f"{i}:{v}" for i, v in enumerate(scalar["feasible_owner_histogram"])))
    print("scalar support survivor histogram:", fmt(scalar["support_histogram"]))
    print("scalar survivors:", len(scalar["survivors"]))
    print("scalar survivor supports:",
          len({support for support, _word, _cap in scalar["survivors"]}))
    print("scalar multiplicity profiles:", len(scalar["multiplicities"]))
    print("scalar capacity vectors:", len(scalar["capacity_vectors"]))
    print("literal scalar-survivor state words:", scalar["literal_survivor_words"])
    print("scalar digest:", scalar["digest"])
    print("Z/8 anchor orders:", tuple(sorted(ANCHORS)))
    print("anchor-bank size histogram:", fmt(flag["bank_size_histogram"]))
    print("relaxed owner-bound histogram:", fmt(flag["bound_histogram"]))
    print("live owners/context:", fmt(flag["live_owner_histogram"]))
    print("flag digest:", flag["digest"])
    print("context multiplication orbit sizes:", fmt(context_orbits))
    print("support multiplication orbit sizes:", fmt(support_orbits))
    print("tournament vertices: six owner obligations (not runners or sheets)")
    print("tournament observable: (relaxed bound, scalar capacity)")
    print("tournament switch/tie path: harder key wins; coordinate order breaks ties")
    print("tournament fingerprint classes:", len(flag["tournament_histogram"]))
    print("tournament tie/flip histogram:", fmt(Counter(
        (key[0], key[1]) for key, count in flag["tournament_histogram"].items()
        for _ in range(count)
    )))
    print("tournament scores: always (0,1,2,3,4,5); triangles=0; SCCs=1^6")
    print("challenged assumption: proof-obligation vertices replace runner vertices")
    print("preserved: owner-local impossibility under a sound upper relaxation")
    print("destroyed: shared unit word and all nonanchor overlaps (telemetry only)")
    print("THEOREM: every scalar row has at least four owners with U_8<32")
    print("VERDICT: primitive proper AP-centred common-scale-32 H6 face EMPTY")


if __name__ == "__main__":
    main()
