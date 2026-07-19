#!/usr/bin/env python3
"""Exact c=36 AP-centred Hamming-six complementary-fibre certificate.

Every labelled six-support and every hereditary effective-order word is
visited.  Scalar owner capacity leaves a finite residual.  On that residual
we use two sound upper relaxations.  The first retains every provider whose
order divides four, and independently maximizes all other masks outside the
exact Z/4 anchor union.  The second does the same for orders dividing nine.
No scalar row is live at all six owners in both gauges.

Tournament Analysis uses the six owner obligations as vertices.  The
observable is the owner-local tuple (live flag, relaxed ceiling, scalar
capacity, anchor-bank size), with coordinate order as the tie Hamiltonian
path.  The tournament is telemetry; the proof-bearing predicate is the
absolute failed-owner ceiling in at least one of the two gauges.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 36
ORDERS = (1, 2, 3, 4, 6, 9, 12, 18, 36)
UNITS = {
    order: (
        (0,)
        if order == 1
        else tuple(unit for unit in range(1, order) if gcd(unit, order) == 1)
    )
    for order in ORDERS
}
FULL = (1 << C) - 1
FLAGS = {
    "Z4": frozenset(order for order in ORDERS if 4 % order == 0),
    "Z9": frozenset(order for order in ORDERS if 9 % order == 0),
}

EXPECTED_DIGESTS = {
    "grammar": "e95823dec1c339d191803fdb8c1571a61a0ef10d5c0239bcf3c4f63627f94119",
    "base": "6525f0e6a79291883867d8992fc79e1c971381a4f98010baa73d60edd957d81b",
    "mask": "a7de515af8812c5c85505feb373bfece7f8986b661bc4328c9e80fe05a3199b0",
    "scalar": "8ef7886a1e87d4ffdcf45fd3a57d3b1c20d55581ea9a4ad50441da759e05a3ad",
    "Z4": "42fc4b8e19246037e83d1b806d50e22b50e680a6b4f9bef0b743195db3e9ab5d",
    "Z9": "604f8092b9c7717bbaa8f0080046d129b4d2b490e1ec7940ecfd9c5aa0d557b7",
    "joint": "de941b205a284043c3511f1205516a4224c114aa86950e4c1751fc0ffee599f8",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fmt(counter: Counter) -> str:
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label: int, order: int, unit: int) -> int:
    if order == 1:
        return label
    step = ((unit - order * label) * pow(P, -1, order)) % order
    return (order * label + P * step) % (P * order)


def literal_crt_base(label: int, order: int, unit: int) -> int:
    answers = tuple(
        value
        for value in range(P * order)
        if value % P == order * label % P and value % order == unit % order
    )
    require(len(answers) == 1, "literal CRT representative is not unique")
    return answers[0]


def local_mask(label: int, order: int, unit: int, owner: int) -> int:
    base = crt_base(label, order, unit)
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
    hits = sum(value % P == target for value in range(-order + 1, order + 1))
    return (C // order) * hits


def build_tables():
    cards = np.zeros((P, len(ORDERS), P), dtype=np.int16)
    choices = {}
    rows = {}
    option_sizes = Counter()
    base_digest = sha256()
    mask_digest = sha256()

    for label in range(1, P):
        for order_index, order in enumerate(ORDERS):
            for unit in UNITS[order]:
                algebraic = crt_base(label, order, unit)
                literal = literal_crt_base(label, order, unit)
                require(algebraic == literal, "algebraic/literal CRT mismatch")
                base_digest.update(bytes((label, order, unit)))
                base_digest.update(algebraic.to_bytes(2, "little"))

            for owner in range(1, P):
                card = analytic_cardinality(label, order, owner)
                exact_masks = tuple(
                    local_mask(label, order, unit, owner)
                    for unit in UNITS[order]
                )
                masks = tuple(sorted(set(exact_masks)))
                require(all(mask.bit_count() == card for mask in masks),
                        "mask/cardinality mismatch")
                option_sizes[(order, len(masks))] += 1
                for unit, mask in zip(UNITS[order], exact_masks):
                    mask_digest.update(bytes((label, order, unit, owner)))
                    mask_digest.update(mask.to_bytes(8, "little"))
                    for sheet in range(C):
                        require(
                            ((mask >> sheet) & 1)
                            == ((mask >> ((sheet + order) % C)) & 1),
                            "effective-order periodicity mismatch",
                        )
                        if order in FLAGS["Z4"]:
                            require(
                                ((mask >> sheet) & 1)
                                == ((mask >> ((sheet + 4) % C)) & 1),
                                "Z/4 anchor periodicity mismatch",
                            )
                        if order in FLAGS["Z9"]:
                            require(
                                ((mask >> sheet) & 1)
                                == ((mask >> ((sheet + 9) % C)) & 1),
                                "Z/9 anchor periodicity mismatch",
                            )
                cards[label, order_index, owner] = card
                choices[label, order, owner] = masks

        rows[label] = tuple(
            analytic_cardinality(label, order, 1) for order in ORDERS
        )

    for label in range(1, P):
        for owner in range(1, P):
            ratio = label * pow(owner, -1, P) % P
            require(
                tuple(int(cards[label, index, owner])
                      for index in range(len(ORDERS))) == rows[ratio],
                "ratio covariance mismatch",
            )

    return {
        "cards": cards,
        "choices": choices,
        "rows": rows,
        "option_sizes": option_sizes,
        "base_digest": base_digest.hexdigest(),
        "mask_digest": mask_digest.hexdigest(),
    }


def hereditary(index_word: tuple[int, ...]) -> bool:
    word = tuple(ORDERS[index] for index in index_word)
    via_lcm = all(
        lcm(*(word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )
    via_carriers = (
        sum(order % 4 == 0 for order in word) >= 2
        and sum(order % 9 == 0 for order in word) >= 2
    )
    require(via_lcm == via_carriers, "lcm/carrier grammar mismatch")
    return via_lcm


def build_grammar():
    words = []
    weighted_states = 0
    carrier_histogram = Counter()
    digest = sha256()
    for index_word in product(range(len(ORDERS)), repeat=6):
        if not hereditary(index_word):
            continue
        word = tuple(ORDERS[index] for index in index_word)
        weight = prod(len(UNITS[order]) for order in word)
        words.append(index_word)
        weighted_states += weight
        carrier_histogram[
            (sum(order % 4 == 0 for order in word),
             sum(order % 9 == 0 for order in word))
        ] += 1
        digest.update(bytes(index_word))
        digest.update(weight.to_bytes(8, "little"))
    return {
        "words": np.asarray(words, dtype=np.int8),
        "weighted_states": weighted_states,
        "carrier_histogram": carrier_histogram,
        "digest": digest.hexdigest(),
    }


def scalar_census(words: np.ndarray, cards: np.ndarray):
    histogram = np.zeros(7, dtype=np.int64)
    support_histogram = Counter()
    survivors = []
    profile_histogram = Counter()
    capacity_vectors = set()
    literal_survivor_states = 0
    digest = sha256()

    for support_tuple in combinations(range(1, P), 6):
        support = np.asarray(support_tuple, dtype=np.int8)
        capacities = np.zeros((len(words), 6), dtype=np.int16)
        for provider, label in enumerate(support):
            capacities += cards[
                label, words[:, provider, None], support[None, :]
            ]
        feasible = (capacities >= C).sum(axis=1)
        histogram += np.bincount(feasible, minlength=7)
        indices = np.flatnonzero(feasible == 6)
        support_histogram[len(indices)] += 1

        for row_index in indices:
            word = tuple(ORDERS[int(index)] for index in words[row_index])
            capacity = tuple(int(value) for value in capacities[row_index])
            survivors.append((support_tuple, word, capacity))
            profile_histogram[
                tuple(word.count(order) for order in ORDERS)
            ] += 1
            capacity_vectors.add(capacity)
            literal_survivor_states += prod(len(UNITS[order]) for order in word)
            digest.update(bytes(support_tuple))
            digest.update(bytes(word))
            for value in capacity:
                digest.update(value.to_bytes(2, "little"))

    return {
        "histogram": tuple(int(value) for value in histogram),
        "support_histogram": support_histogram,
        "survivors": tuple(survivors),
        "profile_histogram": profile_histogram,
        "capacity_vectors": capacity_vectors,
        "literal_survivor_states": literal_survivor_states,
        "digest": digest.hexdigest(),
    }


def relaxation_bound(support, word, owner, choices, anchors):
    bank = frozenset((0,))
    for label, order in zip(support, word):
        if order not in anchors:
            continue
        bank = frozenset(
            partial | option
            for partial in bank
            for option in choices[label, order, owner]
        )

    best = 0
    for anchor_union in bank:
        outside = FULL ^ anchor_union
        bound = anchor_union.bit_count()
        for label, order in zip(support, word):
            if order in anchors:
                continue
            bound += max(
                (mask & outside).bit_count()
                for mask in choices[label, order, owner]
            )
        best = max(best, bound)
    return best, len(bank)


def relaxation_census(survivors, choices):
    bounds = {name: Counter() for name in FLAGS}
    banks = {name: Counter() for name in FLAGS}
    live = {name: Counter() for name in FLAGS}
    own_order_bound = {name: Counter() for name in FLAGS}
    joint_live = Counter()
    tournament = Counter()
    digests = {name: sha256() for name in FLAGS}
    joint_digest = sha256()

    for support, word, capacities in survivors:
        row_bounds = {name: [] for name in FLAGS}
        row_banks = {name: [] for name in FLAGS}
        for owner_index, owner in enumerate(support):
            for name, anchors in FLAGS.items():
                bound, bank_size = relaxation_bound(
                    support, word, owner, choices, anchors
                )
                row_bounds[name].append(bound)
                row_banks[name].append(bank_size)
                bounds[name][bound] += 1
                banks[name][bank_size] += 1
                own_order_bound[name][word[owner_index], bound] += 1
                digests[name].update(bytes(support))
                digests[name].update(bytes(word))
                digests[name].update(bytes((owner, bound)))
                digests[name].update(bank_size.to_bytes(2, "little"))

        live_counts = {
            name: sum(value >= C for value in row_bounds[name])
            for name in FLAGS
        }
        for name in FLAGS:
            live[name][live_counts[name]] += 1
        joint_live[live_counts["Z4"], live_counts["Z9"]] += 1
        joint_digest.update(bytes(support))
        joint_digest.update(bytes(word))
        joint_digest.update(bytes((live_counts["Z4"], live_counts["Z9"])))

        gauge_edges = {}
        gauge_ties = {}
        for name in FLAGS:
            edges = set()
            ties = 0
            for left, right in combinations(range(6), 2):
                left_key = (
                    row_bounds[name][left] >= C,
                    row_bounds[name][left],
                    capacities[left],
                    row_banks[name][left],
                )
                right_key = (
                    row_bounds[name][right] >= C,
                    row_bounds[name][right],
                    capacities[right],
                    row_banks[name][right],
                )
                if left_key == right_key:
                    ties += 1
                    winner, loser = left, right
                elif left_key > right_key:
                    winner, loser = left, right
                else:
                    winner, loser = right, left
                edges.add((winner, loser))
            scores = sorted(
                sum(source == vertex for source, _target in edges)
                for vertex in range(6)
            )
            require(scores == list(range(6)), "owner tournament not transitive")
            gauge_edges[name] = edges
            gauge_ties[name] = ties

        flips = sum(
            ((left, right) in gauge_edges["Z4"])
            != ((left, right) in gauge_edges["Z9"])
            for left, right in combinations(range(6), 2)
        )
        tournament[gauge_ties["Z4"], gauge_ties["Z9"], flips] += 1

    require(joint_live[6, 6] == 0,
            "a scalar row survives both complementary flags")
    tournament_digest = sha256()
    for fingerprint, count in sorted(tournament.items()):
        tournament_digest.update(bytes(fingerprint))
        tournament_digest.update(count.to_bytes(8, "little"))
    return {
        "bounds": bounds,
        "banks": banks,
        "live": live,
        "own_order_bound": own_order_bound,
        "joint_live": joint_live,
        "tournament": tournament,
        "tournament_digest": tournament_digest.hexdigest(),
        "digests": {name: digest.hexdigest() for name, digest in digests.items()},
        "joint_digest": joint_digest.hexdigest(),
    }


def main() -> None:
    tables = build_tables()
    grammar = build_grammar()
    scalar = scalar_census(grammar["words"], tables["cards"])
    flags = relaxation_census(scalar["survivors"], tables["choices"])

    require(len(grammar["words"]) == 223_729, "hereditary word count")
    require(grammar["weighted_states"] == 1_904_124_672,
            "unit-weighted state count")
    require(scalar["histogram"] == (
        2_400_018, 49_776_708, 77_652_252, 56_649_096,
        17_848_230, 2_316_960, 82_332,
    ), "scalar feasible-owner histogram")
    require(len(scalar["survivors"]) == 82_332, "scalar survivor count")
    require(sum(count for rows, count in scalar["support_histogram"].items()
                if rows > 0) == 860, "scalar survivor support count")
    require(flags["live"]["Z4"] == Counter({
        0: 44_946, 1: 11_208, 2: 10_722, 3: 5_172,
        4: 2_280, 5: 180, 6: 7_824,
    }), "Z/4 live-owner histogram")
    require(flags["live"]["Z9"] == Counter({
        0: 45_600, 1: 10_752, 2: 4_428, 3: 5_436,
        4: 2_550, 5: 1_344, 6: 12_222,
    }), "Z/9 live-owner histogram")
    require(flags["joint_live"][6, 6] == 0,
            "complementary terminality")
    require(sum(flags["tournament"].values()) == 82_332,
            "owner tournament count")
    require(grammar["digest"] == EXPECTED_DIGESTS["grammar"],
            "grammar digest")
    require(tables["base_digest"] == EXPECTED_DIGESTS["base"],
            "CRT-base digest")
    require(tables["mask_digest"] == EXPECTED_DIGESTS["mask"],
            "literal-mask digest")
    require(scalar["digest"] == EXPECTED_DIGESTS["scalar"],
            "scalar digest")
    require(flags["digests"]["Z4"] == EXPECTED_DIGESTS["Z4"],
            "Z/4 relaxation digest")
    require(flags["digests"]["Z9"] == EXPECTED_DIGESTS["Z9"],
            "Z/9 relaxation digest")
    require(flags["joint_digest"] == EXPECTED_DIGESTS["joint"],
            "joint-live digest")

    print("### c=36 AP-centred H6 complementary Z/4-Z/9 primary ###")
    print("orders:", ORDERS)
    print("unit counts:", tuple(len(UNITS[order]) for order in ORDERS))
    print("grammar: at least two full 4-carriers and at least two full 9-carriers")
    print("hereditary words:", len(grammar["words"]))
    print("labelled support/order contexts:", 924 * len(grammar["words"]))
    print("literal state words/support:", grammar["weighted_states"])
    print("raw labelled state words:", 924 * grammar["weighted_states"])
    print("carrier histogram:", fmt(grammar["carrier_histogram"]))
    print("grammar digest:", grammar["digest"])
    print("CRT-base digest:", tables["base_digest"])
    print("literal-mask digest:", tables["mask_digest"])
    print("distinct option sizes (order,size):", fmt(tables["option_sizes"]))
    print("cardinality vectors by ratio; columns", ORDERS)
    for ratio in range(1, P):
        print(f"  r={ratio:2d}: {tables['rows'][ratio]}")
    print("scalar feasible-owner histogram:",
          " ".join(f"{i}:{value}" for i, value in enumerate(scalar["histogram"])))
    print("scalar support-survivor histogram:", fmt(scalar["support_histogram"]))
    print("scalar survivors:", len(scalar["survivors"]))
    print("scalar survivor supports:",
          sum(count for rows, count in scalar["support_histogram"].items()
              if rows > 0))
    print("scalar multiplicity profiles:", len(scalar["profile_histogram"]))
    print("scalar capacity vectors:", len(scalar["capacity_vectors"]))
    print("literal scalar-survivor state words:",
          scalar["literal_survivor_states"])
    print("scalar digest:", scalar["digest"])
    for name in ("Z4", "Z9"):
        print(f"{name} anchors:", tuple(sorted(FLAGS[name])))
        print(f"{name} bank sizes:", fmt(flags["banks"][name]))
        print(f"{name} bounds:", fmt(flags["bounds"][name]))
        print(f"{name} live owners/context:", fmt(flags["live"][name]))
        print(f"{name} digest:", flags["digests"][name])
    print("joint (Z4-live,Z9-live) histogram:", fmt(flags["joint_live"]))
    print("joint digest:", flags["joint_digest"])
    print("owner tournament fingerprints:", len(flags["tournament"]))
    print("owner tournament Z4 tie histogram:", fmt(Counter({
        ties: sum(count for (left, _right, _flips), count
                  in flags["tournament"].items() if left == ties)
        for ties in {key[0] for key in flags["tournament"]}
    })))
    print("owner tournament Z9 tie histogram:", fmt(Counter({
        ties: sum(count for (_left, right, _flips), count
                  in flags["tournament"].items() if right == ties)
        for ties in {key[1] for key in flags["tournament"]}
    })))
    print("owner tournament gauge-flip histogram:", fmt(Counter({
        flips: sum(count for (_left, _right, value), count
                   in flags["tournament"].items() if value == flips)
        for flips in {key[2] for key in flags["tournament"]}
    })))
    print("owner tournament digest:", flags["tournament_digest"])
    print("owner tournaments: scores (0,1,2,3,4,5), triangles 0, SCCs 1^6, HP 1")
    print("tournament vertices: owner obligations; gauge: Z/4 versus Z/9 proof cost")
    print("preserved: absolute owner-local upper bound")
    print("destroyed: shared nonanchor units and mutual nonanchor overlaps")
    print("complementary terminal cell (6,6):", flags["joint_live"][6, 6])
    print("VERDICT: primitive proper AP-centred common-scale-36 H6 face EMPTY")


if __name__ == "__main__":
    main()
