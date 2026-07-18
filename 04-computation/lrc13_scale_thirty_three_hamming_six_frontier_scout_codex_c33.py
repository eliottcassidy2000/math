#!/usr/bin/env python3
"""Exact primary for the c=33 AP-centred Hamming-six three-fibre obstruction.

All labelled supports and hereditary words are visited.  At the proof-facing
gate, masks of order 1 or 3 are retained exactly as thick fibres of
Z/33 -> Z/3; every order-11 or order-33 provider is then allowed to choose
its unit independently outside that anchor union.  The resulting sound upper
relaxation is below 33 at every owner whose own order carries 11.  Hereditary
lcm forces at least two such owners.

Tournament Analysis uses owner obligations, not runners, as its main vertex
set.  The pairwise cost observable is compared in the Z/3 and Z/11 gauges,
with coordinate order as the tie Hamiltonian path.  A second three-vertex
tournament compares quotient flags by their absolute owner-local bounds.
Both are telemetry: the numerical Z/3 threshold is the proof object.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 33
ORDERS = (1, 3, 11, 33)
UNITS = {
    d: ((0,) if d == 1 else tuple(u for u in range(1, d) if gcd(u, d) == 1))
    for d in ORDERS
}
FULL = (1 << C) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label: int, order: int, unit: int) -> int:
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
    require(len(answers) == 1, "literal CRT uniqueness")
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
    interval_count = sum(
        value % P == target for value in range(-order + 1, order + 1)
    )
    return (C // order) * interval_count


def build_tables():
    cards = np.zeros((P, len(ORDERS), P), dtype=np.int16)
    choices = {}
    rows = {}
    base_digest = sha256()
    mask_digest = sha256()
    for label in range(1, P):
        for order_index, order in enumerate(ORDERS):
            for unit in UNITS[order]:
                algebraic = crt_base(label, order, unit)
                literal = crt_base_literal(label, order, unit)
                require(algebraic == literal, "algebraic/literal CRT mismatch")
                base_digest.update(bytes((label, order, unit)))
                base_digest.update(algebraic.to_bytes(2, "little"))
            for owner in range(1, P):
                card = analytic_cardinality(label, order, owner)
                masks = tuple(sorted({
                    local_mask(label, order, unit, owner)
                    for unit in UNITS[order]
                }))
                require(all(mask.bit_count() == card for mask in masks),
                        "mask/cardinality mismatch")
                for mask in masks:
                    mask_digest.update(bytes((label, order, owner)))
                    mask_digest.update(mask.to_bytes(8, "little"))
                    for sheet in range(C):
                        require(
                            ((mask >> sheet) & 1)
                            == ((mask >> ((sheet + order) % C)) & 1),
                            "order periodicity mismatch",
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
                tuple(int(cards[label, i, owner]) for i in range(len(ORDERS)))
                == rows[ratio],
                "ratio covariance of scalar table",
            )
    return {
        "cards": cards,
        "choices": choices,
        "rows": rows,
        "base_digest": base_digest.hexdigest(),
        "mask_digest": mask_digest.hexdigest(),
    }


def hereditary(index_word) -> bool:
    word = tuple(ORDERS[i] for i in index_word)
    via_lcm = all(
        lcm(*(word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )
    via_carriers = (
        sum(order % 3 == 0 for order in word) >= 2
        and sum(order % 11 == 0 for order in word) >= 2
    )
    require(via_lcm == via_carriers, "hereditary grammar mismatch")
    return via_lcm


def build_grammar():
    words = []
    weighted = 0
    digest = sha256()
    for index_word in product(range(len(ORDERS)), repeat=6):
        if hereditary(index_word):
            word = tuple(ORDERS[i] for i in index_word)
            words.append(index_word)
            weight = prod(len(UNITS[d]) for d in word)
            weighted += weight
            digest.update(bytes(index_word))
            digest.update(weight.to_bytes(8, "little"))
    require(len(words) == 57**2 == 3249, "grammar count")
    require(weighted == (3**6 - 1 - 6*2) * (11**6 - 1 - 6*10),
            "weighted grammar count")
    return np.asarray(words, dtype=np.int8), weighted, digest.hexdigest()


def scalar_census(words, cards):
    histogram = np.zeros(7, dtype=np.int64)
    support_histogram = Counter()
    survivors = []
    multiplicities = Counter()
    capacity_vectors = set()
    digest = sha256()
    for support_tuple in combinations(range(1, P), 6):
        support = np.asarray(support_tuple, dtype=np.int8)
        capacities = np.zeros((len(words), 6), dtype=np.int16)
        for provider, label in enumerate(support):
            capacities += cards[label, words[:, provider, None], support[None, :]]
        feasible = (capacities >= C).sum(axis=1)
        histogram += np.bincount(feasible, minlength=7)
        indices = np.flatnonzero(feasible == 6)
        support_histogram[len(indices)] += 1
        for row_index in indices:
            word = tuple(ORDERS[int(i)] for i in words[row_index])
            capacity = tuple(int(v) for v in capacities[row_index])
            survivors.append((support_tuple, word, capacity))
            multiplicities[tuple(word.count(d) for d in ORDERS)] += 1
            capacity_vectors.add(capacity)
            digest.update(bytes(support_tuple))
            digest.update(bytes(word))
            digest.update(bytes(capacity))
    literal_survivors = sum(
        count * prod(
            len(UNITS[order]) ** multiplicity
            for order, multiplicity in zip(ORDERS, profile)
        )
        for profile, count in multiplicities.items()
    )
    return {
        "histogram": histogram,
        "support_histogram": support_histogram,
        "survivors": tuple(survivors),
        "multiplicities": multiplicities,
        "capacity_vectors": capacity_vectors,
        "literal_survivors": literal_survivors,
        "digest": digest.hexdigest(),
    }


def anchor_union_bank(support, word, owner, choices, anchors):
    bank = frozenset((0,))
    for label, order in zip(support, word):
        if order not in anchors:
            continue
        bank = frozenset(
            partial | option
            for partial in bank
            for option in choices[label, order, owner]
        )
    return bank


def relaxation_bound(support, word, owner, choices, anchors):
    bank = anchor_union_bank(support, word, owner, choices, anchors)
    best = 0
    for union in bank:
        outside = FULL ^ union
        bound = union.bit_count()
        for label, order in zip(support, word):
            if order in anchors:
                continue
            bound += max(
                (mask & outside).bit_count()
                for mask in choices[label, order, owner]
            )
        best = max(best, bound)
    return best, len(bank)


def flag_census(survivors, choices):
    flags = {
        "Z3": frozenset((1, 3)),
        "Z11": frozenset((1, 11)),
        "mixed": frozenset((1, 3, 11)),
    }
    bound_histograms = {name: Counter() for name in flags}
    bank_histograms = {name: Counter() for name in flags}
    live_histograms = {name: Counter() for name in flags}
    combined_live = Counter()
    z3_order_bound = Counter()
    z3_profile_live = Counter()
    tournament_histogram = Counter()
    flag_rank_histogram = Counter()
    flag_pair_wins = Counter()
    flag_digests = {name: sha256() for name in flags}
    residual = []
    for support, word, capacities in survivors:
        row_bounds = {name: [] for name in flags}
        row_banks = {name: [] for name in flags}
        for owner in support:
            for name, anchors in flags.items():
                bound, bank_size = relaxation_bound(
                    support, word, owner, choices, anchors
                )
                row_bounds[name].append(bound)
                row_banks[name].append(bank_size)
                bound_histograms[name][bound] += 1
                bank_histograms[name][bank_size] += 1
                flag_digests[name].update(bytes(support))
                flag_digests[name].update(bytes(word))
                flag_digests[name].update(bytes((owner, bound)))
                flag_digests[name].update(bank_size.to_bytes(4, "little"))
        for name in flags:
            live_histograms[name][sum(v >= C for v in row_bounds[name])] += 1
        profile = tuple(word.count(d) for d in ORDERS)
        z3_profile_live[(profile, sum(v >= C for v in row_bounds["Z3"]))] += 1
        for owner_index, bound in enumerate(row_bounds["Z3"]):
            z3_order_bound[(word[owner_index], bound)] += 1
            flag_names = tuple(flags)
            # A smaller sound upper bound is the stronger quotient flag.
            ordered = tuple(sorted(
                range(3), key=lambda i: (row_bounds[flag_names[i]][owner_index], i)
            ))
            flag_rank_histogram[tuple(flag_names[i] for i in ordered)] += 1
            for left, right in combinations(range(3), 2):
                left_key = (row_bounds[flag_names[left]][owner_index], left)
                right_key = (row_bounds[flag_names[right]][owner_index], right)
                winner = left if left_key < right_key else right
                loser = right if winner == left else left
                flag_pair_wins[(flag_names[winner], flag_names[loser])] += 1

        def owner_edges(name):
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
            require(
                sorted(sum((a, b).count(vertex) == 1 and a == vertex
                           for a, b in edges) for vertex in range(6))
                == [0, 1, 2, 3, 4, 5],
                "owner tournament score mismatch",
            )
            return edges, ties

        edges3, ties3 = owner_edges("Z3")
        edges11, ties11 = owner_edges("Z11")
        flips = sum(
            ((left, right) in edges3) != ((left, right) in edges11)
            for left, right in combinations(range(6), 2)
        )
        tournament_histogram[(ties3, ties11, flips)] += 1
        live = sum(
            min(row_bounds[name][owner] for name in flags) >= C
            for owner in range(6)
        )
        combined_live[live] += 1
        if live == 6:
            residual.append((support, word, capacities, row_bounds))
    return {
        "bound_histograms": bound_histograms,
        "bank_histograms": bank_histograms,
        "live_histograms": live_histograms,
        "combined_live": combined_live,
        "residual": residual,
        "z3_order_bound": z3_order_bound,
        "z3_profile_live": z3_profile_live,
        "tournament_histogram": tournament_histogram,
        "flag_rank_histogram": flag_rank_histogram,
        "flag_pair_wins": flag_pair_wins,
        "flag_digests": {
            name: digest.hexdigest() for name, digest in flag_digests.items()
        },
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
        require(orbit <= contexts, "context orbit escapes scalar bank")
        context_histogram[len(orbit)] += 1
        remaining -= orbit
    supports = {support for support, _word, _capacity in survivors}
    remaining_supports = set(supports)
    support_histogram = Counter()
    while remaining_supports:
        support = min(remaining_supports)
        orbit = {
            tuple(sorted(multiplier * label % P for label in support))
            for multiplier in range(1, P)
        }
        require(orbit <= supports, "support orbit escapes scalar bank")
        support_histogram[len(orbit)] += 1
        remaining_supports -= orbit
    return context_histogram, support_histogram


def fmt(counter) -> str:
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main() -> None:
    tables = build_tables()
    words, weighted, grammar_digest = build_grammar()
    scalar = scalar_census(words, tables["cards"])
    flags = flag_census(scalar["survivors"], tables["choices"])
    context_orbits, support_orbits = orbit_census(scalar["survivors"])

    # Complete theorem-facing census assertions.
    require(tuple(int(v) for v in scalar["histogram"])
            == (16_086, 446_568, 1_446_252, 891_684, 171_594, 28_548, 1_344),
            "scalar feasible-owner histogram")
    require(scalar["support_histogram"]
            == Counter({0: 774, 2: 48, 3: 24, 4: 24, 9: 24, 18: 12, 36: 18}),
            "scalar support histogram")
    require(len(scalar["survivors"]) == 1_344, "scalar survivor count")
    require(len({row[0] for row in scalar["survivors"]}) == 150,
            "scalar survivor support count")
    require(len(scalar["multiplicities"]) == 11, "multiplicity profile count")
    require(len(scalar["capacity_vectors"]) == 600, "capacity vector count")
    require(scalar["literal_survivors"] == 40_022_400,
            "literal scalar-survivor count")
    require(grammar_digest
            == "2900299ec2c09d25e41e3281575ca5de5fa5077388a34b8694bc0eda40d067c3",
            "grammar digest")
    require(tables["base_digest"]
            == "dc9e41ee276778a8b43d911c2e823109050c771285998307131d8fd39c148859",
            "CRT-base digest")
    require(tables["mask_digest"]
            == "bc85dd8ecfc3d5a178b7da29be292b203472d2a5d1dfbf1b3cbefda8f42e0554",
            "literal-mask digest")
    require(scalar["digest"]
            == "a1d26583a73929cd8f3bd66563741e09b550100dfe414d8334d71989cef280a4",
            "scalar-bank digest")
    require(flags["bank_histograms"]["Z3"]
            == Counter({1: 444, 2: 3_336, 3: 4_284}), "Z3 bank histogram")
    require(flags["bound_histograms"]["Z3"] == Counter({
        25: 996, 26: 4_248, 27: 384, 28: 288,
        29: 228, 30: 96, 33: 1_824,
    }), "Z3 bound histogram")
    require(flags["live_histograms"]["Z3"]
            == Counter({0: 306, 1: 360, 2: 624, 4: 54}),
            "Z3 live-owner histogram")
    require(flags["bank_histograms"]["Z11"] == Counter({
        1: 348, 10: 4_512, 55: 720, 80: 1_176, 90: 432,
        100: 624, 175: 72, 265: 24, 315: 48, 345: 24,
        355: 48, 385: 12, 675: 12, 715: 12,
    }), "Z11 bank histogram")
    require(flags["bound_histograms"]["Z11"] == Counter({
        26: 984, 27: 144, 28: 564, 29: 2_484, 30: 768,
        31: 156, 32: 48, 33: 768, 34: 156, 35: 348,
        36: 60, 37: 24, 38: 1_296, 39: 144, 43: 108, 48: 12,
    }), "Z11 bound histogram")
    require(flags["live_histograms"]["Z11"]
            == Counter({0: 30, 1: 192, 2: 912, 3: 72, 4: 72, 6: 66}),
            "Z11 live-owner histogram")
    require(flags["bound_histograms"]["mixed"] == Counter({
        25: 1_188, 26: 4_152, 27: 288, 28: 288,
        29: 228, 30: 96, 33: 1_824,
    }), "mixed bound histogram")
    require(flags["live_histograms"]["mixed"]
            == Counter({0: 306, 1: 360, 2: 624, 4: 54}),
            "mixed live-owner histogram")
    require(flags["flag_digests"] == {
        "Z3": "af6ad015c829d6b8dd6ffad44752297f51ac6f2535548f1b3dea3b7b3b900564",
        "Z11": "eb3204cee4eaec9f13eaf9f7e44d625423a50b1c066b2af05ff424fd8bd7da5f",
        "mixed": "6ec6a536a05eee349c05e0a0cb351a36754397ec352cf309d5e9e6c3a0162725",
    }, "flag digests")
    require(max(
        bound for (order, bound), count in flags["z3_order_bound"].items()
        if count and order % 11 == 0
    ) == 29, "11-carrier Z3 ceiling")
    require(all(
        sum(order % 11 == 0 for order in word) >= 2
        for _support, word, _capacity in scalar["survivors"]
    ), "hereditary 11-carrier count")
    require(not flags["residual"], "all-owner flag residual")
    require(context_orbits == Counter({6: 4, 12: 110}),
            "context orbit histogram")
    require(support_orbits == Counter({6: 1, 12: 12}),
            "support orbit histogram")
    require(sum(flags["tournament_histogram"].values()) == 1_344,
            "owner tournament count")
    require(max(key[2] for key in flags["tournament_histogram"]) == 9,
            "owner gauge-flip maximum")
    require(flags["flag_rank_histogram"] == Counter({
        ("Z3", "mixed", "Z11"): 6_228,
        ("Z3", "Z11", "mixed"): 1_644,
        ("mixed", "Z3", "Z11"): 192,
    }), "flag ranking histogram")

    print("### c=33 AP-centred H6 exact three-fibre primary ###")
    print("orders:", ORDERS)
    print("unit counts:", tuple(len(UNITS[d]) for d in ORDERS))
    print("grammar: at least two 3-carriers and at least two 11-carriers")
    print("hereditary words:", len(words))
    print("labelled support/order contexts:", 924 * len(words))
    print("literal state words/support:", weighted)
    print("raw labelled state words:", 924 * weighted)
    print("grammar digest:", grammar_digest)
    print("CRT-base digest:", tables["base_digest"])
    print("literal-mask digest:", tables["mask_digest"])
    print("cardinality vectors by ratio; columns D=1,3,11,33:")
    for ratio in range(1, P):
        print(f"  r={ratio:2d}: {tables['rows'][ratio]}")
    print("scalar feasible-owner histogram:",
          " ".join(f"{i}:{int(v)}" for i, v in enumerate(scalar["histogram"])))
    print("scalar support-survivor histogram:", fmt(scalar["support_histogram"]))
    print("scalar survivors:", len(scalar["survivors"]))
    print("scalar survivor supports:", len({row[0] for row in scalar["survivors"]}))
    print("scalar multiplicity profiles:", fmt(scalar["multiplicities"]))
    print("scalar capacity vectors:", len(scalar["capacity_vectors"]))
    print("literal scalar-survivor state words:", scalar["literal_survivors"])
    print("scalar digest:", scalar["digest"])
    for name in ("Z3", "Z11", "mixed"):
        print(f"{name} bank sizes:", fmt(flags["bank_histograms"][name]))
        print(f"{name} bounds:", fmt(flags["bound_histograms"][name]))
        print(f"{name} live owners/context:",
              fmt(flags["live_histograms"][name]))
        print(f"{name} digest:", flags["flag_digests"][name])
    print("best-of-flags live owners/context:", fmt(flags["combined_live"]))
    print("Z3 owner-order/bound:", fmt(flags["z3_order_bound"]))
    print("Z3 profile/live:", fmt(flags["z3_profile_live"]))
    print("all-owner residual:", len(flags["residual"]))
    print("residual multiplicities:", fmt(Counter(
        tuple(word.count(d) for d in ORDERS)
        for _support, word, _capacities, _bounds in flags["residual"]
    )))
    print("context multiplication orbit sizes:", fmt(context_orbits))
    print("support multiplication orbit sizes:", fmt(support_orbits))
    print("owner tournament (Z3 ties,Z11 ties,gauge flips):",
          fmt(flags["tournament_histogram"]))
    print("owner tournaments: score (0,1,2,3,4,5), triangles 0, SCCs 1^6, HP 1")
    print("flag strength rankings:", fmt(flags["flag_rank_histogram"]))
    print("flag pair wins:", fmt(flags["flag_pair_wins"]))
    print("flag tournaments: score (0,1,2), triangles 0, SCCs 1^3, HP 1")
    print("tournament vertices: owner obligations; alternate vertices: quotient flags")
    print("observable/switch: proof cost in Z3 versus Z11 gauge; coordinate tie path")
    print("preserved: owner-local absolute upper bound; destroyed: shared units/needle overlaps")
    print("11-carrier owner Z3 ceiling: 29 < 33")
    print("THEOREM: hereditary grammar supplies at least two terminal 11-carrier owners")
    print("VERDICT: primitive proper AP-centred common-scale-33 H6 face EMPTY")


if __name__ == "__main__":
    main()
