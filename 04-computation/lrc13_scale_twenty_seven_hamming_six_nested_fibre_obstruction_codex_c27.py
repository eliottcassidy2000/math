#!/usr/bin/env python3
"""Exact primary for the c=27 AP-centred Hamming-six nested-fibre obstruction.

The proof-facing obstruction is a nested 3-adic fibre relaxation.  Rows with
an order-three provider are projected through Z/27 -> Z/3.  The thirty scalar
rows without order three are projected through Z/27 -> Z/9.  At a fixed owner,
the anchor-fibre union is kept literally while every remaining provider is
allowed to maximize its contribution outside that union independently.  This
only enlarges the possible union, so an upper bound below 27 is terminal.

NumPy is used only for exact batching of the 924*1909 scalar contexts.  CRT
bases, literal masks, fibre bounds, immutable-set reachability, orbits, and
tournament telemetry use Python integers and deterministic iteration.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 27
ORDERS = (1, 3, 9, 27)
UNITS = {
    order: (
        (0,)
        if order == 1
        else tuple(unit for unit in range(1, order) if gcd(unit, order) == 1)
    )
    for order in ORDERS
}
FULL = (1 << C) - 1
ORDER_THREE_HIGH = frozenset((1, 4, 5, 8, 9))
ORDER_NINE_HIGH = frozenset((1, 2, 5, 8, 11))


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
        return 27 if ratio == 1 else 0
    if order == 3:
        return 9 if ratio in ORDER_THREE_HIGH else 0
    if order == 9:
        return 6 if ratio in ORDER_NINE_HIGH else 3
    require(order == 27, "unknown order")
    return 5 if ratio == 1 else 4


def build_tables():
    masks = {}
    cards = np.zeros((P, len(ORDERS), P), dtype=np.int16)
    mask_digest = sha256()
    base_digest = sha256()
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
                            "closed cardinality mismatch")
                    require(mask.bit_count() == card, "mask/cardinality mismatch")
                    masks[label, order, unit, owner] = mask
                    cards[label, order_index, owner] = card
                    mask_digest.update(bytes((label, order, unit, owner)))
                    mask_digest.update(mask.to_bytes(4, "little"))

                    if order in (3, 9):
                        for sheet in range(C):
                            require(
                                ((mask >> sheet) & 1)
                                == ((mask >> ((sheet + order) % C)) & 1),
                                "fibre periodicity",
                            )
                        occupied = tuple(
                            residue
                            for residue in range(order)
                            if (mask >> residue) & 1
                        )
                        require(
                            len(occupied) * (C // order) == card,
                            "fibre-cardinality identity",
                        )
    return (
        masks,
        cards,
        mask_digest.hexdigest(),
        base_digest.hexdigest(),
    )


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
        via_lcm = hereditary(index_word)
        via_prime_cube = index_word.count(3) >= 2
        require(via_lcm == via_prime_cube, "prime-cube grammar mismatch")
        if not via_lcm:
            continue
        words.append(index_word)
        weight = prod(len(UNITS[ORDERS[index]]) for index in index_word)
        state_words += weight
        digest.update(bytes(index_word))
        digest.update(weight.to_bytes(8, "little"))
    require(len(words) == 1909, "hereditary word count")
    require(state_words == 380_511_756, "literal state words/support")
    return np.asarray(words, dtype=np.int8), state_words, digest.hexdigest()


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
    expected_histogram = (11452, 391056, 849558, 418988, 83472, 8940, 450)
    require(tuple(int(x) for x in feasible_owner_histogram) == expected_histogram,
            "scalar feasible-owner histogram")
    require(support_histogram == Counter({0: 840, 3: 6, 4: 24, 5: 24,
                                          6: 24, 12: 6}),
            "support survivor histogram")
    require(len(survivors) == 450, "scalar survivor count")
    require(len({support for support, _, _ in survivors}) == 84,
            "scalar support count")
    require(all(1 not in word for _, word, _ in survivors),
            "order-one scalar survivor")
    require(
        multiplicities
        == Counter({
            (0, 2, 2, 2): 294,
            (0, 2, 1, 3): 60,
            (0, 3, 1, 2): 36,
            (0, 4, 0, 2): 18,
            (0, 0, 4, 2): 18,
            (0, 0, 3, 3): 12,
            (0, 3, 0, 3): 12,
        }),
        "scalar multiplicity histogram",
    )
    return {
        "survivors": tuple(survivors),
        "feasible_owner_histogram": tuple(int(x) for x in feasible_owner_histogram),
        "support_histogram": support_histogram,
        "multiplicities": multiplicities,
        "capacity_vectors": capacity_vectors,
        "digest": digest.hexdigest(),
    }


def multiply_context(context, multiplier):
    return tuple(sorted((multiplier * label % P, order) for label, order in context))


def multiplication_orbits(survivors):
    contexts = {
        tuple(zip(support, word))
        for support, word, _capacity in survivors
    }
    remaining = set(contexts)
    context_histogram = Counter()
    multiplicity_orbits = Counter()
    representatives = []
    while remaining:
        context = min(remaining)
        orbit = {multiply_context(context, multiplier) for multiplier in range(1, P)}
        require(orbit <= contexts, "context multiplication orbit escapes")
        remaining -= orbit
        context_histogram[len(orbit)] += 1
        multiplicity = tuple(
            sum(order == candidate for _label, order in context)
            for candidate in ORDERS
        )
        multiplicity_orbits[multiplicity, len(orbit)] += 1
        representatives.append((context, len(orbit)))

    supports = {support for support, _word, _capacity in survivors}
    support_remaining = set(supports)
    support_histogram = Counter()
    support_representatives = []
    while support_remaining:
        support = min(support_remaining)
        orbit = {
            tuple(sorted(multiplier * label % P for label in support))
            for multiplier in range(1, P)
        }
        require(orbit <= supports, "support multiplication orbit escapes")
        support_remaining -= orbit
        support_histogram[len(orbit)] += 1
        support_representatives.append((support, len(orbit)))

    require(context_histogram == Counter({6: 3, 12: 36}), "context orbits")
    require(support_histogram == Counter({6: 2, 12: 6}), "support orbits")
    return {
        "context_histogram": context_histogram,
        "multiplicity_orbits": multiplicity_orbits,
        "representatives": tuple(representatives),
        "support_histogram": support_histogram,
        "support_representatives": tuple(support_representatives),
    }


def fibre_relaxation_bound(support, word, owner, anchor_order, masks):
    """Sound union upper bound retaining only anchor-fibre incidence."""
    anchor_indices = [index for index, order in enumerate(word) if order == anchor_order]
    require(anchor_indices, "empty anchor family")
    anchor_unions = frozenset((0,))
    for index in anchor_indices:
        label = support[index]
        options = frozenset(
            masks[label, anchor_order, unit, owner]
            for unit in UNITS[anchor_order]
        )
        anchor_unions = frozenset(
            partial | option
            for partial in anchor_unions
            for option in options
        )

    best = 0
    for anchor_union in anchor_unions:
        bound = anchor_union.bit_count()
        for index, (label, order) in enumerate(zip(support, word)):
            if index in anchor_indices:
                continue
            bound += max(
                (masks[label, order, unit, owner] & ~anchor_union).bit_count()
                for unit in UNITS[order]
            )
        best = max(best, bound)
    return best, len(anchor_unions)


def structural_fibre_census(survivors, masks):
    bound_histogram = Counter()
    context_feasible_histogram = Counter()
    carrier_contexts = Counter()
    carrier_owner_bounds = {3: Counter(), 9: Counter()}
    carrier_union_banks = {3: Counter(), 9: Counter()}
    rows = {}
    digest = sha256()
    for support, word, _capacity in survivors:
        carrier = 3 if 3 in word else 9
        carrier_contexts[carrier] += 1
        bounds = []
        for owner in support:
            bound, union_count = fibre_relaxation_bound(
                support, word, owner, carrier, masks
            )
            bounds.append(bound)
            bound_histogram[bound] += 1
            carrier_owner_bounds[carrier][bound] += 1
            carrier_union_banks[carrier][union_count] += 1
            digest.update(bytes(support))
            digest.update(bytes(word))
            digest.update(bytes((owner, carrier, bound)))
            digest.update(union_count.to_bytes(2, "little"))
        feasible_count = sum(bound >= C for bound in bounds)
        context_feasible_histogram[feasible_count] += 1
        rows[support, word] = tuple(bounds)

    require(carrier_contexts == Counter({3: 420, 9: 30}), "carrier split")
    require(
        context_feasible_histogram == Counter({0: 48, 1: 96, 3: 96, 4: 210}),
        "structural feasible-owner histogram",
    )
    require(max(context_feasible_histogram) == 4, "structural all-owner survivor")
    return {
        "bounds": rows,
        "bound_histogram": bound_histogram,
        "context_feasible_histogram": context_feasible_histogram,
        "carrier_contexts": carrier_contexts,
        "carrier_owner_bounds": carrier_owner_bounds,
        "carrier_union_banks": carrier_union_banks,
        "digest": digest.hexdigest(),
    }


def owner_local(support, word, owner, masks):
    reachable = frozenset((0,))
    layers = []
    for label, order in zip(support, word):
        options = frozenset(
            masks[label, order, unit, owner]
            for unit in UNITS[order]
        )
        reachable = frozenset(
            partial | option
            for partial in reachable
            for option in options
        )
        layers.append(len(reachable))
    maximum = max(mask.bit_count() for mask in reachable)
    feasible = FULL in reachable
    require(feasible == (maximum == C), "feasibility/maximum mismatch")
    maximizing = sum(mask.bit_count() == maximum for mask in reachable)
    return feasible, maximum, len(reachable), maximizing, tuple(layers), reachable


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
            vertex
            for vertex in unused
            if reach[root][vertex] and reach[vertex][root]
        }
        require(component, "empty SCC")
        result.append(len(component))
        unused -= component
    return tuple(sorted(result))


def tournament_fingerprint(local_rows, capacities, structural_bounds):
    adjacency = [[False] * 6 for _ in range(6)]
    ties = 0
    flips = 0
    for left, right in combinations(range(6), 2):
        left_key = (
            int(local_rows[left][0]),
            local_rows[left][1],
            capacities[left],
            structural_bounds[left],
        )
        right_key = (
            int(local_rows[right][0]),
            local_rows[right][1],
            capacities[right],
            structural_bounds[right],
        )
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
                if (
                    (previous_mask >> previous) & 1
                    and adjacency[previous][last]
                ):
                    paths[mask][last] += paths[previous_mask][previous]
    hamiltonian_paths = sum(paths[-1])
    components = scc_sizes(adjacency)
    require(scores == (0, 1, 2, 3, 4, 5), "tournament scores")
    require(triangles == 0, "tournament cycle")
    require(components == (1, 1, 1, 1, 1, 1), "tournament SCCs")
    require(hamiltonian_paths == 1, "tournament Hamiltonian paths")
    return ties, flips


def exact_owner_census(survivors, structural, masks):
    feasible_context_histogram = Counter()
    maximum_histogram = Counter()
    bank_size_histogram = Counter()
    maximizing_histogram = Counter()
    layer_profiles = Counter()
    tie_histogram = Counter()
    flip_histogram = Counter()
    total_masks = 0
    maximum_bank = 0
    digest = sha256()
    for support, word, capacities in survivors:
        local_rows = []
        bounds = structural["bounds"][support, word]
        for owner, bound in zip(support, bounds):
            row = owner_local(support, word, owner, masks)
            feasible, maximum, bank_size, maximizing, layers, reachable = row
            require(maximum <= bound, "exact maximum exceeds structural bound")
            local_rows.append(row)
            maximum_histogram[maximum] += 1
            bank_size_histogram[bank_size] += 1
            maximizing_histogram[maximizing] += 1
            layer_profiles[layers] += 1
            total_masks += bank_size
            maximum_bank = max(maximum_bank, bank_size)
            digest.update(bytes(support))
            digest.update(bytes(word))
            digest.update(bytes((owner, int(feasible), maximum)))
            digest.update(bank_size.to_bytes(4, "little"))
            digest.update(maximizing.to_bytes(4, "little"))
            for mask in sorted(reachable):
                digest.update(mask.to_bytes(4, "little"))
        feasible_count = sum(row[0] for row in local_rows)
        feasible_context_histogram[feasible_count] += 1
        ties, flips = tournament_fingerprint(local_rows, capacities, bounds)
        tie_histogram[ties] += 1
        flip_histogram[flips] += 1

    require(
        feasible_context_histogram == Counter({0: 336, 1: 96, 4: 18}),
        "exact feasible-owner contexts",
    )
    require(
        maximum_histogram
        == Counter({20: 120, 21: 336, 22: 192, 23: 336,
                    24: 528, 25: 432, 26: 588, 27: 168}),
        "exact owner maxima",
    )
    require(total_masks == 13_598_160 and maximum_bank == 128_880,
            "reachable-mask totals")
    return {
        "feasible_context_histogram": feasible_context_histogram,
        "maximum_histogram": maximum_histogram,
        "bank_size_histogram": bank_size_histogram,
        "maximizing_histogram": maximizing_histogram,
        "layer_profiles": layer_profiles,
        "tie_histogram": tie_histogram,
        "flip_histogram": flip_histogram,
        "total_masks": total_masks,
        "maximum_bank": maximum_bank,
        "digest": digest.hexdigest(),
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
        cyclic = False
        for first, second, third in (
            indices,
            (indices[0], indices[2], indices[1]),
        ):
            cyclic |= (
                adjacency[first][second]
                and adjacency[second][third]
                and adjacency[third][first]
            )
        directed_triangles += int(cyclic)
    return (
        len(arcs),
        len(symmetric_edges),
        len(reciprocal_edges),
        directed_triangles,
        scc_sizes(adjacency),
    )


def quadratic_cayley_audit():
    """Expose the K6,6 shadow of the order-nine high-ratio switch."""
    group = frozenset(range(1, P))
    quadratic = frozenset(value * value % P for value in group)
    nonquadratic = group - quadratic
    directed = ORDER_NINE_HIGH - {1}
    symmetric = directed | frozenset(pow(value, -1, P) for value in directed)
    require(symmetric == nonquadratic, "order-nine symmetric switch")
    edges = {
        tuple(sorted((left, right)))
        for left, right in combinations(range(1, P), 2)
        if right * pow(left, -1, P) % P in symmetric
    }
    require(len(edges) == 36, "K6,6 edge count")
    require(
        all((left in quadratic) != (right in quadratic) for left, right in edges),
        "K6,6 bipartition",
    )
    return tuple(sorted(quadratic)), tuple(sorted(nonquadratic)), len(edges)


def fmt(counter):
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main():
    masks, cards, mask_digest, base_digest = build_tables()
    words, state_words, grammar_digest = grammar()
    scalar = scalar_census(words, cards)
    orbits = multiplication_orbits(scalar["survivors"])
    structural = structural_fibre_census(scalar["survivors"], masks)
    exact = exact_owner_census(scalar["survivors"], structural, masks)
    relation3 = relation_fingerprint(ORDER_THREE_HIGH - {1})
    relation9 = relation_fingerprint(ORDER_NINE_HIGH - {1})
    quadratic_cayley = quadratic_cayley_audit()

    require(
        grammar_digest == "b35f76c33a4cb79268765c56b058a4e8c0ed0ab4bb813bcdb9fe7032b9aca4c6",
        "grammar digest",
    )
    require(
        base_digest == "21576171f9d407741077323d8519909349d068af06eb9921e6554b0b73822d32",
        "CRT-base digest",
    )
    require(
        mask_digest == "2f16bfe8cb3191e2614c4c8040d94c719eca4feaa8534b299fb5e5008de1a7a4",
        "literal-mask digest",
    )
    require(
        scalar["digest"] == "141ec6d8c551c2e0ebac31dc102d3f38ad257be268dfe24403835afa7613dc05",
        "scalar digest",
    )
    require(
        structural["digest"] == "14f91a7e4dcb828f5508a2d0c6d45a966652dd7abe9e5ae6936f27f05a72ed3a",
        "structural-fibre digest",
    )
    require(
        exact["digest"] == "ea89048d693763714050844ba0748e0f7128897a6b8ac52831235f4419a91497",
        "exact-owner digest",
    )
    require(
        structural["bound_histogram"]
        == Counter({20: 24, 21: 84, 22: 492, 23: 48, 24: 168,
                    25: 528, 26: 132, 27: 1128, 28: 96}),
        "structural owner bounds",
    )
    require(len(exact["bank_size_histogram"]) == 51, "reachable-bank bins")
    require(len(exact["layer_profiles"]) == 1737, "layer profiles")
    require(
        exact["tie_histogram"] == Counter({0: 96, 1: 168, 2: 36,
                                            3: 102, 4: 12, 7: 36}),
        "tournament ties",
    )
    require(
        exact["flip_histogram"]
        == Counter({0: 10, 1: 9, 2: 21, 3: 34, 4: 55, 5: 49,
                    6: 64, 7: 44, 8: 30, 9: 48, 10: 34, 11: 16,
                    12: 20, 13: 11, 14: 4, 15: 1}),
        "tournament flips",
    )
    require(relation3 == (48, 36, 12, 16, (12,)), "order-three Cayley relation")
    require(relation9 == (48, 36, 12, 0, (12,)), "order-nine Cayley relation")

    literal_survivor_words = sum(
        count * prod(len(UNITS[order]) ** multiplicity
                     for order, multiplicity in zip(ORDERS, pattern))
        for pattern, count in scalar["multiplicities"].items()
    )
    require(literal_survivor_words == 46_002_816, "literal survivor words")

    print("scale-twenty-seven AP-centred Hamming-six prime-cube nested-fibre certificate")
    print("orders", ORDERS, "unit counts", tuple(len(UNITS[d]) for d in ORDERS))
    print("numpy exact-batching version", np.__version__)
    print("hereditary words", len(words), "state words/support", state_words,
          "raw labelled states", 924 * state_words)
    print("grammar SHA256", grammar_digest)
    print("CRT-base SHA256", base_digest)
    print("literal-mask SHA256", mask_digest)
    print("closed scalar cards D1 27 at r=1 else 0; D3 9 on", tuple(sorted(ORDER_THREE_HIGH)),
          "else 0; D9 6 on", tuple(sorted(ORDER_NINE_HIGH)),
          "else 3; D27 5 at r=1 else 4")
    print("scalar feasible-owner histogram", " ".join(
        f"{index}:{value}" for index, value in enumerate(scalar["feasible_owner_histogram"])
    ))
    print("scalar supports-by-survivor-count", fmt(scalar["support_histogram"]))
    print("scalar survivors", len(scalar["survivors"]), "on supports",
          len({support for support, _word, _capacity in scalar["survivors"]}))
    print("scalar multiplicity histogram", fmt(scalar["multiplicities"]))
    print("literal unit words across scalar survivors", literal_survivor_words)
    print("distinct scalar capacity vectors", len(scalar["capacity_vectors"]))
    print("scalar SHA256", scalar["digest"])
    print("context multiplication orbit-size histogram", fmt(orbits["context_histogram"]),
          "orbits", len(orbits["representatives"]))
    print("context orbit/multiplicity histogram", fmt(orbits["multiplicity_orbits"]))
    print("support multiplication orbit-size histogram", fmt(orbits["support_histogram"]),
          "orbits", len(orbits["support_representatives"]))
    print("size-six context orbit representatives")
    for context, size in orbits["representatives"]:
        if size == 6:
            print(" ", context)
    print("nested-fibre carrier contexts", fmt(structural["carrier_contexts"]))
    print("mod3 anchor-union bank-size histogram",
          fmt(structural["carrier_union_banks"][3]))
    print("mod9 anchor-union bank-size histogram",
          fmt(structural["carrier_union_banks"][9]))
    print("mod3 owner-bound histogram", fmt(structural["carrier_owner_bounds"][3]))
    print("mod9 owner-bound histogram", fmt(structural["carrier_owner_bounds"][9]))
    print("combined structural owner-bound histogram", fmt(structural["bound_histogram"]))
    print("structural feasible-owner/context histogram",
          fmt(structural["context_feasible_histogram"]))
    print("structural fibre SHA256", structural["digest"])
    print("exact feasible-owner/context histogram", fmt(exact["feasible_context_histogram"]))
    print("exact maximum-union histogram", fmt(exact["maximum_histogram"]))
    print("exact reachable-bank-size histogram", fmt(exact["bank_size_histogram"]))
    print("exact owner rows", 6 * len(scalar["survivors"]),
          "total reachable masks", exact["total_masks"],
          "maximum bank", exact["maximum_bank"],
          "layer profiles", len(exact["layer_profiles"]))
    print("exact owner SHA256", exact["digest"])
    print("order-three high-ratio Cayley fingerprint arcs,symmetric,reciprocal,triangles,SCCs",
          relation3)
    print("order-nine high-ratio Cayley fingerprint arcs,symmetric,reciprocal,triangles,SCCs",
          relation9)
    print("order-nine symmetric high-ratio relation is K6,6 on QR/NQR",
          quadratic_cayley)
    print("tournament pair observable lexicographic (feasible,max-union,capacity,fibre-bound); switch by larger key and coordinate tie Hamiltonian path")
    print("tournament fingerprints all 450 transitive: scores 0,1,2,3,4,5; cycles 0; singleton SCCs; Hamiltonian paths 1")
    print("tournament tie-edge histogram", fmt(exact["tie_histogram"]))
    print("tournament flip-edge histogram", fmt(exact["flip_histogram"]))
    print("3-adic toothpick/Kakeya carrier Z/27->Z/3 retains D3 thick fibres and independently relaxes D9/D27 needles; Z/27->Z/9 resolves the 30 rows with no D3 anchors")
    print("information-loss ledger the fibre quotient forgets within-fibre positions, pair overlaps among relaxed needles, and shared unit choices, but only enlarges unions; literal DP is the sidecar. Owners retain the terminal deficit; runners, isolated residues, raw sheets, Fano points, and completed tournaments lose the absolute threshold or nested incidence")
    print("frontier verdict scalar-empty no; nested-fibre all-six empty yes; exact owner-local all-six empty yes")


if __name__ == "__main__":
    main()
