#!/usr/bin/env python3
"""Exact c=25 AP-centred Hamming-six prime-square obstruction certificate.

The proof-bearing layer is structural.  It derives the hereditary order
grammar and scalar classification, identifies a Cayley ``K6 disjoint-union
K6`` obstruction on the two quadratic classes, and projects the 25 literal
sheets to their five residue fibres.  The resulting Kakeya/toothpick incidence
count bounds every order-five owner union by 22 and every order-twenty-five
owner union by 21.  A full immutable-set union DP independently checks that
both bounds are sharp.  Multiplication orbits and tournaments are telemetry;
no support, order word, owner, or literal unit is quotiented out of a census.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod


P = 13
C = 25
ORDERS = (1, 5, 25)
UNITS = {
    d: ((0,) if d == 1 else tuple(e for e in range(1, d) if gcd(e, d) == 1))
    for d in ORDERS
}
FULL = (1 << C) - 1
GENERATOR_POWERS = tuple(pow(2, exponent, P) for exponent in range(12))
EXPONENT = {value: exponent for exponent, value in enumerate(GENERATOR_POWERS)}
FORBIDDEN_RATIO_5 = frozenset((4, 9, 12))
FORBIDDEN_OWNER_MULTIPLIER = frozenset((3, 10, 12))


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
    inverse = pow(owner, -1, P)
    mask = 0
    for sheet in range(C):
        value = centered(base * (inverse + P * sheet), P * order)
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


def build_tables():
    masks = {}
    cards = {}
    digest = sha256()
    for label in range(1, P):
        for order in ORDERS:
            for unit in UNITS[order]:
                require(
                    crt_base_algebraic(label, order, unit)
                    == crt_base_literal(label, order, unit),
                    "algebraic/literal CRT mismatch",
                )
                for owner in range(1, P):
                    mask = local_mask(label, order, unit, owner)
                    card = analytic_cardinality(label, order, owner)
                    require(mask.bit_count() == card, "mask/card mismatch")
                    masks[label, order, unit, owner] = mask
                    cards[label, order, owner] = card
                    digest.update(bytes((label, order, unit, owner)))
                    digest.update(mask.to_bytes(4, "little"))
    return masks, cards, digest.hexdigest()


def cardinality_formula_audit(cards):
    """Check the three closed scalar formulas independently of CRT units."""
    for label in range(1, P):
        for owner in range(1, P):
            ratio = label * pow(owner, -1, P) % P
            expected = {
                1: 25 if ratio == 1 else 0,
                5: 0 if ratio in FORBIDDEN_RATIO_5 else 5,
                25: 3 if ratio == 12 else 4,
            }
            for order in ORDERS:
                require(
                    cards[label, order, owner] == expected[order],
                    "closed scalar-cardinality formula mismatch",
                )


def cayley_audit():
    """Audit the proof-facing forbidden-ratio Cayley carrier."""
    group = frozenset(range(1, P))
    quadratic = frozenset(value * value % P for value in group)
    nonquadratic = group - quadratic
    inverse_forbidden = frozenset(pow(value, -1, P) for value in FORBIDDEN_RATIO_5)
    symmetric_switch = FORBIDDEN_RATIO_5 | inverse_forbidden
    require(quadratic == frozenset((1, 3, 4, 9, 10, 12)), "quadratic class")
    require(nonquadratic == frozenset((2, 5, 6, 7, 8, 11)), "nonquadratic class")
    require(inverse_forbidden == FORBIDDEN_OWNER_MULTIPLIER, "inverse ratio switch")
    require(symmetric_switch == quadratic - {1}, "Cayley K6 switch")

    directed_arcs = set()
    symmetric_edges = set()
    reciprocal_edges = set()
    for owner in group:
        for provider in group:
            if owner == provider:
                continue
            ratio = provider * pow(owner, -1, P) % P
            if ratio in FORBIDDEN_RATIO_5:
                directed_arcs.add((owner, provider))
            if ratio in symmetric_switch:
                symmetric_edges.add(tuple(sorted((owner, provider))))
            if (
                ratio in FORBIDDEN_RATIO_5
                and pow(ratio, -1, P) in FORBIDDEN_RATIO_5
            ):
                reciprocal_edges.add(tuple(sorted((owner, provider))))
    require(len(directed_arcs) == 36, "directed Cayley arc count")
    require(len(symmetric_edges) == 30, "two K6 edge count")
    require(len(reciprocal_edges) == 6, "antipodal reciprocal matching")
    require(
        all(
            ((left in quadratic) == (right in quadratic))
            for left, right in symmetric_edges
        ),
        "Cayley edge crosses quadratic classes",
    )
    require(
        all(
            tuple(sorted((left, right))) in symmetric_edges
            for component in (quadratic, nonquadratic)
            for left, right in combinations(sorted(component), 2)
        ),
        "Cayley component is not complete",
    )
    return (
        len(directed_arcs),
        len(symmetric_edges),
        len(reciprocal_edges),
        tuple(sorted((len(quadratic), len(nonquadratic)))),
    )


def ratio_mask_signatures(masks):
    """Audit the proof-bearing reduction of literal masks modulo five."""
    signatures = {}
    for order in (5, 25):
        for ratio in range(1, P):
            rows = set()
            for unit in UNITS[order]:
                mask = masks[ratio, order, unit, 1]
                signature = tuple(
                    sum((mask >> sheet) & 1 for sheet in range(residue, C, 5))
                    for residue in range(5)
                )
                rows.add(signature)
            signatures[order, ratio] = frozenset(rows)

    self5 = signatures[5, 1]
    require(self5 == frozenset(((0, 0, 0, 5, 0),)), "self order-five coset")
    expected_moving5 = frozenset(
        tuple(5 * int(index == residue) for index in range(5))
        for residue in (0, 1, 2, 4)
    )
    for ratio in range(2, P):
        if ratio in FORBIDDEN_RATIO_5:
            require(
                signatures[5, ratio] == frozenset(((0, 0, 0, 0, 0),)),
                "forbidden order-five ratio",
            )
        else:
            require(signatures[5, ratio] == expected_moving5, "moving order-five coset")

    for ratio in range(1, P):
        for signature in signatures[25, ratio]:
            if ratio == 12:
                require(
                    signature[3] == 0
                    and sum(signature) == 3
                    and sum(value == 1 for value in signature) == 3,
                    "antipodal order-twenty-five signature",
                )
            elif ratio in (4, 9):
                require(signature == (1, 1, 1, 0, 1), "full noncentral signature")
            else:
                require(
                    signature[3] == 1
                    and sum(signature) == 4
                    and sum(value == 1 for value in signature) == 4,
                    "central order-twenty-five signature",
                )
    return signatures


def hereditary(word):
    return all(
        lcm(*(word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )


def grammar():
    words = []
    state_words = 0
    digest = sha256()
    for word in product(ORDERS, repeat=6):
        by_lcm = hereditary(word)
        by_prime_square = word.count(25) >= 2
        require(by_lcm == by_prime_square, "prime-square grammar mismatch")
        if by_lcm:
            words.append(word)
            weight = prod(len(UNITS[d]) for d in word)
            state_words += weight
            digest.update(bytes(word))
            digest.update(weight.to_bytes(8, "little"))
    require(len(words) == 473, "hereditary grammar size")
    require(state_words == 243_750_000, "literal state-word count")
    return tuple(words), state_words, digest.hexdigest()


def scalar_capacity(support, word, owner, cards):
    return sum(cards[label, order, owner] for label, order in zip(support, word))


def structural_scalar_predicate(support, word):
    """The derived c=25 scalar classification, independent of cards."""
    if word.count(1):
        return False
    if word.count(25) != 4:
        return False
    providers5 = tuple(label for label, order in zip(support, word) if order == 5)
    if len(providers5) != 2:
        return False
    if (EXPONENT[providers5[0]] - EXPONENT[providers5[1]]) % 2 == 0:
        return False
    forbidden = {
        multiplier * provider % P
        for provider in providers5
        for multiplier in FORBIDDEN_OWNER_MULTIPLIER
    }
    return set(support) == set(range(1, P)) - forbidden


def structural_owner_bound(support, word, owner):
    """Bound the local union using only five-coset occupancy signatures."""
    normalized = tuple(
        (label * pow(owner, -1, P) % P, order)
        for label, order in zip(support, word)
    )
    ratios5 = tuple(ratio for ratio, order in normalized if order == 5)
    ratios25 = tuple(ratio for ratio, order in normalized if order == 25)
    owner_order = word[support.index(owner)]
    require(len(ratios5) == 2 and len(ratios25) == 4, "structural owner grammar")
    if owner_order == 5:
        # The self order-five mask occupies Q_3.  The other occupies one of
        # Q_0,Q_1,Q_2,Q_4.  Exactly two order-25 providers have ratios 4,9,
        # hence hit the moving coset, while the two nonresidue providers each
        # hit Q_3.  Four of their sixteen points are already covered.
        require(1 in ratios5, "order-five self ratio")
        require(sum(ratio in (4, 9) for ratio in ratios25) == 2, "two full signatures")
        require(sum(ratio not in FORBIDDEN_RATIO_5 for ratio in ratios25) == 2,
                "two central signatures")
        return 10 + (16 - 4)

    require(owner_order == 25, "unexpected owner order")
    # The two order-five cosets are either equal (giving the stronger bound
    # 5+15=20) or distinct.  In the distinct case, each of the three central
    # four-point order-25 masks and the antipodal three-point mask meets their
    # two-coset union.  Again at least four points are already covered.
    require(all(ratio not in FORBIDDEN_RATIO_5 and ratio != 1 for ratio in ratios5),
            "order-five ratios at order-twenty-five owner")
    require(1 in ratios25 and 12 in ratios25, "self/antipodal pair")
    require(sum(ratio not in FORBIDDEN_RATIO_5 for ratio in ratios25) == 3,
            "three central signatures")
    return 10 + (15 - 4)


def scalar_census(words, cards):
    feasible_owner_histogram = Counter()
    survivors = []
    multiplicity_histogram = Counter()
    support_histogram = Counter()
    digest = sha256()
    for support in combinations(range(1, P), 6):
        count_on_support = 0
        for word in words:
            capacities = tuple(
                scalar_capacity(support, word, owner, cards) for owner in support
            )
            feasible_count = sum(value >= C for value in capacities)
            feasible_owner_histogram[feasible_count] += 1
            exact = feasible_count == 6
            structural = structural_scalar_predicate(support, word)
            require(exact == structural, "structural scalar predicate mismatch")
            if exact:
                survivors.append((support, word, capacities))
                count_on_support += 1
                multiplicity_histogram[
                    tuple(word.count(order) for order in ORDERS)
                ] += 1
                digest.update(bytes(support))
                digest.update(bytes(word))
                digest.update(bytes(capacities))
        support_histogram[count_on_support] += 1
    require(len(survivors) == 36, "scalar survivor count")
    require(len({support for support, _, _ in survivors}) == 36, "support count")
    require(
        multiplicity_histogram == Counter({(0, 2, 4): 36}),
        "survivor multiplicity",
    )
    require(
        feasible_owner_histogram
        == Counter({0: 1156, 1: 149868, 2: 171636, 3: 90884,
                    4: 21864, 5: 1608, 6: 36}),
        "scalar feasible-owner histogram",
    )
    require(
        support_histogram == Counter({0: 888, 1: 36}),
        "scalar support histogram",
    )
    return {
        "survivors": tuple(survivors),
        "feasible_owner_histogram": feasible_owner_histogram,
        "support_histogram": support_histogram,
        "multiplicity_histogram": multiplicity_histogram,
        "digest": digest.hexdigest(),
    }


def owner_local(support, word, owner, masks):
    reachable = frozenset((0,))
    layers = []
    for label, order in zip(support, word):
        options = frozenset(
            masks[label, order, unit, owner] for unit in UNITS[order]
        )
        reachable = frozenset(
            partial | option for partial in reachable for option in options
        )
        layers.append(len(reachable))
    maxima = max(mask.bit_count() for mask in reachable)
    maximizing = sum(mask.bit_count() == maxima for mask in reachable)
    return FULL in reachable, maxima, len(reachable), maximizing, tuple(layers), reachable


def tournament_fingerprint(rows, capacities):
    adjacency = [[False] * 6 for _ in range(6)]
    ties = 0
    flips = 0
    for left, right in combinations(range(6), 2):
        left_key = (rows[left][0], rows[left][1], capacities[left])
        right_key = (rows[right][0], rows[right][1], capacities[right])
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

    reach = [[adjacency[source][target] for target in range(6)] for source in range(6)]
    for vertex in range(6):
        reach[vertex][vertex] = True
    for middle in range(6):
        for source in range(6):
            if reach[source][middle]:
                for target in range(6):
                    reach[source][target] |= reach[middle][target]
    unused = set(range(6))
    scc_sizes = []
    while unused:
        root = min(unused)
        component = {
            vertex
            for vertex in unused
            if reach[root][vertex] and reach[vertex][root]
        }
        require(component, "empty tournament SCC")
        scc_sizes.append(len(component))
        unused -= component

    paths = [[0] * 6 for _ in range(1 << 6)]
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

    require(scores == (0, 1, 2, 3, 4, 5), "non-transitive tournament")
    require(triangles == 0, "tournament triangle")
    require(tuple(sorted(scc_sizes)) == (1, 1, 1, 1, 1, 1), "tournament SCCs")
    require(hamiltonian_paths == 1, "tournament Hamiltonian paths")
    return ties, flips, scores, triangles, tuple(sorted(scc_sizes)), hamiltonian_paths


def support_orbits(survivors):
    rows = {(support, tuple(label for label, d in zip(support, word) if d == 5))
            for support, word, _ in survivors}
    remaining = set(rows)
    histogram = Counter()
    invariants = Counter()
    while remaining:
        support, providers5 = min(remaining)
        orbit = set()
        for multiplier in range(1, P):
            new_support = tuple(sorted(multiplier * x % P for x in support))
            new_providers = tuple(sorted(multiplier * x % P for x in providers5))
            orbit.add((new_support, new_providers))
        require(orbit <= rows, "multiplication orbit escapes")
        remaining -= orbit
        histogram[len(orbit)] += 1
        ratio = providers5[1] * pow(providers5[0], -1, P) % P
        invariants[min(ratio, pow(ratio, -1, P))] += 1
    representatives = []
    for ratio in (2, 5, 6):
        providers5 = (1, ratio)
        forbidden = {
            multiplier * provider % P
            for provider in providers5
            for multiplier in FORBIDDEN_OWNER_MULTIPLIER
        }
        support = tuple(sorted(set(range(1, P)) - forbidden))
        representatives.append((providers5, support))
    require(
        all((support, providers5) in rows for providers5, support in representatives),
        "declared orbit representative missing",
    )
    return histogram, invariants, tuple(representatives)


def fmt(counter):
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main():
    masks, cards, mask_digest = build_tables()
    cardinality_formula_audit(cards)
    ratio_mask_signatures(masks)
    cayley = cayley_audit()
    words, state_words, grammar_digest = grammar()
    scalar = scalar_census(words, cards)

    feasible_histogram = Counter()
    maximum_histogram = Counter()
    bank_size_histogram = Counter()
    maximizing_histogram = Counter()
    layer_histogram = Counter()
    missing_sheet_histogram = Counter()
    tie_histogram = Counter()
    flip_histogram = Counter()
    owner_digest = sha256()
    owner_profiles = Counter()
    total_masks = 0
    maximum_bank = 0
    for support, word, capacities in scalar["survivors"]:
        rows = []
        for owner in support:
            row = owner_local(support, word, owner, masks)
            feasible, maximum, bank_size, maximizing, layers, reachable = row
            rows.append(row)
            proof_bound = structural_owner_bound(support, word, owner)
            require(maximum <= proof_bound, "literal DP exceeds structural proof bound")
            require(proof_bound == (22 if word[support.index(owner)] == 5 else 21),
                    "owner-type bound")
            maximum_histogram[maximum] += 1
            bank_size_histogram[bank_size] += 1
            maximizing_histogram[maximizing] += 1
            layer_histogram[layers] += 1
            total_masks += bank_size
            maximum_bank = max(maximum_bank, bank_size)
            for mask in reachable:
                if mask.bit_count() == maximum:
                    missing_sheet_histogram[FULL ^ mask] += 1
            owner_digest.update(bytes((owner, int(feasible), maximum)))
            owner_digest.update(bank_size.to_bytes(8, "little"))
            owner_digest.update(maximizing.to_bytes(8, "little"))
            for mask in sorted(reachable):
                owner_digest.update(mask.to_bytes(4, "little"))
        feasible_count = sum(row[0] for row in rows)
        feasible_histogram[feasible_count] += 1
        owner_profiles[tuple(row[1] for row in rows)] += 1
        ties, flips, _, _, _, _ = tournament_fingerprint(rows, capacities)
        tie_histogram[ties] += 1
        flip_histogram[flips] += 1

    orbit_histogram, orbit_invariants, orbit_representatives = support_orbits(
        scalar["survivors"]
    )

    owner_hash = owner_digest.hexdigest()
    require(
        mask_digest == "741748b977fd90f0b506a15780e33d87bb04fdd60ae4f844ea1e40349ff8c47d",
        "mask digest",
    )
    require(
        grammar_digest == "7ae50439ddbd7e09d37516d067fe20f35d1f36b7830c25dc7156c900c6fde62f",
        "grammar digest",
    )
    require(
        scalar["digest"] == "ad266f55f820615eb2f7b4e323b6599024842c56352ff65c1b6dcdd117c250f9",
        "scalar digest",
    )
    require(
        owner_hash == "a9064e057fcd7169395a28a3c411476f7b8178ec6137fbf74235aaee6db1f85c",
        "owner digest",
    )
    require(feasible_histogram == Counter({0: 36}), "owner feasible contexts")
    require(maximum_histogram == Counter({21: 144, 22: 72}), "owner maxima")
    require(
        bank_size_histogram
        == Counter({45200: 24, 48380: 24, 48540: 24,
                    133390: 48, 140330: 48, 141430: 48}),
        "reachable-bank sizes",
    )
    require(
        maximizing_histogram == Counter({80: 24, 90: 144, 100: 24, 140: 24}),
        "maximizing-mask counts",
    )
    require(total_masks == 23_338_080 and maximum_bank == 141_430, "mask totals")
    require(len(layer_histogram) == 108, "layer profile count")
    require(len(missing_sheet_histogram) == 3813, "missing-mask count")
    require(sum(missing_sheet_histogram.values()) == 20_640, "maximizer total")
    require(max(missing_sheet_histogram.values()) == 14, "missing-mask multiplicity")
    require(len(owner_profiles) == 13 and sum(owner_profiles.values()) == 36,
            "owner maximum-vector profiles")
    require(
        all(vector.count(21) == 4 and vector.count(22) == 2 for vector in owner_profiles),
        "owner maximum-vector shape",
    )
    require(tie_histogram == Counter({7: 36}), "tournament tie edges")
    require(
        flip_histogram == Counter({0: 1, 1: 2, 2: 5, 3: 2, 4: 16,
                                   5: 2, 6: 5, 7: 2, 8: 1}),
        "tournament flip edges",
    )
    require(orbit_histogram == Counter({12: 3}), "multiplication orbits")
    require(orbit_invariants == Counter({2: 1, 5: 1, 6: 1}), "orbit ratios")

    print("scale-twenty-five AP-centred Hamming-six prime-square obstruction certificate")
    print("orders", ORDERS, "unit counts", tuple(len(UNITS[d]) for d in ORDERS))
    print("generator powers", GENERATOR_POWERS)
    print("order-five forbidden provider/owner ratios", tuple(sorted(FORBIDDEN_RATIO_5)))
    print("inverse forbidden owner multipliers", tuple(sorted(FORBIDDEN_OWNER_MULTIPLIER)))
    print(
        "forbidden Cayley audit directed-arcs", cayley[0],
        "symmetric-edges", cayley[1], "reciprocal-pairs", cayley[2],
        "component-sizes", cayley[3], "graph K6+K6",
    )
    print("hereditary words", len(words), "state words/support", state_words,
          "raw labelled states", 924 * state_words)
    print("grammar SHA256", grammar_digest)
    print("mask SHA256", mask_digest)
    print("scalar feasible-owner histogram", fmt(scalar["feasible_owner_histogram"]))
    print("scalar supports-by-survivor-count", fmt(scalar["support_histogram"]))
    print("scalar survivors", len(scalar["survivors"]), "multiplicities",
          fmt(scalar["multiplicity_histogram"]))
    print("literal unit words across scalar survivors", 36 * 4**2 * 20**4)
    print("scalar SHA256", scalar["digest"])
    print("multiplication orbit-size histogram", fmt(orbit_histogram))
    print("orbit representative ratio invariants", fmt(orbit_invariants))
    for providers5, support in orbit_representatives:
        print("orbit representative B5", providers5, "support", support,
              "C25", tuple(x for x in support if x not in providers5))
    print("proof bounds every order-five owner union <=22; every order-twenty-five owner union <=21")
    print("Kakeya/toothpick quotient pi: Z/25 -> Z/5 turns D5 masks into full parallel fibres and D25 masks into rigid three/four-point needles; every survivor pays four forced needle/fibre incidences")
    print("owner-local feasible-context histogram", fmt(feasible_histogram))
    print("owner maximum-union histogram", fmt(maximum_histogram))
    print("owner reachable-bank-size histogram", fmt(bank_size_histogram))
    print("owner maximizing-mask-count histogram", fmt(maximizing_histogram))
    print("owner maximum-vector/context histogram", fmt(owner_profiles))
    print("owner rows", 6 * len(scalar["survivors"]), "total reachable masks",
          total_masks, "maximum bank", maximum_bank)
    print("distinct layer-size profiles", len(layer_histogram))
    print("distinct maximum missing masks", len(missing_sheet_histogram))
    print("total maximizing masks", sum(missing_sheet_histogram.values()),
          "largest missing-mask multiplicity", max(missing_sheet_histogram.values()))
    print("owner SHA256", owner_hash)
    print("tournament observable (feasible,max-union,capacity), lex switch, coordinate ties")
    print("tournament all transitive: score histogram 0,1,2,3,4,5; cycles 0; singleton SCCs; Hamiltonian paths 1")
    print("tournament tie-edge histogram", fmt(tie_histogram))
    print("tournament flip-edge histogram", fmt(flip_histogram))
    print("alternate carrier audit quadratic classes and the two D5 obligations preserve the scalar K6+K6 obstruction; quotient fibres preserve the four-incidence overlap tax; owner obligations with literal masks preserve terminal feasibility")
    print("lost-coordinate ledger a completed tournament loses the absolute threshold; providers alone, isolated sheets, gaps, wall events, unlabelled residues, Fano points, and chi7 colours destroy shared-unit, owner, or fibre incidence")
    print("frontier verdict scale 25 structurally empty; next legal untreated common scale 27 because THM-860 excludes multiples of 13")


if __name__ == "__main__":
    main()
