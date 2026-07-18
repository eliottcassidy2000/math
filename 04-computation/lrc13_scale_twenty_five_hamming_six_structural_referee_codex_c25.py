#!/usr/bin/env python3
"""Independent exact referee for the scale-25 Hamming-six obstruction.

This standard-library-only program reconstructs every CRT mask by literal
search.  It then checks, separately:

* the complete hereditary divisor grammar and its literal-unit weight;
* the scalar reduction from 437,052 support/order rows to 36 rows;
* the owner-normalized five-coset signatures used by the structural proof;
* the resulting 22/21 owner bounds without a Cartesian unit search; and
* the exact projected union banks on all 216 surviving owner obligations.

Multiplication orbits and completed owner tournaments are telemetry only.  In
particular, neither quotient is used to skip a labelled row or owner.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod


P = 13
C = 25
LABELS = tuple(range(1, P))
ORDERS = (1, 5, 25)
UNITS = {
    order: tuple(
        unit
        for unit in range(order)
        if (order == 1 and unit == 0)
        or (order > 1 and gcd(unit, order) == 1)
    )
    for order in ORDERS
}
FULL = (1 << C) - 1
QUADRATIC = frozenset((1, 3, 4, 9, 10, 12))
NONQUADRATIC = frozenset(LABELS) - QUADRATIC
ZERO_RATIOS = frozenset((4, 9, 12))
FORBIDDEN_MULTIPLIERS = frozenset(pow(ratio, -1, P) for ratio in ZERO_RATIOS)


def require(condition, message):
    """Optimization-stable assertion."""
    if not condition:
        raise RuntimeError(message)


def centered(value, modulus):
    """Return the representative in (-modulus/2, modulus/2]."""
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def literal_crt_base(label, order, unit):
    """Find x = order*label (mod 13), x = unit (mod order) literally."""
    candidates = tuple(
        value
        for value in range(P * order)
        if value % P == order * label % P and value % order == unit % order
    )
    require(len(candidates) == 1, "CRT representative is not unique")
    return candidates[0]


def local_mask(label, order, unit, owner):
    """Construct the exact strict-left/closed-right local sheet mask."""
    base = literal_crt_base(label, order, unit)
    owner_inverse = pow(owner, -1, P)
    mask = 0
    for sheet in range(C):
        value = centered(base * (owner_inverse + P * sheet), P * order)
        if -order < value <= order:
            mask |= 1 << sheet
    return mask


def analytic_cardinality(label, order, owner):
    """Count one order-period and repeat it C/order times."""
    ratio = label * pow(owner, -1, P) % P
    target = order * ratio % P
    one_period = sum(
        value % P == target for value in range(-order + 1, order + 1)
    )
    return (C // order) * one_period


def build_mask_bank():
    masks = {}
    cards = {}
    digest = sha256()
    for label in LABELS:
        for order in ORDERS:
            for owner in LABELS:
                expected = analytic_cardinality(label, order, owner)
                option_cards = set()
                for unit in UNITS[order]:
                    mask = local_mask(label, order, unit, owner)
                    option_cards.add(mask.bit_count())
                    masks[label, order, unit, owner] = mask
                    digest.update(bytes((label, order, unit, owner)))
                    digest.update(mask.to_bytes(4, "little"))
                require(option_cards == {expected}, "unit-dependent cardinality")
                cards[label, order, owner] = expected
    return masks, cards, digest.hexdigest()


def hereditary(word):
    return all(
        lcm(*(word[index] for index in range(6) if index != omitted)) == C
        for omitted in range(6)
    )


def build_grammar():
    words = []
    weighted_states = 0
    digest = sha256()
    for word in product(ORDERS, repeat=6):
        by_lcm = hereditary(word)
        by_prime_power = word.count(25) >= 2
        require(by_lcm == by_prime_power, "hereditary/prime-power mismatch")
        if not by_lcm:
            continue
        words.append(word)
        fibre = prod(len(UNITS[order]) for order in word)
        weighted_states += fibre
        digest.update(bytes(word))
        digest.update(fibre.to_bytes(8, "little"))
    return tuple(words), weighted_states, digest.hexdigest()


def ratio_cardinality_table(cards):
    table = {}
    for order in ORDERS:
        row = tuple(cards[ratio, order, 1] for ratio in LABELS)
        for owner in LABELS:
            translated = tuple(
                cards[ratio * owner % P, order, owner] for ratio in LABELS
            )
            require(translated == row, "cardinality is not ratio-invariant")
        table[order] = row
    require(table[1] == (25,) + (0,) * 11, "bad order-one row")
    require(
        table[5] == tuple(0 if ratio in ZERO_RATIOS else 5 for ratio in LABELS),
        "bad order-five row",
    )
    require(
        table[25] == tuple(3 if ratio == 12 else 4 for ratio in LABELS),
        "bad order-twenty-five row",
    )
    return table


def scalar_bank(words, cards):
    survivors = []
    scalar_digest = sha256()
    d1_witness_maximum = 0
    formula_checks = 0

    for support in combinations(LABELS, 6):
        for word in words:
            capacities = tuple(
                sum(
                    cards[label, order, owner]
                    for label, order in zip(support, word)
                )
                for owner in support
            )

            if 1 in word:
                owner_index = word.index(25)
                d1_witness_maximum = max(
                    d1_witness_maximum, capacities[owner_index]
                )
            else:
                count_25 = word.count(25)
                for owner, capacity in zip(support, capacities):
                    z_owner = sum(
                        order == 5
                        and label * pow(owner, -1, P) % P in ZERO_RATIOS
                        for label, order in zip(support, word)
                    )
                    delta_owner = sum(
                        order == 25
                        and label * pow(owner, -1, P) % P == 12
                        for label, order in zip(support, word)
                    )
                    require(
                        capacity == 30 - count_25 - 5 * z_owner - delta_owner,
                        "scalar formula mismatch",
                    )
                    formula_checks += 1

            if min(capacities) < C:
                continue
            survivors.append((support, word, capacities))
            scalar_digest.update(bytes(support))
            scalar_digest.update(bytes(word))
            scalar_digest.update(bytes(capacities))

    require(d1_witness_maximum == 23, "order-one witness bound changed")

    expected = set()
    for square_label in QUADRATIC:
        for nonsquare_label in NONQUADRATIC:
            order_five = frozenset((square_label, nonsquare_label))
            forbidden = {
                multiplier * label % P
                for label in order_five
                for multiplier in FORBIDDEN_MULTIPLIERS
            }
            require(len(forbidden) == 6, "forbidden triples overlap")
            support = tuple(sorted(set(LABELS) - forbidden))
            require(order_five <= set(support), "order-five label was forbidden")
            word = tuple(5 if label in order_five else 25 for label in support)
            capacities = tuple(
                sum(
                    cards[label, order, owner]
                    for label, order in zip(support, word)
                )
                for owner in support
            )
            expected.add((support, word, capacities))

    require(set(survivors) == expected, "symbolic survivor classification mismatch")
    require(len(survivors) == 36, "unexpected scalar survivor count")
    require(
        all(tuple(word.count(order) for order in ORDERS) == (0, 2, 4)
            for _, word, _ in survivors),
        "unexpected survivor multiplicity",
    )

    return tuple(survivors), d1_witness_maximum, formula_checks, scalar_digest.hexdigest()


def multiplication_orbits(survivors):
    rows = {
        (
            tuple(label for label, order in zip(support, word) if order == 5),
            support,
        )
        for support, word, _ in survivors
    }
    remaining = set(rows)
    representatives = []
    sizes = Counter()
    while remaining:
        seed = min(remaining)
        orbit = {
            (
                tuple(sorted(multiplier * label % P for label in seed[0])),
                tuple(sorted(multiplier * label % P for label in seed[1])),
            )
            for multiplier in LABELS
        }
        require(orbit <= rows, "multiplication orbit left scalar bank")
        representatives.append(min(orbit))
        sizes[len(orbit)] += 1
        remaining -= orbit
    expected = [
        ((1, 2), (1, 2, 4, 5, 8, 9)),
        ((1, 5), (1, 4, 5, 6, 7, 9)),
        ((1, 6), (1, 2, 4, 6, 9, 11)),
    ]
    require(sorted(representatives) == expected, "orbit representatives changed")
    require(sizes == Counter({12: 3}), "orbit sizes changed")
    return tuple(sorted(representatives)), sizes


def occupied_cosets(mask):
    return tuple(
        sum((mask >> sheet) & 1 for sheet in range(residue, C, 5))
        for residue in range(5)
    )


def five_coset_signature_audit(masks):
    """Audit signatures after naming each owner's self coset Q_3."""
    self_cosets = []
    signature_digest = sha256()
    checks = 0

    for owner in LABELS:
        self_options = {
            masks[owner, 5, unit, owner] for unit in UNITS[5]
        }
        require(len(self_options) == 1, "self order-five mask depends on unit")
        self_mask = next(iter(self_options))
        self_signature = occupied_cosets(self_mask)
        require(sorted(self_signature) == [0, 0, 0, 0, 5], "self mask not a coset")
        self_coset = self_signature.index(5)
        self_cosets.append(self_coset)
        other_cosets = set(range(5)) - {self_coset}

        for ratio in LABELS:
            label = ratio * owner % P
            order_five_masks = {
                masks[label, 5, unit, owner] for unit in UNITS[5]
            }
            if ratio in ZERO_RATIOS:
                require(order_five_masks == {0}, "forbidden order-five mask nonempty")
            elif ratio == 1:
                require(order_five_masks == {self_mask}, "self coset changed")
            else:
                signatures = [occupied_cosets(mask) for mask in order_five_masks]
                require(len(order_five_masks) == 4, "moving order-five fibre collapsed")
                require(
                    {signature.index(5) for signature in signatures} == other_cosets,
                    "moving order-five masks do not traverse the other four cosets",
                )
                require(
                    all(sorted(signature) == [0, 0, 0, 0, 5]
                        for signature in signatures),
                    "moving order-five mask is not a complete coset",
                )

            for unit in UNITS[25]:
                mask = masks[label, 25, unit, owner]
                signature = occupied_cosets(mask)
                occupied = {index for index, count in enumerate(signature) if count}
                require(max(signature) <= 1, "order-25 mask repeats a five-coset")
                if ratio in (4, 9):
                    require(
                        occupied == other_cosets,
                        "ratio 4/9 does not occupy all four nonself cosets",
                    )
                elif ratio == 12:
                    require(
                        len(occupied) == 3 and self_coset not in occupied,
                        "antipodal ratio has wrong three-coset signature",
                    )
                else:
                    require(
                        len(occupied) == 4 and self_coset in occupied,
                        "generic order-25 ratio has wrong four-coset signature",
                    )
                signature_digest.update(bytes((owner, ratio, unit, self_coset)))
                signature_digest.update(bytes(signature))
                checks += 1

    # The raw sheet residue of the self coset is owner-dependent.  The theorem
    # may call it Q_3 only after this owner-local relabelling.
    require(tuple(self_cosets) == (3, 1, 2, 0, 4, 3, 1, 0, 4, 2, 3, 1),
            "raw self-coset map changed")
    return tuple(self_cosets), checks, signature_digest.hexdigest()


def structural_owner_bound(support, word, owner, masks):
    owner_index = support.index(owner)
    owner_order = word[owner_index]
    self_mask = masks[owner, 5, UNITS[5][0], owner]
    self_coset = occupied_cosets(self_mask).index(5)

    five_indices = [index for index, order in enumerate(word) if order == 5]
    twenty_five_indices = [index for index, order in enumerate(word) if order == 25]
    ratios_25 = {
        support[index] * pow(owner, -1, P) % P
        for index in twenty_five_indices
    }

    five_options = [
        {
            masks[support[index], 5, unit, owner]
            for unit in UNITS[5]
        }
        for index in five_indices
    ]

    if owner_order == 5:
        require(
            ratios_25 & QUADRATIC == {4, 9}
            and ratios_25 & NONQUADRATIC <= NONQUADRATIC
            and len(ratios_25 & NONQUADRATIC) == 2,
            "order-five owner ratio pattern changed",
        )
        self_index = five_indices.index(owner_index)
        moving_index = 1 - self_index
        require(five_options[self_index] == {self_mask}, "bad self fibre")
        for moving_mask in five_options[moving_index]:
            base_union = self_mask | moving_mask
            require(base_union.bit_count() == 10, "order-five cosets overlap")
            for index in twenty_five_indices:
                ratio = support[index] * pow(owner, -1, P) % P
                for unit in UNITS[25]:
                    mask = masks[support[index], 25, unit, owner]
                    if ratio in (4, 9):
                        require(mask & moving_mask, "ratio 4/9 misses moving coset")
                    else:
                        require(ratio in NONQUADRATIC, "unexpected residual ratio")
                        require(mask & self_mask, "nonresidue misses self coset")
        return 22

    require(owner_order == 25, "unexpected owner order")
    require(
        ratios_25 & QUADRATIC == {1, 12}
        and len(ratios_25 & NONQUADRATIC) == 2,
        "order-25 owner ratio pattern changed",
    )
    for first_mask in five_options[0]:
        for second_mask in five_options[1]:
            base_union = first_mask | second_mask
            if first_mask == second_mask:
                require(base_union.bit_count() == 5, "equal cosets have wrong size")
                require(5 + 15 == 20, "coincident-coset bound changed")
            else:
                require(base_union.bit_count() == 10, "distinct cosets overlap")
                for index in twenty_five_indices:
                    for unit in UNITS[25]:
                        mask = masks[support[index], 25, unit, owner]
                        require(mask & base_union, "order-25 mask misses two cosets")
                require(10 + 15 - 4 == 21, "distinct-coset bound changed")
    return 21


def strongly_connected_component_sizes(out):
    reach = list(out)
    for vertex in range(6):
        reach[vertex] |= 1 << vertex
    for middle in range(6):
        for source in range(6):
            if (reach[source] >> middle) & 1:
                reach[source] |= reach[middle]
    remaining = set(range(6))
    sizes = []
    while remaining:
        root = min(remaining)
        component = {
            vertex
            for vertex in remaining
            if (reach[root] >> vertex) & 1 and (reach[vertex] >> root) & 1
        }
        require(component, "empty tournament component")
        sizes.append(len(component))
        remaining -= component
    return tuple(sorted(sizes))


def tournament_fingerprint(summaries):
    """Complete the lexicographic owner preorder by coordinate order."""
    out = [0] * 6
    ties = 0
    flips = 0
    for left, right in combinations(range(6), 2):
        if summaries[left] == summaries[right]:
            winner = left
            ties += 1
        elif summaries[right] > summaries[left]:
            winner = right
            flips += 1
        else:
            winner = left
        loser = left + right - winner
        out[winner] |= 1 << loser

    scores = tuple(sorted(mask.bit_count() for mask in out))
    triangles = 0
    for a, b, c in combinations(range(6), 3):
        triangles += int(
            ((out[a] >> b) & 1 and (out[b] >> c) & 1 and (out[c] >> a) & 1)
            or ((out[a] >> c) & 1 and (out[c] >> b) & 1 and (out[b] >> a) & 1)
        )

    paths = [[0] * 6 for _ in range(1 << 6)]
    for last in range(6):
        paths[1 << last][last] = 1
    for mask in range(1, 1 << 6):
        for last in range(6):
            if not (mask >> last) & 1:
                continue
            previous_mask = mask ^ (1 << last)
            for previous in range(6):
                if ((previous_mask >> previous) & 1
                        and (out[previous] >> last) & 1):
                    paths[mask][last] += paths[previous_mask][previous]

    return (
        ties,
        flips,
        scores,
        triangles,
        strongly_connected_component_sizes(out),
        sum(paths[-1]),
    )


def exact_owner_banks(survivors, masks):
    maximum_histogram = Counter()
    maximum_by_order = Counter()
    reachable_histogram = Counter()
    owner_profile_histogram = Counter()
    structural_bound_histogram = Counter()
    tie_histogram = Counter()
    flip_histogram = Counter()
    tournament_fingerprints = Counter()
    feasible_owner_rows = 0
    reachable_total = 0
    digest = sha256()

    for support, word, capacities in survivors:
        summaries = []
        maxima = []
        for owner_index, owner in enumerate(support):
            reachable = frozenset((0,))
            for label, order in zip(support, word):
                options = frozenset(
                    masks[label, order, unit, owner] for unit in UNITS[order]
                )
                reachable = frozenset(
                    partial | option
                    for partial in reachable
                    for option in options
                )

            maximum = max(mask.bit_count() for mask in reachable)
            feasible = FULL in reachable
            structural_bound = structural_owner_bound(
                support, word, owner, masks
            )
            require(maximum <= structural_bound, "structural bound violated")
            require(not feasible, "scale-25 owner counterexample")

            owner_order = word[owner_index]
            maximum_histogram[maximum] += 1
            maximum_by_order[owner_order, maximum] += 1
            reachable_histogram[len(reachable)] += 1
            structural_bound_histogram[owner_order, structural_bound] += 1
            feasible_owner_rows += int(feasible)
            reachable_total += len(reachable)
            maxima.append(maximum)
            summaries.append((maximum, capacities[owner_index], len(reachable)))

            digest.update(bytes(support))
            digest.update(bytes(word))
            digest.update(bytes((owner, owner_order, maximum, int(feasible))))
            digest.update(len(reachable).to_bytes(4, "little"))
            for mask in sorted(reachable):
                digest.update(mask.to_bytes(4, "little"))

        owner_profile_histogram[tuple(maxima)] += 1
        fingerprint = tournament_fingerprint(tuple(summaries))
        ties, flips = fingerprint[:2]
        tie_histogram[ties] += 1
        flip_histogram[flips] += 1
        tournament_fingerprints[fingerprint[2:]] += 1

    require(feasible_owner_rows == 0, "feasible owner row found")
    require(maximum_histogram == Counter({21: 144, 22: 72}),
            "owner maximum histogram changed")
    require(maximum_by_order == Counter({(25, 21): 144, (5, 22): 72}),
            "owner/order maxima changed")
    require(structural_bound_histogram == maximum_by_order,
            "structural bound is not exact")
    require(sum(owner_profile_histogram.values()) == 36,
            "owner-profile census changed")
    require(
        tournament_fingerprints
        == Counter({((0, 1, 2, 3, 4, 5), 0, (1, 1, 1, 1, 1, 1), 1): 36}),
        "completed tournament should be transitive",
    )

    return {
        "maximum_histogram": maximum_histogram,
        "maximum_by_order": maximum_by_order,
        "reachable_histogram": reachable_histogram,
        "reachable_total": reachable_total,
        "owner_profile_count": len(owner_profile_histogram),
        "tie_histogram": tie_histogram,
        "flip_histogram": flip_histogram,
        "tournament_fingerprints": tournament_fingerprints,
        "digest": digest.hexdigest(),
    }


def fmt(counter):
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main():
    masks, cards, mask_digest = build_mask_bank()
    ratio_table = ratio_cardinality_table(cards)
    words, weighted_states, grammar_digest = build_grammar()
    survivors, d1_maximum, formula_checks, scalar_digest = scalar_bank(words, cards)
    representatives, orbit_sizes = multiplication_orbits(survivors)
    self_cosets, signature_checks, signature_digest = five_coset_signature_audit(masks)
    exact = exact_owner_banks(survivors, masks)

    capacity_histogram = Counter(
        capacity for _, _, capacities in survivors for capacity in capacities
    )
    multiplicity_histogram = Counter(
        tuple(word.count(order) for order in ORDERS)
        for _, word, _ in survivors
    )

    print("scale=25 p=13 hamming=6 structural-referee=standard-library-literal-crt")
    print(f"orders={ORDERS} unit-counts={tuple(len(UNITS[d]) for d in ORDERS)}")
    print(
        f"grammar words={len(words)} weighted-states/support={weighted_states} "
        f"unquotiented={len(tuple(combinations(LABELS, 6))) * weighted_states}"
    )
    print(f"ratio-cardinality D1={ratio_table[1]}")
    print(f"ratio-cardinality D5={ratio_table[5]}")
    print(f"ratio-cardinality D25={ratio_table[25]}")
    print(
        f"scalar formula-checks={formula_checks} d1-witness-max={d1_maximum} "
        f"survivors={len(survivors)} multiplicities={fmt(multiplicity_histogram)}"
    )
    print(f"scalar owner-capacities={fmt(capacity_histogram)}")
    print(f"multiplication-orbits={fmt(orbit_sizes)}")
    for order_five, support in representatives:
        print(f"orbit-representative B={order_five} support={support}")
    print(f"raw-self-cosets-by-owner={self_cosets}")
    print(
        f"owner-normalized-signatures={signature_checks} "
        "Q3=name-of-self-coset-not-a-global-sheet-residue"
    )
    print(f"exact maximum-union={fmt(exact['maximum_histogram'])}")
    print(f"exact maximum-by-owner-order={fmt(exact['maximum_by_order'])}")
    print(f"reachable-counts={fmt(exact['reachable_histogram'])}")
    print(
        f"reachable-total={exact['reachable_total']} "
        f"owner-profile-count={exact['owner_profile_count']} feasible-owner-rows=0"
    )
    print(
        f"completed-owner-tournament ties={fmt(exact['tie_histogram'])} "
        f"flips={fmt(exact['flip_histogram'])} scores=(0,1,2,3,4,5) "
        "triangles=0 scc=(1,1,1,1,1,1) hamiltonian-paths=1"
    )
    print(f"mask-bank-sha256={mask_digest}")
    print(f"grammar-sha256={grammar_digest}")
    print(f"scalar-bank-sha256={scalar_digest}")
    print(f"five-coset-signature-sha256={signature_digest}")
    print(f"owner-bank-sha256={exact['digest']}")
    print("RESULT scale-25 scalar bank has 36 rows and every owner misses >=3 sheets")


if __name__ == "__main__":
    main()
