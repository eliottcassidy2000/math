#!/usr/bin/env python3
"""Exact scale-27 AP-centred Hamming-six owner-deficit certificate.

The calculation is standard-library only.  It reconstructs every CRT mask,
enumerates all labelled support/order rows without a symmetry quotient, and
then computes every owner-local reachable union.  Multiplication orbits,
tournaments, and residue signatures are telemetry; the terminal certificate
is the absence of a row whose six owner projections are all feasible.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod


P = 13
C = 27
ORDERS = (1, 3, 9, 27)
UNITS = {
    order: (
        (0,)
        if order == 1
        else tuple(unit for unit in range(1, order + 1) if gcd(unit, order) == 1)
    )
    for order in ORDERS
}
FULL_MASK = (1 << C) - 1
MASK64 = (1 << 64) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base_algebraic(label: int, order: int, unit: int) -> int:
    if order == 1:
        return label
    step = ((unit - order * label) * pow(P, -1, order)) % order
    return (order * label + P * step) % (P * order)


def crt_base_literal(label: int, order: int, unit: int) -> int:
    for value in range(P * order):
        if value % P == order * label % P and value % order == unit % order:
            return value
    raise RuntimeError("literal CRT search failed")


def local_mask(label: int, order: int, unit: int, owner: int) -> int:
    base = crt_base_algebraic(label, order, unit)
    owner_inverse = pow(owner, -1, P)
    result = 0
    for sheet in range(C):
        value = centered(base * (owner_inverse + P * sheet), P * order)
        if -order < value <= order:
            result |= 1 << sheet
    return result


def analytic_cardinality(label: int, order: int, owner: int) -> int:
    ratio = label * pow(owner, -1, P) % P
    target = order * ratio % P
    period = sum(value % P == target for value in range(-order + 1, order + 1))
    return (C // order) * period


def build_tables():
    masks = {}
    distinct_masks = {}
    cards = {}
    mask_digest = sha256()
    residue_digest = sha256()
    choice_histogram = Counter()
    residue_signature_histogram = Counter()

    for label in range(1, P):
        for order in ORDERS:
            for owner in range(1, P):
                card = analytic_cardinality(label, order, owner)
                cards[label, order, owner] = card
                choices = []
                for unit in UNITS[order]:
                    require(
                        crt_base_algebraic(label, order, unit)
                        == crt_base_literal(label, order, unit),
                        "algebraic/literal CRT mismatch",
                    )
                    mask = local_mask(label, order, unit, owner)
                    require(mask.bit_count() == card, "mask/cardinality mismatch")
                    masks[label, order, unit, owner] = mask
                    choices.append(mask)
                    mask_digest.update(bytes((label, order, unit, owner)))
                    mask_digest.update(mask.to_bytes(4, "little"))

                    mod_three = tuple(
                        sum((mask >> sheet) & 1 for sheet in range(residue, C, 3))
                        for residue in range(3)
                    )
                    mod_nine = tuple(
                        sum((mask >> sheet) & 1 for sheet in range(residue, C, 9))
                        for residue in range(9)
                    )
                    residue_signature_histogram[order, card, mod_three, mod_nine] += 1
                    residue_digest.update(bytes((label, order, unit, owner)))
                    residue_digest.update(bytes(mod_three + mod_nine))

                unique = tuple(sorted(set(choices)))
                distinct_masks[label, order, owner] = unique
                choice_histogram[order, card, len(unique)] += 1

    return {
        "masks": masks,
        "distinct_masks": distinct_masks,
        "cards": cards,
        "mask_digest": mask_digest.hexdigest(),
        "residue_digest": residue_digest.hexdigest(),
        "choice_histogram": choice_histogram,
        "residue_signature_histogram": residue_signature_histogram,
    }


def hereditary(order_word: tuple[int, ...]) -> bool:
    return all(
        lcm(*(order_word[index] for index in range(6) if index != omitted)) == C
        for omitted in range(6)
    )


def hereditary_prime_power(order_word: tuple[int, ...]) -> bool:
    return order_word.count(27) >= 2


def build_order_words():
    words = []
    digest = sha256()
    for word in product(ORDERS, repeat=6):
        via_lcm = hereditary(word)
        require(via_lcm == hereditary_prime_power(word), "hereditary grammar mismatch")
        if via_lcm:
            words.append(word)
            digest.update(bytes(word))
    return words, digest.hexdigest()


def scalar_bank(order_words, cards):
    bank = []
    support_histogram = Counter()
    feasible_owner_histogram = Counter()
    multiplicity_histogram = Counter()
    capacity_histogram = Counter()
    minimum_slack = Counter()
    capacity_vectors = set()
    digest = sha256()

    for support in combinations(range(1, P), 6):
        survivors = 0
        for word in order_words:
            capacities = tuple(
                sum(cards[label, order, owner] for label, order in zip(support, word))
                for owner in support
            )
            feasible_count = sum(value >= C for value in capacities)
            feasible_owner_histogram[feasible_count] += 1
            if feasible_count != 6:
                continue

            survivors += 1
            bank.append((support, word, capacities))
            multiplicity_histogram[
                tuple(word.count(order) for order in ORDERS)
            ] += 1
            capacity_histogram.update(capacities)
            minimum_slack[min(capacities) - C] += 1
            capacity_vectors.add(capacities)
            digest.update(bytes(support))
            digest.update(bytes(word))
            digest.update(bytes(capacities))
        support_histogram[survivors] += 1

    return {
        "bank": bank,
        "support_histogram": support_histogram,
        "feasible_owner_histogram": feasible_owner_histogram,
        "multiplicity_histogram": multiplicity_histogram,
        "capacity_histogram": capacity_histogram,
        "minimum_slack": minimum_slack,
        "capacity_vectors": capacity_vectors,
        "digest": digest.hexdigest(),
    }


def mix64(value: int) -> int:
    value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9 & MASK64
    value = (value ^ (value >> 27)) * 0x94D049BB133111EB & MASK64
    return value ^ (value >> 31)


def owner_local(support, word, owner, distinct_masks):
    reachable = frozenset({0})
    for label, order in zip(support, word):
        choices = distinct_masks[label, order, owner]
        reachable = frozenset(partial | choice for partial in reachable for choice in choices)

    maximum = max(mask.bit_count() for mask in reachable)
    feasible = FULL_MASK in reachable
    require(feasible == (maximum == C), "owner feasibility/maximum mismatch")
    mask_sum = sum(reachable) & MASK64
    mask_xor = 0
    for mask in reachable:
        mask_xor ^= mix64(mask)
    return feasible, maximum, len(reachable), mask_sum, mask_xor


def strongly_connected_sizes(adjacency):
    reach = [row[:] for row in adjacency]
    for vertex in range(6):
        reach[vertex][vertex] = True
    for middle in range(6):
        for left in range(6):
            if reach[left][middle]:
                for right in range(6):
                    reach[left][right] = reach[left][right] or reach[middle][right]

    unseen = set(range(6))
    sizes = []
    while unseen:
        seed = min(unseen)
        component = {vertex for vertex in unseen if reach[seed][vertex] and reach[vertex][seed]}
        sizes.append(len(component))
        unseen -= component
    return tuple(sorted(sizes))


def hamiltonian_path_count(adjacency):
    counts = [[0] * 6 for _ in range(1 << 6)]
    for vertex in range(6):
        counts[1 << vertex][vertex] = 1
    for subset in range(1 << 6):
        for last in range(6):
            ways = counts[subset][last]
            if not ways:
                continue
            for following in range(6):
                if not (subset >> following) & 1 and adjacency[last][following]:
                    counts[subset | (1 << following)][following] += ways
    return sum(counts[(1 << 6) - 1])


def tournament_summary(capacities, local_rows):
    adjacency = [[False] * 6 for _ in range(6)]
    ties = 0
    flips = 0
    for left in range(6):
        for right in range(left + 1, 6):
            left_key = (
                local_rows[left][0],
                local_rows[left][1],
                capacities[left],
                local_rows[left][2],
            )
            right_key = (
                local_rows[right][0],
                local_rows[right][1],
                capacities[right],
                local_rows[right][2],
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
    for left, middle, right in combinations(range(6), 3):
        triangles += int(
            (
                adjacency[left][middle]
                and adjacency[middle][right]
                and adjacency[right][left]
            )
            or (
                adjacency[left][right]
                and adjacency[right][middle]
                and adjacency[middle][left]
            )
        )
    return (
        scores,
        triangles,
        strongly_connected_sizes(adjacency),
        hamiltonian_path_count(adjacency),
        ties,
        flips,
    )


def multiply_row(row, multiplier):
    return tuple(sorted((multiplier * label % P, order) for label, order in row))


def multiplication_orbits(rows):
    remaining = set(rows)
    histogram = Counter()
    while remaining:
        row = min(remaining)
        orbit = {multiply_row(row, multiplier) for multiplier in range(1, P)}
        require(orbit <= rows, "multiplication invariance failed")
        remaining -= orbit
        histogram[len(orbit)] += 1
    return histogram


def fmt(counter: Counter) -> str:
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main() -> None:
    tables = build_tables()
    order_words, order_digest = build_order_words()
    state_words_per_support = sum(
        prod(len(UNITS[order]) for order in word) for word in order_words
    )
    scalar = scalar_bank(order_words, tables["cards"])
    bank = scalar["bank"]

    feasible_contexts = Counter()
    maximum_union = Counter()
    minimum_owner_maximum = Counter()
    reachable_counts = Counter()
    owner_vectors = set()
    feasible_owner_sets = Counter()
    tournament_scores = Counter()
    tournament_triangles = Counter()
    tournament_sccs = Counter()
    tournament_paths = Counter()
    tournament_ties = Counter()
    tournament_flips = Counter()
    owner_digest = sha256()

    for support, word, capacities in bank:
        local_rows = []
        for owner in support:
            local = owner_local(
                support, word, owner, tables["distinct_masks"]
            )
            local_rows.append(local)
            feasible, maximum, count, mask_sum, mask_xor = local
            maximum_union[maximum] += 1
            reachable_counts[count] += 1
            owner_digest.update(bytes(support + word + (owner, int(feasible), maximum)))
            owner_digest.update(count.to_bytes(4, "little"))
            owner_digest.update(mask_sum.to_bytes(8, "little"))
            owner_digest.update(mask_xor.to_bytes(8, "little"))

        feasible_labels = tuple(
            owner for owner, local in zip(support, local_rows) if local[0]
        )
        maximum_vector = tuple(local[1] for local in local_rows)
        feasible_contexts[len(feasible_labels)] += 1
        feasible_owner_sets[feasible_labels] += 1
        minimum_owner_maximum[min(maximum_vector)] += 1
        owner_vectors.add(maximum_vector)

        scores, triangles, sccs, paths, ties, flips = tournament_summary(
            capacities, local_rows
        )
        tournament_scores[scores] += 1
        tournament_triangles[triangles] += 1
        tournament_sccs[sccs] += 1
        tournament_paths[paths] += 1
        tournament_ties[ties] += 1
        tournament_flips[flips] += 1

    scalar_rows = {
        tuple(zip(support, word)) for support, word, _ in bank
    }
    support_rows = {
        tuple((label, 0) for label in support) for support, _, _ in bank
    }

    require(len(order_words) == 1_909, "hereditary order count")
    require(state_words_per_support == 380_511_756, "literal state count")
    require(len(bank) == 450, "scalar bank count")
    require(len({support for support, _, _ in bank}) == 84, "scalar support count")
    require(
        scalar["multiplicity_histogram"]
        == Counter(
            {
                (0, 0, 3, 3): 12,
                (0, 0, 4, 2): 18,
                (0, 2, 1, 3): 60,
                (0, 2, 2, 2): 294,
                (0, 3, 0, 3): 12,
                (0, 3, 1, 2): 36,
                (0, 4, 0, 2): 18,
            }
        ),
        "scalar multiplicity profiles",
    )
    require(
        feasible_contexts == Counter({0: 336, 1: 96, 4: 18}),
        "owner deficit",
    )
    require(
        maximum_union
        == Counter({20: 120, 21: 336, 22: 192, 23: 336, 24: 528, 25: 432, 26: 588, 27: 168}),
        "maximum-union histogram",
    )
    require(max(reachable_counts) == 128_880, "reachable-bank maximum")
    require(
        sum(size * count for size, count in reachable_counts.items()) == 13_598_160,
        "reachable-bank total",
    )
    require(feasible_contexts.get(6, 0) == 0, "all-six owner row unexpectedly feasible")

    print("scale-twenty-seven AP-centred Hamming-six owner-deficit certificate")
    print("standard-library exact arithmetic; no support, order word, unit mask, or owner quotient")
    print("divisor grammar 1,3,9,27; hereditary iff at least two order-27 providers")
    print(
        f"supports 924; hereditary order words {len(order_words)}; "
        f"labelled order contexts {924 * len(order_words)}"
    )
    print(
        f"state words/support {state_words_per_support}; raw labelled states "
        f"{924 * state_words_per_support}"
    )
    print(f"mask SHA256 {tables['mask_digest']}; order SHA256 {order_digest}")
    print(f"residue-signature SHA256 {tables['residue_digest']}")
    print("distinct-mask choice histogram (order,cardinality,choices)", fmt(tables["choice_histogram"]))
    print(f"residue-signature types {len(tables['residue_signature_histogram'])}")
    print(
        f"scalar contexts {len(bank)} on {len({support for support, _, _ in bank})} supports; "
        f"capacity vectors {len(scalar['capacity_vectors'])}"
    )
    print(f"scalar-bank SHA256 {scalar['digest']}")
    print("scalar feasible-owner histogram", fmt(scalar["feasible_owner_histogram"]))
    print("contexts-per-support histogram", fmt(scalar["support_histogram"]))
    print("order multiplicity profiles", fmt(scalar["multiplicity_histogram"]))
    print("owner capacity histogram", fmt(scalar["capacity_histogram"]))
    print("minimum scalar-slack histogram", fmt(scalar["minimum_slack"]))
    print("scalar-row multiplication orbits", fmt(multiplication_orbits(scalar_rows)))
    print("support multiplication orbits", fmt(multiplication_orbits(support_rows)))
    print(f"owner-local rows {6 * len(bank)}")
    print("feasible-owner/context histogram", fmt(feasible_contexts))
    print("maximum reachable sheet-union histogram", fmt(maximum_union))
    print("minimum owner maximum/context histogram", fmt(minimum_owner_maximum))
    print(f"distinct owner maximum vectors {len(owner_vectors)}")
    print(
        f"reachable banks {6 * len(bank)}; total masks "
        f"{sum(size * count for size, count in reachable_counts.items())}; "
        f"bank-size bins {len(reachable_counts)}; maximum bank {max(reachable_counts)}"
    )
    print(f"owner-summary SHA256 {owner_digest.hexdigest()}")
    print(f"distinct feasible-owner label sets {len(feasible_owner_sets)}")
    print("owner-local all-six contexts", feasible_contexts.get(6, 0))
    print("tournament pair observable (feasible,max-union,capacity,bank-size); coordinate tie path")
    print("tournament score-sequence histogram", fmt(tournament_scores))
    print("tournament directed-triangle histogram", fmt(tournament_triangles))
    print("tournament SCC-size histogram", fmt(tournament_sccs))
    print("tournament Hamiltonian-path-count histogram", fmt(tournament_paths))
    print("tournament tie-edge histogram", fmt(tournament_ties))
    print("tournament edge-flip histogram", fmt(tournament_flips))
    print("challenged vertices: mod-3 and mod-9 sheet classes preserve the ramified incidence signatures; owner tournaments preserve the terminal deficit vector but destroy sheet colours, mask intersections, and the absolute 27-sheet threshold")
    print("frontier verdict scalar-empty no; owner-local all-six empty yes; scale 27 common-sheet H6 face closed")


if __name__ == "__main__":
    main()
