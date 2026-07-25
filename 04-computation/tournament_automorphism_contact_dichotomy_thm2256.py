#!/usr/bin/env python3
"""Exact finite audit for THM-2256.

The proof of the bounded/linear dichotomy is uniform.  This companion checks
the interaction identities on every labelled tournament through order five
and computes the exact minimum over all Hall-layer multisets for orders at
most four and margins at most four.
"""

from collections import Counter
from itertools import combinations_with_replacement, permutations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def tournament_from_mask(order: int, mask: int) -> tuple[tuple[bool, ...], ...]:
    adjacency = [[False for _ in range(order)] for _ in range(order)]
    bit = 0
    for left in range(order):
        for right in range(left + 1, order):
            if (mask >> bit) & 1:
                adjacency[left][right] = True
            else:
                adjacency[right][left] = True
            bit += 1
    return tuple(tuple(row) for row in adjacency)


def is_automorphism(adjacency, sigma) -> bool:
    order = len(adjacency)
    return all(
        adjacency[left][right]
        == adjacency[sigma[left]][sigma[right]]
        for left in range(order)
        for right in range(order)
    )


def interaction(adjacency, sigma, tau) -> int:
    order = len(adjacency)
    return sum(
        adjacency[left][right]
        and adjacency[tau[right]][sigma[left]]
        for left in range(order)
        for right in range(order)
    )


def interaction_data(adjacency, full_matrix: bool = True):
    order = len(adjacency)
    sym = tuple(permutations(range(order)))
    diagonal = tuple(
        interaction(adjacency, sigma, sigma) for sigma in sym
    )
    automorphism_indices = tuple(
        index
        for index, sigma in enumerate(sym)
        if is_automorphism(adjacency, sigma)
    )
    require(
        tuple(index for index, value in enumerate(diagonal) if value == 0)
        == automorphism_indices,
        "diagonal-zero/automorphism mismatch",
    )
    for position, left in enumerate(automorphism_indices):
        for right in automorphism_indices[position + 1 :]:
            require(
                interaction(adjacency, sym[left], sym[right])
                + interaction(adjacency, sym[right], sym[left])
                > 0,
                "distinct automorphisms have zero contact",
            )
    nonautomorphism_indices = tuple(
        index
        for index in range(len(sym))
        if index not in set(automorphism_indices)
    )
    contacts = tuple(
        (auto, nonauto)
        for auto in automorphism_indices
        for nonauto in nonautomorphism_indices
        if interaction(adjacency, sym[auto], sym[nonauto])
        + interaction(adjacency, sym[nonauto], sym[auto])
        == 0
    )
    weights = None
    if full_matrix:
        weights = tuple(
            tuple(
                interaction(adjacency, sigma, tau)
                + interaction(adjacency, tau, sigma)
                for tau in sym
            )
            for sigma in sym
        )
    return (
        sym,
        diagonal,
        weights,
        automorphism_indices,
        nonautomorphism_indices,
        contacts,
    )


def energy(counts, diagonal, weights) -> int:
    active = tuple(sorted(counts))
    value = sum(
        diagonal[index] * counts[index] * counts[index]
        for index in active
    )
    value += sum(
        weights[left][right] * counts[left] * counts[right]
        for position, left in enumerate(active)
        for right in active[position + 1 :]
    )
    return value


def exact_layer_minimum(adjacency, margin: int) -> tuple[int, tuple[int, ...]]:
    (
        sym,
        diagonal,
        weights,
        automorphisms,
        _,
        contacts,
    ) = interaction_data(adjacency)
    automorphism_set = set(automorphisms)
    best_value = None
    best_layers = None
    for layers in combinations_with_replacement(range(len(sym)), margin):
        if len(set(layers)) == 1 and layers[0] in automorphism_set:
            continue
        counts = Counter(layers)
        value = energy(counts, diagonal, weights)
        if best_value is None or value < best_value:
            best_value = value
            best_layers = layers
    require(best_value is not None, "empty nonfree layer universe")
    require(best_value >= 1, "positive integer floor failed")
    if contacts:
        contact_ceiling = min(diagonal[nonauto] for _, nonauto in contacts)
        require(best_value <= contact_ceiling, "bounded contact witness failed")
    else:
        require(
            best_value >= margin - 1,
            "linear no-contact lower bound failed",
        )
    return best_value, best_layers


def canonical_transitive(order: int):
    return tuple(
        tuple(left < right for right in range(order))
        for left in range(order)
    )


def main() -> None:
    census = {}
    diagonal_checks = 0
    automorphism_pair_checks = 0
    contact_witnesses = 0
    no_contact_witnesses = 0

    for order in range(2, 6):
        contact_count = 0
        no_contact_count = 0
        for mask in range(1 << (order * (order - 1) // 2)):
            adjacency = tournament_from_mask(order, mask)
            (
                sym,
                diagonal,
                weights,
                automorphisms,
                _,
                contacts,
            ) = interaction_data(adjacency, full_matrix=False)
            diagonal_checks += len(sym)
            automorphism_pair_checks += (
                len(automorphisms) * (len(automorphisms) - 1) // 2
            )
            if contacts:
                contact_count += 1
                contact_witnesses += 1
                auto, nonauto = contacts[0]
                require(
                    diagonal[nonauto] >= 1
                    and interaction(
                        adjacency, sym[auto], sym[nonauto]
                    )
                    + interaction(
                        adjacency, sym[nonauto], sym[auto]
                    )
                    == 0,
                    "invalid bounded-class witness",
                )
            else:
                no_contact_count += 1
                no_contact_witnesses += 1
                require(
                    all(
                        interaction(
                            adjacency, sym[auto], sym[nonauto]
                        )
                        + interaction(
                            adjacency, sym[nonauto], sym[auto]
                        )
                        >= 1
                        for auto in automorphisms
                        for nonauto in range(len(sym))
                        if nonauto not in set(automorphisms)
                    ),
                    "invalid linear-class witness",
                )
        census[order] = (contact_count, no_contact_count)

    exact_minima = {}
    layer_multisets = 0
    for order in range(2, 5):
        rows = []
        for mask in range(1 << (order * (order - 1) // 2)):
            adjacency = tournament_from_mask(order, mask)
            margin_values = []
            for margin in range(1, 5):
                value, _ = exact_layer_minimum(adjacency, margin)
                margin_values.append(value)
                permutation_count = 1
                for numerator in range(
                    len(tuple(permutations(range(order)))),
                    len(tuple(permutations(range(order)))) + margin,
                ):
                    permutation_count *= numerator
                for denominator in range(1, margin + 1):
                    permutation_count //= denominator
                layer_multisets += permutation_count
            rows.append(tuple(margin_values))
        exact_minima[order] = (
            min(value for row in rows for value in row),
            max(value for row in rows for value in row),
        )

    transitive_controls = {}
    for order in range(2, 6):
        adjacency = canonical_transitive(order)
        (
            sym,
            diagonal,
            weights,
            automorphisms,
            _,
            contacts,
        ) = interaction_data(adjacency)
        require(len(automorphisms) == 1, "transitive automorphism drift")
        identity = automorphisms[0]
        adjacent_swap = list(range(order))
        adjacent_swap[-2], adjacent_swap[-1] = (
            adjacent_swap[-1],
            adjacent_swap[-2],
        )
        swap_index = sym.index(tuple(adjacent_swap))
        require(
            diagonal[swap_index] == 1
            and weights[identity][swap_index] == 0
            and (identity, swap_index) in contacts,
            "transitive unit-contact witness failed",
        )
        transitive_controls[order] = (
            diagonal[swap_index],
            weights[identity][swap_index],
        )

    require(
        sum(sum(values) for values in census.values()) == 1098,
        "labelled tournament census drift",
    )

    print("THM-2256 AUTOMORPHISM-CONTACT AUDIT")
    print(f"labelled_tournament_census={census}")
    print(f"diagonal_checks={diagonal_checks}")
    print(f"distinct_automorphism_pair_checks={automorphism_pair_checks}")
    print(
        "class_witness_counts="
        f"bounded:{contact_witnesses},linear:{no_contact_witnesses}"
    )
    print(f"exact_layer_multisets_q2_to_q4_N1_to_N4={layer_multisets}")
    print(f"exact_minimum_range_by_order={exact_minima}")
    print(f"transitive_unit_controls={transitive_controls}")
    print("status=PASS")


if __name__ == "__main__":
    main()
