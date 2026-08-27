#!/usr/bin/env python3
"""Exact finite controls for THM-4275's confluent observer matroids.

The recurrence enumeration concerns abstract F4 sidecars. It does not assert
that every word is geometrically realized by a Keller-map incidence.
"""

from collections import Counter, defaultdict
from itertools import combinations, product
import sys


def require(condition, message="audit failure"):
    if not condition:
        raise RuntimeError(message)


# F4=F2[w]/(w^2+w+1), encoded in the basis (1,w) by two bits.
def add(left, right):
    return left ^ right


def multiply(left, right):
    a0, a1 = left & 1, (left >> 1) & 1
    b0, b1 = right & 1, (right >> 1) & 1
    c0 = (a0 * b0) ^ (a1 * b1)
    c1 = (a0 * b1) ^ (a1 * b0) ^ (a1 * b1)
    return c0 | (c1 << 1)


def f4_power(value, exponent):
    answer = 1
    for _ in range(exponent):
        answer = multiply(answer, value)
    return answer


OMEGA = 2
ALPHA = multiply(OMEGA, OMEGA)
require(ALPHA == 3 and multiply(OMEGA, ALPHA) == 1, "bad F4 encoding")


def ambient_word(c0, c1, c2):
    """THM-4258 cubic word alpha^j(c0+j*c1+C(j,2)*c2)."""
    word = []
    for index in range(12):
        value = c0
        if index & 1:
            value = add(value, c1)
        if (index * (index - 1) // 2) & 1:
            value = add(value, c2)
        word.append(multiply(f4_power(ALPHA, index), value))
    return tuple(word)


def visible_word(even_seed, odd_seed):
    """THM-4264 visible word, with one independent seed per parity."""
    word = []
    for index in range(12):
        parity_step = index // 2
        seed = even_seed if index % 2 == 0 else odd_seed
        word.append(multiply(f4_power(OMEGA, parity_step), seed))
    return tuple(word)


AMBIENT = tuple(ambient_word(*seed) for seed in product(range(4), repeat=3))
VISIBLE = tuple(visible_word(*seed) for seed in product(range(4), repeat=2))
require(len(set(AMBIENT)) == 64, "ambient parametrization collision")
require(len(set(VISIBLE)) == 16, "visible parametrization collision")
require(set(VISIBLE).issubset(set(AMBIENT)), "visible module not ambient")


def observer_profile(words, sample):
    fibres = defaultdict(int)
    for word in words:
        fibres[tuple(word[index] for index in sample)] += 1
    fibre_sizes = tuple(sorted(fibres.values()))
    require(len(set(fibre_sizes)) == 1, ("irregular linear fibres", sample))
    return len(fibres), fibre_sizes[0]


def rank_from_image_size(image_size):
    ranks = {1: 0, 4: 1, 16: 2, 64: 3}
    require(image_size in ranks, "non-F4-linear image size")
    return ranks[image_size]


ambient_subset_rank_counts = Counter()
visible_subset_rank_counts = Counter()
ambient_regular_image_fibre_cells = 0
visible_regular_image_fibre_cells = 0

for mask in range(1 << 12):
    sample = tuple(index for index in range(12) if mask & (1 << index))

    image, fibre = observer_profile(AMBIENT, sample)
    rank = rank_from_image_size(image)
    expected = min(len({index % 4 for index in sample}), 3)
    require(rank == expected, ("ambient sample rank", sample, rank, expected))
    require(image * fibre == 64, "ambient fibre ledger")
    ambient_subset_rank_counts[rank] += 1
    ambient_regular_image_fibre_cells += image

    image, fibre = observer_profile(VISIBLE, sample)
    rank = rank_from_image_size(image)
    expected = len({index % 2 for index in sample})
    require(rank == expected, ("visible sample rank", sample, rank, expected))
    require(image * fibre == 16, "visible fibre ledger")
    visible_subset_rank_counts[rank] += 1
    visible_regular_image_fibre_cells += image


def fixed_size_rank_counts(words, size):
    counts = Counter()
    for sample in combinations(range(12), size):
        image, _ = observer_profile(words, sample)
        counts[rank_from_image_size(image)] += 1
    return counts


ambient_pair_ranks = fixed_size_rank_counts(AMBIENT, 2)
ambient_triple_ranks = fixed_size_rank_counts(AMBIENT, 3)
visible_pair_ranks = fixed_size_rank_counts(VISIBLE, 2)
visible_triple_ranks = fixed_size_rank_counts(VISIBLE, 3)
require(ambient_pair_ranks == Counter({2: 54, 1: 12}), "ambient pair census")
require(
    ambient_triple_ranks == Counter({3: 108, 2: 108, 1: 4}),
    "ambient triple census",
)
require(visible_pair_ranks == Counter({2: 36, 1: 30}), "visible pair census")
require(
    visible_triple_ranks == Counter({2: 180, 1: 40}),
    "visible triple census",
)

# A one-edge quotient has perfectly regular density but cannot recognize the
# singleton zero/full-collapse fibre.
one_edge_zero_hostiles = 0
for index in range(12):
    image, fibre = observer_profile(VISIBLE, (index,))
    require((image, fibre) == (4, 4), "one-edge fibre profile")
    zero_fibre = [word for word in VISIBLE if word[index] == 0]
    nonzero_hostiles = sum(any(value != 0 for value in word) for word in zero_fibre)
    require(len(zero_fibre) == 4 and nonzero_hostiles == 3, "one-edge hostile")
    one_edge_zero_hostiles += nonzero_hostiles

# Every opposite-parity pair is faithful; every same-parity pair has a
# four-word kernel. This includes cyclic gaps 3 and 9, which do not generate
# Z/12 additively.
opposite_pair_checks = 0
same_pair_hostiles = 0
gap_counts = Counter()
for left, right in combinations(range(12), 2):
    image, fibre = observer_profile(VISIBLE, (left, right))
    if (left - right) & 1:
        require((image, fibre) == (16, 1), "opposite-parity collision")
        opposite_pair_checks += 1
        gap = (right - left) % 12
        gap_counts[min(gap, 12 - gap)] += 1
    else:
        require((image, fibre) == (4, 4), "same-parity profile")
        same_pair_hostiles += 1
require(opposite_pair_checks == 36, "opposite-pair count")
require(same_pair_hostiles == 30, "same-pair count")
require(gap_counts == Counter({1: 12, 3: 12, 5: 12}), "odd-gap census")


# Fat-contact model b^12=Lambda. Over F13, 2 has order 12. The generic
# Kummer-evaluation matrix is Fourier--Vandermonde rank 12; at a=0 it has
# rank 1, whereas the coefficient/Hasse observer remains the identity.
def matrix_rank_mod(matrix, prime):
    rows = [list(row) for row in matrix]
    rank = 0
    columns = len(rows[0]) if rows else 0
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, len(rows)) if rows[row][column] % prime),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column] % prime, -1, prime)
        rows[rank] = [(entry * inverse) % prime for entry in rows[rank]]
        for row in range(len(rows)):
            if row == rank or rows[row][column] % prime == 0:
                continue
            scalar = rows[row][column] % prime
            rows[row] = [
                (entry - scalar * pivot_entry) % prime
                for entry, pivot_entry in zip(rows[row], rows[rank])
            ]
        rank += 1
    return rank


PRIME = 13
ZETA = 2
require(
    pow(ZETA, 12, PRIME) == 1
    and all(pow(ZETA, exponent, PRIME) != 1 for exponent in (1, 2, 3, 4, 6)),
    "bad twelfth root",
)
contact_ranks = Counter()
contact_fourier_basis_checks = 0
for parameter in range(PRIME):
    matrix = [
        [
            pow((pow(ZETA, row, PRIME) * parameter) % PRIME, column, PRIME)
            for column in range(12)
        ]
        for row in range(12)
    ]
    rank = matrix_rank_mod(matrix, PRIME)
    contact_ranks[rank] += 1
    require(rank == (1 if parameter == 0 else 12), "contact rank")
    if parameter == 0:
        continue
    for basis_column in range(12):
        values = [matrix[row][basis_column] for row in range(12)]
        for fourier_column in range(12):
            recovered = (
                pow(12, -1, PRIME)
                * sum(
                    pow(ZETA, (-row * fourier_column) % 12, PRIME) * values[row]
                    for row in range(12)
                )
            ) % PRIME
            expected = (
                pow(parameter, fourier_column, PRIME)
                if fourier_column == basis_column
                else 0
            )
            require(recovered == expected, "Fourier recovery")
            contact_fourier_basis_checks += 1

coalesced_b_values = (0,) * 12
hasse_b_coefficients = tuple(1 if degree == 1 else 0 for degree in range(12))
require(set(coalesced_b_values) == {0}, "coalesced value hostile")
require(hasse_b_coefficients[1] == 1, "normal-Hasse hostile")


# Exact finite-shell transverse truncation criterion. Duplicate entries test
# full fibres, and unequal entries realize several pair-contact orders.
series_shell = (
    (0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0),
    (0, 1, 0, 0, 0, 0),
    (0, 1, 1, 0, 0, 0),
    (0, 1, 1, 1, 0, 0),
    (0, 1, 1, 1, 1, 0),
    (0, 1, 1, 1, 1, 1),
    (1, 0, 0, 0, 0, 0),
)


def contact_order(left, right):
    for index, (a, b) in enumerate(zip(left, right)):
        if a != b:
            return index
    return None


orders = [
    contact_order(left, right)
    for left, right in combinations(series_shell, 2)
    if left != right
]
minimum_complete_depth = 1 + max(orders)
truncation_pair_checks = 0
for depth in range(0, len(series_shell[0]) + 2):
    same_fibres = True
    for left in series_shell:
        for right in series_shell:
            full_equal = left == right
            truncated_equal = left[:depth] == right[:depth]
            same_fibres &= full_equal == truncated_equal
            truncation_pair_checks += 1
    require(same_fibres == (depth >= minimum_complete_depth), "jet depth gate")


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    print("jc23_confluent_sample_matroids_thm4275")
    print("scope=ABSTRACT_F4_SIDECARS_NOT_GEOMETRIC_REALIZATION")
    print(f"ambient_words={len(AMBIENT)}")
    print(f"visible_words={len(VISIBLE)}")
    print("ambient_sample_rule=rank=min(number_of_residues_mod_4,3)")
    print("visible_sample_rule=rank=number_of_parities_seen")
    print(f"ambient_subset_rank_counts={dict(sorted(ambient_subset_rank_counts.items()))}")
    print(f"visible_subset_rank_counts={dict(sorted(visible_subset_rank_counts.items()))}")
    print(f"ambient_pair_rank_counts={dict(sorted(ambient_pair_ranks.items()))}")
    print(f"ambient_triple_rank_counts={dict(sorted(ambient_triple_ranks.items()))}")
    print(f"visible_pair_rank_counts={dict(sorted(visible_pair_ranks.items()))}")
    print(f"visible_triple_rank_counts={dict(sorted(visible_triple_ranks.items()))}")
    print("ambient_all_subset_checks=4096")
    print("visible_all_subset_checks=4096")
    print(f"ambient_regular_image_fibre_cells={ambient_regular_image_fibre_cells}")
    print(f"visible_regular_image_fibre_cells={visible_regular_image_fibre_cells}")
    print("visible_one_edge=image_4_fibre_4_zero_fibre_has_3_nonzero_words")
    print(f"visible_one_edge_nonzero_zero_fibre_hostiles={one_edge_zero_hostiles}")
    print(f"visible_opposite_parity_faithful_pairs={opposite_pair_checks}")
    print(f"visible_opposite_parity_gap_counts={dict(sorted(gap_counts.items()))}")
    print(f"visible_same_parity_unfaithful_pairs={same_pair_hostiles}")
    print(f"contact_evaluation_rank_counts_F13={dict(sorted(contact_ranks.items()))}")
    print(f"contact_fourier_basis_checks={contact_fourier_basis_checks}")
    print("contact_a0_value_rank=1")
    print("contact_a0_full_hasse_rank=12")
    print("contact_hostile=b_has_zero_coalesced_values_but_first_hasse_coefficient_1")
    print(f"transverse_series_shell={len(series_shell)}")
    print(f"transverse_minimum_complete_depth={minimum_complete_depth}")
    print(f"transverse_truncation_pair_checks={truncation_pair_checks}")
    print("status=PASS")


if __name__ == "__main__":
    main()
