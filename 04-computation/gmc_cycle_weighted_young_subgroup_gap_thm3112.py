#!/usr/bin/env python3
"""Exact companion for THM-3112.

The universal Young-subgroup decomposition is proved symbolically in the
theorem.  This companion checks its finite group-algebra normalization on
small labelled alphabets and derives the sharp degree-six THM-3110 boundary
from the live response banks using rational arithmetic only.
"""

import importlib.util
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations, product
from math import comb, factorial, prod
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
SPEC = importlib.util.spec_from_file_location("thm3110_schur", UPSTREAM)
THM3110 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(THM3110)


def integer_partitions(total, largest=None):
    if total == 0:
        return ((),)
    if largest is None or largest > total:
        largest = total
    answer = []
    for first in range(largest, 0, -1):
        for tail in integer_partitions(total - first, first):
            answer.append((first,) + tail)
    return tuple(answer)


def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def power_product(alphabet, partition):
    return prod(sum(value**power for value in alphabet)
                for power in partition)


def complete_homogeneous(alphabet, degree):
    complete = [0] * (degree + 1)
    complete[0] = 1
    for value in alphabet:
        for index in range(1, degree + 1):
            complete[index] += value * complete[index - 1]
    return complete[degree]


def monomial_symmetric(alphabet, partition):
    require(len(partition) <= len(alphabet), "monomial shape is too long")
    exponents = partition + (0,) * (len(alphabet) - len(partition))
    return sum(
        prod(value**power for value, power in zip(alphabet, powers))
        for powers in set(permutations(exponents))
    )


def colouring_identity(total, alphabet):
    """Check Omega=sum_f a_f P_H coefficient by coefficient."""

    group = tuple(permutations(range(total)))
    denominator = factorial(total)
    right = defaultdict(Fraction)
    total_weight = Fraction(0)
    for colouring in product(range(len(alphabet)), repeat=total):
        weight = prod(alphabet[colour] for colour in colouring)
        counts = Counter(colouring)
        subgroup_size = prod(factorial(count) for count in counts.values())
        coefficient = Fraction(weight, denominator)
        total_weight += coefficient * subgroup_size
        for permutation in group:
            if all(colouring[index] == colouring[permutation[index]]
                   for index in range(total)):
                right[permutation] += coefficient

    left = {
        permutation: Fraction(
            power_product(alphabet, cycle_type(permutation)), denominator)
        for permutation in group
    }
    require(dict(right) == left, "cycle-colouring group-algebra identity failed")
    require(total_weight == complete_homogeneous(alphabet, total),
            "colouring weights do not sum to h_N")
    return len(group), len(alphabet) ** total


def z_value(partition):
    answer = 1
    for part, multiplicity in Counter(partition).items():
        answer *= part**multiplicity * factorial(multiplicity)
    return answer


def hook_dimension(partition):
    denominator = 1
    for row, length in enumerate(partition):
        for column in range(length):
            denominator *= (
                length - column
                + sum(other > column for other in partition[row + 1:])
            )
    return factorial(sum(partition)) // denominator


def matrix_rank(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    if not work:
        return 0
    row = 0
    for column in range(len(work[0])):
        pivot = next((index for index in range(row, len(work))
                      if work[index][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        scale = work[row][column]
        work[row] = [value / scale for value in work[row]]
        for index in range(len(work)):
            if index == row or not work[index][column]:
                continue
            scale = work[index][column]
            work[index] = [left - scale * right
                           for left, right in zip(work[index], work[row])]
        row += 1
        if row == len(work):
            break
    return row


def q_and_spectra(invariant, first=1, second=2, total=6):
    bank = THM3110.BANKS[invariant]
    dominant = (
        ((0, 3), (2, 0), (2, 0), (1, 0))
        if invariant == 0 else
        ((0, 3), (1, 1), (2, 0), (1, 0), (1, 0))
    )
    dominant_coefficient = (1, 2)[invariant]
    base = THM3110.forced_divisor(invariant, first, second)
    q_roots = THM3110.residual_roots(
        invariant, dominant, first, second)
    q_data = THM3110.complete_and_elementary(q_roots, total)
    q_h = THM3110.schur_value(q_data, (total,))

    atom_data = []
    for coefficient, row in bank:
        roots = THM3110.residual_roots(invariant, row, first, second)
        atom_data.append((coefficient, roots,
                          THM3110.complete_and_elementary(roots, total)))
    phi_h = sum(coefficient * THM3110.schur_value(data, (total,))
                for coefficient, _, data in atom_data)
    rate = Fraction(phi_h, dominant_coefficient * q_h)

    class_values = {}
    for partition in integer_partitions(total):
        phi_power = sum(
            coefficient * power_product(roots, partition)
            for coefficient, roots, _ in atom_data
        )
        q_power = power_product(q_roots, partition)
        class_values[partition] = Fraction(
            phi_power - rate * dominant_coefficient * q_power, base)

    require(sum(class_values[partition] / z_value(partition)
                for partition in class_values) == 0,
            "row-normalized central element has nonzero augmentation")

    spectra = {}
    for partition in integer_partitions(total):
        phi_schur = sum(
            coefficient * THM3110.schur_value(data, partition)
            for coefficient, _, data in atom_data
        )
        q_schur = THM3110.schur_value(q_data, partition)
        spectra[partition] = Fraction(
            phi_schur - rate * dominant_coefficient * q_schur,
            base * hook_dimension(partition),
        )
    require(spectra[(total,)] == 0, "trivial Specht scalar is nonzero")
    require(all(value >= 0 for value in spectra.values()),
            "degree-six row-normalized operator is not PSD")
    require(all((value == 0) == (partition == (total,))
                for partition, value in spectra.items()),
            "unexpected degree-six equality boundary")
    return rate, class_values, spectra


def signed_beta_coefficient_hostile(invariant, first=1, second=2):
    total = 5
    shape = (2, 1, 1, 1)
    bank = THM3110.BANKS[invariant]
    dominant = (
        ((0, 3), (2, 0), (2, 0), (1, 0))
        if invariant == 0 else
        ((0, 3), (1, 1), (2, 0), (1, 0), (1, 0))
    )
    dominant_coefficient = (1, 2)[invariant]
    q_roots = THM3110.residual_roots(
        invariant, dominant, first, second)
    atom_roots = tuple(
        (coefficient,
         THM3110.residual_roots(invariant, row, first, second))
        for coefficient, row in bank
    )
    phi_h = sum(coefficient * complete_homogeneous(roots, total)
                for coefficient, roots in atom_roots)
    q_h = complete_homogeneous(q_roots, total)
    row_rate = Fraction(phi_h, dominant_coefficient * q_h)
    phi_m = sum(coefficient * monomial_symmetric(roots, shape)
                for coefficient, roots in atom_roots)
    q_m = monomial_symmetric(q_roots, shape)
    monomial_rate = Fraction(phi_m, dominant_coefficient * q_m)
    return row_rate, monomial_rate, row_rate - monomial_rate


def central_alpha_matrix(total, classes):
    """Nonidentity coefficients of sum_|A|=k alpha_A, k=2..N."""

    columns = []
    for size in range(2, total + 1):
        column = []
        for partition in classes:
            support = sum(part for part in partition if part > 1)
            coefficient = (
                -Fraction(comb(total - support, size - support),
                          factorial(size - 1))
                if size >= support else Fraction(0)
            )
            column.append(coefficient)
        columns.append(column)
    return [[columns[column][row] for column in range(len(columns))]
            for row in range(len(classes))]


colour_controls = []
for total, alphabet in ((3, (0, 1, 2)), (4, (1, 2)), (5, (0, 1, 3))):
    group_size, colouring_count = colouring_identity(total, alphabet)
    colour_controls.append(
        f"N{total}:alphabet={alphabet}:perms={group_size}:"
        f"colourings={colouring_count}"
    )

records = []
classes = tuple(partition for partition in integer_partitions(6)
                if partition != (1, 1, 1, 1, 1, 1))
alpha_matrix = central_alpha_matrix(6, classes)
alpha_rank = matrix_rank(alpha_matrix)
require(alpha_rank == 5, "central alpha span rank drift")

for invariant in (0, 1):
    rate, class_values, spectra = q_and_spectra(invariant)
    positive_nonidentity = tuple(
        (partition, value)
        for partition, value in class_values.items()
        if partition != (1, 1, 1, 1, 1, 1) and value > 0
    )
    require(tuple(partition for partition, _ in positive_nonidentity)
            == ((2, 1, 1, 1, 1),),
            "degree-six positive class boundary drift")
    support_four = (
        class_values[(4, 1, 1)], class_values[(2, 2, 1, 1)])
    require(support_four[0] != support_four[1],
            "same-support class obstruction disappeared")
    augmented_rank = matrix_rank([
        row + [class_values[partition]]
        for row, partition in zip(alpha_matrix, classes)
    ])
    require(augmented_rank == 6,
            "D6 unexpectedly entered the central alpha span")
    records.append((invariant + 1, rate, positive_nonidentity,
                    support_four, min(value for partition, value in spectra.items()
                                      if partition != (6,))))

beta_hostiles = tuple(
    signed_beta_coefficient_hostile(invariant) for invariant in (0, 1)
)
require(beta_hostiles[0] == (
    Fraction(6, 8447), Fraction(10, 3007),
    Fraction(-66428, 25400129)), "I1 signed-beta hostile drift")
require(beta_hostiles[1][2] == Fraction(-238351, 278093196),
        "I2 signed-beta hostile drift")

print("cycle_colouring=" + ";".join(colour_controls))
print("universal_identity=Omega_N(S)=sum_f a_f*P_Hf;sum_f_a_f=h_N(S)")
print("gap_decomposition=beta_N(S)=sum_f a_f*(identity-P_Hf):direct_PSD")
print("singleton_boundary=N*beta_N((1))=alpha_[N]")
print(f"central_alpha_N6=classes={len(classes)}:generators=5:rank={alpha_rank}")
for invariant, rate, positive, support_four, minimum in records:
    partition, value = positive[0]
    print(
        f"I{invariant}_N6=rate={rate}:positive_nonidentity={partition}:{value}:"
        f"q411={support_four[0]}:q2211={support_four[1]}:"
        f"augmented_rank=6:min_nontrivial_spectrum={minimum}"
    )
print("uniform_octopus_boundary=D6_not_in_conjugacy_averaged_span_of_alpha_A")
print(
    "signed_beta_hostile_I1=mu(2,1,1,1):"
    f"row={beta_hostiles[0][0]}:monomial={beta_hostiles[0][1]}:"
    f"difference={beta_hostiles[0][2]}"
)
print(
    "signed_beta_hostile_I2=mu(2,1,1,1):"
    f"row={beta_hostiles[1][0]}:monomial={beta_hostiles[1][1]}:"
    f"difference={beta_hostiles[1][2]}"
)
print("open_target=bank_level_stochastic_refinement_or_rooted_star_mesh_inequality")
print("all_exact_checks=PASS")
