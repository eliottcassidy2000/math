#!/usr/bin/env python3
"""Exact companion for THM-3128's full-kernel no-preimage boundary.

The script reconstructs the first active pole-prefix Hasse hostile directly
from the THM-3110 bank, then keeps the two relevant gauges separate:

* raw coefficients G_lambda multiply the Young gaps Lbar_lambda;
* factorial coordinates c_lambda=G_lambda/w_lambda multiply
  Kbar_lambda=w_lambda Lbar_lambda and are acted on by labelled deletion A.

It proves that the full affine fibre with the same labelled-deletion image
has a fixed negative coarsening-upset mass, so no point of that fibre is a
nonnegative Hasse boundary.  Only integers and exact fractions are used.
"""

from collections import Counter
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import chain
from math import factorial
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent


def load_companion(name, filename):
    spec = spec_from_file_location(name, HERE / filename)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


g = load_companion(
    "thm3110_companion",
    "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py",
)
f = load_companion(
    "thm3120_companion",
    "gmc_product_gamma_row_pole_flag_thm3120.py",
)


SOURCE_TYPES = (
    (5,), (4, 1), (3, 2), (3, 1, 1),
    (2, 2, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1),
)
TARGET_TYPES = ((4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1))
SOURCE_INDEX = {shape: index for index, shape in enumerate(SOURCE_TYPES)}
TARGET_INDEX = {shape: index for index, shape in enumerate(TARGET_TYPES)}


def weight(shape):
    answer = 1
    for part in shape:
        answer *= factorial(part)
    return answer


def lowered(shape, part):
    answer = list(shape)
    answer.remove(part)
    if part > 1:
        answer.append(part - 1)
    return tuple(sorted(answer, reverse=True))


def labelled_deletion_matrix():
    matrix = [[0] * len(SOURCE_TYPES) for _ in TARGET_TYPES]
    for column, shape in enumerate(SOURCE_TYPES):
        for part, multiplicity in Counter(shape).items():
            matrix[TARGET_INDEX[lowered(shape, part)]][column] += (
                part * multiplicity
            )
    return tuple(tuple(row) for row in matrix)


def matvec(matrix, vector):
    return tuple(sum(entry * value for entry, value in zip(row, vector))
                 for row in matrix)


def rational_rank(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rank = 0
    for column in range(len(work[0]) if work else 0):
        pivot = next((row for row in range(rank, len(work))
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][column]
        work[rank] = [value / scale for value in work[rank]]
        for row in range(len(work)):
            if row == rank or work[row][column] == 0:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right
                         for left, right in zip(work[row], work[rank])]
        rank += 1
    return rank


def determinant(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    answer = Fraction(1)
    for column in range(len(work)):
        pivot = next((row for row in range(column, len(work))
                      if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        value = work[column][column]
        answer *= value
        for index in range(column, len(work)):
            work[column][index] /= value
        for row in range(column + 1, len(work)):
            scale = work[row][column]
            for index in range(column, len(work)):
                work[row][index] -= scale * work[column][index]
    require(answer.denominator == 1, "integral determinant lost integrality")
    return answer.numerator


def all_monomial_values(roots, removed, maximum_degree=5):
    """Evaluate every m_lambda on the virtual alphabet roots-removed."""

    power_sums = tuple(
        sum(root ** degree for root in roots)
        - sum(root ** degree for root in removed)
        for degree in range(1, maximum_degree + 1)
    )
    shapes = tuple(sorted(
        chain.from_iterable(g.partitions(degree)
                            for degree in range(1, maximum_degree + 1)),
        key=lambda shape: (sum(shape), len(shape), shape),
    ))
    values = {(): 1}
    for shape in shapes:
        exponent = shape[-1]
        remainder = shape[:-1]
        value = power_sums[exponent - 1] * values[remainder]
        for old_exponent in set(remainder):
            merged = list(remainder)
            merged.remove(old_exponent)
            merged.append(old_exponent + exponent)
            merged = tuple(sorted(merged, reverse=True))
            value -= Counter(merged)[old_exponent + exponent] * values[merged]
        multiplicity = Counter(shape)[exponent]
        require(value % multiplicity == 0,
                "virtual monomial recurrence lost integrality")
        values[shape] = value // multiplicity
    return values


def active_prefix_current():
    invariant, a, b, prefix_length, degree = 1, 1, 3, 5, 5
    numerator, _, remaining, _, _, _ = f.reduced_row_fraction(
        invariant, a, b)
    polynomial = tuple(numerator[5:])
    poles = tuple(sorted(remaining.elements(), reverse=True))
    removed = poles[:prefix_length]
    require(len(polynomial) - 1 == prefix_length,
            "selected prefix is no longer active-top")

    atom_values = tuple(
        (coefficient, all_monomial_values(
            g.residual_roots(invariant, rows, a, b), removed))
        for coefficient, rows in g.BANKS[invariant]
    )
    dominant_rows = tuple(sorted(
        ((0, 3), (1, 1), (2, 0), (1, 0), (1, 0)),
        reverse=True,
    ))
    quotient_roots = g.residual_roots(
        invariant, dominant_rows, a, b)
    quotient = all_monomial_values(quotient_roots, removed)

    phi = tuple(
        sum(coefficient * values[shape]
            for coefficient, values in atom_values)
        for shape in SOURCE_TYPES
    )
    phi_h = sum(phi)
    quotient_vector = tuple(quotient[shape] for shape in SOURCE_TYPES)
    quotient_h = sum(quotient_vector)
    current = tuple(
        phi_h * q_value - phi_value * quotient_h
        for phi_value, q_value in zip(phi, quotient_vector)
    )
    require(sum(current) == 0, "active current lost zero raw mass")
    require(phi[SOURCE_INDEX[(5,)]] == 0,
            "virtual-prefixed power-sum response did not vanish")
    return (
        invariant, a, b, degree, poles, removed, polynomial,
        quotient_roots, phi_h, quotient_h,
        phi[SOURCE_INDEX[(5,)]], quotient[SOURCE_TYPES[0]], current,
    )


def main():
    (
        invariant, a, b, degree, poles, removed, polynomial,
        quotient_roots, phi_h, quotient_h, phi_power, quotient_power, raw,
    ) = active_prefix_current()
    expected_raw = (
        -2324160, 544320, 2237760, -656640,
        -915840, 972000, 142560,
    )
    require(raw == expected_raw, "active hostile coefficient vector drift")
    require((phi_h, quotient_h, phi_power, quotient_power)
            == (1440, -358, 0, -1614),
            "coarsest component data drift")

    weights = tuple(weight(shape) for shape in SOURCE_TYPES)
    target_weights = tuple(weight(shape) for shape in TARGET_TYPES)
    require(weights == (120, 24, 12, 6, 4, 2, 1),
            "degree-five factorial weights drift")
    factorial_coordinates = tuple(Fraction(value, scalar)
                                  for value, scalar in zip(raw, weights))
    require(all(value.denominator == 1 for value in factorial_coordinates),
            "hostile factorial coordinates became nonintegral")
    factorial_coordinates = tuple(value.numerator
                                  for value in factorial_coordinates)
    require(factorial_coordinates == (
        -19368, 22680, 186480, -109440,
        -228960, 486000, 142560,
    ), "factorial-coordinate hostile drift")

    labelled = labelled_deletion_matrix()
    expected_labelled = (
        (5, 1, 0, 0, 0, 0, 0),
        (0, 4, 2, 2, 0, 0, 0),
        (0, 0, 3, 0, 1, 0, 0),
        (0, 0, 0, 3, 4, 3, 0),
        (0, 0, 0, 0, 0, 2, 5),
    )
    require(labelled == expected_labelled,
            "labelled-deletion matrix drift")
    minor = tuple(tuple(row[column] for column in (0, 1, 2, 3, 5))
                  for row in labelled)
    require(determinant(minor) == 360 and rational_rank(labelled) == 5,
            "labelled-deletion rank certificate failed")

    kernel_1 = (-1, 5, -2, -8, 6, 0, 0)
    kernel_2 = (1, -5, 0, 10, 0, -10, 4)
    require(matvec(labelled, kernel_1) == (0,) * 5
            and matvec(labelled, kernel_2) == (0,) * 5,
            "declared labelled-deletion kernel vector escaped")
    require(kernel_1[0] * kernel_2[2] - kernel_1[2] * kernel_2[0] == 2,
            "kernel basis independence control failed")
    # Rank five and two independent kernel vectors prove that these span the
    # full kernel of the 5-by-7 map.

    raw_kernel_1 = tuple(value * scalar
                         for value, scalar in zip(kernel_1, weights))
    raw_kernel_2 = tuple(value * scalar
                         for value, scalar in zip(kernel_2, weights))
    require(raw_kernel_1 == (-120, 120, -24, -48, 24, 0, 0),
            "first factorial raw-kernel vector drift")
    require(raw_kernel_2 == (120, -120, 0, 60, 0, -20, 4),
            "second factorial raw-kernel vector drift")

    # B=W_4 A W_5^{-1} is raw unweighted block deletion.  Construct it
    # exactly and verify that its (4)-row is the indicator of the upset
    # U={(5),(4,1)}.
    raw_deletion = []
    for target_weight, row in zip(target_weights, labelled):
        raw_row = tuple(Fraction(target_weight * entry, source_weight)
                        for entry, source_weight in zip(row, weights))
        require(all(value.denominator == 1 for value in raw_row),
                "factorial conjugacy lost integral raw deletion")
        raw_deletion.append(tuple(value.numerator for value in raw_row))
    raw_deletion = tuple(raw_deletion)
    expected_raw_deletion = (
        (1, 1, 0, 0, 0, 0, 0),
        (0, 1, 1, 2, 0, 0, 0),
        (0, 0, 1, 0, 1, 0, 0),
        (0, 0, 0, 1, 2, 3, 0),
        (0, 0, 0, 0, 0, 1, 5),
    )
    require(raw_deletion == expected_raw_deletion,
            "raw block-deletion matrix drift")
    require(matvec(raw_deletion, raw_kernel_1) == (0,) * 5
            and matvec(raw_deletion, raw_kernel_2) == (0,) * 5,
            "factorial raw-kernel vector escaped block deletion")

    upset = (0, 1)
    upset_mass = sum(raw[index] for index in upset)
    kernel_upset_masses = (
        sum(raw_kernel_1[index] for index in upset),
        sum(raw_kernel_2[index] for index in upset),
    )
    require(upset_mass == -1779840,
            "fixed negative upset mass drift")
    require(kernel_upset_masses == (0, 0),
            "full factorial kernel changed the blind upset")
    raw_image = matvec(raw_deletion, raw)
    labelled_image = matvec(labelled, factorial_coordinates)
    require(raw_image == tuple(weight_value * image_value
                               for weight_value, image_value
                               in zip(target_weights, labelled_image)),
            "factorial conjugacy failed on the hostile current")
    require(raw_image[0] == upset_mass,
            "upset mass is not the fixed (4)-deletion coordinate")

    print("active_prefix=I2:support=(1,3):N=5:j=d=5")
    print(f"reduced_poles={poles}:removed={removed}:P={polynomial}")
    print(f"dominant_Q={quotient_roots}:Phi_h={phi_h}:Q_h={quotient_h}:"
          f"Phi_p5={phi_power}:Q_p5={quotient_power}")
    print(f"raw_Lbar_current={raw}:sum={sum(raw)}")
    print(f"factorial_weights={weights}")
    print(f"factorial_coordinates={factorial_coordinates}")
    print(f"labelled_A={labelled}:rank=5:minor_01235=360")
    print(f"ker_A=K1:{kernel_1}:K2:{kernel_2}:dimension=2")
    print(f"raw_B={raw_deletion}")
    print(f"ker_B=WK1:{raw_kernel_1}:WK2:{raw_kernel_2}:dimension=2")
    print("blind_upset=((5),(4,1)):"
          f"Gmass={upset_mass}:WK1mass={kernel_upset_masses[0]}:"
          f"WK2mass={kernel_upset_masses[1]}")
    print(f"raw_deletion_image={raw_image}:target_(4)={raw_image[0]}")
    print("same_labelled_deletion_Hasse_preimage=NONE")
    print("scope=fixed_active_prefix;full_kernel;Hasse_cone;not_all_selectors")
    print("status=PASS")


if __name__ == "__main__":
    main()
