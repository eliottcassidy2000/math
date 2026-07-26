#!/usr/bin/env python3
"""Exact companion for THM-2408.

The checks use only integers and Fraction arithmetic.  Roots in the exhaustive
Prony controls are rational test nodes; the proof in the theorem is the same
Vandermonde argument over C.
"""

from fractions import Fraction
from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def norm2(z):
    return z[0] * z[0] + z[1] * z[1]


def add_complex(z, w):
    return (z[0] + w[0], z[1] + w[1])


def poly_from_roots(roots):
    """Coefficients in increasing order."""
    coefficients = [Fraction(1)]
    for root in roots:
        updated = [Fraction(0)] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            updated[degree] -= root * coefficient
            updated[degree + 1] += coefficient
        coefficients = updated
    return coefficients


def poly_value(coefficients, x):
    value = Fraction(0)
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def exponential_value(bank, n):
    return sum(
        (coefficient * (node ** n) for node, coefficient in bank.items()),
        Fraction(0),
    )


def apply_shift_polynomial(coefficients, sequence, n):
    return sum(
        (coefficient * sequence(n + degree)
         for degree, coefficient in enumerate(coefficients)),
        Fraction(0),
    )


def resultant_from_roots(left, right):
    value = Fraction(1)
    for alpha in left:
        for beta in right:
            value *= alpha - beta
    return value


def cyclotomic13_remainder(values, character):
    """Evaluate a rational C_13 profile in the basis 1,zeta,...,zeta^11."""

    coefficients = [Fraction(0)] * 13
    for exponent, value in enumerate(values):
        coefficients[(character * exponent) % 13] += value
    top = coefficients[12]
    return tuple(coefficients[index] - top for index in range(12))


def exact_hilbert_checks():
    values = range(-2, 3)
    gaussian = list(product(values, repeat=2))
    cases = 0
    equality_seen = False
    for owner in gaussian:
        for difference in gaussian:
            union = add_complex(owner, difference)
            owner_norm = norm2(owner)
            difference_norm = norm2(difference)
            union_norm = norm2(union)
            require(
                4 * max(owner_norm, union_norm) >= difference_norm,
                "pointwise branch bound failed",
            )
            require(
                2 * (owner_norm + union_norm) >= difference_norm,
                "summed Hilbert bound failed",
            )
            if (
                difference != (0, 0)
                and 4 * owner_norm == difference_norm
                and 4 * union_norm == difference_norm
            ):
                equality_seen = True
            cases += 1
    require(equality_seen, "sharp midpoint equality not represented")
    return cases


def exhaustive_prony_checks():
    nodes = (Fraction(1), Fraction(2), Fraction(3))
    coefficient_choices = (0, -2, -1, 1, 2)
    banks = []
    for values in product(coefficient_choices, repeat=len(nodes)):
        banks.append(
            {
                node: Fraction(coefficient)
                for node, coefficient in zip(nodes, values)
                if coefficient
            }
        )

    pair_count = 0
    nonzero_residual_count = 0
    transverse_count = 0
    for owner in banks:
        for difference in banks:
            union = {
                node: owner.get(node, 0) + difference.get(node, 0)
                for node in nodes
                if owner.get(node, 0) + difference.get(node, 0)
            }
            residual_nodes = tuple(union)
            residual_count = len(residual_nodes)
            if residual_count:
                nonzero_residual_count += 1
                for start in range(5):
                    window = [
                        exponential_value(union, start + offset)
                        for offset in range(residual_count)
                    ]
                    require(
                        any(value != 0 for value in window),
                        "Vandermonde window failed",
                    )
            else:
                for n in range(7):
                    require(
                        exponential_value(owner, n)
                        == -exponential_value(difference, n),
                        "shared-node cancellation classification failed",
                    )

            owner_polynomial = poly_from_roots(tuple(owner))
            transverse = {
                node: coefficient
                for node, coefficient in difference.items()
                if node not in owner
            }
            if transverse:
                transverse_count += 1

            def union_sequence(n):
                return exponential_value(owner, n) + exponential_value(
                    difference, n
                )

            for n in range(5):
                left = apply_shift_polynomial(
                    owner_polynomial, union_sequence, n
                )
                right = sum(
                    (
                        coefficient
                        * poly_value(owner_polynomial, node)
                        * (node ** n)
                        for node, coefficient in transverse.items()
                    ),
                    Fraction(0),
                )
                require(left == right, "shift-operator identity failed")

            if transverse:
                window_length = len(owner) + len(transverse)
                for start in range(5):
                    require(
                        any(
                            union_sequence(start + offset) != 0
                            for offset in range(window_length)
                        ),
                        "transverse r+t window failed",
                    )

            pair_count += 1

    return pair_count, nonzero_residual_count, transverse_count


def resultant_checks():
    nodes = (Fraction(1), Fraction(2), Fraction(3))
    subsets = []
    for mask in range(1 << len(nodes)):
        subsets.append(
            tuple(
                node for index, node in enumerate(nodes)
                if (mask >> index) & 1
            )
        )
    count = 0
    for left in subsets:
        for right in subsets:
            resultant = resultant_from_roots(left, right)
            require(
                (resultant != 0) == set(left).isdisjoint(right),
                "resultant/common-node equivalence failed",
            )
            count += 1
    return count


def sharp_window_check():
    nodes = tuple(map(Fraction, (1, 2, 4, 7)))
    polynomial = poly_from_roots(nodes)
    coefficients = {}
    for node in nodes:
        derivative = sum(
            degree * polynomial[degree] * (node ** (degree - 1))
            for degree in range(1, len(polynomial))
        )
        require(derivative != 0, "repeated sharpness node")
        coefficients[node] = Fraction(1, 1) / derivative
    values = tuple(exponential_value(coefficients, n) for n in range(4))
    require(values == (0, 0, 0, 1), "sharp Lagrange window failed")
    return values


def positivity_hostile_check():
    p = 13
    epsilon = Fraction(1, 2 * p)
    h_nonzero_mode = epsilon / p
    owner_nonzero_mode = -h_nonzero_mode
    union_nonzero_mode = Fraction(0)
    owner_values = [
        Fraction(1, p) - (epsilon if shift == 0 else 0)
        for shift in range(p)
    ]
    require(all(value > 0 for value in owner_values), "owner not strict")
    require(h_nonzero_mode == Fraction(1, 338), "H hostile mode")
    require(owner_nonzero_mode == Fraction(-1, 338), "O hostile mode")
    require(union_nonzero_mode == 0, "U hostile mode")
    energy = 12 * h_nonzero_mode * h_nonzero_mode
    require(energy == Fraction(3, 28561), "hostile energy")
    return h_nonzero_mode, owner_nonzero_mode, energy


def endpoint_hostile_check():
    """Check the two shared oriented jumps at 0 and 1/13 exactly."""

    h_jumps = [Fraction(0)] * 13
    h_jumps[0] = Fraction(1)
    h_jumps[1] = Fraction(-1)
    owner_jumps = tuple(-value for value in h_jumps)
    checked = 0
    for residue in range(1, 13):
        h_remainder = cyclotomic13_remainder(h_jumps, -residue)
        owner_remainder = cyclotomic13_remainder(owner_jumps, -residue)
        require(any(h_remainder), "endpoint hostile H node sum vanished")
        require(
            owner_remainder == tuple(-value for value in h_remainder),
            "endpoint hostile lost opposite coefficients",
        )
        checked += 1
    require(checked == 12, "endpoint hostile residue count")
    return checked


def centered_bound_check():
    modulus = 13
    length = 5
    lower = -(length // 2)
    upper = (length + 1) // 2
    frequencies = [
        k + modulus * n
        for k in range(1, modulus)
        for n in range(lower, upper)
    ]
    bound = modulus * ((length + 1) // 2) - 1
    require(max(map(abs, frequencies)) == bound, "centered bound not sharp")
    return bound


def main():
    hilbert_cases = exact_hilbert_checks()
    pair_count, residual_count, transverse_count = exhaustive_prony_checks()
    resultant_count = resultant_checks()
    sharp_values = sharp_window_check()
    h_mode, owner_mode, hostile_energy = positivity_hostile_check()
    endpoint_residues = endpoint_hostile_check()
    centered_bound = centered_bound_check()

    thm2403_energy = Fraction(9, 250_994_068) / 4
    thm2403_amplitude_denominator = 228_488 * 2
    require(
        thm2403_energy == Fraction(9, 1_003_976_272),
        "THM-2403 branch energy",
    )
    require(
        thm2403_amplitude_denominator == 456_976,
        "THM-2403 branch amplitude",
    )

    print(f"hilbert_gaussian_cases: {hilbert_cases}")
    print("branch_energy_constant: 1/4")
    print("branch_sum_constant: 1/2")
    print("branch_equality: O=-H/2,U=H/2")
    print(f"prony_bank_pairs: {pair_count}")
    print(f"prony_nonzero_residual_pairs: {residual_count}")
    print(f"prony_transverse_pairs: {transverse_count}")
    print(f"resultant_subset_pairs: {resultant_count}")
    print(
        "sharp_four_node_sequence: "
        + ",".join(str(value) for value in sharp_values)
    )
    print(f"positive_hostile_H_mode: {h_mode}")
    print(f"positive_hostile_O_mode: {owner_mode}")
    print(f"positive_hostile_energy: {hostile_energy}")
    print(f"endpoint_hostile_nonzero_residues: {endpoint_residues}")
    print(f"centered_m13_length5_bound: {centered_bound}")
    print(f"thm2403_one_branch_energy: {thm2403_energy}*rho_R^2")
    print(
        "thm2403_one_branch_amplitude: "
        f"1/{thm2403_amplitude_denominator}*rho_R"
    )


if __name__ == "__main__":
    main()
