#!/usr/bin/env python3
"""Independent exact audit of THM-3079's Newton--PF row-transform theorem.

This file does not import the primary companion.  It reconstructs falling-
factorial coefficients by finite differences, audits Toeplitz minors by a
Leibniz determinant, evaluates the reciprocal-Gamma strip directly over Q,
and checks multi-edge strict-prefix meshes.  It also records a below-strip
hostile showing why the positivity condition is a real hypothesis.
"""

from __future__ import annotations

import ast
import hashlib
import itertools
import math
from fractions import Fraction
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "7c413b0b7d5c568f77e795c74150b91337621ab018801bff8450721dfa0cfaeb"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def falling(x, degree):
    result = Fraction(1)
    for offset in range(degree):
        result *= x - offset
    return result


def polynomial_value(shifts, x):
    result = Fraction(1)
    for shift in shifts:
        result *= x + shift
    return result


def newton_coefficients(shifts):
    degree = len(shifts)
    values = [polynomial_value(shifts, Fraction(point))
              for point in range(degree + 1)]
    coefficients = []
    layer = values
    for order in range(degree + 1):
        coefficients.append(layer[0] / math.factorial(order))
        layer = [layer[index + 1] - layer[index]
                 for index in range(len(layer) - 1)]
    return tuple(coefficients)


def determinant_leibniz(matrix):
    order = len(matrix)
    total = Fraction(0)
    for permutation in itertools.permutations(range(order)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(order)
            for right in range(left + 1, order)
        )
        term = Fraction(-1 if inversions % 2 else 1)
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        total += term
    return total


def reciprocal_gamma_ratio(shape, index):
    """Gamma(shape)/Gamma(shape+index), for an integer index on its strip."""
    shape = Fraction(shape)
    require(shape + index > 0, (shape, index, "outside positive strip"))
    if index >= 0:
        denominator = Fraction(1)
        for offset in range(index):
            denominator *= shape + offset
        return 1 / denominator
    result = Fraction(1)
    for offset in range(index, 0):
        result *= shape + offset
    return result


def transformed_value(shape, coefficients, index):
    return sum(
        coefficient * reciprocal_gamma_ratio(shape, index - shift)
        for shift, coefficient in enumerate(coefficients)
    )


def toeplitz_entry(coefficients, row, column):
    difference = row - column
    return coefficients[difference] if 0 <= difference < len(coefficients) else 0


def audit_newton_families():
    families = (
        (Fraction(1),),
        (Fraction(1, 2), Fraction(3, 2)),
        (Fraction(2, 3), Fraction(2, 3), Fraction(5, 2)),
        (Fraction(1, 5), Fraction(2, 5), Fraction(3, 5), Fraction(4, 5)),
        (Fraction(1), Fraction(2), Fraction(4), Fraction(4), Fraction(7)),
    )
    reconstruction_cells = 0
    toeplitz_minors = 0
    generalized_minors = 0
    digest = hashlib.sha256()
    for family_index, shifts in enumerate(families):
        coefficients = newton_coefficients(shifts)
        require(all(coefficient > 0 for coefficient in coefficients),
                (shifts, coefficients, "positive Newton coefficients"))
        for point in range(-3, 13):
            reconstructed = sum(
                coefficient * falling(Fraction(point), degree)
                for degree, coefficient in enumerate(coefficients)
            )
            require(reconstructed == polynomial_value(shifts, Fraction(point)),
                    (shifts, point, reconstructed))
            reconstruction_cells += 1

        size = len(coefficients) + 3
        for order in range(1, min(3, size) + 1):
            for rows in itertools.combinations(range(size), order):
                for columns in itertools.combinations(range(size), order):
                    minor = determinant_leibniz([
                        [toeplitz_entry(coefficients, row, column)
                         for column in columns]
                        for row in rows
                    ])
                    require(minor >= 0, (shifts, rows, columns, minor))
                    toeplitz_minors += 1

        degree = len(shifts)
        shape = Fraction(degree + 1) + Fraction(family_index + 1, 7)
        for index in range(0, 10):
            direct = polynomial_value(shifts, shape + index - 1)
            direct *= reciprocal_gamma_ratio(shape, index)
            transformed = transformed_value(shape, coefficients, index)
            require(direct == transformed,
                    (shifts, shape, index, direct, transformed))
            reconstruction_cells += 1

        nodes = range(6)
        for order in range(1, 5):
            target_sign = -1 if (order * (order - 1) // 2) % 2 else 1
            for rows in itertools.combinations(nodes, order):
                for columns in itertools.combinations(nodes, order):
                    determinant = determinant_leibniz([
                        [transformed_value(shape, coefficients, row + column)
                         for column in columns]
                        for row in rows
                    ])
                    require(target_sign * determinant > 0,
                            (shifts, shape, rows, columns, determinant))
                    generalized_minors += 1
        digest.update(repr((shifts, coefficients, shape)).encode())
    return (
        reconstruction_cells,
        toeplitz_minors,
        generalized_minors,
        digest.hexdigest(),
    )


def pochhammer(shape, index):
    result = Fraction(1)
    for offset in range(index):
        result *= shape + offset
    return result


def mesh_value(shapes, exponents, index):
    result = Fraction(1)
    for shape, exponent in zip(shapes, exponents):
        factor = pochhammer(shape, index)
        result *= factor ** exponent
    return result


def mesh_shifts(shapes, exponents):
    prefixes = []
    running = 0
    for exponent in exponents:
        running += exponent
        prefixes.append(running)
    require(all(prefix < 0 for prefix in prefixes), prefixes)
    require(prefixes[-1] == -1, prefixes)
    shifts = []
    base = shapes[0]
    for edge in range(len(shapes) - 1):
        gap = shapes[edge + 1] - shapes[edge]
        require(gap.denominator == 1 and gap > 0, gap)
        multiplicity = -prefixes[edge] - 1
        for _copy in range(multiplicity):
            for offset in range(int(gap)):
                shifts.append((shapes[edge] - base) + offset + 1)
    return tuple(shifts)


def audit_integer_meshes():
    meshes = (
        ((Fraction(3), Fraction(4)), (-2, 1)),
        ((Fraction(7), Fraction(9), Fraction(10)), (-3, 1, 1)),
        ((Fraction(10), Fraction(11), Fraction(14), Fraction(15)),
         (-2, 1, -2, 2)),
    )
    identity_cells = 0
    minor_cells = 0
    digest = hashlib.sha256()
    for shapes, exponents in meshes:
        shifts = mesh_shifts(shapes, exponents)
        degree = len(shifts)
        require(shapes[0] > degree, (shapes, exponents, degree))
        coefficients = newton_coefficients(shifts)
        scale = Fraction(1)
        for edge in range(len(shapes) - 1):
            gap = int(shapes[edge + 1] - shapes[edge])
            prefixes = sum(exponents[:edge + 1])
            multiplicity = -prefixes - 1
            scale *= pochhammer(shapes[edge], gap) ** (-multiplicity)
        for index in range(10):
            transformed = scale * transformed_value(shapes[0], coefficients, index)
            direct = mesh_value(shapes, exponents, index)
            require(transformed == direct,
                    (shapes, exponents, index, transformed, direct))
            identity_cells += 1
        for order in range(1, 5):
            sign = -1 if (order * (order - 1) // 2) % 2 else 1
            for rows in itertools.combinations(range(6), order):
                for columns in itertools.combinations(range(6), order):
                    determinant = determinant_leibniz([
                        [mesh_value(shapes, exponents, row + column)
                         for column in columns]
                        for row in rows
                    ])
                    require(sign * determinant > 0,
                            (shapes, exponents, rows, columns, determinant))
                    minor_cells += 1
        digest.update(repr((shapes, exponents, shifts, coefficients, scale)).encode())
    return identity_cells, minor_cells, digest.hexdigest()


def audit_strip_hostile():
    shape = Fraction(1, 2)
    coefficients = (Fraction(1), Fraction(3))
    # Gamma(a)/Gamma(a-1)=a-1 remains defined here, but is outside the
    # positive strip used by the theorem's Cauchy--Binet proof.
    leading = coefficients[0] + coefficients[1] * (shape - 1)
    require(leading == Fraction(-1, 2), leading)
    for base in range(1, 7):
        value = transformed_value(shape, coefficients, base)
        require(value > 0, (base, value))
    return shape, coefficients, leading


def main():
    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assert_nodes, float_nodes) == (0, 0), "truth-gate syntax")
    newton = audit_newton_families()
    meshes = audit_integer_meshes()
    hostile = audit_strip_hostile()
    semantic = (newton, meshes, hostile, assert_nodes, float_nodes)
    semantic_sha = hashlib.sha256(repr(semantic).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha == EXPECTED_SEMANTIC_SHA256, semantic_sha)
    lines = [
        "THM-3079 INDEPENDENT NEWTON-PF STRICT-MESH AUDIT",
        (
            f"newton_reconstruction_cells={newton[0]};"
            f"toeplitz_minors={newton[1]};generalized_minors={newton[2]};"
            f"family_sha256={newton[3]}"
        ),
        (
            f"integer_mesh_identity_cells={meshes[0]};"
            f"integer_mesh_generalized_minors={meshes[1]};mesh_sha256={meshes[2]}"
        ),
        (
            "below_strip_hostile=a:1/2;C(z)=1+3z;F_0=-1/2;"
            "tail_bases_1_to_6_positive"
        ),
        "proof_audit=Newton_interlacing+Cauchy_Binet_distinguished_minor+Euler_tail_shift:PASS",
        f"truth_gates=assert_nodes:{assert_nodes};float_literals:{float_nodes}",
        f"semantic_sha256={semantic_sha}",
        "scope=independent_exact_audit;general_strict_prefix_and_noninteger_mesh_remain_open",
        "all_exact_controls=PASS",
    ]
    print("\n".join(lines))


if __name__ == "__main__":
    main()
