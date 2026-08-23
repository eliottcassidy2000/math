#!/usr/bin/env python3
"""Finite-exact shifted three-Cohn sidecar for THM-3822.

This script does not prove an all-word Cohn theorem.  It exhausts the stated
finite universe and uses a one-variable parity certificate that is sufficient
to disprove the necessary square condition H(h,k) in each cell.
"""

from __future__ import annotations

import ast
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


x = sp.symbols("x")
prime = 1009
zero = sp.Poly(0, x, modulus=prime)
one = sp.Poly(1, x, modulus=prime)


def cohn(a: int, b: int, y_value: int):
    """Cohn matrix C(x+a,y+b), already specialized modulo prime."""
    X = sp.Poly(x + a, x, modulus=prime)
    T = sp.Poly(y_value + b, x, modulus=prime)
    return ((4 * T * T, 2 * X * T - one),
            (one + 2 * X * T, X * X))


def inverse(matrix):
    return ((matrix[1][1], -matrix[0][1]),
            (-matrix[1][0], matrix[0][0]))


def multiply(left, right):
    return tuple(
        tuple(
            sum((left[i][r] * right[r][j] for r in range(2)), zero)
            for j in range(2)
        )
        for i in range(2)
    )


def determinant(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def atlas_discriminant(h, k):
    """The square sidecar H(h,k) from the THM-3811 SL2 reconstruction."""
    return (
        84 * h**7 + 36 * h**6 * k**2 + 196 * h**6 * k
        + 84 * h**5 * k**3 + 36 * h**5 * k**2 + 49 * h**4 * k**4
        + 112 * h**4 * k**3 - 12 * h**3 * k**5 - 14 * h**2 * k**6
        + 12 * h**2 * k**5 + k**8
    )


def certified_nonsquare(polynomial: sp.Poly) -> bool:
    """A sufficient nonsquare certificate over an algebraic closure of F_p.

    If f=g^2 and deg(f)=2d, then gcd(f,f') has degree at least d.  Therefore
    odd degree, or a derivative gcd of degree below half, proves that f is not
    a square even after extending the constant field.  Returning False makes
    no claim.
    """
    if polynomial.is_zero:
        return False
    degree = polynomial.degree()
    if degree % 2:
        return True
    return sp.gcd(polynomial, polynomial.diff()).degree() < degree // 2


# Hostile and positive controls for the deliberately one-sided certificate.
square_control = sp.Poly((x**2 + 3 * x + 1) ** 2, x, modulus=prime)
nonsquare_control = sp.Poly(x**4 + x + 1, x, modulus=prime)
check(not certified_nonsquare(square_control), "square control survives")
check(certified_nonsquare(nonsquare_control), "nonsquare control is rejected")


shifts = (-1, 0, 1)
specializations = (2, 3)
factor_labels = list(itertools.product(shifts, shifts, (1, -1)))
check(len(factor_labels) == 18, "eighteen shifted/inverse factors")

matrices = {}
for y_value in specializations:
    row = []
    for a, b, exponent in factor_labels:
        matrix = cohn(a, b, y_value)
        if exponent == -1:
            matrix = inverse(matrix)
        check(determinant(matrix) == one, "each factor has determinant one")
        row.append(matrix)
    matrices[y_value] = row


pattern_counts: Counter[tuple[bool, bool]] = Counter()
survivors = []
for word in itertools.product(range(len(factor_labels)), repeat=3):
    rejection_pattern = []
    for y_value in specializations:
        first = matrices[y_value][word[0]]
        second = matrices[y_value][word[1]]
        third = matrices[y_value][word[2]]
        product = multiply(multiply(first, second), third)
        check(determinant(product) == one, "each word has determinant one")
        k, h = product[0]
        candidate = atlas_discriminant(h, k)
        rejection_pattern.append(certified_nonsquare(candidate))
    pattern = tuple(rejection_pattern)
    pattern_counts[pattern] += 1
    if not any(pattern):
        survivors.append(word)


total = len(factor_labels) ** 3
check(total == 5832, "finite universe size")
check(pattern_counts == Counter({(True, True): 5832}),
      "both specializations reject every word")
check(not survivors, "no finite-universe survivor")


source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "scope": "FINITE-EXACT only: ordered products of three C(x+a,y+b)^eps factors",
    "shifts": "a,b in {-1,0,1}; eps in {+1,-1}",
    "universe": total,
    "field": "GF(1009)",
    "fibres": list(specializations),
    "criterion": "odd degree or deg gcd(f,f') < deg(f)/2 for H(h,k)",
    "patterns": {str(key): value for key, value in sorted(pattern_counts.items())},
    "survivors": len(survivors),
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3822-finite-exact-shifted-triple-Cohn-sidecar")
print("scope=C(x+a,y+b)^eps;ordered_length=3;a,b=-1,0,1;eps=+-1")
print(f"field=GF({prime});specializations=y=2,3")
print(f"factors={len(factor_labels)};words={total}")
print("rejection=odd_degree_or_derivative_gcd_below_half_degree")
print(f"both_fibres_reject={pattern_counts[(True, True)]}")
print(f"survivors={len(survivors)}")
print("scope_warning=FINITE-EXACT_NOT_ALL_WORDS_NOT_A_JC_THEOREM")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
