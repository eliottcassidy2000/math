#!/usr/bin/env python3
"""Exact companion for THM-4010, the higher-jet consecutive sampler."""

from __future__ import annotations

from hashlib import sha256
import json
from math import comb, factorial
import sys

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


sys.stdout.reconfigure(newline="\n")
GATES = 0
EXPECTED_SEMANTIC_SHA256 = "ac4767298585cf7aaed247899cec0e675837b55d57d279d0295681484ee7c0fa"


def require(condition, label):
    global GATES
    GATES += 1
    if not bool(condition):
        raise RuntimeError(label)


def hasse_matrix(a: int, m: int, k: int) -> sp.Matrix:
    size = (m + 1) * k
    return sp.Matrix([
        [
            0 if degree < order else
            comb(degree, order) * (a + node) ** (degree - order)
            for degree in range(size)
        ]
        for node in range(m + 1)
        for order in range(k)
    ])


def translation_matrix(a: int, size: int) -> sp.Matrix:
    # Coefficients of P(a+Y) from coefficients of P(X).
    return sp.Matrix([
        [0 if degree < power else comb(degree, power) * a ** (degree - power)
         for degree in range(size)]
        for power in range(size)
    ])


def newton_basis_matrix(m: int, k: int) -> sp.Matrix:
    x = sp.symbols("x")
    size = (m + 1) * k
    columns = []
    for node in range(m + 1):
        prior = sp.prod((x - left) ** k for left in range(node))
        for order in range(k):
            poly = sp.Poly((x - node) ** order * prior, x, domain=sp.ZZ)
            columns.append([int(poly.nth(degree)) for degree in range(size)])
    return sp.Matrix(size, size, lambda row, col: columns[col][row])


def invariant_factors(matrix: sp.Matrix) -> tuple[int, ...]:
    diagonal = smith_normal_form(matrix, domain=ZZ)
    return tuple(abs(int(diagonal[i, i])) for i in range(diagonal.rows))


def determinant_formula(m: int, k: int) -> int:
    return sp.prod(factorial(node) ** (k * k) for node in range(1, m + 1))


def naive_block_smith(m: int, k: int) -> tuple[int, ...]:
    return tuple(
        factorial(node) ** k
        for node in range(m + 1)
        for _ in range(k)
    )


atlas = []
first_naive_failure = None

for k in range(1, 5):
    for m in range(0, 6):
        size = (m + 1) * k
        base = hasse_matrix(0, m, k)
        factors = invariant_factors(base)
        determinant = abs(int(base.det()))
        expected = int(determinant_formula(m, k))
        require(determinant == expected, ("determinant formula", m, k))
        require(sp.prod(factors) == expected, ("Smith product", m, k))
        require(all(right % left == 0 for left, right in zip(factors, factors[1:])),
                ("Smith divisibility", m, k))

        basis = newton_basis_matrix(m, k)
        require(abs(int(basis.det())) == 1, ("Newton basis unimodular", m, k))
        triangular = base * basis
        for node in range(m + 1):
            for order in range(k):
                index = node * k + order
                require(triangular[index, index] == factorial(node) ** k,
                        ("Newton diagonal", m, k, node, order))
                for later_node in range(node + 1, m + 1):
                    for later_order in range(k):
                        column = later_node * k + later_order
                        require(triangular[index, column] == 0,
                                ("Newton block triangular", m, k, index, column))

        for shift in (-3, 2, 5):
            translated = hasse_matrix(shift, m, k)
            transform = translation_matrix(shift, size)
            require(abs(int(transform.det())) == 1,
                    ("translation unimodular", shift, m, k))
            require(translated == base * transform,
                    ("translation identity", shift, m, k))
            require(invariant_factors(translated) == factors,
                    ("translation Smith", shift, m, k))

        naive = naive_block_smith(m, k)
        if first_naive_failure is None and factors != naive:
            first_naive_failure = (m, k, naive, factors)
        atlas.append((m, k, size, factors, determinant))

require(first_naive_failure is not None, "naive block Smith hostile exists")
require(first_naive_failure[:2] == (3, 2), "first scan-order naive hostile")
require(first_naive_failure[2] == (1, 1, 1, 1, 4, 4, 36, 36),
        "naive hostile tuple")
require(first_naive_failure[3] == (1, 1, 1, 1, 4, 4, 12, 108),
        "actual hostile tuple")

# Direct finite controls for the principal kernel and optimal evaluation
# modulus. The general proof is monic division plus Hasse Taylor expansion.
x = sp.symbols("x")
kernel_controls = []
for a, m, k in ((-2, 2, 2), (1, 3, 2), (-1, 2, 3), (0, 4, 2)):
    size = (m + 1) * k
    F = sp.prod(x - (a + node) for node in range(m + 1))
    divisor = sp.Poly(F ** k, x, domain=sp.ZZ)
    P = sp.Poly(
        sum(((-1) ** degree) * (degree + 2) * x ** degree
            for degree in range(size + 6)),
        x, domain=sp.ZZ,
    )
    quotient, remainder = sp.div(P, divisor, domain=sp.ZZ)
    require(P == quotient * divisor + remainder, ("monic division", a, m, k))
    require(remainder.degree() < size, ("remainder degree", a, m, k))
    for node in range(m + 1):
        point = a + node
        for order in range(k):
            delta = sp.diff(P.as_expr() - remainder.as_expr(), x, order)
            require(delta.subs(x, point) == 0,
                    ("ordinary/Hasse jet equality", a, m, k, node, order))
    for B in (a - 2, a + m + 2, a + m + 5):
        modulus = int(F.subs(x, B)) ** k
        difference = int(P.eval(B) - remainder.eval(B))
        require(modulus != 0 and difference % modulus == 0,
                ("evaluation congruence", a, m, k, B))
        require(int(divisor.eval(B)) == modulus,
                ("sharp generator", a, m, k, B))
    kernel_controls.append((a, m, k, tuple(int(c) for c in remainder.all_coeffs())))

semantic_record = {
    "atlas": [
        [m, k, size, list(factors), str(determinant)]
        for m, k, size, factors, determinant in atlas
    ],
    "first_naive_failure": [
        first_naive_failure[0], first_naive_failure[1],
        list(first_naive_failure[2]), list(first_naive_failure[3]),
    ],
    "kernel_controls": [
        [a, m, k, list(coefficients)]
        for a, m, k, coefficients in kernel_controls
    ],
    "theorem": "ker(jet)=F^k*Z[X];optimal_modulus=abs(F(B))^k",
    "index": "product_{j=1}^m(j!)^(k^2)",
}
semantic_sha256 = sha256(
    json.dumps(semantic_record, sort_keys=True,
               separators=(",", ":")).encode("ascii")
).hexdigest()
if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

print("CONFLUENT_CONSECUTIVE_HASSE_OBSERVER_THM4010")
print("status=FINITE-EXACT_CANONICAL;THM4000_EXTENSION;NO_UNSTATED_GENERAL_SNF_FORMULA")
print("theorem=equal_Hasse_jets_iff_difference_in_F^k_Z[X]")
print("optimal_new_base_modulus=abs(F(B))^k_for_B_outside_nodes")
print("jet_lattice_index=product_(j=1)^m(j!)^(k^2)")
print("translation_invariance=verified_by_unimodular_binomial_shift")
print("newton_basis=block_triangular_with_(j!)^k_repeated_k_on_block_diagonal")
print("naive_smith_hostile=" + repr(first_naive_failure).replace(" ", ""))
print("atlas=" + repr(tuple(atlas)).replace(" ", ""))
print("semantic_sha256=" + semantic_sha256)
print(f"gates={GATES}")
print("RESULT=EXACT_KERNEL_MODULUS_AND_INDEX;FULL_GENERAL_SNF_OPEN")
