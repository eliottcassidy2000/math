#!/usr/bin/env python3
"""Exact coefficientwise blowup-carry controls for THM-3231."""

from __future__ import annotations

import ast
import hashlib
from fractions import Fraction
from math import comb
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = (
    ROOT
    / "01-canon"
    / "theorems"
    / "THM-3223-fourth-fifth-resonant-prs-primitive-walls-pell-content-and-pivot-resurrection.md"
)
EXPECTED_DEPENDENCY_SHA256 = (
    "2db785e5169b1c8b0eb414b4070765b6b76f47c0d02e92680eb206139c54dfb3"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


def pseudo(left: list[object], right: list[object]) -> list[object]:
    connection = left[0] * right[1] - left[1] * right[0]
    return [
        sp.expand(
            left[0] ** 2 * right[index + 2]
            - left[0] * right[0] * left[index + 2]
            - connection * left[index + 1]
        )
        for index in range(min(len(left), len(right)) - 2)
    ]


def p1(left: list[object], right: list[object]) -> object:
    return sp.expand(left[0] * right[1] - left[1] * right[0])


# Universal coefficientwise exceptional-square divisibility.
w = list(sp.symbols("w0:9"))
t = list(sp.symbols("t0:9"))
x = pseudo(w, t)
y = pseudo(x, w)


def exceptional_formula(index: int) -> object:
    return sp.expand(
        t[0]
        * (
            t[0] * w[1] ** 2 * w[index + 4]
            + (t[1] * w[1] ** 2 - t[0] * w[1] * w[2])
            * w[index + 3]
            + (
                t[0] * (w[2] ** 2 - w[1] * w[3])
                - t[1] * w[1] * w[2]
                + t[2] * w[1] ** 2
            )
            * w[index + 2]
            - t[index + 3] * w[1] ** 3
        )
    )


for index, coefficient in enumerate(y):
    polynomial = sp.Poly(coefficient, w[0])
    require(polynomial.nth(0) == 0, ("exceptional constant", index))
    require(polynomial.nth(1) == 0, ("exceptional tangent", index))
    require(
        sp.expand(polynomial.nth(2) - exceptional_formula(index)) == 0,
        ("exceptional quotient", index),
    )
    determinant = sp.det(
        sp.Matrix(
            [
                [t[0], t[1], t[2], t[index + 3]],
                [w[1], w[2], w[3], w[index + 4]],
                [0, w[1], w[2], w[index + 3]],
                [0, 0, w[1], w[index + 2]],
            ]
        )
    )
    require(
        sp.expand(exceptional_formula(index) - t[0] * determinant) == 0,
        ("osculating determinant", index),
    )

# The normal-cone constant is the next reciprocal quotient jet of the shifted
# wall row v=z^(-1)(w-w_0).
v = w[1:]
ev = pseudo(v, t)
normal_constant = p1(ev, v)
require(
    sp.expand(exceptional_formula(0) - t[0] * normal_constant) == 0,
    "normal-cone Pluecker constant",
)
ratio = [sp.cancel(t[0] / v[0])]
for index in range(1, 8):
    ratio.append(
        sp.cancel(
            (t[index] - sum(v[j] * ratio[index - j] for j in range(1, index + 1)))
            / v[0]
        )
    )
require(
    sp.cancel(normal_constant + v[0] ** 4 * ratio[3]) == 0,
    "third reciprocal quotient jet",
)
for index in range(len(y)):
    convolution = sum(v[j] * ratio[index - j + 3] for j in range(index + 1))
    require(
        sp.cancel(
            exceptional_formula(index)
            + t[0] * v[0] ** 3 * convolution
        )
        == 0,
        ("quotient-tail convolution", index),
    )


# Exact p=43, d=p+2 application through the sixth constant pivot.
P = 43
CAP = 10
MODULUS = P**CAP
D = P + 2
A_DEGREE = D - 2


def reciprocal_top(epsilon: int, codimension: int) -> int:
    total = Fraction(0)
    for ell in range(codimension + 1):
        offset = 2 * epsilon - 2 * codimension + ell
        if offset >= 0:
            ratio = 1
            for j in range(1, offset + 1):
                ratio *= 2 * A_DEGREE + j
        else:
            denominator = 1
            for j in range(0, -offset):
                denominator *= 2 * A_DEGREE - j
            ratio = Fraction(1, denominator)
        total += (
            comb(codimension, ell)
            * D ** (codimension - ell)
            * (-1) ** ell
            * ratio
        )
    value = comb(A_DEGREE + epsilon, codimension) * total
    require(
        value.denominator % P != 0,
        ("nonunit reduced denominator", epsilon, codimension),
    )
    return (
        value.numerator
        * pow(value.denominator, -1, MODULUS)
        % MODULUS
    )


def pseudo_mod(left: list[int], right: list[int]) -> list[int]:
    connection = (left[0] * right[1] - left[1] * right[0]) % MODULUS
    return [
        (
            left[0] ** 2 * right[index + 2]
            - left[0] * right[0] * left[index + 2]
            - connection * left[index + 1]
        )
        % MODULUS
        for index in range(min(len(left), len(right)) - 2)
    ]


def valuation(value: int) -> int:
    if value == 0:
        return CAP
    order = 0
    while value % P == 0:
        value //= P
        order += 1
    return order


left = [reciprocal_top(0, index) for index in range(13)]
right = [reciprocal_top(1, index) for index in range(13)]
rows = [pseudo_mod(left, right)]
rows.append(pseudo_mod(rows[-1], left))
for _ in range(4):
    rows.append(pseudo_mod(rows[-1], rows[-2]))

require(tuple(map(len, (left, *rows))) == (13, 11, 9, 7, 5, 3, 1), "lengths")
constant_valuations = tuple(valuation(row[0]) for row in rows)
require(constant_valuations == (0, 0, 0, 1, 0, 2), constant_valuations)
require(tuple(value % P for value in rows[3][:5]) == (0, 24, 20, 16, 3), "row4")
require(tuple(value % P for value in rows[4][:3]) == (23, 12, 1), "row5")

w0_unit = (rows[3][0] // P) % P
y0_unit = (rows[5][0] // P**2) % P
require(w0_unit == 36 and y0_unit == 19, (w0_unit, y0_unit))


def exceptional_constant_mod(w_row: list[int], t_row: list[int]) -> int:
    w1, w2, w3, w4 = (value % P for value in w_row[1:5])
    t0, t1, t2, t3 = (value % P for value in t_row[:4])
    return (
        t0
        * (
            t0 * w1**2 * w4
            - 2 * t0 * w1 * w2 * w3
            + t0 * w2**3
            + t1 * w1**2 * w3
            - t1 * w1 * w2**2
            + t2 * w1**2 * w2
            - t3 * w1**3
        )
        % P
    )


exceptional_unit = exceptional_constant_mod(rows[3], rows[2])
require(exceptional_unit == 39, exceptional_unit)
require(w0_unit**2 * exceptional_unit % P == y0_unit, "blowup carry")

source = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source)
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert")
require(
    not any(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    ),
    "float literal",
)
require(lf_sha256(DEPENDENCY) == EXPECTED_DEPENDENCY_SHA256, "dependency drift")

print("THM-3231 exact simple-pivot exceptional-square carry")
print("dependency_pin=THM-3223")
print("symbolic_row_coefficients=5 constant_and_linear_w0_terms=0 PASS")
print("exceptional_formula_controls=5 PASS")
print("toeplitz_4x4_determinants=5 quotient_tail_convolutions=5 PASS")
print("normal_cone_constant=t0*P1(E(Sw,t),Sw) third_quotient_jet=PASS")
print("p43_jet_lengths=13,11,9,7,5,3,1")
print("p43_constant_valuations=0,0,0,1,0,2")
print("p43_row4_mod43=0,24,20,16,3 row5_mod43=23,12,1")
print("p43_units=w0_over_43:36 exceptional:39 row6_over_43^2:19")
print("mod43_terminal_row_lifts_to_second_order_unit=PASS")
print("scope=exceptional_pivot_carry_not_full_prs_or_root_selector")
